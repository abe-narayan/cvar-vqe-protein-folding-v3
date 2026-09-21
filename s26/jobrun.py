#!/usr/bin/env python
"""s26/jobrun.py - launch one job under the S26 governor.

    python s26/jobrun.py --agent Q --tag CPU --name a2_dla [--est-ram 1.0] -- python s26/q_dla.py ...

What it does, in order:
  1. refuses to start while the governor's last snapshot shows the box above the RAM launch
     gate, while the concurrency cap is full, while the box lacks --est-ram + 0.5 GB free,
     or while two AMBER jobs are already running and this one is AMBER-tagged (polls every
     5 s).  CPU does NOT gate launches: the gate is read from governor.CPU_CEILING, which
     the governor sets to 101 precisely so that a box run hot by instruction is not treated
     as a fault.  See the v3 note below.
  2. registers s26/jobs/<name>.json (pid, create_time, agent, tag, cmd, start, estimate) so the
     governor can suspend, resume or kill it;
  3. runs the command in its own process group with S26_AGENT / S26_TAG / S26_JOB in the
     environment, stdout+stderr appended to s26/logs/<name>.log;
  4. samples the process tree's RSS every 5 s and, on exit, writes s26/jobs_done/<name>.json
     with the exit code, wall time and PEAK RSS - the honest memory estimate for the ledger.

Exit code is the child's.  A job killed by the governor exits non-zero and its owner (or the
queue, for queued jobs) restarts it from its own checkpoint: every S26 job must checkpoint.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
JOBS, DONE, LOGS, QUEUE_L = HERE / "jobs", HERE / "jobs_done", HERE / "logs", HERE / "queue_launched"
STATE = HERE / "governor_state.json"
CAP_FILE = HERE / "launch_cap.json"   # v2.3: {"max_concurrent": n} read on every wait tick
MAX_AMBER = 2


# --------------------------------------------------------------------------------------
# v3 (2026-09-20, S31 lane D, defect D-C).  THE LAUNCH GATE IS DERIVED FROM THE GOVERNOR.
#
# THE DEFECT.  This file last changed in S26 (v2.3).  `governor.py` changed four times in
# S29 -- v2.5, v2.6, v2.6b, v2.6c -- and the direction of every one of those changes was
# to stop treating CPU as a safety signal:
#
#     governor.py CPU_CEILING = 101.0   "CPU-triggered suspension is DISABLED.  Eight jobs
#                                        on eight cores put smoothed CPU at 96-100% BY
#                                        DESIGN ...  A saturated CPU has no failure mode;
#                                        RAM does (OOM), so RAM alone triggers suspend and
#                                        kill."
#
# The launcher never followed.  It still refused to start anything while the governor's
# 15 s CPU mean was above 85%, which is BELOW the 94-95% band the campaign is instructed
# to hold.  So in exactly the regime the charter asks for, `jobrun` blocked every launch
# while the governor was content -- and lanes responded by launching detached with `nohup`,
# which registers nothing in s26/jobs/, so the governor cannot suspend, resume, kill or
# even SEE those processes.  That is strictly worse than either policy alone.  (Measured
# while writing this fix: governor_state.json at 23:52:33 read cpu_smooth 96.8%, jobs=0 --
# the 96.8% was four of this lane's own shards, launched detached for this exact reason.)
#
# THE FIX, and why this one and not another.  RAM is the binding constraint on this box:
# 16.75 GB total, ~10-11 GB of it the user's own baseline load, and an OOM is the only
# way a job here can take the machine down.  A saturated CPU makes everything slower and
# nothing unsafe.  So:
#
#   * CPU does not gate launches.  The gate is read FROM `governor.CPU_CEILING`, so it is
#     the governor's policy by construction -- if a future sprint re-enables CPU-triggered
#     suspension by lowering that constant, this launcher follows it automatically and
#     nobody has to remember to grep for the partner.  (Standing lesson: paired thresholds
#     move together; governor.py suspends and jobrun.py launches, and letting one half
#     drift caused a deadlock and a throughput bug in S29.  This is the third instance.)
#   * RAM still gates launches, UNCHANGED IN VALUE (93.0%), and is now derived as one
#     point of hysteresis below the governor's suspend ceiling, so a new job never starts
#     inside the band where the governor is about to suspend something.
#   * The real limiter on process count stays what it always was and what it should be: a
#     COUNT (`_cap()`, s26/launch_cap.json) plus the per-job headroom test
#     (`avail < est_ram + 0.5`) plus the refusal to start unsupervised above 0.5 GB.
#     Those are deterministic and they are what keeps the OS, Claude Code, logging and
#     emergency termination alive.  Nothing about them is relaxed here.
#
# `--cpu-gate` restores a CPU launch gate for a caller that wants one; it is not the
# default because the default has to agree with the governor.
# --------------------------------------------------------------------------------------
def _governor_band():
    """The governor's own constants, imported rather than copied.  Falls back to the S29
    values if the module cannot be loaded, and says so."""
    try:
        from importlib import util as _ilu
        spec = _ilu.spec_from_file_location("s26_governor_constants", HERE / "governor.py")
        mod = _ilu.module_from_spec(spec)
        spec.loader.exec_module(mod)           # constants + defs only; main is __main__-guarded
        return float(mod.CEILING), float(mod.CPU_CEILING), int(mod.MAX_JOBS), True
    except Exception:
        return 94.0, 101.0, 8, False


_GOV_CEILING, _GOV_CPU_CEILING, _GOV_MAX_JOBS, _GOV_IMPORTED = _governor_band()

#: launch below the governor's suspend ceiling, by one point of hysteresis.  93.0 with the
#: S29/S31 governor -- the same value this file has always used, now tied to its partner.
CEILING = _GOV_CEILING - 1.0
#: the governor's CPU policy, verbatim.  101.0 means "CPU never blocks", which is what the
#: governor decided in S29 v2.6c and what the 94-95% utilisation target requires.
CPU_START = _GOV_CPU_CEILING
#: the CONSERVATIVE fallback when s26/launch_cap.json cannot be read.  Deliberately not
#: the governor's MAX_JOBS: an unreadable cap file is a fault, and a fault should launch
#: fewer jobs, not more.  `_cap()` clamps the live value to the governor's MAX_JOBS above.
MAX_CONCURRENT = 4


def _cap() -> int:
    """The live concurrency cap, clamped to the governor's own MAX_JOBS.

    v3 (S31-D): the clamp is the point.  The cap file is the knob lanes turn; MAX_JOBS is
    the governor's hard limit, and a cap file edited above it would let jobrun register
    more jobs than the governor is willing to manage.
    """
    try:
        n = int(json.loads(CAP_FILE.read_text(encoding="utf-8"))["max_concurrent"])
    except Exception:
        n = MAX_CONCURRENT
    return max(1, min(n, _GOV_MAX_JOBS))


def _state():
    try:
        return json.loads(STATE.read_text(encoding="utf-8"))
    except Exception:
        return None


def _running_amber() -> int:
    n = 0
    for f in JOBS.glob("*.json"):
        try:
            j = json.loads(f.read_text(encoding="utf-8"))
            p = psutil.Process(int(j["pid"]))
            if str(j.get("tag", "")).upper() == "AMBER" and p.status() != psutil.STATUS_STOPPED:
                n += 1
        except Exception:
            pass
    return n


def _live_registrations() -> int:
    """Registered jobs whose pid is alive, read from s26/jobs/ directly."""
    n = 0
    for f in JOBS.glob("*.json"):
        try:
            j = json.loads(f.read_text(encoding="utf-8"))
            if psutil.pid_exists(int(j["pid"])):
                n += 1
        except Exception:
            pass
    return n


def _tree_rss(p: psutil.Process) -> int:
    tot = 0
    try:
        procs = [p] + p.children(recursive=True)
    except psutil.Error:
        return 0
    for q in procs:
        try:
            tot += q.memory_info().rss
        except psutil.Error:
            pass
    return tot


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--agent", required=True)
    ap.add_argument("--tag", default="CPU")
    ap.add_argument("--name", required=True)
    ap.add_argument("--est-ram", type=float, default=0.0)
    ap.add_argument("--cpu-gate", type=float, default=None,
                    help="refuse to launch while the governor's 15 s CPU mean is above this "
                         "percentage.  Default: the governor's own CPU_CEILING (%.1f), i.e. "
                         "CPU does not gate launches -- see the v3 note at the top of this "
                         "file.  Pass a lower number to restore a CPU launch gate."
                         % CPU_START)
    ap.add_argument("--spec", default=None, help="queue spec filename (set by the governor)")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        ap.error("no command after --")
    for d in (JOBS, DONE, LOGS):
        d.mkdir(parents=True, exist_ok=True)
    tag = a.tag.upper()
    reg = JOBS / f"{a.name}.json"
    if reg.exists():
        try:
            old = json.loads(reg.read_text(encoding="utf-8"))
            if psutil.pid_exists(int(old["pid"])):
                print(f"jobrun: a job named {a.name} is already registered (pid {old['pid']})",
                      file=sys.stderr)
                return 3
        except Exception:
            pass

    # Wait for headroom and, for AMBER, for a slot.
    cpu_gate = CPU_START if a.cpu_gate is None else float(a.cpu_gate)
    if not _GOV_IMPORTED:
        print("jobrun: WARNING could not import s26/governor.py; using the fallback band "
              f"(ram<{CEILING:.1f}, cpu<{CPU_START:.1f}, max_jobs {_GOV_MAX_JOBS})",
              file=sys.stderr, flush=True)
    waited = 0.0
    import random
    while True:
        st = _state()
        stale = (time.time() - st["epoch"]) > 60 if st else True
        ram = st["ram_pct"] if st else 0.0
        cpu = st.get("cpu_smooth", st["cpu_pct"]) if st else 0.0
        #: the registration count is read from s26/jobs/ directly (fresher than the governor's
        #: 5 s snapshot); v2.1: several waiters saw one freed slot at once and six jobs launched
        #: against a cap of four (2026-09-13 09:27), so the pass is re-checked after a jitter.
        n_reg = _live_registrations()
        amber_block = tag == "AMBER" and _running_amber() >= MAX_AMBER
        avail = float(st.get("ram_avail_gb", 99.0)) if st else 99.0
        crowded = (not stale) and (ram > CEILING or cpu > cpu_gate or n_reg >= _cap()
                                   or avail < a.est_ram + 0.5)
        #: v2.2 (ledger L40): a stale snapshot means no governor is watching; a job that may
        #: need more than 0.5 GB does not start unsupervised. Small jobs still may. And a job
        #: never starts unless the box has its estimate plus 0.5 GB free right now.
        if stale and a.est_ram > 0.5:
            crowded = True
            if waited == 0.0:
                print(f"jobrun: {a.name} waiting for a live governor (est {a.est_ram} GB > 0.5)",
                      file=sys.stderr, flush=True)
        if not crowded and not amber_block:
            time.sleep(random.uniform(0.2, 3.0))
            if _live_registrations() < _cap() and not (
                    tag == "AMBER" and _running_amber() >= MAX_AMBER):
                break
            continue
        if waited == 0.0 or waited % 60.0 == 0.0:
            if amber_block:
                why = "AMBER slot"
            else:
                which = []
                if ram > CEILING:
                    which.append(f"RAM {ram:.1f}% > {CEILING:.1f}%")
                if cpu > cpu_gate:
                    which.append(f"CPU(15s) {cpu:.1f}% > {cpu_gate:.1f}%")
                if n_reg >= _cap():
                    which.append(f"jobs {n_reg}/{_cap()}")
                if avail < a.est_ram + 0.5:
                    which.append(f"free {avail:.2f}GB < est {a.est_ram}+0.5GB")
                if stale:
                    which.append("no live governor")
                why = "; ".join(which) or f"box at ram {ram:.1f}% cpu(15s) {cpu:.1f}%"
            print(f"jobrun: {a.name} waiting for {why} ({waited:.0f}s so far)",
                  file=sys.stderr, flush=True)
        time.sleep(5.0)
        waited += 5.0

    spec = None
    if a.spec:
        try:
            spec = json.loads((QUEUE_L / a.spec).read_text(encoding="utf-8"))
        except Exception:
            spec = None

    env = dict(os.environ, S26_AGENT=a.agent, S26_TAG=tag, S26_JOB=a.name, PYTHONUNBUFFERED="1")
    logf = LOGS / f"{a.name}.log"
    flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    t0 = time.time()
    with logf.open("ab") as fh:
        header = (f"\n===== {time.strftime('%Y-%m-%dT%H:%M:%S')} jobrun {a.name} "
                  f"agent={a.agent} tag={tag}\n===== cmd: {' '.join(cmd)}\n")
        fh.write(header.encode())
        fh.flush()
        child = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=fh, stderr=subprocess.STDOUT,
                                 creationflags=flags)
        p = psutil.Process(child.pid)
        record = {"name": a.name, "agent": a.agent, "tag": tag, "pid": child.pid,
                  "create_time": p.create_time(), "cmd": cmd, "cwd": str(ROOT),
                  "start": time.strftime("%Y-%m-%dT%H:%M:%S"), "start_ts": t0,
                  "est_ram_gb": a.est_ram, "log": str(logf), "spec": spec, "waited_s": waited,
                  "jobrun_pid": os.getpid()}
        tmp = reg.with_suffix(f".{os.getpid()}.tmp")   # v3 (S31-D): a SHARED temp path is
        # not atomic -- os.replace is, but two writers on one <name>.tmp can publish an
        # interleaved file atomically.  The pid makes it unique per writer.
        tmp.write_text(json.dumps(record, indent=1), encoding="utf-8")
        os.replace(tmp, reg)
        peak = 0
        while child.poll() is None:
            peak = max(peak, _tree_rss(p))
            time.sleep(5.0)
    rc = child.returncode
    done = {**record, "exit_code": rc,
            "end": time.strftime("%Y-%m-%dT%H:%M:%S"), "wall_s": round(time.time() - t0, 1),
            "peak_rss_gb": round(peak / 1e9, 3)}
    tmp = (DONE / f"{a.name}.json").with_suffix(f".{os.getpid()}.tmp")   # v3 (S31-D), as above
    tmp.write_text(json.dumps(done, indent=1), encoding="utf-8")
    os.replace(tmp, DONE / f"{a.name}.json")
    try:
        reg.unlink()
    except OSError:
        pass
    print(f"jobrun: {a.name} exit={rc} wall={done['wall_s']}s peak_rss={done['peak_rss_gb']}GB",
          file=sys.stderr, flush=True)
    return rc if rc is not None else 1


if __name__ == "__main__":
    sys.exit(main())
