#!/usr/bin/env python
"""s26/jobrun.py - launch one job under the S26 governor.

    python s26/jobrun.py --agent Q --tag CPU --name a2_dla [--est-ram 1.0] -- python s26/q_dla.py ...

What it does, in order:
  1. refuses to start while the governor's last snapshot shows the box above the ceiling,
     or while two AMBER jobs are already running and this one is AMBER-tagged (polls every 5 s);
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
CEILING, MAX_AMBER = 93.0, 2
#: v2 (2026-09-13): a new job also waits while the governor's 15 s CPU mean is above
#: CPU_START, or while MAX_CONCURRENT jobs are already registered. Five lanes launching at
#: once put the 8-core box at 98.7% CPU on 2026-09-13 00:37 and the governor thrashed.
CPU_START, MAX_CONCURRENT = 85.0, 4


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
        crowded = (not stale) and (ram > CEILING or cpu > CPU_START or n_reg >= MAX_CONCURRENT)
        if not crowded and not amber_block:
            time.sleep(random.uniform(0.2, 3.0))
            if _live_registrations() < MAX_CONCURRENT and not (
                    tag == "AMBER" and _running_amber() >= MAX_AMBER):
                break
            continue
        if waited == 0.0 or waited % 60.0 == 0.0:
            why = ("AMBER slot" if amber_block else
                   f"box at ram {ram:.1f}% cpu(15s) {cpu:.1f}% jobs {n_reg}/{MAX_CONCURRENT}")
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
        tmp = reg.with_suffix(".tmp")
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
    tmp = (DONE / f"{a.name}.json").with_suffix(".tmp")
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
