#!/usr/bin/env python
"""s26/governor.py - the Sprint 26 resource governor.

Keeps this 16.75 GB / 8-core box inside the band the S26 campaign requires:

    * never above 93% RAM or CPU      -> suspend the newest registered job
    * above 95% for more than 15 s    -> kill the newest job (CTRL_BREAK first, so it can
                                         checkpoint; terminate the tree after a grace period)
                                         and put its spec back on the queue
    * below 88% for more than 60 s    -> launch the next queued job, if any
    * never more than two AMBER jobs  -> a third AMBER-tagged job is suspended on sight
                                         until a slot frees, and queued AMBER specs wait

Every action is appended to s26/governor.log with a timestamp.  A snapshot is written to
s26/governor_state.json every sample so agents (and s26/jobrun.py) can read the live band.

A "job" is any process registered in s26/jobs/<name>.json.  Two ways to register one:

    python s26/jobrun.py --agent Q --tag CPU --name a2_dla -- python s26/q_dla.py     (direct)
    python s26/enqueue.py --agent Q --tag AMBER --name c3_relax -- python s26/c3.py   (queued)

Tags: CPU (default), AMBER (OpenMM; capped at 2), ESM (loads the ESM bank), TEST (pytest).

Run it:

    python s26/governor.py            # foreground, Ctrl-C stops it
    python s26/governor.py --once     # one sample, print it, exit
    python s26/governor.py --status   # print the last snapshot and the registered jobs

Lanes (Claude agents) are not processes this script can start; the coordinator records
them in s26/lanes.json and the governor logs a warning when fewer than 4 or more than 8
are recorded, so the log carries the whole picture.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
JOBS = HERE / "jobs"
DONE = HERE / "jobs_done"
QUEUE = HERE / "queue"
LAUNCHED = HERE / "queue_launched"
LOG = HERE / "governor.log"
STATE = HERE / "governor_state.json"
LANES = HERE / "lanes.json"

CEILING = 94.0        # % RAM (v2.5, S29: the user asked for 94-95%), or % CPU smoothed over CPU_WINDOW samples: suspend the newest job
HARD = 95.5           # % RAM sustained HARD_SECONDS (v2.5): kill the newest job (RAM only, v2)
HARD_SECONDS = 15.0
RESUME_BELOW = 92.0   # RAM band is 92-94 (v2.5); resume only once RAM is back under 90 ...
CPU_CEILING = 99.0    # v2.6: CPU gets its own ceiling; the box is meant to run hot
CPU_RESUME = 99.0     # v2.6 (S29): the box is deliberately run at 94-95% CPU on the user's
                      # instruction, so a resume gated on CPU < 80 could NEVER fire and starved
                      # five jobs across four lanes for 40 min. RAM is the real constraint here.
CPU_WINDOW = 3        # samples in the CPU rolling mean (3 x 5 s = 15 s)
MIN_SUSPEND = 20.0    # seconds a suspended job stays suspended before it may resume (v2)
STALL_SECONDS = 180.0 # RAM parked in [RESUME_BELOW, CEILING) with jobs suspended: kill the fattest (v2.2)
LOW = 91.0            # below this for LOW_SECONDS (v2.5): launch the next queued job
LOW_SECONDS = 60.0
SAMPLE = 5.0
MAX_AMBER = 2
MAX_JOBS = 8
MIN_LANES, MAX_LANES = 4, 8
KILL_GRACE = 25.0     # seconds between CTRL_BREAK and terminate

# v2 (2026-09-13 08:50). The 00:37 log showed the v1 band oscillating: with five jobs the raw
# CPU sample crossed 93% and fell under 90% on alternate 5 s ticks, so the newest job was
# suspended and resumed every 5 s, and a sustained CPU spike would have KILLED a job that was
# doing nothing wrong (CPU overload slows the box; it cannot crash it). v2: CPU is a 15 s
# rolling mean; CPU above the ceiling suspends but never kills; only RAM kills; a suspended
# job stays down at least MIN_SUSPEND seconds; resumption needs RAM < 90 AND CPU < 80.


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def log(kind: str, msg: str) -> None:
    line = f"{now()} {kind:8s} {msg}"
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def read_json(p: Path):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def write_json_atomic(p: Path, obj) -> None:
    """Atomic replace, retried: on Windows `os.replace` raises PermissionError while another
    process (a jobrun waiter, a lane's `--status`) holds the target open for reading. That
    killed governor v2 once (19:36, ledger L59); a snapshot lost to a reader is skipped, not fatal."""
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1), encoding="utf-8")
    for attempt in range(8):
        try:
            os.replace(tmp, p)
            return
        except PermissionError:
            time.sleep(0.05 * (attempt + 1))
    log("WARN", f"could not replace {p.name} after 8 attempts (a reader holds it); snapshot skipped")


# ----------------------------------------------------------------------------- jobs

def _proc(job: dict):
    """The psutil.Process for a job, or None if its pid is gone or reused."""
    try:
        p = psutil.Process(int(job["pid"]))
        ct = job.get("create_time")
        if ct is not None and abs(p.create_time() - float(ct)) > 2.0:
            return None
        return p
    except (psutil.NoSuchProcess, psutil.AccessDenied, ValueError, KeyError):
        return None


def tree(p: psutil.Process) -> list:
    try:
        return [p] + p.children(recursive=True)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return []


def tree_rss(p: psutil.Process) -> int:
    total = 0
    for q in tree(p):
        try:
            total += q.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total


def is_suspended(p: psutil.Process) -> bool:
    try:
        return p.status() == psutil.STATUS_STOPPED
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return False


def load_jobs() -> list:
    """Registered jobs whose process is alive.  Dead registrations are reaped."""
    jobs = []
    for f in sorted(JOBS.glob("*.json")):
        job = read_json(f)
        if not job:
            continue
        p = _proc(job)
        if p is None:
            log("REAP", f"{job.get('name')} (pid {job.get('pid')}) is gone; registration removed")
            try:
                done = DONE / f.name
                if not done.exists():
                    job["reaped_by_governor"] = now()
                    write_json_atomic(done, job)
                f.unlink()
            except OSError:
                pass
            continue
        job["_proc"] = p
        job["_file"] = f
        job["_rss"] = tree_rss(p)
        job["_suspended"] = (job.get("name") in SUSPENDED) or is_suspended(p)
        jobs.append(job)
    jobs.sort(key=lambda j: float(j.get("start_ts", 0.0)))   # oldest first
    return jobs


SUSPENDED: set = set()     # v2.4 (ledger L129): the governor's own memory of what it suspended.
                           # Windows counts suspends per thread, and a root blocked in a wait
                           # reads "running" after suspend(), so psutil's status is not the
                           # truth; this set is. Every SUSPEND is matched by exactly one RESUME.


def suspend(job: dict) -> None:
    if job["name"] in SUSPENDED:
        return
    for q in tree(job["_proc"]):
        try:
            q.suspend()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    SUSPENDED.add(job["name"])
    log("SUSPEND", f"{job['name']} agent={job.get('agent')} tag={job.get('tag')} "
                   f"rss={job['_rss'] / 1e9:.2f}GB")


def resume(job: dict) -> None:
    for q in reversed(tree(job["_proc"])):
        try:
            q.resume()
            #: a process suspended more than once (by an earlier governor) needs as many
            #: resumes; keep resuming while psutil still reports it stopped, bounded.
            for _ in range(32):
                if not is_suspended(q):
                    break
                q.resume()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    SUSPENDED.discard(job["name"])
    log("RESUME", f"{job['name']} agent={job.get('agent')} tag={job.get('tag')}")


def kill(job: dict, reason: str) -> None:
    p = job["_proc"]
    log("KILL", f"{job['name']} agent={job.get('agent')} tag={job.get('tag')} "
                f"rss={job['_rss'] / 1e9:.2f}GB reason={reason}; sending CTRL_BREAK, "
                f"grace {KILL_GRACE:.0f}s for a checkpoint")
    try:
        if job.get("_suspended"):
            resume(job)                      # a stopped process cannot checkpoint
        if os.name == "nt":
            os.kill(p.pid, signal.CTRL_BREAK_EVENT)
        else:
            p.send_signal(signal.SIGINT)
    except Exception as e:                    # noqa: BLE001
        log("WARN", f"signal to {job['name']} failed: {e!r}")
    deadline = time.time() + KILL_GRACE
    while time.time() < deadline and p.is_running():
        time.sleep(1.0)
    for q in reversed(tree(p)):
        try:
            q.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    log("KILLED", f"{job['name']} terminated")
    spec = job.get("spec")
    if spec:
        spec = dict(spec)
        spec["requeued_at"] = now()
        spec["requeue_reason"] = reason
        spec["requeue_count"] = int(spec.get("requeue_count", 0)) + 1
        target = QUEUE / f"{int(spec.get('priority', 50)):03d}_{spec['name']}.json"
        write_json_atomic(target, spec)
        log("REQUEUE", f"{job['name']} -> {target.name} (count {spec['requeue_count']})")
    else:
        log("WARN", f"{job['name']} was a direct job with no spec; its owner must restart "
                    f"it from its checkpoint")
    try:
        job["_file"].unlink()
    except OSError:
        pass


def n_amber(jobs: list) -> int:
    return sum(1 for j in jobs if str(j.get("tag", "")).upper() == "AMBER" and not j["_suspended"])


def launch_next(jobs: list) -> bool:
    specs = sorted(QUEUE.glob("*.json"))
    amber_running = n_amber(jobs)
    for f in specs:
        spec = read_json(f)
        if not spec:
            continue
        if str(spec.get("tag", "CPU")).upper() == "AMBER" and amber_running >= MAX_AMBER:
            continue
        cmd = [sys.executable, str(HERE / "jobrun.py"),
               "--agent", str(spec.get("agent", "?")),
               "--tag", str(spec.get("tag", "CPU")),
               "--name", str(spec["name"]),
               "--est-ram", str(spec.get("est_ram_gb", 0.0)),
               "--spec", str(f.name), "--"] + list(spec["cmd"])
        cwd = spec.get("cwd") or str(ROOT)
        flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
        subprocess.Popen(cmd, cwd=cwd, creationflags=flags,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.replace(f, LAUNCHED / f.name)
        shown = " ".join(map(str, spec["cmd"]))[:160]
        log("LAUNCH", f"{spec['name']} agent={spec.get('agent')} tag={spec.get('tag')} "
                      f"est_ram={spec.get('est_ram_gb', '?')}GB cmd={shown}")
        return True
    return False


# ----------------------------------------------------------------------------- loop

_CPU_HIST: list = []


def sample(jobs: list) -> dict:
    vm = psutil.virtual_memory()
    cpu = psutil.cpu_percent(interval=None)
    _CPU_HIST.append(cpu)
    del _CPU_HIST[:-CPU_WINDOW]
    cpu_smooth = sum(_CPU_HIST) / len(_CPU_HIST)
    lanes = read_json(LANES) or {}
    return {
        "ts": now(), "epoch": time.time(),
        "ram_pct": vm.percent, "ram_used_gb": round(vm.used / 1e9, 3),
        "ram_avail_gb": round(vm.available / 1e9, 3), "ram_total_gb": round(vm.total / 1e9, 3),
        "cpu_pct": cpu, "cpu_smooth": round(cpu_smooth, 1),
        "jobs": [{"name": j["name"], "agent": j.get("agent"), "tag": j.get("tag"),
                  "pid": j["pid"], "rss_gb": round(j["_rss"] / 1e9, 3),
                  "suspended": j["_suspended"], "start": j.get("start")} for j in jobs],
        "n_jobs": len(jobs), "n_amber": n_amber(jobs),
        "n_queued": len(list(QUEUE.glob("*.json"))),
        "lanes": lanes.get("active", []), "n_lanes": len(lanes.get("active", [])),
        "band": {"ceiling": CEILING, "hard": HARD, "resume_below": RESUME_BELOW,
                 "cpu_resume": CPU_RESUME, "low": LOW, "version": 2},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()
    for d in (JOBS, DONE, QUEUE, LAUNCHED):
        d.mkdir(parents=True, exist_ok=True)

    if args.status:
        print(json.dumps(read_json(STATE), indent=1))
        return

    psutil.cpu_percent(interval=None)
    if args.once:
        time.sleep(1.0)
        jobs = load_jobs()
        s = sample(jobs)
        print(json.dumps(s, indent=1))
        return

    log("START", f"governor v2 pid {os.getpid()} ceiling={CEILING}% hard(ram)={HARD}%/{HARD_SECONDS:.0f}s "
                 f"cpu_resume={CPU_RESUME}% cpu_window={CPU_WINDOW} min_suspend={MIN_SUSPEND:.0f}s "
                 f"low={LOW}%/{LOW_SECONDS:.0f}s sample={SAMPLE:.0f}s max_amber={MAX_AMBER}")
    breach_since = None
    hard_since = None
    low_since = None
    suspended_stack: list = []      # names, in the order we suspended them (LIFO resume)
    suspended_at: dict = {}         # name -> time of suspension (MIN_SUSPEND)
    stall_since = None              # v2.2 stall breaker
    last_sample_log = 0.0
    last_lane_warn = 0.0
    while True:
        t0 = time.time()
        try:
            jobs = load_jobs()
            s = sample(jobs)
            write_json_atomic(STATE, s)
        except Exception as exc:          # a transient file race must not stop supervision
            log("WARN", f"sample failed ({type(exc).__name__}: {exc}); retrying next tick")
            time.sleep(SAMPLE)
            continue
        ram, cpu = s["ram_pct"], s["cpu_smooth"]
        # v2.6 (S29): RAM and CPU get SEPARATE ceilings. The user asked for the box to run at
        # 94-95% CPU, so a single CEILING shared with RAM made the governor suspend a job every
        # time the instruction was being followed, then resume it seconds later (thrash).
        # CPU is a throughput signal here, not a failure mode; RAM is the real constraint.
        hot = ram if ram >= CEILING else (cpu if cpu > CPU_CEILING else 0.0)
        which = "ram" if ram >= CEILING else "cpu"

        if t0 - last_sample_log >= 60.0:
            log("SAMPLE", f"ram={ram:.1f}% ({s['ram_used_gb']:.2f}GB used, {s['ram_avail_gb']:.2f}GB free) "
                          f"cpu={s['cpu_pct']:.1f}% (15s mean {cpu:.1f}%) jobs={s['n_jobs']} "
                          f"amber={s['n_amber']} queued={s['n_queued']} lanes={s['n_lanes']}")
            last_sample_log = t0

        if (s["n_lanes"] < MIN_LANES or s["n_lanes"] > MAX_LANES) and t0 - last_lane_warn >= 600:
            log("WARN", f"{s['n_lanes']} lanes recorded in lanes.json; the campaign requires "
                        f"{MIN_LANES}-{MAX_LANES}")
            last_lane_warn = t0

        # Adopt suspensions this governor did not make (v2.3, ledger L83): after a restart the
        # stack is empty, and a job the previous governor suspended would otherwise never
        # be resumed. Any stopped job not on the stack is put on it now.
        for j in jobs:
            if j["_suspended"] and j["name"] not in suspended_stack:
                suspended_stack.append(j["name"])
                suspended_at.setdefault(j["name"], 0.0)
                log("ADOPT", f"{j['name']} was suspended before this governor started; adopted")

        # AMBER cap on direct-launched jobs: suspend any beyond the second, newest first.
        amber_live = [j for j in jobs if str(j.get("tag", "")).upper() == "AMBER" and not j["_suspended"]]
        while len(amber_live) > MAX_AMBER:
            j = amber_live.pop()             # newest
            suspend(j)
            j["_suspended"] = True
            suspended_stack.append(j["name"])
            log("WARN", f"AMBER cap: {j['name']} suspended until a slot frees")

        # Hard breach (RAM only, v2): > HARD for HARD_SECONDS -> kill the newest job.
        if ram > HARD:
            hard_since = hard_since or t0
            if t0 - hard_since > HARD_SECONDS and jobs:
                victim = jobs[-1]
                kill(victim, f"ram {ram:.1f}% > {HARD}% for {t0 - hard_since:.0f}s")
                if victim["name"] in suspended_stack:
                    suspended_stack.remove(victim["name"])
                suspended_at.pop(victim["name"], None)
                hard_since = None
        else:
            hard_since = None

        # Soft breach: > CEILING (RAM, or smoothed CPU) -> suspend the newest running job.
        if hot > CEILING:
            breach_since = breach_since or t0
            running = [j for j in jobs if not j["_suspended"]]
            if running:
                j = running[-1]
                suspend(j)
                suspended_stack.append(j["name"])
                suspended_at[j["name"]] = t0
            elif not jobs:
                if t0 - breach_since > 60 and int(t0 - breach_since) % 60 < SAMPLE:
                    log("WARN", f"{which} at {hot:.1f}% with no registered jobs: the pressure "
                                f"is outside the campaign's processes")
        else:
            breach_since = None

        # Band restored (RAM under 90 AND smoothed CPU under 80): resume, one per sample, the
        # most recently suspended job first, after it has been down at least MIN_SUSPEND
        # seconds, unless it is an AMBER job with no free slot.
        if ram < RESUME_BELOW and cpu < CPU_RESUME and suspended_stack:
            by_name = {j["name"]: j for j in jobs}
            for name in reversed(list(suspended_stack)):
                j = by_name.get(name)
                if j is None or not j["_suspended"]:
                    suspended_stack.remove(name)
                    suspended_at.pop(name, None)
                    continue
                if t0 - suspended_at.get(name, 0.0) < MIN_SUSPEND:
                    continue
                if str(j.get("tag", "")).upper() == "AMBER" and n_amber(jobs) >= MAX_AMBER:
                    continue
                resume(j)
                suspended_stack.remove(name)
                suspended_at.pop(name, None)
                break

        # Stall breaker (v2.2, ledger L74): RAM parked between RESUME_BELOW and CEILING with
        # jobs suspended means the pressure is the user's own load and nothing of ours will
        # ever resume. After STALL_SECONDS in that state, kill the suspended job holding the
        # most memory (CTRL_BREAK first, so it checkpoints); its owner relaunches it later.
        if ram >= RESUME_BELOW and suspended_stack and [j for j in jobs if j["_suspended"]]:
            stall_since = stall_since or t0
            if t0 - stall_since > STALL_SECONDS:
                victim = max((j for j in jobs if j["_suspended"]), key=lambda j: j["_rss"])
                kill(victim, f"stall: ram {ram:.1f}% >= {RESUME_BELOW}% for {t0 - stall_since:.0f}s "
                             f"with {len(suspended_stack)} suspended; freeing {victim['_rss'] / 1e9:.2f}GB")
                if victim["name"] in suspended_stack:
                    suspended_stack.remove(victim["name"])
                suspended_at.pop(victim["name"], None)
                stall_since = None
        else:
            stall_since = None

        # Low water: < LOW for LOW_SECONDS -> launch the next queued job.
        if hot < LOW and not suspended_stack:
            low_since = low_since or t0
            if t0 - low_since > LOW_SECONDS and len(jobs) < MAX_JOBS:
                if launch_next(jobs):
                    low_since = None
        else:
            low_since = None

        time.sleep(max(0.5, SAMPLE - (time.time() - t0)))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("STOP", "governor stopped by Ctrl-C")
