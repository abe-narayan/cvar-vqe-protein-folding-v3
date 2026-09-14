#!/usr/bin/env python
"""s26/i_tree_watchdog.py -- keep a lane's job trees consistent with the governor's view.

Found 02:26 to 02:54 on 2026-09-14 (ledger, lane I): the governor decides "suspended" from the
ROOT process's status.  A root that is blocked in `subprocess.run` waiting for its child reads
"running" to psutil on Windows even after `suspend()`, so the governor re-suspended the whole
tree on every hot sample (18 times, 0 resumes) and never resumed it; on Windows every
`suspend()` increments a per-thread suspend count, so the grandchild needed 16 `resume()` calls
before it ran again.  The job sat for 28 minutes with 13 CPU-seconds done.

This watchdog acts only on jobs of ONE agent (`--agent`, default I) and only in the state the
governor cannot see: when a job's root is RUNNING but a descendant is STOPPED, it resumes the
descendant (repeatedly, until it runs).  It never touches a job whose root is stopped, which is
the governor's deliberate suspension, and never touches other agents' jobs.

    python s26/i_tree_watchdog.py [--agent I] [--interval 10]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import time

import psutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOBS = os.path.join(ROOT, "s26", "jobs")
LOG = os.path.join(ROOT, "s26", "logs", "i_tree_watchdog.log")


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def tick(agent):
    for f in sorted(glob.glob(os.path.join(JOBS, "*.json"))):
        try:
            job = json.loads(open(f, encoding="utf-8").read())
        except Exception:                                                # noqa: BLE001
            continue
        if str(job.get("agent")) != agent:
            continue
        try:
            root = psutil.Process(int(job["pid"]))
            if abs(root.create_time() - float(job.get("create_time", root.create_time()))) > 2.0:
                continue
            if root.status() == psutil.STATUS_STOPPED:
                continue                     # the governor's own suspension: respect it
            for p in root.children(recursive=True):
                try:
                    n = 0
                    while p.status() == psutil.STATUS_STOPPED and n < 200:
                        p.resume()
                        n += 1
                        time.sleep(0.02)
                    if n:
                        log(f"{job['name']}: root {root.pid} running, descendant {p.pid} was stopped; "
                            f"{n} resume call(s) -> {p.status()}")
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--agent", default="I")
    ap.add_argument("--interval", type=float, default=10.0)
    a = ap.parse_args(argv)
    log(f"watchdog start agent={a.agent} interval={a.interval}s")
    while True:
        tick(a.agent)
        time.sleep(a.interval)


if __name__ == "__main__":
    main()
