#!/usr/bin/env python
"""s26/enqueue.py - put a job on the governor's queue.

    python s26/enqueue.py --agent P --tag AMBER --name c3_relax --priority 20 --est-ram 2.5 \
        -- python s26/c3_relax.py --resume

Lower priority number launches first.  The governor launches queued jobs only while the box has
been under 88% for 60 s and no job is suspended; AMBER specs wait for a free AMBER slot.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUEUE = HERE / "queue"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--agent", required=True)
    ap.add_argument("--tag", default="CPU")
    ap.add_argument("--name", required=True)
    ap.add_argument("--priority", type=int, default=50)
    ap.add_argument("--est-ram", type=float, default=0.0)
    ap.add_argument("--note", default="")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        ap.error("no command after --")
    QUEUE.mkdir(parents=True, exist_ok=True)
    spec = {"name": a.name, "agent": a.agent, "tag": a.tag.upper(), "priority": a.priority,
            "est_ram_gb": a.est_ram, "cmd": cmd, "cwd": str(HERE.parent), "note": a.note,
            "enqueued": time.strftime("%Y-%m-%dT%H:%M:%S")}
    target = QUEUE / f"{a.priority:03d}_{a.name}.json"
    tmp = target.with_suffix(".tmp")
    tmp.write_text(json.dumps(spec, indent=1), encoding="utf-8")
    os.replace(tmp, target)
    print(f"enqueued {target.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
