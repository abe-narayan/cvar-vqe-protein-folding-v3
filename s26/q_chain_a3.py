"""s26/q_chain_a3.py -- launch the two A3 build shards through jobrun once both A1 shards are done.

The coordinator's ordering (L33 message): A3 runs the same way as A1, after A1's shards finish.
This waits for s26/jobs_done/a1_build_s0.json and a1_build_s1.json, then starts the two A3
shards through s26/jobrun.py (each with the thread variables at 1, est-ram 0.4, per-target
checkpoints under s26/results/a3/), and waits for them.  It reads nothing but job records.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PY = sys.executable
DONE = os.path.join(HERE, "jobs_done")
WAIT_FOR = ("a1_build_s0.json", "a1_build_s1.json")

A3_ARGS = ["--build", "--tag", "a3", "--variants", "zrank,zraw,asinh,soft",
           "--adapt-variants", "zrank,zraw", "--pools", "L2", "--optimisers", "adam_best"]


def main():
    t0 = time.time()
    while not all(os.path.exists(os.path.join(DONE, f)) for f in WAIT_FOR):
        time.sleep(60)
    print(f"A1 shards done after {time.time() - t0:.0f} s of waiting; launching A3", flush=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    procs = []
    for i in range(2):
        cmd = [PY, os.path.join(HERE, "jobrun.py"), "--agent", "Q", "--tag", "CPU",
               "--name", f"a3_build_s{i}", "--est-ram", "0.4", "--",
               PY, os.path.join(HERE, "q_adapt.py")] + A3_ARGS + ["--shard", f"{i}/2"]
        procs.append(subprocess.Popen(cmd, cwd=ROOT, env=env))
        time.sleep(5)
    rcs = [p.wait() for p in procs]
    print(f"A3 shards exited with {rcs}", flush=True)
    return max(rcs)


if __name__ == "__main__":
    sys.exit(main())
