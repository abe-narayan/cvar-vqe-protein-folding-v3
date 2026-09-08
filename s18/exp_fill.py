"""s18/exp_fill.py -- a coefficient-cache filler for MATH's `Target`, nothing else.

The 126-target experiment's cost is entirely MATH's conditional-moment fit, which is cached on
disk by `(pdb, mu, S, grid)`.  MATH's own `build` populated the first ~51 targets and then
exited; this fills the rest so `s18/exp_run.py main` can consume them.

It computes NOTHING scientific: it calls MATH's `Target.fit` with MATH's declared parameters and
writes the same `.npz` MATH's `build` writes, with MATH's config hash inside it.  A `stride` and
`offset` let several fillers share the list without recomputing each other's targets, and the
main run's own dynamic ordering takes uncached targets from the far end, so the workers converge
instead of colliding.  Every fit is deterministic in `stable_rng(pdb, mu, S, grid, 'mufreeze')`,
so which worker computed a file cannot change its contents.

    python -m s18.exp_fill [offset] [stride]
"""
from __future__ import annotations

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                # noqa: E402
from s18 import math_anova as MA               # noqa: E402
from s18 import math_lib as ML                 # noqa: E402


def run(offset=0, stride=1, mu="pool", S=None, grid=None):
    S = S or MA.NSAMP
    grid = grid or MA.GRID
    tg = I.targets()
    deb = MA.debias_map(tg)
    todo = tg[offset::stride]
    t0 = time.time()
    made = 0
    for t in todo:
        pdb = t["pdb"]
        path = os.path.join(ML.CACHE, f"anova_{pdb}_{mu}_{S}_{grid}.npz")
        if os.path.exists(path):
            continue
        d = MA.gather_one(t, deb[int(t["fold"])])
        #: write to a worker-private temp name and rename, so a concurrent reader can never
        #: see a half-written coefficient file and mistake it for a complete one
        tmp = path + f".w{offset}.tmp.npz"
        #: batch: see the note in `s18/exp_run.py`.  Chunk size only; coefficients identical.
        ob = MA.Target(d["pdb"], d["seq"], d["n"], d["fold"],
                       d["dhat"], d["sd"], d["i"], d["j"]).fit(mu, S, grid, batch=4096)
        ob.save(tmp)
        try:
            os.replace(tmp, path)
        except OSError:
            pass
        made += 1
        print(f"  filled {pdb} ({made}, {time.time()-t0:.0f}s)", flush=True)
    print(f"DONE fill offset={offset} stride={stride}: {made} new in "
          f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 0,
        int(sys.argv[2]) if len(sys.argv) > 2 else 1)
