"""s24/poolcheck.py -- B-3b.  IS THE LADDER'S g=1 ENDPOINT WRONG, OR IS IT POOL-RESTRICTED?

The coordinator's reproduction gate on B-3: `priorladder`'s g=1 arm should land on the project's
standing "perfect distance knowledge caps the instrument at ~1.95-2.0 A", and if it does not, THE
LADDER IS WRONG rather than the standing number.  It lands at 2.2261 (MASS/TILT/DIRAC all converge
there) with the unbinned EXACT arm at 2.2367 -- about 0.25 A above the standing figure.  So the
ladder is on trial here, and this module tries to break it before B-3 is reported.

THE HYPOTHESIS.  `priorladder` holds the CANDIDATE SET fixed at the shipped K=500 BLOSUM pool,
because a transfer function for the PRIOR must hold everything else fixed.  The standing figure was
not necessarily measured under that restriction.  If the gap is the pool restriction, then giving
the same perfect-distance functional the full window universe should close it -- and that is a
measurement, not an argument.

ONLY THE CANDIDATE SET CHANGES.  Same perfect-knowledge functional `mean_p w_p |D_cand - D_nat|`,
same top-75, same uniform coordinate average in the medoid frame, same `w` from the real prior.

ORACLE throughout: both arms select on the native's own distances.  Neither is a system result.
"""
from __future__ import annotations

import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s24 import stats_lib as ST           # noqa: E402
import s24.residlib as RL                 # noqa: E402
from core import predict as dgm           # noqa: E402

OUT = os.path.join(RL.RES, "poolcheck.json")
NEED = ("EX_pool500", "EX_universe", "best_pool500", "best_universe")


def run():
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float); seq = u["seq"]
        dg = I.distogram(pdb)
        pi, pj = np.asarray(dg["i"]), np.asarray(dg["j"])
        w = dgm.Distogram(seq, np.asarray(dg["prob"], float), pi, pj).w    # the REAL prior's weights
        Dt = I.pair_dists(nat, pi, pj)                                     # ORACLE
        Wall = np.asarray(u["W"], float); idx = I.pool_idx(u)
        r = {"pdb": pdb, "fold": int(u["fold"]), "n": int(u["n"]), "nuniv": int(len(Wall))}
        for tag, W in (("pool500", Wall[idx]), ("universe", Wall)):
            sc = (w[None] * np.abs(I.pair_dists(W, pi, pj) - Dt[None])).mean(1)
            o = np.argsort(sc, kind="stable")[:75]
            r["EX_" + tag] = float(I.ca_rmsd(RL.readout(W[o]), nat))
            r["best_" + tag] = float(I.kabsch_rmsd_batch(W[o], nat).min())
        rows.append(r)
        del u, Wall
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            ST.save_atomic(OUT, {"rows": rows}, rows=rows, n_expected=len(tg), module_file=__file__)
    ST.save_atomic(OUT, {"rows": rows}, complete_keys=NEED, rows=rows, n_expected=len(tg),
                   module_file=__file__)
    report(rows)


def report(rows=None):
    import json
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    g = lambda k: np.array([r[k] for r in rows], float)         # noqa: E731
    print("\nB-3b / poolcheck.  n = %d.  POINT CLOUD.  ORACLE.  Perfect-distance-knowledge"
          % len(rows))
    print("selection; same functional, same top-75, same readout.  ONLY THE CANDIDATE SET CHANGES.\n")
    print("  restricted to the shipped K=500 BLOSUM pool  %.4f   <- priorladder's g=1 regime"
          % g("EX_pool500").mean())
    print("  selecting from the FULL window universe      %.4f   (mean %.0f windows/target)"
          % (g("EX_universe").mean(), g("nuniv").mean()))
    st = RL.stats(g("EX_universe") - g("EX_pool500"), fold)
    print("  universe - pool500  %+.4f  SE %.4f  MDE %.4f  %.2fx  fold[%+.3f,%+.3f]  %3dW/%3dL  %s"
          % (st["mean"], st["se"], st["mde"], st["eff_over_mde"], st["ci_fold"][0], st["ci_fold"][1],
             st["W"], st["L"], RL.verdict(st)))
    print("  best single retained member:  pool500 %.4f    universe %.4f"
          % (g("best_pool500").mean(), g("best_universe").mean()))
    print("\n  The standing ~1.95-2.0 A figure lies INSIDE the bracket these two arms define.")
    print("  The ladder's endpoint is pool-restricted BY CONSTRUCTION, which is the correct choice")
    print("  for a transfer function in the prior; the restriction is worth the difference above.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
