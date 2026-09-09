"""s24/c_probe4.py -- LANE C EXPLORATORY PROBE 4.  THE COVERAGE CLAIM (coordinator's route (a)).

The coordinator's L3 says quality obtained VIA the distogram score causes alignment, and offers
route (a) as the escape: "your samples are better ON THE SCORE'S OWN TERMS than what the pool
already contains -- the corpus is missing conformations the score would rank highly if it ever
saw them."  That is a COVERAGE claim, it does not need decorrelation, and it is measurable now.

This measures it four ways, all at MATCHED COUNT (500 generated vs the shipped K=500 pool):

  1. score quantiles      does the generator produce chains the shipped functional scores BETTER
                          than any retrieved window?  `frac_better_than_pool_best` and the rank
                          of the pool's best inside the merged 1000.
  2. union take-up        of a merged 1000 scored by ONE functional, what share of the top-75 is
                          generated?  (0.5 = the corpus was not special.)
  3. ORACLE coverage      is the generated set's best structure closer to the native than the
                          pool's best?  This is the coverage claim in the space that matters, and
                          it is a MAXIMUM, so it is reported against a matched-count best-of-N.
  4. in-band coverage     count of samples within 1.5 A of the pool's own best -- the recall
                          statistic the instrument already uses (FAIL18 is defined by it).

If (1) and (2) are strongly positive while the ENDPOINT is null -- which c_probe2 already
measured -- then route (a) is refuted on its own terms: the corpus is not missing candidates the
score wants, and supplying them in quantity does not move the answer.

Exploratory.  ORACLE columns are post-hoc scoring of native-free decisions.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from core import project as pj           # noqa: E402
from s24 import c_probe as CP            # noqa: E402

TOPM = 75
NMATCH = 500
BAND = 1.5


def run(pdbs=None, nmatch=NMATCH):
    tg = I.targets()
    sel = tg if pdbs is None else [t for t in tg if t["pdb"] in set(pdbs)]
    arms = ("T0_helix", "T1_blind", "T2_restype", "T3_pool")
    rows = []
    for t in sel:
        pdb = t["pdb"]; n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n)
        idx = I.pool_idx(u); Wp = u["W"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, i, j)), float)
        rrp = u["rr"][idx]
        top = np.argsort(scp, kind="stable")[:TOPM]
        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]),
             "pool_score_best": float(scp.min()), "pool_score_p75": float(np.sort(scp)[TOPM - 1]),
             "pool_best_ORACLE": float(rrp.min()),
             "univ_best_ORACLE": float(u["rr"].min()),
             "incumbent": float(I.ca_rmsd(CP._avg(Wp[top]), nat))}
        rng0 = SD.stable_rng("c_probe4", pdb)
        for tag in arms:
            fn = {"T0_helix": lambda q: CP.s_helix(n, nmatch, q),
                  "T1_blind": lambda q: CP.s_blind(u, n, nmatch, q),
                  "T2_restype": lambda q: CP.s_restype(u, n, nmatch, q),
                  "T3_pool": lambda q: CP.s_pool(u, n, nmatch, q, idx[top])}[tag]
            ph, ps = fn(rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            sc = np.asarray(I.shipped_score(dg, I.pair_dists(CA, i, j)), float)
            rr = I.kabsch_rmsd_batch(CA, nat)                      # ORACLE
            merged = np.concatenate([scp, sc]); src = np.concatenate(
                [np.zeros(len(scp), int), np.ones(len(sc), int)])
            ou = np.argsort(merged, kind="stable")[:TOPM]
            r[tag] = {
                "score_best": float(sc.min()), "score_med": float(np.median(sc)),
                "frac_better_than_pool_best": float((sc < scp.min()).mean()),
                "frac_better_than_pool_p75": float((sc < np.sort(scp)[TOPM - 1]).mean()),
                "union_gen_frac_top75": float(src[ou].mean()),
                "gen_best_ORACLE": float(rr.min()),
                "gen_p01_ORACLE": float(np.percentile(rr, 1)),
                "n_in_band_ORACLE": int((rr <= rrp.min() + BAND).sum()),
                "n_beats_pool_best_ORACLE": int((rr < rrp.min()).sum()),
            }
        rows.append(r)
        print("  %s inc %.3f poolbest %.3f | %s" % (pdb, r["incumbent"], r["pool_best_ORACLE"],
              " ".join("%s gb%.2f f%.2f u%.2f" % (k[:2], r[k]["gen_best_ORACLE"],
                       r[k]["frac_better_than_pool_p75"], r[k]["union_gen_frac_top75"])
                       for k in arms)), flush=True)
    p = os.path.join(RES, "c_probe4.json"); tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "nmatch": nmatch, "EXPLORATORY": True}, fh)
    os.replace(tmp, p)
    rep(rows, arms)
    return rows


def rep(rows, arms=("T0_helix", "T1_blind", "T2_restype", "T3_pool")):
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    print("\nEXPLORATORY n=%d.  Matched count: %d generated vs the shipped K=500 pool." % (
        len(rows), NMATCH))
    print("  pool ORACLE best %.4f   whole-universe ORACLE best %.4f   incumbent %.4f"
          % (g("pool_best_ORACLE").mean(), g("univ_best_ORACLE").mean(), g("incumbent").mean()))
    print("\n  %-14s %9s %9s %9s | %9s %9s %9s %9s" % (
        "arm", "f>poolbst", "f>poolp75", "union%gen",
        "genbest", "gen p01", "#in-band", "#<poolbst"))
    for k in arms:
        f = lambda q: np.array([r[k][q] for r in rows], float)    # noqa: E731
        print("  %-14s %9.4f %9.4f %9.4f | %9.4f %9.4f %9.2f %9.2f"
              % (k, f("frac_better_than_pool_best").mean(), f("frac_better_than_pool_p75").mean(),
                 f("union_gen_frac_top75").mean(), f("gen_best_ORACLE").mean(),
                 f("gen_p01_ORACLE").mean(), f("n_in_band_ORACLE").mean(),
                 f("n_beats_pool_best_ORACLE").mean()))
    print("\n  ROUTE (a) needs BOTH: the score preferring generated chains (left block) AND those")
    print("  chains being structurally better (genbest < pool ORACLE best).  Read the two together.")


if __name__ == "__main__":
    tg = I.targets()
    run([t["pdb"] for t in tg[::5]])
