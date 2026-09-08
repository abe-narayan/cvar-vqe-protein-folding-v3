"""s17/phys_gate.py -- THE GATE'S RECALL/BEST PARADOX, resolved.

E3 found something that looks self-contradictory and is not:

    the Legacy gate keeps MORE sub-2 A candidates than a matched random gate of the same
    count, and simultaneously the set BEST after gating is SIGNIFICANTLY WORSE.

Two reasons that could be an artefact rather than a finding, both tested here:

  (a) DIFFERENT DENOMINATORS.  Near-native recall is only defined on targets that have a
      sub-2 A candidate at all; set best is defined on all 126.  So the two statistics are
      not measured on the same targets.  Everything below is recomputed on the SUB-2A
      STRATUM only, so both statistics share one denominator.
  (b) THE STATISTIC.  Recall is a MEAN over the near-native class; set best is its MINIMUM.
      A gate can raise the mean of a class while dropping its extreme member.  Measured
      directly here: the rank (within the pool's Legacy order) of the pool's single best
      candidate, against the rank of the median sub-2 A candidate.

Also runs the coordinator's question: **does a Legacy-informed shortlist keep a near-native
member that a score-ranked shortlist discards?**  That is a COVERAGE measurement, priced
against a matched random shortlist of the same size, on the same targets.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I                     # noqa: E402
from s15.seed import stable_rng                     # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s17 import phys_lib as P                       # noqa: E402
from s17.phys_report import load, BAND, NEAR, GATE_F, N_RAND    # noqa: E402

SHORT = (25, 75, 150)


def main(name="phys_ident_leg.json"):
    rows = load(name)
    folds = np.array([r["fold"] for r in rows])
    has = np.array([bool((np.asarray(r["d"], float) < NEAR).any()) for r in rows])
    print(f"n = {len(rows)}; targets WITH a sub-{NEAR} A candidate in K=500: "
          f"{int(has.sum())}/{len(rows)}.  Every statistic below is on that stratum, so "
          f"recall and set-best share one denominator.")
    out = {"n": int(len(rows)), "n_stratum": int(has.sum())}

    for gate_col in ("s_legacy", "s_leg_torsion", "s_leg_steric"):
        print(f"\n{'='*112}\nGATE = {gate_col}   (reject the worst f by this score, keep the "
              f"rest)\n{'='*112}")
        print(f"  {'f':>6}{'near recall':>13}{'rand recall':>13}"
              f"{'  recall gate-rand [CI] W/L':>36}{'set best':>10}{'rand best':>11}"
              f"{'  best gate-rand [CI] W/L':>34}")
        g = {}
        for f in GATE_F:
            rec, rrec, bst, rbst, prank, mrank = [], [], [], [], [], []
            for r, ok in zip(rows, has):
                if not ok:
                    continue
                s = np.asarray(r[gate_col], float)
                d = np.asarray(r["d"], float)
                k = len(d)
                nkeep = int(round(k * (1.0 - f)))
                near = d < NEAR
                srt = np.argsort(s, kind="stable")[:nkeep]
                rec.append(float(near[srt].sum() / near.sum()))
                bst.append(float(d[srt].min()))
                rr, rb = [], []
                for t in range(N_RAND):
                    rng = stable_rng("s17phys", "gate2", gate_col, f, r["pdb"], t)
                    pick = rng.permutation(k)[:nkeep]
                    rr.append(float(near[pick].sum() / near.sum()))
                    rb.append(float(d[pick].min()))
                rrec.append(np.mean(rr)); rbst.append(np.mean(rb))
                #: (b) the statistic: where does the SINGLE BEST candidate sit in the gate's
                #: own order, against where the MEDIAN near-native sits?  Both in [0, 1];
                #: the gate keeps the lowest fraction (1 - f).
                rk = EL._avg_rank(s) / k
                prank.append(float(rk[int(np.argmin(d))]))
                mrank.append(float(np.median(rk[near])))
            pr = P.paired(rec, rrec, folds=folds[has])
            pb = P.paired(bst, rbst, folds=folds[has])
            print(f"  {f:>6.2f}{np.mean(rec):>13.3f}{np.mean(rrec):>13.3f}"
                  f"   {-pr['mean']:+.4f} [{-pr['ci'][1]:+.4f},{-pr['ci'][0]:+.4f}] "
                  f"{pr['L']:>3}W/{pr['W']:<3}L"
                  f"{np.mean(bst):>10.3f}{np.mean(rbst):>11.3f}"
                  f"   {pb['mean']:+.4f} [{pb['ci'][0]:+.4f},{pb['ci'][1]:+.4f}] "
                  f"{pb['W']:>3}W/{pb['L']:<3}L")
            g[f] = {"recall": float(np.mean(rec)), "rand_recall": float(np.mean(rrec)),
                    "recall_vs_rand": pr, "best": float(np.mean(bst)),
                    "rand_best": float(np.mean(rbst)), "best_vs_rand": pb,
                    "rank_of_pool_best": float(np.mean(prank)),
                    "rank_of_median_near": float(np.mean(mrank))}
        pr = g[GATE_F[0]]
        print(f"\n  THE MECHANISM.  Normalised rank in the gate's OWN order (0 = kept first, "
              f"1 = rejected first):")
        print(f"     the pool's single BEST candidate sits at    "
              f"{pr['rank_of_pool_best']:.3f}")
        print(f"     the MEDIAN sub-{NEAR} A candidate sits at        "
              f"{pr['rank_of_median_near']:.3f}")
        print(f"     a random candidate sits at                  0.500")
        print(f"  If the best sits HIGHER than the median near-native, the gate rejects the "
              f"extreme member\n  of the class it otherwise protects -- which is exactly the "
              f"set-mean trap, at the level of one candidate.")
        out[gate_col] = {str(k): v for k, v in g.items()}

    # ---------------------------------------------------------------- coverage
    print(f"\n{'='*112}\nCOVERAGE -- can a Legacy-informed shortlist keep a near-native the "
          f"SCORE-ranked shortlist discards?\n{'='*112}")
    print("  Shortlists of size M, all native-free.  `dist` = top-M by the shipped distance")
    print("  score.  `mix` = top-M/2 by distance UNION top-M/2 by Legacy.  `rand` = a random")
    print("  M from the pool -- the matched random control of the same size.")
    print(f"  {'M':>5}{'dist has<2A':>13}{'mix has<2A':>12}{'rand has<2A':>13}"
          f"{'  mix - dist best (A) [CI] W/L':>38}{'dist best':>11}{'mix best':>10}")
    cov = {}
    for M in SHORT:
        hd, hm, hr, bd, bm = [], [], [], [], []
        for r, ok in zip(rows, has):
            if not ok:
                continue
            d = np.asarray(r["d"], float)
            dist = np.asarray(r["s_dist"], float)
            leg = np.asarray(r["s_legacy"], float)
            k = len(d)
            sd_ = np.argsort(dist, kind="stable")[:M]
            half = max(1, M // 2)
            mix = np.unique(np.concatenate([np.argsort(dist, kind="stable")[:half],
                                            np.argsort(leg, kind="stable")[:half]]))
            rng = stable_rng("s17phys", "cover", M, r["pdb"])
            rp = rng.permutation(k)[:M]
            hd.append(float((d[sd_] < NEAR).any()))
            hm.append(float((d[mix] < NEAR).any()))
            hr.append(float((d[rp] < NEAR).any()))
            bd.append(float(d[sd_].min())); bm.append(float(d[mix].min()))
        pb = P.paired(bm, bd, folds=folds[has])
        print(f"  {M:>5}{np.mean(hd):>13.3f}{np.mean(hm):>12.3f}{np.mean(hr):>13.3f}"
              f"   {pb['mean']:+.4f} [{pb['ci'][0]:+.4f},{pb['ci'][1]:+.4f}] "
              f"{pb['W']:>3}W/{pb['L']:<3}L{np.mean(bd):>11.3f}{np.mean(bm):>10.3f}")
        cov[M] = {"dist_has": float(np.mean(hd)), "mix_has": float(np.mean(hm)),
                  "rand_has": float(np.mean(hr)), "dist_best": float(np.mean(bd)),
                  "mix_best": float(np.mean(bm)), "mix_vs_dist_best": pb}
    out["coverage"] = {str(k): v for k, v in cov.items()}
    json.dump(out, open(os.path.join(RESULTS, "phys_gate.json"), "w"),
              default=lambda x: float(x) if isinstance(x, np.floating) else str(x))
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phys_ident_leg.json")
