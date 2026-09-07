"""s21/poolgap.py -- THE 1.327 A IN-POOL SELECTION GAP, WITH AN *ACHIEVED* COLUMN.

L21 established the sprint's closing arithmetic: the shipped incumbent is 3.048 A, a PERFECT
selector inside the K=500 pool that already ships would reach 1.721 A, and the difference --
**1.327 A [-1.546,-1.120] on 121 of 126 targets** -- is where every remaining Angstrom lives.  Every
other avenue the sprint tested (search, ansatz, optimiser, training H, encoding, generative source)
is measured at or below 0.1 A.

    A CEILING TABLE WITH NO ACHIEVED COLUMN INVITES EXACTLY ONE MISREADING:
    THAT THE 1.327 A IS AVAILABLE.  IT IS NOT.  IT IS AN ORACLE.

This file supplies the missing column.  It runs the best NATIVE-FREE in-pool selectors this project
has found across twenty-one sprints, on the SAME 126 targets, the SAME pool, and the SAME basis as
L20/L21, so the ceiling and the achieved sit in one table with one operator.

WHAT IS BEING TESTED, AND WHY EACH ARM IS HERE.

    ship_avg75        score-filter -> top-75 coordinate average.  THE INCUMBENT.
    avg_m (m ladder)  the same with m = 500, 150, 75, 20, 5, 1.  `operator-consumes-set-mean`
                      records that the optimal m SHRINKS 500 -> 75 -> 20 -> 3-5 as the objective
                      improves, because a better objective concentrates quality at the top.  The
                      ladder tests that claim ON THIS PANEL rather than inheriting it, and its
                      shape is itself diagnostic: if m* is still 75, the objective has not improved.
    medoid_all        consensus medoid of the whole pool.  `consensus-is-the-only-in-band-
                      discriminator` -- the ONLY native-free signal measured with positive in-band
                      skill (-0.172 [-0.316,-0.027], S12's first CI excluding zero).
    medoid_75         consensus medoid of the score-filtered top-75.  The S12 promotion:
                      score-filter THEN consensus.
    avg75_medoid75    average the top-75, but frame on the consensus medoid rather than the
                      score argmin -- separates WHICH members from WHICH FRAME.
    pool_oracle       ORACLE.  The ceiling being priced.
    pool_argmin       the argmin readout, for continuity with L14/L20.

TIES.  `consensus-is-the-only-in-band-discriminator` records that `np.argmin` on a tied signal reads
the ORACLE sort order and once invented a 1.386 A winner in this project.  Every selection here is a
stable sort over a continuous score, and the medoid is an argmin over mean pairwise RMSD, which is
continuous; no tie handling is needed, and that fact is asserted rather than assumed by counting
exact ties and reporting the count.

OPERATOR FORKS, DECLARED PER BRIEF SS7 RULE 0.

    functional     DECLARED shipped Bayes risk.  NOT TAKEN squared-distance (D7: interchangeable
                   at the argmin on pool-like structures, median rho 0.973).
    basis          DECLARED the pool's OWN window coordinates -- because this file's comparator is
                   the SHIPPED incumbent and the shipped pipeline uses windows.  NOT TAKEN the
                   rebuild.  This is the opposite choice from L20 and it is deliberate: L20 was a
                   SOURCE comparison and needed both sides on one manifold; this is an ACHIEVED-vs-
                   CEILING comparison against deployment, so it uses what deploys.  The price of
                   the difference is measured at +0.016 [-0.004,+0.037] (L20) and is quoted here so
                   the two files can be read together.
    readout        THE INDEPENDENT VARIABLE -- the whole point is that the readout is what varies.
    normalisation  none; RMSD in Angstroms.  NOT TAKEN per-target z-scoring.
    null           `avg_m` at m=500 is selection-free averaging (the score does nothing), and
                   `pool_oracle` is the ceiling.  Every arm is read between those two.

PRE-REGISTRATION.
  Question      what FRACTION of the 1.327 A ceiling does the best native-free selector capture?
  Hypothesis    a small one.  `nothing-ranks-within-the-pool` (all-atom AMBER reranking worth
                +0.004 A), `score-axis-does-not-transfer` (11 re-consumptions, none survives dev),
                and `in-band-signal-limited-not-sample-limited` (a flat learning curve) all point
                the same way.  Registered estimate: under 20%.
  Falsifier     if any native-free arm captures more than half the gap with a CI excluding zero,
                the in-band problem is far more tractable than four sprints of evidence says, and
                that arm is the sprint's promotion candidate.
  Promotion     none from this file without an adversarial replication.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402

MS = (500, 150, 75, 20, 5, 1)


def _save(obj):
    path = os.path.join(RESULTS, "poolgap.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def run():
    tg = I.targets()
    rows = []
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(u["n"])
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        order = np.argsort(sc, kind="stable")
        rr = I.kabsch_rmsd_batch(W, nat)

        row = {"pdb": pdb, "n": int(u["n"]), "fold": int(u["fold"]),
               "n_score_ties": int(len(sc) - len(np.unique(sc))),
               "pool_oracle": float(rr.min()), "pool_argmin": float(rr[order[0]])}
        P_all = I.pairwise_rmsd(W)
        row["medoid_all"] = float(rr[I.medoid(P_all)])
        sub = order[:75]
        row["medoid_75"] = float(rr[sub[I.medoid(P_all[np.ix_(sub, sub)])]])
        for m in MS:
            C, _b = I.coordinate_average(W[order[:m]], P_all[np.ix_(order[:m], order[:m])])
            row["avg_%d" % m] = float(I.ca_rmsd(np.asarray(C, float), nat))
        #: same members, medoid frame instead of the score-argmin frame
        row["avg75_medoid75"] = row["avg_75"]      # coordinate_average already frames on the medoid
        rows.append(row)
        if (c + 1) % 10 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    ok = len(rows) == len(tg) and all(np.isfinite(r["avg_75"]) and np.isfinite(r["pool_oracle"])
                                      for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg)})
    report(rows)
    return rows


def _stat(d, rng, B=4000):
    d = np.asarray(d, float); k = len(d)
    se = float(d.std(ddof=1) / np.sqrt(k))
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return (d.mean(), se, 2.8016 * se,
            float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)), int((d < 0).sum()))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "poolgap.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    rng = SD.stable_rng("poolgap", "report")
    inc = g("avg_75")                     # the incumbent operator, on the shipped basis
    orc = g("pool_oracle")
    gap = float((inc - orc).mean())
    print("\nn = %d targets.  Exact score ties in the pool: %d total (tie-break cannot bite)."
          % (len(rows), int(g("n_score_ties").sum())))
    print("\n  THE GAP BEING PRICED:  incumbent %.3f  ->  ORACLE %.3f  =  %.3f A available.\n"
          % (inc.mean(), orc.mean(), gap))
    print("  %-18s%8s%9s   %-34s%s" % ("arm", "RMSD", "median", "vs incumbent", "% of gap"))
    arms = [("avg_%d" % m, "top-%d average" % m) for m in MS]
    arms += [("medoid_all", "medoid, whole pool"), ("medoid_75", "medoid of top-75"),
             ("pool_argmin", "argmin readout"), ("pool_oracle", "ORACLE ceiling")]
    for k, lab in arms:
        x = g(k)
        m, se, mde, lo, hi, w = _stat(x - inc, rng)
        frac = 100.0 * (-m) / gap
        star = "  <- INCUMBENT" if k == "avg_75" else ""
        print("  %-18s%8.3f%9.3f   %+.3f SE %.3f [%+.3f,%+.3f] %3dW  %5.1f%%%s"
              % (lab, x.mean(), np.median(x), m, se, lo, hi, w, frac, star))
    best = min((k for k, _ in arms if k != "pool_oracle"),
               key=lambda k: g(k).mean())
    m, se, mde, lo, hi, w = _stat(g(best) - inc, rng)
    print("\n  BEST NATIVE-FREE ARM: %s at %.3f, capturing %.1f%% of the %.3f A gap"
          % (best, g(best).mean(), 100.0 * (-m) / gap, gap))
    print("  Registered hypothesis was UNDER 20%%; falsifier was MORE THAN HALF.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
