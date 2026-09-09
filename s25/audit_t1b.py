"""s25/audit_t1b.py -- AUDIT part B: L2's mechanism, the weight/shape split, and the best-of-K null.

Independent of s25/temper.py.  Reads only its persisted artefact s25/results/temper.json and the
cached distograms.  No training.
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

from s12 import instrument as I           # noqa: E402
from s24 import stats_lib as ST           # noqa: E402
from s25.audit_t1 import widen_sd, temper, moments, discrete_median  # noqa: E402

CACHE = os.path.join(ROOT, "s12", "cache")
SEPBANDS = [(2, 2), (3, 3), (4, 5), (6, 8), (9, 15)]


def main():
    A = json.load(open(os.path.join(RES, "temper.json")))
    rows = A["rows"]
    FG, TG = A["fgrid"], A["tgrid"]
    pdbs = [r["pdb"] for r in rows]
    fold = np.array([r["fold"] for r in rows], int)
    R = np.array([r["SD"] for r in rows], float)
    RF = np.array([r["SDFIXW"] for r in rows], float)
    RT = np.array([r["TEMP"] for r in rows], float)
    base = R[:, 0]
    n = len(rows)

    print("=" * 100)
    print("B2  WEIGHT CHANNEL vs SHAPE CHANNEL.  L2 says width is 'very nearly irrelevant'.")
    print("=" * 100)
    print("  The shipped per-pair weight is  w_p = shell_p / (sd_p + 0.5)^g  with shell = 1 and g = 1")
    print("  (core/predict.py:420, score_weights.json intentionally absent).  shell is CONSTANT, so")
    print("  the ONLY per-pair weight in the shipped score is a pure function of the posterior WIDTH.")
    print()
    print("  %8s%12s%12s%12s%14s%14s" % ("f", "SD", "SDFIXW", "base", "total effect", "shape-only"))
    for a, f in enumerate(FG):
        tot = (R[:, a] - base).mean(); shp = (RF[:, a] - base).mean()
        print("  %8.2f%12.4f%12.4f%12.4f%+14.4f%+14.4f" % (f, R[:, a].mean(), RF[:, a].mean(),
                                                           base.mean(), tot, shp))
    print()
    for a, f in enumerate(FG):
        if f == 1.0:
            continue
        tot = (R[:, a] - base).mean(); shp = (RF[:, a] - base).mean()
        print("    f=%-5.2f  total %+.4f = shape %+.4f (%.0f%%) + weight %+.4f (%.0f%%)"
              % (f, tot, shp, 100 * shp / tot, tot - shp, 100 * (tot - shp) / tot))

    print("\n  The three headline contrasts, through stats_lib (negative = better than shipped):")
    for lab, arm in (("SD f=1.5  vs shipped", R[:, FG.index(1.5)]),
                     ("SDFIXW f=1.5 (shape only) vs shipped", RF[:, FG.index(1.5)]),
                     ("SD f=3.0  vs shipped", R[:, FG.index(3.0)]),
                     ("SDFIXW f=3.0 (shape only) vs shipped", RF[:, FG.index(3.0)]),
                     ("TEMP T=5.0 vs shipped", RT[:, TG.index(5.0)])):
        print(ST.fmt(ST.compare(arm, base, fold, names=pdbs, label=lab)))

    print("\n" + "=" * 100)
    print("B3  THE ORACLE ROW'S BEST-OF-8 NULL.  temper.py's null re-samples the TARGET'S OWN")
    print("    8 columns WITH REPLACEMENT, so it is close to self-referential.")
    print("=" * 100)
    dev = R - R.mean(1, keepdims=True)
    obs = float((R.min(1) - R.mean(1)).mean())
    rng = np.random.default_rng(20250908)
    K = R.shape[1]

    # (a) temper.py's null, reimplemented: min of K WITH-REPLACEMENT draws from the SAME row.
    a_null = np.mean([np.minimum.reduce([dev[np.arange(n), rng.integers(0, K, n)]
                                         for _ in range(K)]).mean() for _ in range(400)])
    # (a') its analytic ceiling: min-of-K-with-replacement from a row's OWN K values can only
    #      ever return one of those K values, and returns the true min with prob 1-(1-1/K)^K.
    p_incl = 1.0 - (1.0 - 1.0 / K) ** K
    # (b) EXCHANGEABLE-COLUMN null: deviations resampled ACROSS targets, so any per-target
    #     structure in WHICH f wins is destroyed but the marginal spread is preserved.
    flat = dev.ravel()
    b_null = np.mean([flat[rng.integers(0, len(flat), (n, K))].min(1).mean() for _ in range(400)])
    # (c) SPLIT-HALF TRANSFER, the brief's preferred construction: choose f on one random half
    #     of targets, score it on the other half.  Nulls itself.
    sh = []
    for _ in range(400):
        perm = rng.permutation(n); h1, h2 = perm[:n // 2], perm[n // 2:]
        sh.append(0.5 * (dev[h2][:, int(np.argmin(dev[h1].mean(0)))].mean()
                         + dev[h1][:, int(np.argmin(dev[h2].mean(0)))].mean()))
    shm = float(np.mean(sh))
    k_eff = float(np.exp(np.log(K)))

    print("  observed per-target ORACLE gain vs the row mean:            %+.4f A" % obs)
    print("  (a) temper.py's null  min of %d WITH-REPLACEMENT draws from the SAME row"
          "\n      -> %+.4f   share_accounted %.0f%%   residual %+.4f" % (K, a_null,
                                                                         100 * a_null / obs,
                                                                         obs - a_null))
    print("      P(the row's true minimum is drawn at least once) = 1-(1-1/%d)^%d = %.3f."
          % (K, K, p_incl))
    print("      This null can ONLY return one of the row's own 8 values.  It is bounded below by")
    print("      the observed statistic BY CONSTRUCTION and its share is a near-constant of K,")
    print("      not a measurement.  It tests nothing about whether the oracle carries signal.")
    print("  (b) exchangeable-column null  deviations resampled ACROSS targets"
          "\n      -> %+.4f   share_accounted %.0f%%   residual %+.4f" % (b_null,
                                                                         100 * b_null / obs,
                                                                         obs - b_null))
    print("  (c) SPLIT-HALF TRANSFER (the brief's preferred construction; nulls itself)"
          "\n      -> %+.4f A of the %+.4f oracle transfers to held-out targets  (%.0f%% real)"
          % (shm, obs, 100 * shm / obs))
    print("  k_eff (distinct grid points) %.2f" % k_eff)

    print("\n" + "=" * 100)
    print("B4  IS THE SD OPERATOR'S LOCATION SHIFT SEPARATION-GRADED?")
    print("=" * 100)
    print("  s = sd_p * sqrt(f^2-1), so wide (long-range) pairs are convolved harder.  L1 says the")
    print("  needed location correction ALSO grows with separation (-0.048 -> -0.589 A).  If the two")
    print("  match, the SD arm silently APPLIED the per-separation location correction that BRIEF")
    print("  SS3 lists as the leading open candidate -- and it cost RMSD.")
    acc = {k: [] for k in SEPBANDS}
    need = {k: [] for k in SEPBANDS}
    for r in rows:
        z = np.load(os.path.join(CACHE, "disto_%s.npz" % r["pdb"]))
        P = np.asarray(z["prob"], float); C = np.asarray(z["centres"], float)
        i, j = np.asarray(z["i"], int), np.asarray(z["j"], int)
        exp = np.asarray(z["expected"], float)
        u = I.load_univ(r["pdb"]); nat = np.asarray(u["nat_ca"], float)
        Dt = I.pair_dists(nat[None], i, j)[0]
        sep = (j - i).astype(int)
        m0, _ = moments(P, C)
        m15, _ = moments(widen_sd(P, C, 1.5), C)
        for a, b in SEPBANDS:
            m = (sep >= a) & (sep <= b)
            if m.sum():
                acc[(a, b)].append(float((m15 - m0)[m].mean()))
                need[(a, b)].append(float((Dt - exp)[m].mean()))
    print("  %-8s%18s%18s" % ("sep", "SD(1.5) shift", "L1 signed error"))
    for k in SEPBANDS:
        print("  %-8s%+18.4f%+18.4f" % ("%d-%d" % k, np.mean(acc[k]), np.mean(need[k])))
    xs = np.array([np.mean(acc[k]) for k in SEPBANDS])
    ys = np.array([np.mean(need[k]) for k in SEPBANDS])
    print("  correlation across the 5 bands: %.4f   ratio (shift/needed) %s"
          % (np.corrcoef(xs, ys)[0, 1], np.round(xs / ys, 3).tolist()))


if __name__ == "__main__":
    main()
