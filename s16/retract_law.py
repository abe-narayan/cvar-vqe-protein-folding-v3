"""SPRINT 16 / RETRACT -- TASK A part 1: the fusion law recomputed with the CORRECT mean.

Sprint 15 `s15/coherence.py` line 123 computes the law's `r` as the ARITHMETIC mean of the two
channels' RMSDs:

    r_bar = 0.5 * (rd + rp)
    pred  = sqrt(max(r_bar**2 - (s/2)**2, 0))

The identity it is quoting -- the two-member case of the Krogh-Vedelsby (1995) ambiguity
decomposition, `d_avg^2 = (1/B) sum_b d_b^2 - (1/B) sum_b ||f_b - f_bar||^2` -- requires the
QUADRATIC mean, sqrt((r1^2 + r2^2)/2).  The two agree only when r1 == r2.

This module needs NO refitting: `s15/results/coherence.json` persists rmsd_disto, rmsd_pool,
disagreement and rmsd_coordavg per target, which is everything the law consumes.  It recomputes
both predictions, reports the residual of each, the paired difference with an i.i.d.-target
bootstrap AND a fold-clustered bootstrap (TARGET / FOLD as the unit), and the distribution of
|r1 - r2| so a reader can see when the substitution bites.

Run:  python -m s16.retract_law
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402

RES = os.path.join(ROOT, "s16", "results")
os.makedirs(RES, exist_ok=True)


def fold_cluster_boot(d, folds, n_boot=20000, seed=17):
    """Two-stage cluster bootstrap: resample FOLDS with replacement, then TARGETS within.

    `I.paired` resamples targets i.i.d., which ignores fold clustering.  BRIEF 3.4 asks for a
    fold-clustered interval, so both are reported.  With only 5 folds this interval is coarse and
    that is stated, not hidden.
    """
    d = np.asarray(d, float); folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = {f: np.flatnonzero(folds == f) for f in uf}
    rng = SD.stable_rng("retract_law", f"foldboot{seed}")
    out = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, len(uf), len(uf))
        acc = []
        for k in pick:
            ii = idx[uf[k]]
            acc.append(d[ii[rng.integers(0, len(ii), len(ii))]])
        out[b] = np.concatenate(acc).mean()
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


def law(r1, r2, s, mean="quad"):
    r1 = np.asarray(r1, float); r2 = np.asarray(r2, float); s = np.asarray(s, float)
    if mean == "quad":
        r2bar = 0.5 * (r1 ** 2 + r2 ** 2)
    elif mean == "arith":
        r2bar = (0.5 * (r1 + r2)) ** 2
    else:
        raise ValueError(mean)
    return np.sqrt(np.maximum(r2bar - (s / 2.0) ** 2, 0.0))


def run():
    src = os.path.join(ROOT, "s15", "results", "coherence.json")
    rows = json.load(open(src))["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([int(r["fold"]) for r in rows], int)
    r1 = np.asarray([r["rmsd_disto"] for r in rows], float)
    r2 = np.asarray([r["rmsd_pool"] for r in rows], float)
    s = np.asarray([r["disagreement"] for r in rows], float)
    meas = np.asarray([r["rmsd_coordavg"] for r in rows], float)
    orig = np.asarray([r["law_prediction"] for r in rows], float)

    p_ar = law(r1, r2, s, "arith")
    p_qu = law(r1, r2, s, "quad")
    assert np.allclose(p_ar, orig, atol=1e-9), "arithmetic reconstruction must match the record"

    e_ar = meas - p_ar          # positive = the law UNDER-predicts the measured coordinate average
    e_qu = meas - p_qu
    a_ar, a_qu = np.abs(e_ar), np.abs(e_qu)
    gap = np.abs(r1 - r2)

    out = {
        "n": len(rows),
        "source": "s15/results/coherence.json (persisted per-target rows; NO refit)",
        "reconstruction_matches_record": True,
        "measured_coordavg_mean": float(meas.mean()),
        "measured_coordavg_median": float(np.median(meas)),
        "arith": {
            "pred_mean": float(p_ar.mean()), "pred_median": float(np.median(p_ar)),
            "signed_err_mean": float(e_ar.mean()), "signed_err_median": float(np.median(e_ar)),
            "abs_err_mean": float(a_ar.mean()), "abs_err_median": float(np.median(a_ar)),
            "abs_err_sd": float(a_ar.std(ddof=1)),
            "n_underpredict": int((e_ar > 0).sum()), "n_overpredict": int((e_ar < 0).sum()),
        },
        "quad": {
            "pred_mean": float(p_qu.mean()), "pred_median": float(np.median(p_qu)),
            "signed_err_mean": float(e_qu.mean()), "signed_err_median": float(np.median(e_qu)),
            "abs_err_mean": float(a_qu.mean()), "abs_err_median": float(np.median(a_qu)),
            "abs_err_sd": float(a_qu.std(ddof=1)),
            "n_underpredict": int((e_qu > 0).sum()), "n_overpredict": int((e_qu < 0).sum()),
        },
        "abs_err_quad_minus_arith": I.paired(a_qu, a_ar, folds=folds, names=pdbs),
        "abs_err_quad_minus_arith_foldclustered_ci95": fold_cluster_boot(a_qu - a_ar, folds),
        "pred_quad_minus_arith": I.paired(p_qu, p_ar, folds=folds, names=pdbs),
        "pred_quad_minus_arith_foldclustered_ci95": fold_cluster_boot(p_qu - p_ar, folds),
        "channel_gap_abs_r1_minus_r2": I.summary(gap),
        "channel_gap_quantiles": {q: float(np.quantile(gap, q / 100.0))
                                  for q in (0, 10, 25, 50, 75, 90, 95, 99, 100)},
        "substitution_size": I.summary(p_qu - p_ar),
        "corr_gap_vs_substitution": float(np.corrcoef(gap, p_qu - p_ar)[0, 1]),
        "rows": [{"pdb": p, "fold": int(f), "r_disto": float(x1), "r_pool": float(x2),
                  "gap": float(g), "s": float(ss), "measured": float(m),
                  "pred_arith": float(pa), "pred_quad": float(pq),
                  "err_arith": float(ea), "err_quad": float(eq)}
                 for p, f, x1, x2, g, ss, m, pa, pq, ea, eq
                 in zip(pdbs, folds, r1, r2, gap, s, meas, p_ar, p_qu, e_ar, e_qu)],
    }
    with open(os.path.join(RES, "retract_law.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"n = {out['n']} targets, measured coordinate average "
          f"{out['measured_coordavg_mean']:.4f} (median {out['measured_coordavg_median']:.4f})\n")
    for k in ("arith", "quad"):
        b = out[k]
        print(f"  {k:5s} mean pred {b['pred_mean']:.4f} | signed err {b['signed_err_mean']:+.4f} "
              f"(median {b['signed_err_median']:+.4f}) | |err| {b['abs_err_mean']:.4f} "
              f"(median {b['abs_err_median']:.4f}) | under/over {b['n_underpredict']}/{b['n_overpredict']}")
    w = out["abs_err_quad_minus_arith"]
    print(f"\n  |err| quad - arith  {w['mean_diff']:+.4f} [{w['ci95'][0]:+.4f},{w['ci95'][1]:+.4f}] "
          f"i.i.d.-target | fold-clustered "
          f"[{out['abs_err_quad_minus_arith_foldclustered_ci95'][0]:+.4f},"
          f"{out['abs_err_quad_minus_arith_foldclustered_ci95'][1]:+.4f}] "
          f"median {w['median_diff']:+.4f}  W/L {w['n_better']}/{w['n_worse']}")
    g = out["channel_gap_abs_r1_minus_r2"]
    print(f"\n  |r1 - r2| mean {g['mean']:.4f} median {g['median']:.4f} "
          f"max {g['max']:.4f}; quantiles "
          + ", ".join(f"p{q}={v:.3f}" for q, v in out["channel_gap_quantiles"].items()))
    sub = out["substitution_size"]
    print(f"  pred_quad - pred_arith mean {sub['mean']:+.4f} median {sub['median']:+.4f} "
          f"max {sub['max']:+.4f}; corr with the gap {out['corr_gap_vs_substitution']:+.4f}")
    return out


if __name__ == "__main__":
    run()
