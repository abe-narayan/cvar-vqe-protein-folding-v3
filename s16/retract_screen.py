"""SPRINT 16 / RETRACT -- TASK A part 3: is `s` an a-priori NATIVE-FREE fusion screen?

LIT (`s16/lit_FINDINGS.md` B.5) says two things survive the Krogh-Vedelsby attribution:
  (a) that the identity holds through Kabsch superposition on 126 real targets  -> `s16/retract_exact.py`
  (b) that `s` is native-free, so the law works as an *a priori* screen for whether fusion is worth
      running.
The task is to TEST (b), not repeat it.

THE HOSTILE READING, stated before the numbers.  The law is

    d_avg = sqrt( rbar_quad^2 - (s/2)^2 ),   rbar_quad = sqrt((r1^2+r2^2)/2)

and `rbar_quad` is a mean of RMSDs TO THE NATIVE.  It is not native-free.  So the law's PREDICTION
is not native-free; only its INPUT `s` is.  Two separate claims must therefore be separated:

  B1  can `s` alone predict the fusion gain?  (the absolute gain needs `r`, so this asks how much of
      the per-target gain variance `s` carries on its own)
  B2  can `s` decide the only question a practitioner actually faces -- **does the coordinate average
      beat the better single channel?**  The identity guarantees d_avg <= rbar_quad for ANY s > 0
      (that IS the ambiguity decomposition: the ensemble never beats the mean member by less than the
      diversity).  It says NOTHING about beating the BEST member.  If `s` cannot separate the targets
      where fusion wins from those where it loses, the screen claim is empty.

Controls: the trivial screens `always fuse` and `never fuse`, and the ORACLE per-target pick.
Nothing here refits anything; it reads `s15/results/coherence.json` per-target rows.

Run:  python -m s16.retract_screen
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
from s16.retract_law import fold_cluster_boot   # noqa: E402

RES = os.path.join(ROOT, "s16", "results")


def _spearman(a, b):
    from scipy.stats import rankdata
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def run():
    rows = json.load(open(os.path.join(ROOT, "s15", "results", "coherence.json")))["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([int(r["fold"]) for r in rows], int)
    r1 = np.asarray([r["rmsd_disto"] for r in rows], float)     # the globally better channel
    r2 = np.asarray([r["rmsd_pool"] for r in rows], float)
    s = np.asarray([r["disagreement"] for r in rows], float)
    avg = np.asarray([r["rmsd_coordavg"] for r in rows], float)
    n = np.asarray([r["n"] for r in rows], float)
    rq = np.sqrt(0.5 * (r1 ** 2 + r2 ** 2))

    gain_vs_mean = rq - avg                # the identity's quantity; positive by construction
    gain_vs_better = r1 - avg              # the practitioner's quantity (disto is better on average)
    gain_vs_oracle = np.minimum(r1, r2) - avg

    out = {"n": len(rows),
           "which_channel_is_better_globally": "rmsd_disto (%.4f vs %.4f)" % (r1.mean(), r2.mean())}

    # ---- B1: how much of the gain does the native-free s carry on its own? -------------------
    out["B1_gain_vs_mean"] = {
        "mean": float(gain_vs_mean.mean()), "median": float(np.median(gain_vs_mean)),
        "n_positive": int((gain_vs_mean > 0).sum()),
        "pearson_with_s": float(np.corrcoef(s, gain_vs_mean)[0, 1]),
        "spearman_with_s": _spearman(s, gain_vs_mean),
        "pearson_with_s_over_rq_ORACLE": float(np.corrcoef(s / rq, gain_vs_mean)[0, 1]),
        "pearson_with_chain_length": float(np.corrcoef(n, gain_vs_mean)[0, 1]),
    }
    out["B1_gain_vs_better_channel"] = {
        "mean": float(gain_vs_better.mean()), "median": float(np.median(gain_vs_better)),
        "n_positive_fusion_wins": int((gain_vs_better > 0).sum()),
        "n_negative_fusion_loses": int((gain_vs_better < 0).sum()),
        "pearson_with_s": float(np.corrcoef(s, gain_vs_better)[0, 1]),
        "spearman_with_s": _spearman(s, gain_vs_better),
        "pearson_with_s_over_rq_ORACLE": float(np.corrcoef(s / rq, gain_vs_better)[0, 1]),
        "spearman_with_s_over_rq_ORACLE": _spearman(s / rq, gain_vs_better),
    }

    # ---- B2: can a NATIVE-FREE threshold on s beat "always fuse"? ----------------------------
    # Leave-fold-out: the threshold is chosen on 4 folds by minimising mean emitted RMSD, applied
    # to the held-out fold.  The choice uses native labels (it minimises RMSD) -- FOLD-HONEST
    # SUPERVISED, exactly the status of the step size in s16/steer.py.  Reported as such.
    grid = np.quantile(s, np.linspace(0.0, 1.0, 41))
    emitted_lfo = np.empty(len(rows))
    chosen = {}
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        best_t, best_v = None, None
        for t in grid:
            v = float(np.where(s[tr] > t, avg[tr], r1[tr]).mean())
            if best_v is None or v < best_v:
                best_v, best_t = v, float(t)
        chosen[int(f)] = best_t
        emitted_lfo[te] = np.where(s[te] > best_t, avg[te], r1[te])

    always = avg.copy()
    never = r1.copy()
    oracle = np.minimum(r1, avg)
    out["B2_screen"] = {
        "lfo_thresholds_per_fold_SUPERVISED": chosen,
        "emitted_screen_lfo_mean": float(emitted_lfo.mean()),
        "emitted_screen_lfo_median": float(np.median(emitted_lfo)),
        "always_fuse_mean": float(always.mean()), "always_fuse_median": float(np.median(always)),
        "never_fuse_mean": float(never.mean()), "never_fuse_median": float(np.median(never)),
        "ORACLE_per_target_pick_mean": float(oracle.mean()),
        "screen_vs_always": I.paired(emitted_lfo, always, folds=folds, names=pdbs),
        "screen_vs_always_foldclustered_ci95": fold_cluster_boot(emitted_lfo - always, folds),
        "screen_vs_never": I.paired(emitted_lfo, never, folds=folds, names=pdbs),
        "screen_vs_never_foldclustered_ci95": fold_cluster_boot(emitted_lfo - never, folds),
        "ORACLE_headroom_always_minus_oracle": float(always.mean() - oracle.mean()),
    }

    # ---- B2b: is `s` even correlated with the SIGN of the fusion decision? -------------------
    win = (gain_vs_better > 0).astype(float)
    # AUC of s as a ranker of "fusion wins here"
    order = np.argsort(s)
    ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    npos, nneg = win.sum(), (1 - win).sum()
    auc = float((ranks[win == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))
    out["B2b_sign_discrimination"] = {
        "base_rate_fusion_wins": float(win.mean()),
        "AUC_of_s_for_fusion_wins": auc,
        "AUC_null": 0.5,
        "point_biserial_s_vs_win": float(np.corrcoef(s, win)[0, 1]),
        "mean_s_when_fusion_wins": float(s[win == 1].mean()),
        "mean_s_when_fusion_loses": float(s[win == 0].mean()),
    }

    json.dump(out, open(os.path.join(RES, "retract_screen.json"), "w"), indent=1)

    b = out["B1_gain_vs_mean"]
    print(f"n = {out['n']}.  better channel globally: {out['which_channel_is_better_globally']}\n")
    print("B1  gain of the coordinate average over ...")
    print(f"  the QUADRATIC MEAN member (the identity's quantity): mean {b['mean']:+.4f} "
          f"median {b['median']:+.4f}, positive on {b['n_positive']}/{out['n']}")
    print(f"       corr with s {b['pearson_with_s']:+.3f} (spearman {b['spearman_with_s']:+.3f}); "
          f"with s/rbar (NOT native-free) {b['pearson_with_s_over_rq_ORACLE']:+.3f}")
    c = out["B1_gain_vs_better_channel"]
    print(f"  the BETTER SINGLE CHANNEL (the practitioner's quantity): mean {c['mean']:+.4f} "
          f"median {c['median']:+.4f}, wins {c['n_positive_fusion_wins']}/"
          f"{c['n_negative_fusion_loses']}")
    print(f"       corr with s {c['pearson_with_s']:+.3f} (spearman {c['spearman_with_s']:+.3f}); "
          f"with s/rbar (NOT native-free) {c['pearson_with_s_over_rq_ORACLE']:+.3f}")
    d = out["B2_screen"]
    print(f"\nB2  screen 'fuse iff s > t', t chosen leave-fold-out (SUPERVISED, minimises RMSD)")
    print(f"  thresholds {d['lfo_thresholds_per_fold_SUPERVISED']}")
    print(f"  screen {d['emitted_screen_lfo_mean']:.4f} | always fuse {d['always_fuse_mean']:.4f} "
          f"| never fuse {d['never_fuse_mean']:.4f} | ORACLE pick {d['ORACLE_per_target_pick_mean']:.4f}")
    w = d["screen_vs_always"]
    print(f"  screen - always {w['mean_diff']:+.4f} [{w['ci95'][0]:+.4f},{w['ci95'][1]:+.4f}] iid | "
          f"fold-clustered [{d['screen_vs_always_foldclustered_ci95'][0]:+.4f},"
          f"{d['screen_vs_always_foldclustered_ci95'][1]:+.4f}] median {w['median_diff']:+.4f} "
          f"W/L {w['n_better']}/{w['n_worse']}")
    w = d["screen_vs_never"]
    print(f"  screen - never  {w['mean_diff']:+.4f} [{w['ci95'][0]:+.4f},{w['ci95'][1]:+.4f}] iid | "
          f"fold-clustered [{d['screen_vs_never_foldclustered_ci95'][0]:+.4f},"
          f"{d['screen_vs_never_foldclustered_ci95'][1]:+.4f}] median {w['median_diff']:+.4f} "
          f"W/L {w['n_better']}/{w['n_worse']}")
    e = out["B2b_sign_discrimination"]
    print(f"\nB2b s as a ranker of 'fusion wins here': AUC {e['AUC_of_s_for_fusion_wins']:.3f} "
          f"(null 0.500), base rate {e['base_rate_fusion_wins']:.3f}; "
          f"mean s win {e['mean_s_when_fusion_wins']:.3f} vs lose {e['mean_s_when_fusion_loses']:.3f}")
    return out


if __name__ == "__main__":
    run()
