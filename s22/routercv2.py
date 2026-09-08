"""s22/routercv2.py -- ATTACKING routercv.py's NEGATIVE RESULT before reporting it.

routercv.py's naive full-switch router LOSES to the incumbent on every one of five arm sets
(one significantly: A_plus_latent, CI [+0.005,+0.146]).  Two obvious objections before that goes in
the findings as a disposition:

  (1) IS IT AN IN-FOLD-WINS / OUT-OF-FOLD-LOSES STORY (the exact hazard BRIEF names)?  Fit-and-
      predict on the SAME data (no held-out fold at all) as the upper bound of what this exact
      method could ever achieve with zero generalisation gap.
  (2) DOES A GATED (confidence-thresholded) ROUTER RESCUE IT?  Route away from the incumbent only
      when the predicted margin exceeds a threshold selected by NESTED cross-validation (never on
      the true held-out fold) -- protects against acting on noisy small-margin predictions, which is
      the natural first fix for a router that loses by switching too often.

Both run on the ALL (16-arm) set with the full feature set, the richest and worst-losing
configuration, so the attack is not run only on a set already flattered.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from sklearn.linear_model import RidgeCV

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s15 import seed as SD                  # noqa: E402
from s22.routercv import (SETS, ALL_FEATS, INCUMBENT, ALPHAS, _load, _stat, _fmt)  # noqa: E402

ARMS = SETS["ALL"]
MARGINS = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0, np.inf]


def _fit_predict_all(Xtr, Ytr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1.0
    Xtr_s = (Xtr - mu) / sd; Xte_s = (Xte - mu) / sd
    preds = np.zeros((len(Xte), Ytr.shape[1]))
    for ai in range(Ytr.shape[1]):
        m = RidgeCV(alphas=ALPHAS)
        m.fit(Xtr_s, Ytr[:, ai])
        preds[:, ai] = m.predict(Xte_s)
    return preds


def in_fold(rows):
    """Fit and predict on the SAME 126 targets -- the upper bound with zero generalisation gap."""
    X = np.array([[r[f] for f in ALL_FEATS] for r in rows], float)
    Y = np.array([[r["arms"][a] for a in ARMS] for r in rows], float)
    inc_i = ARMS.index(INCUMBENT)
    preds = _fit_predict_all(X, Y, X)
    pick = np.argmin(preds, axis=1)
    chosen = Y[np.arange(len(rows)), pick]
    inc = Y[:, inc_i]
    return chosen, inc


def gated_router(rows, folds_arr, rng):
    """Route away from the incumbent only if the predicted margin beats a threshold chosen by
    NESTED CV (an inner split of the training folds, never the true held-out fold)."""
    X = np.array([[r[f] for f in ALL_FEATS] for r in rows], float)
    Y = np.array([[r["arms"][a] for a in ARMS] for r in rows], float)
    inc_i = ARMS.index(INCUMBENT)
    n = len(rows)
    chosen_rmsd = np.full(n, np.nan)
    chosen_thr = {}
    for f in sorted(set(folds_arr.tolist())):
        tr = np.where(folds_arr != f)[0]
        te = np.where(folds_arr == f)[0]
        # inner split of the training targets: 70/30, stable per fold
        rr = SD.stable_rng("s22", "gated_inner", str(f))
        perm = rr.permutation(len(tr))
        cut = int(0.7 * len(tr))
        inner_tr, inner_val = tr[perm[:cut]], tr[perm[cut:]]
        preds_val = _fit_predict_all(X[inner_tr], Y[inner_tr], X[inner_val])
        best_thr, best_val_mean = None, np.inf
        for thr in MARGINS:
            margin = preds_val[:, inc_i][:, None] - preds_val   # positive = arm predicted better than incumbent
            beats = margin.copy()
            beats[:, inc_i] = -1.0    # never gate out the incumbent itself
            ok = beats >= thr
            pick = np.where(ok.any(1), np.where(ok, preds_val, np.inf).argmin(1), inc_i)
            val_actual = Y[inner_val][np.arange(len(inner_val)), pick].mean()
            if val_actual < best_val_mean:
                best_val_mean, best_thr = val_actual, thr
        chosen_thr[f] = float(best_thr)
        # refit on FULL training folds with the nested-selected threshold, apply to true held-out fold
        preds_te = _fit_predict_all(X[tr], Y[tr], X[te])
        margin = preds_te[:, inc_i][:, None] - preds_te
        beats = margin.copy(); beats[:, inc_i] = -1.0
        ok = beats >= best_thr
        pick = np.where(ok.any(1), np.where(ok, preds_te, np.inf).argmin(1), inc_i)
        chosen_rmsd[te] = Y[te][np.arange(len(te)), pick]
    inc = Y[:, inc_i]
    return chosen_rmsd, inc, chosen_thr


def main():
    rows = _load()
    folds = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("s22", "routercv2", "report")

    print("ATTACK 1 -- IN-FOLD (fit==predict set, zero generalisation gap, NOT a valid estimate,")
    print("            reported only to show whether the method wins with unlimited leakage):")
    chosen, inc = in_fold(rows)
    d = chosen - inc
    s = _stat(d, rng)
    print("  in-fold router vs incumbent: %.4f vs %.4f   %s" % (chosen.mean(), inc.mean(), _fmt(s)))

    print("\nATTACK 2 -- GATED ROUTER, margin threshold chosen by NESTED CV (never the true held-out fold):")
    chosen_g, inc_g, thr = gated_router(rows, folds, rng)
    dg = chosen_g - inc_g
    sg = _stat(dg, rng)
    print("  gated router vs incumbent:   %.4f vs %.4f   %s" % (chosen_g.mean(), inc_g.mean(), _fmt(sg)))
    print("  nested-selected thresholds per outer fold: %s" % thr)

    out = {"in_fold": {"mean": float(chosen.mean()), "stat": s},
           "gated": {"mean": float(chosen_g.mean()), "stat": sg, "thresholds": thr},
           "complete": True, "note": "diagnostic/attack file, no per-target row list; n=126 both attacks"}
    path = os.path.join(RESULTS, "routercv2.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1)
    os.replace(tmp, path)
    print("\nwrote", path)


if __name__ == "__main__":
    main()
