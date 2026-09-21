#!/usr/bin/env python
"""s31/s31_E4_predictability.py -- S31 lane E: what is PREDICTABLE vs what is VALUABLE, same axis.

Refit S30's leave-fold-out ridge on the SAME features and the SAME pinned folds, but with the
response replaced by each half of the orthogonal decomposition of the ORACLE prior error:

    y = y_along + y_perp     (split against the ORACLE common mode mu, per target)

Out-of-fold R^2 for each response answers "which half can the pool + distogram predict at all?"
The endpoint arms (`s31_E2_agg.py`) answer "which half is worth removing?".  If the two answers
name different halves, that is the bottleneck, stated exactly rather than by analogy.

Controls, matched to the operator's space:
  * matched-dimension random features (same ridge, same folds)
  * the SHARED-REFERENT FLOOR: `y_along` is a multiple of `mu`, and `mu` shares `-d_nat` with `y`,
    so a non-zero R^2 is expected before any information.  The random-feature arm measures it.

ORACLE / NOT DEPLOYABLE: the responses read `d_nat` twice.  This is a diagnostic, not a method.

    python s31/s31_E4_predictability.py
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s31 import s31_E_lib as L                     # noqa: E402
from s31.s31_E2_deltas import p_along, p_perp      # noqa: E402

OUT = os.path.join(HERE, "results", "s31_E4_predictability.json")
SEED = 310504


def main():
    d = L.load()
    X, y, mu, mh = d["X"], d["y"], d["mu"], d["mu_hat"]
    tgt, fold = d["tgt"], d["fold"]
    feats = [str(f) for f in d["feats"]]
    cols = [feats.index(f) for f in L.NEST["N3_plus_pool"]]
    rng = np.random.default_rng(SEED)

    ya = np.empty_like(y); yp = np.empty_like(y)
    for _, m in L.per_target(y, tgt):
        ya[m] = p_along(y[m], mu[m]); yp[m] = p_perp(y[m], mu[m])
    ss = lambda a: float((a ** 2).sum())

    Xr = rng.normal(size=X[:, cols].shape)
    out = {"seed": SEED, "note": "ORACLE / NOT DEPLOYABLE -- the responses read d_nat",
           "energy_share_pooled": {"along": ss(ya) / ss(y), "perp": ss(yp) / ss(y)},
           "responses": {}}
    for nm, resp in (("y_FULL", y), ("y_ALONG_mu", ya), ("y_PERP_mu", yp)):
        yh = L.lfo_ridge(X[:, cols], resp, fold)
        yhr = L.lfo_ridge(Xr, resp, fold)
        out["responses"][nm] = {
            "R2_oof": float(1 - ss(resp - yh) / ss(resp - resp.mean())),
            "R2_oof_randfeat_control": float(1 - ss(resp - yhr) / ss(resp - resp.mean())),
            "R2_oof_vs_zero": float(1 - ss(resp - yh) / ss(resp)),
            "sd_response": float(resp.std()),
            "corr_pred_truth_pooled": float(np.corrcoef(yh, resp)[0, 1]),
            "share_of_full_energy": ss(resp) / ss(y)}
        out["responses"][nm]["excess_over_randfeat"] = (
            out["responses"][nm]["R2_oof"] - out["responses"][nm]["R2_oof_randfeat_control"])
    # and how much of the FULL corrector's prediction lands in each half
    yhf = L.lfo_ridge(X[:, cols], y, fold)
    ha = np.empty_like(y); hp = np.empty_like(y)
    for _, m in L.per_target(y, tgt):
        ha[m] = p_along(yhf[m], mu[m]); hp[m] = p_perp(yhf[m], mu[m])
    out["full_corrector_spends"] = {
        "R2_of_its_ALONG_part_against_y_along": float(1 - ss(ya - ha) / ss(ya - ya.mean())),
        "R2_of_its_PERP_part_against_y_perp": float(1 - ss(yp - hp) / ss(yp - yp.mean())),
        "energy_share_of_yhat_along": ss(ha) / ss(yhf), "energy_share_of_yhat_perp": ss(hp) / ss(yhf)}

    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, OUT)
    print("pooled energy share of y:  along %.4f  perp %.4f"
          % (out["energy_share_pooled"]["along"], out["energy_share_pooled"]["perp"]))
    for nm, r in out["responses"].items():
        print("%-12s  R2_oof %+.4f   randfeat ctrl %+.4f   excess %+.4f   (energy share %.4f)"
              % (nm, r["R2_oof"], r["R2_oof_randfeat_control"], r["excess_over_randfeat"],
                 r["share_of_full_energy"]))
    print("full corrector: R2 on the ALONG half %+.4f, on the PERP half %+.4f; its own energy "
          "%.4f along / %.4f perp"
          % (out["full_corrector_spends"]["R2_of_its_ALONG_part_against_y_along"],
             out["full_corrector_spends"]["R2_of_its_PERP_part_against_y_perp"],
             out["full_corrector_spends"]["energy_share_of_yhat_along"],
             out["full_corrector_spends"]["energy_share_of_yhat_perp"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
