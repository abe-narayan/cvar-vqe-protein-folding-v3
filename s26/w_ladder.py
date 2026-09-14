#!/usr/bin/env python
"""s26/w_ladder.py -- the own-native (leaked) fold models placed on the S24 prior ladder:
gam_eff and cos of each leaked posterior toward the native, against the measured endpoint
deltas of L44 Part C.  Lane W, Sprint 26.  ORACLE DIAGNOSTIC throughout (reads native
distances).  Pre-registration: s26/PREREG_memorisation_on_the_ladder.md.

    python s26/jobrun.py --agent W --tag CPU --name w_ladder --est-ram 0.4 -- python s26/w_ladder.py
"""
from __future__ import annotations

import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import predict as PR                       # noqa: E402
import w_selfcopy as W                               # noqa: E402
import p_ladder as PL                                # noqa: E402

SLOPE = -2.1496                                      #: S24 L13, MASS0.1 - MASS0.0 over 0.1 (cloud basis)


def main():
    W.require_signoff("w_ladder")
    E = json.load(open(os.path.join(W.RES, "w_selfcopy_endpoint.json")))
    crow = {r["pdb"]: r for r in E["C"]["rows"]}
    tg = I.targets()
    rows = []
    t0 = time.time()
    for t in tg:
        pdb, seq, fold = t["pdb"], t["seq"], t["fold"]
        u = I.load_univ(pdb); nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb); P0 = np.asarray(dg["prob"], float); i, j = np.asarray(dg["i"]), np.asarray(dg["j"])
        Dt = I.pair_dists(nat, i, j)
        C = np.asarray(PR.CENTRES, float)
        nb = np.abs(C[None, :] - Dt[:, None]).argmin(1); O = np.eye(PR.NBINS)[nb]
        E0 = (P0 * C[None]).sum(1)
        for jm in range(5):
            if jm == fold:
                continue
            prob, ii, jj = W.pinned_posterior(seq, jm)
            assert np.array_equal(ii, i) and np.array_equal(jj, j)
            EQ = (prob * C[None]).sum(1)
            pr = PL.progress(prob, P0, O, EQ, E0, Dt)
            cr = crow[pdb]["leaked"][str(jm)]
            rows.append({"pdb": pdb, "fold": fold, "model": jm, **pr,
                         "mae_leaked": float(np.abs(EQ - Dt).mean()), "mae_clean": float(np.abs(E0 - Dt).mean()),
                         "delta_cloud": float(cr["cloud"] - crow[pdb]["clean"]["cloud"]),
                         "delta_arm": float(cr["arm"] - crow[pdb]["clean"]["arm"]),
                         "pred_naive": float(SLOPE * pr["gam_prob"]),
                         "pred_discounted": float(SLOPE * pr["gam_prob"] * pr["cos_prob"]) if np.isfinite(pr["cos_prob"]) else float("nan")})
        del u
        if len(rows) % 40 == 0:
            print("  %d cells (%.0f s)" % (len(rows), time.time() - t0), flush=True)
    g = lambda k: np.array([r[k] for r in rows], float)
    summ = {"n_cells": len(rows),
            "gam_prob_mean": float(np.nanmean(g("gam_prob"))), "cos_prob_mean": float(np.nanmean(g("cos_prob"))), "amp_prob_mean": float(np.nanmean(g("amp_prob"))),
            "gam_loc_mean": float(np.nanmean(g("gam_loc"))), "cos_loc_mean": float(np.nanmean(g("cos_loc"))),
            "mae_leaked_mean": float(g("mae_leaked").mean()), "mae_clean_mean": float(g("mae_clean").mean()),
            "delta_cloud_mean": float(g("delta_cloud").mean()), "delta_arm_mean": float(g("delta_arm").mean()),
            "pred_naive_mean": float(np.nanmean(g("pred_naive"))), "pred_discounted_mean": float(np.nanmean(g("pred_discounted"))),
            "corr_pred_discounted_vs_delta_cloud": float(np.corrcoef(np.nan_to_num(g("pred_discounted")), g("delta_cloud"))[0, 1]),
            "corr_gam_prob_vs_delta_cloud": float(np.corrcoef(np.nan_to_num(g("gam_prob")), g("delta_cloud"))[0, 1]),
            "corr_mae_change_vs_delta_cloud": float(np.corrcoef(g("mae_leaked") - g("mae_clean"), g("delta_cloud"))[0, 1]),
            "mde_cloud_L44": 0.2437, "mde_arm_L44": 0.2463}
    summ["disagreement_discounted_minus_measured_cloud"] = summ["pred_discounted_mean"] - summ["delta_cloud_mean"]
    summ["disagreement_naive_minus_measured_cloud"] = summ["pred_naive_mean"] - summ["delta_cloud_mean"]
    summ["falsifier_disagree_beyond_mde"] = bool(abs(summ["disagreement_discounted_minus_measured_cloud"]) > summ["mde_cloud_L44"])
    # per-target (four-model mean) rows for a fold-clustered CI on the disagreement
    per = {}
    for r in rows:
        d = per.setdefault(r["pdb"], {"fold": r["fold"], "pred": [], "meas": []})
        d["pred"].append(r["pred_discounted"]); d["meas"].append(r["delta_cloud"])
    pdbs = sorted(per); folds = ST.pinned_folds(pdbs)
    pred = np.array([np.nanmean(per[p]["pred"]) for p in pdbs]); meas = np.array([np.mean(per[p]["meas"]) for p in pdbs])
    res = ST.compare(pred, meas, folds, names=pdbs, label="ladder prediction (discounted, -2.1496 x gam x cos) minus MEASURED cloud delta, per target (ORACLE)")
    print(ST.fmt(res))
    print(json.dumps(summ, indent=1))
    out = {"label": "own-native models on the S24 ladder (ORACLE DIAGNOSTIC)", "prereg": "s26/PREREG_memorisation_on_the_ladder.md",
           "summary": summ, "disagreement_stats": res, "rows": rows}
    p = W.save("ladder", out, rows=rows, n_expected=len(tg) * 4, complete_keys=("pdb", "model", "gam_prob", "delta_cloud"))
    print("  ->", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
