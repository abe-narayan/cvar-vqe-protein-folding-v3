#!/usr/bin/env python
"""s31/s31_E3_decomp.py -- S31 lane E: where does the prize live, and where does the corrector spend?

Per target, decompose the ORACLE ideal correction `y = expected - d_nat` orthogonally against the
ORACLE common-mode direction `mu = pool75_mean - d_nat`:

    y = y_along + y_perp,     y_along = <y,mu>/<mu,mu> * mu,     y_perp = y - y_along

Then decompose the fitted corrector's SUM-OF-SQUARES REDUCTION the same way.  Because the split is
orthogonal within a target,

    ||y - yhat||^2  =  ||y_along - yhat_along||^2  +  ||y_perp - yhat_perp||^2

so "which component does the corrector actually predict?" is answered exactly, not by a correlation.

This settles a mechanism question S30's §12 stated the other way round.  Everything here reads
`d_nat`: **ORACLE / NOT DEPLOYABLE**, a diagnostic and a price, never a method.

    python s31/s31_E3_decomp.py
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

from s24 import stats_lib as ST            # noqa: E402
from s31 import s31_E_lib as L             # noqa: E402
from s31.s31_E2_deltas import p_along, p_perp   # noqa: E402

OUT = os.path.join(HERE, "results", "s31_E3_decomp.json")


def main():
    d = L.load()
    X, y, mu, mh = d["X"], d["y"], d["mu"], d["mu_hat"]
    tgt, fold = d["tgt"], d["fold"]
    pdbs = [str(p) for p in d["pdbs"]]
    feats = [str(f) for f in d["feats"]]
    tfold = np.array([fold[tgt == q][0] for q in range(len(pdbs))], int)
    yhat = L.lfo_ridge(X[:, [feats.index(f) for f in L.NEST["N3_plus_pool"]]], y, fold)

    out = {"note": "ORACLE / NOT DEPLOYABLE throughout -- every quantity reads d_nat",
           "directions": {}}
    for dname, u in (("mu_ORACLE", mu), ("mu_hat_NATIVE_FREE", mh)):
        rec = {k: [] for k in ("frac_energy_along", "frac_energy_perp", "sse_red_along",
                               "sse_red_perp", "sse_red_total", "share_of_reduction_along",
                               "share_of_reduction_perp", "yhat_frac_energy_along",
                               "ceiling_ratio_perp_over_full")}
        for q, m in L.per_target(y, tgt):
            yy, hh, uu = y[m], yhat[m], u[m]
            ya, yp = p_along(yy, uu), p_perp(yy, uu)
            ha, hp = p_along(hh, uu), p_perp(hh, uu)
            e = float(yy @ yy)
            ra = float(ya @ ya - (ya - ha) @ (ya - ha))
            rp = float(yp @ yp - (yp - hp) @ (yp - hp))
            rec["frac_energy_along"].append(float(ya @ ya) / e)
            rec["frac_energy_perp"].append(float(yp @ yp) / e)
            rec["sse_red_along"].append(ra); rec["sse_red_perp"].append(rp)
            rec["sse_red_total"].append(ra + rp)
            tot = ra + rp
            rec["share_of_reduction_along"].append(ra / tot if abs(tot) > 1e-12 else np.nan)
            rec["share_of_reduction_perp"].append(rp / tot if abs(tot) > 1e-12 else np.nan)
            nh = float(hh @ hh)
            rec["yhat_frac_energy_along"].append(float(ha @ ha) / nh if nh > 1e-12 else np.nan)
            rec["ceiling_ratio_perp_over_full"].append(float(yp @ yp) / e)
        V = {k: np.array(v, float) for k, v in rec.items()}
        blk = {}
        for k, a in V.items():
            c = ST.compare(a, np.zeros(len(a)), folds=tfold, names=pdbs, label="E3|%s|%s" % (dname, k))
            blk[k] = {"mean": float(np.nanmean(a)), "median": float(np.nanmedian(a)),
                      "ci95_fold": c["ci95_fold"], "se": c["se"]}
        # pooled (energy-weighted) versions -- the honest aggregate for an energy share
        ea = np.array(rec["frac_energy_along"]); ep = np.array(rec["frac_energy_perp"])
        tote = np.array([float(y[tgt == q] @ y[tgt == q]) for q in range(len(pdbs))])
        blk["POOLED_frac_energy_along"] = float((ea * tote).sum() / tote.sum())
        blk["POOLED_frac_energy_perp"] = float((ep * tote).sum() / tote.sum())
        blk["POOLED_share_of_reduction_along"] = float(np.nansum(V["sse_red_along"])
                                                       / np.nansum(V["sse_red_total"]))
        blk["POOLED_share_of_reduction_perp"] = float(np.nansum(V["sse_red_perp"])
                                                      / np.nansum(V["sse_red_total"]))
        blk["per_target"] = {pdbs[i]: {k: float(V[k][i]) for k in V} for i in range(len(pdbs))}
        out["directions"][dname] = blk

    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, OUT)

    for dn, b in out["directions"].items():
        print("=== direction %s   (ORACLE / NOT DEPLOYABLE)" % dn)
        for k in ("frac_energy_along", "frac_energy_perp", "share_of_reduction_along",
                  "share_of_reduction_perp", "yhat_frac_energy_along"):
            q = b[k]
            print("  %-28s mean %+.4f  median %+.4f  fold95 [%+.4f,%+.4f]"
                  % (k, q["mean"], q["median"], q["ci95_fold"][0], q["ci95_fold"][1]))
        print("  POOLED energy along %.4f / perp %.4f | POOLED reduction along %.4f / perp %.4f"
              % (b["POOLED_frac_energy_along"], b["POOLED_frac_energy_perp"],
                 b["POOLED_share_of_reduction_along"], b["POOLED_share_of_reduction_perp"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
