#!/usr/bin/env python
"""s31/s31_E6_selection.py -- S31 lane E: does a correction act through SELECTION or through the
average?

The endpoint arms say what each correction is worth.  This says HOW.  For each arm, per target:

    score_m = shipped Bayes-risk score of pool member m against the CORRECTED posterior
    rho     = Spearman(score, the member's own true CA RMSD)        <- selection skill, ORACLE
    top75   = mean true RMSD of the 75 members the corrected score retains
    recall  = how many of the ORACLE-best 75 the corrected score retains

`u["rr"]` is the stored per-member RMSD, so nothing is recomputed.  Cheap: no projection.
Every quantity here reads the native: **ORACLE / NOT DEPLOYABLE**, a mechanism diagnostic.

    python s31/s31_E6_selection.py [--deltas <npz>]
"""
from __future__ import annotations

import argparse
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

from s12 import instrument as I                                  # noqa: E402
from s24 import stats_lib as ST                                  # noqa: E402
from s30.s30_P_prior import score_shift, POOL_K, M, MIN_SEP      # noqa: E402

DELTAS = os.path.join(HERE, "results", "s31_E2_deltas.npz")
OUT = os.path.join(HERE, "results", "s31_E6_selection.json")


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra @ ra) * (rb @ rb))
    return float(ra @ rb / den) if den > 1e-12 else np.nan


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--deltas", default=DELTAS)
    a = ap.parse_args()
    z = np.load(a.deltas, allow_pickle=False)
    off = z["off"]; names = [str(s) for s in z["names"]]
    D_ = {k: z["delta_" + k] for k in names}
    tg = I.targets()
    rec = {k: {"rho": [], "top75": [], "recall": []} for k in names}
    folds, pdbs = [], []
    for k, t in enumerate(tg):
        pdb = t["pdb"]; n = int(t["n"])
        u = I.load_univ(pdb); dg = I.distogram(pdb)
        ii, jj = I.pair_index(n, MIN_SEP)
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        rr = np.asarray(u["rr"], float)[pool]                    # ORACLE per-member RMSD
        Dm = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        best75 = set(np.argsort(rr, kind="stable")[:M].tolist())
        s, e = int(off[k]), int(off[k + 1])
        for nm in names:
            sc = score_shift(dg, Dm, -D_[nm][s:e])
            top = np.argsort(sc, kind="stable")[:M]
            rec[nm]["rho"].append(spearman(sc, rr))
            rec[nm]["top75"].append(float(rr[top].mean()))
            rec[nm]["recall"].append(len(best75 & set(top.tolist())) / float(M))
        folds.append(int(t["fold"])); pdbs.append(pdb)
    folds = np.array(folds)
    base = {k: np.array(rec["PROD"][k], float) for k in ("rho", "top75", "recall")}
    out = {"status": "ORACLE / NOT DEPLOYABLE -- every quantity reads the native",
           "n": len(pdbs), "arms": {}}
    for nm in sorted(names):
        d = {}
        for q in ("rho", "top75", "recall"):
            v = np.array(rec[nm][q], float)
            c = ST.compare(v, base[q], folds=folds, names=pdbs, label="E6|%s|%s" % (nm, q))
            d[q] = {"mean": float(v.mean()), "delta_vs_PROD": c["effect"], "mde": c["mde"],
                    "x_mde": c["effect_over_mde"], "ci95_fold": c["ci95_fold"],
                    "folds_same_sign": c["folds_same_sign"]}
        out["arms"][nm] = d
    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, OUT)
    print("%-22s %8s %8s | %8s %8s | %8s" % ("arm", "rho", "d_rho", "top75", "d_top75", "recall"))
    for nm in sorted(names):
        d = out["arms"][nm]
        print("%-22s %+8.4f %+8.4f | %8.4f %+8.4f | %8.4f"
              % (nm, d["rho"]["mean"], d["rho"]["delta_vs_PROD"], d["top75"]["mean"],
                 d["top75"]["delta_vs_PROD"], d["recall"]["mean"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
