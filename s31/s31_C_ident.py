#!/usr/bin/env python
"""s31/s31_C_ident.py -- lane C, F-C1c (identifiability of the ORACLE readout weights from
native-free features) and the CROSS-TERM ACCOUNTING that prices the whole native-free class.

The exact decomposition (see `s31_C_ladder.py`): for any `sum(w) = 1`,

    ||X(w) - nat||^2 / n  =  c2  +  2<c, D w>/n  +  ||D w||^2 / n

`c2` is fixed by the candidate set and NO readout can touch it.  A readout beats the set mean only
through the CROSS TERM.  So "can a native-free readout work" is exactly "can a native-free rule
produce a negative cross term", and that is a scalar this file measures against the ORACLE's.

F-C1c as registered: fit the ORACLE weight vector from native-free per-candidate features OUT OF
FOLD on the pinned folds; the bar is 5% out-of-fold R^2 and the fitted weights are APPLIED and
measured (contract rule 21 -- never quote an implied conversion alone).

    python s31/s31_C_ident.py run
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                    # noqa: E402
from s24 import stats_lib as ST                    # noqa: E402
from s31 import s31_C_cache as CA                  # noqa: E402
from s31 import s31_C_ladder as LD                 # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "s31_C_ident.json")
TOP = LD.TOP
FEATS = ("dis_z", "rank_z", "typ_z", "medoid_z", "rg_z", "bond_z", "dis_gap", "const")


def features_for(W, idx, dis, n):
    """NATIVE-FREE per-candidate features over the fixed top-`TOP` set.  Nothing here reads the
    native; every one is computable at inference."""
    Wt = np.asarray(W, float)[idx]                                   # (TOP, n, 3)
    d = np.asarray(dis, float)[idx]
    P = I.pairwise_rmsd(Wt)
    typ = P.mean(1)                                                  # typicality (mean pairwise)
    med = int(np.argmin(typ))
    dmed = P[med]
    rg = np.sqrt(((Wt - Wt.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
    bond = np.linalg.norm(Wt[:, 1:] - Wt[:, :-1], axis=-1).mean(-1)
    rank = np.argsort(np.argsort(d)).astype(float)

    def z(x):
        x = np.asarray(x, float)
        return (x - x.mean()) / max(float(x.std()), 1e-12)

    gap = d - d.min()
    return np.vstack([z(d), z(rank), z(typ), z(dmed), z(rg), z(bond), z(gap),
                      np.ones(len(idx))]).T                          # (TOP, 8)


def main():
    pdbs = CA.all_pdbs()
    t0 = time.time()
    Xs, Ys, F, meta = [], [], [], []
    per = {}
    for i, p in enumerate(pdbs):
        from s27 import s28_A_amp as A28
        z = CA.load(p)
        W = np.asarray(z["W"], float)
        order = np.asarray(z["order"], int)
        dis = np.asarray(z["dis"], float)
        nat = np.asarray(z["nat"], float)                            # ORACLE
        top75 = np.asarray(z["top75"], int)
        n, fold = int(z["n"]), int(z["fold"])
        frame = A28.Frame(W, np.sort(top75))
        Wf = frame.Wf
        natp = I.superpose_batch(nat[None], frame.ref)[0].ravel()    # ORACLE
        idx = order[:TOP]
        Wt = Wf[idx]
        A = Wt.T
        w_or = LD.oracle_convex(A, natp)                             # ORACLE target of the fit
        Xs.append(features_for(W, idx, dis, n))
        Ys.append(w_or - 1.0 / TOP)                                  # the zero-sum part
        F.append(fold)
        meta.append(dict(pdb=p, n=n, fold=fold))
        per[p] = dict(Wt=Wt, natp=natp, n=n, nat=nat, w_or=w_or)
        if i % 20 == 0:
            print("[%3d/%3d] %s %.1fs" % (i + 1, len(pdbs), p, time.time() - t0), flush=True)

    F = np.array(F)
    folds = sorted(set(F.tolist()))
    out = {"n": len(pdbs), "top": TOP, "feats": list(FEATS),
           "provenance": ST.provenance(__file__)}

    # ---------------------------------------------------------------- F-C1c: out-of-fold fit
    lam = 1e-3
    oof_pred = [None] * len(pdbs)
    for q in folds:
        tr = [j for j in range(len(pdbs)) if F[j] != q]
        te = [j for j in range(len(pdbs)) if F[j] == q]
        Xtr = np.vstack([Xs[j] for j in tr]); Ytr = np.concatenate([Ys[j] for j in tr])
        G = Xtr.T @ Xtr + lam * len(Xtr) * np.eye(Xtr.shape[1])
        beta = np.linalg.solve(G, Xtr.T @ Ytr)
        for j in te:
            oof_pred[j] = Xs[j] @ beta

    yall = np.concatenate(Ys); pall = np.concatenate(oof_pred)
    ss_res = float(((yall - pall) ** 2).sum())
    ss_tot = float(((yall - yall.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot
    out["F_C1c"] = {"oof_r2": r2, "bar": 0.05, "fires": bool(r2 >= 0.05),
                    "corr": float(np.corrcoef(yall, pall)[0, 1]),
                    "note": "target is the ORACLE convex weight minus uniform; features are "
                            "native-free per-candidate quantities; leave-fold-out on pinned folds"}

    # ---- APPLY the fitted weights (rule 21) ------------------------------------------------
    r_fit, r_prod, r_or = [], [], []
    for j, m in enumerate(meta):
        d = per[m["pdb"]]
        w = 1.0 / TOP + oof_pred[j]
        w = np.clip(w, 0.0, None)
        w = w / max(w.sum(), 1e-12)                                  # renormalise: still sum-1
        r_fit.append(float(I.ca_rmsd(LD.emit(d["Wt"], w, d["n"]), d["nat"])))
        wu = np.zeros(TOP); wu[:75] = 1.0 / 75
        r_prod.append(float(I.ca_rmsd(LD.emit(d["Wt"], wu, d["n"]), d["nat"])))
        r_or.append(float(I.ca_rmsd(LD.emit(d["Wt"], d["w_or"], d["n"]), d["nat"])))
    r_fit = np.array(r_fit); r_prod = np.array(r_prod); r_or = np.array(r_or)
    out["applied"] = ST.compare(r_fit, r_prod, folds=F, names=pdbs,
                                label="oof_fitted_weights_vs_prefix75",
                                seed_parts=("s31C", "ident"))
    out["applied"]["mean_oracle_convex"] = float(r_or.mean())

    # ---------------------------------------------------------------- cross-term accounting
    rows = {json.loads(l)["pdb"]: json.loads(l) for l in open(LD.ROWS)}
    rows = [rows[p] for p in pdbs if p in rows]

    def cr(a):
        return np.array([r["arms"][a]["cross"] for r in rows if "cross" in r["arms"][a]])

    nf = {a: float(cr(a).mean()) for a in
          ("unif_prefix_m75", "unif_prefix_m128", "soft_T0.50", "soft_T1.00", "soft_T2.00",
           "soft_T4.00", "typ_T1.00", "typ_T2.00") if a in rows[0]["arms"]}
    orc = {a: float(cr(a).mean()) for a in
           ("ORACLE_argmin128", "ORACLE_convex128", "ORACLE_affine128")}
    best_nf = min(nf, key=lambda a: nf[a])
    out["cross_accounting"] = {
        "native_free": nf, "oracle": orc,
        "best_native_free_arm": best_nf,
        "capture_of_convex_oracle": float(nf[best_nf] / orc["ORACLE_convex128"]),
        "capture_of_argmin_oracle": float(nf[best_nf] / orc["ORACLE_argmin128"]),
        "note": "the cross term 2<c,Dw>/n is the ONLY channel a sum-to-one readout has; c2 is "
                "fixed by the candidate set.  Units are A^2.",
    }

    # what fraction of the ORACLE cross would be needed to reach 3.00 / 2.50 A on the mean
    c2 = np.array([rows[i]["arms"]["unif_prefix_m128"]["c2"] for i in range(len(rows))])
    d2_or = np.array([rows[i]["arms"]["ORACLE_convex128"]["d2"] for i in range(len(rows))])
    cross_or = cr("ORACLE_convex128")
    fr = np.arange(0.0, 1.001, 0.02)
    curve = []
    for f in fr:
        v = c2 + f * cross_or + (f ** 2) * d2_or                     # scale the ORACLE deviation
        curve.append(float(np.sqrt(np.maximum(v, 0.0)).mean()))
    out["what_would_it_take"] = {
        "fraction_grid": [float(x) for x in fr], "mean_rmsd": curve,
        "model": "w(f) = uniform + f*(w_oracle_convex - uniform); f=0 is the 128-mean, f=1 the "
                 "ORACLE convex optimum.  ORACLE / NOT DEPLOYABLE -- it prices the axis, it is "
                 "not an operator.",
        "f_for_3.00": next((float(fr[i]) for i, v in enumerate(curve) if v <= 3.00), None),
        "f_for_2.50": next((float(fr[i]) for i, v in enumerate(curve) if v <= 2.50), None),
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    print("\nF-C1c -- can the ORACLE readout weights be predicted from native-free features?")
    print("  out-of-fold R2 %.5f against a registered 5%% bar  ->  %s" %
          (out["F_C1c"]["oof_r2"], "FIRES" if out["F_C1c"]["fires"] else "REFUTED"))
    print("  applied: fitted %.4f vs score-prefix-75 %.4f  (%+.4f, %.2fx MDE, %dW/%dL); "
          "ORACLE convex on the same set %.4f" %
          (out["applied"]["mean_a"], out["applied"]["mean_b"], out["applied"]["effect"],
           out["applied"]["effect_over_mde"], out["applied"]["n_better"],
           out["applied"]["n_worse"], out["applied"]["mean_oracle_convex"]))

    print("\nCROSS-TERM ACCOUNTING (A^2; the only channel a sum-to-one readout has)")
    for a, v in sorted(nf.items(), key=lambda kv: kv[1]):
        print("  native-free  %-18s %+8.4f" % (a, v))
    for a, v in orc.items():
        print("  ORACLE       %-18s %+8.4f" % (a, v))
    print("  best native-free captures %.3f%% of the ORACLE convex cross term" %
          (100 * out["cross_accounting"]["capture_of_convex_oracle"]))

    ww = out["what_would_it_take"]
    print("\nWHAT WOULD IT TAKE (ORACLE / NOT DEPLOYABLE; f = fraction of the ORACLE convex step)")
    print("  f:    " + " ".join("%6.2f" % x for x in ww["fraction_grid"][::5]))
    print("  rmsd: " + " ".join("%6.4f" % x for x in ww["mean_rmsd"][::5]))
    print("  f needed for 3.00 A: %s   for 2.50 A: %s" % (ww["f_for_3.00"], ww["f_for_2.50"]))
    print("\n  wrote %s" % OUT)


if __name__ == "__main__":
    main()
