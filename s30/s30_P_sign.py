#!/usr/bin/env python
"""s30/s30_P_sign.py -- S30 lane P, third measurement: CAN ANYTHING SUPPLY THE FIVE SIGNS?

`s30_P_prior.py` reduced the prior route to five signs per target: ORACLE signs with a
leave-fold-out magnitude buy -0.3567 A (point cloud, 2.29x MDE, 5/5 folds), ORACLE magnitudes with
a leave-fold-out sign buy -0.0733 (NOT MEASURED), and the matched-accuracy corruption ladder says
~0.8 sign accuracy is needed for -0.27 A while ~0.6 is worth nothing.

This names and measures every native-free sign channel available, including the distogram's
UNCOLLAPSED posterior -- which the lane's feature regression probed only through two scalars, and
which project memory flags ("do not collapse the posterior", -2.253 A).  Each channel replaces ONLY
the sign in the -0.3567 arm; the magnitude stays leave-fold-out.  So every row is directly
comparable to the ORACLE-sign row and to the accuracy ladder.

    python s30/s30_P_sign.py
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s30.s30_P_prior import bin_of, score_shift, emit, POOL_K, M, MIN_SEP   # noqa: E402

OUT = os.path.join(HERE, "results", "s30_P_sign.json")
CH = ["POSTSKEW", "POSTMODE", "POSTSD", "POOLDIS75", "POOLDIS500", "RGDIS", "LFO_LOGIT"]


def main():
    tg = I.targets()
    R = []
    for t in tg:
        pdb = t["pdb"]; n = int(t["n"])
        u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
        ii, jj = I.pair_index(n, MIN_SEP)
        b = bin_of(jj - ii)
        exp = np.asarray(dg["expected"], float)
        P = np.asarray(dg["prob"], float); cen = np.asarray(dg["centres"], float)
        sd = np.asarray(dg["sd"], float)
        cdf = np.cumsum(P, 1)
        med = cen[np.argmax(cdf >= 0.5, axis=1)]                    # posterior median
        mode = cen[np.argmax(P, 1)]
        nat = np.asarray(u["nat_ca"], float)
        d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)           # ORACLE
        err = exp - d_nat                                           # ORACLE
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        sub = np.asarray(rec["sub"], int)
        rg_pred = float(np.sqrt((np.sum(exp ** 2) * 2 + (n - 1) * 2 * 3.8 ** 2) / (2.0 * n * n)))
        rg_pool = float(np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1)).mean())
        bm = lambda v: np.array([v[b == q].mean() if (b == q).any() else 0.0 for q in range(5)])
        R.append(dict(pdb=pdb, fold=int(t["fold"]), n=n, b=b, exp=exp, D=D, W=W, nat=nat, dg=dg,
                      prof=bm(err),                                             # ORACLE
                      POSTSKEW=bm(exp - med), POSTMODE=bm(exp - mode),
                      POSTSD=bm(sd - sd.mean()),
                      POOLDIS75=bm(exp - D[sub].mean(0)), POOLDIS500=bm(exp - D.mean(0)),
                      RGDIS=np.full(5, rg_pred - rg_pool),
                      pool_mean=float(np.asarray(u["rr"], float)[pool].mean()),
                      pool_best=float(np.asarray(u["rr"], float)[pool].min())))
    folds = np.array([r["fold"] for r in R]); pdbs = [r["pdb"] for r in R]
    F = sorted(set(folds.tolist()))
    PR = np.array([r["prof"] for r in R])
    gmag = {f: np.abs(PR[folds != f]).mean(0) for f in F}

    # LFO logistic-style channel: a leave-fold-out linear model over the six native-free
    # per-bin signals predicting sign(prof).  Trained on distogram+pool features only, so it is
    # the corrector the charter warns about; reported with the matched-accuracy ladder beside it.
    feats = ["POSTSKEW", "POSTMODE", "POSTSD", "POOLDIS75", "POOLDIS500", "RGDIS"]
    X = np.array([[np.concatenate([[r[f][q]] for f in feats] + [[1.0 * (q == k) for k in range(5)]])
                   for q in range(5)] for r in R])                     # (nt, 5, 6+5)
    Y = np.sign(PR)
    lg = np.zeros_like(Y)
    for f in F:
        tr = folds != f
        A = X[tr].reshape(-1, X.shape[2]); yv = Y[tr].reshape(-1)
        mu, s = A.mean(0), np.where(A.std(0) > 1e-12, A.std(0), 1.0)
        A = (A - mu) / s
        w = np.linalg.solve(A.T @ A + 10.0 * np.eye(A.shape[1]), A.T @ yv)
        Bx = (X[folds == f].reshape(-1, X.shape[2]) - mu) / s
        lg[folds == f] = np.sign(Bx @ w).reshape(-1, 5)

    def run(delta_of):
        v = []
        for r in R:
            sc = score_shift(r["dg"], r["D"], -delta_of(r))
            v.append(emit(r["W"], sc, r["nat"])[0])
        return np.array(v)

    base = run(lambda r: np.zeros(len(r["exp"])))
    arms = {}
    acc = {}
    W_ = np.abs(PR)
    for name in CH:
        S = lg if name == "LFO_LOGIT" else np.sign(np.array([r[name] for r in R]))
        agree = (S == Y)
        acc[name] = {"sign_acc_all": float(agree.mean()),
                     "sign_acc_damage_weighted": float((agree * W_).sum() / W_.sum()),
                     "sign_acc_bin4": float(agree[:, 4].mean()),
                     "per_bin": {}}
        # (A) per-bin accuracy against the PLAUSIBLE zero-information baseline -- "always
        # predict the sign of the leave-fold-out global profile", which is positive in every
        # bin, not a coin flip (`zero-information-control-must-be-plausible`).
        for q in range(5):
            bl = (np.array([np.sign(gmag[f][q]) if False else 1.0 for f in folds]) == Y[:, q])
            c = ST.compare(-agree[:, q].astype(float), -bl.astype(float), folds=folds,
                           names=pdbs, label="%s|bin%d" % (name, q))
            acc[name]["per_bin"][q] = {
                "acc": float(agree[:, q].mean()), "baseline_always_plus": float(bl.mean()),
                "excess": float(-c["effect"]), "mde": c["mde"],
                "x_mde": abs(c["effect_over_mde"]), "ci95_fold": [-c["ci95_fold"][1],
                                                                 -c["ci95_fold"][0]],
                "folds_same_sign": c["folds_same_sign"]}
        idx = {r["pdb"]: k for k, r in enumerate(R)}
        arms[name] = run(lambda r, S=S: (S[idx[r["pdb"]]] * gmag[r["fold"]])[r["b"]])
    arms["ORACLEsign"] = run(lambda r: (np.sign(r["prof"]) * gmag[r["fold"]])[r["b"]])
    acc["ORACLEsign"] = {"sign_acc_all": 1.0, "sign_acc_damage_weighted": 1.0,
                         "sign_acc_bin4": 1.0, "per_bin": {}}
    # MANDATORY NULL (`error-coherence-decides-correctors`): the ORACLE sign corrupted at
    # random to each channel's OWN accuracy, same leave-fold-out magnitude.  i.i.d. mistakes
    # at matched accuracy, against the channels' coherent ones.
    rng = np.random.default_rng(300203)
    idx = {r["pdb"]: k for k, r in enumerate(R)}
    for a in (0.50, 0.60, 0.64, 0.70, 0.80):
        Sr = np.sign(PR) * np.where(rng.random(PR.shape) < a, 1.0, -1.0)
        nm = "RANDsign_acc%.2f" % a
        arms[nm] = run(lambda r, S=Sr: (S[idx[r["pdb"]]] * gmag[r["fold"]])[r["b"]])
        ag = (Sr == Y)
        acc[nm] = {"sign_acc_all": float(ag.mean()),
                   "sign_acc_damage_weighted": float((ag * W_).sum() / W_.sum()),
                   "sign_acc_bin4": float(ag[:, 4].mean()), "per_bin": {}}

    nt = len(R)
    pm = np.array([r["pool_mean"] for r in R]); pb = np.array([r["pool_best"] for r in R])
    tA = np.zeros(nt, bool); tA[np.argsort(-pm)[:18]] = True
    tB = np.zeros(nt, bool); tB[np.argsort(-pb)[:18]] = True
    res = {"prod_cloud": float(base.mean()),
           "note": "POINT CLOUD, an INTERMEDIATE.  Magnitude is leave-fold-out in every row; "
                   "only the SIGN channel changes.  Comparable to the accuracy ladder in "
                   "s30_P_prior.json (acc 0.8 -> -0.2744, 0.6 -> -0.0747).",
           "arms": {}}
    for name, v in arms.items():
        c = ST.compare(v, base, folds=folds, names=pdbs, label=name)
        res["arms"][name] = dict(acc[name], mean=float(v.mean()), delta=c["effect"],
                                 se=c["se"], mde=c["mde"], x_mde=c["effect_over_mde"],
                                 ci95_fold=c["ci95_fold"], folds_same_sign=c["folds_same_sign"],
                                 W=c["n_better"], L=c["n_worse"], verdict=c["verdict"],
                                 tailA=float((v - base)[tA].mean()),
                                 tailB=float((v - base)[tB].mean()))
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("production point cloud %.4f" % base.mean())
    for k, a in res["arms"].items():
        print("%-12s acc %.3f (dmg-wt %.3f, bin4 %.3f)  %.4f  delta %+.4f  %.2fxMDE  %dW/%dL  %s"
              % (k, a["sign_acc_all"], a["sign_acc_damage_weighted"], a["sign_acc_bin4"],
                 a["mean"], a["delta"], abs(a["x_mde"]), a["W"], a["L"], a["verdict"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
