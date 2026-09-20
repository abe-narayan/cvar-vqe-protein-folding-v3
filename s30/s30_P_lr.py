#!/usr/bin/env python
"""s30/s30_P_lr.py -- S30 lane P, measurement (B): the regression RESTRICTED TO WHERE THE PRIZE IS.

The coordinator's objection, which is correct: the aggregate continuous R^2 of `s30_P_r2.py` is the
wrong statistic in two ways.  It is diluted by the short-range pairs that carry only ~32% of the
prize, and it prices MAGNITUDE while the prize is SIGN (ORACLE-sign 62.1% against ORACLE-magnitude
12.8%).  So: the same question, restricted to |i-j| >= 7, in the pair space where the prize was
localised, out of fold on the pinned folds, with two matched controls.

PRE-REGISTERED BARS, written here before the first number of this file existed:
  * B2's ceiling  R^2 = 1.96 %,  primary 12.82 %,  ambitious 39.44 %  (same as PREREG_S30_P.md)
  * sign accuracy: the bar is the lane's own ladder -- 0.80 buys -0.274 A, and the PLAUSIBLE
    zero-information baseline ("always positive", the leave-fold-out global profile's sign) is
    already 0.627 in the long-range bin, so the excess over THAT is the statistic.

    python s30/s30_P_lr.py
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

OUT = os.path.join(HERE, "results", "s30_P_lr.json")
POOL_K, M, MIN_SEP = 500, 75, 2
FEATS = ["expected", "sd", "entropy", "exp_minus_median", "exp_minus_mode", "sep", "n",
         "d_prod", "pool75_mean", "pool500_mean", "pool75_sd", "exp_minus_pool75",
         "exp_minus_pool500", "rg_disagree", "post_sd_mean", "score_sd", "sep_over_n",
         "expected_over_n", "prob_tail_hi", "prob_tail_lo"]
DIST_ONLY = {"expected", "sd", "entropy", "exp_minus_median", "exp_minus_mode", "post_sd_mean",
             "prob_tail_hi", "prob_tail_lo", "expected_over_n"}


def rows():
    X, Y, S, T, FD = [], [], [], [], []
    for t in I.targets():
        pdb = t["pdb"]; n = int(t["n"])
        u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
        ii, jj = I.pair_index(n, MIN_SEP); sep = (jj - ii).astype(float)
        exp = np.asarray(dg["expected"], float); sd = np.asarray(dg["sd"], float)
        P = np.asarray(dg["prob"], float); cen = np.asarray(dg["centres"], float)
        ent = -(P * np.log(np.clip(P, 1e-12, None))).sum(1)
        med = cen[np.argmax(np.cumsum(P, 1) >= 0.5, axis=1)]
        mode = cen[np.argmax(P, 1)]
        nat = np.asarray(u["nat_ca"], float)
        d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)              # ORACLE
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        sub = np.asarray(rec["sub"], int)
        C0 = np.asarray(rec["avg_ca"], float)
        d_prod = np.linalg.norm(C0[ii] - C0[jj], axis=1)
        rgp = float(np.sqrt((np.sum(exp ** 2) * 2 + (n - 1) * 2 * 3.8 ** 2) / (2.0 * n * n)))
        rgo = float(np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1)).mean())
        sc = I.shipped_score(dg, D)
        col = dict(expected=exp, sd=sd, entropy=ent, exp_minus_median=exp - med,
                   exp_minus_mode=exp - mode, sep=sep, n=np.full(len(exp), n), d_prod=d_prod,
                   pool75_mean=D[sub].mean(0), pool500_mean=D.mean(0), pool75_sd=D[sub].std(0),
                   exp_minus_pool75=exp - D[sub].mean(0), exp_minus_pool500=exp - D.mean(0),
                   rg_disagree=np.full(len(exp), rgp - rgo),
                   post_sd_mean=np.full(len(exp), sd.mean()),
                   score_sd=np.full(len(exp), sc.std()), sep_over_n=sep / n,
                   expected_over_n=exp / n, prob_tail_hi=P[:, -3:].sum(1),
                   prob_tail_lo=P[:, :3].sum(1))
        X.append(np.column_stack([col[f] for f in FEATS]))
        Y.append(exp - d_nat)                                          # ORACLE
        S.append(sep); T.append(np.full(len(exp), len(T))); FD.append(np.full(len(exp), t["fold"]))
    return (np.concatenate(X), np.concatenate(Y), np.concatenate(S),
            np.concatenate(T).astype(int), np.concatenate(FD).astype(int))


def lfo_r2(X, y, fd, tgt, alphas=(1e0, 1e1, 1e2, 1e3, 1e4)):
    yh = np.zeros_like(y)
    F = sorted(set(fd.tolist()))
    for f in F:
        tr, te = fd != f, fd == f
        best, ba = np.inf, alphas[0]
        for a in alphas:
            err = 0.0
            for g in [q for q in F if q != f]:
                itr, ite = tr & (fd != g), fd == g
                mu, s = X[itr].mean(0), np.where(X[itr].std(0) > 1e-12, X[itr].std(0), 1.0)
                A = (X[itr] - mu) / s; Bx = (X[ite] - mu) / s
                A1 = np.column_stack([A, np.ones(len(A))]); B1 = np.column_stack([Bx, np.ones(len(Bx))])
                w = np.linalg.solve(A1.T @ A1 + a * np.eye(A1.shape[1]), A1.T @ y[itr])
                err += float(((y[ite] - B1 @ w) ** 2).sum())
            if err < best:
                best, ba = err, a
        mu, s = X[tr].mean(0), np.where(X[tr].std(0) > 1e-12, X[tr].std(0), 1.0)
        A = (X[tr] - mu) / s; Bx = (X[te] - mu) / s
        A1 = np.column_stack([A, np.ones(len(A))]); B1 = np.column_stack([Bx, np.ones(len(Bx))])
        w = np.linalg.solve(A1.T @ A1 + ba * np.eye(A1.shape[1]), A1.T @ y[tr])
        yh[te] = B1 @ w
    return yh


def block(X, y, fd, tgt, label, rng):
    out = {}
    yh = lfo_r2(X, y, fd, tgt)
    ss = lambda a: float((a ** 2).sum())
    out["R2_oof"] = 1 - ss(y - yh) / ss(y - y.mean())
    out["R2_oof_vs_zero"] = 1 - ss(y - yh) / ss(y)
    # control 1: matched-dimension random features
    Xr = rng.normal(size=X.shape)
    out["R2_ctrl_randfeat"] = 1 - ss(y - lfo_r2(Xr, y, fd, tgt)) / ss(y - y.mean())
    # control 2: features permuted across targets (keeps marginals, breaks the link)
    perm = rng.permutation(len(y))
    out["R2_ctrl_permrows"] = 1 - ss(y - lfo_r2(X[perm], y, fd, tgt)) / ss(y - y.mean())
    out["excess_over_randfeat"] = out["R2_oof"] - out["R2_ctrl_randfeat"]
    # sign accuracy against the plausible zero-information baseline (always predict +)
    out["sign_acc"] = float((np.sign(yh) == np.sign(y)).mean())
    out["sign_acc_baseline_always_plus"] = float((np.sign(y) > 0).mean())
    out["sign_acc_excess"] = out["sign_acc"] - out["sign_acc_baseline_always_plus"]
    out["n_pairs"] = int(len(y)); out["label"] = label
    return out


def main():
    X, y, sep, tgt, fd = rows()
    rng = np.random.default_rng(300204)
    ind = [k for k, f in enumerate(FEATS) if f not in DIST_ONLY]
    masks = {"all_pairs": np.ones(len(y), bool), "LONGRANGE_sep_ge7": sep >= 7,
             "shortrange_sep_le4": sep <= 4, "mid_sep_5_6": (sep >= 5) & (sep <= 6)}
    res = {"bars": {"B2": 0.0196, "primary": 0.1282, "ambitious": 0.3944},
           "features": FEATS, "no_distogram_features": [FEATS[k] for k in ind], "blocks": {}}
    for nm, mk in masks.items():
        res["blocks"][nm] = block(X[mk], y[mk], fd[mk], tgt[mk], nm, rng)
        res["blocks"][nm + "|nodist"] = block(X[mk][:, ind], y[mk], fd[mk], tgt[mk],
                                              nm + "|nodist", rng)
        res["blocks"][nm + "|poolonly"] = block(
            X[mk][:, [FEATS.index(f) for f in ("exp_minus_pool75", "exp_minus_pool500",
                                               "pool75_sd", "rg_disagree", "sep", "n")]],
            y[mk], fd[mk], tgt[mk], nm + "|poolonly", rng)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for k, b in res["blocks"].items():
        print("%-32s n=%6d  R2 %+.4f  randfeat %+.4f  permrows %+.4f  excess %+.4f | "
              "sign %.3f vs always+ %.3f (%+.3f)"
              % (k, b["n_pairs"], b["R2_oof"], b["R2_ctrl_randfeat"], b["R2_ctrl_permrows"],
                 b["excess_over_randfeat"], b["sign_acc"], b["sign_acc_baseline_always_plus"],
                 b["sign_acc_excess"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
