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


# NESTED blocks.  The response is `expected - d_nat` and several features contain `expected`,
# so a raw R^2 carries a SHARED-REFERENT FLOOR (`shared-referent-floor`) plus the distogram's own
# calibration/shrinkage, which the record has already closed as a route
# (`error-shape-not-mae-decides-ranking`: slope +0.376, "post-hoc correction is closed"; S25-L2:
# calibrating the posterior makes RMSD worse).  Only the INCREMENT over N1 can be new information.
NEST = {
    "N0_separation_prior": ["sep", "n", "sep_over_n"],
    "N1_plus_calibration": ["sep", "n", "sep_over_n", "expected", "expected_over_n"],
    "N2_plus_posterior_shape": ["sep", "n", "sep_over_n", "expected", "expected_over_n", "sd",
                                "entropy", "exp_minus_median", "exp_minus_mode", "prob_tail_hi",
                                "prob_tail_lo", "post_sd_mean"],
    "N3_plus_pool": FEATS,
    "N1_plus_pool_only": ["sep", "n", "sep_over_n", "expected", "expected_over_n", "pool75_mean",
                          "pool500_mean", "pool75_sd", "exp_minus_pool75", "exp_minus_pool500",
                          "rg_disagree", "d_prod", "score_sd"],
}


def main():
    X, y, sep, tgt, fd = rows()
    rng = np.random.default_rng(300204)
    ind = [k for k, f in enumerate(FEATS) if f not in DIST_ONLY]
    masks = {"all_pairs": np.ones(len(y), bool), "LONGRANGE_sep_ge7": sep >= 7,
             "shortrange_sep_le4": sep <= 4, "mid_sep_5_6": (sep >= 5) & (sep <= 6)}
    res = {"bars": {"B2": 0.0196, "primary": 0.1282, "ambitious": 0.3944},
           "features": FEATS, "no_distogram_features": [FEATS[k] for k in ind],
           "blocks": {}, "nested": {}}
    for nm, mk in masks.items():
        prev = 0.0
        for bn, fl in NEST.items():
            cols = [FEATS.index(f) for f in fl]
            b = block(X[mk][:, cols], y[mk], fd[mk], tgt[mk], nm + "|" + bn, rng)
            b["increment_over_previous"] = b["R2_oof"] - prev
            if bn != "N1_plus_pool_only":
                prev = b["R2_oof"]
            res["nested"].setdefault(nm, {})[bn] = b
        res["nested"][nm]["N1_plus_pool_only"]["increment_over_N1"] = (
            res["nested"][nm]["N1_plus_pool_only"]["R2_oof"]
            - res["nested"][nm]["N1_plus_calibration"]["R2_oof"])
        res["nested"][nm]["N3_plus_pool"]["increment_over_N1"] = (
            res["nested"][nm]["N3_plus_pool"]["R2_oof"]
            - res["nested"][nm]["N1_plus_calibration"]["R2_oof"])
        res["blocks"][nm] = block(X[mk], y[mk], fd[mk], tgt[mk], nm, rng)
        res["blocks"][nm + "|nodist"] = block(X[mk][:, ind], y[mk], fd[mk], tgt[mk],
                                              nm + "|nodist", rng)
        res["blocks"][nm + "|poolonly"] = block(
            X[mk][:, [FEATS.index(f) for f in ("exp_minus_pool75", "exp_minus_pool500",
                                               "pool75_sd", "rg_disagree", "sep", "n")]],
            y[mk], fd[mk], tgt[mk], nm + "|poolonly", rng)
    # ---- THE APPLIED ARM: take the leave-fold-out per-pair prediction as the correction and
    #      measure the EMITTED structure.  An R^2 that does not survive being applied is a
    #      statistic about the prior, not a route to the endpoint.
    from s30.s30_P_prior import score_shift, emit
    res["applied"] = {}
    preds = {}
    for bn in ("N1_plus_calibration", "N2_plus_posterior_shape", "N3_plus_pool"):
        preds[bn] = lfo_r2(X[:, [FEATS.index(f) for f in NEST[bn]]], y, fd, tgt)
    # MANDATORY NULL (`error-coherence-decides-correctors`): a corrector with the SAME out-of-fold
    # R^2 whose mistakes are i.i.d. rather than inherited from the pool.  Same R^2, same applied
    # pipeline; if the fitted one hurts and this one helps, the discriminator is COHERENCE.
    for bn in ("N3_plus_pool", "N1_plus_calibration"):
        r2 = float(1 - ((y - preds[bn]) ** 2).sum() / ((y - y.mean()) ** 2).sum())
        r2 = max(min(r2, 0.999), 0.0)
        sy = float(y.std())
        preds["IIDmatched_R2_%.3f" % r2] = (r2 * (y - y.mean())
                                            + np.sqrt(r2 * (1 - r2)) * sy
                                            * rng.standard_normal(len(y)) + y.mean())
    # ---- THE INCOHERENCE STATISTIC, so "incoherent" is a number a candidate channel can be
    #      tested against rather than an adjective.  mu_p = the POOL's COMMON-MODE pair error
    #      (mean over production's 75 of d_m,p - d_nat,p) -- the quantity S30-L7 proved
    #      non-identifiable from within the pool.  A corrector's RESIDUAL r = y - yhat is
    #      COHERENT to the extent it correlates with mu.
    mu = X[:, FEATS.index("pool75_mean")] - X[:, FEATS.index("expected")] + y

    def coh(r):
        out = []
        for q in sorted(set(tgt.tolist())):
            m = tgt == q
            if r[m].std() > 1e-12 and mu[m].std() > 1e-12:
                out.append(np.corrcoef(r[m], mu[m])[0, 1])
        return float(np.mean(out)), float(np.corrcoef(r, mu)[0, 1])

    res["incoherence"] = {}
    res["incoherence"]["UNCORRECTED_production"] = dict(
        zip(("coh_within_target", "coh_pooled"), coh(y)))
    for bn, yh in preds.items():
        res["incoherence"][bn] = dict(zip(("coh_within_target", "coh_pooled"), coh(y - yh)),
                                      R2_oof=float(1 - ((y - yh) ** 2).sum()
                                                   / ((y - y.mean()) ** 2).sum()))
    for bn, yh in preds.items():
        off, vals, base = 0, [], []
        for t in I.targets():
            pdb = t["pdb"]; n = int(t["n"])
            u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
            ii, jj = I.pair_index(n, MIN_SEP); m = len(ii)
            pool = np.asarray(u["order"], int)[:POOL_K]
            W = np.asarray(u["W"], float)[pool]
            nat = np.asarray(u["nat_ca"], float)
            D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
            vals.append(emit(W, score_shift(dg, D, -yh[off:off + m]), nat)[0])
            base.append(emit(W, score_shift(dg, D, np.zeros(m)), nat)[0])
            off += m
        c = ST.compare(np.array(vals), np.array(base),
                       folds=np.array([t["fold"] for t in I.targets()]),
                       names=[t["pdb"] for t in I.targets()], label="applied|" + bn)
        res["applied"][bn] = {"mean_cloud": float(np.mean(vals)), "prod": float(np.mean(base)),
                              "delta": c["effect"], "mde": c["mde"],
                              "x_mde": c["effect_over_mde"], "ci95_fold": c["ci95_fold"],
                              "W": c["n_better"], "L": c["n_worse"], "verdict": c["verdict"]}
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for nm, d in res["nested"].items():
        print("--", nm)
        for bn, b in d.items():
            print("   %-26s R2 %+.4f  increment %+.4f  sign %.3f (always+ %.3f)"
                  % (bn, b["R2_oof"], b.get("increment_over_N1", b["increment_over_previous"]),
                     b["sign_acc"], b["sign_acc_baseline_always_plus"]))
    for bn, a in res["applied"].items():
        c = res["incoherence"].get(bn, {})
        print("APPLIED %-26s %.4f vs prod %.4f  delta %+.4f  %.2fxMDE  %dW/%dL  coh %+.4f  %s"
              % (bn, a["mean_cloud"], a["prod"], a["delta"], abs(a["x_mde"]), a["W"], a["L"],
                 c.get("coh_within_target", float("nan")), a["verdict"]))
    print("coherence of the UNCORRECTED error: %+.4f within target, %+.4f pooled"
          % (res["incoherence"]["UNCORRECTED_production"]["coh_within_target"],
             res["incoherence"]["UNCORRECTED_production"]["coh_pooled"]))
    for k, b in res["blocks"].items():
        print("%-32s n=%6d  R2 %+.4f  randfeat %+.4f  permrows %+.4f  excess %+.4f | "
              "sign %.3f vs always+ %.3f (%+.3f)"
              % (k, b["n_pairs"], b["R2_oof"], b["R2_ctrl_randfeat"], b["R2_ctrl_permrows"],
                 b["excess_over_randfeat"], b["sign_acc"], b["sign_acc_baseline_always_plus"],
                 b["sign_acc_excess"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
