#!/usr/bin/env python
"""s30/s30_F_score.py -- S30 lane F, experiment F3: what is the distogram confidently wrong
about on the pools where its own filter inverts?

Pre-registration `s30/PREREG_S30_F3.md` (committed eb719398, before this file existed).

Measures the shipped Bayes-risk score's PER-TARGET IN-POOL RANKING SKILL over all 500 BLOSUM
pool members, and asks what predicts its collapse.  Builds NO router and NO detector: every
predictor is ORACLE by design, so the question is "what is the score wrong about", not "can we
detect it".  Correlations with the seven closed router feature families are reported regardless.

REPRODUCTION GATE: the recomputed score's top-75 must equal the production record's `sub`
set-wise on all 126 targets.  Carries the S29-L19 caveat that deep order below the top-75 cut
may differ ~2-in-126 from a fresh posterior recomputation.

    python s30/s30_F_score.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s30_F_score.json")
POOL_K, M, MIN_SEP = 500, 75, 2
N_NULL, SEED_NULL = 20000, 30002


def rg(X):
    X = np.asarray(X, float)
    return float(np.sqrt(((X - X.mean(-2, keepdims=True)) ** 2).sum(-1).mean(-1)))


def rows():
    out = []
    for t in I.targets():
        pdb = t["pdb"]
        u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
        pool = np.asarray(u["order"], int)[:POOL_K]
        W = np.asarray(u["W"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        ii, jj = I.pair_index(int(u["n"]), min_sep=MIN_SEP)
        assert np.array_equal(np.asarray(dg["i"]), ii) and np.array_equal(np.asarray(dg["j"]), jj), pdb
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)              # (500, npairs)
        sc = I.shipped_score(dg, D)
        order = np.argsort(sc, kind="stable")
        sub = np.asarray(rec["sub"], int)
        if set(order[:M].tolist()) != set(sub.tolist()):
            raise RuntimeError("REPRODUCTION GATE FAILED on %s" % pdb)
        rr = np.asarray(u["rr"], float)[pool]

        # ---- the headline: in-pool ranking skill (positive = lower score means lower RMSD)
        rho_pool = float(spearmanr(sc, rr).statistic)
        top = order[:M]
        rho_top75 = float(spearmanr(sc[top], rr[top]).statistic)

        # ---- ORACLE: the distogram's own error against the NATIVE's true distances
        d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)
        exp = np.asarray(dg["expected"], float)
        e = exp - d_nat
        scale_add = float(e.mean())                        # signed additive offset  = SCALE error
        shape_err = float(np.abs(e - e.mean()).mean())     # residual after removing it = SHAPE error
        err_mae = float(np.abs(e).mean())
        scale_mult = float(exp.mean() / d_nat.mean())      # multiplicative form (lane L's framing)

        # ---- confidence (posterior concentration) -- F3 family, diagnostic only
        P = np.asarray(dg["prob"], float)
        ent = float((-(P * np.log(np.clip(P, 1e-12, None))).sum(1)).mean())
        sd_post = float(np.asarray(dg["sd"], float).mean())

        # ---- closed-family features, for the mandatory correlation table
        rgs = np.array([rg(w) for w in W])
        out.append(dict(
            pdb=pdb, n=int(t["n"]), fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
            rho_pool=rho_pool, rho_top75=rho_top75,
            err_mae=err_mae, scale_add=scale_add, abs_scale_add=abs(scale_add),
            shape_err=shape_err, scale_mult=scale_mult, abs_log_scale=abs(float(np.log(scale_mult))),
            post_entropy=ent, post_sd=sd_post,
            # F1 length is `n`;  F2 score distribution;  F4 pool geometry;  F5 rg disagreement
            F2_score_sd=float(sc.std()), F2_score_gap=float(sc[order[M]] - sc[order[0]]),
            F4_pool_rg_sd=float(rgs.std()), F4_top75_rg_sd=float(rgs[top].std()),
            F5_rg_disagree=float(exp.mean() - D.mean()),   # predicted vs REALISED (native-free: F5)
            pool_mean=float(rr.mean()), pool_best=float(rr.min()),
            set_mean_benefit=float(rr.mean() - rr[top].mean()),   # the S30-L3 mechanism quantity
            prod=float(rec["rmsd_avg"]),
        ))
    return out


def resid(y, x):
    """y with x regressed out (both ranked), for the length-residualised rows."""
    from scipy.stats import rankdata
    ry, rx = rankdata(y), rankdata(x)
    b = np.polyfit(rx, ry, 1)
    return ry - np.polyval(b, rx)


def main():
    R = rows()
    g = lambda k: np.array([r[k] for r in R], float)
    pdbs = [r["pdb"] for r in R]; folds = g("fold")
    isf = np.array([r["fail18"] for r in R])
    pm, bp = g("pool_mean"), g("pool_best")
    a1 = np.zeros(126, bool); a1[np.argsort(-pm)[:18]] = True
    a2 = np.zeros(126, bool); a2[np.argsort(-bp)[:18]] = True
    strata = dict(all=np.ones(126, bool), other108=~isf, FAIL18=isf,
                  worst18_poolmean=a1, worst18_bestpool=a2)
    rp, rt = g("rho_pool"), g("rho_top75")

    out = dict(prereg="s30/PREREG_S30_F3.md", n=len(R),
               note="ALL PREDICTORS ORACLE; no router, no detector, no deployable object")

    # -------- F3a: does in-pool skill collapse on the tail?
    out["F3a"] = {nm: dict(rho_pool=float(rp[s].mean()), median=float(np.median(rp[s])),
                           rho_top75=float(rt[s].mean()), n=int(s.sum()))
                  for nm, s in strata.items()}
    c = ST.compare(rp, np.zeros(126), folds=folds, names=pdbs, label="rho_pool_vs_zero",
                   seed_parts=("s30F3",))
    co = ST.compare(rp[~isf], np.zeros(108), folds=folds[~isf],
                    names=[p for p, q in zip(pdbs, ~isf) if q], label="rho_pool_108", seed_parts=("s30F3",))
    cf = ST.compare(rp[isf], np.zeros(18), folds=folds[isf],
                    names=[p for p, q in zip(pdbs, isf) if q], label="rho_pool_F18", seed_parts=("s30F3",))
    rng = np.random.default_rng(SEED_NULL); idx = np.arange(126)
    nl = np.array([rp[rng.choice(idx, 18, replace=False)].mean() for _ in range(N_NULL)])
    obs = float(rp[isf].mean())
    out["F3a"]["stats"] = dict(
        all=dict(mean=c["mean_a"], ci95_fold=c["ci95_fold"], folds_same_sign=c["folds_same_sign"]),
        other108=dict(mean=co["mean_a"], ci95_fold=co["ci95_fold"], folds_same_sign=co["folds_same_sign"]),
        FAIL18=dict(mean=cf["mean_a"], ci95_fold=cf["ci95_fold"], folds_same_sign=cf["folds_same_sign"]),
        random18_null=dict(null_mean=float(nl.mean()),
                           ci95=[float(np.percentile(nl, 2.5)), float(np.percentile(nl, 97.5))],
                           p_two_sided=float(2 * min((nl >= obs).mean(), (nl <= obs).mean()))))
    out["F3a"]["fires"] = bool(co["ci95_fold"][0] > 0 and obs <= 0.5 * rp[~isf].mean()
                               and (obs < out["F3a"]["stats"]["random18_null"]["ci95"][0]
                                    or obs > out["F3a"]["stats"]["random18_null"]["ci95"][1]))

    # -------- F3b/F3c/F3d: what predicts the collapse?
    preds = ("err_mae", "abs_scale_add", "shape_err", "abs_log_scale", "post_entropy", "post_sd",
             "n", "F2_score_sd", "F2_score_gap", "F4_pool_rg_sd", "F4_top75_rg_sd", "F5_rg_disagree")
    tab = {}
    for k in preds:
        x = g(k)
        rs = spearmanr(x, rp)
        rr_ = spearmanr(resid(rp, g("n")), resid(x, g("n")))    # length-residualised
        # fold-clustered CI on the Spearman
        rngf = np.random.default_rng(7)
        F = sorted(set(folds.tolist()))
        bs = []
        for _ in range(2000):
            sel = np.concatenate([np.where(folds == q)[0] for q in rngf.choice(F, len(F), replace=True)])
            bs.append(spearmanr(x[sel], rp[sel]).statistic)
        tab[k] = dict(spearman=float(rs.statistic), p=float(rs.pvalue),
                      spearman_resid_n=float(rr_.statistic),
                      ci95_fold=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                      mean_all=float(x.mean()), mean_F18=float(x[isf].mean()),
                      mean_108=float(x[~isf].mean()))
    out["predictors_of_rho_pool"] = tab
    out["F3b"] = dict(stat=tab["err_mae"]["spearman"], bar=-0.40,
                      fires=bool(tab["err_mae"]["spearman"] <= -0.40
                                 and tab["err_mae"]["ci95_fold"][1] < 0))

    # F3c: partial Spearman of rho_pool with SCALE error controlling for SHAPE error
    pr = spearmanr(resid(rp, g("shape_err")), resid(g("abs_scale_add"), g("shape_err")))
    pr2 = spearmanr(resid(rp, g("abs_scale_add")), resid(g("shape_err"), g("abs_scale_add")))
    out["F3c"] = dict(partial_scale_given_shape=float(pr.statistic), p=float(pr.pvalue),
                      partial_shape_given_scale=float(pr2.statistic), p2=float(pr2.pvalue),
                      bar=-0.30, fires=bool(pr.statistic <= -0.30),
                      caveat="statistic guessed by lane F; lane L was unreachable to specify it")

    # F3d: does confidence*error beat error alone?
    inter = g("err_mae") * (1.0 / np.clip(g("post_sd"), 1e-9, None))
    s_i = float(spearmanr(inter, rp).statistic)
    out["F3d"] = dict(error_alone=tab["err_mae"]["spearman"], error_times_confidence=s_i,
                      gain=float(abs(s_i) - abs(tab["err_mae"]["spearman"])), bar=0.10,
                      fires=bool(abs(s_i) - abs(tab["err_mae"]["spearman"]) >= 0.10),
                      note="uses F3-family posterior material; DIAGNOSTIC ONLY, cannot become a detector")

    # -------- does in-pool skill explain the S30-L3 mechanism quantity?
    smb = g("set_mean_benefit")
    out["link_to_S30_L3"] = dict(
        spearman_rho_pool_vs_set_mean_benefit=float(spearmanr(rp, smb).statistic),
        pearson=float(np.corrcoef(rp, smb)[0, 1]))

    out["per_target"] = [{k: r[k] for k in r} for r in R]
    out["provenance"] = ST.provenance(__file__)
    ST.save_atomic(OUT, out)

    print("=== F3a  IN-POOL RANKING SKILL OF THE SHIPPED SCORE (ORACLE), by stratum ===")
    print("%-20s %5s %10s %10s" % ("stratum", "n", "rho_pool", "rho_top75"))
    for nm in ("all", "other108", "FAIL18", "worst18_poolmean", "worst18_bestpool"):
        v = out["F3a"][nm]
        print("%-20s %5d %+10.4f %+10.4f" % (nm, v["n"], v["rho_pool"], v["rho_top75"]))
    st = out["F3a"]["stats"]
    print("  108 fold CI [%+.4f,%+.4f] %d/5 folds | FAIL18 fold CI [%+.4f,%+.4f] %d/%d folds"
          % (st["other108"]["ci95_fold"][0], st["other108"]["ci95_fold"][1], st["other108"]["folds_same_sign"],
             st["FAIL18"]["ci95_fold"][0], st["FAIL18"]["ci95_fold"][1], st["FAIL18"]["folds_same_sign"], 4))
    print("  random-18 null: obs %+.4f  null mean %+.4f  CI [%+.4f,%+.4f]  p=%.4g"
          % (obs, st["random18_null"]["null_mean"], st["random18_null"]["ci95"][0],
             st["random18_null"]["ci95"][1], st["random18_null"]["p_two_sided"]))
    print("  F3a fires:", out["F3a"]["fires"])
    print("\n=== WHAT PREDICTS rho_pool (Spearman over 126; ORACLE predictors marked *) ===")
    print("%-18s %9s %9s %-20s %9s %9s" % ("predictor", "spearman", "resid(n)", "fold CI", "F18 mean", "108 mean"))
    lab = {"err_mae": "*", "abs_scale_add": "*", "shape_err": "*", "abs_log_scale": "*"}
    for k in preds:
        t = tab[k]
        print("%-18s %+9.3f %+9.3f [%+.3f,%+.3f] %9.3f %9.3f %s"
              % (k, t["spearman"], t["spearman_resid_n"], t["ci95_fold"][0], t["ci95_fold"][1],
                 t["mean_F18"], t["mean_108"], lab.get(k, "")))
    print("\n  F3b (error) fires:", out["F3b"]["fires"], " stat", round(out["F3b"]["stat"], 3))
    print("  F3c partial scale|shape %+.3f (p=%.2g) ; shape|scale %+.3f (p=%.2g) ; fires %s"
          % (out["F3c"]["partial_scale_given_shape"], out["F3c"]["p"],
             out["F3c"]["partial_shape_given_scale"], out["F3c"]["p2"], out["F3c"]["fires"]))
    print("  F3d error*confidence %+.3f vs error alone %+.3f ; gain %+.3f ; fires %s"
          % (out["F3d"]["error_times_confidence"], out["F3d"]["error_alone"],
             out["F3d"]["gain"], out["F3d"]["fires"]))
    print("\n  link to S30-L3: spearman(rho_pool, filter set-mean benefit) = %+.3f"
          % out["link_to_S30_L3"]["spearman_rho_pool_vs_set_mean_benefit"])
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
