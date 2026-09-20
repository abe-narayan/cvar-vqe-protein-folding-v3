#!/usr/bin/env python
"""s30/s30_F_stagegap.py -- S30 lane F, experiment F1.

Pre-registration: `s30/PREREG_S30_F1.md` (committed 73626721, before this file existed).

Is the FAIL18 tail POOL-LIMITED (nothing good in the candidate set) or SELECTION-LIMITED (good
members present and discarded)?  Decomposes the production -> ORACLE gap across the three
narrowing stages of the shipped path and controls each against a random subset drawn in the
SAME operator space (contract rule 6), because best-of-k is an order statistic.

EVERY RMSD IN THIS FILE IS AN ORACLE QUANTITY except `prod` (the deployable emitted structure).
`rr` is a POINT CLOUD RMSD; the built-chain counterparts are quoted from S29 lane O's ladder
table and are never mixed with these in one row.

    python s30/s30_F_stagegap.py            # writes s30/results/s30_F_stagegap.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.special import gammaln

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
OUT = os.path.join(RESULTS, "s30_F_stagegap.json")

K, M = 500, 75
N_DRAW = 2000          # prereg: 2,000 draws per target, per control
SEED_CTRL = 30001      # prereg
SEED_NULL = 30002      # prereg
N_NULL = 20000         # prereg: 20,000-draw random-18 null
CAP = 3.00             # prereg F1a threshold


def rows():
    out = []
    for t in I.targets():
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        rec = I.shipped_record(pdb)
        rr = np.asarray(u["rr"], float)
        order = np.asarray(u["order"], int)
        pool = order[:K]
        sub = np.asarray(rec["sub"], int)           # indices INTO THE POOL (verified)
        assert sub.min() >= 0 and sub.max() < K, pdb
        rr_pool = rr[pool]
        rr_top = rr_pool[sub]
        out.append(dict(
            pdb=pdb, n=int(t["n"]), fold=int(t["fold"]),
            fail18=bool(pdb in I.FAIL18),
            n_win=int(len(rr)),
            best_univ=float(rr.min()),
            best_pool=float(rr_pool.min()),
            best_top=float(rr_top.min()),
            mean_pool=float(rr_pool.mean()),
            mean_top=float(rr_top.mean()),
            mean_univ=float(rr.mean()),
            prod=float(rec["rmsd_avg"]),             # emitted, point cloud
            rr_pool=rr_pool, rr_univ=rr,             # kept in memory only
        ))
    return out


def _log_nCk(n, k):
    """log C(n,k), vectorised over n; -inf where k > n."""
    n = np.asarray(n, float)
    out = np.full(n.shape, -np.inf)
    ok = n >= k
    out[ok] = (gammaln(n[ok] + 1.0) - gammaln(k + 1.0) - gammaln(n[ok] - k + 1.0))
    return out


def ctrl_best_of_subset_exact(values, k, obs):
    """EXACT distribution of min over a subset of size k drawn WITHOUT replacement -- the
    filter's own operator space (it retains k DISTINCT members).

    DEVIATION FROM PREREG, declared: the prereg registered a 2,000-draw sampler.  The min of a
    random k-subset depends only on the RANK of the smallest index drawn, whose law is exact and
    closed-form:  P(min rank >= i) = C(n-i, k) / C(n, k).  This computes the same estimand the
    registered sampler approximates, with zero Monte-Carlo error, and cannot be tuned (it has no
    free parameter).  A 2,000-draw sampler is run alongside on the k=75 arm as a cross-check.

    Returns (expected_best, percentile_of_obs) where percentile is P(random best < obs) +
    0.5 * P(random best == obs); 0.5 = no skill, > 0.5 = worse than random.
    """
    v = np.sort(np.asarray(values, float))
    n = len(v)
    if k >= n:
        return float(v[0]), 0.5
    i = np.arange(0, n - k + 1, dtype=float)              # possible ranks of the minimum
    log_tot = _log_nCk(np.array([float(n)]), k)[0]
    surv = np.exp(_log_nCk(n - i, k) - log_tot)           # P(min rank >= i), surv[0] == 1
    pmf = np.empty_like(surv)
    pmf[:-1] = surv[:-1] - surv[1:]
    pmf[-1] = surv[-1]
    exp_best = float((v[:len(pmf)] * pmf).sum())
    lt = int(np.searchsorted(v, obs, side="left"))        # number of values strictly < obs
    eq = int(np.searchsorted(v, obs, side="right")) - lt
    p_ge_obs = float(np.exp(_log_nCk(np.array([float(n - lt)]), k)[0] - log_tot))
    p_gt_obs = float(np.exp(_log_nCk(np.array([float(n - lt - eq)]), k)[0] - log_tot))
    p_lt = 1.0 - p_ge_obs
    p_eq = p_ge_obs - p_gt_obs
    return exp_best, float(p_lt + 0.5 * p_eq)


def ctrl_best_of_subset_sampled(values, k, n_draw, rng):
    """The prereg's registered 2,000-draw sampler, kept as a cross-check on the k=75 arm."""
    v = np.asarray(values, float)
    n = len(v)
    if k >= n:
        return float(v.min())
    best = np.empty(n_draw, float)
    for i in range(n_draw):
        best[i] = v[rng.choice(n, k, replace=False)].min()
    return float(best.mean())


def main():
    R = rows()
    pdbs = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    isf = np.array([r["fail18"] for r in R])
    assert isf.sum() == 18, isf.sum()

    # ---------------------------------------------------------------- controls
    rng = np.random.default_rng(SEED_CTRL)
    for r in R:
        e75, p75 = ctrl_best_of_subset_exact(r["rr_pool"], M, r["best_top"])
        r["ctrl_best_r75"], r["pct_filt"] = e75, p75
        r["ctrl_best_r75_sampled"] = ctrl_best_of_subset_sampled(r["rr_pool"], M, N_DRAW, rng)
        e500, p500 = ctrl_best_of_subset_exact(r["rr_univ"], K, r["best_pool"])
        r["ctrl_best_r500"], r["pct_retr"] = e500, p500

    A = {k: np.array([r[k] for r in R], float) for k in
         ("best_univ", "best_pool", "best_top", "prod", "mean_pool", "mean_top", "mean_univ",
          "ctrl_best_r75", "ctrl_best_r75_sampled", "ctrl_best_r500", "pct_filt", "pct_retr")}

    G_retr = A["best_pool"] - A["best_univ"]
    G_filt = A["best_top"] - A["best_pool"]
    G_read = A["prod"] - A["best_top"]
    G_tot = A["prod"] - A["best_pool"]

    def strat(x):
        return dict(all=float(x.mean()), fail18=float(x[isf].mean()), other108=float(x[~isf].mean()),
                    median_all=float(np.median(x)), median_fail18=float(np.median(x[isf])),
                    median_other108=float(np.median(x[~isf])))

    out = dict(
        prereg="s30/PREREG_S30_F1.md",
        basis="POINT CLOUD (rr); built-chain counterparts quoted from s29 lane O, never mixed",
        n=len(R), n_fail18=int(isf.sum()),
        stages={k: strat(A[k]) for k in
                ("best_univ", "best_pool", "best_top", "prod", "mean_pool", "mean_top",
                 "ctrl_best_r75", "ctrl_best_r75_sampled", "ctrl_best_r500")},
        gaps=dict(G_retr=strat(G_retr), G_filt=strat(G_filt), G_read=strat(G_read),
                  G_tot=strat(G_tot)),
    )

    # ------------------------------------------------- F1a: is the pool good enough on the tail
    out["F1a"] = dict(
        claim="ORACLE best-of-pool on FAIL18 is below the 3.00 A cap",
        cap=CAP,
        oracle_best_pool_fail18=float(A["best_pool"][isf].mean()),
        oracle_best_pool_other108=float(A["best_pool"][~isf].mean()),
        n_fail18_under_cap=int((A["best_pool"][isf] < CAP).sum()),
        n_other108_under_cap=int((A["best_pool"][~isf] < CAP).sum()),
        worst_fail18_best_pool=float(A["best_pool"][isf].max()),
        fires=bool(A["best_pool"][isf].mean() < CAP),
    )

    # ------------------------------------------------- F1b: is the filter cost distinctively tail
    sf_f = float(G_filt[isf].sum() / G_tot[isf].sum())
    sf_o = float(G_filt[~isf].sum() / G_tot[~isf].sum())
    rng2 = np.random.default_rng(SEED_NULL)
    idx = np.arange(len(R))
    null_gfilt = np.empty(N_NULL); null_share = np.empty(N_NULL)
    for i in range(N_NULL):
        s = rng2.choice(idx, 18, replace=False)
        null_gfilt[i] = G_filt[s].mean()
        null_share[i] = G_filt[s].sum() / G_tot[s].sum()
    obs_g = float(G_filt[isf].mean())
    out["F1b"] = dict(
        share_filt_fail18=sf_f, share_filt_other108=sf_o,
        share_ratio=float(sf_f / sf_o) if sf_o else None,
        clause1_share_ratio_ge_2=bool(sf_f >= 2.0 * sf_o),
        random18_null=dict(
            n_draw=N_NULL, seed=SEED_NULL,
            obs_G_filt_fail18=obs_g,
            null_mean=float(null_gfilt.mean()),
            null_ci95=[float(np.percentile(null_gfilt, 2.5)), float(np.percentile(null_gfilt, 97.5))],
            p_two_sided=float(2 * min((null_gfilt >= obs_g).mean(), (null_gfilt <= obs_g).mean())),
            obs_share=sf_f,
            null_share_ci95=[float(np.percentile(null_share, 2.5)),
                             float(np.percentile(null_share, 97.5))],
            p_share_two_sided=float(2 * min((null_share >= sf_f).mean(), (null_share <= sf_f).mean())),
        ),
    )
    out["F1b"]["clause2_outside_null"] = bool(
        obs_g > out["F1b"]["random18_null"]["null_ci95"][1] or
        obs_g < out["F1b"]["random18_null"]["null_ci95"][0])

    # difficulty control: regress G_filt on production RMSD over the 108, extrapolate to FAIL18
    x = A["prod"][~isf]; y = G_filt[~isf]
    b1, b0 = np.polyfit(x, y, 1)
    resid = y - (b0 + b1 * x)
    s_res = float(resid.std(ddof=2))
    xf = A["prod"][isf]; yf = G_filt[isf]
    pred = b0 + b1 * xf
    n108 = len(x); xbar = x.mean(); sxx = float(((x - xbar) ** 2).sum())
    se_pred = s_res * np.sqrt(1.0 + 1.0 / n108 + (xf - xbar) ** 2 / sxx)
    z = (yf - pred) / se_pred
    out["F1b"]["difficulty_control"] = dict(
        slope=float(b1), intercept=float(b0), resid_sd=s_res,
        fail18_mean_obs=float(yf.mean()), fail18_mean_pred=float(pred.mean()),
        fail18_mean_residual=float((yf - pred).mean()),
        n_fail18_above_95_band=int((z > 1.96).sum()),
        n_fail18_below_95_band=int((z < -1.96).sum()),
        mean_z=float(z.mean()),
    )
    out["F1b"]["clause3_survives_difficulty"] = bool((yf - pred).mean() > 0 and (z > 1.96).sum() >= 1)
    out["F1b"]["fires"] = bool(out["F1b"]["clause1_share_ratio_ge_2"] and
                               out["F1b"]["clause2_outside_null"] and
                               out["F1b"]["clause3_survives_difficulty"])

    # ------------------------------------------------- F1c: does the filter have negative skill
    cmp_all = ST.compare(A["best_top"], A["ctrl_best_r75"], folds=folds, names=pdbs,
                         label="F1c_filter_vs_random75_all", seed_parts=("s30F",))
    cmp_f = ST.compare(A["best_top"][isf], A["ctrl_best_r75"][isf], folds=folds[isf],
                       names=[p for p, f in zip(pdbs, isf) if f],
                       label="F1c_filter_vs_random75_FAIL18", seed_parts=("s30F",))
    cmp_o = ST.compare(A["best_top"][~isf], A["ctrl_best_r75"][~isf], folds=folds[~isf],
                       names=[p for p, f in zip(pdbs, isf) if not f],
                       label="F1c_filter_vs_random75_other108", seed_parts=("s30F",))
    cmp_r = ST.compare(A["best_pool"], A["ctrl_best_r500"], folds=folds, names=pdbs,
                       label="F1c_retrieval_vs_random500_all", seed_parts=("s30F",))
    cmp_rf = ST.compare(A["best_pool"][isf], A["ctrl_best_r500"][isf], folds=folds[isf],
                        names=[p for p, f in zip(pdbs, isf) if f],
                        label="F1c_retrieval_vs_random500_FAIL18", seed_parts=("s30F",))
    cmp_ro = ST.compare(A["best_pool"][~isf], A["ctrl_best_r500"][~isf], folds=folds[~isf],
                        names=[p for p, f in zip(pdbs, isf) if not f],
                        label="F1c_retrieval_vs_random500_other108", seed_parts=("s30F",))
    out["F1c"] = dict(
        filter_all=cmp_all, filter_fail18=cmp_f, filter_other108=cmp_o,
        retrieval_all=cmp_r, retrieval_fail18=cmp_rf, retrieval_other108=cmp_ro,
        pct_filt=dict(median_all=float(np.median(A["pct_filt"])),
                      median_fail18=float(np.median(A["pct_filt"][isf])),
                      median_other108=float(np.median(A["pct_filt"][~isf]))),
        pct_retr=dict(median_all=float(np.median(A["pct_retr"])),
                      median_fail18=float(np.median(A["pct_retr"][isf])),
                      median_other108=float(np.median(A["pct_retr"][~isf]))),
    )
    out["F1c"]["fires"] = bool(cmp_f["effect"] > 0 and
                               np.median(A["pct_filt"][isf]) > 0.5 and
                               cmp_f["ci95_fold"] is not None and
                               cmp_f["ci95_fold"][0] * cmp_f["ci95_fold"][1] > 0)

    out["per_target"] = [
        {k: r[k] for k in ("pdb", "n", "fold", "fail18", "n_win", "best_univ", "best_pool",
                           "best_top", "mean_pool", "mean_top", "prod",
                           "ctrl_best_r75", "ctrl_best_r75_sampled", "ctrl_best_r500",
                           "pct_filt", "pct_retr")}
        for r in R]
    out["provenance"] = ST.provenance(__file__)
    ST.save_atomic(OUT, out)
    print(json.dumps({k: out[k] for k in ("stages", "gaps", "F1a")}, indent=1)[:4000])
    print("\nF1b", json.dumps(out["F1b"], indent=1)[:2500])
    print("\nF1c fires:", out["F1c"]["fires"])
    for nm in ("filter_all", "filter_fail18", "filter_other108",
               "retrieval_all", "retrieval_fail18", "retrieval_other108"):
        c = out["F1c"][nm]
        print("  %-24s eff %+.4f  xMDE %+.2f  foldCI [%+.4f,%+.4f]  %d/%d folds  %dW/%dL  %s" % (
            nm, c["effect"], c["effect_over_mde"], c["ci95_fold"][0], c["ci95_fold"][1],
            c["folds_same_sign"], c["n_folds"], c["n_better"], c["n_worse"], c["verdict"]))
    print("  pct_filt medians", out["F1c"]["pct_filt"])
    print("  pct_retr medians", out["F1c"]["pct_retr"])
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
