"""forensics: "what would have let us recognise the failure before it happened?"

`forensics_causal` asked this as CLASSIFICATION of FAIL18 membership and got AUC 0.600 inside a
label-shuffled null.  This asks the better-posed question: leave-fold-out REGRESSION of the
production answer's CA-RMSD on deployable, inference-time-computable features, plus the single-
feature LFO correlations.  A confidence/abstention signal is useful even if it cannot name the 18.
Null control: the same pipeline with the target labels permuted.
"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import spearmanr
from s12 import instrument as I
from s12 import forensics_lib as L

def main():
    d = json.load(open(os.path.join(I.RESULTS, "forensics_part1.json")))["rows"]
    alpha = {x["pdb"]: x for x in json.load(open(os.path.join(I.RESULTS, "forensics_alphabet.json")))["rows"]}
    rgz = np.load(os.path.join(I.CACHE, "forensics_rgpred.npz"), allow_pickle=True)
    rgi = {str(p): k for k, p in enumerate(rgz["pdb"])}
    names = [r["pdb"] for r in d]; folds = np.array([r["fold"] for r in d])
    F = np.array([r["fail18"] for r in d]); y = np.array([r["rmsd_arm"] for r in d])
    feats = {
        "top75_pair_rmsd": [r["top75"]["pair_rmsd_mean"] for r in d],
        "top75_rg_sd": [r["top75"]["rg_sd"] for r in d],
        "top75_nclusters": [r["top75"]["nclusters"] for r in d],
        "pool_nclusters": [r["pool_nclusters"] for r in d],
        "disto_rg": [r["disto_rg"] for r in d],
        "disto_E": [np.mean([c == "E" for c in r["disto_ss"]]) for r in d],
        "disto_H": [np.mean([c == "H" for c in r["disto_ss"]]) for r in d],
        "esm_ss_E": [np.mean([c == "E" for c in alpha[p]["pred_ss"]]) for p in names],
        "esm_ss_H": [np.mean([c == "H" for c in alpha[p]["pred_ss"]]) for p in names],
        "alpha_ll": [alpha[p]["mean_ll_univ"] for p in names],
        "n": [r["n"] for r in d],
        "nw": [np.log10(r["nw"]) for r in d],
        "org_frac_top75": [r["top75"]["org_frac"] for r in d],
        "rg_pred_gbt": [float(rgz["gbt"][rgi[p]]) for p in names],
        "rg_gap_disto_pool": [r["disto_rg"] - r["top75"]["rg_mean"] for r in d],
        "top75_rg_mean": [r["top75"]["rg_mean"] for r in d],
        "sim_pool_min": [r["band_univ"].get("sim_pool_min", np.nan) for r in d],
    }
    X = np.array([feats[k] for k in feats], float).T
    X = np.where(np.isfinite(X), X, np.nanmean(np.where(np.isfinite(X), X, np.nan), 0))
    out = {"features": list(feats), "n": len(y)}
    # single-feature LFO spearman (fit is rank-only, so LFO == in-sample for a monotone feature;
    # report the in-sample rank correlation and a permutation p)
    rng = np.random.default_rng(0)
    single = {}
    for k in feats:
        v = np.array(feats[k], float); m = np.isfinite(v)
        rho = spearmanr(v[m], y[m]).correlation
        null = np.array([abs(spearmanr(rng.permutation(v[m]), y[m]).correlation) for _ in range(2000)])
        single[k] = dict(rho=float(rho), p_perm=float((null >= abs(rho)).mean()),
                         rho_fail18=float(spearmanr(v[m & F], y[m & F]).correlation) if (m & F).sum() > 3 else None,
                         rho_other=float(spearmanr(v[m & ~F], y[m & ~F]).correlation))
    out["single_feature"] = single
    # LFO regression
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    def lfo(yy, seed=0):
        pr, pg = np.zeros(len(yy)), np.zeros(len(yy))
        for f in range(5):
            tr, te = folds != f, folds == f
            sc = StandardScaler().fit(X[tr]); A = sc.transform(X[tr]); B = sc.transform(X[te])
            pr[te] = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(A, yy[tr]).predict(B)
            pg[te] = GradientBoostingRegressor(n_estimators=250, max_depth=2, learning_rate=0.05, subsample=0.8, random_state=seed).fit(A, yy[tr]).predict(B)
        return pr, pg
    pr, pg = lfo(y)
    for nm, p in (("ridge", pr), ("gbt", pg)):
        o = np.argsort(-p)
        out[nm] = dict(spearman=float(spearmanr(p, y).correlation), mae=float(np.abs(p - y).mean()),
                       rmsd_of_worst18_flagged=float(y[o[:18]].mean()), n_fail18_in_worst18=int(F[o[:18]].sum()),
                       rmsd_of_rest=float(y[o[18:]].mean()),
                       # abstention curve: mean RMSD after dropping the k flagged-worst targets
                       abstain=[float(y[o[k:]].mean()) for k in (0, 6, 12, 18, 25)])
    # null: permuted labels
    nulls = []
    for s in range(30):
        yp = np.random.default_rng(s).permutation(y)
        _, pgn = lfo(yp, seed=s)
        nulls.append(spearmanr(pgn, yp).correlation)
    out["null_spearman_mean"] = float(np.mean(nulls)); out["null_spearman_p95"] = float(np.percentile(nulls, 95))
    # what does perfect abstention buy? (ORACLE reference)
    o = np.argsort(-y)
    out["oracle_abstain"] = [float(y[o[k:]].mean()) for k in (0, 6, 12, 18, 25)]
    I.write("forensics_warning", out)
    print("single-feature Spearman with the answer RMSD (deployable features):")
    for k, v in sorted(single.items(), key=lambda kv: -abs(kv[1]["rho"])):
        print(f"  {k:20s} rho={v['rho']:+.3f} p_perm={v['p_perm']:.4f}  F18 {str(v['rho_fail18'])[:6]:>6}  O108 {v['rho_other']:+.3f}")
    print(json.dumps({k: out[k] for k in ("ridge", "gbt", "null_spearman_mean", "null_spearman_p95", "oracle_abstain")}, indent=1))

if __name__ == "__main__":
    main()
