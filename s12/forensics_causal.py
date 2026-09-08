"""forensics Part 1b: find the COMMON CAUSAL STRUCTURE behind FAIL18.

Tests, on all 126 targets, which target-level properties predict failure, and whether the
candidate causes (native strand content, lasso/fibril metadata, library SS coverage,
pool multimodality) are the same fact or separate ones.  Also produces the per-target
FAIL18 detail table and a deployable-only early-warning check.
"""
import os, sys, json, itertools, collections
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact, pearsonr
from s12 import instrument as I
from s12 import forensics_lib as L

def get(r, path, default=np.nan):
    x = r
    for k in path.split("."):
        try:
            x = x[k]
        except Exception:
            return default
    return float(x) if not isinstance(x, bool) else float(x)

def main():
    d = json.load(open(os.path.join(I.RESULTS, "forensics_part1.json")))
    rows = d["rows"]; names = [r["pdb"] for r in rows]
    F = np.array([r["fail18"] for r in rows])
    alpha = {x["pdb"]: x for x in json.load(open(os.path.join(I.RESULTS, "forensics_alphabet.json")))["rows"]}
    out = {}

    # ---------------------------------------------------------------- 1. candidate causes
    feats = {
        "o_native_E": [get(r, "native_ss_comp.E") for r in rows],
        "o_native_H": [get(r, "native_ss_comp.H") for r in rows],
        "o_native_C": [get(r, "native_ss_comp.C") for r in rows],
        "o_pool_best": [get(r, "o_pool_best") for r in rows],
        "o_univ_best": [get(r, "o_univ_best") for r in rows],
        "o_cov20": [np.log10(1 + get(r, "o_cov20")) for r in rows],
        "o_lib_ss_agree_pool": [get(r, "bg.ss_agree_pool") for r in rows],
        "o_band_score_pct": [get(r, "band_pool.score_pct_mean") for r in rows],
        "rho_pool": [get(r, "rho_pool") for r in rows],
        "meta_lasso_or_fibril": [float(r["header"]["lasso"] or r["header"]["fibril"]) for r in rows],
        "meta_any_flag": [float(r["header"]["any_flag"]) for r in rows],
        "dep_disto_rg": [get(r, "disto_rg") for r in rows],
        "dep_disto_E": [float(L.SSINV.index("E") == 1) * np.mean([c == "E" for c in r["disto_ss"]]) for r in rows],
        "dep_disto_H": [np.mean([c == "H" for c in r["disto_ss"]]) for r in rows],
        "dep_pool_nclusters": [get(r, "pool_nclusters") for r in rows],
        "dep_top75_pair_rmsd": [get(r, "top75.pair_rmsd_mean") for r in rows],
        "dep_top75_rg_sd": [get(r, "top75.rg_sd") for r in rows],
        "dep_n": [get(r, "n") for r in rows],
        "dep_alpha_ll_band": [alpha[p].get("mean_ll_univ", np.nan) for p in names],
        "dep_ss_esm_E": [np.mean([c == "E" for c in alpha[p]["pred_ss"]]) for p in names],
        "dep_ss_esm_H": [np.mean([c == "H" for c in alpha[p]["pred_ss"]]) for p in names],
    }
    y = np.array([get(r, "rmsd_arm") for r in rows])
    tab = {}
    for k, v in feats.items():
        v = np.array(v, float)
        m = np.isfinite(v)
        tab[k] = dict(fail18_mean=float(np.nanmean(v[F])), other_mean=float(np.nanmean(v[~F])),
                      mwu_p=float(mannwhitneyu(v[F & m], v[~F & m]).pvalue),
                      spearman_rmsd=float(spearmanr(v[m], y[m]).correlation),
                      auc_fail=float((mannwhitneyu(v[F & m], v[~F & m]).statistic) / (F.sum() * (~F).sum())))
    out["causes"] = tab

    # ---------------------------------------------------------------- 2. is E-content the same fact as lasso/fibril?
    E = np.array(feats["o_native_E"]); LF = np.array(feats["meta_lasso_or_fibril"]).astype(bool)
    out["strand_vs_meta"] = dict(
        E_mean_lassofib=float(E[LF].mean()), E_mean_other=float(E[~LF].mean()),
        fail_rate_LF=f"{int(F[LF].sum())}/{int(LF.sum())}", fail_rate_notLF=f"{int(F[~LF].sum())}/{int((~LF).sum())}",
        fail_rate_highE_notLF=f"{int(F[(E >= 0.35) & ~LF].sum())}/{int(((E >= 0.35) & ~LF).sum())}",
        fail_rate_lowE_notLF=f"{int(F[(E < 0.35) & ~LF].sum())}/{int(((E < 0.35) & ~LF).sum())}",
        fail_rate_highE_LF=f"{int(F[(E >= 0.35) & LF].sum())}/{int(((E >= 0.35) & LF).sum())}",
        fisher_highE_within_notLF=float(fisher_exact([[int(F[(E >= 0.35) & ~LF].sum()), int(((E >= 0.35) & ~LF).sum() - F[(E >= 0.35) & ~LF].sum())],
                                                      [int(F[(E < 0.35) & ~LF].sum()), int(((E < 0.35) & ~LF).sum() - F[(E < 0.35) & ~LF].sum())]])[1]),
        rmsd_highE=float(y[E >= 0.35].mean()), rmsd_lowE=float(y[E < 0.35].mean()),
        rmsd_LF=float(y[LF].mean()), rmsd_notLF=float(y[~LF].mean()))

    # ---------------------------------------------------------------- 3. the filter's SS bias
    # what SS does the score prefer, conditioned on the native's SS?
    hi = E >= 0.35
    out["ss_bias"] = dict(
        top75_E_when_native_E=float(np.mean([get(r, "top75.ss_comp.E") for r, m in zip(rows, hi) if m])),
        top75_E_when_native_notE=float(np.mean([get(r, "top75.ss_comp.E") for r, m in zip(rows, hi) if not m])),
        top75_H_when_native_E=float(np.mean([get(r, "top75.ss_comp.H") for r, m in zip(rows, hi) if m])),
        pool_E_frac_when_native_E=None, n_highE=int(hi.sum()))
    # pool-level E fraction requires the window SS: compute
    pe_hi, pe_lo, ue_hi, ue_lo = [], [], [], []
    for r, m in zip(rows, hi):
        u = I.load_univ(r["pdb"]); cass = L.window_cass(r["pdb"], u); pool = I.pool_idx(u)
        pf = float((cass[pool] == 1).mean()); uf = float((cass == 1).mean())
        (pe_hi if m else pe_lo).append(pf); (ue_hi if m else ue_lo).append(uf)
    out["ss_bias"]["pool_E_frac_when_native_E"] = float(np.mean(pe_hi))
    out["ss_bias"]["pool_E_frac_when_native_notE"] = float(np.mean(pe_lo))
    out["ss_bias"]["univ_E_frac_when_native_E"] = float(np.mean(ue_hi))
    out["ss_bias"]["univ_E_frac_when_native_notE"] = float(np.mean(ue_lo))

    # ---------------------------------------------------------------- 4. deployable early warning
    # can we PREDICT failure from deployable quantities only?  LFO logistic regression.
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    dep = [k for k in feats if k.startswith("dep_")]
    X = np.array([feats[k] for k in dep], float).T
    X = np.nan_to_num(X, nan=np.nanmean(X))
    folds = np.array([r["fold"] for r in rows])
    p_hat = np.zeros(len(rows))
    for f in range(5):
        tr, te = folds != f, folds == f
        sc = StandardScaler().fit(X[tr])
        m = LogisticRegression(max_iter=2000, C=0.5).fit(sc.transform(X[tr]), F[tr])
        p_hat[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
    from sklearn.metrics import roc_auc_score
    out["early_warning"] = dict(features=dep, auc_lfo=float(roc_auc_score(F, p_hat)),
                                spearman_with_rmsd=float(spearmanr(p_hat, y).correlation),
                                top18_precision=float(F[np.argsort(-p_hat)[:18]].mean()),
                                rmsd_top18_flagged=float(y[np.argsort(-p_hat)[:18]].mean()),
                                rmsd_rest=float(y[np.argsort(-p_hat)[18:]].mean()))
    # null control: shuffle labels
    rng = np.random.default_rng(0); aucs = []
    for _ in range(200):
        Fp = rng.permutation(F); ph = np.zeros(len(rows))
        for f in range(5):
            tr, te = folds != f, folds == f
            if len(set(Fp[tr])) < 2:
                ph[te] = 0.5; continue
            sc = StandardScaler().fit(X[tr]); m = LogisticRegression(max_iter=500, C=0.5).fit(sc.transform(X[tr]), Fp[tr])
            ph[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
        aucs.append(roc_auc_score(Fp, ph))
    out["early_warning"]["null_auc_mean"] = float(np.mean(aucs)); out["early_warning"]["null_auc_p95"] = float(np.percentile(aucs, 95))

    # ---------------------------------------------------------------- 5. per-target FAIL18 detail
    detail = []
    for r in rows:
        if not r["fail18"]:
            continue
        a = alpha[r["pdb"]]
        detail.append(dict(pdb=r["pdb"], n=r["n"], seq=r["seq"], o_native_ss=r["o_native_ss"], disto_ss=r["disto_ss"],
                           pred_ss_esm=a["pred_ss"], o_pool_best=r["o_pool_best"], o_univ_best=r["o_univ_best"],
                           shipped=r["shipped"], rmsd_arm=r["rmsd_arm"], o_cov20=r["o_cov20"], o_cov20_pool=r["o_cov20_pool"],
                           n_band_pool=r["n_band_pool"], band_score_pct=get(r, "band_pool.score_pct_mean"),
                           band_score_pct_min=get(r, "band_pool.score_pct_min"),
                           band_org_frac=get(r, "band_univ.org_frac"), band_rank_med=get(r, "band_univ.rank_med"),
                           band_in_pool_frac=get(r, "band_univ.in_pool_frac"), band_ident=get(r, "band_univ.ident_mean"),
                           band_ss_agree=get(r, "band_univ.ss_agree_mean"), lib_ss_agree=get(r, "bg.ss_agree_pool"),
                           top75_ss_agree=get(r, "top75.ss_agree"), top75_ss_mode=r["top75"]["ss_mode"],
                           top75_o_rr_min=get(r, "top75.o_rr_min"), top75_in_near_cluster=r["top75_in_near_cluster"],
                           pool_nclusters=r["pool_nclusters"], rho_pool=r["rho_pool"], rg_native=r["o_native_rg"],
                           rg_top75=get(r, "top75.rg_mean"), disto_rg=r["disto_rg"], diag=r["diag"],
                           band_parents=r["band_univ"].get("parents", [])[:3], header={k: r["header"][k] for k in ("lasso", "fibril", "xray", "membrane", "cosolvent", "bound", "designed", "ssbond", "expdta")}))
    out["fail18_detail"] = detail
    I.write("forensics_causal", out)

    print("=== candidate causes (auc_fail = P(FAIL ranks higher)) ===")
    for k, v in sorted(out["causes"].items(), key=lambda kv: -abs(kv[1]["auc_fail"] - 0.5)):
        print(f"{k:24s} F18 {v['fail18_mean']:8.3f}  O108 {v['other_mean']:8.3f}  p={v['mwu_p']:9.2e}  auc={v['auc_fail']:.3f}  rho_rmsd={v['spearman_rmsd']:+.3f}")
    print("\n=== strand vs metadata ===");  print(json.dumps(out["strand_vs_meta"], indent=1))
    print("\n=== SS bias ===");  print(json.dumps(out["ss_bias"], indent=1))
    print("\n=== deployable early warning ===");  print(json.dumps(out["early_warning"], indent=1))

if __name__ == "__main__":
    main()
