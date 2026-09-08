"""forensics Part 1c: the compactness hypothesis.

H: the distogram's shrinkage (calibration slope 0.376) makes its expected distances too
short for EXTENDED targets, so the Bayes-risk score prefers compact windows; the FAIL18
natives are extended, so their near-native band is scored last.  Controls: (a) is the
score-vs-rg relation the same on the 108? (b) does the distogram's own predicted rg track
the native rg, and with what slope? (c) an rg-matched ORACLE control: restrict the pool to
windows whose rg is within 0.5 A of the native rg and re-run the chain (DIAGNOSTIC, bounds
what perfect knowledge of a single scalar - the size - would buy).
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
from s12 import instrument as I
from s12 import forensics_lib as L

def per_target(t, do_project=True):
    import torch; torch.set_num_threads(2)
    pdb, n = t["pdb"], t["n"]
    u = I.load_univ(pdb); rr = u["rr"]; nat = u["nat_ca"]
    dg = I.distogram(pdb); i, j = I.pair_index(n)
    pool = I.pool_idx(u); Wp = u["W"][pool]
    D = I.pair_dists(Wp, i, j).astype(np.float32).astype(float)
    sc = I.shipped_score(dg, D)
    rgp = np.array([L.rg(w) for w in Wp]); rgn = L.rg(nat)
    e2e = float(np.linalg.norm(nat[-1] - nat[0])) / (n - 1)
    e2ep = np.linalg.norm(Wp[:, -1] - Wp[:, 0], axis=1) / (n - 1)
    # distogram-implied rg: sqrt(sum E[d_ij]^2 / n^2) using predicted for |i-j|>=2 and ideal else
    assert (dg["i"] == i).all() and (dg["j"] == j).all(), "pair index mismatch"
    ex = dg["expected"]; s2 = 2 * (ex ** 2).sum() + 2 * (n - 1) * 3.8 ** 2
    rg_pred = float(np.sqrt(s2 / (2 * n ** 2)))
    out = dict(pdb=pdb, n=n, fail18=pdb in I.FAIL18, o_rg_nat=rgn, rg_pred=rg_pred,
               o_e2e_nat=e2e, rg_pool_mean=float(rgp.mean()), e2e_pool_mean=float(e2ep.mean()),
               rho_score_rg=float(spearmanr(sc, rgp).correlation),
               rho_score_e2e=float(spearmanr(sc, e2ep).correlation),
               o_rho_rr_rg=float(spearmanr(rr[pool], rgp).correlation),
               o_rg_pct_in_pool=float((rgp < rgn).mean()),
               rg_top75=float(rgp[np.argsort(sc, kind="stable")[:I.M]].mean()),
               o_rg_band=float(rgp[rr[pool] <= rr[pool].min() + I.BAND].mean()))
    # oracle rg-matched pool (DIAGNOSTIC): among the K=500 BLOSUM pool, keep windows with |rg - rg_nat| < 0.5
    keep = np.where(np.abs(rgp - rgn) < 0.5)[0]
    out["n_rgmatch"] = int(len(keep))
    if len(keep) >= I.M:
        r = L.chain(u, pool[keep], dg, t["seq"], t["fold"], do_project=do_project)
        out["rgmatch"] = {k: v for k, v in r.items() if k != "top_idx_univ"}
    # oracle rg-matched retrieval over the WHOLE universe (DIAGNOSTIC): BLOSUM top-500 among rg-matched
    rgu = np.linalg.norm(u["W"] - u["W"].mean(1, keepdims=True), axis=2)
    rgu = np.sqrt((rgu ** 2).mean(1))
    keepu = np.where(np.abs(rgu - rgn) < 0.5)[0]
    out["n_rgmatch_univ"] = int(len(keepu))
    if len(keepu) >= I.K:
        sel = keepu[np.argsort(-u["sim"][keepu], kind="stable")[:I.K]]
        r = L.chain(u, sel, dg, t["seq"], t["fold"], do_project=do_project)
        out["rgmatch_univ"] = {k: v for k, v in r.items() if k != "top_idx_univ"}
    return out

if __name__ == "__main__":
    tg = I.targets(); rows = []; t0 = time.time()
    for k, t in enumerate(tg):
        rows.append(per_target(t))
        print(f"[{k+1:3d}/126] {t['pdb']} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
    F = np.array([r["fail18"] for r in rows]); names = [r["pdb"] for r in rows]
    folds = [t["fold"] for t in tg]
    agg = {}
    for m in ["o_rg_nat", "rg_pred", "o_e2e_nat", "rg_pool_mean", "rho_score_rg", "rho_score_e2e", "o_rho_rr_rg",
              "o_rg_pct_in_pool", "rg_top75", "o_rg_band", "n_rgmatch", "n_rgmatch_univ"]:
        v = np.array([r[m] for r in rows], float)
        agg[m] = dict(all=float(v.mean()), fail18=float(v[F].mean()), other=float(v[~F].mean()),
                      p=float(mannwhitneyu(v[F], v[~F]).pvalue))
    rgn = np.array([r["o_rg_nat"] for r in rows]); rgp = np.array([r["rg_pred"] for r in rows])
    sl = np.polyfit(rgn, rgp, 1)
    agg["rg_calibration"] = dict(slope=float(sl[0]), intercept=float(sl[1]), pearson=float(pearsonr(rgn, rgp)[0]),
                                 pred_minus_nat=float((rgp - rgn).mean()))
    base = json.load(open(os.path.join(I.RESULTS, "forensics_part1.json")))["rows"]
    b = {r["pdb"]: r for r in base}
    for arm in ("rgmatch", "rgmatch_univ"):
        got = [r for r in rows if arm in r]
        if not got:
            continue
        a = np.array([r[arm]["rmsd_arm"] for r in got]); bb = np.array([b[r["pdb"]]["rmsd_arm"] for r in got])
        f = np.array([r["fail18"] for r in got])
        agg[arm] = dict(n=len(got), paired=I.paired(a, bb, folds=[folds[names.index(r["pdb"])] for r in got], names=[r["pdb"] for r in got]),
                        fail18=I.paired(a[f], bb[f]) if f.sum() > 1 else None, other=I.paired(a[~f], bb[~f]),
                        pool_best=dict(arm=float(np.mean([r[arm]["pool_best"] for r in got])),
                                       base=float(np.mean([b[r["pdb"]]["o_pool_best"] for r in got]))))
        for q in ("fail18", "other"):
            if agg[arm][q]:
                agg[arm][q] = {x: agg[arm][q][x] for x in ("mean_a", "mean_b", "mean_diff", "ci95", "n_better", "n_worse")}
        agg[arm]["paired"].pop("top10_targets", None)
    I.write("forensics_compact", dict(rows=rows, aggregate=agg))
    print(json.dumps(agg, indent=1, default=str))
