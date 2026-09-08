"""forensics Part 1: per-target forensic table for tuning126, FAIL18 vs other-108.

Per target: header metadata class, native rg / CA-SS (labels), universe coverage
(1.0/1.5/2.0 A), provenance of the near-native band (org, BLOSUM rank, identity, ESM-key
percentile, SS agreement, hydrophobic-pattern agreement, in-pool score percentile),
pool / top-75 cluster structure, distogram-implied SS vs native SS, diagnosis flags.
Everything oracle-derived is labelled `o_`.
"""
import os, sys, time, json, collections
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
from s12 import instrument as I
from s12 import forensics_lib as L

def pct_rank(x, v):
    """percentile (0 best .. 1 worst) of value v among x (lower better)."""
    x = np.asarray(x, float); return float((x < v).mean())

def analyse(t):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    u = I.load_univ(pdb); rr = u["rr"]; nat = u["nat_ca"]
    pool = I.pool_idx(u); rec = I.shipped_record(pdb); sub = np.asarray(rec["sub"], int)
    dg = I.distogram(pdb); i, j = I.pair_index(n)
    Wp = u["W"][pool]; D = I.pair_dists(Wp, i, j).astype(np.float32).astype(float)
    sc = I.shipped_score(dg, D)
    pid, start = L.parent_map(pdb, u); P = L.parents()
    cass = L.window_cass(pdb, u); esm = L.esm_key_cached(pdb, u)
    rank = np.empty(len(u["order"]), int); rank[u["order"]] = np.arange(len(u["order"]))
    o_nss = L.ca_ss(nat); o_nrg = L.rg(nat)
    dss = L.distogram_ss(dg, n)
    hyd = L.hydro_pattern(seq)
    codes = u["S"]; tcode = np.array([I.ALPHABET.index(c) for c in seq])
    ident = (codes == tcode[None, :]).mean(1)
    esm_pct = np.empty(len(esm)); esm_pct[np.argsort(-esm, kind="stable")] = np.arange(len(esm)) / len(esm)
    pool_best = float(rr[pool].min()); thr = pool_best + I.BAND
    band_u = np.where(rr <= thr)[0]; band_p = np.where(rr[pool] <= thr)[0]
    ubest = np.argsort(rr, kind="stable")[:10]
    hyd_w = np.isin(codes, [I.ALPHABET.index(c) for c in L.HYDRO])
    out = dict(pdb=pdb, n=n, fold=fold, seq=seq, fail18=pdb in I.FAIL18, nw=int(len(rr)), org_frac=float(u["org"].mean()),
               header=L.header_class(pdb), o_native_ss=L.ss_str(o_nss), o_native_rg=o_nrg,
               disto_ss=L.ss_str(dss), disto_ss_agree=float((dss == o_nss).mean()),
               o_univ_best=float(rr.min()), o_pool_best=pool_best,
               o_cov10=int((rr <= 1.0).sum()), o_cov15=int((rr <= 1.5).sum()), o_cov20=int((rr <= 2.0).sum()),
               o_cov10_pool=int((rr[pool] <= 1.0).sum()), o_cov15_pool=int((rr[pool] <= 1.5).sum()), o_cov20_pool=int((rr[pool] <= 2.0).sum()),
               shipped=float(rr[pool[int(np.argmin(sc))]]), rmsd_arm=rec["rmsd_arm"], rmsd_avg=rec["rmsd_avg"],
               n_band_univ=int(len(band_u)), n_band_pool=int(len(band_p)))
    # universe's best-10 and band provenance
    def prov(idx):
        if len(idx) == 0:
            return {}
        return dict(n=int(len(idx)), org_frac=float(u["org"][idx].mean()), rank_med=float(np.median(rank[idx])),
                    rank_min=int(rank[idx].min()), in_pool_frac=float((rank[idx] < I.K).mean()),
                    sim_mean=float(u["sim"][idx].mean()), sim_pool_min=float(u["sim"][pool].min()),
                    ident_mean=float(ident[idx].mean()), esm_pct_med=float(np.median(esm_pct[idx])),
                    esm_key_mean=float(esm[idx].mean()),
                    ss_agree_mean=float((cass[idx] == o_nss[None, :]).mean()),
                    hyd_agree_mean=float((hyd_w[idx] == hyd[None, :]).mean()),
                    parents=[(str(P["pdb"][pid[k]]), int(start[k]), float(rr[k]), int(rank[k])) for k in idx[np.argsort(rr[idx])][:5]],
                    ss=[L.ss_str(cass[k]) for k in idx[np.argsort(rr[idx])][:3]])
    out["band_univ"] = prov(band_u); out["univ_top10"] = prov(ubest)
    # universe background for the same keys
    out["bg"] = dict(ident_mean=float(ident.mean()), ss_agree_mean=float((cass == o_nss[None, :]).mean()),
                     hyd_agree_mean=float((hyd_w == hyd[None, :]).mean()), esm_key_mean=float(esm.mean()),
                     ident_pool=float(ident[pool].mean()), ss_agree_pool=float((cass[pool] == o_nss[None, :]).mean()),
                     rg_pool_mean=float(np.mean([L.rg(w) for w in Wp])), rg_pool_sd=float(np.std([L.rg(w) for w in Wp])),
                     rg_univ_pct_native=float(np.mean(np.array([L.rg(w) for w in u["W"][::max(1, len(rr)//3000)]]) < o_nrg)))
    # in-pool band: score percentile
    if len(band_p):
        out["band_pool"] = dict(n=int(len(band_p)), score_pct_mean=float(np.mean([pct_rank(sc, sc[k]) for k in band_p])),
                                score_pct_min=float(min(pct_rank(sc, sc[k]) for k in band_p)),
                                in_top75=int(np.isin(band_p, sub).sum()),
                                ss_agree_mean=float((cass[pool[band_p]] == o_nss[None, :]).mean()),
                                org_frac=float(u["org"][pool[band_p]].mean()))
    rho = spearmanr(sc, rr[pool]).correlation
    out["rho_pool"] = float(rho)
    # pool clusters
    Ppair = I.pairwise_rmsd(Wp)
    ncl, lab = L.cluster_count(Ppair, 2.5)
    cl_stats = []
    for c in range(1, ncl + 1):
        m = lab == c
        if m.sum() < 5:
            continue
        cl_stats.append(dict(size=int(m.sum()), o_rr_min=float(rr[pool][m].min()), o_rr_mean=float(rr[pool][m].mean()),
                             score_mean=float(sc[m].mean()), n_top75=int(np.isin(np.where(m)[0], sub).sum()),
                             rg=float(np.mean([L.rg(w) for w in Wp[m]])),
                             ss=collections.Counter(L.ss_str(cass[pool[m]]).replace("", "")).most_common(3)))
    cl_stats.sort(key=lambda d: -d["size"])
    out["pool_nclusters"] = ncl; out["pool_clusters_ge5"] = cl_stats[:8]
    near = [c for c in cl_stats if c["o_rr_min"] <= thr]
    out["pool_has_near_cluster"] = bool(near)
    out["top75_in_near_cluster"] = int(sum(c["n_top75"] for c in near))
    # top-75
    Wt = Wp[sub]; Pt = Ppair[np.ix_(sub, sub)]
    ss_t = cass[pool[sub]]
    comp = {k: float((ss_t == v).mean()) for k, v in L.SSMAP.items()}
    ncl_t, _ = L.cluster_count(Pt, 2.5)
    rg_t = np.array([L.rg(w) for w in Wt])
    out["top75"] = dict(ss_comp=comp, ss_agree=float((ss_t == o_nss[None, :]).mean()), pair_rmsd_mean=float(Pt[np.triu_indices(len(sub), 1)].mean()),
                        nclusters=ncl_t, rg_mean=float(rg_t.mean()), rg_sd=float(rg_t.std()), o_rr_min=float(rr[pool[sub]].min()),
                        o_rr_mean=float(rr[pool[sub]].mean()), org_frac=float(u["org"][pool[sub]].mean()),
                        ss_mode=collections.Counter(L.ss_str(s) for s in ss_t).most_common(1)[0][0],
                        rg_z_native=float((o_nrg - rg_t.mean()) / max(rg_t.std(), 1e-6)))
    out["native_ss_comp"] = {k: float((o_nss == v).mean()) for k, v in L.SSMAP.items()}
    # distogram predicted rg
    ex = dg["expected"]; d2 = (ex ** 2).sum() + (n - 1) * 3.8 ** 2
    out["disto_rg"] = float(np.sqrt(d2 / n ** 2))
    # diagnosis
    bu = out["band_univ"]
    out["diag"] = dict(
        query=bool(len(band_u) and (rr.min() <= pool_best - 0.5) and bu["in_pool_frac"] < 0.2),
        filter=bool(len(band_p) and out["band_pool"]["score_pct_mean"] > 0.5),
        multimodal=bool(out["pool_has_near_cluster"] and ncl >= 2 and out["top75_in_near_cluster"] < 10),
        ss_mispred=bool(out["disto_ss_agree"] < 0.5),
        rg_mismatch=bool(abs(out["top75"]["rg_z_native"]) > 2.0),
        reference=bool(out["header"]["lasso"] or out["header"]["cyclic"] or out["header"]["fibril"] or out["header"]["xray"]
                       or out["header"]["cosolvent"] or out["header"]["membrane"] or out["header"]["bound"]),
        covalent=bool(out["header"]["covalent"]))
    return out

def aggregate(rows):
    F = [r for r in rows if r["fail18"]]; O = [r for r in rows if not r["fail18"]]
    def col(rs, path):
        v = []
        for r in rs:
            x = r
            try:
                for k in path.split("."):
                    x = x[k]
                v.append(float(x))
            except Exception:
                pass
        return np.array(v)
    keys = ["n", "nw", "o_native_rg", "o_univ_best", "o_pool_best", "o_cov10", "o_cov15", "o_cov20", "o_cov10_pool", "o_cov15_pool", "o_cov20_pool",
            "n_band_univ", "n_band_pool", "shipped", "rmsd_arm", "rho_pool", "disto_ss_agree", "disto_rg",
            "band_univ.org_frac", "band_univ.rank_med", "band_univ.in_pool_frac", "band_univ.ident_mean", "band_univ.esm_pct_med",
            "band_univ.ss_agree_mean", "band_univ.hyd_agree_mean", "univ_top10.rank_med", "univ_top10.in_pool_frac", "univ_top10.org_frac",
            "bg.ident_mean", "bg.ss_agree_mean", "bg.ss_agree_pool", "bg.rg_univ_pct_native", "band_pool.score_pct_mean", "band_pool.score_pct_min",
            "band_pool.ss_agree_mean", "pool_nclusters", "top75_in_near_cluster", "top75.ss_agree", "top75.pair_rmsd_mean", "top75.nclusters",
            "top75.rg_z_native", "top75.o_rr_min", "top75.org_frac", "native_ss_comp.H", "native_ss_comp.E", "native_ss_comp.C",
            "top75.ss_comp.H", "top75.ss_comp.E", "top75.ss_comp.C"]
    table = {}
    for k in keys:
        a, b = col(F, k), col(O, k)
        if len(a) and len(b):
            p = mannwhitneyu(a, b).pvalue
            table[k] = dict(fail18_mean=float(a.mean()), fail18_med=float(np.median(a)), other_mean=float(b.mean()), other_med=float(np.median(b)), p=float(p), n=[len(a), len(b)])
    flags = {}
    for fl in ["lasso", "cyclic", "fibril", "xray", "membrane", "cosolvent", "bound", "designed", "covalent", "nonaqueous", "any_flag"]:
        a = sum(r["header"][fl] for r in F); b = sum(r["header"][fl] for r in O)
        flags[fl] = dict(fail18=f"{a}/18", other=f"{b}/{len(O)}", odds_p=float(fisher_exact([[a, 18 - a], [b, len(O) - b]])[1]))
    for fl in ["query", "filter", "multimodal", "ss_mispred", "rg_mismatch", "reference", "covalent"]:
        a = sum(r["diag"][fl] for r in F); b = sum(r["diag"][fl] for r in O)
        flags["diag_" + fl] = dict(fail18=f"{a}/18", other=f"{b}/{len(O)}", odds_p=float(fisher_exact([[a, 18 - a], [b, len(O) - b]])[1]))
    return dict(columns=table, flags=flags)

if __name__ == "__main__":
    tg = I.targets(); rows = []; t0 = time.time()
    for k, t in enumerate(tg):
        rows.append(analyse(t))
        print(f"[{k+1:3d}/126] {t['pdb']} {'FAIL' if rows[-1]['fail18'] else '    '} {time.time()-t0:.0f}s", flush=True)
    agg = aggregate(rows)
    I.write("forensics_part1", dict(rows=rows, aggregate=agg))
    print(json.dumps(agg["flags"], indent=1))
    for k, v in agg["columns"].items():
        print(f"{k:32s} F18 {v['fail18_mean']:8.3f} ({v['fail18_med']:7.3f})  O108 {v['other_mean']:8.3f} ({v['other_med']:7.3f})  p={v['p']:.3g}")
