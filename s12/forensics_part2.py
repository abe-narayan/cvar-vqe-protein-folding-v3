"""forensics Part 2: retrieval keys at matched K=500 with the full downstream chain.

For each key: universe ranking -> top-500 pool -> shipped score top-75 -> coordinate
average -> project.  Evaluation: pool_best, pool_mean, in-band count (band = SHIPPED pool
best + 1.5 A, a fixed per-target reference), top-75 best, rmsd_avg, rmsd_arm.

Keys (deployable unless marked ORACLE):
  blosum        shipped BLOSUM62 sum (reproduces production)
  esm           mean per-position cosine of ESM-2 pca128 (target residue vs window residue)
  ss_pred       agreement of window CA-SS with the shipped distogram's implied SS (random tie-break)
  hyb1/3/10     blosum + lam * (number of residues whose CA-SS agrees with the distogram SS)
  oss           ORACLE agreement with the native CA-SS (random tie-break)
  ohyb3         ORACLE blosum + 3 * agreement with native CA-SS
  score_univ    shipped distogram score over the WHOLE universe (lower better)
  alpha_pred    predicted structural-alphabet string agreement (from forensics_alphabet cache)
  ahyb3         blosum + 3 * alpha_pred agreement
  oalpha        ORACLE native alphabet agreement
  oahyb3        ORACLE blosum + 3 * native alphabet agreement
Usage: python -m s12.forensics_part2 [--keys k1,k2,...] [--workers 2] [--no-project]
"""
import os, sys, time, json, argparse
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from s12 import instrument as I
from s12 import forensics_lib as L

ALL_KEYS = ["blosum", "esm", "ss_pred", "hyb1", "hyb3", "hyb10", "oss", "ohyb3", "score_univ",
            "alpha_pred", "ahyb3", "alpha_ll", "allhyb3", "ss_esm", "sshyb3", "oalpha", "oahyb3",
            "union", "rrf"]
OUT = os.path.join(I.RESULTS, "forensics_part2_rows")
os.makedirs(OUT, exist_ok=True)

def jitter(n, seed):
    return np.random.default_rng(seed).random(n) * 1e-3

def keys_for(u, dg, pdb, want):
    n = u["n"]; nw = len(u["W"]); sim = u["sim"].astype(float)
    cass = L.window_cass(pdb, u)
    out = {}
    if "blosum" in want:
        out["blosum"] = sim
    if any(k in want for k in ("esm", "union", "rrf")):
        e = L.esm_key_cached(pdb, u).astype(float)
        out["esm"] = e
        rb = np.empty(nw); rb[np.argsort(-sim, kind="stable")] = np.arange(nw)
        re_ = np.empty(nw); re_[np.argsort(-e, kind="stable")] = np.arange(nw)
        out["rrf"] = 1.0 / (60 + rb) + 1.0 / (60 + re_)          # reciprocal rank fusion
        un = np.zeros(nw); un[rb < I.K // 2] = 2.0; un[re_ < I.K // 2] += 1.0
        out["union"] = un - 1e-6 * np.minimum(rb, re_)           # interleaved union of the two halves
    if any(k in want for k in ("ss_pred", "hyb1", "hyb3", "hyb10")):
        dss = L.distogram_ss(dg, n)
        agree = (cass == dss[None, :]).sum(1).astype(float)
        out["ss_pred"] = agree + jitter(nw, 1)
        for lam in (1, 3, 10):
            out[f"hyb{lam}"] = sim + lam * agree
    if any(k in want for k in ("oss", "ohyb3")):
        o_nss = L.ca_ss(u["nat_ca"])                     # ORACLE
        agree = (cass == o_nss[None, :]).sum(1).astype(float)
        out["oss"] = agree + jitter(nw, 2)
        out["ohyb3"] = sim + 3 * agree
    if "score_univ" in want:
        i, j = I.pair_index(n)
        sc = np.empty(nw)
        for a in range(0, nw, 4000):
            D = I.pair_dists(u["W"][a:a + 4000], i, j).astype(np.float32).astype(float)
            sc[a:a + 4000] = I.shipped_score(dg, D)
        out["score_univ"] = -sc
    if any(k in want for k in ("alpha_pred", "ahyb3", "oalpha", "oahyb3", "alpha_ll", "allhyb3", "ss_esm", "sshyb3")):
        path = os.path.join(I.CACHE, f"forensics_alpha_{pdb}.npz")
        if os.path.exists(path):
            z = np.load(path)
            ap = z["agree_pred"].astype(float); ao = z["agree_oracle"].astype(float); ll = z["ll_norm"].astype(float); se = z["ss_esm_agree"].astype(float)
            out["alpha_pred"] = ap + jitter(nw, 3); out["ahyb3"] = sim + 3 * ap
            out["oalpha"] = ao + jitter(nw, 4); out["oahyb3"] = sim + 3 * ao
            out["alpha_ll"] = ll; out["allhyb3"] = sim + 3 * ll
            out["ss_esm"] = se + jitter(nw, 5); out["sshyb3"] = sim + 3 * se
    return {k: v for k, v in out.items() if k in want}

def auc(key, pos):
    """AUC of key (higher = retrieved earlier) for positives `pos` (bool)."""
    from scipy.stats import rankdata
    r = rankdata(key); npos = pos.sum(); nneg = len(pos) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    return float((r[pos].sum() - npos * (npos + 1) / 2) / (npos * nneg))

def run_target(args):
    t, want, do_project = args
    import torch; torch.set_num_threads(2)
    pdb = t["pdb"]
    path = os.path.join(OUT, f"{pdb}.json")
    done = json.load(open(path)) if os.path.exists(path) else {}
    need = [k for k in want if k not in done or ("rmsd_arm" not in done[k] and do_project)]
    if not need:
        return pdb, done
    u = I.load_univ(pdb); dg = I.distogram(pdb); rr = u["rr"]
    thr = float(rr[I.pool_idx(u)].min()) + I.BAND
    keys = keys_for(u, dg, pdb, need)
    for k in need:
        if k not in keys:
            continue
        pool = L.rank_key(keys[k])
        r = L.chain(u, pool, dg, t["seq"], t["fold"], do_project=do_project)
        r["in_band"] = int((rr[pool] <= thr).sum()); r["thr"] = thr
        r["auc_band"] = auc(keys[k], rr <= thr); r["auc_20"] = auc(keys[k], rr <= 2.0)
        band = np.where(rr <= thr)[0]
        rank = np.empty(len(rr), int); rank[np.argsort(-keys[k], kind="stable")] = np.arange(len(rr))
        r["band_rank_med"] = float(np.median(rank[band])) if len(band) else None
        r["org_frac"] = float(u["org"][pool].mean())
        r["overlap_blosum"] = float(np.isin(pool, I.pool_idx(u)).mean())
        del r["top_idx_univ"]
        done[k] = r
    with open(path, "w") as fh:
        json.dump(done, fh)
    return pdb, done

def aggregate(rows, want):
    tg = I.targets(); names = [t["pdb"] for t in tg]; folds = [t["fold"] for t in tg]
    F = np.array([p in I.FAIL18 for p in names])
    res = {"keys": {}, "paired_vs_blosum": {}, "error_corr": {}}
    metrics = ["pool_best", "pool_mean", "in_band", "top_best", "top_mean", "rmsd_avg", "rmsd_arm", "auc_band", "auc_20", "org_frac", "overlap_blosum"]
    vec = {}
    for k in want:
        if not all(k in rows[p] for p in names):
            continue
        d = {m: np.array([rows[p][k].get(m, np.nan) for p in names], float) for m in metrics}
        vec[k] = d
        res["keys"][k] = {m: dict(all=float(np.nanmean(v)), fail18=float(np.nanmean(v[F])), other=float(np.nanmean(v[~F]))) for m, v in d.items()}
    if "blosum" in vec:
        for k in vec:
            if k == "blosum":
                continue
            res["paired_vs_blosum"][k] = {}
            for m in ["pool_best", "pool_mean", "in_band", "top_best", "rmsd_avg", "rmsd_arm"]:
                a, b = vec[k][m], vec["blosum"][m]
                if np.isnan(a).any() or np.isnan(b).any():
                    continue
                pr = I.paired(a, b, folds=folds, names=names)
                pr["fail18"] = I.paired(a[F], b[F]); pr["other"] = I.paired(a[~F], b[~F])
                for q in ("fail18", "other"):
                    pr[q] = {x: pr[q][x] for x in ("mean_a", "mean_b", "mean_diff", "ci95", "n_better", "n_worse")}
                pr.pop("top10_targets", None)
                res["paired_vs_blosum"][k][m] = pr
        ks = [k for k in vec if not np.isnan(vec[k]["rmsd_arm"]).any()]
        E = np.array([vec[k]["rmsd_arm"] for k in ks])
        if len(ks) > 1:
            res["error_corr"] = {"keys": ks, "pearson": np.corrcoef(E).tolist(),
                                 "spearman": __import__("scipy.stats", fromlist=["spearmanr"]).spearmanr(E, axis=1).correlation.tolist(),
                                 "min_over_keys_mean": float(E.min(0).mean()),
                                 "min_over_keys_fail18": float(E[:, F].min(0).mean())}
    return res

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--keys", default=",".join(ALL_KEYS)); ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--no-project", action="store_true"); ap.add_argument("--tag", default="")
    a = ap.parse_args(); want = a.keys.split(",")
    tg = I.targets(); rows = {}; t0 = time.time()
    if a.workers <= 1:
        for k, t in enumerate(tg):
            pdb, r = run_target((t, want, not a.no_project)); rows[pdb] = r
            print(f"[{k+1:3d}/126] {pdb} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
    else:
        import multiprocessing as mp
        with mp.Pool(a.workers) as pool:
            for k, (pdb, r) in enumerate(pool.imap_unordered(run_target, [(t, want, not a.no_project) for t in tg])):
                rows[pdb] = r
                print(f"[{k+1:3d}/126] {pdb} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
    agg = aggregate(rows, want)
    I.write("forensics_part2" + a.tag, agg)
    for k, v in agg["keys"].items():
        print(f"{k:12s} " + "  ".join(f"{m}={v[m]['all']:.3f}/{v[m]['fail18']:.3f}/{v[m]['other']:.3f}" for m in ("pool_best", "pool_mean", "in_band", "top_best", "rmsd_avg", "rmsd_arm")))
