"""E5 -- ranking diagnostics INSIDE the assembled candidate space, and the diversity question.

For each target, enumerate the same assembled candidate space as E3 (key `combo` or `random`,
k pieces per position, all 2- and 3-piece structures) and measure:
  rho_global : Spearman(shipped distogram score, true CA-RMSD) over all assembled candidates
  rho_inband : the same restricted to candidates within (best + 1.5 A)  -- "in-band" metric
  native_pct : percentile of the NATIVE's own distogram score among the assembled candidates
  recall75   : does the distogram top-75 of the assembled pool contain a band member?
  best_k     : best rr in the assembled pool as a function of shortlist size k (4, 8, 20)
and compares them against the same quantities computed on the retrieval K=500 pool.
All rr values are ORACLE evaluation.
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_bank as AB
from s12 import assembly_common as AC
from s12 import assembly_e2 as E2
from s12 import assembly_e3 as E3

KEYS = ["combo", "esm", "random"]
KS = [4, 8, 20]


def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def run_target(t):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    out = os.path.join(I.CACHE, f"asm_e5_{pdb}.json")
    if os.path.exists(out):
        return json.load(open(out))
    u = I.load_univ(pdb); nat = u["nat_ca"]; p = I.pool_idx(u)
    dg = I.distogram(pdb); i, j = I.pair_index(n)
    K = E2.target_keys(t)
    import peptide_db as db
    q = db.by_pdb(pdb)
    nat_ideal = I.build_ca(q.phi[None], q.psi[None])[0]
    nat_score = float(I.shipped_score(dg, I.pair_dists(nat_ideal[None], i, j))[0])
    res = {"pdb": pdb, "n": n, "fold": fold, "native_score": nat_score, "keys": {}}
    # retrieval control
    Wp = u["W"][p]; sp = I.shipped_score(dg, I.pair_dists(Wp, i, j).astype(np.float32).astype(float)); rp = u["rr"][p]
    bandp = rp <= rp.min() + I.BAND
    res["retrieval"] = {"rho_global": spearman(sp, rp), "rho_inband": spearman(sp[bandp], rp[bandp]) if bandp.sum() > 2 else None,
                        "native_pct": float((sp < nat_score).mean() * 100), "band_size": int(bandp.sum()),
                        "recall75": bool(np.isin(np.argsort(sp, kind="stable")[:I.M], np.where(bandp)[0]).any()),
                        "best": float(rp.min())}
    for key in KEYS:
        E3.KEY = key; E3._TOP = {}
        best_by_k = {}
        for k in KS:
            allPHI, allPSI, tag, meta = [], [], [], []
            for comp in E3.structures(n):
                PHI, PSI, grids, tops, junc = E3.enumerate_structure(t, K, comp, k)
                allPHI.append(PHI); allPSI.append(PSI); tag += [len(meta)] * len(PHI); meta.append(comp)
            PHI = np.concatenate(allPHI); PSI = np.concatenate(allPSI); tag = np.array(tag)
            CA = AC.build_many(PHI, PSI)
            sc = I.shipped_score(dg, I.pair_dists(CA, i, j)); rr = I.kabsch_rmsd_batch(CA, nat)
            best_by_k[k] = float(rr.min())
            if k == max(KS):
                band = rr <= rr.min() + I.BAND
                res["keys"][key] = {"N": int(len(rr)), "rho_global": spearman(sc, rr),
                                    "rho_inband": spearman(sc[band], rr[band]) if band.sum() > 2 else None,
                                    "native_pct": float((sc < nat_score).mean() * 100), "band_size": int(band.sum()),
                                    "recall75": bool(np.isin(np.argsort(sc, kind="stable")[:I.M], np.where(band)[0]).any()),
                                    "best": float(rr.min()), "argmin_rr": float(rr[int(np.argmin(sc))]),
                                    "top75_best": float(rr[np.argsort(sc, kind="stable")[:I.M]].min()),
                                    "median": float(np.median(rr)),
                                    "diversity": float(np.mean(I.pairwise_rmsd(CA[np.argsort(sc, kind="stable")[:40]])))}
        res["keys"][key]["best_by_k"] = best_by_k
    json.dump(res, open(out, "w"), default=float)
    return json.load(open(out))


def _init():
    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        import torch; torch.set_num_threads(1)
    except Exception:
        pass


if __name__ == "__main__":
    import multiprocessing as mp
    tg = I.targets(); rows = []
    print("free GB", I.free_gb(), flush=True)
    with mp.Pool(2, initializer=_init) as pool:
        for c, r in enumerate(pool.imap_unordered(run_target, tg)):
            rows.append(r)
            if c % 10 == 0:
                print(f"[{c+1}/126] {r['pdb']} free={I.free_gb():.1f}", flush=True)
    rows.sort(key=lambda r: r["pdb"]); pdbs = [r["pdb"] for r in rows]
    G = lambda v: AC.group_means(v, pdbs)
    agg = {"n": len(rows)}
    fl = lambda vals: [v for v in vals if v is not None and np.isfinite(v)]
    agg["retrieval"] = {"rho_global": float(np.mean(fl([r["retrieval"]["rho_global"] for r in rows]))),
                        "rho_inband": float(np.mean(fl([r["retrieval"]["rho_inband"] for r in rows]))),
                        "native_pct": G([r["retrieval"]["native_pct"] for r in rows]),
                        "recall75": float(np.mean([r["retrieval"]["recall75"] for r in rows])),
                        "best": G([r["retrieval"]["best"] for r in rows])}
    for key in KEYS:
        d = [r["keys"][key] for r in rows]
        agg[key] = {"N_mean": float(np.mean([x["N"] for x in d])),
                    "rho_global": float(np.mean(fl([x["rho_global"] for x in d]))),
                    "rho_inband": float(np.mean(fl([x["rho_inband"] for x in d]))),
                    "native_pct": G([x["native_pct"] for x in d]),
                    "recall75": float(np.mean([x["recall75"] for x in d])),
                    "best": G([x["best"] for x in d]), "argmin_rr": G([x["argmin_rr"] for x in d]),
                    "top75_best": G([x["top75_best"] for x in d]), "median": G([x["median"] for x in d]),
                    "diversity": G([x["diversity"] for x in d]),
                    "best_by_k": {k: G([x["best_by_k"][str(k)] for x in d]) for k in KS}}
        agg[f"paired_{key}_best_vs_retrieval_best"] = I.paired([x["best"] for x in d], [r["retrieval"]["best"] for r in rows],
                                                               folds=[r["fold"] for r in rows], names=pdbs)
    agg["per_target"] = [{"pdb": r["pdb"], "fail18": r["pdb"] in I.FAIL18, "retr_best": r["retrieval"]["best"],
                          **{f"{k}_best": r["keys"][k]["best"] for k in KEYS},
                          **{f"{k}_argmin": r["keys"][k]["argmin_rr"] for k in KEYS}} for r in rows]
    print(I.write("assembly_e5_ranking", agg))
    print("retrieval", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in agg["retrieval"].items()})
    for key in KEYS:
        print(key, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in agg[key].items() if k != "best_by_k"})
        print("   best_by_k", {k: round(v["all"], 3) for k, v in agg[key]["best_by_k"].items()})
