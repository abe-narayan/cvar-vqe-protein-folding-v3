"""forensics Part 2b': is the structural CLASS a usable retrieval key?

Section 4 showed the near-native windows are 18.8x enriched in same-class peptide parents for
lasso targets and 6.6x for fibril targets.  This asks whether that is worth anything through the
full chain, at matched K=500.

Arms (all keep K=500 and the shipped score/top-75/average/project downstream):
  blosum        control (production)
  class_meta    BLOSUM + bonus for windows whose PEPTIDE parent shares a header class flag with
                the target.  The target's own class comes from its PDB header text (metadata:
                available for a real target from the experiment, not from its coordinates).
                Labelled METADATA, not fully deployable.
  class_only    rank by (shares class, then BLOSUM) -- the class as the primary key.
  o_shape       ORACLE/DIAGNOSTIC: BLOSUM + 3*(CA-SS agreement with the native) + 3*(|rg-rg_nat|<1)
                -- bounds what knowing the target's SHAPE CLASS (SS string + size) is worth AS A KEY.
  o_shape_only  ORACLE/DIAGNOSTIC: rank by shape agreement alone.
Leakage: `class_meta` uses PDB header TEXT of the target and of library peptides only.  `o_*` arms
are diagnostic bounds and are labelled as such everywhere.
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from s12 import instrument as I
from s12 import forensics_lib as L
from s12.forensics_class import CLASSES, parent_headers

ARMS = ["blosum", "class_meta", "class_only", "o_shape", "o_shape_only"]

def run(t, PH):
    import torch; torch.set_num_threads(2)
    pdb, n = t["pdb"], t["n"]
    u = I.load_univ(pdb); dg = I.distogram(pdb); rr = u["rr"]
    pid, start = L.parent_map(pdb, u); P = L.parents()
    sim = u["sim"].astype(float); nw = len(sim)
    tcls = L.header_class(pdb); mine = {c for c in CLASSES if tcls.get(c)}
    # class match per window (peptide parents only; fragments have no peptide header)
    share = np.zeros(nw)
    ppdb = np.array([str(P["pdb"][k]) for k in pid], dtype=object)
    ispep = np.array([bool(P["is_pep"][k]) for k in pid])
    cache = {}
    for w in range(nw):
        if not ispep[w]:
            continue
        q = ppdb[w]
        if q not in cache:
            h = PH.get(q, {})
            cache[q] = len(mine & {c for c in CLASSES if h.get(c)})
        share[w] = cache[q]
    cass = L.window_cass(pdb, u)
    o_nss = L.ca_ss(u["nat_ca"]); o_rgn = L.rg(u["nat_ca"])            # ORACLE
    ssag = (cass == o_nss[None, :]).sum(1).astype(float)
    rgw = np.sqrt(((u["W"] - u["W"].mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    shape = ssag + 3.0 * (np.abs(rgw - o_rgn) < 1.0)
    keys = {"blosum": sim, "class_meta": sim + 4.0 * share, "class_only": share * 1e4 + sim,
            "o_shape": sim + 3.0 * shape, "o_shape_only": shape * 1e3 + sim}
    thr = float(rr[I.pool_idx(u)].min()) + I.BAND
    out = {"n_class_windows": int((share > 0).sum()), "tclass": sorted(mine)}
    for k, v in keys.items():
        pool = L.rank_key(v)
        r = L.chain(u, pool, dg, t["seq"], t["fold"], do_project=False)
        r["in_band"] = int((rr[pool] <= thr).sum()); r.pop("top_idx_univ")
        out[k] = r
    return pdb, out

if __name__ == "__main__":
    PH = parent_headers(); tg = I.targets(); rows = {}; t0 = time.time()
    path = os.path.join(I.RESULTS, "forensics_classkey.json")
    if os.path.exists(path):
        rows = json.load(open(path)).get("rows", {})
    for k, t in enumerate(tg):
        if t["pdb"] in rows:
            continue
        pdb, r = run(t, PH); rows[pdb] = r
        print(f"[{k+1:3d}/126] {pdb} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
    names = [t["pdb"] for t in tg]; folds = [t["fold"] for t in tg]
    F = np.array([p in I.FAIL18 for p in names])
    agg = {"arms": {}, "paired_vs_blosum": {}}
    for a in ARMS:
        d = {m: np.array([rows[p][a][m] for p in names], float) for m in ("pool_best", "pool_mean", "in_band", "top_best", "rmsd_avg")}
        agg["arms"][a] = {m: dict(all=float(v.mean()), fail18=float(v[F].mean()), other=float(v[~F].mean())) for m, v in d.items()}
        if a == "blosum":
            continue
        agg["paired_vs_blosum"][a] = {}
        for m in ("pool_best", "in_band", "rmsd_avg"):
            x = np.array([rows[p][a][m] for p in names], float); y = np.array([rows[p]["blosum"][m] for p in names], float)
            pr = I.paired(x, y, folds=folds, names=names); pr.pop("top10_targets", None)
            pr["fail18"] = {q: I.paired(x[F], y[F])[q] for q in ("mean_a", "mean_b", "mean_diff", "ci95", "n_better", "n_worse")}
            pr["other"] = {q: I.paired(x[~F], y[~F])[q] for q in ("mean_a", "mean_b", "mean_diff", "ci95", "n_better", "n_worse")}
            agg["paired_vs_blosum"][a][m] = pr
    I.write("forensics_classkey", dict(rows=rows, aggregate=agg))
    for a in ARMS:
        v = agg["arms"][a]
        print(f"{a:14s} " + "  ".join(f"{m}={v[m]['all']:.3f}/{v[m]['fail18']:.3f}/{v[m]['other']:.3f}" for m in ("pool_best", "in_band", "rmsd_avg")))
    for a, d in agg["paired_vs_blosum"].items():
        for m, pr in d.items():
            print(f"{a:14s} {m:10s} d={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']} | F18 {pr['fail18']['mean_diff']:+.3f} | O108 {pr['other']['mean_diff']:+.3f} | drop10 {pr['drop_top10_mean_diff']:+.3f}")
