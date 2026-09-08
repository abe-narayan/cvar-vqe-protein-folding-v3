"""forensics Part 2d stage 2: the size gate with a DEPLOYABLE predicted size.

Arms (all: gate the WHOLE universe by |rg_window - rg_hat| < w, then BLOSUM top-500, then the
shipped chain).  A target with fewer than K=500 surviving windows falls back to the ungated
BLOSUM pool (reported).
  base        production BLOSUM top-500 (control)
  o_gate05    ORACLE rg_hat = native rg, w = 0.5   (the diagnostic bound from forensics_compact)
  gbt_gate10  DEPLOYABLE rg_hat = LFO GBT prediction, w = 1.0 (= its own MAE)
  gbt_gate15  DEPLOYABLE, w = 1.5
  disto_gate10 DEPLOYABLE rg_hat = shipped distogram's implied rg, w = 1.0
  const_gate10 CONTROL rg_hat = training-set mean rg, w = 1.0  (isolates "any size gate helps")
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from s12 import instrument as I
from s12 import forensics_lib as L

ARMS = [("base", None, None), ("o_gate05", "o_true", 0.5), ("gbt_gate10", "gbt", 1.0),
        ("gbt_gate15", "gbt", 1.5), ("disto_gate10", "disto", 1.0), ("const_gate10", "const", 1.0),
        # width-matched decomposition: at the SAME w=1.0, how much of the gain is target-specific?
        ("o_gate10", "o_true", 1.0), ("shuf_gate10", "shuf", 1.0)]

def main():
    z = np.load(os.path.join(I.CACHE, "forensics_rgpred.npz"), allow_pickle=True)
    P = {k: z[k] for k in z.files}
    idx = {str(p): k for k, p in enumerate(P["pdb"])}
    # NULL control: each target gets ANOTHER target's true rg (derangement) -- same marginal
    # distribution of gate centres, zero target-specific information.
    rng = np.random.default_rng(0); n = len(P["o_true"]); perm = rng.permutation(n)
    while (perm == np.arange(n)).any():
        perm = rng.permutation(n)
    P["shuf"] = P["o_true"][perm]
    tg = I.targets(); rows = {}; t0 = time.time()
    path = os.path.join(I.RESULTS, "forensics_rggate.json")
    if os.path.exists(path):
        rows = json.load(open(path)).get("rows", {})
    import torch; torch.set_num_threads(2)
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in rows:
            continue
        u = I.load_univ(pdb); dg = I.distogram(pdb); rr = u["rr"]
        rgw = np.sqrt(((u["W"] - u["W"].mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        thr = float(rr[I.pool_idx(u)].min()) + I.BAND
        r = {}
        for name, src, w in ARMS:
            if src is None:
                pool = I.pool_idx(u); fell = False
            else:
                c = float(P[src][idx[pdb]])
                keep = np.where(np.abs(rgw - c) < w)[0]
                fell = len(keep) < I.K
                pool = I.pool_idx(u) if fell else keep[np.argsort(-u["sim"][keep], kind="stable")[:I.K]]
            o = L.chain(u, pool, dg, t["seq"], t["fold"], do_project=False); o.pop("top_idx_univ")
            o["in_band"] = int((rr[pool] <= thr).sum()); o["fellback"] = bool(fell)
            r[name] = o
        rows[pdb] = r
        print(f"[{k+1:3d}/126] {pdb} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
        I.write("forensics_rggate", dict(rows=rows))
    names = [t["pdb"] for t in tg]; folds = [t["fold"] for t in tg]
    F = np.array([p in I.FAIL18 for p in names])
    agg = {"arms": {}, "paired_vs_base": {}}
    for name, _, _ in ARMS:
        d = {m: np.array([rows[p][name][m] for p in names], float) for m in ("pool_best", "pool_mean", "in_band", "top_best", "rmsd_avg")}
        agg["arms"][name] = {m: dict(all=float(v.mean()), fail18=float(v[F].mean()), other=float(v[~F].mean())) for m, v in d.items()}
        agg["arms"][name]["n_fellback"] = int(sum(rows[p][name]["fellback"] for p in names))
        if name == "base":
            continue
        agg["paired_vs_base"][name] = {}
        for m in ("pool_best", "in_band", "rmsd_avg"):
            x = np.array([rows[p][name][m] for p in names], float); y = np.array([rows[p]["base"][m] for p in names], float)
            pr = I.paired(x, y, folds=folds, names=names); pr.pop("top10_targets", None)
            pr["fail18"] = {q: I.paired(x[F], y[F])[q] for q in ("mean_a", "mean_b", "mean_diff", "ci95", "n_better", "n_worse")}
            pr["other"] = {q: I.paired(x[~F], y[~F])[q] for q in ("mean_a", "mean_b", "mean_diff", "ci95", "n_better", "n_worse")}
            agg["paired_vs_base"][name][m] = pr
    I.write("forensics_rggate", dict(rows=rows, aggregate=agg))
    for name, _, _ in ARMS:
        v = agg["arms"][name]
        print(f"{name:14s} fb={v['n_fellback']:3d} " + "  ".join(f"{m}={v[m]['all']:.3f}/{v[m]['fail18']:.3f}/{v[m]['other']:.3f}" for m in ("pool_best", "in_band", "rmsd_avg")))
    for name, d in agg["paired_vs_base"].items():
        for m in ("rmsd_avg",):
            pr = d[m]
            print(f"{name:14s} {m} d={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']} drop10={pr['drop_top10_mean_diff']:+.3f} | F18 {pr['fail18']['mean_diff']:+.3f} CI[{pr['fail18']['ci95'][0]:+.3f},{pr['fail18']['ci95'][1]:+.3f}] | O108 {pr['other']['mean_diff']:+.3f}")

if __name__ == "__main__":
    main()
