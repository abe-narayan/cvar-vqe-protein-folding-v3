"""EXPERIMENT 1 -- the terminal-operator ladder.

Stage `raw`  : builds, for every target, a candidate structure C for every operator arm and
               records the RAW CA-RMSD of C to the native.  Dumps the clouds so the (slow)
               projection can be run separately / in parallel.
Stage `proj` : projects a list of arms and records the EMITTED CA-RMSD (lam=0 fit_ca arm,
               exactly the production synthesis path that gives 3.2041 for `avg75`).

Usage:
    python -m s12.agg_ladder raw
    python -m s12.agg_ladder proj <arm1,arm2,...> <shard> <nshard>
    python -m s12.agg_ladder collect
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_common as A

CLOUD = os.path.join(ROOT, "s12", "cache", "agg_clouds")
PROJ = os.path.join(ROOT, "s12", "cache", "agg_proj")
os.makedirs(CLOUD, exist_ok=True)
os.makedirs(PROJ, exist_ok=True)


def _z(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / (s if s > 1e-12 else 1.0)


def simplex_lsq(A, T, iters=400):
    """min ||sum_k w_k A_k - T||_F  s.t. w >= 0, sum w = 1.  Projected gradient (ORACLE)."""
    m = len(A)
    X = A.reshape(m, -1)
    t = T.reshape(-1)
    G = X @ X.T
    c = X @ t
    L = np.linalg.eigvalsh(G).max() + 1e-9
    w = np.ones(m) / m
    for _ in range(iters):
        g = G @ w - c
        w = _proj_simplex(w - g / L)
    return w


def _proj_simplex(v):
    u = np.sort(v)[::-1]
    css = np.cumsum(u)
    rho = np.nonzero(u * np.arange(1, len(v) + 1) > (css - 1))[0][-1]
    theta = (css[rho] - 1.0) / (rho + 1.0)
    return np.maximum(v - theta, 0.0)


def cluster_labels(P, k):
    """Average-linkage agglomerative clustering on the pairwise-RMSD matrix."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    if len(P) <= k:
        return np.arange(len(P))
    Z = linkage(squareform(P, checks=False), method="average")
    return fcluster(Z, t=k, criterion="maxclust") - 1


def arms_for(d, pdb, seq, n):
    """Return {arm_name: C (n,3)} plus diagnostics."""
    Wp = d["Wp"]; rr = d["rr"].astype(float); sc = d["sc"].astype(float); sub = d["sub"].astype(int)
    out = {}
    diag = {}

    W75 = Wp[sub]; rr75 = rr[sub]; sc75 = sc[sub]
    P75 = I.pairwise_rmsd(W75)
    med75 = int(np.argmin(P75.mean(1)))
    ref = W75[med75]
    A75 = I.superpose_batch(W75, ref)
    cons75 = P75.mean(1)

    # ---- singletons
    out["argmin75"] = W75[int(np.argmin(sc75))]
    out["medoid75"] = W75[med75]
    # ---- plain averages
    out["avg75"] = A75.mean(0)
    out["iter_avg75"] = A.iterative_average(W75, iters=4)
    for f in (0.1, 0.2, 0.35, 0.5):
        C, keep = A.trimmed_average(W75, frac=f)
        out[f"trim{int(f*100)}_75"] = C
    out["geomed75"] = A.geometric_median(W75)
    # ---- score weightings
    zs = _z(sc75)
    for T in (0.25, 0.5, 1.0, 2.0):
        w = np.exp(-zs / T); out[f"wscore_T{T}_75"] = (A75 * (w / w.sum())[:, None, None]).sum(0)
    rk = np.argsort(np.argsort(sc75)).astype(float)
    for r0 in (5.0, 20.0):
        w = 1.0 / (rk + r0); out[f"wrank_r{int(r0)}_75"] = (A75 * (w / w.sum())[:, None, None]).sum(0)
    # ---- consensus weightings
    zc = _z(cons75)
    for T in (0.5, 1.0):
        w = np.exp(-zc / T); out[f"wcons_T{T}_75"] = (A75 * (w / w.sum())[:, None, None]).sum(0)
    for T in (0.5, 1.0):
        w = np.exp(-(zs + zc) / T); out[f"wboth_T{T}_75"] = (A75 * (w / w.sum())[:, None, None]).sum(0)
    # ---- cluster-then-average (largest mode)
    for k in (2, 3):
        lab = cluster_labels(P75, k)
        sizes = np.bincount(lab, minlength=k)
        big = int(np.argmax(sizes))
        idx = np.where(lab == big)[0]
        out[f"bigmode{k}_75"] = I.superpose_batch(W75[idx], W75[idx][int(np.argmin(P75[np.ix_(idx, idx)].mean(1)))]).mean(0)
        diag[f"modes{k}_sizes"] = sizes.tolist()
        diag[f"modes{k}_labels"] = lab.tolist()
        # per-mode averages + oracle best-of-k
        bo, mrm = [], []
        for c in range(k):
            ix = np.where(lab == c)[0]
            if len(ix) == 0:
                continue
            Cm = I.superpose_batch(W75[ix], W75[ix][int(np.argmin(P75[np.ix_(ix, ix)].mean(1)))]).mean(0)
            out[f"mode{k}_{c}_75"] = Cm
            bo.append(I.ca_rmsd(Cm, d["nat"]))
        diag[f"oracle_bestof{k}_raw"] = float(min(bo))

    # ---- pool-level cardinality ladder (score order within the K=500 pool)
    order = np.argsort(sc)
    for m in (10, 25, 50, 150, 300, 500):
        idx = order[:m]
        Wm = Wp[idx]
        Pm = I.pairwise_rmsd(Wm)
        rm = Wm[int(np.argmin(Pm.mean(1)))]
        out[f"avg{m}_pool"] = I.superpose_batch(Wm, rm).mean(0)

    # ---- ORACLE arms
    out["ORC_best75"] = W75[int(np.argmin(rr75))]
    zo = _z(rr75)
    for T in (0.25, 0.5, 1.0):
        w = np.exp(-zo / T); out[f"ORC_wrmsd_T{T}_75"] = (A75 * (w / w.sum())[:, None, None]).sum(0)
    ordo = np.argsort(rr75)
    for m in (10, 25, 40):
        idx = ordo[:m]
        out[f"ORC_top{m}avg_75"] = I.superpose_batch(W75[idx], W75[idx][0]).mean(0)
    # oracle simplex weights on native-superposed top-75  (reproduces "oracle re-weighting 2.044")
    An = I.superpose_batch(W75, d["nat"])
    w = simplex_lsq(An, d["nat"])
    out["ORC_simplex75"] = (An * w[:, None, None]).sum(0)
    diag["ORC_simplex75_support"] = int((w > 1e-4).sum())
    # oracle simplex on the whole K=500 pool (the 0.853 hull number)
    Anp = I.superpose_batch(Wp, d["nat"])
    wp = simplex_lsq(Anp, d["nat"], iters=600)
    out["ORC_simplex500"] = (Anp * wp[:, None, None]).sum(0)
    # oracle pool best member
    out["ORC_bestpool"] = Wp[int(np.argmin(rr))]
    # oracle top-25 of the POOL averaged (C1's 1.644)
    ordp = np.argsort(rr)
    out["ORC_top25avg_pool"] = I.superpose_batch(Wp[ordp[:25]], Wp[ordp[0]]).mean(0)
    out["ORC_top75avg_pool"] = I.superpose_batch(Wp[ordp[:75]], Wp[ordp[0]]).mean(0)

    diag["cons75"] = cons75.tolist()
    diag["med75"] = med75
    return out, diag


def stage_raw():
    tg = I.targets()
    rows = {}
    diags = {}
    t0 = time.time()
    for k, t in enumerate(tg):
        d = A.load(t["pdb"], mem=False)
        arms, diag = arms_for(d, t["pdb"], t["seq"], t["n"])
        names = sorted(arms)
        Cs = np.stack([arms[a] for a in names]).astype(np.float32)
        np.savez_compressed(os.path.join(CLOUD, f"{t['pdb']}.npz"), C=Cs, names=np.array(names))
        rows[t["pdb"]] = {a: float(I.ca_rmsd(arms[a], d["nat"])) for a in names}
        diags[t["pdb"]] = diag
        if k % 20 == 0:
            print(f"{k} {t['pdb']} {time.time()-t0:.0f}s", flush=True)
    I.write("agg_ladder_raw", {"raw": rows})
    I.write("agg_ladder_diag", diags)
    print("done", time.time() - t0)


def stage_proj(arms, shard, nshard):
    tg = I.targets()
    tg = [t for k, t in enumerate(tg) if k % nshard == shard]
    for t in tg:
        path = os.path.join(PROJ, f"{t['pdb']}.json")
        cur = json.load(open(path)) if os.path.exists(path) else {}
        z = np.load(os.path.join(CLOUD, f"{t['pdb']}.npz"))
        names = list(z["names"]); C = z["C"].astype(float)
        d = A.load(t["pdb"], mem=False)
        changed = False
        for a in arms:
            if a in cur or a not in names:
                continue
            out = I.project(C[names.index(a)], t["seq"], t["fold"])
            cur[a] = {"fit": float(I.ca_rmsd(out["fit_ca"], d["nat"])),
                      "lam": float(I.ca_rmsd(out["ca"], d["nat"]))}
            changed = True
        if changed:
            json.dump(cur, open(path, "w"))
        print(t["pdb"], flush=True)


def stage_collect():
    tg = I.targets()
    out = {}
    for t in tg:
        path = os.path.join(PROJ, f"{t['pdb']}.json")
        out[t["pdb"]] = json.load(open(path)) if os.path.exists(path) else {}
    I.write("agg_ladder_proj", out)
    print("collected")


if __name__ == "__main__":
    if sys.argv[1] == "raw":
        stage_raw()
    elif sys.argv[1] == "proj":
        stage_proj(sys.argv[2].split(","), int(sys.argv[3]), int(sys.argv[4]))
    else:
        stage_collect()
