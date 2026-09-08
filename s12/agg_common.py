"""Shared cache + helpers for the SET-AGGREGATION agent (s12/agg_*).

Builds, once, a compact per-target npz holding everything the operator ladder and the
learned set decoder need, so no experiment ever reloads the 7k-40k window universes.

    s12/cache/agg/<pdb>.npz
        Wp    (500, n, 3) float32   pool coordinates (shipped K=500 BLOSUM pool)
        rr    (500,)      float32   ORACLE CA-RMSD of each pool member to the native
        sc    (500,)      float32   shipped distogram Bayes-risk score (lower = better)
        sub   (75,)       int32     shipped top-75 indices INTO THE POOL
        nat   (n, 3)      float32   ORACLE native CA
        exp   (npairs,)   float32   distogram expected distance
        sd    (npairs,)   float32   distogram sd
        pi,pj (npairs,)   int32     pair index (min_sep 2)
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

AGG = os.path.join(ROOT, "s12", "cache", "agg")
os.makedirs(AGG, exist_ok=True)


def build(pdb=None, verbose=True):
    tg = I.targets()
    if pdb is not None:
        tg = [t for t in tg if t["pdb"] == pdb]
    for t in tg:
        path = os.path.join(AGG, f"{t['pdb']}.npz")
        if os.path.exists(path):
            continue
        u = I.load_univ(t["pdb"])
        p = I.pool_idx(u)
        Wp = u["W"][p]
        rr = u["rr"][p]
        dg = I.distogram(t["pdb"], t["seq"], t["fold"])
        pi, pj = I.pair_index(t["n"])
        D = I.pair_dists(Wp, pi, pj)
        sc = I.shipped_score(dg, D.astype(np.float32).astype(float))
        rec = I.shipped_record(t["pdb"])
        sub = np.asarray(rec["sub"], np.int32)
        np.savez_compressed(path, Wp=Wp.astype(np.float32), rr=rr.astype(np.float32),
                            sc=sc.astype(np.float32), sub=sub, nat=u["nat_ca"].astype(np.float32),
                            exp=dg["expected"].astype(np.float32), sd=dg["sd"].astype(np.float32),
                            pi=pi.astype(np.int32), pj=pj.astype(np.int32))
        if verbose:
            print(t["pdb"], flush=True)
        del u


_MEM = {}


def load(pdb, mem=True):
    if mem and pdb in _MEM:
        return _MEM[pdb]
    z = np.load(os.path.join(AGG, f"{pdb}.npz"))
    d = {k: z[k] for k in z.files}
    d["Wp"] = d["Wp"].astype(np.float64)
    d["nat"] = d["nat"].astype(np.float64)
    if mem:
        _MEM[pdb] = d
    return d


def risk_of(pdb, D):
    """Shipped Bayes-risk score of pair-distance rows D (b, npairs).  Lower = better."""
    dg = I.distogram(pdb)
    return I.shipped_score(dg, np.asarray(D, float))


# ------------------------------------------------------------------ set operators
def superpose_to(W, T):
    return I.superpose_batch(W, T)


def weighted_average(W, w, ref=None, P=None):
    """Superpose W on a reference (medoid by default), then weighted mean."""
    W = np.asarray(W, float)
    w = np.asarray(w, float)
    w = w / w.sum()
    if ref is None:
        if P is None:
            P = I.pairwise_rmsd(W)
        ref = W[int(np.argmin(P.mean(1)))]
    A = I.superpose_batch(W, ref)
    return (A * w[:, None, None]).sum(0)


def iterative_average(W, w=None, iters=3):
    """Superpose-average-resuperpose until stable (a better estimator than one pass)."""
    W = np.asarray(W, float)
    if w is None:
        w = np.ones(len(W))
    w = np.asarray(w, float); w = w / w.sum()
    P = I.pairwise_rmsd(W)
    ref = W[int(np.argmin(P.mean(1)))]
    C = ref
    for _ in range(iters):
        A = I.superpose_batch(W, C)
        C = (A * w[:, None, None]).sum(0)
    return C


def geometric_median(W, ref=None, iters=25, eps=1e-6):
    """Weiszfeld on superposed coordinates, re-superposing each round."""
    W = np.asarray(W, float)
    P = I.pairwise_rmsd(W)
    C = W[int(np.argmin(P.mean(1)))] if ref is None else ref
    for _ in range(iters):
        A = I.superpose_batch(W, C)
        d = np.sqrt(((A - C) ** 2).sum((1, 2)))
        w = 1.0 / np.maximum(d, eps)
        w = w / w.sum()
        Cn = (A * w[:, None, None]).sum(0)
        if np.abs(Cn - C).max() < 1e-5:
            C = Cn
            break
        C = Cn
    return C


def trimmed_average(W, frac=0.2, iters=3):
    """Drop the `frac` most-outlying members (by distance to the running average)."""
    W = np.asarray(W, float)
    P = I.pairwise_rmsd(W)
    C = W[int(np.argmin(P.mean(1)))]
    keep = np.arange(len(W))
    for _ in range(iters):
        A = I.superpose_batch(W, C)
        d = np.sqrt(((A - C) ** 2).sum((1, 2)) / W.shape[1])
        k = max(2, int(round(len(W) * (1 - frac))))
        keep = np.argsort(d)[:k]
        C = A[keep].mean(0)
    return C, keep


if __name__ == "__main__":
    build()
    print("cache ok", len(os.listdir(AGG)))
