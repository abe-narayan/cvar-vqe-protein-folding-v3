"""EXPERIMENT 5 (coordinator request) -- the OPERATOR x OBJECTIVE-QUALITY surface.

Question: is the optimal terminal operator a FUNCTION OF OBJECTIVE QUALITY?  With the
shipped (bad) objective, averaging beats argmin by 0.406 A raw.  With a perfect objective,
argmin is claimed to beat averaging.  If so, "the production average is the native-free
optimum" is true only CONDITIONAL on the current objective.

Objective axis (reuses the objective agent's `interp` corruption family, i.e. the REAL
distogram error structure scaled in amplitude -- not i.i.d. noise):

    dtarget(a) = max( dtrue + a * (exp - dtrue), 1.5 )      a = 0 -> native
                                                            a = 1 -> shipped expected dists
                                                            a > 1 -> worse than shipped
    score(c)   = mean_p | D[c,p] - dtarget(a)_p |           (L1, the objective agent's scorer)

plus the actual SHIPPED Bayes-risk score as a labelled reference point.

Operator axis: top-m coordinate average for m in {1,2,3,5,8,12,20,35,50,75,150,300,500}
(m = 1 is argmin), each averaged after superposing the SELECTED SUBSET on its OWN medoid
(the production operator), plus trimmed means and score-weighted averages at headline m.

Outcome: RAW CA-RMSD over the whole surface (cheap); headline cells confirmed through the
real L-BFGS projection by `agg_surface_proj.py`.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_common as A

MS = [1, 2, 3, 5, 8, 12, 20, 35, 50, 75, 110, 150, 220, 300, 500]
AS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]
PCACHE = os.path.join(ROOT, "s12", "cache", "agg_P500")
os.makedirs(PCACHE, exist_ok=True)


def pool_P(pdb, W):
    f = os.path.join(PCACHE, f"{pdb}.npy")
    if os.path.exists(f):
        return np.load(f)
    P = I.pairwise_rmsd(W).astype(np.float32)
    np.save(f, P)
    return P


def subset_avg(W, P, idx):
    """Production operator on a subset: superpose on the subset's own medoid, mean."""
    idx = np.asarray(idx, int)
    if len(idx) == 1:
        return W[idx[0]]
    sub = P[np.ix_(idx, idx)]
    b = idx[int(np.argmin(sub.mean(1)))]
    return I.superpose_batch(W[idx], W[b]).mean(0)


def run(out="agg_surface"):
    tg = I.targets()
    rows = {}
    t0 = time.time()
    for q, t in enumerate(tg):
        pdb = t["pdb"]
        d = A.load(pdb, mem=False)
        W = d["Wp"]; nat = d["nat"]; rr = d["rr"].astype(float)
        pi, pj = d["pi"].astype(int), d["pj"].astype(int)
        D = I.pair_dists(W, pi, pj)
        dtrue = np.linalg.norm(nat[pi] - nat[pj], axis=-1)
        e = d["exp"].astype(float)
        P = pool_P(pdb, W)
        r = {"n": t["n"], "fold": t["fold"], "pool_best": float(rr.min())}
        # --- shipped Bayes-risk objective (the real one)
        objs = {"shipped": d["sc"].astype(float)}
        maes = {"shipped": float(np.abs(e - dtrue).mean())}
        for a in AS:
            dt = np.maximum(dtrue + a * (e - dtrue), 1.5)
            objs[f"a{a}"] = np.abs(D - dt[None, :]).mean(1)
            maes[f"a{a}"] = float(np.abs(dt - dtrue).mean())
        r["mae"] = maes
        for nm, sc in objs.items():
            order = np.argsort(sc, kind="stable")
            for m in MS:
                idx = order[:m]
                r[f"{nm}|m{m}"] = float(I.ca_rmsd(subset_avg(W, P, idx), nat))
                r[f"{nm}|m{m}|best"] = float(rr[idx].min())
            # trimmed + weighted at the headline cardinalities
            for m in (25, 75):
                idx = order[:m]
                C, keep = A.trimmed_average(W[idx], frac=0.25)
                r[f"{nm}|trim25_m{m}"] = float(I.ca_rmsd(C, nat))
                z = sc[idx]; z = (z - z.mean()) / (z.std() + 1e-9)
                for T in (0.5, 2.0):
                    w = np.exp(-z / T); w /= w.sum()
                    sub = P[np.ix_(idx, idx)]
                    b = idx[int(np.argmin(sub.mean(1)))]
                    Cw = (I.superpose_batch(W[idx], W[b]) * w[:, None, None]).sum(0)
                    r[f"{nm}|w{T}_m{m}"] = float(I.ca_rmsd(Cw, nat))
        rows[pdb] = r
        if q % 10 == 0:
            print(q, pdb, f"{time.time()-t0:.0f}s", flush=True)
    I.write(out, rows)
    print("done", time.time() - t0)


if __name__ == "__main__":
    run()
