"""ADVERSARIAL AUDIT 1b -- "the objective's own optimum is 4.62 A, worse than the 3.20 A the
pipeline emits, therefore the pipeline beats its own objective".

That number is a CLASSICAL MDS embedding of E[d] (`obj_errstruct.py`).  Classical MDS is a
weak decoder: it eigendecomposes a doubly-centred squared-distance matrix that is not even
Euclidean, it ignores the score's shell weights, and above all it minimises the WRONG
FUNCTIONAL -- strain on E[d], not the shipped Bayes-risk score the pipeline actually ranks
with.  If a proper decoder gets much below 4.62 A, the "pipeline beats its own objective"
claim weakens or dies.

Four decoders, same distogram, same 126 targets:
  A. classical MDS on E[d]                       (reproduce 4.62)
  B. classical MDS output pushed through the PROJECT operator (the pipeline's own manifold)
  C. SMACOF stress majorisation on E[d], shell-weighted, multi-start   (a real MDS)
  D. **direct minimisation of the SHIPPED BAYES-RISK SCORE over torsion space**, i.e. the
     honest argmin of the objective inside the class of structures the pipeline can emit.
     L-BFGS-B over (phi, psi) with batched finite differences through `build_ca_exact`,
     8 multi-starts (pool members' own torsions + random Ramachandran draws).

D is the decisive arm.  Also reports the SCORE each decoder attains, so we can see whether
a decoder that reaches a LOWER score reaches a WORSE structure -- which is the real content
of the claim.
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I


def classical_mds(D):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1][:3]
    return V[:, idx] * np.sqrt(np.maximum(w[idx], 0))[None, :]


def smacof(D, W, n_iter=400, seed=0, n_start=5):
    n = D.shape[0]
    rng = np.random.default_rng(seed)
    best, bestx = np.inf, None
    for s in range(n_start):
        X = classical_mds(D) if s == 0 else rng.normal(0, 5, (n, 3))
        Wsum = W.sum(1)
        for _ in range(n_iter):
            d = np.linalg.norm(X[:, None] - X[None, :], axis=-1)
            np.fill_diagonal(d, 1.0)
            Bm = -W * D / d
            np.fill_diagonal(Bm, 0.0)
            np.fill_diagonal(Bm, -Bm.sum(1))
            X = (Bm @ X) / np.maximum(Wsum, 1e-9)[:, None]
        d = np.linalg.norm(X[:, None] - X[None, :], axis=-1)
        st = float((W * (d - D) ** 2).sum())
        if st < best:
            best, bestx = st, X
    return bestx


def make_score(dg, n):
    """Batched shipped Bayes-risk score of (B, n, 3) coordinate sets."""
    grid = np.asarray(dg["grid"], float)
    risk = np.asarray(dg["risk"], float)
    i, j = I.pair_index(n)
    rows = np.arange(risk.shape[0])[None, :]

    def f(X):
        X = np.atleast_3d(np.asarray(X, float))
        d = np.linalg.norm(X[:, i, :] - X[:, j, :], axis=-1)
        g = np.clip(((d - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
        return risk[rows, g].mean(1)
    return f, i, j


def main(limit=None, n_start=8, maxiter=250, seed=0):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    rows = []
    t0 = time.time()
    for t in tg:
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        W = u["W"][p]; rr = u["rr"][p]; nat = u["nat_ca"]
        dg = I.distogram(pdb, t["seq"], t["fold"])
        sfun, ii, jj = make_score(dg, n)

        # dense E[d] matrix + shell weights implied by the score
        Dm = np.zeros((n, n)); Dm[ii, jj] = dg["expected"]; Dm = Dm + Dm.T
        for s in range(1, 2):                       # |i-j| < 2 has no prediction: ideal CA
            for a in range(n - s):
                Dm[a, a + s] = Dm[a + s, a] = 3.80
        Wm = np.zeros((n, n)); Wm[ii, jj] = 1.0; Wm = Wm + Wm.T
        for a in range(n - 1):
            Wm[a, a + 1] = Wm[a + 1, a] = 1.0

        r = dict(pdb=pdb, n=n, fold=t["fold"], fail=int(pdb in I.FAIL18))
        # ---- A classical MDS
        Xa = classical_mds(Dm)
        r["mds_rmsd"] = float(I.ca_rmsd(Xa, nat)); r["mds_score"] = float(sfun(Xa[None])[0])
        # ---- B MDS -> production projection
        outb = I.project(Xa, t["seq"], t["fold"])
        r["mdsproj_rmsd"] = float(I.ca_rmsd(outb["fit_ca"], nat))
        r["mdsproj_score"] = float(sfun(outb["fit_ca"][None])[0])
        # ---- C SMACOF
        Xc = smacof(Dm, Wm, seed=seed)
        r["smacof_rmsd"] = float(I.ca_rmsd(Xc, nat)); r["smacof_score"] = float(sfun(Xc[None])[0])
        outc = I.project(Xc, t["seq"], t["fold"])
        r["smacofproj_rmsd"] = float(I.ca_rmsd(outc["fit_ca"], nat))
        r["smacofproj_score"] = float(sfun(outc["fit_ca"][None])[0])

        # ---- D direct argmin of the SHIPPED SCORE over torsion space
        rng = np.random.default_rng(seed)
        sc_pool = sfun(W)
        starts = [np.concatenate([u["PHI"][p][k], u["PSI"][p][k]])
                  for k in np.argsort(sc_pool, kind="stable")[:max(1, n_start // 2)]]
        while len(starts) < n_start:
            starts.append(rng.uniform(-np.pi, np.pi, 2 * n))

        def obj(x):
            return float(sfun(I.build_ca(x[:n], x[n:])[None])[0])

        def grad(x, h=1e-4):
            B = np.repeat(x[None], 2 * len(x) + 1, 0)
            for k in range(len(x)):
                B[2 * k + 1, k] += h; B[2 * k + 2, k] -= h
            v = sfun(I.build_ca(B[:, :n], B[:, n:]))
            return (v[1::2] - v[2::2]) / (2 * h)

        bestv, bestx = np.inf, None
        for x0 in starts:
            res = minimize(obj, np.asarray(x0, float), jac=grad, method="L-BFGS-B",
                           options=dict(maxiter=maxiter))
            if res.fun < bestv:
                bestv, bestx = float(res.fun), res.x
        Xd = I.build_ca(bestx[:n], bestx[n:])
        r["tors_rmsd"] = float(I.ca_rmsd(Xd, nat)); r["tors_score"] = bestv

        # ---- references on the SAME score scale
        rec = I.shipped_record(pdb)
        r["emit_rmsd"] = float(rec["rmsd_fit"])
        r["emit_score"] = float(sfun(np.asarray(rec["fit_ca"], float)[None])[0])
        r["nat_score"] = float(sfun(nat[None])[0])
        r["nat_rmsd"] = 0.0
        k = int(np.argmin(sc_pool))
        r["argmin_rmsd"] = float(rr[k]); r["argmin_score"] = float(sc_pool[k])
        r["poolbest_rmsd"] = float(rr.min())
        rows.append(r)
        print(json.dumps({k2: (round(v, 4) if isinstance(v, float) else v)
                          for k2, v in r.items()}), flush=True)

    ks = ["mds", "mdsproj", "smacof", "smacofproj", "tors", "emit", "argmin"]
    fail = np.array([x["fail"] for x in rows], bool)
    agg = {"n": len(rows), "elapsed": time.time() - t0}
    for k in ks:
        a = np.array([x[f"{k}_rmsd"] for x in rows])
        s = np.array([x[f"{k}_score"] for x in rows])
        agg[k] = dict(rmsd=float(a.mean()), rmsd_fail=float(a[fail].mean()),
                      rmsd_other=float(a[~fail].mean()), score=float(s.mean()))
    agg["native"] = dict(rmsd=0.0, score=float(np.mean([x["nat_score"] for x in rows])))
    agg["pool_best"] = float(np.mean([x["poolbest_rmsd"] for x in rows]))
    # does a LOWER score buy a WORSE structure?  rank the decoders both ways.
    agg["ordering"] = sorted(((k, agg[k]["score"], agg[k]["rmsd"]) for k in ks), key=lambda z: z[1])
    print(json.dumps(agg, indent=1))
    I.write("adv_decode", dict(agg=agg, rows=rows))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
