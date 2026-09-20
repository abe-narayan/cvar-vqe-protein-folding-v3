"""S30 lane T -- M4: the ORACLE ceiling of the quadric (second-moment) tail class.

Registered in s30/PREREG_S30_T.md addendum 1. Budget-matched: HALFSPACE, QUADRIC and FREE each
get K = 5,000 ORACLE-scored draws, so the comparison is not a best-of-K artefact
(`grid-oracles-are-order-statistics`). PREFIX is exactly one set at fixed m.

Every arm is ORACLE. This measures a CLASS CEILING, not a deployable operator.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s30.s30_T_bits import kabsch_to, pairwise_rmsd_medoid, pair_index  # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s30", "results", "s30_T_quadric.json")
K = 5000
M = 75
KPC = 6


def rmsd_batch(A, b):
    """CA-RMSD of each A[k] (n,3) to b (n,3), Kabsch, reflection corrected."""
    n = b.shape[0]
    bc = b - b.mean(0)
    Ac = A - A.mean(1, keepdims=True)
    C = np.einsum("kni,nj->kij", Ac, bc)
    U, S, Vt = np.linalg.svd(C)
    det = np.sign(np.linalg.det(np.einsum("kij,kjl->kil", U, Vt)))
    Ssum = S[:, 0] + S[:, 1] + det * S[:, 2]
    sq = (Ac ** 2).sum((1, 2)) + (bc ** 2).sum()
    return np.sqrt(np.maximum(sq - 2 * Ssum, 0.0) / n)


def bottom_m_means(scores, Wflat, m):
    """scores (K,500) -> the coordinate mean of the m lowest, as (K, 3n)."""
    idx = np.argpartition(scores, m - 1, axis=1)[:, :m]
    return Wflat[idx].mean(1)


def main(limit=None):
    rows = []
    files = UNIV if limit is None else UNIV[:limit]
    for fi, f in enumerate(files):
        z = np.load(f)
        pdb = str(z["pdb"])
        n = int(z["n"])
        pool = z["order"].astype(np.int64)[:500]
        nat = z["nat_ca"].astype(np.float64)
        W = z["W"][pool].astype(np.float64)
        b = pairwise_rmsd_medoid(W)
        Ws = kabsch_to(W, W[b])                     # EXOGENOUS frame, fixed once
        Wflat = Ws.reshape(500, -1)

        i, j = pair_index(n, 2)
        P = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=2)
        P = P - P.mean(0, keepdims=True)
        U, S, Vt = np.linalg.svd(P, full_matrices=False)
        k = min(KPC, S.shape[0])
        Z = U[:, :k] * S[:k]
        Z = Z / (Z.std(0, keepdims=True) + 1e-12)   # (500, k) standardised PC scores

        rng = np.random.default_rng(abs(hash(("s30T-M4", pdb))) % (2 ** 32))

        # ---- PREFIX: exactly one set at fixed m
        E = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
        key = rng.random(500)
        pref = np.lexsort((key, E))[:M]
        r_prefix = float(rmsd_batch(Wflat[pref].mean(0).reshape(1, n, 3), nat)[0])

        # ---- HALFSPACE: K random directions in R^k
        G = rng.standard_normal((k, K))
        G /= np.linalg.norm(G, axis=0, keepdims=True)
        sc = (Z @ G).T.astype(np.float64)                     # (K, 500)
        mh = bottom_m_means(sc, Wflat, M)
        r_half = rmsd_batch(mh.reshape(K, n, 3), nat)

        # ---- QUADRIC: K random (G, g); the quadric order is LINEAR in the Veronese map
        iu, ju = np.triu_indices(k)
        V = Z[:, iu] * Z[:, ju]                                # (500, k(k+1)/2)
        V = V / (V.std(0, keepdims=True) + 1e-12)
        A = rng.standard_normal((V.shape[1], K))
        Bq = rng.standard_normal((k, K))
        nrm = np.sqrt((A ** 2).sum(0) + (Bq ** 2).sum(0))
        A /= nrm
        Bq /= nrm
        scq = (V @ A + Z @ Bq).T.astype(np.float64)
        mq = bottom_m_means(scq, Wflat, M)
        r_quad = rmsd_batch(mq.reshape(K, n, 3), nat)

        # ---- FREE: K uniformly random 75-subsets (the matched best-of-K null)
        scr = rng.random((K, 500))
        mr = bottom_m_means(scr, Wflat, M)
        r_free = rmsd_batch(mr.reshape(K, n, 3), nat)

        rows.append(dict(pdb=pdb, n=n, k=k,
                         prefix=r_prefix,
                         half_best=float(r_half.min()), half_med=float(np.median(r_half)),
                         quad_best=float(r_quad.min()), quad_med=float(np.median(r_quad)),
                         free_best=float(r_free.min()), free_med=float(np.median(r_free)),
                         free_sd=float(r_free.std())))
        if (fi + 1) % 10 == 0:
            print(f"  {fi+1}/{len(files)}", flush=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(rows=rows, K=K, M=M, KPC=KPC), open(OUT, "w"))
    print("wrote", OUT, len(rows))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
