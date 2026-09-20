"""S30 lane T -- M1/M2/M3 of s30/PREREG_S30_T.md.

M1  the ORACLE order-statistic ladder inside the deployed BLOSUM top-500 pool, and the
    value-of-a-bit law fitted to it.
M2  what the 3,252 bits of retrieval choice bought, priced in search-equivalent bits.
M3  the rank-collapse pre-check for any tail-then-aggregate lift (Theorem T1b).

Everything from `s8/generate_univ/<pdb>.npz`, the 126 pinned universes. `rr` is ORACLE.
No deployable parameter is chosen here.
"""
import json
import sys
import glob
import os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s30", "results", "s30_T_bits.json")

RS = [0, 1, 2, 3, 4, 5, 6, 7, 8]          # N = 1..256
NS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 500]
NDRAW = 20                                 # random-500 control draws per target


def kabsch_to(P, Q):
    """Superpose each P[k] (n,3) onto Q (n,3); return superposed coords, reflection forbidden."""
    Qc = Q - Q.mean(0)
    Pc = P - P.mean(1, keepdims=True)
    C = np.einsum("kni,nj->kij", Pc, Qc)
    U, S, Vt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(np.einsum("kij,kjl->kil", U, Vt)))
    D = np.zeros_like(U)
    D[:, 0, 0] = 1.0
    D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    R = np.einsum("kij,kjl,klm->kim", U, D, Vt)
    return np.einsum("kni,kij->knj", Pc, R) + Qc.mean(0)


def pairwise_rmsd_medoid(W):
    """Index of the medoid of W (k,n,3) under CA-RMSD, exactly as s12/instrument.py:115."""
    k, n, _ = W.shape
    Wc = W - W.mean(1, keepdims=True)
    G = np.einsum("kni,lnj->klij", Wc, Wc)
    U, S, Vt = np.linalg.svd(G.reshape(k * k, 3, 3))
    det = np.sign(np.linalg.det(np.einsum("mij,mjl->mil", U, Vt)))
    Ssum = S[:, 0] + S[:, 1] + det * S[:, 2]
    sq = (Wc ** 2).sum((1, 2))
    msd = (sq[:, None] + sq[None, :] - 2 * Ssum.reshape(k, k)) / n
    msd = np.maximum(msd, 0.0)
    return int(np.argmin(np.sqrt(msd).mean(1)))


def stable_rank(A):
    """||A||_F^2 / ||A||_2^2 on the column-centred A."""
    Ac = A - A.mean(0, keepdims=True)
    s = np.linalg.svd(Ac, compute_uv=False)
    if s[0] <= 0:
        return 1.0, 1.0
    lam = s ** 2
    return float(lam.sum() / lam[0]), float(lam[0] / lam.sum())


def pair_index(n, min_sep=2):
    i, j = np.triu_indices(n, k=min_sep)
    return i, j


def main(limit=None, do_m3=True):
    rows = []
    files = UNIV if limit is None else UNIV[:limit]
    for fi, f in enumerate(files):
        z = np.load(f)
        pdb = str(z["pdb"])
        n = int(z["n"])
        rr = z["rr"].astype(np.float64)
        order = z["order"].astype(np.int64)
        nw = rr.shape[0]
        pool = order[:500]
        rr_pool = rr[pool]

        # ---- M1: prefix ladder inside the deployed pool, in the deployed order
        ladder = [float(rr_pool[:min(N, 500)].min()) for N in NS]

        # ---- M2: random-500 control, and a random ladder over the whole universe
        rng = np.random.default_rng(abs(hash(("s30T", pdb))) % (2 ** 32))
        rand500 = []
        for _ in range(NDRAW):
            sub = rng.choice(nw, size=min(500, nw), replace=False)
            rand500.append(float(rr[sub].min()))
        # random ladder: mean best-of-N over NDRAW draws, N up to nw
        rand_ladder = {}
        for N in [1, 2, 4, 8, 16, 32, 64, 128, 256, 500, 1000, 2000, 4000, 8000, 16000]:
            if N > nw:
                break
            v = []
            for _ in range(NDRAW):
                sub = rng.choice(nw, size=N, replace=False)
                v.append(float(rr[sub].min()))
            rand_ladder[N] = float(np.mean(v))

        row = dict(pdb=pdb, n=n, nw=nw,
                   ladder=ladder,
                   rr_pool_mean=float(rr_pool.mean()), rr_pool_sd=float(rr_pool.std()),
                   rr_univ_mean=float(rr.mean()), rr_univ_sd=float(rr.std()),
                   rr_univ_min=float(rr.min()),
                   rand500=float(np.mean(rand500)), rand500_sd=float(np.std(rand500)),
                   rand_ladder={str(k): v for k, v in rand_ladder.items()})

        # ---- M3: rank collapse of the candidate feature matrices
        if do_m3:
            W = z["W"][pool].astype(np.float64)
            b = pairwise_rmsd_medoid(W)
            Ws = kabsch_to(W, W[b])
            A_coord = Ws.reshape(len(pool), -1)
            rs_c, top_c = stable_rank(A_coord)
            i, j = pair_index(n, 2)
            D = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=2)
            rs_d, top_d = stable_rank(D)
            row.update(rs_coord=rs_c, top_coord=top_c, rs_dist=rs_d, top_dist=top_d,
                       npairs=int(len(i)), dcoord=int(3 * n))
        rows.append(row)
        if (fi + 1) % 20 == 0:
            print(f"  {fi+1}/{len(files)}", flush=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(dict(rows=rows, NS=NS, NDRAW=NDRAW), fh)
    print("wrote", OUT, len(rows), "targets")


if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else None
    main(lim)
