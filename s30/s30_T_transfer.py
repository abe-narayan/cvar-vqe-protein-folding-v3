"""S30 lane T -- M6: price the K = 5,000 direction search as the best-of-K it is.

Lane Q's request (s30/briefs/Q_to_T_quadric_and_rank.md section 2), endorsed by the coordinator.
`best_of_k_within` needs (n_targets, K) with column k the SAME SETTING on every target, so the
direction bank must be COMMON across targets. It is made comparable by expressing every direction
in the target's own top-6 principal basis with a deterministic, native-free sign convention
(each PC's largest-|loading| entry made positive; PCs already ordered by eigenvalue).

Direction k therefore means "the same linear combination of this target's own sign-fixed top-6
pair-distance PCs" -- a transferable RULE, not a per-target fit. That is exactly what split-half
transfer should be asked about.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "s24"))
from s30.s30_T_bits import kabsch_to, pairwise_rmsd_medoid, pair_index  # noqa: E402
from s30.s30_T_quadric import rmsd_batch, bottom_m_means  # noqa: E402
from s24 import stats_lib as ST  # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s30", "results", "s30_T_transfer.json")
K = 5000
M = 75
KPC = 6


def main(limit=None):
    # ONE common direction bank, drawn once, used on every target.
    rb = np.random.default_rng(20300920)
    G = rb.standard_normal((KPC, K))
    G /= np.linalg.norm(G, axis=0, keepdims=True)
    nq = KPC * (KPC + 1) // 2
    A = rb.standard_normal((nq, K))
    Bq = rb.standard_normal((KPC, K))
    nrm = np.sqrt((A ** 2).sum(0) + (Bq ** 2).sum(0))
    A /= nrm
    Bq /= nrm

    rows, MH, MQ, pdbs = [], [], [], []
    files = UNIV if limit is None else UNIV[:limit]
    for fi, f in enumerate(files):
        z = np.load(f)
        pdb = str(z["pdb"])
        n = int(z["n"])
        pool = z["order"].astype(np.int64)[:500]
        nat = z["nat_ca"].astype(np.float64)
        W = z["W"][pool].astype(np.float64)
        b = pairwise_rmsd_medoid(W)
        Ws = kabsch_to(W, W[b])
        Wflat = Ws.reshape(500, -1)

        i, j = pair_index(n, 2)
        P = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=2)
        P = P - P.mean(0, keepdims=True)
        U, S, Vt = np.linalg.svd(P, full_matrices=False)
        k = min(KPC, S.shape[0])
        # DETERMINISTIC, NATIVE-FREE SIGN CONVENTION: largest-|loading| entry of each PC positive
        sgn = np.sign(Vt[np.arange(k), np.argmax(np.abs(Vt[:k]), axis=1)])
        sgn[sgn == 0] = 1.0
        Z = (U[:, :k] * S[:k]) * sgn
        Z = Z / (Z.std(0, keepdims=True) + 1e-12)
        if k < KPC:                                    # pad (never fires at n >= 9)
            Z = np.hstack([Z, np.zeros((500, KPC - k))])

        iu, ju = np.triu_indices(KPC)
        V = Z[:, iu] * Z[:, ju]
        V = V / (V.std(0, keepdims=True) + 1e-12)

        sh = (Z @ G).T.astype(np.float64)
        rh = rmsd_batch(bottom_m_means(sh, Wflat, M).reshape(K, n, 3), nat)
        sq = (V @ A + Z @ Bq).T.astype(np.float64)
        rq = rmsd_batch(bottom_m_means(sq, Wflat, M).reshape(K, n, 3), nat)

        MH.append(rh.astype(np.float32))
        MQ.append(rq.astype(np.float32))
        pdbs.append(pdb)
        rows.append(dict(pdb=pdb, half_best=float(rh.min()), half_mean=float(rh.mean()),
                         quad_best=float(rq.min()), quad_mean=float(rq.mean())))
        if (fi + 1) % 20 == 0:
            print(f"  {fi+1}/{len(files)}", flush=True)

    MH = np.array(MH)
    MQ = np.array(MQ)
    res = {}
    for nm, Mx in [("HALFSPACE", MH), ("QUADRIC", MQ)]:
        o = ST.best_of_k_within(Mx, n_boot=400, seed_parts=("s30T", nm))
        res[nm] = o
        print("\n" + "=" * 78)
        print(f"M6  best_of_k_within -- {nm}, COMMON direction bank, K={K}, n={len(pdbs)}")
        print("=" * 78)
        print(f"   observed per-target gain vs its own K-mean : {o['observed_gain']:+.4f} A")
        print(f"   VALID across-target null                   : {o['null_across_targets']:+.4f} A"
              f"   ({100*o['share_accounted']:.0f}% accounted)")
        print(f"   residual                                   : {o['residual']:+.4f} A")
        print(f"   SPLIT-HALF TRANSFER (the number to quote)  : {o['split_half']:+.4f} A"
              f"   ({100*o['split_half_frac']:.0f}% of the oracle)")
        print(f"   k_eff {o['k_eff']:.1f} of {K}")
        print(f"   VERDICT: {o['verdict']}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(rows=rows, K=K, M=M, KPC=KPC, bok=res,
                   note="common direction bank; PC signs fixed by largest-|loading|-positive"),
              open(OUT, "w"))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
