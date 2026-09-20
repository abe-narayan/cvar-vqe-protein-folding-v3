"""S30 lane T -- M5: the COMBINATION CEILING, predicted from the pool's own geometry.

Any native-free operator that reweights or selects pool members emits a displacement
    u = sum_x w_x (W_x - c),   sum_x w_x = 1,
which lies in the span of the pool's DEVIATIONS about the production average c. So the ORACLE-best
combination of ANY number of such fields is capped by how much of the oracle error e = t - c lies
in the WELL-CONDITIONED part of that span -- the top-k principal directions, which a bounded
reweighting can move along with O(1) coefficients.

This predicts lane D's Gram-rank / combination-rho answer from a different instrument.
Every arm is ORACLE. It is a ceiling, not an operator.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s30.s30_T_bits import kabsch_to, pairwise_rmsd_medoid  # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s30", "results", "s30_T_combo.json")
KS = [1, 2, 3, 6, 10, 15, 21, 33]
M = 75


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
        Ws = kabsch_to(W, W[b])                              # exogenous frame

        E = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
        rng = np.random.default_rng(abs(hash(("s30T-M5", pdb))) % (2 ** 32))
        top = np.lexsort((rng.random(500), E))[:M]
        c = Ws[top].mean(0)                                  # production average, same frame

        # put the native in the SAME frame (Kabsch onto c), then project out rigid motion
        t = kabsch_to(nat[None], c)[0]
        e = (t - c).ravel()

        Dv = (Ws - c).reshape(500, -1)                       # pool deviations about c
        # the deviation span's principal directions (the well-conditioned part)
        U, S, Vt = np.linalg.svd(Dv, full_matrices=False)
        share = (S ** 2) / (S ** 2).sum()
        proj = {}
        for k in KS:
            kk = min(k, Vt.shape[0])
            P = Vt[:kk]
            proj[str(k)] = float((P @ e) @ (P @ e) / (e @ e))
        rows.append(dict(pdb=pdb, n=n, d=3 * n, rank=int(Vt.shape[0]),
                         proj=proj, share1=float(share[0]),
                         rs_dev=float((S ** 2).sum() / S[0] ** 2),
                         enorm=float(np.linalg.norm(e))))
        if (fi + 1) % 40 == 0:
            print(f"  {fi+1}/{len(files)}", flush=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(rows=rows, KS=KS, M=M), open(OUT, "w"))

    print("\n" + "=" * 78)
    print("M5  THE COMBINATION CEILING FROM THE POOL's OWN GEOMETRY (n=%d, ORACLE)" % len(rows))
    print("=" * 78)
    d = np.mean([r["d"] for r in rows])
    print(f"   mean 3n = {d:.1f}   deviation-matrix stable rank {np.mean([r['rs_dev'] for r in rows]):.3f}"
          f"   PC1 share {np.mean([r['share1'] for r in rows]):.3f}")
    print("\n     k   ||Pi_k e||^2/||e||^2   rho_ceiling = sqrt(.)   isotropic null k/d")
    for k in KS:
        v = np.array([r["proj"][str(k)] for r in rows])
        print(f"   {k:3d}      {v.mean():10.4f}            {np.sqrt(v.mean()):8.4f}"
              f"            {k/d:8.4f}")
    print("\n   'rho_ceiling' is the LARGEST cosine any convex reweighting of the pool can reach")
    print("   while staying inside the top-k well-conditioned directions. It is ORACLE and it is")
    print("   an UPPER bound on the ORACLE-optimal combination of ANY set of such fields.")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
