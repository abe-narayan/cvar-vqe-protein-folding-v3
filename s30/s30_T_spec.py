"""S30 lane T -- the spectral supplement to M3 and the rho<->bits dictionary.

(1) eps-effective dimension of the candidate feature matrices (how many principal directions
    a halfspace cut can actually USE), which is what Theorem T1b's cap needs.
(2) the rho <-> bits dictionary and the codebook cross-check.
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
OUT = os.path.join(ROOT, "s30", "results", "s30_T_spec.json")


def ceff(lam, frac):
    c = np.cumsum(lam) / lam.sum()
    return int(np.searchsorted(c, frac) + 1)


def main():
    rows = []
    for fi, f in enumerate(UNIV):
        z = np.load(f)
        n = int(z["n"])
        pool = z["order"].astype(np.int64)[:500]
        W = z["W"][pool].astype(np.float64)
        b = pairwise_rmsd_medoid(W)
        Ws = kabsch_to(W, W[b])
        A = Ws.reshape(len(pool), -1)
        A = A - A.mean(0, keepdims=True)
        lc = np.linalg.svd(A, compute_uv=False) ** 2
        i, j = pair_index(n, 2)
        D = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=2)
        D = D - D.mean(0, keepdims=True)
        ld = np.linalg.svd(D, compute_uv=False) ** 2
        rows.append(dict(pdb=str(z["pdb"]), n=n, npairs=int(len(i)),
                         k50c=ceff(lc, .50), k90c=ceff(lc, .90), k99c=ceff(lc, .99),
                         k50d=ceff(ld, .50), k90d=ceff(ld, .90), k99d=ceff(ld, .99),
                         lam_d=[float(x) for x in ld[:12]], lam_d_sum=float(ld.sum()),
                         lam_c=[float(x) for x in lc[:12]], lam_c_sum=float(lc.sum())))
        if (fi + 1) % 40 == 0:
            print(f"  {fi+1}/{len(UNIV)}", flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(rows=rows), open(OUT, "w"))

    ns = np.array([r["n"] for r in rows])
    print("\n" + "=" * 78)
    print("M3b  eps-EFFECTIVE DIMENSION OF THE CANDIDATE FEATURE MATRICES (n=126 pools of 500)")
    print("=" * 78)
    print(f"   mean chain length n = {ns.mean():.2f}   mean 3n-6 = {3*ns.mean()-6:.1f}"
          f"   mean npairs = {np.mean([r['npairs'] for r in rows]):.1f}")
    for tag, ks in [("A_coord", ["k50c", "k90c", "k99c"]), ("A_dist ", ["k50d", "k90d", "k99d"])]:
        v = [np.array([r[k] for r in rows]) for k in ks]
        print(f"   {tag}: k50 {v[0].mean():5.2f}   k90 {v[1].mean():5.2f}   k99 {v[2].mean():5.2f}"
              f"   (k90 max {v[1].max()}, k99 max {v[2].max()})")
    # variance share of the leading distance-map directions
    P = np.array([np.array(r["lam_d"][:6]) / r["lam_d_sum"] for r in rows])
    print("   A_dist variance share of PC1..PC6: " + "  ".join(f"{x:.3f}" for x in P.mean(0)))

    # ---- rho <-> bits dictionary
    print("\n" + "=" * 78)
    print("THE rho <-> BITS DICTIONARY   I = -(d/2) log2(1 - rho^2),  d = 3n-6")
    print("=" * 78)
    d = 3 * ns.mean() - 6
    print(f"   d = {d:.2f} (mean over the 126)")
    for lbl, rho in [("every field ever built (|rho|<=0.04)", 0.04),
                     ("B2's ceiling", 0.140),
                     ("3.00 A", 0.358), ("2.50 A (charter)", 0.628),
                     ("2.31 A", 0.695), ("1.71 A", 0.847)]:
        I = -d / 2 * np.log2(1 - rho ** 2)
        print(f"   rho {rho:5.3f}  ({lbl:36s}) -> I = {I:8.3f} bits")

    print("\n   CODEBOOK CROSS-CHECK (the point of the accounting)")
    bits = json.load(open(os.path.join(ROOT, "s30", "results", "s30_T_bits.json")))
    L = np.array([r["ladder"] for r in bits["rows"]])
    D0, D7, D9 = L[:, 0].mean(), L[:, 7].mean(), L[:, 9].mean()
    for R, Dv, nm in [(7.0, D7, "best-of-128 (the readout's own 7 bits)"),
                      (8.966, D9, "best-of-500 (the whole pool)")]:
        rho = np.sqrt(max(0.0, 1 - (Dv / D0) ** 2))
        I = -d / 2 * np.log2(1 - rho ** 2)
        print(f"   R = {R:5.3f} index bits: {D0:.4f} -> {Dv:.4f} A  =>  rho = {rho:.4f}"
              f"  =>  I = {I:7.2f} displacement bits   (ratio {I/R:5.2f}x)")
    print("   The ratio is the CODEBOOK: the index names a structure the pool already supplies.")


if __name__ == "__main__":
    main()
