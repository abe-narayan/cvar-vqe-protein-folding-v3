"""ASM-9: the OVERLAPPING-PIECE join — a physically realisable alternative to torsion
concatenation, and a second measurement of the placement cost.

Piece 1 covers native residues [0, c+o), piece 2 covers [c, n); they share o residues.
Piece 2 is placed by superposing its FIRST o residues onto piece 1's LAST o residues, and
the emitted chain is piece 1's own coordinates for [0, c+o) followed by the transformed
piece 2 for [c+o, n).  Both pieces keep their REAL deposited geometry; the only freedom
consumed is the o-residue superposition, so the join needs no native and is deployable.
o >= 3 makes the join fully determined.

The oracle here restricts each side to the top-M pieces by its OWN rigid fit to its native
segment and then evaluates every pair exactly.

Usage: python -m s12.asm_overlap [--M 150] [--o 3]
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402

LMIN = 4


def kabsch_pairs(B, Aref):
    """R, t with R @ b + t ~= a, for batches B (m,o,3) -> Aref (m,o,3)."""
    bc = B.mean(1, keepdims=True); ac = Aref.mean(1, keepdims=True)
    Bc = B - bc; Ac = Aref - ac
    H = np.einsum("mki,mkj->mij", Bc, Ac)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("mji,mkj->mik", Vt, U)))
    D = np.tile(np.eye(3), (len(B), 1, 1)); D[:, 2, 2] = d
    R = np.einsum("mji,mjk,mlk->mil", Vt, D, U)     # (V D U^T)
    return R, ac[:, 0] - np.einsum("mij,mj->mi", R, bc[:, 0])


def overlap_assembly(t, nat, M=150, o=3):
    n, fold = t["n"], t["fold"]
    best = (np.inf, None, None)
    for c in range(LMIN, n - LMIN + 1):
        L1, L2 = c + o, n - c
        if L1 > n or L2 < o + 1 or L1 < o + 1:
            continue
        b1 = A.bank(fold, L1); b2 = A.bank(fold, L2)
        s1 = A.msd_many(b1["Wc"], b1["n2"], nat[:L1][None])[:, 0]
        s2 = A.msd_many(b2["Wc"], b2["n2"], nat[c:][None])[:, 0]
        k1 = np.argsort(s1)[:M]; k2 = np.argsort(s2)[:M]
        W1 = b1["W"][k1]; W2 = b2["W"][k2]
        ii, jj = np.meshgrid(np.arange(len(k1)), np.arange(len(k2)), indexing="ij")
        ii = ii.ravel(); jj = jj.ravel()
        P1 = W1[ii]                                  # (m, L1, 3)
        P2 = W2[jj]                                  # (m, L2, 3)
        R, tr = kabsch_pairs(P2[:, :o], P1[:, L1 - o:])
        P2m = np.einsum("mij,mkj->mki", R, P2) + tr[:, None, :]
        X = np.concatenate([P1, P2m[:, o:]], 1)      # (m, n, 3)
        r = A.rmsd_to(X, nat)
        k = int(np.argmin(r))
        if r[k] < best[0]:
            bond = np.linalg.norm(np.diff(X[k], axis=0), axis=1)
            best = (float(r[k]), c, (float(bond.min()), float(bond.max())))
    return best


def main():
    av = sys.argv
    M = int(av[av.index("--M") + 1]) if "--M" in av else 150
    o = int(av[av.index("--o") + 1]) if "--o" in av else 3
    tg = sorted(I.targets(), key=lambda t: (t["fold"], t["pdb"]))
    res = {}; cur = None; t0 = time.time()
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); cur = t["fold"]
        u = I.load_univ(t["pdb"]); nat = u["nat_ca"]; del u
        r, c, bond = overlap_assembly(t, nat, M=M, o=o)
        res[t["pdb"]] = {"n": t["n"], "fold": t["fold"], "k2_overlap": r, "cut": c,
                         "bond_min": bond[0] if bond else None,
                         "bond_max": bond[1] if bond else None}
        if (q + 1) % 20 == 0:
            print(f"{q+1}/{len(tg)} {time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)
    I.write(f"asm_overlap_o{o}_M{M}", res)
    v = np.array([res[t["pdb"]]["k2_overlap"] for t in tg], float)
    f18 = np.array([t["pdb"] in I.FAIL18 for t in tg])
    bm = np.array([res[t["pdb"]]["bond_min"] for t in tg], float)
    bx = np.array([res[t["pdb"]]["bond_max"] for t in tg], float)
    print(f"\n2-piece OVERLAP join (o={o}, M={M}): {np.nanmean(v):.3f}  "
          f"FAIL18 {np.nanmean(v[f18]):.3f}  other {np.nanmean(v[~f18]):.3f}")
    print(f"emitted virtual CA-CA bond range: min {np.nanmin(bm):.2f}  max {np.nanmax(bx):.2f} A")
    return res


if __name__ == "__main__":
    main()
