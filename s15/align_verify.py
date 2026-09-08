"""SPRINT 15, ALIGN -- independent verification of the analytic Jacobian.

`s15/align_lib.sup_jacobian` is an ANALYTIC derivative built from the rigid-suffix-rotation
identity, with the six rigid-body modes projected out algebraically.  `s15/info_null.jacobian`
is a CENTRAL-DIFFERENCE derivative at 1 degree with each perturbed trace Kabsch-superposed
onto the base trace.  They share no code.  The brief requires an independent implementation
for every headline number; this is it.

Also checks, on the same targets:
  * the four inert torsions have EXACTLY zero (or rigid-only) columns;
  * `||J e|| / sqrt(n)` predicts the true Kabsch RMSD of a small perturbation;
  * `pair_jacobian` matches central differences of the pair distances.

    python -m s15.align_verify
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import info_null as NL              # noqa: E402
from s15 import seed as SD                   # noqa: E402
import peptide_db as pdb                     # noqa: E402


def main(n_targets=5):
    tg = I.targets()[:n_targets]
    print("ANALYTIC vs CENTRAL-DIFFERENCE JACOBIAN (independent implementations)\n")
    print(f"{'pdb':<7}{'n':>4}{'cosine':>20}{'max|dJ|':>12}{'relerr':>11}"
          f"{'inert cols':>28}")
    worst = 0.0
    for t in tg:
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        phi, psi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        n = len(phi)
        Ja, _, CA = A.sup_jacobian(phi, psi)
        Jf = NL.jacobian(phi, psi)
        cos = float((Ja * Jf).sum() /
                    (np.linalg.norm(Ja) * np.linalg.norm(Jf)))
        mx = float(np.abs(Ja - Jf).max())
        rel = mx / float(np.abs(Jf).max())
        worst = max(worst, rel)
        cn = np.linalg.norm(Ja, axis=0)
        inert = [cn[0], cn[n - 1], cn[n], cn[2 * n - 1]]      # phi0, phi_{n-1}, psi0, psi_{n-1}
        print(f"{p:<7}{n:>4}{cos:>20.15f}{mx:>12.2e}{rel:>11.2e}"
              f"   {inert[0]:.1e} {inert[1]:.1e} {inert[2]:.1e} {inert[3]:.1e}")

    # ---- the linearisation itself: does ||J e||/sqrt(n) predict the real RMSD?
    print("\nFIRST-ORDER RMSD CHECK  ||J e||/sqrt(n) vs true Kabsch RMSD of the perturbation")
    print(f"{'pdb':<7}{'|e| (deg)':>11}{'first-order':>14}{'true':>10}{'ratio':>9}")
    for t in tg[:3]:
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        phi, psi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        n = len(phi)
        J, _, _ = A.sup_jacobian(phi, psi)
        base = I.build_ca(phi, psi)
        rng = SD.stable_rng(p, 11)
        for scale in (0.2, 1.0, 5.0):
            e = rng.normal(0, np.radians(scale), 2 * n)
            fo = float(np.linalg.norm(J @ e) / np.sqrt(n))
            tr = float(I.ca_rmsd(I.build_ca(phi + e[:n], psi + e[n:]), base))
            print(f"{p:<7}{scale:>11.2f}{fo:>14.5f}{tr:>10.5f}{fo / max(tr, 1e-12):>9.4f}")

    # ---- pair jacobian
    print("\nPAIR JACOBIAN vs central differences")
    for t in tg[:3]:
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        phi, psi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        n = len(phi)
        i, j = I.pair_index(n)
        Jr, CA = A.raw_jacobian(phi, psi)
        G, d0 = A.pair_jacobian(Jr, CA, i, j)
        h = 1e-6
        num = np.zeros_like(G)
        for k in range(2 * n):
            dp = np.zeros(n)
            ds = np.zeros(n)
            (dp if k < n else ds)[k % n] = h
            ca1 = np.asarray(I.build_ca(phi + dp, psi + ds), float)
            ca2 = np.asarray(I.build_ca(phi - dp, psi - ds), float)
            d1 = np.sqrt(((ca1[i] - ca1[j]) ** 2).sum(1))
            d2 = np.sqrt(((ca2[i] - ca2[j]) ** 2).sum(1))
            num[:, k] = (d1 - d2) / (2 * h)
        cos = float((G * num).sum() / (np.linalg.norm(G) * np.linalg.norm(num)))
        print(f"  {p:<7} cosine {cos:.15f}   max|dG| {np.abs(G - num).max():.2e}")

    print(f"\nworst relative Jacobian disagreement over {len(tg)} targets: {worst:.2e}")
    return worst


if __name__ == "__main__":
    main()
