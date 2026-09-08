"""SPRINT 14 / VQE -- the live-qubit audit and the signal-family calibration.

Two prerequisites for everything else, both measured rather than assumed.

1. LIVE QUBITS.  The brief records that `core/project.py` leaves `phi[0]`, `psi[n-1]` and
   `phi[n-1]` inert for the CA trace, and concludes "log2(k) qubits per chain are dead".
   That is measured here by brute force -- flip each qubit on 4,000 random configurations
   and ask whether the tabulated CA-RMSD ever changes.

2. THE SIGNAL FAMILY.  `blend_objective` mixes a realistic objective toward the true RMSD
   in rank space.  The map from the nominal mixing weight to the REALISED rank correlation
   is strongly non-linear, so it is measured on every target at every setting and the
   realised rho is what every downstream table reports.

    python -m s14.vqe_calib
"""
from __future__ import annotations

import numpy as np

from s14 import vqe_lib as V

SIGNALS = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0)
BASES = ("legacy", "prior")


def main():
    E = V.load_all()
    print("=" * 78)
    print("1. LIVE-QUBIT AUDIT (measured: flip qubit q, does the tabulated CA-RMSD move?)")
    print("=" * 78)
    print(f"{'pdb':6s} {'n':>3s} {'nominal':>8s} {'live':>5s} {'dead':>5s}  dead qubits")
    rows = {}
    for p, e in E.items():
        a = e.live_qubits()
        rows[p] = a
        print(f"{p:6s} {e.n:3d} {a['n_qubits_nominal']:8d} {a['n_live']:5d} "
              f"{len(a['dead']):5d}  {a['dead']}")
    nl = [rows[p]["n_live"] for p in E]
    print(f"\nmean live qubits {np.mean(nl):.2f} of {E['1CS9'].n_qubits} nominal "
          f"({np.mean(nl) / E['1CS9'].n_qubits:.1%})")
    print("dead-qubit blocks are residue 0 and residue n-1 on every target: "
          f"{all(rows[p]['dead'] == [0, 1, 2 * E[p].n - 2, 2 * E[p].n - 1] for p in E)}")

    print()
    print("=" * 78)
    print("2. SIGNAL-FAMILY CALIBRATION -- realised rho(objective, true RMSD)")
    print("   ORACLE DIAGNOSTIC (the family reads the native by construction)")
    print("=" * 78)
    out = {"live_qubits": rows, "signals": {}}
    for base in BASES:
        print(f"\n--- base = {base} ---")
        hdr = "  ".join(f"{s:>6.2f}" for s in SIGNALS)
        print(f"{'pdb':6s} {hdr}")
        for p, e in E.items():
            b = getattr(e, base)
            vals = []
            for s in SIGNALS:
                o = V.blend_objective(b, e.rmsd, s)
                vals.append(V.spearman(o, e.rmsd))
                out["signals"].setdefault(base, {}).setdefault(p, {})[str(s)] = {
                    "rho": vals[-1],
                    "rmsd_at_argmin": float(e.rmsd[int(np.argmin(o))]),
                    "argmin": int(np.argmin(o)),
                }
            print(f"{p:6s} " + "  ".join(f"{v:+6.3f}" for v in vals))
        print(f"{'MEAN':6s} " + "  ".join(
            f"{np.mean([out['signals'][base][p][str(s)]['rho'] for p in E]):+6.3f}"
            for s in SIGNALS))
        print(f"{'RMSD@argmin':6s}")
        print(f"{'  mean':6s} " + "  ".join(
            f"{np.mean([out['signals'][base][p][str(s)]['rmsd_at_argmin'] for p in E]):6.3f}"
            for s in SIGNALS))

    print()
    print("reference points (mean over the nine targets):")
    print(f"  best in space        {np.mean([e.rmsd.min() for e in E.values()]):.3f}")
    print(f"  uniform random draw  {np.mean([e.rmsd.mean() for e in E.values()]):.3f}")
    print(f"  worst in space       {np.mean([e.rmsd.max() for e in E.values()]):.3f}")
    V.write("vqe_calib", out)


if __name__ == "__main__":
    main()
