"""SPRINT 14 / VQE -- the collapse trap, resolved into its actual governing variable.

The brief records: "at CVaR alpha <= 0.25 two AMBER variants became literally the same
objective."  `s14/test_vqe.py` proves an identity that produces exactly that symptom, and
`s14/vqe_cvar.py` section D measures it on the real AMBER pairs -- and the measurement does
NOT reproduce "small alpha collapses them".  On a BROAD distribution, shrinking alpha makes
two genuinely different objectives MORE distinguishable, not less.

So the governing variable is not alpha alone.  From the identity, the collapse condition is
``alpha <= p(x*)`` -- a joint condition on alpha AND on how concentrated the distribution
is.  This module varies both, on the real AMBER and Legacy pairs, and reports the surface.

    python -m s14.vqe_collapse
"""
from __future__ import annotations

import numpy as np

from core import quantum as Q
from s14 import vqe_lib as V

ALPHAS = (1.0, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01)
CONC = (2.0, 0.6, 0.1, 0.02, 0.005, 0.001)      # Dirichlet concentration: broad -> peaked


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else float("nan")


def pairs(pdb="1CS9", masked=True):
    """The AMBER pairs, on the CLEAN subset only.

    THE MASK MATTERS.  The cached AMBER subset in `qarch_enum_*.npz` is 40% ORACLE
    CONDITIONED -- its populations were not drawn uniformly and it is ~0.49 A better than
    the space it came from (1CS9: subset mean RMSD 3.573 vs full-enumeration 4.064).  Only
    ``amber_kind == 0`` is a uniform draw (1CS9: 1,193 rows, mean RMSD 4.061 against the
    space's 4.064).  `masked=False` reproduces the contaminated version so the size of the
    difference can be reported rather than asserted.  The `rmsd`, `legacy`, `prior` and
    `leg_*` columns are full-enumeration and need no mask.
    """
    e = V.Enum(pdb)
    d = np.load(f"{V.ENUM_DIR}/qarch_enum_{pdb}.npz")
    keep = (np.asarray(d["amber_kind"]) == 0) if masked \
        else np.ones(len(d["amber_kind"]), bool)
    amb = {c[4:]: np.asarray(d[c], float)[keep] for c in d.files if c.startswith("amb_")}
    tot = np.asarray(d["amber_total"], float)[keep]
    idx = np.asarray(d["amber_idx"])[keep]
    return {
        "AMBER total vs AMBER-minus-solvation":
            (V.uniformise(tot), V.uniformise(tot - amb["solvation"])),
        "AMBER total vs AMBER nonbonded+solvation":
            (V.uniformise(tot), V.uniformise(amb["nonbonded"] + amb["solvation"])),
        "AMBER total vs AMBER torsion-only":
            (V.uniformise(tot), V.uniformise(amb["torsion"])),
        "Legacy vs Legacy-minus-steric":
            (V.uniformise(e.legacy[idx]), V.uniformise((e.legacy - e.leg["steric"])[idx])),
        "Legacy vs AMBER total":
            (V.uniformise(e.legacy[idx]), V.uniformise(tot)),
    }


def contamination_check(pdb="1CS9"):
    """How much does the oracle-conditioned AMBER subset move the rank correlations?"""
    print("=" * 100)
    print("0. AMBER MASK CHECK -- the cached subset is 40% oracle-conditioned")
    print("=" * 100)
    a = pairs(pdb, masked=True)
    b = pairs(pdb, masked=False)
    print(f"{'pair':45s} {'rho (kind==0)':>14s} {'rho (all rows)':>15s} {'shift':>8s}")
    for k in a:
        ra = V.spearman(*a[k])
        rb = V.spearman(*b[k])
        print(f"{k[:43]:45s} {ra:+14.4f} {rb:+15.4f} {ra-rb:+8.4f}")
    print(f"\n  n rows: clean {len(next(iter(a.values()))[0])}, "
          f"contaminated {len(next(iter(b.values()))[0])}")
    print("  EVERY AMBER number below uses the clean mask.")


def main():
    contamination_check()
    print()
    P = pairs(masked=True)
    M = len(next(iter(P.values()))[0])
    rng = np.random.default_rng(0)
    print("=" * 100)
    print("THE COLLAPSE SURFACE: cosine of dCVaR/dp between two objectives")
    print("  rows = distribution concentration (Dirichlet alpha_D), cols = CVaR alpha")
    print("  1.000 means the two objectives are LITERALLY the same CVaR function there.")
    print(f"  subset size {M}")
    print("=" * 100)
    out = {}
    for name, (E1, E2) in P.items():
        rho = V.spearman(E1, E2)
        print(f"\n--- {name}   (rho between the objectives = {rho:+.4f}) ---")
        print(f"{'Dir a_D':>9s} {'max p':>9s} " + "".join(f"{a:>9.3g}" for a in ALPHAS))
        for cD in CONC:
            p = rng.dirichlet(np.full(M, cD))
            row = []
            for a in ALPHAS:
                d1 = Q.cvar_exact(E1, p, a)[2]
                d2 = Q.cvar_exact(E2, p, a)[2]
                row.append(_cos(d1, d2))
            out.setdefault(name, {})[str(cD)] = {"max_p": float(p.max()),
                                                 "alphas": list(ALPHAS), "cos": row}
            print(f"{cD:9.3g} {p.max():9.4f} "
                  + "".join(f"{v:9.3f}" if np.isfinite(v) else f"{'flat':>9s}"
                            for v in row))
        out[name]["rho"] = rho
    print()
    print("=" * 100)
    print("READING")
    print("=" * 100)
    print("* Down a column (fixed alpha, distribution concentrating) the cosine rises to")
    print("  1.000 and then the gradient goes FLAT entirely -- that is the identity")
    print("  alpha <= p(x*) from the derivation, and it is the real collapse.")
    print("* Along a row (fixed distribution, alpha shrinking) on a BROAD distribution the")
    print("  cosine FALLS: small alpha makes genuinely different objectives MORE")
    print("  distinguishable, not less.")
    print("* A pair whose rank correlation is already ~1.000 is the same objective at")
    print("  EVERY alpha and at every concentration.  That is a property of the pair, not")
    print("  of CVaR, and it is the likeliest explanation of the recorded observation.")
    V.write("vqe_collapse", out)
    print("\nwritten -> s14/results/vqe_collapse.json")


if __name__ == "__main__":
    main()
