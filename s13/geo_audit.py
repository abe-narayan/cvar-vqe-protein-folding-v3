"""SPRINT 13 GEO -- EXPERIMENT 1: reproduce the instrument before measuring with it.

Nothing downstream is trustworthy unless four things are true, and all four are MEASURED
here rather than asserted:

A. **The state-derivative rule is right.**  Every metric in this study is built from
   ``d psi / d theta_i = (1/2) psi(theta + pi e_i)``.  Checked against central finite
   differences on the statevector itself.

B. **The exact CVaR gradient is right.**  ``grad_cvar_paramshift`` (projector shift rule
   chained with the envelope-theorem ``dCVaR/dp``) against ``grad_cvar_fd`` (central finite
   differences on the CVaR value) -- independent machinery, so it is a real check.

C. **The CVaR-gradient defect is present, measurable, and not used.**  The sampled score
   estimator is run with ``baseline="const"`` (the fix) and ``baseline="tail"`` (the
   recorded defect).  The tail form must reproduce a cosine near +0.58-0.66 at ~0.56-0.76x
   norm against the exact gradient.  If it does not reproduce, the audit is broken; if it
   does, every result in this study uses ``"const"`` and the defect is quantified.

D. **The cached energy table IS the model.**  Random rows re-derived from the live
   `FoldingHamiltonian` / `core.amber.refine` and compared.

E. **Objective validity comes before optimisation** (BRIEF section 5).  Because the whole
   register is enumerated, the ranking diagnostics are EXHAUSTIVE rather than sampled:
   Spearman rho(E, CA-RMSD) over all ``k**L`` states, the RMSD percentile of the energy
   argmin, and the energy percentile of the best representable structure (ORACLE).

    python -m s13.geo_audit
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from core.quantum import (cvar_exact, grad_cvar_paramshift, grad_cvar_fd,
                          grad_cvar_score)             # noqa: E402

CELLS = [("1A13", 5, 4), ("2BFI", 5, 4), ("1A13", 10, 2), ("2BFI", 10, 2)]
ALPHAS = (1.0, 0.25, 0.05)


def cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else float("nan")


def audit_A(circ, theta, h=1e-6):
    """State derivative: pi-shift rule vs central finite differences."""
    _, D = G.dstates(circ, theta)
    P = theta.size
    Gd = np.repeat(theta[None, :], 2 * P, axis=0)
    r = np.arange(P)
    Gd[2 * r, r] += h
    Gd[2 * r + 1, r] -= h
    S = circ.states_batch(Gd)
    Dfd = (S[0::2] - S[1::2]) / (2 * h)
    return {"max_abs_diff": float(np.abs(D - Dfd).max()),
            "cos": cos(D.ravel(), Dfd.ravel()),
            "norm_ratio": float(np.linalg.norm(D) / np.linalg.norm(Dfd))}


def audit_B(circ, theta, E, alpha):
    g = grad_cvar_paramshift(circ, theta, E, alpha)
    f = grad_cvar_fd(circ, theta, E, alpha, h=1e-5)
    return {"cos_paramshift_vs_fd": cos(g, f),
            "norm_ratio": float(np.linalg.norm(g) / max(np.linalg.norm(f), 1e-300)),
            "max_abs_diff": float(np.abs(g - f).max()),
            "grad_norm": float(np.linalg.norm(g))}


def audit_C(circ, theta, E, alpha, shots=16384, seed=0):
    g = grad_cvar_paramshift(circ, theta, E, alpha)
    out = {}
    for b in ("const", "tail"):
        gs, _ = grad_cvar_score(circ, theta, E, alpha, shots=shots,
                                rng=np.random.default_rng(seed), baseline=b)
        out[b] = {"cos": cos(gs, g),
                  "norm_ratio": float(np.linalg.norm(gs) / max(np.linalg.norm(g), 1e-300))}
    return out


def berry(circ, theta):
    psi, D = G.dstates(circ, theta)
    return float(np.abs(D @ psi).max())


def validity(pdb, L, k):
    """BRIEF section 5, exhaustive: does either objective rank the native at all?"""
    R = G.rmsd_table(pdb, L, k)
    out = {"pdb": pdb, "L": L, "k": k, "N": int(len(R)),
           "rmsd_best": float(R.min()), "rmsd_mean": float(R.mean()),
           "rmsd_median": float(np.median(R))}
    for m in G.MODELS:
        E = G.energy_table(pdb, L, k, m)
        fin = np.isfinite(E)
        i_e = int(np.argmin(np.where(fin, E, np.inf)))
        i_r = int(np.argmin(R))                                  # ORACLE
        out[m] = {
            "spearman_E_rmsd": G.spearman(E[fin], R[fin]),
            "pearson_E_rmsd": float(np.corrcoef(E[fin], R[fin])[0, 1]),
            "argmin_E_rmsd": float(R[i_e]),
            "argmin_E_rmsd_percentile": float((R < R[i_e]).mean() * 100),
            "ORACLE_best_rmsd_energy_percentile": float((E[fin] < E[i_r]).mean() * 100),
            "E_min": float(E[fin].min()), "E_median": float(np.median(E[fin])),
            "E_max": float(E[fin].max()), "E_sd": float(E[fin].std()),
            "E_sd_winsor99": float(G.winsor_hi(E[fin]).std()),
            "n_nonfinite": int((~fin).sum()),
            "frac_above_1e4": float((E[fin] > 1e4).mean()),
        }
    return out


def main():
    rows = {"what": "s13 geo experiment 1 -- instrument reproduction and objective validity",
            "audit": [], "validity": [], "table_verify": []}
    for pdb, L, k in CELLS:
        for m in G.MODELS:
            rows["table_verify"].append(
                dict(pdb=pdb, L=L, k=k, model=m, **G.verify_table(pdb, L, k, m, n=16)))
    for pdb, L, k in CELLS:
        n = int(round(np.log2(k))) * L
        circ = G.circuit(n, layers=3)
        for seed in (0, 1, 2):
            theta = G.init_theta(circ.n_params(), seed)
            a = {"pdb": pdb, "L": L, "k": k, "n_qubits": n, "seed": seed,
                 "state_deriv": audit_A(circ, theta),
                 "berry_max_abs": berry(circ, theta)}
            for m in G.MODELS:
                E = G.energy_table(pdb, L, k, m)
                Ez, mu, sd = G.zscore(G.winsor_hi(E))
                a[m] = {"z_mu": mu, "z_sd": sd}
                for al in ALPHAS:
                    a[m][f"alpha{al}"] = {"paramshift_vs_fd": audit_B(circ, theta, Ez, al),
                                          "score_baselines": audit_C(circ, theta, Ez, al)}
            rows["audit"].append(a)
        rows["validity"].append(validity(pdb, L, k))
    G.write("geo_audit", rows)

    # ---- console summary
    print("\nA/B  state-derivative and exact-gradient checks (worst over all cells):")
    print("  state pi-shift vs FD   max|d| = %.3e   cos = %.12f"
          % (max(r["state_deriv"]["max_abs_diff"] for r in rows["audit"]),
             min(r["state_deriv"]["cos"] for r in rows["audit"])))
    print("  Berry term <d_i psi|psi>  max|.| = %.3e"
          % max(r["berry_max_abs"] for r in rows["audit"]))
    cbs = [r[m][f"alpha{al}"]["paramshift_vs_fd"]["cos_paramshift_vs_fd"]
           for r in rows["audit"] for m in G.MODELS for al in ALPHAS]
    print("  paramshift vs finite differences: min cos = %.12f" % min(cbs))
    print("\nC  sampled score estimator vs the exact gradient (mean over cells):")
    for b in ("const", "tail"):
        cc = [r[m][f"alpha{al}"]["score_baselines"][b]["cos"]
              for r in rows["audit"] for m in G.MODELS for al in ALPHAS if al < 1.0]
        nn = [r[m][f"alpha{al}"]["score_baselines"][b]["norm_ratio"]
              for r in rows["audit"] for m in G.MODELS for al in ALPHAS if al < 1.0]
        print(f"  baseline={b:6s} cos {np.mean(cc):+.4f}   |g|/|g_exact| {np.mean(nn):.3f}")
    print("\nD  table == live model:  max abs diff = %.3e"
          % max(r["max_abs_diff"] for r in rows["table_verify"]))
    print("\nE  objective validity (exhaustive over the register):")
    print("  cell            model   rho(E,RMSD)  argminE RMSD  (pctile)  best-RMSD E-pctile")
    for v in rows["validity"]:
        for m in G.MODELS:
            print("  %-5s L=%-2d k=%d %-7s %+8.3f %11.3f %9.1f %14.1f"
                  % (v["pdb"], v["L"], v["k"], m, v[m]["spearman_E_rmsd"],
                     v[m]["argmin_E_rmsd"], v[m]["argmin_E_rmsd_percentile"],
                     v[m]["ORACLE_best_rmsd_energy_percentile"]))


if __name__ == "__main__":
    main()
