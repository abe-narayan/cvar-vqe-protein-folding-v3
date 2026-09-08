"""SPRINT 13 GEO -- EXPERIMENT 3: gradient statistics and the barren-plateau question.

DEFINITIONS USED (and they are the established ones, not "small gradient = plateau")
-----------------------------------------------------------------------------------
For a cost ``C(theta)`` and ``theta ~ U(-pi, pi)^P``:

    * ``Var_theta[ dC/dtheta_i ]`` -- the barren-plateau order parameter.  A plateau is
      the statement that ``Var`` decays EXPONENTIALLY in the qubit count n, i.e.
      ``Var ~ b^n`` with ``b < 1``; polynomial decay is not a plateau (McClean et al. 2018;
      Cerezo et al. 2021).  Both fits are reported with their R^2 and the winner is named.
    * ``E[ ||grad C|| ]`` -- reported too, but it is NOT the order parameter: it mixes P
      (which grows with n) into the statement.
    * anisotropy ``max_i Var_i / mean_i Var_i`` and the per-layer profile.
    * concentration: share of ``||grad||^2`` carried by the largest 10% of coordinates.

THE LOCALITY CONTROL (BRIEF: "test whether ... explained by the cost-locality and causal-
cone framework")
----------------------------------------------------------------------------------------
Cerezo et al.'s result is that for this ansatz family a GLOBAL cost has exponentially
vanishing gradients at any depth while a LOCAL cost does not until the depth grows.  The
control is run on the SAME ansatz, the SAME theta samples and the SAME batched simulation,
by contracting the same probability Jacobian against different diagonal observables:

    ``Z_0``            1-local
    ``Z_0 Z_1``        2-local
    ``|0..0><0..0|``   maximally global
    ``E(x)``           the energy expectation (alpha = 1 CVaR) -- global by construction
    ``CVaR_alpha``     the actual objective

If the CVaR objective decays like the global projector and unlike ``Z_0``, locality
explains it.  If it does not, that is a finding and it is reported as one.

SCALE.  ``CVaR(aE+b) = a CVaR(E) + b`` so every gradient scales with the energy units.
Legacy and AMBER are therefore compared on energies standardised by the register's
INTERQUARTILE RANGE (median removed, IQR divided out).  The plain sd is unusable for AMBER:
it is set by a handful of hard-sphere clash configurations at 1e10-1e18 kcal/mol which
``CVaR_alpha`` discards for every alpha < 1.  Every scale considered is recorded per cell in
``scales`` so any normalisation can be re-derived post hoc.

    python -m s13.geo_grad
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from core.quantum import cvar_exact                    # noqa: E402

LAYERS = (1, 3, 5)
ALPHAS = (1.0, 0.25, 0.05)

#: (pdb, L, k) ladders.  k=4 varies peptide length at 2 bits/residue; k=2 reaches the same
#: qubit counts at DOUBLE the peptide length; k=8 holds the peptide fixed and adds states.
LADDER_K4 = [(p, L, 4) for p in ("1A13", "1A1P", "2BFI") for L in (4, 5, 6, 7, 8)]
LADDER_K2 = [(p, L, 2) for p in ("1A13", "1A1P", "2BFI") for L in (8, 10, 12)]
LADDER_K8 = [(p, L, 8) for p in ("1A13", "1A1P", "2BFI") for L in (4, 5)]


def n_samples(n_qubits: int) -> int:
    return {8: 128, 10: 128, 12: 96, 14: 64, 15: 48, 16: 32}.get(n_qubits, 32)


def jacobian(circ, theta):
    """``(P, dim)`` exact d p(x) / d theta_i, one batched simulation."""
    PR = circ.probs_batch(circ._shift_grid(np.asarray(theta, float), np.pi / 2))
    return (PR[0::2] - PR[1::2]) / 2.0


def observables(n: int, dim: int):
    j = np.arange(dim)
    z0 = 1.0 - 2.0 * ((j >> (n - 1)) & 1)
    z1 = 1.0 - 2.0 * ((j >> (n - 2)) & 1)
    proj = np.zeros(dim); proj[0] = 1.0
    return {"Z0": z0.astype(float), "Z0Z1": (z0 * z1).astype(float), "proj0": proj}


def summarise(Gd: np.ndarray) -> dict:
    """`Gd` is (n_samples, P) gradients."""
    var_i = Gd.var(0)
    nrm = np.linalg.norm(Gd, axis=1)
    P = Gd.shape[1]
    top = max(1, int(round(0.1 * P)))
    sq = np.sort(Gd ** 2, axis=1)[:, ::-1]
    conc = (sq[:, :top].sum(1) / np.maximum(sq.sum(1), 1e-300)).mean()
    return {"var_mean": float(var_i.mean()), "var_max": float(var_i.max()),
            "var_min": float(var_i.min()),
            "anisotropy": float(var_i.max() / max(var_i.mean(), 1e-300)),
            "grad_norm_mean": float(nrm.mean()), "grad_norm_sd": float(nrm.std()),
            "grad_norm_median": float(np.median(nrm)),
            "conc_top10pct": float(conc),
            "var_by_param": [float(x) for x in var_i]}


def cell(pdb, L, k, layers, alphas=ALPHAS, seed=0):
    n = int(round(np.log2(k))) * L
    circ = G.circuit(n, layers=layers)
    P = circ.n_params()
    ns = n_samples(n)
    rng = np.random.default_rng(1000 * seed + n * 10 + layers)
    TH = rng.uniform(-np.pi, np.pi, size=(ns, P))          # SAME thetas for both models
    dim = circ.dim
    obs = observables(n, dim)

    Ez = {}
    sc = {}
    for m in G.VARIANTS:
        E = G.variant_table(pdb, L, k, m)
        s = G.scales(E)
        Ez[m] = (E - np.median(E)) / (s['iqr'] if s['iqr'] > 0 else 1.0)
        sc[m] = s

    acc = {m: {f"cvar{a}": [] for a in alphas} for m in Ez}
    for m in Ez:
        acc[m]["Emean"] = []
    ctrl = {o: [] for o in obs}
    for s in range(ns):
        J = jacobian(circ, TH[s])
        p = circ.probs(TH[s])
        for o, v in obs.items():
            ctrl[o].append(J @ v)
        for m in Ez:
            acc[m]["Emean"].append(J @ Ez[m])
            for a in alphas:
                _, _q, dp = cvar_exact(Ez[m], p, a)
                acc[m][f"cvar{a}"].append(J @ dp)

    out = {"pdb": pdb, "L": L, "k": k, "layers": layers, "n_qubits": n, "P": P,
           "n_samples": ns, "control": {o: summarise(np.array(v)) for o, v in ctrl.items()}}
    for m in Ez:
        out[m] = {"scales": sc[m]}
        for key, v in acc[m].items():
            out[m][key] = summarise(np.array(v))
    return out


def main(which="all"):
    cells = []
    if which in ("k4", "all"):
        cells += LADDER_K4
    if which in ("k2", "all"):
        cells += LADDER_K2
    if which in ("k8", "all"):
        cells += LADDER_K8
    rows = []
    for pdb, L, k in cells:
        if not all(os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, m) + ".npy"))
                   for m in G.MODELS):
            print(f"  skip {pdb} L={L} k={k} (tables not built)", flush=True)
            continue
        for lay in LAYERS:
            t0 = time.time()
            G.gate(tries=1)
            r = cell(pdb, L, k, lay)
            r["wall_s"] = time.time() - t0
            rows.append(r)
            print("  %-5s L=%-2d k=%d n=%-2d lay=%d  Var[dC] legacy %.3e  amber %.3e  "
                  "(Z0 %.3e proj0 %.3e)  %.0fs"
                  % (pdb, L, k, r["n_qubits"], lay,
                     r["legacy"]["cvar0.25"]["var_mean"],
                     r["amber"]["cvar0.25"]["var_mean"],
                     r["control"]["Z0"]["var_mean"],
                     r["control"]["proj0"]["var_mean"], r["wall_s"]), flush=True)
            G.write("geo_grad", {"what": "gradient statistics and locality control",
                                 "layers": list(LAYERS), "alphas": list(ALPHAS),
                                 "cells": rows})
    G.write("geo_grad", {"what": "gradient statistics and locality control",
                         "layers": list(LAYERS), "alphas": list(ALPHAS), "cells": rows})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
