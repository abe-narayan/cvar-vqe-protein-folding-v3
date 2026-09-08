"""SPRINT 13 GEO -- EXPERIMENT 5: does the energy model's PAULI-WEIGHT SPECTRUM predict
its gradient scaling?

WHY THIS IS THE RIGHT OBJECT
============================
Both energy models are DIAGONAL in the computational basis, so each is exactly a weighted
sum of Pauli-Z strings -- a Walsh-Hadamard expansion, computable exactly by one FWHT over
the enumerated register:

    E(x) = sum_S c_S chi_S(x),      chi_S(x) = (-1)^{|S AND x|},      <Z_S> = sum_x p(x) chi_S(x)

with ``sum_{S != 0} c_S^2 = Var(E)`` exactly (Parseval -- asserted, not assumed).  The
**weight spectrum** ``W(w) = sum_{|S| = w} c_S^2 / Var(E)`` is precisely the "locality" the
BRIEF's section 4 describes in words, expressed in the language a VQA audience uses.

THE PREDICTION, AND WHY IT IS A REAL TEST RATHER THAN A RESTATEMENT
===================================================================
Cerezo et al. (2021) turn a cost function's locality into a gradient-scaling statement.
Here that becomes an *empirically closed loop* with no free parameters:

  1. **Kernel.**  For a fixed (n, layers) and theta ~ U(-pi,pi), measure
     ``v(w) = mean over |S| = w and over i of Var_theta[ d<Z_S>/dtheta_i ]``.
     This is a property of the ANSATZ alone -- no energy model in it.
  2. **Spectrum.**  FWHT the energy table -> ``c_S``, hence ``W(w)``.
     This is a property of the ENERGY MODEL alone -- no circuit in it.
  3. **Prediction.**  ``Var_pred[dC/dtheta_i] = sum_S c_S^2 v(|S|)``, i.e. spectrum x kernel,
     valid if the cross-covariances ``Cov(d<Z_S>, d<Z_S'>)`` are negligible.
  4. **Measurement.**  ``Var_meas[d<E>/dtheta_i]`` from the same theta samples.
  5. The ratio ``Var_meas / Var_pred`` is the test.  It also isolates *which* assumption
     fails if it fails: ``Var_exact = sum_S c_S^2 Var_theta[d<Z_S>/dtheta_i]`` (no
     weight-averaging, still no cross terms) is computed too, so "weight is not the right
     summary" and "cross terms matter" are separated.

The whole array ``d<Z_S>/dtheta_i`` for EVERY S is obtained for free: it is the
Walsh-Hadamard transform of the probability Jacobian's rows.  No sampling over Pauli
strings is needed at all -- the spectrum, the kernel and the prediction are exact over the
complete Pauli basis.

Note the scope: steps 3-4 are exact for the ``alpha = 1`` cost ``<E>``, which is a genuine
linear observable.  ``CVaR_alpha`` for ``alpha < 1`` is a piecewise-linear functional of
``p`` and is NOT a fixed observable, so its deviation from the prediction is reported as a
measurement rather than folded into the theory.

    python -m s13.geo_pauli
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from core.quantum import cvar_exact                    # noqa: E402

LAYERS = (1, 2, 3, 5)
CELLS = ([(p, L, 4) for p in ("1A13", "1A1P", "2BFI") for L in (4, 5, 6, 7, 8)]
         + [(p, L, 2) for p in ("1A13", "1A1P", "2BFI") for L in (8, 10, 12)]
         + [(p, L, 8) for p in ("1A13", "1A1P", "2BFI") for L in (4, 5)])


def fwht(a: np.ndarray) -> np.ndarray:
    """Sylvester-ordered fast Walsh-Hadamard transform along the last axis (unnormalised).

    ``fwht(f)[S] = sum_x f[x] (-1)^{popcount(S AND x)}``.
    """
    a = np.array(a, dtype=float, copy=True)
    N = a.shape[-1]
    h = 1
    while h < N:
        b = a.reshape(*a.shape[:-1], N // (2 * h), 2, h)
        x = b[..., 0, :].copy()
        y = b[..., 1, :].copy()
        b[..., 0, :] = x + y
        b[..., 1, :] = x - y
        h *= 2
    return a


def popcount(n: int) -> np.ndarray:
    j = np.arange(1 << n, dtype=np.int64)
    c = np.zeros(1 << n, dtype=np.int16)
    for b in range(n):
        c += ((j >> b) & 1).astype(np.int16)
    return c


def spectrum(E: np.ndarray, n: int, w: np.ndarray):
    """Walsh coefficients and the Pauli-weight variance spectrum of a diagonal energy."""
    N = 1 << n
    c = fwht(np.asarray(E, float)) / N
    c2 = c ** 2
    var = float(np.var(E))
    tot = float(c2[1:].sum())
    Ww = np.bincount(w, weights=c2, minlength=n + 1)
    Ww[0] = 0.0
    share = Ww / max(tot, 1e-300)
    cum = np.cumsum(share)
    # Delta-spike reference: a constant plus a spike at ONE configuration has weight
    # spectrum exactly Binomial(n, 1/2) and mean weight exactly n/2, whatever the physics.
    # `l1_from_binomial` near 0 with `spike.top10` near 1 means the spectrum is measuring
    # a handful of steric clashes, not the force field's interaction structure.
    from math import comb
    binom = np.array([comb(n, ww) for ww in range(n + 1)], float)
    binom[0] = 0.0
    binom = binom / binom.sum()
    return {
        "delta_reference_mean_weight": n / 2.0,
        "l1_from_binomial": float(np.abs(share - binom).sum()),
        "spike": G.spike_share(E),
        "parseval_var": var, "parseval_sum_c2": tot,
        "parseval_rel_err": abs(var - tot) / max(abs(var), 1e-300),
        "weight_share": [float(x) for x in share],
        "mean_weight": float((np.arange(n + 1) * share).sum()),
        "weight_50pct": int(np.searchsorted(cum, 0.50) ),
        "weight_90pct": int(np.searchsorted(cum, 0.90)),
        "weight_95pct": int(np.searchsorted(cum, 0.95)),
        "weight_99pct": int(np.searchsorted(cum, 0.99)),
        "share_w1": float(share[1]), "share_w2": float(share[2]),
        "share_w_le2": float(share[1] + share[2]),
        "share_w_ge5": float(share[5:].sum()) if n >= 5 else 0.0,
    }, c2


def cell(pdb, L, k, layers, n_theta=None, alphas=(1.0, 0.25, 0.05), seed=0):
    n = int(round(np.log2(k))) * L
    N = 1 << n
    circ = G.circuit(n, layers=layers)
    P = circ.n_params()
    if n_theta is None:
        n_theta = {8: 128, 10: 96, 12: 64, 14: 40, 15: 32, 16: 24}.get(n, 24)
    w = popcount(n)

    Es, spec, c2 = {}, {}, {}
    for v in G.VARIANTS:
        E = G.variant_table(pdb, L, k, v)
        Es[v] = E
        spec[v], c2[v] = spectrum(E, n, w)

    rng = np.random.default_rng(7717 * n + 31 * layers + seed)
    TH = rng.uniform(-np.pi, np.pi, size=(n_theta, P))

    s1 = np.zeros((P, N)); s2 = np.zeros((P, N))
    gm = {v: {f"a{a}": np.zeros((n_theta, P)) for a in alphas} for v in Es}
    for t in range(n_theta):
        PR = circ.probs_batch(circ._shift_grid(TH[t], np.pi / 2))
        J = (PR[0::2] - PR[1::2]) / 2.0                     # (P, N) dp(x)/dtheta_i
        Jh = fwht(J)                                        # (P, N) d<Z_S>/dtheta_i
        s1 += Jh; s2 += Jh ** 2
        p = circ.probs(TH[t])
        for v, E in Es.items():
            for a in alphas:
                if a == 1.0:
                    gm[v]["a1.0"][t] = J @ E
                else:
                    _, _q, dp = cvar_exact(E, p, a)
                    gm[v][f"a{a}"][t] = J @ dp
    VarS = s2 / n_theta - (s1 / n_theta) ** 2               # (P, N) Var_theta d<Z_S>/dth_i

    # ---- ansatz kernel v(w): mean over parameters and over strings of that weight
    vw = np.zeros(n + 1)
    m = VarS.mean(0)
    for ww in range(n + 1):
        sel = (w == ww)
        vw[ww] = float(m[sel].mean()) if sel.any() else 0.0

    out = {"pdb": pdb, "L": L, "k": k, "n_qubits": n, "layers": layers, "P": P,
           "n_theta": n_theta,
           "kernel_v_of_weight": [float(x) for x in vw],
           "kernel_decay": G.loglog_fit(np.arange(1, n + 1), vw[1:])}
    for v in Es:
        meas = {f"a{a}": float(gm[v][f"a{a}"].var(0).mean()) for a in alphas}
        pred_w = float((c2[v] * vw[w]).sum())               # spectrum x kernel
        pred_x = float((c2[v][None, :] * VarS).sum(1).mean())   # exact per-string, no cross
        out[v] = {
            "spectrum": spec[v],
            "var_grad_measured": meas,
            "var_grad_pred_weight_kernel": pred_w,
            "var_grad_pred_exact_per_string": pred_x,
            "ratio_meas_over_pred_weight": meas["a1.0"] / max(pred_w, 1e-300),
            "ratio_meas_over_pred_exact": meas["a1.0"] / max(pred_x, 1e-300),
            "scales": G.scales(Es[v]),
        }
    return out


def main():
    rows = []
    done = set()
    p = os.path.join(G.RESULTS, "geo_pauli.json")
    if os.path.exists(p):                       # resume: energy tables land incrementally
        import json
        rows = json.load(open(p)).get("cells", [])
        done = {(r["pdb"], r["L"], r["k"], r["layers"]) for r in rows}
        print(f"  resuming with {len(rows)} cells already computed", flush=True)
    for pdb, L, k in CELLS:
        if not os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")):
            continue
        for lay in LAYERS:
            if (pdb, L, k, lay) in done:
                continue
            t0 = time.time()
            r = cell(pdb, L, k, lay)
            r["wall_s"] = time.time() - t0
            rows.append(r)
            print("  %-5s L=%-2d k=%d n=%-2d lay=%d | mean Pauli weight  leg %.2f  amb %.2f "
                  "| w95 leg %d amb %d | meas/pred(exact) leg %.2f amb %.2f | %.0fs"
                  % (pdb, L, k, r["n_qubits"], lay,
                     r["legacy"]["spectrum"]["mean_weight"],
                     r["amber"]["spectrum"]["mean_weight"],
                     r["legacy"]["spectrum"]["weight_95pct"],
                     r["amber"]["spectrum"]["weight_95pct"],
                     r["legacy"]["ratio_meas_over_pred_exact"],
                     r["amber"]["ratio_meas_over_pred_exact"], r["wall_s"]), flush=True)
            G.write("geo_pauli", {"what": "Pauli-weight (Walsh-Hadamard) spectra of the "
                                         "energy models and the gradient-scaling "
                                         "prediction they imply", "cells": rows})
    G.write("geo_pauli", {"what": "Pauli-weight (Walsh-Hadamard) spectra of the energy "
                                 "models and the gradient-scaling prediction they imply",
                          "cells": rows})


if __name__ == "__main__":
    main()
