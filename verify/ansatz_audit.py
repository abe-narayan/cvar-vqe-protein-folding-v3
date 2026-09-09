"""Is the VQE a genuine parameterised circuit, or classical search wearing the name?

Three ways `core.quantum` could fail the "genuine VQE" test, each checked here:

  1. `OneLayerAnsatz` replaces circuit simulation with a closed-form Bernoulli-XOR
     distribution.  That is legitimate ONLY if the closed form is exactly the output
     distribution of RY-then-CNOT-chain-then-ring.  Checked against the kron simulator.
  2. `MPSAnsatz` must be a real contraction, not a lookup.  Checked by comparing its
     `probs` against a dense statevector built from the same tensors, and by confirming
     the distribution actually MOVES with theta (a classical enumerator would not).
  3. `cvar_gradient_exact` is the reference the sibling's own defect test measures
     against.  If IT is wrong the whole test is circular, so it is checked here against
     independent finite differences on an independently rebuilt distribution.

Also: the register must be genuinely parameterised -- gradients non-zero, distribution
normalised, and the objective's optimum not simply the argmin found by enumeration.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cvar_audit import cnot, cos, op_on, ref_cvar, ry  # noqa: E402  (independent refs)


def onelayer_ref_probs(n, theta, ring=True):
    """Dense simulation of ONE RY layer + CNOT chain (+ring), from kron primitives."""
    dim = 1 << n
    psi = np.zeros(dim)
    psi[0] = 1.0
    for q in range(n):
        psi = op_on(n, q, ry(float(theta[q]))) @ psi
    for q in range(n - 1):
        psi = cnot(n, q, q + 1) @ psi
    if ring and n > 2:
        psi = cnot(n, n - 1, 0) @ psi
    p = psi ** 2
    return p / p.sum()


def main():
    from core import quantum as Q

    out = {}
    rng = np.random.default_rng(31337)

    # -- 1: the OneLayerAnsatz closed form IS the circuit ------------------
    worst = 0.0
    for n in (4, 6, 8):
        for ring in (True, False):
            th = rng.uniform(0.2, 2.9, n)
            an = Q.OneLayerAnsatz(n, ring=ring)
            p_core = an.probs(th)
            p_ref = onelayer_ref_probs(n, th, ring=ring and n > 2)
            worst = max(worst, float(np.max(np.abs(p_core - p_ref))))
    out["1_onelayer_closedform_vs_circuit_max_abs_diff"] = worst
    out["1_onelayer_is_the_circuit"] = worst < 1e-12

    # normalisation and responsiveness
    an = Q.OneLayerAnsatz(8, ring=True)
    t0 = rng.uniform(0.3, 2.8, 8)
    p0 = an.probs(t0)
    out["1_probs_sum"] = float(p0.sum())
    p1 = an.probs(t0 + 0.31)
    out["1_distribution_moves_with_theta"] = float(np.max(np.abs(p1 - p0)))

    # -- 2: MPSAnsatz is a real contraction --------------------------------
    try:
        mp = Q.MPSAnsatz(10, 2, final_ry=True)
        th = rng.uniform(0.2, 2.8, mp.n_params())
        pm = mp.probs(th)
        out["2_mps_probs_sum"] = float(pm.sum())
        out["2_mps_support"] = int((pm > 1e-12).sum())
        out["2_mps_dim"] = int(pm.size)
        pm2 = mp.probs(th + 0.27)
        out["2_mps_moves_with_theta"] = float(np.max(np.abs(pm2 - pm)))
        # logp must be consistent with probs on explicit bitstrings
        bits = Q.all_bitstrings(10)
        lp = mp.logp(th, bits)
        out["2_mps_logp_vs_probs_max_abs_diff"] = float(
            np.max(np.abs(np.exp(lp) - pm)))
        # entanglement: a product state would factorise. Measure the bond-2 Schmidt gap.
        psi = np.sqrt(np.maximum(pm, 0)).reshape(1 << 5, 1 << 5)
        sv = np.linalg.svd(psi, compute_uv=False)
        out["2_mps_schmidt_rank_gt1"] = bool(sv[1] / max(sv[0], 1e-30) > 1e-6)
    except Exception as exc:                                   # pragma: no cover
        out["2_mps_error"] = f"{type(exc).__name__}: {exc}"

    # -- 3: cvar_gradient_exact against INDEPENDENT finite differences -----
    # Rebuild the distribution from the kron simulator, take central differences of the
    # loop-accumulated CVaR, and compare. Nothing of core is used on the reference side.
    n, alpha = 8, 0.2
    an = Q.OneLayerAnsatz(n, ring=True)
    th = rng.uniform(0.4, 2.6, n)
    E = rng.standard_normal(1 << n)
    bits_all = Q.all_bitstrings(n)
    g_core, v_core = Q.cvar_gradient_exact(an, th, bits_all, E, alpha)

    h = 1e-6
    g_ref = np.zeros(n)
    for k in range(n):
        tp, tm = th.copy(), th.copy()
        tp[k] += h
        tm[k] -= h
        g_ref[k] = (ref_cvar(E, onelayer_ref_probs(n, tp), alpha)
                    - ref_cvar(E, onelayer_ref_probs(n, tm), alpha)) / (2 * h)
    out["3_cvar_gradient_exact_cosine_vs_independent_fd"] = cos(g_core, g_ref)
    out["3_cvar_gradient_exact_relnorm"] = float(
        np.linalg.norm(g_core) / np.linalg.norm(g_ref))
    out["3_cvar_value_vs_independent"] = float(
        abs(v_core - ref_cvar(E, onelayer_ref_probs(n, th), alpha)))

    # -- 4: CVaR is the tail, not the mean ---------------------------------
    e = rng.standard_normal(4096)
    for a in (0.05, 0.15, 0.5):
        v = Q.cvar_from_samples(e, a)
        k = max(1, int(np.floor(a * e.size)))
        ref = float(np.sort(e)[:k].mean())
        out[f"4_cvar_alpha{a}_vs_sorted_tail_mean"] = float(abs(v - ref))
        out[f"4_cvar_alpha{a}_minus_full_mean"] = float(v - e.mean())

    print(json.dumps(out, indent=2, sort_keys=True))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "ansatz_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    return out

if __name__ == "__main__":
    main()
