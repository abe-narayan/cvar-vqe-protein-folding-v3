"""Independent audit of the CVaR / VQE machinery in `core.quantum`.

NOTHING in the reference half of this file imports `core`.  The simulator is rebuilt from
`np.kron` on 2x2 gate matrices, the CVaR from an explicit sort-and-accumulate, and the
gradient from central differences on that CVaR.  If `core.quantum` and this file agree,
two independent implementations agree; if the sibling's own test file agreed with
`core.quantum` that would only prove the sibling is self-consistent.

Checks, in order:
  A  the statevector matches a kron-built reference
  B  parameter-shift == exact finite differences   (target cosine 1.000000)
  C  the constant baseline is unbiased             (target cosine ~ +0.994 sampled)
  D  the tail-only baseline is the recorded defect (+0.6556 at 0.758x norm, NO noise)
  E  the deployed default is `const`, on every call site that has one
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ------------------------------------------------------------------ reference
I2 = np.eye(2)
X = np.array([[0.0, 1.0], [1.0, 0.0]])


def ry(t):
    c, s = np.cos(t / 2.0), np.sin(t / 2.0)
    return np.array([[c, -s], [s, c]])


def op_on(n, q, g):
    """`g` (2x2) acting on wire `q` of `n`, MSB-first ordering."""
    m = np.array([[1.0]])
    for w in range(n):
        m = np.kron(m, g if w == q else I2)
    return m


def cnot(n, ctrl, tgt):
    """Full 2**n CNOT built by explicit basis relabelling -- no shortcuts shared with core."""
    dim = 1 << n
    m = np.zeros((dim, dim))
    for j in range(dim):
        cb = (j >> (n - 1 - ctrl)) & 1
        k = j ^ (1 << (n - 1 - tgt)) if cb else j
        m[k, j] = 1.0
    return m


def ref_state(n, layers, theta, ring=True):
    """Naive dense-unitary simulation of layers x (RY all wires, CNOT chain, ring)."""
    dim = 1 << n
    psi = np.zeros(dim)
    psi[0] = 1.0
    th = np.asarray(theta, float).reshape(layers, n)
    for L in range(layers):
        for q in range(n):
            psi = op_on(n, q, ry(th[L, q])) @ psi
        for q in range(n - 1):
            psi = cnot(n, q, q + 1) @ psi
        if ring and n > 2:
            psi = cnot(n, n - 1, 0) @ psi
    return psi


def ref_probs(n, layers, theta, ring=True):
    p = ref_state(n, layers, theta, ring) ** 2
    return p / p.sum()


def ref_cvar(E, p, alpha):
    """CVaR of the lower alpha tail, accumulated by an explicit loop (no argsort tricks)."""
    order = sorted(range(len(E)), key=lambda i: (E[i], i))
    acc_mass = 0.0
    acc_val = 0.0
    for i in order:
        take = min(p[i], alpha - acc_mass)
        if take <= 0:
            break
        acc_val += take * E[i]
        acc_mass += take
        if acc_mass >= alpha:
            break
    return acc_val / alpha


def ref_grad_fd(n, layers, E, theta, alpha, h=1e-6, ring=True):
    """Central differences on `ref_cvar(ref_probs(...))` -- fully independent of core."""
    th = np.asarray(theta, float)
    g = np.zeros(th.size)
    for k in range(th.size):
        tp, tm = th.copy(), th.copy()
        tp[k] += h
        tm[k] -= h
        g[k] = (ref_cvar(E, ref_probs(n, layers, tp, ring), alpha)
                - ref_cvar(E, ref_probs(n, layers, tm, ring), alpha)) / (2 * h)
    return g


def cos(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d else float("nan")


# ---------------------------------------------------------------------- audit
def main():
    from core import quantum as Q

    out = {}
    rng = np.random.default_rng(20260904)
    n, layers, alpha = 7, 3, 0.15
    theta = rng.normal(0, 0.8, n * layers)
    E = rng.normal(0, 1.0, 1 << n)

    # -- A: statevector against the kron reference ------------------------
    circ = Q.StatevectorCircuit(n, layers, ring=True)
    p_core = circ.probs(theta)
    p_ref = ref_probs(n, layers, theta, ring=True)
    out["A_probs_max_abs_diff"] = float(np.max(np.abs(p_core - p_ref)))
    out["A_probs_bit_identical"] = bool(np.array_equal(p_core, p_ref))
    s_core, s_ref = circ.state(theta), ref_state(n, layers, theta)
    out["A_state_max_abs_diff"] = float(np.max(np.abs(s_core - s_ref)))

    # core's own CVaR value against the loop reference
    out["A_cvar_value_diff"] = float(abs(Q.cvar_exact(E, p_core, alpha)[0]
                                         - ref_cvar(E, p_ref, alpha)))

    # -- B: parameter shift vs exact finite differences -------------------
    g_ps = Q.grad_cvar_paramshift(circ, theta, E, alpha)
    g_fd_core = Q.grad_cvar_fd(circ, theta, E, alpha, h=1e-5)
    g_fd_ref = ref_grad_fd(n, layers, E, theta, alpha, h=1e-6)
    out["B_paramshift_vs_core_fd_cosine"] = cos(g_ps, g_fd_core)
    out["B_paramshift_vs_INDEPENDENT_fd_cosine"] = cos(g_ps, g_fd_ref)
    out["B_paramshift_vs_INDEPENDENT_fd_relnorm"] = float(
        np.linalg.norm(g_ps) / np.linalg.norm(g_fd_ref))
    out["B_paramshift_max_abs_diff_vs_ref_fd"] = float(np.max(np.abs(g_ps - g_fd_ref)))

    # -- C/D: the baselines, with ZERO sampling noise ---------------------
    # Exact-expectation form of the score estimator: enumerate the register and weight each
    # bitstring by its true probability instead of drawing shots.  Any cosine below 1 here
    # is BIAS, because there is no variance left to blame.
    dim = 1 << n
    q_quantile = Q.cvar_exact(E, p_core, alpha)[1]
    PR = circ.probs_batch(circ._shift_grid(theta, np.pi / 2))
    G = ((PR[0::2] - PR[1::2]) / 2.0 / np.maximum(p_core, 1e-15)[None, :]).T  # (dim, P)

    w = np.where(E < q_quantile, (E - q_quantile) / alpha, 0.0)

    def exact_expect(weights):
        return (p_core[:, None] * weights[:, None] * G).sum(0)

    g_none = exact_expect(w)
    g_const = exact_expect(w - (p_core * w).sum())          # constant in x
    wt = w.copy()
    m = E < q_quantile
    wt[m] -= (p_core[m] * wt[m]).sum() / p_core[m].sum()    # tail-only: the DEFECT
    wt[~m] = 0.0
    g_tail = exact_expect(wt)

    out["C_const_baseline_cosine_vs_paramshift"] = cos(g_const, g_ps)
    out["C_const_baseline_relnorm"] = float(np.linalg.norm(g_const) / np.linalg.norm(g_ps))
    out["C_no_baseline_cosine_vs_paramshift"] = cos(g_none, g_ps)
    out["D_tail_baseline_cosine_vs_paramshift"] = cos(g_tail, g_ps)
    out["D_tail_baseline_relnorm"] = float(np.linalg.norm(g_tail) / np.linalg.norm(g_ps))

    # sampled estimator, the deployable one, at high shot count
    g_s_const, _ = Q.grad_cvar_score(circ, theta, E, alpha, shots=200000,
                                     rng=np.random.default_rng(7), baseline="const")
    g_s_tail, _ = Q.grad_cvar_score(circ, theta, E, alpha, shots=200000,
                                    rng=np.random.default_rng(7), baseline="tail")
    out["C_sampled_const_cosine"] = cos(g_s_const, g_ps)
    out["D_sampled_tail_cosine"] = cos(g_s_tail, g_ps)

    # -- E: the deployed default -----------------------------------------
    import inspect
    defaults = {}
    for name in ("cvar_gradient", "grad_cvar_score", "optimise", "run", "run_cvar_vqe"):
        fn = getattr(Q, name, None)
        if fn is None or not callable(fn):
            continue
        sig = inspect.signature(fn)
        if "baseline" in sig.parameters:
            defaults[name] = sig.parameters["baseline"].default
    out["E_baseline_defaults"] = defaults
    out["E_all_defaults_const"] = all(v == "const" for v in defaults.values())

    print(json.dumps(out, indent=2, sort_keys=True))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "cvar_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    return out

if __name__ == "__main__":
    main()
