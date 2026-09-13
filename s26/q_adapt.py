"""s26/q_adapt.py -- qubit-ADAPT-VQE on the DEPLOYED CVaR free energy.  Lane Q, Sprint 26.

Pre-registered in s26/PREREG_A1.md (endpoint) and s26/PREREG_A3.md (Hamiltonian variants).
Property half (operator sequences, free energies, divergences) reads no native; the endpoint
half (`--label`) reads the native and is GATED on "PHASE 0 SIGNED OFF" in s26/LEDGER.md.

WHAT IS HERE
============
1. A small exact statevector simulator for circuits made of PAULI-STRING ROTATIONS
   exp(-i theta P / 2), written on numpy from the definition of the Pauli matrices and
   NOT from core.quantum.  qubit q <-> index bit (n-1-q), so qubit 0 is the most
   significant bit, exactly the convention of core.quantum.StatevectorCircuit and of
   pennylane's qml.probs.  The deployed RY/CNOT-chain-plus-ring ansatz is ALSO built here
   (RY(theta) = exp(-i theta Y/2) is a Pauli rotation; the entangler is a basis permutation
   composed independently) so it can be checked against core.quantum to 1e-12
   (s26/q_tests.py) -- the harness itself calls core.quantum.run_cvar_vqe for the fixed arm.

2. The three operator pools, each of Pauli strings with an ODD number of Y's so that
   -iP is a real signed permutation and the amplitudes stay real (the deployed ansatz is
   real; Tang et al. PRX Quantum 2, 020310 (2021) use the same restriction):
      "V"      Tang et al.'s recursive minimal complete pool, 2n-2 strings:
               {V}_2 = {Z2 Y1, Y2};  {V}_n = { Z_n {V}_{n-1}, Y_n, Y_{n-1} }.
      "G"      Tang et al.'s local minimal complete pool, 2n-2 strings:
               {Z_{k+1} Y_k : k = 1..n-1} u {Y_k : k = 2..n}.
      "L2"     every 1- and 2-local Pauli string on ANY pair of qubits with an odd number
               of Y's: Y_q (n) and Y X, X Y, Y Z, Z Y on every pair (2 n (n-1)); 2n^2 - n.
      "L2C"    the COMPLEX option: every 1- and 2-local Pauli string (3n + 9 n(n-1)/2),
               amplitudes complex; reported only when asked for.
   Tang's qubit k (1-based) is wire k-1 here.

3. The ADAPT loop.  Objective: F = CVaR_alpha(E; p_theta) - T H(p_theta), the deployed
   objective, through core.quantum.cvar_exact and the same entropy clip (1e-15).  Growth:
   the derivative of p(x) with respect to a NEW rotation exp(-i phi P/2) appended at the
   end of the circuit, at phi = 0, is exactly (p(+pi/2) - p(-pi/2))/2 = Re(conj(psi) (-iP psi))
   elementwise, so dF/dphi|0 = sum_x d_x psi_x (M psi)_x with d = dCVaR/dp - T dH/dp and
   M = -iP.  The operator with the largest |dF/dphi| is appended (ties: lowest pool index,
   counted), every angle is re-optimised with Adam (lr 0.15, the same arithmetic as
   core.quantum.run_cvar_vqe, fresh optimiser state) for `adam_steps` steps, and growth
   stops at ||g_pool|| < eps or at `max_params` parameters.  The start is ONE layer of RY
   with theta ~ N(0, 0.6^2) (the deployed initialisation law) because from |0...0> every
   first-order ADAPT gradient is exactly zero: p_x(phi) = cos^2(phi/2) delta_{x0}
   + sin^2(phi/2) |<x|P|0>|^2 has zero derivative at phi = 0.

4. The A1/A3 harness (`--build`, `--label`, `--stats`).  `--build` reproduces the
   production pool through core.pipeline (retrieve, score, filter_pool, average, project),
   runs every arm and stores COORDINATES per target (atomic JSON, resumable).  It reads
   no native.  `--label` opens the native through core.pipeline.label and is gated.
   `--stats` runs s24.stats_lib.compare on the pre-registered contrasts.

Every array the harness emits for the fixed arm is asserted bit-for-bit against
bench_results/cache/464a0ddb5f283e04/<pdb>.json (the 126-target Config(quantum=True,
legacy=True) run): `q_ca` and `ca` must match to the last bit (`--check-cache`).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
QCACHE = os.path.join(ROOT, "bench_results", "cache", "464a0ddb5f283e04")

#: the deployed learning rate and initialisation law (core.quantum.run_cvar_vqe)
LR = 0.15
INIT_SD = 0.6
ADAM_B1, ADAM_B2, ADAM_EPS = 0.9, 0.999, 1e-8


# =============================================================== Pauli strings
def popcount(a):
    a = np.asarray(a, np.uint64)
    return np.bitwise_count(a).astype(np.int64)


class Pauli:
    """An n-qubit Pauli string as (x, z) bit masks over the basis-index bits.

    qubit q is index bit (n-1-q).  ``x`` marks X or Y on a qubit, ``z`` marks Z or Y.
    The overall phase is not tracked: rotations use -iP, whose phase is fixed by n_y.
    """

    __slots__ = ("n", "x", "z")

    def __init__(self, n, x, z):
        self.n, self.x, self.z = int(n), int(x), int(z)

    @classmethod
    def from_word(cls, word):
        """``word[q]`` is the letter on qubit q, e.g. "IYZ" = Y on qubit 1, Z on qubit 2."""
        n = len(word)
        x = z = 0
        for q, ch in enumerate(word.upper()):
            b = 1 << (n - 1 - q)
            if ch in "XY":
                x |= b
            if ch in "ZY":
                z |= b
            if ch not in "IXYZ":
                raise ValueError(word)
        return cls(n, x, z)

    def word(self):
        out = []
        for q in range(self.n):
            b = 1 << (self.n - 1 - q)
            xb, zb = bool(self.x & b), bool(self.z & b)
            out.append("Y" if (xb and zb) else "X" if xb else "Z" if zb else "I")
        return "".join(out)

    def n_y(self):
        return int(bin(self.x & self.z).count("1"))

    def weight(self):
        return int(bin(self.x | self.z).count("1"))

    def key(self):
        return (self.x << self.n) | self.z

    def commutes(self, other):
        s = bin(self.x & other.z).count("1") + bin(self.z & other.x).count("1")
        return s % 2 == 0

    def times(self, other):
        """The product string, phase dropped."""
        return Pauli(self.n, self.x ^ other.x, self.z ^ other.z)

    def __repr__(self):
        return f"Pauli({self.word()})"

    def __eq__(self, other):
        return isinstance(other, Pauli) and (self.n, self.x, self.z) == (other.n, other.x, other.z)

    def __hash__(self):
        return hash((self.n, self.x, self.z))


def pauli_action(P: Pauli):
    """``(perm, coef)`` with ``(M psi)[k] = coef[k] * psi[perm[k]]`` for ``M = -iP``.

    From the definition of the Pauli matrices on a basis state |j>:  X|b> = |1-b>,
    Z|b> = (-1)^b |b>, Y|b> = i (-1)^b |1-b>.  So P|j> = i^{n_y} (-1)^{popcount(j & z)} |j ^ x>,
    hence (M psi)[k] = -i i^{n_y} (-1)^{popcount((k ^ x) & z)} psi[k ^ x].  For an odd
    number of Y's the scalar -i i^{n_y} is +-1 and M is a REAL signed permutation.
    """
    dim = 1 << P.n
    k = np.arange(dim, dtype=np.int64)
    perm = k ^ P.x
    par = popcount(perm & P.z) & 1
    c = (-1j) * (1j) ** P.n_y()
    sign = np.where(par == 1, -1.0, 1.0)
    if abs(c.imag) < 1e-15:
        return perm, sign * float(c.real)
    return perm, sign * complex(c)


# ================================================================== the pools
def pool_V(n):
    """Tang et al. (2021) recursive minimal complete pool, 2n-2 strings, qubit k = wire k-1."""
    if n < 2:
        raise ValueError("n >= 2")

    def word(letters):                 # letters: dict qubit(1-based) -> letter
        return "".join(letters.get(k, "I") for k in range(1, n + 1))

    pool = [{2: "Z", 1: "Y"}, {2: "Y"}]            # {V}_2 = {Z2 Y1, Y2}
    for m in range(3, n + 1):
        new = []
        for d in pool:                             # Z_m {V}_{m-1}
            e = dict(d)
            e[m] = "Z"
            new.append(e)
        new.append({m: "Y"})                       # Y_m
        new.append({m - 1: "Y"})                   # Y_{m-1}
        pool = new
    ops = [Pauli.from_word(word(d)) for d in pool]
    assert len(ops) == 2 * n - 2
    return ops


def pool_G(n):
    """Tang et al. (2021) local minimal complete pool: Z_{k+1} Y_k (k=1..n-1) and Y_k (k=2..n)."""
    ops = []
    for k in range(1, n):
        w = ["I"] * n
        w[k - 1] = "Y"
        w[k] = "Z"
        ops.append(Pauli.from_word("".join(w)))
    for k in range(2, n + 1):
        w = ["I"] * n
        w[k - 1] = "Y"
        ops.append(Pauli.from_word("".join(w)))
    assert len(ops) == 2 * n - 2
    return ops


def pool_L2(n, complex_pool=False):
    """All 1- and 2-local strings on any pair; real (odd-Y) unless `complex_pool`."""
    ops = []
    singles = "XYZ" if complex_pool else "Y"
    for q in range(n):
        for a in singles:
            w = ["I"] * n
            w[q] = a
            ops.append(Pauli.from_word("".join(w)))
    pairs = [(a, b) for a in "XYZ" for b in "XYZ"]
    if not complex_pool:
        pairs = [(a, b) for a, b in pairs if (a == "Y") != (b == "Y")]   # exactly one Y
    for q in range(n):
        for r in range(q + 1, n):
            for a, b in pairs:
                w = ["I"] * n
                w[q], w[r] = a, b
                ops.append(Pauli.from_word("".join(w)))
    if complex_pool:
        assert len(ops) == 3 * n + 9 * n * (n - 1) // 2
    else:
        assert len(ops) == 2 * n * n - n
    return ops


def make_pool(name, n):
    if name == "V":
        return pool_V(n)
    if name == "G":
        return pool_G(n)
    if name == "L2":
        return pool_L2(n, complex_pool=False)
    if name == "L2C":
        return pool_L2(n, complex_pool=True)
    raise KeyError(name)


POOLS = ("V", "G", "L2", "L2C")


# ============================================================ the simulator
class PauliCircuit:
    """A product of Pauli-string rotations on |0...0>, batched over parameter vectors.

    ``ops[k]`` is applied k-th (``ops[0]`` first); ``theta[k]`` is its angle.
    ``states_batch(TH)`` with TH of shape (B, P) returns (B, 2^n) amplitudes; real when
    every operator has an odd number of Y's, complex otherwise.
    """

    def __init__(self, n, ops):
        self.n = int(n)
        self.dim = 1 << self.n
        self.ops = list(ops)
        self._act = [pauli_action(P) for P in self.ops]
        self.real = all(np.isrealobj(c) for _, c in self._act)

    def n_params(self):
        return len(self.ops)

    def append(self, P: Pauli):
        self.ops.append(P)
        self._act.append(pauli_action(P))
        self.real = self.real and np.isrealobj(self._act[-1][1])

    def states_batch(self, TH):
        TH = np.atleast_2d(np.asarray(TH, float))
        B = TH.shape[0]
        if TH.shape[1] != len(self.ops):
            raise ValueError(f"{TH.shape[1]} angles for {len(self.ops)} operators")
        psi = np.zeros((B, self.dim), float if self.real else complex)
        psi[:, 0] = 1.0
        C = np.cos(TH / 2.0)
        S = np.sin(TH / 2.0)
        for k, (perm, coef) in enumerate(self._act):
            mpsi = coef[None, :] * psi[:, perm]
            psi = C[:, k][:, None] * psi + S[:, k][:, None] * mpsi
        return psi

    def state(self, theta):
        return self.states_batch(np.asarray(theta, float)[None, :])[0]

    def probs_batch(self, TH):
        psi = self.states_batch(TH)
        p = (psi * np.conj(psi)).real if not self.real else psi * psi
        return p / p.sum(1, keepdims=True)

    def probs(self, theta):
        return self.probs_batch(np.asarray(theta, float)[None, :])[0]

    def apply_M(self, psi, P: Pauli):
        """``(-iP) psi`` for a single state."""
        perm, coef = pauli_action(P)
        return coef * psi[perm]

    def shift_grid(self, theta, h):
        th = np.asarray(theta, float)
        P = th.size
        TH = np.repeat(th[None, :], 2 * P, axis=0)
        r = np.arange(P)
        TH[2 * r, r] += h
        TH[2 * r + 1, r] -= h
        return TH


def ry_layer_ops(n):
    """Y_q for q = 0..n-1: the initial RY layer as Pauli rotations (RY = exp(-i theta Y/2))."""
    ops = []
    for q in range(n):
        w = ["I"] * n
        w[q] = "Y"
        ops.append(Pauli.from_word("".join(w)))
    return ops


def entangler_perm(n, ring=True):
    """The CNOT chain (+ ring) as one basis-index gather, composed here independently."""
    dim = 1 << n
    pairs = [(q, q + 1) for q in range(n - 1)]
    if ring and n > 2:
        pairs.append((n - 1, 0))
    cur = np.arange(dim, dtype=np.int64)
    idx = np.arange(dim, dtype=np.int64)
    for c, t in pairs:
        cb = (idx >> (n - 1 - c)) & 1
        f = np.where(cb == 1, idx ^ (1 << (n - 1 - t)), idx)
        cur = cur[f]
    return cur


class FixedAnsatzAlt:
    """The deployed RY/CNOT-chain-plus-ring ansatz, rebuilt on the Pauli simulator.

    Used ONLY to check the simulator against core.quantum.StatevectorCircuit; the
    harness's fixed arm is core.quantum.run_cvar_vqe itself.
    """

    def __init__(self, n, layers=3, ring=True):
        self.n, self.layers = int(n), int(layers)
        self.dim = 1 << self.n
        self._perm = entangler_perm(self.n, ring)
        self._ry = PauliCircuit(self.n, ry_layer_ops(self.n))

    def n_params(self):
        return self.n * self.layers

    def states_batch(self, TH):
        TH = np.atleast_2d(np.asarray(TH, float)).reshape(-1, self.layers, self.n)
        B = TH.shape[0]
        psi = np.zeros((B, self.dim))
        psi[:, 0] = 1.0
        C = np.cos(TH / 2.0)
        S = np.sin(TH / 2.0)
        for L in range(self.layers):
            for q, (perm, coef) in enumerate(self._ry._act):
                psi = C[:, L, q][:, None] * psi + S[:, L, q][:, None] * (coef[None, :] * psi[:, perm])
            psi = psi[:, self._perm]
        return psi

    def probs(self, theta):
        psi = self.states_batch(np.asarray(theta, float)[None, :])[0]
        p = psi * psi
        return p / p.sum()


# ============================================================ the objective
def objective(E, p, alpha, T):
    """``(F, d)`` for F = CVaR_alpha(E; p) - T H(p), d = dF/dp, exactly as core.quantum.free_energy."""
    from core import quantum as Q
    v, q, dp = Q.cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-15))
    H = float(-(p * lp).sum())
    dH = -(lp + 1.0)
    return float(v - T * H), dp - T * dH, float(v), H


def free_energy_grad(circ: PauliCircuit, theta, E, alpha, T):
    """``(F, grad F, p)`` by the exact two-term parameter-shift rule on every basis probability."""
    theta = np.asarray(theta, float)
    p = circ.probs(theta)
    F, d, v, H = objective(E, p, alpha, T)
    PR = circ.probs_batch(circ.shift_grid(theta, np.pi / 2))
    g = (PR[0::2] - PR[1::2]) @ d / 2.0
    return F, g, p, v, H


def pool_gradients(circ: PauliCircuit, theta, E, alpha, T, pool):
    """dF/dphi at phi = 0 for exp(-i phi P/2) appended at the END, for every P in the pool.

    Exact: (p(+pi/2) - p(-pi/2))/2 = Re(conj(psi) (M psi)) elementwise, M = -iP.
    """
    psi = circ.state(theta)
    p = (psi * np.conj(psi)).real
    p = p / p.sum()
    _, d, _, _ = objective(E, p, alpha, T)
    out = np.empty(len(pool))
    for k, P in enumerate(pool):
        mpsi = circ.apply_M(psi, P)
        out[k] = float((d * (np.conj(psi) * mpsi).real).sum())
    return out


def pool_gradients_shift(circ: PauliCircuit, theta, E, alpha, T, pool):
    """The same quantity by literally simulating the +-pi/2 shifted circuits (a check)."""
    theta = np.asarray(theta, float)
    p0 = circ.probs(theta)
    _, d, _, _ = objective(E, p0, alpha, T)
    out = np.empty(len(pool))
    for k, P in enumerate(pool):
        c2 = PauliCircuit(circ.n, circ.ops + [P])
        pp = c2.probs(np.r_[theta, np.pi / 2])
        pm = c2.probs(np.r_[theta, -np.pi / 2])
        out[k] = float(((pp - pm) / 2.0) @ d)
    return out


def adam_optimise(circ: PauliCircuit, theta, E, alpha, T, steps, lr=LR, best=False):
    """Adam on the exact gradient, the arithmetic of core.quantum.run_cvar_vqe, fresh state.

    ``best=False`` returns the LAST iterate (the deployed driver's convention).
    ``best=True`` returns the iterate with the lowest F among the start and every step.
    Measured on the deployed energy shape (s26/agentQ_FINDINGS.md section 0): a fresh Adam at
    lr 0.15 moves every angle by 0.15 rad on its first step whatever the gradient, so F rises
    before it falls and the last iterate is not stationary (|g| ~ 7e-2); the best iterate
    makes growth monotone by construction and is the pre-registered ADAPT re-optimiser.
    """
    th = np.asarray(theta, float).copy()
    m = np.zeros_like(th)
    v = np.zeros_like(th)
    best_F, best_th = np.inf, th.copy()
    for t in range(1, int(steps) + 1):
        F, g, _, _, _ = free_energy_grad(circ, th, E, alpha, T)
        if best and F < best_F:
            best_F, best_th = F, th.copy()
        m = ADAM_B1 * m + (1 - ADAM_B1) * g
        v = ADAM_B2 * v + (1 - ADAM_B2) * g * g
        th = th - lr * (m / (1 - ADAM_B1 ** t)) / (np.sqrt(v / (1 - ADAM_B2 ** t)) + ADAM_EPS)
    if not best:
        return th
    F = free_energy_grad(circ, th, E, alpha, T)[0]
    if F < best_F:
        best_F, best_th = F, th.copy()
    return best_th


def lbfgs_optimise(circ: PauliCircuit, theta, E, alpha, T, maxiter=300, gtol=1e-7):
    """L-BFGS-B on the exact gradient, the ADAPT literature's re-optimiser (Grimsley 2019,
    Tang 2021 use BFGS).  Converges to |g| ~ 1e-7 where F is smooth (alpha = 1) and stops at
    a CVaR kink otherwise; the returned point is the best F seen by the line search."""
    from scipy.optimize import minimize

    def fg(x):
        F, g, _, _, _ = free_energy_grad(circ, x, E, alpha, T)
        return F, g

    x0 = np.asarray(theta, float)
    res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                   options=dict(maxiter=int(maxiter), gtol=gtol))
    F0 = fg(x0)[0]
    return res.x if res.fun <= F0 else x0


OPTIMISERS = ("adam_best", "adam_last", "lbfgs")


def reoptimise(circ, theta, E, alpha, T, optimiser, steps, lr=LR):
    if optimiser == "adam_best":
        return adam_optimise(circ, theta, E, alpha, T, steps, lr=lr, best=True)
    if optimiser == "adam_last":
        return adam_optimise(circ, theta, E, alpha, T, steps, lr=lr, best=False)
    if optimiser == "lbfgs":
        return lbfgs_optimise(circ, theta, E, alpha, T, maxiter=max(int(steps), 300))
    raise KeyError(optimiser)


# ================================================================ ADAPT-VQE
def run_adapt(E, alpha, T, n, pool_name, max_params=None, adam_steps=50, seed=0,
              eps=1e-3, record_at=(7, 14, 21), lr=LR, init_steps=None, optimiser="adam_best"):
    """qubit-ADAPT-VQE on F = CVaR_alpha(E; p) - T H(p).  Returns a dict (see below).

    Start: one RY layer, theta ~ N(0, INIT_SD^2) from default_rng(seed) (the deployed
    law), re-optimised for `init_steps` (default `adam_steps`) steps -> the P = n state.
    Then grow until ||g_pool|| < eps or len(theta) == max_params (default 3n, i.e.
    parameter-matched to the deployed layers=3 ansatz).  `optimiser` is one of
    OPTIMISERS; every re-optimisation starts from the previous angles with the new angle
    at 0 and a fresh optimiser state.

    If growth stops early (eps), `snapshots` carries the final state under every
    remaining key of `record_at`, each flagged `stopped_early=True` with its actual P.
    """
    E = np.asarray(E, float)
    n = int(n)
    max_params = 3 * n if max_params is None else int(max_params)
    init_steps = adam_steps if init_steps is None else int(init_steps)
    if optimiser not in OPTIMISERS:
        raise KeyError(optimiser)
    pool = make_pool(pool_name, n)
    rng = np.random.default_rng(seed)
    theta = rng.normal(0.0, INIT_SD, n)
    circ = PauliCircuit(n, ry_layer_ops(n))
    theta = reoptimise(circ, theta, E, alpha, T, optimiser, init_steps, lr=lr)
    trace, seq, snaps = [], [], {}
    total_adam = int(init_steps)
    tie_events = 0
    consecutive_repeats = 0
    stopped = "max_params"

    def snapshot():
        F, g, p, v, H = free_energy_grad(circ, theta, E, alpha, T)
        return dict(P=int(len(theta)), F=F, cvar=v, H_nats=H, H_bits=H / math.log(2.0),
                    grad_norm=float(np.linalg.norm(g)), p=p, stopped_early=False)

    s0 = snapshot()
    trace.append(dict(step=0, P=s0["P"], F=s0["F"], H_bits=s0["H_bits"], op=None,
                      pool_grad=None, pool_grad_norm=None))
    if s0["P"] in record_at:
        snaps[s0["P"]] = s0
    step = 0
    while len(theta) < max_params:
        gp = pool_gradients(circ, theta, E, alpha, T, pool)
        norm = float(np.linalg.norm(gp))
        if norm < eps:
            stopped = "eps"
            break
        a = np.abs(gp)
        best = float(a.max())
        ties = np.flatnonzero(np.isclose(a, best, rtol=1e-12, atol=0.0))
        if len(ties) > 1:
            tie_events += 1
        k = int(ties[0])                                  # lowest pool index on a tie
        if seq and pool[k].word() == seq[-1]:
            consecutive_repeats += 1
        circ.append(pool[k])
        theta = np.r_[theta, 0.0]
        theta = reoptimise(circ, theta, E, alpha, T, optimiser, adam_steps, lr=lr)
        total_adam += int(adam_steps)
        step += 1
        s = snapshot()
        seq.append(pool[k].word())
        trace.append(dict(step=step, P=s["P"], F=s["F"], H_bits=s["H_bits"], op=pool[k].word(),
                          pool_grad=float(gp[k]), pool_grad_norm=norm, n_ties=int(len(ties))))
        if s["P"] in record_at:
            snaps[s["P"]] = s
    final = snapshot()
    for P in record_at:
        if P not in snaps and P > len(theta):
            snaps[P] = dict(final, stopped_early=True)
    return dict(pool=pool_name, n=n, alpha=float(alpha), T=float(T), seed=int(seed),
                optimiser=optimiser, eps=float(eps), adam_steps=int(adam_steps),
                init_steps=int(init_steps), total_adam_steps=total_adam,
                max_params=max_params, stopped=stopped, tie_events=tie_events,
                consecutive_repeats=consecutive_repeats,
                n_distinct_ops=len(set(seq)), ops=[P.word() for P in circ.ops], sequence=seq,
                theta=theta.tolist(), trace=trace, snapshots=snaps, final=final,
                pool_size=len(pool))


def product_of_marginals(p, n):
    """The product of the n single-bit marginals of p (qubit 0 = MSB)."""
    p = np.asarray(p, float) / np.sum(p)
    idx = np.arange(1 << n)
    out = np.ones(1 << n)
    for q in range(n):
        bit = (idx >> (n - 1 - q)) & 1
        m1 = float(p[bit == 1].sum())
        out = out * np.where(bit == 1, m1, 1.0 - m1)
    return out


def product_diagnostics(E, T, n):
    """How far the Gibbs state of E at temperature T is from a product state.

    For the DEPLOYED E (standardised ranks; index = rank) E is affine in the register index
    j = sum_q 2^(n-1-q) b_q, so exp(-E_j/T) factorises over the bits and the Gibbs state IS a
    product state: KL(Gibbs || product of its marginals) = 0 up to tie-averaging.  Measured
    here, not assumed.  Returned in nats.
    """
    g = gibbs(E, T)
    pm = product_of_marginals(g, n)
    return dict(kl_gibbs_to_product=kl(g, pm), kl_product_to_gibbs=kl(pm, g),
                H_gibbs_bits=float(-(g * np.log2(np.maximum(g, 1e-300))).sum()),
                H_product_bits=float(-(pm * np.log2(np.maximum(pm, 1e-300))).sum()))


def gibbs(E, T):
    w = -np.asarray(E, float) / T
    w = w - w.max()
    g = np.exp(w)
    return g / g.sum()


def entropy_nats(p):
    p = np.asarray(p, float) / np.sum(p)
    return float(-(p * np.log(np.maximum(p, 1e-300))).sum())


def match_temperature(E, H_target_nats, lo=1e-3, hi=1e3, tol=1e-10):
    """The T at which the Gibbs entropy of E equals H_target (monotone in T; bisection)."""
    E = np.asarray(E, float)
    if entropy_nats(gibbs(E, hi)) < H_target_nats:
        return float(hi)
    if entropy_nats(gibbs(E, lo)) > H_target_nats:
        return float(lo)
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if entropy_nats(gibbs(E, mid)) < H_target_nats:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1.0 < tol:
            break
    return float(math.sqrt(lo * hi))


def random_weights_matched_entropy(H_nats, dim, seed, tol=1e-9):
    """THE MATCHED RANDOM CONTROL: random weights over the SAME candidates at the SAME entropy.

    w_c = softmax(c * g) with g ~ N(0, 1) a fixed draw per (target seed) and c >= 0 bisected
    so that H(w_c) = H_nats.  H(w_c) falls monotonically from log(dim) at c = 0 to 0, so the
    match is unique.  S25 (q_alpha.json, entropy_curve) showed the readout responds to the
    entropy of the weights and to nothing else that was measured; this control carries that
    entropy and no information about the candidates.
    """
    rng = np.random.default_rng(int(seed))
    g = rng.normal(0.0, 1.0, int(dim))
    H_nats = float(min(max(H_nats, 0.0), math.log(dim)))

    def w_of(c):
        x = c * g
        x = x - x.max()
        w = np.exp(x)
        return w / w.sum()

    lo, hi = 0.0, 1.0
    while entropy_nats(w_of(hi)) > H_nats and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if entropy_nats(w_of(mid)) > H_nats:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return w_of(0.5 * (lo + hi))


def kl(p, q):
    p = np.asarray(p, float) / np.sum(p)
    q = np.asarray(q, float) / np.sum(q)
    return float((p * np.log(np.maximum(p, 1e-300) / np.maximum(q, 1e-300))).sum())


def sym_kl(p, q):
    return max(kl(p, q), kl(q, p))


# ============================================================ energy variants (A3)
def energy_variants(sc_top):
    """Monotone (order-preserving) standardisations of the score prefix sc_top (ascending).

    zrank   the DEPLOYED E: standardised ranks, ties averaged (core.pipeline._zrank).
    zraw    moment standardisation of the raw Bayes-risk score: (s - mean)/sd.  Monotone
            affine; the gaps between candidates survive, so E varies between targets.
    asinh   S21 L33's `Nt`: asinh of a robust z-score, then re-standardised (mean 0, sd 1).
    soft    S13 rank-preserving soft compression of the upper tail (geo_common.soft_compress),
            then re-standardised.
    Every variant is a strictly increasing function of the score wherever the score has no
    exact ties, so the CVaR tail is the same prefix of the same order (the S25 theorem).
    """
    from core.pipeline import _zrank
    s = np.asarray(sc_top, float)
    out = {"zrank": _zrank(s)}
    sd = s.std()
    out["zraw"] = (s - s.mean()) / max(sd, 1e-12)
    med = np.median(s)
    mad = np.median(np.abs(s - med))
    scale = 1.4826 * mad if mad > 0 else max(sd, 1e-12)
    a = np.arcsinh((s - med) / scale)
    out["asinh"] = (a - a.mean()) / max(a.std(), 1e-12)
    q = float(np.median(s))
    iqr = float(np.quantile(s, 0.75) - np.quantile(s, 0.25))
    iqr = iqr if iqr > 0 else 1.0
    c = s.copy()
    hi = s > q
    c[hi] = q + iqr * np.log1p((s[hi] - q) / iqr)
    out["soft"] = (c - c.mean()) / max(c.std(), 1e-12)
    for k, v in out.items():
        d = np.diff(v)
        if (d < -1e-12).any():
            raise RuntimeError(f"variant {k} is not monotone in the score")
    return out


# ================================================================== harness
ARM_SPEC = None      # filled by build_arm_spec()


def build_arm_spec(pools=("V", "L2"), record_at=(7, 14, 21), fixed_iters=(50, 750),
                   variants=("zrank",), adapt_variants=("zrank",), optimisers=("adam_best", "lbfgs")):
    """Names of every state the harness produces per target; the prereg lists them."""
    arms = []
    for ev in variants:
        for it in fixed_iters:
            arms.append(f"fixed_{ev}_it{it}")
    for ev in adapt_variants:
        for pn in pools:
            for opt in optimisers:
                for P in record_at:
                    arms.append(f"adapt{pn}_{opt}_{ev}_P{P}")
    arms.append("uniform128")
    arms.append("gibbs_T")
    arms.append("randH_fixed")
    if pools and optimisers:
        arms.append("randH_adaptL2")
    for ev in variants:
        if ev != "zrank":
            arms.append(f"fixed_{ev}_Tmatch_it{int(fixed_iters[0])}")
    return arms


def distinct_states(P, thresh_nats=0.01):
    """Greedy count of distinct distributions: a state is new if its symmetric KL to every
    representative so far exceeds `thresh_nats` (representatives in the given order)."""
    reps = []
    assign = []
    for p in P:
        for j, q in enumerate(reps):
            if sym_kl(p, q) <= thresh_nats:
                assign.append(j)
                break
        else:
            reps.append(np.asarray(p, float))
            assign.append(len(reps) - 1)
    return len(reps), assign


def _targets(cfg, limit=None, order="fold", pdbs=None):
    import core
    from core import pipeline as pl
    db = core.backend("data")
    tg = pl.manifest("tuning126")
    folds = db.folds(cfg.n_folds)
    if pdbs:
        want = set(pdbs)
        tg = [p for p in tg if p.pdb in want]
    if order == "fold":
        tg = sorted(tg, key=lambda p: (folds[p.seq], p.n, p.pdb))
    elif order == "fold_rev":
        tg = sorted(tg, key=lambda p: (-folds[p.seq], -p.n, p.pdb))
    else:
        tg = sorted(tg, key=lambda p: p.pdb)
    if limit:
        tg = tg[:int(limit)]
    return tg, folds


def build_one(target, fold, cfg, pools, record_at, fixed_iters, variants, adapt_variants,
              seed, adam_steps, eps, check_cache=True, verbose=True,
              optimisers=("adam_best", "lbfgs")):
    """Everything native-free for one target.  Returns the record (coordinates, states)."""
    import core
    from core import pipeline as pl
    from core import quantum as Q
    clk = pl.Clock()
    t0 = time.perf_counter()
    pool = pl.retrieve(target, fold, cfg, clk)
    pool = pl.score(pool, target.seq, fold, target.n, cfg, clk)
    sub, Wsub, Ps, top, Pt = pl.filter_pool(pool, cfg, clk)
    C, b = pl.average(Wsub, Ps, clk)
    ca, phi, psi, fit = pl.project(C, target.seq, fold, cfg, clk)
    dim = 1 << cfg.vqe_qubits
    if len(top) < dim:
        raise ValueError(f"{target.pdb}: pool has {len(top)} < {dim}")
    o = np.asarray(top[:dim], int)
    sc = np.asarray(pool["sc"], float)[o]
    block = np.asarray(Pt, float)[:dim, :dim]
    Wo = pool["W"][o]
    alpha, T = pl.VQE_LFO[int(fold) % len(pl.VQE_LFO)]
    EV = energy_variants(sc)
    rec = {"pdb": target.pdb, "n": int(target.n), "fold": int(fold), "seq": target.seq,
           "alpha": float(alpha), "T": float(T), "seed": int(seed),
           "sub": np.asarray(sub, int).tolist(), "o": o.tolist(),
           "sc_top": sc.tolist(), "E_variants": {k: v.tolist() for k, v in EV.items()},
           "avg_ca": np.asarray(C, float).tolist(), "fit_ca": np.asarray(fit, float).tolist(),
           "ca": np.asarray(ca, float).tolist(), "arms": {}, "adapt": {}}
    states = {}
    # ---- fixed ansatz arms (core.quantum, the deployed driver) ------------------
    for ev in variants:
        for it in fixed_iters:
            p, cv, H, _ = Q.run_cvar_vqe(EV[ev], alpha, T, n=cfg.vqe_qubits, layers=cfg.vqe_layers,
                                         iters=int(it), seed=int(seed))
            states[f"fixed_{ev}_it{it}"] = dict(p=np.asarray(p, float), cvar=float(cv), H_nats=float(H))
    # ---- ADAPT arms ---------------------------------------------------------------
    for ev in adapt_variants:
        for pn in pools:
            for opt in optimisers:
                r = run_adapt(EV[ev], alpha, T, cfg.vqe_qubits, pn, max_params=3 * cfg.vqe_qubits,
                              adam_steps=adam_steps, seed=int(seed), eps=eps, record_at=record_at,
                              optimiser=opt)
                key = f"{pn}_{opt}_{ev}"
                rec["adapt"][key] = {k: v for k, v in r.items() if k != "snapshots"}
                rec["adapt"][key]["final"] = {k: v for k, v in r["final"].items() if k != "p"}
                for P, s in r["snapshots"].items():
                    states[f"adapt{pn}_{opt}_{ev}_P{P}"] = dict(
                        p=np.asarray(s["p"], float), cvar=s["cvar"], H_nats=s["H_nats"], F=s["F"],
                        grad_norm=s["grad_norm"], P_actual=s["P"], stopped_early=s["stopped_early"])
    states["uniform128"] = dict(p=np.ones(dim) / dim, cvar=None, H_nats=math.log(dim))
    gT = gibbs(EV["zrank"], T)
    states["gibbs_T"] = dict(p=gT, cvar=None, H_nats=float(-(gT * np.log(gT)).sum()))
    # ---- the matched random controls: same 128 candidates, same entropy, no information --
    import zlib
    tseed = zlib.crc32(target.pdb.encode()) ^ (int(seed) * 7919)
    primary_opt = "adam_best" if "adam_best" in optimisers else (optimisers[0] if optimisers else None)
    for ref_name, tag in (("fixed_zrank_it50", "randH_fixed"),
                          (f"adaptL2_{primary_opt}_zrank_P21" if primary_opt else None, "randH_adaptL2")):
        if ref_name in states:
            w = random_weights_matched_entropy(states[ref_name]["H_nats"], dim, tseed)
            states[tag] = dict(p=w, cvar=None, H_nats=entropy_nats(w), matched_to=ref_name)
    # ---- A3: the fixed ansatz on each variant at a TEMPERATURE that matches the deployed
    # Gibbs entropy (native-free bisection), so target-dependent gaps are separated from
    # sharper weights (S25: the readout responds to the entropy of the weights) ---------
    rec["T_matched"] = {}
    H_dep = entropy_nats(gibbs(EV["zrank"], T))
    for ev in variants:
        if ev == "zrank":
            continue
        T_ev = match_temperature(EV[ev], H_dep)
        rec["T_matched"][ev] = T_ev
        p, cv, H, _ = Q.run_cvar_vqe(EV[ev], alpha, T_ev, n=cfg.vqe_qubits, layers=cfg.vqe_layers,
                                     iters=int(fixed_iters[0]), seed=int(seed))
        states[f"fixed_{ev}_Tmatch_it{int(fixed_iters[0])}"] = dict(
            p=np.asarray(p, float), cvar=float(cv), H_nats=float(H), T_used=float(T_ev))
    # ---- how far each variant's Gibbs state is from a product state (native-free) ------
    rec["product_diagnostics"] = {ev: product_diagnostics(EV[ev], T, cfg.vqe_qubits) for ev in EV}
    rec["product_diagnostics"]["ideal_ladder"] = product_diagnostics(
        (lambda r: (r - r.mean()) / r.std())(np.arange(1, dim + 1, dtype=float)), T, cfg.vqe_qubits)
    # ---- readouts: selection medoid and the built chain, exactly as production ------
    for name, st in states.items():
        p = st["p"]
        local = pl.consensus_medoid(block, p)
        Cq = pl.average_weighted(Wo, block, p, clk)
        qa, qp, qz, _ = pl.project(Cq, target.seq, fold, cfg, clk)
        arm = {"p": p.tolist(), "sel": int(o[local]), "sel_local": int(local),
               "q_avg_ca": np.asarray(Cq, float).tolist(), "q_ca": np.asarray(qa, float).tolist(),
               "entropy_bits": float(st["H_nats"]) / math.log(2.0), "cvar": st["cvar"],
               "kl_to_gibbs": kl(p, gT), "kl_from_gibbs": kl(gT, p),
               "kl_to_product": kl(p, product_of_marginals(p, cfg.vqe_qubits))}
        for k in ("F", "grad_norm", "P_actual", "stopped_early", "matched_to", "T_used"):
            if k in st:
                arm[k] = st[k] if isinstance(st[k], (bool, str)) else float(st[k])
        rec["arms"][name] = arm
    rec["arms"]["uniform128"]["sel_uniform"] = int(o[pl.consensus_medoid(block)])
    # ---- the bit-for-bit check against the production quantum cache ----------------
    if check_cache:
        path = os.path.join(QCACHE, f"{target.pdb}.json")
        if os.path.exists(path):
            with open(path) as fh:
                ref = json.load(fh)
            d_ca = float(np.abs(np.asarray(ref["ca"]) - np.asarray(rec["ca"])).max())
            d_q = float(np.abs(np.asarray(ref["q_ca"]) - np.asarray(rec["arms"]["fixed_zrank_it50"]["q_ca"])).max())
            same_sel = int(ref["quantum"]["sel"]) == rec["arms"]["fixed_zrank_it50"]["sel"]
            rec["cache_check"] = {"ca_max_abs": d_ca, "q_ca_max_abs": d_q, "sel_equal": bool(same_sel),
                                  "alpha_T_equal": (float(ref["quantum"]["alpha"]), float(ref["quantum"]["T"])) == (alpha, T),
                                  "entropy_bits_ref": ref["quantum"]["entropy_bits"],
                                  "entropy_bits_here": rec["arms"]["fixed_zrank_it50"]["entropy_bits"]}
        else:
            rec["cache_check"] = None
    rec["timings"] = {k: round(v, 4) for k, v in clk.t.items()}
    rec["wall"] = round(time.perf_counter() - t0, 3)
    if verbose:
        cc = rec.get("cache_check") or {}
        print(f"  {target.pdb} n={target.n} fold={fold} a={alpha} T={T} wall {rec['wall']:.1f}s "
              f"cache ca {cc.get('ca_max_abs', 'n/a')} q_ca {cc.get('q_ca_max_abs', 'n/a')} "
              f"sel_equal {cc.get('sel_equal', 'n/a')}", flush=True)
    return rec


def label_one(rec, target, cfg):
    """ORACLE: the native is read here and nowhere else.  Gated on PHASE 0 sign-off."""
    import core
    from core import pipeline as pl
    aud = core.backend("numerics")
    from s8.inband import sel_of
    # the pool's per-window RMSD, exactly as core.pipeline.label builds it
    from core import pipeline as _pl
    clk = pl.Clock()
    pool = pl.retrieve(target, rec["fold"], cfg, clk)
    pool = pl.score(pool, target.seq, rec["fold"], target.n, cfg, clk)
    rr = pl._q(aud.kabsch_rmsd_batch(pool["W64"], np.asarray(target.ca, float)), cfg)
    nat = pl._q(np.asarray(target.ca, float), cfg)
    sub = np.asarray(rec["sub"], int)
    out = {"shipped": float(sel_of(pool["sc"], rr)), "pool_best": float(rr.min()),
           "top_m_best": float(rr[sub].min()),
           "rmsd_avg": float(aud.kabsch_rmsd_batch(np.asarray(rec["avg_ca"])[None], nat)[0]),
           "rmsd_arm": float(aud.kabsch_rmsd_batch(np.asarray(rec["ca"])[None], nat)[0]),
           "arms": {}}
    for name, a in rec["arms"].items():
        out["arms"][name] = {
            "rmsd_sel": float(rr[int(a["sel"])]),
            "rmsd_q_avg": float(aud.kabsch_rmsd_batch(np.asarray(a["q_avg_ca"])[None], nat)[0]),
            "rmsd_q_synth": float(aud.kabsch_rmsd_batch(np.asarray(a["q_ca"])[None], nat)[0])}
    out["arms"]["uniform128"]["rmsd_sel_uniform"] = float(rr[int(rec["arms"]["uniform128"]["sel_uniform"])])
    return out


def _atomic_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def _signed_off():
    try:
        with open(os.path.join(HERE, "LEDGER.md"), encoding="utf-8") as fh:
            return "PHASE 0 SIGNED OFF" in fh.read()
    except OSError:
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", action="store_true", help="native-free build of every arm per target")
    ap.add_argument("--label", action="store_true", help="ORACLE labelling (gated on sign-off)")
    ap.add_argument("--stats", action="store_true", help="the pre-registered contrasts")
    ap.add_argument("--tag", default="a1", help="results sub-directory / merged file tag")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--order", default="fold", choices=("fold", "fold_rev", "pdb"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--pdbs", default=None, help="comma-separated pdb ids")
    ap.add_argument("--pools", default="V,L2")
    ap.add_argument("--variants", default="zrank", help="energy variants for the FIXED ansatz")
    ap.add_argument("--adapt-variants", default="zrank", help="energy variants for ADAPT")
    ap.add_argument("--fixed-iters", default="50,750")
    ap.add_argument("--optimisers", default="adam_best,lbfgs")
    ap.add_argument("--adam-steps", type=int, default=50)
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--record-at", default="7,14,21")
    ap.add_argument("--no-check-cache", action="store_true")
    ap.add_argument("--i-am-signed-off", action="store_true",
                    help="override the ledger check (coordinator instruction only)")
    a = ap.parse_args(argv)

    import core
    from core import pipeline as pl
    cfg = replace(pl.PROD, quantum=True, amber=False, legacy=False)
    pools = tuple(x for x in a.pools.split(",") if x)
    variants = tuple(x for x in a.variants.split(",") if x)
    adapt_variants = tuple(x for x in a.adapt_variants.split(",") if x)
    fixed_iters = tuple(int(x) for x in a.fixed_iters.split(",") if x)
    optimisers = tuple(x for x in a.optimisers.split(",") if x)
    record_at = tuple(int(x) for x in a.record_at.split(",") if x)
    outdir = os.path.join(RESULTS, a.tag)
    os.makedirs(outdir, exist_ok=True)
    pdbs = a.pdbs.split(",") if a.pdbs else None
    tg, folds = _targets(cfg, limit=a.limit, order=a.order, pdbs=pdbs)
    pl.guard_esm([p.seq for p in pl.manifest("tuning126")])
    pl._LIB.n_folds = cfg.n_folds
    pl._LIB.cap = cfg.window_cache
    pl._LIB.cap_members = cfg.member_cache

    if a.build:
        print(f"BUILD tag={a.tag} seed={a.seed} order={a.order} n_targets={len(tg)} pools={pools} "
              f"variants={variants} adapt_variants={adapt_variants} fixed_iters={fixed_iters}", flush=True)
        for k, t in enumerate(tg):
            path = os.path.join(outdir, f"{t.pdb}.json")
            if os.path.exists(path):
                try:
                    with open(path) as fh:
                        done = json.load(fh)
                    if done.get("seed") == a.seed and set(done.get("arms", {})) >= set(
                            build_arm_spec(pools, record_at, fixed_iters, variants, adapt_variants,
                                           optimisers)):
                        print(f"  [{k + 1}/{len(tg)}] {t.pdb} cached", flush=True)
                        continue
                except Exception:                                        # noqa: BLE001
                    pass
            if not pl.wait_for_memory(cfg.stop_pct, where=t.pdb):
                print("memory gate: stopping", flush=True)
                return 2
            rec = build_one(t, folds[t.seq], cfg, pools, record_at, fixed_iters, variants,
                            adapt_variants, a.seed, a.adam_steps, a.eps,
                            check_cache=not a.no_check_cache, optimisers=optimisers)
            rec["provenance_seed"] = a.seed
            _atomic_json(path, rec)
            print(f"  [{k + 1}/{len(tg)}] {t.pdb} written", flush=True)
        return 0

    if a.label:
        if not (_signed_off() or a.i_am_signed_off):
            print("REFUSED: --label reads the native; s26/LEDGER.md does not contain "
                  "'PHASE 0 SIGNED OFF'.", flush=True)
            return 3
        db = core.backend("data")
        allt = {p.pdb: p for p in db.load()}
        for k, t in enumerate(tg):
            path = os.path.join(outdir, f"{t.pdb}.json")
            if not os.path.exists(path):
                print(f"  {t.pdb}: no build record, skipped", flush=True)
                continue
            with open(path) as fh:
                rec = json.load(fh)
            if rec.get("label") is not None and rec["label"].get("arms") and \
                    set(rec["label"]["arms"]) >= set(rec["arms"]):
                continue
            rec["label"] = label_one(rec, allt[t.pdb], cfg)
            _atomic_json(path, rec)
            print(f"  [{k + 1}/{len(tg)}] {t.pdb} labelled", flush=True)
        return 0

    if a.stats:
        return stats(outdir, a.tag, cfg)
    ap.print_help()
    return 1


def stats(outdir, tag, cfg):
    """The pre-registered contrasts, through s24.stats_lib.compare, fold-clustered."""
    from s24 import stats_lib as ST
    rows = []
    for f in sorted(os.listdir(outdir)):
        if f.endswith(".json"):
            with open(os.path.join(outdir, f)) as fh:
                rows.append(json.load(fh))
    lab = [r for r in rows if r.get("label")]
    print(f"{len(rows)} build records, {len(lab)} labelled")
    # ---- property summary (no native): free energies, KL to Gibbs, distinct states and
    # operator sequences per (alpha, T) cell, cache checks -------------------------------
    prop = {"tag": tag, "n": len(rows), "rows": [], "distinct_states": {}, "distinct_sequences": {}}
    for r in rows:
        row = {"pdb": r["pdb"], "fold": r["fold"], "alpha": r["alpha"], "T": r["T"],
               "cache_check": r.get("cache_check"), "wall": r.get("wall"),
               "product_diagnostics": r.get("product_diagnostics"), "T_matched": r.get("T_matched")}
        for name, a in r["arms"].items():
            row[f"{name}:H_bits"] = a["entropy_bits"]
            row[f"{name}:kl_to_gibbs"] = a["kl_to_gibbs"]
            row[f"{name}:kl_to_product"] = a.get("kl_to_product")
            if "F" in a:
                row[f"{name}:F"] = a["F"]
        for k, ad in r.get("adapt", {}).items():
            row[f"adapt_{k}:sequence"] = ad["sequence"]
            row[f"adapt_{k}:stopped"] = ad["stopped"]
            row[f"adapt_{k}:F_final"] = ad["final"]["F"]
            row[f"adapt_{k}:n_distinct_ops"] = ad["n_distinct_ops"]
        prop["rows"].append(row)
    cells = sorted({(r["alpha"], r["T"]) for r in rows})
    arm_names = sorted(set.intersection(*[set(r["arms"]) for r in rows])) if rows else []
    for name in arm_names:
        per_cell = {}
        for cell in cells:
            P = [np.asarray(r["arms"][name]["p"], float) for r in rows if (r["alpha"], r["T"]) == cell]
            k, assign = distinct_states(P, 0.01)
            per_cell[f"a{cell[0]}_T{cell[1]}"] = dict(n_targets=len(P), n_distinct=k)
        prop["distinct_states"][name] = per_cell
    adapt_keys = sorted(set.intersection(*[set(r.get("adapt", {})) for r in rows])) if rows else []
    for k in adapt_keys:
        per_cell = {}
        for cell in cells:
            seqs = [tuple(r["adapt"][k]["sequence"]) for r in rows if (r["alpha"], r["T"]) == cell]
            per_cell[f"a{cell[0]}_T{cell[1]}"] = dict(n_targets=len(seqs), n_distinct=len(set(seqs)))
        prop["distinct_sequences"][k] = per_cell
    ST.save_atomic(os.path.join(RESULTS, f"{tag}_property.json"), prop, module_file=__file__)
    print("wrote", os.path.join(RESULTS, f"{tag}_property.json"))
    for name, pc in prop["distinct_states"].items():
        print(f"  distinct states {name:34s} " + "  ".join(f"{c}: {v['n_distinct']}/{v['n_targets']}" for c, v in pc.items()))
    for k, pc in prop["distinct_sequences"].items():
        print(f"  distinct sequences {k:30s} " + "  ".join(f"{c}: {v['n_distinct']}/{v['n_targets']}" for c, v in pc.items()))
    if not lab:
        return 0
    pdbs = [r["pdb"] for r in lab]
    folds = np.array([int(r["fold"]) for r in lab])
    names = sorted(lab[0]["label"]["arms"])

    def col(arm, key):
        return np.array([r["label"]["arms"][arm][key] for r in lab], float)

    ref_arm = "fixed_zrank_it50"
    res = {"tag": tag, "n": len(lab), "pdbs": pdbs, "contrasts": {}, "means": {}}
    for arm in names:
        res["means"][arm] = {k: float(col(arm, k).mean()) for k in ("rmsd_q_synth", "rmsd_sel", "rmsd_q_avg")}
    res["means"]["rmsd_arm"] = float(np.mean([r["label"]["rmsd_arm"] for r in lab]))
    res["means"]["shipped"] = float(np.mean([r["label"]["shipped"] for r in lab]))
    for basis in ("rmsd_q_synth", "rmsd_sel"):
        for arm in names:
            if arm == ref_arm:
                continue
            r = ST.compare(col(arm, basis), col(ref_arm, basis), folds, names=pdbs,
                           label=f"{basis}: {arm} - {ref_arm}")
            print(ST.fmt(r))
            res["contrasts"][f"{basis}:{arm}-{ref_arm}"] = r
    ST.save_atomic(os.path.join(RESULTS, f"{tag}_stats.json"), res, module_file=__file__)
    print("wrote", os.path.join(RESULTS, f"{tag}_stats.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
