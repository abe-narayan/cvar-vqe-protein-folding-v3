"""s12 QUANTUM-ROLE -- the CVaR-VQE arm, run on the project's own genuine machinery.

`core.quantum.run_cvar_vqe` is the shipped driver: a real `layers x (RY on every wire, CNOT
chain + ring)` circuit, simulated exactly, optimised by Adam on the EXACT parameter-shift
CVaR gradient (`grad_cvar_paramshift`), with the entropy term of `free_energy` making the
objective a free energy `F = CVaR_alpha - T H`.  Nothing here reimplements any of it.

The recorded gradient DEFECT (`cvar_gradient(baseline="tail")`, cosine +0.656 at 0.758x
norm) is never used: `run_cvar_vqe` differentiates through `free_energy`, whose gradient is
the exact parameter-shift one, and the sampled-estimator audit below is run with
`baseline="const"`.  `vqe_gradient_audit` re-measures both on THESE Hamiltonians so the
finding is not merely quoted.

A `lightning.qubit` cross-check (`pennylane_probs_match`) confirms the exact statevector
used for speed is the same circuit the pinned device simulates.
"""
from __future__ import annotations
import os, sys, time, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def _q():
    from core import quantum as Q
    return Q


def run(E, alpha, T, n, layers=3, iters=80, restarts=1, seed=0, lr=0.15):
    """Wrap `core.quantum.run_cvar_vqe`; return p, cvar, entropy, wall, circuit."""
    Q = _q()
    t0 = time.perf_counter()
    p, cv, H, circ = Q.run_cvar_vqe(np.asarray(E, float), alpha=alpha, T=T, n=n,
                                    layers=layers, iters=iters, restarts=restarts,
                                    seed=seed, lr=lr)
    return dict(p=p, cvar=float(cv), entropy_nats=float(H),
                entropy_bits=float(H / np.log(2.0)),
                wall=float(time.perf_counter() - t0),
                evals=int(iters * restarts * (1 + 2 * n * layers)), circ=circ)


def readouts(E, p, cfgs=None):
    """Point and distribution readouts from the optimised state.

    argmax_p      the most probable basis state (the VQE's point answer)
    best_in_supp  the lowest-energy state with p above 1/(10 * dim) (what a sampler sees)
    p_weighted    the whole distribution, for a weighted coordinate readout
    """
    E = np.asarray(E, float); p = np.asarray(p, float)
    dim = len(p)
    top = int(np.argmax(p))
    supp = np.where(p > 1.0 / (10.0 * dim))[0]
    bis = int(supp[np.argmin(E[supp])]) if len(supp) else top
    return dict(argmax_state=top, argmax_energy=float(E[top]),
                best_in_support_state=bis, best_in_support_energy=float(E[bis]),
                support_size=int(len(supp)),
                p_top=float(p[top]),
                mean_energy=float((p * E).sum()))


def sample_states(p, shots, rng):
    return rng.choice(len(p), size=int(shots), p=np.asarray(p, float) / np.sum(p))


# -------------------------------------------------------------------- audits
def pennylane_probs_match(n=8, layers=3, seed=0):
    """`lightning.qubit` vs the exact statevector on the same angles."""
    Q = _q()
    rng = np.random.default_rng(seed)
    th = rng.normal(0, 0.6, n * layers)
    circ = Q.build_global_circuit(n, layers, ring=True, device="lightning.qubit")
    pl = np.asarray(circ(th), float).ravel()
    sv = Q.StatevectorCircuit(n, layers, ring=True).probs(th)
    return dict(n=n, layers=layers, max_abs_diff=float(np.abs(pl - sv).max()),
                total_variation=float(0.5 * np.abs(pl - sv).sum()))


def vqe_gradient_audit(E, n, layers=3, alpha=0.2, shots=8192, seed=0):
    """cos(estimator, exact parameter-shift gradient) for both baselines, on THIS E."""
    Q = _q()
    rng = np.random.default_rng(seed)
    circ = Q.StatevectorCircuit(n, layers)
    th = rng.normal(0, 0.6, circ.n_params())
    g_ex = Q.grad_cvar_paramshift(circ, th, E, alpha)
    g_fd = Q.grad_cvar_fd(circ, th, E, alpha)
    out = dict(paramshift_vs_fd_cos=float(g_ex @ g_fd / (np.linalg.norm(g_ex) * np.linalg.norm(g_fd))))
    for bl in ("const", "tail"):
        g, _ = Q.grad_cvar_score(circ, th, E, alpha, shots=shots,
                                 rng=np.random.default_rng(seed + 1), baseline=bl)
        out[f"score_{bl}_cos"] = float(g @ g_ex / (np.linalg.norm(g) * np.linalg.norm(g_ex)))
        out[f"score_{bl}_normratio"] = float(np.linalg.norm(g) / np.linalg.norm(g_ex))
    return out


# ------------------------------------------------- shared-budget device-realistic arm
_HAM_CLS = {}


def qubo_hamiltonian(E, n_qubits, eval_budget=None, cache_limit=1 << 22):
    """A `budget.BudgetedEnergyModel` wrapped around a VQ energy vector.

    `run_global_cvar_vqe` and the classical searches below are handed the SAME object, so
    they compete under the project's own accounting: one charge per UNIQUE bitstring whose
    energy had to be computed, cached thereafter.  That is the sprint's shared-budget
    protocol, not one invented here.
    """
    if "cls" not in _HAM_CLS:
        from budget import BudgetedEnergyModel

        class _QuboHam(BudgetedEnergyModel):
            def __init__(self, E, n, eval_budget, cache_limit):
                self._E = np.asarray(E, float)
                self._n = int(n)
                self._init_budget(cache_limit, eval_budget)

            @property
            def n_qubits(self):
                return self._n

            def energy(self, bitstring):
                bs = str(bitstring)
                if bs in self._cache:
                    return self._cache[bs]
                self._charge()
                e = float(self._E[int(bs, 2)])
                if len(self._cache) < self._cache_limit:
                    self._cache[bs] = e
                return e

        _HAM_CLS["cls"] = _QuboHam
    return _HAM_CLS["cls"](E, n_qubits, eval_budget, cache_limit)


def _stall_cap(ham):
    """Iterations without a fresh charge after which a search is declared saturated.

    A cached bitstring costs nothing, so a local search that has enumerated its whole
    reachable neighbourhood would loop forever against a budget it can no longer spend.
    """
    return 200 * max(ham.n_qubits, 1)


def budgeted_random(ham, seed=0):
    """Uniform random search under the shared budget."""
    from budget import BudgetExhausted
    rng = np.random.default_rng(seed)
    fmt = f"0{ham.n_qubits}b"
    best, bb = np.inf, None
    t0 = time.perf_counter()
    stall, cap = 0, _stall_cap(ham)
    try:
        while stall < cap:
            spent = ham.n_energy_evaluations
            bs = format(int(rng.integers(0, 1 << ham.n_qubits)), fmt)
            e = ham.energy(bs)
            stall = 0 if ham.n_energy_evaluations > spent else stall + 1
            if e < best:
                best, bb = e, bs
    except BudgetExhausted:
        pass
    return dict(best_energy=float(best), best_bitstring=bb,
                n_energy_evaluations=int(ham.n_energy_evaluations),
                wall=float(time.perf_counter() - t0))


def budgeted_anneal(ham, seed=0, T0=None, T1=None):
    """Single-bit-flip simulated annealing under the shared budget."""
    from budget import BudgetExhausted
    rng = np.random.default_rng(seed)
    nq = ham.n_qubits
    t0 = time.perf_counter()
    budget = ham.budget_remaining
    best, bb = np.inf, None
    stall, cap = 0, _stall_cap(ham)
    try:
        x = rng.integers(0, 2, nq)
        cur = ham.energy("".join(map(str, x)))
        best, bb = cur, "".join(map(str, x))
        probe = [ham.energy("".join(map(str, rng.integers(0, 2, nq)))) for _ in range(32)]
        t_hi = (float(np.std(probe)) + 1e-9) if T0 is None else T0
        t_lo = t_hi * 1e-3 if T1 is None else T1
        while stall < cap:
            spent = ham.n_energy_evaluations
            frac = min(1.0, ham.n_energy_evaluations / max(budget, 1.0))
            T = t_hi * (t_lo / t_hi) ** frac
            y = x.copy()
            y[int(rng.integers(0, nq))] ^= 1
            bs = "".join(map(str, y))
            e = ham.energy(bs)
            stall = 0 if ham.n_energy_evaluations > spent else stall + 1
            if e < cur or rng.random() < math.exp(-(e - cur) / max(T, 1e-12)):
                x, cur = y, e
            if cur < best:
                best, bb = cur, bs
    except BudgetExhausted:
        pass
    return dict(best_energy=float(best), best_bitstring=bb,
                n_energy_evaluations=int(ham.n_energy_evaluations),
                wall=float(time.perf_counter() - t0))


def budgeted_greedy_ls(ham, seed=0):
    """Repeated 1-flip steepest descent from random starts, under the shared budget."""
    from budget import BudgetExhausted
    rng = np.random.default_rng(seed)
    nq = ham.n_qubits
    t0 = time.perf_counter()
    best, bb = np.inf, None
    stall, cap = 0, _stall_cap(ham)
    try:
        while stall < cap:
            spent0 = ham.n_energy_evaluations
            x = rng.integers(0, 2, nq)
            cur = ham.energy("".join(map(str, x)))
            if cur < best:
                best, bb = cur, "".join(map(str, x))
            improved = True
            while improved:
                improved = False
                for i in range(nq):
                    y = x.copy(); y[i] ^= 1
                    bs = "".join(map(str, y))
                    e = ham.energy(bs)
                    if e < cur - 1e-12:
                        x, cur, improved = y, e, True
                        if cur < best:
                            best, bb = cur, bs
            stall = 0 if ham.n_energy_evaluations > spent0 else stall + 1
    except BudgetExhausted:
        pass
    return dict(best_energy=float(best), best_bitstring=bb,
                n_energy_evaluations=int(ham.n_energy_evaluations),
                wall=float(time.perf_counter() - t0))
