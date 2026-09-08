"""SPRINT 14 / VQE -- a budget-HONEST CVaR-VQE, and the accounting that makes it honest.

THE ACCOUNTING RULE, and why it is the whole experiment
------------------------------------------------------
An objective evaluation is one read of the energy of one configuration.  Simulating the
circuit is free (it touches no energy); sampling is free; the parameter-shift
log-derivative of the circuit is free.  What costs is looking up ``E(x)``.

This rule is not a convenience, it is what makes the comparison legitimate:

  * `core.quantum.run_cvar_vqe` optimises the EXACT parameter-shift gradient, which needs
    ``cvar_exact(E, p_theta, alpha)`` -- and that reads the ENTIRE energy array, all
    ``4^n`` entries, on EVERY iteration.  An 80-iteration run of it therefore consumes
    ``80 x 262,144 = 21 million`` objective evaluations to search a space that exhaustive
    enumeration solves exactly in 262,144.  It is a beautiful diagnostic and it is not a
    search algorithm.  It is included below as `arm="vqe_exactgrad"` and reported with that
    evaluation count attached, never at nominal parity.
  * The deployable estimator samples ``shots`` bitstrings, reads ``E`` at those and nowhere
    else, and forms the score-function gradient.  ``budget = shots x iters``.  That is the
    arm that can be put next to simulated annealing at matched budget.

READOUT.  Two different structures come out of a VQE and they are not the same object, so
both are reported always:
  ``best_seen``  -- the lowest-energy configuration ever sampled.  This is what the shipped
                    pipeline returns, and it is a best-of-N selector: it makes the VQE a
                    (bad) random search dressed in a circuit unless the distribution
                    actually concentrates.
  ``mode``       -- the argmax of the final distribution.  This is what the variational
                    state actually learned, and it is the only readout for which "the VQE
                    did something" is a meaningful claim.

INITIALISATION.  `init_rmsd` is recorded on every run.  An initialisation that already sits
at the answer is not a VQE result, and the only way to see that is to print it next to the
outcome.
"""
from __future__ import annotations

import numpy as np

from core import quantum as Q
from s14 import vqe_lib as V


# ------------------------------------------------------------------ ansatz factory
def make_ansatz(kind: str, n_qubits: int):
    """Every ansatz here is an exactly-simulated genuine circuit -- no truncation."""
    if kind == "onelayer":
        return Q.OneLayerAnsatz(n_qubits, ring=True)
    if kind.startswith("mps"):
        # "mps3" / "mps3f" / "mps3n" -> layers 3, +final RY, entangler none
        body = kind[3:]
        layers = int(body[0]) if body and body[0].isdigit() else 2
        rest = body[1:] if body and body[0].isdigit() else body
        return Q.MPSAnsatz(n_qubits, layers=layers,
                           final_ry=("f" in rest),
                           entangler="none" if "n" in rest else "cnot")
    raise ValueError(f"unknown ansatz {kind!r}")


def init_theta(ansatz, rng, scale=0.8, mode="random", prior=None, n_res=None,
               bits_per_res=None):
    """Initial angles.  `mode="prior"` warm-starts from a per-residue distribution.

    The warm start is a genuine, native-free initialisation when `prior` is the leakage-safe
    empirical prior; it is an ORACLE initialisation when `prior` came from `ORACLE_prior`.
    Callers must label which, and `run` records `init_rmsd` either way.
    """
    P = ansatz.n_params()
    if mode == "random":
        return rng.uniform(0.0, np.pi, P) if scale is None else \
            np.pi / 2 + rng.normal(0.0, scale, P)
    if mode == "zero":
        return np.full(P, np.pi / 2) + rng.normal(0, 1e-3, P)
    if mode == "prior":
        # target per-qubit marginals from the per-residue state distribution, then invert
        # the RY -> Bernoulli map on the LAST rotation block (the free-marginal one).
        m = _qubit_marginals(prior, bits_per_res)
        th = np.pi / 2 + rng.normal(0, 0.05, P)
        nb = P // len(m)
        th = th.reshape(nb, len(m))
        th[-1] = 2.0 * np.arcsin(np.sqrt(np.clip(m, 1e-6, 1 - 1e-6)))
        th[:-1] = 1e-3 * rng.standard_normal(th[:-1].shape)
        return th.reshape(-1)
    raise ValueError(mode)


def _qubit_marginals(P: np.ndarray, bits_per_res: int) -> np.ndarray:
    """(n, k) state distribution -> (n * bits_per_res,) P(qubit = 1), MSB first."""
    P = np.asarray(P, float)
    n, k = P.shape
    out = np.empty(n * bits_per_res)
    idx = np.arange(k)
    for i in range(n):
        for b in range(bits_per_res):
            hit = (idx >> (bits_per_res - 1 - b)) & 1
            out[i * bits_per_res + b] = float(P[i][hit == 1].sum() / P[i].sum())
    return out


# --------------------------------------------------------------------- the driver
def run(E: np.ndarray, n_qubits: int, alpha, budget: int, shots: int = 512,
        ansatz: str = "mps2f", seed: int = 0, lr: float = 0.15,
        baseline: str = "const", rmsd: np.ndarray = None,
        init: str = "random", init_scale: float = 0.8, prior=None,
        bits_per_res: int = 2, exact_dist: bool = True,
        alpha_anneal: bool = False, keep_seen: bool = False) -> dict:
    """Sampled CVaR-VQE with a hard objective-evaluation budget.

    `alpha` may be a float or, with `alpha_anneal`, the (a0, a1) pair fed to
    `core.quantum.alpha_schedule`.  Returns both readouts, both axes, the trajectory in
    parameter space, and the final distribution's entropy and near-native mass.
    """
    rng = np.random.default_rng(seed)
    an = make_ansatz(ansatz, n_qubits)
    th = init_theta(an, rng, init_scale, init, prior, bits_per_res=bits_per_res)
    th0 = th.copy()
    opt = Q.Adam(an.n_params(), lr=lr)
    c = V.Counter(E, budget)
    iters = max(1, budget // shots)
    traj = [th.copy()]
    cvals, gnorms = [], []
    a0, a1 = (alpha if isinstance(alpha, (tuple, list)) else (alpha, alpha))

    # the initialisation's own structure, before any optimisation
    b0 = an.sample(th0, shots, np.random.default_rng(seed + 991))
    i0 = _bits_to_index(b0)
    init_rmsd = float(rmsd[i0].mean()) if rmsd is not None else float("nan")
    init_best_rmsd = float(rmsd[i0].min()) if rmsd is not None else float("nan")

    for t in range(iters):
        if c.left < shots:
            break
        a = Q.alpha_schedule(t / max(1, iters - 1), a0, a1) if alpha_anneal else a0
        bits = an.sample(th, shots, rng)
        idx = _bits_to_index(bits)
        e = c(idx)
        if e.size < shots:
            break
        g, val = Q.cvar_gradient(an, th, bits, e, a, baseline=baseline)
        th = opt.step(th, g)
        traj.append(th.copy())
        cvals.append(float(val))
        gnorms.append(float(np.linalg.norm(g)))

    out = {
        "arm": f"vqe_{ansatz}_a{a0}" + ("_anneal" if alpha_anneal else ""),
        "ansatz": ansatz, "alpha": a0, "alpha_end": a1, "anneal": bool(alpha_anneal),
        "shots": shots, "iters": len(cvals), "budget": budget, "evals": int(c.used),
        "seed": seed, "lr": lr, "baseline": baseline, "init": init,
        "init_scale": init_scale, "n_qubits": n_qubits,
        "n_params": int(an.n_params()),
        "cvar_trajectory": cvals, "grad_norms": gnorms,
        "grad_var": float(np.var(gnorms)) if gnorms else float("nan"),
        "theta0": th0.tolist(), "theta": th.tolist(),
        "param_path_length": float(sum(np.linalg.norm(traj[i + 1] - traj[i])
                                       for i in range(len(traj) - 1))),
        "param_displacement": float(np.linalg.norm(th - th0)),
        "best_e": float(c.best_e), "best_i": int(c.best_i),
        "init_rmsd_mean": init_rmsd, "init_rmsd_best": init_best_rmsd,
    }
    # -- the two readouts, and the distribution ---------------------------
    if exact_dist:
        p = an.probs(th)
        p = np.maximum(p, 0.0); p = p / p.sum()
        out["mode_i"] = int(np.argmax(p))
        out["entropy_bits"] = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
        out["max_prob"] = float(p.max())
        out["eff_support"] = float(2.0 ** out["entropy_bits"])
        out["mean_energy"] = float(p @ np.asarray(E, float))
        if rmsd is not None:
            out["mode_rmsd"] = float(rmsd[out["mode_i"]])
            out["mean_rmsd_under_p"] = float(p @ rmsd)
            for thr in (2.0, 2.5, 3.0):
                out[f"pmass_below_{thr}"] = float(p[rmsd < thr].sum())
    if rmsd is not None:
        seen = c.all_seen()
        out["rmsd_returned"] = float(rmsd[c.best_i])
        out["rmsd_mean_seen"] = float(rmsd[seen].mean())
        out["rmsd_best_seen"] = float(rmsd[seen].min())
        out["n_distinct"] = int(np.unique(seen).size)
        out["diversity"] = float(np.unique(seen).size / max(1, seen.size))
    if keep_seen:
        out["seen"] = c.all_seen()
    return out


def _bits_to_index(bits: np.ndarray) -> np.ndarray:
    b = np.asarray(bits, np.int64)
    nq = b.shape[1]
    return b @ (1 << np.arange(nq - 1, -1, -1))


def run_exactgrad(E, n_qubits, alpha, iters=60, layers=2, seed=0, lr=0.15,
                  rmsd=None, T=0.0):
    """The EXACT-gradient arm.  Reads all 4^n energies per iteration -- reported as such.

    This is `core.quantum.free_energy` / `run_cvar_vqe`'s machinery, kept because it is the
    noise-free upper bound on what the ansatz can do, but its evaluation count is
    ``iters * 4^n`` and it is never placed at nominal budget parity with a sampled arm.
    """
    circ = Q.StatevectorCircuit(n_qubits, layers)
    rng = np.random.default_rng(seed)
    th = rng.normal(0.0, 0.8, circ.n_params())
    opt = Q.Adam(circ.n_params(), lr=lr)
    cv = []
    for _ in range(iters):
        f, g, p, v, H = Q.free_energy(circ, th, E, alpha, T)
        th = opt.step(th, g)
        cv.append(float(v))
    _, _, p, v, H = Q.free_energy(circ, th, E, alpha, T)
    out = {"arm": f"vqe_exactgrad_a{alpha}", "alpha": alpha, "layers": layers,
           "iters": iters, "seed": seed, "T": T,
           "evals": int(iters * len(E)), "cvar": float(v), "entropy_nats": float(H),
           "cvar_trajectory": cv, "mode_i": int(np.argmax(p)),
           "max_prob": float(p.max()), "mean_energy": float(p @ np.asarray(E, float))}
    if rmsd is not None:
        out["mode_rmsd"] = float(rmsd[out["mode_i"]])
        out["mean_rmsd_under_p"] = float(p @ rmsd)
        for thr in (2.0, 2.5, 3.0):
            out[f"pmass_below_{thr}"] = float(p[rmsd < thr].sum())
    return out
