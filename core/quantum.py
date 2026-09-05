"""The VQE/CVaR subsystem: one module.

CONSOLIDATED FROM
=================
``qansatz.py`` (324)     CVaR, the RY/CNOT ansatz families, Adam, the CVaR gradient
``vqe.py`` (570)         the PennyLane global CVaR-VQE, SPSA, the shared-budget driver
``foldvqe.py`` (683)     warm start, reservoir, refinement, basin selection, ``fold``
``objective.py`` (150)   ``FoldObjective`` -- the batched composite folding score
``hamiltonian.py`` (134) ``FoldingHamiltonian`` -- the knowledge-based bitstring energy

Nothing here is a re-export: the five modules are still on disk and
``tests/test_quantum.py`` uses them as the equivalence oracle for every routine that was
rewritten.  They are the reference, not the implementation.


WHAT IS GENUINE HERE, AND WHY THAT MATTERS
==========================================
This is a real variational quantum eigensolver, not enumeration wearing its name:

* a **parameterised circuit** -- ``layers x (RY on every wire, CNOT chain, optional ring
  closure)`` on ``|0...0>`` -- in three simulators that are each exact for their regime:
  PennyLane ``lightning.qubit`` (``build_global_circuit``, n <= 30), a batched statevector
  (``StatevectorCircuit``, n <= ~20), and an exact matrix-product state (``MPSAnsatz``,
  no qubit ceiling; the chain topology caps the bond dimension at ``2 ** layers``, so
  *nothing is truncated* -- there is no approximation to trade off).
* **real expectation values**: probabilities come out of the simulated state, bitstrings
  are drawn from that distribution, and the objective is evaluated on what was drawn.
* **genuine CVaR** -- the conditional value at risk of the objective distribution over the
  measured bitstrings, i.e. the mean of its lower ``alpha`` tail, not the mean.

``all_bitstrings`` / ``cvar_gradient_exact`` / ``grad_cvar_fd`` enumerate the register.
They are **verification instruments**, they are named so, and they are never the search.

A measured fact worth knowing that is NOT licence to weaken any of this: at these problem
sizes the VQE ties uniform random sampling and loses to annealing (s9 ``stage_search``).
The component is mandatory and stays genuine; the honest number is reported rather than
engineered away.


THE CVaR GRADIENT, AND THE DEFECT THAT MUST NOT COME BACK
=========================================================
The objective is written down and differentiated:

    CVaR_alpha(theta) = max_t  t - (1/alpha) E_{x~p_theta} (t - E(x))_+

At the optimal ``t`` (the alpha-quantile ``q``) the envelope theorem kills the ``t``
dependence and leaves the score-function form

    grad CVaR = -(1/alpha) E_p[ (q - E(x))_+ grad log p_theta(x) ]
              =  E_p[ f(x) grad log p_theta(x) ],   f(x) = -(q - E(x))_+ / alpha .

A baseline ``b`` may be subtracted from ``f`` **if and only if it is constant in x**,
because the correction term is ``b * E_p[grad log p] = 0``.  The shipped
``qansatz.cvar_gradient`` subtracts the TAIL MEAN from the tail entries and leaves the
non-tail entries at zero -- that is ``b(x) = m * 1[x in tail]``, a *function of x*, and
the identity it relies on does not hold on a data-dependent subset.  It is a bias, not
extra variance.

Priced against an exact classical reference (36 checks at 10 qubits, i.e. 1024 amplitudes
enumerated exactly -- ``s9/refine.stage_grad``):

    estimator                                    cos with exact grad   |g| / |g_exact|
    -------------------------------------------  --------------------  ---------------
    exact expectation, CONSTANT baseline               +1.000000            1.000
    exact expectation, shipped TAIL-ONLY baseline      +0.655634            0.758
    sampled, constant baseline                         +0.994               ~1
    parameter-shift vs exact finite differences        +1.000000            1.000

The exact-expectation rows carry *zero sampling noise*, which is what proves the tail-only
number is bias.  ``cvar_gradient`` here defaults to ``baseline="const"`` -- the fix -- and
``baseline="tail"`` reproduces the shipped defect verbatim (to 5.6e-17), so the comparison
stays runnable and so ``tests/test_quantum.py`` can assert the defect is gone by
*measuring* it rather than by reading the source.  That regression test is the single most
important thing in the test file.


WHAT WAS MADE FAST, AND WHAT WAS NOT TOUCHED
============================================
Shots, iterations, restarts, the ansatz family, the Hamiltonian, the objective and every
seed path are unchanged.  The speed came from removing Python/framework overhead, which
the profile said was ~100% of the cost at these register sizes (a 60-site MPS contraction
moves 4x4 matrices; the flops are not the problem, the 1,000 framework calls are):

* ``cvar`` selects the tail with ``np.argpartition`` (O(N)) instead of ``np.argsort``
  (O(N log N)).  For a tail fraction alpha the two are *mathematically identical*, ties
  included -- see ``tail_indices`` for the proof and the tests for the assertion.  Both
  paths ship: below a measured crossover of ~450 samples the sort's constant factor wins,
  and above it the selection is 5.9x at 2,048 samples and 12-25x beyond.
* ``MPSAnsatz`` keeps every site in ONE padded ``(n, chi, 2, chi)`` array instead of a
  Python list of ragged tensors, so a rotation layer is one ``einsum`` and an entangling
  layer six array ops -- rather than ``2n`` and ``7n``.  The CNOT chain is applied as the
  exact bond-dimension-2 MPO it already was.
* the MPS norm and amplitudes are contracted by a **log-depth matrix-chain reduction**
  (associativity), so 60 sites take 6 batched matmuls instead of 60 sequential ones.  This
  is the one place the reduction ORDER changed; measured against the shipped
  left-to-right contraction the difference is ~1e-13 on log p and ~1e-14 relative on the
  gradient (``t_mps_logp_matches_legacy``, ``t_mps_grad_logp_matches_legacy``).
* right environments for sampling come from a Hillis-Steele suffix scan over the transfer
  matrices: 6 batched matmuls instead of 60 sequential ones, and the autoregressive loop
  that follows (which is genuinely sequential -- each bit conditions the next) trades two
  ``np.einsum`` parses per site for two matmuls on contiguous blocks.  Bit-identical bits.
* ``StatevectorCircuit`` evaluates **all 2P shifted parameter vectors in one batched
  simulation** instead of re-running the whole stack per shifted parameter, and the
  entangling layer's basis permutation -- a topology invariant -- is composed once in
  ``__init__`` rather than rebuilt from ``np.arange`` on every gate of every call.  The
  arithmetic is elementwise identical, so parameter-shift results are BIT-identical.
* the PennyLane path builds its tape once and rebinds parameters instead of letting the
  QNode reconstruct and re-transform the tape on every objective call.  Bit-identical.
* ``fold`` can run its restarts concurrently -- they are genuinely independent, and their
  reservoir offers are recorded per restart and replayed in restart order, so the pool is
  bit-identical to the sequential run.  It is OFF by default, because it was measured and
  it does not pay: at 60 qubits and 384 shots the inner loop moves kilobyte arrays, so the
  cost is interpreter dispatch rather than the BLAS calls threads could overlap.  The
  numbers are in `fold`'s docstring.  A negative result, reported rather than buried.

Every one of those claims is a test in ``tests/test_quantum.py``.
"""
import ctypes
import itertools
import math
import os
import threading
import time
import warnings
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from budget import (BudgetExhausted, BudgetedEnergyModel, resolve_maxiter,
                    check_optimizer_budget)

__all__ = [
    "cvar", "cvar_from_samples", "cvar_from_distribution", "cvar_from_probs",
    "alpha_schedule", "tail_indices", "BASELINES",
    "OneLayerAnsatz", "MPSAnsatz", "Adam", "make_ansatz",
    "cvar_gradient", "cvar_gradient_exact", "grad_logp_weighted", "all_bitstrings",
    "StatevectorCircuit", "cvar_exact", "grad_cvar_paramshift", "grad_cvar_fd",
    "grad_cvar_score", "free_energy", "run_cvar_vqe", "entropy_onelayer",
    "build_global_circuit", "n_parameters", "BestSeenTracker",
    "run_global_cvar_vqe", "MIN_USEFUL_SPSA_ITERS",
    "FoldObjective", "FoldingHamiltonian",
    "marginals_to_angles", "state_marginals_to_bits", "Reservoir", "refine",
    "snap_to_states", "basin_average", "basin_select", "build_distogram",
    "build_target", "seed_states", "run_vqe", "fold", "fold_multi",
    "fold_iterated", "consensus_distances",
    "memory_load", "mem_gate", "limit_threads",
]


# ============================================================== resource hygiene
class _MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


def memory_load() -> float:
    """System memory load in percent, via ``GlobalMemoryStatusEx``.

    Deliberately NOT ``Get-CimInstance``: under load that costs minutes per call on this
    box, which makes the gate more expensive than the job it gates.  This is one kernel32
    call and costs microseconds.  Returns ``nan`` off Windows, which reads as "proceed".
    """
    try:
        m = _MEMORYSTATUSEX()
        m.dwLength = ctypes.sizeof(_MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return float(m.dwMemoryLoad)
    except Exception:
        return float("nan")


def mem_gate(tag: str = "", limit: float = 92.0, wait: float = 5.0,
             tries: int = 60) -> float:
    """Block while system memory load is above `limit`. Returns the load it last saw."""
    load = memory_load()
    for _ in range(int(tries)):
        if not (load == load) or load <= limit:
            return load
        time.sleep(wait)
        load = memory_load()
    warnings.warn(f"mem_gate({tag}): still at {load:.0f}% after waiting", RuntimeWarning)
    return load


def limit_threads(n: int = 2) -> None:
    """Pin torch (and BLAS where the environment still allows it) to `n` threads.

    Sibling jobs share this box, but this is a performance fix as much as a courtesy: a
    torch default of one thread per core turns a 4x4 matmul into a thread-pool round trip
    and makes the whole MPS stack slower.
    """
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(v, str(int(n)))
    try:
        import torch
        torch.set_num_threads(int(n))
    except Exception:
        pass


# ========================================================================== CVaR
#: Baselines accepted by `cvar_gradient`. "const" is the fix; "tail" is the recorded
#: defect, kept ONLY so the audit and the regression test can measure it.
BASELINES = ("const", "tail", "none")

#: Below this many samples a full ``np.argsort`` beats the O(N) selection, because the
#: selection needs four passes (partition, ``e < q``, a count, and the tie fix-up) against
#: introsort's one very well optimised one. MEASURED on this box, alpha=0.1, in
#: microseconds per call:
#:
#:     N          64    128    256    384    512   1024   2048   4096  16384  1048576
#:     sort     3.98   5.25   7.61   9.09  13.93  26.01  86.18  242.5   1295   149565
#:     select   7.18   8.51  10.02  10.45  10.90  12.25  14.56   21.3     50     8599
#:
#: The crossover is ~450, so both paths are kept and the cheaper one is taken. They return
#: the identical mask, quantile and value -- ``t_sort_and_partition_*`` asserts that on
#: both sides of the cutoff and under heavy ties -- so this is a cost decision only.
#: `foldvqe` runs 384 shots (sort), the global driver runs 2048 (select, 5.9x) and the
#: gradient audits run 4096 (select, 11.4x).
CVAR_SORT_CUTOFF = 448


def tail_indices(energies: np.ndarray, alpha: float) -> Tuple[int, float, np.ndarray]:
    """``(k, quantile, tail mask)`` for the lower `alpha` tail, in O(N).

    THE SORT/PARTITION EQUIVALENCE
    ------------------------------
    The shipped routine took ``order = np.argsort(e, kind="stable")`` and kept
    ``order[:k]``.  CVaR needs a SELECTION, not an ordering, and ``np.argpartition`` does
    selection in O(N) where a sort is O(N log N).  The two give the identical mask, ties
    included, and here is why: a stable sort places equal elements in index order, so
    ``order[:k]`` is exactly

        {x : e[x] < q}  union  {the (k - #below) lowest INDICES with e[x] == q}

    with ``q`` the k-th smallest value.  That is what this function constructs, so the
    mask, the quantile, and therefore the CVaR value and every gradient weight, are the
    same object rather than an approximation of it.
    """
    e = np.asarray(energies, float)
    if e.size == 0:
        raise ValueError("cvar received no samples")
    if not (0.0 < alpha <= 1.0):
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    n = e.size
    k = max(1, int(math.ceil(alpha * n)))
    if k >= n:
        return n, float(e.max()), np.ones(n, bool)
    if n <= CVAR_SORT_CUTOFF:
        order = np.argsort(e, kind="stable")       # cheaper at this size; same answer
        tail = np.zeros(n, bool)
        tail[order[:k]] = True
        return k, float(e[order[k - 1]]), tail
    part = np.argpartition(e, k - 1)
    q = float(e[part[k - 1]])
    tail = e < q
    below = int(np.count_nonzero(tail))
    need = k - below
    # `need` is 1 whenever the quantile value is unique, which is the common case for
    # continuous energies: exactly k-1 samples are strictly below the k-th smallest.
    if need == 1:
        tail[int(np.argmax(e == q))] = True        # the LOWEST index among the ties
    elif need > 1:
        ties = np.flatnonzero(e == q)              # already in ascending index order
        tail[ties[:need]] = True
    return k, q, tail


def cvar(energies: np.ndarray, alpha: float) -> Tuple[float, float, np.ndarray]:
    """``(cvar, quantile, tail mask)`` for the LOWER alpha tail of `energies`.

    Drop-in for ``qansatz.cvar``: identical outputs, O(N) instead of O(N log N).
    """
    e = np.asarray(energies, float)
    _, q, tail = tail_indices(e, alpha)
    return float(e[tail].mean()), q, tail


def cvar_from_samples(energies: Sequence[float], alpha: float) -> float:
    """Sampled CVaR: the mean of the lowest ``ceil(alpha * N)`` energies.

    This is the objective the global SPSA driver evaluates, so it is worth two details.

    First, order: the shipped ``vqe.cvar_from_samples`` summed the tail in ASCENDING order
    (``np.sort(e)[:k].mean()``), and floating-point addition is not associative, so summing
    the same k values in index order lands one ulp away -- over thousands of noisy SPSA
    evaluations, enough to fork a trajectory.  So the k values are sorted before the mean,
    which makes the result BIT-identical: both arrays hold the same multiset in ascending
    order, so they are equal element by element.

    Second, this one does NOT go through `tail_indices`.  It needs the tail's VALUES, not
    its positions, so ``np.partition`` beats ``np.argpartition`` -- no index array, no
    boolean mask, no gather -- and the ordering rule for ties is then irrelevant, because
    tied values are indistinguishable in a sum.  ``O(N + k log k)``.
    """
    e = np.asarray(energies, dtype=float)
    if e.size == 0:
        raise ValueError("cvar received no samples")
    if not (0.0 < alpha <= 1.0):
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    k = max(1, int(math.ceil(alpha * e.size)))
    if k >= e.size or e.size <= CVAR_SORT_CUTOFF:
        return float(np.sort(e)[:k].mean())
    return float(np.sort(np.partition(e, k - 1)[:k]).mean())


def cvar_from_distribution(energies: np.ndarray, probs: np.ndarray,
                           alpha: float) -> float:
    """CVaR of an explicit discrete distribution, boundary mass split fractionally."""
    energies = np.asarray(energies, dtype=float)
    probs = np.clip(np.asarray(probs, dtype=float), 0.0, None)
    tot = probs.sum()
    if tot <= 0:
        return float(np.min(energies))
    probs = probs / tot
    order = np.argsort(energies)
    acc = esum = 0.0
    for k in order:
        p = probs[k]
        if acc + p < alpha:
            esum += p * energies[k]
            acc += p
        else:
            esum += (alpha - acc) * energies[k]
            acc = alpha
            break
    return esum / acc if acc > 0 else float(energies[order[0]])


def cvar_from_probs(energies: np.ndarray, probs: np.ndarray, alpha: float
                    ) -> Tuple[float, float, np.ndarray]:
    """EXACT CVaR of an explicit distribution: ``(value, quantile, tail mass)``.

    Vectorised form of `cvar_from_distribution` that also returns the per-state mass the
    tail used, which is what makes the value a differentiable function of `probs` and
    therefore what the finite-difference tests differentiate.
    """
    e = np.asarray(energies, float)
    p = np.asarray(probs, float)
    if not (0.0 < alpha <= 1.0):
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    tot = p.sum()
    if tot <= 0:
        raise ValueError("probabilities sum to zero")
    p = p / tot
    order = np.argsort(e, kind="stable")
    cum = np.cumsum(p[order])
    take = np.clip(alpha - (cum - p[order]), 0.0, p[order])
    mass = np.zeros_like(p)
    mass[order] = take
    m = int(np.searchsorted(cum, alpha, side="left"))
    q = float(e[order[min(m, len(order) - 1)]])
    return float((mass * e).sum() / alpha), q, mass


def alpha_schedule(progress: float, a0: float = 0.5, a1: float = 0.05) -> float:
    """Geometric anneal of the CVaR level from `a0` to `a1` over a run.

    Early on a narrow tail concentrates the distribution before the landscape has been
    sampled -- the premature-concentration failure seen at ``init_scale=0.25``.  Late on a
    broad tail is just the mean, and the mean is not what we want to minimise.
    """
    p = min(1.0, max(0.0, float(progress)))
    return float(a0 * (a1 / a0) ** p)


# ======================================================================= ansatz
class OneLayerAnsatz:
    """RY on every wire, then a CNOT chain (and optional ring closure).

    The prepared distribution is exactly: draw independent bits ``y_q ~
    Bern(sin^2(theta_q / 2))``, then output the cumulative XOR.  Sampling is a vectorised
    prefix-XOR and both the log-probability and its gradient are closed form, so nothing
    on this path is a simulation approximation -- it is the exact output distribution of
    the circuit, obtained analytically instead of by contracting a state.

    Unchanged from ``qansatz.OneLayerAnsatz`` apart from hoisting the per-call constants:
    at 96 shots it already costs 25-40 microseconds and is not where the time goes.
    """

    layers = 1

    def __init__(self, n_qubits: int, ring: bool = True):
        self.n = int(n_qubits)
        self.ring = bool(ring) and self.n > 2

    def n_params(self) -> int:
        return self.n

    def _p1(self, theta: np.ndarray) -> np.ndarray:
        return np.sin(np.asarray(theta, float) / 2.0) ** 2

    def sample(self, theta: np.ndarray, shots: int, rng) -> np.ndarray:
        y = (rng.random((int(shots), self.n)) < self._p1(theta)).astype(np.uint8)
        b = np.bitwise_xor.accumulate(y, axis=1)
        if self.ring:
            b = b.copy()
            b[:, 0] ^= b[:, -1]
        return b

    def _pre_bits(self, bits: np.ndarray) -> np.ndarray:
        """Invert the CNOT network: recover the independent Bernoulli bits."""
        b = np.asarray(bits, np.uint8)
        if self.ring:
            b = b.copy()
            # b0' = b0 ^ b_{n-1}; the chain fixes b_{n-1} from the *outputs*, so undo the
            # ring first using the untouched last output bit.
            b[:, 0] ^= b[:, -1]
        y = b.copy()
        y[:, 1:] ^= b[:, :-1]
        return y

    def logp(self, theta: np.ndarray, bits: np.ndarray) -> np.ndarray:
        s = np.clip(self._p1(theta), 1e-12, 1 - 1e-12)
        y = self._pre_bits(bits)
        return (y * np.log(s) + (1 - y) * np.log1p(-s)).sum(1)

    def grad_logp(self, theta: np.ndarray, bits: np.ndarray) -> np.ndarray:
        """``(B, n_params)`` gradient of log p at each sampled bitstring."""
        th = np.asarray(theta, float)
        s = self._p1(th)
        y = self._pre_bits(bits).astype(float)
        sin = np.sin(th)
        sin = np.where(np.abs(sin) < 1e-6, np.sign(sin + 1e-12) * 1e-6, sin)
        return 2.0 * (y - s[None, :]) / sin[None, :]

    def probs(self, theta: np.ndarray) -> np.ndarray:
        """Exact distribution over all ``2**n`` bitstrings. Verification only."""
        return np.exp(self.logp(theta, all_bitstrings(self.n)))


def entropy_onelayer(ansatz: "OneLayerAnsatz", theta: np.ndarray
                     ) -> Tuple[float, np.ndarray]:
    """Exact entropy (nats) and its gradient for `OneLayerAnsatz`.

    The CNOT chain (and the ring closure) is a BIJECTION on bitstrings, so the output
    distribution's entropy equals that of the independent Bernoulli layer: a closed form,
    not an estimate.  Used by the free-energy objective in `run_vqe_freeenergy`, where
    minimising CVaR alone is degenerate -- the minimiser concentrates on the lowest-energy
    basis states and the readout collapses to the argmin -- and the entropy term is what
    makes the optimum interior.
    """
    th = np.asarray(theta, float)
    s = np.clip(np.sin(th / 2.0) ** 2, 1e-12, 1 - 1e-12)
    H = float(-(s * np.log(s) + (1 - s) * np.log1p(-s)).sum())
    g = np.log((1 - s) / s) * (np.sin(th) / 2.0)
    return H, g


_EYE2 = np.eye(2)


class MPSAnsatz:
    """``layers`` repetitions of (RY on every wire, CNOT chain), simulated EXACTLY.

    WHY THIS IS EXACT AND NOT AN APPROXIMATION
    ------------------------------------------
    A CNOT is ``|0><0| (x) I + |1><1| (x) X``, an exact bond-dimension-2 MPO, so applying
    the chain doubles one bond with no truncation and no SVD.  After ``L`` entangling
    layers every interior bond is exactly ``2 ** L`` (2 or 4 in practice) and the two
    boundary bonds are 1.  Nothing is discarded, so the 30-qubit statevector wall is gone:
    cost is O(n chi^3), linear in the number of qubits, and a 64- or 96-qubit register is
    routine.

    ``final_ry`` appends one more RY layer with no entangler after it, and it is not
    cosmetic.  The plain chain cannot express an arbitrary set of per-qubit marginals at
    all: ``b_q`` is the prefix XOR of independent Bernoulli bits, so ``|1 - 2 P(b_q = 1)|``
    is a running product and must be non-increasing along the wire.  Measured against the
    torsion prior's own marginals on 12 targets, 53% of chain steps violate that ordering
    and the best achievable initialisation is 0.13 off per qubit on average (0.29 at
    worst).  Directly optimising the angles to match a random marginal vector leaves a
    residual of 0.21 for the chain and 0.000 with one trailing RY layer.

    ``entangler="none"`` drops the CNOTs entirely, giving a product distribution -- free
    marginals, no correlations.  It is the control that says whether the entanglement is
    doing anything.

    The ring closure is not applied here: it is a long-range gate for an open MPS and
    would cost a swap network for a correlation the chain already carries at depth >= 2.

    HOW THE REWRITE WORKS (and why it is the same arithmetic)
    ---------------------------------------------------------
    The shipped implementation kept a Python LIST of ragged tensors and touched each site
    individually: ``2n`` framework calls for a rotation layer and ``7n`` for an entangling
    layer, then ``n`` sequential einsums for the right environments and ``6n`` more inside
    the log-probability loop.  At n=60 that is ~1,000 torch operations on 4x4 tensors --
    entirely dispatch overhead.

    Here every site lives in ONE padded array ``A`` of shape ``(n, chi, 2, chi)``:

    * a rotation layer is a single ``einsum("nij,najb->naib", G, A)``;
    * an entangling layer is the MPO applied to every site at once.  Per site the CNOT
      chain does a TARGET expansion (new left index ``a*2 + c_in``, carrying the incoming
      control bit) and then a CONTROL expansion (new right index ``r*2 + c_out``, with
      ``c_out`` pinned to this wire's post-target physical value).  Target-then-control is
      the order the sequential loop realises and it is not symmetric: qubit ``q`` is the
      target of CNOT(q-1,q) before it is the control of CNOT(q,q+1), so the control reads
      the *flipped* value.  Site 0 is never a target (its ``c_in`` is masked to 0) and
      site n-1 is never a control (it keeps its unexpanded right bond, zero-padded).

    The index conventions are the shipped ones exactly -- ``c_out`` fast on the right bond,
    ``c_in`` fast on the left -- so the padded array contracts to the identical state.
    ``t_mps_padded_build_matches_legacy`` checks tensor-by-tensor contraction to 0.0.
    """

    def __init__(self, n_qubits: int, layers: int = 2, final_ry: bool = False,
                 entangler: str = "cnot"):
        import torch
        self.torch = torch
        self.n = int(n_qubits)
        self.layers = int(layers)
        self.final_ry = bool(final_ry)
        self.entangler = entangler
        self.ring = False
        self.nblocks = self.layers + int(self.final_ry)
        self._entangling = (self.entangler != "none")
        # ---- topology invariants, hoisted out of every evaluation -----------
        # Bond dimension entering each block, and the masks the two boundary sites need.
        # These depend only on (n, layers, entangler) -- rebuilding them per call was
        # ~15% of the build and 100% avoidable.
        self._Din: List[int] = []
        self._m0: List[np.ndarray] = []           # site-0 left mask, per entangling layer
        self._sm: Optional[np.ndarray] = None     # last-site selector, shape (n,1,1,1)
        D = 1
        for l in range(self.nblocks):
            self._Din.append(D)
            if l >= self.layers or not self._entangling:
                continue
            m = np.ones((self.n, 2 * D, 1, 1))
            m[0, 1:] = 0.0
            self._m0.append(m)
            D *= 2
        self.chi = D
        if self._m0:
            sm = np.zeros((self.n, 1, 1, 1))
            sm[self.n - 1] = 1.0
            self._sm = sm
        self._t_m0: List = []
        self._t_sm = None
        self._t_eye2 = None
        # numpy build cache, keyed on the parameter bytes: `sample` and any repeated
        # evaluation at one theta reuse the contraction instead of redoing it.
        self._np_cache: Tuple = (None, None, None)

    # -- shape/plumbing ---------------------------------------------------
    def n_params(self) -> int:
        return self.n * self.nblocks

    def _torch_consts(self):
        t = self.torch
        if self._t_eye2 is None:
            self._t_eye2 = t.as_tensor(_EYE2)
            self._t_m0 = [t.as_tensor(m) for m in self._m0]
            self._t_sm = None if self._sm is None else t.as_tensor(self._sm)
        return self._t_eye2, self._t_m0, self._t_sm

    # -- construction -----------------------------------------------------
    def build(self, theta):
        """Padded MPS ``(n, chi, 2, chi)`` for torch (autograd) or numpy input."""
        istorch = not isinstance(theta, np.ndarray)
        t = self.torch if istorch else np
        n, nb = self.n, self.nblocks
        th = theta.reshape(nb, n)
        if istorch:
            eye2, m0s, sm = self._torch_consts()
            A = self.torch.zeros(n, 1, 2, 1, dtype=self.torch.float64)
        else:
            eye2, m0s, sm = _EYE2, self._m0, self._sm
            A = np.zeros((n, 1, 2, 1))
        A[:, 0, 0, 0] = 1.0
        li = 0
        for l in range(nb):
            D = self._Din[l]
            c = t.cos(th[l] / 2.0)
            s = t.sin(th[l] / 2.0)
            G = t.stack([t.stack([c, -s], -1), t.stack([s, c], -1)], -2)   # (n, 2, 2)
            A = t.einsum("nij,najb->naib", G, A)
            if l >= self.layers or not self._entangling:
                continue
            # TARGET expansion: new left index a*2 + c_in, value A[a, s ^ c_in, r].
            if istorch:
                T = self.torch.stack([A, A.flip(2)], 2)
            else:
                T = np.stack([A, A[:, :, ::-1]], 2)
            T = T.reshape(n, 2 * D, 2, D) * m0s[li]        # site 0 has no incoming CNOT
            # CONTROL expansion: new right index r*2 + c_out, with c_out == this wire's
            # (post-target) physical value.
            C = (T[..., None] * eye2[None, None, :, None, :]).reshape(n, 2 * D, 2, 2 * D)
            # site n-1 is never a control: keep its unexpanded right bond, zero-padded so
            # the array stays rectangular.  Only right index 0 is ever non-zero there.
            if istorch:
                Tp = self.torch.cat([T, self.torch.zeros_like(T)], -1)
            else:
                Tp = np.concatenate([T, np.zeros_like(T)], -1)
            A = C * (1.0 - sm) + Tp * sm
            li += 1
        return A

    def _build_cached(self, theta: np.ndarray) -> np.ndarray:
        th = np.ascontiguousarray(np.asarray(theta, float))
        key = th.tobytes()
        if self._np_cache[0] == key:
            return self._np_cache[1]
        A = self.build(th)
        self._np_cache = (key, A, None)
        return A

    # -- contraction ------------------------------------------------------
    @staticmethod
    def _chain_reduce(M, istorch, cat):
        """Product ``M[0] @ M[1] @ ... @ M[k-1]`` by pairwise (log-depth) reduction.

        Matrix multiplication is associative, so this is the same product the shipped
        left-to-right loop computed; only the ORDER of the floating-point reductions moved.
        That is the single arithmetic change in this class and it is measured, not
        asserted: ~1e-13 on log p, ~1e-14 relative on the gradient.
        """
        while M.shape[0] > 1:
            m = M.shape[0]
            h = m // 2
            P = M[:2 * h:2] @ M[1:2 * h:2]
            if m % 2:
                P = cat([P, M[-1:]], 0)
            M = P
        return M[0]

    def norm(self, A):
        """``<psi|psi>``, i.e. the shipped ``R[0][0, 0]``.

        Contracts the transfer matrices
        ``M[q][(a,d),(b,c)] = sum_s A[q][a,s,b] A[q][d,s,c]``
        with the right boundary ``R[n] = E_00`` (the last site's right bond is size 1 in
        the unpadded picture, and index 0 in the padded one) and reads off ``(a,d)=(0,0)``.
        """
        istorch = not isinstance(A, np.ndarray)
        t = self.torch if istorch else np
        n, D = self.n, self.chi
        cat = (lambda xs, d: self.torch.cat(xs, d)) if istorch else \
              (lambda xs, d: np.concatenate(xs, d))
        M = t.einsum("qasb,qdsc->qadbc", A, A).reshape(n, D * D, D * D)
        P = self._chain_reduce(M, istorch, cat)
        return P[0, 0]

    def logp_from(self, A, bits, norm):
        """log p for a batch of bitstrings, given a built MPS and its norm.

        ``amp(x) = e_0^T (prod_q A[q][:, x_q, :]) e_0``.  The product is a matrix chain, so
        it reduces pairwise in log depth; each level is rescaled by its own max magnitude
        and the logs are accumulated, which keeps the 60-site product away from underflow
        exactly as the shipped per-site normalisation did.
        """
        istorch = not isinstance(A, np.ndarray)
        n, B = self.n, bits.shape[0]
        if istorch:
            t = self.torch
            idx = bits.to(t.int64).T if hasattr(bits, "to") else \
                t.as_tensor(np.asarray(bits, np.int64)).T
            Ms = A.permute(0, 2, 1, 3)[t.arange(n)[:, None], idx].permute(1, 0, 2, 3)
            ls = t.zeros(B, dtype=t.float64)
            cat, absf, logf = t.cat, (lambda x: x.abs()), t.log
            amax = (lambda x: x.abs().amax((2, 3)))
            clamp = (lambda x, m: t.clamp(x, min=m))
        else:
            b = np.asarray(bits, np.int64)
            Ms = np.transpose(np.transpose(A, (0, 2, 1, 3))[np.arange(n)[:, None], b.T],
                              (1, 0, 2, 3))
            ls = np.zeros(B)
            cat, absf, logf = np.concatenate, np.abs, np.log
            amax = (lambda x: np.abs(x).max((2, 3)))
            clamp = (lambda x, m: np.maximum(x, m))
        while Ms.shape[1] > 1:
            m = Ms.shape[1]
            h = m // 2
            P = Ms[:, :2 * h:2] @ Ms[:, 1:2 * h:2]
            if m % 2:
                P = cat([P, Ms[:, -1:]], 1)
            sc = clamp(amax(P), 1e-300)                    # (B, k)
            P = P / sc[:, :, None, None]
            ls = ls + logf(sc).sum(1)
            Ms = P
        amp = absf(Ms[:, 0, 0, :]).sum(-1)
        return 2.0 * (ls + logf(clamp(amp, 1e-300))) - logf(clamp(norm, 1e-300))

    def logp(self, theta: np.ndarray, bits: np.ndarray) -> np.ndarray:
        """Exact log p(bits) in numpy. No autograd, no torch."""
        A = self._build_cached(theta)
        return self.logp_from(A, np.asarray(bits, np.int64), self.norm(A))

    def logp_torch(self, theta, bits):
        """Differentiable log p for a batch of bitstrings (torch tensors in and out)."""
        A = self.build(theta)
        return self.logp_from(A, bits, self.norm(A))

    def probs(self, theta: np.ndarray) -> np.ndarray:
        """Exact distribution over all ``2**n`` bitstrings. Verification only."""
        return np.exp(self.logp(theta, all_bitstrings(self.n)))

    # -- sampling ---------------------------------------------------------
    def _right_envs(self, A):
        """``R[1] .. R[n]`` as one ``(n, chi, chi)`` array, by a suffix scan.

        ``R[q] = M[q] R[q+1]`` is a matrix-chain suffix product, and a Hillis-Steele
        doubling scan computes ALL suffixes in ceil(log2 n) batched matmuls where the
        shipped right-to-left loop needed n sequential ones.  Same products, same
        associativity argument as `_chain_reduce`.
        """
        n, D = self.n, self.chi
        M = np.einsum("qasb,qdsc->qadbc", A, A).reshape(n, D * D, D * D)
        S = M
        step = 1
        while step < n:
            S = np.concatenate([S[:n - step] @ S[step:], S[n - step:]], axis=0)
            step *= 2
        vr = np.zeros(D * D)
        vr[0] = 1.0                                        # R[n] = E_00
        R = np.empty((n + 1, D, D))
        R[:n] = (S @ vr).reshape(n, D, D)
        R[n] = 0.0
        R[n, 0, 0] = 1.0
        return R

    def sample(self, theta: np.ndarray, shots: int, rng) -> np.ndarray:
        """Exact ancestral sampling from the circuit's output distribution.

        Left-to-right conditional decomposition with right environments -- the standard
        exact MPS sampler.  The per-site loop is genuinely sequential (each bit conditions
        the next), so what was hoisted out of it is the environment scan and the
        ``A R A`` contraction, both precomputed for every site in one batched einsum.
        """
        if not self._entangling:
            # Product state: the marginals are the RY angles of all rotation layers
            # composed along the wire, and there is nothing to contract.
            th = np.asarray(theta, float).reshape(self.nblocks, self.n).sum(0)
            return (rng.random((int(shots), self.n))
                    < np.sin(th / 2.0) ** 2).astype(np.uint8)
        A = self._build_cached(theta)
        R = self._right_envs(A)
        n, D = self.n, self.chi
        # Both contractions become matmuls on contiguous blocks. The shipped loop spent
        # two `np.einsum` parses per site -- ~10 us each at n=60, which is most of the
        # sampler -- for arithmetic that is a (B, chi) x (chi, 2 chi) product and a
        # quadratic form. Same operands, same order, bit-identical bits out.
        Af = np.ascontiguousarray(A.reshape(n, D, 2 * D))
        B = int(shots)
        v = np.zeros((B, D))
        v[:, 0] = 1.0
        out = np.empty((B, n), np.uint8)
        u = rng.random((B, n))
        rows = np.arange(B)
        for q in range(n):
            w = (v @ Af[q]).reshape(B, 2, D)
            w2 = w.reshape(B * 2, D)
            p = np.einsum("ij,ij->i", w2 @ R[q + 1], w2).reshape(B, 2)
            np.maximum(p, 0.0, out=p)
            p /= np.maximum(p.sum(1, keepdims=True), 1e-300)
            bit = (u[:, q] > p[:, 0]).astype(np.int64)
            out[:, q] = bit
            v = w[rows, bit]
            v /= np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-300)
        return out

    # -- gradient ---------------------------------------------------------
    def grad_logp(self, theta: np.ndarray, bits: np.ndarray,
                  weights: np.ndarray) -> np.ndarray:
        """``sum_i weights[i] * grad log p(bits[i])`` -- the score-function accumulator."""
        t = self.torch
        th = t.as_tensor(np.asarray(theta, float), dtype=t.float64).requires_grad_(True)
        lp = self.logp_torch(th, t.as_tensor(np.asarray(bits, np.int64)))
        (lp * t.as_tensor(np.asarray(weights, float))).sum().backward()
        return th.grad.detach().numpy()


# ==================================================================== optimiser
class Adam:
    """Adam, verbatim from ``qansatz.Adam``. The optimiser arithmetic is untouched."""

    def __init__(self, n: int, lr: float = 0.12, b1: float = 0.9, b2: float = 0.999):
        self.lr, self.b1, self.b2 = lr, b1, b2
        self.m = np.zeros(n)
        self.v = np.zeros(n)
        self.t = 0

    def step(self, x: np.ndarray, g: np.ndarray) -> np.ndarray:
        self.t += 1
        self.m = self.b1 * self.m + (1 - self.b1) * g
        self.v = self.b2 * self.v + (1 - self.b2) * g * g
        mh = self.m / (1 - self.b1 ** self.t)
        vh = self.v / (1 - self.b2 ** self.t)
        return x - self.lr * mh / (np.sqrt(vh) + 1e-8)


# ============================================================== CVaR gradients
def all_bitstrings(n_qubits: int) -> np.ndarray:
    """``(2**n, n)`` every bitstring, MSB first. VERIFICATION ONLY -- exponential."""
    n = int(n_qubits)
    if n > 20:
        raise ValueError("refusing to enumerate more than 2**20 bitstrings")
    j = np.arange(1 << n, dtype=np.int64)
    return ((j[:, None] >> np.arange(n - 1, -1, -1)) & 1).astype(np.uint8)


def grad_logp_weighted(ansatz, theta: np.ndarray, bits: np.ndarray,
                       weights: np.ndarray) -> np.ndarray:
    """``sum_i w_i grad log p(bits_i)`` for either ansatz family."""
    if isinstance(ansatz, OneLayerAnsatz) or hasattr(ansatz, "_pre_bits"):
        return np.asarray(weights, float) @ ansatz.grad_logp(theta, bits)
    return ansatz.grad_logp(theta, bits, weights)


def cvar_gradient(ansatz, theta: np.ndarray, bits: np.ndarray,
                  energies: np.ndarray, alpha: float,
                  baseline: str = "const") -> Tuple[np.ndarray, float]:
    """Score-function gradient of ``CVaR_alpha`` and the objective value.

    ``baseline``
      ``"const"``  (default, THE FIX) subtract the mean weight from EVERY sample and
                   accumulate over all of them.  ``E_p[grad log p] = 0`` holds over the
                   whole distribution, so this is a genuine control variate: variance
                   down, bias zero.  Measured cos +1.000000 against the exact-expectation
                   reference, and +0.994 as a sampled estimator.
      ``"tail"``   THE RECORDED DEFECT, reproduced verbatim so it stays measurable:
                   centre the weights within the tail and accumulate over the tail only.
                   ``b(x) = m * 1[x in tail]`` is a function of x, the identity does not
                   hold on a data-dependent subset, and the estimator is biased -- cos
                   +0.655634 at 0.758x norm with ZERO sampling noise.  Reproduces
                   ``qansatz.cvar_gradient`` to 5.6e-17.
      ``"none"``   no baseline: unbiased, higher variance.

    The value returned alongside is the sampled CVaR, unchanged in either mode.

    Note the cost asymmetry, and that it is not a reason to prefer the defect: the tail
    form runs the autograd pass on ``alpha * shots`` rows and the constant form on all
    ``shots`` rows.  With the padded MPS the whole batch is one reduction, so the
    corrected gradient is still ~4.4x faster than the defective one was.
    """
    if baseline not in BASELINES:
        raise ValueError(f"baseline must be one of {BASELINES}, got {baseline!r}")
    e = np.asarray(energies, float)
    val, q, tail = cvar(e, alpha)
    N = len(e)
    w = np.zeros(N)
    # (q - E)_+ / alpha, mean-normalised by N: the envelope-theorem weight on the tail.
    w[tail] = -(q - e[tail]) / (alpha * N)
    if baseline == "tail":
        w[tail] -= w[tail].mean()
        # Mirror `qansatz.cvar_gradient` branch for branch, including WHICH rows the
        # accumulation runs over: the closed-form family accumulates over all `shots` rows
        # (the non-tail weights are exactly zero) and only the autograd family is
        # restricted to the tail. Summing k terms and summing N terms of which N-k are zero
        # are the same number in exact arithmetic and one ulp apart in floating point, and
        # this routine exists to reproduce the shipped defect exactly, not approximately.
        if isinstance(ansatz, OneLayerAnsatz) or hasattr(ansatz, "_pre_bits"):
            return w @ ansatz.grad_logp(theta, bits), val
        return ansatz.grad_logp(theta, bits[tail], w[tail]), val
    if baseline == "const":
        w = w - w.mean()
    return grad_logp_weighted(ansatz, theta, bits, w), val


def cvar_gradient_exact(ansatz, theta: np.ndarray, bits_all: np.ndarray,
                        energies: np.ndarray, alpha: float
                        ) -> Tuple[np.ndarray, float]:
    """Exact expectation gradient of ``CVaR_alpha`` -- no sampling noise at all.

    ``grad = -(1/alpha) sum_x (q - E(x))_+ p(x) grad log p(x)``: the envelope-theorem form
    evaluated against the TRUE distribution instead of a sample.  This is the reference the
    finite-difference test checks and the reference the Monte Carlo estimators are measured
    against.  VERIFICATION ONLY -- it enumerates the register.
    """
    p = np.exp(np.asarray(ansatz.logp(theta, bits_all), float))
    p = p / p.sum()
    val, q, _ = cvar_from_probs(energies, p, alpha)
    e = np.asarray(energies, float)
    w = -np.clip(q - e, 0.0, None) * p / alpha
    return grad_logp_weighted(ansatz, theta, bits_all, w), val


# ================================================== exact statevector circuit
class StatevectorCircuit:
    """``layers x (RY on every wire, CNOT chain, ring closure)``, exact real statevector.

    Real-amplitude RY+CNOT is the standard hardware-efficient ansatz and is expressive
    enough for a multimodal distribution over basis states at depth >= 2, which is what a
    distribution over 2-3 populated structural hypotheses actually requires.  The circuit
    is simulated EXACTLY (``2**n`` amplitudes, n <= ~20), so ``p_theta`` and every gradient
    below are analytic rather than sampled -- which is what makes the gradient audit a real
    verification instead of a noise comparison.

    TWO INVARIANTS HOISTED, ONE LOOP BATCHED
    ----------------------------------------
    1. The entangling layer is a fixed PERMUTATION of basis indices and does not depend on
       the parameters at all.  The shipped version rebuilt ``np.arange(dim)``, a shift, a
       comparison and a ``np.where`` for every CNOT of every layer of every call; here the
       whole chain (ring included) is composed into ONE index array in ``__init__`` and
       applied as a single gather per layer.  ``n`` gathers become 1.
    2. ``probs_batch`` simulates B parameter vectors at once, so parameter shift --
       inherently ``2P`` evaluations of the same circuit -- costs ONE pass over the gate
       list instead of ``2P`` passes.  At n=7, layers=3 that is 24 numpy calls in place of
       1,806.

    Neither changes the arithmetic: the RY update is the same two elementwise expressions
    on the same slices, and a composed permutation is exact.  ``probs`` is BIT-identical to
    the shipped loop (measured max |diff| 0.0 at n=7 and n=10); the parameter-shift gradient
    agrees to 2e-16, the only difference being that the ``dCVaR/dp`` contraction is now one
    matmul instead of ``P`` separate dot products.
    ``t_statevector_probs_bit_identical_to_legacy`` pins both.
    """

    def __init__(self, n: int, layers: int = 3, ring: bool = True):
        self.n, self.layers, self.ring = int(n), int(layers), bool(ring)
        self.dim = 1 << self.n
        self._perm = self._entangler_permutation()

    def _entangler_permutation(self) -> np.ndarray:
        """Compose the whole CNOT chain (+ ring) into one basis-index gather."""
        n, dim = self.n, self.dim
        cur = np.arange(dim, dtype=np.int64)
        pairs = [(q, q + 1) for q in range(n - 1)]
        if self.ring and n > 2:
            pairs.append((n - 1, 0))
        idx = np.arange(dim, dtype=np.int64)
        for ctrl, tgt in pairs:
            cb = (idx >> (n - 1 - ctrl)) & 1
            j = np.where(cb == 1, idx ^ (1 << (n - 1 - tgt)), idx)
            cur = cur[j]
        return cur

    def n_params(self) -> int:
        return self.n * self.layers

    # -- batched simulation ------------------------------------------------
    def states_batch(self, thetas: np.ndarray) -> np.ndarray:
        """``(B, n_params)`` angles -> ``(B, dim)`` real amplitudes."""
        TH = np.atleast_2d(np.asarray(thetas, float))
        B = TH.shape[0]
        n, dim = self.n, self.dim
        TH = TH.reshape(B, self.layers, n)
        # Every rotation angle's cosine and sine in two calls rather than 2 * layers * n
        # calls on scalars -- the values are identical, only the dispatch count moves.
        CO = np.cos(TH / 2.0)
        SI = np.sin(TH / 2.0)
        psi = np.zeros((B, dim))
        psi[:, 0] = 1.0
        for L in range(self.layers):
            for q in range(n):
                c = CO[:, L, q][:, None, None]
                s = SI[:, L, q][:, None, None]
                v = psi.reshape(B, 1 << q, 2, 1 << (n - q - 1))
                a = v[:, :, 0, :].copy()
                b = v[:, :, 1, :]                # a view; index 0 is written first
                v[:, :, 0, :] = c * a - s * b
                v[:, :, 1, :] = s * a + c * b
                psi = v.reshape(B, dim)
            psi = psi[:, self._perm]
        return psi

    def probs_batch(self, thetas: np.ndarray) -> np.ndarray:
        p = self.states_batch(thetas) ** 2
        return p / p.sum(1, keepdims=True)

    def state(self, theta: np.ndarray) -> np.ndarray:
        return self.states_batch(np.asarray(theta, float)[None, :])[0]

    def probs(self, theta: np.ndarray) -> np.ndarray:
        return self.probs_batch(np.asarray(theta, float)[None, :])[0]

    # -- the shift grid, built once per gradient ---------------------------
    def _shift_grid(self, theta: np.ndarray, h: float) -> np.ndarray:
        """``(2P, P)`` stack of ``theta`` with +h / -h on each coordinate in turn."""
        th = np.asarray(theta, float)
        P = th.size
        TH = np.repeat(th[None, :], 2 * P, axis=0)
        r = np.arange(P)
        TH[2 * r, r] += h
        TH[2 * r + 1, r] -= h
        return TH


def cvar_exact(E: np.ndarray, p: np.ndarray, alpha: float
               ) -> Tuple[float, float, np.ndarray]:
    """CVaR of the LOWER `alpha` tail of the energy distribution, from exact `p`.

    ``CVaR_alpha = (1/alpha)[ sum_{E(x)<q} p(x) E(x) + (alpha - P(E<q)) q ]`` with ``q`` the
    alpha-quantile.  Returns ``(value, q, dCVaR/dp)``.  The derivative uses the envelope
    theorem -- q's own dependence on p cancels -- so ``dp`` is exactly ``(E - q)/alpha`` on
    the strict tail and 0 elsewhere, with the boundary state carrying its partial mass.
    """
    E, p = np.asarray(E, float), np.asarray(p, float)
    o = np.argsort(E, kind="stable")
    c = np.cumsum(p[o])
    k = int(np.searchsorted(c, alpha, side="left"))
    k = min(k, len(E) - 1)
    q = E[o[k]]
    below = o[:k]
    mass = c[k - 1] if k > 0 else 0.0
    val = (float((p[below] * E[below]).sum()) + (alpha - mass) * q) / alpha
    dp = np.zeros_like(p)
    dp[below] = (E[below] - q) / alpha
    return float(val), float(q), dp


def grad_cvar_paramshift(circ: StatevectorCircuit, theta: np.ndarray,
                         E: np.ndarray, alpha: float) -> np.ndarray:
    """EXACT gradient by the parameter-shift rule applied to each basis probability.

    ``p(x) = <psi|Pi_x|psi>`` is the expectation of a projector, so every ``p(x)`` obeys the
    exact two-term shift rule for an RY generator.  Chaining it with ``dCVaR/dp`` gives the
    exact analytic gradient with no sampling and no finite-difference step.  This is the
    reference the deployable estimator is audited against.

    All ``2P`` shifted circuits are simulated in one batched pass; the shift rule itself is
    untouched, and the result is bit-identical to evaluating them one at a time.
    """
    theta = np.asarray(theta, float)
    _, _, dp = cvar_exact(E, circ.probs(theta), alpha)
    PR = circ.probs_batch(circ._shift_grid(theta, np.pi / 2))
    return (PR[0::2] - PR[1::2]) @ dp / 2.0


def grad_cvar_fd(circ: StatevectorCircuit, theta: np.ndarray, E: np.ndarray,
                 alpha: float, h: float = 1e-5) -> np.ndarray:
    """Central finite differences on the exact CVaR value.

    Independent machinery from the shift rule -- no projector algebra, no envelope theorem
    -- which is exactly what makes it a real check on `grad_cvar_paramshift` rather than a
    restatement of it.
    """
    theta = np.asarray(theta, float)
    PR = circ.probs_batch(circ._shift_grid(theta, h))
    vals = np.array([cvar_exact(E, PR[i], alpha)[0] for i in range(PR.shape[0])])
    return (vals[0::2] - vals[1::2]) / (2 * h)


def grad_cvar_score(circ: StatevectorCircuit, theta: np.ndarray, E: np.ndarray,
                    alpha: float, shots: int = 8192, rng=None,
                    baseline: str = "const") -> Tuple[np.ndarray, float]:
    """The SAMPLED score-function estimator -- what a real device would compute.

    ``baseline`` is the whole point of this function; see `cvar_gradient`.  ``"const"``
    subtracts a constant from every weight and is unbiased; ``"tail"`` reproduces the
    recorded defect verbatim (subtract the tail mean from the tail entries only), so the
    audit MEASURES the bias rather than asserting it.
    """
    if baseline not in BASELINES:
        raise ValueError(f"baseline must be one of {BASELINES}, got {baseline!r}")
    rng = np.random.default_rng(0) if rng is None else rng
    theta = np.asarray(theta, float)
    p = circ.probs(theta)
    x = rng.choice(len(p), size=int(shots), p=p)
    en = np.asarray(E, float)[x]
    val, q, _ = cvar_exact(E, p, alpha)
    w = np.where(en < q, (en - q) / alpha, 0.0)
    if baseline == "const":
        w = w - w.mean()                    # constant in x -> unbiased
    elif baseline == "tail":
        m = en < q                          # the recorded DEFECT, reproduced verbatim
        if m.any():
            w[m] -= w[m].mean()
    # d log p(x) / d theta_k by the same shift rule, evaluated once per basis state and
    # for all parameters in ONE batched simulation.
    PR = circ.probs_batch(circ._shift_grid(theta, np.pi / 2))
    G = ((PR[0::2] - PR[1::2]) / 2.0 / np.maximum(p, 1e-15)[None, :]).T
    return (w[:, None] * G[x]).mean(0), val


def free_energy(circ: StatevectorCircuit, theta: np.ndarray, E: np.ndarray,
                alpha: float, T: float):
    """``F = CVaR_alpha(E; p_theta) - T H(p_theta)``, and its exact gradient.

    WHY AN ENTROPY TERM, AND WHY IT IS NOT A FUDGE
    ----------------------------------------------
    Minimising CVaR alone is degenerate for this task: for ANY alpha the minimiser
    concentrates p on the lowest-energy basis states, so the readout collapses back to the
    argmin -- the shipped selector, 3.454 A.  That was measured before this term existed
    (every ``vqe_a*`` arm returned the argmin's structure, state entropy 0.01 bits at
    alpha=1) and it is a property of CVaR, not of the optimiser.

    The quantity the consensus readout needs is an ENSEMBLE: concentrated on good
    hypotheses but still broad enough to have a centre.  That is a free energy, and it is
    the same object the project's ensemble-selection result rests on.  ``alpha`` says which
    part of the energy distribution is scored and ``T`` says how broad the ensemble is;
    both are measured, neither is tuned toward a native.

    Returns ``(F, grad F, p, CVaR, H)``.
    """
    p = circ.probs(theta)
    v, q, dp = cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-15))
    H = float(-(p * lp).sum())
    dH = -(lp + 1.0)
    d = dp - T * dH
    PR = circ.probs_batch(circ._shift_grid(np.asarray(theta, float), np.pi / 2))
    g = (PR[0::2] - PR[1::2]) @ d / 2.0
    return float(v - T * H), g, p, v, H


def run_cvar_vqe(E: np.ndarray, alpha: float, T: float = 0.0, n: int = 7,
                 layers: int = 3, iters: int = 80, restarts: int = 1,
                 seed: int = 0, lr: float = 0.15):
    """Optimise the CVaR free energy of the hypothesis distribution; return ``p_theta``.

    Adam on the exact parameter-shift gradient.  The gradient is the audited one; the
    sampled estimator agrees with it at cosine 0.994, so using the exact one here removes
    shot noise from a selection experiment without changing what is being optimised.

    Iterations, restarts, learning rate, seed handling and the ansatz are unchanged from
    the shipped driver; only the per-gradient simulation got batched.
    """
    circ = StatevectorCircuit(n, layers)
    rng = np.random.default_rng(seed)
    best, best_f = None, np.inf
    for r in range(restarts):
        th = rng.normal(0.0, 0.6, circ.n_params())
        m = np.zeros_like(th)
        v = np.zeros_like(th)
        for t in range(1, iters + 1):
            f, g, _, _, _ = free_energy(circ, th, E, alpha, T)
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g * g
            th = th - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
        f, _, p, cv, H = free_energy(circ, th, E, alpha, T)
        if f < best_f:
            best_f, best = f, (p, cv, H)
    return best[0], float(best[1]), float(best[2]), circ


# ============================================== PennyLane global CVaR-VQE
def n_parameters(n_qubits: int, layers: int) -> int:
    return layers * n_qubits


def build_global_circuit(n_qubits: int, layers: int, ring: bool = True,
                         device: str = "lightning.qubit") -> Callable:
    """One circuit over ALL n_qubits. Returns probs over the full register.

    A REAL DEVICE, EXECUTED ONCE PER TAPE INSTEAD OF REBUILT PER CALL
    ----------------------------------------------------------------
    The gate list, the wire map and the measurement are identical on every objective call
    -- only the ``layers * n_qubits`` rotation angles change.  The shipped version wrapped
    a QNode, which reconstructs the tape from the Python function and re-runs the transform
    program on every single evaluation; at n=12, layers=4 that was 3.5 ms of which 2.5 ms
    was construction.  Here the tape is built once and ``bind_new_parameters`` swaps the
    angles in, then the device executes it directly.

    This is still PennyLane and still ``lightning.qubit``: same simulator, same gates, same
    measurement.  The outputs are BIT-identical (measured max |diff| 0.0 at n=12, 16 and
    20), and at n=20 -- where actual simulation, not construction, dominates -- the two
    paths cost the same, which is exactly the signature of an overhead-only fix.
    """
    import pennylane as qml
    try:
        dev = qml.device(device, wires=n_qubits)
    except Exception:
        dev = qml.device("default.qubit", wires=n_qubits)

    n_par = n_parameters(n_qubits, layers)
    try:
        ops = []
        for _l in range(layers):
            for q in range(n_qubits):
                ops.append(qml.RY(0.0, wires=q))
            for q in range(n_qubits - 1):
                ops.append(qml.CNOT(wires=[q, q + 1]))
            if ring and n_qubits > 2:
                ops.append(qml.CNOT(wires=[n_qubits - 1, 0]))
        tape = qml.tape.QuantumScript(ops, [qml.probs(wires=range(n_qubits))])
        slots = list(range(n_par))
        # one probe execution: if anything about this device/tape combination is unusual,
        # fall back rather than fail deep inside an optimiser
        _probe = np.linspace(0.1, 0.2, n_par)
        dev.execute(tape.bind_new_parameters(list(_probe), slots))

        def circuit(params):
            p = np.asarray(params, dtype=float).ravel()
            return dev.execute(tape.bind_new_parameters(list(p), slots))

        return circuit
    except Exception:
        @qml.qnode(dev)
        def circuit_qnode(params):
            p = np.reshape(np.asarray(params, dtype=float), (layers, n_qubits))
            for l in range(layers):
                for q in range(n_qubits):
                    qml.RY(float(p[l][q]), wires=q)
                for q in range(n_qubits - 1):
                    qml.CNOT(wires=[q, q + 1])
                if ring and n_qubits > 2:
                    qml.CNOT(wires=[n_qubits - 1, 0])
            return qml.probs(wires=range(n_qubits))

        return circuit_qnode


class BestSeenTracker:
    """Records the lowest-energy bitstring seen anywhere during a run."""

    def __init__(self):
        self.best_energy = float("inf")
        self.best_bitstring: Optional[str] = None
        self.n_lookups = 0

    def offer(self, bitstring: str, energy: float) -> None:
        self.n_lookups += 1
        if energy < self.best_energy:
            self.best_energy = energy
            self.best_bitstring = bitstring


class _SPSAResult:
    def __init__(self, x, fun):
        self.x = x
        self.fun = fun


def _spsa(objective, x0, n_iter, rng, a=0.25, c=0.15,
          set_sample_tag=None, progress=None):
    """Simultaneous Perturbation Stochastic Approximation (Spall). Arithmetic untouched.

    This is the right family of optimizer for this objective, and the reason is structural
    rather than empirical.  The CVaR objective is a *stochastic* estimate: it draws
    ``shots`` bitstrings from the circuit's distribution, so evaluating the same parameters
    twice returns different values.  COBYLA is a derivative-free **trust-region** method
    that builds a linear model from n+1 points and shrinks its region when predicted and
    actual improvement disagree -- which sampling noise guarantees at every step.  SPSA's
    convergence theory is *for* noisy evaluations, and its cost per iteration is 2
    evaluations regardless of dimension, against COBYLA's n+1 just to build an initial
    simplex (89 evaluations at 88 parameters).

    Three fixes over the pre-2024 implementation, all of which mattered:

    1. **Progress-driven schedule.**  ``A`` and the decay exponents were calibrated against
       ``n_iter``, but the runner deliberately sets ``maxiter`` far above what the shared
       evaluation budget affords, so the budget terminates the run long before ``k``
       approaches ``n_iter``.  With ``n_iter = 50 x n_params = 4000`` that made ``A = 400``
       and a first step of ``a / 401**0.602 = 0.007`` -- five times too small -- and the
       decay never completed.  ``progress()`` returns the fraction of the budget slice
       consumed and the schedule is traversed against that.
    2. **Common random numbers.**  Both arms of the finite difference are evaluated against
       the *same* bitstring sample stream, so the sampling noise largely cancels in
       ``fp - fm``.
    3. **Two evaluations per iteration, not three.**  The midpoint ``(fp + fm) / 2`` is an
       unbiased estimate of ``f(x)`` with *lower* variance than a fresh evaluation, so the
       old third call was pure waste: 1 + 3n evaluations became 2n.
    """
    x = np.array(x0, dtype=float)
    A = max(1, n_iter // 10)
    best_x, best_f = x.copy(), float("inf")
    for k in range(n_iter):
        # Effective index along the intended schedule. Budget-driven when a budget exists,
        # so the decay completes exactly as the slice runs out; otherwise the raw index.
        p = None if progress is None else progress()
        kk = k if p is None else p * n_iter
        ak = a / ((kk + 1.0 + A) ** 0.602)
        ck = c / ((kk + 1.0) ** 0.101)

        d = rng.choice([-1.0, 1.0], size=x.size)
        if set_sample_tag is not None:
            set_sample_tag(k)               # common random numbers across the +/- pair
        fp = objective(x + ck * d)
        fm = objective(x - ck * d)

        f_mid = 0.5 * (fp + fm)
        if f_mid < best_f:
            best_f, best_x = f_mid, x.copy()

        x = x - ak * (fp - fm) / (2.0 * ck) * d
    if set_sample_tag is not None:
        set_sample_tag(None)                # restore per-call streams
    return _SPSAResult(best_x, best_f)


#: Fewest SPSA iterations at which the optimiser has been observed to do useful work on
#: this objective. Below this the run returns something close to its initialisation.
#: Measured: the shipped configuration (shots 2048, budget 20000, restarts 4) affords about
#: 33, and returned a bit-identical answer across 3 seeds and 5 hyperparameter settings --
#: the signature of an optimiser that never moved.
MIN_USEFUL_SPSA_ITERS = 100


def _warn_if_underoptimised(hamiltonian, shots: int, restarts: int,
                            optimizer: str, verbose: bool) -> Optional[str]:
    """Warn when the shared budget cannot afford enough optimiser steps.

    ``hamiltonian.energy`` charges per UNIQUE bitstring, so while the circuit distribution
    is still broad, one objective call costs ~``shots`` of the budget.  The affordable
    number of SPSA iterations is therefore about ``budget / (shots * restarts * 2)`` early
    on -- it improves as the distribution concentrates and samples start hitting the cache,
    but the early phase is exactly when the optimiser needs to move.

    A warning rather than an error because a small run is legitimately allowed to be
    under-optimised; what is not acceptable is for a headline comparison to be
    under-optimised without saying so.
    """
    budget = getattr(hamiltonian, "eval_budget", None)
    if not budget or optimizer.upper() != "SPSA":
        return None
    affordable = budget / max(1, shots * restarts * 2)
    if affordable >= MIN_USEFUL_SPSA_ITERS:
        return None
    msg = (f"VQE is budget-starved: shots={shots} x restarts={restarts} against a shared "
           f"budget of {budget} affords about {affordable:.0f} SPSA iterations early on, "
           f"below the {MIN_USEFUL_SPSA_ITERS} at which this objective has been seen to "
           f"make progress. The result will be close to the initialisation and will look "
           f"deceptively stable across seeds. Reduce shots (try "
           f"{max(32, int(budget / (MIN_USEFUL_SPSA_ITERS * restarts * 2)))}), reduce "
           f"restarts, or raise the budget.")
    warnings.warn(msg, RuntimeWarning, stacklevel=3)
    if verbose:
        print(f"    !! {msg}")
    return msg


def _run_single(hamiltonian, circuit, n_qubits: int, layers: int,
                alpha: float, shots: int, maxiter: int, seed: int,
                optimizer: str, tracker: BestSeenTracker,
                init_scale: float, verbose: bool, sampler=None,
                energy_batch: Optional[Callable] = None,
                spsa_a: float = 0.25, spsa_c: float = 0.15,
                trace_interval: int = 250) -> Dict:
    from scipy.optimize import minimize

    fmt = f"0{n_qubits}b"
    n_par = n_parameters(n_qubits, layers)

    init_rng = np.random.default_rng(np.random.SeedSequence([seed, 0xC0FFEE]))
    params0 = (math.pi / 2.0) + init_rng.normal(0.0, init_scale, size=n_par)

    history: List[float] = []
    trace: List[Dict] = []
    interval = max(1, int(trace_interval))
    spent_at_start = int(hamiltonian.n_energy_evaluations)
    next_trace = {"n": ((spent_at_start // interval) + 1) * interval}
    eval_counter = {"n": 0}
    # Which bitstring-sample stream the next objective call draws from. `None` means "a
    # fresh stream per call", which is what a deterministic optimizer like COBYLA gets.
    # SPSA pins it so both arms of its finite difference share one stream.
    sample_tag = {"v": None}

    def objective(params: np.ndarray) -> float:
        k = eval_counter["n"]
        eval_counter["n"] = k + 1
        tag = k if sample_tag["v"] is None else sample_tag["v"]
        rng = np.random.default_rng(np.random.SeedSequence([seed, tag]))

        if sampler is None:
            probs = np.asarray(circuit(params), dtype=float)
            probs = np.clip(probs, 0.0, None)
            s = probs.sum()
            if s <= 0:
                probs = np.full_like(probs, 1.0 / probs.size)
            else:
                probs = probs / s
            idx = rng.choice(probs.size, size=shots, p=probs)
        else:
            idx = sampler.draw_indices(params, shots, rng)

        uniq, inverse = np.unique(idx, return_inverse=True)
        bslist = [format(int(u), fmt) for u in uniq]
        emap = energy_batch(bslist) if energy_batch is not None else None
        uniq_energies = np.empty(uniq.size, dtype=float)
        for m, bs in enumerate(bslist):
            # A budget-aware batch may return only the prefix that fit. Falling back to
            # .energy() for a missing entry preserves normal BudgetExhausted termination.
            e = (emap[bs] if emap is not None and bs in emap
                 else hamiltonian.energy(bs))
            uniq_energies[m] = e
            tracker.offer(bs, e)
        sample_energies = uniq_energies[inverse]   # length == shots

        val = cvar_from_samples(sample_energies, alpha)
        history.append(val)
        if hamiltonian.n_energy_evaluations >= next_trace["n"]:
            trace.append({
                "n_energy_evaluations": int(hamiltonian.n_energy_evaluations),
                "objective_eval": int(k + 1),
                "cvar": float(val),
                "best_energy": float(tracker.best_energy),
                "best_bitstring": tracker.best_bitstring,
                "elapsed_s": float(time.time() - t0),
            })
            next_trace["n"] += interval
        if verbose and (k % 25 == 0):
            print(f"      eval {k:4d} | CVaR {val:9.3f} | "
                  f"best seen {tracker.best_energy:9.3f}")
        return val

    # Best parameters seen by the objective itself, so an interrupted run still has an
    # answer. Without this, exhausting the budget mid-COBYLA would leave no `res.x`.
    best = {"f": float("inf"), "x": np.array(params0, dtype=float)}
    _raw_objective = objective

    def objective(params):                          # noqa: F811 -- wraps the above
        val = _raw_objective(params)
        if val < best["f"]:
            best["f"], best["x"] = val, np.array(params, dtype=float)
        return val

    # Fraction of THIS restart's budget slice consumed, for SPSA's step schedule. Captured
    # at entry because `n_energy_evaluations` is cumulative across restarts.
    _span = max(1.0, float(hamiltonian.budget_remaining)) \
        if getattr(hamiltonian, "eval_budget", None) else None

    def progress() -> Optional[float]:
        if _span is None:
            return None
        return min(1.0, max(0.0, 1.0 - float(hamiltonian.budget_remaining) / _span))

    t0 = time.time()
    terminated_by = "optimizer"
    try:
        if optimizer.upper() == "COBYLA":
            res = minimize(objective, params0, method="COBYLA",
                           options={"maxiter": maxiter, "rhobeg": 0.4})
        elif optimizer.upper() == "SPSA":
            res = _spsa(objective, params0, maxiter,
                        np.random.default_rng(np.random.SeedSequence([seed, 7])),
                        set_sample_tag=lambda t: sample_tag.__setitem__("v", t),
                        progress=progress, a=spsa_a, c=spsa_c)
        else:
            raise ValueError(f"unknown optimizer {optimizer!r}; use COBYLA or SPSA")
        final_params = np.asarray(res.x if hasattr(res, "x") else res, dtype=float)
        final_objective = (float(res.fun) if hasattr(res, "fun")
                           else float(objective(final_params)))
    except BudgetExhausted:
        # A normal termination condition, not a failure: the arm spent its share of the
        # shared evaluation budget. Return the best point the objective actually saw.
        terminated_by = "budget"
        final_params, final_objective = best["x"], best["f"]
    runtime = time.time() - t0

    return {
        "final_params": final_params,
        "final_objective": final_objective,
        "history": history,
        "trace": trace,
        "n_objective_evals": eval_counter["n"],
        "terminated_by": terminated_by,
        "runtime": runtime,
    }


def run_global_cvar_vqe(hamiltonian, layers: int = 4, alpha: float = 0.15,
                        shots: int = 2048, maxiter: Optional[int] = None,
                        restarts: int = 4, seed: int = 0,
                        optimizer: str = "SPSA", ring: bool = True,
                        final_shots: int = 8192, init_scale: float = 1.0,
                        device: str = "lightning.qubit",
                        readout_reserve_frac: float = 0.01,
                        sampler=None, energy_batch: Optional[Callable] = None,
                        spsa_a: float = 0.25, spsa_c: float = 0.15,
                        trace_interval: int = 250,
                        verbose: bool = False) -> Dict:
    """Global CVaR-VQE over the entire protein configuration register.

    Cost is governed by ``hamiltonian.eval_budget`` -- the number of *unique* structures
    whose energy had to be computed -- not by ``maxiter``.  Every search method is handed
    the same budget through the Hamiltonian it shares, so whichever finds the lowest energy
    inside that budget wins and no method can buy a better answer with a larger optimizer
    allowance.  ``maxiter=None`` resolves deliberately high so the shared budget binds.

    Restarts are NOT parallelised here, and that is deliberate: the budget is sliced
    cumulatively across them (restart r may spend its equal share plus whatever earlier
    restarts left), so they are sequentially dependent by construction.  ``fold`` below,
    whose restarts share nothing but a reservoir, is where the free parallelism lives.

    ``init_scale`` and the shots/budget ratio are LOAD-BEARING TOGETHER, and this is the
    single most consequential thing to know about this function.  Measured on 1UAO at 400
    SPSA iterations, changing nothing but the initial spread:

        init_scale 0.25   ->  E -11.018   CA-RMSD 5.27 A   (the alpha-helix)
        init_scale 1.00   ->  E -11.350   CA-RMSD 1.96 A   (the global minimum's backbone)

    Neither variable does anything on its own, which is why this went undiagnosed.  Wide
    init with too few iterations never concentrates -- the distribution is still
    near-uniform (17-20 bits of a possible 22) when the run ends, i.e. random sampling with
    extra steps.  Enough iterations from a narrow init converges reliably into the SAME
    basin every time: at ``init_scale=0.25`` every parameter starts within 0.25 rad of
    pi/2, so all restarts and all seeds begin at effectively one point and reach one answer.
    That is why 5 hyperparameter settings and 3 seeds all returned a bit-identical -11.018.

    ``sampler`` selects an alternative exact sampler without forking the optimizer logic.
    ``energy_batch`` may evaluate the independent unique structures of one objective call
    concurrently; the Hamiltonian remains the sole owner of the shared hard budget.
    """
    _warn_if_underoptimised(hamiltonian, shots, restarts, optimizer, verbose)

    n_qubits = hamiltonian.n_qubits
    if sampler is None and n_qubits > 30:
        raise MemoryError(
            f"n_qubits={n_qubits} requires ~{2**n_qubits * 16 / 1e9:.0f} GB "
            "for a statevector. A genuine full-system VQE is not simulable "
            "at this size. Reduce protein length or state count.")

    n_par = n_parameters(n_qubits, layers)
    maxiter = resolve_maxiter(maxiter, n_par)
    check_optimizer_budget(maxiter, n_par, optimizer,
                           getattr(hamiltonian, "eval_budget", None))

    circuit = (None if sampler is not None else
               build_global_circuit(n_qubits, layers, ring=ring, device=device))
    tracker = BestSeenTracker()

    hamiltonian.reset_counters()

    # Withhold a slice for the final read-out so a run can always report its answer after
    # the optimizer has spent everything. The read-out scans most-probable-first and those
    # states are nearly all cache hits after a converged search, so this needs to be small
    # -- reserving 10% measurably penalised the VQE's search relative to the classical
    # arms, which is precisely the unfairness the shared budget exists to remove.
    budget = getattr(hamiltonian, "eval_budget", None)
    readout_reserve = 0 if budget is None else max(1, int(budget * readout_reserve_frac))

    restart_seeds = [int(s.generate_state(1)[0])
                     for s in np.random.SeedSequence(seed).spawn(restarts)]

    t0 = time.time()
    runs = []
    best_run = None
    restarts_completed = 0
    for r, rseed in enumerate(restart_seeds):
        if verbose:
            print(f"    --- restart {r + 1}/{restarts} (seed {rseed}) ---")
        if budget is not None:
            # Slice the pool cumulatively: restart r may spend up to its equal share plus
            # anything earlier restarts left unspent. Without this, restart 1 consumes the
            # whole pool and a 4-restart algorithm silently becomes a 1-restart one.
            searchable = budget - readout_reserve
            cum_limit = int(round(searchable * (r + 1) / restarts))
            hamiltonian.reserve(budget - cum_limit)
            if hamiltonian.budget_remaining <= 0:
                if verbose:
                    print(f"    (budget exhausted; {restarts - r} restarts not run)")
                break
        run = _run_single(hamiltonian, circuit, n_qubits, layers, alpha,
                          shots, maxiter, rseed, optimizer, tracker,
                          init_scale, verbose, sampler=sampler,
                          energy_batch=energy_batch, spsa_a=spsa_a,
                          spsa_c=spsa_c, trace_interval=trace_interval)
        if (not run["trace"] or
                run["trace"][-1]["n_energy_evaluations"] !=
                hamiltonian.n_energy_evaluations):
            run["trace"].append({
                "n_energy_evaluations": int(hamiltonian.n_energy_evaluations),
                "objective_eval": int(run["n_objective_evals"]),
                "cvar": (float(run["history"][-1]) if run["history"] else float("nan")),
                "best_energy": float(tracker.best_energy),
                "best_bitstring": tracker.best_bitstring,
                "elapsed_s": float(run["runtime"]),
            })
        for point in run["trace"]:
            point["restart"] = r
        runs.append(run)
        restarts_completed += 1
        if best_run is None or run["final_objective"] < best_run["final_objective"]:
            best_run = run
    total_runtime = time.time() - t0

    if best_run is None:
        raise BudgetExhausted(budget or 0, hamiltonian.n_energy_evaluations)

    hamiltonian.release()               # hand the reserve back for the read-out

    fmt = f"0{n_qubits}b"
    final_rng = np.random.default_rng(np.random.SeedSequence([seed, 0xF1A1]))
    if sampler is None:
        final_probs = np.asarray(circuit(best_run["final_params"]), dtype=float)
        final_probs = np.clip(final_probs, 0.0, None)
        final_probs = final_probs / final_probs.sum()
        final_idx = final_rng.choice(final_probs.size, size=final_shots, p=final_probs)
        modal_bits = format(int(np.argmax(final_probs)), fmt)
        p_sorted = np.sort(final_probs)[::-1]
        top1 = float(p_sorted[0])
        top16 = float(p_sorted[:16].sum())
        nz = final_probs[final_probs > 1e-15]
        entropy = float(-np.sum(nz * np.log2(nz)))
    else:
        final_idx = sampler.draw_indices(best_run["final_params"], final_shots, final_rng)
        vals, counts = np.unique(final_idx, return_counts=True)
        freq = counts.astype(float) / float(final_shots)
        modal_bits = format(int(vals[np.argmax(counts)]), fmt)
        top1 = float(freq.max())
        top16 = float(np.sort(freq)[::-1][:16].sum())
        entropy = float(-np.sum(freq * np.log2(freq)))
    # Scan most-probable-first. np.unique returns basis indices in *sorted* order, so if the
    # budget truncates the read-out a plain unique() scan keeps an arbitrary low-index
    # subset rather than the states the circuit actually favours.
    uniq, counts = np.unique(final_idx, return_counts=True)
    uniq = uniq[np.argsort(-counts, kind="stable")]
    readout_bits = [format(int(u), fmt) for u in uniq]
    readout_energies = (energy_batch(readout_bits)
                        if energy_batch is not None else None)
    vqe_bits, vqe_energy = None, float("inf")
    for bs in readout_bits:
        try:
            e = (readout_energies[bs]
                 if readout_energies is not None and bs in readout_energies
                 else hamiltonian.energy(bs))
        except BudgetExhausted:
            break
        if e < vqe_energy:
            vqe_energy, vqe_bits = e, bs

    try:
        modal_energy = hamiltonian.energy(modal_bits)
    except BudgetExhausted:
        modal_energy = float("nan")

    objective_trace = []
    for restart, run in enumerate(runs):
        for point in run["trace"]:
            q = dict(point)
            q["restart"] = restart
            objective_trace.append(q)
    return {
        "vqe_bitstring": vqe_bits,
        "vqe_energy": float(vqe_energy),
        "vqe_modal_bitstring": modal_bits,
        "vqe_modal_energy": float(modal_energy),
        "best_seen_bitstring": tracker.best_bitstring,
        "best_seen_energy": float(tracker.best_energy),
        "final_objective": best_run["final_objective"],
        "history": best_run["history"],
        "objective_trace": objective_trace,
        "distribution_top1_prob": top1,
        "distribution_top16_mass": top16,
        "distribution_entropy_bits": entropy,
        "max_entropy_bits": float(n_qubits),
        "n_qubits": n_qubits,
        "n_parameters": n_parameters(n_qubits, layers),
        "layers": layers,
        "alpha": alpha,
        "shots_per_eval": shots,
        "final_shots": final_shots,
        "restarts": restarts,
        "optimizer": optimizer,
        "spsa_a": float(spsa_a),
        "spsa_c": float(spsa_c),
        "backend": (getattr(sampler, "backend_name", "sampler")
                    if sampler is not None else "statevector"),
        "n_objective_evals_total": sum(r["n_objective_evals"] for r in runs),
        "n_objective_evals_best_run": best_run["n_objective_evals"],
        "n_spsa_iterations_total": (sum(r["n_objective_evals"] for r in runs) // 2
                                    if optimizer.upper() == "SPSA" else 0),
        "n_energy_evaluations": hamiltonian.n_energy_evaluations,
        "n_unique_structures_cached": hamiltonian.cache_size(),
        "maxiter_resolved": maxiter,
        "maxiter_over_n_params": maxiter / max(1, n_par),
        "eval_budget": budget,
        "budget_exhausted": bool(getattr(hamiltonian, "budget_exhausted", False)),
        "terminated_by": ("budget" if any(r["terminated_by"] == "budget" for r in runs)
                          else "optimizer"),
        "restarts_completed": restarts_completed,
        "runtime": total_runtime,
        "seed": seed,
    }


# ================================================================== objective
# Geometry and the classical energy terms come through the package's backend switch, not
# by importing the root modules, so that a consolidated `core.geometry` / `core.energy`
# going live is the act of committing the file -- and so `core.bench` records which
# implementation this module actually ran against in its cache key. `backend()` falls back
# to the root module whenever the consolidated one is absent or incomplete, and it is
# called here at import time rather than per call because these are module handles, not
# work. There is no import cycle: `core/__init__` imports nothing at module level, and
# neither `core.geometry` nor `core.energy` imports `core.quantum`.
#
# `sidechains` has no backend entry (its consolidation lives inside `core.amber`, which
# this module does not use), so it is imported directly.
import core as _core                                               # noqa: E402
import sidechains as sc                                            # noqa: E402

geo = _core.backend("geometry")
et = _core.backend("energy")


class FoldObjective:
    """Batched composite folding score over discrete torsion-state assignments.

    Four terms, each earning its place from a measurement rather than from physical
    plausibility:

    ``distogram``  weighted L1 Bayes risk against the predicted CA-CA distance
                   distribution. This is the only term with real near-native
                   discrimination.
    ``mrf``        the local torsion prior: per-residue state marginals plus
                   nearest-neighbour couplings, estimated from the held-out database. It is
                   a *local* term, so on its own it has no fold in it -- its job is to keep
                   the chain in populated Ramachandran basins and to supply the VQE's warm
                   start, not to choose a tertiary arrangement.
    ``clash``      excluded volume on CA and CB. A feasibility filter, not a
                   discriminator: flat across everything the search should be choosing
                   between, and steeply positive on the structures it must not return.
    ``legacy``     the coarse knowledge-based potential, standardised against a random
                   sample from this target's own state library and blended in units of the
                   distance term's own spread on that same sample, so one weight transfers
                   across targets and chain lengths. Off by default.

    Everything is a function of the ``(B, n)`` state matrix, and every step from states to
    score is batched, so one call scores a whole VQE shot budget.  Verbatim from
    ``objective.FoldObjective``: this is the energy the VQE minimises and it does not move.
    """

    def __init__(self, sequence: str, rep, dist=None, mrf=None, catr=None,
                 w_dist: float = 1.0, w_mrf: float = 0.35, w_clash: float = 1.0,
                 w_catrace: float = 0.0, legacy=None, w_legacy: float = 0.0,
                 clash_ca: float = 4.1, clash_cb: float = 3.6):
        self._cons_d = None
        self._cons_w = None
        self.w_cons = 0.0
        self.seq = sequence
        self.rep = rep
        self.n = rep.n_residues
        self.dist = dist
        self.mrf = mrf
        # CA-trace (theta, tau) prior. The distance term is invariant under reflection, so
        # without this the objective literally cannot tell a structure from its mirror --
        # and 15-24% of the pool on 1LE0 is mirrored, with the mirror selected on 2 of 3
        # seeds.
        self.catr = catr
        self.w_catrace = w_catrace
        self.legacy = legacy
        self.w_legacy = float(w_legacy)
        self.w_dist, self.w_mrf, self.w_clash = w_dist, w_mrf, w_clash
        self.clash_ca, self.clash_cb = clash_ca, clash_cb
        self.rows = np.arange(self.n)
        self._iu = np.triu_indices(self.n, k=3)
        self.n_calls = 0
        self.n_structures = 0
        # `fold(workers=...)` scores restarts on several threads; these two counters are
        # the only shared mutable state on this object, so they get a lock rather than a
        # racy read-modify-write. Nothing else here is stateful per call.
        self._count_lock = threading.Lock()

    def _count(self, k: int) -> None:
        with self._count_lock:
            self.n_calls += 1
            self.n_structures += int(k)

    # -- geometry ---------------------------------------------------------
    def build(self, states: np.ndarray) -> Dict[str, np.ndarray]:
        S = np.asarray(states, int)
        return geo.build_backbone_batch(self.rep._phi[self.rows[None], S],
                                        self.rep._psi[self.rows[None], S])

    def clash(self, coords: Dict[str, np.ndarray]) -> np.ndarray:
        i, j = self._iu
        out = np.zeros(len(coords["CA"]))
        for key, cut in (("CA", self.clash_ca), ("CB", self.clash_cb)):
            if key not in coords:
                continue
            d = np.linalg.norm(coords[key][:, i] - coords[key][:, j], axis=-1)
            out += np.maximum(cut - d, 0.0).sum(1)
        return out

    def set_consensus(self, distances: Optional[np.ndarray],
                      weights: Optional[np.ndarray] = None,
                      w_cons: float = 0.0) -> None:
        """Add an L1 restraint toward a consensus distance matrix from a first pass.

        The search's own candidate pool carries distance information the prior does not:
        the pool's best member is far closer to native than the one the prior selects. A
        second pass restrained toward the pool's score-weighted mean distances is therefore
        not merely echoing the prior. It is also the obvious way to amplify the prior's
        mistakes, so `w_cons` is small and is fitted on the development peptides.
        """
        self._cons_d = None if distances is None else np.asarray(distances, float)
        self._cons_w = (None if weights is None else
                        np.asarray(weights, float) / max(np.mean(weights), 1e-12))
        self.w_cons = float(w_cons)

    def _consensus(self, ca: np.ndarray) -> np.ndarray:
        d = self.dist
        dd = np.linalg.norm(ca[:, d.i, :] - ca[:, d.j, :], axis=-1)
        err = np.abs(dd - self._cons_d[None, :])
        if self._cons_w is not None:
            err = err * self._cons_w[None, :]
        return err.mean(1)

    # -- score ------------------------------------------------------------
    def score_from(self, states: np.ndarray, coords: Dict[str, np.ndarray]) -> np.ndarray:
        e = np.zeros(len(coords["CA"]))
        if self.dist is not None and self.w_dist:
            e += self.w_dist * self.dist.score(coords["CA"])
        if self.mrf is not None and self.w_mrf:
            e += self.w_mrf * self.mrf.score(states)
        if self.w_clash:
            e += self.w_clash * self.clash(coords)
        if self.catr is not None and self.w_catrace:
            e += self.w_catrace * self.catr.score(coords["CA"])
        if self.legacy is not None and self.w_legacy:
            e += self.w_legacy * self.legacy.score_from_coords(coords, states=states)
        if self.w_cons and self._cons_d is not None:
            e += self.w_cons * self._consensus(coords["CA"])
        return e

    def __call__(self, states: np.ndarray):
        """``states`` (B, n) -> ``(scores (B,), CA (B, n, 3))``."""
        coords = self.build(states)
        self._count(len(coords["CA"]))
        return self.score_from(states, coords), coords["CA"]

    # -- continuous form, for refinement ----------------------------------
    def score_angles(self, phi: np.ndarray, psi: np.ndarray):
        coords = geo.build_backbone_batch(phi, psi)
        e = np.zeros(len(coords["CA"]))
        if self.dist is not None and self.w_dist:
            e += self.w_dist * self.dist.score(coords["CA"])
        if self.w_clash:
            e += self.w_clash * self.clash(coords)
        if self.catr is not None and self.w_catrace:
            e += self.w_catrace * self.catr.score(coords["CA"])
        if self.legacy is not None and self.w_legacy:
            e += self.w_legacy * self.legacy.score_from_coords(coords, phi=phi, psi=psi)
        if self.w_cons and self._cons_d is not None:
            e += self.w_cons * self._consensus(coords["CA"])
        self._count(len(coords["CA"]))
        return e, coords["CA"]


# ================================================================ Hamiltonian
#: Cap on the chi1 rotamer combinations scanned when scoring a *native* structure. The
#: native's chi1 is not available from the PDB (`protein_geometry` reads N, CA and C only),
#: so the fair comparison gives the native the same rotamer freedom the encoding gives a
#: prediction rather than pinning it to one rotamer -- otherwise the native is handicapped
#: by a degree of freedom the prediction gets to optimise. 2**6 = 64 evaluations is a
#: generous ceiling; beyond it the scan falls back to rotamer 0.
MAX_NATIVE_CHI_SCAN_BITS = 6


class FoldingHamiltonian(BudgetedEnergyModel):
    """Knowledge-based Hamiltonian over a discrete structure encoding.

    Verbatim from ``hamiltonian.FoldingHamiltonian``.  This is the energy the global
    CVaR-VQE minimises, it owns the shared evaluation budget, and it does not move.
    """

    def __init__(self, sequence: str, representation,
                 weights: Optional[Dict[str, float]] = None,
                 use_corrected_mj: bool = True,
                 backtracking_penalty: float = 5.0,
                 cache_limit: int = 500_000,
                 eval_budget: Optional[int] = None):
        self.sequence = sequence.strip().upper()
        self.rep = representation
        if len(self.sequence) != self.rep.n_residues:
            raise ValueError(
                f"sequence length {len(self.sequence)} != representation "
                f"n_residues {self.rep.n_residues}")
        self.weights = dict(et.DEFAULT_WEIGHTS if weights is None else weights)
        self.use_corrected_mj = bool(use_corrected_mj)
        self.backtracking_penalty = float(backtracking_penalty)
        self._init_budget(cache_limit, eval_budget)

    @property
    def n_qubits(self) -> int:
        return self.rep.n_qubits

    @property
    def n_bits(self) -> int:
        return self.rep.n_bits

    # -- aromatic rings ---------------------------------------------------
    def _ring_residues(self):
        """Residue indices whose aromatic ring can be built, with their 3-letter code."""
        return [(i, geo.ONE_TO_THREE[aa]) for i, aa in enumerate(self.sequence)
                if geo.ONE_TO_THREE[aa] in sc.AROMATIC_RING_RESIDUES]

    def _rings_from_coords(self, coords: Dict[str, np.ndarray],
                           chi1: Optional[Dict[int, float]] = None
                           ) -> Dict[int, Dict[str, np.ndarray]]:
        """Aromatic ring atoms from a full backbone, at the given chi1.

        Used for structures that did not come from a bitstring -- natives, mainly -- so
        their aromatic term is computed with the *same* functional form as a prediction's.
        Without this the native would silently fall back to the CB proxy while predictions
        used real ring geometry, and the two energies would not be comparable at all.
        """
        if not all(k in coords for k in ("N", "CA", "C", "CB")):
            return {}
        chi1 = chi1 or {}
        out: Dict[int, Dict[str, np.ndarray]] = {}
        for i, key in self._ring_residues():
            names = sc.ring_atom_names(key)
            if not names:
                continue
            atoms = sc.build_sidechain(key, coords["N"][i], coords["CA"][i],
                                       coords["C"][i], coords["CB"][i],
                                       chi1=chi1.get(i))
            out[i] = {n: atoms[n] for n in names if n in atoms}
        return out

    # -- energy -----------------------------------------------------------
    def components(self, bitstring: str) -> Dict[str, float]:
        # One pass through the representation: validate once, decode once, build the
        # backbone once, and build the aromatic rings from that same backbone. Calling
        # decode() and build_coords() separately validated the bitstring four times and
        # decoded it twice, on a path that runs 20,000 times per search arm.
        phi, psi, coords, rings = self.rep.build_all(bitstring)
        comp = et.energy_components(self.sequence, coords, phi, psi,
                                    use_corrected_mj=self.use_corrected_mj,
                                    rings=rings)
        comp["backtracking"] = (self.backtracking_penalty
                                * et.backtracking_term(self.rep, bitstring))
        return comp

    def energy(self, bitstring: str) -> float:
        hit = self._cache.get(bitstring)
        if hit is not None:
            return hit                  # cache hits are free -- revisiting costs nothing
        self._charge()                  # raises BudgetExhausted before doing the work
        comp = self.components(bitstring)
        e = et.total_from_components(comp, self.weights) + comp["backtracking"]
        if len(self._cache) < self._cache_limit:
            self._cache[bitstring] = e
        return e

    def energy_from_coords(self, coords: Dict[str, np.ndarray],
                           phi=None, psi=None,
                           chi1: Optional[Dict[int, float]] = None,
                           scan_chi: bool = True) -> float:
        """Energy of an arbitrary structure, typically a native.

        `scan_chi=True` minimises over the encoded chi1 rotamers, giving this structure the
        same sidechain freedom a bitstring-encoded prediction has. Pass an explicit `chi1`,
        or `scan_chi=False`, to pin the rotamers instead.
        """
        def evaluate(chi):
            comp = et.energy_components(
                self.sequence, coords, phi, psi,
                use_corrected_mj=self.use_corrected_mj,
                rings=self._rings_from_coords(coords, chi) or None)
            return et.total_from_components(comp, self.weights)

        if chi1 is not None or not scan_chi:
            return evaluate(chi1)

        residues = [i for i, _ in self._ring_residues()]
        if not residues or len(residues) > MAX_NATIVE_CHI_SCAN_BITS:
            return evaluate(None)
        rot = {i: sc.CHI1_ROTAMERS[geo.ONE_TO_THREE[self.sequence[i]]]
               for i in residues}
        best = float("inf")
        for combo in itertools.product(*(range(len(rot[i])) for i in residues)):
            best = min(best, evaluate({i: rot[i][k]
                                       for i, k in zip(residues, combo)}))
        return best


# ============================================================== folding driver
# The redesigned pipeline: prior -> warm-started CVaR-VQE -> refine -> select.
#
#     sequence
#        |
#        +-- peptide_db.holdout ------------------------------------------+
#        |      (target and everything >0.6 identity to it removed)       |
#        v                                                                |
#     per-residue torsion library (torsion_lib2)  <---------- context stats
#        |                                                                |
#        +--> torsion MRF prior  --------+                                |
#        |        (h_i, J_i,i+1)         |                                |
#        +--> ESM-2 + sequence distogram-+                                |
#                 |                      |                                |
#                 v                      v                                |
#           FoldObjective         warm-start angles ----+                 |
#                 |                                     |                 |
#                 v                                     v                 |
#         CVaR-VQE  (exact-sampled RY/CNOT ansatz, analytic CVaR gradient) |
#                 |                                                       |
#                 +--> diverse reservoir of candidates                    |
#                              |                                          |
#                              v                                          |
#                   continuous torsion refinement                         |
#                              |                                          |
#                              v                                          |
#               basin selection (population + score)                      |
#                              |                                          |
#                              v                                          |
#               optional physics rescoring: Model A / Amber ff14SB+GB ----+
#
# *Warm start.* The RY angles of the one-layer ansatz set the per-qubit marginals exactly,
# and the marginals of a torsion-state prior are exactly what we know before searching.
# The inversion is closed form. The previous driver initialised at pi/2 plus noise -- the
# uniform distribution -- and spent its budget rediscovering the Ramachandran statistics of
# a 780-peptide database.
#
# *Selection by basin, not by score.* The scorer's own ranking of its top candidates is
# weak; the population of a basin under a decent scorer is a much better signal, and it is
# what guards the worst case.
def marginals_to_angles(bit_marginals: np.ndarray, ring: bool = False) -> np.ndarray:
    """RY angles whose one-layer ansatz reproduces the given per-qubit marginals.

    ``b_q`` is the prefix XOR of independent Bernoulli(s_t), so
    ``1 - 2 P(b_q = 1) = prod_{t <= q} (1 - 2 s_t)``.  Dividing consecutive running products
    recovers each ``s_q`` exactly.  The ring closure is ignored here: it couples the first
    and last qubits, and matching marginals up to that one coupling is enough for an
    initialisation.
    """
    m = np.clip(np.asarray(bit_marginals, float), 0.02, 0.98)
    A = 1.0 - 2.0 * m
    prev = np.concatenate([[1.0], A[:-1]])
    ratio = np.clip(A / prev, -0.999, 0.999)
    s = np.clip((1.0 - ratio) / 2.0, 1e-3, 1 - 1e-3)
    return 2.0 * np.arcsin(np.sqrt(s))


def state_marginals_to_bits(pi: np.ndarray, bits_per_residue: int) -> np.ndarray:
    """``pi`` (n, k) state probabilities -> (n * bits_per_residue,) qubit marginals."""
    n, k = pi.shape
    codes = ((np.arange(k)[:, None] >> np.arange(bits_per_residue - 1, -1, -1)) & 1)
    return (pi @ codes).reshape(-1)


def make_ansatz(kind: str, n_qubits: int, layers: int, rep, mrf, warm_start: bool):
    """``(ansatz, initial angles)``.

    ``chain``     the shipped one-layer RY + CNOT chain + ring, sampled in closed form.
    ``chain_ry``  the same chain with one trailing RY layer, which is what makes the
                  prior's marginals reachable at all (see `MPSAnsatz`).
    ``deep``      ``layers`` chain blocks plus a trailing RY layer, exact MPS.
    ``product``   RY only, no entangler -- the control for whether entanglement helps.
    """
    m = (state_marginals_to_bits(mrf.marginals(), rep.bits_per_residue)
         if warm_start else None)
    warm = marginals_to_angles(m) if warm_start else np.full(n_qubits, np.pi / 2.0)
    if kind == "chain":
        return OneLayerAnsatz(n_qubits, ring=True), warm
    if kind == "product":
        return MPSAnsatz(n_qubits, 1, final_ry=False, entangler="none"), warm
    nblocks = 1 if kind == "chain_ry" else max(1, layers)
    a = MPSAnsatz(n_qubits, nblocks, final_ry=True)
    # Warm start goes in the LAST rotation layer with the entangling blocks at identity:
    # RY(0) is the identity and a CNOT on |0...0> does nothing, so the prepared state is
    # exactly the product distribution with the prior's marginals. Putting it in the first
    # layer instead would push it through the CNOT chain and lose it.
    base = np.concatenate([np.zeros(n_qubits * nblocks), warm])
    return a, base


class Reservoir:
    """Score-ordered, diversity-preserving candidate store."""

    def __init__(self, radius: float = 1.2, capacity: int = 64):
        self.radius = radius
        self.capacity = capacity
        self.scores: List[float] = []
        self.cas: List[np.ndarray] = []
        self.states: List[np.ndarray] = []

    def offer_batch(self, scores: np.ndarray, cas: np.ndarray,
                    states: np.ndarray) -> None:
        for s, c, st in zip(scores, cas, states):
            self._offer(float(s), c, st)

    def _offer(self, score: float, ca: np.ndarray, st: np.ndarray) -> None:
        if not self.cas:
            self.scores.append(score)
            self.cas.append(ca)
            self.states.append(st)
            return
        d = geo.ca_rmsd_batch(np.stack(self.cas), ca)
        j = int(np.argmin(d))
        if d[j] > self.radius:
            self.scores.append(score)
            self.cas.append(ca)
            self.states.append(st)
            if len(self.scores) > self.capacity:
                worst = int(np.argmax(self.scores))
                for lst in (self.scores, self.cas, self.states):
                    lst.pop(worst)
        elif score < self.scores[j]:
            self.scores[j], self.cas[j], self.states[j] = score, ca, st

    def items(self):
        return list(zip(self.scores, self.cas, self.states))


class _OfferLog:
    """Records reservoir offers instead of applying them.

    The Reservoir is order-dependent (an offer inside the diversity radius overwrites its
    neighbour), so concurrent restarts cannot share one.  They do not need to: a restart's
    offers depend only on that restart's own sampling, so each thread records its offers
    and `fold` replays them into the shared reservoir in restart order.  The result is the
    same sequence of `_offer` calls the sequential run made, hence a bit-identical
    reservoir -- which is what `t_fold_workers_reservoir_identical` asserts.
    """

    __slots__ = ("records",)

    def __init__(self):
        self.records = []

    def offer_batch(self, scores, cas, states):
        self.records.append((np.array(scores, float), np.array(cas, float),
                             np.array(states)))

    def replay(self, reservoir: Reservoir) -> None:
        for s, c, st in self.records:
            reservoir.offer_batch(s, c, st)


def refine(objective, states: np.ndarray, sigma0: float = 12.0, steps: int = 40,
           pop: int = 24, seed: int = 0, restarts: int = 2
           ) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Continuous local relaxation of (phi, psi) away from the discrete grid.

    A (1 + pop) evolution strategy in torsion space with a geometrically shrinking step,
    which reuses the batched builder the search already runs at ~4,000 structures per
    second.  The discrete grid costs ~0.5-1.5 A of representation floor; this is what
    recovers the part of that gap the objective can actually see.
    """
    rows = np.arange(objective.n)
    phi0 = objective.rep._phi[rows, states].copy()
    psi0 = objective.rep._psi[rows, states].copy()
    e0, ca0 = objective.score_angles(phi0[None], psi0[None])
    gbest, gphi, gpsi, gca = float(e0[0]), phi0, psi0, ca0[0]
    # Several short anneals from the same discrete start beat one long one: the surface off
    # the grid is rough enough that a single (1+lambda) run stalls at the first basin it
    # shrinks into.
    for r in range(max(1, restarts)):
        rng = np.random.default_rng(seed * 977 + r)
        phi, psi = phi0.copy(), psi0.copy()
        best, best_ca = float(e0[0]), ca0[0]
        sigma = np.radians(sigma0)
        for t in range(steps):
            cand_phi = phi[None] + rng.normal(0, sigma, (pop, objective.n))
            cand_psi = psi[None] + rng.normal(0, sigma, (pop, objective.n))
            e, cas = objective.score_angles(cand_phi, cand_psi)
            k = int(np.argmin(e))
            if e[k] < best:
                best, phi, psi, best_ca = float(e[k]), cand_phi[k], cand_psi[k], cas[k]
            else:
                sigma *= 0.88
            if sigma < np.radians(0.4):
                break
        if best < gbest:
            gbest, gphi, gpsi, gca = best, phi, psi, best_ca
    return gphi, gpsi, gbest, gca


def snap_to_states(rep, ca: np.ndarray, seed: int = 0, restarts: int = 6):
    """Torsion states whose CA trace best matches `ca`, and that trace."""
    import floor
    rng = np.random.default_rng(seed)
    best, best_v = None, float("inf")
    for r in range(restarts):
        start = (np.zeros(rep.n_residues, int) if r == 0
                 else rng.integers(0, rep.n_states, rep.n_residues))
        st, v = floor.descend(rep, np.asarray(ca, float), start)
        if v < best_v:
            best, best_v = st, v
    rows = np.arange(rep.n_residues)
    out = geo.build_backbone_batch(rep._phi[rows[None], best[None]],
                                   rep._psi[rows[None], best[None]])["CA"][0]
    return best, out


def basin_average(cas: Sequence[np.ndarray], weights: Sequence[float],
                  members: Sequence[int], centre: int) -> np.ndarray:
    """Boltzmann-weighted mean of a basin's members, superposed on its medoid.

    Each member carries independent error; averaging them cancels the part of it that is
    random rather than systematic.  The result is not a lattice structure -- it does not
    correspond to any torsion assignment -- so it is only ever the *reported* answer, never
    a search candidate.
    """
    ref = np.asarray(cas[centre], float)
    w = np.asarray([weights[i] for i in members], float)
    w = w / max(w.sum(), 1e-12)
    stack = np.stack([geo.kabsch_superpose(np.asarray(cas[i], float), ref)
                      for i in members])
    return (stack * w[:, None, None]).sum(0)


def basin_select(scores: Sequence[float], cas: Sequence[np.ndarray],
                 radius: float = 2.0, kT: float = 1.0) -> Tuple[int, Dict]:
    """Medoid of the basin with the most Boltzmann-weighted population.

    Population, not best score: on this objective the top candidates cannot be ranked
    against each other reliably, but a basin that many independent restarts fall into is
    much more often the right one.  `kT` is in units of the objective, whose scale is set
    by the distogram term (mean per-pair distance error in Angstroms).
    """
    C = np.stack(cas)
    s = np.asarray(scores, float)
    D = np.stack([geo.ca_rmsd_batch(C, c) for c in C])
    w = np.exp(-(s - s.min()) / kT)
    mass = (D < radius) @ w
    members = np.where(D[int(np.argmax(mass))] < radius)[0]
    sub = D[np.ix_(members, members)]
    pick = int(members[np.argmin((sub * w[members][None, :]).sum(1))])
    return pick, {"n_basins": int(len(set(map(tuple, (D < radius).astype(int))))),
                  "basin_size": int(len(members)),
                  "basin_mass": float(mass.max() / w.sum()),
                  "members": members.tolist(), "weights": w.tolist()}


def build_distogram(sequence: str, rep, use_esm: bool = True, fragments: bool = True,
                    models: str = "auto", n_ref: int = 512, seed: int = 0):
    """The distance prior for a target: one model, or both combined and calibrated.

    ``models="both"`` averages the per-pair MLP and the pair-tensor network after
    standardising each against a fixed sample of random structures from this target's own
    state library.  The sample is drawn from the representation, never from the native.
    """
    import distogram as dgm
    if models == "auto":
        # "both" is better (in-band +0.454 against +0.451 and +0.397 alone), but only when
        # the pair-network fold models actually exist. Training one on demand inside a
        # target build is a 25-minute stall, so availability is checked rather than assumed.
        import glob
        import pairnet as _pn
        models = "both" if glob.glob(os.path.join(_pn.MODEL_DIR, "fold*.pt")) else "mlp"
    parts = []
    if models in ("mlp", "both"):
        parts.append(dgm.Distogram.for_target(sequence, use_esm=use_esm,
                                              fragments=fragments))
    if models in ("pairnet", "both"):
        try:
            import pairnet
            parts.append(pairnet.distogram_for(sequence))
        except Exception:
            pass
    if len(parts) == 1:
        return parts[0]
    rng = np.random.default_rng(seed)
    S = rng.integers(0, rep.n_states, size=(n_ref, rep.n_residues))
    rows = np.arange(rep.n_residues)
    ref = geo.build_backbone_batch(rep._phi[rows[None], S],
                                   rep._psi[rows[None], S])["CA"]
    return dgm.CombinedDistogram(parts, ref_ca=ref)


def build_target(sequence: str, k: int = 8, use_esm: bool = True,
                 exclude_self: bool = True, w_mrf: float = 0.35,
                 w_clash: float = 1.0, fragments: bool = True,
                 models: str = "auto", w_catrace: float = 0.0,
                 w_legacy: float = 0.0,
                 legacy_weights: Optional[Dict[str, float]] = None,
                 legacy_use: Optional[Sequence[str]] = None) -> Dict:
    """Everything that depends only on the sequence: library, priors, objective."""
    import catrace
    import peptide_db as pdb
    import priors
    import torsion_lib2 as tl2
    ex = sequence if exclude_self else ""
    entries = pdb.holdout(sequence) if exclude_self else list(pdb.load())
    tab = tl2.library_for(sequence, k, ex)
    rep = tl2.PerResidueTorsion(sequence, tab, chi_bits=False)
    dist = build_distogram(sequence, rep, use_esm=use_esm, fragments=fragments,
                           models=models)
    mrf = priors.TorsionMRF.fit(sequence, tab, entries)
    catr = catrace.CATracePrior(exclude_seq=ex)
    leg = None
    if w_legacy:
        import legacy_field
        leg = legacy_field.LegacyField(sequence, rep, weights=legacy_weights,
                                       use=legacy_use, dist=dist)
    o = FoldObjective(sequence, rep, dist=dist, mrf=mrf, catr=catr,
                      w_mrf=w_mrf, w_clash=w_clash, w_catrace=w_catrace,
                      legacy=leg, w_legacy=w_legacy)
    return {"rep": rep, "dist": dist, "mrf": mrf, "catrace": catr, "legacy": leg,
            "objective": o, "table": tab}


def seed_states(target: Dict, restarts: int = 4, seed: int = 0) -> np.ndarray:
    """Torsion states whose CA trace best matches the distogram's own embedded structure.

    Distance geometry on the *predicted* matrix gives a coordinate trace; snapping that
    trace onto the state library is the discrete assignment the search should start from.
    Nothing here touches the native -- the target's structure is not read at any point --
    so this is a prior-driven initialisation, not leakage.
    """
    import floor
    rep = target["rep"]
    X = target["dist"].realize(seed=seed, restarts=restarts)
    rng = np.random.default_rng(seed)
    best, best_v = None, float("inf")
    for r in range(6):
        start = (np.zeros(rep.n_residues, int) if r == 0
                 else rng.integers(0, rep.n_states, rep.n_residues))
        st, v = floor.descend(rep, X, start)
        if v < best_v:
            best, best_v = st, v
    return best


def run_vqe(objective, ansatz, theta0: np.ndarray, iters: int, shots: int,
            rng, lr: float = 0.12, a0: float = 0.5, a1: float = 0.05,
            reservoir=None, elite_frac: float = 0.1,
            random_frac: float = 0.0, update: str = "cvar_grad",
            baseline: str = "const", trace: Optional[List] = None) -> Dict:
    """One CVaR-VQE restart. Returns the final angles and the best structure seen.

    ``update="cvar_grad"`` (default) takes a genuine score-function gradient of the CVaR
    objective; ``update="moment"`` keeps the shipped driver's cross-entropy update
    available so the comparison stays measurable -- it is not a gradient of any objective
    and it is not the default.

    THE ONE DEFAULT THAT CHANGED, STATED PLAINLY
    --------------------------------------------
    ``baseline`` defaults to ``"const"``, the CORRECTED control variate, where
    `foldvqe.run_vqe` used the biased tail-only form unconditionally.  So this function
    does NOT reproduce `foldvqe.run_vqe` at its defaults, by design: the shipped estimator
    points ~49 degrees off the true CVaR gradient and is 24% short in magnitude even with
    infinite shots, and carrying that forward to keep a benchmark bit-identical would be
    the wrong trade.  Pass ``baseline="tail"`` for exact reproduction -- the test suite
    does, and gets a bit-identical trajectory on the chain ansatz.

    Nothing downstream is silently affected: `core.pipeline` and `core.bench` do not use
    the quantum backend, and `bench.py`, `validation.py` and `refine2.py` import the root
    `foldvqe` directly rather than through `core.backend`.

    Shots, iterations, the alpha anneal, the Adam settings and the rng stream are exactly
    the shipped ones.
    """
    theta = np.array(theta0, float)
    adam = Adam(len(theta), lr=lr)
    rep = objective.rep
    w, n_res, k = rep.bits_per_residue, rep.n_residues, rep.n_states
    powers = (1 << np.arange(w - 1, -1, -1)).astype(np.int64)
    n_elite = max(1, int(elite_frac * shots))
    n_extra = max(1, int(random_frac * shots)) if random_frac > 0 else 0
    best_e, best_ca, best_st = float("inf"), None, None
    for it in range(iters):
        alpha = alpha_schedule(it / max(1, iters - 1), a0, a1)
        bits = ansatz.sample(theta, shots, rng)
        states = (bits[:, :w * n_res].reshape(shots, n_res, w) @ powers) % k
        e, ca = objective(states)
        m = int(np.argmin(e))
        if e[m] < best_e:
            best_e, best_ca, best_st = float(e[m]), ca[m], states[m]
        if reservoir is not None:
            sel = np.argsort(e)[:n_elite]
            if n_extra:
                # Admission by SCORE is admission by a weak ranker. Instrumentation shows
                # the best structure the objective ever scores averages 1.06 A while the
                # best one that survives the elite filter averages 1.23 A -- the filter
                # itself throws away 0.17 A using a score whose in-band Spearman is 0.45.
                # Offering a score-blind random sample as well lets the diversity radius,
                # rather than the score, decide what the pool covers.
                extra = rng.choice(shots, size=n_extra, replace=False)
                sel = np.unique(np.concatenate([sel, extra]))
            reservoir.offer_batch(e[sel], ca[sel], states[sel])
        if update == "cvar_grad":
            grad, val = cvar_gradient(ansatz, theta, bits, e, alpha, baseline=baseline)
            theta = adam.step(theta, grad)
        else:
            # The shipped driver's update, kept so the comparison is measurable: pull the
            # angles toward the mean bit pattern of the CVaR tail. Not a gradient of any
            # objective -- it is the cross-entropy method.
            val = cvar(e, alpha)[0]
            tail = np.argsort(e)[:max(1, int(np.ceil(alpha * shots)))]
            frac = np.clip(bits[tail].mean(0), 1e-6, 1 - 1e-6)
            tgt = 2.0 * np.arcsin(np.sqrt(frac))
            step_lr = lr * 2.9 * (1.0 - it / max(1, iters - 1))
            theta = theta + step_lr * (tgt - theta) + rng.normal(0, 0.05, theta.shape)
        if trace is not None:
            trace.append({"iter": it, "alpha": alpha, "cvar": val,
                          "min": float(e.min()), "mean": float(e.mean()),
                          "best": best_e})
    return {"theta": theta, "best_energy": best_e, "best_ca": best_ca,
            "best_states": best_st}


def fold(sequence: str, k: int = 8, layers: int = 1, restarts: int = 4,
         iters: int = 60, shots: int = 384, seed: int = 0, use_esm: bool = True,
         warm_start: bool = True, refine_top: int = 8, radius: float = 1.2,
         fragments: bool = True, select_radius: float = 2.0,
         target: Optional[Dict] = None, exclude_self: bool = True, lr: float = 0.12,
         alpha0: float = 0.5, alpha1: float = 0.05, update: str = "cvar_grad",
         baseline: str = "const", refine_kw: Optional[Dict] = None,
         random_frac: float = 0.0, reservoir_capacity: int = 64,
         select: str = "basin", ansatz_kind: str = "chain",
         seed_frac: float = 0.6, seed_restarts: int = 1, use_seed: bool = True,
         workers: int = 1, verbose: bool = False) -> Dict:
    """Fold `sequence`.

    ``exclude_self=True`` is the benchmark setting: the target's own structure and every
    homolog above 0.6 identity are removed from every fitted component.

    ``workers > 1`` runs the VQE restarts on that many threads.  They are genuinely
    independent -- different rng streams, no shared optimiser state -- and their reservoir
    offers are recorded per restart and replayed in restart order, so the pool is
    bit-identical to the sequential run (`test_reservoir_replay_is_identical_to_sequential
    _offering`).

    IT DEFAULTS TO 1 BECAUSE IT WAS MEASURED AND IT DOES NOT PAY.  Four restarts, 30
    iterations, 384 shots, 20 residues at k=8, against the real batched backbone builder,
    wall-clock in ms:

        ansatz     sequential   2 threads   4 threads   8 threads
        chain           2592       3090        3762        3782
        chain_ry        4274       3812        4567        4527

    Free parallelism is only free when the work releases the GIL, and this work does not:
    at 60 qubits and 384 shots every array in the inner loop is a few kilobytes, so the
    cost is numpy and torch DISPATCH -- Python bytecode holding the interpreter lock -- and
    not the BLAS calls that would let threads overlap.  The one arm that gains anything is
    the MPS ansatz at 2 threads (1.12x), because torch's backward pass does release the
    lock, and even that is inside the noise of a box running three sibling jobs.

    The mechanism is kept, defaulted off, because it is proven exactly equivalent and the
    conclusion is a property of THIS objective: an Amber-scored or OpenMM-scored objective
    spends its time outside the interpreter and would scale.  The honest summary is that
    the restarts are parallelisable and the profile says there is nothing to parallelise.
    """
    t_start = time.perf_counter()
    T = target if target is not None else build_target(
        sequence, k=k, use_esm=use_esm, exclude_self=exclude_self, fragments=fragments)
    rep, o, mrf = T["rep"], T["objective"], T["mrf"]
    o.n_calls = o.n_structures = 0      # per-fold accounting, not cumulative
    n_qubits = rep.bits_per_residue * rep.n_residues

    # Distance-geometry seed: embed the predicted distogram, snap the trace onto the state
    # library, and mix that assignment into the warm-start marginals. The VQE still does
    # the search -- this only says where to start looking, and it is derived entirely from
    # the prediction, never from the target's structure.
    seed_st = None
    if use_seed:
        try:
            seed_st = seed_states(T, seed=seed)
        except Exception:
            seed_st = None

    ansatz, base = make_ansatz(ansatz_kind, n_qubits, layers, rep, mrf, warm_start)
    seeded_base = base
    if seed_st is not None and seed_frac > 0 and warm_start:
        import priors
        pi = mrf.marginals()
        onehot = np.zeros_like(pi)
        onehot[np.arange(len(seed_st)), seed_st] = 1.0
        blended = (1.0 - seed_frac) * pi + seed_frac * onehot
        seeded_base = make_ansatz(
            ansatz_kind, n_qubits, layers, rep,
            priors.TorsionMRF(-np.log(np.clip(blended, 1e-6, None)), mrf.J), True)[1]

    res = Reservoir(radius=radius, capacity=reservoir_capacity)
    if seed_st is not None:
        e_seed, ca_seed = o(seed_st[None])
        res.offer_batch(e_seed, ca_seed, seed_st[None])

    # The first `seed_restarts` restarts start from the distance-geometry seed's basin; the
    # rest start from the sequence prior. Seeding every restart measurably narrows the
    # search (it wins on some targets and loses on more), so the two initialisations are
    # run side by side inside one fold instead of chosen.
    def _one(r: int, sink):
        rng = np.random.default_rng([seed, r])
        start = seeded_base if r < seed_restarts else base
        theta0 = start + rng.normal(0, 0.35, size=ansatz.n_params())
        tr: List = []
        run = run_vqe(o, ansatz, theta0, iters, shots, rng, lr=lr, a0=alpha0, a1=alpha1,
                      reservoir=sink, update=update, baseline=baseline,
                      random_frac=random_frac, trace=tr)
        return run, tr

    traces: List = [None] * restarts
    runs: List = [None] * restarts
    if int(workers) > 1 and restarts > 1:
        from concurrent.futures import ThreadPoolExecutor
        logs = [_OfferLog() for _ in range(restarts)]
        with ThreadPoolExecutor(max_workers=int(workers)) as pool:
            for r, (run, tr) in enumerate(pool.map(lambda i: _one(i, logs[i]),
                                                   range(restarts))):
                runs[r], traces[r] = run, tr
        for lg in logs:                          # replay in restart order
            lg.replay(res)
        if verbose:
            for r, run in enumerate(runs):
                print(f"    restart {r}: best {run['best_energy']:.4f}", flush=True)
    else:
        for r in range(restarts):
            run, tr = _one(r, res)
            runs[r], traces[r] = run, tr
            if verbose:
                print(f"    restart {r}: best {run['best_energy']:.4f}, "
                      f"reservoir {len(res.scores)}", flush=True)

    items = res.items()
    scores = [s for s, _, _ in items]
    cas = [c for _, c, _ in items]
    states = [st for _, _, st in items]
    # Torsion angles of every pool member, kept so a *resolution-independent* rescore can
    # use terms that need more than the CA trace. A refined candidate no longer sits on its
    # state library, so its states are not a description of it; its angles are.
    _rows = np.arange(rep.n_residues)
    _S = np.asarray(states, int)
    phis = list(rep._phi[_rows[None, :], _S]) if len(states) else []
    psis = list(rep._psi[_rows[None, :], _S]) if len(states) else []

    # Refinement. Two corrections over the previous version, both from measurement:
    #
    # 1. The old code compared a refined `dist + clash` value against an unrefined
    #    `dist + MRF + clash` one, so a refinement was accepted partly because the MRF term
    #    silently vanished. That mismatch was acting as an accidental conservative accept
    #    filter -- putting both on one basis measures 0.14 A WORSE. Rather than keep an
    #    accident, the refined structure is now ADDED to the pool alongside the unrefined
    #    one and everything is re-scored on the same basis, so selection sees both and no
    #    hidden filter decides for it.
    # 2. `refine2.refine_pool` relaxes the whole reservoir with an exact analytic torsion
    #    gradient at 0.33 s/structure -- 2.7x cheaper than the (1+lambda) search it
    #    replaces, and it reaches a lower objective. Accuracy is the same, because the
    #    objective, not the operator, is the constraint: across nine refinement settings the
    #    within-target Spearman between objective reached and CA-RMSD is +0.014.
    refined = None
    if refine_top:
        order = list(np.argsort(scores)[:refine_top])
        new_cas, new_states, new_scores = [], [], []
        new_phi, new_psi = [], []
        try:
            import refine2
            phi_r, psi_r, e_r, ca_r = refine2.refine_pool(
                o, [states[i] for i in order], seed=seed, **(refine_kw or {}))
            for m, idx in enumerate(order):
                new_cas.append(ca_r[m])
                new_states.append(states[idx])
                new_scores.append(float(e_r[m]))
                new_phi.append(phi_r[m])
                new_psi.append(psi_r[m])
        except Exception:
            for idx in order:
                phi, psi, e, ca_r = refine(o, states[idx], seed=seed + int(idx))
                new_cas.append(ca_r)
                new_states.append(states[idx])
                new_scores.append(float(e))
                new_phi.append(phi)
                new_psi.append(psi)
        # Re-score refined structures on the same basis the unrefined ones carry, so the
        # two populations are comparable.
        if new_cas:
            base_mrf = (o.w_mrf * o.mrf.score(np.stack(new_states))
                        if (o.mrf is not None and o.w_mrf) else 0.0)
            cas.extend(new_cas)
            states.extend(new_states)
            scores.extend(list(np.asarray(new_scores) + base_mrf))
            phis.extend(new_phi)
            psis.extend(new_psi)
        refined = int(len(order))

    if select == "bag":
        # Rank the pool by agreement across random pair subsets, then take the medoid of
        # the best-ranked basin. Combines the two selection signals that were measured to
        # matter: basin population and robustness of the score to which pairs it uses.
        rk = T["dist"].bagged_rank(np.stack(cas), n_bags=16, seed=seed)
        pick, info = basin_select(list(rk), cas, radius=select_radius, kT=6.0)
        info["mode"] = "bag"
    elif select == "best":
        pick, info = int(np.argmin(scores)), {"mode": "best_score"}
    elif select == "consensus":
        C = np.stack(cas)
        D = np.stack([geo.ca_rmsd_batch(C, c) for c in C])
        pick, info = int(np.argmin(D.mean(1))), {"mode": "consensus"}
    else:
        pick, info = basin_select(scores, cas, radius=select_radius)
        info["mode"] = "basin"
    wall = time.perf_counter() - t_start
    rows = np.arange(rep.n_residues)
    _st = np.asarray(states[pick], int)
    sel_coords = geo.build_backbone_batch(rep._phi[rows[None], _st[None]],
                                          rep._psi[rows[None], _st[None]])
    sel_coords = {key: v[0] for key, v in sel_coords.items()}
    return {"ca": cas[pick], "score": scores[pick], "states": states[pick],
            "coords": sel_coords,
            "dist": T["dist"], "objective": o, "rep": rep,
            "pool_k": [k] * len(cas), "reps": {k: rep},
            "clash_of": (lambda C: o.w_clash * o.clash({"CA": np.asarray(C, float)})),
            "pool_scores": scores, "pool_cas": cas, "pool_states": states,
            "pool_phi": phis, "pool_psi": psis,
            "best_score_index": int(np.argmin(scores)),
            "selection": info, "n_candidates": len(items),
            "n_structures": o.n_structures, "wall": wall,
            "structures_per_s": o.n_structures / max(wall, 1e-9),
            "qubits": n_qubits, "layers": layers, "k": k,
            "traces": traces, "refined": refined,
            "vqe_best_energy": float(min(r["best_energy"] for r in runs))}


def fold_multi(sequence: str, ks: Sequence[int] = (8, 16), seed: int = 0,
               targets: Optional[Dict[int, Dict]] = None,
               select_radius: float = 2.0,
               w_legacy_select: Optional[float] = None, **kw) -> Dict:
    """Run the pipeline at several state resolutions and select across the union.

    The ablation says k=8 and k=16 are not ordered: k=16 gives a better median and more
    sub-2 A results (its representation floor is 0.69 A against 0.93 A) and a slightly
    worse mean, because the larger register is harder for the same shot budget.  Running
    both and selecting over the merged reservoir takes the better of the two per target
    without having to pick in advance.

    Candidates are re-scored on a resolution-independent basis before selection -- the
    distance prior plus the clash term, both functions of the CA trace alone -- because the
    MRF term's scale depends on the library and would otherwise bias the merge toward
    whichever k produced it.
    """
    pools_ca: List[np.ndarray] = []
    pools_st: List[np.ndarray] = []
    pools_k: List[int] = []
    pools_phi: List[np.ndarray] = []
    pools_psi: List[np.ndarray] = []
    per_k: Dict[int, Dict] = {}
    ref = None
    for k in ks:
        T = (targets or {}).get(k)
        r = fold(sequence, k=k, seed=seed, target=T, select_radius=select_radius, **kw)
        per_k[k] = r
        pools_ca.extend(r["pool_cas"])
        pools_st.extend(r["pool_states"])
        pools_phi.extend(r.get("pool_phi", []))
        pools_psi.extend(r.get("pool_psi", []))
        pools_k.extend([k] * len(r["pool_cas"]))
        if ref is None:
            ref = r
    o = ref["objective"] if "objective" in ref else None
    dist = per_k[ks[0]]["dist"]
    C = np.stack(pools_ca)
    scores = list(dist.score(C) + ref["clash_of"](C))
    # The Legacy field, if the search used one, is added on the same resolution-independent
    # basis. It is a function of the torsion ANGLES, not of the state library, so it merges
    # across k exactly as the distance term does -- and unlike the distance term it can see
    # backbone hydrogen bonding, which is why the merged pool can contain a structure the
    # distance score alone has no way to prefer.
    #
    # `w_legacy_select` defaults to the weight the SEARCH used. Setting it to 0 separates
    # the two effects cleanly: the pool is whatever generation produced either way, and only
    # the pick changes.
    leg = getattr(o, "legacy", None) if o is not None else None
    wl = (getattr(o, "w_legacy", 0.0) if w_legacy_select is None
          else float(w_legacy_select))
    if leg is not None and wl and len(pools_phi) == len(pools_ca):
        PHI, PSI = np.stack(pools_phi), np.stack(pools_psi)
        coords = geo.build_backbone_batch(PHI, PSI)
        scores = list(np.asarray(scores)
                      + wl * leg.score_from_coords(coords, phi=PHI, psi=PSI))
    pick, info = basin_select(scores, pools_ca, radius=select_radius)
    wall = sum(r["wall"] for r in per_k.values())
    out = dict(per_k[ks[0]])
    out.update({"ca": pools_ca[pick], "score": scores[pick], "states": pools_st[pick],
                "pool_scores": scores, "pool_cas": pools_ca, "pool_states": pools_st,
                "pool_phi": pools_phi, "pool_psi": pools_psi,
                "pool_k": pools_k, "reps": {k: per_k[k]["rep"] for k in ks},
                "n_candidates": len(pools_ca), "selection": info, "wall": wall,
                "n_structures": sum(r["n_structures"] for r in per_k.values()),
                "structures_per_s": (sum(r["n_structures"] for r in per_k.values())
                                     / max(wall, 1e-9)),
                "ks": list(ks)})
    return out


def consensus_distances(pool_cas: Sequence[np.ndarray], pool_scores: Sequence[float],
                        dist, kT: float = 1.0) -> np.ndarray:
    """Score-weighted mean CA-CA distance over the pool, indexed like ``dist.i/j``.

    Measured over 16 targets, this beats the PRIOR's own predicted distances: correlation
    with the native matrix 0.698 -> 0.731 and mean absolute error 2.317 -> 2.217 A.  Uniform
    and Boltzmann weighting are within 0.01 of each other and a 50/50 blend with the prior
    is slightly worse than the pure ensemble, so the pool is not merely echoing the prior --
    it holds distance information the prior does not.
    """
    C = np.stack(pool_cas)
    s = np.asarray(pool_scores, float)
    w = np.exp(-(s - s.min()) / max(kT, 1e-9))
    w = w / w.sum()
    d = np.linalg.norm(C[:, dist.i] - C[:, dist.j], axis=-1)
    return (d * w[:, None]).sum(0)


def fold_iterated(sequence: str, ks: Sequence[int] = (8, 16), seed: int = 0,
                  targets: Optional[Dict[int, Dict]] = None, w_cons: float = 0.35,
                  passes: int = 2, kT: float = 1.0, **kw) -> Dict:
    """Search, build a consensus distance matrix from the pool, then search again.

    The second pass adds an L1 restraint toward the first pass's consensus distances.  The
    restraint weight is deliberately small: this is also the obvious way to amplify the
    prior's own mistakes, and `w_cons` is fitted on the development peptides, never on a
    benchmark target.
    """
    T = targets or {k: build_target(sequence, k=k) for k in ks}
    dist = T[ks[0]]["dist"]
    out = None
    for it in range(max(1, passes)):
        out = fold_multi(sequence, ks=ks, seed=seed + 977 * it, targets=T, **kw)
        if it == passes - 1:
            break
        cons = consensus_distances(out["pool_cas"], out["pool_scores"], dist, kT=kT)
        for k in ks:
            T[k]["objective"].set_consensus(cons, dist.w, w_cons)
    for k in ks:
        T[k]["objective"].set_consensus(None)
    out["passes"] = passes
    out["w_cons"] = w_cons
    return out
