"""SPRINT 22 / WORKSTREAM A -- the CANDIDATE-STATE quantum encoding, H|i> = E_i|i>.

Pre-registered in `s22/PREREG_A.md` before this file existed (A1/A2/A3).

WHAT THIS IS, AND WHY IT IS A DIFFERENT OBJECT FROM EVERY EARLIER CVaR-VQE ARM
------------------------------------------------------------------------------
Every previous quantum arm in this programme (s17-s20) put the circuit's register over a LATENT
that then had to be DECODED into a structure: a k-state lattice code (s17/s18) or a per-residue
von Mises basin selector with a continuously-drawn angle (s19/s20). The circuit's bitstring was
never "the candidate" -- it was an intermediate object one step removed from a structure.

Here the register IS the pool.  `n_qubits = ceil(log2(K))` computational basis states name the
`K` real retrieval-pool candidates directly (plus `2**n_qubits - K` unavoidable padding states,
which get a fixed penalty energy and are therefore never selected by any order-based readout).
`H|i> = E_i|i>` is diagonal by construction: there is no decoder, no build step, no basin.  CVaR-
VQE trained on this register is, unambiguously, a selector over a finite pool -- which is exactly
what s21 L3's exact theorem is ABOUT (a pool-restricted, same-H, order-based CVaR-VQE ties the
pool argmin).  That theorem is reproduced here as a soundness gate, not claimed as new.

What IS new is the GAUGE.  The candidate<->qubit-index bijection is arbitrary, and the ansatz's
entangling structure (a CNOT chain / MPS of bond dimension chi) imposes a locality bias over
qubit-index HAMMING distance.  A relabelling of which candidate sits at which computational basis
index scrambles that bias arbitrarily.  The EXACT optimum cannot see the label (it is a lookup);
the ACHIEVED, finite-iteration optimum might.  Sprint 18 had exactly this failure once already
(gauge percentile 0.00 on the k=4 lattice register) -- this module tests for it rather than
assuming a diagonal-basis encoding is safe by default.

EXACT, NOT SAMPLED
------------------
K <= 500 means n_qubits <= 9, i.e. the register is EXACTLY enumerable (<=512 basis states).  Every
quantity here -- the trained distribution, its CVaR value, its alpha-tail mass, its entropy, its
gradient -- is computed by full enumeration (`core.quantum.cvar_gradient_exact`,
`all_bitstrings`), not sampled.  There is no shot noise anywhere in this module.  This is a
genuine simplification the candidate-basis register buys for free: the s19/s20 continuous
encoding could never do this because its Hamiltonian was a physics call, not a lookup table.

NO AMBER.  Every energy here is either the shipped native-free distogram Bayes-risk score
(`s12.instrument.shipped_score`) or the genuine Legacy potential (`s16.energy_lib`), both cheap,
CPU-only, no OpenMM.  AMBER is a declared scope limit in `s22/PREREG_A.md` -- this module never
imports `core.amber` and never touches the serialised OpenMM context.
"""
from __future__ import annotations

import json
import math
import os
import time
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
import sys
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I               # noqa: E402
from s16 import energy_lib as EL              # noqa: E402
from s15 import seed as SD                    # noqa: E402
from core import quantum as Q                 # noqa: E402

SALT = "s22qcand"
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

K_POOL = 500


# ============================================================ pool + Hamiltonians
def build_pool(pdb: str, k: int = K_POOL) -> Dict:
    """The target's K-member retrieval pool, both native-free scores and the ORACLE label.

    `score_dist` -- the shipped Bayes-risk distogram score (native-free, lower is better).
    `score_leg`  -- genuine Legacy total (native-free, lower is better).
    `rr`         -- ORACLE Ca-RMSD to the native. LABEL ONLY: never enters a training decision.
    `W`          -- (k, n, 3) candidate CA coordinates (point-cloud basis).
    """
    t = {x["pdb"]: x for x in I.targets()}[pdb]
    u = I.load_univ(pdb)
    p = I.pool_idx(u, k=k)
    W = u["W"][p]
    PHI, PSI = u["PHI"][p], u["PSI"][p]
    rr = u["rr"][p]
    dg = I.distogram(pdb, t["seq"], t["fold"])
    i, j = I.pair_index(t["n"])
    D = I.pair_dists(W, i, j)
    score_dist = I.shipped_score(dg, D.astype(np.float32).astype(float))
    comp = EL.legacy_components_of_windows(t["seq"], PHI, PSI)
    score_leg = EL.legacy_total_from(comp)
    # native-free "typicality": mean CA-RMSD of each candidate to the rest of ITS OWN pool.
    # Used only as the A2 structural term (H_struct) and as a native-free basin-coverage
    # diagnostic -- never as a training label.
    typ = _typicality(W)
    return dict(pdb=pdb, n=int(t["n"]), seq=t["seq"], fold=int(t["fold"]), k=int(k),
                W=W, rr=rr, score_dist=score_dist, score_leg=score_leg, typicality=typ,
                nat_ca=u["nat_ca"])


def _typicality(W: np.ndarray, sample: int = 120, seed: int = 0) -> np.ndarray:
    """Mean CA-RMSD of every candidate to a fixed random subsample of its own pool.

    Native-free, O(k * sample) instead of O(k^2) -- the k=500 pairwise matrix is not needed for
    a mean-distance proxy and this keeps `build_pool` cheap enough to call per (target, seed).
    """
    k = len(W)
    rng = np.random.default_rng(seed)
    idx = rng.choice(k, size=min(sample, k), replace=False)
    d = I.kabsch_rmsd_batch(W, W[idx[0]])  # placeholder shape probe
    out = np.zeros(k)
    for ii in idx:
        out += I.kabsch_rmsd_batch(W, W[ii])
    return out / len(idx)


# ============================================================ the candidate-basis encoding
class Encoding:
    """One (target, Hamiltonian, label) instance of `H|i> = E_i|i>`.

    `label[s]` is the qubit computational-basis index holding canonical slot `s`; canonical
    slots `0..K-1` are the K real candidates in the pool's own order, `K..2**n-1` are padding.
    `cand_of_bit[b]` inverts this: the candidate index at bit-index `b`, or -1 for padding.
    """

    def __init__(self, scores: np.ndarray, label: Optional[np.ndarray] = None,
                 pad_margin: float = 10.0):
        self.K = int(len(scores))
        self.n_qubits = int(math.ceil(math.log2(max(2, self.K))))
        self.dim = 1 << self.n_qubits
        self.scores = np.asarray(scores, float)
        sd = float(self.scores.std()) if self.K > 1 else 1.0
        self.pad_energy = float(self.scores.max() + pad_margin * max(sd, 1e-9))
        self.label = (np.arange(self.dim) if label is None else np.asarray(label, int))
        assert self.label.shape == (self.dim,) and len(set(self.label.tolist())) == self.dim, \
            "label must be a permutation of range(2**n_qubits)"
        self.cand_of_bit = np.full(self.dim, -1, int)
        self.cand_of_bit[self.label[:self.K]] = np.arange(self.K)
        self.E = np.full(self.dim, self.pad_energy)
        self.E[self.label[:self.K]] = self.scores

    def bit_index_of_candidate(self, c: int) -> int:
        return int(self.label[c])

    def energies_of_bits(self, bit_idx: np.ndarray) -> np.ndarray:
        return self.E[np.asarray(bit_idx, int)]


def random_label(dim: int, pdb: str, tag: str, k: int, salt: str = SALT) -> np.ndarray:
    return SD.stable_rng(pdb, tag, "qcand_perm", k, salt=salt).permutation(dim)


def bits_matrix_to_int(bits: np.ndarray) -> np.ndarray:
    """`(B, n)` MSB-first 0/1 rows -> `(B,)` integers. Matches `core.quantum.all_bitstrings`."""
    b = np.asarray(bits, np.int64)
    n = b.shape[1]
    w = (1 << np.arange(n - 1, -1, -1)).astype(np.int64)
    return b @ w


# ============================================================ alpha schedules
def alpha_of(schedule, it: int, iters: int, entropy_bits: Optional[float] = None,
            max_entropy_bits: Optional[float] = None) -> float:
    """Resolve one of the six pre-registered alpha arms. `schedule` is a string or a float."""
    if isinstance(schedule, (int, float)):
        return float(schedule)
    prog = it / max(1, iters)
    if schedule == "annealed":
        return Q.alpha_schedule(prog, a0=0.5, a1=0.05)
    if schedule == "adaptive":
        a0, a_min = 0.5, 0.05
        if entropy_bits is None or not max_entropy_bits:
            return a0
        frac = max(0.0, min(1.0, entropy_bits / max_entropy_bits))
        return float(np.clip(frac * a0, a_min, a0))
    raise ValueError(f"unknown alpha schedule {schedule!r}")


# ============================================================ training (EXACT, no sampling)
#
# PRE-REGISTRATION ADDENDUM, dated 2026-09-07 (PREREG_A.md A2 is frozen; this is the declared,
# dated deviation, not an edit in place -- see the module docstring's pilot finding).
#
# PILOT FINDING THAT MOTIVATES IT.  On the first smoke run (`1CS9`, alpha=0.15, 200 Adam steps,
# unregularised CVaR-over-p), the trained distribution collapsed to ESS=1.0 -- a single delta on
# the global argmin, tail_size 1 against a nominal tail of 77.  This is NOT a bug: it is the exact
# fact `core.quantum.free_energy`'s own docstring already states -- "for ANY alpha the minimiser
# concentrates p on the lowest-energy basis states, so the readout collapses back to the argmin".
# For CONTINUOUS, generically-distinct candidate energies there is no numerical tie at the argmin,
# so the theoretical CVaR-minimising "face" (s21 L3c) is a single point unless the objective is
# regularised. The PREREG's beta-weighted typicality term is therefore REPLACED by the project's
# own established mechanism for keeping the ensemble interior -- the entropy-regularised free
# energy `F_alpha,T(theta) = CVaR_alpha(theta) - T * H(p_theta)` (`core.quantum.free_energy`,
# lifted here from the exact StatevectorCircuit form to the exact MPSAnsatz score-function form,
# since both are diagonal-Hamiltonian gradients and the derivation is identical). `T` is the new
# axis that was `struct_beta`; it is reported explicitly wherever it is nonzero, and `T=0`
# reproduces the pure-CVaR collapse as the baseline case, not a hidden control.
def _cvar_and_entropy_weights(p: np.ndarray, E: np.ndarray, alpha: float
                              ) -> Tuple[np.ndarray, float, np.ndarray, float]:
    """Exact score-function weights for CVaR and for the Shannon entropy of `p`, both already
    p(x)-weighted so `grad_logp_weighted(ansatz, theta, allbits, w)` returns the exact gradient
    with NO baseline subtraction (a constant baseline is a no-op only when multiplying grad_logp
    BEFORE the p-weighting, per `cvar_gradient_exact`'s own convention -- applying it to an
    already p-weighted `w` would be wrong, since `sum_x grad_logp(x)` without the p(x) factor is
    not itself zero.  This module never subtracts a baseline for that reason.)."""
    val, q, _ = Q.cvar_from_probs(E, p, alpha)
    w_cvar = -np.clip(q - E, 0.0, None) * p / alpha
    logp_safe = np.log(np.clip(p, 1e-300, None))
    H = float(-(p * logp_safe).sum())
    w_ent = -p * (1.0 + logp_safe)
    return w_cvar, val, w_ent, H


def train(enc: Encoding, alpha_schedule, seed: int, iters: int = 300, lr: float = 0.12,
         entangler: str = "cnot", T_schedule=0.0, layers: int = 2,
         checkpoints: Optional[List[int]] = None,
         theta_init: Optional[np.ndarray] = None,
         opt_state: Optional[Q.Adam] = None) -> Dict:
    """Exact-gradient CVaR-VQE (optionally entropy-regularised) on the candidate-basis register.

    No sampling anywhere: `E` is a lookup table and the full `2**n`-state distribution is
    enumerated every iteration, so every quantity is exact (see the module docstring).
    `T_schedule` is a constant, or one of the same string schedules `alpha_of` accepts, applied
    to the temperature instead of alpha (`T=0` is the pure-CVaR baseline). `checkpoints`, if
    given, stashes `(p, theta)` at those iteration counts so the collapse DYNAMICS -- not only
    the converged endpoint -- can be read off without retraining (PREREG_A A2).

    `theta_init` / `opt_state`, if given, WARM-START from a previous stage's parameters and
    Adam moments instead of a fresh random draw -- the multi-stage arm (PREREG_A A3). The random
    initial theta is still recorded (`theta0`) as whatever was used to seed the FIRST stage, so a
    staged run and its single-stage control can be compared from an identical origin.
    """
    n = enc.n_qubits
    an = Q.MPSAnsatz(n, layers=layers, final_ry=True, entangler=entangler)
    rng = np.random.default_rng(seed)
    if theta_init is None:
        th = (math.pi / 2.0) + rng.normal(0.0, 0.8, an.n_params())
    else:
        th = np.asarray(theta_init, float).copy()
    th0 = th.copy()
    opt = opt_state if opt_state is not None else Q.Adam(an.n_params(), lr=lr)
    allb = Q.all_bitstrings(n)
    E_eff = enc.E.copy()
    max_ent = float(n)  # bits
    ckset = set(checkpoints or [])
    snaps: Dict[int, Dict] = {}
    hist_val: List[float] = []
    hist_alpha: List[float] = []
    hist_T: List[float] = []
    hist_entropy: List[float] = []
    grad_norms: List[float] = []
    for it in range(1, iters + 1):
        p = np.exp(np.asarray(an.logp(th, allb), float))
        p = np.maximum(p, 0.0); p = p / p.sum()
        ent_bits = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
        a_t = alpha_of(alpha_schedule, it, iters, ent_bits, max_ent)
        T_t = alpha_of(T_schedule, it, iters, ent_bits, max_ent) if T_schedule else 0.0
        w_cvar, val, w_ent, H = _cvar_and_entropy_weights(p, E_eff, a_t)
        w = w_cvar - T_t * w_ent if T_t else w_cvar
        g = Q.grad_logp_weighted(an, th, allb, w)
        grad_norms.append(float(np.linalg.norm(g)))
        th = opt.step(th, g)
        hist_val.append(val)
        hist_alpha.append(a_t)
        hist_T.append(T_t)
        hist_entropy.append(H)
        if it in ckset:
            snaps[it] = dict(theta=th.copy(), p=p.copy(), alpha=a_t, T=T_t, cvar=val,
                             entropy_bits=ent_bits)
    p_final = np.exp(np.asarray(an.logp(th, allb), float))
    p_final = np.maximum(p_final, 0.0); p_final = p_final / p_final.sum()
    return dict(theta=th, theta0=th0, ansatz=an, allbits=allb, p_final=p_final,
                E_eff=E_eff, hist_val=hist_val, hist_alpha=hist_alpha, hist_T=hist_T,
                hist_entropy=hist_entropy, snapshots=snaps, opt=opt,
                grad_norm_mean=float(np.mean(grad_norms)), grad_norm_sd=float(np.std(grad_norms)),
                iters=iters, n_qubits=n, entangler=entangler, seed=seed)


def untrained_distribution(n_qubits: int, entangler: str, seed: int,
                           layers: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """The ansatz's distribution at its RANDOM INITIAL theta -- never an initialisation mean.

    Returns `(p, allbits)` at the init the Adam run above would have started from with the SAME
    seed, so `best_of_N` from the untrained circuit is a like-for-like control.
    """
    an = Q.MPSAnsatz(n_qubits, layers=layers, final_ry=True, entangler=entangler)
    rng = np.random.default_rng(seed)
    th0 = (math.pi / 2.0) + rng.normal(0.0, 0.8, an.n_params())
    allb = Q.all_bitstrings(n_qubits)
    p = np.exp(np.asarray(an.logp(th0, allb), float))
    p = np.maximum(p, 0.0); p = p / p.sum()
    return p, allb


# ============================================================ exact readouts (all native-FREE
# except the RMSD scoring, which is post-hoc ORACLE only)
def exact_face(E: np.ndarray, p: np.ndarray, alpha: float) -> Dict:
    """The CVaR-optimal FACE, exactly: which bit-indices carry positive tail mass, and how much.

    `cvar_from_probs` orders states by ENERGY (not by probability) and accumulates probability
    mass bottom-up until `alpha` is reached -- so `mass > 0` is exactly the alpha-tail support of
    the distribution `p` against the fixed energies `E`. This is the FACE (s21 L3c / BRIEF item 2)
    read off exactly, not sampled.
    """
    val, q, mass = Q.cvar_from_probs(E, p, alpha)
    tail = mass > 0
    ess = float(alpha ** 2 / max(1e-300, (mass ** 2).sum()))
    return dict(cvar=val, quantile=q, mass=mass, tail=tail, tail_size=int(tail.sum()),
               ess=ess, nominal_tail_size=int(math.ceil(alpha * len(E))))


#: TIE-BREAKING HAZARD, caught in a pilot run before any A1 number was read for the record
#: (target `5H1H`: top-2 distogram scores are EXACTLY tied at 2.1144738, gap 0.0). A naive
#: `np.argmin` over bit-index order picks whichever tied candidate happens to land FIRST under
#: the label currently in force -- so a random relabelling can flip which member of an exact tie
#: is "the" argmin, manufacturing an apparent gauge failure that is actually project memory's
#: documented trap (`consensus-is-the-only-in-band-discriminator`'s "tie-breaking leaks the pool
#: order"). Both readouts below return the FULL tied set and average the ORACLE RMSD over it,
#: which is gauge-invariant by construction (the tied set is a set of CANDIDATE indices, not of
#: bit indices, so no label can reorder it) and is the convention project memory records.
_TIE_RTOL = 1e-9


def tied_indices(values: np.ndarray, at: float) -> np.ndarray:
    v = np.asarray(values, float)
    tol = max(_TIE_RTOL, _TIE_RTOL * abs(at))
    return np.flatnonzero(np.abs(v - at) <= tol)


def argmin_readout(enc: Encoding, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Order-based readout: the FULL SET of candidates tied at the argmin energy among states
    with `p>0`. Returns `(bit_indices, candidate_indices)`, both possibly length > 1 under a
    numerical tie.

    Reproduces s21 L3's theorem: for a full-support ansatz (generic RY angles put non-zero mass
    everywhere) this set equals the GLOBAL argmin SET regardless of training quality or label --
    which is exactly the soundness gate A1 registers, not a novel claim.
    """
    support = p > 1e-15
    if not support.any():
        support = np.ones_like(p, bool)
    idx = np.arange(len(p))[support]
    emin = float(enc.E[idx].min())
    tied = idx[tied_indices(enc.E[idx], emin)]
    cands = enc.cand_of_bit[tied]
    cands = cands[cands >= 0]
    return tied, cands


def tail_candidates(enc: Encoding, tail_mask: np.ndarray) -> np.ndarray:
    """Real (non-padding) candidate indices inside the tail-support set."""
    bits = np.flatnonzero(tail_mask)
    cands = enc.cand_of_bit[bits]
    return cands[cands >= 0]


def tail_average_rmsd(pool: Dict, cand_idx: np.ndarray) -> Tuple[float, np.ndarray]:
    """Coordinate-average of the given candidates, ORACLE-scored against the native (post-hoc)."""
    if len(cand_idx) == 0:
        return float("nan"), None
    C, _b = I.coordinate_average(pool["W"][cand_idx])
    return float(I.ca_rmsd(C, pool["nat_ca"])), C


def repaired_rmsd(pool: Dict, C: np.ndarray) -> float:
    """`R_repair`: the ideal-geometry projection of a coordinate average (s12.instrument.project).
    ORACLE-scored, post-hoc."""
    pr = I.project(np.asarray(C, float), pool["seq"], pool["fold"])
    return float(I.ca_rmsd(np.asarray(pr["ca"], float), pool["nat_ca"]))


# ============================================================ classical controls
def classical_exact_sort(pool: Dict, scores: np.ndarray, alpha: float
                         ) -> Dict:
    """The GREEDY / EXHAUSTIVE classical control: exact score sort, deterministic.

    In a <=512-state register, exhaustive scoring costs K evaluations -- cheaper than any
    training budget used anywhere in this programme -- so "greedy" and "exact argmin" coincide
    here and that identity is stated rather than hidden (PREREG_A A1).
    """
    order = np.argsort(scores, kind="stable")
    m = max(1, int(math.ceil(alpha * len(scores))))
    top = order[:m]
    emin = float(scores.min())
    argmin_set = tied_indices(scores, emin)          # gauge-invariant: a set of CANDIDATES
    r_argmin, _ = tail_average_rmsd(pool, argmin_set)
    r_tailavg, C = tail_average_rmsd(pool, top)
    return dict(argmin_cand=int(order[0]), argmin_set=argmin_set, r_argmin=r_argmin,
               tail_cands=top, r_tailavg=r_tailavg, coord_avg=C, m=int(m))


def random_subset_control(pool: Dict, m: int, n_draws: int = 16, seed: int = 0
                          ) -> Dict:
    """Matched-count, zero-information null: `n_draws` random subsets of size `m`, averaged."""
    rng = np.random.default_rng(seed)
    k = len(pool["score_dist"])
    r_single, r_avg = [], []
    for d in range(n_draws):
        idx = rng.choice(k, size=m, replace=False)
        r_avg.append(tail_average_rmsd(pool, idx)[0])
        r_single.append(pool["rr"][int(idx[0])])
    return dict(r_argmin_mean=float(np.mean(r_single)), r_tailavg_mean=float(np.mean(r_avg)),
               r_argmin_all=r_single, r_tailavg_all=r_avg)


def simulated_annealing_control(scores: np.ndarray, seed: int, n_evals: Optional[int] = None
                                ) -> Dict:
    """Matched-budget classical SA over candidate identity (discrete Metropolis on the index).

    Near-vacuous by construction here (PREREG_A A1): the register is exhaustively scorable in
    K evaluations, so SA is run at that SAME budget for completeness under the Hard Requirements,
    not because search is expected to be the binding constraint.
    """
    rng = np.random.default_rng(seed)
    k = len(scores)
    n_evals = n_evals or k
    cur = int(rng.integers(k))
    best = cur
    for t in range(n_evals):
        T = max(1e-6, 1.0 * (1.0 - t / n_evals))
        cand = int(rng.integers(k))
        de = scores[cand] - scores[cur]
        if de <= 0 or rng.random() < math.exp(-de / T):
            cur = cand
        if scores[cur] < scores[best]:
            best = cur
    return dict(best_cand=best, found_global=bool(best == int(np.argmin(scores))))


# ============================================================ io
def write(name: str, obj: Dict, required_keys: Optional[List[str]] = None) -> str:
    """Atomic write + a `_COMPLETE` sidecar demanding the FULL key set (BRIEF S3 rules 5-6)."""
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, path)
            break
        except PermissionError:
            time.sleep(0.4)
    rows = obj.get("rows", [])
    ok = bool(required_keys) and len(rows) > 0 and all(
        all(k in r for k in required_keys) for r in rows)
    if required_keys is not None:
        with open(path.replace(".json", "") + "_COMPLETE", "w") as fh:
            fh.write(json.dumps({"complete": ok, "n_rows": len(rows),
                                 "required_keys": required_keys,
                                 "time": time.strftime("%Y-%m-%d %H:%M:%S")}))
    return path


def paired_ci(a, b, n_boot=4000, seed=0):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    d = a - b
    n = len(d)
    if n < 3:
        return {"n": int(n), "mean": float(d.mean()) if n else float("nan")}
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    se = float(d.std(ddof=1) / math.sqrt(n))
    return {"n": int(n), "mean": float(d.mean()), "median": float(np.median(d)), "se": se,
           "mde": float(2.8016 * se), "ci95": [float(np.percentile(bs, 2.5)),
                                                float(np.percentile(bs, 97.5))],
           "W": int((d < 0).sum()), "L": int((d > 0).sum())}
