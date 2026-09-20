#!/usr/bin/env python
"""s29/s29_X_config.py -- S29 LANE X: THE CONFIGURATION-SPACE CVaR-VQE (fragment recombination).

Pre-registered in `s29/PREREG_S29_X.md` (read it first; the nine questions, the falsifiers and
the registered prior live there, not here).

WHAT IS NEW RELATIVE TO THE RECORD
==================================
Every state space this project has run a CVaR-VQE over is either (a) the CANDIDATE-IDENTITY
register (a basis state is a pool member; H is diagonal; the tail is a classical prefix to
1e-13, S28-L21, so the quantum stage can be deleted), (b) the per-residue TORSION-BIN lattice
(S13; k = 4 Ramachandran bins; exhaustively enumerated, no native-free objective finds the good
configurations), or (c) the per-residue BASIN latent (S19-S21; one qubit per residue, von Mises
draws from a parametric prior FITTED to the pool's own marginals; exhaustively enumerated at
n = 126, the exact argmin ties a zero-evaluation pool).

Here a basis state is a CHIMERA: the chain is cut into S contiguous segments and each segment
independently takes its backbone torsions from one of F = 8 retrieved pool members (the DIS
top-8 of the shipped K = 500 pool).  The space is the combinatorial RECOMBINATION of real
retrieved fragments: it CONTAINS its eight parents, it is not in the span of the pool average,
and every configuration is built with the production ideal-geometry builder.  log2(F) = 3 qubits
per segment, q = 3 S qubits, and 2**q = F**S exactly, so the register is padding-free and the
whole space is exactly enumerable -- which makes every classical control EXACT rather than
sampled.

THE HAMILTONIAN
===============
    H_diag(x) = E_pair(x) + E_rama(x)          nats
    E_pair    = -sum_{j-i>=2} log( p_ij[bin(d_ij(x))] + eps )   the shipped leave-fold-out
                distogram's 17-bin posterior, consumed JOINTLY as a log score on the built
                configuration (not through its median, its mean or the L1 Bayes risk)
    E_rama    = -sum_i log P_fold(phi_i, psi_i | aa_i)          the per-fold leakage-safe
                Ramachandran table, the S27 RAMA channel's binning
    H_mix     = -Gamma sum_k X_k               the transverse field: transitions between
                chimeras differing in one bit of one segment's fragment index

    F(theta)  = CVaR_alpha(H_diag; p_theta) - T S(p_theta) - Gamma <psi|sum_k X_k|psi>

At Gamma = 0 the objective and its gradient are `core.quantum.free_energy` LINE FOR LINE and the
training loop is `core.quantum.run_cvar_vqe` (same RNG, same Adam, same settings), so Gamma = 0
reproduces the deployed driver bit for bit (`tests/test_s29_X.py`).  For Gamma > 0 the objective
is NOT a function of p alone -- it depends on the relative SIGNS of the amplitudes on
Hamming-neighbour chimeras -- and the set-equality theorem's premise (H diagonal) fails.

THE READOUT
===========
The CVaR tail's mass vector comes from `core.quantum.cvar_from_probs` (the deployed extraction).
    R1  uniform coordinate average of the tail SET          the DEPLOYED operator
    R2  tail-MASS-weighted average of the same set          the CVaR tail ENSEMBLE (primary)
    R3  probability-weighted average of the state's top-512 configurations   diagnostic
then `s12.instrument.project` -- the production projection -- for the built chain, which is the
reporting basis.  R2 - R1 on the identical set is the whole of what the quantum state can
contribute beyond the tail size m.

NATIVE-FREE / ORACLE SPLIT
==========================
Nothing outside `oracle_*` reads a native.  `--poison` replaces `nat_ca` with NaN and asserts
every emitted structure is bit-identical (the contract's NaN-poison test, also in the test file).

USAGE
=====
    python s29/s29_X_config.py --run [--limit N] [--targets 1A13,...]   resumable, per target
    python s29/s29_X_config.py --analyse                                ST.fmt blocks
    python s29/s29_X_config.py --selftest                               fast invariants
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                          # noqa: E402  (after the thread caps)

from core import predict as PRD             # noqa: E402
from core import project as PJ              # noqa: E402
from core import quantum as Q               # noqa: E402
from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402
from s24 import d_harness as H              # noqa: E402
from s24 import stats_lib as ST             # noqa: E402
from s25 import phys_lib as P               # noqa: E402
from s27 import ham_lib as HL               # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

SALT = "s29X"
# ---- every constant below is pre-registered in s29/PREREG_S29_X.md section 8; none is
# ---- chosen on an RMSD, and none may be changed after a result exists.
F_MEMBERS = 8               # fragment alternatives per segment -> 3 qubits
SEG_LEN = 3                 # Rosetta 3-mers
SEG_MAX = 5                 # register cap: q = 3*S <= 15 (addendum 1: measured cost, not taste)
EPS_P = 1e-4                # posterior floor: one pair costs at most 9.21 nats
ALPHA = 0.15                # production's tail fraction
TEMP = 1.0                  # nats: the posterior's own temperature. NEVER tuned.
LAYERS = 3
ITERS = 80
LR = 0.15
SEEDS = (0, 1)
R3_CAP = 512
MEDOID_EXACT_MAX = 600      # above this the medoid is estimated on a seeded 256-subsample
MEDOID_SUBSAMPLE = 256
SA_COOL = (20.0, 1.0)       # nats, geometric
#: ADDENDUM 2 (S29-L15, lane T's falsifier ladder): the mixer coupling grid, in units of the
#: native-free `gamma_gap`.  Gate 1 is TV(p_Gamma, p_0) > TV_GATE on the SAMPLED distribution;
#: below it the readout provably cannot resolve the difference (S25 L15) and the cell is empty.
GAMMA_GRID = (0.0, 0.5, 1.0, 2.0)
TV_GATE = 0.45

#: S27's trainability set (`s27/run_trainability.py`: `P.targets()[::11][:12]`).
PROBE_TARGETS = ["1A13", "1I6Y", "1M02", "2BFI", "2LWS", "2MP9",
                 "2P5H", "5Z5W", "6MBM", "7JGX", "8HVS", "9KAR"]


# ============================================================ the configuration space
def segments_of(n: int) -> List[np.ndarray]:
    """Contiguous near-equal segments: S = min(ceil(n/SEG_LEN), SEG_MAX) blocks."""
    s = min(int(math.ceil(n / SEG_LEN)), SEG_MAX)
    return [np.asarray(a, int) for a in np.array_split(np.arange(n), s)]


def seg_of_residue(n: int) -> np.ndarray:
    segs = segments_of(n)
    out = np.empty(n, int)
    for s, ix in enumerate(segs):
        out[ix] = s
    return out


def config_digits(S: int, F: int = F_MEMBERS) -> np.ndarray:
    """(F**S, S) base-F digits of every basis index, segment 0 = MOST significant.

    The bit order matches `core.quantum.StatevectorCircuit`, whose qubit q of basis index
    `idx` is `(idx >> (n-1-q)) & 1`: qubit 0 is the most significant bit.  Segment s owns
    qubits [3s, 3s+3), so segment 0 owns the top three bits.
    """
    M = F ** S
    x = np.arange(M, dtype=np.int64)
    return np.stack([(x // (F ** (S - 1 - s))) % F for s in range(S)], axis=1)


def member_indices(pdb: str, cand, dis: np.ndarray, F: int = F_MEMBERS) -> np.ndarray:
    """The F parents: the DIS top-F of the shipped pool, ties broken by a stable random key.

    `s27.run_pool.topm`'s rule (contract item 12: ties never break by array order).
    """
    key = SD.stable_rng(pdb, "members", salt=SALT).random(len(dis))
    order = np.lexsort((key, np.asarray(dis, float)))
    return np.asarray(order[:F], int)


class Space:
    """The chimera space of one target, fully enumerated.  Native-free by construction.

    `scramble=True` builds the MATCHED NULL lane D required (S29-L4 hole (a)): the identical
    cardinality, the identical 8 parents and the identical marginal fragment content, with each
    member's segment blocks permuted ACROSS segment positions so that no configuration is a
    compatible recombination.  Its ORACLE best is the order-statistic control for D1: a minimum
    over 8**S structures is smaller than a minimum over 8 whatever the space contains
    (`grid-oracles-are-order-statistics`), and only the gap to THIS null is recombination.
    """

    def __init__(self, pdb: str, cand, dg, rama_cnt, F: int = F_MEMBERS, scramble: bool = False):
        self.pdb, self.cand, self.dg = pdb, cand, dg
        self.n, self.seq, self.fold = cand.n, cand.seq, cand.fold
        self.F = int(F)
        self.dis = H.score_shipped(cand)
        self.members = member_indices(pdb, cand, self.dis, F)
        self.segs = segments_of(self.n)
        self.S = len(self.segs)
        self.q = self.S * int(round(math.log2(self.F)))
        self.M = self.F ** self.S
        assert self.M == 1 << self.q, "register must be padding-free"
        self.digits = config_digits(self.S, self.F)                    # (M, S)
        self.sor = seg_of_residue(self.n)
        self.scramble = bool(scramble)
        PHI_m = np.asarray(cand.PHI, float)[self.members]              # (F, n)
        PSI_m = np.asarray(cand.PSI, float)[self.members]
        if scramble:
            PHI_b, PSI_b = self._scrambled_blocks(PHI_m, PSI_m)
        else:
            PHI_b = np.stack([PHI_m] * 1)[0]                           # (F, n), unchanged
            PSI_b = PSI_m
        idx = self.digits[:, self.sor]                                 # (M, n) member per res
        rr = np.arange(self.n)[None, :]
        self.PHI = PHI_b[idx, rr]
        self.PSI = PSI_b[idx, rr]
        PHI_m, PSI_m = PHI_b, PSI_b
        self.CA = np.asarray(PJ.build_ca_exact(self.PHI, self.PSI), float)   # (M, n, 3)
        self.i, self.j = I.pair_index(self.n, 2)
        assert np.array_equal(self.i, np.asarray(dg["i"])) and \
               np.array_equal(self.j, np.asarray(dg["j"])), "pair index must match the distogram"
        self.rama_cnt = rama_cnt
        self.E_pair = self._pair_energy(self.CA)
        self.E_rama_tab = self._rama_table(PHI_m, PSI_m)               # (S, F)
        self.E_rama = self.E_rama_tab[np.arange(self.S)[None, :], self.digits].sum(1)
        self.E = self.E_pair + self.E_rama
        self.gamma_gap = self._gamma_gap()

    def _scrambled_blocks(self, PHI_m: np.ndarray, PSI_m: np.ndarray):
        """Each member's segment blocks permuted across segment POSITIONS (the matched null).

        Cardinality, parents and marginal fragment content are identical; what is destroyed is
        that a member's block sits where that member put it.  A source block shorter than its
        destination is cycled, so every destination residue receives a real torsion pair.
        """
        rng = SD.stable_rng(self.pdb, "scramble", salt=SALT)
        PH, PS = PHI_m.copy(), PSI_m.copy()
        for f in range(self.F):
            perm = rng.permutation(self.S)
            for s, ix in enumerate(self.segs):
                src = self.segs[perm[s]]
                take = np.array([src[t % len(src)] for t in range(len(ix))], int)
                PH[f, ix] = PHI_m[f, take]
                PS[f, ix] = PSI_m[f, take]
        return PH, PS

    # -- the two energy terms ---------------------------------------------
    def _pair_energy(self, CA: np.ndarray, chunk: int = 8192) -> np.ndarray:
        prob = np.asarray(self.dg["prob"], float)                      # (npairs, 17)
        lp = -np.log(prob + EPS_P)
        out = np.empty(len(CA))
        for a in range(0, len(CA), chunk):
            D = I.pair_dists(CA[a:a + chunk], self.i, self.j)
            b = np.digitize(D, PRD.BIN_EDGES)
            out[a:a + chunk] = lp[np.arange(lp.shape[0])[None, :], b].sum(1)
        return out

    def _rama_table(self, PHI_m: np.ndarray, PSI_m: np.ndarray) -> np.ndarray:
        """(S, F) 1-body Ramachandran cost.  Residue n-1 is inert for the CA trace (S13
        correction 10) and is excluded; every other residue contributes its (phi, psi) cell."""
        cnt = np.asarray(self.rama_cnt, float) + HL.PSEUDO             # (20, 36, 36)
        Pm = cnt / cnt.sum(axis=(1, 2), keepdims=True)
        logp = np.log(Pm)
        aa = HL.codes(self.seq)
        b = 36
        pb = np.clip(((PHI_m + math.pi) / (2 * math.pi) * b).astype(int), 0, b - 1)
        sb = np.clip(((PSI_m + math.pi) / (2 * math.pi) * b).astype(int), 0, b - 1)
        per = -logp[aa[None, :].repeat(self.F, 0), pb, sb]             # (F, n)
        live = np.ones(self.n, bool)
        live[self.n - 1] = False
        tab = np.zeros((self.S, self.F))
        for s, ix in enumerate(self.segs):
            k = ix[live[ix]]
            if len(k):
                tab[s] = per[:, k].sum(1)
        return tab

    def _gamma_gap(self) -> float:
        """Median |H(x) - H(x xor e_k)| over the H top-75 configurations and the q bits.

        Native-free: it is a property of H alone.  It is the scale at which a first-order
        neighbour amplitude Gamma / dE is order one, i.e. the only non-arbitrary Gamma.
        """
        top = np.argsort(self.E, kind="stable")[:min(75, self.M)]
        g = []
        for k in range(self.q):
            nb = top ^ (1 << (self.q - 1 - k))
            g.append(np.abs(self.E[top] - self.E[nb]))
        return float(np.median(np.concatenate(g)))

    # -- diagnostics -------------------------------------------------------
    def parents_rows(self) -> Dict:
        """Which configuration is parent f?  All S digits equal f."""
        base = sum(f * self.F ** s for s, f in [(s, 1) for s in range(self.S)])
        return {int(f): int(f * (self.M - 1) // (self.F - 1)) for f in range(self.F)}

    def oracle_rmsd_all(self) -> np.ndarray:
        """ORACLE: CA-RMSD of every chimera to the native.  Diagnostic only."""
        return I.kabsch_rmsd_batch(self.CA, self.cand.nat_ca)


# =================================================================== the quantum stage
class ProductCircuit(Q.StatevectorCircuit):
    """The PRODUCT-STATE RESTRICTION: the identical circuit with the entangler deleted."""

    def _entangler_permutation(self) -> np.ndarray:
        return np.arange(self.dim, dtype=np.int64)


def xor_index(dim: int, q: int, k: int) -> np.ndarray:
    return np.arange(dim, dtype=np.int64) ^ (1 << (q - 1 - k))


def mixer_value(psi: np.ndarray, xors: Sequence[np.ndarray]) -> float:
    """<psi| sum_k X_k |psi> for a real, normalised statevector."""
    return float(sum(float(psi @ psi[x]) for x in xors))


def f_and_grad(circ, theta: np.ndarray, E: np.ndarray, alpha: float, T: float,
               gamma: float, xors: Sequence[np.ndarray]):
    """F = CVaR_alpha(E; p) - T S(p) - gamma <psi|sum X_k|psi>, and its EXACT gradient.

    Every term is an expectation of a projector or of a Pauli string, so the two-term
    parameter-shift rule is exact for all of them; the shifted circuits are simulated ONCE
    and both the probability and the amplitude quantities are read off the same pass.

    At gamma = 0 this is `core.quantum.free_energy` line for line (asserted in the tests).
    """
    psi = circ.state(theta)
    p = psi ** 2
    p = p / p.sum()
    v, qv, dp = Q.cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-15))
    S = float(-(p * lp).sum())
    dS = -(lp + 1.0)
    Mv = mixer_value(psi / max(np.linalg.norm(psi), 1e-300), xors) if gamma else 0.0
    A = circ.states_batch(circ._shift_grid(np.asarray(theta, float), np.pi / 2))
    PR = A ** 2
    PR = PR / PR.sum(1, keepdims=True)
    d = dp - T * dS
    g = (PR[0::2] - PR[1::2]) @ d / 2.0
    if gamma:
        nrm = np.maximum(np.linalg.norm(A, axis=1), 1e-300)[:, None]
        An = A / nrm
        Mb = np.zeros(len(A))
        for x in xors:
            Mb += (An * An[:, x]).sum(1)
        g = g - gamma * (Mb[0::2] - Mb[1::2]) / 2.0
    return float(v - T * S - gamma * Mv), g, p, psi, float(v), S, float(Mv)


def train(circ, E: np.ndarray, alpha: float, T: float, gamma: float,
          xors: Sequence[np.ndarray], seed: int, iters: int = ITERS, lr: float = LR):
    """Adam on the exact gradient.  Identical to `core.quantum.run_cvar_vqe` at gamma = 0."""
    rng = np.random.default_rng(int(seed))
    th = rng.normal(0.0, 0.6, circ.n_params())
    th0 = th.copy()
    m = np.zeros_like(th)
    v = np.zeros_like(th)
    hist = []
    for t in range(1, int(iters) + 1):
        f, g, _, _, _, _, _ = f_and_grad(circ, th, E, alpha, T, gamma, xors)
        hist.append(f)
        m = 0.9 * m + 0.1 * g
        v = 0.999 * v + 0.001 * g * g
        th = th - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
    f, _, p, psi, cv, S, Mv = f_and_grad(circ, th, E, alpha, T, gamma, xors)
    return dict(theta=th, theta0=th0, F=f, p=p, psi=psi, cvar=cv, entropy=S, mixer=Mv,
                hist=hist, seed=int(seed), gamma=float(gamma))


def ground_state(E: np.ndarray, gamma: float, q: int, xors: Sequence[np.ndarray]):
    """The EXACT ground state of diag(E) - gamma sum_k X_k (the diagonalised equivalent)."""
    from scipy.sparse.linalg import LinearOperator, eigsh
    dim = len(E)

    def mv(v):
        out = E * v
        if gamma:
            acc = np.zeros_like(v)
            for x in xors:
                acc += v[x]
            out = out - gamma * acc
        return out

    op = LinearOperator((dim, dim), matvec=mv, dtype=float)
    w, V = eigsh(op, k=1, which="SA", tol=1e-9, maxiter=20000)
    g = np.asarray(V[:, 0], float)
    g = g / np.linalg.norm(g)
    return float(w[0]), g


def cvar_optimal_law(E: np.ndarray, alpha: float, T: float) -> Tuple[np.ndarray, float]:
    """The EXACT minimiser of CVaR_alpha(E; p) - T S(p) over the simplex (lane T, S29-L15 Q1).

    By Rockafellar-Uryasev plus Sion's minimax theorem the optimum is
    `p*(x) ~ exp((t* - E_x)_+ / (alpha T))` with `t*` the maximiser of
    `t - T log sum_x exp((t - E_x)_+ / (alpha T))`.  This is lane T's M6 control in the
    strongest available form: the exact optimum of the SAME objective with no circuit, no
    optimiser and no per-target training -- what the quantum stage must beat to contribute.
    """
    E = np.asarray(E, float)
    s = max(alpha * T, 1e-12)

    def obj(t):
        z = np.maximum(t - E, 0.0) / s
        return t - T * (float(np.max(z)) + math.log(float(np.exp(z - z.max()).sum())))
    lo, hi = float(E.min()) - 1.0, float(E.max()) + 1.0
    for _ in range(200):                                # golden-section on a concave function
        a = lo + 0.381966 * (hi - lo)
        b = hi - 0.381966 * (hi - lo)
        if obj(a) < obj(b):
            lo = a
        else:
            hi = b
    t = 0.5 * (lo + hi)
    z = np.maximum(t - E, 0.0) / s
    p = np.exp(z - z.max())
    return p / p.sum(), float(t)


def total_variation(p: np.ndarray, q: np.ndarray) -> float:
    return float(0.5 * np.abs(np.asarray(p, float) - np.asarray(q, float)).sum())


def gibbs(E: np.ndarray, T: float) -> np.ndarray:
    z = -(np.asarray(E, float) - E.min()) / max(float(T), 1e-12)
    p = np.exp(z - z.max())
    return p / p.sum()


def gibbs_matched_entropy(E: np.ndarray, target_S: float) -> Tuple[np.ndarray, float]:
    """The thermal ensemble whose Shannon entropy equals the trained state's (bisection)."""
    def ent(T):
        p = gibbs(E, T)
        pp = p[p > 0]
        return float(-(pp * np.log(pp)).sum())
    lo, hi = 1e-3, 1e3
    if ent(hi) < target_S:
        return gibbs(E, hi), hi
    if ent(lo) > target_S:
        return gibbs(E, lo), lo
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if ent(mid) < target_S:
            lo = mid
        else:
            hi = mid
    T = math.sqrt(lo * hi)
    return gibbs(E, T), T


def anneal(E: np.ndarray, q: int, budget: int, seed: int,
           cool: Tuple[float, float] = SA_COOL) -> Dict:
    """Single-flip Metropolis over the q configuration bits at a matched evaluation budget.

    Returns the visit counts, so the SA tail is read by the IDENTICAL CVaR extraction the
    quantum arm uses (the control matched in the operator's space, contract rule 7).
    """
    rng = SD.stable_rng("sa", int(seed), int(budget), salt=SALT)
    x = int(rng.integers(0, len(E)))
    ex = float(E[x])
    counts = {}
    T0, T1 = cool
    for t in range(int(budget)):
        T = T0 * (T1 / T0) ** (t / max(int(budget) - 1, 1))
        k = int(rng.integers(0, q))
        y = x ^ (1 << (q - 1 - k))
        ey = float(E[y])
        if ey <= ex or rng.random() < math.exp(-(ey - ex) / max(T, 1e-12)):
            x, ex = y, ey
        counts[x] = counts.get(x, 0) + 1
    p = np.zeros(len(E))
    for k, v in counts.items():
        p[k] = v
    p = p / p.sum()
    return dict(p=p, n_distinct=int((p > 0).sum()), budget=int(budget), seed=int(seed))


# ======================================================================== the readouts
def _medoid_frame(Wset: np.ndarray, pdb: str, tag: str) -> int:
    """The medoid of the set.  Exact for small sets; a seeded 256-reference estimate above
    `MEDOID_EXACT_MAX`, because the medoid criterion is a MEAN over references and a
    subsample of references is an unbiased estimate of it.  Native-free either way."""
    m = len(Wset)
    if m <= MEDOID_EXACT_MAX:
        return I.medoid(I.pairwise_rmsd(Wset))
    rng = SD.stable_rng(pdb, "medoid", tag, m, salt=SALT)
    ref = rng.choice(m, size=MEDOID_SUBSAMPLE, replace=False)
    acc = np.zeros(m)
    for r in ref:
        acc += I.kabsch_rmsd_batch(Wset, Wset[r])
    return int(np.argmin(acc))


def weighted_average(Wset: np.ndarray, w: Optional[np.ndarray], pdb: str, tag: str) -> np.ndarray:
    """S8-11's operator with weights: superpose on the medoid, weighted mean.

    With uniform weights this is `s12.instrument.coordinate_average` exactly (tested)."""
    Wset = np.asarray(Wset, float)
    b = _medoid_frame(Wset, pdb, tag)
    Ws = I.superpose_batch(Wset, Wset[b])
    if w is None:
        return Ws.mean(0)
    w = np.asarray(w, float)
    return (w[:, None, None] * Ws).sum(0) / w.sum()


def tail_readouts(space: Space, p: np.ndarray, tag: str) -> Dict:
    """R1 / R2 / R3 from a distribution over configurations, plus the tail's description."""
    val, qv, mass = Q.cvar_from_probs(space.E, p, ALPHA)
    idx = np.nonzero(mass > 0)[0]
    m = int(len(idx))
    W = space.CA[idx]
    R1 = weighted_average(W, None, space.pdb, tag + "|R1")
    R2 = weighted_average(W, mass[idx], space.pdb, tag + "|R2")
    top = np.argsort(-p, kind="stable")[:R3_CAP]
    R3 = weighted_average(space.CA[top], p[top], space.pdb, tag + "|R3")
    pp = p[p > 0]
    w = mass[idx] / mass[idx].sum()
    return dict(R1=R1, R2=R2, R3=R3, m=m, tail_idx=idx, cvar=float(val),
                entropy=float(-(pp * np.log(pp)).sum()),
                ess=float(1.0 / np.square(p).sum()),
                # S29-L4 hole (c): R2 - R1 over the same set is mechanically a difference in
                # EFFECTIVE set size, so the weights' participation ratio is printed beside it.
                pr=float(1.0 / np.square(w).sum()),
                # S29-L4 hole (f): R3 is a top-512 readout unless its captured mass is stated.
                r3_mass=float(p[top].sum()), r3_k=int(len(top)),
                tail_mass_top=float(mass[idx].max() / ALPHA) if m else float("nan"))


def random_weight_control(space: Space, idx: np.ndarray, pr: float, tag: str,
                          n_draws: int = 8) -> List[np.ndarray]:
    """S29-L4 hole (c): Dirichlet weights on the SAME set at the SAME participation ratio.

    If R2 - R1 is reproduced by random weights of matched concentration, the quantum stage
    contributed CONCENTRATION, not information.  `Dirichlet(a, ..., a)` on m components has
    E[PR] ~ m(ma + 1)/(m a + m) ... solved numerically here by bisection on a, which is exact
    enough because the realised PR of each draw is recorded and averaged.
    """
    m = len(idx)
    if m < 2:
        return []

    def pr_of(a, rng):
        w = rng.dirichlet(np.full(m, a))
        return 1.0 / np.square(w).sum()
    lo, hi = 1e-3, 1e4
    for _ in range(40):
        mid = math.sqrt(lo * hi)
        rg = SD.stable_rng(space.pdb, "prsolve", tag, salt=SALT)
        v = float(np.mean([pr_of(mid, rg) for _ in range(6)]))
        if v < pr:
            lo = mid
        else:
            hi = mid
    a = math.sqrt(lo * hi)
    rng = SD.stable_rng(space.pdb, "prctrl", tag, salt=SALT)
    out = []
    for _ in range(int(n_draws)):
        w = rng.dirichlet(np.full(m, a))
        out.append(weighted_average(space.CA[idx], w, space.pdb, tag + "|RW"))
    return out


def shape_of(C: np.ndarray) -> Dict[str, float]:
    """Contract addendum 20(c): the emitted structure's Rg and mean virtual bond."""
    C = np.asarray(C, float)
    g = C - C.mean(0, keepdims=True)
    return dict(rg=float(np.sqrt((g ** 2).sum(1).mean())),
                bond=float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean()))


def uniform_topm(space: Space, m: int, tag: str) -> np.ndarray:
    key = SD.stable_rng(space.pdb, "topm", int(m), salt=SALT).random(space.M)
    order = np.lexsort((key, space.E))
    return weighted_average(space.CA[order[:max(1, int(m))]], None, space.pdb, tag)


# ========================================================================= the driver
def target_path(pdb: str) -> str:
    return os.path.join(RESULTS, f"s29_X_probe_{pdb}.json")


def load_target(pdb: str):
    cand = H.Candidates.from_universe(pdb, k=500)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    rama = np.load(os.path.join(ROOT, "s8", "generate_rama.npz"))["cnt"][cand.fold]
    return cand, dg, rama


def emit(space: Space, name: str, C: np.ndarray, out: Dict, chain: bool = True,
         extra: Optional[Dict] = None) -> None:
    """Score one emitted structure on BOTH bases and record it.  The ONLY native read."""
    row = dict(arm=name, rmsd_cloud=float(I.ca_rmsd(C, space.cand.nat_ca)), **shape_of(C))
    if chain:
        pr = I.project(C, space.seq, space.fold)
        row["rmsd_chain"] = float(I.ca_rmsd(pr["ca"], space.cand.nat_ca))
    if extra:
        row.update(extra)
    out["arms"].append(row)


def run_target(pdb: str, chain: bool = True, verbose: bool = True) -> Dict:
    t0 = time.time()
    cand, dg, rama = load_target(pdb)
    sp = Space(pdb, cand, dg, rama)
    xors = [xor_index(sp.M, sp.q, k) for k in range(sp.q)]
    out = dict(pdb=pdb, n=sp.n, fold=int(sp.fold), S=sp.S, q=sp.q, M=sp.M,
               F=sp.F, alpha=ALPHA, T=TEMP, layers=LAYERS, iters=ITERS, lr=LR,
               gamma_gap=sp.gamma_gap, members=sp.members.tolist(), arms=[],
               seg_sizes=[int(len(s)) for s in sp.segs])

    # ---------------- ORACLE diagnostics (labelled; they choose nothing) -------------
    rr_all = sp.oracle_rmsd_all()
    parent_cfg = np.array([f * (sp.M - 1) // (sp.F - 1) for f in range(sp.F)], int)
    assert np.all(sp.digits[parent_cfg] == np.arange(sp.F)[:, None]), "parent indexing"
    nat_pair = float(sp._pair_energy(np.asarray(cand.nat_ca, float)[None])[0])
    nat_pr = I.project(np.asarray(cand.nat_ca, float), sp.seq, sp.fold)
    nat_pair_pr = float(sp._pair_energy(np.asarray(nat_pr["ca"], float)[None])[0])
    ordE = np.argsort(sp.E, kind="stable")
    rank_of_best = int(np.nonzero(ordE == int(np.argmin(rr_all)))[0][0])
    # S29-L4 hole (a): the matched order-statistic null for D1.
    spn = Space(pdb, cand, dg, rama, scramble=True)
    rr_null = spn.oracle_rmsd_all()
    out["oracle"] = dict(
        scrambled_best_rmsd_cloud=float(rr_null.min()),
        scrambled_mean_rmsd=float(rr_null.mean()),
        scrambled_argminE_rmsd=float(rr_null[int(np.argmin(spn.E))]),
        native_pair_E_projected=nat_pair_pr,
        native_pair_E_projected_pctile=float((sp.E_pair < nat_pair_pr).mean()),
        native_projected_rmsd=float(I.ca_rmsd(np.asarray(nat_pr["ca"], float), cand.nat_ca)),
        best_chimera_rmsd_cloud=float(rr_all.min()),
        best_parent_rmsd_cloud=float(rr_all[parent_cfg].min()),
        parents_rmsd_cloud=[float(v) for v in rr_all[parent_cfg]],
        pool_best_rmsd=float(np.nanmin(np.asarray(cand.oracle_rr, float))),
        space_mean_rmsd=float(rr_all.mean()),
        best_chimera_E_pctile=float(rank_of_best / max(sp.M - 1, 1)),
        native_pair_E=nat_pair,
        native_pair_E_pctile=float((sp.E_pair < nat_pair).mean()),
        E_argmin_rmsd=float(rr_all[int(ordE[0])]),
        corr_E_rmsd=float(np.corrcoef(sp.E, rr_all)[0, 1]),
    )
    emit(sp, "ORACLE_best_chimera", sp.CA[int(np.argmin(rr_all))], out, chain)

    # ---------------- exact classical ladder over H_diag ----------------------------
    emit(sp, "EXACT_argmin", sp.CA[int(ordE[0])], out, chain, dict(m=1))
    for m in (8, 75):
        emit(sp, f"EXACT_top{m}", uniform_topm(sp, m, f"top{m}"), out, chain, dict(m=m))

    # ---------------- the VQE arms ---------------------------------------------------
    circ = Q.StatevectorCircuit(sp.q, LAYERS)
    # ADDENDUM 2 / S29-L15: the Gamma grid comes FIRST and carries lane T's gate-1 quantity,
    # the total variation of the sampled distribution against the Gamma = 0 twin.  Full
    # readouts for the two registered couplings (0 and 1 x gamma_gap); R1/R2 only for the
    # two extra grid points, which exist to place the gate rather than to emit an endpoint.
    vqe_cells = [("VQE_g0_s0", circ, 0.0, 0),
                 ("VQE_g1_s0", circ, sp.gamma_gap, 0),
                 ("VQE_g05_s0", circ, 0.5 * sp.gamma_gap, 0),
                 ("VQE_g2_s0", circ, 2.0 * sp.gamma_gap, 0),
                 ("VQE_g1_s1", circ, sp.gamma_gap, 1),
                 ("VQE_prod_s0", ProductCircuit(sp.q, LAYERS), sp.gamma_gap, 0)]
    m_ref, psi_ref, p_g0 = None, None, None
    out["tv_gate"] = {}
    for name, cc, gam, sd in vqe_cells:
        tt = time.time()
        res = train(cc, sp.E, ALPHA, TEMP, gam, xors, sd)
        ro = tail_readouts(sp, res["p"], name)
        gate = H.gate_set_equality(sp.E, ro["tail_idx"], ro["m"])
        if name == "VQE_g0_s0" and sd == 0:
            p_g0 = res["p"]
        if p_g0 is not None and sd == 0:
            out["tv_gate"][name] = dict(gamma=float(gam),
                                        tv_vs_g0=total_variation(res["p"], p_g0),
                                        m=ro["m"], entropy=ro["entropy"], pr=ro["pr"])
        if name == "VQE_g1_s0":
            m_ref = ro["m"]
            psi_ref = res["psi"]
            out["p_vqe_entropy"] = ro["entropy"]
        meta = dict(m=ro["m"], cvar=ro["cvar"], entropy=ro["entropy"], ess=ro["ess"],
                    pr=ro["pr"], r3_mass=ro["r3_mass"], r3_k=ro["r3_k"],
                    gamma=float(gam), seed=int(sd), F_final=res["F"], F_first=res["hist"][0],
                    mixer=res["mixer"], tail_is_prefix=bool(gate["subset_of_energy_prefix"]),
                    tail_equals_topm=bool(gate["equality"]),
                    # lane T's M5: the CVaR term is exactly constant along every simplex
                    # direction above the VaR; the mixer term is not a function of p at all.
                    cvar_flat_frac=float((sp.M - ro["m"] - 1) / max(sp.M - 1, 1)),
                    secs=float(time.time() - tt))
        for R in (("R1", "R2", "R3") if name in ("VQE_g0_s0", "VQE_g1_s0", "VQE_g1_s1",
                                                 "VQE_prod_s0") else ("R1", "R2")):
            emit(sp, f"{name}|{R}", ro[R], out, chain, dict(meta, readout=R))
        if name == "VQE_g1_s0":
            # S29-L4 hole (c): random weights on the SAME set at the SAME participation ratio
            for d, Cw in enumerate(random_weight_control(sp, ro["tail_idx"], ro["pr"], name)):
                emit(sp, f"RANDW_s0|d{d}", Cw, out, chain,
                     dict(m=ro["m"], pr=ro["pr"], draw=d, readout="R2rand"))
        if verbose:
            print(f"    {name}: m={ro['m']} F {res['hist'][0]:.3f}->{res['F']:.3f} "
                  f"({time.time() - tt:.0f}s)", flush=True)

    # untrained: theta_0 of seed 0, no training (the mandatory control)
    rng = np.random.default_rng(0)
    th0 = rng.normal(0.0, 0.6, circ.n_params())
    p0 = circ.probs(th0)
    ro = tail_readouts(sp, p0, "UNTRAINED")
    for R in ("R1", "R2", "R3"):
        emit(sp, f"UNTRAINED_s0|{R}", ro[R], out, chain,
             dict(m=ro["m"], entropy=ro["entropy"], readout=R))

    # ---------------- the diagonalised equivalent -----------------------------------
    for tag, gam in (("GS_g1", sp.gamma_gap), ("GS_g0", 0.0)):
        w0, g = ground_state(sp.E, gam, sp.q, xors)
        ro = tail_readouts(sp, g ** 2, tag)
        for R in ("R1", "R2"):
            emit(sp, f"{tag}|{R}", ro[R], out, chain,
                 dict(m=ro["m"], gs_energy=w0, gamma=float(gam), readout=R))
        out[f"{tag}_energy"] = w0
    # reachability: overlap of the trained state with the exact ground state (no retraining)
    _, g1 = ground_state(sp.E, sp.gamma_gap, sp.q, xors)
    out["overlap_vqe_gs"] = float((psi_ref / np.linalg.norm(psi_ref) @ g1) ** 2)

    # ---------------- classical ensembles at matched entropy ------------------------
    # lane T's M6 (S29-L15 Q1), in its strongest form: the EXACT minimiser of the same
    # objective at Gamma = 0, with no circuit and no optimiser.  A quantum stage that cannot
    # beat this contributes nothing at Gamma = 0 by construction.
    pstar, tstar = cvar_optimal_law(sp.E, ALPHA, TEMP)
    ro = tail_readouts(sp, pstar, "CVAROPT")
    out["cvaropt"] = dict(t_star=tstar, m=ro["m"], entropy=ro["entropy"], pr=ro["pr"],
                          tv_vs_vqe_g0=total_variation(pstar, p_g0) if p_g0 is not None else None)
    for R in ("R1", "R2"):
        emit(sp, f"CVAROPT|{R}", ro[R], out, chain, dict(m=ro["m"], readout=R))

    pg = gibbs(sp.E, TEMP)
    ro = tail_readouts(sp, pg, "GIBBS_T1")
    for R in ("R1", "R2"):
        emit(sp, f"GIBBS_T1|{R}", ro[R], out, chain, dict(m=ro["m"], readout=R))
    pm, Tm = gibbs_matched_entropy(sp.E, out["p_vqe_entropy"])
    ro = tail_readouts(sp, pm, "GIBBS_match")
    out["gibbs_matched_T"] = float(Tm)
    for R in ("R1", "R2"):
        emit(sp, f"GIBBS_match|{R}", ro[R], out, chain,
             dict(m=ro["m"], T=float(Tm), readout=R))

    # ---------------- simulated annealing at a matched budget -----------------------
    sa = anneal(sp.E, sp.q, sp.M, 0)
    ro = tail_readouts(sp, sa["p"], "SA")
    out["sa_distinct"] = sa["n_distinct"]
    for R in ("R1", "R2"):
        emit(sp, f"SA_s0|{R}", ro[R], out, chain,
             dict(m=ro["m"], budget=sa["budget"], n_distinct=sa["n_distinct"], readout=R))

    # ---------------- the permuted posterior ----------------------------------------
    perm_E = permuted_energy(sp)
    # S29-L4 hole (e): the permuted arm's spectrum beside the real one, so "PERM did not
    # reproduce it" cannot be confounded by PERM having a flatter energy landscape.
    out["spectrum"] = dict(real_sd=float(sp.E.std()), real_range=float(np.ptp(sp.E)),
                           perm_sd=float(perm_E.std()), perm_range=float(np.ptp(perm_E)),
                           real_gap75=float(np.ptp(np.sort(sp.E)[:75])),
                           perm_gap75=float(np.ptp(np.sort(perm_E)[:75])),
                           perm_seed="stable_rng(pdb,'perm',salt='s29X')")
    ordP = np.argsort(perm_E, kind="stable")
    emit(sp, "PERM_EXACT_top75", weighted_average(sp.CA[ordP[:75]], None, pdb, "permtop"),
         out, chain, dict(m=75))
    spp = _ShadowSpace(sp, perm_E)
    resP = train(circ, perm_E, ALPHA, TEMP, sp.gamma_gap, xors, 0)
    roP = tail_readouts(spp, resP["p"], "PERM_VQE")
    for R in ("R1", "R2"):
        emit(sp, f"PERM_VQE_g1_s0|{R}", roP[R], out, chain, dict(m=roP["m"], readout=R))

    # ---------------- the size-matched classical control ----------------------------
    if m_ref:
        emit(sp, "EXACT_topm_matched", uniform_topm(sp, m_ref, "topm"), out, chain,
             dict(m=int(m_ref)))

    out["secs"] = float(time.time() - t0)
    return out


class _ShadowSpace:
    """The same geometry with a different energy vector (for the permuted-posterior arm)."""

    def __init__(self, sp: Space, E: np.ndarray):
        self.pdb, self.CA, self.E, self.M = sp.pdb, sp.CA, np.asarray(E, float), sp.M


def permuted_energy(sp: Space) -> np.ndarray:
    """H_diag recomputed with the posterior's rows permuted WITHIN sequence separation.

    Native-free, structure-preserving in the marginal sense: each pair keeps a real predicted
    distribution, but not its own.  The control for "is the specific posterior doing the work".
    """
    rng = SD.stable_rng(sp.pdb, "perm", salt=SALT)
    prob = np.asarray(sp.dg["prob"], float).copy()
    sep = sp.j - sp.i
    for s in np.unique(sep):
        k = np.nonzero(sep == s)[0]
        if len(k) > 1:
            prob[k] = prob[k[rng.permutation(len(k))]]
    lp = -np.log(prob + EPS_P)
    out = np.empty(sp.M)
    for a in range(0, sp.M, 8192):
        D = I.pair_dists(sp.CA[a:a + 8192], sp.i, sp.j)
        b = np.digitize(D, PRD.BIN_EDGES)
        out[a:a + 8192] = lp[np.arange(lp.shape[0])[None, :], b].sum(1)
    return out + sp.E_rama


# ===================================================================== lane D's meter
def cost_nll(W, ctx):
    """The pair half of H_diag as a cost function on CA clouds, for `s29/s29_D_cost_audit.py`.

    Lower is better.  The Ramachandran half is a function of torsions and cannot be evaluated
    from a cloud, so the meter sees the pair term only; the ledger entry says so.
    """
    W = np.asarray(W, float)
    if W.ndim == 2:
        W = W[None]
    dg = ctx.dg
    i, j = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
    lp = -np.log(np.asarray(dg["prob"], float) + EPS_P)
    D = I.pair_dists(W, i, j)
    b = np.digitize(D, PRD.BIN_EDGES)
    return lp[np.arange(lp.shape[0])[None, :], b].sum(1)


# ========================================================================== analysis
def analyse(pdbs: Optional[Sequence[str]] = None) -> Dict:
    pdbs = list(pdbs or PROBE_TARGETS)
    rows = []
    for p in pdbs:
        f = target_path(p)
        if os.path.exists(f):
            rows.append(json.load(open(f)))
    if not rows:
        return dict(n=0)
    names = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(names)
    # D1, with the matched order-statistic null: recombination's ORACLE value is the gap to
    # the SCRAMBLED space of identical cardinality, never the gap to the eight parents.
    d1 = {}
    for a, b, lab in (("best_chimera_rmsd_cloud", "scrambled_best_rmsd_cloud",
                       "D1 chimera ORACLE best - SCRAMBLED ORACLE best (cloud)"),
                      ("best_chimera_rmsd_cloud", "best_parent_rmsd_cloud",
                       "D1u chimera ORACLE best - best PARENT (cloud, UNMATCHED order statistic)"),
                      ("best_chimera_rmsd_cloud", "pool_best_rmsd",
                       "D1p chimera ORACLE best - POOL ORACLE best K=500 (cloud)")):
        va = np.asarray([r["oracle"][a] for r in rows], float)
        vb = np.asarray([r["oracle"][b] for r in rows], float)
        r = ST.compare(va, vb, folds=folds, names=names, label=lab, seed_parts=("s29X",))
        d1[lab] = ST.fmt(r)
    prod = production_chain(names)
    idx = {r["pdb"]: {a["arm"]: a for a in r["arms"]} for r in rows}
    arms = sorted({a for r in rows for a in idx[r["pdb"]]})
    out = dict(n=len(rows), pdbs=names, arms={}, fmt={}, oracle={}, meta={})
    for basis in ("rmsd_cloud", "rmsd_chain"):
        for a in arms:
            v = [idx[p][a].get(basis) for p in names]
            if any(x is None for x in v):
                continue
            out["arms"][f"{a}|{basis}"] = dict(mean=float(np.mean(v)),
                                               median=float(np.median(v)),
                                               per_target={p: float(x) for p, x in zip(names, v)})
    for k in ("best_chimera_rmsd_cloud", "best_parent_rmsd_cloud", "pool_best_rmsd",
              "space_mean_rmsd", "best_chimera_E_pctile", "native_pair_E_pctile",
              "corr_E_rmsd", "E_argmin_rmsd", "scrambled_best_rmsd_cloud",
              "scrambled_mean_rmsd", "scrambled_argminE_rmsd",
              "native_pair_E_projected_pctile", "native_projected_rmsd"):
        v = [r["oracle"][k] for r in rows]
        out["oracle"][k] = dict(mean=float(np.mean(v)), median=float(np.median(v)),
                                per_target={p: float(x) for p, x in zip(names, v)})
    for k in ("q", "M", "gamma_gap", "overlap_vqe_gs", "gibbs_matched_T", "sa_distinct",
              "p_vqe_entropy", "secs"):
        v = [r.get(k) for r in rows if r.get(k) is not None]
        if v:
            out["meta"][k] = dict(mean=float(np.mean(v)), min=float(np.min(v)),
                                  max=float(np.max(v)))
    # GATE 1 (lane T, S29-L15): TV of the sampled distribution against the Gamma = 0 twin.
    out["tv_gate"] = {}
    for nm in ("VQE_g05_s0", "VQE_g1_s0", "VQE_g2_s0"):
        v = [r["tv_gate"][nm]["tv_vs_g0"] for r in rows if nm in r.get("tv_gate", {})]
        if v:
            out["tv_gate"][nm] = dict(mean=float(np.mean(v)), median=float(np.median(v)),
                                      min=float(np.min(v)), max=float(np.max(v)),
                                      n_above_gate=int(sum(x > TV_GATE for x in v)),
                                      n=len(v), gate=TV_GATE,
                                      per_target={p: float(x) for p, x in zip(names, v)})

    def cmp(a, b, basis, label):
        va = [idx[p][a].get(basis) for p in names]
        vb = ([prod[p][basis] for p in names] if b == "PRODUCTION"
              else [idx[p][b].get(basis) for p in names])
        if any(x is None for x in va) or any(x is None for x in vb):
            return
        r = ST.compare(np.asarray(va, float), np.asarray(vb, float), folds=folds,
                       names=names, label=label, seed_parts=("s29X",))
        out["fmt"][label] = ST.fmt(r)
        out.setdefault("compare", {})[label] = {k: r[k] for k in
                                                ("effect", "median_effect", "se", "mde",
                                                 "effect_over_mde", "ci95_iid", "ci95_fold",
                                                 "n_better", "n_worse", "verdict")}

    for basis in ("rmsd_cloud", "rmsd_chain"):
        b = basis.split("_")[1]
        for a in arms:
            cmp(a, "PRODUCTION", basis, f"P1 {a} - PRODUCTION ({b})")
        cmp("VQE_g1_s0|R2", "SA_s0|R2", basis, f"P2 VQE R2 - SA R2 ({b})")
        cmp("VQE_g1_s0|R2", "VQE_g1_s0|R1", basis, f"P3 VQE R2 - VQE R1 ({b})")
        # S29-L4 hole (c): the PR-matched random-weight control for P3, averaged over draws.
        rw = [a for a in arms if a.startswith("RANDW_s0|")]
        if rw:
            va = [idx[p]["VQE_g1_s0|R2"].get(basis) for p in names]
            vb = [float(np.mean([idx[p][a][basis] for a in rw if a in idx[p]])) for p in names]
            if not any(x is None for x in va):
                r = ST.compare(np.asarray(va, float), np.asarray(vb, float), folds=folds,
                               names=names, label=f"P3c VQE R2 - PR-matched random weights ({b})",
                               seed_parts=("s29X",))
                out["fmt"][r["label"]] = ST.fmt(r)
                out.setdefault("compare", {})[r["label"]] = {
                    k: r[k] for k in ("effect", "se", "mde", "effect_over_mde",
                                      "ci95_fold", "verdict")}
        cmp("VQE_g1_s0|R2", "GIBBS_match|R2", basis, f"P4 VQE R2 - GIBBS(matched) R2 ({b})")
        cmp("VQE_g1_s0|R2", "VQE_g0_s0|R2", basis, f"P5 VQE gamma - VQE gamma=0, R2 ({b})")
        cmp("VQE_g1_s0|R2", "UNTRAINED_s0|R2", basis, f"C-untrained VQE R2 - UNTRAINED R2 ({b})")
        cmp("VQE_g1_s0|R2", "VQE_prod_s0|R2", basis, f"C-product VQE R2 - PRODUCT R2 ({b})")
        cmp("VQE_g1_s0|R2", "VQE_g1_s1|R2", basis, f"C-seed VQE R2 s0 - s1 ({b})")
        cmp("VQE_g1_s0|R2", "GS_g1|R2", basis, f"C-diagonalised VQE R2 - GS R2 ({b})")
        cmp("VQE_g0_s0|R2", "CVAROPT|R2", basis, f"M6 VQE(gamma=0) R2 - CVaR-OPTIMAL LAW R2 ({b})")
        cmp("VQE_g1_s0|R2", "CVAROPT|R2", basis, f"M6g VQE(gamma) R2 - CVaR-OPTIMAL LAW R2 ({b})")
        cmp("VQE_g1_s0|R2", "PERM_VQE_g1_s0|R2", basis, f"C-perm VQE R2 - PERM VQE R2 ({b})")
        cmp("VQE_g1_s0|R2", "EXACT_topm_matched", basis, f"C-ordstat VQE R2 - EXACT top-m ({b})")
        cmp("EXACT_top75", "PRODUCTION", basis, f"D2 EXACT top-75 - PRODUCTION ({b})")
    out["fmt"].update(d1)
    for k in ("pr", "r3_mass", "m"):
        v = [idx[p]["VQE_g1_s0|R2"].get(k) for p in names if k in idx[p]["VQE_g1_s0|R2"]]
        if v:
            out["meta"][f"vqe_{k}"] = dict(mean=float(np.mean(v)), min=float(np.min(v)),
                                           max=float(np.max(v)))
    for nm in ("VQE_g1_s0|R1", "VQE_g1_s0|R2", "EXACT_top75", "ORACLE_best_chimera"):
        for s in ("rg", "bond"):
            v = [idx[p][nm][s] for p in names if nm in idx[p] and s in idx[p][nm]]
            if v:
                out["meta"][f"{nm}|{s}"] = float(np.mean(v))
    path = os.path.join(RESULTS, "s29_X_probe.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return out


def production_chain(pdbs: Sequence[str]) -> Dict[str, Dict[str, float]]:
    """Production's own rows: `s27/results/chain_rows.jsonl :: DIS`, seed 0 (the anchor)."""
    out = {}
    with open(os.path.join(ROOT, "s27", "results", "chain_rows.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("config") == "DIS" and int(r.get("seed", 0)) == 0 and r["pdb"] in pdbs:
                out[r["pdb"]] = dict(rmsd_cloud=float(r["rmsd_cloud"]),
                                     rmsd_chain=float(r["rmsd_chain"]))
    return out


# ========================================================================== selftest
def selftest(pdb: str = "2P5H") -> None:
    cand, dg, rama = load_target(pdb)
    sp = Space(pdb, cand, dg, rama)
    print(f"  {pdb}: n={sp.n} S={sp.S} q={sp.q} M={sp.M} gamma_gap={sp.gamma_gap:.4f}")
    # 1. a parent configuration reproduces its member's rebuild exactly
    for f in range(sp.F):
        c = f * (sp.M - 1) // (sp.F - 1)
        rb = I.build_ca(cand.PHI[sp.members[f]], cand.PSI[sp.members[f]])
        assert np.allclose(sp.CA[c], rb, atol=0, rtol=0), f"parent {f} mismatch"
    print("  parents reproduce their members' rebuilds bit-exactly: OK")
    # 2. gamma = 0 objective equals core.quantum.free_energy
    circ = Q.StatevectorCircuit(sp.q, LAYERS)
    th = np.random.default_rng(3).normal(0, 0.6, circ.n_params())
    f0, g0, _, _, _ = Q.free_energy(circ, th, sp.E, ALPHA, TEMP)
    f1, g1, *_ = f_and_grad(circ, th, sp.E, ALPHA, TEMP, 0.0, [])
    assert abs(f0 - f1) < 1e-12 and np.max(np.abs(g0 - g1)) < 1e-12
    print(f"  gamma=0 == core.quantum.free_energy: |dF| {abs(f0-f1):.2e} "
          f"|dg|max {np.max(np.abs(g0-g1)):.2e}")
    # 3. uniform weighted_average == I.coordinate_average
    W = sp.CA[:40]
    C0, _ = I.coordinate_average(W)
    C1 = weighted_average(W, None, pdb, "selftest")
    assert np.max(np.abs(C0 - C1)) < 1e-12
    print(f"  uniform weighted_average == coordinate_average: {np.max(np.abs(C0-C1)):.2e}")
    # 4. mixer gradient against finite differences
    xors = [xor_index(sp.M, sp.q, k) for k in range(sp.q)]
    gam = sp.gamma_gap
    f, g, *_ = f_and_grad(circ, th, sp.E, ALPHA, TEMP, gam, xors)
    h = 1e-6
    fd = np.empty(3)
    for k in range(3):
        tp = th.copy(); tp[k] += h
        tm = th.copy(); tm[k] -= h
        fd[k] = (f_and_grad(circ, tp, sp.E, ALPHA, TEMP, gam, xors)[0]
                 - f_and_grad(circ, tm, sp.E, ALPHA, TEMP, gam, xors)[0]) / (2 * h)
    print(f"  exact vs FD gradient (3 coords): max |diff| {np.max(np.abs(fd - g[:3])):.2e}")
    assert np.max(np.abs(fd - g[:3])) < 1e-5


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--targets", type=str, default="")
    ap.add_argument("--no-chain", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    pdbs = [s for s in a.targets.split(",") if s] or list(PROBE_TARGETS)
    if not a.targets:
        # lane T's ladder is run CHEAPEST FIRST (S29-L15): the register size sets the cost,
        # so the gate-1 quantity arrives on the small targets within minutes.
        tg = {t["pdb"]: t["n"] for t in I.targets()}
        pdbs = sorted(pdbs, key=lambda p: (min(math.ceil(tg[p] / SEG_LEN), SEG_MAX), p))
    if a.limit:
        pdbs = pdbs[:a.limit]
    if a.run:
        for pdb in pdbs:
            f = target_path(pdb)
            if os.path.exists(f) and not a.force:
                print(f"  {pdb}: cached", flush=True)
                continue
            t0 = time.time()
            out = run_target(pdb, chain=not a.no_chain)
            tmp = f + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(out, fh, indent=1,
                          default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
            os.replace(tmp, f)
            print(f"  {pdb}: q={out['q']} M={out['M']} arms={len(out['arms'])} "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
    if a.analyse:
        o = analyse(pdbs)
        for k in sorted(o.get("fmt", {})):
            print(o["fmt"][k])
            print()


if __name__ == "__main__":
    main()
