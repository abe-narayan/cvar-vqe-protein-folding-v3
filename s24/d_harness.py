"""S24 / LANE D -- THE INTEGRATION HARNESS. Any candidate pool -> CVaR-VQE + classical control.

WHAT THIS IS FOR
================
Lanes B and C are building candidate GENERATORS. Their pools must not be evaluated by
whatever ad-hoc loop each lane happens to write, because the single most repeated error in
this project's memory is `control-must-match-the-operators-space` -- three instances in two
sprints. This module is the one path both pools go through, and it makes the control
STRUCTURALLY impossible to get wrong: the quantum arm and the classical arm are handed the
SAME `Candidates` object, the SAME energy vector, the SAME readout operator and the SAME
frame. Nothing in the API lets a caller vary one arm's inputs without varying the other's.

    from s24.d_harness import Candidates, score_shipped, run_both, aggregate
    cand = Candidates.from_universe("1L2Y")            # or .from_arrays(W, PHI, PSI, ...)
    row  = run_both(cand, score_shipped(cand), alpha=0.15, T=0.5, seed=0)
    agg  = aggregate(rows, folds)                       # SE, MDE, iid + fold CI, W/L, worst

THE THREE PILLARS, HELD (BRIEF s0)
==================================
* The quantum arm is a GENUINE CVaR-VQE. `core.quantum.run_cvar_vqe` -> exact
  `StatevectorCircuit` (layers x [RY on every wire, CNOT chain, ring closure]), Adam on the
  EXACT parameter-shift gradient, genuine CVaR of the trained distribution against a genuine
  diagonal Hamiltonian `H|i> = E_i|i>` over the candidate-identity register (`s22.qcand_lib`
  Encoding, reused as-is, not reimplemented). K <= 2**n candidates means the register is
  exactly enumerable, so there is no shot noise and no sampling approximation anywhere.
* `H_Legacy` is the genuine 11-term potential at `DEFAULT_WEIGHTS` via
  `s16.energy_lib.legacy_components_of_windows` + `legacy_total_from`. Never fitted.
* `H_AMBER` is genuine ff14SB/GBn2 through the existing OpenMM path (`s20.qb2_lib.AmberSP`,
  which asserts bit-exactness against `core.amber.refine_coords` before returning a number).
  It lives behind `LOCK_AMBER` and it is NEVER a default stage -- see `d_hamiltonians.py`.

THE STANDING DETERMINATION THIS HARNESS ENCODES (see `d_setequality_proof.py`)
=============================================================================
It is a THEOREM ABOUT THE ALGORITHM, not a property of this pool, and the exact statement
matters (the coordinator corrected a looser first draft of mine; this is the agreed wording):

  > The realised CVaR tail's support is always a SUBSET of an initial prefix of the energy
  > order, and equals that prefix exactly when every state in the prefix carries positive
  > probability. p_theta can delete a member; it can never add one outside the classical
  > top-m.

`cvar_from_probs` (core/quantum.py:352-354) allocates the alpha mass along the ENERGY order
via a non-decreasing exclusive prefix sum, and `take`'s clip has a SECOND conjunct
`p[order] > 0` that punches zero-probability states out of that prefix; `tail_indices`
(core/quantum.py:229-272) takes no probability vector at all. Measured over 3,888 adversarial
cells: 0 subset-hood violations, 0 holes that were not exactly zero-probability states, and
2268/2268 exact EQUALITY on every full-support cell. A different candidate manifold cannot
break any of it. Therefore:

  * `gate_set_equality` below is a SOUNDNESS GATE ON THE IMPLEMENTATION -- it catches a
    mis-wired label, a padding leak or an energy mismatch between the arms. It is NOT
    evidence about the pool and it is never reported as a finding. It asserts the THEOREM
    (subset-hood) and reports EQUALITY separately, because a trained RY/CNOT state has
    generic angles and therefore full support, which is why s22 saw equality on 2016/2016
    cells -- an empirical regime, not the theorem.
  * The trained state's ONLY remaining degrees of freedom over the SET are `m`, the rung of
    the classical top-m ladder its alpha mass reaches, and (at exact zeros only) WHICH
    prefix members it DELETES. So the harness reports the classical control at BOTH the VQE's
    realised m (the size-matched bar, s23's convention) and at the production rung m=75, and
    the difference between those two IS the whole of what the circuit can contribute to
    selection. Stating that up front is what stops a future reader mistaking an m-ladder
    movement for a quantum effect.

DISCIPLINE BUILT IN, NOT LEFT TO THE CALLER
===========================================
* BASIS. Every RMSD here is POINT-CLOUD (`I.ca_rmsd` of a coordinate average against
  `nat_ca`), tagged `basis="point_cloud"` in every row. Built-chain numbers are produced only
  by `readout_projected` and are tagged `basis="built_chain"`. `aggregate` REFUSES to combine
  rows of different basis. The gap is 0.156 A of pure operator choice.
* LABELS. Every RMSD carries `label` in {ORACLE, ACHIEVABLE, PRODUCTION}. Scores are
  native-free; `nat_ca` and `rr` enter only through `_oracle_*` fields, which are named so.
* NULLS. `arm_random` is the zero-information control matched in the OPERATOR'S space (random
  m-subsets put through the identical coordinate average), not a degenerate one. Best-of-K
  quantities are reported as the distribution of the MAXIMUM.
* SEEDS. `s15.seed.stable_rng` only.
* WRITES. `write()` is atomic (tmp + os.replace) and the `_COMPLETE` sidecar demands the FULL
  key set, not a row count.

SCOPE LIMITS, DECLARED
======================
* No AMBER on this path. `run_both` never imports `core.amber`; the Hamiltonian-disagreement
  work is in `d_hamiltonians.py`, which takes `LOCK_AMBER` and announces it.
* No refinement/repair arm. s23 L11 measured 17 of 17 restrained/unrestrained settings at or
  worse than no repair, bottoming out at the no-op. `readout_projected` exists as a labelled
  diagnostic and is off by default; any arm that uses it must be benchmarked against NO
  projection on matched candidates, and `run_both` will emit both when asked.
* This module runs NO directional experiment of its own and therefore registers no Rule-0
  fork list. The lanes that USE it must pre-register theirs.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s16 import energy_lib as EL           # noqa: E402
from s22 import qcand_lib as QC            # noqa: E402
from core import quantum as Q              # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

SALT = "s24laneD"
#: production rung. `s21/results/poolgap.json` avg_75 = 3.0483380938795324 = the incumbent.
M_PROD = 75
#: the register ceiling. StatevectorCircuit is exact but parameter shift costs 2*P full
#: simulations of 2**n amplitudes; beyond this the harness refuses rather than silently
#: switching to an approximate path.
MAX_QUBITS = 13

BASES = ("point_cloud", "built_chain")
LABELS = ("ORACLE", "ACHIEVABLE", "PRODUCTION")


# ============================================================== the candidate container
@dataclass
class Candidates:
    """An arbitrary candidate pool for one target, plus the ORACLE label kept apart.

    `W` is the ONLY thing the readout touches. `PHI`/`PSI` are needed for Legacy and AMBER
    and may be None for a coordinate-only generator -- in which case `score_legacy` and the
    AMBER path raise rather than silently substituting something.

    `nat_ca` and `oracle_rr` are LABELS. Nothing in this module may read them outside a
    function whose name begins `_oracle` or whose returned key begins `oracle_`.
    """
    pdb: str
    n: int
    seq: str
    fold: int
    W: np.ndarray                                   # (k, n, 3) CA, point-cloud basis
    PHI: Optional[np.ndarray] = None                # (k, n) radians
    PSI: Optional[np.ndarray] = None
    nat_ca: Optional[np.ndarray] = None             # (n, 3) ORACLE
    oracle_rr: Optional[np.ndarray] = None          # (k,)   ORACLE
    source: str = "unspecified"                     # provenance: which lane made these
    meta: Dict = field(default_factory=dict)

    def __post_init__(self):
        self.W = np.asarray(self.W, float)
        if self.W.ndim != 3 or self.W.shape[1] != self.n or self.W.shape[2] != 3:
            raise ValueError(f"W must be (k, {self.n}, 3), got {self.W.shape}")
        if len(self.W) < 2:
            raise ValueError("a pool needs at least 2 candidates")
        if not np.isfinite(self.W).all():
            raise ValueError("W contains non-finite coordinates")
        for nm in ("PHI", "PSI"):
            v = getattr(self, nm)
            if v is not None:
                v = np.asarray(v, float)
                if v.shape != (len(self.W), self.n):
                    raise ValueError(f"{nm} must be ({len(self.W)}, {self.n}), got {v.shape}")
                setattr(self, nm, v)
        if self.nat_ca is not None:
            self.nat_ca = np.asarray(self.nat_ca, float)
        if self.oracle_rr is not None:
            self.oracle_rr = np.asarray(self.oracle_rr, float)

    @property
    def k(self) -> int:
        return int(len(self.W))

    @property
    def has_torsions(self) -> bool:
        return self.PHI is not None and self.PSI is not None

    # -- constructors -------------------------------------------------------
    @classmethod
    def from_universe(cls, pdb: str, k: int = 500, source: str = "incumbent_blosum_top500"):
        """The SHIPPED K=500 BLOSUM retrieval pool. The reference manifold."""
        t = {x["pdb"]: x for x in I.targets()}[pdb]
        u = I.load_univ(pdb)
        p = I.pool_idx(u, k=k)
        return cls(pdb=pdb, n=int(t["n"]), seq=t["seq"], fold=int(t["fold"]),
                   W=u["W"][p], PHI=u["PHI"][p], PSI=u["PSI"][p],
                   nat_ca=u["nat_ca"], oracle_rr=u["rr"][p], source=source,
                   meta=dict(universe_idx=p.tolist()))

    @classmethod
    def from_uniform_library(cls, pdb: str, k: int = 500, seed: int = 0,
                             source: str = "uniform_library_draw"):
        """A blind uniform draw from the target's leakage-safe window universe.

        The zero-cost stand-in generator the coordinator used for L2(d). Its role here is as
        a MATERIALLY DIFFERENT MANIFOLD that already exists: it lets the harness be exercised
        end-to-end, and the set-equality gate be checked on a pool with a visibly different
        energy structure, before either real generator ships.
        """
        t = {x["pdb"]: x for x in I.targets()}[pdb]
        u = I.load_univ(pdb)
        nw = len(u["W"])
        rng = SD.stable_rng(pdb, "unifpool", int(seed), int(k), salt=SALT)
        p = rng.choice(nw, size=min(k, nw), replace=False)
        return cls(pdb=pdb, n=int(t["n"]), seq=t["seq"], fold=int(t["fold"]),
                   W=u["W"][p], PHI=u["PHI"][p], PSI=u["PSI"][p],
                   nat_ca=u["nat_ca"], oracle_rr=u["rr"][p], source=source,
                   meta=dict(universe_idx=p.tolist(), seed=int(seed)))

    @classmethod
    def from_arrays(cls, pdb: str, W, PHI=None, PSI=None, source: str = "generator",
                    meta: Optional[Dict] = None):
        """THE ENTRY POINT FOR LANES B AND C. Target metadata and the native label come from
        the instrument, so a generator cannot accidentally supply its own native."""
        t = {x["pdb"]: x for x in I.targets()}[pdb]
        u = I.load_univ(pdb)
        return cls(pdb=pdb, n=int(t["n"]), seq=t["seq"], fold=int(t["fold"]),
                   W=np.asarray(W, float), PHI=PHI, PSI=PSI, nat_ca=u["nat_ca"],
                   oracle_rr=None, source=source, meta=dict(meta or {}))

    def merge(self, other: "Candidates", source: Optional[str] = None) -> "Candidates":
        """Concatenate two pools for the SAME target -- the §17 union operator.

        Used when a lane wants "incumbent + generated" scored by one functional. The result
        carries a `split` array so the share of the retained set taken from each half can be
        reported, which is the quantity L2(d) found to be 0.355 for a free stand-in source.
        """
        if other.pdb != self.pdb or other.n != self.n:
            raise ValueError("merge requires the same target")
        both_t = self.has_torsions and other.has_torsions
        return Candidates(
            pdb=self.pdb, n=self.n, seq=self.seq, fold=self.fold,
            W=np.concatenate([self.W, other.W], 0),
            PHI=np.concatenate([self.PHI, other.PHI], 0) if both_t else None,
            PSI=np.concatenate([self.PSI, other.PSI], 0) if both_t else None,
            nat_ca=self.nat_ca,
            oracle_rr=(np.concatenate([self.oracle_rr, other.oracle_rr])
                       if self.oracle_rr is not None and other.oracle_rr is not None else None),
            source=source or f"merge({self.source}+{other.source})",
            meta=dict(split=np.r_[np.zeros(self.k, int), np.ones(other.k, int)].tolist(),
                      sources=[self.source, other.source], sizes=[self.k, other.k]))


# ==================================================================== the energy channel
def score_shipped(cand: Candidates) -> np.ndarray:
    """The shipped distogram Bayes-risk score. NATIVE-FREE. Lower is better.

    This is the functional the incumbent 3.0483 A pipeline ranks with, so it is the default
    H for both arms and the thing a generator's samples must actually be good under.
    """
    dg = I.distogram(cand.pdb, cand.seq, cand.fold)
    i, j = I.pair_index(cand.n)
    D = I.pair_dists(cand.W, i, j)
    return np.asarray(I.shipped_score(dg, D.astype(np.float32).astype(float)), float)


def score_legacy(cand: Candidates) -> np.ndarray:
    """GENUINE `H_Legacy`: the 11-term potential at `DEFAULT_WEIGHTS`. NATIVE-FREE.

    `legacy_total_from(comp)` with no `w` uses `DEFAULT_WEIGHTS` verbatim -- never fitted,
    never a regression target (BRIEF pillar 2). Requires torsions.
    """
    if not cand.has_torsions:
        raise ValueError("score_legacy needs PHI/PSI; this pool is coordinate-only")
    comp = EL.legacy_components_of_windows(cand.seq, cand.PHI, cand.PSI)
    return np.asarray(EL.legacy_total_from(comp), float)


def zrank(x) -> np.ndarray:
    """Standardised rank -- the currency `core.pipeline.quantum_stage` puts E in.

    Rank-standardising makes the CVaR temperature comparable across targets whose raw score
    scales differ, which is exactly why the deployed stage does it. It is a MONOTONE map, so
    it changes no ordering and therefore (by the set-equality theorem) no tail membership --
    only the energy GAPS the entropy term trades against.
    """
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(r.std(), 1e-12)


# ======================================================================== the readouts
def readout_uniform(cand: Candidates, idx: Sequence[int]) -> Tuple[np.ndarray, int]:
    """THE DEPLOYED READOUT: uniform coordinate average in the retained set's own medoid
    frame (`s12.instrument.coordinate_average`, S8-11's operator). Point-cloud basis.

    Uniform is not a default of convenience: s23 L8/L5 and s22 established it as the
    minimum-variance combination of exchangeable estimators with i.i.d. errors, and every
    departure measured (sharper weights, clustering, smaller m, amplitude weighting) attacks
    the variance reduction that makes averaging work.
    """
    idx = np.asarray(idx, int)
    if idx.size == 0:
        raise ValueError("empty retained set")
    C, b = I.coordinate_average(cand.W[idx])
    return C, int(b)


def readout_projected(cand: Candidates, C: np.ndarray) -> np.ndarray:
    """Ideal-geometry projection of a coordinate average. BUILT-CHAIN basis, DIAGNOSTIC ONLY.

    Off by default. Never compared against a point-cloud number.
    """
    pr = I.project(np.asarray(C, float), cand.seq, cand.fold)
    return np.asarray(pr["ca"], float)


def _oracle_rmsd(cand: Candidates, C: np.ndarray) -> float:
    """ORACLE scoring of an emitted cloud, post-hoc. The ONLY place `nat_ca` is read."""
    if cand.nat_ca is None:
        return float("nan")
    return float(I.ca_rmsd(np.asarray(C, float), cand.nat_ca))


# =========================================================================== the arms
def arm_vqe(cand: Candidates, E: np.ndarray, alpha: float = 0.15, T: float = 0.5,
            layers: int = 3, iters: int = 80, lr: float = 0.15, restarts: int = 1,
            seed: int = 0, label_perm: bool = False) -> Dict:
    """GENUINE CVaR-VQE selection over the candidate-identity register.

    `H|i> = E_i|i>` diagonal on the SUPPLIED energies -- the harness never lets this arm see
    a different functional from the control. The circuit is the exact `StatevectorCircuit`;
    the optimiser is Adam on the EXACT parameter-shift gradient (`run_cvar_vqe`), so the
    historical `baseline="tail"` CVaR-gradient defect cannot arise on this path. The tail is
    read EXACTLY off `cvar_from_probs`, not sampled.

    `label_perm=True` scrambles the candidate<->basis-index bijection with a stable RNG. The
    exact optimum cannot see the label; a finite-iteration one might (s18 had exactly this
    failure once). It is the GAUGE control and it is cheap, so it is offered rather than
    assumed away.
    """
    E = np.asarray(E, float)
    if len(E) != cand.k:
        raise ValueError(f"energies ({len(E)}) must match candidates ({cand.k})")
    n_qubits = int(math.ceil(math.log2(max(2, cand.k))))
    if n_qubits > MAX_QUBITS:
        raise ValueError(f"pool of {cand.k} needs {n_qubits} qubits, ceiling is {MAX_QUBITS}; "
                         f"prefilter the pool and SAY SO, do not silently approximate")
    lab = (QC.random_label(1 << n_qubits, cand.pdb, f"gauge{seed}", cand.k, salt=SALT)
           if label_perm else None)
    enc = QC.Encoding(E, label=lab)
    t0 = time.time()
    p, cv, Hent, circ = Q.run_cvar_vqe(enc.E, alpha, T, n=enc.n_qubits, layers=layers,
                                       iters=iters, restarts=restarts, seed=int(seed), lr=lr)
    face = QC.exact_face(enc.E, np.asarray(p, float), alpha)
    cands = QC.tail_candidates(enc, face["tail"])
    n_pad = int(face["tail_size"] - len(cands))
    ent_bits = float(-(np.asarray(p)[np.asarray(p) > 0]
                       * np.log2(np.asarray(p)[np.asarray(p) > 0])).sum())
    return dict(arm="vqe", m=int(len(cands)), cands=np.sort(cands).astype(int),
                tail_size=int(face["tail_size"]), n_pad_in_tail=n_pad,
                cvar=float(cv), entropy_bits=ent_bits, entropy_nats=float(Hent),
                ess=float(face["ess"]), alpha=float(alpha), T=float(T),
                n_qubits=int(enc.n_qubits), layers=int(layers), iters=int(iters),
                seed=int(seed), gauge_permuted=bool(label_perm),
                secs=float(time.time() - t0))


def arm_classical(cand: Candidates, E: np.ndarray, m: int) -> Dict:
    """THE CONTROL: exact classical rank, top-m, on the IDENTICAL candidates and energies.

    In a register this size exhaustive scoring costs k evaluations -- cheaper than any
    training budget -- so "greedy" and "exact argmin" coincide and that identity is stated
    rather than hidden.
    """
    E = np.asarray(E, float)
    m = int(max(1, min(int(m), len(E))))
    top = np.argsort(E, kind="stable")[:m]
    return dict(arm="classical", m=m, cands=np.sort(top).astype(int))


def arm_random(cand: Candidates, m: int, n_draws: int = 16, seed: int = 0) -> Dict:
    """ZERO-INFORMATION control, MATCHED IN THE OPERATOR'S SPACE.

    `n_draws` random m-subsets, each put through the IDENTICAL uniform coordinate average.
    Not a degenerate control (a uniform-on-the-torus draw would be a worse measure, not an
    uninformative one) -- these are real candidates from the same pool, so the only thing
    removed is the score's ordering information.

    Reports the MEAN over draws as the null, and the MAX/MIN over draws separately so a
    best-of-K comparison can be priced against the distribution of the MAXIMUM.
    """
    rng = SD.stable_rng(cand.pdb, "randnull", int(seed), int(m), salt=SALT)
    m = int(max(1, min(int(m), cand.k)))
    rs = []
    for _ in range(int(n_draws)):
        pick = rng.choice(cand.k, size=m, replace=False)
        C, _ = readout_uniform(cand, pick)
        rs.append(_oracle_rmsd(cand, C))
    rs = np.asarray(rs, float)
    return dict(arm="random", m=m, n_draws=int(n_draws),
                oracle_rmsd_mean=float(rs.mean()), oracle_rmsd_sd=float(rs.std(ddof=1)),
                oracle_rmsd_best=float(rs.min()), oracle_rmsd_worst=float(rs.max()))


# ==================================================================== the soundness gate
def gate_set_equality(E: np.ndarray, vqe_cands: np.ndarray, m: int) -> Dict:
    """SOUNDNESS GATE ON THE IMPLEMENTATION, not evidence about the pool.

    THE ASSERTION is the theorem: `SUBSET`-hood. The VQE tail must live inside the classical
    energy-ordered prefix that reaches its own worst member -- p_theta can delete, never add.
    A failure here means the harness is mis-wired: a scrambled label, a padding state that
    leaked into the tail, or the two arms seeing different energies. That is the ONLY thing
    this gate is for, and it is never quoted as a result.

    `equality` is REPORTED, not asserted. A trained RY/CNOT state has generic angles and
    therefore full support, so equality holds in practice (s22: 2016/2016 cells) -- but that
    is an empirical regime, not the theorem, and conflating the two is exactly the error the
    coordinator corrected in my first draft. If `equality` ever goes false while `pass_`
    stays true, the circuit has produced exact zeros and the harness should say so rather
    than fail.

    Comparison is VALUE-BASED, so exact ties cannot manufacture a spurious failure (project
    memory: tie-breaking leaks the pool order). Target `5H1H` has an exact top-2 tie.
    """
    E = np.asarray(E, float)
    v = np.asarray(vqe_cands, int)
    order = np.argsort(E, kind="stable")
    # the prefix that reaches the tail's own worst member -- value-based so ties are safe
    e_cut = float(E[v].max()) if v.size else float("-inf")
    prefix = order[E[order] <= e_cut + 1e-12]
    subset = bool(set(v.tolist()) <= set(prefix.tolist()))
    cls = order[:int(m)]
    a, b = np.sort(E[v]), np.sort(E[cls])
    idx_equal = bool(set(v.tolist()) == set(cls.tolist()))
    val_equal = bool(len(a) == len(b) and np.allclose(a, b, rtol=0, atol=1e-12))
    return dict(pass_=subset, subset_of_energy_prefix=subset,
                equality=val_equal, index_identical=idx_equal, value_identical=val_equal,
                n_holes=int(len(prefix) - len(v)), m=int(m),
                max_abs_energy_gap=float(np.max(np.abs(a - b))) if val_equal else None)


# ============================================================ the harness: BOTH arms
def run_both(cand: Candidates, E: np.ndarray, alpha: float = 0.15, T: float = 0.5,
             m_fixed: int = M_PROD, seed: int = 0, n_random_draws: int = 16,
             layers: int = 3, iters: int = 80, label_perm: bool = False,
             with_projection: bool = False) -> Dict:
    """ONE CELL: one pool, one energy, the VQE arm and the classical control side by side.

    Returns a flat row. Both arms consume the same `cand` and the same `E` by construction --
    that is the whole point of this function existing. Every RMSD is POINT-CLOUD and ORACLE
    (an oracle number is never a system result; it is the evaluation of an ACHIEVABLE
    selection, and the selection itself is native-free).

    The three classical columns and what each is for:
      `c_matched`  top-m at the VQE's REALISED m -- the size-matched bar (s23's convention).
                   By the theorem this is the same SET as the VQE's, so any difference between
                   `q_rmsd` and `c_matched_rmsd` is numerical, and a non-zero value is a BUG.
      `c_fixed`    top-`m_fixed` (=75, the production rung) -- what the pipeline actually ships.
      `random`     the zero-information null at the VQE's realised m.
    The VQE's entire selection contribution is therefore `c_matched_rmsd - c_fixed_rmsd`: a
    movement along the classical m-ladder, and nothing else is available to it.
    """
    E = np.asarray(E, float)
    q = arm_vqe(cand, E, alpha=alpha, T=T, layers=layers, iters=iters, seed=seed,
                label_perm=label_perm)
    m_q = int(q["m"])
    if m_q == 0:
        raise RuntimeError(f"{cand.pdb}: VQE tail held no real candidate (padding leak)")
    gate = gate_set_equality(E, q["cands"], m_q)
    c_m = arm_classical(cand, E, m_q)
    c_f = arm_classical(cand, E, m_fixed)
    rnd = arm_random(cand, m_q, n_draws=n_random_draws, seed=seed)

    Cq, bq = readout_uniform(cand, q["cands"])
    Cm, _ = readout_uniform(cand, c_m["cands"])
    Cf, _ = readout_uniform(cand, c_f["cands"])

    row = dict(
        pdb=cand.pdb, n=cand.n, fold=cand.fold, k=cand.k, source=cand.source,
        basis="point_cloud", label="ACHIEVABLE", seed=int(seed),
        alpha=float(alpha), T=float(T), m_fixed=int(m_fixed),
        # --- quantum arm
        q_m=m_q, q_tail_size=int(q["tail_size"]), q_pad_in_tail=int(q["n_pad_in_tail"]),
        q_cvar=q["cvar"], q_entropy_bits=q["entropy_bits"], q_ess=q["ess"],
        q_n_qubits=q["n_qubits"], q_layers=q["layers"], q_iters=q["iters"],
        q_gauge_permuted=q["gauge_permuted"], q_secs=q["secs"],
        q_rmsd=_oracle_rmsd(cand, Cq),
        # --- classical controls on IDENTICAL candidates and IDENTICAL energies
        c_matched_m=int(c_m["m"]), c_matched_rmsd=_oracle_rmsd(cand, Cm),
        c_fixed_m=int(c_f["m"]), c_fixed_rmsd=_oracle_rmsd(cand, Cf),
        # --- zero-information null, matched in the operator's space
        rand_rmsd_mean=rnd["oracle_rmsd_mean"], rand_rmsd_sd=rnd["oracle_rmsd_sd"],
        rand_rmsd_best=rnd["oracle_rmsd_best"], rand_rmsd_worst=rnd["oracle_rmsd_worst"],
        rand_n_draws=int(rnd["n_draws"]),
        # --- soundness
        gate_pass=bool(gate["pass_"]), gate_equality=bool(gate["equality"]),
        gate_index_identical=bool(gate["index_identical"]), gate_n_holes=int(gate["n_holes"]),
        # --- pool descriptors, native-free
        E_mean=float(E.mean()), E_sd=float(E.std()), E_min=float(E.min()),
        medoid_of_q=int(bq),
    )
    # provenance of the retained set when the pool is a merge -- L2(d)'s 0.355 statistic
    split = cand.meta.get("split")
    if split is not None:
        s = np.asarray(split, int)
        row["q_share_second_source"] = float(s[q["cands"]].mean())
        row["c_fixed_share_second_source"] = float(s[c_f["cands"]].mean())
    # ORACLE pool descriptors -- LABEL ONLY, named so, never an input
    if cand.oracle_rr is not None:
        row["oracle_pool_best"] = float(np.min(cand.oracle_rr))
        row["oracle_pool_mean"] = float(np.mean(cand.oracle_rr))
    if with_projection:
        row["proj_q_rmsd"] = _oracle_rmsd(cand, readout_projected(cand, Cq))
        row["proj_c_fixed_rmsd"] = _oracle_rmsd(cand, readout_projected(cand, Cf))
        row["proj_basis"] = "built_chain"
    return row


#: the FULL key set a completed `run_both` row must carry. `write` checks every one of these
#: on every row -- a row count is not a completion flag.
ROW_KEYS = ("pdb", "n", "fold", "k", "source", "basis", "label", "seed", "alpha", "T",
            "q_m", "q_rmsd", "q_cvar", "q_entropy_bits", "q_n_qubits",
            "c_matched_m", "c_matched_rmsd", "c_fixed_m", "c_fixed_rmsd",
            "rand_rmsd_mean", "rand_rmsd_best", "rand_rmsd_worst",
            "gate_pass", "E_mean", "E_sd")


# =============================================================================== stats
def paired_stats(a, b, folds=None, n_boot=4000, seed=0, name_a="a", name_b="b") -> Dict:
    """Paired contrast with everything BRIEF s4 demands. NEGATIVE means `a` is BETTER.

    SE and the effect as a multiple of its OWN MDE (2.8016*SE, per comparison -- never a
    quoted constant), iid bootstrap CI BESIDE a fold-clustered one, W/L, median beside the
    mean (the free early-warning on concentration), and the worst-target degradation.
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    d = a - b
    nn = len(d)
    if nn < 3:
        return dict(n=int(nn), mean=float(d.mean()) if nn else float("nan"),
                    verdict="NOT MEASURED")
    rng = np.random.default_rng(int(seed))
    se = float(d.std(ddof=1) / math.sqrt(nn))
    mde = 2.8016 * se
    bs = np.array([d[rng.integers(0, nn, nn)].mean() for _ in range(int(n_boot))])
    out = dict(name_a=name_a, name_b=name_b, n=int(nn),
               mean_a=float(a.mean()), mean_b=float(b.mean()),
               mean=float(d.mean()), median=float(np.median(d)), se=se, mde=float(mde),
               eff_over_mde=float(abs(d.mean()) / mde) if mde > 0 else float("inf"),
               ci_iid=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
               W=int((d < 0).sum()), L=int((d > 0).sum()), ties=int((d == 0).sum()),
               worst_degradation=float(d.max()), best_improvement=float(d.min()))
    if folds is not None:
        f = np.asarray(folds)[ok]
        uf = np.unique(f)
        per = {int(u): float(d[f == u].mean()) for u in uf}
        # cluster bootstrap: resample FOLDS, not targets
        cb = []
        for _ in range(int(n_boot)):
            pick = rng.integers(0, len(uf), len(uf))
            cb.append(np.mean(np.concatenate([d[f == uf[i]] for i in pick])))
        out["per_fold"] = per
        out["ci_fold"] = [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))]
        s = np.sign(d.mean())
        out["folds_same_sign"] = int(sum(1 for v in per.values() if np.sign(v) == s))
        out["n_folds"] = int(len(uf))
    lo, hi = out.get("ci_fold", out["ci_iid"])
    r = out["eff_over_mde"]
    # DEFECT FIXED 2026-09-08, found by D2's own sign-flip control. The previous rule was
    #   NULL if CI includes 0, else TYPE-M if 0.7<=r<=1.3, else MEASURED
    # which labelled an effect at 0.39x its own MDE "MEASURED" purely because a 5-fold cluster
    # bootstrap CI happened to exclude zero. With only 5 clusters that CI is unstable, and an
    # effect below its own MDE is by definition one the design could not reliably detect --
    # calling it MEASURED is over-claiming in exactly the direction a lane wants to be flattered.
    # Anything under 1.0x MDE is now named as underpowered rather than promoted.
    out["verdict"] = ("NULL" if lo * hi <= 0 else
                      "TYPE-M ZONE (0.7-1.3x MDE, NOT A RESULT)" if 0.7 <= r <= 1.3 else
                      "UNDERPOWERED (<0.7x MDE, NOT A RESULT)" if r < 0.7 else
                      "MEASURED")
    return out


def aggregate(rows: List[Dict], seed: int = 0) -> Dict:
    """Roll up `run_both` rows. REFUSES to mix bases -- the 0.156 A operator gap is real."""
    if not rows:
        raise ValueError("no rows")
    bases = {r["basis"] for r in rows}
    if len(bases) > 1:
        raise ValueError(f"refusing to aggregate mixed bases {bases}: point-cloud and "
                         f"built-chain RMSD are never compared")
    folds = np.array([r["fold"] for r in rows])
    qr = np.array([r["q_rmsd"] for r in rows], float)
    cm = np.array([r["c_matched_rmsd"] for r in rows], float)
    cf = np.array([r["c_fixed_rmsd"] for r in rows], float)
    rn = np.array([r["rand_rmsd_mean"] for r in rows], float)
    rb = np.array([r["rand_rmsd_best"] for r in rows], float)
    gates = np.array([bool(r["gate_pass"]) for r in rows])
    return dict(
        n_rows=len(rows), basis=rows[0]["basis"], sources=sorted({r["source"] for r in rows}),
        gate_pass_rate=float(gates.mean()), gate_failures=[r["pdb"] for r in rows
                                                           if not r["gate_pass"]],
        # equality is the empirical (full-support) regime, reported apart from the theorem
        gate_equality_rate=float(np.mean([bool(r.get("gate_equality", True)) for r in rows])),
        mean_q=float(np.nanmean(qr)), mean_c_matched=float(np.nanmean(cm)),
        mean_c_fixed=float(np.nanmean(cf)), mean_random=float(np.nanmean(rn)),
        mean_random_bestofK=float(np.nanmean(rb)),
        mean_q_m=float(np.mean([r["q_m"] for r in rows])),
        mean_entropy_bits=float(np.mean([r["q_entropy_bits"] for r in rows])),
        # THE contrast the pillar requires: quantum selection vs the classical control
        vqe_vs_matched=paired_stats(qr, cm, folds, seed=seed,
                                    name_a="vqe", name_b="classical_top_m_matched"),
        vqe_vs_fixed=paired_stats(qr, cf, folds, seed=seed,
                                  name_a="vqe", name_b="classical_top75_PRODUCTION"),
        # the m-ladder movement, which the theorem says is the ONLY channel the circuit has
        matched_vs_fixed=paired_stats(cm, cf, folds, seed=seed,
                                      name_a="classical_at_vqe_m", name_b="classical_top75"),
        # against the zero-information control, matched in the operator's space
        vqe_vs_random=paired_stats(qr, rn, folds, seed=seed,
                                   name_a="vqe", name_b="random_m_subset_MEAN"),
        vqe_vs_random_bestofK=paired_stats(qr, rb, folds, seed=seed,
                                           name_a="vqe",
                                           name_b="random_m_subset_BEST_OF_K_max_dist"),
    )


# =============================================================================== io
def write(name: str, obj: Dict, required_keys: Optional[Sequence[str]] = ROW_KEYS) -> str:
    """Atomic write (tmp + os.replace) + a `_COMPLETE` sidecar on the FULL KEY SET."""
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, path)
            break
        except PermissionError:
            time.sleep(0.4)
    if required_keys is not None:
        rows = obj.get("rows", [])
        missing = sorted({k for r in rows for k in required_keys if k not in r})
        ok = bool(rows) and not missing
        with open(path.replace(".json", "") + "_COMPLETE", "w") as fh:
            json.dump(dict(complete=ok, n_rows=len(rows), missing_keys=missing,
                           required_keys=list(required_keys),
                           time=time.strftime("%Y-%m-%d %H:%M:%S")), fh)
    return path


# ============================================================================ locking
class AmberLock:
    """Exclusive `s24/results/LOCK_AMBER`. The OpenMM lane is serialised and it is Lane D's.

    `os.open(..., O_CREAT|O_EXCL)` -- an atomic create that FAILS if the file exists. Holds
    for bounded periods only; the context manager always releases, including on exception.
    """

    def __init__(self, tag: str = "", wait: float = 15.0, tries: int = 240):
        self.path = os.path.join(RESULTS, "LOCK_AMBER")
        self.tag, self.wait, self.tries, self.fd = tag, float(wait), int(tries), None

    def __enter__(self):
        for _ in range(self.tries):
            try:
                self.fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(self.fd, json.dumps(
                    dict(pid=os.getpid(), tag=self.tag, lane="D",
                         opened=time.strftime("%Y-%m-%d %H:%M:%S"))).encode())
                os.close(self.fd)
                self.fd = 1
                print(f"[LOCK_AMBER] OPENED by lane D pid={os.getpid()} tag={self.tag}",
                      flush=True)
                return self
            except FileExistsError:
                time.sleep(self.wait)
        raise TimeoutError(f"LOCK_AMBER held by another process after "
                           f"{self.tries * self.wait:.0f}s")

    def __exit__(self, *exc):
        if self.fd is not None:
            try:
                os.remove(self.path)
            except OSError:
                pass
            print(f"[LOCK_AMBER] RELEASED by lane D pid={os.getpid()} tag={self.tag}",
                  flush=True)
        return False
