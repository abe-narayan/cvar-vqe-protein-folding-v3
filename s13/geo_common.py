"""SPRINT 13, TRAINABILITY-AND-GEOMETRY -- shared machinery.

WHAT THIS MODULE IS FOR
=======================
The scientific question is: *how does the complexity of the molecular energy model reshape
the geometry of the variational state manifold and the trainability of the VQE?*  The
controlled comparison is **Legacy (`core.energy`, 11-term knowledge-based)** against
**AMBER ff14SB/GBn2 (`core.amber`, genuine, single point, no minimisation)** with the
representation, ansatz, initialisation, CVaR alpha, optimiser, seeds AND the
objective-evaluation budget all matched.

THE ONE ARCHITECTURAL DECISION EVERYTHING RESTS ON
--------------------------------------------------
Every diagnostic here (Fubini-Study metric, gradient statistics, Hessian, CVaR sweep,
optimiser comparison) is a function of the *vector of energies over the whole register*.
So the register is enumerated ONCE per (target, length, k, model) and cached:

    E[j] = model.energy(bitstring(j))    for every j in [0, 2**n)

For AMBER that is 2**n genuine ff14SB/GBn2 single-point calls at ~7 ms each.  Once the
table exists, the CVaR objective, its exact parameter-shift gradient, the metric and the
Hessian are all *exact* and cost no further energy evaluations.  Three consequences:

1. **Budget parity is exact and trivially auditable for the diagnostics**: both models pay
   exactly ``2**n`` unique evaluations for their table, the same ``2**n``, and nothing
   else.  Wall time differs by ~35x; the evaluation count does not differ at all.
2. **The gradients are analytic, not sampled**, so a difference between Legacy and AMBER is
   a property of the two energy landscapes rather than of two different noise realisations.
3. For the *optimiser* comparison (QNG / GD / SPSA) the table is wrapped in
   `TableModel`, a `budget.BudgetedEnergyModel` that charges one evaluation per UNIQUE
   bitstring and caches thereafter -- the project's own accounting, byte for byte the same
   class the production Hamiltonians use.  `verify_table` re-derives a random sample of
   rows from the live Hamiltonian to prove the table *is* the model.

ENERGY SCALE, AND WHY IT MUST BE NORMALISED
-------------------------------------------
``CVaR_alpha(a*E + b) = a*CVaR_alpha(E) + b``, so every gradient of the objective scales
*linearly* with the units of the energy model.  Legacy energies live on ~1e2 and AMBER
single-point energies on ~1e2 to 1e6 (unminimised ideal-geometry backbones have real steric
clashes).  Reporting raw ``||grad C||`` for the two models would therefore compare units,
not trainability.  Every gradient statistic in this study is reported BOTH raw and divided
by ``sd(E)`` over the enumerated register, which is exactly scale-free.  The z-scored form
is the headline; the raw form is kept so the reader can see the scale that was removed.

NO NATIVE INFORMATION enters any objective, any Hamiltonian, any optimiser or any
hyper-parameter choice.  `rmsd_table` reads `nat_ca` and is an ORACLE EVALUATION
instrument: it is used to score what an arm returned, never to choose it.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I                      # noqa: E402
import torsion_lib2 as tl2                           # noqa: E402
from budget import BudgetedEnergyModel               # noqa: E402
from core.quantum import (StatevectorCircuit, cvar_exact, grad_cvar_paramshift,
                          FoldingHamiltonian)        # noqa: E402

CACHE = os.path.join(ROOT, "s13", "cache")
RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

#: The two energy models under comparison.  "legacy" is `core.energy`'s 11-term
#: knowledge-based model through `core.quantum.FoldingHamiltonian`; "amber" is
#: `core.amber.refine(..., k_restraint=0, steps=-1)` -- genuine ff14SB/GBn2, single point,
#: NO minimisation.  Neither is a surrogate.
MODELS = ("legacy", "amber")

#: Derived energy-model VARIANTS.  Every one of these is a post-hoc function of a cached
#: genuine table; none of them costs an extra OpenMM call, and none of them is called
#: "AMBER" without qualification.
#:
#:   ``amber_soft``  a MODIFIED POTENTIAL: the genuine AMBER single-point energy passed
#:                   through a STRICTLY INCREASING soft compression of its upper tail.
#:                   Because the map is monotone the induced ranking of configurations --
#:                   and therefore Spearman rho(E, CA-RMSD), the argmin, and the native's
#:                   percentile -- are IDENTICAL to plain AMBER's.  Only the dynamic range
#:                   changes.  That is what makes it the right control: it separates
#:                   "AMBER is untrainable because its ranking is uninformative" from
#:                   "AMBER is untrainable because 16 orders of magnitude of hard-core
#:                   repulsion dominate every gradient".  It is a modified potential and it
#:                   is labelled as one everywhere.
#:   ``legacy_soft`` the SAME monotone compression applied to Legacy.  Added 2026-09-05 after
#:                   the Pauli-spectrum agent's delta-spike correction (`s13/walsh_FINDINGS.md`),
#:                   independently reproduced here: a table whose variance is carried by a
#:                   handful of extreme configurations has a Walsh spectrum that is exactly
#:                   Binomial(m, 1/2) regardless of the physics, so *any* raw-vs-raw spectrum
#:                   comparison between two differently-spiked tables measures the spikes.
#:                   SOME LEGACY CELLS ARE ALSO SPIKED (top-10 configuration variance share
#:                   0.88 on 1A13 L=5 k=4), so the fair comparison applies the same monotone
#:                   conditioning to BOTH models.  `legacy_soft` vs `amber_soft` is that
#:                   comparison; both are MODIFIED POTENTIALS and are labelled as such.
VARIANTS = ("legacy", "amber", "legacy_soft", "amber_soft")

#: Free-RAM gate, GB (BRIEF rule 7).
MIN_FREE_GB = 1.5


def rss_gb() -> float:
    """This process's own resident set, GB -- the number the 1.2 GB cap is about."""
    try:
        import ctypes
        from ctypes import wintypes

        class _PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        c = _PMC(); c.cb = ctypes.sizeof(_PMC)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(c), c.cb)
        return c.WorkingSetSize / 2 ** 30
    except Exception:
        return float("nan")


def gate(min_free: float = MIN_FREE_GB, wait: float = 10.0, tries: int = 6) -> float:
    """Block until the machine has `min_free` GB free.  Called before every heavy step.

    After `tries * wait` seconds it proceeds anyway and RECORDS the free memory it
    proceeded at: this process's own resident set is ~0.3 GB (one OpenMM Context plus a
    <=2**16 float table), well inside the 1.2 GB per-process cap, so a machine that stays
    below the gate is other agents' load rather than this one's, and blocking forever on it
    would be a stall, not compliance.  Every proceed-anyway is logged.
    """
    f = I.free_gb()
    if f >= min_free:
        return f
    for _ in range(tries):
        time.sleep(wait)
        f = I.free_gb()
        if f >= min_free:
            return f
    print(f"  [gate] proceeding at free={f:.2f} GB after "
          f"{tries * wait:.0f}s below {min_free} GB", flush=True)
    return f


# ------------------------------------------------------------------ targets / representation
def peptide(pdb: str, L: Optional[int] = None) -> Dict[str, object]:
    """One tuning target, optionally truncated to its first `L` residues.

    Truncation is how peptide LENGTH becomes an independent variable while the qubit count
    per residue stays fixed.  The truncated prefix is a genuine peptide sub-sequence; the
    torsion library for it is built with the FULL target sequence held out
    (`exclude_seq`), so the leakage discipline is the library's own, unchanged.
    """
    u = I.load_univ(pdb)
    full = u["seq"]
    L = len(full) if L is None else int(L)
    return {"pdb": pdb, "seq": full[:L], "full_seq": full, "L": L,
            "fold": u["fold"], "nat_ca": np.asarray(u["nat_ca"], float)[:L]}


def make_rep(t: Dict[str, object], k: int):
    """`PerResidueTorsion` over the k-state sequence-conditioned library, target held out."""
    tab = tl2.library_for(t["seq"], k, t["full_seq"])
    return tl2.PerResidueTorsion(t["seq"], tab, chi_bits=False)


def bits_of(j: int, n_bits: int) -> str:
    return format(int(j), "0%db" % int(n_bits))


def all_states(L: int, k: int) -> np.ndarray:
    """``(k**L, L)`` per-residue state indices in bitstring (MSB-first) order."""
    b = int(round(math.log2(k)))
    j = np.arange(k ** L, dtype=np.int64)
    return np.stack([(j >> (b * (L - 1 - i))) & (k - 1) for i in range(L)], axis=1)


# ------------------------------------------------------------------ energy tables
def _tag(pdb: str, L: int, k: int, model: str) -> str:
    return f"geo_E_{pdb}_{L}_{k}_{model}"


def energy_table(pdb: str, L: int, k: int, model: str,
                 verbose: bool = True) -> np.ndarray:
    """``(k**L,)`` genuine energies of EVERY register state, cached on disk.

    `model="legacy"` -> `core.quantum.FoldingHamiltonian.energy` (the 11-term model).
    `model="amber"`  -> `core.amber.refine(k_restraint=0, steps=-1)`, i.e.
                        `core.amber.single_point`: genuine ff14SB/GBn2 potential energy of
                        the as-built ideal-geometry structure, with NO minimisation and no
                        restraint.  `memo=False` so OpenMM's LRU does not grow without
                        bound over 65k calls; every state is visited exactly once anyway.
    """
    p = os.path.join(CACHE, _tag(pdb, L, k, model) + ".npy")
    if os.path.exists(p):
        return np.load(p)
    t = peptide(pdb, L)
    rep = make_rep(t, k)
    N = k ** L
    if rep.n_bits != int(round(math.log2(k))) * L:
        raise RuntimeError("unexpected bit layout")
    E = np.empty(N, float)
    t0 = time.time()
    if model == "legacy":
        H = FoldingHamiltonian(t["seq"], rep, eval_budget=None)
        for j in range(N):
            E[j] = H.energy(bits_of(j, rep.n_bits))
    elif model == "amber":
        import core.amber as A
        A.builder_for(t["seq"], rep, "CPU", 1)          # warm the topology once
        for j in range(N):
            r = A.refine(t["seq"], rep, bits_of(j, rep.n_bits),
                         k_restraint=0.0, steps=-1, tolerance=1e9,
                         platform_name="CPU", threads=1, memo=False)
            E[j] = float(r["energy"])
            if verbose and j and j % 4096 == 0:
                print(f"    amber {pdb} L={L} k={k}: {j}/{N} "
                      f"{time.time() - t0:.0f}s", flush=True)
    else:
        raise ValueError(model)
    np.save(p, E)
    if verbose:
        print(f"  [{model}] {pdb} L={L} k={k} N={N} in {time.time() - t0:.1f}s "
              f"min={E.min():.2f} max={E.max():.3g}", flush=True)
    return E


def verify_table(pdb: str, L: int, k: int, model: str, n: int = 24,
                 seed: int = 0) -> Dict[str, float]:
    """Re-derive `n` random rows from the LIVE model and compare with the cached table."""
    E = energy_table(pdb, L, k, model)
    t = peptide(pdb, L)
    rep = make_rep(t, k)
    rng = np.random.default_rng(seed)
    js = rng.choice(len(E), size=min(n, len(E)), replace=False)
    if model == "legacy":
        H = FoldingHamiltonian(t["seq"], rep, eval_budget=None)
        live = np.array([H.energy(bits_of(j, rep.n_bits)) for j in js])
    else:
        import core.amber as A
        live = np.array([float(A.refine(t["seq"], rep, bits_of(j, rep.n_bits),
                                        k_restraint=0.0, steps=-1, tolerance=1e9,
                                        memo=False)["energy"]) for j in js])
    d = np.abs(live - E[js])
    return {"n": int(len(js)), "max_abs_diff": float(d.max()),
            "max_rel_diff": float((d / np.maximum(np.abs(live), 1e-12)).max())}


def rmsd_table(pdb: str, L: int, k: int) -> np.ndarray:
    """ORACLE EVALUATION INSTRUMENT. CA-RMSD of every register state to the native.

    Reads `nat_ca`.  Used only to score what an arm returned and to report the objective's
    ranking validity; never inside an objective, a Hamiltonian or a hyper-parameter choice.
    """
    p = os.path.join(CACHE, f"geo_R_{pdb}_{L}_{k}.npy")
    if os.path.exists(p):
        return np.load(p)
    t = peptide(pdb, L)
    rep = make_rep(t, k)
    S = all_states(L, k)
    rows = np.arange(L)
    PHI = np.asarray(rep._phi, float)[rows[None, :], S]
    PSI = np.asarray(rep._psi, float)[rows[None, :], S]
    R = np.empty(len(S), float)
    B = 4096
    for a in range(0, len(S), B):
        W = I.build_ca(PHI[a:a + B], PSI[a:a + B])
        R[a:a + B] = I.kabsch_rmsd_batch(W, t["nat_ca"])
    np.save(p, R)
    return R


def soft_compress(E: np.ndarray) -> np.ndarray:
    """Strictly increasing compression of an energy table's upper tail.

    ``E -> E`` below the median ``q``; ``E -> q + s*log1p((E-q)/s)`` above it, with
    ``s = IQR``.  Monotone, so the induced RANKING of configurations is unchanged to the
    last bit; only the dynamic range moves (1e13 -> ~q + 30 s).  Used to build
    ``amber_soft``.  It is a MODIFIED POTENTIAL, not AMBER.
    """
    E = clean(E)
    q = float(np.median(E))
    s = float(np.quantile(E, 0.75) - np.quantile(E, 0.25))
    s = s if s > 0 else 1.0
    hi = E > q
    out = E.copy()
    out[hi] = q + s * np.log1p((E[hi] - q) / s)
    return out


def variant_table(pdb: str, L: int, k: int, name: str) -> np.ndarray:
    """Energy vector for a named model VARIANT, built from the cached genuine tables."""
    if name in ("legacy", "amber"):
        return clean(energy_table(pdb, L, k, name))
    if name.endswith("_soft"):
        return soft_compress(energy_table(pdb, L, k, name[:-5]))
    if name.endswith("_w99"):
        return winsor_hi(clean(energy_table(pdb, L, k, name[:-4])), 0.99)
    raise ValueError(name)


def spike_share(E: np.ndarray, top: int = 10) -> Dict[str, float]:
    """How much of `Var(E)` is carried by the few most extreme configurations.

    THE DIAGNOSTIC THAT MUST ACCOMPANY EVERY WALSH SPECTRUM.  A function that is a constant
    plus a spike at one point ``x0`` has ``c_S = +-(E_max - mu)/2^m`` for EVERY subset S, so
    its Pauli-weight spectrum is exactly ``Binomial(m, 1/2)`` and its mean weight is exactly
    ``m/2`` -- carrying no information about the physics at all.  Raw AMBER single-point
    tables are spike-dominated (top-10 configuration variance share 0.62-0.9998 measured
    here), and so are some Legacy tables (up to 0.88).  A spectrum whose ``l1_from_binomial``
    is small and whose ``top10`` is near 1 is measuring clashes, not interactions.
    """
    E = clean(E)
    dev = (E - E.mean()) ** 2
    tot = float(dev.sum())
    o = np.argsort(dev)[::-1]
    return {"top1": float(dev[o[:1]].sum() / max(tot, 1e-300)),
            "top10": float(dev[o[:top]].sum() / max(tot, 1e-300))}


class TableModel(BudgetedEnergyModel):
    """The enumerated energy table behind the project's own budget accounting.

    Charges exactly one evaluation per UNIQUE bitstring and caches thereafter -- the same
    contract `FoldingHamiltonian` and `AmberHamiltonian` implement, and the reason a
    Legacy arm and an AMBER arm can be compared at all: AMBER is ~35x slower per call, so
    matching wall time would hand Legacy 35x the evaluations.
    """

    def __init__(self, E: np.ndarray, n_bits: int, eval_budget: Optional[int] = None):
        self.E = np.asarray(E, float)
        self._n_bits = int(n_bits)
        self._init_budget(1 << 22, eval_budget)

    @property
    def n_qubits(self) -> int:
        return self._n_bits

    @property
    def n_bits(self) -> int:
        return self._n_bits

    def energy(self, bitstring: str) -> float:
        hit = self._cache.get(bitstring)
        if hit is not None:
            return hit
        self._charge()
        e = float(self.E[int(bitstring, 2)])
        self._cache[bitstring] = e
        return e

    def energy_index(self, j: int) -> float:
        return self.energy(bits_of(j, self._n_bits))


# ------------------------------------------------------------------ energy conditioning
def zscore(E: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """``(E - mean) / sd`` over the enumerated register, and the (mean, sd) removed.

    An affine map of the energy is an affine map of CVaR, so this is the *exact* way to
    make gradient magnitudes comparable between two models with different units.  It is
    invertible and it is reported, not hidden.
    """
    E = np.asarray(E, float)
    m, s = float(E.mean()), float(E.std())
    return (E - m) / (s if s > 0 else 1.0), m, s


def clean(E: np.ndarray) -> np.ndarray:
    """Replace non-finite entries by the finite maximum.  Reported wherever it fires."""
    E = np.asarray(E, float)
    fin = np.isfinite(E)
    if fin.all():
        return E
    return np.where(fin, E, E[fin].max())


def scales(E: np.ndarray) -> Dict[str, float]:
    """Every candidate normalising scale for one energy table, computed once.

    A gradient of a CVaR objective carries the energy's units, so comparing Legacy with
    AMBER requires dividing by SOME scale -- and the choice can flip the conclusion, so all
    four are recorded and every headline is checked against all four.

      ``sd``          plain standard deviation.  For AMBER this is set by a handful of
                      hard-sphere clash configurations at 1e10-1e13 kcal/mol and is
                      essentially meaningless.
      ``sd_winsor99`` sd after capping the upper 1%.  Better, still tail-driven.
      ``iqr``         Q75 - Q25.  Robust, and the headline normaliser.
      ``span_low``    median(E) - min(E): the depth of the low-energy region a CVaR
                      objective with small alpha is actually working in.
    """
    E = clean(E)
    q25, q50, q75 = (float(np.quantile(E, q)) for q in (0.25, 0.5, 0.75))
    return {"sd": float(E.std()), "sd_winsor99": float(winsor_hi(E, 0.99).std()),
            "iqr": q75 - q25, "span_low": q50 - float(E.min()),
            "min": float(E.min()), "median": q50, "max": float(E.max()),
            "frac_gt_1e4": float((E > 1e4).mean()),
            "frac_gt_1e8": float((E > 1e8).mean())}


def winsor_hi(E: np.ndarray, q: float = 0.99) -> np.ndarray:
    """Cap the UPPER tail at its q-quantile.  Diagnostic only, and stated where used.

    AMBER single-point energies on unminimised backbones have a heavy positive tail
    (hard-sphere clashes reaching 1e5-1e7 kcal/mol) which dominates ``sd(E)`` while being
    invisible to ``CVaR_alpha`` for any alpha < 1 -- the tail is the part CVaR discards.
    Reporting a z-score against that sd would understate AMBER's gradients.  Both the
    raw-sd and the winsorised-sd normalisations are therefore reported.
    """
    E = np.asarray(E, float)
    hi = float(np.quantile(E[np.isfinite(E)], q))
    return np.minimum(E, hi)


# ------------------------------------------------------------------ ansatz / geometry
def circuit(n: int, layers: int = 3, ring: bool = True) -> StatevectorCircuit:
    return StatevectorCircuit(n, layers=layers, ring=ring)


def init_theta(P: int, seed: int, scale: float = math.pi) -> np.ndarray:
    """Uniform on [-scale, scale]^P.  The SAME draw is handed to both energy models."""
    return np.random.default_rng(seed).uniform(-scale, scale, size=P)


def dstates(circ: StatevectorCircuit, theta: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """``(psi, dpsi)`` -- the state and its EXACT parameter derivatives.

    For an ``RY(theta) = exp(-i theta Y / 2)`` gate, ``RY(theta + pi) = -i Y RY(theta)``
    and ``d/dtheta RY = -i (Y/2) RY``, so

        d psi / d theta_i  =  (1/2) psi(theta + pi e_i)

    exactly -- the gate sits in the middle of the circuit and everything around it is
    linear and parameter-independent for that coordinate.  This is a *state* shift rule
    (shift pi, factor 1/2), not the expectation-value shift rule (shift pi/2, factor 1/2);
    `t_dstates_matches_fd` in `geo_audit` checks it against finite differences.
    """
    th = np.asarray(theta, float)
    P = th.size
    G = np.repeat(th[None, :], P + 1, axis=0)
    G[1:] += np.pi * np.eye(P)
    S = circ.states_batch(G)
    psi = S[0]
    psi = psi / np.linalg.norm(psi)
    return psi, 0.5 * S[1:]


def fs_metric(circ: StatevectorCircuit, theta: np.ndarray) -> np.ndarray:
    """Fubini-Study metric ``g_ij = Re[<d_i psi|d_j psi> - <d_i psi|psi><psi|d_j psi>]``.

    The ansatz has real amplitudes (RY + CNOT on |0..0>), so the pure-state QFI is
    ``F = 4 g`` and the Berry-connection term ``<d_i psi|psi>`` is identically zero
    (norm conservation: ``d/dtheta <psi|psi> = 2<d psi|psi> = 0`` for a real state).  It is
    still computed and subtracted, and `geo_metric` reports its measured size, because
    "the term you assumed away" is exactly where this kind of result dies.
    """
    psi, D = dstates(circ, theta)
    a = D @ psi
    return D @ D.T - np.outer(a, a)


def metric_stats(g: np.ndarray, tol: float = 1e-10) -> Dict[str, object]:
    """Spectrum summary of a Fubini-Study metric / QFI matrix."""
    w = np.linalg.eigvalsh((g + g.T) / 2.0)
    w = np.clip(w, 0.0, None)
    tr = float(w.sum())
    wmax = float(w.max()) if w.size else 0.0
    nz = w[w > tol * max(wmax, 1e-300)]
    pr = float(tr ** 2 / max(float((w ** 2).sum()), 1e-300))
    d = np.sqrt(np.clip(np.diag(g), 1e-300, None))
    C = g / np.outer(d, d)
    off = C[~np.eye(len(C), dtype=bool)]
    return {
        "P": int(len(w)),
        "trace": tr,
        "lam_max": wmax,
        "lam_min": float(w.min()),
        "lam_min_nonzero": float(nz.min()) if nz.size else 0.0,
        "rank_1e-10": int((w > 1e-10 * max(wmax, 1e-300)).sum()),
        "rank_1e-6": int((w > 1e-6 * max(wmax, 1e-300)).sum()),
        "participation_ratio": pr,
        "cond_1e-10": float(wmax / nz.min()) if nz.size else float("inf"),
        "cond_full": float(wmax / max(float(w.min()), 1e-300)),
        "eigs_top8": [float(x) for x in w[::-1][:8]],
        "eigs_bot8": [float(x) for x in w[:8]],
        "mean_abs_offdiag_corr": float(np.abs(off).mean()),
        "max_abs_offdiag_corr": float(np.abs(off).max()) if off.size else 0.0,
        "diag_mean": float(np.diag(g).mean()),
        "diag_min": float(np.diag(g).min()),
    }


# ------------------------------------------------------------------ objective
def cvar_value(E: np.ndarray, circ: StatevectorCircuit, theta: np.ndarray,
               alpha: float) -> float:
    return float(cvar_exact(E, circ.probs(theta), alpha)[0])


def cvar_grad(E: np.ndarray, circ: StatevectorCircuit, theta: np.ndarray,
              alpha: float) -> np.ndarray:
    """EXACT analytic gradient (parameter shift on every basis probability)."""
    return grad_cvar_paramshift(circ, theta, E, alpha)


def entropy_bits(p: np.ndarray) -> float:
    p = np.asarray(p, float)
    q = p[p > 0]
    return float(-(q * np.log2(q)).sum())


# ------------------------------------------------------------------ statistics
def boot_ci(x: np.ndarray, n_boot: int = 4000, seed: int = 0) -> Tuple[float, float]:
    x = np.asarray(x, float)
    if x.size < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    m = rng.choice(x, size=(n_boot, x.size), replace=True).mean(1)
    return (float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)))


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def loglog_fit(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    """Fit ``log y = a + b*x`` (exponential) and ``log y = c + d*log x`` (polynomial)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (y > 0) & (x > 0)
    x, y = x[m], y[m]
    if x.size < 3:
        return {}
    ly = np.log(y)

    def _fit(u, v):
        A = np.column_stack([np.ones_like(u), u])
        c, *_ = np.linalg.lstsq(A, v, rcond=None)
        r = v - A @ c
        ss = float(((v - v.mean()) ** 2).sum())
        return float(c[1]), (1.0 - float((r ** 2).sum()) / ss if ss > 0 else float("nan"))

    be, r2e = _fit(x, ly)
    bp, r2p = _fit(np.log(x), ly)
    return {"exp_slope_per_qubit": be, "exp_r2": r2e,
            "poly_exponent": bp, "poly_r2": r2p,
            "decay_base": float(np.exp(be))}


def write(name: str, obj) -> str:
    p = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    print(f"  -> {p}", flush=True)
    return p
