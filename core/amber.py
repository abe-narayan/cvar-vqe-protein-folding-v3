"""Consolidated AMBER physics: ff14SB + GBn2 through OpenMM.

This is `amber_refine` + `amber_hamiltonian` + `budget` + `sidechains` in one module,
with the driver layer rewritten around a measurement rather than an assumption.

Where the time actually goes
----------------------------
Measured on 1A13 (INWKGIAAMAKKLL, n=14, 237 atoms) under the production recipe
``refine_coords(k_restraint=K_MODERATE, steps=0, tolerance=1.0, components=True)``,
CPU platform, Threads=1:

    refine_coords                       5928 ms / structure
      LocalEnergyMinimizer.minimize     5900 ms      99.5%
      everything else in Python           28 ms       0.5%
        _assemble (frozen H frames)        5.8 ms
        _heavy_positions (sidechains)      0.7 ms
        getState x7, marshalling          19   ms

and inside the minimiser:

    one full energy+force evaluation     7.08 ms
      CustomGBForce  (GBn2 solvation)    7.13 ms     ~98%
      NonbondedForce                     0.109 ms
      PeriodicTorsionForce               0.085 ms
      HarmonicAngleForce                 0.058 ms
      HarmonicBondForce                  0.032 ms
    evaluations to converge at tol=1.0   708-847

So the cost of this module is, to within a rounding error,

    (number of L-BFGS evaluations) x (cost of one GBn2 evaluation)

Two consequences, both of which shaped what this file does and does not do.

**The System/Context rebuild hypothesis is false.** `amber_refine` already cached one
builder per (sequence, representation shape, platform) in `_BUILDERS`, so a pool of
candidates for one target already shared one ForceField, one System and one Context and
only pushed new positions. There was no per-candidate rebuild to remove. The entire
Python/OpenMM boundary -- coordinate marshalling, `Quantity` construction, unit
conversion, every `getState` -- is 0.5% of the call. Optimising it perfectly would buy
1.005x.

**Therefore the hot arithmetic in this file is copied verbatim and deliberately not
"optimised".** `_frame`, `_assemble`, `_place_ring` and the sidechain templates are
byte-for-byte what they were. A vectorised `_frame` is algebraically identical and *not*
bit-identical -- `np.linalg.norm` scales before squaring and `np.dot` accumulates in its
own order, so a batched rewrite disagrees in the last ulp on a third of random 3-vectors
(see `_frame`'s own docstring, where this was already measured). Trading exact
reproducibility for 0.1% is a bad trade, and the pinned golden energies and `budget.py`'s
cache invariant both depend on it not being made.

What this module does instead
-----------------------------
Three changes, in descending order of measured value. All three are exactly equivalent:
they return the same floats as the serial implementation, not merely close ones.

1. `refine_many` -- process-level parallelism across candidates. Candidates are
   independent, so N worker processes each hold one Context at Threads=1 and run the
   identical code path. Bit-identical per structure, near-linear in cores.

2. `_MEMO` -- a content-addressed cache on the exact input coordinates. `_run` is a pure
   function of (builder, heavy positions, k_restraint, steps, tolerance), which is the
   same fact that makes Context reuse valid in the first place. Repeated candidates in a
   retrieval pool cost one hash instead of 5.9 s.

3. `_BUILDERS` is now a bounded LRU with deliberate teardown. It previously grew without
   limit: a 126-target sweep retained 126 live OpenMM Contexts.

A fourth lever exists, is exposed as `threads=N`, and should be treated as a **different
experiment rather than a faster one**. The CPU platform threads CustomGBForce, and on a
single force evaluation the effect is exactly what you would expect from reduction order:
2082.5538003250876 at 1 thread against 2082.55372507416 at 4, a relative shift of 3.6e-8.
That number is misleading, and measuring only it would have been the mistake here.

Pushed through a converged minimisation it does not stay small. Each L-BFGS step is
steered by forces that differ in the last few bits, and ~800 steps descending from
10^4-10^6 kcal/mol of builder strain land in a *different local minimum*. Measured on
1A13, single point against converged:

    single point   relative 6.5e-08   positions identical
    converged      relative 3.0e-04   0.147 kcal/mol, atoms moved 0.065 A

a 4600x amplification, for 1.9-2.7x of wall. 0.147 kcal/mol is small but it is not
nothing against a retrieval pool whose interaction energies span ~24 kcal/mol.

Finding that required fixing something first. `_calibrate_hydrogens` minimises hydrogens
against a scratch System, and it inherited the run's thread count, so at Threads=4 it
landed in a different minimum and *every* structure afterwards inherited different frozen
hydrogen frames -- a 144 kcal/mol shift at the single point, i.e. a different builder
rather than a different rounding. That calibration is now pinned to one thread whatever
the run uses (no effect on the default path, which is Threads=1 already), so `threads=N`
now diverges only through the minimiser. `tests/test_amber.py` asserts both halves.

Use `refine_many` when you want this faster. Use `threads=N` only when you have accepted
that the numbers will not match the pinned set and are not comparing against anything that
was computed at Threads=1.

Hard invariants preserved
-------------------------
Genuine ff14SB + GBn2: `app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")`
with `CustomGBForce` in the built System (asserted in `tests/test_amber.py`). AMBER
parameters, the restraint constant, the tolerance and the `bond+angle > 1000` strain
rejection are untouched. `refine_coords` is the coordinate entry point; `refine` still
projects onto the discrete state library and is still the wrong call for real geometry.

Reference: on 1A13 the native interaction energy (nonbonded + solvation, k=10, steps=0,
tolerance=1.0) is -489.9138948277905 kcal/mol. This module reproduces it bit-for-bit.

THE PLATFORM QUESTION, SETTLED
------------------------------
Everything below is measured, and the raw tables are in `verify/amber_platform.json`
(reproduce with `verify/amber_platform.py`, `amber_affinity.py`, `amber_threads.py`).
This box reports `Reference 1.0 / CPU 10.0 / OpenCL 50.0` and essentially all the cost is
in one `CustomGBForce`, so the GPU looked like the obvious lever. It is not.

**OpenCL: a structural negative.** It initialises (Intel Arc 140V) and it does honour
`Precision=double`. The disqualifying fact is that the OpenCL platform does not expose
`DeterministicForces`, which the CPU platform does and which this module sets to true --
that setting is what makes the golden above reproducible. Bitwise reproducibility is
therefore unavailable on OpenCL by construction, not by a switch we failed to find.
Measured: three runs of the SAME input gave -489.94942 / -489.84657 / -489.83861, a
0.11 kcal/mol spread and 7.2e-3 A of CA movement. A tolerance presupposes that the same
input gives the same answer twice. And the speed is not there either -- six targets gave
3.33x, 0.93x, 0.80x, 1.18x, 0.41x, 1.20x, geometric mean ~1.06x, because these systems are
~230 atoms and kernel launch overhead swamps the GBn2 kernel. Two traps make a naive
OpenCL arm look plausible while being wrong twice over: this module sets no `Precision`
property, so `platform_name="OpenCL"` silently runs SINGLE precision (8.5 kcal/mol off);
and building the Hamiltonian on OpenCL re-solves `_calibrate_hydrogens`, which is itself a
`LocalEnergyMinimizer` call, moving the hydrogen frames 0.188 nm and the energy
2.93 kcal/mol before any production work happens.

**THIS BOX HAS FOUR FAST CORES, NOT EIGHT.** Ranked with the AMBER workload itself rather
than with a proxy kernel -- a GBn2 minimisation is bandwidth-bound and does not rank cores
the way a spin loop does -- the per-core cost of one structure is 4.375 / 4.260 / 4.557 /
4.828 s on cores 0-3 and 6.661 / 6.146 / 7.099 / 6.073 s on cores 4-7. A clean 4/4 split
at 1.67x, which is the Core Ultra 7 256V's P-core / LP-E-core boundary. Summed, the machine
delivers **6.43 fastest-core-equivalents, not 8**, and any feasibility bound of the form
`wall >= CPU_work / 8` is 24% optimistic. It also explains the worker sweep saturating at
~3.9 effective workers without invoking memory at all: that is the P-core count.

**Affinity does not help, and pinning to the fast cores hurts.** One process on the fastest
core is the baseline; four unpinned workers cost +14.4% CPU per structure (reproducing the
sweep's +14.9% in a pure AMBER harness whose footprint is nowhere near the RAM gate, which
rules out page pressure as the cause); four workers PINNED to the four fast cores cost
+24% CPU and 3.3x the wall, reproduced twice. And one process pinned to the SLOWEST core,
with exactly zero contention, costs +40.7% -- larger than the entire w4 regression. So the
w4 inflation is dominated by placement onto slow cores, not by L3 contention, and there is
no cheap bit-identical affinity fix. Every arm returned identical energies, which verifies
rather than assumes that placement cannot change arithmetic.

THE INFLATION FLOOR DEPENDS ON THE SCHEDULE, and an earlier version of this docstring
got that wrong by calling `mean(1/rel)` "the floor". It is the prediction for ONE regime.
With N cores of unequal speed there are two:

    equal COUNT of structures per core   mean(1/rel) = 1.2908   static equal-count chunks
    equal WALL time per core             N/sum(rel)  = 1.2449   any dynamic dispatch

Under dynamic dispatch every core stays busy the same wall time, so the FAST cores complete
more structures and the per-structure CPU average is weighted toward them -- which is why
the dynamic floor is 3.7% below the static one, worth 59 CPU-s on this workload. Both are
now measured against real 8-worker runs and each lands on its own regime's prediction: the
old contiguous equal-count dispatcher gave 1.2913 (+0.03%) and the dynamic longest-first
one gives 1.2533 (+0.67%). The core map is therefore confirmed twice, from two schedules
that disagree with each other, and the residual over the dynamic floor is ~8.6 CPU-s of
1611.59 -- there is no memory-bandwidth contention left in AMBER at all.

A METHODOLOGICAL NOTE, because two agents including this one nearly got it wrong: a
float64 spin loop ranked these cores at 1.24x with no 4/4 pattern, and a cache-resident
matmul at 1.7x with no pattern.  Both led to "heterogeneity is not the constraint".  The
proxy did not add error, it REMOVED THE SIGNAL -- a cache-resident kernel and a
bandwidth-bound one load different parts of the machine and rank the cores by a different
quantity.  Rank hardware with the workload whose cost you are trying to explain.

**Threads: t2 is a regression.** One process on the fast cores: t2 is 0.908x the wall of
t1 while burning 1.35x the CPU-seconds; t4 is 1.969x wall for 1.652x CPU. Since the
feasibility bound counts CPU-seconds, threading makes the bound WORSE and only pays where
cores would otherwise be idle. Deviations from the golden are -3.4e-2 (t2) and -1.6e-2
(t4) kcal/mol -- small, deterministic, but not zero, so `threads=1` remains the reference
path and any threaded arm is a labelled one.  The decisive objection is stronger than the
wall numbers, though: the thread count changes the RESULT, because `_calibrate_hydrogens`
inherits it and settles in a different minimum -- ~144 kcal/mol.  `w8 x t1` is both faster
and bit-identical, so the trade was never live.
"""
from __future__ import annotations

import ctypes
import hashlib
import math
import os
import random
import sys
import time
import warnings
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import openmm
from openmm import app, unit

import protein_geometry as geo

__all__ = [
    # refine API
    "refine", "refine_coords", "refine_ca", "single_point", "refine_many",
    "builder_for", "clear_cache", "cache_stats",
    "K_WEAK", "K_MODERATE", "K_STRONG", "BACKBONE_ATOMS",
    "CONVERGE_MAX_KCAL", "convergence_flags",
    # hamiltonian
    "AmberHamiltonian", "AMBER_TERMS", "KJ_PER_KCAL", "KCAL_PER_KJ",
    # sidechains
    "build_sidechain", "build_full_structure", "sidechain_atom_names",
    "residue_bonds", "heavy_atom_count", "write_full_pdb", "ring_atom_names",
    "SUPPORTED_RESIDUES", "CHI_ANGLES", "CHI1_ROTAMERS", "AROMATIC_RING_RESIDUES",
    "NotImplementedResidueError",
    # budget
    "BudgetExhausted", "BudgetedEnergyModel", "resolve_maxiter",
    "check_optimizer_budget",
    # resources
    "memory_percent", "memory_guard",
]


# ---------------------------------------------------------------- resources
class _MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


def memory_percent() -> float:
    """Physical memory in use, 0-100.

    `GlobalMemoryStatusEx` deliberately, NOT `Get-CimInstance Win32_OperatingSystem`:
    the WMI call costs whole minutes on a loaded box, which makes it useless as a gate
    on the thing that is loading the box. This is a single syscall, microseconds.
    Returns 0.0 off Windows so the guard is a no-op rather than an error.
    """
    if os.name != "nt":
        return 0.0
    s = _MEMORYSTATUSEX()
    s.dwLength = ctypes.sizeof(_MEMORYSTATUSEX)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s))
    return float(s.dwMemoryLoad)


#: Refuse to start more OpenMM workers above this. Sibling processes and an untouchable
#: browser share this box; an OOM here kills someone else's run, not just this one.
MEMORY_LIMIT_PERCENT = 92.0


def memory_guard(limit: float = MEMORY_LIMIT_PERCENT) -> float:
    """Current memory percent, raising `MemoryError` above `limit`."""
    p = memory_percent()
    if p > limit:
        raise MemoryError(
            f"physical memory at {p:.0f}% exceeds the {limit:.0f}% ceiling; "
            "refusing to start more OpenMM contexts")
    return p


# ============================================================================
# SECTION 1 -- evaluation budget: RE-EXPORTED from `budget.py`, not redefined
# ============================================================================
# These four symbols used to be copied into this module, byte for byte.  Two copies of
# an EXCEPTION CLASS is not a style problem, it is a live trap: `except
# budget.BudgetExhausted` does not catch `core.amber.BudgetExhausted`, and the two look
# identical in a diff.  Today nothing is broken, because `core.quantum.FoldingHamiltonian`
# subclasses the root class and the three `except BudgetExhausted` handlers in
# `core.quantum` catch that same root class -- consistent by accident of which copy each
# line happened to reach.  A partial repoint (base class moved, `except` clauses not, or
# the reverse) would send budget exhaustion straight through three handlers that still
# read correctly, turning a graceful stop into a crash inside the VQE path.
#
# THE DEPENDENCY POINTS core -> root, AND THAT DIRECTION IS DELIBERATE.  `budget.py` is on
# the LEGACY arm: it is what `CORE_BACKENDS=legacy` runs.  If it imported from
# `core.amber`, the baseline arm would silently be executing consolidated code, and every
# equivalence proof in this sprint is measured against that arm being pure.  So the
# consolidated module imports from the root module and never the other way round.
#
# `tests/test_amber.py` pins the identity, so a re-introduced copy fails instead of
# passing review.
from budget import (                                                    # noqa: E402
    MAXITER_PER_PARAM,
    MIN_MAXITER_PER_PARAM,
    BudgetExhausted,
    BudgetedEnergyModel,
    check_optimizer_budget,
    resolve_maxiter,
)


# ============================================================================
# SECTION 2 -- all-20-residue sidechain builder (was sidechains.py)
# ============================================================================
class NotImplementedResidueError(NotImplementedError):
    """Raised for residue types this module cannot build."""


THREE = dict(geo.ONE_TO_THREE)
ONE = dict(geo.THREE_TO_ONE)

#: 2026-09-01: the eight remaining standard residues were added. Before that this module
#: built twelve types, which meant energy Model B (Amber ff14SB) could be evaluated on
#: only 6 of the 35 benchmark peptides -- the model was effectively unavailable on 83% of
#: the targets it is supposed to be compared on. The new entries are ideal-geometry NeRF
#: placements at Engh & Huber bond lengths and angles with the modal rotamer as the
#: default chi, exactly the construction the original twelve use.
SUPPORTED_RESIDUES = ("GLY", "ALA", "SER", "THR", "ASP", "GLU",
                      "ASN", "LYS", "PRO", "PHE", "TYR", "TRP",
                      "VAL", "LEU", "ILE", "MET", "CYS", "GLN", "ARG", "HIS")


CHI_ANGLES: Dict[str, Tuple[float, ...]] = {
    "SER": (62.0,),
    "THR": (62.0,),
    "ASP": (-70.0, -15.0),
    "ASN": (-65.0, -20.0),
    "GLU": (-67.0, 180.0, -10.0),
    "LYS": (-67.0, 180.0, 180.0, 180.0),
    "PHE": (-65.0, -85.0),
    "TYR": (-65.0, -85.0),
    "TRP": (-65.0, -90.0),
    "PRO": (),
    "ALA": (),
    "GLY": (),
    "VAL": (175.0,),
    "LEU": (-65.0, 175.0),
    "ILE": (-65.0, 170.0),
    "MET": (-67.0, 180.0, 75.0),
    "CYS": (-65.0,),
    "GLN": (-67.0, 180.0, -25.0),
    "ARG": (-67.0, 180.0, 180.0, 85.0),
    "HIS": (-65.0, -75.0),
}

#: The two dominant chi1 rotamers for aromatic sidechains: m (gauche-, ~50% of observed
#: rotamers) and t (trans, ~35%). One bit per aromatic residue selects between them.
#:
#: Index 0 reproduces `CHI_ANGLES` exactly, so a representation with chi encoding disabled
#: -- or with every chi bit zero -- builds the identical structure it always did. That is
#: what keeps the pinned golden runs reproducible.
CHI1_ROTAMERS: Dict[str, Tuple[float, ...]] = {
    "PHE": (-65.0, 180.0),
    "TYR": (-65.0, 180.0),
    "TRP": (-65.0, 180.0),
}

#: Aromatic residues whose ring this module can build, so real ring geometry (centroid,
#: normal) is available to the energy. Histidine is aromatic but its imidazole template is
#: not implemented, so His falls back to the CB proxy in `energy_terms.aromatic_term`.
AROMATIC_RING_RESIDUES = ("PHE", "TYR", "TRP")


_SPECS: Dict[str, List[Tuple[str, Tuple[str, str, str], float, float, object]]] = {
    "GLY": [],
    "ALA": [],

    "SER": [
        ("OG",  ("N", "CA", "CB"), 1.417, 110.8, ("chi", 0, 0.0)),
    ],


    "THR": [
        ("OG1", ("N", "CA", "CB"), 1.420, 110.1, ("chi", 0, 0.0)),
        ("CG2", ("N", "CA", "CB"), 1.530, 109.2, ("chi", 0, -120.0)),
    ],

    "ASP": [
        ("CG",  ("N", "CA", "CB"),  1.522, 112.6, ("chi", 0, 0.0)),
        ("OD1", ("CA", "CB", "CG"), 1.250, 118.4, ("chi", 1, 0.0)),
        ("OD2", ("CA", "CB", "CG"), 1.250, 118.4, ("chi", 1, 180.0)),
    ],

    "ASN": [
        ("CG",  ("N", "CA", "CB"),  1.521, 112.6, ("chi", 0, 0.0)),
        ("OD1", ("CA", "CB", "CG"), 1.231, 120.8, ("chi", 1, 0.0)),
        ("ND2", ("CA", "CB", "CG"), 1.328, 116.4, ("chi", 1, 180.0)),
    ],

    "GLU": [
        ("CG",  ("N", "CA", "CB"),  1.529, 113.6, ("chi", 0, 0.0)),
        ("CD",  ("CA", "CB", "CG"), 1.523, 112.6, ("chi", 1, 0.0)),
        ("OE1", ("CB", "CG", "CD"), 1.250, 118.4, ("chi", 2, 0.0)),
        ("OE2", ("CB", "CG", "CD"), 1.250, 118.4, ("chi", 2, 180.0)),
    ],

    "LYS": [
        ("CG", ("N",  "CA", "CB"), 1.531, 113.5, ("chi", 0, 0.0)),
        ("CD", ("CA", "CB", "CG"), 1.531, 113.3, ("chi", 1, 0.0)),
        ("CE", ("CB", "CG", "CD"), 1.531, 113.4, ("chi", 2, 0.0)),
        ("NZ", ("CG", "CD", "CE"), 1.486, 110.4, ("chi", 3, 0.0)),
    ],

    "VAL": [
        ("CG1", ("N", "CA", "CB"), 1.521, 110.5, ("chi", 0, 0.0)),
        ("CG2", ("N", "CA", "CB"), 1.521, 110.5, ("chi", 0, 122.0)),
    ],

    "LEU": [
        ("CG",  ("N", "CA", "CB"),  1.530, 116.3, ("chi", 0, 0.0)),
        ("CD1", ("CA", "CB", "CG"), 1.521, 110.7, ("chi", 1, 0.0)),
        ("CD2", ("CA", "CB", "CG"), 1.521, 110.7, ("chi", 1, 122.6)),
    ],

    "ILE": [
        ("CG1", ("N", "CA", "CB"),   1.530, 110.4, ("chi", 0, 0.0)),
        ("CG2", ("N", "CA", "CB"),   1.521, 110.5, ("chi", 0, -122.6)),
        ("CD1", ("CA", "CB", "CG1"), 1.513, 113.8, ("chi", 1, 0.0)),
    ],

    "MET": [
        ("CG", ("N", "CA", "CB"),  1.520, 114.1, ("chi", 0, 0.0)),
        ("SD", ("CA", "CB", "CG"), 1.803, 112.7, ("chi", 1, 0.0)),
        ("CE", ("CB", "CG", "SD"), 1.791, 100.9, ("chi", 2, 0.0)),
    ],

    "CYS": [
        ("SG", ("N", "CA", "CB"), 1.808, 114.4, ("chi", 0, 0.0)),
    ],

    "GLN": [
        ("CG",  ("N", "CA", "CB"),  1.520, 114.1, ("chi", 0, 0.0)),
        ("CD",  ("CA", "CB", "CG"), 1.516, 112.6, ("chi", 1, 0.0)),
        ("OE1", ("CB", "CG", "CD"), 1.231, 120.8, ("chi", 2, 0.0)),
        ("NE2", ("CB", "CG", "CD"), 1.328, 116.4, ("chi", 2, 180.0)),
    ],

    "ARG": [
        ("CG",  ("N",  "CA", "CB"), 1.520, 114.1, ("chi", 0, 0.0)),
        ("CD",  ("CA", "CB", "CG"), 1.520, 111.3, ("chi", 1, 0.0)),
        ("NE",  ("CB", "CG", "CD"), 1.461, 112.0, ("chi", 2, 0.0)),
        ("CZ",  ("CG", "CD", "NE"), 1.329, 124.2, ("chi", 3, 0.0)),
        ("NH1", ("CD", "NE", "CZ"), 1.326, 120.0, 0.0),
        ("NH2", ("CD", "NE", "CZ"), 1.326, 120.0, 180.0),
    ],

    # The imidazole is planar, so the two ring closures are NeRF placements at torsion
    # 180 within the CB-CG-ND1 and CB-CG-CD2 planes rather than a docked template.
    "HIS": [
        ("CG",  ("N",  "CA", "CB"),  1.497, 113.8, ("chi", 0, 0.0)),
        ("ND1", ("CA", "CB", "CG"),  1.378, 122.7, ("chi", 1, 0.0)),
        ("CD2", ("CA", "CB", "CG"),  1.354, 131.0, ("chi", 1, 180.0)),
        ("CE1", ("CB", "CG", "ND1"), 1.321, 109.3, 180.0),
        ("NE2", ("CB", "CG", "CD2"), 1.374, 107.2, 180.0),
    ],
}


_RING_TEMPLATES: Dict[str, Dict[str, Tuple[float, float]]] = {
    "TYR": {
        "CB":  (-2.069719, -0.516542),
        "CG":  (-0.560501, -0.424851),
        "CD1": (0.077557, 0.808854),
        "CD2": (0.222218, -1.572237),
        "CE1": (1.456804, 0.896036),
        "CE2": (1.601875, -1.491804),
        "CZ":  (2.212445, -0.256383),
        "OH":  (3.585912, -0.172939),
    },
    "TRP": {
        "CB":  (-1.827828, -0.599403),
        "CG":  (-0.442130, -0.030330),
        "CD1": (-0.097941, 1.290441),
        "CD2": (0.786942, -0.767086),
        "NE1": (1.270715, 1.421830),
        "CE2": (1.836583, 0.172965),
        "CE3": (1.100358, -2.129541),
        "CZ2": (3.177670, -0.207145),
        "CZ3": (2.430771, -2.503923),
        "CH2": (3.452873, -1.547140),
    },
}

# Phenylalanine is tyrosine's ring without the para hydroxyl. Reusing the same template
# rather than writing new coordinates means it inherits TYR's verified planarity (RMS
# deviation < 1e-6 A) and bond lengths (< 0.02 A of ideal) for free.
_RING_TEMPLATES["PHE"] = {k: v for k, v in _RING_TEMPLATES["TYR"].items() if k != "OH"}

_RING_ANCHOR = {
    #        bond CB-CG, angle CA-CB-CG, bond CG-CD1, angle CB-CG-CD1
    "PHE": (1.502, 113.8, 1.389, 120.8),
    "TYR": (1.512, 113.8, 1.389, 120.8),
    "TRP": (1.498, 113.6, 1.365, 126.9),
}

#: Ring atoms, in template order, for the residues whose rings can be built. These are the
#: atoms `energy_terms` uses to compute a ring centroid and normal, which is what makes a
#: real pi-stacking geometry available instead of a CB-CB proxy.
_RING_ATOMS: Dict[str, Tuple[str, ...]] = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    # The indole's two fused rings are coplanar, so the whole nine-atom system defines one
    # plane and one centroid -- which is also how indole stacking is usually described.
    "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
}


def ring_atom_names(resname: str) -> Tuple[str, ...]:
    """Ring atom names for an aromatic residue, or () if its ring is not implemented."""
    key = str(resname).strip().upper()
    if len(key) == 1:
        key = THREE.get(key, key)
    return _RING_ATOMS.get(key, ())


PRO_RING = {
    "b_CB_CG": 1.526,
    "a_CA_CB_CG": 102.286,
    "t_chi1": 11.675,      # Cg-endo pucker
    "b_CG_CD": 1.526,
    "a_CB_CG_CD": 106.700,
    "t_chi2": -22.549,
}


_SIDECHAIN_BONDS: Dict[str, Tuple[Tuple[str, str], ...]] = {
    "GLY": (),
    "ALA": (),
    "SER": (("CB", "OG"),),
    "THR": (("CB", "OG1"), ("CB", "CG2")),
    "ASP": (("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")),
    "ASN": (("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")),
    "GLU": (("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")),
    "LYS": (("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")),
    "PRO": (("CB", "CG"), ("CG", "CD"), ("CD", "N")),
    "PHE": (("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"),
            ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")),
    "TYR": (("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"),
            ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CZ", "OH")),
    "TRP": (("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "NE1"),
            ("NE1", "CE2"), ("CE2", "CD2"), ("CD2", "CE3"), ("CE3", "CZ3"),
            ("CZ3", "CH2"), ("CH2", "CZ2"), ("CZ2", "CE2")),
    "VAL": (("CB", "CG1"), ("CB", "CG2")),
    "LEU": (("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")),
    "ILE": (("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")),
    "MET": (("CB", "CG"), ("CG", "SD"), ("SD", "CE")),
    "CYS": (("CB", "SG"),),
    "GLN": (("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")),
    "ARG": (("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"),
            ("CZ", "NH1"), ("CZ", "NH2")),
    "HIS": (("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"),
            ("CD2", "NE2"), ("CE1", "NE2")),
}

_ATOM_ORDER: Dict[str, Tuple[str, ...]] = {
    "GLY": ("N", "CA", "C", "O"),
    "ALA": ("N", "CA", "C", "O", "CB"),
    "SER": ("N", "CA", "C", "O", "CB", "OG"),
    "THR": ("N", "CA", "C", "O", "CB", "OG1", "CG2"),
    "PRO": ("N", "CA", "C", "O", "CB", "CG", "CD"),
    "ASP": ("N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"),
    "ASN": ("N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"),
    "GLU": ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"),
    "LYS": ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"),
    "PHE": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2",
            "CZ", "OH"),
    "TRP": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2",
            "CE3", "CZ2", "CZ3", "CH2"),
    "VAL": ("N", "CA", "C", "O", "CB", "CG1", "CG2"),
    "LEU": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"),
    "ILE": ("N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"),
    "MET": ("N", "CA", "C", "O", "CB", "CG", "SD", "CE"),
    "CYS": ("N", "CA", "C", "O", "CB", "SG"),
    "GLN": ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"),
    "ARG": ("N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"),
    "HIS": ("N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"),
}

_ELEMENT = {"N": "N", "O": "O", "C": "C", "S": "S"}


def _resolve(resname: str) -> str:
    """Accept 'W' or 'TRP' (any case) -> three-letter code, or raise."""
    key = str(resname).strip().upper()
    if len(key) == 1:
        if key not in THREE:
            raise NotImplementedResidueError(
                f"sidechain construction not implemented for residue "
                f"{key!r}: unknown one-letter code")
        key = THREE[key]
    if key not in SUPPORTED_RESIDUES:
        raise NotImplementedResidueError(
            f"sidechain construction not implemented for residue {key} "
            f"({ONE.get(key, '?')}); supported residues are "
            f"{', '.join(SUPPORTED_RESIDUES)}")
    return key


def _nerf(a, b, c, length: float, angle_deg: float, torsion_deg: float):
    return np.array(geo._place_atom(a, b, c, float(length),
                                    math.radians(angle_deg),
                                    math.radians(torsion_deg)), dtype=float)


def _torsion_value(spec, chis: Sequence[float]) -> float:
    if isinstance(spec, tuple):
        _, index, offset = spec
        return float(chis[index]) + float(offset)
    return float(spec)


def _frame(origin, x_ref, plane_ref):
    """Orthonormal right-handed frame: e1 towards x_ref, e2 in-plane.

    Deliberately still `np.linalg.norm` / `np.cross` and NOT scalar arithmetic, even
    though profiling puts `np.cross` at ~15% of the energy hot path. A scalar rewrite is
    algebraically the same but not bit-identical: `np.linalg.norm` scales before squaring
    to avoid overflow, and `np.dot` accumulates in its own order, so the two disagree in
    the last ulp on ~10-30% of random 3-vectors. That propagates into ring coordinates and
    then into energies (measured: 133 of 1604 energies moved by up to 1.1e-13), which
    would break the exact reproducibility the golden runs and `budget.py`'s cache
    invariant are pinned to. The speedup here comes from `_ring_dock` instead, which
    removes calls rather than changing arithmetic.
    """
    e1 = np.asarray(x_ref, float) - np.asarray(origin, float)
    e1 = e1 / np.linalg.norm(e1)
    v = np.asarray(plane_ref, float) - np.asarray(origin, float)
    e2 = v - np.dot(v, e1) * e1
    e2 = e2 / np.linalg.norm(e2)
    return e1, e2, np.cross(e1, e2)


#: resname -> (atom names, (k, 3) projection coefficients onto the template frame).
#: The ring template, its frame, and every atom's coordinates in that frame are fixed
#: properties of the residue type, but `_place_ring` recomputed all three on every call --
#: once per aromatic per energy evaluation. Caching them removes one of the two `_frame`
#: calls, the template dict rebuild, and every `np.dot` in the docking loop. The cached
#: coefficients are the same floats the projections produced, so the docking arithmetic
#: below is unchanged.
_RING_DOCK_CACHE: Dict[str, Tuple[Tuple[str, ...], np.ndarray]] = {}


def _ring_dock(resname: str) -> Tuple[Tuple[str, ...], np.ndarray]:
    hit = _RING_DOCK_CACHE.get(resname)
    if hit is not None:
        return hit
    tmpl = _RING_TEMPLATES[resname]
    t3 = {k: np.array([v[0], v[1], 0.0]) for k, v in tmpl.items()}
    te1, te2, te3 = _frame(t3["CG"], t3["CB"], t3["CD1"])
    origin = t3["CG"]
    names = tuple(t3.keys())
    coeffs = np.array([[np.dot(t3[k] - origin, te1),
                        np.dot(t3[k] - origin, te2),
                        np.dot(t3[k] - origin, te3)] for k in names], dtype=float)
    _RING_DOCK_CACHE[resname] = (names, coeffs)
    return names, coeffs


def _place_ring(resname: str, N, CA, CB, chis) -> Dict[str, np.ndarray]:
    """Dock a rigid planar ring template onto the CB-CG-CD1 frame.

    CG and CD1 are placed by NeRF from chi1 and chi2, so the two real degrees
    of freedom are handled exactly as for acyclic sidechains. The template
    then supplies the rigid remainder. Every template atom is coplanar with
    CB, CG and CD1, so the docking has no reflection ambiguity.
    """
    b_cg, a_cg, b_cd1, a_cd1 = _RING_ANCHOR[resname]
    CG = _nerf(N, CA, CB, b_cg, a_cg, chis[0])
    CD1 = _nerf(CA, CB, CG, b_cd1, a_cd1, chis[1])

    names, coeffs = _ring_dock(resname)          # template-only; see `_ring_dock`
    re1, re2, re3 = _frame(CG, CB, CD1)

    out: Dict[str, np.ndarray] = {}
    for name, (c1, c2, c3) in zip(names, coeffs):
        out[name] = CG + c1 * re1 + c2 * re2 + c3 * re3
    return out


def sidechain_atom_names(resname: str) -> Tuple[str, ...]:
    """Sidechain heavy-atom names (CB onwards) in PDB order."""
    key = _resolve(resname)
    return tuple(a for a in _ATOM_ORDER[key] if a not in ("N", "CA", "C", "O"))


def heavy_atom_count(resname: str) -> int:
    """Total heavy atoms for the residue, backbone included (no OXT)."""
    return len(_ATOM_ORDER[_resolve(resname)])


def residue_bonds(resname: str) -> Tuple[Tuple[str, str], ...]:
    """All intra-residue heavy-atom bonds, backbone included."""
    key = _resolve(resname)
    bb = [("N", "CA"), ("CA", "C"), ("C", "O")]
    if key != "GLY":
        bb.append(("CA", "CB"))
    return tuple(bb) + _SIDECHAIN_BONDS[key]


def build_sidechain(resname: str, N, CA, C, CB,
                    chi1: Optional[float] = None) -> Dict[str, np.ndarray]:
    """Heavy-atom sidechain coordinates keyed by PDB atom name.

    Returns CB and everything beyond it. Glycine returns an empty dict: it has
    no CB, and emitting one would break the ff14SB GLY template.

    `C` is accepted for interface completeness; only N, CA and CB determine
    the sidechain frame.

    `chi1` overrides the first chi angle in degrees. `None` uses `CHI_ANGLES`, so every
    existing caller builds exactly the structure it always did. The override is what lets
    `representations` put chi1 in the bitstring: for an aromatic residue the ring's
    orientation about CA-CB is the degree of freedom that decides whether a stacked pair
    is geometrically reachable, and freezing it can put the correct arrangement out of
    reach no matter what the backbone does.
    """
    key = _resolve(resname)
    N = np.asarray(N, dtype=float)
    CA = np.asarray(CA, dtype=float)
    CB = np.asarray(CB, dtype=float)

    if key == "GLY":
        return {}

    atoms: Dict[str, np.ndarray] = {"CB": CB.copy()}

    if key == "PRO":
        # Proline's chi1 is fixed by the ring pucker; an override is meaningless.
        p = PRO_RING
        atoms["CG"] = _nerf(N, CA, CB, p["b_CB_CG"], p["a_CA_CB_CG"],
                            p["t_chi1"])
        atoms["CD"] = _nerf(CA, CB, atoms["CG"], p["b_CG_CD"],
                            p["a_CB_CG_CD"], p["t_chi2"])
        return atoms

    chis = list(CHI_ANGLES[key])
    if chi1 is not None and chis:
        chis[0] = float(chi1)

    if key in _RING_TEMPLATES:
        for name, xyz in _place_ring(key, N, CA, CB, chis).items():
            if name != "CB":
                atoms[name] = xyz
        return atoms

    for name, (an, bn, cn), length, angle, tspec in _SPECS[key]:
        ref = {"N": N, "CA": CA, "CB": CB}
        ref.update(atoms)
        atoms[name] = _nerf(ref[an], ref[bn], ref[cn], length, angle,
                            _torsion_value(tspec, chis))
    return atoms


def build_full_structure(sequence: str, backbone: Dict[str, np.ndarray],
                         add_oxt: bool = True,
                         chi1: Optional[Dict[int, float]] = None
                         ) -> Dict[str, object]:
    """Backbone dict from build_backbone() -> all heavy atoms.

    Returns {"sequence", "residues", "n_atoms"} where "residues" is a list of
    {atom_name: (3,) array}, one per residue, in PDB atom order.

    `chi1` maps residue index -> chi1 in degrees, for residues whose chi1 is encoded.
    Absent entries use `CHI_ANGLES`, so `chi1=None` reproduces the previous behaviour
    exactly -- which is what keeps `run_8state_seed0.py`'s pinned Amber energies valid.
    """
    seq = str(sequence).strip().upper()
    for k in ("N", "CA", "C", "O", "CB"):
        if k not in backbone:
            raise ValueError(f"backbone dict is missing key {k!r}")
    n_res = len(backbone["CA"])
    if len(seq) != n_res:
        raise ValueError(f"sequence length {len(seq)} != {n_res} residues "
                         f"in backbone")
    chi1 = chi1 or {}

    residues: List[Dict[str, np.ndarray]] = []
    for i, aa in enumerate(seq):
        key = _resolve(aa)
        Ni = np.asarray(backbone["N"][i], float)
        CAi = np.asarray(backbone["CA"][i], float)
        Ci = np.asarray(backbone["C"][i], float)
        Oi = np.asarray(backbone["O"][i], float)
        CBi = np.asarray(backbone["CB"][i], float)
        res: Dict[str, np.ndarray] = {"N": Ni, "CA": CAi, "C": Ci, "O": Oi}
        res.update(build_sidechain(key, Ni, CAi, Ci, CBi, chi1=chi1.get(i)))
        ordered = {a: res[a] for a in _ATOM_ORDER[key] if a in res}
        if add_oxt and i == n_res - 1:
            # OXT is the carboxylate oxygen opposite O about the CA-C axis.
            # The ff14SB C-terminal template requires it, and Modeller adds
            # hydrogens but never heavy atoms.
            tor = math.degrees(geo.dihedral(Ni, CAi, Ci, Oi))
            ordered["OXT"] = _nerf(Ni, CAi, Ci, 1.250, 117.0, tor + 180.0)
        residues.append(ordered)
    return {"sequence": seq, "residues": residues,
            "n_atoms": sum(len(r) for r in residues)}


def write_full_pdb(path: str, structure: Dict[str, object],
                   remark: str = "ideal-geometry heavy atoms") -> None:
    """Write an all-heavy-atom PDB. Day 2 feeds this to OpenMM PDBFile."""
    seq = str(structure["sequence"])
    residues = structure["residues"]
    with open(path, "w") as fh:
        fh.write(f"REMARK  {remark}\n")
        serial = 1
        for i, (aa, res) in enumerate(zip(seq, residues)):
            rn = _resolve(aa)
            for name, xyz in res.items():
                el = _ELEMENT.get(name[0], "C")
                nm = name if len(name) >= 4 else " " + name
                fh.write(
                    f"ATOM  {serial:>5d} {nm:<4s} {rn:>3s} A{i + 1:>4d}    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                    f"  1.00  0.00          {el:>2s}\n")
                serial += 1
        fh.write("TER\nEND\n")


# ============================================================================
# SECTION 3 -- ff14SB/GBn2 builder (was amber_hamiltonian.py)
# ============================================================================
KJ_PER_KCAL = 4.184
KCAL_PER_KJ = 1.0 / KJ_PER_KCAL

AMBER_TERMS = ("bond", "angle", "torsion", "nonbonded", "solvation")
_GROUP_OF = {"HarmonicBondForce": 1, "HarmonicAngleForce": 2,
             "PeriodicTorsionForce": 3, "CMAPTorsionForce": 3,
             "NonbondedForce": 4}
_TERM_GROUP = {"bond": 1, "angle": 2, "torsion": 3,
               "nonbonded": 4, "solvation": 5}
_RESTRAINT_GROUP = 6

RESTRAINED_BACKBONE = ("N", "CA", "C")
_NON_FRAME_ATOMS = ("O", "OXT")


class AmberHamiltonian(BudgetedEnergyModel):


    def __init__(self, sequence: str, representation,
                 weights: Optional[Dict[str, float]] = None,
                 restraint_k: float = 100.0,
                 minimization_steps: int = 50,
                 minimization_tolerance: float = 2.0,
                 cache_limit: int = 500_000,
                 platform_name: str = "CPU",
                 threads: int = 1,
                 reference_states: Optional[Sequence[int]] = None,
                 collapse_floor: float = -600.0,
                 collapse_mode: str = "geometry",
                 collapse_min_contact: float = 1.05,
                 eval_budget: Optional[int] = None):
        self.sequence = sequence.strip().upper()
        self.rep = representation
        if len(self.sequence) != self.rep.n_residues:
            raise ValueError(
                f"sequence length {len(self.sequence)} != representation "
                f"n_residues {self.rep.n_residues}")
        if getattr(self.rep, "is_lattice", False):
            raise ValueError(
                "AmberHamiltonian requires a full-backbone representation; "
                "the tetrahedral lattice provides CA only, which cannot be "
                "given an ff14SB topology. Use TorsionStateRepresentation.")

        self.restraint_k = float(restraint_k)
        self.threads = int(threads)
        self.minimization_steps = int(minimization_steps)
        self.minimization_tolerance = float(minimization_tolerance)
        self.collapse_floor = float(collapse_floor)
        # 2026-09-01: the energy-magnitude floor is a defect, not a guard. GB self-energy
        # scales as q^2, so a peptide with net charge +3 or more sits far below -600
        # kcal/mol while being perfectly ordinary geometry. Scanned over the 60-target
        # benchmark with the floor disabled, FIVE targets had their NEAR-NATIVE structure
        # returned as +inf (7OB2 q=+5 at -748.8, 1Y5C +5 at -819.0, 2FBS +4 at -851.9,
        # 2JUT +3 at -890.2, 6RQS +10 at -801.8), and on 2JUT 17 of 20 random decoys were
        # masked. Their radii of gyration match native to within 0.4 A and their closest
        # heavy-heavy contact is ~1.20 A -- an ordinary C=O bond, not a collapse. A real
        # Coulomb collapse has Rg of 2-3 A and two atoms on top of each other.
        #
        # `collapse_mode="geometry"` therefore tests what the guard was always meant to
        # test: a heavy-atom pair closer than `collapse_min_contact`, which no bonded pair
        # reaches. "energy" restores the old magnitude floor for reproducing the pinned
        # golden runs; "none" disables the guard entirely.
        self.collapse_mode = str(collapse_mode)
        self.collapse_min_contact = float(collapse_min_contact)
        self.n_collapsed = 0

        self._init_budget(cache_limit, eval_budget)
        self.t_build = 0.0
        self.t_minimize = 0.0
        self.t_energy = 0.0

        t0 = time.time()
        self._build_topology()
        self._build_system(platform_name)
        self._calibrate_hydrogens(reference_states)
        self.setup_time = time.time() - t0
        self._memo_tag = (
            self.sequence, int(self.n_atoms), str(platform_name),
            int(self.threads),
            hashlib.blake2b(
                np.ascontiguousarray(
                    np.array([c for *_, c in self._h_local],
                             dtype=np.float64).ravel()).tobytes(),
                digest_size=8).digest())

        self.weights = {t: 1.0 for t in AMBER_TERMS}
        if weights:
            self.weights.update(weights)

    @property
    def n_qubits(self) -> int:
        return self.rep.n_qubits

    @property
    def n_bits(self) -> int:
        return self.rep.n_bits

    def _reference_backbone(self) -> Dict[str, np.ndarray]:
        n = len(self.sequence)
        phi = np.full(n, math.radians(-120.0))
        psi = np.full(n, math.radians(130.0))
        phi[0] = geo.DEFAULT_PHI
        psi[n - 1] = geo.DEFAULT_PSI
        return geo.build_backbone(phi, psi)

    def _build_topology(self) -> None:
        seq = self.sequence
        ref = build_full_structure(seq, self._reference_backbone())

        top = app.Topology()
        chain = top.addChain("A")
        self._heavy_names: List[Tuple[int, str]] = []
        positions: List[np.ndarray] = []
        prev_C = None
        for i, aa in enumerate(seq):
            rn = geo.ONE_TO_THREE[aa]
            res = top.addResidue(rn, chain)
            objs = {}
            for nm, xyz in ref["residues"][i].items():
                objs[nm] = top.addAtom(
                    nm, app.element.Element.getBySymbol(nm[0]), res)
                self._heavy_names.append((i, nm))
                positions.append(np.asarray(xyz, dtype=float) * 0.1)  # A -> nm
            for a_, b_ in residue_bonds(rn):
                top.addBond(objs[a_], objs[b_])
            if "OXT" in objs:
                top.addBond(objs["C"], objs["OXT"])
            if prev_C is not None:
                top.addBond(prev_C, objs["N"])
            prev_C = objs["C"]

        self.forcefield = app.ForceField("amber14/protein.ff14SB.xml",
                                         "implicit/gbn2.xml")
        self._reference_platform = openmm.Platform.getPlatformByName("Reference")
        random.seed(0)
        modeller = app.Modeller(top, np.array(positions) * unit.nanometer)
        modeller.addHydrogens(self.forcefield,
                              platform=self._reference_platform)
        self.topology = modeller.topology
        self.n_atoms = self.topology.getNumAtoms()
        self._modeller_positions = np.array(
            modeller.positions.value_in_unit(unit.nanometer), dtype=float)

        self._index_topology()

    def _index_topology(self) -> None:
        self._elem = {a.index: (a.element.symbol if a.element else "X")
                      for a in self.topology.atoms()}
        self._resof = {a.index: a.residue.index for a in self.topology.atoms()}
        self._nameof = {a.index: a.name for a in self.topology.atoms()}

        self._heavy_index: Dict[Tuple[int, str], int] = {}
        self._hydrogens: List[int] = []
        for a in self.topology.atoms():
            if self._elem[a.index] == "H":
                self._hydrogens.append(a.index)
            else:
                self._heavy_index[(a.residue.index, a.name)] = a.index

        self._heavy_order = [self._heavy_index[k] for k in self._heavy_names]
        self._restraint_idx = [self._heavy_index[(i, nm)]
                               for i in range(len(self.sequence))
                               for nm in RESTRAINED_BACKBONE]

        adj: Dict[int, List[int]] = {}
        parent: Dict[int, int] = {}
        for bond in self.topology.bonds():
            i, j = bond[0].index, bond[1].index
            ei, ej = self._elem[i], self._elem[j]
            if ei == "H" and ej != "H":
                parent[i] = j
            elif ej == "H" and ei != "H":
                parent[j] = i
            elif ei != "H" and ej != "H":
                adj.setdefault(i, []).append(j)
                adj.setdefault(j, []).append(i)
        for k in adj:
            adj[k].sort()
        self._adj = adj
        self._hparent = parent
        missing = [h for h in self._hydrogens if h not in parent]
        if missing:
            raise RuntimeError(
                f"{len(missing)} hydrogens have no heavy parent in the "
                "topology; cannot build deterministic frames")

    def _frame_atoms(self, p: int) -> Tuple[int, int]:
        ri = self._resof[p]
        if self._nameof[p] == "N" and ri > 0:
            return self._heavy_index[(ri, "CA")], self._heavy_index[(ri - 1, "C")]
        same = [q for q in self._adj.get(p, [])
                if self._resof[q] == ri and self._nameof[q] not in _NON_FRAME_ATOMS]
        if len(same) >= 2:
            return same[0], same[1]
        if len(same) == 1:
            a = same[0]
            second = [q for q in self._adj.get(a, [])
                      if q != p and self._resof[q] == ri
                      and self._nameof[q] not in _NON_FRAME_ATOMS]
            if second:
                return a, second[0]
        raise RuntimeError(
            f"no rigid frame for {self._nameof[p]} in residue {ri}")

    @staticmethod
    def _frame(p: np.ndarray, a: np.ndarray, b: np.ndarray):
        e1 = a - p
        e1 = e1 / np.linalg.norm(e1)
        v = b - p
        e2 = v - np.dot(v, e1) * e1
        e2 = e2 / np.linalg.norm(e2)
        return e1, e2, np.cross(e1, e2)

    def _build_system(self, platform_name: str) -> None:
        self.system = self.forcefield.createSystem(
            self.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            rigidWater=False,
            removeCMMotion=False,
            implicitSolventKappa=0.0 / unit.nanometer)
        for f in self.system.getForces():
            f.setForceGroup(_GROUP_OF.get(f.__class__.__name__, 5))


        rest = openmm.CustomExternalForce(
            "0.5*k_rest*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        rest.addGlobalParameter("k_rest", 0.0)
        for nm in ("x0", "y0", "z0"):
            rest.addPerParticleParameter(nm)
        for idx in self._restraint_idx:
            rest.addParticle(idx, [0.0, 0.0, 0.0])
        rest.setForceGroup(_RESTRAINT_GROUP)
        self._rest_force = rest
        self.system.addForce(rest)
        self._k_rest_internal = self.restraint_k * KJ_PER_KCAL * 100.0

        self.platform = openmm.Platform.getPlatformByName(platform_name)
        names = set(self.platform.getPropertyNames())
        props = {}
        if "Threads" in names:
            props["Threads"] = str(int(self.threads))
        if "DeterministicForces" in names:
            props["DeterministicForces"] = "true"
        self.platform_properties = props
        self.context = openmm.Context(
            self.system,
            openmm.VerletIntegrator(0.001 * unit.picoseconds),
            self.platform, props)

    def _calibrate_hydrogens(self, reference_states) -> None:
        n = len(self.sequence)
        states = [1] * n if reference_states is None else list(reference_states)
        bits = self.rep.bitstring_from_states(states)
        heavy = self._heavy_positions(self.rep.build_coords(bits),
                                     chi1=self._chi1_of(bits))

        pos = np.array(self._modeller_positions, dtype=float)
        for k, idx in enumerate(self._heavy_order):
            pos[idx] = heavy[k]

        scratch = self.forcefield.createSystem(
            self.topology, nonbondedMethod=app.NoCutoff, constraints=None,
            rigidWater=False, removeCMMotion=False,
            implicitSolventKappa=0.0 / unit.nanometer)
        for i in range(self.n_atoms):
            if self._elem[i] != "H":
                scratch.setParticleMass(i, 0.0)
        calib_props = dict(self.platform_properties)
        if "Threads" in calib_props:
            calib_props["Threads"] = "1"
        ctx = openmm.Context(scratch, openmm.VerletIntegrator(0.001),
                             self.platform, calib_props)
        ctx.setPositions(pos * unit.nanometer)
        openmm.LocalEnergyMinimizer.minimize(ctx, 1e-4, 2000)
        relaxed = np.array(
            ctx.getState(getPositions=True).getPositions(asNumpy=True)
            .value_in_unit(unit.nanometer), dtype=float)
        del ctx

        self._h_local: List[Tuple[int, int, int, int, np.ndarray]] = []
        for h in self._hydrogens:
            p = self._hparent[h]
            a, b = self._frame_atoms(p)
            e1, e2, e3 = self._frame(relaxed[p], relaxed[a], relaxed[b])
            d = relaxed[h] - relaxed[p]
            self._h_local.append(
                (h, p, a, b,
                 np.array([np.dot(d, e1), np.dot(d, e2), np.dot(d, e3)])))

    def _chi1_of(self, bitstring: str) -> Optional[Dict[int, float]]:
        """Encoded chi1 angles, or None when the representation does not encode chi1.

        The atom *set* is independent of chi1 -- only positions move -- so `_heavy_names`,
        `_heavy_order` and the frozen hydrogen frames built in `_calibrate_hydrogens` stay
        valid. That is what makes putting chi1 in the bitstring safe here.
        """
        fn = getattr(self.rep, "chi1_degrees", None)
        if fn is None:
            return None
        return fn(bitstring) or None

    def _heavy_positions(self, backbone: Dict[str, np.ndarray],
                         chi1: Optional[Dict[int, float]] = None) -> np.ndarray:
        full = build_full_structure(self.sequence, backbone, chi1=chi1)
        out = np.empty((len(self._heavy_names), 3), dtype=float)
        for k, (i, nm) in enumerate(self._heavy_names):
            out[k] = full["residues"][i][nm]
        return out * 0.1

    def _assemble(self, heavy_nm: np.ndarray) -> np.ndarray:
        pos = np.empty((self.n_atoms, 3), dtype=float)
        for k, idx in enumerate(self._heavy_order):
            pos[idx] = heavy_nm[k]
        for h, p, a, b, c in self._h_local:
            e1, e2, e3 = self._frame(pos[p], pos[a], pos[b])
            pos[h] = pos[p] + c[0] * e1 + c[1] * e2 + c[2] * e3
        return pos

    def _evaluate(self, heavy_nm: np.ndarray, want_components: bool = False,
                  want_positions: bool = False):
        t0 = time.time()
        pos = self._assemble(heavy_nm)
        self.t_build += time.time() - t0

        ctx = self.context
        t0 = time.time()
        ctx.setPositions(pos * unit.nanometer)
        for k, idx in enumerate(self._restraint_idx):
            self._rest_force.setParticleParameters(k, idx, pos[idx].tolist())
        self._rest_force.updateParametersInContext(ctx)
        ctx.setParameter("k_rest", self._k_rest_internal)
        openmm.LocalEnergyMinimizer.minimize(
            ctx, self.minimization_tolerance, self.minimization_steps)
        self.t_minimize += time.time() - t0

        t0 = time.time()
        ctx.setParameter("k_rest", 0.0)
        # The geometry guard has to see the MINIMISED structure: the as-built ideal-geometry
        # structure routinely has close contacts that minimisation removes, so testing the
        # input coordinates would flag almost everything.
        need_pos = want_positions or self.collapse_mode == "geometry"
        state = ctx.getState(getEnergy=True, getPositions=need_pos)
        energy = state.getPotentialEnergy().value_in_unit(
            unit.kilocalorie_per_mole)
        final_pos = (np.array(state.getPositions(asNumpy=True)
                              .value_in_unit(unit.nanometer), dtype=float)
                     if need_pos else pos)
        if self._is_collapsed(final_pos, energy):
            # Return +inf so CVaR and best-seen can never select it. Applied uniformly to
            # native and predicted structures. See `collapse_mode` in __init__.
            self.n_collapsed += 1
            energy = float("inf")
        comp = None
        if want_components:
            comp = {}
            for term in AMBER_TERMS:
                s = ctx.getState(getEnergy=True, groups={_TERM_GROUP[term]})
                comp[term] = s.getPotentialEnergy().value_in_unit(
                    unit.kilocalorie_per_mole)
        out_pos = final_pos if want_positions else None
        self.t_energy += time.time() - t0
        return energy, comp, out_pos


    def _is_collapsed(self, pos: np.ndarray, energy: float) -> bool:
        """Is this structure a steric collapse rather than a low-energy fold?"""
        if self.collapse_mode == "none":
            return False
        if self.collapse_mode == "energy":
            return energy < self.collapse_floor
        if not np.isfinite(energy):
            return True
        heavy = pos[self._heavy_order] * 10.0                  # nm -> A
        d = np.linalg.norm(heavy[:, None, :] - heavy[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        return bool(d.min() < self.collapse_min_contact)

    def energy(self, bitstring: str) -> float:
        if len(bitstring) != self.rep.n_bits:
            raise ValueError(
                f"bitstring length {len(bitstring)} != n_bits {self.rep.n_bits}")
        hit = self._cache.get(bitstring)
        if hit is not None:
            return hit                  # cache hits are free
        self._charge()                  # raises BudgetExhausted before minimizing
        heavy = self._heavy_positions(self.rep.build_coords(bitstring),
                                     chi1=self._chi1_of(bitstring))
        e, _, _ = self._evaluate(heavy)
        if len(self._cache) < self._cache_limit:
            self._cache[bitstring] = e
        return e

    def components(self, bitstring: str) -> Dict[str, float]:

        if len(bitstring) != self.rep.n_bits:
            raise ValueError(
                f"bitstring length {len(bitstring)} != n_bits {self.rep.n_bits}")
        heavy = self._heavy_positions(self.rep.build_coords(bitstring),
                                     chi1=self._chi1_of(bitstring))
        e, comp, _ = self._evaluate(heavy, want_components=True)
        comp["total"] = e
        return comp

    def energy_from_coords(self, coords: Dict[str, np.ndarray],
                           phi=None, psi=None,
                           chi1: Optional[Dict[int, float]] = None) -> float:
        # Unlike the legacy model this does not scan chi1 rotamers: an OpenMM minimization
        # per rotamer combination is far too expensive, and the minimizer relaxes the
        # (unrestrained) sidechains anyway, so the starting rotamer matters much less here.
        heavy = self._heavy_positions(coords, chi1=chi1)
        e, _, _ = self._evaluate(heavy)
        return e

    def minimized_ca(self, bitstring: str) -> Tuple[np.ndarray, np.ndarray]:

        heavy = self._heavy_positions(self.rep.build_coords(bitstring),
                                     chi1=self._chi1_of(bitstring))
        pre = np.array([heavy[k] for k, (_, nm) in enumerate(self._heavy_names)
                        if nm == "CA"]) * 10.0
        _, _, pos = self._evaluate(heavy, want_positions=True)
        post = np.array([pos[self._heavy_index[(i, "CA")]]
                         for i in range(len(self.sequence))]) * 10.0
        return pre, post

    def reset_counters(self) -> None:
        super().reset_counters()
        self.t_build = self.t_minimize = self.t_energy = 0.0

    def describe(self) -> Dict[str, object]:
        return {
            "forcefield": "amber14/protein.ff14SB.xml + implicit/gbn2.xml",
            "units": "kcal/mol",
            "n_atoms": self.n_atoms,
            "n_hydrogens": len(self._hydrogens),
            "platform": self.platform.getName(),
            "platform_properties": dict(self.platform_properties),
            "restraint_k_kcal_per_mol_A2": self.restraint_k,
            "collapse_floor_kcal": self.collapse_floor,
            "minimization_steps": self.minimization_steps,
            "hydrogen_strategy": "frozen local frames (relaxed once)",
            "setup_time_s": self.setup_time,
        }


# ============================================================================
# SECTION 4 -- restrained refinement driver (was amber_refine.py)
# ============================================================================
#
# `k_restraint` is in kcal/mol/A^2 per restrained atom; the potential is
# 0.5*k*|r - r0|^2, so an atom pulled by an internal force F settles at d = F/k.
#
#     k =   1   WEAK      d ~ 0.1 - 1.0 A   stereochemistry repaired, fold free to drift
#     k =  10   MODERATE  d ~ 0.03 - 0.3 A  the default: relieves bond/angle strain,
#                                           far too little for a torsion to flip
#     k = 100   STRONG    d ~ 0.01 - 0.1 A  backbone pinned, sidechains/H only
#
# `steps=0` means "minimise until converged" (OpenMM's `maxIterations` convention), with
# `tolerance` in kJ/mol/nm on the RMS force. Determinism comes from the CPU platform at
# Threads=1 with DeterministicForces, and from the starting structure being a pure
# function of its inputs: no random seeding, no thermostat, no MD.

#: Restraint strengths in kcal/mol/A^2.
K_WEAK = 1.0
K_MODERATE = 10.0
K_STRONG = 100.0

#: Atoms held by the positional restraint. The carbonyl O is deliberately free: it carries
#: real ideal-geometry strain and is not needed to pin the fold.
BACKBONE_ATOMS = ("N", "CA", "C")

#: kcal/mol/A^2 -> kJ/mol/nm^2
_K_SCALE = KJ_PER_KCAL * 100.0

# ---------------------------------------------------------------- convergence gate
#: SPRINT 16.  `LocalEnergyMinimizer.minimize` is asked to run to convergence and NEVER
#: told whether it got there.  On the 126-target tuning instrument four targets (1D6X,
#: 1MF6, 2NB7, 7BX2) end minimisation above 1000 kcal/mol -- 1MF6 at 8.9e8 -- and were
#: silently scored into every published mean.  Those four alone move the exact rotated-
#: frame null (which is ZERO by construction) from -0.0005 to +0.0117 kcal-free angstrom.
#:
#: THE RULE, DECLARED BEFORE IT WAS APPLIED (s16/energy_gate.py records the declaration
#: and its timestamp).  A minimisation is CONVERGED iff, with the restraint switched off:
#:
#:     final potential energy <= CONVERGE_MAX_KCAL      (default 1000.0 kcal/mol)
#:
#: One threshold, no target-specific tuning, no native-derived quantity, and the same
#: 1000 kcal/mol number the codebase already uses for its bond+angle strain gate. It is
#: reported, never silently applied: `refine_coords` always returns the structure, and
#: sets `converged` / `converge_reason` so the CONSUMER decides.  Any statistic quoted
#: from a gated set must print the excluded count.
CONVERGE_MAX_KCAL = 1000.0


def convergence_flags(energy: float, energy_initial: float,
                      max_kcal: float = CONVERGE_MAX_KCAL) -> Dict[str, object]:
    """The pre-declared gate.  Pure function of the returned energies."""
    e = float(energy)
    ok = np.isfinite(e) and e <= float(max_kcal)
    if not np.isfinite(e):
        why = "non-finite final energy"
    elif e > float(max_kcal):
        why = "final energy %.4g > %.4g kcal/mol" % (e, max_kcal)
    else:
        why = ""
    return {"converged": bool(ok), "converge_reason": why,
            "converge_max_kcal": float(max_kcal),
            "energy_drop": float(energy_initial) - e}


# ---------------------------------------------------------------- builder LRU
#: Bounded, unlike the original dict. A 126-target sweep used to retain 126 live OpenMM
#: Contexts plus their Systems and topologies; on a 16 GB box shared with three sibling
#: runs that is the difference between finishing and swapping. Eviction tears the Context
#: down deliberately rather than waiting for the collector.
BUILDER_CACHE_SIZE = 8

_BUILDERS: "OrderedDict[tuple, AmberHamiltonian]" = __import__(
    "collections").OrderedDict()


def _drop_builder(h: "AmberHamiltonian") -> None:
    """Release an OpenMM Context and its System deliberately."""
    for attr in ("context", "system", "_rest_force", "platform",
                 "_reference_platform", "forcefield"):
        try:
            delattr(h, attr)
        except AttributeError:
            pass


def builder_for(sequence: str, rep, platform_name: str = "CPU",
                threads: int = 1) -> "AmberHamiltonian":
    """Cached ff14SB/GBn2 builder for one sequence + representation shape.

    Construction costs ~2.8 s (topology, hydrogen addition, and a one-off hydrogen-frame
    calibration), so it is cached; every later candidate for the same target is then a
    minimisation and nothing else. This cache is why there is no per-candidate System
    rebuild to remove -- see the module docstring.

    `threads` > 1 lets the CPU platform parallelise CustomGBForce (2.2x at 4 threads) at
    the cost of bit-exactness; see the module docstring before using it.
    """
    key = (sequence.strip().upper(), int(rep.n_bits), int(rep.n_states),
           platform_name, int(threads))
    h = _BUILDERS.get(key)
    if h is not None:
        _BUILDERS.move_to_end(key)
        return h
    memory_guard()
    h = AmberHamiltonian(sequence, rep, restraint_k=K_MODERATE,
                         minimization_steps=0, platform_name=platform_name,
                         collapse_floor=float("-inf"), threads=threads)
    _BUILDERS[key] = h
    while len(_BUILDERS) > BUILDER_CACHE_SIZE:
        _, evicted = _BUILDERS.popitem(last=False)
        _drop_builder(evicted)
    return h


# ---------------------------------------------------------------- result memo
#: `_run` is a pure function of (builder, heavy positions, k_restraint, steps, tolerance)
#: -- the same fact that makes Context reuse valid. Retrieval pools repeat candidates
#: (identical fragment windows recur across targets and across k), and a repeat currently
#: costs a full 5.9 s converged minimisation to rediscover a number already computed.
#: Keyed on the exact bytes of the input coordinates, so a hit is exact, never approximate.
MEMO_SIZE = 512
_MEMO: "OrderedDict[tuple, dict]" = __import__("collections").OrderedDict()
_MEMO_HITS = 0
_MEMO_MISSES = 0


def _memo_key(H: "AmberHamiltonian", heavy_nm: np.ndarray, k_restraint: float,
              steps: int, tolerance: float, components: bool) -> tuple:
    a = np.ascontiguousarray(heavy_nm, dtype=np.float64)
    return (H._memo_tag, hashlib.blake2b(a.tobytes(), digest_size=16).digest(),
            float(k_restraint), int(steps), float(tolerance), bool(components))


def _memo_copy(out: dict) -> dict:
    """A caller may mutate what it gets back; the cache must not see that."""
    c = dict(out)
    for k, v in c.items():
        if isinstance(v, np.ndarray):
            c[k] = v.copy()
        elif isinstance(v, dict):
            c[k] = {kk: (vv.copy() if isinstance(vv, np.ndarray) else vv)
                    for kk, vv in v.items()}
        elif isinstance(v, list):
            c[k] = list(v)
    return c


def cache_stats() -> Dict[str, int]:
    """Builder and memo occupancy, and the memo hit/miss counts."""
    return {"builders": len(_BUILDERS), "memo": len(_MEMO),
            "memo_hits": _MEMO_HITS, "memo_misses": _MEMO_MISSES}


def clear_cache() -> None:
    """Tear down every cached Context and empty the memo."""
    global _MEMO_HITS, _MEMO_MISSES
    for h in _BUILDERS.values():
        _drop_builder(h)
    _BUILDERS.clear()
    _MEMO.clear()
    _MEMO_HITS = _MEMO_MISSES = 0


# ---------------------------------------------------------------- the run itself
def _split(H: "AmberHamiltonian", pos_nm: np.ndarray) -> Dict[str, object]:
    """Full-system positions (nm) -> heavy-atom views in angstrom."""
    n = len(H.sequence)
    ang = pos_nm * 10.0
    bb: Dict[str, np.ndarray] = {}
    for nm in ("N", "CA", "C", "O"):
        bb[nm] = np.array([ang[H._heavy_index[(i, nm)]] for i in range(n)])
    heavy_names = list(H._heavy_names)
    heavy = np.array([ang[H._heavy_index[k]] for k in heavy_names])
    return {"ca": bb["CA"], "backbone": bb, "heavy": heavy,
            "heavy_names": heavy_names}


def _energy_kcal(H: "AmberHamiltonian") -> float:
    """Potential energy with the restraint switched off, in kcal/mol."""
    H.context.setParameter("k_rest", 0.0)
    return H.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)


def _run(H: "AmberHamiltonian", heavy_nm: np.ndarray, k_restraint: float,
         steps: int, tolerance: float, components: bool) -> Dict[str, object]:
    """One restrained minimisation. Byte-for-byte the original arithmetic.

    Nothing here is reordered or fused. The 28 ms of Python around a 5900 ms C++ call is
    not worth the reproducibility risk of touching it.
    """
    t0 = time.time()
    pos = H._assemble(heavy_nm)                      # heavy + frozen-frame hydrogens
    ctx = H.context
    ctx.setPositions(pos * unit.nanometer)

    # Anchor every restrained atom at its own input position.
    for j, idx in enumerate(H._restraint_idx):
        H._rest_force.setParticleParameters(j, idx, pos[idx].tolist())
    H._rest_force.updateParametersInContext(ctx)

    e0 = _energy_kcal(H)
    if steps >= 0:                      # steps < 0 => single point, no minimisation
        ctx.setParameter("k_rest", float(k_restraint) * _K_SCALE)
        openmm.LocalEnergyMinimizer.minimize(ctx, float(tolerance), int(steps))

    ctx.setParameter("k_rest", 0.0)
    state = ctx.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    out_nm = np.asarray(
        state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)

    out = _split(H, out_nm)
    out.update({"energy": float(energy), "energy_initial": float(e0),
                "k_restraint": float(k_restraint), "steps": int(steps),
                "tolerance": float(tolerance), "n_atoms": int(H.n_atoms),
                "positions_nm": out_nm})
    # Displacement of the restrained set: how far the fold actually moved.
    ridx = H._restraint_idx
    d = (out_nm[ridx] - pos[ridx]) * 10.0
    out["restraint_rmsd"] = float(np.sqrt((d ** 2).sum(1).mean()))
    out["restraint_max"] = float(np.sqrt((d ** 2).sum(1)).max())
    # SPRINT 16 convergence gate.  Reported, never silently applied; adds no arithmetic
    # to the minimisation itself, so every prior number reproduces bit-for-bit.
    out.update(convergence_flags(energy, e0))
    if components:
        out["components"] = {
            t: ctx.getState(getEnergy=True, groups={_TERM_GROUP[t]})
                  .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            for t in AMBER_TERMS}
    out["wall"] = time.time() - t0
    return out


def _run_memo(H: "AmberHamiltonian", heavy_nm: np.ndarray, k_restraint: float,
              steps: int, tolerance: float, components: bool,
              memo: bool = True) -> Dict[str, object]:
    global _MEMO_HITS, _MEMO_MISSES
    if not memo or MEMO_SIZE <= 0:
        return _run(H, heavy_nm, k_restraint, steps, tolerance, components)
    t0 = time.time()
    key = _memo_key(H, heavy_nm, k_restraint, steps, tolerance, components)
    hit = _MEMO.get(key)
    if hit is not None:
        _MEMO.move_to_end(key)
        _MEMO_HITS += 1
        out = _memo_copy(hit)
        # `wall` must be THIS call's cost. Returning the original 5.9 s for a 20 us
        # cache hit would silently corrupt every timing roll-up downstream.
        out["wall"] = time.time() - t0
        out["memo_hit"] = True
        return out
    _MEMO_MISSES += 1
    out = _run(H, heavy_nm, k_restraint, steps, tolerance, components)
    out["memo_hit"] = False
    _MEMO[key] = _memo_copy(out)
    while len(_MEMO) > MEMO_SIZE:
        _MEMO.popitem(last=False)
    return out


# ---------------------------------------------------------------- public API
def refine_coords(sequence: str, rep, coords: Dict[str, np.ndarray],
                  k_restraint: float = K_MODERATE, steps: int = 0,
                  tolerance: float = 1.0, chi1: Optional[Dict[int, float]] = None,
                  platform_name: str = "CPU",
                  components: bool = False,
                  threads: int = 1, memo: bool = True) -> Dict[str, object]:
    """Refine from a backbone dict (N, CA, C, O, CB in angstrom).

    This is the coordinate entry point and the one the pipeline uses. `refine` projects
    continuous torsions onto the discrete state library first, which changes the
    experiment; do not substitute it here.
    """
    H = builder_for(sequence, rep, platform_name, threads)
    return _run_memo(H, H._heavy_positions(coords, chi1=chi1),
                     k_restraint, steps, tolerance, components, memo)


def refine(sequence: str, rep, states: Sequence[int],
           k_restraint: float = K_MODERATE, steps: int = 0,
           tolerance: float = 1.0, platform_name: str = "CPU",
           components: bool = False,
           threads: int = 1, memo: bool = True) -> Dict[str, object]:
    """Restrained all-atom ff14SB/GBn2 refinement of one DISCRETE candidate.

    Returns a dict with `ca`, `backbone`, `heavy`, `heavy_names`, `positions_nm`,
    `energy` (restraint excluded), `energy_initial`, `restraint_rmsd`, `restraint_max`
    and `wall`; `components` adds the five-term decomposition.
    """
    if isinstance(states, str):
        bits = states
    else:
        bits = rep.bitstring_from_states(np.asarray(states, dtype=int))
    H = builder_for(sequence, rep, platform_name, threads)
    heavy = H._heavy_positions(rep.build_coords(bits), chi1=H._chi1_of(bits))
    return _run_memo(H, heavy, k_restraint, steps, tolerance, components, memo)


def refine_ca(sequence: str, rep, ca: np.ndarray,
              k_restraint: float = K_MODERATE, steps: int = 0,
              tolerance: float = 1.0, platform_name: str = "CPU",
              sweeps: int = 12, seed: int = 0,
              components: bool = False, threads: int = 1) -> Dict[str, object]:
    """Refine from a bare CA trace; the discrete states are fitted to it first."""
    import floor as _floor
    ca = np.asarray(ca, float)
    n = rep.n_residues
    rng = np.random.default_rng(seed)
    best_s, best = _floor.descend(rep, ca, rng.integers(0, rep.n_states, n),
                                  sweeps=sweeps)
    out = refine(sequence, rep, best_s, k_restraint, steps, tolerance,
                 platform_name, components, threads)
    out["states"] = np.asarray(best_s, int)
    out["fit_rmsd"] = float(best)
    return out


def single_point(sequence: str, rep, states, platform_name: str = "CPU",
                 components: bool = False, threads: int = 1) -> Dict[str, object]:
    """ff14SB + GBn2 energy of the built structure with NO minimisation at all."""
    return refine(sequence, rep, states, k_restraint=0.0, steps=-1,
                  tolerance=1e9, platform_name=platform_name,
                  components=components, threads=threads)


# ---------------------------------------------------------------- parallel map
#: Per-worker resident set, measured: ~230 MB for interpreter + numpy + OpenMM + one
#: 237-atom Context. Used to size the pool against free physical memory.
WORKER_MB = 260.0

_W: Dict[str, object] = {}


def _worker_init(sequence, rep, platform_name, threads):
    """One Context per process, built once and reused for every chunk."""
    # Windows spawns rather than forks, so the child re-imports from scratch. Without the
    # repo root on its path `import core.amber` fails, `refine_many` falls back to serial,
    # and the only symptom is that the parallelism silently did not happen.
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root not in sys.path:
        sys.path.insert(0, root)
    import core.amber as A
    _W["mod"] = A
    _W["args"] = (sequence, rep, platform_name, threads)
    A.builder_for(sequence, rep, platform_name, threads)


def _worker_run(payload):
    A = _W["mod"]
    seq, rep, platform_name, threads = _W["args"]
    idx, coords_list, k_restraint, steps, tolerance, components, chi1 = payload
    out = []
    for j, c in zip(idx, coords_list):
        r = A.refine_coords(seq, rep, c, k_restraint=k_restraint, steps=steps,
                            tolerance=tolerance, chi1=chi1, components=components,
                            platform_name=platform_name, threads=threads)
        r.pop("heavy_names", None)          # rebuilt on the parent side; big and constant
        out.append((j, r))
    return out


def plan_workers(n_tasks: int, workers: Optional[int] = None,
                 reserve_cores: int = 1, limit: float = MEMORY_LIMIT_PERCENT
                 ) -> int:
    """How many OpenMM worker processes this box can actually afford right now.

    Bounded by three things and the smallest wins: the number of candidates, the cores
    left after `reserve_cores` (sibling agents are running), and free physical memory at
    `WORKER_MB` each with `limit` as the ceiling.
    """
    if n_tasks <= 1:
        return 1
    if workers is not None:
        return max(1, min(int(workers), n_tasks))
    cores = max(1, (os.cpu_count() or 2) - max(0, int(reserve_cores)))
    if os.name == "nt":
        s = _MEMORYSTATUSEX()
        s.dwLength = ctypes.sizeof(_MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s))
        total_mb = s.ullTotalPhys / 1e6
        used_mb = total_mb * memory_percent() / 100.0
        headroom_mb = max(0.0, total_mb * limit / 100.0 - used_mb)
        by_mem = int(headroom_mb // WORKER_MB)
    else:
        by_mem = cores
    return max(1, min(n_tasks, cores, by_mem))


def refine_many(sequence: str, rep, coords_list: Sequence[Dict[str, np.ndarray]],
                k_restraint: float = K_MODERATE, steps: int = 0,
                tolerance: float = 1.0, chi1: Optional[Dict[int, float]] = None,
                platform_name: str = "CPU", components: bool = False,
                workers: Optional[int] = None, threads: int = 1,
                reserve_cores: int = 1) -> List[Dict[str, object]]:
    """Refine a whole candidate pool. Results in input order, bit-identical to serial.

    This is the one large win the profile actually supports. 99.5% of `refine_coords` is
    one `LocalEnergyMinimizer` call and candidates are independent, so the only thing
    left to exploit is that there are many of them. Each worker holds its own Context at
    Threads=1 and runs the identical code path, so every returned float is the float the
    serial loop would have produced -- this parallelises the work, it does not approximate
    it.

    Falls back to an in-process loop for one candidate, when only one worker is
    affordable, or if the pool cannot be started (a raised worker count on a box already
    at the memory ceiling is a reason to run serially, not to fail).

    `heavy_names` is dropped from worker results and restored from the local builder: it
    is a constant list of ~230 tuples per structure and pickling it per candidate costs
    more than the refinement saves.
    """
    coords_list = list(coords_list)
    n = len(coords_list)
    if n == 0:
        return []
    # An explicit `workers=` overrides the core and memory *planning*, but not the
    # ceiling: a measured run of this benchmark left the box at 96% because the caller
    # asked for six workers and nothing re-checked. Sibling processes share this machine,
    # so the gate is unconditional and running serially is the correct fallback.
    nw = plan_workers(n, workers, reserve_cores)
    if nw > 1:
        try:
            memory_guard()
        except MemoryError as exc:
            warnings.warn(f"{exc}; running serially", RuntimeWarning)
            nw = 1
    if nw <= 1:
        return [refine_coords(sequence, rep, c, k_restraint, steps, tolerance,
                              chi1, platform_name, components, threads)
                for c in coords_list]

    # Chunk so each worker pays the ~2.8 s builder cost once and then amortises it.
    chunks: List[tuple] = []
    per = max(1, (n + nw - 1) // nw)
    for a in range(0, n, per):
        idx = list(range(a, min(a + per, n)))
        chunks.append((idx, [coords_list[j] for j in idx], k_restraint, steps,
                       tolerance, components, chi1))

    results: List[Optional[dict]] = [None] * n
    try:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(
                max_workers=nw, initializer=_worker_init,
                initargs=(sequence, rep, platform_name, threads)) as ex:
            for got in ex.map(_worker_run, chunks):
                for j, r in got:
                    results[j] = r
    except Exception as exc:                                        # noqa: BLE001
        warnings.warn(f"parallel refine failed ({type(exc).__name__}: {exc}); "
                      "falling back to serial", RuntimeWarning)
        return [refine_coords(sequence, rep, c, k_restraint, steps, tolerance,
                              chi1, platform_name, components, threads)
                for c in coords_list]

    H = builder_for(sequence, rep, platform_name, threads)
    names = list(H._heavy_names)
    for r in results:
        if r is not None:
            r["heavy_names"] = list(names)
    return [r for r in results]                                     # type: ignore[misc]
