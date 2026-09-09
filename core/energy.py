"""Consolidated Legacy physics: the knowledge-based 11-term energy and its field.

This is `energy_terms` + `legacy_field` + `torsion_lib2` + `priors` in one module. The
eleven terms, their `DEFAULT_WEIGHTS`/`FITTED_WEIGHTS` combination, the MJ contact
potential, the DSSP H-bond model, the cooperativity terms and the Ramachandran table are
lifted byte-for-byte; `energy_components` returns exactly what it returned before.

Where the time goes
Measured on real BLOSUM-pool structures for 1A13 (n=14), per structure:

    energy_components            0.306 ms
      hbond_terms                0.106 ms   (dssp_energy_matrix + greedy match)
      steric_term                0.073 ms
      contact                    0.019 ms
      solvation                  0.018 ms
      compactness                0.015 ms
      electrostatic              0.012 ms
      coop_helix + coop_sheet    0.003 ms
      aromatic                   0.000 ms
      sum of terms               0.246 ms
      per-call overhead          0.060 ms   (d_cb, asarray, dict construction)

For scale, one `core.amber.refine_coords` on the same structure costs 5928 ms. Legacy is
19,000x cheaper than the AMBER half and is not the pipeline's bottleneck. It is still
worth batching, because the VQE objective evaluates it ~200k times per fold, where 0.3 ms
becomes a minute per target.

The brief's premise -- "per-pair interaction arithmetic currently written with Python
loops" -- was already stale when this consolidation started. `steric_term` caches its pair
list, radii and exemption mask in `_steric_layout` and costs one gather, one norm and one
dot; `contact_term` is a single `np.dot`; `solvation_term` accumulates coordination with
`np.bincount` off the triangle `energy_components` has already computed. The only real
Python loop left is the greedy one-donor-one-acceptor match in `hbond_terms`, and it is
sequential by definition: each accepted bond removes a donor and an acceptor from
contention. It cannot be vectorised without changing the tie-break, which is why
`BatchLegacy._hbonds` -- which does exactly that -- is documented as a different
tie-break rather than an optimisation of the same one.

What is actually new here
`components_batch` scores a whole pool in one call. It batches every term whose reduction
order it can preserve and loops only the sequential ones, so it is not an approximation of
`energy_components`: measured over real pools, all eleven terms agree to 0.0 absolute
(see `tests/test_energy.py`, which reports the number rather than assuming it). The win is
the per-call overhead and the numpy dispatch, which is most of the cost at these sizes.

`BatchLegacy.terms_from_coords` is the *other* batched path and it is deliberately not the
same thing: it substitutes a vectorised greedy H-bond match and the CB-CB aromatic
fallback so it can drive generation. Its differences from the scalar model are documented
per-term and checked by `verify`. Use `components_batch` when you need the scalar model's
numbers and `BatchLegacy` when you need a generation field.
"""
import math
import os
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import peptide_db as pdb
import protein_geometry as geo
import representations as reps

__all__ = [
    # the 11-term model
    "energy_components", "components_batch", "total_from_components",
    "TERM_NAMES", "DEFAULT_WEIGHTS", "FITTED_WEIGHTS", "WEIGHT_ORIGIN",
    "FREE_WEIGHTS", "COUPLED_TERMS", "SEPARABLE_TERMS", "FILTER_TERMS",
    "steric_term", "contact_term", "hbond_terms", "coop_helix_term",
    "coop_sheet_term", "solvation_term", "electrostatic_term", "aromatic_term",
    "torsion_term", "compactness_term", "backtracking_term", "ring_frame",
    "rama_penalty", "switch", "sequence_arrays", "pair_index",
    "aromatic_indices", "clear_sequence_cache",
    "MJ_CORRECTED", "MJ_RAW", "BURIAL", "CHARGE", "AROMATIC",
    # generation field
    "BatchLegacy", "LegacyField", "TERMS", "SUBSETS", "verify",
    # torsion libraries
    "ContextLibrary", "PerResidueTorsion", "library_for", "make",
    # priors
    "DistogramPrior", "TorsionMRF", "pair_features", "fit_distogram_model",
    "BIN_EDGES", "NBINS",
]

# SECTION 1 -- the Legacy 11-term energy (was energy_terms.py)
MJ_ORDER = ["C", "M", "F", "I", "L", "V", "W", "Y", "A", "G",
            "T", "S", "N", "Q", "D", "E", "H", "R", "K", "P"]

_MJ_RAW = [
    [-5.44, -4.99, -5.80, -5.50, -5.83, -4.96, -4.95, -4.16, -3.57, -3.16, -3.11, -2.86, -2.59, -2.85, -2.41, -2.27, -3.60, -2.57, -1.95, -3.07],
    [-4.99, -5.46, -6.56, -6.02, -6.41, -5.32, -5.55, -4.91, -3.94, -3.39, -3.51, -3.03, -2.95, -3.30, -2.57, -2.89, -3.98, -3.12, -2.48, -3.45],
    [-5.80, -6.56, -7.26, -6.84, -7.28, -6.29, -6.16, -5.66, -4.81, -4.13, -4.28, -4.02, -3.75, -4.10, -3.48, -3.56, -4.77, -3.98, -3.36, -4.25],
    [-5.50, -6.02, -6.84, -6.54, -7.04, -6.05, -5.78, -5.25, -4.58, -3.78, -4.03, -3.75, -3.52, -3.67, -3.17, -3.27, -4.14, -3.63, -3.01, -3.76],
    [-5.83, -6.41, -7.28, -7.04, -7.37, -6.48, -6.14, -5.67, -4.91, -4.16, -4.34, -4.08, -3.75, -4.04, -3.40, -3.59, -4.54, -4.03, -3.37, -4.20],
    [-4.96, -5.32, -6.29, -6.05, -6.48, -5.52, -5.18, -4.62, -4.04, -3.38, -3.46, -3.30, -3.07, -3.28, -2.83, -2.90, -3.58, -3.07, -2.49, -3.32],
    [-4.95, -5.55, -6.16, -5.78, -6.14, -5.18, -5.06, -4.66, -3.82, -3.42, -3.22, -3.07, -3.07, -3.11, -2.84, -2.99, -3.98, -3.41, -2.69, -3.73],
    [-4.16, -4.91, -5.66, -5.25, -5.67, -4.62, -4.66, -4.17, -3.36, -3.01, -3.01, -2.78, -2.83, -2.97, -2.76, -2.79, -3.52, -3.16, -2.60, -3.19],
    [-3.57, -3.94, -4.81, -4.58, -4.91, -4.04, -3.82, -3.36, -2.72, -2.31, -2.32, -2.01, -1.84, -1.89, -1.70, -1.51, -2.41, -1.83, -1.31, -2.03],
    [-3.16, -3.39, -4.13, -3.78, -4.16, -3.38, -3.42, -3.01, -2.31, -2.24, -2.08, -1.82, -1.74, -1.66, -1.59, -1.22, -2.15, -1.72, -1.15, -1.87],
    [-3.11, -3.51, -4.28, -4.03, -4.34, -3.46, -3.22, -3.01, -2.32, -2.08, -2.12, -1.96, -1.88, -1.90, -1.80, -1.74, -2.42, -1.90, -1.31, -1.90],
    [-2.86, -3.03, -4.02, -3.75, -4.08, -3.30, -3.07, -2.78, -2.01, -1.82, -1.96, -1.67, -1.58, -1.49, -1.63, -1.48, -2.11, -1.62, -1.05, -1.57],
    [-2.59, -2.95, -3.75, -3.52, -3.75, -3.07, -3.07, -2.83, -1.84, -1.74, -1.88, -1.58, -1.68, -1.71, -1.68, -1.51, -2.08, -1.64, -1.21, -1.53],
    [-2.85, -3.30, -4.10, -3.67, -4.04, -3.28, -3.11, -2.97, -1.89, -1.66, -1.90, -1.49, -1.71, -1.54, -1.46, -1.42, -1.98, -1.80, -1.29, -1.73],
    [-2.41, -2.57, -3.48, -3.17, -3.40, -2.83, -2.84, -2.76, -1.70, -1.59, -1.80, -1.63, -1.68, -1.46, -1.21, -1.02, -2.32, -2.29, -1.68, -1.33],
    [-2.27, -2.89, -3.56, -3.27, -3.59, -2.90, -2.99, -2.79, -1.51, -1.22, -1.74, -1.48, -1.51, -1.42, -1.02, -0.91, -2.15, -2.27, -1.80, -1.26],
    [-3.60, -3.98, -4.77, -4.14, -4.54, -3.58, -3.98, -3.52, -2.41, -2.15, -2.42, -2.11, -2.08, -1.98, -2.32, -2.15, -3.05, -2.16, -1.35, -2.25],
    [-2.57, -3.12, -3.98, -3.63, -4.03, -3.07, -3.41, -3.16, -1.83, -1.72, -1.90, -1.62, -1.64, -1.80, -2.29, -2.27, -2.16, -1.55, -0.59, -1.70],
    [-1.95, -2.48, -3.36, -3.01, -3.37, -2.49, -2.69, -2.60, -1.31, -1.15, -1.31, -1.05, -1.21, -1.29, -1.68, -1.80, -1.35, -0.59, -0.12, -0.97],
    [-3.07, -3.45, -4.25, -3.76, -4.20, -3.32, -3.73, -3.19, -2.03, -1.87, -1.90, -1.57, -1.53, -1.73, -1.33, -1.26, -2.25, -1.70, -0.97, -1.75],
]


def _build_mj_corrected() -> Dict[Tuple[str, str], float]:
    idx = {aa: i for i, aa in enumerate(MJ_ORDER)}
    self_e = {aa: _MJ_RAW[idx[aa]][idx[aa]] for aa in MJ_ORDER}
    return {(a, b): _MJ_RAW[idx[a]][idx[b]] - 0.5 * (self_e[a] + self_e[b])
            for a in MJ_ORDER for b in MJ_ORDER}

MJ_CORRECTED = _build_mj_corrected()
MJ_RAW = {(a, b): _MJ_RAW[MJ_ORDER.index(a)][MJ_ORDER.index(b)]
          for a in MJ_ORDER for b in MJ_ORDER}

# Burial scale
#: Fauchere-Pliska octanol/water pi (1983), all twenty residues.
#:
#: This replaces Kyte-Doolittle, and the replacement is the single largest accuracy fix
#: in this rewrite. KD is a membrane-spanning propensity scale and places the aromatics
#: on the *hydrophilic* side (W -0.9, Y -1.3, F +2.8 only by comparison). Any peptide
#: whose hydrophobic core is an aromatic cluster is then penalised for forming it, and
#: for a peptide with no aliphatic core every KD value can be negative -- in which case
#: `solvation_term` becomes purely expansive and actively drives the search to an
#: extended chain. Fauchere-Pliska puts W and Y at +2.25 and +0.96, so burying an
#: aromatic cluster is rewarded, which is what the physics says.
#:
#: This is a property of the scale, not of any one sequence: it changes the sign of the
#: solvation term for every aromatic-core peptide in `dataset.CANDIDATE_PDB_IDS`.
BURIAL: Dict[str, float] = {
    "W": 2.25, "I": 1.80, "F": 1.79, "L": 1.70, "C": 1.54, "M": 1.23, "V": 1.22,
    "Y": 0.96, "P": 0.72, "A": 0.31, "T": 0.26, "H": 0.13, "G": 0.00, "S": -0.04,
    "Q": -0.22, "N": -0.60, "E": -0.64, "D": -0.77, "K": -0.99, "R": -1.01,
}
#: Trp, so the normalised scale spans [-0.45, 1.0] and is dimensionless.
BURIAL_NORM = 2.25

CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.5}

#: Aromatic residues, as a class. Ring-ring stacking is real physics that the MJ
#: contact term represents only weakly (its weighted spread is 0.29 against torsion's
#: 2.27), and it is the dominant tertiary interaction in short hairpin peptides
#: generally -- not in any one of them specifically.
AROMATIC = frozenset("FYWH")

#: Coulomb constant, kcal A / (mol e^2). The previous electrostatic term omitted this
#: entirely, so a salt bridge at 4 A scored -0.076 kcal/mol against a physical -1 to -3
#: and the term had the smallest spread of the seven (0.011). It was inert.
COULOMB = 332.0637
#: Effective dielectric for a solvent-exposed peptide. Combined with the existing
#: exp(-d/8) screening this is a coarse Debye-Huckel treatment, so the term's weight is
#: still a calibration target -- but its *magnitude* is now physical.
DIELECTRIC = 40.0

#: van der Waals radii, A. Used by the all-atom backbone clash test.
VDW_RADIUS = {"N": 1.55, "C": 1.70, "O": 1.52, "S": 1.80}
#: Soft-sphere prefactor. Contacts are penalised below SOFTNESS * (r_i + r_j); 0.80
#: leaves room for real close contacts (a C-C pair may legitimately reach 3.4 A) while
#: still catching interpenetration.
SOFTNESS = 0.80

#: |i-j| at or above which an H-bond counts as long-range. A definition, not a weight.
HB_LONGRANGE_SEP = 5

#: Maximum run count `coop_helix_term` may return. NOT a weight -- a weight scales the
#: reward everywhere, this bounds where it stops accumulating, and no setting of
#: `calibrate_weights` reproduces it.
#:
#: Why it is needed. `coop_helix_term` advances by one residue, so its ceiling is n-5 --
#: one n-turn per residue along the whole chain, growing without bound. `coop_sheet_term`
#: steps by two AND is capped by strand length, so its ceiling is ~n/4. That is 5 vs 2 at
#: n=10 and 15 vs 8 at n=20, which makes the docstring below (before this change) false:
#: the two terms are symmetric in FORM and roughly 5:1 asymmetric in MAGNITUDE. Capping the
#: helix run makes the code do what it already claimed to do.
#:
#: Why 2, and why the helix term only. Measured on the 12 natives in `CANDIDATE_PDB_IDS`:
#: coop_sheet reaches at most 3 (1LE1, 1LE3) and coop_helix reaches 8 (1V4Z). A cap of 2 on
#: the helix brings the two into practical parity without touching the sheet term -- capping
#: both at 2 costs the natives 26.0 kcal/mol against 22.0 for helix-only, and the extra 4.0
#: comes entirely off real hairpin natives, which is the opposite of the intent.
#:
#: Measured effect, chignolin, `diagnose_energy_model.py --protein 1UAO --states 4`
#: (all 4,194,304 structures -- these are the numbers that tool reprints):
#:     shipped   global min -17.018 @ 5.27 A   top-100 mean RMSD 5.07   enrichment -0.27
#:     cap = 2   global min -13.125 @ 1.96 A   top-100 mean RMSD 2.47   enrichment +2.33
#: Restricted to the 2,574,772 clash-free structures the numbers are 5.124 -> 2.567 and
#: 0.000 -> +2.557; the two differ only because the full set includes clashing structures.
#:
#: Note Spearman(energy, RMSD) is +0.1629 either way, unchanged to four decimals. The cap
#: acts only on structures carrying a helical run of 3+, which is 1.56% of the space, so it
#: reorders the TOP of the ranking without moving the global rank correlation. Anyone
#: judging this change by Spearman will conclude, wrongly, that it does nothing.
#: No native regresses: the helical natives (1V4Z, 1L2Y, 2JOF) already rank at percentile
#: 0.00 among sampled feasible structures and stay there, because nothing sampled is close.
#:
#: WARNING for anyone evaluating this. `--energy-quality` CANNOT see this change. Random
#: sampling essentially never produces a helical run of 3 or more (0.00% of 4000 draws on
#: chignolin, 0.03-0.18% on the others), so sampled Spearman and enrichment are identical
#: with and without the cap. Judge it on `diagnose_energy_model.py` or a real SA/VQE run,
#: where the search actually reaches the helix basin.
COOP_RUN_CAP = 2

#: Re-exported from `protein_geometry` so the two copies of the DSSP form cannot drift.
#: Below HB_MIN_ON the donor and acceptor heavy atoms are interpenetrating, and nothing else
#: in the model says so: `_steric_layout` exempts every hetero N/O pair, and the amide H is
#: not in the steric layout at all, so the `-1/dOH` term is unbounded below. As dOH -> 0.5
#: (the old numerical guard) that term alone reaches -55.8 kcal/mol against a maximum
#: possible weighted steric of 27.0, so the clash can never win at any weight. Measured: the
#: SA arm at default budget drove into it on 2 of 3 seeds, scoring -58 against an ideal
#: helix's -28.
HB_MIN_ON = geo.HB_MIN_ON
HB_E_FLOOR = geo.HB_E_FLOOR

#: Aromatic ring-ring geometry. Two wells rather than one, because pi-stacking has two
#: distinct favourable arrangements and a single distance well cannot tell them apart:
#:
#:   parallel-displaced  ring planes roughly parallel, centroids ~3.6-4.2 A apart
#:   T-shaped (edge-to-face)  planes roughly perpendicular, centroids ~5.0-5.5 A
#:
#: `theta` is the angle between ring normals folded into [0, 90] degrees, since a normal's
#: sign is arbitrary. The relative depths are a modelling choice under the term's single
#: calibrated weight; the geometry is not.
COOP_LADDER_MIN_SEP =2
AROM_PD_DIST = 4.00
AROM_PD_DEPTH = 1.00
AROM_T_DIST = 5.20
AROM_T_DEPTH = 0.70
AROM_DIST_WIDTH = 1.10
AROM_ANGLE_WIDTH = 35.0
#: Fallback CB-CB well, used for aromatics whose ring this repo cannot build (His) or when
#: chi1 is not encoded so no ring geometry exists.
AROM_CB_DIST = 5.50
AROM_CB_WIDTH = 1.80

# Weights
#: Starting weights. `physical` entries are fixed by the units of the term they scale;
#: `empirical` entries are free and should be set by
#: `energy_quality.calibrate_weights` over a train split of sequences.
DEFAULT_WEIGHTS: Dict[str, float] = {
    "steric": 4.0,
    "contact": 1.0,
    "hbond_local": 1.0,
    "hbond_longrange": 3.0,
    "coop_helix": 2.0,
    "coop_sheet": 2.0,
    "solvation": 0.5,
    "electrostatic": 1.0,
    "aromatic": 0.8,
    "torsion": 0.15,
    "compactness": 0.4,
}

TERM_NAMES = list(DEFAULT_WEIGHTS.keys())

#: Which weights `calibrate_weights` is allowed to move. The rest are pinned by the units
#: of the term they scale.
#:
#: `hbond_longrange` is here rather than fixed at 3.0 deliberately. Weighting the
#: long-range ladder above local bonds is right for a beta hairpin and gives a helix
#: nothing -- every helical H-bond is i->i+4 -- so a hardcoded multiplier is a standing bet
#: that the targets are hairpins. Most of `dataset.CANDIDATE_PDB_IDS` are, but 1DU1 is a
#: charged helix and the trp-cages are mixed, so the multiplier has to be measured on a
#: train split and reported on held-out clusters like any other free parameter. Splitting
#: `hbond` into two terms also means `experiment_energy_ablation` produces separate
#: `no_hbond_local` and `no_hbond_longrange` arms, which is how the hairpin bias becomes
#: visible instead of assumed.
FREE_WEIGHTS = ("aromatic", "torsion", "compactness", "solvation",
                "hbond_longrange", "coop_helix", "coop_sheet")

#: Terms that couple residues: they cannot be evaluated without a three-dimensional
#: arrangement, so they are the only ones carrying tertiary-structure information.
COUPLED_TERMS = ("contact", "hbond_local", "hbond_longrange", "coop_helix",
                 "coop_sheet", "solvation", "electrostatic", "aromatic",
                 "compactness")

#: Terms that are a sum of independent per-residue contributions. Their minimum is found
#: residue by residue, so a landscape they dominate has no fold in it -- which is exactly
#: how the pre-2026-07-28 model ended up anti-correlated with correctness. The invariant
#: worth defending is COUPLED variance > SEPARABLE variance among feasible structures;
#: `energy_quality` measures it and `validation` gates on it.
SEPARABLE_TERMS = ("torsion",)

#: `steric` is in neither list: it is a feasibility filter, not a discriminator. It
#: carried 87% of the variance across *all* random structures while being exactly zero for
#: both the native fold and the global minimum, so among the clash-free structures a
#: search actually explores it is flat. Any variance accounting that does not exclude it
#: balances against a term that does not vary where it matters.
FILTER_TERMS = ("steric",)

WEIGHT_ORIGIN = {
    "steric": "physical  (must dominate; hard-core overlap is forbidden)",
    "contact": "reference (MJ-corrected potential used at unit weight)",
    "hbond_local": "reference (DSSP electrostatic model, kcal/mol, at unit weight)",
    "hbond_longrange": "empirical (hairpin-ladder emphasis; MUST be calibrated -- a fixed "
                       "multiplier is a bet that the targets are hairpins)",
    "coop_helix": "empirical (consecutive n-turns; MUST be calibrated)",
    "coop_sheet": "empirical (consecutive ladder rungs; MUST be calibrated)",
    "solvation": "empirical (burial scale is dimensionless; sets burial vs contact)",
    "electrostatic": "physical  (screened Coulomb, with the 332 prefactor)",
    "aromatic": "empirical (ring-ring well depth, kcal/mol-ish)",
    "torsion": "empirical (Ramachandran basin depth; small because the state library "
               "is already restricted to favourable basins -- see docstring)",
    "compactness": "empirical (one-sided Rg restraint)",
}

# Ramachandran
_RAMA_BASINS = [
    (-63.0, -42.0, 28.0, 1.00),
    (-120.0, 130.0, 40.0, 0.90),
    (-75.0, 145.0, 40.0, 0.75),
    (-85.0, 100.0, 32.0, 0.70),
    (-100.0, -15.0, 30.0, 0.45),
    (75.0, 35.0, 30.0, 0.40),
]
_HELIX_FORMERS = set("AELMQKRH")
_SHEET_FORMERS = set("VIFYTWC")


def _ang_diff(a: float, b: float) -> float:
    return ((a - b + 180.0) % 360.0) - 180.0


@lru_cache(maxsize=100_000)
def rama_penalty(aa: str, phi_rad: float, psi_rad: float) -> float:
    """Ramachandran basin penalty for one residue.

    Memoized: the torsion representations draw phi/psi from a fixed library of discrete
    values, so across a whole search this has only n_residues x n_states distinct
    arguments but is called once per residue per structure. Pure function of its
    arguments, so the cache returns bit-identical values.

    NOTE on the weight this term carries. `representations` restricts phi/psi to a
    curated library of *favourable* basins, so a term whose job is to keep torsions
    physical has little left to do -- while being a sum of independent per-residue
    penalties, so its minimum is reached by choosing each residue's best basin in
    isolation. That minimum is locally ideal everywhere and globally meaningless, and at
    weight 1.0 it accounted for +6.329 of native chignolin's +7.963 penalty relative to
    the global minimum. Hence `DEFAULT_WEIGHTS["torsion"] = 0.15`: it survives as a
    tiebreaker carrying amino-acid specificity, and the *library* (per-residue-class in
    `representations`) now carries the rest.
    """
    pd, sd = math.degrees(phi_rad), math.degrees(psi_rad)
    score = 0.0
    for k, (pc, sc, sig, depth) in enumerate(_RAMA_BASINS):
        d2 = _ang_diff(pd, pc) ** 2 + _ang_diff(sd, sc) ** 2
        w = depth
        if k == 0 and aa in _HELIX_FORMERS:
            w += 0.5
        if k == 1 and aa in _SHEET_FORMERS:
            w += 0.5
        if k == 3 and aa != "G":
            w *= 0.2
        score += w * math.exp(-d2 / (2.0 * sig * sig))
    e = 1.0 - score
    if aa == "P":
        if pd < -90.0:
            e += 0.03 * (-90.0 - pd)
        elif pd > -50.0:
            e += 0.03 * (pd - (-50.0))
    if aa == "G":
        e -= 0.2
    return e


def switch(d: np.ndarray, d0: float, dc: float) -> np.ndarray:
    """Cosine switch: 1 below d0, 0 above dc, smooth in between."""
    d = np.asarray(d, dtype=float)
    s = np.zeros_like(d)
    s[d <= d0] = 1.0
    mid = (d > d0) & (d < dc)
    if np.any(mid):
        s[mid] = 0.5 * (1.0 + np.cos(math.pi * (d[mid] - d0) / (dc - d0)))
    return s

# Per-sequence caches
_SEQ_CACHE: Dict[Tuple[str, bool], Tuple] = {}


def sequence_arrays(sequence: str, use_corrected_mj: bool = True):
    """Cache burial, charge, and pairwise MJ arrays for a sequence.

    Returns ``(burial, charge, mj)``. The first element was Kyte-Doolittle hydropathy
    before this rewrite; it is now the Fauchere-Pliska scale (see `BURIAL`).
    """
    key = (sequence, use_corrected_mj)
    hit = _SEQ_CACHE.get(key)
    if hit is not None:
        return hit
    table = MJ_CORRECTED if use_corrected_mj else MJ_RAW
    burial = np.array([BURIAL.get(a, 0.0) for a in sequence])
    q = np.array([CHARGE.get(a, 0.0) for a in sequence])
    mj = np.array([[table.get((a, b), 0.0) for b in sequence] for a in sequence])
    _SEQ_CACHE[key] = (burial, q, mj)
    return burial, q, mj


@lru_cache(maxsize=4096)
def aromatic_indices(sequence: str) -> Tuple[int, ...]:
    return tuple(i for i, a in enumerate(sequence) if a in AROMATIC)


@lru_cache(maxsize=64)
def pair_index(n: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cached upper-triangle pair indices and separations for `n` residues.

    ``np.triu_indices`` was being rebuilt on every energy evaluation. The arrays are
    returned from a cache, so treat them as read-only -- nothing here mutates them.
    """
    di, dj = np.triu_indices(n, 1)
    return di, dj, (dj - di)


def clear_sequence_cache() -> None:
    _SEQ_CACHE.clear()
    _STERIC_LAYOUT_CACHE.clear()
    rama_penalty.cache_clear()
    aromatic_indices.cache_clear()
    pair_index.cache_clear()


# Terms
def steric_term(coords: Dict[str, np.ndarray], sequence: str,
                rings: Optional[Dict[int, Dict[str, np.ndarray]]] = None,
                min_sep: int = 2) -> float:
    """Soft-sphere overlap over every available heavy backbone atom, plus CB and rings.

    Replaces the previous CA-CA (< 3.8 A) and CB-CB (< 3.4 A) pair test, which had two
    problems. It never looked at N, C or O, so backbone carbonyls could interpenetrate
    for free -- the main geometric constraint that should force an aromatic pair into a
    stacked rather than an overlapping arrangement. And its 3.8 A CA threshold is below
    the geometric floor for |i-j| = 2 (2 * 3.8 * sin(40 deg) ~ 4.9 A), so in the torsion
    representation the sep-2 test could never fire at all.

    N/O pairs are exempt: those are exactly the H-bond donor/acceptor pairs, which
    legitimately approach 2.8-3.0 A, and penalising them would put this term in direct
    conflict with `hbond_terms`.

    `rings`, when supplied, adds aromatic ring atoms. This is what stops two rings simply
    passing through each other -- previously nothing did, because `sidechains` built the
    rings and the energy never looked at them. Note the soft-sphere limit for a carbon pair
    is 0.80 * (1.70 + 1.70) = 2.72 A while parallel-displaced stacking puts the closest
    ring atoms at ~3.5 A, so this forbids interpenetration without forbidding the stack.

    The pair list, element radii and exemption mask depend only on the *layout* -- which
    atom names are present, the residue count, and which ring atoms exist -- never on the
    coordinates. They are cached by `_steric_layout`, so a structure costs one gather, one
    norm and one dot product. Rebuilding them per call made this the most expensive term in
    the model once ring atoms took the pair count from 45 to ~2000.
    """
    names = tuple(k for k in ("N", "CA", "C", "O", "CB") if k in coords)
    if not names:
        return 0.0
    n_res = len(coords[names[0]])
    if n_res != len(sequence):
        raise ValueError(f"coordinate count {n_res} != sequence length {len(sequence)}")

    # Insertion order, which `sidechains.ring_atom_names` fixes to template order, so the
    # gather below lines up with the cached layout. Deterministic, and avoids a per-call
    # sort of the atom names.
    ring_keys = tuple((i, nm) for i, atoms_i in (rings or {}).items() for nm in atoms_i)
    ii, jj, limit = _steric_layout(names, n_res, ring_keys, min_sep)
    if ii.size == 0:
        return 0.0

    blocks = [np.asarray(coords[k], dtype=float) for k in names]
    if ring_keys:
        blocks.append(np.array([rings[i][nm] for i, nm in ring_keys], dtype=float))
    atoms = np.concatenate(blocks, axis=0)

    over = np.maximum(0.0, limit - np.linalg.norm(atoms[ii] - atoms[jj], axis=1))
    return float(over @ over)

_STERIC_LAYOUT_CACHE: Dict[tuple, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}


def _steric_layout(names: Tuple[str, ...], n_res: int,
                   ring_keys: Tuple[Tuple[int, str], ...], min_sep: int):
    """Cached ``(ii, jj, contact_limit)`` for one atom layout. Read-only."""
    key = (names, n_res, ring_keys, min_sep)
    hit = _STERIC_LAYOUT_CACHE.get(key)
    if hit is not None:
        return hit

    elem: List[str] = []
    res_of: List[int] = []
    for k in names:
        # CB is carbon; Gly has no real CB but carries a virtual one, which is the right
        # excluded volume for a CA-only pseudo-atom anyway.
        elem.extend([("N" if k == "N" else "O" if k == "O" else "C")] * n_res)
        res_of.extend(range(n_res))
    for i, nm in ring_keys:
        elem.append(nm[0] if nm[0] in VDW_RADIUS else "C")
        res_of.append(i)

    elem_arr = np.asarray(elem)
    radius = np.array([VDW_RADIUS[e] for e in elem])
    is_no = np.isin(elem_arr, ("N", "O"))
    res_arr = np.asarray(res_of)

    ii, jj = np.triu_indices(len(elem), 1)
    keep = np.abs(res_arr[ii] - res_arr[jj]) >= min_sep
    # exempt N...O and O...N (H-bond partners)
    keep &= ~(is_no[ii] & is_no[jj] & (elem_arr[ii] != elem_arr[jj]))
    ii, jj = ii[keep], jj[keep]
    out = (ii, jj, SOFTNESS * (radius[ii] + radius[jj]))
    _STERIC_LAYOUT_CACHE[key] = out
    return out


def contact_term(mj: np.ndarray, sep: np.ndarray, di: np.ndarray,
                 dj: np.ndarray, d_cb: np.ndarray) -> float:
    """MJ-corrected contact energy over CB pairs at |i-j| >= 3."""
    m3 = sep >= 3
    return float(np.dot(mj[di[m3], dj[m3]], switch(d_cb[m3], 4.5, 8.5)))


def hbond_terms(coords: Dict[str, np.ndarray],
                desolvation_cost: float = 1.0
                ) -> Tuple[float, float, Tuple[Tuple[int, int], ...]]:
    """DSSP electrostatic H-bond energies, split into local and long-range sums.

    Returns ``(local, longrange, matched_pairs)`` where each pair is ``(donor, acceptor)``.
    The pair list is what makes the cooperativity terms possible -- see
    `coop_helix_term` / `coop_sheet_term`.

    Vectorised greedy one-donor-one-acceptor matching. `energy_components` reports the two
    sums as the separate weighted terms `hbond_local` and `hbond_longrange`, and -- the
    part that matters -- does *not* divide either by the residue count. `HB_LONGRANGE_SEP`
    is the boundary.
    """
    if not all(k in coords for k in ("N", "C", "O")):
        return 0.0, 0.0, ()      # representation has no backbone (e.g. the lattice)

    n = len(coords["N"])
    # The DSSP form itself lives in `protein_geometry.dssp_energy_matrix` -- ONE copy,
    # shared with `dssp_hbonds`. It was written out here as well until this was factored;
    # see that function for why a second copy is a correctness hazard rather than a
    # duplication nit.
    E, ok = geo.dssp_energy_matrix(coords, min_sep=2, cutoff=-0.5)

    donors, acceptors = np.where(ok)
    if len(donors) == 0:
        return 0.0, 0.0, ()

    energies = E[donors, acceptors]
    # Stable: HB_E_FLOOR manufactures exact ties at -4.0 and this greedy match is
    # order-dependent, so an unstable sort could map one bitstring to two energies across
    # platforms -- which is exactly the invariant budget.py's cache relies on.
    order = np.argsort(energies, kind="stable")
    donor_used = np.zeros(n, dtype=bool)
    acc_used = np.zeros(n, dtype=bool)
    local = lr = 0.0
    matched = []
    for k in order:
        i, j = int(donors[k]), int(acceptors[k])
        if donor_used[i] or acc_used[j]:
            continue
        donor_used[i] = True
        acc_used[j] = True
        e = float(energies[k]) + desolvation_cost
        if e >= 0.0:
            continue
        matched.append((i, j))
        if abs(i - j) < HB_LONGRANGE_SEP:
            local += e
        else:
            lr += e
    return local, lr, tuple(matched)


def coop_helix_term(pairs: Tuple[Tuple[int, int], ...]) -> float:
    """Reward *consecutive* n-turns, which is what a helix actually is.

    `hbond_local` is a sum over independent matched pairs, so two scattered i->i+4 bonds
    score exactly the same as two adjacent ones. That degeneracy is the reason an additive
    potential cannot tell a helix from a compact tangle with the same bond count. DSSP
    itself defines a helix by *repeating* n-turns, not by their number.

    In DSSP's convention the donor is C-terminal to the acceptor, so a helical bond has
    ``donor - acceptor`` in {3, 4}; the run condition is that ``(d+1, a+1)`` is also bonded.
    An n-residue run scores n-1, against 0 for the same bonds scattered.

    The count is capped at `COOP_RUN_CAP`. Uncapped it grows as n-5 -- one n-turn per
    residue, along the whole chain -- while `coop_sheet_term` is bounded near n/4 by strand
    length, so the longer the peptide the more the objective pays for being a helix for
    reasons of counting rather than physics. See `COOP_RUN_CAP` for the measurements. The
    cap is on the run count, not on which bonds count as helical, so a long helix still
    scores the full `hbond_local` reward for every one of its bonds; what it stops earning
    is additional *cooperativity* bonus past the point where the helix is established.
    """
    s = {(d, a) for d, a in pairs if 3 <= d - a <= 4}
    runs = sum(1 for (d, a) in s if (d + 1, a + 1) in s)
    return -float(min(runs, COOP_RUN_CAP))


def coop_sheet_term(pairs: Tuple[Tuple[int, int], ...],
                    min_sep: int = 2) -> float:
    """Reward *consecutive* rungs of an antiparallel ladder.

    The same additivity problem as `coop_helix_term`, and the one that most likely bounds
    RMSD on a hairpin: `hbond_longrange` cannot distinguish two scattered long-range bonds
    from two adjacent rungs forming a correct register. In an antiparallel ladder the
    H-bonded pairs step as ``(i, j) -> (i+2, j-2)``, so that offset is the run condition.

    Purely topological -- no amino acid, no distance scale, no parameter. Requiring a
    ladder also constrains the turn at the small-|i-j| end, which is the other half of
    register error -- which is why `min_sep` is 2 and NOT `HB_LONGRANGE_SEP`. That
    threshold classifies bond *energy* into local vs long-range; used here it discards
    every pair with |d - a| < 5, i.e. the rungs nearest the turn, truncating each ladder
    from the inside. A 3-rung ladder then scored -1 instead of -2, and on a 10-mer the
    term had almost nothing left to reward. 2 is the floor `hbond_terms` already enforces
    when matching, so every matched pair is admissible and the (i+2, j-2) step does all
    the discrimination.

    Only antiparallel is modelled: hairpins are antiparallel, and parallel sheets do not
    occur in peptides of this length.
    """
    s = {(d, a) for d, a in pairs if abs(d - a) >= min_sep}
    return -float(sum(1 for (d, a) in s if (d + 2, a - 2) in s))


def solvation_term(burial: np.ndarray, d_cb: np.ndarray,
                   di: np.ndarray, dj: np.ndarray, n: int) -> float:
    """Burial reward, proportional to the Fauchere-Pliska value of each residue.

    `burial` used to be Kyte-Doolittle. See the `BURIAL` docstring for why that made
    aromatic-core peptides unfoldable by construction.

    Coordination numbers are accumulated from the upper-triangle distances
    `energy_components` has already computed, rather than from a fresh (n, n) matrix plus a
    `fill_diagonal`. Identical result -- the switch is symmetric and the triangle excludes
    the diagonal -- at half the distance work and no square allocation.
    """
    s = switch(d_cb, 6.0, 10.0)
    coord = (np.bincount(di, weights=s, minlength=n)
             + np.bincount(dj, weights=s, minlength=n))
    return float(-(burial / BURIAL_NORM) @ coord)


def ring_frame(ring_xyz) -> Tuple[np.ndarray, np.ndarray]:
    """Centroid and unit normal of a planar ring.

    The normal is the least-variance direction of the atom cloud, which for a ring that is
    planar to 1e-6 A (as `sidechains` builds them, pinned by
    `validation.test_sidechain_rings_planar`) is the ring normal to machine precision.

    Uses `eigh` on the 3x3 scatter matrix rather than an SVD of the (k, 3) coordinates:
    same eigenvector, order-independent, and several times cheaper -- which matters because
    this runs once per aromatic residue per energy evaluation. `eigh` returns eigenvalues
    ascending, so column 0 is the smallest-variance direction.
    """
    P = np.asarray(ring_xyz, dtype=float)
    c = P.mean(axis=0)          # computed once; this used to be evaluated twice
    M = P - c
    _, vecs = np.linalg.eigh(M.T @ M)
    return c, vecs[:, 0]


def aromatic_term(sequence: str, coords: Dict[str, np.ndarray],
                  rings: Optional[Dict[int, Dict[str, np.ndarray]]] = None,
                  min_sep: int = 3) -> float:
    """Ring-ring attraction between aromatic sidechains (F, Y, W, H).

    Two regimes, and which one applies is a property of the *encoding*:

    * **Real geometry** when `rings` supplies ring atoms for both partners -- available
      once chi1 is encoded (`representations.TorsionStateRepresentation.build_rings`). The
      potential is then a function of centroid distance *and* the angle between ring
      normals, with separate wells for parallel-displaced and T-shaped stacking. This is
      the point of the chi1 bits: three terms in this model (burial, long-range H-bond, and
      this one) all say "bring the aromatics together", but only this one can say *in what
      orientation*, and orientation is what separates a 2 A prediction from a 4 A one.
    * **CB-CB fallback** otherwise -- His, whose imidazole template is not implemented, or
      any representation without chi1 bits. A single soft distance well: crude, but honest
      about the fact that without chi1 the ring's direction is not determined.

    This is the one *new* term rather than a correction, and the one most at risk of
    encoding a preconception about what these peptides look like. It is therefore in
    `TERM_NAMES`, so `experiment_energy_ablation` automatically produces a `no_aromatic`
    arm, and its weight is a `FREE_WEIGHTS` entry to be set on a train split and reported
    on held-out clusters. If it does not survive that, delete it.
    """
    idx = aromatic_indices(sequence)
    if len(idx) < 2:
        return 0.0
    rings = rings or {}
    CB = np.asarray(coords.get("CB", coords["CA"]), dtype=float)

    # dict insertion order is the template order from `sidechains.ring_atom_names`, so this
    # is deterministic without a per-call sort; `ring_frame` is order-independent anyway.
    frames = {i: ring_frame(np.array(list(rings[i].values()), dtype=float))
              for i in idx if i in rings and len(rings[i]) >= 3}

    total = 0.0
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            if j - i < min_sep:
                continue
            if i in frames and j in frames:
                ci, ni = frames[i]
                cj, nj = frames[j]
                d = float(np.linalg.norm(ci - cj))
                # Normals are sign-ambiguous, so fold the angle into [0, 90] degrees.
                cos = min(1.0, abs(float(np.dot(ni, nj))))
                theta = math.degrees(math.acos(cos))
                pd = (AROM_PD_DEPTH
                      * math.exp(-((d - AROM_PD_DIST) / AROM_DIST_WIDTH) ** 2)
                      * math.exp(-(theta / AROM_ANGLE_WIDTH) ** 2))
                t = (AROM_T_DEPTH
                     * math.exp(-((d - AROM_T_DIST) / AROM_DIST_WIDTH) ** 2)
                     * math.exp(-((theta - 90.0) / AROM_ANGLE_WIDTH) ** 2))
                total -= pd + t
            else:
                d = float(np.linalg.norm(CB[i] - CB[j]))
                total -= math.exp(-((d - AROM_CB_DIST) / AROM_CB_WIDTH) ** 2)
    return float(total)


def electrostatic_term(q: np.ndarray, sep: np.ndarray, di: np.ndarray,
                       dj: np.ndarray, d_cb: np.ndarray) -> float:
    """Screened Coulomb between formal sidechain charges.

    Now carries the `COULOMB / DIELECTRIC` prefactor the previous version omitted, which
    made it 13-40x too weak and gave it the smallest spread of any term.

    Masked before the exponential rather than after: most peptides have a handful of charged
    pairs out of n(n-1)/2, and the old form evaluated `exp` over every pair to then discard
    almost all of it.
    """
    qq = q[di] * q[dj]
    mask = (sep >= 2) & (qq != 0.0)
    if not np.any(mask):
        return 0.0
    d = d_cb[mask]
    # The 2.0 A floor guards the 1/r singularity only; the screening exponential keeps the
    # true distance, exactly as before. Clamping both would silently change the term for
    # sub-2 A pairs -- which are severe clashes, but changing them is still a change.
    return float(np.sum((COULOMB / DIELECTRIC) * (qq[mask] / np.maximum(d, 2.0))
                        * np.exp(-d / 8.0)))


def torsion_term(sequence: str, phi: np.ndarray, psi: np.ndarray) -> float:
    """Sum of per-residue Ramachandran penalties. Separable by construction -- see the
    weight note in `rama_penalty`."""
    return float(sum(rama_penalty(sequence[i], phi[i], psi[i])
                     for i in range(len(sequence))))


def compactness_term(CA: np.ndarray, n: int) -> float:
    """One-sided radius-of-gyration restraint.

    Penalises expansion past the expected Rg for a folded chain of this length but not
    over-compaction, which is `steric_term`'s job. A two-sided well made the two terms
    fight, and at the old weight of 0.05 an extended chain paid only 0.69 for being ~4 A
    too large while gaining ~6 in torsion -- the one term that knew the answer should be
    compact was the second-weakest in the model.
    """
    rg = geo.radius_of_gyration(CA)
    target = 2.2 * (n ** 0.38)
    return float(max(0.0, rg - target) ** 2)


def backtracking_term(rep, bitstring: str) -> float:
    if not getattr(rep, "is_lattice", False):
        return 0.0
    from representations import LATTICE_DIRECTIONS
    dirs = [np.array(LATTICE_DIRECTIONS[b]) for b in rep.bond_directions(bitstring)]
    return float(sum(1.0 for a, b in zip(dirs, dirs[1:])
                     if float(np.dot(a, b)) < -2.5))


def energy_components(sequence: str,
                     coords: Dict[str, np.ndarray],
                     phi: Optional[np.ndarray] = None,
                     psi: Optional[np.ndarray] = None,
                     use_corrected_mj: bool = True,
                     rings: Optional[Dict[int, Dict[str, np.ndarray]]] = None
                     ) -> Dict[str, float]:
    """Unweighted term values for one structure.

    All coordinates must be in ANGSTROMS. `representations` guarantees this for both
    representations; before this rewrite the tetrahedral lattice returned coordinates in
    lattice units (bond norm sqrt(3) ~ 1.73 where CA-CA should be 3.80), which made
    every threshold here wrong by 2.19x on that representation.

    `rings` maps residue index -> aromatic ring atoms, from
    `representations.TorsionStateRepresentation.build_rings`. Supplying it turns the
    aromatic term into a real orientation-dependent pi-stacking potential and gives ring
    atoms excluded volume; omitting it falls back to CB proxies.

    Note `hbond_local` and `hbond_longrange` are separate weighted terms. They used to be
    combined -- and divided by the chain length -- which suppressed the only long-range
    structural signal by 10x at N=10 and let a helix's local bonds substitute for a
    hairpin's ladder. Keeping them separate also makes the long-range emphasis a calibrated
    parameter rather than a hardcoded bet that every target is a hairpin.

    `coop_helix` and `coop_sheet` then add the one thing an additive pairwise potential
    structurally cannot express: cooperativity. They are rewarded by the same mechanism at
    the same default weight, but that alone does NOT make them topology-neutral, and this
    docstring used to claim it did. Their ceilings differ by construction: the helix run
    advances one residue at a time over the whole chain (max n-5) while the ladder steps two
    at a time and is bounded by strand length (max ~n/4), i.e. 5 vs 2 at n=10 and 15 vs 8 at
    n=20. Equal weights on unequal ceilings is a preference. `COOP_RUN_CAP` bounds the helix
    run so the two are comparable in practice; see that constant for the measurements.
    """
    n = len(sequence)
    CA = np.asarray(coords["CA"], dtype=float)
    CB = np.asarray(coords.get("CB", coords["CA"]), dtype=float)
    if len(CA) != n:
        raise ValueError(f"coordinate count {len(CA)} != sequence length {n}")

    burial, q, mj = sequence_arrays(sequence, use_corrected_mj)

    di, dj, sep = pair_index(n)          # cached; read-only
    d_cb = np.linalg.norm(CB[di] - CB[dj], axis=1)

    hb_local, hb_lr, hb_pairs = hbond_terms(coords)

    return {
        "steric": steric_term(coords, sequence, rings=rings),
        "contact": contact_term(mj, sep, di, dj, d_cb),
        "hbond_local": hb_local,
        "hbond_longrange": hb_lr,
        "coop_helix": coop_helix_term(hb_pairs),
        "coop_sheet": coop_sheet_term(hb_pairs),
        "solvation": solvation_term(burial, d_cb, di, dj, n),
        "electrostatic": electrostatic_term(q, sep, di, dj, d_cb),
        "aromatic": aromatic_term(sequence, coords, rings=rings),
        "torsion": (0.0 if phi is None else torsion_term(sequence, phi, psi)),
        "compactness": compactness_term(CA, n),
    }


def total_from_components(components: Dict[str, float],
                         weights: Dict[str, float]) -> float:
    return float(sum(weights.get(k, 0.0) * components.get(k, 0.0)
                     for k in TERM_NAMES))

# SECTION 2 -- batched generation field (was legacy_field.py)
#: Term order. Same names as `energy_terms.TERM_NAMES` so weights are interchangeable.
TERMS = list(TERM_NAMES)

#: Named term families, for `LegacyField(use=...)`. `all` is the shipped Legacy model;
#: the rest isolate the families that could plausibly be independent of a CA-CA distance
#: prior. `hb` is pure backbone geometry with no amino-acid identity in it at all, which
#: makes it the cleanest test of independence: nothing in it can be a restatement of what
#: a sequence-conditioned distance predictor already knows.
SUBSETS = {
    "all": None,
    "hb": ("hbond_local", "hbond_longrange", "coop_helix", "coop_sheet"),
    "burial": ("contact", "solvation"),
    "packing": ("steric", "contact", "solvation", "compactness"),
}


class BatchLegacy:
    """Per-sequence tables built once; `terms(states)` scores a whole batch."""

    def __init__(self, sequence: str, representation, use_corrected_mj: bool = True):
        self.seq = sequence
        self.rep = representation
        self.n = n = len(sequence)
        self.burial, self.q, self.mj = sequence_arrays(sequence, use_corrected_mj)
        di, dj, sep = pair_index(n)
        self.di, self.dj, self.sep = di, dj, sep.astype(float)
        self.mj_pair = self.mj[di, dj]
        self.qq = self.q[di] * self.q[dj]
        self.elec_mask = (sep >= 2) & (self.qq != 0.0)
        self.m3 = sep >= 3
        self.rows = np.arange(n)
        self.idx_sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        # steric layout, shared with the scalar model
        self.st_names = ("N", "CA", "C", "O", "CB")
        self.st_ii, self.st_jj, self.st_lim = _steric_layout(
            self.st_names, n, (), 2)
        # aromatic CB pairs at |i-j| >= 3
        idx = aromatic_indices(sequence)
        pairs = [(i, j) for a, i in enumerate(idx) for j in idx[a + 1:] if j - i >= 3]
        self.ar_i = np.array([p[0] for p in pairs], int)
        self.ar_j = np.array([p[1] for p in pairs], int)
        # rama table: per-residue, per-state, exact for a fixed library
        self.rama = np.array(
            [[rama_penalty(sequence[i], representation._phi[i, s],
                              representation._psi[i, s])
              for s in range(representation.n_states)] for i in range(n)])
        self.rg_target = 2.2 * (n ** 0.38)

    # ---------------------------------------------------------------- pieces
    def _steric(self, coords) -> np.ndarray:
        A = np.concatenate([coords[k] for k in self.st_names], axis=1)   # (B, 5n, 3)
        d = np.linalg.norm(A[:, self.st_ii] - A[:, self.st_jj], axis=2)
        over = np.maximum(0.0, self.st_lim[None, :] - d)
        return (over ** 2).sum(1)

    def _hbonds(self, coords):
        """Greedy one-donor-one-acceptor DSSP match, batched.

        Returns ``(local, longrange, bond)`` with `bond` the (B, n, n) boolean matrix of
        matched (donor, acceptor) pairs -- what the cooperativity terms read.
        """
        N, C, O = coords["N"], coords["C"], coords["O"]
        B, n = N.shape[0], self.n
        H = np.full((B, n, 3), np.nan)
        if n >= 2:
            d = C[:, :-1] - O[:, :-1]
            nd = np.linalg.norm(d, axis=2)
            ok0 = nd > 1e-6
            sh = N[:, 1:] + np.divide(d, nd[..., None], out=np.zeros_like(d),
                                      where=ok0[..., None])
            H[:, 1:][ok0] = sh[ok0]
        valid = np.isfinite(H).all(axis=2)
        dON = np.linalg.norm(N[:, :, None, :] - O[:, None, :, :], axis=3)
        dCH = np.linalg.norm(C[:, None, :, :] - H[:, :, None, :], axis=3)
        dOH = np.linalg.norm(H[:, :, None, :] - O[:, None, :, :], axis=3)
        dCN = np.linalg.norm(N[:, :, None, :] - C[:, None, :, :], axis=3)
        with np.errstate(divide="ignore", invalid="ignore"):
            E = 0.084 * 332.0 * (1.0 / dON + 1.0 / dCH - 1.0 / dOH - 1.0 / dCN)
        E = np.maximum(np.nan_to_num(E, nan=0.0, posinf=0.0, neginf=geo.HB_E_FLOOR),
                       geo.HB_E_FLOOR)
        ok = valid[:, :, None] & (self.idx_sep[None] >= 2)
        ok &= (dON > geo.HB_MIN_ON) & (dCH > 0.5) & (dOH > 0.5) & (dCN > 0.5)
        ok &= np.isfinite(E) & (E < -0.5)

        # greedy: repeatedly take the lowest still-admissible pair, block its donor and
        # acceptor. `n` rounds suffice -- each round consumes one donor.
        big = 1e9
        avail = ok.copy()
        bond = np.zeros_like(ok)
        Ef = np.where(ok, E, big)
        bi = np.arange(B)
        for _ in range(n):
            flat = np.where(avail, Ef, big).reshape(B, n * n)
            k = flat.argmin(1)
            v = flat[bi, k]
            live = v < big
            if not live.any():
                break
            d_i, a_j = k // n, k % n
            bond[bi[live], d_i[live], a_j[live]] = True
            avail[bi[live], d_i[live], :] = False
            avail[bi[live], :, a_j[live]] = False
        # desolvation cost of +1.0, matching `energy_terms.hbond_terms`
        e = np.where(bond, E + 1.0, 0.0)
        e = np.where(e < 0.0, e, 0.0)
        bond = bond & (E + 1.0 < 0.0)
        sepm = self.idx_sep[None]
        local = np.where(sepm < HB_LONGRANGE_SEP, e, 0.0).sum((1, 2))
        lr = np.where(sepm >= HB_LONGRANGE_SEP, e, 0.0).sum((1, 2))
        return local, lr, bond

    def _coop(self, bond):
        """Consecutive n-turns and consecutive antiparallel ladder rungs, batched."""
        n = self.n
        d, a = np.indices((n, n))
        hel = bond & (d - a >= 3)[None] & (d - a <= 4)[None]
        run = np.zeros(bond.shape[0])
        if n >= 2:
            nxt = np.zeros_like(hel)
            nxt[:, :-1, :-1] = hel[:, 1:, 1:]
            run = (hel & nxt).sum((1, 2)).astype(float)
        helix = -np.minimum(run, COOP_RUN_CAP)
        lad = bond & (np.abs(d - a) >= 2)[None]
        nxt = np.zeros_like(lad)
        if n >= 3:
            nxt[:, :-2, 2:] = lad[:, 2:, :-2]
        sheet = -(lad & nxt).sum((1, 2)).astype(float)
        return helix, sheet

    # ---------------------------------------------------------------- public
    def terms_from_coords(self, coords: Dict[str, np.ndarray],
                          states: Optional[np.ndarray] = None,
                          phi: Optional[np.ndarray] = None,
                          psi: Optional[np.ndarray] = None) -> np.ndarray:
        """``(B, 11)`` unweighted term values, columns ordered as `TERMS`."""
        CA = coords["CA"]
        CB = coords.get("CB", CA)
        B, n = CA.shape[0], self.n
        dcb = np.linalg.norm(CB[:, self.di, :] - CB[:, self.dj, :], axis=2)
        out = np.zeros((B, len(TERMS)))
        col = {t: k for k, t in enumerate(TERMS)}

        out[:, col["steric"]] = self._steric(coords)
        out[:, col["contact"]] = (switch(dcb[:, self.m3], 4.5, 8.5)
                                  * self.mj_pair[None, self.m3]).sum(1)
        local, lr, bond = self._hbonds(coords)
        out[:, col["hbond_local"]] = local
        out[:, col["hbond_longrange"]] = lr
        h, s = self._coop(bond)
        out[:, col["coop_helix"]] = h
        out[:, col["coop_sheet"]] = s
        # solvation: burial-weighted coordination over CB pairs
        sw = switch(dcb, 6.0, 10.0)
        coord = np.zeros((B, n))
        np.add.at(coord.T, self.di, sw.T)
        np.add.at(coord.T, self.dj, sw.T)
        out[:, col["solvation"]] = -(coord @ (self.burial / BURIAL_NORM))
        # electrostatic
        if self.elec_mask.any():
            de = np.maximum(dcb[:, self.elec_mask], 2.0)
            dr = dcb[:, self.elec_mask]
            out[:, col["electrostatic"]] = (
                (COULOMB / DIELECTRIC)
                * (self.qq[None, self.elec_mask] / de) * np.exp(-dr / 8.0)).sum(1)
        # aromatic, CB proxy
        if len(self.ar_i):
            da = np.linalg.norm(CB[:, self.ar_i, :] - CB[:, self.ar_j, :], axis=2)
            out[:, col["aromatic"]] = -np.exp(
                -((da - AROM_CB_DIST) / AROM_CB_WIDTH) ** 2).sum(1)
        # torsion
        if states is not None:
            S = np.asarray(states, int)
            out[:, col["torsion"]] = self.rama[self.rows[None, :], S].sum(1)
        elif phi is not None:
            out[:, col["torsion"]] = np.array(
                [torsion_term(self.seq, phi[b], psi[b]) for b in range(B)])
        # compactness
        rg = np.sqrt(((CA - CA.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        out[:, col["compactness"]] = np.maximum(0.0, rg - self.rg_target) ** 2
        return out

    def terms(self, states: np.ndarray) -> np.ndarray:
        S = np.asarray(states, int)
        phi = self.rep._phi[self.rows[None, :], S]
        psi = self.rep._psi[self.rows[None, :], S]
        coords = geo.build_backbone_batch(phi, psi)
        return self.terms_from_coords(coords, states=S)

# Weighting
#: Per-term weights for the *generation* field. Fitted on real candidate pools of the
#: 24 development peptides (`work/pools`, dev split), never on a benchmark target and
#: never on the synthetic decoy bank. See `work/fit_legacy.py`.
FITTED_WEIGHTS: Dict[str, float] = dict(DEFAULT_WEIGHTS)


class LegacyField:
    """A single scalar per structure: standardised, weighted Legacy energy.

    `standardise` draws `n_ref` random state assignments from the target's own library and
    records each term's mean and spread on that population. The field returned is then
    dimensionless and comparable across targets and across chain lengths, which is what
    lets one blend weight `w_legacy` transfer. Nothing in the reference sample depends on
    the native structure.
    """

    def __init__(self, sequence: str, representation,
                 weights: Optional[Dict[str, float]] = None,
                 n_ref: int = 512, seed: int = 0,
                 use: Optional[Sequence[str]] = None, dist=None):
        self.bl = BatchLegacy(sequence, representation)
        w = dict(FITTED_WEIGHTS if weights is None else weights)
        self.w = np.array([w.get(t, 0.0) for t in TERMS], float)
        if use is not None:
            keep = set(use)
            self.w = np.array([self.w[k] if t in keep else 0.0
                               for k, t in enumerate(TERMS)])
        self.mu = np.zeros(len(TERMS))
        self.sd = np.ones(len(TERMS))
        self._standardise(representation, n_ref, seed, dist)

    def _standardise(self, rep, n_ref: int, seed: int, dist) -> None:
        rng = np.random.default_rng(seed)
        S = rng.integers(0, rep.n_states, size=(n_ref, rep.n_residues))
        rows = np.arange(rep.n_residues)
        coords = geo.build_backbone_batch(rep._phi[rows[None], S],
                                          rep._psi[rows[None], S])
        T = self.bl.terms_from_coords(coords, states=S)
        self.mu = T.mean(0)
        self.sd = np.maximum(T.std(0), 1e-6)
        # the weighted sum's own scale on the same population, so `w_legacy` is a
        # multiple of "one standard deviation of the random-structure population"
        z = (T - self.mu) / self.sd
        v = z @ self.w
        self.centre = float(v.mean())
        self.scale = float(max(v.std(), 1e-6))
        # ... and, when the distance term is supplied, expressed in ITS units on the same
        # sample. Without this `w_legacy` would mean something different on every target,
        # because the Bayes-risk score's spread depends on how confident the prior is.
        self.unit = 1.0
        if dist is not None:
            d = np.asarray(dist.score(coords["CA"]), float)
            self.unit = float(max(d.std(), 1e-9))

    def _reduce(self, T: np.ndarray) -> np.ndarray:
        z = (T - self.mu[None, :]) / self.sd[None, :]
        return self.unit * (z @ self.w - self.centre) / self.scale

    def score(self, states: np.ndarray) -> np.ndarray:
        return self._reduce(self.bl.terms(states))

    def score_from_coords(self, coords, states=None, phi=None, psi=None) -> np.ndarray:
        return self._reduce(self.bl.terms_from_coords(coords, states=states,
                                                      phi=phi, psi=psi))


def verify(sequence: str, representation, batch: int = 24, seed: int = 0):
    """Max per-term deviation of `BatchLegacy` from `energy_terms.energy_components`."""
    bl = BatchLegacy(sequence, representation)
    rng = np.random.default_rng(seed)
    S = rng.integers(0, representation.n_states, size=(batch, len(sequence)))
    T = bl.terms(S)
    rows = np.arange(len(sequence))
    worst = {t: 0.0 for t in TERMS}
    for k in range(batch):
        phi = representation._phi[rows, S[k]]
        psi = representation._psi[rows, S[k]]
        coords = geo.build_backbone(phi, psi)
        ref = energy_components(sequence, coords, phi, psi)
        for c, t in enumerate(TERMS):
            worst[t] = max(worst[t], abs(ref[t] - T[k, c]))
    return worst

# SECTION 3 -- per-residue torsion libraries (was torsion_lib2.py)
_CLASSES = (reps.CLASS_GENERAL, reps.CLASS_GLY, reps.CLASS_PRO,
            reps.CLASS_PRE_PRO)
MIN_OBS = 40
#: Weight of the backed-off pool relative to the specific one. Not zero: with 60
#: observations, k=16 clusters would otherwise be fitted to ~4 points each.
BACKOFF_W = 0.25


def _circ_kmeans(X4: np.ndarray, w: np.ndarray, k: int, seed: int = 1,
                 iters: int = 60) -> np.ndarray:
    """Weighted k-means on the 4-D (cos phi, sin phi, cos psi, sin psi) embedding."""
    m = len(X4)
    k = min(k, m)
    rng = np.random.default_rng(seed)
    # k-means++ on the weights, so rare-but-real basins are not missed by a uniform draw.
    idx = [int(rng.choice(m, p=w / w.sum()))]
    d2 = ((X4 - X4[idx[0]]) ** 2).sum(1)
    for _ in range(k - 1):
        p = d2 * w
        s = p.sum()
        idx.append(int(rng.choice(m, p=p / s)) if s > 0 else int(rng.integers(m)))
        d2 = np.minimum(d2, ((X4 - X4[idx[-1]]) ** 2).sum(1))
    C = X4[idx].copy()
    for _ in range(iters):
        lab = ((X4[:, None, :] - C[None]) ** 2).sum(-1).argmin(1)
        newC = C.copy()
        for j in range(k):
            sel = lab == j
            if sel.any():
                ww = w[sel][:, None]
                newC[j] = (X4[sel] * ww).sum(0) / ww.sum()
        if np.allclose(newC, C):
            break
        C = newC
    return C


def _angles_from_centres(C: np.ndarray) -> np.ndarray:
    return np.column_stack([np.arctan2(C[:, 1], C[:, 0]),
                            np.arctan2(C[:, 3], C[:, 2])])


class ContextLibrary:
    """Observation pools keyed by sequence context, built once per holdout database."""

    def __init__(self, entries: Sequence[pdb.Peptide]):
        self.pools: Dict[tuple, List[Tuple[float, float]]] = {}
        for p in entries:
            cls = reps.residue_classes(p.seq, p.n)
            for i in range(p.n):
                a = p.seq[i]
                l = p.seq[i - 1] if i > 0 else "^"
                r = p.seq[i + 1] if i + 1 < p.n else "$"
                t = (float(p.phi[i]), float(p.psi[i]))
                for key in ((l, a, r), (a, r), (l, a), (a,), (cls[i],)):
                    self.pools.setdefault(key, []).append(t)
        self._cache: Dict[tuple, np.ndarray] = {}

    def _keys(self, seq: str, i: int, cls: str) -> List[tuple]:
        a = seq[i]
        l = seq[i - 1] if i > 0 else "^"
        r = seq[i + 1] if i + 1 < len(seq) else "$"
        return [(l, a, r), (a, r), (l, a), (a,), (cls,)]

    def states(self, seq: str, i: int, cls: str, k: int, seed: int = 1) -> np.ndarray:
        """(k, 2) array of (phi, psi) in radians for residue `i` of `seq`."""
        key = (seq[max(0, i - 1):i + 2], seq[i], cls, k, i == 0, i == len(seq) - 1)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        chain = self._keys(seq, i, cls)
        pts: List[Tuple[float, float]] = []
        wts: List[float] = []
        used = False
        for depth, kk in enumerate(chain):
            pool = self.pools.get(kk)
            if not pool:
                continue
            if not used and len(pool) >= MIN_OBS:
                pts.extend(pool)
                wts.extend([1.0] * len(pool))
                used = True
            elif used:
                pts.extend(pool)
                wts.extend([BACKOFF_W] * len(pool))
                break
        if not used:                      # fall back to the class pool outright
            pool = self.pools.get((cls,)) or self.pools.get(("GENERAL",)) or []
            pts, wts = list(pool), [1.0] * len(pool)
        P = np.asarray(pts, float)
        w = np.asarray(wts, float)
        X4 = np.column_stack([np.cos(P[:, 0]), np.sin(P[:, 0]),
                              np.cos(P[:, 1]), np.sin(P[:, 1])])
        ang = _angles_from_centres(_circ_kmeans(X4, w, k, seed=seed))
        if len(ang) < k:                  # pad by repeating (degenerate tiny pool)
            ang = np.vstack([ang, np.repeat(ang[:1], k - len(ang), axis=0)])
        self._cache[key] = ang
        return ang


@lru_cache(maxsize=32)
def _class_library(k: int, exclude_seq: str, seed: int) -> tuple:
    """k states per residue CLASS, clustered over the held-out database."""
    entries = pdb.holdout(exclude_seq) if exclude_seq else pdb.load()
    lib = ContextLibrary(entries)
    out = {}
    for c in _CLASSES:
        pool = np.asarray(lib.pools.get((c,))
                          or lib.pools[(reps.CLASS_GENERAL,)], float)
        X4 = np.column_stack([np.cos(pool[:, 0]), np.sin(pool[:, 0]),
                              np.cos(pool[:, 1]), np.sin(pool[:, 1])])
        out[c] = _angles_from_centres(
            _circ_kmeans(X4, np.ones(len(X4)), k, seed=seed))
    return tuple(sorted(out.items()))


@lru_cache(maxsize=32)
def library_for(sequence: str, k: int, exclude_seq: str = "", seed: int = 1,
                mode: str = "class") -> np.ndarray:
    """``(n, k, 2)`` per-residue (phi, psi) table in radians, holding out `exclude_seq`.

    ``mode="class"`` gives every residue the k states of its class (GENERAL/GLY/PRO/
    PRE_PRO); ``mode="ctx"`` gives each residue its own states from the back-off context
    chain above. The context form was measured against the class form on 35 targets at
    k = 4, 8, 16, 32 and is NOT better -- mean floor 1.39/1.02/0.76/0.58 A against
    1.38/0.93/0.69/0.54. Conditioning on the triplet fragments the observation pool faster
    than it sharpens the distribution, so the class form is the default and the context
    form is kept because the comparison is the evidence for that choice.
    """
    cls = reps.residue_classes(sequence, len(sequence))
    if mode == "class":
        lib = dict(_class_library(k, exclude_seq, seed))
        return np.stack([lib[c] for c in cls])
    entries = pdb.holdout(exclude_seq) if exclude_seq else pdb.load()
    ctx = ContextLibrary(entries)
    return np.stack([ctx.states(sequence, i, cls[i], k, seed=seed)
                     for i in range(len(sequence))])


class PerResidueTorsion(reps.TorsionStateRepresentation):
    """`TorsionStateRepresentation` with an explicit per-residue state table.

    Everything downstream -- bitstring layout, chi1 bits, `build_all`, the Hamiltonians,
    `batch_features` -- reads `_phi` / `_psi`, so supplying the table is all it takes.
    """

    name = "torsion_ctx"

    def __init__(self, sequence: str, table: np.ndarray, chi_bits: bool = True):
        n, k, _ = table.shape
        if n != len(sequence):
            raise ValueError(f"table has {n} residues, sequence has {len(sequence)}")
        reps.STATE_LIBRARIES.setdefault(
            k, {c: [(-63.0, -42.0)] * k for c in
                (reps.CLASS_GENERAL, reps.CLASS_GLY, reps.CLASS_PRO,
                 reps.CLASS_PRE_PRO)})
        super().__init__(n, n_states=k, sequence=sequence, chi_bits=chi_bits)
        self._phi = np.ascontiguousarray(table[:, :, 0])
        self._psi = np.ascontiguousarray(table[:, :, 1])
        self.libraries = [[(float(np.degrees(a)), float(np.degrees(b)))
                           for a, b in table[i]] for i in range(n)]


def make(sequence: str, k: int = 8, exclude_seq: Optional[str] = None,
         chi_bits: bool = False, seed: int = 1,
         mode: str = "class") -> PerResidueTorsion:
    tab = library_for(sequence, k, exclude_seq or "", seed, mode)
    return PerResidueTorsion(sequence, tab, chi_bits=chi_bits)

# SECTION 4 -- sequence priors (was priors.py)
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --------------------------------------------------------------------- residue features
#: Chou-Fasman helix/sheet propensity, Kyte-Doolittle hydropathy, formal charge at pH 7,
#: and side-chain volume (A^3). Five numbers per residue type -- enough to let the model
#: generalise across amino acids instead of memorising 20 identities from 780 peptides.
_PROPS = {
    #      helix sheet  hydro  charge  volume
    "A": (1.42, 0.83,  1.8,  0.0,  88.6), "C": (0.70, 1.19,  2.5,  0.0, 108.5),
    "D": (1.01, 0.54, -3.5, -1.0, 111.1), "E": (1.51, 0.37, -3.5, -1.0, 138.4),
    "F": (1.13, 1.38,  2.8,  0.0, 189.9), "G": (0.57, 0.75, -0.4,  0.0,  60.1),
    "H": (1.00, 0.87, -3.2,  0.1, 153.2), "I": (1.08, 1.60,  4.5,  0.0, 166.7),
    "K": (1.16, 0.74, -3.9,  1.0, 168.6), "L": (1.21, 1.30,  3.8,  0.0, 166.7),
    "M": (1.45, 1.05,  1.9,  0.0, 162.9), "N": (0.67, 0.89, -3.5,  0.0, 114.1),
    "P": (0.57, 0.55, -1.6,  0.0, 112.7), "Q": (1.11, 1.10, -3.5,  0.0, 143.8),
    "R": (0.98, 0.93, -4.5,  1.0, 173.4), "S": (0.77, 0.75, -0.8,  0.0,  89.0),
    "T": (0.83, 1.19, -0.7,  0.0, 116.1), "V": (1.06, 1.70,  4.2,  0.0, 140.0),
    "W": (1.08, 1.37, -0.9,  0.0, 227.8), "Y": (0.69, 1.47, -1.3,  0.0, 193.6),
}
_DEFAULT = tuple(float(np.mean([v[k] for v in _PROPS.values()])) for k in range(5))


def _prop(seq: str) -> np.ndarray:
    return np.array([_PROPS.get(a, _DEFAULT) for a in seq], float)

#: Distance bin edges, A. Fine where CA-CA distances are structurally informative
#: (5-14 A spans an i,i+3 helical turn through a hairpin cross-strand pair) and coarse
#: beyond, where "far apart" is all the information there is.
BIN_EDGES = np.array([4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0,
                      11.0, 12.5, 14.0, 16.0, 19.0, 23.0])
NBINS = len(BIN_EDGES) + 1


def pair_features(seq: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(features, i, j)`` for every pair with ``j - i >= 2`` of one sequence."""
    n = len(seq)
    P = _prop(seq)
    # Local windows: mean propensity over +-2 residues, the classic secondary-structure
    # signal, computed by stacking shifted copies so it costs nothing.
    pad = np.vstack([P[:1], P[:1], P, P[-1:], P[-1:]])
    W = np.stack([pad[k:k + n] for k in range(5)]).mean(0)
    i, j = np.triu_indices(n, k=2)
    s = (j - i).astype(float)
    between = np.stack([P[a + 1:b].mean(0) if b > a + 1 else P[a]
                        for a, b in zip(i, j)])
    f = np.column_stack([
        s, np.log(s), s / n, np.full(i.shape, float(n)),
        i / n, j / n, np.minimum(i, n - 1 - j).astype(float),
        P[i], P[j], W[i], W[j],
        P[i] * P[j], np.abs(P[i] - P[j]),
        between,
    ])
    return f.astype(np.float32), i, j


def fit_distogram_model(entries: Sequence[pdb.Peptide], seed: int = 0,
                        max_iter: int = 200):
    from sklearn.ensemble import HistGradientBoostingClassifier
    X, Y = [], []
    for p in entries:
        f, i, j = pair_features(p.seq)
        if not len(f):
            continue
        d = np.linalg.norm(p.ca[i] - p.ca[j], axis=1)
        X.append(f)
        Y.append(np.digitize(d, BIN_EDGES))
    X = np.vstack(X)
    Y = np.concatenate(Y)
    m = HistGradientBoostingClassifier(max_iter=max_iter, learning_rate=0.1, max_depth=7,
                                       l2_regularization=1.0, random_state=seed,
                                       early_stopping=False)
    m.fit(X, Y)
    return m


class DistogramPrior:
    """Per-target log P(distance bin) table, shape ``(n, n, NBINS)``."""

    def __init__(self, sequence: str, logp: np.ndarray):
        self.seq = sequence
        self.n = len(sequence)
        self.logp = logp
        self._iu = np.triu_indices(self.n, k=2)
        self._flat = logp[self._iu[0], self._iu[1]]          # (npairs, NBINS)

    @classmethod
    def fit(cls, sequence: str, entries: Optional[Sequence[pdb.Peptide]] = None,
            model=None, seed: int = 0) -> "DistogramPrior":
        if model is None:
            model = fit_distogram_model(
                entries if entries is not None else pdb.holdout(sequence), seed)
        n = len(sequence)
        f, i, j = pair_features(sequence)
        prob = model.predict_proba(f)
        cols = np.asarray(model.classes_, int)
        tab = np.full((len(f), NBINS), 1e-4)
        tab[:, cols] = np.maximum(prob, 1e-4)
        tab /= tab.sum(1, keepdims=True)
        full = np.full((n, n, NBINS), 1.0 / NBINS)
        full[i, j] = tab
        full[j, i] = tab
        return cls(sequence, np.log(full).astype(np.float32))

    def score(self, ca: np.ndarray) -> np.ndarray:
        """Mean per-pair negative log-likelihood for a batch ``(B, n, 3)`` of CA traces."""
        arr = np.asarray(ca, float)
        single = arr.ndim == 2
        if single:
            arr = arr[None]
        i, j = self._iu
        d = np.linalg.norm(arr[:, i, :] - arr[:, j, :], axis=-1)
        b = np.searchsorted(BIN_EDGES, d)
        lp = np.take_along_axis(self._flat[None], b[:, :, None], axis=2)[:, :, 0]
        out = -lp.mean(1)
        return float(out[0]) if single else out

    def expected_distances(self) -> np.ndarray:
        """Posterior-mean CA-CA distance per pair, for diagnostics and restraints."""
        centres = np.concatenate([[BIN_EDGES[0] - 0.5],
                                  0.5 * (BIN_EDGES[1:] + BIN_EDGES[:-1]),
                                  [BIN_EDGES[-1] + 2.0]])
        return (np.exp(self.logp) * centres).sum(-1)


# --------------------------------------------------------------------- torsion MRF
class TorsionMRF:
    """Per-residue state marginals plus nearest-neighbour state couplings.

    Estimated by projecting held-out database residues onto the target's per-residue
    library: for slot ``i``, every database residue of the same type at a comparable
    fractional chain position contributes its nearest state. Couplings use consecutive
    database pairs whose residue types match ``(seq[i], seq[i+1])``.
    """

    def __init__(self, h: np.ndarray, J: np.ndarray):
        self.h = h                # (n, k)      -log marginal
        self.J = J                # (n-1, k, k) -log (joint / marginal product)

    @classmethod
    def fit(cls, sequence: str, table: np.ndarray,
            entries: Sequence[pdb.Peptide], smooth: float = 1.0,
            pos_tol: float = 0.34) -> "TorsionMRF":
        n, k, _ = table.shape
        counts = np.full((n, k), smooth)
        pair = np.full((n - 1, k, k), smooth / k)
        cphi, sphi = np.cos(table[:, :, 0]), np.sin(table[:, :, 0])
        cpsi, spsi = np.cos(table[:, :, 1]), np.sin(table[:, :, 1])

        def nearest(i: int, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
            d = ((np.cos(phi)[:, None] - cphi[i]) ** 2
                 + (np.sin(phi)[:, None] - sphi[i]) ** 2
                 + (np.cos(psi)[:, None] - cpsi[i]) ** 2
                 + (np.sin(psi)[:, None] - spsi[i]) ** 2)
            return d.argmin(1)

        by_aa: Dict[str, List[Tuple[float, float, float, str, float, float]]] = {}
        for p in entries:
            for t in range(p.n):
                by_aa.setdefault(p.seq[t], []).append(
                    (t / max(1, p.n - 1), float(p.phi[t]), float(p.psi[t]),
                     p.seq[t + 1] if t + 1 < p.n else "",
                     float(p.phi[t + 1]) if t + 1 < p.n else np.nan,
                     float(p.psi[t + 1]) if t + 1 < p.n else np.nan))

        for i in range(n):
            frac = i / max(1, n - 1)
            rows = [r for r in by_aa.get(sequence[i], ())
                    if abs(r[0] - frac) <= pos_tol]
            if rows:
                s = nearest(i, np.array([r[1] for r in rows]),
                            np.array([r[2] for r in rows]))
                np.add.at(counts[i], s, 1.0)
            if i < n - 1:
                nxt = [r for r in rows if r[3] == sequence[i + 1]]
                if nxt:
                    s0 = nearest(i, np.array([r[1] for r in nxt]),
                                 np.array([r[2] for r in nxt]))
                    s1 = nearest(i + 1, np.array([r[4] for r in nxt]),
                                 np.array([r[5] for r in nxt]))
                    np.add.at(pair[i], (s0, s1), 1.0)

        marg = counts / counts.sum(1, keepdims=True)
        h = -np.log(marg)
        J = np.zeros((max(0, n - 1), k, k))
        for i in range(n - 1):
            pj = pair[i] / pair[i].sum()
            J[i] = (-np.log(np.maximum(pj, 1e-12))
                    + np.log(marg[i])[:, None] + np.log(marg[i + 1])[None, :])
        return cls(h.astype(np.float32), J.astype(np.float32))

    def score(self, states: np.ndarray) -> np.ndarray:
        """``states`` (B, n) of state indices -> (B,) mean per-residue MRF energy."""
        S = np.atleast_2d(np.asarray(states, int))
        n = S.shape[1]
        e = self.h[np.arange(n)[None], S].sum(1)
        if n > 1 and len(self.J):
            e = e + self.J[np.arange(n - 1)[None], S[:, :-1], S[:, 1:]].sum(1)
        return e / n

    def marginals(self) -> np.ndarray:
        p = np.exp(-self.h)
        return p / p.sum(1, keepdims=True)

# SECTION 5 -- batched scalar-exact scoring of a whole candidate pool
"""`components_batch` is `energy_components` over B structures at once, exactly.

The design rule here was found by measuring, and it is the opposite of the obvious one.

The obvious version batches the *reductions* -- one `einsum` for steric over the whole
pool, one matrix-vector product for contact, `np.add.at` for solvation. That was built
first. It ran at **1.05x** and it was **not bit-exact**: reassociating the sums moved
steric by 2.8e-17, contact by 8.9e-16 and solvation by 7.1e-15. Both halves of that are
bad. The speedup was absent because these terms are not where the time is, and the
exactness was spent for nothing.

Where the time actually is, at n=14, is numpy *dispatch* on tiny arrays -- and almost all
of it is geometry, not reduction:

    hbond_terms   0.106 ms   of which dssp_energy_matrix builds four (n, n) distance
                             matrices from four (n, n, 3) temporaries
    steric_term   0.073 ms   one gather and one norm over ~2400 pairs
    everything else, all eleven reductions together, ~0.086 ms

So this batches the geometry, where one dispatch replaces B, and then applies the
**unmodified scalar reductions** per structure -- `np.dot` for steric and contact,
`np.bincount` for solvation, the same greedy one-donor-one-acceptor loop for H-bonds.
Those per-structure calls are BLAS on vectors that are already materialised; they cost
almost nothing, and they are the reason every term below agrees with `energy_components`
to exactly 0.0 rather than to 1e-15.

The greedy H-bond match stays sequential on purpose. Each accepted bond removes a donor
and an acceptor from contention, so vectorising it means changing the tie-break -- which
is what `BatchLegacy._hbonds` does, and why that is documented as a different model for
generation rather than an optimisation of this one.
"""


class _PoolTables:
    """Everything about a sequence and a chain length that does not depend on coordinates.

    The scalar path recomputes (or re-looks-up) all of this on every call. Held here it is
    built once per pool. All read-only.
    """

    __slots__ = ("seq", "n", "burial", "q", "mj", "di", "dj", "sep", "mj_m3",
                 "m3", "elec_mask", "qq_m", "st_names", "st_ii", "st_jj",
                 "st_lim", "rg_target", "burial_norm")

    def __init__(self, sequence: str, coords_keys, use_corrected_mj: bool = True):
        self.seq = sequence
        self.n = n = len(sequence)
        self.burial, self.q, self.mj = sequence_arrays(sequence, use_corrected_mj)
        self.burial_norm = self.burial / BURIAL_NORM
        di, dj, sep = pair_index(n)
        self.di, self.dj, self.sep = di, dj, sep
        self.m3 = sep >= 3
        self.mj_m3 = self.mj[di[self.m3], dj[self.m3]]
        qq = self.q[di] * self.q[dj]
        self.elec_mask = (sep >= 2) & (qq != 0.0)
        self.qq_m = qq[self.elec_mask]
        self.st_names = tuple(k for k in ("N", "CA", "C", "O", "CB")
                              if k in coords_keys)
        self.st_ii, self.st_jj, self.st_lim = _steric_layout(
            self.st_names, n, (), 2)
        self.rg_target = 2.2 * (n ** 0.38)

_POOL_TABLES: Dict[tuple, _PoolTables] = {}


def _tables(sequence: str, coords_keys, use_corrected_mj: bool = True) -> _PoolTables:
    key = (sequence, tuple(sorted(coords_keys)), bool(use_corrected_mj))
    hit = _POOL_TABLES.get(key)
    if hit is None:
        hit = _PoolTables(sequence, coords_keys, use_corrected_mj)
        if len(_POOL_TABLES) > 64:
            _POOL_TABLES.clear()
        _POOL_TABLES[key] = hit
    return hit


def clear_pool_cache() -> None:
    _POOL_TABLES.clear()


def components_batch(sequence: str,
                     coords: Dict[str, np.ndarray],
                     phi: Optional[np.ndarray] = None,
                     psi: Optional[np.ndarray] = None,
                     use_corrected_mj: bool = True,
                     rings: Optional[Sequence[Optional[Dict]]] = None
                     ) -> Dict[str, np.ndarray]:
    """`energy_components` for a whole pool. Same numbers, one call.

    Parameters
    coords : dict of (B, n, 3) arrays in ANGSTROMS -- e.g. straight out of
             `protein_geometry.build_backbone_batch`. Must carry at least "CA";
             "N", "C", "O" enable the H-bond terms and "CB" the CB-pair terms,
             exactly as in the scalar path.
    phi, psi : (B, n) or None. None leaves `torsion` at 0.0, as the scalar path does.
    rings  : optional per-structure ring dicts, one per candidate, for the real
             orientation-dependent aromatic term.

    Returns
    dict of term name -> (B,) float array, keys and order exactly `TERM_NAMES`.
    """
    CA = np.asarray(coords["CA"], dtype=float)
    if CA.ndim != 3:
        raise ValueError(f"components_batch wants (B, n, 3) coordinates, got {CA.shape}")
    B, n = CA.shape[0], CA.shape[1]
    if n != len(sequence):
        raise ValueError(f"coordinate count {n} != sequence length {len(sequence)}")
    T = _tables(sequence, coords.keys(), use_corrected_mj)
    CB = np.asarray(coords.get("CB", CA), dtype=float)
    out = {k: np.zeros(B) for k in TERM_NAMES}

    # ================= batched geometry: one dispatch instead of B ==========
    # CB-CB pair distances. The scalar path recomputes these per structure, and again
    # inside `steric_term`'s superset, and again per aromatic pair.
    d_cb = np.linalg.norm(CB[:, T.di, :] - CB[:, T.dj, :], axis=2)     # (B, P)

    over = None
    if T.st_names and T.st_ii.size:
        A = np.concatenate([np.asarray(coords[k], float) for k in T.st_names], axis=1)
        over = np.maximum(0.0, T.st_lim[None, :]
                          - np.linalg.norm(A[:, T.st_ii] - A[:, T.st_jj], axis=2))

    have_bb = all(k in coords for k in ("N", "C", "O"))
    Eb = okb = None
    if have_bb:
        Eb, okb = _dssp_matrices_batch(np.asarray(coords["N"], float),
                                       np.asarray(coords["C"], float),
                                       np.asarray(coords["O"], float))

    sw_contact = switch(d_cb[:, T.m3], 4.5, 8.5)
    sw_solv = switch(d_cb, 6.0, 10.0)
    d_elec = (np.maximum(d_cb[:, T.elec_mask], 2.0), d_cb[:, T.elec_mask]) \
        if T.elec_mask.any() else None
    ca_c = CA - CA.mean(axis=1, keepdims=True)

    # ================= exact scalar reductions, per structure ===============
    # Every line below is the reduction the scalar term function performs, on data that
    # is already materialised. Same primitive, same order, same float.
    aromatic_here = len(aromatic_indices(sequence)) >= 2
    for b in range(B):
        # `.copy()` is load-bearing, not defensive. A row of a (B, M) block and a freshly
        # allocated M-vector hold the same floats, but BLAS `ddot` picks its SIMD path on
        # the buffer's alignment, and a row at an odd offset reduces in a different order.
        # Measured on a real 1A13 pool: 37 of 100 contact energies moved by up to 4.4e-16
        # without this. The copy costs one memcpy of M doubles and buys exact agreement
        # with `contact_term`.
        if over is not None:
            o = over[b].copy()
            out["steric"][b] = float(o @ o)
        out["contact"][b] = float(np.dot(T.mj_m3, sw_contact[b].copy()))
        coord = (np.bincount(T.di, weights=sw_solv[b], minlength=n)
                 + np.bincount(T.dj, weights=sw_solv[b], minlength=n))
        out["solvation"][b] = float(-T.burial_norm @ coord)
        if d_elec is not None:
            out["electrostatic"][b] = float(np.sum(
                (COULOMB / DIELECTRIC) * (T.qq_m / d_elec[0][b])
                * np.exp(-d_elec[1][b] / 8.0)))
        out["compactness"][b] = float(
            max(0.0, float(np.sqrt(np.mean(np.sum(ca_c[b] ** 2, axis=1))))
                - T.rg_target) ** 2)
        if have_bb:
            hl, hr, pairs = _hbond_match(Eb[b], okb[b], n)
            out["hbond_local"][b] = hl
            out["hbond_longrange"][b] = hr
            out["coop_helix"][b] = coop_helix_term(pairs)
            out["coop_sheet"][b] = coop_sheet_term(pairs)
        rb = rings[b] if rings is not None else None
        if aromatic_here:
            out["aromatic"][b] = aromatic_term(
                sequence, {"CA": CA[b], "CB": CB[b]}, rings=rb)
        if phi is not None:
            out["torsion"][b] = torsion_term(sequence, phi[b], psi[b])
    return out


def _dssp_matrices_batch(N: np.ndarray, C: np.ndarray, O: np.ndarray):
    """`protein_geometry.dssp_energy_matrix` over a whole pool.

    The arithmetic is line-for-line the scalar function's, with one leading batch axis.
    Every reduction here is `np.linalg.norm` over the last axis of three contiguous
    floats, which is the same three-element `add.reduce` the scalar path performs, so the
    distances -- and therefore E and ok -- come out bit-identical. `tests/test_energy.py`
    checks that against `protein_geometry.dssp_energy_matrix` directly rather than
    trusting the argument.
    """
    B, n = N.shape[0], N.shape[1]
    H = np.full((B, n, 3), np.nan)
    if n >= 2:
        d = C[:, :-1] - O[:, :-1]
        nd = np.linalg.norm(d, axis=2)
        ok0 = nd > 1e-6
        # `H[1:][ok] = N[1:][ok] + d[ok]/nd[ok][:, None]` with a batch axis in front.
        shifted = N[:, 1:] + np.divide(d, nd[..., None],
                                       out=np.zeros_like(d), where=ok0[..., None])
        tgt = H[:, 1:]
        tgt[ok0] = shifted[ok0]
        H[:, 1:] = tgt
    valid = np.isfinite(H).all(axis=2)

    dON = np.linalg.norm(N[:, :, None, :] - O[:, None, :, :], axis=3)
    dCH = np.linalg.norm(C[:, None, :, :] - H[:, :, None, :], axis=3)
    dOH = np.linalg.norm(H[:, :, None, :] - O[:, None, :, :], axis=3)
    dCN = np.linalg.norm(N[:, :, None, :] - C[:, None, :, :], axis=3)

    with np.errstate(divide="ignore", invalid="ignore"):
        E = 0.084 * 332.0 * (1.0 / dON + 1.0 / dCH - 1.0 / dOH - 1.0 / dCN)
    E = np.maximum(E, geo.HB_E_FLOOR)

    idx = np.arange(n)
    sep = np.abs(idx[:, None] - idx[None, :])
    ok = valid[:, :, None] & (sep >= 2)[None]
    ok &= (dON > geo.HB_MIN_ON) & (dCH > 0.5) & (dOH > 0.5) & (dCN > 0.5)
    ok &= np.isfinite(E) & (E < -0.5)
    return E, ok


def _hbond_match(E: np.ndarray, ok: np.ndarray, n: int,
                 desolvation_cost: float = 1.0):
    """The greedy one-donor-one-acceptor match, lifted verbatim from `hbond_terms`.

    Split out only so it can run on a precomputed (E, ok) pair. The stable sort matters:
    `HB_E_FLOOR` manufactures exact ties at -4.0 and this match is order-dependent, so an
    unstable sort could map one structure to two energies across platforms.
    """
    donors, acceptors = np.where(ok)
    if len(donors) == 0:
        return 0.0, 0.0, ()
    energies = E[donors, acceptors]
    order = np.argsort(energies, kind="stable")
    donor_used = np.zeros(n, dtype=bool)
    acc_used = np.zeros(n, dtype=bool)
    local = lr = 0.0
    matched = []
    for k in order:
        i, j = int(donors[k]), int(acceptors[k])
        if donor_used[i] or acc_used[j]:
            continue
        donor_used[i] = True
        acc_used[j] = True
        e = float(energies[k]) + desolvation_cost
        if e >= 0.0:
            continue
        matched.append((i, j))
        if abs(i - j) < HB_LONGRANGE_SEP:
            local += e
        else:
            lr += e
    return local, lr, tuple(matched)


def totals_batch(components: Dict[str, np.ndarray],
                 weights: Optional[Dict[str, float]] = None) -> np.ndarray:
    """Weighted Legacy energy per structure, from `components_batch` output.

    Defaults to `FITTED_WEIGHTS` -- the weighting fitted on real candidate pools rather
    than on a synthetic decoy bank, which is the combination the shipped Legacy model
    scores with.

    Accumulating this as eleven vectorised `tot += w * term` steps looks equivalent and is
    not: with numpy scalars in the sum, 45 of 100 totals on a real 1A13 pool came out up
    to 3.6e-15 from `total_from_components` on term values that were themselves bit-
    identical. So this calls the shipped function per structure on Python floats. It is
    B x 11 scalar operations -- microseconds against the geometry above -- and it is exact
    by construction rather than by argument.
    """
    w = FITTED_WEIGHTS if weights is None else weights
    B = len(next(iter(components.values())))
    return np.array([
        total_from_components({k: float(components[k][b]) for k in TERM_NAMES}, w)
        for b in range(B)])
