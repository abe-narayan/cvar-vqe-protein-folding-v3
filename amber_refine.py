"""Restrained all-atom AMBER refinement of a discrete-torsion candidate.

Why this module exists
`torsion_lib2` builds structures from a *discrete* library of (phi, psi) states with ideal
bond lengths and angles and a fixed rotamer per sidechain. Such a structure is
stereochemically crude: after `floor.floor()` picks the states closest to a native, the
resulting all-atom model still carries tens to hundreds of kcal/mol of bond/angle/torsion
strain (measured: 640 kcal/mol of angle strain on 6RQS) purely because ideal-geometry
NeRF placement plus frozen hydrogen frames do not satisfy ff14SB's internal coordinates.

This module fixes that strain *without letting the fold move*, by running an OpenMM
`LocalEnergyMinimizer` with a harmonic positional restraint on every backbone N, CA and C.
Everything else -- carbonyl O, OXT, all sidechain heavy atoms, all hydrogens -- is free.

The restraint choice
`k_restraint` is in **kcal/mol/A^2** per restrained atom, matching
`AmberHamiltonian.restraint_k`. The potential is 0.5*k*|r - r0|^2, so an atom pulled by an
internal force F settles at a displacement d = F/k. Typical residual internal forces in an
ideal-geometry model are 10-100 kcal/mol/A, which gives:

    k =   1  kcal/mol/A^2  ->  d ~ 0.1 - 1.0 A   WEAK. Stereochemistry is fully repaired
                                                 but the backbone is free to drift into a
                                                 different (usually more compact) basin.
    k =  10  kcal/mol/A^2  ->  d ~ 0.03 - 0.3 A  MODERATE. This is the default. Enough
                                                 give for the minimiser to relieve bond
                                                 and angle strain, far too little for a
                                                 torsion to flip -- a phi/psi flip needs
                                                 several A of CA motion.
    k = 100  kcal/mol/A^2  ->  d ~ 0.01 - 0.1 A  STRONG. The backbone frame is effectively
                                                 pinned; only sidechains and hydrogens
                                                 relax. This is what `AmberHamiltonian`
                                                 uses for its ranking energy.

The default `steps=0` means "minimise until converged" (OpenMM's convention for
`maxIterations`), with `tolerance` in kJ/mol/nm on the RMS force. Determinism comes from
the CPU platform with `Threads=1` and `DeterministicForces=true`, and from the fact that
the starting structure is a pure function of (sequence, states): no random seeding, no
thermostat, no MD. Two calls with identical arguments return bit-identical coordinates.

**Cost / accuracy setting.** Converged minimisation costs 2.4-24 s CPU per structure.
`steps=200, tolerance=5.0` costs 0.9-2.6 s CPU (3-10x less) and lands within 1-9 kcal/mol
of converged with the restrained set displaced to within 0.02 A of the same place -- i.e.
the same geometry, a slightly less-relaxed energy. Use it when refining many candidates;
use the default when the energy itself is the quantity of interest.

**Measured effect on accuracy** (12 benchmark targets, oracle `floor.floor` starts, see
`work/measure_refine.py`): refinement is RMSD-*neutral* at k=10 and k=100 (mean dCA-RMSD
-0.002 and -0.008 A) and mildly *harmful* at k=1 (+0.083 A). It repairs stereochemistry --
the built structure sits ~19,000 kcal/mol above its own minimum -- but it does not move a
structure toward native. Treat it as an energy-cleanup operator, not an accuracy operator.

Relationship to `AmberHamiltonian`
This module does *not* modify `amber_hamiltonian.py`. It composes with it: an
`AmberHamiltonian` instance is used purely as a *builder* (topology with hydrogens, ff14SB
+ GBn2 System, calibrated hydrogen local frames, OpenMM Context), and this module then
drives its own minimisation with its own restraint constant. Two of `AmberHamiltonian`'s
behaviours are deliberately switched off for refinement:

  * `collapse_floor` is set to -inf. The shipped default of -600 kcal/mol converts genuine
    energies of highly charged peptides into +inf (see `work/amber_audit.md`); a refiner
    must report the number it actually computed.
  * `minimization_steps` on the builder is irrelevant here -- `refine()` never calls
    `AmberHamiltonian.energy()`.

API
---
    refine(sequence, rep, states, k_restraint=10.0, steps=0) -> dict
    refine_coords(sequence, rep, coords, ...)   # a backbone dict, e.g. rep.build_coords()
    refine_ca(sequence, rep, ca, ...)           # CA trace; states fitted first
    single_point(sequence, rep, states, ...)    # energy with no minimisation
"""
from __future__ import annotations

import time
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import openmm
from openmm import unit

import floor as _floor
from amber_hamiltonian import AmberHamiltonian, KJ_PER_KCAL

__all__ = ["refine", "refine_coords", "refine_ca", "single_point",
           "builder_for", "clear_cache",
           "K_WEAK", "K_MODERATE", "K_STRONG", "BACKBONE_ATOMS"]

#: Restraint strengths in kcal/mol/A^2 (see module docstring).
K_WEAK = 1.0
K_MODERATE = 10.0
K_STRONG = 100.0

#: Atoms held by the positional restraint. The carbonyl O is deliberately free: it carries
#: real ideal-geometry strain and is not needed to pin the fold.
BACKBONE_ATOMS = ("N", "CA", "C")

#: kcal/mol/A^2 -> kJ/mol/nm^2
_K_SCALE = KJ_PER_KCAL * 100.0

_BUILDERS: Dict[Tuple[str, int, int, str], AmberHamiltonian] = {}


def builder_for(sequence: str, rep, platform_name: str = "CPU") -> AmberHamiltonian:
    """Cached ff14SB/GBn2 builder for one sequence + representation shape.

    Construction costs ~1-3 s (topology, hydrogen addition, and a one-off hydrogen-frame
    calibration), so it is cached; refinement of a second candidate for the same target is
    then just a minimisation.
    """
    key = (sequence.strip().upper(), int(rep.n_bits), int(rep.n_states), platform_name)
    h = _BUILDERS.get(key)
    if h is None:
        h = AmberHamiltonian(sequence, rep, restraint_k=K_MODERATE,
                             minimization_steps=0, platform_name=platform_name,
                             collapse_floor=float("-inf"))
        _BUILDERS[key] = h
    return h


def clear_cache() -> None:
    _BUILDERS.clear()


def _split(H: AmberHamiltonian, pos_nm: np.ndarray) -> Dict[str, object]:
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


def _energy_kcal(H: AmberHamiltonian) -> float:
    """Potential energy with the restraint switched off, in kcal/mol."""
    H.context.setParameter("k_rest", 0.0)
    return H.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)


def _run(H: AmberHamiltonian, heavy_nm: np.ndarray, k_restraint: float,
         steps: int, tolerance: float, components: bool) -> Dict[str, object]:
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
    if components:
        from amber_hamiltonian import AMBER_TERMS, _TERM_GROUP
        out["components"] = {
            t: ctx.getState(getEnergy=True, groups={_TERM_GROUP[t]})
                  .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            for t in AMBER_TERMS}
    out["wall"] = time.time() - t0
    return out


def refine_coords(sequence: str, rep, coords: Dict[str, np.ndarray],
                  k_restraint: float = K_MODERATE, steps: int = 0,
                  tolerance: float = 1.0, chi1: Optional[Dict[int, float]] = None,
                  platform_name: str = "CPU",
                  components: bool = False) -> Dict[str, object]:
    """Refine from a backbone dict (N, CA, C, O, CB in angstrom), e.g. `rep.build_coords`."""
    H = builder_for(sequence, rep, platform_name)
    return _run(H, H._heavy_positions(coords, chi1=chi1),
                k_restraint, steps, tolerance, components)


def refine(sequence: str, rep, states: Sequence[int],
           k_restraint: float = K_MODERATE, steps: int = 0,
           tolerance: float = 1.0, platform_name: str = "CPU",
           components: bool = False) -> Dict[str, object]:
    """Restrained all-atom ff14SB/GBn2 refinement of one discrete candidate.

    Parameters
    sequence : one-letter sequence.
    rep      : a `torsion_lib2.PerResidueTorsion` (or any full-backbone representation).
    states   : per-residue torsion state indices, or a bitstring for `rep`.
    k_restraint : harmonic positional restraint on backbone N/CA/C, kcal/mol/A^2.
                  See the module docstring; K_WEAK / K_MODERATE / K_STRONG.
    steps    : max minimiser iterations; 0 = run to `tolerance` (OpenMM convention).
    tolerance: RMS force convergence threshold, kJ/mol/nm.

    Returns
    dict with
        ca              (n, 3) refined CA coordinates, angstrom
        backbone        {"N","CA","C","O"} -> (n, 3), angstrom
        heavy           (m, 3) all heavy atoms, angstrom, ordered as `heavy_names`
        heavy_names     [(residue_index, atom_name)] for `heavy`
        positions_nm    (n_atoms, 3) full all-atom positions in nm (hydrogens included)
        energy          ff14SB + GBn2 potential energy of the refined structure,
                        kcal/mol, with the restraint term excluded
        energy_initial  the same quantity for the *input* structure
        restraint_rmsd  RMSD (A) of the restrained N/CA/C set, input -> refined
        restraint_max   largest single restrained-atom displacement, A
        wall            seconds
    """
    if isinstance(states, str):
        bits = states
    else:
        bits = rep.bitstring_from_states(np.asarray(states, dtype=int))
    H = builder_for(sequence, rep, platform_name)
    heavy = H._heavy_positions(rep.build_coords(bits), chi1=H._chi1_of(bits))
    return _run(H, heavy, k_restraint, steps, tolerance, components)


def refine_ca(sequence: str, rep, ca: np.ndarray,
              k_restraint: float = K_MODERATE, steps: int = 0,
              tolerance: float = 1.0, platform_name: str = "CPU",
              sweeps: int = 12, seed: int = 0,
              components: bool = False) -> Dict[str, object]:
    """Refine from a bare CA trace.

    A CA trace does not define an all-atom structure, so the representation's discrete
    states are first fitted to it by the same coordinate descent `floor` uses
    (`floor.descend`), and the resulting full backbone is what gets refined. The returned
    dict carries the fitted `states` and the CA-RMSD of that fit (`fit_rmsd`) so the
    projection error is never hidden inside the refinement number.
    """
    ca = np.asarray(ca, float)
    n = rep.n_residues
    rng = np.random.default_rng(seed)
    best_s, best = _floor.descend(rep, ca, rng.integers(0, rep.n_states, n),
                                  sweeps=sweeps)
    out = refine(sequence, rep, best_s, k_restraint, steps, tolerance,
                 platform_name, components)
    out["states"] = np.asarray(best_s, int)
    out["fit_rmsd"] = float(best)
    return out


def single_point(sequence: str, rep, states, platform_name: str = "CPU",
                 components: bool = False) -> Dict[str, object]:
    """ff14SB + GBn2 energy of the built structure with NO minimisation at all.

    This is the honest single-point cost/number, unlike `AmberHamiltonian.energy()`, which
    always minimises first.
    """
    return refine(sequence, rep, states, k_restraint=0.0, steps=-1,
                  tolerance=1e9, platform_name=platform_name, components=components)
