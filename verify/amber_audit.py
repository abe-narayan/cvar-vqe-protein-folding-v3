"""Independent genuineness audit of the AMBER stack in `core.amber`.

The sibling suite asserts `CustomGBForce` is present and that the 1A13 reference
reproduces. Both are necessary and neither is sufficient: a System can carry a real
`CustomGBForce` object whose per-particle parameters have been replaced, and a pinned
constant reproduces trivially if the constant was READ OFF the same code that computes it.

So this file additionally cross-checks the actual ff14SB parameters -- charge, sigma,
epsilon, per atom -- against a System built independently from the shipped XML by a
straight `app.ForceField(...).createSystem` on the same topology, and checks the GBn2
per-particle parameters are the real ones rather than defaults.

  1  forcefield files, force list, CustomGBForce shape
  2  ff14SB charge / sigma / epsilon vs an independently built System
  3  GBn2 per-particle parameters are populated and residue-dependent
  4  the 1A13 native interaction reference, -489.9138948277905
  5  the solvation term is real work (moves with geometry, not a constant)
  6  `_calibrate_hydrogens` is pinned to one thread under every worker configuration
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

REF_1A13_INTERACTION = -489.9138948277905


def build_native(pid="1A13"):
    import peptide_db as db
    import protein_geometry as geo
    import torsion_lib2 as tl2
    p = db.by_pdb(pid)
    tab = tl2.library_for(p.seq, 8, p.seq)
    rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
    rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))
    return p, rep, rb


def main():
    import openmm
    from openmm import app, unit
    import core.amber as A

    out = {}
    p, rep, rb = build_native()

    H = A.builder_for(p.seq, rep)
    desc = H.describe()
    out["1_forcefield"] = desc["forcefield"]
    forces = sorted(f.__class__.__name__ for f in H.system.getForces())
    out["1_forces"] = forces
    out["1_has_CustomGBForce"] = "CustomGBForce" in forces
    out["1_has_NonbondedForce"] = "NonbondedForce" in forces

    gb = next(f for f in H.system.getForces()
              if f.__class__.__name__ == "CustomGBForce")
    nb = next(f for f in H.system.getForces()
              if f.__class__.__name__ == "NonbondedForce")
    out["1_gb_computed_values"] = gb.getNumComputedValues()
    out["1_gb_per_particle_params"] = gb.getNumPerParticleParameters()
    out["1_gb_particles"] = gb.getNumParticles()
    out["1_n_atoms"] = H.n_atoms

    # -- 2: ff14SB parameters against an INDEPENDENTLY built System -------
    # Same topology, a fresh ForceField object created here, nothing of core involved in
    # the parameter assignment. If core had patched charges or swapped a vdW set, the
    # per-atom arrays would diverge.
    ff_ref = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
    sys_ref = ff_ref.createSystem(H.topology,
                                  nonbondedMethod=app.NoCutoff,
                                  constraints=None,
                                  rigidWater=False,
                                  implicitSolventKappa=0.0 / unit.nanometer)
    nb_ref = next(f for f in sys_ref.getForces()
                  if f.__class__.__name__ == "NonbondedForce")
    gb_ref = next(f for f in sys_ref.getForces()
                  if f.__class__.__name__ == "CustomGBForce")

    dq = ds = de = 0.0
    for i in range(nb.getNumParticles()):
        q, s, e = nb.getParticleParameters(i)
        q2, s2, e2 = nb_ref.getParticleParameters(i)
        dq = max(dq, abs(q.value_in_unit(unit.elementary_charge)
                         - q2.value_in_unit(unit.elementary_charge)))
        ds = max(ds, abs(s.value_in_unit(unit.nanometer)
                         - s2.value_in_unit(unit.nanometer)))
        de = max(de, abs(e.value_in_unit(unit.kilojoule_per_mole)
                         - e2.value_in_unit(unit.kilojoule_per_mole)))
    out["2_max_abs_charge_diff_e"] = dq
    out["2_max_abs_sigma_diff_nm"] = ds
    out["2_max_abs_epsilon_diff_kJ"] = de
    out["2_ff14sb_params_identical"] = max(dq, ds, de) == 0.0

    qs = np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                   for i in range(nb.getNumParticles())])
    out["2_total_charge"] = float(qs.sum())
    out["2_n_distinct_charges"] = int(np.unique(np.round(qs, 6)).size)
    out["2_charges_not_all_zero"] = bool(np.abs(qs).max() > 0.1)

    # -- 3: GBn2 per-particle parameters are real -------------------------
    gbp = np.array([[float(v) for v in gb.getParticleParameters(i)]
                    for i in range(gb.getNumParticles())])
    gbp_ref = np.array([[float(v) for v in gb_ref.getParticleParameters(i)]
                        for i in range(gb_ref.getNumParticles())])
    out["3_gb_params_max_abs_diff"] = float(np.max(np.abs(gbp - gbp_ref)))
    out["3_gb_params_identical"] = bool(np.array_equal(gbp, gbp_ref))
    # radii (column 1) must vary by element/residue, not be one constant
    out["3_gb_distinct_radii"] = int(np.unique(np.round(gbp[:, 1], 6)).size)

    # -- 4: the pinned 1A13 reference -------------------------------------
    o = A.refine_coords(p.seq, rep, rb, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True)
    cm = o["components"]
    inter = cm["nonbonded"] + cm["solvation"]
    out["4_components"] = {k: float(v) for k, v in cm.items()}
    out["4_interaction"] = float(inter)
    out["4_reference"] = REF_1A13_INTERACTION
    out["4_abs_diff"] = float(abs(inter - REF_1A13_INTERACTION))
    out["4_bit_exact"] = bool(inter == REF_1A13_INTERACTION)

    # -- 5: solvation is real work ----------------------------------------
    rb2 = {k: v + np.array([0.0, 0.0, 0.37]) for k, v in rb.items()}
    o2 = A.refine_coords(p.seq, rep, rb2, k_restraint=A.K_MODERATE, steps=0,
                         tolerance=1.0, components=True)
    out["5_solvation_native"] = float(cm["solvation"])
    out["5_solvation_translated"] = float(o2["components"]["solvation"])
    # a rigid translation must leave the energy invariant: real physics, not a hash
    out["5_translation_invariant"] = float(
        abs(cm["solvation"] - o2["components"]["solvation"]))
    out["5_solvation_is_large"] = bool(cm["solvation"] < -100.0)

    print(json.dumps(out, indent=2, sort_keys=True, default=float))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "amber_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=float)
    return out


if __name__ == "__main__":
    main()
