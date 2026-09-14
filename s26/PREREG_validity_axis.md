# PREREG -- THE HEAVY-ATOM VALIDITY AXIS OF THE PRODUCTION RELAXATION (lane PH, S26; coordinator's item under L77)

Written 2026-09-13 23:10, before `s26/ph_validity.py` produced a number. Never edited after.
Code: `s26/ph_validity.py`. Results: `s26/results/ph_validity.json`, cells under
`s26/results/ph_validity_cells/`. Tag AMBER (one relaxation per target).

## 1. What is measured and why

"Refine with physics is kept as a validity step" (C3, L39, L46) rests so far on energies
(e0 above 1e4 on 58.7% of targets, e1 below 1000 on 125 of 126) and on CA-level geometry
(virtual bond, L24). The heavy-atom axis the record uses for validity (`s16.energy_lib.panel`:
`n_clash_2A`, `n_clash_2p6A`, `min_heavy`, `bond_strain`, `angle_strain`, `rama_favoured`,
`rama_outlier`, `cis_frac`, `chirality_L_frac`) has never been printed for the production
emission before and after its relaxation. The production cache holds the relaxed CA, phi and
psi but not the relaxed N, C, O, so the relaxed heavy-atom backbone must be regenerated.

## 2. Operator, gate, and what is compared

Per target: the built chain from the cached production `phi`, `psi` through
`core.geometry.build_backbone` (ideal geometry, N/CA/C/O/CB), panel BEFORE; the production
relaxation re-run exactly as `core.pipeline._relax_inner` does it (`torsion_lib2.library_for(seq,
8, seq)`, `PerResidueTorsion`, `refine_coords(k_restraint=10, steps=0, components=True)`),
panel AFTER on the returned N/CA/C/O backbone (CB from `place_cb`). GATE, per target: the
re-run's CA must equal the cached `amber_ca` to max |dCA| < 1e-6 A and its `e1` the cached
`amber_e1` to 1e-6 kcal/mol; a failure is reported by name and the target is kept with a flag
(the relaxation is deterministic at threads = 1 and S16 verified bit-identity across
interpreters; a mismatch would mean the cache was produced under different settings, which is
a finding about the cache, not about the operator).

Contrasts, paired per target, `ST.compare` with folds, before minus after, on every panel axis;
the zero-information reference for validity is the constant alpha-helix (S16 L27: rama 1.000
and zero clashes by construction), printed beside the built chain and the relaxed chain so no
validity statistic is credited to the force field that a constant helix also achieves. The
"validity step" sentence is written from: the number of targets with any heavy-atom pair below
2.0 A before and after; the mean clash count before and after; the bond and angle strain before
(ideal by construction, expected 0.0 on the built chain) and after; and the two broken-bond
targets (2BP4, 9KAR) by name.

## 3. Falsifier

None in the accuracy sense; this is a descriptive measurement. The claim "the relaxation
removes the builder's clashes" is FALSE if the mean `n_clash_2A` after relaxation is not below
half of the mean before with the fold CI excluding zero, or if `bond_strain` after exceeds 0.05
(5% relative) on more than 10 targets, in which case "validity step" is reworded to say what
the numbers say.

## 4. Cost

126 x (2.8 s builder + 6 to 12 s minimisation) = about 25 min, one process, est-ram 0.4 GB
(the branch probe's peak 0.272 rounded up), `--tag AMBER`, per-target cells, resumable. Runs
after `ph_branch_relax2` releases the AMBER slot (one AMBER job of mine at a time). Probe: the
first target's wall and the gate result are printed before the loop continues.
