# S26 LANE BRIEF -- PH (Physics lane)

You are lane **PH** of Sprint 26. Read `s26/LANE_CONTRACT.md` in full first; it binds you. Then
`docs/STATE_BRIEF_2026-09-12.md` in full. Then the files in section 2 below, in that order. Do not
start any code before the reading is done; a claim you did not read is a claim you will re-test.

The coordinator spawned you at 2026-09-13 00:20. Four lanes are already running (E examiner, I
infrastructure, Q quantum, P prior/learning). Phase 0 is NOT yet signed off: until the coordinator
posts `PHASE 0 SIGNED OFF` in `s26/LEDGER.md` you may read, write PREREGs and IDEAs, write and
unit-test code on synthetic data, and run 1-target memory probes under `s26/jobrun.py`. You may
also run the cis-peptide CENSUS (section 3.2, part 1) before the gate, because it computes no RMSD
and reads only omega angles of dev natives; the coordinator logs that reading in the ledger (L5).
Everything that reads an RMSD to a native waits for the gate.

## 1. What you own

1. **AMBER as a steric REJECT filter, not a ranker** (mandatory tournament direction, Part 4.2 of
   the campaign prompt). Idea, prereg, experiment, findings.
2. **The cis-peptide gap** (mandatory tournament direction). Census of the 126 dev natives, the
   projection floor on cis targets, and a DESIGN (not necessarily an implementation) of a
   two-bond-length projection if the cost is material.
3. **The steric singularity** for the report (Part IV). A short, sourced explainer plus any cheap
   measurement that sharpens it.
4. **C3, the matched-random control for AMBER refinement** (Proposal C, Part 3.3 of the prompt).
   The relaxation's gain over a matched-magnitude random displacement, fold-clustered, against MDE,
   first on the PRODUCTION built chain (artefacts already exist), later on the best C2 rung's
   output when lane P delivers it (the coordinator will message you).
5. At least **three tournament ideas** of your own (`s26/IDEA_<name>.md`), beyond the two mandatory
   ones. Physics ideas nobody has tested. Every idea gets a falsifier before it gets compute.
6. `s26/agentPH_FINDINGS.md` in the S12 to S25 format, `s26/PH_PART_IV_NOTES.md` for the report
   writer (lane E), ledger entries, hourly STATUS lines under `## PH`, commits of your own files.

## 2. Read, in this order (after the contract and the state brief)

- `s25/agentPHYS_FINDINGS.md` (617 lines) and `s25/PREREG_PHYS.md` (316). The seven-configuration
  suite; both physics energies worse than a random 75-subset (+0.330 / +0.455 A, 5/5 folds); the
  rank-permuted control; the landscape table (section 2) with every physics number the report will
  quote; the non-monotone raw-moment trap (section 3); FORM 5 closed (section 5).
- `s24/LEDGER.md` entries L10, L13, L16 (AMBER bias direction; the prior ladder; the fitted score
  closed) and the S24 D1 workstream entries on the filter/partition forms (search the ledger for
  "D1-C" and "filter"). **The functional lever already failed as a FILTER (s24 D1-C, 5/5 harmful).
  Before you run the steric reject you must state, in the prereg, exactly how a binary clash
  reject at a physical threshold differs from the S24 filter, or the Adversary will kill it.**
- `s16/LEDGER.md` L27 and `s16/repair_FINDINGS.md`, plus `s16/repair.py` (the matched-magnitude
  random displacement control and its construction; `random_rigid` at line 153). The production
  AMBER step costs +0.0207 A [+0.0143, +0.0276] on the built chain (state brief section 4); S16
  found a random displacement of matched magnitude at least as accurate as AMBER's, and AMBER's
  displacement pointing slightly AWAY from truth (cos -0.052 [-0.092, -0.013]).
- `s20/LEDGER.md` L6, L8, L12 (the shared Relax operator "H = E o Relax50", the collapse census,
  the steric singularity: participation ratio 0.0745 vs 0.4221, anisotropy 18.8 vs 5.8, 72% of
  the damage is move size) and `s20/agentC_FINDINGS.md`.
- `s13/walsh_FINDINGS.md` headline and section 3 (raw AMBER's Pauli spectrum is a delta spike; the
  nonbonded term owns all of AMBER's variance; `pauli-spectrum-delta-spike-artefact`).
- `s21/LEDGER.md` around line 2319 (a validity axis with "clash" and "cis" costs already measured
  once: read what was measured and cite it in the cis prereg).
- `s23/LEDGER.md` L5 (removing geometric outliers from the average COSTS +0.142 A) and L11 (the
  rule that forbids minimisation inside a scoring functional).
- Code: `core/amber.py` (`refine_coords` at line 1499, `memory_guard` at 252, `memory_percent`,
  `builder_for`, `AmberHamiltonian`, `convergence_flags`), `core/pipeline.py` (`relax` at 989,
  `_relax_inner` at 1014, `label` at 1109, `filter_pool` at 739, `Config`), `core/geometry.py`
  (constants lines 55 to 80: `BOND_*`, `ANGLE_*`, `OMEGA_TRANS = pi`; `build_backbone` at 146,
  `native_coords_from_pdb`), `core/project.py` (module docstring; `build_ca` at 428; the constant
  virtual bond), `core/data.py` (`Peptide` at 376 with the `rebuild` field; `DIRS` at 363 =
  `pdbs/` and `pdbs_ext/`; `_scan` at 389), `s20/qb2_lib.py` (`AmberSP` at 98, bit-exact single
  points), `s12/instrument.py`, `s24/stats_lib.py`, `s24/cache_amber/<pdb>.npz` (keys `pdb, n,
  fold, k, e_amber(500), score_dist(500), universe_idx(500), amber_verify_max_rel`).

## 3. The experiments

### 3.1 AMBER as a steric reject filter (`s26/PREREG_amber_reject.md`, `s26/ph_reject.py`)

The claim to test: both physics energies rank worse than noise, but AMBER's difficulty is a steric
singularity, so a BINARY clash reject on the production top-75 might help where a score cannot.

- Energies: read `s24/cache_amber/<pdb>.npz` `e_amber` aligned to `universe_idx`; assert on every
  target that `universe_idx` equals the instrument's production pool (`I.pool_idx`), bit-identical,
  before any number is read. No new AMBER single points are needed for the production K=500.
- Arms (all native-free): (a) reject candidates in the shipped top-75 whose AMBER single point
  exceeds a FIXED physical threshold (pre-register the threshold family, e.g. 1e3, 1e4, 1e5
  kcal/mol; S25 measured 58.6% of every pool above 1e4, so state what fraction of the top-75 each
  threshold rejects BEFORE reading any RMSD); (b) refill from the next-ranked survivors so the set
  stays at m=75; (c) do not refill, so the set shrinks. Report both.
- Controls, matched in the operator's space: reject the SAME NUMBER of candidates at random from
  the top-75 (16 draws, mean), and reject by a permuted energy vector (marginal kept, correspondence
  destroyed). Zero-information anchor: the shipped top-75 itself (3.0483 point cloud, 3.2148 built).
- Endpoint: point-cloud `rmsd_avg` AND built-chain `rmsd_arm` through `I.project` (state the basis
  on both sides). Paired, `ST.compare`, fold CI beside iid, MDE, W/L, concentration null,
  `best_of_k_within` if you sweep the threshold (a threshold sweep is a grid oracle: quote the
  split-half transfer, never the raw minimum).
- Falsifier: fold-clustered gain above MDE with 5/5 folds. State the power: at n=126 the SE of a
  paired difference on the built chain is about 0.01 to 0.04 A depending on how many targets move.
- Expected outcome, stated before the run: S24 D1-C says filters are harmful; S23 L5 says removing
  members from the average costs; S25 says AMBER's top-75 is +1.10 A more expanded than the pool.
  So the prior is "null or harmful". Say so and run it anyway; that is the point of the tournament.

### 3.2 The cis-peptide gap (`s26/PREREG_cis.md`, `s26/ph_cis.py`)

Part 1, the CENSUS (allowed before the gate): for each of the 126 dev targets
(`[t["pdb"] for t in s12.instrument.targets()]`, natives under `core.data.DIRS`), parse N, CA, C of
every residue and compute omega_i = dihedral(CA_i, C_i, N_{i+1}, CA_{i+1}); a bond is cis if
|omega| < 30 degrees. Cross-check with the CA-CA distance (cis about 2.9 A, trans 3.80 A). Also
census the retrieval pool windows: consecutive CA-CA distance in `s8/generate_univ/<pdb>.npz` `W`
(CA only, float32) below 3.3 A. Report: number of dev targets with at least one cis bond, which
residue (proline or not), and the pool's cis rate. **Never open a benchmark PDB**: iterate over the
126 dev ids only; do not call anything that enumerates `results/benchmark_manifest.json`.

Part 2 (after the gate): on cis targets, the projection floor. `core.data.Peptide.rebuild` is the
RMSD of rebuilding the native from its own phi/psi with ideal trans geometry (verify that reading in
`core/data.py` before using it; if it is something else, compute the floor yourself: project the
native through the production projection and score against the native, ORACLE DIAGNOSTIC). Compare
cis targets against non-cis targets, and compare each cis target's production `rmsd_arm` against
its `rmsd_avg` (the built-chain cost, which is where a wrong virtual bond would show). If the
floor is material (state the number and the MDE), write the DESIGN of a two-bond-length projection:
where omega enters `core/project.py` (`frames` at 305 takes a scalar omega; `build_ca` at 428),
what a per-bond omega would change, what native-free signal would choose cis (proline-preceding
positions are the only defensible one; state the base rate), and what it would cost to test.
Do not implement it on the production path.

### 3.3 The steric singularity (`s26/PH_PART_IV_NOTES.md`)

For the report's Part IV, write two pages in plain sentences for a reader who knows no protein
physics: Legacy's 11 terms (`core/energy.py`, `DEFAULT_WEIGHTS`), AMBER ff14SB/GBn2 through OpenMM
(`core/amber.py`), what "H = E o Relax50" means (s20 L6), why Legacy is a compactness model (top-75
Rg -0.758 A, rho(E,Rg) +0.60) and AMBER has the opposite sign (+1.103 A, -0.27), why both rank worse
than a random subset, what the steric singularity is (unrelaxed ideal-geometry windows evaluated
inside the r^-12 wall; 7% of modes carry all the curvature; 99.6% of the Walsh variance in ten
configurations), and the non-monotone standardisation trap (40 of 126 targets, tie blocks up to 462
of 500, retrieval order leaking through ties). Every number with its artefact path
(`s25/results/phys_landscape.json`, `s25/results/phys_suite.json`, `s20/results/...`,
`s13/results/walsh_xval.json`). One cheap new measurement is welcome if it sharpens the story
(for example the distribution of the minimum heavy-atom distance in the shipped top-75 versus in
the relaxed emission), but it must be pre-registered and it is optional.

### 3.4 C3, the matched-random control (`s26/PREREG_c3_control.md`, `s26/ph_c3.py`)

Stage 1 (after the gate; NO AMBER compute needed): the production cache
`bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` holds `ca` (built chain, `rmsd_arm` 3.2148) and
`amber_ca` (relaxed, `rmsd_full` 3.2355) for 126 targets. Per target, the AMBER displacement is
`amber_ca - ca` after superposition; its magnitude is the RMSD between the two. Build the S16
matched-magnitude random displacement control (read `s16/repair.py` and reproduce its construction,
or state exactly how yours differs): a random direction of the same per-target magnitude, applied
to `ca`, 16 draws, mean; also a matched-magnitude move along the geodesic toward another pool
member (S20's `c_land_null` form) if cheap. Report AMBER minus random, fold-clustered CI, MDE, W/L.
Also report the validity axis: clash count and bond-length deviation of `ca` versus `amber_ca`
(this is what "keep AMBER as a validity step" would rest on; measure it so the presentation can
say it with a number).

Stage 2 (when the coordinator messages you with lane P's best C2 rung output): re-run the
production relaxation (`core.pipeline.relax` semantics: `refine_coords(..., k_restraint=10,
steps=0)`, converged) on that output, 126 targets, under `jobrun --tag AMBER`, checkpointing per
target to `s26/results/ph_c3_cells/<pdb>.json` so a governor kill loses nothing, then the same
control. Probe one target first and quote the peak RSS. Never more than two AMBER jobs at once
(the governor enforces it; you also must not launch a third).

Decision rule (from the prompt): if the relaxation's gain over the random control is below MDE,
the presentation drops "refine with physics" as an accuracy step and keeps it as a validity step.
Write that sentence, with the numbers, in `s26/C3_RESULT.md` so lane P can cite it.

## 4. Rules that bite hardest in this lane

- Basis discipline. `rmsd_avg` is a point cloud, `rmsd_arm` a built chain, `rmsd_full` a relaxed
  chain. Never subtract across bases without saying so. The production result is the built chain.
- Controls must match the operator's space (the project's most repeated error). A reject filter's
  control rejects the same count; a displacement's control moves the same distance.
- ORACLE labelling on anything that reads a native (the census, the floor, the cos to truth).
- Ties: `ST.argmin_tied`, never `np.argmin` on a tied signal (`consensus-is-the-only-in-band-
  discriminator` memory: an argmin on ties read the retrieval order and invented a winner).
- Every AMBER job: `python s26/jobrun.py --agent PH --tag AMBER --name <name> --est-ram <GB> --
  python s26/ph_<x>.py ...`. Probe 1 target first. Checkpoint every 10 minutes.
- Do not edit `core/`, `s16/`, `s20/`, `s24/`, `s25/`. New code under `s26/` only.
- No stock words, no em dashes in anything the presenter or the report will read.
- If you finish everything, take the highest-ranked orphaned tournament idea (ask the coordinator
  by finishing your turn with the question).

## 5. Deliverables checklist

- `s26/IDEA_amber_reject.md`, `s26/IDEA_cis_peptide.md`, plus at least three more `s26/IDEA_*.md`.
- `s26/PREREG_amber_reject.md`, `s26/PREREG_cis.md`, `s26/PREREG_c3_control.md` (before each run).
- `s26/ph_reject.py`, `s26/ph_cis.py`, `s26/ph_c3.py`; results in `s26/results/ph_*.json` with
  `ST.save_atomic` provenance.
- `s26/agentPH_FINDINGS.md`, `s26/PH_PART_IV_NOTES.md`, `s26/C3_RESULT.md`.
- Ledger entries `## L<n> -- TITLE (date, PH)` for every finding; STATUS lines hourly.
- Commits of your own files on branch `s26` with the trailer lines from the contract.

Finish your turn with: what you found, what you did not do and why, and any question for the
coordinator. Do not predict results you have not measured.
