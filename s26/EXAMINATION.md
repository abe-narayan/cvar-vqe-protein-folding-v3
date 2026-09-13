# S26 EXAMINATION (lane E, Examiner / Librarian)

Written 2026-09-13 (reading 00:03 to 00:46, assembly 08:40 onward) on branch `s26`, HEAD `5dc7a3a6`
at assembly time. Every number below carries the artefact it was read from. The scripts that
produced the artefacts are committed (`e3570aed`): `s26/e_module_map.py`, `s26/e_hashes.py`,
`s26/e_reproduce.py`, `s26/e_trace.py`, `s26/e_claims.py`. `python s26/examine.py` (lane I)
regenerates the module map through `e_module_map.py` and, with `--search`, the claim search.

Reading order followed as mandated: LANE_CONTRACT, STATE_BRIEF, README, ARCHITECTURE,
professor_brief, s25/LEDGER, s25/QUANTUM, s25/agentQ_FINDINGS, docs/FINDINGS (1-300, every
heading, S9-10, S10-4, S10-5, S11 in full), CONDENSED_REPORT, every `core/*.py`, the 24 root
modules, `s5/ s7/ s8/ s9/`, `tests/*.py`, `verify/*.py`, `s25/resultslab/*.py`,
`s12/instrument.py`, `s24/stats_lib.py`. Nothing was skimmed; every file was opened in full.

Two incidents, recorded here and in `s26/agentE_FINDINGS.md`: (1) the first run of `e_claims.py`
parsed `s9/final_report.json` before the benchmark exclusion existed and printed one aggregate
leaf (`dist/shipped/mean`), a number already published in `README.md:14`; the exclusion was
added, the artefacts were regenerated, and no per-target benchmark value was printed or read.
(2) The session was cut at 00:46 by the API limit; nothing was lost.

---

## A. MODULE MAP

Source: `s26/results/module_map.json` (699 modules under the repository root; `.git`,
`__pycache__`, `_archive`, `QUARANTINE*` and `distogram_models_large.STALE*` skipped). For each
module the JSON holds path, line count, sha256, first docstring line, every repository import
(module level and lazy, with line numbers), the reverse graph `imported_by`, and the backend
names the module answers to. The two backend reports were produced in two subprocesses by
`core.backend_report()`:

```
CORE_BACKENDS unset  {"amber": "core.amber", "data": "core.data", "energy": "core.energy",
                      "geometry": "core.geometry", "legacy": "core.energy",
                      "numerics": "core.geometry", "predict": "core.predict",
                      "project": "core.project", "quantum": "core.quantum"}
CORE_BACKENDS=legacy {"amber": "amber_refine", "data": "peptide_db", "energy": "energy_terms",
                      "geometry": "protein_geometry", "legacy": "legacy_field",
                      "numerics": "s7.audit", "predict": "distogram", "project": "s8.project",
                      "quantum": "foldvqe"}
```

"importers" is the count of repository modules that import the module (full lists in the JSON).

### A.1 `core/` (11 modules, the production path)

| module | lines | computes | repo imports | importers | backend |
|---|---:|---|---|---:|---|
| `core/__init__.py` | 131 | the backend switch: `CONTRACT`, `REPLACES`, `backend(name)`, `backend_name`, `backend_report`; `CORE_BACKENDS=legacy` forces the root modules | none | 14 | (switch) |
| `core/amber.py` | 1697 | ff14SB + GBn2 through OpenMM: `builder_for`, `refine_coords` (k_restraint, steps, tolerance, components), `refine_many`, `memory_guard` (92% ceiling), `convergence_flags` (`CONVERGE_MAX_KCAL` = 1000), bounded builder cache, all-20 sidechain builder | protein_geometry, budget, floor | 36 | amber |
| `core/bench.py` | 1090 | the speed/science harness behind `bench_results/*.json`: arms, per-stage timings, occupancy and contention verdicts, `MACHINE` core map (6.426 core-equivalents), `ABLATIONS`, `("rmsd_avg", "raw average (illegal)")` at line 677 | core, core.pipeline, s9.final, s7.audit, s7.debias, s7.poolsize, s8.inband, protein_geometry, torsion_lib2 (12 lazy) | 1 | none |
| `core/cache.py` | 355 | deterministic npz cache with parameter sidecars; `key`, `store`, `load`, `cached`, `CacheCollision`, in-process `mem_pct` | none | 7 | none |
| `core/data.py` | 1409 | peptides and windows, `identity` (longer normalisation, line 182), `containment` (216), `identity_many` (321), pinned `clusters`/`folds`, `benchmark`/`dev_set`, encodings and representations, ESM accessors (`esm_raw`, `esm_embed`, `esm_contacts`), `fold_fragments` | core.cache, core.geometry, sidechains | 12 | data |
| `core/energy.py` | 1762 | the Legacy 11-term energy: `energy_components`, `components_batch`, `totals_batch`, `BatchLegacy`, `LegacyField` (off in production) | peptide_db, protein_geometry, representations | 19 | energy, legacy |
| `core/geometry.py` | 1099 | `build_backbone(_batch)` (constants lines 57-65, `OMEGA_TRANS` = pi at 65), `kabsch_rmsd_batch`/`ca_rmsd_batch` (singular-value residual, reflections forbidden), PDB IO with cache and access log, DSSP, `CATracePrior`, `representation_floor` | core.cache, core.data | 41 | geometry, numerics |
| `core/pipeline.py` | 1683 | `Config` (118) and `PROD` (241); manifests; stages `retrieve` (672), `score` (711), `filter_pool` (739), `_zrank` (788), `quantum_stage` (806), `average` (924), `project` (940), `run_target` (1049), `label` (1109); `_q` (666) float32 round trip; `guard_esm` (535); `VQE_LFO` (113); `run` with workers, resume and canonical-order summaries | core, s7.debias, esm_features, core.data, s7.poolsize, s8.consensus2, torsion_lib2, s8.inband (7 lazy) | 20 | none |
| `core/predict.py` | 853 | the distogram: `BIN_EDGES` (54, 17 bins), `_fold_fragments` (224), `_model_path` (339), `train_fold` (345), `Distogram` (396) with `_risk` (npairs, 760) on the 2.0..40.0/0.05 grid, `w = shell/(sd+0.5)^g`, `for_target` (426) | core.cache, core.data | 19 | predict |
| `core/project.py` | 1598 | STAGE 3b: geometry constants (149-155), scan-builder `frames`/`build_ca`/`build_ca_exact`, gradients `exact`/`analytic`/`fd` (`GRAD = "exact"`, 170), penalties `rama`/`rama20`/`vm`/`phip`/`ramah` (`make_penalty` 759), `STARTS` (814), `_emit` (817), `fit_prior`/`fit_multi`/`lam_path` (933), `geometry_report` (1045), `l_signature`; harvest/equiv/exactness stages | core, s8.project, core.pipeline (2 lazy) | 40 | project |
| `core/quantum.py` | 2613 | `StatevectorCircuit`, `OneLayerAnsatz`, `MPSAnsatz` (exact, chi = 2^layers), `cvar`/`cvar_from_samples`/`cvar_from_probs`/`cvar_exact`, `tail_indices`, `grad_cvar_paramshift`/`_fd`/`_score` (default baseline `const`), `free_energy`, `run_cvar_vqe`, `run_vqe` (SPSA / cvar_grad), `Reservoir`, `basin_select`, `FoldObjective`, `FoldingHamiltonian`, PennyLane `build_global_circuit`, `all_bitstrings` (capped at 20 qubits) | budget, core, sidechains, floor, distogram, pairnet, catrace, peptide_db, priors (9 lazy) | 53 | quantum |

### A.2 The 24 root reference modules (the legacy arm)

| module | lines | computes | repo imports | importers | backend |
|---|---:|---|---|---:|---|
| `amber_hamiltonian.py` | 443 | OpenMM System builder for discrete-torsion candidates (legacy arm of `core.amber`) | protein_geometry, sidechains, budget | 2 | none |
| `amber_refine.py` | 269 | restrained all-atom refinement `refine_coords`, `builder_for`, `K_MODERATE` = 10 | floor, amber_hamiltonian | 5 | amber (legacy) |
| `budget.py` | 193 | `BudgetExhausted`, `BudgetedEnergyModel`, `resolve_maxiter`, `check_optimizer_budget`; single definition asserted by `tests/test_amber.py` | none | 9 | none |
| `catrace.py` | 158 | CA-trace prior: `pseudo_angles`, handedness, mirror penalty (`catrace_prior.npz`) | peptide_db | 3 | none |
| `distogram.py` | 491 | legacy distogram: MLP, `_fold_fragments`, `Distogram.for_target`, Bayes-risk `score` | peptide_db, priors, fragment_db, esm_features | 22 | predict (legacy) |
| `energy_terms.py` | 783 | Legacy 11 terms, `TERM_NAMES`, `DEFAULT_WEIGHTS`, MJ table, `sequence_arrays` | protein_geometry, representations | 7 | energy (legacy) |
| `esm_features.py` | 147 | ESM-2 650M per-residue features and contact head; `compute` (tokens only), caches | none | 9 | none |
| `floor.py` | 73 | representation floor `floor`, `project` states | protein_geometry | 5 | none |
| `foldvqe.py` | 681 | prior, warm-started CVaR-VQE `run_vqe`, `Reservoir`, `basin_select`, `marginals_to_angles` | distogram, objective, peptide_db, priors, protein_geometry, qansatz, torsion_lib2, floor, pairnet | 1 | quantum (legacy) |
| `fragment_db.py` | 143 | protein fragment corpus from `prots/` (`fragment_db.npz`) | peptide_db, protein_geometry | 13 | none |
| `hamiltonian.py` | 133 | `FoldingHamiltonian` over discrete encodings with the budget | energy_terms, protein_geometry, sidechains, budget | 1 | none |
| `legacy_field.py` | 315 | `BatchLegacy` generation field, `FITTED_WEIGHTS`, `SUBSETS`, `TERMS` | energy_terms, protein_geometry | 7 | legacy (legacy) |
| `legacy_refine.py` | 537 | batched continuous Legacy energy, `refine_angles` (bounded pattern search) | energy_terms, protein_geometry | 1 | none |
| `objective.py` | 150 | `FoldObjective` composite (distogram, clash, MRF) in batch | protein_geometry | 2 | none |
| `pairnet.py` | 349 | pair-tensor distance predictor (`pairnet_models/`) | distogram, peptide_db, priors, fragment_db, esm_features | 2 | none |
| `peptide_db.py` | 429 | legacy database: `Peptide`, `identity` (119, `shared / max(len(a), len(b))` at 170), `holdout`, `clusters`, `folds`, `benchmark`, `dev_set`, `by_pdb` | protein_geometry | 86 | data (legacy) |
| `priors.py` | 238 | pair features, `BIN_EDGES`, `_PROPS`, `DistogramPrior`, `TorsionMRF` | peptide_db, representations | 8 | none |
| `protein_geometry.py` | 618 | `build_backbone(_batch)`, Kabsch, PDB IO, `extract_torsions`, DSSP, SS assignment | none | 43 | geometry (legacy) |
| `qansatz.py` | 324 | `OneLayerAnsatz`, `MPSAnsatz`, `cvar`, `cvar_gradient` (the tail-only baseline, the recorded defect) | none | 3 | none |
| `refine2.py` | 752 | analytic-gradient continuous refinement of (phi, psi) | protein_geometry, distogram | 2 | none |
| `representations.py` | 626 | bitstring to structure encodings (torsion libraries, lattice), `make_representation` | protein_geometry, sidechains | 10 | none |
| `sidechains.py` | 521 | all-20 sidechain builder `build_full_structure`, residue bond tables | protein_geometry | 6 | none |
| `torsion_lib2.py` | 199 | per-residue torsion state libraries `library_for`, `PerResidueTorsion` | peptide_db, representations | 45 | none |
| `vqe.py` | 568 | SPSA driver `_spsa`, `cvar_from_samples`/`_distribution`, `n_parameters`, `MIN_USEFUL_SPSA_ITERS` | budget | 1 | none |

### A.3 `s5/ s7/ s8/ s9/`

| module | lines | computes | repo imports | importers | backend |
|---|---:|---|---|---:|---|
| `s5/lib.py` | 114 | `B62`, `encode`, `windows_full`, `kabsch_rmsd_batch` (retyped by hand; pinned equal to `core.data.BLOSUM62`) | distogram, peptide_db | 8 | none |
| `s5/esm32.py` | 62 | compact 32-d ESM cache (`s5/esm32.npz`) built in a throwaway subprocess | distogram, esm_features, fragment_db, peptide_db (lazy) | 2 | none |
| `s5/esmraw.py` | 61 | float16 raw 1280-d ESM cache (`s5/esmraw.npz`) | esm_features, fragment_db, peptide_db (lazy) | 1 | none |
| `s5/torsion.py` | 314 | sequence-conditioned von Mises mixture, `for_target`, `mixture_logpdf` (the `vm` penalty) | distogram, peptide_db, s5.esm32, esm_features | 1 | none |
| `s7/audit.py` | 828 | the tuning instrument's primitives: `kabsch_rmsd_batch`, `pair_index`, `pair_dists`, `encode`, `B62`, `build_target`, `build_pool_members`, `windows_of`, `KMAX` | distogram, fragment_db, peptide_db, protein_geometry | 15 | numerics (legacy) |
| `s7/debias.py` | 607 | `score_risk`, `rank_stats`, `tuning_targets` (126), `dev_targets` (24), `paired`, `K` = 500, `s7/debias_cache` | distogram, peptide_db, s7.audit | 14 | none |
| `s7/poolsize.py` | 693 | `GRID` (2.0..40.0 step 0.05), `load_pred` (the cached shipped prediction per target) | distogram, peptide_db, s7.audit, s7.debias, esm_features | 9 | none |
| `s7/amber_native.py` | 237 | `pool_for`, `amber_energy` (interaction-only, converged under k = 10) | distogram, legacy_field, peptide_db, protein_geometry, torsion_lib2, s5.lib, amber_refine | 3 | none |
| `s7/repr_select.py` | 922 | representation ablation on selection (onehot/phys/pca32/pca128/raw/escon); `train_entries` (the production training rule) | distogram, peptide_db, priors, s7.audit, s7.debias, s5.esm32, s5.esmraw, esm_features, fragment_db | 1 | none |
| `s8/consensus2.py` | 1293 | filter x size x operator sweep; `filter_orders`, `op_value`, `_avg_struct`, `fit_torsions`, `fit_consensus` (`FIT_STARTS`), `DEV_ARM` (sc, 75, fit) | protein_geometry, s7.audit, s7.debias, s7.poolsize, s8.generate, s8.inband, legacy_refine | 4 | none |
| `s8/generate.py` | 1118 | the window universes (`s8/generate_univ`, `stage_univ` asserts against `s7/debias_cache`), generator arms, `shipped_score`, `pair_D`, `rg_of`, `rama_logp`, the law `sel = 0.257 best + 0.311 mean + 1.673` | distogram, peptide_db, protein_geometry, s7.audit, s7.debias, s7.poolsize | 5 | none |
| `s8/inband.py` | 1387 | 38 native-free signals (`signal_block`), `sel_of` (tie-averaged), in-band skill, `stage_two` (score filter + consensus medoid), `stage_devpool` | energy_terms, legacy_field, peptide_db, protein_geometry, s7.audit, s7.debias, s7.poolsize, s8.generate, distogram | 6 | none |
| `s8/integrate.py` | 1928 | channel error-correlation matrix, `Fuser`, filter x q grid, `Circuit`, `cvar_exact`, `grad_cvar_paramshift`/`_fd`/`_score`, `free_energy`, `run_cvar_vqe`, `stage_vqe` (`s8/integrate_vqe.json`), roles, refine, k25, natcons | protein_geometry, s8.generate, s8.inband, s8.invfold, amber_refine, torsion_lib2, peptide_db, s7.poolsize, s5.lib | 0 | none |
| `s8/invfold.py` | 1391 | inverse-folding self-consistency scorer (`features`, `InvFoldNet`, `SelfConsistency`), `build_pool` with torsions, leakage audit | peptide_db, protein_geometry, distogram, s5.lib, s7.debias | 1 | none |
| `s8/project.py` | 1732 | torsion-constrained projection: `RamaPenalty`/`RamaHingePenalty`/`VMPenalty`/`PhiPosPenalty`, `logp_tables`, `hinge_thresholds`, `fit_prior`, `fit_multi`, `lam_path`, the reference `core.project` is bit-identical to | protein_geometry, s7.audit, s7.debias, s7.poolsize, s8.consensus2, s8.generate, s8.inband, distogram, peptide_db | 3 | project (legacy) |
| `s9/final.py` | 1231 | the 60-target pass (spent): `K` 500, `M` 75, `PEN` ramah, `LAM` 0.3, `AMBER_K` 10, `AMBER_STEPS` 0; `library_members`, `windows_all`, `build_pool`, `synthesise`; `s9/final_report.json` (never opened here) | distogram, peptide_db, protein_geometry, s7.audit, s7.debias, s7.poolsize, s8.consensus2, s8.project, amber_refine | 3 | none |

### A.4 Results lab, instrument, statistics

| module | lines | computes | repo imports | importers | backend |
|---|---:|---|---|---:|---|
| `s25/resultslab/__init__.py` | 6 | package docstring | none | 0 | none |
| `s25/resultslab/exportlib.py` | 985 | PDB export with REMARK 999 provenance, `quantise`, `read_pdb`, `native_ca` (from the universe), `pool_best` (ORACLE, gate only), `pool_oracle_gate` (598), `difficulty_gate` (653), `ca_bond_stats`, `chain_from_ca` (lam 0.3 through `core.project.lam_path`), `export_pair` | s12.instrument, s24.stats_lib, core.geometry, core.project | 6 | none |
| `s25/resultslab/providers.py` | 229 | provider registry; `incumbent` replays `bench_results/cache/1fc9f2dcf489e2fb`; `synthetic` (fixture-only); `from_json`/`from_npz` (provenance UNKNOWN by default) | s12.instrument, s25.resultslab.exportlib | 3 | none |
| `s25/resultslab/schema.py` | 526 (510 at reading; lane I commit `eb89c165`) | records and leaderboard (`RECORD_KEYS`, `LEADERBOARD_KEYS`), `collect`, `leaderboard` (gate at 309), `build` (three refusals), `fmt_leaderboard` | s12.instrument, s24.stats_lib, exportlib, providers | 1 | none |
| `s25/resultslab/site.py` | 140 | static site renderer (`results/site`), structures inlined | s24.stats_lib, exportlib | 1 | none |
| `s25/resultslab/build.py` | 178 | `--mode selftest` (sandbox, fixtures) and `--mode frozen --spec` | exportlib, providers, schema, site | 0 | none |
| `s25/resultslab/test_export.py` | 499 | the export layer's tests (round trip, gates, provenance, basis) | s12.instrument, exportlib, providers, core.geometry | 0 | none |
| `s12/instrument.py` | 289 | `ca_rmsd`, `kabsch_rmsd_batch`, `targets()` (sorted glob of `s8/generate_univ`), `load_univ`, `pool_idx`, `shipped_record` (`PROD_KEY` = 1fc9f2dcf489e2fb), `shipped_score`, `project`, `selfcheck`, `paired`, `write`; `K` 500, `M` 75, `BAND` 1.5, `FAIL18` | core.project, core.geometry, core.pipeline, core.predict (lazy) | 469 | none |
| `s24/stats_lib.py` | 519 | `compare` (SE, MDE = 2.8016 SE, iid and fold CIs, W/L, concentration null, verdict), `fmt`, `pinned_folds`, `best_of_k_within`, `split_half_transfer`, `argmin_tied`, `save_atomic`, `provenance` | s15.seed, s12.instrument | 27 | none |

Observations from the graph: `s8.integrate` has zero importers (its numbers live in
`s8/integrate_vqe.json`, which `core/pipeline.VQE_LFO` and `verify/vqe_lfo_audit.py` read);
`s12.instrument` is imported by 469 modules and is the ruler every sprint since S12 reports with;
`peptide_db` (86 importers) and `torsion_lib2` (45) are the most-imported root modules, both
still imported directly by consolidated modules (`core.energy`, `core.quantum`, `core.pipeline`),
so the "legacy arm" is not a separable layer.

---

## B. DATAFLOW TRACE

Script `s26/e_trace.py`, run under the governor (`s26/jobs_done/e_trace_1S9Z.json`: wall 15.1 s,
exit 0, peak RSS 0.552 GB; `e_trace_9KAR.json`: 15.0 s, exit 0, 0.553 GB; in-process peak RSS
535.9 and 535.2 MB). Artefacts `s26/results/e_trace_1S9Z.json`, `s26/results/e_trace_9KAR.json`.
The trace mirrors `core.pipeline.run_target` stage by stage with `cfg = pl.PROD` and, for the
selector only, `dataclasses.replace(pl.PROD, quantum=True)` (production runs `quantum=False`, so
the selector does not execute in the 3.2148 arm; it is traced because the brief asks for it).
Backends in effect: the default report above. Both targets are in the hot ESM cache
(`esm_small.npz`); the 1.5 GB bank was not loaded.

T030 is `1S9Z` (`results/summary/target_map.json`; the 0.182 A target). The hard target is
`9KAR` (rmsd_arm 7.4377, the worst pool_best-to-emission gap in the cache: pool best 1.3997,
top-75 best 5.3436), chosen because it shows the failure mode that dominates the FAIL18 stratum.

| stage | quantity | 1S9Z (T030) | 9KAR |
|---|---|---|---|
| sequence | `seq`, n, fold (pinned) | SIRELEARIRELELRI, 16, fold 1 | GGWGTVPDWFFNMNW, 15, fold 3 |
| ESM-2 | cache; raw shape; contacts; pca32; pair features | `esm_small.npz`; (16, 1280); (16, 16); (16, 32); (105, 183) | `esm_small.npz`; (15, 1280); (15, 15); (15, 32); (91, 183) |
| distogram | model file (sha256 first 16); bins; centres | `distogram_models/fold1_esm_frag.pt` (23691cc531435766); 17; 4.0, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.5, 9.5, 10.5, 11.75, 13.25, 15.0, 17.5, 21.0, 25.0 | `distogram_models/fold3_esm_frag.pt` (25ba58715917945b); 17; same centres |
| distogram | `prob` shape; expected min/max; sd min/max; `_risk` shape | (105, 17); 5.0475 / 21.0692; 0.2318 / 2.6985; (105, 760) | (91, 17); 5.6155 / 23.4415; 0.2325 / 7.5805; (91, 760) |
| retrieve | n_windows; K; W shape; sim min/max; from peptides; distinct parents | 7193; 500; (500, 16, 3); -1.0 / 19.0; 144; 420 | 9814; 500; (500, 15, 3); -6.0 / 32.0; 160; 414 |
| score | sc shape; min / max; n_pairs; argmin; distinct scores | (500,); 0.692609 / 7.268827; 105; 401; 488 | (500,); 2.242029 / 5.179138; 91; 324; 483 |
| filter | m; sub shape; `top` shape (prod / quantum); Pt shape; sub first 10; matches record | 75; (75,); (75,) / (128,); (75, 75) / (128, 128); 401, 98, 193, 161, 284, 260, 372, 373, 381, 382; True | 75; (75,); (75,) / (128,); (75, 75) / (128, 128); 324, 201, 45, 73, 89, 137, 248, 485, 416, 124; True |
| Hamiltonian | H = diag(E), E = `_zrank(sc[top[:128]])`; E first 5 | -1.718577, -1.691512, -1.664448, -1.637384, -1.610320 | -1.718581, -1.691517, -1.664453, -1.637389, -1.610324 |
| Hamiltonian | E last 5 | 1.610320, 1.637384, 1.664448, 1.691512, 1.718577 | 1.610324, 1.637389, 1.664453, 1.691517, 1.718581 |
| Hamiltonian | E range; max deviation from the standardised ladder 1..128; as fraction of range; ties in prefix | 3.437153; 0.013536; 0.394%; 2 | 3.437163; 0.013540; 0.394%; 4 |
| circuit | n; layers; params; iters; seed; (alpha, T) from `VQE_LFO[fold]` | 7; 3; 21; 50; 0; (0.25, 0.3) | 7; 3; 21; 50; 0; (1.0, 0.3) |
| CVaR | CVaR value; entropy (bits); collapsed; p_max; realised tail size; prefix length at alpha | -1.673974; 6.2858; False; 0.0691; 8; 8 | -1.006803; 5.6576; False; 0.0344; 128; 128 |
| readout | selected pool index (p-weighted / uniform) | 449 / 449 | 449 / 449 |
| average | avg_ca shape; medoid (local / pool index); virtual bond mean / min / max; max diff vs record | (16, 3); 44 / 449; 3.7529 / 3.3581 / 3.7974; 0.0 | (15, 3); 34 / 449; 1.9696 / 1.1556 / 2.6674; 0.0 |
| projection | ca shape; bond mean / min / max; max diff vs record (ca, phi, psi, fit_ca) | (16, 3); 3.803955 / 3.803955 / 3.803955; 0.0 / 0.0 / 0.0 / 0.0 | (15, 3); 3.803955 / 3.803955 / 3.803955; 0.0 / 0.0 / 0.0 / 0.0 |
| AMBER (from the record; not re-run, RAM below the 85% threshold rule was not needed) | e0; e1; moved; strain after; amber_ca bond mean | 5370.04; -1290.59; 0.1913; 51.99; 3.8619 | 6.020e12; 1262.41; 0.6097; 1166.13; 3.9338 |
| reporting (`s12.instrument.ca_rmsd` vs `pl._q(native)`) | rmsd_avg; rmsd_fit; rmsd_arm; rmsd_full (fresh = record, diff) | 0.195888; 0.181943; 0.181981; 0.317379 (all diff 0.0) | 7.060979; 7.441044; 7.437696; 7.543710 (all diff 0.0) |
| record | shipped; pool_best; top_m_best | 0.302322; 0.159596; 0.159596 | 7.230870; 1.399733; 5.343629 |
| fresh `run_target` + `label` (amber=False) | shipped / pool_best / top_m_best / rmsd_avg / rmsd_fit / rmsd_arm diff; `sub` identical; max diff ca | all 0.0; True; 0.0 | all 0.0; True; 0.0 |
| timings (s) | retrieval / windows / distogram / pairwise / projection / quantum | 0.001 / 0.367 / 0.012 / 0.090 / 4.701 / 0.047 | 0.002 / 0.359 / 0.010 / 0.077 / 4.180 / 0.041 |

What the trace shows. The Hamiltonian is the standardised rank ladder of the top-128 scores to
within tie-averaging on both targets (0.394% of range; S25 L17 quotes a worst case of 1.18% over
8 targets, `s25/results/q_gibbs.json::results/spectrum_target_independence/worst_frac_of_range`
= 0.011821117271971639). On fold 3 the deployed table sets alpha = 1.0, so the CVaR tail is
the whole register (realised tail 128 of 128) and the selector reduces to a free-energy
minimisation with no tail constraint; `s25/results/q_verify.json::results/"folds running alpha
== 1.0 (NO tail constraint)"` = [0, 3, 4], and `s25/results/q_alpha.json::results/
share_of_targets_with_no_tail_constraint` = 0.6190476190476191. On 9KAR the retrieval holds a
1.40 A window but the top-75 score filter keeps nothing under 5.34 A; the coordinate average
of that set has a mean virtual bond of 1.97 A (minimum 1.16 A), the projection rebuilds it at
3.804 A, and the AMBER step starts from 6.0e12 kcal/mol and ends at 1262 kcal/mol, above the
1000 kcal/mol convergence gate of `core.amber.CONVERGE_MAX_KCAL`.

---

## C. CLAIM LEDGER

Method: `s26/e_claims.py` (three passes: a literal grep of every md/py/csv/txt/json under the
repository, a walk of every JSON leaf under the results directories for values that round to
the claim at the quoted precision with the key path recorded, and targeted lookups in the citing
sprint's own results directory). Output `s26/results/claim_search.json` and the readable
`s26/results/claim_search.txt`. `s9/final_report.json`, `s9/final_synth.json`, `s9/final_cache`
and `results/benchmark_manifest.json` are excluded (Rule 1). Status vocabulary: SOURCED (a
results artefact holds the number; file and key given), SOURCED BY TEST (an assertion in
`tests/` pins the number and passed in lane I's governed run), DERIVED (computed from an
artefact's rows; the arithmetic is stated), **UNSOURCED (document-only)** (no artefact holds it).
Stored values are quoted to stored precision.

| id | claim | artefact (file :: key) | stored value | status |
|---|---|---|---|---|
| C01 | 3.2148 built chain (`rmsd_arm`) | `bench_results/baseline_tuning126.json :: science/rmsd_arm/mean`; `optimised_tuning126_w8.json :: science/rmsd_arm/mean`; `s26/results/e_reproduce.json :: summary/rmsd_arm/mean_fresh` | 3.214765154210998; 3.2147651542109985; 3.2147651542109985 | SOURCED |
| C02 | 3.0483 point cloud (`rmsd_avg`) | `bench_results/baseline_tuning126.json :: science/rmsd_avg/mean`; `e_reproduce.json :: summary/rmsd_avg/mean_fresh` | 3.048338093879532; 3.048338093879532 | SOURCED |
| C03 | 3.236 / 3.2355 repaired emission (`rmsd_full`) | `bench_results/baseline_tuning126.json :: science/rmsd_full/mean`; `e_reproduce.json :: summary/rmsd_full/mean_fresh` | 3.2354598538973844 | SOURCED; note that 1 of the 126 relaxations (9KAR, e1 1262.4103583126937) ends above the 1000 kcal/mol gate of `core.amber.CONVERGE_MAX_KCAL` and a second (2BP4, e1 845.3467373422843) converges by that gate but carries bond+angle strain 1172.6967039637698, above the S8 strain-rejection rule (`s8/integrate.py:292`); both are inside this mean (`bench_results/cache/1fc9f2dcf489e2fb/{2BP4,9KAR}.json`; lane PH's L24: 125 of 126 converge) |
| C04 | 2.9610 benchmark full system | `s9/final_report.json :: dist/full/mean` (not opened, Rule 1); asserted within 5e-4 by `tests/test_pipeline.py:338-347` and `tests/test_equivalence.py:330-342`, both passed in `s26/TEST_RUN.md` job `pytest_core` | not read | SOURCED BY TEST |
| C05 | 2.9507 benchmark shipped baseline | same file, `dist/shipped/mean`; same tests | not read (the first `e_claims.py` run printed 2.9507235775391263 before the exclusion; incident recorded) | SOURCED BY TEST |
| C06 | +0.0103 paired, CI [-0.1596, +0.1803], 31W/29L | `s9/final_report.json` (not opened); cited `docs/FINDINGS.md:4505` with paired SE 0.0867, `docs/FINDINGS.md:36`, `README.md:14` | not read | artefact named, not verified by E (Rule 1) |
| C07 | 1.711 / 1.7108 pool best | `bench_results/baseline_tuning126.json :: science/pool_best/mean`; `results/summary/pool_best.json :: mean`; `e_reproduce.json :: summary/stored_pool_best/mean` | 1.7108244199364904 (all three) | SOURCED |
| C08 | 3.425 / 3.4251 random 75-subset null | `s25/results/phys_suite.json :: random_null_rank/mean` (and `random_null_moment/mean`) | 3.4251256499338507 | SOURCED (point-cloud basis) |
| C09 | 4.065 constant alpha-helix | `s12/results/s14_ladder.json :: rows/L0_constant_helix/mean`; `s12/results/s14_coherence.json :: real_emitters/L0_constant_helix/rmsd` | 4.064753929494389 | SOURCED |
| C10 | 3.770 sequence-only torsion predictor | `s13/results/tors_eval.json :: reports/a_pepPos/build/mean` | 3.7704579282796167 | SOURCED |
| C11 | 2.226 / 2.2261 ORACLE perfect prior (gamma = 1) | `s24/results/priorladder.json :: rows[*]/MASS1.0` (mean over 126 rows; TILT1.0 and DIRAC1.0 identical; MASSFIXW1.0 2.233186) | 2.226080 (mean of the stored per-row values) | DERIVED |
| C12 | 1.313 whole-library best | `bench_results/recon_library_saturation.json :: universe_best/small` | 1.3134468768690186 | SOURCED |
| C13 | 1.594 ORACLE torsion-space ceiling (k = 4) | `s13/results/ceiling_report.json :: k/4/descent/mean`; `s13/results/adv_ceiling_k4_report.json :: arms/real/summary/mean` | 1.5940389994797128 | SOURCED |
| C14 | 0.611 ORACLE distance geometry from true distances | `s15/results/distgeo.json :: arms/ORACLE_true_distances/mean` | 0.6108221854859448 | SOURCED |
| C15 | -2.15 A per unit gamma / -2.1496 | `s24/results/priorladder.json`: (mean MASS0.1 - mean MASS0.0) / 0.1 = (2.833382 - 3.048338) / 0.1 | -2.149562 (first-rung slope) | DERIVED (no leaf stores the slope) |
| C16 | gamma = 0.0225 reaches 3.0 A | same rows, linear interpolation on the first rung: 0.0 + (3.048338 - 3.0) / (3.048338 - 2.833382) x 0.1 | 0.022487 | DERIVED |
| C17 | 68% common-mode / 0.676 | `s23/results/errdecomp.json :: rows[*]/f_common` (mean over 126; median 0.673693) | 0.675770 | DERIVED (per-row leaves; the mean is not stored) |
| C18 | 0.902 nats (trained state to its Gibbs optimum at T = 0.3) | `s25/results/q_gibbs.json :: results/divergence_ladder/1/ladder/"50 Adam steps (DEPLOYED)"/kl_nats` | 0.9019158325027764 | SOURCED |
| C19 | 0.24x MDE (circuit minus exact Boltzmann at T = 0.3) | `s25/results/q_alpha.json :: results/expressivity_vs_boltzmann/"circuit - exact Boltzmann @T=0.3"` (mean, mde, eff_over_mde) | -0.03024176839325163; 0.12605272002164683; 0.2399136519073786 | SOURCED |
| C20 | 1.18% worst deviation from the rank ladder | `s25/results/q_gibbs.json :: results/spectrum_target_independence/{worst, E_range, worst_frac_of_range}` | 0.04063087303875024; 3.4371432161567106; 0.011821117271971639 | SOURCED (also reproduced here at 0.394% on 1S9Z and 9KAR) |
| C21 | 5.6e-17 statevector vs independent simulator | `s25/results/q_verify.json :: results/"max \|p_statevector - p_dense\|"` (stored as a string) | "5.551e-17" | SOURCED |
| C22 | cos 1.000000000 parameter shift vs finite differences | `q_verify.json :: results/"cos(param-shift, exact FD), min over 12"` and `"cos(free-energy grad, FD), min over 4 (a,T)"` | "1.000000000"; "1.000000000" | SOURCED |
| C23 | 4.6e-10 relative error | `q_verify.json :: results/"rel \|g_ps - g_fd\| / \|g_fd\|, max over 12"`; `"rel err of free-energy grad, max over 4"` | "4.597e-10"; "4.663e-10" | SOURCED |
| C24 | 0.656 / 0.655634 / 0.524 / 0.566586 (tail-only baseline cosines) | 0.566586: `q_verify.json :: results/"cos(exact-expectation, TAIL baseline)"` with `"\|g_tail\| / \|g_exact\|"` = "0.519" (S25 re-verification, n = 7, layers 3, alpha 0.15). 0.655634: `core/quantum.py:42` (docstring, S9 audit at 10 qubits, 0.758x norm) and `tests/test_quantum.py:60 REC_TAIL_COS`, reproduced live by `test_cvar_gradient_baseline_is_constant_not_tail_only` (passed, job `pytest_core`). 0.524: `docs/FINDINGS.md:3161` (S8-9 sampled estimator); `s8/integrate_cvarcheck.json` is not on disk | "0.566586"; 0.655634 (code constant); 0.524 none | 0.566586 SOURCED; 0.655634 SOURCED BY TEST; **0.524 UNSOURCED (document-only)** |
| C25 | -0.649 / -0.252 / -0.047 / -0.311 log2 Var per qubit | `s25/results/q_plateau.json :: results/{linear_alpha1_T0, cvar_alpha025_T0, cvar_alpha01_T0, deployed_a1_T03}/log2_slope_per_qubit` | -0.6491874127317976; -0.2522012761595383; -0.047158929969141276; -0.31053500868274697 | SOURCED |
| C26 | 36.1 deg vs 36.4 deg (phi MAE, full sequence context vs sequence-blind marginal) | cited `s13/SPRINT13_DOSSIER.md:261-263`, `s13/tors_FINDINGS.md:373, 385`, `s24/LEDGER.md:550`; no JSON leaf in `s13/results/tors_*.json` holds 36.1 or 36.4 under a phi/mae/marg key; `s13/cache/tors_rows.npz` was not opened (npz) | none found | **UNSOURCED (document-only)** |
| C27 | +0.0004 (tuning) / +0.0030 (benchmark) identity-leak prices | cited `docs/FINDINGS.md:231, 4573, 5026` (S10-4 table: synthesis fit 3.2005 to 3.2009, +0.0004 [-0.0004, +0.0013]; benchmark +0.0030 [-0.0002, +0.0061]); the cited artefacts `s10/idaudit_*.json` (`docs/FINDINGS.md:4948`) do not exist on disk (no `s10/` directory) | none found | **UNSOURCED (artefact absent)** |
| C28 | 2/60 and 4/126 self-copies | 4/126: `tests/test_data.py:72-74` (`ok_canon >= 4`, passed), `docs/FINDINGS.md:2660`; both counts: `s26/results/i_identity_audit.json` (lane I, L15: the four dev targets 1CEK, 2FBU, 2P5H, 6B9K with shorter identity 1.0; benchmark 2/60 as a count only) | counts | SOURCED (4/126 by test and audit; 2/60 count-only, unquantified by design) |
| C29 | seven-configuration suite 3.058 / 3.132 / 3.215 / 3.253 / 3.674 / 3.755 / 3.881 | `s25/results/phys_suite.json :: configs_rank/{3,5,4,7,6,1,2}/mean` and `results/summary/leaderboard.json :: rows[*]/mean_secondary` | 3.058016219425174; 3.131686249878013; 3.215116303141773; 3.2530157270854128; 3.6742023446060776; 3.7553415439547573; 3.880496706249537 | SOURCED, on the POINT-CLOUD basis (see section I.2: the built-chain means of the same rows are 3.2187 / 3.3100 / 3.3732 / 3.4221 / 3.8248 / 3.8844 / 4.1015) |
| C30 | Legacy +0.330 and AMBER +0.455 vs a random 75-subset, 5/5 folds | `phys_suite.json :: controls_rank/1/vs_random/{effect, mde, effect_over_mde, ci95_fold, folds_same_sign}`; `controls_rank/2/vs_random/...` | 0.3302158940209069, 0.1662403218173309, 1.9863766528541527, [0.23753999909000342, 0.4158172931555582], 5; 0.45537105631568636, 0.18834102878507575, 2.4178006207841753, [0.328398404491368, 0.5355771582778813], 5 | SOURCED (point-cloud basis) |
| C31 | rho -0.7423 (readout entropy vs mean RMSD over the 9 VQE cells) | `s25/results/q_alpha.json :: results/cell_correlations/H` (alpha 0.26999350922052445, T -0.02344255091539687) | -0.7422836835050723 | SOURCED |
| C32 | 2,592 cells, 0 violations | `q_verify.json :: results/n_cells`; `results/"SUBSET-hood violations  <- THE THEOREM"`; `"full-support cells with EXACT value equality"` | 2592; 0; "972 / 972" | SOURCED |
| C33 | 78-89% of the free-energy gap closed | `q_gibbs.json :: results/training_control/{0,1,2}/frac_gap_closed` (T = 0.1, 0.3, 1.0) | 0.8914955163069174; 0.782951322968052; 0.7833283350457241 | SOURCED |
| C34 | mean \|z_moment\| per channel 0.7529 / 0.8013 / 0.1127 | cited `s25/LEDGER.md:1097`, `s25/agentPHYS_FINDINGS.md:115`; no leaf in `s25/results/*.json` rounds to any of the three under a moment key (`phys_suite.json` rows carry no z; `phys_landscape.json` carries `*_frac_absz_lt_0p1` only) | none found | **UNSOURCED (document-only)** |
| C35 | 355 passed / 13 skipped | cited `docs/STATE_BRIEF_2026-09-12.md:274` only (a live run with no junit file in the tree); reconciliation: the non-AMBER suite is 350 tests (337 passed, 13 skipped) plus `test_amber.py` 15 plus `test_amber_frame_invariance.py` 3 = 368 = 355 + 13; the governed S26 run (`s26/TEST_RUN.md`) is 369 tests, 356 passed, 13 skipped after lane I added one test | none for the S25 figure | **UNSOURCED (document-only)**; superseded by `s26/results/test_run.json` |

Count: 4 claims UNSOURCED (C26, C27, C34, C35) and one component of C24 (0.524). C06 is
named to a file I did not open. The rest are sourced to an artefact leaf, a passing test, or
a stated derivation from stored rows. Ledger entry: `s26/LEDGER.md` L28 (lane E).

Basis notes that matter for the presenter: C08, C29 and C30 are point-cloud numbers (the
physics suite and the leaderboard's secondary column). The state brief lines 57 and 190-191 and
`results/summary/professor_brief.md:115-125` quote them next to the built-chain 3.2148 without
saying so.

---

## D. TEST SUITE

Not run by this lane (contract section 3). Lane I's governed run is the record:
`s26/TEST_RUN.md` and `s26/results/test_run.json` (rendered 2026-09-13 00:41 by
`s26/i_test_report.py`), junit files `s26/results/pytest_core.xml`, `pytest_amber.xml`,
`pytest_amber_frame.xml`, `pytest_nanpoison_post_core_edit.xml`.

- `pytest_core` (tree `a4db170c`, non-AMBER files): 350 tests, 337 passed, 0 failed, 0 errors,
  13 skipped, none from the memory guard; wall 195.3 s, peak RSS 1.694 GB. The 13 skips are 11
  `VERIFY_SLOW=1` gates (3 in `test_equivalence.py`, 8 in `test_integration.py`) and 2 absent
  artefacts (`bench_results/optimised_tuning126_w6.json`, `smoke_opt_test.json`).
- `pytest_amber` (tree `37bddbbb`): 16 passed (15 + lane I's new memory-guard message test),
  260.5 s, 0.872 GB. `pytest_amber_frame`: 3 passed, 255.7 s, 0.324 GB.
- Combined: 369 tests, 356 passed, 0 failed, 0 errors, 13 skipped.

Read here for section C: the two benchmark-record tests
(`tests/test_pipeline.py::test_the_committed_benchmark_report_still_holds_its_reference_numbers`,
`tests/test_equivalence.py::test_benchmark60_cached_record_reproduces_and_is_not_rerun`) are in
the passed set, as are the NaN-poison guard
(`tests/test_pipeline.py::test_nan_poisoning_the_native_changes_nothing_deployable`), the
production-config pin (`test_the_default_configuration_is_the_preregistration`), the four-component
reproduction of `s8/integrate_vqe.json` (`rmsd_vqe_sel` 3.3135, `medoid128` 3.3443), and every
`core.quantum` equivalence test including the tail-baseline regression test.

---

## E. DECLARED DEFECTS, restated

1. **Identity is normalised by the longer sequence.** `core/data.py:182` `identity(a, b, gap,
   norm="longer")` returns `_nw_matches(a, [b])[0] / max(n, m)`; `peptide_db.py:119-170`
   (`shared / max(len(a), len(b))` at line 170) is the legacy copy. A 13-residue peptide sitting
   verbatim inside a 25-residue member scores 0.520 and passes the 0.6 filter
   (`tests/test_data.py:194-211`, 1CEK in 1A11). `containment` (`core/data.py:216`) divides by
   the shorter length and reports 1.0. Consequence: 4 of 126 dev targets and 2 of 60 benchmark
   targets have a verbatim copy in their own fold model's training set (section C, C28). Lane I
   added `norm="shorter"` behind the default (commit `37bddbbb`, ledger L15); nothing on the
   production path passes it, and the pinned clusters and folds are untouched (hashes in F).
2. **The projection builds a constant 3.804 A virtual bond and only trans peptide bonds.**
   `core/project.py:149-155` (`BOND_N_CA` 1.458, `BOND_CA_C` 1.525, `BOND_C_N` 1.329, angles
   111.0 / 116.2 / 121.7 deg, `OMEGA_TRANS = pi`), identical constants in
   `core/geometry.py:57-65`; `frames` takes a scalar omega. Measured on all 126 production
   chains: `s26/results/e_reproduce.json :: summary/virtual_bond/chain_mean_bond` =
   3.803954938363982 with `chain_bond_sd_max` 1.802075856660821e-15, against a native mean of
   3.812176007269521. A cis bond (CA-CA about 2.9 A) is not representable. Cost on the
   instrument: unmeasured; lane PH's `s26/PREREG_cis.md` predicts zero cis natives on model 1
   because `core/data.py` drops any peptide with a step outside 3.5 to 4.1 A.
3. **`pool_gate = WARN` with 2 per-target violations.** `s25/resultslab/exportlib.py:598-639`
   `pool_oracle_gate`: the aggregate (mean below the mean ORACLE pool best 1.7108) is hard and
   raises; a per-target violation only sets WARN. `results/summary/leaderboard.json` production
   row: `pool_gate` "WARN", `n_pool_violations` 2, `pool_margin` 1.501796830456039. The two
   targets are 1D6X and 1KWE (lane I, L7 corrected in L9): the emitted chain is a coordinate
   average projected to ideal geometry, not a pool member, so on a single target it can land
   nearer the native than any window. `schema.py:309` applies the gate; lane I's `eb89c165`
   makes `fmt_leaderboard` name the WARN rows.
4. **The fold models are not independent.** `core/predict.py:224` `_fold_fragments(fold,
   n_folds)` and `:345` `train_fold`: the training set of fold f is every peptide whose pinned
   fold is not f plus the identity-filtered fragments, so each dev native is a training label
   for four of the five fold models (S24 L8: each model saw about 100 of the other 125 dev
   natives). The five files `distogram_models/fold{0..4}_esm_frag.pt` (hashes in F) therefore
   share most of their training data, and the per-target errors of different folds are not
   independent draws; the fold-clustered CI of `s24.stats_lib.compare` is the deciding interval.
5. **Moment standardisation is non-monotone on AMBER energies.** `s25/phys_lib.py:98` `zmoment`
   (an `(x - mean) / sd` on raw energies whose range spans 1e28) can reorder candidates under
   float64; production uses ranks: `core/pipeline.py:788` `_zrank` (`scipy.stats.rankdata`,
   then standardise), applied at line 838 inside `quantum_stage` to the distogram score of the
   top 2^n candidates, and to nothing else. Confirmed by the trace (E first/last 5 above) and by
   lane I's audit (L8): no production path applies a moment z-score to an AMBER energy; the
   production `Config` has `quantum=False`, so in the 3.2148 arm the selector does not run.

Two further observations recorded for the Adversary (not declared defects): the deployed
`VQE_LFO` table runs alpha = 1.0 on folds 0, 3 and 4 (61.9% of targets, no tail constraint;
section B); and one production relaxation (9KAR) ends above `core.amber.CONVERGE_MAX_KCAL` while a
second (2BP4) ends above the S8 bond+angle strain rule (section C, C03); both are counted in
`rmsd_full`.

---

## F. PINNED-FILE HASH TABLE

`s26/e_hashes.py` writes `s26/results/pinned_hashes.json` (written 2026-09-13T00:09:00) and
`--check` reports drift. Bytes only; `results/benchmark_manifest.json` is hashed as a byte
stream and never parsed; `esm_cache.npz` is stat-only (1,509,204,583 bytes, never opened).

| file | sha256 (first 16) | bytes |
|---|---|---:|
| `results/benchmark_manifest.json` | a40581ad01cfd2b7 (equals the S20 record a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d) | 8002 |
| `peptide_folds.json` | 89bc58edc59c6d79 | 18956 |
| `peptide_clusters.json` | 9ae3d8f3afc281df | 20217 |
| `catrace_prior.npz` | d93a5c627d5ca4df | 7172 |
| `results/monomer_manifest.json` | 5a10cf4eac3d095c | 28719 |
| `esm_small.npz` | 5396757c2b119c61 | 22062826 |
| `esm_pca.npz` | a0e2b62fb99d5118 | 169950 |
| `peptide_db.npz` | 3e0b447284e73d1b | 406080 |
| `fragment_db.npz` | 678dce94558460e1 | 1874081 |
| `fragment_db_large.npz` | 10e9feb16814276d | 4082944 |
| `distogram_models/fold0_esm_frag.pt` (production) | 0ce6625ef3fd541b | 1497455 |
| `distogram_models/fold1_esm_frag.pt` (production) | 23691cc531435766 | 1497519 |
| `distogram_models/fold2_esm_frag.pt` (production) | 57214b7a209ff986 | 1497455 |
| `distogram_models/fold3_esm_frag.pt` (production) | 25ba58715917945b | 1497519 |
| `distogram_models/fold4_esm_frag.pt` (production) | 4563414b9874f496 | 1497455 |
| the other 25 files of `distogram_models/` (`fold*_esm.pt`, `*.d20`, `*_s1.pt.d20`, `*_sw-lin.pt`, `*_sw-none.pt`, `fold1_seq.pt`) | in the JSON | 30 files, 44.7 MB |
| `pdbs/` | 61 files, each hashed in the JSON | 31.0 MB |

Verdict recorded in the JSON: `benchmark_manifest_matches_S20_record` = true. Lane I's identity
audit (L15) re-hashed `peptide_folds.json` and `peptide_clusters.json` before and after its
in-memory re-clustering and found them equal to this table.

---

## G. NOT IN GIT

Measured with `git ls-files` and `git check-ignore` on 2026-09-13 (sizes from the file system).
"Breaks" is read from the module that consumes the artefact, not guessed.

| artefact | size / files | tracked | what breaks without it |
|---|---|---|---|
| `bench_results/cache/1fc9f2dcf489e2fb/` (the 126 production records) | 0.7 MB / 126 | no (gitignored) | `s12.instrument.shipped_record`, `s25/resultslab/providers.incumbent`, every S12 to S26 reproduction including section H; `core.pipeline.run` would rebuild it only with AMBER (about 300 s of CPU per the harness) and only if every input below is present. This is the one artefact the whole record stands on that is not versioned. |
| `esm_cache.npz` | 1509.2 MB / 1 | no | `core.data.esm_raw` for any sequence outside `esm_small.npz` (`guard_esm` refuses to load it without a probe); training new distogram inputs (lane P uses `s5/esmraw.npz` instead) |
| `esm_small.npz` | 22.1 MB / 1 | no | the hot ESM cache (360 sequences, all 126 dev targets): without it every `run_target` falls through to the 1.5 GB bank; `tests/test_data.py::test_distogram_predictions_are_identical_end_to_end` skips |
| `esm_pca.npz` | 0.2 MB / 1 | no | `core.data.esm_embed` (the 32-d projection): the 183-column feature stack of every production distogram cannot be built; `run_target` fails on all targets |
| `peptide_db.npz` | 0.4 MB / 1 | no | the database cache; rebuilt by `core.data` scanning `pdbs/` and `pdbs_ext/` |
| `pdbs_ext/` | 578.2 MB / 1463 | no | the database rebuild above (`pdbs/` alone holds 61 files); without it `peptide_db.npz` cannot be regenerated and the retrieval library, natives and folds are unrecoverable |
| `prots/` | 5842.9 MB / 13751 | no | `fragment_db.py` (`fragment_db.npz`, `fragment_db_large.npz`) and `s8/invfold.py` `stage_scan` |
| `fragment_db.npz`, `fragment_db_large.npz` | 1.9 + 4.1 MB | no | `_fold_fragments`: the retrieval universe loses its fragment windows (7193 windows for 1S9Z include them), pools change, fold models cannot be retrained |
| `catrace_prior.npz` | 7 KB | no (pinned) | `core.geometry.CATracePrior` re-derives it on first use; pinned, so regeneration is forbidden |
| `distogram_models/` | 44.7 MB / 30 | no (pinned) | `core.pipeline.fold_model` and `guard_esm` raise; no scoring at all |
| `s8/generate_univ/` | 117.4 MB / 126 | no | the instrument itself: `s12.instrument.targets()` (target order, T-numbers), `load_univ`, `native_ca` for the results lab, `pool_best`, every S12 to S26 script; `results/summary/target_map.json` raises on a changed order |
| `s7/debias_cache/` | 58.9 MB / 131 | no | `s7.debias._load`, `s8.generate.stage_univ`'s reproduction assertion, `s8.invfold.stage_verify` |
| `s8/consensus2_cache/` | 132.7 MB / 252 | no | `s8.consensus2`, `s8.project` (the S8 instrument), `tests/test_geometry.py::test_pool_best_and_shipped_selected_reproduce_the_instrument` (skips) |
| `s8/integrate_dmat/` | 47.0 MB / 126 | no | `s8.integrate.load_dmat` (consensus medoid over the 500 x 500 matrix) |
| `s8/inband_cache/`, `s8/integrate_chan/`, `s8/integrate_amber/`, `s7/audit_cache/`, `s7/repr_models/`, `s8/invfold_scan.npz`, `s8/invfold_windows.npz` | ABSENT | no | `s8.inband.stage_*`, `s8.integrate.stage_chan/corr/fuse/...`, `s8.inband.stage_devpool` (asserts against `s7/audit_cache`), `s8.inband.stage_ens`; S7-11's per-target artefact `s7/repr_tune.json` is lost (L11) |
| `s8/invfold_models/` | 8.4 MB / 3 | no | `s8.invfold.load_model` (the `iv_*` channels of `s8.integrate`) |
| `s5/esm32.npz`, `s5/esmraw.npz`, `s7/repr_cache/` | 11.4 + 243.5 + 7.7 MB | no | `s7.repr_select`, lane P's ladder (`s26/p_ladder.py` reads only these compact tables, L11) |
| `s9/final_cache/` (sealed benchmark pools) | 56.5 MB / 60 | no | `tests/test_pipeline.py` sections 3 and 4 (bit-for-bit reproduction of the committed benchmark synthesis) skip |
| `s9/synth_cache/` | 5.9 MB / 1386 (126 tracked: `full_sc_75_*.json`) | partly | the tracked 126 are the S9 reference `tests/test_equivalence.py` and `verify/headline_audit.py` read; the rest are probe arms |
| `s8/generate_rama.npz`, `s8/inband_ss.npz`, `s8/integrate_native.npz` | < 0.2 MB | no | `s8.generate.rama_logp`, `s8.inband.ss_logp`, `s8.integrate.load_native` |
| `s12/results/` | 56.9 MB / 466 (369 tracked) | partly | untracked entries are `.npz`/`.log` by the global ignore |
| `s24/results/` | 6.7 MB / 50 (36 tracked) | partly | `e_ceff_obs.npy`, `e_q_obs.npy` and logs untracked |
| `s25/results/` | 8.2 MB / 158 (142 tracked) | partly | logs untracked; `suite_cells/` tracked |
| `score_weights.json` | ABSENT by design | n/a | shell weights are 1 (`s8/inband.py:238`, `core.predict`); its absence is the production configuration |
| `s26/models/` | empty | n/a | reserved |

Tracked and pinned: `peptide_folds.json`, `peptide_clusters.json`, `results/benchmark_manifest.json`,
`results/monomer_manifest.json`, `pdbs/` (61), `s8/integrate_vqe.json`, `s8/project_prior.json`,
`s8/project_devarm.json`, `s9/final_report.json`, `s9/final_synth.json`, `results/summary/` (7),
`results/structures/` (2142), `verify/project_{equiv,inputs,exactness}.json`, `bench_results/*.json`.

---

## H. REPRODUCTION (Part 2.3)

`s26/e_reproduce.py` under the governor (`s26/jobs_done/e_reproduce.json`: exit 0, wall 5.0 s;
the 5-second RSS sampler read 0.004 GB because the job finished inside one sampling interval;
the in-process peak, `pl.proc_rss()`, was 46.09 MB). Artefact `s26/results/e_reproduce.json`
(`complete` true, 126 rows, provenance module `e_reproduce.py`, source sha 4bd3111ea1007fb7,
git 74e073e26d5c dirty). Method: every stored structure of the 126 records (`avg_ca`, `fit_ca`,
`ca`, `amber_ca`) re-scored against its native through `s12.instrument.ca_rmsd`, the natives
from `core.backend("data").load()` put through `core.pipeline._q` (the reference's float32
round trip, `Config.reference_precision=True`). Tolerance 2.6e-04 A.

| quantity | fresh mean | stored mean | max per-target abs disagreement | n |
|---|---|---|---|---|
| `rmsd_avg` (point cloud) | 3.048338093879532 | 3.048338093879532 | 0.0 | 126 |
| `rmsd_fit` (lam = 0 chain) | 3.2040761603809194 | 3.2040761603809194 | 0.0 | 126 |
| `rmsd_arm` (built chain, production) | 3.2147651542109985 | 3.2147651542109985 | 0.0 | 126 |
| `rmsd_full` (AMBER emission) | 3.2354598538973844 | 3.2354598538973844 | 0.0 | 126 |

T030 = 1S9Z: `rmsd_arm` fresh 0.18198112330908295, stored 0.18198112330908295 (avg 0.195888,
fit 0.181943, full 0.317379). Stored aggregates re-read from the records: shipped
3.4540004952559396, pool_best 1.7108244199364904, top_m_best 2.3061526409453816. Fold labels of
all 126 records equal the pinned folds (`fold_mismatch_vs_pinned_folds` = []); every record
carries `cfg_key` 1fc9f2dcf489e2fb and `dev_mode` ""; no `amber_err`.

Fresh single-target runs (`e_trace.py`, section B): `pl.run_target(target, fold, Config(amber=False))`
followed by `pl.label` reproduces shipped, pool_best, top_m_best, rmsd_avg, rmsd_fit and rmsd_arm
to 0.0 on 1S9Z and on 9KAR, with the filtered set `sub` identical and the emitted `ca` identical
to the record (max abs diff 0.0). AMBER was not re-run (the record's `amber_ca` is scored instead).

Secondary measurements from the same rows (`summary/virtual_bond`, `summary/projection_gap_arm_minus_avg`):
point cloud mean virtual bond 2.9613938229671057 A against native 3.812176007269521 (22.317%
contracted), worst single bond 0.6491799200244281 A, 80 of 126 targets under 3.4 A; built chain
3.803954938363982 (sd at most 1.8e-15); AMBER emission 3.8665519730270517. Projection gap
(arm minus avg): mean +0.1664270603314659, median +0.0977, IQR [+0.0049, +0.3379], min -0.3483,
max +0.6529, chain better on 16 of 126. These are the S25 L4/L8/L11 figures.

The 3.2126 versus 3.2148 pair: `results/summary/leaderboard.json` production row `mean`
3.2126212503925293 (basis `built_chain_bb`, `mean_secondary` 3.048329101217154) against the
cache's 3.2147651542109985. The leaderboard re-projects the stored point cloud through
`exportlib.chain_from_ca` and round-trips through PDB `%8.3f` quantisation; lane P measured the
per-target gap at up to 0.171 A (L14). The two are the same arm on two bases of emission and
are not a discrepancy; presenter documents should quote one and name it.

Verdict: exact. No reproduction failure; nothing to ledger under Part 2.3. The coordinator
recorded the same in L16b.

---

## I. Examiner's observations for the Adversary

1. **Production-path edits during Phase 0.** Between my traces (git `74e073e2`, dirty) and this
   assembly, lane I committed `37bddbbb`: `core/data.py` (identity `norm` keyword, dead helper
   removed), unused-import removals in `core/amber.py`, `core/bench.py`, `core/cache.py`,
   `core/predict.py`, `verify/vqe_lfo_audit.py`. My AST check (docstrings stripped) found all
   seven working files non-identical to their parent at 00:40; lane I's own `i_ast_check.py`
   verified each diff is exactly the removed names or the new keyword, the changes are ledgered
   (L10, L15), and `tests/test_pipeline.py` + `test_data.py` + `test_amber.py` passed on the
   edited tree under the governor. My section H numbers were produced on the tree before those
   edits and the identity path they touch is not on `run_target`. The Adversary should still
   re-run `python s26/e_reproduce.py` on HEAD; it costs 5 s.
2. **The seven-configuration suite is quoted on the point-cloud basis** next to the built-chain
   production number in `docs/STATE_BRIEF_2026-09-12.md:57, 190-191` and
   `results/summary/professor_brief.md:115-125`. Built-chain means for the same configurations
   (`results/summary/leaderboard.json :: rows[*]/mean`): distogram 3.218708055115966,
   amber_distogram 3.309961759812365, legacy_distogram 3.3731520305329474,
   legacy_amber_distogram 3.4220825100678876, legacy_amber 3.8248241126995124, legacy
   3.8843706034949403, amber 4.101517931701575; the random-75 null (3.4251) has no built-chain
   counterpart on disk. The ranking is unchanged; the magnitudes differ by 0.16 to 0.22 A.
3. **The deployed selector is unconstrained on 78 of 126 targets** (alpha = 1.0 on folds 0, 3,
   4). On those targets "CVaR tail" means the whole register and the readout is a free-energy
   weighting only (section B, 9KAR).
4. **One unconverged AMBER relaxation (9KAR) and one strain-flagged relaxation (2BP4) sit inside
   `rmsd_full`** (C03; lane PH's L24 counts 125 of 126 converged by the e1 gate). The record
   carries `converged` flags only since S16; the production cache predates the gate.
5. **The production cache is not versioned** (section G, first row).
6. **Four presenter-facing numbers have no artefact** (C26, C27, C34, C35) and one benchmark
   number is named to a file this lane did not open (C06).
