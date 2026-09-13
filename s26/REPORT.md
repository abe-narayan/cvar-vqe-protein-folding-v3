# THE PEPTIDE STRUCTURE PROGRAMME, SPRINTS 5 TO 26: THE REPORT

Report writer: lane E, Sprint 26. Built incrementally on branch `s26` from 2026-09-13; every
section is committed as it lands and every number carries an artefact path (Appendix B). The
report builds on, and cites rather than replaces, `docs/FINDINGS.md` (the S5 to S11 record with
its corrections ledger), `docs/CONDENSED_REPORT.md` (S5 to S13), the sprint ledgers `s12/` to
`s25/`, `ARCHITECTURE.md` (the frozen specification), `s25/QUANTUM.md`, `s26/EXAMINATION.md`,
`s26/PH_PART_IV_NOTES.md` and the S26 ledger. A reader who knows some quantum computing and
nothing about proteins should meet no undefined term; Appendix A holds the glossary.

Conventions. Distances are in Angstrom (A, 1e-10 m). Every RMSD contrast names its basis on
both sides: POINT CLOUD (the raw coordinate average, `rmsd_avg`), BUILT CHAIN (the ideal-geometry
chain the pipeline emits, `rmsd_arm`, the production result), RELAXED CHAIN (after the AMBER
step, `rmsd_full`), or SINGLE WINDOW (one retrieved candidate scored on its own). An arm labelled
ORACLE reads the native structure and is a diagnostic, never a result. "Closed" means the
direction was measured to its own ceiling or refuted by a pre-registered falsifier; "open" means
neither has happened. Test names are given as `file::test`; the record of their last run is
`s26/TEST_RUN.md`.

Status of the sections: Parts I to VI, IX and Appendix A are drafted from the record as it
stands; Parts VII and VIII and Appendix C carry slots that fill as S26 verdicts land; Appendix B
lists the numbers of the drafted parts and grows with them.

---

## PART I. THE PROBLEM

### I.1 What is being predicted

A protein is a chain of amino acids (residues). Each residue contributes one alpha carbon (CA)
to the backbone, so a peptide of n residues has a CA trace of n points in space. The chain's
shape is set by two rotatable angles per residue, phi and psi (the third backbone angle, omega,
is almost always fixed near 180 degrees, the trans configuration), and by the fixed geometry of
the covalent bonds, which places consecutive CA atoms about 3.8 A apart. Given the amino-acid
sequence, the problem is to predict the CA trace.

This programme restricts itself to short peptides, 9 to 16 residues long, taken from the
Protein Data Bank (deposited experimental structures, most of them solved by NMR and carrying
several models; the instrument scores against model 1). Two properties of that regime decide
what the work can and cannot do.

First, the answer is scored by CA-RMSD after superposition: the predicted CA trace is rotated
and translated to fit the deposited one as well as possible (the Kabsch algorithm, reflections
forbidden), and the root-mean-square distance between corresponding CA atoms is reported
(`s12/instrument.py:91` `ca_rmsd`). A mirror image is not a match. Any two 12-residue chains
superpose to within a few Angstrom, so the scale of the problem is small: a constant
alpha-helix, which knows nothing about the target, scores 4.0648 A built chain on this
instrument (`s12/results/s14_ladder.json :: rows/L0_constant_helix/mean`), and the production
pipeline scores 3.2148 A built chain.

Second, at this length the sequence says little about the shape. The programme measured the
sequence-to-torsion channel directly: a model that reads the whole 15-residue context predicts
phi no better than a model that reads no sequence at all (the record quotes 36.1 against 36.4
degrees of mean absolute error, `docs/CONDENSED_REPORT.md:104-116`, `s13/SPRINT13_DOSSIER.md`
section 5). The artefact behind those two numbers is not on disk (`s26/EXAMINATION.md` C26),
so they are cited here to the documents that carry them and are not in Appendix B. Short
peptides in isolation are also often flexible, and the deposited model is one member of an
ensemble. Both facts matter for what "accuracy" can mean.

### I.2 Why this is a good testbed

The peptide regime is small enough that the whole candidate space can sometimes be enumerated
(Sprint 21 enumerated every configuration of a 2^n discrete latent space on all 126 targets,
`s21/LEDGER.md` L14, L17), small enough that a 7-qubit register can index a candidate set, and
large enough that retrieval from a library of real protein fragments produces candidates within
1 to 2 A of the answer on most targets. It is therefore a place where a variational quantum
algorithm can be put next to its classical ceiling and its classical controls with nothing
hidden: the exact optimum of every objective can be computed, the ORACLE ceiling of every stage
can be measured against the deposited structure, and every claimed improvement can be paired
against a matched control on the same 126 targets.

### I.3 The shape of the answer, stated first

The programme did not reach its accuracy targets. On the 60-target sealed benchmark, opened
once, the synthesis architecture scored 2.9610 A against the shipped baseline's 2.9507 A (both
single-window selection arms on the S9 instrument), a paired difference of +0.0103 A with 95% CI
[-0.1596, +0.1803], 31 wins to 29 losses (`docs/FINDINGS.md:4479-4607`, S9-10; the two means are
asserted within 5e-4 A of the committed record `s9/final_report.json` by
`tests/test_pipeline.py::test_the_committed_benchmark_report_still_holds_its_reference_numbers`,
passed in `s26/TEST_RUN.md`; this lane did not open the record). On the 126-target development
instrument the production result is 3.2148 A built chain
(`bench_results/baseline_tuning126.json :: science/rmsd_arm/mean`), against an ORACLE best pool
member of 1.7108 A single window (`science/pool_best/mean`). The findings, which are the actual
output, are what Parts IV to VIII are about.

---

## PART II. THE INSTRUMENT

This part is written so that a reader could design a valid experiment on this repository
without further help. The rules are the ones the programme learned by breaking them; the
retraction count that taught them is 23 corrections in `docs/FINDINGS.md:60-118` plus five S25
claims retracted by internal audit (`s25/LEDGER.md` L5, L7, L9, L11, L15).

### II.1 The targets

- **tuning126.** 126 peptides of 9 to 16 residues (`core.pipeline.manifest("tuning126")`, order
  pinned by PDB id, `tests/test_pipeline.py::test_manifest_order_is_pinned`). Every number in
  Sprints 12 to 26 is measured here. The window universe of each target (every candidate window
  the library can supply, with the target's native CA trace as a reporting label) is cached in
  `s8/generate_univ/<pdb>.npz`; `s12.instrument.targets()` is the sorted list of those files,
  and the T-numbers of the results lab (T001 to T126) follow that order
  (`results/summary/target_map.json`).
- **Five pinned folds.** Every sequence in the peptide database is assigned to an identity
  cluster (470 clusters, `peptide_clusters.json`, asserted by `tests/test_data.py:305`) and
  every cluster to one of five folds (`peptide_folds.json`). Fold f's distogram model is trained
  on the other four folds plus the identity-filtered protein fragments (`core/predict.py`
  `train_fold`, `_fold_fragments`), so no target is scored by a model that saw its own sequence.
  The two files are write-on-first-use caches: if either is deleted it is silently re-derived,
  and doing that once moved 13 benchmark targets and invalidated every model (`README.md`, "The
  trap worth knowing about"). Their hashes are in `s26/results/pinned_hashes.json`.
- **dev24.** A 24-target development split, cluster-disjoint from the benchmark (`README.md`,
  "The instrument discipline"), used by Sprints 8 and 9 for single pre-registered passes.
- **benchmark60, sealed.** 60 targets in `results/benchmark_manifest.json`, cluster-disjoint
  from tuning126 (`tests/test_pipeline.py::test_tuning_instrument_is_disjoint_from_the_benchmark`).
  It was opened exactly once, by S9-10, on pre-registered constants
  (`docs/FINDINGS.md:4479-4607`); the CLI refuses to run it without
  `--i-am-spending-the-benchmark`
  (`tests/test_equivalence.py::test_the_benchmark_is_refused_without_the_explicit_flag`), and
  `tests/test_equivalence.py::test_benchmark60_cached_record_reproduces_and_is_not_rerun`
  asserts the cached record is reproduced rather than re-run. Its manifest is hashed as bytes
  and never parsed by the examination (sha256
  a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d, 8002 bytes, the S20 record,
  `s26/results/pinned_hashes.json`).
- **The declared leak.** Identity between two sequences is normalised by the longer one
  (`core/data.py` `identity`), so a peptide sitting verbatim inside a longer library member can
  score below the 0.6 threshold and land in another fold. Four of 126 dev targets and two of 60
  benchmark targets carry such a self-copy in their own fold model's training set
  (`s24/LEDGER.md` L4; `s26/results/i_identity_audit.json`, S26 ledger L15). The dev effect was
  priced in S10-4; the benchmark effect is declared and unquantified by design, because
  quantifying it would open the benchmark.

### II.2 The endpoint and the reporting bases

The endpoint is the mean over targets of the full-chain CA-RMSD to model 1 after superposition.
The production pipeline emits three objects per target and the record reports all three:

| basis | what it is | production mean, 126 targets | artefact |
|---|---|---|---|
| POINT CLOUD | the uniform coordinate average of the retained windows; not a molecule (mean virtual CA-CA bond 2.9614 A against a native 3.8122, 22.3% contracted, worst bond 0.649 A, 80 of 126 targets under 3.4 A) | 3.0483 | `bench_results/baseline_tuning126.json :: science/rmsd_avg/mean`; bonds `s26/results/e_reproduce.json :: summary/virtual_bond` |
| BUILT CHAIN | the ideal-geometry backbone nearest that cloud under a torsion prior; THE PRODUCTION RESULT | 3.2148 | `science/rmsd_arm/mean` |
| RELAXED CHAIN | the built chain after a restrained AMBER minimisation | 3.2355 | `science/rmsd_full/mean` |

The basis ruling (`s25/LEDGER.md` L4, L8, L11) is that the built chain is the result and the
point cloud is a labelled intermediate; the projection costs +0.1664 A on average (median
+0.0977, IQR [+0.0049, +0.3379]; the chain is better than the cloud on 16 of 126 targets;
`s26/results/e_reproduce.json :: summary/projection_gap_arm_minus_avg`). A contrast that mixes
bases is not a contrast. Two further built-chain numbers appear in the record and are the same
arm on two emission paths: 3.2148 from the production cache and 3.2126 from the results lab,
which re-projects the stored cloud and round-trips through PDB coordinate quantisation, a
per-target gap of up to 0.171 A (`results/summary/leaderboard.json :: rows[0]/mean`; S26 ledger
L14). A fourth number, 3.2041 A, is the chain fitted with no torsion prior (`rmsd_fit`,
`science/rmsd_fit/mean`) and is not a production arm.

The S25 ruling has a corollary the S26 examination had to restate (S26 ledger L28 item 3, L29
item 2): the seven-configuration physics suite, the random-75 null and the "worse than noise"
verdicts of Part IV are POINT-CLOUD numbers (`s25/results/phys_suite.json`), and the built-chain
means of the same rows are in `results/summary/leaderboard.json :: rows[*]/mean` (distogram
3.2187 to AMBER 4.1015). Wherever this report quotes one it names the basis.

### II.3 The statistics

Every comparison is paired per target through `s24.stats_lib.compare(a, b, folds, names)`:

- **Effect** is the mean of the per-target differences; the median is printed beside it.
- **SE** is the standard error of that paired difference (sample sd with ddof = 1 over the
  square root of n), never an arm's own SE
  (`tests/test_instrument.py::test_paired_se_is_the_one_the_mde_rule_consumes`).
- **MDE = 2.8016 x SE**, per comparison. It is the smallest effect a paired test at 5%
  two-sided size and 80% power resolves at this n (1.96 + 0.84 standard errors). A result must
  clear its own MDE; 0.7 to 1.3x MDE is the Type-M zone (the magnitude, if real, is inflated by
  selection) and is not a result; below 0.7x is "underpowered", not "null". The earlier
  "MDE = 0.084 A" was a single constant quoted as a property of the instrument and was wrong by
  up to 84x in both directions (`s21/LEDGER.md` L8). For scale: a built-chain contrast between
  the production pipeline and a synthetic quantum-selected variant of it has paired SE 0.0342 A
  and MDE 0.0958 A; the deployed selector against the plain argmin on the S8 selection
  instrument has MDE 0.2051 A (`s26/results/q_mde_reference.json`).
- **Two confidence intervals.** The iid bootstrap CI over targets, and the fold-clustered CI
  (targets grouped by their pinned fold, five clusters). The fold-clustered interval decides,
  because the five distogram models are not independent: each saw most of the other 125 dev
  natives (`s24/LEDGER.md` L8), so per-target errors are correlated within and across folds.
  Until S20 every CI in the record was iid over targets (`s20/LEDGER.md` L9).
- **Wins, losses, ties**, and the per-fold effects; a positive needs 5 of 5 folds in sign.
- **The concentration null.** A mean carried by a handful of targets is suspicious, but a raw
  drop-top-10 threshold is not a test (it misfired once); the mean without the ten largest
  contributions is compared against a uniform-effect null (`s24/stats_lib.compare`;
  `s25/QUANTUM.md` section 6.4 shows one such check landing at the 64.7th percentile of the
  null).
- **A verdict string** (measured, not measured, Type-M zone, underpowered, worse) printed with
  the numbers, never instead of them.

### II.4 Controls

- **Controls match the operator's space.** A displacement is compared to a random displacement
  of the same magnitude in the same coordinates; a reject filter is compared to a random reject
  of the same count; a selector is compared to a random subset of the same size. The programme
  recorded this rule after breaking it three times in two sprints (`s16/LEDGER.md` L27; S26
  `PREREG_c3_control.md`, `PREREG_amber_reject.md`).
- **The zero-information control is plausible, not uniform.** Uniform-on-the-torus torsions are
  a worse structure, not an uninformative one; the standing control is a constant alpha-helix
  (phi = -63, psi = -42 degrees), 4.0648 A built chain on tuning126
  (`s12/results/s14_ladder.json :: rows/L0_constant_helix/mean`).
- **Best-of-N, never an initialisation mean.** An optimiser is compared to the best of N draws
  from its own untrained state (`s25/results/q_gibbs.json :: results/training_control`,
  N = 200).
- **Tied argmins are averaged.** `np.argmin` returns the first index of a tie, and the candidate
  order is the retrieval order, which is not neutral; a tied signal once "won" by 1.386 A by
  reading that order (`docs/FINDINGS.md:2751`, "A methodological correction that changed the
  answer"; `s24.stats_lib.argmin_tied`). A selector's outcome is the mean RMSD over its tied
  argmin set.
- **Grid oracles are order statistics.** A per-target minimum over K arms is mostly best-of-K;
  `s24.stats_lib.best_of_k_within` reports the across-target null and `k_eff`, and only the
  split-half transfer (choose on one half of the pool, score on the other) counts
  (`s24/LEDGER.md` L15, L15-A).
- **Natives never enter a deployable operator.**
  `tests/test_pipeline.py::test_nan_poisoning_the_native_changes_nothing_deployable` NaN-poisons
  every native quantity and asserts the emitted coordinates are bit-identical.

### II.5 The discipline around a run

Every experiment since S15 is pre-registered (`s26/PREREG_<name>.md`: hypothesis, exact
falsifier, comparison arm, basis, expected effect against a computed MDE, memory, agent-hours)
before its first result exists and is never edited afterwards; addenda are appended. Artefacts
are written with `s24.stats_lib.save_atomic`, which stamps the module hash, git commit and dirty
flag and gates `complete: true` on a full key set over the expected number of rows. A positive
result is re-run at a second seed and with the fold processing order reversed and must land
inside its own CI (`s26/BRIEF.md`). Nothing superseded is deleted; retractions are appended.
Every emitted structure in the results lab reproduces its own RMSD through the instrument to
within PDB quantisation, worst 2.6e-4 A (`s25/LEDGER.md` L10).

---

## PART III. THE PIPELINE, STAGE BY STAGE

Frozen specification: `ARCHITECTURE.md` (S25). Code: `core.pipeline.run_target`.
Configuration: `core.pipeline.PROD` (K = 500, m = 75, penalty `ramah`, lam = 0.3, multi-start,
stable tie-break, AMBER k = 10, steps = 0, n_folds = 5, quantum off;
`tests/test_pipeline.py::test_the_default_configuration_is_the_preregistration`). The shapes
and numbers below are from `s26/results/e_trace_1S9Z.json` (T030, the 0.182 A target) and
`s26/results/e_trace_9KAR.json` (a 7.4 A target); every stage of both traces reproduces the
production record at 0.0, and `s26/results/e_reproduce.json` re-scores all 126 stored
structures to 0.0 on all four bases (`s26/EXAMINATION.md` H).

### III.1 Sequence and ESM-2 features

Input: the amino-acid sequence (1S9Z: SIRELEARIRELELRI, n = 16, fold 1; 9KAR:
GGWGTVPDWFFNMNW, n = 15, fold 3). ESM-2, a protein language model (`fair-esm`), gives one
vector per residue plus a contact map; the pipeline uses a 32-dimensional PCA of the vectors
(`esm_pca32`) with contact-map and physicochemical pair features, assembled into a
(n_pairs, 183) matrix over all residue pairs at separation >= 2 (105 pairs at n = 16, 91 at
n = 15). The features come from the hot cache `esm_small.npz` (all 126 targets); the 1.5 GB
bank `esm_cache.npz` is never loaded on the production path (`core.pipeline.guard_esm`).

What it discards: everything about the target that is not a function of its sequence. Where it
can be wrong: the language model was trained on evolutionary sequence statistics, not on
peptide structures; the sequence channel is weak at this length (Part I). What it is worth:
ESM features against one-hot residue identity buy 0.288 A on selection
(`docs/FINDINGS.md:1970-2037`, S7 finding 11; `docs/STATE_BRIEF_2026-09-12.md` 5.7).

### III.2 The distogram (the prior)

`core.predict.Distogram.for_target` loads the target's fold model
(`distogram_models/fold<f>_esm_frag.pt`, sha256 in `s26/results/pinned_hashes.json`) and emits,
for every pair, a probability distribution over 17 distance bins with centres 4.0, 4.75, 5.25,
5.75, 6.25, 6.75, 7.25, 7.75, 8.5, 9.5, 10.5, 11.75, 13.25, 15.0, 17.5, 21.0, 25.0 A (`prob`
shape (105, 17) for 1S9Z). From it the score table `_risk` of shape (n_pairs, 760) is
precomputed on a grid from 2.0 to 40.0 A in 0.05 A steps: for pair p and candidate distance t,
`risk_p(t) = w_p * sum_c prob[p, c] * |t - centre_c|` with `w_p = 1 / (sd_p + 0.5)` (the shell
weights file is absent by design, so the shell factor is 1; `s8/inband.py`). The posterior mean
spans 5.05 to 21.07 A on 1S9Z with per-pair sd 0.23 to 2.70 A; on 9KAR 5.62 to 23.44 A with sd
up to 7.58 A.

What it assumes: that a pair's distance is predictable from sequence alone at 17-bin
resolution. Where it can be wrong: the posterior is over-confident by about a factor of two
(`s25/LEDGER.md` L1); its L1 risk's per-pair minimiser is the posterior median, which sits on
one of 17 atoms, so the per-pair target is quantised (`s25/LEDGER.md` L6); the fold models share
most of their training data (Part II). The ORACLE ceiling above it: with the prior replaced by
the native's own distances, everything else fixed, the point cloud reaches 2.2261 A point cloud
against 3.0483 point cloud (`s24/results/priorladder.json`, mean of column MASS1.0 over 126
rows against MASS0.0), a gain of 0.822 A; the slope at the origin is -2.1496 A per unit of
interpolation toward truth (first rung MASS0.1, mean 2.8334), but only along the native's own
direction: moving the score's target 25% toward truth by any achievable route does not move the
endpoint (`s25/LEDGER.md` L12).

### III.3 Retrieval

`core.pipeline.retrieve`: every window of length n from the library of the target's fold (the
out-of-fold peptides, 787 chains in the database, plus the fold's identity-filtered fragments
from 6,003 protein fragments; `ARCHITECTURE.md` 2.1) is scored by BLOSUM62 similarity to the
target sequence; the top K = 500 by a stable argsort are the pool. 1S9Z: 7,193 windows,
similarity -1 to 19, 144 of the 500 from peptides, 420 distinct parent chains; 9KAR: 9,814
windows, similarity -6 to 32. Pool array `W` of shape (500, n, 3), float32 round trip as the
reference did (`core.pipeline._q`).

What it discards: every window outside the top 500, and every conformation not present in the
library. Where it can be wrong: BLOSUM similarity is a weak key; fixing retrieval outright, by
purchasing perfect retrieval with the native, was worth 0.016 A (`docs/FINDINGS.md:2379-2611`,
S8-6), and structure and sequence are decoupled at this length (`docs/FINDINGS.md:477`, S5
finding 9). The stable sort matters because similarity sums are small integers with large tie
sets, and an unstable order moves up to 47 of the 500 members (`README.md`, "Layout"). ORACLE
ceilings, single window: the best window in the pool averages 1.7108 A (`science/pool_best/mean`),
the best window anywhere in the universe 1.3134 A
(`bench_results/recon_library_saturation.json :: universe_best/small`).

### III.4 Score and filter

`core.pipeline.score` looks every candidate's pair distances up in the risk table and sums
(`sc` shape (500,), 0.6926 to 7.2688 on 1S9Z, 2.2420 to 5.1791 on 9KAR; 488 and 483 distinct
values). `filter_pool` keeps the top m = 75 by stable argsort (`sub`, 75 indices; 1S9Z's first
ten are 401, 98, 193, 161, 284, 260, 372, 373, 381, 382) and computes the 75 x 75 pairwise
superposition RMSD matrix among them (the sub-block is bit-identical to slicing the full
500 x 500 matrix,
`tests/test_pipeline.py::test_subblock_pairwise_matrix_is_bit_identical_to_the_full_one`).

What it assumes: that the score's argmin neighbourhood contains the answer. Where it can be
wrong: this is where the programme's accuracy is lost. The argmin alone scores 3.4540 A single
window (`s26/results/e_reproduce.json :: summary/stored_shipped/mean`); the best member of the
top-75 averages 2.3062 A single window (`summary/stored_top_m_best/mean`), so the filter
discards, on average, the 1.7108 A pool member for a set whose best is 0.6 A worse, and on 9KAR
it keeps nothing under 5.34 A while the pool holds a 1.40 A window. A perfect ranker inside the
shipped top-25 would return 2.6087 A single window (`s17/results/inband.json :: cells/25`,
band_best; reproduced in `s26/PREREG_C1.md`); inside the near-native band the score itself is
worse than a coin flip (`docs/FINDINGS.md:2739-2987`, S8-8).

### III.5 The selector (off in production)

With `quantum=True` the filter keeps the top 128 (`top`, 128 indices, 128 x 128 matrix) and
`quantum_stage` (`core/pipeline.py:818-865`) forms `E = _zrank(sc[top])`: the 128 scores are
replaced by their ranks, standardised to zero mean and unit sd, giving a ladder from -1.7186 to
+1.7186 that is the same on every target up to tie-averaging (max deviation 0.394% of range on
both traced targets; worst 1.18% over the S25 sample,
`s25/results/q_gibbs.json :: results/spectrum_target_independence`). `H = diag(E)` is the
Hamiltonian of a 7-qubit register whose basis states index the 128 candidates; a 3-layer,
21-parameter circuit is trained for 50 Adam steps on `F = CVaR_alpha(E; p) - T H(p)` with
(alpha, T) from `VQE_LFO[fold]` (alpha 0.25 on folds 1 and 2, 1.0 on folds 0, 3, 4; T = 0.3
everywhere; `core/pipeline.py:113-118`), and the CVaR tail of the trained distribution is read
out. On 1S9Z (alpha 0.25) the realised tail is 8 states; on 9KAR (alpha 1.0) it is all 128. Part
V is about this stage. Production runs with `quantum=False`, so the 3.2148 A result never
executes it; the four-component arm that does (`bench_results/fourcomponent_tuning126_w8.json`)
reproduces `s8/integrate_vqe.json` (`rmsd_vqe_sel` 3.3135, `medoid128` 3.3443, both
single-window selection on the S8 instrument;
`tests/test_pipeline.py::test_the_four_components_all_execute_and_reproduce_their_published_numbers`).

### III.6 The coordinate average (the point cloud)

`core.pipeline.average`: every retained window is superposed onto the medoid of the set (the
member with the smallest mean RMSD to the others; local index 44, pool index 449 on 1S9Z) and
the arithmetic mean of the 75 traces is taken (`avg_ca`, (n, 3)). Uniform weights were tested
against rank-power and distance weightings with full leakage and are the optimum; removing the
four most deviant members costs +0.142 A (`ARCHITECTURE.md` 2.5).

What it discards: the individual windows; the mean of points scattered about a curve lies
inside it, so the cloud is contracted (virtual bonds 3.75 A mean on 1S9Z, 1.97 A on 9KAR with a
minimum of 1.16 A) and is not a chain. What it assumes: that the 75 errors are independent
enough to cancel. They are not: 68% of the pool's squared error is a common mode shared by every
member (`s23/results/errdecomp.json :: rows[*]/f_common`, mean 0.6758), and averaging removes
only the other 32% (`s23/LEDGER.md` L9). The point-cloud endpoint is 3.0483 A.

### III.7 The projection (the built chain)

`core.pipeline.project` calls `core.project.lam_path`: from four generic starts (extended,
alpha-helix, beta-strand, polyproline II) L-BFGS-B minimises the CA-RMSD between an
ideal-geometry chain built from (phi, psi) and the point cloud, first with no prior, then with a
hinged Ramachandran penalty at lam = 0.3 (`ramah`, calibrated so that a residue as plausible as
95% of real residues feels no force). The emitted chain has every virtual CA-CA bond at
3.803955 A (sd at most 1.8e-15 over the 126 chains,
`s26/results/e_reproduce.json :: summary/virtual_bond`), every bond length and angle at the
constants of `core/geometry.py`, and omega fixed at 180 degrees, so a cis peptide bond (CA-CA
about 2.9 A) cannot be represented. The shipped mode is bit-identical to the reference
`s8.project` on all 126 targets (`verify/project_equiv.json`).

What it assumes: that the nearest ideal chain to a contracted cloud is a good chain. Where it
can be wrong: the projection is degenerate (a CA trace admits torsion solutions at near-equal
distance, one Ramachandran-plausible and one not), which is why the multi-start is science, not
tuning (`tests/test_equivalence.py::test_the_projection_selects_a_different_degenerate_branch`).
The cost of the step is +0.1664 A built chain against point cloud; the chain is better than the
cloud on 16 targets. Recovering the contraction by other means does not help: the Frechet mean
recovers nothing (`s20/LEDGER.md` L11), iterated Procrustes does not help and the contraction
is intrinsic (`s23/LEDGER.md` L4), and the per-target scale that would fix it is real but not
predictable native-free (`s23/LEDGER.md` L2, L3, L6). The built-chain endpoint is 3.2148 A.

### III.8 The relaxation (the relaxed chain)

`core.pipeline.relax` builds every side chain onto the chain with a fixed rotamer, adds
hydrogens, and minimises the AMBER ff14SB force field with the GBn2 implicit solvent under a
10 kcal/mol/A^2 restraint on N, CA and C to convergence (`core.amber.refine_coords`,
k_restraint 10, steps 0). Per target the record keeps the energy before and after, the
restraint RMSD and the residual bond plus angle strain (1S9Z: 5370 to -1290.6 kcal/mol, moved
0.19 A, strain 52.0; 9KAR: 6.0e12 to +1262.4, moved 0.61 A, strain 1166.1;
`s26/results/e_trace_{1S9Z,9KAR}.json`, relax stage). Across the 126: the CA trace moves
0.220 A RMS, the relaxed energy averages -560 kcal/mol, 125 of 126 converge under the
`core.amber.CONVERGE_MAX_KCAL` = 1000 gate (`core/amber.py:1282`; 9KAR does not), and the
virtual bond stretches to 3.867 A on average, past 4.0 A on two targets (2BP4 5.38 A, 9KAR
4.86 A; `s26/results/ph_c3_nativefree.json`, S26 ledger L24). 2BP4 converges by the energy gate
but carries strain 1172.7, above the S8 strain-rejection rule of 1000 kcal/mol
(`s8/integrate.py:292`).

What it assumes: that the force field's local minimum near the built chain is nearer the
native. The record says the step costs +0.0207 A [+0.0143, +0.0276] relaxed chain against built
chain (`docs/STATE_BRIEF_2026-09-12.md` section 4) and that a random displacement of matched
size is at least as accurate (`s16/LEDGER.md` L27); the S26 control experiment C3 measures it
again against a random move of the same size (Part VII). What it buys is validity: exact bond
geometry, strain removed. The relaxed-chain endpoint is 3.2355 A.

### III.9 The ceilings, on one line each

| ceiling (ORACLE unless stated) | mean, 126 targets | basis | artefact |
|---|---|---|---|
| best window in the whole library universe | 1.3134 | single window | `bench_results/recon_library_saturation.json :: universe_best/small` |
| best window in the K = 500 pool | 1.7108 | single window | `bench_results/baseline_tuning126.json :: science/pool_best/mean` |
| best window in the shipped top-75 | 2.3062 | single window | `s26/results/e_reproduce.json :: summary/stored_top_m_best/mean` |
| perfect ranker inside the shipped top-25 | 2.6087 | single window | `s17/results/inband.json :: cells/25` (band_best) |
| the prior replaced by the native's distances, same readout | 2.2261 | point cloud | `s24/results/priorladder.json :: MASS1.0` (mean) |
| shipped argmin (the S8 baseline, not ORACLE) | 3.4540 | single window | `s26/results/e_reproduce.json :: summary/stored_shipped/mean` |
| production point cloud | 3.0483 | point cloud | `science/rmsd_avg/mean` |
| production built chain | 3.2148 | built chain | `science/rmsd_arm/mean` |
| production relaxed chain | 3.2355 | relaxed chain | `science/rmsd_full/mean` |
| constant alpha-helix (zero information, not ORACLE) | 4.0648 | built chain | `s12/results/s14_ladder.json :: rows/L0_constant_helix/mean` |

The single-window ceilings and the averaged results are on different bases and are shown
together only to place the stages; the one contrast the programme makes across that line is
the gap between what the pool holds (1.71 A, single window) and what the readout returns
(3.21 A, built chain), which is the recognition problem every sprint from S7 to S25 worked on.

## PART IV. THE TWO PHYSICS HAMILTONIANS

Source: `s26/PH_PART_IV_NOTES.md` (lane PH, S26), which carries the artefact path of every
number below; the underlying files are `s25/results/phys_landscape.json`,
`s25/results/phys_suite.json`, `s13/results/walsh_xval.json`, `s13/results/walsh_amber.json`,
`s20/results/c_land_report.txt`, `s20/results/c_land_null.json`,
`s26/results/ph_reject_census.json` and `s26/results/ph_c3_nativefree.json`. Basis is named at
every RMSD.

### IV.1 What the two scores are

The programme has two ways of asking "is this a plausible peptide?" that do not read the native.

**Legacy** is an eleven-term score (`core/energy.py`, `DEFAULT_WEIGHTS`, never fitted): steric
4.0, contact 1.0, local hydrogen bond 1.0, long-range hydrogen bond 3.0, helix cooperativity
2.0, sheet cooperativity 2.0, solvation 0.5, electrostatic 1.0, aromatic 0.8, torsion 0.15,
compactness 0.4. Each term reads the backbone and CB atoms only (the CB is the first side-chain
atom) and adds a number; the weights are hand-set constants. It costs 0.3 ms per structure.

**AMBER** is a molecular force field, ff14SB with the GBn2 implicit-solvent model, evaluated by
OpenMM (`core/amber.py`). Every heavy atom of every side chain is built onto the backbone by a
deterministic builder (a fixed rotamer, no scan), hydrogens are added, and the energy is the
sum of bond, angle, torsion, non-bonded (Lennard-Jones plus Coulomb) and solvation terms. A
single evaluation costs about 9 ms; a restrained minimisation to convergence 6 to 12 s
(`s16/repair_FINDINGS.md` section 2.4). The record's earlier "6 ms" was a memoised re-evaluation
and is corrected to 28 ms for the single point it timed (`docs/CONDENSED_REPORT.md:185`).

### IV.2 What "H = E after 50 relaxation steps" means

When AMBER is used as a Hamiltonian for a quantum selector, the structure it is handed is an
ideal-geometry chain that has never been relaxed, and its energy there is not a finite number on
42% of a lattice register (`s20/LEDGER.md` L6: finite at Relax_1 on 112 of 192 states, at
Relax_50 on 186 of 192). The deployable object is therefore "run 50 steps of restrained
minimisation, then read the energy". That operator, not the raw force field, is the
Hamiltonian; the relaxation is what makes the objective defined, and the Spearman correlation
between the energy after 50 steps and after 1 step is 0.36 to 0.88 across three targets, so the
ordering is still moving when the cap stops it.

### IV.3 Legacy is a compactness model; AMBER has the opposite sign

Measured on the shipped K = 500 pools of all 126 targets
(`s25/results/phys_landscape.json :: summary`): the 75 candidates Legacy likes best are 0.758 A
more compact in radius of gyration than the pool (SE 0.023) and Legacy's rank correlation with
Rg is +0.60. The 75 AMBER likes best are 1.103 A more expanded (SE 0.042), rho(AMBER, Rg)
= -0.27. They disagree with each other about ordering: rho(Legacy, AMBER) = -0.090 (SE 0.018),
reproduced on three instruments (`s20/results/c_q1.json` -0.0886; S24 -0.0829). Legacy's
whole-pool correlation with the true RMSD is +0.307, but a score that reads only Rg gets +0.279,
so almost all of Legacy's apparent skill is "prefer compact". AMBER's is -0.027 whole-pool and
+0.047 with Rg removed. (`s20/LEDGER.md` L8 is the first statement of this.)

### IV.4 Why both rank worse than a random subset

Used as the selector over the pool, with the identical coordinate-average readout and a CVaR-VQE
for all seven configurations (`s25/results/phys_suite.json :: configs_rank, random_null_rank`;
`s25/LEDGER.md` L16), the endpoints are, ALL POINT CLOUD: Distogram 3.058, AMBER+Distogram
3.132, Legacy+Distogram 3.215, all three 3.253, random 75-subset 3.425, Legacy+AMBER 3.674,
Legacy 3.755, AMBER 3.881 A. Legacy alone is +0.330 A point cloud worse than picking 75 at
random (point cloud), 1.99x its MDE; AMBER +0.455 A, 2.42x MDE; 5 of 5 folds each
(`s25/agentPHYS_FINDINGS.md` section 1.4). The built-chain means of the same seven rows are
3.2187, 3.3100, 3.3732, 3.4221, 3.8248, 3.8844 and 4.1015 A (`results/summary/leaderboard.json
:: rows[*]/mean`; the random-75 built-chain row is not in that table, so the "worse than
random" contrast is stated on the point cloud only). Permuting a physics channel while keeping
its distribution improves the endpoint (Legacy +0.327 point cloud worse than its own noise,
AMBER +0.471; `s25/agentPHYS_FINDINGS.md` section 1.5). The reason is IV.3: the pool's dominant
axis is compactness, the distogram already selects on it, and each energy pushes along that
axis in a direction unrelated to which candidate is right. S16 had found the same for Legacy as
a ranker inside the architecture (+0.068 [+0.012, +0.125] worse than a matched random drop at
m = 5, `s16/LEDGER.md` L17), S17 had found the decisive AMBER-against-Legacy comparison at full
scale and concluded "the objective is not weak, it is wrong" (`s17/LEDGER.md` L29, L30), and S18
found every physics filter losing to a random gate (`s18/LEDGER.md` L17).

### IV.5 The steric singularity, and where it lives

An unrelaxed ideal-geometry rebuild of a retrieved window puts atoms on top of each other, and
the Lennard-Jones r^-12 wall turns one overlap into an energy of 1e4 to 1e29 kcal/mol (a relaxed
peptide sits at -1170 to -500). Measured: 58.6% of every pool is above 1e4 kcal/mol; the
energies span 15.3 decades; 97.0% of a pool lands inside |z| < 0.1 of a moment z-score; the ten
worst candidates carry 99.66% of the variance
(`s25/results/phys_landscape.json :: AMB_frac_absz_lt_0p1, AMB_decades, AMB_top10_var_share`).

In the Pauli basis the same fact reads: over the 25 fully enumerated AMBER tables of S13, the ten
most extreme configurations out of 4,096 carry a median 99.56% of raw AMBER's Walsh variance
(range 60.7% to 99.98%; the single worst configuration alone a median 46.5%), against a median
26.6% for Legacy over 41 tables (`s13/results/walsh_xval.json :: concentration`,
`var_share_top10`, `var_share_top1` by `model`). A function that is a constant plus one spike
has a Pauli-weight spectrum of exactly Binomial(m, 1/2), so the spike makes raw AMBER look
maximally non-local for arithmetic reasons that say nothing about the physics
(`s13/walsh_FINDINGS.md`; Part V.7). The spike is one term: in the same enumerations the
non-bonded term owns the whole variance, covariance share 1.000 on every cell, while bond and
angle terms have variance 0.0 and torsion and solvation have variances of order 3 and 900
(kcal/mol)^2 against the non-bonded term's 7e30 to 2e31
(`s13/results/walsh_amber.json :: exact[*].terms`; `s13/walsh_FINDINGS.md` section 3.1). Remove
that one term and AMBER's spectrum falls below Legacy's on 6 of 6 targets (section 3.2 there).

In the torsion-space Hessian the same fact reads: 7.45% of AMBER's modes carry all of its
curvature (participation ratio 0.0745 against Legacy's 0.4221, 30W/0L; anisotropy 18.8 against
5.8, 0W/30L; condition number three orders larger; `s20/results/c_land_report.txt` section 1).
When AMBER is minimised from a pool member it moves the chain 0.577 rad per coordinate, five
times Legacy's 0.104, because its first steps are clash relief; 72% of the RMSD damage that
minimisation does is the size of that move, not its direction
(`s20/results/c_land_null.json`, toward-member null +0.444 of +0.620; `s20/LEDGER.md` L12).

New in S26 (`s26/results/ph_reject_census.json :: singularity`; S26 ledger L23): on the 9,450
shipped top-75 rebuilds, 40.7 of every 75 have two heavy atoms closer than 2.0 A when every side
chain is built, but only 2.6 of 75 do on the backbone plus CB (S19 measured 2.66). Of the 5,057
rebuilds above 1e4 kcal/mol, 96.8% have their closest contact on a side-chain atom (bb-sc 2,627,
sc-sc 2,267, bb-bb 163). Within a top-75 the AMBER single point tracks the minimum heavy-atom
distance at Spearman -0.74. So the singularity is mostly the builder's: the fixed rotamer that
places each side chain onto a backbone it did not see. The backbones themselves are almost
always physically possible.

### IV.6 The non-monotone standardisation trap

Turning an energy into a score by subtracting its mean and dividing by its standard deviation
looks harmless, but with one 1e28 outlier setting the standard deviation, every candidate below
about 1e12 maps to the same double-precision number. On 40 of 126 targets the moment z-score
changes AMBER's own ordering, with exact tie blocks up to 462 of 500; `np.argsort` then returns
those tied candidates in array order, which is the BLOSUM retrieval order, and that order is not
neutral (rho with the true RMSD +0.054, seven standard errors from zero). Every Angstrom of the
moment-z "advantage" on the AMBER configuration lived on exactly those 40 targets
(`s25/agentPHYS_FINDINGS.md` section 3; `s25/results/phys_suite.json :: normalisation_fork`).
Rank standardisation is monotone and has no such failure; it is what production uses
(`core.pipeline._zrank`; `ARCHITECTURE.md` section 4).

### IV.7 The relaxation on the production chain

The last production stage relaxes the built chain with a restraint of 10 kcal/mol/A^2 on N, CA
and C. Measured on the 126 emissions (`s26/results/ph_c3_nativefree.json`; S26 ledger L24):
the built chain's own AMBER energy is above 1e4 kcal/mol on 58.7% of targets before relaxation
(median 8.6e4), the relaxation brings it to -560 on average and converges on 125 of 126 (9KAR
ends at +1262), moves the CA trace 0.220 A RMS, and stretches the virtual CA-CA bond from 3.804
to 3.867 A on average, breaking it beyond 4.0 A on two targets (2BP4 5.38 A, 9KAR 4.86 A). What
that step does to accuracy against a random move of the same size is C3 (`s26/C3_RESULT.md`,
after the gate; Part VII); on the record so far the step costs +0.0207 A [+0.0143, +0.0276]
relaxed chain against built chain (`docs/STATE_BRIEF_2026-09-12.md` section 4) and S16 found a
random displacement of matched size at least as accurate (`s16/LEDGER.md` L27). S17's "steric
validity is free" was quoted against the wrong baseline and corrected by the workstream that
made it (`s17/LEDGER.md` L27), and S17 explained S16's cis-peptide defect in a way that
exonerates the force field: the ideal-geometry builder cannot represent the cis bond, and the
minimiser is asked to repair what the builder could not draw (`s17/LEDGER.md` L28).

### IV.8 What the physics is for

The record's settled role for both energies is a validity stage, not an accuracy stage:
`docs/FINDINGS.md:3663` (S8-12, "AMBER is a validity stage, not an accuracy stage");
`s17/LEDGER.md` L24 (the validity frontier splits in two, and the incumbent's restraint set
wins); `s18/LEDGER.md` L15, L17 (`leg_contact` fails as a term in the objective and is harmful
on the deployed one); `s21/LEDGER.md` L24, L33 (the Legacy-to-AMBER continuation is degenerate
in raw units, preconditioning does not work, and the physics half is worse than chance at
n = 126); `s24/LEDGER.md` L10, L10-A, L16 (two different physics functionals do not escape the
shared referent; AMBER selects along a non-parallel direction and cannot pay for it; the
functional lever is closed in all three forms); `s25/LEDGER.md` L16, L18 (the seven-configuration
suite; closed in all five forms). The sentence lane PH wrote for the presenter is the
programme's position: the two physics scores are real, they disagree with each other about what
a good peptide looks like, and on this pool both are worse than choosing at random because
neither is looking at the thing that decides accuracy; the all-atom force field's number on an
unrelaxed candidate is almost entirely the distance between its two closest atoms, and that
distance is set by how the side chains were placed, not by the backbone. It belongs in the
pipeline as a validity check, which is what it is, and not as a judge of which candidate is
right, which it is not (`s26/PH_PART_IV_NOTES.md` section 9).

## PART V. THE QUANTUM COMPONENT, IN FULL

Sources: `s25/QUANTUM.md` (the S25 account, sections 0 to 8, with its artefacts
`s25/results/q_verify.json`, `q_alpha.json`, `q_gibbs.json`, `q_plateau.json`),
`s13/SPRINT13_DOSSIER.md` sections 6 to 11 (the trainability record), `s26/LEDGER.md` L27 with
`s26/results/q_dla.json` (the dynamical Lie algebra), and the sprint ledgers named inline. The
selector's place in the pipeline is Part III.5.

### V.1 Two quantum paths, and which one is deployed

There are two distinct quantum objects in the codebase, and conflating them is the first
mistake available (`s25/QUANTUM.md` section 1):

| | the selector | the generation lane |
|---|---|---|
| ansatz | `StatevectorCircuit` (`core/quantum.py:885-982`) | `MPSAnsatz` (`core/quantum.py:459-776`) |
| register | candidate identity, n = 7, dim 128 | torsion or basin latent, n up to 96 |
| deployed by | `core.pipeline.quantum_stage` | the S19 to S21 lanes |
| simulation | exact dense statevector | exact matrix product state, bond dimension chi = 2^layers |

`core/pipeline.py` does not use `MPSAnsatz`. Statements about chi apply to the generation lane;
statements about the deployed selector apply to `StatevectorCircuit`. Production runs with
`quantum=False`; the selector is executed by the four-component arm and by the S25 analysis
scripts.

### V.2 The ansatz, exactly

The selector prepares `|psi(theta)> = prod_{l=1..L} [U_ent (x)_{q} RY(theta_{l,q})] |0...0>`
with `RY(t) = exp(-i t Y/2)` and an entangler that is a CNOT chain plus a ring closure,
`U_ent = CNOT(n-1,0) CNOT(n-2,n-1) ... CNOT(0,1)`. Deployed: n = 7 qubits, 128 basis states,
L = 3 layers, P = 21 parameters, 21 RY and 21 CNOT gates (6 chain + 1 ring per layer), two-qubit
depth 21, total depth about 24, 50 Adam steps, one restart, seed 0
(`core/pipeline.py:186-189`). RY and CNOT are real matrices and the initial state is real, so
the amplitudes are real exactly and the reachable manifold lies in SO(2^n), not SU(2^n);
dim so(128) = 8128 against 21 parameters. The simulation carries all 2^n amplitudes; the CNOT
chain is composed into one basis permutation, and `probs_batch` simulates all 2P shifted
parameter vectors in one pass. Verified in S25 against a dense Kronecker-product simulator:
max |p - p_dense| = 5.6e-17 (`s25/results/q_verify.json`).

The generation-lane `MPSAnsatz` is layers of (RY on every wire, CNOT chain), no ring, with an
optional trailing RY layer. A CNOT is an exact bond-dimension-2 matrix product operator, so
after L entangling layers every interior bond is exactly 2^L and nothing is discarded: there is
no SVD and no truncation threshold anywhere in the class (verified: no truncation primitive;
chi by layers {1: 2, 2: 4, 3: 8, 4: 16}; chi independent of n at n = 64; max |p_MPS - p_dense|
3.331e-16 over 24 configurations; max |<psi|psi> - 1| 6.661e-16). Cost is O(n chi^3), linear
in n, which is why 64- and 96-qubit registers were routine in S19 to S21. The bond dimension is
saturated, not merely bounded: the Schmidt spectrum across the middle cut at n = 6, L = 2 is
[0.879145, 0.461538, 0.105625, 0.054141, 0, 0, 0, 0], rank exactly 4. The product-state control
(`entangler="none"`) at identical angles moves a probability by up to 0.334, so the entangler
does change the state (`s25/QUANTUM.md` 2.2, 7.4).

### V.3 The objective, exactly

`E = zrank(score[top[:128]])`, `H = diag(E)` (`core/pipeline.py:852`). `zrank` is the
standardised rank: strictly monotone, so it changes no ordering, no argmin, no level set and
(by V.6) no tail membership; it changes only the gaps, which is what the entropy term trades
against. It is not cosmetic: a moment z-score of raw energies is outlier-dominated (Part IV.6).

**The Hamiltonian barely changes between targets** (`s25/LEDGER.md` L17). `top = argsort(sc)`
(`core/pipeline.py:768`), so `sc[top[:128]]` is already ascending and `_zrank` of an ascending
vector returns the standardised ranks 1 to 128. E is therefore very nearly the same vector on
every one of the 126 targets; it differs only through tie-averaging. On 8 real targets the
worst |E - standardised ranks| is 4.06e-2 against an E range of 3.4371, 1.18% of range, and E
is exactly the standardised ranks on 0 of 8 (`s25/results/q_gibbs.json ::
results/spectrum_target_independence`; the two S26 traces show 0.394%, Part III.5). All of the
per-target information enters through which candidate occupies which rank, and essentially
none through the spectrum of H. Three consequences: the deployed selector solves very nearly
the same variational problem on all 126 targets, so there are effectively two trained states in
the whole deployment, one per (alpha, T) cell, not 126; the quantum stage is insensitive to the
target; and a more expressive state has nothing target-specific to be expressive about, which
is the mechanism behind five sprints of "deeper or wider orders nothing" (V.11).

**CVaR.** For alpha in (0, 1] and q the alpha-quantile of E under p,
`CVaR_alpha = (1/alpha) [sum_{E(x) < q} p(x) E(x) + (alpha - P(E < q)) q]`: the mean of the
lowest-alpha mass of the energy distribution, the boundary state contributing fractionally. Two
implementations, `cvar_exact` (`core/quantum.py:985-1005`, the deployed driver) and
`cvar_from_probs` (`core/quantum.py:336-359`, which also returns the per-state mass the tail
used), are mutually checked; `tests/test_cvar.py` (90 tests) asserts the definition on both
sort-cutoff branches.

**The deployed objective is a free energy**, `F(theta) = CVaR_alpha(E; p_theta) - T H(p_theta)`
with `H(p) = -sum p log p` (`core/quantum.py:1072-1100`). Minimising CVaR alone is degenerate
for a selection task: the minimiser concentrates p on the lowest-energy states and a consensus
readout collapses to the argmin; measured at alpha = 1, T = 0.1 the state carries 0.0761 bits of
a possible 7 and the arm is identical, target by target, to the plain argmin selector. The
(alpha, T) pair comes from a leave-fold-out table (`core/pipeline.py:113-118`):
`VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}`.
T = 0.3 on all five folds; alpha = 1.0 on three of them, where CVaR is the full mean and there
is no tail constraint at all, 78 of 126 targets (share 0.6190,
`s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint`;
`s25/LEDGER.md` L3).

**What is differentiable.** CVaR as a function of p is concave and piecewise linear; by the
envelope theorem `dCVaR/dp(x) = (E(x) - q)/alpha` on the strict tail and 0 elsewhere, with the
boundary state carrying its partial mass; it is non-differentiable exactly where the quantile
crosses a state or where energies tie, a measure-zero set. The entropy is smooth on the
interior (p clipped at 1e-15 inside the log). `p_theta` obeys an exact shift rule (V.4). The
readout is not differentiable at all and is never differentiated: the consensus medoid is an
argmin over candidates.

### V.4 The gradient

Each parameter enters as RY, whose generator Y/2 has eigenvalues +-1/2, so every basis
probability obeys the two-term parameter-shift rule exactly:
`dp(x)/dtheta_k = [p(x; theta_k + pi/2) - p(x; theta_k - pi/2)]/2`. Chained with `dCVaR/dp`
and `dH/dp` this gives `dF/dtheta` exactly, at 2P = 42 circuit evaluations per gradient, batched
into one pass (`grad_cvar_paramshift`, `core/quantum.py:1008-1023`). Verified in S25 against
central finite differences on the exact objective (independent machinery): cosine
1.000000000 (minimum over 12 cells, alpha in {0.1, 0.25, 1}), relative error 4.597e-10 (maximum
over 12); for the free-energy gradient including the entropy term, cosine 1.000000000 and
relative error 4.663e-10 over 4 (alpha, T) cells including the deployed one
(`s25/results/q_verify.json`). `run_cvar_vqe` (`core/quantum.py:1103-1130`) is Adam
(beta1 0.9, beta2 0.999, lr 0.15) on that gradient; it samples nothing, the only random draw
being the initial angles theta ~ N(0, 0.6^2).

### V.5 The CVaR gradient defect

The sampled (device-realisable) estimator uses the score-function form,
`grad CVaR = E_p[f(x) grad log p_theta(x)]` with `f(x) = -(q - E(x))_+/alpha`. A baseline b may
be subtracted from f if and only if it is constant in x, because the correction
`b E_p[grad log p]` vanishes only then. The historically shipped `qansatz.cvar_gradient`
subtracted the tail mean from the tail entries and left the rest at zero, `b(x) = m 1[x in
tail]`, a function of x. That is a bias, not extra variance; it does not shrink with shots.
Both `cvar_gradient` (`core/quantum.py:816-863`) and `grad_cvar_score`
(`core/quantum.py:1040-1069`) now default to `baseline="const"`; `baseline="tail"` reproduces
the defect verbatim so that the regression test in `tests/test_quantum.py` measures its absence
rather than asserting it (39 tests passed, `s26/TEST_RUN.md`). The magnitude is
instrument-dependent and is always cited with its instrument, all three with zero sampling
noise, which is what proves the disagreement is bias (`s25/QUANTUM.md` 4.2):

| instrument | cosine with the exact gradient | norm ratio |
|---|---|---|
| S9, 10 qubits, 1024 amplitudes enumerated (`core/quantum.py:66`) | +0.655634 | 0.758 |
| S25 re-verification, n = 7, L = 3, alpha = 0.15, deployed energy shape | +0.566586 | 0.519 |
| any instrument, constant baseline | +1.000000 | 1.000 |
| sampled estimator, constant baseline | +0.994 | about 1 |

A third figure, +0.524 over 36 checks on the S8 instrument, lives only in project memory and is
not in Appendix B (`s26/EXAMINATION.md` C24). There is no universal constant; what reproduces is
the sign, the order of magnitude and the mechanism. The defect was first found and priced in
S5 (`docs/FINDINGS.md:368`, "A defect in the shipped CVaR gradient") and re-audited against two
references in S8-9 (`docs/FINDINGS.md:3147`) and S9-5 (`docs/FINDINGS.md:4197`).

### V.6 The tail-subset-of-prefix theorem

Statement (`s25/QUANTUM.md` section 5): the realised CVaR tail's support is always a subset of
an initial prefix of the energy order, and equals that prefix exactly when every state in the
prefix carries positive probability. The trained state can delete a member; it can never add
one outside the classical top-m.

Proof, from the source. `cvar_from_probs` computes `order = argsort(e)`, `cum = cumsum(p[order])`
and `take = clip(alpha - (cum - p[order]), 0, p[order])` (`core/quantum.py:352-354`).
`take[j] > 0` needs both conjuncts of the clip: (i) the exclusive prefix sum in the energy order
is below alpha, and that sum is non-decreasing because p >= 0, so it crosses alpha exactly once
and (i) alone defines an initial prefix; (ii) `p[order][j] > 0`, which punches out any
zero-probability state inside the prefix. So the support is a prefix with holes, every hole a
zero-probability state. On the other path, `tail_indices` (`core/quantum.py:229-272`) takes no
probability vector and no parameters: `k = ceil(alpha n)` depends on (n, alpha) only.

Consequences. The trained state moves exactly two things: where the prefix cuts (a rung m on
the classical top-m ladder, plus which prefix members it deletes where it has exact zeros) and
the weights inside the prefix. Membership is bounded above by `argsort(E)[:m]` as an identity,
independent of the landscape's shape, degeneracy, multimodality and scale, and therefore of the
pool that produced it. Subset-hood is the theorem; equality is the empirical regime (a trained
RY/CNOT state has generic angles and full support), and an earlier draft that asserted
prefix-hood without conjunct (ii) was wrong.

Verification, S25, families rebuilt with exact zeros: 6 register sizes x 9 energy structures x
8 probability structures x 6 alphas = 2,592 cells; 1,620 contain an exact zero probability, so
the assertion can fire; subset-hood violations 0; holes that were not exactly zero-probability
0; prefix-hood violations (the wrong, looser claim) 1,424 = 54.9%; full-support cells 972, with
exact value equality on 972 of 972 (`s25/results/q_verify.json`). S24's independent run: 3,888
cells, 0 violations, 29.9% hole rate; the coordinator's: 17,574 trials, 0 violations, 58.8%. The
hole rates differ because the family mixes differ; the three assertions agree exactly.

The honest statement about S22's Gate 1: a set-equality gate reported passing on 2,016 of 2,016
cells was read at the time as evidence about the candidate pool. It was the code being read
back. The same applies to the harness test that shows the quantum arm matching its classical
control to 1e-9. A 100% pass rate on an algebraic identity carries no information about the
problem.

### V.7 The contribution, with the measurements

Basis notice. Two instruments appear in this section and never share a column: the S8
instrument (`s8/integrate_vqe.json`), consensus-medoid selection of one member from a
128-candidate filtered set, SINGLE WINDOW, whose no-circuit rungs are 3.4540 (argmin), 3.3414,
3.2835 and 3.3135 (the recorded VQE arm); and the 126-target production instrument, a
score-filtered uniform top-75 coordinate average, 3.0483 POINT CLOUD and 3.2148 BUILT CHAIN.
Everything measured below is on the S8 instrument at n = 126 over the five pinned folds, with
MDE = 2.8016 x SE per comparison (`s25/q_alpha.py` -> `s25/results/q_alpha.json`).

**The alpha effect does not survive the temperature** (`s25/QUANTUM.md` 6.1). Differences of
mean single-window RMSD, alpha = 0.25 minus alpha = 1.0: at T = 0.1, -0.1067 (SE 0.0644,
0.59x MDE, underpowered); at the deployed T = 0.3, -0.0011 (SE 0.0477, 0.01x MDE, 44W/47L,
null); at T = 1.0, +0.0039 (SE 0.0240, 0.06x MDE, null). Alpha = 0.10 minus 1.0: -0.1126 (SE
0.0792, 0.51x) at T = 0.1, +0.0279 (SE 0.0518, 0.19x) at T = 0.3, +0.0126 (SE 0.0268, 0.17x) at
T = 1.0. The effect appears only at T = 0.1 and changes sign at the other two temperatures. The
"+0.113 A CVaR contribution" of the earlier record was measured at a temperature that does not
ship and was withdrawn (`s25/LEDGER.md` L3, L5).

**One curve: the readout's entropy is the variable** (6.2). All 18 arms in the artefact pass
through the same operator, `consensus_medoid(D, o, w)`, and differ only in the weight vector w,
whose entropy is computable in closed form for the non-circuit arms (argmin 0 bits, top-fraction
log2 k, Boltzmann H(exp(-E/T)), medoid over 128 = 7 bits). Over the nine VQE cells,
corr(H_readout, mean RMSD) = -0.7423, corr(alpha, mean RMSD) = +0.2700, corr(T, mean RMSD)
= -0.0234; alpha's marginal share of the variance left after H and H^2 is 0.032; a quadratic in
H over all 18 arms has R^2 = 0.7045. Fitted on the nine no-circuit arms only, the curve scores
the nine circuit arms with mean residual +0.0090 A (sd 0.0342) against the fit's own residual
sd 0.0268 A: the circuit sits on a curve fitted without it.

**The circuit against its own analytic optimum** (6.3). At alpha = 1 the objective is
`mean_p(E) - T H(p)`, whose minimiser over the simplex is exactly the Gibbs distribution
`p*(x) = exp(-E(x)/T)/Z`, so `boltz_T` is the optimum the 21-parameter state approximates. At
the endpoint, circuit minus exact Boltzmann, single window: +0.0499 (SE 0.0427, 0.42x MDE) at
T = 0.1; -0.0302 (SE 0.0450, 0.24x MDE, 43W/44L) at T = 0.3; +0.0445 (SE 0.0511, 0.31x) at
T = 1.0. In distribution (`s25/q_gibbs.py` -> `s25/results/q_gibbs.json`), using the identity
`F(p) - F(p*) = T KL(p || p*)`, asserted at runtime to < 1e-9 on every row, at T = 0.3:

| state | F | KL nats | KL bits | TV | H bits |
|---|---|---|---|---|---|
| random theta (untrained) | -1.275443 | 3.927427 | 5.66608 | 0.78052 | 4.4562 |
| uniform over 128 | -1.455609 | 3.326874 | 4.79966 | 0.70154 | 7.0000 |
| point mass at the argmin | -1.718572 | 2.450332 | 3.53508 | 0.91374 | 0.0000 |
| 5 Adam steps | -1.644194 | 2.698258 | 3.89276 | 0.70994 | 5.1019 |
| 15 Adam steps | -1.908150 | 1.818405 | 2.62340 | 0.63160 | 5.4014 |
| 50 Adam steps (deployed) | -2.183096 | 0.901916 | 1.30119 | 0.45308 | 5.6706 |
| the Gibbs optimum | -2.453671 | 0 | 0 | 0 | 4.9135 |

KL(trained || optimum) is 1.449 nats at T = 0.1, 0.902 at T = 0.3, 0.373 at T = 1.0; total
variation 0.761, 0.453, 0.351. The mandatory control, best of 200 draws from the untrained
circuit: F_trained -2.183096 against F_init_best -1.483528 and F_gibbs -2.453671 at T = 0.3,
gap closed 0.783, beats best-of-200 at all three temperatures (0.891 at T = 0.1, 0.783 at
T = 1.0). Three readings, in order: the optimiser trains (it beats best-of-200 and closes 78 to
89% of the free-energy gap); it does not reach the optimum (0.902 nats away at the deployed
temperature, disagreeing with it on 45% of its mass, and broader than optimal at 5.67 bits
against 4.91); and the endpoint cannot tell (0.24x MDE). The readout is insensitive to a
distributional difference of nearly half the mass (`s25/LEDGER.md` L15), which is the mechanism
behind the entropy curve and the sharpest fact in the quantum record. Two details: at T = 0.1
the trained state's free energy (-1.7176) is marginally worse than the point mass at the argmin
(-1.7186) with entropy 0.075 bits, the collapse; at T = 1.0 the trained state (KL 0.373) is only
modestly better than the uniform distribution (KL 0.458).

**Against the no-circuit arms** (6.4), single window: VQE_LFO minus argmin -0.1405 (SE 0.0732,
0.68x MDE, 66W/48L, 5 of 5 folds in sign; underpowered, not a result); minus Boltzmann at
T = 0.3, -0.0002 (0.00x MDE); minus uniform top-64, +0.0262 (SE 0.0279, 0.34x); minus uniform
top-128, -0.0308 (SE 0.0590, 0.19x). The medians are near zero (VQE_LFO minus argmin median
-0.0081 against mean -0.1405; 12 of 126 targets tie exactly; 5 targets carry 45.3% of the
effect), and the drop-top-5 statistic sits at the 64.7th percentile of a uniform-effect null,
so concentration is suggested and not established.

**The alpha = 1 folds are a claims problem, not a performance problem** (6.5). Forcing alpha < 1
on every fold, an ORACLE counterfactual not applied, is worth -0.0311 A at 0.27x MDE
(alpha = 0.25) or -0.0021 A (alpha = 0.10). Nothing about the deployed configuration changes on
the strength of this; what changes is what is claimed about it.

**Deleting the entangler.** The standing project result from removing the CNOTs and changing
nothing else is -0.013 A [-0.095, +0.077], an interval containing zero and effects in both
directions; and chi orders nothing (`s21/LEDGER.md` L34). Whether the entanglement changes the
answer is therefore not measured and specifically not refuted (`s25/QUANTUM.md` 7.4).

### V.8 The width and depth sweeps, and what they do and do not say

`s25/q_plateau.py` -> `s25/results/q_plateau.json`, exact parameter-shift gradients (no shot
noise), theta ~ N(0, 0.6^2), the deployed initialisation law, L = 3 throughout the width sweep.
`Var_theta[dF/dtheta_0]`:

| n | dim | P | alpha = 1, T = 0 (linear cost) | alpha = 0.25, T = 0 | alpha = 0.10, T = 0 | alpha = 1, T = 0.3 | alpha = 0.25, T = 0.3 |
|---|---|---|---|---|---|---|---|
| 4 | 16 | 12 | 5.358e-02 | 1.453e-02 | 4.699e-03 | 5.491e-02 | 2.088e-02 |
| 6 | 64 | 18 | 1.699e-02 | 1.140e-02 | 3.214e-03 | 1.188e-02 | 1.559e-02 |
| 7 | 128 | 21 | 2.062e-02 | 1.064e-02 | 1.922e-03 | 1.526e-02 | 1.121e-02 |
| 8 | 256 | 24 | 4.358e-03 | 6.340e-03 | 2.284e-03 | 8.389e-03 | 7.837e-03 |
| 10 | 1024 | 30 | 2.713e-03 | 5.517e-03 | 3.175e-03 | 5.731e-03 | 6.121e-03 |
| 12 | 4096 | 36 | 1.218e-03 | 3.222e-03 | 2.524e-03 | 7.741e-03 | 6.376e-03 |
| 13 | 8192 | 39 | 1.072e-03 | 3.600e-03 | 2.869e-03 | 4.802e-03 | 4.113e-03 |
| fitted log2 Var per qubit | | | -0.6492 | -0.2522 | -0.0472 | -0.3105 | -0.2429 |

No exponential plateau is present at any width measured, in any column, at this depth. The
steepest decay is the linear cost at -0.649 log2 per qubit, about 1.6x per qubit added, against
the factor of 2 a 2-design would give. Draws per row are 250 (n <= 8), 200 (n = 10), 120
(n = 12), 80 (n = 13); the relative SE of a variance estimate is sqrt(2/(m - 1)), about 9% at
m = 250 and 16% at m = 80. Depth sweep at n = 7, alpha = 0.25, T = 0.3, 250 draws: L = 1
2.381e-02, L = 2 1.658e-02, L = 3 9.197e-03 (deployed), L = 4 7.999e-03, L = 6 7.644e-03, L = 8
7.440e-03, L = 12 8.213e-03: depth costs a factor of about 3 and then saturates.

**The CVaR non-linearity's effect on the picture** (7.2). At alpha = 1 the CVaR reduces
identically to the linear cost `<psi| diag(E) |psi>`, so the linear control is the same
expression with the non-linearity switched off, every confound held fixed by construction. The
ratio Var[grad CVaR_alpha] / Var[grad MEAN] at matched n: alpha = 0.25 gives 0.2712, 0.6707,
0.5159 (n = 7, deployed), 1.4549, 2.0336, 2.6447, 3.3585 for n = 4, 6, 7, 8, 10, 12, 13;
alpha = 0.10 gives 0.0877, 0.1892, 0.0932, 0.5241, 1.1704, 2.0715, 2.6765. In magnitude at the
deployed width the non-linearity shrinks the gradient (a tail objective's weight is supported on
an alpha fraction of the states); in scaling it flattens the decay, from -0.649 log2 per qubit
to -0.252 (alpha = 0.25) and -0.047 (alpha = 0.10), so the ratio crosses 1 near n = 8. The
entropy term also flattens the decay (-0.311 and -0.243 in the T = 0.3 columns). A mechanism
was offered as a reading, not tested.

**Scope conditions**, so this is not over-read. The deployed register is n = 7; every larger n
is an extrapolation instrument. P = 3n against dim so(2^n) = 2^(n-1)(2^n - 1), 21 against 8128
at n = 7, so the circuit is nowhere near a 2-design at any width measured, and the observed
decay rate is a property of this shallow, structured ansatz; the classic exponential-in-n
plateau is a statement about 2-design circuits and this is not one. Not measured in S25: the
dynamical Lie algebra (now measured, V.9); gradient variance along the optimisation trajectory
rather than at random initialisations; anything on hardware; whether the same ordering holds
for a different spectral shape. Everything reported is computed exactly and is
shot-count-independent; on hardware the deployed gradient needs the full 128-outcome
distribution resolved, which is distribution reconstruction rather than a few Pauli
expectations, and the two device routes in the codebase (sampled CVaR with SPSA,
`run_global_cvar_vqe` at `core/quantum.py:1448`; the score-function gradient with a constant
baseline) produced none of the numbers above (7.3).

### V.9 The dynamical Lie algebra (S26, lane Q, ledger L27)

`s26/q_dla.py` -> `s26/results/q_dla.json` (complete; 85.3 s, peak RSS 0.479 GB), pre-registered
in `s26/PREREG_A2.md` before the run; a property measurement with no native and no score. The
closure is computed exactly on Pauli strings as a set and cross-checked by dense SVD at
n = 4, 5 (12 of 12 cells agree; the two conjugation conventions agree on 12 of 12). Every
closure lies in the odd-Y (real) set.

| ansatz or pool | dim(DLA) | artefact leaf |
|---|---|---|
| fixed RY/CNOT chain + ring, L = 1 | n (abelian) at every n = 4 to 11 | `results/fixed/n<n>_L1/dim` |
| fixed, n = 4, 5, 7, 8, 10, 11, L >= 2 | dim so(2^n): 120, 496, 8128, 32640, 523776, 2096128 | `results/fixed/n<n>_L2/dim`, `frac_of_so` = 1.0 |
| fixed, n = 6 | 510 / 1023 / 2016 = so(64) at L = 2 / 3 / 4 | `results/fixed/n6_L{2,3,4}/dim` |
| fixed, n = 9 | 32766 / 65535 / 130816 = so(512) at L = 2 / 3 / 4 | `results/fixed/n9_L{2,3,4}/dim` |
| pools V and G (Tang 2021, 2n - 2 strings), n = 4 to 9 | 36, 136, 528, 2080, 8256, 32896 = dim so(2^(n-1) + 1) | `results/pools/V_n<n>/dim`, `G_n<n>/dim` |
| pool L2 (all 1- and 2-local odd-Y strings), n = 4 to 9 | so(2^n) | `results/pools/L2_n<n>/dim` |
| ADAPT-selected sets, n = 7, alpha = 1 | abelian (7) at every step, both pools, both optimisers; L-BFGS stops at P = 7 with no operator selected | `results/adapt_sets/*a1.0*` |
| ADAPT-selected sets, n = 7, alpha = 0.25 | V: 7 -> 16; L2: 7 -> 1025 at P = 21 (12.6% of so(128)) | `results/adapt_sets/*a0.25*` |

The pre-registered prediction H2b, that dim(DLA) at the deployed (n = 7, L = 3) is below 8128
(guess 4095), is falsified: it is 8128, the full so(128), already at L = 2. H2a (depth 1
abelian), H2c (the pools: 2080, 8256, 32896 predicted and measured), H2d (odd-Y) and H2e (ADAPT
sets abelian at alpha = 1) held. That n = 6 and n = 9 need depth 4, with dimensions at L = 2 and
3 equal to 2 dim su(2^(n-2)) and dim su(2^(n-1)), is an observation from two cases labelled
HYPOTHESIS, not explained.

What it means for V.8: the algebra at the deployed cell is maximal, so nothing in the algebra
protects the ansatz from an exponential plateau. Once a circuit is a 2-design over exp(g), the
variance scales as 1/dim(g) with dim(g) = 8128 at n = 7 and growing as 4^n/2; S13 measured the
decay base approaching 0.504 per qubit at depth 8 (`s13/results/geo_kernel.json`), the 2-design
rate. The S25 result "no exponential plateau at n = 4 to 13" is therefore a statement about
depth 3 (P = 3n against dim so(2^n)) and is quoted with "at depth 3" attached. It is not
evidence of a favourable algebra, and small-DLA simulability arguments do not apply: the circuit
is simulable because n = 7, not because of its structure (`s21/LEDGER.md` L4). On the pools: a
"complete" pool in Tang's sense (overlap-matrix rank 2^n - 1) generates so(2^(n-1) + 1), which
acts transitively on the real sphere and has about a quarter of the dimension of so(2^n);
completeness is weaker than controllability, so the ADAPT arm using pool V in A1 is restricted
to that subalgebra by construction and the L2 arm is not.

### V.10 The S13 locality theorem and the Pauli spectrum

This is the trainability record of the torsion-space generation lane
(`s13/SPRINT13_DOSSIER.md` sections 6 to 11; `docs/CONDENSED_REPORT.md`, "What survives").

**An exact locality theorem in torsion space.** Under an ideal-geometry backbone builder, the
CA-CA distance d_ij depends on exactly the j - i - 1 residues strictly between i and j,
contiguous; proven on more than 5,000 (pair, variable) cells with agreement 1.0000 and zero
counterexamples; non-supporting variables move it by exactly 0.000 A; it holds at chain
termini, on glycine (74 targets) and proline (52), and under cis-omega. Mechanism: d_ij is an
internal coordinate of the CA_i to CA_j sub-chain and the torsion-independent virtual bond
removes the two end residues. The all-atom case follows four exact rules (CA: i < m < j; N:
i <= m < j; C/O: i < m <= j; CB: i <= m <= j), max |delta d| outside support 1.4e-13 A against a
minimum inside of 4.6e-2 A (`s13/qarch_FINDINGS.md` section 1). Consequences: a separation-8
pair is a 14-qubit interaction at k = 4 bits per residue, so no 2-local Ising form of a
distance-based molecular objective exists in this encoding; and the union of supports over all
pairs is the whole chain, so both energies are full-register and "AMBER is less local than
Legacy" is a category error, recorded as refuted.

**The artefact, and the correction.** The first Pauli-spectrum measurement reported AMBER's
mean Pauli weight above Legacy's on 14 of 14 cells. Over 141 fully enumerated tables the top-10
configurations of 4,096 carry a median 99.6% of raw AMBER's Walsh variance; a constant plus a
single spike has Walsh weight spectrum exactly Binomial(m, 1/2) with mean m/2, measured 6.001
against predicted 6.001, L1 distance 0.0003. A delta spike is maximally global for arithmetic
reasons; the measurement was of steric clashes (Part IV.5). Some Legacy cells are spiked too
(top-10 share 0.878 on one), and 99th-percentile winsorisation is not sufficient conditioning
(a winsorised AMBER table still carried 0.986 of its variance in ten configurations); only
monotone rank-preserving conditioning works (`s13/walsh_FINDINGS.md`).

**The corrected result.** With rank-preserving soft compression applied identically to both
models: mean Pauli weight Legacy 2.236, AMBER 3.015, AMBER higher on 79 of 79 cells; share of
variance at weight <= 2, 0.641 against 0.392; >= 4-body Sobol share 0.069 against 0.122 (12 of
13 cells); cumulative 1 + 2-body share 0.778 against 0.629. Neither model is 2-local, and the
non-locality is in interaction order, not sequence range. The chain closed with no free
parameter: predicted gradient variance from `sum_S c_S^2 Var_theta[d<Z_S>/dtheta_i]` against
measured, median ratio 1.006 (Legacy) and 1.001 (AMBER, conditioned) on the exact per-string
form over 95 cells; the coarse weight-kernel form (1.278) is withdrawn, and raw AMBER's ratio
(0.913) is not a converging estimator (its across-theta sd rises with the number of samples,
0.15 to 4.18). Per term, AMBER's non-bonded term has covariance share 1.000 on all 26 component
cells and Legacy's steric term 0.955 at m = 18: both models are a steric potential plus
rounding error. The spectrum saturates: Legacy's mean weight goes 3.89, 4.10, 4.17, 4.27 across
m = 12 to 18 with the tail beyond weight 3 flat at 0.60 to 0.64, an effective interaction order
of about 3 residues independent of length. The residue-order spectrum is encoding-invariant to
1e-16, but the same Legacy energy has mean qubit weight 2.70 / 10.00 / 1.27 under binary /
one-hot-mean / one-hot-penalty encodings, a factor of 7.9: a one-hot Hamiltonian's Pauli
spectrum is a free parameter of the implementer.

**Trainability, the honest negative** (section 10). The measured ansatz kernel v(w) is flat in
Pauli weight (median v(w_max)/v(1) = 1.00) and decays in n as 2^(-0.47n) to 2^(-0.86n); the
standard barren-plateau theorem needs local 2-design blocks this ansatz does not have, so it
licenses no prediction at 6 to 18 qubits. Exponential and polynomial decay both fit at R^2 0.93
to 1.00 over n = 6 to 14 and are not discriminable. The decay base falls monotonically with
depth onto the 2-design limit, 0.823, 0.768, 0.648, 0.537, 0.504 at depth 8; at depth 2 the
kernel is not flat (0.862 to 0.078 over n = 6 to 14), so cost-locality would explain behaviour
at depth 2 and does not at the depth the project uses.

**Three clean results** (section 11). The metric (quantum Fisher information) contains no
Hamiltonian: bit-identical across energy models at matched theta (0.000e+00), now a unit test;
the energy model selects which region the optimiser visits, it does not reshape the manifold.
Quantum natural gradient has nothing to fix: the metric is full rank at every theta, n and
depth, g_ii = 0.2500 exactly, off-diagonal correlations 0.008 to 0.037 and shrinking with n,
exactly I/4 at depth 1. CVaR at small alpha is exactly a steric clash filter: for alpha <= 0.25
`amber` and `amber_soft` are the same objective to every printed digit while differing by
1.7e16x in gradient variance at alpha = 1. And SPSA optimises the AMBER objective best of five
arms (percentile 0.088) while returning the worst structure (+0.333 A against random): better
optimisation of a misaligned objective produces worse physics.

### V.11 The quantum record across the sprints

- S5: the CVaR gradient defect found and priced (`docs/FINDINGS.md:368`); a VQE over segment
  choices with its encoding and limits stated (`docs/FINDINGS.md:404`).
- S6: VQE/CVaR assembly gives a narrow pool whose narrowness cannot be exploited
  (`docs/FINDINGS.md:1021`).
- S8-9: CVaR-VQE over the discrete hypothesis set, and the integrated four-component system
  priced component by component (`docs/FINDINGS.md:3240`, `:3299`); the gradient audited against
  two references (`:3147`).
- S9-5: refinement's failure located in selection; the four-component ablation conclusive; the
  gradient defect fixed and priced against an exact reference (`docs/FINDINGS.md:4154-4213`).
- S12: the mandated quantum component, "a verified problem, and a classical win"
  (`s12/SPRINT12_DOSSIER.md` section XVIII): 20 of 20 alpha x T arms lose, and the VQE's eight
  most probable states never beat the eight best classical ones.
- S13: the locality theorem, the Walsh artefact and its correction, the Pauli-spectrum to
  gradient-variance chain, the metric, QNG and CVaR results (V.10).
- S14: a VQE redesign around a structural objective whose optimum is in the right place while
  the energies' is not (`s14/LEDGER.md`; `docs/STATE_BRIEF_2026-09-12.md` section 7).
- S15: the paper-driven programme; `s15/VQE_CVaR.md`, `s15/qens_FINDINGS.md`,
  `s15/qgeom_FINDINGS.md`, `s15/qrestraint_FINDINGS.md`.
- S16: every CVaR-VQE arm loses to uniform random sampling at matched budget (+0.209 to
  +0.301 single window, CIs excluding zero) on the coordinate-average readout as well
  (`s16/LEDGER.md` L17); there is no phase boundary and the axis is not rho (L25); the
  programme's last standing positive quantum result, the unranked ensemble, is retired: both
  classical searches beat the VQE on the ensemble readout too, +0.27 to +0.62 A (L29).
- S17: every matched-budget comparison in three sprints was run past classical saturation
  (`s17/LEDGER.md` L22).
- S18: the quantum branch of the degree-1 objective is closed, the answer forced before the arms
  were run (`s18/LEDGER.md` L6, L16).
- S20: the quantum sampler hypothesis is refuted (`s20/LEDGER.md` L1); every circuit-side
  landscape metric is a difficulty proxy and the CVaR tail is worth nothing (L4); the AMBER
  objective is defined only after relaxation (L6).
- S21: the encoding lever is confounded by step count and closes (`s21/LEDGER.md` L2, L13, L30,
  L37, the last showing both of S20's "significant" cells null); no bond dimension breaks
  classical simulability because the register is too small (L4); the exhaustive latent argmin
  does not beat a zero-evaluation pool (L14, L17); the ansatz ladder on AMBER is won by
  best-of-N and chi orders nothing (L34); MDE is per comparison (L8).
- S22: the entropy-regularised tail does not beat the classical selector (`s22/LEDGER.md` L13);
  the readout-H effect replicates on a different substrate (L14) and is the same fact as the
  A2b null (L16); multi-stage VQE helps the tail a little and beats no classical bar (L15).
- S23: the probability-weighted readout fails; the last quantum-side door is closed
  (`s23/LEDGER.md` L8).
- S25: the +0.113 A CVaR contribution withdrawn (L3, L5); the readout insensitive to nearly half
  the mass (L15); the Hamiltonian barely changes between targets (L17); `s25/QUANTUM.md`.
- S26: A2, the dynamical Lie algebra (V.9); the rest of Proposal A fills Part VII as it lands.

The position the record supports: the deployed selector is an exactly simulated, exactly
differentiated, correctly implemented CVaR free-energy optimiser over a 128-state diagonal
Hamiltonian that is nearly the same on every target; it trains, does not reach its cheap
classical optimum, and is read out by an operator that cannot tell the difference; no
configuration of it has beaten a matched classical control on any readout; and the part of the
quantum work that survives as a result is the trainability measurement of V.10, with V.9's
correction to its scope.

<!-- PART V S26 ADDITIONS -->

## PART VI. THE TWENTY-TWO SPRINTS AS A STORY

Sprints 1 to 4 predate the record kept here; their tree is in git history at `5fa05cd` and their
one surviving claim (that physics ranks real geometry) was re-measured and overturned in S5 and
S7. From S5 the record is `docs/FINDINGS.md` (S5 to S11, with its corrections ledger at lines
60 to 118 and its index at 119 to 211), the sprint dossiers `s12/` to `s15/`, and the live
ledgers `s14/LEDGER.md` to `s26/LEDGER.md`. `docs/CONDENSED_REPORT.md` condenses S5 to S13 and
`docs/STATE_BRIEF_2026-09-12.md` section 7 gives one line per sprint. Each section below has the
same six fields. "Falsifier" names the pre-registered condition where one existed; formal
pre-registration began in S15, and before that the field names the test that decided the
sprint. Numbers are quoted from the named section or ledger entry, which names its artefact;
Appendix B lists them as cited.

### VI.1 Sprint 5: the wall is located (`docs/FINDINGS.md:123`, sections at 248 to 780)

Question. A K = 500 retrieval pool holds a 1.71 A member on average and the shipped selector
returns 3.45 A (both single window): where between the two is the accuracy lost, and can a
larger pool, more search, or a better sequence key recover it?

Falsifier. Paired arms on the 126 targets against the shipped baseline; an arm that did not move
the mean by more than its interval was a null.

Result. The pool was never the problem (section 1) and search is not the constraint at any
budget (2); the loss is a 3.3 A wall (3) which is not "no signal" but no resolution at the top
of the ranking, the shortlist ceiling (4). The shipped CVaR gradient carries a defect (6, Part
V.5). In-band skill is bounded by the prior's own accuracy (7); winner's curse is real and small
(8); structure and sequence are decoupled at this length, so retrieval by sequence cannot work
as a fix (9); a torsion prior is the sixth signal to hit the same wall (10); the causal chain
closes and the prior cannot be improved from sequence, all four arms (11); retrieval by
predicted structure buys better selection and a worse pool with no net gain (13); the synthesis
is null too and its interim numbers were a mirage (14).

Closed. Pool size, search budget, sequence-key retrieval, retrieval by predicted structure.

Retracted. Section 12, "physics does rank real geometry, Sprint 4's measurement was
confounded", was itself overturned in S7 (finding 10) on matched pools with a radius-of-gyration
control: AMBER is worse, not better.

Open. In-band discrimination, carried to S6.

### VI.2 Sprint 6: selection is measured out (`docs/FINDINGS.md:140`, sections at 781 to 1331)

Question. Can any ranker, learned or physical, order the near-native band of a pool?

Falsifier. Properly powered paired comparisons at n = 126 with a correctness audit of every arm
before it was reported.

Result. The benchmark does not have a floor; the ceiling is ours (1). A learned ranker fixes
the global ranking and not the in-band one (2); a correctness bug in the pairwise arm was found
before it was reported (3); properly powered, the learned rankers are exactly null (4); building
a structure instead of picking one is null (5); secondary structure is not the discriminating
axis (6); tail aggregation of pair violations is null (7); VQE/CVaR assembly gives a narrow pool
whose narrowness cannot be exploited (8); refinement's two hypotheses are answered and the
second explains the first (9); the full ablation yields one usable confidence signal,
inter-generator agreement (10); anchored selection is null (11); decomposed all-atom physics in
a learned combiner is null and its +0.320 did not replicate (12); the distance channel is
sufficient, the predictor is not (13).

Closed. Learned in-band rankers; tail aggregation; anchored selection; physics in a learned
combiner.

Retracted. The +0.320 combiner effect (section 12).

Open. The inter-generator agreement lead (closed in S7 as a difficulty thermometer).

### VI.3 Sprint 7: the predictor is the constraint (`docs/FINDINGS.md:155`, sections at 1332 to 2088)

Question. Is the distance predictor improvable, and does anything rank the truth first?

Falsifier. Each lead carried a control that would show it to be difficulty rather than skill.

Result. The agreement lead is a difficulty thermometer (1); thirty times more structural data
makes the predictor monotonically worse (2); the predictor is not shrunk but uncorrelated, so
post-hoc de-biasing is dead (3); chirality is a red herring on the generative path (4); the
objective is misspecified, the native is not its minimum (5); pool size has an interior optimum
and selection tracks the pool mean, not its best (6); training against ranking makes it worse,
and the predictor is approximately typicality (8); the triangle operator improves the
diagnostics and not the outcome (9); no scoring channel ranks the truth first, AMBER is worse
not better, final at n = 70 (10); ESM does beat one-hot on selection, the recorded ablation
having judged by the wrong metric (11); BLOSUM retrieval beats random fragments, the recorded
claim refuted (12).

Closed. Post-hoc de-biasing; more training data; chirality; physics as a ranker.

Retracted. "ESM adds nothing" and "random fragments beat retrieval" (sections 11 and 12).

Open. The aggregation metric, carried to S8.

### VI.4 Sprint 8: the full architectural attack (`docs/FINDINGS.md:169`, sections S8-1 to S8-14 at 2089 to 3899)

Question. With the ranker line closed, does any architectural change (representation, pool
construction, retrieval key, generation, synthesis, physics) move the endpoint?

Falsifier. Every arm paired at n = 126 with an ORACLE ceiling computed first; a single
pre-registered dev24 pass for the survivor.

Result. The representation is not the bottleneck, the aggregation metric is (S8-1);
multimodality is modest, two to three basins (S8-2); the 3.454 A mean is four classes and a 28%
opportunity, of which the coordinator corrected that class contributions are not recoverable
fractions (S8-4); pool construction is capped at 2.406 A single window at any filter skill, and
the law that says why (S8-5; later corrected: measured through a weak terminal, averaging the
same pool gives 1.925, `docs/CONDENSED_REPORT.md:176`); retrieval is fixable and fixing it is
worth 0.016 A (S8-6); inverse folding is the least misspecified objective and still misspecified
(S8-7); inside the near-native band the score is worse than a coin flip and consensus is the
only discriminator, giving the first arm in three sprints with a CI excluding zero, a score
filter plus consensus medoid (S8-8, `docs/FINDINGS.md:2893`); the channels' errors are
decorrelated and it buys 0.01 A (S8-9); generation is worth nothing as candidates but exposes a
second law (S8-10, corrected by S8-13: the transfer law is the identity map where we operate);
synthesis beats selection at 3.204 A, replicated on dev (S8-11); AMBER is a validity stage, not
an accuracy stage (S8-12); the training target is misspecified by 0.5 A per pair and fixing it
is worth 0.044 A (S8-14).

Closed. Pool construction, retrieval key, generation as candidates, learned combinations of
signals, K = 25.

Retracted. The S8-10 transfer law (S8-13); the class-contribution arithmetic (coordinator
correction at 2612).

Open. Synthesis, taken to the benchmark in S9.

### VI.5 Sprint 9: the benchmark pass (`docs/FINDINGS.md:186`, sections S9-1 to S9-10 at 3900 to 4607)

Question. Does the synthesis architecture's tuning gain transfer to the sealed benchmark?

Falsifier. Pre-registered constants; one pass; the paired difference's CI.

Result. The shared bias is real, one interpretable mode, and not an error (S9-1); a torsion
prior buys physical validity for free and no accuracy (S9-2); the relational channel is
typicality wearing a tournament (S9-3); the loop is a contraction toward its own fixed point
(S9-4); refinement fails and the failure is all selection (S9-5, with the gradient defect fixed
and priced); synthesis is at its ceiling and the binding constraint is the distance target
(S9-6); evolutionary information exists, is retrievable, and does not help (S9-7); the in-band
ceiling is informational, not architectural (S9-8); consensus is the wrong operator for a
refinement trajectory (S9-9); and the final benchmark pass: 2.9610 against 2.9507 A single
window, +0.0103 [-0.1596, +0.1803], 31W/29L; the gain does not replicate, the concentration was
the tell, and the leakage audit found a real defect (S9-10).

Closed. The synthesis line as an accuracy claim; the benchmark (spent).

Retracted. The tuning gain as a transferable effect.

Open. Why the non-replication (S10).

### VI.6 Sprint 10: adjudication and bounds (`docs/FINDINGS.md:198`, sections S10-1 to S10-5 at 4608 to 5159)

Question. Was the non-replication contamination or concentration, and where is the remaining
headroom?

Falsifier. An attrition ledger over the 126 targets and a cross-tabulation of the identity
audits.

Result. The loss is set purity and the failure is a minority class (S10-1; its own identity
audit withdrawn); the projection path is closed because the required matrix is 2.7x better than
anything achievable (S10-2); physics does not supply decorrelated error, the force field agrees
with retrieval target-specifically (S10-3); the identity-audit conflict is adjudicated as an
alphabet permutation and the leak priced (S10-4; the artefact it cites is absent,
`s26/EXAMINATION.md` C27, so the price is not quoted here); the reconciled bound ladder shows the
geometry is not the barrier (S10-5).

Closed. Contamination as the explanation; the projection path; physics as a decorrelated
channel.

Retracted. S10-1's identity audit (void, S10-4).

Open. Recognition, named as the barrier.

### VI.7 Sprint 11: performance engineering and consolidation (`docs/FINDINGS.md:205`, S11-1 at 5160, S11-2 at 5219)

Question. Do the four mandated components run, and can the codebase be consolidated without
changing a number?

Falsifier. A paired ablation of all four components at n = 126, and a two-arm equivalence test
between the consolidated `core/` and the reference modules.

Result. The four components run and three buy nothing measurable (S11-1); the projection is not
a function of its input at Angstrom resolution, the reference disagreeing with itself by 1.62 A,
which is the degenerate-branch fact of Part III.7 (S11-2). The consolidation into `core/` kept
the reference arm on disk as the equivalence oracle (`README.md`, "Where everything is";
`docs/consolidation-2026-09-04.json`).

Closed. The performance question; the reference arm frozen.

Open. The last untested ranker class, learned aggregation over the deviation map (S12).

### VI.8 Sprint 12: adversarial, aggregation, assembly (`s12/SPRINT12_DOSSIER.md`; thirteen findings files in `s12/`)

Question. Is in-band discrimination signal-limited or sample-limited, and what does the terminal
operator consume?

Falsifier. A set-transformer over the signed (75 x n_pairs x 18) deviation tensor with the 18
hardest targets held out; a leaked-label harness as the positive control.

Result. The learned set decoder emits 3.184 against production 3.203 (d = -0.019 [-0.058,
+0.020], 63W/63L), with a flat learning curve 3.043 to 3.026 from n = 8 to 75 while the leaked
label reaches 2.534 at n = 8: signal-limited, not sample-limited (section I). The terminal
operator consumes the set mean, not the set best: d_out = 1.16 d_set_mean + 0.04 d_set_best,
R^2 0.893 over 882 perturbations, so a perfect rank-1 decision is worth -1.743 A through argmin
and -0.029 A through the m = 75 average (section III; `docs/CONDENSED_REPORT.md:61`). Sequence
conditioning is worth 1.004 A on the 108 ordinary targets and negative on the 18 hardest, where
a pipeline given random windows emits 5.425 A against the shipped 6.019 A; three routers are
null (section I). Ten of the eighteen failures are fibril segments and lasso peptides, 10/18
against 6/108, Fisher p = 1.2e-6 (section IV). Multi-piece assembly is closed. The mandated
quantum component is a verified problem and a classical win: 20 of 20 alpha x T arms
(section XVIII). Section XVII priced one direction that reaches the target, chemical-shift
torsion restraints at 1.486 A.

Closed. Learned aggregation; multi-piece assembly; the ranking line at every level.

Retracted (later). The 1.486 A restraint figure and its coverage gate (S14 C14).

Open. Torsion space and the quantum architecture (S13); chemical shifts (S14).

### VI.9 Sprint 13: torsion space and the quantum architecture (`s13/SPRINT13_DOSSIER.md` sections 1 to 16)

Question. Does a torsion-space representation with a molecular energy give a VQE something to
optimise, and what does the energy model do to trainability?

Falsifier. Full enumeration of nine targets at 262,144 configurations each; the certified global
optimum against random sampling; leave-fold-out torsion prediction.

Result. The representation is not the barrier: at about 24 live qubits the space contains a
1.594 A answer, but 88% of what the library buys is generic Ramachandran and 0.388 A of the
ceiling is the privileged oracle start (1). No native-free objective finds it (2); neither
energy ranks the native, Legacy's in-decile rank correlation +0.043 and raw AMBER's -0.088, the
native at the 32nd to 40th percentile, and Legacy's certified optimum +0.139 A worse than random
(3). Optimising harder makes the structure worse, 3.764 at 10 evaluations, 3.667 at 300, 3.920
at the certified optimum (4). Sequence-only torsion prediction is dead: full context predicts
phi no better than a sequence-blind marginal, the whole channel is 10.4 degrees of psi, and the
best sequence-only builder emits 3.770 A (5). The locality theorem, the Walsh artefact and its
correction, the corrected Pauli spectrum with its no-free-parameter gradient-variance chain, and
the metric, QNG and CVaR results are Part V.10 (6 to 11).

Closed. Torsion-space search on either energy; sequence-only torsion prediction; QNG; the
"AMBER is less local" framing.

Retracted. "AMBER's weight exceeds Legacy's on 14/14 cells" (the delta-spike artefact);
"optimising Legacy is worse than not optimising" as a general statement (reverses under
stronger optimisers); the energy-ordered Gray-coding gain, corrected to -10.9% against its own
null.

Open. A native-free structural objective (S14); chemical shifts (S14).

### VI.10 Sprint 14: chemical shifts and the VQE redesign (`s14/LEDGER.md`; `s14/SPRINT14_DOSSIER.md`)

Question. Can chemical-shift-derived torsion restraints reach 2.0 A at real coverage, and does a
structural (non-energy) Hamiltonian have a good low-energy region?

Falsifier. Coverage measured from deposits before any structure code (C14); certified
enumeration of the structural objective (C13).

Result. The chemical-shift route is closed by arithmetic: 54 of 126 targets are runnable and
ORACLE-perfect torsions on all of them still leave the instrument at 2.021 A, 55 being needed
(C14); the durable outputs are that a bimodal shift posterior is one qubit with a physical
justification (-2.253 A [-2.642, -1.865] for a search over the top-8 support against the
argmax, C14b) and that confidently wrong costs two to three times what absent costs (C14c). The
structural objective's certified global optimum is 0.885 A better than random where Legacy's is
0.139 A worse (C13); the distogram orders the bulk and the torsion prior places the optimum
(C13b); the selection gap survives certification across disjoint target sets (C13c) and a
0.68 A gap survives infinite budget (V2). Searching a good objective harder buys nothing (C8);
aggregating a low-energy set in coordinate space recovers 0.29 A and still loses to the
incumbent (+0.110 [+0.004, +0.214], C9); the ablation ladder closes with aggregation worth 3.4x
the objective and the objective's whole contribution 0.171 A (C11). `leg_torsion`, the sprint's
strongest open lead, is refuted at n = 126 (+1.729 A, C12). Legacy is a clash gate (98.8% of its
variance is `steric`, E1); AMBER is a validator and an inert refiner (+0.014 A [-0.001, +0.028],
E2); no physical objective exceeds 0.511 pairwise accuracy below a 0.25 A quality gap (E3). The
VQE is the worst optimiser tested on its own axis (V3); two new CVaR defects (V4); plain
expectation-value VQE returns the best structure, low alpha buys entropy not accuracy (V5).

Closed. Chemical shifts as a route to 2.0 A; `leg_torsion`; term reweighting; AMBER refinement
as rescue; the log-qubit continuous-torsion encoding as compression (it is a
reparameterisation, C5c).

Retracted. The ledger's own REFUTED table lists twenty-two items, among them the coordinator's
"eight times better ordered" headline (corrected under a uniform proposal), the ENERGY lane's
`leg_torsion` and reweighting claims, the VQE lane's Gray-coding claim, and an 86% bond
contraction (a frame bug; corrected to 16 to 19%).

Open. N1 to N6 of the ledger: the per-target skill and compactness discrepancy; a
pre-registered restraint constant for the -0.022 A; anything supplying the per-target sign of
the in-band ordering; a forward chemical-shift predictor; QNG at depth >= 2; the
generator-shift test.

### VI.11 Sprint 15: the paper-driven programme (`s15/LEDGER.md`; `s15/FINAL_DOSSIER.md`; `s15/PAPER_DRAFT.md`)

Question. Can a generative torsion distance-geometry architecture, with the instrument frozen
and audited first, beat the incumbent; and what do the mechanism experiments say about the
objective?

Falsifier. Phase 0 audit rows 0.1 to 0.15 as a gate; every claim traced by an adversarial
consistency audit; pre-registration from here on.

Result. The gate passed with one retraction, one correction and one new defect: the peptide
database and all 126 window universes rebuild bit-identically; `pool_best` reproduces from PDB
files as 1.7108244199364904 exactly; an independent RMSD implementation agrees to 2.04e-13 A over
63,000 structures; AMBER parameters are bit-exact against an independently built force field;
AMBER's single point costs 8.3 to 23.3 ms with memoisation defeated, not 28; `BAND = 1.5 A` has
no derivation; 16 of 126 targets carry verbatim own-fold windows with zero measured impact on
`pool_best`. The generative architecture: torsion distance geometry reaches an ORACLE 0.611 A
and a predicted 3.644 A, +0.440 [+0.290, +0.592] against the incumbent, a loss (1.3); the
distogram has MAE 2.386 A, bias +0.509 rising to +1.492 at separation 11 to 15, z-sd 2.633
(1.4); the retrieval pool is a second distance channel with the opposite bias sign (1.9).
Mechanism: the native sits at the objective's 34.7th percentile and is its argmin on 4 of 126,
yet the argmin is 0.949 A [0.753, 1.147] better than random (2.1, 2.2); the full cascade does
not beat the incumbent (+0.117 [+0.050, +0.188]) because filtering the ensemble by the objective
makes aggregation worse (2.6, 2.7). Workstreams: the literature survey (67 papers) found two
claims already taken and one control without precedent; the quantum positive was demoted (27
"cells" were 9 targets x 3 seeds; null at target level); the AMBER positive was weakened on
accuracy and confirmed on validity; the adversarial audit traced 178 claims, 141 matched, 21
mismatched, 16 untraceable, 13 blockers applied. The state brief's line: 3.321 A (worse); error
shape beats magnitude; a 0.611 A true distance-geometry floor.

Closed. The generative torsion architecture as an accuracy route; the pool channel as an
objective (worst arm, 2.4); Family B feasibility (2.5).

Retracted. The 8-target sign of the distance-geometry result (superseded at n = 126); the
"1.06 A quantisation ceiling" (an artefact of a uniform-grid lookup on a non-uniform grid, 1.6,
1.7); `hash(pdb)` seeding (salted per process, 1.10); the running tally is in the ledger.

Open. The mechanism experiments in flight at close (2.8 to 2.12); native-free steering (S16).

### VI.12 Sprint 16: native-free steering (`s16/LEDGER.md` L1 to L29)

Question. Can the pool's error direction be estimated without the native and steered against?

Falsifier. The flagship's rules were written before any measurement (L1, L2) and decided it
(L9); the pivot carried a stronger falsifier (L13).

Result. The flagship is refuted at n = 126 (L9) and the coordinate-space pivot is refuted too
(L13); the audit found the flagship's premise about five times weaker than the record said
(L14); the fusion law is a theorem and the claimed 0.162 A was two of our own errors (L16);
integration found Legacy contributing nothing, every CVaR-VQE arm losing to uniform random
sampling at matched budget, and an AMBER contribution that was then withdrawn (L17, L21, L23);
diversity-aware selection is a null (L20); VERIFY found three of the coordinator's five claims
failing the same way (L21); QPHASE found no phase boundary and the axis is not rho (L25);
REVIEWER found the report breaking its own rule inside the section that states it (L26); REPAIR
removed the last accuracy claim in the programme and recorded the rule that a control must
match the operator's space (L27); the last standing positive quantum result, the unranked
ensemble, is retired (L29).

Closed. Native-free steering in both spaces; every arm chose "do nothing"; the AMBER relaxation
gain dissolves under a matched-random control.

Retracted. The 0.162 A fusion effect; the AMBER-as-ranker interval; the unranked-ensemble
quantum positive.

Open. The averaging-to-selection readout (S17).

### VI.13 Sprint 17: averaging to selection (`s17/LEDGER.md` L1 to L30)

Question. Is the shortlist or the ranker what binds, and does widening K help?

Falsifier. The ORACLE map at n = 126 (L15), coverage-preserving shortlists (L14), and matched
budgets audited for saturation (L22).

Result. The K = 500 "retrieval ceiling" is a truncation artefact (L1); the distance objective
cannot be re-engineered into a better selector (L5); the in-band problem is information-limited
but the information is present and mis-routed (L6); consensus reproduces across a decade of
instruments and typicality does not (L9); the in-band result at n = 126 (L10); widening K hurts
because the shortlist loses the answer (L12); coverage-preserving shortlists recover the ceiling
and the readout cannot use it (L14); the averaging readout has no hidden headroom (L17), except
that at n = 126 the averaging ceiling does beat the best-member ceiling (L18); every
matched-budget comparison in three sprints ran past classical saturation (L22); the last feature
class is closed and the shortlist binds (L23); the validity frontier splits in two and the
incumbent's restraint set wins (L24); the one target-level signal that works is the trivial one
(L25); the decisive AMBER-against-Legacy comparison at full scale (L29); the objective is not
weak, it is wrong (L30).

Closed. Selection at four levels; widening K; feature classes.

Retracted. L3 (the coordinator's own error, caught by SELECT); the "routing ceiling" as a
min-of-N statistic below its noise floor (L8); "steric validity is free" quoted against the
wrong baseline (L27).

Open. The degree-1 objective (S18); why the errors are maladaptive (S19).

### VI.14 Sprint 18: the degree-1 objective (`s18/LEDGER.md` L1 to L17)

Question. Is a degree-1 (single-pair) truncation of the distance objective as good as the full
one, and does a quantum branch on it have anything to optimise?

Falsifier. Five pre-registered falsifiers F1 to F5; F2, F4 and F5 fired (L5).

Result. The 19-target result reproduces exactly and was never significant (L1); the functional
form is sound and the distogram's error is worse than noise (L2); both confound controls survive
(L3); the degree-1 branch is closed (L5) and the quantum branch with it, the answer forced
before the arms ran (L6); the degree-1 object is the torsion prior (L7); the prior arm is refuted
with its falsifier fired (L9); 1/sd^2 weighting is validated native-free (L10, amending L4);
the lambda ladder at n = 126 shows degree-1 decisively worse than the full objective (L11); the
deployed tempering exponent is at the optimum (L14); `leg_contact` fails as a term (L15) and is
harmful, and every physics filter loses to a random gate (L17).

Closed. Degree-1 truncation; the quantum branch on it; `leg_contact`; physics filters.

Retracted. `shuf_paired` priced the weights, not the error assignment (L4); the helix-mu
control is not a match and is NOT MEASURED (L13, correcting L12).

Open. The mechanism of the maladaptive error (S19).

### VI.15 Sprint 19: why the errors are maladaptive (`s19/LEDGER.md` L1 to L20)

Question. What is the pool's systematic error, and can it be estimated native-free?

Falsifier. The opening hypothesis carried its mechanism; the P5 headroom test with an oracle
(L7); the soft-consumption falsifier (L16).

Result. The opening hypothesis is refuted in its mechanism (L1); the harm is the residual's
cross-pair sign correlation (L3); the in-manifold mechanism is measured and deliberately not
promoted (L6); the branch has no headroom even with an oracle (L7); the harmful coherent
component is shared and two thirds of it is reproduced by zero-information references (L11);
score gates damage an averaged set by alignment, not diversity (L12); the harmful mode is named
and is not estimable native-free, an identifiability wall (L14); soft consumption's falsifier
fires (L16); a scalar validity score would have promoted a broken structure (L18); "best built
3.048 A" is a contracted point cloud, not a structure (L20).

Closed. Native-free estimation of the harmful mode; the distribution-shape direction; soft
consumption.

Retracted. The "pathological AMBER minimiser" was CPU starvation (L17); a clip-induced magnitude
confound in the headline (L15); L2's "the modes carry real positional information" downgraded
(L5).

Open. Landscape geometry and Legacy against AMBER (S20).

### VI.16 Sprint 20: landscape geometry (`s20/LEDGER.md` L1 to L15, L-B)

Question. What is each energy's landscape like, and is a quantum sampler worth anything on it?

Falsifier. F-D2 (the relaxation clause, L6) and F-C2 (the validity score, L15) both fired.

Result. The quantum sampler hypothesis is refuted (L1; the two positives that survived it were
nulled by S21 L37); the +0.164 A "AMBER tax" is the projection's (L2); every circuit-side
landscape metric is a difficulty proxy and the CVaR tail is worth nothing (L4); the
shared-referent floor turns "two thirds sequence-independent" into one fifth (L5); the
relaxation is what makes the AMBER objective defined (L6); Legacy is a compactness model and the
two potentials disagree about ordering (L8); the programme's paired CI was iid over targets, not
fold-clustered (L9); the Frechet mean recovers nothing (L11); AMBER's difficulty is a steric
singularity and 72% of its apparent damage is move size (L12); corpus change does not
decorrelate, the carrier is the selector (L13); the whole optimisation question is priced at
zero for the first priority (L-B, L14).

Closed. The sampler hypothesis; the Frechet route; the optimisation question.

Retracted. Five live claims damaged by workstream D (L7); the "two thirds" framing (L5).

Open. CVaR-VQE as a selector on the discrete latent (S21).

### VI.17 Sprint 21: CVaR-VQE as selector (`s21/LEDGER.md` L1 to L37)

Question. Does the encoding, the bond dimension or the Hamiltonian give the selector anything,
and what is the exhaustive latent worth?

Falsifier. The exhaustive latent control (L12, L14, L17); D1's direct replication (L35); the
gauge sweep (L30).

Result. The encoding lever is confounded by step count (L2) and closes (L13, L37, both of S20's
"significant" cells null); "which Hamiltonian is best" is not well defined without a readout
(L3); no bond dimension breaks classical simulability because the register is too small (L4);
"MDE = 0.084 A" is not a property of the instrument (L8); Legacy carries no in-band rank
information beyond the distogram (L9); the budget exceeds the entire latent on 60% of targets
(L11) and the exhaustive argmin does not beat a zero-evaluation pool (L14, complete at n = 126,
L17); the answer is in the objective's top 5% and a top-512 readout ceiling is 1.986 A (L18);
the source is not the lever, the retrieval pool beats the latent through the same operator
(L20); all remaining leverage is in-pool selection, worth 1.33 A (L21), and the 1.338 A
selection gap is entirely uncaptured (L23); AMBER is harder, the continuation is degenerate in
raw units and preconditioning does not work (L24, L33); the sprint's first positive, a
compactness channel clearing where eleven functionals did not (L27), survives attack with the
mechanism being disagreement (L28); the "encoding" is a step-size reparameterisation and the bug
costs 0.015 A (L30); the ansatz ladder on AMBER is won by best-of-N and chi orders nothing
(L34).

Closed. The encoding question; bond dimension; the exhaustive latent; the Legacy-to-AMBER
continuation.

Retracted. R1 diversity-maximising set selection and the mechanism the lane had asserted (L7);
the bimodal route as demonstrated (L25); the coordinator's fork review carried two unchecked
operators (L32).

Open. Routing to a per-target m (S22).

### VI.18 Sprint 22: routing (`s22/LEDGER.md` L1 to L16)

Question. Is a per-target shortlist size real, and can anything native-free choose it?

Falsifier. The lane's own falsifiers on per-target m (L4) and hedging (L9) were stated and
fired in the directions recorded.

Result. The routing ceiling is 0.482 A and not concentrated (L1); the latent and the pool are
complementary with difficulty as the crossover (L2), but two native-free difficulty proxies fail
with the wrong sign (L3); per-target m is real and transferable, not noise-fitting (L4, L5, L6);
four independent router constructions fail, one significantly harmful (L7); hedging across m is
refuted (L9); sample size, not features, is the defensible explanation (L10); a third group's
external corroboration, and the tension not hidden (L11); in the extended ablation distance
alone wins and the torsion prior is the worst term (L12); the quantum lane closes, the
entropy-regularised tail not beating the classical selector (L13); the readout-H effect
replicates on a different substrate (L14), multi-stage VQE beats no classical bar (L15), and the
A2b null and the readout-H "lever" are the same fact (L16).

Closed. Routers (four); hedging; the quantum lane's tail regularisation.

Retracted. None recorded beyond the fired falsifiers.

Open. Driving the development instrument below 3.0 A (S23).

### VI.19 Sprint 23: drive dev below 3.0 (`s23/LEDGER.md` L1 to L11)

Question. Can the contraction be undone, the per-target scale predicted, or the readout
re-weighted to reach 3.0 A built chain?

Falsifier. Workstream A's three falsifiers (L5) and H5's setting sweep (L11), all fired.

Result. Kabsch superposition makes the incumbent's contracted output cost RMSD (L1); the
per-target scale is real and universal, but the global constant is worth nothing (L2), no
native-free feature predicts it (L3), and it transfers almost perfectly while being unreachable
in principle (L6); iterated Procrustes does not help and the contraction is intrinsic (L4); a
fifth router fails, and the signal that fails is a new kind (L7); the probability-weighted
readout fails and the last quantum-side door is closed (L8); the pool's error is 68% shared,
which explains every null in the sprint (L9); the mechanism's forward prediction has the right
direction and no measurable size (L10); H3 is closed by source audit and H5 repair strictly
costs RMSD at every setting (L11).

Closed. Scale correction; Procrustes; routers (five); probability-weighted readout.

Retracted. None.

Open. Learned candidate generation (S24).

### VI.20 Sprint 24: learned candidate generation (`s24/LEDGER.md` L1 to L16 with amendments)

Question. Can a learned generator supply candidates whose biases are independent of the
retrieval pool's, and what is the prior worth?

Falsifier. An ORACLE upper bound computed before the generator was trained (L7); the corpus
census (L8).

Result. The directive's screening statistic is broken, use the cosine (L1); the biases are
partially independent and quality still dominates (L2); the shared bias is the score's, not the
corpus's and not BLOSUM's (L3); a real leakage defect reaches the sealed benchmark, declared and
not quantified (L4); workstream B closes its own lane by an oracle upper bound (L7); the corpus
is not repetitive, it is small, with almost no beta sheet (L8, L8-A); the generation ladder at
full n is null-to-worse in every union and selection's lift is universal, +0.166 (L9, L9-A);
the functional lever closes in all three forms, AMBER selecting along a non-parallel direction
it cannot pay for (L10, L10-A, L16); the prior-attribution ladder is the only steep lever in the
project, -2.15 A per unit toward a perfect prior (L13, L13-A); chain correlation is 7.2% and the
signed-bias tension dissolves (L14); grid oracles are order statistics and only transfer
oracles survive (L15, L15-A).

Closed. Learned candidate generation; the functional lever; per-target grid oracles.

Retracted. L5's two primaries (L11, the mechanism claim withdrawn); L6's "fifth per-target
oracle" is an order statistic (L12); the +0.20 selection lift corrected to +0.166 (L9-A).

Open. The prior itself, and the basis question (S25).

### VI.21 Sprint 25: the endgame (`s25/LEDGER.md` L1 to L18; `ARCHITECTURE.md`; `s25/QUANTUM.md`)

Question. What is the production number, on what basis, and what does the record support about
the prior, the physics and the quantum component?

Falsifier. The round-trip gate on every exported structure (L10); the Phase I prior ladder
(L12); the seven-configuration suite against a random 75-subset (L16).

Result. The distance posterior is over-confident by a factor of two, measured (L1); calibrating
it makes RMSD worse (L2, later retracted in its mechanism, L7); the "+0.113 A CVaR contribution"
was measured at a temperature that does not ship, with alpha = 1 on three of five folds, and
was never a measured effect (L3, L5); the headline number described an object that is not a
protein structure and the basis decision is forced: built chain (L4, L8, L11); the shipped
score's per-pair target is quantised to 17 values (L6); the "86% accounted" null carries no
information (L9); the round-trip gate passes and the frozen build has a native-free leakage gate
(L10); a synthetic native-plus-noise leaderboard was found on disk and quarantined (L11); Phase
I closes: the score's target can be moved 25% toward truth and the endpoint does not move
(L12); three sprints of record exist only on this disk (L13); the test suite is green (L14);
the readout is insensitive to a distributional difference of nearly half the mass (L15); both
physics energies are measurably worse than noise (L16); the Hamiltonian barely changes between
targets (L17); the functional lever is closed in all five forms (L18).

Closed. The basis question; the physics lever; the prior-calibration lever; the quantum
contribution as a claim.

Retracted. L2's mechanism and magnitudes (L7); L5 and L9 as stated; the S25 audit retracted
five claims in all (L5, L7, L9, L11, L15).

Open. Whether a better distance predictor is obtainable; the 2/60 benchmark self-copy leak;
publishing the trainability half; the dynamical Lie algebra (measured in S26).

### VI.22 Sprint 26: examination, proposals and tournament (`s26/LEDGER.md`; `s26/BRIEF.md`)

Question. With the record examined end to end, do Proposals A (the quantum selector), B (the
scaled generation lane) and C (the classical ladder) or any tournament entry move the built
chain below 3.2148 A, and if not, what precisely closes?

Falsifier. `s26/BRIEF.md` states one-sentence falsifiers for each; the pre-registrations are
`s26/PREREG_*.md`; the phase gate held every endpoint experiment until Phase 0 was signed off.

Result so far. Phase 0: the examination reproduces every stored structure on all four bases at
0.0 and the two traced targets end to end (L28; Part III); the module map, pinned hashes,
not-in-git census and claim ledger are `s26/results/module_map.json`, `pinned_hashes.json`,
`claim_search.json` (`s26/EXAMINATION.md`); the governed test run is 370 tests, 357 passed, 13
real skips, 0 failed (lane I, L6 to L10, L15, L16, L21; `s26/TEST_RUN.md`); the identity flag
ships dark and the self-copy audit is `s26/results/i_identity_audit.json` (L15). Lane Q: the
deployed ansatz's dynamical Lie algebra is the full so(128) from depth 2, falsifying the
pre-registered proper-subalgebra prediction (L27; Part V.9). Lane PH: the steric singularity is
mostly the side-chain builder's (L23), and the production relaxation's native-free census (L24;
Part IV.7). Lane P: B1 is infeasible at its stated scale (L13); the repo-native stand-in (L14);
C1 reproduced; the C2 ladder training (L26). One Rule-1 incident, logged with no consequence
(L28 item 5, L29). The rest fills Part VII as the verdicts land.

Closed, retracted, open. Slots; see Parts VII and VIII.

## PART VII. WHAT S26 TESTED AND FOUND

Slots, filled as verdicts land in `s26/LEDGER.md` (the phase gate was still closed when this
part was opened; nothing here is a result yet):

- VII.1 Proposal A (lane Q): A1 the null test of the deployed selector; A2 the dynamical Lie
  algebra (landed, `s26/LEDGER.md` L27, written into Part V.9); A3 a target-dependent
  Hamiltonian; A4 gradient variance, grown against fixed ansatz. [pending]
- VII.2 Proposal B (lane P): B1 infeasible at the stated scale (L13); B2 the feasible-scale
  rungs; B3 the repo-native stand-in (L14). [pending endpoint]
- VII.3 Proposal C (lanes P, PH): C1 the reproduction (`s26/PREREG_C1.md`); C2 the ladder
  (training running, L26); C3 the matched-random control of the relaxation (native-free half
  landed, L24; the accuracy half waits on the gate); C4 routers; C5 the common-mode prediction.
  [pending]
- VII.4 Tournament entries: the AMBER reject filter (census L23), the cis-peptide gap, the
  identity null (L15), physics branch selection, strain as difficulty, rotamer relief,
  coherence-penalised training, window ensembling. [pending]

<!-- PART VIII -->

<!-- PART IX -->

<!-- APPENDIX A -->

## APPENDIX B. EVERY NUMBER IN THIS REPORT, WITH ITS ARTEFACT

Rows are added as parts land. A number whose only source is a document is not in the report;
the four document-only numbers of `s26/EXAMINATION.md` section C (36.1/36.4 deg, +0.0004/+0.0030,
the |z_moment| triple 0.7529/0.8013/0.1127, 355/13) are mentioned only as absent, and the 0.524
sampled tail-only cosine inside C24 is replaced by the three instrument-specific values of
`s25/QUANTUM.md` section 4.2. "as stored" means the leaf is the number to the precision printed;
"as cited" means the number is quoted from the named document section, which names its own
artefact; "as asserted" means a passing test pins it.

| number | where used | artefact (file :: key, or file:line) | stored value |
|---|---|---|---|
| 3.2148 | I, II, III | `bench_results/baseline_tuning126.json :: science/rmsd_arm/mean` | 3.214765154210998 |
| 3.0483 | II, III | `bench_results/baseline_tuning126.json :: science/rmsd_avg/mean` | 3.048338093879532 |
| 3.2355 | II, III | `bench_results/baseline_tuning126.json :: science/rmsd_full/mean` | 3.2354598538973844 |
| 3.2041 | II | `bench_results/baseline_tuning126.json :: science/rmsd_fit/mean` | 3.2040761603809194 |
| 3.2126, 0.171 | II | `results/summary/leaderboard.json :: rows[0]/mean`; S26 ledger L14 | 3.2126212503925293 |
| 3.2187 ... 4.1015 | II | `results/summary/leaderboard.json :: rows[*]/mean` | as stored |
| 1.7108 | I, III | `bench_results/baseline_tuning126.json :: science/pool_best/mean` | 1.7108244199364904 |
| 1.3134 | III | `bench_results/recon_library_saturation.json :: universe_best/small` | 1.3134468768690186 |
| 2.3062 | III | `s26/results/e_reproduce.json :: summary/stored_top_m_best/mean` | 2.3061526409453816 |
| 3.4540 | III | `s26/results/e_reproduce.json :: summary/stored_shipped/mean` | 3.4540004952559396 |
| 2.6087 | III | `s17/results/inband.json :: cells/25` band_best mean (reproduced `s26/PREREG_C1.md`) | 2.6087 |
| 2.2261, 2.8334, -2.1496, 0.822 | III | `s24/results/priorladder.json :: rows[*]/MASS1.0, MASS0.1, MASS0.0` (means 2.226080, 2.833382, 3.048338) | derived |
| 4.0648 | I, II, III | `s12/results/s14_ladder.json :: rows/L0_constant_helix/mean` | 4.064753929494389 |
| 2.9610, 2.9507 | I | `s9/final_report.json :: dist/full/mean, dist/shipped/mean`, asserted within 5e-4 by `tests/test_pipeline.py::test_the_committed_benchmark_report_still_holds_its_reference_numbers` (not opened by this lane) | as asserted |
| +0.0103 [-0.1596, +0.1803], 31W/29L | I | `docs/FINDINGS.md:4479-4607` (S9-10), naming `s9/final_report.json` | as cited |
| 2.9614, 3.8122, 22.3%, 0.649, 80 of 126, 3.803955, 1.8e-15, 3.867 | II, III | `s26/results/e_reproduce.json :: summary/virtual_bond` | as stored |
| +0.1664, +0.0977, [+0.0049, +0.3379], 16 of 126 | II, III | `s26/results/e_reproduce.json :: summary/projection_gap_arm_minus_avg` | as stored |
| 0.0342, 0.0958, 0.2051 | II | `s26/results/q_mde_reference.json` | as stored |
| 0.394%, 1.18% | III | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` (hamiltonian stage); `s25/results/q_gibbs.json :: results/spectrum_target_independence` | as stored |
| 0.6758 (68%) | III | `s23/results/errdecomp.json :: rows[*]/f_common` mean | 0.675770 |
| 0.220, -560, 125/126, 5.38, 4.86 | III | `s26/results/ph_c3_nativefree.json` (S26 ledger L24) | as stored |
| 5370, -1290.6, 0.19, 52.0; 6.0e12, 1262.4, 0.61, 1166.1; 1172.7 (2BP4) | III | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` (relax stage); `bench_results/cache/1fc9f2dcf489e2fb/2BP4.json` | as stored |
| 7193, 9814, -1..19, -6..32, 144, 420, 105, 91, 183, 17 centres, 760, 488, 483, 0.6926..7.2688, 2.2420..5.1791, 401 98 193 161 284 260 372 373 381 382, 44/449, 3.75, 1.97, 1.16, 8, 128, 5.05..21.07, 0.23..2.70, 5.62..23.44, 7.58, 5.34, 1.40 | III | `s26/results/e_trace_1S9Z.json`, `s26/results/e_trace_9KAR.json` | as stored |
| 3.3135, 3.3443 | III | `s8/integrate_vqe.json :: rmsd_vqe_sel, medoid128` (reproduced by `tests/test_pipeline.py::test_the_four_components_all_execute_and_reproduce_their_published_numbers`) | as stored |
| 0.288 | III | `docs/FINDINGS.md:1970-2037` (S7 finding 11) | as cited |
| 0.016 | III | `docs/FINDINGS.md:2379` (S8-6 heading) | as cited |
| +0.142 | III | `ARCHITECTURE.md` section 2.5 | as cited |
| +0.0207 [+0.0143, +0.0276] | III | `docs/STATE_BRIEF_2026-09-12.md` section 4 | as cited |
| 1.386 | II | `docs/FINDINGS.md:2751` (S8-8) | as cited |
| 470, 13, 47 of 500 | II, III | `tests/test_data.py:305`; `README.md` (fold repin; layout) | as asserted |
| 23, 5 | II | `docs/FINDINGS.md:60-118`; `s25/LEDGER.md` L5, L7, L9, L11, L15 | as cited |
| 4/126, 2/60 | II | `s24/LEDGER.md` L4; `s26/results/i_identity_audit.json` | as cited |
| a40581ad...422d, 8002 | II | `s26/results/pinned_hashes.json` | as stored |
| 1000 | III | `core/amber.py:1282`; `s8/integrate.py:292` | as in source |
| 2.6e-4 | II | `s25/LEDGER.md` L10 | as cited |
| 64.7th percentile | II | `s25/results/q_alpha.json` (concentration check, `s25/QUANTUM.md` 6.4) | as cited |
| Legacy weights 4.0, 1.0, 1.0, 3.0, 2.0, 2.0, 0.5, 1.0, 0.8, 0.15, 0.4 | IV | `core/energy.py` `DEFAULT_WEIGHTS` | as in source |
| 0.3 ms, 9 ms, 6 to 12 s | IV | `s16/repair_FINDINGS.md` section 2.4 | as cited |
| 28 ms (corrects 6 ms) | IV | `docs/CONDENSED_REPORT.md:185` (corrections table) | as cited |
| 112/192, 186/192, 0.36 to 0.88 | IV | `s20/LEDGER.md` L6 | as cited |
| -0.758 (SE 0.023), +0.60, +1.103 (SE 0.042), -0.27, -0.090 (SE 0.018), +0.307, +0.279, -0.027, +0.047 | IV | `s25/results/phys_landscape.json :: summary` | as stored |
| -0.0886, -0.0829 | IV | `s20/results/c_q1.json`; S24 (via `s26/PH_PART_IV_NOTES.md` section 3) | as cited |
| 3.058, 3.132, 3.215, 3.253, 3.425, 3.674, 3.755, 3.881 (point cloud) | IV | `s25/results/phys_suite.json :: configs_rank, random_null_rank` | as stored |
| +0.330 (1.99x MDE), +0.455 (2.42x MDE), +0.327, +0.471 | IV | `s25/agentPHYS_FINDINGS.md` sections 1.4, 1.5 (from `s25/results/phys_suite.json`) | as cited |
| 3.2187, 3.3100, 3.3732, 3.4221, 3.8248, 3.8844, 4.1015 (built chain) | IV | `results/summary/leaderboard.json :: rows[*]/mean` | as stored |
| +0.068 [+0.012, +0.125] | IV | `s16/LEDGER.md` L17 (`s16/integrate.py`; interval flagged provisional there) | as cited |
| 58.6%, 15.3, 97.0%, 99.66% | IV | `s25/results/phys_landscape.json :: AMB_frac_absz_lt_0p1, AMB_decades, AMB_top10_var_share` | as stored |
| 99.56% (60.7% to 99.98%), 46.5%, 26.6%, 25, 41 | IV | `s13/results/walsh_xval.json :: concentration` | as stored |
| 1.000, 0.0, 3, 900, 7e30 to 2e31, 6 of 6 | IV | `s13/results/walsh_amber.json :: exact[*].terms`; `s13/walsh_FINDINGS.md` 3.1, 3.2 | as stored |
| 0.0745, 0.4221, 30W/0L, 18.8, 5.8, 0W/30L | IV | `s20/results/c_land_report.txt` section 1 | as stored |
| 0.577, 0.104, +0.444 of +0.620 (72%) | IV | `s20/results/c_land_null.json` | as stored |
| 9450, 40.7/75, 2.6/75, 2.66, 5057, 96.8%, 2627, 2267, 163, -0.74 | IV | `s26/results/ph_reject_census.json :: singularity` | as stored |
| 40/126, 462/500, +0.054 | IV | `s25/results/phys_suite.json :: normalisation_fork`; `s25/agentPHYS_FINDINGS.md` section 3 | as stored |
| 58.7%, 8.6e4, -560, 125/126, +1262, 0.220, 3.804 to 3.867, 5.38, 4.86 | IV | `s26/results/ph_c3_nativefree.json` | as stored |
| n = 7, 128, L = 3, P = 21, 21 RY, 21 CNOT, depth 21 / ~24, 50 steps, seed 0, 8128 | V | `core/pipeline.py:186-189`; `core/quantum.py:885-982`; dim so(128) by arithmetic | as in source |
| 5.6e-17, 3.331e-16, 6.661e-16, {1: 2, 2: 4, 3: 8, 4: 16}, [0.879145, 0.461538, 0.105625, 0.054141, 0, 0, 0, 0], 0.334 | V | `s25/results/q_verify.json` | as stored |
| 4.06e-2, 3.4371, 1.18%, 0 of 8 | V | `s25/results/q_gibbs.json :: results/spectrum_target_independence` | as stored |
| 0.0761 bits | V | `s25/results/q_alpha.json` (alpha = 1, T = 0.1 cell) | as cited (`s25/QUANTUM.md` 3.4) |
| VQE_LFO table, 78 of 126, 0.6190 | V | `core/pipeline.py:113-118`; `s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint` | 0.6190476190476191 |
| 1.000000000, 4.597e-10, 4.663e-10, 42 | V | `s25/results/q_verify.json` | as stored |
| +0.655634 / 0.758, +0.566586 / 0.519, +1.000000, +0.994 | V | `core/quantum.py:66` (S9 instrument); `s25/results/q_verify.json` (S25 re-verification) | as stored |
| 2592, 1620, 0, 0, 1424 (54.9%), 972 / 972; 3888, 29.9%; 17574, 58.8%; 2016 / 2016 | V | `s25/results/q_verify.json`; S24 harness and coordinator runs as cited in `s25/QUANTUM.md` section 5; `s22/LEDGER.md` (Gate 1) | as stored / as cited |
| 3.4540, 3.3414, 3.2835, 3.3135 | V | `s8/integrate_vqe.json` (S8 instrument rungs) | as stored |
| alpha-effect table: -0.1126 (0.0792), -0.1067 (0.0644), +0.0279 (0.0518), -0.0011 (0.0477), +0.0126 (0.0268), +0.0039 (0.0240) | V | `s25/results/q_alpha.json` | as stored |
| -0.7423, +0.2700, -0.0234, 0.032, 0.7045, +0.0090 (sd 0.0342), 0.0268 | V | `s25/results/q_alpha.json` (entropy curve) | as stored |
| +0.0499 (0.0427), -0.0302 (0.0450), +0.0445 (0.0511) | V | `s25/results/q_alpha.json` (circuit minus Boltzmann) | as stored |
| Gibbs ladder at T = 0.3 (seven rows), 1.449 / 0.902 / 0.373, 0.761 / 0.453 / 0.351, -1.483528, -2.453671, 0.891 / 0.783 / 0.783, -1.7176 / -1.7186 / 0.075, 0.458 | V | `s25/results/q_gibbs.json :: results/*, results/training_control` | as stored |
| -0.1405 (0.0732, 66W/48L), -0.0002, +0.0262 (0.0279), -0.0308 (0.0590), -0.0081, 12, 45.3%, 64.7th | V | `s25/results/q_alpha.json` | as stored |
| -0.0311 (0.27x MDE), -0.0021 | V | `s25/results/q_alpha.json` (forced-alpha counterfactual) | as stored |
| -0.013 [-0.095, +0.077] | V | standing project result cited in `s25/QUANTUM.md` 7.4 (no S25 artefact; not re-derived) | as cited |
| width-sweep table (35 variances), slopes -0.6492 / -0.2522 / -0.0472 / -0.3105 / -0.2429, depth sweep (7 values), CVaR ratio tables (14 values), draws 250 / 200 / 120 / 80 | V | `s25/results/q_plateau.json` | as stored |
| DLA table: 120, 496, 510, 1023, 2016, 8128, 32640, 32766, 65535, 130816, 523776, 2096128; pools 36, 136, 528, 2080, 8256, 32896; ADAPT 7, 16, 1025 (12.6%); 12 of 12; 85.3 s; 0.479 GB | V | `s26/results/q_dla.json :: results/fixed, results/pools, results/adapt_sets, results/numeric`; `s26/jobs_done/a2_dla.json` | as stored |
| 0.504 (and 0.823, 0.768, 0.648, 0.537) | V | `s13/results/geo_kernel.json` (via `s13/SPRINT13_DOSSIER.md` section 10) | as cited |
| 1.0000, 5,000+, 74, 52, 1.4e-13, 4.6e-2 | V | `s13/qarch_FINDINGS.md` section 1 (`s13/SPRINT13_DOSSIER.md` section 6) | as cited |
| 141, 99.6%, 6.001 / 6.001, 0.0003, 0.878, 0.986 | V | `s13/walsh_FINDINGS.md` (`s13/SPRINT13_DOSSIER.md` section 7; raw sweep `geo_pauli_v1_rawonly.json`) | as cited |
| 2.236, 3.015, 79 / 79, 0.641, 0.392, 0.069, 0.122, 12 / 13, 0.778, 0.629, 1.006, 1.001, 1.278, 0.913, 95, 0.15 to 4.18 | V | `s13/SPRINT13_DOSSIER.md` section 8 (Pauli tables of `s13/results/`) | as cited |
| 1.000 on 26 cells, 0.955, 3.89 / 4.10 / 4.17 / 4.27, 0.60 to 0.64, 2.70 / 10.00 / 1.27, 7.9, 1e-16 | V | `s13/SPRINT13_DOSSIER.md` section 9 (`s13/results/walsh_amber.json`) | as cited |
| 1.00, 2^(-0.47n) to 2^(-0.86n), 0.93 to 1.00, 0.862 to 0.078 | V | `s13/SPRINT13_DOSSIER.md` section 10 (`s13/results/geo_kernel.json`) | as cited |
| 0.000e+00, 0.2500, 0.008 to 0.037, 1.7e16, 0.088, +0.333 | V | `s13/SPRINT13_DOSSIER.md` section 11 | as cited |
| +0.209 to +0.301, +0.27 to +0.62 | V | `s16/LEDGER.md` L17, L29 (`s16/integrate.py`; `s16/qphase_FINDINGS.md` section 3) | as cited |
| 20 of 20 | V | `s12/SPRINT12_DOSSIER.md` section XVIII | as cited |
| 2.406, 1.925 | VI | `docs/FINDINGS.md:2287` (S8-5); `docs/CONDENSED_REPORT.md:176` (correction) | as cited |
| 3.204 | VI | `docs/FINDINGS.md:3603` (S8-11) | as cited |
| 0.5 per pair, 0.044 | VI | `docs/FINDINGS.md:3810` (S8-14) | as cited |
| 2.7x | VI | `docs/FINDINGS.md:4726` (S10-2) | as cited |
| 1.62 | VI | `docs/FINDINGS.md:5219` (S11-2) | as cited |
| 3.184, 3.203, -0.019 [-0.058, +0.020], 63W/63L, 3.043 to 3.026, 2.534, 1.004, 5.425, 6.019, 10/18, 6/108, 1.2e-6, 1.486 | VI | `s12/SPRINT12_DOSSIER.md` sections I, IV, XVII | as cited |
| 1.16, 0.04, 0.893, 882, -1.743, -0.029 | VI | `s12/SPRINT12_DOSSIER.md` section III; `docs/CONDENSED_REPORT.md:61-65` | as cited |
| 1.594, 88%, 0.388, +0.043, -0.088, 32nd to 40th, +0.139, 3.764 / 3.667 / 3.920, 10.4, 3.770, 262,144, 9 | VI | `s13/SPRINT13_DOSSIER.md` sections 1 to 5; `docs/CONDENSED_REPORT.md:84-116` | as cited |
| -10.9% | VI | `s13/SPRINT13_DOSSIER.md` section 9 | as cited |
| 54/126, 2.021, 55, -2.253 [-2.642, -1.865], 0.885, 0.139, 0.68, +0.110 [+0.004, +0.214], 0.29, 3.4x, 0.171, +1.729, 98.8%, +0.014 [-0.001, +0.028], 0.511, 0.25, 16 to 19% | VI | `s14/LEDGER.md` tables DEMONSTRATED, REFUTED, CLOSED (C8, C9, C11, C12, C13, C14, E1, E2, E3, V2) | as cited |
| 787/787, 1.7108244199364904, 2.04e-13, 63,000, 8.3 to 23.3, 16/126, 0.611, 3.644, +0.440 [+0.290, +0.592], 2.386, +0.509, +1.492, 2.633, 34.7th, 4/126, -0.949 [-1.147, -0.753], +0.117 [+0.050, +0.188], 67, 178, 141, 21, 16, 13, 3.321 | VI | `s15/LEDGER.md` Phase 0 rows 0.1 to 0.15, rows 1.3, 1.4, 2.1, 2.2, 2.6, workstreams 3.2, 3.9 to 3.11; `docs/STATE_BRIEF_2026-09-12.md` section 7 | as cited |
| 0.162 | VI | `s16/LEDGER.md` L16 | as cited |
| 1.986, 1.33, 1.338, 0.015, 60% | VI | `s21/LEDGER.md` L18, L21, L23, L30, L11 | as cited |
| 0.482 | VI | `s22/LEDGER.md` L1 | as cited |
| +0.166, 7.2%, -2.15 | VI | `s24/LEDGER.md` L9-A, L14, L13 | as cited |
| +0.164 | VI | `s20/LEDGER.md` L2 | as cited |
| 370, 357, 13, 0 | VI, IX | `s26/TEST_RUN.md`; `s26/results/test_run.json` | as stored |
<!-- APPENDIX B ROWS -->

<!-- APPENDIX C -->
