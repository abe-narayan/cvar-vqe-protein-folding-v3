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

<!-- PART V -->

<!-- PART VI -->

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
<!-- APPENDIX B ROWS -->

<!-- APPENDIX C -->
