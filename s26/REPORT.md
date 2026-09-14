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

Status: final, at the close of Sprint 26 (ledger L142, 2026-09-14 04:14). Parts I to VI, IX and
Appendix A are written from the record; Parts VII and VIII carry every S26 verdict with its
ledger entry and artefact; Appendix B lists every number in the report, and
`python s26/e_report_check.py` (or `python s26/examine.py --report-check`) re-reads each one from
its artefact; Appendix D reconciles this report with `docs/REPORT_S26.md`; Appendix C reproduces
`s26/LEDGER.md` verbatim and is the last thing in the file.

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
phi with a mean absolute error of 36.133 degrees against 36.416 for a model that reads no
sequence at all (`s26/results/a_c26_phi_mae.json`, re-derived by the S26 Adversary from
`s13/cache/tors_rows.npz`, ledger L31, L32; the record's 36.1 against 36.4 is
`docs/CONDENSED_REPORT.md:104-116`, `s13/SPRINT13_DOSSIER.md` section 5). Short peptides in
isolation are also often flexible, and the deposited model is one member of an ensemble. Both
facts matter for what "accuracy" can mean.

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
[-0.1596, +0.1803], 31 wins to 29 losses (`docs/FINDINGS.md:4503-4631`, S9-10; the two means are
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
retraction count that taught them is 23 corrections in `docs/FINDINGS.md:60-142` plus five S25
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
  (`docs/FINDINGS.md:4503-4631`); the CLI refuses to run it without
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
  (`s24/LEDGER.md` L4; `s26/results/i_identity_audit.json`, S26 ledger L15). At the looser
  >= 0.6 identity criterion 13 of 60 benchmark targets have such a pool window
  (`docs/FINDINGS.md:4595`, S9-10), priced at +0.0030 A in S10-4. The dev price, re-derived in S26
  with an artefact: clean minus production +0.0004 A on the lam = 0 chain (fold CI [-0.0001,
  +0.0010], MDE 0.0012) and +0.0018 A on the built chain [-0.0003, +0.0039], neither clearing
  its MDE (`s26/results/w_selfcopy_endpoint.json`; S26 ledger L44, L50). The benchmark effect is
  never measured, because measuring it would open the benchmark; it is bounded from the dev
  proxy without opening it: dev-proxy price 0.008 A (all four dev self-copies measured on both
  channels, L108); own-native envelope 0.028 A (the mean's
  fold-CI limit) to 0.151 A (the worst single target) on the built chain, under the envelope's
  assumption A2; MINOR under every reading on the built chain and on the paired gain; on the
  selection basis the worst-target reading is 0.194 A; it cannot move the benchmark verdict
  either way (`s26/results/w_selfcopy_bound.json`; L44, L55, L58, L108; Part VII.4).

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
`science/rmsd_fit/mean`) and is not a production arm. The spread behind the built-chain mean, on
the results-lab rebuild (`results/summary/leaderboard.json :: rows[0]`): median 2.9661 A, best
0.1824 (1S9Z), worst 8.2342 (2MQ2), 28.6% of targets under 2 A and 50.8% under 3 A, correlation
with the ORACLE pool best 0.7150, and a rebuild projection gap of 0.1643 A against the cache's
0.1664 (Appendix D).

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
  reading that order (`docs/FINDINGS.md:2775`, "A methodological correction that changed the
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
within PDB quantisation: worst absolute error 2.239e-4 A against a quantisation bound of
8.66e-4 A (`s25/LEDGER.md` L10).

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
(`docs/FINDINGS.md:1994-2061`, S7 finding 11; `docs/STATE_BRIEF_2026-09-12.md` 5.7).

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
(standardised error sd 1.9962 where a calibrated model gives 1; 90% intervals cover 61.6% and 50%
intervals 28.2%; 24.1% of pairs multimodal; `s25/results/calib.json`, means over 126 rows;
`s25/LEDGER.md` L1); its L1 risk's per-pair minimiser is the posterior median, which sits on one
of 17 atoms, so the per-pair target is quantised, up to 1.75 A of location error at long
separation from binning alone (`s25/LEDGER.md` L6); the fold models share most of their training
data (Part II). The ORACLE ceiling above it: with the prior replaced by
the native's own distances, everything else fixed, the point cloud reaches 2.2261 A point cloud
against 3.0483 point cloud (`s24/results/priorladder.json`, mean of column MASS1.0 over 126
rows against MASS0.0), a gain of 0.822 A; the slope at the origin is -2.1496 A per unit of
interpolation toward truth (first rung MASS0.1, mean 2.8334 point cloud, -0.2150 A at 2.4x its
MDE, 113W/13L; `s24/LEDGER.md` L13), but only along the native's own direction: across 51
achievable arms that move the score's target up to 25% toward truth, the correlation between
the distance moved and the endpoint change is +0.054, and the endpoint does not move
(`s25/LEDGER.md` L12).

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
purchasing perfect retrieval with the native, was worth 0.016 A (`docs/FINDINGS.md:2403-2635`,
S8-6), and structure and sequence are decoupled at this length (`docs/FINDINGS.md:501`, S5
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
worse than a coin flip (`docs/FINDINGS.md:2763-3011`, S8-8).

### III.5 The selector (off in production)

With `quantum=True` the filter keeps the top 128 (`top`, 128 indices, 128 x 128 matrix) and
`quantum_stage` (`core/pipeline.py:806`) forms `E = _zrank(sc[top])`: the 128 scores are
replaced by their ranks, standardised to zero mean and unit sd, giving a ladder from -1.7186 to
+1.7186 that is the same on every target up to tie-averaging (max deviation 0.394% of range on
both traced targets; worst 1.18% over the S25 sample,
`s25/results/q_gibbs.json :: results/spectrum_target_independence`). `H = diag(E)` is the
Hamiltonian of a 7-qubit register whose basis states index the 128 candidates; a 3-layer,
21-parameter circuit is trained for 50 Adam steps on `F = CVaR_alpha(E; p) - T H(p)` with
(alpha, T) from `VQE_LFO[fold]` (alpha 0.25 on folds 1 and 2, 1.0 on folds 0, 3, 4; T = 0.3
everywhere; `core/pipeline.py:113`), and the CVaR tail of the trained distribution is read
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
cloud on 16 targets. The representation itself costs almost nothing: the native CA trace
projected through this same operator lands 0.083 A from itself (0.043 A with the prior off; max
0.44; S26 ledger L89, ORACLE DIAGNOSTIC), so the cost of the step is the operator's displacement
of a non-native cloud, not the ideal geometry or the constant omega (Part VII.4). Recovering the contraction by other means does not help: the Frechet mean
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
native. On the means the step costs +0.0207 A, relaxed chain 3.2355 against built chain 3.2148
(`bench_results/baseline_tuning126.json :: science/rmsd_full/mean, science/rmsd_arm/mean`);
paired, +0.0207 with fold CI [+0.0154, +0.0290], 2.16x MDE, 40W/86L (S26 ledger L39). S16 found
a random displacement of matched size at least as accurate (`s16/LEDGER.md` L27), and S26's C3
stage 1 measured it: the relaxed chain is worse than the built chain displaced at random by the
same 0.220 A (+0.0111, fold CI [+0.0062, +0.0171], Type-M zone) and worse than the built chain
moved the same distance toward a random pool member (+0.0385 [+0.0298, +0.0479]); the AMBER
displacement points away from the native (ORACLE cosine -0.049). The step is dropped as an
accuracy step and kept as a validity step on 124 of 126 targets (Part VII.3; L39, L46). What it buys is validity: exact bond
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
and was corrected to 28 ms (`docs/CONDENSED_REPORT.md:185`); the S15 audit with memoisation
defeated then measured 8.3 to 23.3 ms and recorded that budget matching had over-charged AMBER
two to three times (`s15/LEDGER.md` row 0.8).

### IV.2 What "H = E after 50 relaxation steps" means

When AMBER is used as a Hamiltonian for a quantum selector, the structure it is handed is an
ideal-geometry chain that has never been relaxed, and its energy there is not a finite number on
42% of a lattice register (`s20/LEDGER.md` L6: finite at Relax_1 on 112 of 192 states, at
Relax_50 on 186 of 192). The deployable object is therefore "run 50 steps of restrained
minimisation, then read the energy". That operator, not the raw force field, is the
Hamiltonian; the relaxation is what makes the objective defined, and the Spearman correlation
between the energy after 50 steps and after 1 step is 0.358, 0.827 and 0.882 on three targets,
so the ordering is still moving when the cap stops it.

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
random" contrast is stated on the point cloud only). Against the production point cloud itself
(`s25/results/phys_suite.json :: vs_incumbent`, point cloud on both sides): Distogram +0.0097
(0.43x MDE), AMBER+Distogram +0.0833 (0.70x), Legacy+Distogram +0.1668 (0.91x), all three
+0.2047 (1.07x, Type-M zone), Legacy+AMBER +0.6259 (1.89x), Legacy +0.7070 (2.14x), AMBER
+0.8322 (2.67x); on the built chain the distogram-selected configuration is +0.0061 A over
production, 0.19x MDE, 61W/65L (`results/summary/leaderboard.json :: rows[1]`). Permuting a
physics channel while keeping its distribution improves the endpoint (Legacy +0.327 point cloud worse than its own noise,
AMBER +0.471; `s25/agentPHYS_FINDINGS.md` section 1.5). The reason is IV.3: the pool's dominant
axis is compactness, the distogram already selects on it, and each energy pushes along that
axis in a direction unrelated to which candidate is right. S16 had found the same for Legacy as
a ranker inside the architecture (+0.068 [+0.012, +0.125] worse than a matched random drop at
m = 5, same basis on both sides, the S16 integrate readout; `s16/LEDGER.md` L17), S17 had found the decisive AMBER-against-Legacy comparison at full
scale and concluded "the objective is not weak, it is wrong" (`s17/LEDGER.md` L29, L30), and S18
found every physics filter losing to a random gate (`s18/LEDGER.md` L17).

### IV.5 The steric singularity, and where it lives

An unrelaxed ideal-geometry rebuild of a retrieved window puts atoms on top of each other, and
the Lennard-Jones r^-12 wall turns one overlap into an energy of 1e4 to 1e29 kcal/mol (a relaxed
peptide sits at -1170 to -500). Measured: 58.6% of every pool is above 1e4 kcal/mol; the
energies span 15.3 decades; 97.0% of a pool lands inside |z| < 0.1 of a moment z-score; the ten
worst candidates carry 99.66% of the variance
(the fraction above 1e4 from `s26/results/ph_reject_census.json ::
summary/per_threshold/1e4/pool_frac_over`; the rest from `s25/results/phys_landscape.json ::
summary`).

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
chain is built, but only 2.6 of 75 do on the backbone plus CB (S19 measured 2.66 per 75,
`s19/LEDGER.md` L12). Of the 5,057
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
(`s25/agentPHYS_FINDINGS.md:57-58`, `:357`, `:379`, its section 3, computed from
`s25/results/phys_suite.json :: normalisation_fork`).
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
after the gate; Part VII.3): the relaxed chain is +0.0207 A worse than the built chain (fold CI
[+0.0154, +0.0290]), +0.0111 worse than the built chain displaced at random by the same
magnitude and +0.0385 worse than a move of the same size toward a random pool member (S26
ledger L39, L46); S16 had found the random displacement at least as accurate (`s16/LEDGER.md`
L27). S17's "steric
validity is free" was quoted against the wrong baseline and corrected by the workstream that
made it (`s17/LEDGER.md` L27), and S17 explained S16's cis-peptide defect in a way that
exonerates the force field: the ideal-geometry builder cannot represent the cis bond, and the
minimiser is asked to repair what the builder could not draw (`s17/LEDGER.md` L28).

### IV.8 What the physics is for

The record's settled role for both energies is a validity stage, not an accuracy stage:
`docs/FINDINGS.md:3687` (S8-12, "AMBER is a validity stage, not an accuracy stage");
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
(`core/pipeline.py:181-184`, `Config`). RY and CNOT are real matrices and the initial state is real, so
the amplitudes are real exactly and the reachable manifold lies in SO(2^n), not SU(2^n);
dim so(128) = 8128 against 21 parameters. The simulation carries all 2^n amplitudes; the CNOT
chain is composed into one basis permutation, and `probs_batch` simulates all 2P shifted
parameter vectors in one pass. Verified in S25 against a dense Kronecker-product simulator:
max |p - p_dense| = 5.551e-17 (`s25/results/q_verify.json`).

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
(alpha, T) pair comes from a leave-fold-out table (`core/pipeline.py:113`):
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
| S5, the first measurement (`docs/FINDINGS.md:392-427`) | -0.023 | not recorded |
| S9, 10 qubits, 1024 amplitudes enumerated (`core/quantum.py:42`) | +0.655634 | 0.758 |
| the consolidation audit, exact / sampled (`verify/cvar_audit.json`) | +0.6847 / +0.685 | not recorded |
| S25 re-verification, n = 7, L = 3, alpha = 0.15, deployed energy shape | +0.566586 | 0.519 |
| any instrument, constant baseline | +1.000000 | 1.000 |
| sampled estimator, constant baseline | +0.994 | about 1 |

A third figure, +0.524 over 36 checks on the S8 instrument, lives only in project memory and is
not in Appendix B (`s26/EXAMINATION.md` C24). There is no universal constant; what reproduces is
the sign, the order of magnitude and the mechanism. The defect was first found and priced in
S5 (`docs/FINDINGS.md:392`, "A defect in the shipped CVaR gradient") and re-audited against two
references in S8-9 (`docs/FINDINGS.md:3171`) and S9-5 (`docs/FINDINGS.md:4221`).

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
sd 0.0268 A, single window on both sides: the circuit sits on a curve fitted without it.

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
variation 0.761, 0.453, 0.351; the pre-registered falsifier (KL below 0.1 nats, "the circuit
attains its own optimum") fired against the S25 draft (`s25/results/q_gibbs.json ::
results/falsifier`). The mandatory control, best of 200 draws from the untrained
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
nothing else (the S20 continuous encoding, same CVaR estimator, same seed) is -0.013 A [-0.095,
+0.077], an interval containing zero and effects in both directions (`s20/LEDGER.md` L1); and
chi orders nothing (`s21/LEDGER.md` L34). Whether the entanglement changes the
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
arms (percentile 0.088) while returning the worst structure (+0.333 A against random, both
sides built chains from their bitstrings): better optimisation of a misaligned objective
produces worse physics.

### V.11 The quantum record across the sprints

- S5: the CVaR gradient defect found and priced (`docs/FINDINGS.md:392`); a VQE over segment
  choices with its encoding and limits stated (`docs/FINDINGS.md:428`).
- S6: VQE/CVaR assembly gives a narrow pool whose narrowness cannot be exploited
  (`docs/FINDINGS.md:1045`).
- S8-9: CVaR-VQE over the discrete hypothesis set, and the integrated four-component system
  priced component by component (`docs/FINDINGS.md:3264`, `:3299`); the gradient audited against
  two references (`:3147`).
- S9-5: refinement's failure located in selection; the four-component ablation conclusive; the
  gradient defect fixed and priced against an exact reference (`docs/FINDINGS.md:4178-4237`).
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

### V.12 S26 additions as they land

**A2 checked (L45).** The Adversary's independent re-derivation `s26/results/a_dla_check.json`
(no shared closure code, incremental Gram-Schmidt rank) agrees with `q_dla.json` at every cell
at n = 4 to 7, L = 1 to 4, for all four conjugation and gate-order conventions; the 8128 / 1025
reconciliation (the fixed ansatz's algebra against the strings ADAPT selected, 1025 / 8128 =
0.126) is correct. The "S25 slopes are a statement about depth 3" reading is an inference from
the 2-design literature, not a measurement here. Per growth step on the real targets (L138): the
algebra of a grown set tracks the appended strings, not what the circuit does; the inert sets of
the alpha = 1 targets reach closures of up to 530 dimensions and no grown set reaches so(128).

**A4, grown against fixed circuits (L35, checked L47).** Part VII.1 carries the table. In one
sentence for this part: on the deployed Hamiltonian an operator-growing (ADAPT) construction
produces product circuits at alpha = 1, whose gradient variance does not decay with width because
there is nothing entangled to train, and at alpha = 0.25 a 2-local family that decays like the
fixed ansatz (-0.302 against -0.243 log2 per qubit, a difference of -0.056 with a bootstrap 95% CI
of [-0.182, +0.087], L119); no width-scaling argument for the selector comes out of it
(`s26/results/q_var.json :: results/slopes`; `s26/results/q_var_boot.json`).

**A1, the ADAPT endpoint (L68, L70, L75).** Part VII.1 carries the numbers. For this part: on the
deployed Hamiltonian the optimum at alpha = 1 is a product state on every real target
(KL(Gibbs || product of marginals) at most 7.9e-4 nats), a 7-parameter RY layer reaches it, the
21-parameter fixed circuit stops 0.90 nats short, an adaptive ansatz free to entangle appends
operators worth less than 1e-3 nats and leaves a product state, and reaching the optimum exactly
moves the built chain by -0.014 to -0.022 A at seed 0 and -0.045 to -0.052 at seed 1, a third to
four fifths of what the comparison resolves (MDE 0.059 to 0.066): not measured on either seed,
the difference between the seeds being the fixed comparator's own initialisation (3.228 against
3.261 A, L139, L140). Proposal A's verdict is REPLACE (L69), on seed-independent facts. **A3, a target-dependent Hamiltonian (L125).** A
raw-score Hamiltonian at the deployed entropy makes the trained states target-dependent (124 of
126 distinct, a median 1.34 nats apart against 0.010 under the rank ladder) and the emitted
structure does not change (+0.0034 A, 0.04x MDE); without the entropy match the sharper states
are worse by +0.09 to +0.11 A at 0.7 to 0.8x the MDE: the readout responds to the entropy of the
weights and to nothing else that was varied.

## PART VI. THE TWENTY-TWO SPRINTS AS A STORY

Sprints 1 to 4 predate the record kept here; their tree is in git history at `5fa05cd` and their
one surviving claim (that physics ranks real geometry) was re-measured and overturned in S5 and
S7. From S5 the record is `docs/FINDINGS.md` (S5 to S11, with its corrections ledger at lines
60 to 142 (the S26 block appended at L130) and its index at 143 to 235), the sprint dossiers `s12/` to `s15/`, and the live
ledgers `s14/LEDGER.md` to `s26/LEDGER.md`. `docs/CONDENSED_REPORT.md` condenses S5 to S13 and
`docs/STATE_BRIEF_2026-09-12.md` section 7 gives one line per sprint. Each section below has the
same six fields. "Falsifier" names the pre-registered condition where one existed; formal
pre-registration began in S15, and before that the field names the test that decided the
sprint. Numbers are quoted from the named section or ledger entry, which names its artefact;
Appendix B lists them as cited.

### VI.1 Sprint 5: the wall is located (`docs/FINDINGS.md:147`, sections at 248 to 780)

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

### VI.2 Sprint 6: selection is measured out (`docs/FINDINGS.md:164`, sections at 781 to 1331)

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

### VI.3 Sprint 7: the predictor is the constraint (`docs/FINDINGS.md:179`, sections at 1332 to 2088)

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

### VI.4 Sprint 8: the full architectural attack (`docs/FINDINGS.md:193`, sections S8-1 to S8-14 at 2089 to 3899)

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
filter plus consensus medoid at -0.172 A [-0.316, -0.027] single window, 74W/47L (S8-8,
`docs/FINDINGS.md:2917-2959`); the channels' errors are
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

### VI.5 Sprint 9: the benchmark pass (`docs/FINDINGS.md:210`, sections S9-1 to S9-10 at 3900 to 4607)

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

### VI.6 Sprint 10: adjudication and bounds (`docs/FINDINGS.md:222`, sections S10-1 to S10-5 at 4608 to 5159)

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

### VI.7 Sprint 11: performance engineering and consolidation (`docs/FINDINGS.md:229`, S11-1 at 5160, S11-2 at 5219)

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

Result. Every S12 contrast below is on the same basis on both sides, the S12 dossier's emitted
structure. The learned set decoder emits 3.184 against production 3.203 (d = -0.019 [-0.058,
+0.020], 63W/63L), with a flat learning curve 3.043 to 3.026 from n = 8 to 75 while the leaked
label reaches 2.534 at n = 8: signal-limited, not sample-limited (section I). The terminal
operator consumes the set mean, not the set best: d_out = 1.16 d_set_mean + 0.04 d_set_best,
R^2 0.893 over 882 perturbations, so a perfect rank-1 decision is worth -1.743 A through argmin
and -0.029 A through the m = 75 average (section III; `docs/CONDENSED_REPORT.md:61`). Sequence
conditioning is worth 0.776 A over all 126 (blind pipeline 3.989 against shipped 3.213) and
1.004 A on the 108 ordinary targets, and is negative on the 18 hardest, where a pipeline given
random windows emits 5.425 A against the shipped 6.019 A; three routers are null (section I). Ten of the eighteen failures are fibril segments and lasso peptides, 10/18
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

Result. Every S13 number below is a built chain from torsions, the same basis on both sides of
each contrast. The representation is not the barrier: at about 24 live qubits the space contains
a 1.594 A answer, but 88% of what the library buys is generic Ramachandran and 0.388 A of the
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

Result. Every S14 contrast below is on the same basis on both sides, the S14 ladder's emitted
structure (built chains from torsions), ORACLE arms labelled. The chemical-shift route is closed
by arithmetic: 54 of 126 targets are runnable and
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
`pool_best`. The generative architecture (built chains from torsions against the incumbent's built chain,
the same basis on both sides): torsion distance geometry reaches an ORACLE 0.611 A and a
predicted 3.644 A, +0.440 [+0.290, +0.592] against the incumbent, a loss (1.3); the
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
(L11) and the exhaustive argmin does not beat a zero-evaluation pool, the same basis on both sides through
one readout operator (L14, complete at n = 126, L17); the answer is in the objective's top 5% and a top-512 readout ceiling is 1.986 A (L18);
the source is not the lever, the retrieval pool beats the latent through the same operator
(L20); all remaining headroom is in-pool selection, worth 1.33 A (L21), and the 1.338 A
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
C1 reproduced; the C2 ladder training (L26). Lane W: the 2/60 proxy bound pre-registered and a
native-free census of the four dev self-windows (L30). One Rule-1 incident, logged with no
consequence (L28 item 5, L29). The Adversary's audit of the examination: reproduction exact on
HEAD, one material documentary item (the test total, corrected by addenda to 370 / 357 / 13),
C26 re-derived with an artefact, C27 located in git history (L31, L32, L34); Phase 0 signed off
at 09:05 (L33). After the gate: A4 (L35), the cis floor (L38), C3 stage 1 (L39), the steric
reject on the point cloud (L43), the 2/60 leak bound (L44) and the Adversary's checks of each
(L45 to L50), all in Part VII. The campaign paused 09:30 to 19:14 (L37, L40 to L42).

Closed, retracted, open. Slots; see Parts VII and VIII.

## PART VII. WHAT S26 TESTED AND FOUND

Filled from `s26/LEDGER.md` from L33 (PHASE 0 SIGNED OFF, 09:05) to L142 (SPRINT 26 CLOSE,
04:14 the next morning); each item names its ledger entry, its artefact, its pre-registration
and the Adversary's check. Nothing was written here before its ledger entry existed. Bases are
named on both sides of every contrast. What did not run is recorded as such with its reason.

### VII.0 How S26 was run

The sprint was a campaign of eight agent lanes under one coordinator, run on the 16.75 GB,
8-core box the record was made on and controlled by four things.

**The lanes.** E (examiner and librarian, then report writer: `s26/EXAMINATION.md`,
`s26/BRIEF.md`, this report), I (infrastructure: the governed test run, `s26/examine.py`, the
identity flag, the results-lab rebuild), Q (Proposal A: A1 to A4), P (Proposals B and C: the
feasibility of B1, the stand-in B3, the prior ladder C2, the routers C4 and the common-mode
prediction C5), PH (physics: the cis and steric censuses, C3, the reject filter, strain, branch
selection), W (wildcard: the 2/60 bound, the identity floor, the tie-break floor, ensembling),
A (the Adversary: the Phase 0 audit, a check of every ledger entry within the hour,
`s26/RETRACTIONS.md`, `s26/DELIVERABLES_CHECK.md`, the tournament ranking) and PR (the deck,
`vqe_research_overview.pptx`, every number read from an artefact and registered in
`s26/pr_values.json`). Each lane has a brief in `s26/briefs/`, a findings file
`s26/agent<lane>_FINDINGS.md` that names what it did not do, and two STATUS lines per hour in
`s26/STATUS.md`. The coordinator rules, spawns, and writes the ledger entries that decide.

**The ledger.** `s26/LEDGER.md` is append-only and numbered; five lanes append concurrently, so
the rule is re-read the tail immediately before appending and suffix on a collision; nothing
earlier is edited, corrections are appended (L57 corrects L56, L75 corrects L68). Every
endpoint number quotes its artefact, its job record (`s26/jobs_done/<name>.json`: exit code,
wall, peak RSS) and its pre-registration; the Adversary's check of each entry is its own entry
(L31, L34, L45 to L49, L54, L55, L70, L71).

**The phase gate.** No endpoint experiment (nothing that reads a native to score a native-free
operator) ran until the examination of Phase 0 had been audited: the examination landed as L28
(reproduction exact on all four bases, the claim ledger, one Rule-1 incident logged in L29), the
Adversary's audit found one material documentary item and re-derived one document-only number
(L31), the coordinator applied the addenda (L32), and PHASE 0 SIGNED OFF was posted at 09:05
(L33, re-checked L34). Gated scripts refuse to run before that string exists in the ledger
(lane W's `endpoint` and `floor`, lane Q's `label`). Every experiment has a `s26/PREREG_*.md`
filed before compute with its falsifier, comparison arm, basis and MDE; a pre-registration is
never edited, only appended to; the tournament of ideas (`s26/TOURNAMENT.md`, L51) was ranked by
the Adversary before any idea ran.

**The governor.** `s26/governor.py` samples every 5 s and holds the band 88 / 90 / 93 / 95%
(launch queued work under 88% for 60 s; resume under 90%; suspend the newest job above 93% RAM
or CPU; kill and requeue the newest above 95% for 15 s; at most two AMBER-tagged jobs); every
job runs through `s26/jobrun.py`, which registers it, waits for headroom and records exit, wall
and peak RSS. It grew through the night: v2.1 re-checks the launch cap after a jitter after six
jobs launched against a cap of four (L36); v2.2 refuses a stale snapshot for any job above
0.5 GB and any launch without est-ram + 0.5 GB free (L41); the governor itself gained a
retrying atomic write after a Windows file-replace race killed it (L59) and a stall breaker
that kills the fattest suspended job when the box sits at or above 90% for 180 s with every job
suspended (L74); its first kill (L78) showed that the CTRL_BREAK grace signal does not cross
consoles on Windows, so a governor kill is a hard kill after 25 s and the loss is bounded by
each job's own checkpoint granularity (per target or per fold in every S26 job); v2.3 adopts
suspensions made before its start and reads the launch cap from a file, raised from four to
six (L83), after which seven waiters that had read the old cap at import were found starved
for two hours and relaunched (L92); lane I was starved for 75 minutes by other lanes' relaunch
cadence (jobrun has no priority) and the cap went to seven with non-critical launches held
(L111); and v2.4, the last change of the sprint, keeps the governor's own record of what it
suspended after lane I found that a job whose root blocks in a subprocess reads "running" to
psutil after a suspend, so v2.3 re-suspended one audit's tree 18 times and never resumed it,
and Windows counts suspends per thread so the child then needed 16 resume calls (L128, L129).
The box was shared with the user's own load all night, so memory, not CPU, set the parallelism
(headroom 1.4 GB after the pause, L41).

**The interruptions.** Four, none of which lost a number: the API session limit cut every lane
from 00:46 to 08:40 (L17); the user paused the campaign at 09:30 at a usage limit, the host
killed the governor for low memory at 10:10 and every lane was cut by the session limit until
the resume at 19:14 (L37, L40, L41, L42); the governor died on a file-replace race at 19:36 and
was restarted within a minute (L59); and at 21:50 the box stalled at 91.9% RAM with every job
suspended until the stall breaker and a hand stop of the `raw` rung's training cleared it
(L74). Each time the lanes resumed from their STATUS lines and the per-target checkpoints on
disk; no finished work was restarted. At 22:02 the user extended the sprint to a close at
about 04:30 Pacific on 2026-09-14, with this report due at about 04:45 (L77): every deferred
run, test and diagnostic goes ahead under the same rules.

### VII.1 Proposal A: the quantum selector (lane Q)

**A2, the dynamical Lie algebra (L27; Adversary L45: STANDS).** Written in full in Part V.9.
The pre-registered proper-subalgebra prediction was falsified: the deployed RY/CNOT chain-plus-ring
ansatz generates the full so(2^n) from depth 2 at n = 7 (dim 8128), so the algebra offers no
protection against a 2-design plateau and the S25 slopes are a statement about depth 3. The
Adversary re-derived every cell independently (`s26/results/a_dla_check.json`, no shared
closure code, a different numeric rank rule): 8128 at n = 7 for all four conjugation and gate-order
conventions; n = 4, 5 numeric and symbolic counts agree 6 of 6. The per-growth-step item the
prompt asked for under A2 landed as L138 (`s26/q_dla_a1.py` -> `s26/results/q_dla_a1.json`, the
exact Lie closure of every ADAPT run at every growth step over the 126 A1 records;
`s26/PREREG_A2.md` addendum 2; property measurement, no native, no score; the job's wrapper lost
the child before the summary write and the summary was recomputed in-process from the stored
ladders). Final closure dimension, median [min, max]: pool V, Adam, alpha = 1: 7 [7, 82] with dim
7 on 48 of 78 targets (exactly those without a multi-qubit string appended); pool L2, Adam,
alpha = 1: 11 [7, 139], dim 7 on 15; pool V, L-BFGS, alpha = 1: 58 [7, 530], dim 7 on exactly the
18 targets where nothing was appended; pool L2, L-BFGS, alpha = 1: 37 [7, 513], dim 7 on the 10
with nothing appended; at alpha = 0.25 pool L2 reaches a median 1025 [513, 2017] under Adam and
1025 [258, 2017] under L-BFGS, pool V 16; no grown set on any target reaches so(128) = 8128, and
pool V never exceeds 530 against its whole-pool 2080. All four predictions (P2a to P2d) held.
Reading: the dynamical Lie algebra is a property of the generator set; on the 78 alpha = 1
targets the grown sets whose appended strings are inert have closures of up to 530 dimensions
while the fixed ansatz at the same parameter count has 8128 and does no more with it, so for a
grown circuit the DLA dimension counts what could be done if the appended angles were not near
zero, not what the circuit does; the only sets with a large algebra and non-trivial angles are
the alpha = 0.25 L2 sets (1025 of 8128), and A1 measured their endpoint at +0.0000 A against the
fixed circuit on those 48 targets. This closes the A2 deliverable.

**A4, gradient variance of grown against fixed circuits (L35; Adversary L47: STANDS WITH
CAVEAT).** `s26/q_var.py` -> `s26/results/q_var.json` (complete; 2,765 s, peak RSS 0.08 GB),
pre-registered in `s26/PREREG_A4.md`; energies are the deployed shape (standardised ranks 1 to
2^n), no target, no native. Gate: the n = 7 row of all five S25 cells reproduces bit-identically
(`results/reproduction/passed`, `worst_rel` 0.0). Fitted log2 Var[dF/dtheta_0] per qubit over
n = 4 to 13 (`results/slopes`; ADAPT growth by Adam best iterate, 50 steps, pool V = the Tang
2n - 2 strings, pool L2 = all 2-local odd-Y strings):

| cell | fixed (S25) | grown, pool V | grown, pool L2 |
|---|---|---|---|
| alpha = 1, T = 0 | -0.649 | -0.079 | +0.006 |
| alpha = 0.25, T = 0 | -0.252 | degenerate | degenerate |
| alpha = 0.10, T = 0 | -0.047 | degenerate | degenerate |
| alpha = 1, T = 0.3 | -0.311 | +0.024 | -0.008 |
| alpha = 0.25, T = 0.3 | -0.243 | -0.239 | -0.302 |

(The ledger's table prints the P = 3n-only fits, +0.035 and -0.246 for the two pool-V cells; the
stored slopes above use every row.) What the grown circuits are: at alpha = 1 they hold 1 to 3
distinct operators, the other 12 to 25 selections being consecutive repeats of one single-qubit
Y, so they are product circuits with the parameter count of the fixed ansatz, which is why their
variance does not decay with n (a large gradient from a product circuit is the trivial regime,
not trainability). At alpha = 0.25, T = 0.3, pool L2 grows 4 to 21 distinct 2-local strings and
decays at -0.302 per qubit against the fixed ansatz's -0.243. At T = 0 with alpha < 1 growth
stops at P = n on every width (all first-order gradients vanish at a basis state), degenerate by
construction as pre-registered. Grown-to-fixed variance ratio at matched n: 1.63 at n = 4
(alpha = 1, T = 0.3) and 2.5 to 370 everywhere else, 31 of 32 matched rows above 2. Predictions:
H4a held on three of four cells (two sit 0.035 and 0.002 outside the registered band, inside the
fit's error); the falsifier (a slope below -0.5) did not fire; H4b failed on 1 of 32 rows; H4c
held on P < 3n (14 of 14) and failed on Var < 1e-3 for 3 of 7 rows; H4d held. Adversary caveats:
the product-circuit reading rests on the ADAPT closures in `s26/results/q_dla.json ::
results/adapt_sets`, not on `q_var.json`; and "-0.302 is the fixed ansatz's -0.243 within the
sampling error" is asserted without a persisted slope CI (the qualitative conclusion, both slopes
of order -0.25 to -0.30 and far from a 2-design's -1.0, is safe). The bootstrap the caveat asked
for landed as L119 (`s26/q_var_boot.py` -> `s26/results/q_var_boot.json`; 3,928 s, peak RSS
0.09 GB; `s26/PREREG_A4.md` addendum 2): every A4 row re-measured with the same draws reproduces
`q_var.json` at relative deviation 0.0 (35 fixed and 70 grown rows), and 2,000 percentile
resamples of the draws within each width give the slopes their intervals: fixed -0.649 [-0.714,
-0.596], -0.252 [-0.318, -0.187], -0.047 [-0.171, +0.201], -0.311 [-0.356, -0.267], -0.243
[-0.294, -0.193]; grown minus fixed at alpha = 1 +0.571 [+0.500, +0.642] and +0.658 [+0.596,
+0.722] at T = 0, +0.346 [+0.265, +0.428] and +0.301 [+0.232, +0.369] at T = 0.3 (every interval
excluding zero: the grown circuits do not decay); at alpha = 0.25, T = 0.3 the grown-L2 minus
fixed difference is -0.056 [-0.182, +0.087] and grown-V +0.006 [-0.097, +0.123], so "decays like
the fixed ansatz" is now a stored interval that includes zero. Predictions P4a and P4c held;
P4b failed on one of five cells (the alpha = 0.10, T = 0 fixed slope, whose interval is wide
because that gradient is supported on 13 of 128 states). Nothing in L35's direction changes.
Verdict: A4 gives Proposal A no width-scaling argument.

**A1, qubit-ADAPT-VQE in place of the fixed ansatz (L68; corrected by L75; Adversary L70:
STANDS WITH CAVEAT; verdict L69).** Pre-registered in `s26/PREREG_A1.md` (09:00, before any
endpoint existed). Build `s26/q_adapt.py --build --tag a1` (two shards killed by the host at
10:10 after 33 targets; one governed process resumed from the checkpoints, `s26/jobs_done/
a1_build.json`, 7,782 s, peak RSS 0.383 GB); label after the gate (`a1_label`, 25 s); stats
`s26/results/a1_stats.json`; per-target records `s26/results/a1/<pdb>.json` (126), each
reproducing the production quantum arm bit-for-bit (`ca`, `q_ca` and the selection index against
`bench_results/cache/464a0ddb5f283e04/`, 126 of 126). Basis: `rmsd_q_synth`, the production
projection of the weighted average over the 128 candidates, a BUILT CHAIN, on both sides of
every contrast; the comparator is the deployed selector `fixed_zrank_it50` (3.2280 A built chain,
3.3135 A single-window selection; the production top-75 arm is 3.2148 and the shipped argmin
3.4540). The two primaries, ADAPT at the same 21-parameter budget minus the fixed circuit:
pool L2 -0.0138 A (SE 0.0210, MDE 0.0588, 0.23x, fold CI [-0.0705, +0.0442], 58W/68L); pool V
-0.0222 (SE 0.0217, MDE 0.0608, 0.36x, [-0.0854, +0.0424], 61W/65L); both NOT MEASURED, 3 of 5
folds, concentration at the null's 50th and 51st percentile. The registered falsifier "ADAPT is
null" (both primaries inside 0.5x their MDE) fires; "ADAPT helps" does not. Single-window
secondaries -0.0437 and -0.0448 (0.47x MDE, fold CIs excluding zero, 45 to 59 exact ties). All
twelve ADAPT arms (two pools, Adam and L-BFGS, P = 7, 14, 21) are negative on the built chain,
-0.0128 to -0.0245 at 0.21x to 0.40x; controls: the fixed circuit at 750 steps -0.0035, the
exact Gibbs state +0.0088, uniform over 128 +0.0135, a matched-entropy random Hamiltonian
+0.0230 (fixed) and +0.0337 (ADAPT, fold CI [+0.0122, +0.0581], 5/5 folds, 0.34x). Power: the
comparison resolves 0.059 to 0.061 A; the whole quantum synthesis against the classical top-75
arm is +0.0133 A (`s26/results/q_mde_reference.json`), so the resolution is four times the
component's own footprint. The property half, no native: on the 78 alpha = 1 targets the fixed
circuit's KL to the Gibbs state is 0.9027 (max 0.984; S25's 0.902 reproduced) and ADAPT's 0.0002
(max 0.0009); KL(Gibbs || product of its marginals) is 1.4e-4 mean, 7.9e-4 max over 126 targets:
the deployed objective's optimum is a product state on every real target, a 7-parameter RY layer
reaches it, the 21-parameter fixed circuit stops 0.90 nats short, and reaching it moves the
emitted structure by -0.02 A, a third of the MDE (S25's readout-slack finding reproduced by a
second ansatz family). Correction L75: L68's "L-BFGS growth stops with no operator selected on
78 of 78 targets" is RETRACTED; the records (`adapt/*/sequence`, `trace`, `theta`) show operators
appended on 60 of 78 (pool V) and 68 of 78 (pool L2) targets under L-BFGS and on 78 of 78 under
Adam, all multi-qubit, buying at most 1.2e-4 nats (L-BFGS) or 8.6e-4 (Adam) with angles below
0.02 rad and leaving the state a product state to 4.1e-4 nats: inert growth, not absent growth
(on the ideal ladder of L27 and L35 the appended count is exactly zero, and those statements
stand). Adversary caveats (L70): the twelve arms' per-target delta vectors correlate at 0.955 on
average (first principal component 96.0% of their variance; the best-of-13 oracle -0.094 A
transfers +0.006 split-half), so "twelve of twelve" is one observation; the 21 parameters are a
budget, realised 7 to 21; the Gibbs control's +0.0088 is -0.0271 on the 78 alpha = 1 targets and
+0.0671 on the 48 alpha = 0.25 targets; and 29 of 78 targets differ by more than 0.01 A (max
0.206) between two states 2.6e-4 nats apart, the projection amplifying sub-milli-nat differences
(the mechanism the tie-break floor measures, VII.4); after L75 the Adversary corrected its own
caveat 2 (L82; `s26/RETRACTIONS.md` R5): L-BFGS does append operators at alpha = 1, they are
inert, and nothing in the null depends on the count. The registered replication (seed 1,
reversed fold order; `s26/PREREG_A1.md` addendum 2; `s26/results/a1s1_stats.json`,
`s26/results/a1s1/<pdb>.json`, L139) did not reproduce the seed-0 null: at seed 1 the two
primaries are -0.0449 (SE 0.0227, MDE 0.0635, 0.71x, fold CI [-0.0960, -0.0037], 4/5 folds,
72W/54L) and -0.0523 (SE 0.0235, MDE 0.0659, 0.79x, [-0.1087, -0.0054], 77W/49L), the Type-M
zone; the seed-0 points lie inside the seed-1 iid CIs, and neither seed clears its MDE. The
mechanism is in the arm means: the ADAPT arms are seed-stable (3.2142 to 3.2161 and 3.2059 to
3.2087; they converge to the same product Gibbs state whatever the start) while the fixed
21-parameter comparator, which stops 0.90 nats short of that optimum in a seed-dependent place,
moved from 3.2280 to 3.2610 (+0.033 A) with the seed, so the contrast grew because the
comparator landed worse, not because ADAPT emitted a better structure; the production top-75
arm is seed-free (3.2148 on both). A1 is therefore NOT MEASURED on either seed: underpowered at
seed 0 (0.23x / 0.36x), Type-M at seed 1 (0.71x / 0.79x), direction consistent, magnitude set by
the deployed circuit's initialisation; L68's "null at the registered threshold" is a seed-0
statement (`s26/RETRACTIONS.md` R11, a scope correction; `s26/PROPOSAL_A.md` addendum 3). The
Adversary's reading (L140): two seeds are one ADAPT value against two draws of the comparator,
so no pooling across seeds is licensed and by S20's rule a variational arm needs four seeds
before its endpoint is called; the comparator's 0.033 A seed variance is a fact about the
deployed selector worth a qualifier wherever its 3.228 is quoted; and REPLACE stands on
seed-independent facts (the product-state optimum, the full so(2^n) algebra from depth 2, no
width-scaling argument, inert growth) and on an endpoint not measured on either seed. What seed
1 adds is a sentence about the deployed circuit, not about growth: a better-optimised state of
the same product target (the 7-parameter RY layer, -0.050 at seed 1, or the exact Gibbs state,
-0.024) sits 0.02 to 0.05 A nearer at 0.3 to 0.8x MDE because the fixed circuit is a
seed-sensitive under-optimiser; that supports "replace", not "grow". A four-seed run is the
next step and was not run. Verdict for Proposal A (`s26/PROPOSAL_A.md`,
accepted L69): REPLACE, not the expected keep-with-edits, because the mechanism is absent rather
than weak: the Hamiltonian is a constant ladder whose optimum is a product state, the fixed
ansatz's algebra is already the full so(2^n), the grown circuits give no width-scaling argument,
and the endpoint is not measured on either seed at a stated resolution of 0.06 A.

**A3, a target-dependent Hamiltonian (L125; the S25 L17 question).** `s26/PREREG_A3.md` (filed
09:00); build `s26/q_adapt.py --build --tag a3` with the variants `zrank`, `zraw`, `asinh`, `soft`
and ADAPT on `zrank` and `zraw` (`s26/jobs_done/a3_build.json`, 16,373 s, peak RSS 0.352 GB);
label after the gate (45 s); stats `s26/results/a3_stats.json`; records `s26/results/a3/<pdb>.json`
(126, each bit-for-bit against the production quantum cache). Every variant is a strictly
increasing function of the shipped score over the same 128 candidates, re-standardised, so the
classical order and the CVaR prefix are the same on every arm; `Tmatch` runs the fixed circuit
at the per-target temperature at which the variant's Gibbs entropy equals the deployed 4.914
bits, chosen native-free. Built chain on both sides against the deployed selector
`fixed_zrank_it50` (3.2280): entropy-matched, `zraw` +0.0034 (SE 0.0273, 0.04x MDE, fold CI
[-0.0446, +0.0337]), `asinh` +0.0145 (0.18x), `soft` -0.0037 (0.05x); unmatched at T = 0.3, the
sharper states are worse on all three variants, `zraw` +0.1053 (0.77x, fold CI [+0.0367, +0.2112],
5/5 folds), `asinh` +0.0864 (0.69x), `soft` +0.1116 (0.81x), the Type-M zone; ADAPT on `zraw`
+0.0345 to +0.0378 (0.30x to 0.32x); single-window `zraw` +0.0885 (0.55x). The registered
"null" fired for every entropy-matched arm, "helps" did not fire, "harmful" did not fire for
the unmatched arms (three variants, one direction, none past its MDE). The property half, no
native: under the rank ladder the trained states are a tight family (14 of 78 distinct at
alpha = 1 by a 0.01-nat criterion, pairwise median 0.010 nats; lane Q's prediction of at most 3
distinct states was falsified, and "two trained states" holds for the spectrum and the endpoint
but not at a 0.01-nat resolution), while under the entropy-matched raw score 124 of 126 are
distinct, a median 1.34 nats apart; the Gibbs entropy at T = 0.3 is 4.914 bits on every target
under `zrank` and 3.10 (range 0.01 to 5.99) under `zraw`; the matched temperatures average
0.565; KL(Gibbs || product) is at most 8e-4 under `zrank` and 0.050 mean (max 0.263) under
`zraw`; ADAPT's operator sets are 57 distinct on 78 targets under `zrank` (the commonest, {Y_0,
Y_3}, on 10) because the argmax is taken over pool gradients of order 1e-3. The L17 answer: the
Hamiltonian can be made target-dependent while preserving the selection semantics, and at the
deployed entropy that makes the trained states target-dependent (124 of 126 distinct) and the
emitted structure does not change (+0.003 A, 0.04x MDE); without the entropy match the same
gaps sharpen the state (2.8 bits against 5.9) and the built chain is worse by +0.09 to +0.11 A
with fold CIs above zero on 5/5 folds at 0.7 to 0.8x the MDE. The readout responds to the
entropy of the weights (S25 section 6.2) and to nothing else varied here; the spectrum was the
last untried lever on the selector's side and it moves the answer only through that entropy.
IDEA_l17 is closed: measured, not helpful; no cell of `VQE_LFO` changes.

### VII.2 Proposal B: the scaled generation lane (lane P)

**B1 is infeasible on this box (L13).** `s26/results/b1_feasibility.json`: ESMFold v1 needs at
least 8.76 GB resident (the 3B language model halved plus an fp32 trunk; 7.4 GB with everything
in fp16) against 4.4 GB of campaign headroom under the 93% ceiling (64.8% of 16.75 GB in use at
00:19); its `esmfold` module needs `omegaconf` and `openfold`, both absent, and `openfold` needs
`nvcc` and Python <= 3.9 on a Python 3.13 CPU-only box; the checkpoints are not on disk
(esmfold_3B_v1.pt 2,771,653,574 bytes and esm2_t36_3B_UR50D.pt 5,678,116,398 bytes, 8.45 GB to
download). Nothing was downloaded. B1 stops.

**B3, the repo-native stand-in (L14).** Its per-target artefacts exist: the sequence-only torsion
predictor (`s13/cache/tors_rows.npz['a_pepPos']`, `rmsd_build` mean 3.7705 A, built chain from
predicted torsions) and the constant helix (`s14/results/ladder.json`, 4.0648 A built chain).
Paired against the production built chain 3.2148: torsion predictor +0.5557 A (median +0.2314,
SE 0.1318, MDE 0.3694, 37W/89L); helix +0.8500 (SE 0.1494, 32W/94L); torsion predictor against
helix -0.2943 (SE 0.0807, 69W/57L); on FAIL18 the three emit 5.569 / 6.032 / 5.887 A. The
arm-choice ORACLE over {pipeline, torsion predictor, helix} is an order statistic (L14). The B3
question proper, whether the set of targets where the pipeline beats sequence-only is
characterisable native-free (L106, corrected by L107): `s26/PREREG_B3.md`; `s26/p_b3.py run` ->
`s26/results/p_b3.json` (job `p_b3_run`, 100 s, peak RSS 0.053 GB); 45 native-free features
(`s26/results/p_b3_features.json`: length, composition, the top-75 members' helix / strand / coil
content, distogram entropy, ESM contact-map and retrieval-score statistics), a nested
leave-fold-out ridge and a 300-draw label-permutation null; built chains on both sides
(pipeline 3.2148, torsion predictor 3.7705, constant helix 4.0648). The sign classifier is at its
permutation null against both comparators (held-out balanced accuracy 0.522 against a null 95th
percentile of 0.578 for the torsion predictor, p = 0.263; 0.557 against 0.566 for the helix,
p = 0.093); the regression of the gain's size reduces the held-out squared error by -0.5366 (SE
0.2309, 0.83x MDE, NOT MEASURED) against the torsion predictor and by -1.1420 (SE 0.3415, 1.19x
MDE, fold CI [-1.495, -0.700], 5/5, 94W/32L, Type-M zone) against the helix. The falsifier
required both halves; the size of the gain over a helix is partly predictable and the sign is
not. The feature that carries the size is the pool's strand content (standardised ridge weight
-0.612, rho -0.638 with the gain over the helix; L107 struck L106's "helix content and the
length"): a strand-like retrieval pool means a large gain over a constant helix, a helical one a
small gain, S17 L25's trivial target-level signal re-found. ORACLE stratum for the record: on
FAIL18 the torsion predictor beats the pipeline by 0.463 A and on the other 108 the pipeline
wins by 0.725 A, and nothing native-free in this feature set locates that stratum. **B3 verdict:
the falsifier fires; Proposal B's verdict is REPLACE** (the feasible-scale rungs B2 are the C2
ladder of VII.3, none of which clears its MDE in the helpful direction).

**B2, the feasible-scale rungs.** The ladder's rungs are evaluated under C2 (VII.3); the `raw`
rung's training was stopped by hand at the 21:50 stall (L74; fold checkpoints under
`s26/models/p_ladder/`) and resumes on the coordinator's word when the governor shows 2.5 GB
free; `esm8m`, `mix` and `pairnet` follow (L77). The B verdict waits on them.

### VII.3 Proposal C: the classical ladder (lanes P and PH)

**C1, the reproduction (`s26/PREREG_C1.md`; lane P findings section 7).** The in-band ceiling of
S17 reproduces from `s17/results/inband.json`: a perfect ranker inside the shipped top-25 returns
2.6087 A single window (Part III.4).

**C2, the prior ladder.** Training rungs as under B2. The anchor landed (L56, corrected by
L57): the shipped posterior pushed through the ladder's own path (`s26/p_ladder.py`,
`s26/results/p_ladder_shipped_s0.json`, 126 rows; job `p_eval_shipped2`, 416 s, peak RSS
0.13 GB) reproduces the production cache per target at 0.00e+00 on the single-window selection
(3.4540) and 7.25e-14 on the point cloud (3.0483), while its built chain lands at 3.2126, the
leaderboard-rebuild basis, with 120 of 126 targets differing from `rmsd_arm` by more than 1e-6
(mean abs 0.0107, max 0.1711 at 7JS6): the cloud enters the projection 1e-14 away from the
persisted `avg_ca` and the multi-start L-BFGS lands in a different local optimum on most
targets. Every rung and the anchor go through the same code path, so the ladder's paired
contrasts are on one basis, the rebuild built chain 3.2126, stated in every C2 entry; the
0.0022 A offset to the production figure is a basis difference, not an effect. The lam = 0
chain lands at 3.2052 against the production 3.2041 for the same reason. Rung evaluations landed so far, each through the same path, anchor `s26/results/
p_ladder_shipped_s0.json`, artefacts `s26/results/p_ladder_<rung>_s0.json` and
`p_ladder_report_<rung>_s0.json`, BUILT CHAIN on the rebuild basis (3.2126) and single-window
selection (3.4540) on both sides, rung minus shipped, n = 126:

| rung (what changes in the prior) | built chain: effect, SE, x MDE, fold CI, W/L, folds | selection: effect, x MDE, fold CI | verdict | ledger |
|---|---|---|---|---|
| `pca32` (the shipped prior retrained) | +0.0000 on 126/126, all bases | +0.0000 | reproduction gate 2 passes | L65 |
| `noesm` (the ESM block removed, 42-d physicochemical pair features) | +0.2078, SE 0.0733, 1.01x, [+0.0986, +0.3942], 47W/79L, 5/5 | +0.3299, 1.27x, [+0.2122, +0.4475] | WORSE, Type-M zone (S7-11's lost -0.288 to -0.34 reproduced with the opposite sign convention) | L62 |
| `conly` (the 13 contact-head columns, no embedding block) | +0.1223, SE 0.0708, 0.62x, [-0.0350, +0.2674], 53W/73L, 4/5 | +0.1116, 0.42x | NOT MEASURED; against `noesm` -0.086 built chain (0.51x) and -0.218 selection (1.00x) | L63 |
| `wide` (width 768, depth 4, 2.7x the parameters, same inputs) | +0.0317, SE 0.0458, 0.25x, [-0.0300, +0.1279], 54W/72L, 3/5 | +0.0972, 0.48x, [+0.0487, +0.1587], 5/5 | NOT MEASURED, null-to-worse: more capacity does not buy a better prior | L66 |
| `pca32f` (the 32-PCA refitted per fold, no global-PCA leak) | +0.0180, SE 0.0491, 0.13x, [-0.0725, +0.0970], 70W/56L, 4/5 | -0.0316, 0.17x | NOT MEASURED: the shipped global PCA was not a leak worth anything | L67 |
| `pca128` (128 ESM components instead of 32) | +0.0759, SE 0.0630, 0.43x, [-0.0691, +0.2157], 57W/69L, 3/5 | +0.0249, 0.13x | NOT MEASURED; S7-11's "indistinguishable" stands, direction if anything worse | L72 |
| `esm8m` (the 8M-parameter ESM-2 in place of the 650M) | +0.2418, SE 0.0740, 1.17x, [+0.1134, +0.3853], 5/5 | +0.2246, 0.88x, [+0.1221, +0.3046] (NOT MEASURED; L101 struck a wrongly typed "-0.30") | WORSE, Type-M zone; no better than no ESM at all (+0.034 against `noesm`, 0.15x): the channel is the 650M model | L99, L101 |
| `pairnet` (S19's triangle-update PairNet on the deployed ESM inputs, the pinned `pairnet_models/`) | +0.0407, SE 0.0486, 0.30x, [-0.0442, +0.1143], 3/5 | +0.0077, 0.03x | NOT MEASURED; the best MAE of the ladder (2.079 against the shipped 2.339 per pair) and the highest gamma-equivalent (+0.1285 at cosine 0.287) and still no endpoint gain | L103 |
| `mix` (the shipped posterior mixed with the pool's own 17-bin histogram, weight lam chosen leave-fold-out over eight values) | lam = 0 chosen on 5/5 folds, so identical to shipped (126 ties); the mean built chain by lam: 3.2126, 3.2286, 3.2295, 3.2555, 3.2871, 3.3451, 3.6556, 3.9455 at lam 0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1 | +0.0385, 0.21x | the shipped posterior is its own leave-fold-out optimum; typicality (lam = 1) costs +0.73; the per-target ORACLE over eight lams (-0.3988) is 183% accounted for by the best-of-8 null, split-half -0.1647, k_eff 5.81 | L93 |

Adversary check of the rungs (L79: all stand; L66 with a caveat): the phase gate is enforced in
code, natives enter only after selection, the top-75 is a stable argsort, `sel` is tie-averaged,
every rung goes through one path against one anchor whose built chain is the rebuild basis;
`noesm`'s +0.208 and +0.330 are Type-M magnitudes whose sign is measured and whose harm is
tail-carried (median +0.044, worst +3.65 at 8T61); `conly`'s "two thirds of the ESM channel is in
the contact head" is a point-estimate split; `wide`'s selection-basis "real, small degradation" is
not licensed at 0.48x MDE with the iid CI spanning zero (suggestive, not measured); and any
eventual best rung or per-target minimum over rungs is an order statistic to be read through
`best_of_k_within`, never as the raw minimum.

Every rung carries the S25 L12 caveat on its gamma-equivalent (`noesm` has gam_eff +0.1161 at
cos +0.247 while being worse: the product -2.1496 x gam_eff is redeemable only at cos = 1).
The ladder's decision (L112): no rung beats the shipped posterior on the built chain (`noesm`
+0.208 [+0.099, +0.394] and `esm8m` +0.242 [+0.113, +0.385] WORSE at 5/5 folds; `conly` +0.122,
`pca128` +0.076, `pairnet` +0.041, `wide` +0.032, `pca32f` +0.018 with fold CIs spanning zero;
`pca32` and `mix` the identity; `raw` pending), so the best rung is the shipped one and the C3
stage-2 delivery file `s26/results/p_best_rung_chains.json` (126/126, complete) carries the
production emission: `phi`, `psi` in radians, unwrapped (values up to 72.7 occur; wrap before any
use that assumes (-pi, pi], L114), `ca`, `rmsd_arm` (mean 3.2148); stage 2 is therefore the
stage-1 replication PH already ran (L87), verified in L118: `s26/ph_c3_verify.py` ->
`s26/results/ph_c3_stage2_verify.json` (126/126, 5 s) compares the delivery with the production
cache and finds `ca`, the wrapped phi and psi and `rmsd_arm` identical at 0.0 on 126 of 126
(mean 3.2147652), so stage 2 reduces to stage 1 by identity, no relaxation is run, and the
decision rule's outcome stands unchanged (`s26/C3_RESULT.md` addendum 4). The lane's own
re-projection of the same rung
(`p_deliver_shipped`, 567 s) is saved beside it as `p_best_rung_chains_rebuild_basis.json`
(mean 3.2126, the L57 basis) as the cross-check that the ladder's anchor and the production
chains are the same object up to the projection's multi-start sensitivity; it overwrote the
delivery file for six minutes and was moved (L116). The `raw` rung (the 1280-d ESM-2 embedding
with a learned projection) is recorded as NOT RUN (L133): four of five fold models exist
(`s26/models/p_ladder/raw_fold{0,1,2,3}_s0.pt`, 4,477 to 4,900 s per fold at peak RSS 1.93 to 1.96
GB; the first attempt was terminated at 19:11 with the governor dead and fold 3 queued 40 minutes
behind the cap), fold 4 is untrained at the 04:30 close and a rung needs all five; its inputs
are bracketed by `pca32f` (+0.018, 0.13x) and `pca128` (+0.076, 0.43x), both null, and S7-11
measured raw on selection as indistinguishable from both. The deferred
`coherence_penalised_training` (1.25 GB)
and `attn` (3 to 3.5 GB, only if the box empties); length- and FAIL18-stratified tables for every
rung; a replication for any rung clearing its MDE; the best-rung delivery file for C3 stage 2;
the final `PROPOSAL_B.md` and `PROPOSAL_C.md` (L77).

**C3, the relaxation against matched controls.** Native-free half, L24 (Adversary L48: STANDS,
every number recomputed: displacement 0.2198 A, SE 0.0084; e0 median 8.59e4 kcal/mol, 58.73%
above 1e4; e1 mean -559.8, max 1262.4 at 9KAR; 125 of 126 converged; relaxed bond mean 3.867 A,
6 of 126 outside [3.6, 4.0]): Part IV.7. Accuracy half, stage 1, L39 (Adversary L46: STANDS WITH
CAVEAT). `s26/ph_c3.py stage1` -> `s26/results/ph_c3_stage1.json` (complete 126/126; 10 s, peak
RSS 0.041 GB), pre-registered in `s26/PREREG_c3_control.md` sections 2 and 3 and addendum 1
(prediction: AMBER worse than the matched random move by about +0.013 A, cosine negative). No
AMBER compute: the production coordinates `ca` and `amber_ca` from the cache. Bases: the arm's
input is the BUILT CHAIN (3.2148) and its output the RELAXED CHAIN (3.2355), both reproduced
from the stored coordinates to 1e-6; the controls displace the built chain's CA trace by a
vector of the same per-atom RMS magnitude as AMBER's displacement (0.220 A mean, measured after
superposition), either isotropic with the six rigid components projected out (S16's
construction, 16 draws) or along the line toward a random member of the shipped K = 500 pool (16
draws); a displaced trace is not an ideal-geometry chain, so the controls match the operator's
magnitude, not its validity, exactly as S16 did. All contrasts are full-chain CA-RMSD to the
native, n = 126:

| contrast (all on the displaced-built-chain basis) | effect | fold CI | x MDE | W/L | verdict |
|---|---|---|---|---|---|
| AMBER (relaxed chain) minus do-nothing (built chain) | +0.0207 | [+0.0154, +0.0290] | 2.16 | 40/86 | WORSE |
| AMBER minus matched-magnitude random | +0.0111 | [+0.0062, +0.0171] | 1.12 | 52/74 | WORSE, Type-M zone |
| AMBER minus matched-magnitude toward-member | +0.0385 | [+0.0298, +0.0479] | 2.27 | 29/97 | WORSE |
| random minus do-nothing | +0.0096 | [+0.0077, +0.0122] | 2.28 | 29/97 | WORSE |
| toward-member minus do-nothing | -0.0178 | [-0.0255, -0.0124] | 1.49 | 90/36 | BETTER (provisional) |
| ORACLE cosine, AMBER minus random | -0.0503 | [-0.0735, -0.0253] | 1.17 | 74/52 | Type-M zone |
| ORACLE cosine, AMBER minus toward-member | -0.1719 | [-0.2104, -0.1388] | 2.65 | 94/32 | measured |

Every contrast is 5 of 5 folds in sign. The ORACLE cosine between AMBER's displacement and the
true residual is -0.049 (SE 0.015), positive on 36.5% of targets; the S16 identity predicts
+0.0107 for a move of this size orthogonal to the residual. Reading: about half the production
step's cost is the size of the move (a random move of the same 0.220 A costs +0.0096) and the
other half is its direction (AMBER is worse than its own random twin by +0.0111, sign and fold
CI clean, magnitude inflated about 1.07x); the displacement points away from the native and
reproduces S16 L27 to the fourth decimal (S16 cosine -0.0521). A move of the same size toward a
random pool member improves the built chain by -0.0178: a zero-information control, not a
proposal, provisional until its registered replication (`ph_c3_stage1_rep`, new seeds, reversed
order) lands; its mechanism is on the record (the projection moved the chain 0.166 A away from
the cloud with a negative cosine, and the pool members surround the cloud). Decision rule
(`s26/PREREG_c3_control.md` section 3): AMBER does not beat either matched control, so "refine
with physics" is dropped as an accuracy step and kept as a validity step, with the Adversary's
wording: a validity step on 124 of 126 targets; on 2BP4 (relaxed CA-CA 5.38 A) and 9KAR (4.86 A,
not converged) it breaks a virtual bond (L46). `s26/C3_RESULT.md` carries the sentence for lane P
and the presentation. Stage 2 repeats this on lane P's best C2 rung when it is delivered. Replication (L87;
`s26/results/ph_c3_stage1_rep.json`, new draw seeds, reversed target order, registered in
`s26/PREREG_c3_control.md` addendum 2): every contrast lands inside the first run's fold CI,
all 5/5 folds: toward-member minus do-nothing -0.0207 [-0.0250, -0.0166] (1.77x MDE, 94W/32L),
AMBER minus random +0.0100 [+0.0059, +0.0175] (1.02x, Type-M unchanged), AMBER minus
toward-member +0.0414 [+0.0331, +0.0501], random minus do-nothing +0.0107 [+0.0093, +0.0120],
the ORACLE cosine contrast -0.0459 [-0.0734, -0.0212]. The PROVISIONAL label on the toward-member
line is lifted: a zero-information move of the relaxation's own size toward a random pool member
improves the built chain by 0.018 to 0.021 A where the relaxation's move worsens it by 0.021; the
reading (a property of the projection's displacement, not of physics; not a proposal) stands.
Adversary check, L95: STANDS WITH CAVEAT; the toward-member gain is by the discipline a measured,
replicated, native-free 0.02 A, and the mechanism in numbers is that the built chain sits 0.166 A
further from the native than the cloud it was projected from, the pool members surround that
cloud, and a 0.22 A step toward one of them reverses about an eighth of the projection's own
displacement; the natural comparator (the same step toward the cloud itself) was not run, and a
moved chain is no longer an ideal-geometry chain. Not a lever; a measurement of what the
projection costs.

**The validity axis of the relaxation (lane PH, an L77 extension, L100).**
`s26/PREREG_validity_axis.md`; `s26/ph_validity.py run | report` -> `s26/results/ph_validity.json`
(126/126; AMBER job `ph_validity`, 2,568 s, peak RSS 0.274 GB). Gate: the production relaxation
re-run here reproduces the cached `amber_ca` and `amber_e1` on 126 of 126 targets to 0.0 A and
0.0 kcal/mol, so the panel describes the deployed emission before and after its own relaxation.
Native-free; the constant helix (zero clashes, ideal geometry by construction) printed beside
both. Built chain to relaxed chain, per target, n = 126: heavy-atom pairs below 2.0 A 0.4444 to
0.0079 (34 targets affected to 1; fold CI [-0.5546, -0.3172], 2.05x MDE, 34W/0L/92T); pairs below
2.6 A 3.4524 to 0.1190 (3.30x, 81W/0L); closest heavy-atom pair 2.3629 to 2.7846 A (4.62x, 115 of
126); bond strain 0.0 to 1.3% relative (12% on 2BP4), angle strain 0.0 to 2.5%, omega
non-planarity 0.0 to 6.61 degrees mean (48 on 1D6X; the cis fraction rises to 0.33%, two
targets); Ramachandran favoured 0.9302 to 0.9114 (0.84x MDE) and outliers 0.0242 to 0.0327, both
NOT MEASURED; chirality unchanged. The registered falsifier (clashes not halved) does not fire.
The sentence "validity step" rests on: the relaxation turns 34 emissions with a sub-2 A
heavy-atom overlap into 1 and 125 of 126 energies above the 1000 kcal/mol gate into converged
ones, at the price of 0.021 A of accuracy, 1.3% bond and 2.5% angle strain, 6.6 degrees of
peptide-bond non-planarity and a broken virtual bond on 2BP4 and 9KAR; it buys no Ramachandran
(the k = 10 restraint on N, CA and C holds the torsions and lets the force field fix the packing
by bending bonds and angles, S16 L27), and the constant helix beats the relaxed chain on every
validity axis, so a validity statistic a zero-information reference maximises is not evidence
about the force field.

**C4, the routers (lane P, L110, L115).** `s26/PREREG_C4.md`; `s26/p_c4.py run` ->
`s26/results/p_c4.json` (job `p_c4_run`, 2,839 s, peak RSS 0.114 GB). Anchor: the rebuilt m = 75
cloud reproduces `s12/results/agg_surface.json` and `s23/results/errdecomp.json` per target at
0.00e+00, the rebuilt built chain 3.2126 (the L57 basis). Labels: the per-m point-cloud RMSD of
the 15-rung surface and the per-target scale s* of S23; twenty-five new native-free features in
four blocks (the principal-axis spread of the top-75 cloud, retrieval-score entropy and gaps,
posterior bin entropy, ESM contact-map statistics) beside the S22 set, six blocks in all; a
nested leave-fold-out ridge, two m-router constructions (A: multi-output over the rungs then
argmin; B: ridge on log m* then the nearest rung) and an s-router (a global scale about the
cloud's centroid), label-permutation nulls (200 draws for m, 100 for s). The m* oracle over the
grid is -0.408 A and an order statistic (best-of-k share 1.43, split-half -0.099, k_eff 12.7).
Twelve m routers, point cloud against the fixed m = 75 (positive = worse): all new features,
router A +0.0684 (SE 0.0269, 0.91x MDE, worse than 96% of the permutation null; +0.0707 on the
built chain, 0.83x), router B +0.0267 (0.56x; +0.0237 built chain); the S22 features reproduce
their recorded null through the same harness (+0.0244, +0.0206); the other blocks +0.0594 to
-0.0012, the one negative number (retrieval-score entropy, router B, -0.0012 at 0.02x its MDE)
being "not harmful", not "helpful". Six s routers all predict s* with the wrong sign (rho -0.10 to
-0.21) and cost +0.007 to +0.027 against s = 1, as S23 L6 and L9 say they must (s* is a function
of the invisible common mode). C4 is closed: the sixth through seventeenth router constructions
land where the first five did (S22 L7, S23 L7), and the S22 L10 bound (sample size, not
features) stands as the explanation.

**C5, the common-mode prediction (L132).** `s26/PREREG_C5.md` (alpha grid reduced to 0.5, declared
in the commit); `s26/p_c5.py run` -> `s26/results/p_c5.json` (complete, 126/126; job `p_c5_run`,
5,284 s, peak RSS 0.09 GB), statistics `p_c5_stats.json`. Construction: per target the shipped
top-75 cloud c and the common mode c - t after one rigid fit (ORACLE, training folds only inside
the fits); R1 = five shell means of D(c) - D(t) in distance space, applied by stress descent; R2
= the cloud's own principal-axis frame (native-free); GLOBAL = the training-fold mean correction
per length band; RIDGE = a nested leave-fold-out ridge from the 45 B3 features to the five R1
shell means; RANDOM = a random direction of the same magnitude; ORACLE = the true correction.
Every arm through the production projection; built chain (rebuild basis 3.2126) on both sides,
point cloud carried; negative = better than the incumbent. GLOBAL R2 -0.0090 (SE 0.0140, 0.23x
MDE, fold CI [-0.0331, +0.0185]); GLOBAL R1 +0.0312 (0.60x, fold CI above zero at 5/5 but under
the MDE gate); RIDGE R1 +0.1643 (SE 0.0328, 1.79x, fold CI [+0.1176, +0.2253], 5/5, 40W/86L,
WORSE) against its magnitude-matched random control +0.1351 (1.33x): the predicted shell profile
is indistinguishable from a random one of the same size (the nested alphas went to the top of
the grid on 3 of 5 folds, so the ridge emitted a near-constant correction); RANDOM R2 +0.0123
(0.45x). The ceiling: subtracting the true common mode gives -1.8749 A in distance space
(120W/6L) and -3.1293 in the coordinate frame (126W/0L; the R2 oracle is the native by
construction, so its number is the size of the mode plus the projection's cost, not a target).
Point cloud: GLOBAL R1 +0.0360 (0.84x), GLOBAL R2 +0.0164, RIDGE R1 +0.1492 (1.91x, WORSE),
ORACLE R2 -3.0483. The pre-registered falsifier fires for both GLOBAL and RIDGE; C5 is refuted at
one agent-day, as S16, S19 L14 and S24 L7 predicted, now with the ceiling measured beside it on
the same operator: the common mode is 1.9 to 3.1 A of built-chain error and nothing native-free
in these two representations touches it (`s26/PROPOSAL_C.md` addendum 2).

### VII.4 The tournament entries

`s26/TOURNAMENT.md` (the Adversary, L51) scored the sixteen ideas the lanes filed as `s26/IDEA_*.md`
by information x plausibility / agent-hours, found no clean same-operator kill on the record for
any of them (every closer in the record measured a different operator or space), ranked the
survivors and set the run order. The outcome of every idea at the time of writing (the rows are
updated at the close):

| idea (lane) | tournament rank and plausibility | outcome | ledger | basis |
|---|---|---|---|---|
| product_state_optimum (Q) | 1, 0.80 | MEASURED: the deployed objective's optimum is a product state on every real target (KL to the product of marginals at most 7.9e-4 nats); reaching it exactly moves the built chain by -0.014 to -0.022 A at seed 0 and -0.045 to -0.052 at seed 1, not measured on either seed | L68, L70, L75, L139, L140 | built chain (`rmsd_q_synth`) |
| conformational_identity_floor (W) | 2, 0.70 | MEASURED, ORACLE: the same sequence in another deposit sits a median 2.908 A from the native; +0.97 A worse than the pool's best window, -1.44 better than its mean, both Type-M at n = 18 | L52, L80 | single window |
| tiebreak_noise_floor (W) | 3, 0.85 | MEASURED: the pipeline's own convention noise is 0.0039 A on the 126-mean and 0.0236 A as the paired MDE between two tie-break conventions | L64, L71 | built chain (rebuild) |
| strain_difficulty (PH) | 4, 0.50 | MEASURED as a calibration flag and RESTATED: the pool's own disagreement (the top-75's pairwise spread) predicts the emitted chain's error at partial rho +0.452; the relaxation's displacement (+0.433) is its proxy at rho 0.76 and adds +0.08 given it; the "first native-free quantity above 0.4" sentence retracted; not a lever | L53, L81, L121, L123 | built chain (ORACLE label) |
| branch_select (PH) | 5, 0.20 | CLOSED as an accuracy step: +0.0055 against the production choice (0.18x MDE); the relaxed energy picks the branch as well as the objective and no better (-0.102 against random) | L88, L96 | built chain |
| window_ensembling (P, run by W) | 6, 0.15 | CLOSED in the fixed-K form: -0.0005 against shipped (0.02x), +0.0042 against a zero-information resample | L84, L94 | built chain (rebuild) |
| l17_target_dependent_hamiltonian = A3 (Q) | 7, 0.15 | CLOSED, measured and not helpful: a target-dependent Hamiltonian at the deployed entropy makes the trained states target-dependent (124 of 126 distinct) and the readout does not notice (+0.0034 A, 0.04x MDE); the sharper unmatched states are worse (+0.09 to +0.11, 0.7 to 0.8x) | L125 | built chain (`rmsd_q_synth`) |
| rotamer_relief B (PH) | 8, 0.15 | CLOSED: relieving the builder's chi1 halves the catastrophic fraction (0.535 to 0.233 above 1e4) and leaves an energy that ranks nothing (rho -0.002) and still rejects harmfully (+0.030 against the anchor, +0.022 against random, 5/5 folds); Part A is L23 | L23, L131 | point cloud (reject); single window (ranking, ORACLE) |
| window_provenance (W) | 9, 0.15 | census MEASURED (73% of every pool and top-75 is fragment windows, 0.6% whole peptides); the ORACLE class gap 0.42 A in the pool; the achievable readout H_P3 CLOSED (+0.0066 against uniform, 0.15x) | L85, L94, L109 | single window (census contrast); built chain (H_P3) |
| amber_prior_partner (P, run by W) | 10, 0.10 | CLOSED: the leave-fold-out choice is lam = 0 on every fold; every mixture cell worse; the last AMBER form in the record | L105 | point cloud (cells); built chain (chosen arm) |
| amber_reject (PH) | delivered on the point cloud before the ranking | CLOSED on both bases: the physical-threshold reject with refill is +0.228 (point cloud) and +0.248 (built chain) worse than the shipped top-75, harmful at 1e3 and 1e4, underpowered null at 1e5 and 1e6 | L23, L43, L54, L86, L95 | point cloud; built chain |
| cis_peptide (PH) | delivered before the ranking | CLOSED as a lever: 0 cis on the instrument; the own-torsion floor 0.347 A is an upper bound and the tight floor is 0.083 A (0.043 with the prior off) | L22, L38, L48, L89, L97 | ideal chain against the native (ORACLE) |
| selfcopy_proxy_bound (W) | delivered before the ranking | DELIVERED: the 2/60 benchmark leak bounded MINOR without opening the benchmark (dev proxy 0.008 A; envelope 0.028 A mean CI, 0.151 A worst target, built chain) | L30, L44, L49, L55, L58, L108 | built chain and the others, each named |
| better_prior_inputs (P) | deferred on memory (`attn` 3 to 3.5 GB) | the C2 ladder's other inputs measured flat (VII.3); `attn` NOT RUN; `raw` NOT RUN (four of five folds trained, L133) | L112, L133 | built chain (rebuild) |
| coherence_penalised_training (P) | deferred on memory (1.25 GB) | NOT RUN; recorded as such | L51, L77 | |
| trainability_paper (Q) | not an experiment | the outline delivered (`s26/TRAINABILITY_PAPER_OUTLINE.md`); the direction the evidence favours (VII.6) | L117 | |
| partial_recall_gradient (W, an L77 extension) | added under the extension | no gradient at the pre-registered MDE; a weak proximity effect of order rho -0.2 suggested, confounded with retrieval | L90, L97 | built chain (ORACLE label) |
| memorisation_on_the_ladder (W, an L77 extension) | added under the extension | ORACLE calibration of the prior-ladder currency: the discounted currency under-prices a large trained move | L91, L97 | point cloud and built chain, ORACLE |
| validity_axis (PH, an L77 extension) | added under the extension | MEASURED: the relaxation removes the builder's clashes (34 targets to 1) at 1.3% bond and 2.5% angle strain and 6.6 degrees of omega, no Ramachandran gain | L100 | descriptive, built and relaxed chains |

**The AMBER reject filter.** Census, L23 (Adversary L48: STANDS; Part IV.5). Endpoint on the point
cloud, L43: `s26/ph_reject.py cloud` and `report` -> `s26/results/ph_reject_cloud.json` (complete
126/126), `s26/results/ph_reject_report.json`; pre-registered in `s26/PREREG_amber_reject.md`
sections 3 to 5 and addendum 1. BASIS: POINT CLOUD on both sides (the coordinate average of the
retained windows; the anchor reproduces the production 3.0483 to 1e-6). Arms: R rejects members
of the shipped top-75 with AMBER single-point energy above the threshold and refills from the
next-ranked survivors to m = 75; S rejects without refill; matched controls reject the same count
at random (RANDR, RANDS, 16 draws) or apply the threshold to a permuted energy vector (PERMR,
PERMS). At the primary threshold 1e4 kcal/mol, point cloud against point cloud, n = 126:

| contrast | effect | fold CI | x MDE | W/L/T | verdict |
|---|---|---|---|---|---|
| R minus the shipped top-75 (anchor 3.0483) | +0.2276 | [+0.1466, +0.3232] | 1.17 | 49/67/10 | WORSE, Type-M zone |
| R minus RANDR | +0.1672 | [+0.0699, +0.2793] | 1.26 | 47/77/2 | WORSE, Type-M zone |
| R minus PERMR | +0.1053 | [+0.0187, +0.2000] | 0.98 | 52/66/8 | NOT MEASURED |
| S minus anchor | +0.1079 | [+0.0647, +0.1528] | 1.34 | 38/75/13 | WORSE |
| S minus RANDS | +0.0709 | [+0.0317, +0.1087] | 1.00 | 41/72/13 | NOT MEASURED |
| RANDR minus anchor | +0.0605 | [+0.0206, +0.0962] | 0.57 | 62/62/2 | NOT MEASURED |
| RANDS minus anchor | +0.0371 | [+0.0229, +0.0495] | 1.47 | 35/78/13 | WORSE |

The falsifier fired the other way. The dose is monotone: at 1e3 (54.5 of 75 rejected) R costs
+0.532 A; at 1e4 (40.1) +0.228; at 1e5 (30.3) +0.093 (0.66x MDE); at 1e6 (23.3) +0.066 (0.51x);
the less the reject does, the less it costs, and its limit is the anchor. Rejecting the same count
at random and refilling is nearly free (+0.061, not measured), so the energy's choice of what to
reject costs the +0.167; the retained set at 1e4 is +0.254 A more expanded in Rg than the anchor,
S25's "AMBER selects expansion" reproduced as a filter. The threshold sweep is an order statistic
(`best_of_k_within`: R oracle minimum -0.3368, split-half transfer -0.1470, k_eff 3.74). Power:
every arm-against-anchor contrast at 1e3 and 1e4 clears its MDE with the fold CI excluding zero;
at 1e5 and 1e6 R is 0.66x and 0.51x its MDE, so "the mildest threshold is null" is underpowered
below a 0.13 A gain and measured against any harm above it. Disposition, pending the built-chain
twin (`ph_reject_chain`, running): the steric reject at a physical threshold, with refill, is
closed on the point cloud, with a sign; the filter form of the functional lever (S24 D1-C, S19 Q3)
stays closed and gains its sixth instrument. Adversary check, L54: STANDS WITH CAVEAT. The
operators are native-free (the native enters only in the ORACLE scoring), ties and the fixed
primary threshold are clean, and the conclusion is established; but the harm is tail-carried
(R at 1e4 against the anchor has median +0.0027 against mean +0.2276, p90 +1.28, worst +4.15 at
8T61: near zero on the median target, large on the minority whose whole pool has no survivor
and on the deep refills), the two headline magnitudes +0.228 and +0.167 are Type-M-zone numbers
whose direction rests on the measured 1e3 contrasts (+0.532, +0.196), S at 1e4 (+0.108) and RANDS
(+0.037), and the refill-cost / choice-cost split rests on an underpowered +0.061.

Built chain, L86 (the falsifier's last leg and the cross-basis replication): `s26/ph_reject.py
chain` -> `s26/results/ph_reject_chain.json` (complete 126/126), `ph_reject_report.json` (both
bases); job `ph_reject_chain` 11,976 s (3.3 h, 22.8 projections per target), peak RSS 0.09 GB;
`s26/PREREG_amber_reject.md` addendum 2. BASIS: BUILT CHAIN on both sides; every retained set is
averaged and projected through the production call, and the anchor is the shipped top-75
re-projected through the same call, 3.2126 against the stored 3.2148 (the rebuild basis of L57:
the production path round-trips the cloud through float32 before projecting and this instrument
does not; per-target difference mean -0.0021, max 0.171; CA deviation above 0.5 A on 6 targets,
the projection's own branch degeneracy). At 1e4, built chain against built chain, n = 126: R minus
anchor +0.2483 (fold CI [+0.1197, +0.3479], 1.17x MDE, 4/5 folds, 49W/67L/10T, WORSE, Type-M); R
minus RANDR +0.1574 ([+0.0480, +0.2811], 1.15x, 4/5, WORSE, Type-M); R minus PERMR +0.1021 (0.88x,
NOT MEASURED); S minus anchor +0.1037 ([+0.0373, +0.1662], 1.11x, 5/5, WORSE, Type-M); S minus
RANDS +0.0551 (0.73x, NOT MEASURED); RANDR minus anchor +0.0909 (0.75x); RANDS minus anchor
+0.0486 (0.99x). The dose is again monotone: 1e3 R +0.5631 ([+0.4413, +0.7037], 1.77x, 5/5,
MEASURED) and S +0.1921 (1.44x, 5/5); 1e4 R +0.248, S +0.104; 1e5 R +0.1190 (0.78x), S +0.0833
(0.95x); 1e6 R +0.0915 (0.65x), S +0.0511 (0.76x); the threshold sweep is an order statistic (R
oracle minimum -0.3442, split-half transfer -0.1473, k_eff 3.86). One difference from the point
cloud reported rather than smoothed: fold 1 has the opposite sign on R at 1e4 (-0.002 against
+0.267 / +0.226 / +0.421 / +0.310), so R is 4/5 folds here and 5/5 there; the fold CI still
excludes zero. The built chain is noisier by construction (the projection's branch noise on every
arm), so the matched-control contrasts that were Type-M on the point cloud are underpowered
here, and the L54 caveats transfer unchanged: the harm is tail-carried (median +0.002 against
mean +0.248, p90 +1.39, worst +4.20 at 8T61) and the refill-cost / choice-cost split is a point
estimate on both bases. Power: "the mildest threshold is null" is underpowered for a gain below
about 0.14 A on this basis and measured against any harm above it. DISPOSITION: AMBER as a steric
reject filter at a physical threshold, with or without refill, is CLOSED on both bases in the
two forms the record had not measured (physics-set count, refill), with a harmful sign at 1e3
and 1e4 and an underpowered null at 1e5 and 1e6; the functional lever's filter form keeps its
closure and gains a sixth and seventh instrument. Why, from L23: the members the threshold
condemns are condemned by the builder's side-chain placement (96.8%), not by their backbones,
and removing them removes compact members whose error cancelled in the average (S23 L5).
Adversary check, L95: STANDS WITH CAVEAT; on the built chain the R arm holds on 4 of 5 folds
(fold 1 -0.002), so it does not meet the standing rule's fold clause on this basis, while the S
arm (5/5) and the 1e3 arms on both bases carry "closed on both bases, with a sign".

**The cis-peptide gap.** Census, L22 (Adversary L48: STANDS, with stale line numbers for the step
gate, now `core/data.py:439` and `:725`): `s26/results/ph_cis_census.json`. Over the 1,507
peptide bonds of the 126 model-1 natives the omega deviation from planarity has mean 1.91
degrees, median 0.34, p90 5.8, p99 17.4, max 42.8 at 9UV5; 0.53% of bonds are beyond 20 degrees
and 0.13% beyond 30; 0 of 126 natives, 0 of 1,966 ensemble models and 0 of 2,352,893 library
windows are cis by either criterion (the universe's minimum CA step is 3.5045 A); the registered
ensemble-cis prediction is falsified. The two-bond-length projection is therefore worth 0.000 A on
this instrument (L22). Floor, L38 (ORACLE DIAGNOSTIC; Adversary L48: STANDS):
`s26/results/ph_cis_floor.json` (complete 126/126), `s26/PREREG_cis.md` part 2. The ideal-trans
rebuild of the native's own phi/psi sits 0.347 A from the native (SE 0.029, median 0.272, p90
0.752, max 1.474 at 1ID6; N/CA/C rebuild 0.340), above 0.5 A on 32 of 126; it tracks omega
non-planarity (Spearman +0.828 with the maximum deviation) and does not track the production
chain cost of 0.166 A (rho +0.083; chain cost against omega deviation -0.036). The paired contrast
chain cost minus floor is -0.1804 (fold CI [-0.2266, -0.1022], 86W/40L), meaning only that the
production projection pays less than the own-torsion rebuild floor; the floor is an upper bound
on the manifold floor (the tight `floor2` is registered), and the two Spearman nulls exclude
|rho| above about 0.25 at n = 126. Nothing here is a proposal. Tight form, L89 (ORACLE
DIAGNOSTIC; `s26/results/ph_cis_floor2.json`, 126/126; job `ph_cis_floor2` 2,813 s, peak RSS
0.104 GB; `s26/PREREG_cis.md` addendum 2; every RMSD in this paragraph is an ideal-geometry
chain against the native CA trace, the same basis on both sides): the native CA trace itself
projected through the production projection sits 0.0833 A from itself at the production rung (SE 0.0075, median 0.0596,
p90 0.215, max 0.435 at 6MBM, none above 0.5) and 0.0429 A with the torsion prior off; against
L38's own-torsion rebuild floor the tight floor is -0.2636 [fold -0.2998, -0.2289] (113W/13L, 3.46x
MDE), because rebuilding from the native's own phi/psi lets every 2-degree omega error accumulate
down the chain while fitting phi/psi to the trace absorbs them. So the constant omega and the
ideal bond geometry cost the projection at most a few hundredths on this instrument (the two
gaps are priced: cis 0.000 A, non-planarity about 0.04 A); the ramah prior at 0.3 pulls the
projection 0.0404 A [+0.0318, +0.0492] away from a native input (23W/103L), the prior's price on
the one input where it has nothing to fix; and the production chain cost exceeds the tight floor
by +0.0832 [+0.0440, +0.1331] (5/5 folds, 1.51x MDE): what the projection pays on the production
input is the operator's displacement of a non-native cloud (S16 L27), not representation.
L38's 0.347 A is retained as the own-torsion upper bound, and every quotation of it now carries
"own-torsion upper bound; the tight floor is 0.08 A, and the projection's 0.166 A cost is the
operator's displacement, not representation" (Adversary L97: STANDS).

**The identity null and the 2/60 benchmark leak.** The dev audit is L15 (Part II.1). Lane W's
bound, pre-registered as `s26/PREREG_selfcopy_bound.md` (L30; Adversary L49: sound, no benchmark
sequence, native or RMSD read anywhere), landed as L44 with the coordinator's ruling L50 on how
the report quotes it: artefacts `s26/results/w_selfcopy_endpoint.json`, `w_selfcopy_bound.json`,
`w_selfcopy_floor.json`, `w_selfcopy_retrieval.json`, `w_selfcopy_envelope.json`; job
`w_endpoint_report` (160 s, peak RSS 0.302 GB). Reproduction gate: the "with" arm equals the
production emission on all 126 (top-75 identical, cloud max abs 1.4e-14; means 3.4540 single
window, 3.0483 point cloud, 3.2126 built chain on the re-projection rebuild, 3.2052 lam = 0
chain). Bases named on every line; a positive delta means the leaked emission is worse.

- Channel A, the four dev self-copies, dropping the self-window from the pool and refilling
  (S10-4's operator): built chain 2P5H -0.0683, 6B9K +0.0071, 1CEK and 2FBU 0.0000; single
  window 0.0000 on all four and on all 122 controls (the BLOSUM rank-0 window is never the
  argmin). The matched control population (the rank-0 window dropped on the other 122 targets):
  |delta| p50 0.000, p90 0.078, p95 0.147, max 0.511; the four sit at its 61st to 89th percentile.
- C27 re-derived (S10-4's operator over the 13 targets with a >= 0.6 window): clean minus
  production +0.0004 A on the lam = 0 chain (fold CI [-0.0001, +0.0010], MDE 0.0012, 117 ties),
  +0.0018 A on the built chain [-0.0003, +0.0039], +0.0003 on the cloud, 0.0000 single window.
  The report quotes both, each with its basis (L50): the +0.0004 is S10-4's figure on the lam = 0
  chain, now with an artefact; the +0.0018 is the same operator on the production basis; neither
  clears its MDE and both say the dev leak is immaterial.
- Channel C, the ORACLE envelope: a fold model that had the target's own native among its
  training labels emits a built chain 0.6980 A nearer to it (fold CI [-0.8312, -0.5764], 112W/14L,
  5/5 folds; point cloud -0.7094; single window -1.2081 [-1.4377, -0.9745]; lam = 0 chain
  -0.6932), changes half the top-75 (overlap 0.47) and the argmin on 97% of targets; the paired
  gain of the architecture over the argmin is biased against the architecture by +0.5101 per
  leaked target. This is an ORACLE fact about leaked training, not a price of the leak that
  exists.
- Channel B, the carrier chain in the training labels, complete on the four carrier-out models
  (L108; 1A11 out of fold 2, 2LMF and 2P5J out of fold 4, 1U6V out of fold 0, each through lane
  P's exact pca32 training path with 276 / 231 / 120 / 120 pairs removed, 893 to 1,063 s and
  1.24 to 1.25 GB each; the reference models reproduce the pinned emission at 0.000; gated
  re-run `w_endpoint_report2`, 145 s, 0.309 GB). Reference minus carrier-out, built chain
  (positive = the carrier helped): 1CEK +0.0114, 2FBU -0.0021, 2P5H -0.2460, 6B9K -0.0265;
  with both channels removed, production minus clean +0.0114 / -0.0021 / -0.2460 / -0.0243.
  F3 is falsified (|delta| on 2P5H 0.246 above the registered 0.10 line; its control clause,
  carrier-out change against two control-out chains per fold, is measured on one control,
  9BAF out of fold 0; the five remaining control-out models and a second-seed retrain of the
  one large value were withdrawn from the queue by the starvation ruling (L111) and the hold was
  lifted 25 minutes before the close, so the control clause is OPEN and would need about two
  hours of wall on this box, L137). The sign follows the cross-deposit distance of the
  copy (L52): the 2P5J copy sits 2.33 A from 2P5H's native and training on it pulls the
  posterior toward the carrier's geometry (mean |dE[d]| 1.38 A per pair, the largest of the
  four) and the built chain 0.25 A away from the native; on 1CEK, whose carrier is 0.60 A from
  the native, it helps by 0.011. The leak is not a gift to the leaked target but a bias toward
  another deposit's conformation.
- Part E (ORACLE DIAGNOSTIC): the same sequence in a different deposit, 18 dev targets and 22
  verbatim partners, sits a median 2.908 A from the native (min 0.317, max 5.502; 18% under 1.0 A;
  the natives' own ensemble spread is 1.044 A). A verbatim copy is not a near-native answer at
  this length.
- The bound (2/60 times the per-target quantity, assumptions A1 to A5): dev-4 channel A, signed
  max, 0.0023 built chain, IMMATERIAL; the native-free triangle max 0.0067, IMMATERIAL; both
  channels removed, now n = 4, 0.0082 built chain (2P5H), 0.0082 on the paired gain, 0.0054
  cloud, 0.0002 selection, IMMATERIAL (L108; 0.0004 at n = 1); the own-native envelope at its fold-CI limit
  0.0277 built chain / 0.0479 single window / 0.0228 paired gain, MINOR. Pre-registered verdict:
  MINOR, bounded at 0.028 A built chain, 0.048 single window and 0.023 on the paired gain because
  the envelope clause fires; every direct dev-proxy measurement is 0.002 A or less. Against the
  benchmark's own CI half-width of 0.170 A the bound cannot change the benchmark verdict in either
  direction. The caveat that attaches to every benchmark figure now reads: 2/60 self-copies,
  bounded at 0.028 A (built chain) by the own-native envelope and 0.002 A by the dev proxy
  (`s26/results/w_selfcopy_bound.json`). Adversary check, L55: STANDS WITH CAVEAT; every row of
  the bound recomputes ((2/60) x 0.0683 = 0.0023, x 0.2016 = 0.0067, x 0.8312 = 0.0277, x 1.4377
  = 0.0479, x 0.6838 = 0.0228); the artefact's `gain` verdict disagreed with the ledger (it lacked
  the envelope's gain row) and the envelope limit is the fold-CI limit of a mean effect, not a
  per-target bound. Both applied in L58 (`s26/w_bound_addendum.py`; the artefact now carries the
  gain row and the per-target readings): the class MINOR holds under every reading of the
  envelope on the built chain (mean-CI limit 0.0277, p95 target 0.0830, worst target 0.1512 at
  2BP4), the point cloud (0.0280 / 0.0824 / 0.1587) and the paired gain (0.0228 / 0.0816 /
  0.1232); on the selection basis (0.0479 / 0.1154 / 0.1944 at 9KAR) the worst-single-target
  reading crosses the pre-registered 0.170 line, so for the benchmark's argmin mean alone a
  single self-copy behaving like 9KAR could move that mean by up to 0.19 A under A2; the paired
  gain, which is the benchmark's verdict, is MINOR under every reading (worst 0.123). The
  wording to carry: "2/60 self-copies; dev-proxy price 0.002 A; own-native envelope 0.028 A
  (mean CI) to 0.151 A (worst target) on the built chain, under A2; MINOR under every reading
  on the built chain and the paired gain; the selection basis reaches 0.194 at the worst
  target; cannot move the benchmark verdict either way", with "dev-proxy price 0.008 A" in
  place of 0.002 after L108; where the direct measurement is not zero it is harmful to the
  leaked target's built chain and neutral to its argmin, so the un-leaked benchmark paired gain
  would if anything be slightly more favourable to the architecture than the +0.0103 reported,
  by at most (2/60) x 0.246 = 0.008 A under A2.
- Measured on the way, native-free: a one-member change of the 75 can flip the medoid frame
  (cloud moves 1.29 A on 5H1H, 0.74 on 6EY3) and the projection branch (chain moves 0.85 to
  2.11 A on six targets at cloud moves of 0.04 to 0.09 A), recorded for the tie-break floor idea.

**The tournament ranking (L51).** `s26/TOURNAMENT.md`, the Adversary's ranking, sets the run
order under the evening's headroom: product_state_optimum rides A1; conformational_identity_floor,
strain_difficulty and tiebreak_noise_floor next; A3 chained behind A1; branch_select when 1.5 GB
is free; window_ensembling, rotamer_relief, window_provenance and amber_prior_partner as capacity
allows; better_prior_inputs `attn` (3 to 3.5 GB) and coherence_penalised_training (1.25 GB)
deferred on memory.

**Conformational identity floor (tournament item 2, lane W, L52; ORACLE DIAGNOSTIC).** Part E of
`s26/PREREG_selfcopy_bound.md` (falsifier F5: median above 1.5 A) and `s26/PREREG_identity_floor.md`
for the two paired contrasts; `s26/results/w_selfcopy_floor.json` (22 pairs) and
`w_identity_floor.json`. The same sequence in a different deposit sits a median 2.908 A from the
target's native (mean 2.811, min 0.317 at 5V5B, max 5.502 at 7S3O; 4 of 22 below 1.0 A, 6 of 22
below 1.5 A; per target over 18, median 3.055); the four cross-fold self-copies reproduce S24 L4
(1CEK 0.595, 2FBU 3.278, 2P5H 2.334, 6B9K 4.126). Single window against the native on both
sides, n = 18 targets: the copy is worse than the pool's own best window by +0.9677 A (fold CI
[+0.7247, +1.2809], 3W/15L, 1.06x MDE, Type-M zone) and better than the pool's mean window by
-1.4428 A (fold CI [-1.7172, -1.1896], 13W/5L, 1.14x MDE, Type-M zone); F5 holds. Sequence
identity buys a window about 2.9 A from the native, no nearer than the built chain's 3.21 A mean
by a margin the instrument can call at n = 18 and 1 A worse than the best window the pool
already holds; the self-copy leak is small where it is measured because the copy is not the
native (L44). Not done: a length-matched non-verbatim control population. Adversary check, L80:
STANDS WITH CAVEAT; both paired contrasts are Type-M, the number for the report is the median
2.9 A (a descriptive statistic at n = 22), and the sentence about the 2.0 A target is a statement
about these 18 targets, not a general law.

**Strain as difficulty (tournament item 4, lane PH, L53; a calibration flag, not a lever).**
`s26/PREREG_strain_difficulty.md` (committed before the run); `s26/ph_strain.py` ->
`s26/results/ph_strain.json` and the registered replication `ph_strain_rep.json` (new bootstrap
seed, reversed target order; 35 s each, peak RSS 0.1 GB); reads the 126 production records only;
`rmsd_arm` (built chain) is read only as the ORACLE label of the emitted structure. Four
native-free signals the production relaxation emits for free: the built chain's AMBER energy
before relaxation (`log_e0`), the energy removed (`log_drop`), how far the chain moved
(`moved`, the restraint RMSD on N/CA/C) and the strain left (`strain_after`). Partial Spearman
with the built-chain error given n and Rg: `moved` +0.433 (fold CI [+0.247, +0.588], 5/5 folds,
permutation p < 0.00025 against a Bonferroni bar of 0.0125); `log_e0` +0.241 and `log_drop`
+0.250 (4/5 folds, at the 0.25 bar); `strain_after` +0.133 (fold CI spanning zero). The
falsifier cleared on `moved` with margin and the replication lands inside the first run's fold
CI on every quantity (`moved` +0.433, fold CI [+0.248, +0.581]). Quartiles of `moved` (32
targets each), mean built-chain error: below 0.159 A moved, 2.286 A; 0.159 to 0.218, 2.936;
0.218 to 0.288, 3.758; 0.288 and above, 3.923. Without 9KAR and 2BP4 rho is +0.41; the top eight
by `moved` hold easy and hard targets alike. The FAIL18 Fisher tests are null (best one-sided p
0.113). It is the first native-free quantity in the record with a correlation above 0.4 to the
per-target error of the emitted structure (the routers of S22 L7 and S23 L7 reached 0.24 to
0.37 on the per-target sign); the pre-registration forbids turning it into a selector or a
weight, and the record says every such conversion fails held out. Adversary check, L81: STANDS
WITH CAVEAT, provisional: "replicated" means the bootstrap and permutation draws (the Spearman is
a deterministic function of 126 fixed rows), so say "CIs replicated"; and the pre-registration
carried no pool-disagreement control. That control ran (L121, `s26/a_strain_vs_spread.py` ->
`s26/results/a_strain_vs_spread.json`; ORACLE label, native-free signals): the shipped top-75's
own pairwise CA-RMSD spread, available before the relaxation runs, has partial rho +0.452 with
the built-chain error given n and Rg (fold CI [+0.280, +0.609], permutation p below 0.0005, 5/5
folds; the medoid's minimum +0.451); `moved` correlates with that spread at +0.756 and, given
it, adds +0.082 (iid CI [-0.103, +0.259], permutation p 0.39). Lane PH accepted the veto in full
(L123): the sentence "the first native-free quantity in this programme's record with a
correlation above 0.4" is RETRACTED, together with the framings "the physics reports when the
answer is untrustworthy" and "how far the relaxation moves the chain predicts its error"; the
finding is restated as: the pool's own disagreement is the native-free quantity that predicts
the emitted chain's error (rho +0.45 partial on n and Rg, 5/5 folds, zero cost, measurable
without AMBER), and the relaxation's displacement is its proxy (rho 0.76) adding nothing given
it. What stands from L53: the phenomenon, the quartile table as its presentable form, and "not
a lever"; the physics adds nothing the pool did not already say. The pre-registration lesson,
PH's own: the confound list of a difficulty signal must include the pool's own statistics
before the operator's.

**The tie-break noise floor (tournament item 4, lane W, L64; Adversary L71: STANDS WITH
CAVEAT).** `s26/PREREG_tiebreak_floor.md`; native-free half `w_tiebreak_draws` (5,166 s under a
four-job load, peak RSS 0.111 GB; `s26/results/w_selfcopy_tiebreak_draws.json`, 126/126,
production gate top-75 == `sub` on 126/126); gated half `w_tiebreak_endpoint` (25 s, 0.062 GB;
`w_selfcopy_tiebreak_endpoint.json`); composed in `s26/results/w_tiebreak_report.json`. Operator:
on each target the boundary tie class (every window sharing the 500th window's BLOSUM62 sum;
median 115 windows, 54.5 of them inside the production pool) is re-drawn uniformly 8 times,
everything downstream identical; the production chain here is the rebuild basis 3.2126. A
random tie-break replaces 3.57 of the 75 averaged members on average (median 3, max 15) and
leaves the argmin unchanged on 91.5% of cells; the chain moves 0.168 A per draw in the median
(p90 1.08), mostly orthogonally to the native. The floor: the sd of the 126-mean over the 8
draws is 0.0039 A on the built chain (0.0017 point cloud, 0.0032 lam = 0 chain, 0.0136 single
window); the paired MDE between two draws is 0.0236 A median over 28 pairs (max 0.0321) on the
built chain; per target s_tie median 0.0227, p90 0.129, above 0.1 A on 21 of 126, above 0.5 A on
none, range over the 8 draws up to 0.713 at 1D6X (a projection-branch flip). Falsifier (m_tie
below 0.002 and paired MDE below 0.010): not fired; registered predictions m_tie 0.003 to 0.010
and paired MDE 0.02 to 0.05 held, "about 8 of 75 replaced" was an over-estimate (3.6). The
production convention against the mean of the 8 draws: +0.0026 built chain (0.15x MDE), +0.0048
point cloud (0.55x), +0.0201 single window (0.75x, fold CI excluding zero, 107 ties; suggestive
of S25's non-neutral retrieval order, not measured). What it means, with the Adversary's
qualification: a hundredths-level effect is real only as a paired contrast with the tie-break
held fixed; quoted across runs or against another instrument it is inside the pipeline's own
convention noise (0.024 A at n = 126). Every hundredths-level effect on the record (C3's +0.0207
and +0.0111, S24's -0.022, the ORACLE weight 0.015, the +0.004 reranking, C27's two prices, L19's
0.012, L24's +0.013) was measured paired with the tie-break fixed and is not invalidated; S9's
0.004 to 0.005 A residual between the stable and the plain argsort is this floor's effect on the
mean. Not done: 16 draws; the relaxed basis.

**Branch selection (tournament item 5, lane PH, L88).** `s26/PREREG_branch_select.md`;
`s26/ph_branch.py solutions | relax | report` -> `s26/results/ph_branch_solutions.json` (126/126,
the production choice reconstructed to 0.0 A on every target), `ph_branch_relax.json` (126/126;
AMBER jobs `ph_branch_relax` and `ph_branch_relax2`, 41 + 85 cells, peak RSS 0.3 GB, about 65 s
per target), `ph_branch_report.json`. The production projection chooses among five solutions at
its lam = 0.3 rung (the warm start from lam = 0 and four generic starts) by the lowest objective;
the arm relaxes all five built chains with the production operator, reads the converged energy
with the restraint off and emits the BUILT chain of the lowest (ties averaged; none occurred). All
126 targets have at least two distinct solutions (4.77 on average). BUILT CHAIN on both sides:
the energy pick minus the production choice +0.0055 (SE 0.0109, MDE 0.0306, 0.18x, 39W/45L/42T,
NOT MEASURED; a gain of 0.03 A would have been seen); the energy pick minus a uniformly random
pick among the five -0.1021 [fold -0.1270, -0.0750] (2.29x, 91W/35L, 5/5, MEASURED); the random
pick minus production +0.1076 (2.65x); the raw single point's pick minus production +0.0922
(1.10x, Type-M, WORSE); on the relaxed-chain basis the energy pick is -0.0036 (0.09x). Means:
production 3.2126, energy pick 3.2181, raw-energy pick 3.3048, random pick 3.3202, ORACLE minimum
over the five 3.1301 (the -0.19 order statistic is 130% accounted for by the best-of-5 null,
split-half 56%, k_eff 4.63; each chooser finds the per-target best on about 30% of targets against
20% by chance). The measured sentence: the converged force-field energy discriminates among the
projection's own branches exactly as well as the torsion prior already does, and the unrelaxed
energy does not (the relaxation is what makes the energy usable, S20 L6). CLOSED as an accuracy
lever. Adversary check, L96: STANDS; the 56% split-half transfer is measured against the column
mean (the random pick, 3.320 A), so a fixed best start lands at about 3.21 A, the production
objective's own choice; the ORACLE 0.083 A between production and the per-target best branch is
the order statistic L88 says it is.

**Window ensembling (tournament item 6, the mandatory test-time-ensembling direction, lane W,
L84).** `s26/PREREG_window_ensembling.md` (20:05, before the probe); native-free half
`w_ensemble_clouds` (4,351 s under load with two suspensions, peak RSS 0.113 GB;
`s26/results/w_selfcopy_ensemble_clouds.json`, 126/126, production gate top-75 == `sub`, cloud max
abs 1.8e-15); gated half `w_ensemble_endpoint` (10 s; `w_selfcopy_ensemble_endpoint.json`). How
it differs from widening K: widening K (S17 L12) draws one shortlist from a wider pool and the
extra plausible windows displace near-native members out of that single top-75 (its ORACLE best
2.104 to 2.572 A from K = 75 to the full universe); fixed-K ensembling keeps three K = 500
shortlists (BLOSUM45, 62, 80), each cut to its own top-75 by the shipped score (the BLOSUM62 one
is the production set), and averages the three clouds in the production frame, so no shortlist
is widened and only the parallelism of the three clouds' errors can act. Geometry, medians: the
BLOSUM45 / 80 pools overlap the BLOSUM62 pool at 0.83 / 0.89 and their top-75s the production
top-75 at 0.83 / 0.88 (93 distinct members in the union); the clouds sit 0.14 to 0.20 A apart; the
ensemble moves the built chain 0.185 A from the production chain, the zero-information resample
0.249. Endpoints, built chain on both sides (rebuild basis), n = 126: the three-key ensemble
minus shipped -0.0005 (SE 0.0100, MDE 0.0281, 0.02x, fold CI [-0.0194, +0.0162]); minus its
matched zero-information control (a bootstrap resample of the production set, `boot3`) +0.0042
(0.14x); on the point cloud +0.0016 and -0.0008. Secondaries: BLOSUM45 alone -0.0039, BLOSUM80
alone +0.0237 (0.70x, fold CI above zero, Type-M: if anything worse), the K-variant ensemble
-0.0153 (0.40x; on the cloud -0.0116 at 0.52x with the fold CI below zero, suggestive and the S17
L12 direction), the union of three -0.0021, `boot3` -0.0047. Verdict: refuted at this instrument;
the design excludes a gain of 0.028 A or more on the built chain and 0.019 on the cloud; parallel
errors average to themselves (S23 L9, S24 L3). Test-time ensembling is closed in the fixed-K
form with a power statement; the widening-K form was closed by S17 L12. Adversary check, L94:
STANDS.

**Window provenance (tournament item 9, lane W, L85).** `s26/PREREG_window_provenance.md`;
census `w_provenance_census` (5 s; `s26/results/w_selfcopy_provenance_census.json`, 126/126, 0
unresolved windows); ORACLE contrast `w_provenance_oracle` (10 s;
`w_selfcopy_provenance_oracle.json`). Census, means over 126: whole peptides are 0.6% of the K =
500 pool and 0.7% of the top-75 (30 of 126 targets have one in the top-75), terminal peptide
windows 8.0% / 7.2%, interior 18.7% / 18.3%, fragment windows 72.7% / 73.8%; the score keeps the
class mix almost unchanged. ORACLE contrast, single window against the native on both sides,
per-target class means: in the pool any peptide-derived window is nearer the native than a
fragment window by -0.4161 [fold -0.4472, -0.3758] (3.18x MDE, 5/5 folds; whole+terminal -0.421,
interior -0.411; position in the parent does not matter, -0.011 at 0.13x); inside the shipped
top-75 the gap shrinks to -0.086 (0.96x, n = 126) for any peptide window and -0.1313 [-0.1646,
-0.1096] (1.13x, Type-M, n = 114) for whole+terminal windows: the score already removes most of
the class difference. By the pre-registration's rule the achievable readout test H_P3 (a
fragment-class relative weight chosen leave-fold-out against the uniform readout and a
permuted-weight control; expected 0.00 to -0.02 A against an MDE of about 0.05) therefore runs.
Adversary check, L94: STANDS WITH CAVEAT; the pool contrast (-0.416, 101W/25L, not concentrated)
is a clean ORACLE fact and the expectation for H_P3 was null (the uniform readout was optimal
over two weighting families in S25; 68% of the error is common-mode). H_P3 landed (L109): job
`w_provenance_readout` (9,321 s under the six-job load, peak RSS 0.110 GB;
`s26/results/w_selfcopy_provenance_readout.json`, 126/126; the uniform-weight cloud equal to the
production `avg_ca` on every target). The shipped top-75 averaged with fragment-class members at
relative weight w in {0, 0.5, 1, 2}, w chosen on four training folds' built-chain mean and
applied to the fifth (w = 0.5 on folds 0, 3, 4 and 0 on folds 1, 2: it always down-weights),
against the uniform readout and the same weights permuted across the 75 members (4 seeded
draws). Built chain on both sides: routed minus uniform +0.0066 (SE 0.0157, MDE 0.0440, 0.15x,
fold CI [-0.0155, +0.0289]); routed minus permuted -0.0179 (0.38x, fold CI [-0.0316, -0.0057],
4/5); on the point cloud +0.0001 and -0.0284 (0.68x, 5/5). Fixed weights over all folds, built
chain minus uniform: w = 0 -0.0043, w = 0.5 -0.0117 (0.39x), w = 2 +0.0183 (0.68x), the L85
direction and none of it measurable at n = 126. H_P3 refuted: the class information is worth
something against noise with the same weights and nothing against doing nothing; the routed
readout excludes a gain of 0.044 A or more; the 0.13 A single-window gap acts on 8% of the
members of a 75-member mean and the members' change moves the chain mostly orthogonally to the
native (L44, L64). What the census leaves for the report: the pool is three-quarters fragment
windows, a peptide-derived window is 0.42 A nearer the native before the score and 0.09 to
0.13 A after it, and re-weighting by provenance cannot convert that into an emission gain at
this n.

**Partial recall gradient (lane W, an L77 extension, L90).** `s26/PREREG_partial_recall_gradient.md`;
covariates `w_recall_cov` (15 s; `s26/results/w_selfcopy_recall_covariates.json`, 126/126): for
each dev target against its own fold model's training corpus, the maximum pinned identity (mean
0.470, range 0.333 to 0.588, all below 0.6), the maximum containment (mean 0.657; the 1.00s are
the four self-copies) and the longest shared substring (mean 4.4 residues, 3 to 13); endpoint
`w_recall_endpoint` (45 s; `w_selfcopy_recall_endpoint.json`). Spearman with the built-chain
error: identity -0.205 [fold -0.315, -0.099] (0.81x the registered MDE of 0.253; permutation p
0.010, Bonferroni alpha 0.017), containment -0.157, substring -0.088; on the single-window basis
identity -0.230 (0.91x). Verdict by the pre-registered rule: no gradient at the MDE (nothing
reaches rho -0.25), a clean bill for the 0.6 threshold at that resolution, with the honest addition
that a weak sequence-proximity effect of order rho -0.2 is suggested on every basis; a confound
the pre-registration did not name is that the fold model's training corpus is the retrieval
library, so retrieval alone predicts a negative rho (S7-12). Adversary check, L97: STANDS.

**Memorisation on the ladder (lane W, an L77 extension, L91; ORACLE DIAGNOSTIC).**
`s26/PREREG_memorisation_on_the_ladder.md`; job `w_ladder` (40 s, peak 0.274 GB;
`s26/results/w_selfcopy_ladder.json`, 504 (target, model) cells). The four own-native models of
L44 Part C, put through the C2 ladder's own currency: in probability space they travel 28.1% of
the way to the native's distances at cosine 0.584 (amplitude 0.473); in location space 77.3% at
cosine 0.911; the posterior's MAE against the native falls from 2.339 to 0.825 A per pair; the
measured endpoint deltas are -0.709 point cloud and -0.698 built chain (L44). The naive ladder
(-2.1496 times gam_prob) predicts -0.603, inside the MDE of the measured cloud delta; the
direction-discounted currency of S25 L12 (times the cosine) predicts -0.364, disagreeing by
+0.346 (1.42x its MDE, 5/5 folds; the registered falsifier fires). Over the 504 cells the
discounted prediction correlates +0.414 with the measured delta, gam_prob -0.436 and the MAE
change +0.721. Reading, for Proposal C's arithmetic only: for a trained posterior the endpoint
responds to the full probability-space move, so the -2.15 A per unit gamma currency is usable
for a trained prior when gamma is measured in probability space and quoted without the cosine
discount but with the cosine reported; the L12 discount was derived from re-readings that move
small distances at low cosine. Not deployable by construction. Adversary check, L97: STANDS; the
S25 L12 currency is calibrated on small re-readings only and must not price a large off-axis
move, as L62 and L79 already apply.

**AMBER as a partner to the prior (tournament item 10, lane W, L105).**
`s26/PREREG_amber_prior_partner.md` (22:15, before any real-target run); stage 1 native-free
`w_amberprior_clouds` (135 s, peak 0.163 GB; `s26/results/w_selfcopy_amberprior_clouds.json`,
126/126; lam = 0 asserted bit-exact against the shipped risk table), stage 2 gated
`w_amberprior_endpoint` (4,786 s under a six-job load, peak 0.131 GB;
`w_selfcopy_amberprior_endpoint.json`; per-cell contrasts `w_amberprior_cells_cloud.json`).
Operator: per pair, the K = 500 pool's 17-bin distance histogram weighted by exp(-beta x
zrank(E_AMBER)) from the 63,000 cached single points, mixed into the shipped posterior at weight
lam through a real `Distogram` risk table and the production path; lam in {0, 0.05, 0.1, 0.2,
0.35, 0.5} x beta in {0, 0.5, 1, 2, 4}, chosen leave-fold-out on the point cloud, verdict on the
built chain. The leave-fold-out choice is (0, 0) on every fold, so the deployable arm is the
incumbent bit-exactly on 126 of 126 and every registered contrast is an exact tie. Cell by cell
on the point cloud every mixture is worse than the shipped posterior, monotonically in lam
(+0.016 [+0.002, +0.030] at lam 0.05 to +0.145 [+0.001, +0.260] at lam 0.5 with beta 0; 22 of 29
cells with the fold CI above zero, every cell NOT MEASURED by the rule), and beta moves a cell
by at most 0.02 A with no consistent sign: AMBER's ordering carries nothing through the prior
that the unweighted histogram does not (S24 L16's rank-permuted null seen from the prior side).
The per-target ORACLE over the 30 cells (-0.290 on the cloud) is 224% accounted for by the
best-of-k null (k_eff 9.7), split-half transfer -0.022. The last AMBER form in the record, a
distribution inside the prior, is closed on this instrument; registered expectation (0.00 +/-
0.03, lam = 0 on most folds) held.

**Rotamer relief B (tournament item 8, lane PH, L131).** `s26/PREREG_rotamer_relief.md`;
`s26/ph_relief.py run | report` and `ph_relief_reject.py` -> `s26/results/ph_relief_run.json`
(126/126; AMBER job `ph_relief_run`, 6,890 s, peak RSS 0.313 GB), `ph_relief_report.json`,
`ph_relief_reject.json`. The operator, native-free: one greedy sweep over every residue with a
chi1 (all but glycine, alanine and proline), options {the builder's default, 60, 180, 300
degrees}, the lowest ff14SB/GBn2 single point kept per residue, backbone fixed, no minimisation;
41 single points and 0.67 s per member; the relieved energy is at most the raw one on all 9,450
members. Census: the fraction of the shipped top-75 above 1e4 kcal/mol falls from 0.535 (SE
0.026) to 0.233 (SE 0.024), and above 1e6 from 0.311 to 0.057; the registered falsifier (not
below half of the raw fraction, 0.27) does not fire: more than half of the catastrophic single
points were the builder's chi1 choice, as L23 predicted, and about a quarter of every pool stays
above 1e4 (backbone-owned or multi-rotamer clashes a one-pass sweep cannot reach). Does the
relieved energy rank (ORACLE evaluation)? No, and neither did the raw one: Spearman with the true
RMSD within the top-75, raw +0.0000 (SE 0.0202), relieved -0.0018; in-band +0.0043 and +0.0055;
the difference -0.0019 at 0.04x MDE. Does it reject helpfully? No: at 1e4 the relieved threshold
rejects 17.5 of 75 (raw 40.1), and the shrink arm on the relieved energy is +0.0297 worse than
the anchor on the point cloud (fold CI [+0.0188, +0.0404], 5/5 folds, 0.74x MDE: the sign
measured, the size underpowered), +0.0220 worse than rejecting the same count at random
([+0.0106, +0.0318], 5/5, 0.61x), and less harmful than the raw reject by -0.0782 (1.12x, Type-M)
only because it rejects a third as many members; at 1e6 +0.0188 against the anchor (0.59x). The
built-chain half was not run: a point-cloud arm worse than the anchor with a fold CI excluding
zero would need the projection to reverse its sign, which no arm in L86 did; the refill arm was
not run, as registered. Closed: the relief removes the builder's half of the singularity and
leaves an energy that ranks nothing and rejects harmfully, the last route by which an all-atom
single point on retrieved windows could have entered the pipeline. For the report: half of what
the force field called impossible was the way the side chains had been placed; once that is
fixed, the force field still cannot tell a good backbone from a bad one.

**Coherence-penalised training and `better_prior_inputs attn`.** Not run: deferred on memory all
night (1.25 GB and 3 to 3.5 GB against a headroom that never freed) and recorded as such.

Lane W's closing entry (L137) lists its twelve ledger entries and 32 governed jobs (31 exit 0;
the one exit 1 the synthetic test that caught a bound error before any real run; largest peak
RSS 1.250 GB, longest wall 9,321 s), confirms that no benchmark file, sequence, name, native or
RMSD was read at any point, and names what it did not run with the time each would need: the
Part B control-out models (about 2 h), the second-seed 2P5J retrain (about 21 min), the tie-break
floor at 16 draws and on the relaxed basis, the ensembling and provenance replications (each
the pre-declared replication for a positive that did not occur), the recall gradient's gated
mechanism check, the triangle bound on the benchmark itself (a coordinator's option, closed for
this sprint), and the two ideas deferred on memory.

### VII.5 Operational findings of the sprint (for the record)

The governor v2 band held through the sprint until the host killed it for low memory at 10:10
(L40); jobrun v2.1 re-checks the launch cap after a jitter after six jobs launched against a cap
of four (L36); v2.2 refuses a stale governor snapshot for any job above 0.5 GB and any launch
without est-ram + 0.5 GB free (L41). The campaign paused at the user's request from 09:30 to
19:14 (L37, L42); a second report, `docs/REPORT_S26.md`, appeared during the pause and is
reconciled in Appendix D. The governor v2 died again at 19:36 on a Windows file-replace race and v2.1 retries the
snapshot write (L59); at 21:50 the box stalled at 91.9% RAM with every job suspended and v2.2
gained the stall breaker (L74). The deck: lane PR built `vqe_research_overview.pptx` from
artefacts (11 slides; 237 registered numbers, then 296 and 302, each with its path; four typed
from the claim ledger, seven kept off the slides; L61), rebuilt slide 8 in proposal form with the
L70 caveats and flagged the one L68 wording the A1 records state differently (L73), and applied
L75's final wording (L76; verification 0 dashes, 0 banned words). The Adversary's checks of L27,
L39, L35, L22 to L24, L38, L30, L43, L44, L68 and L64 all stand (L45 to L49, L54, L55, L70, L71),
six with caveats recorded above, and the Adversary's checks of the six C2 rungs (L79), the
identity floor (L80), strain (L81, then L121), ensembling and provenance (L94), the reject on
the built chain and the C3 replication (L95), branch selection (L96), the tight floor, the
recall gradient, the memorisation ladder and the mix rung (L97), and the two verdict files and
the ruling (L120) likewise. The deck's final build (L126, L134): slides 9 and 10 from the final
`PROPOSAL_B.md` and `PROPOSAL_C.md`, slide 11 with L117 as amended by L122, L123's strain
wording, L119's intervals and one A3 sentence on slide 8, C5's closure and the `raw` rung's
status on slide 10; 520 registered numbers each read from an artefact at build time; nothing
live, DRAFT or pending; verification 0 dashes, 0 banned words, spoken words 249 / 248 / 249 on
slides 8 / 9 / 10 against a limit of 250. The coordinator appended a "Sprint 26 additions"
table to `docs/FINDINGS.md`'s corrections ledger (L130): R2's scope, S7-11's lost artefact
re-measured (L62), S10-4's leak prices re-sourced with the 2/60 bound, the suite's basis named,
S16's relaxation finding confirmed and sharpened, the routers extended, and S26's own corrections
(R1, R5, L123, L101, L107, L57, L9); no sprint body edited. The Adversary's final deliverables
pass (L135, 04:00 to 04:12; `s26/DELIVERABLES_CHECK.md`): 13 of 14 rows OK with this report's
row PARTIAL at 04:00 and re-checked at 04:40; 29 pre-registrations each filed before compute;
18 ideas and the tournament; all eight findings files; 137 ledger entries; R6 to R10 added
beside R1 to R5 with the prior-sprint disposition table; the repository fixes (items 1 to 5,
defects 6a to 6d, `s24/cache_amber` tracked, the opt-in tier 11 of 11, the frozen rebuild 2016
of 2016); the deck with the Adversary's qualifiers on the slides; the paper outline with one
minor (the literal "at depth 3" in a caption); Rule 13 at 0 / 0 / 0 on every presenter-facing
file and the deck dump. The two checks the Adversary had listed as not done landed as L136: C4's eighteen routers
(nested, permutation-nulled, both bases, no order statistic taken) and the C3 delivery file
(the production emission, cross-checked by L116) both STAND, and L135's "not done" line is
closed. Lane I's re-run of the standalone audits under the governor
(`s26/i_verify_rerun.py`, records `s26/results/verify/`, table `s26/results/verify/REPORT.md`,
each audit executed in-process with its writes redirected and the tracked JSON's sha256
asserted unchanged; the 02:59 render, twelve audits): `vqe_lfo_audit`, `cvar_audit`,
`ansatz_audit`, `legacy_audit` and `amber_audit` IDENTICAL leaf for leaf (32, 18, 19, 30 and 31
leaves); `headline_audit` identical on 123 leaves with 8 added; `leak_audit` differing on one
leaf, the projection backend's name (`s8.project` tracked, `core.project` fresh);
`grad_key_collision` differing because the collision the tracked audit recorded (`COLLISION`
true, the config had no field for the gradient) is fixed on the current tree (false; the field
exists); `equiv_compare` and `project_arms` differing because the fresh runs compare the
production cache `1fc9f2dcf489e2fb` on 126 targets against the 8-target keys the tracked files
were built from (the fresh comparison bit-identical on 126 of 126 for `ca`, `fit_ca`, `phi` and
`avg_ca`, worst 0.0, where the tracked 8-target run had 0 of 8 at up to 22 A on the pre-repair
arm); `amber_platform` and `projection_divergence` with no fresh output (the second raising a
shape mismatch); two not re-runnable (`recon_containment_audit.json`, written by a script
deleted in the consolidation, and `project_inputs.json`, the harvested input the projection
audits read). The final table (L141; runner `s26/i_verify_rerun.py`, deletions and move-outs under `verify/`
refused, the tracked sha256 asserted unchanged after every run): 20 of 20 audits ran, 0 not run;
5 IDENTICAL, 11 DIFFERS each explained, 3 with no tracked JSON to diff (`determinism_audit`,
`hazard_audit`, `project_selfcheck`, whose fresh outputs are recorded: bit-identical across
processes and under 4 threads, no key collisions, the stable argsort, builder against reference
1.732e-13 A), 1 ERROR by construction (`projection_divergence`, whose arm caches are gone and
whose question `project_arms` answers). Every audit whose tracked JSON records a comparable
experiment reproduces its science: `project_exactness` 126 of 126 bit-identical (worst 0.0; the
one differing leaf a timing), `project_equiv` on every science leaf of its four-arm table
(synthesis 3.2148 / 3.2148 / 3.2145 / 3.2057 for the reference, exact, finite-difference and
analytic arms), `project_stability` identical on all 469 science leaves, `equiv_compare`
bit-identical on 8 of 8 including every AMBER quantity; two audits record the repair of what
they found (`grad_key_collision` true to false, the mode now in the cache key; `project_arms` 0
of 8 to 126 of 126 once the exact projection replaced the scan builder); two are not like for
like per target because the module default changed under them (`project_iters`,
`project_degeneracy`) and say so; `amber_platform`'s tracked JSON is a hand-assembled summary
with no writer. The sprint closed at 04:14 with L142: what stands and where it lives, in seven
items (Phase 0 signed off with the production numbers reproduced to 0.0 on all four bases; the
three verdicts; the tournament with every ranked survivor run or recorded as not run; the
repository items and the suite at 370 / 368 / 2; the deck, the paper outline, this report, the
retractions R1 to R11, the deliverables check and the `docs/FINDINGS.md` block; the not-run list
with reasons; 184 governed job records, every one with an exit code and a peak RSS, the box
never above the 93% ceiling under the governor's control). The findings, not the folding, remain
the output: no proposal survived as stated, every mandatory direction has a measured answer, and
the presenter walks in knowing which questions are closed and by which sprint. After the close
the Adversary's report row (L144, on the 04:12 text) passed with two items this final pass
carries (the two-seed A1 statement; the check re-run on the final commit) and one quoted word
rephrased; lane PR applied L139, L140 and L138 to the deck (L143): slide 8 quotes both A1 seeds and says
not measured on either, with the deployed circuit's 0.033 A seed variance beside it; slide 5
carries the same qualifier; slide 4 states the 3.2148 A production number as seed-free; 550
registered numbers, verification 0 dashes, 0 banned words, spoken words 248 / 248 / 249. Lane I's extension
found two things about the test suite itself: the relaunched slow-test jobs had lost
`VERIFY_SLOW=1` because it lived in the launching shell, not in the command (one null run set
aside, both tiers relaunched with the flag inside the command, L98), and the opt-in equivalence
tier had never been able to run: its summary parser took the last brace of an indented JSON
dump and turned the parse failure into a skip, so the three gated equivalence tests were among
the "13 skips" of every recorded run; fixed in the test file (commit `7be8e4b0`, a parse
failure is now a failure) and relaunched (L102); once it could run the tier passed 3 of 3
(`pytest_slow_equivalence3`, 40.2 s, peak RSS 0.313 GB, L104): the legacy and consolidated arms
are bit-identical up to the projection on all 8 smoke8 targets and the projection's divergence
stays inside its pinned band, served from the two on-disk caches (so it certifies the
equivalence of the two recorded arms, the same fact as `compare_tuning126.json`'s
`science_delta` of 0), and the suite record became 370 tests, 360 passed, 0 failed, 10 skipped
with the tier in (357 / 13 without it); the integration tier (8 items) was still queued. One
housekeeping fact: eight lanes share one git index, and a lane's commit once swept four files
another lane had just staged (L104). The final AST gate against `ae86a124` passed:
50 of 55 production files AST-identical with docstrings stripped and the 5 that differ exactly
the ledgered edits (`s26/results/ast_gate_ae86a124.txt`, L98); the frozen results-lab rebuild is
queued (`resultslab_rebuild`) with its comparator self-checked on the untouched tree (2016 of
2016 RMSD values exact, 2142 of 2142 PDBs byte-identical). The rebuild ran (L127; 4,472.9 s under
seven concurrent jobs after 4,765 s in the queue, peak RSS 0.11 GB; the documented command
unchanged): REPRODUCED. Every per-target RMSD on both bases, 2016 of 2016 values, exactly equal;
every leaderboard mean, median, secondary mean, gate, correlation, paired effect, MDE, fold CI,
W/L and verdict identical for all eight configurations (production 3.2126 built chain / 3.0483
point cloud, pool gate WARN with 2 violations, difficulty gate PASS); every one of the 2,142 PDBs
identical in its ATOM records, differing only in the git-commit and module-hash remarks; the
one finding, that the tracked leaderboard's three descriptive columns (`selector`,
`hamiltonians`, `distogram_used`) were written by uncommitted code and rebuild as null (Part
IX.6). `results/summary/results.json`, `results.csv` and `leaderboard.json` were committed
(`083c9b95`); `results/structures` restored to the tracked bytes. Lane I was then starved for 75
minutes by the other lanes' relaunch cadence (jobrun has no priority); the cap went to seven and
non-critical launches were held (L111), and the integration tier ran: 8 passed, 0 failed
(165.6 s, peak RSS 1.139 GB): the 11 Legacy terms on a real pool column for column, the AMBER
System's parameters against one built from the same XML, the pinned 1A13 native interaction
energy bit-exact, single-point invariance under translation, NaN-poisoning through `run_target`
with the relaxation on, and bit-identity across processes and thread counts. With it the whole
opt-in tier is 11 run, 11 passed, and the suite record is 370 tests, 368 passed, 0 failed, 0
errors, 2 skipped (the two absent artefacts; commit `7ad4ef68`, L113). The extension (L77) adds, by lane: I's opt-in test tier
(`VERIFY_SLOW=1`), every `verify/` audit re-run, the results-lab frozen rebuild and the final AST
gate; Q's A3, the per-step DLA of the grown circuits, bootstrap CIs on the A4 slopes and a
second-seed replication of A1's null; P's remaining rungs and stratified tables; PH's reject
chain, branch selection, rotamer relief, the C3 replication and stage 2; W's ensembling,
provenance and Part B retrains; A's checks within the hour and the deliverables pass at about
04:00; this report's final pass at 04:15 for a 04:45 close.

### VII.6 The three verdicts and the direction the evidence favours (L117)

Each proposal was pre-registered with a falsifier before its first endpoint existed and judged
by the campaign's rule that a proposal is not softened to survive; the coordinator's ruling
(L117, 02:02) accepted the three lanes' verdicts with the evidence each rests on.

| proposal | verdict | file | what decided it |
|---|---|---|---|
| A: qubit-ADAPT-VQE in place of the fixed ansatz | REPLACE (the prompt expected keep with edits) | `s26/PROPOSAL_A.md` (addendum 3) | the deployed objective's optimum is a product state on every real target (KL to the product of marginals at most 7.9e-4 nats); an adaptive ansatz appends only inert operators; the fixed circuit's Lie algebra is already the full so(2^n) from depth 2; the grown circuits give no width-scaling argument; the endpoint is not measured on either seed (0.23x / 0.36x MDE at seed 0, 0.71x / 0.79x at seed 1, the difference the fixed comparator's own seed) at a 0.06 A resolution (L27, L35, L45, L47, L68, L69, L70, L75, L82, L119, L138, L139, L140) |
| B: a learned folding model as the prior | REPLACE | `s26/PROPOSAL_B.md`; the replacement `s26/PROPOSAL_B_REPLACEMENT.md` (the trainability paper) | ESMFold cannot run on this box (L13); no feasible-scale model input beats the shipped prior (`noesm` +0.208, `esm8m` +0.242, `conly` +0.122, L62, L63, L99); the set where the pipeline beats sequence-only is not characterisable native-free (L106, L107) |
| C: learn the ranking, then refine with physics | KEEP WITH EDITS | `s26/PROPOSAL_C.md`; `s26/C3_RESULT.md` | the edits: the model should learn a better prior, not a better ranker (the prior's derivative is steep, -2.15 A per unit toward truth, and its inputs on this machine are flat: nine rungs, none beats the shipped prior, L56 to L67, L72, L93, L99, L103, L105, L112); AMBER is a validity step only, worse than a random move of its own size (L39, L46, L87, L100); the routers are closed for the seventh and eighth time (L110, L115); C5 did not complete before the close and keeps its pre-registration |

The direction the evidence favours, in the presenter's words (the slide 11 line, verbatim from
L117): "If I could do one thing next, I would publish the trainability work first. Every figure
in it already exists as a measured artefact, it needs no new machine, and it is the one part of
this project whose result is positive and complete. The only open accuracy lever is the
distance prior, and the honest next step there is a larger language model than this laptop can
hold, so that comes second and needs a bigger machine. I would not spend more time on the
circuit for accuracy: we now know why it cannot matter here." Why this order: the paper's
inputs are all in hand (the S13 locality theorem and Pauli-spectrum prediction, the S25 width
sweep scoped to depth 3 by A2, A2, A4, the product-state fact) and the Adversary has checked
each; the prior lever is real (S24 L13) but every input this machine can compute is measured
flat, so it is a resourcing decision rather than an experiment the presenter can run next week;
and the circuit is closed as an accuracy lever by three independent facts. The Adversary's check
of the two verdict files and of the ruling (L120): PROPOSAL_B STANDS (every number in its notes
at its artefact; two minor sourcing gaps closed, the isolation contrasts now persisted as
`s26/results/a_ladder_isolations.json`: `conly` minus `noesm` -0.2183 single window at 1.00x MDE
and -0.0856 built chain at 0.51x, `esm8m` minus `noesm` +0.0340 built chain at 0.15x; the 8.76 GB
of B1 stated as its derivation, 12.0 / 2 for the fp16 language model plus 2.76 for the fp32
trunk); PROPOSAL_C STANDS after one material document fix (the file had carried a future
"FINAL 03:30" stamp and a past-tense C5 outcome before C5 had finished, a clock error of about
80 minutes in lane P's estimated times; corrected by a real-time header, a present-tense C5
section and addendum 1, L124); and L117 STANDS WITH three wording caveats, accepted by the
coordinator (L122): in the presenter's line "positive and complete" reads "exact and complete"
(the trainability chain is a measured identity, not a gain); "the Adversary has checked each"
applies to the S26 additions (A2, A4 with L119's intervals, the product-state fact) and not to
the S13 inputs, which the paper cites from their artefacts without an S26 re-check; and the S13
Pauli mean weights 2.236 / 3.015, which lane PR could not reproduce from `s13/results/
geo_pauli.json` on the night, stay off every slide and are quoted in the paper outline only with
their artefact path and that note (this report's Part V.10 quotes them from
`s13/SPRINT13_DOSSIER.md` section 8 as cited, with the same note). The line's "the one part
whose result is complete" is fair as a ranking of publishability, not as a statement that
nothing else finished.

### VII.7 S26 retractions (`s26/RETRACTIONS.md`)

Kept by the Adversary, one entry per retraction, the claim verbatim with the artefact that
contradicts it; nothing superseded is deleted anywhere.

| entry | claim | now stands | where |
|---|---|---|---|
| R1 | the governed test count 369 / 356 (`s26/EXAMINATION.md` D, C35; L28 item 2) | 370 tests, 357 passed, 0 failed, 13 skipped, 0 memory-guard skips (`s26/results/test_run.json`) | L31, L32, L34 |
| R2 | S25's "no barren plateau at any width measured" (the slopes of `s25/results/q_plateau.json`) | the slopes stand; the claim is scoped to depth 3, because the algebra is the full so(2^n) (A2); for `docs/FINDINGS.md` at the close | L27, L45; Part V.9 |
| R3 | lane Q's H2b, dim(DLA) at (n = 7, L = 3) below 8128, guess 4095 (`s26/PREREG_A2.md`) | 8128 = so(128) from depth 2; falsified by lane Q itself | L27; Part V.9 |
| R4 | lane PH's registered expectation of cis bonds in other ensemble models (`s26/PREREG_cis.md`) | 0 of 1,966 ensemble models, 0 of 126 natives, 0 of 2,352,893 windows | L22, L48; Part VII.4 |
| L123 | L53's "the first native-free quantity in this programme's record with a correlation above 0.4 to the per-target error of the emitted structure", and the framings "the physics reports when the answer is untrustworthy" / "how far the relaxation moves the chain predicts its error" | the pool's own disagreement predicts the error (partial rho +0.452, 5/5 folds); the relaxation's displacement is its proxy (rho 0.76) and adds +0.08 given it; the phenomenon and the quartile table stand | L121, L123; Part VII.4 |
| R5 (L75, L82) | L68's "L-BFGS growth at alpha = 1 stops with no operator selected on 78 of 78 targets" (and the same sentence in `s26/PROPOSAL_A.md` sections 3 and 5); the Adversary's L70 caveat 2 "L-BFGS grows nothing" | operators are appended on 60 to 78 of 78 targets and are inert (at most 1.2e-4 nats, angles below 0.02 rad, product state to 4.1e-4 nats); the verdict does not move | L73, L75, L76, L82; Part VII.1 |
| R6 (L121, L123) | L53's "the first native-free quantity above 0.4" (the same row as the L123 line above; the Adversary's entry for it) | the pool's own disagreement is the quantity; `moved` is its proxy | Part VII.4 |
| R7 (L57) | L56's "bit-exact on all four bases" for the C2 anchor | bit-exact on selection and the point cloud; the built chain lands on the rebuild basis 3.2126, 120 of 126 targets differing from the cache by more than 1e-6 | Part VII.3 |
| R8 (L101) | L99's "-0.30" for the `esm8m` selection slope | +0.2246 (0.88x MDE, NOT MEASURED) | Part VII.3 |
| R9 (L107) | L106's "helix content and the length" as the feature carrying the size of the gain over a helix | the pool's strand content (standardised ridge weight -0.612, correlation -0.638) | Part VII.2 |
| R10 (L9) | L7's point-cloud values 1.5921 (1D6X) and 2.2812 (1KWE), typed before the query returned | 2.0900 and 2.7313 (`bench_results/cache/1fc9f2dcf489e2fb/{1D6X,1KWE}.json :: rmsd_avg`); the L7 mechanism stands on the corrected numbers | S26 L7, L9 |
| R11 (L139, L140) | L68's "null at the registered threshold" for A1 | a seed-0 statement: at seed 1 the primaries are -0.045 / -0.052 A at 0.71x / 0.79x MDE with fold CIs excluding zero; A1 is not measured on either seed; a scope correction, no number withdrawn | Part VII.1 |

`s26/RETRACTIONS.md` also carries the disposition table of prior-sprint claims S26 contradicts,
scopes or re-sources, mirrored in the "Sprint 26 additions" block of `docs/FINDINGS.md`'s
corrections ledger (L130): the S25 plateau claim scoped to depth 3 (R2); the S13 Pauli mean
weights 2.236 / 3.015 not re-derived in S26 and kept off the slides; S7-11's lost -0.288 replaced
by the measured -0.330 single window and -0.208 built chain (L62); S10-4's leak prices
re-derived on the lam = 0 chain and priced on the built chain with the 2/60 bound (L44, L58,
L108); the seven-configuration suite's basis named as point cloud with its built-chain twins
(L28, L29); S16 L27's relaxation finding confirmed and sharpened (L39, L87, L100); the routers
extended by twelve m and six s constructions (L110, L115); and L38's own-torsion floor scoped as
an upper bound by the tight floor (L89).

Not retractions, recorded there as sourcing corrections: C27's +0.0004 / +0.0030 (re-derived
exactly on the lam = 0 chain by lane W, L44; the S10 artefacts are in git history at `5fa05cd`).
No prior-sprint accuracy claim has been contradicted by an S26 endpoint result to date; L39, L43
and L44 confirm the record.

## PART VIII. CLOSED AND OPEN

Two tables. The rule: an open item is never moved to closed on one failed experiment; it moves
when it has been measured to its own ceiling or refuted by a pre-registered falsifier, and the
row names which. The standing pre-S26 list is `docs/STATE_BRIEF_2026-09-12.md` 5.6 and
`docs/CONDENSED_REPORT.md`, "Closed - do not re-fund"; Part VI carries every closure with its
sprint. This part lists the S26 movements and the items that remain open at the close.

### VIII.1 Closed

| direction | closed by | where |
|---|---|---|
| Further in-band rankers over retrieval pools; post-hoc correction of the distance objective; growing the fragment library; multi-piece assembly; key fusion; either energy in a search objective | measured to ceiling or refuted, S5 to S13 | `docs/CONDENSED_REPORT.md`, "Closed"; Part VI.1 to VI.9 |
| Chemical-shift torsion restraints; sequence-only torsion prediction; degree-1 truncation; native-free error steering; learned candidate generation; Legacy-to-AMBER continuation; probability-weighted readout; per-target scale correction; the encoding lever; bond dimension; routers (five) | pre-registered falsifiers fired or ORACLE ceilings, S14 to S25 | `docs/STATE_BRIEF_2026-09-12.md` 5.6; Part VI.10 to VI.21 |
| The physics energies as selectors, filters or objectives (the functional lever, all five forms) | S17 L29/L30, S18 L17, S21 L33, S24 L16, S25 L16/L18 | Part IV.8 |
| The CVaR-VQE selector as an accuracy contribution on any readout | S16 L29, S21, S22 L13, S23 L8, S25 L3/L5/L15 | Part V.7, V.11 |
| The basis question (built chain is the result) | S25 L4, L8, L11 | Part II.2 |
| The ansatz's dynamical Lie algebra as an unmeasured property | measured: full so(2^n) from depth 2 at n = 7; independently re-derived | S26 L27, L45; Part V.9 |
| A width-scaling argument for Proposal A from grown (ADAPT) circuits | A4: product circuits at alpha = 1, fixed-ansatz decay at alpha = 0.25 | S26 L35, L47; Part VII.1 |
| "Refine with physics" (the production AMBER relaxation) as an accuracy step | C3 stage 1: worse than do-nothing, worse than a matched random move, worse than a matched move toward a pool member; replicated (L87); stage 2 reduces to stage 1 because the best C2 rung is the shipped prior (L118); kept as a validity step on 124 of 126 with its validity axis measured (L100) | S26 L39, L46, L87, L100, L118; Part VII.3 |
| The steric reject at a physical threshold, with or without refill, on both bases | falsifier fired the other way on the point cloud and replicated on the built chain (harmful at 1e3 and 1e4, underpowered null at 1e5 and 1e6); the sixth and seventh instruments for the filter form of the functional lever | S26 L43, L54, L86; Part VII.4 |
| The cis-peptide gap and the omega non-planarity gap as levers on this instrument | 0 cis natives, ensemble models or windows; the two-bond projection worth 0.000 A; the tight representation floor 0.083 A (non-planarity about 0.04 A) | S26 L22, L38, L48, L89; Part VII.4 |
| Branch selection by the relaxed energy as an accuracy step | +0.0055 against the production choice (0.18x MDE; a 0.03 A gain would have been seen); the energy discriminates among branches as well as the objective and no better | S26 L88; Part VII.4 |
| Test-time ensembling of retrieval keys (fixed K) | -0.0005 against shipped (0.02x MDE), +0.0042 against a zero-information resample; a 0.028 A gain excluded | S26 L84; Part VII.4 |
| A provenance-weighted readout (down-weighting fragment windows in the average) | +0.0066 against uniform (0.15x MDE; a 0.044 A gain excluded), -0.018 against the permuted-weight control; the ORACLE class gap (0.42 A in the pool) cannot be converted into an emission gain at this n | S26 L85, L94, L109; Part VII.4 |
| Routers for a per-target shortlist size or scale (C4: twelve m routers and six s routers on six native-free feature blocks) | none clears its MDE, all but one point the harmful way, the s routers predict the wrong sign; the sixth to seventeenth constructions land where the first five did | S26 L110, L115; Part VII.3 |
| ESMFold (B1) on this box | infeasible on three independent grounds | S26 L13; Part VII.2 |
| Proposal B: a native-free characterisation of the set where the pipeline beats sequence-only (B3) | the sign classifier is at its permutation null against both comparators; only the size of the gain over a helix is partly predictable, by the pool's strand content; verdict REPLACE | S26 L14, L106, L107; Part VII.2 |
| AMBER as a distribution inside the prior (the last AMBER form in the record) | the leave-fold-out choice is lam = 0 on every fold; every mixture cell worse, beta without a consistent sign | S26 L105; Part VII.4 |
| Proposal A: an adaptive (ADAPT) ansatz in place of the fixed one, as an accuracy or trainability lever | A1 not measured on either seed (0.23x / 0.36x MDE at seed 0, 0.71x / 0.79x at seed 1, the difference the fixed comparator's own seed; resolution 0.06 A; a four-seed run not done); the optimum is a product state on every real target; the fixed algebra is already maximal (A2, per growth step L138); no width-scaling argument (A4, with intervals L119); verdict REPLACE on the seed-independent facts | S26 L68, L69, L70, L75, L119, L138, L139, L140; Part VII.1 |
| A target-dependent Hamiltonian for the selector (A3, the S25 L17 question) | at matched entropy the trained states become target-dependent and the emitted structure does not change (+0.0034 A, 0.04x MDE); the sharper unmatched states are worse; the readout responds to the weights' entropy and nothing else | S26 L125; Part VII.1 |
| Native-free prediction and subtraction of the pool's common mode (C5) | global corrections null (0.23x and 0.60x MDE), the learned correction harmful (+0.164 A, 1.79x, 5/5 folds) and indistinguishable from a random correction of the same size; the ceiling 1.9 to 3.1 A measured beside it | S26 L132; Part VII.3 |
| An all-atom single point on retrieved windows after relieving the builder's rotamers (rotamer relief B) | the relief halves the catastrophic fraction and leaves an energy that ranks nothing and rejects harmfully; the last route for an all-atom single point on windows | S26 L131; Part VII.4 |
| More prior capacity (`wide`), a per-fold PCA (`pca32f`), 128 ESM components (`pca128`), a pool-histogram mixture (`mix`), a triangle-update PairNet on the same inputs (`pairnet`) as prior levers; a smaller language model (`esm8m`) or none (`noesm`) as alternatives | null-to-worse on the built chain (`wide` +0.032, `pca32f` +0.018, `pca128` +0.076, `pairnet` +0.041, all under 0.5x MDE; `mix` chooses lam = 0 on 5/5 folds); removing or shrinking the ESM channel is worse (`noesm` +0.208, `esm8m` +0.242, both 5/5 folds, Type-M zone) | S26 L62, L63, L66, L67, L72, L93, L99, L103; Part VII.3 |

### VIII.2 Open

| item | what would close it | where |
|---|---|---|
| Whether a better distance predictor is obtainable (the only steep lever, -2.15 A per unit toward truth); every input this machine can compute is measured flat (nine rungs) | a larger language model than this machine can hold, on a bigger machine; Proposal C's kept form (L117, L122) | `s24/LEDGER.md` L13; `s26/PROPOSAL_C.md`; Part VII.3, VII.6 |
| The 2/60 benchmark self-copy leak, now bounded MINOR (dev proxy 0.008 A with both channels measured on all four dev self-copies; own-native envelope 0.028 A mean CI to 0.151 A worst target on the built chain; 0.194 at the worst target on the selection basis); F3's control clause open (one of six control-out models) | by design only a fresh benchmark, which does not exist; the five control-out models and a second seed of the 2P5H retrain, about two hours on this box (L137) | S26 L44, L49, L50, L55, L58, L108, L137; Part VII.4 |
| Where the target-specific third of the pool's coherent error comes from, and whether any native-free proxy is strong enough to act on it | a native-free proxy reaching the in-band ordering 2 A needs | `s19/LEDGER.md` L11, L14; `s17/LEDGER.md` L23 |
| Publishing the trainability half | a manuscript from Part V.10 with V.9's scope correction | `s13/`, `s25/QUANTUM.md`; S26 L27 |
| Not run at the close, recorded as such: the `raw` rung (four of five folds trained), `coherence_penalised_training` and `better_prior_inputs attn` (deferred on memory), a four-seed A1, lane W's Part B controls, `project_stability` beyond L141, the C5 alpha grid beyond 0.5 | their pre-registrations' falsifiers, on a machine with the headroom | `s26/TOURNAMENT.md`; `s26/PREREG_*.md`; S26 L133, L137, L139, L142; Part VII.3, VII.4 |
| The tie-break noise floor: measured (0.004 A on the 126-mean, 0.024 A paired MDE between conventions, built chain); not a lever | nothing; it is the floor every cross-run hundredths-level claim is read against | S26 L64, L71; Part VII.4 |
| A native-free difficulty flag: the pool's own disagreement predicts the emitted chain's error (partial rho +0.452 given n and Rg, 5/5 folds); the relaxation's displacement is its proxy; forbidden as a lever by the pre-registration and by the record (every conversion of a difficulty signal into a selector failed held out) | a use as a confidence label only | S26 L53, L81, L121, L123; Part VII.4 |
| Sequence proximity to the training corpus below the 0.6 threshold: no gradient at the pre-registered MDE (rho -0.25), a weak effect of order -0.2 suggested and confounded with retrieval | a design that separates recall from retrieval; none registered | S26 L90; Part VII.4 |

## PART IX. OPERATING MANUAL

Sources: `README.md` (layout, running, setup, the governor, the traps), `s26/TEST_RUN.md`,
`s26/EXAMINATION.md` sections G and H, `s26/examine.py`, `s25/resultslab/build.py`. Commands
are given for the Windows box the record was made on; the Python is
`C:/Users/abena/miniforge_3/python.exe` and the working directory is the repository root.

### IX.1 Environment

- `pip install -e .` installs the nine runtime dependencies named in `pyproject.toml` (numpy,
  scipy, torch, openmm, pennylane, pennylane-lightning, fair-esm, scikit-learn, biopython);
  `pip install -e ".[test]"` adds pytest. `lightning.qubit` is the device `core.quantum` builds
  circuits for and `tests/test_quantum.py` asserts probabilities against it with `==`, so it is
  not interchangeable.
- The machine: 16.75 GB RAM and 8 cores (4 fast and 4 slow, about 6.43 core-equivalents,
  `docs/STATE_BRIEF_2026-09-12.md` section 8); the baseline load with no campaign work is about
  67% RAM; the working rules are never above 93% RAM or CPU, never more than two OpenMM jobs at
  once (`README.md`, "The S26 resource governor").
- Every `.py` on the production path (`core/`, the 24 root modules, `s5/ s7/ s8/ s9/`) is the
  reference arm for an equivalence claim; edit none of them without re-running
  `tests/test_equivalence.py` and `tests/test_pipeline.py`. The full pre-consolidation tree is
  at commit `5fa05cd`.

### IX.2 Data and caches: which are pinned, which are rebuildable

Pinned, never delete, move or regenerate (`s26/results/pinned_hashes.json` holds their sha256;
`python s26/e_hashes.py --check` reports drift):

| artefact | what it pins | why regeneration is dangerous |
|---|---|---|
| `peptide_folds.json`, `peptide_clusters.json` | the five folds and 470 clusters | write-on-first-use: a missing file is silently re-derived, which once moved 13 benchmark targets and invalidated every model |
| `catrace_prior.npz` | the torsion prior | same write-on-first-use behaviour |
| `results/benchmark_manifest.json` | the sealed 60 | hashed as bytes; never parsed by any S26 script |
| `distogram_models/fold<f>_esm_frag.pt` | the five fold models | valid only for the fold definition committed alongside them; `distogram_models_large.STALE-PRE-FOLD-REPIN-DO-NOT-USE/` is quarantined by rename so that `FRAG_LARGE=1` fails loudly |
| `pdbs/` (61 files) | the reporting set | globbed and sorted by three modules; adding, removing or renaming a file re-orders BLOSUM tie sets and moves up to 47 of 500 pool members |

Rebuildable, ignored by git, owned by the module that rebuilds it (`.gitignore` names each
owner): `prots/` (5.5 GB) and `pdbs_ext/` (555 MB) from the RCSB fetch scripts kept in git
history; `esm_cache.npz` (1.5 GB, never load without a probe; the hot subset `esm_small.npz`
serves the 126 targets) from `esm_features.compute`; `peptide_db.npz` and `fragment_db*.npz`
from their `build` entry points; `*_models/` from `train_fold`. Caches that are the only copy of
a scientific record and must be preserved: `bench_results/cache/1fc9f2dcf489e2fb/` (the
production run's 126 checkpoints, gitignored), `s8/generate_univ/` (the 126 window universes),
`s24/cache_amber/*.npz` (63k AMBER single points; the state brief's outstanding item is to
whitelist and commit it), `s13/cache/` and the sprint `results/` directories. The full
not-in-git census is `s26/EXAMINATION.md` section G.

### IX.3 Running the production path on one target and reproducing its RMSD

The record for a target is `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`; it carries the
emitted `ca`, `avg_ca`, `fit_ca`, `amber_ca`, the shortlist `sub` and the labels. To reproduce
it from scratch, without AMBER:

    python s26/e_trace.py 1S9Z            # every stage, shapes and values, to s26/results/e_trace_1S9Z.json
    python s26/e_trace.py 9KAR --amber    # also the relaxation (tag AMBER; run through jobrun)

In code the trace does what production does:

    from core import pipeline as pl
    target, fold = ...                                   # from core.backend("data").load() and the pinned folds
    cfg = pl.Config(amber=False)                         # PROD otherwise; reference_precision=True
    rec, pool, clk = pl.run_target(target, fold, cfg, clk)
    lab = pl.label(rec, pool, target, clk, cfg)          # the ONLY function that opens the native

`lab` holds `rmsd_avg`, `rmsd_fit`, `rmsd_arm`, `shipped`, `pool_best`, `top_m_best`; on 1S9Z
`rmsd_arm` is 0.18198112330908295 and every value equals the stored record at 0.0
(`s26/EXAMINATION.md` H). To re-score every stored structure of all 126 records through the
instrument in five seconds, `python s26/e_reproduce.py` (writes `s26/results/e_reproduce.json`;
the four bases reproduce at 0.0). The instrument itself is `s12.instrument.ca_rmsd`; the natives
go through `core.pipeline._q` because that is what `label()` scored.

### IX.4 The full run, resuming, and the config key

    python -m core.pipeline run --manifest tuning126 --workers 6
    python -m core.pipeline stats --manifest tuning126

Manifests: `smoke8`, `smoke24`, `tuning126`, `dev24`, `benchmark60`. Every target writes one
atomic checkpoint under `bench_results/cache/<config-key>/<pdb>.json`; the key is a SHA-1 over
every parameter that can change a number (K, M, penalty, lam, the AMBER schedule, the fold
count, the pair separation, the optimiser cap, the tie-break rule, the live backend set), so
resuming is re-running the same command, a run with a changed parameter cannot inherit the old
answer, and `stats` reports how many targets were resumed. `core.bench --arm baseline|optimised
--fresh` measures both arms cold, `--compare` prints the speedup (the recorded cold runs:
2537.568 s at 1 worker against 303.137 s at 8, 8.371x end to end, with every science key
identical to 4.4e-16, `bench_results/compare_tuning126.json`), `--components` runs the
four-component system. `benchmark60` is refused without `--i-am-spending-the-benchmark`; do not
pass it (`README.md`, "The instrument discipline"). `CORE_BACKENDS=legacy` forces the reference
modules (`amber_refine, peptide_db, energy_terms, protein_geometry, legacy_field, s7.audit,
distogram, foldvqe, s8.project`) for an equivalence run.

### IX.5 The test suite under the governor

Never run `pytest tests/` bare on this box; the AMBER files can push it past 92% RAM. The suite
runs in three governed jobs (`s26/TEST_RUN.md` records each command verbatim):

    python s26/jobrun.py --agent I --tag TEST  --name pytest_core        --est-ram 2.0 -- python -m pytest tests/ -q -rs -p no:cacheprovider --deselect tests/test_amber.py --deselect tests/test_amber_frame_invariance.py --ignore=tests/test_amber.py --ignore=tests/test_amber_frame_invariance.py --junitxml=s26/results/pytest_core.xml
    python s26/jobrun.py --agent I --tag AMBER --name pytest_amber       --est-ram 2.0 -- python -m pytest tests/test_amber.py -q -rs -p no:cacheprovider --junitxml=s26/results/pytest_amber.xml
    python s26/jobrun.py --agent I --tag AMBER --name pytest_amber_frame --est-ram 1.5 -- python -m pytest tests/test_amber_frame_invariance.py -q -rs -p no:cacheprovider --junitxml=s26/results/pytest_amber_frame.xml
    python s26/i_test_report.py     # renders s26/TEST_RUN.md and s26/results/test_run.json from the junit files

The last governed run (branch `s26`, commit `601a39c7`): 370 tests, 357 passed, 0 failed, 0
errors, 13 skipped, none from the memory guard; the 13 skips are 11 `VERIFY_SLOW=1` opt-ins
(3 in `test_equivalence.py`, 8 in `test_integration.py`) and 2 absent artefacts; per-file
counts, wall times and peak RSS (1.692, 0.872, 0.324 GB) are in `s26/TEST_RUN.md`. Set
`VERIFY_SLOW=1` to run the full pipeline arms and the OpenMM checks, and set it inside the job's
command, not in the launching shell (a relaunch dropped it once, S26 ledger L98). The three gated
equivalence tests could never run before S26: their summary parser took the last brace of an
indented JSON dump and turned the failure into a skip; the test file was fixed at `7be8e4b0`
(L102) and the slow tiers re-run under the governor: the equivalence tier passed 3 of 3 from the
on-disk caches (L104) and the integration tier 8 of 8 (165.6 s, peak RSS 1.139 GB, L113), so the
record with the whole opt-in tier in is 370 tests, 368 passed, 0 failed, 0 errors, 2 skipped, the
two skips the absent artefacts (`s26/TEST_RUN.md`, `s26/results/test_run.json`, commit
`7ad4ef68`). Both AMBER files carry an
autouse memory guard that skips, with a message naming the ceiling and the governor's last
reading, when the box is above `core.amber`'s 92% ceiling.

### IX.6 The results lab

`s25/resultslab/` is one command:

    python -m s25.resultslab.build --mode selftest --limit 12     # sandbox smoke test, test-fixture labels only
    python -m s25.resultslab.build --mode frozen --spec s25/results/frozen_configs.json

`--mode frozen` is the only mode that writes to `results/`; the spec file, not the module,
decides what a configuration is, and each configuration's path is a JSON of emitted point
clouds keyed by PDB id, which the lab projects to the built chain itself. Four gates run on
every exported set (pool-oracle, difficulty-correlation, per-target sd floor, provenance REMARK
in every PDB header; the build refuses files without one), `providers.register` refuses
synthetic bindings under real configuration names, and every exported structure reproduces its
own RMSD through the instrument to within PDB quantisation (`s25/LEDGER.md` L10, L11). Outputs:
`results/summary/leaderboard.json` (rows with `mean` on `built_chain_bb` and `mean_secondary`
on the point cloud), `results/summary/target_map.json` (T001 to T126), `results/site/`. The S26
frozen rebuild (governed job `resultslab_rebuild`, L127) reproduced every number: 2016 of 2016
per-target RMSDs exactly, every leaderboard statistic and verdict for all eight configurations,
and all 2,142 PDBs in their ATOM records; the three descriptive leaderboard columns `selector`,
`hamiltonians` and `distogram_used` are not produced by any committed `schema.py` and rebuild as
null (the per-record values are intact in `results.json` and the spec), so a rebuild's overview
table shows them empty until `schema.leaderboard()` copies them from the records; not
hand-patched.

### IX.7 Adding a sprint

1. Make `s27/` with a `BRIEF.md`, a `LEDGER.md` and a `STATUS.md`; nothing in `core/` imports
   a sprint directory, and no sprint script edits the production path.
2. Write `s27/PREREG_<name>.md` before any result exists: hypothesis, the exact falsifier, the
   comparison arm, the basis of every RMSD, the expected effect against an MDE computed from
   persisted per-target arrays (`s26/PREREG_C2.md` is the template), memory and agent-hours.
   Never edit it afterwards; append addenda.
3. Run heavy work through `s26/jobrun.py` (or `s26/enqueue.py` when the box is busy) with the
   right tag; quote the peak RSS from `s26/jobs_done/<name>.json`.
4. Write artefacts with `s24.stats_lib.save_atomic` (module hash, git commit, dirty flag,
   `complete` gated on the expected rows) and compare with `s24.stats_lib.compare` (paired
   effect, SE, MDE, iid and fold-clustered CI, W/L, per-fold, concentration null); average tied
   argmins with `argmin_tied`; report grid minima through `best_of_k_within` and
   `split_half_transfer`.
5. Append to the ledger: re-read the tail immediately before appending, number the entry one
   past the last, suffix (`L30b`) on a collision, never edit an earlier entry; two STATUS lines
   per hour under the lane's heading.
6. Register every number the sprint will quote in `s26/results/claims.json` (path, key, value)
   so that `python s26/examine.py` re-reads it; `--search` runs the tree-wide search that finds
   where each number is cited and stored.
7. A positive result is re-run at a second seed and with the fold order reversed and must land
   inside its own CI before it is called a result; a negative one is stated with its MDE.

### IX.8 Running the governor

    python s26/governor.py                    # foreground; samples every 5 s; Ctrl-C stops it
    python s26/governor.py --once             # one sample
    python s26/governor.py --status           # the last snapshot (also s26/governor_state.json)
    python s26/jobrun.py --agent E --tag CPU --name probe --est-ram 0.5 -- python some_script.py
    python s26/enqueue.py --agent P --tag AMBER --name relax --priority 20 --est-ram 2.5 -- python s26/p_relax.py

The band: launch queued work when the box has sat under 88% for 60 s; resume suspended work
under 90%; suspend the newest registered job above 93% RAM or CPU; above 95% for 15 s, kill the
newest (CTRL_BREAK first, 25 s grace) and requeue it; at most two AMBER-tagged jobs at once.
Every action is appended to `s26/governor.log`. Read `s26/governor_state.json` before launching
anything heavy.

### IX.9 The examination

    python s26/examine.py             # module map, pinned-hash check, claim ledger (OK / MISMATCH / ABSENT)
    python s26/examine.py --search    # plus the tree-wide claim search (slower)
    python s26/e_module_map.py        # s26/results/module_map.json (787 modules on the final tree; the count grows with the sprint's scripts)
    python s26/e_hashes.py            # s26/results/pinned_hashes.json (--check to compare)
    python s26/e_claims.py            # s26/results/claim_search.{json,txt}; benchmark files excluded

`examine.sh` and `examine.bat` at the root call `examine.py`; there is no Makefile. Exit status
is non-zero if any step reports a problem. The standalone audits under `verify/` re-run one at a
time as `python s26/i_verify_rerun.py run <name>` (writes redirected to `s26/results/verify/`, the
tracked JSON compared leaf by leaf; table `s26/results/verify/REPORT.md`); at the close 20 of 20
had run with the science reproducing in every audit that has a comparable tracked record (L141;
Part VII.5).

## APPENDIX A. GLOSSARY

Terms a reader who knows some quantum computing and nothing about proteins will meet, in the
order a newcomer meets them. Where a term is a name in the code it is given in backticks.

**Protein, peptide, residue.** A protein is a chain of amino acids; each amino acid in the chain
is a residue. A peptide is a short chain; here 9 to 16 residues.

**Backbone, side chain, CA, CB, N, C, O.** Every residue contributes three backbone atoms, N
(nitrogen), CA (the alpha carbon) and C (the carbonyl carbon, with its oxygen O), to the chain;
what distinguishes one amino acid from another is the side chain attached to CA, whose first
atom is CB. Glycine has no CB.

**CA trace, virtual bond.** The sequence of CA positions, one per residue; consecutive CAs are
about 3.80 A apart in a real chain (3.8040 in the ideal geometry used here). That distance is
the virtual bond. It is 2.9 A across a cis peptide bond.

**phi, psi, omega; trans and cis; Ramachandran.** The three backbone torsion angles per residue.
Omega, across the peptide bond, is nearly always 180 degrees (trans); cis (0 degrees) is rare.
Phi and psi are the free ones, and the joint distribution real residues occupy is the
Ramachandran distribution; a Ramachandran penalty (`ramah`) pushes a chain toward it.

**Rotamer.** One of the discrete preferred conformations of a side chain. The builder here uses
one fixed rotamer per residue type, without scanning.

**Alpha-helix, beta-strand, polyproline II, extended.** Standard backbone conformations; a
constant alpha-helix is phi = -63, psi = -42 degrees at every residue and is the programme's
zero-information control.

**Radius of gyration (Rg).** The root-mean-square distance of the atoms from their centroid; a
compactness measure.

**PDB, deposit, model 1, NMR.** The Protein Data Bank holds deposited experimental structures,
each with a four-character id (1S9Z, 9KAR). NMR structures come as an ensemble of models; the
instrument scores model 1.

**Native, target.** The deposited structure of the peptide being predicted (the native) and the
peptide itself (the target). The native is a reporting label only; no deployable operator reads
it.

**Superposition, Kabsch, CA-RMSD.** To compare two CA traces, one is rotated and translated
onto the other to minimise the root-mean-square distance between corresponding atoms (the
Kabsch algorithm; reflections are not allowed). The minimised value is the CA-RMSD. The endpoint
of the programme is its mean over targets.

**Basis (of an RMSD).** Which object is scored: SINGLE WINDOW, POINT CLOUD (`rmsd_avg`), BUILT
CHAIN (`rmsd_arm`), RELAXED CHAIN (`rmsd_full`), or the lam = 0 chain (`rmsd_fit`). Part II.2.

**Library, fragment, window, pool, universe.** The library is 787 peptide chains plus 6,003
fragments cut from larger proteins. A window is a stretch of n consecutive residues of a
library member with its CA coordinates. The universe is every window of the target's length in
the target's fold; the pool is the top K = 500 windows by BLOSUM similarity.

**BLOSUM62, similarity, identity.** BLOSUM62 is a standard substitution matrix scoring how
alike two amino acids are; a window's similarity to the target is the sum over aligned
positions. Identity is the fraction of identical residues, normalised here by the longer
sequence.

**Cluster, fold (of the cross-validation), leakage, self-copy.** Library sequences are grouped
into identity clusters (470) and clusters into five folds; a target's fold model is trained on
the other four. Leakage is a training sequence that is the target; a self-copy is a target
sitting verbatim inside a longer training sequence that passed the identity threshold.

**Distogram, prior, posterior, bin, over-confident.** The distogram is a neural network that
predicts, for each residue pair, a probability distribution over 17 distance bins from the
sequence; "prior" and "posterior" both refer to that distribution (prior to seeing any
structure; posterior to the sequence). Over-confident means its stated uncertainty is smaller
than its error.

**ESM-2, `esm_pca32`, one-hot.** ESM-2 is a protein language model that maps a sequence to
per-residue vectors; the pipeline uses a 32-component PCA of them. One-hot is the trivial
encoding of residue identity as a 20-way indicator.

**Score, Bayes risk, L1 risk, `sc`.** The score of a candidate window is the sum over pairs of
the expected absolute distance error under the distogram, weighted by 1/(sd + 0.5): a Bayes
risk under the L1 loss. Lower is better. `sc` is the vector of scores over the pool.

**Argmin, shipped selector, top-m, shortlist, in-band.** The argmin is the single lowest-score
candidate (the S8 baseline). Top-m is the m lowest; the shortlist. In-band means among
candidates already close to the native, where ranking skill would have to act.

**Ranker, selector, readout, terminal operator.** A ranker orders candidates; a selector chooses
a subset; the readout (terminal operator) turns the chosen set into one structure. Production's
readout is the uniform coordinate average of the top 75.

**Medoid, consensus, typicality.** The medoid of a set is the member with the smallest mean RMSD
to the others; consensus selection picks it; typicality is a candidate's closeness to the rest
of the pool.

**Coordinate average, point cloud, contraction.** Superposing the retained windows onto the
medoid and averaging their coordinates gives a point cloud, which is contracted (its virtual
bonds are short) because a mean of scattered points lies inside them.

**Projection, built chain, `lam`, multi-start.** Fitting an ideal-geometry chain (fixed bond
lengths and angles, omega = 180) to the point cloud by choosing phi and psi; lam is the weight
of the Ramachandran penalty; multi-start runs the fit from several generic starts because the
problem is degenerate.

**Common mode.** The part of the pool's error that every member shares; 68% of the squared error.
Averaging cannot remove it.

**AMBER ff14SB, GBn2, OpenMM, kcal/mol.** A molecular force field (the set of functions and
parameters giving a molecule's potential energy) and an implicit-solvent model, evaluated by the
OpenMM library. Energies are in kilocalories per mole; a relaxed peptide sits around -1000, a
clashing one at 1e4 and far above.

**Lennard-Jones, steric clash, singularity.** The non-bonded term whose r^-12 wall makes two
overlapping atoms cost an enormous energy; a clash is such an overlap.

**Restraint, minimisation, relaxation, convergence gate, strain.** A restraint is a harmonic
penalty holding named atoms near their starting positions (k = 10 kcal/mol/A^2 on N, CA, C);
minimisation follows the energy downhill; the relaxation is the production minimisation;
`CONVERGE_MAX_KCAL` = 1000 is the energy below which a relaxation counts as converged; strain is
the residual bond-plus-angle energy.

**Legacy energy.** The programme's own eleven-term coarse score over backbone and CB atoms with
hand-set weights (`core/energy.py`).

**ORACLE.** Any arm that reads the native. A diagnostic ceiling, never a result.

**Zero-information control, matched control, best-of-N.** A control that knows nothing about
the target but is a plausible structure (the constant helix); a control drawn from the same
space as the operator with the same magnitude or count; the best of N draws from an untrained
optimiser, the control an optimiser must beat.

**Paired comparison, SE, CI, bootstrap, fold-clustered, W/L, concentration, MDE, Type-M.** Part
II.3. MDE is the minimum detectable effect, 2.8016 times the paired SE; a Type-M error is an
effect whose magnitude is inflated by having been selected for significance.

**Pre-registration, falsifier, artefact, `save_atomic`, complete.** A written plan with a
stated condition that would refute the hypothesis, filed before the result exists; the JSON
file a result is read from; the writer that stamps provenance and marks the file complete only
when every expected row is present.

**tuning126, dev24, benchmark60, FAIL18.** The 126-target development instrument; a 24-target
development split; the 60-target sealed benchmark; the 18 hardest development targets.

**Config key, manifest, cache.** A hash over every parameter that can change a number, naming the
directory of per-target checkpoints; a named list of targets; the checkpoint directory.

**Hamiltonian, diagonal, spectrum, `zrank`.** Here the Hamiltonian is a diagonal 128 x 128
matrix whose entries are the standardised ranks of the 128 shortlisted scores; its spectrum is
that ladder; `zrank` is the rank-then-standardise map.

**Qubit register, basis state, candidate index.** Seven qubits have 128 basis states; each
indexes one shortlisted candidate.

**Ansatz, RY, CNOT, chain, ring, layer, parameter.** The parameterised circuit; a single-qubit
rotation about Y; the two-qubit controlled-NOT; a CNOT on each neighbouring pair in order; the
closing CNOT from the last qubit to the first; one round of rotations plus entangler; one
rotation angle.

**Statevector, exact simulation, matrix product state (MPS), bond dimension chi, Schmidt
rank.** The full vector of 2^n amplitudes; computing it without approximation; a tensor
factorisation of the state whose internal index size is chi; the number of non-zero Schmidt
coefficients across a cut, which chi bounds.

**Parameter shift.** The exact rule giving the derivative of an expectation with respect to a
rotation angle from two evaluations shifted by plus and minus pi/2.

**CVaR, alpha, quantile, tail.** The conditional value at risk: the mean of the lowest-alpha
fraction of the energy distribution; the quantile is the cut; the tail is the set of states
below it.

**Free energy, temperature T, entropy, Gibbs (Boltzmann) distribution.** F = CVaR - T H(p); T
weights the entropy H; at alpha = 1 the minimiser over all distributions is
p(x) proportional to exp(-E(x)/T), the Gibbs distribution.

**KL divergence, total variation.** Two measures of the difference between distributions; KL in
nats or bits; TV as the largest probability mass on which they disagree.

**Adam, L-BFGS-B, SPSA, greedy, annealing.** Optimisers: a first-order stochastic method with
momentum; a quasi-Newton method with bounds; a two-evaluation stochastic gradient
approximation; a best-improvement local search; a temperature-scheduled random search.

**Barren plateau, 2-design, gradient variance.** A regime in which gradients vanish
exponentially in the qubit count; a circuit family matching the Haar distribution up to second
moments, for which that regime is proven; the quantity measured to test for it.

**Dynamical Lie algebra (DLA), so(2^n), su(2^n), controllable.** The Lie algebra generated by the
circuit's generators under nested commutators; the orthogonal and unitary algebras on 2^n
dimensions; an ansatz whose DLA is the full algebra can reach any state in the corresponding
group.

**Pauli string, Pauli weight, Walsh (Fourier) spectrum, Sobol share.** A tensor product of I, X,
Y, Z over the qubits; the number of non-identity factors; the decomposition of a diagonal
function over Z-strings; the fraction of variance carried by interactions of a given order.

**Quantum Fisher information, metric, QNG.** The Fubini-Study metric on the ansatz manifold; the
optimiser that preconditions by it.

**ADAPT, operator pool, Tang pools V and G, L2.** A procedure that grows an ansatz one operator
at a time from a pool; the minimal 2n - 2 element pools of Tang et al. 2021; the pool of all
1- and 2-local odd-Y strings.

**Torsion space, latent, lattice, encoding, k bits per residue.** Representing a chain by its
torsions rather than coordinates; a discrete code for them; the discrete grid of allowed values;
the map from bits to torsions; its resolution.

**Locality theorem.** Part V.10: the CA-CA distance d_ij depends on exactly the residues
strictly between i and j.

**Governor, `jobrun`, peak RSS, AMBER tag.** The S26 resource scheduler; its job wrapper; the
largest resident memory a job reached; the label for jobs that load OpenMM (at most two at
once).

**Results lab.** `s25/resultslab/`: the build that exports every structure set with provenance
and checks it through four gates.

**Ledger, STATUS, lane, sprint, coordinator, Adversary.** The append-only record of a sprint's
findings; the hourly status file; one agent's assignment; one campaign of work; the agent that
plans and rules; the agent that audits.

## APPENDIX B. EVERY NUMBER IN THIS REPORT, WITH ITS ARTEFACT

Rows are added as parts land. A number whose only source is a document is not in the report;
of the four document-only numbers of `s26/EXAMINATION.md` section C, C26 (36.1/36.4 deg) was
re-derived in S26 and now has a row; C27 (+0.0004/+0.0030; its artefact is in git history, L31)
and C34 (the |z_moment| triple 0.7529/0.8013/0.1127) are mentioned only as absent; C35 (355/13)
is superseded by the 370/357/13 of `s26/TEST_RUN.md`; and the 0.524
sampled tail-only cosine inside C24 is replaced by the three instrument-specific values of
`s25/QUANTUM.md` section 4.2. "as stored" means the leaf is the number to the precision printed;
"as cited" means the number is quoted from the named document section, which names its own
artefact; "as asserted" means a passing test pins it. Line numbers into `docs/FINDINGS.md` are
those of commit `73d82db5` (the S26 corrections block, L130, moved every body line by 24).

| number | where used | artefact (file :: key, or file:line) | stored value |
|---|---|---|---|
| 36.133, 36.416 | I | `s26/results/a_c26_phi_mae.json :: summary/arms/p_grid/pooled_mae_phi_deg, summary/arms/n_marg/pooled_mae_phi_deg` (Adversary re-derivation from `s13/cache/tors_rows.npz`, L31; pooled MAE over 1507 residues, 126 targets) | 36.132758, 36.415619 |
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
| 2.6087 | III | `s17/results/inband.json :: rows[*]/cells/25/band_best` (mean over 126 rows; reproduced `s26/PREREG_C1.md`) | 2.608684 |
| 2.2261, 2.8334, -2.1496, 0.822 | III | `s24/results/priorladder.json :: rows[*]/MASS1.0, rows[*]/MASS0.1, rows[*]/MASS0.0` (means over 126 rows) | derived: (2.833382 - 3.048338)/0.1 = -2.1496; 3.048338 - 2.226080 = 0.8223 |
| 4.0648 | I, II, III | `s12/results/s14_ladder.json :: rows/L0_constant_helix/mean` | 4.064753929494389 |
| 2.9610, 2.9507 | I | `s9/final_report.json :: dist/full/mean, dist/shipped/mean`, asserted within 5e-4 by `tests/test_pipeline.py::test_the_committed_benchmark_report_still_holds_its_reference_numbers` (not opened by this lane) | as asserted |
| +0.0103 [-0.1596, +0.1803], 31W/29L | I | `docs/FINDINGS.md:4503-4631` (S9-10), naming `s9/final_report.json` | as cited |
| 2.9614, 3.8122, 22.3%, 0.649, 80 of 126, 3.803955, 1.8e-15, 3.867 | II, III | `s26/results/e_reproduce.json :: summary/virtual_bond, n_rows` | as stored |
| +0.1664, +0.0977, [+0.0049, +0.3379], 16 of 126 | II, III | `s26/results/e_reproduce.json :: summary/projection_gap_arm_minus_avg, n_rows` | as stored |
| 0.0342, 0.0958, 0.2051 | II | `s26/results/q_mde_reference.json` | as stored |
| 0.394%, 1.18% | III | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` (hamiltonian stage); `s25/results/q_gibbs.json :: results/spectrum_target_independence` | as stored |
| 0.6758 (68%) | III | `s23/results/errdecomp.json :: rows[*]/f_common` mean | 0.675770 |
| 0.220, -560, 125/126, 5.38, 4.86 | III | `s26/results/ph_c3_nativefree.json` (S26 ledger L24) | as stored |
| 5370, -1290.6, 0.19, 52.0; 6.0e12, 1262.4, 0.61, 1166.1; 1172.7 (2BP4) | III | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` (relax stage); `bench_results/cache/1fc9f2dcf489e2fb/2BP4.json` | as stored |
| 7193, 9814, -1..19, -6..32, 144, 420, 105, 91, 183, 17 centres, 760, 488, 483, 0.6926..7.2688, 2.2420..5.1791, 401 98 193 161 284 260 372 373 381 382, 44/449, 3.75, 1.97, 1.16, 8, 128, 5.05..21.07, 0.23..2.70, 5.62..23.44, 7.58, 5.34, 1.40 | III | `s26/results/e_trace_1S9Z.json`, `s26/results/e_trace_9KAR.json` | as stored |
| 3.3135, 3.3443 | III | `s8/integrate_vqe.json :: arms/vqe_LFO/sel, arms/medoid128/sel` (reproduced by `tests/test_pipeline.py::test_the_four_components_all_execute_and_reproduce_their_published_numbers`) | as stored |
| 0.288 | III | `docs/FINDINGS.md:1994-2061` (S7 finding 11; its per-target artefact, s7/repr_tune.json, is lost, S26 ledger L11) | as cited |
| -0.172 [-0.316, -0.027], 74W/47L | VI | `docs/FINDINGS.md:2925` (S8-8 table) | as cited |
| 8.3 to 23.3 | IV | `s15/LEDGER.md` row 0.8 (`s15/results/audit_amber_cost_sweep.json`) | as cited |
| 0.016 | III | `docs/FINDINGS.md:2403` (S8-6 heading) | as cited |
| +0.142 | III | `ARCHITECTURE.md` section 2.5 | as cited |
| +0.0207 | III, IV | `bench_results/baseline_tuning126.json :: science/rmsd_full/mean, science/rmsd_arm/mean` | derived: 3.2354598538973844 - 3.214765154210998 = 0.0207 |
| 1.386 | II | `docs/FINDINGS.md:2775-2784` (S8-8) | as cited |
| 470, 13, 47 of 500 | II, III | `tests/test_data.py:305`; `README.md` (fold repin; layout) | as asserted |
| 23, 5 | II | `docs/FINDINGS.md:60-142`; `s25/LEDGER.md` L5, L7, L9, L11, L15 | as cited |
| 4/126, 2/60 | II | `s24/LEDGER.md` L4; `s26/results/i_identity_audit.json` | as cited |
| a40581ad...422d, 8002 | II | `s26/results/pinned_hashes.json` | as stored |
| 1000 | III | `core/amber.py:1282`; `s8/integrate.py:292` | as in source |
| 2.239e-4, 8.66e-4 | II | `s25/LEDGER.md` L10 | as cited |
| 64.7th percentile | II | `s25/agentQ_FINDINGS.md:375-376` (`s25/QUANTUM.md` 6.4) | as cited |
| Legacy weights 4.0, 1.0, 1.0, 3.0, 2.0, 2.0, 0.5, 1.0, 0.8, 0.15, 0.4 | IV | `core/energy.py` `DEFAULT_WEIGHTS` | as in source |
| 0.3 ms, 9 ms, 6 to 12 s | IV | `s16/repair_FINDINGS.md` section 2.4 | as cited |
| 28 ms (corrects 6 ms) | IV | `docs/CONDENSED_REPORT.md:185` (corrections table) | as cited |
| 112/192, 186/192, 0.358, 0.827, 0.882 | IV | `s20/LEDGER.md` L6 | as cited |
| -0.758 (SE 0.023), +0.60, +1.103 (SE 0.042), -0.27, -0.090 (SE 0.018), +0.307, +0.279, -0.027, +0.047 | IV | `s25/results/phys_landscape.json :: summary` | as stored |
| -0.0886, -0.0829 | IV | `s20/results/c_q1.json`; S24 (via `s26/PH_PART_IV_NOTES.md` section 3) | as cited |
| 3.058, 3.132, 3.215, 3.253, 3.425, 3.674, 3.755, 3.881 (point cloud) | IV | `s25/results/phys_suite.json :: configs_rank, random_null_rank` | as stored |
| +0.330 (1.99x MDE), +0.455 (2.42x MDE), +0.327, +0.471 | IV | `s25/agentPHYS_FINDINGS.md` sections 1.4, 1.5 (from `s25/results/phys_suite.json`) | as cited |
| 3.2187, 3.3100, 3.3732, 3.4221, 3.8248, 3.8844, 4.1015 (built chain) | IV | `results/summary/leaderboard.json :: rows[*]/mean` | as stored |
| +0.068 [+0.012, +0.125] | IV | `s16/LEDGER.md` L17 (`s16/integrate.py`; interval flagged provisional there) | as cited |
| 15.3, 97.0%, 99.66% | IV | `s25/results/phys_landscape.json :: summary/AMB_decades/mean, summary/AMB_frac_absz_lt_0p1/mean, summary/AMB_top10_var_share/mean` | 15.2556, 0.96989, 0.99656 |
| 58.6% | IV | `s26/results/ph_reject_census.json :: summary/per_threshold/1e4/pool_frac_over/mean` | 0.58554 |
| 99.56% (60.7% to 99.98%), 46.5%, 26.6%, 25, 41 | IV | `s13/results/walsh_xval.json :: concentration[model=amber], concentration[model=legacy]` (medians over the 25 and 41 entries) | as stored |
| 1.000, 0.0, 3, 900, 7e30 to 2e31, 6 of 6 | IV | `s13/results/walsh_amber.json :: exact[*].terms`; `s13/walsh_FINDINGS.md` 3.1, 3.2 | as stored |
| 0.0745, 0.4221, 30W/0L, 18.8, 5.8, 0W/30L | IV | `s20/results/c_land_report.txt` section 1 | as stored |
| 0.577, 0.104, +0.444 of +0.620 (72%) | IV | `s20/results/c_land_null.json :: rows[*]/null/amber/theta_moved, rows[*]/null/legacy/theta_moved` (means over 30 rows); `s20/LEDGER.md` L12 | as stored |
| 9450, 40.7/75, 2.6/75, 5057, 96.8%, 2627, 2267, 163, -0.74 | IV | `s26/results/ph_reject_census.json :: summary/singularity, n_rows`; `s25/results/phys_landscape.json :: m_prod` | derived: 126*75 = 9450; 163+2627+2267 = 5057; (2627+2267)/5057 = 0.968 |
| 2.66 | IV | `s19/LEDGER.md` L12 | as cited |
| 40/126, 462/500, +0.054 | IV | `s25/agentPHYS_FINDINGS.md:57-58`, `s25/agentPHYS_FINDINGS.md:357`, `s25/agentPHYS_FINDINGS.md:379` (its section 3, computed from `s25/results/phys_suite.json :: normalisation_fork`) | as cited |
| 58.7%, 8.6e4, -560, 125/126, +1262, 0.220, 3.804 to 3.867, 5.38, 4.86 | IV | `s26/results/ph_c3_nativefree.json` | as stored |
| n = 7, 128, L = 3, P = 21, 21 RY, 21 CNOT, depth 21 / ~24, 50 steps, seed 0, 8128 | V | `s25/results/q_verify.json :: results` (deployed register, parameter count); `s25/QUANTUM.md:180-192` (gate counts and depth); `core/pipeline.py:181-184` (`Config`) | derived: 2^6*(2^7-1) = 8128 |
| 5.6e-17, 3.331e-16, 6.661e-16, {1: 2, 2: 4, 3: 8, 4: 16}, [0.879145, 0.461538, 0.105625, 0.054141, 0, 0, 0, 0], 0.334 | V | `s25/results/q_verify.json` | as stored |
| 4.06e-2, 3.4371, 1.18%, 0 of 8 | V | `s25/results/q_gibbs.json :: results/spectrum_target_independence` | as stored |
| 0.0761 bits | V | `s25/results/q_alpha.json` (alpha = 1, T = 0.1 cell) | as cited (`s25/QUANTUM.md` 3.4) |
| VQE_LFO table, 78 of 126, 0.6190 | V | `core/pipeline.py:113`; `s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint, results/n` | derived: 0.6190476190476191*126 = 78 |
| 1.000000000, 4.597e-10, 4.663e-10, 42 | V | `s25/results/q_verify.json :: results` | derived: 2*21 = 42 |
| +0.655634 / 0.758, +0.566586 / 0.519, +1.000000, +0.994 | V | `core/quantum.py:42` (S9 instrument); `s25/results/q_verify.json :: results` (S25 re-verification); `docs/FINDINGS.md:3184` (the sampled estimator, S8-9) | as stored |
| 2592, 1620, 0, 0, 1424 (54.9%), 972 / 972; 3888, 29.9%; 17574, 58.8%; 2016 / 2016 | V | `s25/results/q_verify.json`; S24 harness and coordinator runs as cited in `s25/QUANTUM.md` section 5; `s22/LEDGER.md` (Gate 1) | as stored / as cited |
| 3.4540, 3.3414, 3.2835, 3.3135 | V | `s8/integrate_vqe.json` (S8 instrument rungs) | as stored |
| alpha-effect table: -0.1126 (0.0792), -0.1067 (0.0644), +0.0279 (0.0518), -0.0011 (0.0477), +0.0126 (0.0268), +0.0039 (0.0240) | V | `s25/results/q_alpha.json` | as stored |
| -0.7423, +0.2700, -0.0234, 0.032, 0.7045, +0.0090 (sd 0.0342), 0.0268 | V | `s25/results/q_alpha.json` (entropy curve) | as stored |
| +0.0499 (0.0427), -0.0302 (0.0450), +0.0445 (0.0511) | V | `s25/results/q_alpha.json` (circuit minus Boltzmann) | as stored |
| Gibbs ladder at T = 0.3 (seven rows), 1.449 / 0.902 / 0.373, 0.761 / 0.453 / 0.351, -1.483528, -2.453671, 0.891 / 0.783 / 0.783, -1.7176 / -1.7186 / 0.075, 0.458 | V | `s25/results/q_gibbs.json :: results/*, results/training_control` | as stored |
| -0.1405 (0.0732, 66W/48L), -0.0002, +0.0262 (0.0279), -0.0308 (0.0590), -0.0081, 12 | V | `s25/results/q_alpha.json :: results/vs_no_circuit` | as stored |
| 45.3%, 64.7th | V | `s25/agentQ_FINDINGS.md:375-376` (drop-top-5 against a uniform-effect null; `s25/LEDGER.md` L5) | as cited |
| -0.0311 (0.27x MDE), -0.0021 | V | `s25/results/q_alpha.json` (forced-alpha counterfactual) | as stored |
| -0.013 [-0.095, +0.077] | V | `s20/LEDGER.md` L1 (CNOT deletion on the S20 continuous encoding; cited by `s25/QUANTUM.md` 7.4) | as cited |
| width-sweep table (seven widths, five columns), slopes -0.6492 / -0.2522 / -0.0472 / -0.3105 / -0.2429, depth sweep (seven depths), CVaR ratio tables (two columns of seven), draws 250 / 200 / 120 / 80 | V | `s25/results/q_plateau.json` | as stored |
| DLA table: 120, 496, 510, 1023, 2016, 8128, 32640, 32766, 65535, 130816, 523776, 2096128; pools 36, 136, 528, 2080, 8256, 32896; ADAPT 7, 16, 1025 (12.6%); 12 of 12; 85.3 s; 0.479 GB | V | `s26/results/q_dla.json :: results/fixed, results/pools, results/adapt_sets, results/numeric`; `s26/jobs_done/a2_dla.json` | derived: 1025/8128 = 0.1261 |
| 0.504 (and 0.823, 0.768, 0.648, 0.537) | V | `s13/results/geo_kernel.json` (via `s13/SPRINT13_DOSSIER.md` section 10) | as cited |
| 1.0000, 5,000+, 74, 52, 1.4e-13, 4.6e-2 | V | `s13/qarch_FINDINGS.md` section 1 (`s13/SPRINT13_DOSSIER.md` section 6) | as cited |
| 141, 99.6%, 6.001 / 6.001, 0.0003, 0.878, 0.986 | V | `s13/walsh_FINDINGS.md` (`s13/SPRINT13_DOSSIER.md` section 7; raw sweep `s13/results/geo_pauli_v1_rawonly.json`) | as cited |
| 2.236, 3.015, 79 / 79, 0.641, 0.392, 0.069, 0.122, 12 / 13, 0.778, 0.629, 1.006, 1.001, 1.278, 0.913, 95, 0.15 to 4.18 | V | `s13/SPRINT13_DOSSIER.md` section 8 (Pauli tables of `s13/results/`) | as cited |
| 1.000 on 26 cells, 0.955, 3.89 / 4.10 / 4.17 / 4.27, 0.60 to 0.64, 2.70 / 10.00 / 1.27, 7.9, 1e-16 | V | `s13/SPRINT13_DOSSIER.md` section 9 (`s13/results/walsh_amber.json`) | as cited |
| 1.00, 2^(-0.47n) to 2^(-0.86n), 0.93 to 1.00, 0.862 to 0.078 | V | `s13/SPRINT13_DOSSIER.md` section 10 (`s13/results/geo_kernel.json`) | as cited |
| 0.000e+00, 0.2500, 0.008 to 0.037, 1.7e16, 0.088, +0.333 | V | `s13/SPRINT13_DOSSIER.md` section 11 | as cited |
| +0.209 to +0.301, +0.27 to +0.62 | V | `s16/LEDGER.md` L17, L29 (`s16/integrate.py`; `s16/qphase_FINDINGS.md` section 3) | as cited |
| 20 of 20 | V | `s12/SPRINT12_DOSSIER.md` section XVIII | as cited |
| 2.406, 1.925 | VI | `docs/FINDINGS.md:2311` (S8-5); `docs/CONDENSED_REPORT.md:176` (correction) | as cited |
| 3.204 | VI | `docs/FINDINGS.md:3627` (S8-11) | as cited |
| 0.5 per pair, 0.044 | VI | `docs/FINDINGS.md:3834` (S8-14) | as cited |
| 2.7x | VI | `docs/FINDINGS.md:4750` (S10-2) | as cited |
| 1.62 | VI | `docs/FINDINGS.md:5243` (S11-2) | as cited |
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
| 16.75 GB, 8 cores, 6.43, 67%, 93%, 92%, 88 / 90 / 93 / 95%, 60 s, 15 s, 25 s, 5 s | IX | `README.md` ("The S26 resource governor"); `docs/STATE_BRIEF_2026-09-12.md` section 8; `s26/governor.py` | as in source |
| 5.5 GB, 555 MB, 1.5 GB, 61, 63k | IX | `README.md` ("Data this repository does not carry"); `docs/STATE_BRIEF_2026-09-12.md` section 8 | as cited |
| 0.18198112330908295 | IX | `s26/results/e_reproduce.json` (row 1S9Z, `rmsd_arm`); `s26/results/e_trace_1S9Z.json` | as stored |
| 1.692, 0.872, 0.324, 11, 3, 8, 2, 601a39c7 | IX | `s26/TEST_RUN.md`; `s26/results/test_run.json` | as stored |
| 787 | IX | `s26/results/module_map.json :: n_modules` (regenerated by lane I on the final tree, d0ec972c) | 787 |
| 2.9661, 0.1824, 8.2342, 28.6%, 50.8%, 0.7150, 0.1643 | II, D | `results/summary/leaderboard.json :: rows[0]/median, rows[0]/best, rows[0]/worst, rows[0]/frac_under_2, rows[0]/frac_under_3, rows[0]/corr_with_pool_best, rows[0]/basis_delta_mean` (pointed to by `docs/REPORT_S26.md` B.1) | as stored |
| +0.0061, 0.19x, 61W/65L | IV, D | `results/summary/leaderboard.json :: rows[1]/paired_effect, rows[1]/effect_over_mde, rows[1]/wins, rows[1]/losses` | as stored |
| +0.7070 (2.14x), +0.8322 (2.67x), +0.0097 (0.43x), +0.1668 (0.91x), +0.0833 (0.70x), +0.6259 (1.89x), +0.2047 (1.07x) | IV, D | `s25/results/phys_suite.json :: vs_incumbent` (point cloud against point cloud) | as stored |
| 1.9962, 0.6159, 0.2824, 0.2412 | III, D | `s25/results/calib.json :: rows[*]/z_sd, rows[*]/cov90, rows[*]/cov50, rows[*]/multimodal_frac` (means over 126 rows) | as stored |
| 2.8334, -0.2150, 2.4x, 113W/13L | III, D | `s24/LEDGER.md` L13 | as cited |
| 1.75, 51, +0.054 | III, D | `s25/LEDGER.md` L6, L12 | as cited |
| 5.551e-17 | V, D | `s25/results/q_verify.json :: results` | as stored |
| 0.1 nats, fired | V, D | `s25/results/q_gibbs.json :: results/falsifier` | as stored |
| -0.023 | V, D | `docs/FINDINGS.md:392-427` (S5 section 6) | as cited |
| 0.6847, 0.685 | V, D | `verify/cvar_audit.json :: D_tail_baseline_cosine_vs_paramshift, D_sampled_tail_cosine` | as stored |
| 8.371, 2537.568, 303.137 | IX, D | `bench_results/compare_tuning126.json :: end_to_end_speedup, wall_s` | as stored |
| 3.989, 3.213 | VI, D | `s12/SPRINT12_DOSSIER.md` section I | as cited |
| 13 of 60, +0.0030 | II, D | `docs/FINDINGS.md:4595` (S9-10); S10-4 (artefact in git history at `5fa05cd`, S26 ledger L31) | as cited |
| A4 slopes -0.649 / -0.079 / +0.006, -0.252, -0.047, -0.311 / +0.024 / -0.008, -0.243 / -0.239 / -0.302 | V, VII | `s26/results/q_var.json :: results/slopes` | as stored |
| +0.035, -0.246, 1.63, 2.5 to 370, 31 of 32, 1 to 3, 12 to 25, 4 to 21, 2.4e-2, 14 of 14, 3 of 7, 2,765 s, 0.08 GB | VII | `s26/LEDGER.md` L35; `s26/jobs_done/a4_var.json` | as cited |
| 6 of 6, 1025 / 8128 = 0.126 | V, VII | `s26/results/a_dla_check.json`; `s26/LEDGER.md` L45 | as cited |
| 8.76 GB, 7.4 GB, 4.4 GB, 64.8%, 16.75 GB, 2,771,653,574, 5,678,116,398, 8.45 GB, 2.60 GB, 30 MB | VII | `s26/results/b1_feasibility.json`; `s26/LEDGER.md` L13 | as cited |
| 3.7705, +0.5557, +0.2314, 0.1318, 0.3694, 37W/89L, +0.8500, 0.1494, 32W/94L, -0.2943, 0.0807, 69W/57L, 5.569, 6.032, 5.887 | VII | `s26/LEDGER.md` L14 (`s13/cache/tors_rows.npz`, `s14/results/ladder.json`) | as cited |
| 0.2198, 0.0084, 8.59e4, 58.73%, -559.8, 6 of 126 | VII | `s26/LEDGER.md` L48 (recomputation of L24 from `bench_results/cache/1fc9f2dcf489e2fb`) | as cited |
| C3 stage 1 table: +0.0207 [+0.0154, +0.0290] 2.16, +0.0111 [+0.0062, +0.0171] 1.12, +0.0385 [+0.0298, +0.0479] 2.27, +0.0096 [+0.0077, +0.0122] 2.28, -0.0178 [-0.0255, -0.0124] 1.49, -0.0503 [-0.0735, -0.0253] 1.17, -0.1719 [-0.2104, -0.1388] 2.65; 40/86, 52/74, 29/97, 90/36, 74/52, 94/32; -0.049 (0.015), 36.5%, +0.0107, -0.0521, 1.07x, 0.220, 16, 1e-6, 10 s, 0.041 GB, 124 of 126, 5.38, 4.86 | VII | `s26/results/ph_c3_stage1.json`; `s26/LEDGER.md` L39, L46; `s26/jobs_done/ph_c3_stage1.json` | as stored / as cited |
| steric reject table: +0.2276 [+0.1466, +0.3232] 1.17 49/67/10; +0.1672 [+0.0699, +0.2793] 1.26 47/77/2; +0.1053 [+0.0187, +0.2000] 0.98 52/66/8; +0.1079 [+0.0647, +0.1528] 1.34 38/75/13; +0.0709 [+0.0317, +0.1087] 1.00 41/72/13; +0.0605 [+0.0206, +0.0962] 0.57 62/62/2; +0.0371 [+0.0229, +0.0495] 1.47 35/78/13; +0.532, 54.5, +0.228, 40.1, +0.093, 30.3, 0.66x, +0.066, 23.3, 0.51x, +0.061, +0.167, +0.254, -0.3368, -0.1470, 3.74, 0.13 | VII | `s26/results/ph_reject_report.json`; `s26/LEDGER.md` L43 | as stored / as cited |
| 1,507, 1.91, 0.34, 5.8, 17.4, 42.8, 0.53%, 0.13%, 0 of 126, 1,966, 2,352,893, 3.5045, 439, 725 | VII | `s26/results/ph_cis_census.json`; `s26/LEDGER.md` L22, L48 | as stored / as cited |
| 0.347 (0.029), 0.272, 0.752, 1.474, 0.340, 32 of 126, +0.828, +0.083, -0.036, 0.166, -0.1804 [-0.2266, -0.1022], 86W/40L, 0.25 | VII | `s26/results/ph_cis_floor.json`; `s26/LEDGER.md` L38, L48 | as stored / as cited |
| 160 s, 0.302 GB, 1.4e-14, 3.4540, 3.0483, 3.2126, 3.2052; -0.0683, +0.0071, 0.0000, 122, 47, 0.078, 0.147, 0.511, 61st to 89th; +0.0004 [-0.0001, +0.0010], 0.0012, 117, +0.0018 [-0.0003, +0.0039], +0.0003; -0.6980 [-0.8312, -0.5764], 112/14, -0.7094, -1.2081 [-1.4377, -0.9745], -0.6932, 0.47, 97%, +0.5101; -0.0114, -0.0204, +0.0062; 18, 22, 2.908, 0.317, 5.502, 18%, 1.044; 0.0023, 0.0067, 0.0004, 0.0277, 0.0479, 0.0228, 0.028, 0.048, 0.023, 0.002, 0.170; 1.29, 0.74, 0.85, 2.11, 0.04, 0.09 | II, VII, VIII | `s26/results/w_selfcopy_endpoint.json`, `w_selfcopy_bound.json`, `w_selfcopy_floor.json`; `s26/LEDGER.md` L44, L50; `s26/jobs_done/w_endpoint_report.json` | as stored / as cited |
| 0.0683, 0.2016, 0.8312, 1.4377, 0.6838, 0.0830, 0.1512, 0.0280, 0.0824, 0.1587, 0.0816, 0.1232, 0.1154, 0.1944, 0.194, 0.151, 0.123, 0.19, 0.05 | II, VII, VIII | `s26/results/w_selfcopy_bound.json :: signed_bounds_gated`; `s26/LEDGER.md` L55, L58 | as stored / as cited |
| 416 s, 0.13 GB, 0.00e+00, 7.25e-14, 3.2126, 3.2052, 120 of 126, 0.0107, 0.1711, 7JS6, 0.0022 | VII | `s26/results/p_ladder_shipped_s0.json`; `s26/LEDGER.md` L56, L57; `s26/jobs_done/p_eval_shipped2.json` | as cited |
| +0.0027, +1.28, +4.15, 8T61, +0.532, +0.196, +0.108, +0.037, +0.061 | VII | `s26/LEDGER.md` L54 (from `s26/results/ph_reject_report.json`) | as cited |
| 2.908, 2.811, 0.317, 5.502, 4 of 22, 6 of 22, 3.055, 0.595, 3.278, 2.334, 4.126, +0.9677 [+0.7247, +1.2809], 3W/15L, 1.06x, -1.4428 [-1.7172, -1.1896], 13W/5L, 1.14x, 18, 22 | VII | `s26/results/w_selfcopy_floor.json`, `s26/results/w_identity_floor.json`; `s26/LEDGER.md` L52 | as stored / as cited |
| +0.433 [+0.247, +0.588], +0.241, +0.250, +0.133, 0.00025, 0.0125, [+0.248, +0.581], 0.159, 0.218, 0.288, 2.286, 2.936, 3.758, 3.923, +0.41, 0.113, 35 s, 0.1 GB, 32 | VII, VIII | `s26/results/ph_strain.json`, `s26/results/ph_strain_rep.json`; `s26/LEDGER.md` L53 | as stored / as cited |
| A1: 7,782 s, 0.383 GB, 25 s, 33, 126 of 126, 3.2280, 3.3135, -0.0138 (0.0210, 0.0588, 0.23x, [-0.0705, +0.0442], 58W/68L), -0.0222 (0.0217, 0.0608, 0.36x, [-0.0854, +0.0424], 61W/65L), 50th, 51st, -0.0437, -0.0448, 0.47x, 45 to 59, -0.0128 to -0.0245, 0.21x to 0.40x, -0.0035, +0.0088, +0.0135, +0.0230, +0.0337 [+0.0122, +0.0581], 0.34x, 0.059 to 0.061, +0.0133, 0.9027, 0.984, 0.0002, 0.0009, 1.4e-4, 7.9e-4, 0.90, -0.02, 78, 48 | V, VII, VIII | `s26/results/a1_stats.json`; `s26/results/a1/<pdb>.json`; `s26/jobs_done/a1_build.json`, `a1_label.json`; `s26/results/q_mde_reference.json`; `s26/LEDGER.md` L68 | as stored / as cited |
| L70: 0.955, 0.919, 96.0%, -0.094, +0.006, 7 to 21, -0.0271, +0.0671, +0.0047, 29 of 78, 0.206, 2.6e-4 | VII | `s26/LEDGER.md` L70 | as cited |
| L75: 60 of 78, 68 of 78, 78 of 78, 235, 393, 1092, 1.2e-4, 8.6e-4, 9.9e-5, 0.02, 0.018, 0.28, 0.11, 4.1e-4, 3.3e-4, 3.0e-4, 2.7e-4, 0.053, 0.126, 0.148, 1e-3 | VII | `s26/results/a1/<pdb>.json :: adapt`; `s26/LEDGER.md` L75 | as cited |
| C2 rung table: noesm +0.2078 (0.0733, 1.01x, [+0.0986, +0.3942], 47W/79L), +0.3299 (1.27x, [+0.2122, +0.4475]); conly +0.1223 (0.0708, 0.62x, [-0.0350, +0.2674], 53W/73L), +0.1116 (0.42x), -0.086, -0.218 (1.00x); wide +0.0317 (0.0458, 0.25x, [-0.0300, +0.1279], 54W/72L), +0.0972 (0.48x, [+0.0487, +0.1587]); pca32f +0.0180 (0.0491, 0.13x, [-0.0725, +0.0970], 70W/56L), -0.0316 (0.17x); pca128 +0.0759 (0.0630, 0.43x, [-0.0691, +0.2157], 57W/69L), +0.0249 (0.13x); pca32 +0.0000; +0.1161, +0.247, 42, 13, 768, 4, 2.7x, 128, 32 | VII, VIII | `s26/results/p_ladder_report_{noesm,conly,pca32,wide,pca32f,pca128}_s0.json`; `s26/LEDGER.md` L62, L63, L65, L66, L67, L72 | as stored / as cited |
| tie-break floor: 5,166 s, 0.111 GB, 25 s, 0.062 GB, 115, 54.5, 8, 3.57, 3, 15, 91.5%, 0.168, 1.08, 0.0039, 0.0017, 0.0032, 0.0136, 0.0236, 28, 0.0321, 0.0227, 0.129, 21 of 126, 0.713, 1D6X, 0.002, 0.010, 0.003 to 0.010, 0.02 to 0.05, 3.6, +0.0026 (0.15x), +0.0048 (0.55x), +0.0201 (0.75x), 107, 0.024, 0.004 to 0.005 | VII, VIII | `s26/results/w_tiebreak_report.json`, `w_selfcopy_tiebreak_draws.json`, `w_selfcopy_tiebreak_endpoint.json`; `s26/LEDGER.md` L64, L71 | as stored / as cited |
| deck: 11, 237, 296, 302, 4, 7, 248, 247, 244, 249 | VII | `s26/pr_values.json`; `s26/LEDGER.md` L61, L73, L76 | as cited |
| 91.9%, 180 s, 5.0 GB, 2.3 GB, 1.7 GB, 0.87 GB, 16%, 2.5 GB | VII | `s26/LEDGER.md` L74, L77 | as cited |
| R1 to R4 as listed | VII | `s26/RETRACTIONS.md` | as cited |
| mix: lam ladder 3.2126, 3.2286, 3.2295, 3.2555, 3.2871, 3.3451, 3.6556, 3.9455; +0.0385, 0.21x, -0.3988, 183%, -0.1647, 5.81, 4,332 s, 0.13 GB, +0.73 | VII | `s26/results/p_ladder_report_mix_s0.json`; `s26/LEDGER.md` L93 | as stored / as cited |
| L79: +0.044, +3.65, 0.48x, [-0.051, +0.240], +0.058, 0.37x, 0.18 | VII | `s26/LEDGER.md` L79 | as cited |
| C3 replication: -0.0207 [-0.0250, -0.0166], 1.77x, 94W/32L, +0.0100 [+0.0059, +0.0175], 1.02x, +0.0414 [+0.0331, +0.0501], +0.0107 [+0.0093, +0.0120], -0.0459 [-0.0734, -0.0212], 0.018 to 0.021 | VII, VIII | `s26/results/ph_c3_stage1_rep.json`; `s26/LEDGER.md` L87 | as stored / as cited |
| steric reject, built chain: 11,976 s, 3.3 h, 22.8, 0.09 GB, -0.0021, 0.171, 6, +0.2483 [+0.1197, +0.3479] 1.17x 49W/67L/10T, +0.1574 [+0.0480, +0.2811] 1.15x, +0.1021 0.88x, +0.1037 [+0.0373, +0.1662] 1.11x, +0.0551 0.73x, +0.0909 0.75x, +0.0486 0.99x, +0.5631 [+0.4413, +0.7037] 1.77x, +0.1921 1.44x, +0.1190 0.78x, +0.0833 0.95x, +0.0915 0.65x, +0.0511 0.76x, -0.3442, -0.1473, 3.86, -0.002, +0.267, +0.226, +0.421, +0.310, +0.002, +1.39, +4.20, 0.14 | VII, VIII | `s26/results/ph_reject_chain.json`, `s26/results/ph_reject_report.json`; `s26/LEDGER.md` L86 | as stored / as cited |
| branch select: 4.77, +0.0055 (0.0109, 0.0306, 0.18x, 39W/45L/42T), -0.1021 [-0.1270, -0.0750] 2.29x 91W/35L, +0.1076 2.65x, +0.0922 1.10x, -0.0036 0.09x, 3.2181, 3.3048, 3.3202, 3.1301, -0.19, 130%, 56%, 4.63, 30%, 20%, 41, 85, 0.3 GB, 65 s, 0.03 | VII, VIII | `s26/results/ph_branch_report.json`, `ph_branch_solutions.json`, `ph_branch_relax.json`; `s26/LEDGER.md` L88 | as stored / as cited |
| cis floor tight: 0.0833 (0.0075), 0.0596, 0.215, 0.435, 6MBM, 0.0429, -0.2636 [-0.2998, -0.2289], 113W/13L, 3.46x, +0.0404 [+0.0318, +0.0492], 23W/103L, +0.0832 [+0.0440, +0.1331], 1.51x, 2,813 s, 0.104 GB, 0.083, 0.043, 0.44, 0.04 | III, VII, VIII | `s26/results/ph_cis_floor2.json`; `s26/LEDGER.md` L89 | as stored / as cited |
| ensembling: 4,351 s, 0.113 GB, 1.8e-15, 10 s, 2.104, 2.572, 0.83, 0.89, 0.88, 93, 0.14 to 0.20, 0.185, 0.249, -0.0005 (0.0100, 0.0281, 0.02x, [-0.0194, +0.0162]), +0.0042 (0.14x), +0.0016, -0.0008, -0.0039, +0.0237 (0.70x), -0.0153 (0.40x), -0.0116 (0.52x), -0.0021, -0.0047, 0.028, 0.019 | VII, VIII | `s26/results/w_selfcopy_ensemble_clouds.json`, `w_selfcopy_ensemble_endpoint.json`; `s26/LEDGER.md` L84 | as stored / as cited |
| provenance: 0.6%, 0.7%, 30 of 126, 8.0%, 7.2%, 18.7%, 18.3%, 72.7%, 73.8%, -0.4161 [-0.4472, -0.3758] 3.18x, -0.421, -0.411, -0.011 0.13x, -0.086 0.96x, -0.1313 [-0.1646, -0.1096] 1.13x, 114, 0.42, 0.09, 0.13, 0.05 | VII, VIII | `s26/results/w_selfcopy_provenance_census.json`, `w_selfcopy_provenance_oracle.json`; `s26/LEDGER.md` L85 | as stored / as cited |
| recall gradient: 0.470, 0.333 to 0.588, 0.657, 4.4, 3 to 13, -0.205 [-0.315, -0.099], 0.81x, 0.253, 0.010, 0.017, -0.157, -0.088, -0.230, 0.91x, -0.25, -0.2, 15 s, 0.292 GB, 45 s | VII, VIII | `s26/results/w_selfcopy_recall_covariates.json`, `w_selfcopy_recall_endpoint.json`; `s26/LEDGER.md` L90 | as stored / as cited |
| memorisation ladder: 504, 28.1%, 0.584, 0.473, 77.3%, 0.911, 2.339, 0.825, -0.709, -0.698, -0.603, -0.364, +0.346, 1.42x, +0.414, -0.436, +0.721, 40 s, 0.274 GB | VII | `s26/results/w_selfcopy_ladder.json`; `s26/LEDGER.md` L91 | as stored / as cited |
| 25 s (grace), six (cap), 22:03, 2.9 (median), 18, 22 | VII | `s26/LEDGER.md` L78, L80, L83, L92 | as cited |
| esm8m: +0.2418 (0.0740, 1.17x, [+0.1134, +0.3853]), +0.2246 (0.0913, 0.88x, [+0.1221, +0.3046]), +0.2111 (1.09x), +0.034 (0.15x), 251 s, 0.445 GB, 110 | VII, VIII | `s26/results/p_ladder_report_esm8m_s0.json`; `s26/LEDGER.md` L99, L101 | as stored / as cited |
| pairnet: +0.0407 (0.0486, 0.1361, 0.30x, [-0.0442, +0.1143]), +0.0360 (0.26x), +0.0077 (0.03x), 2.0791, 2.3386, +0.1285, 0.287, 587 s, 0.301 GB, 64, 4 | VII, VIII | `s26/results/p_ladder_report_pairnet_s0.json`; `s26/LEDGER.md` L103, L93 (the shipped MAE 2.3386) | as stored / as cited |
| validity axis: 2,568 s, 0.274 GB, 0.4444, 0.0079, 34, 1, [-0.5546, -0.3172], 2.05x, 34W/0L/92T, 3.4524, 0.1190, 3.30x, 81W/0L, 2.3629, 2.7846, 4.62x, 115, 1.3%, 12%, 2.5%, 6.61, 48, 0.33%, 0.9302, 0.9114, 0.84x, 0.0242, 0.0327, 0.021, 3.0792, 1000 | VII | `s26/results/ph_validity.json`; `s26/LEDGER.md` L100 | as stored / as cited |
| L95: 0.166, 0.22, eighth; L96: 3.320, 3.21, 0.083; L94: 101W/25L | VII | `s26/LEDGER.md` L94, L95, L96 | as cited |
| 50 of 55, 5, 2016 of 2016, 2142 of 2142, 7be8e4b0, 210.6 s, 0.632 GB, 8 | VII, IX | `s26/results/ast_gate_ae86a124.txt`; `s26/LEDGER.md` L98, L102 | as cited |
| B3: 100 s, 0.053 GB, 45, 300, 0.522, 0.498, 0.578, 0.263, 0.557, 0.566, 0.093, +0.244, -0.5366 (0.2309, 0.6469, 0.83x, [-0.882, -0.185], 84W/42L), +0.404, -1.1420 (0.3415, 0.9566, 1.19x, [-1.495, -0.700], 94W/32L), 0.463, 0.725, -0.612, -0.638, +0.463, -0.128, +0.394, -0.273, -0.266, +0.253, +0.230, -0.42, -0.636, -0.601 | VII, VIII | `s26/results/p_b3.json`, `s26/results/p_b3_features.json`; `s26/LEDGER.md` L106, L107 | as stored / as cited |
| amber_prior_partner: 135 s, 0.163 GB, 4,786 s, 0.131 GB, 63,000, 29, 30, +0.016 [+0.002, +0.030], +0.145 [+0.001, +0.260], 22 of 29, 0.02, 0.006, 0.007, -0.290, 2.804, 224%, 9.7, -0.022, 8%, 0.04, 0.19, -0.0010, 0.08x | VII, VIII | `s26/results/w_selfcopy_amberprior_endpoint.json`, `w_selfcopy_amberprior_clouds.json`, `w_amberprior_cells_cloud.json`; `s26/LEDGER.md` L105 | as stored / as cited |
| 3 of 3, 40.2 s, 0.313 GB, 360, 10, 0.05, 0.25, a6ce3ab6 | VII, IX | `s26/results/pytest_slow_equivalence.xml`; `s26/TEST_RUN.md`; `s26/LEDGER.md` L104 | as cited |
| Part B: 276, 231, 120, 1,049 s, 1,063 s, 893 s, 1.244, 1.250, 1.248, 145 s, 0.309 GB, 1.024, 0.76, 0.791, 0.89, 1.382, 0.77, 1.373, 0.91, +0.0114, +0.0204, +0.0115, -0.0062, +0.0177, -0.0021, -0.0002, -0.0040, +0.0019, -0.2460, -0.1634, -0.2310, -0.0265, -0.0318, -0.0188, -0.0243, 0.246, 0.10, 2.33, 1.38, 0.60, 4.13, 3.28, 0.0082, 0.0054, 0.0002, 0.008, 1,440 s, 1.247 | II, VII, VIII | `s26/results/w_selfcopy_endpoint.json`, `w_selfcopy_bound.json :: signed_bounds_gated/both_removed4`; `s26/LEDGER.md` L108 | as stored / as cited |
| H_P3: 9,321 s, 0.110 GB, +0.0066 (0.0157, 0.0440, 0.15x, [-0.0155, +0.0289]), -0.0179 (0.38x, [-0.0316, -0.0057]), +0.0001, -0.0284 (0.68x), -0.0043, -0.0117 (0.39x), +0.0183 (0.68x), 0.044, 4, 8%, 113, 30 | VII, VIII | `s26/results/w_selfcopy_provenance_readout.json`; `s26/LEDGER.md` L109 | as stored / as cited |
| C4: 2,839 s, 0.114 GB, -0.408, 1.43, -0.099, 12.7, +0.0684 (0.0269, 0.0754, 0.91x, 96%), +0.0707 (0.83x), +0.0267 (0.56x), +0.0237, +0.0244, +0.0206, +0.0594, -0.0012 (0.02x, 0.010), +0.0130, -0.203, +0.007 to +0.027, -0.10 to -0.21, 25, 15, 200, 100 | VII, VIII | `s26/results/p_c4.json`; `s26/LEDGER.md` L110, L115 | as stored / as cited |
| L112 to L116: 72.7, 6.19, 567 s, 0.128 GB, 6ed3b367 | VII | `s26/results/p_best_rung_chains.json`, `p_best_rung_chains_rebuild_basis.json`; `s26/LEDGER.md` L112, L114, L116 | as cited |
| integration tier: 8, 165.6 s, 1.139 GB, 11, 368, 2, -489.9138948277905, 4,620 s, seven, 75, 7ad4ef68 | VII, IX | `s26/results/pytest_slow_integration.xml`; `s26/TEST_RUN.md`; `s26/LEDGER.md` L111, L113 | as cited |
| the three verdicts (A REPLACE, B REPLACE, C KEEP WITH EDITS) and the slide 11 line | VII | `s26/PROPOSAL_A.md`, `s26/PROPOSAL_B.md`, `s26/PROPOSAL_B_REPLACEMENT.md`, `s26/PROPOSAL_C.md`, `s26/C3_RESULT.md`; `s26/LEDGER.md` L117 | as cited |
| A4 bootstrap: 3,928 s, 0.09 GB, 35, 70, 2,000, -0.649 [-0.714, -0.596], -0.252 [-0.318, -0.187], -0.047 [-0.171, +0.201], -0.311 [-0.356, -0.267], -0.243 [-0.294, -0.193], +0.571 [+0.500, +0.642], +0.658 [+0.596, +0.722], +0.346 [+0.265, +0.428], +0.301 [+0.232, +0.369], -0.056 [-0.182, +0.087], +0.006 [-0.097, +0.123], 13 of 128, 0.23 to 0.60 | V, VII | `s26/results/q_var_boot.json`; `s26/jobs_done/a4_var_boot.json`; `s26/LEDGER.md` L119 | as stored / as cited |
| C3 stage 2 verify: 3.2147652, 0.0 on 126 of 126, 5 s | VII, VIII | `s26/results/ph_c3_stage2_verify.json`; `s26/LEDGER.md` L118 | as stored / as cited |
| strain spread control: +0.452 [+0.280, +0.609], +0.451 [+0.286, +0.608], +0.433 [+0.247, +0.581], +0.241 [+0.054, +0.376], 0.007, +0.082 [-0.103, +0.259], [+0.011, +0.171], 0.39, +0.052, +0.756, 0.0005, 2,000 | VII, VIII | `s26/results/a_strain_vs_spread.json`; `s26/LEDGER.md` L121, L123 | as stored / as cited |
| L120: 192, 214, 0.5216, 0.5782, 0.5565, 0.5662, 0.4043, -0.2183 [-0.383, -0.030], -0.0856, +0.0340, 12.0, 2.76, af05d987, 80, 30 of 126 | VII | `s26/results/a_ladder_isolations.json`; `s26/results/p_b3.json`; `s26/LEDGER.md` L120, L124 | as cited |
| A3: 16,373 s, 0.352 GB, 45 s, 0.385 GB, 4.914, 3.2280, +0.0034 (0.0273, 0.04x, [-0.0446, +0.0337]), +0.0145 (0.18x), -0.0037 (0.05x), +0.1053 (0.77x, [+0.0367, +0.2112]), +0.0864 (0.69x), +0.1116 (0.81x), +0.0345, +0.0378, 0.30x, 0.32x, +0.0885 (0.55x), 14 of 78, 0.010, 3, 124 of 126, 1.34, 3.10, 0.01, 5.99, 0.565, 8e-4, 0.050, 0.263, 57, 10, 2.8, 5.9, +0.09 to +0.11, 0.7 to 0.8x, 0.041, 1.18%, 0.54 | V, VII, VIII | `s26/results/a3_stats.json`; `s26/results/a3/<pdb>.json`; `s26/jobs_done/a3_build.json`, `a3_label.json`; `s26/LEDGER.md` L125 | as stored / as cited |
| C5: 5,284 s, 0.09 GB, 45, 0.5, -0.0090 (0.0140, 0.23x, [-0.0331, +0.0185]), +0.0312 (0.60x), +0.1643 (0.0328, 1.79x, [+0.1176, +0.2253], 40W/86L), +0.1351 (1.33x), +0.0123 (0.45x), -1.8749 (120W/6L), -3.1293 (126W/0L), +0.0360 (0.84x), +0.0164, +0.1492 (1.91x), -3.0483, 3 of 5, 1.9 to 3.1 | VII, VIII | `s26/results/p_c5.json`, `p_c5_stats.json`; `s26/LEDGER.md` L132 | as stored / as cited |
| raw: 4,477 to 4,900 s, 1.93 to 1.96 GB, 4 of 5, 80, 10 | VII, VIII | `s26/LEDGER.md` L133 | as cited |
| rotamer relief B: 6,890 s, 0.313 GB, 41, 0.67 s, 9,450, 0.535 (0.026), 0.233 (0.024), 0.311, 0.057, 0.27, +0.0000 (0.0202), -0.0018, +0.0043, +0.0055, -0.0019 (0.04x), 17.5, 40.1, 4.3, 23.3, +0.0297 ([+0.0188, +0.0404], 0.74x), +0.0220 ([+0.0106, +0.0318], 0.61x), -0.0782 (1.12x), +0.0188 (0.59x), +0.017 (0.58x), 1.4 h, 10 of 14 | VII, VIII | `s26/results/ph_relief_run.json`, `ph_relief_report.json`, `ph_relief_reject.json`; `s26/LEDGER.md` L131 | as stored / as cited |
| rebuild: 4,472.9 s, 4,765 s, 0.11 GB, 2016 of 2016, 2,142, 8, 1008, 2, 083c9b95 | VII, IX | `s26/results/resultslab_rebuild/post_verdict.json`; `s26/LEDGER.md` L127 | as cited |
| governor v2.4: 18, 0, 16, 28, 13.3, 46, 32 | VII | `s26/LEDGER.md` L128, L129 | as cited |
| deck final: 520, 249, 248, 249, 250 | VII | `s26/pr_values.json`; `s26/pr_verify.txt`; `s26/LEDGER.md` L126, L134 | as cited |
| R7 to R10: 120 of 126, 1e-6, +0.2246, 0.88x, -0.612, -0.638, 1.5921, 2.2812, 2.0900, 2.7313 | VII | `s26/RETRACTIONS.md` R7 to R10; `bench_results/cache/1fc9f2dcf489e2fb/{1D6X,1KWE}.json :: rmsd_avg`; `s26/LEDGER.md` L9, L57, L101, L107 | as cited / as stored |
| L135: 13 of 14, 29, 18, 137, 11 of 11 | VII | `s26/DELIVERABLES_CHECK.md`; `s26/LEDGER.md` L135 | as cited |
| verify re-run: 32, 18, 19, 30, 31, 123, 8, 1, 126 of 126, 0 of 8, 22, twelve, 02:59 | VII, IX | `s26/results/verify/REPORT.md`, `s26/results/verify/REPORT.json`; `s26/results/verify/*.rerun.json` | as stored |
| A2 per step: 7 [7, 82], 48, 11 [7, 139], 15, 58 [7, 530], 18, 37 [7, 513], 10, 16 [7, 289], 1025 [513, 2017], 1025 [258, 2017], 16 [9, 161], 530, 2080, 30, 63, 60, 68, 36, 48 | V, VII | `s26/results/q_dla_a1.json`; `s26/LEDGER.md` L138 | as stored / as cited |
| lane W close: 12, 32, 31, 1.250, 9,321 s, 2 h, 21 min, 25, 16 | VII, VIII | `s26/agentW_FINDINGS.md`; `s26/LEDGER.md` L137 | as cited |
| A1 seed 1: -0.0449 (0.0227, 0.0635, 0.71x, [-0.0960, -0.0037], 72W/54L), -0.0523 (0.0235, 0.0659, 0.79x, [-0.1087, -0.0054], 77W/49L), 3.2161, 3.2087, 3.2610, 3.2280, 3.2142, 3.2059, +0.033, -0.0498 (0.76x), -0.0241 (0.28x), -0.0194, +0.0696, +0.0045, -0.0561 (0.65x), -0.0615 (0.69x), 1.40, 1.29, 25 s, 0.396 GB, 0.200 | V, VII, VIII | `s26/results/a1s1_stats.json`; `s26/results/a1s1/<pdb>.json`; `s26/jobs_done/a1s1_label.json`; `s26/LEDGER.md` L139, L140 | as stored / as cited |
| verify final: 20 of 20, 5, 11, 3, 1, 0, 1.732e-13, 6.828e-06, 1.540, 1894, 12431, 508, 469, 3.2145, 3.2057, 587.0 s, 1798.9 s, 1308.6 s, 0.275 GB | VII, IX | `s26/results/verify/REPORT.md`; `s26/LEDGER.md` L141 | as stored / as cited |
| close: 184, 18, 29, R1 to R11, 13 of 14; deck 550, 248, 248, 249 | VII | `s26/jobs_done/`; `s26/pr_values.json`; `s26/LEDGER.md` L142, L143 | as cited |
<!-- APPENDIX B ROWS -->

## APPENDIX D. RECONCILIATION WITH docs/REPORT_S26.md

`docs/REPORT_S26.md` (181 KB, Parts I to X, Appendices A to D) and `docs/REPORT_S26_SUMMARY.md`
appeared on disk during the 2026-09-13 pause, written from this branch by another session or by
the user, not by any S26 lane, and are left untouched and uncommitted (`s26/LEDGER.md` L42). This
appendix records what this report took from them, where the two disagree and which artefact
decides, and what was deliberately not taken. Every number taken was re-read from the named
artefact and is a row of Appendix B; the docs file is cited as the pointer.

### D.1 Taken (artefact-sourced numbers this report lacked)

| item | artefact | now in |
|---|---|---|
| production built chain per-target spread: median 2.9661, best 0.1824 (1S9Z), worst 8.2342 (2MQ2), 28.6% under 2 A, 50.8% under 3 A, correlation with the pool best 0.7150, rebuild basis delta 0.1643 | `results/summary/leaderboard.json :: rows[0]` | II.2 |
| the seven configurations against production, built chain, effect and x MDE | `s25/results/phys_suite.json :: vs_incumbent` (point cloud); `results/summary/leaderboard.json :: rows[1]` (distogram row, built chain, +0.0061 at 0.19x MDE) | IV.4 |
| the distogram's calibration: z sd 1.9962, 90% coverage 0.6159, 50% coverage 0.2824, multimodal fraction 0.2412 | `s25/results/calib.json :: rows[*]` (means over 126 rows) | III.2 |
| the prior ladder's first rung: 2.8334 point cloud, -0.2150 at 2.4x MDE, 113W/13L | `s24/LEDGER.md` L13 | III.2 |
| the score's per-pair target quantised to 17 values, up to 1.75 A of location error at long separation; 51 achievable arms, corr(gamma_eff, endpoint delta) +0.054 | `s25/LEDGER.md` L6, L12 | III.2 |
| statevector check 5.551e-17 (this report had rounded to 5.6e-17) | `s25/results/q_verify.json :: results` | V.2 |
| the S25 Gibbs falsifier (KL below 0.1 nats) fired against the draft | `s25/results/q_gibbs.json :: results/falsifier` | V.7 |
| the S5 measurement of the gradient defect, cosine -0.023, and the consolidation audit's exact / sampled tail cosines 0.6847 / 0.685 | `docs/FINDINGS.md:392-427`; `verify/cvar_audit.json` | V.5 |
| the entangler-deletion contrast's source, -0.013 [-0.095, +0.077] | `s20/LEDGER.md` L1 | V.7, Appendix B |
| the optimised arm's end-to-end speed-up, 8.371x (2537.568 s to 303.137 s, both cold, 1 against 8 workers) | `bench_results/compare_tuning126.json :: wall_s, end_to_end_speedup` | IX.4 |
| B1's feasibility numbers | `s26/results/b1_feasibility.json`; `s26/LEDGER.md` L13 | VII.2 |

### D.2 Disagreements, and the artefact that decides

| the two readings | decision |
|---|---|
| 3.2126 (their headline's leaderboard rebuild) against 3.2148 (this report's production cache); 0.1643 against 0.1664 for the projection gap; T030 0.18242 against 0.18198 | Not a disagreement: the same built chain on two emission paths (the results lab re-projects the stored cloud and round-trips through PDB quantisation, up to 0.171 A per target, `s26/LEDGER.md` L14; `s26/EXAMINATION.md` H). Both reports say so. This report quotes the cache number as the production result and names the rebuild beside it. |
| `bench_results/compare_tuning126.json :: science.rmsd_arm` (theirs) against `bench_results/baseline_tuning126.json :: science/rmsd_arm/mean` (this report) | The same values: `compare_tuning126.json :: science_delta` is 0.0 on every science key except two at 4.4e-16. Either path is valid. |
| 5.6e-17 (this report, from `s25/QUANTUM.md`) against 5.551e-17 (theirs) | The artefact says 5.551e-17; corrected here. |
| the benchmark leak count: "13 of 60 targets have a pool window at >= 0.6 identity" (their B.2, from `docs/FINDINGS.md:4595`) against the 2/60 of this report | Two criteria, both true: 13/60 is the >= 0.6 identity count of S9-10 (priced +0.0030 A in S10-4, artefact in git history at `5fa05cd`, L31); 2/60 is the verbatim self-copy defect of S24 L4 that lane W bounded (L44). This report keeps 2/60 as the declared leak and now names the 13/60 criterion beside it (II.1). The two benchmark target ids printed in `docs/FINDINGS.md` and in their table are not repeated here. |
| "each fold model saw ~100 of the other 125 dev natives": their citation `s20/LEDGER.md` L-B against this report's `s24/LEDGER.md` L8 | The sentence is at `s24/LEDGER.md:668` (entry L8); the S20 entry carries the fold-CI discussion. This report's citation stands. |
| the blind pipeline: their 3.989 against shipped 3.213 (all 126) beside this report's 3.749 against 2.745 (the 108 ordinary targets) and 5.425 against 6.019 (FAIL18) | Both from `s12/SPRINT12_DOSSIER.md` section I on different subsets; the 126-target pair is added to VI.8 for completeness. |
| ESM against one-hot, -0.288 A: their Appendix D9 records that the per-target artefact `s7/repr_tune.json` is lost (`s26/LEDGER.md` L11) | Agreed; this report's row cites `docs/FINDINGS.md` S7-11 and now carries the lost-artefact caveat. |
| their "Sprint 26, for context only" and "V.12 ADAPT-VQE (in progress)" | Written at 14:25, before L38 to L50 landed; Part VII of this report supersedes them. |

### D.3 Not taken

- Their Appendix B.7 "External (audience)": four items read from the web (arXiv abstracts, an
  author page) and the supervisors' names and affiliation. No external claim and no personal
  detail enters this report (`s26/LEDGER.md` L42 item 1).
- Their Part X (the presentation guide: the ninety-second statement, slide notes, flashcards,
  the postdoc questions, the must-not-claim list). That is lane PR's deliverable; this report is
  the record it draws on. Their must-not-claim list agrees with this report's Part VIII and the
  S25 retractions in every item checked (the +0.113 A CVaR contribution, the 3.0483 headline, the
  0.084 A MDE constant, "the circuit attains its own optimum").
- Their Appendix D12: numbers quoted from `~/.claude/projects/.../memory/*.md` (the 82.8th
  percentile, the 0.986 / 0.600 in-band figures, the 51st percentile, the 6.43 core-equivalents
  as a measurement). A memory file is not an artefact; none of those numbers is quoted here
  except 6.43, which this report attributes to `docs/STATE_BRIEF_2026-09-12.md` section 8 as an
  operational estimate, not a measurement.
- Their Appendix C (the corrections ledger) restates `docs/FINDINGS.md:60-142` and the sprint
  ledgers; this report's Part VI carries the same corrections sprint by sprint and adds nothing
  from their table.
- Their reading of `docs/FINDINGS.md:4595` names two benchmark targets; not repeated here.

## APPENDIX C. THE S26 LEDGER

`s26/LEDGER.md` reproduced verbatim at the close of the sprint, after the coordinator's
sprint-close entry, as the campaign prompt requires. It is the record this report cites by
entry number throughout. Being verbatim it is exempt from the report's own style rules (the
check in `s26/e_report_check.py` scans the prose before this appendix): it contains no em or en
dash, and two of the six words the report's own prose never uses appear in it as the ledger's
own text (L60 names the word this report removed at 92d559bb; L77 uses the adjective form of
another for a seed re-run); they are quotations of the record, not this report's prose.

Reproduced from `s26/LEDGER.md` at commit `95652f6b` (2026-09-14 04:18, through L144), unchanged; 147 numbered entries
(L1 to L144, with the suffixed collisions). The report's own text ends here; everything below is
the ledger.

---

# SPRINT 26 -- DECISION LEDGER

Append-only. One entry per finding or decision. Nobody edits another lane's entry; a
retraction is a new entry naming the entry it retracts. Every number carries an artefact path.

---

## L0 -- THE GOVERNOR IS RUNNING; THE BASELINE LOAD LEAVES ABOUT 4 GB (2026-09-12, coordinator)

`s26/governor.py` started at 23:50:41 (pid 36196), commit `a4db170c`. Self-test: a 0.3 GB job
registered through `s26/jobrun.py`, ran, and reported peak RSS 0.328 GB
(`s26/jobs_done/selftest_alloc.json`). Every action is logged to `s26/governor.log`.

Measured before any campaign process existed (`python s26/governor.py --once`):

    RAM total 16.75 GB (15.6 GiB)   used 11.16 GB = 66.6%   available 5.59 GB   CPU 10.5%

The 66.6% is the user's own environment (several VS Code and claude sessions, Defender;
`psutil` process table, 2026-09-12 23:47). The campaign's 93% ceiling therefore leaves about
4.4 GB for all lanes together. The "90 to 93% band" instruction is read as: fill up to, never
above, 93%; the governor launches queued work only when the box has sat under 88% for 60 s.

---

## L1 -- OPERATIONAL ITEM 1 CLOSED: `s24/cache_amber/*.npz` IS TRACKED (2026-09-12, coordinator)

Commit `6b494780`, the first on branch `s26`: `.gitignore` whitelist `!s24/cache_amber/*.npz`
placed after the blanket `*.npz` rule; 126 files, 1.5 MB, `git ls-files s24/cache_amber`
returns 126. S25 L13's single-point-of-failure on the 63,000 ff14SB/GBn2 single points is
closed.

---

## L2 -- THE PRESENTATION FILE DOES NOT EXIST ON THIS MACHINE (2026-09-12, coordinator)

`vqe_research_overview.pptx` is named in the campaign prompt as the 11-slide deck to tune. It
is not in the repository, not in git history (`git log --all --diff-filter=A -- '*.pptx'`
returns nothing), and not anywhere under `C:\Users\abena` to depth 7 excluding `AppData`
(`find` over Desktop, Documents, Downloads, OneDrive and the profile root; the only `.pptx`
files are six unrelated lecture decks in Downloads and python-pptx's template).

Decision: Phase 4 will build the deck from the structure the prompt describes (title; six
"what I built" slides; three proposal slides; goal and ask; dark theme) using python-pptx,
and `s26/PRESENTATION_CHANGES.md` will record every slide's content and its artefact so the
edits map onto the real file if the user supplies it. No number will be invented for a slide
that has no artefact.

---

## L3 -- `TEST_RUN_RESULT_PLACEHOLDER` IS NOT IN THE TREE (2026-09-12, coordinator)

The prompt says the state brief contains the literal string `TEST_RUN_RESULT_PLACEHOLDER`.
`grep -rn` over every `.md`, `.py` and `.txt` finds no occurrence. The brief on disk
(`docs/STATE_BRIEF_2026-09-12.md` section 6.1) already carries a live run: "355 passed, 13
skipped, 0 failed, 0 errors, exit 0" dated 2026-09-12. Operational item 5 is therefore "replace
that line with the S26 run under the governor, with per-file counts and the skip reasons", and
the Infrastructure lane owns it.

---

## L4 -- BRANCH AND READING ORDER (2026-09-12, coordinator)

Branch `s26` created off `ae86a124`. Coordinator read, in this order and in full: the state
brief, README, ARCHITECTURE, professor brief, S25 LEDGER (1,270 lines), S25 QUANTUM.md (857),
S25 agentQ_FINDINGS, S20 agentD_FINDINGS (1,006), docs/FINDINGS.md (all 9,081 lines, corrections
ledger first), docs/CONDENSED_REPORT, the S24, S23, S22, S21, S20, S19, S18, S17, S16, S15 and
S14 ledgers, the S13 and S12 dossiers, S13 qarch_FINDINGS part 1, the S25 BRIEF and PREREG_Q,
s12/instrument.py, s24/stats_lib.py, core/pipeline.py, and the deployed-selector region of
core/quantum.py. The remaining code read (core modules, root modules, tests, verify,
resultslab) is assigned to the Examiner and Infrastructure lanes, whose module map is the
Phase 0 deliverable.

---
## L5 -- PHASE-GATE READING FOR THE CIS-PEPTIDE CENSUS, AND LANE PH SPAWNED (2026-09-13, coordinator)

The phase gate forbids any endpoint experiment (anything reading an RMSD to a native) before
"PHASE 0 SIGNED OFF". The cis-peptide census reads omega angles of the 126 dev natives and CA-CA
distances of pool windows; it computes no RMSD, selects nothing and changes nothing. It is an ORACLE
DIAGNOSTIC and is allowed before the gate. The projection-floor cost on cis targets (an RMSD) waits
for the gate. Lane PH (Physics) spawned at 00:20 with brief `s26/briefs/PH.md`: AMBER as a steric
reject filter, the cis-peptide gap, the steric singularity for Part IV, and the C3 matched-random
control (stage 1 on the production cache, stage 2 on the best C2 rung when lane P delivers it).
Lanes active: E, I, Q, P, PH (`s26/lanes.json`).

---

## L6 -- OPERATIONAL ITEM 4 CLOSED: THE S8-13 SOURCE LOG HAS A TRACKED EXCERPT (2026-09-13, lane I)

`docs/FINDINGS.md` S8-13 (line 3711) cites `s8/predictor_report.log`; the file lives at
`_archive/logs/s8/predictor_report.log` (untracked by decision, `_archive/` is gitignored),
sha256 `2647f640efb6d8f8e0398620647d503c366aa9a0209f17742d927366dd0f4740`, 13,588 bytes, 172
lines. The number with no other source is the transfer law `selected = 1.009 * ref + 0.011`
(r = 0.981, 2,520 pairs), log line 92: a search for `1.009 * ref` over every json/md/py/txt/log
file finds only `docs/FINDINGS.md` and that log. Now tracked verbatim, with its context lines and
the other S8-13 tables the log carries (converged X_fit 3.217, the 3.653 / 3.285 / 2.220 filter
ceiling, the ESM ablation), in `docs/sources/s8_predictor_report_excerpt.md` (commit `eb89c165`).
Only the excerpt is tracked; the log directory is not.

Two things the excerpt records that were not known before: S8-13's second regression
(`1.019 * ref - 0.001`, r = 0.984, "over converged references") is NOT in the log and has no
on-disk source at all; and S8-13's X_fit positive-phi rate 0.0551 does not match the converged
row in the log (0.0629). `python s26/examine.py` now checks the transfer law as claim
`s8_13_transfer_law` with the excerpt as fallback.

---
## L7 -- DEFECT 6c: pool_gate = WARN IS TWO TARGETS AND IS NOW EXPLAINED IN THE OUTPUT (2026-09-13, lane I)

`results/summary/results.json` (production rows, primary basis `built_chain_bb`): the two
per-target pool violations are **1D6X** (T007, chain 1.7809 A against its pool's ORACLE best
member 2.4194 A) and **1KWE** (T019, 2.4224 against 2.8625). The production cache confirms the
mechanism is averaging, not leakage (`bench_results/cache/1fc9f2dcf489e2fb/{1D6X,1KWE}.json`):
for 1D6X the raw coordinate average is already at 1.5921 A while the best single member of the
top-75 it averages is 2.4194 A (`top_m_best` = `pool_best` there), and the projected chain
lands at 1.7809; for 1KWE the average is 2.2812 against a top-75 best of 2.8625 and the chain
at 2.4224. The emitted structure is a coordinate average projected onto ideal geometry, not a
pool member, so on a single target it can land nearer the native than any one window; it cannot
do so systematically, which is why the aggregate (mean 3.2126 against 1.7108, margin 1.50 A) is
the release condition (S25 L11). The distogram row shows the same two targets, amber_distogram
one (1KWE), legacy_distogram two (1KWE, 8IL1); the four physics rows show none.

Change (commit `eb89c165`, `s25/resultslab/schema.py`): `fmt_leaderboard` appends one line
under the table naming the WARN rows and the count of targets and stating the mechanism, and
the persisted `pool_gate_rule` string now defines PASS / WARN / FAIL in the same words (the site
shows that string on its overview card). Verified on the on-disk payload
(`s26/logs/precheck_fmt_leaderboard.log`). **The files under `results/summary/` were built
before this change and still carry the old rule string; the explanation appears on the next
`--mode frozen` build. `results/` was not rebuilt.**

---
## L8 -- DEFECT 6d: PRODUCTION NORMALISES BY RANK; EVERY MOMENT Z-SCORE IS RESEARCH-ONLY OR NOT ON AMBER (2026-09-13, lane I)

The production selector currency is `core.pipeline._zrank` (line 788: `rankdata`, then
`(r - mean) / sd`), applied at line 838 to the distogram score of the top 2^n candidates in
`quantum_stage`. It is the only normalisation in `core.pipeline`. Two qualifications: the
production `Config` has `quantum=False`, so in the 3.2148 A run the selector does not execute at
all (it runs under `--components`); and `_zrank` is applied to the distogram score, never to an
AMBER energy, because no AMBER energy enters the deployable selector.

Every remaining `(x - mean) / std` on an energy, with a verdict (`grep` over `core/`,
`s25/phys_*.py`, `s25/resultslab/`, `s24/`):

| path | what is standardised | verdict |
|---|---|---|
| `s25/phys_lib.py:98 zmoment` + `combine(norm="moment")` | raw AMBER / Legacy / distogram channels, no outer re-standardisation | RESEARCH ONLY: the DECLARED SECONDARY of the S25 seven-configuration suite; this is where the non-monotone float64 behaviour was measured (40/126 targets). Imported by `s25/phys_suite.py`, `s25/phys_analyse.py` only |
| `s25/phys_landscape.py:132`, `s25/phys_gate.py:71` | per-target energy vectors `e` | RESEARCH ONLY: S25 diagnostics |
| `core/energy.py:1066 LegacyField._standardise` | the 11 Legacy TERMS against 512 random structures from the target's own library | NOT AMBER, NOT PRODUCTION: the optional `legacy` term of `core.quantum.FoldObjective` (`w_legacy=0.0`, "Off by default"); Legacy terms are bounded so the moment is monotone there; `core.pipeline` never instantiates it |
| `core/quantum.py:2126 build_distogram(models="both")` | two distogram models' scores against random structures | NOT AMBER, NOT PRODUCTION: the generation lane's prior; `core.pipeline` does not call it |
| `core/predict.py:560` | input FEATURES of the distogram model (`sigma[k]`) | not an energy; training-time feature scaling |
| `s25/resultslab/exportlib.py:688` | per-target RMSD sd in the difficulty gate | not an energy |

`core/amber.py`, `core/bench.py`, `s24/d_harness.py` (`zrank`, line 275) and
`s24/d_hamiltonians.py` contain no moment standardisation. Conclusion: no production code path
applies a moment z-score to an AMBER energy; the rank normalisation is the one production uses,
and the S25 statement (ARCHITECTURE section 4) stands.

---

## L9 -- CORRECTION TO L7: TWO POINT-CLOUD NUMBERS WERE WRITTEN BEFORE THE QUERY RETURNED (2026-09-13, lane I)

L7 states the raw coordinate average of 1D6X "is already at 1.5921 A" and of 1KWE "2.2812".
Those two numbers are wrong: I composed the entry in the same turn as the cache query and
typed values that were not yet on screen. The production cache
(`bench_results/cache/1fc9f2dcf489e2fb/1D6X.json`, `1KWE.json`) reads:

    1D6X (n 13, fold 4): shipped 3.1750  pool_best 2.4194  top_m_best 2.4194  top_m_mean 3.5117
                         rmsd_avg 2.0900  rmsd_fit 1.9033  rmsd_arm 1.7814  rmsd_full 1.7473
    1KWE (n 16, fold 2): shipped 3.4830  pool_best 2.8625  top_m_best 2.8625  top_m_mean 3.9980
                         rmsd_avg 2.7313  rmsd_fit 2.4598  rmsd_arm 2.4127  rmsd_full 2.3944

The mechanism stated in L7 stands on the correct numbers: on both targets the point cloud
(2.0900, 2.7313) is already nearer the native than the best single member of the 75 it
averages (2.4194, 2.8625), and the projected chain stays below that member (1.7814, 2.4127 in
the cache; 1.7809, 2.4224 in `results/summary/results.json`, which re-projects and round-trips
through PDB, the same rebuild difference as 3.2148 against 3.2126). Every other number in L7
was read from `results/summary/results.json` or the leaderboard before it was written.

Rule for myself, recorded so it is not repeated: a number goes into the ledger only in a turn
AFTER the tool output that carries it has been read.

---

## L10 -- OPERATIONAL ITEM 3 CLOSED: THE SIX HELD UNUSED-IMPORT PROPOSALS ARE APPLIED (2026-09-13, lane I)

S25 L14 held them "until the results and physics lanes stop importing `core.pipeline`". The
condition, checked by grep: `s25/resultslab/*.py` and `s25/phys_*.py` never imported
`core.pipeline` (they import `core.project`, `core.geometry`, `s12.instrument`);
`s24/d_harness.py:275` and `s24/d_hamiltonians.py:130` name it in docstrings only;
`s23/c1_probweight.py:35` imports it (a finished sprint script); `s12/instrument.py:188`
imports it lazily inside one function. No live lane imports it, and at the moment of the edit
(`s26/governor_state.json`, 00:26:31) the governor had no job registered at all, so the
"never edit a module while a job launched from it is running" rule was met by measurement.

All six applied, commit `37bddbbb`, each verified with `s26/i_ast_check.py` (the working file
against HEAD, docstrings stripped, `ast.unparse` diffed) to differ by EXACTLY the removed names:

| file | removed | evidence it was unused |
|---|---|---|
| `core/amber.py:281-282` | `MAXITER_PER_PARAM`, `MIN_MAXITER_PER_PARAM` from the `budget` import | zero references repo-wide outside `budget.py` and the import; not in `__all__` |
| `core/bench.py:47` | `asdict` | single occurrence was the import |
| `core/cache.py:43` | `Iterable` | single occurrence |
| `core/predict.py:38` | `Dict` | single occurrence (`Dict`, strings included) |
| `core/data.py:440-442` | `_seq_index()` (dead `lru_cache` helper) | zero references in 700 modules |
| `verify/vqe_lfo_audit.py:19` | `defaultdict` | single occurrence; not on the production path |

Nothing changes a number: the removed names were never read. Tests after the edit, all under
the governor: `tests/test_pipeline.py` + `tests/test_data.py` (job
`pytest_nanpoison_post_core_edit`, 79 tests, 0 failed, 2 absent-artefact skips, exit 0) and
`tests/test_amber.py` (job `pytest_amber`, 16 passed, exit 0). The full non-AMBER suite is
re-run on the committed tree as the last act of the lane and recorded in `s26/TEST_RUN.md`.
`peptide_db.py` (the legacy arm) was not touched.

---

## L11 -- LANE P: THE COMPACT ESM TABLES COVER THE TRAINING CORPUS; THE BANK COSTS 1.79 GB; ONE S7 ARTEFACT IS LOST (2026-09-13, lane P)

`s26/results/p_probe_esm.json` (jobrun `p_probe_esm`, peak RSS 0.594 GB, 10 s): `s5/esmraw.npz`
(486.6 MB float32 resident, 1.3 s), `s5/esm32.npz` and `s7/repr_cache/esmcon.npz` each cover all
6,790 unique training sequences (787 peptides + 6,003 fragments) and all 126 targets;
`esm_small.npz` covers 360. Per fold 6,609-6,640 chains, 552,199-556,709 pairs.
`s26/results/p_probe_esmcache.json` (jobrun `p_probe_esmcache`, tag ESM, est 2.5 GB):
`esm_cache.npz` = 22,795 sequences, 1.79 GB resident (process peak_wset 1.793 GB; the 5-s jobrun
sampler read 1.456 GB, an under-read of the transient), 19.7 s. Decision: the C2 ladder reads the
compact tables only and patches `core.data.esm_raw/esm_embed/esm_contacts` so no path reaches the
bank. `s12.instrument.project` rebuilds 1A13's persisted `fit_ca` and `ca` at max abs 0.0 in 2.9 s.
`s7/repr_tune.json` (S7-11's per-target ESM contrast) is not on disk and not in git history
(`git log --all -- s7/repr_tune.json` is empty); only `s7/repr_oracle.json` survives. The S7-11
contrast can only be re-measured (rungs noesm / pca32 of `s26/p_ladder.py`).

---

## L12 -- LANE P: C2 CODE AND ITS PROBE. THE RETRAINED pca32 REPRODUCES THE PINNED FOLD-0 MODEL; THE LADDER COSTS ~38 MIN PER RUNG (2026-09-13, lane P)

`s26/p_ladder.py` (rungs shipped, noesm, conly, pca32, pca32f, pca128, raw, esm8m, wide, pairnet,
mix; endpoints sel / cloud / arm / fit; gam_eff and cos in probability and location space; the
phase gate is enforced in code). Synthetic tests `s26/p_ladder_test.py` (jobrun `p_ladder_test`,
10 s, 0.325 GB): lean trainer == `core.predict.MLP.fit` at max |dprob| 1.0e-7; rebuilt pca32
features vs the shipped construction max abs 9.5e-7; `mix` lam = 0 bit-exact; NaN-poison
bit-identical; the guard refuses unknown sequences. Probe (jobrun `p_probe_train_pca32_f0`):
fold 0, 552,199 pairs, 183-d, 40 epochs, 453 s, peak RSS 1.248 GB. Gate
(`s26/results/p_ladder_gate_pca32_fold0_s0_probe.json`, 25 fold-0 targets, native-free): posterior
max abs 1.5e-5 vs `distogram_models/fold0_esm_frag.pt`, risk max abs 6.8e-4, K = 500 argmin
identical 25/25, top-75 overlap 1.000, Spearman 0.99999998. Not bit-exact (float32 rounding of
the cached embeddings), equivalent in every quantity the pipeline consumes. A note for the record:
S7-11's "tri on ESM input has not been measured" is superseded on one basis by S19 L11 (PairNet on
the deployed inputs, -0.080 [-0.242, +0.082] through the distance-geometry fit, paired sd 0.910,
MDE 0.227, `s19/results/a_models.json`); unmeasured through the pipeline readout, which is rung
pairnet. Question for the coordinator: may the 5-fold TRAINING of the ladder (no native RMSD read)
run before sign-off, so eval can start at sign-off? Not launched pending the answer (rule 14).

---

## L13 -- LANE P: B1 IS INFEASIBLE ON THIS BOX ON THREE INDEPENDENT GROUNDS; B1 STOPS (2026-09-13, lane P)

`s26/results/b1_feasibility.json`. (a) Weights on disk: only esm2_t33_650M (2.60 GB) and
esm2_t6_8M (30 MB); no ESMFold or 3B checkpoint. (b) `esm.pretrained.esmfold_v1` exists but
`import esm.esmfold.v1.pretrained` raises `ModuleNotFoundError: No module named 'omegaconf'`
(esmfold.py:10); `openfold` (esmfold.py:11-13, trunk.py:11) is absent; the fair-esm README's
openfold install "requires nvcc" and "python <= 3.9" (box: Python 3.13, torch 2.13 CPU).
(c) README: esmfold_v1 = 48 (+36) layers, 690M (+3B) parameters; HTTP HEAD (network reachable in
< 10 s): esmfold_3B_v1.pt 2,771,653,574 B, esm2_t36_3B_UR50D.pt 5,678,116,398 B; fair-esm halves
the LM and keeps the trunk fp32 (esmfold.py:46) -> 8.76 GB resident minimum (fp16 everything
7.4 GB) against 4.4 GB of campaign headroom (`s26/governor_state.json` 00:19: 64.8% of 16.75 GB,
ceiling 93%). Nothing downloaded. The 650M contact head is exposed (`esm_features.compute`,
`return_contacts=True`), cached for every training sequence, and already 13 of the shipped prior's
183 input columns; S17 L23 measured it as a ranker (+0.116 [+0.047, +0.184] in-band at top-75,
+0.050 over its ESM-free twin, mostly compactness, no argmin gain, shortlist worse than matched
random). Its standalone value as a PRIOR INPUT is rung conly (PREREG_B2); the size axis is
esm8m vs 650M. Proposal B routes to B2/B3 and `s26/PROPOSAL_B_REPLACEMENT.md`.

---

## L14 -- LANE P: B3's PER-TARGET ARTEFACTS EXIST; THE ARM-CHOICE ORACLE OVER {pipeline, torsion predictor, helix} IS AN ORDER STATISTIC (2026-09-13, lane P)

`s13/cache/tors_rows.npz['a_pepPos']` (126 rows, `rmsd_build` mean 3.7705) and
`s14/results/ladder.json per_target['L0_constant_helix']` (126, mean 4.0648; its incumbent column
reproduces `rmsd_fit` 3.2041 at 0.0000). Paired against the production `rmsd_arm` (3.2148),
persisted artefacts only: tors - arm +0.5557 (median +0.2314), SE 0.1318, MDE 0.3694, 37W/89L,
top-10 share 0.567; helix - arm +0.8500 (median +0.1691), SE 0.1494, MDE 0.4185, 32W/94L;
tors - helix -0.2943, SE 0.0807, 69W/57L. FAIL18: tors 5.569 / arm 6.032 / helix 5.887;
other108: 3.471 / 2.745 / 3.761. Per-target min over the three arms 2.9632 (-0.2515 vs arm) is
92% accounted for by `ST.best_of_k_within`'s across-target null (k_eff 2.34). The S13 `base`
column (3.2126, the leaderboard rebuild) differs from `rmsd_arm` by up to 0.171 A per target;
PREREG_B3 pairs against the production cache. Also noted: the coordinator's L5 defines lane PH's
C3 (steric reject filter, matched-random control); lane P's `PREREG_C3.md` was drafted without
that brief and says so in its addendum; PH's brief governs.

---

## L15 -- DEFECT 6a: THE IDENTITY FLAG SHIPS DARK; THE "CORRECTED" NORMALISATION IS AT THE NULL AS A CLUSTERING CRITERION; THE FOLDS STAY PINNED (2026-09-13, lane I)

**Code.** `core.data.identity(a, b, norm="longer")` and `identity_many(..., norm="longer")`
(commit `37bddbbb`): the default path is statement-for-statement the pinned convention;
`norm="shorter"` is the same match count over the SHORTER sequence behind a verbatim-substring
test. `clusters()` / `folds()` do not accept it; nothing on the production path passes it;
`peptide_db.py` (legacy arm) untouched. Test:
`tests/test_data.py::test_identity_norm_flag_ships_dark_and_the_shorter_form_catches_the_self_copy`
(default bit-identical; 1CEK-in-1A11 scores 1.0 under the flag; scalar and batched forms agree
on both sides of `_BATCH_MIN`).

**Audit, in memory only** (`s26/i_identity_audit.py`, job `i_identity_audit3`, 0.056 GB,
`s26/results/i_identity_audit.json`; `peptide_folds.json` and `peptide_clusters.json`
sha256-identical before and after and equal to lane E's `s26/results/pinned_hashes.json`).
The in-memory copy of `core.data.clusters` reproduces the pinned clusters AND the pinned folds
exactly, so the copy is faithful.

1. The brief's corrected criterion (shorter normalisation + substring, same 0.6 threshold)
   collapses **470 clusters to 166**; 594/787 sequences gain or lose a mate; **71/126 dev
   targets** gain a mate sitting in another pinned fold; 48/60 benchmark sequences (count only).
2. That count is chance, not leak. Pre-registered null (`s26/PREREG_identity_null.md`, written
   before the control ran): a real dev sequence passes the shorter criterion against **0.5%** of
   the other 786 members, a composition-preserving shuffle of it against **0.3%** (ratio 0.62;
   the falsifier was ratio < 0.2). At mean degree ~3 on chance edges, single linkage percolates.
   The pinned longer criterion: 0.07% real, 0.001% shuffled. So `norm="shorter"` at 0.6 is a
   per-pair leak TEST, not a replacement clustering threshold; the docstring now says so.
   (My PREREG expected pass rates of 10 to 30%; they are 0.5%. The at-the-null verdict stands,
   the magnitude I predicted was wrong by 30x.)
3. The minimal fix, pinned criterion OR verbatim substring: 459 clusters; membership changes
   on 23 sequences, 20 of them verbatim containments and 3 carried along by single linkage
   (their pinned cluster mate is one); PREREG H3 as written is falsified by those 3, by
   transitivity. Dev targets that gain a mate in another pinned fold: **5** -- the four declared
   self-copies **1CEK** (fold 2; carrier n=25 in fold 3; longer identity 0.520), **2FBU** (4;
   carrier n=23 in 0; 0.522), **2P5H** (4; carrier n=17 in 2; 0.529), **6B9K** (0; carrier n=17
   in 2; 0.588), each with shorter identity 1.0 and each moved into one fold by the fix, plus
   **8ZG2** (fold 1), whose own verbatim carrier (n=23, longer identity 0.478) is in fold 1 and
   which gains a cross-fold mate only because that carrier also carries a peptide from another
   fold. 14 further dev targets have a verbatim relative that already sits in their own fold.
   Benchmark under the minimal fix: **2/60** gain a cross-fold mate (count only), which is the
   S24 L4 figure.
4. The trap, measured: even the minimal fix, re-derived through `folds()`' shuffle, would
   relabel **647/787** sequences (102/126 dev under the shorter criterion), because cluster ids
   renumber and the seeded permutation changes with them. Any future repair must PATCH the
   pinned fold map (move the copies) and retrain only the affected fold models, never re-derive.

Nothing was written to the pinned files; nothing ships. Sign convention for the record: 4/126
dev and 2/60 benchmark stand as the self-copy counts.

---

## L16 -- DEFECT 6b: THE AMBER MEMORY GUARD NAMES THE CEILING AND QUOTES THE GOVERNOR; BOTH AMBER FILES CARRY IT (2026-09-13, lane I)

`tests/test_amber.py` (commit `eb89c165`): the autouse `_memory_ceiling` fixture keeps its
decision exactly (`_memory_verdict` on `core.amber.memory_percent()`, the syscall
`core.amber.memory_guard` itself uses, so the verdict is identical with or without a governor);
what changed is the text. `_governor_reading()` reads `s26/governor_state.json` and returns
its RAM / CPU / job counts / band ceiling when the snapshot is under 120 s old, else `None`
(absent, stale, malformed). `_ceiling_message()` then says "physical memory N%
(core.amber.memory_percent) is above core.amber's 92% ceiling; the S26 governor read X% RAM /
Y% CPU at <ts> (<age> s ago) with J registered job(s), A AMBER, governor ceiling 93.0%", or
"no fresh s26/governor_state.json, so the governor is not running", followed by the existing
"this suite did not cause it" / "THIS SUITE accounts for it" clause. Unit test
`test_the_ceiling_message_names_the_ceiling_and_the_governor_reading` (synthetic snapshot
files; fresh, stale, absent and malformed cases; both verdicts).

`tests/test_amber_frame_invariance.py` imports that fixture object (`from test_amber import
_memory_ceiling`), which registers it as autouse in that module: the S25 L14-approved change,
shared rather than copied, so memory pressure there is now a skip that names the ceiling instead
of a red `MemoryError` out of `core.amber.memory_guard`.

Under the governor, on the tree of `37bddbbb`: `tests/test_amber.py` job `pytest_amber`
**16 passed** (15 + the new test), 260.5 s, peak RSS 0.872 GB; `tests/test_amber_frame_invariance.py`
job `pytest_amber_frame` **3 passed**, 255.7 s, peak RSS 0.324 GB; the guard fired in neither
(box 65 to 75% throughout, `s26/governor.log`). Record: `s26/TEST_RUN.md`,
`s26/results/test_run.json`, `s26/results/pytest_amber.xml`, `pytest_amber_frame.xml`.

---
## L16b -- PHASE-GATE READING FOR MODEL TRAINING; C3 OWNERSHIP; REPRODUCTION EXACT (2026-09-13, coordinator)

(Numbered L16b: lane I's L16 landed while this entry was being written; rule 4 below applies.)

1. Lane P asked (L12) whether the C2 ladder's 5-fold TRAINING may run before sign-off. Training a
   distogram head on the training folds reads native distances as labels, exactly as the pinned
   production models were trained; it reads no RMSD, selects nothing, and produces no claim. It is
   allowed before the gate under three conditions: one training job at a time through the
   governor (probe peak 1.25 GB, L12); rung order fixed by `s26/PREREG_C2.md` before the first
   job; no evaluation of any rung (no RMSD) until "PHASE 0 SIGNED OFF" is posted.
2. Proposal C's item C3 (AMBER relaxation of the best C2 rung's output against the S16
   matched-magnitude random displacement control) is owned by lane PH (L5, `s26/briefs/PH.md`,
   prereg to be `s26/PREREG_c3_control.md`). Lane P's `s26/PREREG_C3.md` describes a different
   hypothesis (AMBER-weighted pool histogram as a mixture partner of the prior; relaxed training
   labels). It is not deleted: lane P appends an addendum saying it is not C3, and enters the
   cheap arm in the tournament as `s26/IDEA_amber_prior_partner.md`.
3. Lane E's reproduction (`s26/results/e_reproduce.json`, job `e_reproduce`, exit 0): fresh
   versus stored means agree to 0.0 on all four bases (rmsd_avg 3.048338, rmsd_fit 3.204076,
   rmsd_arm 3.214765, rmsd_full 3.235460, n=126); T030 = 1S9Z rmsd_arm 0.181981 fresh and
   stored; fold mismatches against the pinned folds: none. The 2.6e-04 A tolerance is met with
   zero disagreement. Phase 0 still waits for `s26/EXAMINATION.md` and the Adversary's audit.
4. Ledger numbering: five lanes append concurrently. Re-read the tail immediately before
   appending and take the next free number; if two entries share a number, the later one adds a
   suffix (for example L17b) rather than renumbering.

---

## L17 -- SESSION-LIMIT INTERRUPTION 00:46 TO 08:40; NOTHING LOST; LANES RESUMED FROM DISK (2026-09-13, coordinator)

All five lanes (E, I, Q, P, PH) were terminated at about 00:46 by the API session limit
("session limit, resets 03:30 America/Los_Angeles"), not by any error. The governor (pid 36196)
ran throughout: `s26/governor.log` is continuous, with SAMPLE lines every minute from 23:52 to
08:36. No governed job was in flight at the cut: the last registration reaped was
`pytest_core_post` at 00:45:56 and `s26/jobs/` was empty at 08:36; every record in
`s26/jobs_done/` carries an exit code. Every lane's files on disk stand (listing at 08:36:
preregs C1-C5, B1-B3, amber_reject, cis, c3_control, identity_null; eight IDEA files; scripts
e_*, i_*, p_*, ph_*, q_*; results in `s26/results/`); lane commits eb89c165, 37bddbbb, 9bb7f4f6,
601a39c7 are on the branch. At 08:40 each lane was resumed from its own transcript with an
instruction to commit what is on disk and continue, not restart. The phase gate is still closed
(`s26/EXAMINATION.md` and `s26/BRIEF.md` not yet written). Rulings added on resumption: lane Q's
A2 (DLA) and A4 (gradient variance versus width) read no native and may run before the gate;
A1 and A3 endpoints wait.

---

## L18 -- HOW L15's BENCHMARK COUNTS WERE OBTAINED: THE MANIFEST WAS READ THROUGH core.data.benchmark(), SEQUENCES ONLY (2026-09-13, lane I)

The counts in L15 that concern the sealed benchmark (48/60 gain a cross-fold mate under the
shorter criterion; 52/60 would be relabelled by a naive re-derivation; 2/60 under the minimal
pinned-OR-substring fix) were computed in `s26/i_identity_audit.py` (sections 3b and 5) over
the SEQUENCES of the 60 benchmark peptides returned by `core.backend("data").benchmark()`.
That function (`core/data.py:617-633`) opens `results/benchmark_manifest.json`, reads the `pdb`
field of every entry in its `targets` list, and returns the matching `Peptide` records of the
peptide database (`peptide_db.npz`). So the manifest's contents WERE read by the audit process,
by the route the lane brief named as the permitted one ("through core.backend("data").benchmark()
sequences only, print no names, report only the count"); stated plainly so the coordinator can
log it as the conservative path. The S24 L4 count (2/60) was made the same way.

What was and was not touched: the audit used only `.seq` of each returned record. The
`Peptide` records carry native CA / phi / psi for every one of the 787 database entries
(they ARE the database, loaded by every caller of `core.data.load()`), and no such field was
accessed for a benchmark entry. No benchmark PDB file was opened, no native coordinate read,
no RMSD computed, no benchmark name printed, saved or logged: `s26/results/i_identity_audit.json`
holds counts only (`benchmark60`, `substring_or_pinned.benchmark60_with_new_mates_in_other_pinned_folds`),
and non-dev database members appear in the audit output by length alone. Nothing further will be
computed on benchmark targets by this lane.

---

## L19 -- OPERATIONAL ITEM 2 CLOSED: verify/run_equiv2.sh RUNS THE COMPARISON IT WAS WRITTEN FOR; BASELINE AND THE SHIPPED MODE ARE BIT-IDENTICAL (2026-09-13, lane I)

Why it was stale: it exported `PROJECT_GRAD`, which `core/project.py` stopped reading when the
mode moved into `core.pipeline.Config.project_grad` so that it reaches the cache key
(`verify/grad_key_collision.py`). The `core.pipeline` CLI has no flag for the mode, so as
written the script ran three arms in one mode. Repair (commit `eb89c165`): the consolidated
arms go through `s26/i_run_equiv2_arm.py`, which builds the identical `Config` to
`python -m core.pipeline run` (the same `replace(PROD, ...)` call) and sets `project_grad`;
four arms (baseline legacy; `exact`, the shipped default; `fd`; `analytic`); the original
`cfg_key` grep kept; any argument passed to every arm. `s26/i_equiv2_compare.py` compares the
arms per target from their cache directories with `==`.

Run: job `run_equiv2`, `sh verify/run_equiv2.sh --no-amber`, tag AMBER (with `--no-amber` no
OpenMM context is created, but `openmm` is still imported because `core.pipeline` resolves
the amber backend at start, so the conservative tag was the right one), exit 0, 145.3 s, peak
RSS 0.596 GB (`s26/jobs_done/run_equiv2.json`). Per arm on smoke8 (8 targets, no stage 4):
baseline 67.4 s, exact 33.0 s, fd 20.2 s, analytic 8.8 s.

    arm              cfg_key            vs opt_exact (8 targets)
    baseline_legacy  65ec272db3d31f05   ca, fit_ca, phi, psi, avg_ca, rmsd_avg, rmsd_fit,
                                        rmsd_arm, n_windows, n_top: bit-identical 8/8
    opt_exact        66050f6daae4ca07   (the shipped default)
    opt_fd           85faafd84d76a827   avg_ca and rmsd_avg identical; rmsd_fit differs on
                                        8/8, max 6.9e-4 A; rmsd_arm max 1.2e-2 A
    opt_analytic     f4e137a48586bd4b   avg_ca and rmsd_avg identical; rmsd_fit max 1.0e-1 A;
                                        rmsd_arm max 4.8e-2 A, 8/8 targets

Four distinct keys: the mode is in the key and the collision the script was flagged for is
closed. Baseline against the shipped mode is bit-identical on every emitted array and scalar:
the consolidation is faithful, and the `fd` and `analytic` differences are the scan builder and
the analytic gradient, as `core/project.py`'s docstring states. (The raw coordinate max|d| of
22 to 32 A and phi/psi differences near 2 pi between modes are lab-frame and angle-wrap
differences of un-superposed emitted chains, not the science; the RMSD scalars are.)

Artefacts: `s26/results/run_equiv2.log` (whitelisted transcript: the jobrun log, all four arm
logs, the comparison), `s26/results/run_equiv2_compare.json`, `verify/e2_*.log` (ignored).
Not run: the AMBER-inclusive form; stage 4 is downstream of the projection the script compares.

---
## L18b -- GOVERNOR v2: SMOOTHED CPU, RAM-ONLY KILLS, 20 s MINIMUM SUSPENSION, LAUNCH CAP OF 4 (2026-09-13, coordinator)

(Numbered L18b: lane I's L18 landed first; the governor.log RESTART line says "ledger L18" and means this
entry. Lane I's L18 records that the benchmark counts in L15 came from `core.data.benchmark()` sequences
only, with no coordinate, native or RMSD read; the coordinator logs that as the conservative path:
membership and sequence through the project's own helper is allowed for counting, and nothing further
is computed on benchmark targets in S26.)

`s26/governor.log` 00:37:01 to 00:38:51: with five governed jobs registered the raw CPU sample
crossed 93% and fell under 90% on alternate 5 s ticks (SAMPLE lines 00:37:41 cpu 98.7%, 00:38:41
92.8%, RAM 67 to 68% throughout), so the v1 band suspended and resumed `q_dla_smoke` eleven times
in two minutes and suspended two PH census jobs. No job was killed and every job exited 0, but a
CPU spike sustained 15 s would have killed the newest job, which for jobrun-launched jobs has no
queue spec to requeue from. Changes (commit follows this entry), all in `s26/governor.py` and
`s26/jobrun.py`, no production module touched:

- CPU enters the band as a 15 s rolling mean (`cpu_smooth`, three samples), written to
  `s26/governor_state.json` beside the raw value.
- The kill rule (95% for 15 s) now reads RAM only; CPU above 93% suspends the newest running job
  and never kills. CPU overload slows the box; only RAM can crash it.
- A suspended job stays suspended at least 20 s; resumption needs RAM under 90% AND smoothed CPU
  under 80%.
- `jobrun.py` also waits while the smoothed CPU is above 85% or four jobs are already
  registered, so five lanes launching at once are serialised at the door instead of suspended
  after the fact.

The v1 governor (pid 36196) was terminated at 08:40:24 with the RESTART line in the log and v2
started in its place. Two lane-Q jobs (`a2_dla`, `a4_var`) were registered at that moment and ran
unsupervised for a few seconds; v2 picked their registrations up on its first sample. The log is
continuous through the restart.

---

## L20 -- OPERATIONAL ITEM 5 CLOSED: THE BRIEF'S 6.1 CARRIES THE GOVERNED RUN; THE ORIGINAL SENTENCE IS IN HISTORY (2026-09-13, lane I)

`docs/STATE_BRIEF_2026-09-12.md` was committed exactly as received as `6e50ea93`, so the
"355 passed, 13 skipped, 0 failed, 0 errors, exit 0" sentence survives in history; then its
section 6.1 live-run paragraph was replaced (commit `83305549`) by the S26 governed run:
**370 tests, 357 passed, 13 skipped, 0 failed, 0 errors, 0 memory-guard skips** on commit
`601a39c7`, in three jobs under `s26/governor.py` (`pytest_core_post`, the 10 non-AMBER files:
351 tests, 338 passed, 13 skipped, 240 s, peak RSS 1.692 GB; `pytest_amber`: 16 passed, 261 s,
0.872 GB; `pytest_amber_frame`: 3 passed, 256 s, 0.324 GB), with the passed count per file,
the 13 skip reasons (11 `VERIFY_SLOW=1` opt-ins in `test_equivalence.py` and
`test_integration.py`, 2 absent artefacts in `test_pipeline.py`) and the pointer to
`s26/TEST_RUN.md` / `s26/results/test_run.json`. The suite is 370 rather than 368 because S26
added one test to `test_amber.py` and one to `test_data.py`.

The same non-AMBER files run BEFORE any S26 edit (job `pytest_core`, tree `a4db170c`: 337
passed, 13 skipped, 0 failed, identical skip list) are recorded beside it and marked
superseded, so the before/after of the production-path edits (L10, L15) is on file.
`TEST_RUN_RESULT_PLACEHOLDER` never existed in the tree (L3); the live-run sentence held the
placeholder's role.

---

## L21 -- HYGIENE: README GOVERNOR SECTION; `make examine` IS `python s26/examine.py` (WITH examine.sh / examine.bat), WRAPPING THE MODULE MAP, THE PINNED HASHES AND THE CLAIM LEDGER (2026-09-13, lane I)

`README.md` (commit `eb89c165`) gained the section "The S26 resource governor": the three
scripts and what each does, the band (88% low water to launch queued work; resume below 90%;
93% ceiling suspends the newest job; 95% for 15 s kills it with CTRL_BREAK first and requeues
it), the AMBER cap of two, the tags CPU / AMBER / ESM / TEST, where the live state
(`s26/governor_state.json`) and the log (`s26/governor.log`, whitelisted) are, how to start it
(`python s26/governor.py`, `--once`, `--status`), how a job is run or queued, the three-way
split of the test suite, and the examine entry point.

There is no Makefile in the repository and no `make` on this box, so the `make examine`
equivalent is a script: `python s26/examine.py`, with the one-line launchers `examine.sh` and
`examine.bat` at the repository root (no name conflict). It calls, imported rather than copied:
lane E's `s26/e_module_map.py` (rewrites `s26/results/module_map.json`), lane E's
`s26/e_hashes.py --check` (re-hashes every pinned artefact against
`s26/results/pinned_hashes.json`; the benchmark manifest as bytes only), and lane I's
`s26/i_claim_check.py` over `s26/results/claims.json` (21 claims, each re-read from the artefact
it names: the four 126-target means of the production cache, the `compare_tuning126.json`
science deltas and baseline constants, the leaderboard's production row and gate, the pinned
hashes, the S8-13 transfer law with the excerpt as fallback, and the S26 test-run totals);
`--search` adds lane E's `s26/e_claims.py`. Exit status is non-zero on any drift, mismatch or
non-optional absence. Validated: 21/21 OK at 00:37 (`s26/results/claim_check.json`), and the
full run at 08:4x is `s26/logs/i_examine_full.log` (it rewrote lane E's `module_map.json` once;
the map is deterministic apart from the sha256 of the files S26 edited).

Also added for the contract's own rules: `s26/i_ast_check.py` (section 5, AST identity modulo
docstrings against a commit, with a readable diff of any deviation) and `s26/i_test_report.py`
(junit XML plus `s26/jobs_done/` into `s26/TEST_RUN.md` and `s26/results/test_run.json`, with
memory-guard skips counted apart from real skips).

---

## L22 -- CIS CENSUS: NO CIS PEPTIDE BOND EXISTS ANYWHERE ON THE INSTRUMENT; THE DATABASE STEP GATE, NOT THE PROJECTION, SETS THE COST (2026-09-13, PH)

`s26/ph_cis.py census`, `s26/results/ph_cis_census.json` (complete 126/126), job
`s26/jobs_done/ph_cis_census.json` (exit 0, 90 s, peak RSS 0.038 GB). Pre-registered in
`s26/PREREG_cis.md` section 2 with the prediction "zero on model 1 and in the pool, because of
the step gate; some cis bonds in other ensemble models". Allowed before the gate by L5: ORACLE
DIAGNOSTIC on the natives (omega angles only, no RMSD), native-free on the pool.

    model-1 natives with a cis bond, |omega| < 30 deg        0 / 126
    the same by consecutive CA-CA < 3.3 A                    0 / 126   (the two criteria agree 126/126)
    deposited models with a cis bond, every ensemble         0 / 1,966
    universe windows with a CA-CA step < 3.3 A               0 / 2,352,893   minimum step 3.5045 A
    K=500 pools and production top-75 with such a window     0 and 0
    omega non-planarity |180 - |omega||, 1,507 bonds         mean 1.91 deg, median 0.34, p90 5.8, p99 17.4,
                                                             max 42.8 (9UV5, bonds 0 and 6: -137.2 and +144.8 deg);
                                                             0.53% of bonds beyond 20 deg, 0.13% beyond 30 deg

The universe minimum step IS the gate: `core/data.py:406-407` drops any peptide with a
consecutive CA-CA step below 3.5 A and `core/data.py:697-698` drops any fragment window with
one. The 126 are drawn from that database (`s7/debias.py:103-120`) and so is the sealed
benchmark, so neither can contain a cis target and no pool can contain a cis window. The
two-bond-length projection is therefore worth exactly 0.000 A on this instrument by
construction of the database, not by any property of the projection. The residual cost of the
constant omega here is the non-planarity tail (0.5% of bonds beyond 20 deg), which part 2
(gated) prices as the ideal-trans floor on the native's own torsions. The registered prediction
that other ensemble models carry cis bonds was WRONG: 0 of 1,966. The cis-peptide question is a
world-supply question (the 16 containment-fresh targets, 10 amyloid) and the design note in
`s26/agentPH_FINDINGS.md` section 1.4 is written for that supply.

---

## L23 -- STERIC REJECT CENSUS: AT 1e4 kcal/mol THE REJECT REMOVES 40 OF 75, EMPTIES 11 TOP-75 SETS AND 8 WHOLE POOLS, AND 96.8% OF THE CATASTROPHES ARE SIDE-CHAIN CONTACTS (2026-09-13, PH)

`s26/ph_reject.py census`, `s26/results/ph_reject_census.json` (complete 126/126), job
`s26/jobs_done/ph_reject_census.json` (exit 0, 100 s, peak RSS 0.117 GB). Native-free: the
63,000 cached single points (`s24/cache_amber`, pool identity `universe_idx == I.pool_idx` and
production `sub` == score top-75 asserted on every target) and the ideal-geometry rebuilds; no
RMSD, no native. Required by `s26/briefs/PH.md` section 3.1 before any RMSD is read; the prior
was registered in `s26/PREREG_amber_reject.md` section 5.

    threshold     pool frac > T   top-75 rejected   zero-reject   all-75 rejected   no survivor in 500   refill depth (mean/median)   R overlap with anchor
    1e3             0.760           54.5 / 75          1              25                 20                270 / 234                    0.27
    1e4 PRIMARY     0.586           40.1 / 75          2              11                  8                201 / 147                    0.46
    1e5             0.446           30.3 / 75          6               3                  1                167 / 118                    0.60
    1e6             0.345           23.3 / 75          8               2                  0                135 / 103                    0.69

The 0.586 reproduces S25's 58.6% (`s25/results/phys_landscape.json`). At the primary threshold
the operator is not a small surgical reject: it removes more than half of every shipped set,
reaches rank 147 (median) of 500 to refill, empties the whole top-75 on 11 targets (1G89 1ID6
1LB7 2MAI 2NB7 2XL1 5MML 5Z5W 7BX2 8UN8 9S5G) and finds no survivor among all 500 candidates on
8 (1G89 1ID6 2MAI 2NB7 2XL1 5MML 7BX2 8UN8). Both empty cases fall back to the anchor, a rule
added in the prereg's addendum 1 before any RMSD is read. Native-free geometry of the retained
set at 1e4: R is +0.254 A more expanded in Rg than the anchor (SE 0.058) and +0.144 A in the
minimum |i-j| >= 3 CA-CA distance; S is +0.060 and +0.070. Same sign as S25's +1.10 A expansion
of AMBER's top-75, smaller because the distogram's order is kept among the survivors.

WHERE THE SINGULARITY LIVES (the optional Part IV measurement, `singularity` block; every top-75
member rebuilt with all heavy atoms through the reference builder `sidechains.py`, no OpenMM):

    closest heavy-atom contact, residue separation >= 2
      all 9,450 members:              bb-bb 1,054   bb-sc 5,300   sc-sc 3,096
      the 5,057 members above 1e4:    bb-bb   163   bb-sc 2,627   sc-sc 2,267    -> 96.8% side-chain-involving
    members per target with a heavy-atom pair closer than 2.0 A:   all-atom 40.7 of 75 (SE 1.9)
                                                                    backbone+CB 2.61 of 75 (SE 0.34)   [S19 section 3.1: 2.66, reproduced]
    Spearman(e_amber, minimum heavy-atom distance) within a top-75: mean -0.743, median -0.803
    minimum heavy-atom distance: members above 1e4, 1.48 A; members below 1e4, 2.31 A

The AMBER single point on the top-75 is, to rho -0.74, the minimum heavy-atom distance of the
rebuild, and 96.8% of the rebuilds it condemns are condemned by a contact involving a side chain
placed by the deterministic builder (fixed chi1, no rotamer scan). The physically impossible
class on the backbone is 2.6 per 75 (S19); the class the 1e4 reject removes is 40 per 75.
Stated before any RMSD is read: the steric reject at the primary threshold is a reject of the
builder's side-chain placement, not of the pool's backbones. This is Part A of
`s26/IDEA_rotamer_relief.md` and it survives.

---

## L24 -- C3 NATIVE-FREE PART: THE PRODUCTION RELAXATION MOVES THE CA TRACE 0.220 A RMS; 58.7% OF BUILT CHAINS START ABOVE 1e4 kcal/mol; 125 OF 126 CONVERGE (2026-09-13, PH)

`s26/ph_c3.py nativefree`, `s26/results/ph_c3_nativefree.json` (complete 126/126), job
`s26/jobs_done/ph_c3_nativefree.json` (exit 0, 5 s; the RSS sampler polls every 5 s and
under-reads a 5 s process, so no memory number is claimed). Read from
`bench_results/cache/1fc9f2dcf489e2fb` with every RMSD key stripped; nothing here reads a native.

    AMBER displacement of the built chain, per-atom RMS after superposition   0.220 A (SE 0.008, median 0.197, range 0.103 to 0.591)
    `amber_moved` (restraint RMSD on N/CA/C, unsuperposed)                     0.233 A
    the built chain's own AMBER energy e0 before relaxation                    median 8.6e4 kcal/mol, min -473, max 1.3e14; 58.7% above 1e4
    relaxed energy e1                                                          mean -560 (SE 30), max +1262 (9KAR)
    converged (e1 <= CONVERGE_MAX_KCAL = 1000)                                 125 / 126 (9KAR fails)
    bond + angle strain after                                                  mean 60 kcal/mol
    virtual CA-CA bond: built 3.80395 (exact, sd 1e-16) -> relaxed             mean 3.867, min 3.12 (1M02), max 5.38 (2BP4, e1 +845); 4.86 on 9KAR
    targets with a relaxed CA-CA outside [3.6, 4.0]                            6 / 126
    Rg change                                                                  +0.046 A

Two things the validity story must carry: the emission's own strain census (58.7% of built
chains sit above 1e4 kcal/mol before relaxation, the same fraction as the pool), and the two
targets where the relaxation itself breaks a virtual bond (2BP4, 9KAR). Derived prediction for
stage 1, registered in `s26/PREREG_c3_control.md` addendum 1 before the gate: an orthogonal move
of 0.220 A on a 3.21 A chain costs about m^2 / (2 RMSD) = 0.0075 A by the S16 identity; the
production step costs +0.0207; so AMBER is predicted WORSE than the matched random control by
about +0.013 A. Measured after the gate.

---
## L25 -- LANE W (WILDCARD) SPAWNED; THE 2/60 PROXY BOUND IS ITS FIRST MANDATORY DIRECTION (2026-09-13, coordinator)

Six lanes active (E, I, Q, P, PH, W; `s26/lanes.json`). Brief `s26/briefs/W.md`. Of the mandatory
tournament directions, test-time window ensembling is already claimed by lane P
(`s26/IDEA_window_ensembling.md`) and the AMBER reject and cis census by PH (L22, L23); the 2/60
benchmark self-copy leak bounded from a dev-set proxy without opening the benchmark is unclaimed
and goes to W, building on lane I's L15 audit (the four dev self-copies and their carriers) and
S10-4's +0.0004 A dev price. W also owes at least three ideas nobody else has proposed and takes
the highest-ranked orphaned survivor once the Adversary ranks the tournament.

---


## L26 -- LANE P: THE C2 TRAINING CHAIN IS RUNNING, ONE GOVERNOR JOB AT A TIME, IN THE ORDER FIXED BY PREREG_C2 ADDENDUM 2; THE C3 HAND-OFF FORMAT IS IMPLEMENTED (2026-09-13, lane P)

Resumed 08:37 after the 00:46 session cut; no training job had been launched before the cut
(`s26/jobs_done/` held only the probes). Per L16b: order noesm -> conly -> pca32 -> wide -> pca32f
-> pca128 -> featurise-esm8m -> esm8m -> raw, fixed in `s26/PREREG_C2.md` addendum 2 BEFORE the
first job; `s26/p_train_chain.sh` runs them sequentially through `s26/jobrun.py` (est-ram 1.5 GB,
raw 1.8, NT = 2), one `.pt` per fold as the checkpoint, a second `_p2` pass resuming anything
the governor kills. `p_train_noesm` registered 08:52. No `eval`/`report` before "PHASE 0 SIGNED
OFF" (refused in code). Expected machine time 6-10 h (raw last: 5175 x 384 first layer).
While it runs: `s26/PREREG_C3.md` addendum 2 (not Proposal C's C3; C3 is lane PH's,
`s26/PREREG_c3_control.md`), `s26/IDEA_amber_prior_partner.md` (H_C3a entered in the tournament;
H_C3b withdrawn by S8-14's arithmetic), the C1 closure reproduction written into
`s26/agentP_FINDINGS.md` section 7 (S12 learning curve 3.0433 -> 3.0258 real vs 2.5342 leaked at
n = 8; S17 band best 2.6087 vs random 3.4676, all to the third decimal), `s26/p_stats.py`
(nested ridge, selftest OK: planted 0.866 vs null p95 0.558), `s26/p_b3.py`, `s26/p_c4.py`,
`s26/p_c5.py` (native-free feature builds and synthetic selftests queued under jobrun; their
`run` commands are gated). The C3 stage-2 hand-off format (`s26/results/p_best_rung_chains.json`:
rows keyed by pdb with phi/psi in RADIANS as `s12.instrument.project` emits them at lam = 0.3,
`ca`, `rmsd_arm`; `ST.save_atomic` with complete_keys and n_expected = 126) is implemented in
`s26/p_deliver.py`, a separate module so that `p_ladder.py` is not edited while its jobs run; it
re-projects the chosen rung and asserts equality with the eval JSON before writing.

## L27 -- LANE Q, A2: THE DEPLOYED ANSATZ'S DYNAMICAL LIE ALGEBRA IS THE FULL so(2^n) FROM DEPTH 2 AT n = 7; MY PRE-REGISTERED PREDICTION OF A PROPER SUBALGEBRA IS FALSIFIED; THE TANG MINIMAL POOLS GENERATE so(2^(n-1)+1) (2026-09-13, lane Q)

`s26/q_dla.py` -> `s26/results/q_dla.json` (complete; `s26/jobs_done/a2_dla.json`: 85.3 s, peak
RSS 0.479 GB). Pre-registered in `s26/PREREG_A2.md` before the run. Property measurement, no
native, no score: the closure is computed exactly on Pauli strings as a set and cross-checked
by dense SVD (rtol 1e-10) at n = 4, 5: 12 of 12 cells agree; the two conjugation conventions
agree on 12 of 12. Every closure lies in the odd-Y (real) set.

    fixed ansatz (RY / CNOT chain + ring), dim(DLA) by width n and depth L
      L = 1                       n            (abelian) at every n = 4..11
      n = 4, 5, 7, 8, 10, 11      so(2^n)      from L = 2 on   (8128 at the deployed n = 7)
      n = 6                       510 / 1023 / 2016 = so(64)  at L = 2 / 3 / 4
      n = 9                       32766 / 65535 / 130816 = so(512)  at L = 2 / 3 / 4
    pools V and G (Tang 2021, 2n-2 strings)   36, 136, 528, 2080, 8256, 32896  at n = 4..9
                                              = dim so(2^(n-1) + 1) at all six widths
    pool L2 (all 1-, 2-local odd-Y strings)   so(2^n) at n = 4..9
    ADAPT-selected sets, ideal ladder, n = 7  alpha = 1: abelian (7) at every step, both pools,
                                              both optimisers; L-BFGS stops at P = 7 with no
                                              operator selected.  alpha = 0.25: V 7 -> 16;
                                              L2 7 -> 1025 at P = 21 (12.6% of so(128)).

H2b of the prereg said dim(DLA) at (n = 7, L = 3) is below 8128, guess 4095. It is 8128, the
full so(128), already at L = 2. Falsified; the measured value stands. H2a (depth 1 abelian),
H2c (pools: 2080, 8256, 32896 predicted and measured), H2d (odd-Y) and H2e (ADAPT sets abelian
at alpha = 1) held. The only widths where depth 2 is not enough are n = 6 and n = 9, whose
dimensions at L = 2 and 3 equal 2 dim su(2^(n-2)) and dim su(2^(n-1)); that n divisible by 3
is the condition is an observation from two cases and is labelled HYPOTHESIS, not explained.

What it means for the S25 slopes (`s25/results/q_plateau.json`): the algebra at the deployed
cell is maximal, so nothing in the algebra protects the ansatz from an exponential plateau.
Larocca et al. 2022 / Ragone et al. 2024 give Var ~ 1/dim(g) once the circuit is a 2-design
over exp(g); dim(g) = 8128 at n = 7 and grows as 4^n / 2. S13 measured the decay base
approaching 0.504 per qubit at depth 8 (`s13/results/geo_kernel.json`), the 2-design rate.
The S25 result "no exponential plateau at n = 4..13" is therefore a statement about depth 3
(P = 3n against dim so(2^n)), and must be quoted with "at depth 3" attached. It is not
evidence of a favourable algebra, and Cerezo et al. 2025 does not apply through the
small-DLA route; the circuit is simulable because n = 7, not because of its structure
(S21 L4).

The pools: a "complete" pool in Tang's sense (overlap-matrix rank 2^n - 1) generates
so(2^(n-1) + 1), which acts transitively on the real sphere and has about a quarter of the
dimension of so(2^n). Completeness is weaker than controllability; the ADAPT arm using pool V
in A1 is therefore restricted to that subalgebra by construction, and the L2 arm is not.
## L28 -- PHASE 0 EXAMINATION LANDED: REPRODUCTION EXACT; FOUR DOCUMENT-ONLY NUMBERS; THE SEVEN-CONFIGURATION SUITE IS QUOTED ON THE POINT-CLOUD BASIS; ONE RULE-1 INCIDENT (2026-09-13, lane E)

`s26/EXAMINATION.md` (sections A to I), `s26/BRIEF.md`, `s26/agentE_FINDINGS.md`; scripts and
artefacts in commit `e3570aed` (`s26/e_module_map.py`, `e_hashes.py`, `e_reproduce.py`,
`e_trace.py`, `e_claims.py`; `s26/results/module_map.json`, `pinned_hashes.json`,
`e_reproduce.json`, `e_trace_1S9Z.json`, `e_trace_9KAR.json`, `claim_search.json`, `.txt`).

1. Reproduction (Part 2.3): fresh vs stored means 3.048338093879532 / 3.2040761603809194 /
   3.2147651542109985 / 3.2354598538973844, max per-target disagreement 0.0 on all four bases,
   T030 = 1S9Z 0.18198112330908295 (L16b records the same). Fresh `run_target` + `label` on
   1S9Z and 9KAR: every arm 0.0, `sub` identical, emitted `ca` identical. Every dataflow stage
   of both traces matches the record at 0.0; the Hamiltonian is the rank ladder to 0.394% of
   range (S25's worst case 1.18%, `s25/results/q_gibbs.json`).
2. Claim ledger (EXAMINATION C, 35 claims): 30 sourced to an artefact leaf, a passing test or a
   stated derivation from stored rows. **UNSOURCED (document-only):** 36.1 / 36.4 deg phi MAE
   (C26), the +0.0004 / +0.0030 identity-leak prices (C27; the cited `s10/idaudit_*.json` do
   not exist), the |z_moment| triple 0.7529 / 0.8013 / 0.1127 (C34), 355 passed / 13 skipped
   (C35; superseded by `s26/TEST_RUN.md`: 369 / 356 / 13), and the 0.524 sampled tail-only
   cosine inside C24. The benchmark delta +0.0103 (C06) is named to `s9/final_report.json`,
   which this lane did not open; the two means beside it are asserted to 5e-4 by two passing tests.
3. Basis: the seven-configuration suite (3.058 ... 3.881), the random-75 null (3.4251) and the
   Legacy +0.330 / AMBER +0.455 verdicts are point-cloud numbers (`s25/results/phys_suite.json`,
   `results/summary/leaderboard.json :: rows[*]/mean_secondary`); the state brief (lines 57,
   190-191) and `results/summary/professor_brief.md:115-125` quote them beside the built-chain
   3.2148 without saying so. Built-chain means of the same rows: 3.2187 / 3.3100 / 3.3732 /
   3.4221 / 3.8248 / 3.8844 / 4.1015 (`leaderboard.json :: rows[*]/mean`).
4. For the Adversary: `VQE_LFO` runs alpha = 1.0 on folds 0, 3, 4 (78 of 126 targets, no tail
   constraint; `s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint`
   0.6190476190476191); 9KAR's production relaxation ends at 1262.4 kcal/mol, above
   `core.amber.CONVERGE_MAX_KCAL` (L24: 125 of 126 converge) and is inside the 3.2355 mean;
   2BP4 converges by that gate (e1 845.3) but carries bond+angle strain 1172.7, above the S8
   strain-rejection rule; `bench_results/cache/1fc9f2dcf489e2fb/` is gitignored (EXAMINATION G).
5. Incident (Rule 1): the first run of `s26/e_claims.py` parsed `s9/final_report.json` before the
   benchmark exclusion existed and printed one aggregate leaf, `dist/shipped/mean` =
   2.9507235775391263, already published in `README.md:14`. No per-target benchmark value was
   printed or read. The exclusion (`SKIP_FILES`, `SKIP_SUBSTR`) is in the committed script and
   the artefacts were regenerated with it.
6. Provenance note: the traces and the reproduction ran at git `74e073e2` (dirty: lane I's
   later-committed `37bddbbb` edits to `core/data.py` and four unused imports were in the tree).
   The identity path they touch is not on `run_target`; the reproduction costs 5 s and should be
   re-run on HEAD by the Adversary (`python s26/e_reproduce.py`).

---

## L29 -- ADVERSARY SPAWNED ON THE EXAMINATION; THE RULE-1 INCIDENT LOGGED; THE POINT-CLOUD BASIS OF THE SUITE NOTED (2026-09-13, coordinator)

Lane A (Adversary) spawned at 08:55 with brief `s26/briefs/A.md` on `s26/EXAMINATION.md`,
`s26/BRIEF.md` and `s26/agentE_FINDINGS.md` (L28). Phase 0 remains closed until its audit is
clean or every material finding is fixed and re-checked. Six lanes active (E, Q, P, PH, W, A);
lane I finished (L21).

Three things from L28 the coordinator rules on now:

1. The Rule-1 incident (L28 item 5): the first run of `s26/e_claims.py` parsed
   `s9/final_report.json` and printed one aggregate, `dist/shipped/mean` = 2.9507235775391263,
   a number already published in `README.md:14` and the state brief. No per-target benchmark
   value was read or printed; the exclusion is in the committed script and the artefacts were
   regenerated with it. Logged as an incident with no consequence: the benchmark was not spent
   and nothing new was learned about it. Rule for the rest of the sprint: no script parses
   `s9/final_report.json`, `results/benchmark_manifest.json` or any file under `results/`
   naming a benchmark target; lane I's L18 read (sequences through `core.data.benchmark()` for
   a count) and this one are the only two benchmark-adjacent reads of S26.
2. The seven-configuration suite, the random-75 null and the Legacy +0.330 / AMBER +0.455
   verdicts are point-cloud numbers (L28 item 3; `s25/results/phys_suite.json`). The state brief
   and `results/summary/professor_brief.md` quote them beside the built-chain 3.2148 without
   naming the basis. Every S26 deliverable that quotes them (the report, the deck, the proposal
   verdicts) names the basis on both sides; the built-chain means of the same rows are in
   `results/summary/leaderboard.json :: rows[*]/mean` (distogram 3.2187 ... AMBER 4.1015).
3. Four presenter-facing numbers are document-only (C26 36.1/36.4 deg, C27 +0.0004/+0.0030,
   C34 the |z_moment| triple, C35 355/13). The Adversary decides whether each is material; the
   coordinator's default is that a document-only number is not put on a slide unless a lane
   re-derives it from cached rows this sprint (C26 from `s13/cache/tors_rows.npz` if it holds
   the per-residue predictions; C35 is superseded by `s26/TEST_RUN.md` and is simply replaced).

---


## L30 -- LANE W: THE 2/60 PROXY BOUND IS PRE-REGISTERED; NATIVE-FREE CENSUS: TWO OF THE FOUR DEV SELF-WINDOWS NEVER REACH THE EMISSION; EVERY TARGET HAS A BLOSUM TIE AT THE POOL BOUNDARY (2026-09-13, W)

`s26/PREREG_selfcopy_bound.md` (written before any result), `s26/IDEA_selfcopy_proxy_bound.md`,
code `s26/w_selfcopy.py`, synthetic tests `s26/w_selfcopy_test.py` (job `w_selfcopy_test2`, ALL
OK, 5 s; the tests caught one definitional error in the `sel` bound before any real run). Nothing
here reads a native, an RMSD to a native, a benchmark sequence, PDB or the manifest; every
native-free command overwrites `rr` and `nat_ca` with NaN on load and asserts its outputs finite.
The benchmark facts used are the record's only: 2/60 (S24 L4, L15) and the mechanism (S24 L4).

Design, in short: the leak has two channels, (A) the carrier's self-window at BLOSUM rank 0 in
the K = 500 pool, (B) the carrier's native distances in the fold model's training labels. Both
are live, in the same direction, on 1CEK / 2FBU / 2P5H / 6B9K. Part A drops the identity-1.0
window and refills the pool (S10-4's operator, exact copies only) on the 4, with the other 122
targets' BLOSUM rank-0 window dropped as the matched control population, the >= 0.6 variant
beside it (the C27 re-derivation on the production basis, L28), and an ORACLE insertion of the
withheld same-fold carriers on 8 more targets. Part B retrains the three fold models with the
carrier removed (through `s26/p_ladder.py`'s pca32 path; controls: two dev chains per fold).
Part C, the n = 126 envelope: the four pinned fold models that trained on the target's OWN
native, against the clean one, paired, fold-clustered (ORACLE). Part D: the bound
(2/60) x max per-target effect, and (2/60) x the fold-CI limit of the envelope, judged against
0.017 A (one tenth of the benchmark CI half-width). Before any RMSD is read, every part stores
its emissions and the TRIANGLE BOUND: Kabsch CA-RMSD is a metric, so RMSD(leaked emission,
un-leaked emission) bounds |change in RMSD-to-native| without the native. Assumptions A1 to A5
are in the PREREG.

Native-free census (`s26/results/w_selfcopy_census.json`, job `w_selfcopy_census`, 20 s, 0.064 GB):

    target  n  fold  exact self-window: universe idx / BLOSUM rank / pool pos / in shipped top-75   >= 0.6 windows: universe / pool / top-75
    1CEK   13   2          7 / 0 / 0 / NO                                                            8 / 3 / 0
    2FBU   12   4       1602 / 0 / 0 / NO                                                            4 / 1 / 0
    2P5H    9   4       3316 / 0 / 0 / YES                                                           5 / 1 / 1
    6B9K   10   0       1082 / 0 / 0 / YES                                                           7 / 1 / 1

So channel A cannot touch the emitted structure of 1CEK or 2FBU at all (the self-window is
filtered out by the distogram score); on 2P5H and 6B9K it is one member of the 75 averaged.
21/126 targets carry a >= 0.6 window somewhere in the universe (S10-4's 21 reproduced; its
13-in-pool count is checked when the retrieval run completes). Only 31% of targets keep their
BLOSUM rank-0 window in the top-75. Probes: retrieval on 1CEK (`w_selfcopy_retrieval_probe`,
15 s, 0.104 GB): with and without the self-window the cloud, chain and lam = 0 chain are
IDENTICAL (triangle bound 0.000) and the refill window does not enter the top-75; envelope on
1CEK (`w_selfcopy_envelope_probe`, 25 s, 0.298 GB): production gate true (top-75 equals the
cache's `sub`), the four leaked models keep 81 to 87% of the top-75 and move the built chain by
0.079 / 0.124 / 0.137 / 0.122 A (triangle bounds, models 1 / 0 / 3 / 4). Full native-free runs
`w_selfcopy_retrieval` (est 0.5 GB, ~20 min) and `w_selfcopy_envelope` (est 0.6 GB, ~40 min)
are registered; the gated `endpoint`, `floor` and `report` wait for sign-off.

A second census fact, recorded because it is the basis of `s26/IDEA_tiebreak_noise_floor.md`:
every one of the 126 targets has a BLOSUM tie at the K = 500 boundary; the median tie class at
the boundary score holds 115 windows, about 57 inside the pool and 60 outside, so ~11% of every
pool is chosen by corpus order. S17 measured only the ORACLE pool best under random tie-breaks
(sd 0.018); the production endpoint's floor is unmeasured. Three own ideas filed:
`IDEA_tiebreak_noise_floor.md`, `IDEA_conformational_identity_floor.md` (Part E of the
PREREG: the same sequence in a different deposit), `IDEA_window_provenance.md` (census first,
plausibility 0.15).

Question for the coordinator (rule 14; not run until answered): may Part B's 10 retrains (tag
CPU, est-ram 1.5 GB, 456 s each, one at a time, order fixed in the PREREG section 8) run before
sign-off under L16b's conditions, and may they interleave with lane P's chain (two training
jobs resident, ~2.5 GB, against 4.7 GB available at 08:54)? Default if unanswered: the
reference is lane P's `pca32` fold models and Part B trains after sign-off. Noted and NOT done:
the triangle bound is computable on the benchmark itself with no native or RMSD read, but it
would need the two benchmark sequences and universes, which L18b/L29 close for this sprint.

## L31 -- EXAMINATION AUDIT: REPRODUCTION EXACT ON HEAD; ONE MATERIAL FINDING, THE TEST-SUITE TOTAL (369/356 QUOTED, 370/357 IN THE ARTEFACT); C26 RE-DERIVED; C27 IS IN GIT HISTORY (2026-09-13, A)

`s26/EXAMINATION_AUDIT.md`, sections 1 to 9. (1) `s26/e_reproduce.py` re-run on HEAD `76287a33`
through `s26/a_reproduce_head.py` (the Examiner's `main()` unchanged, the write redirected so
lane E's artefact is not overwritten; job `a_reproduce_head`, exit 0, 5.0 s;
`s26/results/a_reproduce_head.json`): all four means and all 126 rows equal
`s26/results/e_reproduce.json` at 0.0; T030 = 1S9Z `rmsd_arm` 0.18198112330908295. (2) 29
claim-ledger leaves opened (seed-26 sample C02, C05, C08, C15, C17, C19, C21, C23, C28, C32, plus
every UNSOURCED and markdown-sourced entry); every sourced value sits at its cited key at stored
precision. (3) The trace is faithful to `core/pipeline.py:836-848` (selector), `:934` (average)
and `core/predict.py:416, 422` (score weight `1/(sd+0.5)`, `_risk`); the identical pool index
449 on both traced targets is verified in `bench_results/cache/1fc9f2dcf489e2fb/{1S9Z,9KAR}.json
:: sub` (positions 44 and 34). (4) The junit files agree with `s26/TEST_RUN.md` and
`test_run.json`: 370 / 357 / 0 / 0 / 13, 0 memory-guard skips, 11 `VERIFY_SLOW` and 2
absent-artefact skips. (5) Every defect file:line resolves at HEAD. (6) Four hashes recomputed,
four equal. (7) Three sizes as listed, all ignored. (8) Rule 1: the grep over `s26/*.py` hits
only `e_claims.py` (the SKIP lists, committed `e3570aed`), `e_hashes.py` (bytes) and
`i_identity_audit.py` (L18's sequence count). L28 item 3: `s25/results/phys_suite.json :: basis`
= "point_cloud"; the artefact declares its own basis.

MATERIAL A1. `s26/EXAMINATION.md` section D ("Combined: 369 tests, 356 passed"), row C35,
`s26/agentE_FINDINGS.md` X1 and L28 item 2 quote 369 / 356 / 13 for `s26/TEST_RUN.md`, which
holds 370 / 357 / 13 (`pytest_core_post.xml` 351 / 338 / 13 plus `pytest_amber` 16 and
`pytest_amber_frame` 3). The count is the 00:41 render, before `pytest_core_post` finished at
00:45:56; section D lists four junit files and omits `pytest_core_post.xml`. L20 and the state
brief 6.1 already carry the right number; nothing in Phase 1 depends on it. Fix: lane E, or the
coordinator if E has finished, appends a correction to `s26/EXAMINATION.md` (D and C35) and
`s26/agentE_FINDINGS.md` and posts a ledger entry naming L28 item 2; A re-checks on the append.

MINOR. C26 (36.1 / 36.4 deg) is now DERIVED: `s13/cache/tors_rows.npz` arms `p_grid` / `n_marg`,
key `err_phi`, pooled over 1,507 residues: 36.133 / 36.416 deg (psi 62.355 / 72.772), matching
`s13/SPRINT13_DOSSIER.md:261-263`. C27's artefacts are not absent: `s10/idaudit_price.json` and
four siblings are in git history at `5fa05cd` (commit `49fc708` exists); not opened, because the
price stage covers the sealed benchmark; the +0.0030 half stays historical, the +0.0004 dev half
is lane W's to re-derive (L25). C34 confirmed document-only (a grep of `s25/` hits only the two
md files). C35's replacement is 370 / 357 / 13. C24's 0.524: `stage_cvarcheck` exists
(`s8/integrate.py:1162`) but its JSON and console log do not; quote 0.655634
(`tests/test_quantum.py:59`) or 0.566586 (`s25/results/q_verify.json`) instead.

Question for the coordinator: `tests/test_pipeline.py:338-347` parses `s9/final_report.json`
aggregates, so any further full `pytest tests/` in S26 does so again; accept, or `--deselect` it
for the rest of the sprint.

AUDIT: MATERIAL FINDINGS 1, LISTED ABOVE

---
## L32 -- L28 ITEM 2 CORRECTED: THE GOVERNED TEST RUN IS 370 / 357 / 13; THE test_pipeline BENCHMARK-AGGREGATE READ IS ACCEPTED; C26 IS DERIVED (2026-09-13, coordinator)

1. The Adversary's one MATERIAL finding (L31 item A1) is fixed by appended addenda in
   `s26/EXAMINATION.md` (sections D and C35) and `s26/agentE_FINDINGS.md` (X1): the governed
   suite is 370 tests, 357 passed, 0 failed, 0 errors, 13 skipped, 0 memory-guard skips
   (`s26/results/test_run.json :: combined`; junit files `pytest_core_post.xml` 351/338/13,
   `pytest_amber.xml` 16/16, `pytest_amber_frame.xml` 3/3). L28 item 2's "369 / 356 / 13" was
   the 00:41 render before `pytest_core_post` finished. Lane E's lines are left as written.
2. Ruling on the Adversary's question: `tests/test_pipeline.py:338-347` reads two published
   aggregates of `s9/final_report.json` (the benchmark means 2.9507 and 2.9610 already in
   `README.md`) and asserts them to 5e-4. It is the repository's own regression test, it
   predates S26, it reads no per-target value, and a pass or fail carries no new information
   about the benchmark. Accepted; it is not deselected, so the S26 test record stays whole.
   No S26 lane adds any benchmark read of its own (L29 item 1 stands).
3. C26 (phi MAE 36.1 vs 36.4 deg) is DERIVED by the Adversary from `s13/cache/tors_rows.npz`
   (arms `p_grid` / `n_marg`, key `err_phi`, 1,507 residues: 36.133 / 36.416 deg). The Adversary
   saves that derivation as `s26/results/a_c26_phi_mae.json` with provenance; with it the number
   may appear on a slide with that artefact path. C27's +0.0004 dev half is lane W's to
   re-derive (L25, L30); C34 and C24's 0.524 stay document-only and off every slide.

---

## L33 -- PHASE 0 SIGNED OFF (2026-09-13 09:05, coordinator)

The Adversary's audit (`s26/EXAMINATION_AUDIT.md`, L31) found the reproduction exact on HEAD on
all four bases and all 126 rows, 29 claim leaves at stored precision, the trace faithful to the
code, the hashes and sizes matching, Rule 1 closed, and one MATERIAL documentary finding, which
L32 fixes with appended addenda; the Adversary re-checks the append on its next turn and may
veto by ledger entry, in which case this sign-off is retracted by a new entry and every result
read in between is provisional. Phase 1 begins: lanes Q (A1, A3), P (ladder evaluation, B3, C4,
C5), PH (cis floor, C3 stage 1, steric reject) and W (the 2/60 proxy bound) may read RMSDs to
natives under their pre-registrations. Standing: every positive result is attacked by the
Adversary before it enters the closed/open tables; basis named on both sides of every contrast;
the phase gate string in this heading is the one the lanes' code checks.

---

## L34 -- ADVERSARY RE-CHECK OF L32: THE ADDENDA STAND; THE L33 SIGN-OFF IS NOT VETOED (2026-09-13, A)

Re-read at commit `195cf93e`: the appended addenda in `s26/EXAMINATION.md` (after section I) and `s26/agentE_FINDINGS.md` (after the unsourced list) both state 370 tests, 357 passed, 0 failed, 0 errors, 13 skipped, 0 memory-guard skips with the artefact paths (`s26/results/test_run.json :: combined`, `pytest_core_post.xml` 351/338/13, `pytest_amber.xml` 16/16, `pytest_amber_frame.xml` 3/3), and L32 item 1 matches them; the original lines are left as written. L31 item A1: STANDS as fixed. No veto.

---

## L35 -- LANE Q, A4: ADAPT-GROWN CIRCUITS ON THE DEPLOYED HAMILTONIAN ARE PRODUCT CIRCUITS AT alpha = 1 (NO DECAY, NOTHING TO TRAIN) AND DECAY LIKE THE FIXED ANSATZ AT alpha = 0.25; THE S25 n = 7 ROWS REPRODUCE EXACTLY (2026-09-13, lane Q)

`s26/q_var.py` -> `s26/results/q_var.json` (complete); `s26/jobs_done/a4_var.json`: 2,765 s,
peak RSS 0.08 GB. Pre-registered in `s26/PREREG_A4.md` before the run. Property measurement:
energies are the deployed SHAPE (standardised ranks 1..2^n), no target, no native. Figure:
`s26/figures/a4_variance_slopes.png` (190 dpi, white), from the JSON.

Gate: `s25.q_plateau.measure(7, 3, alpha, T, 250, seed=1007)` reproduces the n = 7 row of
all five cells of `s25/results/q_plateau.json` at relative deviation 0.0 (bit-identical).

    fitted log2 Var[dF/dtheta_0] per qubit, n = 4..13, theta ~ N(0, 0.6^2), same draws as S25
    (grown: rows with P = 3n only, count in brackets; growth by Adam best-iterate, 50 steps,
     eps 1e-3, pool V = Tang 2n-2, pool L2 = 2-local odd-Y)
      cell                    fixed (S25)   grown V          grown L2
      alpha=1,    T=0         -0.649        -0.079 (7)       +0.006 (7)
      alpha=0.25, T=0         -0.252        degenerate (0)   degenerate (0)
      alpha=0.10, T=0         -0.047        degenerate (0)   degenerate (0)
      alpha=1,    T=0.3       -0.311        +0.035 (5)       -0.008 (7)
      alpha=0.25, T=0.3       -0.243        -0.246 (6)       -0.302 (7)
    grown/fixed variance ratio at matched n: 1.63 (n = 4, alpha=1 T=0.3, both pools) and
    2.5 to 370 everywhere else (31 of 32 matched rows above 2).

What the grown circuits are. At alpha = 1 (T = 0 or 0.3) they hold 1 to 3 DISTINCT operators;
the other 12 to 25 selections are consecutive repeats of the same single-qubit Y (a repeated
rotation merges with the previous one): they are product circuits with the parameter count
of the fixed ansatz, which is why their variance does not decay with n. At alpha = 0.25,
T = 0.3, pool L2 grows 4 to 21 distinct 2-local strings and its decay, -0.302 per qubit, is
the fixed ansatz's -0.243 within the sampling error of a 7-point slope (relative SE of a
variance 9 to 16%). At T = 0 with alpha < 1 growth stops at P = n on every width (the
collapse; all first-order gradients vanish at a basis state, the lemma in s26/q_adapt.py) and
the variance is 0 to 2.4e-2; those cells are degenerate by construction, as pre-registered.

Predictions (PREREG_A4): H4a (T = 0.3 grown slopes in [-0.3, 0]) held for three of four
(V alpha=1 +0.035 and L2 alpha=0.25 -0.302 sit just outside by 0.035 and 0.002, inside the
error of the fit); the falsifier (below -0.5) did not fire. H4b (ratio > 2 at every matched
n) failed on 1 of 32 rows (n = 4, 1.63) and held on 31. H4c held on P < 3n (14 of 14) and
failed on "Var < 1e-3" for 3 of 7 rows at alpha = 0.25 (max 2.4e-2 at n = 4). H4d held
exactly. All four recorded as measured, none softened.

What it means, in the record's own terms. A large gradient from a product circuit is not
trainability; it is the trivial regime (rule 10's mirror: never call a large gradient a
merit). The one non-trivial grown family (alpha = 0.25, L2) decays like the fixed ansatz it
would replace. A4 gives Proposal A no width-scaling argument.
## L36 -- jobrun v2.1: THE LAUNCH CAP IS RE-CHECKED AFTER A JITTER; SIX JOBS HAD LAUNCHED AGAINST A CAP OF FOUR (2026-09-13, coordinator)

`s26/governor_state.json` at 09:28:11 showed six registered jobs (p_train_conly, two W census
passes, a1_build_s0 and s1, p_eval_shipped) against jobrun's cap of four: when A4 finished,
three waiting jobrun processes read the same 5 s snapshot with three jobs and all launched.
The governor's band held (smoothed CPU 85.7%, RAM 76.9%; a1_build_s0 suspended as the newest
job, resumed when CPU fell), so nothing was lost, but the cap was not a cap. `s26/jobrun.py`
now counts live registrations in `s26/jobs/` directly instead of the snapshot, sleeps a random
0.2 to 3 s after passing the gate and re-checks the count before registering. Waiters started
before this edit run the old code until they launch; new waiters use v2.1. No production module
touched; `s26/governor.py` unchanged.

---

## L37 -- CAMPAIGN PAUSED AT THE USER'S REQUEST (2026-09-13 09:30); STATE SAVED; GOVERNED JOBS KEEP RUNNING (2026-09-13, coordinator)

The user asked for a break with the usage limit at 96%. Every lane was told to commit its files,
append a STATUS line and end its turn; the coordinator commits whatever remains as a checkpoint.
Governed python jobs keep running on the box (they cost no API usage and checkpoint to disk):
at 2026-09-13T09:30:16: p_train_conly (P), w_selfcopy_retrieval (W), w_selfcopy_envelope (W), a1_build_s1 (Q), a1_build_s0 (Q, suspended), p_eval_shipped (P); RAM 77.1%, smoothed CPU 93.6%. The governor v2
(pid 25480) stays up. `s26/models/` (11 MB of ladder checkpoints, regenerable) is not committed.

Where the sprint stands: Phase 0 signed off (L33, audit L31/L34). Done: A2 (L27), A4 (L35), the
cis census (L22), the steric-reject census (L23), the C3 native-free part (L24), C1 (lane P
findings section 7), B1 infeasible (L13), all operational items and declared defects (L6 to L10,
L15, L16, L18 to L21), the report draft Parts I to VI, IX, Appendix A (`s26/REPORT.md`), the
Proposal B replacement outline draft, the trainability paper outline, A2/A4 figures. Running or
queued: A1 shards then A3 (Q), the C2 ladder training and evaluation (P), the cis floor, C3
stage 1 and the steric reject (PH), the 2/60 proxy bound (W), the Adversary's checks of L27,
L22 to L24, L30, L35 and the tournament ranking (A), the report's Appendix B check (E).

To resume: read this ledger from L33, `s26/STATUS.md`, `s26/jobs_done/` and `s26/results/`,
then continue each lane from its STATUS line; the lane briefs are in `s26/briefs/` and every
lane's findings file names what it did not do. No lane restarts finished work.

---

## L38 -- CIS FLOOR (ORACLE DIAGNOSTIC): THE IDEAL-TRANS REBUILD OF THE NATIVE'S OWN TORSIONS SITS 0.347 A FROM THE NATIVE (max 1.47), TRACKS OMEGA NON-PLANARITY AT rho 0.83, AND DOES NOT TRACK THE 0.166 A CHAIN COST (rho 0.08) (2026-09-13, PH)

`s26/ph_cis.py floor`, `s26/results/ph_cis_floor.json` (complete 126/126), job
`s26/jobs_done/ph_cis_floor.json` (exit 0, 10 s, peak RSS 0.098 GB, held 345 s at the job cap).
Pre-registered in `s26/PREREG_cis.md` part 2 and addendum 1. ORACLE DIAGNOSTIC: every number
reads the native. Bases, named: `floor_ca` is the CA-RMSD between the model-1 native CA trace
and the chain rebuilt from the native's OWN phi/psi at ideal trans geometry
(`core.geometry.build_backbone`), after superposition; `rebuild_bb` is the same on N/CA/C
(`core.data.Peptide.rebuild`, verified from `core/data.py:410-419`); `chain_cost` is the
production `rmsd_arm - rmsd_avg` (built chain minus point cloud, a difference across bases,
`bench_results/cache/1fc9f2dcf489e2fb`).

    floor_ca   (CA, ideal trans on the native's own torsions)   mean 0.347 A (SE 0.029), median 0.272, p90 0.752, max 1.474 (1ID6)
    rebuild_bb (N/CA/C, `Peptide.rebuild`)                         mean 0.340 A (SE 0.028), median 0.248, max 1.413
    chain_cost (production rmsd_arm - rmsd_avg)                    mean 0.166 A (SE 0.018), median 0.098
    targets with floor_ca > 0.5 A                                  32 / 126;  floor_ca > chain_cost on 86 / 126
    Spearman(floor_ca, max omega deviation)                        +0.828;  with mean omega deviation +0.833
    Spearman(chain_cost, max omega deviation)                      -0.036;  Spearman(chain_cost, floor_ca) +0.083
    cis targets                                                    0 (L22); the cis-vs-non-cis contrast is empty and is not printed
    9UV5 (the one bond beyond 30 deg)                              floor_ca 0.663, rebuild_bb 0.708, chain_cost 0.372

`ST.fmt`, verbatim (the two quantities are both per-target Angstroms but are different
objects; "BETTER" below means only that the production chain cost is smaller than the
own-torsion rebuild floor, not that anything improved):

    production chain cost (rmsd_arm - rmsd_avg) MINUS the ORACLE CA floor of ideal trans geometry on the native's own torsions
      a 0.1664 (med 0.0977)   b 0.3468 (med 0.2722)   n=126
      effect -0.1804   median -0.1261   SE 0.0333   MDE 0.0933   effect/MDE -1.93
      iid  CI95 [-0.2433, -0.1176]
      fold CI95 [-0.2266, -0.1022]   folds same sign 5/5   per-fold 0:-0.025 1:-0.224 2:-0.230 3:-0.223 4:-0.202
      86W/40L/0T   worst degradation +0.5161 (1RSW)   p90 +0.2720   power 1.00  Type-M 1.00
      concentration: drop-top10 -0.1100 vs uniform-effect null p10/p50/p90 -0.1523/-0.1116/-0.0702 -> pctile 0.518
      VERDICT: BETTER

Reading. (1) The constant omega does carry a representation cost on this instrument even
without a single cis bond: the ideal-trans rebuild of the native's own torsions misses the
native by 0.35 A on average and by more than 0.5 A on a quarter of the targets, and that miss is
the omega non-planarity to rho 0.83 (L22: 0.5% of bonds beyond 20 deg, mean deviation 1.9 deg;
small deviations accumulate along the chain). (2) It is NOT what the production projection
pays: the chain cost (0.166 A) does not correlate with the floor (rho +0.08) nor with the omega
deviation (rho -0.04). The projection's cost is a displacement effect of the operator (S16 L27:
the projection moves 0.813 A with a negative cosine), not a representation effect. (3)
`floor_ca` is an UPPER bound on the manifold floor, because the projection fits phi/psi to the
trace rather than rebuilding from the native's torsions; `s26/PREREG_cis.md` addendum 2
registers `floor2` (the native projected through the production projection) as the tight
number, 10 min CPU, to run after the reject jobs. Power: n = 126, SE 0.033, MDE 0.093 on the
contrast; a Spearman at n = 126 has SE about 0.09, so the two nulls (+0.08, -0.04) exclude
|rho| above about 0.25. Nothing here is a proposal; the two-bond-length projection stays worth
0.000 A on this instrument (L22), and the design note in `s26/agentPH_FINDINGS.md` 1.4 stands.

---

## L39 -- C3 STAGE 1: THE PRODUCTION RELAXATION IS WORSE THAN A RANDOM MOVE OF ITS OWN SIZE (+0.0111, fold CI [+0.0062, +0.0171]) AND WORSE THAN A MOVE TOWARD A POOL MEMBER (+0.0385); ITS DISPLACEMENT POINTS AWAY FROM THE NATIVE (cos -0.049); REFINE-WITH-PHYSICS IS A VALIDITY STEP, NOT AN ACCURACY STEP (2026-09-13, PH)

`s26/ph_c3.py stage1`, `s26/results/ph_c3_stage1.json` (complete 126/126), job
`s26/jobs_done/ph_c3_stage1.json` (exit 0, 10 s, peak RSS 0.041 GB, held 345 s at the job cap).
Pre-registered in `s26/PREREG_c3_control.md` sections 2 and 3 and addendum 1 (prediction:
AMBER worse than the matched random move by about +0.013 A, cosine negative). No AMBER compute:
the production coordinates `ca` and `amber_ca` from `bench_results/cache/1fc9f2dcf489e2fb`.

Bases, named on both sides. Arm: input the BUILT CHAIN (`rmsd_arm`, 3.2148), output the
RELAXED CHAIN (`rmsd_full`, 3.2355), both reproduced from the stored coordinates to 1e-6.
Controls: the built chain's CA trace displaced by a vector of the SAME per-atom RMS magnitude
as AMBER's displacement (0.220 A mean, measured after superposing `amber_ca` onto `ca`);
`rand` = S16's construction (isotropic Gaussian, six rigid components projected out on
`rigid_basis`, 16 draws, mean); `member` = the same magnitude along the straight line toward a
random member of the shipped K=500 pool superposed onto `ca` (16 draws, mean; 6 of 2,016 draws
overshot the member). A displaced CA trace is not an ideal-geometry chain: the controls match the
operator's SPACE (magnitude), not its validity, exactly as S16 did. Full-chain CA-RMSD to the
native, ORACLE evaluation of native-free operators. `ST.fmt`, verbatim, all seven contrasts:

    stage1 AMBER (relaxed) minus do-nothing (built chain)
      a 3.2355 (med 2.9757)   b 3.2148 (med 2.9661)   n=126
      effect +0.0207   median +0.0158   SE 0.0034   MDE 0.0096   effect/MDE +2.16
      iid  CI95 [+0.0141, +0.0274]
      fold CI95 [+0.0154, +0.0290]   folds same sign 5/5   per-fold 0:+0.021 1:+0.037 2:+0.017 3:+0.013 4:+0.017
      40W/86L/0T   worst degradation +0.1653 (2MID)   p90 +0.0640   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0258 vs uniform-effect null p10/p50/p90 +0.0213/+0.0256/+0.0302 -> pctile 0.529
      VERDICT: WORSE
    stage1 AMBER minus matched-magnitude RANDOM (16 draws)
      a 3.2355 (med 2.9757)   b 3.2244 (med 2.9669)   n=126
      effect +0.0111   median +0.0070   SE 0.0035   MDE 0.0099   effect/MDE +1.12
      iid  CI95 [+0.0044, +0.0178]
      fold CI95 [+0.0062, +0.0171]   folds same sign 5/5   per-fold 0:+0.011 1:+0.023 2:+0.009 3:+0.003 4:+0.010
      52W/74L/0T   worst degradation +0.1518 (2MID)   p90 +0.0539   power 0.88  Type-M 1.07
      concentration: drop-top10 +0.0168 vs uniform-effect null p10/p50/p90 +0.0122/+0.0168/+0.0216 -> pctile 0.497
      VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.07x]
    stage1 AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)
      a 3.2355 (med 2.9757)   b 3.1969 (med 2.9755)   n=126
      effect +0.0385   median +0.0338   SE 0.0061   MDE 0.0170   effect/MDE +2.27
      iid  CI95 [+0.0272, +0.0506]
      fold CI95 [+0.0298, +0.0479]   folds same sign 5/5   per-fold 0:+0.039 1:+0.052 2:+0.049 3:+0.029 4:+0.027
      29W/97L/0T   worst degradation +0.3497 (2BP4)   p90 +0.1092   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0485 vs uniform-effect null p10/p50/p90 +0.0407/+0.0481/+0.0561 -> pctile 0.521
      VERDICT: WORSE
    stage1 RANDOM minus do-nothing
      a 3.2244 (med 2.9669)   b 3.2148 (med 2.9661)   n=126
      effect +0.0096   median +0.0071   SE 0.0015   MDE 0.0042   effect/MDE +2.28
      iid  CI95 [+0.0070, +0.0127]
      fold CI95 [+0.0077, +0.0122]   folds same sign 5/5   per-fold 0:+0.010 1:+0.015 2:+0.008 3:+0.010 4:+0.007
      29W/97L/0T   worst degradation +0.1147 (1I93)   p90 +0.0276   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0114 vs uniform-effect null p10/p50/p90 +0.0094/+0.0113/+0.0135 -> pctile 0.528
      VERDICT: WORSE
    stage1 TOWARD-MEMBER minus do-nothing
      a 3.1969 (med 2.9755)   b 3.2148 (med 2.9661)   n=126
      effect -0.0178   median -0.0118   SE 0.0043   MDE 0.0120   effect/MDE -1.49
      iid  CI95 [-0.0266, -0.0096]
      fold CI95 [-0.0255, -0.0124]   folds same sign 5/5   per-fold 0:-0.018 1:-0.015 2:-0.031 3:-0.016 4:-0.010
      90W/36L/0T   worst degradation +0.1204 (8T62)   p90 +0.0274   power 0.99  Type-M 1.01
      concentration: drop-top10 -0.0088 vs uniform-effect null p10/p50/p90 -0.0141/-0.0092/-0.0040 -> pctile 0.543
      VERDICT: BETTER
    stage1 ORACLE cos: AMBER minus RANDOM
      a -0.0491 (med -0.0498)   b 0.0013 (med 0.0028)   n=126
      effect -0.0503   median -0.0509   SE 0.0154   MDE 0.0432   effect/MDE -1.17
      iid  CI95 [-0.0799, -0.0211]
      fold CI95 [-0.0735, -0.0253]   folds same sign 5/5   per-fold 0:-0.050 1:-0.088 2:-0.066 3:-0.001 4:-0.046
      74W/52L/0T   worst degradation +0.3845 (5MXS)   p90 +0.1647   power 0.90  Type-M 1.06
      concentration: drop-top10 -0.0225 vs uniform-effect null p10/p50/p90 -0.0433/-0.0232/-0.0028 -> pctile 0.517
      VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.06x]
    stage1 ORACLE cos: AMBER minus TOWARD-MEMBER
      a -0.0491 (med -0.0498)   b 0.1229 (med 0.1090)   n=126
      effect -0.1719   median -0.1733   SE 0.0232   MDE 0.0649   effect/MDE -2.65
      iid  CI95 [-0.2156, -0.1266]
      fold CI95 [-0.2104, -0.1388]   folds same sign 5/5   per-fold 0:-0.168 1:-0.190 2:-0.244 3:-0.121 4:-0.140
      94W/32L/0T   worst degradation +0.5848 (2NDN)   p90 +0.1612   power 1.00  Type-M 1.00
      concentration: drop-top10 -0.1320 vs uniform-effect null p10/p50/p90 -0.1626/-0.1327/-0.1020 -> pctile 0.512
      VERDICT: BETTER

    ORACLE cos(AMBER displacement, true residual)   mean -0.049 (SE 0.015), median -0.050, positive on 36.5% of targets
    orthogonal-move cost predicted by the S16 identity from the magnitude alone   +0.0107 (SE 0.0011), median +0.0079
    identity n RMSD^2 = |r|^2 - 2 v.r + |v|^2 as an upper bound on the superposed RMSD   max violation 0.0062 A (it is a bound, not an equality, under re-superposition)

Reading. (1) The production step reproduces exactly: +0.0207 [+0.0154, +0.0290], 2.16x MDE.
(2) About half of that cost is the SIZE of the move: a random displacement of the same 0.220 A
costs +0.0096, which is what the S16 identity predicts for a move orthogonal to the residual
(+0.0107). The other half is DIRECTION: AMBER is worse than its own random twin by +0.0111
(fold CI excluding zero, 5/5 folds, but 1.12x MDE, so the magnitude is an upper bound), and its
displacement points away from the native, cos -0.049 (36.5% positive), which reproduces S16 L27
to the fourth decimal (S16: cos -0.0521, AMBER minus random -0.0491; here AMBER minus random on
the cosine -0.0503). The registered prediction (+0.013, negative cosine) held. (3) A move of the
same size TOWARD A RANDOM POOL MEMBER improves the built chain by -0.0178 [-0.0255, -0.0124],
90W/36L, 1.49x MDE, 5/5 folds, with ORACLE cos +0.123. This is a zero-information control, not
a proposal, and it is PROVISIONAL until the replication registered in `s26/PREREG_c3_control.md`
addendum 2 lands (job `ph_c3_stage1_rep`: new seeds, reversed order). Its mechanism is on the
record: the projection moved the chain 0.166 A away from the point cloud with a negative cosine
(S16 L27), and the pool members surround the cloud, so any move back toward one of them
recovers part of that displacement. It says nothing about physics and everything about the
projection's cost; it will be handed to the Adversary as such. (4) Power: SEs 0.0015 to 0.0061,
MDEs 0.004 to 0.017; every contrast except AMBER-vs-random and the AMBER-vs-random cosine is
above 1.3x its MDE; those two are in the Type-M zone, direction measured, magnitude inflated
about 1.07x.

DECISION RULE (`s26/PREREG_c3_control.md` section 3): AMBER does not beat either matched
control; it is worse than both. "Refine with physics" is DROPPED as an accuracy step and KEPT
as a validity step. `s26/C3_RESULT.md` carries the sentence with the numbers for lane P and the
presentation. Stage 2 repeats this on lane P's best C2 rung when it is delivered.

---
## L40 -- GOVERNOR KILLED BY THE HOST AT 10:10 FOR LOW MEMORY; AT 19:10 THE BOX WAS AT 91.8% RAM; THE TWO REMAINING LANE-P JOBS WERE STOPPED BY THE COORDINATOR (2026-09-13 19:11, coordinator)

The harness stopped the governor v2 background process (pid 25480) at about 10:10 with the message
"the system is running low on memory" (`s26/governor.log` ends at 10:10:36 with RAM 87.5%, six
jobs registered). Every lane was then cut by the API session limit. At 19:10 the coordinator found:
RAM 91.1 to 91.8% used with 1.5 GB available, no governor alive, a1_build_s0/s1, ph_reject_cloud
and w_train_chain gone (their records, if any, are in `s26/jobs_done/`), `p_eval_shipped` suspended
since 09:25 (0.03 GB; suspended by the governor before it died, never resumed), and `p_train_raw`
running unsupervised since 16:50 at 1.11 GB (launched by lane P's chain driver through jobrun,
which treats a stale governor snapshot as "go"). With the box at the campaign ceiling and the host
already killing processes for memory, the coordinator terminated both lane-P jobs and the chain
driver (per-fold checkpoints exist under `s26/models/p_ladder/`; lane P resumes the `raw` rung and
the shipped evaluation from them) and restarted the governor. Nothing else of the campaign's is
running. Rule added for jobrun: a stale governor snapshot is NOT "go" for a job whose est-ram
exceeds 0.5 GB (implemented at the next resume; recorded here first).

---

## L41 -- THE BOX IS AT 84% WITH NOTHING OF OURS RUNNING; jobrun v2.2; WHAT FINISHED DURING THE PAUSE (2026-09-13 19:13, coordinator)

After L40's stops the box reads 84.2% RAM used (14.1 GB) with no campaign process alive: the
load is the user's own (twelve Chrome renderers at 0.3 to 0.8 GB each, VS Code, three claude
sessions). Available memory is 2.6 GB, so the campaign's working headroom under the 93% ceiling
is now about 1.4 GB, a third of what L0 measured. Every remaining job must be launched one at a
time with a measured peak below 1 GB until the user's load drops; the ladder's `raw` rung (peak
of its sibling `pca32f` 1.85 GB) does not fit and waits.

`s26/jobrun.py` v2.2: (i) a stale governor snapshot is no longer "go" for a job with est-ram
above 0.5 GB (L40); (ii) a job never starts unless the snapshot shows its estimate plus 0.5 GB
available. Lane P's shell chain drivers (`p_train_chain.sh`, `p_eval_chain.sh`) were terminated
because they re-launch the next rung whenever the previous job exits, including on a kill; on
resume lane P launches rungs singly.

Finished during the pause (`s26/jobs_done/`, exit 0): w_selfcopy_retrieval (09:34),
p_train_conly (09:51), w_selfcopy_envelope (09:57), w_tiebreak_probe (09:58), p_train_pca32
(10:33), p_train_wide (14:32), p_train_pca32f (15:13), p_train_pca128 (16:08),
p_featurise_esm8m and p_train_esm8m (16:50). The "_p2" records at 19:11 to 19:12 are the chain
driver re-touching finished rungs (5 s, no work). Killed without a record: a1_build_s0/s1
(per-target checkpoints under `s26/results/a1/`), ph_reject_cloud, w_train_chain,
p_eval_shipped, p_train_raw (fold checkpoints under `s26/models/p_ladder/`), p_eval_noesm.
Each resumes from its checkpoint; none restarts.

---

## L42 -- RESUMING AT 2026-09-13 19:14: A SECOND REPORT WAS FOUND IN docs/ (NOT WRITTEN BY ANY S26 LANE); LANES RELAUNCHED UNDER A 1.4 GB HEADROOM; PR SPAWNED (2026-09-13, coordinator)

1. `docs/REPORT_S26.md` (181 KB, Parts I to X, Appendices A to D, mtime 14:25) and
   `docs/REPORT_S26_SUMMARY.md` (14:24) appeared on disk during the pause, while every S26 lane
   was dead (session limit, 10:10 to 19:10). They were not written by any lane of this campaign
   and are not in git; their header names HEAD `b7f1f280`, so they were written from this
   branch by another session or by the user. They cite external sources (arXiv abstracts, an
   author page) that no S26 lane fetched. They are left untouched and uncommitted by the
   coordinator; the user decides whether they are tracked. The campaign's report deliverable
   remains `s26/REPORT.md`; lane E cross-checks the two, takes any artefact-sourced number or
   presenter material it lacks, cites the docs file when it does, and never copies an external
   claim or a personal detail into a deliverable.
2. Every lane is relaunched from its own transcript with the instruction to resume from its
   STATUS line and the checkpoints on disk. Headroom rule until the user's own load drops: one
   governed job at a time whose measured or probed peak is under 1.0 GB; the `raw` rung
   (sibling peak 1.85 GB) and lane W's Part B retrains (1.3 GB) wait; jobrun v2.2 blocks any
   launch without est-ram + 0.5 GB free.
3. Lane PR (Presentation) is spawned now: `s26/C3_RESULT.md` (L39) fixes Proposal C's AMBER
   sentence, A2/A4 (L27, L35) fix the slide 6/8 figures, and the slides 1 to 7 and 11 depend
   only on the Phase 0 claim ledger. Seven lanes active (E, Q, P, PH, W, A, PR).

---

## L43 -- STERIC REJECT, POINT CLOUD: THE FALSIFIER FIRED THE OTHER WAY. AT 1e4 THE REFILL ARM IS +0.228 A WORSE THAN THE SHIPPED TOP-75 AND +0.167 WORSE THAN REJECTING THE SAME COUNT AT RANDOM (5/5 FOLDS); THE DOSE IS MONOTONE AND ITS LIMIT IS THE ANCHOR (2026-09-13, PH)

`s26/ph_reject.py cloud` then `report`, `s26/results/ph_reject_cloud.json` (complete 126/126),
`s26/results/ph_reject_report.json`, log `s26/results/ph_reject_report_cloud.log`; job
`s26/jobs_done/ph_reject_cloud2.json` (exit 0, 170 s for the last 10 cells, peak RSS 0.041 GB;
the first 116 cells came from `ph_reject_cloud`, killed with the governor by the host at about
10:10, L40; the per-target cells made the restart lossless). Pre-registered in
`s26/PREREG_amber_reject.md` sections 3 to 5 and addendum 1; the moved-subset secondary and the
R-first reading were fixed by the coordinator at 09:10 (findings section 2.3) before any RMSD
was read.

BASIS: POINT CLOUD on both sides (`I.coordinate_average` of the retained windows, scored ORACLE
against `nat_ca`; the anchor reproduces the production `rmsd_avg` 3.0483 to 1e-6 on every
target). The built-chain twin (`rmsd_arm` basis) is job `ph_reject_chain`, running. Arms: R =
reject e_amber > T from the shipped top-75 and refill from the next-ranked survivors to m = 75
(the PRIMARY reading; it empties only when the whole pool has no survivor: 8 targets at 1e4,
which fall back to the anchor and count as ties); S = reject without refill (empties 11 sets at
1e4, same fallback). Controls matched in the operator's space: RANDR / RANDS reject the SAME
COUNT at random (16 draws, mean); PERMR / PERMS apply the same threshold to a permuted energy
vector (16 permutations). Twelve `ST.fmt` blocks verbatim, primary threshold 1e4:

      point_cloud [all, n=126] R@1e4 minus anchor (negative = arm better)
        a 3.2760 (med 3.1167)   b 3.0483 (med 2.8373)   n=126
        effect +0.2276   median +0.0027   SE 0.0693   MDE 0.1941   effect/MDE +1.17
        iid  CI95 [+0.0950, +0.3672]
        fold CI95 [+0.1466, +0.3232]   folds same sign 5/5   per-fold 0:+0.278 1:+0.076 2:+0.195 3:+0.398 4:+0.199
        49W/67L/10T   worst degradation +4.1545 (8T61)   p90 +1.2819   power 0.91  Type-M 1.06
        concentration: drop-top10 +0.3256 vs uniform-effect null p10/p50/p90 +0.2365/+0.3216/+0.4116 -> pctile 0.521
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      point_cloud [all, n=126] R@1e4 minus RANDR@1e4 (matched control)
        a 3.2760 (med 3.1167)   b 3.1088 (med 3.0266)   n=126
        effect +0.1672   median +0.0128   SE 0.0473   MDE 0.1326   effect/MDE +1.26
        iid  CI95 [+0.0799, +0.2639]
        fold CI95 [+0.0699, +0.2793]   folds same sign 5/5   per-fold 0:+0.156 1:+0.007 2:+0.191 3:+0.384 4:+0.113
        47W/77L/2T   worst degradation +3.2867 (8T61)   p90 +0.7875   power 0.94  Type-M 1.03
        concentration: drop-top10 +0.2306 vs uniform-effect null p10/p50/p90 +0.1709/+0.2276/+0.2895 -> pctile 0.523
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.03x]
      point_cloud [all, n=126] R@1e4 minus PERMR@1e4 (matched control)
        a 3.2760 (med 3.1167)   b 3.1707 (med 3.0234)   n=126
        effect +0.1053   median +0.0031   SE 0.0385   MDE 0.1078   effect/MDE +0.98
        iid  CI95 [+0.0350, +0.1821]
        fold CI95 [+0.0187, +0.2000]   folds same sign 4/5   per-fold 0:+0.070 1:-0.063 2:+0.128 3:+0.274 4:+0.115
        52W/66L/8T   worst degradation +2.6065 (9KAR)   p90 +0.6742   power 0.78  Type-M 1.14
        concentration: drop-top10 +0.1610 vs uniform-effect null p10/p50/p90 +0.1131/+0.1602/+0.2120 -> pctile 0.509
        VERDICT: NOT MEASURED (|effect| 0.1053 <= its own MDE 0.1078, 0.98x)
      point_cloud [all, n=126] S@1e4 minus anchor (negative = arm better)
        a 3.1563 (med 3.0432)   b 3.0483 (med 2.8373)   n=126
        effect +0.1079   median +0.0147   SE 0.0287   MDE 0.0805   effect/MDE +1.34
        iid  CI95 [+0.0542, +0.1679]
        fold CI95 [+0.0647, +0.1528]   folds same sign 5/5   per-fold 0:+0.043 1:+0.054 2:+0.133 3:+0.193 4:+0.117
        38W/75L/13T   worst degradation +1.7827 (8T61)   p90 +0.4359   power 0.96  Type-M 1.02
        concentration: drop-top10 +0.1427 vs uniform-effect null p10/p50/p90 +0.1043/+0.1412/+0.1794 -> pctile 0.516
        VERDICT: WORSE
      point_cloud [all, n=126] S@1e4 minus RANDS@1e4 (matched control)
        a 3.1563 (med 3.0432)   b 3.0854 (med 2.8491)   n=126
        effect +0.0709   median +0.0096   SE 0.0254   MDE 0.0711   effect/MDE +1.00
        iid  CI95 [+0.0235, +0.1229]
        fold CI95 [+0.0317, +0.1087]   folds same sign 4/5   per-fold 0:-0.000 1:+0.040 2:+0.105 3:+0.133 4:+0.077
        41W/72L/13T   worst degradation +1.5172 (8T61)   p90 +0.3281   power 0.80  Type-M 1.13
        concentration: drop-top10 +0.1041 vs uniform-effect null p10/p50/p90 +0.0714/+0.1030/+0.1371 -> pctile 0.520
        VERDICT: NOT MEASURED (|effect| 0.0709 <= its own MDE 0.0711, 1.00x)
      point_cloud [all, n=126] S@1e4 minus PERMS@1e4 (matched control)
        a 3.1563 (med 3.0432)   b 3.0809 (med 2.8611)   n=126
        effect +0.0754   median +0.0054   SE 0.0274   MDE 0.0767   effect/MDE +0.98
        iid  CI95 [+0.0253, +0.1312]
        fold CI95 [+0.0311, +0.1255]   folds same sign 5/5   per-fold 0:+0.020 1:+0.023 2:+0.113 3:+0.164 4:+0.063
        49W/69L/8T   worst degradation +1.5869 (8T61)   p90 +0.3941   power 0.79  Type-M 1.13
        concentration: drop-top10 +0.1137 vs uniform-effect null p10/p50/p90 +0.0787/+0.1119/+0.1493 -> pctile 0.526
        VERDICT: NOT MEASURED (|effect| 0.0754 <= its own MDE 0.0767, 0.98x)
      point_cloud [all, n=126] RANDR@1e4 minus anchor (negative = arm better)
        a 3.1088 (med 3.0266)   b 3.0483 (med 2.8373)   n=126
        effect +0.0605   median +0.0000   SE 0.0377   MDE 0.1056   effect/MDE +0.57
        iid  CI95 [-0.0120, +0.1392]
        fold CI95 [+0.0206, +0.0962]   folds same sign 5/5   per-fold 0:+0.121 1:+0.069 2:+0.004 3:+0.014 4:+0.086
        62W/62L/2T   worst degradation +1.6061 (8FLP)   p90 +0.5672   power 0.36  Type-M 1.65
        concentration: drop-top10 +0.1250 vs uniform-effect null p10/p50/p90 +0.0781/+0.1226/+0.1703 -> pctile 0.527
        VERDICT: NOT MEASURED (|effect| 0.0605 <= its own MDE 0.1056, 0.57x)
      point_cloud [all, n=126] RANDS@1e4 minus anchor (negative = arm better)
        a 3.0854 (med 2.8491)   b 3.0483 (med 2.8373)   n=126
        effect +0.0371   median +0.0023   SE 0.0090   MDE 0.0253   effect/MDE +1.47
        iid  CI95 [+0.0214, +0.0556]
        fold CI95 [+0.0229, +0.0495]   folds same sign 5/5   per-fold 0:+0.043 1:+0.014 2:+0.028 3:+0.060 4:+0.040
        35W/78L/13T   worst degradation +0.6973 (9KAR)   p90 +0.1193   power 0.98  Type-M 1.01
        concentration: drop-top10 +0.0450 vs uniform-effect null p10/p50/p90 +0.0330/+0.0445/+0.0570 -> pctile 0.519
        VERDICT: WORSE
      point_cloud [moved, n=116] R@1e4 minus anchor (negative = arm better)
        a 3.2910 (med 3.1277)   b 3.0438 (med 2.8573)   n=116
        effect +0.2473   median +0.0142   SE 0.0750   MDE 0.2101   effect/MDE +1.18
        iid  CI95 [+0.1024, +0.3992]
        fold CI95 [+0.1515, +0.3538]   folds same sign 5/5   per-fold 0:+0.278 1:+0.080 2:+0.244 3:+0.457 4:+0.206
        49W/67L/0T   worst degradation +4.1545 (8T61)   p90 +1.3987   power 0.91  Type-M 1.05
        concentration: drop-top10 +0.3564 vs uniform-effect null p10/p50/p90 +0.2617/+0.3505/+0.4495 -> pctile 0.526
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.05x]
      point_cloud [moved, n=116] R@1e4 minus RANDR@1e4 (matched control)
        a 3.2910 (med 3.1277)   b 3.1149 (med 3.0689)   n=116
        effect +0.1761   median +0.0170   SE 0.0507   MDE 0.1421   effect/MDE +1.24
        iid  CI95 [+0.0825, +0.2822]
        fold CI95 [+0.0809, +0.3007]   folds same sign 5/5   per-fold 0:+0.156 1:+0.017 2:+0.199 3:+0.439 4:+0.117
        43W/73L/0T   worst degradation +3.2867 (8T61)   p90 +0.8178   power 0.93  Type-M 1.04
        concentration: drop-top10 +0.2451 vs uniform-effect null p10/p50/p90 +0.1808/+0.2413/+0.3129 -> pctile 0.531
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.04x]
      point_cloud [moved, n=113] S@1e4 minus anchor (negative = arm better)
        a 3.1525 (med 3.0824)   b 3.0321 (med 2.8563)   n=113
        effect +0.1203   median +0.0233   SE 0.0319   MDE 0.0892   effect/MDE +1.35
        iid  CI95 [+0.0601, +0.1875]
        fold CI95 [+0.0668, +0.1862]   folds same sign 5/5   per-fold 0:+0.043 1:+0.057 2:+0.167 3:+0.246 4:+0.126
        38W/75L/0T   worst degradation +1.7827 (8T61)   p90 +0.4649   power 0.97  Type-M 1.02
        concentration: drop-top10 +0.1607 vs uniform-effect null p10/p50/p90 +0.1199/+0.1586/+0.2023 -> pctile 0.524
        VERDICT: WORSE
      point_cloud [moved, n=113] S@1e4 minus RANDS@1e4 (matched control)
        a 3.1525 (med 3.0824)   b 3.0734 (med 2.8570)   n=113
        effect +0.0790   median +0.0154   SE 0.0282   MDE 0.0790   effect/MDE +1.00
        iid  CI95 [+0.0257, +0.1358]
        fold CI95 [+0.0302, +0.1325]   folds same sign 4/5   per-fold 0:-0.000 1:+0.042 2:+0.132 3:+0.170 4:+0.083
        41W/72L/0T   worst degradation +1.5172 (8T61)   p90 +0.3480   power 0.80  Type-M 1.12
        concentration: drop-top10 +0.1172 vs uniform-effect null p10/p50/p90 +0.0815/+0.1153/+0.1524 -> pctile 0.531
        VERDICT: NOT MEASURED (|effect| 0.0790 <= its own MDE 0.0790, 1.00x)

All four thresholds, `[all, n=126]`, arm minus anchor and arm minus its matched random control:

    1e3  R  vs anchor  +0.5317 (+1.75x MDE 0.3045)  fold [+0.4328, +0.6566]   36W/69L/21T  WORSE         | vs RANDR  +0.4611 (+1.87x)  fold [+0.3723, +0.6036]   WORSE
    1e3  S  vs anchor  +0.1960 (+1.66x MDE 0.1182)  fold [+0.1347, +0.2330]   29W/71L/26T  WORSE         | vs RANDS  +0.1118 (+1.20x)  fold [+0.0583, +0.1424]   WORSE
    1e4  R  vs anchor  +0.2276 (+1.17x MDE 0.1941)  fold [+0.1466, +0.3232]   49W/67L/10T  WORSE         | vs RANDR  +0.1672 (+1.26x)  fold [+0.0699, +0.2793]   WORSE
    1e4  S  vs anchor  +0.1079 (+1.34x MDE 0.0805)  fold [+0.0647, +0.1528]   38W/75L/13T  WORSE         | vs RANDS  +0.0709 (+1.00x)  fold [+0.0317, +0.1087]   NOT MEASURED
    1e5  R  vs anchor  +0.0929 (+0.66x MDE 0.1410)  fold [+0.0332, +0.1433]   58W/61L/7T   NOT MEASURED  | vs RANDR  +0.0416 (+0.44x)  fold [-0.0027, +0.0963]   NOT MEASURED
    1e5  S  vs anchor  +0.0739 (+1.17x MDE 0.0634)  fold [+0.0385, +0.1077]   50W/67L/9T   WORSE         | vs RANDS  +0.0568 (+0.95x)  fold [+0.0264, +0.0871]   NOT MEASURED
    1e6  R  vs anchor  +0.0662 (+0.51x MDE 0.1305)  fold [+0.0196, +0.1013]   56W/62L/8T   NOT MEASURED  | vs RANDR  +0.0297 (+0.30x)  fold [-0.0071, +0.0588]   NOT MEASURED
    1e6  S  vs anchor  +0.0504 (+0.97x MDE 0.0517)  fold [+0.0261, +0.0717]   52W/64L/10T  NOT MEASURED  | vs RANDS  +0.0377 (+0.81x)  fold [+0.0139, +0.0590]   NOT MEASURED

Threshold sweep (an order statistic): `ST.best_of_k_within` over the four thresholds, R: oracle
per-target minimum -0.3368 A, split-half transfer -0.1470 (44%), k_eff 3.74; S: -0.1177,
transfer -0.0556 (47%), k_eff 3.72. The transfer is negative because the sweep contains 1e6,
which rejects least and so damages least: "the least harmful threshold transfers" is not a gain.
R@1e6 is +0.066 against the anchor (0.51x MDE, NOT MEASURED) and S@1e6 +0.050 (0.97x).

READING. (1) The falsifier did not clear; it fired the other way. At the primary threshold the
refill arm R is WORSE than the shipped top-75 by +0.228 A [fold +0.147, +0.323], 49W/67L/10T,
5/5 folds (1.17x MDE, Type-M zone: the size is an upper bound, the sign is measured), and WORSE
than rejecting the same number at random and refilling with no energy judgment by +0.167
[+0.070, +0.279], 5/5 folds (1.26x MDE, Type-M). The shrink arm S is worse than the anchor by
+0.108 [+0.065, +0.153], 1.34x MDE, measured; against its random twin +0.071 at exactly 1.00x
MDE, not measured. Restricting to the targets the operator actually moved (the secondary,
n = 116 / 113) changes nothing: R +0.247 and S +0.120, both WORSE; R vs RANDR +0.176, WORSE.
(2) The dose is monotone: 1e3 (54.5 of 75 rejected) +0.532 A; 1e4 (40.1) +0.228; 1e5 (30.3)
+0.093 (0.66x MDE); 1e6 (23.3) +0.066 (0.51x). The less the reject does, the less it costs, and
its limit is the anchor: S16 L27's k-ladder shape and S23 L5's mechanism (the members a physics
score removes carry error that cancels in the average) through a third operator. (3) The
matched controls say where the harm comes from. Rejecting the same COUNT at random and
refilling (RANDR) costs +0.061 at 1e4 (0.57x MDE, not measured), so refilling from ranks 76 to
147 is nearly free and the energy's CHOICE of what to reject costs the other +0.167. Rejecting
the same count at random without refill (RANDS) costs +0.037 (1.47x MDE, measured): shrinking
a 75-set to 35 at random costs 0.04 A and the energy's choice of which 40 costs 0.07 more. The
permuted-energy controls sit between random and real (PERMR +0.122, PERMS +0.033) because a
permuted single point rejects the same energies from other candidates. (4) Why: the census
(L23) found 96.8% of the members the threshold condemns are condemned by a side-chain contact
placed by the deterministic builder, and the retained set at 1e4 is +0.254 A more expanded in
Rg than the anchor; the reject keeps expanded members and discards compact ones. S25's "AMBER
selects expansion" (+1.10 A Rg on its own top-75) reproduced as a filter.

POWER. Every arm-vs-anchor contrast at 1e3 and 1e4 is above its MDE with the fold CI excluding
zero and 5/5 folds. At 1e5 and 1e6 R is 0.66x and 0.51x its MDE (SE 0.050 and 0.047, so a gain
of 0.13 A or more would have been seen): "the mildest threshold is null" is UNDERPOWERED for a
gain below 0.13 A and MEASURED against any harm above it. Nothing here is positive, so no seed
replication is due; the cross-basis replication is the built-chain run, and the direction must
hold there or the point-cloud result is a basis artefact (S16 showed the two bases can disagree
about a repair). DISPOSITION, pending the built chain: the steric reject at a physical
threshold, with refill, is closed on the point cloud in the two forms the record had not
measured (physics-set count, refill), with a sign, on all 126 and on the moved subset; the
filter form of the functional lever (s24 D1-C, s19 Q3) stays closed and gains its sixth
instrument.

---

## L44 -- LANE W: THE 2/60 BENCHMARK SELF-COPY LEAK, BOUNDED FROM THE DEV PROXY WITHOUT OPENING THE BENCHMARK: 0.002 A BY THE DEV-4 MEASUREMENT, 0.028 A (BUILT CHAIN) / 0.048 A (SELECTION) / 0.023 A (PAIRED GAIN) BY THE OWN-NATIVE ENVELOPE; PRE-REGISTERED CLASS MINOR; C27's +0.0004 DEV PRICE RE-DERIVED EXACTLY; A FOLD MODEL THAT TRAINED ON THE TARGET'S OWN NATIVE EMITS A CHAIN 0.70 A NEARER TO IT (2026-09-13, W)

Pre-registration `s26/PREREG_selfcopy_bound.md` (L30); gated job `w_endpoint_report`
(`s26/w_endpoint_report.py`: posterior, endpoint, report in one process; exit 0, 160 s, peak RSS
0.302 GB, `s26/jobs_done/w_endpoint_report.json`); artefacts `s26/results/w_selfcopy_endpoint.json`
(Parts A, B, C signed), `w_selfcopy_bound.json` (Part D), `w_selfcopy_floor.json` (Part E),
`w_selfcopy_retrieval.json` and `w_selfcopy_envelope.json` (the native-free halves, complete
126/126, exit 0, peak 0.116 / 0.318 GB, L41), `w_selfcopy_posterior.json`. No benchmark file,
sequence, native or RMSD was read; the benchmark facts used are 2/60 and the mechanism (S24 L4).
Reproduction gate: the "with" arm equals the production emission on all 126 (top-75 == `sub`
126/126, cloud max abs 1.4e-14; means sel 3.4540 / cloud 3.0483 / arm 3.2126 / fit 3.2052, the arm
being the re-projection rebuild figure of L9). Basis stated on every line; positive delta = the
leaked emission is WORSE.

**Part A, channel A (the self-window in the K = 500 pool), the four dev self-copies, PRODUCTION
basis.** Signed delta = with minus without (self-window dropped, pool refilled from the next BLOSUM
windows, S10-4's operator on exact copies):

    target   in top-75   arm      cloud    fit      sel      | triangle bound arm | control percentile (arm)
    1CEK        no      +0.0000  +0.0000  +0.0000  +0.0000  |  0.000             |  0.61
    2FBU        no      +0.0000  +0.0000  +0.0000  +0.0000  |  0.000             |  0.61
    2P5H        yes     -0.0683  -0.0114  -0.0193  +0.0000  |  0.202             |  0.89
    6B9K        yes     +0.0071  +0.0069  +0.0046  +0.0000  |  0.034             |  0.71

The self-window is never the score's argmin (sel 0.0000 on all four, and 0.0000 on all 122
controls: the BLOSUM rank-0 window is never the argmin anywhere on the instrument). The matched
control population (the rank-0 window dropped on the other 122 targets, 47 of which had it in
the top-75): |delta arm| p50 0.000, p90 0.078, p95 0.147, max 0.511; the four self-copies sit at
its 61st to 89th percentile. Registered F1 holds (|delta arm| < 0.10 and bound < 0.30 on all four).

**C27's dev half re-derived (S10-4's operator: drop every >= 0.6 window, refill; 13 targets with
such a window in the pool, the same 13 as S10-4, 21 in the universe): CLEAN minus PRODUCTION
+0.0004 on the lam = 0 chain, fold CI [-0.0001, +0.0010], iid [-0.0003, +0.0013], MDE 0.0012, 117
ties; +0.0018 on the built chain [-0.0003, +0.0039]; +0.0003 on the cloud; 0.0000 on sel.**
S10-4's +0.0004 [-0.0004, +0.0013] reproduces on the same basis (its "synthesis fit"), and its "0
on the shipped argmin" reproduces exactly. The number now has an artefact (L28 C27, L31).

**Part C, the envelope (ORACLE by construction): the four pinned fold models that had the target's
OWN native among their training labels, against the clean model, everything else identical.**

    basis   leaked-model mean minus clean   fold CI95            MDE     x MDE   W/L      folds   per model
    arm     -0.6980 (median -0.2663)        [-0.8312, -0.5764]   0.246   2.83    112/14   5/5     -0.749 / -0.646 / -0.702 / -0.694 / -0.700
    cloud   -0.7094                         [-0.8394, -0.5680]   0.244   2.91    110/16   5/5
    sel     -1.2081 (median -0.9239)        [-1.4377, -0.9745]   0.315   3.83    112/12   5/5     -1.324 / -1.102 / -1.254 / -1.163 / -1.201
    fit     -0.6932                         [-0.8234, -0.5615]   0.245   2.83    110/16   5/5
    gain (arm - sel)  +0.5101               [+0.3698, +0.6838]   0.255   2.00    38/88    5/5

The shipped MLP memorises at the pipeline endpoint: a fold model that saw the target's native emits
a built chain 0.70 A nearer to it, changes half the top-75 (overlap 0.47) and the argmin on 97% of
targets, and its argmin selection is 1.21 A nearer. The registered expectation (-0.2 to -0.8 A on
arm) held; concentration at the uniform-effect null's 52nd percentile (not concentrated); the
spread among the four leaked models (sd 0.07 A on arm) is small against the effect, so the effect
is the inclusion of the native, not the corpus fifth. The leak inflates the argmin arm more than
the built chain, so under a leak of this strength the paired gain of the architecture over the
shipped argmin is biased AGAINST the architecture by +0.51 A per leaked target. This is also the
first measurement behind EXAMINATION E4: per-target dev results may never be read across folds.

**Part B (channel B, the carrier chain in the training labels), partial: one of ten retrains
finished before the host killed `w_train_chain` (L40, L41).** 1CEK with 1A11 removed from fold 2's
corpus (`s26/models/w_selfcopy/pca32_fold2_s0_out_1A11.pt`, 276 pairs fewer), reference lane P's
`pca32_fold2_s0.pt` (which reproduces the pinned emission at 0.000): the carrier's presence is worth
-0.0114 on the built chain, -0.0204 cloud, +0.0062 sel (posterior mean |dE[d]| 1.02 A per pair,
top-75 overlap 0.76). Both channels removed on 1CEK: production is 0.0114 A BETTER on arm, 0.0177 on
the gain. 1CEK is the one dev case whose copy is near-native (Part E, 0.595 A); the three others
sit 2.3 to 4.1 A away. The control-out models and the other three carrier-out models wait for
headroom (L42) and for the tournament.

**Part E (ORACLE DIAGNOSTIC): the same sequence in a different deposit, 18 dev targets, 22 verbatim
partners: median 2.908 A from the native (min 0.317, max 5.502; 18% under 1.0 A, 27% under 1.5 A;
the natives' own ensemble spread is 1.044 A; the copy is worse than the target's pool MEAN on 5 of
22 and beats the pool's best on 3 of 22).** F5 holds. A verbatim copy is not a near-native answer at
this length; this is why channels A and B are small where they are measured.

**Part D, the bound, (2 / 60) x the per-target quantity, assumptions A1 to A5 of the PREREG (same
mechanism and direction; the benchmark 2 no worse than the dev 4 or than the envelope; same pipeline;
no length correction; the envelope bounds help, harm is measured on the dev 4):**

    source                                            arm       sel       paired gain   class (0.017 = one tenth of the benchmark CI half-width)
    dev-4 channel A, signed max                        0.0023    0.0000    0.0023         IMMATERIAL
    dev-4 channel A, native-free triangle max          0.0067    0.0000    0.0067         IMMATERIAL
    both channels removed (1CEK only, n = 1)           0.0004    0.0002    0.0006         IMMATERIAL
    own-native envelope, fold-CI limit (n = 126)       0.0277    0.0479    0.0228         MINOR

**Pre-registered verdict: MINOR, bounded at 0.028 A on the built chain, 0.048 A on the selection
basis and 0.023 A on the paired gain, because the envelope clause of Part D fires; every direct
measurement on the dev proxy is IMMATERIAL (0.002 A or less).** The envelope is loose by
construction: it prices a model trained on the target's OWN native, and the only measured carrier
effect (1CEK, the near-native case) is 60x smaller. Against the benchmark's own CI half-width 0.170
(S9-10: +0.0103 [-0.1596, +0.1803]) the bound cannot change the benchmark verdict (no validated
gain) in either direction: the un-leaked paired gain lies in [+0.0103 - 0.023, +0.0103 + 0.023]
under the envelope and within 0.003 of +0.0103 under the dev-proxy measurement. The caveat that
attaches to every benchmark figure now reads: 2/60 self-copies, bounded at 0.028 A (built chain)
by the own-native envelope and 0.002 A by the dev proxy, `s26/results/w_selfcopy_bound.json`.

Also measured on the way, native-free, and recorded for `s26/IDEA_tiebreak_noise_floor.md`: a
one-member change of the 75 flips the medoid frame (cloud moves 1.29 A on 5H1H, 0.74 on 6EY3; the
same set in the original frame is 0.06 to 0.08 A away) and flips the projection branch (chain moves
0.85 to 2.11 A on 1NIZ, 1CS9, 1M02, 1RSW, 2LNG, 2BP4 at cloud moves of 0.04 to 0.09 A); the signed
deltas on those targets are 3 to 10x smaller than the triangle bounds, so the flips move the chain
mostly orthogonally to the native. ORACLE insertion of the withheld same-fold carriers (8 targets,
11 windows): 5 of 11 enter the top-75 (F2's first clause, "at least half", misses by one); |delta
arm| < 0.10 on 7 of 8 targets; the largest move is 8ZG2, -0.277 A from a window that is itself
5.06 A from the native (a branch flip, not the window's geometry).

Deviations from the PREREG, stated: Part B has 1 of 10 models (the rest killed by the host, L40);
the `both_removed4` key of `w_selfcopy_bound.json` therefore holds n = 1; no replication of Part
C at a second seed is needed (deterministic, pinned models); Part A's control population is the
replication of its operator. Question for the coordinator: the S10-4 number is now sourced on the
`fit` basis; does the report quote the built-chain +0.0018 [-0.0003, +0.0039] beside it or in its
place?

## L45 -- ADVERSARY CHECK OF L27 (A2, THE DLA): STANDS (2026-09-13, A)

Independent re-derivation `s26/a_dla_check.py` -> `s26/results/a_dla_check.json` (88 s, run
directly, under 200 MB; written from the construction in `s26/PREREG_A2.md` and
`core/quantum.py`, sharing no closure code with `s26/q_dla.py`, and using a different numeric
rank rule -- incremental Gram-Schmidt, no shared svd tolerance). Symbolic closure at n = 4..7,
L = 1..4 equals `s26/results/q_dla.json :: results/fixed/*` at every cell: n7 dim(DLA) = 8128 =
so(128) from L = 2 (7 at L = 1); n6 510 / 1023 / 2016; n4 120 from L = 2; n5 496 from L = 2. All
four combinations of conjugation convention (right, left) and gate order (forward, reverse) give
8128 at n = 7, L = 2 and L = 3. The independent numeric route agrees with the symbolic count at
n = 4, 5, L = 1..3 (6 of 6); every closure is odd-Y. The 8128 / 1025 reconciliation in
`s26/agentQ_FINDINGS.md` 2.1 is correct: 8128 is `q_dla.json :: results/fixed/n7_L3/dim` (the
fixed ansatz, so(128), 1.0 of so); 1025 is `results/adapt_sets/adapt_L2_*_n7_a0.25_T0.3/
ladder[P=21]/dim` (the strings ADAPT selected), 1025 / 8128 = 0.126. Leakage / ties / iid-vs-fold
CI: not applicable (an exact count, no native, no RMSD). Q's own H2b prediction (dim < 8128,
guess 4095) is correctly recorded as FALSIFIED. Verdict: STANDS. The "S25 slopes are a statement
about depth 3" reading is an inference from Larocca 2022 / Ragone 2024, not a measurement here;
it is a scope note for RETRACTIONS, not a retracted number.

---

## L46 -- ADVERSARY CHECK OF L39 (C3 STAGE 1): STANDS WITH CAVEAT (2026-09-13, A)

The three contrasts reproduce from `s26/results/ph_c3_stage1.json` (complete, 126 rows,
provenance e480fc15): AMBER minus do-nothing +0.02069 (effect/MDE 2.16, fold CI [0.0154,
0.0290], 5/5, WORSE), AMBER minus matched RANDOM +0.01108 (effect/MDE 1.12, fold CI [0.0062,
0.0171], 5/5), AMBER minus TOWARD-MEMBER +0.03851 (effect/MDE 2.27, fold CI [0.0298, 0.0479],
5/5); cos(AMBER, true residual) mean -0.0491 (SE 0.0149, 36.5% positive). Checks:

1. Control construction against S16. `s26/ph_lib.py:random_displacement` reproduces
   `s16/repair.py`'s `rand` exactly: an isotropic Gaussian on the CA trace, the six rigid-body
   directions removed by projection onto `s15/align_lib.py:rigid_basis`, scaled so
   ||g||/sqrt(n) equals AMBER's per-atom RMS displacement. The magnitude matched is `mag_sup`
   0.2198 A (SE 0.0084), measured after superposition, the same definition S16 used. Confirmed.
2. The Type-M reading. AMBER-minus-RANDOM is at 1.12x MDE with `type_m_flag` True (0.7-1.3x is
   the Type-M zone). So the HEADLINE +0.0111 is a Type-M number: its sign and fold CI are clean
   (5/5, CI excludes zero) but its magnitude is inflated ~1.07x and by the discipline it is
   "not a result" as a magnitude. The two clean contrasts are AMBER-minus-do-nothing (+0.0207,
   2.16x) and AMBER-minus-member (+0.0385, 2.27x), both clear of the Type-M zone. Any
   presentation line must quote +0.0111 with the Type-M flag or lean on the two clean contrasts.
3. The decision rule holds regardless: "accuracy step" required AMBER to BEAT both controls;
   it is worse than both, so it is not an accuracy step. Note that TOWARD-MEMBER minus
   do-nothing is -0.0178 (BETTER, fold CI [-0.0255, -0.0124], 5/5): a zero-information move
   toward a random pool member IMPROVES RMSD where AMBER's physics move worsens it.
4. "Validity step" wording. Sound for the 124 targets that converge with a sane virtual bond,
   but on 2BP4 (relaxed CA-CA 5.38 A) and 9KAR (4.86 A, e1 1262 > CONVERGE_MAX_KCAL, not
   converged) the relaxation BREAKS a virtual bond; on those two it is not cleanly a validity
   step either. The presentation should say "a validity step on 124 of 126; on 2 it breaks a
   virtual bond, one of which does not converge."

Verdict: STANDS WITH CAVEAT (the +0.0111 headline is Type-M; "validity step" carries the 2BP4 /
9KAR exception). The finding's direction and its central conclusion are correct.

---

## L47 -- ADVERSARY CHECK OF L35 (A4, GROWN-CIRCUIT VARIANCE): STANDS WITH CAVEAT (2026-09-13, A)

The reproduction gate holds: `s26/results/q_var.json :: results/reproduction/passed` True,
`worst_rel` 0.0 (the S25 n = 7 rows of all five cells reproduce bit-identically). The slopes
reproduce (`results/slopes`): fixed deployed_a1_T03 -0.311, grown V +0.024, grown L2 -0.008;
the one non-trivial family, deployed_a025_T03, fixed -0.243, grown V -0.239, grown L2 -0.302.
Checks:

1. "Product circuit at alpha = 1" is established, but by `s26/results/q_dla.json`, not by
   `q_var.json`. `q_var.json` carries only the variance slopes; the product-circuit reading
   rests on the ADAPT closures in `q_dla.json :: results/adapt_sets` (alpha = 1 selects an
   abelian, dim-7 = n set, `n_distinct_ops` 1 to 2, the ladder dim stays 7), which I
   independently reproduced in `s26/results/a_dla_check.json` (L27 above). A large gradient
   from a product circuit is the trivial regime, not trainability (rule 10's mirror). Sound,
   with that provenance noted.
2. The one non-trivial slope's comparison. L35 says grown L2 -0.302 is "the fixed ansatz's
   -0.243 within the sampling error of a 7-point slope (relative SE of a variance 9 to 16%)."
   `q_var.json` stores the point slope but no CI on it, so "within error" is asserted, not
   computed: the difference is 0.059 in log2-slope-per-qubit and no SE on that difference is on
   disk. The qualitative conclusion is safe (both slopes are order -0.25 to -0.30, both far from
   the -1.0 of a 2-design), but the precise "within error" claim is not backed by a stored CI.
   Also: the ledger table's grown-V value -0.246 for that cell is the P = 3n-only 6-point fit,
   while `q_var.json :: slopes/deployed_a025_T03/grown_V_adam_best` is -0.239; the two differ
   because of the degenerate-row exclusion the ledger states, not a discrepancy.

Verdict: STANDS WITH CAVEAT (the negative conclusion -- ADAPT gives Proposal A no width-scaling
argument -- holds; the "within error" of -0.302 vs -0.243 is asserted without a persisted slope
CI, and "product circuit" is grounded in q_dla.json, not q_var.json).

---

## L48 -- ADVERSARY CHECK OF L22, L23, L24, L38 (PH CENSUSES AND THE CIS FLOOR): ALL STAND; L22 WITH A STALE-LINE-NUMBER CAVEAT (2026-09-13, A)

L22 (cis census, `s26/results/ph_cis_census.json`, complete). Omega statistics recomputed from
its rows: 1,507 bonds (= sum of n-1), mean 1.91, median 0.34, p90 5.8, p99 17.4, max 42.8 deg
at 9UV5 (bonds -137.2 and +144.8), 0.53% beyond 20 deg, 0.13% beyond 30; 0 of 126 model-1
natives cis by either criterion (the two agree 126/126), 0 of 1,966 ensemble models, 0 of
2,352,893 windows, universe minimum step 3.5045 A. Native-free path confirmed (natives through
`ph_lib.native_backbone`, omega/CA only; bank through `univ_nativefree`, which refuses `rr` /
`nat_ca`; record through `prod_record_nativefree`, RMSD keys stripped). CAVEAT: L22 cites the
step gate at `core/data.py:406-407` and `697-698`; at HEAD it is `core/data.py:439`
(`step.min() < 3.5 or step.max() > 4.1`) and `:725` -- the numbering before lane I's `37bddbbb`.
Same code, stale line numbers. Q's registered ensemble-cis prediction is correctly recorded as
FALSIFIED (0 of 1,966). Verdict: STANDS WITH CAVEAT.

L23 (steric reject census, `s26/results/ph_reject_census.json`, complete). The pool-identity
assertion is live code (`s26/ph_reject.py:96-102`: `universe_idx == I.pool_idx`,
`amber_verify_max_rel == 0.0`, production `sub` == score top-75 as a set) and holds on 1A13,
1S9Z, 9KAR (`s24/cache_amber/<pdb>.npz :: universe_idx` equals `order[:500]`; top-75 members
above 1e4 = 12 / 27 / 74, equal to the census rows; pool frac above 1e4 = 0.390 / 0.650 /
0.888). 96.8% = (2,627 + 2,267) / 5,057 side-chain-involving; rho(e, min heavy-atom distance)
within top-75 -0.743 (SE 0.017). Verdict: STANDS.

L24 (C3 native-free part). Every number recomputed from `bench_results/cache/1fc9f2dcf489e2fb`
with my own Kabsch: displacement 0.2198 A (SE 0.0084, median 0.1972, range 0.103-0.591),
`amber_moved` 0.2332, e0 median 8.59e4 / min -473 / max 1.3e14 / 58.73% above 1e4, e1 mean
-559.8 (SE 30.0) / max 1262.4 (9KAR), 125 of 126 converged, strain 60.4, relaxed bond mean
3.867 (min 3.12 at 1M02, max 5.38 at 2BP4, 4.86 at 9KAR), 6 of 126 outside [3.6, 4.0], Rg
+0.0457. All match L24. Verdict: STANDS.

L38 (cis floor, ORACLE DIAGNOSTIC, `s26/results/ph_cis_floor.json`, complete). floor_ca mean
0.3468 A (SE 0.029, median 0.272, max 1.474 at 1ID6); rebuild_bb 0.3400; chain_cost 0.1664;
rho(floor_ca, max omega dev) +0.828, rho(chain_cost, omega dev) -0.036, rho(chain_cost,
floor_ca) +0.083; the paired contrast cost-minus-floor -0.1804 (fold CI [-0.2266, -0.1022],
verdict BETTER). Every quantity reads the native and is ORACLE-labelled; the "BETTER" is
explicitly disarmed in L38 (it means only that the production chain cost is smaller than the
own-torsion rebuild floor). floor_ca is declared an UPPER bound on the manifold floor and the
tight `floor2` is registered. The two Spearman nulls (+0.08, -0.04) exclude |rho| above ~0.25
at n = 126. Verdict: STANDS.

---

## L49 -- ADVERSARY CHECK OF L30 (W's 2/60 PROXY-BOUND PREREG): SOUND, STANDS (2026-09-13, A)

L30 is a pre-registration plus a native-free census, not an endpoint result; the check is on
its soundness and on whether the census reads anything it must not.

Reads. `s26/w_selfcopy.py census` iterates `I.targets()` (the 126 dev targets), loads every
universe through `load_blind` (which overwrites `rr` and `nat_ca` with NaN, lines 105-115), and
reads the production record through `I.shipped_record` (native-free). No benchmark sequence,
name, PDB, native, RMSD or manifest is read anywhere in the census, retrieval, envelope or
posterior commands; the only native reads (`I.load_univ` at lines 686, 821) are inside the
gated `endpoint` and `floor`, which refuse to run before "PHASE 0 SIGNED OFF". The one
benchmark-derived fact used, 2/60, is on the record (S24 L4, lane I L15/L18). Confirmed clean.

Assumption set (A1-A5). Reasonably complete for the quantity claimed (the leak's contribution
to a benchmark mean, bounded in absolute value). A1 (mechanism match, no interaction when one
carrier carries both benchmark targets) and A5 (the envelope bounds channel B in the HELP
direction only; HARM is measured on the dev 4 directly) are the two load-bearing assumptions and
both are stated. The gap a reader would press -- that the benchmark carriers might be MORE
homologous to their targets than the dev-4 carriers, which would make the benchmark effect
larger than "max of four" -- is covered by the ORACLE-insertion arm (8 dev targets with a
same-fold verbatim carrier, longer identity 0.6 to 0.93 against the dev 4's <= 0.52), reported
beside the dev 4. The honest limitations are self-declared: n = 4 supports no quantile ("max of
four" is named as such), and the AMBER stage's second-order contribution to a per-target delta
is an assumption (A3), not a measurement.

Verdict: STANDS as a sound pre-registration. The bound's headline will rest on n = 4 for the
realistic arm and on the n = 126 envelope for the guard; when the gated endpoints land, the
Adversary re-checks the signed deltas, the triangle bounds and the envelope's fold CI against
this prereg before any benchmark caveat text is written.

---

## L50 -- RULING ON L44's QUESTION: THE REPORT QUOTES BOTH LEAK PRICES, EACH WITH ITS BASIS; NEITHER REPLACES THE OTHER (2026-09-13 19:27, coordinator)

Lane W asked whether the report quotes the built-chain +0.0018 A [-0.0003, +0.0039] beside
S10-4's +0.0004 A [-0.0004, +0.0013] or in its place. Both, each with its basis named: the
+0.0004 is the S10-4 figure on the lam = 0 chain (`fit`), now re-derived exactly
(`s26/results/w_selfcopy_endpoint.json`, L44) and is the number the state brief carries; the
+0.0018 is the same operator on the production built chain and is the number that matches the
report's production basis. Neither clears its MDE (0.0012 on `fit`; the built-chain CI includes
zero) and both say the dev leak is immaterial; the envelope (0.028 A built chain) is the bound
the report states for the benchmark's 2/60, labelled MINOR as pre-registered. The 0.70 A
own-native effect (fold CI [-0.83, -0.58], 5/5 folds) is an ORACLE fact about leaked training,
not a price of the leak that exists, and the report says so in the same sentence.

---

## L51 -- TOURNAMENT RANKED (`s26/TOURNAMENT.md`); ASSIGNMENTS AND THE RUN ORDER UNDER TONIGHT'S HEADROOM (2026-09-13 19:28, coordinator)

The Adversary's ranking stands as posted. Assignments, in run order, one governed job at a time
beside the three endpoint jobs already running (a1_build, p_eval_*, ph_reject_chain):

1. product_state_optimum (Q): rides A1; Q reports it with A1's ledger entry.
2. conformational_identity_floor (W): gated, under a minute; W runs it first.
3. strain_difficulty (PH): reads the 126 cached records; PH runs it between chain cells.
4. tiebreak_noise_floor (W): about 50 min CPU, checkpointed; W runs it after item 2.
5. l17_target_dependent_hamiltonian = A3 (Q): already chained behind A1.
   branch_select (PH): AMBER, one job, about 1.0 GB; PH probes one target first and runs it
   only when the governor snapshot shows 1.5 GB free (jobrun v2.2 enforces it).
6. window_ensembling (P's idea, orphaned to W): W runs it after items 2 and 4; it is the
   mandatory test-time-ensembling direction and must state how it differs from widening K.
   rotamer_relief B (PH), window_provenance (W), amber_prior_partner (orphaned to W): as
   capacity allows, in that order, after everything above.
7. Deferred on memory: better_prior_inputs `attn` (3 to 3.5 GB) and coherence_penalised_training
   (1.25 GB); the `ragp` rung folds under C2 in lane P's queue. If the box frees up, W takes
   coherence_penalised_training; `attn` needs the box to itself and is unlikely tonight.

The amber_reject built-chain arm remains PH's (the falsifier's last leg, running as
ph_reject_chain). Every survivor still needs its PREREG on disk before compute (the ones
without one: strain_difficulty, window_ensembling, window_provenance, amber_prior_partner,
product_state_optimum if its arm is not already inside PREREG_A1).

---


## L52 -- TOURNAMENT ITEM 2, conformational_identity_floor (W): THE SAME SEQUENCE IN A DIFFERENT DEPOSIT SITS 2.9 A FROM THE NATIVE (MEDIAN OVER 22 PAIRS; 3.1 A OVER 18 TARGETS), 0.97 A WORSE THAN THE POOL'S BEST WINDOW AND 1.44 A BETTER THAN ITS MEAN; F5 HOLDS; ORACLE DIAGNOSTIC (2026-09-13, W)

Pre-registered as Part E of `s26/PREREG_selfcopy_bound.md` (falsifier F5: median above 1.5 A;
FALSIFIED if below 1.0 A) and, for the two paired contrasts, `s26/PREREG_identity_floor.md`
(written before `s26/w_identity_floor_stats.py` ran). Measurement: job `w_selfcopy_floor` (exit 0,
5 s; `s26/results/w_selfcopy_floor.json`, complete 22/22, provenance e480fc15); statistics: job
`w_identity_floor_stats` (exit 0, 5 s, peak 0.004 GB; `s26/results/w_identity_floor.json`). Every
quantity reads the native (the copy's CA-RMSD to the target's model-1 native at the shared
segment; the pool's ORACLE `rr`): ORACLE DIAGNOSTIC, nothing selects, nothing is deployable.
Population: the 18 dev targets with a verbatim relative in the peptide database (lane I's L15
list; 4 cross-fold carriers, the rest same-fold), 22 partners; benchmark sequences never used.

    cross-deposit CA-RMSD, 22 pairs:  median 2.908   mean 2.811   min 0.317 (5V5B)   max 5.502 (7S3O)
                                       below 1.0 A: 4/22 (18%)     below 1.5 A: 6/22 (27%)
    per target (partners averaged), 18:  median 3.055   below 1.0 A: 4/18   below 1.5 A: 5/18
    the four cross-fold self-copies:  1CEK 0.595   2FBU 3.278   2P5H 2.334   6B9K 4.126   (S24 L4's four, reproduced)
    role: carrier segment (n = 15) median 2.98;  carried whole chain (n = 7) median 2.04;  Spearman(length ratio, RMSD) -0.08 (p 0.72)
    reference scales: the natives' own intra-ensemble spread 1.044 A (record); the S24 L8 "universe best" 1.313 A

F5 holds (2.908 against the 1.5 A line; 18% below 1.0 A against "fewer than a third"). The two
registered paired contrasts (`ST.fmt` verbatim; basis: a single window against the native on both
sides; n = 18 targets, five folds):

  copy_minus_pool_best (ORACLE both sides; single window vs native), n=18 targets
    a 2.7968 (med 3.0548)   b 1.8291 (med 1.9778)   n=18
    effect +0.9677   median +0.9703   SE 0.3269   MDE 0.9158   effect/MDE +1.06
    iid  CI95 [+0.3349, +1.5834]
    fold CI95 [+0.7247, +1.2809]   folds same sign 5/5   per-fold 0:+1.113 1:+1.523 2:+0.736 3:+0.635 4:+0.937
    3W/15L/0T   worst degradation +2.9126 (7S3O)   p90 +2.5643   power 0.84  Type-M 1.10
    concentration: drop-top10 +2.1403 vs uniform-effect null p10/p50/p90 +1.6224/+2.1031/+2.5080 -> pctile 0.544
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.10x]
  copy_minus_pool_mean (ORACLE both sides; single window vs native), n=18 targets
    a 2.7968 (med 3.0548)   b 4.2396 (med 4.1842)   n=18
    effect -1.4428   median -0.9970   SE 0.4504   MDE 1.2618   effect/MDE -1.14
    iid  CI95 [-2.3156, -0.6046]
    fold CI95 [-1.7172, -1.1896]   folds same sign 5/5   per-fold 0:-1.855 1:-1.217 2:-1.773 3:-1.452 4:-1.023
    13W/5L/0T   worst degradation +0.9307 (8ZG2)   p90 +0.4068   power 0.89  Type-M 1.06
    concentration: drop-top10 +0.1234 vs uniform-effect null p10/p50/p90 -0.3330/+0.0854/+0.4138 -> pctile 0.555
    VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.06x]

Reading, with the discipline's labels: the copy is WORSE than the pool's own best window (+0.97 A,
5/5 folds, fold CI excluding zero, but 1.06x MDE: Type-M zone, sign clean, magnitude not a result)
and BETTER than a typical pool window (-1.44 A, 5/5 folds, 1.14x MDE: Type-M zone likewise). Both
were registered as expected in sign; the second was registered as "not measured", and it is
measured in the Type-M sense only. The number for the report: sequence identity buys a window that
is, in the median, 2.9 A from the native, i.e. no nearer than the shipped built chain's 3.21 A mean
by a margin the instrument can call at n = 18, and 1 A worse than the best window the pool already
holds. The 2.0 A target of S15 to S25 is below what a verbatim sequence lookup reaches on this
instrument. Two consequences carried forward: the self-copy leak is small where it is measured
because the copy is not the native (L44); and "containment-fresh" novelty (memory
`no-fresh-benchmark-exists`) is a weaker notion than the record has treated it as, since a
verbatim copy is worth about one pool window. Replication: deterministic (Kabsch on fixed
coordinates; no seed, no fit); the only random element, the bootstrap, is seeded by `stats_lib`.
Not done: a length-matched non-verbatim control population (a random same-fold window against
each target) is the natural third contrast and was not registered; it is the pool mean's role
here and is left for the report's reader as such.

## L53 -- STRAIN DIFFICULTY (tournament 4): HOW FAR THE RELAXATION MOVES THE BUILT CHAIN PREDICTS ITS ERROR (SPEARMAN +0.433 WITH rmsd_arm, FOLD CI [+0.25, +0.59], 5/5 FOLDS, PARTIAL ON n AND Rg, REPLICATED); QUARTILE MEANS 2.29 / 2.94 / 3.76 / 3.92 A; A CALIBRATION FLAG, NOT A LEVER (2026-09-13, PH)

`s26/ph_strain.py`, `s26/results/ph_strain.json` (complete 126/126) and the registered
replication `s26/results/ph_strain_rep.json` (different bootstrap seed, targets in reversed
order); jobs `s26/jobs_done/ph_strain.json` and `ph_strain_rep.json` (exit 0, 35 s each, peak
RSS 0.1 GB). Pre-registered in `s26/PREREG_strain_difficulty.md` (tournament rank 4, L51),
committed at 59d8e934 before the run. Reads the 126 production records
(`bench_results/cache/1fc9f2dcf489e2fb`) and nothing else. CALIBRATION, not an accuracy lever:
nothing is selected, tuned or moved; `rmsd_arm` (BUILT CHAIN basis) is read only as the ORACLE
label of the already-emitted structure. Signals, all native-free and emitted by the production
relaxation for free: `log_e0` (the built chain's own AMBER energy before relaxation),
`log_drop` (energy removed), `moved` (restraint RMSD on N/CA/C, how far the chain moved),
`strain_after` (bond + angle energy left). Confounds n and Rg partialled by residualising both
variables on them.

Spearman with the ORACLE `rmsd_arm`, n = 126, 4000-draw iid and fold-clustered bootstrap CIs,
4000-draw label-permutation p, per-fold signs:

    raw rho(log_e0, rmsd_arm)                     rho +0.224  iid [+0.061, +0.376]  fold [+0.116, +0.328]  perm p 0.0110  folds same sign 5/5
    raw rho(log_drop, rmsd_arm)                   rho +0.226  iid [+0.060, +0.376]  fold [+0.118, +0.328]  perm p 0.0112  folds same sign 4/5
    raw rho(moved, rmsd_arm)                      rho +0.431  iid [+0.267, +0.565]  fold [+0.278, +0.565]  perm p 0.0000  folds same sign 5/5
    raw rho(strain_after, rmsd_arm)               rho +0.175  iid [+0.003, +0.348]  fold [+0.109, +0.236]  perm p 0.0470  folds same sign 5/5
    partial rho(log_e0, rmsd_arm | n, Rg)         rho +0.241  iid [+0.070, +0.391]  fold [+0.047, +0.379]  perm p 0.0075  folds same sign 4/5
    partial rho(log_drop, rmsd_arm | n, Rg)       rho +0.250  iid [+0.081, +0.402]  fold [+0.068, +0.381]  perm p 0.0037  folds same sign 4/5
    partial rho(moved, rmsd_arm | n, Rg)          rho +0.433  iid [+0.272, +0.573]  fold [+0.247, +0.588]  perm p 0.0000  folds same sign 5/5
    partial rho(strain_after, rmsd_arm | n, Rg)   rho +0.133  iid [-0.052, +0.305]  fold [-0.178, +0.452]  perm p 0.1353  folds same sign 3/5

    FAIL18 in the top quartile of log_e0       : 7/31 (18 of 126 overall)  odds 2.23  one-sided p 0.113
    FAIL18 in the top quartile of log_drop     : 7/31 (18 of 126 overall)  odds 2.23  one-sided p 0.113
    FAIL18 in the top quartile of moved        : 6/31 (18 of 126 overall)  odds 1.66  one-sided p 0.257
    FAIL18 in the top quartile of strain_after : 4/31 (18 of 126 overall)  odds 0.86  one-sided p 0.698
    confounds: rho(n, rmsd_arm) +0.114; rho(Rg, rmsd_arm) +0.213; rho(n, log_e0) +0.340; rho(Rg, log_e0) +0.140

THE FALSIFIER CLEARS ON ONE SIGNAL, `moved`, WITH MARGIN: partial rho +0.433 (bar 0.25), fold CI
[+0.247, +0.588] excluding zero, 5/5 folds, permutation p < 0.00025 (bar 0.0125 after Bonferroni
over four signals). REPLICATED as registered (new seed, reversed order):

    partial rho(moved, rmsd_arm | n, Rg)          rho +0.433  iid [+0.272, +0.572]  fold [+0.248, +0.581]  perm p 0.0000  folds same sign 5/5
    partial rho(log_e0, rmsd_arm | n, Rg)         rho +0.241  iid [+0.080, +0.394]  fold [+0.047, +0.375]  perm p 0.0063  folds same sign 4/5
    partial rho(log_drop, rmsd_arm | n, Rg)       rho +0.250  iid [+0.084, +0.413]  fold [+0.068, +0.381]  perm p 0.0037  folds same sign 4/5

The replication lands inside the first run's fold CI on every quantity. `log_e0` and `log_drop`
are positive at +0.24 and +0.25 partial but fail the 5/5-fold sign rule (4/5) and sit at the
0.25 bar; `strain_after` is null (partial +0.13, fold CI spanning zero). The FAIL18 Fisher test is
null for every signal (best one-sided p 0.113), as predicted: FAIL18 is about retrieval recall,
not strain.

The presentable form (built chain, ORACLE labels, quartiles of `moved`, 32 targets each):

    moved (A)     < 0.159      0.159 to 0.218    0.218 to 0.288    >= 0.288
    mean rmsd_arm   2.286         2.936             3.758             3.923

A chain the force field has to move 0.29 A or more to make physical has a mean error 1.6 A larger
than one it moves less than 0.16 A. Mechanism checks (ORACLE, diagnostic): without 9KAR and 2BP4
(the two broken-bond emissions, moved 0.61 and 0.50) rho is +0.41, so it is not two outliers;
rho(moved, log_e0) is +0.41 while `moved` beats `e0` on the label by 0.19 in rho, so it is not the
raw energy in disguise; the top-8 by `moved` contain both easy (1D6X 1.78 A) and hard (9KAR
7.44, 2MFV 6.04) targets, so it is a graded signal and not a flag for a few catastrophes.

What it is and is not. It is the first native-free quantity in this programme's record with a
correlation above 0.4 to the per-target error of the emitted structure (the routers of S22 L7 /
S23 L7 and the compactness proxies of `in-band-ordering-is-per-target` reached 0.24 to 0.37 on
the per-target sign). It is a confidence flag the presentation can attach to every emitted
structure at zero cost ("how far the physics had to move this chain to make it physical").
It is NOT an accuracy lever: the prereg forbids converting it into a selector or a weight, and
the record says every such conversion fails held out (S22 L7, S23 L7). The physical reading is
plain: a coordinate average that the projection turned into a strained chain is one whose pool
members disagreed, and disagreement is error. Power: at n = 126 the design resolves |rho| of
0.25 at about 0.8; the three weaker signals are at or below that bar and are reported as
measured, not as nulls. HYPOTHESIS for lane PR: quote it as a calibration curve, never as a
gain. The Adversary's check is invited.

---

## L54 -- ADVERSARY CHECK OF L43 (STERIC REJECT, POINT CLOUD): STANDS WITH CAVEAT (2026-09-13, A)

A harmful result, so the checklist is applied to the controls and the power, not to a gain.

- Leakage: the operators (`reject_refill`, `reject_shrink`, `random_shrink`, `random_refill`,
  `permuted_energy`, `retained_sets`, `s26/ph_reject.py:109-200`) read `e_amber`, `score_dist`,
  `sub` and `order` only; no `nat_ca` and no RMSD inside them (the `rr` at `:169` is an rng handle,
  not the ORACLE array; `univ_nativefree` refuses `rr`). The native enters only in `cloud_rmsd`,
  the ORACLE scoring of a native-free operator. Clean.
- Tie-breaking: the score order is `argsort_stable` (the production rule); an empty retained set
  falls back to the anchor and is counted as a tie (10 ties at 1e4 R, 13 at S), as pre-registered
  in addendum 1. Clean.
- iid vs fold CI: both quoted verbatim; at 1e4 both exclude zero for R vs anchor and R vs RANDR,
  5/5 folds. Clean.
- Concentration / median-vs-mean: the uniform-effect null is computed and not flagged (pctile
  0.52). But the harm is TAIL-CARRIED: R@1e4 vs anchor has median +0.0027 against mean +0.2276,
  49W/67L/10T, p90 +1.28, worst +4.15 (8T61). The presentation must say "near zero on the median
  target, catastrophic on a minority (the 8 targets whose whole pool has no survivor, and the deep
  refills)", never "+0.228 on every target". CAVEAT 1.
- k_eff for the threshold sweep: `ST.best_of_k_within` applied (k_eff 3.74 R / 3.72 S), and the
  negative split-half transfer is correctly read as "the mildest threshold transfers", not a gain.
  Clean.
- Tuned parameter: the primary 1e4 was fixed before the run (PREREG section 3). Clean.
- Baseline / basis: the anchor is the production `rmsd_avg` (reproduced to 1e-6 on every
  target); controls matched in count, refill and permutation; point cloud on both sides, stated.
  Clean.
- Type-M: the two HEADLINE numbers, R@1e4 vs anchor +0.228 (1.17x MDE) and R vs RANDR +0.167
  (1.26x), are Type-M-zone magnitudes. The direction "harmful" rests on the measured contrasts:
  1e3 R +0.532 (1.75x), 1e3 S +0.196 (1.66x), 1e4 S +0.108 (1.34x), RANDS +0.037 (1.47x), 5/5
  folds throughout, and on the monotone dose. Quote the 1e4 R magnitudes with the flag. CAVEAT 2.
- "Refilling from ranks 76 to 147 is nearly free" rests on RANDR vs anchor +0.061 at 0.57x MDE
  (UNDERPOWERED, iid CI includes zero), so the refill-cost / choice-cost decomposition is a
  point-estimate split, not a measured one. CAVEAT 3.
- Replication: a negative; the cross-basis replication is the built-chain job `ph_reject_chain`,
  running; the direction must hold there. Power stated correctly for the 1e5 / 1e6 nulls.

Verdict: STANDS WITH CAVEAT (the conclusion -- the physical-threshold reject with refill is
harmful on the point cloud, its limit is the anchor -- is established; the +0.228 / +0.167
magnitudes are Type-M, the harm is tail-carried, and the refill-cost split is underpowered).
Disposition "closed on the point cloud in the two unmeasured forms" is accepted pending the
built chain.

---

## L55 -- ADVERSARY CHECK OF L44 (THE 2/60 PROXY BOUND): STANDS WITH CAVEAT; THE CLASS MINOR IS STABLE UNDER EVERY READING OF THE ENVELOPE, THE NUMBER 0.028 IS NOT (2026-09-13, A)

- Leakage / reads: the native-free halves use `load_blind` (`rr`, `nat_ca` NaN-poisoned); the
  gated `endpoint` and `floor` read natives for scoring only; no benchmark file, sequence or name
  anywhere (confirmed in L49 and re-checked on `s26/w_endpoint_report.py`'s inputs). Clean.
- Ties: `emit` uses `argsort(kind="stable")`, the argmin tie set is stored and `sel` is averaged
  over it (`ST.argmin_tied`'s rule); the refill follows the production corpus order. Clean.
- iid vs fold: Part C on arm -0.698, fold CI [-0.831, -0.576], 2.83x MDE, 112W/14L, 5/5;
  concentration at the null's 52nd percentile; the median (-0.266) sits well inside the mean
  (-0.698), a broad but skewed effect. Clean, ORACLE-labelled, and L50 rightly keeps it out of
  the price of the leak that exists.
- Bound arithmetic: every row of `s26/results/w_selfcopy_bound.json` recomputes: (2/60) x
  0.0683 = 0.0023 (2P5H); (2/60) x 0.2016 = 0.0067 (triangle); (2/60) x 0.8312 = 0.0277 (arm
  envelope fold-CI limit); (2/60) x 1.4377 = 0.0479 (sel); (2/60) x 0.6838 = 0.0228 (gain).
  Materiality thresholds pre-registered (0.017 / 0.170). Clean.

Two caveats.

1. **Artefact / ledger disagreement on the paired gain.** `w_selfcopy_bound.json :: verdict/gain`
   reads IMMATERIAL at 0.0006 (the n = 1 `both_removed4` row) because `C_envelope_fold_ci` carries
   no `gain` row and `report()`'s `max(cands)` (`w_selfcopy.py:913`) therefore never saw the
   envelope for that basis; L44's Part D table and its pre-registered verdict say MINOR at 0.023
   (from the gain row +0.5101 [+0.3698, +0.6838]). The ledger's class is the prereg-correct one
   (Part D: IMMATERIAL only if B_real AND B_env are both below 0.017). Lane W should add the
   envelope gain row to the artefact so the report cites a JSON that agrees with the ledger.
2. **B_env is the fold-CI limit of a MEAN effect, not a per-target bound.** The prereg defined it
   so (Part D, assumption A2), and A2 is stated, but a bound on the contribution of two SPECIFIC
   targets is (2/60) x a per-target quantity. From `s26/results/w_selfcopy_endpoint.json :: C/rows`
   (leaked minus clean, arm, mean over the four leaked models): median -0.266, p05 -2.490, worst
   -4.536 A (2BP4); over (target, model) pairs the worst is -4.557. So the same envelope gives
   (2/60) x 4.536 = **0.151 A** as the worst-single-target bound and (2/60) x p95 = 0.083 A;
   the artefact's own native-free triangle rows say the same (`C_envelope_p95_over_targets/arm`
   0.147, `C_envelope_max_over_models/arm` 0.211, the latter loose because branch flips move
   the chain orthogonally to the native, as L44 notes). The class MINOR therefore holds under
   every reading of the envelope (0.028 mean-CI, 0.083 p95, 0.151 worst target; all below
   0.170), which is a stronger statement than L44 makes; but "bounded at 0.028 A" must be quoted
   as "expected contribution 0.023, mean-CI limit 0.028, worst single target 0.151, under A2",
   and the envelope named as the over-bound it is (own-native training, 60x the one measured
   carrier effect).

Also noted, not a caveat: Part B is 1 of 10 retrains (host kill, L40), stated as a deviation;
F2's "at least half" missed by one (5 of 11) and is reported as such; F1, F3 (n = 1), F4, F5
hold. Verdict: STANDS WITH CAVEAT. The benchmark caveat text should read: "2/60 self-copies;
dev-proxy price 0.002 A; own-native envelope 0.028 A (mean CI) to 0.151 A (worst target), MINOR
under every reading; cannot move the benchmark verdict either way."

---


## L56 -- LANE P, C2 RUNG SHIPPED: THE ANCHOR. THE SHIPPED POSTERIOR THROUGH THE LADDER'S OWN PATH IS BIT-EXACT AGAINST THE PRODUCTION CACHE ON ALL FOUR BASES (2026-09-13 19:34, lane P)

Artefacts: `s26/results/p_ladder_shipped_s0.json` (126 rows, complete), `s26/results/p_ladder_report_shipped_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_shipped2` (resumed from a 30-row checkpoint): exit 0, 416 s, peak RSS 0.13 GB (`s26/jobs_done/p_eval_shipped2.json`); the first attempt `p_eval_shipped` was suspended by the governor at 09:25 and terminated at 19:11 (L42). Every later rung is paired against these 126 rows. The per-target production values are `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` (`shipped`, `rmsd_avg`, `rmsd_arm`, `rmsd_fit`).

```
SHIPPED THROUGH THE LADDER'S OWN PATH (the anchor), n=126
  sel   mean 3.4540 (pinned 3.4540)  vs production cache per target: max abs 0.00e+00, mean abs 0.00e+00, n(diff>1e-6)=0
  cloud mean 3.0483 (pinned 3.0483)  vs production cache per target: max abs 7.25e-14, mean abs 2.69e-15, n(diff>1e-6)=0
  arm   mean 3.2126 (pinned 3.2148)  vs production cache per target: max abs 1.71e-01, mean abs 1.07e-02, n(diff>1e-6)=120
  fit   mean 3.2052 (pinned 3.2041)  vs production cache per target: max abs 1.41e-01, mean abs 1.80e-03, n(diff>1e-6)=118
  gam_eff/cos of the shipped posterior against itself: prob +0.0000/nan loc +0.0000/nan (must be 0/nan)  MAE 2.3386
```

---

## L57 -- CORRECTION TO L56: THE ANCHOR'S BUILT CHAIN IS THE LEADERBOARD-REBUILD BASIS (3.2126), NOT THE PRODUCTION CACHE (3.2148); SELECTION AND POINT CLOUD ARE BIT-EXACT (2026-09-13 19:35, lane P)

L56's title says "bit-exact on all four bases"; the block under it says otherwise and the block is
right. Through `s26/p_ladder.py`'s path the shipped posterior reproduces the production cache
per target at 0.00e+00 on selection and 7.25e-14 on the point cloud, but the BUILT CHAIN lands at
mean 3.2126 (the S13 `base` / S24 L9 "leaderboard rebuild" figure; the state brief lists both:
"3.2148 (leaderboard rebuild: 3.2126)"), with 120/126 targets differing from `rmsd_arm` by more
than 1e-6: mean abs 0.0107, max abs 0.1711, signed mean -0.0021, worst
7JS6 0.171, 6QAX 0.148, 2LWU 0.121, 2NDN 0.115, 2LNG 0.101. Mechanism: the cloud enters `s12.instrument.project`
differing from the persisted `avg_ca` by 1e-14, and the multi-start L-BFGS projection lands in a
different local optimum on most targets; the 1A13 probe (L11) matched at 0.0 because it was fed
the persisted `avg_ca` itself. Consequence for the ladder: every rung and the anchor are built
through the SAME process and code path, so the paired contrasts are on one basis, the
rebuild basis 3.2126, and that basis is stated in every C2 entry from here on; the 0.0022 A
offset to the production figure is a basis difference, not an effect. The `fit` (lam = 0) basis
lands at 3.2052 (production 3.2041) for the same reason. L56 stands with its title read as
"bit-exact on selection and point cloud".

---

## L58 -- ADDENDUM TO L44 PER L55: THE BOUND ARTEFACT NOW CARRIES THE GAIN ENVELOPE ROW AND THE PER-TARGET READINGS; 0.028 IS QUOTED WITH 0.151 (WORST TARGET, 2BP4) BESIDE IT; ON THE SELECTION BASIS THE WORST-TARGET READING IS 0.194, ABOVE THE 0.170 LINE (2026-09-13, W)

Both L55 caveats applied by `s26/w_bound_addendum.py` (run directly, 4 s, one JSON read, under
200 MB; `s26/w_selfcopy.py` was not edited because `w_tiebreak_draws`, which imports it, is
running). `s26/results/w_selfcopy_bound.json` is regenerated with this script's provenance and
the original provenance kept under `provenance_original`; every number comes from
`s26/results/w_selfcopy_endpoint.json`.

1. `signed_bounds_gated/C_envelope_fold_ci/gain` added: leaked-model mean of (arm - sel) minus
   clean, +0.5101, fold CI [+0.3698, +0.6838], MDE 0.2548, 5/5 folds; bound (2/60) x 0.6838 =
   0.0228. `verdict/gain` now reads MINOR at 0.0228, as L44's Part D table says; the artefact and
   the ledger agree.
2. `signed_bounds_gated/C_envelope_per_target` added, every basis: (2/60) x |leaked-mean minus
   clean| at the worst single target, at the p95 target, and over (target, model) pairs. The
   Adversary's numbers reproduce exactly:

       basis   mean-CI limit   p95 target   worst target (pdb)     worst (target, model) pair   class under every reading
       arm       0.0277          0.0830       0.1512 (2BP4, 4.54 A)   0.1519                       MINOR
       gain      0.0228          0.0816       0.1232 (2NBC)           0.1261                       MINOR
       cloud     0.0280          0.0824       0.1587 (2MQ2)           0.1604                       MINOR
       sel       0.0479          0.1154       0.1944 (9KAR)           0.1944                       crosses 0.170 at the worst target

   So the class MINOR holds under every reading of the envelope on the built chain, the point
   cloud and the paired gain (the benchmark's headline contrast), and on the SELECTION basis
   (the benchmark's `shipped` argmin arm, the basis most sensitive to memorisation: L44 Part C
   -1.21 A) the worst-single-target reading is 0.194 A, above the pre-registered 0.170 line.
   Stated as the pre-registered rule gives it: for the benchmark's argmin mean alone, under the
   over-bound of own-native training and assumption A2, a single benchmark self-copy behaving
   like 9KAR could move that mean by up to 0.19 A; the expected contribution on every basis stays
   0.002 A (dev proxy) to 0.05 A (mean envelope), and the paired gain, which is the benchmark's
   verdict, is MINOR under every reading (worst 0.123). The benchmark verdict (no validated gain,
   +0.0103 [-0.160, +0.180]) does not move.
3. Wording to carry everywhere, as L55 asks: "2/60 self-copies; dev-proxy price 0.002 A; own-native
   envelope 0.028 A (mean CI) to 0.151 A (worst target) on the built chain, under A2; MINOR under
   every reading on the built chain and the paired gain; the selection basis reaches 0.194 at the
   worst target; cannot move the benchmark verdict either way." Applied in
   `s26/agentW_FINDINGS.md` (sections 0, 2.5, 4, 5), `s26/PREREG_selfcopy_bound.md` addendum 2 and
   `s26/IDEA_selfcopy_proxy_bound.md`. L44 itself is not edited.
## L59 -- GOVERNOR v2 DIED AT 19:36 ON A WINDOWS FILE-REPLACE RACE; v2.1 RETRIES AND NEVER EXITS THE LOOP ON A TRANSIENT ERROR (2026-09-13 19:37, coordinator)

`s26/governor.py` exited with `PermissionError: [WinError 5] Access is denied` from
`os.replace(governor_state.json.tmp, governor_state.json)` at 19:36:xx (its background task's
output): a reader (a jobrun waiter polling the snapshot, or a lane's `--status`) held the target
open at the instant of the replace. Four jobs were registered (a1_build, ph_reject_chain,
w_tiebreak_draws, ph_branch_solutions) and ran unsupervised for about a minute at 75% RAM; none
was harmed. Fix: `write_json_atomic` retries the replace eight times with a short back-off and
logs a skipped snapshot instead of raising; the main loop wraps the sample in a try/except that
logs and retries on the next tick. Restarted as v2.1 at 19:37. No production module touched.

---

## L60 -- REPORT CHECK: EVERY APPENDIX B NUMBER OF s26/REPORT.md RE-READ FROM ITS ARTEFACT (106 ROWS: 79 PASS / 2 MISSING PATH / 25 NOT FOUND ON THE FIRST RUN; 138 / 138 AFTER THE FIXES); STYLE 0 / 0 / 0; PARTS VII AND VIII FILLED TO L59; APPENDIX D RECONCILES docs/REPORT_S26.md (2026-09-13 19:50, lane E)

`s26/e_report_check.py` (committed b981b523, revised through 7db128c0) parses Appendix B of
`s26/REPORT.md` row by row, resolves every backticked path (`file :: key, key`, `file:line`,
`file:a-b`, `tests/x.py::test`, ledger entries named after a ledger path), checks existence, and
checks every number quoted in the row against the artefact: JSON leaves under the stated keys
(list means, container lengths, `[field=value]` filters, means and counts over `rows[*]/x`,
numbers inside string leaves), text lines or ledger entries by literal or by tolerance (half a
unit in the last quoted digit; a percentage as x and x/100), and `derived: expr = value` cells
evaluated with every four-decimal literal itself required to be in the row's artefacts. Rule 1:
`s9/final_report.json` and anything named benchmark are checked for existence only (two rows,
the benchmark means, are recorded as unverifiable behind Rule 1). `--style` scans the whole report
for the six banned words, U+2013 / U+2014, and any RMSD sentence contrasting two numbers with no
basis label on the sentence or its paragraph. Output: the table on stdout and
`s26/results/e_report_check.json` (save_atomic provenance, report sha256). Repeatable as
`python s26/examine.py --report-check` (examine.py gained the flag; a9984267 did not parse and
6f3708dc repairs it).

Runs, all through jobrun (tag CPU, est-ram 0.3): `e_report_check` (55 s including the wait, peak
RSS 0.032 GB) on the draft as committed at fcfbaf1d: 106 rows, 79 PASS, 2 MISSING PATH, 25 NUMBER
NOT FOUND; `e_report_check2` and `e_report_check3` after the fixes: 109 / 109 PASS, style 0 / 0 /
0; `e_report_check4` on the Part VII/VIII/D state is queued behind the four registered jobs (the
job cap) and refreshes the artefact when a slot frees; the iteration runs between them ran
directly to a scratch output (5 s, 9 MB). Current state, 7db128c0: 138 rows, 138 PASS, style 0 /
0 / 0.

What the 27 failures were, and what was done (drop or correct, never a new number without an
artefact): 2 MISSING PATH (a root `README.md` resolved against the previous citation's directory,
parser fixed; `geo_pauli_v1_rawonly.json` cited without its `s13/results/` directory). 25 NUMBER
NOT FOUND: 10 wrong or incomplete key paths (`cells/25` -> `rows[*]/cells/25/band_best`;
`rows[*]/MASS1.0, MASS0.1, MASS0.0` -> three full paths; `rmsd_vqe_sel, medoid128` ->
`arms/vqe_LFO/sel, arms/medoid128/sel`; `singularity` -> `summary/singularity`; the
`phys_landscape` summary keys; `concentration` -> `[model=amber]` / `[model=legacy]` for the 25 /
41 table counts; `q_alpha`, `q_verify` and `c_land_null` narrowed to their keys; `n_rows` added
where a row quoted "of 126"); 6 derived cells declared (9450 = 126 x 75; 5057 = 163 + 2627 + 2267;
96.8%; 78 = 0.6190 x 126; 42 = 2 P; 8128 = dim so(128); 12.6% = 1025 / 8128; the prior-ladder slope
and gain; +0.0207 = 3.2355 - 3.2148); 5 citations moved to the document that carries the number
(`s25/agentPHYS_FINDINGS.md:57-58, :357, :379` for 40/126, 462/500, +0.054;
`s25/agentQ_FINDINGS.md:375-376` for 45.3% and the 64.7th percentile; `docs/FINDINGS.md:3160` for
+0.994; `core/quantum.py:42` for +0.655634; `docs/FINDINGS.md:2751-2760` for 1.386); 4 numbers
replaced by the artefact's own (0.36 to 0.88 -> 0.358 / 0.827 / 0.882, `s20/LEDGER.md` L6; 2.6e-4 ->
2.239e-4 against the 8.66e-4 quantisation bound, `s25/LEDGER.md` L10; 699 -> 725 modules,
`module_map.json :: n_modules`; 58.6% re-sourced to `ph_reject_census.json ::
summary/per_threshold/1e4/pool_frac_over`); 1 number DROPPED: the relaxation cost's interval
"+0.0207 [+0.0143, +0.0276]" quoted from the state brief has no artefact (the brief does not
contain it either); the mean +0.0207 is kept as a derived difference of two stored means and the
paired interval now comes from L39 (+0.0207, fold CI [+0.0154, +0.0290]). One S19 number (2.66 per
75) re-cited to `s19/LEDGER.md` L12. Style: the first, broad detector flagged 33 sentences; 24
were not RMSD contrasts (correlations, free energies, variances, Pauli weights) and the detector
now skips those; 9 real RMSD contrasts lacked a basis on the sentence and are labelled (S12,
S13, S14, S15 and S21 numbers in Part VI, one S16 contrast in IV.4, the S8-instrument residual in
V.7, the SPSA contrast in V.10). Banned words: one ("leverage", VI.17) removed at 92d559bb;
dashes: none.

Appendix D reconciles `docs/REPORT_S26.md` and its summary (L42 item 1; read only, not edited or
committed): eleven artefact-sourced items taken with the same paths and the docs file cited as
the pointer (the leaderboard spread, the seven configurations against the incumbent, the
calibration means, the ladder's first rung, the binning and the 51-arm null, the Gibbs falsifier,
the S5 cosine and the consolidation audit's cosines, the entangler-deletion source, the
end-to-end speed-up, B1's numbers); one disagreement decided by the artefact (5.6e-17 -> 5.551e-17,
mine corrected); five apparent disagreements that are the same quantity on two paths or two
criteria (3.2126 / 3.2148 and 0.18242 / 0.18198, `compare_` / `baseline_tuning126.json`, 13/60 at
>= 0.6 identity against the 2/60 verbatim self-copies, the "100 of 125" citation, the blind
pipeline on 126 against 108); not taken: their external sources (arXiv, an author page), the
supervisors' names and affiliation, their Part X presentation guide (lane PR's), their memory-file
citations, and the two benchmark target ids they print from `docs/FINDINGS.md:4571`.

Parts VII and VIII now carry L35, L38, L39, L43, L44, L52, L53, L56 / L57 and the Adversary's
L45 to L49, L54, L55, L58, each row with its ledger entry and artefact; the C2 anchor's basis
(the rebuild built chain 3.2126, L57) is stated where the ladder is discussed; the leak wording of
L58 is carried in II.1, VII.4 and VIII.2. Next: A1 / A3, B2 / C2 rungs, the steric reject on the
built chain, C3 stage 2, the tie-break floor and the remaining tournament entries as they land;
Appendix C at the close.

---

## L61 -- LANE PR: THE DECK IS BUILT FROM ARTEFACTS (11 SLIDES; SLIDES 1-7, 9, 10, 11 FILLED; SLIDE 8 A PENDING PLACEHOLDER WITH THE A2/A4 FACTS); 237 REGISTERED NUMBERS, EACH WITH ITS PATH; FOUR NUMBERS TYPED FROM THE CLAIM LEDGER, SEVEN KEPT OFF THE SLIDES (2026-09-13 19:55, PR)

`vqe_research_overview.pptx` (repository root) did not exist (L2), so it is BUILT by `s26/pr_build_deck.py`
(python-pptx 1.0.2; dark theme; one accent; figures on white plates at 190 dpi; title, body and notes on every
slide) from `s26/pr_values.py` (every number read from its artefact at build time; the record is
`s26/pr_values.json`, 237 tokens with value, path, basis, status), `s26/pr_notes.md` (the spoken text with
`{TOKEN}` placeholders; each notes frame ends with a SOURCES list) and `s26/pr_figures.py` (the two
ORACLE-superposed CA overlays for 1S9Z = T030 and 9KAR, their RMSDs recomputed through `s12.instrument` and
equal to the record at 1e-9: 0.181981 / 7.437696; the S25 width and depth sweeps rebuilt from
`s25/results/q_plateau.json`; the accuracy ladder with the basis on every bar). Lane Q's A2 and A4 figures sit on
slides 6 and 8. Commits 71efb76d (a), a6738b35 (b), e93fc98c (c), cb90ffb3 (d placeholder). Verification
(`s26/pr_verify.txt`, pasted into `s26/PRESENTATION_CHANGES.md`): 11 slides; 0 U+2014, 0 U+2013; 0 banned words;
spoken words 234 / 248 / 244 on slides 8 / 9 / 10 (limit 250). Build peak RSS 0.12 GB, run directly.

Basis rule (L28 item 3, L29 item 2) applied throughout: the suite, the random-75 3.4251 and the Legacy / AMBER
verdicts are named point cloud; 3.2148 built chain; the C2 anchor the rebuild basis 3.2126 (L57); the S25
quantum contrasts the selection basis. The Adversary's caveats are in the wording: +0.0111 with the Type-M flag
and "a validity step on 124 of 126" (L46); the product circuit grounded in `q_dla.json` and no CI claimed on the
-0.302 vs -0.243 slopes (L47); the steric reject with the flag and "tail-carried" (L54); the strain signal as a
calibration flag, never a gain (L53); the benchmark caveat in L58's words on slide 4's notes.

Not sourced to a results artefact and labelled as such: the benchmark +0.0103 [-0.1596, +0.1803] 31W/29L (typed
from claim C06; `s9/final_report.json` not opened; the two means SOURCED_BY_TEST), the 4.4 GB headroom (L13,
a governor reading), the 11-rung count (the PREREG), the MDE factor (the contract). Kept off every slide: the
+0.0030 benchmark leak price (C27; the caveat uses lane W's bound instead), the z_moment triple (C34), the 0.524
cosine (C24), the 355/13 count (C35), the lost ESM -0.288 (L11), the S13 Pauli mean weights 2.236 / 3.015 (not
reproduced by me from `geo_pauli.json`; the median measured/predicted ratio 0.9969 over 104 cells is what slide
9 quotes), and the S12 operator-law coefficients. The phi statistic is on slide 2 with
`s26/results/a_c26_phi_mae.json` (L32). Ledger numbers recomputed on the way and found identical: L14, L26,
L35, L46, L53. Slide 11's direction line is DRAFT until the coordinator's verdict entry. Slide 8 waits for
`s26/PROPOSAL_A.md`; the builder swaps the placeholder for the proposal slide when the file exists.

---

## L62 -- LANE P, C2 RUNG NOESM: REMOVING THE ESM BLOCK COSTS +0.208 A ON THE BUILT CHAIN (5/5 FOLDS, 1.01x MDE, TYPE-M ZONE) AND +0.330 A ON SELECTION, REPRODUCING S7-11's LOST -0.288 TO -0.34; AND ITS gam_eff IS POSITIVE WHILE IT IS WORSE (2026-09-13 20:01, lane P)

Artefacts: `s26/results/p_ladder_noesm_s0.json` (126 rows, complete), `s26/results/p_ladder_report_noesm_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_noesm2`: exit 0, 521 s, peak RSS 0.305 GB (`s26/jobs_done/p_eval_noesm2.json`); fold models `s26/models/p_ladder/noesm_fold{0..4}_s0.pt` (42-d physicochemical pair features, the shipped MLP and trainer; `p_train_noesm` 1,743 s, peak 0.63 GB). Reading: (1) the S7-11 figure whose artefact is lost (L11) is re-measured on the same instrument: the ESM block is worth -0.330 A [fold -0.448, -0.212] on selection, 5/5 folds, against S7-11's -0.288 [-0.484, -0.092] for pca32 and -0.34 for pca128; (2) on the built chain it is worth -0.208 A, 5/5 folds, at 1.01x its MDE (Type-M zone; the magnitude is uncertain, the sign is not), so the ESM channel survives the pipeline, which S17 L23's 'filter, not discriminator' did not settle for the readout; (3) the medians (+0.044 arm, +0.071 sel) are far below the means and the worst target is +3.65 A (8T61): the loss is concentrated on a few targets, but the uniform-effect null puts the drop-top-10 at its 50th percentile, so it is not a concentration artefact; (4) gam_eff is POSITIVE (+0.116 prob-space at cos 0.25; +0.343 loc-space at cos 0.34, amplitude ~1) for a rung that is WORSE at 5/5 folds and has worse MAE (2.548 vs 2.339): a large move at low cosine projects positively onto the truth direction while adding more orthogonal error than it removes. This is S25 L12's caveat reproduced from the training side, on an achievable rung, and it settles how gam_eff will be read for every rung below: never without the cosine and the amplitude, and never as a prediction of the endpoint.

```
  noesm vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.4205 (med 3.4139)   b 3.2126 (med 2.9661)   n=126
    effect +0.2078   median +0.0440   SE 0.0733   MDE 0.2054   effect/MDE +1.01
    iid  CI95 [+0.0647, +0.3511]
    fold CI95 [+0.0986, +0.3942]   folds same sign 5/5   per-fold 0:+0.115 1:+0.126 2:+0.200 3:+0.573 4:+0.075
    47W/79L/0T   worst degradation +3.6485 (8T61)   p90 +0.9859   power 0.81  Type-M 1.12
    concentration: drop-top10 +0.3348 vs uniform-effect null p10/p50/p90 +0.2551/+0.3322/+0.4153 -> pctile 0.514
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.12x]
  noesm vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.2442 (med 3.1375)   b 3.0483 (med 2.8373)   n=126
    effect +0.1959   median +0.0576   SE 0.0684   MDE 0.1917   effect/MDE +1.02
    iid  CI95 [+0.0554, +0.3310]
    fold CI95 [+0.0775, +0.3814]   folds same sign 5/5   per-fold 0:+0.112 1:+0.193 2:+0.153 3:+0.547 4:+0.035
    47W/79L/0T   worst degradation +3.4270 (8T61)   p90 +1.0224   power 0.82  Type-M 1.11
    concentration: drop-top10 +0.3171 vs uniform-effect null p10/p50/p90 +0.2411/+0.3155/+0.3936 -> pctile 0.509
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.11x]
  noesm vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.7839 (med 3.6764)   b 3.4540 (med 3.4779)   n=126
    effect +0.3299   median +0.0713   SE 0.0924   MDE 0.2589   effect/MDE +1.27
    iid  CI95 [+0.1499, +0.5173]
    fold CI95 [+0.2122, +0.4475]   folds same sign 5/5   per-fold 0:+0.149 1:+0.203 2:+0.521 3:+0.443 4:+0.332
    48W/70L/8T   worst degradation +3.1238 (2LWS)   p90 +1.7459   power 0.95  Type-M 1.03
    concentration: drop-top10 +0.5026 vs uniform-effect null p10/p50/p90 +0.3881/+0.4977/+0.6130 -> pctile 0.524
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.03x]
  gamma-equivalent: gam_eff prob-space +0.1161 at cos +0.247 ; loc-space +0.3430 at cos +0.341 ; MAE 2.5483 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 5/5 ; verdict (arm): WORSE [TYPE-M ZONE: magnitude inflated ~1.12x]
```

---

## L63 -- LANE P, C2 RUNG CONLY: THE CONTACT HEAD ALONE (55-d, NO EMBEDDING BLOCK) SITS BETWEEN noesm AND THE SHIPPED PRIOR: +0.122 A ON THE BUILT CHAIN (0.62x MDE, 4/5 FOLDS, NOT MEASURED); AGAINST noesm IT BUYS -0.218 A ON SELECTION (1.00x MDE, 4/5) AND -0.086 A ON THE BUILT CHAIN (0.51x) (2026-09-13 20:12, lane P)

Artefacts: `s26/results/p_ladder_conly_s0.json` (126 rows, complete), `s26/results/p_ladder_report_conly_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_conly`: exit 0, 621 s, peak RSS 0.305 GB; models `s26/models/p_ladder/conly_fold{0..4}_s0.pt` (physicochemical pair block + the 13 ESM-2 650M contact-head columns; `p_train_conly` peak 0.69 GB). PREREG_B2's first size-axis point. Reading: the rung is UNDERPOWERED against the shipped prior (0.62x MDE on the built chain, 0.42x on selection; power 0.41 / 0.22), not null; its medians (+0.018 arm, 0.000 sel with 15 exact ties) say most targets are unmoved and the mean is carried by a few (worst 1CEK +3.70 A, the target whose verbatim self-window is BLOSUM rank 1, S24 L4). Isolated against noesm (paired, same code path): the contact head alone is worth -0.218 A on selection [fold -0.374, -0.037], 4/5 folds, exactly at its MDE (1.00x, Type-M zone), and -0.086 A on the built chain [fold -0.212, +0.010], 0.51x MDE, 2/5 folds. So of the ESM channel's -0.330 A on selection (L62), about two thirds is in the 13 contact-head columns and the rest in the 32-d embedding block (pca32, next); on the built chain the split is not resolved at n = 126. S17 L23 measured the contact head as a RANKER (filter, not discriminator); as a PRIOR INPUT it carries selection skill that the ranking measurement could not see, and that is consistent with the record: a better filter moves the argmin over K = 500 and barely moves a top-75 average. gam_eff/cos for the record: +0.109 prob at cos 0.24, +0.343 loc at cos 0.37, same shape as noesm (a large low-cosine move), MAE 2.437.

```
  conly vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.3349 (med 3.1110)   b 3.2126 (med 2.9661)   n=126
    effect +0.1223   median +0.0176   SE 0.0708   MDE 0.1983   effect/MDE +0.62
    iid  CI95 [-0.0113, +0.2620]
    fold CI95 [-0.0350, +0.2674]   folds same sign 4/5   per-fold 0:+0.115 1:-0.177 2:+0.221 3:+0.379 4:+0.079
    53W/73L/0T   worst degradation +3.7009 (1CEK)   p90 +0.8198   power 0.41  Type-M 1.55
    concentration: drop-top10 +0.2519 vs uniform-effect null p10/p50/p90 +0.1750/+0.2472/+0.3239 -> pctile 0.529
    VERDICT: NOT MEASURED (|effect| 0.1223 <= its own MDE 0.1983, 0.62x)
  conly vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.1532 (med 2.8753)   b 3.0483 (med 2.8373)   n=126
    effect +0.1049   median +0.0227   SE 0.0638   MDE 0.1788   effect/MDE +0.59
    iid  CI95 [-0.0201, +0.2354]
    fold CI95 [-0.0176, +0.2458]   folds same sign 4/5   per-fold 0:+0.096 1:-0.136 2:+0.196 3:+0.359 4:+0.026
    50W/76L/0T   worst degradation +3.1669 (1CEK)   p90 +0.6442   power 0.38  Type-M 1.61
    concentration: drop-top10 +0.2243 vs uniform-effect null p10/p50/p90 +0.1540/+0.2207/+0.2924 -> pctile 0.523
    VERDICT: NOT MEASURED (|effect| 0.1049 <= its own MDE 0.1788, 0.59x)
  conly vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.5656 (med 3.5772)   b 3.4540 (med 3.4779)   n=126
    effect +0.1116   median +0.0000   SE 0.0950   MDE 0.2660   effect/MDE +0.42
    iid  CI95 [-0.0761, +0.3032]
    fold CI95 [-0.0975, +0.3008]   folds same sign 3/5   per-fold 0:+0.279 1:-0.288 2:+0.373 3:+0.191 4:-0.000
    56W/55L/15T   worst degradation +4.4311 (1CEK)   p90 +1.5913   power 0.22  Type-M 2.15
    concentration: drop-top10 +0.2906 vs uniform-effect null p10/p50/p90 +0.1683/+0.2864/+0.4050 -> pctile 0.520
    VERDICT: NOT MEASURED (|effect| 0.1116 <= its own MDE 0.2660, 0.42x)
  gamma-equivalent: gam_eff prob-space +0.1093 at cos +0.241 ; loc-space +0.3425 at cos +0.368 ; MAE 2.4373 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 4/5 ; verdict (arm): NOT MEASURED (|effect| 0.1223 <= its own MDE 0.1983, 0.62x)
```

---

## L64 -- TOURNAMENT ITEM 4, tiebreak_noise_floor (W): THE PIPELINE'S OWN NOISE FLOOR FROM BLOSUM TIE-BREAKING AT THE K = 500 BOUNDARY IS 0.004 A ON THE 126-MEAN (BUILT CHAIN) AND 0.024 A AS THE PAIRED MDE BETWEEN TWO CONVENTIONS; EVERY HUNDREDTHS-LEVEL EFFECT ON THE RECORD IS INSIDE THE LATTER; THE PRODUCTION CONVENTION IS NOT DISTINGUISHABLE FROM A RANDOM DRAW (2026-09-13, W)

Pre-registered in `s26/PREREG_tiebreak_floor.md` (falsifiers and expectations in section 1, written
before the probe). Native-free half: job `w_tiebreak_draws` (exit 0, 5166 s under a four-job load,
peak RSS 0.111 GB; `s26/results/w_selfcopy_tiebreak_draws.json`, complete 126/126, production gate
top-75 == `sub` on 126/126). Gated half: job `w_tiebreak_endpoint` (exit 0, 25 s, 0.062 GB;
`s26/results/w_selfcopy_tiebreak_endpoint.json`); ledger numbers composed by
`s26/w_tiebreak_report.py` -> `s26/results/w_tiebreak_report.json`. Operator: on each target the
boundary tie class (every window sharing the 500th window's BLOSUM62 sum; median 115 windows, 54.5
of them inside the production pool) is re-drawn uniformly 8 times (`s15.seed.stable_rng`, seeds
("w_tiebreak", pdb, k)), the non-tied prefix unchanged, everything downstream identical (shipped
score, top-75, medoid-frame average, ramah 0.3 multi-start projection). Basis stated on every
line; the production chain here is the rebuild basis 3.2126 (L9, L57). Negative = the first arm is
better.

**Membership (native-free).** A random tie-break replaces 3.57 of the 75 averaged members on
average (median 3, max 15; one target never changes) and leaves the argmin unchanged on 91.5% of
(target, draw) cells: the score's top-75 is drawn mostly from the non-tied prefix. The registered
"about 8 of 75" was an over-estimate. The chain moves 0.168 A per draw in the median (triangle
bound, p90 1.08 A), 7x more than its RMSD-to-native changes (below): the moves are mostly
orthogonal to the native, as in L44.

**The floor, per basis (8 draws x 126 targets):**

    basis   m_tie (sd of the 126-mean over 8 draws)   mean over draws   production   s_tie median / p90   paired MDE between two draws: median / max over 28 pairs   max |paired effect| between draws
    arm     0.0039                                     3.2152            3.2126       0.0227 / 0.1290      0.0236 / 0.0321                                              0.0108
    cloud   0.0017                                     3.0532            3.0483       0.0126 / 0.0408      0.0100 / 0.0124                                              0.0046
    fit     0.0032                                     3.2099            3.2052       0.0148 / 0.0745      0.0157 / 0.0187                                              0.0104
    sel     0.0136                                     3.4741            3.4540       0.0000 / 0.0901      0.0466 / 0.0714                                              0.0398

Per target on the built chain: s_tie above 0.1 A on 21/126, above 0.5 A on 0/126; the range over
the 8 draws median 0.065 A, p90 0.362, max 0.713 (1D6X, the `pool_gate = WARN` target of L7/L9,
whose emission sits on a projection-branch flip). m_tie from 8 draws carries a relative SE of about
27% (sd of 8 numbers).

**Falsifier (PREREG section 1).** The claim is FALSIFIED if m_tie < 0.002 AND the paired MDE < 0.010
on the built chain: m_tie 0.0039 and paired MDE 0.0236, so not falsified. CONFIRMED if any recorded
effect is below the paired MDE between two draws: all nine listed are (`w_tiebreak_report.json ::
recorded_effects_vs_floor`): C3's production relaxation +0.0207 (L39) and its +0.0111 against the
matched random move, the S24 restrained relaxation -0.022, the ORACLE functional weight 0.015 (S24
L16), the all-atom reranking +0.004, C27's +0.0004 / +0.0018 (L44), the fd-vs-exact projection
difference 0.012 (L19), L24's predicted +0.013. Three of them (C27's two prices and the all-atom
reranking) are also below 2 x m_tie = 0.008 A, the spread of the 126-mean that the convention alone
produces. Registered predictions: m_tie 0.003 to 0.010 (measured 0.0039, holds); paired MDE 0.02 to
0.05 (0.0236, holds); median s_tie 0.03 to 0.08 (0.0227, just below); ~8 of 75 replaced (3.6, below).

**What the floor means, and what it does not.** It is the size of effect an unstated convention
produces on its own: two runs of the same pipeline that differ only in how the boundary ties are
broken differ on the 126-mean by an sd of 0.004 A (built chain) and can be told apart at 80% power
only above 0.024 A. It applies to any contrast whose two arms do NOT share the tie-break: a re-run
after a change of corpus order (the pinned `pdbs/` file set, state brief section 3), a different K,
a re-retrieval, a comparison across instruments or sprints (S9's 0.004 to 0.005 A residual between
the stable and the plain argsort, `docs/FINDINGS.md:1526, 7261`, is exactly this floor's effect on
the mean). It does NOT apply to a paired contrast whose two arms share the pool and its tie-break
(C3's relaxation, the functional weight, the reranking were all measured that way, and their
paired SEs already exclude this noise). So the sentence for the report is: "a hundredths-level
effect is real only as a paired contrast with the tie-break held fixed; quoted across runs or
against another instrument it is inside the pipeline's own convention noise (0.024 A at n = 126)."

**The production convention against random draws (registered secondary; `ST.fmt` verbatim):**

  production convention minus mean over 8 draws [arm], n=126
    a 3.2126 (med 2.9661)   b 3.2152 (med 2.9625)   n=126
    effect -0.0026   median -0.0004   SE 0.0062   MDE 0.0174   effect/MDE -0.15
    iid  CI95 [-0.0142, +0.0098]
    fold CI95 [-0.0155, +0.0156]   folds same sign 4/5   per-fold 0:-0.010 1:+0.031 2:-0.001 3:-0.004 4:-0.022
    64W/61L/1T   worst degradation +0.3828 (2LNG)   p90 +0.0577   power 0.07  Type-M 5.64
    concentration: drop-top10 +0.0099 vs uniform-effect null p10/p50/p90 +0.0023/+0.0095/+0.0171 -> pctile 0.534
    VERDICT: NOT MEASURED (|effect| 0.0026 <= its own MDE 0.0174, 0.15x)
  production convention minus mean over 8 draws [cloud], n=126
    a 3.0483 (med 2.8373)   b 3.0532 (med 2.8162)   n=126
    effect -0.0048   median -0.0011   SE 0.0032   MDE 0.0088   effect/MDE -0.55
    iid  CI95 [-0.0112, +0.0011]
    fold CI95 [-0.0089, -0.0006]   folds same sign 4/5   per-fold 0:-0.011 1:-0.007 2:+0.003 3:-0.008 4:-0.003
    70W/55L/1T   worst degradation +0.1433 (7YFS)   p90 +0.0133   power 0.34  Type-M 1.71
    concentration: drop-top10 +0.0019 vs uniform-effect null p10/p50/p90 -0.0013/+0.0017/+0.0050 -> pctile 0.524
    VERDICT: NOT MEASURED (|effect| 0.0048 <= its own MDE 0.0088, 0.55x)
  production convention minus mean over 8 draws [sel], n=126
    a 3.4540 (med 3.4779)   b 3.4741 (med 3.4796)   n=126
    effect -0.0201   median +0.0000   SE 0.0096   MDE 0.0269   effect/MDE -0.75
    iid  CI95 [-0.0392, -0.0033]
    fold CI95 [-0.0392, -0.0018]   folds same sign 4/5   per-fold 0:+0.009 1:-0.000 2:-0.055 3:-0.028 4:-0.024
    13W/6L/107T   worst degradation +0.2561 (6F3V)   p90 +0.0000   power 0.55  Type-M 1.34
    concentration: drop-top10 +0.0065 vs uniform-effect null p10/p50/p90 -0.0009/+0.0052/+0.0106 -> pctile 0.618
    VERDICT: NOT MEASURED (|effect| 0.0201 <= its own MDE 0.0269, 0.75x)

The corpus-order convention is not distinguishable from a random draw on the built chain (0.15x
MDE) or the point cloud (0.55x); on the selection basis it is 0.020 A better than a random draw at
0.75x MDE with the fold CI excluding zero and 107 ties: suggestive of the S25 finding that the
retrieval order is not neutral (rho(pool index, ORACLE RMSD) +0.054), and NOT MEASURED by the rule.
Nothing here is a lever: a tie-break cannot be chosen native-free, and the point of the number is
the floor. One representative draw-to-draw contrast (`ST.fmt` verbatim; the other 27 pairs are
summarised above):

  tie-break draw 0 minus draw 1 [arm], n=126
    a 3.2182 (med 2.9452)   b 3.2143 (med 3.0146)   n=126
    effect +0.0039   median +0.0007   SE 0.0084   MDE 0.0236   effect/MDE +0.17
    iid  CI95 [-0.0123, +0.0203]
    fold CI95 [-0.0041, +0.0193]   folds same sign 1/5   per-fold 0:-0.001 1:-0.001 2:-0.004 3:+0.034 4:-0.005
    54W/69L/3T   worst degradation +0.2871 (8T61)   p90 +0.0708   power 0.07  Type-M 5.16
    concentration: drop-top10 +0.0222 vs uniform-effect null p10/p50/p90 +0.0117/+0.0216/+0.0310 -> pctile 0.538
    VERDICT: NOT MEASURED (|effect| 0.0039 <= its own MDE 0.0236, 0.17x)

Replication: the 8 draws are 8 seeds; the 28 pairwise contrasts give the paired MDE's spread
(0.0236 median, 0.0321 max on the built chain); there is nothing positive to replicate at a second
seed, and the fold order does not enter (no fit). Power: the floor's own numbers ARE the power
statement. Not done: 16 draws (PREREG fork; the 8-draw m_tie carries a 27% relative SE, which the
verdict does not depend on); the AMBER-relaxed basis.

## L65 -- LANE P, C2 RUNG PCA32: REPRODUCTION GATE 2 PASSES AT THE ENDPOINT: THE RETRAINED pca32 EMITS THE SHIPPED PIPELINE'S ANSWER ON 126/126 TARGETS, ALL THREE BASES, 0 A DIFFERENCE (2026-09-13 21:01, lane P)

Artefacts: `s26/results/p_ladder_pca32_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pca32_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pca32`: exit 0, 536 s, peak RSS 0.306 GB; models `s26/models/p_ladder/pca32_fold{0..4}_s0.pt` retrained from the compact tables with the shipped recipe (`p_train_pca32`, peak 0.69 GB). The retrained posterior differs from the pinned one by 1.5e-5 per pair (fold-0 gate, L12) and that difference changes no argmin, no top-75 membership and therefore no emitted cloud or chain on any of the 126 targets: 126 exact ties on selection, point cloud and built chain (the FLAG on the concentration line is the degenerate all-zero case). This closes PREREG_C2's gate 2 in the strongest form: the ladder's trainer IS the production trainer at the endpoint, so any rung difference below is attributable to its inputs or architecture and not to retraining noise. It also fixes the ESM contrast: noesm (L62) is the shipped pipeline minus its ESM block and nothing else. gam_eff 0 at cos 0 (the identity), MAE 2.3386 = the shipped diagnostic.

```
  pca32 vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2126 (med 2.9661)   b 3.2126 (med 2.9661)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  pca32 vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0483 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  pca32 vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4540 (med 3.4779)   b 3.4540 (med 3.4779)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  gamma-equivalent: gam_eff prob-space -0.0000 at cos -0.001 ; loc-space -0.0000 at cos +0.010 ; MAE 2.3386 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 5/5 ; verdict (arm): NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
```

---

## L66 -- LANE P, C2 RUNG WIDE: THE CAPACITY RUNG (WIDTH 768, DEPTH 4, SAME INPUTS) IS NULL-TO-WORSE: +0.032 A ON THE BUILT CHAIN (0.25x MDE), +0.097 A ON SELECTION (0.48x MDE, FOLD CI ABOVE ZERO 5/5); MORE CAPACITY DOES NOT BUY A BETTER PRIOR (2026-09-13 21:11, lane P)

Artefacts: `s26/results/p_ladder_wide_s0.json` (126 rows, complete), `s26/results/p_ladder_report_wide_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_wide`: exit 0, 566 s, peak RSS 0.337 GB; models `s26/models/p_ladder/wide_fold{0..4}_s0.pt` (2.7x the parameters of the shipped 384x3, same 183-d inputs, corpus, epochs, loss; `p_train_wide` peak 1.03 GB). The built-chain contrast is UNDERPOWERED (0.25x MDE, power 0.11), so no gain of 0.13 A or more exists, and the point estimate points the wrong way. On selection the fold-clustered CI excludes zero on the harmful side at 5/5 folds while the effect is 0.48x its MDE: a real, small degradation, worth reporting as 'the wider head selects slightly worse', not as a measured effect size. MAE barely moves (2.348 vs 2.339): the wider net reaches the same conditional mean. This is S9-8's monotone-decreasing capacity curve (ridge beats a 512x3 net and a 1500-tree GBT) reproduced at the PRIOR rather than at the ranker: the information is not in the inputs, and adding capacity spends it on fitting the fragment distribution (S7-2). gam_eff +0.096 prob / +0.291 loc at cos 0.18 / 0.35: the third rung in a row with positive gam_eff and a worse endpoint.

```
  wide vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2444 (med 3.1364)   b 3.2126 (med 2.9661)   n=126
    effect +0.0317   median +0.0099   SE 0.0458   MDE 0.1282   effect/MDE +0.25
    iid  CI95 [-0.0566, +0.1201]
    fold CI95 [-0.0300, +0.1279]   folds same sign 3/5   per-fold 0:+0.018 1:+0.021 2:-0.056 3:+0.230 4:-0.027
    54W/72L/0T   worst degradation +2.4883 (6RRO)   p90 +0.3976   power 0.11  Type-M 3.50
    concentration: drop-top10 +0.1272 vs uniform-effect null p10/p50/p90 +0.0792/+0.1235/+0.1749 -> pctile 0.536
    VERDICT: NOT MEASURED (|effect| 0.0317 <= its own MDE 0.1282, 0.25x)
  wide vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0710 (med 2.8223)   b 3.0483 (med 2.8373)   n=126
    effect +0.0227   median +0.0123   SE 0.0427   MDE 0.1195   effect/MDE +0.19
    iid  CI95 [-0.0650, +0.1031]
    fold CI95 [-0.0500, +0.1161]   folds same sign 3/5   per-fold 0:+0.022 1:+0.060 2:-0.061 3:+0.199 4:-0.071
    56W/70L/0T   worst degradation +2.1217 (6RRO)   p90 +0.3972   power 0.08  Type-M 4.52
    concentration: drop-top10 +0.1131 vs uniform-effect null p10/p50/p90 +0.0692/+0.1096/+0.1545 -> pctile 0.544
    VERDICT: NOT MEASURED (|effect| 0.0227 <= its own MDE 0.1195, 0.19x)
  wide vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.5512 (med 3.5063)   b 3.4540 (med 3.4779)   n=126
    effect +0.0972   median +0.0000   SE 0.0725   MDE 0.2032   effect/MDE +0.48
    iid  CI95 [-0.0508, +0.2397]
    fold CI95 [+0.0487, +0.1587]   folds same sign 5/5   per-fold 0:+0.018 1:+0.099 2:+0.216 3:+0.099 4:+0.062
    55W/55L/16T   worst degradation +2.8395 (5MXS)   p90 +1.1442   power 0.27  Type-M 1.92
    concentration: drop-top10 +0.2280 vs uniform-effect null p10/p50/p90 +0.1348/+0.2230/+0.3121 -> pctile 0.529
    VERDICT: NOT MEASURED (|effect| 0.0972 <= its own MDE 0.2032, 0.48x)
  gamma-equivalent: gam_eff prob-space +0.0958 at cos +0.182 ; loc-space +0.2909 at cos +0.350 ; MAE 2.3480 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 3/5 ; verdict (arm): NOT MEASURED (|effect| 0.0317 <= its own MDE 0.1282, 0.25x)
```

---

## L67 -- LANE P, C2 RUNG PCA32F: REFITTING THE 32-PCA PER FOLD (NO GLOBAL-PCA LEAK) CHANGES NOTHING MEASURABLE: +0.018 A BUILT CHAIN (0.13x MDE), -0.032 A SELECTION (0.17x); THE SHIPPED GLOBAL PCA WAS NOT A LEAK WORTH ANYTHING (2026-09-13 21:34, lane P)

Artefacts: `s26/results/p_ladder_pca32f_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pca32f_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pca32f`: exit 0, 551 s, peak RSS 0.796 GB (the lean trainer's residue table, `s5/esmraw.npz`, is resident at eval); models `s26/models/p_ladder/pca32f_fold{0..4}_s0.pt` with per-fold bases `pca32f_pca_fold{0..4}.npz` fitted on training sequences only (`p_train_pca32f` peak 1.06 GB). Every contrast is far under its MDE (0.11x to 0.17x; power 0.06 to 0.08) with medians at or below 0.01 A: the rung is the shipped prior up to fitting noise. Two things this settles: the shipped `esm_pca.npz` (fitted once over the whole database, S7-11's 'weak leak') is worth nothing to the endpoint, so no ladder result is contaminated by it; and the isolate for pca128 is clean, since pca32f and pca128 differ only in the number of components. MAE 2.317 (the lowest so far) with no endpoint movement: S7-3's rule again. gam_eff +0.111 prob at cos 0.24, +0.348 loc at cos 0.40.

```
  pca32f vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2307 (med 3.0813)   b 3.2126 (med 2.9661)   n=126
    effect +0.0180   median -0.0092   SE 0.0491   MDE 0.1376   effect/MDE +0.13
    iid  CI95 [-0.0806, +0.1111]
    fold CI95 [-0.0725, +0.0970]   folds same sign 4/5   per-fold 0:+0.015 1:+0.103 2:+0.043 3:+0.125 4:-0.147
    70W/56L/0T   worst degradation +1.8708 (9BAF)   p90 +0.6142   power 0.07  Type-M 6.46
    concentration: drop-top10 +0.1163 vs uniform-effect null p10/p50/p90 +0.0582/+0.1142/+0.1710 -> pctile 0.518
    VERDICT: NOT MEASURED (|effect| 0.0180 <= its own MDE 0.1376, 0.13x)
  pca32f vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0347 (med 2.9544)   b 3.0483 (med 2.8373)   n=126
    effect -0.0136   median -0.0119   SE 0.0453   MDE 0.1270   effect/MDE -0.11
    iid  CI95 [-0.1023, +0.0752]
    fold CI95 [-0.1099, +0.0751]   folds same sign 3/5   per-fold 0:-0.000 1:+0.098 2:-0.021 3:+0.093 4:-0.187
    71W/55L/0T   worst degradation +1.3544 (2N9M)   p90 +0.5068   power 0.06  Type-M 7.85
    concentration: drop-top10 +0.0852 vs uniform-effect null p10/p50/p90 +0.0314/+0.0820/+0.1313 -> pctile 0.532
    VERDICT: NOT MEASURED (|effect| 0.0136 <= its own MDE 0.1270, 0.11x)
  pca32f vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4224 (med 3.3752)   b 3.4540 (med 3.4779)   n=126
    effect -0.0316   median +0.0000   SE 0.0649   MDE 0.1819   effect/MDE -0.17
    iid  CI95 [-0.1541, +0.0917]
    fold CI95 [-0.1108, +0.0541]   folds same sign 3/5   per-fold 0:+0.084 1:-0.136 2:-0.092 3:+0.101 4:-0.099
    60W/49L/17T   worst degradation +2.1793 (2MAI)   p90 +0.7372   power 0.08  Type-M 4.91
    concentration: drop-top10 +0.1025 vs uniform-effect null p10/p50/p90 +0.0195/+0.0986/+0.1767 -> pctile 0.523
    VERDICT: NOT MEASURED (|effect| 0.0316 <= its own MDE 0.1819, 0.17x)
  gamma-equivalent: gam_eff prob-space +0.1105 at cos +0.237 ; loc-space +0.3479 at cos +0.400 ; MAE 2.3168 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 4/5 ; verdict (arm): NOT MEASURED (|effect| 0.0180 <= its own MDE 0.1376, 0.13x)
```

---

## L68 -- LANE Q, A1: qubit-ADAPT-VQE IN PLACE OF THE FIXED ANSATZ IS NULL AT THE REGISTERED THRESHOLD ON THE BUILT CHAIN (-0.014 / -0.022 A, 0.23x / 0.36x MDE, FOLD CIs SPAN ZERO); EVERY ADAPT ARM POINTS THE SAME WAY AND NONE CLEARS ITS MDE; ON THE 78 alpha = 1 TARGETS ADAPT REACHES THE PRODUCT GIBBS STATE THE FIXED CIRCUIT MISSES BY 0.90 NATS (2026-09-13, lane Q)

Pre-registered in `s26/PREREG_A1.md` (filed 09:00, before any endpoint existed). Build:
`s26/q_adapt.py --build --tag a1` (two shards killed by the host at 10:10 after 33 targets;
one governed process resumed from those checkpoints: `s26/jobs_done/a1_build.json`, 7,782 s,
peak RSS 0.383 GB). Label: `s26/jobs_done/a1_label.json` (25 s, 0.36 GB), run after L33.
Stats: `s26/results/a1_stats.json`, `s26/logs/a1_stats.log` (every ST.fmt block).
Per-target records: `s26/results/a1/<pdb>.json` (126). Every record reproduces the
production quantum arm bit-for-bit: `ca` and `q_ca` against
`bench_results/cache/464a0ddb5f283e04/<pdb>.json` at max |diff| 0.0 and the selection index
equal, 126 of 126 (key `cache_check`). Basis: built chain (`rmsd_q_synth`, the production
projection of the weighted average) on both sides of every contrast below; the selection
readout (`rmsd_sel`) is the s8-instrument secondary and is named where it appears.
Comparator: `fixed_zrank_it50`, the deployed selector (3.2280 A built chain, 3.3135 A
selection; the production top-75 arm `rmsd_arm` is 3.2148 A, the shipped argmin 3.4540).

### The two PRIMARY contrasts (built chain), their selection-basis secondaries, and the matched random control, verbatim

  rmsd_q_synth: adaptL2_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2142 (med 3.1307)   b 3.2280 (med 2.9926)   n=126
    effect -0.0138   median +0.0066   SE 0.0210   MDE 0.0588   effect/MDE -0.23
    iid  CI95 [-0.0562, +0.0258]
    fold CI95 [-0.0705, +0.0442]   folds same sign 3/5   per-fold 0:-0.116 1:-0.008 2:+0.008 3:+0.096 4:-0.035
    58W/68L/0T   worst degradation +0.5351 (2MK7)   p90 +0.2886   power 0.10  Type-M 3.68
    concentration: drop-top10 +0.0294 vs uniform-effect null p10/p50/p90 +0.0046/+0.0295/+0.0533 -> pctile 0.498
    VERDICT: NOT MEASURED (|effect| 0.0138 <= its own MDE 0.0588, 0.23x)

  rmsd_q_synth: adaptV_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2059 (med 3.0924)   b 3.2280 (med 2.9926)   n=126
    effect -0.0222   median +0.0012   SE 0.0217   MDE 0.0608   effect/MDE -0.36
    iid  CI95 [-0.0665, +0.0194]
    fold CI95 [-0.0854, +0.0424]   folds same sign 3/5   per-fold 0:-0.129 1:-0.068 2:+0.030 3:+0.095 4:-0.031
    61W/65L/0T   worst degradation +0.5816 (2N9M)   p90 +0.2586   power 0.18  Type-M 2.44
    concentration: drop-top10 +0.0229 vs uniform-effect null p10/p50/p90 -0.0018/+0.0224/+0.0469 -> pctile 0.509
    VERDICT: NOT MEASURED (|effect| 0.0222 <= its own MDE 0.0608, 0.36x)

  rmsd_q_synth: randH_adaptL2 - fixed_zrank_it50
    a 3.2617 (med 3.0966)   b 3.2280 (med 2.9926)   n=126
    effect +0.0337   median +0.0214   SE 0.0350   MDE 0.0982   effect/MDE +0.34
    iid  CI95 [-0.0354, +0.1008]
    fold CI95 [+0.0122, +0.0581]   folds same sign 5/5   per-fold 0:+0.010 1:+0.006 2:+0.054 3:+0.075 4:+0.026
    55W/71L/0T   worst degradation +1.1926 (1D6X)   p90 +0.4926   power 0.16  Type-M 2.58
    concentration: drop-top10 +0.1077 vs uniform-effect null p10/p50/p90 +0.0701/+0.1069/+0.1436 -> pctile 0.513
    VERDICT: NOT MEASURED (|effect| 0.0337 <= its own MDE 0.0982, 0.34x)

  rmsd_sel: adaptL2_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2698 (med 3.1903)   b 3.3135 (med 3.1524)   n=126
    effect -0.0437   median +0.0000   SE 0.0334   MDE 0.0936   effect/MDE -0.47
    iid  CI95 [-0.1103, +0.0213]
    fold CI95 [-0.0720, -0.0224]   folds same sign 5/5   per-fold 0:-0.100 1:-0.052 2:-0.028 3:-0.025 4:-0.018
    39W/29L/58T   worst degradation +1.7418 (8HVS)   p90 +0.0706   power 0.26  Type-M 1.96
    concentration: drop-top10 +0.0289 vs uniform-effect null p10/p50/p90 -0.0089/+0.0268/+0.0645 -> pctile 0.527
    VERDICT: NOT MEASURED (|effect| 0.0437 <= its own MDE 0.0936, 0.47x)

  rmsd_sel: adaptV_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2688 (med 3.1903)   b 3.3135 (med 3.1524)   n=126
    effect -0.0448   median +0.0000   SE 0.0343   MDE 0.0961   effect/MDE -0.47
    iid  CI95 [-0.1109, +0.0218]
    fold CI95 [-0.0874, -0.0066]   folds same sign 4/5   per-fold 0:-0.099 1:-0.106 2:+0.016 3:-0.025 4:-0.018
    43W/34L/49T   worst degradation +1.7418 (8HVS)   p90 +0.1779   power 0.26  Type-M 1.96
    concentration: drop-top10 +0.0294 vs uniform-effect null p10/p50/p90 -0.0106/+0.0272/+0.0650 -> pctile 0.531
    VERDICT: NOT MEASURED (|effect| 0.0448 <= its own MDE 0.0961, 0.47x)

### The registered falsifiers

"ADAPT is null" fires on both primaries: -0.0138 A (0.23x its MDE of 0.0588) and -0.0222 A
(0.36x its MDE of 0.0608) are inside +-0.5x MDE, the fold-clustered CIs [-0.071, +0.044] and
[-0.085, +0.042] span zero, 3 of 5 folds agree in sign, W/L is 58/68 and 61/65, and the
drop-top-10 statistic sits at the 50th and 51st percentile of the uniform-effect null (no
concentration). "ADAPT helps" does not fire; no replication is owed.

Power statement. The comparison resolves 0.059 to 0.061 A on the built chain (MDE = 2.8016
SE at n = 126). An effect smaller than that could exist unseen; the observed effects are a
third of it, and Gelman-Carlin power for a true effect equal to the observed is 0.10 and 0.18
(Type-M 3.7 and 2.4). For scale, the whole quantum synthesis against the classical top-75 arm
is +0.0133 A (`s26/results/q_mde_reference.json`), so the resolution is four times the size
of the component's own footprint. Null at the registered threshold; underpowered below 0.06 A.

### Every ADAPT arm, built chain vs fixed_zrank_it50 (effect, x MDE; all fold CIs span zero; 3-4/5 folds)

    adaptL2 adam_best  P7 / P14 / P21   -0.0128 (0.21x) / -0.0192 (0.33x) / -0.0138 (0.23x)
    adaptL2 lbfgs      P7 / P14 / P21   -0.0206 (0.34x) / -0.0187 (0.32x) / -0.0194 (0.33x)
    adaptV  adam_best  P7 / P14 / P21   -0.0128 (0.21x) / -0.0195 (0.32x) / -0.0222 (0.36x)
    adaptV  lbfgs      P7 / P14 / P21   -0.0206 (0.34x) / -0.0244 (0.39x) / -0.0245 (0.40x)
    controls:  fixed 750 steps -0.0035 (0.16x)   gibbs_T +0.0088 (0.11x)   uniform128 +0.0135 (0.14x)
               randH_fixed +0.0230 (0.24x, 4/5)   randH_adaptL2 +0.0337 (0.34x, fold CI [+0.012,+0.058], 5/5)
    alpha subsets: L2 P21 alpha=1 -0.0223 (0.24x, n=78), alpha=0.25 +0.0000 (0.00x, n=48);
                   V  P21 alpha=1 -0.0254 (0.28x),       alpha=0.25 -0.0169 (0.27x)

Twelve of twelve ADAPT arms are negative on the built chain and twelve of twelve on the
selection readout (-0.039 to -0.048 A, 0.42x to 0.47x MDE, fold CIs excluding zero on all
twelve, 45 to 59 exact ties per contrast; 58 of 126 targets emit the identical selected
candidate under ADAPT-L2-P21 and the fixed circuit). By the standing rule (clear the MDE AND
fold CI excluding zero) none is a result; it is the shape of S25's `VQE_LFO - argmin`
(-0.1405 at 0.68x, 5/5 folds): a consistent direction with an unmeasured magnitude. The
matched-entropy random control (same 128 candidates, same entropy, no information) is worse
than the fixed circuit by +0.034 with a fold CI above zero on 5/5 folds and is likewise below
its MDE. The 7-parameter arms (P7, one trained RY layer) carry the same effect as the
21-parameter ones on both bases.

### The property half, on the real targets (no native in any of these numbers)

    KL(p || Gibbs) on the 78 alpha = 1 targets   fixed circuit 0.9027 (max 0.984)   [S25's 0.902 reproduced]
                                                 ADAPT, both pools, both optimisers 0.0002 (max 0.0009)
    KL(Gibbs_zrank || product of its marginals)  mean 1.4e-4, max 7.9e-4 over 126 targets
    L-BFGS growth at alpha = 1                   stops with no operator selected on 78 of 78 targets
    ADAPT sequences, alpha = 0.25, pool L2       mean 13.7 distinct 2-local strings; 38 distinct sequences on 48 targets

The deployed selector's target at alpha = 1 is a product state on every real target; a
7-parameter RY layer reaches it and the 21-parameter fixed circuit stops 0.90 nats short.
Reaching it changes the emitted structure by -0.02 A, a third of the MDE. That is the one
sentence Proposal A earns, and it is negative: a better-optimised state on the deployed
Hamiltonian is not visible through the readout, S25's readout-slack finding reproduced by a
second ansatz family.

Verdict for Proposal A (`s26/PROPOSAL_A.md`): REPLACE. The record now holds both halves: the
Hamiltonian is a constant ladder whose optimum is a product state (nothing to grow), and
growing anyway moves the answer by a third of the MDE.
## L69 -- PROPOSAL A VERDICT ACCEPTED: REPLACE (NOT THE EXPECTED KEEP WITH EDITS); HOW SLIDES 8 AND 9 STAY DISTINCT (2026-09-13 21:40, coordinator)

The campaign prompt expected KEEP WITH EDITS for Proposal A ("ADAPT-VQE is the right tool for
the trainability question and the DLA measurement"). Lane Q's evidence (L27, L35, L68;
`s26/PROPOSAL_A.md`) says the mechanism is absent rather than weak: the deployed objective's
optimum is a product state on every target (KL to the product of marginals at most 7.9e-4
nats), ADAPT offered entangling operators grows none on any of the 78 alpha = 1 targets, the
fixed circuit's DLA is already the full so(2^n) from depth 2 so there is no expressivity gap
for an adaptive ansatz to fill, the grown circuits give no width-scaling argument (L35), and the
endpoint is null at 0.23x and 0.36x its MDE with fold CIs spanning zero and a stated
resolution of 0.06 A. The prompt's own rule governs: a proposal is not softened to survive. The
verdict REPLACE stands, subject to the Adversary's check of L68.

Structure of the deck after this verdict (for lane PR): Proposal A and Proposal B both point
at the trainability paper, so slides 8 and 9 must not say the same thing twice. Slide 8 is
"what we learned by letting the circuit grow" (the ADAPT measurements and the product-state
diagnosis, the reason an ansatz cannot matter here); slide 9 is "what we publish" (the
trainability paper: the S13 locality theorem and Pauli-spectrum prediction, the S25 width sweep
scoped to depth 3 by the DLA result, A2, A4, the product-state fact, with the venue-honest gap
statement). Lane P's Proposal B verdict, when it lands, adds the feasible-scale ladder facts
(B2) and the B3 characterisability result to slide 9's notes as the evidence that the original
B had no testable form on this machine.

---

## L70 -- ADVERSARY CHECK OF L68 (A1, qubit-ADAPT vs THE FIXED ANSATZ): STANDS WITH CAVEAT; THE REPLACE VERDICT (L69) IS SUPPORTED (2026-09-13, A)

Checked against `s26/results/a1/<pdb>.json` (126 records), `s26/results/a1_stats.json`,
`s26/q_adapt.py` and `s26/PREREG_A1.md`.

Confirmed clean.
- Readout: every arm goes through the production functions `pl.consensus_medoid`,
  `pl.average_weighted` and `pl.project` on its own `p` (`s26/q_adapt.py:849-851`); the
  comparator `fixed_zrank_it50` reproduces the S25 quantum-arm cache bit-for-bit (`cache_check`:
  `ca` and `q_ca` max abs 0.0, `sel_equal` True, 126 of 126). Basis built chain both sides.
- Leakage: the native enters only in `label_one` (`:887-911`, RMSD scoring), gated on "PHASE 0
  SIGNED OFF" (`:925`); the build half reads no native. Clean.
- Registered threshold: "ADAPT is null fires if both PRIMARY contrasts lie within +-0.5x their
  own MDE" (`s26/PREREG_A1.md:85`); measured 0.23x and 0.36x. Fires as registered.
- Power statement: complete and correct (MDE 0.0588 / 0.0608, Gelman-Carlin power 0.10 / 0.18,
  Type-M 3.7 / 2.4, "null at the registered threshold; underpowered below 0.06 A"). Fold CIs
  span zero, 3 of 5 folds, concentration at the null's 50th percentile. No replication owed.
- The product-state claim rests on persisted per-target rows: `product_diagnostics/zrank/
  kl_gibbs_to_product` (recomputed over 126: mean 1.40e-4, max 7.90e-4, as quoted) and
  `arms/fixed_zrank_it50/kl_to_gibbs` on the 78 alpha = 1 targets (mean 0.9027, max 0.9836, as
  quoted; ADAPT-L2-P21 2.6e-4 mean, 7.0e-4 max). Sourced.

Three caveats.
1. **"Twelve of twelve arms negative" is one observation, not twelve.** The twelve ADAPT delta
   vectors (built chain, per target) have mean pairwise correlation 0.955 (min 0.919, max 1.000)
   and their first principal component carries 96.0% of their variance; `ST.best_of_k_within`
   over the 13 arms (fixed + 12) gives a per-target oracle of -0.094 A whose split-half transfer
   is +0.006 (-6% of the oracle: NOT A SIGNAL). L68 rightly calls none of them a result; the
   phrase "a consistent direction" must be read with the 0.955 beside it, so that twelve
   correlated draws at 0.2-0.4x MDE are not taken for a trend.
2. **The parameter match is a budget match, not a realised one.** `max_params = 3n = 21` for the
   arms labelled P21, but `len(ops)` runs 14 to 21 for Adam and 0 to 21 for L-BFGS (`adapt/*/ops`
   in the records): on the 78 alpha = 1 targets L-BFGS grows nothing and the "P21" arm is the
   7-parameter RY layer, and Adam grows 7 to 14 strings. Nothing in the null depends on it (P7,
   P14 and P21 agree within 0.01 A on both bases), but "21 parameters both sides" should be
   stated as "the same 21-parameter budget; realised counts 7 to 21".
3. **The `gibbs_T` control's +0.0088 overall masks two opposite subsets.** On the 78 alpha = 1
   targets the exact Gibbs state is -0.0271 against the fixed circuit and ADAPT-L2-P21 is
   -0.0223; the two agree with each other (+0.0047, SE 0.0060, 0.28x MDE, fold CI [-0.0009,
   +0.0098]) as two states 2.6e-4 nats apart should. On the 48 alpha = 0.25 targets `gibbs_T` is
   +0.0671 (the Gibbs state is not the CVaR optimum there) and ADAPT is +0.0000. So the
   product-state diagnosis is internally consistent on the subset it applies to; and, noted for
   the record, 29 of 78 alpha = 1 targets differ by more than 0.01 A between two states 2.6e-4
   nats apart (max 0.206 A): the projection readout amplifies sub-milli-nat differences on some
   targets, the mechanism L64 measures as the pipeline's own noise floor.

Verdict: STANDS WITH CAVEAT. The REPLACE verdict for Proposal A (L69) is supported: the
endpoint is null at the registered threshold with a stated resolution of 0.06 A, the optimum is
a product state on every real target by persisted per-target rows, and no arm family changes
the emitted structure beyond a third of the MDE.

---

## L71 -- ADVERSARY CHECK OF L64 (THE TIE-BREAK NOISE FLOOR): STANDS WITH CAVEAT (2026-09-13, A)

- Operator and reads: the boundary tie class is re-drawn uniformly with `s15.seed.stable_rng`
  (8 seeds), the non-tied prefix untouched, everything downstream production; the draws are
  native-free (`load_blind`) and the endpoint is gated; production gate top-75 == `sub` on
  126 of 126. Clean.
- Numbers: m_tie 0.0039 A (sd of the 126-mean over 8 draws; W states its ~27% relative SE),
  paired MDE between two draws median 0.0236 A over the 28 pairs (max 0.0321), per-target s_tie
  median 0.0227, p90 0.129, above 0.1 A on 21 of 126; 3.57 of 75 members replaced per draw; the
  argmin unchanged on 91.5% of cells. Registered predictions held (m_tie, paired MDE) or were
  over-estimates (8 of 75 replaced; median s_tie 0.03-0.08). `s26/results/w_tiebreak_report.json`.
- The production convention against a random draw: 0.15x MDE on the built chain, 0.55x on the
  cloud, 0.75x on `sel` with 107 ties and the fold CI excluding zero -- correctly NOT MEASURED,
  correctly read as suggestive of S25's non-neutral retrieval order and not a lever.
- Replication: the 8 seeds are the replication; no fit, no fold order. Power: the floor is the
  power statement.

Caveat. The title's sentence "EVERY HUNDREDTHS-LEVEL EFFECT ON THE RECORD IS INSIDE THE LATTER"
is true only as W's body qualifies it: the floor is the noise between two RUNS that do not share
the tie-break (a corpus re-order, a re-retrieval, a cross-instrument or cross-sprint comparison).
Every one of the nine listed effects (C3's +0.0207 and +0.0111, S24's -0.022, the 0.015 ORACLE
weight, the +0.004 reranking, C27's two prices, L19's 0.012, L24's +0.013) was measured as a
PAIRED contrast with the pool and its tie-break held fixed, and their paired SEs (C3: 0.0034)
already exclude this noise; none of them is invalidated. The presenter must quote the report
sentence W wrote ("a hundredths-level effect is real only as a paired contrast with the
tie-break held fixed; quoted across runs or against another instrument it is inside the
pipeline's own convention noise, 0.024 A at n = 126"), never the title alone. Second, minor:
m_tie from 8 draws carries a 27% relative SE, so "0.004" is 0.003 to 0.005; the verdict does
not depend on it.

Verdict: STANDS WITH CAVEAT (the title needs the body's qualifier wherever it is quoted).

---


## L72 -- LANE P, C2 RUNG PCA128: 128 ESM COMPONENTS INSTEAD OF 32: +0.076 A ON THE BUILT CHAIN (0.43x MDE, 3/5 FOLDS), +0.025 A ON SELECTION (0.13x); S7-11's '-0.337 vs -0.288, INDISTINGUISHABLE' STANDS AND THE DIRECTION IS, IF ANYTHING, WORSE (2026-09-13 21:45, lane P)

Artefacts: `s26/results/p_ladder_pca128_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pca128_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pca128`: exit 0, 576 s, peak RSS 0.799 GB; models `s26/models/p_ladder/pca128_fold{0..4}_s0.pt` (567-d inputs, the S7-11 per-fold 128-component bases `s7/repr_cache/pca128_fold{0..4}.npz`; `p_train_pca128` peak 1.41 GB). Isolated against pca32f (same construction, 32 vs 128 components; paired, same path): +0.058 A on the built chain, SE 0.056, MDE 0.156, 0.37x, fold CI [-0.048, +0.176], 3/5 folds, NOT MEASURED. S7-11 ranked pca128 nominally best on selection (3.417 vs pca32's 3.466, inside its CI); on this instrument it is +0.025 A on selection with 12 exact ties and +0.076 on the built chain, both far under their MDEs. Power: the design can see a 0.18 A built-chain gain and none of 0.18 A exists. MAE 2.386 (worse than pca32f's 2.317): more components fit the fragment corpus's directions, not the peptides'. gam_eff +0.117 prob at cos 0.21, +0.374 loc at cos 0.39.

```
  pca128 vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2885 (med 3.2575)   b 3.2126 (med 2.9661)   n=126
    effect +0.0759   median +0.0158   SE 0.0630   MDE 0.1764   effect/MDE +0.43
    iid  CI95 [-0.0410, +0.1932]
    fold CI95 [-0.0691, +0.2157]   folds same sign 3/5   per-fold 0:+0.229 1:-0.173 2:-0.081 3:+0.256 4:+0.131
    57W/69L/0T   worst degradation +2.2800 (6RRO)   p90 +0.7311   power 0.23  Type-M 2.10
    concentration: drop-top10 +0.1958 vs uniform-effect null p10/p50/p90 +0.1236/+0.1917/+0.2650 -> pctile 0.528
    VERDICT: NOT MEASURED (|effect| 0.0759 <= its own MDE 0.1764, 0.43x)
  pca128 vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.1089 (med 2.9577)   b 3.0483 (med 2.8373)   n=126
    effect +0.0606   median +0.0227   SE 0.0570   MDE 0.1598   effect/MDE +0.38
    iid  CI95 [-0.0573, +0.1725]
    fold CI95 [-0.0632, +0.1822]   folds same sign 3/5   per-fold 0:+0.204 1:-0.143 2:-0.088 3:+0.218 4:+0.100
    54W/72L/0T   worst degradation +1.8200 (2LER)   p90 +0.5037   power 0.19  Type-M 2.36
    concentration: drop-top10 +0.1731 vs uniform-effect null p10/p50/p90 +0.1081/+0.1682/+0.2315 -> pctile 0.538
    VERDICT: NOT MEASURED (|effect| 0.0606 <= its own MDE 0.1598, 0.38x)
  pca128 vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4789 (med 3.5336)   b 3.4540 (med 3.4779)   n=126
    effect +0.0249   median +0.0000   SE 0.0710   MDE 0.1990   effect/MDE +0.13
    iid  CI95 [-0.1117, +0.1606]
    fold CI95 [-0.1169, +0.1456]   folds same sign 4/5   per-fold 0:+0.036 1:-0.228 2:+0.037 3:+0.009 4:+0.211
    53W/61L/12T   worst degradation +2.5698 (7JGX)   p90 +0.7897   power 0.06  Type-M 6.76
    concentration: drop-top10 +0.1708 vs uniform-effect null p10/p50/p90 +0.0828/+0.1676/+0.2526 -> pctile 0.515
    VERDICT: NOT MEASURED (|effect| 0.0249 <= its own MDE 0.1990, 0.13x)
  gamma-equivalent: gam_eff prob-space +0.1171 at cos +0.214 ; loc-space +0.3744 at cos +0.387 ; MAE 2.3861 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 3/5 ; verdict (arm): NOT MEASURED (|effect| 0.0759 <= its own MDE 0.1764, 0.43x)
```

---

## L73 -- LANE PR: SLIDE 8 REBUILT IN PROPOSAL FORM FROM PROPOSAL_A.md (VERDICT REPLACE, L68 / L69) WITH THE ADVERSARY'S L70 CAVEATS ON THE SLIDE; SLIDES 8 AND 9 KEPT DISTINCT (DIAGNOSIS / PUBLICATION); ONE WORDING OF PROPOSAL_A.md THAT THE A1 RECORDS STATE DIFFERENTLY (2026-09-13 21:55, PR)

`vqe_research_overview.pptx` rebuilt (`python s26/pr_build_deck.py --no-figures`, then `s26/pr_changes.py`):
slide 8 now carries Proposal A's two-minute script adapted from `s26/PROPOSAL_A.md` section 5 (248 spoken
words) with every number read from `s26/results/a1_stats.json` (primaries: 2-local pool -0.0138 A, SE 0.0210,
MDE 0.0588, 0.23x, fold CI [-0.071, +0.044], 58W/68L; Tang pool -0.0222, SE 0.0217, MDE 0.0608, 0.36x,
[-0.085, +0.042], 61W/65L; NOT MEASURED; resolution 0.06 A; all 12 ADAPT arms -0.0128 to -0.0245 at 0.21x to
0.40x; the exact Gibbs state +0.0088 at 0.11x) and from the 126 records of `s26/results/a1/` (KL(Gibbs ||
product of marginals) mean 1.4e-4, max 7.9e-4; on the 78 alpha = 1 targets the fixed circuit's KL to Gibbs
0.9027, max 0.984, the 7-parameter RY layer's 1.5e-4, max 7.9e-4). Lane Q's A4 figure stays on slide 8; slide 9
is retitled "Direction B, and what we publish: the trainability paper" and names its Sprint 26 additions (the
DLA census, the product-state target, the non-decaying grown circuits, the A1 null), so that 8 is the diagnosis
and 9 the publication as L69 asks. The Adversary's L70 caveats are on slide 8, recomputed from the same records:
the twelve arms' per-target delta vectors correlate at 0.955 on average (min 0.919), one observation not
twelve; the 21 parameters are a budget (the Adam primaries realise 21, L-BFGS 7 to 21); the Gibbs control is
-0.0271 on the 78 alpha = 1 targets and +0.0671 on the 48 alpha = 0.25 targets. A1's basis is named on the
slide: `rmsd_q_synth`, the production projection of the weighted average over the 128 candidates (3.2280 A for
the deployed selector), not the production top-75 arm 3.2148. Verification (`s26/pr_verify.txt`, pasted into
`s26/PRESENTATION_CHANGES.md`): 11 slides, 0 U+2014, 0 U+2013, 0 banned words, spoken words 248 / 247 / 244 on
slides 8 / 9 / 10. Registry 296 tokens (`s26/pr_values.json`).

One statement of `s26/PROPOSAL_A.md` (sections 3 and 5) and of L68 that the records word differently, for lane
Q and the Adversary: "ADAPT ... selects no entangling operator on any of the 78 alpha = 1 targets under
L-BFGS" / "declines them on all 78 targets" / "stops with no operator selected on 78 of 78". In
`s26/results/a1/*.json :: adapt/{V,L2}_lbfgs_zrank`, the growth halts by the eps = 1e-3 pool-gradient criterion
on 156 of 156 (pool, target) cells, but before the halt it ADDS operators on 60 of 78 targets (pool V; 0 added on
18) and 68 of 78 (pool L2; 0 on 10), at most 11, the first of them entangling strings (1A13, pool V: IYIZZZZ at
pool-gradient norm 1.29e-3, just above eps); each such addition lowers the free energy by at most 1.2e-4 nats
(median 9e-6) from a 7-rotation state already within 7.9e-4 nats of the Gibbs optimum. The substance of the
claim holds (the growth is inert: nothing the objective or the readout can see), and the deck says it in the
records' words with the counts; the literal "no operator selected" does not hold. L70's caveat 2 ("on the 78
alpha = 1 targets L-BFGS grows nothing") reads the same way on the records: realised P for the L-BFGS P21 arms
runs 7 to 21. Lane Q's file is not edited. Slide 11's direction line stays DRAFT until the coordinator rules on
all three verdicts; lane P's B and C verdicts will trigger one more rebuild.

---
## L74 -- STALL AT 91.9% RAM WITH EVERY JOB SUSPENDED; GOVERNOR v2.2 ADDS A STALL BREAKER; THE raw RUNG'S TRAINING IS STOPPED BY HAND (21:51, coordinator)

At 21:50 the box read 91.9% RAM (user load: Chrome 5.0 GB, three claude sessions 2.3 GB, VS
Code 1.7 GB) with all four registered jobs suspended (ph_reject_chain, p_train_raw_f1 at
0.87 GB, w_ensemble_clouds, a3_build) and CPU at 16%: the v2 band suspends above 93% and
resumes only below 90%, so with the user's own load parked in between nothing of ours could
resume. v2.2: if RAM sits at or above 90% for 180 s while jobs are suspended, the governor
kills the suspended job holding the most memory (CTRL_BREAK first, so it checkpoints) and logs
it; the owner relaunches from the checkpoint when the box frees. To clear tonight's stall at
once the coordinator stops p_train_raw_f1 by hand (fold checkpoints under
`s26/models/p_ladder/`; lane P resumes it only on the coordinator's word) and restarts the
governor as v2.2. The raw rung is the ladder's last training item and not on the path to any
verdict.

---

## L75 -- LANE Q CORRECTS L68 ON PR's L73 FLAG: AT alpha = 1 ADAPT's GROWTH IS INERT, NOT ABSENT. OPERATORS ARE APPENDED ON 60 TO 78 OF 78 TARGETS, ALL MULTI-QUBIT UNDER L-BFGS, AND THEY BUY AT MOST 1.2e-4 NATS AND LEAVE A PRODUCT STATE (KL <= 4.1e-4); THE VERDICT DOES NOT MOVE (2026-09-13, lane Q)

PR's L73 read `s26/results/a1/` correctly and I did not. Re-read of my own records
(`s26/results/a1/<pdb>.json`, keys `adapt.<pool>_<optimiser>_zrank.sequence`, `.trace`,
`.theta`, `.stopped`, and `arms.*.kl_to_product`), the 78 alpha = 1 targets:

    pool / re-optimiser    stop reason     targets with >= 1     appended strings     F(P=7) - F(final)      max |angle| of      KL(p || product
                                           string appended       total / weight >= 2  mean / max, nats       appended strings    of marginals) max
    V  / L-BFGS-B          eps 78/78       60 / 78               235 / 235            1.5e-5 / 9.9e-5        0.018 rad           3.3e-4
    L2 / L-BFGS-B          eps 78/78       68 / 78               393 / 393            1.7e-5 / 1.2e-4        0.018 rad           4.1e-4
    V  / Adam best         P = 21 78/78    78 / 78               1092 / 304           7.5e-4 / 8.6e-4        0.28 rad (1 target) 3.0e-4
    L2 / Adam best         P = 21 78/78    78 / 78               1092 / 642           7.4e-4 / 8.3e-4        0.11 rad            2.7e-4
    (alpha = 0.25, L2 / Adam best, for contrast: F(P=7) - F(final) 0.053 nats; KL to product 0.126 mean, 0.148 max)

What "selected" means in `s26/q_adapt.py` `run_adapt`: at each growth step the pool operator
with the largest |dF/dphi| at phi = 0 is APPENDED to the circuit and all angles are
re-optimised; every appended operator stays with whatever angle it receives. There is no
retention step. So PR's count (appended) is the right count and my "no operator selected"
was false; the correct statement is that the appended operators are inert.

Why the ideal-ladder runs (L27, L35: L-BFGS appends nothing at alpha = 1; the Adam-grown
sets are abelian) differ from the real targets: on the ideal ladder E is exactly affine in
the register index, the RY layer reaches the exact Gibbs state and every pool gradient is
exactly zero. On a real target tie-averaging makes E not exactly affine (KL(Gibbs || product)
up to 7.9e-4), the RY layer's residual gradient in some multi-qubit direction exceeds
eps = 1e-3, an operator is appended, re-optimisation gives it an angle of order 5e-3 rad,
the free energy moves by order 1e-5 nats, and the criterion then fires. That is inert growth.

CORRECTIONS. In L68, the line "L-BFGS growth at alpha = 1 stops with no operator selected on
78 of 78 targets" and the sentence "ADAPT with either pool selects no entangling operator on
any of the 78 alpha = 1 targets under L-BFGS" are RETRACTED and replaced by the table above.
In `s26/PROPOSAL_A.md` section 3 ("selects no entangling operator on any of the 78 alpha = 1
targets") and section 5 ("it declines them on all 78 targets; there is nothing to entangle")
are corrected by an appended addendum. In `s26/agentQ_FINDINGS.md` sections 0.1 and 1.3 the
same sentences are corrected by an appended note (section 9). The ideal-ladder statements in
L27 and L35 stand as written (they are about the ladder, where the appended count is zero).

What does not change: every endpoint number in L68 (the two primaries, the twelve arms, the
controls, the subsets); the product-state result (KL(Gibbs_zrank || product) <= 7.9e-4 on 126
targets; ADAPT reaches the Gibbs state to KL <= 9e-4 on the 78 alpha = 1 targets; the fixed
circuit stops 0.90 nats short); the A2 and A4 results; the verdict REPLACE. The mechanism
sentence becomes: at alpha = 1 the objective's optimum is a product state, and an adaptive
ansatz that is free to entangle appends operators worth less than 1e-3 nats and leaves the
state a product state; there is nothing for growth to do, and growing anyway moves the
built chain by a third of the MDE.

FINAL WORDING FOR SLIDE 8's NOTES (for lane PR), replacing any "grows nothing" / "declines
them" phrasing:
"At the deployed setting the optimum the circuit is asked to reach is a product state:
seven single-qubit rotations represent it exactly (KL 0.0002 nats). When ADAPT is offered
entangling operators it does append them, on 60 to 78 of the 78 targets, but they are inert:
together they lower the objective by less than 0.001 nats, their angles stay below 0.02
radians under L-BFGS, and the state remains a product state to 0.0004 nats. The deployed
21-parameter circuit stops 0.9 nats short of that same optimum, and reaching it exactly
moves the emitted structure by 0.014 to 0.022 angstroms, a third of what the comparison can
resolve. Sources: s26/results/a1/*.json (adapt.*.sequence, .trace, arms.*.kl_to_product),
s26/results/a1_stats.json, ledger L68 and L75."

## L76 -- LANE PR: L75's FINAL WORDING APPLIED TO SLIDE 8 (NOTES AND TEXT) AND TO SLIDE 9's SPRINT-26 LINE; NO "GROWS NOTHING" / "DECLINES THEM" / "SELECTS NO ENTANGLING OPERATOR" PHRASING REMAINS; EVERY L75 QUANTITY RECOMPUTED FROM THE A1 RECORDS; THE RETRACTION IS IN THE LANE'S CHANGE LOG (2026-09-13 22:25, PR)

Applied L75 (lane Q's reconciliation of this lane's L73 flag; the two L68 sentences retracted there;
`s26/PROPOSAL_A.md` addendum). Slide 8's spoken notes now carry L75's final wording with every number
read from the records: "seven single-qubit rotations represent it exactly, to 0.0002 nats"
(`s26/results/a1/*.json :: arms/adaptV_lbfgs_zrank_P7/kl_to_gibbs`, mean over the 78 alpha = 1 targets
1.5e-4); "it does append them, on 60 to 78 of the 78 targets" (`adapt/V_lbfgs_zrank/sequence` non-empty on
60, `L2_lbfgs_zrank` on 68, `*_adam_best_zrank` on 78); "together they lower the objective by less than
0.001 nats" (`adapt/*/trace`: F(P = 7) - F(final) at most 1.2e-4 under L-BFGS, 8.6e-4 under Adam);
"their angles stay below 0.02 radians under L-BFGS" (`adapt/*_lbfgs_zrank/theta` beyond the first 7:
max 0.0179 rad); "the state remains a product state to 0.0004 nats" (`arms/adapt*_P21/kl_to_product`
max 4.1e-4); "stops 0.90 nats short" (`arms/fixed_zrank_it50/kl_to_gibbs` mean 0.9027). The slide text
says the same in one bullet, and slide 9's Sprint-26 line scopes the A4 non-decay to the ideal ladder
(L27, L35 stand) with "inert growth on the real targets (L75)". The deck dump has no "grows nothing",
"declines them" or "selects no entangling operator" left; the notes cite the retracted phrase once, as
retracted. Verification PASS: 11 slides, 0 U+2014, 0 U+2013, 0 banned words, spoken words 249 / 247 /
244 on slides 8 / 9 / 10. Registry 302 tokens. `s26/PRESENTATION_CHANGES.md` row updated to name the
retraction and its artefact keys; `s26/agentPR_FINDINGS.md` section 7 is a dated change log (the flag
at L73, the retraction at L75, this application) for the Adversary's RETRACTIONS.md. No endpoint
number or verdict changed. Still pending: lane P's B and C verdicts (one more rebuild) and the
coordinator's ruling on all three directions for slide 11's DRAFT line.

---

(Clock correction, PR: L76's header time is 21:56, not 22:25; the same correction applies to the 22:25 lines in `s26/STATUS.md` and `s26/agentPR_FINDINGS.md` section 7. Nothing else in L76 changes.)

---
## L77 -- THE USER EXTENDS THE SPRINT: EVERY DEFERRED RUN, TEST AND DIAGNOSTIC GOES AHEAD; CLOSE AT ABOUT 04:30 PACIFIC, FINAL REPORT AT ABOUT 04:45 (2026-09-13 22:02, coordinator)

The user's instruction at 22:0x: run all the tests and additional diagnostics possible, push the
close to about 04:30 Pacific (2026-09-14), and deliver the final report with all detail at about
04:45. Nothing running is interrupted. The extra scope, by lane, under the same rules (prereg
before compute; governor; one honest est-ram per job; the 93% ceiling stands, so memory, which
is the user's own load tonight, sets the parallelism, and the stall breaker (L74) may kill the
fattest suspended job):

- I (resumed): the opt-in test tier `VERIFY_SLOW=1` (the 11 real skips of `s26/TEST_RUN.md`:
  OpenMM and pipeline checks, AMBER-tagged, two at most); every standalone audit under
  `verify/` re-run with its JSON compared to the tracked one; the results lab `--mode frozen`
  rebuild as a governed job with the four gates and the leaderboard numbers compared to the
  tracked `results/summary/` (3.2126 built chain, 3.0483 point cloud; the L7 WARN explanation
  appears); the final AST gate of every production `.py` against `ae86a124` (docstrings
  stripped); `python s26/examine.py` on the final tree; `s26/TEST_RUN.md` extended.
- Q: A3 to completion; the DLA dimension of the ADAPT-grown circuit at every growth step across
  the 126 A1 records (the A2 item the prompt asked for per growth step); bootstrap CIs on the
  A4 slopes (the L47 caveat); a second-seed and reversed-fold-order replication of A1's
  fixed-versus-ADAPT contrast as a robustness check even though it is a null.
- P: the remaining rungs (esm8m, mix, pairnet), then raw (retrain fold by fold when the
  governor shows 2.5 GB free), then the deferred coherence_penalised_training (1.25 GB) if the
  box allows; B3, C4, C5; length-stratified and FAIL18-stratified tables for every rung;
  replication for any rung clearing its MDE; the best-rung delivery file; the final
  PROPOSAL_B.md and PROPOSAL_C.md. `attn` (3 to 3.5 GB) only if the box empties.
- PH: the reject chain to completion; branch_select relaunched from its 41 cells; rotamer_relief
  B; `stage1 --rep` (the C3 replication) and `floor2`; C3 stage 2 on P's rung file; the
  heavy-atom validity axis of the relaxed emission.
- W: ensembling to completion; window_provenance; amber_prior_partner; Part B's retrains one at a
  time as memory allows; the two IDEA files found on the way get preregs and, if cheap, runs.
- A: every new entry checked within the hour; RETRACTIONS.md and DELIVERABLES_CHECK.md kept
  current; the final deliverables pass at about 04:00.
- E: the report grows with every result; the final pass (Parts VII, VIII, Appendix B, Appendix C
  = this ledger reproduced) starts at about 04:15 for a 04:45 close; the report check re-run last.
- PR: the last rebuild after the B and C verdicts and the coordinator's slide 11 ruling.
- Coordinator: the slide 11 ruling once the three verdicts are in; the `docs/FINDINGS.md`
  corrections-ledger update from RETRACTIONS.md; the sprint-close entry.

Eight lanes active (E, I, Q, P, PH, W, A, PR): the maximum.

---

## L78 -- THE STALL BREAKER'S FIRST KILL (p_eval_esm8m, 22:03); THE CTRL_BREAK GRACE SIGNAL DOES NOT CROSS CONSOLES ON WINDOWS, SO A GOVERNOR KILL IS A HARD KILL AFTER 25 s (2026-09-13 22:08, coordinator)

`s26/governor.log` 22:03:30: the v2.2 stall breaker (L74) fired as designed on the fattest
suspended job while RAM sat above 90%, `p_eval_esm8m`; the CTRL_BREAK it sends first failed
with `OSError(22, 'The parameter is incorrect')` (WinError 87): `GenerateConsoleCtrlEvent` can
only reach a process group attached to the caller's console, and the governor runs in its own
console. The job was terminated at 22:03:55 after the 25 s grace and lane P relaunches it from
its checkpoint. Consequence, stated for every lane: a governor kill is a hard kill after 25 s;
the loss is bounded by the job's own checkpoint granularity (per target or per fold in every
S26 job), and a killed job's partial rows are not a result until the relaunched run completes
with `complete: true`. Not changed tonight: a file-based stop request that jobs would poll is
the right fix and is noted for the next campaign in `s26/agentI_FINDINGS.md`'s hygiene list.

---

## L79 -- ADVERSARY CHECK OF THE C2 RUNGS L62, L63, L65, L66, L67, L72: ALL STAND; L66's "A REAL, SMALL DEGRADATION" ON SELECTION IS NOT LICENSED AT 0.48x MDE; THE LADDER'S EVENTUAL BEST RUNG IS AN ORDER STATISTIC (2026-09-13, A)

Common path checked once (`s26/p_ladder.py`): the phase gate is enforced in code (`:118`
`SIGNOFF`; `:64`); native quantities enter only in the read-after-selection half (`:513-522`,
"ORACLE-READ-AFTER-SELECTION"), after `score_target`; the top-75 is `argsort(kind="stable")`
(`:505`); `sel` is tie-averaged with `ST.argmin_tied` (`:517`); every rung goes through one code
path against the one anchor (`p_ladder_shipped_s0.json`), whose built chain is the rebuild basis
3.2126 (L57) and is bit-exact against the cache on selection and point cloud; gate 2 (L65) shows
the retrained pca32 emits the shipped answer on 126/126, so rung differences are inputs or
architecture, not trainer noise. iid and fold CIs on every contrast; concentration nulls
computed; medians beside means; power stated. Clean on the checklist.

Per rung:
- L62 noesm: WORSE at 5/5 folds on all three bases, but every magnitude is Type-M (1.01x /
  1.02x / 1.27x): the sign is measured, the "+0.208 A" and "+0.330 A" are upper bounds; L62
  says so. The harm is tail-carried (median +0.044, worst +3.65 on 8T61); the uniform-effect null
  does not flag it. "Reproducing S7-11's -0.288 to -0.34" is a reproduction of a Type-M number
  by a Type-M number; quote both with the flag. The gam_eff-positive-while-worse reading (S25
  L12 from the training side) is correct and important: gam_eff is never a prediction of the
  endpoint.
- L63 conly: correctly UNDERPOWERED (0.62x), not null; the "two thirds of the ESM channel is in
  the contact head" is a point-estimate split (-0.218 at exactly 1.00x MDE against -0.330 at
  1.27x, both Type-M); the built-chain split is unresolved and L63 says so. The worst target is
  1CEK (+3.70), the verbatim self-copy target: noted, not read into.
- L65 pca32: 126 exact ties on all bases; the degenerate all-zero concentration FLAG is the
  identity case. The strongest possible gate-2 pass.
- L66 wide: built chain UNDERPOWERED (0.25x). CAVEAT: the selection sentence "a real, small
  degradation, worth reporting as 'the wider head selects slightly worse'" is not licensed:
  +0.097 is 0.48x its MDE, the iid CI spans zero ([-0.051, +0.240]) and only the fold CI excludes
  it; by the standing rule (clear the MDE AND fold CI excluding zero) it is not a result and the
  word is "suggestive, not measured". The conclusion "more capacity does not buy a better prior"
  stands on the built-chain null and the flat MAE.
- L67 pca32f: every contrast at 0.11x to 0.17x; "changes nothing measurable" is correct at the
  stated MDEs (0.13 to 0.18 A) and closes the global-PCA leak question at that resolution.
- L72 pca128: 0.43x / 0.38x / 0.13x, all under MDE; the isolate against pca32f (+0.058, 0.37x)
  is the clean size-axis point; power stated (a 0.18 A gain would be seen).

Multiplicity, stated now for the ladder's close: seven rungs so far and four to come, all
paired against one anchor on one pool, so the per-rung deltas are correlated (the rungs share
the pool, the score and the projection); any "best rung" or per-target minimum over rungs must
go through `ST.best_of_k_within` with its k_eff and split-half transfer (PREREG_C2's own rule),
never the raw minimum. No rung is positive so far, so nothing is owed yet. Verdicts: L62, L63,
L65, L67, L72 STAND; L66 STANDS WITH CAVEAT (the selection-basis wording).

---

## L80 -- ADVERSARY CHECK OF L52 (conformational_identity_floor, ORACLE DIAGNOSTIC): STANDS WITH CAVEAT (n = 18; BOTH PAIRED CONTRASTS ARE TYPE-M) (2026-09-13, A)

Reads: the copy's CA-RMSD to the target's model-1 native at the shared segment and the pool's
ORACLE `rr`; every quantity reads the native and is labelled ORACLE DIAGNOSTIC; nothing selects;
dev targets only (lane I's L15 list, 18 targets, 22 partners; no benchmark sequence, confirmed
in L49 for the same script). Registered falsifier F5 (median above 1.5 A; FALSIFIED below 1.0)
holds at 2.908 A over 22 pairs (18% below 1.0 A, 27% below 1.5 A); the four cross-fold
self-copies reproduce S24 L4's four values. Deterministic (Kabsch on fixed coordinates; the
bootstrap is seeded). Power: n = 18 gives MDEs of 0.92 and 1.26 A on the two paired contrasts.

Caveat. Both paired contrasts are in the Type-M zone (copy minus pool best +0.968 at 1.06x MDE;
copy minus pool mean -1.443 at 1.14x); W labels them so ("sign clean, magnitude not a result").
The number for the report is the median 2.9 A, a descriptive statistic at n = 22, not the paired
effects; and the sentence "the 2.0 A target is below what a verbatim lookup reaches" is a
statement about this instrument's 18 targets with a verbatim relative, not a general law. W's
own "not done" (no length-matched non-verbatim control population) is the right missing arm;
the pool mean plays that role loosely. Verdict: STANDS WITH CAVEAT.

---

## L81 -- ADVERSARY CHECK OF L53 (strain_difficulty): STANDS WITH CAVEAT, PROVISIONAL UNTIL THE POOL-SPREAD CONTROL LANDS; "REPLICATED" MEANS THE CIs, NOT THE ESTIMATE (2026-09-13, A)

Checked against `s26/ph_strain.py` and `s26/results/ph_strain.json` / `ph_strain_rep.json`.
- Leakage: the four signals are the production cache's own relaxation scalars (`amber_moved`,
  `amber_e0`, `amber_e1`, `amber_strain_after`), native-free by construction; `rmsd_arm` enters
  only as the ORACLE label (`:84-90`); gated (`:74`); nothing is selected, weighted or tuned.
  Clean.
- Multiplicity: four pre-registered signals, Bonferroni (p < 0.0125); the falsifier fires on one
  (`moved`, partial rho +0.433, fold CI [+0.247, +0.588], 5/5 folds, permutation p < 0.00025).
  Confounds n and Rg partialled by residualisation; per-fold signs and a 4000-draw permutation
  null. Clean. The two 4/5-fold signals (log_e0, log_drop at +0.24, +0.25) are correctly reported
  as measured-but-failing-the-rule, and `strain_after` as null; the FAIL18 Fisher tests as null.
- Outliers: rho +0.41 without 9KAR and 2BP4 (the broken-bond emissions). Clean.
- Basis: built chain, ORACLE label, stated. Clean.

Two caveats.
1. **"REPLICATED" in the title means the bootstrap and permutation draws, not the estimate.**
   The Spearman is a deterministic function of 126 fixed rows; the registered replication
   (new seed, reversed order) can only re-draw the CIs and the permutation null, and it did
   (identical rho +0.433, fold CI [+0.248, +0.581]). It cannot test whether the signal
   generalises to new targets; that needs targets the record does not have (no fresh
   benchmark exists). Say "CIs replicated", not "result replicated".
2. **The prereg carried no pool-disagreement control.** The obvious native-free difficulty
   proxy is the top-75's own pairwise CA-RMSD spread (the members' disagreement, no native).
   If `moved` is that spread in disguise, L53's "first native-free quantity above 0.4" is a
   re-labelling. The Adversary's control `s26/a_strain_vs_spread.py` (rho of the spread with
   `rmsd_arm`; partial rho of `moved` given n, Rg AND the spread; job `a_strain_vs_spread`,
   queued behind the four-job cap at 22:3x) decides it; the verdict below is provisional until
   it lands and is confirmed or amended in a follow-up entry.

Verdict: STANDS WITH CAVEAT, provisional. As a calibration flag with the L53 wording ("a
calibration curve, never a gain") it may be quoted now; the novelty sentence waits for the
control.

---

## L82 -- CORRECTION TO THE ADVERSARY'S L70 CAVEAT 2 PER L75: L-BFGS DOES APPEND OPERATORS AT alpha = 1; THEY ARE INERT; THE CAVEAT'S CONCLUSION IS UNCHANGED (2026-09-13, A)

L70 caveat 2 said "on the 78 alpha = 1 targets L-BFGS grows nothing and the 'P21' arm is the
7-parameter RY layer". Per L75 (lane Q, from the records' `adapt.*_lbfgs_zrank.sequence`), under
L-BFGS operators are appended on 60 of 78 (pool V) and 68 of 78 (pool L2) targets, all
multi-qubit, worth at most 1.2e-4 nats, with angles at most 0.018 rad, leaving the state a
product to 4.1e-4 nats. My `len(ops)` reading (0 to 21 under L-BFGS) was consistent with that
and I mis-stated it as "grows nothing". Corrected wording of caveat 2: the parameter match is a
budget match (21) and the realised counts run 7 to 21 under both optimisers, with the L-BFGS
appended operators inert; nothing in the A1 null depends on the count (P7, P14, P21 agree within
0.01 A). Recorded in `s26/RETRACTIONS.md` R5 together with lane Q's own corrections.

---

## L83 -- GOVERNOR v2.3 ADOPTS SUSPENSIONS MADE BEFORE ITS START; THE LAUNCH CAP IS NOW A FILE AND IS RAISED TO SIX (2026-09-13 22:21, coordinator)

Two throttles found at 22:20 with the box at 78% RAM and 3.7 GB free: (1) `a3_build` had been
suspended by the v2.1 governor and was never resumed by v2.2, whose suspended stack starts
empty on a restart; v2.3 adopts any stopped registered job onto its stack on every tick and
logs ADOPT, so it resumes under the normal rule. (2) Eight launcher waiters (P esm8m2 and mix,
E report_check5, A strain_vs_spread, I slow_equivalence and slow_integration, Q a4_var_boot and
a2_dla_a1) sat behind jobrun's fixed cap of four with room on the box. `s26/jobrun.py` v2.3
reads the cap from `s26/launch_cap.json` on every wait tick (coordinator-only file); it is set
to 6 now. Waiters started before this edit keep the old fixed cap of four until they launch;
every new waiter uses the file. The memory rules (est-ram + 0.5 GB free; the 93% ceiling; the
stall breaker) are unchanged and still bind before the cap does.

---


## L84 -- TOURNAMENT ITEM 6, window_ensembling (P's idea, orphaned to W; THE MANDATORY TEST-TIME-ENSEMBLING DIRECTION): FIXED-K ENSEMBLING OF THREE BLOSUM KEYS IS -0.0005 A ON THE BUILT CHAIN (0.02x MDE) AND INDISTINGUISHABLE FROM A ZERO-INFORMATION RESAMPLE OF THE PRODUCTION SET (+0.0042, 0.14x); THE THREE CLOUDS SIT 0.14 TO 0.20 A APART AND THEIR ERRORS ARE PARALLEL, AS S23 L9 / S24 L3 PREDICT; A GAIN OF 0.028 A OR MORE IS EXCLUDED (2026-09-13, W)

Pre-registered in `s26/PREREG_window_ensembling.md` (written 20:05, before the probe; falsifier,
arms, controls, the widening-K distinction and the expected 0.00 to +0.03 A all fixed there).
Native-free half: job `w_ensemble_clouds` (exit 0, 4351 s wall under a four-job load with two
governor suspensions, peak RSS 0.113 GB; `s26/results/w_selfcopy_ensemble_clouds.json`, complete
126/126; production gate top-75 == `sub` 126/126, cloud max abs 1.8e-15; BLOSUM62 sums asserted
equal to the universe's `sim` and `core.data.top_k` to the pinned pool on every target). Gated
half: job `w_ensemble_endpoint` (exit 0, 10 s, 0.057 GB; `s26/results/w_selfcopy_ensemble_endpoint.json`).
Probe `w_ensemble_probe` (1A13, 35 s, 0.097 GB). Basis on every line: built chain `arm` PRIMARY
(the rebuild basis 3.2126, L57), point cloud carried (3.0483); the selection basis does not exist
for an ensemble. Negative = the arm is better than the shipped emission.

**How fixed-K ensembling differs from widening K, in one sentence each.** Widening K (S17 L12)
draws ONE shortlist from a wider pool, and the extra plausible windows displace near-native
members out of that single top-75 (its ORACLE best 2.104 -> 2.572 A from K = 75 to the full
universe); fixed-K ensembling keeps three shortlists of K = 500, each cut to its own top-75 by
the shipped score (the BLOSUM62 one IS the production set, its ORACLE best 2.306 untouched), and
averages the three CLOUDS in the production frame, so no shortlist is widened and displacement
cannot act; what can act is only whether the three clouds' errors are parallel, and S23 L9 (68%
of the pool's error is common-mode, shared by every member of the universe) and S24 L3 (score
selection makes the bias parallel across sources, cosine 0.943 against a within-source 0.933)
predicted that they are.

**Native-free geometry (medians over 126).** The BLOSUM45 / BLOSUM80 pools overlap the BLOSUM62
pool at 0.83 / 0.89 and their score-selected top-75s overlap the production top-75 at 0.83 /
0.88 (93 distinct members in the union of three); the three clouds sit 0.18 (45 vs 62), 0.14
(80 vs 62) and 0.20 A (45 vs 80) apart; the ensemble cloud moves the built chain 0.185 A from
the production chain in the median (b45 0.33, b80 0.25, ensK 0.32, union3 0.15, boot3 0.25:
the zero-information resample moves it as far as the ensemble does).

**Endpoints (`ST.fmt` verbatim; the PRIMARY contrast and the two controls the falsifier names):**

  ens3 minus shipped [arm]
    a 3.2121 (med 2.9453)   b 3.2126 (med 2.9661)   n=126
    effect -0.0005   median -0.0002   SE 0.0100   MDE 0.0281   effect/MDE -0.02
    iid  CI95 [-0.0206, +0.0199]
    fold CI95 [-0.0194, +0.0162]   folds same sign 3/5   per-fold 0:-0.021 1:-0.012 2:+0.022 3:-0.021 4:+0.023
    64W/62L/0T   worst degradation +0.5570 (1D6X)   p90 +0.0615   power 0.05  Type-M 48.95
    concentration: drop-top10 +0.0179 vs uniform-effect null p10/p50/p90 +0.0062/+0.0171/+0.0296 -> pctile 0.531
    VERDICT: NOT MEASURED (|effect| 0.0005 <= its own MDE 0.0281, 0.02x)
  ens3 minus shipped [cloud]
    a 3.0499 (med 2.7932)   b 3.0483 (med 2.8373)   n=126
    effect +0.0016   median -0.0014   SE 0.0066   MDE 0.0185   effect/MDE +0.09
    iid  CI95 [-0.0095, +0.0158]
    fold CI95 [-0.0086, +0.0161]   folds same sign 2/5   per-fold 0:-0.005 1:+0.028 2:+0.007 3:-0.011 4:-0.008
    68W/58L/0T   worst degradation +0.6286 (1U62)   p90 +0.0475   power 0.06  Type-M 9.87
    concentration: drop-top10 +0.0112 vs uniform-effect null p10/p50/p90 +0.0028/+0.0105/+0.0198 -> pctile 0.540
    VERDICT: NOT MEASURED (|effect| 0.0016 <= its own MDE 0.0185, 0.09x)
  ens3 minus boot3 (the matched zero-information ensembling control) [arm]
    a 3.2121 (med 2.9453)   b 3.2079 (med 3.0424)   n=126
    effect +0.0042   median +0.0004   SE 0.0106   MDE 0.0297   effect/MDE +0.14
    iid  CI95 [-0.0162, +0.0243]
    fold CI95 [-0.0178, +0.0300]   folds same sign 2/5   per-fold 0:-0.028 1:-0.009 2:+0.048 3:-0.011 4:+0.017
    61W/65L/0T   worst degradation +0.5241 (9BAF)   p90 +0.1299   power 0.07  Type-M 5.96
    concentration: drop-top10 +0.0242 vs uniform-effect null p10/p50/p90 +0.0117/+0.0238/+0.0366 -> pctile 0.514
    VERDICT: NOT MEASURED (|effect| 0.0042 <= its own MDE 0.0297, 0.14x)
  ens3 minus boot3 (the matched zero-information ensembling control) [cloud]
    a 3.0499 (med 2.7932)   b 3.0507 (med 2.8659)   n=126
    effect -0.0008   median -0.0025   SE 0.0064   MDE 0.0178   effect/MDE -0.04
    iid  CI95 [-0.0134, +0.0117]
    fold CI95 [-0.0051, +0.0048]   folds same sign 3/5   per-fold 0:-0.006 1:+0.003 2:+0.009 3:-0.005 4:-0.004
    72W/54L/0T   worst degradation +0.2384 (6CEJ)   p90 +0.0724   power 0.05  Type-M 19.64
    concentration: drop-top10 +0.0125 vs uniform-effect null p10/p50/p90 +0.0047/+0.0119/+0.0197 -> pctile 0.541
    VERDICT: NOT MEASURED (|effect| 0.0008 <= its own MDE 0.0178, 0.04x)

**The secondaries, minus shipped (effect / MDE / fold CI / verdict), both bases:**

    arm      b45 -0.0039 / 0.0405 / [-0.0348, +0.0212] NOT MEASURED   b80 +0.0237 / 0.0337 / [+0.0052, +0.0402] NOT MEASURED (0.70x, fold CI above zero, 4/5 folds, Type-M zone: the BLOSUM80 shortlist alone is if anything WORSE)
             ensK -0.0153 / 0.0386 / [-0.0380, +0.0047] NOT MEASURED   union3 -0.0021 / 0.0191 / [-0.0180, +0.0168] NOT MEASURED   boot3 -0.0047 / 0.0273 / [-0.0161, +0.0045] NOT MEASURED
    cloud    b45 +0.0009 / 0.0374 NOT MEASURED   b80 +0.0136 / 0.0250 / [+0.0024, +0.0258] NOT MEASURED (0.54x)   ensK -0.0116 / 0.0223 / [-0.0179, -0.0045] NOT MEASURED (0.52x, 5/5 folds, fold CI below zero: the K-ensemble's CLOUD is 0.012 A nearer, and its built chain -0.015 at 0.40x is not)
             union3 +0.0013 / 0.0113 NOT MEASURED   boot3 +0.0023 / 0.0209 NOT MEASURED

**Verdict.** H_E is refuted at this instrument: the ensemble neither beats the shipped chain
(-0.0005, 0.02x MDE) nor its zero-information control (+0.0042, 0.14x). Power: the design
excludes a gain of 0.028 A or more on the built chain (the MDE of ens3 minus shipped) and of
0.019 A or more on the cloud; the registered expectation (0.00 to +0.03) held. The mechanism is
the registered one: the three shortlists overlap at 0.83 to 0.88, their clouds sit 0.14 to 0.20 A
apart, and averaging them moves the emission by exactly as much as a resample of the production
set does (0.185 vs 0.249 A median chain move) with the same effect on accuracy (none): parallel
errors average to themselves. The only cell with the fold CI excluding zero in the helpful
direction is the K-variant ensemble on the CLOUD (-0.0116, 0.52x MDE, 5/5 folds), which its own
built chain does not carry (-0.0153, 0.40x); it is quoted as suggestive and nothing more, and it
is the S17 L12 direction (the consensus readout improves with K on the point cloud, and the
built chain does not follow). The mandatory test-time-ensembling direction is closed in the
fixed-K form with a power statement; the widening-K form was closed by S17 L12. Replication:
`boot3` is seeded (`s15.seed.stable_rng`); a second seed and the reversed fold order are the
pre-declared replication for a positive, and there is none to replicate; the ensemble itself
has no random element. Deviations from the PREREG: none; the `sel` basis is undefined for an
ensemble and is not quoted. Artefacts: `s26/PREREG_window_ensembling.md`, `s26/w_ensemble.py`,
`s26/w_ensemble_test.py` (ALL OK), `s26/results/w_selfcopy_ensemble_probe_1A13.json`,
`w_selfcopy_ensemble_clouds.json`, `w_selfcopy_ensemble_endpoint.json`,
`s26/jobs_done/w_ensemble_{probe,clouds,endpoint}.json`.

## L85 -- TOURNAMENT ITEM 9, window_provenance (W), CENSUS AND ORACLE CONTRAST: 73% OF EVERY POOL AND OF EVERY TOP-75 IS FRAGMENT WINDOWS, 0.6% WHOLE PEPTIDES; IN THE POOL A PEPTIDE-DERIVED WINDOW IS 0.42 A NEARER THE NATIVE THAN A FRAGMENT WINDOW (3.2x MDE, 5/5 FOLDS); INSIDE THE SCORE-SELECTED TOP-75 THE GAP SHRINKS TO 0.09 (0.96x MDE) AND WHOLE+TERMINAL VS FRAGMENT TO 0.13 (1.13x, TYPE-M ZONE); THE READOUT TEST H_P3 RUNS NEXT (2026-09-13, W)

Pre-registered in `s26/PREREG_window_provenance.md` (written before the code ran on any target).
Census (native-free, codes and sequences only): job `w_provenance_census` (exit 0, 5 s, peak
0.004 GB; `s26/results/w_selfcopy_provenance_census.json`, complete 126/126, 0 unresolved
windows). ORACLE contrast (single-window basis, the universe's `rr`): job `w_provenance_oracle`
(exit 0, 10 s, 0.055 GB; `s26/results/w_selfcopy_provenance_oracle.json`). Class rule: a pool
window's string looked up among the length-n substrings of the 787 peptides (whole: parent length
n; terminal: offset 0 or the end; interior) and the 6,003 fragments (fragment), the FIRST parent
in database order winning, the universe's `org` flag arbitrating a string present in both banks;
synthetic test `s26/w_provenance_test.py` ALL OK. Multi-parent strings are common (15,990 at
length 9 down to 1,224 at 16, mostly overlapping fragments of one protein) and are counted, not
guessed.

**Census (H_P1), means over 126 targets:**

    class       K = 500 pool   shipped top-75   targets with one in the top-75
    whole          0.6%            0.7%           30 / 126
    terminal       8.0%            7.2%          113 / 126
    interior      18.7%           18.3%
    fragment      72.7%           73.8%
    unresolved     0.0%            0.0%

The score keeps the pool's class mix almost unchanged: it neither favours nor removes peptide
windows. The registered "about 1% whole peptides" holds (0.6 to 0.7%).

**ORACLE contrast (H_P2), per-target mean `rr` of one class minus another, `ST.fmt` verbatim for
the two contrasts the falsifier names (single-window basis on both sides; negative = the first
class is nearer the native):**

  pool: class peptide_any minus class fragment, per-target mean ORACLE rr (single-window basis), n=126 targets with both
    a 4.1721 (med 4.0394)   b 4.5882 (med 4.3129)   n=126
    effect -0.4161   median -0.4358   SE 0.0468   MDE 0.1310   effect/MDE -3.18
    iid  CI95 [-0.5066, -0.3220]
    fold CI95 [-0.4472, -0.3758]   folds same sign 5/5   per-fold 0:-0.341 1:-0.401 2:-0.463 3:-0.439 4:-0.434
    101W/25L/0T   worst degradation +1.1974 (7VI4)   p90 +0.2005   power 1.00  Type-M 1.00
    concentration: drop-top10 -0.3272 vs uniform-effect null p10/p50/p90 -0.3902/-0.3278/-0.2696 -> pctile 0.504
    VERDICT: BETTER
  top75: class whole_or_terminal minus class fragment, per-target mean ORACLE rr (single-window basis), n=114 targets with both
    a 3.5381 (med 3.4172)   b 3.6695 (med 3.7564)   n=114
    effect -0.1313   median -0.0533   SE 0.0414   MDE 0.1159   effect/MDE -1.13
    iid  CI95 [-0.2139, -0.0533]
    fold CI95 [-0.1646, -0.1096]   folds same sign 5/5   per-fold 0:-0.106 1:-0.193 2:-0.123 3:-0.143 4:-0.107
    66W/48L/0T   worst degradation +0.8160 (5MXS)   p90 +0.2992   power 0.89  Type-M 1.07
    concentration: drop-top10 -0.0286 vs uniform-effect null p10/p50/p90 -0.0789/-0.0306/+0.0136 -> pctile 0.520
    VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.07x]

The other contrasts (effect / MDE / fold CI / verdict): pool whole+terminal minus fragment -0.421
/ 0.136 / [-0.453, -0.388] BETTER (3.1x, 5/5); pool interior minus fragment -0.411 / 0.136 /
[-0.447, -0.368] BETTER (3.0x, 5/5); pool whole+terminal minus interior -0.011 / 0.085 NOT
MEASURED (0.13x); top-75 peptide_any minus fragment -0.086 / 0.090 / [-0.116, -0.056] NOT MEASURED
(0.96x, 5/5 folds); top-75 interior minus fragment -0.061 / 0.093 NOT MEASURED (0.66x); top-75
whole+terminal minus interior -0.056 / 0.114 NOT MEASURED (0.49x).

Reading. In the retrieval pool, ANY peptide-derived window (whole, terminal or interior alike) is
0.41 to 0.42 A nearer the native than a fragment window, 5/5 folds, 3x its MDE: S19's "the peptide
corpus carries the sequence-structure channel; fragments are training-only material whose
conformation is held by contacts outside the window" (S24 L8), measured for the first time on the
single-window basis of the shipped pool. The position of the window in its parent peptide does
not matter (whole+terminal minus interior 0.13x MDE). Inside the score-selected top-75 the gap
shrinks to -0.09 (0.96x) for any peptide window and -0.13 (1.13x, Type-M zone) for whole+terminal
windows: the shipped score already removes most of the class difference, and what remains is at
the edge of what n = 114 to 126 can see. The registered expectation for H_P2 (NOT MEASURED inside
the top-75, well under 0.1 A) held for the peptide_any and interior contrasts and was exceeded, in
the Type-M zone, by the whole+terminal one. By the PREREG's rule the ACHIEVABLE readout test H_P3
(fragment-class relative weight in {0, 0.5, 1, 2}, chosen leave-fold-out on the built chain,
against the uniform readout and the permuted-weight control) therefore runs; the registered
expectation for it is 0.00 to -0.02 A on the built chain against an MDE of about 0.05
(NOT MEASURED), because a 0.13 A single-window difference on 8% of the members of a 75-member
average is worth at most a few hundredths through the mean, and the L44/L64 controls show that a
few members' change moves the chain mostly orthogonally to the native. Nothing here is deployable
yet: the class of a window is native-free, the contrast is ORACLE. Replication: deterministic
(the census and the `rr` are fixed); the bootstrap is seeded.

## L86 -- STERIC REJECT, BUILT CHAIN: THE POINT-CLOUD RESULT REPLICATES ACROSS BASES. R@1e4 IS +0.248 A WORSE THAN THE RE-PROJECTED SHIPPED TOP-75 (fold CI [+0.12, +0.35], 4/5 FOLDS, TYPE-M) AND +0.157 WORSE THAN THE SAME COUNT REJECTED AT RANDOM; DOSE MONOTONE; THE FILTER IS CLOSED ON BOTH BASES (2026-09-13, PH)

`s26/ph_reject.py chain` then `report`, `s26/results/ph_reject_chain.json` (complete 126/126),
`s26/results/ph_reject_report.json` (both bases), job `s26/jobs_done/ph_reject_chain.json` (exit
0, 11,976 s = 3.3 h, 22.8 projections per target, peak RSS 0.09 GB). Pre-registered in
`s26/PREREG_amber_reject.md` sections 3 to 5, addendum 1 (fallback and moved subset) and
addendum 2 (the re-projected anchor). This is the falsifier's last leg and the cross-basis
replication of L43.

BASIS: BUILT CHAIN on both sides (`rmsd_arm`): every retained set is coordinate-averaged and
then projected through `s12.instrument.project` (ramah at 0.3, multi-start, exact gradient,
the production call) and the chain is scored ORACLE against `nat_ca`. The ANCHOR is the
shipped top-75 RE-PROJECTED through the same call, not the stored production chain, so all arms
share one instrument. Its deviation from the stored chain, reported as registered: mean
3.2126 A against the production 3.2148 (the "leaderboard rebuild 3.2126" of the
state brief and L57: the production path round-trips the cloud through float32 before
projecting, `Config.reference_precision`, and this instrument does not); per-target RMSD
difference mean -0.0021, max |0.171|; CA deviation between the two chains mean
0.105 A, above 0.05 A on 46 targets and above 0.5 A on 6
(2LWU 1.25, 7JS6 1.54, 6QAX 1.56, 2BP4 1.62), which is the projection's own branch
degeneracy (`core/project.py` docstring: a 1e-13 input difference can route L-BFGS-B into the
other torsion branch, up to 1.6 A). Arms, controls and fallbacks as in L43 (R = reject e_amber
> T and refill to 75, the PRIMARY reading; S = reject, no refill; RANDR/RANDS same count at
random; PERMR/PERMS the same threshold on a permuted energy); the built-chain controls at 1e4
carry 4 draws each (registered in section 3), the other thresholds carry R and S only.
Eleven `ST.fmt` blocks verbatim:

      built_chain [all, n=126] R@1e4 minus anchor (negative = arm better)
        a 3.4609 (med 3.2449)   b 3.2126 (med 2.9661)   n=126
        effect +0.2483   median +0.0016   SE 0.0758   MDE 0.2123   effect/MDE +1.17
        iid  CI95 [+0.1021, +0.4044]
        fold CI95 [+0.1197, +0.3479]   folds same sign 4/5   per-fold 0:+0.267 1:-0.002 2:+0.226 3:+0.421 4:+0.310
        49W/67L/10T   worst degradation +4.2043 (8T61)   p90 +1.3949   power 0.91  Type-M 1.06
        concentration: drop-top10 +0.3567 vs uniform-effect null p10/p50/p90 +0.2615/+0.3545/+0.4565 -> pctile 0.511
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      built_chain [all, n=126] R@1e4 minus RANDR@1e4 (matched control)
        a 3.4609 (med 3.2449)   b 3.3036 (med 3.1840)   n=126
        effect +0.1574   median +0.0085   SE 0.0488   MDE 0.1367   effect/MDE +1.15
        iid  CI95 [+0.0701, +0.2539]
        fold CI95 [+0.0480, +0.2811]   folds same sign 4/5   per-fold 0:+0.120 1:-0.049 2:+0.182 3:+0.376 4:+0.160
        56W/68L/2T   worst degradation +3.2399 (8T61)   p90 +0.8244   power 0.90  Type-M 1.06
        concentration: drop-top10 +0.2204 vs uniform-effect null p10/p50/p90 +0.1582/+0.2166/+0.2815 -> pctile 0.530
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      built_chain [all, n=126] R@1e4 minus PERMR@1e4 (matched control)
        a 3.4609 (med 3.2449)   b 3.3589 (med 3.1940)   n=126
        effect +0.1021   median +0.0000   SE 0.0414   MDE 0.1159   effect/MDE +0.88
        iid  CI95 [+0.0241, +0.1886]
        fold CI95 [-0.0080, +0.2034]   folds same sign 4/5   per-fold 0:+0.017 1:-0.113 2:+0.137 3:+0.277 4:+0.174
        57W/61L/8T   worst degradation +2.6686 (9KAR)   p90 +0.6660   power 0.69  Type-M 1.20
        concentration: drop-top10 +0.1607 vs uniform-effect null p10/p50/p90 +0.1091/+0.1587/+0.2107 -> pctile 0.518
        VERDICT: NOT MEASURED (|effect| 0.1021 <= its own MDE 0.1159, 0.88x)
      built_chain [all, n=126] S@1e4 minus anchor (negative = arm better)
        a 3.3164 (med 3.1608)   b 3.2126 (med 2.9661)   n=126
        effect +0.1037   median +0.0025   SE 0.0333   MDE 0.0933   effect/MDE +1.11
        iid  CI95 [+0.0414, +0.1747]
        fold CI95 [+0.0373, +0.1662]   folds same sign 5/5   per-fold 0:+0.014 1:+0.002 2:+0.172 3:+0.194 4:+0.129
        47W/66L/13T   worst degradation +1.7906 (9BAF)   p90 +0.4638   power 0.88  Type-M 1.08
        concentration: drop-top10 +0.1422 vs uniform-effect null p10/p50/p90 +0.1000/+0.1424/+0.1914 -> pctile 0.498
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.08x]
      built_chain [all, n=126] S@1e4 minus RANDS@1e4 (matched control)
        a 3.3164 (med 3.1608)   b 3.2612 (med 3.0604)   n=126
        effect +0.0551   median +0.0028   SE 0.0269   MDE 0.0754   effect/MDE +0.73
        iid  CI95 [+0.0041, +0.1101]
        fold CI95 [-0.0005, +0.1064]   folds same sign 4/5   per-fold 0:-0.048 1:+0.027 2:+0.118 3:+0.127 4:+0.055
        45W/68L/13T   worst degradation +1.3396 (8T61)   p90 +0.3322   power 0.53  Type-M 1.36
        concentration: drop-top10 +0.0941 vs uniform-effect null p10/p50/p90 +0.0586/+0.0924/+0.1297 -> pctile 0.524
        VERDICT: NOT MEASURED (|effect| 0.0551 <= its own MDE 0.0754, 0.73x)
      built_chain [all, n=126] RANDR@1e4 minus anchor (negative = arm better)
        a 3.3036 (med 3.1840)   b 3.2126 (med 2.9661)   n=126
        effect +0.0909   median +0.0026   SE 0.0432   MDE 0.1211   effect/MDE +0.75
        iid  CI95 [+0.0110, +0.1787]
        fold CI95 [+0.0456, +0.1322]   folds same sign 5/5   per-fold 0:+0.148 1:+0.047 2:+0.045 3:+0.046 4:+0.151
        59W/65L/2T   worst degradation +2.1588 (8FLP)   p90 +0.7406   power 0.56  Type-M 1.34
        concentration: drop-top10 +0.1645 vs uniform-effect null p10/p50/p90 +0.1103/+0.1633/+0.2176 -> pctile 0.513
        VERDICT: NOT MEASURED (|effect| 0.0909 <= its own MDE 0.1211, 0.75x)
      built_chain [all, n=126] RANDS@1e4 minus anchor (negative = arm better)
        a 3.2612 (med 3.0604)   b 3.2126 (med 2.9661)   n=126
        effect +0.0486   median +0.0002   SE 0.0175   MDE 0.0490   effect/MDE +0.99
        iid  CI95 [+0.0167, +0.0836]
        fold CI95 [+0.0121, +0.0698]   folds same sign 4/5   per-fold 0:+0.063 1:-0.025 2:+0.055 3:+0.067 4:+0.074
        48W/65L/13T   worst degradation +1.1826 (9BAF)   p90 +0.2279   power 0.79  Type-M 1.13
        concentration: drop-top10 +0.0695 vs uniform-effect null p10/p50/p90 +0.0477/+0.0684/+0.0920 -> pctile 0.528
        VERDICT: NOT MEASURED (|effect| 0.0486 <= its own MDE 0.0490, 0.99x)
      built_chain [all, n=126] R@1e3 minus anchor (negative = arm better)
        a 3.7757 (med 3.5036)   b 3.2126 (med 2.9661)   n=126
        effect +0.5631   median +0.0248   SE 0.1138   MDE 0.3187   effect/MDE +1.77
        iid  CI95 [+0.3549, +0.7891]
        fold CI95 [+0.4413, +0.7037]   folds same sign 5/5   per-fold 0:+0.541 1:+0.346 2:+0.530 3:+0.822 4:+0.578
        36W/69L/21T   worst degradation +4.6340 (2NDN)   p90 +2.7185   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.6994 vs uniform-effect null p10/p50/p90 +0.5482/+0.6942/+0.8531 -> pctile 0.518
        VERDICT: WORSE
      built_chain [all, n=126] S@1e3 minus anchor (negative = arm better)
        a 3.4048 (med 3.2359)   b 3.2126 (med 2.9661)   n=126
        effect +0.1921   median +0.0007   SE 0.0477   MDE 0.1338   effect/MDE +1.44
        iid  CI95 [+0.1048, +0.2907]
        fold CI95 [+0.1032, +0.2469]   folds same sign 5/5   per-fold 0:+0.231 1:+0.021 2:+0.273 3:+0.212 4:+0.208
        37W/63L/26T   worst degradation +2.4536 (9BAF)   p90 +0.7604   power 0.98  Type-M 1.01
        concentration: drop-top10 +0.2423 vs uniform-effect null p10/p50/p90 +0.1781/+0.2389/+0.3077 -> pctile 0.529
        VERDICT: WORSE
      built_chain [moved, n=116] R@1e4 minus anchor (negative = arm better)
        a 3.4869 (med 3.2751)   b 3.2172 (med 3.0247)   n=116
        effect +0.2697   median +0.0077   SE 0.0820   MDE 0.2298   effect/MDE +1.17
        iid  CI95 [+0.1144, +0.4415]
        fold CI95 [+0.1290, +0.3943]   folds same sign 4/5   per-fold 0:+0.267 1:-0.002 2:+0.283 3:+0.484 4:+0.321
        49W/67L/0T   worst degradation +4.2043 (8T61)   p90 +1.5485   power 0.91  Type-M 1.06
        concentration: drop-top10 +0.3904 vs uniform-effect null p10/p50/p90 +0.2849/+0.3857/+0.4945 -> pctile 0.520
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      built_chain [moved, n=116] R@1e4 minus RANDR@1e4 (matched control)
        a 3.4869 (med 3.2751)   b 3.3178 (med 3.2429)   n=116
        effect +0.1692   median +0.0099   SE 0.0522   MDE 0.1462   effect/MDE +1.16
        iid  CI95 [+0.0731, +0.2704]
        fold CI95 [+0.0519, +0.3172]   folds same sign 4/5   per-fold 0:+0.120 1:-0.043 2:+0.202 3:+0.438 4:+0.165
        51W/65L/0T   worst degradation +3.2399 (8T61)   p90 +0.8244   power 0.90  Type-M 1.06
        concentration: drop-top10 +0.2393 vs uniform-effect null p10/p50/p90 +0.1716/+0.2352/+0.3086 -> pctile 0.533
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]

All four thresholds, `[all, n=126]`, arm minus anchor, built chain:

    1e3  R   +0.5631 (median  +0.0248, +1.77x MDE)  fold [+0.4413, +0.7037]   5/5  36W/69L/21T  WORSE
    1e3  S   +0.1921 (median  +0.0007, +1.44x MDE)  fold [+0.1032, +0.2469]   5/5  37W/63L/26T  WORSE
    1e4  R   +0.2483 (median  +0.0016, +1.17x MDE)  fold [+0.1197, +0.3479]   4/5  49W/67L/10T  WORSE
    1e4  S   +0.1037 (median  +0.0025, +1.11x MDE)  fold [+0.0373, +0.1662]   5/5  47W/66L/13T  WORSE
    1e5  R   +0.1190 (median  +0.0000, +0.78x MDE)  fold [+0.0482, +0.1844]   5/5  57W/62L/7T   NOT MEASURED
    1e5  S   +0.0833 (median  +0.0000, +0.95x MDE)  fold [+0.0020, +0.1560]   3/5  58W/59L/9T   NOT MEASURED
    1e6  R   +0.0915 (median  +0.0000, +0.65x MDE)  fold [+0.0328, +0.1377]   4/5  57W/61L/8T   NOT MEASURED
    1e6  S   +0.0511 (median  +0.0000, +0.76x MDE)  fold [-0.0016, +0.0983]   3/5  61W/55L/10T  NOT MEASURED

Threshold sweep (order statistic): R oracle minimum -0.3442, split-half transfer -0.1473 (43%%),
k_eff 3.86; S -0.1345, transfer -0.0535 (40%%): the mildest rung transfers, not a gain (R@1e6
+0.092 at 0.65x MDE, S@1e6 +0.051 at 0.76x, both NOT MEASURED).

READING. (1) The built chain REPLICATES the point cloud in sign on the primary reading: R@1e4
is +0.248 A worse than the re-projected shipped top-75 (fold CI [+0.120, +0.348], 49W/67L/10T)
and +0.157 worse than rejecting the same count at random with a judgment-free refill (fold CI
[+0.048, +0.281]); L43 had +0.228 and +0.167 on the point cloud. Both are Type-M-zone
magnitudes (1.17x and 1.15x MDE), as on the point cloud, and both are tail-carried (median
+0.002 and +0.009 against means +0.248 and +0.157; p90 +1.39; worst +4.20 on 8T61), so the
Adversary's L54 caveats 1 and 2 transfer unchanged to this basis: near zero on the median
target, catastrophic on the minority whose pool has no survivor or whose refill reaches deep.
One difference from L43 that is reported rather than smoothed: fold 1 has the opposite sign on
R (-0.002 against +0.267 / +0.226 / +0.421 / +0.310 on the other four), so R@1e4 is 4/5 folds
on the built chain where it was 5/5 on the point cloud; the fold CI still excludes zero. S@1e4
is +0.104 worse (5/5 folds, 1.11x MDE, Type-M) against +0.108 on the point cloud. (2) The
dose is again monotone: 1e3 +0.563 (1.77x MDE, 5/5 folds, MEASURED), 1e4 +0.248, 1e5 +0.119
(0.78x), 1e6 +0.092 (0.65x); the limit is the anchor. (3) The built chain is noisier than the
point cloud by construction (the projection adds its own branch noise on every arm), so the
matched-control contrasts that were Type-M on the point cloud are UNDERPOWERED here: S vs
RANDS +0.055 at 0.73x MDE (fold CI [-0.0005, +0.106]), R vs PERMR +0.102 at 0.88x, S vs PERMS
+0.050 at 0.59x; RANDR vs anchor +0.091 at 0.75x, RANDS vs anchor +0.049 at 0.99x. The L54
caveat 3 stands: the refill-cost / choice-cost split is a point estimate on both bases. (4)
The moved subset (secondary, n = 116 / 113) agrees: R +0.270 WORSE (Type-M), R vs RANDR +0.169
WORSE (Type-M), S +0.116 WORSE (Type-M).

POWER. Built-chain SEs 0.03 to 0.11 A, MDEs 0.09 to 0.32; at 1e5 and 1e6 R sits at 0.78x and
0.65x its MDE (SE 0.055 and 0.050), so "the mildest threshold is null" is UNDERPOWERED for a
gain below about 0.14 A on this basis and MEASURED against any harm above it. Nothing is
positive; the replication of the harmful direction is this entry (cross-basis, same sign, same
monotone dose, same tail structure).

DISPOSITION. The falsifier of `s26/PREREG_amber_reject.md` section 4 required R or S at 1e4 to
BEAT the anchor and its matched control on the built chain; both are WORSE than the anchor
with fold CIs excluding zero (R 4/5 folds, S 5/5), and R is worse than its matched control.
AMBER as a steric reject filter at a physical threshold, with or without refill, is CLOSED on
both bases in the two forms the record had not measured (physics-set count, refill), with a
harmful sign at 1e3 and 1e4 and an underpowered null at 1e5 and 1e6. The functional lever's
filter form (s19 Q3, s24 D1-C) keeps its closure and gains a sixth and seventh instrument. Why,
from L23: the members the threshold condemns are condemned by the builder's side-chain
placement (96.8%), not by their backbones; removing them removes compact members whose error
cancelled in the average (S23 L5). The remaining physics question on this pool is
`IDEA_rotamer_relief` part B (does the singularity move when the side chains are relieved),
which runs next under the tournament.

---

## L87 -- C3 STAGE 1 REPLICATION (new seeds, reversed order): EVERY CONTRAST LANDS INSIDE L39's FOLD CI; THE TOWARD-MEMBER CONTROL's -0.018 A IMPROVEMENT IS NO LONGER PROVISIONAL (-0.0207 [-0.0250, -0.0166], 5/5) (2026-09-13, PH)

`s26/ph_c3.py stage1 --rep`, `s26/results/ph_c3_stage1_rep.json` (complete 126/126), job
`s26/jobs_done/ph_c3_stage1_rep.json` (exit 0, 10 s, peak RSS 0.04 GB). Registered in
`s26/PREREG_c3_control.md` addendum 2 before it ran: new draw seeds (salt "rep"), targets
processed in reversed pinned order, fold labels untouched. Basis as in L39: arm input the
BUILT CHAIN (`rmsd_arm` 3.2148), arm output the RELAXED CHAIN (`rmsd_full` 3.2355); the
controls displace the built chain's CA trace by AMBER's own per-atom RMS magnitude. The arm
has no randomness, so only the controls' draws change. `ST.fmt` verbatim, the five contrasts
the replication was registered to judge:

      stage1-REP TOWARD-MEMBER minus do-nothing
        a 3.1941 (med 2.9650)   b 3.2148 (med 2.9661)   n=126
        effect -0.0207   median -0.0158   SE 0.0042   MDE 0.0117   effect/MDE -1.77
        iid  CI95 [-0.0287, -0.0126]
        fold CI95 [-0.0250, -0.0166]   folds same sign 5/5   per-fold 0:-0.026 1:-0.020 2:-0.026 3:-0.015 4:-0.016
        94W/32L/0T   worst degradation +0.1064 (8TXS)   p90 +0.0279   power 1.00  Type-M 1.00
        concentration: drop-top10 -0.0121 vs uniform-effect null p10/p50/p90 -0.0173/-0.0123/-0.0073 -> pctile 0.522
        VERDICT: BETTER
      stage1-REP AMBER minus matched-magnitude RANDOM (16 draws)
        a 3.2355 (med 2.9757)   b 3.2255 (med 2.9698)   n=126
        effect +0.0100   median +0.0094   SE 0.0035   MDE 0.0098   effect/MDE +1.02
        iid  CI95 [+0.0033, +0.0172]
        fold CI95 [+0.0059, +0.0175]   folds same sign 5/5   per-fold 0:+0.008 1:+0.026 2:+0.006 3:+0.005 4:+0.007
        49W/77L/0T   worst degradation +0.1442 (2LM8)   p90 +0.0561   power 0.81  Type-M 1.12
        concentration: drop-top10 +0.0156 vs uniform-effect null p10/p50/p90 +0.0109/+0.0156/+0.0204 -> pctile 0.509
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.12x]
      stage1-REP AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)
        a 3.2355 (med 2.9757)   b 3.1941 (med 2.9650)   n=126
        effect +0.0414   median +0.0387   SE 0.0059   MDE 0.0165   effect/MDE +2.50
        iid  CI95 [+0.0299, +0.0527]
        fold CI95 [+0.0331, +0.0501]   folds same sign 5/5   per-fold 0:+0.047 1:+0.058 2:+0.043 3:+0.028 4:+0.033
        30W/96L/0T   worst degradation +0.3600 (2BP4)   p90 +0.1093   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0512 vs uniform-effect null p10/p50/p90 +0.0439/+0.0513/+0.0588 -> pctile 0.492
        VERDICT: WORSE
      stage1-REP RANDOM minus do-nothing
        a 3.2255 (med 2.9698)   b 3.2148 (med 2.9661)   n=126
        effect +0.0107   median +0.0083   SE 0.0014   MDE 0.0039   effect/MDE +2.79
        iid  CI95 [+0.0082, +0.0135]
        fold CI95 [+0.0093, +0.0120]   folds same sign 5/5   per-fold 0:+0.013 1:+0.012 2:+0.011 3:+0.008 4:+0.010
        28W/98L/0T   worst degradation +0.1035 (1S9Z)   p90 +0.0283   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0124 vs uniform-effect null p10/p50/p90 +0.0106/+0.0124/+0.0143 -> pctile 0.513
        VERDICT: WORSE
      stage1-REP ORACLE cos: AMBER minus RANDOM
        a -0.0491 (med -0.0498)   b -0.0032 (med -0.0030)   n=126
        effect -0.0459   median -0.0401   SE 0.0155   MDE 0.0435   effect/MDE -1.05
        iid  CI95 [-0.0770, -0.0157]
        fold CI95 [-0.0734, -0.0212]   folds same sign 5/5   per-fold 0:-0.042 1:-0.101 2:-0.050 3:-0.003 4:-0.036
        77W/49L/0T   worst degradation +0.4596 (5MXS)   p90 +0.1598   power 0.84  Type-M 1.10
        concentration: drop-top10 -0.0172 vs uniform-effect null p10/p50/p90 -0.0377/-0.0175/+0.0020 -> pctile 0.508
        VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.10x]

    first run (L39)                          replication              first run's fold CI       inside?
    toward-member minus do-nothing  -0.0178  -0.0207 [-0.0250, -0.0166]  [-0.0255, -0.0124]     yes
    AMBER minus random              +0.0111  +0.0100 [+0.0059, +0.0175]  [+0.0062, +0.0171]     yes
    AMBER minus toward-member       +0.0385  +0.0414 [+0.0331, +0.0501]  [+0.0298, +0.0479]     yes
    random minus do-nothing         +0.0096  +0.0107 [+0.0093, +0.0120]  [+0.0077, +0.0122]     yes
    ORACLE cos AMBER minus random   -0.0503  -0.0459 [-0.0734, -0.0212]  [-0.0735, -0.0253]     yes

Every replicated quantity lands inside the first run's fold-clustered CI, all 5/5 folds. The
PROVISIONAL label on the toward-member line of L39 is lifted: a zero-information move of the
same size as the relaxation's, toward a random member of the shipped pool, improves the built
chain by 0.018 to 0.021 A (94W/32L in the replication, 1.77x its MDE), where the relaxation's
own move worsens it by 0.021. The reading in L39 stands unchanged: this is a property of the
projection's displacement (S16 L27: the projection moved the chain 0.166 A away from the point
cloud with a negative cosine, and the pool members surround the cloud), not of physics, and it
is not a proposal. The Type-M status of AMBER-minus-random (1.02x MDE in the replication) is
unchanged and the L46 caveat still applies to any quotation of that number.

---

## L88 -- BRANCH SELECT (tournament 5): THE RELAXED ENERGY PICKS THE PROJECTION BRANCH AS WELL AS THE OBJECTIVE AND NO BETTER (+0.0055 vs production, 0.18x MDE, 42 ties; -0.102 vs a random branch, 2.29x, 5/5); THE RAW SINGLE POINT PICKS WORSE (+0.092); CLOSED AS AN ACCURACY STEP (2026-09-13, PH)

`s26/ph_branch.py solutions | relax | report`, artefacts `s26/results/ph_branch_solutions.json`
(126/126, G1 exact 0.0 on every target), `ph_branch_relax.json` (126/126; jobs `ph_branch_relax`
41 cells + `ph_branch_relax2` 85 cells, AMBER, peak RSS 0.27 GB probe / 0.3 GB run, about 65 s
per target), `ph_branch_report.json`. Pre-registered in `s26/PREREG_branch_select.md` (sections
1 to 7, addendum 1 with the gates); tournament item 5 (L51). Per-target cells under
`s26/results/ph_branch_{sol_,}cells/`.

THE OPERATOR, native-free. The production projection chooses among FIVE solutions at its
lam = 0.3 rung (the warm start from lam = 0 and the four generic starts alpha / beta / PPII /
extended) by the lowest objective; the choice was reconstructed and equals `I.project` to
0.0 A on 126/126. The arm relaxes all five built chains with the production operator
(`refine_coords(k=10, steps=0)`), reads the converged energy with the restraint off (`e1`) and
emits the BUILT chain of the lowest converged `e1` (ties averaged by `ST.argmin_tied`; 0 ties
occurred; on 10 targets fewer than five converged and the pick is among the converged). All 126
targets have at least two distinct solutions (4.77 on average), so the effective n is 126 and
the moved subset equals the whole. BASIS: BUILT CHAIN (`rmsd_arm`) for the primary; the relaxed
chain of the same pick as a secondary. Controls in the operator's space: the production choice
(anchor, the re-projected 3.2126); the exact expectation of a uniformly random pick among the
five; the raw single point `e0` as the declared secondary picker. `ST.fmt` verbatim:

      [all, n=126] pick_e1 minus production (built chain)
        a 3.2181 (med 2.9950)   b 3.2126 (med 2.9661)   n=126
        effect +0.0055   median +0.0000   SE 0.0109   MDE 0.0306   effect/MDE +0.18
        iid  CI95 [-0.0154, +0.0271]
        fold CI95 [+0.0004, +0.0106]   folds same sign 4/5   per-fold 0:+0.008 1:+0.012 2:+0.010 3:-0.003 4:+0.001
        39W/45L/42T   worst degradation +0.6492 (1D6X)   p90 +0.1157   power 0.08  Type-M 4.77
        concentration: drop-top10 +0.0269 vs uniform-effect null p10/p50/p90 +0.0142/+0.0262/+0.0389 -> pctile 0.534
        VERDICT: NOT MEASURED (|effect| 0.0055 <= its own MDE 0.0306, 0.18x)
      [all, n=126] pick_e1 minus random pick (built chain)
        a 3.2181 (med 2.9950)   b 3.3202 (med 3.0808)   n=126
        effect -0.1021   median -0.0752   SE 0.0159   MDE 0.0446   effect/MDE -2.29
        iid  CI95 [-0.1339, -0.0718]
        fold CI95 [-0.1270, -0.0750]   folds same sign 5/5   per-fold 0:-0.064 1:-0.099 2:-0.119 3:-0.074 4:-0.142
        91W/35L/0T   worst degradation +0.2983 (2BP4)   p90 +0.0937   power 1.00  Type-M 1.00
        concentration: drop-top10 -0.0662 vs uniform-effect null p10/p50/p90 -0.0868/-0.0661/-0.0483 -> pctile 0.497
        VERDICT: BETTER
      [all, n=126] random pick minus production (built chain)
        a 3.3202 (med 3.0808)   b 3.2126 (med 2.9661)   n=126
        effect +0.1076   median +0.0641   SE 0.0145   MDE 0.0406   effect/MDE +2.65
        iid  CI95 [+0.0795, +0.1365]
        fold CI95 [+0.0795, +0.1327]   folds same sign 5/5   per-fold 0:+0.073 1:+0.111 2:+0.130 3:+0.071 4:+0.143
        32W/94L/0T   worst degradation +0.6698 (2RUO)   p90 +0.3744   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.1267 vs uniform-effect null p10/p50/p90 +0.1074/+0.1266/+0.1460 -> pctile 0.501
        VERDICT: WORSE
      [all, n=126] pick_e0 minus production (built chain)
        a 3.3048 (med 3.0309)   b 3.2126 (med 2.9661)   n=126
        effect +0.0922   median +0.0000   SE 0.0300   MDE 0.0840   effect/MDE +1.10
        iid  CI95 [+0.0387, +0.1557]
        fold CI95 [+0.0270, +0.1676]   folds same sign 4/5   per-fold 0:-0.004 1:+0.108 2:+0.076 3:+0.034 4:+0.218
        39W/58L/29T   worst degradation +2.2014 (2F3A)   p90 +0.3328   power 0.87  Type-M 1.08
        concentration: drop-top10 +0.1219 vs uniform-effect null p10/p50/p90 +0.0823/+0.1208/+0.1623 -> pctile 0.516
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.08x]
      [all, n=126] pick_e1 minus production (relaxed chain)
        a 3.2400 (med 3.0159)   b 3.2437 (med 2.9873)   n=126
        effect -0.0036   median +0.0000   SE 0.0141   MDE 0.0395   effect/MDE -0.09
        iid  CI95 [-0.0342, +0.0220]
        fold CI95 [-0.0212, +0.0089]   folds same sign 2/5   per-fold 0:+0.014 1:+0.005 2:-0.038 3:+0.000 4:-0.000
        35W/49L/42T   worst degradation +0.6488 (1D6X)   p90 +0.0842   power 0.06  Type-M 9.18
        concentration: drop-top10 +0.0264 vs uniform-effect null p10/p50/p90 +0.0135/+0.0255/+0.0385 -> pctile 0.538
        VERDICT: NOT MEASURED (|effect| 0.0036 <= its own MDE 0.0395, 0.09x)

    means, built chain:  production 3.2126   e1 pick 3.2181   e0 pick 3.3048   random pick 3.3202   ORACLE min over five 3.1301
    e1 pick equals the production choice on 42 of 126 targets (the 42 exact ties in the first block)
    ORACLE DIAGNOSTIC: the e1 pick is the per-target best on 32.5% of targets, the production choice on 29.4%
    ORACLE min over the five (an order statistic): -0.1901; valid across-target null -0.2475 (share 1.30);
      split-half transfer -0.1067 (56%); k_eff 4.63; argmin counts over the five columns [33, 24, 39, 18, 12]

READING. (1) The falsifier (`PREREG_branch_select.md` section 4) required the relaxed-energy
pick to BEAT the production choice past its MDE with the fold CI excluding zero. It does not:
+0.0055 A, 0.18x MDE, 39W/45L/42T, the fold CI [+0.0004, +0.0106] on the wrong side of zero
anyway. NOT MEASURED, and at SE 0.011 an improvement of 0.03 A or more would have been seen:
underpowered below that, null above it. The idea is CLOSED as an accuracy step. (2) What the
energy DOES do, and it is the interesting half: the converged relaxed energy beats a uniformly
random choice among the same five solutions by -0.102 A [fold -0.127, -0.075], 91W/35L, 5/5
folds, 2.29x MDE, MEASURED. The production objective (CA-RMSD to the cloud plus 0.3 x ramah)
beats the random pick by 0.108 (2.65x MDE). So the converged all-atom energy carries the SAME
discriminating power among the projection's branches as the 2D torsion prior, and adds nothing
to it: on 42 targets it picks the identical solution, on the rest it trades one near-equivalent
branch for another (median difference 0.000). This is the first place in the record where an
AMBER quantity ranks a real discrete choice as well as the structural objective does; it is
also the first place where that choice is among only five candidates that all came from the
same cloud. (3) The raw single point `e0` is WORSE than the objective by +0.092 (1.10x MDE,
Type-M, 39W/58L/29T): the unrelaxed energy, which L23 showed is the builder's side-chain clash,
picks the wrong branch where the relaxed energy does not. Relaxation is what makes the energy a
usable discriminator here, the same operator fact as S20 L6 ("the relaxation is what makes the
AMBER objective defined"). (4) On the relaxed-chain basis the e1 pick is a null against the
production choice relaxed (-0.004, 0.09x MDE). (5) The ORACLE minimum over the five solutions is
-0.19 A below the production choice, but the valid best-of-5 null accounts for 130% of it and
the split-half transfer is 56% (-0.107): the branch degeneracy is real (a per-target-consistent
column exists) and it is worth about 0.1 A to a perfect chooser, which neither the objective
nor the energy is (each finds the per-target best on about 30% of targets, chance being 20%).

POWER. n = 126 with 42 exact ties; SE of the primary 0.011, MDE 0.031: a gain of 0.03 A would
have cleared. No positive result on the falsifier, so no seed replication is due (the arm has
no randomness; the relaxation is deterministic at threads = 1). The measured negative-control
contrast (energy beats random by 0.10) is not a proposal and needs no replication to be quoted
as a diagnostic. DISPOSITION: branch_select CLOSED as an accuracy lever; the record gains one
measured sentence for the report: "the converged force-field energy discriminates among the
projection's own branches exactly as well as the torsion prior already does, and the
unrelaxed energy does not."

---

## L89 -- CIS FLOOR, TIGHT FORM (ORACLE DIAGNOSTIC): THE NATIVE PROJECTED THROUGH THE PRODUCTION PROJECTION SITS 0.083 A FROM ITSELF (max 0.44, 0.043 WITH THE PRIOR OFF); L38's 0.347 A WAS THE OWN-TORSION UPPER BOUND; THE CONSTANT OMEGA COSTS ABOUT 0.04 A HERE AND THE PROJECTION'S 0.166 A COST IS NOT REPRESENTATION (2026-09-13, PH)

`s26/ph_cis.py floor2`, `s26/results/ph_cis_floor2.json` (complete 126/126), job
`s26/jobs_done/ph_cis_floor2.json` (exit 0, 2813 s, peak RSS 0.104 GB). Registered in
`s26/PREREG_cis.md` addendum 2 before it ran, as the tight form of L38's floor: the native CA
trace itself is projected through the PRODUCTION projection (`I.project`: ramah at 0.3,
multi-start, exact gradient) and its lam = 0 rung, and each is scored against the native.
ORACLE DIAGNOSTIC (the native is the input). BASIS: CA-RMSD of an ideal-geometry chain
against the native CA trace, on both sides; `chain_cost` is the production
`rmsd_arm - rmsd_avg` as in L38. `ST.fmt` verbatim:

      floor2 (native projected, lam 0.3) MINUS floor_ca (own-torsion rebuild)
        a 0.0833 (med 0.0596)   b 0.3468 (med 0.2722)   n=126
        effect -0.2636   median -0.1713   SE 0.0272   MDE 0.0762   effect/MDE -3.46
        iid  CI95 [-0.3177, -0.2126]
        fold CI95 [-0.2998, -0.2289]   folds same sign 5/5   per-fold 0:-0.203 1:-0.293 2:-0.268 3:-0.324 4:-0.242
        113W/13L/0T   worst degradation +0.2127 (9L1M)   p90 +0.0059   power 1.00  Type-M 1.00
        concentration: drop-top10 -0.1993 vs uniform-effect null p10/p50/p90 -0.2339/-0.1999/-0.1694 -> pctile 0.513
        VERDICT: BETTER
      floor2 lam 0.3 MINUS floor2 lam 0
        a 0.0833 (med 0.0596)   b 0.0429 (med 0.0228)   n=126
        effect +0.0404   median +0.0087   SE 0.0060   MDE 0.0169   effect/MDE +2.39
        iid  CI95 [+0.0292, +0.0522]
        fold CI95 [+0.0318, +0.0492]   folds same sign 5/5   per-fold 0:+0.029 1:+0.042 2:+0.058 3:+0.030 4:+0.042
        23W/103L/0T   worst degradation +0.2741 (6GIJ)   p90 +0.1276   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0452 vs uniform-effect null p10/p50/p90 +0.0370/+0.0449/+0.0533 -> pctile 0.522
        VERDICT: WORSE
      production chain cost MINUS floor2 (lam 0.3)
        a 0.1664 (med 0.0977)   b 0.0833 (med 0.0596)   n=126
        effect +0.0832   median +0.0214   SE 0.0197   MDE 0.0552   effect/MDE +1.51
        iid  CI95 [+0.0440, +0.1214]
        fold CI95 [+0.0440, +0.1331]   folds same sign 5/5   per-fold 0:+0.178 1:+0.069 2:+0.038 3:+0.101 4:+0.040
        57W/69L/0T   worst degradation +0.5706 (1RSW)   p90 +0.3935   power 0.99  Type-M 1.01
        concentration: drop-top10 +0.1171 vs uniform-effect null p10/p50/p90 +0.0913/+0.1167/+0.1418 -> pctile 0.510
        VERDICT: WORSE

    floor2, lam 0.3 (the production rung)   mean 0.0833 A (SE 0.0075), median 0.0596, p90 0.215, max 0.435 (6MBM), none above 0.5
    floor2, lam 0 (no torsion prior)        mean 0.0429 A (SE 0.0053), median 0.0228
    L38's own-torsion rebuild floor          mean 0.3468 A, max 1.474
    Spearman(floor2, max omega deviation) +0.372;  Spearman(floor2, chain cost) +0.114;  floor2 above the rebuild floor on 13 targets

READING. (1) The tight representation floor of the production manifold is SMALL: the nearest
ideal-trans chain to any dev native is 0.083 A away at the production rung (0.043 A with the
prior off), never more than 0.44 A. L38's 0.347 A was the own-torsion rebuild, an upper bound
that the registered addendum predicted would fall, and it fell by 0.264 A [fold -0.300,
-0.229], 113W/13L. The reason is arithmetic: rebuilding from the native's own phi/psi lets
every 2-degree omega error accumulate down the chain, while fitting phi/psi to the trace
absorbs those errors into compensating torsions. (2) The constant omega and the ideal bond
geometry therefore cost the projection at most a few hundredths on this instrument; the
Spearman with omega non-planarity drops from 0.83 (L38) to 0.37. The cis-peptide gap (L22:
zero cis on the instrument) and the non-planarity gap are both priced now: 0.000 A and about
0.04 A. (3) The ramah prior at 0.3 pulls the projection 0.040 A [+0.032, +0.049] away from the
native when the input IS the native (23W/103L): the prior is a bias toward typical torsions,
and on a perfect input it can only cost. On the production input (a coordinate average, not a
native) the prior is what chooses the Ramachandran-plausible branch (S9-2), so this is not a
proposal to remove it; it is the prior's price on the one input where it has nothing to fix.
(4) The production chain cost (0.166 A) exceeds the tight floor by +0.083 [+0.044, +0.133],
57W/69L, 5/5 folds, 1.51x MDE: what the projection pays on the production input is not
representation (0.08 of it is) but the operator's own displacement of a non-native cloud (S16
L27), and floor2's rho with the chain cost is +0.11.

POWER. SE 0.006 to 0.027 A on the three contrasts; the smallest MDE is 0.017 A. Nothing here is
a proposal; the two-bond-length projection stays worth 0.000 A on this instrument (L22) and the
design note in `s26/agentPH_FINDINGS.md` 1.4 stands. L38's number is superseded as "the floor"
by this one and is retained as the own-torsion upper bound; both are in the findings.

---

## L90 -- partial_recall_gradient (W, L77 extension): NO GRADIENT AT THE PRE-REGISTERED MDE (rho <= -0.25); THE STRONGEST COVARIATE, THE MAX PINNED IDENTITY TO THE TRAINING CORPUS, IS rho -0.205 WITH rmsd_arm (0.81x MDE, fold CI [-0.315, -0.099], PERMUTATION p 0.010, -0.176 AFTER PARTIALLING OUT LENGTH): SUGGESTIVE, NOT MEASURED, AND CONFOUNDED WITH RETRIEVAL BY CONSTRUCTION (2026-09-14, W)

Pre-registered in `s26/PREREG_partial_recall_gradient.md` (written before any covariate was
correlated with any RMSD). Covariates, native-free, sequences only: job `w_recall_cov` (exit 0,
15 s, peak 0.292 GB; `s26/results/w_selfcopy_recall_covariates.json`, complete 126/126): for each
dev target, against its own fold model's TRAINING corpus (`p_ladder.train_entries(fold)`, the
list `core.predict.train_fold` uses; the target's own chain asserted absent), I_long = max
pinned longer-normalised identity (mean 0.470, range 0.333 to 0.588, all below 0.6 as the fold
discipline requires), I_short = max containment (mean 0.657, range 0.50 to 1.00; the 1.00s are
the four self-copies of L44), L_kmer = the longest shared exact substring (mean 4.4 residues,
range 3 to 13). Gated: job `w_recall_endpoint` (exit 0, 45 s, 0.111 GB;
`s26/results/w_selfcopy_recall_endpoint.json`), Spearman against the production cache's
`rmsd_arm` (built chain, PRIMARY), `rmsd_avg` and `shipped`, fold-clustered bootstrap CI (5
folds, 4,000 draws), 500-draw label permutation one-sided in the direction of help, Bonferroni
alpha 0.017; the rank-partial correlation given chain length beside each.

    covariate   outcome     rho      fold CI95           perm p (help)   partial | n   verdict at the registered MDE 0.253
    I_long      rmsd_arm   -0.205   [-0.315, -0.099]    0.010           -0.176        no gradient (0.81x MDE)
    I_long      rmsd_avg   -0.210   [-0.312, -0.091]    0.014           -0.169        no gradient
    I_long      shipped    -0.230   [-0.280, -0.156]    0.004           -0.195        no gradient (0.91x)
    I_short     rmsd_arm   -0.157   [-0.254, -0.071]    0.036           -0.138        no gradient
    I_short     rmsd_avg   -0.136   [-0.229, -0.064]    0.068           -0.109        no gradient
    I_short     shipped    -0.175   [-0.252, -0.097]    0.026           -0.152        no gradient
    L_kmer      rmsd_arm   -0.088   [-0.198, +0.006]    0.174           -0.090        no gradient
    L_kmer      rmsd_avg   -0.081   [-0.192, +0.030]    0.178           -0.084        no gradient
    L_kmer      shipped    -0.140   [-0.244, -0.035]    0.056           -0.143        no gradient

Verdict by the pre-registered rule: no covariate reaches rho <= -0.25, so no gradient EXISTS at
the MDE; the registered expectation (all |rho| within 0.15) is exceeded by I_long on all three
bases (-0.205 to -0.230), whose fold CI excludes zero and whose permutation p passes Bonferroni
on `rmsd_arm` and `shipped`, at 0.81 to 0.91x the MDE: the correlation analogue of the Type-M
zone, sign consistent, magnitude not a result. Power: a gradient of |rho| >= 0.25 is excluded
(the MDE); one of 0.2 is not.

A confound the PREREG did not name, stated now: the fold model's training corpus IS the
retrieval library (out-of-fold peptides + fold fragments; `s12/instrument.py`), so a target with
a close relative in the corpus also has a closer BLOSUM retrieval pool, and a negative rho is
expected from retrieval alone (S7-12: BLOSUM has skill; S10-4's table). The covariate cannot
separate "the model recalls a training native" from "the pool holds a better window", and the
registered ORACLE mechanism check (the posterior's MAE against I_short) was conditioned on a
gradient at the MDE and did not run. What stands: a clean bill for the 0.6 threshold at the
pre-registered resolution (no near-copy below 0.6 moves the MLP by an amount the instrument can
call), with the honest addition that a weak sequence-proximity effect of either mechanism, of
order rho -0.2, is suggested on every basis. Replication: the covariates are deterministic; the
permutation and bootstrap are seeded (`s15.seed.stable_rng`). Not deployable: the covariate is
native-free, the outcome is the native; nothing selects.

## L91 -- memorisation_on_the_ladder (W, L77 extension; ORACLE DIAGNOSTIC): THE OWN-NATIVE MODELS TRAVEL 28% OF THE WAY TO THE NATIVE IN PROBABILITY SPACE AT COSINE 0.58 (77% AT COSINE 0.91 IN LOCATION SPACE, MAE 2.34 -> 0.82 A PER PAIR); THE DIRECTION-DISCOUNTED LADDER PREDICTS -0.36 A AGAINST A MEASURED -0.71 A (DISAGREEMENT +0.35, 1.5x MDE, 5/5 FOLDS): THE S25 L12 CURRENCY UNDER-PRICES A TRAINED OPERATOR; THE UNDISCOUNTED -2.15 x gam_eff LANDS WITHIN THE MDE (-0.60) (2026-09-14, W)

Pre-registered in `s26/PREREG_memorisation_on_the_ladder.md` (written before any leaked model's
gam_eff was computed). Job `w_ladder` (exit 0, 40 s, peak 0.274 GB;
`s26/results/w_selfcopy_ladder.json`, complete 504/504 (target, model) cells). ORACLE on every
line: gam_eff, cos and MAE are measured against the native's distances; the endpoint deltas are
L44 Part C's (leaked-model mean minus clean, `s26/results/w_selfcopy_endpoint.json :: C/rows`).
`gam_eff` and `cos` from `s26/p_ladder.progress` (the C2 ladder's code): probability space is
the S24 MASS construction <Q - P0, O - P0> / |O - P0|^2; location space is the S25 loc currency
on E[d].

    quantity (mean over 504 cells)                       value
    gam_prob / cos_prob / amp_prob                       0.281 / 0.584 / 0.473
    gam_loc / cos_loc                                    0.773 / 0.911
    MAE of E[d] vs native: clean -> leaked               2.339 -> 0.825 A per pair
    measured delta, cloud / built chain (L44)            -0.709 / -0.698 A
    ladder, naive  -2.1496 x gam_prob                    -0.603 A     (disagreement with the measured cloud delta +0.106, inside the MDE 0.244)
    ladder, discounted  -2.1496 x gam_prob x cos_prob    -0.364 A     (disagreement +0.346, 1.42x the MDE; the registered falsifier FIRES)
    corr over 504 cells: pred_discounted vs delta_cloud  +0.414   (S25 L12's 51 re-readings: +0.054)
                         gam_prob vs delta_cloud         -0.436
                         MAE change vs delta_cloud       +0.721

The registered per-target contrast (`ST.fmt` verbatim; the ladder's prediction minus the
measured cloud delta, four-model means, n = 126; positive = the ladder predicts LESS gain than
measured):

  ladder prediction (discounted, -2.1496 x gam x cos) minus MEASURED cloud delta, per target (ORACLE)
    a -0.3635 (med -0.3440)   b -0.7094 (med -0.3121)   n=126
    effect +0.3459   median -0.0122   SE 0.0826   MDE 0.2313   effect/MDE +1.50
    iid  CI95 [+0.1905, +0.5148]
    fold CI95 [+0.2058, +0.4719]   folds same sign 5/5   per-fold 0:+0.094 1:+0.575 2:+0.273 3:+0.362 4:+0.429
    65W/61L/0T   worst degradation +4.1162 (2MQ2)   p90 +1.6348   power 0.99  Type-M 1.01
    concentration: drop-top10 +0.4324 vs uniform-effect null p10/p50/p90 +0.3213/+0.4306/+0.5470 -> pctile 0.507
    VERDICT: WORSE

Reading. The registered expectation held: the direction-discounted currency of S25 L12
(gam_eff x cos) under-prices the one trained operator that moved far, by 0.35 A on the mean,
1.5x its MDE, 5/5 folds; the median disagreement is -0.01 (the mean is carried by the targets
where memorisation is large and off-axis: 65W/61L). The undiscounted ladder slope applied to
gam_prob alone lands within the MDE of the measured gain (-0.60 vs -0.71). So for a TRAINED
posterior the endpoint responds to the full probability-space move, not only to its projection
onto the native's direction; the L12 discount was derived from re-readings that move small
distances at low cosine, and it does not transfer to a large trained move. Two consequences,
both ORACLE and both for Proposal C's arithmetic rather than for any deployable arm: (i) the
"-2.15 A per unit gamma" transfer function is a usable currency for a trained prior when gam_eff
is measured in probability space and quoted without the cosine discount, but with the cosine
reported (this operator: cos 0.58); (ii) MAE change tracks the endpoint at +0.72 here, so S7-3's
"MAE does not price selection" is a statement about between-model differences of a few tenths
of an Angstrom per pair, not about a 1.5 A per pair move. What the leaked models do to the
posterior: E[d] travels 77% of the way to the native's distances at cosine 0.91 (they recall
the native's distance matrix nearly along the right direction), while the 17-bin probability
vectors travel only 28% at cosine 0.58 (the recalled mass sits in bins adjacent to the native's,
not on it): the location currency and the probability currency disagree about how far the
operator went, and the endpoint follows the location. Replication: deterministic (pinned
models, fixed natives); nothing to replicate. Deviations from the PREREG: none. Not deployable
by construction (own-native training); the value is the calibration of the currency the
record uses to price a better prior.
## L92 -- SEVEN LAUNCHER WAITERS STARTED BEFORE v2.3 WERE STARVED FOR TWO HOURS BY THE OLD FIXED CAP; RELAUNCHED UNDER THE FILE CAP (2026-09-14 00:10, coordinator)

Waiters created before 22:25 read the cap of four at import and never saw the file cap of six
(L83); with five or six jobs registered under the new cap they could never pass their own
test. Found at 00:10 with no child and a start time of 22:05 to 22:07: p_eval_esm8m2 (P),
pytest_slow_equivalence and pytest_slow_integration (I, AMBER), a2_dla_a1 and a4_var_boot (Q),
a_strain_vs_spread (A), e_report_check5 (E). Each was terminated before it had launched
anything and its exact command re-run under jobrun v2.3 (logs `s26/logs/relaunch_<name>.out`);
the owners' done-file watchers are unaffected because the job names are unchanged. Lesson for
the hygiene list: a launcher's cap must be read at every tick from the day it is written.

---


## L93 -- LANE P, C2 RUNG MIX: THE POOL-HISTOGRAM MIXTURE CHOOSES lam = 0 ON 5/5 FOLDS: THE SHIPPED POSTERIOR IS ITS OWN LEAVE-FOLD-OUT OPTIMUM; EVERY lam > 0 IS MONOTONICALLY WORSE ON THE BUILT CHAIN, lam = 1 (TYPICALITY) BY +0.733 A (2026-09-14 00:10, lane P)

Artefacts: `s26/results/p_ladder_mix_s0.json` (126 rows, complete), `s26/results/p_ladder_report_mix_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_mix`: exit 0, 4,332 s, peak RSS 0.13 GB; eight mixtures P = (1-lam) P_shipped + lam H_pool per target, H_pool the K = 500 pool's own 17-bin per-pair histogram (+1e-3 floor), lam in {0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1}; lam = 0 asserted bit-exact to the incumbent's risk table on every target. Mean built chain by lam: 3.2126 / 3.2286 / 3.2295 / 3.2555 / 3.2871 / 3.3451 / 3.6556 / 3.9455 (monotone). Leave-fold-out lam on the built chain = 0 on all five folds, so the deployable arm IS the incumbent (126 ties; the FLAGs are the degenerate all-zero case). On selection the fold-wise lam is nonzero on three folds and the arm is +0.039 (0.21x MDE, 57 ties). Per-target ORACLE over the 8 lams: -0.399 A, but `ST.best_of_k_within`'s valid across-target null is -0.729 (183% accounted; k_eff 5.8): an order statistic. (The split-half figure of -0.165 measures that lam = 0 is the consistently best COLUMN against the row mean, which includes lam = 1 at +0.73; it is not per-target transfer.) This is the pre-registered prediction (PREREG_C2, S7-3 `gain 0`, S10-2, S19 L14): the pool's own distance statistics add nothing to the prior as a distribution, and the more of them enter, the worse. lam = 1 reproduces S7-3's typicality penalty on the built chain (+0.733 here; S7-3 measured +0.280 on selection over the older instrument).

```
MIX: leave-fold-out lam per fold {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0} ; per-target ORACLE over 8 lams: gain -0.3988, valid null -0.7285 (share 1.83), split-half -0.1647 (41% of oracle), k_eff 5.81
     mean arm by lam: 0:3.2126  0.05:3.2286  0.1:3.2295  0.2:3.2555  0.35:3.2871  0.5:3.3451  0.75:3.6556  1:3.9455
  mix vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2126 (med 2.9661)   b 3.2126 (med 2.9661)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  mix vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0483 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  mix vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4925 (med 3.4817)   b 3.4540 (med 3.4779)   n=126
    effect +0.0385   median +0.0000   SE 0.0643   MDE 0.1801   effect/MDE +0.21
    iid  CI95 [-0.0837, +0.1619]
    fold CI95 [-0.0487, +0.1590]   folds same sign 2/5   per-fold 0:+0.267 1:-0.037 2:-0.052 3:+0.080 4:-0.050
    33W/36L/57T   worst degradation +2.4137 (8T62)   p90 +0.6572   power 0.09  Type-M 4.03
    concentration: drop-top10 +0.1772 vs uniform-effect null p10/p50/p90 +0.1057/+0.1722/+0.2422 -> pctile 0.543
    VERDICT: NOT MEASURED (|effect| 0.0385 <= its own MDE 0.1801, 0.21x)
  strata [arm]: len 9-10 n=20 +0.000 (SE 0.000, med +0.000, 0W/0L) | len 11-12 n=32 +0.000 (SE 0.000, med +0.000, 0W/0L) | len 13-14 n=37 +0.000 (SE 0.000, med +0.000, 0W/0L) | len 15-16 n=37 +0.000 (SE 0.000, med +0.000, 0W/0L) | FAIL18 n=18 +0.000 (SE 0.000, med +0.000, 0W/0L) | other108 n=108 +0.000 (SE 0.000, med +0.000, 0W/0L)
  strata [cloud]: len 9-10 n=20 +0.000 (SE 0.000, med +0.000, 0W/0L) | len 11-12 n=32 +0.000 (SE 0.000, med +0.000, 0W/0L) | len 13-14 n=37 +0.000 (SE 0.000, med +0.000, 0W/0L) | len 15-16 n=37 +0.000 (SE 0.000, med +0.000, 0W/0L) | FAIL18 n=18 +0.000 (SE 0.000, med +0.000, 0W/0L) | other108 n=108 +0.000 (SE 0.000, med +0.000, 0W/0L)
  strata [sel]: len 9-10 n=20 +0.039 (SE 0.059, med +0.000, 6W/6L) | len 11-12 n=32 -0.043 (SE 0.088, med +0.000, 10W/5L) | len 13-14 n=37 +0.066 (SE 0.140, med +0.000, 8W/13L) | len 15-16 n=37 +0.081 (SE 0.150, med +0.000, 9W/12L) | FAIL18 n=18 -0.071 (SE 0.143, med +0.000, 7W/4L) | other108 n=108 +0.057 (SE 0.071, med +0.000, 26W/32L)
  gamma-equivalent: gam_eff prob-space +0.0000 at cos +nan ; loc-space +0.0000 at cos +nan ; MAE 2.3386 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 5/5 ; verdict (arm): NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
```

---

## L94 -- ADVERSARY CHECK OF L84 AND L85 (W: window ensembling, window provenance): L84 STANDS; L85 STANDS WITH CAVEAT (AN ORACLE POOL FACT; THE DEPLOYABLE TEST H_P3 IS STILL TO RUN) (2026-09-14, A)

L84 (fixed-K ensembling of three BLOSUM keys). Pre-registered before the probe; basis stated
(built chain on the rebuild basis, point cloud carried, `sel` undefined for an ensemble and not
quoted); the zero-information control is matched in member count (three 75-resamples of the
production set); the primary is -0.0005 A at 0.02x MDE and the control +0.0042 at 0.14x, both
NOT MEASURED; every secondary is under its MDE, and the one that touches the Type-M edge (the
BLOSUM80 shortlist alone, +0.0237 at 0.70x with the fold CI above zero, 4/5) is flagged "if
anything worse", which is the right wording. Power: the design resolves 0.028 A and sees
nothing; the mandatory direction is closed at that resolution in its fixed-K form, the
widening-K form having been closed by S17 L12. Deterministic, nothing to replicate. STANDS.

L85 (window provenance). The census is native-free (codes and sequences); the contrast is
ORACLE on the single-window basis and labelled so. Pool: peptide-derived minus fragment -0.416
A, 3.18x MDE, fold CI [-0.447, -0.376], 101W/25L, 5/5, not concentrated: a clean ORACLE fact,
and the one S19 already implied (the peptide corpus carries the sequence-structure channel; the
fragments are 73% of every pool). Top-75: whole-or-terminal minus fragment -0.131 at 1.13x
(Type-M, flagged). Two caveats. (1) Nothing deployable is claimed yet: the class is native-free
but the contrast is ORACLE; the pre-registered ACHIEVABLE test H_P3 (class-weighted readout
against uniform and a permuted-weight control) is the only thing that could turn this into a
lever, and by the record (uniform readout optimal over two families, S25; 68% common-mode) the
expectation is null; the Adversary checks H_P3 when it lands. (2) Six class contrasts are
listed; the two the prereg named are the ones that count, and the four secondaries are under
their MDEs. STANDS WITH CAVEAT.

---

## L95 -- ADVERSARY CHECK OF L86 AND L87 (PH: the reject on the built chain; the C3 replication): BOTH STAND WITH CAVEAT; THE TOWARD-MEMBER GAIN IS A MEASURED, REPLICATED 0.02 A THAT UNDOES PART OF THE PROJECTION'S OWN MOVE (2026-09-14, A)

L86 (steric reject, built chain). The cross-basis replication L54 asked for: R@1e4 +0.248 A
against the re-projected shipped top-75 (1.17x MDE, Type-M, fold CI [+0.120, +0.348]) and
+0.157 against the same count rejected at random (1.15x, Type-M), S@1e4 +0.104 (1.11x, 5/5),
S vs RANDS +0.055 (0.73x, NOT MEASURED). Same caveats as L54, plus one: on the built chain the
R-arm primary holds on 4 of 5 folds (fold 1 -0.002), not 5/5, so it does not meet the standing
rule's fold clause on this basis; the S arm does (5/5), the 1e3 arms are measured on both bases
and the dose is monotone, so "closed on both bases, with a sign" stands on those. The harm is
again tail-carried (median +0.0016, worst +4.20 on 8T61). The anchor is the re-projected
shipped set (addendum 2), the right comparator for a re-projected arm. STANDS WITH CAVEAT.

L87 (C3 stage-1 replication). New seeds, reversed order, fold labels fixed: every replicated
quantity lands inside L39's fold CI. The toward-member control minus do-nothing is now -0.0207
[-0.0250, -0.0166], 5/5, 1.77x MDE, 94W/32L, replicated: by the discipline this IS a measured
result, and it is native-free (a move of AMBER's own magnitude toward a random member of the
shipped pool). PH reads it as a property of the projection's displacement, not of physics, and
not a proposal; I agree, with the mechanism stated in numbers so the presenter cannot mistake
it: the built chain sits 0.166 A further from the native than the point cloud it was projected
from (L38 / EXAMINATION H), the pool members surround that cloud, so a 0.22 A step toward a
member is a partial reversal of the projection's own displacement and recovers about an eighth
of it. The natural comparator, the same step toward the cloud itself, was not run and would
say whether "member" matters at all; and a moved chain is no longer an ideal-geometry chain,
so the emitted object changes basis. Not a lever; a measurement of what the projection costs.
The AMBER-minus-random Type-M status (1.02x in the replication) is unchanged. STANDS WITH CAVEAT.

---

## L96 -- ADVERSARY CHECK OF L88 (branch_select): STANDS; THE 56% SPLIT-HALF TRANSFER IS AGAINST THE COLUMN MEAN AND LANDS ON THE PRODUCTION CHOICE, NOT ABOVE IT (2026-09-14, A)

Operator native-free (the converged e1 of each of the five projection solutions; ties by
`ST.argmin_tied`, 0 ties reported); ORACLE evaluation only. The falsifier did not fire: e1 pick
minus production +0.0055 at 0.18x MDE with 42 exact ties (the e1 pick IS the production choice
on a third of the targets), the fold CI on the harmful side, 4/5 folds: NOT MEASURED, closed as
an accuracy step with the power stated. The e1 pick beats a random branch by -0.102 (2.29x,
5/5) and the raw single point picks worse (+0.092, Type-M): the relaxed energy carries the
same discriminating power among branches as the 2D torsion prior, and the raw energy carries
the builder's side-chain clash (L23). Clean on leakage, ties, both CIs, concentration and power.

One clarification, so the 56% is not read as a hidden lever. `ST.best_of_k_within` over the five
columns gives an ORACLE per-target minimum of -0.190 (3.130 A) against the row mean, with a
split-half transfer of -0.107 (56%, k_eff 4.63, argmin counts [33, 24, 39, 18, 12]). A
transfer is measured against the COLUMN MEAN, which here is the random pick (3.320 A); a fixed
best column therefore lands at about 3.21 A, which is the production objective's own choice
(3.213). So "a per-target-consistent column exists" means one start is systematically better
than a random start, and the objective already captures that; a fixed-start rule would at best
tie production. The ORACLE 0.083 A between production and the per-target best branch is the
order statistic L88 says it is. STANDS.

---

## L97 -- ADVERSARY CHECK OF L89, L90, L91, L93 (the cis floor's tight form; the recall gradient; memorisation on the ladder; the mix rung): ALL STAND; L38's 0.347 A IS NOW TO BE QUOTED AS THE OWN-TORSION UPPER BOUND BESIDE THE TIGHT 0.083 A (2026-09-14, A)

L89 (ORACLE DIAGNOSTIC). The native projected through the production projection sits 0.083 A
from itself (0.043 with the prior off; max 0.44); floor2 minus L38's own-torsion floor -0.264
(3.46x, 113W/13L); the prior costs +0.040 when the input is the native (2.39x, 23W/103L); the
chain cost exceeds floor2 by +0.083 (1.51x, 5/5). L38 declared its 0.347 an upper bound and
registered floor2 in advance, so nothing is retracted; every quotation of "0.35 A" now carries
"own-torsion upper bound; the tight floor is 0.08 A, and the projection's 0.166 A cost is the
operator's displacement, not representation" (added to the slide-qualifier table). STANDS.

L90 (native-free covariates against the ORACLE label). No covariate reaches the registered
|rho| >= 0.25; the strongest, max pinned identity to the corpus, is -0.205 (0.81x, fold CI
[-0.315, -0.099], permutation p 0.010; -0.176 after partialling n): the correlation analogue of
the Type-M zone, and W says so. The unnamed confound W adds (the training corpus IS the
retrieval library, so "recall" and "a nearer pool window" are not separable here) is the right
one and is stated before any reading. Power stated (|rho| >= 0.25 excluded at ~0.8). STANDS.

L91 (ORACLE DIAGNOSTIC, own-native models). The direction-discounted ladder (-2.1496 x gam x
cos) over-predicts the measured cloud delta by +0.346 A (1.50x MDE, 5/5, 65W/61L); the
registered falsifier fires, so the S25 L12 currency is calibrated on small re-readings only and
must not be used to price a large off-axis move. Not deployable by construction; the value is
the calibration, which L62 and L79 already apply ("gam_eff is never a prediction of the
endpoint"). STANDS.

L93 (mix rung). Leave-fold-out chooses lam = 0 on 5/5 folds, so the arm is the shipped
posterior on every basis (126 ties, the degenerate FLAG); every lam > 0 is monotonically worse
(+0.733 at lam = 1, typicality). The per-target ORACLE over the eight lams (-0.399) is 1.83x
covered by the valid across-target null and its 41% split-half transfer is the trivial statement
that lam = 0 beats the mean over lams; no lam other than 0 transfers. Nested by construction
(the lam chosen leave-fold-out), so no regression-to-the-mean concern. S7-3 / S19 L14's
foreshadowed null, measured with the falsifier. STANDS.

---

## L98 -- THE RELAUNCHED SLOW-TEST JOBS LOST VERIFY_SLOW=1: ONE NULL RUN SET ASIDE, BOTH TIERS RELAUNCHED WITH THE FLAG INSIDE THE COMMAND; THE AST GATE PASSES; THE REBUILD IS QUEUED (2026-09-14 00:17, lane I)

L92's relaunch re-ran my two slow-test commands exactly, but `VERIFY_SLOW=1` had been in the
launching shell's environment, not in the command, so the relaunched jobs ran without it:
`pytest_slow_integration` (00:13:41) exited 0 in 5 s with all 8 items SKIPPED (a null run, its
records set aside as `s26/results/pytest_slow_integration.NULLRUN_no_VERIFY_SLOW.{json,xml}`,
not counted), and the waiting `pytest_slow_equivalence` would have done the same, so I
terminated that waiter (pid 2664) and my own chain launcher that was gated on those names.
Relaunched at 00:17 as `pytest_slow_equivalence2` (AMBER, est 1.2 GB) and
`pytest_slow_integration2` (AMBER, est 1.8 GB) with the flag set inside the command
(`python -c "os.environ['VERIFY_SLOW']='1'; pytest.main([...])"`), so no relaunch can drop it
again; the verify chain launcher re-gated on the new names. Lesson for the hygiene list beside
L92's: a job's environment belongs in its command.

Also done without the box: the final AST gate, `s26/i_ast_check.py --ref ae86a124` over every
production module (`core/` 11, the 24 root modules, `s5/ s7/ s8/ s9/` 20 = 55 files,
docstrings stripped): **50 AST-identical; 5 differ, and the five are exactly the ledgered edits**
(`core/amber.py`, `core/bench.py`, `core/cache.py`, `core/predict.py`: one unused import each,
L10; `core/data.py`: the `norm` keyword on `identity` / `identity_many` and the removed dead
`_seq_index`, L15). Table: `s26/results/ast_gate_ae86a124.txt` (commit `[see git log]`).
The frozen results-lab rebuild is queued as governed job `resultslab_rebuild` (CPU, est 1.2 GB)
behind `s26/i_resultslab_rebuild.py pre` (snapshot of the tracked `results/summary` and the
sha256 + ATOM-record sha256 of all 2,142 tracked PDBs); its `post` comparator was self-checked on
the untouched tree (2016/2016 RMSD values exact, 2142/2142 PDBs byte-identical).

---

## L99 -- LANE P, C2 RUNG ESM8M: THE 8M LANGUAGE MODEL IN PLACE OF THE 650M IS WORSE THAN THE SHIPPED PRIOR BY +0.242 A ON THE BUILT CHAIN (1.17x MDE, 5/5 FOLDS) AND NO BETTER THAN NO ESM AT ALL (+0.034 vs noesm, 0.15x): THE ESM CHANNEL'S VALUE IS IN THE 650M MODEL, NOT IN 'ANY LANGUAGE MODEL' (2026-09-14 00:19, lane P)

Artefacts: `s26/results/p_ladder_esm8m_s0.json` (126 rows, complete), `s26/results/p_ladder_report_esm8m_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_esm8m` (killed by the stall breaker at 110/126, L74) resumed as `p_eval_esm8m2`: exit 0, 251 s, peak RSS 0.445 GB; the result file completed with `complete: true` on the resumed run and is used only in that form. Models `s26/models/p_ladder/esm8m_fold{0..4}_s0.pt`: esm2_t6_8M_UR50D reps (320-d, PCA-32 refit per fold) and ITS OWN contact head, shipped architecture (`p_train_esm8m` peak 1.85 GB). PREREG_B2's downward size point. Read as the monotone test PREREG_B2 declared: esm8m is worse than pca32 (the 650M) by +0.242 A [fold +0.113, +0.385], 5/5 folds, 1.17x MDE (Type-M zone: the sign is established, the magnitude is inflated ~1.06x), and indistinguishable from noesm (+0.034, 0.15x MDE, median -0.016) and from conly on selection (+0.113, 0.46x). So the size axis is NOT flat: an 8M model carries none of the ESM channel and the 650M carries all of it that this instrument can see; the ladder cannot say whether 3B would carry more (infeasible here, L13), only that the slope from 8M to 650M is about -0.24 A per 1.9 decades of parameters on the built chain and -0.30 on selection. This corrects the reading S7-11's 'any reduction of it' invited: it is any reduction of the 650M representation, not any language model. gam_eff +0.107 prob at cos 0.24 (positive while worse, the sixth instance). Strata: worse in every length band and on other108; on FAIL18 -0.37 (12W/6L, SE 0.28), the same pattern as every other worse rung (findings section 9, H1).

```
  esm8m vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.4544 (med 3.3360)   b 3.2126 (med 2.9661)   n=126
    effect +0.2418   median +0.0442   SE 0.0740   MDE 0.2074   effect/MDE +1.17
    iid  CI95 [+0.1006, +0.3937]
    fold CI95 [+0.1134, +0.3853]   folds same sign 5/5   per-fold 0:+0.165 1:+0.201 2:+0.470 3:+0.385 4:+0.036
    52W/74L/0T   worst degradation +3.9336 (1CEK)   p90 +1.0899   power 0.90  Type-M 1.06
    concentration: drop-top10 +0.3418 vs uniform-effect null p10/p50/p90 +0.2441/+0.3351/+0.4368 -> pctile 0.535
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
  esm8m vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.2595 (med 3.0590)   b 3.0483 (med 2.8373)   n=126
    effect +0.2111   median +0.0407   SE 0.0688   MDE 0.1928   effect/MDE +1.09
    iid  CI95 [+0.0829, +0.3486]
    fold CI95 [+0.0765, +0.3424]   folds same sign 4/5   per-fold 0:+0.175 1:+0.220 2:+0.403 3:+0.328 4:-0.015
    53W/73L/0T   worst degradation +3.4867 (1CEK)   p90 +1.1721   power 0.87  Type-M 1.08
    concentration: drop-top10 +0.3071 vs uniform-effect null p10/p50/p90 +0.2193/+0.3046/+0.3948 -> pctile 0.516
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.08x]
  esm8m vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.6786 (med 3.5757)   b 3.4540 (med 3.4779)   n=126
    effect +0.2246   median +0.0355   SE 0.0913   MDE 0.2557   effect/MDE +0.88
    iid  CI95 [+0.0470, +0.4071]
    fold CI95 [+0.1221, +0.3046]   folds same sign 5/5   per-fold 0:+0.276 1:+0.231 2:+0.330 3:+0.293 4:+0.037
    57W/65L/4T   worst degradation +3.9038 (1CEK)   p90 +1.6314   power 0.69  Type-M 1.21
    concentration: drop-top10 +0.3882 vs uniform-effect null p10/p50/p90 +0.2703/+0.3837/+0.4991 -> pctile 0.519
    VERDICT: NOT MEASURED (|effect| 0.2246 <= its own MDE 0.2557, 0.88x)
  strata [arm]: len 9-10 n=20 +0.136 (SE 0.081, med +0.096, 6W/14L) | len 11-12 n=32 +0.153 (SE 0.137, med -0.007, 16W/16L) | len 13-14 n=37 +0.263 (SE 0.167, med +0.025, 17W/20L) | len 15-16 n=37 +0.355 (SE 0.143, med +0.063, 13W/24L) | FAIL18 n=18 -0.069 (SE 0.109, med -0.018, 10W/8L) | other108 n=108 +0.294 (SE 0.084, med +0.052, 42W/66L)
  strata [cloud]: len 9-10 n=20 +0.092 (SE 0.093, med +0.011, 9W/11L) | len 11-12 n=32 +0.140 (SE 0.131, med +0.025, 14W/18L) | len 13-14 n=37 +0.249 (SE 0.154, med +0.064, 15W/22L) | len 15-16 n=37 +0.299 (SE 0.128, med +0.046, 15W/22L) | FAIL18 n=18 -0.083 (SE 0.101, med -0.016, 9W/9L) | other108 n=108 +0.260 (SE 0.078, med +0.049, 44W/64L)
  strata [sel]: len 9-10 n=20 +0.012 (SE 0.171, med +0.056, 9W/10L) | len 11-12 n=32 +0.210 (SE 0.165, med +0.085, 14W/18L) | len 13-14 n=37 +0.353 (SE 0.191, med +0.000, 17W/18L) | len 15-16 n=37 +0.224 (SE 0.180, med +0.036, 17W/19L) | FAIL18 n=18 +0.128 (SE 0.156, med +0.053, 6W/10L) | other108 n=108 +0.241 (SE 0.103, med +0.034, 51W/55L)
  gamma-equivalent: gam_eff prob-space +0.1139 at cos +0.230 ; loc-space +0.3579 at cos +0.375 ; MAE 2.4666 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 5/5 ; verdict (arm): WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
```

---

## L100 -- THE VALIDITY AXIS OF THE PRODUCTION RELAXATION, MEASURED: HEAVY-ATOM CLASHES 0.44 PER TARGET (34 TARGETS) TO 0.008 (1), CLOSEST PAIR 2.36 TO 2.78 A, AT 1.3% BOND / 2.5% ANGLE STRAIN AND 6.6 DEG OF OMEGA, NO RAMACHANDRAN GAIN; RE-RUN BIT-IDENTICAL TO THE CACHE ON 126/126 (2026-09-14, PH)

`s26/ph_validity.py run | report`, `s26/results/ph_validity.json` (complete 126/126), job
`s26/jobs_done/ph_validity.json` (exit 0, 2568 s, peak RSS 0.274 GB, AMBER; held 620 s at the
cap). Pre-registered in `s26/PREREG_validity_axis.md` before the run. Native-free: the panel
(`s16.energy_lib.panel`) reads no native; folds enter only the CI.

GATE, per target: the production relaxation (`refine_coords(k=10, steps=0)` on the built chain
from the cached `phi`, `psi`, exactly `core.pipeline._relax_inner`) re-run here reproduces the
cached `amber_ca` and `amber_e1` on 126 of 126 targets to max |dCA| = 0.0 A and max |dE| = 0.0
kcal/mol. The panel below therefore describes the DEPLOYED emission before and after its own
relaxation, not a re-implementation. Zero-information reference: the constant alpha-helix
(phi -63, psi -42; rama 1.000 and zero clashes by construction, S16 L27), printed beside both.

    axis              built chain   relaxed chain   constant helix
    n_clash_2A        built   0.4444   relaxed   0.0079   helix   0.0000
    n_clash_2p6A      built   3.4524   relaxed   0.1190   helix   0.0000
    min_heavy         built   2.3629   relaxed   2.7846   helix   3.0792
    bond_strain       built   0.0000   relaxed   0.0130   helix   0.0000
    angle_strain      built   0.0000   relaxed   0.0246   helix   0.0000
    rama_favoured     built   0.9302   relaxed   0.9114   helix   1.0000
    rama_outlier      built   0.0242   relaxed   0.0327   helix   0.0000
    cis_frac          built   0.0000   relaxed   0.0033   helix   0.0000
    chirality_L_frac  built   1.0000   relaxed   1.0000   helix   1.0000
    omega_dev         built   0.0000   relaxed   6.6120   helix   0.0000
    targets with any heavy-atom pair below 2.0 A: built 34, relaxed 1;  relaxed bond_strain above 0.05 on 2 targets (2BP4, 9KAR)

`ST.fmt` verbatim, relaxed minus built (negative = lower after), the seven axes that move:

      n_clash_2A: relaxed minus built (negative = lower after)
        a 0.0079 (med 0.0000)   b 0.4444 (med 0.0000)   n=126
        effect -0.4365   median +0.0000   SE 0.0760   MDE 0.2129   effect/MDE -2.05
        iid  CI95 [-0.5873, -0.2937]
        fold CI95 [-0.5546, -0.3172]   folds same sign 5/5   per-fold 0:-0.560 1:-0.522 2:-0.240 3:-0.565 4:-0.333
        34W/0L/92T   worst degradation +0.0000 (1A13)   p90 +0.0000   power 1.00  Type-M 1.00
        concentration: drop-top10 -0.2586 vs uniform-effect null p10/p50/p90 -0.3448/-0.2586/-0.1724 -> pctile 0.473
        VERDICT: BETTER
      n_clash_2p6A: relaxed minus built (negative = lower after)
        a 0.1190 (med 0.0000)   b 3.4524 (med 1.0000)   n=126
        effect -3.3333   median -1.0000   SE 0.3606   MDE 1.0103   effect/MDE -3.30
        iid  CI95 [-4.0556, -2.6429]
        fold CI95 [-3.9916, -2.7385]   folds same sign 5/5   per-fold 0:-3.880 1:-3.783 2:-2.400 3:-4.217 4:-2.633
        81W/0L/45T   worst degradation +0.0000 (1A1P)   p90 +0.0000   power 1.00  Type-M 1.00
        concentration: drop-top10 -2.5517 vs uniform-effect null p10/p50/p90 -3.0086/-2.5517/-2.1293 -> pctile 0.497
        VERDICT: BETTER
      min_heavy: relaxed minus built (negative = lower after)
        a 2.7846 (med 2.8160)   b 2.3629 (med 2.4744)   n=126
        effect +0.4217   median +0.3380   SE 0.0326   MDE 0.0912   effect/MDE +4.62
        iid  CI95 [+0.3586, +0.4871]
        fold CI95 [+0.3680, +0.4783]   folds same sign 5/5   per-fold 0:+0.489 1:+0.412 2:+0.376 3:+0.508 4:+0.345
        11W/115L/0T   worst degradation +1.8827 (1I93)   p90 +0.8626   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.4631 vs uniform-effect null p10/p50/p90 +0.4197/+0.4625/+0.5075 -> pctile 0.507
        VERDICT: WORSE
      bond_strain: relaxed minus built (negative = lower after)
        a 0.0130 (med 0.0106)   b 0.0000 (med 0.0000)   n=126
        effect +0.0129   median +0.0106   SE 0.0012   MDE 0.0034   effect/MDE +3.78
        iid  CI95 [+0.0110, +0.0157]
        fold CI95 [+0.0112, +0.0150]   folds same sign 5/5   per-fold 0:+0.012 1:+0.016 2:+0.011 3:+0.015 4:+0.011
        0W/126L/0T   worst degradation +0.1206 (2BP4)   p90 +0.0130   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0132 vs uniform-effect null p10/p50/p90 +0.0116/+0.0131/+0.0151 -> pctile 0.545
        VERDICT: WORSE
      angle_strain: relaxed minus built (negative = lower after)
        a 0.0246 (med 0.0234)   b 0.0000 (med 0.0000)   n=126
        effect +0.0246   median +0.0234   SE 0.0005   MDE 0.0013   effect/MDE +18.37
        iid  CI95 [+0.0236, +0.0255]
        fold CI95 [+0.0239, +0.0254]   folds same sign 5/5   per-fold 0:+0.026 1:+0.025 2:+0.023 3:+0.024 4:+0.024
        0W/126L/0T   worst degradation +0.0501 (1D6X)   p90 +0.0309   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0252 vs uniform-effect null p10/p50/p90 +0.0245/+0.0252/+0.0259 -> pctile 0.508
        VERDICT: WORSE
      omega_dev: relaxed minus built (negative = lower after)
        a 6.6120 (med 5.6118)   b 0.0000 (med 0.0000)   n=126
        effect +6.6120   median +5.6118   SE 0.4332   MDE 1.2135   effect/MDE +5.45
        iid  CI95 [+5.8798, +7.5240]
        fold CI95 [+6.0008, +7.2286]   folds same sign 5/5   per-fold 0:+7.583 1:+6.470 2:+5.445 3:+6.420 4:+7.032
        0W/126L/0T   worst degradation +48.0241 (1D6X)   p90 +9.6210   power 1.00  Type-M 1.00
        concentration: drop-top10 +6.9076 vs uniform-effect null p10/p50/p90 +6.3454/+6.8782/+7.5429 -> pctile 0.524
        VERDICT: WORSE
      rama_favoured: relaxed minus built (negative = lower after)
        a 0.9114 (med 1.0000)   b 0.9302 (med 1.0000)   n=126
        effect -0.0188   median +0.0000   SE 0.0080   MDE 0.0225   effect/MDE -0.84
        iid  CI95 [-0.0355, -0.0037]
        fold CI95 [-0.0334, +0.0013]   folds same sign 4/5   per-fold 0:-0.020 1:+0.020 2:-0.035 3:-0.014 4:-0.037
        27W/16L/83T   worst degradation +0.1818 (2LM8)   p90 +0.0833   power 0.65  Type-M 1.24
        concentration: drop-top10 +0.0003 vs uniform-effect null p10/p50/p90 -0.0083/+0.0003/+0.0086 -> pctile 0.500
        VERDICT: NOT MEASURED (|effect| 0.0188 <= its own MDE 0.0225, 0.84x)

READING, the numbers "validity step" rests on. (1) The relaxation removes the builder's
clashes: heavy-atom pairs below 2.0 A fall from 0.444 per target (34 targets affected) to
0.008 (1 target), fold CI [-0.555, -0.317], 34W/0L/92T; contacts below 2.6 A fall from 3.45
to 0.12 per target, 81W/0L/45T; the closest heavy-atom pair moves from 2.36 to 2.78 A, 115 of
126 targets. The registered falsifier ("clashes not halved") does not fire. (2) It pays in
covalent geometry, which the built chain had ideal by construction: bond strain 0.0 to 1.3%
relative (0/126 unchanged; 12% on 2BP4), angle strain 0.0 to 2.5%, omega non-planarity 0.0 to
6.6 degrees mean (48 degrees on 1D6X; the cis fraction rises to 0.3%, two targets). (3) It buys
no Ramachandran: favoured 0.930 to 0.911, outliers 0.024 to 0.033, both NOT MEASURED. This is
the S16 L27 finding on the production input: the restraint k = 10 on N/CA/C holds the backbone
torsions where they were and lets the force field fix the packing by bending bonds and
angles. (4) The constant helix beats the relaxed chain on every axis (0 clashes, 3.08 A minimum
separation, ideal covalent geometry, rama 1.000): a validity statistic a zero-information
reference maximises is not evidence about the force field, so the defensible statement is
the conjunction, as S16 wrote it: "removes the builder's clashes (34 targets to 1) while moving
the CA trace 0.220 A and holding the torsions, at a covalent price of 1.3% bond and 2.5% angle
strain and 6.6 degrees of omega".

THE SENTENCE, for `s26/C3_RESULT.md` and the presentation: the relaxation is a validity step
that turns 34 emissions with a sub-2 A heavy-atom overlap into 1, and 125 of 126 energies above
the 1000 kcal/mol gate into converged ones, at the price of 0.021 A of accuracy, 1.3% bond and
2.5% angle strain, 6.6 degrees of peptide-bond non-planarity, and a broken virtual bond on 2BP4
and 9KAR. POWER: n = 126, SEs 0.0005 to 0.36 on the axes; every claimed change is above 2x its
MDE with 5/5 folds; the three nulls (rama favoured 0.84x, outliers 0.55x, cis 0.46x) are
underpowered for effects below their MDEs (0.023, 0.015, 0.007) and measured against anything
larger. Descriptive; no replication is due.

---

## L101 -- CORRECTION TO L99: THE esm8m SELECTION CONTRAST IS +0.225 A (0.88x MDE, NOT MEASURED), NOT "-0.30" (2026-09-14 00:20, lane P)

L99's note quotes a selection slope of "-0.30" for the 8M-to-650M step; the block under it has
the right number: esm8m vs shipped on selection is +0.2246, SE 0.0913, MDE 0.2557, 0.88x, fold
CI [+0.122, +0.305], 5/5 folds, VERDICT NOT MEASURED (below its MDE; the fold CI excludes zero
but the MDE gate is not met, so no size is claimed). The built-chain figure (+0.242, 1.17x, WORSE)
is the one that carries the reading. Struck: "-0.30 on selection".

---

## L102 -- THE OPT-IN EQUIVALENCE TIER HAD NEVER BEEN ABLE TO RUN: ITS SUMMARY PARSER TOOK THE WRONG BRACE AND SKIPPED ALL THREE ARMS AFTER THE PIPELINES HAD FINISHED (2026-09-14 00:25, lane I)

Job `pytest_slow_equivalence2` (AMBER, `VERIFY_SLOW=1` inside the command, 210.6 s, peak RSS
0.632 GB, exit 0, `s26/logs/pytest_slow_equivalence2.log`): the three gated tests of
`tests/test_equivalence.py` all reported `SKIPPED: could not parse the harness summary`. The
cause is in the test, not the pipeline. `_run_arm` runs `python -m core.pipeline run --manifest
smoke8 --workers 1 [--backends legacy]` as a subprocess, asserts rc 0 (which passed: six pipeline
arms completed, resumed from the production and baseline caches), then locates the summary JSON
with `cp.stdout.rfind("{")`. The harness prints its summary with `indent=2`, so the last brace
in stdout opens a nested dict (`"stage_totals": { "n_top_mean": 75.0, ...`) and `json.loads`
raises; the `except` turned that into a skip. Reproduced offline on `verify/e2_baseline.log`
(the same CLI output): `rfind` fails, a line-anchored top-level brace parses. So the "13 skips"
of every recorded run include three tests that could not have run even with the flag set.

Fix (commit `7be8e4b0`, a test file, not the production path): the parser takes the last line
that is exactly `{`, and a parse failure is now a FAIL, not a skip. Fast tier unchanged (8
passed, 3 skipped without the flag). The tier is relaunched as `pytest_slow_equivalence3`; the
parser-skip run's records are kept aside as `s26/results/pytest_slow_equivalence.PARSER_SKIP_run2.*`
and are not counted.

---

## L103 -- LANE P, C2 RUNG PAIRNET: THE TRIANGLE-UPDATE PAIRNET ON ESM INPUT, THROUGH THE PIPELINE: +0.041 A ON THE BUILT CHAIN (0.30x MDE), +0.008 ON SELECTION (0.03x); THE BEST MAE OF THE LADDER (2.079 vs 2.339) AND THE HIGHEST gam_eff (+0.129 at cos 0.29, +0.441 at cos 0.53 in location space) BUY NOTHING (2026-09-14 00:29, lane P)

Artefacts: `s26/results/p_ladder_pairnet_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pairnet_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pairnet`: exit 0, 587 s, peak RSS 0.301 GB; models the pinned `pairnet_models/fold{0..4}_c64b4_s0.pt` (S19's PairNet, triangle multiplicative update + axial attention, c = 64, 4 blocks, on the deployed ESM inputs; no training this sprint). This closes S7-11's 'tri on ESM input has not been measured' on the pipeline basis: null-to-slightly-worse at every endpoint, 3/5 folds, medians +0.007. It is the sharpest instance of 'better matrix, worse ranking' in the ladder: PairNet cuts the distance MAE by 11% (2.079 against 2.339; S19 L11 measured 2.073) and moves the posterior mean 44% of the way to the truth in location space at cos 0.53 (the highest cosine of any rung; S25 L12's example operator is exactly 'cos 0.5', predicted +0.024 A) and the built chain moves +0.041 +/- 0.049. Joint consistency (S19 L11's closure through the distance-geometry fit, -0.080 [-0.242, +0.082]) is now also closed through the readout the pipeline uses. Strata: FAIL18 -0.10 (10W/8L), other108 +0.06; the 9-10-mers +0.12 on the built chain and -0.24 on selection (n = 20, SE 0.2), unresolved.

```
  pairnet vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2534 (med 3.0884)   b 3.2126 (med 2.9661)   n=126
    effect +0.0407   median +0.0067   SE 0.0486   MDE 0.1361   effect/MDE +0.30
    iid  CI95 [-0.0546, +0.1306]
    fold CI95 [-0.0442, +0.1143]   folds same sign 3/5   per-fold 0:+0.137 1:-0.011 2:+0.005 3:-0.104 4:+0.141
    58W/68L/0T   worst degradation +1.6875 (2N9M)   p90 +0.6843   power 0.13  Type-M 2.93
    concentration: drop-top10 +0.1472 vs uniform-effect null p10/p50/p90 +0.0837/+0.1440/+0.2047 -> pctile 0.527
    VERDICT: NOT MEASURED (|effect| 0.0407 <= its own MDE 0.1361, 0.30x)
  pairnet vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0843 (med 2.9058)   b 3.0483 (med 2.8373)   n=126
    effect +0.0360   median +0.0069   SE 0.0490   MDE 0.1373   effect/MDE +0.26
    iid  CI95 [-0.0606, +0.1295]
    fold CI95 [-0.0439, +0.1059]   folds same sign 4/5   per-fold 0:+0.140 1:+0.021 2:+0.010 3:-0.130 4:+0.110
    59W/67L/0T   worst degradation +1.6203 (2N9M)   p90 +0.6848   power 0.11  Type-M 3.32
    concentration: drop-top10 +0.1455 vs uniform-effect null p10/p50/p90 +0.0843/+0.1430/+0.2003 -> pctile 0.524
    VERDICT: NOT MEASURED (|effect| 0.0360 <= its own MDE 0.1373, 0.26x)
  pairnet vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4617 (med 3.5274)   b 3.4540 (med 3.4779)   n=126
    effect +0.0077   median +0.0067   SE 0.0817   MDE 0.2289   effect/MDE +0.03
    iid  CI95 [-0.1504, +0.1603]
    fold CI95 [-0.2030, +0.1603]   folds same sign 3/5   per-fold 0:+0.247 1:-0.033 2:+0.027 3:-0.410 4:+0.143
    57W/63L/6T   worst degradation +2.4293 (2LWS)   p90 +1.0498   power 0.05  Type-M 24.82
    concentration: drop-top10 +0.1830 vs uniform-effect null p10/p50/p90 +0.0796/+0.1797/+0.2840 -> pctile 0.518
    VERDICT: NOT MEASURED (|effect| 0.0077 <= its own MDE 0.2289, 0.03x)
  strata [arm]: len 9-10 n=20 +0.115 (SE 0.125, med +0.033, 9W/11L) | len 11-12 n=32 +0.023 (SE 0.100, med +0.009, 13W/19L) | len 13-14 n=37 -0.010 (SE 0.075, med +0.016, 17W/20L) | len 15-16 n=37 +0.066 (SE 0.101, med -0.000, 19W/18L) | FAIL18 n=18 -0.096 (SE 0.127, med -0.005, 10W/8L) | other108 n=108 +0.064 (SE 0.052, med +0.010, 48W/60L)
  strata [cloud]: len 9-10 n=20 +0.072 (SE 0.131, med +0.024, 10W/10L) | len 11-12 n=32 +0.036 (SE 0.104, med +0.020, 13W/19L) | len 13-14 n=37 +0.003 (SE 0.076, med +0.011, 17W/20L) | len 15-16 n=37 +0.049 (SE 0.099, med -0.001, 19W/18L) | FAIL18 n=18 -0.142 (SE 0.142, med -0.005, 10W/8L) | other108 n=108 +0.066 (SE 0.052, med +0.012, 49W/59L)
  strata [sel]: len 9-10 n=20 -0.239 (SE 0.205, med -0.190, 11W/8L) | len 11-12 n=32 +0.178 (SE 0.167, med +0.158, 11W/21L) | len 13-14 n=37 +0.072 (SE 0.141, med +0.000, 17W/17L) | len 15-16 n=37 -0.071 (SE 0.156, med +0.000, 18W/17L) | FAIL18 n=18 +0.082 (SE 0.163, med +0.051, 8W/9L) | other108 n=108 -0.005 (SE 0.092, med +0.007, 49W/54L)
  gamma-equivalent: gam_eff prob-space +0.1285 at cos +0.287 ; loc-space +0.4412 at cos +0.533 ; MAE 2.0791 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 3/5 ; verdict (arm): NOT MEASURED (|effect| 0.0407 <= its own MDE 0.1361, 0.30x)
```

---

## L104 -- THE OPT-IN EQUIVALENCE TIER PASSES 3/3 ONCE IT CAN RUN: THE TWO ON-DISK ARMS ARE BIT-IDENTICAL UP TO THE PROJECTION AND THE PROJECTION'S DIVERGENCE STAYS INSIDE ITS PINNED BAND; ONE STAGED-FILE COLLISION BETWEEN LANES (2026-09-14 00:31, lane I)

Job `pytest_slow_equivalence3` (AMBER, `VERIFY_SLOW=1` in the command, tree `7e08b968`):
**3 passed, 0 failed, 40.2 s, peak RSS 0.313 GB**, `s26/results/pytest_slow_equivalence.xml`,
`s26/logs/pytest_slow_equivalence3.log`. What the three assert, now that the L102 parser lets
them run: `test_everything_up_to_the_projection_is_bit_identical` (retrieval, filter membership
and order, the coordinate average `avg_ca` and the scalars `shipped`, `pool_best`, `pool_mean`,
`top_m_best`, `top_m_mean`, `rmsd_avg`, `n_windows` equal with `==` on all 8 smoke8 targets,
legacy arm vs consolidated arm), `test_the_projection_selects_a_different_degenerate_branch`
(upstream exact, and the `rmsd_full` delta between the arms within |mean| < 0.05 and max < 0.25 A),
`test_no_stage_is_skipped_in_the_consolidated_arm` (every stage timed > 0, AMBER lowered the
energy and moved the structure, 75 retained).

Basis of the result, stated plainly: both arms were SERVED FROM THE ON-DISK CACHES
(production `bench_results/cache/1fc9f2dcf489e2fb`, baseline `464a0ddb5f283e04`; record mtimes
2026-09-04, unchanged), which the test's own docstring allows ("uses the harness's own cache"),
so the 40 s is six resumed pipeline runs and the tier certifies the equivalence of the two
recorded arms, the same fact `bench_results/compare_tuning126.json` records as
`science_delta = 0`, not a fresh computation. `verify/run_equiv2.sh` (L19) is the fresh one.

Record: `s26/TEST_RUN.md` / `s26/results/test_run.json` now fold opt-in jobs into the skip count
instead of adding tests: **370 tests, 360 passed, 0 failed, 0 errors, 10 skipped** with the tier in
(357 / 13 without it). The integration tier (8 items) is still queued as `pytest_slow_integration2`.

Housekeeping fact for the hygiene list: at 00:30:53 lane PH's `git commit` swept four files I had
just staged (`s26/i_test_report.py`, `s26/TEST_RUN.md`, `s26/results/test_run.json`,
`pytest_slow_equivalence.xml`) into its commit `a6ce3ab6` because eight lanes share one index;
the content is correct and committed, only the message is PH's. Stage-and-commit in one call.

---

## L105 -- TOURNAMENT ITEM 10, amber_prior_partner (P's H_C3a, orphaned to W): THE LEAVE-FOLD-OUT CHOICE DECLINES TO MIX ON ALL FIVE FOLDS (lam* = 0), SO THE DEPLOYABLE ARM IS THE INCUMBENT BIT-EXACTLY; EVERY ONE OF THE 29 MIXTURE CELLS IS WORSE THAN THE SHIPPED POSTERIOR ON THE POINT CLOUD (+0.010 TO +0.152 A, NONE PAST ITS MDE, 22 OF 29 WITH THE FOLD CI ABOVE ZERO); THE ENERGY WEIGHTING CHANGES NOTHING THE UNWEIGHTED HISTOGRAM DOES NOT; THE PER-TARGET ORACLE OVER 30 CELLS TRANSFERS 8% IN SPLIT HALF (2026-09-14, W)

Pre-registered in `s26/PREREG_amber_prior_partner.md` (written 22:15, before any real-target run;
addendum 1 records the staging refinement before the run). Stage 1 (native-free): job
`w_amberprior_clouds` (exit 0, 135 s, peak 0.163 GB; `s26/results/w_selfcopy_amberprior_clouds.json`,
complete 126/126; gate top-75 == `sub` 126/126; lam = 0 asserted bit-exact against the shipped risk
table and beta = 0 against `p_ladder.pool_histogram` on every target). Stage 2 (gated): job
`w_amberprior_endpoint` (exit 0, 4786 s wall under a six-job load, peak 0.131 GB;
`s26/results/w_selfcopy_amberprior_endpoint.json`, complete 126/126); per-cell cloud contrasts in
`s26/results/w_amberprior_cells_cloud.json`. Probe `w_amberprior_probe` (1A13, 5 s). Operator: per
pair, the K = 500 pool's 17-bin distance histogram weighted by exp(-beta x zrank(E_AMBER)) from
the 63,000 cached ff14SB/GBn2 single points (`s24/cache_amber`, pool identity asserted), mixed
into the shipped posterior at weight lam, through a genuine `core.predict.Distogram` risk table,
the shipped score, the top-75 and the medoid-frame average; grid lam in {0, 0.05, 0.1, 0.2, 0.35,
0.5} x beta in {0, 0.5, 1, 2, 4}; (lam*, beta*) chosen leave-fold-out on the point cloud (the
declared cost fork), verdict on the built chain. Basis on every line; negative = better.

**The leave-fold-out choice is (0, 0) on every fold** (`chosen_cells_by_fold`: 0,0 for folds 0 to
4; the beta = 0 row's own choice likewise 0,0). The deployable arm is therefore the incumbent
bit-exactly on 126/126 targets, and every registered contrast (chosen minus shipped; chosen minus
the unweighted mixture; chosen minus the rank-permuted-AMBER mixture) is an exact tie on both
bases (`ST.fmt`: effect +0.0000, 0W/0L/126T, "NOT MEASURED (|effect| 0.0000 <= its own MDE
0.0000)"). H_AP is refuted in the strongest form the design allows: the choice the rule makes with
four training folds is to leave the prior alone.

**Why, cell by cell (point cloud, all folds, `ST.compare` against the shipped posterior; lam = 0 is
the identity at every beta):**

    lam    beta 0                  beta 0.5                beta 1                  beta 2                  beta 4
    0.05   +0.016 [+0.002,+0.030]  +0.015 [+0.003,+0.024]  +0.019 [+0.005,+0.033]  +0.012 [+0.001,+0.030]  +0.010 [+0.001,+0.020]
    0.1    +0.021 [-0.008,+0.045]  +0.022 [-0.008,+0.053]  +0.021 [-0.010,+0.053]  +0.021 [-0.003,+0.047]  +0.009 [-0.003,+0.024]
    0.2    +0.041 [-0.013,+0.086]  +0.035 [-0.012,+0.086]  +0.029 [-0.013,+0.086]  +0.033 [+0.002,+0.069]  +0.042 [+0.011,+0.078]
    0.35   +0.080 [-0.005,+0.151]  +0.059 [-0.008,+0.131]  +0.057 [-0.007,+0.136]  +0.060 [+0.009,+0.125]  +0.075 [+0.021,+0.127]
    0.5    +0.145 [+0.001,+0.260]  +0.143 [+0.007,+0.260]  +0.125 [+0.031,+0.235]  +0.125 [+0.032,+0.228]  +0.152 [+0.061,+0.243]
    (effect and fold CI95, A; MDEs 0.036 to 0.186; effect/MDE 0.15 to 1.00; W/L from 60/64 to 49/77; every cell NOT MEASURED by the rule, 22 of 29 with the fold CI above zero)

Every mixture cell is worse than the shipped posterior on the point cloud, monotonically in lam
(the typicality direction of S7-3 / S19 L14, reproduced: moving the prior toward the pool's own
histogram costs accuracy), and beta moves a cell by at most 0.02 A in either direction with no
consistent sign (at lam 0.05 the energy weighting helps by 0.006; at lam 0.5 beta 4 hurts by
0.007): AMBER's ordering carries nothing through the prior that the unweighted histogram does not,
which is S24 L16's rank-permuted null (-0.0010, 0.08x MDE) seen from the prior side. The
per-target ORACLE over the 30 cells is -0.290 A on the cloud (2.804 vs 3.048), and
`ST.best_of_k_within` accounts for 224% of it with its valid across-target null (k_eff 9.7), the
split-half transfer is -0.022 (8% of the oracle): "NOT A SIGNAL", as the PREREG required it to be
quoted. Registered expectation (0.00 +/- 0.03, lam* = 0 on most folds): held, with lam* = 0 on all
five. Power: for the deployable arm the question is moot (it is the identity); for the cells, a
gain of 0.04 A or more at lam 0.05 and of 0.19 A at lam 0.5 is excluded on the cloud. The built
chain was projected only for the chosen cell (the identity), the beta = 0 chosen cell (the
identity) and the permuted control (the identity when lam = 0), per the declared staging; no
built-chain number exists for a non-chosen cell, and none is claimed. Replication: the choice is
deterministic given the artefact; the permutation seeds were never consumed (lam* = 0). The last
AMBER form in the record, as a distribution inside the prior, is closed on this instrument.
Deviations from the PREREG: none beyond addendum 1's staging.

## L106 -- B3: THE SET WHERE THE PIPELINE BEATS SEQUENCE-ONLY IS NOT CHARACTERISABLE NATIVE-FREE BY THE PRE-REGISTERED STANDARD; THE SIGN CLASSIFIER IS AT THE PERMUTATION NULL AGAINST BOTH COMPARATORS, AND ONLY THE SIZE OF THE GAIN OVER THE CONSTANT HELIX IS PARTLY PREDICTABLE (R2 0.40), WHICH IS THE HELICITY SIGNAL (2026-09-14 00:55, lane P)

`s26/p_b3.py run`, `s26/results/p_b3.json` (job `p_b3_run`, exit 0, 100 s, peak RSS 0.053 GB),
pre-registered in `s26/PREREG_B3.md`; 45 native-free features (`s26/results/p_b3_features.json`:
length, composition, the top-75 members' H/E/C content, distogram entropy, ESM contact-map and
retrieval-score statistics), nested leave-fold-out ridge (`s26/p_stats.py`), 300-draw label
permutation null. Built chains on both sides: pipeline `rmsd_arm` 3.2148, sequence-only torsion
predictor 3.7705 (`s13/cache/tors_rows.npz`), constant helix 4.0648 (`s14/results/ladder.json`).

    sign(arm - tors): held-out balanced accuracy 0.522  vs permutation null mean 0.498, 95th pct 0.578  (p = 0.263)
    sign(arm - helix): held-out balanced accuracy 0.557  vs null mean 0.498, 95th pct 0.566  (p = 0.093)
    d = arm - tors, ridge regression: held-out R2 +0.244; squared-error reduction vs the training-fold mean
      -0.5366  SE 0.2309  MDE 0.6469  (-0.83x)  fold CI [-0.882, -0.185]  4/5  84W/42L  NOT MEASURED
    d = arm - helix, ridge regression: held-out R2 +0.404; squared-error reduction
      -1.1420  SE 0.3415  MDE 0.9566  (-1.19x)  fold CI [-1.495, -0.700]  5/5  94W/32L  BETTER (Type-M zone)

The falsifier required BOTH the classifier above its permutation null AND the regression's
fold-clustered MSE reduction beyond MDE. Against the torsion predictor both halves fail
(classifier at the null; regression 0.83x MDE). Against the constant helix the classifier fails
(0.557 against a 95th percentile of 0.566) and the regression passes at 1.19x MDE: the SIZE of
the pipeline's gain over a helix is partly predictable, the SIGN is not. The features that carry
the size are the top-75 members' helix content and the length (descriptive ridge weights in
`s26/agentP_FINDINGS.md` section 10): a helical retrieval pool means the constant helix is
already close and the gain is small. That is S17 L25's "the one target-level signal that works,
and it is the trivial one", re-found. ORACLE stratum for the record: on FAIL18 the torsion
predictor beats the pipeline by 0.463 A and on the other 108 the pipeline wins by 0.725 A
(S12's failure-class split), and nothing native-free in this feature set locates that stratum.
**B3 verdict: the falsifier fires; Proposal B's B3 route does not yield a native-free
characterisable set. Per the campaign rule, Proposal B's verdict is REPLACE.**

---

## L107 -- CORRECTION TO L106: THE FEATURE THAT CARRIES THE SIZE OF THE GAIN IS THE POOL'S STRAND CONTENT, NOT ITS HELIX CONTENT OR THE LENGTH (2026-09-14 00:55, lane P)

Descriptive ridge weights and correlations, full data (`s26/agentP_FINDINGS.md` section 10):
for d = arm - helix the largest standardised weight is the top-75 members' E (strand) fraction
(-0.612; rho(ss_E, d) = -0.638), then the distogram's mean sd (+0.394), its max entropy (-0.273),
the cysteine fraction (-0.266), its multimodal fraction (+0.253) and the ESM long-range contact
mass (+0.230; rho -0.42); rho(ss_H, d) = +0.463 and rho(n, d) = -0.128. For d = arm - tors the
same ordering (ss_E -0.636, rho -0.601). So: when the retrieval pool is strand-like the
pipeline's gain over a constant helix is LARGE, and when it is helical the gain is small; the
sign stays unpredictable (L106). "Helix content and the length" in L106 is struck; "the pool's
strand content" is the sentence.

---

## L108 -- ADDENDUM TO L44 / L58: PART B COMPLETE ON THE FOUR CARRIER-OUT MODELS. REMOVING THE CARRIER FROM THE FOLD MODEL'S TRAINING LABELS CHANGES THE BUILT CHAIN BY +0.011 / -0.002 / -0.246 / -0.027 A ON 1CEK / 2FBU / 2P5H / 6B9K: F3 IS FALSIFIED BY 2P5H, IN THE HARMFUL DIRECTION (THE LEAK MAKES THE LEAKED TARGET WORSE); THE DIRECT DEV-PROXY BOUND RISES FROM 0.002 TO 0.008 A, STILL IMMATERIAL; THE CLASS STAYS MINOR BY THE ENVELOPE (2026-09-14, W)

Part B of `s26/PREREG_selfcopy_bound.md` is now measured on all four carrier-out models: 1A11 out
of fold 2 (built before the host kill, L40), 2LMF out of fold 4 (job `w_train_out_2LMF`, 1049 s,
peak 1.244 GB), 2P5J out of fold 4 (`w_train_out_2P5J`, 1063 s, 1.250 GB), 1U6V out of fold 0
(`w_train_out_1U6V`, 893 s, 1.248 GB), each through lane P's exact pca32 training path with the
one chain's pairs removed (276 / 231 / 120 / 120 pairs fewer); the reference for every fold is
lane P's `pca32_fold<f>_s0.pt` (identical function, corpus and seed, nothing removed), which
reproduces the pinned emission at 0.000 on all four targets. Gated re-run: job `w_endpoint_report2`
(posterior, endpoint, report; exit 0, 145 s, peak 0.309 GB); the L55 additions re-applied by
`s26/w_bound_addendum.py`; artefacts `s26/results/w_selfcopy_posterior.json`,
`w_selfcopy_endpoint.json`, `w_selfcopy_bound.json` (regenerated, provenance stamped, the
original kept under `provenance_original`). Parts A, C, D and E are unchanged. Basis on every
line; delta = reference minus carrier-out: positive = the carrier's presence HELPED the target.

    target  carrier  fold   posterior change (mean |dE[d]| per pair; top-75 overlap)   delta arm    cloud     fit       sel       paired gain
    1CEK    1A11     2      1.024 A; 0.76                                              +0.0114     +0.0204   +0.0115   -0.0062   +0.0177
    2FBU    2LMF     4      0.791 A; 0.89                                              -0.0021     -0.0002   -0.0021   -0.0040   +0.0019
    2P5H    2P5J     4      1.382 A; 0.77                                              -0.2460     -0.1634   -0.2310   +0.0000   -0.2460
    6B9K    1U6V     0      1.373 A; 0.91                                              -0.0265     -0.0318   -0.0188   +0.0000   -0.0265

    both channels removed (self-window dropped AND carrier out of the model), production minus clean:
    1CEK +0.0114   2FBU -0.0021   2P5H -0.2460   6B9K -0.0243   (arm)

**F3 is falsified**: |delta arm| on 2P5H is 0.246 A, above the registered 0.10 A line; its
control clause (carrier-out change against two control-out chains per fold) is measured so far on
one control (9BAF out of fold 0, `w_train_out_9BAF_f0`, 1440 s, 1.247 GB) and the other five are
queued; the clause is reported when they exist. The direction is the one L52 makes expected:
2P5J's copy of 2P5H's sequence sits 2.33 A from 2P5H's native, and training the fold model on
that carrier pulls the posterior toward the carrier's geometry (mean |dE[d]| 1.38 A per pair, the
largest of the four) and the built chain 0.25 A AWAY from the native. On 1CEK, whose carrier
segment is 0.60 A from the native, the carrier helps by 0.011; on 6B9K (carrier 4.13 A away) it
hurts by 0.027; on 2FBU (3.28 A away) it is 0.002. So the training channel's sign follows the
cross-deposit distance of the copy, which is what a memorising model does with a copy that is
not the native: the leak is not a gift to the leaked target, it is a bias toward another
deposit's conformation.

**The bound, updated (Part D; `w_selfcopy_bound.json :: signed_bounds_gated/both_removed4`, now
n = 4):** (2/60) x max |delta| with both channels removed = 0.0082 A on the built chain (2P5H),
0.0082 on the paired gain, 0.0054 on the cloud, 0.0002 on selection: IMMATERIAL (below 0.017) on
every basis, up from 0.0004 / 0.0006 at n = 1. The envelope rows are unchanged (mean-CI 0.028
built chain / 0.048 selection / 0.023 gain; worst target 0.151 / 0.194 / 0.123 under A2, L58), so
the pre-registered class stays MINOR by the envelope clause and IMMATERIAL by every direct
measurement, now with all four dev self-copies measured on both channels. The sign matters for
the benchmark verdict: where the direct measurement is not zero it is HARMFUL to the leaked
target's built chain and neutral to its argmin, so the un-leaked benchmark paired gain would, if
anything, be slightly more favourable to the architecture than the +0.0103 reported; by at most
(2/60) x 0.246 = 0.008 A under A2. The caveat text (L55) stands with "dev-proxy price 0.008 A"
in place of 0.002. Replication: not applicable to deterministic retrains at seed 0 (a second seed
of the 2P5H retrain is the natural check of the one large value and is queued behind the
control-out models if time allows). Deviation from the PREREG: Part B's controls are 1 of 6 at
this entry; F3's control clause is open.

## L109 -- window_provenance, H_P3 (THE ACHIEVABLE READOUT): DOWN-WEIGHTING FRAGMENT WINDOWS IN THE TOP-75 AVERAGE, THE WEIGHT CHOSEN LEAVE-FOLD-OUT, IS +0.0066 A ON THE BUILT CHAIN (0.15x MDE) AGAINST THE UNIFORM READOUT AND -0.018 (0.38x) AGAINST ITS PERMUTED-WEIGHT CONTROL; NOT MEASURED; A GAIN OF 0.044 A OR MORE IS EXCLUDED; THE 0.13 A ORACLE CLASS GAP OF L85 DOES NOT REACH THE EMISSION (2026-09-14, W)

Pre-registered in `s26/PREREG_window_provenance.md` (section 1, H_P3; addendum 1 declared the one
reduction, 4 permutation draws instead of 8, before the run). Job `w_provenance_readout` (exit 0,
9321 s wall under the six-job load, peak RSS 0.110 GB; `s26/results/w_selfcopy_provenance_readout.json`,
complete 126/126; the uniform-weight cloud asserted equal to the production `avg_ca` on every
target). Operator: the shipped top-75 (production `sub`) averaged in the production medoid frame
with fragment-class members at relative weight w in {0, 0.5, 1, 2} and peptide-derived members
at 1, one projection per weight; w chosen on the four training folds' built-chain mean and
applied to the fifth; control: the same weights permuted across the 75 members (4 seeded draws,
`s15.seed.stable_rng`) at the chosen weight. Basis on every line; negative = the routed arm is
better.

The leave-fold-out choice: w = 0.5 on folds 0, 3, 4 and w = 0 on folds 1, 2 (the routed arm
down-weights fragments everywhere; it never chooses the uniform readout and never up-weights).

  routed_minus_uniform [arm], fragment weight chosen leave-fold-out on the built chain
    a 3.2193 (med 2.9905)   b 3.2126 (med 2.9661)   n=126
    effect +0.0066   median -0.0017   SE 0.0157   MDE 0.0440   effect/MDE +0.15
    iid  CI95 [-0.0238, +0.0370]
    fold CI95 [-0.0155, +0.0289]   folds same sign 2/5   per-fold 0:-0.002 1:-0.023 2:+0.046 3:-0.016 4:+0.020
    68W/58L/0T   worst degradation +0.7206 (9UV5)   p90 +0.1643   power 0.07  Type-M 5.65
    concentration: drop-top10 +0.0345 vs uniform-effect null p10/p50/p90 +0.0154/+0.0335/+0.0535 -> pctile 0.535
    VERDICT: NOT MEASURED (|effect| 0.0066 <= its own MDE 0.0440, 0.15x)
  routed_minus_permuted [arm], fragment weight chosen leave-fold-out on the built chain
    a 3.2193 (med 2.9905)   b 3.2372 (med 3.0123)   n=126
    effect -0.0179   median -0.0050   SE 0.0168   MDE 0.0470   effect/MDE -0.38
    iid  CI95 [-0.0503, +0.0147]
    fold CI95 [-0.0316, -0.0057]   folds same sign 4/5   per-fold 0:+0.004 1:-0.039 2:-0.032 3:-0.017 4:-0.010
    70W/56L/0T   worst degradation +1.0627 (3SGO)   p90 +0.1302   power 0.19  Type-M 2.35
    concentration: drop-top10 +0.0150 vs uniform-effect null p10/p50/p90 -0.0054/+0.0138/+0.0340 -> pctile 0.525
    VERDICT: NOT MEASURED (|effect| 0.0179 <= its own MDE 0.0470, 0.38x)
  routed_minus_uniform [cloud], fragment weight chosen leave-fold-out on the built chain
    a 3.0484 (med 2.8024)   b 3.0483 (med 2.8373)   n=126
    effect +0.0001   median -0.0042   SE 0.0150   MDE 0.0421   effect/MDE +0.00
    iid  CI95 [-0.0270, +0.0289]
    fold CI95 [-0.0274, +0.0450]   folds same sign 1/5   per-fold 0:-0.019 1:-0.026 2:+0.089 3:-0.034 4:-0.013
    69W/57L/0T   worst degradation +0.9086 (9UV5)   p90 +0.0997   power 0.05  Type-M 495.60
    concentration: drop-top10 +0.0252 vs uniform-effect null p10/p50/p90 +0.0062/+0.0240/+0.0436 -> pctile 0.527
    VERDICT: NOT MEASURED (|effect| 0.0001 <= its own MDE 0.0421, 0.00x)
  routed_minus_permuted [cloud], fragment weight chosen leave-fold-out on the built chain
    a 3.0484 (med 2.8024)   b 3.0768 (med 2.8418)   n=126
    effect -0.0284   median -0.0130   SE 0.0150   MDE 0.0420   effect/MDE -0.68
    iid  CI95 [-0.0565, +0.0029]
    fold CI95 [-0.0480, -0.0167]   folds same sign 5/5   per-fold 0:-0.017 1:-0.065 2:-0.016 3:-0.033 4:-0.017
    85W/41L/0T   worst degradation +1.0753 (3SGO)   p90 +0.0692   power 0.48  Type-M 1.44
    concentration: drop-top10 +0.0010 vs uniform-effect null p10/p50/p90 -0.0171/-0.0002/+0.0188 -> pctile 0.537
    VERDICT: NOT MEASURED (|effect| 0.0284 <= its own MDE 0.0420, 0.68x)

Fixed weights over all folds (diagnostic, NOT the routed arm; built chain minus uniform): w = 0
-0.0043 (0.06x MDE 0.076; fold CI [-0.042, +0.030]); w = 0.5 -0.0117 (0.39x MDE 0.030; fold CI
[-0.041, +0.010], 4/5 folds); w = 2 +0.0183 (0.68x MDE 0.027; fold CI [-0.008, +0.040]). The
direction is consistent with L85 (down-weighting fragments helps slightly, up-weighting hurts
slightly), and none of it is measurable at n = 126.

Verdict: H_P3 refuted; the idea closes with its power statement. The routed readout is +0.007 A
against uniform (0.15x MDE; it excludes a gain of 0.044 A or more) and -0.018 against the
permuted-weight control (0.38x; the fold CI excludes zero at 4/5 folds and 0.38x MDE, i.e. the
class information is worth something against noise-with-the-same-weights but nothing against
doing nothing). Mechanism: the 0.13 A single-window class gap inside the top-75 (L85) acts on 8%
of the members of a 75-member mean, and the L44 / L64 controls already showed that a few
members' change moves the chain mostly orthogonally to the native; the registered expectation
(0.00 to -0.02 against an MDE of about 0.05) held. What the census leaves for the report: the
pool is three-quarters fragment windows; a peptide-derived window is 0.42 A nearer the native
before the score and 0.09 to 0.13 A after it; and re-weighting by provenance cannot convert that
into an emission gain at this n. Replication: for a positive, the second permutation seed and the
reversed fold order were pre-declared; there is no positive to replicate. Deviations from the
PREREG: the 4-draw control (declared in addendum 1); the `fallback` flag (a target with no
peptide-derived member and w = 0 falls back to uniform) fired on 0 targets at the chosen weights
(w = 0.5 or 0 with at least one peptide member on every target: 113 have a terminal, 30 a whole,
every target an interior or terminal member by the census).

## L110 -- C4 FIRST BLOCK (ALL 25 NEW FEATURES): THE SIXTH AND SEVENTH m* ROUTERS AND THE SECOND s* ROUTER ARE NULL-TO-HARMFUL; NEW FEATURE PROVENANCE DOES NOT MATTER, AS S23 L7 SAID (2026-09-14 01:27, lane P)

`s26/p_c4.py run` (job `p_c4_run`, running; `s26/results/p_c4.json` is rewritten after every
block), pre-registered in `s26/PREREG_C4.md`. Anchor: the rebuilt m = 75 cloud reproduces
`s12/results/agg_surface.json` and `s23/results/errdecomp.json` per target at 0.00e+00. Nested
leave-fold-out ridge (`s26/p_stats.py`), 200-draw label permutation null for m, 100 for s.
Features: pax_* (principal-axis spread of the top-75 cloud), sim_* (retrieval-score entropy and
gaps), dg_* (posterior bin entropy), con_* (ESM contact-map statistics), 25 in all, jointly.

```
  C4 m-router A [new_all] routed - fixed m=75 (POINT CLOUD, persisted surface)
    a 3.1168 (med 2.9363)   b 3.0483 (med 2.8373)   n=126
    effect +0.0684   median +0.0000   SE 0.0269   MDE 0.0754   effect/MDE +0.91
    iid  CI95 [+0.0169, +0.1236]
    fold CI95 [+0.0023, +0.1592]   folds same sign 3/5   per-fold 0:+0.105 1:+0.238 2:-0.004 3:+0.057 4:-0.022
    40W/61L/25T   worst degradation +1.6752 (7N2I)   p90 +0.3081   power 0.72  Type-M 1.18
    concentration: drop-top10 +0.1082 vs uniform-effect null p10/p50/p90 +0.0745/+0.1061/+0.1412 -> pctile 0.533
    VERDICT: NOT MEASURED (|effect| 0.0684 <= its own MDE 0.0754, 0.91x)
    perm null mean +0.0265 p05 +0.0003  p_perm 0.960
  C4 m-router A [new_all] routed - fixed m=75 (BUILT CHAIN, same projection both sides)
    a 3.2834 (med 3.1928)   b 3.2126 (med 2.9661)   n=126
    effect +0.0707   median +0.0000   SE 0.0303   MDE 0.0849   effect/MDE +0.83
    iid  CI95 [+0.0133, +0.1340]
    fold CI95 [+0.0141, +0.1418]   folds same sign 4/5   per-fold 0:+0.091 1:+0.216 2:-0.006 3:+0.063 4:+0.012
    40W/61L/25T   worst degradation +2.1542 (7N2I)   p90 +0.4232   power 0.65  Type-M 1.25
    concentration: drop-top10 +0.1164 vs uniform-effect null p10/p50/p90 +0.0790/+0.1138/+0.1556 -> pctile 0.534
    VERDICT: NOT MEASURED (|effect| 0.0707 <= its own MDE 0.0849, 0.83x)
  C4 m-router B [new_all] routed - fixed m=75 (POINT CLOUD, persisted surface)
    a 3.0751 (med 2.9202)   b 3.0483 (med 2.8373)   n=126
    effect +0.0267   median +0.0035   SE 0.0169   MDE 0.0474   effect/MDE +0.56
    iid  CI95 [-0.0064, +0.0594]
    fold CI95 [+0.0120, +0.0457]   folds same sign 5/5   per-fold 0:+0.016 1:+0.036 2:+0.019 3:+0.065 4:+0.005
    44W/64L/18T   worst degradation +0.7631 (8HVS)   p90 +0.2352   power 0.35  Type-M 1.67
    concentration: drop-top10 +0.0623 vs uniform-effect null p10/p50/p90 +0.0419/+0.0611/+0.0805 -> pctile 0.530
    VERDICT: NOT MEASURED (|effect| 0.0267 <= its own MDE 0.0474, 0.56x)
    perm null mean +0.0192 p05 +0.0065  p_perm 0.865
  C4 m-router B [new_all] routed - fixed m=75 (BUILT CHAIN, same projection both sides)
    a 3.2364 (med 3.0593)   b 3.2126 (med 2.9661)   n=126
    effect +0.0237   median +0.0000   SE 0.0193   MDE 0.0542   effect/MDE +0.44
    iid  CI95 [-0.0137, +0.0606]
    fold CI95 [+0.0045, +0.0529]   folds same sign 5/5   per-fold 0:+0.001 1:+0.002 2:+0.025 3:+0.082 4:+0.013
    48W/60L/18T   worst degradation +0.8722 (8HVS)   p90 +0.2487   power 0.23  Type-M 2.07
    concentration: drop-top10 +0.0616 vs uniform-effect null p10/p50/p90 +0.0391/+0.0607/+0.0841 -> pctile 0.520
    VERDICT: NOT MEASURED (|effect| 0.0237 <= its own MDE 0.0542, 0.44x)
  C4 s-router [new_all] routed - s=1 (point cloud)
    a 3.0614 (med 2.8817)   b 3.0483 (med 2.8373)   n=126
    effect +0.0130   median +0.0138   SE 0.0147   MDE 0.0412   effect/MDE +0.32
    iid  CI95 [-0.0150, +0.0409]
    fold CI95 [-0.0161, +0.0529]   folds same sign 3/5   per-fold 0:+0.028 1:+0.085 2:-0.007 3:+0.005 4:-0.031
    58W/68L/0T   worst degradation +0.6981 (5V5B)   p90 +0.2071   power 0.14  Type-M 2.78
    concentration: drop-top10 +0.0385 vs uniform-effect null p10/p50/p90 +0.0196/+0.0380/+0.0576 -> pctile 0.517
    VERDICT: NOT MEASURED (|effect| 0.0130 <= its own MDE 0.0412, 0.32x)
    rho(pred, s*) -0.203 ; perm null mean +0.0089
```

Reading. The m* oracle over the 15-rung ladder is real (S22 L4: -0.239 split-half) and the
routers point the HARMFUL way: router A (multi-output ridge, argmin of the prediction) +0.068 A
on the point cloud at 0.91x its MDE with a median of 0.000 and worse than 96% of the permutation
null (+0.027), and +0.071 on the built chain through the same projection (0.83x); router B
(ridge on log m*, nearest rung) +0.027 (0.56x), null +0.019, p_perm 0.865. The s* router
predicts s* with rho -0.203 (the wrong sign) and costs +0.013 (0.32x MDE) against s = 1.
Neither m router clears its MDE and none is significantly harmful in the fold CI, so
PREREG_C4's stop rule (a fold CI above zero) is not triggered and the remaining blocks run;
their entries follow as the job writes them. The finding is the one the record predicted (S22
L7/L10, S23 L3/L6/L7): a per-target quantity that transfers split-half is invisible to a
linear router on any native-free feature set at n = 126, whatever the features describe.

---
## L111 -- LANE I STARVED FOR 75 MINUTES BY THE OTHER LANES' RELAUNCH CADENCE; CAP RAISED TO SEVEN AND NON-CRITICAL LAUNCHES HELD SO THE INTEGRATION TIER, THE FROZEN REBUILD AND THE VERIFY CHAIN CAN RUN (2026-09-14 01:33, coordinator)

jobrun has no priority: a lane that relaunches within seconds of a job finishing wins the
freed slot, and lane I's three waiters (pytest_slow_integration2, resultslab_rebuild, the
20-audit verify chain gated behind them) lost every race since 00:17. Ruling: the file cap is
raised from six to seven; lane W holds its Part B retrains after the one now training and
launches nothing new until lane I's three items are registered; lane Q holds the A1
replication until then (A3 keeps running); lane PH has nothing new to launch; lane P keeps
its critical path (raw folds, C4, C5, the verdicts). Memory still binds: the integration tier
needs 2.3 GB free and the rebuild 1.7 GB under jobrun's rule, so lane W's retrains (1.3 GB)
stay held until the rebuild has run. If the three items cannot run by 03:30, lane I closes
with them recorded as not run and the reason.

---


## L112 -- THE BEST C2 RUNG IS THE SHIPPED PRIOR; THE C3 STAGE-2 DELIVERY FILE CARRIES THE PRODUCTION EMISSION (`s26/results/p_best_rung_chains.json`, 126/126, complete) (2026-09-14 01:37, lane P)

Decision, with the fold-clustered evidence: no rung of the C2 ladder beats the shipped
posterior on the built chain. noesm +0.208 [fold +0.099, +0.394] 5/5 WORSE (L62); esm8m +0.242
[+0.113, +0.385] 5/5 WORSE (L99); conly +0.122 [-0.035, +0.267] (L63); pca128 +0.076 [-0.069,
+0.216] (L72); pairnet +0.041 [-0.044, +0.114] (L103); wide +0.032 [-0.030, +0.128] (L66);
pca32f +0.018 [-0.073, +0.097] (L67); pca32 and mix are the identity (L65, L93); raw pending.
The best rung is therefore the shipped one, and per the coordinator's rule the delivery file
carries it: rows keyed by pdb with `phi`, `psi` (RADIANS, the production `core.project.lam_path`
emission at lam = 0.3, from `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`), `ca` (n x 3),
`rmsd_arm` (mean 3.2148), `rmsd_cloud`, `rmsd_fit`; provenance-stamped, `complete: true`,
n = 126. The production emission is delivered rather than this lane's re-projection because it
is the chain PH's C3 stage 1 already relaxed (L39), so stage 2 IS the stage-1 replication PH
ran; the rebuild-basis re-projection (`p_deliver_shipped`, queued 47 min behind the job cap) is
recorded when it lands as a cross-check only.

---

## L113 -- THE OPT-IN TIER IS COMPLETE: 11/11 PASS. WITH IT THE SUITE READS 370 TESTS, 368 PASSED, 2 SKIPPED (ABSENT ARTEFACTS), 0 FAILED (2026-09-14 01:37, lane I)

Job `pytest_slow_integration2` (AMBER, `VERIFY_SLOW=1` in the command, tree `7e08b968`,
registered 01:33 after 4,620 s in the queue, the 7th slot of L111): **8 passed, 0 failed,
165.6 s, peak RSS 1.139 GB**, `s26/results/pytest_slow_integration.xml`,
`s26/logs/pytest_slow_integration2.log`. The eight items: the 11 Legacy terms on a real
retrieval pool, core vs root, column for column; the AMBER System's charge / sigma / epsilon and
GBn2 per-particle parameters against a System built here from the same ff14SB + GBn2 XML on the
same topology; the pinned 1A13 native interaction energy (-489.9138948277905) bit-exact;
single-point invariance under a rigid translation (relative < 1e-5); NaN-poisoning of the native
through `run_target` with stage 4 on 1CS9 and 1CB3 (every deployable quantity bit-identical,
no NaN reaches any record); and bit-identity of the same target across processes and under a
4-thread ambient environment (child processes through `verify/determinism_audit.run_child`).

With L104's equivalence tier (3/3) that is the whole opt-in tier: **11 run, 11 passed**. Folded
into the suite record (`s26/TEST_RUN.md`, `s26/results/test_run.json`, commit `7ad4ef68`):
**370 unique tests, 368 passed, 0 failed, 0 errors, 2 skipped**, the two skips being the absent
artefacts `bench_results/optimised_tuning126_w6.json` and "no smoke result on disk"
(`tests/test_pipeline.py`). Peak RSS of every job is in the record; the memory guard fired in none.

---

## L114 -- NOTE TO L112 FOR LANE PH: THE DELIVERED TORSIONS ARE RADIANS BUT UNWRAPPED (2026-09-14 01:37, lane P)

`s26/results/p_best_rung_chains.json` carries the production `phi`/`psi` exactly as
`core.project.lam_path` emitted them: radians, not reduced modulo 2*pi (|value| up to 72.7
occurs; 1A13 psi reaches 6.19). Wrap with ((x + pi) mod 2*pi) - pi before any use that
assumes (-pi, pi]; `ca` is unaffected. The file's `torsion_units` field now says so.

---

## L115 -- C4 COMPLETE: TWELVE m* ROUTERS AND SIX s* ROUTERS ON SIX FEATURE BLOCKS, NONE CLEARS ITS MDE, ALL BUT ONE POINT THE HARMFUL WAY; THE ONE THAT DOES NOT (RETRIEVAL-SCORE ENTROPY, ROUTER B, -0.001 A) IS 0.02x ITS MDE (2026-09-14 01:45, lane P)

`s26/results/p_c4.json` (job `p_c4_run`, exit 0, 2,839 s, peak RSS 0.114 GB), pre-registered
in `s26/PREREG_C4.md`. Anchor: rebuilt m = 75 cloud vs `agg_surface` and `errdecomp` 0.00e+00;
rebuilt m = 75 built chain 3.2126 (the L57 basis). Labels: per-m point-cloud RMSD from
`s12/results/agg_surface.json` (15 rungs), s* from `s23/results/errdecomp.json`. The m* oracle
over the grid: -0.408 A, `ST.best_of_k_within` share 1.43 (an order statistic; split-half
-0.099, k_eff 12.7), consistent with S24 L15/L15-A for grid oracles (the real m* signal is the
S22 L4 split-half transfer, which this grid does not measure).

m routers (nested leave-fold-out ridge; A = multi-output over the 15 rungs then argmin, B = ridge
on log m* then nearest rung; 200-draw label permutation null; positive = worse than fixed m = 75):
    new_all    router A  point cloud +0.0684  SE 0.0269  MDE 0.0754 (0.91x)  fold[+0.002,+0.159]  3/5  40W/61L  perm-null mean +0.0265  p_perm 0.960  | built chain +0.0707 (0.83x)
    new_all    router B  point cloud +0.0267  SE 0.0169  MDE 0.0474 (0.56x)  fold[+0.012,+0.046]  5/5  44W/64L  perm-null mean +0.0192  p_perm 0.865  | built chain +0.0237 (0.44x)
    old_S22    router A  point cloud +0.0244  SE 0.0261  MDE 0.0732 (0.33x)  fold[-0.010,+0.057]  4/5  51W/56L  perm-null mean +0.0212  p_perm 0.675  | built chain +0.0466 (0.55x)
    old_S22    router B  point cloud +0.0206  SE 0.0177  MDE 0.0495 (0.42x)  fold[+0.009,+0.034]  5/5  45W/73L  perm-null mean +0.0198  p_perm 0.595  | built chain +0.0198 (0.34x)
    new_pax    router A  point cloud +0.0594  SE 0.0281  MDE 0.0786 (0.76x)  fold[-0.001,+0.138]  4/5  38W/52L  perm-null mean +0.0235  p_perm 0.945
    new_pax    router B  point cloud +0.0198  SE 0.0194  MDE 0.0543 (0.36x)  fold[+0.011,+0.031]  5/5  44W/68L  perm-null mean +0.0189  p_perm 0.595
    new_sim    router A  point cloud +0.0057  SE 0.0251  MDE 0.0703 (0.08x)  fold[-0.020,+0.033]  2/5  43W/53L  perm-null mean +0.0228  p_perm 0.075
    new_sim    router B  point cloud -0.0012  SE 0.0176  MDE 0.0493 (-0.02x)  fold[-0.023,+0.023]  2/5  46W/57L  perm-null mean +0.0188  p_perm 0.010
    new_dgent  router A  point cloud +0.0187  SE 0.0314  MDE 0.0878 (0.21x)  fold[-0.065,+0.124]  3/5  51W/59L  perm-null mean +0.0212  p_perm 0.520
    new_dgent  router B  point cloud +0.0180  SE 0.0186  MDE 0.0522 (0.34x)  fold[-0.022,+0.062]  3/5  42W/53L  perm-null mean +0.0184  p_perm 0.375
    new_con    router A  point cloud +0.0251  SE 0.0225  MDE 0.0632 (0.40x)  fold[-0.014,+0.068]  3/5  35W/49L  perm-null mean +0.0215  p_perm 0.680
    new_con    router B  point cloud +0.0028  SE 0.0158  MDE 0.0442 (0.06x)  fold[-0.022,+0.031]  2/5  33W/50L  perm-null mean +0.0189  p_perm 0.045

s routers (ridge on s*, applied as a global scale about the cloud's centroid; 100-draw null):
    new_all    routed s - s=1  +0.0130  SE 0.0147  MDE 0.0412 (0.32x)  fold[-0.016,+0.053]  3/5  rho(pred, s*) -0.203  perm-null +0.0089
    old_S22    routed s - s=1  +0.0263  SE 0.0169  MDE 0.0474 (0.56x)  fold[-0.016,+0.076]  3/5  rho(pred, s*) -0.199  perm-null +0.0094
    new_pax    routed s - s=1  +0.0066  SE 0.0137  MDE 0.0383 (0.17x)  fold[-0.018,+0.032]  3/5  rho(pred, s*) -0.180  perm-null +0.0119
    new_sim    routed s - s=1  +0.0111  SE 0.0130  MDE 0.0364 (0.30x)  fold[-0.013,+0.042]  3/5  rho(pred, s*) -0.209  perm-null +0.0072
    new_dgent  routed s - s=1  +0.0228  SE 0.0159  MDE 0.0445 (0.51x)  fold[-0.015,+0.065]  3/5  rho(pred, s*) -0.170  perm-null +0.0050
    new_con    routed s - s=1  +0.0270  SE 0.0155  MDE 0.0435 (0.62x)  fold[-0.014,+0.073]  3/5  rho(pred, s*) -0.096  perm-null +0.0074

Reading. Twelve m routers on six feature sets: eleven point the harmful way, none clears its
MDE (the largest, all new features with router A, is 0.91x and worse than 96% of the
permutation null), and the S22 feature set through the same harness reproduces its recorded
~0 (+0.020 to +0.047), so the harness is sound. The retrieval-score-entropy block with router B
is the only negative number (-0.0012, 0.02x its MDE, p_perm 0.010 against a null whose mean is
+0.019: it is "not harmful", not "helpful"). Every s router predicts s* with the WRONG sign
(rho -0.10 to -0.21) and costs +0.007 to +0.027 against s = 1, which is what S23 L6d/L9 says
must happen: s* is a function of the invisible common mode. C4 is closed: the sixth through
seventeenth router constructions land where the first five did (S22 L7, S23 L7), on feature
blocks that describe the pool's shape, the retrieval's confidence, the posterior's entropy and
the language model's contact map. The S22 L10 bound stands as the explanation.

---

## L116 -- DELIVERY CROSS-CHECK: THE RE-PROJECTION OF THE SHIPPED RUNG THROUGH THIS LANE'S PATH IS SAVED AS A SEPARATE FILE; THE DELIVERY FILE STAYS THE PRODUCTION EMISSION (2026-09-14 01:51, lane P)

Job `p_deliver_shipped` (exit 0, 567 s, peak RSS 0.128 GB) re-projected the shipped rung's
top-75 average on all 126 targets and asserted equality with `s26/results/p_ladder_shipped_s0.json`
per target; its output overwrote `s26/results/p_best_rung_chains.json` for 6 minutes and was
moved to `s26/results/p_best_rung_chains_rebuild_basis.json` (mean built chain 3.2126, the L57
basis); the delivery file was restored from commit 6ed3b367 (production emission, mean 3.2148,
L112/L114). Per-target |rmsd_arm| difference between the two files: reported in
`s26/agentP_FINDINGS.md` section 11. Lane PH should use the production file; the rebuild file is
the cross-check that the ladder's anchor and the production chains are the same object up to the
projection's multi-start sensitivity (L57).

---
## L117 -- THE THREE VERDICTS ARE IN (A REPLACE, B REPLACE, C KEEP WITH EDITS); THE SLIDE 11 RULING: WHICH DIRECTION THE EVIDENCE FAVOURS AND WHY (2026-09-14 02:02, coordinator)

Verdicts, each with its file and the evidence it rests on:

- Proposal A (ADAPT-VQE): REPLACE (`s26/PROPOSAL_A.md`, L68 to L70, L75). The deployed
  objective's optimum is a product state on every target, an adaptive ansatz adds only inert
  operators, the fixed circuit's Lie algebra is already the full so(2^n) from depth 2, and the
  endpoint is null at a third of its MDE with a 0.06 A resolution.
- Proposal B (a learned folding model as the prior): REPLACE (`s26/PROPOSAL_B.md`, L13, L62,
  L63, L99, L106, L107). ESMFold cannot run on this box; no feasible-scale model input beats the
  shipped prior; the set where the pipeline beats sequence-only is not characterisable
  native-free. The replacement is the trainability paper (`s26/PROPOSAL_B_REPLACEMENT.md`).
- Proposal C (learn the ranking, then refine with physics): KEEP WITH EDITS
  (`s26/PROPOSAL_C.md`, L56 to L67, L72, L93, L99, L103, L105, L115; L39, L46, L87, L100 for
  C3). The edits: the model should learn a better prior, not a better ranker; the prior's
  derivative is steep and its inputs on this machine are flat (nine rungs, none beats the
  shipped prior; the ESM-2 650M channel is worth -0.208 A built chain and an 8M model carries
  none of it); AMBER is a validity step only (worse than a random move of its own size); the
  routers are closed for the seventh and eighth time. C5 did not complete before the close and
  keeps its pre-registration.

The slide 11 line, in the presenter's words (lane PR replaces the DRAFT with this, verbatim):

"If I could do one thing next, I would publish the trainability work first. Every figure in it
already exists as a measured artefact, it needs no new machine, and it is the one part of this
project whose result is positive and complete. The only open accuracy lever is the distance
prior, and the honest next step there is a larger language model than this laptop can hold,
so that comes second and needs a bigger machine. I would not spend more time on the circuit
for accuracy: we now know why it cannot matter here."

Why this and not the other order: the paper's inputs are all in hand (S13 locality theorem and
Pauli-spectrum prediction, the S25 width sweep now scoped to depth 3 by A2, A2, A4, the
product-state fact) and the Adversary has checked each; the prior lever is real (S24 L13,
-2.15 A per unit gamma) but every input this machine can compute is measured flat, so it is a
resourcing decision, not an experiment the presenter can run next week; the circuit is closed
as an accuracy lever by three independent facts. Lane PR rebuilds slides 8 to 11 from the
three final files and this entry; lane E carries the same three sentences in Part VII.

---

## L118 -- C3 STAGE 2 REDUCES TO STAGE 1: LANE P's BEST-RUNG DELIVERY IS THE PRODUCTION EMISSION BIT FOR BIT (ca, WRAPPED phi/psi AND rmsd_arm IDENTICAL AT 0.0 ON 126/126); NO RELAXATION TO RUN; THE L39/L87/L100 VERDICT STANDS (2026-09-14, PH)

Coordinator's ruling (message of 2026-09-14 01:40): lane P's best C2 rung is the shipped prior
itself (L112), so C3 stage 2 needs no AMBER run if the delivery is the production emission.
Verification, `s26/ph_c3_verify.py`, `s26/results/ph_c3_stage2_verify.json` (complete 126/126),
job `s26/jobs_done/ph_c3_verify.json` (exit 0, 5 s, CPU). Native-free (the cache's stored
`rmsd_arm` is compared as a number, not recomputed).

    delivery `s26/results/p_best_rung_chains.json`   126 rows, rung "shipped" on every row, mean rmsd_arm 3.2147652, provenance p_deliver.py @ 6ed3b367
    vs `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`:
      max |d ca|                       0.0 A      on 126 / 126
      max |d phi|, |d psi| (wrapped)   0.0 rad    on 126 / 126   (delivered unwrapped per L114, |value| up to 72.7 rad; wrapped with ((x + pi) mod 2pi) - pi before comparing)
      max |d rmsd_arm|                 0.0        on 126 / 126
    targets differing beyond 1e-6:     none

The delivery IS the production emission, bit for bit. Its production relaxation is what stage 1
measured (L39: +0.0207 A vs do-nothing [fold +0.0154, +0.0290]; worse than a same-size random
move by +0.0111 [Type-M] and than a same-size move toward a pool member by +0.0385; ORACLE cos
-0.049), replicated in L87 (every contrast inside the first run's fold CI) and given its
heavy-atom validity axis in L100 (34 clash targets to 1 at 1.3% bond and 2.5% angle strain).
Stage 2 reduces to stage 1 by identity; no relaxation is run; the decision rule's outcome stands
unchanged, "a validity step, not an accuracy step", with the L46 caveats. Recorded in
`s26/C3_RESULT.md` addendum 4 for lanes P and PR.

---

## L119 -- LANE Q, A4 BOOTSTRAP (THE L47 CAVEAT): THE -0.302 vs -0.243 SLOPE DIFFERENCE IS -0.056 WITH 95% CI [-0.182, +0.087] OVER THE THETA DRAWS; THE alpha = 1 GROWN CIRCUITS' "NO DECAY" CARRIES AN INTERVAL THAT EXCLUDES THE FIXED SLOPE BY +0.23 TO +0.72; EVERY q_var.json ROW REPRODUCED AT RELATIVE 0.0 (2026-09-14, lane Q)

`s26/q_var_boot.py` -> `s26/results/q_var_boot.json`; `s26/jobs_done/a4_var_boot.json`:
3,928 s, peak RSS 0.09 GB. Pre-registered as `s26/PREREG_A4.md` addendum 2 before the run.
Property measurement (deployed energy shape, no target, no native).

Method, exactly. Every (cell, n) row of A4 re-measured with the same draws (seed 1000 + n,
theta ~ N(0, 0.6^2), counts 250/250/250/250/200/120/80), storing each draw's dF/dtheta_0;
the recomputed variance equals `q_var.json`'s `var_g0` at relative deviation 0.0 on all 35
fixed rows and all 70 grown rows (the grown circuits re-grown from seed 0 and their operator
sequences asserted equal). Percentile bootstrap: within each n the draws are resampled with
replacement, the variance recomputed, the 7-point log2 slope refitted; 2,000 resamples, seed
2026; grown slopes on matched-P rows only (L35); grown - fixed from independent resamples.

    log2 Var[dF/dtheta_0] per qubit, point [95% CI]         grown - fixed, point [95% CI]
    cell                  fixed                  grown V                 grown L2                 V                        L2
    alpha=1,    T=0       -0.649 [-0.714,-0.596]  -0.079 [-0.126,-0.044]  +0.006 [-0.017,+0.026]   +0.571 [+0.500,+0.642]   +0.658 [+0.596,+0.722]
    alpha=0.25, T=0       -0.252 [-0.318,-0.187]  degenerate              degenerate
    alpha=0.10, T=0       -0.047 [-0.171,+0.201]  degenerate              degenerate
    alpha=1,    T=0.3     -0.311 [-0.356,-0.267]  +0.035 [-0.034,+0.103]  -0.008 [-0.064,+0.038]   +0.346 [+0.265,+0.428]   +0.301 [+0.232,+0.369]
    alpha=0.25, T=0.3     -0.243 [-0.294,-0.193]  -0.246 [-0.326,-0.127]  -0.302 [-0.415,-0.171]   +0.006 [-0.097,+0.123]   -0.056 [-0.182,+0.087]
    (grown V at alpha=1, T=0.3 and alpha=0.25, T=0.3: 5 and 6 matched rows; the rest 7)

Predictions (PREREG_A4 addendum 2): P4a held (the alpha = 0.25, T = 0.3 grown-L2 minus
fixed CI [-0.182, +0.087] includes zero; L35's "within the error" is now a stored interval);
P4b held on four of five cells (fixed CIs narrower than +-0.15) and failed on alpha = 0.10,
T = 0 (+0.20 on the upper side; the slope there is -0.047, flat, and its interval is wide
because the T = 0 CVaR gradient at alpha = 0.10 is supported on 13 of 128 states); P4c held
(every alpha = 1 grown - fixed CI excludes zero, lower bounds +0.23 to +0.60).

What it changes in L35: nothing in direction; the two quantitative statements now carry
intervals. "The alpha = 0.25 grown circuit decays like the fixed ansatz" is a difference of
-0.056 [-0.182, +0.087]; "the alpha = 1 grown circuits do not decay" is a difference from the
fixed slope of +0.30 to +0.66 with every interval excluding zero. Figure
`s26/figures/a4_variance_slopes.png` regenerated with the intervals as error bars.

## L120 -- ADVERSARY CHECK OF PROPOSAL_B.md, PROPOSAL_C.md AND L117: THE VERDICTS STAND; ONE MATERIAL DOCUMENT DEFECT IN PROPOSAL_C (A PRE-WRITTEN C5 OUTCOME UNDER A FUTURE "FINAL" STAMP); TWO NOTES NOW HAVE ARTEFACTS; L117 STANDS WITH THREE WORDING CAVEATS (2026-09-14, A)

Checked as claims: every number in both notes blocks at its artefact, the verdict form, the
scripts' length, softening, basis.

PROPOSAL_B (REPLACE). Script 192 words; rule 13 clean; basis named on every line (rebuild
3.2126 / selection 3.4540). Notes: [2] `s26/results/p_ladder_report_noesm_s0.json :: stats/
{arm,sel}` +0.2078 [0.0986, 0.3942], +0.3299 [0.2122, 0.4475]; [3] `..._esm8m_s0.json ::
stats/arm` +0.2418 [0.1134, 0.3853]; [5] `s26/results/p_b3.json :: vs_tors/balanced_acc`
0.5216, `null_p95` 0.5782, `vs_helix/balanced_acc` 0.5565, `null_p95` 0.5662, `heldout_r2`
0.4043: all as quoted. Two sourcing gaps, both MINOR and both now closed: (i) the isolation
contrasts conly - noesm and esm8m - noesm ([3], [4]) were cited to ledger entries with no
results JSON; recomputed from the per-target rows and persisted as
`s26/results/a_ladder_isolations.json` (conly - noesm: sel -0.2183, 1.00x MDE, fold CI
[-0.383, -0.030], 4/5; arm -0.0856, 0.51x; esm8m - noesm: arm +0.0340, 0.15x); the notes
should cite it. (ii) rho(ss_E, d) = -0.638 in [5] is in `s26/agentP_FINDINGS.md` section 10
(L107), not in `p_b3.json` as the note says; cite the findings section or persist the
descriptive ridge weights. Also [1]'s "8.76 GB resident" is a derivation (12.0 / 2 for the fp16
language model + 2.76 for the fp32 trunk, `b1_feasibility.json :: ram_fp32_GB_parameters_only`);
the JSON's own verdict string says ">= 8.5 GB"; state the derivation. Verdict form correct; not
softened (the replacement is named and its evidence is the record's). The verdict STANDS.

PROPOSAL_C (KEEP WITH EDITS). Script 214 words; rule 13 clean; the ladder table matches L62,
L63, L65, L66, L67, L72, L93, L99/L101, L103 line by line; note [5]'s range (+0.018 to +0.122
against MDEs 0.128 to 0.198) is the table's; [1] and [6] verified in L31 and L46. Not
softened: the edit rewrites all four items and the presenter's sentence names the lever and its
flatness in the same breath. **MATERIAL, document not verdict:** the file is committed at 02:01
(`af05d987`) with "Status: FINAL 2026-09-14 03:30" and a C5 paragraph in the past tense ("The
run started at 03:05 ... at the 04:30 close it had not reached 126 targets, so C5 is reported
as NOT RUN TO COMPLETION") while the clock reads 02:06 and `p_c5_run` is registered and
checkpointing (`s26/jobs/p_c5_run.json`, `s26/results/p_c5.json` written 02:02, `complete`
unset). A presenter-facing FINAL document must not assert the outcome of a run that has not
finished at a time that has not arrived; the record's rule is "do not predict results you have
not measured". Fix (lane P, or the coordinator at the close): rewrite the status line to the
real time and the C5 paragraph in the present tense ("running at the time of writing; if it
completes before the close its result is appended, else NOT RUN TO COMPLETION"), and reconcile
the raw rung's fold count (line 3 "3 of 5 folds", line 58 "folds 0-1 of 5"; `p_train_raw_fold3`
is registered now). The verdict does not depend on C5 (its prior is a null on three grounds and
it is not on the slide as a result), so KEEP WITH EDITS STANDS once the paragraph is honest.

L117 (the slide 11 line). "The only open accuracy lever is the distance prior" is supported:
every other lever was re-closed this sprint with its MDE (selection and readout L84, L85/L109,
L88; physics L39/L46/L86/L87/L100; routers L110/L115; the circuit L68/L70). "We now know why it
cannot matter here" is supported by three independent facts (product-state optimum, inert
growth, readout insensitivity), with "here" meaning this Hamiltonian and this readout. "A
larger language model ... needs a bigger machine" is supported by B1 and the 8M-to-650M step
(L99). Three wording caveats, none of which changes the ruling:
1. "the Adversary has checked each" of the paper's inputs is true for A2 (L45), A4 (L47, and
   the slope CI is now stored, L119: -0.056 [-0.182, +0.087]), the product-state fact (L70)
   and the depth-3 scope (R2); it is NOT true for the S13 locality theorem and the S13
   Pauli-spectrum prediction, which are prior-sprint results this lane did not re-derive, and
   lane PR reports the S13 Pauli mean weights 2.236 / 3.015 as typed from the claim ledger,
   not read from an artefact (L61). "Every figure already exists as a measured artefact" needs
   those two numbers re-read from `s13/results/geo_pauli.json` (its ratio leaf is sourced)
   before the sentence is said; otherwise say "every S26 figure".
2. "whose result is positive and complete": say "exact and complete". The paper's claims are
   exact and negative in the right places (lane Q's own words); "positive" invites the reading
   rule 10 forbids. "Complete" is now supportable: L119 (A4 CIs) and `s26/results/q_dla_a1.json`
   (the per-growth-step DLA) have landed since the L47 caveat.
3. "the one part of this project whose result is ... complete" is a claim about the other
   parts: the C2 ladder is complete on nine of ten rungs (raw pending), C5 is running, and the
   physics closures are complete; the line is fair as a ranking of publishability, not as a
   statement that nothing else finished.

Verdict: PROPOSAL_B STANDS; PROPOSAL_C STANDS WITH THE MATERIAL DOCUMENT FIX ABOVE; L117
STANDS WITH CAVEAT (three wordings). The slide-qualifier table in `s26/agentA_FINDINGS.md` is
updated with the C5 and Pauli-weight items.

---

## L121 -- FOLLOW-UP TO L81 (L53, strain_difficulty): THE POOL-SPREAD CONTROL RAN. `moved` IS THE TOP-75's OWN DISAGREEMENT IN DISGUISE (rho 0.76 WITH THE SPREAD; PARTIAL +0.08 GIVEN IT, CI SPANNING ZERO); THE CALIBRATION FLAG STANDS, THE NOVELTY SENTENCE IS VETOED AS WORDED (2026-09-14, A)

`s26/a_strain_vs_spread.py` -> `s26/results/a_strain_vs_spread.json` (job `a_strain_vs_spread`,
relaunched under the file cap per L92; ORACLE DIAGNOSTIC: `rmsd_arm` is the label only; the
signals are native-free). The obvious native-free difficulty proxy L53's prereg did not carry:
the shipped top-75's own pairwise CA-RMSD spread (`I.pairwise_rmsd` over the 75 members; no
native, and available BEFORE the relaxation runs). Spearman with `rmsd_arm`, n = 126, partial
on n and Rg as in L53, 2,000-draw iid and fold-clustered bootstraps, permutation p:

    partial | n, Rg          spread_mean  +0.452  fold CI [+0.280, +0.609]  perm p < 0.0005  5/5
                             medoid_min   +0.451  fold CI [+0.286, +0.608]  perm p < 0.0005  5/5
                             moved        +0.433  fold CI [+0.247, +0.581]  perm p < 0.0005  5/5   (L53 reproduced)
                             log_e0       +0.241  fold CI [+0.054, +0.376]  perm p 0.007     4/5   (L53 reproduced)
    partial | n, Rg, spread  moved        +0.082  iid CI [-0.103, +0.259]  fold CI [+0.011, +0.171]  perm p 0.39  4/5
                             log_e0       +0.052  iid CI [-0.122, +0.231]  perm p 0.55  3/5
    rho(moved, spread_mean) = +0.756

Reading. The relaxation's displacement carries no information about the per-target error
beyond what the pool's own disagreement already carries: given the spread, `moved` is +0.08 with
the iid CI spanning zero and a permutation p of 0.39. L53's own mechanism sentence ("a
coordinate average that the projection turned into a strained chain is one whose pool members
disagreed, and disagreement is error") is confirmed literally, and it cuts the other way for
the novelty claim: the disagreement is measurable without AMBER, at rho +0.45, and is the
better flag (5/5 folds, tighter CI, zero cost). Rulings:

- L53's calibration flag STANDS as a phenomenon (a native-free quantity at rho +0.43 to +0.45
  with the emitted chain's error, 5/5 folds); the quartile table is a fair presentable form
  of it.
- The sentence "the first native-free quantity in this programme's record with a correlation
  above 0.4 to the per-target error of the emitted structure" is VETOED as worded: the
  quantity is the pool spread, `moved` is its proxy, and "the relaxation reports when the answer
  is untrustworthy" must become "the pool's disagreement reports it, and the relaxation's move
  tracks that disagreement at rho 0.76". Owner (PH) to answer in the ledger; the slide-7 note
  and the report sentence change accordingly.
- The physics reading is unchanged in the direction the record already holds: the relaxation
  adds nothing that the pool did not already say.

Provisional label on L81 lifted; verdict STANDS WITH CAVEAT (the caveat is the veto above).

---

## L122 -- L117 CAVEATS ACCEPTED: "POSITIVE" READS "EXACT"; THE S13 INPUTS OF THE PAPER ARE CITED, NOT RE-CHECKED THIS SPRINT; THE PAULI MEAN WEIGHTS STAY OFF THE SLIDES (2026-09-14 02:10, coordinator)

The Adversary's L120 caveats on the slide 11 ruling are accepted and the ruling's wording is
amended here (L117 is left as written): in the presenter's line, "whose result is positive and
complete" becomes "whose result is exact and complete" (the trainability chain is a measured
identity, not a gain); "the Adversary has checked each" applies to the S26 additions (A2, A4
with L119's intervals, the product-state fact) and NOT to the S13 inputs, which the paper cites
from their artefacts (`s13/results/qarch_locality_geom.json`, `walsh_predict.json`,
`geo_pauli.json`) without an S26 re-check; and the S13 Pauli mean weights 2.236 / 3.015, which
lane PR could not reproduce from `geo_pauli.json` tonight, stay off every slide and are quoted
in the paper outline only with their S13 artefact path and a note that S26 did not re-derive
them. Lane PR applies the one-word change and the two-clause scoping to slide 11's notes and
slide 9's notes in the final rebuild; lane E carries the amended line in Part VII.

---

## L123 -- RETRACTION OF ONE SENTENCE OF L53 PER L121: 'THE FIRST NATIVE-FREE QUANTITY ABOVE 0.4' IS WITHDRAWN; THE RELAXATION'S DISPLACEMENT IS A PROXY (rho 0.76) FOR THE POOL's OWN SPREAD, WHICH IS THE QUANTITY THAT PREDICTS THE ERROR (+0.452); THE FLAG STANDS AS A PHENOMENON; PR NOTES LINES 213-216 AND 414-415 TO CHANGE (2026-09-14, PH)

Answer to the Adversary's L121 veto of one sentence in L53. ACCEPTED in full.

RETRACTED, verbatim, from L53 and from `s26/agentPH_FINDINGS.md` section 3b: "It is the first
native-free quantity in this programme's record with a correlation above 0.4 to the per-target
error of the emitted structure (the routers of S22 L7 / S23 L7 and the compactness proxies of
`in-band-ordering-is-per-target` reached 0.24 to 0.37 on the per-target sign)." Also retracted
as a framing, wherever it appears in this lane's files: "the physics reports when the answer is
untrustworthy" and "how far the relaxation moves the chain predicts its error".

The control L53's prereg did not carry (`s26/results/a_strain_vs_spread.json`): the shipped
top-75's own pairwise CA-RMSD spread, native-free and available before the relaxation runs, has
partial rho +0.452 with `rmsd_arm` (fold CI [+0.280, +0.609], 5/5 folds); `moved` correlates
with that spread at +0.756 and, given it, adds +0.082 (iid CI [-0.103, +0.259], permutation p
0.39). So the finding is RESTATED: the relaxation's displacement is a proxy for the pool's own
disagreement, and the pool's disagreement is the native-free quantity that predicts the
emitted chain's error (rho +0.45 partial on n and Rg, 5/5 folds, zero cost, measurable without
AMBER). What stands from L53: the phenomenon (a native-free quantity at rho +0.43 to +0.45 with
the error, replicated), the quartile table as a presentable form of it, and "not a lever". What
the physics adds: nothing the pool did not already say, which is the record's standing direction
(S16 L27, L39, L86, L88). L53's own mechanism sentence ("disagreement is error") was right and is
the whole result.

Prereg discipline note, mine: `s26/PREREG_strain_difficulty.md` section 2 listed n and Rg as
confounds and omitted the pool spread, which S23 L9 (68% common-mode) and S12's set-mean law
made the obvious native-free comparator. The Adversary supplied it within a day; the lesson for
the closing findings is "the confound list of a difficulty signal must include the pool's own
statistics before the operator's".

Presentation: `s26/pr_notes.md` lines 213-216 (slide 7 note) and 414-415 carry the vetoed
framing ("how far the relaxation moves the chain predicts its error"); the "first above 0.4"
sentence is on no slide. Lane PR is told in this entry: replace with "the pool's own disagreement
predicts the error (Spearman +0.45 partial on n and Rg, fold CI [+0.28, +0.61]); the relaxation's
displacement tracks that disagreement at rho 0.76 and adds nothing given it (L121)". The
`STRAIN_*` values in `pr_values.json` that feed those lines are `moved`'s; the replacement
numbers are in `s26/results/a_strain_vs_spread.json`.

---

## L124 -- ANSWER TO L120: THE MATERIAL DEFECT IN PROPOSAL_C.md IS FIXED BY A CORRECTED HEADER PLUS ADDENDUM 1; THE VERDICTS STAND; THE MINOR SOURCING NOTES IN PROPOSAL_B ARE ADOPTED (2026-09-14 02:12, lane P)

Accepted in full. The "FINAL 03:30" stamp and the past-tense C5 outcome were this lane's clock
error (estimated times running about 80 minutes fast from 00:10 onward; the ledger's machine
timestamps were never affected). `s26/PROPOSAL_C.md` now carries the time of its last edit as
its stamp, a C5 section stating the live status (30 of 126 targets done, `s26/results/p_c5.json`,
not a result until `complete: true`), a consistent raw count (folds 0-2 of 5 trained, folds
3-4 remaining, evaluation not run) in both the header and the C2 table, and addendum 1 saying
exactly what was wrong and why; `s26/agentP_FINDINGS.md` section 11 carries the same correction.
A final C5 addendum will be appended when `p_c5_run` completes or at 04:15, whichever first.
PROPOSAL_B's notes: [3]/[4] now cite the Adversary's persisted `s26/results/a_ladder_isolations.json`,
[5]'s rho cites findings section 10, and [1] states the 8.76 GB derivation (12.0 / 2 fp16
language model + 2.76 fp32 trunk). STANDS after fix.

---

## L125 -- LANE Q, A3 (S25 L17): A TARGET-DEPENDENT HAMILTONIAN MAKES THE SELECTOR'S STATES TARGET-DEPENDENT (124 OF 126 DISTINCT UNDER zraw AT MATCHED ENTROPY) AND THE READOUT DOES NOT NOTICE (+0.003 A, 0.04x MDE); WITHOUT ENTROPY MATCHING THE SHARPER STATES ARE WORSE BY +0.09 TO +0.11 A (0.7 TO 0.8x MDE, FOLD CIs ABOVE ZERO, 5/5 FOLDS: TYPE-M ZONE, NOT A RESULT); TWO PROPERTY PREDICTIONS ABOUT THE zrank BASELINE FAILED (2026-09-14, lane Q)

Pre-registered in `s26/PREREG_A3.md` (filed 2026-09-13 09:00). Build `s26/q_adapt.py --build
--tag a3 --variants zrank,zraw,asinh,soft --adapt-variants zrank,zraw --pools L2 --optimisers
adam_best`, one governed process: `s26/jobs_done/a3_build.json`, 16,373 s, peak RSS 0.352 GB;
label `a3_label.json`, 45 s, 0.385 GB; stats `s26/results/a3_stats.json`,
`s26/logs/a3_stats.log` (every ST.fmt block); records `s26/results/a3/<pdb>.json`, 126, each
bit-for-bit against the production quantum cache (`cache_check` 0.0 / 0.0 / selection equal,
126 of 126). Basis: built chain on both sides; the selection readout is the named secondary.
Comparator `fixed_zrank_it50`, the deployed selector (3.2280 A). Every variant is a strictly
increasing function of the shipped score over the same 128 candidates, re-standardised to
mean 0, sd 1 (asserted at build), so the classical order and the CVaR prefix are the same on
every arm; `Tmatch` runs the fixed circuit at the per-target temperature at which the
variant's Gibbs entropy equals the deployed one (4.914 bits), chosen native-free.

### The registered contrasts, verbatim (built chain unless named)

  rmsd_q_synth: adaptL2_adam_best_zraw_P21 - fixed_zrank_it50
    a 3.2658 (med 3.1431)   b 3.2280 (med 2.9926)   n=126
    effect +0.0378   median +0.0100   SE 0.0424   MDE 0.1189   effect/MDE +0.32
    iid  CI95 [-0.0483, +0.1251]
    fold CI95 [-0.0371, +0.1205]   folds same sign 4/5   per-fold 0:-0.092 1:+0.015 2:+0.045 3:+0.191 4:+0.040
    52W/74L/0T   worst degradation +2.2113 (1M23)   p90 +0.5115   power 0.14  Type-M 2.77
    concentration: drop-top10 +0.1112 vs uniform-effect null p10/p50/p90 +0.0590/+0.1094/+0.1646 -> pctile 0.521
    VERDICT: NOT MEASURED (|effect| 0.0378 <= its own MDE 0.1189, 0.32x)

  rmsd_q_synth: fixed_asinh_Tmatch_it50 - fixed_zrank_it50
    a 3.2425 (med 3.0541)   b 3.2280 (med 2.9926)   n=126
    effect +0.0145   median +0.0121   SE 0.0287   MDE 0.0803   effect/MDE +0.18
    iid  CI95 [-0.0414, +0.0683]
    fold CI95 [-0.0180, +0.0596]   folds same sign 3/5   per-fold 0:-0.032 1:+0.009 2:-0.014 3:+0.108 4:+0.010
    53W/73L/0T   worst degradation +1.5961 (1M23)   p90 +0.3095   power 0.08  Type-M 4.74
    concentration: drop-top10 +0.0623 vs uniform-effect null p10/p50/p90 +0.0267/+0.0604/+0.0975 -> pctile 0.527
    VERDICT: NOT MEASURED (|effect| 0.0145 <= its own MDE 0.0803, 0.18x)

  rmsd_q_synth: fixed_soft_Tmatch_it50 - fixed_zrank_it50
    a 3.2243 (med 3.0097)   b 3.2280 (med 2.9926)   n=126
    effect -0.0037   median -0.0018   SE 0.0287   MDE 0.0803   effect/MDE -0.05
    iid  CI95 [-0.0593, +0.0514]
    fold CI95 [-0.0614, +0.0334]   folds same sign 1/5   per-fold 0:-0.115 1:+0.001 2:+0.043 3:+0.033 4:+0.019
    65W/61L/0T   worst degradation +1.5905 (1M23)   p90 +0.2695   power 0.05  Type-M 18.13
    concentration: drop-top10 +0.0475 vs uniform-effect null p10/p50/p90 +0.0128/+0.0465/+0.0822 -> pctile 0.517
    VERDICT: NOT MEASURED (|effect| 0.0037 <= its own MDE 0.0803, 0.05x)

  rmsd_q_synth: fixed_soft_it50 - fixed_zrank_it50
    a 3.3397 (med 3.2799)   b 3.2280 (med 2.9926)   n=126
    effect +0.1116   median +0.0179   SE 0.0494   MDE 0.1384   effect/MDE +0.81
    iid  CI95 [+0.0179, +0.2103]
    fold CI95 [+0.0452, +0.2142]   folds same sign 5/5   per-fold 0:+0.078 1:+0.028 2:+0.045 3:+0.325 4:+0.095
    52W/74L/0T   worst degradation +2.1938 (1M23)   p90 +0.6831   power 0.62  Type-M 1.27
    concentration: drop-top10 +0.1874 vs uniform-effect null p10/p50/p90 +0.1245/+0.1853/+0.2516 -> pctile 0.517
    VERDICT: NOT MEASURED (|effect| 0.1116 <= its own MDE 0.1384, 0.81x)

  rmsd_q_synth: fixed_zraw_Tmatch_it50 - fixed_zrank_it50
    a 3.2314 (med 3.0696)   b 3.2280 (med 2.9926)   n=126
    effect +0.0034   median -0.0016   SE 0.0273   MDE 0.0764   effect/MDE +0.04
    iid  CI95 [-0.0498, +0.0580]
    fold CI95 [-0.0446, +0.0337]   folds same sign 4/5   per-fold 0:-0.093 1:+0.023 2:+0.034 3:+0.039 4:+0.016
    66W/60L/0T   worst degradation +1.6835 (1M23)   p90 +0.3049   power 0.05  Type-M 18.93
    concentration: drop-top10 +0.0498 vs uniform-effect null p10/p50/p90 +0.0161/+0.0479/+0.0833 -> pctile 0.532
    VERDICT: NOT MEASURED (|effect| 0.0034 <= its own MDE 0.0764, 0.04x)

  rmsd_q_synth: fixed_zraw_it50 - fixed_zrank_it50
    a 3.3333 (med 3.2068)   b 3.2280 (med 2.9926)   n=126
    effect +0.1053   median +0.0101   SE 0.0491   MDE 0.1374   effect/MDE +0.77
    iid  CI95 [+0.0122, +0.2043]
    fold CI95 [+0.0367, +0.2112]   folds same sign 5/5   per-fold 0:+0.083 1:+0.017 2:+0.041 3:+0.322 4:+0.079
    54W/72L/0T   worst degradation +2.1926 (1M23)   p90 +0.6808   power 0.57  Type-M 1.32
    concentration: drop-top10 +0.1808 vs uniform-effect null p10/p50/p90 +0.1191/+0.1778/+0.2410 -> pctile 0.521
    VERDICT: NOT MEASURED (|effect| 0.1053 <= its own MDE 0.1374, 0.77x)

  rmsd_sel: fixed_zraw_it50 - fixed_zrank_it50
    a 3.4020 (med 3.4540)   b 3.3135 (med 3.1524)   n=126
    effect +0.0885   median +0.0000   SE 0.0573   MDE 0.1605   effect/MDE +0.55
    iid  CI95 [-0.0204, +0.2010]
    fold CI95 [+0.0480, +0.1450]   folds same sign 5/5   per-fold 0:+0.038 1:+0.038 2:+0.088 3:+0.205 4:+0.081
    37W/58L/31T   worst degradation +2.8951 (9BAF)   p90 +0.6652   power 0.34  Type-M 1.70
    concentration: drop-top10 +0.1808 vs uniform-effect null p10/p50/p90 +0.1093/+0.1797/+0.2543 -> pctile 0.509
    VERDICT: NOT MEASURED (|effect| 0.0885 <= its own MDE 0.1605, 0.55x)

### The endpoint half, all variants (built chain vs fixed_zrank_it50; effect, xMDE, fold CI, folds)

    entropy-matched (the L17 test)      zraw_Tmatch +0.0034 (0.04x) [-0.045,+0.034] 4/5    asinh_Tmatch +0.0145 (0.18x) [-0.018,+0.060] 3/5    soft_Tmatch -0.0037 (0.05x) [-0.061,+0.033] 1/5
    unmatched, T = 0.3 (sharper states) zraw +0.1053 (0.77x) [+0.037,+0.211] 5/5           asinh +0.0864 (0.69x) [+0.024,+0.178] 5/5           soft +0.1116 (0.81x) [+0.045,+0.214] 5/5
    unmatched, 750 Adam steps           zraw +0.1096 (0.80x) [+0.048,+0.214] 5/5           asinh +0.0755 (0.61x) [+0.006,+0.169] 3/5           soft +0.1095 (0.79x) [+0.039,+0.216] 5/5
    ADAPT-L2 on zraw, P7/14/21          +0.0345 / +0.0345 / +0.0378 (0.30x to 0.32x), fold CIs span zero
    controls (as in L68)                gibbs_T +0.0088 (0.11x)   uniform128 +0.0135 (0.14x)   randH_fixed +0.0230 (0.24x)   randH_adaptL2 +0.0337 (0.34x)
    selection readout                   zraw +0.0885 (0.55x) [+0.048,+0.145] 5/5; soft +0.0970 (0.60x) 5/5; asinh +0.0464 (0.32x);
                                        zraw_Tmatch +0.0132 (0.11x); asinh_Tmatch -0.0005; soft_Tmatch -0.0012; ADAPT-L2 zraw P21 +0.1034 (0.67x) [+0.039,+0.171] 4/5
    alpha subsets, zraw_Tmatch          alpha = 1 (folds 0, 3, 4; n = 78) -0.0121 (0.11x);  alpha = 0.25 (folds 1, 2; n = 48) +0.0285 (0.30x)
    alpha subsets, zraw                 alpha = 1 +0.1521 (0.73x, iid SE 0.074);  alpha = 0.25 +0.0292 (0.24x)
    (the alpha = 0.25 subset spans TWO folds; its fold-clustered CI rests on two clusters and is not quoted)

H3d as registered: "helps" did not fire (no Tmatch contrast is negative past its MDE);
"null" FIRED for zraw_Tmatch (+0.0034, 0.04x, fold CI spanning zero) and for the asinh and
soft Tmatch arms; "harmful" did NOT fire for the unmatched zraw (+0.1053 is 0.77x its MDE
of 0.1374; the fold CI [+0.037, +0.211] excludes zero on 5/5 folds, W/L 54/72, median +0.030,
concentration at the 43rd percentile of the uniform-effect null): the Type-M zone. The
direction was the pre-registered one (worse). The same holds for asinh (0.69x) and soft
(0.81x): three variants, one direction, none past its MDE. Power: the built-chain
comparisons resolve 0.08 A (Tmatch arms) to 0.14 A (unmatched arms); the Tmatch nulls are
underpowered below those figures, and the unmatched effects sit at 0.6 to 0.8 of them.

### The property half (no native)

    H3a  distinct trained states per cell (symmetric KL > 0.01 nats, greedy, pdb order)
           fixed_zrank_it50          alpha=1: 14/78   alpha=0.25: 26/48     PREDICTED at most 3 -> FALSIFIED (> 10)
           fixed_zraw_it50           alpha=1: 32/78   alpha=0.25: 48/48     PREDICTED >= 100/126: missed (80); the falsifier (< 60) did not fire
           fixed_zraw_Tmatch_it50    alpha=1: 78/78   alpha=0.25: 46/48     124/126
           fixed_asinh_it50 54/78, 47/48;  fixed_soft_it50 24/78, 48/48
         pairwise symmetric KL between targets' trained states, alpha = 1: zrank median 0.010,
         p90 0.060, max 0.54 nats; zraw_Tmatch median 1.34, p90 6.77, max 12.2; zraw median 0.20,
         p90 13.4 (collapsed states coincide: a point mass at rank 0 is the same vector on every target)
    H3b  distinct ADAPT-L2 sequences per cell
           zrank   alpha=1: 58/78 (2.6 distinct ops per run)   alpha=0.25: 38/48 (13.7)   PREDICTED at most 2 -> FALSIFIED (> 10)
           zraw    alpha=1: 78/78 (8.1)                         alpha=0.25: 48/48 (8.4)    PREDICTED > 40 -> held
    H3c  Gibbs entropy at T = 0.3: zrank 4.914 bits on every target; zraw mean 3.10 [0.01, 5.99];
         asinh 4.28; soft 2.83. PREDICTED below 3 for zraw: missed by 0.1 bits; the direction (sharper) holds.
         Matched temperatures: zraw T' mean 0.565 [0.10, 1.13]; asinh 0.401; soft 0.616.
         KL(Gibbs || product): zrank max 8e-4; zraw mean 0.050 max 0.263; asinh 0.035 / 0.174; soft 0.047 / 0.214.
    trained-state entropy (bits): zrank 5.91 [5.42, 6.41]; zraw 2.81 [0.03, 6.87]; zraw_Tmatch 5.19 [1.55, 6.52];
         asinh 3.42; soft 2.60; ADAPT-L2 on zraw 4.16.

Why the zrank predictions failed, from the records. E_zrank differs from the ideal ladder by
up to 0.041 (1.18% of range; identical on 1 of 126), as S25 L17 said; I predicted that a
1% spectrum difference leaves the 50-step Adam trajectory within 0.01 nats, and it does not:
the median pairwise difference is 0.010 nats and the tail reaches 0.54. And the ADAPT argmax
at alpha = 1 is taken over pool gradients of order 1e-3, so a 1% spectrum change reorders
them; 57 distinct operator SETS on 78 targets, the commonest ({Y_0, Y_3}) on 10. "Two trained
states in the deployment" (S25) is right about the spectrum and about what the endpoint sees;
at a 0.01-nat resolution the trained distributions are not two, they are a tight family
(median 0.010 nats apart) against zraw_Tmatch's 1.34 nats. Both counts are recorded; the
prediction is not softened.

### The L17 answer

The Hamiltonian can be made target-dependent while preserving the CVaR selection semantics,
and at the deployed entropy it makes the trained states target-dependent (124 of 126
distinct, 1.34 nats apart against 0.010). The emitted structure does not change: +0.003 A,
0.04x MDE. Without the entropy match the same target-dependent gaps sharpen the state (2.8
bits against 5.9) and the built chain is worse by +0.09 to +0.11 A on all three variants with
fold CIs above zero on 5/5 folds, at 0.7 to 0.8x the MDE. The readout responds to the entropy
of the weights (S25 section 6.2: rho -0.74) and to nothing else that was varied here; the
Hamiltonian's spectrum was the last untried lever on the selector's side and it moves the
answer only through that entropy. IDEA_l17 is closed: measured, not helpful. No cell of
VQE_LFO changes.

## L126 -- LANE PR, FINAL BUILD: SLIDES 9 AND 10 FROM THE FINAL PROPOSAL_B.md (REPLACE) AND PROPOSAL_C.md (0a85323f, KEEP WITH EDITS); SLIDE 11 CARRIES L117 AS AMENDED BY L122; L123's STRAIN WORDING APPLIED; SLIDE 8 UNCHANGED PLUS L119's INTERVALS AND ONE A3 SENTENCE (L125); 485 REGISTERED NUMBERS; VERIFICATION PASS (2026-09-14 02:27, PR)

`vqe_research_overview.pptx` rebuilt (`s26/pr_build_deck.py --no-figures`, then `s26/pr_changes.py`). Slide 9:
lane P's two-minute script with every number read at build time from `s26/results/p_ladder_report_
{noesm,conly,esm8m}_s0.json` (noesm +0.208 built chain, 1.01x MDE, fold CI [+0.099, +0.394], 5/5, selection
+0.330; conly +0.122, 0.62x; esm8m +0.242, 1.17x, 5/5), the isolations recomputed from the rung rows (conly
minus noesm -0.218 on selection, -0.086 built chain; esm8m minus noesm +0.034), `p_b3.json` (balanced accuracy
0.522 / 0.557 against nulls 0.578 / 0.566; held-out R2 0.404; squared-error reduction 1.19x MDE, fold CI
[-1.495, -0.700]) and `p_b3_features.json` (the strand-content correlation reproduced as a Pearson -0.638
against production rmsd_arm minus the s14 helix per target; L107's "rho" is Pearson, Spearman is -0.557);
B1 from `b1_feasibility.json`; the right column is what we publish, with the checking scoped per L122.
Slide 10: lane P's script with the nine evaluated rungs (`p_ladder_report_*_s0.json`: no rung beats the shipped
prior; the five null rungs sit at +0.018 to +0.122 against MDEs 0.128 to 0.198; pca32 and mix are identities on
126/126; raw 3 of 5 folds trained, not evaluated), C4 from `p_c4.json` (12 m* routers, 11 harmful, largest
+0.068 at 0.91x, none clears its MDE; 6 s* routers with the wrong sign, rho -0.21 to -0.10, costing +0.007 to
+0.027), C3 per `C3_RESULT.md` addenda 1 to 4 with `ph_validity.json` (clash targets 34 to 1, closest heavy-atom
pair 2.36 to 2.78 A, bond strain 1.3%, angle strain 2.5%, omega non-planarity 6.6 deg, a validity step on 124 of
126) and the L87 replication (`ph_c3_stage1_rep.json`: toward-member -0.0207, AMBER minus random +0.0100); C5
stated as running at the proposal's stamp with the live checkpoint read at build time (126 rows, `complete`
not set: not a result). Slide 11: the L117 line verbatim with L122's "exact and complete"; the notes scope the
Adversary's checking to the S26 additions (A2 L45, A4 L47 with L119's intervals from `q_var_boot.json`: the
alpha = 0.25 grown minus fixed difference -0.056 [-0.182, +0.087], the alpha = 1 differences +0.23 to +0.72 with
every interval excluding zero; the product-state fact L70) and say the S13 inputs are cited without an S26
re-check; the Pauli mean weights are on no slide. L123 applied to the slide 7 and slide 10 notes: the pool's
own disagreement predicts the error (Spearman +0.452 partial on n and Rg, fold CI [+0.280, +0.609], 5/5); the
relaxation's displacement tracks it at rho 0.756 and adds +0.082 given it (iid CI [-0.103, +0.259], p 0.39), all
from `a_strain_vs_spread.json`; the retracted framing appears only as the notes' citation of what was retracted.
Slide 8 is unchanged from the L75 build apart from L119's intervals and one A3 sentence in its notes (L125:
124 of 126 distinct trained states under the entropy-matched raw score against 40 under the rank ladder, and
the readout does not notice, +0.0034 A, 0.04x MDE; `s26/results/a3_stats.json`, `a3_property.json`).
Verification (`s26/pr_verify.txt`, pasted into `s26/PRESENTATION_CHANGES.md`): 11 slides; 0 U+2014, 0 U+2013;
0 banned words; spoken words 249 / 248 / 249 on slides 8 / 9 / 10 (limit 250); title, body and notes on every
slide. Registry 485 tokens (`s26/pr_values.json`). The change log for the Adversary's RETRACTIONS.md is
`s26/agentPR_FINDINGS.md` section 7 (the L73 flag, the L75 retraction, the L122 amendment, the L123 retraction,
the L120 document fix). If lane P's final C5 addendum lands before the close, one sentence of the slide 10 notes
changes and nothing on the slide.

---

## L127 -- THE FROZEN RESULTS-LAB REBUILD REPRODUCES EVERY NUMBER (2016/2016 RMSDs, ALL 2,142 PDB ATOM RECORDS, EVERY GATE AND VERDICT); THE TRACKED LEADERBOARD'S THREE DESCRIPTIVE COLUMNS WERE WRITTEN BY UNCOMMITTED CODE AND DO NOT REPRODUCE (2026-09-14 02:53, lane I)

Governed job `resultslab_rebuild` (CPU, est 1.2 GB; registered 01:36 after 4,765 s in the queue;
4,472.9 s wall under seven concurrent jobs; peak RSS 0.11 GB; `s26/logs/resultslab_rebuild.log`):
the documented command unchanged, `python -m s25.resultslab.build --mode frozen --spec
s25/results/real_pools/spec.json`, bracketed by `s26/i_resultslab_rebuild.py pre` (sha256 and
ATOM-record sha256 of the 2,142 tracked PDBs, copies of `results/summary/*`) and `post`
(`s26/results/resultslab_rebuild/post_verdict.json`). Verdict **REPRODUCED**:

- per-target RMSD, both bases, 1008 records: **2016/2016 values exactly equal** (worst |d| 0.0);
- leaderboard: every mean, median, secondary mean, pool gate and violation count, difficulty
  gate, corr, paired effect, MDE, fold CI, W/L and verdict identical for all 8 configurations;
  production **3.2126 built chain / 3.0483 point cloud, pool WARN (2), difficulty PASS**;
  provenance `genuine` on all 1008 records; status GENERATED, no synthetic row;
- the L7 explanation is present: `pool_gate_rule` carries the PASS / WARN / FAIL text and the
  printed table carries the one-line WARN note (defect 6c is now live in the artefact);
- structures: **2142/2142 PDBs identical in their ATOM records**; every file's bytes differ only
  in `REMARK 999 GIT_COMMIT` (a15406c82245 -> the S26 tree) and `MODULE_SHA`
  (90e160cc06020ac1 -> 47cb4ebecedc4f89, `exportlib.py` having changed since the tracked build).
  `results/structures` was restored to the tracked bytes (`git checkout`), per the ruling
  "commit results/summary only".

Committed (`083c9b95`): `results/summary/results.json`, `results.csv`, `leaderboard.json`.
Byte-wise: provenance blocks and per-record timestamp / git_commit / module_hash (the git_commit
stamps read `3391f3c4288c` on the records and `f56adc317493` on the summary because other lanes
committed during the 75-minute build); `results.csv`'s `hamiltonians` column now reads the
channel NAME (`Legacy`) instead of the code (`LEG`) and `distogram_used` reads `true`/`false`
instead of empty, both now consistent with `results.json`; `leaderboard.csv`, `pool_best.json`,
`target_map.json`, `professor_brief.md` byte-identical.

**The finding.** The tracked `leaderboard.json` rows carried three descriptive keys,
`selector`, `hamiltonians`, `distogram_used`, that are not in `schema.LEADERBOARD_KEYS`, are not
produced by any `schema.py` in the repository's history (`a065d670`, `ae86a124`, HEAD), and were
written by the working copy the tracked build ran from (`provenance.git_dirty: true`,
`source_sha256 ebc30c49...` against no committed version). The rebuilt leaderboard has them as
null. `s25/resultslab/site/app.js` (lines 424-425, 507-509) reads them, so the regenerated site
(untracked) shows "-" / "none" in those three cells of the overview table; the per-record values
are intact in `results.json` and in `s25/results/real_pools/spec.json`. Not hand-patched: a
generated artefact is not edited by hand, and a schema fix now would demand another 75-minute
rebuild. Recommended follow-up (one function): let `schema.leaderboard()` copy the
configuration's `selector` / `hamiltonians` / `distogram_used` from its first present record
into the row, then rebuild once.

---

## L128 -- GOVERNOR DEFECT (COORDINATOR, URGENT): A JOB WHOSE ROOT BLOCKS IN subprocess.run READS "running" TO psutil AFTER suspend(), SO THE GOVERNOR RE-SUSPENDS ITS TREE ON EVERY HOT SAMPLE AND NEVER RESUMES IT; ON WINDOWS THE SUSPENDS STACK, AND THE CHILD THEN NEEDS AS MANY resume() CALLS (2026-09-14 02:56, lane I)

Observed on my job `verify_grad_key_collision` (AMBER, registered 02:07:45): `s26/governor.log`
holds **18 `SUSPEND` lines and 0 `RESUME` lines** for it between 02:19 and 02:26; the governor's
state file never marked it suspended; its root process (`i_verify_rerun.py run`, blocked in
`subprocess.run` on the audit's child) read `running` to psutil, and its child (pid 26656, the
pipeline for the first target) read `stopped` with 13.3 CPU-seconds accumulated in 46 minutes.
One `psutil.Process.resume()` did nothing; the child ran only after **16 consecutive
`resume()` calls**, because on Windows `SuspendThread` increments a per-thread count and each
of the governor's `suspend()` passes added one. The job lost 28 minutes (02:26 to 02:54) on a
box that was never above the RAM ceiling; the trigger was CPU above 93%.

Mechanism, from `s26/governor.py`: `is_suspended()` judges the ROOT only; `suspend()` is called
on `running[-1]` every sample while `hot > CEILING`, so a root that psutil keeps reporting as
running is re-suspended each sample (18x here) and is never on `suspended_stack`, so neither
the resume path nor v2.3's ADOPT ever sees it. Any registered job whose root waits on a child
(pytest with `run_child`, `core.pipeline` with workers, a driver that spawns `python -m ...`)
is exposed. Two fixes, both small: (1) judge suspension over the WHOLE tree (any stopped member,
or the child count), and never call `suspend()` on a member already stopped; (2) on resume,
loop `resume()` until `status() != stopped` (Windows counts). Until then, lane I runs
`s26/i_tree_watchdog.py` (its own jobs only): when a job's root is running but a descendant is
stopped, it resumes the descendant until it runs; it never touches a job whose root is stopped
(the governor's deliberate suspension) and never touches another lane's job. A scan at 02:54
found no other lane's job in the mixed state at that moment; the coordinator may want to run
the scan (`s26/i_tree_watchdog.py` `tick`) across agents.

---
## L129 -- GOVERNOR v2.4: SUSPENSION STATE IS THE GOVERNOR'S OWN MEMORY, NOT psutil's STATUS; RESUME REPEATS UNTIL THE PROCESS RUNS (2026-09-14 02:58, coordinator)

Lane I's L128: a job whose root process blocks in `subprocess.run` reads "running" to psutil
after `suspend()` on Windows, so the v2.3 band re-suspended its tree on every hot tick (18
suspends, 0 resumes on `verify_grad_key_collision`), and because Windows counts suspends per
thread the child then needed 16 `resume()` calls; the audit stalled 28 minutes. v2.4: the
governor keeps its own set of the jobs it has suspended and never suspends one twice; a job is
"suspended" if it is in that set or psutil says stopped; `resume()` repeats until psutil no
longer reports the process stopped (bounded at 32), so a process suspended many times by an
earlier governor is fully released. Restarted at 02:58; the running jobs are adopted (L83). The
last governor change of the sprint; the hygiene list in `s26/agentI_FINDINGS.md` carries the
lesson (the ledger of a supervisor's actions must be the supervisor's own state).

---

## L130 -- docs/FINDINGS.md CORRECTIONS LEDGER UPDATED WITH THE S26 BLOCK (DELIVERABLE 14) (2026-09-14 03:03, coordinator)

A "Sprint 26 additions" table is appended to the front-matter corrections section of
`docs/FINDINGS.md` (after the two structural lessons, before the Index); no sprint body is
edited. Rows: the S25 "no barren plateau" claim scoped to depth 3 (R2); S7-11's lost -0.288
artefact re-measured (L62); S10-4's leak prices reproduced and re-sourced with the 2/60 bound
(L44, L58); the seven-configuration suite's basis named (L28, L29); S16's AMBER-relaxation
finding confirmed and sharpened (L39, L87, L100); the routers extended (L110, L115); and S26's
own corrections R1, R5, L123, L101, L107, L57, L9. If the Adversary's final RETRACTIONS.md adds
an entry before the close, the coordinator appends a matching row. Rule 9 respected: nothing
superseded is deleted anywhere.

---

## L131 -- ROTAMER RELIEF B (tournament 8): RELIEVING THE BUILDER's CHI1 HALVES THE CATASTROPHIC FRACTION (0.535 TO 0.233 ABOVE 1e4) AND LEAVES AN ENERGY THAT RANKS NOTHING (rho -0.002) AND STILL REJECTS HARMFULLY (+0.030 VS ANCHOR, +0.022 VS RANDOM, FOLD CIs EXCLUDING ZERO); CLOSED (2026-09-14, PH)

`s26/ph_relief.py run | report` and `s26/ph_relief_reject.py`; artefacts
`s26/results/ph_relief_run.json` (126/126; job `ph_relief_run`, AMBER, 6890 s, peak RSS 0.313 GB;
gate G1 at 0.0 against `refine_coords` and the cache on 4 members of every target after the
k = 4 fix of addendum 1), `ph_relief_report.json` (B ii, gated), `ph_relief_reject.json` (B iii,
point cloud, 16 draws). Pre-registered in `s26/PREREG_rotamer_relief.md` sections 1 to 7 and
addenda 1 to 2; tournament item 8 (L51). Part A is L23.

THE OPERATOR, native-free: one greedy sweep over every residue with a chi1 (all but G, A, P;
10 of 14 on 1A13), options {builder default, 60, 180, 300} degrees, the lowest genuine
ff14SB/GBn2 single point kept per residue (no minimisation, backbone fixed); 41 single points
and 0.67 s per member. `E_relief <= e_raw` by construction and asserted on 9,450 members.

B(i), THE CENSUS (native-free): fraction of the shipped top-75 above 1e4 kcal/mol, raw
0.535 (SE 0.026) to relieved 0.233 (SE 0.024); above 1e6,
0.311 to 0.057. The registered falsifier ("not below half of the raw fraction, 0.27")
does NOT fire: 0.233 < 0.27. More than half of the catastrophic single points on the shipped
pool were the deterministic builder's chi1 choice, as L23 predicted from the contact classes
(96.8% side-chain-involving). About a quarter of every pool remains above 1e4 after relief:
those are the backbone-owned or multi-rotamer clashes a one-pass greedy sweep cannot reach.

B(ii), DOES THE RELIEVED ENERGY RANK? (ORACLE evaluation, `rr` as the label): NO, and neither
did the raw one. Spearman(E, ORACLE RMSD) within the top-75: raw +0.0000 (SE 0.0202),
relieved -0.0018; in-band (members within 3 A of the best) raw +0.0043, relieved +0.0055.

B(iii), THE REJECT ON THE RELIEVED ENERGY, arm S (reject from the shipped top-75, no refill),
POINT CLOUD on both sides, matched-count random and permuted controls: at 1e4 the relieved
threshold rejects 17.5 of 75 (raw 40.1); at 1e6, 4.3 (raw 23.3). `ST.fmt` verbatim:

      rho(E_relief, ORACLE d) minus rho(E_raw, ORACLE d), whole top-75
        a -0.0018 (med 0.0232)   b 0.0000 (med -0.0079)   n=126
        effect -0.0019   median +0.0044   SE 0.0168   MDE 0.0471   effect/MDE -0.04
        iid  CI95 [-0.0345, +0.0300]
        fold CI95 [-0.0315, +0.0260]   folds same sign 2/5   per-fold 0:-0.038 1:-0.042 2:+0.054 3:+0.002 4:+0.009
        62W/64L/0T   worst degradation +0.5945 (7JGX)   p90 +0.2374   power 0.05  Type-M 21.25
        concentration: drop-top10 +0.0290 vs uniform-effect null p10/p50/p90 +0.0084/+0.0286/+0.0496 -> pctile 0.509
        VERDICT: NOT MEASURED (|effect| 0.0019 <= its own MDE 0.0471, 0.04x)
      in-band rho(E_relief) minus rho(E_raw)
        a 0.0055 (med 0.0347)   b 0.0043 (med -0.0045)   n=126
        effect +0.0012   median +0.0069   SE 0.0167   MDE 0.0467   effect/MDE +0.03
        iid  CI95 [-0.0323, +0.0324]
        fold CI95 [-0.0228, +0.0261]   folds same sign 2/5   per-fold 0:-0.034 1:-0.022 2:+0.049 3:-0.005 4:+0.014
        62W/64L/0T   worst degradation +0.4460 (9BAF)   p90 +0.2436   power 0.05  Type-M 32.72
        concentration: drop-top10 +0.0327 vs uniform-effect null p10/p50/p90 +0.0116/+0.0323/+0.0525 -> pctile 0.509
        VERDICT: NOT MEASURED (|effect| 0.0012 <= its own MDE 0.0467, 0.03x)
      point_cloud [all, n=126] S_relief@1e4 minus anchor
        a 3.0781 (med 2.8640)   b 3.0483 (med 2.8373)   n=126
        effect +0.0297   median +0.0000   SE 0.0144   MDE 0.0403   effect/MDE +0.74
        iid  CI95 [+0.0044, +0.0608]
        fold CI95 [+0.0188, +0.0404]   folds same sign 5/5   per-fold 0:+0.016 1:+0.043 2:+0.022 3:+0.020 4:+0.045
        46W/54L/26T   worst degradation +1.2601 (1U62)   p90 +0.1078   power 0.54  Type-M 1.35
        concentration: drop-top10 +0.0457 vs uniform-effect null p10/p50/p90 +0.0275/+0.0444/+0.0659 -> pctile 0.534
        VERDICT: NOT MEASURED (|effect| 0.0297 <= its own MDE 0.0403, 0.74x)
      point_cloud [all, n=126] S_relief@1e4 minus RANDS@1e4
        a 3.0781 (med 2.8640)   b 3.0561 (med 2.8358)   n=126
        effect +0.0220   median +0.0000   SE 0.0129   MDE 0.0360   effect/MDE +0.61
        iid  CI95 [-0.0017, +0.0484]
        fold CI95 [+0.0106, +0.0318]   folds same sign 5/5   per-fold 0:+0.004 1:+0.037 2:+0.021 3:+0.015 4:+0.032
        46W/54L/26T   worst degradation +1.1635 (1U62)   p90 +0.0985   power 0.40  Type-M 1.56
        concentration: drop-top10 +0.0380 vs uniform-effect null p10/p50/p90 +0.0224/+0.0368/+0.0552 -> pctile 0.540
        VERDICT: NOT MEASURED (|effect| 0.0220 <= its own MDE 0.0360, 0.61x)
      point_cloud [all, n=126] S_relief@1e4 minus S_raw@1e4
        a 3.0781 (med 2.8640)   b 3.1563 (med 3.0432)   n=126
        effect -0.0782   median -0.0166   SE 0.0250   MDE 0.0700   effect/MDE -1.12
        iid  CI95 [-0.1307, -0.0317]
        fold CI95 [-0.1291, -0.0319]   folds same sign 5/5   per-fold 0:-0.027 1:-0.011 2:-0.111 3:-0.173 4:-0.072
        76W/38L/12T   worst degradation +0.9317 (6HVK)   p90 +0.0813   power 0.88  Type-M 1.07
        concentration: drop-top10 -0.0174 vs uniform-effect null p10/p50/p90 -0.0437/-0.0183/+0.0054 -> pctile 0.520
        VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.07x]
      point_cloud [all, n=126] S_relief@1e6 minus anchor
        a 3.0672 (med 2.8422)   b 3.0483 (med 2.8373)   n=126
        effect +0.0188   median +0.0000   SE 0.0114   MDE 0.0320   effect/MDE +0.59
        iid  CI95 [+0.0010, +0.0452]
        fold CI95 [+0.0012, +0.0399]   folds same sign 5/5   per-fold 0:+0.001 1:+0.057 2:+0.002 3:+0.001 4:+0.032
        38W/46L/42T   worst degradation +1.2577 (1U62)   p90 +0.0450   power 0.38  Type-M 1.61
        concentration: drop-top10 +0.0264 vs uniform-effect null p10/p50/p90 +0.0117/+0.0248/+0.0425 -> pctile 0.554
        VERDICT: NOT MEASURED (|effect| 0.0188 <= its own MDE 0.0320, 0.59x)

READING. (1) The relief is real as a census fact: it halves the fraction of the pool the single
point calls catastrophic, so "AMBER measures the builder's rotamers" (L23) is confirmed by
intervention, not only by contact class. (2) It does not make the energy a ranker: rho 0.000 to
-0.002 whole-set, +0.004 to +0.006 in-band, difference 0.04x MDE with 2/5 folds. A single point
with its side chains relieved is still blind to which candidate is right, which is the S25
landscape result (`rho_AMB_rmsd_inband` +0.09 there, on a different band definition) reaching
the same place through a cleaner observable. (3) It does not make the reject helpful: S on the
relieved energy at 1e4 is +0.030 A worse than the anchor (fold CI [+0.019, +0.040], 5/5 folds,
0.74x MDE: the SIGN is measured by the CI, the size is underpowered) and +0.022 worse than
rejecting the same count at random (fold CI [+0.011, +0.032], 5/5, 0.61x MDE, underpowered);
it is less harmful than the raw reject by -0.078 (1.12x MDE, 5/5, Type-M) only because it
rejects a third as many members. At 1e6 it is +0.019 vs the anchor (0.59x) and +0.017 vs random
(0.58x), both underpowered nulls. The falsifier of section 5 (beat the anchor and the matched
control past MDE on the built chain) cannot clear on the point cloud, so the built-chain half
(1.4 h of projection) was NOT run: a point-cloud arm that is worse than the anchor with a fold
CI excluding zero would need the projection to reverse its sign, which no arm in L86 did.
Arm R (refill) was not run, as registered (it needs relieved energies beyond the top-75).

POWER. B(ii): SE 0.017 on the rho difference, MDE 0.047: a ranking gain of 0.05 in rho would
have shown; below that, underpowered. B(iii): SEs 0.010 to 0.014 at 1e4 and 1e6; the arm's harm
vs the anchor at 1e4 is measured in sign (fold CI excludes zero on all of S_relief-vs-anchor,
S_relief-vs-RANDS, S_relief-vs-PERMS) and its size is an upper bound. No positive result, no
replication due; the direction replicates L43 / L86 on a different energy.

DISPOSITION. rotamer_relief B CLOSED: the relief removes the builder's half of the singularity
and leaves an energy that ranks nothing and rejects harmfully, which closes the last route by
which an all-atom single point on retrieved windows could have entered the pipeline. For the
report: "half of what the force field called impossible was the way the side chains had been
placed; once that is fixed, the force field still cannot tell a good backbone from a bad one."

---

## L132 -- C5 COMPLETE: PREDICTING THE COMMON-MODE DIRECTION ON HELD-OUT FOLDS AND SUBTRACTING IT IS NULL (GLOBAL) TO HARMFUL (RIDGE +0.164 A, 1.79x MDE, 5/5 FOLDS); THE ORACLE CEILING IS -1.87 A (DISTANCE SPACE) AND -3.13 A (COORDINATE FRAME); THE RIDGE PREDICTOR IS WORSE THAN A RANDOM MOVE OF ITS OWN SIZE (2026-09-14 03:54, lane P)

`s26/p_c5.py run`, `s26/results/p_c5.json` (complete: true, 126/126; job `p_c5_run`, exit 0, 5,284 s, peak RSS 0.09 GB),
statistics `s26/results/p_c5_stats.json`; pre-registered in `s26/PREREG_C5.md` (alpha grid reduced to 0.5, declared in the
commit of `p_c5.py`). Construction: per target the shipped top-75 cloud c; the common mode ebar = c - t after one rigid
Kabsch fit (ORACLE, training folds only inside the fits). R1 = distance space, five SHELL means of D(c) - D(t), applied
by stress descent initialised at c; R2 = the cloud's own principal-axis frame (native-free; round trip exact, rotation-
invariant to 3e-15 in the selftest). GLOBAL = the training-fold mean correction per length band; RIDGE = nested
leave-fold-out ridge from the 45 B3 features to the five R1 shell means (alpha chosen 100 to 10,000 per fold, i.e.
near-constant predictions); RANDOM = a random direction of the same magnitude through the same operator; ORACLE =
the true correction (the bound). Every arm through `s12.instrument.project`; built chain primary (rebuild basis
3.2126), point cloud carried. Negative = better than the incumbent.

```
  C5 global_R1_a0.5 - incumbent [BUILT CHAIN]
    a 3.2438 (med 2.9603)   b 3.2126 (med 2.9661)   n=126
    effect +0.0312   median +0.0320   SE 0.0185   MDE 0.0519   effect/MDE +0.60
    iid  CI95 [-0.0046, +0.0668]
    fold CI95 [+0.0157, +0.0495]   folds same sign 5/5   per-fold 0:+0.034 1:+0.040 2:+0.062 3:+0.012 4:+0.011
    54W/72L/0T   worst degradation +0.5890 (7BX2)   p90 +0.2540   power 0.39  Type-M 1.58
    concentration: drop-top10 +0.0708 vs uniform-effect null p10/p50/p90 +0.0501/+0.0694/+0.0894 -> pctile 0.539
    VERDICT: NOT MEASURED (|effect| 0.0312 <= its own MDE 0.0519, 0.60x)

  C5 global_R2_a0.5 - incumbent [BUILT CHAIN]
    a 3.2036 (med 2.9685)   b 3.2126 (med 2.9661)   n=126
    effect -0.0090   median -0.0060   SE 0.0140   MDE 0.0392   effect/MDE -0.23
    iid  CI95 [-0.0358, +0.0177]
    fold CI95 [-0.0331, +0.0185]   folds same sign 3/5   per-fold 0:+0.043 1:-0.037 2:+0.007 3:-0.021 4:-0.035
    70W/56L/0T   worst degradation +0.5215 (2MSA)   p90 +0.1934   power 0.10  Type-M 3.75
    concentration: drop-top10 +0.0166 vs uniform-effect null p10/p50/p90 -0.0022/+0.0160/+0.0344 -> pctile 0.516
    VERDICT: NOT MEASURED (|effect| 0.0090 <= its own MDE 0.0392, 0.23x)

  C5 ridge_R1 - incumbent [BUILT CHAIN]
    a 3.3769 (med 3.1985)   b 3.2126 (med 2.9661)   n=126
    effect +0.1643   median +0.1294   SE 0.0328   MDE 0.0920   effect/MDE +1.79
    iid  CI95 [+0.1021, +0.2293]
    fold CI95 [+0.1176, +0.2253]   folds same sign 5/5   per-fold 0:+0.105 1:+0.283 2:+0.186 3:+0.115 4:+0.142
    40W/86L/0T   worst degradation +1.4057 (9BAF)   p90 +0.5774   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.2222 vs uniform-effect null p10/p50/p90 +0.1784/+0.2204/+0.2628 -> pctile 0.527
    VERDICT: WORSE

  C5 random_R1 - incumbent [BUILT CHAIN]
    a 3.3477 (med 3.1257)   b 3.2126 (med 2.9661)   n=126
    effect +0.1351   median +0.1268   SE 0.0363   MDE 0.1018   effect/MDE +1.33
    iid  CI95 [+0.0654, +0.2080]
    fold CI95 [+0.0292, +0.2463]   folds same sign 4/5   per-fold 0:+0.043 1:+0.197 2:+0.329 3:-0.032 4:+0.131
    49W/77L/0T   worst degradation +2.0828 (6MK8)   p90 +0.6066   power 0.96  Type-M 1.02
    concentration: drop-top10 +0.1952 vs uniform-effect null p10/p50/p90 +0.1472/+0.1936/+0.2416 -> pctile 0.515
    VERDICT: WORSE

  C5 random_R2_a0.5 - incumbent [BUILT CHAIN]
    a 3.2249 (med 3.0150)   b 3.2126 (med 2.9661)   n=126
    effect +0.0123   median +0.0034   SE 0.0098   MDE 0.0274   effect/MDE +0.45
    iid  CI95 [-0.0070, +0.0314]
    fold CI95 [-0.0081, +0.0299]   folds same sign 3/5   per-fold 0:+0.023 1:-0.019 2:+0.023 3:-0.014 4:+0.038
    60W/66L/0T   worst degradation +0.3933 (5H1H)   p90 +0.1343   power 0.24  Type-M 2.03
    concentration: drop-top10 +0.0304 vs uniform-effect null p10/p50/p90 +0.0188/+0.0302/+0.0419 -> pctile 0.513
    VERDICT: NOT MEASURED (|effect| 0.0123 <= its own MDE 0.0274, 0.45x)

  C5 oracle_R1 - incumbent [BUILT CHAIN]
    a 1.3377 (med 0.9955)   b 3.2126 (med 2.9661)   n=126
    effect -1.8749   median -1.5004   SE 0.1300   MDE 0.3641   effect/MDE -5.15
    iid  CI95 [-2.1296, -1.6178]
    fold CI95 [-2.0919, -1.6021]   folds same sign 5/5   per-fold 0:-1.331 1:-1.916 2:-1.840 3:-2.171 4:-2.099
    120W/6L/0T   worst degradation +0.2322 (2LER)   p90 -0.3580   power 1.00  Type-M 1.00
    concentration: drop-top10 -1.5777 vs uniform-effect null p10/p50/p90 -1.7452/-1.5851/-1.4257 -> pctile 0.522
    VERDICT: BETTER

  C5 oracle_R2 - incumbent [BUILT CHAIN]
    a 0.0834 (med 0.0616)   b 3.2126 (med 2.9661)   n=126
    effect -3.1293   median -2.8729   SE 0.1534   MDE 0.4297   effect/MDE -7.28
    iid  CI95 [-3.4275, -2.8439]
    fold CI95 [-3.2715, -3.0057]   folds same sign 5/5   per-fold 0:-2.952 1:-3.096 2:-3.018 3:-3.400 4:-3.188
    126W/0L/0T   worst degradation -0.1490 (1S9Z)   p90 -0.9687   power 1.00  Type-M 1.00
    concentration: drop-top10 -2.8210 vs uniform-effect null p10/p50/p90 -3.0202/-2.8272/-2.6380 -> pctile 0.516
    VERDICT: BETTER

  C5 global_R1_a0.5 - incumbent [POINT CLOUD]
    a 3.0843 (med 2.8019)   b 3.0483 (med 2.8373)   n=126
    effect +0.0360   median +0.0513   SE 0.0152   MDE 0.0427   effect/MDE +0.84
    iid  CI95 [+0.0065, +0.0656]
    fold CI95 [+0.0081, +0.0643]   folds same sign 4/5   per-fold 0:+0.048 1:+0.088 2:+0.051 3:+0.006 4:-0.003
    47W/79L/0T   worst degradation +0.4770 (2MSA)   p90 +0.2234   power 0.66  Type-M 1.24
    concentration: drop-top10 +0.0671 vs uniform-effect null p10/p50/p90 +0.0495/+0.0670/+0.0839 -> pctile 0.502
    VERDICT: NOT MEASURED (|effect| 0.0360 <= its own MDE 0.0427, 0.84x)

  C5 global_R2_a0.5 - incumbent [POINT CLOUD]
    a 3.0647 (med 2.8423)   b 3.0483 (med 2.8373)   n=126
    effect +0.0164   median +0.0207   SE 0.0141   MDE 0.0396   effect/MDE +0.41
    iid  CI95 [-0.0112, +0.0438]
    fold CI95 [-0.0172, +0.0549]   folds same sign 3/5   per-fold 0:+0.082 1:+0.028 2:+0.032 3:-0.008 4:-0.042
    55W/71L/0T   worst degradation +0.5433 (2MSA)   p90 +0.2299   power 0.21  Type-M 2.18
    concentration: drop-top10 +0.0409 vs uniform-effect null p10/p50/p90 +0.0224/+0.0404/+0.0594 -> pctile 0.515
    VERDICT: NOT MEASURED (|effect| 0.0164 <= its own MDE 0.0396, 0.41x)

  C5 ridge_R1 - incumbent [POINT CLOUD]
    a 3.1975 (med 3.0171)   b 3.0483 (med 2.8373)   n=126
    effect +0.1492   median +0.1359   SE 0.0278   MDE 0.0780   effect/MDE +1.91
    iid  CI95 [+0.0959, +0.2014]
    fold CI95 [+0.0947, +0.2292]   folds same sign 5/5   per-fold 0:+0.102 1:+0.310 2:+0.152 3:+0.071 4:+0.123
    36W/90L/0T   worst degradation +1.2111 (6BJF)   p90 +0.5080   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.2024 vs uniform-effect null p10/p50/p90 +0.1674/+0.2014/+0.2357 -> pctile 0.515
    VERDICT: WORSE

  C5 random_R2_a0.5 - incumbent [POINT CLOUD]
    a 3.0673 (med 2.8066)   b 3.0483 (med 2.8373)   n=126
    effect +0.0190   median +0.0136   SE 0.0078   MDE 0.0219   effect/MDE +0.87
    iid  CI95 [+0.0042, +0.0346]
    fold CI95 [+0.0030, +0.0350]   folds same sign 4/5   per-fold 0:+0.040 1:+0.040 2:+0.017 3:-0.013 4:+0.013
    50W/76L/0T   worst degradation +0.3575 (5H1H)   p90 +0.1132   power 0.68  Type-M 1.22
    concentration: drop-top10 +0.0312 vs uniform-effect null p10/p50/p90 +0.0207/+0.0307/+0.0417 -> pctile 0.521
    VERDICT: NOT MEASURED (|effect| 0.0190 <= its own MDE 0.0219, 0.87x)

  C5 oracle_R2 - incumbent [POINT CLOUD]
    a 0.0000 (med 0.0000)   b 3.0483 (med 2.8373)   n=126
    effect -3.0483   median -2.8373   SE 0.1467   MDE 0.4110   effect/MDE -7.42
    iid  CI95 [-3.3418, -2.7668]
    fold CI95 [-3.2083, -2.9005]   folds same sign 5/5   per-fold 0:-2.772 1:-3.018 2:-2.979 3:-3.314 4:-3.156
    126W/0L/0T   worst degradation -0.1959 (1S9Z)   p90 -1.1540   power 1.00  Type-M 1.00
    concentration: drop-top10 -2.7422 vs uniform-effect null p10/p50/p90 -2.9345/-2.7430/-2.5646 -> pctile 0.501
    VERDICT: BETTER
```

Reading. (1) The ceiling is large: subtracting the TRUE common mode gives 1.34 A (distance space, -1.87, 120W/6L)
and 0.08 A (the coordinate frame, -3.13, 126W/0L; the R2 oracle is the native by construction, so its number is
the size of the mode plus the projection's cost, not an achievable target). (2) The two achievable arms are null:
GLOBAL R2 -0.009 (0.23x MDE, 3/5 folds, 70W/56L) and GLOBAL R1 +0.031 (0.60x, fold CI above zero at 5/5 but under
the MDE gate). (3) The learned arm is HARMFUL: RIDGE R1 +0.164 A [fold +0.118, +0.225], 1.79x MDE, 5/5 folds,
40W/86L, and its magnitude-matched random control is +0.135 [+0.029, +0.246], 1.33x, so the predicted shell profile
is indistinguishable from a random shell profile of the same size (+0.029, inside the noise). The nested alphas
went to the top of the grid on 3 of 5 folds: the ridge found nothing to fit and emitted a near-constant
correction, which is the GLOBAL arm plus noise. (4) Power: the design resolves 0.04 to 0.05 A on the GLOBAL arms
and 0.09 A on the ridge arm; a real gain of 0.05 A would have shown as a fold CI below zero on GLOBAL R2 and did not.
Verdict: the pre-registered falsifier fires for both GLOBAL and RIDGE; C5 is refuted at one agent-day, as S16,
S19 L14 ("the estimator is made of the bias") and S24 L7 predicted, now with the ceiling measured beside it on the
same operator: the common mode is 1.9 to 3.1 A of built-chain error, and nothing native-free in these two
representations touches it.

---

## L133 -- THE raw RUNG (1280-d ESM-2 WITH A LEARNED PROJECTION) IS NOT RUN: 4 OF 5 FOLDS TRAINED, FOLD 4 UNTRAINED, NO EVALUATION; MEMORY AND TIME (2026-09-14 03:54, lane P)

`s26/models/p_ladder/raw_fold{0,1,2,3}_s0.pt` exist (jobs `p_train_raw_f1`, `p_train_raw_fold1/2/3`: 4,477 to 4,900 s
per fold at peak RSS 1.93 to 1.96 GB; the first attempt was terminated at 19:11 with the governor dead, L42, and
fold 3 was queued 40 min behind the job cap). At the 04:30 close fold 4 is untrained and the rung cannot be
evaluated (a rung needs all five leave-fold-out models). No further launch: a fifth fold is 80 min of a 2 GB job
and the evaluation another 10, past the close. What the record says instead: raw's inputs are bracketed by pca32f
(32 components, +0.018 built chain, 0.13x MDE, L67) and pca128 (128 components, +0.076, 0.43x, L72), both null,
and S7-11 measured raw on selection as indistinguishable from pca32 and pca128. The rung is recorded as NOT RUN,
not as a result; its four checkpoints stay on disk for the next sprint.

---

## L134 -- LANE PR: C5's CLOSURE (L132, PROPOSAL_C.md ADDENDUM 2) AND THE raw RUNG's STATUS (L133) APPLIED TO SLIDE 10 AND ITS NOTES; THE DECK IS FINAL WITH NOTHING LIVE OR PENDING; 520 REGISTERED NUMBERS; VERIFICATION PASS (2026-09-14 03:55, PR)

Rebuilt after `s26/PROPOSAL_C.md` addendum 2 landed at 03:54 (`python s26/pr_build_deck.py --no-figures; python
s26/pr_changes.py`). Slide 10's first bullet and its notes now state C5 as CLOSED at one agent-day, every number read
from `s26/results/p_c5.json :: summary` (complete: true, 126/126): GLOBAL in the coordinate frame -0.0090 A (0.23x MDE,
fold CI [-0.033, +0.019], NOT MEASURED), GLOBAL in distance space +0.0312 (0.60x, NOT MEASURED), RIDGE in distance
space +0.1643 A (1.79x MDE, fold CI [+0.118, +0.225], 5/5 folds, WORSE), the random correction of the same size
+0.1351 (1.33x), the ORACLE ceilings -1.875 A (distance space) and -3.129 A (coordinate frame); the raw rung is stated
as NOT RUN (4 of 5 fold models on disk at build time, no evaluation, L133). Nothing else on the deck changed.
Verification (`s26/pr_verify.txt`, pasted into `s26/PRESENTATION_CHANGES.md`): 11 slides; 0 U+2014, 0 U+2013; 0
banned words; spoken words 249 / 248 / 249 on slides 8 / 9 / 10; title, body and notes on every slide. Registry 520
tokens (`s26/pr_values.json`). `s26/agentPR_FINDINGS.md` section 7 carries the change-log row. The deck is final:
no slide, note or verdict is live, DRAFT or PENDING.

---

## L135 -- DELIVERABLES CHECK (FINAL PASS, 04:00 TO 04:12): 13 OF 14 ROWS OK, THE REPORT ROW RE-CHECKED AT 04:40; RULE 13 CLEAN ON EVERY PRESENTER-FACING FILE AND THE DECK DUMP; RETRACTIONS R6-R10 ADDED; L131 AND L132 STAND (2026-09-14, A)

`s26/DELIVERABLES_CHECK.md` (final). Rows: (1) governor + log OK (v2 to v2.4 ledgered,
log continuous); (2) EXAMINATION + AUDIT OK; (3) BRIEF OK; (4) PREREG_* OK (29, each before
compute); (5) PROPOSAL_A / B / C + B_REPLACEMENT OK (all FINAL; C's MATERIAL document defect
fixed by L124 and C5 final in addendum 2; B's rho -0.638 cites `p_b3.json` but lives in
`agentP_FINDINGS.md` section 10, MINOR, adopted in L124); (6) IDEA_* (18) + TOURNAMENT OK;
(7) agent*_FINDINGS OK, all eight lanes; (8) LEDGER (137 entries) + RETRACTIONS OK: R6 (L123,
the strain sentence), R7 (L57), R8 (L101), R9 (L107), R10 (L9) added at 04:05 beside R1-R5,
with the prior-sprint disposition table (the S25 plateau scoped to depth 3; the S13 Pauli mean
weights not re-derived and kept off the slides; S7-11's -0.288 replaced by L62's -0.330 /
-0.208; S10-4 re-derived; the suite's basis; S16 sharpened; the routers extended; L38 scoped by
L89); (9) repository fixes OK (items 1-5, defects 6a-6d, cache_amber tracked, the opt-in tier
11/11, the frozen rebuild 2016/2016); (10) the deck OK (`vqe_research_overview.pptx` built
03:55; `PRESENTATION_CHANGES.md` and `pr_notes.md` map every number; the Adversary's
qualifiers are on the slides, L46/L54/L55/L70/L75/L122/L123); (11) TRAINABILITY_PAPER_OUTLINE
OK, one MINOR (put the literal "at depth 3" in the F5 caption); (12) REPORT PARTIAL at 04:00
(3,423 lines at `439ce5f8`, rule-13 clean; lane E's final pass 04:15-04:45; re-checked at
04:40, next entry); (13) STATUS OK; (14) docs/FINDINGS.md corrections block OK at line 117
(L130), with the R6-R10 rows owed by the coordinator's own clause.

Rule 13 at 04:02 (banned words; U+2014; U+2013): 0 / 0 / 0 on `REPORT.md`,
`PRESENTATION_CHANGES.md`, `pr_notes.md`, the deck dump `pr_verify_dump.txt`, the four
proposal files, the paper outline, TOURNAMENT, RETRACTIONS, EXAMINATION, EXAMINATION_AUDIT,
BRIEF, agentA_FINDINGS, and the S26 block of `docs/FINDINGS.md`.

L131 (rotamer relief B): the operator is native-free (a greedy chi1 sweep on the builder's
side chains); the census clears its registered falsifier (0.535 to 0.233 above 1e4, below half
the raw fraction); the relieved energy ranks nothing (rho difference -0.002, 0.04x MDE, in-band
+0.001) and its reject is +0.030 against the anchor at 0.74x MDE with the fold CI above zero
5/5 (Type-M zone: the sign matches L43 / L86, the size is not a result). Closed as an accuracy
step with the power stated. STANDS; the singularity is the builder's (L23) and relieving it
changes nothing the record cares about.

L132 (C5): predictors fitted on training folds (nested), a random direction of matched
magnitude through the same operator as the control, the ORACLE ceilings labelled; GLOBAL R1
+0.031 (0.60x, fold CI above zero 5/5: suggestive harm, not measured), R2 -0.009 (0.23x), RIDGE
+0.164 (1.79x, WORSE, 5/5) and +0.135 (1.33x, WORSE, 4/5); the alpha grid's reduction to 0.5 was
declared. Leakage: none in the predictors; the ceilings read the native and say so. STANDS; the
record's prediction (S16, S19 L14, S24 L7) measured with a falsifier, and Proposal C's addendum
2 carries it.

Not done by this lane before the close: the checks of L115 (C4, eighteen routers) and
L112 / L116 (the delivery file) were displaced by the proposals, L117, the spread control and
this pass; both entries are nulls or reproductions (C4 null-to-harmful on all eighteen; the
delivery file is the production emission, cross-checked by L116), so nothing positive stands
unchecked. Recorded in `s26/agentA_FINDINGS.md` "what I did not do".

---

## L136 -- ADVERSARY CHECK OF L115 AND L112 / L116 (C4's EIGHTEEN ROUTERS; THE C3 DELIVERY FILE): BOTH STAND; THE "NOT DONE" LINE OF L135 IS CLOSED (2026-09-14, A)

L115 (C4). Nested leave-fold-out ridge (the router is fitted on training folds only; the m* and
s* labels are ORACLE per-target optima used as labels, never at inference); 200-draw label
permutation nulls per router; the anchor reproduces `agg_surface` and `errdecomp` at 0.00e+00;
the S22 feature set runs through the same harness as the control and reproduces its recorded
~0 (+0.020 to +0.047); both bases stated (point cloud where m* and s* live, built chain
carried). Twelve m routers: eleven harmful in sign, none clears its MDE (largest 0.91x, worse
than 96% of its null); the one negative number (retrieval-score entropy, router B, -0.0012,
0.02x MDE) is read as "not harmful", which is right, since the permutation null's own mean is
+0.019 (any routing away from m = 75 costs about 0.02 A, so p_perm 0.010 against that null is
not a gain). Six s routers predict s* with the wrong sign (rho -0.10 to -0.21) and cost +0.007
to +0.027 against s = 1. No order statistic is taken (each routed value is the ridge's own
pick, not a minimum). Power stated per router. S22 L7 / S23 L7 extended to eighteen
constructions; S22 L10's bound stands. STANDS.

L112 / L116 (the delivery file). `s26/results/p_best_rung_chains.json` is the production
emission (126/126, complete, mean `rmsd_arm` 3.2148, torsions in radians, unwrapped per L114),
delivered because no C2 rung beat the shipped prior on the built chain (the nine fold CIs listed
in L112 all include or exceed zero on the harmful side) and because PH's C3 stage 1 already
relaxed that chain, so stage 2 is L87's replication. The rebuild-basis re-projection is saved
separately (`p_best_rung_chains_rebuild_basis.json`, mean 3.2126, L57) as the cross-check; its
six-minute overwrite of the delivery file was caught and the file restored from commit
6ed3b367, with the per-target difference recorded in `s26/agentP_FINDINGS.md` section 11.
STANDS.

---


## L137 -- LANE W CLOSING ENTRY: WHAT RAN, WHAT DID NOT, AND WHAT EACH OPEN ITEM WOULD NEED (2026-09-14 04:05, W)

Everything of lane W's is committed (last commit `7b2e7e77` before this entry); nothing of lane
W's is registered with the governor; `s26/agentW_FINDINGS.md` sections 0 to 6 are final.

**Ran and posted (12 ledger entries, every number with its artefact path, every RMSD line with
its basis, every null with its power):**

    L30   the 2/60 proxy bound pre-registered; native-free census (two of the four self-windows never reach the emission; a BLOSUM tie at every pool boundary)
    L44   the bound measured: direct dev-proxy price 0.002 A (later 0.008, L108), own-native envelope 0.028 A (built chain) / 0.048 (selection) / 0.023 (paired gain); class MINOR by the envelope; C27's +0.0004 re-derived exactly; a fold model that saw the target's native emits a chain 0.70 A nearer (2.83x MDE, 5/5)
    L52   tournament item 2: the same sequence in a different deposit sits 2.91 A from the native (22 pairs; F5 holds)
    L58   L55's caveats applied: the artefact carries the gain envelope row and the per-target readings (worst target 0.151 A under A2)
    L64   tournament item 4: the tie-break noise floor, 0.004 A on the 126-mean, 0.024 A paired MDE between two conventions; every hundredths-level recorded effect inside it, with the scope stated (L71 caveat)
    L84   tournament item 6 (the orphan, the mandatory test-time-ensembling direction): fixed-K ensembling refuted with power (-0.0005 A, 0.02x MDE; a gain of 0.028 A or more excluded)
    L85   tournament item 9: provenance census (73% fragment windows, 0.6% whole peptides) and the ORACLE class contrast (0.42 A in the pool, 0.09 to 0.13 inside the top-75)
    L90   partial recall gradient: none at the MDE; I_long rho -0.205 suggestive and confounded with retrieval
    L91   the own-native models on the S24 ladder: the direction-discounted currency under-prices the trained operator by 0.35 A (1.5x MDE)
    L105  tournament item 10 (the second orphan): AMBER as a prior partner; the leave-fold-out choice declines to mix on all five folds; every cell worse on the cloud
    L108  Part B on all four carrier-out models: F3 falsified by 2P5H in the harmful direction; direct bound 0.008 A, IMMATERIAL; class unchanged
    L109  provenance readout H_P3: not measured; a gain of 0.044 A or more excluded

Governed jobs run by this lane: 32 (`s26/jobs_done/w_*.json`), 31 exit 0 and one exit 1 (the first synthetic test run, `w_selfcopy_test`, which caught the `sel`-bound error and was fixed before any real run; L30); largest peak RSS
1.250 GB (the four Part B retrains at 1.244 to 1.250; every other job under 0.32 GB); longest
wall 9,321 s (the provenance readout under a six-job load). No benchmark file, sequence, name,
native or RMSD was read at any point; the only benchmark-derived inputs were the record's count
2/60 and the mechanism of S24 L4.

**Did not run, with the reason and the time each would need:**

1. Part B's five remaining control-out models (8TXS/f0, 8TXS/f2, 8T63/f2, 9BAF/f4, 8T63/f4;
   one of six, 9BAF/f0, is built) and therefore F3's CONTROL CLAUSE (whether the carrier-out
   change exceeds the change from removing an unrelated dev chain of the same fold). Reason: the
   coordinator's starvation ruling (L111) withdrew the queued waiter at 01:36 and the hold was
   lifted at 04:05, 25 minutes before the close. Time needed: 5 x 15 to 24 min of training
   (measured 893 to 1,440 s per model at est-ram 1.4 GB, one at a time), then `w_selfcopy.py
   posterior` + `w_endpoint_report.py` + `w_bound_addendum.py` (about 3 min), then a ledger
   addendum: about 2 h of wall on tonight's box. What it would decide: whether 2P5H's -0.246 A
   (the one large value) is larger than the noise of removing any one chain from the corpus;
   until then F3's control clause is OPEN and L108 says so.
2. A second-seed retrain of 2P5J out of fold 4 (the replication of the one large Part B value):
   about 18 min plus the 3-min re-run. Not run for the same reason.
3. The amber prior partner's permuted-AMBER control at a non-identity cell: moot, the
   leave-fold-out choice never left (0, 0); no time was needed and none was spent.
4. The tie-break floor at 16 draws (PREREG fork; about 50 min) and on the AMBER-relaxed basis;
   the ensembling control at a second seed; the provenance readout at 8 permutation draws:
   each was the pre-declared replication or extension for a POSITIVE, and none of the three
   was positive.
5. The recall gradient's ORACLE mechanism check and a retrieval-only covariate to partial out
   (about 5 min): the PREREG gated the check on a gradient at the MDE, which did not exist; the
   addendum says how to run it next time.
6. The triangle bound on the benchmark itself (needs no native, no RMSD, but the two benchmark
   sequences and their universes): closed by L18b / L29 for this sprint; a coordinator's option.
7. `attn` and coherence_penalised_training (deferred on memory by L51 / L77): never reached
   lane W; nothing was started.

**The two sentences this lane leaves for the report and the deck**, each with its artefact:
"The 2/60 benchmark self-copy leak is bounded from the dev proxy at 0.008 A (direct, both
channels, `s26/results/w_selfcopy_bound.json :: signed_bounds_gated/both_removed4`) and at
0.028 A mean-CI / 0.151 A worst target by the own-native envelope under assumption A2 (the same
file, `C_envelope_per_target`); it is MINOR under every reading on the built chain and the paired
gain and cannot move the benchmark verdict either way; where the leak is measurable it HURTS
the leaked target." And: "The pipeline's own convention noise is 0.004 A on the 126-mean and
0.024 A as the paired MDE between two tie-breaks of the pool boundary
(`s26/results/w_tiebreak_report.json`); a hundredths-level effect is real only as a paired
contrast with the tie-break held fixed."

## L138 -- LANE Q, A2 PER GROWTH STEP ON THE 126 A1 RECORDS: THE ALGEBRA OF THE GROWN SET TRACKS THE APPENDED STRINGS, NOT WHAT THE CIRCUIT DOES (dim 7 EXACTLY WHERE NOTHING WAS APPENDED, UP TO 530 WHERE INERT MULTI-QUBIT STRINGS WERE); THE alpha = 0.25 L2 SETS REACH A MEDIAN 1025 OF 8128 AT P = 21; ALL FOUR ADDENDUM-2 PREDICTIONS HELD (2026-09-14, lane Q)

`s26/q_dla_a1.py` -> `s26/results/q_dla_a1.json` (126 records, every ADAPT run's exact Lie
closure at every growth step) and `s26/figures/a2_dla_grown_ladder.png`. Pre-registered as
`s26/PREREG_A2.md` addendum 2 before the run. Property measurement: reads only the operator
lists in `s26/results/a1/<pdb>.json`; no native, no score. The job's `jobs_done` record is
absent: it was relaunched by the coordinator at 00:10 (L92), computed all 126 records (file
provenance 00:14 local), and its wrapper lost the child before the final summary write; the
summary was recomputed from the stored ladders in-process (seconds, no closures) at 04:03.

    per (pool, re-optimiser, cell): final dim of the closure, median [min, max]; targets at dim 7;
    targets with nothing appended; targets with a multi-qubit string appended; targets reaching so(128) = 8128
      V  Adam-best   alpha=1 (78)       7 [7, 82]      dim7 48   nothing 0    multi 30   so128 0
      V  Adam-best   alpha=0.25 (48)   16 [7, 289]     dim7 12   nothing 0    multi 36   so128 0
      L2 Adam-best   alpha=1 (78)      11 [7, 139]     dim7 15   nothing 0    multi 63   so128 0
      L2 Adam-best   alpha=0.25 (48) 1025 [513, 2017]  dim7 0    nothing 0    multi 48   so128 0
      V  L-BFGS-B    alpha=1 (78)      58 [7, 530]     dim7 18   nothing 18   multi 60   so128 0
      V  L-BFGS-B    alpha=0.25 (48)   16 [9, 161]     dim7 0    nothing 0    multi 48   so128 0
      L2 L-BFGS-B    alpha=1 (78)      37 [7, 513]     dim7 10   nothing 10   multi 68   so128 0
      L2 L-BFGS-B    alpha=0.25 (48) 1025 [258, 2017]  dim7 0    nothing 0    multi 48   so128 0

Predictions: P2a held (under L-BFGS at alpha = 1, dim = 7 on exactly the 18 (V) and 10 (L2)
targets where nothing was appended, above 7 on every other); P2b held (alpha = 0.25, L2:
median 1025 at P = 21, inside [100, 2000]; no target reaches 8128); P2c held (pool V never
above 530 against its whole-pool 2080, L27); P2d held (Adam-best at alpha = 1: dim 7 on
exactly the 48 = 78 - 30 (V) and 15 = 78 - 63 (L2) targets without a multi-qubit string).

Reading. The dynamical Lie algebra is a property of the generator set. On the 78 alpha = 1
targets the grown sets whose appended strings are inert (L75: worth <= 1.2e-4 nats, angles
<= 0.018 rad, the state a product to KL <= 4.1e-4) have closures of up to 530 dimensions,
while the fixed ansatz at the same P has 8128 and does no more with it. For a grown circuit
the DLA dimension does not track what the circuit does; it counts what could be done if the
appended angles were not near zero. The only sets with a large algebra AND non-trivial
angles are the alpha = 0.25 L2 sets (13.7 distinct 2-local strings, 1025 of 8128), and A1
measured their endpoint at +0.0000 A against the fixed circuit on those 48 targets (L68).
This is the per-growth-step item the prompt asked for under A2; it adds to
`a2_dla_dimension.png` the real-target ladders (figure `a2_dla_grown_ladder.png`) and closes
the A2 deliverable.

## L139 -- LANE Q, A1 REPLICATION (SEED 1, REVERSED FOLD ORDER): THE SEED-0 "NULL" IS NOT STABLE ACROSS SEEDS. AT SEED 1 THE TWO PRIMARIES ARE -0.045 / -0.052 A AT 0.71x / 0.79x MDE WITH FOLD CIs EXCLUDING ZERO (TYPE-M ZONE); THE SEED-0 POINTS LIE INSIDE THE SEED-1 CIs; THE DIFFERENCE IS THE COMPARATOR'S SEED SENSITIVITY (FIXED 3.228 -> 3.261 A) NOT THE ADAPT ARMS' (3.214 -> 3.216); A1 IS NOT MEASURED ON EITHER SEED; THE VERDICT STANDS (2026-09-14 04:12, lane Q)

Pre-registered as `s26/PREREG_A1.md` addendum 2 (before the run). `s26/q_adapt.py --build --tag
a1s1 --seed 1 --order fold_rev --fixed-iters 50 --optimisers adam_best`, one governed process:
`s26/jobs_done/a1s1_build.json`; label `a1s1_label.json` (25 s, 0.396 GB); stats
`s26/results/a1s1_stats.json`, `s26/logs/a1s1_stats.log`; records `s26/results/a1s1/<pdb>.json`
(126). The bit-for-bit check against the seed-0 production cache holds for the classical arm
(`ca` 0.0 on 126 of 126) and, as it must, not for the seed-1 quantum arm (`q_ca` differs; the
check is seed-0-specific by construction). Basis: built chain on both sides. Seed 1 changes the
fixed circuit's initial angles, ADAPT's initial RY layer and the random controls' draws; the
reversed fold order changes nothing numerical (each target is independent).

### The two PRIMARY contrasts at seed 1, verbatim

  rmsd_q_synth: adaptL2_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2161 (med 3.0983)   b 3.2610 (med 3.2022)   n=126
    effect -0.0449   median -0.0121   SE 0.0227   MDE 0.0635   effect/MDE -0.71
    iid  CI95 [-0.0935, -0.0033]
    fold CI95 [-0.0960, -0.0037]   folds same sign 4/5   per-fold 0:-0.145 1:+0.009 2:-0.009 3:-0.016 4:-0.054
    72W/54L/0T   worst degradation +0.5512 (1NIZ)   p90 +0.1914   power 0.51  Type-M 1.40
    concentration: drop-top10 +0.0069 vs uniform-effect null p10/p50/p90 -0.0139/+0.0059/+0.0250 -> pctile 0.527
    VERDICT: NOT MEASURED (|effect| 0.0449 <= its own MDE 0.0635, 0.71x)

  rmsd_q_synth: adaptV_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2087 (med 3.0874)   b 3.2610 (med 3.2022)   n=126
    effect -0.0523   median -0.0151   SE 0.0235   MDE 0.0659   effect/MDE -0.79
    iid  CI95 [-0.1017, -0.0107]
    fold CI95 [-0.1087, -0.0054]   folds same sign 4/5   per-fold 0:-0.156 1:-0.051 2:+0.026 3:-0.030 4:-0.050
    77W/49L/0T   worst degradation +0.4790 (1NIZ)   p90 +0.2181   power 0.60  Type-M 1.29
    concentration: drop-top10 +0.0019 vs uniform-effect null p10/p50/p90 -0.0208/+0.0005/+0.0212 -> pctile 0.539
    VERDICT: NOT MEASURED (|effect| 0.0523 <= its own MDE 0.0659, 0.79x)

### Beside seed 0

    built chain, arm - fixed_zrank_it50        seed 0 (L68)                         seed 1
      adaptL2 P21                              -0.0138 (0.23x) fold [-0.071,+0.044]   -0.0449 (0.71x) fold [-0.096,-0.004] 4/5  72W/54L
      adaptV  P21                              -0.0222 (0.36x) fold [-0.085,+0.042]   -0.0523 (0.79x) fold [-0.109,-0.005] 4/5  77W/49L
      adaptL2 P7 (the trained RY layer)        -0.0128 (0.21x)                        -0.0498 (0.76x) fold [-0.104,-0.000]
      gibbs_T (exact optimum, no circuit)      +0.0088 (0.11x)                        -0.0241 (0.28x) fold [-0.095,+0.057]
      uniform128                               +0.0135 (0.14x)                        -0.0194 (0.17x)
      randH_fixed / randH_adaptL2              +0.0230 / +0.0337                      +0.0696 (0.51x) / +0.0045
    arm means, built chain                     seed 0            seed 1
      fixed_zrank_it50 (the comparator)        3.2280            3.2610      (+0.033 A with the seed)
      adaptL2_adam_best_zrank_P21              3.2142            3.2161      (+0.002)
      adaptV_adam_best_zrank_P21               3.2059            3.2087      (+0.003)
      rmsd_arm (classical top-75, seed-free)   3.2148            3.2148
    selection readout, seed 1: L2 P21 -0.0561 (0.65x) fold [-0.143,-0.008] 5/5; V P21 -0.0615 (0.69x) 5/5.

### The registered falsifier FIRED, and what the records say it means

Addendum 2 predicted both primaries inside +-0.5x MDE with fold CIs spanning zero; at seed 1
they are at 0.71x and 0.79x with fold CIs excluding zero. The second half of the prediction
held: the seed-0 points (-0.0138, -0.0222) lie inside the seed-1 iid CIs ([-0.094, -0.003] and
[-0.102, -0.011]). Neither seed clears its MDE (0.059 to 0.066 A); seed 0 is below 0.5x and
seed 1 is in the Type-M zone, so by the standing rule A1 is NOT MEASURED on either seed, and
L68's "null at the registered threshold" is a seed-0 statement that does not hold at seed 1.
By S20 L-B's rule (fewer than 4 seeds is NOT MEASURED for a variational arm; ansatz-seed
sensitivity 0.200 A within-target) two seeds do not settle it either way.

The mechanism is in the arm means: the ADAPT arms are seed-stable (they converge to the same
product Gibbs state at alpha = 1 whatever the start, L68 1.3) and moved 0.002 to 0.003 A; the
fixed 21-parameter circuit, which stops 0.90 nats short of that optimum in a seed-dependent
place, moved +0.033 A. The contrast grew because the comparator landed worse at seed 1, not
because ADAPT emitted a better structure. Direction consistent on both seeds (24 of 24 ADAPT
arms negative over the two seeds); magnitude 0.2x to 0.8x MDE and set by the fixed circuit's
initialisation.

What changes: L68's headline "null at the registered threshold" becomes "NOT MEASURED on two
seeds: 0.23x / 0.36x at seed 0, 0.71x / 0.79x at seed 1, direction consistent, magnitude
governed by the comparator's seed". `s26/PROPOSAL_A.md` gets addendum 3 with that sentence.
The verdict REPLACE stands: nothing here says a grown ansatz emits a measurably better
structure; it says the deployed fixed circuit is a seed-sensitive under-optimiser of a
product-state target, which the 7-parameter RY layer (P7, -0.050 A at seed 1, same as P21)
already reaches. A four-seed run would be the next step and is not run tonight.

## L140 -- ADVERSARY CHECK OF L139 AND L138: "NOT MEASURED ON EITHER SEED" IS THE RIGHT READING; REPLACE STANDS ON SEED-INDEPENDENT FACTS; R11 ENTERED FOR L68's "NULL"; ONE QUALIFIER FOR SLIDE 8 (AND 5), ONE SENTENCE FOR SLIDE 4 (2026-09-14, A)

L139 (A1, seed 1, reversed fold order; `s26/results/a1s1_stats.json`, `a1s1/<pdb>.json`).
1. The reading. Seed 0: 0.23x / 0.36x MDE, fold CIs span zero -- UNDERPOWERED (below 0.7x).
   Seed 1: -0.045 / -0.052 A at 0.71x / 0.79x, fold CIs [-0.096, -0.004] and [-0.109, -0.005],
   4/5 folds -- the TYPE-M ZONE (0.7 to 1.3x), which the contract says is not a result: the
   sign is measured, the magnitude is inflated (Type-M 1.40 / 1.29). "NOT MEASURED on either
   seed" is therefore the contract's own vocabulary and is correct; the precise pair of words is
   "underpowered at seed 0, Type-M at seed 1". L68's "null at the registered threshold"
   (PREREG_A1's +-0.5x MDE) was true at seed 0 and does not hold at seed 1; lane Q says so and
   corrects PROPOSAL_A section 3 by addendum 3. Entered as R11 (a scope correction, not a
   withdrawn number).
2. Two seeds are not two replicates of one effect. The ADAPT arms are seed-stable (3.214 to
   3.216, they converge to the same product Gibbs state), so the two seed contrasts are ONE
   ADAPT value against TWO draws of the fixed comparator (3.228 / 3.261); the between-seed
   change of the contrast (+0.031 / +0.030) is the comparator's initialisation, exactly as L139
   diagnoses. So "24 of 24 arms negative over two seeds" is about two draws, not twenty-four
   (the twelve arms within a seed correlate at 0.955, L70), and no pooling across seeds is
   licensed. By S20 L-B's rule a variational arm needs four seeds before its endpoint is
   called; two do not settle it and L139 says so.
3. The comparator's seed variance is a fact about the deployed selector worth a qualifier:
   the fixed 21-parameter circuit's readout moves 0.033 A between two seeds (3.228 -> 3.261)
   while the production top-75 arm is seed-free (3.2148 on both, the selector is off in
   production). It belongs on slide 8 (and on slide 5 if it quotes the fixed circuit's 3.228 or
   the readout insensitivity), not on slide 4, whose 3.2148 does not run the selector.
4. The verdict. REPLACE stands: it rests on seed-independent, exact facts (the optimum is a
   product state on every target, L68 / L75; the DLA is the full so(2^n) from depth 2, L27 /
   L45; the grown circuits give no width-scaling argument, L35 / L119; the appended operators
   are inert, L75 / L138) and on an endpoint that is not measured on either seed. What seed 1
   adds is a sentence about the DEPLOYED circuit, not about growth: a better-optimised state of
   the same product target (reached by the 7-parameter RY layer, P7 -0.050 at seed 1, the same
   as P21, and by the exact Gibbs state, gibbs_T -0.024) sits 0.02 to 0.05 A nearer at 0.3 to
   0.8x MDE, because the fixed circuit is a seed-sensitive under-optimiser. That supports
   "replace", not "grow".
Verdict: STANDS. Owner's correction (addendum 3) accepted.

Line to lane PR (a qualifier must reach a slide): on slide 8, beside the A1 numbers: "A1 is
not measured on two seeds (0.2x to 0.8x of what the comparison resolves); the grown circuits
are seed-stable (3.214 / 3.216 A) and the deployed fixed circuit moves 0.033 A between seeds
(3.228 / 3.261), so the size of the contrast is set by the deployed circuit's seed, not by
growth." On slide 5, if the fixed circuit's readout is quoted: "on two seeds its built chain
is 3.228 and 3.261 A". On slide 4: "3.2148 does not run the selector and is seed-free."
Sources: `s26/results/a1_stats.json`, `a1s1_stats.json`, ledger L68, L139.

L138 (A2 per growth step; `s26/results/q_dla_a1.json`, 126 records). Property measurement
over the stored operator lists, no native, no score; the four pre-registered predictions
(addendum 2) held exactly: dim 7 on precisely the targets where nothing was appended (18 V, 10
L2 under L-BFGS at alpha = 1; 48 and 15 under Adam-best), up to 530 where inert multi-qubit
strings were appended, median 1025 at P = 21 for L2 at alpha = 0.25, never so(128). The reading
is right and matters for the paper: the DLA is a property of the generator set and counts what
could be done, not what the circuit does; for a grown circuit with near-zero appended angles it
over-states expressivity, and the only large-algebra sets with real angles (alpha = 0.25, L2)
measured +0.0000 A at the endpoint (L68). Provenance irregularity, stated by Q: the job's
`jobs_done` record is absent (the wrapper lost the child after the records were written, L92)
and the summary was recomputed in-process from the stored per-target ladders; the 126 record
files carry their own provenance (00:14). Verdict: STANDS WITH THAT NOTE.

---

## L141 -- EVERY STANDALONE AUDIT UNDER verify/ RE-RUN UNDER THE GOVERNOR WITH ITS TRACKED JSON UNTOUCHED: 20 OF 20 RAN; THE SCIENCE REPRODUCES IN EVERY AUDIT THAT HAS A COMPARABLE TRACKED RECORD; THE DIFFERENCES ARE TIMINGS, CACHE KEYS, FIELDS ADDED SINCE, AND TWO AUDITS THAT NOW RECORD THE REPAIR OF THE DEFECT THEY FOUND (2026-09-14 04:13, lane I)

Runner `s26/i_verify_rerun.py` (commits `b0912487`, `2505e575`): each audit executed in-process
with writes under `verify/` and `bench_results/` redirected to `s26/results/verify/`, deletions
and move-outs there refused (six operations probed), the tracked file's sha256 asserted unchanged
after every run, the fresh JSON diffed leaf by leaf. One job per audit through `jobrun` (peak RSS
from `s26/jobs_done/verify_*.json`, logs `s26/logs/verify_*.log`, per-audit records
`s26/results/verify/<audit>.rerun.json`, fresh JSONs `verify__<name>.json`, report
`s26/results/verify/REPORT.md`). `determinism_audit` ran D1/D3/D4/D5 by function; its D2 was not
run because, as written, it moves the last two PRODUCTION cache records to a temp dir, re-runs
smoke8 (which does not contain them) and deletes the backup.

| # | audit | tag | tracked JSON | verdict | leaves identical / different / added / removed | wall | peak RSS | reading |
|---:|---|---|---|---|---|---:|---:|---|
| 1 | `vqe_lfo_audit` | CPU | vqe_lfo_audit.json | **IDENTICAL** | 32 / 0 / 0 / 0 | 0.3 s | 0.002 GB |  |
| 2 | `cvar_audit` | CPU | cvar_audit.json | **IDENTICAL** | 18 / 0 / 0 / 0 | 0.9 s | 0.003 GB |  |
| 3 | `ansatz_audit` | CPU | ansatz_audit.json | **IDENTICAL** | 19 / 0 / 0 / 0 | 3.3 s | 0.005 GB |  |
| 4 | `headline_audit` | CPU | headline_audit.json | **DIFFERS** | 123 / 0 / 8 / 0 | 0.6 s | 0.002 GB | all 123 recorded leaves identical; the 8 added leaves are `harness_cfg.*` fields `Config` gained since the tracked run (quantum, legacy, vqe_*, report_single_start_fit); no value changed |
| 5 | `legacy_audit` | CPU | legacy_audit.json | **IDENTICAL** | 30 / 0 / 0 / 0 | 27.8 s | 0.127 GB |  |
| 6 | `amber_audit` | AMBER | amber_audit.json | **IDENTICAL** | 31 / 0 / 0 / 0 | 42.2 s | 0.265 GB |  |
| 7 | `amber_platform` | AMBER | amber_platform.json | **NO FRESH OUTPUT** |  | 70.4 s | 0.304 GB | the script prints a table and writes no JSON (the tracked file is a hand-assembled summary of three probes); from the printed table the CPU arm reproduces the golden 1A13 interaction bit-exactly (-489.9138948277905, 726 evaluations) and the OpenCL arms differ run to run as the tracked verdict says they must (hybrid_double -489.874 vs the recorded -489.84..-489.95) |
| 8 | `leak_audit` | AMBER | leak_audit.json | **DIFFERS** | 24 / 1 / 0 / 0 | 118.0 s | 0.581 GB | ALL_CLEAN stands on 24/25 leaves; the one difference is `backends.project`: the tracked run predates `core.project` (`s8.project`), the fresh run records `core.project` |
| 9 | `grad_key_collision` | AMBER | grad_key_collision.json | **DIFFERS** | 7 / 3 / 0 / 10 | 2995.1 s | 0.578 GB | differs in the direction the fix implies: tracked `COLLISION: true`, no Config field; fresh `COLLISION: false`, `config_has_a_field_for_the_gradient: true`, and its two arms are bit-identical because `PROJECT_GRAD` in the environment no longer selects anything (both children run `exact`); the audit is stale the way `run_equiv2.sh` was (L19), 26 run-property leaves are timings |
| 10 | `equiv_compare` | CPU | equivalence.json | **DIFFERS** | 924 / 5 / 1 / 0 | 0.8 s | 0.003 GB | same verdict, BIT-IDENTICAL on all 8 targets including every AMBER quantity; the 5 differing leaves are the arm keys and directory names (tracked `29cc..`/`4077..` are gone; fresh compares the current baseline key `44a9..`, 8 records, against the production cache `1fc9..`, 126 records), hence `n_optimised` 8 -> 126 and the added `only_optimised` list |
| 11 | `projection_divergence` | CPU | projection_divergence.json | **ERROR** |  | 1.2 s | 0.002 GB | ERROR by construction: its two smoke8 arm caches (fd vs analytic) no longer exist, and against the baseline vs production keys it trips (`shape mismatch`) on 7 S11 baseline records whose `amber_ca` is null because AMBER declined at the 92% ceiling in that run; its question is answered by the fresh `project_arms` (126/126 bit-identical) |
| 12 | `project_arms` | CPU | project_arms.json | **DIFFERS** | 47 / 124 / 0 / 0 | 0.3 s | 0.008 GB | differs in the informative direction: the tracked S11 run (scan-builder era) had `ca` bit-identical on 0/8 with 22.3 A worst apart; the fresh run, baseline vs the shipped exact mode over the two on-disk caches, has `ca`, `fit_ca`, `avg_ca` bit-identical on 126/126 and `amber_ca` on 119/119 (the 7 null-AMBER baseline records excluded); the `fd` and `analytic` arms were passed as the production key (those caches are gone), so those two comparisons are identical by construction and carry no information |
| 13 | `determinism_audit` | AMBER | (none) | **no tracked JSON** |  | 59.8 s | 0.579 GB | no tracked JSON; fresh: bit-identical across processes and under 4 threads (D1), no key collisions (D3), `CacheCollision` raises (D4), config key complete (D5); D2 not run (it would delete two production-cache records) |
| 14 | `hazard_audit` | AMBER | (none) | **no tracked JSON** |  | 21.7 s | 0.525 GB | no tracked JSON; fresh: H2 thread-env energies identical 3/3, H3 alignment-invariant and pair distances bit-identical, H6 stable argsort, H7 longer-normalised and symmetric; H1 single-point pair differs by 1239 kcal/mol on 2.1e8 of builder strain (relative 5.8e-6, inside the 1e-5 bar of the integration test); H5 flags the alternate alphabet string in core/data.py, which is `ALPHABET_ALT`, present by design and never used to encode |
| 15 | `project_selfcheck` | CPU | (none) | **no tracked JSON** |  | 0.6 s | 0.003 GB | no tracked JSON (prints only): builder vs reference 1.732e-13 A, analytic vs central gradient 6.828e-06 on |g| = 1.540 |
| 16 | `project_exactness` | CPU | project_exactness.json | **DIFFERS** | 1894 / 1 / 2 / 0 | 587.0 s | 0.092 GB | 126/126 bit-identical, worst |d| 0.0, exactly as tracked; the 1 different leaf is `agg.s_per_target` (4.57 -> 4.65 s, a timing), the run-property leaf is `agg.seconds`, the 2 added leaves are `agg.complete` / `agg.requested` (bookkeeping added with the `_outpath` fix) |
| 17 | `project_equiv` | CPU | project_equiv.json | **DIFFERS** | 12431 / 508 / 2 / 0 | 1798.9 s | 0.275 GB | the four-arm table reproduces on every science leaf: agg.ref / ex / fd / an (synthesis 3.2148 / 3.2148 / 3.2145 / 3.2057, fit, projection_cost, fval) and avg_rmsd 3.048338 identical; all 508 differing leaves are the per-target and aggregate TIMINGS (t_ref / t_ex / t_fd / t_an, s_per_target; 504 + 4) and the 10 run-property leaves are seconds / speedups (exact-mode speedup 2.19x tracked -> 2.02x under seven concurrent jobs); 2 added bookkeeping leaves (complete, requested) |
| 18 | `project_degeneracy` | CPU | project_degeneracy.json | **DIFFERS** | 256 / 510 / 2 / 0 | 239.5 s | 0.09 GB | not like for like at the per-target level, same shape in aggregate: the tracked run measured the four starts under the module default of its day (analytic), the fresh run under the shipped exact mode, so no per-target best objective is identical (max |d| 0.179 in objective units); the degeneracy the audit exists to measure reproduces: 43 / 33 targets with a best-to-runner-up gap under 1e-2 / 1e-3 in both runs, 12 same-point ties in both, gap median 0.0667 -> 0.0643, runner-up distance median 0.846 -> 0.877 A, genuine branch ties (gap < 1e-3, structures > 0.5 A apart) 1 -> 3 |
| 19 | `project_iters` | CPU | project_iters.json | **DIFFERS** | 4 / 14 / 7 / 0 | 384.7 s | 0.273 GB | not like for like: the tracked file is a 10-target run with three arms (ref / fd / an, 90 solves each, 10.0-13.3% of solves hitting maxiter=300, nit mean 157-161) from before the exact mode existed; the fresh CLI default is 24 targets with four arms (ref / ex / fd / an, 216 solves each: 18.5 / 18.5 / 18.1 / 17.6% hitting maxiter, nit mean 171-181), written as project_iters_n24.json; the audit's question (is maxiter binding?) reads the same way at both sizes: a minority of starts hit the cap, and the ref and ex arms are identical solve for solve |
| 20 | `project_stability` | CPU | project_stability_partial77.json | **DIFFERS** | 469 / 2 / 0 / 1 | 1308.6 s | 0.273 GB | IDENTICAL on every science leaf (469/469): the reference against itself under a rigid motion on the same 77 targets reproduces the tracked partial exactly (d(arm) median 0.0323 A, worst 1.6246 A, |dfval| max 0.0890, 62 / 46 / 18 targets moved by more than 1e-3 / 0.01 / 0.1 A); the 2 different leaves and 1 removed leaf are bookkeeping: the tracked file was a run of 126 stopped by the RAM gate at 77 (`stopped_by`, `complete: false`, `requested: 126`), the fresh one was asked for 77 and completed |

20 audits: 5 IDENTICAL, 11 DIFFERS (each explained), 3 with no tracked JSON to diff, 1 ERROR by construction, 0 NOT RUN. Not re-runnable: `recon_containment_audit.json` (written by a deleted script), `project_inputs.json` (the harvested input the projection audits read).

Reading of the whole: every audit whose tracked JSON records a comparable experiment reproduces
its science (`vqe_lfo`, `cvar`, `ansatz`, `legacy`, `amber_audit` bit-for-bit; `headline` on all
123 recorded leaves; `leak` ALL_CLEAN; `equiv_compare` BIT-IDENTICAL on 8/8; `project_exactness`
126/126; `project_equiv` on every science leaf of the four-arm table; `project_degeneracy` in its
aggregate shape). Two audits now record the repair of what they found: `grad_key_collision`
(COLLISION true -> false, the mode is in the cache key) and `project_arms` (0/8 -> 126/126
bit-identical once the shipped exact mode replaced the scan builder). Two are not like for like
at the per-target level because the module default changed under them (`project_iters` size and
arms, `project_degeneracy` mode) and say so. `projection_divergence` cannot run on the current
caches (its arms are gone) and its question is answered by `project_arms`. Three audits print
and write no tracked JSON; their fresh outputs are recorded. `amber_platform`'s tracked JSON is a
hand-assembled summary with no writer.

---
## L142 -- SPRINT 26 CLOSE (2026-09-14 04:14, coordinator)

The sprint closes at the user's extended time. What stands, with where it lives:

1. Phase 0: `s26/EXAMINATION.md` (module map, two dataflow traces, the 35-claim ledger, the
   five defects, the pinned hashes, the not-in-git census), audited in
   `s26/EXAMINATION_AUDIT.md` (L31, one documentary item fixed, L32/L34), signed off (L33). The
   production numbers reproduce to 0.0 on all four bases (3.0483 / 3.2041 / 3.2148 / 3.2355,
   n = 126; T030 = 1S9Z 0.181981), re-run on HEAD by the Adversary.
2. The three proposals: A REPLACE (`s26/PROPOSAL_A.md`; L68, L70, L75, L139, L140: not
   measured on either seed; the mechanism is absent: product-state optimum, full so(2^n) DLA
   from depth 2, inert growth); B REPLACE (`s26/PROPOSAL_B.md`, replacement
   `s26/PROPOSAL_B_REPLACEMENT.md`, the trainability paper; L13, L106, L107, L120); C KEEP WITH
   EDITS (`s26/PROPOSAL_C.md`; nine rungs, none beats the shipped prior, the ESM-2 650M channel
   worth -0.208 A built chain and -0.330 A selection at 5/5 folds; AMBER a validity step only;
   the routers closed again; C5 null to harmful; L56 to L67, L72, L93, L99, L103, L105, L115,
   L124, L132). The slide 11 ruling is L117 as amended by L122.
3. The tournament: `s26/TOURNAMENT.md`; 18 IDEA files, 29 PREREG files, every
   ranked survivor run or recorded as not run with its reason (L51, L137, L133): the cis gap
   (0 cis on the instrument, floor 0.083 A, L22/L38/L89), the steric reject (harmful on both
   bases, L43/L86), the 2/60 leak bound (0.028 A mean, 0.151 worst, MINOR, L44/L58), window
   ensembling (null with power, L84), branch selection (L88), strain difficulty (a proxy for
   the pool's own spread, L53/L121/L123), rotamer relief (L131), the identity and tie-break
   floors (L52, L64), provenance (L85/L109), the AMBER prior partner (L105), the recall gradient
   and the memorisation ladder (L90, L91), A3 (L125), the per-step DLA (L138).
4. Repository: operational items 1 to 5 and defects 6a to 6d closed on branch `s26` (L1, L6 to
   L10, L15, L16, L18 to L21); the suite at 370 tests, 368 passed, 2 absent-artefact skips with
   the opt-in tier (L113); the frozen results-lab rebuild reproducing every number and gate
   (L127); 20 of 20 `verify/` audits re-run with the science reproducing (L141); the AST gate
   50/55 identical and the 5 exactly the ledgered edits (L98); the governor v2.4 and jobrun
   v2.3 with their incident history (L18b, L36, L40, L41, L59, L74, L78, L83, L92, L111,
   L128, L129).
5. The deck `vqe_research_overview.pptx` (built from artefacts, 11 slides, L61, L126, L134 and
   the L139/L140 rebuild) with `s26/PRESENTATION_CHANGES.md`; `s26/TRAINABILITY_PAPER_OUTLINE.md`;
   `s26/REPORT.md` (lane E's final pass ends at about 04:45 with this ledger as Appendix C);
   `s26/RETRACTIONS.md` (R1 to R11); `s26/DELIVERABLES_CHECK.md` (13 of 14 rows at 04:12, the
   report row at 04:40); the `docs/FINDINGS.md` corrections block (L130 and its two added rows);
   eight `s26/agent*_FINDINGS.md`; `s26/STATUS.md` hourly.
6. Not run, recorded with reasons: the raw rung's last fold and evaluation (L133), the `attn`
   input (3 to 3.5 GB, never fit), coherence-penalised training, lane W's remaining Part B
   models and F3 (L137), `project_stability` beyond what L141 states, the C5 alpha grid beyond
   0.5. No superseded claim was deleted anywhere; every retraction is appended.
7. Governed jobs: 184 records in `s26/jobs_done/`, every one with an exit code and a peak
   RSS; the box never exceeded the 93% ceiling under the governor's control (the two
   excursions above it, L40 and L74, were the user's own load with the governor dead or every
   job suspended). One job of lane Q (`a2_dla_a1`, a second run) may still be finishing when
   this entry is written; its record lands in `jobs_done` and changes nothing above.

The findings, not the folding, remain the output: no proposal survived as stated, every
mandatory direction has a measured answer, and the presenter walks in knowing which questions
are closed and by which sprint.

---


## L143 -- LANE PR: L139 / L140 / L138 APPLIED TO SLIDES 4, 5 AND 8; A1 READS "NOT MEASURED ON EITHER SEED" WITH BOTH SEEDS' MDE MULTIPLES AND THE COMPARATOR'S 0.033 A SEED VARIANCE; THE DECK IS FINAL; 550 REGISTERED NUMBERS; VERIFICATION PASS (2026-09-14 04:16, PR)

Rebuilt (`python s26/pr_build_deck.py --no-figures; python s26/pr_changes.py`). Slide 8's A1 block now quotes both seeds
from `s26/results/a1_stats.json` and `a1s1_stats.json` (seed 0: -0.0138 / -0.0222 A at 0.23x / 0.36x MDE, fold CIs
spanning zero; seed 1: -0.0449 / -0.0523 A at 0.71x / 0.79x MDE, fold CIs [-0.096, -0.004] / [-0.109, -0.005], 4/5 folds)
and says NOT MEASURED on either seed (underpowered at seed 0, Type-M at seed 1), with the L140 qualifier: the grown
circuits are seed-stable (3.214 / 3.216 A) and the deployed fixed circuit moves 0.033 A between seeds (3.228 / 3.261,
`means/fixed_zrank_it50/rmsd_q_synth` of the two stats files), so the size of the contrast is set by the deployed
circuit's seed, not by growth; the spoken text carries PROPOSAL_A.md addendum 3's corrected sentence ("on two seeds ...
0.014 to 0.052 angstroms closer ... between a quarter and four fifths of what the comparison can resolve ... the spread
between the seeds comes from the deployed circuit"); the verdict line says the endpoint is not measured on either seed
and REPLACE stands on the product-state diagnosis. L138 added as one bullet (`s26/results/q_dla_a1.json`: dim 7 exactly
where nothing was appended, up to 530 where inert strings were, median 1025 of 8128 for the alpha = 0.25 L2 sets, no
grown set reaches so(128)). Slide 5 carries the same seed qualifier beside the deployed circuit (3.228 / 3.261 A on two
seeds; an ADAPT-grown circuit seed-stable at 3.214 / 3.216). Slide 4 states the 3.2148 A production number as seed-free
(Config quantum = False; 3.2148 on both A1 seeds). The earlier builds' "no change" / "null" wording for A1 is superseded
(R11) and recorded in `s26/agentPR_FINDINGS.md` section 7 and `s26/PRESENTATION_CHANGES.md`. Verification
(`s26/pr_verify.txt`): 11 slides; 0 U+2014, 0 U+2013; 0 banned words; spoken words 248 / 248 / 249 on slides 8 / 9 / 10;
title, body and notes on every slide. Registry 550 tokens. Nothing on the deck is live, DRAFT or PENDING.

---

## L144 -- DELIVERABLES CHECK, THE REPORT ROW (04:16, ON THE COMMITTED TEXT ed5e44cd): OK; TWO ITEMS PENDING FOR LANE E'S FINAL PASS (A1 ON TWO SEEDS; THE REPORT CHECK RE-RUN ON THE FINAL SHA); ONE QUOTED BANNED WORD (2026-09-14, A)

`s26/REPORT.md` at `ed5e44cd` (3,644 lines; Parts I-IX, Appendices A-D). Rule 13: one banned
word, line 3584 (Appendix D quotes L77's phrase for the seed re-run; rephrase), no U+2014, no
U+2013. Every qualifier in `s26/agentA_FINDINGS.md`'s table is carried (line numbers in
`s26/DELIVERABLES_CHECK.md`), including L122's "exact and complete" applied after the verbatim
L117 line (2874) and the S13 Pauli mean weights cited to the dossier "as cited" and kept off the
slides (Appendix B row 97; 2907). Pending for lane E's 04:15-04:45 pass: (i) A1 on two seeds
(L139, L140, R11) is not yet in Part VII or Appendix B, since it landed after the 04:04 commit;
(ii) `s26/results/e_report_check.json` (175 of 175) was computed on report sha 45eaa279, earlier
than HEAD's aa528df3, so the check is re-run on the final commit as lane E's last act (L77). With
those two, row 12 is OK and the Part 10 table is 14 of 14. Final state of this lane's files:
`EXAMINATION_AUDIT.md`, `TOURNAMENT.md`, `RETRACTIONS.md` (R1-R11 and the prior-sprint
table), `DELIVERABLES_CHECK.md` (final), `agentA_FINDINGS.md` (with the per-slide qualifier
table), 24 ledger entries (L31, L34, L45-L49, L54, L55, L70, L71, L79-L82, L94-L97, L120, L121,
L135, L136, L140, this), four artefacts (`a_reproduce_head.json`, `a_c26_phi_mae.json`,
`a_dla_check.json`, `a_strain_vs_spread.json`, `a_ladder_isolations.json`).

---

---

End of Appendix C. Entries after L144, if any, are in `s26/LEDGER.md` itself.
