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
stands; Parts VII and VIII fill from the S26 ledger as verdicts land (DRAFT until the sprint
closes) and Appendix C is reproduced at the close; Appendix B lists every number of the drafted
parts, and `python s26/e_report_check.py` (or `python s26/examine.py --report-check`) re-reads
each one from its artefact; Appendix D reconciles this report with `docs/REPORT_S26.md`.

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
  (`s24/LEDGER.md` L4; `s26/results/i_identity_audit.json`, S26 ledger L15). At the looser
  >= 0.6 identity criterion 13 of 60 benchmark targets have such a pool window
  (`docs/FINDINGS.md:4571`, S9-10), priced at +0.0030 A in S10-4. The dev price, re-derived in S26
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
| S5, the first measurement (`docs/FINDINGS.md:368-403`) | -0.023 | not recorded |
| S9, 10 qubits, 1024 amplitudes enumerated (`core/quantum.py:42`) | +0.655634 | 0.758 |
| the consolidation audit, exact / sampled (`verify/cvar_audit.json`) | +0.6847 / +0.685 | not recorded |
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

### V.12 S26 additions as they land

**A2 checked (L45).** The Adversary's independent re-derivation `s26/results/a_dla_check.json`
(no shared closure code, incremental Gram-Schmidt rank) agrees with `q_dla.json` at every cell
at n = 4 to 7, L = 1 to 4, for all four conjugation and gate-order conventions; the 8128 / 1025
reconciliation (the fixed ansatz's algebra against the strings ADAPT selected, 1025 / 8128 =
0.126) is correct. The "S25 slopes are a statement about depth 3" reading is an inference from
the 2-design literature, not a measurement here.

**A4, grown against fixed circuits (L35, checked L47).** Part VII.1 carries the table. In one
sentence for this part: on the deployed Hamiltonian an operator-growing (ADAPT) construction
produces product circuits at alpha = 1, whose gradient variance does not decay with width because
there is nothing entangled to train, and at alpha = 0.25 a 2-local family that decays like the
fixed ansatz (-0.302 against -0.243 log2 per qubit); no width-scaling argument for the selector
comes out of it (`s26/results/q_var.json :: results/slopes`).

**A1, the ADAPT endpoint (L68, L70, L75).** Part VII.1 carries the numbers. For this part: on the
deployed Hamiltonian the optimum at alpha = 1 is a product state on every real target
(KL(Gibbs || product of marginals) at most 7.9e-4 nats), a 7-parameter RY layer reaches it, the
21-parameter fixed circuit stops 0.90 nats short, an adaptive ansatz free to entangle appends
operators worth less than 1e-3 nats and leaves a product state, and reaching the optimum exactly
moves the built chain by -0.014 to -0.022 A, a third of what the comparison resolves (MDE 0.059 to
0.061). Proposal A's verdict is REPLACE (L69). **A3.** Pending.

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
filter plus consensus medoid at -0.172 A [-0.316, -0.027] single window, 74W/47L (S8-8,
`docs/FINDINGS.md:2893-2935`); the channels' errors are
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

Filled from `s26/LEDGER.md` from L33 (PHASE 0 SIGNED OFF, 09:05) onward; each item names its
ledger entry, its artefact, its pre-registration and the Adversary's check where one has landed.
Nothing is written here before its ledger entry exists. Bases are named on both sides of every
contrast. "Pending" means the run or its ledger entry has not landed.

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
for two hours and relaunched (L92). The box was shared with the user's own load all night, so
memory, not CPU, set the parallelism (headroom 1.4 GB after the pause, L41).

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
conventions; n = 4, 5 numeric and symbolic counts agree 6 of 6.

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
of order -0.25 to -0.30 and far from a 2-design's -1.0, is safe). Verdict: A4 gives Proposal A no
width-scaling argument.

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
inert, and nothing in the null depends on the count. Verdict for Proposal A (`s26/PROPOSAL_A.md`,
accepted L69): REPLACE, not the expected keep-with-edits, because the mechanism is absent rather
than weak: the Hamiltonian is a constant ladder whose optimum is a product state, the fixed
ansatz's algebra is already the full so(2^n), the grown circuits give no width-scaling argument,
and the endpoint is null at a stated resolution of 0.06 A.

**A3, a target-dependent Hamiltonian.** Pending (`a3_build`, suspended by v2.1 and never resumed
by v2.2, adopted and resumed by v2.3, L83).

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
stage-1 replication PH already ran (L87). The lane's own re-projection of the same rung
(`p_deliver_shipped`, 567 s) is saved beside it as `p_best_rung_chains_rebuild_basis.json`
(mean 3.2126, the L57 basis) as the cross-check that the ladder's anchor and the production
chains are the same object up to the projection's multi-start sensitivity; it overwrote the
delivery file for six minutes and was moved (L116). Pending: `raw` (its training stopped by hand
at the 21:50 stall, L74), the deferred `coherence_penalised_training` (1.25 GB)
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

**C5, the common-mode prediction.** Did not complete before the close; it keeps its
pre-registration (L117).

### VII.4 The tournament entries

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
  carrier-out change against two control-out chains per fold, is measured on one control so
  far, 9BAF out of fold 0, five queued). The sign follows the cross-deposit distance of the
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
carried no pool-disagreement control (the top-75's own pairwise CA-RMSD spread), so whether
`moved` is that spread in disguise waits on the Adversary's `a_strain_vs_spread` (relaunched,
L92); the calibration-flag wording may be quoted now, the novelty sentence waits.

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

**Rotamer relief, coherence-penalised training, the provenance readout test H_P3, the
pool-spread control for strain difficulty.** Pending.

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
identity floor (L80), strain (L81, provisional), ensembling and provenance (L94), the reject on
the built chain and the C3 replication (L95), branch selection (L96) and the tight floor, the
recall gradient, the memorisation ladder and the mix rung (L97) likewise. Lane I's extension
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
2016 RMSD values exact, 2142 of 2142 PDBs byte-identical). Lane I was then starved for 75
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
| A: qubit-ADAPT-VQE in place of the fixed ansatz | REPLACE (the prompt expected keep with edits) | `s26/PROPOSAL_A.md` | the deployed objective's optimum is a product state on every real target (KL to the product of marginals at most 7.9e-4 nats); an adaptive ansatz appends only inert operators; the fixed circuit's Lie algebra is already the full so(2^n) from depth 2; the grown circuits give no width-scaling argument; the endpoint is null at 0.23x and 0.36x its MDE with a 0.06 A resolution (L27, L35, L45, L47, L68, L69, L70, L75, L82) |
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
and the circuit is closed as an accuracy lever by three independent facts.

### VII.7 S26 retractions (`s26/RETRACTIONS.md`)

Kept by the Adversary, one entry per retraction, the claim verbatim with the artefact that
contradicts it; nothing superseded is deleted anywhere.

| entry | claim | now stands | where |
|---|---|---|---|
| R1 | the governed test count 369 / 356 (`s26/EXAMINATION.md` D, C35; L28 item 2) | 370 tests, 357 passed, 0 failed, 13 skipped, 0 memory-guard skips (`s26/results/test_run.json`) | L31, L32, L34 |
| R2 | S25's "no barren plateau at any width measured" (the slopes of `s25/results/q_plateau.json`) | the slopes stand; the claim is scoped to depth 3, because the algebra is the full so(2^n) (A2); for `docs/FINDINGS.md` at the close | L27, L45; Part V.9 |
| R3 | lane Q's H2b, dim(DLA) at (n = 7, L = 3) below 8128, guess 4095 (`s26/PREREG_A2.md`) | 8128 = so(128) from depth 2; falsified by lane Q itself | L27; Part V.9 |
| R4 | lane PH's registered expectation of cis bonds in other ensemble models (`s26/PREREG_cis.md`) | 0 of 1,966 ensemble models, 0 of 126 natives, 0 of 2,352,893 windows | L22, L48; Part VII.4 |
| R5 (L75, L82) | L68's "L-BFGS growth at alpha = 1 stops with no operator selected on 78 of 78 targets" (and the same sentence in `s26/PROPOSAL_A.md` sections 3 and 5); the Adversary's L70 caveat 2 "L-BFGS grows nothing" | operators are appended on 60 to 78 of 78 targets and are inert (at most 1.2e-4 nats, angles below 0.02 rad, product state to 4.1e-4 nats); the verdict does not move | L73, L75, L76, L82; Part VII.1 |

Not retractions, recorded there as sourcing corrections: C27's +0.0004 / +0.0030 (re-derived
exactly on the lam = 0 chain by lane W, L44; the S10 artefacts are in git history at `5fa05cd`).
No prior-sprint accuracy claim has been contradicted by an S26 endpoint result to date; L39, L43
and L44 confirm the record.

## PART VIII. CLOSED AND OPEN

Two tables. The rule: an open item is never moved to closed on one failed experiment; it moves
when it has been measured to its own ceiling or refuted by a pre-registered falsifier, and the
row names which. The standing pre-S26 list is `docs/STATE_BRIEF_2026-09-12.md` 5.6 and
`docs/CONDENSED_REPORT.md`, "Closed - do not re-fund"; Part VI carries every closure with its
sprint. This part lists the S26 movements and the items that remain open at the time of writing
(DRAFT until the sprint closes).

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
| "Refine with physics" (the production AMBER relaxation) as an accuracy step | C3 stage 1: worse than do-nothing, worse than a matched random move, worse than a matched move toward a pool member; replicates S16 L27; kept as a validity step on 124 of 126 | S26 L39, L46; Part VII.3 |
| The steric reject at a physical threshold, with or without refill, on both bases | falsifier fired the other way on the point cloud and replicated on the built chain (harmful at 1e3 and 1e4, underpowered null at 1e5 and 1e6); the sixth and seventh instruments for the filter form of the functional lever | S26 L43, L54, L86; Part VII.4 |
| The cis-peptide gap and the omega non-planarity gap as levers on this instrument | 0 cis natives, ensemble models or windows; the two-bond projection worth 0.000 A; the tight representation floor 0.083 A (non-planarity about 0.04 A) | S26 L22, L38, L48, L89; Part VII.4 |
| Branch selection by the relaxed energy as an accuracy step | +0.0055 against the production choice (0.18x MDE; a 0.03 A gain would have been seen); the energy discriminates among branches as well as the objective and no better | S26 L88; Part VII.4 |
| Test-time ensembling of retrieval keys (fixed K) | -0.0005 against shipped (0.02x MDE), +0.0042 against a zero-information resample; a 0.028 A gain excluded | S26 L84; Part VII.4 |
| A provenance-weighted readout (down-weighting fragment windows in the average) | +0.0066 against uniform (0.15x MDE; a 0.044 A gain excluded), -0.018 against the permuted-weight control; the ORACLE class gap (0.42 A in the pool) cannot be converted into an emission gain at this n | S26 L85, L94, L109; Part VII.4 |
| Routers for a per-target shortlist size or scale (C4: twelve m routers and six s routers on six native-free feature blocks) | none clears its MDE, all but one point the harmful way, the s routers predict the wrong sign; the sixth to seventeenth constructions land where the first five did | S26 L110, L115; Part VII.3 |
| ESMFold (B1) on this box | infeasible on three independent grounds | S26 L13; Part VII.2 |
| Proposal B: a native-free characterisation of the set where the pipeline beats sequence-only (B3) | the sign classifier is at its permutation null against both comparators; only the size of the gain over a helix is partly predictable, by the pool's strand content; verdict REPLACE | S26 L14, L106, L107; Part VII.2 |
| AMBER as a distribution inside the prior (the last AMBER form in the record) | the leave-fold-out choice is lam = 0 on every fold; every mixture cell worse, beta without a consistent sign | S26 L105; Part VII.4 |
| Proposal A: an adaptive (ADAPT) ansatz in place of the fixed one, as an accuracy or trainability lever | A1 null at the registered threshold (0.23x and 0.36x MDE, fold CIs spanning zero, resolution 0.06 A); the optimum is a product state on every real target; the fixed algebra is already maximal (A2); no width-scaling argument (A4); verdict REPLACE | S26 L68, L69, L70, L75; Part VII.1 |
| More prior capacity (`wide`), a per-fold PCA (`pca32f`), 128 ESM components (`pca128`), a pool-histogram mixture (`mix`), a triangle-update PairNet on the same inputs (`pairnet`) as prior levers; a smaller language model (`esm8m`) or none (`noesm`) as alternatives | null-to-worse on the built chain (`wide` +0.032, `pca32f` +0.018, `pca128` +0.076, `pairnet` +0.041, all under 0.5x MDE; `mix` chooses lam = 0 on 5/5 folds); removing or shrinking the ESM channel is worse (`noesm` +0.208, `esm8m` +0.242, both 5/5 folds, Type-M zone) | S26 L62, L63, L66, L67, L72, L93, L99, L103; Part VII.3 |

### VIII.2 Open

| item | what would close it | where |
|---|---|---|
| Whether a better distance predictor is obtainable (the only steep lever, -2.15 A per unit toward truth); every input this machine can compute is measured flat (nine rungs) | a larger language model than this machine can hold, on a bigger machine; Proposal C's kept form (L117) | `s24/LEDGER.md` L13; `s26/PROPOSAL_C.md`; Part VII.3, VII.6 |
| Proposal A's target-dependent Hamiltonian (A3) | the pre-registered A3 verdict (`a3_build` running) | Part VII.1 |
| The rest of the C2 ladder (`raw` and the deferred rungs; the best rung so far is the shipped prior, L112) | a rung beating the shipped prior on the fold-clustered CI of the built chain, replicated | Part VII.3 |
| C3 stage 2 (the relaxation on the best C2 rung), which is the stage-1 replication already run because the best rung is the shipped prior | nothing further unless `raw` or a deferred rung beats the shipped prior | S26 L39, L87, L112; Part VII.3 |
| The 2/60 benchmark self-copy leak, now bounded MINOR (dev proxy 0.008 A with both channels measured on all four dev self-copies; own-native envelope 0.028 A mean CI to 0.151 A worst target on the built chain; 0.194 at the worst target on the selection basis); F3's control clause open (one of six control-out models) | by design only a fresh benchmark, which does not exist; the control-out models and a second seed of the 2P5H retrain if time allows | S26 L44, L49, L50, L55, L58, L108; Part VII.4 |
| Where the target-specific third of the pool's coherent error comes from, and whether any native-free proxy is strong enough to act on it | a native-free proxy reaching the in-band ordering 2 A needs | `s19/LEDGER.md` L11, L14; `s17/LEDGER.md` L23 |
| Publishing the trainability half | a manuscript from Part V.10 with V.9's scope correction | `s13/`, `s25/QUANTUM.md`; S26 L27 |
| The tournament entries not run or not complete at the close: the common-mode prediction (C5, keeps its pre-registration), rotamer relief, coherence-penalised training | their pre-registrations' falsifiers | `s26/TOURNAMENT.md`; `s26/PREREG_*.md`; Part VII.4 |
| The tie-break noise floor: measured (0.004 A on the 126-mean, 0.024 A paired MDE between conventions, built chain); not a lever | nothing; it is the floor every cross-run hundredths-level claim is read against | S26 L64, L71; Part VII.4 |
| Strain as difficulty: measured as a calibration flag (partial rho +0.433 of `moved` with the built-chain error), forbidden as a lever by its pre-registration; provisional | the pool-spread control `a_strain_vs_spread` (whether `moved` is the top-75's own disagreement in disguise) | S26 L53, L81; Part VII.4 |
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
on the point cloud), `results/summary/target_map.json` (T001 to T126), `results/site/`.

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
    python s26/e_module_map.py        # s26/results/module_map.json (725 modules at the last check; grows with the sprint's scripts)
    python s26/e_hashes.py            # s26/results/pinned_hashes.json (--check to compare)
    python s26/e_claims.py            # s26/results/claim_search.{json,txt}; benchmark files excluded

`examine.sh` and `examine.bat` at the root call `examine.py`; there is no Makefile. Exit status
is non-zero if any step reports a problem.

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
artefact; "as asserted" means a passing test pins it.

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
| +0.0103 [-0.1596, +0.1803], 31W/29L | I | `docs/FINDINGS.md:4479-4607` (S9-10), naming `s9/final_report.json` | as cited |
| 2.9614, 3.8122, 22.3%, 0.649, 80 of 126, 3.803955, 1.8e-15, 3.867 | II, III | `s26/results/e_reproduce.json :: summary/virtual_bond, n_rows` | as stored |
| +0.1664, +0.0977, [+0.0049, +0.3379], 16 of 126 | II, III | `s26/results/e_reproduce.json :: summary/projection_gap_arm_minus_avg, n_rows` | as stored |
| 0.0342, 0.0958, 0.2051 | II | `s26/results/q_mde_reference.json` | as stored |
| 0.394%, 1.18% | III | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` (hamiltonian stage); `s25/results/q_gibbs.json :: results/spectrum_target_independence` | as stored |
| 0.6758 (68%) | III | `s23/results/errdecomp.json :: rows[*]/f_common` mean | 0.675770 |
| 0.220, -560, 125/126, 5.38, 4.86 | III | `s26/results/ph_c3_nativefree.json` (S26 ledger L24) | as stored |
| 5370, -1290.6, 0.19, 52.0; 6.0e12, 1262.4, 0.61, 1166.1; 1172.7 (2BP4) | III | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` (relax stage); `bench_results/cache/1fc9f2dcf489e2fb/2BP4.json` | as stored |
| 7193, 9814, -1..19, -6..32, 144, 420, 105, 91, 183, 17 centres, 760, 488, 483, 0.6926..7.2688, 2.2420..5.1791, 401 98 193 161 284 260 372 373 381 382, 44/449, 3.75, 1.97, 1.16, 8, 128, 5.05..21.07, 0.23..2.70, 5.62..23.44, 7.58, 5.34, 1.40 | III | `s26/results/e_trace_1S9Z.json`, `s26/results/e_trace_9KAR.json` | as stored |
| 3.3135, 3.3443 | III | `s8/integrate_vqe.json :: arms/vqe_LFO/sel, arms/medoid128/sel` (reproduced by `tests/test_pipeline.py::test_the_four_components_all_execute_and_reproduce_their_published_numbers`) | as stored |
| 0.288 | III | `docs/FINDINGS.md:1970-2037` (S7 finding 11; its per-target artefact, s7/repr_tune.json, is lost, S26 ledger L11) | as cited |
| -0.172 [-0.316, -0.027], 74W/47L | VI | `docs/FINDINGS.md:2901` (S8-8 table) | as cited |
| 8.3 to 23.3 | IV | `s15/LEDGER.md` row 0.8 (`s15/results/audit_amber_cost_sweep.json`) | as cited |
| 0.016 | III | `docs/FINDINGS.md:2379` (S8-6 heading) | as cited |
| +0.142 | III | `ARCHITECTURE.md` section 2.5 | as cited |
| +0.0207 | III, IV | `bench_results/baseline_tuning126.json :: science/rmsd_full/mean, science/rmsd_arm/mean` | derived: 3.2354598538973844 - 3.214765154210998 = 0.0207 |
| 1.386 | II | `docs/FINDINGS.md:2751-2760` (S8-8) | as cited |
| 470, 13, 47 of 500 | II, III | `tests/test_data.py:305`; `README.md` (fold repin; layout) | as asserted |
| 23, 5 | II | `docs/FINDINGS.md:60-118`; `s25/LEDGER.md` L5, L7, L9, L11, L15 | as cited |
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
| +0.655634 / 0.758, +0.566586 / 0.519, +1.000000, +0.994 | V | `core/quantum.py:42` (S9 instrument); `s25/results/q_verify.json :: results` (S25 re-verification); `docs/FINDINGS.md:3160` (the sampled estimator, S8-9) | as stored |
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
| 16.75 GB, 8 cores, 6.43, 67%, 93%, 92%, 88 / 90 / 93 / 95%, 60 s, 15 s, 25 s, 5 s | IX | `README.md` ("The S26 resource governor"); `docs/STATE_BRIEF_2026-09-12.md` section 8; `s26/governor.py` | as in source |
| 5.5 GB, 555 MB, 1.5 GB, 61, 63k | IX | `README.md` ("Data this repository does not carry"); `docs/STATE_BRIEF_2026-09-12.md` section 8 | as cited |
| 0.18198112330908295 | IX | `s26/results/e_reproduce.json` (row 1S9Z, `rmsd_arm`); `s26/results/e_trace_1S9Z.json` | as stored |
| 1.692, 0.872, 0.324, 11, 3, 8, 2, 601a39c7 | IX | `s26/TEST_RUN.md`; `s26/results/test_run.json` | as stored |
| 725 | IX | `s26/results/module_map.json :: n_modules` | 725 |
| 2.9661, 0.1824, 8.2342, 28.6%, 50.8%, 0.7150, 0.1643 | II, D | `results/summary/leaderboard.json :: rows[0]/median, rows[0]/best, rows[0]/worst, rows[0]/frac_under_2, rows[0]/frac_under_3, rows[0]/corr_with_pool_best, rows[0]/basis_delta_mean` (pointed to by `docs/REPORT_S26.md` B.1) | as stored |
| +0.0061, 0.19x, 61W/65L | IV, D | `results/summary/leaderboard.json :: rows[1]/paired_effect, rows[1]/effect_over_mde, rows[1]/wins, rows[1]/losses` | as stored |
| +0.7070 (2.14x), +0.8322 (2.67x), +0.0097 (0.43x), +0.1668 (0.91x), +0.0833 (0.70x), +0.6259 (1.89x), +0.2047 (1.07x) | IV, D | `s25/results/phys_suite.json :: vs_incumbent` (point cloud against point cloud) | as stored |
| 1.9962, 0.6159, 0.2824, 0.2412 | III, D | `s25/results/calib.json :: rows[*]/z_sd, rows[*]/cov90, rows[*]/cov50, rows[*]/multimodal_frac` (means over 126 rows) | as stored |
| 2.8334, -0.2150, 2.4x, 113W/13L | III, D | `s24/LEDGER.md` L13 | as cited |
| 1.75, 51, +0.054 | III, D | `s25/LEDGER.md` L6, L12 | as cited |
| 5.551e-17 | V, D | `s25/results/q_verify.json :: results` | as stored |
| 0.1 nats, fired | V, D | `s25/results/q_gibbs.json :: results/falsifier` | as stored |
| -0.023 | V, D | `docs/FINDINGS.md:368-403` (S5 section 6) | as cited |
| 0.6847, 0.685 | V, D | `verify/cvar_audit.json :: D_tail_baseline_cosine_vs_paramshift, D_sampled_tail_cosine` | as stored |
| 8.371, 2537.568, 303.137 | IX, D | `bench_results/compare_tuning126.json :: end_to_end_speedup, wall_s` | as stored |
| 3.989, 3.213 | VI, D | `s12/SPRINT12_DOSSIER.md` section I | as cited |
| 13 of 60, +0.0030 | II, D | `docs/FINDINGS.md:4571` (S9-10); S10-4 (artefact in git history at `5fa05cd`, S26 ledger L31) | as cited |
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
<!-- APPENDIX B ROWS -->

## APPENDIX C. THE S26 LEDGER (DRAFT: reproduced at the close)

Reproduced verbatim from `s26/LEDGER.md` at the close of the sprint; until then this slot
points at the live file, whose tail at the time of the last report commit is noted here.

<!-- APPENDIX C LEDGER (filled at close) -->

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
| the S5 measurement of the gradient defect, cosine -0.023, and the consolidation audit's exact / sampled tail cosines 0.6847 / 0.685 | `docs/FINDINGS.md:368-403`; `verify/cvar_audit.json` | V.5 |
| the entangler-deletion contrast's source, -0.013 [-0.095, +0.077] | `s20/LEDGER.md` L1 | V.7, Appendix B |
| the optimised arm's end-to-end speed-up, 8.371x (2537.568 s to 303.137 s, both cold, 1 against 8 workers) | `bench_results/compare_tuning126.json :: wall_s, end_to_end_speedup` | IX.4 |
| B1's feasibility numbers | `s26/results/b1_feasibility.json`; `s26/LEDGER.md` L13 | VII.2 |

### D.2 Disagreements, and the artefact that decides

| the two readings | decision |
|---|---|
| 3.2126 (their headline's leaderboard rebuild) against 3.2148 (this report's production cache); 0.1643 against 0.1664 for the projection gap; T030 0.18242 against 0.18198 | Not a disagreement: the same built chain on two emission paths (the results lab re-projects the stored cloud and round-trips through PDB quantisation, up to 0.171 A per target, `s26/LEDGER.md` L14; `s26/EXAMINATION.md` H). Both reports say so. This report quotes the cache number as the production result and names the rebuild beside it. |
| `bench_results/compare_tuning126.json :: science.rmsd_arm` (theirs) against `bench_results/baseline_tuning126.json :: science/rmsd_arm/mean` (this report) | The same values: `compare_tuning126.json :: science_delta` is 0.0 on every science key except two at 4.4e-16. Either path is valid. |
| 5.6e-17 (this report, from `s25/QUANTUM.md`) against 5.551e-17 (theirs) | The artefact says 5.551e-17; corrected here. |
| the benchmark leak count: "13 of 60 targets have a pool window at >= 0.6 identity" (their B.2, from `docs/FINDINGS.md:4571`) against the 2/60 of this report | Two criteria, both true: 13/60 is the >= 0.6 identity count of S9-10 (priced +0.0030 A in S10-4, artefact in git history at `5fa05cd`, L31); 2/60 is the verbatim self-copy defect of S24 L4 that lane W bounded (L44). This report keeps 2/60 as the declared leak and now names the 13/60 criterion beside it (II.1). The two benchmark target ids printed in `docs/FINDINGS.md` and in their table are not repeated here. |
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
- Their Appendix C (the corrections ledger) restates `docs/FINDINGS.md:60-118` and the sprint
  ledgers; this report's Part VI carries the same corrections sprint by sprint and adds nothing
  from their table.
- Their reading of `docs/FINDINGS.md:4571` names two benchmark targets; not repeated here.
