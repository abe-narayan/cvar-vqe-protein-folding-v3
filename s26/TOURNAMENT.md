# S26 TOURNAMENT (lane A, the Adversary)

Ranked 2026-09-13 over every `s26/IDEA_*.md` on disk (16 files: `ls s26/IDEA_*.md`). For each
idea: a kill attempt against the record (sprint and ledger line cited; a "Closed by S<n> L<m>"
counts only if the closing experiment measured the SAME operator in the SAME space, and I say
when it did not), then a score of expected information (0-3) x plausibility (0-1) / cost. Cost
is agent-hours plus governor load; memory is priced tonight against the coordinator's ~1 GB
free-headroom figure, and an idea whose peak RSS exceeds ~1 GB is DEFERRED, not killed.

Standing facts used as kill ammunition (memory index + ledgers): the distance prior is the
ceiling and its derivative is the only steep lever (S24 L13); the pool's error is 68%
common-mode (S23 L9); selection is closed at four levels and the terminal operator consumes the
set mean (S12, S17); both physics energies rank worse than a random subset and the functional
lever is closed in five forms (S25 L16/L18, S24 L16); AMBER relaxation's gain dissolves under a
matched random move (S16 L27, reconfirmed on the production input this sprint, L39: +0.0111 over
random, Type-M); consensus is outlier avoidance (S23 L5); routers fail, two harmful held out
(S22 L7, S23 L7); grid oracles are order statistics, quote `best_of_k_within` and the split-half
transfer (S24); MDE is per comparison; tie-breaking leaks the pool order (S24); the sequence
channel is near-blind (phi 36.1 vs 36.4, `s26/results/a_c26_phi_mae.json`); no fresh benchmark
exists (204 clusters spent); controls must match the operator's space and a zero-information
control must be plausible (constant helix, not uniform); concentration is wrong when
discrimination binds; the shared-referent floor; unstated operators align with the builder's
hypothesis.

## 1. Kill attempts

None of the 16 has a clean same-operator kill: in every case the closing experiment the record
offers measured a different operator or a different space. So all 16 survive as scored
hypotheses; the low-plausibility ones are deprioritised, not deleted (rule 11: never move an
open item to closed on one failed experiment, and never soften a proposal to make it survive --
the mirror applies, do not manufacture a kill).

| idea (lane) | strongest closer on the record | why it is not a kill | residual plausibility |
|---|---|---|---|
| amber_prior_partner (P) | functional lever closed in 5 forms, S24 L16 / S25 L18; rank-permuted AMBER -0.0010 at 0.08x MDE | all five forms act on candidate RANKING; this acts on the PRIOR's per-pair location/width, a different operator space | 0.10 |
| amber_reject (PH) | S19 Q3 fixed-count reject null-to-harmful; S24 L10-A fixed-fraction harmful | neither used a physics-set count with refill on the built chain -- but **L43 this sprint fired the falsifier the other way (+0.228 worse, point cloud)**; the built-chain arm remains | 0.05 (all but closed by L43) |
| better_prior_inputs (P) | every re-reading of the posterior closed, S25 L2/L6/L7/L12; S17 L23 ESM contact scalar mostly compactness | S17 L23 used the contact SCALAR, never the attention MAPS; retrieval histogram as a training INPUT (to learn the 68% common mode) is untested | 0.15 |
| branch_select (PH) | S25 Form 5 per-target physics sign closed at 0.14x MDE | Form 5 switched near-identical arms (corr 0.95); the projection branches differ up to 1.6 A -- a materially different discrete choice | 0.20 |
| cis_peptide (PH) | -- (DONE: L22 census, L38 floor) | delivered; 0 cis on the instrument, floor 0.347 A tracks omega non-planarity, worth 0.000 A here | n/a |
| coherence_penalised_training (P) | S19 L14 "the estimator is made of the bias"; S7-8 ranking/anti-shrinkage losses failed | S19 L14 is about native-free inference; a training loss sees the native, and no arm put a coherence term in the loss | 0.15 |
| conformational_identity_floor (W) | S24 L4 (n=4 cross-deposit values) | never generalised past n=4; this is the n=18/22 measurement | 0.70 (descriptive) |
| l17_target_dependent_hamiltonian (Q) | S25 L17 constant spectrum, two trained states | L17 ran no target-dependent variant; this is that variant (order-preserving, matched entropy) | 0.15 |
| product_state_optimum (Q) | S25 measured 0.902 nats divergence, endpoint tie on selection | S25 did not identify WHY (the target is a product state) nor run the 7-parameter arm on the built chain | 0.80 (the fact is proven; the endpoint arm is cheap) |
| rotamer_relief (PH) | -- (Part A DONE, L23: 96.8% side-chain) | Part A survived; Part B (the relief single point as a ranker/reject) is untested | 0.15 accuracy / 0.60 that relief removes the mass |
| selfcopy_proxy_bound (W) | S10-4 dev +0.0004, S9-10 bench +0.0030 (document-only, C27) | never separated exact from partial copies, never on the built chain, never channel B -- **DELIVERED this sprint, L44: MINOR, 0.028 A envelope** | n/a (done) |
| strain_difficulty (PH) | routers failed held out, S22 L7 / S23 L7 | routers never used the relaxation's own post-hoc energies; this is a calibration flag, not an accuracy lever | 0.50 signal / 0 RMSD |
| tiebreak_noise_floor (W) | S17 L628 oracle-best under tie-breaks sd 0.018 | that is the oracle best of 500, not the emitted chain; the production endpoint's tie-break floor is unmeasured | 0.85 (measures an MDE) |
| trainability_paper (Q) | -- (writing; outline done) | open item 3; not an experiment | n/a |
| window_ensembling (P) | S17 L12 widening K hurts | ensembling at FIXED K never widens a shortlist, so S17's displacement mechanism cannot act; but 3 same-provenance keys give parallel clouds (bias cos 0.943, S24 L3) | 0.15 |
| window_provenance (W) | uniform readout optimal over two families, S25; 68% common-mode | provenance is not a function of rank-power or distance-to-medoid; the class census is new | 0.15 |

## 2. Ranked table (survivors that still have something to run)

Score = information x plausibility / cost-in-agent-hours (compute is governed, so agent-hours is
the binding cost); memory gate applied separately. Ideas already DONE or DELIVERED
(cis_peptide L22/L38, selfcopy_proxy_bound L44, the amber_reject point-cloud arm L43,
trainability outline) are listed in section 4, not ranked here.

| rank | idea (lane) | info | plaus | cost (ah / peak RSS) | score | note |
|---:|---|---:|---:|---|---:|---|
| 1 | product_state_optimum (Q) | 2 | 0.80 | 1 / rides A1, ~0.4 GB | 1.60 | the endpoint arm is inside the A1 harness; near-free, high-plausibility presenter sentence |
| 2 | conformational_identity_floor (W) | 2 | 0.70 | 0.5 / <0.1 GB | 2.80 | Part E of the selfcopy prereg, code written, gated; one report sentence ("identity buys X A") |
| 3 | tiebreak_noise_floor (W) | 2.5 | 0.85 | 3 / <0.5 GB | 0.71 | measures the pipeline's own noise floor; contextualises every hundredths-level effect (AMBER relax -0.022, ORACLE weight 0.015, C27 +0.0004) |
| 4 | strain_difficulty (PH) | 1.5 | 0.50 | 1.5 / <0.3 GB | 0.50 | a native-free confidence flag; the cheapest remaining difficulty predictor; reads 126 cached records |
| 5 | branch_select (PH) | 2 | 0.20 | 3 / ~1.0 GB AMBER | 0.13 | first accuracy role for physics on a small discrete choice; AMBER, one job at a time |
| 6 | window_ensembling (P) | 2 | 0.15 | 2 / <0.6 GB | 0.15 | fixed-K ensembling; the one form widening-K did not test |
| 7 | l17_target_dependent_hamiltonian (Q) | 2 | 0.15 | 5 / <0.4 GB | 0.06 | mandatory direction; A3, in Q's queue; expected null by the entropy mechanism |
| 8 | rotamer_relief B (PH) | 2 | 0.15 | 3 / ~0.8 GB AMBER | 0.10 | the relieved single point as a ranker; decides if all-atom scoring of windows is ever worth it |
| 9 | window_provenance (W) | 1 | 0.15 | 2 / <0.3 GB | 0.075 | census-first; a one-line pool-composition fact |
| 10 | amber_prior_partner (P) | 2 | 0.10 | 2 / <0.6 GB | 0.10 | AMBER as a prior mixture partner; the last AMBER form |

DEFERRED on memory (>1 GB tonight; NOT killed, run when the box frees up):
- better_prior_inputs `attn` (P): the 650M model forward with head weights, 3.0-3.5 GB, needs
  the box to itself and a probe first. Its `ragp` rung (<1 GB) is not deferred and folds under C2.
- coherence_penalised_training (P): training peak 1.248 GB (`s26/jobs_done/p_probe_train_pca32_f0.json`),
  just over the 1 GB line; defer until the C2 ladder training frees the box.

## 3. Survivors in run order

Cheap-and-high-plausibility first, then the mandatory directions already in a lane's queue, then
the low-plausibility closers, memory-deferred last:

1. product_state_optimum (Q) -- rides A1, no extra compute.
2. conformational_identity_floor (W) -- gated, <1 min, the moment W's endpoints run.
3. strain_difficulty (PH) -- reads 126 cached records, <1 min, native-free correlation.
4. tiebreak_noise_floor (W) -- 8 tie-break draws x 126, ~50 min CPU, checkpointed.
5. l17_target_dependent_hamiltonian (Q, A3) and branch_select (PH) -- both endpoint, one AMBER
   job at a time for branch_select.
6. window_ensembling (P), rotamer_relief B (PH), window_provenance (W), amber_prior_partner (P)
   -- run only as capacity allows; all low plausibility.
7. Deferred until RAM frees: better_prior_inputs `attn`, coherence_penalised_training.

The amber_reject built-chain arm (PH) still owes its endpoint after the point-cloud arm fired
harmful (L43); it is the falsifier's last leg, not a survivor to schedule fresh.

## 4. Done / delivered (not ranked)

- cis_peptide (PH): L22 (census, 0 cis on the instrument), L38 (floor 0.347 A).
- selfcopy_proxy_bound (W): L44 (MINOR, 0.028 A envelope; C27 +0.0004 re-derived).
- amber_reject point-cloud arm (PH): L43 (falsifier fired the other way, +0.228 worse).
- rotamer_relief Part A (PH): L23 (96.8% side-chain-involving).
- trainability_paper (Q): outline delivered (`s26/TRAINABILITY_PAPER_OUTLINE.md`).

## 5. Orphaned survivors for the Wildcard lane (W)

Lane P is saturated (the C2 ladder is a 6-10 h training chain plus B3/C4/C5), so its
tournament ideas beyond the C-series are orphaned and go to W in this order (highest score
first), memory permitting:

1. window_ensembling (P) -- <0.6 GB, 20 min, 2 ah; the cheapest orphan and the mandatory
   test-time-ensembling direction, so it should be run.
2. amber_prior_partner (P) -- <0.6 GB, low plausibility; run only after window_ensembling.
3. better_prior_inputs `ragp` rung (P) -- <1 GB; the `attn` rung stays deferred on memory.
4. coherence_penalised_training (P) -- DEFERRED on memory (1.25 GB); W picks it up only when
   the box frees up.

W's own three ideas (conformational_identity_floor, tiebreak_noise_floor, window_provenance)
are ranked in section 2 and are not orphans; W runs them before taking a P orphan.
