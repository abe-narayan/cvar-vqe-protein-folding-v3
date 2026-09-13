# PREREG_C2 -- THE PRIOR-INPUT LADDER (lane P, Sprint 26)

Written 2026-09-13 00:40, before any endpoint number exists. Code: `s26/p_ladder.py`
(unit-tested on synthetic data by `s26/p_ladder_test.py`; one-fold training probe under
`s26/jobrun.py`). Results will live in `s26/results/p_ladder_<rung>_s<seed>.json`. Never edited
after the first endpoint run; addenda are appended.

## 0. What the record already knows, so nobody re-discovers it

- The prior's derivative is the only steep lever: -2.1496 A per unit gamma at the origin,
  gamma = 0.0225 reaches 3.0 A (S24 L13, `s24/results/priorladder.json`), but ONLY along the
  native's own direction; a real operator at cos 0.5 travelling 25% is worth +0.024 A (S25 L12).
- Every re-reading of the shipped posterior is closed (S25 L2/L6/L7/L12: width, location, mode,
  functional, metric projection, de-quantisation; 51 arms, corr(progress, endpoint) = +0.054).
- ESM-pca32 buys -0.288 A [-0.484, -0.092] over one-hot on selection (S7-11, n=126, p=0.0046);
  pca128 -0.337, raw -0.301, mutually indistinguishable. The per-target artefact
  (`s7/repr_tune.json`) NO LONGER EXISTS on disk or in git history (verified 2026-09-13; only
  `s7/repr_oracle.json` survives), so those numbers cannot be re-read per target and this ladder
  is also their only possible reproduction.
- More fragment data makes the predictor monotonically worse (S7-2); the ranking-loss family is
  worse (S7-8); the triangle operator on one-hot input buys diagnostics and no selection gain
  (S7-9); PairNet on the SAME inputs as the deployed model is -0.080 [-0.242, +0.082] through a
  distance-geometry fit (S19 L11 `s19/results/a_models.json`: paired sd 0.910, MDE 0.227) -- so
  S7-11's "tri on ESM input has not been measured" is superseded on the deployed-fit basis and
  still unmeasured on the pipeline basis (top-75 average, built chain), which is what this ladder
  measures.
- The capacity curve of learned in-band rankers is monotone decreasing (S9-8); de-noising the
  training labels is worth 0.044 A (S8-14); the pool-histogram/typicality direction is closed on
  selection: ranking by the pool's own mean profile selects 3.734 against 3.454 (S7-3 `gain 0`),
  the per-pair Bayes combination of pool distances with the distogram is monotone the wrong way
  (S10-2), and every "move toward the pool" arm is worse than the incumbent (S19 L14).
- The fold models are not independent: each saw ~100 of the other 125 dev natives (S24 L8), so
  the fold-clustered CI is the deciding interval.

## 1. Hypothesis and exact falsifier

**H_C2.** A distogram trained on inputs the machine can compute today -- the ESM-2 650M
representation at higher dimension, the raw 1280-d representation with a learned projection, a
different model size, more MLP capacity, the triangle update on ESM input, or a calibrated mixture
with the retrieval pool's own histogram -- is a better prior than the shipped one, and the gain
survives the pipeline.

**Falsifier (per rung).** A rung is a rung ONLY if, on the BUILT CHAIN (`arm`,
`s12.instrument.project` of the shipped top-75 average, the 3.2148 A basis), it beats the shipped
posterior scored through the identical code path by more than its own MDE (`ST.compare`,
MDE = 2.8016 x SE), with the fold-clustered CI excluding zero AND 5/5 folds the same sign. If no
rung clears that, H_C2 is refuted at this ladder and Proposal C's "learn a better prior" reduces
to "no input the machine can compute today does it"; the write-up will state the power of every
rung (MDE and effect/MDE), never "null" for an underpowered rung.

**Secondary falsifiers, reported not decided on.** The same test on the point cloud (`cloud`,
3.0483 basis) and on selection (`sel`, argmin over K=500, tie-averaged, 3.4540 basis). A rung that
wins on `sel` and not on `arm` is a better FILTER, not a better prior for this pipeline (S17 L23's
distinction), and is reported as such.

**Reproduction gates (native-free, before any endpoint is read).**
1. `shipped` through this module reproduces 3.4540 / 3.0483 / 3.2148 / 3.2041 (sel / cloud /
   arm / fit) to 4 decimals; `I.project` already reproduces the persisted `fit_ca` and `ca` of
   1A13 at max abs 0.0 (`s26/results/p_probe_esm.json`).
2. Rung `pca32` retrained from the compact tables reproduces the pinned
   `distogram_models/fold<f>_esm_frag.pt` (`p_ladder.py gate`): risk-table max abs diff, argmin
   agreement over the K=500 pool, top-75 overlap. If the retrain is not bit-exact (torch CPU
   nondeterminism, float32 rounding of `s5/esm32.npz`), the gate falls back to endpoint
   equivalence within MDE and the ladder is read against the RETRAINED pca32 as well as the
   pinned posterior, so that a "gain" can never be a retraining artefact.
3. Rung `mix` at lam = 0 reproduces the shipped risk table bit-exactly per target (asserted).

## 2. Rungs, cheapest first, and the recorded prior for each

| rung | what changes | corpus / training | expected (from the record) |
|---|---|---|---|
| shipped | nothing; the pinned posterior | -- | the identity |
| noesm | ESM block removed (42-d physicochemical pair features) | shipped MLP, retrained | WORSE by +0.29 to +0.34 on sel (S7-11); built chain unknown |
| pca32 | shipped inputs from `s5/esm32.npz` + `s7/repr_cache/esmcon.npz` | shipped MLP, retrained | reproduces the pinned model |
| pca32f | 32-PCA refit PER FOLD on training sequences (no global-PCA leak) | lean trainer | null vs pca32 (isolates the fitting set) |
| pca128 | 128-PCA per fold (`s7/repr_cache/pca128_fold<f>.npz`) | lean trainer | null vs pca32 (S7-11: 3.417 vs 3.466, inside MDE 0.28) |
| raw | full 1280-d, first Linear layer is the projection | lean trainer | null vs pca32 (S7-11: 3.453) |
| esm8m | esm2_t6_8M reps + its contact head (the downward size rung) | lean trainer | worse than pca32; sizes the ESM channel |
| wide | width 768, depth 4 on shipped inputs | shipped trainer | null-to-worse (S9-8 capacity curve) |
| pairnet | pinned `pairnet_models/fold<f>_c64b4_s0.pt` on ESM input | none | null (S19 L11 -0.080 [-0.242,+0.082]) |
| mix | (1-lam) shipped + lam pool histogram, lam leave-fold-out | none | lam* = 0; lam = 1 is typicality, +0.28 on sel (S7-3) |

The 3B ESM-2 rung is NOT in the ladder: `esm2_t36_3B_UR50D.pt` is not on disk, is 5,678,116,398
bytes to download (HTTP HEAD 2026-09-13), and needs >= 5.7 GB resident in fp16 against ~4.4 GB of
campaign headroom (`s26/results/b1_feasibility.json`). The 8M model IS on disk
(`~/.cache/torch/hub/checkpoints/esm2_t6_8M_UR50D.pt`, 30 MB) and gives the size axis a second
point below 650M.

## 3. Comparison arm, basis, readout, statistics

- Comparison arm: the shipped posterior (`s12/cache/disto_<pdb>.npz`) scored through the SAME
  `score_target -> endpoint` path as every rung, paired per target.
- Basis: built chain PRIMARY (`arm`); point cloud and selection carried on both sides. Stated on
  every line. `fit` (the lam = 0 chain) is carried and never decided on.
- Readout: the shipped top-75 uniform coordinate average in the medoid frame, m = 75 fixed.
- Statistics: `ST.compare(rung, shipped, folds=ST.pinned_folds(pdbs))`; iid and fold-clustered
  CIs, MDE, W/L/ties, median beside mean, the uniform-effect concentration null; verdict from the
  fold CI and the MDE gate only. The mixture's leave-fold-out lam is chosen on `arm` within the 4
  training folds and applied to the 5th; the per-target oracle over lam is reported beside
  `ST.best_of_k_within` and its split-half transfer, never as a result.
- Prior-side currency: for every rung, `gam_eff` and `cos` toward the native in the S24 MASS
  construction (probability space) AND in location space (S25 loc), reported side by side;
  -2.1496 x gam_eff is never quoted alone. MAE is a diagnostic column only.
- Every positive result is re-run at `--seed 1` and with the fold-processing order reversed
  (`--folds 4,3,2,1,0`; fold labels stay pinned) and must replicate within its own CI.

## 4. Expected effect against the computed MDE

Computed from persisted per-target arrays (this file's own numbers, 2026-09-13):

| endpoint | proxy comparison | paired sd | SE | MDE |
|---|---|---|---|---|
| point cloud, a prior change | `s24/results/priorladder.json` MASS0.1 - MASS0.0 | 0.352 | 0.0313 | 0.0878 |
| point cloud, a small prior change | MASSFIXW0.1 - MASS0.0 | 0.207 | 0.0184 | 0.0517 |
| built chain vs cloud (projection noise) | `bench_results/cache/1fc9f2dcf489e2fb` rmsd_arm - rmsd_avg | 0.205 | 0.0182 | 0.0511 |
| distogram swap, distance-geometry fit | `s19/results/a_models.json` pairnet - deployed | 0.910 | 0.0811 | 0.2271 |
| selection, ESM vs one-hot | S7-11 CI [-0.484, -0.092] | -- | 0.100 | 0.280 |
| selection, oracle matrix vs shipped | `s7/repr_oracle.json` - `shipped` | 1.373 | 0.1224 | 0.3428 |

So a distogram swap should carry an MDE of roughly 0.09 to 0.20 A on the built chain and about
0.28 A on selection. The record's prior for every ESM-input rung is |effect| < 0.05 A (S7-11's
mutual indistinguishability), i.e. well UNDER the MDE: the likely outcome is "not measured", and
the sprint gets the power statement, not a rung. For the noesm control the record predicts a
+0.3 A degradation on selection (about 1x MDE, Type-M zone) and an unknown built-chain effect.
gamma-equivalent expectation: the pinned posterior is the origin; a rung that reached gamma =
0.0225 at cos 1 would be worth -0.048 A on the point cloud, which is BELOW the point-cloud MDE of
a prior change (0.088) -- a rung can be real and invisible at n = 126, and the write-up will say
so with the gam_eff/cos numbers.

## 5. Memory and agent-hours (measured, not guessed)

- Compact-table coverage: `s5/esmraw.npz` (487 MB float32 resident, 1.3 s load),
  `s5/esm32.npz`, `s7/repr_cache/esmcon.npz` cover all 6,790 training sequences and all 126
  targets; `esm_small.npz` covers only 360 (`s26/results/p_probe_esm.json`, jobrun peak 0.594 GB).
  `esm_cache.npz` costs 1.79 GB resident and 19.7 s (`s26/results/p_probe_esmcache.json`) and is
  never loaded by the ladder.
- Training arrays: 552k-557k pairs per fold; materialised 183-d X is 0.40 GB (+0.40 GB normalised
  copy) -> ~1.5 GB job; lean rungs gather from the residue table -> ~1.0 GB. The one-fold training
  probe (`s26/jobs_done/p_probe_train_pca32_f0.json`) supplies the wall time and peak RSS; both are
  appended below when it finishes.
- Eval: `I.project` 2.9 s per target -> ~7 min per rung for 126 targets; mix at 8 lam values ~55 min.
- Agent-hours: ~2-3 h machine time for 6 trainable rungs x 5 folds, ~2 h eval, ~2 h write-up.

## 6. Rule 0 -- six operator forks, each naming the alternative not taken

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | the shipped Bayes-risk score from a GENUINE `core.predict.Distogram` built from each rung's posterior (`w` recomputed self-consistently, as a better distogram would) | re-implementing the risk table; any reweighting of pairs (closed S24 L6, S25 L2); a new functional |
| basis | built chain primary; point cloud and selection carried, stated on both sides | MAE, NLL, calibration slope or in-band rho as the decision metric (S7-3, S7-9: proxies do not price selection) |
| readout | shipped top-75 uniform average, m = 75 fixed for every rung | re-tuning m per rung (S12 says m* shrinks as the objective improves, so a genuinely better prior would be UNDER-read at m = 75; carried as `top75_best`/`top75_mean`, and the m-ladder is a pre-declared FOLLOW-UP only for a rung that clears the falsifier) |
| normalisation | identical corpus, order, epochs, batch, optimiser, schedule, seed and loss for every trainable rung; per-fold PCA for pca32f/pca128/esm8m, the shipped global PCA for pca32 (that is what production uses) | per-rung hyperparameter search; a global PCA for the new rungs (it lets held-out sequences shape the projection) |
| null | the shipped posterior through the identical path (the identity); the pinned-model reproduction gate; `noesm` as the zero-ESM control; `mix` lam = 0 asserted bit-exact | comparing against the pinned constants alone; a permuted-feature control (the endpoint's noise is already the shipped arm's) |
| label | continuous CA-RMSD of the built chain; `sel` and `cloud` beside it; gam_eff with cos beside it | any binarised "wins", MAE, or gam_eff quoted as -2.1496 x gam_eff |

## 7. What I will do if a rung clears the falsifier

Re-run at seed 1 and with the reversed fold order; run the m-ladder (m in {20, 35, 50, 75, 110})
on that rung as a pre-declared secondary; report gam_eff and cos; then hand the rung to the
Examiner for the operator-fork re-enumeration before it is called a result.

## ADDENDUM 1 (2026-09-13 00:45) -- the pre-sign-off probes, appended, nothing above edited

- Synthetic unit tests `s26/p_ladder_test.py` (jobrun `p_ladder_test`, 10 s, peak 0.325 GB): lean
  trainer == `core.predict.MLP.fit` at max |dprob| 1.0e-7; rebuilt pca32 features vs the shipped
  construction on 12 targets max abs 9.5e-7 (13.8% of entries differ by one float32 ulp);
  noesm features exact; mix lam = 0 exact and `risk_table` reproduces the cached shipped risk at 0.0;
  NaN-poison bit-identical; the guard refuses unknown sequences. ALL OK.
- One-fold training probe, rung pca32, fold 0 (jobrun `p_probe_train_pca32_f0`): 6,630 chains,
  552,199 pairs, 183-d, 40 epochs, **453 s wall, peak RSS 1.248 GB**
  (`s26/jobs_done/p_probe_train_pca32_f0.json`). So the exact-path rungs cost ~38 min per rung
  for 5 folds inside a 1.5 GB job.
- Reproduction gate 2 (`s26/results/p_ladder_gate_pca32_fold0_s0_probe.json`, 25 fold-0 targets,
  native-free): input features vs shipped max abs 1.26e-6; posterior max abs 1.5e-5; risk table
  max abs 6.8e-4 (the pinned model vs the s12 cache itself differs by 2.2e-5); **K = 500 argmin
  identical on 25/25; top-75 overlap 1.000; score Spearman 0.99999998.** Not bit-exact
  (float32 rounding of the cached embeddings), equivalent in every quantity the pipeline consumes.
- Rung `conly` (PREREG_B2) added to the ladder after the probe job exited: 55-d, no embedding block.
- Fold-processing order: `--folds 4,3,2,1,0`; seed replication: `--seed 1`.

## ADDENDUM 2 (2026-09-13 08:45) -- TRAINING ORDER FIXED BEFORE THE FIRST JOB, per LEDGER L16b

Session note: the lane was cut by the API session limit at ~00:46 and resumed at 08:37; no
training job had been launched before the cut (`s26/jobs_done/` holds only the probes; the only
model on disk is `models/p_ladder/pca32_fold0_s0_probe.pt`). Nothing is restarted; the probe
model is kept as the probe and the ladder's own pca32 fold 0 is trained afresh under the
declared seed/tag so that the gate compares like with like.

Rung TRAINING order, one governor job at a time, est-ram 1.5 GB (raw 1.8 GB), NT = 2, one .pt
per fold as the checkpoint, each job resumable (existing fold checkpoints are skipped):

    1 noesm  ->  2 conly  ->  3 pca32  ->  4 wide  ->  5 pca32f  ->  6 pca128
    ->  7 featurise-esm8m then esm8m  ->  8 raw (last: its first Linear layer is 5175 x 384)

Fold order inside every job: 0,1,2,3,4 (the replication run uses 4,3,2,1,0 and seed 1). pairnet
and mix train nothing. No `eval`/`report` runs until "PHASE 0 SIGNED OFF" is in `s26/LEDGER.md`
(enforced in code). Evaluation order after sign-off: shipped (the identity gate), then the rungs
in the order their checkpoints completed.
