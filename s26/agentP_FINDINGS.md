# agentP_FINDINGS -- lane P (Prior / Learning), Sprint 26

Started 2026-09-13 00:50. Format: S12-S25 tiers. Every number carries its artefact path. Nothing
here is an endpoint result: the phase gate ("PHASE 0 SIGNED OFF") has not been posted, and no
RMSD to a native has been computed by this lane.

## 0. What this lane owns and where it stands

Proposal C (C1-C5), Proposal B (B1-B3), every learned-model idea in the tournament. Written
before sign-off: `PREREG_C1..C5.md`, `PREREG_B1..B3.md`, `IDEA_better_prior_inputs.md`,
`IDEA_window_ensembling.md`, `IDEA_coherence_penalised_training.md`, `p_ladder.py` (C2 code),
`p_ladder_test.py` (synthetic tests, all pass), `p_probe_esm.py`, `p_probe_esmcache.py`,
`results/b1_feasibility.json`, `results/p_probe_esm.json`, `results/p_probe_esmcache.json`,
`results/p_ladder_gate_pca32_fold0_s0_probe.json`, `models/p_ladder/pca32_fold0_s0_probe.pt`.

## 1. DEMONSTRATED (native-free measurements, complete)

**D1. The compact ESM tables cover the whole training corpus; the 1.5 GB bank is not needed.**
`results/p_probe_esm.json` (jobrun `p_probe_esm`, peak RSS 0.594 GB, 10 s): `s5/esmraw.npz`,
`s5/esm32.npz` and `s7/repr_cache/esmcon.npz` each hold all 6,790 unique training sequences
(787 peptides + 6,003 fragments) and all 126 targets; `esm_small.npz` holds 360. `s5/esmraw.npz`
materialises to 486.6 MB float32 in 1.3 s. Per fold: 6,609-6,640 chains, 552,199-556,709 pairs.
`results/p_probe_esmcache.json` (jobrun `p_probe_esmcache`, tag ESM): `esm_cache.npz` holds
22,795 sequences (1280-d, mean length 14.0), costs 1.79 GB resident (process peak_wset 1.793 GB;
the jobrun 5-s sampler read 1.456 GB, an under-read of the transient) and 19.7 s. The ladder never
loads it (`p_ladder.install_guard`).

**D2. `s12.instrument.project` reproduces the production emission bit-exactly, in 2.9 s.**
1A13: `fit_ca` and `ca` rebuilt from the persisted `avg_ca` differ from the cached emissions by
max abs 0.0 (`results/p_probe_esm.json`). The built-chain endpoint costs ~6 min per rung.

**D3. The C2 trainer reproduces the pinned distogram.** Rung pca32, fold 0, retrained from the
compact tables with `core.predict.MLP` (40 epochs, 552,199 pairs, 183-d): 453 s, peak RSS
1.248 GB (`jobs_done/p_probe_train_pca32_f0.json`). Against `distogram_models/fold0_esm_frag.pt`
on the 25 fold-0 targets (`results/p_ladder_gate_pca32_fold0_s0_probe.json`, native-free): input
features max abs 1.26e-6 (24% of entries differ by one float32 ulp because the cached embedding
is float32), posterior max abs 1.5e-5, risk table max abs 6.8e-4 (the pinned model vs the s12
cache itself: 2.2e-5), K = 500 argmin identical 25/25, top-75 overlap 1.000, score Spearman
0.99999998. Synthetic tests (`p_ladder_test.py`, jobrun 10 s, 0.325 GB): lean trainer ==
`MLP.fit` at 1.0e-7; `mix` lam = 0 bit-exact; NaN-poison bit-identical; the guard refuses unknown
sequences.

**D4. B1 is infeasible on this box (three independent blockers).** `results/b1_feasibility.json`:
no ESMFold/3B weights on disk; `import esm.esmfold.v1.pretrained` raises `ModuleNotFoundError:
omegaconf`, and `openfold` (needed by esmfold.py:11-13 and trunk.py:11) is absent; the README's
openfold path needs nvcc and Python <= 3.9 (box: 3.13, CPU torch 2.13); checkpoints
2,771,653,574 + 5,678,116,398 bytes (HTTP HEAD, network reachable in < 10 s); fair-esm halves the
3B LM and keeps the trunk fp32 (esmfold.py:43-46) -> 8.76 GB resident minimum against 4.4 GB of
campaign headroom (governor: 16.75 GB, 64.8% used, ceiling 93%). Nothing was downloaded.

**D5. B3's per-target artefacts exist.** `s13/cache/tors_rows.npz['a_pepPos']` (126 rows,
`rmsd_build` mean 3.7705) and `s14/results/ladder.json per_target['L0_constant_helix']` (126,
mean 4.0648). Paired against the production `rmsd_arm` (3.2148): tors - arm +0.5557 (median
+0.2314), SE 0.1318, MDE 0.3694, 37W/89L, top-10 share 0.567; helix - arm +0.8500 (median
+0.1691), SE 0.1494, 32W/94L; FAIL18: tors 5.569 / arm 6.032 / helix 5.887; other108: 3.471 /
2.745 / 3.761. The per-target min over {arm, tors, helix} (2.9632, -0.2515 vs arm) is 92%
accounted for by `ST.best_of_k_within`'s across-target null (k_eff 2.34): an order statistic,
as S24 L15-A found for the 19-arm panel. (Read from persisted artefacts; no structure built.)

**D6. One S7 artefact is lost.** `s7/repr_tune.json`, the per-target array behind S7-11's
"ESM-pca32 buys -0.288 A", is not on disk and not in git history (only `s7/repr_oracle.json`
survives: native-matrix ranking 1.9938, pool best 1.7108). The ESM contrast can only be
re-measured, which the ladder's noesm / pca32 rungs do.

## 2. HYPOTHESIS (pre-registered, not run)

PREREG_C2 (the ladder), PREREG_B2 (conly / esm8m size axis), PREREG_B3 (characterisability),
PREREG_C4 (routers on unused features), PREREG_C5 (common-mode subtraction), PREREG_C3 (physics
prior partner; the coordinator's C3 for lane PH is a different item, see the PREREG's addendum).
The record's prior for every one of them is a null; the MDEs are computed in PREREG_C2 section 4.

## 3. Corrections to the record noticed while reading (not new measurements)

- S7-11 says "tri on ESM input has not been measured". S19 L11 measured PairNet with the deployed
  inputs (which include ESM): -0.080 [-0.242, +0.082] through the distance-geometry fit,
  paired sd 0.910, MDE 0.227 (`s19/results/a_models.json`). Still unmeasured through the
  pipeline's readout; that is rung `pairnet`.
- The S13 torsion arm's `base` column (3.2126) is the leaderboard rebuild and differs from the
  production `rmsd_arm` by up to 0.171 A per target; B3 pairs against the production cache.

## 4. What damaged my own expectations

- I expected the pca32 retrain to be bit-exact or chaotic; it was neither: 1.5e-5 in the
  posterior and identical selection. Training here is stable to input perturbations of 1e-6.
- I expected `esm_cache.npz` to cost 3-4 GB; it costs 1.79 GB. The bank was avoidable anyway.
- I expected the arm-choice oracle over {pipeline, sequence-only, helix} to carry per-target
  signal on FAIL18; the across-target null accounts for 92% of it.

## 5. What I did not do and why

- No endpoint experiment (phase gate). No full-ladder training either: the contract allows
  1-target probes before sign-off and full 5-fold training is more than that; the question is
  logged for the coordinator (end of turn).
- Did not run `featurise-esm8m` (a corpus-wide forward pass; same reason).
- Did not download anything (rule: nothing above 100 MB without a probe and the coordinator's OK).
- Did not edit `p_ladder.py` while its probe job was running; `conly` was added after it exited.
- Did not write PREREG_C3 to the coordinator's definition (lane PH's brief was not visible to
  me); the file says so in its addendum.

## 6. Artefact index

`s26/results/p_probe_esm.json` · `s26/results/p_probe_esmcache.json` ·
`s26/results/b1_feasibility.json` · `s26/results/p_ladder_gate_pca32_fold0_s0_probe.json` ·
`s26/jobs_done/{p_probe_esm,p_probe_esmcache,p_ladder_test,p_probe_train_pca32_f0}.json` ·
`s26/logs/p_*.log` · `s26/models/p_ladder/pca32_fold0_s0_probe.pt`

## 7. C1 -- CLOSURE REPRODUCTION FROM PERSISTED ARTEFACTS (PREREG_C1; reads recorded numbers only)

Both closures that Proposal C rests on reproduce from the artefacts on disk to the third
decimal. (a) The set-transformer over the full signed deviation map (S12 agg_FINDINGS section 6)
has a flat learning curve: `s12/results/agg_dec_v2_n8.json` mean_raw 3.0433, `_n16` 3.0491,
`_n32` 3.0456, `_n64` 3.0358, `agg_dec_v2.json` (3 folds, ~75 targets, early-stopped) 3.0258;
the same harness with the RMSD label leaked: `agg_dec_v2_oracle_n8.json` 2.5342,
`_oracle_n32` 2.2004, `agg_dec_v2_oracle.json` 2.1460 (v1 harness: 3.0561 -> 3.0640 real;
2.5335 -> 2.1466 leaked). The real-feature range over an 8x change in training data is 0.023 A,
inside the seed spread (3.0258 vs 3.0337); the leaked label is at 57% of its final gain at n = 8.
Basis: RAW weighted coordinate average (point cloud); avg75 reference 3.0483. (b) The perfect
ranker inside the shipped top-25: `s17/results/inband.json`, 126 rows, cells['25']: band_best
2.6087 (median 2.5155), random-in-band 3.4676, distance argmin 3.4540, consensus 3.3692,
Legacy 3.4149, band mean 3.5016; band_best - random -0.8589, SE 0.0600, MDE 0.1680, 126W/0L.
Basis: selected single candidate (selection). Neither number is a built chain.

| closure | quantity | recorded | reproduced | artefact |
|---|---|---|---|---|
| set-transformer, real features | n=8 / 16 / 32 / 64 / ~75 | 3.043 / 3.049 / 3.046 / 3.036 / 3.026 | 3.0433 / 3.0491 / 3.0456 / 3.0358 / 3.0258 | `s12/results/agg_dec_v2*.json` |
| set-transformer, leaked label | n=8 / 32 / ~75 | 2.534 / 2.200 / 2.146 | 2.5342 / 2.2004 / 2.1460 | `s12/results/agg_dec_v2_oracle*.json` |
| perfect ranker in top-25 | band best / random / dist argmin | 2.609 / 3.468 / 3.454 | 2.6087 / 3.4676 / 3.4540 | `s17/results/inband.json` |

Verdict for PROPOSAL_C: "selection inside the pool is closed" stands as CITED AND REPRODUCED;
the S7-11 ESM number is CITED, ARTEFACT LOST (D6) and is re-measured by the C2 ladder.

## 8. After L16b: the training chain, the code for B3/C4/C5, and the hand-off format (2026-09-13 09:10)

- Session note: cut at ~00:46 by the API session limit, resumed 08:37; nothing had been
  launched, nothing restarted. The probe model `models/p_ladder/pca32_fold0_s0_probe.pt` stays the
  probe; the ladder's own pca32 is trained under the declared tag.
- Training: `s26/p_train_chain.sh` runs the eight trainable rungs one governor job at a time in
  the order fixed by PREREG_C2 addendum 2 (noesm, conly, pca32, wide, pca32f, pca128, esm8m after
  `featurise-esm8m`, raw last), est-ram 1.5 GB (raw 1.8), one `.pt` per fold, a `_p2` pass to
  resume anything the governor kills. `p_train_noesm` registered 08:52. No evaluation before
  sign-off; `p_ladder.py eval/report` refuse in code.
- Code written and self-tested on synthetic data (no RMSD read): `s26/p_stats.py` (closed-form
  ridge, nested leave-fold-out alpha, balanced accuracy, permutation null; selftest: planted
  signal 0.866 held-out balanced accuracy vs null p95 0.558, random labels 0.503, planted
  regression R2 0.961), `s26/p_b3.py` (PREREG_B3: native-free features cached to
  `results/p_b3_features.json`; classifier + regression + nulls; `run` gated), `s26/p_c4.py`
  (PREREG_C4: four new feature blocks, the S22 feature set as the harness positive control,
  routers A/B/S with label-permutation nulls; `run` gated), `s26/p_c5.py` (PREREG_C5: R1
  distance-space and R2 cloud-frame representations, GLOBAL / RIDGE / ORACLE / RANDOM-MATCHED
  arms through the same projection; selftest covers the frame round trip, rotation invariance,
  the realisation and NaN-poison; `run` gated), `s26/p_deliver.py` (the C3 stage-2 hand-off:
  `results/p_best_rung_chains.json`, phi/psi in radians as `I.project` emits them, `ca`,
  `rmsd_arm`, provenance-stamped, equality with the eval JSON asserted).
- `IDEA_amber_prior_partner.md` entered (H_C3a); H_C3b withdrawn (S8-14 arithmetic).

### What I did not do (update)
- Did not evaluate any rung, fit any router or classifier on real labels, or run C5's fit:
  all read native quantities and wait for "PHASE 0 SIGNED OFF".
- Did not edit `p_ladder.py` after `p_train_noesm` registered; the hand-off code lives in
  `p_deliver.py` for that reason.

## 9. C2 LADDER RESULTS SO FAR (2026-09-13 22:50; ledger L56, L57, L62, L63, L65, L66, L67, L72)

All rungs through one path (`s26/p_ladder.py`), paired per target against the shipped posterior;
built chain on the rebuild basis 3.2126 (L57), point cloud 3.0483, selection 3.4540; MDE per
comparison from `ST.compare`; fold CI is the deciding interval.

| rung | inputs / change | arm effect (SE, MDE, x) | fold CI | folds | sel effect (x MDE) | verdict (arm) | gam_eff prob / cos | MAE |
|---|---|---|---|---|---|---|---|---|
| pca32 | retrained shipped recipe | 0.0000 (126 ties) | -- | 5/5 | 0.0000 | IDENTITY (gate 2 passes) | 0 / 0 | 2.339 |
| noesm | no ESM block | +0.208 (0.073, 0.205, 1.01x) | [+0.099, +0.394] | 5/5 | +0.330 (1.27x) | WORSE, Type-M zone | +0.116 / 0.25 | 2.548 |
| conly | contact head only | +0.122 (0.071, 0.198, 0.62x) | [-0.035, +0.267] | 4/5 | +0.112 (0.42x) | NOT MEASURED | +0.109 / 0.24 | 2.437 |
| wide | 768 x 4 | +0.032 (0.046, 0.128, 0.25x) | [-0.030, +0.128] | 3/5 | +0.097 (0.48x; fold CI > 0) | NOT MEASURED | +0.096 / 0.18 | 2.348 |
| pca32f | per-fold 32-PCA | +0.018 (0.049, 0.138, 0.13x) | [-0.073, +0.097] | 4/5 | -0.032 (0.17x) | NOT MEASURED | +0.111 / 0.24 | 2.317 |
| pca128 | per-fold 128-PCA | +0.076 (0.063, 0.176, 0.43x) | [-0.069, +0.216] | 3/5 | +0.025 (0.13x) | NOT MEASURED | +0.117 / 0.21 | 2.386 |

Isolations (paired, same path): conly - noesm (the contact head alone) -0.086 arm (0.51x MDE),
-0.218 sel (1.00x, fold [-0.374, -0.037], 4/5); pca128 - pca32f (components only) +0.058 arm
(0.32x, 3/5).

**DEMONSTRATED.** (D7) The ESM channel is worth -0.208 A on the built chain and -0.330 A on
selection, 5/5 folds each (L62); the lost S7-11 figure is re-measured and extended to the readout.
(D8) The retrained shipped recipe reproduces the pipeline's emission on 126/126 targets (L65):
every rung difference is attributable to inputs or architecture, not retraining noise.
(D9) No rung beats the shipped prior; the closest to a gain (pca32f on the point cloud, -0.014)
is 0.11x its MDE. **The ESM-input axis is flat at n = 126 for every reduction of the 650M model
and for 2.7x the capacity, as S7-11 and S9-8 predicted.**

**HYPOTHESIS (descriptive strata, n = 18 to 48, no verdicts).** (H1) On FAIL18 every worse rung
is BETTER than the shipped prior (noesm -0.406, conly -0.513, wide -0.195, pca32f -0.187, pca128
-0.253; 10-12W of 18 each; SE 0.16-0.28) while every one is worse on other108 (+0.05 to +0.31):
S12's "on the failure class sequence conditioning is harmful" (blind 5.425 vs shipped 6.019)
reproduced at the level of the prior's inputs. FAIL18 is an ORACLE stratum (zero recall) and
this is not a router (S22 L10). (H2) The 9-10-mers carry the ESM channel's largest loss when it
is removed (+0.40 noesm, +0.19 conly) and the 13-14-mers the smallest (+0.03, +0.10).

**What damaged my expectations.** gam_eff is POSITIVE (+0.10 to +0.12 in probability space at
cos 0.18 to 0.25) for five rungs that are null-to-worse: a large move at low cosine projects
onto the truth direction and loses on the orthogonal component. S25 L12's caveat is reproduced
from the training side, on achievable rungs. I expected the sign of gam_eff to track the sign of
the endpoint; it does not.

## 10. B3 result (ledger L106, L107; `s26/results/p_b3.json`)

DEMONSTRATED (D10): the set where the pipeline beats sequence-only is NOT characterisable
native-free by the pre-registered standard. sign(arm - tors): held-out balanced accuracy 0.522
vs permutation null 95th pct 0.578 (p 0.263); sign(arm - helix): 0.557 vs 0.566 (p 0.093).
Regression of the size: R2 0.244 (vs tors; MSE reduction 0.83x MDE, NOT MEASURED) and 0.404
(vs helix; -1.142, 1.19x MDE, fold CI [-1.495, -0.700], 5/5, BETTER, Type-M zone). Descriptive
weights (ridge alpha 10, full data, standardised): ss_E -0.61 (rho -0.64), dg_sd_mean +0.39,
dg_ent_max -0.27, aa_C -0.27, dg_multimodal +0.25, con_lr +0.23; rho(ss_H, d) +0.46, rho(n, d)
-0.13. The pool's strand content sets how much a constant helix loses; nothing sets the sign.
ORACLE stratum: FAIL18 d_tors +0.463 (the torsion predictor wins), other108 -0.725.

Proposal B verdict per the campaign rule: REPLACE (`s26/PROPOSAL_B.md`).

## 11. C4, the delivery file, and the ladder's remaining rungs (2026-09-14 03:20; ledger L110, L112-L116)

- C4 (`s26/results/p_c4.json`): twelve m* routers and six s* routers on six feature blocks;
  none clears its MDE; eleven m routers point the harmful way (largest +0.068, 0.91x, worse
  than 96% of the permutation null); the S22 feature set reproduces ~0 through the same
  harness; every s router predicts s* with the wrong sign (rho -0.10 to -0.21). DEMONSTRATED
  (D11): the seventeenth router construction lands where the first five did.
- Delivery (`s26/results/p_best_rung_chains.json`, production emission, mean rmsd_arm 3.2148,
  126/126 complete, torsions in unwrapped radians): the best rung is the shipped prior.
  Cross-check `s26/results/p_best_rung_chains_rebuild_basis.json` (this lane's re-projection,
  mean 3.2126): per-target |rmsd_arm| difference max 0.171, mean 0.011, 120/126 differ by more
  than 1e-6, the L57 projection sensitivity.
- raw: folds 0-2 of 5 trained (79 min per fold at peak 1.94-1.96 GB); fold 3 queued at 03:00
  behind the job cap; the rung cannot be evaluated before the 04:30 close and is reported as
  NOT RUN (its input ablation is bracketed by pca32f and pca128, both null).
- C5: `p_c5_run` started 03:05 (alpha 0.5 only, declared in the commit); result in the ledger
  when it lands. Coherence rungs (PREREG_coherence, `p_coh.py`): NOT RUN, no slot before close.

### What I did not do and why (final)
- No replication runs (second seed, reversed fold order): no rung cleared its MDE with 5/5
  folds in the improving direction, so the replication rule never triggered.
- No `attn` rung (IDEA_better_prior_inputs): the box never emptied.
- No coherence rungs, no raw evaluation: the job cap (6-7 concurrent jobs from five lanes) and
  79-minute raw folds consumed the window.

### Correction to section 11 (per ledger L120)
The clock times "03:05" and "03:20" in section 11 were estimated, not read, and are about 80
minutes fast: C5 started at about 01:45 and section 11 was written at about 02:00; raw fold 3
was queued at about 01:25. The ledger entries' timestamps are machine-written and correct.

## 12. C5 result and the close (ledger L132, L133; 2026-09-14 03:54)

DEMONSTRATED (D12): predicting the common-mode direction native-free and subtracting it is null (GLOBAL R2
-0.009, 0.23x MDE; GLOBAL R1 +0.031, 0.60x) to harmful (RIDGE R1 +0.164, 1.79x MDE, 5/5 folds, 40W/86L; its
magnitude-matched random control +0.135, 1.33x). ORACLE DIAGNOSTIC: the true common mode subtracted is
worth -1.87 A (distance space, 120W/6L) and -3.13 A (coordinate frame, 126W/0L) on the built chain. Artefacts
`s26/results/p_c5.json`, `s26/results/p_c5_stats.json`. raw: NOT RUN (4 of 5 folds trained; L133).

## 13. What damaged my own expectations (closing)

- I expected the sign of gam_eff to track the sign of the endpoint. It does not: every null-to-worse rung has
  positive gam_eff (+0.10 to +0.13, cos 0.18 to 0.29), and PairNet's cos 0.53 with an 11 percent MAE cut buys +0.041.
- I expected the contact head to be worth little as a prior input (S17 L23 called it a filter); it carries
  two thirds of the ESM channel's selection value (-0.218 vs noesm) and the 8M model carries none of it.
- I expected the retrained recipe to differ from the pinned models by noise; it emits the pipeline's answer on
  126/126 targets, which made the ladder cleaner than planned.
- I expected the arm-choice oracle over {pipeline, sequence-only, helix} to carry per-target signal on FAIL18;
  the across-target null accounts for 92 percent of it (L14), and B3's classifier is at chance.
- I expected my own clock to be right. It ran 80 minutes fast for two hours and put a future stamp on a
  proposal (L120/L124); the fix is in the record.

## 14. What I did not do and why (closing)

- raw rung: not evaluated (fold 4 untrained; 80-minute folds at 2 GB behind a 6-7 job cap; ledger L133).
- Coherence rungs (PREREG_coherence, p_coh.py, synthetic-tested): not run; no slot before the close.
- Replication (second seed, reversed fold order): never triggered; no rung cleared its MDE in the improving
  direction with 5/5 folds.
- `attn` (IDEA_better_prior_inputs): the box never emptied; the 650M forward pass with head weights was not probed.
- C5 alpha grid reduced from three values to one (0.5) for the close, declared before the run.
- No benchmark was opened; no pinned file was written; `esm_cache.npz` was loaded once, under jobrun, as a probe (L11).
