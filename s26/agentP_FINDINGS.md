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
