# PREREG_B3 -- WHERE THE PIPELINE GAINS OVER SEQUENCE-ONLY, AND IS THAT SET CHARACTERISABLE NATIVE-FREE? (lane P, Sprint 26)

Written 2026-09-13, before any classifier is fitted. B3 is the repo-native version of Proposal B:
instead of a learned folding model (B1: infeasible, `s26/results/b1_feasibility.json`), the
project's own sequence-only arms are the comparator.

## The persisted per-target arms (verified 2026-09-13)

- Sequence-only torsion predictor, direct build, full coverage: `s13/cache/tors_rows.npz` key
  `a_pepPos`, 126 rows, field `rmsd_build`, mean 3.7705 (the recorded 3.770, S13 tors_FINDINGS
  section 5.1). Its `base` field (3.2126) is the S13 leaderboard rebuild of the incumbent and
  differs from the production cache's `rmsd_arm` by up to 0.171 A per target; B3 pairs against
  the production cache (`bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`, `rmsd_arm` 3.2148).
- Constant alpha-helix, zero information: `s14/results/ladder.json` `per_target['L0_constant_helix']`,
  126 entries, mean 4.0648 (the recorded 4.065); its incumbent column reproduces `rmsd_fit`
  (3.2041) at max abs 0.0000.
- Read from those artefacts (allowed pre-sign-off; no structure built):
  tors - arm +0.5557 (median +0.2314), SE 0.1318, MDE 0.3694, sd 1.480, 37W/89L, top-10 share 0.567;
  helix - arm +0.8500 (median +0.1691), SE 0.1494, MDE 0.4185, 32W/94L;
  tors - helix -0.2943, SE 0.0807, MDE 0.2260, 69W/57L.
  FAIL18: tors 5.569 / arm 6.032 / helix 5.887. other108: tors 3.471 / arm 2.745 / helix 3.761.
  Per-target min over {arm, tors, helix} 2.9632 (-0.2515 vs arm) is 92% accounted for by
  `ST.best_of_k_within`'s across-target null (k_eff 2.34): the arm-choice oracle is an order
  statistic, exactly as S24 L15-A found for the 19-arm panel. B3 is NOT a routing proposal.

## Hypothesis

**H_B3.** The pipeline's gain over sequence-only, d_i = arm_i - tors_i (negative = pipeline
better; 89 of 126 targets), is concentrated on a subset that is characterisable from native-free
quantities: length, composition, the retrieval pool's secondary-structure class (H/E/C string of
the top-75 windows' torsions through `I.ss_of`), the shipped distogram's per-pair bin entropy,
and ESM-2 contact-map statistics.

## Exact falsifier

A leave-fold-out classifier of sign(d_i) (ridge-logistic, alpha by inner CV within the 4
training folds; pinned folds `ST.pinned_folds`) must reach a held-out balanced accuracy above the
95th percentile of a 500-draw label-permutation null (labels permuted across targets, the same
nested fit each time), AND a leave-fold-out ridge regression of d_i must have held-out R^2 > 0
with the fold-clustered bootstrap CI of the held-out mean-squared-error reduction excluding zero.
If either fails, the set is NOT characterisable native-free at n = 126 and B3 records that with
its power: the permutation null's 95th-percentile accuracy is the MDE-equivalent for the
classifier; for the regression, the MDE on d is 0.369 A (above), so only a subset carrying a
mean gain above that can be resolved.

## Comparison arm, basis

The comparison is arm (built chain, `rmsd_arm`, 3.2148) vs tors (built chain from predicted
torsions, 3.7705) and vs helix (built chain, 4.0648): all three are built chains, so the basis is
matched on both sides. The selection-side arm (`shipped`, 3.4540) is carried as a third comparator
for the presenter's table only.

## Expected effect vs computed MDE

The record predicts a null: three routers over the same regime split were null in S12 (section I:
"no deployable signal tells the two regimes apart"), five router constructions failed for m*
(S22 L7, S23 L7), the finite-sample bound puts a p = 5 linear router at n = 281 for a 0.24 A
effect (S22 L10). One known structure: FAIL18 (14% of targets) carries 57% of the pipeline's
loss and is the set where sequence conditioning is harmful (S12: blind pipeline 5.425 vs shipped
6.019 on FAIL18); FAIL18 is defined by ORACLE zero recall and is not a native-free label. The
classifier's expected accuracy is at the null (0.50 balanced); the honest prior is that the
S13/S14 secondary-structure split (the helix arm is BETTER than the pipeline on FAIL18: 5.887 vs
6.032) is the only signal and it is weak.

## Memory, agent-hours

Reads persisted rows and 126 universe files (1.4 MB each, one at a time); features from
`s12/cache/disto_*.npz` and `s7/repr_cache/esmcon.npz`; < 0.5 GB (probe `p_probe_esm` 0.594 GB
covers a superset). 3 agent-hours.

## Rule 0 forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | sign(d) classifier and d regression, ridge, nested leave-fold-out | a tree ensemble (S22 L7: in-fold -0.276 -> held-out +0.079, the textbook overfit) |
| basis | built chain on both sides (rmsd_arm vs rmsd_build vs helix) | the S13 `base` rebuild (3.2126) as the incumbent |
| readout | per-target d as stored; no re-emission | rebuilding either arm |
| normalisation | features z-scored within training folds only | global z-scoring |
| null | label permutation across targets, 500 draws, nested fit repeated | a permutation of features within a fold |
| label | sign(d) and continuous d; FAIL18 membership reported as an ORACLE stratum, never a feature | routing gain (an order statistic, above) |
