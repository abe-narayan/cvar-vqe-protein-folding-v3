# PREREG_coherence -- RETRAIN THE PRIOR AGAINST ITS NAMED HARMFUL MODE (lane P, Sprint 26; from IDEA_coherence_penalised_training.md, ordered by the coordinator at L77)

Written 2026-09-13 22:40 before any number exists. Code: `s26/p_coh.py` (imports `s26/p_ladder.py`
unchanged; the training jobs in flight use that module and it is not edited).

## Hypothesis

Two training-loss changes to the shipped MLP, everything else identical (corpus, order, 183-d
inputs from the compact tables, 40 epochs, batch 4096, AdamW 2e-3, OneCycle, soft-binned CE,
seed 0):

- **A5b (`sw_util`)**: per-pair loss weight favouring SHORT separations, the mirror of the
  anti-utility `lin` weight that S19 L11 measured at +0.181 [+0.081, +0.285]: w(sep) = max(6 -
  sep/3, 1) normalised to mean 1 (so sep 2 carries 5.3x the weight of sep 15). S19's own
  pre-registration (`s19/PREREG_A.md` A5) was blocked on RAM; the compact tables unblock it.
- **COH (`coh_l<lambda>`)**: minibatches grouped by chain; loss = the shipped per-pair CE + lambda
  x mean over chains in the batch of sum over the 5 `SHELLS` of (mean signed deviation of the
  differentiable posterior mean E[d] from the true distance, within the shell)^2 / (A^2, A = 1).
  This penalises exactly the per-target separation-profile mode that owns 52.5% of the harm
  (S19 L14) and nothing else. lambda in {0.3, 1.0, 3.0}; lambda = 0 is the shipped trainer
  (rung pca32, which reproduces the pipeline on 126/126, L65).

## Exact falsifier

PREREG_C2's, per arm: built chain (rebuild basis 3.2126), effect beyond its own MDE, fold-
clustered CI excluding zero, 5/5 folds. A5b is REFUTED if its CI spans zero or the sign is
positive (S19 A5's wording). COH is additionally scored on the S19 L14 mode: the held-out
per-target shell-profile residual of the posterior mean, |mean_shell(E[d] - D_nat)| averaged
over shells, against the shipped posterior's (an ORACLE diagnostic) so a null can be told
apart from a penalty that did nothing. lambda is chosen leave-fold-out on `arm` only if a
lambda clears the falsifier; otherwise every lambda is reported and none is selected.

## Comparison arm, basis, readout, statistics

As PREREG_C2 section 3: paired against the shipped posterior through the same path; built chain
primary, point cloud and selection carried; ST.fmt verbatim; gam_eff with cos beside it; the
length and FAIL18 strata as descriptive tables.

## Expected effect vs computed MDE

MDE on the built chain 0.13-0.21 A on this instrument (L62-L72: 0.128 to 0.205). Expected: A5b
0.00 +/- 0.05 (S19: the best single MLP was the UNIFORM one; a mirror of a harm is not a
guaranteed gain); COH -0.05 to +0.05 (two thirds of the harmful component is present in a
zero-information reference, S19 L11, so the penalty can act on at most a third of it). Both
likely "not measured"; the mode-residual diagnostic is the informative output either way.

## Memory, agent-hours

The materialised 183-d trainer (pca32's path: 453 s per fold, peak 1.248 GB in the probe; 0.69
GB in the chain run) plus chain grouping; est-ram 1.25 GB per the coordinator; only while the
governor snapshot shows 2.5 GB free; resumable per fold. 4 arms x 5 folds ~ 2.5 h; evals 10 min
each. 3 agent-hours.

## Rule 0 forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | shipped Bayes risk from a genuine `core.predict.Distogram` | any new functional |
| basis | built chain primary; cloud and selection carried | MAE, NLL, the mode residual as the decision |
| readout | shipped top-75 average, m = 75 | re-tuning m |
| normalisation | weights normalised to mean 1 so the effective learning rate is unchanged; lambda on a 3-point log grid declared here | tuning lambda on the endpoint per fold without declaring the grid |
| null | the shipped trainer (pca32 = the pipeline); `sw_lin` from S19 as the anti-utility direction (through the deployed fit, not re-run) | a random per-pair weight |
| label | built-chain CA-RMSD; the shell-profile residual as an ORACLE diagnostic | W/L |
