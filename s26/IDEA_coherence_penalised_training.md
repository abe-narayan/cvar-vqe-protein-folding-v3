# IDEA_coherence_penalised_training -- RETRAIN THE PRIOR AGAINST ITS NAMED HARMFUL MODE, AND RUN THE ONE S19 ARM THAT WAS BLOCKED ON RAM (lane P, Sprint 26; lane P's own entry)

## Hypothesis

The distogram's harm is not its magnitude but its COHERENCE: cross-pair sign correlation of the
residual (S19 L3), of which the per-target separation-profile mode owns 52.5% (ORACLE removal
3.610 -> 3.085, S19 L14 section 1), while the per-residue mode explains the most variance and
owns none of the harm. The shipped loss is a per-pair soft-binned cross-entropy that is blind to
coherence across the pairs of one chain. Two training-loss changes, both native-free at
inference and both scored through PREREG_C2's path:

1. **S19 A5b, unblocked.** Retrain with the loss weighted toward SHORT separations (the utility
   direction; `separation_weights` mirror of `lin`, i.e. `max(6 - sep/3, 1)`-shaped, pre-declared
   in `s19/PREREG_A.md`), against the uniform arm retrained under the identical script and the
   existing anti-utility arm `sw_lin` (+0.181 [+0.081, +0.285] through the deployed fit, S19
   L11). A5b was "NOT RUN, blocked on RAM: 6,429 of 6,788 training sequences absent from the hot
   ESM cache" (S19 L11); the compact tables (`s26/results/p_probe_esm.json`: 6,790/6,790 covered,
   0.6 GB) remove that block.
2. **Coherence penalty.** Minibatches grouped by chain; loss = the shipped per-pair CE plus
   lambda x mean over chains of sum over the 5 `SHELLS` of (mean signed deviation of E[d] within
   the shell)^2, where E[d] is the differentiable posterior mean. This penalises exactly the
   mode that owns the harm and nothing else; lambda on a 3-point grid chosen leave-fold-out on
   the built chain (lambda = 0 is the shipped trainer).

## Why the record does not already close it

S7-8 tried ranking losses (worse) and an anti-shrinkage loss (null): losses aimed at candidate
ORDER or at amplitude, not at cross-pair coherence. S19 L11 tried the long-favouring weight
(harmful at unchanged MAE) and left the short-favouring mirror pre-registered and unrun. No arm
in the record has put a coherence term in the training loss; S19 L14's wall ("the estimator is
made of the bias") is about estimating the mode native-free AT INFERENCE, and a training loss
sees the native. The counter-evidence to weigh: two thirds of the harmful component is present in
a zero-information constant-helix reference (S19 L11), so the penalty can at most act on the
remaining third, and S9-8's capacity curve says learned objectives on this instrument prefer the
simplest form.

## Exact falsifier

PREREG_C2's, per arm: built chain, beyond its own MDE, fold-clustered CI excluding zero, 5/5
folds; the A5 controls (uniform retrained under the same script; `sw_lin` as the anti-utility
direction) and MAE reported and not decided on. A5b is REFUTED if its CI spans zero or its sign is
positive, exactly as `s19/PREREG_A.md` states. The coherence arm is additionally scored on the
S19 L14 mode decomposition (does the separation-profile residual actually shrink on held-out
targets, an ORACLE diagnostic) so a null can be told apart from a penalty that did nothing.

## Expected effect vs computed MDE

MDE 0.09-0.20 A on the built chain (PREREG_C2 section 4; the S19 zoo's sw_lin - sw_none paired
sd 0.600 gives 0.150 on the deployed fit). Expected: A5b -0.05 to +0.05 (the mirror of a +0.181
harm is not a guaranteed gain; the loss weight changes where capacity is spent, and S19 says the
best single MLP was the UNIFORM one); coherence arm -0.10 to +0.05 with the ORACLE bound on the
mode at -0.525 on the deployed-fit basis. Both likely "not measured"; the mode-decomposition
diagnostic is the informative output either way.

## Memory and agent-hours

Lean trainer with chain-grouped batches, 183-d inputs: ~1.3 GB (the pca32 probe: 1.248 GB peak,
453 s per fold); 2 arms x 3 lambda x 5 folds ~ 6 h machine time at worst, 2.5 h for A5b alone;
eval 7 min per arm. 4 agent-hours. Code: `s26/p_ladder.py`'s `LeanTrainer` with a `chain` index
and a `sample_w` / `penalty` hook (not yet written; written only if this entry wins the tournament).
