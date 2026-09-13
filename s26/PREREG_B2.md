# PREREG_B2 -- A LEARNED-FOLDING-MODEL PRIOR AT FEASIBLE SCALE (lane P, Sprint 26)

Written 2026-09-13, before any endpoint number exists. B2 as originally implied by Proposal B
(run the folding model and use its output) is conditional on B1, and B1 is infeasible on this box
(`s26/results/b1_feasibility.json`). B2 therefore becomes the REPLACEMENT question: what part of
a large pretrained protein model CAN this box consume as a prior input, and does any of it move
the endpoint? It is executed inside `s26/p_ladder.py` as three rungs so that every arm is scored
through the identical path and the identical falsifier as C2.

## Hypothesis

**H_B2.** The structure-supervised part of ESM-2 650M (its contact head) and the size of the
language model are inputs the machine can compute, and at least one of the following beats the
shipped posterior on the built chain by more than its MDE, fold-clustered, 5/5 folds:
- `conly`: the shipped architecture on the physicochemical pair block PLUS the 13 contact-head
  columns only (no per-residue embedding block), 55-d -- the contact head's standalone value as a
  PRIOR INPUT, which S17 L23 measured only as a ranking feature;
- `esm8m`: esm2_t6_8M reps (PCA-32 per fold) and its own contact head in place of the 650M
  model's -- the downward point on the model-size axis;
- `pca32` (650M, shipped) is the middle point; the 3B point is infeasible (5.68 GB fp16 on disk,
  not downloaded; >= 5.7 GB resident against 4.4 GB of headroom).

## Exact falsifier

PREREG_C2's: built chain (`arm`), effect beyond its own MDE, fold-clustered CI excluding zero,
5/5 folds. In addition the size axis is read as a MONOTONE test: if arm(esm8m) is not worse than
arm(pca32) by more than the MDE, the ESM channel's value does not scale with model size at this
length and the "bigger language model" direction (the only thing B1 would have added) is closed
by its own downward rung; if it IS worse, the slope per log-parameter is reported with its CI and
extrapolated to 3B ONLY as an upper bound with the cosine caveat.

## Comparison arm, basis, readout, statistics

As PREREG_C2 sections 3-4 (same code, same tables, same MDEs: 0.09-0.20 A on the built chain,
0.28 A on selection). `conly` vs `noesm` isolates the contact head; `pca32` vs `conly` isolates the
embedding block; `esm8m` vs `pca32` isolates model size at fixed architecture.

## Expected effect

The record: ESM buys a FILTER, not a discriminator (S17 L23); the contact head's in-band
increment over its ESM-free twin is +0.050 [+0.005, +0.097] in Spearman and is mostly
compactness. S7-11: the physicochemical block ties one-hot (+0.001). So `conly` is expected
between `noesm` and `pca32` on selection, closer to `noesm`; the built chain is unknown. `esm8m`
is expected worse than `pca32` by an unknown amount; if the whole 0.29 A ESM gain survives at 8M
parameters, the gain is not "scale" but "any language model", which is what S7-11's "any
reduction of it" already hinted.

## Memory, agent-hours

`conly`: 55-d X, 0.12 GB; `esm8m`: one forward pass of the 30 MB model over 6,916 sequences
(minutes, < 0.5 GB) plus the lean trainer (< 1.0 GB). Reported from `s26/jobs_done/*.json` peak
RSS. 2 agent-hours on top of C2's pipeline.

## Rule 0 forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | shipped Bayes risk from a genuine `core.predict.Distogram` | any new functional |
| basis | built chain primary, cloud and selection carried | pLDDT-like confidence metrics |
| readout | shipped top-75 average, m = 75 | re-tuning m |
| normalisation | esm8m's PCA refit per fold on training sequences; its OWN contact head (not the 650M's) so the rung is one model, not a hybrid | mixing the 8M reps with the 650M contact head |
| null | shipped; noesm; pca32 as the middle size point | a randomly initialised 8M ESM (a different question: architecture vs pretraining) |
| label | built-chain CA-RMSD; gam_eff with cos beside it | MAE |
