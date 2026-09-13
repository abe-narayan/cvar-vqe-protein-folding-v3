# PREREG_C5 -- PREDICT THE COMMON-MODE DIRECTION ON HELD-OUT FOLDS AND SUBTRACT IT (lane P, Sprint 26)

Written 2026-09-13, before any number exists. Budget: one agent-day, hard stop.

## What the record already knows

- The pool's error is 68% common-mode (S23 L9, exact identity in the medoid frame,
  `s23/errdecomp.py`): mean_k|e_k|^2 = |ebar|^2 + mean_k|d_k|^2, with |ebar|^2 invisible from
  inside the pool. Averaging removes only the idiosyncratic 32%.
- The harmful coherent component is SHARED across predictor families and two thirds of it is
  reproduced by zero-information references (S19 L11); the mode that owns 52.5% of the harm is
  the per-target separation profile (S19 L14 section 1); the best native-free estimator of it is
  the retrieval pool, which shares 0.87 of it -- "the estimator is made of the bias" (S19 L14
  section 2). Every native-free error-direction steering arm in S16 chose "do nothing".
- A residual conditioned on target-level information makes coherent mistakes by construction,
  and at matched accuracy coherent mistakes cost 1.1-1.4 A more than i.i.d. ones (S24 L7); the
  per-target scale s* (the radial component of ebar) is unreachable in principle (S23 L6d/L9).

So C5 is the direct, cheap, fully controlled test of a direction the record has closed three
times from three sides. It is run because the lane brief asks for it and because the ORACLE
ceiling and the matched-random control make the result interpretable whichever way it goes.

## Construction

Per target: the shipped top-75 members in the medoid frame (`s23/errdecomp.py`'s frame), the
cloud c, the native t (ORACLE, training folds only). The common mode is ebar = c - t after ONE
rigid Kabsch fit of c onto t (rotation only), an (n, 3) vector per target.

A learnable target needs a native-free frame. Two representations, both pre-declared:
- (R1) DISTANCE SPACE, frame-free: e_D = D(c) - D(t) over the score's pair set, decomposed on the
  S19 L14 mode basis (per-separation-shell means over the 5 `SHELLS`, a per-residue additive
  term, offset, stretch). The correction is applied by moving c toward D(c) - e_hat_D with 200
  steps of the weighted stress majorisation already in `core.predict.Distogram.realize`
  (chain restraint on), INITIALISED AT c so the mirror image is never crossed, then `I.project`.
- (R2) COORDINATE SPACE in the cloud's own principal-axis frame (PCA of c; PC1 sign fixed by the
  N->C direction, PC3 by right-handedness, PC2 by the sign of the cloud's chirality signature):
  ebar expressed in that frame, per residue, padded to n = 16.

Arms (every arm through the SAME readout and projection):
- GLOBAL: the mean training-fold common mode (per length band 9-11/12-13/14-16 in R2; per shell
  in R1), scaled by alpha chosen on the training folds; applied to the held-out fold.
- RIDGE: per-target prediction of the mode coefficients from native-free features (length,
  composition, ESM pca32 mean, the pool's per-shell mean distances and their spread, the
  distogram's expected per-shell profile and entropy), ridge with alpha by inner CV.
- ORACLE ceiling: subtract the true ebar (R2) / true e_D (R1): the idiosyncratic-only endpoint.
- RANDOM-MATCHED: a random direction of the SAME magnitude as the predicted correction, same
  operator (S16's discipline; `control-must-match-the-operators-space`).
- SCALE-ONLY: the S23 s* as the one-parameter version, ORACLE, for scale.

## Hypothesis and exact falsifier

**H_C5.** The predicted component of the common mode, fitted on training folds, reduces the
held-out built-chain CA-RMSD.

**Falsifier.** Held-out GLOBAL or RIDGE minus incumbent must be beyond its MDE with the
fold-clustered CI excluding zero and 5/5 folds the same sign, AND must beat RANDOM-MATCHED by the
same standard. If not, C5 is refuted at one agent-day and the S16/S19/S23 closure stands with a
fourth instrument. If the ORACLE ceiling itself is small on the built chain (the projection may
absorb it), that is reported first: it bounds what any predictor could have bought.

## Comparison arm, basis, readout

Incumbent: the shipped cloud through `I.project` (3.2148), paired per target; point cloud carried
(3.0483). Readout: the shipped top-75 uniform average; only the post-average correction differs.

## Expected effect vs MDE

MDE for a point-cloud change of this size: 0.05-0.09 A (PREREG_C2 section 4). Expected: GLOBAL
~0 (it is the length-typical contraction the projection already partly repairs; S24 L14: averaging
contracts the cloud, and the scale half is bimodal, S23 L2, so a global correction cancels);
RIDGE ~0 or harmful (S22 L10 bound; S24 L7 coherence). ORACLE ceiling: large on the point cloud
(the common mode is 68% of the squared error), smaller on the built chain.

## Memory, agent-hours

126 universes one at a time, the projection at 2.9 s per target per arm (5 arms, ~30 min);
< 0.6 GB. One agent-day, including the write-up; if R2's frame convention fails its own
round-trip test (apply the oracle ebar in the canonical frame and recover the native to < 1e-6)
R2 is dropped and R1 alone is reported.

## Rule 0 forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | a post-average additive correction in a native-free frame; alpha and the ridge fitted on training folds only | correcting the members before averaging (that is the residual lane, closed S24 L7) |
| basis | built chain primary; point cloud carried | distance-space MAE of the corrected matrix |
| readout | shipped top-75 average, then correction, then projection | correcting after projection (leaves the ideal-geometry manifold) |
| normalisation | R1 modes on the S19 basis; R2 in the PCA frame with declared sign rules; both reported | a single representation chosen after seeing which works |
| null | incumbent; RANDOM-MATCHED at the predicted magnitude; the ORACLE ceiling as the bound | the constant-helix source (not an operator in this space) |
| label | built-chain CA-RMSD | the fraction of |ebar|^2 removed |
