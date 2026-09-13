# PROPOSAL C -- LEARN A BETTER DISTANCE PRIOR (lane P, Sprint 26)

Status: DRAFT 2026-09-13 09:35. C1 is final. C2 to C5 are filled in rung by rung from
`s26/results/p_ladder_report_<rung>_s0.json`, `s26/results/p_c4.json`, `s26/results/p_c5.json`,
and lane PH's `s26/C3_RESULT.md`. Verdict form per the campaign prompt: KEEP AS STATED / KEEP
WITH EDITS (and the edits) / REPLACE.

## What the proposal says

The pipeline's accuracy is bounded by its distance prior, and the prior can be improved by
learning: better inputs (C2), physics kept honest (C3), routers for the per-target set size and
scale (C4), and a learned correction of the common-mode error (C5). C1 asks that the closures
this rests on be reproduced first.

## What the repository already knows (the honest tuning of the proposal)

Selection inside the pool is closed at four levels: 58 distogram functionals and 29 in-band
signals (S17), 27 target-level calibration features (S17 L23), and the set-transformer with the
flat learning curve (S12). A perfect ranker inside the shipped top-25 returns 2.609 A (S17
L10/L23). The terminal operator consumes the set mean (S12). The pool's error is 68% common-mode
(S23 L9). Relaxation's apparent gain dissolves under a matched-magnitude random displacement
(S16 L27). The prior's derivative is the only steep lever (S24 L13) and only along the native's
direction (S25 L12). So the proposal is not "learn a better ranker"; it is "learn a better PRIOR",
and the open question is whether one is obtainable from inputs this machine can compute.

## C1 -- closure reproduction (final; `s26/PREREG_C1.md`, `s26/agentP_FINDINGS.md` section 7)

Both closures reproduce from the artefacts to the third decimal.

| closure | quantity | recorded | reproduced | artefact |
|---|---|---|---|---|
| set-transformer, real features | n = 8 / 16 / 32 / 64 / ~75 | 3.043 / 3.049 / 3.046 / 3.036 / 3.026 | 3.0433 / 3.0491 / 3.0456 / 3.0358 / 3.0258 | `s12/results/agg_dec_v2*.json` |
| set-transformer, leaked label | n = 8 / 32 / ~75 | 2.534 / 2.200 / 2.146 | 2.5342 / 2.2004 / 2.1460 | `s12/results/agg_dec_v2_oracle*.json` |
| perfect ranker in the top-25 | band best / random / dist argmin | 2.609 / 3.468 / 3.454 | 2.6087 / 3.4676 / 3.4540 | `s17/results/inband.json` |

The real-feature learning curve spans 0.023 A over an 8x change in training data (inside the
seed spread); the leaked label reaches 57% of its final gain at n = 8. The S7-11 ESM number
(-0.288 A on selection) is cited, artefact lost, and re-measured by C2's noesm/pca32 rungs.

## C2 -- the prior-input ladder (pending; `s26/PREREG_C2.md`, `s26/p_ladder.py`)

Rungs: shipped (anchor), noesm, conly, pca32 (reproduces the pinned model: fold-0 gate argmin
25/25, top-75 overlap 1.000), pca32f, pca128, raw, esm8m, wide, pairnet, mix. Falsifier: beats
the shipped posterior on the built chain by more than its own MDE, fold-clustered CI excluding
zero, 5/5 folds. Results table filled rung by rung (built chain primary; point cloud and
selection carried; gam_eff with cos beside it).

## C3 -- AMBER refinement kept honest (lane PH; `s26/PREREG_c3_control.md`, `s26/C3_RESULT.md`)

Filled from PH's result when it exists. Lane P's `PREREG_C3.md` is a different hypothesis
(entered as `IDEA_amber_prior_partner.md`).

## C4 -- routers on a feature set no previous router used (pending; `s26/PREREG_C4.md`)

## C5 -- predict and subtract the common mode (pending; `s26/PREREG_C5.md`, one agent-day)

## Verdict

Pending.
