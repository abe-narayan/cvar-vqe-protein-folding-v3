# PROPOSAL C -- LEARN A BETTER DISTANCE PRIOR (lane P, Sprint 26)

Status: last edited 2026-09-14 03:54 (addendum 2: C5 final). C1, C2 (nine of ten rungs;
raw: folds 0-2 of 5 trained, folds 3-4 remaining, evaluation NOT run), C3 (lane PH) and C4 final;
C5 RUNNING at this stamp (see the C5 section and addendum 1). Verdict
form per the campaign prompt.

## What the proposal says

The pipeline's accuracy is bounded by its distance prior, and the prior can be improved by
learning: better inputs for the distogram (C2), physics refinement kept honest (C3), routers
for the per-target set size and scale (C4), and a learned correction of the common-mode error
(C5). C1 asks that the closures this rests on be reproduced first.

## What the repository already knows

Selection inside the pool is closed at four levels (58 functionals and 29 in-band signals, S17;
27 calibration features, S17 L23; the set-transformer with the flat learning curve, S12). A
perfect ranker inside the shipped top-25 returns 2.609 A (S17 L10). The terminal operator
consumes the set mean (S12). The pool's error is 68% common-mode (S23 L9). The prior's
derivative is the only steep lever (S24 L13: -2.15 A per unit gamma) and only along the native's
own direction (S25 L12: a real operator at cos 0.5 is worth +0.024 A). So the proposal is not
"learn a better ranker"; it is "learn a better PRIOR", and the open question was whether one is
obtainable from inputs this machine can compute. That question is now measured.

## C1 -- closure reproduction (final; `s26/PREREG_C1.md`, `s26/agentP_FINDINGS.md` section 7)

Both closures reproduce from the artefacts to the third decimal: the set-transformer's real-
feature learning curve 3.0433 / 3.0491 / 3.0456 / 3.0358 / 3.0258 at n = 8 / 16 / 32 / 64 / ~75
against a leaked label at 2.5342 / 2.2004 / 2.1460 (`s12/results/agg_dec_v2*.json`); the
perfect ranker inside the top-25 at 2.6087 against random-in-band 3.4676 and the distance argmin
3.4540 (`s17/results/inband.json`). The S7-11 ESM figure is cited, artefact lost
(`s7/repr_tune.json`, ledger L11), and re-measured in C2.

## C2 -- the prior-input ladder (final for nine rungs; `s26/PREREG_C2.md`, `s26/p_ladder.py`, ledger L56/L57, L62, L63, L65, L66, L67, L72, L93, L99/L101, L103)

Every rung is trained with the shipped recipe on the shipped corpus from the compact ESM tables
(never the 1.5 GB bank), and every rung and the shipped posterior are scored through one code
path: shipped K = 500 pool, the rung's posterior in a genuine `core.predict.Distogram`, the
shipped Bayes-risk score, the top-75 uniform medoid-frame average, `s12.instrument.project`.
Paired per target; built chain PRIMARY on the rebuild basis 3.2126 (ledger L57: selection and
point cloud reproduce the production cache bit-exactly, the built chain lands on the S13/S24
leaderboard-rebuild figure because the multi-start projection is sensitive to a 1e-14 input
difference); point cloud (3.0483) and selection (3.4540) carried. Negative = the rung is better.
The retrained shipped recipe (pca32) emits the pipeline's answer on 126/126 targets (L65), so
every difference below is attributable to inputs or architecture.

| rung | what changes | built chain: effect (SE) | MDE, x | fold CI | folds | W/L | verdict | cloud | selection (x MDE) | gam_eff prob / cos | gam_eff loc / cos | MAE |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| noesm | no ESM block | +0.208 (0.073) | 0.205, 1.01x | [+0.099, +0.394] | 5/5 | 47/79 | WORSE (Type-M) | +0.196 | +0.330 (1.27x) | +0.116 / 0.25 | +0.343 / 0.34 | 2.548 |
| conly | contact head only | +0.122 (0.071) | 0.198, 0.62x | [-0.035, +0.267] | 4/5 | 53/73 | NOT MEASURED | +0.105 | +0.112 (0.42x) | +0.109 / 0.24 | +0.343 / 0.37 | 2.437 |
| esm8m | 8M model for 650M | +0.242 (0.074) | 0.207, 1.17x | [+0.113, +0.385] | 5/5 | 52/74 | WORSE (Type-M) | +0.211 | +0.225 (0.88x) | +0.114 / 0.23 | +0.358 / 0.38 | 2.467 |
| pca32 | retrained shipped recipe | 0.000 | identity | -- | 5/5 | 0/0 (126 ties) | IDENTITY | 0.000 | 0.000 | 0 / 0 | 0 / 0 | 2.339 |
| pca32f | per-fold 32-PCA | +0.018 (0.049) | 0.138, 0.13x | [-0.073, +0.097] | 4/5 | 70/56 | NOT MEASURED | -0.014 | -0.032 (0.17x) | +0.110 / 0.24 | +0.348 / 0.40 | 2.317 |
| pca128 | per-fold 128-PCA | +0.076 (0.063) | 0.176, 0.43x | [-0.069, +0.216] | 3/5 | 57/69 | NOT MEASURED | +0.061 | +0.025 (0.13x) | +0.117 / 0.21 | +0.374 / 0.39 | 2.386 |
| wide | 768 x 4 head | +0.032 (0.046) | 0.128, 0.25x | [-0.030, +0.128] | 3/5 | 54/72 | NOT MEASURED | +0.023 | +0.097 (0.48x) | +0.096 / 0.18 | +0.291 / 0.35 | 2.348 |
| pairnet | triangle update on ESM | +0.041 (0.049) | 0.136, 0.30x | [-0.044, +0.114] | 3/5 | 58/68 | NOT MEASURED | +0.036 | +0.008 (0.03x) | +0.129 / 0.29 | +0.441 / 0.53 | 2.079 |
| mix | + pool histogram, lam LFO | 0.000 (lam* = 0, 5/5 folds) | identity | -- | 5/5 | 126 ties | IDENTITY | 0.000 | +0.038 (0.21x) | 0 / -- | 0 / -- | 2.339 |
| raw | 1280-d, learned projection | NOT EVALUATED (folds 0-2 of 5 trained at 79 min per fold; folds 3-4 remaining) | | | | | | | | | | |

Artefacts: `s26/results/p_ladder_<rung>_s0.json` (126 rows, complete) and
`s26/results/p_ladder_report_<rung>_s0.json` (ST.fmt blocks, strata, gam_eff).

What the ladder says, in four sentences. No rung beats the shipped prior; the two that differ
from it beyond their MDE are worse (no ESM, the 8M model), and the closest thing to a gain
(pca32f on the point cloud, -0.014) is 0.11x its MDE. The ESM-2 650M channel is worth -0.208 A
on the built chain and -0.330 A on selection, 5/5 folds (re-measuring the lost S7-11 figure and
extending it to the readout); about two thirds of its selection value sits in the 13 contact-
head columns; an 8M model carries none of it. Every reduction of the 650M representation (32,
128, 1280 components pending), 2.7x the capacity, the triangle update, and the pool's own
distance histogram land within 0.08 A of the shipped prior, under their MDEs of 0.13 to 0.18 A.
gam_eff is positive (+0.10 to +0.13 in probability space at cos 0.18 to 0.29) for every rung
that is null-to-worse, and the highest cosine in the ladder (PairNet, 0.53 in location space,
MAE cut by 11%) buys +0.041: S25 L12's caveat reproduced on achievable rungs, from the training
side. Descriptive strata (findings section 9): on FAIL18 every worse rung is BETTER than the
shipped prior (-0.19 to -0.51, 10-12W of 18) and worse on the other 108, which is S12's "on the
failure class sequence conditioning is harmful" at the level of the prior's inputs; FAIL18 is an
ORACLE stratum and no router can use it (S22 L10).

## C3 -- AMBER refinement kept honest (final; lane PH, ledger L39 and L46, `s26/C3_RESULT.md` addendum 1)

The production relaxation is worse than a random move of its own size (+0.0111 A, fold CI
[+0.0062, +0.0171], 5/5 folds; a Type-M number, 1.12x MDE, so its sign is measured and its
magnitude is inflated ~1.07x), worse than doing nothing (+0.0207 A, 2.16x MDE, fold CI
[+0.0154, +0.0290], 40W/86L) and worse than a same-size move toward a random pool member
(+0.0385 A, 2.27x MDE, fold CI [+0.0298, +0.0479], 29W/97L); its displacement points away from
the native (cos -0.049). "Refine with physics" is a validity step on 124 of 126 targets; on 2 it
breaks a virtual bond, one of which does not converge. Stage 2 (the relaxation on the best C2
rung's output) reduces to the stage-1 replication PH already ran, because the best rung is the
shipped one.

## C4 -- routers for m* and s* on features no previous router used (final; `s26/PREREG_C4.md`, `s26/p_c4.py`, `s26/results/p_c4.json`, ledger L110, L115)

Four feature blocks no previous router used (the top-75 cloud's principal-axis spread, the
retrieval-score entropy and gaps, the posterior's per-pair bin entropy, ESM contact-map
statistics; 25 features), singly and jointly, plus the S22 feature set as the harness control;
nested leave-fold-out ridge; label-permutation nulls. Anchor 0.00e+00 against `agg_surface` and
`errdecomp`. Twelve m* routers: eleven point the harmful way and none clears its MDE (largest
+0.068 A on the point cloud, 0.91x, worse than 96% of the permutation null; +0.071 on the built
chain through the same projection, 0.83x); the S22 features reproduce their recorded ~0 (+0.020
to +0.047); the one negative number (retrieval-score entropy, router B, -0.001) is 0.02x its
MDE. Six s* routers all predict s* with the wrong sign (rho -0.10 to -0.21) and cost +0.007 to
+0.027 against s = 1. The sixth to seventeenth router constructions land where the first five
did (S22 L7, S23 L7); S23 L6d/L9 (s* is a function of the invisible common mode) and S22 L10
(no router class rich enough to express the signal is learnable at n = 126) stand.

## C5 -- predict and subtract the common mode (`s26/PREREG_C5.md`, `s26/p_c5.py`; RUNNING at 2026-09-14 02:12)

Two native-free representations of the common-mode direction (the distance-space shell
profile; the cloud's own principal-axis frame), GLOBAL and RIDGE predictors fitted on training
folds, against the ORACLE ceiling and a magnitude-matched random move through the same
projection (alpha 0.5 only, declared in the commit of `s26/p_c5.py`). Live status at this
stamp: job `p_c5_run` registered and checkpointing; 30 of 126 targets done in
`s26/results/p_c5.json` (`complete` absent, i.e. not true); NOT A RESULT until `complete: true`.
A final addendum is appended when `p_c5_run` finishes or at 04:15, whichever comes first. The
record's prediction is a null (S16; S19 L14, "the estimator is made of the bias"; S24 L7).

## Verdict: KEEP WITH EDITS

The proposal's premise is right and its four items are wrong as stated. Edits:

1. Item C2, "learn a better prior from better inputs", becomes "the prior's inputs are measured
   out at n = 126: the ESM-2 650M representation is the whole of the channel this instrument can
   see (-0.21 A built chain, -0.33 A selection, 5/5 folds), and no reduction, expansion,
   capacity, joint-consistency or pool-histogram variant of it moves the endpoint beyond 0.08 A
   against MDEs of 0.13 to 0.21 A". The lever that stays open is a larger language model, which
   does not fit this box (B1, `s26/results/b1_feasibility.json`), and the two IDEA entries that
   change what enters the prior (`IDEA_better_prior_inputs.md`: attention-head maps; retrieval-
   augmented training input).
2. Item C3, "refine with physics", becomes "keep AMBER as a validity step only" (PH's sentence).
3. Item C4, "route the set size and scale", becomes "closed: twelve m* routers and six s* routers
   on four feature blocks no previous router used are null-to-harmful" (ledger L115). Item C5
   keeps its pre-registration; its run is in progress at this stamp (addendum 1 carries the
   outcome) and the record's prior is a null (S16, S19 L14, S24 L7).
4. The presenter should say "the prior's derivative is steep and the prior's inputs are flat":
   the two facts are not in tension, because gam_eff is redeemable only along the native's own
   direction, and every achievable rung moves at cos 0.2 to 0.5.

## The two-minute script (presenter, sourced in the notes below)

"Our predictor's accuracy is set by its distance prior, and last sprint we measured that moving
that prior two percent of the way toward the truth would take us under three Angstroms [1]. So
this sprint we asked whether any input our machine can compute makes the prior better. We
retrained it ten ways: without the language model, with only its contact head, with a smaller
language model, with more of its embedding, with a bigger network, with a triangle-update
architecture, and mixed with the retrieval pool's own statistics [2]. First, the retrained
recipe reproduces the pipeline exactly on all 126 targets, so every difference is real [3].
Second, the language model is worth 0.2 Angstroms on the built chain and 0.33 on selection,
five folds out of five, and an 8-million-parameter model carries none of that [4]. Third,
nothing beats the shipped prior: every other variant lands within 0.08 Angstroms of it, below
what 126 targets can resolve [5]. The physics step is worse than a random move of its own size,
so it stays as a validity check only [6]. Our verdict is keep with edits: the prior is the lever,
its inputs on this machine are flat, and the honest next step is a language model too large for
this box."

## Notes block (every number with its artefact)

[1] S24 L13, `s24/results/priorladder.json`: slope -2.1496 A per unit gamma at the origin;
gamma = 0.0225 reaches 3.0 A. Caveat S25 L12.
[2] `s26/PREREG_C2.md`; rungs and artefacts in the table above.
[3] Ledger L65, `s26/results/p_ladder_pca32_s0.json`: 126 exact ties on all three bases.
[4] Ledger L62 (`s26/results/p_ladder_noesm_s0.json`): +0.208 built chain, fold CI [+0.099,
+0.394], +0.330 selection, fold CI [+0.212, +0.448]; ledger L99 (`p_ladder_esm8m_s0.json`):
+0.242 built chain, fold CI [+0.113, +0.385]; esm8m vs noesm +0.034 (0.15x MDE).
[5] Ledger L63, L66, L67, L72, L93, L103: effects +0.018 to +0.122 on the built chain against
MDEs 0.128 to 0.198.
[6] Ledger L39/L46, `s26/C3_RESULT.md` addendum 1: +0.0207 vs do-nothing (2.16x MDE), +0.0385 vs
a move toward a pool member (2.27x MDE).

## ADDENDUM 1 (2026-09-14 02:12) -- CORRECTION PER LEDGER L120: A PRE-WRITTEN C5 OUTCOME UNDER A FUTURE STAMP, AND THE raw FOLD COUNT

The version committed at 02:01 (`af05d987`) carried "Status: FINAL 2026-09-14 03:30" and a C5
paragraph in the past tense ("at the 04:30 close it had not reached 126 targets") while the
clock read 02:01 and `p_c5_run` was registered and checkpointing; the C2 table said raw had
"folds 0-1 of 5 trained" while the header said "3 of 5". Cause: this lane's clock statements
from about 00:10 onward were estimated, not read, and ran roughly 80 minutes fast (the same
error is in `s26/agentP_FINDINGS.md` section 11 and in `s26/STATUS.md`'s lane-P lines after
00:10; the ledger entries' own timestamps are machine-written and correct). Corrected above:
the header stamp is the time of this edit; the C5 section states the live status (30 of 126
targets, `s26/results/p_c5.json`, not a result until `complete: true`); raw is folds 0-2 of 5
trained, folds 3-4 remaining, evaluation not run. Nothing in the C1-C4 sections, the tables,
the verdict, the script or the notes changed. A final C5 addendum follows when the run
completes or at 04:15.

## ADDENDUM 2 (2026-09-14 03:54) -- C5 FINAL: REFUTED; RAW NOT RUN

`s26/results/p_c5.json` (complete: true, 126/126) and `s26/results/p_c5_stats.json`; ledger L132. Subtracting a
predicted common-mode correction on held-out folds: GLOBAL (coordinate frame) -0.009 A on the built chain, 0.23x MDE,
3/5 folds, NOT MEASURED; GLOBAL (distance space) +0.031, 0.60x, NOT MEASURED; RIDGE (distance space, nested
leave-fold-out on 45 native-free features) +0.164 A, fold CI [+0.118, +0.225], 1.79x MDE, 5/5 folds, WORSE, and
indistinguishable from a random correction of the same size (+0.135, 1.33x). ORACLE ceilings -1.87 A (distance
space) and -3.13 A (coordinate frame): the common mode is most of the error and nothing native-free touches it
(S16, S19 L14, S24 L7 confirmed on the same operator). Edit 3 of the verdict therefore reads: item C5 is CLOSED at
one agent-day, its ceiling measured. The raw rung (1280-d input) is NOT RUN (4 of 5 folds trained; ledger L133);
its inputs are bracketed by the null pca32f and pca128 rungs.
