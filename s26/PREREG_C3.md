# PREREG_C3 -- PHYSICS AS A PRIOR PARTNER, NOT AS A RANKER (drafted by lane P for the Physics lane, Sprint 26)

Written 2026-09-13 so the Physics lane can run it without re-deriving the discipline. **Lane P
did not have the text of Proposal C's item C3; this file states the physics-lane version of "learn
a better prior" that the record leaves formally open, and the Physics lane must confirm the
hypothesis matches the proposal before the first run.** Nothing here has run.

## What the record already knows

- Neither physics energy ranks nativeness on this pool: Legacy +0.330 and AMBER +0.455 WORSE than
  a random 75-subset, 5/5 folds (S25 L16); the functional lever is closed in all five forms
  (filter, partition, additive score, per-target sign, ORACLE weight worth 0.015 A: S24 L16,
  S25 L18); nested CV anti-weights AMBER (w = -0.5 on 4 of 5 folds, S24 L16).
- AMBER selects along a genuinely non-parallel direction (bias cosine 0.569 vs the incumbent,
  -0.215 against its floor, S24 L10-A) and cannot pay for it on quality (4.134 A standalone).
- 63,000 genuine ff14SB/GBn2 single points for the 126 x 500 pool are cached
  (`s24/cache_amber/*.npz`, tracked since S26 L1), so any pool-level AMBER question is answerable
  without OpenMM.
- Relaxation's apparent gain dissolved under a matched-magnitude random displacement (S16 L27).

## Hypothesis (two arms, cheapest first)

**H_C3a (cheap; no OpenMM).** An energy-weighted POOL HISTOGRAM -- per pair, the 17-bin
histogram of the K = 500 pool's distances weighted by exp(-beta * zrank(E_AMBER)) from the cache
-- used as the mixture partner of the shipped posterior (PREREG_C2's `mix` rung with H_pool
replaced by H_AMBER, lam and beta chosen leave-fold-out on the built chain) beats the shipped
posterior. This is the one place AMBER's non-parallel direction has NOT been tried: as a
DISTRIBUTION inside the prior rather than as a score, filter or partition.

**H_C3b (expensive; two OpenMM slots).** Retraining the shipped MLP on AMBER-RELAXED model-1
targets (restrained ff14SB/GBn2 relaxation of each training chain, k = 10, the production
relaxation) instead of raw model-1 coordinates gives a better prior. The record predicts a null
by arithmetic: label noise is 0.5 A per pair and de-noising is worth 0.044 A (S8-14); a
relaxation moves a chain by 0.16 A (production `amber_moved` on 1A13) and the relaxed target is
still one member of a ~1 A ensemble.

## Exact falsifier

PREREG_C2's: built chain, effect beyond its own MDE, fold-clustered CI excluding zero, 5/5 folds;
`lam = 0` and `beta = 0` are the incumbent bit-exactly and are asserted per target. For H_C3a the
per-target oracle over (lam, beta) is reported beside `ST.best_of_k_within` and its split-half
transfer. For H_C3b the pinned-model reproduction gate of PREREG_C2 (rung pca32) must pass first,
so that "relaxed labels" is compared against a retrained, not a pinned, control.

## Comparison arm, basis, readout, MDE, memory, agent-hours

As PREREG_C2 (shipped posterior through the same path; built chain primary; MDE 0.09-0.20 A).
H_C3a: reads `s24/cache_amber` and the universes; < 0.6 GB; 2 agent-hours. H_C3b: 6,790
restrained relaxations at ~10 s each = ~19 h of OpenMM at two slots (~10 h wall) -- the Physics
lane decides whether that is worth it AFTER H_C3a and the C2 ladder have reported; lane P's
recommendation is that it is not, given S8-14.

## Rule 0 forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | shipped Bayes risk from a genuine `core.predict.Distogram` built from the mixed posterior; rank-standardised AMBER (the only monotone normalisation, S24 L10) | moment-standardised AMBER (the 1e28-outlier artefact, S25 L16); AMBER as an additive score (closed) |
| basis | built chain primary; cloud and selection carried | AMBER energy of the emission |
| readout | shipped top-75 average, m = 75 | AMBER relaxation of the emission (S16 L27) |
| normalisation | beta on the zrank scale in {0.5, 1, 2, 4}; lam on PREREG_C2's grid; both leave-fold-out | per-target beta |
| null | lam = 0 (incumbent); the UNWEIGHTED pool histogram (PREREG_C2's `mix`) so AMBER's contribution is isolated from the histogram's; a rank-permuted AMBER (preserves the marginal, destroys the correspondence, S24 L16) | a random energy |
| label | built-chain CA-RMSD | native energy percentile (S7-10: AMBER puts the native at the 54th percentile) |

## ADDENDUM 1 (2026-09-13 00:45) -- the coordinator's definition of C3 differs from this draft

`s26/LEDGER.md` L5 (coordinator) defines lane PH's C3 as "the C3 matched-random control (stage 1
on the production cache, stage 2 on the best C2 rung when lane P delivers it)" attached to the
AMBER steric reject filter. This file was drafted by lane P without that brief and describes a
different, physics-as-prior-partner experiment (H_C3a/H_C3b). **Lane PH's brief governs.** The
draft above is offered to PH as an optional arm; it is not the C3 of record. Lane P will hand PH
the best C2 rung's per-target posteriors and endpoints (`s26/results/p_ladder_<rung>_s0.json`)
the moment they exist after sign-off, which is what stage 2 needs.

## ADDENDUM 2 (2026-09-13 08:45) -- NOT PROPOSAL C's C3, per the coordinator (LEDGER L16b)

Proposal C's item C3 is "AMBER refinement kept honest": re-run the restrained relaxation on the
best C2 rung's output against the S16 matched-magnitude random displacement control, owned by
lane PH (`s26/briefs/PH.md`, pre-registered in `s26/PREREG_c3_control.md`, result to be
`s26/C3_RESULT.md`). This file's H_C3a/H_C3b are a DIFFERENT hypothesis (physics as a prior
partner). Per L16b the file is neither deleted nor retitled; its cheap arm H_C3a is entered in
the tournament as `s26/IDEA_amber_prior_partner.md` and the Adversary decides whether it runs.
H_C3b (relaxed training labels) is withdrawn from the tournament: S8-14's arithmetic (label noise
0.5 A per pair, de-noising worth 0.044 A) already bounds it below every MDE in PREREG_C2.
