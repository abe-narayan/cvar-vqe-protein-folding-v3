# IDEA -- AMBER as a steric REJECT filter, not a ranker (lane PH; mandatory direction)

## Hypothesis
The genuine ff14SB/GBn2 single point cannot rank candidates (worse than a random 75-subset at
2.42x MDE, 5/5 folds, `s25/results/phys_suite.json`), but it does measure one thing reliably:
whether an unrelaxed ideal-geometry rebuild contains a hard steric overlap (58.6% of every pool
above 1e4 kcal/mol; ten of 500 candidates carry 99.7% of its variance,
`s25/results/phys_landscape.json`). A binary reject at a fixed physical threshold, with the count
set by the physics rather than by rank and the retained set refilled from the next-ranked
survivors, uses that one bit and nothing else. If clash-ridden members carry error that does not
cancel in the average, rejecting them lowers the built-chain RMSD.

## Why the record does not already close it
- S19 Q3 (`s19/agentC_FINDINGS.md` section 3): rejected a FIXED COUNT r in {2,5,10,19} of the
  top-75 by rank on Legacy steric, min-heavy and the AMBER single point; no refill; point cloud
  only. AMBER at r=10: +0.0109 [-0.0049, +0.0285] vs matched-random, NOT MEASURED. Steric and
  min-heavy worse than random with a sign. Closed the fixed-count form.
- S24 D1-C (`s24/LEDGER.md` L10-A): physics pre-filters as fixed-fraction partitions of the pool,
  n=30, all harmful in sign, none past MDE. Closed the fixed-fraction form at n=30.
- S25 (functional lever closed in five forms): filter (that S24 one), partition, fitted score,
  equal-weight score, per-target audit. None is a physical-threshold reject with refill.
- Not measured anywhere: the count chosen by a threshold (zero on clean targets), refill to
  m=75, and the built-chain basis.

## Exact falsifier
On the built chain at the primary threshold 1e4 kcal/mol, `R_1e4` (reject + refill) or `S_1e4`
(reject, shrink) beats the shipped top-75 by more than its own MDE with the fold-clustered CI
excluding zero and 5/5 folds, AND beats its matched control (same count rejected at random, with
or without refill respectively) by more than that comparison's own MDE with the fold CI excluding
zero. Otherwise closed on this instrument. Full pre-registration: `s26/PREREG_amber_reject.md`.

## Expected effect against the MDE
Prior: null or harmful (S19 Q3, S24 D1-C, S23 L5 "removing deviant members costs +0.142 A",
S25 "AMBER's top-75 is +1.10 A more expanded"). The MDE on the built chain at n=126 will be
0.03 to 0.11 A depending on how many targets the threshold touches; the effect I expect is
within +-0.03 A. If the falsifier clears, it changes the closed table: "AMBER has no accuracy
role" becomes "AMBER has a reject role", and the deployed filter stage gains one step.

## Cost
Census (native-free, before the gate): 3 min CPU, < 0.5 GB. Point cloud: 5 min. Built chain:
up to 3.2 h CPU (126 targets x <= 21 projections x 4.4 s), < 0.7 GB, checkpointed per target.
Agent-hours: 3. No AMBER compute (the 63,000 cached single points are re-used).

## Information value
Either way the answer closes the last untested form of "physics as a filter" and gives the
report a sentence with a number: how many of the shipped top-75 are physically impossible, and
what removing them is worth on the production basis.
