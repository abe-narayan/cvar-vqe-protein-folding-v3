# IDEA_tiebreak_noise_floor -- THE PIPELINE'S OWN NOISE FLOOR: THE K = 500 POOL BOUNDARY IS SET BY BLOSUM TIES BROKEN BY CORPUS ORDER (lane W, Sprint 26; own idea 1)

## Hypothesis

The retrieval pool is the top 500 windows by BLOSUM62 sum, an integer, with ties broken by the
stable argsort over the universe's corpus order (which is why the `pdbs/` file set is pinned:
"sort order feeds BLOSUM tie-breaking", state brief section 3). Measured native-free this turn
(`s26/results/w_selfcopy_census.json`, job `w_selfcopy_census`, 20 s, 0.064 GB): EVERY one of
the 126 targets has a tie at the boundary; the median tie class at the boundary score holds
115 windows, of which about 57 are inside the pool and 60 outside. So on every target roughly
11% of the pool is chosen by an arbitrary convention, not by the score. The hypothesis is that
re-drawing the boundary members at random (a uniformly random subset of the tie class of the
same size, everything downstream unchanged: shipped score, top-75, average, projection) moves
the built chain by a per-target sd s_tie and moves the 126-target mean by m_tie, and that m_tie
and the MDE of a paired contrast between two random tie-breaks are the pipeline's own noise
floor: no paired effect below that floor is interpretable, whatever its CI says, because the
same effect size is produced by an unstated convention. Several effects the record carries at
the hundredths level would then have to be quoted beside it: the restrained AMBER relaxation
-0.022 [-0.036, -0.009] (S25 / memory `averaging-space-beats-the-objective`), the ORACLE
functional weight 0.015 (S24 L16), the all-atom reranking +0.004, C27's +0.0004, and the 1e-2
differences between projection gradient modes (L19).

## Why the record does not already close it

- S17 (`s17/LEDGER.md` L628) measured the pool's ORACLE BEST member under random tie-breaks:
  1.708 A against 1.711, sd 0.018, and called the 11.4% of the pool set by corpus order
  "robust". That is the oracle best of 500, not the emitted structure.
- S8-12 / S10 (`docs/FINDINGS.md:3668, 6106`) note that "a different tie-break gives a
  different pool" (1A1P pool best 3.270 vs 3.336) and S9 (`docs/FINDINGS.md:1526, 7261`)
  attributes 0.004 to 0.005 A residuals on 126-target MEANS between the stable and the plain
  argsort to tie-break noise. Those are two conventions compared once, at the mean.
- No sprint measured the per-target sd of the production endpoint (built chain) over random
  tie-breaks, the sd of its mean, or the paired MDE between two tie-breaks, and no small effect
  in the record is quoted beside that floor. The rule "never break a tie by array order"
  (contract rule 8) governs the argmin over scores; the pool boundary is the same trap one
  stage earlier, and it is inside production.

## Exact falsifier

Native-free first (before the gate): 8 seeded random tie-breaks per target (the tie class at the
boundary score re-drawn uniformly, `s15.seed.stable_rng`), the emissions stored, and the
triangle bounds between each draw's built chain and the production chain (RMSD between the two
chains, a rigorous bound on the change in RMSD-to-native, `s26/w_selfcopy.py`'s `triangle`).
Gated: per draw, the 126-target mean of the built chain; s_tie (per-target sd over draws),
m_tie (sd of the 126-mean over draws), and for each pair of draws `ST.compare(draw_a, draw_b,
folds)`: the MDE of that contrast is the floor. Registered predictions: m_tie between 0.003 and
0.010 A; per-target s_tie median between 0.03 and 0.08 A; the paired MDE between two draws
0.02 to 0.05 A. The idea's claim is FALSIFIED if m_tie < 0.002 A and the paired MDE < 0.01 A
(then the convention is harmless and nothing in the record needs a caveat); it is CONFIRMED if
at least one of the record's effects listed above is below the paired MDE between two
tie-breaks. Secondary, reported not decided on: whether the top-75 membership changes with the
tie-break (fraction of the 75 replaced, expected about 75 x 11% = 8 members) and whether the
production convention is systematically better or worse than a random draw (the corpus order
is not neutral, S21 L1459: it is `pdbs/` alphabetical, and a sign here would be a finding about
the pinned file order).

## Expected effect against the computed MDE

This idea measures an MDE; it has no effect to clear. Its deliverable is a sentence for the
report with three numbers (m_tie, median s_tie, the tie-break paired MDE) and the list of
recorded effects that sit below them.

## Memory and agent-hours

One universe at a time plus the s12 distogram cache: the retrieval probe of the sibling script
ran at 0.104 GB (`s26/jobs_done/w_selfcopy_retrieval_probe.json`); 126 targets x 8 draws x one
projection (~3 s) = 50 min CPU, checkpointed per target; est-ram 0.5 GB. Gated part: seconds.
Agent-hours: 1 code (the emission path already exists in `s26/w_selfcopy.py`), 1 run, 1 write-up.

## Information value

High for the report and the presenter whatever the number: either "the pipeline's own
arbitrary tie-breaking moves the answer by X A, and these N recorded effects are inside that",
or "the floor is 0.00X A and every recorded effect stands clear of it". It also converts
contract rule 8 from a rule about argmins into a measured property of the production pipeline.
