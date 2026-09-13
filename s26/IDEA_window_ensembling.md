# IDEA_window_ensembling -- TEST-TIME ENSEMBLING OVER RETRIEVAL VARIANTS AT FIXED K (lane P, Sprint 26)

## Hypothesis

Retrieve the K = 500 pool three times with different keys (BLOSUM45, BLOSUM62, BLOSUM80; a
second variant over K in {250, 500, 1000} at the shipped BLOSUM62), select each pool's top-75 with
the shipped distogram score, average each top-75 in its own medoid frame, superpose the three
clouds and average them, then project once (`s12.instrument.project`, ramah 0.3). The ensemble's
built chain beats the shipped built chain.

## How this differs from widening K, and why it might not just reproduce S17 L12

S17 L12 widened ONE shortlist: from K = 75 to the full universe the top-75's mean improves
(4.287 -> 3.557) and its best is destroyed (2.104 -> 2.572) because extra plausible windows
DISPLACE near-native members from the single top-75; the consensus readout is best at K = 500
(3.282) and worse at 2000 / full (3.344 / 3.461). Ensembling at FIXED K never widens any
shortlist: each of the three top-75s keeps its own near-native members (top-75 ORACLE best
2.306 at K = 500 is untouched in the BLOSUM62 member), and only the OUTPUTS are combined.
Mechanically it is a 225-member average with three score-selected shortlists, not a 75-member
average from one shortlist over 1,500 candidates. If the three shortlists' errors were
independent the average would cancel part of them; the S17 mechanism (displacement) cannot act.

## Why the record probably closes it anyway (stated before running)

- The pool's error is 68% common-mode (S23 L9) and two score-selected sources of DIFFERENT
  provenance have bias cosine 0.943, above the within-source control 0.933 (S24 L3): what the
  shipped score selects lands on the same error direction whatever the candidate source. Three
  BLOSUM keys are three sources of the same provenance; their clouds are predicted parallel, and
  averaging parallel errors removes nothing.
- The m-ladder on the shipped row is flat over m in [35, 110] and m = 75 is its argmin (S12
  agg_FINDINGS 5a: m = 110 3.050, m = 150 3.062): a 225-member average of same-quality members
  is predicted slightly WORSE, not better.
- Retrieval is real but nearly saturated as a lever: fixing it is worth 0.016 A (S8-6); BLOSUM62
  beats random at every K and the key stops mattering by K = 2000 (S7-12).
So the recorded expectation is null-to-slightly-worse, and a gain would have to come from
BLOSUM45/80 retrieving BETTER shortlists on average, which S8-6 bounds at ~0.02 A.

## Exact falsifier

Built chain, paired per target: ensemble minus shipped beyond its own MDE, fold-clustered CI
excluding zero, 5/5 folds. Two controls, same operator: (i) three bootstrap 75-subsets of the
SAME BLOSUM62 top-75 averaged the same way (a zero-information ensembling control matched in
member count), and (ii) the single BLOSUM45 and BLOSUM80 clouds alone (are the variants better
shortlists at all). The point cloud is carried; the selection basis is not relevant (there is no
single selected candidate).

## Expected effect vs computed MDE

The variant clouds differ from the shipped one by a re-ranking of the same universe, so the
paired sd should resemble a small prior change: MASSFIXW0.1 - MASS0.0 sd 0.207, MDE 0.052 on the
point cloud (PREREG_C2 section 4); on the built chain ~0.05-0.06. Expected effect 0.00 to +0.03:
below the MDE, i.e. the design can only refute a gain of 0.05 or more.

## Why it is cheap enough to run anyway

Everything comes from the cached universes: `s8/generate_univ/<pdb>.npz` stores every window's
codes `S` (int8, ALPHABET order), so BLOSUM45/80 sums are `M[q, S].sum(-1)` with Biopython's
`substitution_matrices` (installed; verified 2026-09-13) re-ordered to `core.data.ALPHABET`; the
K = 500 and top-75 rules are `core.data.top_k`'s stable argsort and the shipped score
(`I.shipped_score`); the readout is `s24/residlib.readout`; one projection per target per arm
(2.9 s). ~20 minutes of machine time, < 0.6 GB, 2 agent-hours including the write-up.
