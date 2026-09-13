# IDEA_conformational_identity_floor -- THE SAME SEQUENCE IN A DIFFERENT DEPOSIT: HOW FAR A VERBATIM COPY SITS FROM THE NATIVE (lane W, Sprint 26; own idea 2; Part E of PREREG_selfcopy_bound)

## Hypothesis

At 9 to 16 residues the same sequence in a different deposit does not adopt the same
conformation: the CA-RMSD between a dev target's native (model 1) and the verbatim copy of its
sequence inside another database peptide (the carrier's segment at the substring position, or
the shorter relative's whole chain against the target's segment) is, in the median, above
1.5 A, i.e. above the intra-ensemble spread of the natives themselves (1.044 A mean to model 1,
`docs/FINDINGS.md:2184`) and not far below the pool's best member. Two consequences follow.
First, the mechanism behind a small self-copy leak bound (`s26/IDEA_selfcopy_proxy_bound.md`):
a verbatim training or retrieval copy carries limited structural information. Second, a
statement about the endpoint itself: a perfect sequence-identical LOOKUP is not a 2.0 A
predictor on this instrument, so the 2.0 A target set in S15 to S25 was below the floor that
sequence identity alone can reach for these peptides, and "containment-fresh" (memory
`no-fresh-benchmark-exists`) is a weaker notion of novelty than the record has treated it as.

## Why the record does not already close it

- S24 L4 measured four cross-deposit values on the way to declaring the leak (0.595, 3.278,
  2.334, 4.126 A for 1CEK, 2FBU, 2P5H, 6B9K) and drew the sentence "the same peptide sequence
  in a different deposit adopts a different conformation" from them. n = 4, never generalised.
- The intra-ensemble spread (S7, n = 111, 1.044 A mean; S12: it does not predict the selection
  gap) is WITHIN one deposit; S8-14's de-noised labels are ensemble means within a deposit.
  No sprint measured the BETWEEN-deposit, same-sequence distance systematically.
- Lane I's audit (L15, `s26/results/i_identity_audit.json`) found 18 dev targets with a
  verbatim relative (4 in another fold, 14 in their own), which is the whole measurable
  population without touching a benchmark target; the 20 verbatim containments in the database
  include pairs involving benchmark sequences and are deliberately NOT used.

## Exact falsifier

Gated (reads natives): for the 18 dev targets and their 22 verbatim partners, the cross-deposit
CA-RMSD at the shared segment (`s26/w_selfcopy.py floor`), beside the target's pool mean and
pool best RMSD and the record's 1.044 A ensemble spread. Registered prediction: median above
1.5 A, fewer than a third below 1.0 A. FALSIFIED if the median is below 1.0 A: then verbatim
copies ARE near-native at this length, the dev self-copy bound rests on luck, and only the
n = 126 envelope of the sibling PREREG is quoted for the benchmark. Also reported: the
dependence on partner length ratio (a 13-mer inside a 25-mer against a 12-mer inside a 13-mer)
and on role (carrier or carried), without a claim at n = 22.

## Expected effect against the computed MDE

No paired effect; a descriptive distribution with n = 22, quoted with its min, median, max and
the fraction below 1.0 and 1.5 A. The MDE concept does not apply; the falsifier is a threshold
on the median.

## Memory and agent-hours

Under 0.1 GB, one minute (22 Kabsch superpositions and 18 universe loads for the pool
reference). 0.5 agent-hours. Runs the minute "PHASE 0 SIGNED OFF" lands.

## Information value

One sentence for the report that the record does not yet have: "sequence identity buys X A on
this instrument", with the 2.0 A target placed against it, and the mechanism for why the
declared leak is small stated with a number instead of an anecdote.
