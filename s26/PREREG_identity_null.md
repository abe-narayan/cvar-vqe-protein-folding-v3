# PREREG: is the "corrected" identity criterion discriminative, or at the null?

Lane I, Sprint 26, 2026-09-13 00:33. Written before the control is run. Diagnostic only: no
RMSD, no native coordinate, nothing selected; it is allowed before the phase gate.

## Context

`s26/i_identity_audit.py` (first run, `s26/logs/i_identity_audit.log`, job
`i_identity_audit`, 15 s, peak RSS 0.052 GB) reproduces the pinned clusters and folds exactly
from an in-memory copy of `core.data.clusters`, then re-clusters under the brief's "corrected"
criterion: the same Needleman-Wunsch match count normalised by the SHORTER sequence, plus a
verbatim-substring test, at the same 0.6 threshold. Result: 470 pinned clusters become 166;
594/787 sequences gain or lose a mate; 71/126 dev targets gain a mate that sits in another
pinned fold; 48/60 benchmark sequences likewise (counts only).

That is far more than the four declared self-copies, and the project memory
(`containment-threshold-is-at-the-null.md`, S24) already records that a 0.6 containment
threshold is at the null for peptides against longer members: a global alignment with free
end gaps reports the longest common subsequence, and LCS(9, 26) over 20 letters is typically
5 to 6, i.e. 0.56 to 0.67 of the shorter length. So the count above may be chance, not leak.

## Hypothesis

H1: the shorter-normalised criterion at 0.6 is at the null: a SHUFFLED dev-target sequence
(same length and composition, no homology) will pass the criterion against about as many
database members as the real sequence does.

H2: the pinned longer-normalised criterion at 0.6 is discriminative: shuffled sequences pass
it against almost no member.

H3: the minimal fix (pinned criterion OR verbatim substring) moves exactly the sequences that
carry or are carried verbatim, and nothing else.

## Exact falsifier

For each of the 126 dev targets: the fraction of the other 786 database members that pass
the criterion, for the real sequence and for 5 seeded shuffles (`np.random.default_rng(seed)`
with seed = 1000 + target index + 126 * k). Rates averaged over targets.

- H1 is FALSIFIED if the shuffled pass rate under the shorter criterion is below 20% of the
  real pass rate (the criterion would then carry real signal and 71/126 would stand as a
  leak count).
- H2 is FALSIFIED if the shuffled pass rate under the longer criterion exceeds 0.5% of pairs.
- H3 is FALSIFIED if the substring-OR-pinned clustering changes the cluster membership of any
  sequence that is neither a verbatim substring of another member nor contains one.

## Expected

Shorter criterion: real and shuffled pass rates within a factor of 2 of each other, both of
order 10 to 30% of members for the short (9 to 10 residue) targets. Longer criterion: shuffled
rate under 0.1%. Substring-OR-pinned: a handful of merges (the four known pairs plus whatever
other verbatim containments exist in the database), no other change.

## Cost

In memory, read-only, about 30 s and under 0.1 GB. Job `i_identity_audit2` under `jobrun`.
Pinned files hashed before and after, compared with `s26/results/pinned_hashes.json`.
