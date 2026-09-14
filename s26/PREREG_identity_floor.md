# PREREG_identity_floor -- THE SAME SEQUENCE IN A DIFFERENT DEPOSIT (lane W, Sprint 26; tournament rank 2)

Written 2026-09-13 19:55, before the statistics script runs. The measurement itself was
pre-registered as Part E of `s26/PREREG_selfcopy_bound.md` (section 3, falsifier F5) and ran as
job `w_selfcopy_floor` (exit 0, 5 s; `s26/results/w_selfcopy_floor.json`, complete 22/22,
provenance e480fc15); its median is already on the record in L44. This file registers, before
they are computed, the two paired contrasts the idea file said would be reported beside the
median, so that the ledger entry can carry `ST.fmt` verbatim. ORACLE DIAGNOSTIC throughout: every
quantity reads the native; nothing selects; nothing is deployable.

## Hypothesis and the already-registered falsifier

H: at 9 to 16 residues the same sequence in a different deposit is not a near-native answer:
median cross-deposit CA-RMSD above 1.5 A, fewer than a third of pairs below 1.0 A. F5 (already
registered): FALSIFIED if the median is below 1.0 A. Status: measured, median 2.908 A, 18% below
1.0 A (L44). F5 did not fire.

## The two paired contrasts registered here (reported, not decided on)

Unit of analysis: the target (18 targets with a verbatim relative; where a target has several
partners the copy's RMSD is averaged over its partners first, so n = 18, one row per target).
Both quantities on each side are CA-RMSD to the target's native after superposition; basis
"single window against the native", the same basis on both sides.

1. copy minus POOL BEST: the copy's RMSD minus the ORACLE best member of the target's K = 500
   pool. Registered expectation: positive (the copy is worse than the best window), median above
   +1.0 A.
2. copy minus POOL MEAN: the copy's RMSD minus the mean RMSD of the K = 500 pool. Registered
   expectation: negative (the copy is better than a typical window) but with the fold CI
   including zero at n = 18, i.e. NOT MEASURED; power stated.

`ST.compare(copy, other, folds=ST.pinned_folds(pdbs), names=pdbs)`; `ST.fmt` verbatim in the
ledger; the fold-clustered CI decides; MDE = 2.8016 x SE. Both contrasts are ORACLE on both sides
and are context for the median, not results about the pipeline.

## Cost

`s26/w_identity_floor_stats.py` reads one JSON; under 0.1 GB, seconds; run through jobrun as
`w_identity_floor_stats`. No native is re-read (the per-row RMSDs are already in the artefact).
