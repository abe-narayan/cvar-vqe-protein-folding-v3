# IDEA_window_provenance -- WHOLE-PEPTIDE, TERMINAL AND INTERIOR WINDOWS ARE NOT THE SAME EVIDENCE (lane W, Sprint 26; own idea 3; lowest plausibility of the three, census first)

## Hypothesis

A retrieved window is a segment of a parent chain: a whole short peptide (parent length equals
the target length, real termini on both ends), a terminal segment of a longer peptide (one real
terminus), an interior segment of a peptide, or a segment of a protein fragment. The dev
targets are free peptides with two real termini. Measured native-free on one target this turn
(1CEK, in-session, from the universe's codes matched back to the database): of its 500 pool
members 341 are fragment windows, 6 whole peptides, 38 terminal and 115 interior peptide
segments. The hypothesis is that, inside the shipped top-75, whole-peptide and terminal windows
sit nearer the native than interior and fragment windows of the same score rank, and that a
provenance-stratified readout (class weights chosen leave-fold-out on the built chain, uniform
as the incumbent) moves the built chain.

## Why the record does not already close it

- Uniform readout weights were shown optimal over two families, rank-power and
  distance-to-medoid (ARCHITECTURE 2.5, S25); provenance is not a function of either.
- S19 measured that the peptide corpus carries about 7x the sequence-structure channel of the
  protein fragments "that are 80% of every pool" (memory `sequence-conditioning-hurts-the-
  failures`), which is a statement about the two corpora as wholes; within the peptide corpus
  the three positional classes were never separated, and the universe files do not even store
  the parent or the offset (`s8/generate_univ/*.npz` keeps `org`, a peptide/fragment flag, and
  nothing else about provenance); they can be recovered from the window codes.
- The honest prior is null: the pool's error is 68% common-mode (S23 L9) and a per-class weight
  is a within-pool operator; S25's leakage-oracle grid found nothing but uniform.

## Exact falsifier

Native-free census (before the gate): recover (parent length, offset, class) for every pool
member of every target from the codes; the class composition of the K = 500 pool and of the
shipped top-75 per target. Gated: within the top-75, the ORACLE mean RMSD-to-native per class
paired within target (a class is present on a target or not; paired only where both classes
are present), fold-clustered; then the leave-fold-out class-weighted readout against uniform on
the built chain, beyond its own MDE with the fold CI excluding zero, 5/5 folds, with the S25
random-weight control (weights permuted across members, marginal kept). The idea is REFUTED if
the ORACLE per-class difference within the top-75 is inside its MDE (then no weighting can act)
or if the weighted readout does not clear the falsifier. Plausibility 0.15.

## Expected effect against the computed MDE

MDE for a readout re-weighting on the built chain is about 0.05 A (PREREG_C2 section 4: a small
prior change, sd 0.207 on the cloud). Expected 0.00 to -0.02 A: below the MDE, so the likely
deliverable is the census and the ORACLE per-class table, which is the first provenance
decomposition of the pool on the record.

## Memory and agent-hours

Census: 126 targets x up to 40,000 windows matched by substring against 787 peptides and
6,003 fragments, about 5 min, < 0.3 GB. Gated readout: 126 x 5 weights x one projection = 30
min. 2 agent-hours.

## Information value

Low unless the ORACLE table shows a class difference; the census alone is a one-line fact for
the report ("32% of a pool is peptide-derived, 1% whole peptides").
