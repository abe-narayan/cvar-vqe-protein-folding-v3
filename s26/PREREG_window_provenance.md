# PREREG_window_provenance -- WHOLE-PEPTIDE, TERMINAL, INTERIOR AND FRAGMENT WINDOWS IN THE POOL AND THE TOP-75 (lane W, Sprint 26; tournament rank 9; census first)

Written 2026-09-13 22:10, before any code ran on a real target. Idea file:
`s26/IDEA_window_provenance.md`. Code: `s26/w_provenance.py` (synthetic tests inside
`s26/w_provenance_test.py`). Results: `s26/results/w_provenance_*.json`. Native-free half
(`census`) reads sequences and codes only; gated halves (`oracle`, `readout`) read natives and
refuse to run unless "PHASE 0 SIGNED OFF" is in the ledger.

## 0. What the record knows

- The universe files store a window's origin only as `org` (peptide database or protein
  fragment), never the parent chain or the offset; a window's class must be recovered from its
  codes against the database and fragment sequences (a length-n window of a peptide of length
  n is a whole peptide; offset 0 or n_parent - n is terminal; otherwise interior). In-session on
  1CEK (L30): of 500 pool members 341 fragments, 6 whole peptides, 38 terminal, 115 interior.
- Uniform readout weights are optimal over two families (rank-power, distance-to-medoid; S25,
  ARCHITECTURE 2.5); removing four deviant members costs +0.142 A (S23 L5). Provenance is not a
  function of either family.
- The peptide corpus carries about 7x the sequence-structure channel of the fragments that are
  80% of every pool (S19, memory `sequence-conditioning-hurts-the-failures`); within the peptide
  corpus the positional classes were never separated. S24 L8: `fragment_db`'s docstring declares
  fragments training-only, their conformation held by contacts outside the window.
- The pool's error is 68% common-mode (S23 L9); any per-class weighting is a within-pool
  operator and the honest prior is null.

## 1. Hypotheses and exact falsifiers

**H_P1 (census, native-free).** The K = 500 pool and the shipped top-75 are mostly fragment
windows; whole-peptide windows are rare (about 1%). Descriptive; no falsifier; the numbers are
the deliverable ("X% of a pool is peptide-derived, Y% whole peptides, Z% terminal").

**H_P2 (ORACLE, gated).** Inside the shipped top-75, whole-peptide and terminal windows sit
nearer the native than interior and fragment windows. Test: per target, the mean ORACLE
RMSD-to-native (`rr`) of each class present in the top-75; paired within target for every pair
of classes both present; `ST.compare` over the targets where both are present, fold-clustered;
the falsifier is that the whole+terminal class beats the fragment class by more than its MDE
with the fold CI excluding zero on 5/5 folds. Registered expectation: NOT MEASURED (a class
difference inside a score-selected top-75 of well under 0.1 A), because the score already
equalises the classes on what it can see. Also reported: the same contrast at the K = 500 pool
level (before the score), where a class difference is expected (fragments are training-only
material; S24 L8).

**H_P3 (ACHIEVABLE readout, gated; runs only if H_P2's ORACLE difference exceeds its MDE).**
A class-weighted uniform readout (weights w_class in {0, 0.5, 1, 2} for the fragment class
relative to 1 for peptide windows, chosen leave-fold-out on the built chain) beats the uniform
readout on the built chain beyond its MDE with the fold CI excluding zero on 5/5 folds, AND beats
the S25 random-weight control (the same weights permuted across the 75 members, 8 draws). If
H_P2 is NOT MEASURED, H_P3 does not run and the idea closes with H_P2's power statement (no
weighting can act on a difference that is not there).

## 2. Operator, exactly

Census: for every universe window in the K = 500 pool (and the top-75 = production `sub`),
decode the codes to a string, look it up in the set of length-n substrings of the 787 peptide
chains (whole / terminal / interior by parent length and offset; a window string found in more
than one parent takes the class of the FIRST parent in database order, and the count of
multi-parent windows is reported) and, for `org == False`, in the 6,003 fragment chains
(fragment class). Unresolved windows are counted, never guessed.

ORACLE contrast: `rr` per member from the universe, averaged per class per target.

Readout (if run): the shipped top-75 with per-member weights by class, weighted mean in the
production medoid frame, one projection (`s12.instrument.project`, ramah 0.3); weights chosen
on the four training folds, applied to the fifth.

## 3. Statistics

`ST.compare` with `folds=ST.pinned_folds`, `ST.fmt` verbatim, both bases where an RMSD-of-emission
is involved (arm PRIMARY, cloud carried); the ORACLE class contrast is on the single-window
basis (`rr`), stated as such. Replication for a positive: the permuted-weight control at a second
seed and the reversed fold order. Power for a null: the MDE of the ORACLE class contrast at the
number of targets where both classes are present.

## 4. Expected effect against the computed MDE

Class contrast inside the top-75: per-target class means over 10 to 70 members each; expected
MDE 0.05 to 0.10 A on the single-window basis; expected difference within +-0.05 (NOT MEASURED).
Readout re-weighting MDE about 0.05 A on the built chain (PREREG_C2 section 4); expected 0.00 to
-0.02, below the MDE.

## 5. Memory, time

Census: 126 targets x 500 windows matched against a hash set of all length-n substrings of
787 + 6,003 chains (built once per length; about 700k substrings of length 9 to 16): under 0.3
GB, about 5 minutes. ORACLE contrast: seconds. Readout: 126 x 4 weights x one projection = 30
minutes if it runs. One-target probe first. Agent-hours: 1.5.

## 6. Operator forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| class definition | whole (parent length == n), terminal (offset 0 or end), interior, fragment; unresolved counted | finer offsets; parent length bands |
| multi-parent windows | first parent in database order; count reported | dropping them |
| ORACLE unit | per-target class mean, paired within target where both classes present | pooling members across targets (pseudo-replication) |
| readout weights | fragment class relative weight in {0, 0.5, 1, 2}, leave-fold-out | continuous weights fitted (S25: leakage-oracle grids found uniform) |
