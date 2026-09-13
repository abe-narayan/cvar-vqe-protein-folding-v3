
---
## ADDENDUM 2 (2026-09-13 09:32) -- THE REPLICATION RUN, DECLARED BEFORE IT RUNS

Stage 1 (`s26/results/ph_c3_stage1.json`) returned one positive contrast that is not the arm
under test: the zero-information toward-member displacement BEATS do-nothing on the built chain
(-0.0178, fold CI [-0.0255, -0.0124], 90W/36L). The contract requires every positive result to be
re-run on a different seed and a different fold-processing order and to land inside its own CI.
`python s26/ph_c3.py stage1 --rep` does exactly that: draw seeds salted "rep"
(`stable_rng(pdb, "c3randrep")`, `stable_rng(pdb, "c3memberrep")`), targets processed in reversed
pinned order, fold labels untouched, artefact `s26/results/ph_c3_stage1_rep.json`. The replication
is judged on toward-member-minus-do-nothing and on AMBER-minus-random (the Type-M-zone contrast):
each must land inside the first run's fold CI. The arm itself (AMBER's displacement) has no
randomness, so its contrasts move only through the controls' draws.
