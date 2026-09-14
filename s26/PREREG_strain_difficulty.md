# PREREG -- THE RELAXATION'S OWN STRAIN AS A NATIVE-FREE DIFFICULTY SIGNAL (lane PH, S26; tournament rank 4)

Written 2026-09-13 19:38, before `s26/ph_strain.py` produced a number. Never edited after;
addenda appended. Idea: `s26/IDEA_strain_difficulty.md`. Code: `s26/ph_strain.py`. Result:
`s26/results/ph_strain.json`. Reads the 126 production records
(`bench_results/cache/1fc9f2dcf489e2fb`) and nothing else; under a minute; < 0.3 GB.

## 1. Question

Do the native-free scalars the production relaxation emits for free predict the per-target
error of the built chain? This is a CALIBRATION question ("does the physics report when the
answer is untrustworthy"), not an accuracy lever: nothing is selected, tuned or moved. The
endpoint `rmsd_arm` is read only as the ORACLE label of a fixed, already-emitted structure.

## 2. Signals (native-free, all from the cached record) and confounds

    s1  log10(max(e0, 1))      the built chain's own AMBER energy before relaxation (its strain)
    s2  log10(max(e0 - e1, 1)) the energy the relaxation removed
    s3  amber_moved             restraint RMSD on N/CA/C, how far the chain had to move
    s4  amber_strain_after      bond + angle energy left after relaxation
    confounds: chain length n; Rg of the built chain (a longer or more compact chain has more atoms
    and more contacts). Each rho is reported raw and partial on (n, Rg), the partial by residualising
    both variables on the two confounds by least squares before the Spearman.

## 3. Statistics and falsifier

Spearman across the 126 targets between each signal and `rmsd_arm`, with (a) a 4000-draw
fold-clustered bootstrap CI (resample the five pinned folds with replacement, recompute rho)
beside the iid bootstrap CI, (b) a 4000-draw label-permutation null (rmsd_arm shuffled across
targets) giving a two-sided p, (c) the sign of rho on each of the five folds computed within
fold. Practical form: the top quartile of each signal (32 targets) against membership of
`I.FAIL18`, one-sided Fisher exact test.

FALSIFIER. A signal exists iff at least one of s1..s4 has |rho| >= 0.25 with the fold-clustered
CI excluding zero, the same sign on 5/5 folds, and the permutation p below 0.05 / 4 (Bonferroni
over the four signals), on the PARTIAL correlation (confounds removed). Otherwise the relaxation's
strain does not predict the error and the idea is closed on this instrument. The raw correlation
is reported beside it and a raw-only positive is a confound, not a signal.

## 4. Expected outcome, stated before the run

A raw rho of 0.2 to 0.3 for s1 driven by chain length (longer chains have more atoms in the
builder's fixed rotamers and higher RMSD), falling below 0.2 once n and Rg are partialled. The
Fisher test on FAIL18 is expected null (FAIL18 is about retrieval recall, not about strain).
Prior for the falsifier clearing: 0.3. At n = 126 a Spearman has SE about 0.09, so the design
resolves |rho| of 0.25 with power about 0.8; anything below 0.2 is UNDERPOWERED, not null.

## 5. Cost

One process, < 0.3 GB, under a minute, tag CPU, no checkpoint needed (a single pass over 126
JSON files). Agent-hours: 1.

## 6. Replication

A positive result is re-run with a different bootstrap seed (`stable_rng` salt "rep") and the
folds processed in reversed order; the CI and the fold signs must reproduce.
