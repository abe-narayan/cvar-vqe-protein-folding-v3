# WORKSTREAM C -- SPRINT 23 PRE-REGISTRATION (H_C1, H_C2)

Filed BEFORE any RMSD from either experiment is read, per BRIEF Rule 0. Both experiments build
directly on established, cited facts rather than rediscovering them (see refs at each item).

Instrument: `s12.instrument` (targets/pool/distogram/score, the shared 126-target dev
instrument), `s21/results/poolgap.json` (cached full-pool classical m-ladder, 126/126, verified
`avg_75` mean = 3.0483, bit-matches the incumbent), `core.quantum.run_cvar_vqe` /
`StatevectorCircuit` / `cvar_from_probs` (the genuine CVaR-VQE machinery, unmodified). No AMBER,
no OpenMM lock needed. Benchmark stays sealed; only the 126 dev targets are touched.

---

## H_C1 -- PROBABILITY-WEIGHTED CVaR-VQE READOUT

**Object.** Candidate-identity register, `H|i> = E_i|i>` diagonal on the shipped distogram
Bayes-risk score, `n_qubits = 9` (K<=500 -> dim=512, 12 padding states at a penalty energy),
exact `StatevectorCircuit`, exact `cvar_from_probs`. Continues Sprint 22 Workstream A's exact
setup (its D8b names this as "the one door the theorem does not close" -- untried, not the
shipped `average_weighted` operator, which weights the FULL support rather than the CVaR tail).

**Readout under test.** `X_out = sum_{i in tail} (m_i / sum m) * X_i`, `m_i` = the per-state CVaR
tail mass from `cvar_from_probs(E, p_theta, alpha)` (fractional at the boundary state, by that
function's own convention). Compared against, in the order of importance BRIEF specifies:

  a. **THE BAR** -- classical top-m coordinate average, uniform weights, m = the TRAINED tail's
     own realised size (size-matched per target, not a fixed m=75; s22 A2's own convention).
  b. unweighted average over the identical trained-tail SET (uniform weights, same members).
     By s22 D8's set-equality theorem this is expected to be numerically IDENTICAL to (a) --
     verified as Gate 1 below, not assumed.
  c. **THE CONTROL THAT MATTERS MOST** -- the untrained circuit (same seed, same architecture,
     `iters=0` through the identical `run_cvar_vqe` code path) weighted the same way, at ITS OWN
     realised tail (generally a different size from the trained tail).

**Hyperparameters, fixed before any RMSD is read.** `n_qubits=9`, `layers=3`, `iters=80`,
`lr=0.15`, `restarts=1` (`core.quantum.run_cvar_vqe` defaults). `alpha=0.15` (Sprint 21/22's
canonical "incumbent-analogue" alpha). `T`: PRIMARY `T=0.5` (continues s22 D5's own T-sweep --
its richest, least-collapsed cell, chosen for statistical power to see a weighting margin, not
for RMSD outcome); SECONDARY `T in {0.1, 0.2}` reported regardless of sign. 4 seeds per target
(0-3), matched between trained and untrained arms; per-target value = mean over seeds. n=126 dev
targets, all 5 pinned folds. Basis: POINT CLOUD throughout (never mixed with built-chain).

**Six-axis operator fork, per Rule 0.**

| axis | DECLARED | NOT TAKEN |
|---|---|---|
| functional | diagonal H = shipped distogram Bayes-risk score | Legacy/AMBER H (BRIEF names the shipped score; AMBER needs no lock this way) |
| basis | point-cloud CA of the pool's own windows | built-chain (adds a second, projection, confound) |
| readout | tail-mass-weighted average vs (a)/(b)/(c) above | full-support weighting (the shipped `average_weighted`, a different, already-existing operator -- not the untested object) |
| normalisation | raw Angstrom RMSD | per-target z-scored RMSD |
| null | the three matched controls (a)/(b)/(c) themselves | a size-matched random-tail null (already closed, s21 L18/s22 A1; would dilute the specific weighting question) |
| THE LABEL | Cα-RMSD to native, point-cloud, Kabsch | the CVaR value itself (BRIEF: a lower CVaR that doesn't lower RMSD is not a result) |

**Gate 1 (soundness, reported before any RMSD).** Trained-tail candidate SET ==
`argsort(scores)[:m_trained]` classical SET, exact equality, every (target, seed, T) cell. This
reproduces s22 D8's theorem; any miss is a bug to find, not a result.

**Falsifier.** Primary = `W_trained - avg_matched_trained` (arm a). NOT MEASURED / null if
`|mean diff| < MDE` or the CI includes zero. If it clears AND `W_untrained` also ties its own
bar (c vs its own a-analogue), the effect is weighting-not-training -- reported as the clean
negative BRIEF invites ("If (c) ties (a) ... training is not [doing the work]").

**Registered expectation, stated before running.** PESSIMISTIC. Every other lever on this exact
register (s22 D1-D8) is either a-priori closed or an empirical null. A linear (weighted-mean)
functional over an already energy-ordering-determined SET has no obvious channel to move the
point estimate unless the trained amplitude PATTERN within the tail carries RMSD-relevant
information the rank order does not -- which nothing in this project's closure history predicts.
Expect trained ~= untrained ~= classical bar, i.e. a clean null on both counts. Reported either
way.

---

## H_C2 -- rg_z-CONDITIONED AGGREGATION WIDTH

**Object.** The classical top-m coordinate-average LADDER, m in {500,150,75,20,5,1}, exactly as
cached in `s21/results/poolgap.json` (avg_75 = incumbent, bit-verified). BRIEF treats
"conditioning alpha" and "conditioning tail width" as the same lever for this purpose: `m ~=
ceil(alpha*K)`, and s22 D8 proved the CVaR tail's MEMBERSHIP is always a classical top-m prefix
at whatever alpha training realises -- so conditioning the classical ladder directly is the
size-matched, mechanism-identical, and vastly cheaper operationalisation of conditioning alpha
inside the trained VQE. No retraining is needed to test this lever; retraining would test the
same object at higher cost, not a different one.

**Signal under test.** `rg_z = (rg_disto - rg_pool_mean) / rg_pool_sd`, freshly computed for all
126 dev targets (the `s21/rgsign.py` recipe; that file's own cached artefact only covers n=75,
D's n<=13-residue panel, so this is a genuine extension to the full instrument, not a reread).
`rg_z` is L27/L28's PRIMARY arm: partial rho 0.34-0.38 against ORACLE ordering-quality, survives
length stratification (impossible-by-construction per D's structural argument), a native-free
difficulty control (`pool_spread`), AND the pool's own realised extension -- the strongest
native-free signal this project has produced, and specifically NOT a difficulty proxy
(`r=+0.126` with true difficulty, per BRIEF). `rg_gap` (survives the same controls) is carried as
a secondary/robustness arm, computed identically.

**Router class, deliberately minimal per BRIEF's non-negotiable bound (L10, s22).** ONE fitted
scalar: a threshold `tau` on `rg_z`. Two DECLARED-A-PRIORI rungs, NOT fit: `m_lo = 75` (the
already-established global optimum, s22 L5/A2) for `rg_z < tau`, `m_hi = 150` (next rung up) for
`rg_z >= tau`. Direction fixed in advance by L27-L28's mechanism: positive rg_z means the
distogram predicts a MORE EXTENDED structure than the retrieval pool actually realises, which
L19's concentrated-wrong-region finding says is where the objective's top window is most likely
to be a confidently wrong region -- the hedge is to WIDEN, not narrow. This is a 2-arm, 1-global-
threshold router, i.e. the single simplest class in L10's own table (comparable to its "1 global
threshold, K=13 arms" row, gap 0.39 A at n~100/fold -- here even simpler, K=2).

**Nested CV.** Leave-one-pinned-fold-out (5 folds; `fold` field already on every target record).
On each training set (4 folds), grid-search `tau` over training targets' own rg_z values
(midpoints of the sorted training values), MINIMISING training-fold mean RMSD via the m_lo/m_hi
lookup, with the m_lo<->low, m_hi<->high DIRECTION FIXED (never fit). Apply the fitted `tau` to
the held-out fold only; every target is scored on a threshold it never influenced. Aggregate all
5 held-out folds into one paired comparison against fixed m=75.

**Six-axis operator fork, per Rule 0.**

| axis | DECLARED | NOT TAKEN |
|---|---|---|
| functional | the classical top-m ladder (shipped averaging operator) | retraining the VQE's own alpha per threshold (same object by s22 D8, strictly more expensive) |
| basis | point-cloud, matches poolgap.json and the incumbent | built-chain |
| readout | the independent variable -- m assignment as f(rg_z) | -- |
| normalisation | raw Angstrom RMSD, mean over held-out targets | per-target normalised gain |
| null | FIXED incumbent m=75 (already the established global optimum) | RANDOM m per target (already closed at -0.412 the WRONG way, s22 L5; not the right bar for "does rg_z beat the optimum", only for "does routing beat noise") |
| THE LABEL | `rg_z` (L27/L28's primary, strongest-surviving arm) | `rg_disto` / `rg_pool_mean` (both DEMOTED under L28's realised-extension control) |

**Secondary, reported regardless of sign.** (i) `rg_gap` in place of `rg_z`, identical procedure.
(ii) The REVERSED direction (m_hi for low rg_z, m_lo for high) as an explicit sanity/adversarial
control -- if the reversed direction "accidentally" wins, that is evidence of direction
cherry-picking rather than a real signal, and will be reported as such rather than quietly
dropped.

**Falsifier.** If the routed arm does not beat fixed m=75 by more than its own MDE, with BOTH the
iid and fold-clustered CI excluding zero, DOES NOT FIRE -- NOT DEMONSTRATED, joining L7's four
prior router failures (not escaping them).

**Registered expectation, stated before running.** PESSIMISTIC. L10's own bound: a single global
threshold's generalisation gap (0.39 A) is comparable to the ENTIRE 0.482 A routing ceiling at
n~100/fold. Three prior routers using comparably-sized native-free correlations (r=0.36-0.43)
already failed or reversed sign (s22 L7). rg_z's correlation (0.34-0.38) is not larger than those.
Expect NULL or NOT MEASURED. This is registered as the a-priori likely outcome, not discovered
after the fact -- the reason to still run it is that rg_z is mechanistically different (a
compactness-disagreement signal, not a functional of the objective's own score distribution) and
has never been tried as a router over m specifically, so a positive result, if it occurred, would
be new; a null extends L7/L10 to a fifth independent construction rather than repeating them.

---

## SHARED

- MDE = 2.8016 x SE, per comparison, reported beside every effect.
- Paired bootstrap CI (iid, 4000 resamples) AND fold-clustered CI (resample the 5 fold-clusters
  with replacement, 4000 resamples) beside each other for every primary.
- W/L counts and worst-target degradation reported for every primary.
- ORACLE vs ACHIEVED vs FLOOR labelled at every appearance; point-cloud/built-chain never mixed.
- Completion flags require the full key set (n=126, all folds, all declared T/seed cells), not a
  row count.
- Preserving Pillar 1: both experiments keep the CVaR-VQE (H_C1 directly; H_C2 via the proven
  membership-equivalence to the trained tail, s22 D8) as the selector under test. What is
  QUANTUM in H_C1: the trained amplitude pattern that sets the intra-tail WEIGHTS (untested by
  any prior closure). What is NOT: the tail's MEMBERSHIP (s22 D8, closed a priori) and, in H_C2,
  the routing threshold itself (a classical scalar fit on a native-free signal, exactly the
  BRIEF-mandated router class -- the VQE enters only insofar as its tail membership is what m
  operationalises).
