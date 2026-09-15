# S28 LANE A FINDINGS -- THE AMPLITUDE READOUT (a CVaR-VQE read as a SIGNED combination)

Sprint 28, 2026-09-14. Pre-registration `s27/PREREG_S28_A.md` (addenda 1 to 3, all written
before the runs they govern). Code `s27/s28_A_amp.py` (the readout, the surrogate, the exact
gradients, the circuit and the classical families, the ORACLE ceilings, the three phases),
`s27/s28_A_analyse.py` (every contrast through `s24.stats_lib.compare`, printed with `ST.fmt`),
`s27/s28_A_objdiag.py` (the objective diagnostic). Tests `tests/test_s28_A.py` (23 pass).
Ledger entries S28-L1b, L12, L17, and the entries that follow them. Every number below carries
its artefact path; nothing is quoted from memory.

STATUS OF THIS FILE: the ORACLE half (section 2) is complete on 126 targets on the point cloud
and PARTIAL on the built chain; the recognition half (section 3) is PARTIAL (43/126, point
cloud) while the governed jobs queue behind the coordinator's launch hold. Sections marked
PARTIAL carry no verdict. This file is rewritten when the runs land.

## 1. The question

Every operator the project ships is convex. The real amplitudes of the deployed 9-qubit,
depth-3, 27-parameter RY/CNOT circuit are signed, and read as affine weights
w_i = psi_i / sum_j psi_j over the 500 pool windows (posed on the DIS-top-75 medoid, the deployed
average's own frame) they give a SIGNED combination C(theta) = sum_i w_i W_i. Two halves,
measured apart: EXPRESSIVITY (ORACLE: the family's ceiling, theta chosen against the native) and
RECOGNITION (deployable: the native-free objective F_lam = CVaR_alpha(E; psi^2) - T H(psi^2)
+ lam S~(C), with S~ the shipped distogram Bayes risk of the emitted structure's own distance
map, read by linear interpolation of the shipped risk table).

## 2. ORACLE DIAGNOSTIC -- expressivity: the family reaches 0.29 A; the affine ceilings are vacuous

(Every number here is ORACLE. `s27/results/s28_A_oracle_rows.jsonl`, 126 rows;
`s27/results/s28_A_summary.json :: oracle_cloud`; job `s26/jobs_done/s28A_oracle_126.json`,
693 s, peak RSS 0.322 GB.)

Soundness: uniform psi on the DIS top-75 reproduces the S25/S27 point-cloud anchor 3.048338 on
all 126 targets, max |dev| 5.7e-14 (`s28_A_summary.json :: prod_frame_max_dev`).

| family (ORACLE unless marked) | mean | median | < 2 A |
|---|---|---|---|
| production uniform top-75 average (deployable; for scale) | 3.0483 | 2.8373 | 0.294 |
| ORACLE best single member of the top-75 | 2.3062 | | |
| ORACLE best single member of the pool | 1.7108 | | |
| ORACLE random 27-dim affine subspace of the 500 windows, MEDIAN of 8 | 0.6030 | 0.5281 | |
| ORACLE random 27-dim affine subspace, mean of 8 | 0.6083 | 0.5235 | 1.000 |
| ORACLE circuit family, MEAN over 5 starts | 0.7441 | 0.6655 | 1.000 |
| ORACLE circuit family, MEDIAN over 5 starts | 0.4468 | 0.3635 | |
| ORACLE circuit family, best of 5 starts (an order statistic) | 0.2884 | 0.2364 | 1.000 |
| ORACLE affine hull of the top-75, least squares | 0.0000 | | |
| ORACLE affine hull of all 500, least squares | 0.0000 | | |

- The affine-hull rows are zero by DIMENSION COUNTING: a structure has 3n <= 48 coordinates,
  75 generic windows span the whole space, so any affine family with at least 3n + 1 members
  contains every structure exactly. That is what S10-5's "affine span 0.035 A" was. The
  expressivity question for a signed readout is therefore only ever about a LOW-PARAMETER
  family; the 500-weight and 75-weight affine families are complete and say nothing.
- The 27-parameter circuit family reaches 0.288 A (best of 5 local optimisations, Adam 300
  iterations) and 0.447 on the median start. Like for like (lane D's S28-L13, reproduced here:
  `ST.compare` on the per-start and per-subspace columns): the circuit's MEDIAN-of-5 beats the
  random linear subspace's MEDIAN-of-8 exact optimum by -0.156 A (2.27x MDE, fold CI
  [-0.179, -0.130], 5/5 folds, 94W/32L); on the MEAN the circuit is worse by +0.136 because start
  seed 4 sits in a 1.97 A basin on every target (per-start means 0.385 / 0.584 / 0.416 / 0.362 /
  1.973). The best-of-5 figure carries a best-of-K premium (`best_of_k_within`: 116% accounted,
  k_eff 3.38). The comparative claim is the median one: a curved 27-parameter family is more
  expressive than a linear one of the same count by about 0.16 A, and the circuit's optimum is a
  local one against the subspace's exact one, so the contrast is conservative in the circuit's
  favour.
- Per length n = 9..16: circuit 0.120 / 0.114 / 0.197 / 0.209 / 0.227 / 0.373 / 0.468 / 0.455;
  random 27-subspace 0.000 / 0.000 / 0.193 / 0.496 / 0.634 / 0.878 / 1.089 / 0.972 (at n <= 10
  the 26 affine degrees of freedom exceed the 3n - 6 shape degrees of freedom, and both are
  complete).
- The ORACLE optimum is a heavily cancelling signed combination: 49% of the weights negative,
  negative mass sum |w_-| 10.6, effective members 1 / sum w^2 = 0.78, denominator sum psi 0.77.
  It is NOT contracted: Rg 6.59 (pool 6.80, production 6.21), mean virtual CA-CA bond 3.77
  (pool 3.81, production 2.96). The convex average shrinks the chain by 22% (S25); the signed
  combination that reaches the native keeps native bond lengths.
- Built chain, PARTIAL 22/126 (`s27/results/s28_A_chain_rows.jsonl`): the ORACLE optimum
  projects at 0.2934 (point cloud 0.2865 on the same 22; price +0.007), 22/22 under 2 A;
  production on the same 22 projects at 2.9325 (point cloud 2.8272; price +0.105). The complete
  126-target emitted ceiling is posted when the chain job lands.

Reading: EXPRESSIVITY IS NOT THE BARRIER. The family contains structures 0.3 to 0.4 A from the
native on every target, and they survive the projection at a price of 0.007 A.

## 3. RECOGNITION -- does the native-free objective find them?  (PARTIAL 43/126, point cloud; no verdict)

The registered reason to expect failure, measured before any recognition arm was read (S28-L1b):
the objective term S~ at the ORACLE optimum is 2.094 against 1.674 at production, and the ORACLE
structure is preferred on only 25/126 targets (`s28_A_summary.json :: S_at_oracle_vs_prod`).
The objective ranks a 0.29 A structure behind the 3.05 A average on 80% of targets.

Partial per-target picture on 43 targets (`s27/results/s28_A_recog_rows.jsonl`, sorted-pdb
order 1A13..2LNG; mean paired d vs production on the same 43; no CI, no MDE, NOT A RESULT):
circuit lam 0 +84.9 (median +23.0; the deployed state's signs are the initialisation's and the
affine sum blows up); circuit lam 0.3 / 1 / 3 at 80 iterations +0.31 / +0.35 / +0.40; lam 1 at
400 iterations +0.42; a500 rand matched / converged (lam 1) +0.56 / +0.65; a75 rand +0.43 /
+0.52; simplex rand matched +0.03; simplex from production, converged, +0.003; S-only a500
converged +0.72; untrained draw 0 +12.2 (median +5.6).

Scale of the terms (lane D's caveat (a); `parts0` / `parts` in the rows; mean over 43): at
theta0 CVaR -1.38, T H 2.01, lam S~ 1.02 / 3.41 / 10.23 for lam 0.3 / 1 / 3; at the optimum
CVaR -1.53, T H 2.9, lam S~ 0.42 / 1.27 / 3.74. The grid spans minor to dominant.

The mechanism visible already: at lam 1 the optimiser drives S~ to 1.27, BELOW production's
1.67, while the RMSD rises 0.35 A; the ORACLE structure sits at 2.09. The objective orders
(circuit optimum) < (average) < (near-native). It is minimised; its minimiser is away from the
native. That is S8-9's 37th percentile seen from the other side.

Verdicts, the FAIL18 split, the nested lam choice, the F2 controls on the chain, the seed-1
replication if anything is positive, and the NaN-poison confirmation are written when 126/126
and the chain land.

## 4. What damaged my own expectations

- I expected the affine-500 ceiling to be a measurement; it is a tautology (3n + 1 <= 49 < 75).
  The bound ladder's 0.035 A was dimension counting, and I should have seen it before running.
- I expected the circuit family's best-of-5 to be the number to quote; lane D's like-for-like
  check (S28-L13) was right that the median is, and that the mean is worse than a random
  subspace because one initialisation basin is bad on every target.

## 5. What I did not do and why

- No AMBER; no touch of `core/`; no parameter chosen on a native quantity.
- The affine-500 emitted ceiling was dropped from the chain job after 22 targets: its point
  cloud is 0.000 by construction and S10-5 already reports 0.064 emitted; under the user's
  load each projection costs 40 s and the budget goes to the arms that decide something.
- The target-order reversal is a determinism check only (D's caveat (b)); seed 1 is the
  replication and is run only if a positive appears.
