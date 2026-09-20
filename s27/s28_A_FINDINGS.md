# S28 LANE A FINDINGS -- THE AMPLITUDE READOUT (a CVaR-VQE read as a SIGNED combination), AND A2 (the objective's local behaviour at the production point)

Sprint 28, 2026-09-14, paused 2026-09-14 23:30 and completed 2026-09-19 (the A2 built chain).
Pre-registration `s27/PREREG_S28_A.md` (addenda 1 to 4, every one written before the run it
governs). Code `s27/s28_A_amp.py` (the readout, the surrogate, the exact gradients, the circuit
and the classical families, the ORACLE ceilings, the three phases), `s27/s28_A_analyse.py`,
`s27/s28_A_objdiag.py`, `s27/s28_A2_local.py`, `s27/s28_A2_analyse.py`, `s27/s28_A2_prodcheck.py`.
Tests `tests/test_s28_A.py` (23) and `tests/test_s28_A2.py` (5), all passing. Ledger entries
S28-L1b, L12, L17, L18b, L19, L23b, L26b, L27, L30, L39; lane D's checks S28-L1, L13, L20, L23,
L27b, L31 (every one STANDS or STANDS WITH CAVEAT, the caveats answered in the entries named;
the check of S28-L39 is lane D's next). Every number carries its artefact path; nothing is
quoted from memory. Basis is stated on both sides of every contrast; the built chain is the
verdict basis.

STATUS: complete. Lanes A and A2 close with S28-L39 (section 4.3).

## 1. The question

Every operator the project ships is convex. The real amplitudes of the deployed 9-qubit,
depth-3, 27-parameter RY/CNOT circuit are signed, and read as affine weights
w_i = psi_i / sum_j psi_j over the 500 pool windows (posed on the DIS-top-75 medoid, the
deployed average's own frame) they give a SIGNED combination C(theta) = sum_i w_i W_i. Two
halves, measured apart and labelled apart: EXPRESSIVITY (ORACLE: the family's ceiling, theta
chosen against the native) and RECOGNITION (deployable: the native-free objective
F_lam = CVaR_alpha(E; psi^2) - T H(psi^2) + lam S~(C), with S~ the shipped distogram Bayes risk
of the emitted structure's own distance map, read by linear interpolation of the shipped risk
table and extended linearly beyond the grid, addendum 1). A2 then asks the LOCAL question at
the pipeline's own output: does the objective's gradient at production point toward the native?

## 2. DEMONSTRATED (deployable, built chain, 126/126)

**The amplitude readout is worse than production at every lam.** `s27/results/s28_A_chain_rows.jsonl`,
`s27/results/s28_A_summary.json :: chain_means, contrasts` (ST.fmt verbatim in S28-L26b):

| arm (built chain, same job, same code path as production) | mean | effect vs production | x MDE | fold CI | folds | W/L | verdict |
|---|---|---|---|---|---|---|---|
| production, re-projected (point cloud 3.0483) | 3.2071 | | | | | | comparator |
| circuit lam 0.3, 80 iterations | 3.4432 | +0.236 | 1.48 | [+0.128, +0.327] | 5/5 | 52/74 | WORSE |
| circuit lam 1 | 3.4330 | +0.226 | 1.29 | [+0.090, +0.356] | 4/5 | 48/78 | WORSE, size in the Type-M zone |
| circuit lam 3 | 3.4708 | +0.264 | 1.59 | [+0.115, +0.369] | 5/5 | 40/86 | WORSE |
| circuit lam 0 (the deployed state read as signed weights) | 6.1366 | +2.930 | 4.36 | [+2.67, +3.18] | 5/5 | 16/110 | WORSE |
| lam chosen leave-fold-out (built chain; folds chose 0.3/1/1/1/0.3) | 3.4726 | +0.266 | 1.54 | [+0.148, +0.374] | 5/5 | 50/76 | WORSE |

Worse at every lam; the size at lam 1 is in the Type-M zone (lane D's wording, S28-L27b). The
per-target minimum over the lam grid is 90% accounted for by its order-statistic null (k_eff
3.00, split-half transfer -7%: NOT A SIGNAL, `chain_lam_grid_bok`). Falsifier F1 does not fire;
F2 is moot; no seed-1 run is owed. Point cloud beside it (intermediate; S28-L18b): +0.335 /
+0.337 / +0.390 at 2.0 to 2.4x MDE, 5/5 folds; lam 0 +84 A (median +11.6); the untrained
circuit read as signed weights +19.1 A (mean of 16 draws; `s28_A_recog_rows.jsonl`).

FAIL18 / 108 on the chain (`ST.compare` per stratum): lam 0.3 is -0.073 on FAIL18 (SE 0.163,
0.16x MDE, 10W/8L, NOT MEASURED) and +0.288 on the 108 (1.73x, fold CI [+0.224, +0.352], 5/5,
WORSE); lam 1 -0.073 (0.13x) / +0.276 (1.54x); lam 3 -0.026 (0.05x) / +0.312 (1.84x). The
harm is carried by the 108 targets where the distogram is right; on the 18 where it is wrong
the signed readout is a coin toss (S27 L9's regime pattern; no native-free switch, S22/S23/S26).

Soundness: uniform psi on the DIS top-75 reproduces the S25/S27 point-cloud anchor 3.048338 on
all 126, max |dev| 5.7e-14; at lam = 0 the loop reproduces `core.quantum.run_cvar_vqe`'s `p`
bit-for-bit on 126/126 (max |dp| 0.0); undefined readouts 0 in every arm (0 of 2016 untrained
draws), so pre-registered rule (e) was moot; NaN-poison: `tests/test_s28_A.py` and lane D's
end-to-end poison of the recognition phase (S28-L10 item 2).

**The classical controls at the same objective (point cloud, `s28_A_summary.json :: cloud`).**
Every SIGNED family is worse than production and worse in proportion to how well it minimises
S~: unconstrained a500 rand-init at matched budget (80 evaluations) 3.5755 / converged (2000)
3.6869 (S~ 1.13 / 1.08); a75 3.4880 / 3.5387; the random 27-dim subspace, mean of 8, 3.5445;
the brief's S-only a500 converged 3.7530 (S~ 1.06). The CONVEX simplex under the same objective
lands at production: rand-init matched 3.0563 (+0.008), from production converged 3.0522; a75
from production at lam 0.3 is 3.0460 (-0.002). Among signed families the circuit is the least
bad (vs a500 -0.19 at matched budget, -0.30 converged, budgets differing 25x; vs the random
27-subspace mean-of-8 -0.16 at 1.13x MDE, Type-M, 4/5 folds, sign only) and it is WORSE than
the convex simplex (+0.329, 2.07x MDE, 5/5). The parameter count does not do it and the family
does not do it: the sign freedom itself is what the objective misuses.

## 3. ORACLE DIAGNOSTIC -- expressivity is not the barrier; the objective is

(Every number here is ORACLE unless marked; none is a result.)

**3.1 The family's ceiling.** `s27/results/s28_A_oracle_rows.jsonl`, `s28_A_chain_rows.jsonl :: oracle_circ`:
| family (ORACLE) | point cloud mean | built chain mean |
|---|---|---|
| production uniform top-75 average (deployable; for scale) | 3.0483 | 3.2071 |
| ORACLE best single member of the top-75 / of the pool | 2.3062 / 1.7108 | |
| ORACLE random 27-dim affine subspace, MEDIAN of 8 (mean of 8) | 0.6030 (0.6083) | |
| ORACLE circuit family, MEDIAN over 5 starts (mean 0.7441; best of 5, an order statistic, 0.2884) | 0.4468 | best of 5: 0.2516, 126/126 under 2 A, worst 0.970 |
| ORACLE affine hull of the top-75 / of all 500 (least squares) | 0.0000 / 0.0000 | affine-500: 0.1378 on 22 |

- The affine-hull rows are zero by DIMENSION COUNTING (3n <= 48 coordinates, 75 generic windows
  span the whole space): S10-5's "affine span 0.035 A" was this tautology, and a signed readout's
  expressivity question is only ever about a LOW-PARAMETER family.
- The 27-parameter circuit family reaches 0.288 A (best of 5 local optimisations) and 0.447 on
  the median start; like for like (lane D's S28-L13, reproduced), the circuit's MEDIAN-of-5
  beats the random linear subspace's MEDIAN-of-8 exact optimum by -0.156 A (2.27x MDE, fold CI
  [-0.179, -0.130], 5/5, 94W/32L); on the MEAN it is worse by +0.136 because start seed 4 sits
  in a 1.97 A basin on every target (per-start means 0.385 / 0.584 / 0.416 / 0.362 / 1.973;
  `best_of_k_within` on the starts: 116% accounted, k_eff 3.38). The best-of-5 figure carries a
  best-of-K premium of about 0.1 A.
- The ORACLE optimum survives the projection at a GAIN: 0.288 point cloud -> 0.252 built chain
  (-0.037), 126/126 under 2 A; production pays +0.159 for the same projection. The signed
  optimum is not contracted (Rg 6.59, mean virtual bond 3.77; pool 6.80 / 3.81; production
  6.21 / 2.96); the recognition arms' signed optima are not contracted either (bond 3.40 to
  3.54) and pay only +0.03 to +0.06 for the projection. The signed readout escapes the averaging
  contraction (S25's 22%). It does not escape the objective.

**3.2 The objective's ordering (`s27/results/s28_A_objdiag.json`; the native read here only).**
S~ at the native posed in the frame 2.089; at the ORACLE circuit optimum 2.094 (a 0.29 A
structure scores like the native); at production 1.674; at the pool's best-scoring member
1.423; the pool mean 2.487. The native sits at the 36.9th percentile of its own pool under the
shipped objective (S8-9's "37th percentile", reproduced to the decimal on the surrogate).
Production scores BELOW the native on 99/126 targets and below the pool's best member on 6/126.
The ordering is (pool's best-scoring member) < (production) < (native): any optimiser of this
objective, given the freedom of signs, walks away from the native and toward the pool's
best-scoring member. The recognition rows show it happening: at lam 1 the circuit's S~ falls
from 3.41 at theta0 to 1.34 at the optimum, below production's 1.67 and far below the native's
2.09, while the RMSD rises; 400 iterations lower S~ to 1.26 and raise the RMSD further (+0.483
point cloud, 2.7x MDE); the converged unconstrained optimum has the lowest S~ (1.06 to 1.08)
and the worst RMSD (3.69 to 3.77).

**3.3 A2.1, the local statement (`s27/results/s28_A2_cosine_rows.jsonl`, S28-L23b).** At the
production point the cosine between the shipped objective's steepest-descent direction and the
direction to the native (both with rigid-body components removed) is -0.034 (SE 0.021),
median -0.043, 56/126 positive (sign test p 0.25); a random shape field gives mean |cos| 0.140.
FAIL18 -0.143 (SE 0.066, 6/18 positive), the 108 -0.016 (SE 0.022); Spearman(cos, production
RMSD) -0.372: the worse the target, the more the descent direction opposes the native. Not one
S27 channel is locally informative either (RG_LAW +0.034, DISTPOT -0.002, CONTACT +0.030, ENV
-0.006, CAGEO +0.010; the last 15/18 positive on FAIL18 at a mean of +0.04, one of 21 sign
tests, a third of the random reference: recorded, not a result). The distogram's first step
from its own output is uninformative, and where the distogram is wrong it points away.

## 4. A2.2 -- the deployable step ladder (native-free)

**4.1 Point cloud (intermediate; `s27/results/s28_A2_ladder_rows.jsonl`, `s28_A2_summary.json`, S28-L30).**
C(e) = C0 - e g / rms(g), e in {0.1, 0.3, 1.0} A: +0.006 (0.97x MDE, Type-M) / +0.031 (1.7x,
WORSE) / +0.237 (4.1x, WORSE) vs production, 5/5 folds each; against a random direction of the
same RMS displacement (mean of 8 draws) +0.003 / +0.009 / +0.034, under MDE at every e with the
fold CI above zero at every e. S~ falls on 126/126 targets at e = 0.1 and 113/126 at 0.3. The
circuit family passes within 0.114 A RMS of production (its nearest point reproduces production
at -0.001, 0.4x MDE) and its one steepest-descent step degrades like the raw step (+0.004 /
+0.023 / +0.202, the last at 2.2x MDE, WORSE, 5/5). FAIL18 / 108: +0.015 / +0.050 / +0.215 and
+0.004 / +0.028 / +0.240 (raw means). Prior held on every count.

**4.2 What it means.** From production, the shipped objective's descent direction is worth what
a random direction is worth, or a little less; a trust region around production has nothing to
follow. Together with 3.2 and 3.3 this is the sprint's mechanism sentence for lane A: the
objective's minimum is away from the native (S15's 3.321, S28-L18b), its first step from the
pipeline's own output is blind (cosine -0.03), and its ordering puts the native behind the
average on 99/126 targets; a signed readout, which can reach 0.25 A on every target, therefore
walks the wrong way, and so does every convex re-weighting that leaves production.

**4.3 Built chain (the verdict basis; 126/126; S28-L39).** `s27/results/s28_A_chain_rows.jsonl`
(the 13 A2 arms per target; jobs `s28A2_chain` 94/126 on 2026-09-14 and `s28A2_chain_r2` the
remaining 32 on 2026-09-19, 2221.6 s wall with three governor suspensions, peak RSS 0.313 GB),
`s27/results/s28_A2_summary.json :: contrasts (chain:*), chain_strata, chain_projection_price,
chain_rand_best_of_2, chain_circP, chain_step_e_grid_bok, chain_circ_e_grid_bok,
falsifier_A22_chain`. The comparator is production's chain from the primary chain job (S28-L26b,
3.2071): the identical stored cloud through the identical code path, re-projected today on
2LNG / 1A13 / 9BFL / 6RRO bit-for-bit (4/4, `s28_A2_prodcheck.json`), so the S28-L18 / L27b
floor does not apply. ST.fmt verbatim for every contrast in S28-L39.

| arm (built chain) | mean | vs production | x MDE | fold CI | folds | W/L | verdict |
|---|---|---|---|---|---|---|---|
| production (S28-L26b) | 3.2071 | | | | | | comparator |
| gradient step e = 0.1 A | 3.2034 | -0.0037 | 0.21 | [-0.016, +0.007] | 2/5 | 67/59 | NOT MEASURED |
| gradient step e = 0.3 | 3.2249 | +0.0178 | 0.74 | [+0.009, +0.027] | 5/5 | 57/69 | NOT MEASURED (sign only; Type-M zone) |
| gradient step e = 1 | 3.3116 | +0.1045 | 2.10 | [+0.081, +0.128] | 5/5 | 39/87 | WORSE |
| random direction, mean of the SAME two projected draws, e = 0.1 / 0.3 / 1 | 3.2084 / 3.2250 / 3.3255 | +0.0014 / +0.0180 / +0.1185 | 0.08 / 0.99 / 2.64 | | 2/5, 4/5, 5/5 | | NOT MEASURED / NOT MEASURED / WORSE |
| step vs that random mean, e = 0.1 / 0.3 / 1 | | -0.0051 / -0.0001 / -0.0139 | 0.22 / 0.00 / 0.26 | all include zero | | 70/56, 72/54, 70/56 | NOT MEASURED at every e |
| circP (the family's nearest point to production; residual 0.114 A RMS) | 3.2143 | +0.0072 | 0.38 | [-0.010, +0.023] | 3/5 | 65/61 | NOT MEASURED |
| circuit one step vs circP, e = 0.1 / 0.3 / 1 | 3.2143 / 3.2220 / 3.3364 | +0.0000 / +0.0078 / +0.1222 | 0.00 / 0.21 / 1.34 | e = 1: [+0.068, +0.162] | e = 1: 5/5 | e = 1: 45/81 | NOT MEASURED / NOT MEASURED / WORSE |

- The falsifier (some e beats production beyond MDE, fold CI, 5/5, AND beats the random mean
  beyond MDE) does not fire at any e, decided from the stored contrasts. The prior held on every
  measurable count: the step degrades at e = 0.3 (sign) and e = 1 (2.1x MDE); nothing is measured
  at e = 0.1 on either side; the random direction degrades about equally; the circuit's one step
  degrades from its own baseline at e = 1 (1.34x, WORSE).
- Like for like (lane D's S28-L23 (c)): against the mean of the SAME two projected draws the step
  is 0.0x to 0.3x MDE with the fold CI including zero at every e. On the point cloud (4.1) the
  same contrast carried a sign (fold CI above zero, under MDE); the projection removed it. The
  random best-of-2 is 72 to 83% accounted by the across-target null, split-half NOT A SIGNAL.
- The e grid (D's (b)): the per-target minimum over e is 93% (step) / 94% (circuit) accounted,
  k_eff 2.7; the split-half "residual" (-0.043 / -0.038) is the column-mean effect of the
  smallest step (its column mean sits 0.0433 below the average of the three), i.e. the choice
  "e = 0.1", whose own contrast with production is the null -0.0037. No e is chosen.
- FAIL18 / 108 (`chain_strata`; FAIL18 spans 4 folds): at e = 1 the step is worse in both strata
  (FAIL18 +0.134, 1.01x, Type-M, 4/4, 4W/14L; the 108 +0.100, 1.85x, 5/5, WORSE); at e <= 0.3
  nothing is measured in either; the circuit step at e = 0.1 / 0.3 costs +0.041 / +0.117 on
  FAIL18 (1.08 to 1.09x, Type-M, sign only) and 0.05x on the 108: where the distogram is wrong
  its descent costs more (3.3's -0.14 cosine, in Angstroms), and no regime is helped.
- Through the projection: the step's harm at e = 1 shrinks from +0.237 (cloud) to +0.105
  (chain) and the same two random draws' from +0.212 to +0.119, because a 1 A displacement
  de-contracts the average (projection price +0.159 for production, +0.026 for the step, +0.06
  to +0.07 for the random draws; emitted Rg 6.48 production, 6.88 step e = 1, 6.56 random). The
  objective's descent from production is largely an expansion of the contracted average, not
  toward the native (cosine -0.03), and the emitted chain ends where a random step's does. S~
  falls as asked: 1.674 -> 1.619 on 126/126 at e = 0.1, 1.541 on 113/126 at 0.3, 1.533 on
  79/126 at 1.
- Verdict: REFUTED for accuracy ("the objective is locally informative at the production
  point"). 0.1 A from production is free and worthless, 0.3 A costs the sign, 1 A costs 0.10 A,
  and at every size the price is a random step's. A trust region around production has nothing
  to follow. No seed-1 run is owed. Lanes A and A2 close here.

## 5. HYPOTHESIS / REFUTED / OPEN

- REFUTED: "a signed-amplitude readout of the CVaR-VQE, under the deployed objective plus the
  emitted structure's own distogram score, moves the built-chain RMSD" (F1 silent, every arm
  worse; S28-L26b, D's S28-L27b STANDS).
- REFUTED: "the circuit's inductive bias regularises the objective's failure" (F2 moot; the
  circuit is the least bad signed family and worse than the convex simplex; A2's one-step arm
  degrades like the raw step).
- REFUTED (A2 prior held, built chain, S28-L39): "the objective is locally informative at the
  production point": the deployable step is null at 0.1 A, worse in sign at 0.3, worse at 1
  (2.1x MDE), and at every size worth a random step of the same size (0.0x to 0.3x MDE
  against the mean of the same two projected draws); the circuit family adds nothing.
- ORACLE DIAGNOSTIC, standing: the 27-parameter signed family contains a structure within
  0.29 A (point cloud) / 0.25 A (emitted) of the native on every target; a curved 27-parameter
  family is more expressive than a linear one of the same count by 0.16 A on the median.
- HYPOTHESIS (not tested here): an objective that scores the native below the average on most
  targets would have to be built on something other than the shipped per-pair posterior
  medians; S24 L13's prior-derivative ladder says how much a better prior is worth and S28 A
  says the readout is ready for it (expressivity 0.25 A emitted). Nothing in this sprint shows
  such an objective exists (S27: 203 configurations, none).
- OPEN: nothing quantum. The amplitude readout is the first non-diagonal USE of the state (its
  signs) on record; the state's signs at lam = 0 are the initialisation's (S28-L18b, D's (d)),
  and at lam > 0 they are the objective's, which is classical and wrong.

## 6. What damaged my own expectations

- I expected the affine-500 ceiling to be a measurement; it is a tautology (3n + 1 <= 49 < 75).
  The bound ladder's 0.035 A was dimension counting; I should have seen it before running.
- I expected the best-of-5 ceiling to be the number to quote; lane D's like-for-like check
  (S28-L13) was right that the median is, and that one initialisation basin is bad everywhere.
- I expected the circuit family NOT to contain production; it passes within 0.11 A RMS of it on
  every target, so the "family regularises" hypothesis had no room to work at the local scale.
- Three ledger numbers collided with other lanes' entries (L1b, L18b, L26b) and one append was
  lost to a write race (L23b): four lanes appending to one file need read-and-append in one
  process, which I adopted from S28-L27 on.
- Times in S28-L12 and prereg addendum 3 were stamped ahead of the clock (about 19:45, not
  20:05 / 20:10); nothing else in them changes.

## 7. What I did not do and why

- No AMBER; no touch of `core/`; no parameter chosen on a native quantity; no seed-1 run
  (no positive to replicate; the target-order reversal is a determinism check only, D's (b)).
- The classical controls (a500, a75, the subspaces, the simplex) were not projected to the
  built chain: every one is worse than production beyond its MDE on the point cloud and F2 is
  moot (prereg section 4 and addendum 2's rule: an arm that clears 0.7x MDE in the WORSE
  direction has nothing to decide on the chain; the projection has never reversed a sign on
  this record, and the 756 projections would have cost 1.5 to 8 hours under the user's load).
- The affine-500 emitted ceiling was dropped after 22 targets (0.138 on 22; S10-5 has 0.064 on
  126): its point cloud is 0.000 by construction and the budget went to the arms that decide.
- The 400-iteration circuit cells were run at lam 1 only, as a convergence diagnostic (addendum
  2): S~ falls further and the structure gets worse, which is the whole point.
- On the A2 built chain only random draws 0 and 1 per e are projected (16 projections per
  target would be the budget's limit; D's S28-L23 caveat (c): the step is paired against the
  mean of the SAME two draws).
- Production was not re-projected inside the A2 chain job (the arm list is the 13 A2 arms);
  the comparator is the primary chain job's projection of the identical stored cloud by the
  identical code path, and instead of asserting determinism I measured it on four targets
  including the branch-flip target 2LNG (bit-identical 4/4, `s28_A2_prodcheck.json`, D's (d)).
- No seed-1 replication of the random control's draws (no positive to replicate); no e was
  chosen (the grid's split-half transfer is the column-mean effect of the smallest step).

## 8. Resources

Jobs (all under `s26/jobrun.py --agent S28A`, tag CPU; peak RSS from `s26/jobs_done/*.json`):
`s28A_oracle_126` 693 s / 0.322 GB; `s28A_recog_126` + `_r2` 713 + 1420 s / 0.336 GB (one
governor kill at 43/126, resumed from the per-target checkpoint); `s28A_chain_oracle` + `_r2`
(killed at 20/126 and 22/126 by the governor at 95 to 96% RAM from the user's load; both
resumed); `s28A_chain_primary_1` 4262 s / 0.311 GB (one suspension); `s28A2_cosine_126` 60 s /
0.298 GB; `s28A2_ladder_126` 436 s / 0.26 GB; `s28A_objdiag` 30 s / 0.26 GB; `s28A2_chain`
(2026-09-14, 94/126, 13 projections per target, died with the harness's low-memory kill of
its session at the pause) + `s28A2_chain_r2` (2026-09-19, the remaining 32 targets, 2221.6 s
wall of which three governor suspensions under the user's 93 to 95% box load, peak 0.313 GB);
`s28A2_prodcheck` 25 s / 0.308 GB. Every job
checkpoints per target and is resumable; every result file is provenance-stamped.
