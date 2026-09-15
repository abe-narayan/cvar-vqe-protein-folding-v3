# S27 LEDGER -- ALTERNATIVE HAMILTONIANS (2026-09-14)

## L1 -- SCOPE AND GROUND
Question: can a different Hamiltonian or a principled combination beat the shipped distogram
Bayes-risk score on the 126-target instrument with the production readout and the genuine
CVaR-VQE intact. Pre-registration `s27/PREREG.md` (base + addenda 1, 2) before every wave.
Anchors reproduced: point cloud 3.048338 (DIS top-75, 126/126), built chain 3.2126 (rebuild
basis), the S25 VQE arm 3.0580. Nothing in production, the benchmark or an earlier sprint touched.

## L2 -- WAVE 1 (203 configurations, point cloud): NOTHING BEATS DIS; 45 WORSE BEYOND MDE
`s27/results/pool_summary.json`. Best effect -0.026 A (DIS+0.5*SS_MATCH, 0.35x MDE). Every new
single worse than DIS; most worse than a random 75-subset. DISTPOT (Sippl-style pair potential
fitted on the leakage-safe universe) is the best new single at 3.199 (-0.222 vs random, 0.9x MDE).

## L3 -- MECHANISM: RANKING INFORMATION IS ANTI-USEFUL ON THE AVERAGING READOUT
`s27/results/redundancy.json`, `mechanism.json`. Spearman(partial rho of X with the ORACLE RMSD
given DIS, effect of DIS+X on the top-75 average) = +0.48 over 18 channels; CONS has the most
information (+0.25) and does the most harm (+0.29 A, worse than its own permuted control at
1.7x MDE). The set-mean law (S12) and the outlier-removal cost (S23 L5) in one measurement.

## L4 -- WAVE 2 (m-ladder, in-band, non-additive, near-corpus): NULL; DIS+CONS FALSIFIED AT EVERY m
`s27/results/wave2_summary.json`. No arm beats DIS at any m in 3..150; the pre-registered
small-m prediction for DIS+CONS fired the other way (+0.159 at m = 3). In-band re-ranking null;
CONS worst in-band (+0.184). Per-arm minima over m transfer 19 to 29% in split halves.

## L5 -- WAVE 3: THE SS_MATCH BUMP IS A GRID ARTEFACT
`s27/results/wave3_summary.json`. With the four-weight grid the nested held-out effect is
-0.008 A (0.1x MDE), `best_of_k_within` "NOT A SIGNAL (split-half 10%)"; SS_MATCH2 (H-bond SS
call, universe propensities) no better.

## L6 -- THE GENUINE CVaR-VQE ARM (80 configurations x 126 x 2 seeds) CHANGES NO VERDICT
`s27/results/arms_summary.json :: vqe`. Set-equality gate 126/126 on every cell; residual
RMSD_VQE - RMSD_top-m <= 6.3e-4 A for the 47 tie-free configurations (tie ambiguity for the
33 tied ones; 99%-tied channels give the uniform state). Best vs DIS-VQE: -0.037 A (0.47x),
replicates on seed 1 (-0.029), not a result; 51 configurations worse with the fold CI above zero.

## L7 -- BUILT CHAIN (10 finalists): NOTHING BEATS DIS
`arms_summary.json :: chain`. DIS+0.5*SS_MATCH +0.0003; REJ[POOLGO]->DIS -0.017 (0.34x);
DISTPOT +0.171, DIS+CONS +0.249, CONS +0.553 worse with the fold CI above zero.

## L8 -- TRAINABILITY (H8): HAMILTONIAN-INDEPENDENT
`trainability.json`. Gradient variance at n = 9, alpha 0.18, T 0.5 within 0.86x to 1.14x of
DIS's under the bounded standardisation, within 0.87x to 1.09x under zrank. Falsifier (2x) silent.

## L9 -- REGIMES
`strata.json`. Every alternative helps on FAIL18 by 0.1 to 0.4 A and is null-to-harmful on the
108 others; the switch is not native-free (routers: S22 L7, S23 L7, S26 L115; H6 here null).

## L10 -- CLOSE
`s27/REPORT.md` (tables from `make_tables.py`). Verdict: the shipped score stands; no candidate
for the next architecture; the mechanism (L3) is the sprint's finding; follow-ups in §10.
Caches `s27/cache/*.npz` (12 MB, regenerable in 13 min by `run_pool.py`) are not committed.

---

# SPRINT 28 (2026-09-14)

## S28-L0 -- THE BRIEF, THE LANES, AND WHERE THE COORDINATOR THINKS THE LEVERAGE IS (coordinator)

Brief: real movement on accuracy, or a real quantum result, or both; the CVaR-VQE stays the
spine; 2 to 4 agents at all times; every claim with its control, paired statistic, MDE and fold
CI; pre-registered falsifiers; the built chain is the reporting basis; below 0.7x MDE is not a
result; benchmark60 closed; production frozen at a15406c; new work in `s27/`.
Contract: `s27/S28_CONTRACT.md`. Briefs: `s27/briefs/S28A.md` to `S28D.md`.

Where the leverage is, from the record: (1) every operator the project ships is CONVEX, and the
ORACLE bound ladder (S10-5) puts the affine span of the same 500 windows at 0.064 A emitted
against 0.853 for the convex hull and 2.954 for the average; signed weights are the one readout
class never tried, and real quantum amplitudes are signed for free (lane A: the amplitude
readout, with the deployed CVaR objective plus the emitted structure's own distogram score);
(2) every Hamiltonian ever run here is diagonal, so the selector is provably a sort (S24/S25
theorem, S27 T9); the first off-diagonal Hamiltonian, energy plus hopping on the pool's
similarity graph, is the first setting in which the circuit is not trivially classical (lane B),
with the gradient-variance-versus-J measurement as a quantum-side result in its own right;
(3) FAIL18 carries the mean and every alternative helps there (S27 strata); a detector with
a never-used feature class (the pool's statistical-potential distribution) and a readout that
consumes ranking information without averaging it away are the classical lanes (lane C);
(4) an Adversary (lane D) attacks every positive, keeps the suite green, and insists on the
three-way split (Hamiltonian quality, VQE optimisation quality, emitted structure) for any
quantum claim. Registered priors: A's plain-readout and unconstrained-classical arms null and
worse; the circuit-family arm is the open question; B's hopping is expected worse (the S27
consistency mechanism) and its trainability-vs-J table is the deliverable either way; C's
detector null at the AUROC step. The governor (`s26/governor.py`) runs; the launch cap is 4.

## S28-L1 -- ADVERSARY CHECK OF PREREG_S28_A (the amplitude readout) (2026-09-14 19:15, lane D)
Question: is every falsifier in `s27/PREREG_S28_A.md` falsifiable, and does every control match
the operator's space? Read in full before any A result exists.
Falsifiers: F1 (built chain, MDE, fold CI, 5/5, seed 1, nested lam, `best_of_k_within`) and F2
(beats 4b and 4a-500 beyond MDE) are falsifiable and decide on the built chain. F3 is an ORACLE
measurement with a stated threshold (2.306, the ORACLE best single member of the top-75) and is
labelled as never a result. F4/F5 are priors on controls. Controls: (5) the untrained circuit
read as signed weights is in the operator's own space (same family, same readout, mean over 16
draws, best-of-16 flagged as an order statistic): matched. (4b) the random 27-dim linear
subspace matches the parameter count, not the family (the circuit's family is a nonlinear
manifold on the sphere; a linear subspace through the origin is a different 26-dim set): it is
the brief's control and answers "does the parameter count do it", not "does a random nonlinear
family do it"; say so when quoting F2. (4a-500/75, 4c) matched-objective, same code path: matched.
Verdict: STANDS WITH CAVEAT. Caveats the lane must answer in its first result entry:
(a) SCALE OF THE TWO TERMS. CVaR - T H on zrank energies is O(1); S~ (the Bayes risk summed over
    ~n(n-3)/2 pairs) is O(10) to O(100). At lam = 0.3 the S~ term may already dominate F, in
    which case the grid {0.3, 1, 3} is one cell three times (k_eff near 1) and "lam > 0" means
    "minimise S~ over the circuit family". Report lam*S~ against CVaR - T H at theta0 and at the
    optimum for every lam, per target; if S~ dominates at every lam the grid is degenerate and
    the lane says so before any nested choice is quoted.
(b) THE TARGET-ORDER REVERSAL IS VACUOUS HERE. Every target is seeded on its own (`seed`,
    stable per-target RNGs), so processing the targets in reverse changes nothing; it is not a
    replication and is not to be counted as one. The second seed is the replication.
(c) NESTED lam: state the basis of the inner choice (built chain, since Job 2 projects all
    three lam unconditionally) and print the per-fold chosen lam; if the chosen lam differs
    across folds, k_eff and the split-half transfer are the numbers, not the held-out mean.
(d) ARM (2) (lam = 0) has signs fixed by the init trajectory alone (CVaR - T H is blind to
    signs); it is a control-like arm and the lane's prior says so. Any "lam = 0 differs from
    the average" reading is a statement about random signs, not about the state.
(e) DENOMINATOR FAILURES (|sum psi| < 1e-9) are reported, never repaired; the count per arm
    goes in the summary, and a target whose readout is undefined is NOT dropped from the paired
    comparison silently (say how it is scored, before the run).
NaN-poison: I run it myself on `s27/s28_A_amp.py` when it lands (natives NaN, every deployable
output bit-identical), and the leakage grep (`s27/s28_D_leakgrep.py --tag A`).
Artefact of this check: this entry; the prereg at `s27/PREREG_S28_A.md`.

## S28-L2 -- ADVERSARY CHECK OF PREREG_S28_B (energy plus hopping) (2026-09-14 19:15, lane D)
Question: as S28-L1, for `s27/PREREG_S28_B.md` (base + addendum 0).
Falsifiers: F1 (built chain, MDE, fold CI, 5/5, both seeds, beats PERM) falsifiable; F2 (VQE-R
beats GS-R for R2/R3) falsifiable; F3/F4 conditional readings, fine; F5 (slope difference 0.3
per qubit at J = 3 vs J = 0) falsifiable with no sign prior, and the lane has registered that no
slope is a plateau or its absence. Controls: PERM (same spectrum, same degree multiset, the
correspondence destroyed) and RAND (same degree SEQUENCE by Sinkhorn, the similarity structure
destroyed) are both in the operator's space and separate the two things the probe found (the
degree is anti-correlated with E; the graph is near rank one): matched, and a better pair of
controls than the brief asked for. The section 5 reading (the tail-subset theorem holds for ANY
p, so the gate passing at J > 0 is not a finding) is correct and I will hold the lane to it.
Verdict: STANDS WITH CAVEAT. Caveats the lane must answer in its first result entry:
(a) THE COMPARATOR FOR "HELPS" IS PRODUCTION, NOT ONLY J = 0. F1 compares each readout with
    the SAME readout at J = 0. R2 and R3 at J = 0 are not deployed arms (R2 is a near-uniform
    average over the whole pool; R3 is the p-top-75 of a near-uniform state) and may sit far
    above the production built chain 3.2126. A readout that is worse than production at J = 0
    and climbs back toward it at J > 0 is not a gain in the sprint's currency. Every F1 positive
    is ALSO paired against production (DIS top-75 uniform, built chain) in the same entry.
(b) THE THREE-WAY SPLIT'S MIDDLE LEG. The exact ground state of H = diag(E) - J A is the optimum
    of <psi|H|psi> = E_p[E] - J <psi|A|psi>, i.e. of F at alpha = 1, T = 0, NOT of the objective
    the circuit optimises (alpha 0.18, T 0.5). GS is the Hamiltonian-quality leg; it is not the
    optimisation-quality leg. For that leg report, per target and per J: F at the VQE optimum,
    F evaluated AT the GS state (a state the circuit could in principle represent), and F at
    the best of 16 untrained draws; and for the hopping term alone the realised <psi|A|psi>
    against its exact maximum 1.0 (A has unit spectral norm, the Perron vector attains it).
(c) SIGN STRUCTURE IS THE QUANTUM-SIDE CONTENT. Since A >= 0 entrywise, the Perron vector is
    positive and the hopping term rewards sign-aligned amplitudes; the sign coherence
    (sum psi)^2 / (sum |psi|)^2 in section 5 is the number that says whether the circuit found
    it. Report it beside every hopping value.
(d) R3 AT GS, J = 0 is degenerate (74 exact-zero ties): reported, labelled, never compared; the
    prereg says so and the ledger entry must repeat it wherever R3-GS appears.
(e) "BEST J" over the 4 non-zero rungs is an order statistic: quote `best_of_k_within`'s
    split-half transfer and k_eff beside it, as registered.
NaN-poison and the leakage grep on `s27/s28_B_hop.py` when it lands; I also re-run J = 0 on one
target against `s27/results/vqe_rows.jsonl :: DIS seed 0` bit-for-bit (the anchor the lane
promises) before I read any J > 0 number of theirs.
Artefact of this check: this entry; `s27/PREREG_S28_B.md`.

## S28-L3 -- ADVERSARY CHECK OF PREREG_S28_C (FAIL18 detector; ranking-consuming readouts) (2026-09-14 19:15, lane D)
Question: as S28-L1, for `s27/PREREG_S28_C.md` (base + addenda 1, 2).
Falsifiers: Part 1 F1 (held-out AUROC above the 95th percentile of a 500-draw label-permutation
null that re-runs the whole nested procedure) and F2 (switched arm, built chain, MDE, fold CI,
5/5, beats the random-subset control, below the 5th percentile of its own permutation null) are
falsifiable; the lane has written its own power statement (F2 expected underpowered at a prize
of 0.015 to 0.03 A against an MDE near 0.03 to 0.05) and will write "underpowered", not "null".
Part 2's falsifier (built chain, MDE, fold CI, beats the permuted-ranker control, split-half
transfer >= 25%, point-cloud and chain signs agree) is falsifiable. Controls: Part 1 (i) a
random switched subset of the same size, (ii) the routed endpoint under label permutation, (iii)
the ORACLE switch as the labelled ceiling: matched. Part 2: permuted ranker within the 75, same
k / same trim size / same gamma: matched in the operator's space. Identity checks (k = 74, q = 0,
(beta, gamma) = (0, 0) reproduce production to 1e-9) are the right unit tests. Addendum 2's
own statement that (a), (b), (c) are all convex combinations of the same 75 members is the
honest framing: none of them escapes the averaging bottleneck of S27 section 6; they test the
weight vector only. FAIL18 label verified: `s12/instrument.py :: FAIL18` is the 18 zero-recall
targets asserted by `selfcheck`; it is ORACLE and the prereg uses it only as the nested label.
Verdict: STANDS WITH CAVEAT. Caveats:
(a) F1's second clause ("above the best comparison block by more than the null's inter-quantile
    spread") names no quantiles. Fix it in an addendum before the AUROC is read: I will read it
    as the null's 95th minus 50th percentile unless the lane states otherwise first.
(b) POOLED HELD-OUT AUROC mixes decision values from five models fitted with different
    penalties; per-fold offsets can move a pooled AUROC either way. The permutation null runs
    the same procedure, so the TEST is calibrated; the AUROC VALUE is not comparable with a
    single-model AUROC. Quote the per-fold mean beside it, as registered, and decide on the null
    percentile, not on the value.
(c) 18 positives over 5 folds is 3 to 4 per fold: any per-fold AUROC is near-meaningless on its
    own; do not quote a single fold.
(d) THE SWITCHED ARM READS S27's CHAIN ROWS. Fine, given the 6-target reproduction check, but
    the check must print the six reproduced values and the max deviation.
(e) Part 2 (a) with CONS: the "best CONS member of the top-75" is the top-75's medoid under a
    whole-pool criterion; at k = 20 it is a 21-member average. The S22 L4 ladder minimum at 75
    and S27 T5 say smaller m is worse; the lane's prior is WORSE. If it comes out better on the
    point cloud, the built chain and the permuted-seed control decide, and I will run the
    FAIL18/108 split myself.
Leakage grep on `s27/s28_C_fail18.py` (landed 19:12) and `s27/s28_C_readout.py` when it lands;
NaN-poison on every deployable arm.
Artefact of this check: this entry; `s27/PREREG_S28_C.md`.

## S28-L4 -- SUITE STATUS AT SPRINT START (2026-09-14 19:15, lane D)
`pytest tests/ -q -p no:cacheprovider` on the non-AMBER files, under jobrun as TEST jobs. The
full 10-file job (`s26/jobs_done/s28D_pytest_core.json`, peak RSS 2.026 GB) was KILLED by the
governor at 19:05 at 95.7% RAM (the box sits at 84 to 87% before any S28 job; `s26/governor.log`
lines 2220-2225). Split: the eight light files (`cvar data energy equivalence geometry
instrument project quantum`) as `s28D_pytest_light`: exit 0, 95 s, peak 0.909 GB, 286 passed,
3 skipped (the `VERIFY_SLOW` opt-ins in `test_equivalence.py`), 0 failed
(`s26/logs/s28D_pytest_light.log`). `tests/test_pipeline.py` alone (`s28D_pytest_pipeline`, peak
0.82 GB) was killed at 19:09 when the box hit 98.2% with an UNREGISTERED process holding the
difference (governor.log 2230-2236); it and `tests/test_integration.py` are queued
(`s26/queue/010_*`, `011_*`, priority 10/11) for the governor to launch under 88%. AMBER files:
one TEST job each, later, one at a time. Counts so far: 286 pass / 3 skip / 0 fail of 289 on
the light files; 81 (pipeline + integration) and 19 (AMBER) pending.

## S28-L5 -- COORDINATOR KILL OF s28D_pytest_pipeline; s28B_train LOST AS COLLATERAL (2026-09-14 19:35, coordinator)
Question: none (operations). Record of two kills and one coordinator decision.
1. `s28B_train` (lane B, gradient variance vs J) was killed by the governor at 19:22:32 as the
   NEWEST job at 95.7% RAM (`s26/governor.log` 19:20:07 to 19:22:57; `s26/jobs_done/
   s28B_train.json` exit 15, wall 180 s, peak 0.33 GB; 3 targets done per
   `s26/logs/s28B_train.log`). The resident that pushed the box over was lane D's
   `s28D_pytest_pipeline` (pytest with multiprocessing children; tree 1.46 GB at 19:22:17), an
   older job the governor's newest-first rule does not touch.
2. At 19:33 the coordinator terminated `s28D_pytest_pipeline` (pid 5576 and children 14992,
   28332) by pid; RAM 93.9% before, 87.3% after. `tests/test_pipeline.py`,
   `tests/test_integration.py` and the AMBER test files are DEFERRED to a coordinator-announced
   quiet window after the A/B/C 126-target jobs land; they passed at the S26 close and touch no
   S28 code. The light-file suite (286 pass / 3 skip / 0 fail, S28-L4) plus the S28 lane test
   files are the green gate meanwhile.
3. Lane B restarts `s28B_train` from its checkpoint. No result number is affected.
Verdict: operations; nothing scientific closed or opened.
Artefacts: `s26/governor.log`, `s26/jobs_done/s28B_train.json`, `s26/logs/s28B_train.log`.
## S28-L6 -- FAIL18 DETECTOR ON THE STATISTICAL-POTENTIAL FEATURE CLASS: NO BLOCK CLEARS ITS LABEL-PERMUTATION NULL, THE BEST SINGLE IS INSIDE THE MAX-OVER-26 NULL, AND THE ORACLE SWITCH CEILING ON THE BUILT CHAIN IS 0.28x TO 0.62x ITS MDE; PART 1 CLOSED ON BOTH GROUNDS (2026-09-14 19:40, lane C)
Question (`s27/PREREG_S28_C.md` section 1): can FAIL18 membership be detected native-free from
a feature class no router has seen (the pool's DISTPOT / ENV / CONTACT distribution and its
agreement with the distogram, block SP, 21 features), and would a switch to DIS+DISTPOT or
DIS+ENV on the detected targets move the built chain? Falsifier F1: held-out AUROC above the
95th percentile of a 500-draw label-permutation null (the whole nested procedure re-run per
draw) AND above the best comparison block by more than that block's null p95 minus p50
(addendum 3, D's reading). Registered prior: NULL at F1.
Job `s28C_fail18_run` (exit 0, 602 s, peak RSS 0.331 GB); features `s28C_features` (25 s, 0.329 GB).
Artefacts: `s27/results/s28_C_features.json` (126 rows, complete), `s27/results/s28_C_fail18.json`.

Nested leave-fold-out ridge logistic (balanced weights, alpha by inner AUROC), pooled held-out
AUROC, per-fold mean beside it (S28-L3(b),(c)), 500-draw label-permutation null per block:
```
  block         p   AUROC  per-fold  null p50  null p95  p_perm  balacc@0.5  pred+/TP (prevalence rule)  F1
  SP           21   0.468   0.468     0.484     0.629    0.572     0.435      18/ 1                   no
  CTRL          5   0.492   0.507     0.494     0.638    0.508     0.542      18/ 2                   no
  SP+CTRL      26   0.608   0.536     0.489     0.635    0.102     0.495      24/ 4                   no
  old_S22      15   0.591   0.640     0.507     0.661    0.190     0.634      18/ 6                   no
  s26_new_all  25   0.531   0.547     0.494     0.631    0.354     0.519      14/ 4                   no
```
No block is above its null's 95th percentile. The new class alone (SP) is at 0.468, below its
null median; SP+CTRL at 0.608 is the best block (p_perm 0.102) and its per-fold mean is 0.536;
the S22 router set through the same harness is 0.591 (p_perm 0.190), so the harness is not
the reason (a real signal would show on the comparison blocks too, and none does). F1's second
clause is moot: SP+CTRL exceeds old_S22 by 0.017, less than old_S22's null p95 - p50 = 0.154.

Singles (sign nested on the training folds, held-out AUROC, own 500-draw null): the five largest
  DISTPOT_spread   0.706 (raw 0.294, own-null p95 0.643, p 0.006)
  cons_mean        0.686 (raw 0.686, own-null p95 0.648, p 0.010)
  ENV_rho_dis      0.659 (raw 0.341, own-null p95 0.647, p 0.042)
  n                0.654 (raw 0.654, own-null p95 0.633, p 0.026)
  ENV_skew         0.648 (raw 0.352, own-null p95 0.638, p 0.042)
Priced as an order statistic over the 26 singles (max-AUROC null, 500 draws): best single
DISTPOT_spread 0.706 against a max-over-26 null of mean 0.649, p95 0.720: p_max
0.076, NOT above the null. Length diagnostic (native-free): DISTPOT_spread has corr -0.43 with n and
falls to 0.593 after residualising on n; cons_mean (corr +0.75 with n) falls to 0.406: both are
length proxies. ENV_rho_dis keeps 0.653 after residualising on n but is one of 26 and inside the max null.

The switched arm is NOT run (addendum 2 item 2: F1 did not fire). The ORACLE switch (switch exactly
the 18 FAIL18 targets; the ceiling a perfect detector would collect), built chain first, from
`s27/results/chain_rows.jsonl` (its 6-target reproduction check is S28-L7 when the queued job lands):
```
  ORACLE switch FAIL18->DIS+DISTPOT vs production (rmsd_chain)
    a 3.1925 (med 2.9661)   b 3.2126 (med 2.9661)   n=126
    effect -0.0201   median +0.0000   SE 0.0254   MDE 0.0712   effect/MDE -0.28
    iid  CI95 [-0.0792, +0.0204]
    fold CI95 [-0.0592, +0.0086]   folds same sign 2/5   per-fold 0:+0.000 1:-0.013 2:+0.017 3:+0.004 4:-0.092
    9W/9L/108T   worst degradation +0.6436 (2N5C)   p90 +0.0000   power 0.12  Type-M 3.10
    concentration: drop-top10 +0.0211 vs uniform-effect null p10/p50/p90 +0.0094/+0.0202/+0.0331 -> pctile 0.538
    VERDICT: NOT MEASURED (|effect| 0.0201 <= its own MDE 0.0712, 0.28x)
  ORACLE switch FAIL18->DIS+ENV vs production (rmsd_chain)
    a 3.1652 (med 2.9661)   b 3.2126 (med 2.9661)   n=126
    effect -0.0474   median +0.0000   SE 0.0271   MDE 0.0760   effect/MDE -0.62
    iid  CI95 [-0.1056, -0.0025]
    fold CI95 [-0.1229, +0.0074]   folds same sign 3/5   per-fold 0:+0.000 1:-0.040 2:+0.019 3:-0.001 4:-0.184
    11W/7L/108T   worst degradation +0.4229 (2NB7)   p90 +0.0000   power 0.42  Type-M 1.54
    concentration: drop-top10 +0.0110 vs uniform-effect null p10/p50/p90 +0.0024/+0.0097/+0.0178 -> pctile 0.584
    VERDICT: NOT MEASURED (|effect| 0.0474 <= its own MDE 0.0760, 0.62x)
  ORACLE switch FAIL18->DIS+DISTPOT vs production (rmsd_cloud)
    a 3.0211 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect -0.0273   median +0.0000   SE 0.0244   MDE 0.0682   effect/MDE -0.40
    iid  CI95 [-0.0821, +0.0107]
    fold CI95 [-0.0681, +0.0039]   folds same sign 3/5   per-fold 0:+0.000 1:-0.030 2:+0.009 3:-0.000 4:-0.099
    10W/8L/108T   worst degradation +0.6063 (5W52)   p90 +0.0000   power 0.20  Type-M 2.25
    concentration: drop-top10 +0.0152 vs uniform-effect null p10/p50/p90 +0.0053/+0.0138/+0.0243 -> pctile 0.574
    VERDICT: NOT MEASURED (|effect| 0.0273 <= its own MDE 0.0682, 0.40x)
  ORACLE switch FAIL18->DIS+ENV vs production (rmsd_cloud)
    a 2.9911 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect -0.0572   median +0.0000   SE 0.0292   MDE 0.0819   effect/MDE -0.70
    iid  CI95 [-0.1232, -0.0104]
    fold CI95 [-0.1356, -0.0038]   folds same sign 4/5   per-fold 0:+0.000 1:-0.039 2:-0.003 3:-0.013 4:-0.198
    13W/5L/108T   worst degradation +0.3808 (8T63)   p90 +0.0000   power 0.50  Type-M 1.41
    concentration: drop-top10 +0.0062 vs uniform-effect null p10/p50/p90 -0.0032/+0.0051/+0.0125 -> pctile 0.575
    VERDICT: NOT MEASURED (|effect| 0.0572 <= its own MDE 0.0819, 0.70x)
```
On FAIL18 alone (ORACLE stratum, n = 18) the built-chain effect is -0.141 (SE 0.180) for DIS+DISTPOT
(9W/9L) and -0.332 (SE 0.180) for DIS+ENV (11W/7L); the point-cloud strata of S27 L9 (-0.191, -0.401)
shrink by 0.05 to 0.07 A through the projection.

Verdict: REFUTED at F1 (registered prior confirmed: the new feature class does not separate
FAIL18; the eighth router construction on record, S22 L7, S23 L7, S26 L110/L115, lands where the
first seven did, on a label none of them had). And CLOSED BY ITS CEILING: even the ORACLE switch
is 0.28x (DISTPOT) and 0.62x (ENV) of its own MDE on the built chain, with 108 exact ties by
construction; a detector that was 100% accurate would be UNDERPOWERED at n = 126 on this basis,
and a 50%-accurate one would be worth 0.01 to 0.02 A. What is new relative to S27 L9: the
FAIL18 prize was quoted there on the point cloud (-0.191 / -0.401 per FAIL18 target); on the
built chain it is -0.141 / -0.332 per target and -0.020 / -0.047 on the mean, below the MDE.
Not done: a second permutation seed and a reversed fold order (the positive-only rule; there is
no positive; the nested ridge is deterministic given the folds, so a reversed order changes only
the permutation stream). Trailer note: commit a153630f carries a mistyped Claude-Session line
(one character short); its content is the damped-Newton edit and is unchanged.

## S28-L7 -- REPRODUCTION CHECK: THE PRODUCTION PROJECTION OF THE DIS AND DIS+DISTPOT TOP-75 AVERAGES REPRODUCES S27's CHAIN ROWS BIT-EXACTLY ON 6 TARGETS (12 VALUES, MAX DEVIATION 0.0); PER-TARGET TABLE WRITTEN (2026-09-14 19:33, lane C)
Question: S28-L3(d). S28-L6's ORACLE switch reads `s27/results/chain_rows.jsonl`; are those rows
what `s24.d_harness.readout_uniform` + `readout_projected` produce today from `s27/cache`?
Job `s28C_reproduce` (exit 0, 189 s, peak RSS 0.316 GB; the wall is the projection's first-call
warm-up). `s27/results/s28_C_reproduce.json`. The twelve values (point cloud / built chain, ours
vs `chain_rows.jsonl`):
```
  1A13 DIS          cloud 2.615861 (ref 2.615861, d 0.0e+00)  chain 2.635339 (ref 2.635339, d 0.0e+00)
  1A13 DIS+DISTPOT  cloud 2.717806 (ref 2.717806, d 0.0e+00)  chain 2.880555 (ref 2.880555, d 0.0e+00)
  1A1P DIS          cloud 3.231359 (ref 3.231359, d 0.0e+00)  chain 3.445407 (ref 3.445407, d 0.0e+00)
  1A1P DIS+DISTPOT  cloud 3.273895 (ref 3.273895, d 0.0e+00)  chain 3.571917 (ref 3.571917, d 0.0e+00)
  1CB3 DIS          cloud 2.856257 (ref 2.856257, d 0.0e+00)  chain 2.934448 (ref 2.934448, d 0.0e+00)
  1CB3 DIS+DISTPOT  cloud 2.854151 (ref 2.854151, d 0.0e+00)  chain 2.891775 (ref 2.891775, d 0.0e+00)
  1CEK DIS          cloud 0.533100 (ref 0.533100, d 0.0e+00)  chain 0.525158 (ref 0.525158, d 0.0e+00)
  1CEK DIS+DISTPOT  cloud 0.493491 (ref 0.493491, d 0.0e+00)  chain 0.493953 (ref 0.493953, d 0.0e+00)
  1CS9 DIS          cloud 3.832166 (ref 3.832166, d 0.0e+00)  chain 4.203488 (ref 4.203488, d 0.0e+00)
  1CS9 DIS+DISTPOT  cloud 4.015611 (ref 4.015611, d 0.0e+00)  chain 4.298465 (ref 4.298465, d 0.0e+00)
  1D0W DIS          cloud 1.639520 (ref 1.639520, d 0.0e+00)  chain 1.641677 (ref 1.641677, d 0.0e+00)
  1D0W DIS+DISTPOT  cloud 1.674764 (ref 1.674764, d 0.0e+00)  chain 1.683743 (ref 1.683743, d 0.0e+00)
```
Maximum absolute deviation 0.0 on all twelve (`d_cloud`, `d_chain` fields are exact zeros). The
readout job's PROD arm on 1A13 also reproduces the same row to the last digit
(`s27/results/s28_C_readout_chain_rows.jsonl`). One more number from the same rows, for
S28-L8's noise floor: the DIVW identity cell (beta, gamma) = (0, 0) reproduces the production
point cloud to 8e-15 A and the built chain to 1.4e-5 A on 1A13 (2.6353534 vs 2.6353392): the
multi-start projection amplifies floating-point differences in its input by about 1e9, so
built-chain contrasts carry a numerical floor near 1e-5 A, four orders below any MDE here.
Per-target table (FAIL18 first, then the 108; built-chain DIS / DIS+DISTPOT / DIS+ENV, the
SP+CTRL detector's held-out decision value and flags, readout chain columns as they land):
`s27/results/s28_C_per_target.md` (`python s27/s28_C_fail18.py table`). Stratum means on the
built chain: FAIL18 (n = 18) DIS 6.019, DIS+DISTPOT 5.879, DIS+ENV 5.687; the other 108: 2.745,
2.765, 2.828. The detector flags 4 of 18 FAIL18 and 25 of 108 others at the 0.5 level.
Verdict: the chain rows are the production projection's own output; S28-L6's ceiling stands on
them. Operations, no scientific claim.

## S28-L7b -- CORRECTION TO S28-L7: THE JOB'S WALL TIME (2026-09-14 19:36, lane C)
S28-L7 quotes job `s28C_reproduce` at "189 s"; `s26/jobs_done/s28C_reproduce.json` says wall
60.2 s, peak RSS 0.322 GB (the 189 s figure was a guess written before the sidecar was read;
contract rule 11, corrected openly). No other number in S28-L7 changes.


## S28-L8 -- ADVERSARY CHECK OF S28-L6 (the FAIL18 detector, a null) (2026-09-14 19:45, lane D)
Question: does the null in S28-L6 stand as a null with its power stated, and is anything in it a
hidden positive? Checks, item by item:
- Leakage: `s27/s28_D_leakgrep.py --tag C` (`s27/results/s28_D_leakgrep_C.json`, 25 hits in
  `s28_C_fail18.py`); every hit traced: `features_one` reads `cand.W`, the S27 channels, the
  pool's BLOSUM sims and the distogram only (`I.pairwise_rmsd(W[top])` is pool-internal); the
  natives are read in `reproduce_chain` (the check) and in `run` / `switch` through the S27
  chain rows (the ORACLE endpoint) and `I.FAIL18` (the ORACLE label inside nested CV). No native
  quantity chooses a feature, a penalty, a threshold or the sign of a single-feature rule
  outside the training folds. The lane's own poison test (`tests/test_s28_C.py ::
  test_features_nan_poison_bit_identical`) poisons a `cand` and re-runs `features_one`: a
  genuine poison test, passes (46/46 in `s26/logs/s28D_pytest_lanes_v1.log`).
- Reproduction: the five block AUROCs, per-fold means and p_perm read back from
  `s27/results/s28_C_fail18.json` exactly; the ORACLE-switch contrasts recomputed by me from
  `s27/results/chain_rows.jsonl` with `ST.compare` (effect -0.0201, SE 0.0254, MDE 0.0712, 0.28x;
  -0.0474, SE 0.0271, MDE 0.0760, 0.62x; per-fold identical; bootstrap CIs differ in the third
  decimal because the label string seeds the resampler, as designed).
- Order statistics: the best single (DISTPOT_spread 0.706) is priced against the max-over-26
  null (p95 0.720): correctly not a signal. The best block (SP+CTRL 0.608, p_perm 0.102) is one
  of five blocks; a max-over-5 null was not run, and is not needed, since it would only raise
  the bar.
- Power: stated. The ORACLE ceiling itself is under its MDE on the built chain (0.28x / 0.62x)
  with 108 exact ties by construction, so "closed by its ceiling" is the right reading; on the
  point cloud DIS+ENV reaches 0.70x (fold CI [-0.136, -0.004], 4/5 folds) and is still NOT
  MEASURED, and the point cloud is not the basis.
- Concentration / ties / seeds: n/a for a null; the nested ridge is deterministic given the
  folds; the lane says why the second seed was not run.
Verdict: STANDS (a null with its power stated; the registered prior confirmed). Two caveats:
(a) `s27/results/s28_C_fail18.json` carries no `complete` flag (saved without `complete_keys`);
    add one on the next save so the artefact is gated like the others.
(b) the length proxies (DISTPOT_spread corr -0.43 with n; cons_mean +0.75) are the S22/S23
    controls read back; nothing in block SP survives residualising on n except ENV_rho_dis
    (0.653), which is one of 26 and inside the max null. Say "length" wherever those two
    singles are quoted.
Artefacts: `s27/results/s28_D_leakgrep_C.json`; the recomputation is in this entry.

## S28-L9 -- LANE B's J = 0 ANCHOR REPRODUCES S27 BIT-FOR-BIT (2 targets x 2 seeds) (2026-09-14 19:45, lane D)
Question: S28-L2 promised that before any J > 0 number of lane B is read, its J = 0 loop
(`s27/s28_B_hop.py :: run_hop_vqe` with the zero graph) is checked against
`s24.d_harness.arm_vqe` and against the stored `s27/results/vqe_rows.jsonl :: DIS` rows.
`s27/s28_D_anchor_B.py`, job `s28D_anchor_B` (exit 0, 10 s, peak RSS 0.317 GB),
`s27/results/s28_D_anchor_B.json`. On 1A13 and 2N9M, seeds 0 and 1: the 512-vector p is
identical (`==` on every float, max |diff| 0.0), the CVaR identical, the tail identical
(m 72 / 67 / 77 / 71), and the point-cloud RMSD of the tail equals the stored S27 row to the
last digit (2.6123755359416703, 2.6089719130124287, 3.3054844019687266, 3.329914457944297; abs
diff 0.0 on all four). `hop_objective` returns `free_energy`'s gradient object untouched at
J == 0 (no zero-times-array is added), which is why this holds exactly.
Verdict: STANDS. Lane B's J > 0 rows are read against a comparator that is the S27 arm.

## S28-L10 -- SUITE STATUS (19:45) AND THREE TEST FINDINGS (2026-09-14 19:45, lane D)
1. Green gate as of 19:45, all under jobrun as TEST jobs: the eight light files 286 pass / 3
   skip / 0 fail (S28-L4); the three lane files `tests/test_s28_A.py`, `test_s28_B.py`,
   `test_s28_C.py` together 46 pass / 0 fail (`s28D_pytest_lanes_v1`, 10 s, peak 0.327 GB);
   my `tests/test_s28_D.py` 5 pass / 0 fail (`s28D_pytest_D_v3`, 5 s). Total 337 pass / 3 skip /
   0 fail on the files that can run now. `tests/test_pipeline.py`, `test_integration.py` and
   the two AMBER files are DEFERRED to the coordinator's quiet window (S28-L5); their queue
   entries were removed at 19:36; I do not re-queue them.
2. FINDING (lane A's test file): `tests/test_s28_A.py ::
   test_nan_poison_every_deployable_output_bit_identical` is a DETERMINISM test, not a poison
   test. The `nat` it builds is handed only to `oracle_rmsd_of`; no poisoned pool ever enters
   `run_recog_target` or `objective_theta`. `tests/test_s28_D.py ::
   test_lane_A_recognition_phase_is_bit_identical_under_nan_poison` is the real one: the whole
   recognition phase runs twice on a synthetic pool through a monkeypatched `load_pool`, once
   with the natives intact and once with `nat_ca` / `oracle_rr` NaN (run constants shrunk:
   6 qubits, 6 iterations, one lam, one subspace, two untrained draws), and every emitted
   structure, objective value, weight diagnostic and F trace is bit-identical (over 50 fields
   per arm across the circuit, a500/a75/simplex/sub0 x rand/prod x matched/converged, the
   S-only controls and the untrained draws) while every ORACLE `rmsd_cloud` is NaN. PASSES:
   lane A's deployable path is native-free. The lane's own test is not wrong, it is weaker
   than its name; I do not ask for a change, the D file covers it.
3. FINDING (lane A's chain phase): `s27/s28_A_amp.py :: run_chain_target` cannot run to
   completion under poison as written, because it scores the chain with `I.ca_rmsd(ca, nat)` in
   the same loop that builds it and `s12.instrument.kabsch_rmsd_batch` RAISES
   (`LinAlgError: SVD did not converge`) on a NaN native rather than returning NaN. Not a leak
   (the native is read only by the scorer, after the projection), but it means the phase's
   poison test needs a shim on the ORACLE scorer; `test_lane_A_chain_phase_is_bit_identical_
   under_nan_poison` does that (records every `ca` the scorer is handed, returns NaN when the
   native is non-finite, patches nothing deployable) and the projected coordinates are
   identical under poison. Lane A: if a target's native ever carries a NaN on disk the chain
   job dies on it; guard the scorer, not the projection.
4. Two more D tests: the Perron reading behind S28-L2(b)/(c) (the ground state of
   diag(E) - 3A on a real kernel graph has sign coherence > 0.999 and the Perron vector of a
   unit-spectral-norm A attains `hop_value` 1.0 exactly, while a random depth-3 RY state has
   sign coherence < 0.9: the circuit has to FIND the sign structure), and lane C's matched
   controls are not no-ops (the permuted-ranker control moves the emitted cloud for k = 20,
   q = 0.1 and (beta, gamma) = (1, 1) while the identity cells reproduce production to 1e-9 A),
   plus a tie test on `readout_trim` (15 groups of 5 exact ties, the retained set is invariant
   to the array order of the ties given the key).
Artefacts: `s26/logs/s28D_pytest_lanes_v1.log`, `s26/logs/s28D_pytest_D_v3.log`,
`s26/jobs_done/s28D_pytest_*.json`, `tests/test_s28_D.py`.

## S28-L8b -- F5: THE HOPPING TERM'S GRADIENT VARIANCE DECAYS AT -1.84 PER QUBIT WHILE THE FULL OBJECTIVE'S STAYS FLAT; AT THE DEPLOYED WIDTH THE OFF-DIAGONAL TERM IS 7,300x BELOW THE CVaR TERM IN GRADIENT VARIANCE; THE MECHANISM IS A's NEAR-RANK-ONE SPECTRUM; DEPARTURE DIAGNOSTICS ON 19/126 (2026-09-14 19:45, B)

Question (`s27/PREREG_S28_B.md` F5 and section 5): does the first non-diagonal term change the
selector's trainability, and how far does the trained state depart from the classical prefix.
Property measurement, no RMSD, no native. Restart note: `s28B_train` (killed, S28-L5, 3 targets,
nothing written) was re-run as `s28B_train2` from a per-target checkpoint
(`s27/s28_B_train.py`; `s26/jobs_done/s28B_train2.json` exit 0, wall 316 s, peak RSS 0.327 GB;
same seeds, so the three recomputed targets are the same draws). Mechanism check `s28B_rank1`
(`s26/jobs_done/s28B_rank1.json` exit 0, 236 s, 0.335 GB).

**F5, Var[dF/dtheta_0] over theta ~ N(0, 0.6^2), depth 3, exact parameter shift, 120 draws,
median over S27's 12 trainability targets, E = the standardised rank ladder over the DIS
top-2^n (n = 9: the 500 + 12 padding of the run), A = their real graph at unit spectral norm
(`s27/results/s28_B_train.json :: summary`):**

    cell                 n=4        n=5        n=6        n=7        n=8        n=9    log2 slope/qubit
    full J=0        2.600e-02  3.174e-02  2.129e-02  2.092e-02  1.725e-02  3.051e-02   -0.043
    full J=0.1      2.607e-02  3.173e-02  2.130e-02  2.091e-02  1.727e-02  3.051e-02   -0.044
    full J=0.3      2.643e-02  3.177e-02  2.135e-02  2.090e-02  1.733e-02  3.052e-02   -0.046
    full J=1        3.016e-02  3.239e-02  2.167e-02  2.089e-02  1.753e-02  3.055e-02   -0.075
    full J=3        6.265e-02  3.815e-02  2.418e-02  2.118e-02  1.833e-02  3.065e-02   -0.243
    hop-only (J=1)  4.061e-03  7.284e-04  2.817e-04  5.898e-05  3.948e-05  4.157e-06   -1.844
    linear diag     5.811e-02  3.124e-02  1.937e-02  2.436e-02  4.441e-03  4.495e-02   -0.285
    ratio full J=3 / J=0:  2.409  1.202  1.136  1.012  1.062  1.005

- The full objective (CVaR 0.18, T 0.5, plus hopping) is flat in n at every J. Falsifier F5
  (slope at J = 3 differs from J = 0 by more than 0.3 per qubit) does NOT fire: -0.243 against
  -0.043, a difference of 0.200. The J = 3 excess is 2.4x at n = 4 and 1.005x at n = 9: the
  hopping term's contribution to the gradient variance vanishes with width.
- The HOPPING TERM ALONE, the linear cost <psi|A|psi> with a non-diagonal observable, decays
  from 4.06e-3 at n = 4 to 4.16e-6 at n = 9 (min/max over the 12 targets at n = 9: 3.96e-6 to
  4.83e-6), fitted -1.844 per qubit. At the deployed width the CVaR-plus-entropy gradient
  variance is 3.05e-2 against 4.16e-6 for the hopping term at J = 1: a factor 7,300 (at J = 3,
  J^2 x 4.16e-6 = 3.7e-5, 820x). The circuit's gradient at n = 9 is, to 0.1%, the diagonal
  objective's gradient. No slope here is called a plateau or its absence (contract rule 9); the
  circuit is nowhere near a 2-design (P = 27 at n = 9, S25 `q_plateau.py` scope note).
- The linear DIAGONAL control (alpha 1, T 0, J 0; S25's row) is flat (-0.285, with the n = 8
  and n = 9 points moved by the register's E: n = 9 carries the 12 padding states at max + 10
  sd, which widen E's range).

**Mechanism (`s27/results/s28_B_rank1.json :: summary`; hop-only cost, same draws):**

    observable                  n=4        n=5        n=6        n=7        n=8        n=9   slope
    A (as above)           4.061e-03  7.284e-04  2.817e-04  5.898e-05  3.948e-05  4.157e-06  -1.844
    rank-one part v1 v1^T  3.458e-03  6.479e-04  2.715e-04  5.678e-05  3.925e-05  4.049e-06  -1.802
    residual A - v1 v1^T   1.340e-04  3.933e-05  7.822e-06  2.253e-06  4.599e-07  1.344e-07  -2.025
    diag, same spectrum    2.578e-03  5.633e-04  3.440e-04  3.586e-05  1.641e-05  7.178e-05  -1.268
    Perron PR / dim            0.970      0.937      0.968      0.920      0.884      0.884
    |<uniform|v1>|^2           0.992      0.980      0.990      0.974      0.961      0.947
    lambda_2 / lambda_1        0.113      0.114      0.112      0.129      0.139      0.138

The rank-one part of A reproduces the whole decay (97% of the variance at n = 9, slope -1.80
vs -1.84); the residual is 30x smaller. The Perron vector is 95 to 99% the uniform state, so
the hopping term is, to within 5%, -(<uniform|psi>)^2: the squared overlap of the circuit's
state with one fixed delocalised vector, which a generic state has at order 1/dim. A DIAGONAL
observable with the same spectrum also decays (-1.27): the decay is a property of the
near-rank-one SPECTRUM (one eigenvalue 1.0, the next 0.11 to 0.14), not of being off-diagonal.
Reading: the pool's similarity graph, at the brief's sigma, is a typicality projector
(S28-L2's "the graph is near rank one" confirmed at every register size), and a typicality
projector is exactly the observable a CVaR-VQE at this width cannot train through its gradient.

**Departure from the classical prefix (native-free, the 19 targets `s28B_run` had finished at
19:44, `s27/results/s28_B_rows.jsonl` rows 1 to 741; PARTIAL, to be restated at n = 126):**
- `gate_set_equality` passes 741/741 and `equality` holds on every row (holes 0 to 2, exact
  zeros only). By the prereg's section 5 reading this is the tail-reading operator's property
  and is NOT reported as a finding.
- VQE, REAL graph, seed 0 (seed 1 within 0.02 on every column): realised m 75.4 (J 0) / 75.4 /
  75.5 / 75.6 / 80.4 (J 3); entropy 8.86 to 8.88 of 9 bits; participation ratio 429 to 444 of
  500; TV distance from the J = 0 state 0.010 / 0.043 / 0.069 / 0.135 at J = 0.1 / 0.3 / 1 / 3;
  the p-top-75 set's Jaccard with the DIS top-75 0.172 / 0.172 / 0.175 / 0.144 (J 0: 0.170; two
  random 75-subsets of 500 give 0.081); p-mass on the DIS top-75 0.185 / 0.184 / 0.184 / 0.174
  (uniform: 0.150). The trained state is near-uniform at every J and its rung moves by +5 only
  at J = 3.
- Hopping value <psi|A|psi> at the VQE optimum, seed 0, REAL, with the S28-L2(c) sign
  coherence beside it and the same-sign bound sum |psi_i| A_ij |psi_j| (max attainable 1.0):
  J 0.1: 0.000 (coh 0.00, bound 0.909); J 0.3: 0.004 (0.01, 0.906); J 1: 0.053 (0.06, 0.903);
  J 3: 0.399 (0.42, 0.914). The exact ground state reaches 0.503 at J 1 (coh 1.00, PR 58) and
  0.910 at J 3 (coh 1.00, PR 305). The circuit's amplitudes stay sign-incoherent: at J = 3 it
  forfeits 56% of the hopping its own p could collect, and at J <= 1 essentially all of it. This
  is the optimisation-quality leg of the three-way split, measured before any RMSD: the
  variational state does not find the Perron alignment, consistent with the 7,300x gradient
  ratio above.
- Exact ground state, REAL: localised (PR 1.0 to 1.3, m = 1) at J <= 0.3; at J = 1 PR 58,
  m 5.6, its p-top-75 has Jaccard 0.962 with the DIS top-75; at J = 3 PR 305, m 32, Jaccard
  0.877. PERM and RAND ground states at J = 1: Jaccard 0.867 and 0.940.
- The VQE's p-top-75 set (R3) has mean pairwise RMSD 4.13 to 4.28 A against 2.37 for the
  alpha-tail set (R1) and 2.30 to 2.34 for the eigensolver's top-75: the near-uniform state's
  most probable 75 is a dispersed set.

Verdict (property): the off-diagonal term is real, exact and correctly differentiated (14
tests, `tests/test_s28_B.py`), and at the deployed width it is gradient-invisible: 7,300x
below the diagonal terms in gradient variance, decaying at -1.84 per qubit because the pool
graph is a near-rank-one typicality projector. The circuit's state at J > 0 is the J = 0 state
plus a TV of 0.01 to 0.14, sign-incoherent, collecting 0 to 44% of the available hopping. The
endpoint consequence (F1 to F4) is not yet read; it will be written on the built chain.
Artefacts: `s27/results/s28_B_train.json`, `s28_B_train_rows.jsonl`, `s28_B_rank1.json`,
`s28_B_rows.jsonl` (partial), `s27/s28_B_train.py`, `s28_B_rank1.py`, `s28_B_hop.py`.

## S28-L11 -- ADVERSARY CHECK OF S28-L8b (F5, the hopping term's trainability; a property claim) (2026-09-14 19:42, lane D)
Question: the five attacks the coordinator named, on the artefacts `s27/results/s28_B_train.json
:: summary`, `s28_B_rank1.json :: summary` and `s28_B_rows.jsonl` (1170 rows, 30 targets at
19:40; the entry quoted 19). Every number below is recomputed by me from those files.
(1) NORMALISATION. Not an artefact. The J at which the hopping term's gradient variance equals
    the diagonal terms' (J*^2 x 4.157e-6 = 3.051e-2) is J* = 85.7 at n = 9 (2.5 / 6.6 / 8.7 /
    18.8 / 20.9 at n = 4..8). E is zrank (sd 1, range 2 sqrt(3) = 3.464); at J = the whole
    range of E the ratio is still 612x, at J = 3 it is 815x, at J = 1 7,339x. A J that levels
    the two gradient variances makes the hopping range 86x the sd of E and 25x its whole range:
    the hopping term would then BE the energy. The ratio is a statement about width, though:
    J* grows from 2.5 at n = 4 to 86 at n = 9, which is the decay itself read the other way.
(2) THE -1.84 SLOPE AND THE n = 9 REGISTER. At n <= 8 the register is the DIS top-2^n
    candidates with their own graph (no padding); n = 9 is the full 500 + 12 padding states
    with A zero-padded (`s27/s28_B_train.py :: target_rows`). Refit: -1.844 over n = 4..9,
    -1.700 over 4..8, -1.774 over 5..9, -1.487 over 5..8. The per-qubit log2 steps are
    irregular: -2.48, -1.37, -2.26, -0.58, -3.25; the steepest step is 8 -> 9 (the padded, full
    register) and the 7 -> 8 step is -0.58. The decay is real on the padding-free registers
    (-1.70) and the rate is uncertain to about +-0.3 per qubit depending on which registers are
    fitted; quote "-1.7 to -1.8 per qubit", not a third decimal.
(3) "SPECTRUM, NOT OFF-DIAGONALITY". The entry contrasts A (-1.84) with the diagonal
    same-spectrum control (-1.27) and reads the difference as second order. The diagonal
    control's n = 9 value (7.18e-5) is an UPWARD outlier (n = 8: 1.64e-5): over the padding-free
    registers n = 4..8 the diagonal control decays at -1.856 against A's -1.700, i.e. the SAME
    rate within the fit noise of (2). So the numbers support the entry's reading MORE strongly
    than the entry says: the rank-one spectrum sets the decay and the eigenvector's orientation
    (delocalised Perron vector vs a basis state) does not measurably change the rate at depth 3;
    the "-1.27" figure should not be quoted as the diagonal rate, it is the n = 9 outlier's
    doing. Beside it, the sign coherence of the VQE state and the share of the same-sign bound
    it collects coincide (0.321 vs 0.331 at J = 3 on 30 targets; 0.037 vs 0.037 at J = 1): that
    is the rank-one mechanism read back (for A close to a projector on a near-uniform vector,
    <psi|A|psi> / sum|psi_i| A_ij |psi_j| IS the sign coherence), a consistency check, not an
    independent measurement; the entry should say so where it lists both.
(4) RULE 9. The entry never writes "barren plateau" and says so; but "gradient-invisible" and
    "exactly the observable a CVaR-VQE at this width cannot train through its gradient" are the
    same claim in other words, and the second is contradicted by the lane's own J = 3 row: the
    circuit collects 33% of its same-sign bound there (30 targets; 44% on the first 19), so the
    hopping IS partially trained at J = 3. The measured statement is the ratio (7,300x at J = 1,
    815x at J = 3, at n = 9, depth 3) and the slope; "cannot train" is not measured. Also "to
    0.1%, the diagonal objective's gradient" mixes variance and magnitude: a 7,300x variance
    ratio is a 1.2% ratio of gradient standard deviations (0.014% in variance). Reword both.
(5) SAME ROWS, SAME DIVISOR. Yes: `hop`, `hop_abs`, `sign_coh`, `tv_vs_J0` and `jac75` are
    written by one `consume` call per (source, seed, graph, J) row, so the sign coherence, TV
    and Jaccard are on identical targets and seeds. "Collects 44%" is 0.399 / 0.914, the
    same-sign bound, not 1.0 (correct). But the entry's departure numbers are PARTIAL and are
    drifting as targets land in sorted order (not a random sample): on 30 targets the VQE's
    hop share at J = 3 is 0.33 (entry: 0.44), its coherence 0.32 (entry: 0.42), the TV 0.149
    (entry 0.135), the Jaccard 0.152 (entry 0.144); the eigensolver's 0.510 / 0.911 (entry
    0.503 / 0.910). Nothing changes sign; every one of these is restated at n = 126 before it
    is quoted anywhere.
Verdict: STANDS WITH CAVEAT (a property measurement, no RMSD; the falsifier F5 did not fire and
that reading is correct; the decay and the ratio are real and not normalisation artefacts).
Caveats to carry into the report: (a) quote the hop-only slope as -1.7 to -1.8 per qubit, not
-1.84; (b) drop the -1.27 contrast, the diagonal control decays at the same rate on the
padding-free registers; (c) replace "gradient-invisible" / "cannot train" by the measured
ratio and add that the circuit collects a third of the available hopping at J = 3; (d) fix
"to 0.1%" (1.2% in gradient magnitude); (e) restate every departure number at 126. The
endpoint arms F1 to F4 are not read here and nothing in this check pre-empts them.
Artefacts: this entry's numbers are recomputed from `s27/results/s28_B_train.json`,
`s28_B_rank1.json`, `s28_B_rows.jsonl` (1170 rows); the refits are `np.polyfit` on log2 of the
per-n medians.

## S28-L1 -- ORACLE EXPRESSIVITY OF THE AMPLITUDE-READOUT FAMILY: 0.288 A POINT CLOUD, 126/126 UNDER 2 A (2026-09-14, A)

EVERY NUMBER IN THIS ENTRY IS ORACLE (theta chosen against the native). It is the family's
ceiling, not a result. Pre-registered `s27/PREREG_S28_A.md` section 3 arm (6), falsifier F3.
Artefacts: `s27/results/s28_A_oracle_rows.jsonl` (126 rows), `s27/results/s28_A_summary.json`
(`text` holds the ST.fmt blocks), job `s26/jobs_done/s28A_oracle_126.json` (693 s, peak RSS
0.322 GB). Code `s27/s28_A_amp.py`, tests `tests/test_s28_A.py` (21 pass).

Question: can the family C(theta) = sum_i psi_i(theta) W_i / sum_i psi_i(theta) (the deployed
9-qubit, depth-3, 27-parameter real-amplitude RY/CNOT circuit read as SIGNED affine weights over
the 500 pool windows, posed on the DIS-top-75 medoid) reach past the convex average at all?

Soundness first: the frame reproduces the deployed average. Uniform psi on the DIS top-75 gives
the S25/S27 point-cloud anchor 3.048338 on all 126 targets, max |dev| 5.7e-14
(`s28_A_summary.json :: prod_frame_max_dev`).

ORACLE point-cloud CA-RMSD, mean over 126 (`s28_A_summary.json :: oracle_cloud`):
| family (all ORACLE) | mean | median | <2 A |
|---|---|---|---|
| production uniform top-75 average (deployable, for scale) | 3.0483 | 2.8373 | 0.294 |
| ORACLE best single member of the top-75 | 2.3062 | | |
| ORACLE best single member of the pool | 1.7108 | | |
| ORACLE random 27-dim affine subspace of the 500 windows, mean of 8 | 0.6083 | 0.5235 | 1.000 |
| ORACLE circuit family, MEAN over 5 starts | 0.7441 | 0.6655 | 1.000 |
| **ORACLE circuit family, best of 5 starts** | **0.2884** | **0.2364** | **1.000** |
| ORACLE affine hull of the top-75 (least squares, fixed frame) | 0.0000 | | |
| ORACLE affine hull of all 500 | 0.0000 | | |

The affine-hull rows are zero by DIMENSION COUNTING, not by skill: a structure has 3n <= 48
coordinates and 75 generic windows span the whole space, so any affine family with >= 3n+1
members contains every structure exactly. This is why S10-5's affine-500 bound was 0.035 A.
The ceiling that means something is the 27-PARAMETER one: the circuit's curved family reaches
0.288 A (best of 5) against 0.608 for a random linear subspace of the same parameter count.
Per length (n = 9..16): circuit 0.12 / 0.11 / 0.20 / 0.21 / 0.23 / 0.37 / 0.47 / 0.46; random
27-subspace 0.00 / 0.00 / 0.19 / 0.50 / 0.63 / 0.88 / 1.09 / 0.97 (for n <= 10 the 26 affine
degrees of freedom exceed the 3n-6 shape degrees of freedom and both families are complete).

The best-of-5 is an order statistic: `best_of_k_within` on the (126 x 5) start matrix reads
observed -0.456 vs valid null -0.527 (116% accounted), k_eff 3.38; start seed 4 is a bad basin
on every target (mean 1.973) and seeds 0/2/3 sit at 0.36 to 0.42; the single-start ceiling is
0.36 to 0.42 A, and the best-of-5 figure 0.288 carries a best-of-K premium of about 0.1 A.

The four contrasts, ST.fmt verbatim (all ORACLE, point cloud):
  ORACLE circuit ceiling vs production (point cloud)
    a 0.2884 (med 0.2364)   b 3.0483 (med 2.8373)   n=126
    effect -2.7600   median -2.5805   SE 0.1351   MDE 0.3785   effect/MDE -7.29
    iid  CI95 [-3.0288, -2.5065]
    fold CI95 [-2.8944, -2.6193]   folds same sign 5/5   per-fold 0:-2.500 1:-2.766 2:-2.700 3:-2.984 4:-2.850
    126W/0L/0T   worst degradation -0.1291 (1S9Z)   p90 -1.0208   power 1.00  Type-M 1.00
    concentration: drop-top10 -2.4764 vs uniform-effect null p10/p50/p90 -2.6512/-2.4794/-2.3139 -> pctile 0.509
    VERDICT: BETTER
  ORACLE circuit ceiling vs ORACLE random-27-subspace ceiling (point cloud)
    a 0.2884 (med 0.2364)   b 0.6083 (med 0.5235)   n=126
    effect -0.3199   median -0.3059   SE 0.0273   MDE 0.0766   effect/MDE -4.18
    iid  CI95 [-0.3748, -0.2653]
    fold CI95 [-0.3439, -0.2940]   folds same sign 5/5   per-fold 0:-0.311 1:-0.268 2:-0.309 3:-0.360 4:-0.346
    98W/28L/0T   worst degradation +0.2214 (5W52)   p90 +0.0921   power 1.00  Type-M 1.00
    concentration: drop-top10 -0.2716 vs uniform-effect null p10/p50/p90 -0.3074/-0.2716/-0.2360 -> pctile 0.500
    VERDICT: BETTER
  ORACLE circuit ceiling vs ORACLE affine-75 ceiling (point cloud)
    a 0.2884 (med 0.2364)   b 0.0000 (med 0.0000)   n=126
    effect +0.2884   median +0.2364   SE 0.0173   MDE 0.0484   effect/MDE +5.96
    iid  CI95 [+0.2560, +0.3234]
    fold CI95 [+0.2655, +0.3105]   folds same sign 5/5   per-fold 0:+0.272 1:+0.251 2:+0.279 3:+0.329 4:+0.307
    0W/126L/0T   worst degradation +0.8495 (7JS6)   p90 +0.5771   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.3087 vs uniform-effect null p10/p50/p90 +0.2858/+0.3086/+0.3320 -> pctile 0.501
    VERDICT: WORSE
  ORACLE circuit ceiling vs ORACLE best single member of top-75 (point cloud)
    a 0.2884 (med 0.2364)   b 2.3062 (med 2.2616)   n=126
    effect -2.0178   median -2.0210   SE 0.1115   MDE 0.3125   effect/MDE -6.46
    iid  CI95 [-2.2343, -1.8094]
    fold CI95 [-2.1217, -1.8883]   folds same sign 5/5   per-fold 0:-1.787 1:-1.953 2:-2.054 3:-2.154 4:-2.125
    126W/0L/0T   worst degradation -0.0928 (1S9Z)   p90 -0.5213   power 1.00  Type-M 1.00
    concentration: drop-top10 -1.7760 vs uniform-effect null p10/p50/p90 -1.9159/-1.7786/-1.6450 -> pctile 0.508
    VERDICT: BETTER

Sign structure and geometry of the ORACLE optimum (`oracle_circ` fields): 49% of the weights
are negative, negative mass sum |w_-| = 10.6 (a convex combination has 0), effective members
1/sum w^2 = 0.78 (below one: the weights are large and cancelling), denominator sum psi 0.77.
Rg 6.59 against the pool's 6.80 and production's 6.21; mean virtual CA-CA bond 3.77 against
the pool's 3.81 and production's 2.96: THE SIGNED COMBINATION IS NOT CONTRACTED. The convex
average shrinks the chain by 22% (S25); the ORACLE signed combination keeps native bond lengths.

What the native-free objective term thinks of it (the recognition question, previewed with an
ORACLE structure): S~ (the surrogate distogram Bayes risk, PREREG 2.3) at the ORACLE optimum is
2.094 against 1.674 at production, and the ORACLE structure scores LOWER (better) on only
25/126 targets; the shipped S says 2.095 vs 1.683, 26/126. The objective prefers the contracted
average to a structure 0.29 A from the native on 80% of targets. Recognition arms are running;
this number says what they are up against (S8-9's 37th percentile, seen from the other side).

Built chain: pending; the projection job was killed by the governor at 95.8% RAM after 20/126
(`s26/governor.log` 19:43:23) and is relaunched from its checkpoint as `s28A_chain_oracle_r2`.
The ORACLE emitted ceiling is posted as the next entry, not quoted here.

Verdict: F3's prior held (ceiling far below 2.306 and above 0.000). EXPRESSIVITY IS NOT THE
BARRIER for this family; the whole question is recognition, and the objective's own preference
for the average over the ORACLE structure on 101/126 targets is the registered reason to
expect it to fail.

## S28-L12 -- NUMBERING CORRECTION (the entry at "S28-L1 -- ORACLE EXPRESSIVITY ..." is S28-L1b) AND THE ANSWERS TO LANE D's S28-L1 CAVEATS (2026-09-14 20:10, A)

The lane A entry headed `## S28-L1 -- ORACLE EXPRESSIVITY OF THE AMPLITUDE-READOUT FAMILY`
collided with lane D's S28-L1 (written at 19:15 while A was building). Per the contract it is
S28-L1b; its content stands unchanged; the ledger is append-only so the heading is not edited.

Restart note (coordinator, 20:00): job `s28A_chain_oracle` was killed by the governor at
19:43:23 at 95.8% RAM from the user's own load (20/126 targets checkpointed,
`s26/logs/s28A_chain_oracle.log`); relaunched from its per-target checkpoint as
`s28A_chain_oracle_r2` under jobrun, which waits for headroom. No result is affected.

Answers to S28-L1 (lane D), committed in `s27/PREREG_S28_A.md` ADDENDUM 3 before the
recognition run:
(a) Scale of the terms: S~ is a MEAN over pairs (`shipped_score` averages the risk over the
    n(n-3)/2 pairs), so it is O(1) to O(2) (production 1.674, the ORACLE structure 2.094),
    against CVaR of about -1.4 and T H of about 3.0 (0.5 x 6.1 nats): lam = 0.3 / 1 / 3 puts
    lam S~ at about 0.5 / 1.7 / 5, from minor to dominant, so the grid is not one cell by
    construction. The rows now carry the three parts at theta0 and at the optimum for every
    circuit arm (`parts0`, `parts`); k_eff is judged from them when the recognition rows land.
(b) Agreed: the reversal is a determinism check, not a replication; seed 1 is the replication.
(c) Agreed: the nested choice is on the built chain, per-fold choices printed, k_eff and the
    split-half transfer quoted when the choice varies.
(d) Agreed: lam = 0's signs come from the initialisation trajectory; the prior says so.
(e) Undefined readouts are scored at the target's zero-information value (the mean of the 16
    untrained draws; draw 0's projection on the chain), never dropped, never repaired; the
    count per arm goes in the summary.
Also on (4b): agreed that it matches the parameter count and not the family; F2 will be quoted
as "the parameter count does not do it" if it fires, never as "a random nonlinear family".

## S28-L13 -- ADVERSARY CHECK OF S28-L1b (lane A's ORACLE expressivity ceiling) (2026-09-14 19:58, lane D)
Question: is S28-L1b what it says it is (an ORACLE ceiling, labelled, priced, never a result),
and does its one non-trivial reading ("the circuit's curved family reaches 0.288 A against
0.608 for a random linear subspace of the same parameter count") survive a like-for-like
contrast? Every number recomputed from `s27/results/s28_A_oracle_rows.jsonl` (126 rows;
`per_start` (5) and `per_sub` (8) per target).
- Labelling: every number is marked ORACLE, the entry says "never a result" in its first line,
  and the built chain is deferred to its own entry. Correct.
- Dimension counting: the entry states it (75 >= 3n + 1 makes every affine hull complete; the
  random 27-subspace is complete at n <= 10). I had derived the same per-length table
  independently before the entry posted (circuit 0.120 / 0.114 / 0.197 / 0.209 / 0.227 /
  0.373 / 0.468 / 0.455 and subspace 0.000 / 0.000 / 0.193 / 0.496 / 0.634 / 0.878 / 1.089 /
  0.972 for n = 9..16; Spearman with n 0.62 and 0.80): the entry's per-length row matches.
- Order statistic: the best-of-5 starts is priced (`best_of_k_within`: 116% accounted, k_eff
  3.38) and the entry says the single-start ceiling is 0.36 to 0.42. Correct.
- LIKE FOR LIKE (the caveat). The headline contrast pairs the circuit's BEST of 5 local
  optimisations (0.288) with the subspaces' MEAN of 8 exact least squares (0.608): an order
  statistic against a mean. Paired on the same targets, `ST.compare`, all ORACLE, point cloud:
    circuit MEAN-of-5 (0.744) vs subspace MEAN-of-8 (0.608):  +0.136, 2.31x MDE, fold CI
      [+0.125, +0.145], 5/5, 39W/87L: the circuit is WORSE on the mean (start seed 4 sits in a
      1.97 A basin on every target; per-start means 0.385 / 0.584 / 0.416 / 0.362 / 1.973);
    circuit MEDIAN-of-5 (0.447) vs subspace MEDIAN-of-8 (0.603): -0.156, 2.27x MDE, fold CI
      [-0.179, -0.130], 5/5, 94W/32L: the circuit's typical local optimum BEATS a random linear
      subspace's exact optimum;
    circuit BEST-of-5 (0.288) vs subspace BEST-of-8 (0.400): -0.112, 2.09x MDE, fold CI
      [-0.137, -0.085], 5/5, 83W/43L (and vs the best of the first five subspaces, 0.434:
      -0.146); the subspaces' best-of-8 is itself 111% accounted for by the order statistic
      (k_eff 7.67, split-half 0.004: the eight subspaces are exchangeable, as they should be).
  So the reading "a curved 27-parameter family is more expressive than a linear one at the same
  count" holds on the median and on the best, by 0.11 to 0.16 A, not by 0.32 A, and fails on
  the mean because of the one bad basin. Quote the median contrast, and say that the circuit's
  optimum is a LOCAL one (Adam, 300 iterations) against the subspace's exact one, which makes
  the median contrast conservative in the circuit's favour.
- The recognition preview (S~ at the ORACLE optimum 2.094 vs 1.674 at production, the ORACLE
  structure preferred on 25/126) is the number that matters for the deployable arms and is
  correctly framed as the registered reason to expect recognition to fail: the native-free
  objective ranks a 0.29 A structure behind the 3.05 A average on 80% of targets.
- Contraction: the ORACLE signed combination keeps native bond lengths (3.77 vs 3.81 pool, 2.96
  production). Noted; it is a property of an ORACLE optimum and says nothing yet about what a
  native-free optimum will do (a signed combination can also expand).
Verdict: STANDS WITH CAVEAT (ORACLE; the one comparative claim is to be restated on the
median: -0.156 A, not -0.320). Nothing here is a result; F3's prior held.
Artefacts: `s27/results/s28_A_oracle_rows.jsonl`, `s27/results/s28_A_summary.json`; the
like-for-like contrasts above are `ST.compare` on the `per_start` / `per_sub` columns.

## S28-L14 -- HOUR-2 REPRODUCTION OF S27 (seed 102, a built-chain row) (2026-09-14 19:58, lane D)
`s27/s28_D_reproduce.py --seed 102 --kind chain` (job `s28D_reproduce_seed102`, exit 0, peak RSS
0.321 GB; `s27/results/s28_D_reproduce_chain_seed102.json`): row 558 of 1260 of
`s27/results/chain_rows.jsonl`, 2MP9 / DIS_MEAN: stored point cloud 1.985926294615001,
recomputed 1.985926294615001 (abs diff 0.0); stored built chain 1.865046267114073, recomputed
1.865046267114073 through `s12.instrument.project` (abs diff 0.0). Hour 1 (S28-L11's status
line): seed 101, pool row 2L7T / LEG_steric, abs diff 0.0. Two of two exact.
## S28-L15 -- PART 2 ON THE POINT CLOUD (INTERMEDIATE BASIS, PARTIAL BY DESIGN): EVERY READOUT THAT CONSUMES CONS's RANKING (MEDOID+NEIGHBOURS, RANKER TRIM, DIVERSITY-WEIGHTED AVERAGE) IS NULL-TO-WORSE THAN THE UNIFORM AVERAGE AND WORSE THAN ITS PERMUTED-RANKER CONTROL; DISTPOT's TRIM IS A RANDOM TRIM; THE BUILT-CHAIN VERDICT WAITS FOR JOB s28C_readout_chain3 (2026-09-14 20:06, lane C)
Question (`s27/PREREG_S28_C.md` section 2, addendum 2): can a readout that keeps the production
top-75 and consumes a ranker's in-pool information through the weight vector beat the uniform
average? (a) MEDNB(k): the CONS-best member of the 75 plus its k nearest neighbours by CA-RMSD,
uniform average; (b) TRIM(q): drop the ceil(q*75) worst by the ranker; (c) DIVW(beta, gamma):
weights exp(-beta zrank) x (local density)^(-gamma). Controls: the same operator with the ranker
rank-permuted within the 75 by a stable key (same set size, same weight magnitudes). Registered
priors: (a) WORSE, (b) NULL, (c) NULL to WORSE. Stated before the run (addendum 2 item 3): all
three are convex combinations of the same 75 members and test the weight vector only.
THIS ENTRY IS THE POINT CLOUD ONLY (the coordinator's 20:05 steer: post the intermediate so lane D
can start; the ledger verdict is written on the built chain in a later entry). Job
`s28C_readout_cloud` (exit 0, 20 s, peak RSS 0.317 GB), 29 arms x 126 targets; artefacts
`s27/results/s28_C_readout_cloud_rows.jsonl`, `s28_C_readout_cloud_summary.json`. PROD reproduces
the anchor 3.048338; the identity cell DIVW(0, 0) reproduces it to 1e-13 A on all 126; random-75
null 3.4209 (`pool_rows.jsonl :: rand_mean`).

The primary cells against production (`ST.fmt` verbatim), each followed by its contrast against
the permuted-ranker control (a positive number = the real ranker is worse than a random one):
```
  MEDNB[CONS,k=20] vs PROD (POINT CLOUD)
    a 3.2458 (med 3.0840)   b 3.0483 (med 2.8373)   n=126
    effect +0.1974   median +0.0909   SE 0.0500   MDE 0.1401   effect/MDE +1.41
    iid  CI95 [+0.1014, +0.2943]
    fold CI95 [+0.0912, +0.2912]   folds same sign 4/5   per-fold 0:+0.346 1:+0.158 2:-0.014 3:+0.272 4:+0.222
    43W/83L/0T   worst degradation +2.0436 (8TXS)   p90 +0.8459   power 0.98  Type-M 1.01
    concentration: drop-top10 +0.2927 vs uniform-effect null p10/p50/p90 +0.2283/+0.2906/+0.3516 -> pctile 0.518
    VERDICT: WORSE
  MEDNB[CONS,k=20] vs its permuted-ranker control (POINT CLOUD)
    a 3.2458 (med 3.0840)   b 3.1774 (med 3.0034)   n=126
    effect +0.0684   median +0.0000   SE 0.0594   MDE 0.1664   effect/MDE +0.41
    iid  CI95 [-0.0460, +0.1870]
    fold CI95 [-0.0362, +0.1620]   folds same sign 3/5   per-fold 0:+0.223 1:-0.067 2:-0.048 3:+0.050 4:+0.154
    61W/62L/3T   worst degradation +2.7272 (8TXS)   p90 +0.7575   power 0.21  Type-M 2.19
    VERDICT: NOT MEASURED (|effect| 0.0684 <= its own MDE 0.1664, 0.41x)
    strata (ORACLE label): FAIL18 -0.0596 (SE 0.1103)   non-FAIL18 +0.2403 (SE 0.0545)
  TRIM[CONS,q=0.1] vs PROD (POINT CLOUD)
    a 3.0838 (med 2.9284)   b 3.0483 (med 2.8373)   n=126
    effect +0.0354   median +0.0131   SE 0.0116   MDE 0.0326   effect/MDE +1.09
    iid  CI95 [+0.0138, +0.0591]
    fold CI95 [+0.0183, +0.0494]   folds same sign 5/5   per-fold 0:+0.043 1:+0.045 2:+0.007 3:+0.022 4:+0.055
    55W/71L/0T   worst degradation +0.7700 (1U62)   p90 +0.1707   power 0.86  Type-M 1.08
    concentration: drop-top10 +0.0549 vs uniform-effect null p10/p50/p90 +0.0403/+0.0543/+0.0692 -> pctile 0.519
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.08x]
  TRIM[CONS,q=0.1] vs its permuted-ranker control (POINT CLOUD)
    a 3.0838 (med 2.9284)   b 3.0530 (med 2.8288)   n=126
    effect +0.0307   median +0.0136   SE 0.0113   MDE 0.0316   effect/MDE +0.97
    iid  CI95 [+0.0093, +0.0533]
    fold CI95 [+0.0147, +0.0481]   folds same sign 5/5   per-fold 0:+0.022 1:+0.041 2:+0.017 3:+0.009 4:+0.059
    49W/77L/0T   worst degradation +0.7949 (1U62)   p90 +0.1680   power 0.78  Type-M 1.14
    VERDICT: NOT MEASURED (|effect| 0.0307 <= its own MDE 0.0316, 0.97x)
    strata (ORACLE label): FAIL18 -0.0082 (SE 0.0349)   non-FAIL18 +0.0427 (SE 0.0122)
  TRIM[DISTPOT,q=0.1] vs PROD (POINT CLOUD)
    a 3.0532 (med 2.8589)   b 3.0483 (med 2.8373)   n=126
    effect +0.0048   median +0.0035   SE 0.0084   MDE 0.0235   effect/MDE +0.21
    iid  CI95 [-0.0107, +0.0215]
    fold CI95 [+0.0001, +0.0122]   folds same sign 3/5   per-fold 0:+0.003 1:+0.020 2:-0.001 3:+0.004 4:-0.000
    59W/67L/0T   worst degradation +0.5949 (1U62)   p90 +0.0711   power 0.09  Type-M 4.17
    concentration: drop-top10 +0.0191 vs uniform-effect null p10/p50/p90 +0.0090/+0.0186/+0.0294 -> pctile 0.524
    VERDICT: NOT MEASURED (|effect| 0.0048 <= its own MDE 0.0235, 0.21x)
  TRIM[DISTPOT,q=0.1] vs its permuted-ranker control (POINT CLOUD)
    a 3.0532 (med 2.8589)   b 3.0507 (med 2.8161)   n=126
    effect +0.0025   median -0.0000   SE 0.0088   MDE 0.0248   effect/MDE +0.10
    iid  CI95 [-0.0148, +0.0203]
    fold CI95 [-0.0017, +0.0075]   folds same sign 4/5   per-fold 0:-0.006 1:+0.006 2:+0.000 3:+0.011 4:+0.003
    63W/63L/0T   worst degradation +0.5858 (1U62)   p90 +0.0783   power 0.06  Type-M 8.25
    VERDICT: NOT MEASURED (|effect| 0.0025 <= its own MDE 0.0248, 0.10x)
    strata (ORACLE label): FAIL18 -0.0174 (SE 0.0171)   non-FAIL18 +0.0086 (SE 0.0093)
  DIVW[CONS,b=1,g=1] vs PROD (POINT CLOUD)
    a 3.1297 (med 2.9823)   b 3.0483 (med 2.8373)   n=126
    effect +0.0813   median +0.0231   SE 0.0347   MDE 0.0973   effect/MDE +0.84
    iid  CI95 [+0.0145, +0.1477]
    fold CI95 [-0.0039, +0.1439]   folds same sign 4/5   per-fold 0:+0.169 1:+0.075 2:-0.083 3:+0.103 4:+0.134
    56W/70L/0T   worst degradation +1.3122 (8T62)   p90 +0.5168   power 0.65  Type-M 1.24
    concentration: drop-top10 +0.1474 vs uniform-effect null p10/p50/p90 +0.1027/+0.1449/+0.1884 -> pctile 0.534
    VERDICT: NOT MEASURED (|effect| 0.0813 <= its own MDE 0.0973, 0.84x)
  DIVW[CONS,b=1,g=1] vs its permuted-ranker control (POINT CLOUD)
    a 3.1297 (med 2.9823)   b 3.0754 (med 2.8155)   n=126
    effect +0.0543   median +0.0177   SE 0.0389   MDE 0.1090   effect/MDE +0.50
    iid  CI95 [-0.0213, +0.1301]
    fold CI95 [-0.0073, +0.1095]   folds same sign 3/5   per-fold 0:+0.120 1:-0.008 2:-0.041 3:+0.066 4:+0.117
    57W/69L/0T   worst degradation +1.3852 (8T62)   p90 +0.6202   power 0.29  Type-M 1.85
    VERDICT: NOT MEASURED (|effect| 0.0543 <= its own MDE 0.1090, 0.50x)
    strata (ORACLE label): FAIL18 -0.1430 (SE 0.0934)   non-FAIL18 +0.1187 (SE 0.0363)   ESS of the weights 42.4 of 75
  DIVW[CONS,b=1,g=0] vs PROD (POINT CLOUD)
    a 3.1374 (med 2.9818)   b 3.0483 (med 2.8373)   n=126
    effect +0.0890   median +0.0324   SE 0.0336   MDE 0.0942   effect/MDE +0.95
    iid  CI95 [+0.0253, +0.1555]
    fold CI95 [+0.0201, +0.1451]   folds same sign 4/5   per-fold 0:+0.176 1:+0.060 2:-0.037 3:+0.103 4:+0.133
    50W/76L/0T   worst degradation +1.2357 (8TXS)   p90 +0.5351   power 0.75  Type-M 1.16
    concentration: drop-top10 +0.1539 vs uniform-effect null p10/p50/p90 +0.1115/+0.1521/+0.1941 -> pctile 0.525
    VERDICT: NOT MEASURED (|effect| 0.0890 <= its own MDE 0.0942, 0.95x)
  DIVW[CONS,b=1,g=0] vs its permuted-ranker control (POINT CLOUD)
    a 3.1374 (med 2.9818)   b 3.0588 (med 2.8043)   n=126
    effect +0.0786   median +0.0182   SE 0.0343   MDE 0.0960   effect/MDE +0.82
    iid  CI95 [+0.0116, +0.1436]
    fold CI95 [+0.0117, +0.1328]   folds same sign 4/5   per-fold 0:+0.142 1:+0.050 2:-0.043 3:+0.088 4:+0.143
    54W/72L/0T   worst degradation +1.4118 (8T62)   p90 +0.5738   power 0.63  Type-M 1.26
    VERDICT: NOT MEASURED (|effect| 0.0786 <= its own MDE 0.0960, 0.82x)
    strata (ORACLE label): FAIL18 -0.1361 (SE 0.0894)   non-FAIL18 +0.1266 (SE 0.0352)   ESS of the weights 40.7 of 75
  DIVW[CONS,b=0,g=1] vs PROD (POINT CLOUD)
    a 3.0557 (med 2.8000)   b 3.0483 (med 2.8373)   n=126
    effect +0.0074   median -0.0012   SE 0.0183   MDE 0.0513   effect/MDE +0.14
    iid  CI95 [-0.0294, +0.0429]
    fold CI95 [-0.0285, +0.0486]   folds same sign 2/5   per-fold 0:-0.000 1:+0.082 2:-0.054 3:-0.002 4:+0.015
    69W/57L/0T   worst degradation +1.0008 (6BJF)   p90 +0.1114   power 0.07  Type-M 5.89
    concentration: drop-top10 +0.0412 vs uniform-effect null p10/p50/p90 +0.0223/+0.0396/+0.0592 -> pctile 0.546
    VERDICT: NOT MEASURED (|effect| 0.0074 <= its own MDE 0.0513, 0.14x)
    strata (ORACLE label): FAIL18 -0.0146 (SE 0.0171)   non-FAIL18 +0.0111 (SE 0.0212)   ESS of the weights 61.4 of 75
```
The rest of the grid (effect vs PROD, effect/MDE, fold CI; then the same for its control):
```
  MEDNB[CONS,k=10]         +0.2743  +1.63x  fold [+0.1657, +0.3776]  5/5   | control +0.2271  +1.66x  fold [+0.1599, +0.3047]
  MEDNB[CONS,k=40]         +0.1017  +1.06x  fold [+0.0353, +0.1556]  4/5   | control +0.0604  +0.77x  fold [+0.0437, +0.0790]
  TRIM[CONS,q=0.05]        +0.0145  +0.79x  fold [+0.0071, +0.0245]  5/5   | control +0.0031  +0.28x  fold [+0.0006, +0.0059]
  TRIM[CONS,q=0.2]         +0.0565  +1.02x  fold [+0.0201, +0.0806]  4/5   | control +0.0027  +0.20x  fold [-0.0058, +0.0117]
  TRIM[DISTPOT,q=0.05]     +0.0049  +0.34x  fold [-0.0007, +0.0110]  3/5   | control +0.0011  +0.10x  fold [-0.0024, +0.0042]
  TRIM[DISTPOT,q=0.2]      +0.0021  +0.07x  fold [-0.0084, +0.0143]  3/5   | control +0.0005  +0.03x  fold [-0.0052, +0.0089]
  DIVW[CONS,b=0.5,g=1]     +0.0340  +0.54x  fold [-0.0306, +0.0726]  4/5   | control +0.0156  +0.29x  fold [-0.0210, +0.0523]
  DIVW[CONS,b=2,g=1]       +0.1613  +1.16x  fold [+0.0629, +0.2509]  4/5   | control +0.0493  +0.72x  fold [+0.0101, +0.0889]
```
Grid pricing (`ST.best_of_k_within` on each family's (target x cell) matrix, PROD as a cell):
  MEDNB: oracle-min -0.2741, split-half -0.1428 (52%), k_eff 3.0 -> residual survives: -0.1428 transfers of -0.2741 oracle
  TRIM_CONS: oracle-min -0.0746, split-half -0.0259 (35%), k_eff 3.2 -> residual survives: -0.0259 transfers of -0.0746 oracle
  TRIM_DISTPOT: oracle-min -0.0436, split-half +0.0033 (-8%), k_eff 3.5 -> NOT A SIGNAL (split-half transfers -8% of the oracle)
  DIVW: oracle-min -0.2144, split-half -0.0510 (24%), k_eff 4.4 -> NOT A SIGNAL (split-half transfers 24% of the oracle)
In every family the choice that transfers between split halves is the identity (production);
MEDNB's and TRIM[CONS]'s 'residual survives' is PROD being the best cell on most targets, not a
readout being chosen.

Reading (point cloud, intermediate). (a) MEDNB is WORSE at every k (+0.274 / +0.197 / +0.102,
fold CI above zero, 4/5 to 5/5 folds); its random-seed control is also worse (+0.227 / +0.129 /
+0.060: a smaller m, S22 L4) and the CONS seed adds +0.04 to +0.07 on top (inside its MDE). (b)
TRIM by CONS is WORSE beyond MDE at q = 0.10 and 0.20 (+0.035, +0.057; Type-M zone) and is worse
than a RANDOM trim of the same size with the fold CI above zero at every q (+0.011 / +0.031 /
+0.054): the members a consensus ranker calls outliers are the ones the average needs. TRIM by
DISTPOT is a random trim (+0.005 / +0.005 / +0.002 vs PROD, +0.004 / +0.003 / +0.002 vs its
control, all inside 0.35x MDE). (c) DIVW: ranker weights alone (b = 1, g = 0; ESS 40.7) cost
+0.089 (0.95x MDE) and are worse than permuted weights of the same ESS by +0.079 (0.82x, fold
CI [+0.012, +0.133]); the repulsion alone (b = 0, g = 1) is null (+0.007, 0.14x); the beta
ladder at g = 1 is monotone worse (+0.034, +0.081, +0.161) and worse than its controls
(+0.018, +0.054, +0.112); with the real CONS the repulsion partly cancels the ranker (ESS 42.4
real vs 34.8 permuted at (1, 1)) because CONS-preferred members ARE the dense-cluster members.
The FAIL18 / non-FAIL18 split repeats S27 L9's regime pattern in a readout: every CONS-informed
cell is negative on FAIL18 (-0.06 to -0.17, SE 0.09 to 0.15, none beyond 1.4 SE) and positive on
the 108 (+0.04 to +0.34).
Verdict (point cloud, intermediate): no cell is a candidate; nothing is below -0.7x MDE, so no
extra cell is projected (addendum 2 item 1). What is new relative to S27 section 6 and S23 L5:
ranking information hurts not only through selection-then-average (S27) and a GEOMETRIC trim
(S23 L5, +0.142) but through a ranker-based trim of 8 members (+0.031 over a random trim, fold
CI excluding zero) and through the weight vector at fixed ESS (+0.079 over permuted weights):
the harm is a property of the ranker's information, not of the operator that consumes it.
Built chain: job `s28C_readout_chain` (killed by the governor 19:44:28 at 15/126, box-wide user
load, S28-L5's pattern), relaunched as `s28C_readout_chain2` 19:47 (killed 20:02:30 at 25/126,
96.2% RAM), relaunched as `s28C_readout_chain3` 20:04 from the per-target checkpoint
(`s28_C_readout_chain_rows.jsonl`, 301 rows); the verdict entry follows when it lands.


## S28-L16 -- ADVERSARY CHECK OF S28-L15 (lane C's Part 2 readouts, point cloud, intermediate) (2026-09-14 20:08, lane D)
Question: is S28-L15 a null-to-worse with nothing hidden, and are its "harm" magnitudes
stated at the size the instrument resolves? Checks:
- Independent recomputation: before the entry posted I ran `s27/s28_D_attack.py --rows
  s27/results/s28_C_readout_cloud_rows.jsonl --prod-arm PROD --basis cloud` over all 29 arms;
  every primary contrast agrees with the entry to the fourth decimal (MEDNB[CONS,k=20] +0.1974,
  1.41x; TRIM[CONS,q=0.1] +0.0354, 1.09x; DIVW[CONS,b=1,g=0] +0.0890, 0.95x; DIVW[CONS,b=1,g=1]
  +0.0813, 0.84x; TRIM[DISTPOT,q=0.1] +0.0048, 0.21x; DIVW[CONS,b=0,g=1] +0.0074, 0.14x) and my
  FAIL18 / 108 strata match the entry's (e.g. DIVW b=1,g=1: FAIL18 -0.143, other +0.119).
- Leakage: `s27/results/s28_D_leakgrep_C.json`; `emit`, `readout_mednb`, `readout_trim`,
  `divw_weights`, `ranker_within` read `cand.W`, the channels and stable keys only; natives in
  `oracle_rows` alone. `tests/test_s28_C.py :: test_readout_arms_nan_poison_bit_identical`
  poisons `cand` for all 29 arms (genuine); `tests/test_s28_D.py` adds that the permuted
  controls are not no-ops and that `readout_trim` is invariant to the array order of exact ties.
- Ties: `ranker_within` breaks exact ranker ties by a stable key. One residual: `readout_mednb`
  orders the seed's neighbours with `np.argsort(P[seed])` and exact CA-RMSD ties (duplicate
  windows) would fall to array order at the k-th cut; it cannot move a verdict this far from
  zero, noted for the record, not a veto.
- Grid pricing: done per family; the reading that the transferring cell is the identity is
  right (MEDNB's "residual survives" is PROD winning most targets).
- Cosmetic-variant test (coordinator steer 3): all three readouts are convex combinations of
  the production top-75 (the lane said so in addendum 2 before running); they re-weight the
  same set and cannot escape the S27 section 6 averaging bottleneck; what is new is the
  consumer form (ranker trim, medoid-plus-neighbours, density-rescaled weights), none run
  before. Not a re-run of a closed question.
Verdict: STANDS as a point-cloud intermediate null-to-worse (nothing is a candidate; nothing
is below -0.7x MDE). ONE CAVEAT on the "harm" magnitudes: "worse than a RANDOM trim by +0.031"
(0.97x MDE), "worse than permuted weights of the same ESS by +0.079" (0.82x) and "+0.011" at
q = 0.05 (0.28x) are Type-M-zone or under-MDE contrasts; their fold CIs exclude zero and the
signs hold on 5/5 and 4/5 folds, so the DIRECTION (a real ranker hurts more than a random one)
is supported and the SIZE is not (S26 L46's rule: say the sign, not the size). Quote them as
"the sign is measured at 5/5 folds; the magnitude is inside the Type-M zone". The built-chain
entry decides; nothing here pre-empts it.
Artefacts: `s27/results/s28_C_readout_cloud_rows.jsonl`, `s28_C_readout_cloud_summary.json`;
my recomputation printed by `s27/s28_D_attack.py` (not persisted; rerunnable in 20 s).

## S28-L17 -- ORACLE EMITTED CEILING, PARTIAL 22/126 (built chain); RECOGNITION ARMS, PARTIAL 43/126 (point cloud): A PARTIAL NOTE, NOT A RESULT (2026-09-14 20:20, A)

Times in S28-L12 and PREREG addendum 3 were stamped ahead of the clock (written about 19:45);
nothing else changes. Both jobs below were killed by the governor at 95 to 97% RAM from the
user's own load and relaunched from their per-target checkpoints (`s28A_recog_126_r2` queued
behind the cap of 2; the chain job is chained to start when it finishes, one job per lane).

ORACLE (theta chosen against the native; the family's ceiling, never a result), BUILT CHAIN,
PARTIAL 22/126 (`s27/results/s28_A_chain_rows.jsonl`, sorted-pdb order 1A13..1KZ2):
  ORACLE circuit optimum through `s12.instrument.project`: mean 0.2934 A (point cloud 0.2865 on
  the same 22; projection price +0.0068), median 0.267, worst 0.955, 22/22 under 2 A, mean
  virtual bond of the emitted chain 3.804 (ideal geometry, as it must be). Production on the same
  22: chain 2.9325 (point cloud 2.8272, projection price +0.1052; S27's DIS chain rows agree on
  these 22 to max |diff| 0.0186 A, so the paired comparator is re-projected here, not reused).
  ORACLE affine-500 emitted: 0.1378 on 22 (S10-5: 0.064 on 126 in the joint-transform frame).
  Reading (PARTIAL): the projection costs the signed combination 0.007 A and the convex average
  0.105 A, because the signed combination is not contracted (S28-L1b). The complete 126-target
  emitted ceiling is posted when the chain job lands.

RECOGNITION (deployable, native-free), POINT CLOUD, PARTIAL 43/126
(`s27/results/s28_A_recog_rows.jsonl`, sorted-pdb order 1A13..2LNG), mean paired d vs
production on the same 43, no CI, no MDE, NOT A RESULT:
  circuit lam 0 (the deployed CVaR-VQE state read as signed weights): +84.9 (median 23.0 A):
    the signs are the initialisation's and the affine sum blows up, as registered (F1 prior);
  circuit lam 0.3 / 1 / 3 at 80 iterations: +0.31 / +0.35 / +0.40 (means 3.24 / 3.27 / 3.33);
    lam 1 at 400 iterations: +0.42 (convergence does not help);
  a500 rand matched / converged (lam 1): +0.56 / +0.65; a75 rand: +0.43 / +0.52;
  simplex rand matched: +0.03; simplex from production, converged: +0.003 (the objective does
    not move production anywhere);
  S-only a500 converged: +0.72; untrained draw 0: +12.2 (median 5.6).
  Lane D's caveat (a), the scale of the terms, mean over 43: at theta0 CVaR -1.38, T H 2.01 and
  lam S~ 1.02 / 3.41 / 10.23 for lam 0.3 / 1 / 3; at the optimum CVaR -1.53, T H 2.9 and lam S~
  0.42 / 1.27 / 3.74. The grid spans minor to dominant; it is not one cell.
  The mechanism, visible already: at lam 1 the optimiser drives S~ to 1.27, BELOW production's
  1.67, while the RMSD rises 0.35 A; and the ORACLE structure sits at S~ 2.09 (S28-L1b). The
  objective orders (circuit optimum) < (average) < (near-native): it is minimised, and its
  minimiser is away from the native. This is S8-9's 37th percentile from the other side and is
  the registered reason F1 was expected not to fire. Verdicts wait for 126/126 and the chain.

## S28-L18 -- INSTRUMENT: THE BUILT CHAIN'S NUMERICAL FLOOR IS 0.02 A PER TARGET (0.0025 A MEAN), NOT 1e-5, WHEN THE INPUT CLOUD DIFFERS AT 1e-13 (2026-09-14 20:28, lane D)
Question: S28-L7 (lane C) states "built-chain contrasts carry a numerical floor near 1e-5 A"
from one target (1A13, input difference 8e-15 A, chain difference 1.4e-5). Lane A's S28-L17
re-projects the production average in its own frame (input identical to the deployed average
to `prod_frame_max_dev` 5.7e-14, S28-L1b) and notes S27's DIS chain rows "agree on these 22 to
max |diff| 0.0186 A". Which floor is right?
From `s27/results/s28_A_chain_rows.jsonl` (22 targets, arm `prod`) against
`s27/results/chain_rows.jsonl :: DIS` (the same targets), point clouds equal to 6 decimals on
all 22 (e.g. 1CS9 3.832166 both): the built chains differ by +0.01865 (1CS9), +0.01036 (1I8E),
-0.00556 (1FUV), -0.00358 (1A13), +0.00332, +0.00316, -0.00278, -0.00237, ...; max |diff|
0.0186, mean |diff| 0.0025, 17 of 22 above 1e-4.
Reading: the production projection is a multi-start optimisation (`s12.instrument.project`,
`core.project.lam_path`, ramah 0.3) with near-tied branches (S26 L88: the branch degeneracy is
worth 0.08 A to a perfect chooser); a 1e-13 difference in the input cloud flips the chosen
branch on some targets and the emitted chain moves by 1e-3 to 2e-2 A. Lane C's 1e-5 was a
target on which no branch flipped. Consequences for this sprint:
(a) a built-chain contrast between two arms projected from BIT-IDENTICAL clouds (lane C's
    readouts against S27's rows, S28-L7's reproduction at 0.0) has no floor; a contrast whose
    two sides were projected from clouds that differ at floating-point level (lane A's frame
    against S27's rows; any re-implementation of the average) carries a per-target floor of up
    to 0.02 A and a mean floor of about 0.003 A. Lane A's choice to re-project its own comparator
    (S28-L17) is the right one; the two sides of every chain contrast must come from the same
    code path, and the entry must say which.
(b) an effect of 0.003 A on the built chain is inside this floor whatever its CI says; the
    sprint's MDEs (0.03 to 0.10 A) are above it, so no verdict so far is touched.
(c) the identity cells (DIVW(0,0), k = 74, q = 0) reproduce production to 1e-13 on the point
    cloud and to 1e-5 on the chain ONLY when no branch flips; a lane that reports an identity
    check on the chain must report the max over targets, not one target.
Verdict: an instrument property, recorded; S28-L7's "1e-5" is superseded by "1e-5 to 2e-2 per
target, 0.003 mean" (`s27/RETRACTIONS_S28.md` R1, a scope correction, not a retraction of a
result). No S28 verdict changes.
Artefacts: `s27/results/s28_A_chain_rows.jsonl`, `s27/results/chain_rows.jsonl`; the per-target
table is in this entry.

## S28-L18 -- RECOGNITION ON THE POINT CLOUD (INTERMEDIATE BASIS, 126/126, NO VERDICT): EVERY SIGNED ARM IS WORSE THAN PRODUCTION; THE OBJECTIVE IS MINIMISED AND ITS MINIMISER IS AWAY FROM THE NATIVE; THE BUILT-CHAIN VERDICT WAITS FOR JOB s28A_chain_primary (2026-09-14 20:50, A)

Artefacts: `s27/results/s28_A_recog_rows.jsonl` (126 rows, every arm of PREREG section 3),
`s27/results/s28_A_summary.json` (`cloud`, `contrasts`, `text`), jobs
`s26/jobs_done/s28A_recog_126.json` (killed at 43/126, 713 s) and `s28A_recog_126_r2.json`
(resumed to 126/126, 1420 s, peak RSS 0.336 GB). The lam = 0 soundness gate: my loop's `p`
equals `core.quantum.run_cvar_vqe`'s bit-for-bit on all 126 targets (max |dp| 0.00e+00).
Undefined readouts (|sum psi| < 1e-9): 0 in every arm, 0 of 2016 untrained draws; rule (e) moot.
Basis: POINT CLOUD, an intermediate; no arm is called anything here; the verdict is written on
the built chain (S28-L1b's ORACLE emitted ceiling and the circuit arms are in job
`s28A_chain_primary_1`, running).

Means (point cloud, `s28_A_summary.json :: cloud`): production 3.0483; circuit lam 0.3 / 1 / 3
at 80 iterations 3.3832 / 3.3850 / 3.4382; circuit lam 0 87.03 (median 13.94); lam 1 at 400
iterations 3.5315; untrained circuit mean of 16 draws 22.16 (best-of-16, an order statistic,
3.088); a500 rand matched / converged (lam 1) 3.5755 / 3.6869; a75 rand 3.4880 / 3.5387;
random 27-subspace mean of 8 (rand, matched, lam 1) 3.5445; simplex rand matched (lam 1) 3.0563;
simplex from production, converged (lam 1) 3.0522; a75 from production (lam 0.3) 3.0460;
S-only a500 converged 3.7530; S-only simplex from production 3.2383.

  circ_l0.3_i80 vs production (point cloud)
    a 3.3832 (med 3.1282)   b 3.0483 (med 2.8373)   n=126
    effect +0.3349   median +0.2074   SE 0.0535   MDE 0.1500   effect/MDE +2.23
    iid  CI95 [+0.2353, +0.4403]
    fold CI95 [+0.2058, +0.4537]   folds same sign 5/5   per-fold 0:+0.469 1:+0.471 2:+0.392 3:+0.280 4:+0.114
    32W/94L/0T   worst degradation +2.1461 (9BAF)   p90 +1.0105   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.4246 vs uniform-effect null p10/p50/p90 +0.3590/+0.4220/+0.4898 -> pctile 0.519
    VERDICT: WORSE

  circ_l1_i80 vs production (point cloud)
    a 3.3850 (med 3.1263)   b 3.0483 (med 2.8373)   n=126
    effect +0.3367   median +0.2331   SE 0.0606   MDE 0.1697   effect/MDE +1.98
    iid  CI95 [+0.2198, +0.4556]
    fold CI95 [+0.1739, +0.4771]   folds same sign 5/5   per-fold 0:+0.395 1:+0.539 2:+0.443 3:+0.328 4:+0.051
    36W/90L/0T   worst degradation +2.4460 (5MXS)   p90 +1.2211   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.4401 vs uniform-effect null p10/p50/p90 +0.3678/+0.4393/+0.5104 -> pctile 0.508
    VERDICT: WORSE

  circ_l3_i80 vs production (point cloud)
    a 3.4382 (med 3.2587)   b 3.0483 (med 2.8373)   n=126
    effect +0.3899   median +0.3201   SE 0.0583   MDE 0.1633   effect/MDE +2.39
    iid  CI95 [+0.2764, +0.5028]
    fold CI95 [+0.2176, +0.5300]   folds same sign 5/5   per-fold 0:+0.562 1:+0.509 2:+0.505 3:+0.353 4:+0.087
    35W/91L/0T   worst degradation +2.6900 (2RUO)   p90 +1.1795   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.4855 vs uniform-effect null p10/p50/p90 +0.4093/+0.4829/+0.5602 -> pctile 0.520
    VERDICT: WORSE

  circ_l0_i80 vs production (point cloud)
    a 87.0262 (med 13.9399)   b 3.0483 (med 2.8373)   n=126
    effect +83.9779   median +11.6371   SE 24.7247   MDE 69.2687   effect/MDE +1.21
    iid  CI95 [+44.8255, +142.4829]
    fold CI95 [+55.3432, +115.7913]   folds same sign 5/5   per-fold 0:+139.418 1:+48.708 2:+59.327 3:+59.954 4:+103.779
    4W/122L/0T   worst degradation +2703.3080 (6B9K)   p90 +188.7409   power 0.92  Type-M 1.05
    concentration: drop-top10 +91.1882 vs uniform-effect null p10/p50/p90 +58.5632/+88.4532/+129.5828 -> pctile 0.541
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.05x]

  diag_circ_l1_i400 vs production (point cloud)
    a 3.5315 (med 3.4231)   b 3.0483 (med 2.8373)   n=126
    effect +0.4832   median +0.3641   SE 0.0636   MDE 0.1781   effect/MDE +2.71
    iid  CI95 [+0.3600, +0.6087]
    fold CI95 [+0.3485, +0.5945]   folds same sign 5/5   per-fold 0:+0.629 1:+0.568 2:+0.572 3:+0.456 4:+0.243
    29W/97L/0T   worst degradation +2.5526 (2MK7)   p90 +1.4327   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.5852 vs uniform-effect null p10/p50/p90 +0.5062/+0.5829/+0.6660 -> pctile 0.515
    VERDICT: WORSE

  untrained circuit (mean of 16 draws) vs production (point cloud)
    a 22.1615 (med 13.5413)   b 3.0483 (med 2.8373)   n=126
    effect +19.1131   median +10.8237   SE 2.9939   MDE 8.3876   effect/MDE +2.28
    iid  CI95 [+14.2443, +25.3897]
    fold CI95 [+15.6080, +22.4128]   folds same sign 5/5   per-fold 0:+13.999 1:+21.211 2:+14.749 3:+21.787 4:+23.353
    0W/126L/0T   worst degradation +313.6326 (2MLQ)   p90 +33.1158   power 1.00  Type-M 1.00
    concentration: drop-top10 +20.4835 vs uniform-effect null p10/p50/p90 +16.6128/+20.0990/+24.4410 -> pctile 0.547
    VERDICT: WORSE

  circ_l1_i80 vs random-27-subspace control (mean of 8, matched budget) (point cloud)
    a 3.3850 (med 3.1263)   b 3.5445 (med 3.4406)   n=126
    effect -0.1595   median -0.0661   SE 0.0503   MDE 0.1410   effect/MDE -1.13
    iid  CI95 [-0.2589, -0.0626]
    fold CI95 [-0.2500, -0.0583]   folds same sign 4/5   per-fold 0:-0.271 1:+0.010 2:-0.093 3:-0.119 4:-0.283
    74W/52L/0T   worst degradation +0.9635 (8T61)   p90 +0.3167   power 0.89  Type-M 1.07
    concentration: drop-top10 -0.0399 vs uniform-effect null p10/p50/p90 -0.0934/-0.0421/+0.0096 -> pctile 0.524
    VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.07x]

  circ_l1_i80 vs a500_rand_c_l1 (point cloud)
    a 3.3850 (med 3.1263)   b 3.6869 (med 3.8446)   n=126
    effect -0.3019   median -0.1372   SE 0.0633   MDE 0.1772   effect/MDE -1.70
    iid  CI95 [-0.4294, -0.1850]
    fold CI95 [-0.4069, -0.1845]   folds same sign 5/5   per-fold 0:-0.436 1:-0.131 2:-0.168 3:-0.327 4:-0.413
    85W/41L/0T   worst degradation +1.8996 (2EFZ)   p90 +0.3955   power 1.00  Type-M 1.00
    concentration: drop-top10 -0.1612 vs uniform-effect null p10/p50/p90 -0.2467/-0.1627/-0.0864 -> pctile 0.507
    VERDICT: BETTER

  circ_l1_i80 vs simplex_rand_m_l1 (point cloud)
    a 3.3850 (med 3.1263)   b 3.0563 (med 2.7179)   n=126
    effect +0.3287   median +0.2082   SE 0.0566   MDE 0.1587   effect/MDE +2.07
    iid  CI95 [+0.2230, +0.4394]
    fold CI95 [+0.2066, +0.4489]   folds same sign 5/5   per-fold 0:+0.309 1:+0.481 2:+0.471 3:+0.335 4:+0.105
    35W/91L/0T   worst degradation +3.1255 (7BX2)   p90 +1.0374   power 1.00  Type-M 1.00
    concentration: drop-top10 +0.4187 vs uniform-effect null p10/p50/p90 +0.3423/+0.4142/+0.4888 -> pctile 0.527
    VERDICT: WORSE

Reading (intermediate basis, no verdict):
1. Every signed arm is worse than production. The circuit at lam 0.3 / 1 / 3 is +0.335 / +0.337
   / +0.390 A at 2.0 to 2.4x MDE, fold CI above zero, 5/5 folds; lam 0 (the deployed state read
   as signed weights) is +84 A because its signs are the initialisation's (D's caveat (d)); the
   untrained circuit is +19 A. F1's registered prior held on this basis.
2. The optimiser works and that is the problem. At lam 1 the surrogate falls from 3.41 at theta0
   to 1.34 at the optimum (`parts`), BELOW production's 1.674 and far below the ORACLE
   structure's 2.094 (S28-L1b); 400 iterations lower S~ further (1.263) and RAISE the RMSD to
   3.53 (+0.483, 2.7x MDE). Across the classical families the ordering is the same: the
   converged unconstrained affine optimum has the lowest S~ (1.06 to 1.08) and the worst RMSD
   (3.69 to 3.77). The objective term is minimised, and its minimiser is not near the native:
   S15's "optimising a structure against the distogram objective gives 3.321" and S8-9's 37th
   percentile, reproduced for the emitted structure's own score, on the amplitude family.
3. The grid is not degenerate (D's caveat (a)): lam S~ at the optimum is 0.42 / 1.27 / 3.74
   against CVaR -1.53 and T H 2.7 to 2.9; the three cells give three different structures
   (S~ 1.475 / 1.338 / 1.296) and the same verdict.
4. F2 preview (moot unless F1 fires on the chain): among SIGNED families at the same objective
   and budget the circuit is the least bad (vs the random 27-subspace mean-of-8 -0.16 A at
   1.13x MDE, Type-M zone; vs unconstrained a500 converged -0.30 at 1.7x MDE; vs S-only a500
   -0.37 at 2.0x). It is WORSE than the CONVEX simplex under the same objective (+0.329 at 2.07x
   MDE, 5/5): the convex family, started at random, lands at 3.056, production's value. So the
   parameter count does not do it (D's S28-L1 wording), and neither does the circuit's family:
   the sign freedom itself is what the objective misuses.
5. The prod-init controls show the objective does not move production anywhere useful: a75 from
   production at lam 0.3 is 3.0460 (-0.002), the simplex from production 3.052 to 3.115.
6. Geometry: every signed optimum is uncontracted (mean virtual bond 3.40 to 3.54, Rg 6.6 to
   6.7, against production's 2.96 / 6.21 and the pool's 3.81 / 6.80), with 47% negative weights
   and negative mass 5.7 to 6.5. The signed readout escapes the averaging contraction, as the
   ORACLE optimum did; it does not escape the objective.
The chain projections of these arms are running; nothing above is a verdict.

## S28-L19 -- NUMBERING CORRECTION: lane A's entry headed "S28-L18 -- RECOGNITION ON THE POINT CLOUD ..." is S28-L18b (2026-09-14 20:55, A)

It collided with lane D's S28-L18 (the built chain's numerical floor, posted 20:28 while A was
writing). Content unchanged. Lane D's S28-L18 is accepted as it stands: every chain contrast in
lane A is between arms projected by the SAME code path in the SAME job (`s28_A_amp.py chain`,
production re-projected from A's own frame, never S27's rows), and the per-target floor of up
to 0.02 A (mean 0.003) is below every MDE A will quote on the chain.

## S28-L20 -- ADVERSARY CHECK OF S28-L18b (lane A's recognition arms, point cloud, intermediate) (2026-09-14 20:58, lane D)
Question: does S28-L18b's "every signed arm is worse" stand as an intermediate negative, and what
must the built-chain verdict entry carry? Checks (all recomputed from
`s27/results/s28_A_recog_rows.jsonl` and `s28_A_oracle_rows.jsonl :: prod`):
- Reproduction: circ_l1_i80 vs production +0.3367 (SE 0.0606, MDE 0.1697, 1.98x, 36W/90L, 5/5),
  circ_l0.3_i80 +0.3349 (2.23x), circ vs simplex_rand_m_l1 +0.3287 (2.07x): identical to the
  entry to four decimals (bootstrap CIs differ in the third decimal by the label seeding).
- Leakage / poison: the deployable path is native-free end to end (`tests/test_s28_D.py`,
  S28-L10 item 2); the lam = 0 gate passes 126/126 (max |dp| 0.0); undefined readouts 0/2016.
- Ties, concentration, seeds: no arm is positive, so none apply; the concentration null sits at
  the 50th percentile on every WORSE contrast (the harm is spread, not carried by ten targets).
- FAIL18 / other-108 split (mine; the entry has none): circ_l1_i80 is WORSE on the 108 (+0.386,
  SE 0.061, 2.26x MDE, fold CI [+0.295, +0.458], 5/5, 27W/81L) and NULL on FAIL18 (+0.040,
  SE 0.206, 0.07x, 9W/9L); circ_l0.3_i80 the same (+0.390 / +0.004). The convex simplex from a
  random start at matched budget (simplex_rand_m_l1, 3.0563 overall, 0.09x) is -0.203 on
  FAIL18 (SE 0.100, 0.72x MDE, fold CI [-0.275, -0.065], 3/4 folds, 14W/4L) and +0.043 on the
  108 (0.50x): S27 L9's regime pattern (every alternative helps the 18 and costs the 108) read
  back in a re-weighting of the whole pool; Type-M zone on n = 18, not a result, recorded so it
  is not rediscovered.
- Cosmetic-variant test (coordinator steer 3): the signed arms are NOT re-parameterisations of
  the convex top-75 average: 47% negative weights, negative mass 5.7 to 6.5, mean virtual bond
  3.40 to 3.54 against production's 2.96; they leave the hull and escape the contraction, and
  the objective sends them the wrong way. The prod-init convex arms (a75 from production
  -0.002; simplex from production +0.004 to +0.066) ARE re-parameterisations and sit at
  production. New relative to S26/S27: yes, the first signed readout on record (S10-5's affine
  bound was ORACLE only); it goes beyond S27 section 6 by testing a non-convex consumer.
- The three-way split (coordinator steer 4), read for lane A: (i) OBJECTIVE quality: the
  native-free S~ orders (circuit optimum 1.34) < (production 1.67) < (ORACLE structure 2.09),
  so the objective's minimum is away from the native; this is the failing leg. (ii)
  OPTIMISATION quality: the circuit reaches S~ 1.34 at 80 iterations and 1.26 at 400, above
  the converged unconstrained optimum (1.06 to 1.08) and below production; the optimiser
  works, and working harder makes the structure worse (+0.483 at 400 iterations, 2.7x MDE).
  (iii) the emitted structure: worse at every lam. A gain would not have been a quantum result
  in any case: the convex simplex under the same objective lands at production, and the
  signed families are worse in proportion to how well they minimise S~.
Verdict: STANDS as an intermediate negative on the point cloud (F1's prior held on this basis).
Two caveats for the chain entry: (a) quote the matched-budget a500 (3.5755, -0.19 in the
circuit's favour) beside the converged one (3.6869, -0.30) when saying the circuit is the least
bad signed family; the budgets differ 25x; (b) the "least bad vs the random 27-subspace"
contrast is Type-M (1.13x, 4/5 folds): say the sign. The built-chain entry decides F1 and is
attacked in its turn.
Artefacts: `s27/results/s28_A_recog_rows.jsonl`, `s28_A_summary.json`; the strata above are
`ST.compare` on the FAIL18 mask (`s12.instrument.FAIL18`).

## S28-L21 -- HOPPING ON THE POINT CLOUD (INTERMEDIATE BASIS, 126/126, NO VERDICT): THE DEPLOYED READOUT MOVES BY AT MOST 0.013 A AT ANY J; THE TWO OTHER READOUTS ARE +0.24 TO +0.33 A ABOVE PRODUCTION AT EVERY J; THE EIGENSOLVER AT J = 1 RE-SELECTS THE DIS TOP-75; THE DEPARTURE TABLE AT n = 126 AND THE S28-L11 REWORDING OF F5; THE BUILT-CHAIN VERDICT WAITS FOR JOB s28B_chain (2026-09-14 21:20, B)

Point cloud is the SCREENING basis (contract addendum 1, item 13); nothing here is a verdict.
Run: `s28B_run` (64 targets) + `s28B_run2` (1) + `s28B_run3` (61), all resumed from
`s27/results/s28_B_rows.jsonl` (4,914 rows = 126 targets x 39 arms, no duplicates at the
64/65/66 boundaries); kills at 20:04 and 20:10 were the governor's stall breaker under the
user's own load (`s26/governor.log`), peak RSS 0.35 GB, 27 to 45 s per target. Statistics
`s27/s28_B_analyse.py` -> `s27/results/s28_B_summary.json` (every ST.fmt block under `fmt`).

Anchors: J = 0 R1 seed 0 = 3.058016, seed 1 = 3.056181, bit-identical to
`s27/results/vqe_rows.jsonl :: DIS` on 126/126 targets both seeds (`summary.json ::
anchors.s27_bit_identical_s*`; lane D's own check S28-L9); DIS top-75 3.048338. The
eigensolver at J <= 0.3 is the one-hot argmin and its R1/R2 read 3.4540, the S25 argmin
selector's number.

**Means (point cloud, A; `summary.json :: arms`), seed 0 / seed 1 for the VQE:**

    arm            R1 tail          R2 weighted      R3 p-top-75      m       PR
    VQE J=0        3.0580 / 3.0562  3.3697 / 3.3636  3.3096 / 3.3207  74 / 71  428 / 412
    VQE REAL 0.1   3.0546 / 3.0558  3.3665 / 3.3643  3.3097 / 3.2999  74 / 71  428 / 414
    VQE REAL 0.3   3.0564 / 3.0558  3.3657 / 3.3635  3.2876 / 3.2992  74 / 72  427 / 414
    VQE REAL 1     3.0541 / 3.0533  3.3671 / 3.3661  3.2906 / 3.3215  74 / 72  427 / 415
    VQE REAL 3     3.0594 / 3.0623  3.3813 / 3.3656  3.3207 / 3.3055  79 / 76  436 / 425
    GS   REAL 0.1  3.4540 (m=1)     3.4520           3.0446            1        1.0
    GS   REAL 1    3.2185 (m=5.5)   3.0447           3.0597            5.5      62
    GS   REAL 3    3.0697 (m=32)    3.1645           3.0829           32.5     305
    (PERM and RAND: within 0.01 of REAL on every VQE cell; GS PERM/RAND J = 1 R3 3.0452 /
     3.0583; production DIS top-75 3.0483; random-75 null 3.4209, S27 `pool_rows.jsonl`)

**F1 (each VQE arm vs the SAME readout at J = 0, same seed; 72 contrasts):** none clears its
MDE with the fold CI excluding zero on both seeds. Nine cells sit at 0.7x MDE or better and
all nine go to the chain (contract item 13). The largest, ST.fmt verbatim:

  F1 vqe|s1|REAL|J0.1 R3 - vqe|s1|NONE|J0 R3 (point cloud)
    a 3.2999 (med 3.1363)   b 3.3207 (med 3.1716)   n=126
    effect -0.0208   median +0.0000   SE 0.0071   MDE 0.0199   effect/MDE -1.05
    iid  CI95 [-0.0358, -0.0082]
    fold CI95 [-0.0338, -0.0089]   folds same sign 5/5   per-fold 0:-0.045 1:-0.005 2:-0.026 3:-0.023 4:-0.006
    58W/31L/37T   worst degradation +0.1038 (6A5J)   p90 +0.0321   power 0.83  Type-M 1.10
    concentration: drop-top10 -0.0016 vs uniform-effect null p10/p50/p90 -0.0071/-0.0018/+0.0029 -> pctile 0.514
    VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.10x]

  and its two killers, same basis:
  F1 vqe|s0|REAL|J0.1 R3 - vqe|s0|NONE|J0 R3: effect +0.0001, MDE 0.0133, 0.01x, fold CI
    [-0.0076, +0.0086], 52W/39L/35T (the seed-0 twin: does not replicate);
  F1 vqe|s1|PERM|J0.1 R3 - vqe|s1|NONE|J0 R3: effect -0.0219, 0.82x, fold CI [-0.0408,
    -0.0082], 46W/34L/46T (the PERM control, correspondence destroyed, shows the same movement:
    F3 says the specific graph contributed nothing).
  The deployed readout's largest cell, F1 vqe|s0|REAL|J0.1 R1: effect -0.0034, MDE 0.0046,
    0.74x, 22W/10L/94T (94 identical sets), seed 1 -0.0004 (0.17x): not a result.

**S28-L2(a), every readout against PRODUCTION (DIS top-75 uniform, point cloud):** R1 at
J = 0: +0.0097 (0.43x, fold CI [-0.0005, +0.0250], 60W/64L/2T); R1 at every J > 0: +0.005 to
+0.014, all inside MDE. R2 at J = 0: +0.3214 (1.41x, fold [+0.1661, +0.4629], 48W/78L, WORSE);
R3 at J = 0: +0.2613 (1.33x, fold [+0.1315, +0.3879], 46W/80L, WORSE); at every J > 0, graph
and seed R2 is +0.313 to +0.333 and R3 +0.226 to +0.282, all WORSE with the fold CI above zero.
The R3 "gain" above is a readout that is 0.26 A above production climbing 0.02 A toward it.

**F2 (VQE vs the exact ground state, same J, R2/R3):** at J = 1 and 3 the eigensolver is
BETTER than the circuit on both readouts beyond MDE (R3 J = 1: +0.2309, 1.20x, fold [+0.0994,
+0.3575], 47W/79L; R2 J = 3: +0.2168, 1.87x); at J <= 0.3 GS-R2 is the argmin (3.45) and the
circuit's R2 is inside MDE of it. The registered prior (eigensolver wins or ties) held.

**The eigensolver as a selector (GS vs the deployed J = 0 R1, and vs production):** GS REAL
J = 1 R3 = 3.0597: +0.0017 vs J = 0 R1 (0.07x) and +0.0114 vs production (0.45x, fold CI
[+0.0011, +0.0223], 49W/65L/12T); its p-top-75 has Jaccard 0.960 with the DIS top-75, i.e. the
hopping ground state at J = 1 RE-SELECTS the distogram's own top-75 (the probe's 0.92 to 0.97
reading back at n = 126). GS REAL J = 1 R2 (p_gs-weighted over the whole pool) = 3.0447: -0.0037
vs production (0.09x, fold CI [-0.0221, +0.0176], 79W/47L) and -0.0134 vs J = 0 R1 (0.31x):
inside MDE. GS J = 3 R1 (m = 32.5) 3.0697, +0.012 vs J = 0 R1 (0.19x). No eigensolver readout
is beyond its MDE against production.

**F3/F4 (REAL vs PERM / RAND, same J, seed, readout):** 72 contrasts, none beyond MDE; the
largest is R1 J = 1 seed 0 REAL - PERM = +0.0092 (0.57x, PERM better). The graph's specific
similarity structure and its degree sequence are indistinguishable at this instrument.

**S28-L2(e), the best J as an order statistic (`ST.best_of_k_within` over the 4 non-zero
rungs, 18 (graph, seed, readout) cells):** every REAL and RAND cell reads "NOT A SIGNAL"
(split-half transfers -14% to +15% of the oracle; k_eff 3.2 to 4.0); the two cells with a
surviving residual are PERM R1/R2 seed 0 (-0.007 / -0.006 transfers, in a CONTROL). The best J
for R1 by the point-cloud mean is J = 1 on both seeds (3.0541 / 3.0533; -0.0039 vs J = 0 at
0.29x, fold [-0.0084, +0.0001], 35W/38L/53T); it goes to the chain as registered.

**The departure table at n = 126 (S28-L11 caveat (e); native-free; `s28_B_rows.jsonl`; REAL,
seed 0 / seed 1; `share` = <psi|A|psi> / sum |psi_i| A_ij |psi_j|, which for a near-rank-one A
coincides with the sign coherence, a consistency check and not an independent number
(S28-L11(3)):**

    J     m            holes  PR        Jac(ptop75,DIS75)  mass on DIS75  TV vs J=0      hop            bound   share / coh
    0     74.1 / 71.2  <=2    428/412   0.177 / 0.180      0.187 / 0.193  -              0              -       - / 0.00
    0.1   74.3 / 71.2  <=2    428/414   0.177 / 0.180      0.186 / 0.193  0.025 / 0.027  -0.001/-0.001  0.90    0.00 / 0.00
    0.3   74.1 / 71.9  <=2    427/414   0.180 / 0.179      0.187 / 0.191  0.050 / 0.047  -0.000/-0.001  0.90    0.00 / 0.00
    1     73.7 / 71.8  <=1    427/415   0.181 / 0.177      0.188 / 0.191  0.098 / 0.088  0.007 / 0.029  0.89    0.01 / 0.03
    3     79.3 / 76.3  <=2    436/425   0.163 / 0.164      0.177 / 0.182  0.155 / 0.152  0.318 / 0.238  0.91    0.34 / 0.25
    GS 1  5.5          <=2    62        0.960              0.811          0.956          0.499          0.499   1.00 / 1.00
    GS 3  32.5         <=2    305       0.842              0.373          0.994          0.911          0.911   1.00 / 1.00
    (random-75 vs random-75 Jaccard 0.082; uniform mass on 75 of 500 = 0.150; `gate_set_equality`
     passes 4,914/4,914 with `equality` true on every row, the operator's property, not a finding)

The S28-L8b partial numbers restated at 126: the circuit collects 34% / 25% (seeds 0 / 1) of
its same-sign bound at J = 3 (entry: 44% on 19 targets; S28-L11: 33% on 30), 0.7% / 3.2% at
J = 1 (entry: "essentially none"); TV 0.155 / 0.152 at J = 3 (entry 0.135); Jaccard 0.163 /
0.164 (entry 0.144); the eigensolver 0.499 / 0.911 at J = 1 / 3 (entry 0.503 / 0.910). Signs
and orderings unchanged. The retained sets: the VQE's alpha-tail (R1) has mean pairwise RMSD
2.48 A, the eigensolver's top-75 2.36 to 2.45, the VQE's p-top-75 (R3) 4.06 to 4.12: the
near-uniform state's most probable 75 is a dispersed set, which is why R3 is 0.26 A worse than
production and the eigensolver's R3 is not.

**S28-L11 caveats (a) to (d), accepted and reworded here for the record:** (a) the hop-only
slope is -1.7 to -1.8 per qubit (-1.844 over n = 4..9, -1.700 over the padding-free 4..8);
(b) the "-1.27" diagonal-control rate is withdrawn: on the padding-free registers the
diagonal same-spectrum control decays at -1.856, the same rate as A within the fit noise, so
the decay is the near-rank-one spectrum's and the eigenvector's orientation does not
measurably change it; (c) "gradient-invisible" and "cannot train" are withdrawn: the measured
statements are the variance ratio (7,300x at J = 1, 820x at J = 3, n = 9, depth 3) and the
slope, and the circuit collects a quarter to a third of its same-sign bound at J = 3;
(d) "to 0.1%" is withdrawn: a 7,300x variance ratio is 1.2% in gradient magnitude.

Built chain: job `s28B_chain` (18 arms x 126: J = 0 R1 both seeds; best-J R1 J = 1 both
seeds; the nine 0.7x cells; the seed-0 twin and the J = 0 R3 comparators of the largest cell;
GS REAL J = 1 R3 and R2) is running under jobrun with per-(arm, target) checkpoints
(`s27/results/s28_B_chain_rows.jsonl`); the verdict entry is written on it. `s28B_split` (the
three-way split's middle leg, S28-L2(b)) and `s28B_mladder` (the T9 decomposition) run beside
it. No verdict here.
Artefacts: `s27/results/s28_B_rows.jsonl`, `s28_B_summary.json`, `s28_B_chain_arms.txt`,
`s26/jobs_done/s28B_run.json`, `s28B_run2.json`, `s28B_run3.json`.

## S28-L22 -- ADVERSARY CHECK OF S28-L21 (lane B's hopping arms, point cloud, intermediate) (2026-09-14 21:12, lane D)
Question: does the intermediate stand, and is the most promising cell on the board dead?
Every number recomputed from `s27/results/s28_B_rows.jsonl` (4,914 rows) with
`s27/s28_D_attack.py --arm-expr "source|seed|graph|J"`.
- The most promising number on the board at 21:00 was F1's largest cell, R3 (p-top-75) at
  seed 1, REAL, J = 0.1 vs its J = 0 twin: -0.0208, SE 0.0071, MDE 0.0199, 1.05x, fold CI
  [-0.0336, -0.0089], 5/5, 58W/31L/37T, Type-M zone (reproduced to four decimals). The one
  experiment most likely to kill it was already in the rows and I confirm both killers: the
  seed-0 twin is +0.0001 (0.01x, MDE 0.0133) and the PERM control at the same seed and J is
  -0.0219 (0.82x, fold CI [-0.0408, -0.0082], 5/5): a control with the correspondence destroyed
  moves as much, so the graph did nothing (F3); RAND is -0.0139 (0.42x). Against PRODUCTION the
  same arm is +0.2516 (1.22x MDE, WORSE): a readout 0.26 A above production climbing 0.02 A
  toward it, exactly S28-L2(a)'s trap, and the entry says so. DEAD on the point cloud; the
  chain projection of the cell (registered) cannot revive it because the comparator is
  production.
- The deployed readout (R1): every J within 0.014 A of production and inside MDE; the largest
  F1 cell (seed 0, J = 0.1, -0.0034, 0.74x) has 94 identical retained sets of 126 and a seed-1
  twin at 0.17x. Not a result, and the entry does not call it one.
- The eigensolver: GS REAL J = 1 R2 vs production -0.0037 (0.09x, MDE 0.0400, fold CI
  [-0.0221, +0.0176]) with a LOPSIDED 79W/47L and median -0.018 against a mean of -0.004: W/L
  and the median-vs-mean gap cannot diagnose a result (project memory: a raw W/L or drop-top
  read is not a test; the concentration null sits at the 52nd percentile); FAIL18 -0.064
  (0.55x), other-108 +0.006 (0.15x). Under MDE on every basis. It is in the chain list; the
  chain decides and the W/L is not to be quoted as evidence.
- Order statistics: the best J is priced ("NOT A SIGNAL" on every REAL and RAND cell); the
  NINE F1 cells at 0.7x or better that go to the chain are nine of 72 correlated contrasts
  (same seeds and readouts across J; k_eff 3.2 to 4.0 per family): a chain positive among them
  is a best-of-nine and must be priced with `best_of_k_within` over the nine before it is
  called anything.
- Three-way split: (i) Hamiltonian quality: the exact ground state at J = 1 re-selects the DIS
  top-75 (Jaccard 0.960) and at J = 3 a 305-state delocalised typicality mode (Jaccard 0.842,
  mass on DIS75 0.373); its readouts are inside MDE of production; (ii) optimisation quality:
  the circuit's state departs from J = 0 by TV 0.03 to 0.16 and collects 0 to 34% of its
  same-sign hopping bound (the `s28B_split` F-values are pending); (iii) emitted structure:
  no readout beyond MDE of its J = 0 twin on both seeds. No leg carries a gain; the F2 prior
  (eigensolver wins or ties) held. Nothing quantum is claimed and nothing could be.
- S28-L11's caveats (a) to (e) are all accepted and reworded in the entry; the departure
  numbers are restated at 126 with signs and orderings unchanged.
Verdict: STANDS as an intermediate negative (F1, F3, F4 silent; F2's prior held). The chain
entry decides, with the caveats above (production as the comparator; best-of-nine pricing;
no W/L quoted).
Artefacts: `s27/results/s28_B_rows.jsonl`, `s28_B_summary.json`; my recomputations are
rerunnable in under a minute with the command line above.

## S28-L23 -- ADVERSARY CHECK OF THE SECOND-WAVE PREREGS: PREREG_S28_A addendum 4 (A2, the objective's local behaviour) and PREREG_S28_B addendum 1 (B2, a spread-spectrum hopping graph) (2026-09-14 21:22, lane D)
Question: as S28-L1 to L3, for the two addenda written before any second-wave number (A2 at
21:30, B2 at 21:40, both after the point-cloud intermediates and before the chain verdicts).
A2 (addendum 4). Falsifier (some e beats production on the built chain beyond MDE, fold CI,
5/5, AND beats the random-direction mean beyond MDE): falsifiable. Controls: (i) a random
direction at the same RMS displacement, rigid-body removed, is the matched control in the
operator's space; (ii) the circuit-family one-step arm against its own projected baseline is
the right pairing for "does the family's bias help at the local scale". The ORACLE cosine is
labelled in every sentence and never enters the ladder; the NaN-poison of the ladder is
registered as a test. Verdict: STANDS WITH CAVEAT. Caveats before the first A2 number:
(a) THE COSINE NEEDS ITS NULL. In 3n - 6 shape dimensions a random unit direction has
    E|cos| about sqrt(2 / (pi (3n - 6))) = 0.12 to 0.18 for n = 9..16, not zero; "cosine
    near zero" is only readable against the 8 random directions' cosines on the same target.
    Report the ORACLE cosine beside the mean and the 95th percentile of the random directions'
    |cos| per target, and the FAIL18 / 108 split of both.
(b) e in {0.1, 0.3, 1.0} is a grid of three: any "best e" is priced with `best_of_k_within`
    (k_eff, split-half), and the same for the circuit one-step arm's e.
(c) LIKE FOR LIKE ON THE CHAIN. The random control is 8 draws on the point cloud and draws
    0 and 1 on the chain; the chain contrast must pair the step against the mean of the SAME
    two draws (not the point-cloud eight), and say so; a best-of-2 is an order statistic.
(d) The step is taken from C0 in the point cloud and projected; production is the projection
    of C0 itself in the same job (the S28-L18 floor applies only if the code paths differ; the
    addendum re-projects production in the same job, which is right).
(e) The circuit one-step arm's baseline C(theta_P) is the family's nearest point to C0 with a
    non-zero residual (the family need not contain C0); quote that residual's RMSD to
    production so "degrades from its own baseline" is not "starts below production".
B2 (addendum 1). Condition (b) of the brief (the mechanism is the near-rank-one spectrum) is
what S28-L11 found, and the addendum cites it correctly; condition (a) waits for the S28B chain
verdict and gates only the endpoint arms, as it should. Falsifier for "the spectrum was the
mechanism" (kNN hop-only slope shallower than -1.0 per qubit AND the circuit collects more than
half its same-sign bound at J = 1) is concrete and two-sided (the addendum says what a slope at
or below -1.7 would mean). Controls: PERM, RAND (Sinkhorn on the binary kNN degrees), J = 0,
both seeds, production as the comparator: matched. Verdict: STANDS WITH CAVEAT. Caveats:
(a) THE PERRON VECTOR IS STILL NEAR-UNIFORM. For D^-1/2 A D^-1/2 the top eigenvector is
    proportional to sqrt(degree); on a near-regular kNN graph it is again 0.95+ overlapped with
    the uniform state, so the "typicality projector" component does not go away, it is joined
    by a spread remainder. Report the top eigenvector's uniform overlap beside lambda_2 /
    lambda_1 (registered), and decompose the hop-only variance into the rank-one part and the
    remainder as `s28_B_rank1.py` did for the Gaussian graph, so the reading is "the remainder
    carries X% of the variance" and not "the spectrum is spread".
(b) DISCONNECTED GRAPHS. A symmetric kNN graph at k = 5 on 500 windows can have more than one
    component; then lambda_1 = 1 is degenerate (one Perron vector per component) and `eigh`'s
    "ground state" at large J is an arbitrary combination. Count the components (registered)
    and label every GS row on a disconnected graph DEGENERATE, as R3-GS at J = 0 is.
(c) The endpoint prior (WORSE or null, the consistency mechanism) is the right one; the F1
    comparator is production, per S28-L2(a), and the addendum says so.
(d) The B2 endpoint and chain jobs are gated on the S28B chain verdict AND on my check of it
    (contract addendum 1); no B2 endpoint number is read before both are posted.
Artefacts of this check: this entry; `s27/PREREG_S28_A.md` (addendum 4), `s27/PREREG_S28_B.md`
(addendum 1). The C2 prereg (`s27/PREREG_S28_C2.md`) is checked when it lands.

## S28-L23b -- A2.1 ORACLE DIAGNOSTIC: AT THE PRODUCTION POINT THE SHIPPED OBJECTIVE'S STEEPEST-DESCENT DIRECTION HAS COSINE -0.03 (SE 0.02) WITH THE DIRECTION TO THE NATIVE, 56/126 POSITIVE; -0.14 ON FAIL18; NO S27 CHANNEL DOES BETTER (2026-09-14 21:40, A2; re-appended 21:58 after a write race with S28-L23 lost the first copy)

EVERY COSINE HERE IS ORACLE (u is the direction from the production cloud to the native).
Pre-registered `s27/PREREG_S28_A.md` ADDENDUM 4, A2.1. Artefacts:
`s27/results/s28_A2_cosine_rows.jsonl` (126 rows), `s27/results/s28_A2_summary.json`
(`cosine`, `text`), job `s26/jobs_done/s28A2_cosine_126.json` (60 s, peak RSS 0.298 GB). Code
`s27/s28_A2_local.py`, tests `tests/test_s28_A2.py` (5 pass: analytic gradients vs finite
differences, rigid-body removal, the circuit step's displacement, NaN-poison of the ladder).

Setup: C0 = the DIS top-75 uniform average in its medoid frame (anchored 3.048338); g = dS~/dC
(the shipped Bayes risk read by linear interpolation, analytic); u = (native Kabsch-aligned
onto C0) - C0; rigid-body components (3 translations, 3 infinitesimal rotations about the
centroid) projected out of both; cos(-g, u) per target. RG_LAW and EXVOL analytic; DISTPOT,
CONTACT, ENV, CAGEO are histogram lookups (zero gradient almost everywhere) and get a smoothed
central difference at h = 0.5 A, labelled so; NaN where the smoothed gradient is zero (EXVOL has
active pairs on 16/126, CONTACT on 77, ENV on 92).

  DIS (analytic)
  all 126                                  mean -0.034  SE 0.021  median -0.043  positive 56/126  sign-test p 0.247
  FAIL18                                   mean -0.143  SE 0.066  median -0.164  positive 6/18  sign-test p 0.238
  other 108                                mean -0.016  SE 0.022  median -0.022  positive 50/108  sign-test p 0.501

  RG_LAW (analytic)
  all 126                                  mean +0.034  SE 0.033  median +0.021  positive 66/126  sign-test p 0.656
  FAIL18                                   mean +0.138  SE 0.143  median +0.286  positive 12/18  sign-test p 0.238
  other 108                                mean +0.017  SE 0.030  median +0.004  positive 54/108  sign-test p 1

  EXVOL (analytic)
  all 126                                  mean -0.035  SE 0.030  median -0.038  positive 7/16  sign-test p 0.804
  FAIL18                                   mean -0.089  SE 0.046  median -0.101  positive 1/4  sign-test p 0.625
  other 108                                mean -0.016  SE 0.036  median +0.006  positive 6/12  sign-test p 1

  DISTPOT (smoothed FD, step function)
  all 126                                  mean -0.002  SE 0.015  median -0.011  positive 60/126  sign-test p 0.656
  FAIL18                                   mean -0.052  SE 0.049  median -0.057  positive 6/18  sign-test p 0.238
  other 108                                mean +0.007  SE 0.015  median +0.001  positive 54/108  sign-test p 1

  CONTACT (smoothed FD, step function)
  all 126                                  mean +0.030  SE 0.017  median +0.048  positive 44/77  sign-test p 0.254
  FAIL18                                   mean +0.052  SE 0.044  median +0.065  positive 7/13  sign-test p 1
  other 108                                mean +0.026  SE 0.019  median +0.046  positive 37/64  sign-test p 0.26

  ENV (smoothed FD, step function)
  all 126                                  mean -0.006  SE 0.020  median -0.006  positive 45/92  sign-test p 0.917
  FAIL18                                   mean -0.082  SE 0.048  median -0.068  positive 4/15  sign-test p 0.118
  other 108                                mean +0.009  SE 0.021  median +0.017  positive 41/77  sign-test p 0.649

  CAGEO (smoothed FD, step function)
  all 126                                  mean +0.010  SE 0.013  median +0.014  positive 74/126  sign-test p 0.0609
  FAIL18                                   mean +0.043  SE 0.022  median +0.030  positive 15/18  sign-test p 0.00754
  other 108                                mean +0.004  SE 0.015  median +0.010  positive 59/108  sign-test p 0.387

  random-direction reference: mean |cos| of a random shape field with u = 0.140 (the scale of a meaningless cosine at 3n-6 dof)
  |grad S~| RMS per atom: mean 0.0488 A^-1; |u| RMS (distance to the native in the frame): mean 3.048 A
  Spearman(cos_DIS, production RMSD) = -0.372 (n=126)

Lane D's S28-L23 caveat (a), the cosine's null: 16 random shape fields per target (rigid-body
removed) give mean |cos| 0.140 with u over 126 (D's estimate sqrt(2/(pi(3n-6))) = 0.12 to
0.18); the ORACLE cosine of the distogram, -0.034 with SE 0.021, is inside that null in
magnitude and below it in sign. The per-target random |cos| values (16 per target) are in
`cos_random_ref` in every row, so the 95th percentile and the FAIL18 / 108 split of the null
are recomputable from the artefact.

Reading (ORACLE diagnostic):
1. The distogram's gradient at the pipeline's own output carries no directional information
   about the native: mean cosine -0.034 (SE 0.021), median -0.043, 56/126 positive (sign test
   p 0.25). Steepest descent on the shipped objective from production is, on average, a random
   direction with a slight lean AWAY from the native. This is the LOCAL statement behind
   S28-L18b's global one (the objective's minimiser is away from the native) and behind S15's
   3.321 and S8-9's 37th percentile: not only is the optimum in the wrong place, the first step
   is uninformative.
2. The regime split, as the prior said to check: on FAIL18 the cosine is -0.143 (SE 0.066,
   6/18 positive), on the 108 -0.016 (SE 0.022, 50/108). Where the distogram is wrong its
   gradient points away from the native; where it is right the gradient is orthogonal.
   Spearman(cos_DIS, production RMSD) = -0.372 over 126: the worse the target, the more the
   objective's descent direction opposes the native. That is the sequence-conditioning sign
   flip (memory: blind beats shipped on FAIL18) seen from the objective's side, and it is NOT
   "positive on the 108": the 108 are at zero, not above it.
3. No S27 channel is locally informative either: RG_LAW +0.034 (SE 0.033), DISTPOT -0.002,
   CONTACT +0.030, ENV -0.006, CAGEO +0.010 (SE 0.013; 74/126 positive, sign test p 0.06;
   15/18 on FAIL18, p 0.008, at a mean of +0.04, one of 21 sign tests in this table and a
   magnitude a third of the random reference: recorded, not a result). The distogram is the
   most informative ranker on this pool (S27) and its gradient is as blind as the rest.
4. What this predicts for A2.2 (deployable): a step of e A along -g should cost about what a
   random direction of the same size costs, and more on FAIL18; the prior stands.
Answers to the rest of S28-L23's caveats, before any A2.2 number: (b) any best e over the grid
of three is priced with `best_of_k_within`, for the step arm and the circuit arm; (c) on the
chain the step is paired against the mean of the SAME two random draws (0 and 1) that are
projected, said so in the entry, and a best-of-2 is an order statistic; (d) production is
re-projected in the same job; (e) the circuit arm's baseline C(theta_P) is quoted with its
residual RMS to C0 and its own RMSD to the native, and its step is read against that baseline.
A2.2 (point cloud, job `s28A2_ladder_126`) is running; its built chain runs after the primary
chain job (`s28A_chain_primary_1`, suspended by the governor at 64/126 at 21:28, resumed) has
landed and its verdict is posted.

## S28-L24 -- ADVERSARY CHECK OF S28-L23b (A2.1, the ORACLE cosine at the production point) AND SUITE STATUS 21:40 (2026-09-14 21:40, lane D)
Question: is S28-L23b an ORACLE diagnostic with its null, and does its regime reading rest on
anything the shared-referent floor could manufacture?
- Labelling: every cosine is ORACLE, said in the first line and in every reading. Correct.
- The null (S28-L23 caveat a): supplied, 16 random shape fields per target, mean |cos| 0.140,
  inside my estimate 0.12 to 0.18; the distogram's -0.034 (SE 0.021) is inside the null in
  magnitude. The per-target null values are in the rows (`cos_random_ref`), so the 95th
  percentile is recomputable. Correct.
- Multiple comparisons: 21 sign tests in the table; the one at p = 0.008 (CAGEO on FAIL18,
  mean +0.04, a third of the random reference) is flagged as one of 21 and not called a result.
  Correct.
- The shared-referent floor (project memory: two quantities measured against a common
  reference correlate by construction). Spearman(cos_DIS, production RMSD) = -0.372: both sides
  read the native. The floor here is zero by symmetry: the cosine is scale-free in u, and for a
  random shape field the cosine with u is independent of |u|, so a common referent alone gives
  no correlation between cos and RMSD; the -0.372 is a property of the distogram's gradient
  (it opposes u more where u is long). Say "not a floor artefact, by symmetry of a random
  direction" when it is quoted. It is still ORACLE and still a correlation of n = 126.
- The regime reading ("negative on FAIL18, zero on the 108, not positive on the 108") is what
  the numbers say (FAIL18 -0.143 SE 0.066, 6/18; other -0.016 SE 0.022, 50/108); the prior's
  "positive on the 108" half did not happen and the entry says so.
- The prediction for A2.2 (a step along -g costs what a random step costs, more on FAIL18) is
  registered before A2.2's numbers; I hold A2 to it.
Verdict: STANDS (an ORACLE diagnostic with its null; nothing deployable claimed).
SUITE STATUS 21:40: the six S28 test files (`tests/test_s28_A.py`, `test_s28_A2.py`,
`test_s28_B.py`, `test_s28_B2.py`, `test_s28_C.py`, `test_s28_D.py`) as one TEST job
`s28D_pytest_lanes_v2`: exit 0, 64 passed, 0 failed, 10 s, peak RSS 0.328 GB
(`s26/logs/s28D_pytest_lanes_v2.log`). Green gate: 286 pass / 3 skip (light files, S28-L4) +
64 pass (S28 files) = 350 pass / 3 skip / 0 fail on the files that can run; pipeline,
integration and the two AMBER files remain deferred (S28-L5).
Artefacts: `s27/results/s28_A2_cosine_rows.jsonl`, `s28_A2_summary.json`;
`s26/jobs_done/s28D_pytest_lanes_v2.json`.

## S28-L25 -- F5-B2: ON THE kNN SPREAD-SPECTRUM GRAPH THE HOPPING TERM'S GRADIENT VARIANCE DECAYS AT -1.0 PER QUBIT (GAUSSIAN: -1.7 TO -1.8) AND IS 30 TO 60x LARGER AT THE DEPLOYED WIDTH; THE PRE-REGISTERED BRIGHT LINE (-1.0) IS STRADDLED, NOT CLEARED; THE SECOND CLAUSE (HOPPING SHARE AT J = 1) IS RUNNING; THREE-WAY SPLIT MIDDLE LEG FOR S28B (2026-09-14 22:20, B2)

Question (`s27/PREREG_S28_B.md` addendum 1, F5-B2): was the -1.7 to -1.8 per qubit decay of the
hopping term's gradient variance (S28-L8b, S28-L11) the Gaussian graph's near-rank-one SPECTRUM?
Property measurement, no RMSD, no native. Job `s28B2_train` (killed once by the governor at
96.6% box RAM after 9 of 12 targets, re-queued by the governor itself, resumed from its
per-target checkpoint; `s26/jobs_done/s28B2_train.json` exit 0, 131 s for the last 3 targets,
peak RSS 0.33 GB; `s27/results/s28_B2_train_rows.jsonl` 720 rows, `s28_B2_train.json`).
Code `s27/s28_B2_knn.py` (6 tests, `tests/test_s28_B2.py`).

**The graph.** Symmetric binary kNN on the same pairwise CA-RMSD, degree-normalised
D^-1/2 A D^-1/2, unit spectral norm, k = 10 (arm) and k = 5 (replication). Spectral report,
median over the 12 trainability targets, n = 9 (the 500-window pool):
lambda_2 / lambda_1 = 0.982 (k 10), 0.991 (k 5) against the Gaussian graph's 0.138; one
connected component on every target; degrees 10..30 (k 10), 5..16 (k 5); the top eigenvector's
participation ratio 0.94 / 0.92 of dim and overlap with the uniform state 0.985 / 0.980 (the
Gaussian graph's 0.884 and 0.947); corr(degree, E) -0.27 / -0.18 (Gaussian -0.47 to -0.90). The
spectrum is spread, as designed, and the graph is much less a typicality projector.

**F5-B2, Var[dF/dtheta_0], the same draws, settings and targets as F5
(`s28_B2_train.json :: summary`; Gaussian row from `s28_B_train.json`):**

    cell                    n=4        n=5        n=6        n=7        n=8        n=9    slope 4..9   4..8
    k10 hop-only (J=1)  4.729e-03  1.981e-03  1.123e-03  5.183e-04  2.674e-04  1.224e-04   -1.033   -1.022
    k5  hop-only (J=1)  7.806e-03  4.107e-03  2.082e-03  1.054e-03  5.118e-04  2.478e-04   -0.997   -0.982
    Gaussian hop-only   4.061e-03  7.284e-04  2.817e-04  5.898e-05  3.948e-05  4.157e-06   -1.844   -1.700
    k10 full J=0        2.600e-02  3.174e-02  2.129e-02  2.092e-02  1.725e-02  3.051e-02   -0.043
    k10 full J=1        3.062e-02  3.331e-02  2.233e-02  2.135e-02  1.793e-02  3.070e-02   -0.078
    k10 full J=3        6.860e-02  4.695e-02  3.017e-02  2.521e-02  2.095e-02  3.190e-02   -0.265
    k5  full J=3        9.696e-02  6.542e-02  4.026e-02  2.938e-02  2.358e-02  3.264e-02   -0.364
    ratio full J=3 / J=0 at n = 9:  k10 1.046   k5 1.070   (Gaussian 1.005)

- The hop-only decay rate HALVES: -1.03 (k 10) / -1.00 (k 5) per qubit against -1.84 (-1.70 on
  the padding-free registers) for the Gaussian graph, and the per-qubit steps are regular
  (k 10: -1.26, -0.82, -1.12, -0.95, -1.13). At n = 9 the hopping term's gradient variance is
  1.22e-4 (k 10) / 2.48e-4 (k 5): 29x / 60x the Gaussian graph's 4.16e-6. Against the diagonal
  terms' 3.05e-2 the ratio is now 250x (k 10) / 123x (k 5) at J = 1 and 28x / 14x at J = 3
  (Gaussian: 7,300x / 820x).
- THE PRE-REGISTERED BRIGHT LINE IS STRADDLED, NOT CLEARED. Addendum 1's first clause was
  "shallower than -1.0 per qubit": k 5 reads -0.997 (4..9) / -0.982 (4..8), shallower by 0.003 /
  0.018; k 10 reads -1.033 / -1.022, steeper by 0.033 / 0.022. With S28-L11's +-0.3 fit
  uncertainty both are AT the line. What is measured cleanly is the change: the decay rate
  halved and the n = 9 variance rose 30 to 60x when lambda_2 / lambda_1 went from 0.14 to 0.98
  with everything else identical (same draws, same E, same targets, same circuit). That
  supports the spectrum reading of S28-L11 in direction and size; it does not clear the
  registered threshold, and I do not claim it does. No slope here is a plateau or its absence;
  a -1.0 per qubit decay on a spread-spectrum observable at depth 3 is reported as a number.
- The full objective's variance is still flat in n at every J (the CVaR and entropy terms
  dominate at n = 9: the hopping contributes 4.6% (k 10) / 7.0% (k 5) of the variance at J = 3).
- The second clause (the circuit collects more than half of its same-sign bound at J = 1) is
  measured by job `s28B2_share` (queued; native-free, no readout), and will be appended as its
  own entry with the Gaussian graph's share on the same targets beside it.

**The three-way split's middle leg for S28B (S28-L2(b); `s27/results/s28_B_split.json`, job
`s28B_split` + `s28B_split2` resumed from checkpoint, exit 0, peak RSS 0.441 GB; 126 targets x
13 cells):** F = CVaR_0.18 - 0.5 H - J <A> (the circuit's objective) evaluated at the VQE
optimum, at the exact ground state of H, and at the best of 16 untrained draws (draw 1 is the
VQE's own initial theta), means over 126, REAL graph:

    J     F_vqe s0   F_vqe s1   F_gs       F_untrained16   VQE < GS (frac, s0)   VQE < untrained   hop_vqe s0   hop_gs   coh s0
    0     -4.5610    -4.5561    -1.7283    -4.0629         1.000                 1.000             0.000        0.000    0.001
    0.1   -4.5608    -4.5565    -1.7757    -4.0627         1.000                 1.000            -0.001        0.002    0.001
    0.3   -4.5611    -4.5585    -2.0654    -4.0625         1.000                 1.000            -0.000        0.021    0.002
    1     -4.5663    -4.5854    -4.5779    -4.0615         0.468                 1.000             0.007        0.499    0.009
    3     -5.5089    -5.2652    -7.3114    -4.0588         0.286                 1.000             0.318        0.911    0.339
    (PERM and RAND within 0.1 of REAL on every column)

Reading: at J <= 0.3 the exact ground state is one-hot and the entropy term makes it a poor
point for the circuit's objective (the VQE state is lower on 126/126). At J = 1 the two tie
(-4.566 / -4.585 vs -4.578; VQE lower on 47%). At J = 3 the exact ground state (PR 305, 8.5
bits, hopping 0.911) has F = -7.31 against the circuit's -5.51 / -5.27: from the untrained
best-of-16 (-4.06) the circuit closes 45% (seed 0) / 37% (seed 1) of the gap to a state with far
lower F, and collects 0.32 of the 0.91 hopping the ground state collects, with sign coherence
0.34. The VQE beats its own best-of-16 untrained draws on 126/126 targets at every J (it
trains); it does not find the sign-coherent delocalised state at J = 3. Whether that state is
representable by the 27-parameter circuit is not measured; what is measured is that the
optimiser's reached F is 1.8 above it. This is the optimisation-quality leg for the S28B
verdict entry; the Hamiltonian-quality and emitted-structure legs are written there.

S28-L23 caveats: (b) every kNN graph on the 12 targets is CONNECTED at both k and every n
(`s28_B2_train.json :: spectral.n_components` = 1 throughout), so lambda_1 = 1 is
non-degenerate and no GS row needs the DEGENERATE label on this set; the count is carried per
target into any endpoint run. (a) the top eigenvector's uniform overlap is 0.985 (k 10) /
0.980 (k 5) at n = 9, so the typicality component is still present and is joined by a spread
remainder; the rank-one / remainder decomposition of the kNN hop-only variance (`s28_B2_rank1.py`,
job `s28B2_rank1`, queued) is reported in the next B2 entry so the reading is "the remainder
carries X% of the variance", as asked.

Verdict (property): the Gaussian graph's decay was the spectrum's (the rate halves and the
n = 9 variance rises 30 to 60x on a spread spectrum, same draws), but the registered -1.0 line
is straddled (k 5 -0.997, k 10 -1.033) and the first clause of F5-B2 is NOT cleanly passed;
the second clause is pending. The B2 endpoint stays gated on the S28B built-chain verdict and
lane D's check.
Artefacts: `s27/results/s28_B2_train.json`, `s28_B2_train_rows.jsonl`, `s28_B_split.json`,
`s28_B_split_rows.jsonl`, `s26/jobs_done/s28B2_train.json`, `s28B_split2.json`.

## S28-L26 -- ADVERSARY CHECK OF S28-L25 (F5-B2, the spread-spectrum graph's trainability; the S28B split's middle leg) (2026-09-14 21:55, lane D)
Question: do the numbers say what the entry says, and is the "straddled, not cleared" reading
the honest one? Recomputed from `s27/results/s28_B2_train.json :: summary` and
`s28_B_train.json`:
- Slopes: k5 hop-only -0.997 (4..9) / -0.982 (4..8), steps -0.93, -0.98, -0.98, -1.04, -1.05;
  k10 -1.033 / -1.022, steps -1.26, -0.82, -1.12, -0.95, -1.13; n = 9 variance 59.6x (k5) and
  29.4x (k10) the Gaussian graph's 4.157e-6. All as the entry states. The steps are regular
  (no n = 9 excess as in the Gaussian row), so the fit uncertainty here is smaller than the
  +-0.3 I attached to the Gaussian slope; the two k straddle -1.0 by 0.003 and 0.033. The
  registered clause was a bright line and the entry does not claim it; correct.
- What is cleanly measured: same draws, same E, same targets, same circuit, only the graph
  changed (lambda_2 / lambda_1 from 0.14 to 0.98), and the hop-only decay rate halved while
  the n = 9 variance rose 30 to 60x. That supports S28-L11's "the decay is the spectrum's" in
  direction and size, as the entry says. The full objective stays flat at every J (hopping is
  5 to 7% of the variance at J = 3): the hopping term is still a minor part of the gradient at
  the deployed width, so the endpoint prior (null or worse) is unchanged.
- The one thing the entry cannot yet say, and should not until the rank-one decomposition
  lands (job `s28B2_rank1`, queued): the kNN top eigenvector is MORE uniform than the
  Gaussian's (overlap 0.985 / 0.980 vs 0.947), so the typicality projector is still there;
  the 30 to 60x must come from the remainder. If the remainder carries most of the variance the
  reading is "a spread remainder is trainable where a rank-one projector is not"; if it does
  not, the reading is wrong. Hold "the spectrum was the mechanism" until that number posts.
- The middle leg (S28-L2(b)) for S28B, from `s28_B_split.json`: the VQE beats its own
  best-of-16 untrained draws on 126/126 at every J (it trains); at J <= 0.3 the one-hot ground
  state is a poor point for the circuit's OWN objective (entropy), so "VQE < GS on 126/126"
  there is the entropy term, not skill; at J = 1 they tie; at J = 3 the exact ground state
  (F -7.31, hopping 0.911, coherent) is 1.8 below the circuit's reached F (-5.51 / -5.27) and
  the circuit collects 0.32 of the 0.91. The entry says "whether that state is representable by
  the 27-parameter circuit is not measured". THAT IS THE MISSING PIECE OF THE LEG, and it is
  cheap and native-free: fit theta to maximise |<GS|psi(theta)>|^2 (Adam, 16 starts, the same
  loop) at J = 3 on the 12 trainability targets and report the best overlap and F at the
  fitted state. Overlap near 1 makes the 1.8 gap an OPTIMISATION shortfall; overlap well below
  1 makes it EXPRESSIVITY, and then the circuit is not failing to find a state it holds. Without
  it "optimisation quality" is a gap of unknown provenance. Requested before the S28B verdict
  entry quotes the leg.
- Rule 9: "a -1.0 per qubit decay ... is reported as a number": compliant; no plateau wording.
- Second clause (share at J = 1 > 0.5): pending (`s28B2_share`); the Gaussian's 1% at J = 1 is
  the comparator on the same targets, as registered.
Verdict: STANDS WITH CAVEAT (a property measurement; the bright line is straddled and the
entry says so; the mechanism reading waits for the rank-one decomposition; the middle leg needs
the representability fit before it is quoted as "optimisation quality").
Artefacts: `s27/results/s28_B2_train.json`, `s28_B_train.json`, `s28_B_split.json`; my refits
are `np.polyfit` on log2 of the per-n medians.

## S28-L26 -- VERDICT ON THE BUILT CHAIN (126/126): THE AMPLITUDE READOUT IS WORSE THAN PRODUCTION AT EVERY lam (+0.23 TO +0.26 A, 1.3x TO 1.6x MDE, FOLD CI ABOVE ZERO); F1 DOES NOT FIRE, F2 IS MOOT; THE ORACLE CEILING OF THE SAME FAMILY IS 0.252 A EMITTED; EXPRESSIVITY IS NOT THE BARRIER, RECOGNITION IS, AND THE OBJECTIVE IS THE FAILING LEG (2026-09-14 22:40, A)

Basis: BUILT CHAIN, `s12.instrument.project`, every arm and PRODUCTION projected by the same
code path in the same job (`s28A_chain_primary_1`, 4262 s, peak RSS 0.311 GB, resumed after
one governor suspension). Artefacts: `s27/results/s28_A_chain_rows.jsonl` (126 rows: prod,
oracle_circ, circ_l0/0.3/1/3_i80; oracle_aff500 on 22), `s27/results/s28_A_summary.json`
(`chain_means`, `contrasts`, `chain_lam_grid_bok`, `text`). Pre-registration
`s27/PREREG_S28_A.md` sections 3 to 5, addenda 1 to 3.

Anchor: production re-projected here is 3.2071 over 126 (point cloud 3.0483; projection price
+0.159). S27's DIS chain rows give 3.2126 on the same clouds (equal to 5.7e-14): mean
difference -0.0056, 12/126 targets differ by more than 0.02 A and 3 by more than 0.1
(2LNG -0.513, 9BFL -0.160, 2NDN -0.115): a branch flip of the multi-start projection under a
1e-14 input perturbation, lane D's S28-L18 floor with a longer tail than its 22-target sample
showed (max 0.51, not 0.02). It does not touch any contrast here (same code path both sides),
and it is why S27's rows were not reused as the comparator.

Means on the built chain (`chain_means`): production 3.2071; circuit lam 0.3 / 1 / 3 at 80
iterations 3.4432 / 3.4330 / 3.4708; lam 0 (the deployed state read as signed weights) 6.1366;
ORACLE circuit optimum 0.2516 (126/126 under 2 A, worst 0.970). Projection prices (chain minus
cloud): production +0.159, circuit lam arms +0.033 to +0.060, the ORACLE optimum -0.037 (the
projection moves a near-native signed combination TOWARD the native, S10-5's hull result seen
again), lam 0 -80.9 (a 87 A cloud projected onto an ideal chain lands at 6.1 A).

  circ_l0.3_i80 vs production (BUILT CHAIN)
    a 3.4432 (med 3.2660)   b 3.2071 (med 2.9661)   n=126
    effect +0.2362   median +0.1570   SE 0.0568   MDE 0.1591   effect/MDE +1.48
    iid  CI95 [+0.1292, +0.3461]
    fold CI95 [+0.1279, +0.3270]   folds same sign 5/5   per-fold 0:+0.314 1:+0.350 2:+0.310 3:+0.195 4:+0.054
    52W/74L/0T   worst degradation +2.6150 (9BAF)   p90 +0.9704   power 0.99  Type-M 1.01
    concentration: drop-top10 +0.3216 vs uniform-effect null p10/p50/p90 +0.2506/+0.3192/+0.3937 -> pctile 0.515
    VERDICT: WORSE

  circ_l1_i80 vs production (BUILT CHAIN)
    a 3.4330 (med 3.2831)   b 3.2071 (med 2.9661)   n=126
    effect +0.2260   median +0.1141   SE 0.0624   MDE 0.1749   effect/MDE +1.29
    iid  CI95 [+0.0982, +0.3490]
    fold CI95 [+0.0898, +0.3557]   folds same sign 4/5   per-fold 0:+0.206 1:+0.418 2:+0.349 3:+0.245 4:-0.022
    48W/78L/0T   worst degradation +2.6684 (5MXS)   p90 +1.1138   power 0.95  Type-M 1.03
    concentration: drop-top10 +0.3370 vs uniform-effect null p10/p50/p90 +0.2586/+0.3344/+0.4084 -> pctile 0.518
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.03x]

  circ_l3_i80 vs production (BUILT CHAIN)
    a 3.4708 (med 3.3342)   b 3.2071 (med 2.9661)   n=126
    effect +0.2638   median +0.2035   SE 0.0590   MDE 0.1654   effect/MDE +1.59
    iid  CI95 [+0.1495, +0.3757]
    fold CI95 [+0.1146, +0.3692]   folds same sign 5/5   per-fold 0:+0.346 1:+0.390 2:+0.361 3:+0.285 4:+0.001
    40W/86L/0T   worst degradation +2.8741 (2RUO)   p90 +1.1084   power 0.99  Type-M 1.00
    concentration: drop-top10 +0.3626 vs uniform-effect null p10/p50/p90 +0.2828/+0.3595/+0.4367 -> pctile 0.525
    VERDICT: WORSE

  circ_l0_i80 vs production (BUILT CHAIN)
    a 6.1366 (med 5.7725)   b 3.2071 (med 2.9661)   n=126
    effect +2.9295   median +2.9339   SE 0.2400   MDE 0.6724   effect/MDE +4.36
    iid  CI95 [+2.4729, +3.4112]
    fold CI95 [+2.6722, +3.1849]   folds same sign 5/5   per-fold 0:+2.820 1:+2.433 2:+3.320 3:+2.831 4:+3.152
    16W/110L/0T   worst degradation +11.8605 (8TXS)   p90 +5.9532   power 1.00  Type-M 1.00
    concentration: drop-top10 +3.3145 vs uniform-effect null p10/p50/p90 +2.9936/+3.3130/+3.6346 -> pctile 0.503
    VERDICT: WORSE

  circuit, lam chosen leave-fold-out, HELD-OUT vs production (BUILT CHAIN)
    a 3.4726 (med 3.2705)   b 3.2071 (med 2.9661)   n=126
    effect +0.2656   median +0.0970   SE 0.0615   MDE 0.1724   effect/MDE +1.54
    iid  CI95 [+0.1449, +0.3868]
    fold CI95 [+0.1479, +0.3743]   folds same sign 5/5   per-fold 0:+0.314 1:+0.418 2:+0.349 3:+0.245 4:+0.054
    50W/76L/0T   worst degradation +2.6684 (5MXS)   p90 +1.0648   power 0.99  Type-M 1.01
    concentration: drop-top10 +0.3687 vs uniform-effect null p10/p50/p90 +0.2895/+0.3674/+0.4467 -> pctile 0.507
    VERDICT: WORSE

  ORACLE oracle_circ vs production (BUILT CHAIN)
    a 0.2516 (med 0.2240)   b 3.2071 (med 2.9661)   n=126
    effect -2.9555   median -2.7272   SE 0.1450   MDE 0.4064   effect/MDE -7.27
    iid  CI95 [-3.2421, -2.6803]
    fold CI95 [-3.0662, -2.8496]   folds same sign 5/5   per-fold 0:-2.798 1:-2.940 2:-2.859 3:-3.174 4:-3.011
    126W/0L/0T   worst degradation -0.1218 (1S9Z)   p90 -0.9489   power 1.00  Type-M 1.00
    concentration: drop-top10 -2.6615 vs uniform-effect null p10/p50/p90 -2.8536/-2.6621/-2.4811 -> pctile 0.501
    VERDICT: BETTER

F1 ("the formulation moves accuracy"): DOES NOT FIRE. Every lam arm is WORSE than production
on the built chain with the fold CI above zero (lam 0.3: +0.236, 1.48x MDE, 5/5; lam 1: +0.226,
1.29x, 4/5, Type-M zone; lam 3: +0.264, 1.59x, 5/5). The nested leave-fold-out choice (inner
basis: the built chain; per-fold choices 0.3 / 1 / 1 / 1 / 0.3) is +0.266 at 1.54x MDE, 5/5,
WORSE; the per-target minimum over the lam grid is 90% accounted for by its order-statistic
null, k_eff 3.00, split-half transfer -7%: NOT A SIGNAL (`chain_lam_grid_bok`). No seed-1 run
is owed (no positive). F2 ("the circuit matters") is MOOT; its point-cloud preview stands as
S28-L18b item 4 with lane D's S28-L20 caveats answered here: (a) the circuit's advantage over
the unconstrained a500 family is -0.19 at matched budget (80 evaluations) and -0.30 converged
(2000), the budgets differ 25x and both are quoted; (b) the circuit-vs-random-27-subspace
contrast is Type-M (1.13x MDE, 4/5 folds) and its SIGN is that the circuit is the less bad of
the two signed families, nothing more.

The coordinator's three questions for a positive, answered for the negative:
(a) Through the projection: the harm persists (point cloud +0.335 / +0.337 / +0.390, chain
    +0.236 / +0.226 / +0.264; the projection removes 0.1 A of it because the signed clouds are
    uncontracted, mean virtual bond 3.40 to 3.54 before projection against production's 2.96,
    so production pays +0.159 for the projection and the signed arms +0.03 to +0.06). FAIL18 /
    108 on the chain (`ST.compare` per stratum; artefact rows): lam 0.3 is -0.073 on FAIL18
    (SE 0.163, 0.16x MDE, 10W/8L, NOT MEASURED) and +0.288 on the 108 (1.73x MDE, fold CI
    [+0.224, +0.352], 5/5, WORSE); lam 1: -0.073 (0.13x) / +0.276 (1.54x, WORSE); lam 3:
    -0.026 (0.05x) / +0.312 (1.84x, WORSE). The harm is carried by the 108 targets where the
    distogram is right; on the 18 where it is wrong the signed readout is a coin toss. S27 L9's
    regime pattern, again, and again with no native-free switch.
(b) New relative to S26/S27: the first SIGNED (non-convex) consumer of the pool on record
    (S10-5's affine bound was ORACLE only; every S26/S27 readout was convex). It goes beyond
    S27 section 6 (ranking information is anti-useful on an averaging readout) by showing that
    the readout's convexity was not the binding constraint: given sign freedom, the objective
    uses it to move AWAY from the native (S28-L18b item 2: S~ at the circuit optimum 1.34 <
    production 1.67 < the ORACLE structure 2.09). And the family's ORACLE ceiling, 0.252 A
    emitted on 126/126, is 2.955 A below production: EXPRESSIVITY IS NOT THE BARRIER.
(c) The three-way split: (i) OBJECTIVE quality: the failing leg; the native-free S~ prefers
    the circuit's optimum to production and production to a 0.29 A structure on 80% of
    targets. (ii) OPTIMISATION quality: the optimiser works (lam = 0 reproduces
    `run_cvar_vqe` bit-for-bit; S~ falls from 3.4 to 1.3 in 80 iterations; 400 iterations
    lower it further and make the structure worse, +0.483 point cloud). (iii) The EMITTED
    structure: worse at every lam on the chain, above. A positive here would not have been a
    quantum result in any case: the convex simplex under the same objective lands at
    production (3.056 point cloud), and every signed family is worse in proportion to how
    well it minimises S~.

What the objective would have to know (PREREG section 6, the brief's last question): the
ORACLE structure exists in the family on every target and survives the projection; to pick
it, an objective would have to score a 0.29 A signed combination below the contracted average
and below the objective's own signed optimum. The shipped Bayes risk does the opposite on
101/126 targets (S28-L1b). The objective's own error is the common-mode error of the pool
(S23 L9): a per-pair posterior median that the average already sits at.

Verdict on the built chain: REFUTED for accuracy (F1 silent, every arm worse beyond its MDE
with the fold CI above zero except lam 1 in the Type-M zone, same sign); the ORACLE ceiling
(0.252 A emitted, 126/126 under 2 A) is a DIAGNOSTIC and the sprint's mechanism statement is
S28-L18b item 2 plus S28-L23b (the gradient at production is blind). Lane D attacks this
entry next; nothing is built on it. The A2 chain job (`s28A2_chain`, 13 projections per
target) is queued behind the governor and reports in its own entry.

## S28-L27 -- NUMBERING CORRECTION: lane A's entry headed "S28-L26 -- VERDICT ON THE BUILT CHAIN ..." is S28-L26b (2026-09-14 22:45, A)

It collided with lane D's S28-L26 (the check of S28-L25), which was appended between A's read
of the tail and A's append. Content unchanged. From here on lane A reads the tail and appends
in one process.

## S28-L27 -- ADVERSARY CHECK OF LANE A's BUILT-CHAIN VERDICT (the entry headed "S28-L26 -- VERDICT ON THE BUILT CHAIN", a number collision with my S28-L26; read as S28-L26b) AND A SCOPE CORRECTION OF MY S28-L18 (2026-09-14 22:12, lane D)
Question: does the refutation stand, and is anything in it a hidden positive? Every number
recomputed from `s27/results/s28_A_chain_rows.jsonl` (126 rows) with `ST.compare`,
`s12.instrument.FAIL18` and `ST.best_of_k_within`.
- Contrasts on the built chain, production re-projected in the same job: circ lam 0.3 +0.2362
  (SE 0.0568, MDE 0.1591, 1.48x, fold CI [+0.128, +0.327], 5/5, 52W/74L), lam 1 +0.2260 (1.29x,
  4/5, Type-M zone), lam 3 +0.2638 (1.59x, 5/5): identical to the entry to four decimals.
- FAIL18 / 108 (mine agrees with the entry): lam 0.3 -0.073 on FAIL18 (0.16x, 10W/8L, NOT
  MEASURED) and +0.288 on the 108 (1.73x, fold CI [+0.224, +0.352], WORSE); lam 1 -0.073 /
  +0.276 (1.54x); lam 3 -0.026 / +0.312 (1.84x). The harm is carried by the 108; the 18 are a
  coin toss. S27 L9's regime pattern, no switch.
- The lam grid as an order statistic: observed -0.303, null -0.273, 90% accounted, k_eff 3.00,
  split-half -7%: NOT A SIGNAL. The nested held-out arm is WORSE (+0.266, 1.54x, 5/5). Correct.
- Leakage / poison: the chain phase projects stored clouds and reads the native only in the
  scorer (S28-L10 item 3; `tests/test_s28_D.py`); the recognition phase is native-free end to
  end (S28-L10 item 2). The ORACLE arm is labelled in every line it appears.
- Cosmetic-variant test: not a re-parameterisation of the convex average (S28-L20); new
  relative to S26/S27; the entry's (b) is the right statement.
- Three-way split: the entry's (c) matches S28-L20's reading; the failing leg is the objective.
- Type-M: lam 1 at 1.29x is in the zone with the same sign as the other two at 1.48x and
  1.59x; the entry flags it. Say "worse at every lam; the size at lam 1 is in the Type-M zone".
Verdict: STANDS (REFUTED for accuracy; F1 silent; F2 moot; the ORACLE ceiling 0.252 A emitted
is a diagnostic and labelled so). Nothing positive is hidden in the strata, the grid or the
controls. No seed-1 run is owed.

SCOPE CORRECTION OF MY S28-L18 (the built chain's numerical floor). S28-L18 measured the
production cloud re-projected by lane A against S27's chain rows on the 22 targets then
available: max |diff| 0.0186, mean 0.0025. On all 126 (`s28_A_chain_rows.jsonl :: prod` vs
`chain_rows.jsonl :: DIS`, clouds equal to 5.7e-14): mean difference -0.0056, 12 targets above
0.02, 3 above 0.1, and 2LNG at -0.513, 9BFL -0.160, 2NDN -0.115. The floor's TAIL is 0.5 A on
one target in 126, not 0.02; the mean of 126 moves by 0.006 A (3.2071 vs 3.2126). The reading
of S28-L18 stands and sharpens: the multi-start projection's branch choice flips under
floating-point input differences, and on a few targets the branches are 0.1 to 0.5 A apart
(S26 L88's degeneracy seen from the reproducibility side). Consequences: (a) any built-chain
contrast whose two sides come from different code paths carries a per-target floor with a
0.5 A tail and a mean floor near 0.006 A; (b) the production anchor itself is 3.2126 (S27's
rows) or 3.2071 (lane A's re-projection) depending on the path, and an entry must say which;
(c) both are inside every MDE quoted in this sprint, so no verdict changes. Recorded as
`s27/RETRACTIONS_S28.md` R2 (my own scope correction, 22 -> 126 targets).
Artefacts: `s27/results/s28_A_chain_rows.jsonl`, `s27/results/chain_rows.jsonl`; the per-target
tail is in this entry.

## S28-L28 -- NUMBERING CORRECTION: lane D's entry headed "S28-L27 -- ADVERSARY CHECK OF LANE A's BUILT-CHAIN VERDICT ..." is S28-L27b (2026-09-14 22:16, lane D)
It was appended seconds after lane A's "S28-L27 -- NUMBERING CORRECTION" (both read the tail
at S28-L26b). Content unchanged; the ledger is append-only so the heading is not edited. Lane
A's chain verdict is S28-L26b; my check of it is S28-L27b; my S28-L18 scope correction
(RETRACTIONS_S28 R2) is in S28-L27b.
