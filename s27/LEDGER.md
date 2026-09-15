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

