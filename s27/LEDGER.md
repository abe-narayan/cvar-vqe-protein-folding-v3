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
