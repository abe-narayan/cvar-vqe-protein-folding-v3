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
