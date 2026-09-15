# S28 LANE C FINDINGS -- WHERE THE ERROR LIVES: FAIL18 DETECTION ON A NEW FEATURE CLASS, AND READOUTS THAT CONSUME RANKING INFORMATION

Lane C, Sprint 28, 2026-09-14. Pre-registration `s27/PREREG_S28_C.md` (base + addenda 1 to 3,
all written before the number they govern). Code `s27/s28_C_fail18.py`, `s27/s28_C_readout.py`;
tests `tests/test_s28_C.py` (11, all pass). Ledger entries S28-L6 (Part 1), S28-L7/L7b (reproduction
check), S28-L15 (Part 2, point cloud, intermediate), and the built-chain entry that follows it. Every number below carries its artefact path. Bases are
named on both sides of every contrast; the BUILT CHAIN (`s12.instrument.project`) decides.

## 0. One paragraph

Part 1: the pool's own statistical-potential distribution (DISTPOT / ENV / CONTACT spread, skew,
rank agreement and top-75 overlap with the distogram, the top-75's mean rank under each
potential, the pool's consensus spread and contact entropy: 21 features no router had seen) does
not separate FAIL18 from the other 108 targets. Held-out AUROC of the nested ridge-logistic
classifier is 0.468 on the new class, 0.608 with the control class added, both inside their
500-draw label-permutation nulls (p 0.57 and 0.10); the S22 router feature set through the same
harness is 0.591 (p 0.19). The best single feature (DISTPOT_spread, 0.706) is inside the
max-over-26-singles null (p 0.076) and is 43% a length proxy. The switched arm was therefore
not run. Independently of any detector, the ORACLE switch (route exactly the 18 FAIL18 targets
to DIS+DISTPOT or DIS+ENV) moves the built chain by -0.020 A (0.28x its MDE) and -0.047 A
(0.62x): a perfect detector would be underpowered on this instrument. Part 1 is closed on both
grounds. Part 2: three readouts that keep the production top-75 and consume CONS's (or
DISTPOT's) ranking through the weight vector (medoid-plus-k-neighbours, a ranker-based trim, a
diversity-preserving weighted average) are all null-to-WORSE than the uniform average on the
point cloud, and worse than their permuted-ranker controls wherever the contrast is
measurable; the built-chain verdicts are in section 3 (job `s28C_readout_chain`).

## 1. DEMONSTRATED

### 1.1 The FAIL18 detector on the new feature class is null at the AUROC step (S28-L6)
`s27/results/s28_C_fail18.json` (job `s28C_fail18_run`, exit 0, 602 s, peak RSS 0.331 GB).
Label: `s12.instrument.FAIL18` (ORACLE; the 18 zero-recall targets), used only as the nested
training label. Features: `s27/results/s28_C_features.json` (126 rows, `complete: true`,
native-free; the NaN-poison test `test_features_nan_poison_bit_identical` asserts bit-identical
features with `nat_ca` and `oracle_rr` set to NaN).

| block | p | held-out AUROC | per-fold mean | null p50 | null p95 | p_perm |
|---|---:|---:|---:|---:|---:|---:|
| SP (new class) | 21 | 0.468 | 0.468 | 0.484 | 0.629 | 0.572 |
| CTRL (n, dg entropy, sim entropy, sim gap, top-75 spread) | 5 | 0.492 | 0.507 | 0.494 | 0.638 | 0.508 |
| SP+CTRL | 26 | 0.608 | 0.536 | 0.489 | 0.635 | 0.102 |
| old_S22 (S22 router set, comparison) | 15 | 0.591 | 0.640 | 0.507 | 0.661 | 0.190 |
| s26_new_all (S26 C4 blocks, comparison) | 25 | 0.531 | 0.547 | 0.494 | 0.631 | 0.354 |

F1 does not fire on any block. Balanced accuracy at the 0.5 decision level: 0.435 (SP), 0.495
(SP+CTRL); at the prevalence-matched threshold the SP block flags 18 targets of which 1 is
FAIL18, SP+CTRL flags 24 of which 4 are.

Singles (sign nested; own null): DISTPOT_spread 0.706 (p 0.006), cons_mean 0.686 (p 0.010),
ENV_rho_dis 0.659 (p 0.042), n 0.654 (p 0.026), ENV_skew 0.648 (p 0.042). Priced as an order
statistic over the 26 singles (max-AUROC null, 500 draws): max-null mean 0.649, p95 0.720, p_max
0.076 for the best single. Length diagnostic (native-free): DISTPOT_spread has corr -0.43 with
n and falls to 0.593 after residualising on n; cons_mean (corr +0.75) falls to 0.406.

### 1.2 The ORACLE switch ceiling is below its MDE on the built chain (S28-L6)
`s27/results/s28_C_fail18.json :: oracle_switch`, per-target endpoints from
`s27/results/chain_rows.jsonl` (reproduction check: section 1.3). Switching exactly the 18
FAIL18 targets (ORACLE) to the alternative's top-75, everything else production:

| ORACLE switch | basis | effect | SE | MDE | effect/MDE | fold CI | FAIL18 stratum (n=18) |
|---|---|---:|---:|---:|---:|---|---:|
| FAIL18 -> DIS+DISTPOT | built chain | -0.0201 | 0.0254 | 0.0712 | -0.28 | [-0.059, +0.009] | -0.141 (SE 0.180), 9W/9L |
| FAIL18 -> DIS+ENV | built chain | -0.0474 | 0.0271 | 0.0760 | -0.62 | [-0.123, +0.007] | -0.332 (SE 0.180), 11W/7L |
| FAIL18 -> DIS+DISTPOT | point cloud | -0.0273 | 0.0244 | 0.0682 | -0.40 | [-0.068, +0.004] | -0.191 (SE 0.169) |
| FAIL18 -> DIS+ENV | point cloud | -0.0572 | 0.0292 | 0.0819 | -0.70 | [-0.136, -0.004] | -0.401 (SE 0.189) |

The point-cloud strata reproduce S27 L9 (`s27/results/strata.json`: -0.191 / -0.401). Through
the projection the per-FAIL18-target gain shrinks by 0.05 to 0.07 A and the mean effect of a
PERFECT switch is 0.28x to 0.62x its own MDE. A detector that is 50% accurate on FAIL18 (and
harms the false positives by the non-FAIL18 stratum of S27 L9, +0.010 / +0.064 per target)
would be worth 0.01 to 0.02 A on the mean.

### 1.3 The production projection reproduces S27's chain rows (S28-L7)
`s27/results/s28_C_reproduce.json` (job `s28C_reproduce`, exit 0, 60 s, peak RSS 0.322 GB): the
DIS and DIS+DISTPOT top-75 averages of the first six targets (1A13, 1A1P, 1CB3, 1CEK, 1CS9,
1D0W), re-emitted through `readout_uniform` + `readout_projected`, reproduce `chain_rows.jsonl`
with a maximum absolute deviation of 0.0 on all twelve point-cloud and twelve built-chain values
(S28-L7 prints them). The readout job's PROD arm reproduces `chain_rows.jsonl :: DIS` on 1A13 to
the last printed digit (2.6158612331525704 / 2.635339216442019 on both files,
`s27/results/s28_C_readout_chain_rows.jsonl`). Noise floor: the identity cell DIVW(0, 0) differs
from PROD by 8e-15 A on the point cloud and 1.4e-5 A on the built chain (the multi-start
projection amplifies floating-point input differences by about 1e9); built-chain contrasts
therefore carry a numerical floor near 1e-5 A, four orders below any MDE here.

### 1.4 Ranking-consuming readouts on the point cloud (intermediate basis, S28-L15; the built-chain entry decides)
`s27/results/s28_C_readout_cloud_rows.jsonl` (job `s28C_readout_cloud`, 29 arms x 126 targets)
and `s27/results/s28_C_readout_cloud_summary.json`. Every arm keeps the production top-75 and
changes only the weight vector (addendum 2). Production point cloud 3.0483; random-75 null 3.4209.

| arm (ranker CONS unless stated) | mean | effect vs PROD | effect/MDE | fold CI | verdict | vs permuted-ranker control |
|---|---:|---:|---:|---|---|---|
| MEDNB k=10 (11-member average) | 3.3226 | +0.274 | 1.63 | [+0.166, +0.378] | WORSE 5/5 | +0.047 (0.25x) |
| MEDNB k=20 | 3.2458 | +0.197 | 1.41 | [+0.091, +0.291] | WORSE | +0.068 (0.41x) |
| MEDNB k=40 | 3.1501 | +0.102 | 1.06 | [+0.035, +0.156] | WORSE, Type-M | +0.041 (0.37x) |
| TRIM q=0.05 (4 dropped) | 3.0628 | +0.015 | 0.79 | [+0.007, +0.025] | not measured (Type-M zone) | +0.011, fold CI [+0.001, +0.023] |
| TRIM q=0.10 (8 dropped) | 3.0838 | +0.035 | 1.09 | [+0.018, +0.049] | WORSE, Type-M | +0.031, fold CI [+0.015, +0.048] (0.97x) |
| TRIM q=0.20 (15 dropped) | 3.1048 | +0.057 | 1.02 | [+0.020, +0.081] | WORSE, Type-M | +0.054, fold CI [+0.018, +0.081] (0.96x) |
| TRIM[DISTPOT] q=0.05 / 0.10 / 0.20 | 3.053 / 3.053 / 3.050 | +0.005 / +0.005 / +0.002 | 0.34 / 0.21 / 0.07 | all spanning zero | null | +0.004 / +0.003 / +0.002 |
| DIVW b=1, g=0 (ranker weights only) | 3.1374 | +0.089 | 0.95 | [+0.020, +0.145] | not measured (0.95x) | +0.079, fold CI [+0.012, +0.133] |
| DIVW b=0, g=1 (repulsion only, no ranker) | 3.0557 | +0.007 | 0.14 | [-0.029, +0.049] | null | (control is PROD) |
| DIVW b=0.5, g=1 | 3.0823 | +0.034 | 0.54 | [-0.031, +0.073] | null | +0.018 (0.30x) |
| DIVW b=1, g=1 | 3.1297 | +0.081 | 0.84 | [-0.004, +0.144] | not measured | +0.054 (0.50x) |
| DIVW b=2, g=1 | 3.2096 | +0.161 | 1.16 | [+0.063, +0.251] | WORSE, Type-M | +0.112, fold CI [+0.028, +0.194] |
| DIVW b=0, g=0 (identity) | 3.0483 | -0.0000 (max 1e-13) | | | reproduces production | |

The permuted-ranker controls (same operator, CONS rank-permuted within the 75) are all less
harmful than the real-ranker arms: MEDNB~perm +0.227 / +0.129 / +0.060, TRIM~perm +0.003 /
+0.005 / +0.003, DIVW~perm +0.011 / +0.016 / +0.027 / +0.049. Grid pricing
(`ST.best_of_k_within` on each family's (target x cell) matrix including PROD): the
transferable choice is the identity in every family (MEDNB split-half -0.143 of a -0.274
oracle, k_eff 3.0; TRIM[CONS] -0.026 of -0.075; TRIM[DISTPOT] "NOT A SIGNAL", -8%; DIVW "NOT A
SIGNAL", 24%). Strata: every CONS-informed arm is negative on FAIL18 (-0.06 to -0.17, SE 0.09
to 0.15, none beyond 1.4 SE) and positive on the 108 (+0.04 to +0.34), the S27 L9 regime
pattern in a readout rather than a Hamiltonian.

## 2. ORACLE DIAGNOSTIC
- The ORACLE switch ceiling (1.2). It reads FAIL18 membership; it is the prize, never a result.
- The FAIL18 / non-FAIL18 strata of every readout arm (1.4, section 3).

## 3. Built chain (the reporting basis)
[Filled from `s27/results/s28_C_readout_chain_rows.jsonl` and
`s27/results/s28_C_readout_chain_summary.json` when job `s28C_readout_chain` lands.]

## 4. REFUTED
- H(Part 1): "the pool's statistical-potential distribution separates FAIL18". AUROC inside the
  permutation null on every block; the best single inside the max-over-26 null; the registered
  prior (null) confirmed. The eighth router construction on record lands where the seven before
  it did (S22 L7, S23 L7, S26 L110/L115), now on the FAIL18 label itself.
- H(Part 2a): "the CONS medoid of the top-75 plus its k nearest neighbours beats the uniform
  average". WORSE at every k on the point cloud, beyond MDE, 4/5 to 5/5 folds. (Built chain:
  section 3.)
- H(Part 2b): "a ranker-based outlier trim preserves the average's variance reduction". Trimming
  the CONS-worst members is WORSE than trimming random members (fold CI above zero at q = 0.05,
  0.10, 0.20): the CONS outliers are useful to the average. Trimming by DISTPOT is
  indistinguishable from a random trim. (Built chain: section 3.)
- H(Part 2c): "a diversity-preserving weighted average consumes ranking information". Ranker
  weights alone (b = 1, g = 0; ESS 40.7 of 75) cost +0.089 (0.95x MDE, fold CI above zero, the
  Type-M zone) and are worse than permuted weights of the same magnitudes (same ESS 40.7:
  +0.079, 0.82x MDE, fold CI [+0.012, +0.133]); the repulsion term alone is null (+0.007,
  0.14x); adding the repulsion to the ranker weights does not rescue them (b = 1, g = 1:
  +0.081, 0.84x). At b = 2 the cell is WORSE beyond its MDE (+0.161, 1.16x) and worse than its
  control (+0.112, fold CI [+0.028, +0.193]). (Built chain: section 3.)

## 5. HYPOTHESIS
- None advanced. The one direction the data points to (ENV_rho_dis, 0.659 held-out, unchanged
  after residualising on length) is one of 26 singles and inside the max null; it is recorded so
  it is not rediscovered, not proposed.

## 6. OPEN
- Whether a FAIL18 detector could matter on a longer instrument: the ORACLE ceiling (-0.047 A
  for DIS+ENV) would need n of roughly 126 x (0.076 / 0.047)^2 = 330 targets to be a measurable
  perfect-switch effect, and a real detector would be well below that. No fresh benchmark exists
  (`no-fresh-benchmark-exists`).

## 7. What damaged my own expectations
- I expected the 26-feature block to do at least as well as its best single; the ridge logistic
  on 18 positives with alpha chosen by inner AUROC is worse than a one-feature rule (0.608 vs
  0.706), and the single itself is inside the max-over-26 null. The capacity that S22 L10 said
  is not learnable at n = 126 is not learnable here either.
- I expected the ORACLE switch to be a clean, measurable ceiling that the detector would then
  fail to reach. It is not measurable itself on the built chain (0.28x to 0.62x MDE). The
  question "can FAIL18 be detected" was the wrong bottleneck: the prize is below the instrument's
  resolution once the projection is applied.
- I expected the ranker-based trim (b) to be null. Trimming the CONS-worst 8 of 75 costs +0.035
  on the point cloud and is worse than a random trim of the same size with the fold CI above
  zero: S27 section 6's mechanism holds at the level of 8 members, and the "outliers" a
  consensus ranker names are the members whose removal moves the average toward the mode.

## 8. What I did not do and why
- The switched arm (prereg 1.6): gated on F1 by addendum 2; F1 did not fire. Its ORACLE ceiling
  is reported instead.
- A second permutation seed and a reversed fold order: the positive-only replication rule; there
  is no positive.
- The DIVW cell (1, 0.5) named in prereg 2.3: replaced by (1, 1) under addendum 2 (one beta
  ladder at gamma = 1 plus the two one-factor references and the identity); the code's
  `PRIMARY_CHAIN` matches addendum 2.
- Any AMBER: none in this lane. The S26 `moved` signal (partial rho +0.43 with rmsd_arm, S26
  L1863) is the one native-free FAIL18-adjacent quantity not in the feature table; it needs an
  OpenMM relaxation per target and was outside the brief's feature list.
- Extending the grid beyond addendum 2 (k, q, beta): every cell's transferable choice is the
  identity; more cells would be cosmetic variants of a weighted convex average (coordinator
  steer, addendum 2 item 3).
