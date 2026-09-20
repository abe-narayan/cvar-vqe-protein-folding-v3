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
measurable; on the built chain (the reporting basis, S28-L32) every cell is on the harmful side at
0.29x to 0.99x its MDE, none clears it, and the two ranker trims are worse than random trims with the
fold CI above zero. Part 2 is refuted as registered.

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
from PROD by 8e-15 A on the point cloud and 1.4e-5 A on the built chain on 1A13; lane D's
S28-L18 measured the same floor across 22 targets at up to 0.02 A per target (0.003 A mean)
when the multi-start projection flips a near-tied branch, and supersedes S28-L7's "1e-5"
(`s27/RETRACTIONS_S28.md` R1, a scope correction). It does not touch the contrasts here: every
readout arm and its PROD comparator are projected in the same job through the same code path
(`s28_C_readout.py :: oracle_rows`), and PROD reproduces `chain_rows.jsonl :: DIS` bit-exactly
on every target it has processed (max |diff| 0.0 on 49/49 at 20:45; the final count is in the
built-chain entry).

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

## 3. Built chain (the reporting basis; S28-L32)
`s27/results/s28_C_readout_chain_rows.jsonl` (jobs `s28C_readout_chain` to `_chain4`; three governor
kills at 15, 25 and 25 targets under box-wide user load, resumed from the per-target checkpoint; 1,513
rows, no duplicate (pdb, arm)) and `s27/results/s28_C_readout_chain_summary.json`. Both sides of every
contrast are projected in the same job through the same code path; PROD reproduces
`chain_rows.jsonl :: DIS` bit-exactly on 126/126. Production built chain 3.2126.

| cell | mean | effect vs PROD | effect/MDE | fold CI | folds | its control vs PROD | cell vs control | FAIL18 (n=18) | other 108 |
|---|---:|---:|---:|---|---|---|---|---:|---:|
| MEDNB[CONS,k=20] | 3.3613 | +0.1487 | 0.99 | [+0.0342, +0.2367] | 4/5 | +0.1072 (0.80x) | +0.0415 (0.24x, fold [-0.062, +0.134]) | -0.117 (SE 0.123) | +0.193 (SE 0.058) |
| TRIM[CONS,q=0.1] | 3.2434 | +0.0308 | 0.75 | [+0.0108, +0.0461] | 4/5 | +0.0038 (0.17x) | +0.0270 (0.66x, fold [+0.019, +0.035]) | -0.016 (SE 0.033) | +0.039 (SE 0.016) |
| TRIM[DISTPOT,q=0.1] | 3.2343 | +0.0216 | 0.62 | [+0.0029, +0.0397] | 4/5 | +0.0014 (0.06x) | +0.0203 (0.52x, fold [+0.010, +0.036]) | -0.010 (SE 0.021) | +0.027 (SE 0.014) |
| DIVW[CONS,b=1,g=1] | 3.2833 | +0.0707 | 0.66 | [-0.0237, +0.1426] | 4/5 | +0.0517 (0.80x) | +0.0190 (0.16x, fold [-0.075, +0.086]) | -0.146 (SE 0.098) | +0.107 (SE 0.041) |
| DIVW[CONS,b=1,g=0] | 3.2862 | +0.0736 | 0.70 | [-0.0174, +0.1427] | 4/5 | +0.0188 (0.47x) | +0.0548 (0.53x, fold [-0.048, +0.123]) | -0.152 (SE 0.093) | +0.111 (SE 0.040) |
| DIVW[CONS,b=0,g=1] | 3.2278 | +0.0152 | 0.29 | [-0.0112, +0.0354] | 4/5 | PROD | n/a | -0.009 (SE 0.027) | +0.019 (SE 0.021) |

Every real-ranker cell is on the harmful side and none clears its own MDE (0.29x to 0.99x): the
verdict is NOT MEASURED in the harmful direction on all six, with the point-cloud sign preserved on
6/6 and the magnitude shrunk by 0.01 to 0.05 A through the projection. The permuted-ranker controls
are less harmful than the real cells in every pair; the two ranker trims are worse than a random trim
of the same size with the fold CI above zero (CONS +0.027 at 0.66x MDE; DISTPOT +0.020 at 0.52x): by
S28-L16's rule the sign is measured and the size is not. The repulsion-only cell is null (+0.015,
0.29x). Grid pricing on the chain for the DIVW family ((1,0), (0,1), (1,1), PROD):
oracle-min -0.1760, split-half -0.0284 (16%), k_eff 3.9, NOT A SIGNAL (split-half transfers 16% of the oracle); the k and q grids are
priced on the point cloud (section 1.4), where the transferring cell is the identity.
Strata (ORACLE label): the CONS-informed cells are negative on FAIL18 (-0.12 to -0.15, at most 1.7
SE) and positive on the 108 (+0.04 to +0.19, 2.4 to 3.3 SE); the controls are positive on both.

Verdict: the Part 2 falsifier fails on every cell. A ranker-informed convex weight vector over the
production top-75 does not beat the uniform one on this instrument, and where the contrast is
resolvable it is worse than a random weight vector of the same shape.

## 4. REFUTED
- H(Part 1): "the pool's statistical-potential distribution separates FAIL18". AUROC inside the
  permutation null on every block; the best single inside the max-over-26 null; the registered
  prior (null) confirmed. The eighth router construction on record lands where the seven before
  it did (S22 L7, S23 L7, S26 L110/L115), now on the FAIL18 label itself.
- H(Part 2a): "the CONS medoid of the top-75 plus its k nearest neighbours beats the uniform
  average". WORSE at every k on the point cloud, beyond MDE, 4/5 to 5/5 folds; on the built chain k = 20 is
  +0.149 at 0.99x MDE (fold CI [+0.034, +0.237]), harmful and not measured.
- H(Part 2b): "a ranker-based outlier trim preserves the average's variance reduction". Trimming
  the CONS-worst members is WORSE than trimming random members (fold CI above zero at q = 0.05,
  0.10, 0.20): the CONS outliers are useful to the average. Trimming by DISTPOT is
  indistinguishable from a random trim on the point cloud and +0.020 worse than one on the built chain
  (0.52x MDE, fold CI above zero: sign, not size). Built chain for the CONS trim: +0.031 vs PROD
  (0.75x), +0.027 vs a random trim (0.66x, fold CI [+0.019, +0.035], 5/5 folds).
- H(Part 2c): "a diversity-preserving weighted average consumes ranking information". Ranker
  weights alone (b = 1, g = 0; ESS 40.7 of 75) cost +0.089 (0.95x MDE, fold CI above zero, the
  Type-M zone) and are worse than permuted weights of the same magnitudes (same ESS 40.7:
  +0.079, 0.82x MDE, fold CI [+0.012, +0.133]); the repulsion term alone is null (+0.007,
  0.14x); adding the repulsion to the ranker weights does not rescue them (b = 1, g = 1:
  +0.081, 0.84x). At b = 2 the cell is WORSE beyond its MDE (+0.161, 1.16x) and worse than its
  control (+0.112, fold CI [+0.028, +0.193]). Built chain: (1,0) +0.074 (0.70x), (1,1) +0.071
  (0.66x), (0,1) +0.015 (0.29x); all harmful-side, none measured; the DIVW chain grid is NOT A SIGNAL.

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

---

# C2 -- THE RECOGNITION AUDIT (second-wave brief `s27/briefs/S28C2.md`; ORACLE DIAGNOSTIC THROUGHOUT)

Pre-registration `s27/PREREG_S28_C2.md` (base + addendum 1, before the number each governs).
Code `s27/s28_C2_recog_audit.py`; tests `tests/test_s28_C2.py` (6 pass). Ledger S28-L35 (CA level),
S28-L37 (the pool-member veto accepted), S28-L48 (built chain, the verdict). Every structure except production and the controls'
DIRECTIONS is chosen against the native; every number from them is ORACLE; nothing is deployable.

## C2.0 One paragraph
Which native-free scorer in the S27 library prefers the 0.29 A ORACLE structures the amplitude
family expresses (S28-L1b) to the 3.05 A production average? At the CA level, 9 of 15 prefer the
average with the fold CI below 0.5 (DIS on 79%, the consistency scorers on 80 to 98%, DISTPOT,
RG_LAW, SS_MATCH), EXVOL is uninformative (93 ties), ENV / HP / RG_UNIV are coin tosses, CONTACT
clears the registered bar by 0.004 on the best-of-5 only, and ONE scorer, CAGEO (the universe-fitted
CA virtual-angle / virtual-torsion potential), prefers the ORACLE structure on 61% [0.55, 0.69], beats
the random-signed control by +0.19 [+0.14, +0.25] (second seed +0.17), does so on the single start
as well, holds on FAIL18 (11/18) and the 108 (0.61), and sits exactly at the 95th percentile of the
max-over-15 sign-flip null (p 0.050). Head-to-head against a random signed combination at the same
displacement (production absent) DIS prefers the ORACLE structure on 89% and CAGEO on 96%: the
scorers are not blind to structure among signed combinations; they rate the contracted average
above both. The learned linear combination reaches 0.97 held-out but is anti-production (it prefers
random signed combinations on 0.89 and cannot tell ORACLE from random head-to-head, 0.53).
On the BUILT CHAIN (the reporting basis, S28-L48; 31 scorers, one max-over-31 null, the pool-member
control on 27): the closure claim stands. 20 of 31 scorers prefer the projected production average to a
0.25 A ORACLE chain with the fold CI below 0.5 (18 informative, 2 more tie-dominated; DIS on 93%, LEG on 79%),
3 more are tie-dominated, 7 are coin tosses; CAGEO's CA-level preference collapses to anti-recognition (0.611 -> 0.421) once the projection
gives production the same 3.80 A bonds, as lane D registered (S28-L36(a)); the one scorer that clears both
registered clauses on both ORACLE starts, CONTACT@chain (0.583), is the expected maximum of 31 null
scorers (p_max 0.388) and, on the pool-member control, measured in sign (+0.140, lane D's max-over-27 null
p 0.010), Type-M in size (1.23x) and deployably worthless (S27 priced CONTACT as a selector at +0.515 A;
S28-L49/L50): the same marginal cell as at the CA level, not a candidate. No candidate objective is named
for a lane A2 run.

## C2.1 ORACLE DIAGNOSTIC: the CA-level table (S28-L35)
`s27/results/s28_C2_ca_summary.json`, `s28_C2_ca_seed2_summary.json` (jobs `s28C2_ca_v2`, 176 s,
peak RSS 0.348 GB; `s28C2_ca_seed2_v2`, 226 s, 0.343 GB). ORACLE mean RMSD of the ladder: PROD
3.048, sub0 0.620, circ_s0 0.385, circ_best 0.288, NATIVE 0.000; controls RAND_SIGNED 3.871,
GAUSS_0.3 3.064, GAUSS_MATCHED 4.170. Geometry (mean virtual bond / Rg): PROD 2.96 / 6.21, the
ORACLE ladder 3.71 to 3.81 / 6.56 to 6.60, RAND_SIGNED 3.76 / 6.74 (bond-matched to the ORACLE
structures), GAUSS_MATCHED 5.39 / 7.11. Ties count 0.5.

| scorer | pref(circ_best) [fold CI] | circ_s0 | NATIVE | RAND_SIGNED | vs RAND_SIGNED [fold CI] (seed 2) | h2h vs RAND | ladder rho | FAIL18 k/18 | reading |
|---|---:|---:|---:|---:|---|---:|---:|---:|---|
| DIS | 0.206 [0.116, 0.293] | 0.206 | 0.222 | 0.030 | +0.177 [+0.099, +0.250] (+0.167) | 0.887 | -0.18 | 1 | anti |
| DIS_MEAN | 0.151 [0.098, 0.220] | 0.151 | 0.143 | 0.048 | +0.103 [+0.056, +0.163] (+0.097) | 0.752 | -0.31 | 0 | anti |
| CONTACT_LL | 0.321 [0.213, 0.421] | 0.313 | 0.329 | 0.093 | +0.228 [+0.126, +0.309] (+0.234) | 0.717 | -0.10 | 3 | anti |
| DISTPOT | 0.349 [0.264, 0.427] | 0.349 | 0.373 | 0.548 | -0.198 [-0.279, -0.124] (-0.204) | 0.369 | -0.18 | 7 | anti |
| CONTACT | 0.579 [0.504, 0.643] | 0.548 | 0.567 | 0.438 | +0.141 [+0.075, +0.203] (+0.132) | 0.624 | +0.08 | 10 | marginal (best-of-5 only) |
| ENV | 0.532 [0.430, 0.634] | 0.524 | 0.516 | 0.458 | +0.073 [-0.033, +0.166] (+0.085) | 0.563 | +0.02 | 10 | neither |
| HP | 0.544 [0.474, 0.605] | 0.544 | 0.544 | 0.507 | +0.037 [-0.058, +0.119] (+0.045) | 0.557 | +0.01 | 11 | neither |
| RG_LAW | 0.357 [0.222, 0.474] | 0.373 | 0.357 | 0.317 | +0.040 [-0.079, +0.129] (+0.032) | 0.558 | -0.18 | 5 | anti |
| RG_UNIV | 0.468 [0.420, 0.521] | 0.460 | 0.468 | 0.349 | +0.119 [+0.080, +0.162] (+0.103) | 0.593 | -0.00 | 6 | neither |
| EXVOL | 0.464 [0.407, 0.517] | 0.468 | 0.480 | 0.290 | +0.175 [+0.111, +0.241] (+0.196) | 0.685 | +0.03 | 8 | uninformative (93 ties) |
| CAGEO | 0.611 [0.551, 0.691] | 0.587 | 0.619 | 0.421 | +0.190 [+0.136, +0.252] (+0.169) | 0.962 | +0.42 | 11 | BOTH clauses, on circ_s0 too |
| SS_MATCH | 0.393 [0.333, 0.444] | 0.381 | 0.333 | 0.493 | -0.100 [-0.172, -0.056] (-0.065) | 0.373 | -0.27 | 7 | anti |
| CONS_POOL | 0.024 [0.000, 0.057] | 0.024 | 0.024 | 0.137 | -0.113 [-0.182, -0.055] (-0.117) | 0.379 | -0.69 | 1 | anti |
| DMAP_CONS_POOL | 0.135 [0.107, 0.169] | 0.135 | 0.135 | 0.173 | -0.038 [-0.071, +0.019] (-0.050) | 0.627 | -0.38 | 3 | anti |
| POOLGO_POOL | 0.202 [0.137, 0.254] | 0.187 | 0.202 | 0.369 | -0.167 [-0.212, -0.128] (-0.144) | 0.415 | -0.46 | 3 | anti |

Multiplicity: best single CAGEO 0.611 against the max-over-15 sign-flip null (mean 0.570, p95
0.607): p_max 0.050. Linear combination (nested pairwise logistic): 0.968 held-out [0.957, 0.984]
against a sign-flip null of 0.496; the same rule prefers RAND_SIGNED on 0.887 and GAUSS_MATCHED
on 0.895; head-to-head ORACLE vs RAND_SIGNED 0.526 [0.393, 0.626]: anti-production, not recognition.

## C2.2 What the audit says about the objective (mechanism, ORACLE diagnostic)
The NATIVE column tracks the circ_best column within 0.02 on 14/15 scorers: the library's
dislike is of nativeness relative to the average, not of the ORACLE optimisation. DIS tells the
ORACLE structure from a random signed combination on 89% and from a Gaussian perturbation on 91%,
and rates the average above both: the average sits at the posterior's per-pair median by
construction (S23 L9, S28-L26b(c)). What an objective would have to know: that the Bayes-risk
minimiser of a 2x over-confident posterior (S25 L2) is a contracted structure, and that a 0.29 A
structure with native bond lengths scores worse under it than that contraction. CAGEO is the one
scorer in the library that reads local CA geometry rather than pair distances or pool typicality,
which is why it is the one that does not prefer the contraction; whether that survives the
projection (which re-imposes ideal CA geometry on every structure) is the built-chain question.

## C2.2b The pool-member control (S28-L36, lane D; reproduced in S28-L37)
CAGEO's candidate status is VETOED: a random pool member (a real protein trace with no
information about the native) beats the contracted production average under CAGEO as often
(0.618) as the ORACLE structure does (0.611); the paired contrast is -0.007 [-0.027, +0.013],
0.18x MDE (single start -0.031); head-to-head CAGEO rates the ORACLE structure better than a
random pool member on 24% of targets and the NATIVE on 31%. CAGEO's +0.19 over RAND_SIGNED is
protein-like local geometry against scrambled geometry, not nativeness. The closure claim
stands with this control: no scorer prefers the ORACLE structure to production more often than
it prefers an arbitrary pool member to production beyond its MDE (CONTACT_LL +0.137 at 1.28x and
CONTACT +0.137 at 1.18x are the closest, both Type-M, both on a control that was not
registered; CONTACT also fails the single-start clause). `s27/results/s28_C2_ca_summary.json ::
pool_member_control`. S28-L35's "closure falsified at the registered bar" is withdrawn
(`s27/RETRACTIONS_S28.md` R3).

## C2.3 Built chain: the verdict (S28-L48; the reporting basis)
`s27/results/s28_C2_chain_rows.jsonl` (126/126; jobs `s28C2_chain2` to the 09-14 pause at 80, `s28C2_chain3`
resumed at 81 and killed by the governor at 97.3% box RAM after one target, `s28C2_chain4` from the governor
queue, 82 to 126, exit 0, 1,685 s, peak RSS 0.355 GB; per-target checkpoint, no target computed twice) and
`s27/results/s28_C2_chain_summary.json` (`analyse_chain`, job `s28C2_analyse_chain`, 95 s, 0.339 GB). Eight
projections per target through `s12.instrument.project` in one call per target: PROD, ORACLE sub0 / circ_s0 /
circ_best, NATIVE(aff500), and ONE draw of each control (the CA level averaged four). 31 scorers: the 16
backbone scorers (RAMA, DSSPHB, ELEC, TORS_CONS_POOL, LEG and its 11 terms) and the 15 CA scorers re-evaluated
on the projected CA traces (`@chain`); one max-over-31 sign-flip null; the pool-member control (S28-L36) for
the 27 scorers with a pool channel in `s27/cache`.

Instrument. ORACLE mean RMSD of the projected structures: PROD 3.207 (anchor 3.2126; no target bit-identical
to `chain_rows.jsonl :: DIS`, 22 differ by more than 0.01 A, the largest 0.513 A on 2LNG: the projection's
branch-flip floor, S28-L18/L27b, which touches no contrast here because all eight structures of a target share
one `I.project` call), sub0 0.499, circ_s0 0.317, circ_best 0.252 (S28-L26b's emitted ORACLE to the third
decimal), NATIVE 0.084; RAND_SIGNED 3.705, GAUSS_0.3 3.230, GAUSS_MATCHED 3.383 (matched before projection,
not after). Geometry: every projected structure has a 3.80 A mean virtual bond; PROD's Rg rises from 6.21 to
6.48 A against 6.60 to 6.61 for the ORACLE ladder and the NATIVE, so the 22% contraction of the average
becomes a 2% Rg deficit and every geometric preference is between structures of identical bond length.

| scorer (built chain) | pref(circ_best) [fold CI] | circ_s0 | NATIVE | vs RAND_SIGNED [fold CI] | pool-member contrast (x MDE) | h2h vs pool member | FAIL18 k/18 | ties | reading |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| DIS@chain | 0.071 [0.034, 0.109] | 0.071 | 0.079 | +0.063 [+0.033, +0.089] | +0.051 (0.80x) | 0.632 | 0 | 0 | anti (from 0.206) |
| DIS_MEAN@chain | 0.071 [0.049, 0.098] | 0.071 | 0.063 | +0.032 [+0.009, +0.052] | +0.046 (0.73x) | 0.617 | 0 | 0 | anti |
| CONTACT_LL@chain | 0.234 [0.168, 0.311] | 0.238 | 0.242 | +0.087 [+0.012, +0.154] | +0.100 (1.02x) | 0.608 | 1 | 13 | anti |
| DISTPOT@chain | 0.341 [0.319, 0.361] | 0.373 | 0.349 | +0.103 [+0.051, +0.148] | +0.071 (0.63x) | 0.568 | 6 | 0 | anti |
| CONTACT@chain | 0.583 [0.530, 0.645] | 0.575 | 0.575 | +0.107 [+0.009, +0.199] (0.71x, NOT MEASURED) | +0.140 (1.23x, 4/5 folds) | 0.594 | 8 | 9 | both clauses; p_max 0.388; Type-M: marginal |
| ENV@chain | 0.492 [0.408, 0.581] | 0.484 | 0.484 | +0.091 [+0.033, +0.163] | +0.131 (1.19x) | 0.612 | 11 | 16 | coin toss |
| HP@chain | 0.456 [0.361, 0.541] | 0.421 | 0.417 | +0.024 [-0.034, +0.093] | +0.018 (0.15x) | 0.545 | 9 | 13 | coin toss |
| RG_LAW@chain | 0.405 [0.311, 0.481] | 0.413 | 0.421 | +0.071 [-0.009, +0.121] | -0.030 (0.27x) | 0.480 | 6 | 0 | anti |
| RG_UNIV@chain | 0.397 [0.336, 0.451] | 0.389 | 0.349 | +0.048 [-0.042, +0.114] | +0.001 (0.01x) | 0.485 | 6 | 0 | anti |
| EXVOL@chain | 0.417 [0.368, 0.462] | 0.409 | 0.425 | +0.063 [+0.011, +0.118] | -0.076 (1.60x) | 0.425 | 6.5 | 101 | uninformative |
| CAGEO@chain | 0.421 [0.370, 0.492] | 0.381 | 0.452 | +0.238 [+0.143, +0.359] | -0.108 (1.21x, 5/5) | 0.297 | 7 | 0 | anti (from 0.611) |
| SS_MATCH@chain | 0.393 [0.340, 0.443] | 0.381 | 0.377 | -0.155 [-0.228, -0.085] | +0.072 (0.93x) | 0.591 | 5 | 11 | anti |
| CONS_POOL@chain | 0.056 [0.037, 0.089] | 0.056 | 0.063 | -0.127 [-0.185, -0.070] | (pool-relative) | | 2 | 0 | anti |
| DMAP_CONS_POOL@chain | 0.143 [0.109, 0.171] | 0.135 | 0.143 | -0.095 [-0.157, -0.025] | (pool-relative) | | 3 | 0 | anti |
| POOLGO_POOL@chain | 0.345 [0.261, 0.414] | 0.349 | 0.345 | -0.063 [-0.145, +0.007] | (pool-relative) | | 4 | 13 | anti |
| RAMA | 0.413 [0.336, 0.496] | 0.444 | 0.476 | +0.206 [+0.135, +0.290] | -0.377 (3.50x; cross-basis) | 0.117 | 9 | 0 | anti |
| DSSPHB | 0.226 [0.198, 0.256] | 0.210 | 0.258 | -0.171 [-0.207, -0.134] | -0.104 (1.07x; cross-basis) | 0.403 | 3 | 1 | anti |
| ELEC | 0.484 [0.416, 0.560] | 0.484 | 0.492 | -0.063 [-0.119, +0.000] | -0.016 (0.15x) | 0.516 | 8.5 | 18 | coin toss |
| TORS_CONS_POOL | 0.262 [0.226, 0.306] | 0.286 | 0.278 | +0.079 [+0.008, +0.147] | (pool-relative) | | 3 | 0 | anti |
| LEG (total) | 0.206 [0.137, 0.289] | 0.206 | 0.238 | +0.079 [+0.007, +0.168] | -0.221 (2.37x; cross-basis) | 0.340 | 2 | 0 | anti |
| LEG_hbond_local | 0.190 [0.107, 0.269] | 0.218 | 0.250 | -0.016 [-0.073, +0.015] | -0.226 (2.60x; cross-basis) | 0.317 | 3 | 8 | anti |
| LEG_solvation | 0.349 [0.278, 0.436] | 0.365 | 0.365 | +0.008 [-0.094, +0.136] | -0.013 (0.12x) | 0.545 | 9 | 0 | anti |
| LEG_torsion | 0.405 [0.320, 0.496] | 0.333 | 0.452 | +0.254 [+0.180, +0.322] | -0.293 (3.03x; cross-basis) | 0.221 | 8 | 0 | anti |
| LEG_coop_helix | 0.361 [0.305, 0.427] | 0.349 | 0.369 | +0.071 [+0.054, +0.088] | -0.146 (2.50x; cross-basis) | 0.366 | 5 | 79 | anti, tie-dominated |
| LEG_steric | 0.508 [0.384, 0.632] | 0.476 | 0.520 | +0.266 [+0.150, +0.397] | -0.116 (1.41x; cross-basis) | 0.399 | 6 | 54 | coin toss |
| LEG_contact | 0.532 [0.459, 0.615] | 0.540 | 0.516 | -0.012 [-0.068, +0.084] | +0.021 (0.18x) | 0.521 | 10 | 2 | coin toss |
| LEG_electrostatic | 0.492 [0.423, 0.562] | 0.484 | 0.484 | -0.048 [-0.089, +0.009] | -0.008 (0.07x) | 0.515 | 8.5 | 18 | coin toss |
| LEG_compactness | 0.468 [0.393, 0.533] | 0.464 | 0.480 | +0.103 [+0.032, +0.155] | +0.035 (0.31x) | 0.528 | 9 | 14 | coin toss |
| LEG_hbond_longrange | 0.548 [0.533, 0.562] | 0.556 | 0.583 | -0.004 [-0.033, +0.023] | -0.006 (0.17x) | 0.493 | 10.5 | 114 | uninformative (12 decided targets) |
| LEG_coop_sheet | 0.508 [0.500, 0.517] | 0.504 | 0.520 | +0.008 [+0.000, +0.017] | -0.003 (0.17x) | 0.497 | 9 | 124 | uninformative |
| LEG_aromatic | 0.429 [0.357, 0.500] | 0.452 | 0.468 | +0.008 [-0.052, +0.080] | -0.028 (0.35x) | 0.506 | 10.5 | 64 | uninformative |

Multiplicity: best single CONTACT@chain 0.583 against the max-over-31 sign-flip null (mean 0.579, p95
0.627): p_max 0.388 (lane D independently, 4,000 draws: 0.577 / 0.619 / 0.399, `s28_D_c2_chain_null.json`).
Linear combination (31 scorers, nested): 0.960 [0.938, 0.983] held out against a null of 0.502; the same
rule prefers RAND_SIGNED on 0.952 (contrast +0.008 [+0.000, +0.022], one target) and GAUSS_MATCHED on
0.960 (+0.000); head-to-head ORACLE vs RAND_SIGNED 0.492 [0.364, 0.611]: anti-production, purer than at
the CA level (0.887 / 0.526).

The pool-member control on the backbone scorers is CROSS-BASIS (the members' cached values are on their
real torsions, the audited structures on projected torsions) and the projection's torsions are atypical
of real fragments: the projected NATIVE itself is worse than 93% of its pool under RAMA (median), 78%
under LEG_torsion, 65% under LEG. The large negative contrasts on RAMA, LEG_torsion, LEG_hbond_local,
LEG_coop_helix, LEG, DSSPHB, LEG_steric measure that signature, not anti-recognition; on those channels
the same-code-path contrasts (vs RAND_SIGNED / GAUSS) and the ladder are the like-for-like reading. For
the CA scorers on the chain the control is like for like (3.8 A bonds on both sides), which is exactly
why it now resolves CAGEO.

Verdict. The closure claim ("no scorer in the S27 library recognises the ORACLE structures the amplitude
family expresses") STANDS on the built chain over all 31 scorers with the registered multiplicity pricing
and the pool-member control. The literal two-clause falsifier is met by one scorer of 31, CONTACT@chain,
on both ORACLE starts; its pref is the expected maximum of 31 null scorers (p_max 0.388), its clause-2
margin is +0.107 at 0.71x MDE on a single control draw (`ST.compare`: NOT MEASURED), its pool-member
contrast is +0.140 at 1.23x with 4/5 folds (Type-M), its head-to-head against a random pool member is
0.594 (the NATIVE 0.590), and it is below half on FAIL18 (8 of 18): the marginal, Type-M cell S28-L36 named
and S28-L37 accepted, unchanged by the projection, not a candidate. ENV@chain (+0.131, 1.19x) and
CONTACT_LL@chain (+0.100, 1.02x) are in the same zone on the pool-member contrast and fail clause 1. Lane D's
S28-L49 sharpened the reason (adopted in S28-L50): the contact family's pool-member contrast is NOT a chance
maximum (per-target sign-flip null over the 27 scorers, p_max 0.010), so its sign is measured and its size is
Type-M; what it measures is that a native-like structure sits 4 percentile points better than the average
in a pool where the average is a median member under CONTACT (pct 0.443, against DIS's 0.020), a separation
of 0.59 that is below DIS's own 0.63; and S27 already consumed that signal on the same 126 targets
(`s27/results/pool_rows.jsonl`: CONTACT as the selector +0.515 A, ENV +0.419, CONTACT_LL +0.521, every
DIS+CONTACT admixture null-to-worse). "Measured in sign, Type-M in size, deployably worthless" is the
reading; multiplicity of the registered statistic (p_max 0.388) is the decision rule. The
registered prior holds where it was stated; the one CA-level exception, CAGEO, was the contraction of the
average and collapses to anti-recognition on the projected chains (0.611 -> 0.421; pool-member contrast
-0.007 -> -0.108; pct(PROD) median 1.00 -> 0.60), as S28-L36(a) registered. What an objective would have
to know is unchanged and now holds with ideal geometry on both sides: the projected Bayes-risk minimiser of
the over-confident posterior is a better structure than a 0.25 A chain to the objective itself (93% of
targets; the projected average is better than 99% of its own pool under DIS) and to every consistency,
geometric and legacy-energy channel in the library. C2 is closed on the reporting basis with no candidate
objective named for a lane A2 run.

## C2.4 What damaged my own expectations
- I expected CAGEO's CA-level preference to be the one live cell going into the chain and wrote S28-L35's
  verdict as "falsified at the registered bar" before the pool-member control existed. The control that
  killed it was not one of the three I registered: a random signed combination has protein-like bonds but
  scrambled virtual angles, Gaussian noise has neither, and neither could see that a 22%-contracted
  average is off the virtual-angle distribution CAGEO was fitted on. Lane D's control (a real trace with no
  information about the native) saw it in one line, and the projection confirmed it on the reporting basis:
  CAGEO went from 0.611 to 0.421 the moment production had real bonds. "Beats a random signed control" is
  not "recognises the ORACLE structure"; the control has to have the property the scorer measures.
- I expected the backbone scorers to add information the CA scorers lacked. They add none in the
  recognition direction: none of the 16 clears clause 1, the legacy energy prefers the production chain to
  the 0.08 A NATIVE chain on 76% of targets, and five are tie-dominated on this basis.
- I expected the pool-member control to be a single fair yardstick on the chain. For the backbone
  scorers it is cross-basis, and the projected NATIVE at the 93rd RAMA percentile of its own pool says
  the projection's torsion signature is larger than any recognition signal on those channels. Two sides of
  a chain contrast must share a code path (S28-L27b); a control drawn from the cache does not.
- I expected the shipped objective's anti-preference to weaken once the average was de-contracted. It
  strengthened (DIS 0.206 -> 0.071; pct(PROD) 0.05 -> 0.01): the projected average is closer to the
  posterior's mode than the contracted one, and a 0.25 A structure with the same bonds is further from it.

## C2.5 What I did not do and why
- AMB on the projected chains (brief, prereg section 3): deferred throughout; one AMBER process per box
  under the user's own 84 to 98% load was not guaranteed to fit (contract rule 8), and the 30 non-AMBER
  scorers gave the closure verdict with margin; recorded here as not run.
- A second control seed on the chain (4 more draws x 3 controls x 126 projections, about 1.3 h): the
  positive-only replication rule (prereg section 5); no positive fired. The chain's single control draw is
  stated beside every clause-2 number.
- A null for the multiplicity of the pool-member contrast across the 27 scorers: not registered and not run
  by me; lane D ran it (S28-L49: p_max 0.010) and S28-L50 adopts it. The registered pricing (the best single
  pref as an order statistic) is inside its null and decides.
- The seed-2 pool-member control at the CA level WAS regenerated on 09-19 (`s28_C2_ca_seed2_summary.json`);
  every circ_best / circ_s0 contrast is identical to seed 0 to the fourth decimal, as it must be (the
  control draws do not enter that contrast); only the RAND_SIGNED-vs-pool-member head-to-head moves, by at
  most 0.025.
- Re-optimising the circuit under any scorer (lane A2's job if one had passed): none passed.
