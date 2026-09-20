# LANE B -- THE COMPATIBILITY HAMILTONIAN, AND THE TAIL-THEN-AGGREGATE OBJECTIVE. SPRINT 29. FINDINGS.

Pre-registration `s29/PREREG_S29_B.md` (base `1c345f07`; addendum 1 `9d745692`; addendum 2
`d3ccf6b6`; addendum 3 `c89e3106`) -- every one committed before the number it governs. Code
`s29/s29_B_compat.py` (the three matrices, the PERM and SPEC controls, an independent hop-gradient
estimator, the exact ground state, the readouts), `s29/s29_B_tta.py` (the tail-then-aggregate
objective with its envelope gradient, the flatness report, the exhaustive subset search, the
endpoint and the built-chain runner), `s29/s29_B_analyse.py`. Tests `tests/test_s29_B.py`
(19 pass). Results `s29/results/s29_B_*.json|jsonl`. Ledger: **S29-L25** (the set-equality
counterexample, 12 targets), **S29-L27** (the flatness measurement), **S29-L36** (measurement 1),
**S29-L45** (the subset search at 126 with the tie fix), **S29-L49** (measurement 2, the gate),
**S29-L54** (the endpoint on the built chain). Lane D's checks: **S29-L38**. Tiers as in S12 to
S28: DEMONSTRATED / ORACLE DIAGNOSTIC / DERIVED / REFUTED / OPEN.

BASIS NOTICE. Two bases, never in one column: the BUILT CHAIN (`rmsd_chain`, production in this
lane's process 3.2105 -- see the anchor note below) and the POINT CLOUD (production 3.048338 at
n = 126, 3.252928 on the 12 trainability targets). The verdict is on the built chain.

ANCHOR NOTE (carried into every chain number). This lane's production point cloud equals S27's
`chain_rows.jsonl :: DIS` `rmsd_cloud` to **3.7e-14** on 126/126. Its PROJECTION differs by up to
**0.5174 A** on 125/126 targets (means 3.210534 vs 3.212625): the S28-L18 / L27b / L43
input-difference floor of `s12.instrument.project`, reproduced to the digit. It does not enter any
contrast here, because every arm and its comparator are projected in one process from clouds built
by one loader.

ORACLE NOTICE. Every RMSD below is an ORACLE evaluation of an ACHIEVABLE (native-free) selection.
The matrices, ground states, subsets, f, the tie draws, the frame and the readouts never see a
native (`tests/test_s29_B.py :: test_selection_is_native_free_and_nan_poison_is_bit_identical`).

---

## 1. DERIVED -- THE SIGN-MIXING LEMMA (S29-L36, S29-L49)
A_c 1 = H A H 1 = 0 and G 1 = Delta(Delta^T 1)/n_res = 0, so every eigenvector with lambda != 0 is
orthogonal to 1 and **has entries of both signs**. Measured as an exact zero: the top eigenvector's
squared overlap with the uniform state is **0.0000** for A_c and G at every register size, against
0.96 to 0.99 for S28's A. A centered off-diagonal term's content is a signed contrast; every
deployed readout is a function of p and discards the sign. Reached independently of lane T's
S29-L11(d) (my prereg section 3.2 predates reading it), asserted in code.

## 2. DERIVED -- THE STABLE-RANK BOUND (S29-L36)
"Within 30x of the diagonal terms" at D = 512 needs r_stable >= 512^2 x 3.051e-2/30 = **266**.
rank(G) <= 3 n_res - 3 = 45 exactly on posed windows, and **rank(G) saturates at 30 measured while
the register grows 128 -> 256 -> 500**. No deviation Gram can be gradient-visible here, by
construction. Measured r_stable at n = 9: 1.036 (A) / 1.591 (A_c) / 1.675 (G).

## 3. DEMONSTRATED -- MEASUREMENT 1, AND A FOURTH CHECK ON LANE T's LAW (S29-L36)
An independent implementation reproduces lane T's 216 gradient cells to **2.03e-14** (median
1.27e-15) and S28's own `hop_only|J1` row to **4e-16** at n = 4..8 (7.0e-2 at n = 9, the same
register-ordering difference T reported); parameter shift against central differences 1.46e-09.
**The bright line B1 is REFUTED on both clauses**: slopes **-2.305** (A_c) and **-1.900** (G),
steeper than A's -1.830 rather than shallower than -1.0, and n = 9 variances 3.60e-6 / 5.05e-6,
283x and 201x below the 30x threshold. Centering removes the lambda_2/lambda_1 degeneracy
(0.1385 -> 0.4654 -> 0.6342) and makes the decay WORSE. Lane T's T1 and T2 are not falsified.
**The justification that spawned this lane is dead**, and for a structural reason.

## 4. ORACLE DIAGNOSTIC -- MEASUREMENT 2: THE GATE IS CLOSED (S29-L49)
Of 24 REAL (M, J) cells, **zero** clear 0.7x MDE against the same target's own DIS top-75 average
on either readout, so conditions 2 to 5 were never reached and **measurement 3 was not run**. The
signed readout is 0.3 to 4.7 A worse with a **coin-flip sign** (sign_correct 0.40 to 0.75, median
0.50), and lane O's S29-L21 has that family's ORACLE ceiling at exactly 0.0000 A. The pole
cancellation is visible: `pole_imbalance` collapses from -1.00 to -0.26 as the state spreads over
M's principal contrast. **Lateral motion**: the p-readout's emitted cloud moves **0.22 to 1.45 A**
from production while its RMSD moves under 0.1 A -- which falsifies S29-L11 prediction 3 as I
operationalised it (a 0.05 A floor) and confirms the law behind it. At n = 12 with SE 0.5641
nothing here is a measured contrast; the gate is a decision, and it decided no.

## 5. DEMONSTRATED -- THE SET-EQUALITY THEOREM FAILS ON REAL POOLS (S29-L25, S29-L45)
Exhaustive over all C(500, 2) = 124,750 pairs, **126 targets**: the f-optimal pair is **not the
energy prefix on 114/126**, the f-optimal m = 5 subset **not on 124/126**, and the per-state sort
finds the exhaustive optimum on only **12/126**. Mean objective gaps +0.1041 (pair) and +0.1499
(m = 5), clearing lane T's registered 0.10. The optimum is reachable neither by sorting E nor by
sorting f: the first measurement in this project's record of "which set" as a real optimisation
variable. Independently reproduced by lane D from the entry's description and under the exact
shipped risk lookup, with zero sign flips (S29-L38).

## 6. DEMONSTRATED -- AND THE ESCAPE BUYS NO ACCURACY, WITH THE DECOMPOSITION AS THE FINDING (S29-L45)
The f-optimal m = 5 subset is **+0.2451 A worse than production** (1.45x MDE, fold CI
[+0.171, +0.326], 5/5 folds, power 0.98, Type-M 1.01, not concentrated) -- but **+0.1647 of that is
simply m = 5 versus m = 75**, and the part attributable to choosing the f-optimal set rather than
the prefix at the same m is **+0.0804 at 0.73x MDE, NOT MEASURED** (fold CI [+0.028, +0.135], 4/5).
Read before the sign, as lane D's field survey requires: MSET_5 has signed cosine +0.063 and
per-target |cos| 0.271, so a harmful result is the expected value of an unsigned displacement.

## 7. DEMONSTRATED -- THE FLATNESS MEASUREMENT, AND ITS GATE CORRECTED IN ADVANCE (S29-L27)
Registered before measuring: the f term of the tail-then-aggregate objective is flat on **exactly**
the CVaR term's subspace, because dR_alpha/dp_y = 0 above the VaR by the same envelope argument, so
the raw 85.5% cannot fall and the quantity a gate must read is the OVERLAP. Measured on 12 targets:
flat_cvar = flat_f = flat_info_bearing = **0.8552** identically at every lam; the TTA readout's own
flat set is the same; **overlap = 0.0000** -- all 74 directions the readout consumes are seen by
the objective, against the deployed readout, which moves along **0 of 511** continuous directions.
Readout anchor: R_alpha is +0.0271 from production (0.43x, NOT MEASURED) with radius of gyration
0.9995 and mean bond 1.0036 of production's, so the sharper average does not contract.

## 8. REFUTED -- THE TAIL-THEN-AGGREGATE ENDPOINT (S29-L54)
**F5b refuted as registered; the registered prior (mine and the coordinator's) held on both
clauses.** On the built chain, 126/126, 9 arms: no lam beats production or its own lam = 0 on
either basis or either seed; lam = 1 is +0.1042 / +0.1157 above production (0.95x / 1.03x, one cell
WORSE) and lam = 3 is +0.1055 / +0.1028 (0.80x / 0.81x). The mechanism is alive and irrelevant:
as lam rises the objective takes the tail (m 74 -> 39, PR 407 -> 164, Jaccard with the DIS top-75
**0.933 -> 0.518**) while the **deployed tail-SET readout stays flat at every lam** (largest
|effect| 0.046, every cell NOT MEASURED). Half the set is replaced and the emitted accuracy does
not follow -- the terminal-operator law (S29-L44 addendum 2) plus the pool's 68% common-mode error.
M6, the target-independent fixed profile (my implementation reproduces lane T's exact optimum to
F -4.7236 vs -4.7237), is **0.31x MDE from the trained circuit**: the quantum stage does not
measurably beat a profile with no circuit in it. The lam grid does not transfer (split-half
+0.047 / +0.053 of a -0.168 oracle).

## 9. A CLAIM NOT MADE, AND A DEFECT FIXED
- **Not made**: my `tail_is_prefix` column reads 0.91 to 0.99 and would have said set-equality
  broke inside the VQE. It is an artefact -- the value-based gate passes with `equality` on every
  cell tested, and the 37 flagged targets are exactly those with exact E ties (up to 4-way). The
  CVaR tail under TTA remains the energy prefix **by construction**, because the tail's order stays
  a per-state scalar, which is the condition lane T attached to the derivation.
- **Fixed**: the exhaustive pair search took `tie[0]` -- array order -- under a comment claiming the
  opposite, and on a DIS-sorted pool that is biased toward the prefix the experiment tests (the
  project's named failure mode, `tie-breaking-leaks-the-pool-order`). Corrected in `7b2e83e1` to one
  global argmin with a stable RNG tie draw, the same rule in the greedy and swap steps, the tie-set
  size on every row (max 4 at n = 126, so the fix was not cosmetic), and a regression test. Lane D
  verified the defect moved no number in S29-L25; the pre-fix rows are kept, not deleted.

## 10. WHAT THIS LANE LEAVES
The off-diagonal route on the **candidate-index register** is shut, and not for S28's reason. S28
closed it on a near-rank-one spectrum; that reason is retired (section 3). It is shut because every
common-mode-removed operator's informative eigenvectors are sign-mixed, a p-readout cancels them,
a signed readout consumes them and is worse with a coin-flip sign and an ORACLE ceiling of zero,
and what motion survives is lateral. That is a statement about the ENCODING, not about
off-diagonality, and it points at lane T's S29-L15 Q2: a local mixer on a configuration space.
The constructive half is the pair S29-L27 + S29-L54: the flatness was fixable, fixing it changed
nothing, and the barrier is the information content of f.

## 11. NOT DONE (named rather than omitted)
- **Measurement 3** (the compatibility-Hamiltonian CVaR-VQE endpoint): not run, by the gate. This
  is a decision the prereg made, not an omission.
- **G_res** (deviations residualised on the posterior's own predicted deviation): registered as NOT
  RUN in prereg section 6, with the reason (it needs a second construction with its own convenience
  choices, and by rule 21 the marginal class it lives in is bounded).
- **The `cost_tail_avg` meter submission**: prereg B2.5 registered f for lane D's meter before the
  endpoint. The endpoint ran on the coordinator's release after lane D's band measurement
  (S29-L33) came back null, and f is the SHIPPED score, whose four meter numbers lane D already
  published in S29-L2 (ladder rho -0.402 chain / -0.182 CA, cosine -0.034, native percentile 0.369,
  preference 0.071 / 0.206 against a pool-member control at 0.020 / 0.126). No separate submission
  was made, and no new cost function was introduced that would have needed one.
- **The tailset readout on the built chain**: only `C_Ralpha` was projected (9 arms x 126). The
  tail-SET readout's chain numbers are not measured; its point-cloud numbers are flat at every lam,
  and the arms that would have been projected under the two-sided 0.7x rule are the R_alpha ones.
- **A second seed for the subset search**: the exhaustive pair search is deterministic given the
  tie draws; only the tie draws carry a seed, and they were not varied.
