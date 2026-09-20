# LANE B -- THE FIRST NON-DIAGONAL HAMILTONIAN. SPRINT 28. FINDINGS.

Pre-registration `s27/PREREG_S28_B.md` (base + addendum 0 before any endpoint number;
addenda 1 and 2 for B2). Code `s27/s28_B_hop.py` (the module), `s28_B_train.py` (F5,
checkpointed), `s28_B_rank1.py` (mechanism), `s28_B_split.py` (three-way split), `s28_B_mladder.py`
(T9 decomposition), `s28_B_analyse.py` (statistics), `s28_B_chain_part.py` (optional half-split
chain runner, unused), `s28_B_represent.py` (the representability fit, S28-L26), `s28_B_prodcheck.py`
(production re-projected in this lane's process), `s28_B2_knn.py`, `s28_B2_rank1.py` (B2). Tests
`tests/test_s28_B.py` (16 pass), `tests/test_s28_B2.py` (8 pass); lane D's `tests/test_s28_D.py` adds the Perron reading.
Results `s27/results/s28_B_*.json|jsonl`, `s28_B2_*.json|jsonl`. Ledger: S28-L8b (F5), S28-L21
(point cloud, intermediate), S28-L25 (B2 trainability), S28-L29 (B2 share and decomposition),
S28-L41 (THE VERDICT on the built chain, with the representability fit), S28-L44 (answers
S28-L43), S28-L46 (the B2 endpoint, closed); lane D's checks S28-L2, S28-L9, S28-L11, S28-L22,
S28-L23, S28-L26, S28-L31, S28-L38, S28-L43. Tiers as in S12 to S25:
DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN. Every number carries its
artefact path.

BASIS NOTICE. Two bases appear and never share a column: the BUILT CHAIN (`rmsd_chain` in
`s28_B_chain_rows.jsonl`, the production projection `s12.instrument.project` of a coordinate
average; production on this basis 3.2126, S27 `chain_rows.jsonl :: DIS`) and the POINT CLOUD
(`R1/R2/R3.rmsd` in `s28_B_rows.jsonl`; production 3.048338). The verdict is on the built chain
(contract addendum 1, item 13); the point cloud screened which arms were projected.

ORACLE NOTICE. Every RMSD is an ORACLE evaluation of an ACHIEVABLE (native-free) selection.
The selections, weights, graphs, J, sigma and readouts never see a native
(`tests/test_s28_B.py :: test_select_target_is_native_free_and_nan_poison_is_bit_identical`).

---

## 0. WHAT WAS BUILT

H = diag(E) - J A on the 9-qubit candidate register: E = zrank(shipped distogram score) over
the 500 windows plus 12 padding states at max + 10 sd; A the Gaussian kernel exp(-d^2 / 2
sigma^2) on the pairwise CA-RMSD of the windows, sigma = median pairwise RMSD (3.97 to 5.00 A
on the probe), zero diagonal, zero padding rows, scaled to unit spectral norm so the hopping
term is bounded by J in the units of E (sd 1, rank gap 0.0069, range 3.464). The objective
keeps the spine: F(theta) = CVaR_0.18(E; p_theta) - 0.5 H(p_theta) - J <psi_theta|A|psi_theta>;
the first two terms and their gradient are `core.quantum.free_energy` untouched, the hopping
term's gradient is the two-term parameter-shift rule on the batched shifted STATES
(`grad_hop_paramshift`; equals central finite differences to 1e-6, cosine > 1 - 1e-9). The
training loop is `core.quantum.run_cvar_vqe` line for line; at J = 0 every float is identical
to `s24.d_harness.arm_vqe` (126/126 targets, both seeds, `s28_B_summary.json ::
anchors.s27_bit_identical_s*`; lane D's independent check S28-L9).

For J > 0 the objective is not a function of p: A >= 0 entrywise, so two states with the same
p differ in F through their amplitude signs (`test_objective_is_not_a_function_of_p_for_J_
positive`). This is the first objective in the project with that property.

Three readouts of a probability vector (VQE optimum, or the exact ground state of H by
`numpy.linalg.eigh`): R1 the alpha-tail along the E order, uniform average (the deployed
readout); R2 the p-weighted average over all 500 in the p-weighted consensus-medoid frame
(S23 L8's convention on the full p); R3 the 75 most probable candidates, uniform average, ties
by the stable key. Controls: J = 0; PERM (P A P^T, same spectrum and degree multiset,
correspondence with E destroyed); RAND (a random geometric graph Sinkhorn-scaled to the same
degree SEQUENCE, similarity structure destroyed); the random-75 null (S27 `pool_rows.jsonl`,
3.4209). Grid J in {0, 0.1, 0.3, 1, 3}, seeds 0 and 1, 39 arms per target, 126 targets
(`s28_B_rows.jsonl`, 4,914 rows; jobs `s28B_run` / `run2` / `run3`, peak RSS 0.35 GB, 27 to 45
s per target).

---

## 1. THE QUANTUM-SIDE MEASUREMENTS (no RMSD)

### 1.1 The hopping term's gradient variance decays at -1.7 to -1.8 per qubit; the full objective's is flat. DEMONSTRATED (property).

`s27/results/s28_B_train.json :: summary` (S28-L8b; lane D S28-L11 STANDS WITH CAVEAT, the
caveats accepted in S28-L21). Var[dF/dtheta_0], theta ~ N(0, 0.6^2), depth 3, exact parameter
shift, 120 draws, median over S27's 12 trainability targets, E the rank ladder over the DIS
top-2^n, A their real graph:

    cell                 n=4        n=5        n=6        n=7        n=8        n=9    log2 slope/qubit
    full J=0        2.600e-02  3.174e-02  2.129e-02  2.092e-02  1.725e-02  3.051e-02   -0.043
    full J=1        3.016e-02  3.239e-02  2.167e-02  2.089e-02  1.753e-02  3.055e-02   -0.075
    full J=3        6.265e-02  3.815e-02  2.418e-02  2.118e-02  1.833e-02  3.065e-02   -0.243
    hop-only (J=1)  4.061e-03  7.284e-04  2.817e-04  5.898e-05  3.948e-05  4.157e-06   -1.844 (4..9), -1.700 (4..8)
    linear diag     5.811e-02  3.124e-02  1.937e-02  2.436e-02  4.441e-03  4.495e-02   -0.285

At the deployed width the diagonal terms' gradient variance is 7,300x the hopping term's at
J = 1 (820x at J = 3), a 1.2% ratio in gradient magnitude. Falsifier F5 (slope difference
> 0.3 per qubit between J = 3 and J = 0) did not fire (0.200). No slope is a plateau or its
absence; the circuit has 27 parameters and is nowhere near a 2-design.

### 1.2 The mechanism is the graph's near-rank-one spectrum. DEMONSTRATED (property).

`s27/results/s28_B_rank1.json :: summary`. The rank-one part v1 v1^T of A reproduces 97% of
the hop-only variance at n = 9 and its slope (-1.80 vs -1.84); the residual is 30x smaller; a
DIAGONAL observable with A's spectrum decays at the same rate on the padding-free registers
(-1.86 vs -1.70, S28-L11(3)). lambda_2 / lambda_1 = 0.11 to 0.14 at every n; the Perron
vector's participation ratio is 0.88 to 0.97 of dim and its overlap with the uniform state
0.95 to 0.99. The hopping term is, to 5%, -(<uniform|psi>)^2. The pool's similarity graph at
the brief's sigma is a typicality projector (degree anti-correlated with E at -0.47 to -0.90,
`s28_B_rows.jsonl :: deg_E_corr`).

### 1.3 On a spread-spectrum (kNN) graph the decay rate halves and the circuit still does not collect the hopping at J = 1. DEMONSTRATED (property); F5-B2 NOT PASSED as registered.

`s28_B2_train.json`, `s28_B2_share.json`, `s28_B2_rank1.json` (S28-L25, S28-L29). Symmetric
kNN graph, degree-normalised, unit spectral norm, lambda_2 / lambda_1 = 0.98 (k 10) / 0.99
(k 5) at n = 9. Hop-only slope -1.033 (k 10) / -0.997 (k 5) per qubit against -1.844; the
registered line (-1.0) is straddled, not cleared. At n = 9 the hopping variance is 29x / 60x
the Gaussian graph's, still 250x / 123x below the diagonal terms at J = 1. The remainder
A - v1 v1^T carries 96.5% / 98.6% of the kNN variance at n = 9 and decays at -0.66 / -0.85 per
qubit; the rank-one part decays at -1.78 as on the Gaussian graph. Share of the same-sign
bound collected by the trained circuit at J = 1: 3% / 12% (k 10, seeds 0 / 1), 23% / 7% (k 5),
Gaussian 8% / 8% on the same targets: the second clause ("more than half") FAILED. At J = 3,
k = 10 the circuit collects 68% / 44% (Gaussian 33% / 25%). Reading: the hopping is collected
only where its gradient variance is within about 30x of the diagonal terms' (k 10, J 3: 28x).
The k = 5 graph is disconnected on some targets (up to 4 components at n = 9); its GS rows are
degenerate and its endpoint replication is dropped (prereg addendum 2).

### 1.4 The departure from the classical prefix, at n = 126. DEMONSTRATED (property).

`s28_B_rows.jsonl`; S28-L21's table; `s28_B_mladder.json`. The alpha-tail support stays a
subset of an E-prefix at every J (4,914/4,914 rows): the tail-reading operator's property,
not a finding; and RMSD_R1 - RMSD_top-m(E) <= 1.1e-13 A on every arm, so every R1 movement at
J > 0 is the m-ladder term (+0.005 to +0.014 A at m 71 to 79). The VQE state's TV distance
from its J = 0 state: 0.025 / 0.050 / 0.098 / 0.155 at J = 0.1 / 0.3 / 1 / 3 (seed 0); realised
rung 74 / 74 / 74 / 79; its p-top-75's Jaccard with the DIS top-75 0.177 / 0.180 / 0.181 /
0.163 (random pairs 0.082); hopping value 0 / 0 / 0.007 / 0.318 against a same-sign bound of
0.89 to 0.91 (sign coherence 0.01 at J 1, 0.34 / 0.25 at J 3 by seed; the J = 3 means are
MIXTURES of fully coherent and fully incoherent cells, section 2.3). The exact ground state
is localised (PR 1.0 to 1.3, m = 1) at J <= 0.3, PR 62 at J = 1 with its p-top-75 at Jaccard
0.960 with the DIS top-75 (it RE-SELECTS the distogram's own set), PR 305 at J = 3.

### 1.5 The three-way split's middle leg (S28-L2(b)). DEMONSTRATED (property).

`s28_B_split.json` (126 targets x 13 cells). On the circuit's own objective F: at J <= 0.3 the
VQE state beats the exact ground state on 126/126 (the ground state is one-hot; -4.56 vs
-1.73 to -2.07); at J = 1 they tie (-4.566 / -4.585 vs -4.578; VQE lower on 47%); at J = 3 the
ground state has F = -7.31 against the circuit's -5.51 / -5.27 (VQE lower on 29%): from the
best of 16 untrained draws (-4.06) the circuit closes 45% / 37% of the gap and collects 0.32 of
the ground state's 0.91 hopping. The VQE beats its untrained best-of-16 on 126/126 at every J.

---

## 2. THE ENDPOINT

### 2.1 Point cloud (screening basis; S28-L21, lane D S28-L22 STANDS as an intermediate negative). No arm beats its J = 0 counterpart beyond MDE on both seeds; no readout is better than production.

`s28_B_summary.json :: contrasts`. R1 (deployed) at every J, graph, seed: within 0.013 A of
J = 0 (largest 0.74x MDE, 94 identical sets, not replicated) and within 0.014 A of production
(inside MDE). R2 and R3 are +0.24 to +0.33 A above production at every J with the fold CI
above zero (WORSE); the one cell beyond MDE against its own J = 0 (seed-1 R3 J = 0.1, -0.021,
1.05x, Type-M zone) has a seed-0 twin at +0.0001 (0.01x) and a PERM twin at -0.022 (0.82x):
dead. The best J over the four rungs is NOT A SIGNAL on every REAL and RAND cell
(`best_of_k_within`, split-half -14% to +15%, k_eff 3.2 to 4.0). REAL vs PERM and REAL vs
RAND: 72 contrasts, none beyond MDE (F3, F4 silent). F2's prior (eigensolver wins or ties)
held: at J >= 1 the eigensolver's R2/R3 beat the circuit's by 0.20 to 0.33 A beyond MDE. The
eigensolver at J = 1 re-selects the DIS top-75 (R3 = 3.0597, +0.011 vs production, 0.45x).

### 2.2 Built chain (THE VERDICT; S28-L41). The deployed readout emits production at every J; the two readouts that could consume the off-diagonal structure are 0.25 to 0.37 A worse than production at every J, as at J = 0; the eigensolver ties production. REFUTED as registered; the registered prior (WORSE or null) held.

`s28_B_chain_rows.jsonl` (2,268 rows = 126 x 18; `s28B_chain` 103 targets on 2026-09-14,
`s28B_chain_r2` 1 and `s28B_chain_r3` 22 on 2026-09-19, resumed per (arm, target) after two
kills, nothing recomputed), `s28_B_summary.json :: chain`; production = S27 `chain_rows.jsonl
:: DIS` 3.2126, RE-PROJECTED in this lane's own process on 126/126 and identical to the last
digit on every target (`s28_B_prodcheck.json`: 126/126, max |d| 0.0e+00), so every
contrast is on one code path and the S28-L18/L27b branch-flip floor does not enter. Effect =
arm - comparator, negative is better; MDE = 2.8016 x SE; the decision is on the fold CI.

    arm                      chain    | vs production        x/MDE   fold CI            | vs J = 0 R1 (same seed)  | 108 / FAIL18 vs production
    vqe|s0|NONE|J0|R1        3.2187   | +0.0061  +0.19x  [-0.021,+0.030]  | (comparator)             | +0.004 (+0.11x) / +0.018 (+0.26x)
    vqe|s1|NONE|J0|R1        3.2241   | +0.0115  +0.45x  [-0.004,+0.032]  | (comparator)             | +0.012 (+0.41x) / +0.006 (+0.25x)
    vqe|s0|REAL|J1|R1        3.2182   | +0.0056  +0.21x  [-0.016,+0.026]  | -0.0005 -0.03x           | +0.006 (+0.22x) / +0.000 (+0.01x)
    vqe|s1|REAL|J1|R1        3.2251   | +0.0125  +0.46x  [-0.005,+0.025]  | +0.0010 +0.04x           | +0.015 (+0.48x) / -0.003 (-0.09x)
    vqe|s0|REAL|J0.1|R1      3.2190   | +0.0063  +0.20x  [-0.019,+0.033]  | +0.0003 +0.02x           | +0.005 (+0.14x) / +0.014 (+0.21x)
    vqe|s0|REAL|J3|R2        3.5833   | +0.3706  +1.58x  [+0.186,+0.507]  | +0.3645 +1.54x           | +0.510 (+2.10x) / -0.465 (-0.97x)
    vqe|s0|PERM|J0.1|R1      3.2163   | +0.0037  +0.12x  [-0.019,+0.028]  | -0.0024 -0.22x           | +0.002 (+0.05x) / +0.015 (+0.22x)
    vqe|s0|PERM|J0.3|R3      3.4822   | +0.2696  +1.36x  [+0.113,+0.400]  | +0.2635 +1.33x           | +0.381 (+1.86x) / -0.399 (-0.92x)
    vqe|s0|PERM|J1|R3        3.4624   | +0.2498  +1.28x  [+0.122,+0.368]  | +0.2437 +1.24x           | +0.365 (+1.82x) / -0.443 (-1.02x)
    vqe|s0|RAND|J3|R2        3.5823   | +0.3696  +1.58x  [+0.187,+0.506]  | +0.3636 +1.55x           | +0.506 (+2.09x) / -0.449 (-0.92x)
    vqe|s1|REAL|J0.1|R3      3.4921   | +0.2795  +1.35x  [+0.133,+0.369]  | +0.2680 +1.32x           | +0.401 (+1.90x) / -0.449 (-0.89x)
    vqe|s1|PERM|J0.1|R3      3.4987   | +0.2861  +1.40x  [+0.133,+0.389]  | +0.2746 +1.36x           | +0.406 (+1.95x) / -0.436 (-0.89x)
    vqe|s1|PERM|J0.3|R3      3.4824   | +0.2698  +1.34x  [+0.133,+0.363]  | +0.2583 +1.30x           | +0.386 (+1.87x) / -0.430 (-0.91x)
    vqe|s0|REAL|J0.1|R3      3.4949   | +0.2823  +1.43x  [+0.135,+0.419]  | +0.2762 +1.39x           | +0.385 (+1.86x) / -0.336 (-0.77x)
    vqe|s0|NONE|J0|R3        3.4912   | +0.2786  +1.42x  [+0.133,+0.412]  | +0.2725 +1.38x           | +0.383 (+1.87x) / -0.348 (-0.80x)
    vqe|s1|NONE|J0|R3        3.5290   | +0.3164  +1.53x  [+0.173,+0.420]  | +0.3049 +1.51x           | +0.442 (+2.12x) / -0.436 (-0.88x)
    gs|s-1|REAL|J1|R3        3.2100   | -0.0026  -0.09x  [-0.031,+0.022]  | -0.0087 -0.30x           | +0.005 (+0.15x) / -0.049 (-0.84x)
    gs|s-1|REAL|J1|R2        3.2068   | -0.0058  -0.13x  [-0.026,+0.018]  | -0.0119 -0.23x           | +0.004 (+0.08x) / -0.064 (-0.53x)

- R1 (the deployed readout) at J = 0, 0.1 and 1, both seeds, PERM included: within
  0.0125 A of production (0.46x MDE at most) and within 0.0024 A of its J = 0 twin
  (0.22x at most). The tail is an E-prefix at every J (section 1.4), so R1 can only move
  through m, and on the chain it does not measurably.
- R2 and R3 at every projected J, graph and seed: +0.25 to +0.37 A WORSE than production
  (1.28x to 1.58x, fold CI above zero on every cell), exactly as at J = 0
  (+0.28 / +0.32). The near-uniform state's most probable 75 is a dispersed set.
- F1 does not fire. The one cell that beats its own J = 0 twin beyond MDE,
  `vqe|s1|REAL|J0.1|R3` (-0.0369, -1.24x, 5/5), is a readout +0.279 A above
  production climbing 0.04 A toward it (S28-L2(a)); its seed-0 twin reads
  +0.0037 (+0.20x) and its PERM control moves -0.0303 against ITS twin
  (REAL - PERM -0.0067, -0.18x). Both-seeds and PERM clauses fail on their own. F3 / F4 silent.
- Best-of-nine (S28-L22): the nine 0.7x screen cells mix two R1 cells at ~0 with two R2 and
  five R3 cells at +0.24 to +0.36, so `best_of_k_within` over the nine "transfers" the R1-vs-R3
  readout gap (49% of a -0.461 oracle), a WORSE effect already measured at J = 0.
  Within the R1 pair k_eff 1.34 and the split-half transfer is -0.0007 A (0.06x the R1
  MDE); within the seven R2 / R3 cells NOT A SIGNAL (11%). The best cell is a control
  (`vqe|s0|PERM|J0.1|R1`, -0.0024).
- FAIL18 vs the 108 (addendum 1 item 14(a)): R1 nothing on either stratum; R2 / R3 worse on
  the 108 by +0.37 to +0.51 A (1.8x to 2.1x) and on the other side on FAIL18 (-0.34 to -0.47 A,
  0.8x to 1.0x MDE at n = 18): S27 L9's regime read back (a less contracted average is closer
  where production fails by 6 A), a property of the readout present at J = 0, sign only at
  n = 18, and not a native-free switch (S27 L9, S28-L6).
- The eigensolver (Hamiltonian quality) on the chain: GS REAL J = 1 R3 -0.0026 (-0.09x) and R2
  -0.0058 (-0.13x) vs production: it re-selects the distogram's set and emits production.

### 2.3 The three-way split's middle leg, measured: the J = 3 shortfall is basin selection, not expressivity (S28-L41; S28-L26's request). DEMONSTRATED (property).

`s28_B_represent.json` (job `s28B_represent`, 12 trainability targets, the 9-qubit endpoint
register, REAL graph, J = 3; theta fitted to maximise (g . psi(theta))^2 with g the exact
ground state, Adam lr 0.15, 16 starts, 80 and 400 iterations; the re-run VQE F equals
`s28_B_rows.jsonl` bit-for-bit on 24/24 cells). The ground state ITSELF is not representable:
best overlap 0.80 to 0.88 (median 0.851), reached from every start at 80 iterations
and unchanged at 400 (the 27-parameter state stays at PR ~465 against the ground state's
~301). But the cap does not bind on F: the fitted state has F below the ground state's on
10/12 targets (median -7.409 vs -7.327; the ground state optimises <H>, not F) and below
the VQE's reached F on 24/24 cells. The circuit holds a state with F 2.8 below what the
optimiser reached on the failing cells: an optimisation shortfall, and its form is basin
selection. The trained VQE state at J = 3 is BIMODAL: 7/24 cells fully sign-coherent
(coherence 1.000, overlap 0.78 to 0.86 with the ground state, F at or below the fitted
optimum) and 17/24 fully incoherent (coherence <= 0.02, F -4.56, the J = 0 value), none in
between; on the full 126 (`s28_B_rows.jsonl`, REAL J = 3) 43/126 (seed 0) and 32/126 (seed 1)
cells are coherent, no cell in (0.1, 0.5], 10 targets coherent on both seeds, 61 on neither.
The sign-aligned VQE state |psi_vqe| IS representable (overlap median 0.986), so the wall
between the basins is the sign structure. THIS CORRECTS 1.4 and S28-L21: "the circuit
collects a third of its bound at J = 3" is a mixture (a third of the cells collect all of it,
the rest none), not a partial alignment (lane D's RETRACTIONS_S28 R4). Where the coherent
basin is found the endpoint does not measurably move (point cloud, seed 0, R1 at J = 3 minus
J = 0: +0.009 on the 43 coherent cells, 0.26x its MDE, against -0.003 on the 83 others; the
difference inside a random-subset null, p 0.275, S28-L43): finding the hopping ground state's
basin makes no measurable difference to the endpoint.

---

## 3. WHAT DAMAGED MY OWN EXPECTATIONS

- I registered "WORSE" for the hopping arms on the deployed readout, by the S27 consistency
  mechanism. On the point cloud they are not worse; they are identical to J = 0 within
  0.013 A, because the circuit does not train the hopping term at this width (1.2% of the
  gradient magnitude): the consistency mechanism never acts through the VQE. It acts through
  the eigensolver, and there at J = 1 it re-selects the distogram's own top-75 rather than the
  pool's mode.
- I expected the -1.84 slope to be a statement about off-diagonality. It is a statement about
  the spectrum: a diagonal observable with the same spectrum decays at the same rate, and a
  spread-spectrum graph halves the rate.
- I expected a spread spectrum to make the hopping trainable at J = 1. It did not: the rate
  halved but the ratio to the diagonal terms at n = 9 stayed above 100x, and the circuit
  collected 3 to 23% of its bound. The bright line I registered (-1.0 per qubit) was straddled
  by 0.03 either way; I report it as not cleared.
- The J = 3 rung I added as "the delocalised limit" is where the circuit finally collects a
  third of the hopping, and it is the rung where R1 drifts up (m 74 -> 79): hopping, when it
  acts, acts as an extra entropy.
- I read "collects a third of its bound" as a partial alignment. It is a mixture: a third of
  the (target, seed) cells are fully sign-coherent and the rest are the J = 0 state; which
  basin a run lands in is a property of the draw, not the target (10 of 126 targets coherent
  on both seeds, 61 on neither).
- I expected the 1.8 gap to the exact ground state at J = 3 to be at least partly expressivity
  (the circuit cannot hold a PR-300 state: max overlap 0.85). It is not, on the objective that
  matters: a representable state has F below the ground state's on 10/12 targets. The circuit
  holds the answer to its own objective and Adam finds it on a third of its draws.

## 4. WHAT I DID NOT DO AND WHY

- No sigma other than the median (addendum 0: cosmetic variant). No J beyond 3 and no finer
  grid (the grid is already an order statistic that does not transfer).
- No AMBER. No projection of arms below 0.7x MDE other than J = 0, the best J and the GS
  comparators (18 arms, `s28_B_chain_arms.txt`).
- The 7-qubit production selector (top-128, VQE_LFO table) was not re-run: the S25 suite's
  9-qubit setting is the registered instrument and the tail-reading theorem's scope is
  identical.
- The B2 endpoint ran at the addendum-2 scope only (section 5); no other k or J, no RAND
  graph, no chain beyond the two cells the 0.7x rule named.
- The representability fit was run on the 12 trainability targets only (S28-L26's request);
  the bimodality count is on all 126 from the endpoint rows.
- The k = 5 B2 replication is dropped on the endpoint (disconnected graphs on some targets).

---

## 5. B2 -- THE SPREAD-SPECTRUM (kNN) GRAPH'S ENDPOINT AT THE ADDENDUM-2 SCOPE (S28-L46). REFUTED as registered; the prior (WORSE or null) held; B2 CLOSED.

`s28_B2_rows_k10.jsonl` (1,134 rows = 126 x 9; job `s28B2_run`), `s28_B2_chain_rows_k10.jsonl`
(252 = 126 x 2; `s28B2_chain`), `s28_B2_summary.json`. Scope: k = 10 (connected on 126/126, no
DEGENERATE ground state), J in {0, 3}, graphs REAL and PERM, seeds 0 and 1, VQE R1 / R3, GS R3;
the J = 0 rows are bit-identical to the Gaussian run's on 126/126, both seeds. The coherence-class
split was registered as prereg addendum 3 BEFORE the job (S28-L43's condition).

- Deployed readout R1 at J = 3: within 0.0105 A of J = 0 (0.42x MDE at most) and within
  0.0202 A of production (0.80x at most, the PERM cell; REAL 0.37x). Null.
- R3 at J = 3: +0.28 to +0.29 A WORSE than production on the point cloud (1.32x to
  1.42x, fold CI above zero), as at J = 0; on the built chain (the two REAL cells, projected
  because they reached 0.7x MDE in the WORSE direction) +0.3062 / +0.3125 A (1.45x / 1.45x,
  5/5) and within 0.35x MDE of their J = 0 R3 twin.
- REAL vs PERM (F3): -0.0103 (-0.31x), +0.0115 (+0.19x), +0.0052 (+0.18x), +0.0067 (+0.11x): silent.
- The kNN eigensolver is the predicted tight cluster (PR 53, m 14.5, Jaccard 0.577 with the
  DIS top-75) and emits production within MDE (+0.0440, +0.54x); F2 (VQE R3 vs GS R3): the
  eigensolver on the better side by 0.24 A, fold CI above zero, 0.99x / 0.95x (sign only).
- Classes (addendum 3): the circuit lands in the sign-coherent basin on 49 / 33 of 126
  cells (Gaussian 43 / 32), 1 / 3 mixed; the incoherent class is not the J = 0 state on
  this graph (hopping 0.111 at coherence 0.001, F -4.80: share and coherence
  separate on a spread spectrum). No class contrast of R1 leaves its MDE; R3 is the same WORSE
  readout on both classes: landing in the coherent basin makes no measurable difference.
- S28-L29's "collects 44 to 68% of its bound" is a mixture (S28-L43 / R4), stated here as counts.

What is new against S28-L41 and what it does not change: a different H in the one respect the
mechanism named (lambda_2 / lambda_1 0.98 vs 0.14), a ground state that is the opposite object
(PR 53 cluster vs PR 305 typicality mode), a circuit that departs further from J = 0 (TV 0.29 /
0.35 vs 0.16); none of it reaches the endpoint. Nothing further is licensed at this scope.
