# LANE B -- THE FIRST NON-DIAGONAL HAMILTONIAN. SPRINT 28. FINDINGS.

**DRAFT, PAUSED 2026-09-14 23:30 BY THE USER'S PAUSE ORDER. Section 2.2 (the built-chain
verdict) is NOT written: job `s28B_chain` was at 98/126 targets when the pause came and is
still running; nothing in this file is a verdict on the reporting basis. The file is committed
so nothing is lost; it is completed from `s27/RESUME_S28.md :: B` when the sprint resumes.**

Pre-registration `s27/PREREG_S28_B.md` (base + addendum 0 before any endpoint number;
addenda 1 and 2 for B2). Code `s27/s28_B_hop.py` (the module), `s28_B_train.py` (F5,
checkpointed), `s28_B_rank1.py` (mechanism), `s28_B_split.py` (three-way split), `s28_B_mladder.py`
(T9 decomposition), `s28_B_analyse.py` (statistics), `s28_B_chain_part.py` (optional half-split
chain runner, unused), `s28_B2_knn.py`, `s28_B2_rank1.py` (B2). Tests `tests/test_s28_B.py` (14
pass), `tests/test_s28_B2.py` (6 pass); lane D's `tests/test_s28_D.py` adds the Perron reading.
Results `s27/results/s28_B_*.json|jsonl`, `s28_B2_*.json|jsonl`. Ledger: S28-L8b (F5), S28-L21
(point cloud, intermediate), S28-L25 (B2 trainability), S28-L29 (B2 share and decomposition);
lane D's checks S28-L2, S28-L9, S28-L11, S28-L22, S28-L23. Tiers as in S12 to S25:
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
0.89 to 0.91 (sign coherence 0.01 at J 1, 0.34 / 0.25 at J 3 by seed). The exact ground state
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

### 2.2 Built chain (THE VERDICT). NOT YET WRITTEN.

Job `s28B_chain` (18 arms x 126; `s28_B_chain_rows.jsonl`, 98/126 targets at the pause) is
running. When it lands: `python s27/s28_B_analyse.py` writes `s28_B_summary.json :: chain`
(each arm vs J = 0 R1 on the chain, vs production (S27 `chain_rows.jsonl :: DIS`), FAIL18 /
108); the verdict entry answers S28-L2 (a) to (e), S28-L22 (best-of-nine pricing over the nine
0.7x cells with `best_of_k_within`; no W/L quoted as evidence) and contract addendum 1 (built
chain only; FAIL18 vs the 108; what is new vs S26/S27; the three-way split from 1.4 and 1.5).

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

## 4. WHAT I DID NOT DO AND WHY

- No sigma other than the median (addendum 0: cosmetic variant). No J beyond 3 and no finer
  grid (the grid is already an order statistic that does not transfer).
- No AMBER. No projection of arms below 0.7x MDE other than J = 0, the best J and the GS
  comparators (18 arms, `s28_B_chain_arms.txt`).
- The 7-qubit production selector (top-128, VQE_LFO table) was not re-run: the S25 suite's
  9-qubit setting is the registered instrument and the tail-reading theorem's scope is
  identical.
- The B2 endpoint (k = 10, J = 3 only, prereg addendum 2) was not started: it is gated on the
  Gaussian chain verdict and lane D's check, neither of which exists at the pause.
- The k = 5 B2 replication is dropped on the endpoint (disconnected graphs on some targets).
