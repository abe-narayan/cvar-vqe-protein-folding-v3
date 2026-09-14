# LANE Q -- QUANTUM. SPRINT 26. FINDINGS.

Pre-registrations: `s26/PREREG_A1.md` to `PREREG_A4.md`. Tournament entries:
`s26/IDEA_l17_target_dependent_hamiltonian.md`, `s26/IDEA_trainability_paper.md`,
`s26/IDEA_product_state_optimum.md`. Paper outline: `s26/TRAINABILITY_PAPER_OUTLINE.md`.
Code: `s26/q_adapt.py`, `s26/q_dla.py`, `s26/q_var.py`, `s26/q_mde.py`, `s26/q_tests.py`.
Tiers as in S12 to S25: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN.
Every number carries its artefact path. No native was read by this lane before the phase
gate; the `--label` phase of the harness refuses to run without "PHASE 0 SIGNED OFF".

BASIS NOTICE. Two readouts appear and never share a column: the built chain
(`rmsd_q_synth`, the production projection of a weighted coordinate average; the production
result on this basis is 3.2148 A for the top-75 uniform arm and 3.2280 A for the deployed
quantum synthesis, `s26/results/q_mde_reference.json`) and the s8-instrument selection
(`rmsd_vqe_sel`, one candidate; 3.3135 A). The point cloud is recorded and never compared.

---

## 0. PRE-SIGN-OFF OBSERVATIONS (synthetic energies at n <= 7, and the one allowed real-target probe)

### 0.1 The deployed selector's optimisation target is a product state. DEMONSTRATED.

`E = _zrank(sc[top[:128]])` with `top` already in score order, so E is affine in the register
index j (E_j = a (j+1) + b up to tie-averaging), and j = sum_q 2^(6-q) b_q. Hence
exp(-E_j/T) = prod_q exp(-a 2^(6-q) b_q / T) x const: the Gibbs state of the deployed
Hamiltonian at alpha = 1 is a product distribution over the 7 bits.

    KL(Gibbs || product of its marginals), nats
      ideal ladder n = 4, 7, 10; T = 0.1, 0.3, 1.0      < 1e-12 (7.2e-17 at n=7, T=0.3)
      1A13, real E (with ties), T = 0.3                 5.4e-5        s26/results/probe/1A13.json
      1A13, zraw / asinh / soft variants                0.042 / 0.025 / 0.038
    (`s26/q_tests.py` test_gibbs_state_of_the_deployed_ladder_is_a_product_state; the probe record)

A single RY layer (7 parameters, no entangler) reaches it: KL(p || Gibbs) = 1e-4 after
L-BFGS-B (arm `adaptV_lbfgs_zrank_P7`, `s26/results/probe/1A13.json`), against 0.903 nats for
the deployed circuit at 50 Adam steps and 0.884 at 750 (arms `fixed_zrank_it50`,
`fixed_zrank_it750`, same file). S25 measured the 0.902 nats and called the trained state
"broader than optimal"; the cause is that the CNOT layers stand between a product input and a
product target. qubit-ADAPT sees it: at alpha = 1 neither pool selects an entangling operator
and growth stops by its gradient criterion at P = 8 (pool V) or P = 10 (pool L2).

What this does NOT say: nothing about the endpoint (A1 measures it), and nothing about the
alpha = 0.25 folds, where the CVaR objective is not the Gibbs functional.

### 0.2 At alpha = 0.25 ADAPT with the 2-local pool reaches a lower free energy than the fixed ansatz. DEMONSTRATED on the ideal ladder.

    ideal ladder n = 7, alpha = 0.25, T = 0.3, seed 0
      fixed ansatz, 50 Adam steps (deployed)      F = -2.9850
      fixed ansatz, 750 Adam steps                F = -2.9889
      ADAPT L2, 14 distinct 2-local strings       F = -3.0392 (Adam last), -3.0344 (best), -3.0350 (L-BFGS)
      ADAPT V                                     F = -2.9843 to -2.9846
    (console record of the 00:40 diagnostic; to be re-derived into s26/results/q_dla.json by A2's
     ADAPT-set sweep, which stores F_final per cell)

A lower F is an objective-level statement. S25 section 6.3 measured that the endpoint cannot
tell a 45%-of-mass distributional difference; A1 tests whether this one is visible.

### 0.3 The deployed Adam is not a converging re-optimiser for ADAPT; the fix is pre-registered. DEMONSTRATED.

A fresh Adam at lr 0.15 moves every angle by 0.15 rad on its first step whatever the gradient:
on 1A13-shaped synthetic E the free energy jumps from -2.6238 to -2.5227 and oscillates for
30 steps; its last iterate is not stationary (|g| ~ 7e-2), the ADAPT stop criterion
(||g_pool|| < 1e-3) can never fire, and the same single-qubit Y is re-selected up to 13 times in
a row (a repeated rotation merges with the previous one and adds no expressivity).
Measured on the ideal ladder at n = 7, both VQE_LFO cells, both pools:

    re-optimiser        F reached (alpha=1)   consecutive repeats   |g| at the end
    Adam, last iterate  -2.4529 to -2.4531    11 to 12              6.6e-2 to 7.6e-2
    Adam, best iterate  -2.4536               12 to 13              7.4e-3 to 1.2e-2
    L-BFGS-B            -2.4537 (= Gibbs)     0 (stops by eps)      8e-9 to 1.4e-8

The pre-registered primary is Adam with the best iterate (the deployed optimiser's arithmetic,
one documented deviation), and L-BFGS-B is the secondary. The repeat count is recorded per run
and reads as a non-convergence diagnostic.

### 0.4 The dynamical Lie algebra at n <= 6 (synthetic, `s26/results/q_dla_smoke_n456.json`). DEMONSTRATED.

    fixed ansatz  n=4: L=1 4,  L=2 120 = so(16)
                  n=5: L=1 5,  L=2 496 = so(32)
                  n=6: L=1 6,  L=2 510, L=3 1023, L=4 2016 = so(64)
    pools V, G:   36, 136, 528 at n = 4, 5, 6  (= dim so(2^(n-1)+1))
    pool L2:      120, 496, 2016 = so(2^n)
    symbolic == numeric (dense SVD, rtol 1e-10) on 8 of 8 cells at n = 4, 5

Depth 1 is abelian, which is the algebraic form of S13's exact I/4 metric at depth 1. The
n = 7 row is A2 (running); predictions are in `s26/PREREG_A2.md`.

### 0.5 The harness reproduces the production quantum arm bit-for-bit. DEMONSTRATED.

`ca` and `q_ca` from `s26/q_adapt.py --build` equal `bench_results/cache/464a0ddb5f283e04/
<pdb>.json` at max |diff| = 0.0, and the selection index is equal, on 1A13 (fold 4, alpha 1),
1A1P (fold 1, alpha 0.25), 1CS9 (fold 2, alpha 0.25), 2MK7 (fold 4)
(`s26/results/probe/1A13.json`, `s26/results/bitforbit/{1A1P,1CS9,2MK7}.json`, key
`cache_check`). The entropy of the deployed state on 1A13 reproduces to every printed digit
(5.660897219523758 bits).

### 0.6 Memory and time (the honest estimates, from `s26/jobs_done/*.json`)

    q_tests_synthetic5      17/17 tests     10 s     peak RSS 0.132 GB
    q_probe_1A13            16 arms         90 s     peak RSS 0.338 GB   (75.7 s in the harness, 73.2 s projection)
    q_probe_bitforbit3      3 targets       45 s     peak RSS 0.379 GB
    q_dla_smoke             n <= 6         156 s     peak RSS 0.384 GB
    q_var_smoke             n = 4, 6        10 s     peak RSS 0.049 GB

### 0.7 Reference MDEs for the endpoint pre-registrations (`s26/results/q_mde_reference.json`)

    built chain  rmsd_q_synth - rmsd_u_synth        -0.0135  SE 0.0342  MDE 0.0958  62W/64L
    built chain  rmsd_q_synth - rmsd_arm            +0.0133  SE 0.0181  MDE 0.0508
    selection    vqe_a1.0_T0.3 - boltz_T0.3         -0.0302  SE 0.0450  MDE 0.1261
    selection    rmsd_vqe_sel - rmsd_vqe_sel_uniform -0.0308 SE 0.0590  MDE 0.1652
    selection    vqe_LFO - argmin                   -0.1405  SE 0.0732  MDE 0.2051  0.68x (S25's number, reproduced)

A note for the presenter: the production headline 3.2148 A (`bench_results/cache/
1fc9f2dcf489e2fb`, `quantum = None`) never passes through the VQE; the quantum synthesis
is a separate arm at 3.2280 A (`464a0ddb5f283e04`), +0.0133 A, 0.26x MDE, NOT MEASURED.

---

## 1. A1 -- ADAPT-VQE vs THE FIXED ANSATZ (ENDPOINT). COMPLETE. NULL AT THE REGISTERED THRESHOLD.

Pre-registered in `s26/PREREG_A1.md` (addendum records the outcome). Ledger L68 carries the
five ST.fmt blocks verbatim. Artefacts: `s26/results/a1/<pdb>.json` (126, every one
bit-for-bit against the production quantum cache: `cache_check` ca and q_ca 0.0, selection
index equal, 126/126), `s26/results/a1_stats.json`, `s26/logs/a1_stats.log`;
`s26/jobs_done/a1_build.json` (7,782 s for 93 targets, peak RSS 0.383 GB),
`a1_label.json` (25 s, 0.36 GB). Basis: built chain on both sides; the selection readout is
the named secondary.

### 1.1 The primaries

    built chain, vs fixed_zrank_it50 (3.2280 A), n = 126, paired
      adaptL2_adam_best_zrank_P21   -0.0138  SE 0.0210  MDE 0.0588  0.23x  fold [-0.0705,+0.0442]  3/5  58W/68L  conc. pctile 0.498
      adaptV_adam_best_zrank_P21    -0.0222  SE 0.0217  MDE 0.0608  0.36x  fold [-0.0854,+0.0424]  3/5  61W/65L  conc. pctile 0.509
    "ADAPT is null" (within +-0.5x MDE) FIRED on both. "ADAPT helps" did not fire.
    Power: resolution 0.059 to 0.061 A; Gelman-Carlin power at the observed effect 0.10 / 0.18,
    Type-M 3.7 / 2.4; the whole quantum synthesis vs the classical arm is +0.0133 A
    (`s26/results/q_mde_reference.json`), four times smaller than the resolution. NULL at the
    registered threshold, UNDERPOWERED below 0.06 A.

### 1.2 Every arm (built chain; all fold CIs span zero; 3-4/5 folds)

    adaptL2 adam_best  P7 / P14 / P21   -0.0128 (0.21x) / -0.0192 (0.33x) / -0.0138 (0.23x)
    adaptL2 lbfgs      P7 / P14 / P21   -0.0206 (0.34x) / -0.0187 (0.32x) / -0.0194 (0.33x)
    adaptV  adam_best  P7 / P14 / P21   -0.0128 (0.21x) / -0.0195 (0.32x) / -0.0222 (0.36x)
    adaptV  lbfgs      P7 / P14 / P21   -0.0206 (0.34x) / -0.0244 (0.39x) / -0.0245 (0.40x)
    fixed 750 steps -0.0035 (0.16x)   gibbs_T +0.0088 (0.11x)   uniform128 +0.0135 (0.14x)
    randH_fixed +0.0230 (0.24x, 4/5)   randH_adaptL2 +0.0337 (0.34x, fold CI [+0.012,+0.058], 5/5)
    alpha = 1 subset (n=78):  L2 P21 -0.0223 (0.24x)   V P21 -0.0254 (0.28x)
    alpha = 0.25 subset (n=48): L2 P21 +0.0000 (0.00x)  V P21 -0.0169 (0.27x)

Selection readout: all twelve ADAPT arms -0.039 to -0.048 A, 0.42x to 0.47x MDE, fold CIs
excluding zero on all twelve, 45 to 59 exact ties (58 of 126 targets select the identical
candidate under ADAPT-L2-P21 and the fixed circuit). Not a result by the standing rule; the
shape of S25's `VQE_LFO - argmin` (0.68x, 5/5). The 7-parameter arms carry the same effect
as the 21-parameter ones on both bases: whatever the direction is, it is not expressivity.

### 1.3 The property half on the real targets (no native)

    KL(p || Gibbs), 78 alpha = 1 targets   fixed 0.9027 (max 0.984; S25's 0.902)   ADAPT 0.0002 (max 0.0009), both pools, both optimisers
    KL(Gibbs_zrank || product)             mean 1.4e-4, max 7.9e-4 over 126 targets
    L-BFGS growth at alpha = 1              no operator selected, 78 of 78 targets (stopped by eps at P = 7)
    ADAPT L2 at alpha = 0.25                mean 13.7 distinct 2-local strings; 38 distinct sequences on 48 targets
    KL(p || Gibbs), 48 alpha = 0.25 targets fixed 2.67, ADAPT 1.59 (Gibbs is not the CVaR optimum there; recorded, not interpreted)

Section 0.1's product-state finding holds on every real target: the deployed selector's
alpha = 1 target is a product state, seven RY angles reach it, the fixed circuit misses it
by 0.90 nats, and reaching it exactly moves the built chain by -0.02 A, a third of the MDE.

### 1.4 What damaged my expectations here

I expected the alpha = 0.25 subset (where the objective is not the Gibbs functional and
ADAPT-L2 reaches a lower F) to be where any signal lived. It is the one cell that is exactly
zero (+0.0000, fold CI [-0.008, +0.008], n = 48): a lower free energy on the CVaR objective
emits the same structure. The direction that does appear is on the alpha = 1 targets, where
ADAPT differs from the fixed circuit only by reaching the product optimum.

### 1.5 Verdict carried to `s26/PROPOSAL_A.md`: REPLACE.

## 2. A2 -- THE DYNAMICAL LIE ALGEBRA. COMPLETE. DEMONSTRATED (exact; property, no native).

`s26/q_dla.py` -> `s26/results/q_dla.json`; `s26/jobs_done/a2_dla.json`: 85.3 s, peak RSS
0.479 GB. Pre-registered in `s26/PREREG_A2.md`; the addendum there records which
predictions held.

### 2.1 The fixed ansatz: maximal algebra from depth 2 at the deployed width

    dim(DLA), RY / CNOT chain + ring, by n and L         so(2^n)      su(2^n)
      n = 7   L = 1     7                                   8128        16383
              L = 2..6  8128  = so(128)                     (the deployed cell is L = 3)
      n = 4   L = 1  4;  L >= 2  120  = so(16)
      n = 5   L = 1  5;  L >= 2  496  = so(32)
      n = 6   L = 1  6;  L = 2  510;  L = 3  1023;  L >= 4  2016 = so(64)
      n = 8   L = 1  8;  L >= 2  32640 = so(256)
      n = 9   L = 1  9;  L = 2  32766;  L = 3  65535;  L >= 4  130816 = so(512)
      n = 10  L = 1  10; L >= 2  523776 = so(1024)
      n = 11  L = 1  11; L >= 2  2096128 = so(2048)
    symbolic == numeric (dense SVD, rtol 1e-10) at n = 4, 5, L = 1..6: 12 of 12
    left and right conjugation conventions: same dimension, 12 of 12
    every closure inside the odd-Y (real) set: yes, at every cell

Method, stated exactly. Symbolic route: each generator is a Pauli string (x, z) bit pair;
two strings anticommute iff popcount(x & z') + popcount(z & x') is odd; the commutator of two
anticommuting strings is the string (x ^ x', z ^ z') up to a non-zero scalar; the closure is
the breadth-first set of strings reached by commuting generators into already-reached
strings; dim(DLA) is the size of that set (distinct strings are linearly independent and
right-normed nested commutators span the algebra). Numeric route: the same generators as dense
2^n x 2^n real antisymmetric matrices -iP/2; the nested commutators [g, b] of every generator
g with every basis element b are stacked as flattened vectors with the current basis; rank =
the number of singular values above 1e-10 times the largest singular value (relative
tolerance 1e-10, `numpy.linalg.svd`); the basis is rebuilt from the right singular vectors
and the step repeats until the rank stops growing. The two routes agree on all 12 (n, L)
cells at n = 4, 5 (`s26/results/q_dla.json`, key `numeric`, field `agree`).

Two different numbers in `s26/logs/a2_dla.log` are two different generator sets, and a reader
should not average them: 8128 is dim(DLA) of the FIXED ansatz's 21 conjugated RY generators
at n = 7 (L = 3, and already at L = 2); 1025 is dim(DLA) of the 21 strings that ADAPT
SELECTED with pool L2 at alpha = 0.25 (7 initial Y_q plus 14 grown 2-local strings), 12.6%
of so(128). The fixed ansatz reaches the whole real algebra by depth 2; the grown ansatz at
P = 21 spans an eighth of it.

Depth 1 is the abelian algebra of the n commuting RY generators; that is the algebraic form
of S13's exact I/4 metric at depth 1 (`s13/results/geo_metric_arm1.json`). From depth 2 the
algebra is the whole of so(2^n) at every width except n = 6 and n = 9, where it takes depth 4
and the intermediate dimensions equal 2 dim su(2^(n-2)) and dim su(2^(n-1)). That the
exceptions are the multiples of 3 is an observation from two widths: HYPOTHESIS, not
explained.

My prediction H2b (a proper subalgebra at n = 7, L = 3, guessed at 4095) is FALSIFIED. The
deployed ansatz is controllable on the real sphere already at depth 2.

### 2.2 The pools

    V and G (Tang 2021, 2n-2 strings)   36, 136, 528, 2080, 8256, 32896   n = 4..9
                                        = dim so(2^(n-1)+1) at all six widths (predicted for
                                          n = 7, 8, 9 in the prereg: held)
    L2 (1-, 2-local, odd Y)             so(2^n) at n = 4..9

A "complete" pool in Tang's sense (overlap-matrix rank 2^n - 1) generates an algebra
transitive on the real sphere with about a quarter of the dimension of so(2^n). The A1 arm
with pool V is confined to it by construction; the L2 arm is not.

### 2.3 The ADAPT-selected sets on the ideal ladder, n = 7

    alpha = 1,    T = 0.3   both pools, both optimisers: abelian, dim 7 at every step;
                            L-BFGS selects nothing and stops at P = 7 (the product target).
    alpha = 0.25, T = 0.3   V: 7 -> 16 (one string, IIIYZZZ);  L2: 7 -> 1025 at P = 21,
                            12.6% of so(128), for both optimisers.

### 2.4 What it means for the S25 slopes

The algebra at the deployed cell is maximal, so nothing algebraic protects the ansatz from a
plateau. Larocca et al. (2022) and Ragone et al. (2024) give Var ~ 1/dim(g) once the circuit
is a 2-design over exp(g); with dim(g) = 8128 at n = 7 and growing as 4^n / 2, that is the
2^-n rate S13 measured at depth 8 (decay base 0.504, `s13/results/geo_kernel.json`). The S25
"no exponential plateau at n = 4..13" (`s25/results/q_plateau.json`) is a statement about
depth 3, P = 3n against dim so(2^n), and is to be quoted with "at depth 3" attached. It is
not evidence of a favourable algebra. The small-DLA route of Cerezo et al. 2025 does not apply;
the circuit is simulable because n = 7 (S21 L4), not because of its structure.

## 3. A3 -- target-dependent Hamiltonians. PENDING the phase gate (endpoint half).

Pre-registered in `s26/PREREG_A3.md`.

## 4. A4 -- GRADIENT VARIANCE OF GROWN CIRCUITS. COMPLETE. DEMONSTRATED (property, no native).

`s26/q_var.py` -> `s26/results/q_var.json`; `s26/jobs_done/a4_var.json`: 2,765 s, peak RSS
0.08 GB. Pre-registered in `s26/PREREG_A4.md` (addendum records which predictions held).
Figure `s26/figures/a4_variance_slopes.png`. Ledger L35.

### 4.1 The gate: S25 reproduced exactly

`s25.q_plateau.measure(7, 3, alpha, T, 250, seed=1007, init_sd=0.6)` returns the n = 7
`var_g0` of every cell of `s25/results/q_plateau.json` at relative deviation 0.0 (bit-identical
on all five cells).

### 4.2 The side-by-side table (matched P = 3n rows only for the grown circuits)

    fitted log2 Var[dF/dtheta_0] per qubit, n = 4..13, theta ~ N(0, 0.6^2), draws 250/250/250/250/200/120/80
      cell                    fixed (S25)    grown V (rows)     grown L2 (rows)
      alpha=1,    T=0         -0.649         -0.079 (7)         +0.006 (7)
      alpha=0.25, T=0         -0.252         degenerate (0)     degenerate (0)
      alpha=0.10, T=0         -0.047         degenerate (0)     degenerate (0)
      alpha=1,    T=0.3       -0.311         +0.035 (5)         -0.008 (7)
      alpha=0.25, T=0.3       -0.243         -0.246 (6)         -0.302 (7)
    unmatched rows (eps stop before 3n): V alpha=1 T=0.3 at n = 12 (P 13), 13 (P 15);
                                          V alpha=0.25 T=0.3 at n = 12 (P 13)
    grown/fixed variance ratio at matched n: 1.63 at n = 4 (alpha=1, T=0.3, both pools),
                                             2.5 to 370 on the other 31 matched rows

The slopes printed at the end of `s26/logs/a4_var.log` include the early-stopped rows and
differ from the table above for the two T = 0 cells (-0.733, +0.065 printed) and slightly for
V at T = 0.3; the matched-only values are the pre-registered quantity.

### 4.3 What the grown circuits are, which is the finding

At alpha = 1 (T = 0 and T = 0.3) the grown circuits hold 1 to 3 distinct operators; the other
12 to 25 selections are consecutive repeats of one single-qubit Y (`n_distinct_ops`,
`consecutive_repeats` per row in the artefact). They are product circuits carrying the fixed
ansatz's parameter count, so their gradient variance does not decay with n (slopes -0.08 to
+0.04) and is 2 to 370 times the fixed ansatz's. That is the trivial regime, and section 0.1
says why: the target at alpha = 1 is a product state. At T = 0 with alpha < 1 growth stops at
P = n on every width (the collapse; all gradients vanish at a basis state) and the variance
is 0 to 2.4e-2: degenerate, as pre-registered. The one non-trivial family, alpha = 0.25,
T = 0.3, pool L2, grows 4 to 21 distinct 2-local strings and decays at -0.302 per qubit
against the fixed ansatz's -0.243, equal within the error of a 7-point slope fitted to
variances with 9 to 16% relative SE; its variance is 4.1 to 20.7 times the fixed one's at
matched n, with no trend. The L2 alpha = 1 curve is erratic in n (0.089, 0.102, 0.372, 0.092,
0.410, 0.097, 0.086) because which operators are grown changes with n; that is a property of
the grown family, not sampling error.

### 4.4 Predictions, as measured

    H4a  T = 0.3 grown slopes in [-0.3, 0]   3 of 4 inside; +0.035 (V, alpha=1) and -0.302 (L2,
                                             alpha=0.25) outside by 0.035 and 0.002; the
                                             falsifier (below -0.5) did not fire
    H4b  grown/fixed ratio > 2 at every matched n   failed on 1 of 32 rows (n = 4, 1.63), held on 31
    H4c  T = 0, alpha < 1 degenerate         P < 3n on 14 of 14 rows; "Var < 1e-3" failed on
                                             3 of 7 rows at alpha = 0.25 (2.4e-2, 1.0e-3, 1.2e-3)
    H4d  S25 n = 7 reproduced                exactly

### 4.5 What it means for Proposal A

A large gradient from a product circuit is not trainability (rule 10's mirror). The only grown
family that is not a product circuit decays like the ansatz it would replace. A4 gives Proposal
A no width-scaling argument, and the record's "no exponential plateau at n <= 13" keeps its
scope: depth 3, this ansatz, this spectrum, with the algebra (A2) offering no protection.

---

## 5. PRE-REGISTRATION SUMMARY

    A1  endpoint   fixed vs ADAPT (V, L2; Adam-best primary, L-BFGS secondary), P = 7/14/21,
                   controls fixed_750, gibbs_T, uniform128, randH; built chain primary.
                   "helps": > MDE, fold CI excludes 0, 5/5, replicates seed 1 and reverse order.
                   "null": within 0.5x MDE. Expected null.
    A2  property   dim(DLA) by n and L vs so/su; pools; ADAPT sets; symbolic == numeric.
    A3  property + endpoint   distinct states / sequences under zrank vs zraw/asinh/soft;
                   endpoint zraw and zraw_Tmatch vs zrank. Expected: zraw harmful, Tmatch null.
    A4  property   grown vs fixed log2 Var per qubit at matched P; T = 0 grown cells degenerate.

## 6. WHAT DAMAGED MY OWN EXPECTATIONS

1. I expected ADAPT to grow entangling operators and to be judged on whether entanglement
   helps. On the deployed Hamiltonian at alpha = 1 there is nothing to entangle: the target is
   a product state (section 0.1, then 126 of 126 real targets in 1.3). That was available from
   the definition of `_zrank` and the register index the whole time.
2. I expected the brief's re-optimiser (Adam, lr 0.15) to be usable as written. Its first
   step is a fixed 0.15 rad in every direction; ADAPT's stop criterion never fires under it
   (section 0.3). The prereg carries the best-iterate rule as a documented deviation.
3. I expected the Tang minimal pools to generate so(2^n) (they are "complete"). They generate
   so(2^(n-1)+1) at n = 4..9 (A2, L27). Completeness (overlap-matrix rank 2^n - 1) is weaker
   than controllability.
4. I predicted a proper subalgebra for the fixed ansatz at n = 7, depth 3 (PREREG_A2 H2b,
   guess 4095). It is the whole so(128) from depth 2. Recorded as falsified in L27.
5. I expected the alpha = 0.25 subset, where ADAPT-L2 reaches a lower free energy than the
   fixed circuit, to be where any endpoint signal lived. It is the one cell that is exactly
   zero (+0.0000, fold CI [-0.008, +0.008], n = 48; section 1.4). A lower objective emitted
   the same structure.
6. I wrote "selects no entangling operator on 78 of 78 targets" into L68 and PROPOSAL_A.md
   from the ideal-ladder run without re-reading the real-target records; PR's L73 re-read them
   and the appended count is 60 to 78 of 78, every string multi-qubit, all inert (<= 1.2e-4
   nats). Corrected in L75 and section 9. The error had the direction of my own hypothesis,
   which is the direction the record says errors are hardest to see.
7. I expected the DLA to be the right diagnostic for "what the grown circuit does". The L75
   table shows why it is not: the algebra is a property of the generator SET, and a set of
   inert multi-qubit strings has a large algebra while the state stays a product. The
   per-growth-step DLA on the real records (PREREG_A2 addendum 2, running) is pre-registered
   to show exactly that.

## 7. WHAT I DID NOT DO AND WHY

- No replication of A1 was owed: both primaries fired the "null" falsifier. The seed-1 and
  reversed-order run (PREREG_A1 addendum 2) is a robustness check under L77's extended scope
  and launches after A3 as one governed process; it is not a contract requirement.
- The complex pool L2C is implemented and tested (pennylane cross-check at n <= 6) but is not
  an arm: the deployed amplitudes are real and a complex pool would confound "grown" with
  "complex".
- The L-BFGS-B arms are secondary, as registered; they carry the same endpoint numbers as the
  Adam arms (section 1.2) and stop early on 78 of 78 alpha = 1 targets after inert additions.
- The squared-risk (posterior-mean) functional is not an A3 variant: it changes the candidate
  ORDER and the brief asked for order-preserving Hamiltonians (IDEA_l17 item 5).
- `pytest tests/` was not run by this lane (lane I owns it); `s26/q_tests.py` (17 tests) was
  run under `jobrun` instead.
- No literature was fetched beyond Tang et al. 2021 (for the exact pool definitions); the
  related-work list is the record's own from S13 to S22.
- The A2 exceptions at n = 6 and n = 9 (depth 4 needed; "multiples of 3") are an observation
  from two widths and were not pursued; HYPOTHESIS.
- Nothing on hardware, noise or shot cost: every number is exact simulation, and the outline
  says what a submission would still need.
- A3's endpoint half runs after this file's last edit; its section is written when it lands.

## 8. ARTEFACTS

    s26/results/q_mde_reference.json      reference MDEs from stored artefacts (q_mde.py)
    s26/results/probe/1A13.json           the 1-target probe (all arms, no native)
    s26/results/bitforbit/*.json          the 3-target fixed-only reproduction check
    s26/results/probe_property.json       property summary of the probe
    s26/results/q_dla_smoke_n456.json     DLA at n <= 6 (synthetic)
    s26/results/q_var_smoke_n46.json      A4 smoke (n = 4, 6; S25 n=7 reproduction passed)
    s26/results/q_dla.json                A2 (complete; ledger L27)  -> s26/figures/a2_dla_dimension.png
    s26/results/q_var.json                A4 (complete; ledger L35)  -> s26/figures/a4_variance_slopes.png
    s26/results/a1/<pdb>.json             A1 per-target records, 126 (built, labelled after L33)
    s26/results/a1_stats.json             A1 contrasts (ledger L68, corrected by L75); s26/logs/a1_stats.log
    s26/results/q_dla_a1.json             per-growth-step DLA on the A1 records (running) -> s26/figures/a2_dla_grown_ladder.png
    s26/results/q_var_boot.json           A4 slope bootstrap CIs (running)
    s26/results/a3/<pdb>.json             A3 per-target records (building)
    s26/results/a1s1/<pdb>.json           A1 seed-1 / reversed-order replication (after A3)
    s26/jobs_done/q_*.json, a2_dla.json, a4_var.json, a1_build.json, a1_label.json   peak RSS and wall per job

## 9. CORRECTION (2026-09-13 22:15, ledger L75): sections 0.1 and 1.3 overstated "no operator selected"

Section 0.1 ("qubit-ADAPT sees it: at alpha = 1 neither pool selects an entangling operator
and growth stops by its gradient criterion at P = 8 (pool V) or P = 10 (pool L2)") and section
1.3 ("L-BFGS growth at alpha = 1: no operator selected, 78 of 78 targets") are RETRACTED. On the
1A13 probe the L-BFGS runs appended IYIZZZZ (V) and three YZ-type strings (L2), all multi-qubit;
on the 78 alpha = 1 targets operators were appended on 60 (V) and 68 (L2) of 78 under L-BFGS,
every one multi-qubit, and on 78 of 78 under Adam. They are inert: F moves by at most 1.2e-4
nats (L-BFGS) or 8.6e-4 (Adam), appended angles are at most 0.018 rad (L-BFGS), and the state
stays a product state to KL <= 4.1e-4 (table in ledger L75). The ideal-ladder statements in
sections 2.3 and 4.3 (nothing appended under L-BFGS; abelian selected sets) stand, and the
difference is tie-averaging: on real targets E is not exactly affine, so the RY layer's
residual gradient in a multi-qubit direction can exceed eps = 1e-3. The endpoint numbers and
the verdict are unchanged. This is section 6's list, one entry longer: I read the ideal-ladder
run into the real-target sentence without re-reading the records; PR did re-read them.
