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

## 1. A1 -- ADAPT-VQE vs the fixed ansatz (endpoint). PENDING the phase gate.

Pre-registered in `s26/PREREG_A1.md`. Expected: null (both primaries within 0.5x MDE).

## 2. A2 -- the dynamical Lie algebra. RUNNING (`s26/jobs/a2_dla.json` -> `s26/results/q_dla.json`).

Pre-registered in `s26/PREREG_A2.md`.

## 3. A3 -- target-dependent Hamiltonians. PENDING the phase gate (endpoint half).

Pre-registered in `s26/PREREG_A3.md`.

## 4. A4 -- gradient variance of grown circuits. RUNNING (`s26/jobs/a4_var.json` -> `s26/results/q_var.json`).

Pre-registered in `s26/PREREG_A4.md`.

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
   a product state. That was available from the definition of `_zrank` and the register index.
2. I expected the brief's re-optimiser (Adam, lr 0.15) to be usable as written. Its first step
   is a fixed 0.15 rad in every direction; ADAPT's stop criterion never fires under it. The
   prereg carries the best-iterate rule as a documented deviation rather than pretending.
3. I expected the Tang minimal pools to generate so(2^n) (they are "complete"). At n = 4, 5, 6
   they generate so(2^(n-1)+1), an algebra transitive on the real sphere but far smaller than
   so(2^n). Completeness in Tang's sense (rank 2^n - 1 of the overlap matrix) is weaker than
   controllability.

## 7. WHAT I DID NOT DO AND WHY

- No endpoint (no RMSD to a native) was computed: the phase gate. The harness's `--label`
  refuses without the ledger line.
- No A1/A3 build beyond the 1-target probe and the 3-target fixed-only reproduction check:
  the contract allows 1-target probes before sign-off; the 3-target check reads no native and
  was the brief's own pre-sign-off assertion.
- The complex pool L2C is implemented and tested (pennylane cross-check) but is not an arm:
  the deployed amplitudes are real and a complex pool would confound "grown" with "complex".
- `pytest tests/` was not run by this lane (lane I owns it); `s26/q_tests.py` was run under
  `jobrun` instead.
- No literature was fetched beyond one paper (Tang et al. 2021, for the exact pool
  definitions); the related-work list is the record's own from S13 to S22.
- The re-derivation of 0.2's console numbers into an artefact is left to A2's ADAPT-set sweep
  (it stores F_final per cell); until it lands they are console numbers and are labelled so.

## 8. ARTEFACTS

    s26/results/q_mde_reference.json      reference MDEs from stored artefacts (q_mde.py)
    s26/results/probe/1A13.json           the 1-target probe (all arms, no native)
    s26/results/bitforbit/*.json          the 3-target fixed-only reproduction check
    s26/results/probe_property.json       property summary of the probe
    s26/results/q_dla_smoke_n456.json     DLA at n <= 6 (synthetic)
    s26/results/q_var_smoke_n46.json      A4 smoke (n = 4, 6; S25 n=7 reproduction passed)
    s26/results/q_dla.json                A2 (running)
    s26/results/q_var.json                A4 (running)
    s26/jobs_done/q_*.json, a2_dla.json, a4_var.json   peak RSS and wall per job
