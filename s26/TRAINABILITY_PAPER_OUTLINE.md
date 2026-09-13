# TRAINABILITY PAPER OUTLINE (lane Q, Sprint 26, deliverable 11)

Working title: "What a molecular force field's Pauli spectrum does and does not predict: an
exact study of variational trainability with the energy model as the independent variable"

Everything below is tied to an artefact in this repository. Where a figure depends on this
sprint's runs (A2, A4) the artefact is named and the row is marked PENDING until it lands.

## 1. Claims the paper makes (each with its scope condition)

C1. Exact locality theorem in torsion space: under an ideal-geometry builder the CA-CA distance
    d_ij depends on exactly the j-i-1 residues strictly between i and j; agreement 1.0000 over
    2,993 (pair, variable) cells, zero counterexamples; all-atom supports are one residue
    wider and populate the s = 0 and s = 1 shells. Scope: ideal geometry, this builder.
    `s13/results/qarch_locality_geom.json`, `s13/results/qarch_locality_amber.json`,
    `s13/results/walsh_predict.json` (outside-support Walsh variance <= 1.1e-29).
C2. The Pauli-weight spectrum of a raw all-atom force field on unrelaxed builds is the
    signature of a delta spike, not of interaction structure: the ten most extreme
    configurations carry a median 99.6% of the Walsh variance and the mean weight lands on
    m/2 to three digits (Binomial(m, 1/2)). Only monotone rank-preserving conditioning
    removes it; 99th-percentile winsorisation does not. `s13/results/walsh_xval.json`,
    `s13/results/geo_pauli_v1_rawonly.json`, `s13/results/geo_anova.json`.
C3. Conditioned identically, AMBER ff14SB/GBn2 is the higher-weight observable: mean Pauli
    weight 3.015 against Legacy 2.236, AMBER higher on 79 of 79 cells; both are a steric term
    plus rounding error (AMBER nonbonded covariance share 1.000 on 26 cells; Legacy steric
    0.955 at m = 18). `s13/results/geo_pauli.json`, `s13/results/walsh_amber.json`,
    `s13/results/walsh_legacy_k4_m18_6t.json`.
C4. Spectrum times ansatz kernel predicts the gradient variance with no free parameter:
    median measured/predicted 0.996 (Legacy), 0.998 (Legacy-soft and AMBER-soft) over 78
    cell x depth combinations; the one approximation is dropping cross-covariances between
    Pauli strings, and the ratio is its test. `s13/results/geo_pauli.json`.
C5. The kernel is flat in Pauli weight at depth >= 3 (v(n)/v(1) 0.88 to 1.06, no trend in n)
    and locality matters only at depth 1 to 2 (0.862 -> 0.078 over n = 6 -> 14 at depth 2);
    so Cerezo et al. 2021's cost-locality mechanism is inactive at the depth in use, and the
    theorem's assumptions (local 2-design blocks, O(log n) depth) do not hold for this ansatz.
    `s13/results/geo_kernel.json`.
C6. The Fubini-Study metric contains no Hamiltonian (bit-identical across energy models at
    matched theta, 0.000e+00), is full rank at every theta, n and depth measured, has
    g_ii = 0.2500 exactly and is I/4 at depth 1; QNG has nothing to correct.
    `s13/results/geo_metric_arm1.json`.
C7. Width sweep of the deployed CVaR free energy on the deployed spectrum, n = 4..13, exact
    parameter-shift gradients, theta ~ N(0, 0.6^2): fitted log2 Var per qubit -0.649 (linear
    cost), -0.252 (alpha 0.25), -0.047 (alpha 0.10), -0.311 and -0.243 (deployed form). The
    CVaR non-linearity shrinks the gradient at n = 7 (0.52x, 0.093x) and flattens its decay so
    the ratio to the linear cost crosses 1 near n = 8 and reaches 3.36x / 2.68x at n = 13.
    Depth saturates at n = 7 from L = 4 (7.4e-3 to 8.2e-3, flat to sampling error). Scope:
    this shallow real-amplitude ansatz, these widths, this spectrum; not a 2-design result.
    `s25/results/q_plateau.json`.
C8. The dynamical Lie algebra of the deployed ansatz is the FULL so(2^n) from depth 2 at
    n = 4, 5, 7, 8, 10, 11 (8128 = so(128) at the deployed n = 7; depth 1 is abelian, dim n;
    n = 6 and n = 9 need depth 4, with 510 / 1023 and 32766 / 65535 at depths 2 / 3). The
    closure lies in the odd-Y (real) subalgebra at every cell. The Tang minimal pools V and G
    generate so(2^(n-1)+1) (36, 136, 528, 2080, 8256, 32896 at n = 4..9); the 2-local odd-Y
    pool generates so(2^n); ADAPT-selected sets at alpha = 1 stay abelian and at alpha = 0.25
    reach 1025 of 8128 at P = 21 with pool L2. Method: exact closure on Pauli strings as a set
    (symplectic anticommutation, bitwise product), agreeing with the numeric rank (dense
    -iP/2, nested commutators, SVD rank at relative tolerance 1e-10) on 12 of 12 cells at
    n = 4, 5. The 8128 and the 1025 are two generator sets (the fixed ansatz's conjugated
    RY generators; the strings ADAPT selected) and are never compared as one quantity.
    `s26/results/q_dla.json`, `s26/results/q_dla_smoke_n456.json`; pre-registered
    prediction of a proper subalgebra at n = 7 falsified (`s26/PREREG_A2.md` addendum).
C9. PENDING (A4): the grown circuits' variance beside the fixed ansatz at matched P.
    `s26/results/q_var.json` (running).
C10. The optimiser trains and the readout cannot tell: beats best-of-200 untrained draws at
    every T, closes 78 to 89% of the free-energy gap, sits 0.902 nats and 45% of mass from its
    Gibbs optimum at the deployed T, and the endpoint difference is 0.24x MDE.
    `s25/results/q_gibbs.json`, `s25/results/q_alpha.json`, `s26/results/q_mde_reference.json`.
C11. Set-equality theorem: the CVaR tail's support is a subset of a prefix of the energy
    order; 2,592 cells, 0 violations; equality on all 972 full-support cells.
    `s25/results/q_verify.json`.
C12. New (S26): the deployed Gibbs target is a product state (E affine in the register
    index); KL(Gibbs || product) 7e-17 on the ideal ladder, 5.4e-5 on a real target; a
    7-parameter RY layer reaches it, the 21-parameter circuit does not.
    `s26/results/probe/1A13.json`, `s26/q_tests.py`.

## 2. Figure list (each tied to a path)

    F1  support size vs sequence separation, CA-CA (exactly s-1) and all-atom (s+0.04);
        `s13/results/qarch_locality_geom.json`, `s13/results/qarch_locality_amber.json`
    F2  Pauli-weight spectra V_w: Legacy, raw AMBER, conditioned AMBER, with the Binomial(m,1/2)
        delta reference; `s13/results/walsh_xval.json`, `s13/results/geo_pauli.json`,
        `s13/results/geo_pauli_v1_rawonly.json`
    F3  measured vs predicted gradient variance, 78 points on the diagonal;
        `s13/results/geo_pauli.json`
    F4  the kernel v(w; n, depth): non-monotone, peak at w = 2, flat at depth >= 3;
        `s13/results/geo_kernel.json`
    F5  log2 Var[dF/dtheta_0] vs n for the five cells, fitted slopes; `s25/results/q_plateau.json`
    F6  Var[grad CVaR_alpha] / Var[grad mean] vs n (the ratio table);
        `s25/results/q_plateau.json` (`nonlinearity_variance_ratio`)
    F7  depth sweep at n = 7; `s25/results/q_plateau.json` (`depth_sweep_n7`)
    F8  dim(DLA) vs n for L = 1..6 against so(2^n), su(2^n), and the pools;
        `s26/results/q_dla.json` (PENDING)
    F9  grown vs fixed slopes, side by side, with actual P per row;
        `s26/results/q_var.json` (PENDING)
    F10 the Gibbs ladder: F, KL, TV, entropy from random theta through 5/15/50 Adam steps to
        the optimum, three temperatures; `s25/results/q_gibbs.json`
    F11 the product-state structure of the deployed Gibbs target: KL(Gibbs || product) by
        energy variant; `s26/results/probe/1A13.json` (one target; full table after A3)
    T1  the set-equality cells (n_cells, exact zeros, violations); `s25/results/q_verify.json`
    T2  the CVaR-vs-linear ratio table (as F6, numbers)
    T3  the Gibbs-ladder table at T = 0.3 (as F10, numbers)
    T4  the metric table by (n, depth): rank, trace, g_ii, off-diagonal correlation, cond;
        `s13/results/geo_metric_arm1.json`

## 3. Venue-honest statement: what would have to be added

For Quantum / PRX Quantum:
- Hardware or a noise model. Everything here is noiseless exact simulation (statevector, and
  exact MPS in the generation lane); no shot noise enters any number. The full-distribution
  readout the deployed gradient needs (128 outcomes) is a distribution-reconstruction cost,
  not a Pauli-expectation cost; a shot-cost analysis in the sense of Qiu et al. 2026
  (estimability) is absent.
- Width beyond n = 13 (statevector limit) with an MPS extension, and a 2-design control
  ansatz, so that "no exponential plateau at n <= 13" can be placed against the regime the
  theorems describe. The S13 exponential-vs-polynomial fits cannot discriminate at 6 to 14
  qubits (both R^2 >= 0.93); the depth-8 decay base 0.504 is the only 2-design-like number.
- A mechanism, or a theorem, for the CVaR flattening (C7); it is currently a reading.
- The DLA (C8) and its relation to the measured slopes through Larocca 2022 / Ragone 2024.
- A treatment of raw force-field conditioning as a modelling decision with its own
  spectrum (C2): the capped and soft potentials are modified observables and are labelled so.
- Multiple seeds and error bars on every variance (S25 used 250/200/120/80 draws; relative
  SE 9 to 16%).
For a QIP or TQC workshop contribution: C1 to C6 and C10 to C12 as they stand, with the
scope conditions in the text, would be enough for a contributed talk; C7 to C9 as
supporting material.

## 4. Two claims the paper must NOT make

1. Any quantum advantage. The register is 7 qubits, every state is simulated exactly, the
   selection is classical by theorem (C11), and Cerezo et al. 2025 place provably BP-free
   landscapes in the classically simulable class; the product-state target (C12) is the
   extreme case.
2. A barren plateau from a small gradient. A barren plateau is a statement about the decay
   of the gradient variance with width under random initialisation; a small gradient at one
   width is not one, and 6 to 14 qubits cannot distinguish exponential from polynomial decay
   (S13 section 6b). Where the paper says "no plateau" it says "at n <= 13, at depth 3, on
   this ansatz".

## 5. Related work (what the S13, S20, S22 literature passes found; see s13/lit_FINDINGS.md,
##    s14/lit_FINDINGS.md, s15/LITERATURE.md, s16/lit_FINDINGS.md, s21/agentD_FINDINGS.md 3.1-3.3)

- Larocca, Czelusta, Cerezo et al., "Diagnosing barren plateaus with tools from quantum
  optimal control", Quantum 6, 824 (2022): the DLA and the 1/dim(g) variance scaling.
- Ragone et al., "A Lie algebraic theory of barren plateaus for deep parameterized quantum
  circuits", Nat. Commun. 15, 7172 (2024); Fontana et al., Nat. Commun. 15, 7171 (2024).
- Larocca et al., "Barren plateaus in variational quantum computing", Nat. Rev. Phys. (2025),
  arXiv:2405.00781: the taxonomy used in s13/lit_FINDINGS.md section 2.
- Cerezo, Sone, Volkoff, Cincio, Coles, "Cost function dependent barren plateaus in shallow
  parametrized quantum circuits", Nat. Commun. 12, 1791 (2021): the locality mechanism whose
  assumptions this ansatz does not satisfy (C5).
- Cerezo et al., "Does provable absence of barren plateaus imply classical simulability? Or,
  why we need to rethink variational quantum computing", Nat. Commun. (2025), arXiv:2312.09121.
- McClean et al. 2018 (2-design plateaus); Holmes et al. 2022 (expressibility); Anschuetz and
  Kiani 2022 (traps without plateaus); Ortiz Marrero et al. 2021; Wang et al. 2021 (noise).
- Qiu, Lumbreras, Li, Rebentrost, "Quantum tilted loss in variational optimization",
  arXiv:2605.02850 (2026): CVaR as a tilted loss; does not remove plateaus; estimability.
- Barkoutsos et al., "Improving variational quantum optimization using CVaR", Quantum 4, 256
  (2020): the CVaR-VQE objective.
- Grimsley, Economou, Barnes, Mayhall, "An adaptive variational algorithm for exact molecular
  simulations on a quantum computer", Nat. Commun. 10, 3007 (2019); Tang et al.,
  "qubit-ADAPT-VQE", PRX Quantum 2, 020310 (2021): the pools and the growth rule used in A1.
- Wiersema, Kökcü, Kemper, Bakalov, "Classification of dynamical Lie algebras for
  translation-invariant 2-local spin systems in one dimension", arXiv:2309.05690: the
  reference for reading the A2 closures.
- QuPepFold, Uttarkar et al., PLoS One 21(2): e0342012 (2026): CVaR-VQE for peptide folding
  is a known combination; no structural metric reported there.
- Cumbo et al., "Logarithmic-scale variational quantum eigensolver for off-lattice protein
  structure prediction in continuous torsional angle space", arXiv:2609.02113 (2026):
  independent report that energy ranking fails on sampled landscapes.
- Zhang et al., arXiv:2510.06413 (2025): per-residue binned torsions, binary index encoding,
  IBM 127-qubit, 75 fragments, 4.89 A mean; the representational precedent.
