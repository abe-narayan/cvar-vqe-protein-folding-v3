# PROPOSAL B, REPLACEMENT FORM: PUBLISH THE TRAINABILITY RESULT (lane Q, Sprint 26)

This is the paper-outline form the campaign prompt requires if Proposal B is replaced. Lane P's
`s26/PROPOSAL_B.md` carries the B1 to B3 evidence and points here. It is built from
`s26/TRAINABILITY_PAPER_OUTLINE.md` with this sprint's A2 and A4, both measured
(`s26/results/q_dla.json`, `s26/results/q_var.json`; ledger L27, L35). Every figure is tied to
an artefact path. Nothing here is a claim of advantage.

## 1. The replacement in one paragraph (what the presenter says)

We measured, exactly and with no free parameter, how a real molecular force field's structure
reaches the gradients of a variational quantum circuit. The chain runs: chain geometry fixes
which residues a distance term can depend on (an exact theorem); that fixes the energy's
Pauli-weight spectrum (an exact transform, once the raw force field's clash spike is
conditioned away); the spectrum times the circuit's own gradient kernel predicts the measured
gradient variance to 0.4%. At the depth we use, the kernel is flat in Pauli weight, so the
force field's locality is not what limits training; the circuit's width is. The circuit's
Lie algebra is the whole real algebra from depth 2, so nothing algebraic protects it from a
plateau at scale; what we observe at 4 to 13 qubits is the shallow regime. The state the
optimiser is asked to reach is a product state, which seven parameters represent exactly,
and an adaptive ansatz grown on that target grows no entanglement: its gradients do not decay
because there is nothing to train. The paper is about what a spectrum does and does not
predict. It does not claim a quantum advantage, does not call a small gradient a barren
plateau, and does not call a large gradient trainability.

## 2. Claims, each with scope and artefact

    C1  exact locality theorem in torsion space              s13/results/qarch_locality_geom.json,
        (support = the j-i-1 residues between; all-atom     qarch_locality_amber.json,
        one residue wider; no 2-local Ising form)           walsh_predict.json
    C2  raw AMBER's Pauli spectrum is a delta spike          s13/results/walsh_xval.json,
        (top-10 configurations carry 99.6% of the variance; geo_pauli_v1_rawonly.json, geo_anova.json
        mean weight = m/2 exactly; winsorisation fails)
    C3  conditioned identically, AMBER is the higher-weight  s13/results/geo_pauli.json,
        observable: 3.015 vs 2.236, 79/79 cells; both       walsh_amber.json, walsh_legacy_k4_m18_6t.json
        models are a steric term plus rounding error
    C4  spectrum x kernel predicts gradient variance,        s13/results/geo_pauli.json
        median measured/predicted 0.996 and 0.998, 78 cells
    C5  the kernel is flat in Pauli weight at depth >= 3;    s13/results/geo_kernel.json
        Cerezo 2021's mechanism is inactive here
    C6  the metric has no Hamiltonian in it; full rank;      s13/results/geo_metric_arm1.json
        g_ii = 0.25; I/4 at depth 1
    C7  width and depth sweeps of the deployed CVaR free     s25/results/q_plateau.json
        energy, n = 4..13: slopes -0.649 / -0.252 / -0.047 /
        -0.311 / -0.243 log2 per qubit; CVaR/linear ratio
        crosses 1 near n = 8; depth saturates from L = 4
    C8  MEASURED THIS SPRINT (A2): the DLA is the full       s26/results/q_dla.json,
        so(2^n) from depth 2 at n = 7 (8128); depth 1        s26/results/q_dla_smoke_n456.json
        abelian; n = 6, 9 need depth 4; Tang pools generate
        so(2^(n-1)+1); 2-local odd-Y pool generates so(2^n);
        ADAPT sets abelian at alpha = 1, 1025 of 8128 at
        alpha = 0.25; symbolic == SVD rank (rtol 1e-10) 12/12
    C9  MEASURED THIS SPRINT (A4): ADAPT-grown circuits at   s26/results/q_var.json,
        matched P = 3n are product circuits at alpha = 1     s26/PREREG_A4.md addendum
        (no decay: -0.079 / +0.006 at T = 0, +0.035 / -0.008
        at T = 0.3; variance 2 to 370x the fixed ansatz's);
        the alpha = 0.25, T = 0.3 L2-grown circuit decays at
        -0.302 vs fixed -0.243 (equal within error); T = 0
        grown cells degenerate (collapse); S25's n = 7 rows
        reproduced at relative deviation 0.0
    C10 the optimiser trains and the readout cannot tell:    s25/results/q_gibbs.json,
        78 to 89% of the gap closed; 0.902 nats and 45% of   s25/results/q_alpha.json,
        mass from the optimum; endpoint 0.24x MDE            s26/results/q_mde_reference.json
    C11 set-equality theorem: 2,592 cells, 0 violations      s25/results/q_verify.json
    C13 NEW (A3): a target-dependent, order-preserving          s26/results/a3_stats.json, ledger L125
        Hamiltonian makes the trained states target-dependent
        (124/126 distinct) and the emitted structure does not
        move (+0.003 A, 0.04x MDE) once the entropy is matched;
        unmatched, sharper states are worse (+0.09 to +0.11 A,
        0.7 to 0.8x MDE, Type-M zone)
    C12 NEW: the deployed Gibbs target is a product state;   s26/results/probe/1A13.json,
        KL 7e-17 ideal, 5.4e-5 real; a 7-parameter RY layer  s26/q_tests.py
        reaches it, the 21-parameter circuit does not

## 3. Figure list (slides 6 and 8 use F8 and F9)

    F1  support size vs separation, CA-CA and all-atom           s13/results/qarch_locality_geom.json, qarch_locality_amber.json
    F2  Pauli-weight spectra with the Binomial(m,1/2) reference   s13/results/walsh_xval.json, geo_pauli.json, geo_pauli_v1_rawonly.json
    F3  measured vs predicted gradient variance (78 points)       s13/results/geo_pauli.json
    F4  the kernel v(w; n, depth)                                 s13/results/geo_kernel.json
    F5  log2 Var vs n, five cells, fitted slopes                  s25/results/q_plateau.json
    F6  Var[grad CVaR]/Var[grad mean] vs n                        s25/results/q_plateau.json
    F7  depth sweep at n = 7                                      s25/results/q_plateau.json
    F8  dim(DLA) vs n: fixed ansatz by depth, pools, ADAPT sets,  s26/results/q_dla.json  ->  s26/figures/a2_dla_dimension.png
        against so(2^n) and su(2^n)
    F9  fitted log2 Var per qubit, fixed vs grown, per cell      s26/results/q_var.json  ->  s26/figures/a4_variance_slopes.png
        (matched-P rows), with the deployed cells' Var-vs-n curves
    F10 the Gibbs ladder (F, KL, TV, entropy along training)      s25/results/q_gibbs.json
    F11 KL(Gibbs || product) by energy variant                    s26/results/probe/1A13.json (one target; full table after A3)
    F12 dim of the Lie closure of the ADAPT-grown set at every growth step on the 126 A1
        records (median and min..max per cell): the algebra of the SET grows while the
        angles are inert at alpha = 1; s26/results/q_dla_a1.json -> s26/figures/a2_dla_grown_ladder.png (PENDING)
    F9b F9 now carries 95% bootstrap CIs on every slope (error bars) and on grown - fixed per
        cell; s26/results/q_var_boot.json (grown L2 - fixed at alpha=0.25, T=0.3: -0.056 [-0.182, +0.087])
    T1  set-equality cells                                        s25/results/q_verify.json
    T2  the CVaR-vs-linear ratio table                            s25/results/q_plateau.json
    T3  the Gibbs-ladder table at T = 0.3                         s25/results/q_gibbs.json
    T4  the metric by (n, depth)                                  s13/results/geo_metric_arm1.json

## 4. Venue-honest statement: what a submission would still need

Quantum / PRX Quantum:
- Hardware runs or a noise model. Every number is noiseless exact simulation (statevector;
  exact MPS in the generation lane). The deployed gradient needs the full 128-outcome
  distribution, a reconstruction cost; the estimability analysis in the sense of Qiu et al.
  2026 (arXiv:2605.02850) is absent.
- Width beyond n = 13 with an MPS extension and a 2-design control ansatz, so that "no
  exponential plateau at n <= 13" is placed against the regime the theorems describe. At
  6 to 14 qubits exponential and polynomial fits are indistinguishable (S13 section 6b, both
  R^2 >= 0.93); the depth-8 decay base 0.504 is the only 2-design-like number.
- A mechanism, or a theorem, for the CVaR flattening (C7); today it is a reading.
- The link from C8's maximal algebra to C7's slopes stated as a depth statement, with the
  parameter count against dim(g) at every n (P = 3n against 2^(n-1)(2^n - 1)).
- A treatment of raw force-field conditioning as a modelling decision (C2): the capped and
  soft potentials are modified observables and must be named so throughout.
- Seeds and error bars on every variance (relative SE 9 to 16% at 250 to 80 draws).
QIP / TQC workshop: C1 to C6 and C10 to C12 as they stand, with scope conditions in the
text, are a contributed talk; C7 to C9 as supporting material.

## 5. Two claims that must not be made

1. Any quantum advantage: 7 qubits, exact simulation, selection classical by theorem (C11),
   the optimisation target a product state (C12), and Cerezo et al. 2025 on BP-free
   landscapes being classically simulable.
2. A barren plateau from a small gradient: a plateau is the decay of the variance with width
   under random initialisation; "no plateau" is said only as "at n <= 13, depth 3, this
   ansatz", and C8 says the algebra offers no protection at scale.

## 6. Related work

As in `s26/TRAINABILITY_PAPER_OUTLINE.md` section 5: Larocca 2022 (Quantum 6, 824); Ragone
2024 and Fontana 2024 (Nat. Commun. 15); Larocca 2025 review (arXiv:2405.00781); Cerezo 2021
(Nat. Commun. 12, 1791); Cerezo 2025 (arXiv:2312.09121); McClean 2018; Holmes 2022; Anschuetz
and Kiani 2022; Qiu et al. 2026 (arXiv:2605.02850); Barkoutsos 2020 (Quantum 4, 256);
Grimsley 2019 (Nat. Commun. 10, 3007); Tang 2021 (PRX Quantum 2, 020310); Wiersema et al.
(arXiv:2309.05690); QuPepFold (PLoS One 21(2): e0342012); Cumbo et al. (arXiv:2609.02113);
Zhang et al. (arXiv:2510.06413).

## 7. Notes block for the presenter (every number's source)

    "no free parameter, 0.4%"        median measured/predicted 0.996, s13/results/geo_pauli.json
    "flat in Pauli weight"           v(n)/v(1) = 0.88 to 1.06 at depth >= 3, s13/results/geo_kernel.json
    "whole real algebra from depth 2" dim 8128 = so(128) at n = 7, L = 2, s26/results/q_dla.json
    "4 to 13 qubits, shallow regime" slopes -0.649 to -0.047 log2 per qubit, s25/results/q_plateau.json;
                                     2-design rate base 0.504 at depth 8, s13/results/geo_kernel.json
    "a product state, seven params"  KL(Gibbs || product) 7.2e-17 (ideal) / 5.4e-5 (1A13),
                                     s26/q_tests.py, s26/results/probe/1A13.json
    "grown circuits do not decay"    slopes -0.079 / +0.006 / +0.035 / -0.008 at alpha = 1;
                                     -0.302 vs -0.243 at alpha = 0.25, s26/results/q_var.json
    "clash spike, 99.6%"             top-10 share 0.996 median, s13/results/walsh_xval.json
    "0.902 nats, 45% of mass"        s25/results/q_gibbs.json; endpoint 0.24x MDE, s26/results/q_mde_reference.json
