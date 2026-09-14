# PROPOSAL A: REPLACE THE FIXED ANSATZ BY qubit-ADAPT-VQE. VERDICT: REPLACE.

Lane Q, Sprint 26. Every number below is sourced in the notes block at the end.

## 1. What the proposal says

Replace the deployed selector's fixed 3-layer RY/CNOT circuit (7 qubits, 21 parameters) by
an ansatz that qubit-ADAPT-VQE grows one operator at a time on the same CVaR free energy,
choosing at each step the operator with the largest gradient. The claim behind it: a circuit
shaped by the problem should represent the trained distribution better than a fixed one, and
a better-represented distribution should give a better structure.

## 2. What the repository already knew

- The Hamiltonian the selector optimises is the same standardised rank ladder on every target
  (S25 L17: within 1.18% of its range). There are two trained states in the deployment, one
  per `(alpha, T)` cell of `VQE_LFO`.
- The selector's choice of candidates is classical by theorem (the CVaR tail is a prefix of
  the score order; 2,592 adversarial cells, 0 violations).
- The readout cannot tell the trained state from its own analytic optimum, which it misses
  by 0.902 nats and 45% of its mass: the endpoint gap is 0.24x its MDE.
- Deeper circuits and larger bond dimension ordered nothing, five sprints in a row.
- No exponential gradient plateau at 4 to 13 qubits; the algebra of the ansatz was not
  measured.

## 3. What this sprint measured

A2, the algebra (`s26/results/q_dla.json`, ledger L27). The fixed circuit's dynamical Lie
algebra is the whole real algebra so(2^n) from depth 2 at n = 7 (dimension 8128); depth 1 is
abelian. Nothing in the algebra protects the circuit from a plateau at scale; the "no plateau"
of S25 is a depth-3 statement. The Tang minimal pools generate a proper subalgebra
(dimension 2080 at n = 7); the 2-local pool generates all of so(128).

A4, the gradients of grown circuits (`s26/results/q_var.json`, ledger L35). On the deployed
Hamiltonian at alpha = 1, ADAPT grows no entangling operator: the grown circuit is a product
circuit whose gradients do not decay with width because there is nothing to train. At
alpha = 0.25 the grown 2-local circuit decays like the fixed one (-0.302 against -0.243 log2
per qubit).

The product-state fact (`s26/results/probe/1A13.json`, `s26/results/a1/`, ledger L68). The
deployed energy is affine in the register index, so the Gibbs state the alpha = 1 objective
asks for is a product state on every one of the 126 targets (KL to the product of its
marginals at most 7.9e-4 nats). Seven RY angles reach it (KL 0.0002); the 21-parameter fixed
circuit stops 0.90 nats short; ADAPT with either pool selects no entangling operator on any
of the 78 alpha = 1 targets under L-BFGS.

A1, the endpoint (`s26/results/a1_stats.json`, ledger L68). Paired on 126 targets, built
chain through the production projection, fixed ansatz against ADAPT at matched 21 parameters:

    ADAPT, 2-local pool  - fixed   -0.0138 A   SE 0.0210   MDE 0.0588   0.23x   fold CI [-0.071, +0.044]   58W/68L
    ADAPT, Tang pool     - fixed   -0.0222 A   SE 0.0217   MDE 0.0608   0.36x   fold CI [-0.085, +0.042]   61W/65L

Both are inside +-0.5x their own MDE with fold CIs spanning zero: the registered "ADAPT is
null" falsifier fired. All twelve ADAPT arms (two pools, two optimisers, 7/14/21 parameters)
are negative by 0.013 to 0.025 A and none clears its MDE. On the selection readout the twelve
arms are -0.039 to -0.048 A at 0.42x to 0.47x MDE with 45 to 59 exact ties: the same shape as
S25's headline, a consistent direction whose magnitude the instrument cannot measure. Power:
the comparison resolves 0.06 A; underpowered below that; the whole quantum synthesis is worth
+0.013 A against the classical arm, so the resolution is four times the component's footprint.

## 4. Verdict: REPLACE

Not KEEP WITH EDITS, because the proposal's mechanism is absent, not weak: the object the
circuit is asked to represent is a product state, an adaptive ansatz confirms it by growing
nothing, and reaching that state exactly moves the answer by a third of the MDE. A better
ansatz has nothing to be better at.

What replaces it: Proposal B's replacement form (`s26/PROPOSAL_B_REPLACEMENT.md`), the
trainability paper, now with A2, A4 and the product-state result as its Sprint 26 additions.

## 5. What the presenter says (under two minutes)

"We tested the idea of letting the circuit grow itself. qubit-ADAPT-VQE starts from one layer
of single-qubit rotations and adds, one at a time, whichever operator would lower the
objective fastest. We ran it on all 126 development peptides with the same objective, the
same seed and the same readout as the deployed circuit, and compared the emitted structures
pairwise.

The answer is no change: the grown circuit's structures are 0.014 to 0.022 angstroms closer
to the experimental ones on average, a third of what this comparison can resolve, with the
confidence interval across the five folds straddling zero. Every variant we tried points the
same way and none of them clears the bar.

We also found out why. The energy the circuit minimises is the same rank ladder on every
peptide, and its optimum at the deployed setting is a product state: seven independent
rotations represent it exactly. When ADAPT is offered entangling operators it declines them
on all 78 targets at that setting; there is nothing to entangle. The deployed 21-parameter
circuit actually stops 0.9 nats short of that optimum, and reaching it exactly changes the
structure by a third of the resolution. We also computed the circuit's Lie algebra: it is the
full real algebra from depth 2, so its trainability at 4 to 13 qubits is a property of its
shallowness, not of any special structure.

So Proposal A is replaced. The adaptive circuit is a correct, well-tested tool, and it
confirmed the diagnosis instead of curing it: the selector's Hamiltonian carries almost no
per-target information, and the readout cannot see the difference between two well-trained
states. What we publish instead is the trainability work, with these measurements added."

## 6. Notes block (every number and its artefact)

    "126 development peptides, paired"          s26/results/a1_stats.json, n = 126, folds pinned
    "0.014 to 0.022 angstroms closer"           built chain, adaptL2/adaptV P21 - fixed_zrank_it50: -0.0138 / -0.0222 A
    "a third of what it can resolve"            0.23x / 0.36x of MDE 0.0588 / 0.0608 A (ST.compare, MDE = 2.8016 SE)
    "interval straddling zero"                  fold-clustered CI [-0.071, +0.044] / [-0.085, +0.042]
    "every variant points the same way"         12 of 12 ADAPT arms negative, -0.013 to -0.025 A, 0.21x to 0.40x MDE, s26/logs/a1_stats.log
    "same rank ladder on every peptide"         S25 L17: worst deviation 1.18% of range, s25/results/q_gibbs.json
    "optimum is a product state"                KL(Gibbs || product) max 7.9e-4 nats over 126 targets, s26/results/a1/*.json product_diagnostics
    "seven rotations represent it exactly"      adaptV_lbfgs_zrank_P7 KL to Gibbs 0.0002 (max 0.0009) on the 78 alpha = 1 targets
    "declines them on all 78 targets"           L-BFGS growth stops with no operator selected, 78 of 78, s26/results/a1/*.json adapt.*.stopped
    "0.9 nats short"                            fixed_zrank_it50 KL to Gibbs 0.9027 on the 78 alpha = 1 targets; S25: 0.902, s25/results/q_gibbs.json
    "full real algebra from depth 2"            dim(DLA) = 8128 = dim so(128) at n = 7, L = 2, s26/results/q_dla.json, ledger L27
    "shallowness, not structure"                S25 slopes -0.65 to -0.05 at depth 3 (s25/results/q_plateau.json) vs the depth-8 decay base 0.504 (s13/results/geo_kernel.json)
    "readout cannot see the difference"         S25: circuit vs exact Boltzmann -0.0302 A, 0.24x MDE, s26/results/q_mde_reference.json; A1: gibbs_T - fixed +0.0088 A, 0.11x
    "well-tested tool"                          s26/q_tests.py 17/17; simulator vs core.quantum and pennylane to 1e-12; production arm reproduced bit-for-bit on 126/126 (cache_check)
    "grown circuits' gradients do not decay"    s26/results/q_var.json slopes -0.079 / +0.006 / +0.035 / -0.008 at alpha = 1; ledger L35
