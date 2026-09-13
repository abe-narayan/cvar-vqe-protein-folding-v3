# IDEA (lane Q) -- THE DEPLOYED SELECTOR'S TARGET IS A PRODUCT STATE; SAY SO, AND PRICE THE
# 7-PARAMETER REPLACEMENT

## Hypothesis

At the deployed cell alpha = 1, T = 0.3 (folds 0, 3, 4; 78 of 126 targets) the objective is
`mean_p(E) - T H(p)`, whose minimiser is the Gibbs state `exp(-E/T)/Z`. Because `E` is affine
in the register index (standardised ranks; the register index IS the rank), `exp(-E/T)`
factorises over the 7 bits: the Gibbs state is a PRODUCT state. A single layer of RY rotations
(7 parameters, no entangler) represents it exactly; the deployed 21-parameter RY/CNOT circuit
sits 0.902 nats from it because its CNOT layers stand between a product input and a product
target. The proposal: (i) state the mechanism in the claims; (ii) measure whether the exact
product Gibbs state, or the 7-parameter RY layer trained by L-BFGS (about 10 gradient
evaluations against 2,100 circuit settings), changes the emitted structure.

## Evidence at pre-registration (all native-free)

- KL(Gibbs || product of its marginals) = 7.2e-17 nats on the ideal ladder at n = 7 and every
  T in {0.1, 0.3, 1.0} (`s26/q_tests.py` `test_gibbs_state_of_the_deployed_ladder_is_a_product_state`),
  5.4e-5 nats on 1A13's real E with its ties (`s26/results/probe/1A13.json`).
- The trained RY layer reaches KL 1e-4 to the Gibbs state (arm `adaptV_lbfgs_zrank_P7`, same
  file); the deployed circuit 0.903 nats at 50 steps, 0.884 at 750.
- ADAPT with either pool selects no entangling operator at alpha = 1 and stops at P = 8 or 10.
- S25 already has the selection-basis endpoint of the exact Gibbs state: `VQE_LFO -
  boltz_T0.3 = -0.0002 A`, 0.00x MDE (`s26/results/q_mde_reference.json`). The built-chain
  basis is not measured; the A1 arm `gibbs_T` measures it.

## Why the record does not already close it

S25 measured the divergence (0.902 nats) and the endpoint equivalence with the Boltzmann
weights on the selection basis; it did not identify WHY the circuit cannot reach its optimum,
and it reported the circuit as "broader than optimal" without the product-state reading.
S13 section 4d found depth 3 to 4 needed to fit a multimodal target; the deployed target is
not multimodal, it is a product, which is a different statement about what the entangler is
for. No sprint ran the 7-parameter product ansatz as an arm.

## Exact falsifier

A1 arm `gibbs_T - fixed_zrank_it50` on the built chain: within +-0.5x MDE is "the replacement
is free"; negative past its MDE with the fold CI excluding zero is "the deployed circuit is
costing accuracy by missing its own optimum"; positive past MDE is "the deployed circuit's
broadness helps", which would be a real finding about the readout's entropy preference.

## Expected effect against the computed MDE

Built-chain reference MDE 0.0958 A (`s26/results/q_mde_reference.json`); expected |effect|
below 0.03 A (the selection-basis analogue is -0.0002 A). The value of the idea is the claim,
not the Angstrom: "the quantum selector's optimisation target is a product state that 7
classical parameters represent exactly" is a sentence the presenter can say in ten seconds
and it replaces "the readout is insensitive to a 45% mass difference" with its cause.

## Memory and agent-hours

Zero extra compute: the arm is inside the A1 harness. Writing 1 agent-hour.
