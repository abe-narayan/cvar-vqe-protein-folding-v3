# IDEA (lane Q) -- A TARGET-DEPENDENT HAMILTONIAN FOR THE CVaR-VQE SELECTOR (S25 L17)

## Hypothesis

The deployed Hamiltonian `H = diag(zrank(sorted scores))` is the same rank ladder on every
target (S25 L17: 1.18% of range, 0 of 8 identical). A Hamiltonian whose spectrum carries the
target's own score GAPS, while every candidate keeps its rank, gives the variational state
target-specific structure to represent. If the readout responds to that structure, the
selector stops being two trained states and becomes 126.

## Candidates, all order-preserving (the CVaR tail stays a prefix of the score order; S25 theorem)

1. `zraw`: the moment z-score of the raw Bayes-risk score over the 128. Cheapest; the gaps
   survive intact. Risk: the Gibbs state at T = 0.3 is much sharper (1.84 bits on 1A13 against
   4.91 for zrank, `s26/results/probe/1A13.json`).
2. `asinh` (S21 L33's `Nt`): asinh of a robust z-score, re-standardised; compresses outliers
   without touching the order. 4.05 bits on 1A13.
3. `soft` (S13 sections 7 and 8, `s13/geo_common.py` `soft_compress`): identity below the
   median, log compression above; re-standardised. 1.72 bits on 1A13.
4. Any of 1 to 3 with a per-target temperature chosen native-free so that the Gibbs entropy
   equals the deployed one (4.91 bits): keeps the readout's spread fixed and varies only the
   gaps. This is the version that tests the L17 mechanism cleanly.
5. NOT order-preserving, listed for completeness: the posterior-MEAN risk (squared loss)
   instead of the median (L1). It changes which candidate is where; S21 L14 measured the
   squared-vs-Bayes functional at +0.114 A SE 0.069 on the latent and +0.050 on pool windows.
   Outside the "preserve selection semantics" brief and not run here.

## Why the record does not already close it

- S25 L17 identified the constant spectrum and stated the consequence (two trained states)
  but ran no target-dependent variant.
- S25 L16 measured that raw moment standardisation of AMBER is not monotone in float64; that
  is a property of 1e28 outliers in raw force-field energies and does not apply to a
  distogram score of order 1 A. Monotonicity is asserted at build time.
- S23 L8 closed the probability-WEIGHTED readout for the deployed Hamiltonian at two
  temperatures; it did not vary the Hamiltonian.
- S22 L16 and S25 L15 close the SET (membership) channel by theorem; the weights channel
  with a target-dependent spectrum was never opened.
- New this sprint (`s26/agentQ_FINDINGS.md` section 0): the deployed Gibbs state is a product
  state, which is why the entangled ansatz has nothing to represent. `zraw` on 1A13 moves the
  Gibbs state 0.042 nats from a product state; `asinh` 0.025; `soft` 0.038. Small, but not
  zero: the first non-product target the selector has ever been given.

## Exact falsifier (PREREG_A3)

On the built chain, paired at n = 126, fold-clustered: `fixed_zraw_Tmatch - fixed_zrank`
negative past its MDE with the fold CI excluding zero, 5/5 folds, replicated on seed 1 and in
reverse order, is "helps". Within +-0.5x MDE is "null". Positive past MDE is "harmful", and is
a result about the entropy mechanism. Property falsifiers: distinct trained states under zraw
at least 100 of 126 (else the variant did not make the states target-dependent).

## Expected effect against the computed MDE

Reference built-chain MDE 0.0958 A (weighted vs uniform over the same 128,
`s26/results/q_mde_reference.json`). Expected: `zraw` alone POSITIVE (worse) by +0.03 to +0.10
on the selection basis via S25's entropy curve; `zraw_Tmatch` within 0.5x MDE (null). The
honest prior is that the readout responds to the entropy of the weights and to nothing else
(S25 section 6.2: rho -0.74 with entropy; alpha adds 3.2% of the variance), so a
target-dependent spectrum with matched entropy is expected to be invisible at the endpoint.

## Memory and agent-hours

Same harness as A1: peak RSS under 0.4 GB (`s26/jobs_done/q_probe_1A13.json`), about 95 s per
target with 21 projections, 3.3 h wall for 126 targets on one process. Agent-hours 4 to 6
including the write-up. Runs after "PHASE 0 SIGNED OFF" (endpoint half); the property half
(distinct states, sequences, product diagnostics) reads no native.
