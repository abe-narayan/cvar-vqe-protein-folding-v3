# IDEA_amber_prior_partner -- AMBER AS A DISTRIBUTION INSIDE THE PRIOR, NOT AS A RANKER (lane P, Sprint 26; entered per LEDGER L16b from PREREG_C3.md's H_C3a)

## Hypothesis

Per pair, weight the K = 500 pool's 17-bin distance histogram by exp(-beta * zrank(E_AMBER))
using the 63,000 cached ff14SB/GBn2 single points (`s24/cache_amber/*.npz`, one per target,
energies aligned with the shipped pool indices), and use that energy-weighted histogram as the
mixture partner of the shipped posterior: P = (1 - lam) P_shipped + lam H_AMBER(beta), lam and
beta chosen leave-fold-out on the built chain. The mixed posterior goes through the shipped
Bayes-risk functional unchanged (`s26/p_ladder.py` rung `mix` with H_pool replaced by H_AMBER).
The claim: AMBER's non-parallel selection direction (bias cosine 0.569 vs the incumbent, -0.215
against its own floor, S24 L10-A) becomes usable when it enters as a distribution the prior is
mixed with, rather than as a score, filter or partition.

## Why the record does not already close it

The functional lever is closed in all five forms (S24 L16, S25 L18): AMBER as a filter
(+0.26 to +0.55, S24 L10-A), as a partition (every partition worse), as an additive score
(nested CV picks w = -0.5; ORACLE w worth 0.015 A), as a per-target sign, and as an ORACLE
weight. In every one of those AMBER acts on the CANDIDATES' ranking. Here it acts on the PRIOR's
location and width per pair, and the shipped score then ranks candidates against that prior.
The unweighted version (the plain pool histogram) is rung `mix` of PREREG_C2 and is itself
foreshadowed null (S7-3 `gain 0`, S10-2, S19 L14); the energy weighting is the only new element,
and it is separable from the histogram by construction. The honest prior is a null: S24 L16's
rank-permuted AMBER control (-0.0010, 0.08x MDE) says AMBER's ordering carries no accuracy
through the score, and a histogram weighted by an ordering that carries no accuracy is the
unweighted histogram plus noise.

## Exact falsifier

PREREG_C2's, on the built chain: (lam*, beta*) chosen leave-fold-out must beat the shipped
posterior by more than its own MDE with the fold-clustered CI excluding zero and 5/5 folds,
AND must beat the UNWEIGHTED mixture (rung `mix`) at its own leave-fold-out lam by the same
standard, AND must beat a rank-permuted-AMBER mixture (marginal preserved, correspondence
destroyed, S24 L16's control). lam = 0 and beta = 0 are asserted bit-exact to the incumbent per
target. The per-target oracle over the (lam, beta) grid is reported beside `ST.best_of_k_within`
and its split-half transfer, never as a result.

## Expected effect vs computed MDE

MDE on the built chain 0.09-0.20 A (PREREG_C2 section 4). Expected effect: 0.00 +/- 0.03
(indistinguishable from `mix`); the design can refute a 0.1 A gain and cannot see a 0.03 A one.

## Memory and agent-hours

Reads `s24/cache_amber/<pdb>.npz` (1.5 MB total), the universes one at a time, and the s12
distogram cache; < 0.6 GB; 8 lam x 4 beta x 126 targets x one projection (2.9 s) = 3.2 h of
machine time if every cell is projected, 25 min if only the leave-fold-out winner and the three
controls are. 2 agent-hours. No OpenMM call: every energy is already on disk (S24 L16).
