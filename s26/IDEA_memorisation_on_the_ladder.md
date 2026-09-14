# IDEA_memorisation_on_the_ladder -- PLACE THE OWN-NATIVE MODEL ON THE S24 PRIOR LADDER: THE FIRST REAL TRAINED OPERATOR WITH A MEASURED gam_eff, cos AND ENDPOINT (lane W, Sprint 26; found on the way to L44)

## Hypothesis

The project's only steep lever is the prior's derivative: interpolating the shipped posterior
toward the native moves the endpoint -2.15 A per unit gamma at the origin (S24 L13), but that
ladder's gamma is redeemable only along the native's own direction (cos = 1 by construction) and
S25 L12 caps a real operator travelling 25% at cos 0.5 at +0.024 A. Every real operator measured
since (51 re-readings of the posterior, S25 L12; the C2 ladder's rungs) sits at gam_eff within a
few hundredths. L44 Part C supplies a real, TRAINED posterior that moved a long way toward the
native and bought 0.70 A on the built chain: the four pinned models that had the target's native
among their labels. Nobody has measured where on the ladder that operator sits. The hypothesis
is that its (gam_eff, cos) pair, entered into the ladder's transfer function, predicts its measured
-0.70 A: i.e. the ladder is a valid transfer function for a trained operator, not only for an
interpolation. If it does not (for example gam_eff 0.3 at cos 0.6 predicting -0.3 A against a
measured -0.70), the ladder under-prices what a better trained prior buys, and the "0.0225 reaches
3.0 A" arithmetic (S24 L13) is the wrong currency for Proposal C.

## Why the record does not already close it

- S24 L13 built the ladder from interpolations (MASS, TILT, DIRAC, MASSFIXW constructions), all
  at cos = 1. S25 L12 measured 51 re-readings at small gam_eff and found corr(gam_eff, endpoint) =
  +0.054, i.e. the ladder does not predict SMALL moves off the native direction. No arm in the
  record has a LARGE gam_eff at cos < 1 with an endpoint: the interpolations have cos 1, the
  re-readings have gam_eff ~ 0. The own-native models are the missing quadrant.
- Lane P's C2 ladder measures gam_eff and cos for every rung (`s26/p_ladder.py progress`), but
  every rung is a legitimate leave-fold-out model with gam_eff expected within 0.05; none is
  designed to be far from the origin.
- L44 stored the leaked models' emissions and their posterior differences from the clean one
  (max/mean |dprob|, mean |dE[d]|: 0.70 to 0.79 A per pair on 1CEK), not gam_eff or cos; the
  posteriors themselves are recomputed in seconds from the five pinned models.

## Exact falsifier

ORACLE DIAGNOSTIC throughout (the operator saw the native). For each target and each of its four
leaked models: gam_eff and cos in probability space and in location space (the S24 MASS
construction and the S25 loc currency, both already coded in `s26/p_ladder.py progress`), the
posterior's MAE against the native, and the measured built-chain delta from
`s26/results/w_selfcopy_endpoint.json :: C/rows`. The ladder's prediction for the delta is
-2.1496 x gam_eff x g(cos) where g is the S25 L12 direction discount; the registered test is
whether the 504 (target, model) points and their mean fall on the ladder: the falsifier is that
the ladder's prediction and the measured -0.698 A disagree by more than the measured effect's
own MDE (0.246 A on the built chain). Falsified either way informs: agreement validates the
ladder as a transfer function for trained priors; a ladder prediction much SMALLER than -0.70
(the registered expectation, because a trained model's move is off the native's direction, cos
0.5 to 0.7) says the ladder under-prices trained priors and the Proposal C arithmetic should be
restated in terms of the measured operator.

## Expected effect against the computed MDE

Expected gam_eff (probability space) 0.15 to 0.35 at cos 0.4 to 0.7 for the leaked models
(they change half the top-75); the ladder at cos 1 would price gam_eff 0.25 at -0.54 A, and
with the S25 direction discount at about -0.1 to -0.2 A, against the measured -0.70. So the
registered expectation is DISAGREEMENT by 2 to 3 MDE: the ladder under-prices a trained
operator. If instead gam_eff comes out near 0.7 at cos ~ 0.9 (near-recall of the native), the
ladder agrees and the memorisation is a near-interpolation. The measurement decides between
the two readings; both are informative.

## Memory and agent-hours

Five pinned MLPs and `esm_small.npz` (the envelope run loaded them at 0.318 GB peak); 504
posteriors at ~0.1 s, no projection (the endpoints exist): under 0.4 GB, 5 minutes, one
governed job. 1.5 agent-hours. Gated (reads native distances for gam_eff and MAE); needs its
PREREG first; ORACLE labelled on every line.

## Information value

High for Proposal C and for the report's central sentence about the prior: it either validates
the "-2.15 A per unit gamma" transfer function on the one trained operator that moved far, or
shows that the currency under-prices trained priors, which changes what "gamma = 0.0225 reaches
3.0 A" means for a C2-style rung. Plausibility that the ladder disagrees: 0.6.
