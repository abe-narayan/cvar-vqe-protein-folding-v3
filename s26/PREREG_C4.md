# PREREG_C4 -- ROUTERS FOR THE PER-TARGET SCALE s* AND SET SIZE m WITH A FEATURE SET NO PREVIOUS ROUTER USED (lane P, Sprint 26)

Written 2026-09-13, before any router is fitted.

## What the record already knows

- Per-target m* is REAL and transfers split-half (-0.2394, SE 0.0353, 2.42x MDE, 5/5 folds, S22
  L4/L6; re-audited S24 L15-A). Five router constructions fail to predict it native-free: length
  quartile (+0.031, wrong-signed), Ridge over 16 arms (+0.070 [+0.005, +0.146], significantly
  HARMFUL), Random Forest (in-fold -0.276 -> held-out +0.079), Ridge + RF on candidate-set
  geometry (+0.014; RF +0.060 [+0.004, +0.117], significantly HARMFUL) (S22 L7), and rg_z
  (+0.0012, 2% of its MDE, S23 L7). The S22 L10 bound: a p = 5 linear router needs n ~ 281 for a
  0.24 A effect; a one-threshold router needs n ~ 180.
- Per-target s* is REAL, closed-form (s* = <c,t>/|c|^2), transfers 98.9% split-half (-0.3198,
  S23 L6/S24 L15-A) and is UNREACHABLE IN PRINCIPLE: s* = 1 - <c, ebar>/|c|^2 depends only on the
  invisible common-mode error (S23 L9), and a mismatched same-length native serves as well as the
  true one (S23 L6d). Eight native-free features are all inside |rho| <= 0.11 (S23 L3).
- Features already used and therefore NOT new: n; score mean/sd/skew/gaps/IQR; distogram sd
  mean/max and expected mean/sd; consensus mean/sd/median; BLOSUM sim mean/sd (top-75 and pool);
  org fraction; rg of the average, predicted rg, rg_gap, rg_z, rg_sd; pool_spread and the
  spread(m) ladder (S12 `agg_router.features`, S22 `routercv.ALL_FEATS`, S22 `mgeomrouter.GEOM_FEATS`,
  S23 `c2_rgcond`).

## Hypothesis

**H_C4.** A feature set no previous router used -- (1) the principal-axis spread of the top-75
cloud (eigenvalues of the superposed members' coordinate covariance: anisotropy ratios, the
fraction of variance on the first axis), (2) the entropy of the retrieval-score distribution
(softmax over the K = 500 BLOSUM sims at the pinned temperature 1, and the top-1/top-75 sim
gaps normalised by the pool sd), (3) the per-pair 17-bin entropy of the shipped posterior (mean,
max, sd, fraction of pairs with >= 2 modes above 0.02; the information-theoretic cousin of the
S12 `dg_sd_mean`, declared as such), (4) ESM-2 650M contact-head statistics (mean and entropy of
the contact map, long-range (|i-j| >= 6) contact mass, max contact probability, contact-map
anisotropy) -- predicts m* or s* native-free with held-out gain above its MDE.

## Exact falsifier

Nested leave-fold-out (alpha by inner CV within the 4 training folds; pinned folds), routed
endpoint vs the fixed incumbent (m = 75; s = 1), paired, `ST.compare` with the fold-clustered CI.
The router is dead unless the held-out gain is beyond its MDE with the fold CI excluding zero and
5/5 folds the same sign; a label-permutation null (200 draws, same nested fit) is reported
beside it. **If the held-out gain is significantly HARMFUL (fold CI above zero), the file says so
in its first line and C4 stops; two previous routers were.** The s* half is pre-declared as a
control on the method: because s* is provably a function of the invisible common mode, a router
that "predicts" it is predicting the reference, and a positive there would be a leak to be
hunted, not a result.

## Comparison arm, basis, readout

Labels: m* per target from `s12/results/agg_surface.json` (`shipped|m<m>` over the 15 pinned
cardinalities; the argmin over the grid, tie-averaged) and, as the honest alternative, the
split-half construction of S22 `mreal.py` (m chosen on half A of the pool, scored on half B);
s* from `s23/results/errdecomp.json` (`s_star`, closed form). Endpoint: the routed emission
through the shipped readout, point cloud (the basis m and s* were defined on, 3.0483) AND the
built chain through `I.project` (3.2148; note that ideal-geometry projection re-imposes the
3.804 A virtual bond and will undo part of any scale change -- reported, not hidden). The
selection-side arm is carried for the presenter's table only.

## Expected effect vs computed MDE

MDE for the m router: the S23 L7 routed-vs-fixed contrast had SE 0.0196, MDE 0.0550 (point
cloud); the oracle m* gain is -0.244 and its transfer -0.239. MDE for the s router: the errdecomp
s* contrast had SE 0.0561, MDE 0.157; oracle -0.3403. Expected held-out gains: m router 0.00 +/-
0.02 (five routers, all null-to-harmful); s router 0.00 (unreachable in principle). Power: the
design can see a 0.055 A m-routing gain; a router capturing 20% of m* (-0.05) sits at the MDE and
would be Type-M inflated if it appeared.

## Memory, agent-hours

Features from the universe files (one at a time), `s12/cache/disto_*.npz`, `s7/repr_cache/esmcon.npz`;
endpoints from `agg_surface.json` (m ladder already persisted per target; no re-emission needed
for the point-cloud basis) and one `I.project` per routed target for the built chain (2.9 s each,
~6 min). < 0.6 GB. 3 agent-hours.

## Rule 0 forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | ridge (linear) and a depth-2 tree, both nested; the linear one is primary because S22 L10 says nothing richer is learnable at n = 126 | RF, GBT (the overfit signature is already on the record) |
| basis | point cloud (where m*, s* live) and built chain | selection |
| readout | the shipped readout at the routed m or s | any new readout |
| normalisation | features z-scored inside the training folds; the four new blocks entered as four separate routers first, then jointly, all reported | a single joint router only |
| null | fixed m = 75 / s = 1; label permutation; the S22/S23 feature set re-run through the same code as a positive-control-for-the-harness (it must reproduce ~0) | comparing to the oracle m*/s* |
| label | held-out routed CA-RMSD; the oracle and its `best_of_k_within` null beside it | W/L (a best-of-K arm wins by construction) |
