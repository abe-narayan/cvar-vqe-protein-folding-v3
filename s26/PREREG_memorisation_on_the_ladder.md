# PREREG_memorisation_on_the_ladder -- WHERE THE OWN-NATIVE MODELS SIT ON THE S24 PRIOR LADDER (lane W, Sprint 26; L77 extension; ORACLE DIAGNOSTIC)

Written 2026-09-13 22:25, before any gam_eff of a leaked model was computed. Idea file:
`s26/IDEA_memorisation_on_the_ladder.md`. Code: `s26/w_ladder.py`. Results:
`s26/results/w_ladder_*.json`. Every quantity reads the native (gam_eff, cos and MAE are
measured against the native's distances; the endpoint deltas come from L44's gated artefact):
ORACLE DIAGNOSTIC on every line; nothing selects; nothing is deployable.

## 0. What the record knows

S24 L13: interpolating the shipped posterior toward the native (MASS construction, cos = 1 by
construction) moves the point cloud -2.1496 A per unit gamma at the origin (MASS0.1 minus
MASS0.0 = -0.215; `s24/results/priorladder.json`), concave; MASS1.0 reaches 2.226. S25 L12:
across 51 re-readings of the posterior (gam_eff -0.015 to +0.25 at cos up to 0.5),
corr(gam_eff, endpoint delta) = +0.054; a real operator at cos 0.5 travelling 25% is priced at
+0.024 A; "nobody may quote -2.1496 x gam_eff for an achievable operator without the cosine
beside it". L44 Part C: the four pinned models that trained on the target's native move the
built chain -0.698 A (cloud -0.709), 5/5 folds, and their posteriors differ from the clean one
by a mean |dE[d]| of about 0.7 A per pair (1CEK probe). No arm in the record has a large
gam_eff at cos < 1 together with an endpoint.

## 1. Hypothesis and exact falsifier

**H_L.** The own-native models' progress toward the native, (gam_eff, cos) in the S24 MASS
construction (probability space) and in the S25 location currency, predicts their measured
endpoint gain through the ladder's transfer function; i.e. the ladder is a valid currency for
a TRAINED operator, not only for an interpolation.

Quantities per (target, leaked model j != fold), 504 cells, from `s26/p_ladder.progress` (the
same code lane P's C2 rungs use): gam_prob, cos_prob, amp_prob; gam_loc, cos_loc, amp_loc; MAE
of E[d] against the native; and from `s26/results/w_selfcopy_endpoint.json :: C/rows` the
measured deltas (leaked minus clean) on the cloud and the built chain.

Ladder prediction, two readings, both reported: (a) the naive one, -2.1496 x gam_prob (cos = 1
assumed: the reading S25 L12 forbids for achievable operators, listed to show its error); (b)
the direction-discounted one, -2.1496 x gam_prob x cos_prob (the projection of the move onto
the native's direction, the L12 currency). Falsifier: the mean over the 504 cells of prediction
(b) and the measured cloud delta (-0.709) disagree by more than the measured effect's own MDE
(0.244 A, L44). Registered expectation: DISAGREEMENT, with (b) smaller in magnitude than
-0.709 (the trained move lies off the native's direction, cos_prob 0.4 to 0.7, gam_prob 0.15
to 0.35, so (b) is about -0.1 to -0.3): the ladder under-prices a trained operator. If (b)
agrees within the MDE, the ladder is validated as a transfer function for trained priors and
the memorisation is a near-interpolation (gam_prob near 0.7 at cos near 0.9). Secondary,
reported: the per-cell correlation between prediction (b) and the measured delta across the
504 cells (S25 L12's +0.054 is the reference); the same on the built chain; the relation of MAE
to the delta (S7-3: MAE does not price selection).

## 2. Expected effect against the computed MDE

The falsifier's MDE is L44's 0.244 A on the cloud (0.246 built chain). The registered
expectation is a disagreement of 0.4 to 0.6 A, i.e. 2 to 3 MDE.

## 3. Memory, time

Five pinned MLPs and `esm_small.npz` (the envelope run: 0.318 GB peak); 504 posteriors at
~0.1 s each, no projection; the natives' distances from the universes. Under 0.4 GB, about 5
minutes, one governed job. No random element; nothing to replicate.

## 4. Operator forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| currency | S24 MASS construction (probability space) PRIMARY; S25 location space carried | a new currency |
| discount | cos multiplied in (the L12 reading) | the S25 amp ratio; a fitted transfer curve (that would be fitting the ladder to the data) |
| endpoint | the cloud (the ladder's own basis, 3.0483) PRIMARY; built chain carried | selection |
| cells | all 504 (target, model) pairs and their mean | the four-model mean only |
