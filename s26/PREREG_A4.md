# PREREG A4 -- GRADIENT VARIANCE OF ADAPT-GROWN CIRCUITS BESIDE THE FIXED ANSATZ (lane Q)

Filed 2026-09-13 08:45, before `s26/q_var.py` has been run above n = 6. Not edited after;
addenda are appended. Property measurement: reads no native, no RMSD. The energies are the
deployed energy SHAPE (standardised ranks 1..2^n), not any target's score, so no pool is opened.
The coordinator ruled (08:40) that it may run before "PHASE 0 SIGNED OFF".

## Why it exists

`s25/q_plateau.py` measured Var_theta[dF/dtheta_0] for the fixed RY/CNOT ansatz at layers 3,
n = 4..13, theta ~ N(0, 0.6^2), five (alpha, T) cells, and found no exponential plateau
(`s25/results/q_plateau.json`: log2 Var per qubit -0.649 linear cost, -0.252 alpha=0.25,
-0.047 alpha=0.10, -0.311 and -0.243 for the deployed form). Proposal A replaces the fixed
ansatz by an ADAPT-grown one. The same measurement on the grown circuit, at matched parameter
count, says whether growth changes the variance picture, and in which direction.

## Design (extends S25 exactly; nothing in the fixed arm is re-implemented)

1. Reproduction gate. `s25.q_plateau.measure(7, 3, alpha, T, 250, seed=1007, init_sd=0.6)`
   must reproduce the n = 7 row of every cell in `s25/results/q_plateau.json` to relative
   1e-9 before anything else runs (asserted; the smoke run passed it, worst relative
   deviation recorded in `s26/results/q_var_smoke_n46.json`).
2. Fixed ansatz: the S25 code path at n = 4, 6, 7, 8, 10, 12, 13, draws 250/250/250/250/
   200/120/80, seed 1000 + n, so the table can be reprinted from this sprint's artefact.
3. Grown circuits: for each pool in {V (Tang, 2n-2), L2 (2-local odd-Y, 2n^2-n)} and each cell
   in {(1,0), (0.25,0), (0.10,0), (1,0.3), (0.25,0.3)}, `run_adapt` grows on the deployed shape
   under that cell (seed 0, Adam best-iterate lr 0.15, 50 steps per growth step, eps 1e-3,
   max P = 3n). Then theta ~ N(0, 0.6^2) is drawn for the grown circuit's OWN P parameters
   (seed 1000 + n, same draw counts) and Var[dF/dtheta_0] and E|g|^2/P are measured with the
   exact parameter-shift gradient. theta_0 is the first RY angle in both ansaetze.
4. Report fitted log2 Var per qubit for fixed and grown side by side, per cell, and the
   actual P of each grown circuit (early eps stops are labelled, not hidden).

## What the smoke run at n = 4, 6 already showed (`s26/results/q_var_smoke_n46.json`)

- At T = 0 with alpha < 1 ADAPT converges toward the CVaR collapse, a computational basis
  state, where every first-order gradient vanishes exactly (the lemma in `s26/q_adapt.py`:
  p_x(phi) = cos^2(phi/2) delta + sin^2(phi/2)|<x|P|0>|^2). It stops by eps at P = n, and the
  variance of the grown circuit is 0 (alpha = 0.10) or 4e-5 (alpha = 0.25). Those cells are
  degenerate for a grown circuit by construction and are reported with that label.
- At alpha = 1 (T = 0 or 0.3) the grown circuit is single-qubit rotations only (the optimum is
  a product state, see PREREG_A1 section 2) and its variance is 5 to 10 times the fixed
  ansatz's (n = 6: 0.14 vs 0.017 linear; 0.10 vs 0.012 deployed form).

## Hypotheses and predictions

H4a. Deployed-form cells (T = 0.3): the grown circuits' log2 Var per qubit is FLATTER (closer
     to 0) than the fixed ansatz's, because at alpha = 1 they are product circuits (abelian
     DLA, no plateau possible) and at alpha = 0.25 they carry few entangling 2-local strings.
     Predicted: grown slope in [-0.3, 0] for both pools at both T = 0.3 cells; fixed -0.31 /
     -0.24 (S25).
H4b. The grown circuits' variance MAGNITUDE exceeds the fixed ansatz's at every n where P is
     matched, by a factor > 2, for the same reason.
H4c. The T = 0, alpha < 1 grown cells are degenerate (P < 3n, variance < 1e-3) at every n.
H4d. The fixed-ansatz rows reproduce S25 at n = 7 to 1e-9 and re-derive the S25 slopes to
     within the sampling error of a variance estimate (relative SE sqrt(2/(m-1)), 9% at 250
     draws, 16% at 80).

## Falsifiers

- H4a is falsified if either grown slope at a T = 0.3 cell is below -0.5 log2 per qubit with
  P matched at every n (then growth makes the landscape steeper in n, which would need its
  own explanation and would be reported as such).
- H4b is falsified by any matched-P n where the grown variance is below the fixed one.
- H4d failing stops the run: the S25 code path is not being reused and nothing is reported.
- A grown circuit whose P is not matched at some n is excluded from the slope fit and the row
  is listed; it does not count for or against H4a.

## What this can and cannot say (to be repeated in the findings)

A large gradient variance from a product circuit is not trainability in a useful sense; it is
the trivial regime, and by Cerezo et al. 2025 (arXiv:2312.09121) the extreme case of a
classically simulable landscape. "No barren plateau" for a grown circuit will not be written
as a merit. Nothing here is a claim about hardware, noise, shot cost or advantage
(`s26/LANE_CONTRACT.md` rule 10).

## Memory, time, governor

Smoke run (n = 4, 6, pool V): peak RSS 0.049 GB, 10 s (`s26/jobs_done/q_var_smoke.json`).
At n = 13 a shifted batch is 78 x 8192 doubles (5 MB); growth with the L2 pool (325 strings)
costs about 39 growth steps x (325 pool gradients + 50 Adam steps x 78 shifted circuits) of
8192-amplitude circuits with up to 39 rotations; estimate 3 min per (pool, cell) at n = 13,
under 1 h in total; peak RSS estimate < 0.3 GB. Launch: `python s26/jobrun.py --agent Q --tag
CPU --name a4_var --est-ram 0.4 -- python s26/q_var.py`. Checkpoint: `s26/results/q_var.json`
rewritten atomically after every (pool, cell, n) row; resumable. Agent-hours: 1 compute
(waiting), 1 writing.

## Rule 0

A property measurement with no outcome that can favour an endpoint hypothesis; falsifiers
stand in place of the six forks. The one choice that could flatter one arm is the initial
law for the grown circuit's parameters; it is the same N(0, 0.6^2) law S25 used for the fixed
one, and the alternative (drawing only the grown angles and keeping the RY layer at its
trained values) is NOT TAKEN because a barren-plateau statement is about random
initialisation.

## ADDENDUM (2026-09-13 09:40, after the run; the text above is unchanged)

`s26/results/q_var.json` landed (2,765 s, peak RSS 0.08 GB, `s26/jobs_done/a4_var.json`).
H4d held exactly (relative deviation 0.0 on all five n = 7 rows). H4a: matched-P grown slopes
at T = 0.3 are +0.035 (V, 5 rows), -0.008 (L2, 7), -0.246 (V, 6), -0.302 (L2, 7); the
falsifier (below -0.5) did not fire; two values sit just outside the predicted [-0.3, 0] by
0.035 and 0.002. H4b failed on 1 of 32 matched rows (n = 4, alpha = 1, T = 0.3: ratio 1.63)
and held on 31 (2.5 to 370). H4c: P < 3n held on 14 of 14 rows; "Var < 1e-3" failed on 3 of 7
rows at alpha = 0.25, T = 0 (2.4e-2, 1.0e-3, 1.2e-3). The printed side-by-side slopes in
`s26/logs/a4_var.log` include the early-stopped rows; the matched-only slopes above are the
pre-registered quantity and are the ones in the ledger (L35) and the figure. Ledger entry L35.

## ADDENDUM 2 (2026-09-13 22:35, BEFORE the run; bootstrap CIs on the slopes, the L47 caveat)

L47 records that `q_var.json` stores point slopes and no CI, so "-0.302 is -0.243 within
error" is asserted, not computed. Method: re-measure every (cell, n) row of A4, fixed and
grown, with the SAME draws (seed 1000 + n, theta ~ N(0, 0.6^2), the same draw counts) but
STORING every per-draw dF/dtheta_0; assert that the recomputed Var equals `q_var.json`'s
`var_g0` to relative 1e-12 for every row (the reproduction gate of this addendum); then a
percentile bootstrap over the theta draws within each n (resample the draws with replacement
at every n independently, recompute Var, refit the 7-point log2 slope; 2000 resamples, seed
2026), for every fixed and grown slope, and for the DIFFERENCE grown - fixed per cell (the two
circuits' draws are independent, so the difference resamples both). Matched-P rows only for
the grown slopes, as in L35. Code `s26/q_var_boot.py`, artefact `s26/results/q_var_boot.json`
(per-draw g0 arrays, the CIs), figure `s26/figures/a4_variance_slopes.png` regenerated with
error bars if the CIs exist.
Predictions: P4a the 95% CI of (grown L2 - fixed) at alpha = 0.25, T = 0.3 includes zero;
P4b the fixed-ansatz CIs are narrower than +-0.15 log2 per qubit at every cell; P4c every
alpha = 1 grown slope's CI excludes the fixed slope's point value at the same cell (the
"no decay" statement carries an interval). Falsifiers: the stated intervals, measured. The
grown circuits are re-grown deterministically (seed 0) and must reproduce their
`sequence` from `q_var.json` exactly (asserted). Time: the A4 run took 2,765 s; this one
repeats it with storage, about 50 minutes, < 0.3 GB. Launch through jobrun as `a4_var_boot`.

## ADDENDUM 3 (2026-09-14, after the addendum-2 run; ledger L119)

`s26/results/q_var_boot.json` landed (3,928 s, 0.09 GB). Reproduction exact (relative 0.0 on
35 fixed and 70 grown rows). P4a HELD (grown L2 - fixed at alpha = 0.25, T = 0.3: -0.056
[-0.182, +0.087]). P4b held on 4 of 5 cells and FAILED on alpha = 0.10, T = 0 (upper half-width
+0.20). P4c HELD (alpha = 1 differences +0.30 to +0.66, all CIs excluding zero).
