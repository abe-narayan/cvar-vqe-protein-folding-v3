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
