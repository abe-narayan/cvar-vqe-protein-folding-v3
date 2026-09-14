# PREREG -- PHYSICS CHOOSES THE PROJECTION BRANCH (lane PH, S26; own idea, runs only if ranked)

Written 2026-09-13 09:20, before `s26/ph_branch.py` produced a number. Never edited after;
addenda appended. Idea: `s26/IDEA_branch_select.md`. Code: `s26/ph_branch.py`. Results:
`s26/results/ph_branch_{solutions,relax,report}.json`, cells under
`s26/results/ph_branch_sol_cells/` (projection solutions, CPU) and `s26/results/ph_branch_cells/`
(relaxed energies, AMBER).

## 1. The operator, read from the code

`core.project.lam_path(C, pen, (0.0, 0.3), maxiter=300, multi=True, grad="exact")` (what
`s12.instrument.project` and `core.pipeline.project` call) does this at the production rung:
`cur = fit_multi(C, pen=None, lam=0)` (the four generic starts `STARTS` = alpha (-57,-47), beta
(-139,135), PPII (-75,145), extended (-120,130), lowest objective kept); `warm = fit_prior(C,
cur.phi, cur.psi, pen, lam=0.3)`; `alt = fit_multi(C, pen, lam=0.3)` (the four generic starts
again); the emitted chain is `alt` if `alt.fval < warm.fval` (strict) else `warm`. So every
target has FIVE candidate solutions at the production rung (warm plus four generic), and the
production choice among them is the lowest value of CA-RMSD-to-cloud + 0.3 x ramah penalty.
The projection is degenerate (module docstring: two torsion branches at near-equal objective
distance; the reference disagrees with itself by up to 1.6 A on some targets), so on the
targets where it bites the five solutions differ materially and the objective's choice is a
coin toss weighted by a 2D torsion prior.

The arm: relax each of the five built chains with the production relaxation
(`refine_coords(k_restraint=10, steps=0)`, `s26/ph_c3.relax_chain`), read the converged energy
with the restraint off (`e1`), and emit the BUILT chain (basis `rmsd_arm`) of the solution with
the lowest converged `e1`. Ties in `e1` (identical solutions relax to identical energies) are
averaged over the tied outcome (`ST.argmin_tied`), never broken by array order. If no solution
converges (`e1 > 1000`), the arm emits the production choice. Native-free.

## 2. Gates, before any RMSD

G1: on every target the production choice reconstructed here (`prod_pick`) must equal
`I.project(C, seq, fold)["ca"]` to max |dCA| < 1e-9 A (the same call, the same arithmetic). A
failure stops the run.
G2: the five solutions' distinctness (pairwise CA-RMSD > 1e-3 A) is counted; the number of
targets with at least two distinct solutions is the effective n of the arm and is printed
before any RMSD.
G3: the AMBER probe on one target (`python s26/jobrun.py --agent PH --tag AMBER --name
ph_branch_probe --est-ram 1.0 -- python s26/ph_branch.py probe --pdb 1A13`): five relaxations,
wall and peak RSS quoted from `s26/jobs_done/ph_branch_probe.json`, energies finite, at least
one converged.

## 3. Controls, matched in the operator's space (a choice among the same five)

- anchor: the production choice (lowest objective). Zero information beyond the objective.
- random pick: the expected RMSD of a uniformly random choice among the five, computed exactly
  as the mean over the five (no draws needed); also over the distinct solutions only.
- permuted energy: a uniformly random permutation of five energies picks a uniformly random
  solution, so this control is identical in expectation to the random pick and is not run
  separately; stated here so it is not asked for later.
- raw-energy pick: the same arm on `e0` (the unrelaxed single point of each built chain), a
  declared secondary; the record says `e0` measures the worst clash and I expect it to be worse
  than the relaxed pick.
- ORACLE ceiling: the per-target minimum over the five, an order statistic, quoted only through
  `ST.best_of_k_within` (split-half transfer, k_eff), never as the raw minimum.

## 4. Endpoints (gated) and falsifier

Built chain (`rmsd_arm`) of the chosen solution, ORACLE evaluation of a native-free choice; the
relaxed chain (`rmsd_full`) of the same solution as a secondary. `ST.compare`, paired, fold CI
beside iid, MDE, W/L, concentration null, on all 126 and on the targets with at least two
distinct solutions (n printed).

The idea survives iff the relaxed-energy pick beats the anchor by more than its own MDE with
the fold-clustered CI excluding zero and 5/5 folds in sign, AND beats the random pick by more
than that comparison's own MDE with the fold CI excluding zero. Otherwise the last place in the
pipeline where an all-atom energy could act on a small discrete choice is closed, and I say so.
ORACLE diagnostics reported beside it: how often the objective's choice, the energy's choice and
the per-target best coincide.

## 5. Expected outcome

Null, plausibly a small harmful sign: the energies do not rank real geometry once Rg is
controlled (`s25/results/phys_landscape.json`, in-band rho +0.09 for AMBER), and the relaxed
energy prefers expanded chains (S25: +1.10 A Rg). The oracle ceiling will be real on the targets
where the branches differ (up to 1.6 A) and will mostly fail to transfer under split-half. The
MDE will be 0.02 to 0.06 A depending on how many targets have distinct solutions.

## 6. Cost, memory, checkpoints

Solutions (CPU): one projection per target plus the five fits it already performs, about 5 s per
target, 10 min total, < 0.6 GB, cells `s26/results/ph_branch_sol_cells/<pdb>.json` (phi, psi,
CA, fval, d_to_C of each solution; prod index; G1 deviation; n_distinct). Relaxation (AMBER):
five converged minimisations per target at 5 to 12 s each (S16 measured 12.1 s at k=30 under
contention), 126 x ~45 s = about 1.6 h, one process, est-ram 1.0 GB (S16: ~260 MB per OpenMM
context plus the interpreter), cells `s26/results/ph_branch_cells/<pdb>.json` written after
each target, resumable, never a third AMBER job. Report (gated, CPU): 2 min. Agent-hours: 3.

## 7. Replication

A positive result is re-run with the fold processing order reversed and must land inside its own
fold CI (the arm has no random draws; the relaxation is deterministic at threads = 1).

---
## ADDENDUM 1 (2026-09-13 19:55) -- THE GATES PASSED AND THE PROBE, BEFORE THE RUN

G1: the reconstructed production choice equals `I.project` to max |dCA| = 0.0 on 126/126
(`s26/results/ph_branch_solutions.json`, `summary.g1_max` 0.0; job `ph_branch_solutions`,
1102 s, peak 0.088 GB). G2: 4.77 distinct solutions per target on average; ALL 126 targets have
at least two distinct solutions, so the effective n of the arm is 126; the production choice is
the warm start on 60 targets and one of the four generic starts on 66 (counts [60, 6, 34, 14,
12]). G3, the AMBER probe on 1A13 (job `ph_branch_probe`, exit 0, wall 70 s, PEAK RSS 0.272 GB,
held 510 s at the job cap): five relaxations, all converged, e1 -306 to -341 kcal/mol, 11.5 to
16.1 s each, 65 s per target. Extrapolated run: 126 x 65 s = about 2.3 h, est-ram 0.4 GB
(probe peak 0.272 rounded up), one job, `--tag AMBER`, per-target cells. Launched at 19:56 with
3.9 GB free (rule: probe peak + 0.5 GB = 0.77 GB under the free reading).
