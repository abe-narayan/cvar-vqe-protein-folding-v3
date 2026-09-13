# PREREG -- C3, THE MATCHED-RANDOM CONTROL FOR AMBER REFINEMENT (lane PH, S26)

Written 2026-09-13 00:40, before `s26/ph_c3.py` read any RMSD. Never edited after; addenda
appended. Proposal C, Part 3.3 of the campaign prompt. Code: `s26/ph_c3.py`. Results:
`s26/results/ph_c3_{nativefree,stage1,stage2}.json`, cells under `s26/results/ph_c3_cells/`.

## 1. The question

The production pipeline ends with one restrained ff14SB/GBn2 minimisation of the built chain
(`core.pipeline.relax`, `refine_coords(k_restraint=10, steps=0)`, converged). On the 126 dev
targets it moves the built chain (`rmsd_arm` 3.2148) to the relaxed chain (`rmsd_full` 3.2355),
+0.0207 A [+0.0143, +0.0276] (state brief section 4; `bench_results/cache/1fc9f2dcf489e2fb`).
S16 L27 measured, on a different input (the all-atom coordinate average at k=30), that a random
displacement of matched magnitude is at least as accurate as AMBER's (AMBER minus random +0.0225
[-0.0138, +0.0571]) and that AMBER's displacement has ORACLE cos -0.052 [-0.092, -0.013] with the
true residual. C3 asks whether the production step's +0.0207 A is what any move of that size
costs, on the production input and the production operator.

## 2. Stage 1 (no AMBER compute; after the gate)

Per target, from the persisted coordinates: the displacement v = superpose(amber_ca onto ca) -
ca and its per-atom RMS magnitude; the ORACLE residual r = superpose(nat onto ca) - ca and
cos(v, r). Controls, matched in the operator's space:

- `rand`: S16's construction reproduced from `s16/repair_report.py:228-238`: isotropic Gaussian
  direction on the CA trace, six rigid-body components removed by projection onto
  `s15.align_lib.rigid_basis(ca)`, scaled to the same per-atom RMS magnitude; 16 draws, mean.
- `member`: S20's `c_land_null` "toward another pool member" form in coordinate space: the same
  RMS magnitude along the straight line from `ca` toward a randomly chosen member of the shipped
  K=500 pool superposed onto `ca`; 16 draws, mean; the number of draws where the member is closer
  than the magnitude (overshoot) is recorded.
- zero-information anchor: `ca` (do nothing), the production result.

Differences from S16, stated: input is the built chain, not the all-atom average; k = 10 (the
production constant), not 30; 16 draws, not 3; the magnitude is measured on CA after
superposition (S16's `moved_ca`, the same definition).

Endpoint: full-chain CA-RMSD to the native (ORACLE evaluation of native-free operators).
Contrasts, all paired, `ST.compare`, fold CI beside iid, MDE, W/L, concentration null:
AMBER minus do-nothing (must reproduce +0.0207); AMBER minus `rand`; AMBER minus `member`;
`rand` minus do-nothing; `member` minus do-nothing; ORACLE cos(AMBER) minus cos(rand) and minus
cos(member). Also the S16 identity n RMSD_after^2 = |r|^2 - 2 v.r + |v|^2 (an upper bound on the
superposed RMSD), and the orthogonal-move cost it predicts from the magnitude alone.

Validity axis, native-free, from the cache: `amber_e0` (the built chain's own AMBER energy, its
strain before relaxation), `amber_e1`, `amber_strain_after` (bond+angle after), `amber_moved`
(restraint RMSD on N/CA/C), the convergence flag (e1 <= 1000 kcal/mol, `core.amber.
CONVERGE_MAX_KCAL`), and a CA-level axis for `ca` and `amber_ca` (virtual bond mean/min/max, CA
contacts below 4.0 and 3.8 A at |i-j| >= 3, Rg). The heavy-atom axis (clash count, bond
deviation) requires the relaxed backbone, which the cache does not hold; stage 2 measures it
with `s16.energy_lib.panel` on the input and output of `refine_coords`.

## 3. Falsifier and decision rule

"Refine with physics" is an ACCURACY step iff AMBER beats BOTH matched controls (`rand` and
`member`) by more than each comparison's own MDE with the fold-clustered CI excluding zero and
5/5 folds in sign. Otherwise the presentation drops "refine with physics" as an accuracy step
and keeps it as a validity step, and `s26/C3_RESULT.md` says so with the numbers.

Expected outcome, derived rather than guessed: S16's identity says a move of magnitude m
orthogonal to the residual costs about m^2 / (2 RMSD). The production step moves the CA trace
by about 0.15 A (`amber_moved` 0.156 on 1A13 in the cache; the census will give the mean), so an
orthogonal move of that size costs a few thousandths of an Angstrom on a 3.2 A chain. The
production step costs +0.0207. If the census magnitude is ~0.15 A, the prediction is that AMBER
is WORSE than the random control by roughly +0.01 to +0.02 A, i.e. its displacement is not
orthogonal to the residual but points slightly away from it, reproducing S16's negative cosine.
I state this so it can be wrong.

## 4. Stage 2 (AMBER compute; when lane P delivers the best C2 rung)

Input: a JSON with rows `pdb`, `phi`, `psi` (radians, the built chain of the rung). Operator:
`relax_chain` in `s26/ph_c3.py`, which mirrors `core.pipeline._relax_inner` (ideal backbone from
the torsions, `torsion_lib2.library_for(seq, 8, seq)`, `PerResidueTorsion`, `refine_coords(k=10,
steps=0, components=True)`), then the identical stage-1 controls on that input, plus the
heavy-atom validity panel on input and output. Probe one target first under
`python s26/jobrun.py --agent PH --tag AMBER --name ph_c3_probe --est-ram 1.0 -- python
s26/ph_c3.py probe --input <file> --pdb <pdb>` and quote the peak RSS and wall time. Then the
126 under `--tag AMBER`, checkpointed per target to `s26/results/ph_c3_cells/<pdb>.json`,
resumable. Never a third AMBER job.

Memory: S16 measured ~230 to 260 MB per OpenMM context plus the interpreter (`core.amber.
WORKER_MB`); the estimate is 1.0 GB for one process. Wall: S16 measured 12.1 s per converged
k=30 minimisation at threads=1 under contention; expect 126 x (5 to 15 s) = 10 to 30 min.

## 5. Cost

Stage 1: ~2 min CPU, < 0.5 GB. Stage 2: ~30 min AMBER, ~1.0 GB. Agent-hours: 1.5 + 1.

## 6. Replication

A positive result (AMBER beating both controls) is re-run with different draw seeds and the fold
order reversed and must land inside its own fold CI.
