# IDEA -- Physics chooses the projection BRANCH, not the candidate (lane PH; own idea 1)

## Hypothesis
The projection is degenerate: a CA point cloud admits two ideal-geometry torsion solutions at
near-equal objective distance, one Ramachandran-plausible and one not, and a warm-started
optimiser cannot cross between them (`core/project.py` module docstring; "the reference
disagrees with itself by up to 1.6 A on the same targets"). Production resolves this by
multi-start from four generic starts (`core.project.STARTS`, alpha, beta, PPII, extended) at
each rung and keeps the lowest objective (`fit_multi`, `lam_path`). The objective is
CA-RMSD-to-cloud plus 0.3 x a hinged Ramachandran prior; the choice among near-degenerate
branches is therefore made by a 2D torsion prior, never by an all-atom energy. The hypothesis is
that the CONVERGED restrained AMBER energy of each branch's built chain (E after
`refine_coords(k=10, steps=0)`, the operator S20 L6 showed is what makes the AMBER objective
defined) selects the physically plausible branch more often than the objective does, and that
this lowers the built-chain RMSD.

This is a different operator space from every physics use on the record: the physics never
ranks pool candidates (closed, S25), never scores the average (closed, S24/S25), never moves
atoms toward accuracy (closed, S16). It picks ONE of four to eight discrete solutions of the
projection on a single target, where the candidates differ by up to 1.6 A and share the same
input.

## Why the record does not already close it
- S13/S14 (`structural-objective-beats-the-energies`): in torsion space the energies' global
  optimum is in the wrong place and the structural objective's is in the right place. That is a
  statement about optimising the energy, not about choosing among solutions of the structural
  objective.
- S16 L27: AMBER's displacement is worse than random. Here AMBER moves nothing; the emitted
  structure is one of the projection's own solutions.
- S25 Form 5: physics supplies a per-target sign to switch between two configurations; closed at
  0.14x MDE. That switched between two near-identical arms (corr 0.95, k_eff 1.02). The branch
  solutions here differ by up to 1.6 A on the targets where the degeneracy bites, so the ORACLE
  ceiling is not a noise artefact by construction; it must still be quoted through
  `best_of_k_within` because it is a per-target minimum.
- No sprint has looked at the multi-start solutions individually; `fit_multi` returns the argmin
  and discards the rest.

## Exact falsifier
Native-free arm: for each target, take the lam=0.3 solutions from the four generic starts plus
the warm start (five candidates, all emitted by `fit_prior` with the production settings),
relax each with the production operator, pick the lowest converged energy (ties by
`ST.argmin_tied`), emit its BUILT chain (not the relaxed one, so the basis is `rmsd_arm`).
Controls: (a) the production choice (lowest objective), the anchor; (b) the same choice made on
a permuted energy vector across the five (marginal kept), 16 draws; (c) a random pick among the
five, 16 draws. ORACLE ceiling: the per-target minimum over the five, quoted with
`ST.best_of_k_within` and the split-half transfer. The idea survives iff the energy-chosen arm
beats the anchor by more than its own MDE with the fold CI excluding zero and beats the random
pick by more than that comparison's MDE. Secondary diagnostic (ORACLE): how often the
objective's choice and the energy's choice coincide, and how often each picks the per-target best.

## Expected effect against the MDE
The multi-start finds a strictly lower objective on some targets and the branches are within
noise on most, so the number of targets that move will be well under 126 and the MDE will be
0.03 to 0.08 A. Expected effect: null (the energies do not rank real geometry once Rg is
controlled, `s25/results/phys_landscape.json` in-band rho +0.09), with a plausible small
harmful sign because the relaxed energy prefers the more expanded branch. Plausibility 0.2.

## Cost
Five projections' solutions per target are produced by the existing code path at no extra cost
beyond one projection (~4.4 s). Five converged relaxations per target at ~6 to 12 s each: 126 x
5 x 9 s = 1.6 h under `--tag AMBER`, ~1.0 GB, checkpointed per target. Agent-hours: 3.
Probe one target first.

## Information value
High if positive (a first accuracy role for the all-atom energy, on a mechanism the record
predicts and has not tested); if negative, it closes the last place in the pipeline where a
physics energy could act on a small discrete choice rather than a ranking.
