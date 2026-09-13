
---
## ADDENDUM 2 (2026-09-13 09:32) -- `floor2`: THE TRUE MANIFOLD FLOOR, DECLARED BEFORE IT RUNS

Part 2 (`s26/results/ph_cis_floor.json`) measured `floor_ca` = CA-RMSD between the native and
the ideal-trans chain rebuilt from the native's OWN phi/psi: 0.347 A mean, max 1.47 A (1ID6),
Spearman +0.83 with the per-target maximum omega deviation. That number is an UPPER bound on
the representation floor of the production projection: the projection does not rebuild from
the native's torsions, it fits phi/psi to a CA trace, and the nearest ideal-trans chain to the
native is closer than the native's own-torsion rebuild. The brief's second form ("project the
native through the production projection and score against the native, ORACLE DIAGNOSTIC")
is therefore registered here as `python s26/ph_cis.py floor2`: for each target,
`I.project(nat_ca, seq, fold)` (ramah at 0.3, multi-start, exact gradient, the production
call) and its lam=0 rung `fit_ca`, each scored against `nat_ca`; artefact
`s26/results/ph_cis_floor2.json`; compared paired to `floor_ca` and to the chain cost with
`ST.compare`. Expected: `floor2` below `floor_ca` on every target (it is the minimum over the
manifold, of which the rebuild is one point, up to L-BFGS-B convergence), and the lam=0 rung
below the lam=0.3 rung. Cost 126 x ~4.4 s = 10 min CPU, < 0.5 GB, runs after the reject jobs.
