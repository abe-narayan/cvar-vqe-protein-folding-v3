# PREREG_amber_prior_partner -- AMBER AS A DISTRIBUTION INSIDE THE PRIOR: THE ENERGY-WEIGHTED POOL HISTOGRAM AS THE MIXTURE PARTNER (lane P's H_C3a, orphaned to lane W; tournament rank 10)

Written 2026-09-13 22:15, before any code ran on a real target. Idea file:
`s26/IDEA_amber_prior_partner.md`; parent draft `s26/PREREG_C3.md` (lane P; H_C3a only, H_C3b
withdrawn there by S8-14's arithmetic). Code to be written: `s26/w_amberprior.py` (reuses
`s26/p_ladder.py`'s `pool_histogram`, `risk_table`, `posterior("mix")` semantics through
`s26/w_selfcopy.emit`). Results: `s26/results/w_amberprior_*.json`. The native-free half
(`clouds`) stores the emissions for every (lam, beta) cell; the gated half (`endpoint`) reads
natives, chooses lam and beta leave-fold-out on the built chain, and runs `ST.compare`.

## 0. What the record knows

- The functional lever is closed in all five forms in which AMBER acts on the CANDIDATES'
  ranking (filter, partition, additive score, per-target sign, ORACLE weight worth 0.015 A: S24
  L16, S25 L18); AMBER's rank-permuted control is -0.0010 at 0.08x MDE (S24 L16): its ordering
  carries no accuracy through the score; AMBER's top-75 is +1.10 A more expanded than the
  distogram's (S25 L16). Here AMBER acts on the PRIOR's per-pair location and width, and the
  shipped score then ranks candidates against that prior: a different operator space, never
  tried (lane P's idea file).
- The unweighted pool histogram as a mixture partner is rung `mix` of PREREG_C2 (lane P), itself
  foreshadowed null (S7-3 `gain 0`, S10-2, S19 L14: every move toward the pool is worse). If
  lane P's `mix` result is on the ledger before this runs, it is the primary comparator; if not,
  this file's beta = 0 arm IS the unweighted mixture and serves as that comparator.
- The 63,000 ff14SB/GBn2 single points are cached (`s24/cache_amber/<pdb>.npz`: `e_amber` (500,)
  aligned with `universe_idx` == `I.pool_idx`, asserted per target; verified fresh single points
  at max_rel 0.0, S25 L16 F0), so no OpenMM call is needed.

## 1. Hypothesis and exact falsifier

**H_AP.** With H_beta(pair, bin) = the K = 500 pool's per-pair 17-bin distance histogram weighted
by w_c = exp(-beta * zrank(E_AMBER)_c) (rank-standardised, the only monotone normalisation, S24
L10 / S25 L16), the mixed posterior P = (1 - lam) P_shipped + lam H_beta, through the shipped
Bayes-risk functional, top-75, medoid-frame average and projection, beats the shipped posterior
on the built chain.

**Falsifier (PREREG_C2's, as the idea file states).** (lam*, beta*) chosen leave-fold-out on the
built chain (the four training folds pick the cell, the fifth is scored; `arm` PRIMARY) must beat
the shipped posterior by more than its own MDE with the fold-clustered CI excluding zero and
5/5 folds the same sign, AND beat the UNWEIGHTED mixture (beta = 0 at its own leave-fold-out
lam) by the same standard, AND beat the rank-PERMUTED-AMBER mixture (the energies permuted across
the 500 candidates within each target, marginal preserved, correspondence destroyed; 4 seeded
draws, `s15.seed.stable_rng("w_amberprior", pdb, k)`) at the same (lam*, beta*) by the same
standard. lam = 0 (any beta) and beta = 0 (any lam) are asserted bit-exact against the incumbent
risk table and against `p_ladder`'s `pool_histogram` respectively, per target. Otherwise H_AP is
refuted at this instrument with its power stated; the per-target oracle over the (lam, beta)
grid is reported beside `ST.best_of_k_within` and its split-half transfer, never as a result.
Registered expectation (lane P's): 0.00 +/- 0.03, indistinguishable from the unweighted mixture;
the design can refute a 0.1 A gain and cannot see a 0.03 A one.

## 2. Grid, arms, operator

    lam  in {0, 0.05, 0.1, 0.2, 0.35, 0.5}      (PREREG_C2's MIX_LAMS without 0.75 and 1.0, which are typicality: S7-3)
    beta in {0, 0.5, 1, 2, 4}                    (PREREG_C3's grid, on the zrank scale)
    H_beta(pair, c) = sum_k w_k [bin(d_pair(k)) == c] / sum_k w_k + eps, normalised per pair (eps 1e-3, p_ladder's convention)
    arms per target: 30 cells; the permuted-AMBER control at every cell with beta > 0 for 4 draws
    emission per cell: shipped score from a genuine core.predict.Distogram (risk table), top-75, cloud, one projection

Leave-fold-out choice: for held-out fold f, (lam*, beta*) = argmin over the grid of the mean
built-chain RMSD on the other four folds (ties: the cell nearest (0, 0) in grid order, stated);
the held-out targets take that cell. The same for the beta = 0 row alone (lam*_0).

## 3. Statistics

`ST.compare(arm, shipped, folds=ST.pinned_folds(pdbs), names=pdbs)` for: the LFO cell vs shipped;
the LFO cell vs the beta = 0 LFO cell; the LFO cell vs the permuted control at the same cell
(mean over 4 draws); `ST.fmt` verbatim; arm PRIMARY, cloud and sel carried on both sides; the
per-target oracle over 30 cells through `ST.best_of_k_within` with k_eff and the split-half
transfer. Replication for a positive: a second permutation seed set and the reversed
fold-processing order. Power for a null: the MDE of the LFO cell vs shipped.

## 4. Expected effect against the computed MDE

MDE 0.05 to 0.09 A on the built chain for a small prior change (PREREG_C2 section 4: sd 0.207
on the cloud for MASSFIXW0.1 minus MASS0.0, plus projection noise); expected effect 0.00 +/- 0.03;
the deliverable is the power statement plus the (lam*, beta*) the LFO choice makes (the
registered expectation is lam* = 0 on most folds, i.e. the choice declines to mix).

## 5. Memory, time

Per target: 500 x npairs distances (already in `I.pair_dists`), 30 histograms and risk tables,
30 + 4 x 24 = 126 emissions, of which the PROJECTION is the cost (3.5 s each): 126 x 126 x 3.5 s
= 15 h. Too much. Declared reduction, before running: project only the cells the falsifier needs,
by staging: stage 1 (native-free) stores every cell's CLOUD (no projection; seconds per target)
and the top-75 membership; stage 2 (gated) chooses (lam*, beta*) leave-fold-out on the POINT
CLOUD, then projects only the chosen cell, the beta = 0 chosen cell and the 4 permuted draws at
the chosen cell per target (6 projections, ~21 s per target, 45 min for 126). The built chain
stays PRIMARY for the verdict; the LFO choice is made on the cloud, which is stated as a fork
(NOT TAKEN: choosing on the built chain at 15 h). est-ram 0.5 GB (the same loads as
`w_selfcopy.retrieval`, 0.116 GB measured, plus one 500 x npairs x 17 histogram tensor per cell).
One-target probe first. Agent-hours: 2 code, 1.5 run, 1 write-up.

## 6. Operator forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| energy normalisation | zrank of the raw single points (monotone) | moment z-score (non-monotone in float64, S25 L16); raw Boltzmann weights (a 1e28 range) |
| histogram | p_ladder's `pool_histogram` with weights, same bins, same eps | kernel-smoothed histogram |
| choice of (lam, beta) | leave-fold-out on the point cloud (cost), verdict on the built chain | on the built chain (15 h); per-target beta |
| control | rank-permuted AMBER at the chosen cell, 4 draws; the beta = 0 mixture | S25's random-75 null (a readout-level control, not a prior-level one) |
| basis | built chain PRIMARY; cloud and sel carried | AMBER energy of the emission |

## ADDENDUM 1 (2026-09-13 22:50) -- staging refined before the first real-target run; nothing above edited

- Stage 1 stores the 30 cells' clouds only; the medoid frame (the cost, 75 x 75 Kabsch) is cached
  by top-75 set within a target, and `cell_cloud` is checked on synthetic data to reproduce
  `w_selfcopy.emit`'s cloud exactly. The rank-permuted-AMBER control (4 seeded draws) is computed
  in stage 2 at the chosen (lam*, beta*) cell only, which is the only place the falsifier reads it.
  Nothing else in sections 1 to 6 changes.
