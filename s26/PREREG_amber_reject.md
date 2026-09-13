# PREREG -- AMBER AS A STERIC REJECT FILTER, NOT A RANKER (lane PH, S26)

Written 2026-09-13 00:40, before any RMSD in `s26/ph_reject.py` was read. Never edited after;
addenda are appended at the end. Mandatory tournament direction (campaign prompt Part 4.2).
Code: `s26/ph_reject.py`. Results: `s26/results/ph_reject_{census,cloud,chain,report}.json`.

## 1. Hypothesis

Both physics energies rank worse than noise when used as a selector (S25 PHYS section 1.4 and
1.5: Legacy +0.330 A and AMBER +0.455 A against a random 75-subset, 5/5 folds; AMBER inside a
distogram-led score is indistinguishable from its own rank-permuted twin at 0.43x MDE). What the
AMBER single point does measure is its worst clash: 58.6% of every pool is above 1e4 kcal/mol
(`s25/results/phys_landscape.json`, `AMB_frac_absz_lt_0p1` 0.970, `AMB_top10_var_share` 0.9966),
and on unrelaxed ideal-geometry rebuilds the energy is finite-but-meaningless above ~1e4 (S20 L8,
`s20/results/c_q1.json`: median +16,062 kcal/mol, 53.5% above 1e4).

The hypothesis is that a BINARY reject at a FIXED PHYSICAL THRESHOLD, which uses only the one
bit the energy is known to carry ("this candidate contains a hard overlap") and none of the
ordering it is known not to carry, improves the production endpoint.

## 2. How this differs from the filters already closed (the Adversary's question)

| earlier test | what it rejected | count set by | refill | basis | n |
|---|---|---|---|---|---|
| S19 Q3 `s19/agentC_reject.py` (`s19/agentC_FINDINGS.md` section 3) | the r worst of the top-75 by Legacy steric / min-heavy / AMBER single point, r in {2,5,10,19} | rank (fixed r on every target) | no | point cloud | 126 |
| S24 D1-C (`s24/LEDGER.md` L10-A, "pre-filter arm") | a fixed fraction of the K=500 pool by which energy PREFERS a candidate (five partitions) | rank (fixed fraction) | n/a (the score selected inside the partition) | point cloud | 30 |
| S25 suite C2/C5 (`s25/results/phys_suite.json`) | nothing; AMBER as a score/rank channel | n/a | n/a | point cloud | 126 |
| **this test** | every top-75 member above a physical energy threshold | **physics** (zero on a clean target, large on a clash-ridden one) | **both arms**: refill to m=75 from the next-ranked survivors (R), and no refill (S) | **point cloud AND built chain** | 126 |

S19's AMBER row at r=10 was +0.0109 [-0.0049, +0.0285] against matched-random (not measured) and
+0.0138 against no rejection; its `steric` and `minheavy` rows were measurably worse than random.
The record therefore says that rejecting a fixed rank-count is null-to-harmful. It does not say
what happens when the count is chosen by the physics and the retained set is refilled, and it
was never measured on the built chain. That is the whole of the difference, and if the
difference does not matter the result will say so.

## 3. Operators (all native-free), thresholds, controls

Energies: `s24/cache_amber/<pdb>.npz` `e_amber` aligned to `universe_idx`; asserted on every
target equal to `I.pool_idx(u)` bit-for-bit and `amber_verify_max_rel == 0.0`; the production
`sub` asserted equal, as a set, to the first 75 of the stable argsort of the cached distogram
score. No new AMBER single points. Ties: the score order is `argsort(kind="stable")`, the
production tie-break (`core.pipeline.Config.tie_break`); the threshold test has no ties beyond
byte-identical duplicate windows, which are rejected or kept together.

Threshold family: T in {1e3, 1e4, 1e5, 1e6} kcal/mol. **PRIMARY T = 1e4**, chosen before the
run because it is the value S20 declared before its own run as the boundary of the "physical
subset" (three orders above the relaxed range of -1170 to -500 kcal/mol) and because S25
measured 58.6% of every pool above it. The other three are shape; any claim about "the best
threshold" is an order statistic and is quoted through `ST.best_of_k_within` (split-half
transfer, k_eff), never as the raw minimum.

Arms per T:
- `R_T`  reject e > T, refill from the next-ranked survivors until m = 75 (identically, the
  top-75 of the pool's survivors in the shipped score order); `R_short` recorded if the pool
  has fewer than 75 survivors.
- `S_T`  reject e > T from the shipped top-75, no refill; an empty set falls back to the anchor
  and is counted.
- `RANDS_T`  control matched to S: reject the SAME COUNT at random from the top-75, 16 draws.
- `RANDR_T`  control matched to R: reject the same count at random, refill with the next-ranked
  candidates in score order with no energy judgment, 16 draws.
- `PERMS_T`, `PERMR_T`  control: the energy vector permuted within the pool (marginal kept,
  correspondence destroyed), same threshold, same two operators, 16 permutations.
- anchor: the shipped top-75 (3.0483 point cloud, 3.2148 built chain), zero information.

Readout: uniform coordinate average in the retained set's own medoid frame
(`I.coordinate_average`), then the production projection (`I.project`, ramah at 0.3,
multi-start, exact gradient) for the built chain. The anchor must reproduce `rmsd_avg` to 1e-6
and the anchor's projected chain must reproduce the production `ca` (`anchor_vs_prod` is the CA
RMSD between them; the max over targets is reported).

Built-chain coverage: anchor, R_T and S_T at every T; the four controls at the PRIMARY T with
4 draws each (a projection costs ~4.4 s; 16 draws x 4 controls x 126 targets is 12 hours of
CPU and is not affordable under a shared 4 GB box). A retained set identical to another set on
the same target reuses its projection. This is stated here, not discovered later.

## 4. Falsifier

The steric reject is a live accuracy step iff, on the BUILT CHAIN at the primary threshold,
`R_1e4` or `S_1e4` beats the anchor by more than its own MDE (2.8016 x SE, computed per
comparison by `ST.compare`) with the fold-clustered CI excluding zero and 5/5 folds in sign,
AND beats its matched control (`RANDR_1e4` or `RANDS_1e4` respectively) by more than that
comparison's own MDE with the fold CI excluding zero.

If it does not, the direction is closed on this instrument in the two forms not previously
measured (physics-set count, refill), on both bases, and I say so. A point-cloud gain without
a built-chain gain is reported as a point-cloud gain and is not the result.

Secondary, with multiplicity stated: the same test at the three other thresholds. A threshold
other than 1e4 that clears the bar is reported as shape and is believed only if the split-half
transfer of the R (or S) grid is at least 50% of its oracle gain.

## 5. Expected outcome, stated before the run

Prior: **null or harmful**. Reasons on the record: S24 D1-C (five physics pre-filters, all
harmful in sign); S19 Q3 (fixed-count steric rejection worse than random with a sign); S23 L5
(removing four geometrically deviant members from the average costs +0.142 A, because the
outliers carry error that cancels); S25 (AMBER's top-75 is +1.103 A more expanded in Rg than the
pool, so a reject that keeps low-energy members will select expanded ones). Two things could
break the prior: the physics-set count is zero on the targets where rejection would have hurt,
and refilling keeps the averaging mechanism intact. I run it because that is the point of the
tournament, and I expect to close it.

Registered prior on the census (to be replaced by measurement in Addendum 1 before any RMSD is
read): at T = 1e4 about half of every top-75 is rejected (S20: 53.5% of the top-75 rebuilds are
above 1e4; S25: 58.6% of the pool); at 1e6 about a third; at 1e3 more than half; at 1e5
somewhere between. Refill at 1e4 will reach deep into the pool (rank 150 to 200).

## 6. Power

At n = 126 the SE of a paired built-chain difference is 0.01 to 0.04 A depending on how many
targets move (a threshold that rejects nothing on many targets produces exact ties). The MDE is
therefore 0.03 to 0.11 A. The production relaxation's +0.0207 A was measured with SE 0.0034, so
an operator that moves every target can be resolved at ~0.01 A; one that moves 30 targets cannot.
Every negative is reported with its MDE; below 0.7x MDE the word is "underpowered", not "null".

## 7. Cost and memory

Census: ~3 min CPU, < 0.5 GB (126 x 75 all-atom rebuilds through `sidechains.py`, no OpenMM).
Point cloud: ~5 min, < 0.6 GB. Built chain: 126 targets x up to 21 projections x ~4.4 s = up to
3.2 h CPU, < 0.7 GB, one process, checkpointed per target under `s26/results/ph_reject_chain_cells/`.
Agent-hours: 2 to write, 1 to analyse and report. Tags: CPU. Runs after `PHASE 0 SIGNED OFF`.

## 8. Replication

A positive result is re-run with the draw seeds changed (`stable_rng` salt "s26ph" to "s26ph-rep")
and the fold processing order reversed, and must land inside its own fold CI, or it is not a result.

---
## ADDENDUM 1 (2026-09-13 08:41) -- THE CENSUS, WRITTEN BEFORE ANY RMSD IS READ

`s26/results/ph_reject_census.json`, ledger L23. Measured what section 5 guessed: at 1e4 the
reject removes 40.1 of 75 (53.5%; prior "about half"), refill reaches rank 147 (median; prior
150 to 200), 11 targets lose the whole top-75 and 8 have no survivor in the whole pool (the prior
did not anticipate whole-pool emptiness). At 1e3 it is 54.5 of 75 with 25 empty sets and 20
empty pools; at 1e5 30.3 of 75; at 1e6 23.3 of 75.

Rules added before any endpoint is read: (1) an EMPTY retained set, S or R, falls back to the
anchor (do nothing), the only native-free deployable choice; the number of fallbacks per arm and
threshold, and the number of targets actually moved, are reported beside every MDE. (2) Every
contrast is shown on all 126 and on the moved subset. (3) The operator at the primary threshold
is now known to be a reject of the builder's side-chain placement (96.8% of condemned members
have a side-chain-involving closest contact; rho(e, min heavy-atom distance) = -0.74) and not of
the pool's backbones (2.6 of 75 have a backbone+CB contact below 2.0 A; S19 reproduced). That
does not change the falsifier; it changes what a positive result would mean, and it is stated
now so it cannot be discovered afterwards. Prior unchanged: null or harmful.

The census artefact was stamped by `ph_reject.py` at source sha ba028f7717d38725, before the
fallback lines (`R_empty`, `R_eff`) were added to `retained_sets`; the census computation does
not read those fields. The committed file (5dc7a3a6) is the one the endpoint runs will stamp.
