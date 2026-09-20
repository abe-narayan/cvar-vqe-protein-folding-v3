# PREREG_S29_O -- THE ORACLE CEILING LADDER BY OPERATOR CLASS, AND THE TYPICALITY-AXIS PROBE

Lane O, Sprint 29. Written 2026-09-19 23:43 Pacific, before the first number. Brief:
`s29/briefs/S29O.md`. Contract: `s29/S29_CONTRACT.md` (inherits `s27/S28_CONTRACT.md`).
Code: `s29/s29_O_ladder.py`; tests `tests/test_s29_O.py`; results `s29/results/s29_O_*.json`,
`s29/results/s29_O_*_rows.jsonl`, `s29/results/s29_O_structs/<pdb>.npz`; findings
`s29/s29_O_FINDINGS.md`. Addenda are appended, never edited, once a result exists.

## 0. Labels, once and for all
EVERY NUMBER THIS LANE PRODUCES IS ORACLE (it reads `nat_ca` or `oracle_rr`) EXCEPT the
leave-fold-out step of rung 6, which is the lane's one DEPLOYABLE arm and its one paired
contrast. Every sentence carrying an ORACLE number says ORACLE in that sentence. Nothing here
tunes anything; nothing here is a result; the ladder is a ceiling table by operator class.

## 1. Findings engaged (contract rule 13)
The ladder attacks nothing: it MEASURES charter findings 8 (recognition) and 11 (common-mode
error) by operator class, and accepts 1 to 7, 9 and 10 as binding. Rung 6 is the coordinator's
H1 (`s29/STATE.md`), which attacks 3, 4, 5, 8, 11 and accepts 1, 2, 6, 7, 9; the registered
prior for it, from S24 L2/L3's geometry (`s24/LEDGER.md` lines 34 to 178), is that it FAILS:
cos near zero or negative, the step null. What is new relative to S16 (native-free
error-direction steering) and S24 L2 (mixtures with the blind source): the EXTRAPOLATION
(negative weight on the blind average, t > 0 below) was never run; S24 only ran interior
mixtures (t in [-1, 0]). What is new relative to S22 L4 (per-target m transfers) is nothing:
rung 2 re-reads it as the ceiling of the router class and prices its order statistic.

## 2. The instrument (reused, not re-implemented)
- Pool and production set: `s27.s28_A_amp.load_pool(pdb)` = `s27.run_pool.channels_for` (the
  shipped K=500 BLOSUM pool, `Candidates.from_universe`), the DIS channel, the production
  top-75 by `run_pool.topm` with the stable tie key. The DIS ORDER for rung 2 is
  `np.lexsort((key, DIS))` (the same tie key; ties never broken by array order).
- Point cloud: `s12.instrument.coordinate_average` (medoid frame, uniform mean) and
  `s12.instrument.ca_rmsd` (Kabsch). Built chain: `s12.instrument.project` (lam 0.3,
  multi-start, maxiter 300, the production projection), scored by `ca_rmsd`.
- Statistics: `s24.stats_lib.compare` / `fmt` / `best_of_k_within` / `split_half_transfer`;
  folds `s24.stats_lib.pinned_folds`; FAIL18 = `s12.instrument.FAIL18`.
- Blind clouds (rung 6): regenerated from S24's OWN stable RNGs, not new draws:
  LIB75 = `s24/biasalign.py` C1 (`SD.stable_rng("biasalign", pdb).choice(N, 75)`, the first
  draw, a uniform 75-window draw from the target's leakage-safe universe: no retrieval, no
  score, the "typical peptide of this length"); BPRIME = `s24/qmatch.py` B'
  (`SD.stable_rng("qmatch", pdb)`, 2,000 uniform windows with the BLOSUM pool excluded, the
  top-75 by the shipped score). Neither cloud is stored on disk (S24 stored only RMSDs and
  cosines), so both are REBUILT exactly as S24 built them and GATED: LIB75's RMSD must equal
  `s24/results/biasalign.json :: mix.m0_75` and its native-frame bias cosine
  `cos_A_C`; BPRIME's must equal `s24/results/qmatch.json :: rmsd_B` and `cos_A_B`; each to
  1e-6 on every target, or the run stops. The S12 blind pipeline (`s12/coord_null.py`,
  3.989 A) seeded with `hash(pdb)` and is not reproducible, so it is not used.

## 3. The ladder (all ORACLE; point cloud AND built chain; 126 targets)
Per target, ONE row of point-cloud values and ONE stored structure per chain item.
1. Best single member: of the top-75 (`oracle_rr[top].min()`) and of the pool
   (`oracle_rr.min()`). Anchors to reproduce exactly: S28-L1b 2.3062 / 1.7108. Chain items:
   `best1_top75`, `best1_pool` (the argmin window; ties reported, first taken, since tied
   `rr` at float32 come from duplicate windows).
2. Best m-subset of the DIS order, m = 1..500: `X_m` = the uniform medoid-frame average of the
   first m of the DIS order (m = 75 is production; m = 1 the shipped argmin). The full curve
   r_i(m) is stored; ORACLE per-target best m (the router ceiling, S22 L4's object at full
   pool), priced as an order statistic over 500 columns with `best_of_k_within` and the
   split-half transfer (prior: about 65% transfers, S22 L4); ORACLE best GLOBAL m (one m for
   all targets, the min of the mean curve) reported beside it. Chain item: `bestm` (the
   per-target ORACLE m).
3. Best basin average: for S in {top-75, top-500} and k in {2, 3, 4, 6, 8}: agglomerative
   AVERAGE-linkage clustering (`scipy.cluster.hierarchy.linkage(method="average")`,
   `fcluster(criterion="maxclust")`) on the pairwise CA-RMSD matrix
   (`s12.instrument.pairwise_rmsd`, which is `kabsch_rmsd_batch` row by row); every basin of
   size >= 2 is averaged (medoid frame, uniform mean); a singleton basin is rung 1's class
   and is not eligible; the ORACLE best basin per (S, k) is kept with its size. Chain items:
   `basin_<S>_k<k>` (10). The winning (S, k) and its basin size are reported per target and
   as counts; the per-target best over the 10 cells is an order statistic and is priced.
4. Best sparse convex combination with s in {2, 3, 5, 10, 20} members, over S in {top-75,
   top-500}: greedy forward support selection with a non-negative, sum-to-one least-squares
   refit (NNLS with the sum-to-one constraint as a weighted augmented row, weights
   renormalised exactly) against the native posed on the S28 frame (`s27.s28_A_amp.Frame`,
   all windows posed on the production medoid); the native's pose is re-solved on the current
   combination (frame alternation, 3 rounds per step) so the fit is the joint frame-and-weight
   problem of S10-5's `cf` convention. The greedy path is nested, so one pass gives every s.
   Reported RMSD is after Kabsch re-superposition (`ca_rmsd`). Chain items:
   `sparse_<S>_s<s>` (10).
5. Best convex combination over the top-75 and over 500: the same solver at full support
   (5 alternation rounds). Anchors: S10-5 raw 1.802 / 0.953, emitted 1.840 / 0.853 (Frank-
   Wolfe, certified). Gate: within 0.05 A of the raw anchors on the mean; if my solver is
   WORSE by more than 0.05 it is under-converged and the entry says so; if BETTER the S10-5
   convention differed and the entry says which. Chain items: `hull_top75`, `hull_pool`.
6. The typicality-axis probe (below). Chain items: `prod` (production re-projected in this
   job, the same-code-path comparator, S28-L26b's reason), `lfo_LIB75`, `lfo_BPRIME`.
7. For every rung: mean, median, FAIL18 mean, 108 mean, fraction under 2 A, on both bases,
   and the projection price (chain minus cloud) with its fold CI.
Projection budget: 28 structures per target (2 + 1 + 10 + 10 + 2 + 3), 3,528 projections at
3 to 6 s = 3 to 6 h in ONE resumable job (`--agent S29O --tag CPU --est-ram 0.5`), per-(item,
target) checkpoint (one jsonl line per projected structure), rung 6's three items for ALL
targets first, then rungs 1, 2, 5, then 4, then 3. One target is probed first and its peak
RSS quoted.

## 4. Rung 6: definitions, falsifiers, the one deployable arm
Per target, in the production cloud's own frame (C = the medoid-frame top-75 average):
- B_sup = the blind average Kabsch-superposed on C; N_sup = the native superposed on C
  (ORACLE); u = C - B_sup; v = N_sup - C. Both are projected onto the orthogonal complement
  of the 6-dim rigid-body tangent space at C (3 translations, 3 infinitesimal rotations about
  C's centroid), so each lives in the 3n - 6 shape degrees of freedom.
- cos(u, v) per target (ORACLE). Random reference: 16 Gaussian fields per target in the same
  3n - 6 space (`SD.stable_rng(pdb, "s29O_randfield", j, salt="s29O")`), cos(g_j, v): the
  signed mean (0 by symmetry) and the mean |cos| (analytic sqrt(2 / (pi (3n - 6))), about
  0.14 at n = 12) are both reported, measured beside analytic.
- The step: X(t) = C + t u for t on the grid T = {-1.0, -0.9, ..., +2.0} (31 values, fixed
  now; t = -1 is the blind average itself, t = 0 production, t > 0 the extrapolation S24
  never ran). r_i(t) = ca_rmsd(X_i(t), nat_i) (ORACLE). ORACLE global step t* = argmin over
  T of the mean curve (one scalar); ORACLE per-target step t_i* (an order statistic over 31,
  priced with `best_of_k_within` and split-half). The LEAVE-FOLD-OUT step: for each pinned
  fold f, t_f = argmin over T of the mean of r_i(t) over i NOT in f, applied to every i in f.
  This arm is DEPLOYABLE (a parameter fitted on other targets' natives, as any trained
  model); its NaN-poison test is in `tests/test_s29_O.py` (the held-out target's native set
  to NaN leaves its emitted X bit-identical).
- PRIMARY blind definition: LIB75 (the sequence-blind object H1 is about). SECONDARY: BPRIME
  (the brief's named fallback; S24 L3 measured its bias cosine with production at 0.943, so
  u is expected small and the step expected null; run because the brief names it, counted in
  the multiplicity).
- The one paired contrast: `lfo_LIB75` vs `prod`, point cloud first (free), then the BUILT
  CHAIN (`s12.instrument.project` on both sides in the same job), `ST.compare` with the
  pinned folds, `ST.fmt` verbatim, FAIL18 / 108 strata, and for any stratum claim the
  random-18 null (20,000 random 18-subsets of the 126 without replacement, one-sided p in the
  observed direction, S28-L40's construction).
Falsifiers (H1 dies if EITHER fires; the registered prior is that BOTH fire):
- F6a (the cosine): the fold-clustered 95% CI of the mean cos(u, v) over 126 (LIB75) includes
  the random reference's signed mean (0), OR the mean cos is at or below the measured mean
  |cos| of the random fields (about 0.14) with the fold CI including that value.
- F6b (the step): `lfo_LIB75` vs `prod` on the point cloud does not clear 0.7x its own MDE
  with the fold CI excluding zero in the improving direction. If it clears on the cloud, the
  built-chain contrast decides; a chain effect below 0.7x MDE is NOT A RESULT.
- Regime clause: a FAIL18-only sign (conditioning was harmful there, S12 coord_null) is a
  regime claim only if the random-18 null gives p < 0.05; otherwise it is the whole-set effect
  read on 18 targets and is not quoted as a stratum result.
The registered prior (S24 L2: the blind source's bias cosine with production 0.647, q 1.231;
S24 L3: the score-selected retrieval-free source at cos 0.943, mixtures on a straight line):
u has a component ALONG the shared error and a large component that is minus the blind
cloud's orthogonal error, so cos(u, v) is near zero or negative and every t > 0 is worse than
production. STATE.md's prior for H1 (cos 0.2 to 0.4 on the 108) is the alternative. Either way
is stated as measured.

## 5. Reproduction gates (the run stops if one fails)
- Production point cloud: mean 3.048338 over 126 and per-target max |dev| <= 1e-9 against
  `s27/results/s28_A_structs/<pdb>.npz :: prod` (S28-L1b's frame).
- Rung 1 anchors: 2.3062 / 1.7108 (S28-L1b) to 1e-4.
- LIB75 and BPRIME gates of section 2 (S24's own artefacts) to 1e-6 per target.
- Rung 2 at m = 75 equals production on every target to 1e-9.
- Chain anchor (NOT a gate, a documented comparison): `prod` re-projected here against
  S28-L26b's 3.2071 (`s27/results/s28_A_chain_rows.jsonl :: prod`) and S27's 3.2126
  (`s27/results/chain_rows.jsonl :: DIS`); the branch-flip floor (S28-L18/L26b: 12/126
  targets differ by > 0.02 A, 3 by > 0.1) is expected and reported, never a contrast.
- Unit tests (`tests/test_s29_O.py`, synthetic pools, no target, no cache): rung 2 at m = k
  and rung 3 at k = 1 and rung 4/5 at uniform-optimal data reduce to the uniform average;
  the sparse solver at s = k reproduces the full-support hull solve; the rigid-body removal
  annihilates every rigid displacement and leaves shape displacements unchanged; cos of a
  vector with itself is 1 and with its rigid image is 0 after removal; the LFO step's
  NaN-poison; the random-18 null on a synthetic vector.

## 6. Multiplicity
Endpoint (deployable) comparisons this lane will run: 4 = {LIB75, BPRIME} x {point cloud,
built chain}; the PRIMARY is `lfo_LIB75` vs `prod` on the built chain; the point-cloud
contrast is its gate, the BPRIME pair is secondary. Stratum tests: 2 (FAIL18 for each
definition on the chain), each against the random-18 null. Every ORACLE contrast in the
ladder is a diagnostic and is not an endpoint comparison. Reported per entry for the
coordinator's sprint total.

## 7. What this lane will not do
No tuning of any deployable parameter on a native; no chain of any structure not listed in
section 3; no second seed for a null; no re-projection of S28's stored structures (their
chain rows are reused as anchors, never recomputed); no benchmark60; no AMBER.
