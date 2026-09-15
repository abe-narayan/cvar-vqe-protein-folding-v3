# S28 LANE C PRE-REGISTRATION -- FAIL18 DETECTION ON A NEW FEATURE CLASS, AND READOUTS THAT CONSUME RANKING INFORMATION WITHOUT AVERAGING IT AWAY

Written 2026-09-14 19:10 before any endpoint number of this lane was read. Binding contract:
`s27/S28_CONTRACT.md`. Brief: `s27/briefs/S28C.md`. Never edited after a result exists;
addenda are appended below with a timestamp.

## 0. What is held fixed

| | value |
|---|---|
| targets | the 126 dev targets, `sorted(pdb)` (`s25.phys_lib.targets`), all of them |
| folds | the 5 pinned folds (`s24.stats_lib.pinned_folds`) |
| pool | the shipped K = 500 BLOSUM pool per target (`s24.d_harness.Candidates.from_universe` through `s27.run_pool.channels_for`, which also asserts the S25 cache identity) |
| channels | `s27/cache/<pdb>.npz` (the S27 channels, native-free, regenerable by `s27/run_pool.py`) |
| production set | the shipped score's top-75 (`s27.run_pool.topm(ch["DIS"], 75, key)`, ties by the stable per-target key), whose uniform average reproduces the point-cloud anchor 3.048338 and whose projection reproduces the built-chain anchor 3.2126 (`s27/results/chain_rows.jsonl :: DIS`) |
| readout (reference) | the deployed uniform coordinate average in the retained set's medoid frame (`s24.d_harness.readout_uniform`) |
| projection | `s12.instrument.project` (ramah 0.3, multi-start), the production projection |
| bases | POINT CLOUD (intermediate, stated) and BUILT CHAIN (the reporting basis, decides) |
| statistics | `s24.stats_lib.compare`: paired per target, SE, MDE = 2.8016 x SE, iid CI beside the fold-clustered CI, W/L/ties, median, concentration null; decided on the fold-clustered CI |
| grids | any minimum over a grid is priced with `ST.best_of_k_within` (split-half transfer, k_eff); the nested leave-fold-out choice is what is reported |
| ties | never by array order (`s27.run_pool.topm` with the stable key; `ST.argmin_tied`) |
| ORACLE | `nat_ca`, `oracle_rr`, FAIL18 membership, any RMSD. FAIL18 membership is the training LABEL of Part 1 inside nested CV and is never read at inference. Every ORACLE quantity is labelled ORACLE in the sentence that uses it |

Reference numbers already on disk (no new endpoint read): point cloud DIS top-75 3.0483
(`s27/results/pool_summary.json :: dis_mean`), built chain DIS 3.2126 and DIS+DISTPOT 3.2100,
DIS+ENV 3.2361 (`s27/results/chain_rows.jsonl`), point-cloud FAIL18 strata DIS+DISTPOT -0.191
(SE 0.169, n = 18) / +0.010 (SE 0.034, n = 108), DIS+ENV -0.401 (SE 0.189) / +0.064 (SE 0.046)
(`s27/results/strata.json`).

## 1. Part 1: the FAIL18 detector

### 1.1 The new angle, stated against what closed
Seven router constructions (S22 L7: four; S23 L7: rg_z; S26 L110/L115: twelve m routers and six
s routers on six feature blocks) failed to predict a per-target quantity (m*, s*) native-free. None
of them (a) had FAIL18 membership as the label, and (b) saw the pool's own statistical-potential
distribution or its agreement with the distogram. Both are new here. S26 L110/L115's feature
blocks (retrieval-score entropy, posterior entropy, ESM contact statistics, principal-axis spread,
and the S22 set) are carried as COMPARISON blocks through the same harness, so "new feature class"
is a measured contrast, not a claim.

### 1.2 Label (ORACLE)
`s12.instrument.FAIL18`: the 18 targets whose shipped top-75 contains no member within 1.5 A of the
pool's best member (zero recall, S10-1). y = 1 for FAIL18, 0 otherwise. 18 positives / 108
negatives.

### 1.3 Features (native-free, `s27/s28_C_fail18.py :: features_one`), all read from the pool and
the S27 channel cache; none reads `nat_ca`, `oracle_rr` or any RMSD.

Block SP (the new class; for each X in {DISTPOT, ENV, CONTACT}, the three universe-fitted
statistical potentials of `s27/ham_lib.py`):
- `X_spread`: IQR of X over the 500, divided by the channel's natural size (number of pairs with
  |i-j| >= 3 for DISTPOT and CONTACT, n for ENV);
- `X_skew`: sample skewness of X over the 500;
- `X_rho_dis`: Spearman rank correlation of X with DIS over the 500;
- `X_overlap`: |top-75(X) intersect top-75(DIS)| / 75 (ties by the stable key);
- `X_top75_z`: mean of zrank(X) over DIS's top-75 (negative = DIS's set looks good under X);
- `X_top75_sd`: sd of zrank(X) over DIS's top-75;
plus three pool-consistency scalars: `cons_mean` (the pool's mean pairwise CA-RMSD, = mean of the
CONS channel), `cons_top75_gap` (mean CONS over DIS's top-75 minus `cons_mean`, in units of the
CONS sd), `pool_contact_ent` (mean binary entropy of the pool's 8 A contact frequencies,
|i-j| >= 3). 21 features.

Block CTRL (the classes previous routers used, recomputed here): `n`, `dg_ent_mean` (the
distogram's mean per-pair entropy, nats), `sim_ent` (softmax entropy of the pool's BLOSUM sims,
normalised), `sim_gap75`, `top75_spread` (mean pairwise CA-RMSD inside DIS's top-75). 5 features.

Blocks S26 (read from `s26/results/p_c4_features.json`, complete: true): `old_S22` (15
features, the S22 router set) and `s26_new_all` (25 features: sim_, dg_, con_, pax_). Comparison
blocks only.

Models run on: SP, CTRL, SP+CTRL, old_S22, s26_new_all, and every single feature.

### 1.4 Models, all nested leave-fold-out
- Ridge logistic regression (numpy Newton iterations, L2 penalty on standardised columns,
  balanced class weights 108/18), penalty alpha chosen by inner leave-one-fold-out over the four
  training folds maximising inner held-out AUROC, grid `np.logspace(-2, 4, 13)`; refit on the
  four folds; the held-out fold's decision values are recorded. Standardisation statistics are
  computed on the training rows only.
- One-threshold rule per feature: the direction (sign) is chosen on the training folds by AUROC;
  the held-out decision value is the signed feature. (For a single feature AUROC is
  threshold-free; the nested part is the sign.)

### 1.5 Metric and null
Primary: held-out AUROC of the pooled held-out decision values over the 126 targets (the
per-fold mean AUROC is reported beside it). Null: 500 label permutations across targets, the
identical nested procedure, per block. p_perm = share of null AUROCs >= observed.
Balanced accuracy at the 0.5 decision level (balanced weights) is reported beside AUROC.

### 1.6 The switched arm (built chain, the reporting basis; point cloud beside it)
For the SP+CTRL block (primary) and the SP block (secondary): predicted-FAIL (held-out
probability >= 0.5 under balanced weights) -> the DIS+DISTPOT top-75 (primary alternative;
DIS+ENV secondary); predicted-OK -> the DIS top-75. Per-target endpoints are read from
`s27/results/chain_rows.jsonl` (`rmsd_chain`, `rmsd_cloud`; the sets are deterministic
functions of the cached channels and the projection is deterministic), after a reproduction
check: the projection of the DIS and DIS+DISTPOT top-75 averages is re-run through
`s24.d_harness.readout_uniform` + `readout_projected` on the first 6 targets in sorted order and
must reproduce `rmsd_chain` to 1e-6 A. Paired against production (DIS built chain, 3.2126).
Controls: (i) the same switch applied to a random subset of the same size as the predicted-FAIL
set, 200 draws (mean and the distribution); (ii) the switched-arm endpoint under 200 label
permutations (the routed endpoint's own null, as S26 L110); (iii) the ORACLE switch (switch
exactly FAIL18), labelled ORACLE, which is the prize a perfect detector would collect and is
the ceiling of this part.

### 1.7 Falsifier and prior
The detector is a result only if BOTH hold: (F1) held-out AUROC of the SP or SP+CTRL block is
above the 95th percentile of its 500-draw label-permutation null AND above the best comparison
block's AUROC by more than the null's inter-quantile spread (so the class, not the harness,
carries it); (F2) the switched arm beats production on the built chain beyond its own MDE with
the fold-clustered CI excluding zero on 5/5 folds, beats control (i), and its permutation null
(ii) has the observed effect below its 5th percentile. Registered prior: NULL at F1 (AUROC
inside the permutation null, as the seven routers before it; S26 L1863 measured the FAIL18
Fisher test null for four strain signals). If F1 fires and F2 does not, the detector is an ORACLE
DIAGNOSTIC of where the error lives, not a lever. A positive is re-run with a second permutation
seed and the folds processed in reversed order before it is believed.

Expected effect if a perfect switch existed: on the point cloud the S27 strata give
18/126 x (-0.191) = -0.027 A (DISTPOT) and 18/126 x (-0.401) = -0.057 A (ENV); the built-chain
value is the ORACLE switch of 1.6(iii) and will be smaller or larger, unknown before the run.
The production-vs-production MDE at n = 126 on the built chain for a switch touching ~18
targets is expected near 0.03 to 0.05 A (S26 L110: 0.054 to 0.085 for full routers), so a
50%-accurate detector (prize ~0.015 to 0.03) is BELOW the MDE: Part 1 is expected to be
underpowered at F2 even if F1 fires, and that will be written as "underpowered", not "null".

## 2. Part 2: readouts that consume ranking information

### 2.1 The angle
S27 section 6 (`s27/results/mechanism.json`): a channel's correct ranking information is
anti-useful on the uniform top-75 average (Spearman +0.48 between a channel's partial rho with
the ORACLE candidate RMSD given DIS and the harm it does when added). CONS has the most (partial
rho +0.246) and hurts most (+0.286). Every readout below keeps the PRODUCTION retained set (DIS
top-75) and changes only how the set is consumed, so the selection lever (closed, S17/S21/S27)
is not re-opened; what is tested is the consumer. S12's medoid readout (-0.172 A on the S12
instrument, superseded on this basis) and S23 L5 (a geometric trim costs +0.142 A) are the
record; the new angle is a ranker-based (not geometry-based) trim, a medoid-plus-neighbours
hybrid, and a weighted average with an explicit diversity term, none of which has been run.

### 2.2 The readouts (`s27/s28_C_readout.py`), each a native-free function of the DIS top-75
(W75, the 75 x 75 pairwise CA-RMSD matrix P, and a ranker r over the 75)
(a) MEDOID+NEIGHBOURS(k): the member of the top-75 with the best (lowest) ranker value, plus
    its k nearest neighbours by CA-RMSD inside the top-75; uniform average of those k+1 members
    in their own medoid frame (`readout_uniform`). k in {10, 20, 40}. Identity: k = 74
    reproduces the production average exactly (same members, same function).
(b) TRIM(q): drop the ceil(q x 75) members with the worst ranker value; uniform average of the
    rest. q in {0.05, 0.10, 0.20} (4, 8, 15 dropped). Identity: q = 0.
(c) DIVERSITY-WEIGHTED(beta, gamma): all 75 members superposed on the production medoid
    (the same frame as production), then the weighted mean with
    w_i proportional to exp(-beta z_i) x (sum_j exp(-P_ij^2 / (2 sigma^2)))^(-gamma),
    z_i = zrank of the ranker within the 75, sigma = median of the off-diagonal P (native-free);
    the second factor is an inverse local density: members in a dense cluster share weight,
    isolated members keep theirs. beta in {0, 0.5, 1, 2}, gamma in {0, 0.5, 1}. Identity:
    (beta, gamma) = (0, 0) reproduces the production average to floating-point (asserted at
    1e-9 A).
Ranker: CONS (primary, the channel S27 names). DISTPOT as the secondary ranker for (b) only
(S27 section 10 item 1: the one non-redundant channel that does not collapse diversity).
Controls, matched in the operator's space: the same operator with the ranker rank-permuted
within the 75 by a stable per-target key (for (a): a random seed member; for (b): a random
trim of the same size; for (c): permuted z_i, gamma unchanged). Plus the random-75 null
(`s27/results/pool_rows.jsonl :: rand_mean`, point cloud) and production itself.
Every arm passes the NaN-poison test: `nat_ca` and `oracle_rr` replaced by NaN, the emitted
coordinates bit-identical.

### 2.3 Endpoints and budget
Point cloud for every cell (cheap). Built chain (3 s per target per arm) for the primary cells
and their matched controls: (a) k = 20 (the S27 section 10 item 3 hybrid) and its control;
(b) q = 0.10 with CONS and with DISTPOT, and their controls; (c) (1, 0.5) and its control,
(1, 0) and its control, (0, 1) (no ranker; its control is production). 10 built-chain arms.
Any further cell whose point-cloud effect against production is below -0.7 x its MDE is also
projected before it is discussed. Checkpointed per (arm, target) in
`s27/results/s28_C_readout_rows.jsonl`, resumable, under `s26/jobrun.py --agent S28C`.

### 2.4 Falsifier and prior
A readout is a result only if: its built-chain effect against production clears its own MDE
with the fold CI excluding zero on 5/5 folds; it beats its matched-random control by the same
standard (so the ranker, not the operator's geometry, carries it); its grid position is not an
order statistic (`ST.best_of_k_within` on the (target x cell) matrix of the family: split-half
transfer >= 25% of the oracle); and the point-cloud and built-chain signs agree. Registered
priors: (a) WORSE (it is a smaller m with a consensus seed: S22 L4's ladder has its minimum at
75 and S27 T5's DIS+CONS is worse at every m); (b) NULL (8 of 75 members change; the operator
is close to S23 L5's trim, which cost +0.142); (c) NULL to WORSE for beta > 0 (sharper weights
attack the variance reduction, S23 L8), NULL for gamma > 0 alone. The value of the part is to
close S27's "consumer of ranking information" question by measurement; a positive would be the
first readout on record that uses a ranker's in-pool information.

## 3. Rule-0 forks not taken
- A detector trained on the switched-arm gain (regression on the per-target difference) rather
  than FAIL18 membership: that is S22's router construction with a new label; not run, so the
  classification result is clean.
- Re-selecting the retained set with CONS before the readout: that is S27 H2/H10, closed.
- AMBER anywhere: no AMBER compute in this lane.

## 4. Memory and time
All classical. Channels are 500-vectors; the pool is (500, n, 3). Expected peak RSS below 0.6
GB per job (S27's `run_vqe_chain --chain` was the same shape). Probe one target first under
jobrun and quote the peak. Part 1: features 126 targets x ~1 s; nested logistic with 500
permutations x 5 blocks: minutes. Part 2: point cloud ~2 min; built chain ~65 min for the 10
primary arms.

## ADDENDUM 1 (2026-09-14 19:20, before any run)
Part 1 decision rule for the switched arm: PRIMARY as written (held-out probability >= 0.5 under
balanced class weights). SECONDARY, reported beside it: a prevalence-matched threshold, the
(1 - prevalence of the training folds) quantile of the training folds' decision values, so the
predicted-positive rate on the held-out fold is about the training prevalence (14%). Both are
nested; neither reads a held-out label. The single-feature rules use the same two thresholds on
the signed feature.
