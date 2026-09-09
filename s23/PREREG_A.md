# PREREG_A — WORKSTREAM A, ENSEMBLE AGGREGATION (Sprint 23)

Never edited after data is read. Addenda are dated and appended, never in-place edits.

**Instrument.** `s12/instrument.py`: 126 cluster-disjoint dev targets, 5 pinned folds, K=500
BLOSUM pool, shipped Bayes-risk distogram score, full-chain Cα-RMSD after Kabsch.
**Comparator, all experiments: `avg_75` = `I.coordinate_average(W[order[:75]])`, the incumbent,
3.048 Å.** Point-cloud basis throughout — every arm here is a point cloud on the pool's own window
coordinates, compared only against other point clouds, per BRIEF hard rule.

**Statistics, every table.** Paired bootstrap (4000 resamples) iid CI, fold-clustered CI (resample
whole folds), SE, MDE = 2.8016×SE, effect reported as a multiple of its own MDE, W/L, worst-target
delta, mean virtual Cα–Cα bond length per arm (physical 3.805 Å, incumbent already 2.961 — contracted
geometry is a cost, not a free win). Nested CV (leave-one-fold-out, matching `s22/routercv.py`'s
convention) for every fitted parameter: select on the 4 training folds, apply to the held-out fold,
concatenate the 5 held-out columns into the achieved arm. In-sample ("oracle" over the grid, all 126
targets) is reported beside it and labelled a **leakage ceiling, not achievable**.

---

## EXPERIMENT 1 — CLUSTER-RESTRICTED AVERAGING (priority 1)

**Mechanism.** If the top-m spans two structural basins, the coordinate average sits between both
and resembles neither. `medoid_75` (pick one, 3.282) and `avg_75` (average all, 3.048) are both
measured; "cluster the top-m, average within the dominant cluster" is not.

**Construction.** Primary at m=75 (matches the incumbent's own width — the cleanest matched
comparison). Pairwise Cα-RMSD among the top-m (already the instrument's own metric). Average-linkage
hierarchical clustering (`scipy.cluster.hierarchy.linkage(method="average")`) on that distance
matrix, cut to k clusters via `fcluster(criterion="maxclust")`, k ∈ {2,3,4,5}. Within the chosen
cluster, `I.coordinate_average` (superpose on the cluster's own medoid, mean).

**Cluster-choice rules**, all computed native-free from the pool's own score and geometry:
- `size` — largest cluster by member count.
- `score` — cluster with the lowest mean shipped score.
- `combined` — lowest (size-rank + score-rank) sum, ranks 0-indexed best-first on each axis
  independently (size descending, score ascending). Declared before any RMSD is read.
- `random` — uniformly random cluster choice (control: does the CHOICE rule matter, or does any
  restriction to one cluster help regardless of which).
- `random_same_size_partition` — at k=3 only: shuffle pool members into groups matching the REAL
  clustering's size distribution (same sizes, random membership), apply `combined`'s pick logic.
  Isolates whether genuine geometric similarity is doing the work or whether any same-sized subset
  restriction would do as well.
- `oracle` — CEILING: the cluster (of the k found) whose coordinate average has the lowest RMSD to
  native. Reported beside every achievable arm, never itself an achievable arm.

**Robustness.** The k=3, `combined` cell is repeated at m ∈ {150, 250} to test whether the effect
strengthens as the pool widens (more room for a second basin to appear).

**Fork alternative named and not taken:** k-medoids/PAM on the same distance matrix is the natural
alternative to average-linkage hierarchical clustering and is NOT run here (budget); if this
experiment shows a live effect, k-medoids is the first replication to run before promotion.

**Six-axis operator fork (Rule 0), declared before any RMSD is read:**
| axis | DECLARED | NOT TAKEN, and why |
|---|---|---|
| functional | shipped Bayes-risk distogram score, for the top-m filter stage (identical to the incumbent's own filter) | squared-distance functional — S21 D7 found it interchangeable at the argmin on pool-like structures, not tested here for budget |
| basis | pool's own window coordinates (point cloud), matching the incumbent's own basis | rebuilt/built-chain manifold — would confound the point-cloud/built-chain comparison the BRIEF forbids mixing |
| readout | average-linkage hierarchical clustering, cut to k, coordinate-average within the chosen cluster | k-medoids (named above, not run); DBSCAN/density clustering (no natural k, harder to make native-free) |
| normalisation | none; raw Å Cα-RMSD; cluster ranks (size, score) computed on raw counts/scores | per-target z-scoring of cluster sizes or scores before ranking |
| null | k=1 clustering degenerates exactly to `avg_75`/`avg_m` (the existing incumbent/ladder); `random_same_size_partition` is the second null, isolating geometry-informativeness from restriction-per-se | a fully random single-structure baseline — already priced elsewhere (`pool_argmin`, `medoid_all`) |
| THE LABEL | continuous point-cloud Cα-RMSD | any binarised "did clustering separate correctly" label, whose threshold would be chosen after seeing data |

**Hypothesis.** At least one achievable (non-oracle, non-random) cluster-choice rule beats `avg_75`
by more than its own MDE with a CI excluding zero, at some (m,k).

**Falsifier.** No achievable rule beats `avg_75` past its own MDE with a CI excluding zero at ANY
(m,k) tested. Recorded as a full refutation of the reachable mechanism — separately, if `oracle`
ALSO fails to beat `avg_75` by more than half its own MDE, the mechanism itself (not just its
reachability) is refuted: there is no basin-splitting headroom to capture on this pool.

**Promotion.** None from this file without an adversarial replication (k-medoids check, at minimum).

---

## EXPERIMENT 2 — ENSEMBLE-WIDTH SWEEP, DONE PROPERLY (priority 2)

**Scope, explicit.** Per-target m selection is CLOSED (S22 L4/L5/L7/L10: real, transferable,
0.243 Å, but unreachable by any native-free router this project's sample size supports). This
experiment does NOT build a router. It re-runs the GLOBAL m-ladder — m ∈ {20,30,50,75,100,150,250,
500} — under two aggregation operators to see whether the global optimum moves:
1. `avg_m` — plain coordinate average (extends the existing {500,150,75,20,5,1} ladder to the full
   8-point grid on this instrument's own pool for self-consistency).
2. `clust_m` — the k=3, `combined`-rule cluster-restricted average from Experiment 1, applied at
   every m in the grid. k=3 and the `combined` rule are FIXED here (not re-optimised per m), declared
   before running, precisely to avoid a post-hoc "pick whichever cluster config looks best at each m"
   leak.

**Reported per m, both operators:** mean, median, WORST-target RMSD, variance, failure rate
(pre-registered threshold: fraction of targets with RMSD > 5.0 Å — a fixed absolute bar chosen before
any RMSD is read, not a data-dependent percentile).

**Fork axes**, abbreviated since they inherit Experiment 1's readout declaration:
functional/basis/normalisation/THE LABEL identical to Experiment 1. **null**: `avg_m` at each m is
itself the null for `clust_m` at that m (does clustering help AT THIS WIDTH, holding width fixed).
**readout**: DECLARED the width grid above; NOT TAKEN a continuous/interpolated m.

**Hypothesis.** The optimal m under `clust_m` differs from 75 (the confirmed `avg_m` optimum, S22
L5), because clustering removes the basin-mixing cost that penalises wide pools under plain
averaging, so a wider m may become optimal once mixing is controlled.

**Falsifier.** `clust_m`'s optimal m (by mean RMSD) is 75, matching `avg_m`'s, at every m the
`clust_m` mean is within its own MDE of `avg_m`'s mean — i.e., clustering does not change the shape
of the curve, only (at best) its level.

**Promotion.** Descriptive; no new arm is promoted from this experiment alone.

---

## EXPERIMENT 3 — WEIGHTED AVERAGING (priority 3)

**Construction, at m=75 (matches the incumbent width exactly — isolates the weighting decision from
the width decision).** Three weight families, each with ONE free parameter, declared BEFORE any RMSD
is read:

1. **Score-weighted (softmax).** `w_i ∝ exp(-score_i / τ)`, `τ = c · sd(score in top-75)` (native-free
   normalisation — score scale varies by target/pool). Grid: `c ∈ {0.1, 0.25, 0.5, 1, 2, 5, 1e6}`
   (`c=1e6` ≈ uniform = the incumbent, the null end of the grid).
2. **Rank-weighted.** `w_i ∝ (75 − rank_i)^p`, rank 0-indexed best-first. Grid:
   `p ∈ {0, 0.5, 1, 2, 4, 8}` (`p=0` = uniform = the incumbent).
3. **Distance-to-medoid-weighted.** `d_i` = pairwise Cα-RMSD from member i to the top-75's own
   medoid (already computed for `avg_75`). `w_i ∝ exp(-d_i / σ)`, `σ = c · sd(d)`. Grid:
   `c ∈ {0.1, 0.25, 0.5, 1, 2, 5, 1e6}` (`c=1e6` ≈ uniform = the incumbent).

All three superpose every member onto the SAME frame `avg_75` already uses (the top-75's own medoid)
before taking the weighted mean — the frame is held fixed so only the weight decision varies.

**Fitting.** Nested CV, leave-one-fold-out (`s22/routercv.py` convention): for held-out fold f,
select the grid value minimising mean RMSD on the other 4 folds; apply to fold f; concatenate the 5
held-out columns as the ACHIEVED arm. In-sample-oracle grid value (fit and scored on all 126) is
reported beside it as a leakage ceiling, never an achievable number.

**Six-axis operator fork:**
| axis | DECLARED | NOT TAKEN |
|---|---|---|
| functional | shipped Bayes-risk score (for the weighting itself, where applicable) | squared-distance functional |
| basis | pool's own window coordinates | rebuilt manifold |
| readout | weighted mean, frame = `avg_75`'s own medoid (frame held fixed across all weight arms) | weighted mean with a per-arm-optimal frame (would confound weight choice with frame choice) |
| normalisation | per-target sd-based bandwidth (`c · sd(score)`, `c · sd(d)`) — native-free, scale-matched to each target's own pool | a single global absolute τ/σ in Å or score-units, which would not transfer across targets with different score/RMSD scales |
| null | `c=1e6` / `p=0` (uniform weights) inside each family — algebraically `avg_75` itself | a size-matched random-weight control |
| THE LABEL | continuous Cα-RMSD | any win/loss binarisation |

**Hypothesis.** At least one family's nested-CV achieved arm beats `avg_75` (uniform) past its own
MDE with a CI excluding zero.

**Falsifier.** None of the three families' achieved (held-out) arms beats `avg_75` past its own MDE
with a CI excluding zero — weighted averaging is refuted as a reachable lever, independent of whether
the in-sample oracle shows headroom (which would then be reported as an unreachable ceiling, in the
same register as H1's scale correction and S22's per-target m).

**Promotion.** None from this file without an ablation across the frame choice.

---

## LOGGING (BRIEF §5)

experiment id `agentA_exp{1,2,3}` · git commit at run time · seed via `s15.seed.stable_rng` ·
target ids = the pinned 126 · config = grids above, frozen at this file's timestamp · Hamiltonian =
none (structural aggregation only, downstream of the shipped distogram score) · ensemble size = m
per arm · post-processing = none (no repair/relax in this workstream; point-cloud only) · mean/
median/worst RMSD, CI, per-target delta vs `avg_75` reported for every arm.

Dated 2026-09-08, before any experiment script is run.
