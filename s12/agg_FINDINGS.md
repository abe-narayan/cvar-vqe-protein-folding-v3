# s12 SET-AGGREGATION — findings

**Question.** Given the candidate set the pipeline already has, how much accuracy is
recoverable by a better TERMINAL OPERATOR — and is learned aggregation over the full
candidate-vs-objective deviation structure the missing piece that scalar-feature rankers
could not be?

## 0. Instrument validation (run first)

`python -m s12.instrument` reproduced exactly:

| quantity | expected | observed |
|---|---|---|
| shipped distogram argmin | 3.4540 | 3.4540005 |
| K=500 pool best (oracle) | 1.7108 | 1.7108244 |
| top-75 best (oracle) | 2.3062 | 2.3061526 |
| synthesis (avg -> project, lam=0) | 3.2041 | 3.2040762 |
| zero-recall targets | 18 | 18 (== FAIL18) |

I additionally re-derived the production synthesis path from my own cache
(`s12/agg_common.py`): medoid superposition -> uniform average of the shipped top-75 ->
`I.project(...)["fit_ca"]` reproduces the shipped `fit_ca` to **< 0.006 A** on every target
checked (typically 2e-4). So the emitted-value harness is the production one.

Code: `s12/agg_common.py` (cache + operators), `s12/agg_ladder.py` (experiment 1),
`s12/agg_meta.py` (pool metadata).

---

## 1. The operator ladder (RAW, i.e. before projection) — `s12/agg_ladder.py`

**Hypothesis.** The terminal operator is worth a lot (C1), so a better native-free operator
on the candidate set the pipeline already has should recover accuracy.

44 operator arms, all 126 targets, RAW CA-RMSD of the emitted point cloud to the native.
Baseline for the paired column is `avg75` = the production operator (medoid superposition +
uniform coordinate average of the shipped top-75).

### 1a. Native-free arms (deployable)

| arm | mean | FAIL18 | other108 | frac<2 | d vs avg75 | CI95 | W/L | drop10 |
|---|---|---|---|---|---|---|---|---|
| wscore_T2.0_75 | 3.0445 | 5.822 | 2.582 | 0.294 | -0.0038 | [-0.027,+0.019] | 65/61 | +0.021 |
| wrank_r20_75 | 3.0459 | 5.821 | 2.583 | 0.294 | -0.0025 | [-0.021,+0.015] | 63/63 | +0.017 |
| **avg75 (production)** | **3.0483** | 5.832 | 2.584 | 0.294 | 0 | — | — | — |
| iter_avg75 | 3.0523 | 5.821 | 2.591 | 0.286 | +0.0040 | [-0.008,+0.015] | 51/75 | +0.016 |
| trim10_75 | 3.0632 | 5.862 | 2.597 | 0.286 | +0.0148 | [-0.001,+0.030] | 58/68 | +0.029 |
| geomed75 | 3.0682 | 5.839 | 2.606 | 0.278 | +0.0199 | [-0.003,+0.043] | 47/79 | +0.042 |
| avg50_pool | 3.0676 | 5.827 | 2.608 | 0.286 | +0.0193 | [-0.010,+0.048] | 53/73 | +0.050 |
| avg150_pool | 3.0716 | 5.718 | 2.631 | 0.286 | +0.0233 | [-0.033,+0.080] | 66/60 | +0.079 |
| trim20/35_75 | 3.069/3.079 | | | | +0.021/+0.031 | | 53/73, 54/72 | |
| bigmode3_75 (largest structural mode only) | 3.0951 | 5.895 | 2.629 | 0.325 | +0.0468 | [-0.018,+0.110] | 53/73 | +0.114 |
| wcons_T1.0_75 (consensus weighting) | 3.0848 | 5.903 | 2.615 | 0.278 | +0.0364 | [+0.002,+0.069] | 47/79 | +0.073 |
| avg300_pool | 3.1715 | 5.348 | 2.809 | 0.278 | +0.1232 | [+0.009,+0.243] | 53/73 | +0.239 |
| medoid75 | 3.2822 | 6.100 | 2.813 | 0.278 | +0.2339 | [+0.162,+0.305] | 34/92 | +0.299 |
| avg500_pool | 3.3961 | 5.220 | 3.092 | 0.238 | +0.3478 | [+0.188,+0.517] | 49/77 | +0.484 |
| argmin75 (shipped baseline) | 3.4540 | 6.008 | 3.028 | 0.214 | +0.4057 | [+0.279,+0.538] | 31/95 | +0.516 |

**Result: the production operator is already at the native-free optimum.** Of 30 deployable
operator arms — trimmed means at four trim fractions, the geometric median, iterated
superposition-averaging, four score-weighting temperatures, two rank weightings, two
consensus weightings, two joint weightings, cluster-then-average at k=2 and k=3, and a
six-point cardinality ladder over the K=500 pool — **not one beats the uniform average of
the shipped top-75 by more than 0.004 A, and no CI excludes zero on the good side.** The
best two (`wscore_T2.0`, `wrank_r20`) are near-uniform weightings, i.e. perturbations of
`avg75` toward `avg75`; their advantage reverses when the top 10 targets are dropped.

The *operator form* is therefore not the missing piece. Argmin -> average is worth
0.406 A raw (C1 is confirmed: the old argmin cap was an artefact), and that 0.406 A is
**already banked in the shipped pipeline**. There is no second helping.

### 1b. Oracle arms — where the headroom actually is

| arm | mean raw | FAIL18 | other108 | frac<2 |
|---|---|---|---|---|
| ORC_simplex500 (affine/convex fit of the whole K=500, oracle) | 1.374 | 2.013 | 1.267 | 0.873 |
| ORC_top25avg_pool (uniform avg of the 25 truly-best pool members) | 1.608 | 2.329 | 1.488 | 0.738 |
| ORC_bestpool (single best pool member) | 1.711 | 2.284 | 1.615 | 0.579 |
| **ORC_simplex75 (optimal convex weights over the shipped top-75)** | **1.936** | 4.008* | 1.663* | 0.619 |
| ORC_wrmsd_T0.25_75 (softmax on true RMSD) | 2.156 | 4.571 | 1.754 | 0.516 |
| ORC_best75 (single best of the top-75) | 2.306 | 4.677 | 1.911 | 0.437 |
| ORC_top25avg_75 (uniform avg of the 25 truly-best of the 75) | 2.437 | 5.049 | 2.001 | 0.460 |
| avg75 (deployable) | 3.048 | 5.832 | 2.584 | 0.294 |

*(the FAIL18/other108 split for ORC_simplex75 is from the 400-iteration run, mean 1.998;
the 800-iteration value 1.936 is the converged one. C1's 1.644 for oracle top-25 pool
averaging and the record's 2.044 for oracle top-75 re-weighting are both reproduced within
the projection/convergence tolerance.)*

**Decomposition of the 1.11 A oracle headroom on the top-75** (raw, converged simplex):

| step | value | share |
|---|---|---|
| avg75 (uniform, all 75) | 3.048 | — |
| uniform average over the ORACLE SIMPLEX SUPPORT (mean support 40 of 75) | 2.441 | 0.607 A = **55%** |
| optimally-weighted average over the same set | 1.936 | 0.505 A = 45% |

So **subset SELECTION and continuous WEIGHTING each carry about half** of the recoverable
accuracy. Selection alone (choosing ~40 of 75 and averaging uniformly) is worth 0.61 A.
This is what makes the combinatorial formulation of section 4 worth posing.

### 1c. Why no deployable weighting finds it

Correlations of the oracle optimal weights with the deployable signals, averaged over 126:

| signal | corr with oracle simplex weight |
|---|---|
| true CA-RMSD of the candidate (ORACLE) | **-0.579** |
| shipped distogram Bayes-risk score | **+0.033** |
| mean pairwise RMSD to the set (consensus) | +0.232 |

The optimal weights are strongly (but far from perfectly) anti-correlated with true RMSD, and
**essentially orthogonal to the shipped score (+0.03)**. Note the sign on consensus: the
optimal weighting puts *more* weight on the more outlying members, the opposite of what
consensus/medoid heuristics do — consistent with the record's "consensus is outlier
avoidance, not a nativeness signal". The information needed to weight is present in the set
(it exists, -0.58 against the label) and is not in the objective the pipeline scores with.


---

## 2. MULTI-HYPOTHESIS OUTPUT -- `s12/agg_modes.py` (`s12/results/agg_modes.json`)

**Hypothesis.** Emitting 2-5 structural modes instead of one average opens oracle headroom;
with the set-level information a native-free rule may reach some of it. Also: a MODE-AWARE
average (equal weight per mode, de-biasing the average for mode population) may beat the
global average even when the mode cannot be picked.

Clusterings: agglomerative-average, Ward, k-medoids on the 75x75 CA-RMSD matrix, k = 2..5.
Rules (13 per clustering): largest, tightest, best mean score, best min score, best risk of
the mode average, best risk of the mode mean pair distances, contains-the-argmin, closest
to the global average, best BLOSUM sim, most peptide-derived (org), closest rg, best mean
rank, smallest member spread; plus two soft mixtures and the mode-aware average.
**Ties are resolved by averaging the outcome over the whole tied argmin set** (the
`consensus-is-the-only-in-band-discriminator` trap: `np.argmin` on a tied signal reads the
pool sort order and invents a winner).

| clustering | ORACLE best-of-k | ORACLE worst-of-k | mode-aware avg | best native-free rule |
|---|---|---|---|---|
| average, k=2 | 2.791 | 3.795 | 3.115 | largest 3.072 |
| average, k=3 | **2.681** | 4.102 | 3.153 | largest 3.088 |
| average, k=5 | 2.584 | 4.341 | 3.119 | largest 3.130 |
| ward, k=3 | 2.693 | 3.879 | 3.077 | nearglobal 3.119 |
| ward, k=5 | 2.591 | 4.081 | **3.061** | nearglobal 3.169 |
| kmed, k=3 | 2.674 | 3.931 | 3.108 | largest 3.098 |
| -- global average (avg75) | | | **3.048** | |

Best 4 of the **192** native-free mode arms, paired against avg75:

| arm | mean | FAIL18 | other108 | d | CI95 | W/L | drop10 |
|---|---|---|---|---|---|---|---|
| ward5_modeavg | 3.0614 | 5.841 | 2.598 | +0.0131 | [-0.022,+0.045] | 62/64 | +0.048 |
| average2_pick_largest | 3.0717 | 5.899 | 2.601 | +0.0234 | [-0.021,+0.067] | 50/76 | +0.069 |
| ward2_modeavg | 3.0761 | 5.846 | 2.614 | +0.0278 | [-0.016,+0.072] | 59/67 | +0.071 |
| ward3_modeavg | 3.0774 | 5.838 | 2.617 | +0.0290 | [-0.012,+0.071] | 66/60 | +0.067 |
| (worst of 192) average5_pick_tightest | 3.8237 | | | +0.775 | | | |

**Result: the multi-hypothesis channel stays closed.** Best-of-3 oracle headroom is 0.367 A
(the record ~0.43 is reproduced within clustering choice); **all 192 native-free arms are
worse than the single global average**, the best by +0.013 A, and that is the maximum over
192 comparisons, so the honest best is if anything worse. The MODE-AWARE average does not
help either (+0.013 at best): the population imbalance of the modes (mean sizes 57/13/5 at
k=3) is not a bias worth correcting. Picking by "tightest" is actively catastrophic
(+0.775) -- the outlier-avoidance failure mode again.

---

## 3. SET SELECTION AS A COMBINATORIAL PROBLEM -- `s12/agg_subset.py`

### 3a. Exact objective on the full 75 (`s12/results/agg_subset_exact.json`)

    F(S) = Risk(S) + mu * Cons(S)
    Risk_dist(S)   = mean_p risk_p( mean_{c in S} d_{p,c} )      (cardinality-neutral)
    Risk_struct(S) = mean_p risk_p( d_p( mean_{c in S} A_c ) )   (shipped score of the
                                                                  averaged structure)
    Cons(S)        = mean_{c != c'} P_{c c'}

`risk` is the SHIPPED leave-fold-out distogram Bayes risk (bit-identical to
`I.shipped_score`). Solvers: forward greedy, greedy + 1-opt local search (drop/add/swap to
convergence), simulated annealing (3000 steps, 2 seeds), and greedy at fixed k.

| arm | mean | FAIL18 | other108 | d vs avg75 | CI95 | W/L | drop10 |
|---|---|---|---|---|---|---|---|
| **ORC_greedy** (greedy on TRUE RMSD, mean k=3.8) | **2.0615** | 4.572 | 1.643 | -0.9869 | [-1.104,-0.877] | 126/0 | -0.865 |
| avg75 (take everything) | 3.0483 | 5.832 | 2.584 | 0 | -- | -- | -- |
| greedy fixed k=40 | 3.0521 | 5.832 | 2.589 | +0.0037 | [-0.033,+0.041] | 64/62 | +0.039 |
| SA, dist, mu=0.03 | 3.0633 | 5.783 | 2.610 | +0.0150 | [-0.030,+0.061] | 66/60 | +0.060 |
| greedy fixed k=25 | 3.0774 | 5.804 | 2.623 | +0.0291 | [-0.021,+0.080] | 57/69 | +0.076 |
| greedy fixed k=10 | 3.1275 | 5.787 | 2.684 | +0.0791 | [+0.009,+0.154] | 55/71 | +0.144 |
| greedy+LS, dist, mu=0 (mean k=5.5) | 3.1601 | 5.839 | 2.714 | +0.1118 | [+0.028,+0.198] | 52/74 | +0.186 |
| greedy+LS, struct, mu=0 (mean k=4.3) | 3.2840 | 5.941 | 2.841 | +0.2357 | [+0.131,+0.350] | 42/84 | +0.324 |
| argmin (single best-scoring) | 3.4540 | 6.008 | 3.028 | +0.4057 | [+0.279,+0.538] | 31/95 | +0.516 |

**Result: the combinatorial problem is well-posed and has 0.99 A of ORACLE headroom (2.062,
with a mean of only 3.8 selected members) -- and the native-free objective captures none of
it.** Every solved arm is >= avg75, and the harder the selection (smaller k), the worse the
answer, monotonically: k=40 ~ avg75, k=10 +0.08, unconstrained greedy (k~5) +0.11, `struct`
variant +0.24, argmin +0.41. A clean instance of "optimise the objective harder, get
further from the native".

### 3b. EXACT optima: the QUBO instances (`s12/results/agg_subset_instances.json`)

Hypothesis set: the 16 medoids of a 16-way agglomerative clustering of the shipped top-75
(so the 16 candidates span the structural diversity of the set rather than duplicating it).
16 binary variables, so all 2^16 = 65536 subsets are enumerated exactly, for **all 126
targets**.

Objective dumped (quadratic, ready for an Ising / CVaR-VQE encoding):

    H(x) = (1/k^2) x' M x  -  (2/k) b' x  +  const ,   k = sum(x) >= 2 ,  x in {0,1}^16
    M     = Q + mu * P   (mu = 0.01)
    Q_cc' = (1/npairs) sum_p d_{p,c} d_{p,c'} / sd_p^2
    b_c   = (1/npairs) sum_p e_p d_{p,c} / sd_p^2
    const = (1/npairs) sum_p e_p^2 / sd_p^2

i.e. the sd-weighted squared deviation of the mean pair distances of the subset from the
distogram predicted distances, plus a consistency term (P = candidate-vs-candidate CA-RMSD).
Emission = uniform coordinate average of the selected medoid-superposed candidates.
Each dumped instance carries `M`, `b`, `const`, `hypo_idx_into_top75`, the certified
exhaustive optimum, the greedy solution, and the oracle best subset.

| arm | mean | FAIL18 | other108 | frac<2 | d vs uniform-16 | CI95 | W/L |
|---|---|---|---|---|---|---|---|
| avg75 (all 75, no selection) | 3.0483 | 5.832 | 2.584 | 0.294 | -0.0241 | [-0.071,+0.024] | 66/60 |
| uniform over all 16 hypotheses | 3.0724 | 5.685 | 2.637 | 0.270 | 0 | -- | -- |
| greedy on H | 3.1768 | 5.806 | 2.739 | 0.246 | +0.1044 | [+0.045,+0.167] | 52/74 |
| **EXACT optimum of H (2^16 enumeration)** | **3.1990** | 5.815 | 2.763 | 0.246 | +0.1266 | [+0.067,+0.186] | 46/80 |
| exact optimum of the shipped-risk F | 3.1911 | 5.837 | 2.750 | 0.254 | +0.1186 | [+0.061,+0.176] | 43/83 |
| ORACLE best subset of the 16 (mean k=2.7) | 2.2810 | 4.749 | 1.870 | 0.476 | -0.7914 | [-0.881,-0.708] | 126/0 |
| ORACLE best single member of the 16 | 2.4968 | 4.816 | 2.110 | 0.365 | -0.5756 | [-0.683,-0.471] | 113/13 |

**Solution-quality audit for the quantum arm:**

| quantity | value |
|---|---|
| greedy objective gap to the exact optimum | 0.0319 mean; greedy is EXACTLY optimal on **40.5%** of instances |
| accuracy cost of closing that gap (greedy -> exact) | **+0.0222 A WORSE** |
| corr(H(S), true RMSD(S)) over all 65534 subsets | **+0.077** (median +0.099, positive on 58%) |
| corr(shipped-risk F(S), true RMSD(S)) | +0.069 |
| RMSD percentile of the EXACT optimum among all subsets | **57.8th** (worse than a coin flip) |
| mean k: exact-H optimum / exact-F optimum / ORACLE | 4.7 / 5.1 / 2.7 |

**This is the honest brief for the CVaR-VQE arm.** The instance family is real, non-trivial
and genuinely combinatorial (0.79 A of oracle headroom inside a 16-bit search space), but
(i) the classical optimiser already solves it essentially exactly and cheaply, and (ii) the
objective is anti-aligned with what we want: its exact minimiser sits at the 58th percentile
of achievable RMSD, so **any optimiser that finds a better H finds a worse structure**. A
quantum arm can be scored honestly against `exhaustive_surrogate.H` (the certified optimum
is in the dump for every instance) -- but a win on H is not a win on accuracy. The blocker
is the objective, not the optimiser.

---

## 4. THE LEARNED SET DECODER — the central experiment (`s12/agg_decoder.py`)

**Hypothesis (C2's untested class).** Previous in-band rankers saw only scalar poolings of
the candidate-vs-objective comparison and reached rank-correlation ~+0.20 for ~0 A. A
permutation-invariant model that reads the FULL per-pair deviation map might do better.

**Inputs (`s12/agg_features.py`).** Per target: a tensor `X (m=75, npairs, F=16)` — for
candidate c and pair (i,j): the candidate's distance; the distogram's expected distance and
sd; the raw and sd-normalised deviation; the per-pair Bayes risk of the candidate's
distance; the set mean / sd / median at that pair; the candidate's deviation from the set
mean, raw and normalised; expected-minus-set-mean; sequence separation; the distogram's
per-pair entropy; the risk of the set mean; and the candidate's rank among the set at that
pair. Plus `G (75, 12)` candidate scalars (z-scored shipped score, consensus, rank, org
flag, BLOSUM sim, RMSD to medoid and to the set average, rg, distogram MAE and correlation,
mean risk, fraction of pairs where it is the best candidate) and the full `P (75,75)`
candidate-vs-candidate CA-RMSD matrix.

**Architecture.** Per-(candidate,pair) MLP -> mean/max/std pooling over pairs -> concat the
candidate scalars -> 2 self-attention blocks over the 75 candidates with a learned
additive attention bias `-softplus(alpha) * P` (this is how the candidate-vs-candidate
geometry enters) -> per-candidate logit -> softmax weights. 51k parameters, CPU, 1 thread.

**Objective.** `L = CA-RMSD( sum_c w_c A_c , native )`, differentiable through the
coordinate average (the candidates are pre-superposed on the set medoid, so the average is
linear in w) with the optimal rotation detached — by the envelope theorem this is the exact
gradient of the aligned RMSD. NOT differentiable through the projection, as specified.
Comparison heads: `listwise` (CE against softmax(-rr/0.3)) and `regress` (MSE on z(rr) —
the objective class every previous scalar ranker used).

**Discipline.** 5 pinned folds; for test fold f the model trains on 3 folds and early-stops
on a 4th (never on f); normalisation statistics from the training targets only.

### 4a. Results (RAW, paired against avg75 = 3.0483)

| arm | mean | FAIL18 | other108 | d | CI95 | W/L | drop10 | rank skill* |
|---|---|---|---|---|---|---|---|---|
| **ORACLE positive control** (true RMSD leaked as 1 input) | **2.1466** | 4.620 | 1.734 | -0.9018 | [-1.016,-0.794] | **126/0** | -0.782 | **+0.755** |
| ORACLE, trained on 32 targets | 2.2087 | 4.726 | 1.789 | -0.8396 | [-0.949,-0.738] | 126/0 | -0.732 | +0.766 |
| ORACLE, trained on 8 targets | 2.5335 | 5.079 | 2.109 | -0.5149 | [-0.603,-0.431] | 125/1 | -0.416 | +0.583 |
| weight_noattn (no set attention) | 3.0324 | 5.841 | 2.564 | -0.0160 | [-0.037,+0.004] | 69/57 | +0.011 | +0.038 |
| **null (features permuted across candidates)** | **3.0335** | 5.881 | 2.559 | -0.0149 | [-0.036,+0.006] | 63/63 | +0.008 | -0.008 |
| weight, seed 1 | 3.0363 | 5.816 | 2.573 | -0.0121 | [-0.030,+0.006] | 68/58 | +0.009 | -0.027 |
| weight_hi (H=64, 3 layers, 150 ep) | 3.0461 | 5.808 | 2.586 | -0.0022 | [-0.034,+0.027] | 64/62 | +0.031 | +0.076 |
| **avg75 (production operator)** | **3.0483** | 5.832 | 2.584 | 0 | — | — | — | — |
| listwise head | 3.0615 | 5.883 | 2.591 | +0.0132 | [-0.008,+0.037] | 60/66 | +0.033 | +0.015 |
| **weight (primary arm)** | **3.0640** | 5.973 | 2.579 | +0.0156 | [-0.018,+0.055] | 61/65 | +0.043 | -0.035 |
| regress head (the old ranker objective) | 3.0839 | 5.976 | 2.602 | +0.0356 | [-0.026,+0.097] | 51/75 | +0.104 | -0.053 |
| — reference: the SHIPPED distogram score | | | | | | | | +0.029 |

*rank skill = Spearman(-model weight, true CA-RMSD) restricted to the near-native band
(pool_best + 1.5 A) **within the top-75** — a harder set than the K=500 in-band problem the
+0.20 record refers to, so the shipped score's own +0.029 here is the right yardstick.

**EMITTED** (through the real L-BFGS projection, `s12/results/agg_decproj_weight.json`):
decoder 3.2282 vs the 3.2028 incumbent, **d = +0.0253 [-0.010,+0.065], 61W/65L** — no gain.

### 4b. The learning curve — the decisive control

| training targets | learned decoder (real features) | ORACLE positive control |
|---|---|---|
| 8 | 3.0561 | **2.5335** |
| 16 | 3.0482 | — |
| 32 | 3.0824 | **2.2087** |
| 64 | 3.0554 | — |
| ~75 (3 folds, early stopped) | 3.0640 | **2.1466** |
| ~101 (all 4 folds, fixed 40 epochs, no early stop) | 3.0793 | — |
| — avg75 reference | 3.0483 | 3.0483 |

**The learning curve is FLAT — not rising.** Range 3.048-3.082 with no monotone trend, all
within the seed-to-seed spread (seed 0 = 3.0640, seed 1 = 3.0363). Meanwhile the ORACLE
control on the same harness, same folds, same optimiser is **already at half its final gain
with 8 training targets** and saturates by 32. So the harness converts a per-candidate
nativeness signal into accuracy at very small n; it simply is not finding one.

**Interpretation. C2's remaining hypothesis is now tested and it is a NEGATIVE.** The
untested class — learned aggregation over the full n x n candidate-vs-objective deviation
map — was given exactly the information it was said to be missing (per-pair deviations, set
statistics per pair, per-pair Bayes risk, the full candidate-vs-candidate RMSD matrix via
attention), a set-level objective that directly minimises the emitted RMSD proxy, and a
positive control proving the pipeline works. It reaches rank skill ~0.00 (the shipped
score's own in-band skill is +0.029) and 0 A. This is not sample-limitation on this axis:
the sample-size curve of a signal that DOES exist is steep and saturates at n=32.

**Nulls.** Permuting the per-candidate features (breaking the feature-to-candidate
correspondence while leaving the geometry channel intact) gives 3.0335 — indistinguishable
from, and nominally *better* than, the real-feature model. Deleting the attention blocks
(`weight_noattn`) gives 3.0324, also indistinguishable. Every "real" arm lies inside the
null band. The only arm outside it is the oracle control.

### 4c. Leakage audit

- `rr` and `nat_ca` enter only (i) the training loss, (ii) evaluation, (iii) arms explicitly
  labelled ORACLE. Every deployable input is a function of pool coordinates, the shipped
  leave-fold-out distogram, BLOSUM sim and the peptide/fragment origin flag.
- Normalisation statistics are computed on training targets only, inside the fold loop.
- The candidate superposition reference is the set medoid (native-free). The native is used
  in the loss only to align the *emitted* structure, never to align the candidates.
- Early stopping uses an inner fold disjoint from the test fold.
- The tie-break trap: every "pick the argmin of a native-free signal" measurement in
  sections 2 and 3 averages the outcome over the whole tied argmin set.

---

## 5. THE OPERATOR x OBJECTIVE-QUALITY SURFACE (coordinator request)
`s12/agg_surface.py`, `s12/agg_surface_proj.py`, `s12/agg_router.py`

**Question.** Is the optimal terminal operator a FUNCTION of objective quality — so that
"the production average is at the native-free optimum" holds only CONDITIONAL on the
current objective?

**Design.** Objective axis: the objective agent's `interp` corruption family, i.e. the REAL
distogram error vector scaled in amplitude (not i.i.d. noise), so the error STRUCTURE is the
predictor's own at every quality level:

    dtarget(a) = max( dtrue + a*(exp - dtrue), 1.5 )   a=0 native, a=1 shipped expectations
    score(c)   = mean_p | D[c,p] - dtarget(a)_p |      (the objective agent's L1 scorer)

plus the actual SHIPPED Bayes-risk objective as a labelled reference. Operator axis:
top-m coordinate average, m in {1,2,3,5,8,12,20,35,50,75,110,150,220,300,500}, each
superposing the SELECTED SUBSET on its OWN medoid (the production operator; m=1 is argmin).
All 126 targets, K=500 pool.

**Proxy validation (mine, independent).** Across 20 projected operator arms x 126 targets,
RAW vs EMITTED CA-RMSD: **r = 0.9956** per structure, mean offset **+0.132 A** (sd 0.170),
arm-ordering Spearman **0.949**. The proxy is sound; headline cells were re-run through the
real L-BFGS projection anyway.

### 5a. The surface (RAW mean CA-RMSD over 126)

| objective | MAE | m=1 | m=3 | m=5 | m=12 | m=20 | m=50 | m=75 | m=150 | m=500 | **best m** | best | cost of holding m=75 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| native | 0.00 | 1.994 | 1.837 | **1.832** | 1.891 | 1.930 | 2.099 | 2.218 | 2.495 | 3.396 | **5** | 1.832 | **+0.386** |
| a=0.1 | 0.23 | 2.014 | **1.862** | 1.864 | 1.912 | 1.951 | 2.124 | 2.234 | 2.506 | 3.396 | 3 | 1.862 | +0.372 |
| a=0.2 | 0.47 | 2.082 | 1.916 | **1.898** | 1.953 | 2.002 | 2.149 | 2.256 | 2.524 | 3.396 | 5 | 1.898 | +0.358 |
| a=0.3 | 0.70 | 2.203 | 2.026 | **1.978** | 2.002 | 2.065 | 2.195 | 2.286 | 2.551 | 3.396 | 5 | 1.978 | +0.309 |
| a=0.4 | 0.94 | 2.378 | 2.157 | 2.101 | 2.100 | 2.138 | 2.256 | 2.343 | 2.594 | 3.396 | 8 | 2.068 | +0.275 |
| a=0.5 | 1.17 | 2.531 | 2.297 | 2.276 | 2.272 | **2.261** | 2.351 | 2.432 | 2.641 | 3.396 | 20 | 2.261 | +0.172 |
| a=0.6 | 1.40 | 2.714 | 2.522 | 2.440 | 2.419 | **2.413** | 2.473 | 2.531 | 2.717 | 3.396 | 20 | 2.413 | +0.119 |
| a=0.7 | 1.64 | 2.954 | 2.724 | 2.654 | 2.604 | 2.593 | 2.612 | 2.666 | 2.798 | 3.396 | 35 | 2.587 | +0.079 |
| a=0.8 | 1.87 | 3.168 | 2.976 | 2.892 | 2.787 | 2.777 | **2.760** | 2.799 | 2.884 | 3.396 | 50 | 2.760 | +0.039 |
| a=0.9 | 2.10 | 3.331 | 3.140 | 3.081 | 2.999 | 2.972 | **2.928** | 2.931 | 2.997 | 3.396 | 50 | 2.928 | +0.004 |
| a=1.0 | 2.34 | 3.509 | 3.329 | 3.256 | 3.197 | 3.136 | 3.086 | 3.078 | 3.093 | 3.396 | 110 | 3.075 | +0.003 |
| a=1.2 | 2.81 | 3.911 | 3.657 | 3.574 | 3.481 | 3.407 | 3.346 | 3.308 | **3.252** | 3.396 | 150 | 3.252 | +0.055 |
| a=1.5 | 3.51 | 4.248 | 4.049 | 3.932 | 3.838 | 3.793 | 3.643 | 3.574 | 3.451 | **3.396** | 500 | 3.396 | +0.178 |
| **SHIPPED (Bayes risk)** | 2.34 | 3.454 | 3.274 | 3.231 | 3.132 | 3.091 | 3.068 | **3.048** | 3.072 | 3.396 | **75** | 3.048 | **0.000** |

**The coordinator's interaction is REAL and now quantified — but the direction was stated
slightly wrong.** Three findings:

1. **The optimal cardinality shrinks monotonically as the objective improves**, from
   m=500 (worse-than-shipped objective: do not select at all) through m=75-110 at shipped
   quality to **m=3-5 at perfect quality**. That curve is the deliverable and it is smooth
   and monotone in MAE with no plateau.
2. **Argmin never becomes optimal.** Even with the exact native distance matrix, m=1 gives
   1.994 while m=5 gives **1.832**. The apparent "argmin beats averaging at perfect
   objective" is a comparison of m=1 against the *fixed* m=75 (1.994 vs 2.218), and both
   lose to a small average. The reconciliation the coordinator asked about:
   S10-5's "perfect distance ranker over K=500 = 1.994" is exactly my **(a=0, m=1) = 1.994**;
   the objective agent's "native matrix, top-75, average = 2.395 emitted" is my
   **(a=0, m=75) = 2.218 raw + 0.13-0.16 projection = ~2.37**. Both numbers are right; they
   are two cells of one surface, and **neither is the optimum of that surface**.
3. **Where the shipped objective sits.** The m=1/m=75 crossover is at **MAE ~= 0.85**
   (a=0.3: 2.203 vs 2.286, argmin wins; a=0.4: 2.378 vs 2.343, average wins). The shipped
   objective is at **MAE 2.34, i.e. 2.7x worse than the crossover** — far on the averaging
   side. And within the shipped row the surface is FLAT over m in [35,110] (3.069, 3.068,
   3.048, 3.050): **m=75 is not just safe, it is the argmin of the shipped row.** So my
   section-1 result stands unchanged at the current operating point.

**But the interaction bites the moment the objective improves.** Holding m=75 costs
+0.004 A at MAE 2.10, +0.079 at 1.64, +0.172 at 1.17 and +0.386 at MAE 0. Since the record
says the projection path needs MAE ~0.87 to reach 2 A, **any objective improvement that
matters lands squarely in the region where m must move** — at MAE 0.87 roughly 0.28 A of
the win would be thrown away by keeping m=75. The operator is not a live tuning axis today
and *is* a live axis for any future objective work. Operator and objective cannot be
optimised independently; the correct statement is "m* = m*(objective quality)".

### 5b. FAIL18 vs other-108 (coordinator item 4)

| objective | subgroup | m=1 | m=5 | m=75 | best m | best |
|---|---|---|---|---|---|---|
| native | FAIL18 | 2.624 | **2.485** | 3.381 | 5 | 2.485 |
| native | other108 | 1.889 | 1.723 | 2.024 | 3 | **1.720** |
| a=0.5 (MAE 1.17) | FAIL18 | 4.182 | 3.706 | 3.922 | 35 | 3.657 |
| a=0.5 | other108 | 2.256 | 2.038 | 2.184 | 20 | 2.026 |
| SHIPPED | FAIL18 | 6.008 | 5.860 | 5.832 | 500 | 5.220 |
| SHIPPED | other108 | 3.028 | 2.793 | 2.584 | 75 | 2.584 |

**Yes — the FAIL18 "operator gap" is mostly an operator gap, and a different operator closes
most of it once the objective is good.** The forensics agent measured a perfect distogram
still emitting 3.640 A on FAIL18 against a 2.285 A single-best floor, a 1.36 A gap
attributed to the operator. My surface says: with a perfect objective, **m=75 gives 3.381
raw (~3.51 emitted, matching their 3.640) but m=5 gives 2.485 raw (~2.62 emitted)**. So
**0.90 A of that 1.36 A gap is recovered simply by moving m from 75 to 5**, leaving ~0.20 A
against the pool-best floor. The gap is not "irreducible operator loss"; it is the cost of a
cardinality tuned for a bad objective being applied to a good one. Note also that FAIL18
prefers a *larger* m than other-108 at every intermediate quality (m=35 vs m=20 at MAE 1.17,
m=500 vs m=75 at shipped) — the modes are further apart there so the average has to be
broader to stay safe.

### 5c. Is there a DEPLOYABLE per-target router for m? (coordinator item 3) — **NO**

Labels: the per-target best m on the shipped row (training labels only, leave-fold-out).
Features: 19 deployable per-target statistics (length; pool score mean/sd/skew and the
rank-1-to-75 and 75-to-500 score gaps; distogram sd mean/max and expected mean/sd; top-75
consensus mean/sd/median; BLOSUM sim mean/sd over the top-75 and over the pool; org
fraction; rg spread). Multi-output ridge on the per-cardinality RMSD, argmin of the
prediction, 5 pinned folds, label-permutation null.

| grid {20,75,300} | mean | d vs const m=75 | CI95 | W/L | drop10 |
|---|---|---|---|---|---|
| constant m=75 | 3.0483 | 0 | — | — | — |
| **LFO ridge router** | **3.0462** | **-0.0021** | **[-0.062,+0.055]** | **39/47** | +0.066 |
| label-permutation null | 3.1437 | +0.0954 | [+0.038,+0.159] | 23/41 | +0.138 |
| ORACLE router (ceiling) | 2.8037 | -0.2446 | [-0.314,-0.183] | 83/0 | -0.158 |

| grid {1,75} | mean | d | CI95 | W/L |
|---|---|---|---|---|
| constant m=75 | 3.0483 | 0 | — | — |
| LFO ridge router | 3.0588 | +0.0104 | [+0.001,+0.028] | 1/6 |
| ORACLE router | 2.9527 | -0.0956 | [-0.154,-0.051] | 31/0 |

**Negative, as the record predicted.** The oracle router is worth 0.245 A over a 3-point
grid; the leave-fold-out router captures **0.9% of it** (-0.0021 A), its CI spans zero, it
loses more targets than it wins (39/47), and the effect reverses when the top 10 targets are
dropped. The label distribution is nearly uniform (35/43/48 across {20,75,300}), i.e. the
per-target preference is real but unpredictable from any deployable summary I could build.
This is the same null as the record's LFO routing AUC 0.558 vs permutation 0.511.

### 5d. Adjudication of the disagreement the coordinator flagged

**The coordinator's reading is correct, and I confirm it with one refinement.**

The objective agent's "0.68 A the terminal operator loses between the best pool member
(1.711) and the emitted average (~2.39 with a perfect objective)" and my "no native-free
operator beats the production average by more than 0.004 A" are both true and are not in
tension, because they are statements about different variables:

- **My statement is about operator FORM at fixed selection quality.** Given the members the
  current objective actually puts in the top-75, no trimming, weighting, geometric-median,
  clustering or subset-selection rule improves on the uniform mean. Section 3 proves the
  strong version: even the *exact* combinatorial optimum of a native-free selection
  objective, certified by 2^16 enumeration on all 126 targets, is 0.13 A WORSE than not
  selecting at all, and sits at the 58th percentile of achievable RMSD.
- **Their statement is about the gap to the best member**, which is a statement about
  RANKING. Closing it requires knowing *which* members are good. Sections 1c, 3 and 4 all
  say that information is not available: the optimal weights correlate +0.03 with the
  shipped score, the exact subset optimum is anti-aligned with RMSD, and a set decoder with
  the full deviation map reaches rank skill ~0.00 while the same harness reaches +0.755 when
  the label is leaked.

So: **the 0.68 A is real but it is not an operator budget — it is the ranking budget seen
through the operator.** The recommendation "spend effort on the terminal operator" is not
actionable as stated. **The refinement is section 5a**: there IS one operator degree of
freedom that is genuinely coupled to objective quality — the cardinality m — and it is worth
0 A today and up to 0.386 A if the objective ever improves. The actionable form of their
recommendation is therefore: *do not spend effort on the operator now; spend it on the
objective, and re-tune m the moment the objective moves.*

### 5e. EMITTED confirmation of the headline cells (real L-BFGS projection, all 126)
`s12/results/agg_surface_proj.json`

| cell | RAW | **EMITTED** | projection cost | FAIL18 (emitted) | other108 (emitted) | d vs 3.2052 | CI95 | W/L |
|---|---|---|---|---|---|---|---|---|
| perfect objective, m=1 (argmin) | 1.9938 | 1.9918 | **-0.002** | 2.623 | 1.887 | -1.2134 | [-1.448,-0.995] | 116/10 |
| **perfect objective, m=5** | 1.8319 | **1.9483** | +0.116 | 2.655 | 1.830 | **-1.2569** | [-1.481,-1.049] | **120/6** |
| perfect objective, m=25 | 1.9679 | 2.1231 | +0.155 | 3.039 | 1.970 | -1.0821 | [-1.291,-0.888] | 119/7 |
| perfect objective, m=75 | 2.2183 | 2.3946 | +0.176 | 3.652 | 2.185 | -0.8106 | [-0.994,-0.640] | 116/10 |
| MAE 1.17 objective, m=20 | 2.2606 | 2.4202 | +0.160 | 3.872 | 2.178 | -0.7850 | [-0.945,-0.640] | 114/12 |
| MAE 1.17 objective, m=75 | 2.4323 | 2.6002 | +0.168 | 4.144 | 2.343 | -0.6050 | [-0.758,-0.466] | 109/17 |
| shipped, m=5 | 3.2312 | 3.3540 | +0.123 | 5.981 | 2.916 | +0.1488 | [+0.048,+0.251] | 44/82 |
| shipped, m=20 | 3.0907 | 3.2507 | +0.160 | 5.975 | 2.797 | +0.0455 | [-0.022,+0.111] | 53/73 |
| **shipped, m=75 (production)** | 3.0483 | **3.2052** | +0.157 | 6.034 | 2.734 | 0 | — | — |

Three things the projection adds to the raw surface:

1. **The projection penalty is monotone in m** (0.00 at m=1, 0.116 at m=5, 0.155 at m=25,
   0.176 at m=75). Averaging more members moves the cloud further off the ideal-geometry
   manifold. This shifts the optimal cardinality DOWN relative to the raw surface, and it
   compresses the m=5-over-m=1 margin from 0.162 A raw to 0.043 A emitted. Argmin still
   does not win, but emitted it is much closer than raw suggests.
2. **A perfect distance objective with the right cardinality emits 1.948 A** — below the
   sprint's 2.0 A target, on the pool the pipeline already retrieves, using nothing but a
   different m. With the production m=75 the same perfect objective emits 2.395.
3. **FAIL18, emitted:** perfect objective at m=75 gives **3.652** (the forensics agent's
   3.640, reproduced), and at m=5 gives **2.655**. So **0.997 A of the 1.36 A FAIL18
   "operator gap" is closed by cardinality alone**, confirmed through the real projection.

### 5f. A correction the coordinator should carry: "perfect filter" is ambiguous

The forensics number "perfect filtering at m=25 emits 1.644" is a **perfect-RMSD filter**
(rank the pool by TRUE CA-RMSD). That is my `ORC_top25avg_pool`: raw 1.608, emitted 1.720.
It is **not** what a perfect *distance objective* achieves: scoring the same pool with the
exact native distance matrix and taking the top 25 emits **2.123** (raw 1.968). The 0.40 A
between them is the part of the filter budget that no distance-space objective can reach,
because L1-in-distance-space rank is not CA-RMSD rank.

| "perfect filter", m=25 | raw | emitted |
|---|---|---|
| rank by TRUE CA-RMSD (the 1.644 number) | 1.608 | 1.720 |
| rank by the EXACT NATIVE DISTANCE MATRIX (the best any distogram can be) | 1.968 | 2.123 |
| best cell of the whole perfect-objective row (m=5) | 1.832 | **1.948** |

So the filter budget quoted as 1.40 A is really **~1.26 A reachable by a perfect distogram
(3.205 -> 1.948, at the right m)** plus ~0.17 A that is only reachable with RMSD knowledge.
That does not change the conclusion — the filter is still where the accuracy is, and a
perfect distogram still crosses 2.0 A — but the target for distogram work is 1.95, not 1.64.

---

## 6. THE DECODER RE-RUN TO THE COORDINATOR'S SPECIFICATION (v2)
`s12/agg_features.py` (v2 tensor), `s12/agg_decoder.py`, `s12/agg_dec_run2.py`,
`s12/agg_dec_score.py`, `s12/results/agg_dec_v2*.json`, `s12/results/agg_dec_terminals.json`

Four changes were requested and all four are in place:

| requested | status |
|---|---|
| optimise the RMSD of the AGGREGATE, not a per-candidate rank/regression | **already the primary head** — `L = CA-RMSD(sum_c w_c A_c, native)`, exact gradient through the coordinate average (rotation detached, envelope theorem). `regress` and `listwise` are the *comparison* arms and both do worse. |
| SIGNED per-pair deviation map, no magnitude-only pooling | X keeps 8 signed columns; v2 adds two more that express COHERENT DIRECTION explicitly: the candidate's mean signed deviation **within each \|i-j\| shell**, and its **global mean signed deviation**. G adds absolute rg, the **distogram-implied rg**, **rg - rg_implied** (signed), n, and mean(d)-mean(e). F: 16 -> 18, GF: 12 -> 18. |
| FAIL18 held OUT OF TRAINING entirely | `no_fail18_train`: for test fold f the model trains on `{fold != f} \ FAIL18`. FAIL18 is a pure held-out group everywhere below. |
| score every arm through argmin, m=25 and m=75; report retained-set MEAN | `s12/agg_dec_score.py`, tables 6b/6c. |

### 6a. v2 results (RAW weighted average, paired vs avg75 = 3.0483)

| arm | mean | FAIL18 (pure held-out) | other108 | d | CI95 | W/L | drop10 | drop20 |
|---|---|---|---|---|---|---|---|---|
| **v2_oracle** (label leaked) | **2.1460** | 4.629 | 1.732 | -0.9024 | [-1.017,-0.798] | **125/1** | -0.786 | -0.695 |
| v2_oracle, 32 train targets | 2.2004 | — | — | -0.848 | [-0.959,-0.744] | 126/0 | — | — |
| v2_oracle, 8 train targets | 2.5342 | — | — | -0.514 | [-0.595,-0.437] | 125/1 | — | — |
| **v2 (signed, seed 0)** | 3.0258 | 5.869 | 2.552 | -0.0225 | [-0.059,+0.014] | 69/57 | **+0.025** | +0.042 |
| v2, seed 1 | 3.0337 | 5.889 | 2.558 | -0.0147 | [-0.051,+0.023] | 61/65 | +0.027 | +0.046 |
| v2, 64 train targets | 3.0358 | 5.915 | 2.556 | -0.0126 | [-0.051,+0.028] | 64/62 | +0.034 | +0.052 |
| v2, 8 train targets | 3.0433 | 5.801 | 2.584 | -0.0050 | [-0.014,+0.002] | 63/63 | +0.004 | +0.005 |
| v2, 32 train targets | 3.0456 | 5.825 | 2.583 | -0.0027 | [-0.030,+0.024] | 62/64 | +0.030 | +0.041 |
| **avg75 (production)** | **3.0483** | 5.832 | 2.584 | 0 | — | — | — | — |
| v2, 16 train targets | 3.0491 | 5.819 | 2.588 | +0.0008 | [-0.012,+0.013] | 55/71 | +0.015 | +0.020 |
| **v2_null (labels/features permuted)** | 3.0499 | 5.875 | 2.579 | +0.0016 | [-0.027,+0.031] | 65/61 | +0.029 | +0.044 |
| **v2_abs (MAGNITUDE-ONLY ablation)** | 3.0624 | 5.974 | 2.577 | +0.0141 | [-0.030,+0.060] | 58/68 | +0.061 | +0.078 |
| v2_hi (H=64, 3 layers, 150 ep) | 3.0957 | — | — | +0.047 | — | — | — | — |
| v2, all 101 train targets (fixed 40 ep) | 3.1471 | 6.084 | 2.658 | +0.0988 | [+0.008,+0.191] | 53/73 | +0.184 | +0.236 |

**EMITTED** (real L-BFGS projection, all 126, `s12/results/agg_decproj_v2.json`):

| | emitted | d vs 3.2028 | CI95 | W/L | drop10 | drop20 | per-fold |
|---|---|---|---|---|---|---|---|
| v2 | **3.1840** | **-0.0189** | [-0.058,+0.020] | **63/63** | +0.031 | +0.052 | -.047/+.027/+.026/-.093/-.012 |
| — FAIL18 (pure held-out, n=18) | 6.0711 | **+0.0379** | [-0.108,+0.180] | 7/11 | | | |
| — other-108 | 2.7028 | -0.0283 | [-0.069,+0.009] | 56/52 | | | |

**Verdict: negative on every criterion in the brief.** |d| = 0.019 < the 0.03 threshold; the
CI spans zero; W/L is exactly 63/63; the sign REVERSES when the top 10 targets are dropped;
two of five folds are positive; and it is **worse on the FAIL18**, the group that carries
the effect. The apparent -0.023 raw gain is inside the seed-to-seed spread (seed 0 -0.023,
seed 1 -0.015) and inside the null band (v2_null +0.002).

### 6b. The learning curve — flat, against a steep oracle curve

| training targets | v2 (signed features) | v2_oracle (label leaked) |
|---|---|---|
| 8 | 3.0433 | **2.5342** |
| 16 | 3.0491 | — |
| 32 | 3.0456 | **2.2004** |
| 64 | 3.0358 | — |
| ~75 (3 folds, early stopped) | 3.0258 | **2.1460** |
| ~101 (4 folds, fixed 40 ep, no early stop) | 3.1471 | — |
| avg75 reference | 3.0483 | 3.0483 |

**The v2 curve is FLAT** — total range 3.026-3.049 over an 8x change in training data, with
no monotone trend, and the two extreme points (8 and 75) differ by 0.018 A, less than the
seed spread. The n=101 arm is worse, but it is confounded (no early stopping), so I do not
read a decline from it. **The same harness on the same folds with a signal that definitely
exists is already at 57% of its final gain at n=8 and at 92% by n=32.** So the answer to
the question C2 left open is: **flat, not still rising** — this hypothesis class is not
sample-limited, it is signal-limited.

### 6c. DOUBLE-SCORING through all three terminals (the sprint-wide rule)

Every arm read off the SAME weight vector: `argmin` = highest-weight candidate; `avgM` =
uniform coordinate average of the M highest-weight candidates; `wavg75` = the trained
weighted average. `setmeanM` = mean TRUE CA-RMSD of the retained set (the adversarial
agent's high-power surrogate). Baseline = the shipped distogram score used the same way.

| arm | argmin | avg10 | avg25 | wavg75 | setmean10 | setmean25 | setbest25 |
|---|---|---|---|---|---|---|---|
| SHIPPED score (baseline) | 3.4540 | 3.1372 | 3.0793 | 3.0483 | 3.4779 | 3.5016 | 2.6087 |
| **v2_oracle** (label leaked) | **2.3634** | **2.2820** | **2.5214** | **2.1460** | **2.7386** | **3.0276** | **2.3062** |
| v2_oracle, 8 train targets | 2.5710 | 2.3860 | 2.5698 | 2.5342 | 2.8409 | 3.0820 | 2.3096 |
| v2 (signed) | 3.6581 | 3.1089 | 3.0355 | 3.0258 | 3.5271 | 3.5125 | 2.5379 |
| v2_null | 3.4493 | 3.0629 | 3.0419 | 3.0499 | 3.5264 | 3.5363 | 2.4757 |
| v2_abs (magnitude only) | 3.6098 | 3.1483 | 3.0872 | 3.0624 | 3.5599 | 3.5567 | 2.5911 |
| weight (v1) | 3.6030 | 3.1855 | 3.1240 | 3.0640 | 3.5926 | 3.5788 | 2.6055 |
| regress head | 3.4043 | 3.1439 | 3.0998 | 3.0839 | 3.4516 | 3.4587 | 2.6540 |

Paired vs the shipped score, per terminal:

| arm | argmin | avg25 | wavg75 |
|---|---|---|---|
| **v2_oracle** | **-1.091 [-1.250,-0.945] 124/1** | **-0.558 [-0.644,-0.478] 123/3** | **-0.902 [-1.017,-0.798] 125/1** |
| v2_oracle, n=8 | -0.883 [-1.031,-0.740] 113/8 | -0.509 [-0.598,-0.428] 122/4 | -0.514 [-0.595,-0.437] 125/1 |
| **v2 (signed)** | **+0.204 [+0.040,+0.359] 51/72** | -0.044 [-0.122,+0.034] 64/62 | -0.022 [-0.059,+0.014] 69/57 |
| v2_null | -0.005 [-0.169,+0.160] 56/67 | -0.037 [-0.087,+0.013] 70/56 | +0.002 [-0.027,+0.031] 65/61 |
| v2_abs | +0.156 [-0.006,+0.312] 56/65 | +0.008 [-0.077,+0.095] 58/68 | +0.014 [-0.030,+0.059] 58/68 |
| weight (v1) | +0.149 [-0.013,+0.316] 57/66 | +0.045 [-0.027,+0.129] 62/64 | +0.016 [-0.018,+0.055] 61/65 |
| regress | -0.050 [-0.198,+0.096] 60/62 | +0.021 [-0.062,+0.103] 61/65 | +0.036 [-0.026,+0.097] 51/75 |

**This closes the "a real gain would be invisible through the average terminal" escape
route.** The oracle control is loudly visible through ALL THREE terminals — **-1.09 through
argmin, -0.56 through avg25, -0.90 through the weighted average**, 123-125 wins out of 126
in each. The argmin terminal, which the adversarial law says is the *most* sensitive to
ranking quality, is exactly where the real-feature decoder is at its WORST: v2 is
**+0.204 A worse than the shipped score as a ranker** (51W/72L, CI excludes zero). So the
decoder is not a good ranker whose gain is being absorbed by the operator; it is a worse
ranker than the objective it was given, and the operator is hiding that fact, not a gain.

**On the retained-set MEAN — the coordinator's primary surrogate — the answer is the same
and it is higher-powered.** The shipped score already delivers setmean25 = 3.5016; v2
delivers **3.5125 (worse)**; the magnitude-only ablation 3.5567; the null 3.5363. The
oracle delivers **3.0276** at m=25 and 2.7386 at m=10. **No learned arm moves the set mean
at all.** I independently reproduce the adversarial agent's transfer law from my own arms:
across (shipped m=10 -> oracle m=10) the implied slope is
(3.1372-2.2820)/(3.4779-2.7386) = **1.157**, against their 1.162.

### 6d. Does the SIGN carry the information? (the sharp null)

| comparison | d | CI95 | W/L | drop10 |
|---|---|---|---|---|
| v2 (signed) vs v2_abs (magnitude only) | **-0.0366** | **[-0.065,-0.011]** | **82/44** | -0.003 |
| v2 (signed) vs v2_null (features permuted) | -0.0240 | [-0.074,+0.023] | 61/65 | +0.035 |
| v2 (signed) vs avg75 (no model) | -0.0225 | [-0.059,+0.014] | 69/57 | +0.025 |

Read honestly: the only comparison whose CI excludes zero is **signed vs magnitude-only**,
and it excludes zero because **the magnitude-only model is worse than doing nothing**
(3.0624 vs avg75 3.0483), not because the signed model is better than doing nothing
(3.0258, CI spans zero, 69/57). Destroying the sign makes the model actively harmful;
keeping it makes the model harmless. **That is not evidence that the sign carries an
exploitable signal** — it is evidence that a magnitude-only view of the deviation map is a
misleading input. The coherent-direction hypothesis is therefore *consistent* with what I
measure but is not *supported* by it: giving the model the shell-wise and global signed
deviations and the rg-vs-implied-rg gap (exactly the coherent-direction coordinates) does
not move emitted RMSD, the retained-set mean, or rank skill.

---

## 7. NEGATIVES, LEAKAGE AUDIT, AND WHAT THE COORDINATOR SHOULD DO NEXT

### 7a. Complete list of negatives (all reported at full strength)

| # | thing tried | arms | result |
|---|---|---|---|
| N1 | operator FORM on the shipped top-75 (trimmed x4, geometric median, iterated averaging, score/rank/consensus/joint weightings x10, cluster-then-average) | 30 | none beats uniform averaging by > 0.004 A; no CI excludes zero on the good side |
| N2 | cardinality ladder over the K=500 pool at shipped objective quality (m = 10..500) | 6 | m=75 is the argmin; the row is flat over m in [35,110] |
| N3 | multi-hypothesis modes: 3 clusterings x k=2..5 x 13 native-free picking rules + mode-aware averaging + soft mixtures | 192 | ALL worse than the single global average; best +0.013 A (and that is the max over 192) |
| N4 | combinatorial subset selection on the exact shipped objective (greedy, greedy+1-opt, SA x2 seeds, fixed-k greedy) | 24 | all >= avg75; harder selection is monotonically worse |
| N5 | EXACT 2^16 optimum of the quadratic subset objective, all 126 targets | 2 | +0.127 A worse than not selecting; optimum sits at the 58th RMSD percentile |
| N6 | learned set decoder over the full deviation map: weight / listwise / regress heads, +/- attention, H=32/64, 2/3 layers, 2 seeds, 6 training-set sizes, v1 and v2 feature sets | 28 | emitted 3.184 vs 3.203, d -0.019 [-0.058,+0.020], 63W/63L, reverses on drop-top-10, worse on FAIL18 |
| N7 | the magnitude-only ablation | 1 | worse than doing nothing (3.062) — sign matters only in that removing it hurts |
| N8 | deployable per-target router for the averaging cardinality m (19 features, ridge, LFO, 4 grids) | 8 | captures 0.9% of a 0.245 A oracle; CI spans zero; 39W/47L; reverses on drop-top-10 |
| N9 | learning curve of the decoder over 8/16/32/64/75/101 training targets | 6 | FLAT; the same harness with a leaked label is at 57% of its gain by n=8 |

### 7b. Positives

| # | finding | value |
|---|---|---|
| P1 | the optimal averaging cardinality m* is a monotone function of objective quality | m* = 500 (MAE 3.5) -> 75-110 (MAE 2.34, today) -> 20 (MAE 1.2) -> 3-5 (MAE 0) |
| P2 | a PERFECT distance objective at m=5 emits **1.948 A** — under the sprint target — on the pool the pipeline already retrieves | vs 2.395 at the production m=75 |
| P3 | the FAIL18 "operator gap" is mostly a CARDINALITY gap | perfect objective: m=75 emits 3.652, m=5 emits 2.655 (0.997 A of the 1.36 A recovered), confirmed through the real projection |
| P4 | "perfect filter" is ambiguous and the sprint should stop conflating the two readings | perfect-RMSD filter at m=25 emits 1.720; perfect-DISTANCE-objective filter at m=25 emits 2.123 |
| P5 | the oracle optimal weights are orthogonal to the shipped score (+0.03) and POSITIVELY correlated with outlyingness (+0.23) | the consensus/medoid heuristic has the wrong sign |
| P6 | the subset problem splits ~50/50 into selection and weighting | uniform average over the oracle support 2.441; optimally weighted 1.936; from 3.048 |
| P7 | independent reproduction of the adversarial operator transfer law from my own arms | implied slope 1.157 vs their 1.162 |

### 7c. Leakage audit

- `rr` / `nat_ca` appear only in (i) evaluation, (ii) the training loss of the decoder under
  leave-fold-out discipline, (iii) arms explicitly prefixed `ORC_`/`ORACLE`. No inference-time
  decision anywhere reads them.
- Every deployable feature is a function of pool CA coordinates, the shipped leave-fold-out
  distogram (expected / sd / prob / risk), BLOSUM similarity, the peptide-vs-fragment origin
  flag and sequence length.
- Fold discipline: 5 pinned folds from `instrument.targets()`; the decoder trains on 3 folds
  and early-stops on a 4th, never the test fold; normalisation statistics come from training
  targets only, computed inside the fold loop. v2 additionally removes FAIL18 from training.
- Superposition references are always native-free (the set or subset medoid). The native is
  used only to align the FINAL emitted structure for scoring.
- Tie-breaking: every "argmin of a native-free signal" measurement averages the outcome over
  the whole tied argmin set (the pool-order leak documented in the memory record).
- benchmark60 untouched: nothing in `s12/agg_*` reads `results/benchmark_manifest.json`,
  `s9/final_cache/*` or `bench_results/` except `bench_results/cache/<PROD_KEY>/<pdb>.json`
  for the 126 tuning targets, via `instrument.shipped_record`. dev24 untouched.
- Multiple comparisons: sections 2 and 3 report the MAXIMUM over 192 and 24 arms
  respectively and are read as upper bounds, not as discoveries.

### 7d. Ranked recommendations

1. **Stop work on the terminal operator as a source of accuracy, and record the reason.**
   The operator FORM is exhausted (30 arms, 192 mode arms, 24 combinatorial arms, an exact
   2^16 optimum, and a set decoder). It is at its conditional optimum and the conditional is
   the objective.
2. **Adopt m* = m*(objective quality) as a standing rule.** m is the one operator degree of
   freedom that is genuinely coupled to the objective; it is worth 0 A today and up to
   0.39 A raw / 0.45 A emitted the moment the distogram improves. Any future objective
   result must be reported at its own best m, or it will understate itself. Concretely: at
   MAE 1.2 use m~20, at MAE <0.5 use m~5.
3. **Retarget the distogram at 1.95 A, not 1.64 A.** A perfect distance objective cannot
   reach the perfect-RMSD-filter number; 1.948 emitted (m=5) is the true ceiling of the
   distogram line, and it does clear 2.0 A.
4. **Do not give the CVaR-VQE arm this subset problem expecting accuracy.** Instances are
   dumped and certified (`agg_subset_instances.json`), so a quantum arm CAN be scored
   honestly on objective value — but the exact classical optimum is already known, greedy
   is optimal on 40% of instances, and finding a better optimum makes the structure worse.
   If the mandated component must run, run it as an OPTIMISER BENCHMARK against the
   certified optimum and say so; do not report it as an accuracy method.
5. **Where the remaining information must come from.** Three independent measurements in
   this report say the missing channel is not in the candidate set: the oracle weights are
   orthogonal to the objective (+0.03); the exact optimum of the native-free objective is at
   the 58th RMSD percentile; and a set model with the full signed deviation map, the shell-
   wise coherent-direction coordinates and the rg-vs-implied-rg gap reaches rank skill 0.00
   while the same harness reaches +0.755 on a leaked label. The information needed to rank
   *within* the retained set is not present in the deployable view. It has to be added, and
   the only lever this report can see is a better distance/shape prior — which is exactly
   where P1-P3 say the accuracy is.
6. **If one more decoder experiment is run**, run it on the K=500 pool rather than the
   top-75 (the top-75 is already the filtered set, so the decoder is being asked to re-rank
   a set the same objective already ranked). Feature tensors for `topk200` are built and
   cached (`s12/cache/agg_feat/*_topk200.npz`); this is the one variant I did not have time
   to train. I would still expect a negative, but it is the cheapest remaining shot and it
   is the only version of the experiment that could in principle rescue members the shipped
   filter discarded — which is where the FAIL18 loss actually is.

---

## 8. Deliverables index

Code (all under `s12/`, nothing outside `s12/agg_*` was modified):
`agg_common.py` cache+operators | `agg_meta.py` pool metadata | `agg_ladder.py` operator ladder
`agg_modes.py` multi-hypothesis | `agg_subset.py` combinatorial selection + QUBO instances
`agg_features.py` (m x npairs x F) deviation tensor | `agg_decoder.py` set decoder
`agg_dec_run.py` / `agg_dec_run2.py` drivers | `agg_dec_proj.py` decoder projection
`agg_dec_score.py` three-terminal double-scoring | `agg_surface.py` + `agg_surface_proj.py`
operator x objective surface | `agg_router.py` cardinality router | `agg_report.py` stats

Results (`s12/results/`):
`agg_ladder_raw.json` / `agg_ladder_proj.json` � 44 operator arms raw, 20 emitted
`agg_ladder_diag.json` � mode sizes, oracle simplex support, consensus
`agg_modes.json` � 192 native-free multi-hypothesis arms
`agg_subset_exact.json` � greedy / local-search / SA / fixed-k on the full 75
**`agg_subset_instances.json`** � 126 certified 16-bit QUBO instances for the quantum arm
`agg_dec_*.json` (28 arms) + `agg_decproj_weight.json` / `agg_decproj_v2.json` (emitted)
**`agg_dec_terminals.json`** � every arm through argmin / avg10 / avg25 / avg75 + set means
`agg_surface.json` (14 objective levels x 15 cardinalities x 126) / `agg_surface_proj.json`
`agg_router.json` � the cardinality router and its null

Note for the coordinator: the operator x objective-quality surface (section 5) was already
complete when it was reassigned to the adversarial agent. It is reported here rather than
duplicated; sections 5a-5f and `agg_surface*.json` are available for their use.
