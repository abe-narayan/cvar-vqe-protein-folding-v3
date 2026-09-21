# PREREG — S32 LANE P (POOL, FILTER, SELECTION)

**Committed before the first number exists.** Branch `s26`. Endpoint: mean **built-chain** Cα RMSD,
`tuning126`, n = 126, production **3.2105 Å**. Cloud (3.0483) and set mean (3.5507) are different
objects and are never differenced against the chain.

Charter sections served: §10 (decompose the pool), §27/§28 (candidate generation), §49 (FAIL18),
§56 (earliest irreversible loss). Contract rules that bind hardest here: **9** (price best-of-K),
**10** (a random control needs its own draw distribution), **12** (FAIL18 is circular), **16** (the
projection price is a property of the object), **17** (the pool is not the bottleneck until shown).

---

## 0. What this lane must not assume

The brief's ladder attributes `+0.4357 Å` to "the filter K=500 → 128". Two facts about the code
change what that sentence can mean, and both are established by reading `core/pipeline.py` before
any measurement:

1. **The 500 → 128 cut is the DISTOGRAM BAYES-RISK score, not BLOSUM.** `retrieve()` builds the
   K = 500 pool by BLOSUM62 sum over the whole leakage-safe window universe
   (`core/pipeline.py:698-700`); `filter_pool()` then sorts by `pool["sc"]`, the distogram risk
   (`core/pipeline.py:755-760`). The *BLOSUM* filter is **universe → 500**, a rung nobody has
   priced, and it is strictly earlier than every rung on the S29 ladder.
2. **The 128 window does not exist on the deployable path.** `want = max(cfg.m, 128 if cfg.quantum
   else 0)`; production is `quantum = False`, so the deployable filter is **500 → 75**. A ladder
   rung measured on the 128-window is a statement about the quantum register, not about production.

Therefore this lane measures the filter at **both** cut points and names which one each number is
about, and it opens the rung **above** K=500.

---

## H-P1 — The six-way pool decomposition, measured INDEPENDENTLY

**HYPOTHESIS.** The endpoint loss can be attributed to six stages that are separately measurable,
and the attribution is *not* the cumulative ladder: a cumulative ladder charges each rung with the
residual of every earlier one.

**MECHANISM.** Each stage is idealised **one at a time**, all others at production. The six:
`POOL QUALITY` (ORACLE best member of the set actually consumed), `POOL DIVERSITY` (within-set
dispersion and effective rank), `POOL COMMON-MODE BIAS` (`f = |ē|² / mean|e_k|²`, S23's exact
identity), `WITHIN-POOL RANKING` (in-band ρ and the deployable-argmin-vs-oracle-best gap), `READOUT`
(uniform average vs ORACLE convex over the SAME set), `RECONSTRUCTION` (chain − cloud, per arm).

**PREDICTION.** The six one-at-a-time idealisations sum to **less** than the cumulative 2.0966 Å,
and the largest single one is READOUT-given-a-perfect-quality-estimate, not POOL QUALITY.

**FALSIFIER.** If idealising POOL QUALITY alone (best member of the consumed 75, everything else
production) already exceeds the READOUT idealisation, generation outranks selection and rule 17's
burden is discharged in generation's favour.

**CONTROL.** Every idealisation is ORACLE and labelled. The no-op arm (production reproduced
exactly) must return 3.2105 ± the A2 floor.

**STATISTICAL RULE.** Descriptive with SE beside every mean; no verdict claimed from a decomposition.

**DEPLOYMENT CONDITION.** None. Diagnostic.

---

## H-P2a — The filter's ORACLE ladder, with best-of-K priced at every rung

**HYPOTHESIS.** Most of the "filter loss" on the best-member axis is an **order statistic**, not
retrieval skill, and the earliest rung (universe → 500 by BLOSUM) has never been priced.

**MECHANISM.** `best-of-K` falls with K for any ordering, including a random one. A filter has
*skill* only relative to a random subset of the same size drawn from the same parent set.

**PREDICTION (registered, directional).**
- `best(BLOSUM-500 of universe)` is **better** than the mean of random-500-of-universe draws
  (BLOSUM retrieval has real skill at the best-member axis).
- `best(score-128 of 500)` is **better** than random-128-of-500 draws.
- Both gaps are **smaller** than the raw rung losses, i.e. ≥ 40% of each rung loss is order
  statistic.

**FALSIFIER.** If score-128 ≈ random-128 (within the draw sd), the 500 → 128 filter has **zero**
best-member retention skill and the 0.4357 Å is entirely order statistic plus noise. If BLOSUM-500
is no better than random-500, retrieval itself is the earliest irreversible loss.

**CONTROL.** ≥ 5 independent draws at every random rung; report draw mean **and** draw-to-draw sd
(rule 10). Seeds derived from `(pdb, rung, draw)`, pinned in the artefact.

**STATISTICAL RULE.** `stats_lib.compare` with pinned folds for arm-to-arm; the random rungs are
reported as a distribution, never as their best draw.

**DEPLOYMENT CONDITION.** None — every arm here is ORACLE / NOT DEPLOYABLE.

---

## H-P2b — Can a native-free rule retain more REACHABLE quality? (the deployable arm)

**HYPOTHESIS.** The deployed terminal operator is a uniform average over 75 members. What such an
operator can spend is the **set mean** and the set's common mode; the set's *best member* is nearly
unreachable through it. Therefore an ORACLE best-member loss at the filter is **not** evidence that
the filter destroys reachable quality.

**MECHANISM.** `operator-consumes-set-mean` (S18): `d_out ≈ 1.16·d_set_mean + 0.04·d_set_best`.
That law was fitted under **random** gates and is recorded as BREAKING under score-based gates, so
it is **re-measured here**, not quoted.

**PREDICTION (registered).**
1. Across the P2b arms, a regression of built-chain `d_out` on (`set_mean`, `set_best`) gives a
   `set_mean` coefficient at least **5×** the `set_best` coefficient.
2. **No native-free set-selection rule beats production by ≥ 1 MDE on the built chain.**

**ARMS** (the 75-member set fed to the unchanged medoid-average + projection):
`PROD` score top-75 · `DEDUP` score top-75 over byte-distinct windows · `WIDE_CONS` score top-150
then the 75 nearest the top-150 medoid · `MMR` score with a diversity penalty (λ fixed a priori at
0.5, **not tuned**) · `BLOSUM75` the 75 highest BLOSUM sums · `RAND75` 3 draws from the 500.

**FALSIFIER.** Either prediction failing is the result: a ≥ 1 MDE native-free win would show the
filter *does* destroy reachable quality; a `set_best` coefficient comparable to `set_mean` would
overturn the reachability argument and re-open the best-member ladder as a deployable target.

**CONTROL.** `RAND75` carries its own draw distribution. All arms projected **in the same job** from
clouds built in that job (rule 3); the shared `PROD` rows are checked bit-identical against the
stored production record and the check is printed.

**STATISTICAL RULE.** Contract §1 in full. `< 0.7×` NOT A RESULT, `0.7–1.0×` NOT MEASURED.

**DEPLOYMENT CONDITION.** A rule is a deployment candidate only if ≥ 1 MDE better on the built
chain with ≥ 4/5 folds agreeing and no ORACLE input anywhere in its definition, including λ.

---

## H-P3 — The hard wall: what in-band skill IS, and what it costs

**HYPOTHESIS.** Positive in-band skill is not merely undiscovered. In the frame the readout actually
works in, the in-band ordering is **algebraically dominated by a term no pool-internal statistic can
observe**, and that term's unknown is exactly the common mode `μ`. The two bottlenecks S31 §20.3
left open — *"a per-candidate quality estimate"* and *"the pool's common mode"* — are **one
requirement**.

**MECHANISM (derivation, to be verified numerically).** In S31's own frame (all 128 superposed on
the uniform medoid, native superposed likewise), with `c` the pool centroid, `μ = c − t` and
`d_k = x_k − c`:

```
a_k  =  |x_k − t|²  =  |μ|²  +  2⟨μ, d_k⟩  +  |d_k|²
         ^^^^^^        ^^^^^^^^^^^^^^^^      ^^^^^^^
         constant      UNOBSERVABLE (U_k)    OBSERVABLE (V_k)
         in k          needs μ's direction   a pool statistic
```

`|μ|²` is constant across `k`, so it cancels from every within-pool ordering: **the entire in-band
problem is ordering `U_k + V_k`.** Every native-free signal ever tested here (consensus, typicality,
medoid criterion, dispersion) is a function of `V_k` alone. Two consequences:

1. **The negative in-band ρ of a typicality signal is forced.** The band is defined as the smallest
   `a_k`, i.e. by conditioning on `U + V`. Conditioning on a sum is a **collider**: inside the band
   `U` and `V` acquire a negative correlation, so a `V`-only signal must lose sign there. The
   magnitude is predicted by `Var(U)/Var(V)`.
2. **Perfect in-band ordering ⟺ knowing μ.** Forward: `μ` gives `U_k` and `V_k` is observable, so
   `a_k` follows exactly. Converse: `{a_k}` with observable `{d_k}` gives `3n + 1` unknowns
   (`μ`, `|μ|²`) against `K` linear equations `a_k − |d_k|² = |μ|² + 2⟨μ,d_k⟩`, over-determined at
   K = 128 ≫ 3n ≤ 48, so `μ` is **recoverable by least squares**.

**PREDICTIONS (registered, in order of how decisively each can fail).**

| id | prediction | falsifier |
|---|---|---|
| P3-1 | `Var(U)/Var(V)` over the top-128, median across targets, is **> 1** | ≤ 1 ⟹ the unobservable term is not dominant and the mechanism is wrong |
| P3-2 | The in-band Spearman of `−V` (the consensus direction) is **negative** on ≥ 80% of targets, and its mean lies within ±0.10 of the value predicted from the per-target `(U,V)` covariance inside the band | a positive mean, or a prediction error > 0.10, falsifies the collider account |
| P3-3 | Least-squares recovery of `μ` from `{a_k, d_k}` has relative residual **< 1e-6** on ≥ 120/126 targets | a larger residual falsifies the equivalence |
| P3-4 | **ZERO-μ CONTROL (ORACLE):** re-defining the label as `a'_k = |d_k|²` (i.e. the same pool with the common mode removed) flips the in-band ρ of the consensus signal to **> +0.5** | no flip ⟹ the sign is not caused by `μ` and the mechanism is wrong |
| P3-5 | **THE PRICE CURVE:** for `μ̂ = cos·μ̂_dir + noise` at controlled `cos(μ̂, μ)`, in-band ρ of `2⟨μ̂,d⟩ + V` crosses zero at some `cos* > 0`; `cos*` is reported as **the price of in-band skill** | if in-band ρ is positive even at `cos = 0`, the requirement is not μ |
| P3-6 | **THE DEPLOYABLE TEST:** the distogram's implied correction direction `μ̂_DIS` has `cos(μ̂_DIS, μ) < cos*`, i.e. the only external channel in the pipeline is **below the price** | `cos(μ̂_DIS, μ) > cos*` ⟹ a native-free signal WITH positive in-band skill exists and must be carried to the endpoint |

**CONTROLS.** (a) random-direction `μ_rand` matched in **norm** to `μ` (rule 6: matched to the
operator's own space) — must reproduce the same negative in-band ρ, showing the sign is driven by
`|μ|` and not by any structure in `μ`; (b) the oracle signal `U+V` must give in-band ρ = 1 (a
verification that can fail, rule 5); (c) a shuffled-`d` control for the recovery in P3-3.

**STATISTICAL RULE.** Per-target ρ aggregated with fold-clustered CIs; any claim of positive in-band
skill needs the fold CI to exclude zero **and** ≥ 4/5 folds agreeing, and must then be carried to the
built chain before it counts.

**DEPLOYMENT CONDITION.** P3-6 firing is the only deployable outcome in H-P3; it would be carried
through the convex readout to the built chain in a single paired job.

---

## H-P4 — Generation, only against the 2.10 Å bar (conditional arm)

**HYPOTHESIS.** Adding candidates cannot help while the terminal operator consumes the set mean and
the pool's error is 68% common-mode: new candidates drawn from the same library inherit the same
common mode, so they raise diversity without moving `μ`.

**MECHANISM.** Averaging removes only the idiosyncratic 32%. A wider or more diverse pool changes
`mean|d_k|²` and leaves `|μ|²` where it is.

**PREDICTION.** Widening the retrieval (K = 500 → 2000, same BLOSUM key) lowers the ORACLE best
member **substantially** (an order statistic) and leaves the common-mode fraction's numerator
`|μ|²` statistically unchanged, and therefore leaves the deployable endpoint unchanged or worse.

**FALSIFIER.** `|μ|²` falling with K by more than its draw sd would show generation moves the common
mode and would discharge rule 17's burden.

**CONTROL.** Matched-K random subsets; `|μ|²` measured on equal-size sets so the comparison is not
confounded by set size.

**STATISTICAL RULE / DEPLOYMENT.** As H-P2b. This arm runs only if H-P2b or H-P3 leaves room for it.

---

## H-P5 — Strata: mechanism only, filter-independent first

FAIL18 is outcome-defined (contract §12). Every stratified statement is computed on
**`worst18_poolmean`** and **`worst18_bestpool`** first (filter-independent), with FAIL18 reported
beside them as a labelled diagnostic. **PREDICTION:** `Var(U)/Var(V)` is **higher** on the
filter-independent tail than on the other 108 — i.e. the tail is where the unobservable term
dominates most, which is S31's *"the tail's pools are displaced together, not less varied"*
(`S/B` 0.776 → 0.423 at unchanged `n_distinct`) restated as a measurable quantity.
**FALSIFIER:** equal or lower `Var(U)/Var(V)` on the tail.

---

## Multiplicity

Every emitted comparison is appended to `s32/MULTIPLICITY.md` as it is emitted, marked
REGISTERED (an arm named above) or EXPLORATORY. The registered count at commit time is:
H-P1 six descriptive quantities (no verdicts), H-P2a four rungs × (1 ordered + 5 random draws),
H-P2b six arms vs PROD, H-P3 six predictions + three controls, H-P4 two arms, H-P5 two strata × one
quantity.

## Artefacts

`s32/results/s32_P_*.json` and `s32/results/s32_P_*_rows.jsonl`, stable keys, written through
`stats_lib.save_atomic` where a row count exists.
