# Sprint 12 — ADVERSARIAL AUDIT

Role: try to destroy the sprint's conclusions. A claim that survives here is worth more
than one that has not met this file. Everything below is a measurement, not an opinion.

Code: `s12/adv_*.py`. Results: `s12/results/adv_*.json`.

---

## A. THE INSTRUMENT (`s12/instrument.py`) — audited FIRST, because a defect here
## contaminates every agent at once.

Script: `s12/adv_instr.py`, `s12/adv_proj.py`. Results: `s12/results/adv_instrument.json`,
`s12/results/adv_projection.json`.

### A1. `shipped_score` vs production (`core/pipeline.py:743` → `s7/debias.py:242`)

The instrument stores `risk` as **float32** (`instrument.py:176`) and gathers/means in
float32; production upcasts to float64 (`core/pipeline.py:650`
`np.asarray(d._risk, float)`). I recomputed the whole 126-target filter under three
precisions and compared the top-75 set to the production record's `sub`.

| path | targets whose top-75 == production `sub` | mean members differing |
|---|---|---|
| float64 risk (production) | **126 / 126** | 0.00 |
| float32 risk (instrument `shipped_score`) | **126 / 126** | 0.00 |
| float32 risk + float32 round-trip on D (`selfcheck` line 238) | **126 / 126** | 0.00 |

- max |score(f32) − score(f64)| over all 126×500 = **3.4e-06**; argmin moves on **0/126**.
- shipped argmin mean: 3.4540004952559396 under *both* precisions. top-75 best 2.3061526
  under production `sub`, my f64 recompute and my f32 recompute — **identical to 7 d.p.**

**VERDICT: no defect.** The float32/float64 discrepancy exists but is provably inert on
this instrument. Every agent's filter is the production filter.

**Side finding worth knowing (tie-break exposure).** The shipped score has on average
**35.2 exactly-tied values per 500-member pool** (duplicate windows retrieved through two
parents; production reports `n_distinct` = 66 of 75 on 1A13). Production resolves them
with `argsort(kind="stable")`, i.e. by BLOSUM/universe position — deployable and correct.
But any agent that ran `np.argmin`/`np.argsort` on a signal derived from a tied score and
then read an oracle-sorted array has the exact leak the project record already burned once
(the 1.386 Å phantom). Audited per-agent in §E.

### A2. `coordinate_average` vs production `avg_ca` (`core/pipeline.py:939`)

126/126 targets: max CA-RMSD between `I.coordinate_average(u["W"][pool][sub])[0]` and the
production `avg_ca` = **2.13e-07 Å**, mean 4.6e-08 Å. Medoid choice agrees on all 126
(a medoid disagreement would show as ~1 Å, not 1e-7). **VERDICT: no defect.**

### A3. `project` vs production (`core/pipeline.py:955`) — ONE REAL, SMALL DEFECT

8 targets sampled (seed 0). Instrument vs production record:

| arm | max CA-RMSD(mine, production) | max |Δ RMSD-to-native| |
|---|---|---|
| `fit_ca` (λ=0 — this is the 3.204 synthesis arm everyone quotes) | **7.0e-04 Å** | **1.2e-04 Å** |
| `ca` (λ=0.3 — the shipped arm) | **0.103 Å** (1FUV) | **0.021 Å** (2NB7) |

Cause is documented in `core/pipeline.py:145-152`: the λ projection is **degenerate** — a
CA trace admits two ideal-geometry torsion branches, and *"a 1e-13 Å difference in the
forward map routes some L-BFGS-B trajectories into the other torsion branch"*. The
instrument's input `C` differs from production's by ~1e-7 Å (production averages the
float32-round-tripped `pool["W"]`, `core/pipeline.py:795`; the instrument averages the
float64 upcast of the same float32 npz, and `cc.superpose_batch` vs `I.superpose_batch`
differ in the last ULP). That 1e-7 amplifies to 0.10 Å of structure.

**Practical consequence and the number to use: any result computed through
`I.project(...)["ca"]` (the λ=0.3 arm) carries an irreducible ~0.02 Å per-target
reproducibility floor. The brief's own ignore-threshold is 0.03 Å, so this does not
invalidate anything, but a λ-arm effect reported at 0.02–0.05 Å is inside the noise of
the tool that measured it.** The λ=0 arm (`fit_ca`), which is what almost every agent
actually uses, reproduces to 1e-4 Å and is safe.

### A4. Are `s12/cache/disto_*.npz` really the leave-fold-out models for the right folds?

This is checked *for free and with maximum power* by A1: the cached distogram of every one
of the 126 targets reproduces production's top-75 **exactly, on all 126 targets**. A cache
poisoned with the wrong fold's model (or another target's distogram) could not select the
same 75 of 500 windows on a single target, let alone 126. `risk.shape[0] == len(triu(n,2))`
on 126/126. `I.pair_index(n, min_sep=2)` is bit-identical to `core/predict.py:102`.
**VERDICT: no defect.**

### A5. Instrument summary

| check | verdict |
|---|---|
| `shipped_score` ≡ production filter | PASS (126/126 exact) |
| `coordinate_average` ≡ production | PASS (2e-7 Å) |
| `project`, λ=0 arm | PASS (7e-4 Å) |
| `project`, λ=0.3 arm | **~0.02 Å noise floor** — usable, but quote no λ-arm effect below 0.03 Å |
| distogram cache = right LFO fold model | PASS (implied at full strength by A1) |
| pair index / min_sep | PASS |
| float32 round-trip placement | harmless (measured, not assumed) |

**The instrument is not the problem. No agent's numbers are contaminated by it.**

---

## B. CLAIM 3 — the fibril/lasso enrichment (10/18 vs 6/108, Fisher p = 1.2e-06)
## **SURVIVES. I attacked it three ways and it did not move.**

Script `s12/adv_meta.py` -> `s12/results/adv_meta.json`.

### B1. Multiplicity — the composite was chosen post hoc, so I priced the whole search space

`fibril|lasso` is a UNION of two flags picked out of the 11-flag family in
`s12/forensics_part1.py:151` (lasso, cyclic, fibril, xray, membrane, cosolvent, bound,
designed, covalent, nonaqueous, any_flag), plus 7 diagnostic flags at `:154` and a further
15 in `s12/fail_headers.py:98-108`. The honest family for a post-hoc *pair union* is
11 singles + all 55 pairwise unions = **66 hypotheses**. I enumerated all 66 and corrected:

| correction | p |
|---|---|
| reported (uncorrected) | 1.22e-06 |
| Bonferroni over the 11 singles | 1.35e-05 |
| **Bonferroni over all 66 (11 singles + 55 unions)** | **8.08e-05** |
| **max-statistic permutation over the same 66, 4000 label shuffles** | **0 / 4000, p_FWER < 2.5e-04** |

`fibril|lasso` is also the family *maximum* — the most enriched of all 66 — which is exactly
the statistic the permutation null controls. It clears every correction by >= 2 orders of
magnitude. Runners-up: `fibril` alone 5.5e-04 (6/18 vs 4/108), `lasso|xray` 1.2e-03,
`lasso` alone 3.8e-03 (4/18 vs 2/108). The union beats either half, and both halves are
individually significant — this is not one flag carrying a passenger.

### B2. Leave-one-out — not one target, not one class

Dropping each of the 126 targets in turn: **the worst p over all 126 drops is 5.38e-06**
(dropping 2BFI). Dropping any one of the 10 flagged FAIL18 members gives p = 5.38e-06 in
every case. There is no influential point.

### B3. I re-read the headers myself and re-derived the labels independently

`s12/adv_meta.py:26-56` implements my own regexes over `pdbs_ext/`+`pdbs/`
HEADER/TITLE/COMPND/KEYWDS/EXPDTA and REMARK 210/245/280. Diffed against
`fail_headers.json`:

- **fibril: 0 disagreements / 126. lasso: 0 / 126. cyclic: 0 / 126.**
- 16 disagreements, all in `membrane` (11) and `bound` (5), all attributable to two
  deliberate differences in my own regexes (I removed TFE/trifluoroethanol from `membrane`
  — a cosolvent is not a membrane mimic — and added ANTAGONIST|INHIBITOR to `bound`).
  Neither flag is part of the claim.

Reading the 18 titles by hand I agree with every fibril and lasso call: 2BFI/3SGO/5W52/9L1M
are explicitly amyloid or steric-zipper segments; 2BP4 is A-beta(1-16) in 80 % TFE and 2JN5 an
alpha-synuclein dodecapeptide (both amyloid by keyword); 2N5C/7JS6/7LCW/9KAR are all titled
LASSO PEPTIDE with a LINK record for the macrolactam.

### B4. The honest caveat that should travel with the claim

The flag is neither necessary nor sufficient: **6 of the 16 flagged targets do not fail, and
8 of the 18 failures are not flagged.** 10/16 = 62.5 % failure rate (exact 95 % CI
[0.35, 0.85]) against 8/110 = 7.3 % ([0.032, 0.138]). The intervals are far apart, but the
flagged group is n = 16 and the point estimate is soft.

---

## C. CLAIM 2 — is the FAIL18 a stable object, or a threshold artefact?
## **SPLIT VERDICT. Stable under the relative-band family; the enrichment DIES under an
## absolute-band definition; and the headline "they return 6.0 A" is near-tautological.**

Script `s12/adv_fail18.py` -> `s12/results/adv_fail18.json`. 99 grid cells:
BAND in {0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0} x K in {250,500,1000} x M in {50,75,100}, plus
an ABSOLUTE band (rr <= 2.0/2.5/3.0 A, which never references `pool_best`).

**There is no random seed anywhere in this chain** — retrieval is
`argsort(-sim, kind="stable")` (`core/pipeline.py:711`) and the filter a deterministic
gather — so "does it survive a fresh seed" has no content here; the only instability
available is in the three constants. Verified by recomputing the whole selfcheck twice:
bit-identical.

### C1. Membership under the constants

| definition | n | Jaccard vs FAIL18 | kept of 18 | added |
|---|---|---|---|---|
| **band 1.5, K 500, M 75 (production)** | 18 | 1.000 | 18 | 0 |
| band 1.25 | 20 | 0.900 | **18** | 2 |
| band 1.75 | 15 | 0.833 | 15 | **0** |
| band 1.0 | 24 | 0.750 | **18** | 6 |
| band 2.0 | 13 | 0.722 | 13 | **0** |
| band 0.75 | 34 | 0.529 | **18** | 16 |
| M = 50 | 19 | 0.947 | **18** | 1 |
| M = 100 | 16 | 0.889 | 16 | **0** |
| K = 250 | 10 | 0.556 | 10 | **0** |
| K = 1000 | 21 | 0.857 | **18** | 3 |
| **absolute band 2.5 A** | 26 | **0.294** | **10** | **16** |
| **absolute band 2.0 A** | 18 | **0.241** | **7** | **11** |
| absolute band 3.0 A | 22 | 0.481 | 13 | 9 |

**The relative-band family is NESTED, not swapping.** In every cell of the relative family
the set is a subset or a superset of FAIL18 — not once does a perturbation eject a core
member *and* introduce a stranger. That is much stronger stability than the mean Jaccard
(0.498, dragged down by the extreme cells and the absolute family) suggests, and it is the
property that matters: a "FAIL18 vs other-108" contrast cannot be a swap artefact.

**And the criterion is not a knife-edge cut through a continuum.** The distribution of
recall (band members kept by the top-75) is *bimodal*: 18 targets at 0, **only 4 targets at
1-3**, 104 at >= 4 and 94 at >= 10. There is a real gap at the threshold.

### C2. Shrinkage of the subgroup contrasts

| definition | n | emitted gap (A) | rho gap | fibril/lasso Fisher p |
|---|---|---|---|---|
| **production (band 1.5, K500, M75)** | 18 | **+3.292** | **-0.538** | **1.22e-06** |
| band 1.25 | 20 | +3.246 | -0.575 | 4.66e-06 |
| band 1.75 | 15 | +3.277 | -0.568 | 4.23e-05 |
| band 1.0 | 24 | +2.831 | -0.512 | 4.01e-05 |
| band 2.0 | 13 | +3.416 | -0.625 | 9.62e-06 |
| band 0.75 | 34 | +2.446 | -0.519 | 2.32e-05 |
| M = 50 | 19 | +3.209 | -0.559 | 2.45e-06 |
| M = 100 | 16 | +3.340 | -0.561 | 5.35e-06 |
| K = 250 | 10 | +2.950 | -0.521 | 2.53e-04 |
| K = 1000 | 21 | +2.938 | -0.500 | 8.45e-06 |
| **absolute 2.5 A** | 26 | **+1.775** | **-0.425** | **0.097 (n.s.)** |
| **absolute 2.0 A** | 18 | **+1.103** | **-0.233** | **0.244 (n.s.)** |
| absolute 3.0 A | 22 | +2.586 | -0.552 | 1.3e-03 |

Over all 93 evaluable cells the emitted gap runs 1.01-4.47 A (mean 2.72) and the
fibril/lasso p is < 0.05 in **83 %** of them.

**WEAKENING, stated precisely.** Every "FAIL18 vs other-108" contrast in the sprint is a
contrast about *filter* failure measured **relative to each target's own pool**. Re-pose the
same question in absolute terms — "did the filter discard everything within 2.5 A of native?"
— and you get a 26-target set sharing only 10 members with FAIL18, whose emitted gap has
shrunk by 46 % and **whose fibril/lasso enrichment is no longer significant (p = 0.097; at
2.0 A, p = 0.244)**. The reason is structural, not statistical: the relative band conditions
on `pool_best`, so it selects targets the *filter* failed; the absolute band mixes in targets
whose *pool* failed, and those are a different population. Both questions are legitimate. The
sprint has answered only the first, and its causal story ("the 18 are peptides whose deposited
conformation is not the free-solution conformation") is a story about the first.

### C3. The circularity, quantified

Zero recall is *logically equivalent* to `band_score_pct_min > M/K = 0.15`. So the reported
`o_band_score_pct_min` contrast (**0.391** on FAIL18 vs 0.009 on the other 108,
`s12/forensics_FINDINGS.md:69`) has a hard definitional floor of 0.15 built into the FAIL18
side: even a FAIL18 that only *just* fails must score >= 0.15. **Only 0.391 - 0.15 = 0.241 of
the 0.382 gap (63 %) is non-definitional.** The same applies, more weakly, to
`o_band_score_pct` (0.687 vs 0.287) and to the headline "**the 18 return ~6.0 A**": zero
recall forces all 75 averaged members to have rr > pool_best + 1.5 = 3.78 A on this group, so
a large emitted RMSD is very nearly entailed by the definition and carries little independent
information. Cite it as a description of the group, never as a discovery about it.

The control the sprint never ran is the **low-recall group** (1-3 band members kept), subject
to almost the same force but not in FAIL18:

| group | n | band_pct_mean | band_pct_min | emitted | rho(score, rr) | fibril/lasso |
|---|---|---|---|---|---|---|
| zero recall (FAIL18) | 18 | 0.687 | 0.391 | 6.026 | +0.107 | 0.556 |
| **low recall (1-3)** | **4** | **0.506** | **0.061** | **3.945** | **+0.145** | **0.000** |
| high recall (>=4) | 104 | 0.278 | 0.007 | 2.687 | +0.664 | 0.058 |

The low-recall group has FAIL18-like **rho** (+0.145 vs +0.107) with *none* of the
lasso/fibril character and an intermediate emitted RMSD. Consistent with the sprint's story
(anti-ranking is continuous; the out-of-distribution-reference cause is what pushes a target
all the way to zero) but n = 4 and it settles nothing on its own.

**The findings that DO carry independent information about the 18** — because the definition
does not constrain them — are: the header metadata enrichment (B, survives), the native
strand fraction (0.470 vs 0.166), and the *negative* result that BLOSUM retrieval is equally
good on both groups (in-pool band fraction 0.079 vs 0.073, p = 0.96). Those stand.

---

## D. THE CENTRAL CLAIM, PART 1c — "generation is wasted because the OBJECTIVE cannot rank"
## **REFRAMED, and one of its four evidentiary legs is REFUTED. The binding constraint at
## the terminal stage is the OPERATOR'S SET-AVERAGING, not the objective's ranking.**

Script `s12/adv_operator.py` -> `s12/results/adv_operator.json`. All 126 targets.

### D1. The design

I gave the objective a *perfect ranking oracle* on the production pool and asked what each
terminal operator does with it. The top-75 is built as the shipped top-75 with its `j`
worst-scoring members replaced by the pool's `j` ORACLE-BEST members. `j = 0` is exactly
production; `j >= 1` means the pool's single best structure is present *and at rank 1* —
which is what a perfect objective would have handed the operator, on the same candidate set.

ORACLE/DIAGNOSTIC: `rr` chooses the inserted members. This is an attribution experiment,
not a method.

### D2. The ladder

| j | set best | set mean | **argmin emits** | **avg emits** | projected `fit` |
|---|---|---|---|---|---|
| **0 (production)** | 2.3062 | 3.5507 | **3.4540** | **3.0483** | **3.2052** |
| 1 | **1.7108** | 3.5327 | **1.7108** | 3.0230 | 3.1761 |
| 2 | 1.7108 | 3.5159 | 1.7108 | 3.0004 | — |
| 5 | 1.7108 | 3.4698 | 1.7108 | 2.9406 | 3.1033 |
| 10 | 1.7108 | 3.4016 | 1.7108 | 2.8497 | — |
| 25 | 1.7108 | 3.2192 | 1.7108 | 2.5801 | — |
| 50 | 1.7108 | 2.9374 | 1.7108 | 2.1953 | — |
| 75 (perfect filter) | 1.7108 | 2.6855 | 1.7108 | **1.9630** | **2.1018** |

Cross-check: `j=75` avg = **1.9630** reproduces `forensics_FINDINGS.md:302`'s "perfect
FILTER, m=75 = 1.963" **exactly**, independently derived. The instrument and forensics agree.

### D3. The result, in one line

| operator | value of a PERFECT rank-1 decision (j=0 -> j=1) |
|---|---|
| argmin (the legacy arm) | **-1.7432 A** |
| **coordinate average -> project (what the system emits)** | **-0.0291 A** [-0.0382, -0.0211], 103W/20L |

The set now contains a structure 0.595 A better than anything it had, *at rank 1*, and the
deployed terminal operator converts that into **0.029 A** — below the brief's own 0.03 A
ignore threshold, and 1/60th of what the same information is worth to argmin.

### D4. The law: the operator's output is a function of the set MEAN, not the set BEST

Pooled OLS over all 882 (target x j) perturbations, in per-target differences from j=0:

    d_avg = 0.0385 * d_set_best  +  1.1620 * d_set_mean  -  0.0277        R^2 = 0.893

| model | R^2 | slope |
|---|---|---|
| **set_mean alone** | **0.891** | **1.188** |
| set_best alone | 0.187 | 0.333 |

Adding `set_best` to a model that already has `set_mean` moves R^2 by 0.002 and gives it a
coefficient of 0.039. And the ratio d_avg / d_set_mean is nearly constant down the whole
ladder — 1.409, 1.378, 1.332, 1.332, 1.412, 1.391, 1.254 — i.e. the operator is a
near-linear map of the set mean with gain ~1.2-1.4.

The "fixed quantile" hypothesis is confirmed for realistic perturbations: the emitted
structure lands at quantile **q = 0.249, 0.224, 0.213, 0.196, 0.201, 0.230** of its own
input set's RMSD distribution for j = 0, 1, 2, 5, 10, 25 (it only breaks to 0.151 / 0.068
at j = 50 / 75). The operator emits something better than ~78 % of what it is given,
**essentially regardless of what it is given**, until you replace the majority of the set.

### D5. What this does to the sprint's central claim

The claim assembled by the coordinator is: *the binding constraint is the objective's lack
of pair-specific information, and every generation-side improvement is wasted because the
objective cannot rank within the improved set.* Four legs were offered. This measurement
breaks the load-bearing one:

* **REFUTED as evidence:** the assembly agent's "the deployable assembly improved the oracle
  pool floor 1.711 -> 1.175 while the emitted answer got WORSE by +0.218" does **not**
  demonstrate that the objective cannot rank. **A perfect ranking would not have helped
  either**, because improving a set's *floor* is not the currency the terminal operator
  spends. Under the fitted law, a 0.536 A improvement in the set's best member with the set
  mean unchanged buys 0.536 x 0.039 = **0.021 A**. The observed asm result (+0.218 A) is
  fully explained by the assembly top-75's set mean being ~0.18 A worse than the pool's,
  with no ranking claim required.
* **SURVIVES:** the objective agent's r_sep and iso-MAE work (see F), and the aggregation
  agent's "the operator is at its native-free optimum" — that claim is about the operator's
  *form* at a fixed set, and my result is about the *set*, so they are consistent, not in
  tension. `agg_FINDINGS.md:110` already says selection is worth 0.61 A of the 1.11 A oracle
  headroom on the top-75; my law is the mechanism behind that number.
* **REFRAMED:** the correct statement is *the terminal operator consumes a 75-member mean,
  so any generation- or ranking-side improvement that does not move the BULK of the top-75 is
  invisible to it by construction.* That is a statement about the architecture, not about the
  distogram's information content, and it changes what to optimise: **the objective's job in
  this pipeline is to raise the mean quality of the 75, not to find the best one.**
* **Practical corollary the sprint should act on:** the system already contains a second
  operator (argmin) that *does* cash ranking improvements at 1.74 A per perfect top-1. A
  ranking improvement worth ~0 A through the average can be worth ~1 A through argmin. Any
  future ranking work must be scored through BOTH terminals or it will be mis-priced.

---

## E. THE CENTRAL CLAIM, PART 1b — "the objective's own optimum is 4.62 A, worse than the
## 3.20 A the pipeline emits"
## **The NUMBER IS REFUTED (4.62 -> 3.53). The CONCLUSION SURVIVES, and my attack makes it
## much stronger than the version that was reported.**

Script `s12/adv_decode.py` -> `s12/results/adv_decode.json`. All 126 targets, ~31 min.

4.62 A is a **classical-MDS** embedding of `E[d]` (`s12/obj_errstruct.py`). Classical MDS is
a weak decoder: it eigendecomposes a doubly-centred squared-distance matrix that need not be
Euclidean, it ignores the score's shell weights, and it minimises strain on `E[d]` rather
than the shipped Bayes-risk score the pipeline actually ranks with. I ran four decoders:

| decoder | mean shipped SCORE | mean CA-RMSD | FAIL18 | other-108 |
|---|---|---|---|---|
| **A. classical MDS on E[d] (as published)** | 1.6706 | **4.621** | 6.515 | 4.305 |
| B. classical MDS -> the production projection | 1.4627 | 4.033 | 5.953 | 3.713 |
| C. SMACOF stress majorisation, shell-weighted, multi-start | 1.4267 | 4.243 | 6.265 | 3.906 |
| D. SMACOF -> the production projection | **1.3662** | 4.026 | 6.089 | 3.682 |
| **E. direct L-BFGS argmin of the SHIPPED SCORE over torsion space, 8 starts** | 1.4205 | **3.532** | 6.008 | 3.120 |
| the pool's best-scoring member (argmin arm) | 1.4238 | 3.454 | 6.008 | 3.028 |
| **what the pipeline emits (`fit_ca`)** | 1.5089 | **3.204** | 6.026 | 2.734 |
| **the NATIVE structure** | **2.0900** | **0.000** | — | — |

Arm A reproduces the published 4.621 exactly.

### E1. The number is wrong by 1.09 A

Arm E — the honest argmin of the actual objective inside the class of structures the
pipeline can emit — reaches **3.532 A**, not 4.621. Paired: **-1.089 A [-1.256, -0.923],
108W/18L, drop-top-10 -0.932**. *"The objective's own geometric optimum is 4.62 A"* should
not be cited; the correct figure is 3.53 A.

### E2. The conclusion survives, at a smaller margin

The pipeline still beats the honest argmin of its own objective: **+0.328 A
[+0.213, +0.448], 38W/88L, drop-top-10 +0.419**. The claimed margin was 1.42 A; the true
margin is 0.33 A. It is real, it is robust to dropping the top 10, and it is not close to
zero — but it is a quarter of what was reported.

Notably, free torsion-space optimisation of the objective is worth nothing over simply
picking the pool's best-scoring window: **+0.078 A [-0.007, +0.160]**. Ten thousand
gradient steps in a 2n-dimensional space buy no accuracy over a library lookup.

### E3. My attack STRENGTHENS the underlying point, on the score axis

Order the seven decoders plus the native by the score they attain (lower = better under the
objective) and read off their RMSD:

| rank on the OBJECTIVE | arm | score | RMSD |
|---|---|---|---|
| 1 (best) | SMACOF -> project | 1.3662 | 4.026 |
| 2 | torsion argmin | 1.4205 | 3.532 |
| 3 | pool argmin | 1.4238 | 3.454 |
| 4 | SMACOF | 1.4267 | 4.243 |
| 5 | MDS -> project | 1.4627 | 4.033 |
| 6 | **the pipeline's emitted structure** | 1.5089 | **3.204** |
| 7 | MDS | 1.6706 | 4.621 |
| 8 (worst) | **the NATIVE** | **2.0900** | **0.000** |

The arm that minimises the objective best (SMACOF->project, 1.3662) is the second-worst
structure; the arm that scores worst (the native, 2.0900) is the correct answer. **The
native scores worse than the emitted structure on 117/126 targets and worse than the
torsion argmin on 122/126.** Within-target Spearman(arm score, arm RMSD) across these eight
structure classes is **-0.042** (negative on 69/126); excluding the native it is +0.262.

This is the cleanest demonstration in the sprint that *minimising the shipped score harder
produces worse structures*, and it is now established against a real optimiser rather than
against a weak decoder — which is exactly the objection the 4.62 A number was vulnerable to.

---

## F. THE CENTRAL CLAIM, PART 1a — "r_sep = 0.196: the distogram carries almost no
## pair-specific information"
## **The NUMBER is estimator-dependent by a factor of 2 (0.196 -> 0.37-0.39 under a pooled,
## cross-validated partialling). The CONTRAST it rests on SURVIVES every estimator.**

Script `s12/adv_rsep.py` -> `s12/results/adv_rsep.json`. All 126 targets.

`s12/obj_errstruct.py:84-90` partials out separation by subtracting **per-shell means
within each target**, then averages the 126 within-target correlations. For a length-n
target, shell s holds n-s pairs: the longest shells hold 2 pairs and 1 pair, and a 1-pair
shell contributes residuals of exactly (0, 0) to both vectors. Mean 67.8 pairs across 11.0
shells = **17.3 % of the degrees of freedom consumed by the partialling**. Six estimators
of the same quantity, with the matched i.i.d. null (the arm the claim is contrasted against)
computed identically:

| estimator | distogram | i.i.d. null at the same MAE | ratio |
|---|---|---|---|
| **E1 per-shell means, all shells (as published)** | **0.1959** | 0.4803 | 0.41 |
| E2 drop shells with < 3 pairs | 0.2033 | 0.4820 | 0.42 |
| E3 drop shells with < 5 pairs | 0.2093 | 0.4738 | 0.44 |
| **E5 pooled across targets, LFO regression on |i-j|** | **0.3675** | 0.7147 | 0.51 |
| **E6 pooled, within-target z-scored first, LFO** | **0.3909** | 0.6508 | 0.60 |

E1 reproduces the published 0.196 to 4 d.p. (0.1959).

**Finding 1 — the small-shell degeneracy is NOT the artefact.** Dropping every shell with
fewer than 3 (or 5) pairs moves r_sep by **+0.007 / +0.013**, i.e. slightly *up*. The
concern that per-shell means on 1-2 pairs absorb real signal is measurable and it is worth
0.01, not 0.2.

**Finding 2 — the estimator choice is worth a factor of two, and this matters for how the
number is quoted.** A pooled, leave-fold-out partialling on a flexible basis of |i-j|
(one-hot shells + n + sep/n + log sep) gives **0.367-0.391**, not 0.196. The reason is
substantive, not technical: a *within-target* partialling removes each target's own shell
profile, which is precisely the thing the distogram predicts reasonably well; a *pooled*
partialling removes only the corpus-average profile and leaves the per-target compactness
signal in. **0.37 is not "almost no pair-specific information".** The obj agent is aware of
both regimes and says so (`obj_FINDINGS.md:129-131`: "across targets the model does track
overall compactness; within a target it barely orders the pairs"), and my E5 (0.367)
recovers the record's pooled slope of 0.376 — so this is a *quoting* hazard, not an error.

**Finding 3 — the within-target estimate is under-powered per target but fine in aggregate.**
dof-corrected t-statistics: mean |t| = 1.97, and only **37 % of the 126 targets** reach
|t| > 2. Any per-target use of r_sep is noise; the 126-target mean is not.

**Finding 4 — the claim survives.** What the argument actually needs is that the distogram
retains *less* within-shell signal than a matched-MAE corruption of the native matrix. It
does, under **every** estimator: 41 %, 42 %, 44 %, 51 %, 60 % of the null's value. The
published framing ("2-3x more within-shell signal at the same MAE") becomes "1.7-2.4x" under
the pooled estimator. Weaker, same sign, same conclusion. And nothing in the obj agent's
causal chain depends on r_sep's absolute value: the load-bearing results are the iso-MAE
sweep and the pair-shuffle (`obj_FINDINGS.md` s9), which are measured on emitted RMSD and
are untouched by this.

---

## G. CLAIM 4 — the coordinator's own audits
## **G1 and G2 SURVIVE. I found no defect.**

### G1. `s12/coord_contain.py` — the 4 verbatim training-set leaks

Verified independently. The containment cache and the fold assignment are consistent, and
the four pairs (`1CEK ⊂ 1A11`, `2FBU ⊂ 2LMF`, `2P5H ⊂ 2P5J`, `6B9K ⊂ 1U6V`) reproduce.
The coordinator's own reading is the right one and it is already the conservative one: the
FAIL18 have *lower* training containment than the rest (0.301 vs 0.417), i.e. the leak runs
*against* the failure story, and rho(containment, emitted) = -0.129 (p = 0.15) on n = 126.
I add one point the coordinator did not make explicitly: the containment audit's power is
capped by C2's own null (a random sequence scores 0.56-0.63 containment against the fragment
bank), so the "74/126 at >= 0.6 containment" figure is **not interpretable** and only the
verbatim-substring subset (4/126, p ~ 0 under the null at n >= 9) is evidence. The
coordinator states this in C2 and then still reports the 74 in the C1 table; the two should
be merged so the 74 is never read as a leak count.

### G2. `s12/coord_benchN.py` — "no adequately-powered fresh benchmark exists"
### **SURVIVES. The gates are not too strict; I checked all 93 rejections individually.**

Script `s12/adv_gate.py` -> `s12/results/adv_gate.json`. I re-parsed every rejected PDB
directly from the text, independently of `geo.native_coords_from_pdb`, enumerating every
chain in MODEL 1 and its standard-residue CA count.

| rejection reason | n | my independent verdict |
|---|---|---|
| **length** | **67** | **CORRECT on all 67.** Every one is multi-chain, the parser took the first chain, and **every chain in every one of the 67 files has the same length and it is 1-8 residues**: 30 files of 6-mers, 15 of 8-mers, 10 of 7-mers, 5 of 5-mers, 3 of 4-mers, 2 of 1-mers, 1 of 3-mers, 1 of 2-mers. **Not one of the 67 contains any chain of length 9-16.** No wrong-chain selection, no truncation. |
| **CA step** | **16** | **CORRECT on all 16.** All are 9-12mers, but every one has **non-contiguous residue numbering** — genuinely missing residues. Examples: 4E0M chain A is resseq 1-7, 9, 10, 13 (16.07 A CA-CA jump); 6UG3 is 1,2,4,6,8,10,12 (three gaps); 8GJC has breaks at 8-10 and 12-14. A structure with missing interior residues cannot serve as a contiguous 9-16mer benchmark native. |
| sequence already in corpus | 8 | correct dedup |
| parse error (1S4A) / rebuild 2.67 (9RD7) | 2 | correct |

**The only defensible relaxation I can find** is to take the longest *contiguous* run from a
gapped chain. I computed it for all 16: it reaches 9-16 residues in exactly **three** cases
(6B17 res 2-12 = 11; 3SGM res 3-11 = 9; 6B79 res 3-11 = 9). That takes the structural-gate
pass from 40 to 43, and after the unchanged dedup / leakage / one-per-cluster stages would
add at most 2-3 to the final 16. **It does not reach 40, and it is not close.** The
coordinator's conclusion — no fresh, adequately-powered, distribution-matched benchmark
exists — stands, and the three validation options it names (dev24, a length-extrapolation
benchmark, declared nested CV) remain the only ones.

---

## H. THE ASSEMBLY LEG, FINISHED — the measurement `asm_deploy.py` did not make, and the
## control that kills it.
## **The assembly candidate set offers NO real headroom: a bank of pure Ramachandran noise
## of the same size BEATS it. asm's reading is correct in form and vacuous in content.**

Scripts `s12/adv_asmoracle.py`, `s12/adv_asmnull.py` -> `s12/results/adv_asmoracle.json`,
`s12/results/adv_asmnull.json`. 24 targets (the first 24 by (fold, pdb) — **none of them
FAIL18**, which is a real limitation of this sub-experiment).

### H1. The missing measurement: does assembly improve the achievable SET, or only the floor?

`asm_FINDINGS.md:305` reports that assembly improves the oracle *best member* 1.711 -> 1.175
while the emitted answer worsens by +0.218, and reads this as "the objective cannot rank".
Section D shows the floor is not the currency the operator spends, so the reading needed a
control: what does the assembly bank's **best possible 75-member set** emit?

| arm (24 targets) | set best | set mean | emits (avg) | emits (projected `fit`) |
|---|---|---|---|---|
| pool, shipped-score top-75 | 2.041 | 3.396 | 2.788 | 2.995 |
| pool, ORACLE top-75 | 1.824 | 2.732 | 1.930 | 2.073 |
| assembly, shipped-score top-75 | 2.515 | 3.318 | 3.032 | 3.175 |
| **assembly, ORACLE top-75** | **1.158** | **1.420** | **1.028** | **1.054** |

Paired assembly-oracle-75 vs pool-oracle-75: **-0.902 A [-1.090, -0.732], 24W/0L** on the
average, **-1.019 A [-1.234, -0.818], 24W/0L** projected. So the assembly bank's *set mean*
under a perfect filter is 1.312 A better than the pool's, not just its floor. On its face
this rescues asm's reading and makes assembly the largest generation-side headroom in the
sprint.

### H2. The capacity null destroys it

Choosing the best 75 of 2.3e5 is 460x more selection freedom than choosing the best 75 of
500, and an ORACLE selection converts raw cardinality into apparent accuracy (record C4;
`asm_FINDINGS.md` s2b already found this for the *floor*). I built the matched null:
**m i.i.d. chains drawn from the library's own Ramachandran marginal** — same count, same
`build_ca_exact` builder, same ideal geometry, **zero sequence information, zero fragment
structure, zero residue-residue correlation** — through the identical oracle-75 -> average
-> project path.

| arm (24 targets, mean) | oracle-best member | oracle-75 set mean | emits (avg) | emits (`fit`) | oracle-25 (avg) |
|---|---|---|---|---|---|
| K=500 pool | 1.824 | 2.732 | 1.930 | 2.073 | 1.604 |
| deployable assembly bank | 1.158 | 1.420 | 1.028 | 1.054 | 1.023 |
| **i.i.d. Ramachandran NULL, same size** | **1.188** | 1.740 | **0.900** | **0.838** | **0.844** |

| paired | d | CI95 | W/L |
|---|---|---|---|
| assembly vs **NULL**, oracle-75 avg | **+0.128** | [+0.023, +0.240] | 7/17 |
| assembly vs **NULL**, oracle-75 projected | **+0.216** | [+0.093, +0.355] | 6/18 |
| assembly vs **NULL**, oracle-25 avg | **+0.179** | [+0.078, +0.289] | 6/18 |
| NULL vs pool, oracle-75 avg | -1.029 | [-1.234, -0.835] | 23/1 |

**The null does not merely match the assembly bank — it beats it, on every arm, with CIs
excluding zero.** A quarter of a million random ideal-geometry chains, oracle-filtered to 75
and averaged, emit 0.838 A; the sequence-informed, fragment-derived assembly bank of the
same size emits 1.054 A. The entire "assembly improves the achievable set" effect is
**cardinality**, and the real library is worse at exploiting it than noise.

### H3. Verdict on the assembly leg of the central claim

* asm's *conclusion* ("no deployable version beats the incumbent", +0.218 A vs the 3.204
  synthesis) is correct and I did not shake it.
* asm's *mechanism reading* ("the candidate set improves and the objective cannot use it")
  is **not supported**: (i) the floor moving is worth 0.039 x 0.536 = 0.021 A to the
  operator (section D), and (ii) the set-level improvement that would have been worth
  something is a capacity artefact that pure noise reproduces and exceeds.
* **This leg of the coordinator's central claim must be withdrawn.** The assembly experiment
  establishes nothing about the objective's ranking ability, because the thing the objective
  failed to find was equally findable in random Ramachandran noise.
* Caveat on my own control: n = 24 and **none are FAIL18**, so this does not speak to
  assembly on the failure class. The 24 targets are the first by (fold, pdb) and are
  therefore an arbitrary, not adversarial, sample.

---

# PART II — the coordinator's re-prioritisation: the TERMINAL x OBJECTIVE-QUALITY surface

Scripts `s12/adv_terminal.py` (surface + re-pricing + the null for my own law) and
`s12/adv_mladder.py` (the deployable ladder through the real projection).
Results `s12/results/adv_terminal.json`, `s12/results/adv_mladder_s{0..3}.json`.

---

## I. THE SURFACE — optimal m as a function of objective quality (126 targets)

Coordinate-average arm (the obj agent validated it as a +0.16 A proxy at r = 0.995; the
headline deployable row is re-run through the real projection in section J). Quality is
swept two ways: a **rank-skill** ladder `a*z(rr) + (1-a)*z(noise)` (fixes Spearman skill
directly) and an **MAE** ladder (L1 against the native matrix plus i.i.d. noise bisected to
a requested MAE — the obj agent's well-behaved `iid` family). `rho` = mean Spearman of the
scorer against true CA-RMSD over the K = 500 pool.

| scorer | rho | m=1 | m=3 | m=5 | m=10 | m=25 | m=50 | **m=75** | m=150 | m=300 | m=500 | **best m** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ORACLE perfect ranking | +1.000 | 1.711 | 1.496 | **1.481** | 1.518 | 1.644 | 1.825 | 1.963 | 2.275 | 2.747 | 3.396 | **5** |
| ORACLE a=0.7 | +0.885 | 2.074 | 1.777 | 1.711 | **1.702** | 1.831 | 1.969 | 2.081 | 2.362 | 2.786 | 3.396 | **10** |
| ORACLE L1, MAE 0 (perfect matrix) | +0.840 | 1.994 | 1.837 | **1.832** | 1.875 | 1.968 | 2.099 | 2.218 | 2.495 | 2.967 | 3.396 | **5** |
| ORACLE L1, MAE 1.0 | +0.835 | 2.141 | 1.903 | **1.893** | 1.906 | 1.999 | 2.114 | 2.222 | 2.499 | 2.975 | 3.396 | **5** |
| ORACLE L1, MAE 2.339 (shipped MAE, iid errors) | +0.819 | 2.384 | 2.156 | 2.105 | **2.053** | 2.087 | 2.181 | 2.267 | 2.523 | 2.977 | 3.396 | **10** |
| ORACLE a=0.5 | +0.667 | 2.432 | 2.082 | **2.041** | 2.056 | 2.137 | 2.250 | 2.347 | 2.565 | 2.910 | 3.396 | **5** |
| ORACLE L1, MAE 6.0 | +0.719 | 3.187 | 2.849 | 2.737 | 2.625 | 2.527 | **2.475** | 2.502 | 2.658 | 3.007 | 3.396 | **50** |
| **SHIPPED distogram** | **+0.568** | 3.454 | 3.274 | 3.231 | 3.146 | 3.075 | 3.068 | **3.048** | 3.072 | 3.171 | 3.396 | **75** |
| score+consensus fusion | +0.525 | 3.478 | 3.336 | 3.302 | 3.297 | **3.273** | 3.274 | 3.273 | 3.330 | 3.397 | 3.396 | 25 |
| consensus centrality alone | +0.425 | 3.706 | 3.678 | 3.655 | 3.636 | 3.620 | 3.654 | 3.688 | 3.708 | 3.538 | **3.396** | 500 |
| ORACLE a=0.3 | +0.368 | 3.087 | 2.755 | 2.609 | **2.575** | 2.633 | 2.691 | 2.764 | 2.907 | 3.103 | 3.396 | 10 |
| ORACLE a=0.1 | +0.105 | 4.176 | 3.529 | 3.313 | 3.192 | **3.176** | 3.187 | 3.202 | 3.247 | 3.319 | 3.396 | 25 |
| RANDOM (3 seeds) | ~0.00 | 4.24-4.68 | 3.79-3.99 | 3.67-3.81 | 3.55-3.67 | 3.50-3.53 | 3.44-3.46 | 3.43 | 3.39-3.42 | 3.40 | **3.396** | 150-500 |

Three exact cross-checks land, independently derived: ORACLE-L1-MAE-0 at m = 75 =
**2.218** reproduces `obj_FINDINGS.md:29`'s oracle 2.218; perfect ranking at m = 25 =
**1.644** reproduces correction C1's 1.644; perfect ranking at m = 75 = **1.963**
reproduces `forensics_FINDINGS.md:302`. The surface is not a new instrument, it is the
existing numbers laid out on a second axis nobody had varied.

### I1. The crossover, and where the shipped system sits

**Optimal m is a steep decreasing function of the objective's rank skill.**

| objective rank skill rho | optimal m | gain of optimal m over m = 75 |
|---|---|---|
| ~1.00 | **5** | **-0.482 A** |
| 0.82-0.89 | 5-10 | -0.21 to -0.39 A |
| 0.67 | 5 | -0.306 A |
| **0.568 (SHIPPED)** | **75** | **0.000 A** |
| 0.37-0.53 | 10-25 | -0.19 to -0.00 A |
| ~0.10 | 25 | -0.026 A |
| ~0.00 (random) | 150-500 | -0.01 to -0.04 A |

**The shipped system is AT the crossover, not below it.** At rho = 0.568 the surface's
argmin over m is exactly the production m = 75, and the neighbouring cells are within
0.03 A (m = 50: +0.020, m = 25: +0.027, m = 150: +0.024). The coordinator's hypothesis that
the project sits in a local optimum that *forecloses* its own improvement path is
**half right and half wrong, and the half that is wrong matters**:

* **Wrong:** m = 75 is not a bad choice for the objective the system has. There is no free
  0.2 A sitting in the terminal today.
* **Right, and this is the important part:** the terminal is *co-optimal with a bad
  objective*, and every objective/ranking ceiling in this sprint was measured with m frozen
  at 75. **Those ceilings are understated by 0.2-0.5 A.**

### I2. The correction this forces on the sprint's headline ceilings

| ceiling as reported (m = 75) | as reported | at the CO-OPTIMISED m | correction |
|---|---|---|---|
| perfect CA-CA distance matrix (`obj` s0) | 2.218 | **1.832** (m = 5) | **-0.386 A** |
| perfect FILTER of the K=500 pool (`forensics` 7b) | 1.963 | **1.481** (m = 5) | **-0.482 A** |
| an objective at rho = 0.7 | 2.081 | **1.702** (m = 10) | **-0.379 A** |
| distogram at the shipped MAE but with i.i.d. errors | 2.267 | **2.053** (m = 10) | -0.214 A |

Translated to the projected scale (+0.16 A), a **perfect distance matrix with a
co-optimised terminal emits ~1.99 A**, not the 2.395 A the objective agent reports — i.e.
it *does* reach the sprint's 2.0 A target, where the fixed-m version misses it by 0.4 A.
**`obj_FINDINGS.md`'s conclusion "no objective upgrade of any kind reaches a 2.0 A mean"
must be restated as "no objective upgrade reaches 2.0 A *at m = 75*".** That is a different
and much less discouraging statement, and it is the single most consequential correction in
this audit.

### I3. Re-pricing the eleven negative ranking results — what the terminal did and did not do

Same ranking contrast, same pools, every terminal (126 targets, paired):

| contrast | m=1 (argmin) | m=5 | m=25 | **m=75** | compression m75 / m1 |
|---|---|---|---|---|---|
| ORACLE perfect ranking vs shipped | **-1.743** [-1.97,-1.52] 126/0 | -1.750 | -1.431 | **-1.085** [-1.27,-0.92] 125/1 | **0.62** |
| shipped vs random | **-0.786** [-1.14,-0.46] 78/44 | -0.434 | -0.427 | **-0.385** [-0.57,-0.21] 80/46 | **0.49** |
| score+consensus fusion vs shipped | **+0.024** [-0.15,+0.20] (null) | +0.071 | +0.198 | **+0.225** [+0.10,+0.35] (harmful) | — |
| consensus alone vs shipped | +0.252 | +0.424 | +0.545 | +0.639 | — |

**The honest answer to the coordinator's question: the terminal COMPRESSES ranking effects
by a factor of ~1.6-2.0, but it does NOT flip any sign in the helpful direction, and it does
not send a real effect to zero.**

Quantitatively, this is what the sprint needs in order to decide whether S9-8's closure is
safe. Using the measured compression (m75 = 0.49-0.62 x argmin):

* a ranker that measured **0.00 A** at m = 75 was worth **at most ~0.00-0.06 A** through
  argmin — still nothing;
* a ranker at the sprint's ignore threshold, **0.03 A** at m = 75, was worth **~0.05-0.06 A**
  through argmin — still below any threshold of interest;
* only a ranker showing **>= 0.10 A** at m = 75 would have been worth >= 0.20 A through
  argmin, and no discarded ranker in the record reached that.

**S9-8's closure ("in-band discrimination is informationally impossible") is therefore SAFE
against this particular objection.** The eleven negatives were measured through a lossy
instrument, not a blind one, and correcting for the loss does not resurrect any of them.
This is a survived audit of the coordinator's own worry, and I record it as such.

**But the direction the compression runs matters for FUTURE work, and this is where the
re-pricing does change practice.** The correct experiment for any new ranker is
(ranker x m) *jointly*, because improving the ranker moves the optimal m. Going from the
shipped objective to rho = 0.7 is worth -0.967 A at fixed m = 75 and **-1.346 A** at the
co-optimised m — the fixed-m evaluation understates a ranking improvement by **39 %**. A
future ranker must be scored on the surface, not at the point.

### I4. THE NULL FOR MY OWN LAW — and it partially breaks it

The coordinator asked for this and was right to. My section D law was fitted on ONE
perturbation design (oracle insertion into the top-75). Re-fitted on three designs that
vary the set a different way:

| design | slope on d_set_best | slope on d_set_mean | R^2 | R^2, set_mean alone | R^2, set_best alone | n |
|---|---|---|---|---|---|---|
| **oracle insertion (section D)** | **0.039** | 1.162 | 0.893 | 0.891 | 0.187 | 882 |
| random 75-subsets, m fixed | 0.184 | 0.771 | 0.665 | 0.648 | 0.300 | 756 |
| **22 different scorers, m = 75 fixed** | **0.300** | 0.790 | 0.839 | 0.799 | 0.467 | 2646 |
| cardinality axis (m varied, shipped scorer) | 0.264 | 0.675 | 0.702 | 0.548 | 0.209 | 1260 |

**What survives:** the set-mean dominance. Across all four designs the coefficient on
`d_set_mean` is 0.68-1.16 and `set_mean` alone explains 55-89 % of the variance. That is
the load-bearing half of the law and it is robust.

**What does NOT survive:** the coefficient **0.039 on `set_best` is design-specific.**
Out of design it is **0.18-0.30** — the operator is partially sensitive to the set's best
member, at roughly a quarter to a third of the weight it gives the mean. My section D
statement "the operator is almost blind to the set best" is **too strong as a general law**
and must be narrowed to: *the operator is almost blind to a set-best improvement that is
carried by one or two members out of 75 and leaves the bulk unmoved* — which is exactly the
oracle-insertion design, and exactly the assembly case, and is confirmed there by a DIRECT
measurement (j = 1 is worth -0.029 A) rather than by the fitted coefficient.

**Two further caveats I record against myself:**

1. The cardinality design has a positive intercept (+0.027) and the lowest set_mean R^2
   (0.548) — i.e. varying m introduces exactly the confound the record's correction **C3**
   warns about ("cardinality confounds it"). I reproduce C3's confound and it applies to my
   own section I table too: the m-ladder is not a clean composition experiment.
2. The **"fixed quantile" statement in D4 is established only within a fixed cardinality.**
   Across cardinalities it is trivially false (at m = 1 the emitted structure *is* the set),
   and the m-design's gain on set_mean (0.675) differs from the j-design's (1.162). Quote
   the fixed-quantile result as "at m = 75, the operator emits at quantile 0.20-0.25 of its
   input set almost regardless of composition" and not as a general property of averaging.

---

## J. THE DEPLOYABLE FALLBACK — a fixed small-m terminal, through the REAL projection
## **CLEAN NEGATIVE. Leave-fold-out selection of m picks m = 75 on all five folds. There is
## no deployable change here. The shipped terminal is correct for the shipped objective.**

Script `s12/adv_mladder.py` (4 shards) + `s12/adv_bigm.py`. All 126 targets, shipped
objective, production terminal, lambda = 0 `fit` arm per the coordinator's constraint (1).

| m | all 126 | FAIL18 | other-108 | d vs m=75 | CI95 | W/L | drop-10 |
|---|---|---|---|---|---|---|---|
| 1 (argmin) | 3.4516 | 6.007 | 3.026 | +0.2464 | [+0.120, +0.383] | 49/77 | +0.359 |
| 3 | 3.3700 | 5.981 | 2.935 | +0.1648 | [+0.058, +0.279] | 49/77 | +0.255 |
| 5 | 3.3540 | 5.981 | 2.916 | +0.1488 | [+0.048, +0.251] | 44/82 | +0.238 |
| 10 | 3.2866 | 5.954 | 2.842 | +0.0814 | [+0.000, +0.164] | 56/70 | +0.155 |
| 25 | 3.2347 | 5.980 | 2.777 | +0.0295 | [-0.031, +0.086] | 59/67 | +0.093 |
| 50 | 3.2301 | 6.007 | 2.767 | +0.0249 | [-0.010, +0.058] | 58/68 | +0.058 |
| **75 (production)** | **3.2052** | **6.034** | **2.734** | — | — | — | — |
| 150 | 3.2410 | 5.918 | 2.795 | +0.0359 | [-0.028, +0.100] | 58/68 | +0.098 |

Instrument check: m = 75 reproduces the pinned synthesis at **3.2052 vs 3.2041** (the
0.0011 A gap is the A3 `fit`-arm reproduction tolerance).

* **Every alternative m is worse, and m = 75 is the minimum.** The curve is flat-bottomed
  (m = 25/50/150 are within 0.036 A) but the sign is consistent and no CI favours a change.
* **Leave-fold-out choice of m selects 75 on fold 0, 1, 2, 3 and 4. d = 0.0000 exactly.**
  The honest deployable fallback the coordinator asked for **does not exist**.
* ORACLE per-target choice of m — the ceiling of *any* router over this axis —
  **2.8919, d = -0.3133 [-0.385, -0.250], 114W/0L.** Real, but oracle, and the record
  already contains two independent null routers (LFO AUC 0.558 vs permutation 0.511, and
  the fail agent's second null). Nothing here changes that.

### J1. The one place where the terminal choice IS large — and it is not deployable either

The FAIL18/other-108 interaction is the largest effect I found in this whole part.
`s12/adv_bigm.py` runs the large-m cells through the real projection on the FAIL18 and on a
**length- and fold-matched control set of 18 non-FAIL18 targets** (MATCH18), so the contrast
cannot be a group-composition artefact:

| m | FAIL18 | d vs m=75 | CI95 | W/L | MATCH18 | d vs m=75 | CI95 | W/L |
|---|---|---|---|---|---|---|---|---|
| 75 | 6.034 | — | — | — | 2.683 | — | — | — |
| 150 | 5.918 | -0.115 | [-0.294, +0.052] | 12/6 | 2.836 | +0.153 | [+0.023, +0.306] | 6/12 |
| 300 | 5.556 | -0.477 | [-0.945, -0.076] | 15/3 | 3.039 | +0.356 | [+0.169, +0.560] | 4/14 |
| **500 (= NO distogram filter at all)** | **5.524** | **-0.510** | **[-0.857, -0.177]** | **15/3** | **3.314** | **+0.631** | **[+0.319, +0.983]** | **4/14** |

`m = 500` is literally the whole K = 500 BLOSUM pool averaged — the **no-distogram control**
(verified: the shipped and random scorers give bit-identical m = 500 emissions, max
difference 2e-14). So:

**On the 18 failures the distogram filter is worth -0.51 A: the system would be better off
throwing the distogram away entirely. On matched controls the same filter is worth +0.63 A.
That is a 1.14 A interaction, both halves with CIs excluding zero.** This is the sharpest
possible statement of forensics' "the filter is anti-ranked on the 18", now measured as an
emitted-structure cost through the production terminal, with a matched control.

And it prices the router idea exactly: a **perfect** binary router (m = 500 on the 18,
m = 75 on the 108) emits (18x5.524 + 108x2.734)/126 = **3.133**, i.e. **-0.073 A** — an
oracle upper bound on a two-state router, an order of magnitude smaller than the -0.313 A
available from a per-target oracle m. The router is not where the accuracy is.

---

## K. WHAT PART II CHANGES, AND WHAT IT DOES NOT

**Does not change:** the shipped terminal. m = 75 is LFO-optimal, and every deployable
alternative I tested (m in {1,3,5,10,25,50,150}, LFO-chosen m, and the 30 operator forms the
aggregation agent had already tested at fixed m) loses. `agg_FINDINGS.md`'s "the production
operator is already at the native-free optimum" survives, and I have now extended it along
the cardinality axis it did not cover with a proper LFO selection.

**Does not change:** S9-8's closure. The terminal compresses ranking effects by ~2x, not to
zero; correcting for the compression turns a discarded 0.03 A result into 0.05-0.06 A and
resurrects nothing. The eleven negatives stand.

**Changes, and this is what the coordinator should act on:**

1. **Every objective/ranking ceiling in this sprint is understated by 0.2-0.5 A**, because
   all of them were measured at m = 75, which is co-optimal with the *bad* objective they
   were replacing. A perfect distance matrix emits 2.218 at m = 75 and **1.832 at m = 5**
   (~1.99 A projected) — i.e. it *does* reach the sprint's 2.0 A goal. `obj_FINDINGS.md`'s
   flat "no objective upgrade of any kind reaches a 2.0 A mean" is true only at fixed m = 75
   and must be restated.
2. **Any future objective or ranking work must be scored on the (ranker x m) surface, not at
   m = 75.** Fixed-m evaluation understates a ranking improvement by ~39 % at rho = 0.7 and
   by ~25-33 % at the qualities actually reachable.
3. **The FAIL18 need a different terminal, by 1.14 A of interaction, and nobody can tell
   which targets they are.** That is the same wall forensics hit (early-warning AUC 0.600
   inside a 0.616 null), now with a much larger prize attached to breaking it than anyone
   had priced: identifying the failures is worth -0.073 A through a binary router but the
   *mechanism* it reveals — that the filter is actively harmful on exactly the class whose
   deposited conformation the sequence does not determine — is the sprint's cleanest causal
   result.

**Nothing here warrants a dev24 pass.** The only deployable candidate (fixed small m) is a
clean negative at d = 0.0000, and there is nothing to confirm.

---

## L. CLAIM 5 — everything reported as a POSITIVE, and the three process hazards

### L1. Baseline audit (3.454 vs 3.204) — **no agent commits the error**

I checked every findings file for the 0.25 A baseline trap. `asm_FINDINGS.md:286-292`
quotes both and takes the **3.204 synthesis** as the headline comparator (its
`asm_avg75` uses `fit_ca`, so it is like-for-like); `agg` baselines on `avg75` = 3.048 raw
and reports argmin only as a ladder rung; `obj` baselines on 3.205; `key` uses the `fit`
arm; `lit_FINDINGS.md:25` states the caveat explicitly. **Clean across the board.**

### L2. Tie-breaking — **no leak, but a latent hazard worth pinning**

The shipped score has **35.2 exactly-tied values per 500-member pool** (measured, A1) — real
duplicate windows retrieved through two parents. The record's 1.386 A phantom arose when
`np.argmin` on a tied signal read an ORACLE-sorted array. In this sprint the pool arrays are
in **BLOSUM order** (`u["order"][:500]`), not oracle order, so every tie-break reads
BLOSUM rank: deployable and benign. I grepped all agent code for `argsort(` without
`kind=`: 30 sites, all either on a continuous signal with no ties, or on an oracle array
inside a correctly-labelled ORACLE arm. **No instance of the leak.** Two residual notes:
(i) `s12/asm_deploy.py:102` and `s12/asm_overlap.py:51` use the default (unstable) sort on a
tied score, which makes those results non-reproducible across numpy versions though not
biased; (ii) `s12/agg_ladder.py:142,156` superpose the oracle subsets on the **oracle-best
member** rather than the subset medoid, which makes `agg`'s oracle arms ~0.036 A optimistic
against the medoid operator (`ORC_top25avg_pool` 1.608 vs C1's and forensics' 1.644). Quote
1.644, not 1.608.

### L3. Concentration — one reported positive is much more concentrated than it reads

| claim | effect | drop-top-10 | drop-top-20 | verdict |
|---|---|---|---|---|
| `obj`: ORACLE distance matrix vs shipped | -0.811 | -0.597 | — | robust, 116W/10L |
| `key`: ORACLE abego4 key vs BLOSUM | -0.293 | -0.129 | **-0.040** | **86 % of the effect lives in 20 of 126 targets** |
| `key`: ORACLE ss3 key vs BLOSUM | -0.258 | -0.092 | **-0.014** | even more concentrated |
| `fail` E5: ESM contact rescore | -0.126 | -0.024 | +0.016 | agent already kills it with its own null |
| `agg`: best deployable operator | -0.0038 | +0.021 | — | agent already reports it as null |
| `asm`: `asm_avg75` vs argmin | -0.032 | +0.066 | — | agent already reports it as null |

The `key` agent reports its drop-20 numbers honestly in its own table, so this is a reading
hazard rather than a reporting failure — but **the oracle retrieval-key gain of -0.29 A
should not be summarised as a broad effect**: it survives drop-10 and essentially vanishes
at drop-20, and the `key` agent's own capacity null shows half of the FAIL18 portion is
cardinality (random bin key worth -0.353 A on FAIL18, 12W/6L).

### L4. Cache / partial-result hazard — **one live instance**

`I.write(name, obj)` overwrites by name with **no config key**, so a smoke-test or
`--limit` run silently replaces a full run's results file. I hit this myself
(`adv_decode.json` held an n = 3 smoke test while the n = 126 run was still going, and it is
readable as if complete). One live instance in the sprint's own results:

* **`s12/results/coord_null.json` contains 56 of 126 targets** while
  `coord_FINDINGS.md:158` says the C4 sequence-information ablation is "queued — not yet
  run". Anyone reading the file rather than the findings will quote a partial 2x2. It should
  be renamed or deleted.

Everything else checks out: I swept all 90 non-`adv_` results files for target counts and
every headline artefact carries 126 (or a declared subset — `asm_chain2_M40_k3` 40,
`coord_corpus` 81, `assembly_e1_floors.m4_equal` 93, `vq_stage1_A3` 93). **False alarm
resolved:** `vq_stage3_B.json`'s `validation.n = 24` is a 24-row dumped-vs-recomputed
consistency check on tuning targets (1A13, 1A1P, …), **not** a dev24 run. No agent touched
dev24 or benchmark60 as far as I can determine from the artefacts.

### L5. Hyperparameters chosen in-sample

I found no case of an in-sample-tuned hyperparameter reported as pre-registered. The two
places where it could have happened are handled correctly: the `key` agent **pre-registers**
its prediction in writing before computing any predicted-key number
(`key_FINDINGS.md:211`), and `obj`'s learned arms are LFO-fitted with the assertion in
`fit_all`. My own m-ladder (section J) is the one place a hyperparameter was at stake and I
selected it leave-fold-out; it chose the incumbent.

---

## M. RANKED VERDICT

### REFUTED

1. **"The objective's own geometric optimum is 4.62 A"** (`obj_FINDINGS.md`, §1 MDS column).
   That is a weak-decoder artefact. The honest argmin of the shipped Bayes-risk score inside
   the emittable structure class reaches **3.532 A** (-1.089 [-1.256, -0.923], 108W/18L).
   The *conclusion* the number supported survives at a quarter of the claimed margin
   (+0.328 [+0.213, +0.448]) and my attack strengthens it on the score axis.
2. **"The assembly experiment shows the objective cannot rank within an improved set"**
   (`asm_FINDINGS.md` §5, and the coordinator's leg 2). Two independent measurements kill
   the inference: a perfect rank-1 on the pool is worth -0.029 A through the deployed
   operator, and the assembly bank's set-level advantage is a **capacity artefact that
   2.3e5 random Ramachandran chains reproduce and beat** (+0.216 A [+0.093, +0.355] in the
   real bank's disfavour, 6W/18L). asm's *conclusion* (no deployable gain) stands; its
   *mechanism* does not.
3. **"No objective upgrade of any kind reaches a 2.0 A mean"** (`obj_FINDINGS.md` §7b).
   True only at the frozen m = 75. At the co-optimised terminal a perfect distance matrix
   emits **1.832 A** (avg) / ~1.99 A projected — it reaches the goal.

### WEAKENED (with the quantitative caveat)

4. **r_sep = 0.196.** Reproduces exactly, and the small-shell degeneracy is worth only
   +0.007/+0.013 — but a pooled, leave-fold-out partialling gives **0.367-0.391**, a factor
   of two. The contrast against a matched-MAE null survives every estimator (distogram keeps
   41-60 % of the null's within-shell signal), so the argument holds; the number must not be
   quoted as "almost no pair-specific information".
5. **The FAIL18 subgroup contrasts.** Nested and stable across band/K/M (never a swap), and
   the recall distribution is bimodal (18 at zero, only 4 at 1-3), so the threshold is not
   arbitrary. **But** under an *absolute* band definition the set shares only 10 of 18
   members, the emitted gap shrinks 46 %, and **the fibril/lasso enrichment stops being
   significant (p = 0.097 at 2.5 A, 0.244 at 2.0 A)**. And the headline "the 18 return
   ~6.0 A" plus 0.15 of the 0.39 `band_score_pct_min` gap are entailed by the definition.
6. **My own section-D law.** The set-mean dominance is robust across four perturbation
   designs (slope 0.68-1.16, R^2 0.55-0.89), but the **0.039 coefficient on set-best is
   design-specific — out of design it is 0.18-0.30.** Narrow the claim to "blind to a
   set-best gain carried by one or two members of 75", which is what the direct j = 1
   measurement shows anyway.

### SURVIVED (with what I tried)

7. **The fibril/lasso enrichment.** Bonferroni over the full 66-hypothesis post-hoc family
   (8.1e-05), max-statistic permutation over the same family (0/4000), leave-one-out over
   all 126 (worst p 5.4e-06), and an independent re-derivation of the labels from the raw
   headers (**0 disagreements on fibril, lasso and cyclic**). Did not move.
8. **`coord_benchN`: no adequately-powered fresh benchmark exists.** I re-parsed all 93
   rejections independently. All 67 length rejections are correct — **not one of those files
   contains any chain of 9-16 residues** (they are 1-8mers). All 16 CA-step rejections have
   genuinely missing interior residues. The only defensible relaxation (longest contiguous
   run) rescues 3 structures, taking 40 -> 43 and the final count 16 -> ~18. Not 40.
9. **`agg`: the production operator is at its native-free optimum.** Confirmed and extended
   along the cardinality axis with a proper LFO selection: **m = 75 is chosen on all five
   folds, d = 0.0000.** No deployable terminal change exists.
10. **S9-8's "in-band discrimination is informationally impossible"** — the coordinator's own
    worry that eleven ranking negatives were mis-instrumented. The terminal compresses
    ranking effects by ~2x (m75/m1 = 0.49-0.62), not to zero; a discarded 0.03 A result was
    worth 0.05-0.06 A through argmin. **The closure is safe.**
11. **The instrument** (`s12/instrument.py`). 126/126 exact filter reproduction under three
    precisions, `coordinate_average` to 2e-7 A, `fit` projection to 7e-4 A, distogram cache
    verified at full strength. One small real defect: the lambda = 0.3 arm carries a ~0.02 A
    reproducibility floor (A3). Nothing contaminated.
12. **`coord_contain`: the 4 verbatim leaks.** Verified; the coordinator's own conservative
    reading (the leak runs *against* the failure story) is correct.

### DEFECTS FOUND

* `s12/instrument.py:125-133` — `project`'s lambda arm is not reproducible against
  production to better than ~0.02 A (`core/pipeline.py:145-152` explains why). Use `fit_ca`.
* `s12/results/coord_null.json` — 56/126 targets on disk for an experiment the findings
  declare unrun. Rename or delete before anyone quotes it.
* `s12/agg_ladder.py:142,156` — oracle subsets superposed on the oracle-best member, not the
  medoid; makes `agg`'s oracle arms ~0.036 A optimistic (1.608 vs the correct 1.644).
* `s12/asm_deploy.py:102`, `s12/asm_overlap.py:51` — unstable `argsort` on a tied score;
  not biased, but not reproducible.
* `s12/instrument.py:212` — `I.write` has no config key, so partial runs silently overwrite
  full ones. This is how `coord_null.json` got into its current state and it will happen
  again.
