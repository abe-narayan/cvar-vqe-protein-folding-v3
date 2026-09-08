# SPRINT 18 — EXPERIMENT WORKSTREAM — PRE-REGISTRATION

Written **before any arm was run**. n = 126 tuning targets, pinned order, `s12.instrument.targets()`.
Start structure for every arm is **exactly the same** coordinate average of the shipped top-75
(Control A), so no arm can win by starting somewhere better. Native read only to score.
MDE at 80% power on this instrument = **0.084 Å**; a null smaller than that is uninformative.

Metric: full-chain Cα-RMSD, frozen implementation, all residues, proper rotations, model 1.
Seeding: `s15.seed.stable_rng` throughout. Sealed 60-target benchmark: not read, not probed.

---

## AMENDMENT — the object is MATH's, and it is RESIDUE-additive

MATH landed `s18/math_anova.py` after this pre-registration was written and before any arm was
read. Their deliverable is consumed as-is; `s18/exp_obj.py` adds only analytic gradients (verified
against central differences of *their* callables) and the convex λ mix. Three consequences,
recorded before results:

1. **MATH's `E_le1` is the RESIDUE-additive object** — the ANOVA variable is `θ_r = (φ_r, ψ_r)`
   on the 2-torus. So arms B2 and B3 of the table above are the same arm, and my provisional
   coordinate-additive object is retired. MATH's separate `harmonic_truncate(order=1)` is the
   "degree-1 in angles" object and is run as its own arm (`E_h1`).
2. **`f_0 = f_{n−1} = 0` EXACTLY** (MATH's locality corollary): the degree-1 object exerts *no
   force at all* on the two terminal residues, which the frozen metric scores. MATH's shipped
   `argmin_le1` resolves the tie by taking `argmin` of an all-zero table, which returns mesh cell
   (0,0) — a fixed torus corner. That is an arbitrary structural commitment on 2 of n residues.
   **I run both**: MATH's as shipped, and `argmin_le1_hold`, which holds undetermined residues at
   the start value. The difference prices the terminal patch, and I pre-commit to reporting it.
3. **MATH revised their own declaration mid-sprint** (G 24 → 16, S 2048 → 512, a moment-based
   estimator, no self-recentring of the mesh, and a new `E_le1_ang`). I re-synced before running
   anything, read `GRID`/`NSAMP` dynamically from their module rather than pinning my own, and
   persist their `CONFIG` into every artefact. Two arms follow from the new API and are added
   here **before results**:
   - `refine_ang` / `argmin_ang` — MATH's **sub-residue (angle-additive)** object
     `E_0 + Σ a_r(φ_r) + Σ b_r(ψ_r)`. Reported as *the closest continuous relative* of the
     lattice's per-qubit truncation, never as the same object.
   - **F5 has a specific resolution from MATH, and it is not "the mapping is invalid":** the two
     lattice bits index k-means clusters of the joint (φ,ψ) library and do not correspond to φ
     and ψ at all, so **no continuous object equals the strict Walsh weight-≤1 projection.** My
     job is therefore to measure the residue-additive and angle-additive objects and report the
     gap between them, not to claim either *is* the lattice object. MATH owns the verdict on F5.
4. **Quadrature.** MATH originally declared S = 2048, G = 24. That is ~190 s/target/fit; with the controls
   and the μ swap it does not fit the compute envelope. I will **measure** the S-sensitivity of
   the *downstream outcome* (S ∈ {256, 512, 2048} on a stratified subset) and run the 126-target
   experiment at the smallest S whose outcome matches S = 2048 by **less than the 0.084 Å MDE**.
   If no reduced S matches, I run at 2048 and drop an arm instead. The chosen S is reported with
   its evidence, never assumed.

## P0 — VALIDITY GATE (run first, no conclusions drawn from it)

**P0.a** λ = 1 of the ladder must reproduce s17 `refine_full` to < 0.01 Å mean.
**P0.b** `E_le1` + `E_ge2` must equal `E_full` pointwise to machine precision by construction.
**P0.c** Monte-Carlo noise in the ANOVA tables must be small against the table's own range:
report `mc_se / range(f_i)` per target. Gate: median ratio < 0.05.
**P0.d** The trigonometric interpolant must reproduce its own grid values to < 1e-8.
*If P0.a fails the harness is wrong and nothing downstream is reportable.*

---

## ARM LIST (all from the identical start)

| id | arm | native-free? |
|---|---|---|
| A0 | `avg` — coordinate average (**Control A**, the incumbent mechanism, ~3.048) | yes |
| A1 | `proj` — its ideal-geometry projection (the torsion start point) | yes |
| B1 | `refine_full` — L-BFGS on `E_full` (λ=1) | yes |
| B2 | `refine_le1` — L-BFGS on `E_le1` (λ=0) | yes |
| B3 | `refine_res` — L-BFGS on the residue-additive `E_res` | yes |
| B4 | `argmin_le1` — **exact separable global argmin** of `E_le1`, no search | yes |
| B5 | `argmin_res` — exact separable global argmin of `E_res` | yes |
| C1 | `rand_obj_le1` — B2 against a **shuffled distogram** (zero-information) | control |
| C2 | `rand_move_le1` — random torsion move of **matched displacement** to B2 | control |
| C3 | `const_obj` — refinement toward a **zero-information constant** objective | control |
| C4 | `rand_obj_full`, `rand_move_full` — s17's two controls, re-run | control |
| D1–D5 | **λ ladder** `E_λ = E_le1 + λ E_ge2`, λ ∈ {0, 0.25, 0.5, 0.75, 1} | yes |
| E1 | `greedy_full` — greedy 1-opt at matched evaluation budget | yes |
| E2 | `greedy_le1` — greedy 1-opt on `E_le1`, same budget | yes |
| E3 | `rls` — random local search, same budget | yes |
| E4 | `anneal` — Metropolis, same budget | yes |
| F1 | μ-sensitivity: `E_le1` under **uniform** μ instead of the pool marginal | yes |

**λ is pre-registered as exactly {0, 0.25, 0.5, 0.75, 1} and will not be extended.**

---

## HYPOTHESES, EXPECTATIONS, CONTROLS, CRITERIA, FALSIFIERS

### H1 — the sprint's question. Degree-1 refinement beats full refinement.
- **Expected (honest prior):** Phase 0 says the 19-target effect had median 0 and a CI spanning
  zero, and that degree-1 is a *worse* global correlate of RMSD (ρ +0.153 vs +0.264). I expect
  **B2 ≈ B1, both worse than A0.** I expect F1/F2/F4 to fire.
- **Strongest control:** C1 (shuffled distogram through the *identical* ANOVA machinery). If
  degree-1 gains as much with a shuffled distogram, the gain is the smoothing an additive
  surrogate applies, not the objective's information.
- **Success:** B2 (or B4) beats A0 at n=126 with a paired fold-clustered CI excluding zero, by
  more than 0.084 Å, **and** beats C1 and C2 at matched displacement.
- **Falsifier F1:** B2/B4 land ≥ 3.5 Å. **F2:** B2 − B1 CI includes zero. **F3:** the gain does
  not exceed C1/C2. **F4:** the effect is absent at n=126 though present on the 19.

### H2 — the λ decomposition (Phase 4, Control E). *The most important single experiment.*
- **Hypothesis:** if the higher-order component `E_ge2` is systematically **mis-specified**, RMSD
  is monotone increasing in λ, best at λ=0.
- **Three pre-declared readings, committed now:**
  - **monotone increasing in λ** → higher-order terms are mis-specified; degree-1 is the right
    truncation; this is the sprint's positive result.
  - **interior minimum** → neither pure object is right; the truncation level is a tunable and
    the "degree-1" framing is wrong. *I will report the interior λ but will NOT tune on it.*
  - **flat / monotone decreasing** → the truncation buys nothing; **F2 fires**.
  - **all five λ worse than A0** → the whole objective family is mis-aimed regardless of
    truncation; the branch closes and the finding is about the objective, not its degree.
- **Reported as:** full curve, per-target ΔRMSD distributions, medians, W/L, fold-clustered CIs
  at every λ — never endpoints alone.

### H2b — THE COORDINATOR'S HYPOTHESIS (added before results were read; `s18/COORD_FINDING.md`, `s18/results/objceil.json`, ledger L2)

The coordinator held the functional form fixed and varied only the distances: `d_α = (1−α)d̂ + α·d_true`,
from the identical coordinate-average start. Two results change what my λ curve is read against:

- **The functional form is sound.** At α = 1 the identical objective reaches **1.152 Å**, 80.2% of
  targets under 2.5. So degree-1 cannot be framed as "fixing a broken functional form" — I withdraw
  any such framing from H1's motivation. What is broken is the *distances*, not the *form*.
- **The distogram's error is worse than random error of the same size.** Destroying only the
  *direction* of the residuals while keeping their magnitudes gives **2.609** against α=0's
  **3.610** (−1.001 [−1.226, −0.784]); matched isotropic noise gives **2.573**. The error is
  structured and the structure is actively harmful.

> **PRE-REGISTERED HYPOTHESIS H2b.** If the *pairwise* component of the distogram's error is the
> mis-directed part, then a degree-1 objective — which averages over pairs — is **more robust to
> structured distogram bias**, and should close part of the 1.0 Å gap between α = 0 (3.610) and
> the shuffled control (2.609).

**Pre-declared readings, committed now:**
- degree-1 lands ≈ **3.6** → it has done nothing; the averaging does not launder the bias. **F1/F2.**
- degree-1 lands ≈ **2.6** → it recovers exactly what randomising the error direction recovers.
  A strong, interpretable positive: the truncation is a bias-robustness operator.
- degree-1 lands **below 2.6** → it does something the shuffle does not; the strongest outcome.
- degree-1 lands **between 2.6 and 3.6** → partial; report the fraction of the gap closed with its CI.

**Reference arms carried on the same targets and the same starts**, merged from `objceil.json`
rather than rebuilt: `shuffled` (2.609) and `isotropic` (2.573).

**Caveat I commit to stating whenever I cite the α curve** (the coordinator's own, and it cuts
against their result): α is a *blend toward truth*, not a model of how a better predictor would
err. A real improved distogram would carry its own error structure, not interpolate toward the
native. The α ladder is a requirement in the "how much residual must go away" sense, **not** a
prediction about any achievable predictor. Every α > 0 arm is **ORACLE**, a ceiling, never a result.

### H7 — WHAT LIMITS DEGREE-1: the distances, or the truncation? (added before results were read)

MATH's `Target.with_dhat` rebuilds `E_0`, `f_r`, `a_r`, `b_r` **exactly** for any distance vector
from the already-computed conditional moment tables, with **no new chain builds**. Three arms
therefore cost nothing and are added now, before any 126-target number was read:

| arm | distances fed to the degree-1 object | label |
|---|---|---|
| `shufobj_*` | the deployed `d̂` with its **pairs permuted** (weights `w` left in place, so only information is destroyed, not the weighting) | native-free CONTROL |
| `orc_le1`, `orc_argmin` | the **native's own** Cα–Cα distances | **ORACLE** |
| `orcshuf_le1` | the coordinator's construction: residual **magnitudes kept, direction destroyed** | **ORACLE** |

**Hypothesis.** The coordinator showed the *full* objective with perfect distances reaches 1.152 Å.
If the degree-1 object is limited by the **distances**, `orc_le1` should improve by a comparable
amount. If it is limited by the **truncation**, `orc_le1` stays poor even with perfect distances.

**Pre-declared readings, committed now:**
- `orc_le1` ≈ 1.2 Å → degree-1's failure is entirely the distogram's; the truncation is innocent
  and the branch's problem is the same problem the whole programme has.
- `orc_le1` stays ≳ 3 Å → **the truncation itself destroys the information**, and no distogram
  improvement could rescue a degree-1 objective. That closes the branch for a reason the λ curve
  alone cannot establish.
- intermediate → report the fraction of the full objective's oracle gain that degree-1 retains.

`orc_*` arms are ORACLE ceilings and are never reported as results, never mixed with the realized
numbers, and never used to select anything.

### H8 — GEOMETRIC-LEVERAGE RE-WEIGHTING (new arm, pre-registered before it was written)

> **CORRECTION, entered after the arm was pre-registered and before its results were read.**
> The coordinator has **withdrawn** the motivating ceiling. ADVERSARIAL audited `s18/objceil.py`
> and found line 163 passes `sd[pi]` — a *permuted weight vector* — into the fit, so the
> `shuf_paired` arm is not the deployed functional at all. Against `shuffled` (2.609, residual
> permuted, weights kept in place) it differs **only** in whether `sd` moves, so the extra
> 0.537 Å is attributable to permuting the **weights**, not to where the errors land. The
> sentence "the gap from 3.610 to 2.072 is entirely error *assignment*" **is withdrawn**, and
> the honest split is now unmeasured.
>
> **The arm stays pre-registered and is run as written**, at the coordinator's explicit request,
> so the record shows an arm that **outlived its original motivation**. It is a legitimate
> native-free hypothesis on its own terms. It is **not** sized or tuned against 2.072, and
> 2.072 is not quoted as its ceiling anywhere in the findings.
>
> The arms that still stand as reference are `shuffled` 2.609, `shuf_strat` 2.560, `isotropic`
> 2.573 — all with correct weights.
>
> *Not mine to run:* ADVERSARIAL owns `wperm_only` (real d̂, permuted weights) and `wflat`
> (real d̂, uniform weights). **My λ ladder has no weighting axis** — λ mixes two objectives at
> fixed `w` — so there is no duplication. H8's own control permutes the *leverage* factor while
> leaving `w = 1/sd²` in place, which is deliberately a different question from `wperm_only`.

The original motivation, recorded as it stood (`s18/COORD_FINDING.md`, ledger L3, n = 126):

    shuffled     plain permutation of residuals            2.609  −1.001 [−1.226,−0.784] vs α=0
    shuf_strat   permuted WITHIN separation bins           2.560  −1.050 [−1.257,−0.860]
    shuf_paired  residual AND its 1/sd² weight moved       2.072   ** RETRACTED, confounded **
    isotropic    matched-RMS Gaussian                      2.573  −1.037 [−1.273,−0.802]

**Arm.** A native-free re-weighting of the refinement objective by **geometric leverage** — how
much a given pair's *target* distance can move the emitted structure. At the Control-A start,
with `J` the superposed Cα Jacobian, `∇d_p` the pair-distance gradient and `H = 2 Gᵀ diag(w) G`
the Gauss–Newton Hessian of the objective, implicit differentiation of the stationarity condition
gives, exactly,

    dθ*/dd̂_p = 2 w_p H⁺ ∇d_p        lev_p = ‖J · dθ*/dd̂_p‖ / √n   (Å of structure per Å of d̂_p)

and the arm minimises the same objective with `w'_p = w_p · (lev_p / median_p lev)^(−2β)`.

**Pre-registered ladder, NOT to be extended:** β ∈ {0, 0.25, 0.5, 0.75, 1}. β = 0 is exactly the
deployed weighting and must reproduce `refine_full`; that is this arm's validity gate.

**Controls (both mandatory, both run).**
- *zero-information*: `lev_p` **permuted across pairs**, so the weight distribution is identical
  and only the leverage→pair assignment is destroyed. This is the exact analogue of the
  coordinator's own shuffle and is the control that decides whether leverage carries information.
- *matched-random*: a random re-weighting drawn to match the realised `w'` distribution.

**Success criterion.** Some β beats `refine_full` at n = 126 by more than 0.084 Å with a
fold-clustered CI excluding zero, **and** beats the permuted-leverage control by more than the MDE.

**Falsifier.** No β beats the permuted control → leverage carries no usable assignment
information and the branch closes with the coordinator's ceiling unclaimed.

**Two cautions I adopt verbatim from the coordinator and will repeat in the findings.**
1. *This is adjacent to a refuted experiment but is not it.* Sprint 17 refuted 58 functionals of
   the distogram **including uncertainty weighting** — for **SELECTION**. This is **REFINEMENT**:
   a different use and a different functional. It inherits neither that refutation nor the α
   ceiling, and is reported on its own evidence.
2. *2.072 Å is ORACLE and is not a target for this arm.* `shuf_paired` knows where the errors are
   because it constructs them; leverage does not. The honest prediction is "some fraction of the
   1.54 Å", and **a null here is a clean result, not a surprise.**

**Honesty condition on the ladder, stated now.** There is **no native-free rule that selects β**.
So if an interior β is best, that is a DIAGNOSTIC, not a deployable result, and it will be
labelled as one — the recorded failure mode of picking a restraint on dev-set RMSD.

### H3 — separability makes degree-1 search-trivial.
- **Hypothesis (near-certain, stated so it is not claimed as a discovery):** `E_le1` is a sum of
  univariate functions, so its global optimum is EXACT at O(n·G) table lookups. **EXACT**, a
  theorem, not a result.
- **Consequence to test:** s17 found greedy 1-opt certifies the *full* objective's optimum at
  1,024 evaluations. **Prediction:** budget is not binding for either objective. If neither is
  budget-bound, **search quality is not the sprint's problem and Control B is a null by
  construction** — I will say so rather than dress it up.

### H4 — the four quantities dissociate (BRIEF §5).
For every objective I report separately: (i) objective reduction %, (ii) RMSD change,
(iii) alignment on **four named axes** — pairwise ordering accuracy, Spearman ρ, Pearson r, and
in-band (≤ 1.5 Å band) Spearman — computed over a per-target sampled ensemble, and (iv) search
quality: evaluations to convergence and whether the certified/exact optimum is reached.
- **Expected:** a large objective reduction with a positive RMSD change, repeating s17. If
  degree-1 reduces less objective but more RMSD, that is the alignment axis moving, and it must
  be shown on the alignment measurement, not inferred from the argmin.

### H5 — Phase 9, protecting the coordinate average.
For every objective, three routes from the identical pool: **pool → coordinate average** (A0),
**pool → objective refinement**, **pool → selection by that objective (argmin over the 75)**.
- **Hypothesis:** the coordinate average survives all three for both objectives.
- **Architectural question answered by the numbers, declared now:** the new objective belongs
  *before averaging* (rerank/filter the 75), *after averaging* (refine), *on individuals*,
  *as a generator*, or **nowhere**. I commit to reporting "nowhere" if that is what the data say.

### H6 — μ-sensitivity (BRIEF §4 point 2).
- The reference measure is a modelling decision. **Success criterion for reporting the result at
  all:** the sign of the B2−A0 and B2−B1 comparisons must be the same under the pool marginal and
  under uniform μ. If the sign flips with μ, the degree-1 object is **not identified** and any
  conclusion about it is an artefact of the measure — that is a fourth way F5 can fire, and I
  will hand it to MATH.

---

## WHAT WOULD MAKE ME REPORT A POSITIVE
All of: (a) B2 or B4 beats A0 by > 0.084 Å with a fold-clustered CI excluding zero; (b) it beats
C1 **and** C2; (c) the λ curve is monotone with its minimum at λ=0; (d) the sign survives the μ
swap; (e) the median moves, not only the mean, and the effect is not carried by < 10 targets
(drop-top-10 reported against a uniform-effect null, per the standing warning).

## WHAT I WILL NOT DO
Tune weights, tune against RMSD, extend the λ ladder, re-pick μ after seeing RMSD, report an
alternative metric as primary, read the sealed benchmark, or soften a falsification.
