# SPRINT 15 — FINAL DOSSIER

The complete record of a research programme that set out to build a better peptide-folding
architecture, failed to, and in failing measured why — precisely enough to say what would be needed
instead.

**Bottom line, stated first.** The incumbent pipeline emits **3.204 Å**. The architecture built this
sprint emits **3.321 Å**, a significant loss (+0.117 [+0.050, +0.188]). Nothing measured here reaches
2.5 Å, and the sprint's own mechanism explains why the obvious routes cannot: **they correct error
magnitudes, and the obstacle is an error direction.**

What the sprint delivers instead is a set of measurements that change where effort should go — and,
on the quantum side, the finding that **every apparent effect dissolves under a classical control
using the same circuit's own samples at 1/200th of the budget**.

---

## 1. THE RESULTS THAT MATTER

Ordered by how much they change what someone should do next, not by how good they look.

### 1.1 The accuracy ceiling was a property of the library, not the information

Solving for a conformation by torsion-space distance geometry, instead of selecting one from a
retrieval library, reaches **0.611 Å mean / 0.055 Å median (ORACLE, from true distances)** on all 126
targets, where two standing project results put the ceiling at ≈1.95 Å. Those were measured with the
library in the loop.

**Perfect distance knowledge is worth 0.611 Å, not 1.95 Å.** Every prior estimate of the value of
improving distance prediction was computed against the wrong floor and roughly **tripled the floor** (1.95 / 0.611 = 3.19×) and **doubled the headroom** a
distance-prediction improvement can address ((3.204 − 0.611)/(3.204 − 1.95) = 2.07×). *(The counterweight travels with it: from predicted distances the same machinery gives
3.644 Å and loses to the incumbent by +0.440 [+0.290, +0.592].)*

### 1.2 Error shape beats error magnitude — in both directions, by more than an ångström

On a controlled phase diagram where the error's structure is **imposed** rather than observed:

| channel | effective error RMS | emitted RMSD |
|---|---|---|
| outlier-shaped (5% of pairs badly wrong) | **4.12 Å** | **1.993 Å** |
| i.i.d. Gaussian | 3.00 Å | 2.561 Å |
| **the real distogram** | **3.70 Å** | **3.644 Å** |

A *larger* error concentrated in a few pairs is far *cheaper* than a smaller one spread evenly, and
both are far cheaper than the real predictor at comparable magnitude. Making the error grow with
sequence separation costs essentially nothing once magnitude is matched (2.95 against 2.87 Å for the
2.5 Å target) — a clean null that retires an obvious hypothesis.

The same phenomenon runs the **other way** on the torsion channel: real torsion errors are **0.6–1.7 Å cheaper** than i.i.d. errors of the same magnitude across
**native-free** arms, and up to **2.9 Å cheaper for the ORACLE best-matching pool window**. Alignment is
**0.565 (ORACLE window)** and **0.665 (the incumbent's projected torsions)** against a **0.945**
random-direction null — but the control that matters is that a **zero-information constant α-helix
already scores 0.709**, so only those two arms clear it. Sign-flipping the ORACLE window's errors costs
**+1.596 Å [+1.385, +1.813]**; for the native-free retrieval mean the same operation is **null**
(−0.150 [−0.366, +0.070]).

### 1.3 The mechanism: the distogram's errors describe a coherent WRONG STRUCTURE

| distance set handed to the fit | unrealizable part |
|---|---|
| the **real** predicted distances | **0.977 Å** |
| the retrieval pool's distances | **0.458 Å** |
| ORACLE i.i.d. Gaussian of the same RMS — the null | **2.369 Å** |
| ORACLE the same errors permuted across pairs | 2.081 Å |
| ORACLE the same magnitudes with random signs | 2.070 Å |
| **ratio real/null** | **0.413**, paired **−1.391 [−1.623, −1.184]** |

**The real prediction is 2.4× closer to being realizable by an actual conformation than
matched-magnitude noise**, and both structure-destroying surrogates move it back to within 13% of the
null. About **59% of the error is geometrically realizable**: the predicted distance matrix is a good
description of a wrong structure, and a faithful optimiser builds it. The negative correlation
completes it — the part the fit *cannot* satisfy is the part anti-correlated with truth, so the
optimiser **discards the incoherent component and follows the coherent one**.

### 1.4 This unifies every negative the programme has produced

Every intervention tried acts on the error's **magnitude profile**; a realizable error is a
**direction** and lies in the null space of all of them. Six independent confirmations:

| intervention | worth |
|---|---|
| recalibrating the `sd` (2.633× over-confident) | **provably zero** — a weighted least-squares argmin is invariant under uniform rescaling of the weights |
| leave-fold-out separation debiasing | little on both axes; **negative on ranking**, at the noise floor on the fit |
| outlier clipping / robust losses | 0.12 Å (ORACLE) |
| per-target rescaling of the restraints | **0.1 Å ORACLE ceiling**; the pool estimates the scale at r = +0.035 |
| rescaling the consensus before projection | **0.275 Å ORACLE ceiling**, and every native-free estimator is uncorrelated with the truth (|r| ≤ 0.106) |
| channel fusion | 0.17–0.22 Å, and **predicted in advance** by §1.5 |

### 1.5 ~~A parameter-free, native-free law for what fusion is worth~~ — WITHDRAWN, Sprint 16

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT workstream. This subsection is withdrawn as a
> novelty claim. The original text is preserved below, struck through.**
>
> 1. **Prior art.** The "law" is the two-member **Krogh–Vedelsby (1995) ambiguity decomposition**
>    (`s16/lit_FINDINGS.md` §B.5). It must be cited, not claimed.
> 2. **Wrong mean.** `r` was computed as the **arithmetic** mean of the two channels' RMSDs; the
>    identity needs the **quadratic** mean. Corrected: prediction ~~3.205~~ → **3.292 Å**, error
>    ~~+0.162~~ → **+0.075 Å** (median +0.019); paired **−0.087 [−0.119, −0.061]** i.i.d.-target,
>    **[−0.129, −0.056]** fold-clustered, **W/L 126/0**, n = 126, TARGET as the unit
>    (`s16/retract_law.py`). **54% of the published prediction error was arithmetic, not physics.**
> 3. **The residual is a frame convention, not a prediction error.** In a single common frame with
>    the quadratic mean, the identity is exact to **2.65e−15** (max 2.52e−14, sign balanced 55/61)
>    on all 126 targets (`s16/retract_exact.py`, which reproduces the Sprint 15 structures
>    bit-identically). The surviving +0.075 Å is **+0.174** (averaging in the medoid frame rather
>    than the target's) **−0.099** (`s` Kabsch-minimised rather than measured in that same frame).
> 4. **The native-free screen is refuted, not merely unproven.** `s` as a ranker of "does fusion win
>    on this target" has **AUC 0.401 against a 0.500 null** — it points the wrong way; a supervised
>    leave-fold-out threshold degenerates to "always fuse" on 5/5 folds (`s16/retract_screen.py`).
>
> **What survives:** the *measurement* that two channels whose per-pair errors correlate at +0.602
> build structures 3.138 Å apart, and the `s/r` scaling statement below (which is algebra and is
> correct). Full account: `s16/retract_FINDINGS.md`.

Two channels whose per-pair errors correlate at **+0.602** nonetheless build wrong structures **3.138 Å apart**. Error correlation is the wrong statistic; what matters is whether the realizable components
point at the same wrong structure.

> ~~**d_avg ≈ √(r² − (s/2)²)**  — predicted 3.205 Å, measured 3.367 Å, **error +0.162 Å**, n = 126, no free parameters~~
>
> **CORRECTED:** `d_avg = √((r₁²+r₂²)/2 − (s/2)²)`, exact — predicted **3.292 Å**, measured 3.367 Å, residual **+0.075 Å** (a frame convention, n = 126)

~~`s` is measurable **without labels**, so the fusion decision becomes a calculation on unlabelled data.~~
**REFUTED** — the *prediction* needs `r`, which is a distance to the native; and `s` alone has AUC 0.401 as a fusion screen.
The gain is small whenever `s ≪ r` and large only as `s → 2r`; here `s/r = 0.88`. **An ångström from
fusion needs channels that disagree about twice as much as they are wrong**, and no pair in this
project comes close.

### 1.6 The variational comparison is not robust to how it is reported

Against best-of-N sampled from the **same circuit at its untrained parameters** — a control the
literature survey found **no paper running**, and strictly stronger than the uniform-random sampling
the nearest published comparison uses — the sign of the result depends on choices routinely left
unstated:

| dimension | what changes |
|---|---|
| **budget** | VQE **wins** at 8,192 (−0.202, 13/19, CI excluding zero); **0 of 6 cells significant at 32,768** |
| **readout** | wins on argmin; **loses significantly** on diversity-respecting set aggregation (+0.372, +0.392 at α = 0.25) |
| **n** | three of the workstream's own reads reversed between n = 1, n = 9 and n = 19 |
| **objective** | the win is **objective-independent** — its strongest cell is the *selective* objective with the worst measured in-tail ordering, so it is not exploiting ordering skill at all |
| **classical control** | **0 wins in 16 against greedy**, losing 3 |

> **And the quantity that dominates every arm — a selection gap of 1.47–2.14 Å — is identical for the
> quantum and control arms, flat across a 16× budget range, and collapses to 0.29–0.44 Å under ORACLE
> restraints. It is restraint error, not the sampler.**

Sprint 14's mechanism reproduces and is **not** a selective-class pathology: `tail − bulk` ordering
skill is negative for every objective, including all four generative ones and the ORACLE.

**Competing explanations eliminated by measurement, not assumption:** not barren plateaus
(`Var ~ 2^(−0.18 to −0.36 n)`, every CI excluding −1.0); not ill-conditioning (the gradient's
null-space share is 0 to 1e-30 and it avoids the small eigendirections); not a broken estimator (CVaR
verified to 4.2e-14); not shot noise.

**And the recorded CVaR gradient defect *helps*.** The biased tail-centred baseline gives −0.266
where the corrected one gives −0.145, because its gradient norm is **2.4× smaller** — a de-facto
step-size reduction. **A weaker optimiser is a better sampler**, because it concentrates less, and
concentration is what costs the distinct samples a best-of-N readout is paid in. A stated liability
inverts into a mechanism.

**What concentration is, measured directly:** at α = 1 the drawn set is **0.847 Å better on the mean
and 1.092 Å worse on the best** — two significant effects of opposite sign, both scaled by α. α is a
**diversity dial** (2.7 → 408 distinct configurations; entropy 0.73 → 6.58 bits), not an accuracy
dial.

### 1.7 Alignment is a description, not a lever — and it closes the route to 2.0 Å

The information-channel work showed that a channel's cost is set by **where its error points**
relative to the map consuming it. The obvious follow-up is whether that direction can be *engineered*
rather than inherited. It was tested across 22 arms on 126 targets.

| | |
|---|---|
| **ORACLE**: rotate the fit's own true error out of the loud half at **fixed magnitude** | **1.855 Å** — −1.822 [−2.037, −1.613], W/L **121/5** |
| best **native-free** leave-fold-out arm, of 22 | **−0.007 [−0.074, +0.055]** |
| alignment actually moved by any arm | ≤ **0.104**, against the **0.387** required |
| pooled partial corr(ΔRMSD, Δalignment given Δraw error), 2,142 pairs | **+0.087**, against β = +0.586 for distance MAE |

> **The information is present in the fit's own error vector — an ORACLE rotation clears 2.0 Å from
> the *same* restraints at the *same* error magnitude — and the best native-free rotation captures 4%
> of it.**

So **2.0 Å remains unsupported, but the reason has changed from information to control.** Family D,
declared unallocated and to be filled only if this worked, **stays empty**.

Three things fell out of it, each refuting something the sprint believed:

- **The geometry was more interesting than expected.** There are **five** exactly-null torsion
  directions, not four — the fifth a *distributed whole-chain crank* (Σδφ ≈ −Σδψ), which no
  per-torsion perturbation could find, and the count is structural: `2n − 4` parameters onto a
  `2n − 5` shape space. The Cα trace responds to a participation ratio of **3.47 effective directions
  out of ~25.9**.
- **The terminal-gap lever does not transfer.** *Perfect native* terminal torsions are worth
  **−0.000 Å [−0.068, +0.066]** while cutting torsion RMS by 17°; dropping terminal restraints
  **costs +0.210**. The 1.79 Å spread is a property of channels with holes, not of a solver.
- **Robust losses are refuted at n = 126** (−0.011 and −0.025), firing the module's own pre-declared
  falsification condition and closing Family B for the second time. Most tellingly, the robust arms
  **do** reach far better restraint residuals (median |z| 0.364 against 0.469) and **convert none of
  it into accuracy**.

**And one genuinely open lead.** A native-free error-direction surrogate **exists**: `θ_fit − θ_pool`
has |cos| **0.390** with the true error against a 0.157 null, and the quiet subspace is itself
identifiable native-free (subspace overlap cos² **0.827** against a 0.499 null). Both ingredients an
alignment intervention needs are separately estimable without the native. **It is untested as a
steering signal**, and it is the cleanest thing this sprint leaves behind.

---

### 1.8 One thing that does work: per-target confidence, native-free

**The programme cannot make these predictions better, and it can tell you which of them to believe.**
Leave-fold-out across the five pinned folds, using only quantities the pipeline already computes
(mean predicted `sd`, retrieval pool spread, top-75 BLOSUM similarity):

| | |
|---|---|
| ρ with the fit's own RMSD | **+0.665** (Pearson +0.629) |
| **best quartile** (n = 32) | **2.354 Å** |
| **worst quartile** (n = 32) | **4.884 Å** |
| **spread** | **2.530 Å** |

Even the mean predicted `sd` **alone** gives ρ = +0.617 and a 2.338 Å quartile spread. It improves no
prediction — it sorts them — but it is free, and it is the natural companion to everything above.

**And the control that shaped it.** I proposed two *new* signals from this sprint's own machinery —
the fit's restraint residual (ρ +0.604) and the channel disagreement (ρ +0.607) — and ran the
comparison before writing them up. Adding both to the baseline moves ρ from **+0.665 to +0.664**.
**They add nothing**, and the single simplest quantity already beats either. The proposed
contribution is refuted; the capability survives.

**Not to be confused with in-band selection.** Ordering *targets* by difficulty is a far easier
problem than ordering *candidates within a target's pool*, which fifteen sprints have failed at. The
predictors here are properties of the **restraint set**, available before any structure is chosen.
Nothing here ranks candidates.

---

### 1.9 The quantum question, closed: the effect is real, and it is not quantum

An apparent positive — VQE beating best-of-N when its ensemble is consumed **unranked** — was
replicated at power (**19 targets, 8 seeds, 1,672 new cells**, after reproducing all 108 original
cells **bit-for-bit at maximum difference 0.000e+00**).

**The effect is real and replicates out of sample.** At the target-level unit it is
**−0.385 [−0.748, −0.035]** at α = 1 and **−0.295 [−0.495, −0.069]** at α = 0.25 — where the
original found nothing — and it is **stronger on the ten targets never used to find it** than on the
original nine. A 10-point α map is smooth, monotone and significant at 9 of 10, so the
non-monotonicity I had cited as a reason to disbelieve was **two noisy three-seed cells**.

**And it is not a quantum effect.**

> **A classical Boltzmann reweighting of the *same untrained circuit* by the *same objective* at
> matched entropy, using 1/200th of the objective budget, beats the CVaR-VQE by +0.25 to +0.48 Å on
> every readout at every α** (W/L up to 1/18). Simulated annealing at **1/400th** of the budget
> matches or beats it. The VQE beats only uniform random — which proves nothing here, since a
> zero-information constant α-helix beats that control too.

On the optimisation axis the gap is starker: **simulated annealing finds the certified global optimum
in 100% of cells at 2,048 evaluations**, while the VQE puts **5.9% of its mass** on it using
**819,200**.

**What the effect is.** Location, entirely: `d_coordavg = 0.94 · d_set_mean` at R² = 0.75, with no
residual shape term. And the per-target win correlates **+0.52 / +0.64 / +0.69** with where the
objective's own certified argmin sits in the ORACLE RMSD distribution — the method is rewarded on
targets where the *objective* happens to be right, and loses by +1.25 Å where it is not. Matched
diversity removes what remains, to at or below the 0.08 Å floor.

**The sub-2 Å enrichment result collapses identically.** It replicates exactly (5.20×) — and
classical annealing at the *same* budget gives 5.34×, null against it and significantly better on the
ORACLE objective. The headline is also a ratio of means: the **median per-target ratio is 1.13–2.66×,
with only 5 of 9 targets showing any enrichment**.

**A prediction of mine, refuted.** I proposed that the win should track distinct-configuration count,
since concentration costs the distinct samples a best-of-N readout is paid in. Within-α
**ρ = +0.004** across 1,216 pooled cells, and the aggregate direction is the *opposite*.

> **There is no quantum result here, positive or negative, that survives its classical control. Every
> apparent effect dissolves under a control using the same circuit's samples, the same objective, and
> two to three orders of magnitude less budget. The distribution is worth something; the optimiser is
> not the cheapest way to produce it.**

---

## 2. THE ARCHITECTURE, AND WHERE IT FAILED

| stage | mean | vs incumbent |
|---|---|---|
| **G** best in the generated ensemble (**ORACLE**) | 2.760 | **−0.444 [−0.573, −0.325]** |
| **S** the objective's argmin | 3.511 | +0.307 [+0.178, +0.439] |
| **A** coordinate consensus | 3.151 | −0.053 [−0.125, +0.020] |
| **F** after projection — *the like-for-like number* | **3.321** | **+0.117 [+0.050, +0.188]** |

**The ensemble is good and the readout is bad.** Better structures than the incumbent emits are
generated on most targets (49% of targets hold a member under 2.5 Å), and the entire deficit is in
reading them out: selection costs **+0.751**, aggregation recovers **−0.360**, projection costs
**+0.170**. Pre-filtering the ensemble by the objective makes aggregation *worse*, closing the
obvious response.

Family status: **A** (distance-restrained torsion VQE) alive, quantum arm in flight. **B**
(constrained/augmented-Lagrangian) **falsified as declared** — at a one-sigma band fewer than half the
restraints admit the native — and redirected to robust potentials. **C** (conditional ensemble
generation) alive and measured above. **D** (reserve) **deliberately unfilled**, pending the alignment
result.

---

## 3. THE POSITIVE CONTROL

Without any search, the restraint objective's argmin over a 500-member pool beats a random member by
**−0.949 Å [−1.147, −0.753]** on 79% of targets, ρ(objective, RMSD) = +0.562 positive on 88.1%. The
ORACLE control verifies the machinery: the native at percentile **0.000 on 126/126**.

At the same chain lengths where a published lattice contact potential's minimum-cost conformation is
**worse** than random, a distance-restraint objective's minimum is about an ångström **better**. That
pathology is a property of contact potentials, not of short-chain prediction.

The native itself nonetheless sits at the **34.7th percentile** and is the argmin on 4 of 126 targets
— independently reproducing a standing result (36.8th percentile, 3 of 126) through different code.

---

## 4. WHAT WAS WRONG, AND WHO FOUND IT

The programme's retraction rate is reported because it calibrates everything else.

| # | claim | fate |
|---|---|---|
| 1 | an 8-target read putting the predicted arm at 2.792 Å | superseded — the full instrument gives **3.644 Å and reverses the sign** |
| 2 | a "1.06 Å quantisation ceiling" | **RETRACTED same day** — a uniform-grid bin lookup on a non-uniform grid, my own bug |
| 3 | "coordinate averaging contracts the backbone by 25.8%" | **WITHDRAWN** — measured at **3.5%**; I reused it in four places without checking |
| 4 | multi-start reproducibility | **broken** — `hash()` is salted per process; found by an agent, outside its brief |
| 5 | every AMBER row in one Sprint 14 artefact | **RETRACTED** — no stratification mask, yet flagged complete |
| 6 | the brief's oracle-conditioning share | **transposed** — 21.7% oracle, 40% unbiased |
| 7 | `amber_kind == 0` is an unbiased stratum | **not in its tails** — an ORACLE index is force-included |
| 8 | `I.FAIL18` as an object | **a threshold artefact** — only 1 of 18 targets survives every threshold |
| 9 | Sprint 14's metric condition numbers | **a seed-averaging artefact** — per point the median is 735, max 1.5e6 |
| 10 | six of the quantum workstream's own claims, one of the info workstream's headline | **self-refuted** |

**Two of the three coordinator retractions were caught by checking machinery against its own stated
assumptions, not by scientific intuition — and in both cases the wrong result was more plausible than
the right one.** That is the single most transferable lesson here.

A methodological trap worth the same weight: **the budget convention decides the QNG conclusion** —
equal iterations gives 5/5 wins, equal hardware cost gives 0/7, with clean sign flips. *"Reporting
only convention A would have produced a false positive."*

And one that invalidates a law this project uses constantly: **the set-mean proxy is out of domain on
a collapsed set.** A method that concentrates improves the proxy and degrades the output.

---

## 5. THE REQUIREMENT, STATED

| target | requirement on the distance channel |
|---|---|
| 2.5 Å | effective RMS ≤ **2.87 Å** for an i.i.d. channel — a **22% reduction** from the current 3.70 Å |
| 2.0 Å | effective RMS ≤ **1.87 Å** — roughly **halving** the current error |

And ours is **~1.1 Å worse than i.i.d. at equal magnitude**, so the real requirement is stiffer.

**It is not "make the distogram better."** It is: **reduce its RMSE by 22%, or make its errors
incoherent.** §1.2 shows an error that is *larger* but unstructured in the right way is worth more
than a smaller one — and nothing in this project, or as far as the literature survey found, anywhere
else, has tried to shape a predictor's error rather than shrink it.

### 5.1 And the second option — measured, and harder than the first read suggested

*(n = 40, a partial of the full-instrument run; supersedes an n = 8 read that gave materially
different numbers and one claim that has been withdrawn)*

Interpolating between the real error and a sign-destroyed version of it, **at constant error
magnitude**, traces the dose–response:

| | t = 0 | 0.25 | 0.5 | 0.75 | t = 1 |
|---|---|---|---|---|---|
| **constant-magnitude ladder** — the experiment | **3.411** | 3.415 | 2.875 | 2.368 | **2.127** |
| *magnitude allowed to shrink* — the confound | *3.411* | *2.899* | *2.401* | *2.181* | *2.127* |

**Destroying the coherence is worth 1.284 Å with no improvement whatsoever in the error's size.** The
response is **not linear**: it is flat to t = 0.25 and then falls steeply, so partial decorrelation
buys nothing until a threshold is passed.

> **Reaching 2.5 Å requires destroying about 69% of the error's coherence** (t ≈ 0.69), at zero
> improvement in error magnitude.

**The confound is large and must be controlled** — 0.474 Å at t = 0.5 — because mixing two vectors
shrinks the result, and plain error reduction is exactly what this experiment is trying not to
measure.

**A withdrawn claim.** An n = 8 read of this ladder put the requirement at 30% and reported that the
uncontrolled ladder dipped *below both its endpoints* at t = 0.5 — "geometrically impossible as a
decorrelation effect", which I wrote up as a third instance of an experiment inverting without a
control. **At n = 40 the uncontrolled ladder is monotone and never dips below its endpoint.** That
reading does not exist in the data and is withdrawn; the control is still the right design and the
confound is still large, but the rhetorical point it supported was an artefact of eight targets.

---

## 6. NOVELTY, HONESTLY BOUNDED

**Ceded:** certified-optimum-is-worse-than-random is published for lattice contact potentials
(Roget et al.); in-loop CVaR for peptide folding is published (QuPepFold); tail-restricted
discrimination has a framework (BEDROC). Four Sprint 14 claims are removed from the story entirely.

**Survives, ranked:** the untrained-circuit best-of-N control (no precedent found); the CVaR estimator
defect triple; the realizability mechanism and its null control; the parameter-free fusion law; the
error-shape phase diagram; the ceiling reframe.

**Not claimed:** any quantum advantage; any comparison to prior predictive work (the nearest
comparator is 4.89 Å on a different target set — our 3.204 Å is **not** a like-for-like win); any
experimental validation.

---

## 7. STATUS AND WHAT REMAINS

**Complete.** The Phase 0 audit (gate PASS, with one retraction, one correction to the brief, one new
defect and four disclosures). The literature and novelty boundary (67 papers). The information-channel
audit. The quantum geometry and CVaR mathematics. The representation ladder. The generative-objective
quantum experiment. Alignment engineering. The replication of the claimed quantum positive — which
**demoted it**. An adversarial consistency audit that traced 178 quantitative claims and raised
**13 blockers, all applied**. The ceiling, mechanism, feasibility, cascade, phase-diagram, projection
and coherence experiments. **26 deliverable documents; 13 of 16 figures rendered**, the rest blocked
on two runs in flight and refused by an `n ≥ 100` guard with three declared exceptions, each printing
its own n.

**The AMBER replication, landed — and it splits the claim in two.** The **accuracy** half is
**WEAKENED**: it replicates in sign and magnitude (5/5 draws, sd 0.0076 Å) and is **not** a
multi-start artefact — that pipeline has no RNG and reproduces bit-identically across interpreters,
so **my application of the 0.08 Å floor to it was wrong** — but it is **absent on the
frame-reproducible 69% of the instrument** (a *losing* 38W/46L, three frames), **erased by an exact
rotation null that must be zero** (+0.0117, max 1.51 Å), and four targets are **silently
non-converged** at up to 8.9 × 10⁸ kcal/mol because `refine_coords` has no convergence gate. It may
no longer be quoted as "−0.022 Å at valid geometry". The **validity** half is **CONFIRMED and was
understated**, once two Sprint 14 reporting defects were fixed: **0.466 → 0.874
Ramachandran-favoured (116W/2L)** and **1.397 → 0.000 clashes (63W/0L)**. **AMBER's defensible role
is stereochemical repair, not accuracy.**

**Stopped deliberately, and recorded as such.** The full `robust.py` run (its hypothesis was refuted
at n = 126 inside another workstream, firing the module's own falsification condition) and the full
`scale.py` run (its ORACLE ceiling had already closed the direction at ~0.1 Å, and an independent
full-instrument measurement of the same idea agreed). Both modules are committed and runnable; both
smoke reads are labelled, and the figure module refuses to render them.

**Not started, by design.** The 60-target benchmark. It has **not been read** — independently
confirmed by the consistency audit, which found one module hashing its manifest's bytes without
loading any content, and no other access anywhere. The protocol and its decision rule are
pre-registered, **including the rule for whether to spend the benchmark at all**: if no intervention
qualifies, the frozen protocol is the incumbent unchanged, a confirmatory run could only return the
pre-registered `null`, and the instrument stays sealed for a future sprint with a real candidate.

**Nothing qualified, and the protocol is frozen as the incumbent unchanged.** Thirteen candidates
were tested against their own controls on the full instrument. Twelve failed outright — the cascade
(+0.117), consensus rescaling (+0.080 to +0.123, with every native-free scale estimator uncorrelated
with the truth), robust losses (−0.011, refuted at its own falsification condition), alignment
engineering (−0.007), terminal relaxation (+0.210), pool augmentation (+0.004, with the leave-fold-out
procedure choosing **w = 0 on 4 of 5 folds** and its premise falsified: the fits are 0.108 Å *worse*
than the set they would join), and the quantum ensemble readout (demoted, then beaten by a classical
reweighting at 1/200th the budget). Two cleared *their own* control — `tik_LFO` at −0.216 Å and
coordinate-level fusion at −0.255 Å — and both improve a **losing** arm without reaching the
incumbent.

**So the 60-target benchmark is not spent.** By a rule fixed in advance, a confirmatory run on an
unmodified pipeline could only return the pre-registered `null`, and the instrument stays sealed for a
sprint with a real candidate. The nearest miss, `tik_LFO`, needs re-basing on the incumbent and a
converged λ grid — a defined next experiment, not a hope.

---

## 8. THE HONEST SUMMARY

This sprint did not produce a better folding method. It produced an explanation for why fifteen
sprints of corrections bought tenths of an ångström, and the explanation is specific, measured, and
falsifiable: **the predictor's errors are mutually consistent — they describe a plausible wrong
structure, and every tool applied to them so far corrects magnitudes rather than directions.**

That reframes the engineering target from shrinking a predictor's error to **decorrelating** it, and
it explains the quantum result in the same breath: a variational optimiser concentrates, and
concentration only helps when the objective can order what it concentrates onto. Here it cannot — so
concentrating is strictly harmful, and the untrained circuit wins.
