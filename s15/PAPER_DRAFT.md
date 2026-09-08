# Restraint coherence, not restraint accuracy, limits short-peptide structure prediction — and the variational-optimiser comparison is not robust to how it is reported

*Draft. Every number is on the full 126-target instrument unless marked *(in flight)*, which report a
smoke read with its n inline. Every ORACLE quantity is labelled as such wherever it appears.*

---

## ABSTRACT

We study backbone structure prediction for peptides of 9–16 residues, a regime where the statistical
potentials and scoring functions calibrated on globular proteins sit outside their fitted domain. We
report six results, five of them negative.

**First**, we show that every apparent effect of a CVaR variational quantum eigensolver on this
problem **dissolves under a classical control that uses the same circuit's samples, the same
objective, and two to three orders of magnitude less budget**. A Boltzmann reweighting of the
*untrained* circuit at matched entropy, costing 1/200th of the objective evaluations, beats the VQE
by **+0.25 to +0.48 Å on every readout at every α**; simulated annealing at 1/400th matches or beats
it, and finds the certified global optimum in **100% of cells at 2,048 evaluations** where the VQE
places 5.9% of its mass using **819,200**. An apparent VQE advantage does replicate out of sample
(−0.385 [−0.748, −0.035] at 19 targets and 8 seeds, stronger on the ten targets never used to find
it) — but it is pure **location**, it is explained by *where the objective's own optimum happens to
lie* (ρ = +0.52 to +0.69 per target), and it vanishes under matched diversity.

We also show that whether a VQE beats **best-of-N sampled from
the same circuit at its untrained parameters** — a control we could not find in the literature, and
strictly stronger than the uniform-random sampling the nearest published comparison uses — depends on
choices that are routinely left unstated. It **wins** at a budget of 8,192 evaluations (−0.202,
13/19, CI excluding zero) and the effect is **gone at 32,768** (0 of 6 cells). The win is
**objective-independent**: its strongest cell is on the *selective* objective with the worst measured
in-tail ordering, so it is not exploiting the objective at all. On diversity-respecting set readouts
it **loses** significantly. And it **never beats a classical greedy search** (0 of 16). Across every
arm the dominant quantity — a selection gap of 1.47–2.14 Å — is **identical for the quantum and control
arms** and flat across a 16× budget range, collapsing to 0.29–0.44 Å under oracle restraints: it is
restraint error, not the sampler. We eliminate barren plateaus, metric ill-conditioning and estimator
error as explanations by direct measurement, and find that a **recorded defect in the CVaR gradient
baseline *helps*** — its gradient norm is 2.4× smaller, so it acts as a step-size reduction, and a
weaker optimiser is a better sampler.

**Second**, we show that the accuracy ceiling long assumed for this problem was a property of the
*retrieval library*, not of the information channel. Solving for a conformation by torsion-space
distance geometry, rather than selecting one from a library, reaches **0.611 Å mean / 0.055 Å median
(ORACLE, from true distances)** where the library-mediated ceiling was ≈1.95 Å. That **triples the floor** (3.19×) and **doubles the addressable headroom** (2.07× against the
incumbent).

**Third**, and centrally, we show that **the shape of a predictor's error matters more than its
magnitude, and in both directions**. On a controlled phase diagram, outlier-shaped distance error at
4.12 Å RMS emits 1.993 Å while i.i.d. error at 3.00 Å RMS emits 2.561 Å — a *larger* error is
*cheaper*. A learned distogram at 3.70 Å RMSE emits 3.644 Å, worse than either, because **59% of its
error is geometrically realizable**: the predicted distance matrix is a good description of a *wrong
structure* rather than a noisy description of the right one, and it is 2.4× closer to realizable than
matched-magnitude noise (n = 126, paired −1.391 [−1.623, −1.184]). The optimiser reproduces the wrong structure faithfully. The same phenomenon
runs the *other* way on the torsion channel, where **native-free** real errors are 0.6–1.7 Å
**cheaper** than i.i.d. (up to 2.9 Å for an ORACLE arm) because they lie in the Jacobian's
low-response subspace — though a constant α-helix already scores 0.709 on the alignment statistic, so
only the ORACLE window (0.565) and the incumbent's projection (0.665) clear that control.

~~**Fourth**, we give a parameter-free, **native-free** law for what combining two channels is worth:
for fits at mean distance `r` from the truth and `s` from each other, their average sits at
`√(r² − (s/2)²)`, predicting the measured value to **0.162 Å across 126 targets**.~~
**⚠ WITHDRAWN — Sprint 16, 2026-09-06 (RETRACT). This sentence must not appear in any draft.** The
expression is the two-member **Krogh–Vedelsby (1995) ambiguity decomposition** (prior art); it was
computed with the arithmetic rather than the quadratic mean of `r` (corrected residual **+0.075 Å**,
not 0.162); it predicts nothing unseen, being exact to **2.65e−15** in a common frame on all 126
targets; and `s` as a native-free fusion screen has **AUC 0.401 against a 0.500 null**. See §6 and
`s16/retract_FINDINGS.md`. **Fourth**, what survives is a *measurement*: two channels whose per-pair
errors correlate at +0.602 nonetheless build wrong structures 3.14 Å apart, so error correlation is
the wrong statistic for the fusion decision.

**Fifth**, we calibrate an **empirical false-positive floor** by running a comparison that is zero by
construction: it returned **+0.081 Å [+0.014, +0.169]** on one random restart draw and −0.003 on
another — a 95% interval excluding zero on a null effect. Absolute arm means vary with sd 0.132 Å
across draws. We apply the floor retrospectively to three of our own results — and record where it does **not**
transfer: a deterministic pipeline with no random restart is unreachable by a floor whose mechanism
is a per-process-salted multi-start, a limit we established only after misapplying it once.

**Finally, one capability that does work.** Per-target confidence is achievable **native-free**: a
leave-fold-out model over quantities the pipeline already computes reaches ρ = **+0.665** with the
method's own error, separating a best quartile at **2.354 Å** from a worst quartile at **4.884 Å**.
It improves no prediction — it sorts them. We report it with the control that shaped it: two signals
newly derived from our realizability analysis correlate at +0.604 and +0.607 individually and, added
to the baseline, move ρ from +0.665 to **+0.664**. They add nothing.

The architecture built to exploit these findings **does not beat the retrieval baseline** (3.321 Å
against 3.204 Å, +0.117 [+0.050, +0.188]). We report the requirement instead: reaching 2.5 Å needs a
22% reduction in distance-prediction RMSE for an i.i.d. channel — or, from a dose–response curve at
*constant* error magnitude, the removal of about **69% of the error's coherence**, which is a
different and more tractable objective than any predictor is currently trained for.

---

## 1. INTRODUCTION

Predicting the backbone of a short peptide looks like a small version of protein structure
prediction. It is not, and the ways it differs are the subject of this paper.

At 9–16 residues a peptide is conformationally heterogeneous, often lacks a single dominant fold, and
— decisively — sits outside the fitted domain of nearly every scoring function the field uses. Roget
et al. make this quantitative for lattice contact potentials: enumerating every peptide up to length
15 across 12,446 PDB structures, they find the **minimum-cost conformation has, on average, a larger
RMSD than a randomly selected feasible one**, because the Miyazawa–Jernigan potential was
parameterised on proteins above 50 residues. Optimising the objective moves away from the answer.

That result frames the problem but does not settle it, because it is about one potential class. We
find the opposite sign for a distance-restraint objective at the same chain lengths: its argmin over a
500-member candidate pool is **−0.949 Å [−1.147, −0.753]** *better* than a random member, winning on
79% of targets. So the pathology is a property of contact potentials, not of short-chain prediction —
and the interesting question is what *does* limit the problem.

**Two questions are usually conflated, and separating them is this paper's main contribution.** The
first is how much information a channel carries. The second is what *shape* that information's error
has. The field measures the first, routinely, as an accuracy: an RMSE, an MAE, a calibration curve.
We show the second dominates.

Concretely, on a controlled phase diagram in which the error's structure is imposed rather than
observed, a distance channel with an **outlier-shaped** error of 4.12 Å RMS produces structures at
1.993 Å, while an **i.i.d.** channel with a *smaller* 3.00 Å RMS produces 2.561 Å. A real learned
predictor at 3.70 Å RMSE produces 3.644 Å — worse than either. The ordering
`outliers ≪ i.i.d. < real` at comparable magnitude is the whole argument: a fit can ignore a few
catastrophic restraints, cannot ignore many mildly wrong ones, and is actively misled by errors that
agree with one another.

The mechanism for that last case is measurable. About **59% of the learned predictor's error is
geometrically realizable** — its predicted distance matrix is a good description of a *wrong
structure* rather than a noisy description of the right one, and it sits **2.4× closer to realizable**
than matched-magnitude noise. A faithful optimiser therefore builds the wrong structure. This
unifies, from one measurement, why every correction we and others have tried buys tenths of an
ångström: recalibration, debiasing, robust losses, outlier trimming and rescaling all act on the
error's **magnitude profile**, and the expensive part is a **direction**.

We also revisit a ceiling. The belief that perfect distance knowledge caps this problem near 1.95 Å
turns out to be a property of the *retrieval library* through which distances were consumed, not of
the distance channel: solving directly for a conformation reaches **0.611 Å** from true distances.
That triples the floor (3.19×) and doubles the addressable headroom (2.07×) — though we are careful
that this is an oracle ceiling, and our own predicted arm loses to the baseline.

Finally, we treat the quantum component as a scientific variable rather than a framing device. The
controlling comparison for any variational optimiser is **best-of-N sampled from the same circuit at
its untrained parameters** at matched budget — a control we could not find in the literature, and
strictly stronger than the uniform-random sampling the nearest published comparison uses. Running it
across a grid rather than at a point, we find the answer is **not a sign but a dependence**: the
result flips with budget, with readout, and with n, and never survives comparison to a classical
greedy search. The quantity that dominates every arm is identical in the quantum and control arms.

**This paper reports a method that does not work.** Our architecture emits 3.321 Å where the baseline
emits 3.204 Å. We report it because the measurements that establish *why* are more useful than the
method would have been, and because three of them — an empirically calibrated false-positive floor, a
~~parameter-free law for what channel fusion is worth~~ (**WITHDRAWN 2026-09-06, Sprint 16 RETRACT — prior art, wrong mean; see §6**), and a reporting-dependence in the quantum
comparison — are, we think, transferable well beyond peptides.

---

## 2. METHODS

*[summarised; full detail in METHODS.md]*

126 tuning targets, five pinned folds, a 60-target benchmark **held out and unread** until the
protocol is frozen. RMSD frozen as Cα-only, full chain including both termini, proper rotations only,
model 1 of the native — independently reimplemented by Horn's quaternion method and agreeing with the
production Kabsch to **2.04e-13 Å over 63,000 structures**. Full-chain versus `[1:-1]` differs by up
to 2.4 Å on one structure here, so the convention is not cosmetic.

Conformations are continuous backbone torsions under ideal geometry, with an exact **O(n)
reverse-mode gradient** verified at **cosine 1.000000000000000** against central differences. Four
torsions per chain are inert for the Cα trace, so the live parameter count is `2n − 4` and the live
qubit register is `(n − 2)·log₂k`.

Every comparison is paired across targets with bootstrap 95% intervals, fold-aware. Multi-start
selection uses the **objective only**; RMSD is read post hoc. ORACLE quantities — anything derived
from native coordinates — are labelled in every table row and never appear as predictive results.

---

## 3. THE VARIATIONAL COMPARISON IS NOT ROBUST TO HOW IT IS REPORTED

### 3.1 The control

Best-of-N sampled from the ansatz at its **untrained** parameters, at a matched hard budget of
objective evaluations. This is stronger than the field's usual control and, to our knowledge,
unreported.

### 3.2 The result depends on the budget, the readout, and n

On 19 enumerated targets with a matched hard budget, **positive means VQE is worse**:

| arm, budget 8,192 | paired difference | wins |
|---|---|---|
| `E_ml_pred` (generative) | **−0.202** | 13/19 |
| `E_ls_pool` (generative) | **−0.184** | 14/19 |
| **`S14_disto_bayes` (SELECTIVE control)** | **−0.192 — the strongest cell** | **15/19** |

| budget | cells significant, of 6 |
|---|---|
| 2,048 | 1 |
| 8,192 | 3 |
| **32,768** | **0** |

Two things follow immediately. **The win is objective-independent** — its strongest cell is the
objective with the *worst* in-tail ordering measured anywhere here, so it is not exploiting ordering
skill. And it is **gone at the largest budget**. A single-budget report would have supported any of
three different conclusions.

**Against a classical greedy search at matched budget, VQE wins 0 of 16 and loses 3.**

**And the mechanism reproduces regardless.** `tail − bulk` ordering skill is negative for *every*
objective — all four generative, the selective control, and the ORACLE. The generative objectives do
have better **in-tail** ρ (+0.127, z = 13.1, against −0.042 for the selective energy), and it still
does not become *positive* relative to the bulk. **Globally the two classes are close** — ρ = +0.523
for `E_combined` against +0.485 for the selective control — so the generative class does not have a
global-ordering advantage worth the name.

> **The binding quantity is a selection gap of 1.47–2.14 Å, identical for the quantum and control
> arms, flat across a 16× budget range, and collapsing to 0.29–0.44 Å under ORACLE restraints. It is
> restraint error, not the sampler.**

### 3.2b The set-level experiment, which goes the other way

With **exact statevector gradients** (the most favourable setting available), 12-qubit sub-registers
of nine enumerated targets, 4,096 draws, paired over 27 cells per α — positive means VQE is worse:

| readout | α = 1.0 | α = 0.25 | α = 0.05 | α = 0.01 |
|---|---|---|---|---|
| argmin by objective | −0.036 | **+0.251** | +0.068 | +0.048 |
| top-20 coordinate average | +0.068 | **+0.372** | +0.197 | +0.177 |
| top-75 coordinate average | +0.147 | **+0.392** | **+0.276** | +0.164 |
| drawn-set mean | **−0.847** | **−0.372** | **−0.191** | −0.066 |
| drawn-set best | **+1.092** | **+0.345** | **+0.155** | +0.087 |

Bold = significant. **No entry is a significant win for VQE on a structural readout.**

### 3.3 What concentration is

At α = 1 the drawn set is **0.847 Å better on the mean and 1.092 Å worse on the best** — two
significant effects of opposite sign and near-equal size, both shrinking monotonically to zero as α
falls. α is a **diversity dial**: 2.7 → 408 distinct configurations out of 4,096 draws, entropy
0.73 → 6.58 bits.

### 3.4 Explanations eliminated, not invoked

- **Not barren plateaus** — gradient variance and scaling measured directly.
- **Not ill-conditioning** — the gradient's null-space share is 0 to 1e-30 (a theorem, confirmed) and
  the bottom eigenvalue decile carries 0.0016–0.0043 of its squared norm where uniform would be 0.10.
  The ill-conditioning is real and *the gradient does not point into it*.
- **Not a broken estimator** — CVaR verified to 4.2e-14 against Rockafellar–Uryasev.
- **Not shot noise** — exact gradients.
- **And the one recorded estimator defect *helps*.** The biased tail-centred baseline gives −0.266
  where the corrected one gives −0.145, because its gradient norm is **2.4× smaller**: it is a
  de-facto step-size reduction, and a weaker optimiser is a better sampler. A second recorded defect
  is inert (0.000 across all runs).

We also correct two published-adjacent facts: **QNG and classical natural gradient are the same
algorithm here** (classical Fisher = 4 × Fubini–Study to 3.3e-16), and previously reported metric
condition numbers of 2.2–3.7 are a **seed-averaging artefact** — per point the median is 735 and the
maximum 1.5e6, and one common ansatz is *exactly rank-deficient* at depth ≥ 2.

### 3.5 A methodological trap

The standing law that a pipeline's terminal operator consumes the *set mean* is **out of domain on a
collapsed set**. At α = 1 the VQE returns two distinct structures out of 4,096 draws, so its "set
mean" looks 0.847 Å better while its actual coordinate average is worse. **A method that concentrates
improves the proxy and degrades the output.** Every set-aggregating arm in this work reports its
diversity alongside its RMSD.

---

## 4. THE CEILING WAS THE LIBRARY

Fitting torsions to the **true** distance matrix (ORACLE) through machinery identical to the
predictive arm:

| arm | mean | median | <2 Å | vs baseline |
|---|---|---|---|---|
| **ORACLE**, true distances | **0.611** | **0.055** | 0.86 | −2.593 [−2.901, −2.295] |
| predicted, 1/sd² weighted | 3.644 | 3.466 | 0.16 | +0.440 [+0.290, +0.592] |

A 0.055 Å median means ideal-geometry torsions reproduce the native trace from its own distances to a
rounding error on half of targets. Since both arms share optimiser, starts, selection rule and
manifold, **the 3.03 Å between them is entirely restraint error**.

Prior estimates put this ceiling at ≈1.95 Å — but those were measured with a retrieval library in the
loop, using distances to *select* among library members. That is the ceiling of selection through a
finite library. **Removing the library moves it to 0.611 Å.**

---

## 5. ERROR SHAPE BEATS ERROR MAGNITUDE

### 5.1 A controlled phase diagram

True distances corrupted under three error models at eight severities, on one comparable axis (the
effective RMS of the injected error):

| channel | effective RMS | emitted RMSD |
|---|---|---|
| outlier-shaped (5% badly wrong) | **4.12 Å** | **1.993 Å** |
| i.i.d. Gaussian | 3.00 Å | 2.561 Å |
| **the real distogram** | **3.70 Å** | **3.644 Å** |

**A larger error concentrated in a few pairs is far cheaper than a smaller error spread evenly, and
both are far cheaper than the real predictor at comparable magnitude.** Requirements: 2.5 Å needs an
i.i.d. RMS of **2.87 Å**; 2.0 Å needs **1.87 Å**. Making the error grow with sequence separation
costs essentially nothing once magnitude is matched (2.95 versus 2.87) — a clean null that retires an
obvious hypothesis.

### 5.2 The mechanism: realizable error

| distance set handed to the fit | unrealizable part |
|---|---|
| the **real** predicted distances | **0.977 Å** |
| the retrieval pool's distances | **0.458 Å** |
| ORACLE i.i.d. Gaussian of the same RMS — the null | **2.369 Å** |
| ORACLE the same errors permuted across pairs | 2.081 Å |
| ORACLE the same magnitudes with random signs | 2.070 Å |
| **ratio real/null** | **0.413**, paired **−1.391 [−1.623, −1.184]** |

**The real prediction is 2.4× closer to being realizable by an actual conformation than
matched-magnitude noise.** Both structure-destroying surrogates move it back to within 13% of the
null. The negative correlation completes it: the part the fit cannot satisfy is the part
anti-correlated with truth, so **the optimiser discards the incoherent component and follows the
coherent one**.

Surrogate destruction on the real error isolates the expensive property *(n = 6, in flight)*:
randomising the error **signs** while keeping every magnitude is worth **−1.007 Å [−1.314, −0.664]**,
against 0.33 for matching a Gaussian, 0.26 for permuting across pairs, and 0.12 for clipping outliers.

### 5.3 The same phenomenon, opposite sign, on the torsion channel

Real torsion errors are **0.6–1.7 Å cheaper** than i.i.d. errors of the same magnitude across
**native-free** arms, and up to **2.9 Å cheaper for the ORACLE best-matching pool window**. Alignment is
**0.565 (ORACLE window)** and **0.665 (the incumbent's projected torsions)** against a **0.945**
random-direction null — but the control that matters is that a **zero-information constant α-helix
already scores 0.709**, so only those two arms clear it. Sign-flipping the ORACLE window's errors costs
**+1.596 Å [+1.385, +1.813]**; for the native-free retrieval mean the same operation is **null**
(−0.150 [−0.366, +0.070]). The 64°-RMS channel emitting **1.77 Å** where i.i.d. emits 4.71 Å is the ORACLE window. None of this is visible in pairwise correlations (|r| ≤ 0.16),
and synthetic coherence moves the result the *opposite* way — **correlated is not aligned**.

The Cα trace of an ideal-geometry chain responds to a participation ratio of only **3.47 effective
directions out of ~25.9**, with five exactly-null directions — the four inert torsions plus a fifth,
distributed, target-specific combination not previously recorded *(in flight)*.

### 5.4 Why every correction under-delivers

All of them act on the error's **magnitude profile**; the expensive part is a **direction**.
Recalibration is provably worth zero (a weighted least-squares argmin is invariant under uniform
rescaling of the weights, and the standardised residual's spread is near-constant across separation
at 2.19–2.98). Separation debiasing is worth little — negative on ranking, at the noise floor on the fit. Outlier clipping is worth
0.12 Å. Per-target rescaling has a 0.1 Å ORACLE ceiling. Consensus rescaling has a 0.275 Å ORACLE
ceiling that the distogram overshoots because its own bias corrupts the estimate.

---

## 6. ~~A NATIVE-FREE LAW FOR FUSION~~ — THE AMBIGUITY DECOMPOSITION, AND WHAT IT IS NOT

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT workstream. This section may not be drafted as a
> contribution.** Original text preserved below, struck through.
>
> | | |
> |---|---|
> | the expression | the **two-member Krogh–Vedelsby (1995) ambiguity decomposition** — ensemble error = mean member error − diversity. 31 years old, standard ensemble learning. **Cite it.** |
> | the mean | Sprint 15 used the **arithmetic** mean of `r₁, r₂`; the identity needs the **quadratic** mean. Channels differ by mean **0.980 Å** per target. Corrected residual **+0.075 Å**, not +0.162; paired **−0.087 [−0.119, −0.061]** i.i.d.-target, **[−0.129, −0.056]** fold-clustered, **W/L 126/0** |
> | the "prediction error" | **not a prediction error.** In one common frame with the correct mean the identity is exact to **2.65e−15** on all 126 real targets. The +0.075 Å decomposes exactly into **+0.174** (medoid-frame averaging) **−0.099** (Kabsch-minimised `s`) |
> | "`s` is native-free, so the fusion decision is a calculation on unlabelled data" | **REFUTED.** The prediction needs `r`, a distance to the native. `s` alone ranks "does fusion win here" at **AUC 0.401** against a 0.500 null, and a supervised leave-fold-out threshold degenerates to "always fuse" on 5/5 folds |
>
> **What may be drafted:** (a) the *measurement* that per-pair error correlation +0.602 coexists with
> 3.138 Å structural disagreement; (b) the `s/r` scaling statement, which is correct algebra;
> (c) the negative result in the last row, which is a genuine finding about consensus screens.
> `s16/retract_law.py`, `s16/retract_exact.py`, `s16/retract_screen.py`, `s16/retract_FINDINGS.md`.

Two channels' per-pair errors correlate at **+0.602**, and the structures they imply sit **3.138 Å
apart** — nearly as far from each other as either is from the native. Error correlation is the wrong
statistic; what matters is whether the realizable components point at the same wrong structure.

> ~~**d_avg ≈ √(r² − (s/2)²)**~~ → **d_avg = √( (r₁²+r₂²)/2 − (s/2)² )**, exact

| | ~~published~~ | **corrected** |
|---|---|---|
| `r` mean single-channel error | ~~3.641 Å~~ (arithmetic) | quadratic mean, per target |
| `s` structural disagreement (**native-free**) | 3.138 Å (Kabsch-minimised) | 3.506 Å (common frame) |
| predicted | ~~**3.205 Å**~~ | **3.292 Å** |
| measured | **3.367 Å** | 3.367 Å |
| ~~error~~ residual | ~~**+0.162 Å**~~ | **+0.075 Å**, and **2.65e−15** in a common frame |

~~No free parameters, and `s` needs no labels — so the fusion decision becomes a calculation on
unlabelled data.~~ **REFUTED — see the correction box.** The gain is small whenever `s ≪ r` and large
only as `s → 2r`; here `s/r = 0.88`, worth ~0.2 Å.

---

## 7. WHAT THE ARCHITECTURE ACHIEVES, AND WHAT IT DOES NOT

| stage | mean | vs baseline |
|---|---|---|
| **G** best in ensemble (**ORACLE**) | 2.760 | −0.444 [−0.573, −0.325] |
| **S** objective's argmin | 3.511 | +0.307 [+0.178, +0.439] |
| **A** coordinate consensus | 3.151 | −0.053 [−0.125, +0.020] |
| **F** after projection — like-for-like | **3.321** | **+0.117 [+0.050, +0.188]** |

**The method loses.** But the ensemble *contains* better structures than the baseline emits (G, −0.444
with an interval nowhere near zero; 49% of targets hold a member under 2.5 Å), so the entire deficit
is in the readout: selection costs +0.751, aggregation recovers −0.360, projection costs +0.170. And
pre-filtering the ensemble by the objective makes aggregation **worse**, closing the obvious response.

**A positive control on the objective.** Without any search, the restraint objective's argmin over a
500-member pool beats a random member by **−0.949 Å [−1.147, −0.753]** on 79% of targets, with
ρ(objective, RMSD) = +0.562 positive on 88.1%. At the same chain lengths where a contact potential's
minimum is *worse* than random, a distance-restraint objective's minimum is about an ångström
*better* — so that pathology is a property of contact potentials, not of short-chain prediction. The
native itself nonetheless sits at the **34.7th percentile** and is the argmin on 4 of 126 targets.

---

## 8. RELATION TO PRIOR WORK, INCLUDING WHAT WE CEDE

- **Certified-optimum-is-worse-than-random is published** for lattice contact potentials (Roget et
  al.). Our version survives only in its all-atom, continuous-torsion form and cites theirs as prior
  art for the phenomenon.
- **In-loop CVaR for peptide folding is published** (QuPepFold). We claim only the estimator defects
  and the diversity measurement.
- **Tail-restricted discrimination has a framework** (BEDROC, Truchon & Bayly). We cite it and
  position our metric as an application.
- **Corroboration:** AlphaFold2's rank-0 model is the lowest-RMSD one only 13% of the time against a
  20% null — state-of-the-art ranking performing worse than chance inside its own pool.
- **No predictive comparison is claimed.** The nearest comparator reports 4.89 Å on 75 fragments, a
  different target set; our 3.204 Å is **not** a like-for-like win.

---

## 9. LIMITATIONS

No experimental validation. Two retrieval constants were fitted on the reported targets (one is
load-bearing, one is not; sensitivities published). Targets share a library, folds and a predictor,
so intervals are fold-aware but dependence is not fully addressed. Quantum results are exact
simulation of a 262,144-configuration space at n = 9, k = 4. Several mechanism measurements are
currently n = 6–24, labelled inline, with full-instrument runs in flight. **An empirically measured
false-positive floor of ~0.08 Å applies to single-draw paired comparisons in this machinery**, and
every claim at or below it is marked as unresolved rather than reported as an effect. **No family-wise error
control** across the programme's arms. The defence was a 60-target benchmark held out and unread; the
protocol has now been frozen as **the baseline unchanged**, because none of thirteen candidates
cleared its control, and by a rule fixed in advance the benchmark is therefore **left unspent** rather
than used to re-measure an unmodified system. Everything reported here is **exploratory**, and that is
the honest limit of all of it.

---

## 10. CONCLUSION

The obstacle in this regime is not the conformational representation, not the search, and not the raw
accuracy of the distance channel. It is that the channel's errors are **mutually consistent**: they
describe a plausible wrong structure, and a faithful optimiser builds it. That reframes the
engineering target from *shrinking* a predictor's error to *decorrelating* it — a direction with, as
far as we can tell, no prior attempts — and it explains why a decade of magnitude-directed
corrections, and every one we tried, buy tenths of an ångström.

**We can also say how much decorrelation would be needed, and it is a lot.** At constant error
magnitude, destroying the coherence is worth 1.284 Å, and reaching 2.5 Å requires removing about
**69%** of it. The response is threshold-like rather than linear: partial decorrelation buys nothing
until most of the coherence is gone. That is a harder target than the first small-sample read
suggested, and it is stated at the n it was measured on.

**And we can say why the obvious way to act on it does not work.** The mechanism is real and enormous
— an oracle rotation of the fit's own error out of the structure's responsive subspace, at fixed
magnitude, clears 2.0 Å from the same restraints. But **no native-free intervention we could build
captures more than 4% of it**, across 22 arms, and pooled over 2,142 (target, arm) pairs the partial
correlation of ΔRMSD with Δalignment is +0.087 against β = +0.586 for plain error magnitude. Among
achievable interventions, **alignment is a description of what good channels happen to have, not a
lever**. So 2.0 Å is not information-limited here; it is **control-limited**, which is a different and
more tractable-sounding problem that we nonetheless could not solve.

The quantum result resolves the same way. A variational optimiser concentrates its distribution, and
concentration only helps when the objective can order what it concentrates onto. Here it cannot, so
concentrating is at best neutral — and whether it looks harmful or helpful depends on the budget, the
readout and the number of targets, none of which the field routinely reports. **The one apparent
positive we found did not survive its own replication**: the unit of analysis overstated n threefold,
the surviving cell failed its concentration check, and the control had been charged 2,048 objective
evaluations against the arm's 819,200.

**What we would hand to the next attempt.** A native-free error-direction surrogate exists —
`θ_fit − θ_pool` has |cos| 0.390 with the true error against a 0.157 null — and the structure's quiet
subspace is separately identifiable without the native (subspace overlap cos² 0.827 against a 0.499
null). Both ingredients an alignment intervention needs are estimable; putting them together as a
steering signal is untested, and it is the cleanest thing this work leaves behind.

---

## FIGURES

Rendered by `python -m s15.figures`; provenance and caveats in `s15/figures/README.md`. ORACLE arms
are drawn hatched and labelled; the incumbent is a reference rule on every accuracy panel; an
`n ≥ 100` guard refuses smoke reads, with three declared exceptions that each print their own n.

| | panel | where it appears |
|---|---|---|
| **1** | the generative cascade's four stages against the baseline | §7 |
| **2** | **error shape beats error magnitude** — three error models on one comparable axis | §5.1 |
| **3** | the constrained family's feasible set excludes the native at every usable ε | §7 |
| **4** | where each objective ranks the truth, and whether it orders the pool at all | §7 |
| **5** | the predictor's systematic over-prediction, growing with sequence separation | §5.4 |
| **6** | perfect distance knowledge is worth 0.611 Å, not ≈1.95 Å | §4 |
| **9** | undoing the averaging contraction — every native-free scale reference fails | §5.4 |
| **11** | **the errors describe a consistent wrong structure** | §5.2 |
| **12** | ~~a parameter-free law predicting fusion gain from native-free disagreement~~ **SUPERSEDED 2026-09-06 (Sprint 16 RETRACT): the two-member Krogh-Vedelsby ambiguity decomposition, exact; the plotted residual is a frame convention and the wrong mean** | §6 |
| **14** | no available channel reaches the accuracy 2.0 Å requires | §5.3 |
| **15** | concentration buys the set mean and pays the set best; α is a diversity dial | §3.2b |
| **16** | where the ångströms go — selection and projection lose ground, aggregation recovers part |§7 |

*(7, 8, 10 and 13 are held back: their runs were stopped as superseded or are still in flight, and the
figure module refuses to render them from smoke reads.)*
