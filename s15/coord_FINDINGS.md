# SPRINT 15 — coordinator findings

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **LITERATURE-SUPPORTED** / **HYPOTHESIS** /
**REFUTED**.

Instrument verified at sprint start: `python -m s12.instrument` reproduces
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`.

---

## HOW TO READ THIS FILE

Sections appear in the order they were **written**, not the order that makes them intelligible,
because the sequence of belief is part of the record — three findings here were retracted by their
author and one was corrected within the hour. Nothing has been reordered or edited away. This map
gives the logical reading order.

**The thesis and the frame**
- PROVISIONAL THESIS — written before the major experiments, timestamped so it cannot be retrofitted
- THE FOUR ARCHITECTURE FAMILIES — declared before implementation

**What is not the barrier** (read these first; they license everything after)
- **K0** the analytic torsion gradient is exact; four torsions per chain are inert
- **K1** → **K1-CORRECTED** the representation is not the barrier — ORACLE 0.611 Å — but the full
  instrument REVERSES the predicted arm's sign against the 8-target read
- **K5** the project's long-standing ~1.95 Å ceiling was a property of the LIBRARY, not the channel

**What the barrier is**
- **K2** the distogram is biased and 2.633× over-confident
- **K12** the accuracy requirement, and the controlled demonstration that error SHAPE beats MAGNITUDE
- **K9** the mechanism: the errors are geometrically REALIZABLE — a coherent wrong structure
- **K10** a native-free law for what channel fusion is worth
- **K3** the retrieval pool as a second distance channel
- **I1** the information-channel audit — the same phenomenon on the torsion channel, opposite sign

**What was tried and what happened**
- **K6** the objective's argmin beats random — the rebuttal to Roget et al.
- **K7** Family B falsified as declared, redirected to robust losses
- **K8** **the central negative: the full cascade does not beat the incumbent**
- **K11** the projection gap is worth 0.275 Å and the distogram cannot close it

**The quantum component**
- **Q1** CVaR-VQE loses to its own untrained initialisation on every diversity-respecting readout
- **Q2** the variational geometry, corrected

**Retractions and corrections — read these to calibrate how much to trust the rest**
- **K4** → **K4-RETRACTED** a "quantisation ceiling" that was my own indexing bug
- **K11** withdraws a 25.8% figure I reused in four places without checking (it is 3.5%)
- **A2** the multi-start draw was not reproducible across processes
- **A1** the Phase 0 audit — one retraction, one correction to the brief, one new defect

**Positioning**
- **L1** the novelty boundary — two claims taken by the literature, one that survives and leads

---

## PROVISIONAL THESIS, written before the major experiments

Required by the brief's paper-first principle, and to be updated or discarded the moment
evidence invalidates it. Recorded here with a timestamp so it cannot be retrofitted.

> **Peptide structure prediction in this regime is limited by the ACCURACY OF PREDICTED
> DISTANCE RESTRAINTS, not by conformational search, not by the conformational
> representation, and not by the physical energy model.**
>
> Supporting sub-claims, each independently falsifiable:
>
> 1. Torsion-space distance geometry recovers a structure from its true distance matrix
>    essentially exactly, so the representation is not the barrier. *(measured; see K1)*
> 2. The predicted restraints carry large, systematic, separation-dependent defects that are
>    **removable** rather than information limits. *(measured; see K2)*
> 3. Framing the problem as **solving** for a restraint-consistent conformation, rather than
>    generating-and-selecting candidates, bypasses the discrimination bottleneck that
>    defeated every Sprint 14 arm.
> 4. Within that frame, VQE/CVaR has a defensible role in **conditional ensemble
>    generation** — producing the bundle of restraint-consistent conformers that NMR
>    structure determination reports — rather than in single-structure energy minimisation,
>    where Sprint 14 showed it is worse than doing nothing.

**What would falsify it.** If the predicted-restraint fit cannot be brought below the
incumbent's 3.204 Å once its measured defects are corrected, sub-claim 2 is false and the
restraints are information-limited rather than defect-limited. If the corrected fit improves
but an equally-informed classical sampler matches the quantum ensemble at matched budget,
sub-claim 4 is false and the thesis reduces to a classical distance-geometry result.

---

## K0. The analytic torsion gradient is exact, and all four terminal torsions are inert

**DEMONSTRATED.** `core.project.frames` + `_torsion_grad` give an exact O(n) reverse-mode
gradient of any CA-function with respect to `(phi, psi)`. Checked against central differences
on a real target: **cosine 1.000000000000000**, and the only entries disagreeing are those
that are numerically zero.

Sprint 14's corrected claim that **four** torsions per chain are inert for the CA trace is
confirmed by two independent routes. Direct perturbation of ±0.7 rad moves the Kabsch RMSD by
**0 to 1.4e-07 Å** for `phi[0]`, `psi[0]`, `phi[n-1]`, `psi[n-1]` across six targets. And the
analytic gradient gives `|g|` of exactly 0, **2.98e-14**, 0 and 0 at those four entries
against a maximum of 3.31e+02 and an interior median of 6.76e+01 — a ratio of 9e-17.

This matters for qubit accounting in any paper: the live count is `(n-2) * log2(k)`, not
`n * log2(k)`.

## K1. Torsion-space distance geometry recovers a structure from its true distances essentially exactly

**ORACLE DIAGNOSTIC** for the ceiling arm; the predicted arm is **DEMONSTRATED** and
native-free. `s15/distgeo.py`, first read on 8 targets (full 126 running).

The method is new to this project and is a **generation** method rather than a selection one:
fit continuous `(phi, psi)` by L-BFGS to minimise
`sum_{|i-j|>=2} w_ij (d_ij(phi,psi) - dhat_ij)^2`, using the exact analytic gradient, from six
native-free starts, **selecting among starts by the objective only and never by RMSD**.

| arm | mean | median | <2 Å |
|---|---|---|---|
| **ORACLE, fit to the TRUE distance matrix** | **0.697** | **0.077** | **0.88** |
| predicted distances, inverse-variance weighted | 2.792 | 2.810 | 0.12 |
| predicted distances, unweighted | 3.102 | 3.129 | 0.00 |
| the retrieval circular-mean start it was given | 3.247 | — | — |

Three things follow, and the first is the most important the sprint has produced so far.

**The representation is emphatically not the barrier.** A median of **0.077 Å** means that on
half these targets, ideal-geometry torsions reproduce the native CA trace from its own
distance matrix to within a rounding error. This is a strictly stronger statement than the
Sprint 13 k=4 descent ceiling of 1.982 Å, because it is continuous, solved rather than
searched, and does not depend on a library.

**The per-pair uncertainty is worth 0.31 Å**, the first time this project has used the
distogram's `sd` for anything at all. It has been sitting unused in `I.distogram` for
fourteen sprints.

**And the fit improves on its own start by 0.455 Å** (3.247 → 2.792), so the optimisation is
doing real work rather than returning what it was handed.

Caveat carried forward: on these 8 (easy) targets the incumbent is 2.277, so the predicted arm
**loses** by +0.515 [-0.049, +1.128]. The full-instrument number is what counts.

## K2. The distogram is systematically BIASED and its uncertainties are 2.6x over-confident

**DEMONSTRATED** (post-hoc evaluation of a leave-fold-out predictor; not leakage).
All 126 targets, 8,549 pairs.

| separation | n pairs | MAE | **bias** | RMSE | mean sd | **z-sd** |
|---|---|---|---|---|---|---|
| 2-3 | 2636 | 0.974 | **+0.113** | 1.463 | — | 2.186 |
| 4-5 | 2132 | 2.148 | **+0.323** | 2.963 | — | 2.983 |
| 6-7 | 1628 | 2.818 | **+0.615** | 4.001 | — | 2.699 |
| 8-10 | 1506 | 3.745 | **+0.931** | 5.216 | — | 2.701 |
| 11-15 | 647 | 4.666 | **+1.492** | 6.328 | — | 2.671 |
| **global** | 8549 | 2.386 | **+0.509** | 3.704 | — | **2.633** |

**Defect 1 — the `sd` is 2.633x over-confident.** A calibrated predictor gives a standardised
residual with sd 1.0. This is a genuine reporting defect that must be declared in any paper
using this distogram.

**But it barely matters for the fit, and that is worth knowing before anyone spends a sprint
on it.** The z-sd is nearly *constant* across separation (2.19, 2.98, 2.70, 2.70, 2.67), so
the `sd` captures the error's **shape** correctly and is wrong only by a near-uniform factor —
and a weighted least-squares argmin is **invariant under uniform rescaling of all weights**.
Recalibration alone cannot change a single structure.

**Defect 2 — a systematic positive bias that grows monotonically with separation**, from
+0.113 Å at separation 2-3 to **+1.492 Å at separation 11-15**. The distogram systematically
**over-predicts** distance, so a least-squares fit to it produces a systematically
**over-extended** structure. This is a pure, removable, additive error and **nobody has ever
corrected it.**

It is also neatly opposed to a Sprint 14 result: coordinate averaging **contracts** the
backbone. (The 25.8% figure quoted here was reused without checking and is WITHDRAWN -- see
K11; measured directly, the production consensus is contracted by 3.5% against the true
distances. The direction stands, the magnitude does not.) The two stages of the incumbent pipeline have **opposite** geometric
pathologies, which is a striking observation in its own right and may explain part of why the
incumbent works as well as it does.

**Correcting the bias is legitimate and native-free at inference**, because the bias is a
property of the *predictor*, estimable from training folds alone. `s15/distcal.py` fits
`b(sep)` on four folds and applies it to the held-out fifth, with an in-fold ORACLE arm as the
ceiling — the control that separates a transfer failure from a functional-form failure, a
distinction Sprint 14 showed is decisive.

---

## THE FOUR ARCHITECTURE FAMILIES

Declared before implementation, per the brief. Two architectures count as distinct only if an
experiment exists whose predicted outcome differs between them.

**A — DISTANCE-RESTRAINED TORSION VQE.** *Hypothesis:* predicted distance restraints define a
many-body diagonal Hamiltonian over discrete torsion states, and VQE/CVaR can sample
restraint-consistent conformations. *Why it differs:* the objective is **generative**
(restraint satisfaction) rather than **selective** (ranking candidates), so the tail-ordering
pathology that defeats every Sprint 14 search arm does not arise. *Distinguishing experiment:*
does the quantum arm beat the continuous L-BFGS fit on ensemble diversity at matched restraint
satisfaction? *Falsified if:* the discrete arm cannot reach the continuous fit's restraint
residual at any budget.

**B — CONSTRAINED / AUGMENTED-LAGRANGIAN.** *Hypothesis:* restraints should **constrain** the
physical manifold, not compete with the energy in a weighted sum — `min E_AMBER` subject to
`C_restraint <= epsilon`. *Why it differs:* this is literally restrained refinement as
practised in NMR structure determination, it gives AMBER a legitimate first-class role, and
its feasible-set analysis directly answers whether the restraint set even *contains*
sub-2.5 Å structures. *Distinguishing experiment:* the epsilon phase diagram — a weighted sum
has no analogue of the feasibility boundary. *Falsified if:* the feasible set is empty, or
contains no structure better than the unconstrained fit.

**C — CONDITIONAL ENSEMBLE GENERATION.** *Hypothesis:* the right output is `P(x | S, O)`, a
bundle of restraint-consistent conformers, not a single argmin — which is what NMR
determination actually reports. *Why it differs:* it is evaluated on ensemble statistics
(coverage, diversity, calibration) that A and B do not produce. *Distinguishing experiment:*
does the ensemble's coordinate consensus beat its own best single member? Sprint 14 showed
aggregation is worth 3.4x the objective, so this is where that leverage lives.
*Falsified if:* consensus over the ensemble is no better than the single fit.

**D — RESERVE.** Unallocated by design. To be filled only if Phase 1 shows A-C cannot attack
the measured bottleneck. Not filled speculatively.

**The four mandatory components have measured, legitimate roles in all three**: VQE samples or
optimises the torsion configuration; CVaR shapes the ensemble tail (Family C is where it can
finally earn its place); **Legacy is a clash gate**, which Sprint 14 measured as its only
honest role; and **AMBER is a validator and refiner** in A and C and the **objective** in B.

---

## K1-CORRECTED. The full instrument REVERSES the sign of the predicted arm

**DEMONSTRATED, and it is a correction to K1 above, which was written on 8 targets.**
`s15/results/distgeo.json`, all 126 targets. The 8-target read is left standing above,
unedited, so the size of the small-sample error stays visible.

| arm | mean | median | sd | <2 Å | FAIL18 | vs incumbent (paired) |
|---|---|---|---|---|---|---|
| **ORACLE, fit to the TRUE distance matrix** | **0.611** | **0.055** | 1.143 | **0.86** | 1.156 | **−2.593 [−2.901, −2.295]** |
| predicted distances, inverse-variance weighted | 3.644 | 3.466 | 1.580 | 0.16 | 5.972 | **+0.440 [+0.290, +0.592]** |
| predicted distances, unweighted | 3.798 | 3.625 | 1.530 | 0.12 | 6.120 | +0.594 [+0.430, +0.761] |
| the retrieval circular-mean start it was handed | 4.072 | — | — | — | — | +0.868 |

**What survives from K1, and gets stronger.** On the full 126 the oracle fit is **0.611 Å mean,
0.055 Å median, 86% under 2 Å**, beating the incumbent by −2.593 with an interval nowhere near
zero. Torsion-space distance geometry recovers a native CA trace from its own distance matrix
essentially exactly. The conformational representation, the ideal-geometry constraint and the
search are **definitively not the barrier**. This is the cleanest ceiling the project has
measured.

The `sd` is still worth real accuracy (3.644 vs 3.798 unweighted, +0.154 Å here against the
0.31 Å seen on 8 targets), and the fit still improves on the start it was handed by **0.428 Å**
(4.072 → 3.644), so the optimiser is doing real work rather than returning its input.

**What is REFUTED.** The 8-target table put the predicted arm at 2.792 and said it loses by
+0.515 "on easy targets — the full-instrument number is what counts." The full-instrument
number is **+0.440 [+0.290, +0.592]** at a mean of **3.644**, not 2.792. So:

> **Naive least-squares distance geometry from the raw predicted distogram is WORSE than the
> incumbent pipeline, on the full instrument, with a confidence interval excluding zero.**

That is a negative result for the simplest form of the thesis and it is recorded as one. The
provisional thesis is *not* yet falsified — its stated falsification condition was "if the
predicted-restraint fit cannot be brought below 3.204 Å **once its measured defects are
corrected**", and the raw arm corrects nothing. But the burden has moved: bias removal (K2), the
pool channel (K3) and aggregation must now find **0.44 Å before they find anything at all**, and
each has to be measured rather than assumed.

**Why this is the more informative failure.** The 3.03 Å between the oracle fit (0.611) and the
predicted fit (3.644) is now measured to be **entirely restraint error** — both arms share the
identical optimiser, starts, selection rule and geometry, so nothing else can be carrying it. No
earlier sprint had a decomposition this clean. Sprint 13 could only say a k=4 torsion descent
bottoms out at 1.982 Å; this says the continuous solve bottoms out at 0.611 Å and every
remaining ångström belongs to the distogram.

## K3. The retrieval pool is a genuinely independent second distance channel

**DEMONSTRATED** (error profile ORACLE-scored post hoc; the channel itself is native-free).

Nobody in fourteen sprints has read the **distances** out of the K=500 retrieval pool, though
its coordinates have been in the universe cache the whole time. Sprint 14 used the pool for
torsions and found its untrained prior beats the trained leave-fold-out sequence model. The
distance analogue:

| | pool median | distogram |
|---|---|---|
| global MAE | **2.249** | 2.386 |
| global bias | −0.128 | **+0.509** |

The pool is slightly *more* accurate with **no training at all**, and — the part that matters —
its bias has the **opposite sign at short and medium separation**. Two channels whose errors
point in opposite directions is the precondition for fusion to do anything, and it is the first
time this project has had one.

The caution travelling with it is Sprint 14's fusion arithmetic: gain goes as the **square** of
the weaker channel's skill, which is why Legacy+AMBER fusion bought +0.004 to +0.011 Å despite a
genuine +0.096 truth-partialled decorrelation. The `combined` arm must therefore be *measured*,
and a fusion that buys nothing is the expected outcome, not a surprise.

---

## K12. THE ACCURACY REQUIREMENT, and the clean demonstration that SHAPE beats MAGNITUDE

**DEMONSTRATED.** `s15/distacc.py`, 40 targets, three realistic error models at eight severities each.
This converts the whole programme from an aspiration into a stated engineering requirement, and it
supplies the cleanest possible demonstration of K9's mechanism — because here the error shape is
*controlled* rather than observed.

### The requirement

True distances are corrupted under a named error model at a controlled severity and refitted through
the identical machinery. The three models inject different *shapes* at different magnitudes, so the
severities are put on one comparable axis — the **effective RMS of the injected error** (the
`sep_scaled` model's exact moment factor is `√E[sep²]/E[sep] = 1.1331`, the `outliers` model's is
`√(0.09 + 0.05·36) = 1.375`).

| error model | required effective RMS for 3.0 Å | for 2.5 Å | for 2.0 Å |
|---|---|---|---|
| `abs_gauss` — i.i.d. Gaussian | > 3.00 | **2.87** | **1.87** |
| `sep_scaled` — error growing with separation | > 3.40 | **2.95** | **1.85** |
| `outliers` — 5% of pairs badly wrong, rest nearly right | **> 4.12** | **> 4.12** | **> 4.12** |

**The distogram's actual RMSE is 3.704 Å, and it emits 3.644 Å.**

> **To reach 2.5 Å the distance channel would have to improve from 3.70 Å RMSE to about 2.87 Å — a
> 22% reduction — and that is the requirement for an i.i.d. channel. Ours is not i.i.d.**

For 2.0 Å the requirement is an RMSE of **1.87 Å**, which is roughly **halving** the current error.

### The clean result: separation-scaling costs nothing, outlier-shaping is nearly free

**`sep_scaled` and `abs_gauss` need the same accuracy** — 2.95 against 2.87 for 2.5 Å, and 1.85
against 1.87 for 2.0 Å. Once magnitude is matched, making the error grow with sequence separation
costs essentially nothing. That is a clean null, and it retires a natural hypothesis: the distogram's
separation-dependent error profile (MAE 0.974 Å at separation 2–3 rising to 4.666 Å at 11–15) is
**not** what makes it expensive.

**The `outliers` row is the striking one.** At an effective RMS of **4.12 Å** — larger than the
distogram's own 3.70 Å — an outlier-shaped channel still emits **1.993 Å**, better than the 2.5 Å
target and far better than the 2.561 Å that i.i.d. noise of only 3.00 Å RMS produces.

| channel | effective error RMS | emitted RMSD |
|---|---|---|
| outlier-shaped | **4.12 Å** | **1.993 Å** |
| i.i.d. Gaussian | 3.00 Å | 2.561 Å |
| **the real distogram** | **3.70 Å** | **3.644 Å** |

> **A larger error concentrated in a few pairs is far cheaper than a smaller error spread evenly, and
> both are far cheaper than the real distogram's error at comparable magnitude. Error shape dominates
> error magnitude, in both directions, by more than an ångström.**

### Why this settles the mechanism rather than merely illustrating it

Every earlier statement of K9 rested on *observing* the error's structure and inferring what it cost.
Here the structure is **imposed**, so the causal direction is not in question, and the three rows
bracket the real channel from both sides:

- the fit tolerates a **few catastrophic** restraints easily — it can ignore them, which is why
  outlier-shaped error at 4.12 Å RMS is nearly free, and why winsorising the real error is worth only
  0.12 Å;
- the fit tolerates **many mildly wrong** restraints less well — i.i.d. noise at 3.00 Å RMS already
  costs 2.561 Å, because every restraint pulls a little and none can be ignored;
- and the real distogram, at 3.70 Å RMS, is **worse than either** at 3.644 Å, because its errors are
  not merely spread out but **mutually consistent** — they describe a coherent wrong structure that
  the fit reproduces faithfully (K9: 56% realizable, 2.3× closer to realizable than matched noise).

The ordering `outliers ≪ i.i.d. < real` at comparable magnitude is the whole argument in one line.

### What this means for the sprint's targets

**2.5 Å is not reachable by any means measured in this sprint.** It requires a 22% reduction in
distance RMSE even for an i.i.d. channel, and our channel is worse than i.i.d. at equal magnitude by
about 1.1 Å. No intervention tested — recalibration (provably zero), separation debiasing (worth
little; negative on ranking, at the noise floor on the fit), robust losses (**refuted at n = 126**), per-target rescaling
(0.1 Å ceiling), consensus rescaling (0.275 Å ceiling, unreachable native-free), or channel fusion
(0.17–0.22 Å, and priced in advance by K10's law) — comes close to closing 0.8 Å.

**And the requirement is now specific enough to be actionable.** It is not "make the distogram
better." It is: **reduce its RMSE by 22%, or make its errors incoherent.** The second is the
interesting half, because K12 shows an error that is *larger* but *unstructured in the right way* is
worth more than a smaller one — and nothing in this project has ever tried to shape a predictor's
error rather than shrink it.

---

## K9/K10-CONFIRMED. The full instrument: realizability strengthens, the fusion law holds, and one of my controls was mislabelled

**DEMONSTRATED, all 126 targets.** `s15/coherence.py`. This is the run that answers the sharpest
criticism the adversarial reviews produced — that the sprint's interpretive core rested on n = 6, 8
and 24 in a project whose history includes an 8-target read that reversed a sign.

### Realizability — CONFIRMED, and stronger than the smoke read

| distance set handed to the fit | unrealizable part |
|---|---|
| **the REAL predicted distances** | **0.977 Å** |
| the retrieval pool's distances | **0.458 Å** |
| **ORACLE i.i.d. Gaussian, matched RMS — the null** | **2.369 Å** |
| ORACLE errors permuted across pairs | 2.081 Å |
| ORACLE magnitudes kept, signs randomised | 2.070 Å |
| *(mean \|error\| for reference)* | *2.339 Å* |

**ratio real/null = 0.413**, paired **−1.391 [−1.623, −1.184]**.

The effect is **larger** at n = 126 than at n = 24 (ratio 0.413 against 0.427), with an interval
nowhere near zero. The real prediction is **2.4× closer to realizable** than matched-magnitude noise,
and both structure-destroying surrogates land within 13% of the i.i.d. null — so the property being
destroyed is exactly the property that made it realizable.

**A new observation: the retrieval pool is even more coherent than the distogram** (0.458 against
0.977, against a 2.339 Å mean error). Its distances describe an *even better* wrong structure. That
sharpens K6's result — the pool is a better estimator and a worse objective — into something more
specific: the pool is confidently, self-consistently wrong.

### The fusion law — holds, with its prediction error growing slightly

> ### ⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT workstream. THIS SUBSECTION IS SUPERSEDED.
>
> **The law is prior art.** It is the two-member case of the **Krogh–Vedelsby (1995) ambiguity
> decomposition** (`s16/lit_FINDINGS.md` §B.5). Cite it; do not claim it.
>
> **It was computed with the wrong mean.** `s15/coherence.py` used the ARITHMETIC mean of the two
> channels' RMSDs; the identity requires the **QUADRATIC** mean `√((r₁²+r₂²)/2)`. The two channels
> differ by mean **0.980 Å** per target (median 0.567, max 4.798), so the substitution is not small.
>
> **The corrected numbers** (`s16/retract_law.py`, n = 126, no refit — recomputed from the
> per-target rows of `coherence.json`):
>
> | | ~~arithmetic (published)~~ | **quadratic (correct)** |
> |---|---|---|
> | law prediction | ~~3.205~~ | **3.292** |
> | measured coordinate average | 3.367 | 3.367 |
> | **prediction error** | ~~**+0.162**~~ | **+0.075** (median **+0.019**) |
>
> Paired, TARGET as the unit: **−0.087 [−0.119, −0.061]** i.i.d.-target bootstrap,
> **[−0.129, −0.056]** fold-clustered, median −0.015, **W/L 126/0**. **54% of the published
> "prediction error" was an arithmetic artefact, not physics.**
>
> **And the remaining 0.075 Å is a frame convention, not physics either.**
> `s16/retract_exact.py` reproduces this table's `Xd`/`Xp` **bit-identically** (max |Δ| = 0.000e+00
> on all four quantities) and shows that with both structures superposed onto the native and the
> quadratic mean, the identity is **exact to |err| mean 2.65e−15, max 2.52e−14, sign balanced
> 55/61** — i.e. exactly floating-point noise. The surviving +0.075 Å decomposes exactly as
> **+0.174** (`I.coordinate_average` averages in the **medoid** frame, not the target's) **−0.099**
> (`s` is Kabsch-**minimised**, 3.138 Å, where the identity wants the same-frame disagreement,
> 3.506 Å). The final Kabsch of the average onto the native contributes **exactly zero**
> (mean = median = max = 0.000000 on 126/126).
>
> **Therefore "a parameter-free identity predicting an unseen quantity to 0.162 Å across 126
> targets" may no longer be quoted.** The identity predicts nothing unseen: it is an exact algebraic
> rearrangement, and this subsection measured a frame mismatch. Full account:
> `s16/retract_FINDINGS.md`.
>
> **What is NOT retracted, and is unaffected:** the *dissociation* below — per-pair errors correlate
> at +0.602 while the structures they imply sit 3.138 Å apart — is a measurement, not an identity,
> and it stands. So does the whole realizability half of this module.

**ORIGINAL TEXT, PRESERVED:**

| | n = 24 | **n = 126** |
|---|---|---|
| distogram fit | 3.270 | **3.622** |
| pool fit | 3.469 | **3.659** |
| structural disagreement `s` (**native-free**) | 2.980 | **3.138** (median 3.420) |
| raw per-pair error correlation | +0.627 | **+0.602** |
| ~~**law prediction**~~ | ~~2.970~~ | ~~**3.205**~~ → **3.292** |
| **measured coordinate average** | 3.105 | **3.367** |
| ~~**prediction error**~~ | ~~+0.135~~ | ~~**+0.162**~~ → **+0.075** |

~~A **parameter-free** identity predicting an unseen quantity to **0.162 Å** across 126 targets. It
cannot be rescued by refitting because there is nothing to fit.~~ The dissociation that motivates it
also holds: per-pair errors correlate at **+0.602** while the structures they imply sit **3.138 Å
apart** — correlated errors, different wrong structures.

### THE CORRECTION: my own control was an ORACLE, and it flips the sign

`s15/coherence.py` reported a field named `coordavg_vs_best_single`, computed as
`np.minimum(rmsd_disto, rmsd_pool)` per target. **That is a per-target ORACLE selection** — it
requires knowing which channel won on each target — and naming it "the best single channel" makes
fusion look worse than it is. On the full instrument the two comparisons have **opposite signs**:

| comparison | coordinate average | restraint-level fusion |
|---|---|---|
| **vs the better channel chosen globally — NATIVE-FREE** | **−0.255 [−0.370, −0.145]** | −0.071 [−0.156, +0.014] |
| vs a per-target ORACLE pick — *a ceiling, not a control* | +0.216 [+0.133, +0.302] | +0.400 [+0.262, +0.546] |

**The honest, native-free result is that coordinate averaging beats the better single channel by
−0.255 Å with an interval excluding zero.** Fusion works. Reported against the oracle control it
would have looked like a significant *loss* of +0.216, and I would have written that down.

The field is renamed and both comparisons are now emitted, labelled for what they are. This is the
fourth control in this sprint that inverts its own conclusion if mislabelled — after the shrinkage
confound in the decorrelation ladder, the projected-versus-unprojected comparison in the cascade, and
the budget convention in the quantum-natural-gradient work.

### A second sign flip from the small sample, worth recording

At n = 24, **restraint-level fusion beat coordinate averaging** (3.046 against 3.105) and I recorded
it as "a direction rather than a result, since n = 24". At n = 126 it is **reversed**: coordinate
averaging 3.367 beats restraint-level fusion 3.551, and only the coordinate average clears its
native-free control significantly (−0.255 [−0.370, −0.145] against −0.071 [−0.156, +0.014]).

Combining two channels at the **structure** level beats combining them at the **restraint** level.
That is consistent with the realizability picture: a single fit to blended restraints is forced to be
realizable throughout, so it must pick *one* coherent wrong structure; averaging two separately-fitted
structures is not so constrained and can land between two different wrong structures, which is
precisely where the law says the gain lives.

Caution retained on my earlier phrasing: I wrote at n = 24 that restraint-level fusion "supports the
cascade's use of the combined likelihood over a two-stage alternative." **That is now refuted** — the
full instrument says the two-stage alternative is better.

---

## K13. The dose–response curve for DECORRELATING the error, and an actionable target

**ORACLE DIAGNOSTIC, n = 8 smoke read; the full 126-target run is queued.** Recorded now because the
*design* — specifically its confound control — is the contribution, and because the number it
produces is the first genuinely actionable engineering target the sprint has generated.

### The gap this closes

K9 says the distogram's errors are coherent and that magnitude-directed corrections cannot touch
them. K12 says reaching 2.5 Å needs a **22% reduction in RMSE** for an i.i.d. channel. Together they
imply the lever is **decorrelation, not shrinkage** — but "decorrelate it" is not a number, and
nobody can act on it without knowing *how much*. A development effort that must halve the coherence
is a very different project from one that must eliminate it.

### The experiment, and the confound that would have destroyed it

Interpolate continuously between the real error and a sign-destroyed version of itself,
`e(t) = (1 − t)·e_real + t·e_destroyed`, and measure the emitted structure along the way.

**Mixing two vectors shrinks the result.** `e(t)` has smaller RMS than either endpoint whenever they
are not parallel, because the components partially cancel — so an uncontrolled interpolation
confounds decorrelation with plain error reduction, and plain error reduction is precisely the thing
K9 says is *not* the lever. Every arm is therefore run twice: **`rescaled`**, renormalised to the
real error's RMS, which is the experiment; and **`raw`**, which is the confound, reported beside it.

| ladder | t = 0 | t = 0.5 | t = 1.0 |
|---|---|---|---|
| **`rescaled` — constant error magnitude** | **2.749** | **2.337** | **2.012** |
| `raw` — magnitude allowed to shrink | 2.749 | **1.709** | 2.012 |

**The confound is large and it is not subtle.** At t = 0.5 the `raw` ladder reads 1.709 Å — better
than *either* endpoint, which is geometrically impossible as a decorrelation effect and is entirely
the shrinkage. Reporting the `raw` ladder alone would have produced a spectacular and completely
spurious result. This is the third time in this sprint that an obvious experiment needed a control to
avoid inverting its own conclusion.

### The result

**At constant error magnitude, destroying the coherence takes 2.749 → 2.012 Å — a gain of 0.737 Å
with no improvement whatsoever in the error's size.** The response is close to linear in `t`.

Inverting the ladder:

> **Reaching 2.5 Å requires destroying only about 30% of the error's coherence (t = 0.302), at zero
> improvement in error magnitude.**

### Why this reframes the engineering target

K12's requirement — a 22% RMSE reduction — is a demand for a materially better predictor, and this
project has been trying and failing to produce one for fifteen sprints. K13's requirement is a demand
for a predictor whose errors are **less mutually consistent**, which is a different objective
entirely, and one that:

- **is not what any predictor is trained for.** Distogram training minimises expected per-pair loss,
  which is a magnitude objective; nothing in it penalises errors that agree with one another. A
  training objective that penalised the *realizability* of the residual would be a genuinely new
  thing, and K9's measurement gives it a loss to minimise;
- **is cheaper than it sounds.** 30% of the coherence, not all of it;
- **is measurable natively.** The realizability ratio (`|d_fit − d̂|` against `|d̂ − d_true|`) needs the
  native, but its native-free half — the fit residual `|d_fit − d̂|` — is computable at inference and
  is a candidate regulariser and a candidate confidence signal in its own right.

### K13-FINAL, all 126 targets

**DEMONSTRATED.** Every arm except t = 0 is an **ORACLE DIAGNOSTIC**: interpolating away from the true
error requires knowing it, so this measures what decorrelation would be *worth*, not a method.

| ladder | t = 0 | 0.25 | 0.5 | 0.75 | t = 1 |
|---|---|---|---|---|---|
| **`rescaled` / signs** — constant magnitude, **the experiment** | **3.624** | 3.606 | 3.086 | 2.635 | **2.320** |
| `rescaled` / permute — constant magnitude | 3.624 | 3.540 | 3.325 | 3.102 | 2.818 |
| *`raw` / signs* — magnitude allowed to shrink, **the confound** | *3.624* | *3.147* | *2.580* | *2.436* | *2.320* |
| *`raw` / permute* | *3.624* | *3.232* | *2.940* | *2.897* | *2.818* |

**Destroying the sign coherence at constant error magnitude is worth 1.304 Å** (3.624 → 2.320). The
response is **threshold-like**: nothing happens by t = 0.25, then it falls steeply.

**And the requirement is far stiffer than either smaller read suggested:**

| target | fraction of the coherence that must be destroyed |
|---|---|
| the incumbent's 3.204 Å | **t = 0.443** |
| **2.5 Å** | **t = 0.857** |
| **2.0 Å** | **not reachable** — complete decorrelation stops at 2.320 Å |

> **Even PERFECT decorrelation of this predictor's error does not reach 2.0 Å.** At the current error
> magnitude the fully-decorrelated ceiling is 2.320 Å, so 2.0 Å requires *both* decorrelation **and**
> a smaller error — the two levers are not interchangeable and neither alone suffices.

**The sign structure matters more than the pair structure.** Destroying signs reaches 2.320 Å where
permuting across pairs reaches only 2.818 Å, on identical magnitudes — consistent with K9's
realizability account, in which what costs is a coherent *direction* rather than a scrambled
*assignment*.

**The confound is large and the control earns its place**: 0.506 Å at t = 0.5 (raw 2.580 against
rescaled 3.086). **And the raw ladder is monotone** (3.624 → 3.147 → 2.580 → 2.436 → 2.320),
confirming at full scale the withdrawal recorded below — the "dips below both endpoints, geometrically
impossible" reading was an eight-target artefact and does not exist.

**The requirement moved twice**: 30% at n = 8, 69% at n = 40, **86% at n = 126**, with 2.0 Å going
from "reachable" to "not reachable by this lever at all". A fourth instance of a small-n read moving
materially on the full instrument.

### K13-CORRECTED at n = 40, and one claim withdrawn

**The n = 8 read above is superseded.** At n = 40 the constant-magnitude ladder is
**3.411 / 3.415 / 2.875 / 2.368 / 2.127** across t = 0 / 0.25 / 0.5 / 0.75 / 1.

| | n = 8 | **n = 40** |
|---|---|---|
| gain from full decorrelation | 0.737 Å | **1.284 Å** |
| shape | "close to linear" | **flat to t = 0.25, then steep** — there is a threshold |
| coherence that must go to reach 2.5 Å | 30% | **69%** |

**And the confound story is WITHDRAWN.** I wrote that the uncontrolled ladder read 1.709 Å at
t = 0.5 — below *both* endpoints, "geometrically impossible as a decorrelation effect" — and used it
as a third example of an experiment inverting without its control. **At n = 40 the uncontrolled
ladder is monotone** (3.411 → 2.899 → 2.401 → 2.181 → 2.127) and never dips below its endpoint. The
confound is real and large (0.474 Å at t = 0.5) and the control is still the right design — but the
specific inversion I described was an artefact of eight targets, and the paragraph claiming it is
withdrawn.

That makes this the **third** n ≤ 8 read in the sprint to move materially on the full instrument,
after K1 (2.792 → 3.644, sign reversed) and the robust-loss "redescending signature" (refuted
outright at n = 126).

### Scope, stated plainly

n = 8, on the first eight targets, where the incumbent is 2.277 Å rather than 3.204 Å — an easier
slice. **Every arm except t = 0 is an ORACLE DIAGNOSTIC**: interpolating away from the true error
requires knowing it, so this measures what decorrelation would be *worth*, and is not a method. The
full-instrument run will move the magnitudes; the claim being made is the *shape* — a near-linear
response, with the 2.5 Å target reachable at partial decorrelation — and the demonstration that the
shrinkage confound must be controlled or the experiment inverts.

Two things must follow before this can be built on: the 126-target ladder, and the separate question
of whether a predictor **can** be trained to have incoherent errors, which this experiment is
designed to justify or forestall rather than answer.

---

## K9. THE DISTOGRAM'S ERRORS ARE GEOMETRICALLY REALIZABLE — they describe a consistent WRONG STRUCTURE

**DEMONSTRATED.** This is the sprint's unifying mechanism. It explains, from one measurement, why
every correction attempted in fifteen sprints has bought so little, and it makes a sharp prediction
about which interventions can possibly work.

### The chain of measurements that produced it

**Step 1 — the real error costs more than i.i.d. noise of the same size.** The phase diagram
corrupts *true* distances with i.i.d. Gaussian noise and refits:

| σ (Å) | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|
| emitted RMSD | 0.554 | 0.709 | 0.939 | 1.306 | 1.450 | 1.787 | 2.075 | **2.561** |

The distogram's MAE of 2.386 Å corresponds to a Gaussian σ ≈ 2.99, where the surface says **2.56 Å**.
The real distogram emits **3.644 Å**. The real error costs roughly **1.1 Å more than i.i.d. noise of
its own magnitude**.

**Step 2 — surrogate destruction says the expensive property is the SIGN pattern.** Taking the
measured error vector `e = d̂ − d_true`, destroying exactly one property, and refitting (n = 6;
every surrogate arm is an **ORACLE DIAGNOSTIC**, since destroying a property requires knowing it):

| surrogate | mean | vs real |
|---|---|---|
| `real` (native-free) | 2.832 | — |
| **`ORACLE_shuffle_signs`** — every magnitude kept, signs randomised | **1.825** | **−1.007 [−1.314, −0.664]** |
| `ORACLE_gauss_matched` — i.i.d. Gaussian at the target's own error RMS | 2.501 | −0.332 [−0.868, +0.195] |
| `ORACLE_shuffle_pairs` — permuted across pairs | 2.577 | −0.256 [−0.884, +0.431] |
| `ORACLE_winsor_z3` — outlier tail clipped | 2.711 | −0.121 [−0.322, +0.057] |
| `debias_sep` (native-free) — the existing leave-fold-out correction | 3.012 | **+0.179** |

Keeping every error magnitude exactly and randomising only the **signs** is worth about an ångström,
with an interval excluding zero at n = 6. Outliers are worth 0.12; separation structure 0.26;
matching an i.i.d. Gaussian 0.33. **The expensive property is that the errors agree in direction.**

**Step 3 — it is not a scale error, and the pool cannot estimate it.** The obvious reading of
"agree in direction" is that the structure is predicted at the wrong overall size. Tested with per
target rescaling, native-free from the retrieval pool, against ORACLE ceilings (n = 8):

| arm | mean | vs real |
|---|---|---|
| `ORACLE_scale` — the best possible isotropic rescale | 2.673 | **−0.090 [−0.284, +0.033]** |
| `ORACLE_affine` — the best possible affine map | 2.647 | −0.116 [−0.448, +0.194] |
| `pool_scale`, `pool_shift`, `pool_affine`, `pool_scale_rg` | 2.79–2.91 | **all +0.02 to +0.15, i.e. worse** |

**The ceiling of *any* per-target isotropic rescaling is about 0.1 Å**, not the 1.0 Å that sign
destruction buys. And the pool's estimate of the per-target scale correlates **+0.035** with the true
one — indistinguishable from nothing. So the coherent component is **not a scale**, and this
direction is closed by its own ceiling in a single run.

**Step 4 — the direct test, and the answer.** If the errors agree in direction because they describe
a *consistent alternative structure*, then the predicted distance matrix should be nearly
**realizable**: there should exist a real 3-D conformation that reproduces it closely. Measured on 24
targets by fitting to `d̂` and comparing the converged fit's own residual against the prediction's
true error:

| quantity | value |
|---|---|
| mean \|d_fit − d̂\| — the **unrealizable** part of the prediction | **0.879 Å** |
| mean \|d̂ − d_true\| — how **wrong** the prediction is | **2.009 Å** |
| ratio unrealizable / wrong | **0.437** |
| corr(fit residual, true error) | **−0.513** |

**Step 4b — the null control, which is what makes step 4 sound.** A residual of 0.879 Å is only
meaningful against what the *same fit* leaves on a distance set of the same magnitude with no
structure. The fit is constrained — `2n − 4` degrees of freedom against `n(n−1)/2 − (n−1)` pairs — so
it cannot match an arbitrary distance matrix, and some residual is guaranteed. Measured on the same
24 targets:

| distance set handed to the fit | unrealizable part |
|---|---|
| **the REAL predicted distances** | **0.879 Å** |
| i.i.d. Gaussian of the same RMS — **the null** | **2.058 Å** |
| the same errors permuted across pairs | 1.824 Å |
| the same magnitudes with random signs | 1.790 Å |
| *(mean \|error\| for reference)* | *2.009 Å* |

**The real prediction is 2.3× closer to realizable than matched-magnitude noise** (ratio 0.427), and
both surrogates that destroy its structure move it back to within 13% of the i.i.d. null. So the low
residual is a property of the prediction, not of the fit's flexibility, and the two surrogate arms
independently confirm that the structure being destroyed is exactly the structure that made it
realizable.

> **About 56% of the distogram's error is geometrically realizable. The predicted distance matrix is
> not a noisy version of the right structure — it is a good description of a WRONG structure, and the
> fit reproduces that wrong structure faithfully.**

The negative correlation of −0.513 completes the picture: the part of the error the fit *cannot*
satisfy is the part anti-correlated with the truth. **The optimiser discards the incoherent component
and follows the coherent one** — precisely backwards from what would be useful.

### Why this unifies the sprint's negatives

Every intervention attempted on this channel operates on the error's **magnitude profile**:
reweighting, recalibration, separation debiasing, robust and redescending losses, outlier trimming.
A realizable error is a **direction**, not a magnitude, and it lies in the null space of all of them.
So the following stop being separate disappointments and become one prediction, confirmed five times:

- **recalibrating the `sd` changes nothing** — provably, since a weighted least-squares argmin is
  invariant under uniform rescaling of the weights (K2);
- **separation debiasing is worth little** — negative on ranking (N5); on the fit axis the +0.179 quoted above is from a pre-seed-fix run and does not reproduce (the current artefact gives −0.065), so it sits at the noise floor;
- **robust and redescending losses move the mean little**, because the errors are not outliers —
  outlier clipping is worth 0.12 Å (K7's reformulation is real but bounded);
- **the per-target scale correction has a 0.1 Å ceiling** (step 3);
- **fusing the pool channel buys little** despite genuinely decorrelated errors, because the gain
  goes as the square of the weaker channel's skill (N10, K6).

It also explains the ORACLE arm's extraordinary 0.611 Å: with zero error there is no wrong structure
to follow, and the representation does the rest.

### The prediction this makes, which is falsifiable and should be tested

> **No reweighting, loss-shaping, calibration, debiasing, or robust potential applied to *this*
> distogram can recover more than a few tenths of an ångström.** The only interventions that can work
> are (a) a genuinely more accurate predictor, or (b) an independent channel whose errors are
> realizable as a *different* wrong structure, so that the two disagree about which wrong structure
> to build.

Point (b) is the one live route and it reframes what fusion is for. Fusion has been justified in this
project by error decorrelation, and dismissed by the squared-skill arithmetic. The realizability
result gives a different and stronger reason to want a second channel: **not to average away noise,
but to break the coherence of a confident mistake.** Whether the retrieval pool's distances are
realizable as a *different* wrong structure than the distogram's is a direct, cheap measurement, and
it is the obvious next experiment.

### Honest scope

Step 2 is n = 6 and step 3 is n = 8; both are smoke reads and the full-instrument runs are queued.
Step 4 is n = 24. The interval on the decisive sign-destruction effect excludes zero even at n = 6,
and the direction of every arm is consistent across steps, but the magnitudes will move. The claim
being made is qualitative — that the dominant error mode is coherent and realizable rather than
noisy — and it is the magnitudes, not the direction, that the full runs will refine.

---

## K10. A NATIVE-FREE LAW FOR WHAT FUSION IS WORTH, from structural disagreement alone

**DEMONSTRATED**, 24 targets. This follows directly from K9 and it is, I think, the most immediately
useful thing the sprint has produced: it lets you predict how much combining two channels will buy
**before running the combination, and without the native**.

### The two channels build genuinely different wrong structures

K9 established that the distogram's errors are largely realizable — they describe a consistent wrong
structure. The obvious question is whether the retrieval pool's distances describe the **same** wrong
structure or a different one. Fitting each channel separately, on 24 targets:

| | |
|---|---|
| distogram fit vs native | 3.270 Å |
| pool fit vs native | 3.469 Å |
| **the two fits disagree with each other by** | **2.980 Å** (median 2.996) |
| pool's unrealizable / wrong ratio | **0.255** — the pool is *more* realizable than the distogram |
| raw per-pair error correlation between the channels | **+0.627** |

**The dissociation is the finding.** The two channels' raw per-pair errors correlate at **+0.627** —
strongly, which by the usual reading would say they share their mistakes and fusion is pointless. Yet
the structures they imply sit **2.98 Å apart**, nearly as far from each other as either is from the
native. Correlated errors, different wrong structures.

This is why error correlation has been the wrong statistic for this decision throughout the project.
What matters is not whether two channels err on the same *pairs*, but whether their **realizable
components point at the same wrong structure**. Those are different questions and here they give
opposite answers.

### The law

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT.** Three things in this subsection are wrong.
> (i) The law is **not ours**: it is the two-member **Krogh–Vedelsby (1995) ambiguity
> decomposition** (`s16/lit_FINDINGS.md` §B.5). (ii) The stated caveats — *"exact when the native
> is equidistant from both and the three points are coplanar"* — are **unnecessary**. The identity
> is exact in **any** inner-product space with **no** equidistance and **no** coplanarity
> assumption, provided `r` is the **quadratic** mean. Verified here to 2.65e−15 on all 126 real
> targets (`s16/retract_exact.py`). (iii) `r` below is the **arithmetic** mean, which is the wrong
> mean; the n = 126 correction is +0.162 → **+0.075 Å**. The n = 24 numbers in the table below were
> never recomputed and are superseded by the n = 126 line above; they are preserved, not endorsed.

If two fits sit at mean distance `r` from the native and `s` from each other, ~~elementary geometry~~
the two-member ambiguity decomposition puts their mean at

> ~~**d_avg ≈ √(r² − (s/2)²)**~~ &nbsp;→&nbsp; **d_avg = √( (r₁²+r₂²)/2 − (s/2)² )**, exactly

~~which is exact when the native is equidistant from both and the three points are coplanar.~~
Tested against the measured coordinate average:

| | value |
|---|---|
| mean single-channel error `r` | 3.370 Å |
| structural disagreement `s` | 2.980 Å |
| **predicted average** | **2.970 Å** |
| **measured coordinate average** | **3.105 Å** |
| prediction error | **+0.135 Å** |
| restraint-level fusion (inverse-variance weighted `d̂`) | 3.046 Å |
| gain over the better single channel | **−0.166 Å** |

**The law predicts the measured fusion gain to within 0.135 Å from the disagreement alone.**

### Why this is worth having

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT. THIS CLAIM IS REFUTED, and it was tested rather
> than repeated** (`s16/retract_screen.py`, n = 126, no refit).
>
> `s` is native-free, but the law's *prediction* is not: it needs `r`, a mean RMSD **to the native**.
> The screen therefore cannot be run on unlabelled data as written. Tested anyway, three ways:
>
> | | value |
> |---|---|
> | corr(`s`, gain of the average over the **quadratic-mean member**) | **+0.777** — but this is *algebra*, not prediction: the identity makes that gain a deterministic function of `s` and `r` |
> | corr(`s`, gain over the **better single channel**) — the only decision a practitioner faces | **+0.112** (Spearman **+0.040**) |
> | AUC of `s` as a ranker of "does fusion win on this target" | **0.401** against a 0.500 null — `s` points the **wrong way** (mean `s` = 2.947 where fusion wins, 3.493 where it loses) |
> | leave-fold-out threshold "fuse iff `s > t`", *t* chosen **supervised** on 4 folds | every fold chose the grid minimum, i.e. **"always fuse"**; +0.0002 [+0.0000, +0.0006] over always-fuse |
>
> **Mechanism.** The ambiguity decomposition guarantees the ensemble beats the *average* member by
> the diversity, for **any** `s > 0`. It is completely silent on beating the **best** member — which
> is the only question a screen is for. So "`s` tells you whether fusion is worth running" is not a
> weak claim, it is the wrong claim: the identity answers a question nobody asks.
> (An ORACLE per-target pick would buy 3.367 → 3.242 Å, so headroom exists; `s` does not reach it.)
>
> **ORIGINAL TEXT, PRESERVED:**

~~**`s` is native-free.** The disagreement between two fits is measured without any reference to the
truth. So the law converts "should we fuse these channels?" from a question requiring a labelled
experiment into a **one-line calculation on unlabelled data**, and it says what the answer depends
on: not accuracy, not error correlation, but **structural disagreement relative to accuracy**.~~

It also states the condition for fusion to be transformative rather than marginal. The gain is
`r − √(r² − (s/2)²)`, so it is small whenever `s ≪ r` and only becomes large as `s → 2r` — two
channels making *opposite* errors, whose midpoint lands on the answer. Here `s/r = 0.88`, which buys
about 0.2 Å and nothing more. **To get an ångström from fusion you would need channels that disagree
about twice as much as they are wrong**, and no pair of channels in this project comes close.

That is a much sharper statement than the standing heuristic that fusion gain goes as the square of
the weaker channel's skill. The heuristic is a correlation-space argument and it gets the right
answer here for the wrong reason; the geometric law says *why*, is measurable in advance, and makes a
falsifiable numerical prediction rather than a qualitative one.

### A secondary result worth noting

**Restraint-level fusion (3.046 Å) slightly beats coordinate-level averaging (3.105 Å).** Combining
the two channels' distances and fitting once is a little better than fitting twice and averaging the
structures — presumably because the single fit is constrained to be realizable throughout, whereas an
average of two realizable structures need not be. The difference is small and n = 24, so it is
recorded as a direction rather than a result, but it supports the cascade's use of the combined
likelihood over a two-stage alternative.

### Scope

n = 24 targets, one channel pair, single fits rather than ensembles. The law's functional form is
geometry and should hold generally; its calibration here rests on one pair of channels and the
coplanarity assumption will not hold exactly. It should be checked on a second channel pair before
being quoted as general, and the full-instrument `errstruct` and `scale` runs will supply the wider
error-structure context it sits in.

---

## K11. The projection gap is real and worth 0.275 Å — but the distogram cannot be used to close it, and my 25.8% figure was wrong

**DEMONSTRATED**, all 126 targets, `s15/expand.py`. Two results and one correction to something I
asserted repeatedly earlier in this sprint.

### The correction, first

I wrote, in `GEOMETRY.md`, `ARCHITECTURE.md`, K8 and the `expand.py` docstring, that **coordinate
averaging contracts the backbone by 25.8%**, reusing a Sprint 14 figure. Measured directly on the
production consensus (`avg_ca` from every shipped record) against the true distances:

> **The production consensus is contracted by 3.5% against the true distances, and by 9.5% against
> the predicted ones — not 25.8%.**

I do not know what the 25.8% referred to; it may have been a different stage, a different set, or a
different measure. **What I do know is that I reused it without checking it, in four places, and used
it to motivate an experiment.** The number is withdrawn wherever I used it. The *direction* survives —
the consensus is genuinely contracted, and the projection genuinely pays to re-expand it — but the
magnitude was overstated by a factor of seven, and the mechanism story I built on it was
correspondingly overconfident.

### The result

| arm | mean | median | <2 Å | <2.5 Å | vs `project_raw` |
|---|---|---|---|---|---|
| `project_raw` — **is** the production `fit_ca` | **3.204** | 2.966 | 0.28 | 0.35 | — |
| `project_scaled` — native-free, distogram scale | 3.327 | 3.053 | 0.27 | 0.32 | **+0.123 [+0.065, +0.184]** |
| `project_scaled_deb` — native-free, debiased | 3.289 | 3.067 | 0.26 | 0.33 | **+0.085 [+0.040, +0.131]** |
| **`ORACLE_scale_true`** — the true distance scale | **2.929** | 2.811 | 0.28 | 0.37 | **−0.275 [−0.397, −0.174]** |

`project_raw` reproduces the incumbent's 3.204 Å exactly, which validates the whole construction —
the arm read from the cached `fit_ca` and the live projection agree.

**There is 0.275 Å on the table.** Rescaling the consensus by the correct isotropic factor before
projecting improves the *production* pipeline by −0.275 [−0.397, −0.174], with an interval nowhere
near zero. This is not an experimental arm; it is the shipped pipeline's own final stage.

**And the distogram cannot supply the factor.** The scale the truth wants is **1.0412**; the
distogram asks for **1.0997**. It over-expands by more than twice the required amount, because its
own **+0.509 Å expansion bias** is baked into the estimate. Applied, it costs **+0.123 Å**. Debiasing
halves the damage (+0.085) without reversing it, exactly as expected if the bias is the cause and the
leave-fold-out profile removes only its population-average part.

### Why this is K9 again, in a fourth place

This is the **fourth** distinct measurement in which the distogram's coherent expansion bias
contaminates a derived quantity and turns a sound intervention negative:

1. maximum likelihood under its own distribution loses to least squares on its mean;
2. the leave-fold-out separation debias is worth little on both axes;
3. the per-target scale correction has a 0.1 Å ceiling and the pool cannot estimate it;
4. **the consensus rescaling has a 0.275 Å ORACLE ceiling and the distogram overshoots it.**

Each was a reasonable intervention, each was defeated by the same coherent error, and each was
defeated in a way that a magnitude-based diagnostic would not have predicted. That is the practical
content of K9: **the error is a direction, and every correction this project has tried is a
magnitude.**

### K11-COMPLETE. Every native-free scale reference fails, and the reason is that the scale is unpredictable

**DEMONSTRATED**, all 126 targets. The pool was the one remaining candidate — K3 measured its
distances as having the **opposite bias sign** to the distogram's, so if the distogram over-expands
the scale estimate the pool should under-expand it and their average should land near the truth.

| arm | mean | median | <2 Å | <2.5 Å | vs `project_raw` |
|---|---|---|---|---|---|
| `project_raw` — **is** the production `fit_ca` | **3.204** | 2.966 | 0.28 | 0.35 | — |
| `project_scaled` distogram | 3.327 | 3.053 | 0.27 | 0.32 | +0.123 [+0.065, +0.184] |
| `project_scaled_deb` debiased distogram | 3.289 | 3.067 | 0.26 | 0.33 | +0.085 [+0.040, +0.131] |
| `project_scaled_pool` **the pool** | 3.299 | 3.005 | 0.25 | 0.33 | **+0.095 [+0.033, +0.158]** |
| `project_scaled_both` **pool + debiased distogram** | 3.284 | 3.069 | 0.25 | 0.33 | **+0.080 [+0.035, +0.129]** ‡ |
| **`ORACLE_scale_true`** | **2.929** | 2.811 | 0.28 | 0.37 | **−0.275 [−0.397, −0.174]** |

**All four native-free arms are on the wrong side of zero**, every interval excluding it. The pool
does not rescue the estimate, and neither does averaging it with the debiased distogram.

**Why, and it is not a subtle reason.** Every native-free estimator of the per-target scale is
**uncorrelated with the truth**:

| estimator | correlation with the true scale | MAE | mean scale it asks for |
|---|---|---|---|
| distogram | **+0.027** | 0.164 | 1.0997 |
| debiased distogram | **−0.001** | 0.163 | 1.0687–1.09 |
| retrieval pool | **+0.106** | 0.167 | 1.0799 |
| pool + debiased distogram | **+0.070** | 0.159 | 1.0687 |
| *(the truth)* | — | — | **1.0412** |

The per-target scale carries real information — correcting it with the true value is worth
**−0.275 Å** — and **no channel available to this project can estimate it at all**. Every estimator
also over-expands, so applying any of them makes things worse in a consistent direction.

**Method validated.** The cached-projection shortcut is exact: a live projection reproduces the
stored `fit_ca` at **max absolute coordinate error 0.0** on all four checked targets, so
`project_raw` genuinely *is* the production pipeline. And the isotropic-rescaling ceiling measured
without any projection is **0.216 Å** (3.048 → 2.833), consistent with the 0.275 Å measured through
it.

**Verdict: the projection gap is real, its ceiling is 0.275 Å, and it is closed.** This is the fifth
place where the distogram's coherent error defeats a sound intervention, and the first where the
failure is not "the estimate is biased" but "**the quantity is not predictable**" — which is a
stronger and more final form of the same obstacle.

---

## K8. THE FULL CASCADE DOES NOT BEAT THE INCUMBENT. The central architecture, measured

**DEMONSTRATED, all 126 targets.** `s15/cascade.py`, combined channel, leave-fold-out separation
debias, 16 native-free starts per target, clash gate at 2% of pairs below 3.6 Å. Two structures
excluded by the gate across the whole instrument.

This is the sprint's central table and it is a negative. It is reported first and without softening.

| stage | mean | median | <2 Å | <2.5 Å | FAIL18 | vs incumbent (paired) |
|---|---|---|---|---|---|---|
| **G** best in ensemble — **ORACLE** | 2.760 | 2.547 | 0.35 | 0.49 | 5.155 | **−0.444 [−0.573, −0.325]** |
| **S** objective's argmin | 3.511 | 3.327 | 0.24 | 0.29 | 6.217 | **+0.307 [+0.178, +0.439]** |
| **A** coordinate consensus | 3.151 | 2.946 | 0.27 | 0.37 | 5.838 | −0.053 [−0.125, **+0.020**] |
| **A** over objective-best half | 3.206 | 2.962 | 0.27 | 0.38 | 5.957 | +0.002 [−0.086, +0.088] |
| **F** after projection — *the like-for-like number* | 3.321 | 3.156 | 0.26 | 0.35 | 6.042 | **+0.117 [+0.050, +0.188]** |

gaps: **G→S +0.751**, **S→A −0.360**, **A→F +0.170**

### The verdict, stated so it cannot be misread

**The honest headline comparison is F, and F loses by +0.117 [+0.050, +0.188].** The incumbent's
3.204 Å is measured *after* the ideal-geometry projection, so only F is like-for-like. `A` at 3.151
is an unprojected coordinate average — it is not a physically valid structure and quoting it against
the incumbent would be comparing a projected structure to an unprojected one. Even taken at face
value its interval **includes zero** (+0.020 at the upper end), so it is parity, not a win.

Restraint-driven generation with the best channel combination, leave-fold-out debiasing, sixteen
diverse starts and coordinate consensus **does not beat a BLOSUM retrieval pipeline.** Sub-claim 3
of the provisional thesis — that framing the problem as *solving* rather than *generating and
selecting* bypasses the discrimination bottleneck — is **not supported** in this form.

### What the decomposition says, which is more useful than the verdict

**The ensemble is good and the readout is bad.** G = 2.760 beats the incumbent by −0.444 with an
interval nowhere near zero. Structures better than what the incumbent emits **are being generated on
most targets**; 49% of targets have an ensemble member under 2.5 Å and 35% under 2.0 Å. Every
ångström of the loss is in reading that ensemble out.

**Selection is actively harmful and aggregation partly rescues it.** G→S is **+0.751 Å** — choosing
the objective's argmin from sixteen fits throws away three quarters of an ångström. Aggregation
recovers **−0.360 Å** of it. This independently reproduces the standing law that aggregation is worth
several times any objective, now measured inside a generative architecture rather than a retrieval
one.

**Using the objective to pre-filter the ensemble makes aggregation worse.** Averaging the
objective-best half gives 3.206 against 3.151 for averaging everything. The objective's within-
ensemble ordering is worse than useless: discarding the members it ranks lowest *removes useful
structures*. This is the in-band blindness of K6 (native at the 34.7th percentile) appearing in a
second, independent place, and it closes the obvious "just filter the ensemble" response.

**The projection costs +0.170 Å**, and that is now the largest single remaining gap in the
cascade — larger than the whole S→A recovery. The consensus structure is off-manifold and paying to
be put back on it. Sprint 14 measured a related price and initially misattributed it; the standing
correction is that it is the cost of **re-expanding a contracted structure** (the contraction is
3.5% against true distances -- K11 withdraws the 25.8% figure I had been reusing). That is a specific, attackable mechanism rather than a
tax.

### The collapse check, run because Q1.3 would otherwise invalidate the A stage

The QGEOM workstream measured a CVaR-VQE whose "set mean" looked 0.85 Å better than its control
while its actual coordinate average was worse, because it had collapsed to **2.7 distinct structures
out of 4,096 draws**. The set-mean law is out of domain on a collapsed set, and K8's A stage is an
average over sixteen fits, so the check is mandatory rather than optional.

**It passes.** The fit ensemble's mean pairwise RMSD is **0.527 Å** (median 0.477), and only
**17 of 126** targets fall below a 0.25 Å spread. The ensemble mean RMSD is 3.528 Å against a best
of 2.760 Å, so there is real spread between members and genuine material to average. The aggregation
result — S→A −0.360 Å — is inside the law's domain and stands.

### Why this is a coherent negative rather than a puzzling one

The incumbent averages **75 retrieved windows**; the cascade averages **16 fitted structures**. The
retrieval pool simply contains better material: its top-75 best is **2.306 Å** against the fit
ensemble's ORACLE best of **2.760 Å**. Generation is not producing structures the library cannot
match — on this evidence, it is producing slightly worse ones, more slowly.

That points at the combination rather than the replacement, which was already queued before this
result landed: the fits are **0.43 Å better than the retrieval set mean** they would join
(3.644 against 4.072), and the terminal operator consumes the set **mean**, not the set best. Mixing
solved conformers into the retrieved pool needs only that inequality to hold, not the generative arm
to win. Its stated failure mode remains that both share the distogram, so their errors may not be
independent.

### What is now falsified, and what is not

**Falsified:** that restraint-driven generation, on its own, improves on retrieval in this regime.

**Not falsified, and strengthened:** that the barrier is restraint accuracy. G = 2.760 from
predicted restraints against 0.611 Å from true ones (K1-CORRECTED) puts the entire deficit in the
restraints, measured now through two different readouts on the same machinery.

**Newly sharpened:** the terminal operator, not the objective, is where the remaining leverage is.
Three of this sprint's measurements now say the same thing from different directions — S→A −0.360,
the best-half arm's +0.055 penalty, and A→F +0.170.

---

## K6. The restraint objective's minimum is BETTER than random — the direct rebuttal to Roget et al.

**DEMONSTRATED**, all 126 targets, K=500 pool members each, no optimisation involved.
`s15/feasible.py`, `s15/results/feasible.json`.

The experiment needs no search at all: evaluate every candidate objective at the native CA trace
and at all 500 retrieved windows, then ask where the native ranks and what the objective's own
argmin is worth. `native pct` is the native's rank among pool members **on the objective**, so 0.5
means the objective is blind to nativeness and above 0.5 means it actively disfavours the native.

| arm | native pct | median | is argmin | argmin RMSD | vs random | wins | ρ(f, RMSD) |
|---|---|---|---|---|---|---|---|
| `ls_pred` distogram, 1/sd² | 0.347 | 0.312 | 0.032 | 3.504 | **−0.949** | 0.79 | **+0.562** |
| `ls_debias` + separation debias | 0.359 | 0.324 | 0.032 | 3.559 | −0.895 | 0.76 | +0.557 |
| `ml_pred` full 17-bin likelihood | 0.394 | 0.330 | 0.016 | 3.450 | −1.003 | 0.83 | +0.532 |
| `ls_pool` pool median distances | 0.564 | 0.577 | 0.000 | 3.757 | −0.696 | 0.71 | +0.401 |
| `ml_pool` pool histogram | **0.665** | 0.722 | 0.000 | 4.026 | −0.427 | 0.56 | +0.365 |
| `ml_comb` distogram × pool | 0.463 | 0.434 | 0.000 | **3.421** | **−1.032** | **0.85** | +0.517 |
| **`ORACLE_true`** true distances | **0.000** | 0.000 | **1.000** | 2.006 | −2.447 | 1.00 | +0.862 |

**The machinery check passes exactly.** `ORACLE_true` places the native at percentile **0.000 on
126/126 targets**, is its own argmin on 126/126, and has ρ = +0.862 with a positive ρ on **100%**
of targets. An evaluation that could not do this would be measuring its own bugs.

**RESULT 1 — the objective's argmin beats a random pool member by −0.949 Å [−1.147, −0.753]**,
winning on 79% of targets, with ρ(f, RMSD) positive on **88.1%** of them. Fusion does slightly
better: `ml_comb` gives **−1.032 Å [−1.213, −0.859]** and wins on 85%.

This is the **direct rebuttal to Roget et al.** (arXiv:2606.21241), and it is like-for-like in the
one way that matters. They enumerate every peptide up to length 15 across 12,446 PDB structures
and report that **the minimum-cost conformation has on average a *larger* RMSD than a randomly
selected feasible one** — for a Miyazawa–Jernigan lattice contact potential, whose parameterisation
on proteins above 50 residues they identify as the cause. We are at the same chain lengths, on the
same kind of question, and our restraint objective's minimum is **about an ångström *better* than
random with a confidence interval nowhere near zero**. Their pathology is a property of contact
potentials at short chain length, not of short-chain structure prediction. That distinction is
exactly what makes our claim survivable next to theirs, and it is now measured rather than argued.

**RESULT 2 — but the native itself sits at the 34.7th percentile, and is the argmin on only 4 of
126 targets.** So the objective orders the pool *broadly* well while failing to put the truth at
the top. It is a global signal, not an in-band one — the same shape this project has found in every
channel it has ever measured, and it independently reproduces the standing result that had the
native at the **36.8th percentile and the argmin on 3/126** through completely different code. Two
implementations agreeing to within a percentile and one target is the best kind of confirmation.

The practical reading: mechanism (A) from the module docstring, not (B). Optimising the restraint
objective does **not** move away from the answer — it moves toward it, on average, by about an
ångström. What it cannot do is find the answer, because the native is not where the objective's
minimum is. Search is worth having and is not enough, which is a materially different and more
hopeful situation than Sprint 14's selective objectives, where search was strictly harmful.

**RESULT 3 — the pool channel is a BETTER estimator and a WORSE objective.** K3 measured the pool's
distances as more accurate than the distogram's (MAE 2.249 vs 2.386) with no training. As a
*ranking objective* it is the worst arm in the table: it puts the native at the **66.5th
percentile** — above chance, so it actively disfavours the native — with ρ of only +0.365 and its
argmin barely better than random (−0.427, 56% wins). This is a clean new instance of the standing
law that **MAE does not price selected RMSD**, and it is a useful warning: K3's "second independent
channel" is real as information and must not be assumed to be usable as an objective.

Fusion nonetheless helps: `ml_comb` has the best argmin of any native-free arm despite containing
the worst single channel. The gain is small, as Sprint 14's fusion arithmetic predicts.

**RESULT 4 — separation debiasing helps the FIT and slightly hurts the RANKING.** `ls_debias` is
marginally worse than `ls_pred` on every column here (native pct 0.359 vs 0.347; argmin 3.559 vs
3.504). A uniform additive shift of the restraint targets barely moves the *ordering* of candidates
while it does move the *location* of the minimum, so a correction that is worth having in
generation need not be worth having in selection. Both must be measured on their own axis; neither
result licenses the other.

## K7. Family B's feasible set does not contain the native at any usable ε

**DEMONSTRATED**, same run, 126 targets. This falsifies the natural formulation of Family B before
a line of it is written, which is the cheapest possible way to learn it.

Family B was declared as `min E_AMBER subject to C_restraint ≤ ε` — restrained refinement as
practised in NMR structure determination. Its feasibility question is whether the native satisfies
the predicted restraint bands at all. Measuring `z = |d_native − dhat| / sd` over all 8,549 pairs:

| ε (in units of the predictor's own sd) | fraction of pairs containing the native |
|---|---|
| 0.5 | 0.248 |
| 1.0 | 0.458 |
| 2.0 | 0.712 |
| 3.0 | 0.840 |
| median z | **1.336** |
| max z | **7.237** |

**At a one-sigma band, fewer than half the restraints admit the native.** Even at three sigma —
already so loose that the restraints barely constrain anything — one pair in six still excludes it,
and the worst pair is off by more than seven of the predictor's own standard deviations.

Two things follow, and they pull in opposite directions, so both must be said.

**Against Family B as declared:** a uniform-ε hard-constraint formulation has an empty or
truth-excluding feasible set at every ε that meaningfully constrains the search. The augmented
Lagrangian would drive the structure toward a region the native does not occupy. As specified, the
family is falsified.

**In its favour, and this is the more interesting half:** the failure is *concentrated*. The median
z is 1.336, so a typical restraint is only mildly inconsistent with the truth; the damage is done
by a minority of badly wrong pairs, one of them at z = 7.2. That is the signature of an **outlier**
problem, not a calibration problem — and outliers are exactly what restrained refinement in NMR is
built to handle, through per-restraint robust potentials and violation-aware reweighting rather
than a uniform band. So Family B is not dead; it is **redirected**, from a uniform hard constraint
to a robust one, and it now has a specific, measured quantity to attack. It also explains something
K2 could not: recalibrating the `sd` by a uniform 2.633× factor changes no structure, because a
uniform rescaling cannot separate a mildly-wrong restraint from a catastrophically-wrong one.

---

## K5. The project's long-standing accuracy ceiling was a property of the LIBRARY, not of the distance channel

**DEMONSTRATED**, by combining K1-CORRECTED with two standing project results. This is the most
consequential structural finding of the sprint so far, and it is a *reframing* rather than a new
measurement — which is exactly why it went unnoticed for so long.

Two results have been on the books for several sprints:

- *"The distance prior is the ceiling"* — oracle distances give **0.36 Å pool-best** and
  **0.98 Å selected**, through the retrieval library.
- *"The sequence signal is the ceiling"* — perfect distance knowledge caps the delivered answer at
  **~1.95–2.0 Å**.

Both were measured with the library in the loop: distances were used to **choose among, filter, and
weight retrieved windows**. The 1.95–2.0 Å figure is therefore the ceiling of *selection through a
finite library*, not the ceiling of the distance channel itself. The distinction was never drawn,
and the 2.0 Å number has been treated as a physical limit of the information available.

K1-CORRECTED measures the same channel without the library: continuous torsion-space distance
geometry, solving for a conformation rather than selecting one, reaches **0.611 Å mean, 0.055 Å
median, 86% under 2 Å** on all 126 targets from the true distance matrix.

> **Perfect distance knowledge is worth ~0.61 Å, not ~1.95 Å. The library, not the distance
> channel, was carrying the ceiling — and generation removes it.**

**What this changes.** Every prior estimate of the value of improving distance prediction was
computed against a 2.0 Å floor, which made the whole channel look nearly exhausted: an improvement
that halved distance error could at best recover a fraction of a 1.2 Å gap. Against a 0.61 Å floor
the same improvement addresses roughly **2.07× the headroom**, against a floor 3.19× lower. Fourteen sprints of "the prior is the
ceiling, stop working on it" rested on a measurement of the wrong operator.

**What this does NOT change, and must not be overclaimed.** It does not say the distance channel can
be improved, only that improving it is now worth much more. On the full instrument the raw predicted
arm is **3.644 Å and loses to the incumbent by +0.440 [+0.290, +0.592]** (K1-CORRECTED), so today
the generative route is *behind*, not ahead. The headroom is real; the means to reach it is not yet
demonstrated. Anyone quoting the 0.611 Å without the 3.644 Å beside it is quoting an oracle as a
result, which the brief forbids.

**The one experiment that prices this.** What matters operationally is not the ceiling but the
*transfer function*: how much predicted-distance accuracy buys how much RMSD. `s15/distacc.py`'s
phase diagram measures exactly this by corrupting true distances under realistic error models
(`abs_gauss`, `rel_gauss`, `sep_scaled`, `biased`, `outliers`) at a ladder of severities, and
inverting it gives **the distance accuracy required to reach 2.5 Å and 2.0 Å**. That converts the
whole research programme into a stated, falsifiable engineering requirement, and it is queued to run
the moment the cascade releases the CPU.

## K4 [RETRACTED -- see K4-RETRACTED below. Left unedited, per the brief.] The distogram's 17-BIN OUTPUT REPRESENTATION appears to cost ~1.06 Å, and nobody has ever measured it

**ORACLE DIAGNOSTIC, 8 targets, matched.** Handed to the RESTRAINT agent for a proper 40–126 target
ladder; recorded here at the moment of discovery so the provenance of the idea is dated.

The distogram's bins are **non-uniform** and coarsen severely with distance:

```
centres  4.00 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.50 9.50 10.50 11.75 13.25 15.00 17.50 21.00 25.00
width      —  0.75 0.50 0.50 0.50 0.50 0.50 0.50 0.75 1.00  1.00  1.25  1.50  1.75  2.50  3.50  4.00
```

There is **no bin below 4.0 Å**, and beyond 13 Å the resolution is 2.5–4.0 Å per bin. Long-range
pairs — precisely the pairs that fix the global fold — are therefore only resolvable to about ±1.5 Å
**by construction**, before the predictor makes a single error.

Two oracles that both know the true distance matrix exactly, on the same 8 targets, with the same
optimiser, starts, selection rule and geometry, disagree by 1.06 Å:

| ORACLE arm | representation | mean |
|---|---|---|
| `ORACLE_true_distances` (`s15/distgeo.py`) | continuous weighted least squares | **0.697** |
| `ORACLE_ml_true` (`s15/distml.py`) | ML against a 0.6 Å Gaussian **sampled at the 17 bin centres** | **1.756** |

Nothing differs between them except the output representation. So on this read the binning alone
costs **~1.06 Å of achievable accuracy — six times the entire measured contribution of any objective
in Sprint 14 (0.171 Å)**, and it has never been measured, discussed, or budgeted for in fourteen
sprints of work on this predictor.

**It also explains an otherwise puzzling result.** Maximum likelihood under the full distribution is
theoretically *strictly more informative* than least squares on its mean — least squares is the
Gaussian special case. Yet on predicted restraints ML is **worse**: `ml_full` 3.021 against
`ls_mean_sd` 2.877. The proposed mechanism is that the ML objective's minimum sits on the bin grid,
so ML inherits the binning as a hard resolution floor, while least squares consumes the
distribution's continuous *mean* and does not. ML's extra information is real but is being spent
paying a quantisation tax larger than the information is worth.

**HYPOTHESIS, now under test by the RESTRAINT agent:** a continuous density fitted to the histogram —
kernel-smoothed, a differentiated monotone CDF interpolant, or a small Gaussian mixture — keeps the
multimodality that motivates ML while restoring sub-bin resolution, and should let ML finally beat
least squares. *Falsified if* a continuous density does not recover the ORACLE gap, or if it recovers
the oracle gap but the predictive arms do not move — which would mean the predictor's own error
already swamps the bin width and the quantisation ceiling is inert in practice.

That last case is a live possibility and must not be glossed: the measured distogram MAE at
separation 11–15 is **4.666 Å** against bins 2.5–4.0 Å wide, so at long range the prediction error
and the bin width are the same order. The ceiling is only worth attacking where prediction error is
*below* it, which on this evidence is short and medium separation. The agent has been told to report
the short-range and long-range decomposition separately for exactly this reason.


---

## K4-RETRACTED. The "quantisation ceiling" was my own indexing bug, caught within the hour

**RETRACTED, same day it was written.** K4 above claimed the distogram's 17-bin output
representation costs ~1.06 Å, inferred from `ORACLE_ml_true` (1.756 Å) sitting 1.06 Å above
`ORACLE_true_distances` (0.697 Å) on matched targets. **That inference is void.** The number was
produced through a defective lookup that I wrote.

**The defect.** `s15/distml.py`'s `LogPTable` computed the bin index as
`floor((d − c[0]) / (c[1] − c[0]))` — a **uniform-grid** interpolation. The distogram's centres are
not uniform: they run 0.5 Å apart at short range and up to 4.0 Å apart at long range. So every
query above roughly 5 Å read the wrong bin, and not by a little:

| query | bin it read | that bin's centre |
|---|---|---|
| 6.25 Å | 3 | 5.75 Å |
| 8.50 Å | 6 | 7.25 Å |
| 11.75 Å | 10 | 10.50 Å |
| 15.00 Å | 14 | **17.50 Å** |
| 25.00 Å | 15 | **21.00 Å** |

Everything above 15.25 Å collapsed into a single bin. The likelihood was being evaluated against a
systematically wrong part of the distribution, worsening with distance — which is exactly the
signature I had just attributed to quantisation.

**What is void.** Every `ml_*` arm computed before 2026-09-05: `ml_full` (3.021),
`ml_full_debias` (2.814), `ml_plus_prior` (2.625), `ORACLE_ml_true` (1.756) in
`s15/results/distml.json`; the `pool_hist`, `pool_sim_weighted` and `combined` arms of
`s15/pooldist.py`, which wrap the same class; and the in-flight `s15/cascade.py` run, which was
**killed at 40/126 and relaunched** rather than allowed to finish and be quoted.

**What is NOT void.** `s15/distgeo.py` and `s15/distcal.py` never touch `LogPTable` — they are
weighted least squares on `dhat` and `sd`. So K1-CORRECTED (0.611 Å oracle, 3.644 Å predicted,
+0.440 vs incumbent), K2 (the bias and z-sd profile), K3 (the pool error profile) and K5 (the
ceiling reframe) all stand. The bug touches the maximum-likelihood branch only.

**The fix, and its check.** The lookup is now a `searchsorted` on the real centres with
per-interval widths, correct for any monotone grid, plus a vectorised `rowsum` for scoring whole
configuration spaces. The analytic derivative is re-verified against central differences:
**cosine 1.000000000000000, maximum absolute error 5e-09** — the same check K0 applied to the
torsion gradient, which is why it was cheap to run and should have been run on this class first.

**RESOLVED, by the RESTRAINT workstream (see R1).** The binned representation costs an ORACLE
**+0.381 Å [+0.156, +0.614]** — 111/126 targets worse, same sign in 5/5 folds, reproduced at +0.416
and +0.486 on independent start draws. **So 75% of the 1.06 Å I claimed was my bug**, and the
remaining quarter is real. It is bin **placement**, not count: a uniform 0.5 Å grid costs +0.243 and a
uniform 0.25 Å grid +0.085 with a CI crossing zero. A continuous density recovers 45% of it, and
maximum likelihood **still** does not beat least squares — the mechanism is calibration, not binning
(R1.4). The paragraph below was written before that measurement and is left as it stood.

**The open question as it stood before R1.** How much the binned representation actually costs was
*unmeasured*. It may be substantial, it may be nearly zero, and the honest prior is weak in both
directions: the bins are undeniably coarse at long range (2.5–4.0 Å), but the distogram's own MAE
at separation 11–15 is **4.666 Å**, so prediction error and bin width are the same order there and
the ceiling may simply be inert. The RESTRAINT agent has been given the corrected class and the
ORACLE representation ladder — continuous, snapped to bin centres, snapped to uniform 0.5 Å,
snapped to uniform 0.25 Å, and a Gaussian sampled at centres versus evaluated continuously — which
isolates bin **width** from bin **placement** from interpolation, reported separately for short and
long separation.

**The methodological lesson, which belongs in the paper's negative-results section.** A
representation defect produced a clean, plausible, mechanistically satisfying finding: a coarse
binning explaining why maximum likelihood underperformed least squares despite being strictly more
informative. The story was coherent, it fit prior results, and it was wrong. Nothing about its
internal plausibility flagged it; only checking the lookup against the grid it claimed to use did.
This is the second time in two sprints that a satisfying mechanism has turned out to be an
artefact of the instrument rather than a property of the problem.

---

## L1. THE NOVELTY POSITION HAS CHANGED MATERIALLY

Relayed from the LIT agent (`s15/LITERATURE.md`, 961 lines, 67 papers, sources tiered by
verification level). These are not stylistic notes: two of them **take** claims Sprint 14
intended to lead with, so the framing has to change now rather than at writing time.

**L1.1 — Certified-optimum-is-worse-than-random is TAKEN, for lattice contact potentials.**
Roget et al. (arXiv:2606.21241, 19 Jun 2026) enumerate every peptide up to length 15 across
12,446 PDB peptides and state that the minimum-cost conformation has, on average, a *larger*
RMSD than a randomly selected feasible one. They also supply the mechanism: the
Miyazawa–Jernigan contact potential was parameterised on proteins above 50 residues, and at
9–16 residues ρ(cost, RMSD) is negative. **What this leaves us:** theirs is a lattice contact
potential; ours is a **certified optimum of an all-atom force field (AMBER ff14SB) plus a
Legacy statistical potential, in continuous torsion space**. Their mechanism cannot cover ours —
ff14SB is not length-parameterised on large proteins. The claim survives only in its all-atom
form, must cite Roget as prior art for the phenomenon, and must be stated as an *extension to a
different potential class*, never as a discovery.

**L1.2 — In-loop CVaR for peptide folding is TAKEN.** Uttarkar et al., PLoS One 21(2):e0342012
(2026), "QuPepFold", uses CVaR inside the loop for peptide folding. Sprint 14's CVaR work cannot
be positioned as "CVaR applied to peptide folding." **What remains ours** is the **CVaR defect
triple** — the gradient-baseline defect (the baseline is centred on the tail only; cosine −0.023
with the true gradient) and its two companions in the Sprint 14 dossier. Those are mechanism
findings about the *estimator*, not applications of it, and no source was found describing them.

**L1.3 — The untrained-circuit best-of-N control is the strongest surviving novelty. Lead with
it.** The LIT agent found no paper anywhere that runs it. The nearest is Boulebnane et al.
(npj QI 9:70, 2023), which uses **uniform random sampling** — a strictly weaker control, because
it does not share the ansatz's support — and reports *parity*. We find **0 wins out of 12**
against best-of-N from the untrained circuit, at +0.65 to +1.32 Å. Stronger control, stronger
negative, no precedent. This is the paper's lead claim.

**L1.4 — "Nothing ranks within the pool" has an excellent citation.** Gulsevin & Meiler
(bioRxiv 2022.02.17.480937; Structure 31(1)) report that AlphaFold2's rank-0 model was the
lowest-RMSD one only **13%** of the time against a **20%** null — state-of-the-art ranking doing
worse than chance inside its own pool. This corroborates rather than takes.

**L1.5 — Tail-restricted discrimination is partially TAKEN as a framework.** Truchon & Bayly
(J Chem Inf Model 47:488, 2007) introduced BEDROC for exactly the early-recognition problem our
in-band metric addresses. We cite it and position our metric as an application, not an
invention.

**L1.6 — The only predictive comparator worth quoting** is Zhang et al. (arXiv:2510.06413) at
**4.89 Å mean on 75 fragments**. Our incumbent's 3.204 Å is on a different target set and is
*not* a like-for-like win; it must be stated as such.

**Four Sprint 14 claims must be removed from the story** (LITERATURE.md §5.4). The LIT agent
also self-corrected one of its own errors in-file: the AMBER/implicit-solvent literature does
**not** document ff14SB failing on peptides in the way an earlier draft asserted. Highest
priority unfetched source: PCCP 2018 / PMID 29480910.

---

## K15. POOL AUGMENTATION FAILS — and its premise was wrong, because I compared against the wrong set

**DEMONSTRATED, all 126 targets.** `s15/augment.py`. This was the last live candidate for the frozen
protocol, and it fails at its own precondition, which the module checked and printed before running a
single arm.

### The premise, and my error

The idea was that the terminal operator consumes the set **mean** (`d_out = 1.16·d_set_mean +
0.04·d_set_best`, R² = 0.89), so mixing conformers *better than the set mean* into the retrieval pool
must lower the emitted structure — **without needing the generative arm to beat the incumbent at
all**. I wrote that the fits were "0.43 Å better than the retrieval set mean they would join",
comparing the fits' 3.644 Å against **4.072 Å**.

**4.072 Å is the retrieval circular-mean *start* the fits are handed. It is not the set they join.**
The set they join is the production **top-75**, whose mean is **3.551 Å**. Measured directly:

```
top-75 set mean 3.551   fit mean 3.659   (the hypothesis needs fit mean < set mean: NO)
```

**The fits are 0.108 Å WORSE than the set they would be mixed into, not 0.43 Å better.** The premise
is false, and the sign of my motivating comparison was wrong because I compared against the wrong
quantity.

### The result, which follows

| arm | mean | median | <2 Å | <2.5 Å | vs `incumbent_avg` |
|---|---|---|---|---|---|
| `incumbent_avg` | **3.205** | 2.966 | 0.28 | 0.35 | — |
| `fits_only` | 3.605 | 3.352 | 0.14 | 0.22 | **+0.400 [+0.250, +0.560]** |
| `augment_w0.25` | 3.208 | 2.989 | 0.28 | 0.35 | +0.003 [−0.009, +0.015] |
| `augment_w1.0` | 3.224 | 3.079 | 0.27 | 0.36 | +0.019 [−0.001, +0.040] |
| `augment_w4.0` | 3.275 | 3.145 | 0.25 | 0.33 | +0.070 [+0.022, +0.119] |
| `augment_replace50` | 3.289 | 3.056 | 0.26 | 0.34 | +0.084 [+0.010, +0.156] |
| **`augment_LFO`** — the honest headline | **3.209** | 2.988 | 0.28 | 0.35 | **+0.004 [−0.001, +0.011]** |
| `augment_replace_LFO` | 3.224 | 3.079 | 0.27 | 0.36 | +0.019 [−0.001, +0.040] |
| **`ORACLE_bestfit`** — ceiling, never a prediction | 3.181 | 3.000 | 0.28 | 0.37 | −0.024 [−0.044, −0.005] |

**Every weight is monotonically harmful**, and the leave-fold-out procedure **correctly selected
w = 0 on four of five folds** (0.25 on the fifth) — i.e. asked to choose how much of the generative
ensemble to add, cross-validation chose *none*. That is the honest outcome and the procedure working
as intended.

**And even the ORACLE ceiling is small and below the floor.** Adding the single best-RMSD fit at high
weight is worth **−0.024 Å**, which sits below the sprint's measured ~0.08 Å false-positive floor. So
there was never more than a few hundredths available here, even with oracle knowledge of which fit to
add.

### The diversity check passes, so the law was in domain

Fit ensemble mean pairwise RMSD **1.397 Å** against the retrieval pool's 2.486 Å, with **56.3%** of
fit pairs separated by more than 0.25 Å. The ensemble is genuinely diverse, so the set-mean law was
applicable and the negative is not the collapsed-set artefact that invalidated an earlier quantum
reading.

### What this closes

**The last candidate for the frozen protocol has failed.** Nothing this sprint produced clears
§2.1 of the pre-registered rule against the incumbent, so the frozen protocol is the **incumbent
pipeline, unchanged**, and by the pre-registered §2.6 the benchmark is **not spent** — a confirmatory
run on an unmodified system could only return the `null` that was written down in advance.

**And the methodological lesson is the sprint's most frequent one, in a new place.** The experiment
was well designed: it checked its own precondition, reported diversity alongside RMSD, selected its
hyperparameter leave-fold-out, and included an ORACLE ceiling. **All of that was built on a
comparison I had made against the wrong baseline** — and the module printed the correct comparison in
its first line of output. A good design does not protect a wrong premise; it just makes the premise's
wrongness visible, which is what happened here.

---

## K14. PER-TARGET CONFIDENCE WORKS NATIVE-FREE — and my proposed new signal adds nothing to it

**DEMONSTRATED, all 126 targets, leave-fold-out.** No new computation; every quantity was already on
disk. **Recorded together with the control that refuted the version I was about to write**, because
the sequence is the point.

### What I set out to claim

K9 noted that the realizability measure has a **native-free half** — the converged fit's own restraint
residual `|d_fit − d̂|`, computable at inference — and K10 established the same of the **channel
disagreement** `s`. Neither had been tested as a signal. Both turned out to correlate with the fit's
own accuracy:

| signal | ρ with the fit's RMSD |
|---|---|
| fit residual `\|d_fit − d̂\|` | **+0.604** |
| channel disagreement `s` | **+0.607** |
| the two jointly, leave-fold-out | **+0.619** |

I drafted this as the sprint's first clearly positive native-free result.

### The control that killed it

Before recording it, I ran the comparison a reviewer would demand: does it add anything over the
signals the project already has?

| signal | ρ with the fit's RMSD | |
|---|---|---|
| **mean predicted `sd`** | **+0.635** | **already available; beats both new signals** |
| retrieval pool spread | +0.506 | already available |
| top-75 BLOSUM similarity | −0.424 | already available |
| fit residual (proposed new) | +0.604 | |
| channel disagreement (proposed new) | +0.607 | |

| leave-fold-out joint model | ρ |
|---|---|
| baseline: `sd` + pool spread + similarity | **+0.665** |
| the two new signals only | +0.619 |
| **baseline + the two new signals** | **+0.664** |

> **Adding both new signals to the baseline changes ρ from +0.665 to +0.664 — it adds nothing.** And
> the single simplest quantity, the mean predicted `sd`, already beats either new signal on its own.

**So the proposed contribution is REFUTED.** The realizability residual and the channel disagreement
are genuine signals, and they are redundant with a column the distogram has been emitting all along.

### What survives, and it is worth having

**Per-target confidence is achievable native-free, using only quantities the pipeline already
computes.** Leave-fold-out across the five pinned folds:

| | |
|---|---|
| ρ (Spearman) | **+0.665**, Pearson +0.629 |
| **best quartile** (n = 32) | **2.354 Å** |
| **worst quartile** (n = 32) | **4.884 Å** |
| **spread** | **2.530 Å** |

Even the mean predicted `sd` **alone** gives ρ = +0.617 and a 2.338 Å quartile spread.

**This does not improve a single prediction** — the emitted structures are exactly the ones the method
already produces; all it does is sort them. But it is a usable triage signal, it costs nothing to
deploy, and it is the natural companion to the sprint's central negative: *we cannot yet make these
predictions better, and we can tell you which of them to believe.*

### Why this works where in-band selection fails, and the distinction matters

**This is per-target confidence, not in-band selection.** Ordering *targets* by how well a method did
on them is a different and much easier problem than ordering *candidate structures within one
target's pool*, which fifteen sprints have failed at — the native sits at the 34.7th percentile of
every objective tried, and in-band accuracy tops out at 0.566 against a 0.638 requirement.

The two must never be conflated. What makes the per-target problem tractable is that its predictors
are properties of the **restraint set**, available before any structure is chosen; nothing here ranks
candidates.

### The methodological note

This is the fourth time in the sprint that a control inverted a conclusion — after the non-uniform
bin lookup, the per-process hash, and the per-target-minimum "control" that was an oracle. It is the
first time the control ran **before** the claim was written down rather than after. The difference
was one command and about two minutes.

---

## AL1. ALIGNMENT IS A DESCRIPTION, NOT A LEVER — and Family D stays empty

Relayed from the ALIGN workstream (`s15/align_FINDINGS.md`, 760 lines, artefacts in
`s15/results/align_{jac,fit,axes,free,reg}.json`). The instrument reproduced at start and end. **This
closes the sprint's last live hope for reaching 2.0 Å, refutes two of my own claims and one of the
information workstream's, and supersedes an experiment I had running.**

### AL1.1 There are FIVE exactly-null torsion directions, not four — and the fifth is derivable

**DEMONSTRATED on 126/126 targets**, with a clean gap of 1e-15 against 2e-3. An exact analytic
Jacobian `J = ∂(superposed CA)/∂(φ,ψ)` was built and verified at cosine 1.000000000000000 against
pair-distance central differences; it **independently reproduces every alignment figure in the
information workstream to three decimals** (0.565 / 0.665 / 0.693 / 0.709 / 0.738).

K0's four inert coordinate torsions are confirmed exactly (column norms 0, 0, 5e-14, 0). **The fifth
is new**: a *distributed whole-chain crank* with `Σδφ ≈ −Σδψ`, not a coordinate, which is why no
per-torsion perturbation test could have found it. And the count is not empirical but structural:
**`2n − 4` effective parameters mapping onto a `2n − 5`-dimensional CA shape space.**

The Cα trace responds to a **participation ratio of 3.47 effective directions out of 2n ≈ 25.9**;
90% of `trace(JᵀJ)` lives in **4.5** directions. An error placed in the quiet half tolerates
**15.7× more RMS angle** at the same RMSD cost.

### AL1.2 The mechanism is real and enormous — and unreachable

| | |
|---|---|
| **ORACLE**: rotate the fit's own true error out of the loud half at **fixed magnitude** | **1.855 Å**, i.e. **−1.822 [−2.037, −1.613]**, W/L **121/5**, 5/5 folds |
| best **native-free** leave-fold-out alignment arm, of 22 tried | **−0.007 [−0.074, +0.055]** |
| alignment actually moved by any arm | ≤ **0.104**, against the **0.387** required |

> **The information is present in the fit's own error vector — an oracle rotation clears 2.0 Å from
> the *same* restraints at the *same* error magnitude — and the best native-free rotation captures
> 4% of it.**

Pooled over 2,142 (target, arm) pairs: **partial corr(ΔRMSD, Δalignment | Δraw error) = +0.087**,
against **β = +0.586** for distance MAE. **Among achievable interventions, alignment is a description
of what good channels happen to have, not a lever that can be pulled.**

**Family D therefore stays empty.** It was declared unallocated, to be filled only if alignment
engineering worked. It did not, and the slot is left empty rather than filled with the arm that came
closest.

### AL1.3 The agent refuted its own headline, with the right control

Its geometry-aware metrics *appeared* to help — until a **matched-shrinkage** control
(`trace(A) = 2n`) was imposed. Then **isotropic Tikhonov beats both geometry-aware metrics at every
λ**, and beats the alignment-aware one significantly at the selected λ (**+0.124 [+0.025, +0.228]**).

So the workstream's one native-free gain — **`tik_LFO` −0.216 Å [−0.345, −0.094], W/L 80/46, 5/5
folds** — is **ordinary regularisation**, not alignment. It is above the sprint's ~0.08 Å
false-positive floor and is a real effect. **It does not qualify for the frozen protocol**, because
its control is the generative fit (3.677 Å), not the incumbent: it improves a losing arm to ≈3.46 Å,
still well short of 3.204. Its **λ grid is also not converged** — λ = 30 was chosen at the boundary on
5/5 folds — so even its magnitude is unsettled.

### AL1.4 THE TERMINAL-GAP LEVER DOES NOT TRANSFER — refuting the information workstream's largest lever

The information workstream measured that **where coverage gaps fall beats how many there are**:
terminal gaps at coverage 0.7 give 1.716 Å against 3.505 Å for mid-chain, a **1.79 Å** spread, and
called it "the largest lever in the diagram, and no channel exploits it." I relayed that as a live
opportunity and briefed an agent on it.

**It does not survive contact with a solver:**

| arm | effect |
|---|---|
| `ORACLE_term_native` — *perfect native* terminal torsions | **−0.000 Å [−0.068, +0.066]**, while cutting torsion RMS by 17° |
| dropping terminal restraints | **+0.210 [+0.055, +0.369]** — actively harmful |

**The 1.79 Å spread is a property of channels with holes, not of a solver.** A method that *fits*
rather than *selects* does not benefit from where its ignorance is located, and paying perfect
attention to the termini buys exactly nothing.

### AL1.5 ROBUST LOSSES FAIL ON THE FULL INSTRUMENT — K7's outlier reading is REFUTED

**`CAUCHY_lfo` −0.011 [−0.102, +0.077]; `gnc_cauchy` −0.025.** Both fire `s15/robust.py`'s own
pre-declared falsification condition, which read: *"If no robust loss beats `squared` on the full
instrument, then the outlier reading of K7 is wrong: the restraint errors are diffuse rather than
concentrated."*

**They do not, so it is.** K7's reformulation of Family B — from a uniform-ε hard constraint to a
robust potential — was the right *diagnosis* of the feasibility data (median z = 1.336, max z = 7.237)
and the wrong *prescription*. `robust.json` was only ever an 8-target smoke test, and the 8-target
"redescending signature" I wrote up (better median, more sub-2 Å, worse mean) does not survive.

**And the most informative detail:** the robust arms **do** reach far better restraint residuals
(median |z| **0.364** against 0.469) — and **convert none of it into accuracy**. Satisfying the
restraints better does not produce a better structure, which is K9 stated from yet another direction.

**I stopped my own full `robust.py` run on receiving this**, since the same arms had already been
measured at n = 126 inside `align_fit`. That decision and its reason are recorded in the scope notes.

### AL1.6 A native-free error-direction surrogate EXISTS — the one genuinely open lead

The agent expected none and found one: **`θ_fit − θ_pool` has |cos| 0.390 with the true error vector
against a 0.157 null**, with alignment correlation **+0.499 [+0.33, +0.65]**. The quiet subspace is
also identifiable native-free — the bottom-half subspace of `J(emitted)` overlaps `J(native)`'s at
mean cos² **0.827** against a 0.499 random null, with spectra correlated **+0.980**.

So both ingredients an alignment intervention needs — *which directions are quiet*, and *which way the
error points* — are separately estimable without the native. **It is untested as a steering signal**,
and it is the single cleanest open lead the sprint leaves behind.

### AL1.7 An incidental result that undercuts K2

**The debiased substrate emits 3.677 Å where undebiased `distgeo` emits 3.644 Å.** This is the first
full-instrument evidence that K2's separation debiasing **does not buy fit accuracy either** — it was
already known not to help ranking (N5). The difference is at the noise floor and not significant on
its own, but it points the same way as the verification workstream's independent finding that the
`debias_sep` surrogate row reverses sign on re-running.

### AL1.8 The phase diagram, closed

3.0 Å **supported**; 2.5 Å **still marginal**; **2.0 Å still not supported predictively — but the
reason has changed from information to control.** The information is demonstrably present (the ORACLE
rotation clears 2.0 Å from the same restraints at the same error magnitude); what is missing is any
native-free way to act on it. The alignment correction explains **60%** of the i.i.d. surface's
over-pricing, independently reproducing the information workstream's "two-thirds linear".

---

## R1. THE RESTRAINT WORKSTREAM: an empirical FALSE-POSITIVE FLOOR that binds the whole sprint

Relayed from `s15/restraint_FINDINGS.md` (675 lines, 14 JSONs, `s15/quant.py` at 1,047 lines). The
instrument reproduced the pinned constants at start and end; all derivatives verified at cosine
1.000000000000 including the torsion path.

**The most important item is §R1.5. It is a constraint on how every small effect in this sprint may
be reported, including several of mine.**

### R1.1 My retraction was right, and the corrected number is a quarter of what I claimed

The agent found the `LogPTable` bin bug independently, before running anything, and found a **second
defect I had missed**: above 15.25 Å the table collapsed everything into one interval *and returned a
non-zero derivative for a flat value* — cosine **0.8099** against central differences, with every
disagreeing sample above 16.13 Å. Its independent `LogPGrid` and my corrected class agree to
**max |Δ| = 0.0e+00**, which is genuine replication rather than agreement by construction.

The corrected measurement, within one start-matched run:

> **The pre-fix table cost +2.081 Å where the same information through a correct table costs +0.528.
> So 75% of the "1.06 Å quantisation ceiling" was my bug.** And the recorded "ML loses to least
> squares" (+0.230 [+0.093, +0.365]) **was the defect, not a finding**.

### R1.2 The real quantisation cost, and what caps any fix

**Binning costs an ORACLE +0.381 Å [+0.156, +0.614]**, 111/126 targets worse, same sign in 5/5 folds,
reproduced at two further independent start draws as **+0.416** and **+0.486**.

**It is bin PLACEMENT, not bin count.** A uniform 0.5 Å grid costs +0.243; a uniform 0.25 Å grid
costs +0.085 with a CI crossing zero. And the ceiling on any de-binning fix is set before a structure
is ever fitted: bin-snapping adds **0.604 Å RMSE against the predictor's own 3.704 Å — 2.7% of the
error variance**.

**Recommendation, adopted:** respace the bins uniformly at ~0.5 Å *if* the distogram is retrained;
otherwise keep consuming mean + sd, and spend the effort on calibration.

### R1.3 A continuous density recovers 45% of it — and ML still does not beat least squares

Kernel smoothing takes binned ML from +0.528 to **+0.291**; mass-exact PCHIP recovers **nothing**
(+0.504). **It is the smoothing, not the mass fidelity.** Both are ORACLE.

On predicted restraints, **every corrected ML arm is −0.041 to −0.102 with a confidence interval
touching or crossing zero**, and the honest leave-fold-out arm that also selects the density family
is the **weakest** of all (−0.040, 60 wins to 66, median **+0.025**). So the answer to "does a
continuous density let ML beat least squares" is **no**.

### R1.4 The mechanism is calibration, exactly as the information-channel workstream predicted

This is the cleanest mechanism confirmation in the sprint, and it was a pre-stated prediction rather
than a post-hoc fit.

ρ(predictive information, ML − LS) = **−0.292 (p = 0.001) / −0.217 (p = 0.014) / −0.186 (p = 0.037)**
across three arms, all outside their permutation nulls. And decisively:

> **ML beats least squares by −0.13 to −0.18 Å with confidence intervals excluding zero on exactly
> the 52 targets where the distogram's cross-entropy is positive, and by nothing on the other 74.**

The agent's reimplementation reproduces the information workstream's −0.818 bits/pair and its
52-of-126 split exactly. (The split is ORACLE, so this is a mechanism, not a method.)

**And the mechanism in the brief I wrote for this agent is FALSE.** I hypothesised that ML's
advantage would come from multimodality — contact versus non-contact structure. Strict bimodality is
**8.1%** of pairs and is *lower* at long separation (5.9%) than short (9.9%) — the opposite of what I
predicted. Per-target ML advantage is uncorrelated with bimodality (ρ = +0.058) and skewness
(−0.107), both inside permutation nulls, on all four arms.

### R1.5 A DEMONSTRATED FALSE POSITIVE, AND THE FLOOR IT ESTABLISHES

**ML against a constant-width Gaussian *is* least squares** — the comparison is provably zero by
construction. Measured:

| run | result |
|---|---|
| one start draw | **+0.081 [+0.014, +0.169]** — a confidence interval **excluding zero** |
| another start draw | −0.003 [−0.071, +0.060] |

> **A paired confidence interval excluded zero on a quantity that is exactly zero, purely from the
> multi-start draw.** The null-calibrated concentration test independently flags that row — and only
> that row — as **FAIL (p = 0.024)**, which is the check earning its place.

**This establishes an empirical false-positive floor of roughly 0.08 Å for single-draw paired
comparisons in this machinery.** Every predictive effect in the agent's own ML table (−0.040 to
−0.102) sits at or below it, which is why the honest answer to R1.3 is "no" rather than "a small yes".

**A LIMIT ON THE FLOOR, established afterwards and recorded here.** The floor's mechanism is a
**per-process-salted random multi-start**. It transfers only to pipelines that *have* one. The AMBER
k = 30 pipeline does not — it contains no RNG and reproduces bit-identically across interpreters — so
applying the floor to it, as I did below, was wrong (P1.1). **A floor measured on one pipeline does
not transfer to another by size alone; it transfers only if the mechanism does.** Check the mechanism
before importing the number.

**It binds my results too, and I am applying it rather than noting it.** Effects at or below the
floor, which must now be reported as indistinguishable from start-draw noise unless replicated across
draws:

| result | effect | status under the floor |
|---|---|---|
| `project_scaled_both` (K11) | +0.080 [+0.035, +0.129] | **at the floor** — the direction is corroborated by three sibling arms, but the magnitude is not resolvable |
| restraint-level fusion vs distogram (K10) | −0.071 [−0.156, +0.014] | already non-significant; now doubly so |
| cascade `A` vs incumbent (K8) | −0.053 [−0.125, +0.020] | already non-significant |
| AMBER relaxation at k = 30 (Sprint 14) | −0.022 [−0.036, −0.009] | **THE FLOOR DOES NOT APPLY — see P1.1.** This pipeline contains **no RNG**; it reproduces bit-identically across interpreters (max \|Δ\| = 0.000e+00, n = 126). The floor's mechanism is a per-process-salted random multi-start, which cannot act here. Its own constructed floor is ≈0.004 Å gated. *(The result is nonetheless WEAKENED, for unrelated reasons — P1.)* |

Effects comfortably **above** the floor and unaffected: the cascade's G→S (+0.751), S→A (−0.360),
A→F (+0.170) and F (+0.117); coordinate averaging vs the better channel (−0.255); realizability
(−1.391); the ceiling (−2.593); the quantisation cost (+0.381).

### R1.6 The start-draw variance, quantified

**sd 0.132 Å on an arm's mean across four draws**, with a worst-target spread of **5.15 Å**. So
**K1's headline 0.611 Å is one draw from that spread** — the agent's own draw returned **0.569 Å**.
Paired *differences* are far more stable (0.381 / 0.416 / 0.486 for the same quantity), which is why
the sprint's comparative claims survive and its absolute constants need an error bar.

**K1's ORACLE ceiling should be quoted as ≈0.6 Å with a start-draw sd of 0.13 Å, not as 0.611.**

And the agent's framing of why this was missed is exactly right: **it was invisible to every
determinism check the sprint had, because none of them starts a second interpreter.**

### R1.7 A correction to one of my own analytical claims

I wrote that long-range pairs carry nearly all of the quantisation loss. **Wrong in magnitude.**
Separation bands cost +0.232 / +0.196 / +0.358 *alone* but **sum to +0.785 against a joint value of
+0.381** — strongly **sub-additive**, because partially corrupting a restraint set leaves the bad
restraints fighting exact ones, which is a harder problem than corrupting all of them.

> **Band decompositions of geometric fits must print the joint value. Summing single-band effects
> overstates the total by a factor of two here.**

### R1.8 The agent's own retraction, preserved

It attributed a ladder gap to scipy's relative `ftol` and the Gaussian normalisation constant, and
built a dedicated control (`GaussDensityNoConst`) to prove it. **Its own control killed the
explanation** (−0.019 Å). The real cause was the start draw — which is how it found R1.5.

---

## P1. THE AMBER RESULT: WEAKENED on accuracy, CONFIRMED and understated on validity — and my floor did not apply

Relayed from the PHYS workstream (`s15/phys_FINDINGS.md`, module `s15/phys_repl.py`, 16 artefacts).
Instrument re-verified exact at the end. **This corrects an error of mine and finds three defects
nobody had recorded.**

### P1.1 MY OWN ERROR FIRST: the 0.081 Å floor does not reach this pipeline

I applied R1.5's empirical false-positive floor to the AMBER k = 30 result and wrote that −0.022 Å
"sits below the floor" and "must be replicated before it is quoted." **That was wrong, and the reason
is mechanical.**

> **The k = 30 pipeline contains no RNG at all.** `core.project.STARTS` is a fixed 4-tuple — `multi =
> True` means *all four fixed starts*, not random ones — AMBER runs at `Threads = 1` with
> `DeterministicForces`, and the only seeds are constants. Re-run in a **fresh interpreter** it
> reproduces the Sprint 14 artefact **bit-identically on both arms, n = 126, max |Δ| = 0.000e+00**.

The floor's mechanism is a **per-process-salted `hash()` driving a random multi-start**. There is no
such draw here, so the floor **cannot act on this pipeline and must not be imported into it**. A
floor measured on one pipeline does not transfer to another by size alone — it transfers only if the
mechanism does, and I did not check that before applying it in four documents.

**The pipeline's own floor had to be built**, and it is a different number for a different reason:
**≈0.004 Å gated, up to 0.038 Å ungated** (P1.3).

### P1.2 Replication holds

Five bootstrap redraws of the 75-window input: **−0.0145 / −0.0235 / −0.0207 / −0.0226 / −0.0354**.

| | |
|---|---|
| mean | **−0.0233 Å**, sd **0.0076 Å** |
| same sign | **5/5** (8/8 including the frame draws) |
| CIs excluding zero | 4/5 |
| against the pre-registered threshold | sd 0.0076 < **0.0085**; effect/SE = **2.25** |

The threshold was derived and written down **before any draw returned**. **Caveat, also
pre-registered and then confirmed:** r(projection cost, effect) = **−0.826** across draws, so the
bootstrap axis **flatters** the effect.

### P1.3 But the exact null fails, and it exposes a new defect

**Null A: relax the same structure in a rotated lab frame.** Zero by construction — ff14SB/GBn2, the
restraint and every RMSD are rigid-invariant. Measured: **+0.0117 [−0.0034, +0.0383], sd 0.137, max
1.51 Å**. On that frame the effect becomes **−0.0102 [−0.0332, +0.0215] — a confidence interval
including zero.**

**NEW DEFECT: `refine_coords` has no convergence gate.** Four targets — 1D6X, 1MF6, 2NB7, 7BX2 — end
minimisation **above 1000 kcal/mol**, with **1MF6 at 8.9 × 10⁸ kcal/mol** and 2 clashes, and are
**silently scored into the published mean**. Nothing anywhere in the pipeline catches it.

**Those four targets alone move the null from −0.0005 to +0.0117.** Gated, the null collapses to
**−0.0005 [−0.0046, +0.0035]** — a floor of ≈0.004 Å — and the effect reproduces across frames to
**0.0003 Å**.

**Null B: permute `fit_multi`'s four fixed starts.** Returns **exactly 0.000e+00** on 126 targets ×
3 permutations × 2 λ arms. **The agent's own hypothesis — that start permutation would flip optimiser
basins and reproduce the 0.08 Å mechanism — is refuted.** All irreproducibility lives in the AMBER
minimiser.

### P1.4 The damaging part: the effect is absent on 69% of the instrument

On the **84 targets (69%) where the minimiser is frame-reproducible**, the effect is
**−0.009 [−0.024, +0.006] with a *losing* 38W/46L** — replicated on **three independent frames**
(−0.009 / −0.005 / −0.008; 38/46, 36/49, 38/47).

> **The whole effect sits on the 31% of targets where the AMBER minimiser is not even
> frame-reproducible.**

And Sprint 14's "5/5 folds same sign" evidence **does not fully replicate** — bootstrap draw 1 gives
3/5.

### P1.5 The validity half is CONFIRMED, and it was badly understated — because of two Sprint 14 defects

Separating accuracy from validity turned up two reporting defects in the original:

- **Arm A's Ramachandran fraction was scored on a DIFFERENT STRUCTURE from its RMSD** — the λ = 0.3
  arm's geometry against the λ = 0 arm's RMSD.
- **Arm A's clash rate was a hard-coded placeholder**, `[0.0] * len`, not a measurement.

Corrected:

| quantity | measured against the correct baseline |
|---|---|
| Ramachandran-favoured | **0.466 → 0.874**, **116W / 2L** *(not 0.734 → 0.874)* |
| heavy-atom clashes | **1.397 → 0.000**, **63W / 0L** |
| minimum heavy separation | **+0.843 Å**, 3W / 123L |

**That half is nowhere near any noise floor.** The stereochemistry claim is robust and larger than
recorded; the accuracy claim is not.

### P1.6 Cost, and what the effect mechanically is

**No budget-matched comparison enters the −0.022 Å figure at all** — it compares two one-shot
post-processing operators. The audit's 28 ms → 8–14 ms correction moves nothing here, and the
single-point figure is the wrong number anyway: **this arm is a full minimisation at 12.57 s mean**,
about **1,500× a single point** and **2.8× the projection**.

Mechanically, the agent's reading is that the effect is **a price of validity, not an accuracy gain**.

### P1.7 The standing claim, restated

> **The accuracy claim may no longer be quoted as "−0.022 Å at valid geometry."** It replicates in
> sign and magnitude, it is not a multi-start artefact, and it is **absent on the frame-reproducible
> 69% of the instrument**, **erased by an exact null that must be zero**, and **gated behind a
> convergence check the pipeline does not perform**.
>
> **The validity claim is CONFIRMED, and larger than recorded:** 0.466 → 0.874 Ramachandran-favoured
> (116W/2L) and 1.397 → 0.000 clashes (63W/0L).

### P1.8 And a familiar trap, in a new place

The agent nearly cut the frame block to two draws after a **3-target probe** gave |Δ| ~ 1e−4. On the
full 126 the same null contains a **1.51 Å outlier** and a mean of +0.0117 Å. **A 3-target probe of a
heavy-tailed quantity measures the wrong thing** — the project's "n ≤ 4 has reversed a conclusion
four times" trap, appearing for the fifth time this sprint.

---

## Q6. THE QUANTUM QUESTION, CLOSED: the effect is real, replicates out of sample, and is not quantum

Relayed from the QENS replication workstream (`s15/qens_FINDINGS.md`). **Verdict: CONFIRMED as an
effect at power; REFUTED as evidence of anything the circuit does.** This supersedes Q5's demotion,
which was correct about the original inference and incomplete about the effect.

### Power, and an exact reproduction first

| | original | replication |
|---|---|---|
| targets | 9 | **19** (nine n = 9 plus ten n = 10) |
| seeds | 3 | **8** |
| cells per α | 27 | **152** |
| α values | 4 | 3 primary **+ a 10-point map** |
| new cells | — | **1,672** |

**All 108 original cells reproduced bit-for-bit — maximum absolute difference 0.000e+00 on every
readout.** The machinery is the original machinery, so nothing below turns on an implementation
difference.

### The effect is real, and it is stronger out of sample

At the **target-level** unit (n = 19), concentration check DIFFUSE:

| α | rand-5 coordinate average |
|---|---|
| 1.0 | **−0.385 [−0.748, −0.035]** |
| 0.25 | **−0.295 [−0.495, −0.069]** — where the original found *nothing* |
| 0.05 | −0.161 [−0.328, +0.008] — null |

**It is stronger on the ten targets never used to find it** (−0.413 / −0.357 / −0.206, all
significant) **than on the original nine** (all null). That is the out-of-sample confirmation that
matters, and it is why Q5's demotion was incomplete: the original *inference* was wrong, the
*effect* is not.

**Both defects Q5 identified are confirmed and corrected**: the 27 cells were 9 targets × 3 seeds
treated as independent, and `set_coordavg_rmsd` subsampled `np.unique(idx)[:400]` — **sorted, not
random** — biasing arms unequally. An unbiased readout was added.

### And the α non-monotonicity — my stated "reason to disbelieve" — is REFUTED

Ten α values × 19 targets: **−0.442, −0.565, −0.464, −0.310, −0.334, −0.302, −0.273, −0.295, −0.188,
−0.158.** Smooth, monotone, significant at 9 of 10. **The reported non-monotonicity was two noisy
three-seed cells.** Per-seed spread at α = 1 runs −0.09 to −0.65 (sd ≈ 0.19), so a three-seed study
could have reported almost anything.

### THE DECISIVE CONTROL: a classical reweighting of the same circuit wins, for 1/200th of the budget

> **A classical Boltzmann reweighting of the *same untrained circuit* by the *same objective* at
> matched entropy, using 1/200th of the objective budget, beats the CVaR-VQE by +0.25 to +0.48 Å on
> every readout at every α (W/L up to 1/18). Simulated annealing at 1/400th of the budget matches or
> beats it.**

The VQE beats only uniform random sampling — which this instrument's own history says proves nothing,
since a zero-information constant α-helix beats that control too.

**And on the optimisation axis the gap is starker still:** simulated annealing finds the **certified
global optimum in 100% of cells at 2,048 evaluations**, while the VQE puts **5.9% of its mass** on it
using **819,200**.

### What the effect actually is

**Location, entirely.** `d_coordavg = 0.94 · d_set_mean` at α = 1 with R² = 0.75 (1.01× and 1.20× at
α = 0.25 and 0.05), with **no residual shape term**. The VQE's distribution simply sits lower:
E_p[RMSD] is −0.834 / −0.474 / −0.291 below the control's.

**And the per-target win is explained by the objective, not the sampler.** It correlates
**+0.52 / +0.64 / +0.69** with where the objective's own certified argmin sits in the ORACLE RMSD
distribution. On the two targets where the argmin is at percentile 0.000 the VQE "wins" by −2.02 and
−1.83 Å; on the two where it sits at 0.939 and 0.938 it **loses** by +1.25 and +0.74 Å. The method is
being rewarded for targets on which the objective happens to be right.

**Matched diversity removes it** wherever the arms are genuinely diverse — −0.02 to −0.07 Å at
α = 0.25 and 0.05, at or below the 0.08 Å false-positive floor, reported **UNRESOLVED**. (At α = 1
matching degenerates to two distinct structures and *manufactures* a larger apparent win; the
agent's expected decisive test was the wrong one, and it says so.)

### The enrichment result collapses the same way

The sub-2 Å probability-mass enrichment from the other quantum workstream **replicates exactly**
(5.20 / 6.16 / 0.98 / 6.70 / 18.64×). But:

- **classical simulated annealing at the identical 8,192-evaluation budget gives 5.34 / 8.07 / 0.68 /
  8.00 / 24.31×** — null against the VQE on all four native-free objectives and significantly
  **better** on the ORACLE one;
- an entropy-matched classical use of the same objective roughly **doubles** it;
- and the headline 5.2× is a **ratio of means**: the **median per-target ratio is 1.13–2.66×, with
  only 5 of 9 targets showing any enrichment at all**.

Matched-diversity untrained sampling gives only 1.04–1.11×, so the enrichment is **not** free
diversity — it is the objective being exploited, and classical methods exploit it better.

### My own prediction, refuted

I predicted to this workstream that the win should track the **distinct-configuration count**, since
concentration costs the distinct samples a best-of-N readout is paid in. **Within-α ρ = +0.004 across
1,216 pooled cells**, and the aggregate direction is the *opposite* — the arms that win concentrate
**most**. The mechanism I proposed is wrong; the objective-quality mechanism above is what carries it.

### The closing position on the quantum component

**There is no quantum result, positive or negative, that survives its classical control.** What
survives is a methodological finding, and it is the strongest version of what this sprint has been
saying all along:

> Every apparent quantum effect measured here dissolves under a control that uses **the same
> circuit's samples, the same objective, and two to three orders of magnitude less budget**. The
> distribution the optimiser produces is worth something; the optimiser is not the cheapest way to
> produce it.

---

## Q4. THE QUANTUM ANSWER IS NOT A SIGN, IT IS A DEPENDENCE — and my own framing was overstated

Relayed from the QRESTRAINT workstream (`s15/qrestraint_FINDINGS.md`, 766 lines, five result JSONs).
**This corrects the framing I had been carrying in `PAPER_DRAFT.md`, `FINAL_DOSSIER.md` and
`VQE_CVaR.md`**, which stated the Sprint 14 negative as though it held unconditionally.

The agent extended the handed-over module from 9 to **19 targets**, added a **budget ladder**, a
**random-tail null**, **concentration measurement**, and — decisively — **`S14_disto_bayes`**, Sprint
14's own selective objective, so that the generative and selective classes could be compared *on one
instrument with one control*. It found **no arithmetic bug** in the code I handed over, and verified
the instrument by reproducing Sprint 14's flagship table to three decimals (0.655 / 0.549 / 0.509 /
**0.391** against a recorded 0.390).

### Q4.1 VQE beats the untrained control here — and the win is objective-independent, so it is not what it looks like

| | result |
|---|---|
| `E_ml_pred`, budget 8192 | **−0.202**, 13/19, CI excluding zero |
| `E_ls_pool`, budget 8192 | **−0.184**, 14/19, CI excluding zero |
| **`S14_disto_bayes`** — the *selective* control, budget 8192 | **−0.192, 15/19 — the STRONGEST cell** |

The largest gain comes from **the objective with the worst in-tail ordering measured anywhere in the
programme**. Whatever VQE is doing here, **it is not exploiting the generative objective's ordering
skill**, because it does the same thing — slightly better — on an objective that has none.

**So the answer to the sprint's central quantum question is: the Sprint 14 negative does not
reproduce on this instrument at this budget, and the reason is not the change of objective class.**

### Q4.2 The sign depends on the budget, and vanishes at the largest one

| budget | cells significant (of 6) |
|---|---|
| 2,048 | 1 |
| 8,192 | 3 |
| **32,768** | **0** |

Non-monotone, and **gone at the largest budget**. A single-budget report would have produced any of
three different papers.

### Q4.3 And VQE never beats a classical greedy search

**0 wins in 16 against greedy, losing 3 of 16.** Beating your own untrained initialisation is
necessary and nowhere near sufficient; the classical control is the one that decides whether any of
this matters, and it is not close.

### Q4.4 Sprint 14's mechanism REPRODUCES, and is not a selective-class pathology

`tail − bulk` ordering skill is **negative for every objective** — all four generative ones, the
selective control, **and the ORACLE**. The generative objectives do have better **in-tail** ρ (+0.127,
z = 13.1, against −0.042), but it never turns positive relative to the bulk; and **globally** the two
classes are close (+0.523 for `E_combined` against +0.485 for the selective control). *(An earlier
draft called +0.127 the global ρ. It is the in-tail figure; corrected.)*

> **The binding quantity is the selection gap: 1.47–2.14 Å, identical for VQE and control, flat
> across a 16× budget range, and collapsing to 0.29–0.44 Å under ORACLE restraints. It is restraint
> error, not the sampler.**

That last clause is the sprint's thesis arriving from the quantum side independently: the quantum
apparatus is not the bottleneck, and neither is the objective class. The restraints are.

### Q4.5 THE CVaR GRADIENT DEFECT HELPS — because a weaker optimiser is a better sampler

Sprint 14 recorded the CVaR baseline as **defective** (centred on the tail only; cosine −0.023 with
the true gradient), and I have been reporting it as a liability that weakens every negative result
obtained with it.

Measured directly: the **biased** estimator gives **−0.266** where the **corrected** one gives
**−0.145** on the ORACLE objective. The defect **helps**, and the mechanism is measured — the biased
gradient's norm is **2.4× smaller**, so it acts as a de-facto step-size reduction. **A weaker
optimiser is a better sampler**, because it concentrates less, and concentration is what costs the
distinct samples a best-of-N readout is paid in.

This is consistent with everything else in the sprint and it inverts a stated liability into a
mechanism. A third recorded defect (`dead_start`) is **inert**: 0.000 across all 19 × 3 × 3 × 8 runs.

### Q4.6 One clean positive, correctly priced

The trained state enriches sub-2 Å probability mass by **5.2–6.7×** over its own untrained
initialisation (**18.6×** under ORACLE restraints), with an internal control confirming it is real —
the weakest objective gives 0.98×, i.e. nothing.

**And it is worth ≤ 0.2 Å**, because an argmin readout ignores probability mass. That prices Family C
exactly: the quantity VQE genuinely improves is the one the terminal operator does not consume.

**The gap the agent named itself:** no classical ensemble sampler was run against this enrichment, so
it measures *what the VQE does*, not *what only a VQE can do*. Relayed to the replication workstream.

### Q4.7 What this changes in my own documents

**The abstract's "worse than not running it" is too strong as an unconditional claim and is being
revised.** The defensible statement is narrower and, I think, more interesting:

> Whether a variational optimiser beats its own untrained initialisation depends on the **budget**
> (significant at 8,192, absent at 32,768), on the **readout** (argmin versus diversity-respecting
> aggregation, which the geometry workstream found go opposite ways), and on **n** (three of this
> agent's own reads reversed between n = 1, n = 9 and n = 19). It **never** beats a classical greedy
> search at matched budget. And the selection gap that dominates every arm is **identical for VQE and
> control**, so the sampler is not what is being measured.

**The paper's lead claim therefore changes from "VQE loses" to "the VQE comparison is not robust to
reporting choices that the literature routinely leaves unstated, and the underlying quantity is
restraint error in both arms."** That is a harder claim to make and a better one, and it is the third
budget/convention trap this sprint has found, after the QNG equal-iterations-versus-equal-cost
inversion and the set-mean-on-a-collapsed-set trap.

### Q4.8 The agent refuted three of its own reads, in place

Its n = 1 read (generative objectives improve inside their tail) — refuted at n = 9. Its n = 9 read
(selective objectives have positive ordering terms) — refuted at n = 19. Its n = 9 budget-trap law
(ρ = −0.491) — collapsed to **−0.073** at n = 19. All three are preserved in place in its §5.

It also traced, rather than suppressed, the one tolerance it could not make exact: `all_ca` reproduces
`Enum.rmsd` to 3.86e-07 because the cached label is stored **float32** (spacing 2.4e-7 at 3 Å), and a
pure-float64 rebuild disagrees by the same amount. The tolerance was not weakened to make a check
pass.

---

## Q1. CVaR-VQE loses to its own untrained initialisation on every diversity-respecting readout

Relayed from the QGEOM workstream (`s15/qgeom_FINDINGS.md`). This is the sprint's quantum result at
the *set* level, and it was obtained in the most favourable setting a CVaR-VQE can be given: the
**exact statevector gradient**, no shot noise, no SPSA, with the budget charged only for readout
draws and the control given the same draws. 12-qubit sub-registers of the nine enumerated targets,
200 iterations, 3 seeds, 4,096 draws, paired over 27 cells per α.

**Positive numbers mean VQE is worse than best-of-N from its own untrained initial distribution at
the same budget.**

| readout | α = 1.0 | α = 0.25 | α = 0.05 | α = 0.01 |
|---|---|---|---|---|
| argmin (by objective) | −0.036 null | **+0.251 SIG** | +0.068 null | +0.048 null |
| top-20 mean | −0.177 null | +0.111 null | −0.070 null | −0.091 null |
| **top-20 coordinate average** | +0.068 null | **+0.372 SIG** | +0.197 null | +0.177 null |
| **top-75 coordinate average** | +0.147 null | **+0.392 SIG** | **+0.276 SIG** | +0.164 null |
| drawn-set MEAN | **−0.847 SIG** [22/5] | **−0.372 SIG** | **−0.191 SIG** | −0.066 null |
| drawn-set BEST | **+1.092 SIG** [1/25] | **+0.345 SIG** | **+0.155 SIG** | +0.087 null |

**Q1.1 — What concentration IS, measured directly.** At α = 1 the VQE's drawn set has a mean RMSD
**0.847 Å better** than the control's and a best member **1.092 Å worse** — two significant effects
of opposite sign and near-equal magnitude, both shrinking monotonically to zero as α falls.
**Concentration buys the mean and pays the best, and α scales both.** This is the cleanest statement
of the Sprint 14 mechanism the project has, and it is now measured rather than inferred.

**Q1.2 — On every readout that survives to a structure, VQE loses.** Both coordinate-average
readouts — which are what the pipeline actually uses — are significantly *worse* than the untrained
control at α = 0.25, and top-75 is worse at α = 0.05 as well. Nothing in the table is a significant
win for VQE on a structural readout.

**Q1.3 — A TRAP THAT INVALIDATES A LAW THIS PROJECT USES CONSTANTLY.** The standing law is that the
terminal operator consumes the set **mean** (`d_out = 1.16·d_set_mean + 0.04·d_set_best`, R² 0.89).
Read naively, the `drawn-set MEAN` row says CVaR-VQE wins by 0.85 Å. **The actual coordinate average
says it loses.** The resolution is diversity:

| α | distinct configurations among 4,096 VQE draws | control | final entropy (bits) | max p |
|---|---|---|---|---|
| 1.00 | **2.7** (median **2**) | 510.3 | 0.73 | 0.702 |
| 0.25 | 149.4 | 510.3 | 3.93 | 0.319 |
| 0.05 | 314.9 | 510.3 | 5.79 | 0.155 |
| 0.01 | 408.0 | 510.3 | 6.58 | 0.101 |

**At α = 1 the VQE returns two distinct structures out of 4,096 draws.** Its "top-20 set" is twenty
copies of one structure, so its set mean *is* that structure's RMSD and there is nothing left to
average.

> **The set-mean law was fitted on sets with real diversity and is OUT OF DOMAIN on a collapsed set.
> Set mean is a valid proxy for the terminal operator only at fixed set diversity; a method that
> concentrates improves the proxy and degrades the actual output.**

**This binds my own work and I am recording the check rather than assuming it.** `s15/augment.py` is
motivated by exactly this law — that mixing conformers better than the set mean into the retrieval
pool must lower the emitted structure. That argument is only valid because the augmenting set has
**genuine diversity**: eight independent multi-start fits, not eight copies of one solution. The
augment run must therefore report the **distinct-structure count and the spread of the fit
ensemble** alongside its RMSD numbers, or it is quoting a law outside its domain in precisely the way
this finding warns about. Same for the cascade's `A` stage, which averages 16 fits.

**Q1.4 — α is a diversity dial, not an accuracy dial.** The monotone entropy column (0.73 → 6.58
bits) is the cleanest reading of what α does. That is a legitimate and honest role for CVaR in
Family C — shaping ensemble spread — and it is *not* the role the CVaR-VQE literature claims for it.

## Q2. The variational geometry, corrected

Also from QGEOM, and all three items correct or upgrade a Sprint 14 record.

**Q2.1 — QNG and classical natural gradient are the same algorithm here.** For a real-amplitude
circuit measured in the computational basis, the classical Fisher information of the measured
distribution equals **four times** the Fubini–Study metric, to **3.3e-16**. The factor is absorbed
into the learning rate. This removes an entire arm of the design space rather than answering a
question about it.

**Q2.2 — Sprint 14's condition numbers were a seed-averaging artefact.** They are statistics of a
seed-averaged metric, not the condition number of the metric at any point. Recomputed per point the
spread is enormous: for a chain at L = 3, **median 735, maximum 1.5e6**, against a recorded 2.2–3.7.
The recorded rank column is the same artefact — the `block` ansatz at depth ≥ 2 is **exactly
rank-deficient**, which the seed average hides.

**Q2.3 — The gradient lies in the range of the metric, and avoids its small directions.** The
null-space share of the gradient is **0 to 1e-30** (a theorem, confirmed numerically), and the
bottom eigenvalue decile carries **0.0016–0.0043** of the gradient's squared norm where a uniform
gradient would carry 0.10. **The ill-conditioning is real and the gradient does not point into it.**
This upgrades a Sprint 14 hypothesis to a proof for the exactly-null directions and a measurement
for the rest — and it removes ill-conditioning as an explanation for the VQE's failure, which
strengthens Q1 by eliminating a competing account.

The agent also refuted its own mechanism for why QNG helps on one cost and not another. That
refutation is preserved in its findings file.

---

## Q3. The QGEOM workstream, completed: one positive result, six self-refutations, and a methodological trap about budget conventions

Relayed from `s15/qgeom_FINDINGS.md` (1,197 lines, 14 result JSONs). The instrument reproduced
exactly at the start **and** end of the workstream. Verification discipline throughout: exact
analytic derivatives against Sprint 14's finite differences (2.0e-11), classical Fisher against 4g
(3.3e-16), adjoint against derivative-state gradient (2.4e-15), exact Hessian against central
differences of the exact gradient (5.2e-11), CA rebuild against stored RMSD (1e-8).

### Q3.1 THE PROJECT'S FIRST POSITIVE QUANTUM RESULT — and it is flagged LOW POWER

**Consumed as an ensemble with no ranker anywhere** — no argmin, no top-k, no objective filtering,
just the drawn distribution aggregated — **VQE beats best-of-N by 0.36 to 0.57 Å**, with confidence
intervals excluding zero, at **α = 1 and α = 0.05**: exactly the settings where the distribution
retains diversity. The concentration check **PASSES**.

The agent's own verdict is **LOW POWER — replicate before building on it**, and that verdict is
adopted. This is now under replication at power by a dedicated workstream, with a matched-diversity
control, because two things are simultaneously true and only measurement separates them:

- **A mechanism exists.** The programme has measured repeatedly that *nothing ranks within the pool*,
  and in the generative cascade selecting the objective's argmin **costs +0.751 Å** while aggregation
  recovers only −0.360 Å. If every ranker is useless, a method whose value lies in the *shape of its
  distribution* rather than in any member it can identify would appear exactly here — invisible on
  argmin readouts, visible only when the ensemble is consumed unranked.
- **A reason to disbelieve exists.** At α = 1 the VQE returns **2.7 distinct structures out of 4,096
  draws**, and the win is non-monotone in α (present at 1 and 0.05, absent at 0.25). A non-monotone
  effect with no mechanism is usually noise.

**Until it replicates, this is not quoted as a result.** It is the only positive quantum claim the
project has ever had, which makes it the one most likely to be wrong.

### Q3.2 CVaR has no defensible role as an objective, tail-shaper, or constraint

A plain expectation penalty drives clash-violating probability mass to **exactly zero**; CVaR leaves
**2.8–5.6%** and returns worse structures. Both lose to best-of-N by **+0.41 to +0.70 Å**.

So of the three roles CVaR was given in the architecture families, two are dead and the third — the
diversity dial — is the only one the evidence supports.

### Q3.3 A TRAP: the budget convention decides the QNG conclusion

This is the most transferable methodological finding of the workstream.

| convention | result |
|---|---|
| **equal iterations** | QNG beats plain gradient descent **5/5** where the condition number > 5, **0/2** where < 5 |
| **equal hardware cost** | the full quantum geometric tensor wins **0 of 7** cells, with two clean sign flips |

> **"Reporting only convention A would have produced a false positive."**

Only the *diagonal* preconditioner survives cost accounting. And on the real VQE across 24 cells, QNG
is **indistinguishable from plain gradient descent** where the condition number is 688–2,677
(+0.0016, CI [−0.010, +0.014], 7 wins to 10) and **loses to plain Adam 17/17** (+0.067, CI excluding
zero), by up to 27× at equal cost.

**Why, mechanically:** `F = 4g` exactly, so no second metric exists to exploit; the gradient provably
lies in `range(g)` (confirmed to 1e-30) and *avoids* the small eigendirections (bottom decile
0.0016–0.0043 against 0.10 for a uniform gradient); so the rescaling QNG applies **along the
directions the gradient actually occupies** is only **1.4–3.1×** even at condition number 2,000.

This also refuted the agent's own pre-registered prediction — Sprint 14's depth-1 refutation of QNG
genuinely does not extend to greater depth under the equal-iteration convention.

### Q3.4 Target conditioning does NOT reshape the manifold in any information-bearing way

Preparing a state for the retrieval-conditioned torsion prior raises the Fubini–Study condition
number **3.85× to 36×**, cuts metric volume by 1.8–3.7 log units and effective-rank fraction by
4.1–7.5 points, at **27 wins to 0**. It looks like a large, clean effect.

Against an **entropy-matched scramble of its own prior** — per-residue entropies identical to the
bit, target information destroyed — it is **NULL on all seven statistics, on both ansätze, across 54
paired cells**, with the scramble reproducing the effect to within 3–7%. The fit-quality confound is
checked and null. The one statistic whose interval excluded zero (gradient magnitude) **fails the
mandated null-calibrated concentration check at p = 0.000–0.001**.

> **The chain `target information → manifold → metric → trainability` does not exist. The middle
> link is real and information-blind.**

### Q3.5 There is no barren plateau, and the scaling axis was misidentified

`Var ~ 2^(−0.18 to −0.36 n)`, with **every confidence interval excluding −1.0**. The real axes are
**depth** (`L^−2.6` at n = 16) and **cost locality** (6.05 log₂ at n = 14, against 2.8 log₂ for the
entire width axis) — refuting Sprint 14's "what matters is n" in its operational reading. The 1-local
retrieval prior has **zero** width decay (+0.006).

### Q3.6 Optimising the landscape well and getting a better structure are anti-correlated

Gradient norm falls **313×** while the mode's RMSD moves **0.049 Å**, and only **25%** of runs reach
a positive-definite Hessian. Most starkly: **the ansatz that optimises best reaches a minimum in 0 of
18 runs and *improves* structure by 0.46 Å, while the one that does reach a minimum makes structure
0.37 Å worse.**

That is the two-axis separation — optimisation quality versus structural accuracy — appearing in its
sharpest form yet, and it is the single cleanest illustration of why this programme reports both axes
and never substitutes one for the other.

### Q3.7 What was refuted, including six of the agent's own claims

Refuted: that QNG is neutral; the cost-profile mechanism for where QNG helps; `p(argmin) → α`; a
confounded locality cell; a meaningless control cell; and its own argmin claim, which **failed the
null-calibrated concentration check**. Also refuted: **Sprint 14's metric table** (an averaging
artefact — the condition numbers and ranks are those of a *seed-averaged* metric and do not reproduce
from the recorded code at either width; true per-point spread is far larger, `chain` at L3 having
median 735 and maximum 1.5e6, and `block` being *exactly* rank-deficient), **Sprint 14's locality
claim**, and **the paper's own thesis H2**.

`g_ii = 0.250000` survives.

### Q3.8 Two gaps, named rather than hidden

The agent complied with the resource throttle by serialising, and two cells are incomplete: the last
two rungs of the A3 ladder and A4b beyond 11 cells. **Direction is established in both**; the missing
cells are for completeness. It also verified mechanically that both of my corrections — the
`LogPTable` bin bug and the AMBER `snap_index` rule — are **inapplicable** to its work (`grep` for
`LogPTable|ShiftedLogP|SumLogP|distml` and for `amber|AMBER` across `s15/qgeom_*.py` returns
nothing), and modified no tracked file outside its own.

---

## I1. The information-channel audit: ERROR STRUCTURE, NOT ERROR SIZE, decides what a channel is worth

Relayed from the INFO workstream (`s15/info_FINDINGS.md`, figure at
`s15/figures/info_phase_diagram.png`, JSON in `s15/results/info_*.json`). The instrument selfcheck
reproduces exactly at the start and end of that run. This is the deepest result the sprint has
produced and it reframes the accuracy question a second time, after K5.

### I1.1 The phase diagram verdict, stated as a requirement

Required in-band pairwise accuracy to reach 2.0 Å through a top-100 operator: **0.638**, against a
chance level of 0.500. The best channel measured is **0.566**.

| target | verdict | what it requires |
|---|---|---|
| 3.0 Å | **supported** | coverage 1.0 at ≤ 26°, or 0.8 at ≤ 13°, or **0.5 with terminal gaps at ≤ 11°** |
| 2.5 Å | **marginal** | coverage ≥ 0.9 at ≤ 15°, or ≥ 0.6 with terminal gaps |
| 2.0 Å | **not supported** by any i.i.d.-class channel | 16.3° RMS at full coverage; the best real channel is 67.6° RMS |

Independent verifications of the diagram: a constant-helix control reproduces 4.065 Å exactly; the
incumbent's σ-equivalent measures 28.6° against a recorded 27.1°; the 2.0 Å requirement measures
16.3° against a recorded 15.1°.

### I1.2 The i.i.d. surface OVER-PRICES real channels — by 0.6–1.7 Å native-free, up to 2.9 Å ORACLE

This is the finding that matters most, and it invalidates the way this project has priced channels
for several sprints.

A real fragment channel at **64° RMS error emits 1.77 Å**. An i.i.d. channel of **the same error
magnitude emits 4.71 Å**. Error magnitude does not price a channel. What prices it is **where the
error points**.

The mechanism is measured, not inferred: **subspace alignment**, `||Je|| / (||e|| · s_rms)` =
**0.565** for the best fragment channel against a **0.945** random-direction null, with every
confidence interval excluding the null. Causality is established by surrogate destruction —
**sign-flipping the errors while preserving every `|error|` exactly costs +1.596 Å [+1.385,
+1.813]**. The magnitudes are untouched; only the directions change; the structure falls apart.

Two corollaries that close off natural but wrong responses:

- **None of this is visible in pairwise correlations** (all |r| ≤ 0.16). So Sprint 14's finding that
  real torsion predictors have near-i.i.d. errors is *confirmed and non-informative*: it measured
  the wrong statistic. Lag-1 correlation says nothing about Jacobian alignment.
- **Correlated is not aligned.** Synthetic coherence — AR(1) processes, class bias — moves the
  result the *opposite* way. The measured coherence factor is 1.24, not the 2 that was assumed.

### I1.3 Where the gaps fall beats how many there are

At σ = 0 and coverage 0.7, the emitted structure is **1.716 Å when the gaps are terminal** and
**3.505 Å when they are mid-chain** — a **1.79 Å spread**, three to four times larger than the
0.40–0.50 Å previously recorded.

**Terminal-gap-shaped coverage is the largest single lever in the phase diagram, and no channel
exploits it.** That observation is now the ALIGN workstream's task.

### I1.4 The channel table

126 real pools, in-band defined as `pool_best + 1.5`, ties averaged.

| channel | availability | in-band accuracy | CI95 | emitted @ 24 |
|---|---|---|---|---|
| **distogram Bayes risk, low-`sd` half** (new) | 126/126 | **0.566** | [0.547, 0.586] | 3.536 |
| distogram L2, low-`sd` quartile / half / 1/sd² | 126/126 | 0.560 / 0.559 / 0.548 | — | 3.65 / 3.55 / 3.52 |
| distogram Bayes risk (incumbent) | 126/126 | 0.545 | [0.522, 0.567] | **3.501** |
| ESM-2 contact map | 126/126 | 0.513 | [0.494, 0.531] | 4.384 |
| pool-consensus secondary structure | 126/126 | 0.513 | [0.498, 0.527] | 4.070 |
| radius of gyration | 126/126 | 0.497 | [0.475, 0.517] | 4.566 |
| pool typicality | 126/126 | 0.495 | [0.465, 0.526] | 3.723 |
| random null | — | 0.498 | [0.490, 0.507] | 4.492 |

### I1.5 The `sd` column is a real discriminator that the operator does not consume

Spearman(`sd`, |error|) = **+0.540**, and the RMS error rises **1.34 → 5.93 Å** across `sd`
quintiles. So the never-used uncertainty column carries genuine information about where the
predictor is wrong.

Conditioning the shipped score on it gives **+0.021 [+0.010, +0.034]** in-band accuracy, 74 wins to
51, concentration check PASS — **the best in-band number this project has ever measured** — and
**+0.033 [−0.051, +0.116] *worse* emitted RMSD.**

A channel can be genuinely informative and still not survive the operator that consumes it. **This
does not contradict K1**: `sd` weighting is worth 0.154 Å in the *generative* fit (3.644 vs 3.798
unweighted). Selection and generation are different axes and the same column pays off on one and
not the other — which is now the third independent instance of that pattern in this sprint, after
K6's pool channel and N5's debiasing.

### I1.6 The distogram is a good ranker and a bad probability

Held-out predictive cross-entropy against baseline: **−0.818 bits/pair [−1.265, −0.419], positive
on only 52 of 126 targets.** ESM-2 contacts are far worse at **−2.570 bits/pair, positive on 0/126**.

This bears directly on the maximum-likelihood work and supplies a mechanism independent of binning.
ML consumes the distribution *as a probability*, so the whole density shape enters the objective;
least squares consumes only its first two moments. If the shape is miscalibrated on 74 of 126
targets, ML is faithfully fitting a badly wrong density while least squares discards precisely the
part that is wrong. **Prediction: a continuous density will recover the quantisation cost and ML
will still lose**, because de-discretising an uncalibrated density does not calibrate it. Relayed to
the RESTRAINT agent with a proposed test — do the targets where ML wins coincide with the 52 where
cross-entropy is positive?

### I1.7 The retrieval torsion prior carries no target-specific bits

**−0.11 to −0.23 bits/torsion** against a generic Ramachandran distribution at every bin count. It
*is* target-specific (a target-shuffled null sits at −0.70 to −0.80) — just not *better* than the
generic prior. Two independent cross-checks agree: modal hit rate 0.371 versus 0.360, and the
structural twin gives 4.072 Å against a **4.065 Å constant α-helix**.

The plug-in mutual-information estimator says **+0.40 to +0.88 bits on the same data** — a measured
demonstration, on our own data, of exactly the estimator trap the brief warns about. Use held-out
predictive cross-entropy; never plug-in MI.

### I1.8 A methodological result: the skill-definition discrepancy was five things, not one

The long-standing disagreement between a +0.416 and a +0.909 skill figure is resolved. The **+0.416
replicates exactly**. The **+0.909 is a near-tautology on a structural axis**: ρ(Rg, RMSD) against
native-Rg z is **−0.935 (p = 0.000)**, which is definitional rather than empirical. On the actual
leave-fold-out learned model the same quantity is **+0.154, p = 0.531**.

Two facts were hidden inside the disagreement. The two skill statistics **anti-correlate**
(ρ_global +0.154 versus in-band accuracy −0.409/−0.489, p = 0.032), so "the sign" was never one
quantity. And against **emitted RMSD**, every compactness variant gives |ρ| ≤ 0.12 at p ≥ 0.63 — the
chain breaks at the last link.

### I1.9 The INFO agent refuted its own headline

Its native-free per-target sign proxy picks the sign at **0.947** on the enumerated band and
captures **99.7%** of an ORACLE channel worth −0.339 Å [−0.564, −0.159] at 19 wins to 0. **It does
not transfer.** On the 126 real pools the entire ORACLE channel is worth only 0.038–0.083 Å, the
proxy recovers 0.020 Å [−0.051, +0.010] at permutation p = 0.067, and **outside the band it is
harmful, +0.43 to +0.68 Å**. Tiered REFUTED as a route to 2.5 Å.

This is the decoy-bank-is-not-a-pool-proxy law applying to its own author's headline result, which
is the standard the programme is trying to hold.

### I1.10 Whole windows beat recombination on the mean and lose on the best

Whole retrieved windows beat independent per-residue recombination by **−0.457 Å on the mean** and
**lose by +0.245 Å on best-of-500** (2.048 versus 2.293). Recombination raises the ceiling and
lowers the typical member, so **it helps only if the selector is sharp** — and nothing in this
project ranks within the pool. A direction to revisit only if selection is ever solved.

---

## A1. Phase 0 audit — the gate PASSES, and four things change how everything else is reported

Full record in `s15/audit_FINDINGS.md`; the binding consequences are now written into `s15/BRIEF.md`
so every agent inherits them. The four that change *this* file's tables:

1. **`I.FAIL18` is a threshold artefact.** `BAND = 1.5 Å` has no derivation anywhere, and
   `n_zero_recall` runs 45 → 2 across BAND 0.5 → 3.0. On Sprint 12's own 99-combination sweep,
   **exactly one of the 18 targets (9KAR) survives in every version of the set**; mean Jaccard 0.498.
   Three sprints have treated it as an object. The column stays in every table — it is comparable to
   every prior sprint and dropping it would break continuity — but no argument may rest on set
   membership, and it must never be described as "the hard targets" without this caveat.
2. **`amber_kind == 0` is unbiased on the mean but not in the tails**, because `a_idx` force-includes
   the ORACLE snap index, which is the minimum-RMSD member on 6/19 files and inside the <1.5 Å
   in-band set on 14/19. Binding rule for every tail, in-band, argmin or top-k statistic:
   `amber_kind == 0 AND amber_idx != snap_index`, with per-target n printed.
3. **The brief's oracle-conditioning share was transposed** — 21.7% oracle-conditioned, 40%
   unbiased, not the reverse. The −0.401 Å effect size was right.
4. **Every AMBER row in `s14/results/obj_floor.json` is retracted** (computed with no `amber_kind`
   mask, yet flagged complete). Corrected values live in `s14/obj_FINDINGS.md` §1.3.

Two further items matter for the paper rather than for the tables. **RMSD is now frozen and
independently reimplemented** — Horn's quaternion method agrees with `kabsch_rmsd_batch` to
2.04e-13 Å over 63,000 real structures and passes the mirror test the naive SVD fails; the frozen
definition is full-chain CA including both termini, proper rotations only. Full-chain versus
`[1:-1]` differs by up to **2.4 Å** on one structure, so the choice is not cosmetic and must be
stated. And **AMBER costs 8–14 ms on these targets, not the 28 ms the brief claimed**, so Sprint
14's budget matching over-charged AMBER by 2–3× — which means AMBER was given *less* search than a
fair match, and its negative result is if anything understated.

---

## A2. The multi-start draw was not reproducible across processes. Fixed, with the scope of the damage stated

**DEMONSTRATED, and it is a defect in code I wrote.** Found by the RESTRAINT agent while checking
that its representation ladder shared start points with the coordinator's reference arm — not part
of its brief, and reported anyway, which is exactly the behaviour the programme needs.

Every multi-start module in this sprint seeded its RNG with

```python
rng = np.random.default_rng(hash(pdb) % (2 ** 32))
```

**Python's string hash is salted per process** unless `PYTHONHASHSEED` is set, and the sprint's own
environment lock (`s15/results/audit_env_lock.json`) records it as `null`. Verified directly: two
consecutive interpreters return `hash('1A13') = 217586290314588545` and `2408026022170661001`. So
every run drew a **different set of multi-start initialisations**, and no multi-start number in this
sprint was bit-reproducible on re-running.

**The scope of the damage, stated precisely, because it decides which results survive.**

- **Within one process the defect is harmless.** Every arm in a single run shares one start set, so
  every paired comparison inside a single result file is start-matched and valid. All the
  within-file conclusions of this sprint stand unaltered — K1-CORRECTED's arm ordering, K6's arm
  table, the restraint ladder's decomposition.
- **Across processes it is not harmless.** A figure from one module compared against a figure from
  another compares different start draws. And the published constants — the ORACLE 0.611 Å and the
  predicted 3.644 Å among them — are **not bit-reproducible on re-running**. They are correct
  measurements with a run-to-run start-draw variance that was never quantified and never declared.

**Fixed.** `s15/seed.py` supplies `stable_seed`/`stable_rng` built on `blake2b`, which is stable
across processes, platforms and Python versions. Seven modules — `cascade`, `distacc`, `distcal`,
`distgeo`, `distml`, `pooldist`, `robust` — are patched to use it, all import cleanly, and the seed
is verified identical across two separate interpreters (3025285571 both times).

**What is deliberately NOT being re-run right now, and why.** The in-flight cascade was launched
before the patch and will finish under the old seeding. It is left to finish. Its purpose is to
answer a scientific question — whether restraint-driven generation with aggregation clears the
incumbent — and a within-run paired comparison answers that correctly. **The frozen protocol will
be re-run under stable seeding before the benchmark is touched**, and only those numbers will be
quoted as reproducible constants. Exploratory runs need not be bit-reproducible; a frozen protocol
must be. That distinction is recorded here so it cannot later be presented as a convenience.

**The methodological point worth carrying into the paper.** This defect was invisible to every
check the programme already ran. Determinism was verified — repeatedly, and at three thread counts
— but always *within* a process or against a *cached* artefact, so a per-process salt could never
show up. A reproducibility audit that never starts a second interpreter cannot detect the most
common source of irreproducibility in Python.

---

## V1. THE ADVERSARIAL CONSISTENCY AUDIT: 13 blockers, all applied

Relayed from the VERIFY workstream (`s15/verify_FINDINGS.md`, 744 lines). It traced **178 distinct
quantitative claims**: **141 matched**, **21 mismatched**, **16 untraceable** (each itemised with
why). Every blocker below has been fixed in the documents; this section records what was wrong, so
the corrections are part of the record rather than silently absorbed.

**What reconciled exactly**, so the absence of a finding is informative: the entire cascade table,
the ceiling table, the ranking table, the feasibility table, the phase diagram, the realizability
table, the n = 126 fusion law, the projection gap, and all 20 CVaR readout cells. The 25.8%
withdrawal is complete across all 24 documents. And there is **no second instance** of the
per-target-minimum mislabelling anywhere in `s15/*.py`.

### The thirteen

| # | what was wrong | fix applied |
|---|---|---|
| **B1** | the benchmark is untouched — but `audit_provenance.py` declares its manifest "not touched here" and then **hashes its bytes**. No leak (digest, size, mtimes only); the *comment* was false | comment corrected to say the bytes are read to hash it |
| **B2–B4** | n = 24 fusion numbers (3.270 / 3.469 / **2.980** / **+0.627**) standing **with no n**, sixteen lines before the n = 126 table in the same document; the same stale pair in the dossier and the paper, whose own abstract already gave the right values; and **−0.513** sitting inside n = 126 tables although `coherence.py` computes no such field | all replaced with n = 126 values; −0.513 moved out of the table and labelled n = 24, not yet recomputed |
| **B5** | RESULTS §5.2's surrogate table reproduced at **no** prefix of its artefact, and its **only native-free row flipped sign** (+0.179 worse in the document, −0.065 better in the artefact) | table rewritten from the artefact at n = 10; the derived claim restricted to the **ranking** axis, where it is directly measured |
| **B6** | the decorrelation ladder did not reproduce, and **its most-quoted sentence was no longer true** — the "geometrically impossible" dip below both endpoints does not exist in the data | §5.1 rewritten at n = 40 (the requirement is **69%**, not 30%; the gain is 1.284 Å, not 0.737; the response is **threshold-like, not linear**), and the spurious paragraph **withdrawn** |
| **B7** | "**2,048 draws**" wrong in 20+ places — the experiment drew **4,096** | corrected everywhere, including the figure axis and `relayed.json` |
| **B8** | "+0.127, z = 13.1" described as the **global** ρ; it is the **in-tail** ρ. Globally the classes are close (+0.523 against +0.485), so "generative objectives have better global ρ" was **false as stated** | corrected in three documents |
| **B9** | the torsion-channel headline (0.565, 1.77 Å, +1.596 Å, and the top of "0.6–2.9 Å") is **ORACLE** and was unlabelled; the honest native-free range is **0.6–1.7 Å**; and the docs omitted the control that matters — **a constant α-helix already scores 0.709** against the 0.945 null | labelled ORACLE, native-free range stated, helix control added, and the sign-flip result noted as **null for the native-free arm** |
| **B10** | "**five times the headroom**" — in the paper's abstract — **is not derivable from any number in the sprint**. The supported ratios are **3.19×** (floor) and **2.07×** (headroom) | corrected in eight documents |
| **B11** | the dossier stated a quantum positive the replication workstream had **already demoted** | rewritten as a **demotion**; see Q5 |
| **B12** | the `n ≥ 100` figure guard was **bypassed three ways**, one of them a headline figure, and one relayed panel drew **9 targets** with no n printed | guard now finds counts under any of three key names, `_relayed()` routed through it, and **all three remaining exceptions are declared and print their own n on the panel** |
| **B13** | two captions contradicted their own data — fig16 said "every stage after it loses ground" over an arrow drawn as a **−0.360 gain**, and fig05 called the bias "a pure additive defect" where the sprint measures correcting it as worth almost nothing | both rewritten |

### Q5. The quantum positive is withdrawn

The audit's B11 is the most consequential single item, and it is recorded here rather than only in
the table because it removes a claim the dossier led with.

| the original reading | what the audit established |
|---|---|
| "27 cells, intervals excluding zero" | the 27 cells are **9 targets × 3 seeds**; at the **target level (n = 9) the α = 1 cells are NULL** — −0.567 [−1.234, **+0.085**] and −0.437 [−1.056, **+0.160**] |
| "concentration check PASSED" | at the *cell* level only; at the target level the surviving α = 0.05 cell **FAILS**, carried by two or three of nine targets |
| matched budget | the VQE was charged **819,200** evaluations against the control's **2,048** — a **400× advantage**, and that is the convention the positive was measured under |
| the readout | `np.unique` returns *sorted* indices, so the subsample was the lowest-indexed configurations rather than a random draw, biting unequally across arms |

**The programme has no positive quantum result.** And my stated "reason to disbelieve" — non-monotone
in α, 2.7 distinct structures — was **not** the operative one. The operative ones are the unit of
analysis, the concentration failure, and a 400× budget asymmetry.

### What this audit says about the rest of the sprint

Three of the thirteen (**B5, B6, and the numbers behind B11**) share one root cause: **quoting
results produced before the seeding fix**. `PROTOCOL_FROZEN` §4 already forbade that, and I did it
anyway, three times. Its reassurance that such numbers are "correct but not bit-reproducible" is
**too weak a statement**: re-running flipped a sign in B5 and destroyed a rhetorical claim in B6.

Two more (**B9, B10**) are of a different and more embarrassing kind: a number quoted so often that
its derivation was never re-checked. "Five times the headroom" appeared in eight documents including
the paper's abstract, and the sentence that introduces it contains arithmetic yielding 2.07.

**The pattern across all thirteen is that none of them was found by thinking harder about the
science.** They were found by loading the artefacts and comparing. That is the third time this sprint
that mechanical checking beat intuition, after the non-uniform bin lookup and the per-process hash.

---

## SCOPE DECISIONS, recorded rather than left implicit

**`s15/robust.py`'s full-instrument run was stopped, because it had already been answered.** The
ALIGN workstream measured the same arms at n = 126 inside `align_fit` — `CAUCHY_lfo`
−0.011 [−0.102, +0.077], `gnc_cauchy` −0.025 — firing `robust.py`'s own pre-declared falsification
condition. Continuing my run would have re-measured a refuted hypothesis on a machine at 99% CPU.
**What is lost:** my module's specific arm set (`huber`, `soft_l1`, `welsch`, `trimmed`, and the
graduated-non-convexity schedule) at n = 126; the two redescending families were covered by ALIGN,
the rest were not. **What is gained:** capacity for `decorr` and `augment`, which test live
hypotheses. `s15/robust.py` is committed and runnable, and the 8-target smoke read in
`s15/results/robust.json` is labelled as such wherever it appears — the figure module refuses it.


The brief requires that anything scaled down be stated as such, so:

**`s15/scale.py`'s full-instrument run was stopped, deliberately.** Its ORACLE arms had already closed
the direction at n = 8: the ceiling of *any* per-target isotropic rescaling of the restraints is
**−0.090 [−0.284, +0.033]**, i.e. about 0.1 Å, and every pool-referenced estimator was on the wrong
side of zero with the pool's scale estimate correlating **+0.035** with the truth. K11's independent
full-instrument measurement of the same idea on the *consensus* (rather than the restraints) reached
the same verdict at n = 126, with all four native-free estimators uncorrelated with the truth
(|r| ≤ 0.106). Running the restraint version to 126 would confirm a direction that two measurements
already close, and it was consuming a slot on a machine pinned at 99% CPU with 1.9 GB free.

**What is lost by that decision:** a full-instrument confidence interval on an arm whose ORACLE
ceiling is below the sprint's measured false-positive floor of ~0.08 Å. **What is gained:** capacity
for `robust`, `decorr` and `augment`, which test live hypotheses rather than confirming a closed one.

`s15/scale.py` is committed and runnable; the n = 8 result is in `s15/results/scale.json` and is
labelled as a smoke read wherever it appears. The figure module **refuses** to render it.

---

## WHERE THE SPRINT ENDED

**The architecture failed and the protocol is frozen as the incumbent, unchanged.** Thirteen
candidates were tested against their own controls on the full 126-target instrument. Twelve failed
outright; two cleared *their own* control while improving a losing arm without reaching the
incumbent. By a rule fixed in advance, the **60-target benchmark is left unspent and sealed** —
a confirmatory run on an unmodified pipeline could only return the pre-registered `null`.

**What the sprint established instead**, in the order that makes it intelligible:

1. **The ceiling was the library, not the channel** — perfect distance knowledge is worth ≈0.6 Å
   (sd 0.13 across start draws), not ≈1.95 Å. K1-CORRECTED, K5.
2. **Error SHAPE beats error MAGNITUDE, in both directions, by more than an ångström.** Outlier-shaped
   error at 4.12 Å RMS emits 1.993 Å; i.i.d. at 3.00 Å emits 2.561 Å; the real predictor at 3.70 Å
   emits 3.644 Å. K12.
3. **The mechanism: the errors are geometrically REALIZABLE** — 2.4× closer to realizable than matched
   noise, paired −1.391 [−1.623, −1.184]. The prediction describes a coherent *wrong structure* and
   the optimiser builds it. K9, confirmed at n = 126.
4. ~~**A parameter-free, native-free law for fusion**, predicting to 0.162 Å across 126 targets. K10.~~
   **⚠ WITHDRAWN, Sprint 16, 2026-09-06 (RETRACT).** The law is the Krogh–Vedelsby (1995) ambiguity
   decomposition (prior art); the residual with the correct quadratic mean is +0.075 Å, not 0.162 Å,
   and even that is a frame convention rather than a prediction error — the identity is exact to
   2.65e−15 on all 126 targets. The native-free screen it was supposed to enable has AUC **0.401**
   against a 0.500 null. See §"The fusion law" above and `s16/retract_FINDINGS.md`.
5. **Alignment is a description, not a lever.** The ORACLE rotation is worth −1.822 Å; the best
   native-free arm of 22 is −0.007. AL1. Family D stays empty.
6. **Even perfect decorrelation does not reach 2.0 Å** — it stops at 2.320 Å; 2.5 Å needs 86% of the
   coherence destroyed. K13-FINAL.
7. **Every quantum effect dissolves under a classical control** using the same circuit's samples, the
   same objective, and 1/200th the budget. Q6.
8. **Per-target confidence works native-free** (ρ +0.665; best quartile 2.354 Å against worst 4.884) —
   and the two signals I proposed for it add **nothing** over what the pipeline already emits. K14.

**Eleven retractions and corrections by the coordinator**, listed in `LEDGER.md`. Two were caught by
checking machinery against its own assumptions; seven by an adversarial audit that loaded the
artefacts and compared them to the prose. **None was caught by thinking harder about the science.**

**The one genuinely open lead.** A native-free error-direction surrogate exists — `θ_fit − θ_pool`
has |cos| 0.390 with the true error against a 0.157 null — and the structure's quiet subspace is
separately identifiable native-free (subspace overlap cos² 0.827 against a 0.499 null). Both
ingredients an alignment intervention needs are estimable without the native; putting them together
as a steering signal is untested. AL1.6.

**Still in flight at the time of writing:** the replication of the AMBER k = 30 relaxation, the
programme's one positive physics result, which sits **below** the measured ~0.08 Å false-positive
floor and is quoted only with that caveat until it lands.
