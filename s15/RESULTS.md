# RESULTS

Every result of Sprint 15, in the order that makes them intelligible rather than the order they were
obtained. ORACLE arms are labelled in every row and never appear as predictive results. Sections
marked *(in flight)* have a full-instrument run in progress and report the smoke read with its n
stated inline.

**The headline, stated first and without softening: the architecture this sprint was built to test
does not beat the incumbent.** What it produced instead is a mechanism that explains why, and that
mechanism is the contribution.

---

## 1. THE HEADLINE TABLE

The generative cascade, all 126 targets, combined restraint channel, leave-fold-out separation
debiasing, 16 native-free starts, clash gate at 2% of pairs below 3.6 Å (2 structures excluded across
the whole instrument).

| stage | mean | median | <2 Å | <2.5 Å | FAIL18† | vs incumbent (paired) |
|---|---|---|---|---|---|---|
| **G** best in ensemble — **ORACLE** | 2.760 | 2.547 | 0.35 | 0.49 | 5.155 | **−0.444 [−0.573, −0.325]** |
| **S** objective's argmin | 3.511 | 3.327 | 0.24 | 0.29 | 6.217 | **+0.307 [+0.178, +0.439]** |
| **A** coordinate consensus | 3.151 | 2.946 | 0.27 | 0.37 | 5.838 | −0.053 [−0.125, **+0.020**] |
| **A** over the objective-best half | 3.206 | 2.962 | 0.27 | 0.38 | 5.957 | +0.002 [−0.086, +0.088] |
| **F** after projection — *the like-for-like number* | 3.321 | 3.156 | 0.26 | 0.35 | 6.042 | **+0.117 [+0.050, +0.188]** |

gaps: **G→S +0.751** · **S→A −0.360** · **A→F +0.170**

† `I.FAIL18` is a threshold artefact — `BAND = 1.5 Å` has no derivation, the zero-recall count runs
45 → 2 across BAND 0.5 → 3.0, and only 1 of the 18 targets survives every threshold choice. Reported
for cross-sprint comparability; no argument rests on set membership.

**Reading it.** Only **F** is like-for-like: the incumbent's 3.204 Å is measured after the same
projection, and `A` at 3.151 is an unprojected coordinate average that is not a physically valid
structure. F **loses by +0.117 [+0.050, +0.188]**; even A's interval includes zero.

**But the ensemble is good and the readout is bad.** G = 2.760 beats the incumbent by −0.444 with an
interval nowhere near zero: structures better than what the incumbent emits are generated on most
targets, with 49% of targets holding a member under 2.5 Å and 35% under 2.0 Å. Every ångström of the
loss is in reading that ensemble out.

**Selection is actively harmful.** G→S costs +0.751 Å. Aggregation recovers −0.360 Å of it. And
pre-filtering the ensemble by the objective makes aggregation *worse* (3.206 against 3.151) — the
objective's within-ensemble ordering is worse than useless, which independently reproduces the
in-band blindness of §3 and closes the obvious "just filter the ensemble" response.

**The collapse check.** The set-mean law is out of domain on a collapsed ensemble (§6.3), so it was
verified: the fit ensemble's mean pairwise RMSD is **0.527 Å** (median 0.477) with only 17 of 126
targets below a 0.25 Å spread. The aggregation result is inside the law's domain.

---

## 2. THE CEILING, AND WHAT IT REFRAMES

Fitting continuous torsions to the **true** distance matrix through the identical machinery:

| arm | mean | median | <2 Å | vs incumbent |
|---|---|---|---|---|
| **ORACLE**, true distances | **0.611**† | **0.055** | **0.86** | **−2.593 [−2.901, −2.295]** |
| predicted, 1/sd² weighted | 3.644 | 3.466 | 0.16 | **+0.440 [+0.290, +0.592]** |
| predicted, unweighted | 3.798 | 3.625 | 0.12 | +0.594 [+0.430, +0.761] |
| the retrieval circular-mean start it was handed | 4.072 | — | — | +0.868 |

A median of **0.055 Å** means that on half these targets ideal-geometry torsions reproduce the native
Cα trace from its own distance matrix to a rounding error. **The representation, the ideal-geometry
constraint and the search are not the barrier**, and because both arms share optimiser, starts,
selection rule and manifold, the 3.03 Å between them is **entirely restraint error**.

**The reframing.** Two standing project results said perfect distance knowledge caps this problem at
≈1.95–2.0 Å. Both were measured with the retrieval library in the loop — that is the ceiling of
*selection through a finite library*, not of the distance channel.

> **Perfect distance knowledge is worth 0.611 Å, not 1.95 Å. Every prior estimate of the value of
> better distance prediction was computed against the wrong floor and understated the floor by 3.19× and the addressable headroom by 2.07×.**

The 0.611 Å is an ORACLE ceiling and never appears without the 3.644 Å beside it. The headroom is
real; the means to reach it is not demonstrated.

**Two supporting measurements.** The distogram's `sd` — unused in fourteen sprints — is worth
**0.154 Å** in the fit (3.644 against 3.798). And the fit improves on the start it was handed by
**0.428 Å**, so the optimiser does real work rather than returning its input.

---

## 3. DOES THE OBJECTIVE RANK THE TRUTH? THE REBUTTAL TO ROGET ET AL.

No search involved: evaluate every candidate objective at the native and at all 500 retrieved
windows, 126 targets. `native pct` is the native's rank among pool members **on the objective**;
0.5 is blind, above 0.5 actively disfavours the native.

| arm | native pct | median | is argmin | argmin RMSD | vs random | wins | ρ(f, RMSD) |
|---|---|---|---|---|---|---|---|
| `ls_pred` distogram, 1/sd² | 0.347 | 0.312 | 0.032 | 3.504 | **−0.949** | 0.79 | **+0.562** |
| `ls_debias` + separation debias | 0.359 | 0.324 | 0.032 | 3.559 | −0.895 | 0.76 | +0.557 |
| `ml_pred` full 17-bin likelihood | 0.394 | 0.330 | 0.016 | 3.450 | −1.003 | 0.83 | +0.532 |
| `ls_pool` pool median distances | 0.564 | 0.577 | 0.000 | 3.757 | −0.696 | 0.71 | +0.401 |
| `ml_pool` pool histogram | **0.665** | 0.722 | 0.000 | 4.026 | −0.427 | 0.56 | +0.365 |
| `ml_comb` distogram × pool | 0.463 | 0.434 | 0.000 | **3.421** | **−1.032** | **0.85** | +0.517 |
| **`ORACLE_true`** | **0.000** | 0.000 | **1.000** | 2.006 | −2.447 | 1.00 | +0.862 |

**Machinery verified**: `ORACLE_true` places the native at percentile **0.000 on 126/126**, is its
own argmin on 126/126, ρ = +0.862 positive on 100% of targets.

**The objective's argmin beats a random pool member by −0.949 Å [−1.147, −0.753]**, winning on 79% of
targets with ρ positive on **88.1%**. Fusion does slightly better at −1.032 [−1.213, −0.859] on 85%.

Roget et al. (arXiv:2606.21241) enumerate every peptide up to length 15 across 12,446 PDB structures
and find the minimum-cost conformation has, on average, a **larger** RMSD than a random feasible one —
for a Miyazawa–Jernigan lattice contact potential, whose parameterisation on proteins above 50
residues they identify as the cause. **At the same chain lengths our restraint objective's argmin is
about an ångström better than random.** Their pathology is a property of contact potentials, not of
short-chain structure prediction.

**But the native itself sits at the 34.7th percentile and is the argmin on only 4 of 126 targets.**
The objective orders the pool broadly well and cannot find the truth. This independently reproduces a
standing result (36.8th percentile, argmin on 3/126) through completely different code.

**A better estimator is a worse objective.** The pool channel has lower MAE than the distogram
(2.249 against 2.386) with **no training**, and is the worst arm here — native at the 66.5th
percentile, *above* chance. A new instance of the standing law that MAE does not price selected RMSD.

---

## 4. FAMILY B: FALSIFIED AS DECLARED, REDIRECTED ON EVIDENCE

Family B was `min E_AMBER subject to C_restraint ≤ ε`. Its precondition is that the native satisfies
the predicted bands. Measuring `z = |d_native − d̂| / sd` over 8,549 pairs:

| ε (units of the predictor's own sd) | fraction of pairs admitting the native |
|---|---|
| 0.5 | 0.248 |
| 1.0 | **0.458** |
| 2.0 | 0.712 |
| 3.0 | 0.840 |
| median z | **1.336** |
| max z | **7.237** |

At one sigma **fewer than half** the restraints admit the native; at three sigma one pair in six
still excludes it. A uniform-ε hard constraint has a feasible set that excludes the answer at every ε
that meaningfully constrains.

**The failure is concentrated, not diffuse** — a typical restraint is only mildly inconsistent and a
minority are catastrophically wrong — which is an outlier problem and redirects the family to robust
potentials rather than uniform bands. *(the full robust-loss ladder is in flight; an 8-target smoke
read shows the redescending signature: `welsch` gives a better median 2.466 against 2.777 and three
times as many targets under 2 Å, with a worse mean 2.967 against 2.704)*

---

## 5. THE MECHANISM: ERRORS THAT DESCRIBE A WRONG STRUCTURE

### 5.1 Real distance errors cost more than i.i.d. noise of the same size

Corrupting **true** distances with i.i.d. Gaussian noise and refitting (40 targets):

| σ (Å) | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|
| emitted RMSD | 0.554 | 0.709 | 0.939 | 1.306 | 1.450 | 1.787 | 2.075 | **2.561** |

The distogram's MAE of 2.386 Å is a Gaussian σ ≈ 2.99, where the surface says 2.56 Å. The real
distogram emits **3.644 Å** — about **1.1 Å more expensive than noise of its own magnitude**.

### 5.2 Surrogate destruction: the sign pattern carries the whole effect

*(n = 10, the current state of a full-instrument run; supersedes an n = 6 read whose numbers were
produced **before** the seeding fix and are not reproducible. Every surrogate is an **ORACLE
DIAGNOSTIC** — destroying a property of the error requires knowing the error.)*

| surrogate | mean | vs `real` |
|---|---|---|
| `real` (native-free) | 2.960 | — |
| **`ORACLE_shuffle_signs`** magnitudes kept, signs randomised | **1.842** | **−1.118** |
| `ORACLE_gauss_matched` i.i.d. at the target's own RMS | 2.401 | −0.559 |
| `ORACLE_shuffle_pairs` permuted across pairs | 2.448 | −0.512 |
| `ORACLE_winsor_z3` outlier tail clipped | 2.848 | −0.112 |
| `debias_sep` (native-free) | 2.895 | −0.065 |
| `ORACLE_true` — the ceiling | 0.591 | −2.369 |

**Randomising only the directions, keeping every magnitude, is worth more than an ångström** — far
more than clipping outliers (0.112) or matching a Gaussian (0.559).

**A correction to an earlier version of this table.** It previously showed `debias_sep` at **+0.179
worse**; the current artefact has it **−0.065 better**, and it is better at every n checked. The
earlier numbers predate the seeding fix and are not reproducible, which is exactly what
`PROTOCOL_FROZEN` §4 forbids quoting. The claim that separation debiasing is "worth little and
sometimes negative" is therefore **restricted to the ranking axis**, where it is directly measured
(`ls_debias` native percentile 0.359 against `ls_pred`'s 0.347, argmin RMSD 3.559 against 3.504). On
the **fit** axis it is worth little and mildly positive at this n — the honest statement is that it is
**at the noise floor on both axes**.

### 5.3 The coherent component is not a scale, and the pool cannot estimate it

*(n = 8 smoke read; full run in flight)*

| arm | vs real |
|---|---|
| `ORACLE_scale` — the best possible isotropic rescale | **−0.090 [−0.284, +0.033]** |
| `ORACLE_affine` — the best possible affine map | −0.116 [−0.448, +0.194] |
| `pool_scale` / `pool_shift` / `pool_affine` / `pool_scale_rg` | all **+0.02 to +0.15**, i.e. worse |

**The ceiling of any per-target rescaling is ~0.1 Å**, and the pool's estimate of the true scale
correlates **+0.035** with it. The direction is closed by its own ceiling.

### 5.4 The answer: the errors are geometrically REALIZABLE

Fitting to `d̂` and comparing the converged fit's own residual against the prediction's true error
(n = 24):

*(now on all 126 targets)*

| distance set handed to the fit | unrealizable part |
|---|---|
| **the REAL predicted distances** | **0.977 Å** |
| the retrieval pool's distances | **0.458 Å** |
| **ORACLE i.i.d. Gaussian of the same RMS — the null** | **2.369 Å** |
| ORACLE the same errors permuted across pairs | 2.081 Å |
| ORACLE the same magnitudes with random signs | 2.070 Å |
| *(mean \|error\| for reference)* | *2.339 Å* |
| **ratio real/null** | **0.413**, paired **−1.391 [−1.623, −1.184]** |

> **The real prediction is 2.3× closer to realizable than matched-magnitude noise. The predicted
> distance matrix is not a noisy version of the right structure — it is a good description of a WRONG
> structure, and the fit reproduces it faithfully.**

Both surrogates that destroy the structure move the residual back to within 13% of the i.i.d. null,
so the low residual is a property of the prediction, not of the fit's flexibility.

### 5.5 Why this unifies every negative

Every intervention attempted operates on the error's **magnitude profile**; a realizable error is a
**direction**, and lies in the null space of all of them. Five independent confirmations:
recalibration is provably worth zero (a weighted least-squares argmin is invariant under uniform
rescaling of the weights); separation debiasing is worth little (and negative on the ranking axis); outlier
clipping is worth 0.12 Å; per-target rescaling has a 0.1 Å ceiling; consensus rescaling has a 0.275 Å
ceiling the distogram overshoots (§7).

---

## 6. FUSION: A PARAMETER-FREE LAW

### 6.1 Correlated errors, different wrong structures

| | n = 126 |
|---|---|
| distogram fit vs native | 3.622 Å |
| pool fit vs native | 3.659 Å |
| **the two fits disagree with each other by** | **3.138 Å** (median 3.420) |
| raw per-pair error correlation | **+0.602** |
| the pool's unrealizable part | 0.458 Å against the distogram's 0.977 — *more* realizable |

Error correlation has been the wrong statistic. What matters is whether two channels' realizable
components point at the **same wrong structure**, and here they do not.

### 6.2 ~~The law~~ — the two-member ambiguity decomposition (SUPERSEDED, Sprint 16)

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT workstream.** Three corrections; the original
> table is preserved beneath with the affected rows struck through.
>
> * **Prior art.** This is the two-member **Krogh–Vedelsby (1995) ambiguity decomposition**
>   (`s16/lit_FINDINGS.md` §B.5). Cite; do not claim.
> * **Wrong mean.** `r` was the **arithmetic** mean of the two channels' RMSDs; the identity needs
>   the **quadratic** mean `√((r₁²+r₂²)/2)`. The channels differ by mean 0.980 Å per target
>   (median 0.567, max 4.798). Corrected: prediction **3.292**, error **+0.075** (median +0.019);
>   paired **−0.087 [−0.119, −0.061]** i.i.d.-target, **[−0.129, −0.056]** fold-clustered,
>   **W/L 126/0**, n = 126, TARGET as the unit.
> * **The residual is a frame convention, not a prediction error.** With both structures in one
>   common frame and the correct mean, the identity is **exact to 2.65e−15** (max 2.52e−14) on all
>   126 targets. `s16/retract_law.py`, `s16/retract_exact.py`, `s16/retract_FINDINGS.md`.

> ~~**d_avg ≈ √(r² − (s/2)²)**~~ → **d_avg = √( (r₁²+r₂²)/2 − (s/2)² )**, an exact identity

*(now on all 126 targets)*

| | n = 24 | **n = 126** |
|---|---|---|
| distogram fit | 3.270 | **3.622** |
| pool fit | 3.469 | **3.659** |
| structural disagreement `s` — **native-free** | 2.980 | **3.138** (3.506 in the common frame) |
| raw per-pair error correlation | +0.627 | **+0.602** |
| ~~**law prediction**~~ | ~~2.970~~ | ~~**3.205**~~ → **3.292** |
| **measured coordinate average** | 3.105 | **3.367** |
| ~~**prediction error**~~ | ~~+0.135~~ | ~~**+0.162**~~ → **+0.075** |

~~A **parameter-free** identity predicting an unseen quantity to **0.162 Å** across 126 targets; it
cannot be rescued by refitting because there is nothing to fit.~~ **WITHDRAWN.** It predicts nothing
unseen — it is an exact algebraic rearrangement, and what Sprint 15 measured as its "prediction
error" was a mixture of the wrong mean (54%) and a frame mismatch (46%).

**A control correction that flips a sign.** A field originally named `coordavg_vs_best_single` was
computed as a **per-target minimum** over the two channels — a per-target ORACLE selection, not a
control. The two comparisons have opposite signs:

| comparison | coordinate average | restraint-level fusion |
|---|---|---|
| **vs the better channel chosen globally — NATIVE-FREE** | **−0.255 [−0.370, −0.145]** | −0.071 [−0.156, +0.014] |
| vs a per-target ORACLE pick — *a ceiling, not a control* | +0.216 [+0.133, +0.302] | +0.400 [+0.262, +0.546] |

**The honest result is that coordinate averaging beats the better single channel by −0.255 Å with an
interval excluding zero.** Reported against the oracle it would have read as a significant loss.

**And a second sign flip from the small sample:** at n = 24 restraint-level fusion beat coordinate
averaging (3.046 against 3.105); at n = 126 it is reversed (3.551 against 3.367), and only the
coordinate average clears its native-free control. Combining at the **structure** level beats
combining at the **restraint** level — consistent with realizability, since a single fit to blended
restraints must pick *one* coherent wrong structure while an average of two separately-fitted
structures can land between two different ones, which is where the law says the gain lives.

The law has **no free parameters**, so it cannot be rescued by refitting, and `s` is measurable
without the native — which turns "should we fuse these channels?" into a calculation on unlabelled
data. It also states when fusion matters: the gain is small whenever `s ≪ r` and large only as
`s → 2r`. Here `s/r = 0.88`. **An ångström from fusion needs channels that disagree about twice as
much as they are wrong**, and no pair in this project comes close. *(full-instrument run in flight)*

### 6.3 The diversity trap that invalidates a law this project uses constantly

CVaR-VQE at α = 1 returns **2.7 distinct configurations out of 4,096 draws**. Its "set mean" looks
0.847 Å better than its control while its actual coordinate average is worse.

> **The set-mean law (`d_out = 1.16·d_set_mean + 0.04·d_set_best`) was fitted on sets with real
> diversity and is out of domain on a collapsed set. A method that concentrates improves the proxy
> and degrades the output.**

Every arm in the programme that invokes it now reports its set's diversity alongside its RMSD.

---

## 7. THE PROJECTION GAP, AND A CORRECTION TO MY OWN NUMBER

| arm | mean | vs `project_raw` |
|---|---|---|
| `project_raw` — **is** the production `fit_ca` | **3.204** | — |
| `project_scaled` — native-free, distogram scale | 3.327 | **+0.123 [+0.065, +0.184]** |
| `project_scaled_deb` — native-free, debiased | 3.289 | +0.085 [+0.040, +0.131] |
| **`ORACLE_scale_true`** | **2.929** | **−0.275 [−0.397, −0.174]** |

**There is 0.275 Å on the table** in the *production* pipeline's final stage. **And the distogram
cannot supply it**: the truth wants a scale of 1.0412 and the distogram asks for 1.0997, because its
own +0.509 Å expansion bias is baked into the estimate. Applied, it costs +0.123 Å. *(pool-referenced
arms in flight)*

**Correction.** I wrote in four places that coordinate averaging contracts the backbone by **25.8%**,
reusing a Sprint 14 figure without checking it. Measured directly: the production consensus is
contracted by **3.5%** against the true distances and 9.5% against the predicted ones. The number is
withdrawn; the direction survives, the magnitude was overstated sevenfold.

---

## 8. THE QUANTUM RESULT

Exact statevector gradients, no shot noise, budget charged only for readout draws, control given the
same draws. 12-qubit sub-registers of the nine enumerated targets, 200 iterations, 3 seeds, 4,096
draws, paired over 27 cells per α. **Positive means VQE is worse than best-of-N from its own
untrained initialisation.**

| readout | α = 1.0 | α = 0.25 | α = 0.05 | α = 0.01 |
|---|---|---|---|---|
| argmin (by objective) | −0.036 null | **+0.251 SIG** | +0.068 null | +0.048 null |
| **top-20 coordinate average** | +0.068 null | **+0.372 SIG** | +0.197 null | +0.177 null |
| **top-75 coordinate average** | +0.147 null | **+0.392 SIG** | **+0.276 SIG** | +0.164 null |
| drawn-set MEAN | **−0.847 SIG** | **−0.372 SIG** | **−0.191 SIG** | −0.066 null |
| drawn-set BEST | **+1.092 SIG** | **+0.345 SIG** | **+0.155 SIG** | +0.087 null |

**Nothing in the table is a significant win for VQE on a structural readout.** At α = 1 the drawn set
is 0.847 Å better on the mean and 1.092 Å worse on the best — two significant effects of opposite
sign shrinking monotonically to zero as α falls. **That is what concentration is, measured directly.**

**Competing explanations eliminated rather than invoked:** not barren plateaus (measured); not
ill-conditioning (the gradient's null-space share is 0 to 1e-30, and the bottom eigenvalue decile
carries 0.0016–0.0043 of its squared norm where uniform would be 0.10); not a broken CVaR (4.2e-14
against Rockafellar–Uryasev); not shot noise (exact gradients).

**α is a diversity dial**: 2.7 → 408 distinct configurations, entropy 0.73 → 6.58 bits.

**QNG and classical natural gradient are the same algorithm here** — the classical Fisher information
equals four times the Fubini–Study metric to 3.3e-16, removing an arm of the design space. And the
previously published condition numbers were a **seed-averaging artefact**: recomputed per point the
median is 735 and the maximum 1.5e6 against a recorded 2.2–3.7, and the `block` ansatz at depth ≥ 2 is
**exactly rank-deficient**.

*(whether the negative survives a change to a generative objective is in flight; the instrument is
verified and Sprint 14's flagship table reproduces exactly)*

---

## 9. WHAT WOULD BE NEEDED FOR 2.5 Å

### 9.1 The distance channel: a stated requirement

True distances corrupted under three realistic error models at eight severities each, refitted
through the identical machinery, 40 targets. The models inject different magnitudes at the same
nominal severity, so they are put on one comparable axis — the **effective RMS of the injected
error** (`sep_scaled`'s exact moment factor is √E[sep²]/E[sep] = 1.1331; `outliers`' is
√(0.09 + 0.05·36) = 1.375).

| error model | required effective RMS for 3.0 Å | for 2.5 Å | for 2.0 Å |
|---|---|---|---|
| `abs_gauss` — i.i.d. Gaussian | > 3.00 | **2.87** | **1.87** |
| `sep_scaled` — error growing with separation | > 3.40 | **2.95** | **1.85** |
| `outliers` — 5% of pairs badly wrong, rest nearly right | **> 4.12** | **> 4.12** | **> 4.12** |

**The distogram's actual RMSE is 3.704 Å and it emits 3.644 Å.**

> **Reaching 2.5 Å requires improving the distance channel from 3.70 Å RMSE to about 2.87 Å — a 22%
> reduction — and that is the requirement for an i.i.d. channel. Ours is worse than i.i.d. at equal
> magnitude by about 1.1 Å. Reaching 2.0 Å requires roughly halving the error.**

### 9.2 The controlled demonstration that shape beats magnitude

Two rows settle the mechanism, because here the error's structure is **imposed** rather than
observed, so the causal direction is not in question:

| channel | effective error RMS | emitted RMSD |
|---|---|---|
| outlier-shaped | **4.12 Å** | **1.993 Å** |
| i.i.d. Gaussian | 3.00 Å | 2.561 Å |
| **the real distogram** | **3.70 Å** | **3.644 Å** |

**A larger error concentrated in a few pairs is far cheaper than a smaller error spread evenly, and
both are far cheaper than the real distogram at comparable magnitude.** The ordering
`outliers ≪ i.i.d. < real` is the whole argument in one line: the fit can *ignore* a few
catastrophic restraints, cannot ignore many mildly wrong ones, and is actively *misled* by errors
that are mutually consistent.

**And a clean null.** `sep_scaled` and `abs_gauss` require the same accuracy (2.95 against 2.87 for
2.5 Å; 1.85 against 1.87 for 2.0 Å). Once magnitude is matched, making the error grow with sequence
separation costs essentially nothing — so the distogram's separation-dependent profile is **not**
what makes it expensive.

### 9.3 The torsion channel, where the structure runs the other way

| target | verdict on the i.i.d. surface | requirement |
|---|---|---|
| 3.0 Å | **supported** | coverage 1.0 at ≤ 26°, or 0.8 at ≤ 13°, or 0.5 with terminal gaps at ≤ 11° |
| 2.5 Å | **marginal** | coverage ≥ 0.9 at ≤ 15°, or ≥ 0.6 with terminal gaps |
| 2.0 Å | **not supported** by any i.i.d.-class channel | 16.3° at full coverage; best channel 67.6° |

**But the i.i.d. surface over-prices real torsion channels**: real torsion errors are **0.6–1.7 Å cheaper** than i.i.d. errors of the same magnitude across
**native-free** arms, and up to **2.9 Å cheaper for the ORACLE best-matching pool window**. Alignment is
**0.565 (ORACLE window)** and **0.665 (the incumbent's projected torsions)** against a **0.945**
random-direction null — but the control that matters is that a **zero-information constant α-helix
already scores 0.709**, so only those two arms clear it. Sign-flipping the ORACLE window's errors costs
**+1.596 Å [+1.385, +1.813]**; for the native-free retrieval mean the same operation is **null**
(−0.150 [−0.366, +0.070]). The 64°-RMS channel emitting **1.77 Å** where i.i.d. of the same size emits 4.71 Å is the
**ORACLE** window. And **where the
coverage gaps fall beats how many there are**: terminal gaps at coverage 0.7 give 1.716 Å against
3.505 Å for mid-chain — a **1.79 Å** spread that no channel exploits.

So the two channels' error structures point in opposite directions: alignment makes the torsion
channel **cheaper** than its magnitude implies, and realizability makes the distance channel **more
expensive**. Whether alignment can be *engineered* rather than inherited is the one open question
that could move 2.0 Å from unsupported to reachable without any improvement in raw accuracy.
*(in flight)*

### 9.4 What is actionable

The requirement is now specific enough to act on. It is not "make the distogram better." It is
**reduce its RMSE by 22%, or make its errors incoherent** — and §9.2 shows an error that is *larger*
but unstructured in the right way is worth more than a smaller one. Nothing in this project has ever
tried to shape a predictor's error rather than shrink it.

## 10. THE HONEST BOTTOM LINE

The incumbent emits **3.204 Å**. The generative cascade emits **3.321 Å**. Nothing measured in this
sprint closes the 0.8 Å to 2.5 Å, and the sprint's own mechanism (§5) explains why the obvious
routes cannot: they correct magnitudes and the error is a direction.

What the sprint delivers instead is a set of measurements that change where effort should go — the
ceiling is 0.611 Å rather than 1.95 Å; the barrier is restraint *coherence* rather than restraint
*accuracy*; fusion is priced by a native-free geometric law; and a variational optimiser loses to its
own untrained initialisation on a control the literature does not run.

---

† **Absolute constants carry a start-draw error bar.** The multi-start initialisation was seeded with
`hash(pdb)`, which Python salts per process, so a run's absolute mean varies with the draw: measured
**sd 0.132 Å** across four independent draws, with an independent workstream's draw returning
**0.569 Å** for this same arm. **Quote this ceiling as ≈ 0.6 Å with a start-draw sd of 0.13 Å, not as
0.611.** Paired *differences* are far more stable (the same quantity measured at 0.381 / 0.416 / 0.486
across draws), which is why the comparative claims in this document survive while its absolute
constants need the bar. Seeding is now stable (`s15/seed.py`), and the frozen protocol is re-run under
it before the benchmark is touched.

‡ **An empirical false-positive floor of ≈ 0.08 Å applies to single-draw paired comparisons.**
Demonstrated directly: a comparison that is **exactly zero by construction** (maximum likelihood
against a constant-width Gaussian *is* least squares) returned **+0.081 [+0.014, +0.169]** on one
start draw and −0.003 on another — a confidence interval excluding zero on a null effect. The
null-calibrated concentration test flags that row, and only that row, as FAIL (p = 0.024). **Effects
at or below 0.08 Å are not resolvable here without replication across draws**, and are marked ‡
wherever they appear.
