# SPRINT 15 — RESTRAINT agent findings

The binned distogram output as a representation: what it costs, and whether a continuous
density buys it back.

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Everything here runs as `python -m s15.quant <check|mm|budget|ladder|pred|tolcheck|startnoise|
repro|mech>` from the repo root; results land in `s15/results/quant_*.json`. `python -m
s12.instrument` reproduced `shipped 3.4540004952559396, pool_best 1.7108244199364904,
top75_best 2.3061526409453816, synthesis_fit 3.2040761603809194, n_zero_recall 18` at the start
and at the end of this work. The 60-target protected benchmark was not read, listed or touched.
Every derivative in `s15/quant.py` is verified against central differences at
**cosine 1.000000000000**, including the full torsion-gradient path
(`quant_derivcheck.json`); the single exception is the frozen pre-fix table, deliberately kept
broken, at 0.8099.

---

## HEADLINE, in the four sentences the brief asks for

**(a) How much accuracy does the binned representation cost?**
**+0.381 Å [+0.156, +0.614]**, reproduced at two further independent start draws as +0.416 and
+0.486 (mean **≈ +0.43 Å**), to an oracle that knows the true distance matrix perfectly
(0.569 → 0.950 Å, worse on 111 of 126 targets, same sign on 5/5 folds) — but the cost is bin
**placement**, not bin count, because a uniform 0.5 Å grid over the same range costs only
+0.243 and a uniform 0.25 Å grid +0.085 with a CI crossing zero; and against the *predictor's*
own 3.704 Å error the binning is worth only **2.7% of the error variance**, which caps what any
de-binning fix can buy on real restraints before a single structure is fitted.

**(b) Does a continuous density recover it?**
**Partly, and only in the ORACLE setting**: a kernel density tied to local bin width takes
binned-ML from +0.528 to **+0.291** — 45% recovered, and better than least squares on snapped
bin centres — while a mass-exact PCHIP density recovers essentially nothing (+0.504), so it is
the **smoothing** and not the mass fidelity that helps; on *predicted* restraints the entire
spread across thirteen representations is 0.33 Å against a 3.03 Å oracle-to-predicted gap.

**(c) Does ML then beat least squares?**
**No — but it no longer loses, and the recorded "ML loses" was a bug.** Through the pre-fix
table ML was **+0.230 [+0.093, +0.365]** worse; with the table corrected every ML arm has a
negative point estimate (−0.041 to −0.102) and **every confidence interval touches or crosses
zero**, the leave-fold-out arm that also chooses the density family is the weakest of all
(**−0.040 [−0.153, +0.074]**, W/L 60/66, median +0.025), and a mathematically-identical
comparison in this same harness produced a **false** +0.081 [+0.014, +0.169] from the start
draw alone — so nothing here clears its own measurement floor.

**(d) Practical recommendation for the architecture.**
Do **not** spend a sprint on the output representation: **respace the bins uniformly at ≈0.5 Å
if the distogram is ever retrained** (free, and worth ~0.14 Å of restraint RMSE), keep
consuming the distogram through its **mean and sd** rather than its full density, and spend the
effort instead on **calibration** — ML beats least squares by −0.13 to −0.18 Å with a CI
excluding zero on exactly the 52 of 126 targets where the distogram is a better probability
than a pooled marginal, and by nothing on the other 74, which is where the remaining 3.03 Å
actually lives.

---

## R0. THE 1.06 Å "QUANTISATION CEILING" THAT CREATED THIS TASK WAS AN INDEXING BUG

**DEMONSTRATED. This is a retraction of the premise I was given, found before any experiment
was run, and independently found by the coordinator within the same hour.**

`s15/distml.LogPTable` computed its bin index as `floor((d - c[0]) / (c[1] - c[0]))` — it
assumed the 17 bin centres were a uniform grid of width `c[1]-c[0] = 0.75 Å`. They are not:

    4.0 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.5 9.5 10.5 11.75 13.25 15.0 17.5 21.0 25.0

so the spacing runs 0.5 Å at short range and 4.0 Å at long range. Consequences, verified by
direct evaluation:

| query | bin the table read | bin it should have read |
|---|---|---|
| 5.0 Å | centre 4.75 | centre 4.75 (correct by luck) |
| 8.0 Å | centre **6.75** | centre 7.75 |
| 12.0 Å | centre **10.50** | centre 11.75 |
| 16.0 Å | centre **21.00** | centre 15.00 |
| > 15.25 Å | all collapse into the last interval | — |

There was a **second, independent defect in the same class**, which I have measured and which
had not been stated: above `c[0] + (B-2)*0.75 = 15.25 Å` the fractional coordinate `t` is
clipped to 1, so the returned **value is flat while the returned derivative is non-zero**. The
table therefore handed L-BFGS a gradient inconsistent with its own objective on exactly the
long-range restraints. A central-difference check confirms it: with the kinks of every
piecewise representation excluded, the frozen pre-fix class scores **cosine 0.8099** against
finite differences, with `max|err| 6.15`, and **every** disagreeing sample lies above 16.13 Å.

The corrected class scores cosine **1.000000000000** — as does every representation built in
this module. The defect was in the **index map**, not in the calculus, which is precisely why
it survived whatever gradient checking it received. That is the transferable lesson: *a
gradient check validates a function against itself and cannot see a wrong lookup.*

**Status of the numbers this task was built on.** `ORACLE_ml_true = 1.756` and
`ml_full = 3.021` are **VOID**. The frozen buggy class is preserved as
`s15.quant.LogPTableUniformBUG`, deliberately never to be fixed, so the ladder can measure how
much of the apparent effect was the bug (see R2).

**Independent replication.** `s15.quant.LogPGrid` and the coordinator's corrected
`s15.distml.LogPTable` were written separately from the same specification and agree to
`max|Δ log P| = 0.0e+00`, `max|Δ grad| = 0.0e+00`. Two independent implementations, not two
copies.

---

## R1. THE ARITHMETIC CEILING, BEFORE ANY STRUCTURE IS FITTED

**DEMONSTRATED** (native distances used post hoc to score a representation error; the quantity
is a property of the *binning*, not of any target).
`s15/results/quant_snaperr.json`, `quant_noisebudget.json`. 126 targets, 8,549 pairs.

Snapping a **perfect** distance matrix to its bin centre is the entire quantisation operation.
Its error, measured directly:

| band | n pairs | snap-to-bin RMSE | snap bias | uniform 0.5 Å grid | uniform 0.25 Å grid |
|---|---|---|---|---|---|
| short, \|i−j\| 2–5 | 4768 | **0.260** | +0.010 | 0.143 | 0.073 |
| mid, \|i−j\| 6–7 | 1628 | **0.499** | +0.056 | 0.145 | 0.073 |
| long, \|i−j\| ≥ 8 | 2153 | **1.055** | −0.066 | 0.144 | 0.071 |
| **all** | 8549 | **0.604** | −0.001 | 0.144 | 0.073 |

Two things follow immediately and neither needs a fit.

**The binning is 4.2× worse than a uniform 0.5 Å grid over the same range** (0.604 vs 0.144),
and the whole of that penalty is long-range: at \|i−j\| ≥ 8 the non-uniform grid costs 1.055 Å
against a uniform grid's 0.144 Å, a factor of **7.3**. So the coordinator's prediction that
*bin placement*, not bin count, carries the loss is confirmed at the restraint level.

**But against the predictor's own error it is small.** The distogram's measured RMSE is 3.704 Å
globally (coord K2). Quantisation adds 0.604 Å in quadrature, i.e.
`(0.604/3.704)² = 2.7%` of the error variance; at long range `(1.055/5.574)² = 3.6%`. A
purely analytic bound using each bin's uniform variance `w²/12` gives 1.3%, in the same place.

> **Whatever the binned representation costs an ORACLE, it can cost the PREDICTIVE arms at
> most a few percent of their error variance, because the predictor's own error is six times
> larger than the quantisation it is expressed in.** This caps the hypothesis before a
> structure is fitted, and the predictive experiment (R3) is a test of that cap, not of the
> hypothesis's plausibility.

Note also that the *populated* long-range bins are much finer than the grid's worst: the mean
width of the bin a real CA–CA distance actually falls into is 1.86 Å at \|i−j\| ≥ 8, not 4.0 Å.
Peptides of 9–16 residues rarely reach the 16–23 Å bins at all.

---

## R2. THE ORACLE QUANTISATION LADDER — the binning costs an oracle 0.381 Å

**ORACLE DIAGNOSTIC throughout. Every arm knows the TRUE distance matrix; only the
REPRESENTATION of that knowledge changes. Nothing in this table can ever be a predictive
headline.** `s15/results/quant_ladder.json`, all 126 targets, 6 native-free starts, selection
by objective only, RMSD read post hoc.

| arm | mean | median | <2 Å | <2.5 Å | FAIL18 | vs continuous (paired, 95% CI) | W/L |
|---|---|---|---|---|---|---|---|
| **ORACLE continuous least squares** | **0.569** | 0.062 | 0.87 | 0.88 | 1.033 | — | — |
| ORACLE snapped to **bin centre** | 0.950 | 0.388 | 0.81 | 0.86 | 1.592 | **+0.381 [+0.156, +0.614]** | 15/111 |
| ORACLE snapped to uniform **0.5 Å** | 0.812 | 0.166 | 0.82 | 0.84 | 1.520 | +0.243 [+0.056, +0.428] | 17/109 |
| ORACLE snapped to uniform **0.25 Å** | 0.654 | 0.104 | 0.86 | 0.88 | 1.396 | +0.085 [−0.123, +0.298] | 26/100 |
| ORACLE bin-snap, **short 2–5 only** | 0.801 | 0.203 | 0.83 | 0.86 | 1.159 | +0.232 [+0.032, +0.438] | 23/103 |
| ORACLE bin-snap, **mid 6–7 only** | 0.765 | 0.207 | 0.85 | 0.89 | 1.389 | +0.196 [+0.001, +0.393] | 18/108 |
| ORACLE bin-snap, **long ≥8 only** | 0.927 | 0.290 | 0.79 | 0.83 | 2.102 | +0.358 [+0.113, +0.603] | 18/108 |
| ORACLE ML, Gaussian w=0.3, **sampled at centres** | 1.290 | 0.691 | 0.73 | 0.79 | 2.427 | +0.722 [+0.472, +0.978] | 18/108 |
| ORACLE ML, Gaussian w=0.6, **sampled at centres** | 1.097 | 0.726 | 0.83 | 0.87 | 1.930 | +0.528 [+0.344, +0.718] | 18/108 |
| ORACLE ML, Gaussian w=1.2, **sampled at centres** | 1.299 | 0.743 | 0.71 | 0.83 | 2.188 | +0.730 [+0.488, +0.983] | 20/106 |
| ORACLE ML, Gaussian w=0.6, **evaluated continuously** | 0.650 | 0.076 | 0.84 | 0.87 | 1.206 | +0.081 [+0.014, +0.169] ‡ | 62/64 |
| ORACLE ML, w=0.6 histogram → **KDE density** | 0.860 | 0.313 | 0.83 | 0.86 | 1.675 | +0.291 [+0.047, +0.540] | 23/103 |
| ORACLE ML, w=0.6 histogram → **PCHIP density** | 1.073 | 0.523 | 0.83 | 0.87 | 2.195 | +0.504 [+0.295, +0.722] | 20/106 |
| **ORACLE ml_true through the VOID pre-fix table** | 2.650 | 2.755 | 0.33 | 0.43 | 3.287 | **+2.081 [+1.837, +2.332]** | 7/119 |

‡ **this row is a demonstrated FALSE POSITIVE — see R7.** ML against a continuous Gaussian of
constant width *is* unweighted least squares, so this row's true value is exactly zero.

Null-calibrated concentration verdicts (`s15.info_lib.concentration_verdict`): **PASS on every
row except the ‡ row, which FAILS at p = 0.024** — the one row we can prove is spurious is the
one the concentration test flags. `ORACLE_snap_bin` is positive on **5 of 5 folds**
(+0.851, +0.616, +0.013, +0.051, +0.369), median difference +0.222.

### What the ladder establishes

**(1) The quantisation cost to an oracle is +0.381 Å [+0.156, +0.614], 111 of 126 targets
worse.** That is the number the task asked for, and it is the first valid measurement of it.

**(2) Bin PLACEMENT, not bin count, is most of it.** A **uniform 0.5 Å** grid over the same
range costs +0.243 and a **uniform 0.25 Å** grid costs +0.085 with a CI crossing zero. The
shipped grid has 17 bins and a 0.5 Å floor; a uniform 0.25 Å grid over 3.5–27 Å would need 94.
So *the loss is bought back by spending bins where the distances are, not by having more of
them* — and R1's restraint-level RMSEs say the same thing (0.604 Å for the shipped grid,
0.144 Å for a uniform 0.5 Å grid).

**(3) The coordinator's long-range prediction is CORRECT IN RANK BUT NOT IN MAGNITUDE, and the
decomposition is strongly sub-additive.** Quantising long-range pairs alone costs the most
(+0.358), but quantising *short*-range pairs alone already costs +0.232, and the three bands
sum to **+0.785 against a joint cost of +0.381**. So it is **false** that "long-range pairs
carry nearly all of the loss": long range carries the largest single share, but no band is
negligible and the parts do not add.

> The sub-additivity is itself a mechanism worth stating. Quantising **one** band leaves that
> band's corrupted restraints fighting **exact** restraints elsewhere, and the fit is dragged
> into an inconsistent compromise. Quantising **all** pairs produces errors that are mutually
> consistent — they all come from one snapped matrix — and the fit finds a coherent, if
> shifted, structure. **A partial-corruption arm therefore over-states the damage, and any
> band decomposition of a geometric fit must report the joint value beside the parts.**

**(4) Maximum likelihood through a sampled histogram costs MORE than least squares on snapped
means** (+0.528 for w=0.6 against +0.381 for `snap_bin`), and the width sweep is
**non-monotone** (+0.722 at w=0.3, +0.528 at w=0.6, +0.730 at w=1.2) — a narrow law is
punished for having no mass between bin centres, and a wide one for being uninformative.

**(5) Most of the "1.06 Å ceiling" that created this task was the indexing bug.** Within one
start-matched run, the pre-fix table costs **+2.081** where the same information through a
correct table costs **+0.528**. So **75% of the apparent effect was the defect** and 25% was
representation. The retraction in R0 is quantified here.

**(6) A continuous density recovers a real part of the ORACLE binning cost — and which
reconstruction you use matters more than that you used one.** KDE takes the binned-ML cost
from +0.528 to **+0.291**, recovering 45% and landing *better than least squares on snapped
bin centres* (+0.381). PCHIP recovers almost nothing (+0.504). The two differ in exactly the
way the mechanism predicts: PCHIP reproduces each bin's mass **exactly**, so it faithfully
preserves the histogram's blockiness; the KDE deliberately does not, and it is the smoothing —
not the mass fidelity — that buys the accuracy back.

---

## R3. THE PREDICTIVE ARMS — with the table fixed, ML no longer loses, and does not clearly win

**DEMONSTRATED (native-free; the only ORACLE quantity is the post-hoc RMSD).**
`s15/results/quant_pred.json`, all 126 targets. Incumbent 3.204 Å. Base arm `ls_mean_sd`.

| arm | mean | median | <2 Å | <2.5 Å | FAIL18 | vs least squares (paired, 95% CI) | W/L |
|---|---|---|---|---|---|---|---|
| `ls_mean_sd` (the s15/distgeo predictive arm) | 3.620 | 3.396 | 0.15 | 0.21 | 5.921 | — | — |
| `ml_full` through the **VOID pre-fix table** | 3.850 | 3.865 | 0.13 | 0.17 | 5.828 | **+0.230 [+0.093, +0.365]** | 51/75 |
| `ml_full_fixedgrid_mass` (untuned) | 3.569 | 3.432 | 0.21 | 0.25 | 6.051 | −0.051 [−0.168, +0.068] | 74/52 |
| `ml_full_fixedgrid_density` (untuned) | 3.579 | 3.532 | 0.21 | 0.26 | 6.175 | −0.041 [−0.150, +0.073] | 72/54 |
| `ml_kde_a0.5` *(EXPLORATORY, tuned on these targets)* | 3.528 | 3.512 | 0.21 | 0.25 | 5.990 | −0.092 [−0.192, −0.000] | 67/59 |
| `ml_kde_a1.0` *(EXPLORATORY)* | 3.616 | 3.498 | 0.18 | 0.25 | 6.110 | −0.004 [−0.097, +0.086] | 61/65 |
| `ml_kde_a2.0` *(EXPLORATORY)* | 3.679 | 3.595 | 0.16 | 0.24 | 5.891 | +0.059 [−0.039, +0.155] | 52/74 |
| `ml_gmm_K2` *(EXPLORATORY)* | 3.582 | 3.433 | 0.18 | 0.25 | 6.057 | −0.038 [−0.144, +0.072] | 74/52 |
| `ml_gmm_K3` *(EXPLORATORY)* | 3.555 | 3.475 | 0.19 | 0.27 | 6.018 | −0.065 [−0.159, +0.030] | 74/52 |
| `ml_pchip` (untuned) | 3.518 | 3.403 | 0.21 | 0.25 | 6.027 | −0.102 [−0.215, +0.014] | 69/57 |
| **`ml_kde_LFO`** *(CONFIRMATORY, leave-fold-out)* | 3.528 | 3.512 | 0.21 | 0.25 | 5.990 | **−0.092 [−0.192, −0.000]** | 67/59 |
| **`ml_gmm_LFO`** *(CONFIRMATORY)* | 3.555 | 3.475 | 0.19 | 0.27 | 6.018 | −0.065 [−0.159, +0.030] | 74/52 |
| **`ml_cont_LFO`** *(CONFIRMATORY, family also chosen LFO)* | 3.580 | 3.500 | 0.21 | 0.24 | 6.082 | **−0.040 [−0.153, +0.074]** | 60/66 |

Leave-fold-out choices: `ml_kde_LFO` picked α = 0.5 in **all five** folds; `ml_gmm_LFO` picked
K = 3 in all five; `ml_cont_LFO` picked PCHIP in folds 0/1/3 and KDE α = 0.5 in folds 2/4.
That the bandwidth and component count are unanimous across folds is a good sign for the
*hyperparameter*; the family split is not.

Null-calibrated concentration verdicts: **PASS on every arm**, but every one carries the
"mean/sd small: test has little power" flag (best |mean/sd| = 0.17). Per the brief's own
trap-list, that is the case where a raw drop-top-10 threshold **would have misfired**: the
drop-top-10 differences are +0.019 for `ml_kde_LFO` and +0.085 for `ml_cont_LFO` — a sign flip
— with `top10_share` of 1.19 and 2.97. The null-calibrated test says those shares are **not**
extreme for effects this small relative to their spread. Both statements are reported together,
as one verdict, exactly as the brief requires: **not demonstrably concentrated, and not
demonstrably present.**

### What the predictive table establishes

**(1) The pre-fix table was making ML significantly worse: +0.230 [+0.093, +0.365].** The
sprint's recorded "ML loses to least squares, 3.021 vs 2.877" was the bug, not a finding. That
is the concrete instance of a representation defect masquerading as a scientific result which
the coordinator asked to have on the record.

**(2) With the table corrected, ML is no longer worse — and is not clearly better.** Every
corrected ML arm has a negative point estimate (−0.041 to −0.102) and **every confidence
interval touches or crosses zero.** The single arm whose CI excludes zero does so at
**−0.000**. And the whole range of these point estimates lies **at or below the measurement
floor established in R7**, where a comparison whose true value is provably exactly zero
returned +0.081 [+0.014, +0.169] from the start draw alone. **No arm in this table clears its
own noise floor, and none should be reported as a gain.**

**(3) The honest arm is the weakest, which is the signature of selection noise.** The
exploratory best (`ml_pchip`, −0.102) and the LFO arm that also chooses the family
(`ml_cont_LFO`, **−0.040**, W/L 60/66, median **+0.025** — i.e. the median target is *worse*)
differ by a factor of 2.5. When the honest version of an effect is much smaller than the
best-of-sweep version, the sweep was reading noise.

**(4) Every arm here is worse than the incumbent's 3.204 Å**, so none of this changes the
pipeline. The whole spread across thirteen representations of the same distogram is 0.33 Å,
against a 3.03 Å oracle-to-predicted gap. **The representation of the restraints is not where
the accuracy is.**

---

## R4. MULTIMODALITY IS RARE, AND IT IS NOT CONCENTRATED AT LONG RANGE

**DEMONSTRATED, and it REFUTES the mechanism the ML-over-least-squares argument rests on.**
`s15/results/quant_multimodal_density.json`, all 126 targets, 8,549 pairs.

**Definition, stated so it can be argued with.** Work on the **density** `h_b = p_b / width_b`,
not the raw mass — bins differing eightfold in width are not comparable as masses, and using
mass manufactures a spurious mode at every wide bin. A pair is **strictly bimodal** when it
has two local density maxima at least 2 bins apart, the minimum density strictly between them
is at most **0.5×** the smaller peak (a real valley, not a shoulder), and the smaller mode's
basin — the bins on its side of the valley — carries at least **0.15** of the total
probability. The loose count (any two local maxima, no separation or depth requirement) is
reported beside it so a reader can see how much work the strictness does.

| band | n pairs | ≥ 2 local maxima | **strictly bimodal** | mean n modes |
|---|---|---|---|---|
| short, \|i−j\| 2–5 | 4768 | 0.786 | **0.099** | 1.97 |
| mid, \|i−j\| 6–7 | 1628 | 0.697 | **0.057** | 1.86 |
| long, \|i−j\| ≥ 8 | 2153 | 0.752 | **0.059** | 2.00 |
| **all** | 8549 | 0.761 | **0.081** | 1.96 |

Control on the raw-mass basis (`quant_multimodal_mass.json`): 0.107 / 0.068 / 0.063 / 0.088 —
the answer does not depend on the basis, so the density normalisation is not driving it.

Two readings, and I am obliged to state the second plainly.

**Local maxima are common (76%) but almost all of them are shoulders.** The mean mode count is
1.96, i.e. the typical histogram wiggles once; only **8.1%** of pairs carry two genuinely
separated modes with a real valley and non-trivial mass on both sides.

> **REFUTED: the coordinator's stated mechanism — that long-separation pairs are
> "either in contact or not" and therefore bimodal, which is where ML's extra information over
> a mean would come from — is false for this distogram.** Strict bimodality is *lower* at
> \|i−j\| ≥ 8 (0.059) than at \|i−j\| 2–5 (0.099). The contact/non-contact structure the
> argument assumes is not in the predicted histograms. On 92% of pairs the distribution is a
> single lump plus a shoulder, and a mean-and-sd summary of a single lump loses very little.
>
> **The ML-over-least-squares argument is therefore weak on its own terms, independent of any
> fitting result.** I say so before reporting the fits, so it cannot be read as a
> rationalisation of them.

---

## R4c. THE COMPETING EXPLANATION: THE DISTOGRAM IS A GOOD RANKER AND A BAD PROBABILITY

**DEMONSTRATED, and independently replicated here.** Relayed from the INFO workstream and
re-derived in `s15.quant.per_target_ipred` from `s15/info_mi.distogram_information`'s
definition, without using INFO's code path.

`I_pred = H0 − CE` in bits per CA–CA pair, where `H0` is a leave-one-target-out pooled
distance marginal's cross-entropy on the native bin label and `CE` is the distogram's own.
My reimplementation returns **mean −0.818 bits/pair, median −0.138, positive on 52 of 126
targets** — identical to INFO's reported −0.818 / −0.138 / 52. Two independent
implementations, so the statistic is not an artefact of either.

This matters because it supplies a mechanism for ML losing to least squares that has **nothing
to do with binning**:

> Maximum likelihood consumes the predicted distribution **as a probability** — every feature
> of its shape enters the objective. Least squares consumes only its **first two moments**.
> If the shape is worse than a pooled marginal on 74 of 126 targets, ML is faithfully fitting
> a badly wrong density and least squares is discarding exactly the part that is wrong.

Under that account, a continuous density should **recover the quantisation cost and still
leave ML losing**, because removing the discretisation does not make an uncalibrated density
calibrated. The two questions therefore come apart and are reported apart:

- **Representation question:** does a continuous density close the ORACLE gap? (R2)
- **Calibration question:** does it let ML beat least squares on *predicted* restraints? (R3)

and the falsification test that separates them — does ML's per-target advantage track
`I_pred`, and does ML win on the 52 calibrated targets? — is R6.

---

## R4b. A SECOND, UNRELATED DEFECT FOUND ON THE WAY: THE MULTI-START DRAW IS NOT REPRODUCIBLE

**DEMONSTRATED.** Not part of my brief; found while checking that my ladder shared starts
with the coordinator's reference arm, and reported because it prices every cross-module
comparison in this sprint.

`s15/distgeo.py`, `s15/distml.py`, `s15/distcal.py` and (initially) this module all seed the
multi-start RNG with

    rng = np.random.default_rng(hash(pdb) % (2 ** 32))

Python's **string hash is salted per process** unless `PYTHONHASHSEED` is set, and
`s15/results/audit_env_lock.json` records it as `null`. Verified directly: two consecutive
interpreters return `hash('1A13') = 217586290314588545` and `2408026022170661001`. So **every
run draws a different set of pool starts.**

What this does and does not break:

- **Within one process it is harmless.** All arms in a single run share the same start set, so
  every paired comparison inside `distgeo.json`, `distcal.json`, `distml.json` or
  `quant_*.json` is start-matched and valid. My decomposition below is of this kind.
- **Across processes it is not.** Any figure from one module compared with a figure from
  another is comparing different start draws, and the published constants — the ORACLE
  0.611 Å and the predicted 3.644 Å of coord K1-CORRECTED among them — are **not bit-
  reproducible on re-running**.

`s15/results/quant_startnoise.json` measures the size of that noise on the cleanest arm in the
sprint: **sd of the mean 0.132 Å across four draws, range 0.273 Å, worst single target 5.15 Å**
— full table and the paired-versus-level decomposition in **R7**. The short version is that the
noise is large on an arm's absolute level and small relative to a paired difference, so
K1-CORRECTED's *conclusion* is safe and its *constants* need a spread attached.

**Recommendation, adopted.** The coordinator has added `s15/seed.py` with
`stable_seed`/`stable_rng` on `blake2b` and patched seven modules; `s15/quant.py` now uses
`SD.stable_rng(pdb, "ladder"|"pred"|"startnoise", ...)` too. Verified across two interpreters:
`stable_seed('1A13') = 3025285571` both times.

**Provenance of the runs reported below, stated rather than hidden.** `quant_ladder.json` and
`quant_pred.json` were **launched before the patch** and therefore ran under `hash(pdb)`
seeding. Following the coordinator's sprint-wide policy for the in-flight `cascade` run, they
are left to finish: their job is answered by a *within-run paired* comparison, and every arm
in each file shares one start set, so every comparison reported here is start-matched and
valid. Re-running them under `stable_rng` will move each arm by the start-draw noise measured
in `quant_startnoise.json` and will not move the paired differences by more than that. They
are **exploratory**, not frozen-protocol constants.

**And the part that belongs in the paper's reproducibility section.** This defect was
invisible to every determinism check the sprint already had — `s15/audit_determinism.py`,
`audit_repro_*.json`, `audit_env_lock.json` — because all of them run **within a single
interpreter** or compare against a **cached artefact**. A process-salted hash is identical
every time you ask it inside one process. *A reproducibility audit that never starts a second
interpreter cannot detect the most common source of irreproducibility in Python.* The
minimal fix is one line in any audit: run the check twice, in two subprocesses, and compare.

---

## R5. WHICH ARMS ARE TUNED, AND HOW

The brief forbids reporting a number tuned on the targets it is reported on. Stated
explicitly:

- **Untuned, no free parameters:** `ls_mean_sd`, `ml_full_fixedgrid_mass`,
  `ml_full_fixedgrid_density`, `ml_pchip`, and every arm in the ORACLE ladder except the
  Gaussian-width sweep (whose widths are declared *a priori* as 0.3/0.6/1.2 Å, reported as a
  sweep, and never selected from).
- **EXPLORATORY, tuned on the reported targets, and labelled as such in the table:**
  `ml_kde_a0.5`, `ml_kde_a1.0`, `ml_kde_a2.0`, `ml_gmm_K2`, `ml_gmm_K3`. These are the sweep.
  Their best member is **not** a result and is not quoted as one.
- **CONFIRMATORY, honest:** `ml_kde_LFO`, `ml_gmm_LFO`, `ml_cont_LFO`. For each of the five
  pinned folds, the bandwidth / component count / density family is chosen by mean RMSD on the
  **other four folds only**, then applied to the held-out fold. No target contributes to the
  choice of the hyperparameter under which it is scored. This is the same leave-fold-out
  pattern `s15/distcal.py` uses for the bias correction. The selection criterion is
  training-fold RMSD, which is native-derived — that is permitted by the brief for training
  under leave-fold-out discipline, and it is stated here rather than buried.

The five folds are the pinned `peptide_folds.json` folds. The 60-target protected benchmark
was not read, listed, or touched.

---

## R6. WHAT ML'S ADVANTAGE ACTUALLY TRACKS: NOT SHAPE, CALIBRATION

**ORACLE DIAGNOSTIC** — the splitting variable `I_pred` is scored against native distances, so
this is a mechanism explanation and **not** a usable selector. `s15/results/quant_mechanism_*.json`,
n = 126, Spearman rho against a 2000-draw permutation null.

If ML beats least squares *because* it can represent asymmetry and bimodality that a mean and
an sd cannot, its per-target advantage must track those features. If instead the INFO account
is right, it must track calibration. Three correlations, four arms:

| arm | rho(**strict bimodality**) | rho(**\|skewness\|**) | rho(**I_pred**) |
|---|---|---|---|
| `ml_full_fixedgrid_mass` | +0.104 (p 0.245) | −0.146 (p 0.104) | **−0.292 (p 0.001)** |
| `ml_full_fixedgrid_density` | +0.058 (p 0.519) | −0.107 (p 0.231) | **−0.217 (p 0.014)** |
| `ml_pchip` | +0.068 (p 0.448) | −0.069 (p 0.442) | **−0.186 (p 0.037)** |
| `ml_kde_LFO` | −0.054 (p 0.547) | −0.121 (p 0.176) | −0.149 (p 0.095) |

Permutation nulls are ≈ [−0.18, +0.18] throughout. Negative rho on `I_pred` means: *the better
calibrated the density, the more ML beats least squares.*

Splitting on the sign of `I_pred` (52 calibrated targets against 74 miscalibrated):

| arm | ML − LS where `I_pred` > 0 (n=52) | where `I_pred` ≤ 0 (n=74) |
|---|---|---|
| `ml_full_fixedgrid_mass` | **−0.165 [−0.298, −0.031]** | +0.029 [−0.140, +0.196] |
| `ml_full_fixedgrid_density` | **−0.127 [−0.250, −0.009]** | +0.020 [−0.148, +0.184] |
| `ml_pchip` | **−0.182 [−0.327, −0.034]** | −0.046 [−0.210, +0.109] |
| `ml_kde_LFO` | −0.104 [−0.228, +0.014] | −0.084 [−0.222, +0.043] |

**REFUTED, on its own terms: the shape hypothesis.** Bimodality does not predict where ML wins
on any arm — every rho sits inside its permutation null, and the sign is not even consistent.
Neither does skewness. This is the second, independent refutation of the same mechanism; R4
already showed strict bimodality is rare (8.1%) and *lower* at long separation.

**CONFIRMED: the INFO account.** `I_pred` predicts it on three of four arms with p ≤ 0.037, and
ML beats least squares by −0.13 to −0.18 Å with a CI excluding zero on precisely the 52 targets
where the distogram is a better probability than a pooled marginal, and by nothing at all on
the other 74. My reimplementation of `I_pred` reproduces INFO's −0.818 / −0.138 / 52-of-126
exactly (R4c), so the splitting variable is not an artefact of either implementation.

Two details worth keeping:

- **The effect is strongest for the rawest table and weakest for the most smoothed one**
  (−0.292 for the plain fixed-grid histogram, −0.149 for KDE α=0.5). That is exactly what the
  mechanism predicts: smoothing washes out the density shape, so the arm becomes less sensitive
  to whether the shape is right — and correspondingly less able to profit when it is.
- **The mechanism does not become a method.** The split is ORACLE. Turning it into a predictive
  gain would need a native-free per-target proxy for calibration, which is a new question, and
  Sprint 14's arithmetic (gain goes as the square of the weaker channel's skill) says a weak
  proxy would buy very little.

---

## R7. THE MEASUREMENT FLOOR — AND A DEMONSTRATED FALSE POSITIVE AT THE EFFECT SIZES IN PLAY

**DEMONSTRATED. This is the control that decides how to read R3, and it should be read before
any small paired difference elsewhere in Sprint 15.**

Maximum likelihood against a Gaussian of **constant** width is unweighted least squares:

    −Σ log N(d; d̂, σ)  =  (1/2σ²) Σ (d − d̂)²  +  P·log(σ√(2π))

an affine map with positive slope. Identical argmin, identical minimiser, identical everything.
`s15.quant.tol_check` fits the true distance matrix three ways from one shared start set:

| arm | mean | median | vs `ls_ref` (paired, 95% CI) | W/L |
|---|---|---|---|---|
| `ls_ref` — Σ(d − d_true)² | 0.625 | 0.051 | — | — |
| `ml_gauss_const` — the identical objective through `fit_ml` | 0.622 | 0.051 | −0.003 [−0.071, +0.060] | 64/62 |
| `ml_gauss_noconst` — the same, constant term dropped | 0.607 | 0.050 | −0.019 [−0.057, +0.001] | 59/66 |

So the least-squares/ML machinery agrees with itself to within ±0.019 Å, W/L ≈ 50/50, when the
objective is genuinely the same. **Good news for the harness.** But now compare that row with
the same comparison inside the ladder:

| the same mathematically-identical comparison, two independent runs | result |
|---|---|
| in `quant_ladder` (start draw A) | **+0.081 [+0.014, +0.169]** — CI excludes zero |
| in `quant_tolcheck` (start draw B) | **−0.003 [−0.071, +0.060]** — CI contains zero |

> **A comparison whose true value is provably exactly zero produced a paired mean difference of
> +0.081 Å with a 95% confidence interval excluding zero, from the multi-start draw alone.**
> The null-calibrated concentration test independently flags that same row, and only that row,
> as **FAIL at p = 0.024**. Two diagnostics, one spurious result, caught.

That is the honest measurement floor for this instrument at n = 126 with 6 starts, and it is
the same size as every predictive effect in R3 (−0.040 to −0.102). It is also larger than
several differences the sprint has treated as real, including the −0.022 Å restraint effect
carried over from Sprint 12.

### The other half: an arm's absolute mean is not reproducible either

`s15/results/quant_startnoise.json`. The **same** ORACLE arm — continuous least squares fitted
to the true distance matrix, 126 targets, 6 starts — at four independent start draws:

| seed | mean | median |
|---|---|---|
| 0 | 0.6060 | 0.0530 |
| 1 | 0.7722 | 0.0741 |
| 2 | 0.5183 | 0.0506 |
| 3 | 0.7915 | 0.0619 |

**sd of the mean across draws = 0.132 Å; range = 0.273 Å.** Per target, mean sd = 0.449 Å and
the worst single target swings by **5.15 Å**.

> The coordinator's published ORACLE constant of **0.611 Å** (K1-CORRECTED) is one draw from a
> distribution with a standard deviation of 0.132 Å; my own run of the identical arm returned
> **0.569 Å**, and the four seeds here return 0.518–0.792. The *median* is far more stable
> (0.051–0.074) because the variance lives in a few targets with many local minima. **Nothing
> is wrong with the finding — ideal-geometry torsions do reproduce a native CA trace from its
> own distance matrix — but the mean should be quoted with this spread, or the median quoted
> instead.**

The distinction that matters: this 0.132 Å is noise on an arm's **absolute level**. A **paired**
difference within one run shares its start set and cancels most of it, which is why the
headline +0.381 Å can still be sound. R7b tests exactly that.

### R7b. INDEPENDENT VERIFICATION OF THE HEADLINE — it reproduces

`s15/results/quant_ladder_repro.json`. The two arms that carry the headline, re-fitted at
independent start draws through `stable_rng`:

| run | continuous mean | snap_bin mean | **paired difference (95% CI)** | W/L |
|---|---|---|---|---|
| `quant_ladder` (original) | 0.569 | 0.950 | **+0.381 [+0.156, +0.614]** | 15/111 |
| `ladder_repro` seed 1 | 0.654 | 1.069 | **+0.416 [+0.169, +0.672]** | 19/107 |
| `ladder_repro` seed 2 | 0.619 | 1.106 | **+0.486 [+0.282, +0.695]** | 15/111 |

**Three independent start draws give +0.381, +0.416, +0.486 — all three CIs exclude zero, all
three are 107–111 losses out of 126.** The absolute continuous-arm level moves over
0.569–0.654; the paired quantisation cost moves over 0.381–0.486, i.e. a spread of 0.105 Å
around a mean of **+0.428 Å**.

This is the decomposition the sprint needs: start-draw noise is **large on levels** (sd 0.132 Å,
and the ORACLE mean is not reproducible) and **small relative to the effect on paired
differences** (spread 0.105 Å on an effect of 0.43 Å). So the headline stands, and it should be
quoted as **≈ +0.4 Å, reproduced at three independent start draws**, rather than as +0.381
alone. It also explains why the +0.081 Å false positive in R7 was possible: an effect that size
sits *below* the residual paired noise, while +0.4 Å sits well above it.

**Practical rule this yields, and I would put it in the paper's methods:** at this effect size,
a paired CI excluding zero is **not** sufficient evidence. A difference below roughly 0.1 Å on
this instrument must be reproduced at an independent start draw before it is believed, and the
cheapest sufficient control is to re-run the two arms that carry the claim with a different
seed. That costs minutes and it caught a false positive here on the first try.

---

## R8. WHAT I SET OUT TO CONFIRM AND FAILED TO — the hypothesis I was handed is REFUTED

Recorded plainly, per the brief. The hypothesis in my brief was:

> *"Maximum likelihood against a binned distribution inherits the binning as a hard resolution
> floor… If true, the fix is a continuous density fitted to the histogram that keeps its
> multimodality but restores sub-bin resolution — which should let ML finally beat least
> squares, because ML's extra information (asymmetry, bimodality) becomes usable."*

It fails in three separate places, and each failure is independent of the others:

1. **The premise was a bug.** 75% of the gap the hypothesis was built to explain came from a
   uniform-grid index map, not from binning (R0, R2).
2. **The stated mechanism does not exist in this distogram.** Strict bimodality is 8.1% of
   pairs and is *lower* at long separation than at short (R4); and the per-target ML advantage
   is uncorrelated with bimodality (rho +0.058, inside the permutation null) and with skewness
   (R6). There is no "asymmetry and bimodality" for ML to exploit.
3. **The predicted consequence does not occur.** A continuous density does recover part of the
   ORACLE quantisation cost (R2, item 6), but on predicted restraints ML still does not
   beat least squares by anything that clears the measurement floor (R3, R7).

What replaced it is better than what it proposed: the binning cost is measured for the first
time (≈ +0.43 Å to an oracle, 2.7% of the predictor's error variance), the "ML loses" record is
retracted with its cause identified, and the real discriminator between ML and least squares
turns out to be **calibration**, not shape — which was INFO's prediction, and which my data
confirm on three of four arms.

**Errors and retractions of my own, preserved in place:**

- I initially attributed the ladder's `ORACLE_ml_gauss_w0.6_continuous` gap to `scipy`'s
  **relative `ftol`** interacting with the Gaussian's additive normalisation constant, and
  built `GaussDensityNoConst` to prove it. **That explanation is REFUTED by my own control**:
  `ml_gauss_const` and `ml_gauss_noconst` differ by only −0.019 Å (R7), so the constant is not
  the cause. The correct explanation is plainer and worse — the two runs used different start
  draws, and the gap was the start draw. The refuted class is kept in the module because it is
  the control that killed my own hypothesis.
- My first read of the ladder was on 4 targets and showed a `mid_6_7` arm at 0.876 Å mean with
  a 0.213 Å median — a single outlier at n=4. On 126 the same arm is 0.765/0.207. Another
  instance of the sprint's recurring small-sample reversal; I did not act on the n=4 read.
- I let a chained background job re-run `pred` after `ladder` finished, and it overwrote the
  completed `s15/results/quant_pred.json` with a 20-target partial. The complete result was
  recovered from `s12/results/s15_quant_pred.json`, written by `I.write` — **the exact hazard
  that function's docstring warns about, saving the result rather than destroying it.** The
  redundant process was killed.

---

## R9. OPEN, AND WHAT I WOULD DO NEXT

- **A native-free proxy for per-target calibration.** R6 shows ML wins −0.13 to −0.18 Å exactly
  where `I_pred > 0`, but `I_pred` is ORACLE. Sprint 14's arithmetic says a weak proxy buys
  little (gain goes as the square of the weaker channel's skill), so this is a real but
  probably small lever — and it is testable cheaply.
- **Retraining the distogram on a uniform ≈0.5 Å grid.** R1 and R2 say the placement, not the
  count, is the loss. Nobody has to change the architecture to test it.
- **Calibration of the density itself** (temperature scaling per separation band, fitted
  leave-fold-out) is the untried intervention that R6 points at directly. It is the analogue of
  `s15/distcal.py`'s bias correction applied to the *shape* rather than the *mean*.
- **Not worth doing:** more density families, more bandwidths, or any further work on the
  output representation. The entire spread across thirteen representations is 0.33 Å and the
  honest arms sit inside the measurement floor.

---

## Reproduction

    python -m s12.instrument                 # instrument check, must reproduce the constants
    python -m s15.quant check                # derivative cosines, all representations
    python -m s15.quant mm                   # multimodality census, both bases
    python -m s15.quant budget               # analytic quantisation-noise budget
    python -m s15.quant ladder               # ORACLE representation ladder, 126 targets
    python -m s15.quant pred                 # predictive arms, 126 targets
    python -m s15.quant tolcheck             # least squares IS Gaussian ML: optimiser control
    python -m s15.quant startnoise           # multi-start draw noise, 4 seeds
    python -m s15.quant repro                # headline re-run at 2 independent start draws
    python -m s15.quant mech                 # shape-vs-calibration mechanism test

`s15/results/`: `quant_derivcheck.json`, `quant_multimodal_density.json`,
`quant_multimodal_mass.json`, `quant_snaperr.json`, `quant_noisebudget.json`,
`quant_ladder.json`, `quant_pred.json`, `quant_tolcheck.json`, `quant_startnoise.json`,
`quant_ladder_repro.json`, `quant_mechanism_*.json`.

**Provenance caveat, repeated here so it is not missed:** `quant_ladder.json` and
`quant_pred.json` were produced under the pre-fix `hash(pdb)` seeding (R4b) and will not
re-run bit-identically; every comparison *within* each file is start-matched and valid, and the
headline is separately reproduced under `stable_rng` in `quant_ladder_repro.json`. Everything
else in this file was produced after the fix.
