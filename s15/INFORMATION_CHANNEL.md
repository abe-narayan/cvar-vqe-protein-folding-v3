# THE INFORMATION CHANNEL

What information is available for predicting a 9–16 residue backbone, how much of it there is, and —
the part that turned out to matter most — what **shape** it comes in.

The central result of this document is that **error magnitude does not price a channel**. Two
channels in the same pipeline, measured the same way, have structured errors; the structure makes one
channel **2.9 Å cheaper** than its magnitude implies and the other **1.1 Å more expensive**. Every
accuracy budget this project computed before measuring that was computed on the wrong axis.

---

## 1. THE CHANNELS, AND WHAT EACH IS WORTH

126 real retrieval pools, in-band defined as `pool_best + 1.5 Å`, ties averaged. "Required for 2.0 Å
through a top-100 operator" is **0.638**; chance is 0.500.

| channel | availability | in-band accuracy | CI95 | emitted @ 24 |
|---|---|---|---|---|
| **distogram Bayes risk, low-`sd` half** | 126/126 | **0.566** | [0.547, 0.586] | 3.536 |
| distogram L2, low-`sd` quartile / half / 1/sd² | 126/126 | 0.560 / 0.559 / 0.548 | — | 3.65 / 3.55 / 3.52 |
| distogram Bayes risk (the incumbent's filter) | 126/126 | 0.545 | [0.522, 0.567] | **3.501** |
| ESM-2 contact map | 126/126 | 0.513 | [0.494, 0.531] | 4.384 |
| pool-consensus secondary structure | 126/126 | 0.513 | [0.498, 0.527] | 4.070 |
| radius of gyration | 126/126 | 0.497 | [0.475, 0.517] | 4.566 |
| pool typicality | 126/126 | 0.495 | [0.465, 0.526] | 3.723 |
| random null | — | 0.498 | [0.490, 0.507] | 4.492 |

**The best channel available reaches 0.566 where 0.638 is required.** No channel in this table, or
any combination of them, closes that gap.

### 1.1 Information content, measured properly

Using **held-out predictive cross-entropy**, not plug-in mutual information:

| channel | bits | positive on |
|---|---|---|
| retrieval pool → native torsion | **−0.11 to −0.23 bits/torsion** vs a generic Ramachandran prior | — |
| distogram → native distance | **−0.818 bits/pair** [−1.265, −0.419] | 52/126 targets |
| ESM-2 contacts → native contact | **−2.570 bits/pair** | **0/126** targets |

Two things follow. **The retrieval torsion prior carries no target-specific advantage** — it *is*
target-specific (a target-shuffled null sits at −0.70 to −0.80) but it is not better than a generic
Ramachandran distribution, and two independent cross-checks agree: modal hit rate 0.371 against
0.360, and a structural twin giving 4.072 Å against a **4.065 Å constant α-helix**. And **the
distogram is a good ranker and a bad probability** — negative cross-entropy on 74 of 126 targets,
which is why maximum likelihood under its own distribution loses to least squares on its mean.

> **A measured methodological warning.** The plug-in mutual-information estimator reports **+0.40 to
> +0.88 bits** on the same data where held-out cross-entropy reports **−0.11 to −0.23**. The sign is
> wrong, not just the magnitude. Plug-in MI is unusable at this sample size and this project should
> never quote it again.

---

## 2. THE CENTRAL RESULT: STRUCTURE, NOT SIZE

### 2.1 The torsion channel — real errors are CHEAPER than i.i.d. (0.6–1.7 Å native-free; 2.9 Å ORACLE)

A phase diagram over torsion-channel accuracy, with i.i.d. errors of controlled size, says 2.0 Å
requires **16.3° RMS at full coverage** and the best available channel is **67.6° RMS**. On that
surface the problem is hopeless.

The surface is wrong. The i.i.d. surface **over-prices real channels by 0.6–1.7 Å across
native-free arms**, and by up to **2.9 Å** for the **ORACLE best-matching pool window** — which is the
arm that emits **1.77 Å** where i.i.d. of the same magnitude emits **4.71 Å**. The largest
**native-free** discount is the incumbent's projected torsions at **1.7 Å**.

**The mechanism is subspace alignment**, measured rather than assumed: `||Je|| / (||e||·s_rms)` =
**0.565** for the ORACLE window and **0.665** for the incumbent's projection, against a **0.945**
random-direction null.

> **And the null is not the control that matters.** A **zero-information constant α-helix already
> scores 0.709.** Alignment below 1 is therefore not by itself evidence of a good channel — it must
> clear the helix, and only those two arms do. This is the project's own recurring trap (a constant
> helix beats the random control on this instrument) appearing on a new statistic.

**Causality is established by surrogate destruction — for the ORACLE arm.** Sign-flipping its errors,
preserving every magnitude exactly and changing only directions, costs **+1.596 Å [+1.385, +1.813]**.
**For the native-free retrieval mean the same operation is null** (−0.150 [−0.366, +0.070]), so the
causal demonstration currently rests on an ORACLE arm and must be stated that way.

Two corollaries that close off the natural responses:

- **This is invisible in pairwise correlations** (all |r| ≤ 0.16). Sprint 14's finding that real
  torsion predictors have near-i.i.d. errors is confirmed *and non-informative*: lag-1 correlation
  says nothing about Jacobian alignment.
- **Correlated is not aligned.** Synthetic coherence — AR(1), class bias — moves the result the
  opposite way. The measured coherence factor is 1.24, not the 2 that had been assumed.

### 2.2 The distance channel — real errors are 1.1 Å MORE EXPENSIVE than i.i.d.

The same experiment on distances gives the opposite sign. Corrupting **true** distances with i.i.d.
Gaussian noise and refitting:

| σ (Å) | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|
| emitted RMSD | 0.554 | 0.709 | 0.939 | 1.306 | 1.450 | 1.787 | 2.075 | **2.561** |

The distogram's MAE of 2.386 Å is a Gaussian σ ≈ 2.99, where the surface says **2.56 Å**. The real
distogram emits **3.644 Å**.

Surrogate destruction identifies the expensive property (ORACLE diagnostics; destroying a property
requires knowing it):

| surrogate | mean | vs real |
|---|---|---|
| `real` (native-free) | 2.832 | — |
| **`ORACLE_shuffle_signs`** — magnitudes kept, signs randomised | **1.825** | **−1.007 [−1.314, −0.664]** |
| `ORACLE_gauss_matched` — i.i.d. at the target's own RMS | 2.501 | −0.332 |
| `ORACLE_shuffle_pairs` — permuted across pairs | 2.577 | −0.256 |
| `ORACLE_winsor_z3` — outlier tail clipped | 2.711 | −0.121 |
| `debias_sep` (native-free) | 3.012 | **+0.179** |

**The sign pattern carries the whole effect.** Outliers are worth 0.12; separation structure 0.26;
Gaussian matching 0.33. Randomising the directions while keeping every magnitude is worth **an
ångström**.

### 2.3 What the sign coherence actually is: the errors are REALIZABLE

Fitting to the predicted distances and comparing the converged fit's own residual against the
prediction's true error, on 24 targets:

| quantity | value |
|---|---|
| mean \|d_fit − d̂\| — the **unrealizable** part of the prediction | **0.879 Å** |
| mean \|d̂ − d_true\| — how **wrong** the prediction is | **2.009 Å** |
| ratio | **0.437** |

> **About 56% of the distogram's error is geometrically realizable. The predicted distance matrix is
> not a noisy version of the right structure — it is a good description of a WRONG structure.**

**The null control, which is what makes this sound.** A residual of 0.879 Å means nothing without
knowing what the *same* fit leaves on a distance set of the same magnitude with no structure — the
fit is constrained (`2n − 4` degrees of freedom against `n(n−1)/2 − (n−1)` pairs), so some residual
is guaranteed. On the same targets:

| distance set handed to the fit | unrealizable part |
|---|---|
| **the REAL predicted distances** | **0.879 Å** |
| ORACLE i.i.d. Gaussian of the same RMS — **the null** | **2.058 Å** |
| ORACLE the same errors permuted across pairs | 1.824 Å |
| ORACLE the same magnitudes with random signs | 1.790 Å |

**The real prediction is 2.3× closer to realizable than matched-magnitude noise**, and both
structure-destroying surrogates move it back to within 13% of the null. The low residual is a
property of the prediction, not of the fit's flexibility.

And the negative correlation completes it: the part of the error the fit *cannot* satisfy is
anti-correlated with the truth, so **the optimiser discards the incoherent component and faithfully
follows the coherent one**.

### 2.4 The unified statement

Both channels have structured errors. The structure lies in the Jacobian's **low-response** subspace
for torsions, so the errors are cheap; it lies in the **realizable** subspace for distances, so the
errors are followed. Same phenomenon, opposite sign, one explanation: **what a channel costs is set
by where its error points relative to the map that consumes it.**

This is why every correction attempted in this project has under-delivered. Reweighting,
recalibration, separation debiasing, robust and redescending losses, outlier trimming, per-target
rescaling — **all of them operate on the error's magnitude profile, and the expensive part is a
direction.** Five confirmations, each independently derived: recalibration is provably worth zero;
debiasing is worth little — negative on ranking, at the noise floor on the fit; outlier clipping is worth 0.12 Å; per-target
rescaling has a 0.1 Å ceiling; consensus rescaling has a 0.275 Å ceiling that the distogram
overshoots.

---

## 3. WHERE THE GAPS FALL BEATS HOW MANY THERE ARE

At σ = 0 and coverage 0.7, the emitted structure is **1.716 Å when the gaps are terminal** and
**3.505 Å when they are mid-chain** — a **1.79 Å spread**, three to four times larger than previously
recorded.

**Terminal-gap-shaped coverage is the largest single lever in the phase diagram, and no channel
exploits it.**

---

## 4. FUSION, AND A LAW FOR WHAT IT IS WORTH

The realizability result reframes what a second channel is *for*. Fusion has been justified in this
project by error decorrelation and dismissed by a squared-skill argument. Neither is the right
account.

**Two channels can be correlated in their errors and still disagree about the structure.** Fitting
the distogram and the retrieval pool separately, on 24 targets:

| | |
|---|---|
| distogram fit vs native | 3.622 Å |
| pool fit vs native | 3.659 Å |
| **the two fits disagree with each other by** | **3.138 Å** |
| raw per-pair error correlation | **+0.602** |
| the pool's unrealizable part | 0.458 Å against the distogram's 0.977 — *more* realizable |

Correlated errors, **different wrong structures**. So the useful question is not whether two channels
err on the same pairs but whether their realizable components point at the same wrong structure.

**The law.** For two fits at mean distance `r` from the native and `s` from each other, the midpoint
sits at

> ~~**d_avg ≈ √(r² − (s/2)²)**~~ → **d_avg = √( (r₁²+r₂²)/2 − (s/2)² )**, exact

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT.** This is the two-member **Krogh–Vedelsby (1995)
> ambiguity decomposition** — prior art, cite it. `r` below is the **arithmetic** mean where the
> identity needs the **quadratic** mean; with the correct mean the residual is **+0.075 Å**, not
> +0.162 Å (paired **−0.087 [−0.119, −0.061]** i.i.d.-target, **[−0.129, −0.056]** fold-clustered,
> W/L 126/0, n = 126, target as the unit). In a single common frame with the correct mean the
> identity is **exact to 2.65e−15** on all 126 targets, so what remains is a frame convention, not a
> prediction error. `s16/retract_law.py`, `s16/retract_exact.py`, `s16/retract_FINDINGS.md`.

*(all 126 targets)*

| | ~~published~~ | **corrected** |
|---|---|---|
| `r` (arithmetic) / `r` (quadratic, per target then averaged) | ~~3.641 Å~~ | — |
| `s` (Kabsch-minimised) / `s` in the common frame | 3.138 Å | 3.506 Å |
| predicted average | ~~**3.205 Å**~~ | **3.292 Å** |
| measured coordinate average (medoid frame) | **3.367 Å** | 3.367 Å |
| **prediction error** | ~~**+0.162 Å**~~ | **+0.075 Å** (median +0.019) |
| residual in one common frame, quadratic mean | — | **2.65e−15 Å** |
| restraint-level fusion | 3.551 Å | 3.551 Å |

Coordinate averaging beats the better single channel by **−0.255 [−0.370, −0.145]**, and beats
restraint-level fusion — reversing a 24-target read that had said the opposite. Combining at the
**structure** level beats combining at the **restraint** level, which is consistent with
realizability: a single fit to blended restraints must pick *one* coherent wrong structure, while an
average of two separately-fitted structures can land between two different ones.

~~**`s` is native-free**, so this converts "should we fuse these channels?" into a calculation on
unlabelled data.~~ **⚠ REFUTED, Sprint 16, 2026-09-06 (RETRACT).** The prediction needs `r`, a
distance **to the native**, so it is not a calculation on unlabelled data; and `s` on its own ranks
"does fusion win on this target" at **AUC 0.401 against a 0.500 null** — it points the wrong way. A
supervised leave-fold-out threshold on `s` chose "always fuse" on 5/5 folds (`s16/retract_screen.py`).
The scaling statement that follows is algebra and is correct: the gain is small whenever
`s ≪ r` and only becomes large as `s → 2r`. Here `s/r = 0.88`, worth about 0.2 Å. **An ångström from
fusion would require channels that disagree about twice as much as they are wrong**, and no pair in
this project comes close.

---

## 5. THE `sd` COLUMN: A REAL DISCRIMINATOR THE OPERATOR DOES NOT CONSUME

Spearman(`sd`, |error|) = **+0.540**, and RMS error rises **1.34 → 5.93 Å** across `sd` quintiles.
The column had never been used for anything in fourteen sprints.

Conditioning the shipped score on it gives **+0.021 [+0.010, +0.034]** in-band accuracy, 74 wins to
51, concentration check PASS — **the best in-band number this project has measured** — and
**+0.033 [−0.051, +0.116] *worse* emitted RMSD.**

A channel can be genuinely informative and not survive the operator that consumes it. Note this does
**not** contradict the generative result that `sd` weighting is worth 0.154 Å in the fit: selection
and generation are different axes, and the same column pays off on one and not the other. That is the
third independent instance of the pattern in this sprint, after the pool channel and separation
debiasing.

---

## 6. WHAT THE PHASE DIAGRAM SAYS ABOUT THE TARGETS

| target | verdict | requirement |
|---|---|---|
| 3.0 Å | **supported** | coverage 1.0 at ≤ 26°, or 0.8 at ≤ 13°, or **0.5 with terminal gaps at ≤ 11°** |
| 2.5 Å | **marginal** | coverage ≥ 0.9 at ≤ 15°, or ≥ 0.6 with terminal gaps |
| 2.0 Å | **not supported** by any i.i.d.-class channel | 16.3° at full coverage; best channel 67.6° |

Verified independently: a constant-helix control reproduces 4.065 Å exactly; the incumbent's
σ-equivalent measures 28.6° against a recorded 27.1°; the 2.0 Å requirement measures 16.3° against a
recorded 15.1°.

**These verdicts are computed on the i.i.d. surface, which §2 shows over-prices real channels.** They
are therefore *conservative* for the torsion channel and *optimistic* for the distance channel — and
restating them with alignment accounted for is the one open question that could move 2.0 Å from
unsupported to reachable without any improvement in raw predictor accuracy.

---

## 7. REFUTATIONS, INCLUDING THE WORKSTREAM'S OWN HEADLINE

- **A native-free per-target sign proxy picked the sign at 0.947 on the enumerated band and captured
  99.7% of an ORACLE channel worth −0.339 Å at 19 wins to 0 — and did not transfer.** On the 126 real
  pools the entire ORACLE channel is worth 0.038–0.083 Å, the proxy recovers 0.020 Å [−0.051, +0.010]
  at permutation p = 0.067, and **outside the band it is harmful, +0.43 to +0.68 Å**. Tiered REFUTED.
- **The +0.909 skill figure was a near-tautology**: ρ(Rg, RMSD) against native-Rg z is −0.935
  (p = 0.000), which is definitional. On the actual leave-fold-out model the same quantity is +0.154,
  p = 0.531. The two skill statistics **anti-correlate** (p = 0.032), so "the sign" was never one
  quantity, and against emitted RMSD every compactness variant gives |ρ| ≤ 0.12 at p ≥ 0.63.
- **Per-torsion σ does not price a generation channel** (§2.1).
- **Near-i.i.d. lag-1 does not license i.i.d. modelling** (§2.1).
- **ESM contact maps are unusable at peptide length** — 0/126 targets positive.
- **Whole windows beat per-residue recombination by −0.457 Å on the mean and lose by +0.245 Å on
  best-of-500.** Recombination helps only if the selector is sharp, and nothing here ranks within the
  pool.
