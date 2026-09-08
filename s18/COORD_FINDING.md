# COORDINATOR FINDING — the functional form is SOUND and the distogram's error is WORSE THAN NOISE

`s18/objceil.py`, **n = 126**, complete. Start = the coordinate average of the shipped top-75
(3.048 A). Same functional form throughout -- a weighted sum of squared Ca-pair distance residuals
-- varying **only the distances it is asked to satisfy**:

    d_alpha = (1 - alpha) * dhat + alpha * d_true

| arm | RMSD | median | vs alpha = 0 | vs the average | W/L vs avg |
|---|---|---|---|---|---|
| coordinate average (start) | **3.048** | 2.837 | — | — | — |
| **alpha = 0** (deployed distogram) | **3.610** | 3.370 | +0.000 | +0.561 [+0.410, +0.717] | 31/95 |
| alpha = 0.25 | 3.105 | 3.137 | -0.504 [-0.589, -0.417] | +0.057 [-0.109, +0.214] | 61/65 |
| **alpha = 0.5** | **2.448** | 2.509 | -1.162 [-1.334, -0.992] | **-0.600 [-0.786, -0.418]** | 91/35 |
| alpha = 0.75 | 1.789 | 1.383 | -1.821 [-2.069, -1.568] | -1.260 [-1.507, -1.021] | 104/22 |
| **alpha = 1** (perfect distances) | **1.152** | 0.307 | -2.458 [-2.770, -2.142] | **-1.896 [-2.205, -1.608]** | **114/12** |
| **shuffled** (magnitudes kept, direction destroyed) | **2.609** | 2.536 | **-1.001 [-1.226, -0.784]** | -0.439 [-0.682, -0.206] | 70/56 |
| **isotropic** (matched-RMS Gaussian noise) | **2.573** | 2.533 | **-1.037 [-1.279, -0.801]** | -0.476 [-0.746, -0.210] | 72/54 |

Distogram residual RMS **3.205 A**, MAE 2.347 A. Every alpha > 0 arm is an **ORACLE DIAGNOSTIC** --
a ceiling and a requirement, never a result.

## Two findings, and the second is the important one

**1. THE FUNCTIONAL FORM IS SOUND.** With perfect distances the *identical* objective reaches
**1.152 A** -- 80.2% of targets below 2.5 A, 73.8% below 2.0 A, beating the coordinate average on
**114 of 126**. The sum-of-squared-pair-residuals form is not the problem, and re-engineering it is
very likely a distraction.

**2. THE DISTOGRAM'S ERRORS ARE WORSE THAN RANDOM ERRORS OF THE SAME SIZE.** Preserving the residual
magnitudes and destroying only their **direction** improves the refined structure by
**-1.001 [-1.226, -0.784]**; matched-RMS isotropic noise gives **-1.037 [-1.279, -0.801]**. Both
intervals exclude zero decisively.

> The distogram is not merely inaccurate. Its error is **structured, and the structure is
> actively harmful** -- worse than having no information about the error at all. This is the
> programme's recorded *"full amplitude, wrong direction"* measured at the level of the optimum
> rather than the predictions.

## The requirement curve

    to beat the coordinate average by refinement     alpha > ~0.3
    to reach 2.5 A                                   alpha ~ 0.5   (halve the residual RMS)
    to reach 2.0 A                                   alpha ~ 0.7

## What this means for each workstream

**This does not kill the degree-1 branch -- it sharpens it, and gives it a better hypothesis.**
If the *pairwise* content of the distogram's error is the mis-directed part, then a degree-1
objective, which averages over pairs, should be **more robust to structured distogram bias**. That
is a sharper, more testable framing than "removing harmful inter-residue structure", and it makes a
specific prediction: **degree-1 should close part of the gap between alpha = 0 (3.610) and the
shuffled control (2.609)**, because both are ways of not trusting the error's direction.

The alpha ladder here is the natural companion to the lambda ladder. Read them together.


---

## ADDENDUM — both confound controls survive, and one makes the result sharper (n = 126)

| control | RMSD | vs alpha = 0 | W/L vs the average |
|---|---|---|---|
| `shuffled` (plain permutation) | 2.609 | -1.001 [-1.226, -0.784] | 70/56 |
| **`shuf_strat`** (within separation bins) | **2.560** | **-1.050 [-1.257, -0.860]** | 72/54 |
| **`shuf_paired`** (residual + its weight permuted together) | **2.072** | **-1.537 [-1.758, -1.323]** | **98/28** |
| `isotropic` | 2.573 | -1.037 [-1.273, -0.802] | 72/54 |

The separation confound is refuted (stratified shuffle is marginally *stronger*). The
weight-alignment confound is refuted and the effect nearly doubles.

`shuf_paired` preserves the residual magnitudes **and** the residual-weight pairing and randomises
only **which pair each error lands on** -- and it reaches **2.072 A**, below the sprint's primary
target, where the real distogram reaches 3.610.

> **The distogram's errors are harmful because of WHICH PAIRS they fall on** -- not their size, and
> not their relation to the model's own confidence. They concentrate on the pairs that most
> determine the structure, and the uncertainty estimates do not compensate.

**The lead:** a native-free re-weighting of the objective by **geometric leverage** now has a
measured ceiling behind it. Caution: Sprint 17 refuted 58 distogram functionals including
uncertainty weighting -- but for **selection**, not **refinement**. Different use; pre-register it
as its own experiment rather than assuming either way.

---

## AMENDMENT 2 — the `shuf_paired` row is WITHDRAWN (2026-09-06)

`s18/objceil.py` line 163 passes `sd[pi]` rather than `sd`, so the `shuf_paired` arm permutes the
residual **and** the weight vector, and is not the deployed functional. Since `shuffled` (2.609)
and `shuf_paired` (2.072) differ in exactly one thing — whether `sd` moves — the 0.537 A between
them prices **permuting the weights**, not the error assignment. Caught by the ADVERSARIAL lane
reading the source; confirmed by the coordinator at line 163. See `s18/LEDGER.md` L4.

**Withdrawn**: "it is which pairs the errors fall on, not their size nor their relation to the
model's own confidence", and every downstream use of 2.072 A as a ceiling on error assignment.

**Unaffected**: `shuffled` 2.609, `shuf_strat` 2.560, `isotropic` 2.573 all carry correct weights
and all beat the real distogram's 3.610 by ~1.0 A with intervals excluding zero. Conclusion (2) of
this document — *the distogram's errors are worse than random errors of the same magnitude* — stands
as measured. Conclusion (1) — *the functional form is sound* — is untouched by this and is now under
attack 3 (objective-limited vs optimiser-limited) in the adversarial lane.

**What replaces it**: the isolating arms `wperm_only` (real dhat, permuted weights) and `wflat`
(real dhat, uniform weights). Both are **native-free**, so whatever they show is a deployable
result rather than an oracle ceiling.
