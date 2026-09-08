# SPRINT 19 — AGENT A (distogram / structure prediction)

**Question owned.** Why does the sequence-conditioned distogram produce errors that are
systematically worse than random errors of the same magnitude, and what change to the predictor,
its loss or its output representation moves the error somewhere useful?

Pre-registration: `s19/PREREG_A.md`, written before any arm ran and not edited since.
All artefacts under `s19/results/` with `.COMPLETE` flags. Code: `s19/a_lib.py`, `a_fit.py`,
`a_topo.py`, `a_coh.py`, `a_struct.py`, `a_models.py`, `a_source.py`, `a_subset.py`,
`a_sepfix.py`, `a_audit.py`.

**Instrument.** 126 cluster-disjoint tuning targets, pinned folds, `s12/instrument.py`.
The fit harness (`s19/a_fit.py`) is byte-identical to `s18/objceil.py`'s: start = ideal-geometry
projection of the coordinate average of the shipped top-75; objective = the deployed
`Σ w_ij (d_ij(φ,ψ) − d̂_ij)²` with `w = 1/sd²`; metric = full-chain Cα-RMSD, model 1 of the native.
**Reproduction gate passes exactly on every module: `avg` 3.048, `real`/`a0.0` 3.610.**
The leave-fold-out separation debias is refitted on the FULL 126-target list in every module
(the subset trap, brief §7). BLAS threads capped at 1 throughout.

---

# 0. THE HEADLINE, AND WHAT IT COSTS

> **The Sprint-18 effect is the residual's cross-pair SIGN CORRELATION, and nothing else.**
> Destroying only that — every pair keeping its own magnitude and its own weight — is worth
> **−1.202 Å [−1.408, −1.004], 112W/14L, 5/5 folds**, at **exactly matched residual RMS**
> (§1.0), which is **128 % of the whole 3.610 → 2.671 gap**.
> A **perfectly coherent** error of matched magnitude is the **worst arm on the board**.

And the price, stated first because it damages my own lane hardest:

> **Not one of eight real predictors — including two that cut ORACLE MAE by 0.24–0.27 Å, a
> different architecture with triangle updates, and two ensembles — moves Cα-RMSD by more than
> the 0.084 Å MDE with a CI excluding zero.** The best is `combo` at −0.094 [−0.191, **+0.001**].
> **Accuracy is not the lever, and I have no native-free estimator of coherence.**

The one thing in my lane that *did* move with a CI excluding zero is a **negative**: weighting the
training loss toward long-range pairs costs **+0.181 Å [+0.081, +0.285]** while leaving ORACLE MAE
unchanged (+0.010 [−0.039, +0.062]).

---

# 1. Q2 — IS THE ERROR COHERENT?  **ESTABLISHED. This is the mechanism.**

`s19/a_coh.py` · `s19/results/a_coh.json` (COMPLETE, n=126) · **every arm is ORACLE**

| arm | RMSD | median | resid RMS | κ | EDM defect | tri % | vs `real` |
|---|---|---|---|---|---|---|---|
| `real` (deployed d̂) | 3.610 | 3.370 | 3.205 | 0.809 | 0.286 | 4.09 | — |
| **`signflip`** † | **2.368** | 2.224 | 3.062 | 0.352 | 0.340 | 10.85 | **−1.242 [−1.455, −1.044]** 113W/13L |
| `shuffled` (S18 control) | 2.649 | 2.606 | 3.008 | 0.367 | 0.361 | 11.45 | −0.961 [−1.184, −0.740] 100W/26L |
| `iso` | 2.655 | 2.622 | 3.006 | 0.367 | 0.368 | 12.02 | −0.955 [−1.190, −0.733] 94W/32L |
| **`coherent`** | **3.749** | 3.575 | 3.097 | 0.990 | 0.002 | 0.02 | **+0.139 [+0.023, +0.253]** 56W/70L |
| **`coh_scaled`** | **3.843** | 3.632 | 3.205 | 0.987 | 0.010 | 0.57 | **+0.233 [+0.122, +0.348]** 46W/80L |

† **Superseded — read §1.0 before quoting this row.** `signflip`'s residual RMS is 3.062, not
`real`'s 3.205, because the 2.0 Å clip shrinks a sign-flipped field. The magnitude-exact value is
**−1.202 [−1.408, −1.004]**. (CIs differ in the third digit between modules only through the
bootstrap seed; the realisations are identical.)

**What `signflip` is, read from the source and not from this sentence.** `dtrue + r·ε`, ε = ±1
i.i.d. per pair. Every pair keeps its own residual magnitude and **its own weight** — no `sd[pi]`
orphaning, the exact defect that retracted Sprint-18's G5 — and the spatial pattern of `|r|` is
untouched. Only the cross-pair sign correlation is destroyed.

**What `coherent` is.** The residual is replaced by the residual of a **real alternative
structure**: the top-75 pool member whose residual RMS is closest to the real one. It is a
genuine Euclidean distance matrix — perfectly realisable — at matched magnitude.

### 1.0 **SELF-AUDIT: a defect in my own headline, found and corrected**

`s19/a_audit.py` · `s19/results/a_audit.json` (COMPLETE, n=126).

Every field is clipped at 2.0 Å before the fit, exactly as `s18/objceil.py` does. **Flipping a
residual's sign can push a distance below 2.0, so the clip shrinks the whitened field**:
`signflip` carries residual RMS 3.062 against `real`'s 3.205, a 4.5 % reduction. On the α ladder
a 25 % magnitude cut is worth 0.504 Å, so 4.5 % extrapolates to ≈ 0.091 Å — about **7 % of the
1.242 Å headline**. Small, but this is the sprint's central result.

`*_exact` rescales the field *after* the clip to a fixed point where the realised residual RMS
equals `real`'s to machine precision (3.2049 for every arm):

| arm | RMSD | median | resid RMS | vs `real` | W/L | folds |
|---|---|---|---|---|---|---|
| `real` | 3.610 | 3.370 | 3.2049 | — | | |
| `signflip` (as published) | 2.368 | 2.224 | 3.0621 | −1.242 [−1.458, −1.038] | 113W/13L | 5/5 |
| **`signflip_exact`** | **2.408** | 2.238 | **3.2049** | **−1.202 [−1.408, −1.004]** | 112W/14L | 5/5 |
| `shuffled` | 2.649 | 2.606 | 3.0076 | −0.961 [−1.194, −0.747] | 100W/26L | 5/5 |
| **`shuffled_exact`** | 2.671 | 2.576 | **3.2049** | −0.938 [−1.169, −0.717] | 99W/27L | 5/5 |
| `iso_exact` | 2.733 | 2.718 | 3.2049 | −0.877 [−1.106, −0.652] | 93W/33L | 5/5 |

**`signflip_exact − shuffled_exact` = −0.264 [−0.448, −0.084], 83W/43L, 5/5 folds.**

> **The headline survives at exactly matched magnitude.** The clip was worth 0.040 Å of the
> 1.242, i.e. 3 %, not the 7 % the extrapolation suggested. **Quote −1.202 [−1.408, −1.004] and
> 128 % of the gap.** The published −1.242 / 129 % figures are superseded by 0.040 Å; every
> qualitative conclusion is unchanged.
>
> Note the confound pointed **against** `signflip` in the published `signflip − shuffled`
> comparison: `shuffled` and `iso` were clipped *harder* (3.008, 3.006 — a 6.2 % reduction), so
> they carried *less* error than `signflip` did. Correcting all three leaves the gap at −0.264
> where it was −0.282.

### 1.1 The three things this establishes

1. **Sign coherence is the whole effect and more.** At exactly matched magnitude (§1.0):
   `real − shuffled_exact` = +0.938; `real − signflip_exact` = **+1.202**. `signflip` beats
   `shuffled` by **−0.264 [−0.448, −0.084], 5/5 folds**. So "which pair carries which magnitude"
   is worth 0.26 Å in the *opposite* direction — it is a second-order **cost** of the Sprint-18
   control, not part of the effect.
2. **A positive control points the same way.** A perfectly coherent error of the same size is
   *worse* than the predictor's actual error, with a CI excluding zero on both variants.
3. **At matched magnitude the whole family is monotone in coherence.** Every arm sits at residual
   RMS 3.01–3.21. κ (the fraction of the predicted deviation that some ideal-geometry structure
   *can* realise) orders the outcome: κ 0.35–0.37 → 2.37–2.65 Å; κ 0.809 → 3.610; κ 0.99 →
   3.75–3.84. **Span 1.475 Å at constant error size.** The deployed distogram sits at κ = 0.809:
   four fifths of its error is realisable by a wrong structure.

### 1.2 **A4 is REFUTED, in the direction nobody expected**

`s19/a_topo.py` measured metric realisability on its own terms, as the brief asks.
The predicted field really is badly non-Euclidean:

| quantity | predicted | native (ORACLE) |
|---|---|---|
| triangle-inequality violations | **4.09 %** | 0.02 % |
| rank-3 EDM defect of −½·J·D²·J | **0.286** | 0.0018 |
| non-Euclidean mass | 0.183 | 0.0009 |

**But that is a consequence of error magnitude, not a mechanism.** A magnitude-matched
*incoherent* field is worse on both (0.340 defect, 10.85 % violations) and lands **1.24 Å
better**. Across the arm family, **unrealisability is anti-correlated with harm.**

> **"Independent pairwise marginals are not jointly realisable" is real, measured, and NOT the
> mechanism.** Nobody should spend on EDM projection, triangle repair, or joint-consistency
> enforcement as a route to RMSD. §4.2 confirms this at the architecture level.

*(ρ(EDM defect, RMSD) = +0.511 across targets — but ρ(residual RMS, RMSD) = +0.891, so the defect
adds nothing over error size. It is a magnitude proxy.)*

---

# 2. Q1 — PAIRWISE ERROR TOPOLOGY  **A1 CONFIRMED, and it is NOT the mechanism**

`s19/a_topo.py` · `s19/results/a_topo.json` (COMPLETE) · n = 8,549 pairs over 126 targets.
Feature axis native-free; error axis ORACLE.

### 2.1 The error is anti-correlated with the Sprint-18 utility compass

| | ρ(|r|, x) | pre-registered prediction |
|---|---|---|
| sequence separation | **+0.466** | > +0.30 ✓ |
| confidence 1/sd² | **−0.520** | < −0.30 ✓ (D9 recorded −0.476 — reproduced) |

| stratum | n | MAE | RMS | share of squared error |
|---|---|---|---|---|
| sep 2–3 | 2636 | 0.976 | 1.462 | 0.049 |
| sep 4–5 | 2132 | 2.157 | 2.955 | 0.161 |
| sep 6–7 | 1628 | 2.824 | 3.971 | 0.221 |
| sep 8–10 | 1506 | 3.760 | 5.190 | 0.350 |
| sep 11+ | 647 | 4.661 | 6.272 | 0.219 |

| confidence quartile | MAE | bias | share of squared error |
|---|---|---|---|
| Q1 (most confident, sd < 0.67) | 0.876 | −0.422 | **0.042** |
| Q2 | 1.775 | −0.346 | 0.127 |
| Q3 | 2.692 | −0.113 | 0.267 |
| Q4 (sd > 2.00) | 4.225 | +0.926 | **0.564** |

> **The high-utility cell the Sprint-18 compass names — short (sep ≤ 5) AND confident — holds
> 39.6 % of the pairs and 9.4 % of the squared error.** The predictor's error lives almost
> entirely where improving it converts worst (`fix_long` +0.165, `fix_unconfident` +0.244).

**This is a real anti-alignment and it is NOT the Sprint-18 mechanism.** A1 is about *where the
error is*; §1 is about *what shape it has*. §1 shows shape is worth 1.24 Å at constant magnitude,
and I refuse to let a confirmed A1 be read as the mechanism.

### 2.2 Two further strata worth recording

* **Terminal pairs dominate.** min(i, n−1−j) = 0 is 30.8 % of pairs and **45.5 % of the squared
  error** (MAE 3.039 against 1.619 for the core, min ≥ 3). The metric is full-chain including
  termini, by rule, so this is not trimmable.
* **A large conditional bias survives the deployed debias, and it is a shrinkage signature.**
  The separation debias zeroes bias *by separation bin* (−0.001 to +0.045 above). Conditioned on
  the **predicted distance regime** it does not: predicted contacts (d̂ < 8 Å) carry bias
  **−0.743**, non-contacts **+0.456**. The predictions are *under-dispersed* — regression to the
  mean, as a probabilistic mean readout must be. **This is not a new lever:** it is exactly what
  `s15/distcal.fit_correction(..., "sep_affine")` fits, and post-hoc distogram correction is a
  closed branch (calibration slope +0.376, corr +0.283). Recorded as topology, not as a fix.

### 2.3 A native-free model can rank the error, moderately

A linear model of `|r|` on 19 native-free features reaches R² = 0.282 and **Spearman +0.537**
against the ORACLE error. Largest standardised coefficients: `d̂` +0.844, `sep` +0.623,
`ent` +0.597, `maxp` +0.459, `sd` +0.403. So *detecting* the error is partly possible — §5 asks
whether detecting it is worth anything.

---

# 3. Q4a — THE PREDICTOR ZOO  **NULL, and it damages my own hypothesis**

`s19/a_models.py` · `s19/results/a_models.json` (COMPLETE, n=126). Eight leave-fold-out
predictors, all post-fold-repin, all 183 features, all with a **per-family** leave-fold-out
separation debias refitted on the full target list. Every arm native-free; MAE and κ are labelled
ORACLE diagnostics.

| model | RMSD | median | MAE\* | κ\* | EDM def | tri % | vs deployed |
|---|---|---|---|---|---|---|---|
| `deployed` (incumbent) | 3.610 | 3.370 | 2.347 | 0.809 | 0.286 | 4.09 | — |
| `d20` (dropout 0.2) | 3.574 | 3.490 | 2.115 | 0.873 | 0.201 | 1.36 | −0.036 [−0.159, +0.090] 71W/55L f3/5 |
| `d20_s1` (seed 1) | 3.608 | 3.579 | 2.105 | 0.872 | 0.199 | 1.37 | −0.001 [−0.121, +0.119] 63W/63L |
| `sw_none` (uniform loss) | 3.508 | 3.332 | 2.117 | 0.866 | 0.203 | 1.43 | −0.101 [−0.239, +0.032] 74W/52L f4/5 |
| `sw_lin` (long-weighted loss) | 3.689 | 3.596 | 2.127 | 0.868 | 0.202 | 1.25 | +0.080 [−0.051, +0.206] 58W/68L |
| `pairnet` (triangle + axial) | 3.529 | 3.506 | 2.104 | **0.914** | **0.148** | **0.68** | −0.080 [−0.242, +0.082] 70W/56L |
| `d20_ens` (2-seed ensemble) | 3.592 | 3.506 | **2.073** | 0.881 | 0.178 | 0.99 | −0.018 [−0.141, +0.108] |
| `combo` (MLP + PairNet) | 3.516 | 3.374 | 2.082 | 0.876 | 0.210 | 1.41 | −0.094 [−0.191, **+0.001**] 78W/48L f4/5 |

\* ORACLE diagnostics.

### 3.1 **A 0.27 Å cut in MAE is worth ≤ 0.10 Å in RMSD, and every CI spans zero**

`deployed` MAE 2.347 → `d20_ens` 2.073, a **11.7 % reduction**, buys −0.018 [−0.141, +0.108].
`combo` cuts MAE by 0.265 and buys −0.094 [−0.191, **+0.001**] — the CI touches zero and the
point estimate is above the 0.084 MDE, so it is **NOT MEASURED, never "an improvement"**. Under
brief §9 no arm here is a validated improvement.

> **"Better matrix, worse ranking" reproduces on the objective path.** A lower MAE that does not
> lower Cα-RMSD is a negative result, and eight of eight are negative.

### 3.2 **Joint consistency works as advertised at the geometry level and does not convert**

`pairnet` — the triangle multiplicative update, the pre-registered "joint consistency
enforcement" candidate — does exactly what its architecture claims: it drives triangle violations
from **4.09 % to 0.68 %**, EDM defect from 0.286 to **0.148**, and κ up to **0.914**, the highest
of any model. And it buys **−0.080 [−0.242, +0.082]**, 70W/56L. This is the architecture-level
confirmation of §1.2: **making the field more realisable makes it more coherently wrong.**

**Candidate 3 of the brief's item-4 list ("joint consistency enforcement") is CLOSED.**

### 3.3 **The pre-registered matched pair: the training loss's separation weight**

`sw_lin` and `sw_none` are the same architecture, the same data, the same regularisation, and
differ **only** in `separation_weights`: `sw_lin` weights the loss `min(sep/3, 6)` — up to 6×
toward long-range pairs, the direction `core/predict.py`'s own comment argues for ("an unweighted
loss spends capacity where the model is already right") and the direction the Sprint-18 compass
says converts **worse than uniform**.

| | value | CI | W/L | folds |
|---|---|---|---|---|
| **Cα-RMSD**, `sw_lin − sw_none` | **+0.181** | **[+0.081, +0.285]** | 49W/77L | 4/5 |
| ORACLE MAE, `sw_lin − sw_none` | +0.010 | [−0.039, +0.062] | 71W/55L | — |

> **Weighting the training loss toward long-range pairs costs 0.181 Å of Cα-RMSD while changing
> ORACLE MAE by nothing measurable.** Above the MDE, CI excluding zero, 4/5 folds. This is the
> Sprint-18 compass reproduced *at the predictor*, in the anti-utility direction, and it is the
> cleanest dissociation of MAE from RMSD in my lane.

It is also a **caution about the incumbent**: `sw_none` (3.508) is the numerically best single
MLP on the board, and it is the one whose loss is *uniform*.

### 3.4 A sign that must not be misread

Across the eight *real* models, ρ(κ, RMSD) = **−0.333** — higher κ goes with *lower* RMSD, the
opposite of §1's within-family ordering. There is no contradiction: §1 varies coherence at
**constant** magnitude; across real models κ and MAE move together (the better models are both
more accurate and more self-consistent), so the across-model correlation is confounded by
magnitude. **n = 8 models, descriptive, not an inferential result.** I flag it because it is
exactly the kind of statistic this programme has misread before.

---

# 4. SOURCE ATTRIBUTION — **the coherent error is SHARED, not architectural**

`s19/a_source.py` · `s19/results/a_source.json` (COMPLETE, n=126) · all ORACLE.
Run at the coordinator's direction, ahead of the mode decomposition, because it gates whether
any predictor intervention is worth money.

**The decomposition.** Each family's error is split *exactly*, using that family's own terminal
fit `X_F`:

    r_F = d̂_F − d_true = r_coh_F + r_inc_F ,   r_coh_F = d(X_F) − d_true ,  r_inc_F = d̂_F − d(X_F)

`r_coh_F` is a real structure's distance error, so cross-family alignment of `r_coh` is the
question. **This is an identity, not a finding.** The null for "shared" is two independently
seeded members of the *same* architecture — never zero.

| pair | raw r | **coherent** | incoherent | RMSD(X_A, X_B) | what it is |
|---|---|---|---|---|---|
| `d20 ~ d20_s1` | 0.865 | **0.898** | 0.400 | 1.763 | same architecture, different seed — **the ceiling** |
| `deployed ~ d20` | 0.737 | 0.833 | 0.237 | 2.368 | same architecture, different regularisation |
| `d20 ~ pairnet` | 0.732 | 0.754 | 0.115 | 2.609 | **different architecture**, matched regularisation |
| `deployed ~ pairnet` | 0.630 | 0.731 | 0.097 | 2.690 | **different architecture**, same inputs |
| `deployed ~ poolmean` | 0.709 | **0.784** | 0.090 | 2.522 | **not a network at all** — retrieval |
| `pairnet ~ poolmean` | 0.783 | 0.759 | 0.185 | 2.552 | |
| `deployed ~ sepprior` | 0.574 | **0.575** | −0.001 | 3.936 | **ZERO-INFORMATION** sequence-blind sep prior |
| `deployed ~ helix` | 0.507 | **0.598** | 0.017 | 3.724 | **ZERO-INFORMATION** constant ideal α-helix |
| `deployed ~ its own signflip` | | | | 3.474 | **the structural null** |

As a share of the same-architecture ceiling (0.898): cross-architecture **0.81**, retrieval
**0.87**, sequence-blind separation prior **0.64**, constant α-helix **0.67**.

**The incoherent component does not share at all** (0.09–0.19 cross-family, ≈0.00 against
zero-information). The sharing is specifically in the part that hurts.

### Surgery, magnitude-matched (whitening bound here: `signflip` −1.257)

| arm | RMSD | vs `real` | W/L |
|---|---|---|---|
| `minus_shared_m` — PairNet-aligned component removed | 2.830 | −0.780 [−0.960, −0.612] | 102W/24L |
| `minus_pool_m` — retrieval-aligned removed | **2.623** | **−0.987 [−1.177, −0.801]** | 109W/17L |
| `minus_sep_m` — zero-info-aligned removed | 2.858 | −0.752 [−0.934, −0.578] | 96W/30L |
| `only_shared_m` — **only** the shared component | 3.842 | **+0.232 [+0.114, +0.360]** | 55W/71L |
| `only_pool_m` — **only** the retrieval-aligned | **4.035** | **+0.425 [+0.305, +0.550]** | 39W/87L |
| `only_sep_m` — **only** the zero-info-aligned | 3.961 | +0.351 [+0.200, +0.502] | 46W/80L |

Removing the cross-family-shared component recovers 62–79 % of the whitening bound. Keeping
**only** the shared component, at matched magnitude, is **worse than the predictor's actual
error** — all three, CI excluding zero.

> **VERDICT (ESTABLISHED).** The harmful coherent component is shared across independent
> families, and **two thirds of it is reproduced by references that know nothing about the
> sequence.** Every predictor — conditioned or not — emits a *typical peptide of that length*,
> and the harmful coherent error is the systematic difference between typical and this native.
> That is an **identifiability statement about the inputs**, not a defect of any architecture.

**I will not overstate it.** Conditioning is not worthless: 0.833 (same architecture) and 0.784
(retrieval) sit well above 0.598 (helix). The claim is that the *majority* of the harmful
component survives with no sequence information — not that sequence buys nothing.

**One structural observation.** The single most harmful shared component is the one aligned with
the **retrieval pool** (0.784, `only_pool_m` +0.425) — and the retrieval pool is *both* what
generates the candidate set *and* where the training fragments come from. The pipeline is
self-consistently biased toward the same "typical peptide", at generation and at prediction,
from one source.

---

# 5. WHICH COHERENT MODE — **one mode owns it, and it is a per-target separation profile**

`s19/a_struct.py` · `s19/results/a_struct.json` (COMPLETE, n=126) · all ORACLE.

**Part 1 — weighted variance of the residual explained (nested, per target):**
offset (1 param) 0.141 · stretch (1) 0.226 · separation profile (5) 0.361 · per-residue
additive (13) 0.430.
Standardised-residual correlation, pairs **sharing a residue** +0.087, sharing **none** −0.061,
within-target permutation null −0.017 — the coherence is a real, local, share-a-residue effect.

**Part 2 — surgery through the deployed fit. Magnitude-matched twins are the honest arms**,
because every removal shrinks the residual and an arm could otherwise win by having less error.
The Sprint-18 gap here is `real − shuffled` = **+1.000 Å**.

| mode removed | RMSD (matched) | vs `real` | share of the gap |
|---|---|---|---|
| per-target **offset** | 3.784 | **+0.174 [+0.092, +0.265]** | **−17.4 %** — removing it *hurts* |
| per-target **stretch** | 3.679 | +0.069 [+0.006, +0.135] | −6.9 % — removing it *hurts* |
| per-target **separation profile** | **3.085** | **−0.525 [−0.721, −0.343]** | **+52.5 %** |
| per-residue **additive** | 3.536 | −0.073 [−0.181, +0.036] | +7.3 % — **NOT MEASURED** |

> **The harmful coherent mode is a per-target separation profile: five numbers per target.**
> The model gets *this* target's distance-versus-separation curve wrong in a target-specific
> way. A per-target constant and a global stretch are *not* the problem — removing either makes
> the fit worse. The per-residue additive mode explains the most variance (0.430) and owns none
> of the harm, which is a clean warning that variance explained does not price damage.

---

# 6. **IS THE HARMFUL MODE ESTIMABLE NATIVE-FREE?  NO.** The wall, measured.

`s19/a_sepfix.py` · `s19/results/a_sepfix.json` (COMPLETE, n=126). Five numbers per target is
small enough that a native-free estimator might exist. It does not. For each separation bin the
ORACLE correction is `mean_b(d̂) − mean_b(d_true)`; each arm replaces `mean_b(d_true)` with a
native-free estimate. **Only `ORACLE_no_sep` reads the native.**

| arm | RMSD | median | resid RMS | vs `real` | W/L | folds |
|---|---|---|---|---|---|---|
| **incumbent (coordinate average)** | **3.048** | | | | | |
| `real` (deployed) | 3.610 | 3.370 | 3.205 | — | | |
| `ORACLE_no_sep` (**the ceiling**) | **3.085** | 3.137 | 3.159 | −0.525 [−0.721, −0.343] | 84W/42L | 5/5 |
| `poolfull` (replace d̂ by the pool's distances) | **3.318** | 3.156 | 2.528 | −0.292 [−0.450, −0.131] | 80W/46L | 5/5 |
| `poolmed_sep` | 3.455 | 3.277 | 3.186 | −0.155 [−0.266, −0.048] | 81W/45L | 5/5 |
| `pool_sep` | 3.512 | 3.388 | 3.159 | −0.098 [−0.195, **−0.004**] | 71W/55L | 4/5 |
| `pool_sep_half` | 3.554 | 3.309 | 3.138 | −0.056 [−0.120, +0.008] | 77W/49L | 4/5 |
| `rand_sep` (**matched-random**) | 3.865 | 3.782 | 3.514 | +0.256 [+0.145, +0.372] | 49W/77L | 5/5 |
| `helix_sep` (**zero-information**) | 3.878 | 3.760 | 3.746 | +0.268 [+0.077, +0.466] | 59W/67L | 4/5 |
| `global_sep` (**zero-information**) | 4.349 | 4.253 | 3.709 | +0.739 [+0.537, +0.942] | 32W/94L | 5/5 |

**It clears the controls I was required to carry** — `pool_sep` beats matched-random by
−0.354 [−0.493, −0.214] and the constant-helix zero-information control by −0.366
[−0.536, −0.199], both 5/5 folds. **And it fails the trap I pre-stated in the module header
before running it:** `pool_sep − poolfull` = **+0.193 [+0.070, +0.324], 44W/82L**.

> **The five-parameter profile match is dominated by simply replacing d̂ with the pool's
> distances wholesale, and the whole ladder — 3.610 → 3.554 → 3.512 → 3.455 → 3.318 — is
> monotone interpolation toward `poolfull`. It is "move toward the retrieval pool", not "estimate
> the coherent mode".** Under brief §12 ("reproduced by a simpler classical baseline") this is
> **DOWNGRADED**, not a result. And every arm on the board, including `poolfull`, is worse than
> the 3.048 Å incumbent.

The native-free route recovers **0.098–0.155 Å of the 0.525 Å ORACLE ceiling — 19–30 %** — and
only by regressing toward what the pipeline already does.

> **CONCLUSION: an identifiability wall.** The harmful mode is nameable, localisable, and worth
> 0.525 Å with the native in hand. The best native-free estimator of it available in this
> project is the retrieval pool, which shares the bias being estimated (§4). **There is still no
> native-free estimator of coherence.**

---

# 7. THE INTERVENTION LEDGER

| # | branch | mechanism | primary metric | result | status |
|---|---|---|---|---|---|
| A1 | error concentrated in the low-utility region | topology | ρ(\|r\|, sep) +0.466; ρ(\|r\|, 1/sd²) −0.520; 9.4 % of SSE in the high-utility cell | confirmed, and explicitly **not** the mechanism | **ESTABLISHED** |
| A3 | the error is coherent | sign correlation | **`signflip_exact` −1.202 [−1.408, −1.004], 5/5 folds, 128 % of the gap at exactly matched magnitude** (§1.0); `coherent` +0.139, `coh_scaled` +0.233 | mechanism | **ESTABLISHED** |
| A4 | the field is not jointly realisable | EDM defect | 0.286 vs native 0.0018 — but a matched *incoherent* field is worse (0.340) and lands 1.24 Å better | measured, **anti-correlated with harm** | **REFUTED as a mechanism** |
| A2 | harmful error in a detectable subset | 25 % ORACLE repair vs size-matched random | **no native-free criterion beats random at matched magnitude; all seven are worse.** ORACLE-perfect detection itself buys only −0.205 [−0.364, −0.041] | §8 | **REFUTED** |
| A5a | training loss weighted **toward long range** (anti-utility) | `sw_lin` vs `sw_none` | **+0.181 [+0.081, +0.285]**, 49W/77L, 4/5 folds; ORACLE MAE +0.010 [−0.039, +0.062] | the compass reproduced at the predictor | **ESTABLISHED (negative)** |
| A5b | training loss weighted **toward short range** (utility) | retrain 5 folds | **NOT RUN** — 6,429/6,788 training sequences absent from the hot ESM cache, needs the 1.5 GB bank, free RAM 1.81–1.90 GB. Declined under brief §15 and `machine-fits-two-heavy-jobs` | pre-registration unedited | **OPEN — NOT MEASURED** |
| A6 | joint consistency enforcement (PairNet triangle update) | architecture | tri 4.09 %→0.68 %, defect 0.286→0.148, κ→0.914 (best on the board); RMSD −0.080 [−0.242, +0.082] | works geometrically, does not convert | **CLOSED** |
| A7 | decorrelating families (ensembles) | `d20_ens`, `combo` | −0.018 [−0.141, +0.108]; −0.094 [−0.191, **+0.001**] | CI touches zero | **NOT MEASURED** |
| A8 | conditional weighting law w(sd, sep, SS, pair type) | — | not run; §3–§6 make its premise (accuracy is the lever) false | superseded | **CLOSED by §4** |
| A9 | native-free estimate of the harmful mode | per-target separation profile from the retrieval pool | −0.098 to −0.155 vs real, clears random and zero-info controls, **loses to `poolfull` +0.193 [+0.070, +0.324]**; all arms worse than the 3.048 incumbent | interpolation toward the pool | **DOWNGRADED** |

## What I would tell the next sprint

0. **Quote −1.202 [−1.408, −1.004], not −1.242** — see the self-audit in §1.0, which found and
   removed a 2.0 Å-clip magnitude confound in my own headline arm.
1. **The mechanism is the residual's cross-pair sign correlation.** Anything that changes error
   *magnitude* without changing its *correlation structure* is measured to be worth nothing:
   eight predictors, 11.7 % of MAE, zero CIs excluding zero.
2. **Do not spend on metric realisability.** Making the field more realisable makes it more
   coherently wrong — confirmed by construction (`coherent` +0.139) and by architecture
   (PairNet, κ 0.914, −0.080 [−0.242, +0.082]).
3. **Do not spend on predictor interventions.** The harmful component is shared with a different
   architecture (0.81 of the same-seed ceiling), with a non-network retrieval baseline (0.87),
   and two thirds of it with references that have no sequence information at all.
4. **The one live number is ORACLE.** 2.353–2.368 Å is a ceiling, not a method, and must never
   be quoted as an achieved number.
5. **If anyone reopens this**, the object to attack is the retrieval pool's own bias, because it
   is the largest identified shared source and it enters the pipeline twice — as the candidate
   generator and as the predictor's training fragments.

---

# 8. Q1b — **A2 REFUTED.** The harmful error is not in a detectable subset

`s19/a_subset.py` · `s19/results/a_subset.json` (COMPLETE, n=126). Repair 25 % of the pairs
(set `d̂ := d_true` on them) and ask whether any native-free criterion for choosing them beats a
size-matched **random** choice. The `matched` column rescales the surviving residual so **every
arm carries the same per-target residual RMS** — without it a criterion wins purely by picking
big residuals, which is what the `raw` column shows and why it must not be read alone.

| criterion | RMSD raw | resid RMS | vs random (raw) | **RMSD matched** | **vs random (matched)** |
|---|---|---|---|---|---|
| `nf_model` (LFO native-free error model) | 2.937 | 2.138 | −0.163 [−0.304, −0.028] | 3.348 | **+0.232 [+0.103, +0.356]** 43W/83L |
| `hi_sd` (least confident) | 3.108 | 2.308 | +0.007 [−0.095, +0.112] | 3.436 | **+0.320 [+0.220, +0.423]** |
| `lo_sd` (most confident) | 3.459 | 3.117 | +0.358 [+0.233, +0.489] | 3.277 | **+0.161 [+0.034, +0.285]** |
| `long_sep` | 2.836 | 2.145 | −0.264 [−0.411, −0.119] | 3.210 | +0.093 [−0.053, +0.232] |
| `short_sep` | 3.535 | 3.150 | +0.434 [+0.299, +0.568] | 3.286 | **+0.170 [+0.053, +0.291]** |
| `terminal` | 2.956 | 2.450 | −0.145 [−0.269, −0.017] | 3.182 | +0.066 [−0.042, +0.181] |
| `multimode` (highest entropy) | 3.208 | 2.686 | +0.107 [+0.000, +0.215] | 3.273 | **+0.157 [+0.054, +0.259]** |
| **`random` (the control)** | 3.101 | 2.791 | — | **3.116** | — |
| `ORACLE_absr` (**perfect detection**) | 2.180 | 1.430 | −0.921 [−1.104, −0.741] | 2.912 | **−0.205 [−0.364, −0.041]** |

> **Not one native-free criterion beats the size-matched random repair at matched magnitude —
> all seven are worse, five with CIs excluding zero.** And **ORACLE-perfect detection of the
> largest 25 % of errors is worth only −0.205 Å** once magnitude is held fixed. The `raw` column
> is entirely a magnitude confound: the criteria that "win" there (`long_sep`, `nf_model`,
> `terminal`) are the ones that strip the most squared error (resid RMS 2.14–2.45 against
> random's 2.79).

**A2 is REFUTED, and the §2.3 native-free error model (Spearman +0.537 against the ORACLE error)
does not convert.** Detecting the error is partly possible; repairing what you detect is worth
nothing.

**A note on the apparent tension with Sprint-18 G2e, which is not a contradiction.** G2e
(`fix_confident` −0.308, `fix_short` −0.294 vs uniform) *shrinks* residuals within a class while
keeping every pair; this experiment *zeroes* a class and rescales the complement. Mine
concentrates the surviving residual into a structured subset and then amplifies it — which
**increases** its coherence — so a random repair, which leaves the field's structure intact,
wins. The two results agree with the §1 mechanism and disagree about nothing; they are different
operations and must not be quoted against each other.
