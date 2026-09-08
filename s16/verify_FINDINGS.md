# SPRINT 16 / VERIFY — HOSTILE RE-READING OF SPRINT 16'S OWN CLAIMS

**Modules** `s16/verify_steer.py`, `s16/verify_csteer.py`, `s16/verify_cosmag.py`,
`s16/verify_stats.py`, `s16/verify_ties.py`
**Artefacts** `s16/results/verify_steer.json`, `verify_csteer.json`, `verify_cosmag.json`,
`verify_stats.json`, `verify_ties.json`
**Mandate** assume every Sprint 16 claim is wrong until reproduced from the artefacts; the
coordinator's first. **No new science** — reproduce, control, break.

> **Every number in Sprint 16 that I checked reproduces from its artefact exactly.** As in
> Sprint 15, that is not the question. What fails is, three times over, the same failure RETRACT
> documented on the fusion law: **a quantity is measured correctly and then read as if it were a
> different quantity.** A transcription audit passes all three.

---

## 0. VERDICTS, worst first

| # | claim | verdict |
|---|---|---|
| 1 | `CSTEER_FINDINGS` §2 — "from the best start every native-free direction has a cosine with the true residual indistinguishable from zero"; "the direction does not exist" | **REFUTED.** The module quotes the **signed** mean of a quantity its own law consumes as **c²**. Read as \|c\| against the matched random control the sprint's best channel is **0.305 vs 0.128, +0.177 [+0.142, +0.212], W/L 100/26**, and a **zero-information constant ideal β-strand beats every channel the sprint built** at \|c\| **0.398 (+0.262 [+0.216, +0.310])**. The ladder also has **no negative rung**, so on 45–67% of targets no arm could have used what is there. |
| 2 | `STEER_FINDINGS` §3 — the mechanism: "the fit's torsion errors are *compensating* … moving part-way breaks the cancellation" | **REFUTED as a causal claim** (the numbers reproduce exactly). A native-free matched control — geodesic torsion interpolation between two arbitrary retrieval windows — reproduces the curve and, at matched angular distance, **bulges more** (−0.024 vs +0.014 at s = 0.5). The effect is a function of \|\|Δθ\|\| alone (Spearman **−0.571**). The **geodesic objection VERIFY was asked to test is itself REFUTED**: `A.wrap` gives the exact per-angle geodesic. |
| 3 | `CSTEER_FINDINGS` §1 — "coordinate space is exactly linear, to three decimals at every rung, and the identity survives Kabsch superposition"; "the sprint's cleanest exhibit" | **TAUTOLOGY.** Demonstrated on 400 random point-cloud pairs with no protein in them: max deviation **4.6e−07 Å**. Same species as the fusion law RETRACT retracted **on the same day, in the same sprint**. |
| 4 | `INTEGRATE_FINDINGS` §2 / L17 — "AMBER as a ranker is −0.054 [−0.088, −0.020] at m = 75 … the only positive accuracy effect in Sprint 16" | **DOES NOT SURVIVE.** Target as the unit: **−0.054 [−0.148, +0.040]**, CI **2.79× wider**, contains zero, median −0.045, **5 targets better / 4 worse**. Removing any one of **8 of 9** targets puts the CI over zero; **seed 2 alone gives −0.0003**; **no per-generator CI excludes zero**. The coordinator self-declared the row-bootstrap defect (L17) — this is the recomputation, and the point estimate does not carry the claim. **The tie trap does NOT apply** (censused: zero ties). |
| 5 | `INTEGRATE_FINDINGS` §4 / L18 — "the generators have almost identical mean member error … they differ in diversity, not in conformer quality; the terminal operator's entire gain is the diversity term" | **HALF WRONG.** Panel A is an identity at machine epsilon (correctly labelled). But the table prints the **free-superposition** member error (spread 0.056 Å) while the identity uses the **common-frame** error (spread **0.122 Å**). In the identity's own terms the readout spread is **~50% member error** for three of the four CVaR arms (cvar0.25 **52%**, cvar1.0 51%, cvar0.5 42%). |
| — | `CSTEER_FINDINGS` §3 (the endpoint-substitution trap), §4 (projection costs 0.165 Å); `INTEGRATE_FINDINGS` §5 (coordinate beats torsion readout); §3 (the VQE loses to uniform) | **SURVIVE.** See §6. |

---

## 1. CLAIM 1 — the pivot measured the wrong statistic, and a zero-information reference beats it

### 1.1 What `csteer.py` derives, and what `CSTEER_FINDINGS.md` reports

The module's own docstring derives the law, correctly:

    ||r − t v||² = ||r||² − 2t(v·r) + t²   minimised at   t* = v·r = c||r||
    giving                                                ||r_new|| = ||r|| √(1 − c²)

**`t*` carries the sign of `c`, and the payoff depends on `c²`.** The law is **two-sided**: a
direction at c = −0.4 is worth exactly what one at c = +0.4 is worth — you step backwards along
it. The code implements this correctly (`t_star = cc*nr/nd`, `law_pred = base*sqrt(1-cc*cc)`).

`CSTEER_FINDINGS.md` §2 then tabulates and quotes the **arithmetic mean of the signed `c` over
targets**, and concludes from `0.012 [−0.052, +0.075]` that "the direction does not exist".
Averaging the signed cosine over targets destroys exactly the quantity the law consumes.

### 1.2 The artefact, re-read as the law reads it (`s16/verify_cosmag.py`, C1)

From `avg` (3.048 Å, the best start), n = 126, fold-clustered paired CIs, matched random
control taken from the module's own `rand0..2` arms:

| direction | signed c (as published) | **\|c\|** | random \|c\| | **\|c\| − random, paired [95% CI]** | W/L | ORACLE two-sided ceiling |
|---|---|---|---|---|---|---|
| toward `fit` | +0.012 | **0.305** | 0.128 | **+0.177 [+0.142, +0.212]** | **100/26** | 2.843 Å |
| toward `ens` | −0.004 | **0.301** | 0.128 | **+0.173 [+0.135, +0.210]** | 98/28 | 2.836 Å |
| toward `med` | −0.010 | 0.251 | 0.128 | +0.123 [+0.087, +0.160] | 83/43 | 2.894 Å |
| toward `proj` | −0.068 | 0.187 | 0.128 | +0.059 [+0.033, +0.087] | 70/56 | 2.970 Å |

The random control returns **0.128**, which is `√(2/(π·3n))` at n = 13 — the correct isotropic
null, so the excess is real and not a magnitude artefact of taking an absolute value.

### 1.3 The zero-information reference AUDIT made mandatory — and it fires (`verify_cosmag.py`, C2)

Same start, same residual, same random control. These directions know **nothing** about the
target except its chain length:

| direction | kind | signed c | **\|c\|** | \|c\| − random [95% CI] | ORACLE two-sided ceiling |
|---|---|---|---|---|---|
| **constant ideal β-strand** | **ZERO-INFORMATION** | −0.012 | **0.398** | **+0.262 [+0.216, +0.310]** | **2.491 Å** |
| constant ideal PPII | ZERO-INFORMATION | −0.008 | 0.382 | +0.246 [+0.200, +0.293] | 2.561 Å |
| constant ideal α-helix | ZERO-INFORMATION | −0.027 | 0.335 | +0.199 [+0.155, +0.243] | 2.733 Å |
| an arbitrary retrieval window | ZERO-INFO | −0.015 | 0.276 | +0.140 [+0.103, +0.179] | 2.859 Å |
| pool medoid | native-free channel | −0.010 | 0.251 | +0.115 [+0.077, +0.152] | 2.894 Å |
| **pure isotropic breathing** (the start structure scaled about its own centroid) | **ZERO-INFORMATION** | −0.067 | **0.296** | **+0.167 [+0.123, +0.214]** | 2.708 Å |
| random | control | −0.016 | 0.136 | +0.000 | 3.009 Å |

**A constant ideal β-strand carries more usable direction from the best start than any channel
Sprint 16 built.** A pure breathing mode — scale the starting structure up or down about its own
centroid, an operation with literally no input — matches the best channel (0.296 vs 0.305).

This is the mechanism: the residual `avg → native` is dominated by a **compactness scalar**.
Every extended reference direction (strand, PPII, breathing) aligns with it in **magnitude**, and
**58% of targets need the contracting sign and 42% the expanding one**. That is the standing
project record — *averaging contracts the backbone 25.8%*, and *the only leverage supplies the
per-target SIGN at inference* — reappearing exactly.

### 1.4 The ladder cannot take the step the law prescribes (`verify_cosmag.py`, C3)

`csteer.STEPS = (0.0, 0.1, …, 1.0)`. There is **no negative rung**. Fraction of targets whose
optimal step from `avg` is negative:

| direction | `fit` | `proj` | `med` | `ens` |
|---|---|---|---|---|
| fraction needing a NEGATIVE step | **0.45** | **0.67** | **0.52** | **0.48** |

On roughly half the instrument, **no arm on that ladder could have helped regardless of what the
direction carried.** A leave-fold-out selector choosing a single non-negative step for a
direction whose sign flips per target must return "do nothing", and it did. The observed null is
therefore partly a property of the ladder, not only of the channels.

### 1.5 What survives, restated correctly

**The operational conclusion survives.** Nothing here is deployable: the sign of `c` and the
magnitude `||r||` are both native-derived, so every "ORACLE two-sided ceiling" above is a ceiling
and must never be quoted as an effect. Coordinate-space steering with a leave-fold-out step and
these channels buys nothing, and 2.5 Å is not reached.

**The stated reason is wrong, and the correction matters.** §2 should read:

> From the best start, every native-free direction has a cosine whose **sign is unpredictable**
> across targets (signed mean 0.012 [−0.052, +0.075]) while its **magnitude is 2–3× the isotropic
> null** (\|c\| 0.305 vs 0.128, W/L 100/26) — and a zero-information constant β-strand has more of
> it (0.398) than any channel we built. The missing ingredient is not the direction. It is the
> **per-target sign**, one bit, plus the step magnitude.

The difference is not cosmetic. "The direction does not exist" closes the line. "The direction
exists at \|c\| ≈ 0.4 from a zero-information reference, and one bit per target is missing"
leaves an ORACLE ceiling of **2.491 Å** — below the sprint's own primary goal — sitting behind a
scalar. Whether that bit is obtainable native-free is **not tested here** (it would be new
science) and is the open question this refutation leaves.

---

## 2. CLAIM 2 — the flagship's numbers are right, its mechanism is not, and the geodesic objection is dead

### 2.1 Reproduction: exact (`s16/verify_steer.py`, Q1)

Recomputed from `s16/results/steer.json`, n = 126:

| step s | 0.0 | 0.1 | 0.2 | 0.35 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|---|---|---|
| mean RMSD (Å) | 3.647 | 3.599 | 3.641 | 3.663 | 3.466 | 2.596 | 0.043 |
| frac of gain, ratio-of-means | 0 | 1.33% | 0.16% | −0.43% | **5.02%** | 29.16% | 100% |
| frac of gain, **per target** | 0 | 1.60% | 0.83% | 1.10% | **7.66%** | 32.31% | 100% |
| improved / worsened | — | 82/44 | 75/51 | 73/53 | **77/49** | 104/22 | 126/0 |

Published: 7.7% [−0.0, +14.9], median 16.1%, 77/49. Recomputed: **7.66%**, fold-clustered
**[−0.04%, +14.85%]**, median **16.08%**, **77/49**. **The published table matches the artefact,
the per-target/ratio-of-means switch is stated in the writeup and is correct, and the per-target
aggregation is the right one.** No error here.

### 2.2 The geodesic objection is REFUTED (`verify_steer.py`, Q2)

`θ_fit + s·A.wrap(θ_native − θ_fit)` with `wrap(a) = (a + π) mod 2π − π` **is** the per-angle
geodesic on the torus. Checked against an independently written great-circle construction
(signed angle from `atan2` of the cross and dot products of the two unit vectors), on 200 random
pairs deliberately placed across the branch cut × 7 steps:

    max |wrap(coordinator − independent geodesic)|  =  1.78e−15 rad
    max CA-RMSD between the two builds              =  3.24e−07 Å
    max |wrap(coordinator − NAIVE unwrapped lerp)|  =  3.14e+00 rad   (the interpolant NOT used)

The naive interpolant differs by up to π and was not used. **The non-monotonicity is not a
wrapping artefact.**

### 2.3 The mechanism claim, and the control the flagship does not have (`verify_steer.py`, Q3)

`STEER_FINDINGS` §3 asserts a cause: *"The fit's torsion errors are compensating — wrong in a
correlated way that partly cancels along the chain … Moving part-way breaks the cancellation
without establishing the correct structure."* That is a causal claim about the **error**, made
with no matched control. The competing explanation is that geodesic torsion interpolation
between **any** two distant conformations bulges — a property of the **parameterisation**.

The missing control, entirely native-free: geodesically interpolate between two arbitrary
top-75 retrieval windows (1,510 pairs over the 126 targets), and from a real window toward a
random torsion vector of matched norm (1,512 pairs). Statistic identical to the flagship's,
`frac(s) = 1 − d(x(s), T)/d(x(0), T)`.

**Matched on angular distance** — restricted to the flagship's own \|\|Δθ\|\| interquartile range
[6.30, 8.96] rad:

| arm | n | 0.1 | 0.2 | 0.35 | **0.5** | 0.7 |
|---|---|---|---|---|---|---|
| flagship, fit → native (ORACLE) | 62 | −0.001 | −0.024 | −0.038 | **+0.014** | +0.293 |
| **CONTROL, pool → pool (native-free)** | 493 | −0.007 | −0.035 | −0.069 | **−0.024** | +0.243 |
| CONTROL, pool → random target | 1446 | +0.013 | +0.028 | +0.055 | +0.127 | +0.373 |

**The native-free control is non-monotonic over the first third and sub-linear at half-way, more
so than the flagship arm.** Nothing about the fit's error is required to produce the curve.

The single variable that produces it is the angular distance:

| \|\|Δθ\|\| bin (rad) | fit → native | pool → pool | pool → random |
|---|---|---|---|
| [0.15, 1.30) | +0.554 (n=2) | +0.498 (n=271) | — |
| [1.30, 3.72) | +0.495 (n=8) | +0.472 (n=264) | — |
| [3.72, 5.71) | +0.408 (n=15) | +0.268 (n=258) | +0.227 (n=13) |
| [5.71, 6.94) | +0.242 (n=21) | +0.011 (n=252) | +0.186 (n=284) |
| [6.94, 8.24) | +0.020 (n=30) | −0.032 (n=242) | +0.118 (n=1215) |
| [8.24, 11.02) | **−0.145** (n=50) | −0.052 (n=222) | — |

Spearman(\|\|e\|\|, frac at s = 0.5) = **−0.571**. Small-\|\|e\|\| half: **+0.277**. Large-\|\|e\|\|
half: **−0.124**. **Below ≈ 4 radians torsion interpolation is linear** — for the flagship arm
(+0.50) and for the control (+0.47) alike.

### 2.4 Why the flagship's error sits in the bad regime — and what that makes the finding

The fit's exact native torsion error has, at the median target (n = 13, 26 angles):

    ||e|| = 7.79 rad;  per-angle RMS error = 90.2 degrees
    a uniformly random assignment of every torsion gives ||.|| = 9.25 rad
    ||e|| / random = 0.868 (median);  59% of targets exceed 0.80 of it

**The "exact native torsion error" is 87% of the way to a random angle assignment.** "Half-way
along it" is therefore a walk to a conformation essentially unrelated to either endpoint, and
that it is not worth half the gain is arithmetic about distance, not a discovery about error
structure.

### 2.5 What §3 should say

- **KEEP.** The numbers, the aggregation, the pre-registered rules and which fired, and the
  operational conclusion: at *this instrument's error magnitude* a leave-fold-out partial step in
  torsion space pays essentially nothing, so the flagship is refuted.
- **WITHDRAW.** "The fit's torsion errors are *compensating* … moving part-way breaks the
  cancellation" — refuted by the pool→pool control at matched distance. Also withdraw the link
  to Sprint 15's realizability result, which the control shows is not needed.
- **WITHDRAW as stated.** "Torsion space is not a space in which you may move part-way." It is,
  below ≈ 4 rad, where the curve is linear on this instrument. The correct statement is a
  **magnitude regime**: partial motion in torsion space is linear for small displacements and
  breaks down beyond ≈ 6 rad, and the prediction error is at 7.8 rad.
- **REPHRASE.** "No improvement to the direction estimator would rescue it" is too strong. What
  is shown is that no scheme stepping *a large fraction of a near-maximal error* can work. A
  small step in a good direction stays in the linear regime; it simply pays proportionally
  little, which is a different and weaker statement.

---

## 3. CLAIM 3 — the coordinate-linearity exhibit is a theorem (`s16/verify_csteer.py`, P1)

**Theorem.** Let `N` be Kabsch-superposed onto `S`, so the identity rotation is optimal for the
pair; equivalently the centroids coincide and the cross-covariance `H = Nᶜᵀ Sᶜ` is symmetric PSD.
For `X_s = S + s(N − S)`,

    Xᶜᵀ Nᶜ = (1 − s) Sᶜᵀ Nᶜ + s Nᶜᵀ Nᶜ = (1 − s) Hᵀ + s Nᶜᵀ Nᶜ

is a sum of two symmetric PSD matrices, hence symmetric PSD, so the identity stays Kabsch-optimal
along the whole segment and `rmsd(X_s, N) = (1 − s)·rmsd(S, N)` **exactly, for every s**.

Demonstrated on **400 random point-cloud pairs (n = 5…40), 9 rungs, no protein involved**:

    max | rmsd(S + s·r, N) − (1 − s)·rmsd(S, N) |  =  4.559e−07 Å

**`3.647·(1 − s)` to three decimals at every rung is a check that floating point works.** It is
described in `CSTEER_FINDINGS` §1 as "a real and, as far as the programme's own record goes,
previously unmeasured fact" and in L13 as "the sprint's cleanest exhibit". It is neither: it is
the same PSD-preservation argument RETRACT used, the same day, to retract the fusion law's
"parameter-free validation to 0.162 Å" — and `csteer.py`'s docstring **already writes the algebra
out**, so the module knew.

**What content remains in the torsion/coordinate contrast:** only the torsion half, and §2.3
above shows the torsion half is a **displacement-magnitude regime** rather than a property of the
space. The honest one-line contrast is: *linear interpolation in a linear space is linear by
construction; geodesic interpolation on a torus pushed through a nonlinear embedding is linear
for small displacements and not for large ones, and our error is large.* That is worth stating.
It is not an exhibit.

---

## 4. CLAIM 4 — the integration result does not survive the target as the unit (`s16/verify_stats.py`)

`integrate._boot` and `diversity._boot` both resample **rows**, and a row is a
(target × seed × generator) cell: **162 rows on 9 targets**. BRIEF §3.4 makes the target the unit.
The coordinator declared this defect himself (L17, and at the head of `INTEGRATE_FINDINGS.md`)
before VERIFY reported; this is the recomputation he asked for. **Point estimates are unchanged
by construction** (the per-target means are balanced). The intervals are not.

### 4.1 The recomputation, m = 75

| contrast | ROWS (published) | **TARGET (n = 9)** | median | W/L (targets) | CI widening |
|---|---|---|---|---|---|
| AMBER = amb − ctrl | −0.0536 [−0.0869, −0.0198] | **−0.0536 [−0.1476, +0.0397]** | −0.045 | 4 worse / 5 better | **2.79×** |
| AMBER = amb − rnd | −0.0539 [−0.0878, −0.0188] | **−0.0539 [−0.1489, +0.0409]** | −0.047 | 4/5 | 2.75× |
| LEGACY = leg − ctrl | +0.0052 [−0.0339, +0.0440] | +0.0052 [−0.1095, +0.1242] | −0.037 | 4/5 | 3.00× |
| INTER | −0.0044 [−0.0416, +0.0323] | −0.0044 [−0.1040, +0.0966] | +0.016 | 5/4 | 2.72× |
| **tors − ctrl** | +0.4463 [+0.3557, +0.5359] | **+0.4463 [+0.1934, +0.7057]** | +0.227 | **8/1** | 2.84× |

At m = 20 the AMBER − rnd contrast likewise goes **−0.0513 [−0.0904, −0.0125] → −0.0513
[−0.1442, +0.0504]**, and at m = 5 the "Legacy is significantly worse than random" claim goes
**+0.0684 [+0.0128, +0.1263] → +0.0684 [−0.0389, +0.1757]**.

### 4.2 Fragility, beyond the interval

Per-target AMBER effect at m = 75 (negative = AMBER helps):

    1CS9 +0.190   2MK7 −0.240   2P5H −0.045   6EY3 +0.040   6F3V +0.044
    6S0N −0.270   7N2I −0.185   8IS3 −0.067   9UV5 +0.050

- **Leave one target out:** dropping any of **8 of the 9** leaves a CI containing zero. Only
  dropping 1CS9 — the largest *adverse* target — strengthens it, which is the wrong direction for
  a robustness argument.
- **Per generator** (target as the unit): cvar0.1 −0.053 [−0.161, +0.048]; cvar0.25 −0.080
  [−0.190, +0.026]; cvar0.5 −0.076 [−0.196, +0.053]; cvar1.0 −0.042 [−0.190, +0.101]; bestofN
  −0.044 [−0.091, +0.004]; uniform −0.027 [−0.140, +0.099]. **None excludes zero.**
- **By seed:** −0.093, −0.068, **−0.0003**. It does not replicate in the third seed.

### 4.3 The tie trap does NOT apply (`s16/verify_ties.py`)

`_readouts` calls `np.argsort(legacy, kind="stable")` on an `idx` that `_pick` has sorted by the
objective, so a tie would break toward better objective rank, which correlates with RMSD, and
would hand `leg`/`amb` skill the random control cannot get. Censused by regenerating the
`uniform` ensembles bit-identically (they are `stable_rng`-seeded, no VQE needed), m = 75:

| pdb | prior ties | Legacy ties | AMBER ties | ρ(idx rank, RMSD) |
|---|---|---|---|---|
| 1CS9 | 30 | **0** | **0** | +0.060 |
| 6F3V | 5 | **0** | **0** | +0.241 |
| 9UV5 | 22 | **0** | **0** | −0.134 |

**Zero ties in either energy.** The trap is real in the code path and does not bite in this data.
`en.prior` does tie (5–30 of 75), but `_pick` breaks those by configuration index — enumeration
order, not oracle order — and ρ(idx rank, RMSD) is +0.06 to +0.24, so the leak channel is weak
even where it exists. **Clean survival on this point.**

### 4.4 Verdict

"AMBER contributes … it is the only positive accuracy effect anywhere in Sprint 16" **must be
withdrawn as a claim with an interval.** What is supportable: *the AMBER point estimate is
−0.054 Å at m = 75 with 5 of 9 targets improving and one seed at zero; at n = 9 targets the
interval contains zero and the effect is carried by no stable subset.* Legacy's null stands and
is now merely wider. The **AMBER − Legacy** contrast (−0.059 published) is the one worth re-running
at n = 126 if anything here is to be promoted; on 9 targets it cannot be decided.

---

## 5. CLAIM 5 — the diversity account is a theorem plus a mis-attributed measurement

### 5.1 Panel A is an identity, and it is the parallel-axis theorem (`verify_stats.py`)

`I.coordinate_average(W, P)` returns `superpose_batch(W, W[b]).mean(0)` — the plain mean of the
members **in the medoid frame**. `diversity.py` then measures `err2` and `div2` in that same
frame. For any point `y`, `mean_i ||x_i − y||² = mean_i ||x_i − x̄||² + ||x̄ − y||²`. With
`y = natm` that **is** the Krogh–Vedelsby line, exactly, with no protein content.

    residual (RMSD units)   mean −3.13e−16   max|·| 7.55e−15
    residual (squared form) mean −2.84e−15   max|·| 7.11e−14
    max relative to readout²                 5.50e−15      (machine epsilon ≈ 2.2e−16)

The coordinator labels this "a theorem check, not a discovery" in both L18 and
`INTEGRATE_FINDINGS` §4. **That labelling is correct and is credited.**

### 5.2 The content claim uses a different quantity from the identity — and it changes the answer

The identity's error term is `err2_common_frame` (members superposed to the medoid, native to the
average). The findings table prints `mean_member_quad`, the **free-superposition** per-member
RMSD, which is a different and systematically smaller quantity. At m = 75:

| generator | readout | free-sup RMSq **(published)** | **COMMON-FRAME err (the identity's term)** | diversity |
|---|---|---|---|---|
| cvar0.1 | 3.426 | 3.950 | **4.216** | 2.393 |
| cvar0.25 | 3.515 | 4.006 | **4.301** | 2.418 |
| cvar0.5 | 3.448 | 3.971 | **4.252** | 2.425 |
| cvar1.0 | 3.423 | 3.983 | **4.261** | 2.475 |
| bestofN | 3.311 | 3.966 | **4.180** | 2.501 |
| uniform | 3.214 | 3.966 | **4.179** | 2.626 |
| **spread** | **0.301** | **0.056** | **0.122** | 0.233 |

The published "the mean member error spans 0.056 Å" is true of a quantity the identity does not
contain. The identity's own error term spans **0.122 Å, 2.2× larger.**

Decomposing the readout gap exactly, per target, each generator against `uniform`
(`d(readout²) = d(err²) − d(div²)`, an identity, so the split is exact):

| generator | d(readout²) | from ERROR | from DIVERSITY | **error share** |
|---|---|---|---|---|
| cvar0.25 | +2.155 | +1.131 | +1.024 | **0.52** |
| cvar1.0 | +1.558 | +0.801 | +0.757 | **0.51** |
| cvar0.5 | +1.749 | +0.743 | +1.006 | **0.42** |
| cvar0.1 | +1.615 | +0.461 | +1.154 | 0.29 |
| bestofN | +0.590 | −0.045 | +0.635 | 0.07 |

**For three of the four CVaR arms the readout penalty is about half member error and half
diversity.** Only the untrained-circuit control is "almost entirely diversity".

And the premise itself was never demonstrated: with the target as the unit (n = 9), the
common-frame member error differences against `uniform` are +0.037 [−0.135, +0.199], +0.122
[−0.037, +0.286], +0.073, +0.082, +0.001 — **every CI contains zero in both directions**, so
"almost identical mean member error" is *not distinguishable from zero at n = 9*, which is not the
same as equal.

### 5.3 Verdict

- **KEEP.** The identity, correctly labelled. The diversity column and its ordering. The
  within-target rank correlations (−0.145 to −0.270) and their honest framing as a weak heuristic.
  The divselect null and its "both terms rose and cancelled" reading (which my analysis supports —
  the terms genuinely do carry comparable weight).
- **WITHDRAW.** "They differ in diversity, not in conformer quality"; "the terminal operator's
  entire gain is the diversity term"; "saying which of its two terms the generators differ in is a
  complete causal account, and the answer is diversity". The account is complete only if the
  identity's own error term is used, and then the answer is **roughly half and half**.
- **WITHDRAW.** L18's "This corrects a recorded project law … the set means are equal to 0.056 Å".
  Measured in the frame the identity requires, the set means are **not** equal, and the standing
  law "the terminal operator consumes the set MEAN" retains about half the effect. The correct
  refinement is *the operator consumes the set mean and the set diversity, in comparable measure*,
  not *diversity replaces the mean*.

---

## 6. WHAT SURVIVED EVERYTHING

A clean survival is a result. These were attacked and held:

1. **Every published number reproduces from its artefact.** `steer.json`: 7.66% (published 7.7),
   median 16.08% (16.1), 77/49, ratio-of-means 5.02% (5.0), curvature ratio 2.709, cosine 0.318,
   top-4 share 0.893. `csteer.json`: every cosine in §2 to three decimals, and fold-clustering the
   CIs (which `csteer._boot` does not do) changes nothing material. `integrate.json` and
   `diversity.json`: every point estimate.
2. **The flagship's interpolant is the correct geodesic** (§2.2). The objection VERIFY was
   commissioned to test does not stand.
3. **The flagship's headline aggregation is the right one.** The writeup switches to the
   per-target mean, says so, prints the ratio-of-means beside it, and the published table is the
   per-target one. This is the opposite of the failure mode being hunted.
4. **No tie leak in the integration** (§4.3), on the exact code path the programme was burned on.
5. **The coordinate readout beats the torsion readout**, and it survives the target as the unit:
   **+0.446 [+0.193, +0.706] at m = 75, 8 targets of 9**, +0.297 [+0.088, +0.509] at m = 20,
   +0.214 [+0.062, +0.363] at m = 5. Monotone in ensemble size. **This is the one accuracy-relevant
   claim in `INTEGRATE_FINDINGS` that survives target-level analysis.**
6. **The VQE loses to uniform random at matched budget**, target as the unit, m = 75: cvar0.25
   **+0.301 [+0.115, +0.485]** (7/2), cvar0.5 +0.234 [+0.011, +0.458], cvar0.1 +0.212
   [+0.003, +0.404]; cvar1.0 +0.209 [−0.002, +0.410] now touches zero. Against the untrained
   circuit only cvar0.25 excludes zero (+0.204 [+0.028, +0.384]), which is what the file already
   says. These contrasts are paired within (target, seed), so the row bootstrap was nearly right
   here — the widening that destroyed the AMBER interval does not apply to them.
7. **`CSTEER_FINDINGS` §3 — the endpoint-substitution trap.** The −0.599 Å and −0.164 Å "effects"
   are correctly diagnosed as arms walking to step 1.0 and landing on a structure the pipeline
   already has, and are correctly barred from quotation. Verified: `fit → avg` and `proj → avg`
   both pick step 1.0 on 5/5 folds and both land at 3.048 Å. The general rule adopted from it
   ("a selector on the boundary of its ladder is choosing an endpoint") is sound and fired again
   in `divselect`.
8. **`CSTEER_FINDINGS` §4 — the projection costs 0.165 Å.** `avg` 3.048, `proj` 3.213. Reproduces.
   Note it also means the pre-projection average is **0.156 Å better than the incumbent's 3.204**.

---

## 7. THE PATTERN, AND WHAT IT IMPLIES FOR THE NEXT AUDIT

Three of the five claims fail the same way, and it is the way RETRACT named:

| claim | quantity measured | quantity the argument needs |
|---|---|---|
| csteer §2 | mean of **signed** c over targets | **\|c\|**, because the law is `√(1 − c²)` and `t*` carries the sign |
| diversity §4 | **free-superposition** member RMSD | the **common-frame** member error, because that is the identity's term |
| fusion law (Sprint 15, retracted) | **arithmetic** mean of member errors | the **quadratic** mean |

Each is arithmetically correct, transcribed correctly, and reported against the wrong reference.
**AUDIT's rule — every arm carries a zero-information reference and a matched random control —
is necessary and, on this evidence, not sufficient.** A fourth check is needed and is proposed
here:

> **Before a statistic is quoted as evidence for a law, derive which functional of it the law
> consumes, and quote that functional.** Where a law is quadratic, sign-invariant, or
> frame-dependent, the arithmetic mean of the raw per-target quantity is not it.

The two remaining failures are different species and both are cheap to prevent:

> **A mechanism sentence needs a control of its own.** "Because the errors are compensating"
> (steer §3) required a matched-distance interpolation control and did not have one; the control
> took twenty minutes and reversed the mechanism while leaving the result standing.

> **An exhibit that can be derived on random point clouds is not an exhibit.** Before publishing
> an empirical identity, run it on inputs from outside the domain. `csteer` §1 and Sprint 15's
> fusion law both fail that test in one line.

---

## 8. LIMITATIONS OF THIS AUDIT — stated so it is not over-read

- **The two 126-target flagship modules were not re-run end to end.** `steer.py`'s and
  `csteer.py`'s fits were read from their artefacts. If a fit is wrong, I would not see it. What I
  re-derived independently are the aggregations, the interpolant, the cosine statistics, and every
  control. The `avg`, `med`, helix/strand/PPII/breathing directions and the `avg → native`
  residual in §1.3 **were** recomputed from source and reproduce `csteer.json`'s `med` cosine to
  three decimals, which cross-validates that pipeline.
- **§1's ORACLE ceilings are ceilings.** They require the sign of `c` and `||r||`, both native.
  Nothing in §1 is deployable and none of it contradicts the pivot's operational verdict.
- **§2.3's control uses retrieval windows, not fits.** Both are ideal-geometry-buildable
  conformations of the same peptide and the comparison is matched on angular distance, but they
  are not the same object. The flagship arm's own distance-binned column (which *is* the same
  object) shows the identical monotone trend, so the conclusion does not rest on the control alone.
- **n = 9 remains n = 9.** My recomputation of the integration widens intervals correctly; it
  cannot make nine targets informative. The AMBER question is open, not answered negatively — the
  right response is to measure it on the 126-target instrument, not to relitigate these nine.
- **I did not audit** QPHASE, ENERGY, REPAIR, `audit_align*`, `retract_*`, `divselect`, or the
  literature review, beyond reading their claims for cross-consistency with the five assigned.
- **The sealed 60-target benchmark was not read, listed, or referenced by any module here.**
