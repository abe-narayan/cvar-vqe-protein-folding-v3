# AGENT D — THEORETICAL / ADVERSARIAL / RESEARCH AUDITOR
## Sprint 19 findings

Pre-registered in `s19/PREREG_D.md`, written before any result was inspected and not edited
after. Every claim carries a §10 label. **TARGET is the unit** throughout; pair-level numbers are
labelled diagnostics. The sealed 60-target benchmark was not read, probed, or derived from, and
`results/benchmark_manifest.json` was not opened.

Artefacts, each with a `COMPLETE` flag:

    s19/results/D_P1P2P3_modality/     the census, harm, coherence          n=126
    s19/results/D_P2P3_hardened/       min-of-N nulls, matched cells        n=126
    s19/results/D_P6_extract/          extractability of the mode signal    n=126
    s19/results/D_P5_headroom/         the structural headroom test         n=126
    s19/results/D_T3_tangent/          tangent decomposition, SUPERSEDED    n=126
    s19/results/D_T3_tangent2/         the same, corrected projection       n=126

Code: `s19/d_modality.py`, `s19/d_modality2.py`, `s19/d_extract.py`, `s19/d_headroom.py`,
`s19/d_tangent.py` (superseded), `s19/d_tangent_report.py`, `s19/d_tangent2.py`.

`D_T3_tangent/` is kept rather than deleted: its absolute values are not quotable (see §2.1) and
its supersession is part of the record.

---

# 0. EXECUTIVE SUMMARY

> ## The coordinator's opening hypothesis is refuted in its mechanism, and the refutation localises where the Sprint-18 phenomenon actually lives.
>
> ## The distogram's errors ARE spatially coherent — and they are coherent to the same degree on the pairs where the model is unimodal and confident. Multimodality is not the carrier.

Seven results, in the order they were produced.

| # | statement | label |
|---|---|---|
| D1 | The multimodality census reproduces at n=126 but was overstated ~2.4×: 0.241 by mass-local-max, **0.097 under a prominence gate**, 0.155/0.165 after a 2× re-binning. My bin-edge attack **failed** on the axis I predicted — the mass→density correction changes nothing (0.240). | ESTABLISHED |
| D2 | Multimodality is **not harmful at the aim point**. The naive +0.410 Å penalty is entirely a (separation, sd) confound; under a nearest-neighbour match it is **+0.007 [−0.291, +0.316]**. | ESTABLISHED (the harm REFUTED) |
| D3 | The multimodality is **not coherent across pairs**. Multimodal pairs are scattered (adjacency 0.072 vs a stratified permutation null of 0.070), and the correction they imply has **zero** spatial coherence (+0.003 [−0.001, +0.008]). | REFUTED (the coordinator's mechanism) |
| D4 | **The dissociation.** The real residual IS spatially coherent, +0.0385 [+0.0320, +0.0457] — and equally so on unimodal (+0.0429) and multimodal (+0.0489) pairs. | ESTABLISHED — the sprint lead |
| D5 | The modes carry **−0.445 ± 0.026 Å** of genuine ORACLE information per multimodal pair, above a matched min-of-N null averaged over 25 sign draws — and **none of it transfers** to any native-free rule, in-fold or out, including a split-half fit *inside the same fold*. | ESTABLISHED (the direction CLOSED by mechanism) |
| D6 | Structurally, on the deployed refinement instrument at n=126: a **perfect ORACLE mode-picker** lands at **3.472 Å** against the 3.048 Å the pipeline already builds. Against a min-of-N-matched null the mode information is **−0.053 [−0.153, +0.044] — NOT MEASURED**. Both native-free mode rules are null-to-worse. | ESTABLISHED — the branch has no headroom toward the mission |
| D7 | **Why a realisable error is harmful, with a scale.** The whitened distance Jacobian has rank **exactly 2n − 5**; a generic direction fills that subspace to **0.330** (measured) and the distogram's real error field fills it to **0.695**, against **0.397** for a sign-randomised version of the *same* field: **+0.298 [+0.262, +0.331], 121/5**. The fit denoises the orthogonal component and cannot touch the in-manifold one. | EXACT (the rank) + ESTABLISHED (the fractions) + **SUPPORTED, not promoted** (the design consequence) |

And one identity that voids a coordinator counter-argument:

| T1 | `Distogram._risk[·,x] = Σ_b p_b\|x − c_b\|·w` is a non-negatively-weighted sum of absolute values, hence **convex in x with minimiser the weighted median of `CENTRES`**. The selection path consumes the full distribution in the letter and is **structurally incapable of expressing multimodality**. "Why did it not show in selection?" is not evidence either way. | **EXACT** — a theorem, not a discovery |

---

# 1. PRIORITY 1 — THE ATTACK ON THE COORDINATOR'S HYPOTHESIS

## 1.1 The census (Block P1). My predicted attack failed; a different one lands.

I expected the 21.7% to be a binning artefact. `BIN_EDGES` are non-uniform — 0.5 Å below 8 Å
rising to 4 Å above 19 — so under any smooth density, bin **mass** ∝ bin **width**, and a wide bin
can be a local mass maximum with no density maximum at all. **That attack fails.** At n=126, over
all 8,549 pairs:

| detector | multimodal fraction (target-level mean) | 95% CI |
|---|---|---|
| mass local maximum, ≥2% mass (reconstruction of the coordinator's rule) | **0.241** | [0.219, 0.263] |
| the same on **density** `p/w` | **0.240** | [0.219, 0.261] |
| density + prominence gate (valley < 0.5× the smaller peak, ≥5% mass) | **0.097** | [0.087, 0.107] |
| density modes after merging bins pairwise, offset 0 | 0.155 | [0.139, 0.171] |
| density modes after merging bins pairwise, offset 1 | 0.165 | [0.149, 0.182] |

The n=12 figures reproduce at n=126 (mass-within-1 Å 0.634 vs the reported 0.641; entropy 1.44 vs
1.43; multimodal 0.241 vs 0.217), so the coordinator's measurement is sound. **But 60% of the
flagged pairs are shoulders, not modes.** They do not survive a prominence criterion and a third do
not survive a 2× coarsening of the bins. The defensible figure for pairs carrying a genuinely
separated second mode is **~0.10, not 0.22**. F-D1 as pre-registered (density < half of mass) did
**not** fire; the census is damaged by prominence and re-binning instead, and I record the failed
prediction as a failed prediction.

## 1.2 Is multimodality harmful? (Block P2a) — the harm is a confound.

Multimodal pairs carry mean sd 2.16 against 1.25 and are structurally the harder pairs. Four
controls, in increasing order of how much they are trusted:

| control | multi − uni, \|mean − true\| | reading |
|---|---|---|
| none (the naive comparison) | **+0.410 Å** | confounded |
| matched cells, separation × sd-decile, target-level | −0.174 [−0.408, +0.072] | **NOT MEASURED** (§9) |
| regression adjustment on separation × log sd | −0.422 [−0.595, −0.248] | model-dependent |
| **nearest-neighbour match** — same target, same separation, nearest sd, refusing \|Δsd\| > 0.25 | **+0.007 [−0.291, +0.316]**, n_t=121, W/L 69/52 | **the honest answer** |

**Verdict: once separation and spread are controlled, multimodality has no measurable effect on the
accuracy of the mean.** The harm is REFUTED; the benefit the regression adjustment suggested is NOT
SUPPORTED and I withdrew it from the coordinator's write-up myself. The premise "the deployed
objective aims at a value the model itself considers unlikely" is true as a statement about
probability mass and false as a statement about accuracy.

## 1.3 Is it coherent? (Block P3) — no, and the real error is.

The brief asserted, without evidence, that the collapse error acts *"coherently across correlated
pairs, because a whole region flipping between two conformer families moves many pairs the same way
at once."* Two independent measurements, both at n=126.

**Clustering.** Fraction of multimodal pairs sharing a residue with another multimodal pair, against
a **separation-stratified label permutation** (500 permutations per target, permuting the flag
within separation bins so that the greater multimodality and different connectivity of long-range
pairs are both preserved):

    observed 0.072   null 0.070   z = +0.269 [+0.079, +0.462]

An interval excluding zero at a z of 0.27 is a real but negligible tendency. **The multimodal pairs
are scattered, not clustered into regions.**

**Signed coherence.** If a region flips as a unit, the correction each pair needs must agree in sign
between pairs sharing a residue. Taking sign(nearest-mode-to-truth − mean) — the direction a
conformer flip would pull each pair:

    adjacent 0.540   non-adjacent 0.537   difference +0.003 [−0.001, +0.008]

**Zero.** F-D3 fires on both limbs. The mechanism as stated is REFUTED.

## 1.4 The dissociation — the part that is worth more than the refutation

The same measurement applied to the **real** residual, sign(d_true − mean):

| class of pair | adjacent − non-adjacent sign agreement | 95% CI |
|---|---|---|
| all pairs | **+0.0385** | [+0.0320, +0.0457] |
| both pairs multimodal | +0.0489 | [+0.0268, +0.0731] |
| both pairs unimodal | **+0.0429** | [+0.0336, +0.0537] |
| one of each | +0.0414 | [+0.0313, +0.0524] |

**The distogram's errors are spatially coherent, and they are coherent to the same degree where the
model is unimodal.** The property Sprint 18 found — errors worse than random errors of the same
magnitude — is real, and multimodality is not its carrier. Both limbs matter: a refutation that
only said "no" would have left the coordinator with nowhere to go.

> **The standing question this leaves: what makes the residual spatially coherent on pairs where
> the model is confident and unimodal?** That is where the Sprint-18 mechanism is.

## 1.5 Was the ensemble the cause? (Block P4) — SUPERSEDED, and why

P4 asked whether `Distogram.for_target`'s averaging of a dropout ensemble manufactures modes no
member believes. **I did not run it, and the reason is a measurement, not a shortage of time:**
P2a shows multimodality has no accuracy consequence and P3 shows it has no coherence, so the
provenance of the modes can no longer change any conclusion in this sprint. Recorded as **OPEN —
NOT MEASURED**, with the note that it would become live again only if some future arm found a use
for the mode structure. It is cheap (one fold's members, `predict_proba` per member) and I am
leaving the pre-registration standing rather than deleting it.

## 1.6 Is the information extractable? (Block P6) — no, and this closes the direction

The coordinator asked for the one number that bounds the whole direction. Per-pair
\|aim − d_true\| against the **deployed** debiased mean; target-level bootstrap, n=126, the pinned
5 folds. Features: the full 17-bin vector plus 14 shape summaries (sd, entropy, skew, kurtosis,
mass-within-1 Å, Bayes-median − mean, both mode counts, top-two mode offsets and masses,
separation). Label: the signed correction d_true − mean (ORACLE).

**On multimodal pairs:**

| arm | vs deployed | reading |
|---|---|---|
| ORACLE nearest mode | −1.095 [−1.201, −0.987] | the ceiling of mode selection |
| **min-of-N matched null** (same \|offsets\|, random signs, minimised the same way) | −0.611 [−0.710, −0.515] | 56% of the ceiling is pure selection over k candidates |
| **mode information proper** = ORACLE − null | **−0.484 [−0.570, −0.401]** | what the modes actually know — see the note below on the draw |
| leave-fold-out ridge, full shape | −0.069 [−0.177, +0.050] | spans zero |
| leave-fold-out GBM, full shape | −0.045 [−0.158, +0.077] | spans zero |
| leave-one-**target**-out ridge | −0.096 [−0.207, +0.021] | spans zero |
| **within-fold split-half GBM** | **+0.057 [−0.075, +0.190]** | spans zero, wrong sign |
| leave-fold-out with permuted labels | −0.003 [−0.033, +0.027] | the zero-information control sits at zero, as it must |
| in-fold GBM | −0.982 | **IN-SAMPLE. Memorisation. NOT quoted as a bound.** |

**A discrepancy in my own two scripts, resolved by averaging rather than by picking.**
`d_modality2.py` reported this quantity as **−0.414** and `d_extract.py` as **−0.484**. Both are
correct: the min-of-N null is a **single random-sign draw**, and the two scripts consumed different
positions of the same `stable_rng` stream. Averaged over **25 independent draws** the value is
**−0.445, across-draw sd 0.026, range [−0.480, −0.367]**. **Quote −0.445 ± 0.026.** Neither
single-draw figure should be cited as *the* number, and the across-draw spread (0.026) is the right
measure of how precisely a single-draw min-of-N null pins this quantity.

**Zero of the −0.445 Å transfers.** The split-half arm is the decisive one: fitted and evaluated
inside the *same* fold, on the *same* distribution, out of sample — still null. So this is not
cross-fold distribution shift; the relationship is not there to learn. And the two obvious
native-free rules are *worse* than doing nothing — measured on the multimodal subset in
`D_P2P3_hardened`, against the deployed mean: **highest-density mode +0.315, uniformly random mode
+0.529.**

On **all** pairs, the only thing that transfers is a better separation/mean/sd map:
`lfo_noshape` (separation, mean, sd only) −0.117 [−0.212, −0.020]; adding the entire 17-bin
distribution moves it to −0.152 [−0.245, −0.053] and leave-one-target-out to −0.168 [−0.265,
−0.072]. **Distribution shape is worth ~0.035 Å over no shape at all, inside the noise.** And per
§5 a lower aggregate error metric that does not lower Cα-RMSD is a negative result — this is a
per-pair diagnostic and nothing more.

---

# 2. THEORY — two statements, labelled as identities

**T1 (EXACT).** `Distogram.__init__` builds `risk[·, x] = Σ_b p_b |x − c_b|` and stores
`risk · w`. For fixed non-negative `w`, this is a non-negatively-weighted sum of absolute values,
hence **convex in x**, and its minimiser is the weighted median of `CENTRES` under `p`. Verified
numerically: `grid[argmin_x risk]` is used throughout as the Bayes-L1 aim point.

*Consequence.* The brief's counter-argument — "the selection path already consumes the whole
distribution, so if full-distribution consumption were the answer, why did it not show there?" —
is void. The selection path consumes the distribution only to compute **one robust location
estimate**. It cannot express multimodality, so its silence carries no information about the
hypothesis. This is the §10 hazard applied to the coordinator's own argument rather than a
workstream's.

**T3 — WHY A REALISABLE ERROR IS HARMFUL, AND THE FREE ANALYTIC NULL FOR κ.** *(numbers in §2.1)*

Agent A measured `kappa` — the fraction of a distance-error field that a real ideal-geometry
structure realises — as real 0.809, signflip 0.352, shuffled/iso 0.367, a real alternative
structure 0.990, with RMSD ordering the arms identically (+0.984 across the six arm means;
per-target Spearman > 0 on 116/126). **That statistic has a free analytic reference and nobody had
to run a control for it.**

The weighted fit is, to first order, an orthogonal projection of the whitened residual onto the
**tangent space of the ideal-geometry manifold** at the start point. A **generic** error vector
retains a fraction ≈ dim(tangent)/dim(pair space) of its energy under that projection, carrying no
information at all.

    the fit is a DENOISER against the component of the error ORTHOGONAL to the manifold,
    and has NO power whatsoever against the component INSIDE it.

That is the whole of "the errors are worse than random errors of the same magnitude": the
distogram's error has more of its energy inside the manifold than chance requires, so the fit
*realises* it instead of averaging it away. It also explains Agent A's whitening result — flipping
residual signs pair-by-pair destroys the cross-pair correlation that puts a vector into a
low-dimensional tangent space — and why her "real alternative structure" control is the worst arm
on the board: that residual lies almost entirely inside the manifold by construction.

**§10 boundary, stated before the numbers.** "A realisable error gets realised" is close to
definitional. The non-trivial empirical content is exactly two quantities: the real field's
in-manifold fraction, and the fact that a matched-magnitude incoherent field's fraction equals the
dimension ratio. The dose–response between them is an observation, not a theorem.

**A correction to my own first statement of the null.** I initially quoted D = 2n − 2 from Sprint
18's B5. That **over-counts**: the map (φ, ψ) → CA trace is not injective, since an ideal-geometry
CA trace of n residues has at most 2n − 5 internal coordinates (n − 2 virtual angles and n − 3
virtual dihedrals). The honest null is the **measured** `rank(J_w)/npairs`, which is what §2.1
reports; `(2n−2)/npairs` is quoted beside it only as the naive version and as the figure I sent the
coordinator before I checked it.

## 2.1 T3 measured — the explicit tangent projection

`s19/d_tangent.py`, `s19/d_tangent_report.py`; artefacts `s19/results/D_T3_tangent/`. Per target:
the finite-difference Jacobian `J = ∂d_pair/∂(φ,ψ)` at the deployed start point, whitened by
`1/sd`; the in-tangent energy fraction of each error field by explicit QR projection; the measured
`rank(J_w)`; and `fin_random`, the in-tangent fraction of 20 genuinely random directions in pair
space — **the empirical null, which needs no derivation and cannot fail the way my algebra did.**
No optimisation is run.

**THE DECISION RULE, written before the numbers and to be applied mechanically whichever way they
fall.** If the across-arm-means correlation between the start-point tangent fraction and Agent A's
realised κ is high **and** the per-target Spearman of in-tangent excess against the per-target
real−signflip gap is low, then **the statistic orders arms and not targets, the tangent projection
is a reference and not a model, and T3's design consequence stays at SUPPORTED.** It is promoted
only by per-target evidence.

### The rank is exact, and the measurement confirms it

    mean rank(J_w) = 20.92        mean 2n-5 = 20.92        mean 2n = 25.92

**The whitened distance Jacobian's rank is exactly 2n − 5 on every target** — the corrected
derivation, verified numerically rather than assumed. And the empirical geometric null lands where
that requires:

| quantity | value | 95% CI |
|---|---|---|
| `rand_white` — a direction isotropic in the **whitened** space | **0.330** | [0.318, 0.341] |
| `rank(J_w)/npairs` | 0.329 | [0.318, 0.340] |

Label the rank statement **EXACT**; 0.330 is a measurement that agrees with it.

### Every control sits at its own correct null, and the matched contrast is large

| field | in-tangent energy fraction | 95% CI |
|---|---|---|
| **real** | **0.695** | [0.669, 0.722] |
| signflip — identical per-pair whitened magnitudes, cross-pair signs randomised | **0.397** | [0.373, 0.420] |
| shuffled | 0.586 | [0.565, 0.608] |
| iso | 0.577 | [0.553, 0.601] |
| coherent — a residual from a real alternative structure | 0.846 | [0.822, 0.870] |
| `rand_raw` — isotropic in **raw** distance space, then whitened | 0.581 | [0.567, 0.595] |
| `rand_white` — the geometric null | 0.330 | [0.318, 0.341] |

    real - signflip   +0.298 [+0.262, +0.331]   W/L 121/5
    real - rand_white +0.365 [+0.338, +0.392]   W/L 125/1
    iso  - rand_white +0.247 [+0.226, +0.268]   W/L 121/5   (the whitening anisotropy)

`iso` (0.577) and `shuffled` (0.586) are at **`rand_raw` (0.581)**, not at the geometric null —
because they are isotropic in *raw* distance space, and after whitening by `1/sd` that is not
isotropic: they carry more energy where `sd` is small, which is where the Jacobian is large. Their
correct reference is 0.581. **The distogram's error field puts 0.695 of its whitened energy inside
a subspace that a generic direction fills to 0.330 and a sign-randomised version of the same field
fills to 0.397.** That is the mechanism with a scale, measured rather than derived.

### THE DECISION RULE FIRES. The design consequence is NOT promoted.

    across the 5 ARM MEANS, pearson(start tangent, A's realised kappa) = +0.889      HIGH
    per-target rho(in-tangent(real) - in-tangent(signflip), gap) = +0.134, p = 0.135,
                                                        95% CI [-0.057, +0.301]     LOW
    per-target rho(in-tangent(real) - rand_white,          gap) = +0.160, p = 0.074
    per-target rho(in-tangent(real), raw,                  gap) = +0.157, p = 0.080

**The statistic orders arms and not targets. The tangent projection is a reference and not a
model. T3's design consequence stays at SUPPORTED.** The correction moved the per-target
correlation from +0.143 to +0.134; it did not rescue it, and I would not have promoted on +0.16
either.

### One dissociation I cannot explain — recorded as OPEN

The linear projection tracks Agent A's realised κ per target well (ρ = +0.656 for `real`), yet
**neither predicts the per-target gap the way κ itself does**:

    rho(A's kappa(real),                   gap) = +0.360   p = 3.4e-05
    rho(A's kappa(real) - kappa(signflip), gap) = +0.482   p = 1.1e-08
    rho(my in-tangent(real),               gap) = +0.157   p = 0.080

So the **post-fit** quantity carries per-target information about the outcome that the start-point
geometry does not. Something happens along the ~5-radian trajectory that neither lane is measuring.
**Caveat volunteered before anyone leans on +0.482:** κ and the gap are both derived from the same
two fits, so a common cause — targets where the fit converges cleanly having both a low objective
and a low RMSD — is not excluded. It is suggestive, not clean; a clean test needs a native-free
predictor of κ.

### Two wrong nulls I published before measuring, both struck

I sent the coordinator `(2n−2)/P = 0.381`, then corrected it to `(2n−5)/P = 0.329`, and on the
strength of the first claimed that Agent A's incoherent κ values "sat on the null". Neither
algebraic figure should be quoted: the first over-counted the rank, and the comparison was between
a *post-fit* κ and a *linear* null in the first place. Separately, my own first implementation
projected with `np.linalg.qr(J_w)`, whose reduced form returns 2n orthonormal columns even when
`J_w` is rank-deficient — so every fraction in the first run was computed against a 2n-dimensional
subspace. **The tell was the empirical null coming out at 0.417 ≈ 2n/npairs instead of 0.329**, and
it was the measurement, not the review, that caught it. Contrasts in that run remained valid
because the random control shared the inflated subspace, but no absolute number from
`D_T3_tangent/` is quotable; `D_T3_tangent2/` supersedes it and both are kept.

**T2 (EXACT, and it bounds the whole family).** Any per-pair-separable objective `Σ_pq f_pq(d_pq)`
— which includes the deployed `((d−d̂)/sd)²`, the Bayes risk, and any per-pair `−log p(d)` — cannot
represent the *joint* constraint that different pairs select the same conformer. A non-convex
per-pair `f` can at least let the geometry settle into a self-consistent set of modes, which the
convex forms cannot; that is an empirical question, and P5 bounds it from above with an oracle.
The corollary for design: **if the goal is joint conformer selection, the change must be in the
representation or the decoder, not in the per-pair loss.** This is the direct route to §4's
recommendations.

---

# 3. PRIORITY 3 — AUDIT

## 3.1 Two defects in the coordinator's live arm, `s19/distobj.py`

Reported to the coordinator during the run, not saved for this document.

**(a) The docstring's weighting claim is false for two of four arms.** `_fit`'s docstring states
"Every branch keeps the SAME weights `1/sd^2`". For `moment`/`mode`/`unimodal`, `_target_terms`
returns `r*r` with `r = (d − dhat)*inv`, so `1/sd²` is embedded. For `nll` and `risk` there is no
`inv` anywhere and `fg` does `f = v.sum()`. **Those two arms change the target and the effective
weighting at once**, and Sprint 18 priced the `1/sd²` weighting at +0.153 [+0.066, +0.247] (G6).
They cannot be read as "the objective consumes the full distribution". `mode` is unaffected and the
coordinator's reason for making it the primary arm survives intact. This is the same class of
defect as Sprint 18's `objceil.py` line 163.

**(b) `nll` is a negative log bin MASS, not a negative log density.** With non-uniform `BIN_EDGES`,
`−log p_b = −log density_b − log w_b`, so the loss carries a spurious `−log w_b` reward for wide
bins spanning **2.079 nats** — a systematic outward pull toward the coarse long-distance bins, in
exactly the direction the shipped separation debias exists to remove. One-line fix:
`logP = np.log(np.maximum(P[m] / WIDTHS[None], 1e-9))`.

**(c) The same issue reaches `mode`.** `modality()` takes `mode = CENTRES[np.argmax(P, axis=1)]` —
the modal bin **by mass**. Measured at n=126 over all 8,549 pairs: mass-argmax and density-argmax
disagree on **14.0%** of pairs, and where they disagree the mass mode sits **+2.72 Å further out**;
the mean shift over all pairs is **+0.380 Å outward**. The distogram already over-predicts distance
(+0.509 Å global, +1.492 Å at separation 11–15) and the coordinate average contracts 22–25%, so
this arm carries an extra expansion on top of the target change. If `mode` loses, some of the loss
is this. My P5 uses the density mode, which is the reason the two instruments could disagree, and
I have said so in advance rather than after.

**(d) Minor.** `pk` counts peaks over interior bins `1..15` only, so a mode in bin 0 (<4.5 Å) or
bin 16 (>23 Å) is never counted.

## 3.2 Sprint-wide leakage and seeding sweep — clean

Over every `s19/*.py` at the time of the sweep (lanes A, B, C, D and the coordinator):

* **no bare `hash()`** anywhere (the only matches are inside Agent B's own audit that greps for it);
* **no unseeded `np.random.*`** and **no `default_rng`** — every stream goes through
  `s15/seed.stable_rng`;
* **no reference** to `benchmark_manifest`, `bench60` or `benchmark60` in any lane;
* native quantities (`nat_ca`, `dtrue`, `rr`) appear only in scoring and in arms that are
  **declared ORACLE diagnostics in their own module docstrings** — checked line by line in
  `a_coh.py`, `a_struct.py`, `agentC_lib.py`, `agentC_kv.py`.

Two specific checks against known traps:

* **Lane A's error-structure arms carry their weights correctly.** `a_coh.py`'s `signflip`,
  `shuffled` and `iso` fields modify the *target distance vector* only; `w` is computed from `sd`
  and is never permuted alongside. `signflip` keeps each residual attached to its own pair, which
  is the right construction after Sprint 18's G6a ("orphaning a residual from its weight is bad").
* **Lane A's `signflip` is not exactly magnitude-matched, by 4.5%.** `a_coh.py` applies
  `fld = np.maximum(field, 2.0)` *after* constructing the sign-flipped field, and computes the
  arm's residual RMS after that clip. Real 3.205 Å, signflip 3.062 Å — sign-flipped residuals that
  push a target below the 2.0 Å floor are truncated. Small against a −1.242 Å effect and it moves
  the arm in the *conservative* direction (less residual), but it should be stated rather than
  left for a reviewer to find. `coh_scaled` is exactly matched (3.205 vs 3.205), as designed.
* **Lane A's `kappa` is a LOWER bound, which strengthens rather than weakens her claim.** It is
  `1 − f/E_t` where `f` is `align_lib.fit`'s returned objective from a **single start**, so it
  measures what the fit *did* realise, not what is realisable. The real field's 0.809 is therefore
  a floor.
* **An attribution correction I made and then had corrected.** I flagged an inverted sign in a
  claim about realisability and harm. The inversion was in the coordinator's *relay* only; Agent A's
  original wording ("unrealisability is anti-correlated with harm") was correct and logically
  identical to mine, and `s19/LEDGER.md` L3 never carried the error. Recorded here because a
  mis-attributed catch is itself a defect, and this lane should not be exempt from that.
* **Lane C states its identity as an identity.** `agentC_lib.py` derives
  `readout² = E_mem² − D²` from the parallel-axis theorem, labels it **EXACT**, and explicitly
  warns that how a gate moves the two right-hand terms is *not* a theorem. It also documents the
  free-superposition/common-frame confusion that made a prior sprint's `FRAME²` come out negative.
  This is §10 being enforced by the lane on itself; no action needed.

## 3.3 My own instrument's validity gate

`s19/d_headroom.py`'s `raw` arm reproduces `s17/refine.py`'s `refine_full` **bit-exactly**
(max |Δ| = 0.0000). `sd` is passed identically to every arm. The leave-fold-out separation debias
is fitted on the full `I.targets()` and only then applied per target — §7's subset trap. Every
random stream is `stable_rng(pdb, tag)`.

**Engineering observation for a future lane.** `s19/d_tangent.py` recomputes `I.project` per
target — 2–4 s of multi-start `lam_path` — where Agent A already caches exactly that start in
`s19/cache/start_<pdb>.npz` via `a_fit.start`. That start is deterministic and identical for every
arm and every lane in this sprint; any module needing it should read A's cache rather than rebuild
it. My run cost roughly 8 contended minutes it did not need to.

**Compute discipline, recorded against myself.** I killed my own P5 replicate at 70/126 when the
box was at 91% CPU with 7 python processes, and removed its directory so no partial file could be
mistaken for a result. §15 allows one heavy process per workstream and the replicate was my second.
The replicate's entire content — whether `stable_rng` gives the same stream in a fresh process —
was then tested directly at n=126: **630 arrays compared, 0 mismatches**, checksum
176466.6545295715, and distinct streams per target confirmed (a constant stream would itself have
been the defect). The stronger replication was already on the board and is cross-implementation:
`raw` reproducing `s17/refine.py`'s `refine_full` at max |Δ| = 0.0002 Å on 126/126.

**One process-hygiene note against myself:** the first P5 launch was detached with `nohup … &` and
was killed at 70/126 when its parent shell exited; the partial file was correctly flagged
`complete: false` and was never read as a result. Relaunched under the task manager. Recorded
because §7 says never treat a partial file as complete, and the guard worked.

---

# 4. PRIORITY 2 — LITERATURE

## 4.1 `arXiv:2609.02113` dissected — *Logarithmic-scale VQE for off-lattice protein structure prediction in continuous torsional angle space* (Cumbo, Raubenolt, Puram, Katzenmeyer, Joshi, Blankenberg; Cleveland Clinic; submitted 2 Sept 2026)

I read the paper in full rather than the abstract. Seven findings.

**(i) Its reported RMSDs are NOT comparable to ours, and the gap is not small.** The headline
0.623 Å for chignolin is the **best retained snapshot**, selected by RMSD to the native, out of
**7,332,100 retained snapshots across 2,333 final models** (their Table 1) on **two** targets. Even
their "Snapshot RMSD median" column is a median over replicas of a **per-replica ORACLE best**
("the lowest-effective-RMSD snapshot from each replica"). The only genuinely native-free numbers in
the paper are the final-model medians:

    chignolin (5AWL, 10 res)   2.90 / 4.20 / 3.84 Å   (custom / Rosetta / OpenMM)
    Trp-cage  (2JOF, 20 res)   5.39 / 6.59 / 6.21 Å

Against our incumbent of **3.204 Å** as a mean over **126 cluster-disjoint 9–16mers**, their
native-free median on a single 10-residue peptide is 2.90 Å and on a single 20-residue peptide
5.39 Å. **We are not behind this paper.** Anyone in this programme who quotes 0.623 Å as the
state of the art is quoting an oracle min-of-N with no min-of-N null.

**(ii) There is no cluster-disjointness discipline, and the two targets are the two most
over-studied peptides in existence.** Chignolin and Trp-cage are the standard validation systems
for `amber99sb-ildn` (used for their final GROMACS minimisation), for Rosetta's all-atom energy,
and for a decade of folding-simulation papers. There is no train/test split because there is no
training — but the *force fields* have been tuned on exactly these systems, which is the
force-field analogue of leakage and is not discussed.

**(iii) The quantum component performs neither quantum optimisation nor a quantum energy
evaluation.** The paper is explicit and to its credit transparent: *"Unlike electronic structure
VQE, the total scalar score is not decomposed into quantum observables and measured as the
expectation value of a qubit Hamiltonian; it is computed classically from the reconstructed
Cartesian structure after the circuit readout and minimized by the classical optimizer."* So the
circuit is a **parameterised nonlinear map θ → torsions**; the energy is classical; COBYLA/SLSQP
does the optimising. In statevector mode the map is deterministic (a reparameterisation); on
hardware the CDF decoder makes it a **sampler**. It is a quantum-parameterised generator, not a
VQE in the usual sense. **This is the same architectural question this programme answers
experimentally under §13, and it is worth noting that a September-2026 paper in this area does not
answer it at all.**

**(iv) The O(log N) qubit claim is real but purchased with Θ(N) circuit depth, and the classical
problem is not compressed.** They set `n = max(2, ⌈log₂ N_torsion⌉)` and `d = ⌈N_torsion/n⌉ + 2`
repetition layers of EfficientSU2. Parameter count is `2n(d+1) ≈ 2N + 6n` — **linear in N**, and
entangling depth is `≈ d·n ≈ N`. For chignolin, 44 torsions → 6 qubits and ~10 layers, i.e. ~132
classical parameters for 44 degrees of freedom. So the register width scales logarithmically and
**everything that actually binds — variational parameter count, circuit depth, energy evaluations,
and (they say so themselves) shot count — remains linear or worse.** Their own Hardware Feasibility
section concedes this: *"this intensive shot requirement and the resulting noisy readouts
effectively trade the problem of qubit quantity for a problem of measurement complexity."*

**(v) Their energy model.** A custom hybrid Hamiltonian: Lennard-Jones, Coulomb, a
sigmoid-weighted SASA implicit solvent, explicit hydrogen bonding, local and non-local steric/clash
terms, Ramachandran-shaped harmonic/Gaussian backbone and χ₁ rotamer priors, ω trans priors, and
hard geometric-integrity penalties (L-chirality at Cα, peptide-plane twist, proline ring closure).
Alternates: Rosetta all-atom via PyRosetta, and OpenMM/AMBER-ff14SB. Final models are protonated
and minimised in GROMACS/amber99sb-ildn. **Read against our record:** their "geometric integrity"
term is doing the job our AMBER repair does, and their custom function beats both established force
fields on final-model quality — consistent with our finding that generic physics is a poor
*ranker* and a fine *validity operator*.

**(vi) Their limitations illuminate our position precisely, and this is the paper's real value to
us.** They report, independently and on a completely different architecture:

> *"Across all three energy functions, low-energy structures did not always correspond to the
> lowest-RMSD conformations, highlighting persistent energy-ranking imbalance."*
> *"the scoring objectives do not consistently provide enough resolution near the basin floor to
> select those conformations as final endpoints."*
> Median RMSD recovered by snapshot mining: **2.18 Å (custom), 2.66 Å (OpenMM), 3.26 Å (Rosetta)**.
> *"the primary bottleneck is [not sampling but] final-model selection."*

That is **`search saturates; discrimination binds`, reproduced externally**, with a selection gap
of 2.2–3.3 Å against our measured ~1.9 Å. It is independent corroboration of the programme's
central negative result, from a group that was not trying to prove it, using force fields and a
generator that share nothing with ours. **It should be cited in any write-up of that finding.**

**(vii) THE FINDING TO LEAD WITH — external replication of this programme's central negative.**
Independently, on an architecture sharing nothing with ours — a quantum-circuit-parameterised
torsion generator, three all-atom force fields, no retrieval, no distogram, no machine learning —
this group reports:

> *"Across all three energy functions, low-energy structures did not always correspond to the
> lowest-RMSD conformations, highlighting persistent energy-ranking imbalance in the sampled
> landscapes."*
>
> *"the scoring objectives do not consistently provide enough resolution near the basin floor to
> select those conformations as final endpoints."*
>
> *"the logarithmic quantum generator is highly capable of sampling the correct conformational
> space: the primary bottleneck [is] final-model selection."*

Their measured selection gap — median per-replica improvement from picking the best visited
structure rather than the endpoint — is **2.18 Å (custom), 2.66 Å (OpenMM), 3.26 Å (Rosetta)**,
against this programme's measured ~1.9 Å.

**This is `search saturates; discrimination binds`, reproduced from outside.** It is the strongest
evidence available that the finding is a property of the *problem* and not of our pipeline, because
it comes from a group with no contact with this work, who were not trying to demonstrate it, and
who report it as an obstacle to their own claim rather than as a result.

**A third, independent corroboration of a second finding.** `arXiv:2606.21241` (*Assessing Cost
Hamiltonian Reliability in Quantum Protein Structure Prediction*, June 2026) tests the quantum
folding community's own contact-energy Hamiltonian directly against RMSD and concludes that *"the
energy landscape of the considered cost Hamiltonian is not correlated well enough to the actual
error to provide meaningful predictions"* for small peptides — its ground state does not correspond
to the native. That is this programme's **`the objective does not rank the native`** (native at the
36.8th percentile, argmin on 3/126) arrived at independently, and its stated methodological lesson —
*"the importance of studying the reliability of the cost Hamiltonian independently from the quantum
approach used"* — is §13 of our own brief. Their one positive note, that correlation improves with
larger instances and more interaction shells, is consistent with our record that these signals are
weakest at peptide length.

**(viii) One structural detail worth stealing conceptually, and one to avoid.** Worth noting: their
end-to-end differentiability through NeRF reconstruction, so that Cartesian-space gradients
backpropagate to generator parameters — the same trick as our `pj._torsion_grad`, independently
arrived at. To avoid: their final models are systematically **expanded** (median end-to-end 7.84 Å
against an experimental 5.51 Å for chignolin, 16.78 Å under Rosetta), the mirror image of our
coordinate average's 22–25% **contraction**. Two architectures, opposite compactness pathologies,
both uncorrected — a radius-of-gyration prior is apparently a general unclaimed lever, and this
programme has already measured that one (`forensics_rggate`).

## 4.2 What else is genuinely usable

**External calibration for our length range — the single most useful item.** A 2026 benchmark of
AlphaFold2, RoseTTAFold2, ESMFold, OmegaFold and DMPfold2 on curated experimental short-peptide
structures (10–49 aa) reports that **accuracy improves systematically with length**, that all
models do markedly worse on β-rich and disordered peptides than on helical ones, and — decisively
for us — that **AlphaFold2's training excluded NMR structures and peptides shorter than ~16
residues.** Our 126 targets are 9–16 residues and largely NMR. That is a documented, external
reason why an off-the-shelf predictor is not a trivial upper bound here, and it should be stated in
any paper this programme writes rather than left implicit.

**The architectural lesson from AF2/trRosetta, which bears directly on T2.** trRosetta and its
predecessors did exactly what this programme does — predict inter-residue distance distributions,
convert them to restraints, and minimise. The literature on that route is candid about its failure
mode: gradient descent on a distance potential *"will not generally locate the global minimum"*, is
*"sensitive to initial conformations"*, and *"the common problem of multiple minima persists even
for a reasonably accurate distance potential"*; and distance-only objectives *"struggle to enforce
higher-order geometric consistency"* so that *"a predicted distance map could in principle
correspond to a physically implausible structure"*. **AF2 abandoned this route.** In AF2/AF3 the
distogram survives only as an auxiliary head and training loss; the structure comes from a learned
decoder (IPA + FAPE) over the pair representation. Distance-AF, which retro-fits user distance
restraints onto AF2, does so *inside* the structure module rather than by least squares outside it.
**Our Sprint-17/18 result — that the objective's optimum is worse than the structure the pipeline
already builds — is the published failure mode of the pre-AF2 generation, reproduced.**

**Conformer generation.** Torsional diffusion (Jing et al.) diffuses on the hypertorus of torsion
angles with an extrinsic-to-intrinsic score model, generating conformer *ensembles* rather than a
point estimate, orders of magnitude faster than Cartesian diffusion. Internal-coordinate diffusion
has since been applied to macrocyclic and mixed-chirality cyclic peptide ensembles at
simulation quality (2025–2026). This is the mature version of the generator this programme builds
by retrieval, and it produces a *distribution*, which is what our operator (coordinate averaging)
actually consumes.

**Distance geometry theory.** Recent EDM work is moving away from rank-constrained SDP toward
(a) Riemannian optimisation over PSD Gram matrices, where non-negativity and the triangle
inequality are enforced *implicitly*, and (b) **hierarchical Bayesian priors placed directly on the
latent point set that generates the EDM** (arXiv:2601.22765, Jan 2026). The second is exactly the
structural answer to T2: a prior on the point set is *not* per-pair separable, so it can express
joint constraints that no reweighting of `((d−d̂)/sd)²` can.

**Model quality assessment.** The CASP record is that **consensus methods beat single-model
methods**, and that single-model correlation degrades precisely as the pool gets worse. Our
programme independently found the same shape (`consensus is the only in-band discriminator`;
`consensus is outlier avoidance, not a nativeness signal`). The best 2025 EMA numbers (Pearson
~0.75 on CASP15 with graph transformers) are on *protein-sized* targets with far more structural
signal than a 12-residue peptide offers; nothing in that literature suggests a single-model scorer
would beat our measured in-band ceiling.

**Quantum — the third external corroboration.** `arXiv:2606.21241` is covered in §4.1(vii).
QuPepFold (PLOS One, 2026) packages CVaR-VQE for peptide folding and reports ~30%
faster convergence than expectation-value VQE and >90% fidelity on IonQ Aria-1 — a *convergence*
claim on a lattice-style Hamiltonian, not an accuracy claim, and **without a matched classical
baseline of the kind §13 requires**. The earlier CVaR-VQE-vs-MD comparison claims better global
optimisation than MD, again without matched budget. **Nothing in the current quantum-folding
literature would survive this programme's own control discipline**, which is worth stating plainly:
our closure of the quantum branch (Sprint 18 C1–C5) is, as far as I can find, the most controlled
negative result in the area.

## 4.3 Sources

- Cumbo, Raubenolt, Puram, Katzenmeyer, Joshi, Blankenberg, *Logarithmic-scale variational quantum
  eigensolver for off-lattice protein structure prediction in continuous torsional angle space*,
  [arXiv:2609.02113](https://arxiv.org/abs/2609.02113) (2 Sept 2026) — dissected in §4.1.
- *Assessing Cost Hamiltonian Reliability in Quantum Protein Structure Prediction*,
  [arXiv:2606.21241](https://arxiv.org/pdf/2606.21241) — the cost Hamiltonian does not track RMSD.
- *A Comprehensive Evaluation of Protein Structure Prediction Models for Short Peptides*,
  [bioRxiv 2026.07.02.736085](https://www.biorxiv.org/content/10.64898/2026.07.02.736085v1) —
  10–49 aa, five models, accuracy rising with length, AF2 trained without NMR or <16-residue
  peptides.
- *AlphaFold-based peptide structure prediction: opportunities, limitations, future directions*,
  [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0734975026001837);
  *Benchmarking AlphaFold2 on peptide structure prediction*,
  [Structure 2022](https://www.sciencedirect.com/science/article/pii/S0969212622004798).
- Jing et al., *Torsional Diffusion for Molecular Conformer Generation*,
  [arXiv:2206.01729](https://arxiv.org/abs/2206.01729); internal-coordinate diffusion for
  macrocyclic peptides, [arXiv:2305.19800](https://arxiv.org/pdf/2305.19800); cyclic-peptide
  ensembles via diffusion, [PMC13262345](https://pmc.ncbi.nlm.nih.gov/articles/PMC13262345/).
- *Distance-AF: modifying AlphaFold2 models with user-specified distance constraints*,
  [Commun. Biol. 2025](https://www.nature.com/articles/s42003-025-08783-5) — restraints applied
  *inside* the structure module, not by external least squares.
- *Distance-guided protein folding based on generalized descent direction*,
  [Brief. Bioinform. 22(6)](https://academic.oup.com/bib/article/22/6/bbab296/6341661) — the
  multiple-minima failure mode of distance-potential gradient descent.
- *Bayesian Matrix Completion Under Geometric Constraints*,
  [arXiv:2601.22765](https://arxiv.org/html/2601.22765v1) — priors on the latent point set, the
  structural answer to T2; *Provable Non-Convex Euclidean Distance Matrix Completion*,
  [arXiv:2508.00091](https://arxiv.org/abs/2508.00091).
- CASP model-accuracy-estimation record: consensus over single-model,
  [CASP11 assessment](https://ncbi.nlm.nih.gov/pmc/articles/PMC4781682); GATE at ρ=0.748 on CASP15,
  [bioRxiv 2025.02.04.636562](https://www.biorxiv.org/content/10.1101/2025.02.04.636562v2.full).
- QuPepFold, [PLOS One 2026](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0342012);
  CVaR-VQE vs MD, [Int. J. Biol. Macromol. 2024](https://www.sciencedirect.com/science/article/abs/pii/S0141813024038388);
  utility-level QPSP, [arXiv:2506.22677](https://arxiv.org/abs/2506.22677).

---

# 5. PRIORITY 1, CONTINUED — THE STRUCTURAL HEADROOM TEST (Block P5)

`s19/results/D_P5_headroom/`, n=126, complete. Instrument = `s17/refine.py`'s `refine_full`
exactly. **Only the target vector changes between arms; `sd` is passed bit-identically to every
one.** Validity gate: `raw` reproduces `s17` refine_full at max |Δ| = **0.0002 Å**, 126/126 within
0.01.

| arm | RMSD | median | vs raw [95% CI] | W/L | objective | \|shift\| |
|---|---|---|---|---|---|---|
| avg (the start) | **3.048** | 2.837 | — | — | — | — |
| proj | 3.213 | 2.966 | — | — | — | — |
| **raw** (deployed) | 3.610 | 3.370 | +0.000 | — | 56.5 | 0.000 |
| med (Bayes-L1 median) | 3.627 | 3.535 | +0.018 [−0.030, +0.073] | 68/58 | 66.7 | 0.402 |
| mode1 (density mode) | 3.644 | 3.431 | +0.035 [−0.017, +0.087] | 54/72 | 72.4 | 0.677 |
| **ORACLE_bestmode** | 3.472 | 3.299 | **−0.138 [−0.213, −0.058]** | 83/43 | 67.5 | 0.843 |
| **minofN_null** | 3.525 | 3.339 | −0.085 [−0.160, −0.009] | 80/46 | 68.3 | 0.763 |
| rand_flip | 3.629 | 3.503 | +0.019 [−0.052, +0.091] | 58/68 | 78.7 | 0.841 |
| ORACLE_true | 1.148 | 0.307 | −2.462 [−2.770, −2.146] | 117/9 | 3.0 | 2.347 |

    ORACLE_bestmode - minofN_null  (the structural MODE INFORMATION)
        = -0.053 [-0.153, +0.044]   spans zero, below the 0.084 A MDE   NOT MEASURED

**Reported exactly as pre-registered, including the clause that splits.** F-D5 said: if
`ORACLE_bestmode − raw` does not reach −0.084 Å with a CI excluding zero, no mode-selection scheme
can move this metric. **F-D5 does NOT fire** — −0.138 [−0.213, −0.058]. F-D5b, the control clause,
was written against `rand_flip` and does not fire against it either. Against the **tighter** control
that I also pre-registered as an arm — `minofN_null`, matched on the offset magnitudes *and* on the
min-over-k selection — the mode information is **NOT MEASURED**. I report both, name which clause I
wrote, and state my own view: `rand_flip` under-controls, because the oracle arm is a minimum over
k candidates and `rand_flip` is a single draw.

**The statement that settles it for the mission, and it needs no control at all.** A **perfect**
oracle mode-picker lands at **3.472 Å**. The structure the pipeline already builds is **3.048 Å**.
Even with the native in hand to choose every mode, refinement is **+0.424 Å worse than doing
nothing.** The branch has no headroom toward the mission whether or not the mode information is
real.

**Three cross-checks that come free.**

1. **The two instruments agree and the mode-definition defect was not load-bearing.** My `mode1` is
   the **density** mode; the coordinator's `distobj` used the **mass** mode. Mine gives +0.035
   [−0.017, +0.087]; his gave +0.037 [−0.009, +0.088]. Two definitions that differ on 14% of pairs
   by +2.72 Å each move the outcome by 0.002 Å. The defect was real and worth fixing; it did not
   change the answer.
2. **`ORACLE_true` reaches 1.148 Å with objective 3.0**, reproducing Sprint 18's 1.152 basin and
   G2's "perfect distances make this functional point at the native", on an independently written
   harness.
3. **The Bayes-L1 median — the aim point the *selection* path actually uses, by T1 — is worth
   +0.018 Å in the objective.** Transplanting selection's own consumption of the distribution into
   the objective buys nothing. That closes the last untested variant of the branch.

---

# 6. PRIORITY 4 — TWO RECOMMENDATIONS

Judged by expected information value. Both follow from T2/T3 and from the programme's own record,
not from novelty. Neither is proposed because it is advanced.

**The framing both inherit, from T3.** If the harm is the component of the predictor's error that
lies *inside* the ideal-geometry manifold, then **no change to the loss can fix it**: any objective
whose argmin lies on the manifold realises whatever part of the error lies in the manifold. That
leaves exactly two levers — **(1) reduce the in-manifold component of the predictor's error**, and
**(2) do not consume the argmin at all**. Everything else in the objective family is measuring a
projection whose input is already wrong in the one direction the projection cannot suppress.

## R0 — THE LEAD, HELD PENDING ONE MEASUREMENT: a training loss in the RIGHT BASIS

**The proposal.** The §5 compass says where error lands beats how much there is, and expresses that
in the basis of *pair confidence and separation*. T3 says the basis that actually matters is
**in-manifold versus orthogonal**. A predictor trained with a loss that up-weights the in-manifold
component of its own error — the component a downstream least-squares fit will realise rather than
average away — is optimising the quantity that reaches Cα-RMSD, which MAE is not.

**Why it is NOT recommended for implementation yet, and this is the coordinator's condition, which
I agree with.** If PairNet, a second MLP seed, and the retrieval pool's own error all carry the
*same* in-manifold component, then that component is an **identifiability limit of the inputs**, no
loss in any basis removes it, and the lever is dead before it costs anything. If it is
family-specific, it is the main event. Agent A is measuring exactly that attribution. **The design
work waits for it.** Recorded here so the reasoning exists before the answer does.

## R1 — SAMPLE THE OBJECTIVE'S POSTERIOR; DO NOT MINIMISE IT

**The argument.** The programme has measured, repeatedly and from different directions, that
concentrating is wrong here: *"running the VQE is WORSE than not running it"*; *"searching harder
hurts on a bad objective"*; the deployed objective's own optimum is +0.561 Å worse than the
structure the pipeline already builds (G1); and the largest single durable result in the torsion
work was **"do not collapse the posterior" (−2.253 Å)**. Meanwhile the operator that *does* work —
coordinate averaging — is a crude posterior mean over a candidate set. **The architecture is
already doing approximate posterior inference and then throwing it away at the last stage by
running an argmin.**

**First experiment (cheap; reuses existing gradients).** On the `s17/refine.py` instrument, replace
the single L-BFGS fit with a short Langevin / SGLD sample from `exp(−E/T)` where
`E = Σ (d − d̂)²/sd²`, using `s15/align_lib.fit`'s exact gradient. Calibrate `T` **native-free**, by
requiring that the sampled per-pair distance spread match the distogram's own predicted `sd` — a
temperature chosen on the model's own uncertainty, not on RMSD. Then apply the operator the
programme has already validated: the coordinate average of the sample, AMBER-relaxed at k=30.
Arms: `argmin` (= refine_full, 3.610), `posterior mean at calibrated T`, a **temperature ladder**,
and the coordinate average (3.048).

**Controls, both mandatory.** (i) A **matched-magnitude isotropic thermostat** — the same
displacement distribution with no objective information — which prices "the sample is just
smoothing". (ii) A **constant α-helix start** for the sampler, since §8 records that two properly
powered zero-information controls in Sprint 18 matched their informative arms.

**Falsifier.** If the posterior mean at the calibrated temperature does not beat the argmin by
≥0.084 Å at n=126 with a fold-aware CI excluding zero, **or** if the isotropic thermostat matches
it, close it. Note the honest hazard in advance: sampling then averaging **contracts** the
backbone, and the programme has measured that averaging contracts 22–25% already — so an Rg control
must be reported, not discovered afterwards.

**And the honest limit of R1 under T3, stated so nobody is surprised by it.** If the error is
coherent and in-manifold, the posterior is *concentrated at the wrong place*, and sampling will not
rescue it — a broad posterior helps only against the orthogonal component, which the argmin already
suppresses. So R1's realistic ceiling is small, and its value is diagnostic: it is the one lever in
the "do not consume the argmin" family that nobody in this programme is working on, it is cheap,
and a null on it would close that family cleanly rather than leaving it as a standing "someone
should try sampling".

## R2 — REPLACE THE LEAST-SQUARES INVERSION WITH A LEARNED DECODER

**The argument.** T2 says no per-pair loss can express joint conformer selection. The literature
says the same thing empirically: distance-restraint minimisation is the *pre-AF2* generation's
known failure mode, and AF2/AF3 kept the distogram as an auxiliary head while replacing the
inversion with a learned decoder over the pair representation. This programme has independently
measured that the inversion is where the loss is — *"where you average beats what you rank with"*,
coordinate averaging beating torsion averaging by 1.024 Å while the objective contributes 0.171 Å.
**The stage that is never questioned is the one the field replaced.**

**First experiment.** Train a small decoder — input the full 17-bin distogram plus separation and
sequence features, output (φ, ψ) per residue — on a **coordinate/FAPE-style loss**, leave-fold-out
on the pinned 5 folds, and emit its structure directly.

**Controls.** (i) The same decoder fed a **shuffled** distogram — the zero-information control that
says whether the gain is the distogram's information or the decoder's learned Ramachandran prior;
this matters because §8 records a constant α-helix matching conditioned torsion channels twice.
(ii) A **constant ideal α-helix** output. (iii) The **coordinate average** at 3.048 Å as the number
to beat.

**Falsifier.** If the decoder does not beat 3.048 Å at n=126 with a fold-aware CI excluding zero,
**or** does not beat its shuffled-distogram control, close it.

**What is new here and what is not, stated honestly.** The programme has already measured that a
*sequence-only* torsion head emits 4.151 Å and that φ carries no sequence signal at peptide length.
This proposal is not that. Its content is **distogram conditioning plus a coordinate-space loss** —
the decoder is asked to invert a distance distribution, not to predict torsions from sequence. If
the shuffled-distogram control matches it, the two are the same experiment after all and it closes
immediately.

## Lever (2) — "do not consume the argmin at all": what it would actually mean

Nobody in the programme is on this, and it is worth one paragraph of definition because it is
larger than R1's Langevin sampler, which is only its cheapest instance.

Every stage of this pipeline currently ends in an **argmin or an argmax over a scalar**: the fit
returns `res.x`; selection returns a top-*m*; the operator returns one coordinate average. T3 says
the fit's argmin is exactly where the in-manifold error gets realised, and Sprint 18 says the
argmin's location is worse than the structure the pipeline already builds. **The alternative is to
make the pipeline's output a distribution and to defer collapsing it until the last possible
moment — or to never collapse it at all and report an ensemble.** Concretely, three forms, in
increasing order of how much they change: (i) *keep the fit but consume its posterior* — R1;
(ii) *keep the ensemble downstream of the fit*, so the deliverable is a set of structures with
weights and the metric is an ensemble metric, which is how NMR peptide structures are deposited in
the first place and how this programme's targets are actually defined; (iii) *never form a point
estimate anywhere* — retrieval emits a weighted set, the objective reweights rather than optimises,
and the terminal operator is a weighted barycentre whose weight vector is the only thing ever
chosen.

**The reason it is worth defining rather than dismissing.** The frozen metric is a point metric —
full-chain Cα-RMSD of one structure — so form (ii) cannot be reported as the primary number without
changing the metric, which §7 forbids. That is a real constraint and it is why nobody has done it.
But the programme has already measured the two facts that make the direction interesting: **the
terminal operator consumes the set MEAN, not the set BEST** (`d_out = 1.16·d_set_mean +
0.04·d_set_best`), and **"do not collapse the posterior" was worth −2.253 Å** in the torsion-
restraint work. A pipeline whose every stage optimises a point estimate, feeding an operator that
only reads the set mean, is mis-matched at every joint. Making that mis-match explicit — measuring
what the pipeline would produce if no stage ever argmin'd — is a legitimate diagnostic even under a
point metric, because the ensemble's own mean is a point structure and can be scored.

## What I am NOT recommending, and why

Equivariant networks, pair-formers, mixture-of-experts routing over target difficulty, and
test-time optimisation are all plausible and all fail the same test: on 126 targets of 9–16
residues, with the sequence channel already measured to be the ceiling, they add capacity to a
problem that is **not capacity-limited** — Sprint 14 measured in-band ordering saturating at a
*linear* model, and P6 above measured a boosted tree extracting nothing a ridge did not. Adding
architecture to a signal-limited problem is the failure mode this lane exists to prevent.

---

# 7. LEDGER ENTRIES

See `s19/PREREG_D.md` for the pre-registrations. Statuses:

| block | hypothesis | result | status |
|---|---|---|---|
| P1 | the census is a bin-edge artefact | density correction changes nothing (0.240 vs 0.241); prominence and re-binning cut it to 0.097/0.155 | **F-D1 did NOT fire as pre-registered**; census damaged on other grounds |
| P2a | multimodality is harmful at the aim point | +0.410 unmatched, +0.007 [−0.291, +0.316] NN-matched | **REFUTED** (the harm) |
| P2b | the modes carry no information beyond their offsets | −0.484 [−0.570, −0.401] above the min-of-N null | **F-D2b did NOT fire** — the information is real |
| P3 | the multimodality is coherent across pairs | adjacency z=+0.27; signed coherence +0.003 [−0.001, +0.008] | **REFUTED** |
| P3′ | (unregistered, discovered) the real residual's coherence is carried by multimodality | equal on unimodal (+0.0429) and multimodal (+0.0489) pairs | **REFUTED** — the sprint's lead |
| P4 | the modes are an ensembling artefact | not run; P2a/P3 make its outcome unable to change a conclusion | **OPEN — NOT MEASURED** |
| P6 | some native-free rule extracts the mode information | every out-of-sample arm spans zero, including split-half within fold | **CLOSED by mechanism** |
| T1 | the selection path can express multimodality | the Bayes risk is convex in x | **EXACT — refuted by identity** |
| P5 | mode targets beat moment targets structurally | native-free mode +0.035, median +0.018; ORACLE mode-picker 3.472 vs the 3.048 the pipeline builds; mode information −0.053 [−0.153, +0.044] against a min-of-N-matched null | **F-D5 did NOT fire; F-D5b splits by control** — reported both ways; branch **CLOSED on the mission** |
| T3 | the in-manifold mechanism, and κ's null | rank = 2n−5 EXACT; geometric null 0.330 measured vs 0.329 predicted; real 0.695 vs matched sign-randomised 0.397, +0.298 [+0.262,+0.331], 121/5; per-target excess vs gap +0.134 [−0.057,+0.301] | **decision rule FIRED** — reference not model; design consequence held at **SUPPORTED** |
| T3′ | (unregistered, discovered) post-fit κ predicts the per-target gap where start-point geometry does not (+0.360 / +0.482 vs +0.157) | not explained; κ and the gap share the same two fits, so a common cause is not excluded | **OPEN** |

---

# 8. A METHODOLOGICAL FINDING, stated as a finding and not as an apology

This lane made **five corrections to its own output** in one working day. Four of the five were
caught by **running a measurement**, not by re-reading an argument. That ratio is the finding.

| # | the error | how it was caught |
|---|---|---|
| 1 | A detached `nohup … &` run was killed when its parent shell exited, and its relaunch then contended for the same JSON — two writers, one file. | The `complete: false` flag held. No number from it was ever reported. Fixed by a per-run result directory (`D_P5_TAG`). |
| 2 | "Conditioned on spread and separation, a multimodal mean is **at least as close** to the truth" — leaning on the most model-dependent of three controls, and it was my own significant result. | A **fourth** control: non-parametric nearest-neighbour matching, +0.007 [−0.291, +0.316]. Withdrawn from the coordinator's ledger at my own request. |
| 3 | The tangent-space null quoted as `(2n−2)/P = 0.381`, taken from Sprint 18's B5. | Checking injectivity of (φ,ψ) → CA trace: B5 constrains which torsions matter, not the rank of the distance map. Corrected to 2n−5. |
| 4 | On the strength of (3), the claim that Agent A's incoherent κ values "sat on the null" — comparing a **post-fit** quantity to a **linear** null. | Measuring the linear fractions directly. |
| 5 | `frac_in` projected with `np.linalg.qr(J_w)`, whose reduced form returns 2n orthonormal columns from a rank-deficient matrix, so every fraction was computed against a 2n-dimensional subspace. | **The empirical null refusing to match the algebra**: it came out 0.417 ≈ 2n/npairs instead of 0.329. |

**The sequence in (3)–(5) is the instructive part, and it is not what it looks like.** The algebra
was right the whole way through — rank = 2n − 5, verified at 20.92 against a predicted 20.92, and
the whitened null measured at 0.330 against a predicted 0.329. What was broken was the
*implementation of the test*. Correct reasoning, broken instrument, and **only running it twice
separated those two.** A review of the derivation would have passed it; a review of the code might
have caught the QR; only the disagreeing null caught it immediately and unambiguously.

> **The standing rule this lane would put in the programme's record: when an analytic null and a
> measured null disagree, the measurement is the null.** Publish the measured one, keep the
> derivation as a check on it, and treat the disagreement as a defect in the instrument until
> proven otherwise — because that is what it was, both times, here.

The corollary for §8 of the brief, which already says *"controls are mandatory, and uniform is not a
control"*: **a control must be matched in the space the operator actually works in.** `iso` and
`shuffled` are isotropic in raw distance space and the fit works in whitened space, so their
correct null is 0.581, not 0.330 — a 0.25 discrepancy that looked like a finding until it was
measured. The same defect in a different guise cost Sprint 18 a retracted claim (`objceil.py`
line 163 permuting the weights along with the residuals) and cost this sprint the `modeshuf`
inference (permuting mode offsets across pairs also orphans them from their weights). **Three
instances in two sprints of a control that was matched in the wrong space.** It is the programme's
most repeated error class and it deserves its own line in the next brief.

---

# 9. WHAT I LEAVE OPEN, AND WHAT WOULD CLOSE IT

| # | open question | what would close it |
|---|---|---|
| O1 | **Why does post-fit κ predict the per-target gap when start-point geometry does not?** ρ(κ, gap) = +0.360 and ρ(κ_real − κ_signflip, gap) = +0.482, against +0.134 [−0.057, +0.301] for the linear in-tangent excess — while the linear projection tracks κ itself well per target (ρ ≈ +0.66). **I have no explanation.** Something along the ~5-radian trajectory carries per-target information that the geometry at the start does not. | A **native-free predictor of κ**, which does not currently exist. Until then the +0.482 is confounded: κ and the gap are both derived from the same two fits, so targets where the fit converges cleanly having both a low objective and a low RMSD is not excluded as a common cause. A second route: measure the in-tangent fraction along the trajectory, not only at its start, and see where the per-target signal appears. |
| O2 | **Block P4 — are the modes an ensembling artefact?** Not run; P2a and P3 made its outcome unable to change any conclusion in this sprint. | Per-member `predict_proba` for one fold against the ensemble mean, on the same pairs. Cheap. Only becomes live again if some future arm finds a use for the mode structure. |
| O3 | **Is the in-manifold component shared across predictor families?** This gates lever (1) entirely: if PairNet, a second MLP seed and the retrieval pool's own error carry the *same* in-manifold component, it is an identifiability limit of the inputs and no loss in any basis removes it. | Agent A's cross-family attribution. **Not my lane, and I have not touched it.** |
