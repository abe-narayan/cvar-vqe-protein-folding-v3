# SPRINT 16 — AUDIT workstream findings

Owner: AUDIT (runs first). *Are the two Sprint 15 measurements the flagship rests on as strong as
recorded?* Tiering: **DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED**.

**Every native-derived quantity below is labelled ORACLE, in the tables and in the prose.**
The statistical unit is the **TARGET** (n = 126). Intervals are `I.paired`'s paired bootstrap over
targets with the per-fold means printed; a fold-clustered twin is printed wherever it changes a
verdict.

**Reproduction**

    python -m s12.instrument                 # the five pinned constants
    python -m s16.audit_align                # Parts A/B/C/D, 126 targets   (~5 min)
    python -m s16.audit_align2               # the ORACLE line search        (~6 min)
    python -m s16.audit_align3               # its missing random-direction null (~7 min)

Artefacts: `s16/results/audit_align{,2,3}.json` and `.log`. Threads pinned to 2 throughout; peak
RSS under 200 MB; one heavy job at a time, launched only below 90% CPU.

---

## 0. VERDICT TABLE

| # | claim under audit | recorded | measured here | verdict |
|---|---|---|---|---|
| **1a** | \|cos\|(θ_fit − θ_pool, true error) — **ORACLE** eval | **0.390** | **0.390** (median 0.379) | **RECONFIRMED**, exactly |
| **1b** | its null | **0.157** | isotropic null is right (**0.158** analytic, 0.160 sampled) but the **magnitude-matched** null is **0.216** | **WEAKENED** |
| **1c** | is \|cos\| computed on wrapped or raw differences? | (not stated) | **wrapped on both sides**; the unwrapped twin would read **0.676** | **RECONFIRMED** — the conservative choice was made |
| **1d** | alignment tracking **+0.499 [+0.33, +0.65]** — **ORACLE** eval | +0.499 | **+0.499**, i.i.d. CI **[+0.326, +0.645]**, fold-clustered **[+0.407, +0.580]**, per-fold +0.61/+0.64/+0.45/+0.32/+0.48 | **RECONFIRMED** |
| **1e** | is it confounded with error MAGNITUDE? | not tested | corr(\|cos\|, ‖e‖) = **−0.327** (wrong sign for that confound); partial \| ‖e‖ = **+0.367** | **RECONFIRMED** — not a magnitude artefact, though 0.216 of the 0.390 is per-coordinate magnitude |
| **1f** | *"a native-free error-direction surrogate exists"* — the **pool channel's own** contribution | implied to be the whole 0.390 | a **zero-information constant α-helix reference** gives \|cos\| **0.357** and tracking **+0.441**; the pool adds **+0.033 [−0.007, +0.075]** | **WEAKENED severely** |
| **1g** | operational value as a steering direction | untested (§5 open lead) | ORACLE line search: fit 3.677 → **random direction 2.894** → **α-helix 2.679** → **pool surrogate 2.584** → true error 0.317. The pool channel is **2.8%** of the ceiling | **WEAKENED to near-nil** |
| **2a** | bottom-half subspace overlap cos² **0.827** | 0.827 | **0.827** | **RECONFIRMED**, exactly |
| **2b** | its null **0.499** | 0.499 | correct, and it is exactly **k/m = 0.500** — trivially "half of the space" | **WEAKENED** (right number, near-vacuous baseline) |
| **2c** | the shrinking-subspace ladder | half 0.827, quarter 0.735 | half **0.827**, quarter **0.735**, eighth **0.903**, exactly-null-5 **0.830** — but ≥ 0.800 of the null-5 figure is **automatic** (4 of the 5 are the same fixed coordinate axes) | **RECONFIRMED numerically, REINTERPRETED** |
| **2d** | *"the geometry does not need the answer"* (native-free identifiability) | 0.827 vs 0.499 | **zero-information controls match or beat it**: constant α-helix **0.829**, a random pool window **0.833**, spectra +0.956/+0.978. Paired: indistinguishable at the half, **significantly worse** at the quarter (−0.016 [−0.029, −0.003] and −0.029 [−0.043, −0.016]) | **REFUTED as an identification claim** |
| **2e** | **five** exactly-null directions | 5 on 126/126 | **5 on 126/126**; the 4 coordinate torsions project into the null space at **1.000000** (min over 126×4) | **RECONFIRMED** |
| **2f** | the fifth is a distributed whole-chain crank, Σδφ ≈ −Σδψ | qualitative | deflating the 4 axes leaves **exactly one** residual direction (residual SVs 1.000, 0.000, 0.000); damage ‖J f₅‖ = **5e-15** against s₁ = 19.1; participation **8.41** coordinates, max load 0.503; −Σδψ/Σδφ **median 0.866** (mean 0.522, p10–p90 −0.72 to +1.97) | **RECONFIRMED, with the crank ratio looser than "≈"** |
| **3a** | ORACLE rotation **1.855 Å**, −1.822 [−2.037, −1.613], W/L 121/5 | as recorded | **1.855**, **−1.822 [−2.037, −1.613]**, W/L **121/5** | **RECONFIRMED**, exactly |
| **3b** | what the operation *does* | "rotates the fit's own error at fixed magnitude" | **(b)** — it builds **θ_native + e_rotated** (`align_fit.py:233–242`). Algebraically equal to θ_fit − P Pᵀe for the un-normalised twin, so (a) and (b) coincide there; both need e, hence θ_native | **CONFIRMED as (b)**; labelling was already correct, the **quotation** must change |
| **3c** | how much of the 1.822 Å is the error's **direction** | implied: all of it | its own zero-information control — ‖e‖ × a **random** direction in the **same** quiet half — emits **2.164 Å**. Direction is worth **−0.308 [−0.435, −0.172]**, W/L 97/29 = **16.9%**; **83.1%** is magnitude + subspace + the native anchor | **WEAKENED** |
| **3d** | the arm is norm-preserving only up to wrapping (−6.84°) | recorded as a defect of the arm | the **structure** is built from a perturbation of **exactly** ‖e‖; only the *reported* `tors_rms_deg` is affected by re-wrapping, and torsions are periodic so RMSD is untouched | **CORRECTED in the arm's favour** |
| **4** | quiet motion barely changes the objective **and** barely changes RMSD | untested theory | **RMSD: yes** (quiet/loud = 0.0018 at 5°, 0.0024 at 20°). **Objective: no** (0.0180 at 5°, 0.0306 at 20°; a 20° quiet step raises the objective by **53% of its own value at the optimum**). **Ratio of ratios 10–17**, per-target median 11.4, **> 1 on 99.2% of targets** | **HALF-REFUTED** — the first clause fails, the second holds |
| **5** | leakage and reproducibility | — | 60-target benchmark **0/60** in the tuning set and read by **no** module in `s15/` or `s16/` except `audit_provenance.py`'s byte hash; `s15/seed.py` used at **every** multi-start draw; the five pinned constants reproduce **exactly**. **One rule violation found**: `s15/align_jac.py:136` seeds the subspace null with `hash()` | **RECONFIRMED, with one §3.7 violation** |

---

## 1. THE ERROR-DIRECTION SURROGATE

`python -m s16.audit_align` → `s16/results/audit_align.json`, Part A. n = 126.

### 1.1 The recorded numbers reproduce exactly — DEMONSTRATED

`s15/align_free.py` was re-implemented from its inputs (same starts under `SD.stable_rng(p)`, same
leave-fold-out separation-debiased restraints, same six-start objective-only selection) and every
headline reproduces to three decimals: \|cos\| for `θ_fit − θ_pool` **0.390**, for `θ_fit − θ_start`
**0.354**, the isotropic null **0.158**, the alignment-tracking correlation **+0.499** with i.i.d.
CI **[+0.326, +0.645]** against the recorded [+0.33, +0.65]. Nothing in §5 of
`s15/align_FINDINGS.md` is mis-transcribed.

### 1.2 The wrapping is right, and it is the conservative choice — DEMONSTRATED

`align_free.py:86` builds the true error as `wrap(phi − nphi) ⊕ wrap(psi − npsi)` and `:99–102`
wraps every surrogate. Both sides are wrapped. Removing the wrap entirely gives \|cos\| = **0.676**,
and wrapping only the surrogate gives **0.425**. So angle wrapping *deflates* this statistic by 0.29
and the recorded 0.390 is the conservative reading, not an inflated one. **The audit's concern here
is disposed of in the recorded work's favour.**

### 1.3 The null is right for the model it states, and wrong for the model it needs — WEAKENED

0.157 is `sqrt(2/(π·2n))` at 2n = 25.9, the expected \|cos\| between two **isotropic** directions.
That is arithmetically correct (sampled: 0.160). But the surrogate is not isotropic: its
per-coordinate magnitude profile correlates with the true error's at **+0.447** pooled over
coordinates. A **matched** null that keeps that profile and randomises only the signs
(`v_k → ±|v_k|`) reads **0.216**; a within-block permutation null reads **0.166**.

| null | value | 0.390 − null (paired) | W/L |
|---|---:|---|---|
| isotropic, as recorded | 0.158 | +0.2304 | 100/26 |
| block permutation, MATCHED dimensionality | 0.166 | +0.2241 | 98/28 |
| **sign-flip, MATCHED per-coordinate magnitude** | **0.216** | **+0.1745** [+0.147, +0.216] | **93/33** |

**The honest headline is 0.390 against 0.216, not against 0.157.** A quarter of the recorded gap is
"the surrogate is big in the same coordinates the error is big in", which carries no direction.

### 1.4 It is NOT confounded with error MAGNITUDE — RECONFIRMED

The brief's suspicion was that the surrogate is large exactly where the error is large. The sign is
the opposite: corr(\|cos\|, ‖e‖ **ORACLE**) = **−0.327**, corr(\|cos\|, torsion RMS) = **−0.363**,
corr(\|cos\|, RMSD) = **−0.170**. The surrogate is *better aligned on the targets where the error is
small*. Partialling ‖e‖ out of the tracking correlation takes +0.499 to **+0.367**; partialling out
both ‖e‖ and torsion RMS leaves **+0.350**. **The tracking result survives the magnitude control.**

### 1.5 It survives the target as the unit and a fold-aware interval — RECONFIRMED

n = 126 targets was already the unit. Adding a fold-clustered bootstrap (resampling whole folds)
gives **[+0.407, +0.580]**, and the five per-fold correlations are **+0.606, +0.643, +0.452, +0.315,
+0.479** — all positive, none crossing zero in sign. **5/5 folds the same sign.**

### 1.6 …but almost none of it comes from the POOL — WEAKENED SEVERELY

This is the finding that matters for the sprint. The surrogate `θ_fit − θ_pool` is a **difference
between the fit and a reference**. The project's standing trap (*"a zero-information constant
α-helix beats uniform random"*) says the reference must be controlled. It never was.

| reference direction subtracted from θ_fit | \|cos\| with the true error (**ORACLE**) | tracks the true alignment (**ORACLE**) |
|---|---:|---|
| **the retrieval pool circular mean** (recorded) | **0.390** | **+0.499 [+0.326, +0.645]** |
| the winning multi-start | 0.354 | +0.346 [+0.168, +0.511] |
| **a constant α-helix — ZERO INFORMATION** | **0.357** | **+0.441 [+0.261, +0.606]** |
| **one random pool window — ZERO INFORMATION** | **0.340** | **+0.360 [+0.161, +0.547]** |
| matched sign-flip null | 0.216 | — |

**pool − α-helix on \|cos\| = +0.033 [−0.007, +0.075]**, W/L 70/56, median +0.010, concentration
PASS (p = 0.665, mean/sd +0.14, *low power flagged*). The interval **includes zero** by the
project's standard paired bootstrap. (A fold-clustered twin gives [+0.009, +0.054]; with only five
clusters that interval is not the more trustworthy of the two, and the honest reading is the
i.i.d. one.)

**So the zero-information reference reproduces 0.357 of the 0.390 — 91% of the raw statistic and
81% of its excess over the matched null (0.141 of 0.174).** What `θ_fit − θ_pool` mostly measures is *the fit's own displacement from any smooth
chain-like torsion field*, which is dominated by θ_fit — and θ_fit appears in the true error too.
The "channel disagreement" reading is not supported.

*A structural detail that reinforces this.* `θ_pool` **is** `starts[0]`: `s15/distgeo.starts`
makes the retrieval top-75 circular mean the first multi-start initialisation, and `align_free`
builds its surrogate from the same circular mean. On the **24 of 126** targets where the
objective-argmin start is `retrieval_circmean` (`s16/results/audit_align2.json`), rows 1 and 2 of
`s15/align_FINDINGS.md` §5 are **literally the same vector** and are not two surrogates.

### 1.7 What it is worth as a STEERING direction — the number the flagship needs

`python -m s16.audit_align2` and `audit_align3`. Move the fit along −c·‖e‖·û for the best c on a
61-point grid spanning ±1.5, with **c chosen on the true RMSD** — an **ORACLE** line search that
gives away both the step size and the sign, i.e. strictly more than any native-free steerer can have.

| direction moved along (all with an **ORACLE** line search) | RMSD | step |
|---|---:|---:|
| the `squared` fit itself | 3.677 | — |
| **a RANDOM direction** — the line search's own null | **2.894** | **−0.783** |
| the ZERO-INFORMATION constant-α-helix direction | 2.679 | −0.215 |
| **`θ_fit − θ_pool`, the recorded surrogate** | **2.584** | **−0.095** |
| the TRUE error direction (**ORACLE** ceiling of a 1-D move) | **0.317** | −2.267 |

* **pool − random**: −0.310 [−0.457, −0.171], W/L 78/48, concentration PASS (p = 0.875, mean/sd −0.38).
* **α-helix − random**: −0.215 [−0.370, −0.072], W/L 66/60, PASS (low power).
* **pool − α-helix**: **−0.095 [−0.219, +0.026]**, W/L 68/57 — **not significant**.
* The ORACLE step size along the surrogate has **61% sign consistency** (mean c = +0.29, median
  +0.45). A native-free steerer must supply that sign; `error-coherence-decides-correctors` records
  that a 0.688-accurate sign with **coherent** mistakes emits **+0.31 Å**.

**The pool channel is worth 0.095 Å of a 3.360 Å ORACLE ceiling — 2.8% — with the step and the sign
already given away.** An ORACLE line search along a *random* direction is worth eight times more.

### 1.8 One thing that goes the flagship's way — DEMONSTRATED

The component a steerer must actually cancel is the error's **loud**-half component (Part D shows
quiet-half error is free). Restricted to the loud half of J at the fit:

| | value |
|---|---:|
| share of ‖e‖² in the loud half | 0.440 |
| **\|cos\| of the surrogate with the error, INSIDE the loud half** | **0.497** (median 0.508) |
| its matched sign-flip null | 0.304 |
| its isotropic null | 0.224 |
| paired gap over the **matched** null | **+0.193**, W/L 101/25 |

**The surrogate is better aligned with the RMSD-relevant half of the error (0.497) than with the
error as a whole (0.390).** Removing the exact null space alone takes it to 0.480, and restricting
to the 2n−3 coordinates the builder reads takes it to 0.422. *(I predicted the opposite for that last
one — that the three never-read torsions, where θ_fit is pinned to θ_start on both sides, would
manufacture agreement. They dilute it instead. The prediction is wrong and is preserved here.)*

**Tier: DEMONSTRATED for the existence of the signal; REFUTED for the claim that it is a
pool-channel signal; WEAKENED to ~3% of the ceiling as a steering direction.**

---

## 2. THE QUIET-SUBSPACE ESTIMATE

`python -m s16.audit_align` Part B, n = 126. J at the **native** torsions is **ORACLE**; every other
J is native-free.

### 2.1 0.827 reproduces, and 0.499 is exactly half the space — RECONFIRMED / WEAKENED

The recorded bottom-half overlap **0.827**, bottom-quarter **0.735** and spectral correlation
**+0.980** all reproduce. For two random k-dimensional subspaces of an m-dimensional space the
expected mean cos² of the principal angles is **k/m**; at k = m/2 that is **0.500**, and the sampled
null in `align_jac.py` returns 0.499. **The null is correct and it is near-vacuous**: a bottom-*half*
subspace starts at 0.5 by arithmetic.

### 2.2 The zero-information controls match it — REFUTED as an identification claim

The claim under audit is `s15/align_FINDINGS.md` §2.5: *"The geometry does not need the answer."*
Test it against Jacobians that contain **no information about the target at all**.

| J compared against **J(native)** — **ORACLE** eval | bottom half | bottom quarter | bottom eighth | the null 5 | spectra |
|---|---:|---:|---:|---:|---:|
| **incumbent emitted torsions (as recorded)** | **0.827** | **0.735** | 0.903 | 0.830 | **+0.980** |
| the restraint fit's own structure | 0.819 | 0.724 | 0.901 | 0.821 | +0.982 |
| **a constant α-helix — ZERO INFORMATION** | **0.829** | **0.751** | 0.917 | 0.849 | +0.956 |
| **one random pool window — ZERO INFORMATION** | **0.833** | **0.765** | 0.919 | 0.857 | +0.978 |
| analytic random-subspace null `k/m` | 0.500 | 0.249 | 0.124 | 0.199 |  |

Paired, n = 126 (positive = the emitted structure's subspace is the better estimate):

| paired comparison | bottom half | bottom quarter |
|---|---|---|
| incumbent emitted − constant α-helix | −0.0022 [−0.0105, +0.0066] W/L 60/66 | **−0.0158 [−0.0287, −0.0029]** W/L 55/71 |
| incumbent emitted − random pool window | −0.0057 [−0.0132, +0.0016] W/L 56/70 | **−0.0292 [−0.0431, −0.0160]** W/L 44/82 |

**At the half the emitted structure is statistically indistinguishable from a zero-information one;
at the quarter it is significantly WORSE than both.** The overlap is a property of *any* chain
of that length with ideal virtual geometry — it is not an estimate of anything target-specific. The
recorded number is right; the inference drawn from it is not.

Deflating the four known-inert coordinate axes out of both subspaces (they are shared exactly, by
construction, and inflate every row) does not change the verdict:

| deflated overlap | bottom half | bottom quarter |
|---|---:|---:|
| incumbent emitted | 0.744 | 0.257 |
| the fit | 0.730 | 0.240 |
| **α-helix ZERO-INFO** | **0.745** | **0.317** |
| **random pool window ZERO-INFO** | **0.751** | **0.338** |
| null `(k−4)/(m−4)` | 0.405 | 0.103 |

There is real geometric structure here — every row is well above its null — but **none of it is
information the emitted structure supplies.**

### 2.3 The ladder, and what "bottom-eighth 0.903" actually is

Bottom-eighth ≈ 3 directions at 2n ≈ 26, which sits **inside** the 5-dimensional exact null space.
The 0.903 is therefore not a strengthening of the result at higher resolution; it is the null space
agreeing with itself. Likewise the null-5 figure: four of those five directions are the *same fixed
coordinate axes* on every structure, so mean cos² ≥ 4/5 = **0.800 automatically**. Measured 0.830
implies cos² of the **fifth** direction, native vs emitted, of **0.150** — i.e. the one
target-specific null direction **does not transfer between structures**. At the resolution where the
subspace becomes specific, the overlap collapses.

### 2.4 The five exactly-null directions — RECONFIRMED independently

* **5 exactly-null directions on 126/126 targets** (threshold `s ≤ 1e-8·s₁`; the gap is 1e-15 vs 2e-3).
* The four coordinate torsions `phi[0], phi[n−1], psi[0], psi[n−1]` project into the null space with
  fraction **1.000000** — the minimum over all 126 × 4 measurements.
* Deflating those four out of the 5-dimensional null basis leaves residual singular values
  **1.000, 0.000, 0.000, 0.000, 0.000** — **exactly one** further direction, on every target.
* That fifth direction has damage ‖J f₅‖ = **5.0e-15** against s₁ = 19.1: exactly inert.
* It is **distributed**: participation ratio **8.41 coordinates** out of ~26, largest single loading
  **0.503**. It is not a coordinate.
* The crank signature: **−Σδψ / Σδφ has median 0.866**, mean 0.522, p10–p90 −0.72 to +1.97, and lies
  in [0.5, 1.5] on **57%** of targets. **"Σδφ ≈ −Σδψ" is right on the median target and loose in the
  tails** — it is a whole-chain crank, but the compensation is not tight enough to state as ≈ without
  the spread beside it.

**Tier: DEMONSTRATED for the count, the derivation and the distributed character; RECONFIRMED with a
looser crank ratio than recorded; REFUTED for the native-free *identifiability* inference.**

---

## 3. THE ORACLE CEILING — WHAT THE ROTATION ACTUALLY DOES

### 3.1 The number reproduces exactly — RECONFIRMED

`s15/results/align_fit.json` and my independent recomputation agree to three decimals:
`squared` **3.677**, `ORACLE_kill_loud_0.5` **1.443**, `ORACLE_kill_loud_norm_0.5` **1.855**,
**−1.822 [−2.037, −1.613]**, W/L **121/5**, per-fold −1.82/−2.03/−1.88/−1.96/−1.50.

### 3.2 It is (b): it reconstructs from θ_native + rotated error — CONFIRMED

`s15/align_fit.py:228–242`, verbatim:

```python
e  = np.concatenate([A.wrap(sq_phi - nphi), A.wrap(sq_psi - npsi)])   # ORACLE error
P  = V_sv[:, :r]                                  # loud half of J at the FIT (native-free J)
e2 = e - P @ (P.T @ e)                            # the error's quiet-half component
th2 = np.concatenate([nphi, npsi]) + e2 * (np.linalg.norm(e) / np.linalg.norm(e2))
```

The emitted structure is **θ_native + a rescaled quiet-half component of the error**. Two honest
qualifications in both directions:

* Because θ_fit = θ_native + e (up to wrapping), the un-normalised twin is algebraically identical
  to **θ_fit − P Pᵀe**, so for `ORACLE_kill_loud_0.5` readings (a) and (b) coincide. The normalised
  arm equals θ_fit − P Pᵀe + (κ−1)e₂ with κ = ‖e‖/‖e₂‖, which is still expressible from the fit but
  needs ‖e‖.
* Either way it **needs e, hence θ_native**. This is not a leak — the arm is named `ORACLE_*` and
  `align_FINDINGS.md` labels it throughout. **What must change is the quotation, not the labelling.**

Note also that the operation is **not** "steer into the quiet subspace". −P Pᵀe is a displacement
**inside the loud subspace**, of RMS magnitude √0.440 × 78.65° = **52.2°**. The arm *deletes loud
error*; it does not move anything into the quiet half.

### 3.3 …and 83% of it has nothing to do with the error's DIRECTION — WEAKENED

The arm's own zero-information control: build **θ_native + ‖e‖ × û** where û is a **random** unit
vector in the *same* bottom-half right-singular subspace of J(fit). Same magnitude, same subspace,
same native anchor, **zero** information about which way the fit is wrong.

| arm — every one built as **θ_native + δ**, all **ORACLE** | RMSD | median | vs the `squared` fit |
|---|---:|---:|---|
| the `squared` fit itself | 3.677 | 3.595 | — |
| `ORACLE_kill_loud_0.5` (confounded twin) | 1.443 | 1.391 | −2.234 [−2.462, −2.018] W/L 124/2 |
| **`ORACLE_kill_loud_norm_0.5` — the recorded ceiling** | **1.855** | 1.745 | **−1.822 [−2.037, −1.613]** W/L 121/5 |
| **CONTROL: ‖e‖ × a RANDOM direction in the SAME quiet half** | **2.164** | 2.126 | **−1.514 [−1.716, −1.325]** W/L 116/10 |
| CONTROL: ‖e‖ × a random direction in the **loud** half | 4.991 | 4.751 | +1.314 [+1.084, +1.535] W/L 19/107 |
| CONTROL: ‖e‖ × an isotropic random direction | 4.852 | 4.584 | +1.175 [+0.954, +1.385] W/L 22/104 |
| the same rotation but with the subspace taken from **J(native)** | 1.606 | 1.463 | −2.071 [−2.298, −1.856] W/L 126/0 |

> **`ORACLE_kill_loud_norm_0.5` minus its own zero-information control: −0.308 [−0.435, −0.172],
> W/L 97/29, median −0.413, per-fold −0.50/−0.40/−0.41/−0.16/−0.11, concentration PASS (p = 0.070,
> mean/sd −0.41).**

**16.9% of the recorded 1.822 Å is the error's direction. 83.1% is "a structure that sits ‖e‖ away
from the native entirely inside the quiet half of a chain Jacobian."** The correct quotation is:

> *A structure differing from the native by 72° RMS **placed anywhere in the quiet half** emits
> ≈2.16 Å; using the fit's actual error direction inside that half improves it to 1.855 Å.*

The concentration p = 0.070 is the closest to failing of any check in this file and should be
watched, not waved through.

### 3.4 A correction in the recorded work's favour

`s15/align_FINDINGS.md` §3 records a defect: *"`ORACLE_kill_loud_norm` preserves ‖e‖ before the
torsions are re-wrapped … the arm's measured torsion RMS falls by about 5° rather than by 0°."*
The **structure** is built from `θ_native + e₂·κ`, a perturbation of **exactly** ‖e‖, and
`I.build_ca` is periodic in the torsions, so the emitted geometry and its RMSD are unaffected. The
−6.84° is an artefact of the *reported statistic* in `score_axes`, which re-wraps before taking the
RMS. **The arm is exactly norm-preserving as built.** The recorded caveat is more damaging to itself
than the facts require.

---

## 4. THE THEORY CHECK — QUIET DIRECTIONS, THE OBJECTIVE, AND RMSD

**Hypothesis under test (stated in the audit brief):** quiet torsion directions barely change the Cα
trace; the restraint objective is a function of Cα distances only; therefore movement in quiet
directions barely changes the objective, the fit is nearly unconstrained along them, and RMSD is
nearly insensitive to them too.

`s16/audit_align.py` Part D, at the `squared` fit's own solution, with J and the quiet/loud split
taken from J at the fit (native-free). Objective `F = Σ_p w_p ((d_p − d̂_p)/sd_p)²`; the objective's
gradient is ~0 at the optimum, so the linear-regime statement is a **curvature** statement, and the
finite-step version is an exact rebuild with no linearisation.

| | quiet / loud |
|---|---:|
| Gauss-Newton curvature of the **restraint objective** | **0.0424** (median 0.0243) |
| curvature of **RMSD²** | **0.00251** (median 0.00258) |
| **RATIO OF THE TWO RATIOS** | **16.9** (per-target median **11.4**, p10 3.9, p90 40.4, **> 1 on 99.2% of targets**) |

Finite step, exact rebuild, objective at the fit = 82.27, RMSD at the fit = 3.677 Å:

| step | Δobjective quiet | Δobjective loud | ratio | ΔRMSD quiet | ΔRMSD loud | ratio | **ratio of ratios** |
|---|---:|---:|---:|---:|---:|---:|---:|
| 5° RMS | +2.20 | +122.3 | 0.0180 | **−0.0004 Å** | +0.085 Å | −0.0018 | **−10.0** |
| 20° RMS | **+43.28** | +1413.6 | 0.0306 | **−0.002 Å** | +0.758 Å | −0.0024 | **−12.6** |

**The second clause holds and the first fails.**

* **RMSD is effectively blind to quiet motion.** A 20° RMS displacement of the fit inside the quiet
  half changes Cα-RMSD by **−0.002 Å** — 0.05% of the fit's own 3.677 Å, and *negative*. Per unit
  motion RMSD is **~400×** less responsive to quiet than to loud directions.
* **The objective is NOT blind to it.** The same 20° quiet step raises the restraint objective by
  **43.3, i.e. 53% of its entire value at the optimum**. Relative to what RMSD cares about, the
  objective **over-prices** quiet directions by a factor of **10–17**, on 99.2% of targets.

**Three consequences, and they decide the sprint.**

1. **No native-free operation applied to the fit, inside the quiet subspace, can move RMSD.** The
   response is 0.002 Å per 20° RMS. "Rotate the error into the quiet subspace" cannot be implemented
   as a motion of the emitted structure, because motions there do not change the answer. This is not
   a control problem; it is a geometry problem.
2. **The ceiling's gain is a loud-subspace operation, not a quiet-subspace one** (§3.2): it deletes
   52.2° RMS of error from the directions the restraints *can* see. The reason the fit has that error
   is that the restraints are wrong — the fit's emitted Cα–Cα distances are **2.125 Å MAE** from
   the true ones (`squared.dist_mae`, **ORACLE** eval), which the rotation cuts to **1.340** — not
   that the fit points its error badly. **`the distance prior is the ceiling` is the operative law here, not
   a new alignment law.**
3. The fit's position along quiet directions is **not** unconstrained by the data: the objective
   resists quiet motion 17× more than RMSD warrants. The restraint fit therefore spends real
   optimisation effort pinning directions that do not affect the answer — which is a coherent and
   independent restatement of §4.6's *"a better restraint objective is not a better structure"*.

**Tier: DEMONSTRATED (both halves measured, with the theoretical clause REFUTED and its consequence
CONFIRMED).**

---

## 5. LEAKAGE AND REPRODUCIBILITY

### 5.1 The 60-target benchmark — SEALED, one previously-declared exception, nothing new

* `results/benchmark_manifest.json` (`n_targets: 60`, `peptide_db.benchmark(n=60, seed=11)`) is the
  protected set. Its 60 PDB IDs are **disjoint from the 126 tuning targets** (overlap = **0**).
* `grep` over every `.py` in `s15/` and `s16/` for `benchmark`, `benchN`, `bench60`,
  `benchmark_manifest`: the **only** hit outside comments is `s15/audit_provenance.py:38`, which
  passes the manifest path to `_hash()` and emits a digest, a byte count and two mtimes. That is the
  known and accepted exception, and `s15/verify_FINDINGS.md` B1 already logged it with a suggested
  fix. **Nothing else touches it.**
* The 17 `s16/*.py` modules present at the time of this audit (mine plus four other workstreams')
  contain **no** benchmark read.
* **Disclosed:** this audit itself opened the manifest to read the 60 **PDB IDs only**, to prove
  disjointness. No sequence, no coordinate, no label was loaded, and no ID appears in any artefact I
  wrote. Auditing is an explicitly permitted use under BRIEF §3.1.

### 5.2 Seeding — `s15/seed.py` at every multi-start, with one violation elsewhere

| site | seed | status |
|---|---|---|
| `s15/distgeo.py:142` | `SD.stable_rng(pdb, seed)` | OK |
| `s15/robust.py:233` | `SD.stable_rng(p)` | OK |
| `s15/align_fit.py:145` | `SD.stable_rng(p)` | OK |
| `s15/align_free.py:76` | `SD.stable_rng(p)` | OK |
| `s15/align_reg.py:137` | `SD.stable_rng(p)` | OK |
| **`s15/align_jac.py:136`** | **`np.random.default_rng(abs(hash(p)) % 7 + 5)`** | **VIOLATES BRIEF §3.7** |
| `s15/info_regime.py:59` | `np.random.default_rng(hash(p) % 2**31)` | violates §3.7; draws only the `"one random top-75 window"` channel, which `align_jac`'s `CH` list does not report |

`align_jac.py:136` seeds the **random-subspace null** — the source of the 0.499 — with `hash()`,
which Python salts per process, and then only into a **7-value** seed pool. Every multi-start draw in
the sprint is clean; this one null is not reproducible across interpreters. **The reported value is
nonetheless correct**: §2.1 shows analytically that the quantity is `k/m = 0.500`, and my own
`stable_rng`-seeded recomputation returns the same. No number changes; the rule was broken.

### 5.3 The instrument

`python -m s12.instrument`, threads pinned to 2:

    shipped 3.4540004952559396   pool_best 1.7108244199364904
    top75_best 2.3061526409453816   synthesis_fit 3.2040761603809194   n_zero_recall 18

**All five pinned constants reproduce exactly.** (A cache read, not a pipeline reproduction, per the
Phase 0 audit.)

---

## 6. IS THE FLAGSHIP HYPOTHESIS STILL WORTH THE SPRINT'S MAIN EFFORT?

**No — not in the form the brief states it, and the evidence says what should replace it.** Every
recorded number reproduces exactly; nothing here is a transcription problem. But all three legs of
the flagship's premise are weaker than they read. The **quiet-subspace estimate is not an estimate**:
a constant α-helix recovers the native's bottom-half subspace at cos² 0.829 against the emitted
structure's 0.827, so there is no target-specific information to combine, and at the resolution
where the subspace does become specific — the fifth null direction — the overlap collapses to 0.150.
The **error-direction surrogate is real but is not a pool signal**: a zero-information α-helix
reference reproduces 0.357 of its 0.390 \|cos\| and +0.441 of its +0.499 tracking, and under an
ORACLE line search that already gives away the step size *and* the sign it beats a random direction
by 0.310 Å while beating the α-helix by only 0.095 [−0.219, +0.026] — 2.8% of that ceiling. And the
**1.855 Å ceiling is 83% direction-free**: its own zero-information control, ‖e‖ × a random vector in
the same quiet half, emits 2.164 Å, leaving 0.308 [−0.435, −0.172] for the direction. Worse, the
theory check shows the two ingredients cannot be combined the way the hypothesis assumes: RMSD
responds to quiet-subspace motion at 0.002 Å per 20° RMS, so **no operation performed on the emitted
structure inside the quiet subspace can change the answer at all**, while the restraint objective
over-prices those same directions by 10–17×. The ceiling arm does not steer error into the quiet
half — it *deletes* 52.2° RMS of error from the **loud** half, an operation defined relative to
θ_native. What is left, and what I would put the sprint's main effort into, is the one measurement
that survived every control: **the surrogate's \|cos\| with the LOUD half of the error is 0.497
against a magnitude-matched null of 0.304 (W/L 101/25)** — the loud half is the only half that
changes RMSD, it carries 0.440 of ‖e‖², and cancelling it is exactly what the ceiling arm does. That
reframes the target from "steer into the quiet subspace" (geometrically inert) to **"cancel loud
error with a native-free estimate of it"** — which is a restraint-accuracy problem, is measured
against `the distance prior is the ceiling` rather than against a new alignment law, and must be run
with a **zero-information reference control and a random-direction control in every arm**, because
without those two the recorded evidence would have read as five times stronger than it is.

---

## 7. WHAT THIS AUDIT DID NOT SETTLE

* The concentration check on the ceiling's direction-only effect returns **PASS at p = 0.070** — the
  closest to failure in this file. It is not a fail, and it is not comfortable.
* `pool − α-helix` is measured on \|cos\| and on the ORACLE line search. It is **not** measured on
  an actual native-free steering pipeline, because none exists; the line search is an upper bound,
  not a simulation of one.
* The 0.497 loud-half \|cos\| is measured against a magnitude-matched null but **not** against the
  zero-information α-helix reference. That control is the first thing the flagship should run.
* The fifth null direction's transfer (cos² 0.150) is measured native-vs-incumbent. Whether it
  transfers between two *native-free* structures is untested and would matter for any arm that uses
  it.
* Everything here inherits `s15/align_FINDINGS.md` §9's scope limit: it is all measured inside the
  one continuous torsion-space restraint fit, which loses to the incumbent by +0.473 Å.
