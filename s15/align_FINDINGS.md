# SPRINT 15 — ALIGN workstream findings

Owner: ALIGN (Phase 2). *Can error alignment be ENGINEERED rather than inherited?*
Written continuously; newest sections appended. Tiering:
**DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

**Instrument selfcheck, start of workstream** (`python -m s12.instrument`) reproduces

    shipped 3.4540004952559396   pool_best 1.7108244199364904
    top75_best 2.3061526409453816   synthesis_fit 3.2040761603809194   n_zero_recall 18

(Per the Phase 0 audit this is a **cache read**, not a pipeline reproduction. Threads capped
at 2 everywhere, per the brief; the distogram is not bit-stable across thread counts.)

**Reproduction commands**

    python -m s15.align_verify     # the analytic Jacobian against info_null's central differences
    python -m s15.align_jac        # TASK 1: the low-response subspace, 126 targets   (~4 min)
    python -m s15.align_fit        # TASKS 2/3: the interventions, 126 targets        (~2 h)
    python -m s15.align_axes       # TASK 3: the two-axis honest test, reads align_fit.json
    python -m s15.align_free       # is there a NATIVE-FREE alignment diagnostic?     (~10 min)
    python -m s15.align_reg        # THE DECIDING TEST: same shrinkage, different metric (~50 min)

Shared machinery: `s15/align_lib.py`. All threads capped at 2; peak RSS under 100 MB per job;
one heavy job at a time on a machine that was at 98-100% CPU throughout from other agents.

---

## THE FOUR ANSWERS, ONE SENTENCE EACH

**(a) What does the low-response torsion subspace look like, and how much error can it
absorb?** It is almost the whole space — the CA trace of an ideal-geometry chain responds to a
**participation ratio of 3.47 effective directions out of 2n ≈ 25.9**, exactly **5** of them
are *exactly* null (the four known inert coordinate torsions plus a fifth, distributed,
target-specific combination nobody had recorded), and an angular error placed isotropically in
the quiet half of the spectrum tolerates **15.7× more RMS angle** than a random direction at
the same RMSD cost, so in principle the entire error budget can be absorbed there.
**[DEMONSTRATED, §2]**

**(b) Can alignment be induced native-free, and by which intervention?** **No — by none of
them.** Six intervention families and 22 arms on 126 targets give a best leave-fold-out result
of **−0.007 Å [−0.074, +0.055]**, no family's leave-fold-out arm has a CI excluding zero, every
arm moves the alignment statistic by at most 0.104 against the 0.387 the ORACLE ceiling
requires, and the one arm with a clear effect (`tik_LFO`, −0.216 Å) is **ordinary isotropic
shrinkage that beats both alignment-aware metrics at matched shrinkage**. The mechanism is real
and worth **1.822 Å [1.613, 2.037]** (ORACLE), and nothing native-free reaches it.
**[REFUTED, §4 and §7 — this is my own hypothesis]**

**(c) Are robustness and alignment the same lever or two?** **Neither — on the full instrument
they are two non-levers, and the interaction test has no power to separate them.** Robust
losses do not beat squared error on 126 targets (`CAUCHY_lfo` −0.011 [−0.102, +0.077]), which
fires `s15/robust.py`'s own falsification condition and refutes the outlier reading of K7; all
three combination arms test "additive", but with parts of |effect| ≤ 0.09 Å and interaction CIs
about 0.23 Å wide that verdict is uninformative and must not be quoted. What *is* measurable is
that their per-target gains correlate only **+0.24 to +0.45**, so they act partly on different
targets. **[REFUTED for robustness; the same/two question is UNDERPOWERED, §4.6, §4.7]**

**(d) After alignment engineering, is 2.5 Å supported and is 2.0 Å still refuted?** **2.5 Å is
still marginal and 2.0 Å is still not supported predictively — but the reason has changed from
information to control.** The ORACLE alignment ceiling emits **1.855 Å** from the *same*
restraints at the *same* error magnitude, so the information for sub-2 Å is demonstrably
present in the fit's own error vector; what is absent is any native-free mechanism that rotates
it, and the best measured is 4% of the required rotation. **[§6]**

---

## 1. THE JACOBIAN IS EXACT, AND IT INDEPENDENTLY REPRODUCES THE INFO HEADLINE

**DEMONSTRATED.** `python -m s15.align_verify`.

`s15/info_null.py` computes `J = d(superposed CA)/d(phi, psi)` by **central differences** at
1°, Kabsch-superposing each perturbed trace. `s15/align_lib.sup_jacobian` computes it
**analytically**, from the same rigid-suffix-rotation identity `core.project._torsion_grad`
uses in reverse mode, with the six rigid-body modes projected out algebraically. They share no
code. Over five targets:

| check | result |
|---|---|
| cosine(analytic, central-difference) | **0.999999998 – 0.9999999999** |
| worst relative disagreement | **1.4e-04** (= the O(h²) truncation of the 1° difference) |
| `pair_jacobian` vs central differences of the pair distances | cosine **1.000000000000000**, max 3e-08 |
| `‖J e‖/√n` vs the true Kabsch RMSD of the perturbation | ratio 0.999–1.030 over 0.2°–5° |

And the alignment statistic computed through the analytic Jacobian on all 126 targets
reproduces `s15/results/info_null.json` **to three decimals on every channel**: retrieval
top-75 circular mean 0.738, top-75 medoid 0.693, constant α-helix 0.709, incumbent projection
0.665, ORACLE best pool window **0.565**. Two independent implementations, one number. The
INFO workstream's mechanism claim is verified rather than re-asserted.

---

## 2. TASK 1 — THE LOW-RESPONSE SUBSPACE

`python -m s15.align_jac` → `s15/results/align_jac.json`, `align_jac.log`. 126 targets,
J at the **native** torsions. **ORACLE DIAGNOSTIC** for the spectrum-at-native; the
native-free twin is measured in §2.4 and it agrees.

### 2.1 The spectrum collapses almost immediately — DEMONSTRATED

Normalised singular values `s_i/s_1` against fractional index `i/2n`, mean over 126 targets:

| i/2n | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| s/s₁ | 1.000 | 0.456 | 0.261 | 0.137 | 0.060 | **0.040** | 0.032 | 0.021 | 0.007 | 0.000 | 0.000 |

| statistic | mean | median | p10–p90 |
|---|---:|---:|---|
| directions carrying 90% of `trace(JᵀJ)` | **4.51** | 4.00 | 3–6 |
| participation ratio `(Σsᵢ²)²/Σsᵢ⁴` (effective rank) | **3.47** | 3.39 | 2.16–4.88 |
| directions with `sᵢ ≤ 0.10 s₁` | 17.27 | 16.0 | 13–22.5 |
| `2n` | 25.9 | 26 | 18–32 |
| rank90 as a fraction of `2n` | **0.176** | 0.188 | — |

**The CA trace of a 9–16mer responds to about three and a half effective torsion directions
out of twenty-six.** That is the structural fact underneath the whole alignment story: the map
from torsion space to CA shape is extraordinarily ill-conditioned, so *where* an error points
matters far more than how big it is. This is measured here for the first time in the project.

### 2.2 There are FIVE exactly-null directions, not four — DEMONSTRATED, and it corrects K0

Sprint 14's correction, confirmed by the coordinator as K0, records **four** inert torsions per
chain (`phi[0]`, `psi[0]`, `phi[n-1]`, `psi[n-1]`). That is right and this workstream confirms
it — the four coordinate directions lie **entirely** inside the null space (projection
fraction 1.0000 each; column norms 0.0, 0.0, 4.9e-14, 0.0 against an interior median of 5.44,
a ratio of 9e-15).

But the null space is **five**-dimensional on **126 of 126 targets**, with a clean gap: the
smallest *non*-zero singular value is 2e-3 to 3e-2 of `s₁` while the five null values are
1e-15 of it. The fifth direction is **not** a coordinate direction; after projecting the four
out, a one-dimensional residual remains with singular value exactly 1.

**The count is derivable, which is the check that it is real and not a numerical artefact.**
A chain of `n` CA atoms with fixed virtual bond lengths has `(n−2)` virtual angles and `(n−3)`
virtual dihedrals — `2n−5` internal degrees of freedom. The builder reads `2n−3` torsions
(`phi[0]`, `phi[n-1]`, `psi[n-1]` are never read), of which `psi[0]` is a rigid rotation about
an axis through CA[0] and therefore carries no RMSD. That leaves `2n−4` effective parameters
mapping onto a `2n−5`-dimensional shape space, so the kernel is `3 + 1 + 1 = 5`. Measurement
and arithmetic agree on every target.

The fifth mode is a **distributed, target-specific crank**: interior `phi` and `psi` loadings
of opposite sign whose sums approximately cancel (e.g. 1CEK: Σδφ = −1.934, Σδψ = +1.894;
1A13: −1.908 / +1.738). It is the whole-chain generalisation of the peptide-plane compensating
pair, and it is exactly inert.

> **Correction to the record: `(n−2)·log2(k)` over-counts the live qubits by one *direction*.**
> The four inert torsions are inert *individually*; a fifth *combination* is inert as well. The
> qubit count is unaffected (a combination is not a coordinate and cannot be dropped from a
> register), but any statement of the form "the CA trace has `2n−4` degrees of freedom" is
> wrong by one, and any Fisher/QNG-style metric on this parameterisation is singular in **five**
> directions, not four. `s15/qgeom_*` should read this row.

### 2.3 How much error the quiet subspace can absorb — DEMONSTRATED

An error placed isotropically in a set `S` of singular directions does RMS damage
`sqrt(mean_{i∈S} sᵢ²)` per unit angle, so the tolerated angle scales inversely.

| | value |
|---|---:|
| RMS singular value of the quiet half ÷ that of the loud half | **0.049** |
| angular budget of the quiet half ÷ that of a random direction | **15.69×** (median 14.16, p10–p90 11.2–22.8) |
| dimension of the exactly-null space as a fraction of `2n` | 0.199 |

**In the null directions the budget is unbounded at first order.** So the answer to "how much
of the error budget can be absorbed there" is: *all of it, geometrically*. Nothing in the
geometry stops a channel from being wrong by 70° in the right directions and emitting 0.35 Å.
Whether anything can be *steered* there is Task 2, and it is a different question.

### 2.4 The 0.565, decomposed — ORACLE DIAGNOSTIC

Share of `‖e‖²` by spectral region, per channel. An **isotropic** error gives
0.250 / 0.500 / 0.250 / 0.199 in the four columns by construction, so the columns read
directly against those numbers.

| channel | align | top ¼ | quiet ½ | bottom ¼ | null space | damage in quiet ½ | RMS deg |
|---|---:|---:|---:|---:|---:|---:|---:|
| retrieval top-75 circular mean | 0.738 | 0.162 | 0.560 | 0.383 | 0.327 | 0.009 | 63.6 |
| top-75 medoid window | 0.693 | 0.160 | 0.567 | 0.367 | 0.315 | 0.009 | 68.0 |
| constant α-helix (control) | 0.709 | 0.185 | 0.499 | 0.301 | 0.257 | 0.009 | 72.0 |
| incumbent projected torsions | 0.665 | 0.153 | 0.578 | 0.391 | 0.342 | 0.015 | 75.1 |
| **ORACLE best pool window** | **0.565** | **0.121** | **0.634** | **0.443** | **0.386** | 0.024 | 60.5 |

Three readings, and the third is the one that constrains Task 2.

1. **The 0.565 is a modest displacement, not a dramatic one.** The best fragment puts 0.121 of
   its squared error in the loud quarter against an isotropic 0.250, and 0.386 in the exact
   null space against an isotropic 0.199 — roughly a factor of two in each direction. The
   alignment discount is bought by moving about **one eighth of the error budget** out of the
   loud quarter.
2. **Every channel is already aligned, including the zero-information control.** A constant
   α-helix scores 0.709 and puts 1.29× the isotropic share in the null space. Alignment below
   1 is a property of *any* smooth, chain-like torsion field, not evidence of a good channel.
   This is the same trap as `a zero-information constant alpha-helix beats uniform random`.
3. **The damage column is the honest one and it is brutal.** Even for the best-aligned channel,
   **97.6% of the RMSD damage comes from the loud half**, which carries only 0.366 of the
   error. The quiet half is where the error *is*; the loud half is where the RMSD *is*. Any
   intervention has to move error out of a region that already holds only a third of it.

### 2.5 The quiet subspace is identifiable NATIVE-FREE — DEMONSTRATED

The subspace half of the alignment statistic needs only `J`, and `J` can be evaluated at the
working structure instead of the native one. Comparing the bottom-half right-singular subspace
of `J(native)` with that of `J(incumbent emitted torsions)`, 126 targets:

| | mean cos² of the principal angles |
|---|---:|
| bottom-**half** subspace, native vs emitted | **0.827** (p10–p90 0.774–0.885) |
| bottom-**quarter** subspace, native vs emitted | 0.735 (0.615–0.842) |
| random-subspace null, same dimensions | **0.499** |
| correlation of the two singular spectra | **+0.980** |

**The geometry does not need the answer.** This is the enabling result for Task 2 — every
intervention in `s15/align_fit.py` evaluates `J` at its own working structure and is therefore
fully native-free. The *error direction* is a different matter and is treated in §5.

*(Caveat kept in view: the minimum principal cosine is 0.170 on average, so the two subspaces
share their bulk but not their edges. A rank-truncated arm that keeps only a few directions is
therefore on much weaker ground than one that reweights smoothly.)*

---

## 3. TASK 2 — THE INTERVENTIONS, AND WHAT THEY WERE MEASURED AGAINST

`python -m s15.align_fit` → `s15/results/align_fit.json`, `align_fit.log`;
`python -m s15.align_reg` → `s15/results/align_reg.json`.

**The common substrate, held identical across every arm.** 126 targets. Restraints are the
shipped leave-fold-out distogram's `expected`, corrected by the leave-fold-out
separation-bias profile of `s15/distcal.py` (K2) and floored at 2.0 Å; weights `1/sd²`;
six native-free starts from `s15.distgeo.starts` under `s15.seed.stable_rng` (so the draw is
reproducible across processes — A2); the exact analytic O(n) torsion gradient; L-BFGS-B to
`ftol 1e-12`; and **multi-start selection on the OBJECTIVE only, never on RMSD**. The control
arm `squared` is therefore exactly `s15/robust.py`'s `squared` arm, which had only ever been
run on an 8-target smoke test.

**The arms, and the null each has to beat.**

| arm | what it changes | its null |
|---|---|---|
| `unit_w` | uniform restraint weights | `squared` (1/sd²) |
| `term<b>_t<t>` | restraints touching the outer `b` residues down-weighted by `t` | `squared`; `t = 1` *is* `squared` |
| `jacw_a<a>` | `w *= (‖∂d_p/∂θ‖ / median)^a` — crude response weighting | `squared` and `unit_w` |
| `damp_k<k>` | `w /= (1 + k (a_p/median)²)`, `a_p = ‖J M⁻¹ g_p‖` = coordinate damage per unit restraint error | `squared` |
| `rank_f<f>` | solve inside `θ_start + span(top-r right singular vectors of J)` | the full-space fit |
| `cauchy_s`, `gnc_cauchy` | robust / graduated-non-convexity losses (K7) | `squared` |
| `tik / loud / quiet` | a quadratic penalty `λ δᵀAδ` on the displacement from the start, with `trace(A) = 2n` in all three so the **shrinkage is matched and only the metric differs** | `tik` (isotropic) is the null for `loud` and `quiet` |

Every Jacobian is evaluated at the arm's own **stage-1 working structure**, never at the
native, so every arm above is **native-free**. Hyperparameters (`t`, `a`, `k`, `f`, `s`, `λ`)
are assembled **leave-fold-out**: the value reported on held-out fold *f* is the one with the
best mean RMSD on the other four folds. The in-fold best is reported beside it as
`*_INFOLD_ORACLE`, which separates a transfer failure from a functional-form failure.

**The ORACLE ceilings**, which can never be headline results:

* `ORACLE_kill_loud_<f>` — take the `squared` solution's **true** torsion error and delete the
  component lying in the loud (top-`f`) subspace. Confounded: deleting a component also shrinks
  `‖e‖`.
* `ORACLE_kill_loud_norm_<f>` — the same, **rescaled to the original ‖e‖**. This is the honest
  ceiling: the *same angular error budget*, rotated out of the loud subspace. It is the single
  number that says how much alignment engineering could ever be worth.
* `ORACLE_term_native` — the outer two residues' torsions replaced by the native ones.

*Defect in my own ORACLE arm, recorded rather than hidden.* `ORACLE_kill_loud_norm` preserves
`‖e‖` **before** the torsions are re-wrapped to `(−π, π]`. When the rescale factor pushes a
component past π the wrap folds it back, so the arm's measured torsion RMS is not exactly
unchanged — it falls by about 5° (7%) rather than by 0°. The residual is reported in the
`dTORSION` column of every table below and must be read with the claim; it is far too small to
account for the RMSD effect (on the INFO phase surface, 75° → 70° at coverage 1 is worth about
0.02 Å), but the arm is "nearly norm-preserving", not exactly so.

---

## 4. TASK 2/3 RESULTS — THE MECHANISM IS HUGE AND NOTHING NATIVE-FREE REACHES IT

`s15/results/align_fit.json`, `align_fit.log`, `align_axes.log`. **126 targets**, 108 min,
`squared` control = **3.677 Å** (incumbent 3.204 on the same targets, so this whole family
loses to the incumbent by **+0.473 [+0.298, +0.648]**, reproducing K1-CORRECTED's sign).

### 4.1 The headline table

RMSD, both raw-accuracy axes, and alignment, all paired against `squared`, all n = 126, all
concentration verdicts **PASS** against the simulated uniform-effect null. `dALIGN` positive
= alignment got *worse*.

| arm | RMSD | median | dRMSD vs squared | dTORSION (deg) | dDISTMAE | dALIGN | verdict |
|---|---:|---:|---|---|---:|---:|---|
| **ORACLE_kill_loud_0.5** | **1.443** | 1.391 | **−2.234 [−2.462,−2.018]** W/L 124/2 | −22.01 [−24.34,−19.77] | −1.100 | −0.518 | confounded (‖e‖ shrinks) |
| **ORACLE_kill_loud_norm_0.5** | **1.855** | 1.745 | **−1.822 [−2.037,−1.613]** W/L 121/5 | −6.84 [−8.30,−5.44] | −0.785 | −0.387 | **the ceiling** |
| ORACLE_kill_loud_0.25 | 3.280 | 3.328 | −0.397 [−0.647,−0.143] | −9.71 | −0.008 | −0.203 | confounded |
| ORACLE_kill_loud_norm_0.25 | 3.547 | 3.533 | −0.130 [−0.392,+0.146] | −2.86 | +0.205 | −0.146 | null |
| **ORACLE_term_native** | 3.677 | 3.718 | **−0.000 [−0.068,+0.066]** W/L 65/61 | **−17.03** | −0.002 | +0.170 | **null** |
| `jacw_a-0.5` | 3.599 | 3.470 | **−0.078 [−0.152,−0.012]** W/L 67/59 | +0.18 [−1.40,+1.77] | −0.006 | −0.024 | LEVER? — see 4.4 |
| `cauchy1_x_term` | 3.591 | 3.628 | −0.086 [−0.192,+0.012] | −0.65 | −0.028 | +0.060 | null |
| `damp_k2.0` | 3.616 | 3.610 | −0.061 [−0.143,+0.016] | +0.72 | +0.013 | +0.034 | null |
| `term2_t0.25` | 3.633 | 3.504 | −0.044 [−0.125,+0.033] | +1.42 | −0.022 | +0.054 | null |
| `rank_f0.75` | 3.640 | 3.545 | −0.037 [−0.091,+0.007] | +0.21 | −0.015 | +0.044 | null |
| `cauchy_s2.0` | 3.645 | 3.586 | −0.032 [−0.116,+0.051] | +1.13 | −0.005 | +0.020 | null |
| `gnc_cauchy` | 3.652 | 3.483 | −0.025 [−0.126,+0.073] | +2.23 | −0.008 | +0.061 | null |
| `unit_w` | 3.678 | 3.539 | +0.001 [−0.090,+0.095] | +2.29 | −0.013 | +0.083 | null |
| `rank_f0.5` | 3.712 | 3.639 | +0.034 [−0.042,+0.116] | −1.71 | +0.011 | +0.033 | null |
| **`term2_t0.0`** | 3.888 | 3.841 | **+0.210 [+0.055,+0.369]** W/L 53/73 | −2.00 | +0.301 | +0.104 | **HURTS** |
| `JACW_lfo` (leave-fold-out) | 3.643 | 3.634 | −0.034 [−0.106,+0.031] | — | — | — | null |
| `CAUCHY_lfo` | 3.666 | 3.569 | −0.011 [−0.102,+0.077] | — | — | — | null |
| `TERMINAL_lfo` | 3.663 | 3.551 | −0.014 [−0.091,+0.063] | — | — | — | null |
| `RANK_lfo` | 3.640 | 3.545 | −0.037 [−0.091,+0.007] | — | — | — | null |
| `DAMP_lfo` | 3.616 | 3.610 | −0.061 [−0.143,+0.016] | — | — | — | null |
| **`ALL_ALIGN_lfo`** (pick the best alignment intervention leave-fold-out) | 3.670 | 3.688 | **−0.007 [−0.074,+0.055]** | — | — | — | **null** |

**The cascade, for the `squared` arm** (the brief requires it for every architecture). This
family has no aggregation stage, so `A = F`:

| stage | value | gap |
|---|---:|---|
| **G** generation — the retrieval start the fit is handed (`START_of_winner`, §7) | 3.877 | |
| **S** selection — six starts, objective-only argmin, then the fit | **3.677** | `G−S` = **−0.200** |
| **A/F** aggregation — none in this family | 3.677 | `S−A` = 0.000 |
| *(ORACLE, not part of the cascade)* alignment ceiling | 1.855 | −1.822 |

**The largest gap by an order of magnitude is the one no stage occupies**: the 1.822 Å sitting
between the emitted structure and the same structure with its error rotated. That is this
workstream's statement of "the largest gap is the next research target".

### 4.2 The mechanism is worth 1.82 Å and it is not an artefact of error size — ORACLE DIAGNOSTIC

Rotating the `squared` fit's **own true torsion error** out of the loud half of the spectrum,
**at essentially the same error magnitude**, takes the fit from **3.677 Å to 1.855 Å**:
**−1.822 [−2.037, −1.613]**, W/L **121/5**, median −1.844, all five folds between −1.50 and
−2.03, concentration **PASS** (p = 0.68). The residual magnitude change is −6.8° of torsion
RMS (an artefact of re-wrapping, §3), which on the INFO phase surface at coverage 1 is worth
about **0.02 Å** — three orders of magnitude too small to carry the effect.

**So the alignment mechanism is real, it is enormous, and it is *not* a restatement of
accuracy.** At today's restraint quality, 1.855 Å — below the project's 2.0 Å "exceptional"
line and 1.349 Å better than the incumbent — is sitting inside the fit's own error, purely as
a matter of which direction that error points. This is the strongest ORACLE ceiling this
workstream produces and, per the brief, it can never be a headline predictive result.

### 4.3 …and nothing native-free gets to it — DEMONSTRATED (a negative)

Six intervention families, twenty-two arms, 126 targets:

* **The best leave-fold-out alignment arm is `ALL_ALIGN_lfo` at −0.007 Å [−0.074, +0.055].**
  Choosing among the terminal, Jacobian-weighted, damage-damped and rank-restricted arms by
  leave-fold-out training-fold RMSD is worth **nothing**.
* **No individual family's leave-fold-out arm has a CI excluding zero** (`DAMP_lfo` −0.061
  [−0.143,+0.016] is the closest).
* **Every arm moves alignment by at most 0.104**, against the **0.387** the ORACLE ceiling
  moves it. The interventions reach roughly a quarter of the required displacement, and several
  move it the *wrong way*.
* **The single arm that clearly hurts is the one the phase diagram most strongly predicted
  would help** (`term2_t0.0`, +0.210 [+0.055,+0.369]) — see 4.5.

### 4.4 The one arm that passes the honest test, and why it is not enough

`jacw_a-0.5` — restraint weights `1/sd² × (‖∂d_p/∂θ‖ / median)^(−0.5)`, i.e. down-weight the
restraints the structure responds to most strongly — is the **only** arm in the study that
improves RMSD (**−0.078 [−0.152, −0.012]**) while its two raw-accuracy axes stay flat
(dTORSION +0.18 [−1.40,+1.77]; dDISTMAE −0.006) and its alignment improves (dALIGN −0.024).
Per-fold: −0.067, −0.114, −0.046, −0.041, −0.116 — **5/5 folds the same sign**.

Everything that argues against calling it a result, stated in one place:

1. It is **1 arm out of 22**. At 22 comparisons a 95% interval that barely excludes zero is
   the expected outcome under the null; nothing about it was pre-registered.
2. **W/L is 67/59 with a median of only −0.013 against a mean of −0.078.** By this project's
   own early-warning rule that is the concentration signature; the null-calibrated check
   returns **PASS (p = 0.718)** but with **mean/sd = −0.19**, which the same helper flags as
   *"test has little power"*. The verdict is therefore "not shown to be concentrated", not
   "shown to be uniform".
3. **Its own leave-fold-out twin loses the effect**: `JACW_lfo` = −0.034 [−0.106, +0.031],
   because fold 4's training folds prefer `a = +0.5` and that fold then pays +0.071. The rule
   picks `a = −0.5` on four folds of five; one fold destroys the significance.
4. The effect is **0.078 Å = 2.1%** of the arm's own RMSD and **4.3%** of the ORACLE ceiling
   it is trying to reach.

**Tier: HYPOTHESIS**, explicitly not DEMONSTRATED. The honest sentence is: *response-weighted
restraints are the only intervention that behaves the way the alignment story predicts, and the
effect is 4% of the available mechanism and does not survive its own tuning rule.*

### 4.5 The terminal lever is REFUTED at its own ORACLE ceiling — DEMONSTRATED

This is the result I expected most and it is the cleanest negative in the file.

INFO C.2 measured a **1.79 Å** spread between terminal and mid-chain missingness and called
terminal-gap-shaped coverage "the largest lever in the diagram". In the restraint-fit setting
it is worth **zero**, and this is settled by the ORACLE arm rather than by an intervention
failing:

> **`ORACLE_term_native` — replacing the outer two residues' torsions with the NATIVE ones —
> changes the emitted RMSD by −0.000 Å [−0.068, +0.066], W/L 65/61, while cutting torsion RMS
> error by 17.03° [15.60, 18.53].**

Perfect terminal torsions, free, and they buy nothing. Note also that this arm makes the
*alignment* worse (+0.170): removing the terminal error removes error that was already sitting
in the quiet subspace, leaving a smaller but proportionally louder residual. That is the
alignment statistic behaving exactly as its definition says it should, and it is a useful check
that the statistic is not merely tracking `‖e‖`.

The two down-weighting arms agree: `term2_t0.25` is −0.044 [−0.125, +0.033] (null) and
**dropping the terminal restraints entirely costs +0.210 [+0.055, +0.369]** — the terminal
restraints are load-bearing for the *interior*, because a pair `(1, 7)` constrains residues 2–6.

**Why this does not contradict INFO C.2, and what it corrects.** C.2 varied *coverage of a
torsion channel* with the gaps filled from the retrieval pool; the cheapness of terminal gaps
is a statement about where an emitter may be *ignorant*. A restraint fit is not ignorant at the
termini — it already places them about as well as it places anything, and `ORACLE_term_native`
shows there is no accuracy left there to buy. **The gap-shape lever is a property of channels
with holes, not of a solver, and it does not transfer to the generative frame.** Recorded
beside INFO's entry, not over it.

### 4.6 Robust losses do NOT beat squared on the full instrument — REFUTED, triggering `s15/robust.py`'s own falsification condition

`s15/results/robust.json` is an **8-target smoke test** (n = 8, incumbent 2.277); this is the
first time the robust family has been run on the full instrument, start-matched against the
identical control.

| arm | RMSD | vs `squared` |
|---|---:|---|
| `cauchy_s2.0` | 3.645 | −0.032 [−0.116, +0.051] |
| `cauchy_s3.0` | 3.653 | −0.024 [−0.105, +0.051] |
| `gnc_cauchy` | 3.652 | −0.025 [−0.126, +0.073] |
| `cauchy_s1.0` | 3.703 | +0.026 [−0.070, +0.121] |
| **`CAUCHY_lfo`** (leave-fold-out scale) | 3.666 | **−0.011 [−0.102, +0.077]** |

`s15/robust.py`'s docstring states the condition: *"If no robust loss beats `squared` on the
full instrument, then the outlier reading of K7 is wrong: the restraint errors would be diffuse
rather than concentrated, and the deficit is genuine information loss that no loss function can
recover."* **No robust loss beats `squared`.** The condition is met.

The diagnostic detail that makes this more than a null: the robust arms **do** achieve much
better restraint residuals — median |z| of **0.364** (`gnc_cauchy`) and **0.430**
(`cauchy_s2.0`) against `squared`'s **0.469** — and convert none of it into RMSD. This is the
project's standing law (`nothing ranks within the pool`, `better matrix, worse ranking`)
appearing in a *generative* setting for the first time: **a better restraint objective is not a
better structure.**

### 4.7 THE HONEST TEST (Task 3) — pooled over 2,142 (target, arm) pairs

| | value |
|---|---:|
| marginal corr(ΔRMSD, Δtorsion RMS) | **+0.022** |
| marginal corr(ΔRMSD, Δdistance MAE) | **+0.593** |
| marginal corr(ΔRMSD, Δalignment) | +0.127 |
| standardised β: torsion / distMAE / alignment | +0.024 / **+0.586** / +0.070 |
| **partial corr(ΔRMSD, Δalignment ∣ Δtorsion, ΔdistMAE)** | **+0.087** |
| the same with the native-free (self-)Jacobian | +0.088 |

**Within the space of interventions that are actually achievable, alignment is a description,
not a lever.** Almost all of the variation in what these interventions emit is explained by how
accurate the emitted *distances* are (β = +0.586); alignment adds +0.087 of partial correlation
on top. Torsion RMS error explains essentially nothing (+0.022) — an independent confirmation
of INFO's "per-torsion sigma is not a sufficient statistic".

**And yet the ORACLE arm proves the lever exists.** The reconciliation is quantitative, not
rhetorical: the achievable interventions move alignment by ≤ 0.10 where the ceiling moves it by
0.39, and the RMSD response is roughly linear in that displacement (0.10/0.39 × 1.82 ≈ 0.47 Å
if the response were linear and perfectly captured; the arms deliver 0.0–0.08 Å, so they are
not even on that line).

> **The finding, stated so it cannot be over-read:** *alignment is a real, first-order,
> 1.8 Å mechanism (ORACLE), and it is a description rather than a lever for every native-free
> intervention measured here.*

### 4.8 An incidental correction the run produced

The substrate here is the leave-fold-out **separation-debiased** restraints of K2/`distcal`
(`squared` = `s15/robust.py`'s control arm). It emits **3.677 Å** where `s15/distgeo.py`'s
undebiased `pred_invvar_weighted` emits **3.644 Å** on the same 126 targets with the same
starts and optimiser. The two runs are not perfectly matched (this one also floors `dhat` at
2.0 Å, as `robust.py` does), so this is an observation and not a controlled ablation — but it
is the first full-instrument evidence that **K2's separation-bias correction does not buy
accuracy in the fit**, and it points the same way as K6 RESULT 4 (`ls_debias` is marginally
worse than `ls_pred` on every ranking column). Flagged for the coordinator; not resolved here.

---

## 5. IS THERE A NATIVE-FREE VERSION OF THE ALIGNMENT DIAGNOSTIC?

`python -m s15.align_free` → `s15/results/align_free.json`, `align_free.log`. 126 targets.
The brief requires this for every ORACLE diagnostic, and the statistic splits cleanly in two.

**The subspace half is native-free and it works** (§2.5): `J` at the emitted structure gives a
bottom-half subspace overlapping the native one at mean cos² **0.827** against a **0.499**
random null, and spectra correlated **+0.980**. On the fits actually produced here, the
alignment of the true error measured through the emitted-structure Jacobian is **0.777**
against **0.736** through the native Jacobian, correlated **+0.698** across targets.

**The error-direction half is ORACLE by construction — but a native-free surrogate carries real
signal.** Four surrogates for the unknown error direction, each scored two ways:

| native-free surrogate for the error direction | its own alignment | corr with the TRUE alignment | \|cos\| with the true error |
|---|---:|---|---:|
| `θ_fit − pool circular mean` (disagreement with the retrieval channel) | 0.718 | **+0.499 [+0.33, +0.65]** | **0.390** |
| `θ_fit − θ_start` (where the fit moved) | 0.656 | +0.346 [+0.17, +0.51] | 0.354 |
| the multi-start spread | 1.352 | +0.453 [+0.29, +0.61] | 0.197 |
| the pool's own principal direction of disagreement | 0.644 | +0.345 [+0.18, +0.51] | 0.218 |
| *null*: two independent directions in 2n ≈ 26 dimensions | — | — | **0.157** |

**DEMONSTRATED.** The channel-disagreement direction `θ_fit − θ_pool` points at the true error
with **\|cos\| = 0.390 against a 0.157 random baseline** and its alignment tracks the true
alignment at **+0.499** with a CI far from zero. So a native-free *diagnostic* of alignment
does exist, and this is the first time the project has had one.

Two things it does **not** license, both of which the standing ledger predicts:

* \|cos\| = 0.39 means the surrogate explains about **15% of the error direction's variance**.
  `decorrelated-errors-exist-but-are-unusable.md` records that fusion gain goes as the *square*
  of the weaker channel's skill; the same arithmetic applies to steering by a noisy direction
  estimate.
* The multi-start spread has alignment **1.352** — it lies in the *loud* directions, more so
  than a random direction. Multi-start solutions differ in exactly the directions that change
  the shape, which is intuitive in hindsight and means the spread is a **shape-diversity**
  measure and not an error-direction estimate.

*(Recorded as an open lead, not a result: nothing in this workstream tests whether steering by
`θ_fit − θ_pool` improves RMSD. It was measured as a diagnostic. Given §4.7's partial
correlation of +0.087, the prior on it converting is low.)*

---

## 6. TASK 4 — CLOSING THE PHASE DIAGRAM

### 6.1 The alignment-corrected surface, and how much of the over-pricing it explains

Read a channel off the INFO i.i.d. surface at coverage 1 using an **effective sigma**
`σ_eff = σ · α / 0.945` (0.945 is the random-direction null of the alignment statistic at these
dimensions, so a channel with `α = 0.945` reads off the surface unchanged):

| channel | σ (deg) | α | σ_eff | raw surface | alignment-corrected surface | actually emits |
|---|---:|---:|---:|---:|---:|---:|
| retrieval top-75 circular mean | 67.6 | 0.738 | 52.8 | 4.82 | 4.48 | 4.072 |
| incumbent projected torsions | 78.4 | 0.665 | 55.2 | 4.90 | 4.53 | 3.215 |
| ORACLE best pool window | 64.4 | 0.565 | 38.5 | 4.71 | 3.89 | 1.773 |
| `squared` restraint fit (this work) | 78.6 | 0.736 | 61.2 | 4.90 | 4.65 | 3.677 |
| **ORACLE_kill_loud_norm_0.5** | 71.8 | **0.349** | **26.5** | 4.82 | **3.03** | **1.855** |

For the alignment-engineered arm the first-order correction closes **(4.82 − 3.03) / (4.82 −
1.86) = 60%** of the over-pricing — an independent reproduction of INFO C.3's "about two-thirds
of the discount is linear subspace alignment and one-third is nonlinear", obtained on a
different object (a solved fit rather than a retrieved fragment) with a different Jacobian
implementation. **The alignment correction is the single largest known improvement to the phase
surface's predictive validity, and it still leaves 40% unexplained.**

### 6.2 What each level requires, restated with alignment as a free variable

Raw per-torsion RMS error required at coverage 1, as a function of the channel's alignment:

| alignment α | required σ for 3.0 Å | for 2.5 Å | for 2.0 Å |
|---|---:|---:|---:|
| 0.945 (i.i.d. / random direction) | 26.1° | 20.9° | 16.3° |
| 0.738 (best real generative channel) | 33.4° | 26.8° | 20.9° |
| 0.565 (best fragment, ORACLE-selected) | 43.7° | 34.9° | 27.2° |
| **0.349 (measured ORACLE ceiling of alignment engineering)** | **70.7°** | **56.6°** | **44.0°** |

Available raw σ: **67.6°** (best generative channel), **78.6°** (this restraint fit). So:

* At the alignment engineering ceiling, **3.0 Å is reachable at today's raw accuracy**
  (70.7° required vs 67.6° available).
* **2.5 Å needs 56.6°** — a **16%** cut in raw error even with perfect alignment engineering,
  by the first-order surface.
* **2.0 Å needs 44.0°** — a **35%** cut, again with perfect alignment engineering.

But the surface still under-states the discount by 40% (6.1), so the **empirical** statement is
the one to quote, and it is stronger: the alignment-engineered arm **actually emits 1.855 Å**
at 71.8° raw error. **At today's restraint accuracy, sub-2 Å is already inside the fit's own
error vector.**

### 6.3 The verdict, restated against INFO C.4

| level | INFO's verdict (before this work) | ALIGN's verdict (after) |
|---|---|---|
| **3.0 Å** | supported; needs ~10% better effective channel | **supported, and now with two independent routes**: a 10% effective-channel gain, or alignment engineering, which reaches 3.0 Å at today's raw accuracy (α ≤ 0.6 suffices) |
| **2.5 Å** | "hard, not closed"; marginal, and only through selection | **STILL MARGINAL, and the reason has changed.** It is *not* information-limited: the ORACLE arm clears it by 0.65 Å at today's restraint accuracy. It is **control-limited** — no native-free intervention moves alignment more than a quarter of the way. |
| **2.0 Å** | "not supported by any i.i.d.-class channel"; every route needs an oracle | **STILL NOT SUPPORTED as a predictive result, but REFRAMED.** The claim "the information is not there" is now false for the generative route: 1.855 Å is achievable by rotating the fit's *existing* error. What is absent is any native-free mechanism that performs the rotation — the best measured is 0.078 Å of the required 1.822 Å (4%), and it does not survive its own tuning rule. |

> **The one-line restatement for the paper.** *2.0 Å in this frame is no longer an information
> limit; it is a control limit.* Sprint 14 and INFO both closed 2.0 Å by arithmetic on how noisy
> the channels are. This workstream shows the channel is noisy in directions that mostly do not
> matter, that a 1.82 Å discount is sitting in the existing error, and that the open problem is
> **steering**, not **accuracy**. That is a different research target and it is stated here for
> the first time.

**The counter-caveat, which must travel with it.** A ceiling that requires the true error
vector is exactly the kind of arm the brief calls out — *"a model that reaches 2 Å only under
unrealistically perfect constraints is not a 2 Å solution"*. `ORACLE_kill_loud_norm_0.5` is
that kind of arm and is labelled so everywhere it appears. Nothing in this file is a predictive
sub-3.6 Å result; the best native-free arm here is **3.599 Å**, which loses to the incumbent's
3.204 Å by **+0.395 [+0.213, +0.574]**.

---

## 7. THE DECIDING EXPERIMENT — SAME SHRINKAGE, DIFFERENT METRIC

`python -m s15.align_reg` → `s15/results/align_reg.json`, `align_reg.log`. 126 targets, 50 min.

Everything in §4 confounds *how much* an intervention constrains the fit with *which
directions* it constrains. This experiment removes that confound completely. Add to the
objective a quadratic penalty `λ·δᵀAδ` on the fit's displacement from its own start, with three
different `A`, **each normalised to `trace(A) = 2n`** so that on an isotropic random
displacement all three impose exactly the same expected penalty:

* `tik` — `A = I`. Ordinary Tikhonov. Every torsion direction costs the same. **The null.**
* `loud` — `A ∝ JᵀJ`. Penalises motion in the RMSD-loud directions ⇒ confines the fit to the
  quiet subspace.
* `quiet` — `A ∝ s_max²I − JᵀJ`. Penalises motion in the RMSD-quiet directions ⇒ confines the
  fit to the well-conditioned subspace. The **smooth** version of the brief's rank-projection arm.

| arm | RMSD | median | torsion RMS | align | vs `squared` | vs incumbent |
|---|---:|---:|---:|---:|---|---|
| **`tik_LFO` (isotropic — the NULL)** | **3.461** | 3.515 | 69.2 | 0.747 | **−0.216 [−0.345, −0.094]** W/L 80/46 | +0.257 [+0.118, +0.394] |
| `quiet_LFO` (confine to the well-conditioned subspace) | 3.501 | 3.532 | 69.8 | 0.754 | −0.176 [−0.297, −0.054] W/L 77/49 | +0.297 |
| `loud_LFO` (confine to the quiet subspace) | 3.615 | 3.481 | 78.7 | 0.755 | −0.062 [−0.160, +0.035] W/L 81/45 | +0.411 |
| `START_of_winner` (the `λ → ∞` limit) | 3.877 | 3.630 | 70.5 | 0.733 | +0.200 [−0.061, +0.479] | +0.673 |
| `squared` (`λ = 0`) | 3.677 | 3.595 | 78.6 | 0.736 | — | +0.473 |

Head-to-head **at matched λ**, geometry against isotropic (positive = the geometry-aware metric
is worse):

| λ | `loud` − `tik` | `quiet` − `tik` |
|---|---|---|
| 0.3 | +0.012 [−0.042, +0.073] | −0.014 [−0.046, +0.008] |
| 3.0 | +0.022 [−0.053, +0.099] | +0.010 [−0.025, +0.047] |
| **30.0** | **+0.124 [+0.025, +0.228]** | +0.040 [−0.050, +0.131] |

### 7.1 The verdict — REFUTED, and it is my own hypothesis that dies

**At matched shrinkage the isotropic metric is the best of the three, and the alignment-aware
metric is significantly the worst at the λ the tuning actually selects.** Ordering is
`tik ≥ quiet ≥ loud` at every λ tested. There is no λ at which knowing the Jacobian helps.

So the one intervention in this workstream with a clear native-free effect —
**`tik_LFO`, −0.216 Å [−0.345, −0.094], W/L 80/46, median −0.101, 5/5 folds the same sign
(−0.100 to −0.329), concentration PASS (p = 0.664, mean/sd −0.29, low power flagged)** — works
because it is **ordinary regularisation**, not because it exploits the RMSD-quiet subspace. Its
raw torsion accuracy also improves (78.6° → 69.2°) and its alignment gets slightly *worse*
(0.736 → 0.747), so by §4.7's own rule it is **confounded**: it moves the raw axis, which is
precisely what the alignment story says a real alignment lever would not need to do.

**Mechanistically this is the least exotic possible explanation and it is the right one.** The
distance-geometry fit **overfits its restraints**; shrinking it toward its retrieval start
recovers 0.216 Å; and the geometry of *which* directions to shrink is irrelevant.

### 7.2 Two caveats that limit the positive half

* **The λ ladder is not converged.** `λ = 30` is the top of the tested grid and the
  leave-fold-out rule selects it on **5 folds of 5** in the `tik` and `quiet` families. The
  optimum is at or beyond the boundary, so 3.461 Å is *what this grid delivers*, not the
  family's best. A wider ladder must be run before the number is quoted as the family's value.
  (The `λ → ∞` limit is bounded by `START_of_winner` at 3.877 Å, so the curve does turn — the
  optimum is interior, just not inside this grid.)
* **It still loses to the incumbent** by +0.257 [+0.118, +0.394]. This is a −0.216 Å improvement
  to an arm that was +0.473 Å behind. K1-CORRECTED's deficit is roughly halved and not closed.

---

## 8. WHAT I REFUTED, INCLUDING MY OWN HYPOTHESES

1. **My own headline hypothesis — that alignment can be engineered.** Six intervention
   families, 22 arms, 126 targets: the best leave-fold-out arm is **−0.007 Å [−0.074, +0.055]**
   and no family's leave-fold-out arm has a CI excluding zero. **REFUTED for every intervention
   tested**, against an ORACLE ceiling of −1.822 Å that proves the mechanism is there.
2. **…and refuted a second time, in the controlled version.** At matched shrinkage
   (`trace(A) = 2n` in all three families) the **isotropic** metric beats both geometry-aware
   metrics at every λ, and beats the alignment-aware one significantly at the λ the tuner picks
   (**+0.124 [+0.025, +0.228]**). The only native-free gain in the workstream is ordinary
   regularisation. **REFUTED with the confound removed**, which is the version that counts.
3. **That the terminal-gap lever transfers to a solver.** INFO called terminal-gap-shaped
   coverage "the largest lever in the diagram" (1.79 Å). In a restraint fit, **perfect native
   terminal torsions are worth −0.000 Å [−0.068, +0.066]** and dropping the terminal restraints
   costs **+0.210 [+0.055, +0.369]**. **REFUTED in this frame**; INFO's measurement stands in
   its own (channels with holes).
4. **That K7's concentrated restraint violations make robust losses the right fix.** First
   full-instrument test: `CAUCHY_lfo` −0.011 [−0.102, +0.077], `gnc_cauchy` −0.025, no arm
   beating `squared`. This is `s15/robust.py`'s **own stated falsification condition** and it
   fires. **REFUTED**, and the 8-target smoke test that suggested otherwise is superseded, not
   deleted.
5. **That "the information is not there for 2.0 Å" in the generative frame.** The
   alignment-engineered ORACLE arm emits **1.855 Å** from the *same* restraints and the *same*
   error magnitude. The limit is control, not information. **REFUTED as an information claim**
   — while 2.0 Å remains unreached predictively.
6. **That the alignment statistic is a re-description of error size.** `ORACLE_term_native`
   cuts torsion RMS by 17.0° and makes alignment **worse** (+0.170) while leaving RMSD
   unchanged; `ORACLE_kill_loud_norm_0.5` leaves the magnitude essentially fixed and improves
   RMSD by 1.82 Å. The two axes move independently in both directions. **CONFIRMED as
   independent.**
7. **That "four torsions per chain are inert" is the complete statement.** There are **five**
   exactly-null directions on 126/126 targets, and the count is derivable (`2n−4` effective
   parameters onto a `2n−5`-dimensional shape space). The fifth is a distributed
   combination, not a coordinate. **CORRECTED, not refuted** — K0's four coordinate torsions
   are confirmed exactly.
8. **That a native-free surrogate for the error DIRECTION does not exist.** The
   channel-disagreement direction `θ_fit − θ_pool` has \|cos\| **0.390** with the true error
   against a **0.157** random baseline, and its alignment tracks the true alignment at
   **+0.499 [+0.33, +0.65]**. **REFUTED** — a weak but real native-free diagnostic exists.
   Whether it can steer anything is untested.

---

## 9. LIMITATIONS, STATED RATHER THAN BURIED

* **Everything here is measured inside ONE architecture** — the continuous torsion-space
  restraint fit of `s15/distgeo.py`, which itself loses to the incumbent by +0.44 to +0.47 Å.
  The alignment ceiling is a ceiling *of that fit's error*. Nothing says the same 1.8 Å sits
  inside the incumbent retrieval pipeline's error, and I did not measure it there.
* **`jacw_a-0.5` is one arm of twenty-two** and is reported as HYPOTHESIS for that reason
  (§4.4). No multiplicity correction would leave it standing.
* **The `ORACLE_kill_loud_norm` arm is norm-preserving only up to torsion wrapping** (§3), a
  residual −6.8° of RMS, quantified and shown to be ~0.02 Å on the phase surface.
* **The Jacobian is a first-order object** and real errors are 60–80° in magnitude. The
  linearisation reproduces the true Kabsch RMSD to 0.999–1.030 at 0.2–5° (§1) and explains 60%
  of the discount at 72° (§6.1); the other 40% is nonlinear and unmodelled.
* **The concentration checks all return PASS but several have mean/sd < 0.3**, which the helper
  itself flags as low power. "Not shown to be concentrated" is not "shown to be uniform".
* **FAIL18 columns are reported for continuity only**; per the Phase 0 audit, `BAND = 1.5 Å`
  has no derivation and only 9KAR survives every threshold choice. No argument here rests on
  set membership.
* **The substrate uses K2's leave-fold-out separation debias**, which §4.8 suggests is not
  helping. Every arm shares it, so every paired comparison is valid, but the absolute level
  (3.677 Å) is 0.03 Å above the undebiased `s15/distgeo.py` arm.

---

## 10. FOR THE COORDINATOR — WHAT TO ACT ON, IN ORDER

1. **Run the `tik` λ ladder past 30.** It is the only native-free gain the workstream found
   (−0.216 Å [−0.345, −0.094], 5/5 folds) and the grid's optimum is at its boundary. It halves
   K1-CORRECTED's +0.44 Å deficit and it is one cheap sweep. `s15/align_reg.py` takes a wider
   `LAMS` tuple and nothing else changes. Note what it *is*: shrinkage of the distance-geometry
   fit toward its retrieval start — i.e. **the pool prior and the distogram restraints fused
   inside the fit**, which is a different object from K3's channel fusion at the distance level.
2. **K7's outlier reading is refuted on the full instrument.** Robust losses buy nothing
   (§4.6). Family B's redirection "from a uniform hard constraint to a robust one" needs a new
   justification or should be closed.
3. **Five, not four, exactly-inert directions** (§2.2). Any Fisher/QNG metric on `(phi, psi)` is
   singular in five directions; `s15/qgeom_qng.py`'s condition numbers should be read against
   that. The qubit count is unaffected.
4. **The terminal-gap lever does not transfer to a solver** (§4.5). Do not build a
   coverage-shaping arm on INFO C.2 without re-testing it in the frame it will be used in.
5. **`ORACLE_kill_loud_norm_0.5` = 1.855 Å is the sharpest ceiling the sprint has for the
   generative route** and reframes 2.0 Å from an information limit to a control limit (§6.3).
   It belongs in the paper as a ceiling, labelled ORACLE, never as a result.
6. **Open lead, untested:** steering by the native-free error-direction surrogate
   `θ_fit − θ_pool` (\|cos\| 0.390 vs a 0.157 null, §5). Measured as a diagnostic only. The
   prior on it converting is low (§4.7 partial correlation +0.087) and it is the only untried
   member of the family.

---

## INSTRUMENT SELFCHECK, END OF WORKSTREAM

`python -m s12.instrument` → `s15/results/instrument_end_align.log`:

    shipped 3.4540004952559396   pool_best 1.7108244199364904
    top75_best 2.3061526409453816   synthesis_fit 3.2040761603809194   n_zero_recall 18

Unchanged from the start of the workstream. (Cache read, not a pipeline reproduction — Phase 0.)
