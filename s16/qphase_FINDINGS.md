# Sprint 16 — QPHASE workstream findings

**Agent: QPHASE.** The sprint's decisive quantum question: *at what objective quality does
CVaR-VQE begin to beat its classical controls, and where does the real peptide-folding
objective sit relative to that boundary?*

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Instrument reproduced at start, exactly:
`shipped 3.4540004952559396 · pool_best 1.7108244199364904 · top75_best 2.3061526409453816 ·
synthesis_fit 3.2040761603809194 · n_zero_recall 18` (`s12.instrument.selfcheck`, a cache read).

Everything is measured on the **full enumerated registers of 19 targets** — nine at n=9
(4⁹ = 262,144 configurations) and ten at n=10 (4¹⁰ = 1,048,576) — **1.28 × 10⁷ exactly
labelled structures**, so the certified global optimum of every objective is known exactly
at every rung of the quality ladder. `Enum.rmsd` and every RMSD are ORACLE quantities used
for post-hoc scoring and for the deliberately-labelled quality knobs, never inside an arm.

**The statistical unit is the TARGET** throughout: the three seeds are averaged *within* a
target before any interval is formed. Intervals are paired fold-aware bootstraps via
`I.paired`. The empirical false-positive floor of ≈0.08 Å is printed on every row.

Reproduce: `python -m s16.qphase {verify,ladder,ladder_noise,realarms,real,alpha,lr,doing,sched}`
then `python -m s16.qphase_report report`. Every mode checkpoints after each target.

---

## 0. HEADLINE — ρ\*, and where this problem sits

**1. ρ\* for each control (blend family, argmin readout, unit = target, n = 19):**

| control | cost | ρ\* — the LOWEST-ρ rung where CVaR-VQE wins with a CI excluding zero and an effect above the ≈0.08 Å floor | tier |
|---|---|---|---|
| **best-of-N from the untrained circuit** | 8,192 evals (matched) | **ρ\* = +0.355** (−0.507 [−0.692,−0.335], W/L 18/1, DIFFUSE) | ORACLE DIAGNOSTIC |
| **classical Boltzmann reweight of the same circuit's samples, matched entropy** | 8,192 evals (matched), **0 extra** | argmin: **undefined by construction** — a reweighting of a drawn set has the same argmin. Ensemble: **ρ\* ≤ +0.264**, i.e. already crossed at the REAL objective | ORACLE DIAGNOSTIC |
| **exact Boltzmann tilt of `p0`, matched diversity** | 676,055 evals (**83×**) | argmin: **ρ\* = +0.355**. Ensemble: **NO CROSSING** — NULL at 9 of 10 rungs | ORACLE DIAGNOSTIC |
| **simulated annealing, matched evaluations** | 8,192 evals (matched) | **NO CONSISTENT CROSSING.** Blend: never significantly better, significantly **worse** at 6 of 10 rungs. Noise: null at 8 of 9 rungs with one isolated win at ρ = +0.980 | ORACLE DIAGNOSTIC |
| **simulated annealing at ¼ the budget** | 2,048 evals (**0.25×**) | **NO CONSISTENT CROSSING.** Blend: null at 9 of 10 rungs, every effect ≤ floor; one isolated cell at ρ = +0.962 (−0.085 [−0.165,−0.006], **0.005 Å above the floor**). Noise: one isolated cell at ρ = +0.500 which **FAILS** the null-calibrated concentration check | ORACLE DIAGNOSTIC |
| **greedy 1-opt with restarts** | 8,192 evals (matched) | **NO CONSISTENT CROSSING.** Blend: VQE significantly **worse at 10 of 10 rungs**. Noise: 2 significant losses, 1 isolated win at ρ = +0.500, rest null | ORACLE DIAGNOSTIC |
| uniform random (**WEAK CONTROL, proves nothing**) | 8,192 evals | ρ\* = +0.355 | — |

**2. ρ IS NOT THE GOVERNING AXIS, and the two independent quality knobs prove it.** The
`blend` and `noise` families dissociate ρ from where the objective's optimum sits. They give
**different ρ\* for the same control**: blend says ρ\*(untrained) = +0.355; noise says the
VQE already beats untrained at ρ = **+0.140** (−0.165 [−0.325,−0.012] SIG) and is negative at
ρ = +0.070. **Disagreement between the two parameterisations is itself the finding**, and the
brief pre-registered that reading. A cell-level fit of the paired difference on ρ has
**R² = 0.014** and puts the crossing at **−2.74 (blend) / −1.10 (noise)** — both outside the
possible range of a correlation, and differing by 2.5×.

**3. WHAT DOES GOVERN IT, measured.** The single coordinate on which both families agree is
**how much better the objective's own certified argmin is than what the control already
returns**, `Δ = RMSD(argmin of E) − RMSD(control)`:

| axis | Spearman with the paired diff | R² pooled | R² blend | R² noise |
|---|---|---|---|---|
| ρ global | −0.059 | 0.014 | 0.007 | 0.006 |
| ρ in-tail (best 1%) | −0.301 | 0.066 | 0.020 | 0.008 |
| argmin percentile | +0.408 | 0.038 | 0.070 | 0.000 |
| **`RMSD(argmin E) − RMSD(control)`** | **+0.636** | **0.427** | **0.859** | 0.124 |

On the structured (blend) landscape the fit is `diff = +0.051 + 0.904 · Δ`, R² = **0.859**:
**the CVaR-VQE returns a structure 90% of the way from what the untrained circuit's best-of-N
returns to the objective's certified argmin.** Classical greedy at the same budget returns
the certified argmin *exactly* (objective gap 0.00000, reached in **68–100%** of cells against
the VQE's **0–32%**). So the VQE is a partial, expensive implementation of "find the
objective's optimum", and the classical arms implement it completely.

**4. THE VQE-vs-CLASSICAL GAP DOES NOT DEPEND ON OBJECTIVE QUALITY AT ALL.** Across the whole
ten-rung blend ladder, ρ = +0.264 → +1.000:

| comparison | mean over rungs | sd over rungs | min | max |
|---|---|---|---|---|
| vqe − anneal (matched evals) | **+0.064** | 0.035 | +0.018 | +0.143 |
| vqe − anneal at ¼ budget | **−0.043** | 0.041 | −0.130 | +0.022 |
| vqe − greedy (matched evals) | **+0.098** | 0.037 | +0.052 | +0.171 |

It is a flat line, and it never crosses zero. **Improving the objective does not move the
quantum-versus-classical comparison; it moves both arms together.**

**5. WHERE THE REAL OBJECTIVES SIT** (19 targets, full register, the identical arm table —
these rows are **not** oracle diagnostics; the objectives are native-free):

| objective | class | ρ global | ρ in-tail 1% (z) | argmin pct | **RMSD of its own argmin** | vqe − untrained | vqe − anneal | vqe − greedy |
|---|---|---|---|---|---|---|---|---|
| `S14_disto_bayes` | SELECTIVE | **+0.485** | +0.048 (+8.4) | 0.282 | 3.077 | −0.183\* | −0.005 | −0.010 |
| `E_ls_pred` | GENERATIVE | +0.449 | +0.076 (+9.7) | 0.283 | 3.123 | −0.144 | +0.041 | +0.005 |
| **`hamil`** (deployed) | STRUCTURAL | +0.264 | +0.083 (+8.6) | **0.231** | **2.661** | −0.309 | +0.030 | +0.121\* |
| `amber_total` | SELECTIVE | +0.148 | +0.011 | 0.551 | 3.931 | *(not an arm — 933-row stratum)* | | |
| `legacy` | SELECTIVE | +0.105 | −0.042 (−2.1) | 0.538 | 4.085 | +0.128 | +0.010 | +0.056 |
| `prior` | SELECTIVE | +0.005 | −0.011 (−0.4) | 0.389 | 3.463 | −0.665\* | −0.050 | +0.131\* |
| *space reference* | | | | | best **1.030**, mean **4.003** | | | |

**THE TWO AXES DISAGREE, and the disagreement is the most important row in this report.**
`S14_disto_bayes` has nearly **twice** `hamil`'s global ρ and a **worse** optimum (3.077 vs
2.661 Å). `prior` has ρ = **+0.005** — no ordering skill whatsoever — and yields the
**largest** VQE-over-untrained gain in the table (−0.665 Å). Across all 95
(objective × target) cells the paired difference correlates **−0.038** with ρ global,
**+0.519** with the argmin percentile — reproducing Sprint 15's +0.52…+0.69 on a new
instrument and on real objectives — and **+0.867** (Pearson +0.911) with
`RMSD(argmin) − RMSD(untrained)`.

**So: this problem does not sit on one side of a ρ boundary, because there is no ρ boundary.
It sits at `RMSD(argmin of the best real objective) = 2.661 Å` against a space best of
1.030 Å — and a classical 1-opt with restarts reaches that 2.661 Å exactly, at the same
budget, in 100% of cells.** The binding constraint is the *location of the objective's
optimum*, and no amount of variational machinery can move it.

---

## V. VERIFICATION — before any claim

`python -m s16.qphase verify`, `s16/results/qphase_verify.json`.

| check | result |
|---|---|
| hoisted blend `(1−s)u_base + s·u_truth` vs `V.blend_objective` at s = 0, 0.37, 1 | max abs diff **0.0, 0.0, 0.0** |
| hoisted noise vs `V.noisy_truth` | max abs diff **0.0** |
| `signal = 0` preserves the base ranking / `signal = 1` the truth's | ρ = **1.000 / 1.000** |
| rebuilt CA traces vs the stored `Enum.rmsd` label | max abs err **3.40e-07** (the label's float32 spacing at 3 Å is 2.4e-07) |
| **a reweighting of a drawn set cannot beat its argmin** | asserted TRUE; the resampled support is a subset of the untrained draws |
| exact `p0 = probs(θ₀)` sums to 1, entropy | 1.0, **12.99 bits** |
| closed-form `E[#distinct in B draws]` vs simulation (B = 2,048) | 1526.41 vs **1532.25** |
| `s12` instrument constants | reproduced to the printed digit |
| **independent reproduction of Sprint 15**: per-target ρ of `E_ls_pred` and `S14_disto_bayes` on 1CS9 | **−0.290 / −0.306**, bit-matching `s15/results/qrestraint.json` |
| **independent reproduction of Sprint 15**: `amber_total` ρ global on the binding-rule stratum, 19 targets | **+0.148**, matching qrestraint's +0.148 |

Seeds come from `SD.stable_rng` (blake2b) with salt `s16qphase`; `hash()` appears nowhere.

**Tie discipline.** Every argmin readout averages the returned RMSD over the *full tied
argmin set*, and the tie multiplicity is stored. The ties are real and structured: at
`signal = 1` every objective has **16** tied global minimisers at n=9, which is exactly the
dead-qubit structure `s14/vqe_lib.py` documents (residue n−1 invisible to the CA trace ×
`phi[0]` inert = 4 × 4). Taking `np.argmin` there would have silently selected one of sixteen
by index order.

---

## 1. DESIGN — what changed, what was held fixed

**One variable moves: the objective's quality.** Everything else is Sprint 14/15 apparatus,
unmodified — the same `Enum` index algebra, the same `mps2f` ansatz, the same `Counter` hard
budget, the same `s14.vqe_run.run` CVaR-VQE with the corrected `baseline="const"` estimator,
the same classical searches, the same two-axis reporting.

* budget **8,192** objective evaluations, **512** shots → 16 iterations, α = **0.25**, lr 0.15
* three seeds (0, 1, 2); θ₀ is a function of (n_qubits, seed) only, so **the control arm is
  literally the circuit the VQE starts from**, shared across every objective
* **blend ladder**: 10 rungs, `signal ∈ {0, 0.1, …, 0.9, 1.0}` over `base = hamil`
* **noise ladder**: 9 rungs, `sigma ∈ {4, 2, 1, 0.5, 0.3, 0.2, 0.12, 0.06, 0}` on the
  uniformised truth
* 19 targets × 19 rungs × 3 seeds × 8 arms = **8,664 arm runs**, 1,083 of them CVaR-VQE

**Why the full register and not Sprint 15's 12-qubit sub-register.** On a
4,096-configuration space any budget large enough to train a VQE already enumerates the
space, so every argmin comparison degenerates. On the full register 8,192 evaluations is
3.1% (n=9) or 0.8% (n=10) of the space — a real search problem. This changes the *cost
ranking of the reweighting control*, and that is stated wherever it matters: Sprint 15's
exact tilt was **200× cheaper** than its exact-gradient VQE on 4,096 configurations; here the
exact tilt reads all 262,144 / 1,048,576 entries and is **83× more expensive** than the
sampled VQE. The direction of that comparison is a property of the register size and must
not be carried across sprints.

**Every rung above `signal = 0` / `sigma = 0` is an ORACLE DIAGNOSTIC** — it reads the native
to build the knob. The BOUNDARY is the deliverable. No blended or noised arm is a predictive
result, and none is quoted as one.

**The realised ρ is measured at every rung** and is the only quality axis used. The mixing
weight is never used as a quality axis, and the docstring's warning is confirmed: on 1CS9
`signal = 0.2` buys ρ = −0.027 while on 2MK7 the same weight buys ρ = +0.524.

---

## 2. THE BOUNDARY — blend family, argmin readout — ORACLE DIAGNOSTIC

Unit = target, n = 19, seeds averaged within target. **Positive = the VQE is WORSE.**
`apct` is the true-RMSD percentile of the objective's certified argmin.

### 2a. vs best-of-N from the untrained circuit (matched 8,192 evals)

| signal | ρ global | ρ tail | apct | VQE | control | diff | median | CI95 | W/L | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| **0 (REAL `hamil`)** | +0.264 | +0.083 | 0.231 | 2.782 | 3.091 | −0.309 | −0.232 | [−0.649,+0.037] | 14/5 | **NULL** |
| 0.1 | +0.355 | +0.468 | 0.006 | 1.645 | 2.152 | **−0.507** | −0.465 | [−0.692,−0.335] | 18/1 | SIG, DIFFUSE |
| 0.2 | +0.455 | +0.489 | 0.005 | 1.556 | 2.046 | −0.490 | −0.435 | [−0.661,−0.335] | 18/1 | SIG, DIFFUSE |
| 0.3 | +0.562 | +0.509 | 0.004 | 1.504 | 1.986 | −0.481 | −0.430 | [−0.674,−0.303] | 17/2 | SIG, DIFFUSE |
| 0.4 | +0.675 | +0.521 | 0.003 | 1.477 | 1.918 | −0.441 | −0.304 | [−0.625,−0.279] | 18/1 | SIG, DIFFUSE |
| 0.5 | +0.787 | +0.543 | 0.003 | 1.444 | 1.879 | −0.435 | −0.333 | [−0.607,−0.280] | 18/1 | SIG, DIFFUSE |
| 0.6 | +0.879 | +0.554 | 0.002 | 1.401 | 1.834 | −0.432 | −0.251 | [−0.595,−0.290] | 18/1 | SIG, DIFFUSE |
| 0.75 | +0.962 | +0.567 | 0.001 | 1.366 | 1.760 | −0.394 | −0.282 | [−0.539,−0.261] | 18/1 | SIG, DIFFUSE |
| 0.9 | +0.995 | +0.610 | 0.000 | 1.327 | 1.675 | −0.348 | −0.265 | [−0.485,−0.224] | 18/1 | SIG, DIFFUSE |
| 1.0 (ORACLE) | +1.000 | +1.000 | 0.000 | 1.216 | 1.469 | −0.252 | −0.188 | [−0.380,−0.146] | 16/3 | SIG, DIFFUSE |

**ρ\*(untrained) = +0.355 on this family.** Note the shape: the gain **peaks at ρ ≈ 0.36 and
then falls monotonically** to −0.252 at a perfect objective. The better the objective, the
*less* the VQE adds over simply drawing from its own untrained circuit and keeping the best —
because best-of-N gets better too. That is the opposite of the naive "better objective ⇒
bigger quantum payoff" reading of Sprint 15's mechanism, and it is a genuine correction to
the expectation this experiment was built on.

### 2b. vs simulated annealing and greedy — **NO CROSSING ON THIS FAMILY**

| signal | ρ global | vqe − anneal (8,192) | vqe − anneal_q (2,048) | vqe − greedy (8,192) |
|---|---|---|---|---|
| 0 (REAL) | +0.264 | +0.030 null ≤FLOOR | −0.130 null | **+0.121 SIG** |
| 0.1 | +0.355 | **+0.081 SIG** | −0.032 null ≤FLOOR | **+0.123 SIG** |
| 0.2 | +0.455 | +0.070 SIG ≤FLOOR | −0.053 null ≤FLOOR | **+0.080 SIG** |
| 0.3 | +0.562 | +0.042 null ≤FLOOR | −0.035 null ≤FLOOR | +0.064 SIG ≤FLOOR |
| 0.4 | +0.675 | +0.060 SIG ≤FLOOR | −0.005 null ≤FLOOR | **+0.082 SIG** |
| 0.5 | +0.787 | +0.031 null ≤FLOOR | −0.016 null ≤FLOOR | +0.053 SIG ≤FLOOR |
| 0.6 | +0.879 | +0.018 null ≤FLOOR | −0.032 null ≤FLOOR | +0.052 SIG ≤FLOOR |
| 0.75 | +0.962 | **+0.097 SIG** | −0.085 SIG | **+0.099 SIG** |
| 0.9 | +0.995 | +0.065 SIG ≤FLOOR | −0.067 null ≤FLOOR | **+0.136 SIG** |
| 1.0 (ORACLE) | +1.000 | **+0.143 SIG** | +0.022 null ≤FLOOR | **+0.171 SIG** |

* **greedy**: the VQE is worse at **10 of 10 rungs**, and the loss is *largest where the
  objective is best*. **ρ\*(greedy) does not exist in [0,1] on this family.**
* **annealing at matched cost**: the VQE is never significantly better; significantly worse
  at 6 rungs (4 of them ≤ floor). **ρ\*(anneal) does not exist in [0,1] on this family.**
* **annealing at a QUARTER of the VQE's objective budget** is statistically indistinguishable
  from the VQE at 9 of 10 rungs, every effect at or below the ±0.08 Å floor; the one
  exception, ρ = 0.962, is −0.085 with a CI upper end of −0.006, i.e. **0.005 Å above the
  false-positive floor** and not a boundary. Four times the objective budget, an ansatz, a
  gradient estimator and a CVaR tail parameter buy **nothing measurable** over 2,048
  evaluations of textbook simulated annealing, at **any** objective quality.

**The noise family qualifies this, and it is reported rather than dropped.** There, against
annealing the comparison is null at 8 of 9 rungs with one significant VQE win at ρ = +0.980
(−0.113 [−0.199,−0.014]); against annealing at ¼ budget one significant win at ρ = +0.500
which **FAILS** the null-calibrated concentration check; against greedy, two significant
losses (ρ = 1.000 and 0.980), one significant win (ρ = +0.500) and six nulls. **Scattered
significant cells of both signs with no monotone structure is what "no boundary" looks like
in the presence of noise, and it is what the noise family shows.**

### 2c. THE OBJECTIVE AXIS, reported separately and never substituted

Mean objective gap in rank units (0 = certified global optimum) / % of cells reaching it exactly:

| signal | vqe | untrained | anneal | anneal_q | **greedy** | tilt_exact | random |
|---|---|---|---|---|---|---|---|
| 0 | 0.00000 / **32%** | 0.00128 / 0% | 0.00000 / 21% | 0.00002 / 0% | 0.00000 / **100%** | 0.00093 / 0% | 0.00016 / 0% |
| 0.3 | 0.00101 / 21% | 0.00760 / 0% | 0.00016 / 32% | 0.00095 / 11% | 0.00000 / **100%** | 0.00463 / 0% | 0.00247 / 0% |
| 0.6 | 0.00103 / **0%** | 0.00732 / 0% | 0.00025 / 26% | 0.00095 / 0% | 0.00008 / **95%** | 0.00380 / 0% | 0.00295 / 0% |
| 1.0 | 0.00025 / 5% | 0.00083 / 0% | 0.00001 / 47% | 0.00022 / 5% | 0.00000 / **68%** | 0.00046 / 0% | 0.00011 / 5% |

**Greedy 1-opt with restarts reaches the certified global optimum in 68–100% of cells at
8,192 evaluations; the CVaR-VQE reaches it in 0–32%.** On the optimisation axis there is no
regime, at any objective quality, in which the quantum arm is competitive. This is Sprint
15's "annealing finds the certified optimum in 100% of cells" reproduced on the *full*
register with a *sampled* (deployable) VQE rather than an exact-gradient one, and with greedy
rather than annealing as the strongest arm.

### 2d. DIVERSITY, reported beside every aggregate (mean distinct configurations of 8,192 draws)

| signal | vqe | untrained | anneal | greedy | tilt_exact | tilt_samples |
|---|---|---|---|---|---|---|
| 0 | 4,761 | 4,994 | 7,265 | 2,899 | 4,516 | 2,937 |
| 0.5 | 4,263 | 4,994 | 6,681 | 4,120 | 4,200 | 2,603 |
| 1.0 | 4,449 | 4,994 | 7,110 | 6,054 | 4,314 | 2,787 |

At α = 0.25 the VQE loses only 12–15% of its distinct-configuration count relative to the
untrained circuit. `tilt_exact` is matched to the VQE's distinct count *by construction*
(closed-form expected-distinct bisection), which is the equal-diversity budget convention.
`tilt_samples` is matched on **entropy**, and the table shows the two conventions are not the
same: entropy-matching leaves it at 2,600 distinct against the VQE's 4,300. **Any comparison
against `tilt_samples` is therefore diversity-confounded in the VQE's favour and is reported
with `tilt_exact` beside it.**

### 2e. vs the classical Boltzmann reweighting

Two conventions, and they give different answers, so both are reported.

* **`tilt_samples` — matched cost, matched entropy, ZERO extra objective evaluations.**
  Reweight the untrained circuit's own 8,192 draws by `exp(−E/T)` with `T` solved so the
  reweighted draw entropy equals the VQE's, and resample. **On the argmin readout this arm
  is identical to the untrained control by construction** (verified): a reweighting of a set
  cannot beat that set's best member. So *no boundary against a matched-cost reweighting is
  definable on the argmin readout at all* — a structural fact, recorded because Sprint 15's
  headline "the effect dissolves under reweighting at 1/200th the budget" was measured on an
  ENSEMBLE readout and does not transfer to an argmin readout.
* **`tilt_exact` — matched diversity, 676,055 evaluations (83× the VQE).** Argmin readout:
  the VQE wins from ρ = +0.355 (−0.407 [−0.552,−0.271]) down to −0.123 at ρ = 1. Ensemble
  readout: **NULL at 9 of 10 rungs**, the single exception (ρ = +0.355) failing the
  null-calibrated concentration check (§3).

---

## 3. THE ENSEMBLE READOUT — reported separately, never substituted

`rand5` coordinate average of 5 configurations drawn uniformly from the arm's own visited
set, no ranker anywhere — Sprint 15's E2 readout.

| signal | ρ global | vqe − untrained | vqe − tilt_samples | vqe − tilt_exact | vqe − anneal | vqe − greedy |
|---|---|---|---|---|---|---|
| 0 (REAL) | +0.264 | −0.241 null | **−0.249 SIG** | −0.156 null | +0.139 null | **+0.334 SIG** |
| 0.1 | +0.355 | **−0.478 SIG** (conc FAIL) | −0.343 SIG (conc FAIL) | −0.247 SIG (conc FAIL) | +0.158 null | **+0.422 SIG** |
| 0.3 | +0.562 | **−0.574 SIG** (conc FAIL) | −0.185 null | −0.186 null | **+0.283 SIG** | **+0.505 SIG** |
| 0.6 | +0.879 | **−0.773 SIG** | −0.196 SIG | −0.102 null | **+0.220 SIG** | **+0.284 SIG** |
| 0.9 | +0.995 | **−0.828 SIG** | −0.137 null | −0.079 null ≤FLOOR | **+0.291 SIG** | **+0.272 SIG** |
| 1.0 (ORACLE) | +1.000 | **−0.762 SIG** | −0.207 SIG (conc FAIL) | −0.128 null | +0.153 null | **+0.303 SIG** |

**Three readings.**

1. **The two readouts have opposite quality dependence.** The VQE's ensemble gain over its
   untrained start *grows* with ρ (−0.241 → −0.828) where its argmin gain *shrinks*
   (−0.309 → −0.252). Reporting one alone would have produced a different headline.
2. **Against a diversity-matched classical tilt the ensemble effect is NULL at 9 of 10
   rungs**, the single exception failing the null-calibrated concentration check. The VQE beats the entropy-matched *sample* reweighting (−0.09 to −0.34) and does
   not beat the diversity-matched *exact* tilt. The mechanism is visible in §2d: the
   sample-reweighting cannot leave the 4,994 configurations the untrained circuit drew and
   loses 45% of its distinct count doing the reweighting, while the VQE and the exact tilt
   both move support. **What the VQE contributes over a reweighting is that its support
   moves; a classical tilt over the same space with matched diversity supplies the same
   thing.** That is a *refinement* of Sprint 15's §7a rather than a contradiction of it —
   §7a's tilt was an exact tilt on a 4,096-configuration space, i.e. this arm, and it too
   was not beaten there.
3. **Both classical searches beat the VQE on the ensemble readout as well**, greedy at every
   rung (+0.27 to +0.62). Sprint 15's "the one place a VQE wins — an ensemble consumed
   without a ranker" **does not survive** on the full register once annealing and greedy are
   in the table.

---

## 4. THE CROSS-CHECK — blend vs noise, and the axis that survives it — DEMONSTRATED

The noise family is the independent second parameterisation. Because `sigma` acts on the
uniformised truth, the realised ρ depends on `sigma` alone and is target-independent to
±0.005 — a *calibrated* quality knob. And it dissociates ρ from the argmin's location: the
argmin of `u_truth + σz` over 262,144 configurations is chosen by the largest negative noise
draw, so at σ = 4 the objective has ρ = +0.070 with its argmin at the **39th** percentile,
while a blend at ρ = +0.264 has its argmin at the **23rd**.

**vs untrained, argmin readout, noise family:**

| sigma | ρ global | ρ tail | apct | diff | CI95 | W/L | verdict |
|---|---|---|---|---|---|---|---|
| 4.0 | **+0.070** | +0.016 | 0.392 | −0.122 | [−0.364,+0.136] | 12/7 | null |
| 2.0 | **+0.140** | +0.036 | 0.340 | **−0.165** | [−0.325,−0.012] | 12/7 | **SIG**, DIFFUSE |
| 1.0 | +0.273 | +0.064 | 0.153 | +0.007 | [−0.268,+0.282] | 10/9 | null ≤FLOOR |
| 0.5 | +0.500 | +0.086 | 0.121 | −0.210 | [−0.354,−0.088] | 15/4 | SIG (conc FAIL) |
| 0.3 | +0.701 | +0.093 | 0.079 | −0.190 | [−0.362,−0.017] | 13/6 | SIG, DIFFUSE |
| 0.2 | +0.832 | +0.115 | 0.050 | −0.092 | [−0.260,+0.071] | 10/9 | null |
| 0.12 | +0.929 | +0.121 | 0.025 | −0.179 | [−0.337,−0.053] | 13/5 | SIG (conc FAIL) |
| 0.06 | +0.980 | +0.149 | 0.015 | −0.158 | [−0.227,−0.090] | 16/3 | SIG, DIFFUSE |
| 0.0 | +1.000 | +1.000 | 0.000 | −0.252 | [−0.376,−0.147] | 16/3 | SIG, DIFFUSE |

**The two families do not agree on ρ\*.** Blend: no significant win below ρ = +0.355. Noise:
a significant win at ρ = +0.140 and a negative point estimate at ρ = +0.070. **Under the
brief's pre-registered reading, that disagreement means the boundary is not about objective
QUALITY as ρ measures it.**

**And the two families DO agree on `Δ = RMSD(argmin E) − RMSD(control)`**: pooled over all
361 (family × rung × target) cells, Spearman +0.636, R² 0.427, with both families' fitted
crossings inside the achievable range (blend −0.057 Å, noise +0.453 Å) where every ρ-based
crossing is outside it. The blend family's fit is `diff = +0.051 + 0.904·Δ`, R² **0.859**.

> **The phase boundary exists, and it is not in ρ. It is at Δ = 0: CVaR-VQE beats a fixed
> sampler exactly when the objective's own optimum is better than what that sampler already
> returns, and it captures ~90% of that difference on a structured landscape and ~25% on a
> noise-dominated one. Because a classical 1-opt captures 100% of it at the same budget,
> there is no objective quality at which the quantum arm is the right tool.**

**Why annealing and greedy have no boundary, mechanistically.** Improving the objective moves
`RMSD(argmin E)` down, and the classical arms *are already at the argmin*. The VQE is at 90%
of the way there. A fixed 10% shortfall of a shrinking gap is a shrinking absolute loss — but
it never changes sign. The measured gap is flat at +0.064 ± 0.035 (anneal) and +0.098 ± 0.035
(greedy) across the entire ladder (§0.4).

---

## 5. WHERE THE REAL OBJECTIVES SIT — DEMONSTRATED (native-free arms)

`python -m s16.qphase realarms`. Five real objectives, each rank-uniformised exactly as every
ladder rung is, each through the identical arm table, 19 targets × 3 seeds.

### 5a. The placement table, absolute values (mean returned CA-RMSD, Å)

| objective | ρ global | ρ tail 1% | argmin pct | **RMSD of argmin** | VQE | untrained | anneal | **greedy** |
|---|---|---|---|---|---|---|---|---|
| `prior` | +0.005 | −0.011 | 0.389 | 3.463 | 3.594 | 4.259 | 3.645 | **3.463** |
| `legacy` | +0.105 | −0.042 | 0.538 | 4.085 | 4.034 | 3.906 | 4.024 | 3.978 |
| **`hamil`** | +0.264 | +0.083 | 0.231 | **2.661** | 2.782 | 3.091 | 2.752 | **2.661** |
| `E_ls_pred` | +0.449 | +0.076 | 0.283 | 3.123 | 3.122 | 3.266 | 3.080 | 3.117 |
| `S14_disto_bayes` | +0.485 | +0.048 | 0.282 | 3.077 | 3.079 | 3.262 | 3.084 | 3.089 |
| *space* | | | | best **1.030** | | | mean **4.003** | |

**Greedy lands on the objective's certified argmin to within 0.012 Å on four of the five
objectives** (0.000, 0.000, 0.006, 0.012); on `legacy` it returns **3.978 Å against that
objective's argmin of 4.085 Å**, i.e. it fails to reach the certified optimum and is *better*
for having failed — the clearest single instance in this report that on a bad objective
optimising harder makes the answer worse. The VQE lands within 0.12 Å of the argmin on all
five. The whole quantum-versus-classical question, on the real objectives, is worth
**≤ 0.13 Å** — against a **1.63 Å** gap between the best real objective's optimum (2.661 Å)
and the space's best structure (1.030 Å).

### 5b. The two quality axes disagree, and only one predicts

95 (objective × target) cells, paired VQE-minus-untrained on the argmin readout:

| candidate axis | Spearman | Pearson |
|---|---|---|
| ρ global | **−0.038** | +0.067 |
| ρ in-tail (best 1%) | −0.254 | −0.226 |
| argmin percentile | **+0.519** | +0.388 |
| RMSD of the argmin | +0.546 | +0.542 |
| **`RMSD(argmin) − RMSD(untrained)`** | **+0.867** | **+0.911** |

The +0.519 against the argmin percentile **reproduces Sprint 15's +0.52…+0.69** on a new
instrument (full register instead of a 12-qubit sub-register), a new estimator (sampled
instead of exact-gradient) and *real* objectives instead of a blend knob. **The mechanism
replicates. The ρ framing of it does not.**

Two concrete demonstrations of the disagreement:

* `S14_disto_bayes` (ρ = +0.485, the best global ranker measured here) has a **worse**
  optimum than `hamil` (ρ = +0.264): 3.077 vs 2.661 Å. Global ρ ranks them in the reverse of
  the order that matters.
* `prior` has ρ = **+0.005** — literally no ordering skill — and produces the **largest**
  VQE-over-untrained gain in the study, −0.665 Å [CI excludes zero]. It does so because the
  retrieval torsion prior is a per-residue product whose optimum is a near-constant
  secondary-structure assignment at 3.463 Å, comfortably better than the 4.003 Å space mean,
  while best-of-N *by an objective with ρ = 0* returns an essentially random 4.259 Å draw.
  **A zero-ρ objective can produce a large "VQE win", and this is exactly the brief's
  standing warning that a zero-information constant α-helix beats uniform random.**

### 5c. In-tail versus global — the tail is where the optimiser concentrates

Every real objective's `tail − bulk` is negative or at chance, as Sprint 15 recorded. In-tail
ρ at the best 1% is +0.048 to +0.083 for the three distogram-derived objectives (z = +8.4 to
+9.7 against matched random-tail nulls, so real but small), **−0.042 (z = −2.1) for Legacy**
and **+0.011 for AMBER** on the binding-rule stratum (`amber_kind == 0 AND amber_idx !=
snap_index`, mean n = 933 rows per target, range 698–1,195). The in-tail axis and the global
axis put `S14_disto_bayes` (global 1st, in-tail 4th) and `hamil` (global 3rd, in-tail 1st) on
opposite orderings, and **neither predicts the VQE result**: in-tail ρ correlates −0.254 with
the paired difference.

---

## 6. WHAT α ACTUALLY IS — a temperature, and a weaker one than a classical thermostat

`python -m s16.qphase alpha`. Six α × 19 targets × 3 seeds × two quality rungs. The decisive
test the brief specifies: **does an entropy-matched classical temperature reproduce α's
entire effect on the readout?** Two classical thermostats are run beside every α:
`tiltmatch` (Boltzmann reweight of the untrained circuit's OWN 8,192 draws, `T` solved so the
reweighted draw entropy equals the VQE's, **zero extra objective evaluations**) and `tiltx`
(exact Boltzmann tilt of `p0` over the whole register, `T` solved so the **expected distinct
count** equals the VQE's, 676,055 evaluations).

### 6a. The α map, at the ORACLE truth objective (signal = 1) — ORACLE DIAGNOSTIC

| α | α·shots | VQE distinct | VQE H (bits) | VQE argmin | VQE rand5 | VQE set-mean | tiltx distinct | tiltx argmin | tiltx rand5 | tiltx set-mean |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 10.2 | 5,116 | 11.72 | **1.139** | 2.928 | 3.559 | 4,722 | 1.367 | 3.010 | 3.679 |
| 0.05 | 25.6 | 5,054 | 11.67 | 1.152 | 2.805 | 3.525 | 4,670 | 1.376 | 3.043 | 3.637 |
| 0.10 | 51.2 | 4,918 | 11.58 | 1.198 | 2.823 | 3.458 | 4,602 | 1.368 | 3.022 | 3.635 |
| 0.25 | 128 | 4,449 | 11.19 | 1.216 | 2.677 | 3.277 | 4,314 | 1.339 | 2.721 | 3.349 |
| 0.50 | 256 | 4,184 | 10.90 | 1.228 | 2.601 | 3.175 | 4,072 | 1.314 | 2.648 | 3.188 |
| 1.00 | 512 | 3,793 | 10.59 | 1.256 | **2.551** | **3.090** | 3,756 | **1.245** | **2.385** | **2.934** |
| *untrained* | — | 4,994 | 11.70 | 1.469 | 3.450 | 4.007 | | | | |

**6b. VERDICT: α is a temperature. DEMONSTRATED.** The α-curve of the VQE's readouts is
reproduced by a one-parameter classical thermostat with Pearson

| readout | VQE vs diversity-matched `tiltx` | VQE vs entropy-matched `tiltmatch` | VQE range | classical range |
|---|---|---|---|---|
| `rand5` coordinate average | **+0.927** | **+0.961** | 0.377 Å | **0.657 Å** |
| drawn-set mean | **+0.979** | **+0.994** | 0.469 Å | **0.745 Å** |
| argmin | −0.857 | −0.944 | 0.117 Å | 0.130 Å |

**On both ensemble readouts a classical Boltzmann temperature reproduces α's entire monotone
effect and moves the readout FURTHER than α does** (0.657 vs 0.377 Å; 0.745 vs 0.469 Å). At
α = 1 the diversity-matched classical tilt is **better** than the VQE on both
(rand5 +0.165 [CI excludes zero]; set-mean +0.157 [CI excludes zero]). At the real objective
(signal = 0) the same correlations are +0.638 (rand5) and +0.941 (set-mean).

**6c. The one thing α does that a temperature does not** is on the ARGMIN readout, where the
two curves run in *opposite* directions (Pearson −0.86): the VQE's argmin degrades as α grows
(1.139 → 1.256) while the tilt's improves (1.367 → 1.245). The effect is 0.117 Å, at the
false-positive floor's scale. **The trade Sprint 15 named — concentration buys the set mean
and pays the set best — reproduces exactly**: over the α range the VQE's set mean improves
0.469 Å while its argmin degrades 0.117 Å, and both are monotone.

**6d. A SCOPE CORRECTION TO THE RECORDED α RANGE.** Sprint 15 measured α moving the
distribution from **2.7 distinct configurations out of 4,096 to 408**, entropy 0.73 → 6.58
bits. Here α moves it from **5,116 distinct out of 8,192 draws to 3,793**, entropy
11.72 → 10.59 bits — a **1.35×** range against a 150× one. Both are correct measurements of
different experiments. The difference is **how far the optimiser is run**, not α: Sprint 15's
figure is 200 iterations of an *exact-gradient* VQE on a 4,096-configuration register; this
is 16 iterations of the *sampled* estimator on a 262,144-configuration one. **α's strength as
a diversity dial is a property of the optimisation budget, and the recorded 2.7-distinct
collapse must not be quoted as a property of α.**

**6e. A confound carried forward, not hidden.** Sprint 15's Defect 2 (sampled-CVaR upward
bias at non-integer `α·N`) is active at α = 0.02, 0.05, 0.10 (α·shots = 10.2, 25.6, 51.2) and
inactive at 0.25, 0.5, 1.0. The best argmin cell in the sweep (α = 0.02) is one where the
bias **is** active. As in Sprint 15, the two cannot be separated on this data and that is
stated rather than resolved.

---

## 7. THE WEAK-OPTIMISER TEST — the recorded defect does NOT reproduce, and its recorded MECHANISM is wrong

`python -m s16.qphase lr`. Five learning rates spanning 16× on the **corrected**
`baseline="const"` estimator, plus the recorded defective `baseline="tail"` estimator at the
reference rate, on 19 targets × 3 seeds at two quality rungs.

### 7a. The diagnostics reproduce; the effect does not

| arm | ‖g‖ first iter | ‖g‖ mean | **‖Δθ‖ displacement** | distinct | H (bits) |
|---|---|---|---|---|---|
| `const` lr 0.15 (rung 0) | **0.21** | 0.03 | **6.74** | 4,761 | 11.42 |
| `tail` lr 0.15 (rung 0) | **0.08** | 0.01 | **6.29** | 5,235 | 11.75 |
| `const` lr 0.15 (rung 1) | 0.15 | 0.03 | 6.40 | 4,449 | 11.19 |
| `tail` lr 0.15 (rung 1) | 0.07 | 0.02 | 6.10 | 5,086 | 11.69 |

**The gradient-norm ratio reproduces exactly** — the biased estimator's gradient is
**2.6× / 2.1× smaller**, against Sprint 15's recorded 2.2–2.5× and the brief's 2.4×.

**But the parameter displacement is only 7% and 5% smaller.** **REFUTED: the recorded
mechanism — "a de-facto step-size reduction" — is not what the biased estimator does.** Adam
divides by a running RMS of the gradient's own magnitude, so it is invariant to a uniform
rescaling of the gradient; a 2.6× smaller gradient therefore produces a nearly identical
step. What `baseline="tail"` changes is the *direction* (Sprint 15 measured cosine +0.656
against the exact reference), not the step length. Any account of the defect that rests on
step size is mechanically unavailable under Adam.

### 7b. The effect itself, at the target unit — NULL

`tail` minus `const` at the identical learning rate, paired over 19 targets:

| rung | argmin readout | rand5 ensemble readout |
|---|---|---|
| 0 (REAL `hamil`) | −0.079 [−0.221, +0.039] W/L 9/5 **NULL** | −0.025 [−0.178, +0.123] W/L 9/10 **NULL** |
| 1 (ORACLE truth) | −0.018 [−0.077, +0.050] W/L 13/6 **NULL** ≤FLOOR | **+0.314 [+0.093, +0.532] W/L 4/15 — the DEFECT is significantly WORSE** |

**REFUTED at n = 19 on the target unit.** Sprint 15's `−0.266` vs `−0.145` (the defect
better) does not reproduce: the two estimators are statistically indistinguishable on the
argmin readout at both quality levels, and the corrected one is **significantly better** on
the ensemble readout at a good objective. Sprint 15's cells were (target × seed × objective)
on the full register with three seeds and a different objective family; this is the same
instrument at the target unit with the two readouts separated.

### 7c. The lr ladder — and Sprint 15's "less optimisation is better" is REFUTED here

Each arm minus the untrained control, 19 targets:

| arm | rung 0, argmin | rung 0, rand5 | rung 1, argmin | rung 1, rand5 |
|---|---|---|---|---|
| `const` lr **0.15** | −0.309 [−0.646,+0.049] | **−0.252** [−0.447,−0.057] | **−0.252** [−0.379,−0.140] | **−0.773** [−0.963,−0.589] |
| `const` lr 0.075 | **−0.334** [−0.676,−0.009] | −0.237 [−0.429,−0.050] | −0.217 [−0.324,−0.132] | −0.717 [−0.931,−0.510] |
| `const` lr 0.0375 | −0.288 [−0.528,−0.042] | −0.189 [−0.354,−0.028] | −0.197 [−0.299,−0.118] | −0.468 [−0.714,−0.227] |
| `const` lr 0.019 | −0.166 [−0.352,+0.011] | −0.117 [−0.304,+0.083] | −0.146 [−0.214,−0.085] | −0.218 [−0.443,+0.006] |
| `const` lr 0.0094 | −0.176 [−0.348,−0.011] | −0.029 [−0.222,+0.155] | −0.097 [−0.149,−0.049] | −0.097 [−0.325,+0.119] |
| `tail` lr 0.15 | −0.388 [−0.722,−0.063] | −0.277 [−0.482,−0.093] | −0.270 [−0.375,−0.185] | −0.458 [−0.726,−0.186] |

Dividing the learning rate by 16 costs **+0.156 [+0.070,+0.253]** on the argmin readout and
**+0.676 [+0.487,+0.855]** on the ensemble readout at a good objective, and
**+0.223 [+0.020,+0.422]** on the ensemble readout at the real one. **The response is
monotone in the optimiser's strength, in the direction that MORE optimisation is better.**

> **VERDICT on the weak-optimiser hypothesis: REFUTED on this instrument.** CVaR-VQE at this
> scale is **not** over-stepping; it is under-stepping, and the recorded "biased estimator
> acts as a step-size reduction" is mechanically impossible under Adam and empirically null
> at the target unit. Sprint 15's mechanism table ("in all three axes, less optimisation is
> better") does not extend from its 12-qubit exact-gradient setting to a sampled VQE on the
> full register. **Note the two settings differ in three ways at once** — estimator, register
> size and iteration count — so this is a scope limitation on the recorded claim, not a
> demonstration that it was wrong where it was measured.

---

## 8. WHAT THE VQE IS DOING — four axes, never substituted

`python -m s16.qphase doing`, nine n=9 targets × 3 seeds × 4 quality rungs, `exact_dist`, so
the whole 2¹⁸ distribution is enumerated rather than sampled. Untrained reference:
entropy **12.74 bits**, `E_p0[E percentile] 0.50`, mode RMSD **3.792**, `P(<2 Å) = 0.0126`.

| signal | ρ | objective gap | **p(certified argmin)** | `E_p[E pct]` | H (bits) | distinct/8192 | **mode RMSD** | `E_p[RMSD]` | **P(<2 Å)** | returned | best seen | cert. optimum |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (REAL) | +0.176 | 3.2e−06 | 0.0051 (p0 2.9e−06) | **0.168** | 11.06 | 4,020 | **3.014** | 3.536 | **0.0459** (**3.6×**) | 2.983 | **1.240** | 2.896 |
| 0.3 | +0.495 | 1.6e−03 | 0.0014 | 0.148 | 9.55 | 3,588 | 2.310 | 3.094 | 0.1309 (10.4×) | 1.700 | 1.177 | 1.563 |
| 0.6 | +0.866 | 1.7e−03 | 0.0049 | 0.182 | 9.23 | 3,551 | 1.863 | 2.923 | 0.2075 (16.5×) | 1.499 | 1.153 | 1.403 |
| 1.0 (ORACLE) | +1.000 | 2.9e−04 | 0.0063 | 0.248 | 10.27 | 3,916 | **1.663** | 2.849 | 0.2685 (21.3×) | 1.144 | 1.144 | 0.969 |

**The classification, one axis at a time.**

* **OPTIMISING — partially, and never to completion.** The trained state places
  **0.14 %–0.63 %** of its mass on the objective's certified global minimiser. That is a
  **440×–1,700× enrichment** over the untrained circuit's 3e−06, and it is still not
  convergence: Sprint 15's "the VQE places 5.9 % of its mass on the optimum where annealing
  finds it in 100 % of cells" is reproduced here an order of magnitude weaker on a 64×
  larger register. Greedy reaches the certified optimum outright in 68–100 % of cells (§2c).
* **SAMPLING — yes, and this is the axis where it clearly works.** `E_p[E percentile]` moves
  from the untrained circuit's **0.50** to **0.15–0.25**. The variational state genuinely
  learns to put mass where the objective is low. **This is the one unambiguous positive.**
* **EXPLORING — no. It explores LESS than doing nothing.** Distinct configurations fall from
  the untrained circuit's 4,994 to 3,551–4,020, entropy from 12.74 to 9.2–11.1 bits, and
  simulated annealing at the same budget visits **7,110–7,265**. Whatever CVaR-VQE is, it is
  not an explorer.
* **REPRESENTING — yes, and it is the largest effect in the study.** Near-native probability
  mass is enriched **3.6× at the real objective** and 10×–21× on the blended ones, against
  the untrained circuit; the **mode** — the only readout for which "the variational state
  learned something" is meaningful — improves from 3.792 Å to **3.014 Å** at the real
  objective and to 1.663 Å at the ORACLE one. This reproduces the direction of Sprint 15's
  5–7× enrichment on generative objectives.

**8a. And the representation is worth nothing through the readout the pipeline uses.** At the
real objective the arm *returns* 2.983 Å — within 0.087 Å of the objective's certified
optimum, 2.896 Å — while the best structure it actually drew is **1.240 Å**. **The selection
gap is 1.74 Å and it belongs entirely to the objective.** This is Sprint 15's 1.47–2.14 Å
arm-level selection gap, reproduced to within its own range, on a different instrument.

**8b. The two axes come apart, again, and in the recorded direction.** At signal = 0 the
VQE's objective gap is 3.2e−06 — essentially *at* the objective's optimum — and its returned
structure is 2.983 Å against a space best of **0.969 Å** on these nine targets. Optimising
this objective perfectly would move the answer by **0.087 Å**. **The optimisation axis is
saturated and the structural axis is not, so nothing on the optimiser side of the problem
can matter.** (§8's table is the nine n=9 targets only, so its ρ = +0.176 at signal = 0 and
its 0.969 Å space best differ from §5's nineteen-target +0.264 and 1.030 Å.)

---

## 9. THE SCHEDULE CONTROL — a partial answer, tiered HYPOTHESIS

`python -m s16.qphase sched`. A VQE samples a *sequence* of distributions, broad early and
narrow late; a single Boltzmann tilt does not. `tilt_anneal` gives the classical reweighting
exactly that schedule over exactly the VQE's own untrained support — a geometric temperature
ladder ending at the temperature whose expected distinct count matches the VQE's — and
nothing else. 19 targets × 3 seeds × 4 rungs.

| signal | ρ | VQE | untrained | `tilt_samples` | `tilt_exact` | `tilt_anneal` | anneal |
|---|---|---|---|---|---|---|---|
| *argmin* 0 | +0.264 | 2.782 | 3.091 (−0.309) | 3.128 (−0.346\*) | 3.012 (−0.230) | 3.018 (−0.237) | **2.752** (+0.030) |
| *argmin* 0.6 | +0.879 | 1.401 | 1.834 (−0.432\*) | 1.826 (−0.424\*) | 1.683 (−0.281\*) | 1.740 (−0.339\*) | **1.383** (+0.018) |
| *argmin* 1.0 | +1.000 | 1.216 | 1.469 (−0.252\*) | 1.488 (−0.272\*) | 1.339 (−0.123\*) | 1.400 (−0.184\*) | **1.073** (+0.143\*) |
| *rand5* 0.6 | +0.879 | 2.665 | 3.439 (−0.773\*) | 2.862 (−0.196\*) | **2.767 (−0.102 NULL)** | 3.277 (−0.611\*) | **2.548** (+0.117\*) |
| *rand5* 1.0 | +1.000 | 2.677 | 3.439 (−0.762\*) | 2.884 (−0.207\*) | **2.805 (−0.128 NULL)** | 3.282 (−0.605\*) | **2.601** (+0.076) |

**The annealed independence sampler does NOT reproduce the VQE**, so "the VQE's advantage over
a fixed tilt is simply its schedule" is **not demonstrated**. **But the arm is
mis-calibrated and I say so**: only its *final* temperature is diversity-matched, so it
spends most of its budget hot and ends at 4,985 distinct against the VQE's 4,449 — it
under-concentrates by construction and is a lower bound, not a fair schedule control. The
properly diversity-matched **fixed**-temperature `tilt_exact` is the stronger classical
thermostat and it **ties the VQE on the ensemble readout** (−0.102 and −0.128, both NULL).
**Tier: HYPOTHESIS.** A correctly schedule-and-diversity-matched annealed reweighting is the
obvious next experiment and it is cheap.

---

## 10. WHAT I REFUTED, INCLUDING MY OWN EXPECTATIONS

1. **REFUTED — the boundary is in ρ.** It is not. The two quality parameterisations give
   ρ\*(untrained) = +0.355 and ≤ +0.140 for the same control; a cell-level fit on ρ has
   R² = 0.014 and puts the crossing outside [−1, 1]. My own pre-experiment expectation was
   that a ρ\* would exist and would be comparable across families. It is not.
2. **REFUTED — "better objective ⇒ bigger quantum payoff."** The VQE's argmin gain over its
   untrained start *peaks* at ρ ≈ 0.36 and **falls** monotonically to a perfect objective,
   because best-of-N improves faster than the VQE does.
3. **REFUTED — the CVaR gradient-baseline defect helps.** Null at the target unit on the
   argmin readout at both quality levels; significantly **worse** on the ensemble readout at
   a good objective (+0.314 [+0.093,+0.532]).
4. **REFUTED — the defect is a de-facto step-size reduction.** A 2.6× smaller gradient
   produces a **7% smaller** parameter displacement, because Adam is invariant to a uniform
   rescaling of the gradient. The mechanism recorded in the ledger is unavailable.
5. **REFUTED (in scope) — "less optimisation is better."** Every lr rung says the opposite
   here: cutting the learning rate 16× costs +0.156 Å (argmin) and +0.676 Å (ensemble).
6. **REFUTED — "the one place a VQE wins is an ensemble consumed without a ranker."**
   It survives against the untrained circuit and against an entropy-matched sample
   reweighting; it does **not** survive against a diversity-matched classical tilt (NULL at
   8/10 rungs), and it loses to greedy at every rung and to annealing at most.
7. **CONFIRMED and re-scoped — α is a diversity dial.** But its *strength* as a dial is a
   property of the optimisation budget (1.35× here vs the recorded 150×), and a classical
   thermostat reproduces its whole effect on the ensemble readouts with Pearson +0.93 to
   +0.99 while moving them 1.6–1.7× further.
8. **CONFIRMED — the Sprint 15 mechanism, on a new instrument.** The per-target VQE-minus-
   control difference correlates **+0.519** with the true-RMSD percentile of the objective's
   certified argmin, against Sprint 15's +0.52…+0.69, and **+0.867** with the direct form
   `RMSD(argmin) − RMSD(control)`.

---

## 11. HONEST LIMITATIONS

* **k = 4 and n ≤ 10.** The space's best structure averages **1.030 Å** and the best real
  objective's optimum **2.661 Å**. Nothing here tests anything below the discrete
  instrument's own floor, and none of it speaks to the 126-target pipeline's 3.204 Å.
* **One ansatz (`mps2f`), one budget (8,192), one shot count (512), three seeds.** The
  budget convention is varied only in the classical arms (2,048 / 8,192 / 676,055). Sprint 15
  showed a budget convention can flip a comparison, and the ¼-budget annealing row is the
  place that would bite: it is *null*, and a still smaller classical budget was not run.
* **`tilt_samples` is entropy-matched, not diversity-matched**, and the two conventions differ
  by 1.6× in distinct count here. Every comparison against it is confounded in the VQE's
  favour and is reported with the diversity-matched `tilt_exact` beside it.
* **`tilt_anneal` is mis-calibrated** (§9) and answers its question only as a lower bound.
* **The blend family's ρ and its argmin location are confounded** — argmin percentile
  collapses 0.231 → 0.006 between signal 0 and 0.1. The noise family breaks that confound,
  which is why both are run; but the blend ladder alone would have supported a ρ\* claim that
  the noise ladder refutes.
* **The nine n=9 and ten n=10 targets are different instruments** and are pooled. Per-fold
  means are stored in every `I.paired` record; the pooled conclusions do not turn on it, but
  the split was not analysed separately here.
* **Contention.** The box ran three other agents' jobs concurrently for most of this
  workstream; CPU sat at 88–99 %. Every gate is recorded in the logs and one blend target
  (`2MD2`) waited on the CPU gate. Wall-clock costs quoted anywhere should be read as
  contended.
* **Sprint 15 disagreements are instrument differences, not contradictions.** Three results
  here differ in sign or magnitude from the record (the reweighting control's cost ranking,
  α's diversity range, the baseline defect). In each case the two experiments differ in
  register size, estimator and iteration count simultaneously, so neither supersedes the
  other; both scopes are stated in place.

---

## 12. WHAT THIS MEANS, AND WHAT REMAINS OPEN

**The positive, generalisable statement this workstream set out to produce.** The intended
form was *"a variational sampler is the right tool when the objective's rank correlation
exceeds ρ\*."* **That statement is false**, and the measurement that replaces it is sharper:

> **A variational sampler pays exactly when the objective's own optimum is better than what a
> fixed sampler already returns from that objective, and it captures ~90 % of that difference
> on a structured landscape. Because a classical 1-opt with restarts captures 100 % of it at
> the same budget, and simulated annealing captures it at a quarter of the budget, there is
> no objective quality at which the variational arm is the right tool for an argmin readout.
> Its only distinctive product is a distribution — a 3.6× near-native mass enrichment at the
> real objective — and a diversity-matched classical Boltzmann tilt of the same untrained
> circuit supplies an ensemble of equal value.**

**Where this problem sits.** The deployed native-free structural objective has ρ = +0.264,
in-tail ρ = +0.083, and **its certified global optimum is 2.661 Å against a space best of
1.030 Å**. Legacy's optimum is 4.085 Å (worse than the 4.003 Å space mean); AMBER's is
3.931 Å at the 55th percentile; the distogram consumed generatively 3.123 Å and selectively
3.077 Å. **Every real objective in this project has an optimum in the 2.7–4.1 Å band, and a
classical 1-opt reaches every one of them exactly at 8,192 evaluations.** The binding
constraint is not search, not concentration and not the sampler; it is where the objective's
minimum is, which is the same statement the standing ledger records as *"search saturates;
discrimination binds"* — now measured with the optimum's location as the explicit variable.

**Open.**

1. A correctly schedule- **and** diversity-matched annealed classical reweighting (§9).
2. The boundary was measured only for one ansatz and one budget. The prediction that follows
   from §0.4 — that the VQE-minus-classical gap is flat in objective quality for *any*
   ansatz — is untested.
3. Whether `RMSD(argmin E) − RMSD(control)` is a **native-free** quantity in any usable form.
   It is not: it reads the native twice. But `RMSD(argmin) − RMSD(control)` has a native-free
   *upper bound* through the pool's own diversity, and that has never been priced.
4. Nothing here escapes k = 4. The discrete instrument's own floor is 1.030 Å.

---

## 13. REPRODUCTION

```
python -m s16.qphase verify        # every pre-condition, before any claim
python -m s16.qphase ladder        # T1a  blend ladder, 19 targets x 10 rungs x 3 seeds
python -m s16.qphase ladder_noise  # T1b  noise ladder, 19 targets x  9 rungs x 3 seeds
python -m s16.qphase realarms      # T2   the five REAL objectives through the same arms
python -m s16.qphase real          # T2   rho profiles incl. AMBER on the binding stratum
python -m s16.qphase alpha         # T4   six alphas + two entropy/diversity-matched tilts
python -m s16.qphase lr            # T3   the lr ladder + the tail-baseline defect
python -m s16.qphase doing         # T5   exact 2^18 distributions, four axes
python -m s16.qphase sched         # the schedule control
python -m s16.qphase_report report        # every table above, from the checkpoints
python -m s16.qphase_report report sched  # section 9, the schedule control
python -m s16.qphase_report report curves # sections 6b and 5b, regenerated not typed
```

Artefacts: `s16/results/qphase_{verify,ladder_blend,ladder_noise,realarms,real,alpha,lr,doing,sched}.json`.
Every mode checkpoints after each target; `report` regenerates every table from the artefact,
so no number in this document is typed by hand.

### Leakage audit

`Enum.rmsd`, `nat_ca`, `u_truth` and every RMSD are ORACLE quantities. They enter (i) post-hoc
scoring, (ii) the two deliberately-labelled quality knobs, and nowhere else. The five real
objectives in §5 are native-free: `hamil` is the retrieval prior plus the leave-fold-out
distogram; `E_ls_pred` and `S14_disto_bayes` are two consumptions of the same leave-fold-out
distogram; `legacy` and `prior` read no native. No native-derived quantity enters any arm's
parameters, any stopping rule, any threshold, or any selection. The 60-target benchmark was
not read.

