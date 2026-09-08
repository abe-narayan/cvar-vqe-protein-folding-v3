# SPRINT 16 / FLAGSHIP — JACOBIAN QUIET-SUBSPACE STEERING

**Module** `s16/steer.py` · **artefact** `s16/results/steer.json` · **log** `s16/results/steer.log`
**Pre-registration** `s16/PREREG_steer.md`, written while the n = 126 run was in flight
**n = 126**, target as the unit, paired fold-clustered bootstrap, step size chosen leave-fold-out

---

## The verdict in one line

**The hypothesis is refuted, and the reason is not the one the pre-registration anticipated.**
Steering fails not because the native-free direction is too weak, and not because the loud/quiet
split is wrong, but because **RMSD is a strongly convex function of torsion displacement toward the
truth**: a partial step in the *exactly correct* direction buys almost nothing.

> **§3's causal explanation was corrected on 2026-09-06 by the VERIFY workstream, and the correction
> matters.** The convexity is real and every number in §3 reproduces, but it is **a large-displacement
> regime, not a property of torsion space** — a zero-information control interpolating between two
> arbitrary retrieval windows is *more* non-monotonic than the fit. The displacement is large because
> **the torsion-space fit is 87% of a uniformly random torsion assignment**. See §3 and ledger L21.

---

## 1. What was measured

| quantity | value |
|---|---|
| plain torsion-space fit | **3.647 Å** |
| incumbent (BLOSUM62 → top-75 → coordinate average → projection → AMBER) | 3.204 Å |
| native-free \|cos(u, e)\| on this fit | **0.318 [0.280, 0.359]**, median 0.277 |
| top-4 singular share of the superposed-CA Jacobian | 0.893 |
| first-order prediction ‖J e‖/√n vs realised RMSD | **9.879 Å vs 3.647 Å — ratio 2.709** |

### The decomposition (ORACLE — it reads the native error)

| rank r | ‖e_loud‖/‖e‖ | RMSD share of loud | \|cos(P_r u, e_loud)\| |
|---|---|---|---|
| 2 | 0.197 | 0.716 | **0.712** |
| 4 | 0.303 | **0.875** | **0.536** |
| 6 | 0.402 | 0.943 | 0.457 |
| 10 | 0.552 | 0.989 | 0.438 |
| 13 | 0.655 | 0.996 | 0.417 |

**The projection did exactly what it was designed to do.** Restricting the native-free direction to
the loud subspace raises its alignment with the loud part of the true error from **0.318 to 0.536**
at rank 4 and to **0.712** at rank 2. The mechanism's *first* step works.

---

## 2. The steering ladder

Leave-fold-out step, five folds. `s = 0` is on the ladder, so a training fold can always choose
"do nothing"; a positive (worse) held-out effect is therefore a real generalisation failure, not an
artefact of the selection rule.

| arm | mean | median | vs plain fit | LFO step per fold |
|---|---|---|---|---|
| `full` (native-free) | 3.647 | 3.412 | **+0.000 [+0.000, +0.000]** | 0, 0, 0, 0, 0 |
| `consensus` (native-free) | 3.654 | 3.442 | +0.007 [−0.026, +0.039] | 0.1, 0.1, 0.1, 0, 0 |
| `loud4` | 3.647 | 3.412 | **+0.000 [+0.000, +0.000]** | 0, 0, 0, 0, 0 |
| `loud10` | 3.647 | 3.412 | **+0.000 [+0.000, +0.000]** | 0, 0, 0, 0, 0 |
| `loud13` | 3.647 | 3.412 | +0.000 [+0.000, +0.000] | 0, 0, 0, 0, 0 |
| `quiet_half` | 3.645 | 3.402 | −0.003 [−0.011, +0.007] | 0.35 × 4, 0.5 |
| `quiet_quarter` | 3.652 | 3.409 | +0.005 [+0.000, +0.011] | 0, 0, 0, 0, 0.7 |
| `ORACLE_quiet_half` | 3.631 | 3.389 | −0.016 [−0.048, +0.018] | 0.5 × 4, 0.7 |
| `ORACLE_loud4` | 3.619 | 3.445 | **−0.028 [−0.106, +0.051]** | 0.1 × 5 |
| `ORACLE_true` | 0.043 | 0.023 | **−3.604 [−3.895, −3.343]** | 1.0 × 5 |

**Every native-free arm chose "do nothing" on every fold.** The loud projection raised the cosine
and bought nothing.

### The pre-registered rules, and which fired

| pre-registered rule | fired? |
|---|---|
| `loud*` beats `full` with a CI excluding zero → hypothesis supported | **no** |
| `loud*` does not beat `full`, both ≥ 0 → **refuted at the mechanism level** | **YES** |
| `quiet_half` improves with a CI excluding zero → framing is wrong | **no** — see §4 |
| `ORACLE_loud4` worse than plain fit while rank-4 RMSD share > 0.8 → the first-order decomposition may never be quoted as a finite-step budget | **YES** (−0.028, CI includes zero, against a 0.875 first-order share) |

---

## 3. The mechanism: RMSD is convex in torsion displacement

The `ORACLE_true` curve is the finding. It steps the **exact** native error, so at `s = 1` it lands
on the native torsions by construction (0.043 Å). What it does in between is the point:

| step s along the exact error | 0.0 | 0.1 | 0.2 | 0.35 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|---|---|---|
| mean RMSD (Å) | 3.647 | 3.599 | 3.641 | 3.663 | **3.466** | 2.596 | 0.043 |
| fraction of gain, from the means | 0 | 1.3% | 0.2% | −0.4% | **5.0%** | 29.2% | 100% |
| **fraction of gain, per target** (mean [95% CI]) | 0 | 1.6% [−0.4, +3.4] | 0.8% [−3.2, +4.7] | 1.1% [−5.3, +7.0] | **7.7% [−0.0, +14.9]** | 32.3% [+25.2, +38.7] | 100% |
| targets improved / worsened | — | 82/44 | 75/51 | 73/53 | **77/49** | 104/22 | 126/0 |

The per-target aggregation is the correct one and is the one quoted; the ratio-of-means row is shown
only so the two are not confused. Median per-target fraction at half-way: **16.1%**.

**Half-way to the correct answer, in the correct direction, is worth about 8% of the answer, with a
confidence interval that touches zero — and it makes 49 of 126 targets WORSE than not moving at
all.** A linear map would give 50% on every target. The curve is not merely sub-linear; it is
*non-monotonic* over the first third, and on more than a third of the instrument a half-step in the
exactly correct direction is actively harmful.

This is why the loud projection cannot help and why no improvement to the direction estimator would
rescue it. Steering requires that partial motion toward the truth pays partially. In torsion space,
it does not: the payoff is concentrated in the last 30% of the path.

**Why — CORRECTED 2026-09-06 by the VERIFY workstream (ledger L21).**

~~The fit's torsion errors are *compensating* — wrong in a correlated way that partly cancels along
the chain, which is the same property Sprint 15 recorded as **realizability**. Moving part-way
breaks the cancellation without establishing the correct structure, so the intermediate conformers
are worse than either endpoint.~~ **[RETRACTED. This sentence had no control, and the control
refutes it.]**

VERIFY supplied the missing control: geodesic interpolation between **two arbitrary retrieval
windows** — native-free, no fit involved, matched on ‖Δθ‖ — gives **−0.024 at s = 0.5 against the
flagship's +0.014**. *The zero-information control bulges more than the fit does.* So the
non-monotonicity is not a property of the fit's errors being correlated or compensating; it is not a
property of the fit at all.

**What it actually is: a magnitude regime.** The effect tracks ‖Δθ‖ with Spearman **−0.571**, and
**below about 4 rad both the fit and the control are linear** (+0.50 and +0.47 of the gain at
half-step). The flagship's error simply sits far above that: **90° per-angle RMS, which is 87% of a
uniformly random torsion assignment.**

**The corrected reading is more damaging to the programme, not less.** The convexity is not a special
property of torsion space that a better parameterisation might dodge — it is what happens to *any*
large torsion displacement. And the reason the displacement is large is that **the torsion-space
distance-geometry fit is very nearly a random torsion assignment.** The measured curvature ratio,
**2.709**, is the same fact in the Jacobian's language: the error is 2.7× outside the regime where
the linearisation describes the motion.

**The two things this does not change.** The numbers in the table above all reproduce exactly, and
the per-target aggregation was the right choice (VERIFY confirmed both). And **the geodesic
objection is refuted**: `A.wrap` *is* the exact per-angle geodesic on the torus, agreeing with an
independent great-circle construction to 1.8e−15 rad, so the non-monotonicity is not an artefact of
interpolating wrapped angles off the circle.

**A consequence for the whole programme, which follows from the corrected reading.** If the fit is
87% of a random torsion assignment, then every downstream stage is repairing near-noise, and the
programme's accuracy is set by how much structure the *retrieval and averaging* stages recover, not
by anything the torsion fit contributes. That is consistent with the incumbent (3.204 A) beating the
fit (3.647 A) by a wide margin, and with the pool coordinate average (3.048 A) beating both.

**This also disposes of the Sprint 15 observation that motivated the sprint.** The ORACLE rotation
into the loud subspace was computed in the *tangent space at θ_fit*. The tangent space is a faithful
description of a neighbourhood the prediction error does not live in. The five exactly-null
directions, the participation ratio of 3.47, and the loud/quiet split are all real and all correctly
measured — they are simply **local** facts being asked a **finite-step** question.

**Connection to a standing project law.** This is the same geometry as the recorded result that
*coordinate averaging beats torsion averaging by 1.024 Å*. Both say: **torsion space is not a space
in which you may move part-way.** The flagship is an independent, quantitative confirmation of that
law from the opposite direction.

---

## 4. The n = 8 smoke was wrong about everything that mattered

Recorded here because the programme's most expensive recurring error is reading a small-n probe.

| arm | n = 8 smoke | n = 126 | verdict on the smoke |
|---|---|---|---|
| `full` | +0.382 [−0.043, +0.898] | +0.000 | **wrong** |
| `loud4` | +0.272 [−0.008, +0.659] | +0.000 | **wrong** |
| `loud10` | +0.495 [−0.010, +1.240] | +0.000 | **wrong** |
| `quiet_half` | **−0.106 [−0.185, −0.031]** (CI excluded zero) | **−0.003 [−0.011, +0.007]** | **wrong — a false positive with a CI that excluded zero** |
| `ORACLE_loud4` | +0.320 | −0.028 | wrong in sign |
| curvature ratio | 2.269 | 2.709 | the one thing that held |

The quiet-control anomaly that drove the entire follow-up plan **evaporated at n = 126**. This is the
sixth time in the programme that n ≤ 8 has reversed a conclusion, and the first time one did so with
a confidence interval that excluded zero at the small n. **A CI at n = 8 is not protection.**

---

## 5. What survives, and what may still be quoted

**Survives.**
- The loud/quiet spectral split of the superposed-CA Jacobian, as a **local** description. Top-4
  singular share 0.893; rank-4 carries 0.875 of the RMSD to first order.
- The native-free direction's alignment with the true error, **0.318 [0.280, 0.359]** on this fit.
  Note this is *lower* than the 0.390 Sprint 15 recorded on its own fit; the AUDIT workstream
  independently reproduced 0.390 exactly (`cos_d_fit_pool` mean 0.3902 [0.363, 0.426], n = 126), so
  the gap is a difference between the two fits, not a contradiction, and both are reported.
- The projection raises alignment as designed (0.318 → 0.536 at rank 4).

**Withdrawn.**
- Any reading of the rank-r RMSD share as a budget for a finite step. It is a first-order quantity
  and the error is 2.7× outside the linear regime.
- The n = 8 quiet-control anomaly and everything built on it.

**Never was a result.** Every `ORACLE_*` arm and the whole of the decomposition table. They are
ceilings and diagnostics.
