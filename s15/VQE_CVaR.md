# VQE AND CVaR

The quantum component, treated as a scientific variable rather than a decoration. This document
states what the variational algorithm does on this problem, what CVaR does to the distribution it
produces, what is correct in the implementation, what is defective, and what the whole apparatus is
worth against the one control that matters.

Status: **living**. All sections are complete; the replication of §4.2c is in flight.

---

## 1. THE CONTROL THAT DECIDES EVERYTHING

A quantum optimiser must beat **best-of-N sampled from the same circuit at its untrained
parameters**, at a matched hard budget of objective evaluations, paired across targets, with a
confidence interval excluding zero.

This is a stronger control than the field uses. The nearest published comparison samples **uniformly
at random**, which is strictly weaker because it does not share the ansatz's support — a circuit's
untrained distribution is already structured, and beating uniform noise demonstrates only that the
structure exists, not that training added anything. The literature survey found **no paper anywhere**
that runs the untrained-circuit control.

**Beating uniform random proves nothing on this instrument.** A zero-information constant α-helix
beats the random control here. Any arm that reports only a win over random has reported nothing.

---

## 2. WHAT IS CORRECT IN THE IMPLEMENTATION

Verified before any claim was drawn from it, because a negative result obtained with broken
machinery is worthless.

| quantity | check | result |
|---|---|---|
| CVaR | against the Rockafellar–Uryasev definition | max error **4.2e-14** over 400 cases |
| `dCVaR/dp` | finite differences | median error **0.0** |
| the analytic torsion gradient | central differences | **cosine 1.000000000000000** |
| the metric | classical Fisher vs Fubini–Study | equal to a factor of 4, to **3.3e-16** |

---

## 3. THE VARIATIONAL GEOMETRY

### 3.1 QNG and classical natural gradient are the same algorithm here

For a real-amplitude circuit measured in the computational basis, the **classical** Fisher
information of the measured distribution equals **four times** the Fubini–Study metric, to 3.3e-16.
The factor is absorbed into the learning rate.

This is a structural result rather than an empirical one, and it **removes an entire arm of the
design space**: there is no separate quantum-natural-gradient method to evaluate on this class of
circuit and this measurement basis. Any paper proposing QNG here is proposing natural gradient.

### 3.2 The published condition numbers were a seed-averaging artefact

Sprint 14 recorded metric condition numbers of 2.2–3.7. Those are statistics of a **seed-averaged**
metric, not the condition number of the metric at any point. Recomputed per point, the real spread
is enormous: for a chain at depth L = 3, **median 735, maximum 1.5e6**.

The recorded rank column is the same artefact. The `block` ansatz at depth ≥ 2 is **exactly
rank-deficient**, which the seed average hides completely.

### 3.3 The ill-conditioning is real, and the gradient does not point into it

The null-space share of the gradient is **0 to 1e-30** — a theorem, confirmed numerically. The
bottom eigenvalue decile carries **0.0016–0.0043** of the gradient's squared norm, where a uniformly
distributed gradient would carry 0.10.

This upgrades a Sprint 14 *hypothesis* ("the badly-conditioned directions are the reducible ones") to
a **proof** for the exactly-null directions and a **measurement** for the rest.

It also matters for interpreting the negative results below: **ill-conditioning is eliminated as an
explanation for the optimiser's failure.** A competing account has been ruled out rather than
ignored.

### 3.4 It is not a barren plateau, and the scaling axis was misidentified

`Var ~ 2^(−0.18 to −0.36 n)`, with **every confidence interval excluding −1.0**. Describing this as a
barren plateau would be wrong.

The real axes are **depth** (`L^−2.6` at n = 16) and **cost locality** (6.05 log₂ at n = 14, against
2.8 log₂ for the entire width axis) — which refutes the earlier "what matters is n" in its
operational reading. The 1-local retrieval prior has **zero** width decay (+0.006).

### 3.5 Does QNG help? The budget convention decides — a trap worth the whole section

| convention | result |
|---|---|
| **equal iterations** | QNG beats plain gradient descent **5/5** where the condition number > 5, **0/2** where < 5 |
| **equal hardware cost** | the full quantum geometric tensor wins **0 of 7** cells, with two clean sign flips |

> **Reporting only the equal-iteration convention would have produced a false positive.**

Only the *diagonal* preconditioner survives cost accounting. On the real VQE across 24 cells, QNG is
**indistinguishable from plain gradient descent** where the condition number is 688–2,677 (+0.0016,
CI [−0.010, +0.014], 7 wins to 10) and **loses to plain Adam 17/17** (+0.067, CI excluding zero) — by
up to 27× at equal cost.

The mechanism follows from §3.1 and §3.3: `F = 4g` exactly, so there is no second metric to exploit;
the gradient lies in `range(g)` and avoids its small eigendirections; so the rescaling QNG applies
**along the directions the gradient actually occupies** is only **1.4–3.1×** even at condition number
2,000.

### 3.6 Target conditioning does not reshape the manifold in any information-bearing way

Preparing a state for the retrieval-conditioned torsion prior raises the Fubini–Study condition
number **3.85× to 36×** and cuts metric volume by 1.8–3.7 log units, at **27 wins to 0**. Against an
**entropy-matched scramble of its own prior** — per-residue entropies identical to the bit, target
information destroyed — it is **NULL on all seven statistics, both ansätze, 54 paired cells**, the
scramble reproducing the effect to within 3–7%. The one statistic whose interval excluded zero fails
the null-calibrated concentration check at p = 0.000–0.001.

> **The chain `target information → manifold → metric → trainability` does not exist. The middle link
> is real and information-blind.**

### 3.7 Optimising well and getting a better structure are anti-correlated

Gradient norm falls **313×** while the mode's RMSD moves **0.049 Å**, and only 25% of runs reach a
positive-definite Hessian. Most starkly: **the ansatz that optimises best reaches a minimum in 0 of
18 runs and *improves* structure by 0.46 Å, while the one that does reach a minimum makes structure
0.37 Å worse.** This is the two-axis separation in its sharpest form, and it is why this programme
reports optimisation and accuracy separately and never substitutes one for the other.

---

## 4. WHAT CVaR ACTUALLY DOES

Measured in the most favourable setting a CVaR-VQE can be given: the **exact statevector gradient**,
no shot noise, no SPSA, budget charged only for readout draws, and the control given the same draws.
12-qubit sub-registers of the nine enumerated targets, 200 iterations, 3 seeds, 4,096 draws, paired
over 27 cells per α.

**Positive means VQE is worse than best-of-N from its own untrained initial distribution.**

| readout | α = 1.0 | α = 0.25 | α = 0.05 | α = 0.01 |
|---|---|---|---|---|
| argmin (by objective) | −0.036 null | **+0.251 SIG** | +0.068 null | +0.048 null |
| top-20 mean | −0.177 null | +0.111 null | −0.070 null | −0.091 null |
| **top-20 coordinate average** | +0.068 null | **+0.372 SIG** | +0.197 null | +0.177 null |
| **top-75 coordinate average** | +0.147 null | **+0.392 SIG** | **+0.276 SIG** | +0.164 null |
| drawn-set MEAN | **−0.847 SIG** | **−0.372 SIG** | **−0.191 SIG** | −0.066 null |
| drawn-set BEST | **+1.092 SIG** | **+0.345 SIG** | **+0.155 SIG** | +0.087 null |

### 4.1 What concentration is, in one pair of rows

At α = 1 the VQE's drawn set has a mean RMSD **0.847 Å better** than the control's and a best member
**1.092 Å worse** — two significant effects of opposite sign and near-equal magnitude, both shrinking
monotonically to zero as α falls.

**Concentration buys the mean and pays the best, and α scales both.** This is the Sprint 14 mechanism
measured directly rather than inferred.

### 4.2 On every readout that survives to a structure, VQE loses

Both coordinate-average readouts — which are what the pipeline actually consumes — are significantly
**worse** than the untrained control at α = 0.25, and top-75 is worse at α = 0.05 too. **Nothing in
the table is a significant win for VQE on a structural readout.**

### 4.2b CVaR has no role as a constraint either

A plain expectation penalty drives clash-violating probability mass to **exactly zero**; CVaR leaves
**2.8–5.6%** and returns worse structures. Both lose to best-of-N by **+0.41 to +0.70 Å**. Of the
three roles CVaR was given in the architecture families — objective, tail-shaper, constraint — two
are dead and only the diversity dial survives.

### 4.2c THE ONE POSITIVE RESULT, AND WHY IT IS NOT YET QUOTED

Consumed as an **ensemble with no ranker anywhere**, VQE *appeared* to beat best-of-N by
**0.36–0.57 Å** with intervals excluding zero. **The replication workstream's audit of the same data
withdraws that reading.** The 27 "cells" are **9 targets × 3 seeds**, so the unit of analysis
overstates n threefold; at the **target level (n = 9) the α = 1 cells are NULL** (−0.567 [−1.234,
+0.085] and −0.437 [−1.056, +0.160]); the surviving α = 0.05 cell **FAILS** the null-calibrated
concentration check, carried by two or three of nine targets; the readout subsamples by *sorted*
index rather than at random, biting unequally across arms; and the VQE was charged **819,200
objective evaluations against the control's 2,048** — a **400× budget advantage**. **The programme has
no positive quantum result.**

**A mechanism exists.** The programme has measured repeatedly that *nothing ranks within the pool*,
and in the generative cascade selecting the objective's argmin **costs +0.751 Å**. If every ranker is
useless, a method whose value lies in the *shape of its distribution* rather than in any member it
can identify would appear exactly here — invisible on argmin readouts, visible only unranked.

**A reason to disbelieve exists.** At α = 1 the ensemble holds 2.7 distinct structures out of 4,096
draws, and the win is **non-monotone in α** (present at 1 and 0.05, absent at 0.25). A non-monotone
effect with no mechanism is usually noise.

**Until it replicates it is not a result.** It is the only positive quantum claim in the project's
history, which makes it the one most likely to be wrong.

### 4.3 α is a diversity dial, not an accuracy dial

| α | distinct configurations among 4,096 draws | control | final entropy (bits) | max p |
|---|---|---|---|---|
| 1.00 | **2.7** (median **2**) | 510.3 | 0.73 | 0.702 |
| 0.25 | 149.4 | 510.3 | 3.93 | 0.319 |
| 0.05 | 314.9 | 510.3 | 5.79 | 0.155 |
| 0.01 | 408.0 | 510.3 | 6.58 | 0.101 |

The monotone entropy column is the cleanest reading of what α does. Shaping ensemble spread is a
legitimate and honest role for CVaR in a conditional-ensemble architecture. **It is not the role the
CVaR-VQE literature claims for it.**

### 4.4 A trap that invalidates a law this project uses constantly

The standing law is that the terminal operator consumes the set **mean**
(`d_out = 1.16·d_set_mean + 0.04·d_set_best`, R² = 0.89). Read naively, the `drawn-set MEAN` row says
CVaR-VQE wins by 0.85 Å. **The actual coordinate average says it loses.**

At α = 1 the VQE returns **two distinct structures out of 4,096 draws**. Its "top-20 set" is twenty
copies of one structure, so its set *mean* is that structure's RMSD and there is nothing to average.

> **The set-mean law was fitted on sets with real diversity and is out of domain on a collapsed set.
> Set mean is a valid proxy for the terminal operator only at fixed set diversity; a method that
> concentrates improves the proxy and degrades the actual output.**

Every arm in this programme that invokes the set-mean law now reports the diversity of its set
alongside its RMSD. The generative cascade's ensemble passes the check — mean pairwise RMSD 0.527 Å,
with only 17 of 126 targets below a 0.25 Å spread — so its aggregation result is inside the law's
domain.

---

## 5. THE CVaR DEFECT TRIPLE

The estimator as deployed carries three defects. They are **findings**, not bugs to be quietly
fixed, and they are what remains of this programme's CVaR novelty after the literature survey found
that in-loop CVaR for peptide folding is already published.

The most consequential: **the gradient baseline is centred on the tail only**, giving a cosine of
**−0.023** with the true gradient — very nearly orthogonal to the real one, and slightly anti-aligned.

**And measurement inverts the expected consequence: the defect HELPS.** The biased estimator gives
**−0.266** where the corrected one gives **−0.145** on the ORACLE objective. The mechanism is measured
rather than guessed — the biased gradient's norm is **2.4× smaller**, so it acts as a de-facto
step-size reduction. **A weaker optimiser is a better sampler**, because it concentrates less, and
concentration is what costs the distinct samples a best-of-N readout is paid in.

So the defect is not a liability weakening the negatives above; it is a mechanism consistent with all
of them. A third recorded defect is **inert**: `dead_start` = 0.000 across all 19 × 3 × 3 × 8 runs.
§4's headline was nonetheless obtained with the **exact** statevector gradient, so no reading of it
depends on the estimator at all.

---

## 6. DOES THE NEGATIVE SURVIVE A CHANGE OF OBJECTIVE CLASS? — ANSWERED: THE QUESTION WAS THE WRONG ONE

Every objective in §4 and in Sprint 14 was **selective**. The restraint objective is **generative**,
and its ordering skill is measured to be real without any search: ρ(objective, RMSD) = **+0.562**,
positive on 88.1% of targets, its argmin over a 500-member pool beating a random member by
**−0.949 Å [−1.147, −0.753]**.

The experiment was run on **19 enumerated targets** with a budget ladder, a random-tail null, and —
decisively — **Sprint 14's own selective objective included as an arm**, so that the two classes could
be compared on one instrument with one control. The instrument was verified by reproducing Sprint
14's flagship table to three decimals.

### 6.1 The answer is not a sign

| arm, budget 8,192 | paired difference | wins |
|---|---|---|
| `E_ml_pred` (generative) | **−0.202** | 13/19 |
| `E_ls_pool` (generative) | **−0.184** | 14/19 |
| **`S14_disto_bayes` — the SELECTIVE control** | **−0.192, the strongest cell** | **15/19** |

| budget | cells significant, of 6 |
|---|---|
| 2,048 | 1 |
| 8,192 | 3 |
| **32,768** | **0** |

**VQE beats the untrained control at 8,192 — and the win is objective-independent.** Its strongest
cell is the objective with the *worst* in-tail ordering measured anywhere in the programme. Whatever
it is doing, it is not exploiting the generative objective. And the effect is **gone at the largest
budget**.

**Against a classical greedy search at matched budget: 0 wins in 16, losing 3.**

### 6.2 The mechanism reproduces, and is not a selective-class pathology

`tail − bulk` ordering skill is **negative for every objective** — all four generative ones, the
selective control, and the **ORACLE**. The generative objectives do have better **in-tail** ρ (+0.127,
z = 13.1, against −0.042), but it never becomes positive relative to the bulk; and **globally** the
two classes are close (ρ = +0.523 for `E_combined` against +0.485 for the selective control), so there
is no global-ordering advantage worth the name either.

> **The binding quantity across every arm is a selection gap of 1.47–2.14 Å — identical for VQE and
> control, flat across a 16× budget range, collapsing to 0.29–0.44 Å under ORACLE restraints. It is
> restraint error, not the sampler.**

### 6.3 What this means for how the result should be stated

The defensible claim is **not** "VQE loses". It is:

> Whether a variational optimiser beats its own untrained initialisation depends on the **budget**,
> the **readout** (argmin and diversity-respecting aggregation go opposite ways), and **n** — three of
> the workstream's own reads reversed between n = 1, n = 9 and n = 19. It never beats a classical
> greedy search. And the quantity that dominates every arm is identical in the quantum and control
> arms, so the sampler is not what is being measured.

That is the third reporting trap in this document, after the QNG budget convention (§3.5) and the
set-mean-on-a-collapsed-set trap (§4.4). All three share a shape: **a single cell, reported without
the axis it depends on, supports a conclusion the full grid does not.**

### 6.4 A second positive cell, and its stated gap

The trained state enriches sub-2 Å probability mass by **5.2–6.7×** over its own untrained
initialisation (**18.6×** under ORACLE restraints), with an internal control confirming it is real —
the weakest objective gives 0.98×. **It is worth ≤ 0.2 Å**, because an argmin readout ignores
probability mass. That prices Family C precisely: the quantity VQE genuinely improves is the one the
terminal operator does not consume.

**The gap, named by the workstream itself:** no classical ensemble sampler was run against this
enrichment, so it measures *what the VQE does*, not *what only a VQE can do*. Under test.

---

## 6.5 THE CLOSING CONTROL: every effect dissolves classically

The apparent positive of §4.2c was replicated at power — **19 targets, 8 seeds, 1,672 new cells**,
after reproducing all 108 original cells **bit-for-bit (max difference 0.000e+00)**.

**The effect is real.** Target-level, concentration DIFFUSE: **−0.385 [−0.748, −0.035]** at α = 1 and
**−0.295 [−0.495, −0.069]** at α = 0.25, **stronger on the ten targets never used to find it** than
on the original nine. A 10-point α map is smooth and monotone, significant at 9 of 10 — so the
non-monotonicity previously cited as a reason to disbelieve was two noisy three-seed cells.

**And it is not quantum.**

> **A classical Boltzmann reweighting of the same untrained circuit by the same objective at matched
> entropy, at 1/200th of the budget, beats the CVaR-VQE by +0.25 to +0.48 Å on every readout at every
> α** (W/L up to 1/18). Annealing at 1/400th matches or beats it. **Annealing finds the certified
> global optimum in 100% of cells at 2,048 evaluations; the VQE puts 5.9% of its mass there using
> 819,200.**

**What it is.** Pure location — `d_coordavg = 0.94 · d_set_mean`, R² = 0.75, no residual shape term.
The per-target win correlates **+0.52 / +0.64 / +0.69** with where the objective's certified argmin
sits in the ORACLE RMSD distribution: the method is rewarded where the *objective* is right and loses
by +1.25 Å where it is not. Matched diversity removes the remainder, to at or below the 0.08 Å floor.

**The enrichment collapses identically.** It replicates at 5.20× — and classical annealing at the
same budget gives 5.34× (null against it; significantly better on the ORACLE objective). The headline
is a ratio of means; the **median per-target ratio is 1.13–2.66×, with 5 of 9 targets showing any
enrichment**.

**A coordinator prediction, refuted.** That the win should track distinct-configuration count:
within-α **ρ = +0.004** over 1,216 pooled cells, aggregate direction *opposite*.

---

## 7. QUBIT ACCOUNTING

With per-residue torsion library size k and chain length n, the nominal register is `n · log₂ k`
qubits. **The live count is `(n − 2) · log₂ k`.**

Four torsions per chain — `phi[0]`, `psi[0]`, `phi[n−1]`, `psi[n−1]` — are **inert for the Cα
trace**. Confirmed by two independent routes: direct perturbation of ±0.7 rad moves the Kabsch RMSD
by 0 to 1.4e-07 Å; and the analytic gradient at those entries is 0, 2.98e-14, 0 and 0 against a
maximum of 3.31e+02 and an interior median of 6.76e+01 — a ratio of 9e-17.

Any qubit count that ignores this **overstates the register by four**.

---

## 8. WHAT THIS SECTION DOES NOT CLAIM

- No quantum advantage is claimed, and none survives. Some arms beat the untrained circuit at some
  budgets, but the effect is objective-independent, absent at the largest budget, **beaten by a
  classical reweighting of the same circuit's own samples at 1/200th of the cost**, and beaten by
  greedy search (0 of 16) and by annealing on both axes.
- The negative is not attributed to barren plateaus (§3.4), to ill-conditioning (§3.3), or to a
  broken CVaR implementation (§2) — each was measured and eliminated.
- The results in §4 are on **12-qubit sub-registers of nine enumerated targets at k = 4**. They are
  conclusions about a discrete configuration space of that size, not about the continuous problem.
- In-loop CVaR for peptide folding is **already published** (QuPepFold, PLoS One 21(2):e0342012).
  This programme's CVaR contribution is the defect triple and the diversity measurement, not the
  application.
