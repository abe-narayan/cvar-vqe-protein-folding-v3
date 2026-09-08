# ARCHITECTURE

The four conceptual families declared at the start of Sprint 15, their distinguishing experiments,
and their status as measured. Two families count as distinct only if an experiment exists whose
predicted outcome differs between them; that test is applied below and one family fails it in its
original form.

Status legend: **ALIVE** (under measurement) / **REDIRECTED** (falsified as declared, reformulated
on evidence) / **RESERVE** (deliberately unfilled) / **CLOSED**.

---

## THE COMMON SUBSTRATE

All four families share one departure from every previous sprint. The incumbent pipeline **selects**
among retrieved candidates; these families **solve** for a conformation.

- **Representation.** Continuous backbone torsions `(φ, ψ)` under ideal bond geometry, built with
  the bit-exact production builder, so solved structures lie on the manifold the pipeline emits.
- **Gradient.** Exact O(n) reverse-mode, verified at **cosine 1.000000000000000** against central
  differences. Four torsions per chain are inert for the Cα trace, so the live parameter count is
  `2n − 4` and the live qubit count is `(n − 2)·log₂k`.
- **Ceiling.** From the true distance matrix this machinery reaches **0.611 Å mean, 0.055 Å median,
  86% under 2 Å** on all 126 targets (ORACLE). The representation, the ideal-geometry constraint and
  the search are **definitively not the barrier**.
- **The barrier.** Restraint accuracy, and nothing else. Both the ORACLE and the predicted arms use
  identical optimisers, starts, selection rules and geometry, so the 3.03 Å between 0.611 and 3.644
  is entirely restraint error.

**The reframing this substrate bought (K5).** The project's standing belief that perfect distance
knowledge caps the answer at ~1.95–2.0 Å was a measurement of *selection through a finite library*,
not of the distance channel. Removing the library moves that ceiling to **0.611 Å**. Every prior
estimate of the value of better distance prediction was computed against the wrong floor and
understated the floor by 3.19× and the addressable headroom by 2.07×.

**The honest counterweight.** On the full instrument the raw predicted generative arm is **3.644 Å
and loses to the incumbent by +0.440 [+0.290, +0.592]**. The headroom is real; the means to reach it
is not demonstrated. The 0.611 Å is never quoted without the 3.644 Å beside it.

---

## FAMILY A — DISTANCE-RESTRAINED TORSION VQE

**Status: MEASURED. The distinguishing question was answered, and the answer is that the question was
the wrong one.**

**Hypothesis.** Predicted distance restraints define a many-body diagonal Hamiltonian over discrete
torsion-library states, and VQE with CVaR can sample restraint-consistent conformations. The objective
is **generative** rather than **selective**, so the tail-ordering pathology that defeated every earlier
search arm should not arise.

**The experiment.** One variable changed; the entire earlier apparatus reused unmodified — same index
algebra, same MPS ansatz, same hard evaluation `Counter`, same classical searches, same two-axis
reporting. Extended to **19 enumerated targets**, with a budget ladder, a random-tail null, and —
decisively — **the previous sprint's own selective objective included as an arm**, so the two classes
could be compared on one instrument with one control. The instrument was verified by reproducing that
sprint's flagship table to three decimals.

**The result.** At a budget of 8,192, VQE **beats** best-of-N from its own untrained circuit:
`E_ml_pred` −0.202 (13/19), `E_ls_pool` −0.184 (14/19). But:

| observation | consequence |
|---|---|
| the **strongest** cell is the *selective* control (−0.192, 15/19) — the objective with the **worst** measured in-tail ordering | the gain is **objective-independent**; it is not exploiting the generative objective at all |
| significant cells: 1/6 at 2,048, 3/6 at 8,192, **0/6 at 32,768** | the effect is **absent at the largest budget** |
| against classical greedy at matched budget | **0 wins in 16**, losing 3 |
| `tail − bulk` ordering skill | **negative for every objective**, including all four generative ones **and the ORACLE** |

> **The selection gap — 1.47–2.14 Å — is identical for the quantum and control arms, flat across a 16×
> budget range, and collapses to 0.29–0.44 Å under ORACLE restraints. It is restraint error, not the
> sampler.**

**So Family A's premise is refuted on its own terms.** Changing the objective class does not change
what the optimiser is limited by, because the limit was never the objective class. What the family
*did* produce is a reporting result: **whether a variational optimiser beats its own initialisation
depends on budget, readout and n**, and a single-cell report supports any of three conclusions.

**One measured positive, correctly priced.** The trained state enriches sub-2 Å probability mass
**5.2–6.7×** over its own untrained initialisation (18.6× under ORACLE restraints), with an internal
control — the weakest objective gives 0.98×. **It is worth ≤ 0.2 Å**, because an argmin readout ignores
probability mass. That prices Family C precisely: the quantity VQE genuinely improves is the one the
terminal operator does not consume. Whether a classical ensemble sampler does the same is under test.

**And a stated liability inverted.** The recorded CVaR gradient-baseline defect **helps** — the biased
estimator gives −0.266 where the corrected one gives −0.145, because its gradient norm is **2.4×
smaller**: a de-facto step-size reduction. **A weaker optimiser is a better sampler**, because it
concentrates less, and concentration is what costs the distinct samples a best-of-N readout is paid in.

---

## FAMILY B — CONSTRAINED / AUGMENTED-LAGRANGIAN

**Status: CLOSED.** Falsified twice — as declared, before implementation, by a measurement that cost
almost nothing.

**As declared.** `min E_AMBER subject to C_restraint ≤ ε` — restrained refinement as practised in
NMR structure determination, giving AMBER a first-class role, with a feasibility phase diagram as
its distinguishing experiment (a weighted sum has no analogue of a feasibility boundary).

**The falsification (K7).** The precondition is that the native satisfies the predicted restraint
bands. Measuring `z = |d_native − d̂| / sd` over 8,549 pairs:

| ε | fraction of pairs admitting the native |
|---|---|
| 0.5 | 0.248 |
| 1.0 | **0.458** |
| 2.0 | 0.712 |
| 3.0 | 0.840 |
| median z | 1.336 |
| max z | **7.237** |

At one sigma **fewer than half** the restraints admit the native; at three sigma — already so loose
the restraints barely constrain — one pair in six still excludes it. A uniform-ε hard constraint has
a feasible set that excludes the answer at every ε that meaningfully constrains the search. As
specified, the family is dead.

**Why it is redirected rather than closed.** The failure is **concentrated, not diffuse**. Median z
is 1.336, so a typical restraint is only mildly inconsistent; the damage is done by a minority of
catastrophically wrong pairs. That is an outlier problem — exactly what NMR refinement has handled
for thirty years with robust potentials and violation-aware reweighting, rather than uniform bands.

**The reformulation, now under measurement.** Robust and redescending losses on the standardised
residual (`huber`, `soft_l1`, `cauchy`, `welsch`), explicit violation trimming, and **graduated
non-convexity** — annealing the loss scale from nearly quadratic down to sharp, warm-starting each
stage, so restraints are discarded only once the structure is approximately right.

**And the reformulation is REFUTED on the full instrument.** `CAUCHY_lfo` gives
−0.011 [−0.102, +0.077] and `gnc_cauchy` −0.025 across 126 targets — firing `s15/robust.py`'s own
pre-declared falsification condition, which read: *"If no robust loss beats `squared` on the full
instrument, then the outlier reading of K7 is wrong: the restraint errors are diffuse rather than
concentrated."* They do not, so it is. The 8-target "redescending signature" I wrote up (better
median, three times as many sub-2 Å, worse mean) does not survive.

**The most informative detail is what the robust arms *do* achieve.** They reach far better restraint
residuals — median |z| **0.364** against 0.469 — and **convert none of it into accuracy**. Satisfying
the restraints better does not produce a better structure, which is the realizability result arriving
from another direction.

So Family B is **closed**, not merely redirected: the uniform-ε form is falsified by feasibility, and
the robust form is falsified by measurement.

**Why recalibration is not the fix, and this is worth stating because it looks like one.** The
distogram's `sd` is 2.633× over-confident. But the standardised residual's spread is nearly constant
across separation (2.19, 2.98, 2.70, 2.70, 2.67), so the `sd` has the error's **shape** right and is
wrong by a near-uniform factor — and a weighted least-squares argmin is **invariant under uniform
rescaling of all weights**. Recalibration alone cannot change a single structure. Only the loss can
distinguish a mildly wrong restraint from a catastrophically wrong one.

---

## FAMILY C — CONDITIONAL ENSEMBLE GENERATION

**Status: MEASURED. Does not beat the incumbent.** The full-instrument cascade emits **3.321 Å**
against 3.204 Å (+0.117 [+0.050, +0.188]); its unprojected consensus is at parity with an interval
including zero. The ensemble is good (**G = 2.760**, −0.444 [−0.573, −0.325]) and the readout is bad —
selection costs +0.751, aggregation recovers −0.360, projection costs +0.170, and pre-filtering the
ensemble by the objective makes aggregation *worse*. The collapse check passes (mean pairwise spread
0.527 Å; 17/126 below 0.25 Å), so the aggregation result is inside the set-mean law's domain.

**Hypothesis.** The right output is `P(x | S, O)` — a bundle of restraint-consistent conformers,
which is what NMR structure determination actually reports — rather than a single argmin.

**Why it is distinct.** It is evaluated on ensemble statistics (coverage, diversity, calibration)
that A and B do not produce, and it is where the project's largest measured lever lives.

**The lever, already measured.** Coordinate aggregation is worth **3.4× the entire contribution of
any objective**: averaging buys 0.587 Å where the whole contribution of any ranking is 0.171 Å. And
the terminal operator consumes the set **mean**, not the set **best** — `d_out = 1.16·d_set_mean +
0.04·d_set_best` at R² = 0.89. A perfect rank-1 selection is worth −1.74 Å through an argmin and
**−0.03 Å** through the m = 75 average.

**The four-stage readout this family requires**, and which all families now report:

| | |
|---|---|
| **G** | best structure anywhere in the ensemble — **ORACLE readout** |
| **S** | what the method identifies without the native — the objective's argmin |
| **A** | what an estimator recovers by combining — coordinate consensus |
| **F** | what is finally reported, after the validity gate and refinement |

`S − A` decides whether selection matters at all.

**A cheap corollary being tested separately.** Since the operator consumes the set mean, and a
restraint-fitted conformer (3.644 Å) is **0.43 Å better than the retrieval set it would join**
(4.072 Å), mixing solved conformers into the retrieved pool should lower the emitted structure
*without the generative arm needing to beat the incumbent at all*. This needs only the fits to beat
the set mean, which is already measured to hold. Its most likely failure mode is stated in advance:
the fits share the distogram with the incumbent's filter, so their errors may not be independent of
the retrieval errors.

**What would falsify the family.** Consensus over the ensemble being no better than the single best
fit.

---

## FAMILY D — RESERVE

**Status: RESERVE — TESTED AND LEFT EMPTY.**

The candidate was alignment engineering, and it was measured across 22 arms on 126 targets. **The
mechanism is real and enormous**: an ORACLE rotation of the fit's own true error out of the loud
half, at **fixed magnitude**, gives **1.855 Å** — −1.822 [−2.037, −1.613], W/L 121/5. **Nothing
native-free reaches it**: the best leave-fold-out arm is −0.007 [−0.074, +0.055], and no arm moves
alignment by more than 0.104 against the 0.387 required. Pooled over 2,142 (target, arm) pairs, the
partial correlation of ΔRMSD with Δalignment given Δraw error is **+0.087**, against β = +0.586 for
distance MAE.

> **Among achievable interventions, alignment is a description of what good channels happen to have,
> not a lever that can be pulled. The slot stays empty.**

The workstream's one native-free gain — `tik_LFO` −0.216 Å [−0.345, −0.094] — turned out under a
matched-shrinkage control to be **ordinary isotropic Tikhonov regularisation**, which beats both
geometry-aware metrics at every λ. It does not qualify for the frozen protocol: it improves a losing
arm from 3.677 Å to ≈3.46 Å, still well short of the incumbent's 3.204 Å, and its λ grid is not
converged.

Declared unallocated at the start, to be filled only if A–C proved unable to attack the measured
bottleneck. It is **not** filled speculatively, and it is **not** filled merely because a fourth
family would look tidier in a paper.

**The candidate that has emerged, and the evidence for it.** The information-channel audit found
that error *structure* prices a channel far better than error *magnitude*: a real 64°-RMS channel
emits **1.77 Å** where an i.i.d. channel of the same magnitude emits **4.71 Å**, and the mechanism
is subspace alignment (`||Je||/(||e||·s_rms)` = 0.565 against a 0.945 random-direction null), with
causality established by surrogate destruction — **sign-flipping errors while preserving every
magnitude exactly costs +1.596 Å [+1.385, +1.813]**. It also found that **where coverage gaps fall
beats how many there are**: terminal gaps at coverage 0.7 give 1.716 Å against 3.505 Å for
mid-chain, a 1.79 Å spread, and **no channel exploits this.**

A family built on **engineering alignment rather than inheriting it** — Jacobian-weighted
restraints, terminal-relaxed fitting, fitting restricted to well-conditioned singular directions —
would be genuinely distinct from A–C, because its distinguishing experiment is one none of them can
run: does emitted RMSD improve *without* raw error improving?

That question is under measurement now. **Family D is filled only if the answer is yes.** If
alignment turns out to be a description of what good channels happen to have rather than a lever
that can be pulled, the slot stays empty and that is recorded as the result.

---

## WHERE THE FOUR MANDATED COMPONENTS SIT

Each has a measured role, not a decorative one, and each role is defended by a control.

| component | role | the control that defends it |
|---|---|---|
| **VQE** | samples/optimises the discrete torsion configuration in A | best-of-N from the **untrained** circuit at matched budget |
| **CVaR** | shapes the ensemble tail; α is a diversity dial in C | an α sweep, plus arms with the corrected gradient baseline |
| **Legacy** | **clash gate** — its only measured-honest role | the gate's exclusion count is reported, never hidden |
| **AMBER** | validator and refiner in A and C; the **objective** in B | matched-budget classical search, at the corrected 8–14 ms cost |

Legacy's eleven weights were never fitted — documented as "a variance-balanced starting point, not a
fit" — which is why its role is confined to gating. AMBER's role has been re-measured under a
pre-registered replication, and it splits in two.

**The ACCURACY claim is WEAKENED.** The −0.022 Å effect replicates in sign and magnitude (5/5 draws,
sd 0.0076 Å) and is **not** a multi-start artefact — this pipeline has no RNG and reproduces
bit-identically across interpreters. But it is **absent on the frame-reproducible 69% of the
instrument** (−0.009 [−0.024, +0.006], a *losing* 38W/46L, replicated on three frames), it is
**erased by an exact rotation null that must be zero** (+0.0117, max 1.51 Å), and four targets are
**silently non-converged** — up to 8.9 × 10⁸ kcal/mol — because `refine_coords` has no convergence
gate. **It may no longer be quoted as "−0.022 Å at valid geometry."**

**The VALIDITY claim is CONFIRMED and was understated**, once two Sprint 14 reporting defects were
corrected (Ramachandran scored on a different structure from its RMSD; clash rate a hard-coded
placeholder): **0.466 → 0.874 Ramachandran-favoured, 116W/2L**, and **1.397 → 0.000 clashes,
63W/0L**. AMBER's defensible role is **stereochemical repair**, and on that axis it is larger than
the record said.

---

## WHAT THIS ARCHITECTURE DOES NOT CLAIM

- It does not claim the generative frame currently beats the incumbent. On the full instrument the
  raw arm **loses by +0.440 [+0.290, +0.592]**.
- It does not claim the distance channel can be improved — only that improving it is now worth about
  3.19× on the floor and 2.07× on the addressable headroom.
- It does not claim quantum advantage. The declared bar is beating an untrained circuit at matched
  budget, and no previous sprint has cleared it.
- It does not treat any ORACLE ceiling as an achievable result.
