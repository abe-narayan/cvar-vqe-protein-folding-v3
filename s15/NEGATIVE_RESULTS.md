# NEGATIVE RESULTS, REFUTATIONS, AND RETRACTIONS

A required deliverable of this programme, and the document most likely to be its most useful
output. Every entry is a thing that was believed, tested, and found false — including things
believed by this project's own earlier sprints, and including things believed by the coordinator
earlier the same day.

**Nothing in this file is deleted when it becomes inconvenient.** Entries are amended by appending,
never by editing away. Where a claim was published internally before its control was run, that fact
is recorded, because the sequence matters more than the conclusion.

Status: **living**. Entries marked *(in flight)* have a run in progress and will be completed
before the protocol freeze.

---

## PART I — REFUTATIONS FROM THIS SPRINT

### N1. Naive distance-geometry generation is WORSE than the incumbent it was meant to replace

**REFUTED on the full instrument.** The sprint's opening move was to reframe the problem from
*selecting* among retrieved candidates to *solving* for a restraint-consistent conformation. On the
first 8 targets this looked promising and the write-up said as much, with the caveat that "the
full-instrument number is what counts."

The full-instrument number, all 126 targets:

| arm | mean | median | <2 Å | vs incumbent (paired) |
|---|---|---|---|---|
| **ORACLE**, fit to the TRUE distance matrix | 0.611 | 0.055 | 0.86 | −2.593 [−2.901, −2.295] |
| predicted distances, inverse-variance weighted | **3.644** | 3.466 | 0.16 | **+0.440 [+0.290, +0.592]** |
| predicted distances, unweighted | 3.798 | 3.625 | 0.12 | +0.594 [+0.430, +0.761] |

The 8-target read put the predicted arm at 2.792 Å. The real figure is 3.644 Å, and it **loses to
the incumbent with a confidence interval excluding zero**. The reframing does not, by itself, buy
anything.

The mechanism is nonetheless clean, and this is why the negative is informative rather than merely
disappointing: both arms use the identical optimiser, starts, selection rule and geometry, so the
3.03 Å between them is **entirely restraint error**, with nothing else able to carry it.

### N2. My own "quantisation ceiling" was an indexing bug — retracted within the hour

**RETRACTED.** Believed and written up: that the distogram's 17-bin output representation costs
about 1.06 Å of achievable accuracy, inferred from an ORACLE maximum-likelihood arm (1.756 Å)
sitting that far above an ORACLE least-squares arm (0.697 Å) on matched targets. The story was
coherent, it explained why maximum likelihood underperformed least squares despite being strictly
more informative, and it fit prior results.

It was wrong. The likelihood table computed its bin index as `floor((d − c₀)/(c₁ − c₀))` — a
**uniform-grid** interpolation on a grid that is not uniform. The distogram's centres run 0.5 Å
apart at short range and up to 4.0 Å apart at long range, so:

| query | bin actually read | that bin's centre |
|---|---|---|
| 6.25 Å | 3 | 5.75 Å |
| 8.50 Å | 6 | 7.25 Å |
| 15.00 Å | 14 | **17.50 Å** |
| 25.00 Å | 15 | **21.00 Å** |

Everything above 15.25 Å collapsed into one bin. The error grew with distance — precisely the
signature that had just been attributed to quantisation.

Void as a result: every `ml_*` arm computed before the fix, and the first full cascade run, which
was killed at 40/126 rather than allowed to finish and be quoted. Unaffected: every least-squares
arm, since they never touch the class.

**And the corrected measurement, for the record.** The binned representation does cost something —
an ORACLE **+0.381 Å [+0.156, +0.614]**, 111/126 targets worse, same sign in 5/5 folds. **So 75% of
the claimed 1.06 Å was the bug and 25% was real.** It is bin *placement* rather than count (a uniform
0.5 Å grid costs +0.243; a uniform 0.25 Å grid +0.085 with a CI crossing zero), and a continuous
density recovers 45% of it — while maximum likelihood **still** does not beat least squares, because
the mechanism is calibration rather than binning.

**The lesson, which is the reason this entry exists.** Nothing about the finding's internal
plausibility flagged it. It was mechanistically satisfying, quantitatively specific, and consistent
with everything around it. Only checking the lookup against the grid it claimed to use caught it.
This is the second time in two sprints that a satisfying mechanism has turned out to be a property
of the instrument rather than of the problem.

### N-FLOOR. An empirical false-positive floor of ~0.08 Å, demonstrated on a null effect

**DEMONSTRATED, and it binds every small effect in this programme.**

Maximum likelihood against a **constant-width Gaussian is least squares** — the two are the same
estimator, so the comparison between them is exactly zero by construction. Measured through this
machinery:

| start draw | result |
|---|---|
| one | **+0.081 [+0.014, +0.169]** — a confidence interval **excluding zero** |
| another | −0.003 [−0.071, +0.060] |

**A paired 95% interval excluded zero on a quantity that is provably zero, purely from the multi-start
draw.** The null-calibrated concentration test independently flagged that row — and only that row — as
**FAIL (p = 0.024)**, which is that check earning its keep.

The underlying variance is quantified: **sd 0.132 Å on an arm's absolute mean across four independent
start draws**, worst target 5.15 Å. Paired *differences* are far more stable (the same quantity at
0.381 / 0.416 / 0.486), which is why this programme's comparative claims survive while its absolute
constants need an error bar.

**Consequences, applied rather than noted.**

1. **Absolute constants get a start-draw bar.** The ORACLE ceiling is **≈0.6 Å with sd 0.13 Å**, not
   0.611 — an independent workstream's draw returned 0.569.
2. **Effects at or below 0.08 Å are not resolvable without replication across draws** — *in pipelines
   that have a random multi-start.* In this programme: `project_scaled_both` (+0.080),
   restraint-level fusion vs the distogram (−0.071), and the cascade's `A` stage (−0.053).

   **A LIMIT ON THE FLOOR, established afterwards.** I also applied it to the AMBER k = 30 relaxation
   (−0.022 Å). **That was wrong.** The floor's mechanism is a *per-process-salted random
   multi-start*, and the AMBER pipeline **contains no RNG at all** — `core.project.STARTS` is a fixed
   4-tuple, AMBER runs single-threaded with deterministic forces, and the whole arm reproduces
   **bit-identically in a fresh interpreter (max |Δ| = 0.000e+00, n = 126)**. Its own constructed
   floor is **≈0.004 Å gated**, for an unrelated reason. **A floor measured on one pipeline does not
   transfer to another by size alone; it transfers only if the mechanism does.**
3. **Every effect above the floor is unaffected**: the cascade gaps (+0.751, −0.360, +0.170, +0.117),
   coordinate averaging versus the better channel (−0.255), realizability (−1.391), the ceiling
   (−2.593), the quantisation cost (+0.381).

**Why it was invisible.** Determinism was checked repeatedly — within a process, and against cached
artefacts. `hash()` is salted per process, so **no check that never starts a second interpreter can
see this**. Seeding is now stable, and the frozen protocol is re-run under it before the benchmark.

### N-AMBER. The one positive physics result: accuracy WEAKENED, validity CONFIRMED and understated

**Replicated under a decision rule pre-registered before any draw returned.**

**The accuracy claim is WEAKENED and may no longer be quoted as "−0.022 Å at valid geometry."** It
does replicate — 5/5 draws same sign, mean −0.0233 Å, sd 0.0076 Å, clearing its pre-registered 0.0085
threshold — and it is **not** a multi-start artefact. But it fails three of its own tests:

1. **An exact null that must be zero returns non-zero.** Relaxing the same structure in a **rotated
   lab frame** is zero by construction (ff14SB/GBn2, the restraint and every RMSD are
   rigid-invariant). Measured **+0.0117 [−0.0034, +0.0383]**, max **1.51 Å** — and on that frame the
   effect becomes −0.0102 [−0.0332, +0.0215], **CI including zero**.
2. **The effect is absent on 69% of the instrument.** On the 84 frame-reproducible targets it is
   **−0.009 [−0.024, +0.006] with a *losing* 38W/46L**, replicated on three independent frames. The
   whole effect sits on the 31% where the AMBER minimiser is not frame-reproducible. Sprint 14's
   "5/5 folds same sign" also does not fully replicate.
3. **The convergence defect below.**

**The validity claim is CONFIRMED and was badly understated** — see the two defects below — at
**0.466 → 0.874 Ramachandran-favoured (116W/2L)** and **1.397 → 0.000 clashes (63W/0L)**. That half
is nowhere near any noise floor. **AMBER's defensible role is stereochemical repair, not accuracy.**

### N-DEFECT-1. `refine_coords` has no convergence gate, and four targets are silently non-converged

**NEW, in no prior record.** 1D6X, 1MF6, 2NB7 and 7BX2 end minimisation **above 1000 kcal/mol** —
**1MF6 at 8.9 × 10⁸ kcal/mol** with 2 clashes — and are **silently scored into the published mean**.
No gate anywhere in the pipeline catches it.

**Those four targets alone move the exact null from −0.0005 to +0.0117 Å.** Gated, the null collapses
to **−0.0005 [−0.0046, +0.0035]** and the effect reproduces across frames to **0.0003 Å**. Four
targets out of 126, ungated, were carrying the entire apparent failure of the null.

### N-DEFECT-2. Sprint 14 scored one arm's Ramachandran on a different structure from its RMSD

**Arm A's Ramachandran fraction was measured on the λ = 0.3 arm's geometry while its RMSD was quoted
from the λ = 0 arm** — two different structures reported as one. Correcting it changes the baseline
from 0.734 to **0.466**, which is why the validity gain was understated by nearly half.

### N-DEFECT-3. Sprint 14's arm-A clash rate was a hard-coded placeholder

`[0.0] * len` — not a measurement. The true baseline is **1.397 clashes**, going to **0.000** after
relaxation on 63 targets with no losses. A row of zeros that looked like a measured result was never
measured at all.

### N3. Family B, as declared, is falsified before implementation

**REFUTED.** Family B was declared as constrained refinement: `min E_AMBER subject to
C_restraint ≤ ε` — the NMR structure-determination formulation, which would give AMBER a
first-class role. Its precondition is that the native satisfies the predicted restraint bands.

Measuring `z = |d_native − d̂| / sd` over 126 targets and 8,549 pairs:

| ε (units of the predictor's own sd) | fraction of pairs admitting the native |
|---|---|
| 0.5 | 0.248 |
| 1.0 | **0.458** |
| 2.0 | 0.712 |
| 3.0 | 0.840 |
| median z | 1.336 |
| max z | **7.237** |

At a one-sigma band **fewer than half** the restraints admit the native; at three sigma — already so
loose the restraints barely constrain — one pair in six still excludes it. A uniform-ε hard
constraint has a feasible set that excludes the answer at every ε that meaningfully constrains.

The redirection this implies is in RESULTS (K7): the failure is **concentrated**, not diffuse — a
typical restraint is only mildly wrong and a minority are catastrophically wrong — so the family
survives as a *robust* constrained formulation, not a uniform-band one.

### N4. The retrieval pool is a better distance ESTIMATOR and a worse OBJECTIVE

**Both halves demonstrated; the second refutes the natural expectation from the first.** The K=500
retrieval pool supplies a complete Cα–Cα distance matrix per member and had never been used for its
distances in fourteen sprints. As an estimator it beats the trained predictor with no training at
all: global MAE **2.249 vs 2.386**, with the opposite bias sign at short and medium separation.

As a *ranking objective* it is the worst arm measured. It places the native at the **66.5th
percentile** of its own retrieval pool — above chance, so it **actively disfavours** the native —
with ρ(f, RMSD) of only +0.365 and an argmin barely better than a random pool member (−0.427 Å,
winning on 56% of targets), against the distogram's −0.949 Å on 79%.

This is a new instance of a standing law of this project: **MAE does not price selected RMSD.** A
channel being more accurate on average does not make it more useful for choosing, and the two must
be measured separately every time.

### N5. Separation debiasing helps the fit and slightly hurts the ranking

**Demonstrated, and it constrains how the correction may be used.** The distogram over-predicts
distance systematically, from +0.113 Å at separation 2–3 to +1.492 Å at 11–15, and correcting this
leave-fold-out is legitimate and native-free. On the *ranking* axis it is marginally harmful: native
percentile 0.359 versus 0.347, argmin RMSD 3.559 versus 3.504.

A uniform additive shift of the restraint targets barely reorders candidates while it does move the
location of the minimum. So a correction worth having in **generation** need not be worth having in
**selection**, and neither result licenses the other.

### N6. Recalibrating the distogram's uncertainties cannot change a single structure

**Demonstrated analytically and confirmed numerically.** The distogram's `sd` is **2.633× over-
confident** — a genuine reporting defect that must be declared in any paper using it. It is
tempting to treat recalibration as an obvious improvement.

It is not, and the reason is worth stating because it saves a sprint. The standardised residual's
spread is nearly *constant* across separation (2.19, 2.98, 2.70, 2.70, 2.67), so the `sd` captures
the error's **shape** correctly and is wrong only by a near-uniform factor. **A weighted
least-squares argmin is invariant under uniform rescaling of all weights.** Recalibration alone
changes nothing about the emitted structure.

This also explains why the problem must be attacked through the **loss** rather than the weights: a
uniform rescaling cannot distinguish a mildly wrong restraint from a catastrophically wrong one,
which per N3 is exactly the distinction that matters.

### N7. Euclidean-distance-matrix projection is redundant here — considered and rejected without a run

**Rejected by reasoning, recorded so it is not re-proposed.** Predicted distances are produced per
pair independently, so the predicted matrix generally violates the triangle inequality and is not
realisable by any 3D structure. Projecting it onto the cone of valid Euclidean distance matrices
before fitting is a standard, principled preprocessing step this project has never performed, and
it looks obviously worth doing.

It is not, in this architecture. The torsion fit **already** enforces a stronger constraint than
EDM projection does: EDM projection enforces rank-3 embeddability, while fitting `(φ, ψ)` under
ideal bond geometry enforces rank-3 embeddability *and* fixed bond lengths and angles, reducing the
degrees of freedom from `3n − 6` to `2n − 4`. Under a squared loss with fixed weights, projecting
first and then fitting is very nearly the same operation as fitting directly. The denoising the EDM
step would provide is already being provided.

Recorded as reasoning rather than measurement, and flagged as such: if a robust or redescending
loss changes the effective weighting enough, the equivalence weakens, and the question could become
live again.

---

## PART II — REFUTATIONS INHERITED FROM SPRINT 14, CARRIED FORWARD

These were established by measurement in the previous sprint and constrain what may be claimed now.

### N8. Running the variational optimiser is worse than not running it

**0 wins out of 12** against best-of-N drawn from the **untrained** circuit, at **+0.65 to +1.32 Å**,
at matched budget. The mechanism was traced rather than asserted: optimisation concentrates the
sampling distribution onto a tail in which every objective is at chance, so the argmin sits at the
**41st–48th percentile of its own tail**, and `tail mean − tail best` = 1.93 Å matches the certified
1.876 Å and sampled 2.040 Å selection gaps.

The literature survey found **no precedent for this control anywhere**. The nearest published work
uses uniform random sampling — a strictly weaker control, since it does not share the ansatz's
support — and reports parity where we find 0/12. This is the strongest surviving novelty of the
programme and the paper leads with it.

Whether it survives a change of objective class, from selective to generative, is the sprint's
central quantum question. *(in flight)*

### N9. In-band discrimination is signal-limited, not sample-limited

A set-transformer over the full signed deviation map has a **flat learning curve**, while a leaked
label is loud at n = 8. More data does not fix in-band ranking; the signal is not there. In-band
ordering is learnable to 0.986 *within* a target and collapses to **0.600 across** targets.

### N10. Errors are decorrelated, and that cannot be exploited

Truth-partialled error correlations between channels run 0.04–0.26 — genuinely decorrelated — and
fusion nonetheless buys only **+0.004 to +0.011 Å**, because the fusion gain goes as the **square**
of the weaker channel's skill. Any claim that two channels should be combined must clear this
arithmetic before it is made. (K6's `ml_comb` clears it, but only just, and only on one axis.)

### N11. Consensus is outlier avoidance, not a nativeness signal

The native sits at the **82.8th percentile** of the criterion the consensus medoid minimises. The
mechanism is capped at the pool's mode, so it cannot exceed the pool's typical structure however
much it is refined.

### N12. Nothing ranks within the pool

Dev pool best 2.355 Å against 3.324 Å returned; all-atom AMBER reranking is worth **+0.004 Å**. The
best external corroboration is that AlphaFold2's rank-0 model is the lowest-RMSD one only **13%**
of the time against a **20%** null — state-of-the-art ranking performing *worse than chance* inside
its own pool.

---

## PART III — CLAIMS WITHDRAWN OR CORRECTED BY THE PHASE 0 AUDIT

### N13. Retraction: every AMBER row in `s14/results/obj_floor.json`

Computed on the full labelled subsample with **no `amber_kind` mask** (0.42 Å oracle-conditioned),
containing no stratification string, and flagged `"complete": true`. Corrected values exist in the
Sprint 14 findings; the artefact was never regenerated. Do not read that file.

### N14. Correction: the oracle-conditioning share was transposed

Stated as "40% oracle-conditioned". Measured across all 19 enumerated files: **21.7%**
oracle-conditioned, **40% unbiased**, 39% prior-conditioned. The −0.401 Å effect size was correct.

### N15. New defect: the "unbiased" AMBER stratum is not unbiased in its tails

The index array **force-includes the ORACLE snap index**, which lands in the `kind == 0` stratum on
16 of 19 files. There it is the **minimum-RMSD member on 6/19**, in the RMSD-lowest 1% on 10/16, and
**inside the < 1.5 Å in-band set on 14/19** — touching 40%/25%/17%/15%/11% of all in-band pairs on
five targets whose bands hold only 5, 8, 12, 13 and 19 members.

`kind == 0` is unbiased on the **mean** (+0.0003 Å) and not in the tails. Binding rule for every
tail, in-band, argmin or top-k statistic: `amber_kind == 0 AND amber_idx != snap_index`, with the
per-target n printed. The direction is safe — the headline can only have been flattered — so the
Sprint 14 negative stands and is if anything understated.

### N16. `I.FAIL18` is a threshold artefact that three sprints treated as an object

`BAND = 1.5 Å` has no derivation anywhere in the repository. The zero-recall count runs **45 → 2**
across BAND 0.5 → 3.0, and on a 99-combination sweep **exactly one of the 18 targets (9KAR) appears
in every version of the set** (mean Jaccard 0.498). The column remains in tables for cross-sprint
comparability; **no argument may rest on set membership.**

### N17. Budget matching over-charged AMBER by 2–3×

A single AMBER point costs **8.3 ms at 127 atoms rising to 23.3 ms at 237 atoms** — 8–14 ms on the
enumerated targets — not the 28 ms assumed. AMBER therefore received *less* search than a fair
matched budget would have given it, so its negative result is understated rather than overstated.

---

## PART IV — CLAIMS THE LITERATURE TAKES FROM US

Not refutations of our measurements, but refutations of our **novelty**. Recorded here because
overclaiming priority is a form of error this document exists to prevent.

### N18. "The certified optimum is worse than random" is published — for lattice contact potentials

Roget et al. (arXiv:2606.21241, 19 Jun 2026) enumerate every peptide up to length 15 across 12,446
PDB peptides and report that the minimum-cost conformation has, on average, a **larger** RMSD than a
randomly selected feasible one. They also give the mechanism: the Miyazawa–Jernigan potential was
parameterised on proteins above 50 residues, and at 9–16 residues ρ(cost, RMSD) is negative.

Our version survives only in its all-atom form — a certified optimum of AMBER ff14SB plus a
statistical potential in continuous torsion space, which their length-parameterisation mechanism
cannot cover — and must cite them as prior art for the phenomenon rather than presenting it as a
discovery.

Our K6 result is the useful counterpart: for a *distance-restraint* objective at the same chain
lengths, the argmin is about an ångström **better** than random, not worse. Their pathology is a
property of contact potentials, not of short-chain structure prediction.

### N19. In-loop CVaR for peptide folding is taken

Uttarkar et al., PLoS One 21(2):e0342012 (2026), "QuPepFold". Our CVaR work cannot be positioned as
"CVaR applied to peptide folding". What remains ours is the **defect triple** — findings about the
estimator, not applications of it — of which the gradient-baseline defect (baseline centred on the
tail only; cosine **−0.023** with the true gradient) is the most consequential.

### N20. Tail-restricted discrimination is partially taken as a framework

Truchon & Bayly (J Chem Inf Model 47:488, 2007) introduced BEDROC for exactly the early-recognition
problem our in-band metric addresses. We cite it and position our metric as an application.

### N21. Four Sprint 14 claims must be removed from the story entirely

Enumerated in the literature survey. They are removed rather than reworded.

---

## PART V — ERRORS OF PROCESS, NOT OF MEASUREMENT

Kept separate because they are the ones most likely to recur.

- **A result was published internally on a colleague's word before its control was run.** The
  "26W/24L is the good pattern" claim was circulated before the drop-top curve existed.
- **And the correction was also wrong.** "The concentration curve disqualifies it" was itself
  invalid: a raw drop-top threshold is not a test of concentration, because dropping the largest
  effects shrinks any mean. It needs a uniform-effect null simulation.
- **A sign error propagated into a shared brief.** "Small tail weighting collapses two AMBER
  variants" was backwards; concentration collapses them.
- **An `argmin` over a tied signal silently read the ORACLE sort order** and invented a 1.386 Å
  winner. Where a signal can tie, the outcome must be averaged over the tied argmin set.
- **A partial-overwrite hazard destroyed a control.** A module calling the ladder with a subset of
  arm names overwrote the shared results file with only its own arms, so a figure rendered without
  the control it was drawn to display. Fixed by giving the writer merge semantics.
- **"The first native-free accuracy improvement in fourteen sprints" was false** — an earlier
  sprint's consensus medoid was −0.172 Å, and the claim was made without checking.
- **Results produced before the seeding fix were quoted three times, after a document forbidding it
  had been written.** An adversarial audit found that re-running flipped the sign of a native-free row
  in one table and destroyed a rhetorical claim in another — so the reassurance that such numbers are
  "correct but not bit-reproducible" was too weak.
- **"Five times the headroom" appeared in eight documents including the paper's abstract, and the
  sentence introducing it contains arithmetic yielding 2.07.** The supported ratios are 3.19× on the
  floor and 2.07× on the addressable headroom. A number quoted often enough stopped being checked.
- **An ORACLE arm's numbers were presented as "real errors" in four documents**, and the control that
  matters — a zero-information constant α-helix scoring 0.709 on the same statistic — was dropped in
  relaying. The native-free range is 0.6–1.7 Å, not 0.6–2.9 Å.
- **A draw count was wrong in more than twenty places** (2,048 for 4,096), including a figure axis.
- **A figure guard written to stop smoke reads becoming figures was bypassed three ways**, one of them
  a headline panel, and one relayed panel drew 9 targets with no n printed.
- **Two figure captions asserted the opposite of the data drawn beneath them.**
- **A claimed quantum positive was led with in the dossier after the replication workstream had
  already demoted it** — the unit of analysis overstated n threefold, the surviving cell failed its
  concentration check, and the control was charged 2,048 evaluations against the arm's 819,200.

**None of these was found by thinking harder about the science.** All were found by loading the
artefacts and comparing them to the prose — the third time in this sprint that mechanical checking
beat intuition, after the non-uniform bin lookup and the per-process hash.

- **A conclusion was drawn from 2 targets and reversed on 9.** The "universal discrimination floor"
  claim. Small-n reads from the enumerated set have reversed a conclusion four times.
