# SPRINT 18 — DECISION AND RETRACTION LEDGER

Written as the sprint runs. Sprint 15 closed with eleven coordinator retractions, Sprint 16 with
five, Sprint 17 with six. The same standard applies.

---

## L1 — PHASE 0: the 19-target result reproduces exactly, and was never significant (2026-09-06, coordinator)

Full account in `s18/PHASE0.md`. Recomputed independently from `s17/results/quantum_theory.json`,
19 exhaustively enumerated targets, certified argmins (no search error).

**Every figure reproduces**: full objective argmin **2.661 A**, degree-1 argmin **2.411 A**, space
best 1.030 A, Walsh weight-1 variance fraction 0.613, cumulative <= weight 2 0.930. Phase 0's
reproduction gate is passed.

**But the inference was never statistically supported:**

    deg1 - full   mean -0.249   median +0.000   W/L 7/5   95% CI [-0.650, +0.093]

The interval **includes zero**; the **median is exactly 0.000** because 7 of 19 targets share an
identical argmin; and the mean is carried by three targets (7VI4 3.936 -> 1.396, 1CS9 5.256 ->
3.655, 7T3H 3.156 -> 1.791) against three that move the other way (8HVS +0.85, 1N9U +0.71, 5V5B
+0.59).

**And a tension pointing the other way**: rho(objective, RMSD) is **+0.264 for the full objective
and +0.153 for degree-1**. The truncation is a *worse* global correlate of RMSD even though its
argmin is better on average -- so whatever is happening is a property of the **argmin**, not of
alignment.

**Sprint 17's report is amended**: it quoted "2.411 vs 2.661" as a lead without its interval. That
is this programme's recurring failure mode -- a number quoted without the functional the claim
needs -- committed by the coordinator, and recorded here rather than in a footnote.

**Falsifier F4 is therefore live before the 126-target test begins.** Degree-1 enters Sprint 18 as
a hypothesis at the noise floor of its own instrument, not as a result to be translated.

---

## L2 — THE FUNCTIONAL FORM IS SOUND; THE DISTOGRAM'S ERROR IS WORSE THAN NOISE (2026-09-06, coordinator)

`s18/objceil.py`, n = 126. Full table in `s18/COORD_FINDING.md`.

Holding the functional form fixed and varying only the distances it must satisfy,
`d_alpha = (1-alpha) dhat + alpha d_true`:

- **alpha = 1 reaches 1.152 A** (median 0.307), beating the coordinate average on **114/126**,
  with **80.2%** of targets below 2.5 A and **73.8%** below 2.0 A. **The sum-of-squared-pair-
  residual form is not the problem.**
- **alpha = 0.5 reaches 2.448 A** -- halving the distogram's residual RMS clears the primary target.
- **alpha = 0 is 3.610 A**, 0.561 worse than the coordinate average (Sprint 17's L30, reproduced
  here exactly on an independent code path).

**And the controls invert the usual reading.** Destroying only the DIRECTION of the residuals while
keeping their magnitudes gives **2.609 A (-1.001 [-1.226, -0.784] against alpha = 0)**; matched-RMS
isotropic noise gives **2.573 A (-1.037 [-1.279, -0.801])**.

> **The distogram's errors are worse than random errors of the same magnitude.** Not merely large:
> structured, and the structure is actively harmful. This is the recorded *"full amplitude, wrong
> direction"* law measured at the optimum rather than at the predictions, and it is the sharpest
> form the programme has produced.

**Consequence.** The objective-engineering branch is probably a distraction *as a way of fixing the
form* -- but degree-1 acquires a better and more testable rationale: if the **pairwise** component
of the error is the mis-directed part, a degree-1 objective should be **more robust to structured
distogram bias**, and should close part of the 1.0 A gap between alpha = 0 and the shuffled control.
Communicated to all four workstreams.

---

## L3 — BOTH CONFOUND CONTROLS SURVIVE, AND THE SECOND MAKES THE RESULT SHARPER (2026-09-06, coordinator)

I raised two possible confounds against my own L2 headline when briefing the adversarial lane, then
took them myself because four workstreams were building on it. `s18/objceil.py`, n = 126, complete.

| control | RMSD | vs alpha = 0 | vs the average | W/L vs avg |
|---|---|---|---|---|
| `shuffled` (plain permutation) | 2.609 | -1.001 [-1.226, -0.784] | -0.439 [-0.682, -0.206] | 70/56 |
| **`shuf_strat`** (permuted **within separation bins**) | **2.560** | **-1.050 [-1.257, -0.860]** | -0.489 [-0.705, -0.273] | 72/54 |
| **`shuf_paired`** (residual **and its weight** permuted together) | **2.072** | **-1.537 [-1.758, -1.323]** | **-0.976 [-1.205, -0.766]** | **98/28** |
| `isotropic` (matched-RMS Gaussian) | 2.573 | -1.037 [-1.273, -0.802] | -0.476 [-0.749, -0.219] | 72/54 |

**Confound 1 REFUTED.** The plain shuffle also destroys the residual's correlation with sequence
separation, and the record says the distogram's bias grows with separation -- so the gain might have
been accidental shrinkage of heavily-weighted long-range residuals. Permuting **only within**
separation bins gives **2.560**, marginally *stronger* than the plain shuffle. It is not a
separation artefact.

**Confound 2 REFUTED, and the effect nearly doubles.** Permuting the residual **together with its
`1/sd^2` weight** -- preserving any adverse alignment between large residuals and large weights, and
preserving the model's own uncertainty calibration -- gives **2.072 A**, better than either other
shuffle and beating the coordinate average on **98 of 126 targets**.

### What the second control actually says, and it is the sharpest form yet

`shuf_paired` preserves the residual magnitude distribution **and** the residual-weight pairing, and
randomises only **which pair each error lands on**. It reaches **2.072 A** where the real distogram
reaches **3.610 A**.

> **The distogram's errors are harmful because of WHICH PAIRS they fall on.** Not their size, not
> their relationship to the model's own confidence -- their *assignment*. The errors are
> concentrated on the pairs that most determine the structure, and the model's uncertainty estimates
> do not compensate.

And the magnitude of the opportunity is stated by the same number: **a randomly re-assigned version
of the distogram's own errors already lands at 2.072 A, below the sprint's primary target**, while
the real one lands at 3.610.

### The concrete lead this creates

A native-free re-weighting of the objective by **geometric leverage** -- how much a given pair's
residual can distort the structure -- is now the cheapest intervention with a measured ceiling
behind it. Note the caution: Sprint 17 refuted 58 functionals of the distogram including uncertainty
weighting, but that was for **selection**, and this is **refinement**. Different use, different
functional; it is not a repeat, and it should be pre-registered as its own experiment rather than
assumed.

---

## L4 — RETRACTION: `shuf_paired` priced the WEIGHTS, not the error assignment (2026-09-06, ADVERSARIAL caught it; coordinator's error)

`s18/objceil.py` line 163 is

    p_, q_, _f = A.fit(pw, sd[pi], i, j, phi0, psi0)          # sd[pi], NOT sd

so the arm permutes the residual **and** hands a permuted weight vector to the fit. It is therefore
not the deployed functional sum_p (d_p - dhat_p)^2 / sd_p^2 at all. `shuffled` (residual permuted,
weights kept in place) = 2.609 and `shuf_paired` = 2.072 differ in **exactly one thing** — whether
`sd` moves — so the whole 0.537 A between them is the price of permuting the WEIGHTS.

**RETRACTED**: "the distogram's errors are harmful because of which pairs they fall on, and the gap
from 3.610 to 2.072 is entirely error assignment." That sentence had already been sent to EXPERIMENT
as the justification for a leverage-weighting arm; the retraction was issued to that lane before it
built on it, and the arm now stands or falls on its own pre-registration.

**WHAT SURVIVES, unchanged**: `shuffled` 2.609 [-1.226,-0.784], `shuf_strat` 2.560
[-1.257,-0.860], `isotropic` 2.573 [-1.273,-0.802] — all three carry the CORRECT weights and all
three beat the real distogram's 3.610 by ~1.0 A. So "the distogram's errors are worse than random
errors of the same magnitude" is untouched. Only the *localisation* of the mechanism is withdrawn.
The separation confound (L3, attack 1) and the weight-alignment confound as originally posed are
still refuted.

**WHAT IT OPENS, and it is worth more than what it closed.** The isolating controls are native-free
and therefore deployable, not ceilings: `wperm_only` (real dhat, permuted weights) and `wflat` (real
dhat, UNIFORM weights). Neither reads a native distance. If `wflat` moves the refinement, the
finding is that **the distogram's own 1/sd^2 confidences are anti-informative inside the refinement
objective** — the first native-free lever in three sprints. Owned by ADVERSARIAL
(`s18/q_ceilattack2.py`, n=126); the coordinator is deliberately not duplicating it.

**Process note.** Seventh instance in three sprints of a control pricing a different quantity than
its name, and the FIRST caught before it left the building. The catch came from the adversarial lane
reading the source rather than the table — which is the argument for keeping that lane funded.

---

## L5 — THE DEGREE-1 BRANCH IS CLOSED: F2, F4 and F5 all fire (2026-09-06, ADVERSARIAL + MATH, independently)

The sprint's central question is answered, and the answer is that **the 19-target 2.411 A was never
a property of the objective.** It is an artefact of which 2-bit code names which torsion state.

| | |
|---|---|
| RA (residue-additive, the gauge-INVARIANT object) vs full | **-0.004 A [-0.390, +0.318]**, W/L 3/5 with 11 ties — 5% of the instrument's 0.084 A MDE |
| W1 (strict Walsh <=1) under a RANDOM relabelling of the k=4 states | **+0.191 A** vs full |
| the codebase's own labelling | gauge percentile **0.00**, p < 0.0001 |
| encoding luck: W1(identity) - W1(gauge mean) | **-0.440 [-0.844, -0.133], folds 5/5** |

That last interval is **the only one in the entire degree-1 story that excludes zero, and it
measures the arbitrariness of a bit encoding.** RA is invariant under relabelling; W1 is not; the
shipped labelling sat at the extreme of its own gauge orbit. MATH reached the same numbers
independently (2.862 vs 2.852; RA invariant 19/19).

Phase 0's "better argmin, worse rho" tension resolves **against** degree-1: decomposing
argmin = band mean + within-band draw, 143% of the -0.249 is the within-band min-of-N draw (CI spans
zero), and the only significant component is BAND QUALITY at **+0.108 [+0.011, +0.251], W/L 4/15,
folds 5/5** — pointing the wrong way. Tie trap checked and clean (all 19 argmin_ties = 1).

**COORDINATOR ERROR, recorded as such.** `s18/BRIEF.md` section 4 asserts that under a uniform
lattice measure the first-order functional ANOVA object *is* the Walsh weight-<=1 projection, and
says to verify it numerically. It is **FALSE as written**, and I then leaned on it as though
verified. The exact statement, verified at 1e-12: the residue-additive ANOVA object equals the
projection onto Walsh coefficients supported INSIDE ONE RESIDUE — weight 0, weight 1, **and
intra-residue weight 2**. They differ by 1.4 objective sd. Label: EXACT. The brief carried a
mathematical error into three workstreams; §4 is corrected in place with this note.

---

## L6 — THE QUANTUM BRANCH IS CLOSED, and the answer was forced before the arms were run (2026-09-06, ADVERSARIAL)

Walsh spectrum on the 19 enumerated registers, rank-uniformised (`s18/q_report anova`):

| object | weight-1 | inside 1 residue | INTER-RESIDUE | mean Pauli wt |
|---|---|---|---|---|
| full (deployed) | 0.6133 | 0.9132 | **0.0868** | 1.538 |
| RA (residue-additive) | 0.6706 | 1.0000 | **0.0000** | 1.329 |
| W1 (strict Walsh <=1) | 1.0000 | 1.0000 | **0.0000** | 1.000 |

The truncation **raises** weight-1 above the full objective's 0.613 and drives inter-residue coupling
variance to exactly zero, so the sprint's conditional quantum test cannot fire on this objective —
the prediction the coordinator registered in advance. Confirmed operationally: greedy 1-opt certifies
the truncations' global optimum in **100% of cells at 36 evaluations** (one coordinate pass) against
63.2% on the full objective, and both truncations' optima are closed forms at ~38 table reads.
RA-full at 36 evals is +0.368 [+0.200, +0.529], folds 5/5.

**The entangled-vs-CNOT-free arms were run anyway, with the pre-registration recording that the
answer had already been forced.** Null on the consumed readout at every alpha, on RA, W1 and full.
Q1/Q2/Q3 all fired. That sequence — answer forced first, arms run anyway, order pre-registered — is
the difference between a null and a null that could be accused of having been arranged, and it is
reported that way.

---

## L7 — MATH closes the branch independently, and names a MECHANISM nobody had: the degree-1 object IS the torsion prior (2026-09-06, MATH)

Three findings, in descending order of how much they damage the hypothesis.

**1. The sprint was porting the wrong object.** The brief's own ANOVA bridge derives the
residue-additive (RA) object, not the strict Walsh weight-<=1 (W1) object. RA's certified argmin is
**2.657 vs full 2.661**: **-0.004 [-0.380, +0.307], 3W/5L/11 ties, folds 1/5** — 5% of the 0.084 A
MDE, measured **on the very instrument where the effect was born.**

**2. W1's advantage is a statement about a bit encoding.** Relabelling which 2-bit code names which
torsion state — pure bookkeeping, no structure, energy or ordering changes — moves W1's argmin to
**2.862 +/- 0.135 A** over 24 relabellings, **23/24 worse than the full objective's 2.661**, with the
shipped labelling's 2.411 at the **0th percentile**. RA is **exactly invariant, sd 0.000 on 19/19**.
Independent reproduction of ADVERSARIAL (gauge mean 2.862 vs 2.852). The concrete failure is named:
of the three ways to split 4 states into pairs, the encoding privileges two as weight-1 and discards
the third; only 8 of 24 permutations preserve that.

**3. THE MECHANISM, and it is the sprint's best explanation of why nothing transfers.** `hamil` is
**75% an exactly-additive torsion prior.** Measured exactly, 19/19: prior term **1.0000**
residue-additive, distogram term **0.404**, blend **0.9495** — and
**rho(RA(hamil), the pure torsion prior) = 0.986**, argmins matching (2.657 vs 2.615).

> **The degree-1 object on the 19-target instrument is, to three nines, the retrieval-pool torsion
> prior.** And `s17/refine.py`'s objective contains **no prior term at all** — so the continuum
> truncation has nothing to reduce to. F4 and F5 are not just "it didn't transfer"; the thing that
> was being truncated does not exist in the continuum objective.

**THE CONTROL THAT CUTS THE OTHER WAY, and must be quoted beside the closure.** Zero-info 4.003;
ORACLE ceiling 1.030; **matched-random projection at the SAME retained variance 3.552**, which W1
beats by **-1.140 [-1.714, -0.561]**. Truncating *by degree* is genuinely not truncating at random.
The branch is closed on transfer and on gauge, **not** on "degree structure is meaningless".

**THE THEOREM (EXACT), verified rather than asserted.** RA = projection onto Walsh coefficients
supported inside one residue (weight 0 + all weight 1 + intra-residue weight 2); W1 = qubit-ANOVA
order 1. Residuals **8.9e-16 / 1.1e-15**, orthogonality 5.5e-17, variance budget 5.0e-16, on 19
targets, raw and rank-uniformised. W1 retains **0.613** of the variance, RA **0.913** — W1 discards a
third of the per-residue field.

**A STRUCTURAL COROLLARY WITH TEETH (EXACT).** The objective is a function of the **n-2 interior
residues**: `f_0 = f_{n-1} = 0` identically. **No arm of this objective can place the two termini
that the frozen full-chain metric scores.** This is a standing constraint on every distance-objective
arm in the programme, not a degree-1 result, and it belongs in the "missing angstroms" accounting.

**MU SENSITIVITY, kept as sensitivity rather than rescue.** RA under uniform / pool / rama:
-0.004 / -0.066 / -0.021, every interval covering zero, argmin identical across all three on 14/19.
The largest estimate any mu gives is still below the MDE. And pool and rama **do not factorise over
qubits on any target**, so W1 is not even an ANOVA under any informative mu.

**DELIVERABLE**: `s18/math_iface.py` (contract `s18/MATH_to_EXP.md`) — `E_full/E_res/E_ang/E_ge2/
E_lambda`, analytic gradients (<=2.9e-9), **`E_full` identical to `s17/refine._obj` to 2.3e-13**,
`E_lambda(.,1) - E_full = 0.0`, Control D PASS at 1.0e-15.

**ONE OPEN ITEM POINTING AGAINST THE COORDINATOR.** d(theta*)/d(dhat) computed analytically: at
n=29 (smoke, stopped) **there is no sign of the predicted robustification, and both point estimates
go the other way** (+0.407 [+0.135, +0.684] unweighted, +0.212 [-0.040, +0.479] sd-scaled). Labelled
OPEN, correctly, and recorded as a smoke.

MATH's verdicts: F1 not tested in this lane · **F2 FIRED** · **F3 FIRED on the gauge control, not
the variance control** · **F4 FIRED with the mechanism named** · **F5 SPLIT — the construction is
exact, the transfer of the published result is invalid.**

---

## L8 — ATTACKS 3 AND 4 LAND: the conclusion strengthens, the NUMBERS die (2026-09-06, ADVERSARIAL)

Both attacks were raised by the adversarial lane against the coordinator's own headline and both
were run by that lane. **The qualitative claim survives; two quoted figures do not.**

### Attack 3 — alpha=1's 1.152 A is a BASIN, not a ceiling

The parameterisation floor, measured by projecting the NATIVE into the same ideal-geometry torsion
parameterisation: **mean 0.083 A**, median 0.060, max 0.435. The parameterisation can represent the
native essentially exactly and is not the constraint anywhere.

Multi-start at alpha=1, same machinery, same weights:

| start | start RMSD | end RMSD | median | objective at end | W/L vs avgproj |
|---|---|---|---|---|---|
| avgproj (the coordinator's) | 3.213 | 1.152 | 0.307 | **3.04** | 0/0 |
| member | 3.571 | 1.253 | 0.554 | 3.86 | 55/71 |
| **helix (ZERO-INFORMATION)** | 4.070 | **1.028** | 0.412 | 3.68 | 61/65 |
| randtors | 5.258 | 2.080 | 2.278 | 18.38 | 37/89 |
| nativeproj_ORACLE | 0.083 | **0.082** | 0.033 | **0.30** | 92/33 |

`nativeproj_ORACLE - avgproj = -1.070 [-1.321, -0.829]`, median -0.094, 92W/33L.

**The objective column decides it.** From the ORACLE start the fit reaches objective **0.30**; from
the coordinator's start it stops at **3.04** — ten times higher, same objective, same weights. The
optimum is not at 1.152; **the optimiser stops there.**

**But the median must be read before anything is rewritten**: avgproj's median is 0.307 and the
paired median difference is only **-0.094**. On more than half the targets the deployed start
already reaches the native, and *the 1.152 mean is carried by a minority that fall into a distant
basin*. Honest statement: **alpha=1 reaches the objective's optimum on most targets and a distant
basin on a minority, and the mean reports the minority.**

**CONSEQUENCE 1 — the coordinator's conclusion is CORRECT AND UNDERSTATED.** With true distances
this functional points at the native; the ceiling is **0.083 A, not 1.152**.

**CONSEQUENCE 2 — the alpha ladder is NOT a requirement curve.** Every rung is one refinement from
one start, so each rung mixes "better target distances" with "better-conditioned landscape and
larger basin". The ladder measures objective quality **and** basin reachability jointly and cannot
be read as pricing the former alone. Every downstream use of it as a requirement curve is withdrawn.

**CONSEQUENCE 3 — a control the coordinator's arm did not carry.** A constant ideal alpha-helix
start (zero-information, 4.070 A) reaches **1.028**, as well as the coordinate average:
-0.124 [-0.366, +0.109]. **The coordinate average is not a privileged start at alpha=1.**

### Attack 4 — "halve the residual -> 2.5 A" is one point in a 0.55 A band that straddles 2.5

**First, a correction to the coordinator's own framing of the attack.** The ladder
`d_alpha = d_true + (1-alpha)(dhat - d_true)` **is exactly a uniform rescaling of the residual
vector** — it already preserves the residual's direction and its entire correlation structure. The
coordinator's request to "shrink the residual preserving correlation structure" is what the ladder
already does. The live objection is a different one: **a real predictor does not improve uniformly
across pairs.**

A family of error models **all at identical residual RMS** (0.5x deployed, the 2.5 A rung),
differing only in *where* the improvement lands:

| error model | resid RMS | RMSD | median | vs unif | W/L |
|---|---|---|---|---|---|
| unif (= the coordinator's a0.5) | 1.604 | 2.453 | 2.509 | +0.000 | 0/0 |
| **fix_confident** | 1.604 | **2.145** | 1.865 | **-0.308 [-0.462,-0.150]** | 91/35 |
| **fix_short** | 1.604 | **2.159** | 1.969 | **-0.294 [-0.424,-0.167]** | 87/39 |
| fix_big | 1.604 | 2.324 | 2.157 | -0.129 [-0.267,+0.004] | 72/54 |
| debias_sep | 1.604 | 2.349 | 2.177 | -0.104 [-0.207,+0.002] | 64/62 |
| fix_small | 1.604 | 2.487 | 2.319 | +0.034 [-0.065,+0.143] | 66/60 |
| **fix_long** | 1.604 | 2.618 | 2.643 | **+0.165 [+0.028,+0.303]** | 48/78 |
| **fix_unconfident** | 1.604 | **2.697** | 2.578 | **+0.244 [+0.097,+0.386]** | 35/91 |

**Spread 2.145 to 2.697 A at identical residual RMS — range 0.552 A, straddling 2.5.** Residual RMS
does **not** determine the outcome. The requirement is quoted from here on as *"halving the residual
RMS reaches 2.15-2.70 A depending on which pairs improve"*, never as 2.5 A.

**A USABLE DIRECTION, and it is the opposite of the intuitive one.** Improving **short-range** pairs
(-0.294) and pairs the model was **already confident about** (-0.308) buys more than improving
long-range or low-confidence pairs — both of which are significantly **worse** than a uniform
improvement. A predictor that reduces RMS by fixing its worst, longest-range, least-confident pairs
lands at 2.62-2.70, not 2.45. This is a live instruction for any predictor work and it inverts the
natural engineering instinct.

---

## L9 — THE PRIOR ARM: primary hypothesis REFUTED, pre-registered falsifier FIRED (2026-09-06, coordinator)

`s18/priorfit.py`, n = 126 complete, `s18/results/priorfit.json` (`complete: true`). The arm MATH's
mechanism opened: the deployed distance objective **plus** lam x a native-free per-residue von Mises
torsion prior fitted to the retrieved windows.

**Reproduction gate PASSES at full scale.** `prior0` = **3.610**, identical to `s17/refine.py`'s
`refine_full`. The arm is the same fit plus a term.

| arm | RMSD | median | vs avg [95% CI] | W/L | disp |
|---|---|---|---|---|---|
| coordinate average (start) | **3.048** | 2.837 | -- | -- | -- |
| prior0 (= s17 refine_full) | 3.610 | 3.370 | +0.561 [+0.405,+0.715] | 31/95 | 5.034 |
| prior1 | 3.510 | 3.335 | +0.462 [+0.343,+0.593] | 31/95 | 5.202 |
| prior3 | 3.423 | 3.155 | +0.374 [+0.264,+0.486] | 30/96 | 5.023 |
| **prior10 (best)** | **3.387** | 3.036 | **+0.339 [+0.241,+0.449]** | 34/92 | 4.960 |
| prior30 | 3.399 | 3.110 | +0.351 [+0.234,+0.473] | 37/89 | 4.972 |
| prioronly | 3.886 | 3.632 | +0.838 [+0.567,+1.122] | 36/90 | -- |

**1. THE PRIMARY HYPOTHESIS IS REFUTED. No lam beats the coordinate average.** The best rung
recovers **0.222 A of the 0.561 A refinement penalty — 40% — and still ends +0.339 [+0.241, +0.449]
worse than doing nothing.** Adding the missing prior term does not rescue refinement toward this
objective.

**2. THE PRE-REGISTERED FALSIFIER FIRES.** It was written as: *if the best `prior{lam}` is within
the instrument's 0.084 A MDE of `helix{lam}`, the prior's positional information contributes nothing
and the finding is "the fit needed regularising".*

    lam=10   prior - helix   -0.024 [-0.138, +0.088]     74/52
    lam=3    prior - helix   +0.005 [-0.081, +0.092]     69/57
    lam=1    prior - helix   -0.016 [-0.084, +0.055]     67/59
    lam=0.03 prior - helix   +0.031 [+0.005, +0.063]     65/61   (prior significantly WORSE)

**At every lam the retrieval-pool torsion prior is statistically indistinguishable from a constant
ideal alpha-helix, and at the smallest lam it is significantly worse.** The best torsion channel
measured anywhere in this project -- phi 33.6 deg, psi 59.2 deg, beating the trained leave-fold-out
sequence predictor -- buys **nothing over a zero-information constant** inside this objective. The
0.222 A that the term does recover is **REGULARISATION of a large-displacement fit**, and is
reported as that smaller claim.

**3. ONE REAL SECONDARY EFFECT, and it is a precise statement.** The `shuf` control -- the SAME
per-residue (mu, kappa) tuples permuted across residues -- IS significantly worse than the true
prior at large lam: -0.108 [-0.211, -0.006] at lam=10, -0.138 [-0.261, -0.016] at lam=30.

> **A WRONG positional assignment hurts; the RIGHT one is worth nothing over making no positional
> claim at all.** The prior's information is detectable and unusable in the same measurement.

**4. THE PATTERN THIS COMPLETES, and it is the sprint's most cross-cutting observation.** This is
the **third** independent zero-information control in Sprint 18 to match an informative arm:

  * ADVERSARIAL attack 3 (L8): a constant alpha-helix START reaches 1.028 A at alpha=1 against the
    coordinate average's 1.152 -- **the coordinate average is not a privileged start.**
  * this arm: a constant alpha-helix PRIOR matches the retrieval prior at every lam.
  * and the programme already recorded, before this sprint, that a zero-information constant
    alpha-helix beats the random control in torsion space.

On this instrument, **generic Ramachandran plausibility reproduces most of what the project's
best sequence- and retrieval-conditioned torsion channels deliver through an objective.** Any future
torsion-channel arm must carry a constant-helix control or it is uninterpretable.

**Displacement is NOT the mechanism, and that refutes my own stated reason for expecting a gain.**
The prior recovers 0.222 A while moving displacement only 5.034 -> 4.960 rad (1.5%). I predicted the
gain would come from shrinking a large-displacement fit; the displacement barely moved. Whatever the
0.222 A is, it is not the displacement reduction I pre-registered as its cause.

---

## L10 — AMENDS L4: the weights are NOT anti-informative. 1/sd^2 is VALIDATED, native-free (2026-09-06, ADVERSARIAL)

The isolating arms landed. `s18/results/q_ceilattack2.json`, n = 126.

| arm | reads native? | RMSD | median | vs a0 = 3.610 | W/L |
|---|---|---|---|---|---|
| a0_wkeep (deployed) | no | **3.610** | 3.370 | +0.000 | 0/0 |
| wperm_only | no | 3.701 | 3.547 | +0.091 [-0.017, +0.197] | 49/77 |
| **wflat** | no | **3.763** | 3.654 | **+0.153 [+0.066, +0.247]** | 50/76 |
| shufr_wkeep | ORACLE | 2.635 | 2.570 | -0.975 [-1.170, -0.783] | 102/24 |
| shufr_wperm (= `shuf_paired`) | ORACLE | 2.167 | 2.055 | -1.442 [-1.680, -1.216] | 113/13 |
| shufr_wflat | ORACLE | 2.440 | 2.301 | -1.169 [-1.396, -0.942] | 106/20 |
| a1_wkeep | ORACLE | 1.152 | 0.307 | -2.458 [-2.775, -2.149] | 117/9 |
| a1_wperm | ORACLE | 1.162 | 0.388 | -2.448 [-2.743, -2.153] | 117/9 |

**1. "sd is anti-informative inside the fit" is REFUTED, and the incumbent's design is VINDICATED.**
On the real distogram the deployed 1/sd^2 weighting is the best of the three. Uniform is
**+0.153 [+0.066, +0.247] WORSE**, CI excluding zero. Permuted is +0.091, spanning zero.
**Deleting the distogram's confidence estimates costs 0.153 A.** This is a **native-free positive
result** — the confidences are informative *and* 1/sd^2 is a reasonable functional form for
consuming them. The lane that proposed the hypothesis is the lane that refuted it, and recorded the
error as its own.

**2. L4's WITHDRAWAL STANDS; L4's stated REASON IS AMENDED.** Neither `shuffled` nor `shuf_paired`
isolates *where the errors land* — both move errors to random pairs. But the 0.537 A between them
does not price "permuting the weights is good". Read the shufr rows:

    shufr_wperm - shufr_wkeep = -0.467 [-0.598, -0.338]

`shufr_wkeep` hands pair p the residual r[pi[p]] with the weight sd[p] — residual and weight
**mismatched**. `shufr_wperm` hands it both r[pi[p]] and sd[pi[p]] — residual and weight **travelling
together**. So the difference prices **whether each residual keeps its own weight, and preserving
that alignment HELPS.** The correct sentence is *"orphaning a residual from its weight is bad"*, and
`wperm_only` confirms it on real data (+0.091, orphaning hurts slightly). L4's phrasing
"the weights are anti-informative" is superseded by this entry.

**3. THE COORDINATOR'S CORE HEADLINE IS REPRODUCED INDEPENDENTLY, on a different start.**
`shufr_wkeep` = **2.635** against a0 = 3.610, **-0.975 [-1.170, -0.783], 102W/24L** — the coordinator's
`shuffled` arm (2.609) rebuilt in another lane with the deployed weighting and only the residuals
permuted. **The distogram's errors are worse than random errors of the same magnitude** has now
survived every control thrown at it: separation-stratified, isotropic, weight-permuted,
weight-flattened, and an independent reimplementation.

**4. A FREE INTERNAL CONSISTENCY CHECK.** `a1_wperm` (perfect distances, permuted weights) = 1.162
against `a1_wkeep` = 1.152. **With perfect distances the weighting is irrelevant** — exactly as it
must be, and it validates the whole arm family.

**STILL RUNNING**: a tempering sweep, weight = sd^-p for p in {0, 0.5, 1, 1.5, 2, 3, 4} plus a
matched-random permuted-weight control, n = 126. Extended past p = 2 *because* `wflat` came out
worse: if the response is monotone increasing in p, the deployed exponent is not the optimum in the
direction nobody expected. Native-free and deployable if an interior or higher p wins.

---

## L11 — THE LAMBDA LADDER, n = 126: degree-1 is not equal to the full objective, it is DECISIVELY WORSE (2026-09-06, EXPERIMENT)

| lambda | 0 (degree-1) | 0.25 | 0.5 | 0.75 | 1 (full) | Control A (coordinate average) |
|---|---|---|---|---|---|---|
| mean RMSD | **4.335** | 3.742 | 3.683 | 3.619 | 3.610 | **3.048** |
| vs Control A | +1.287 [+1.124, +1.473] | +0.694 | +0.635 | +0.570 | +0.561 | -- |
| W/L | 18/108 | 23/103 | 27/99 | 32/94 | 31/95 | -- |

**Monotone DECREASING in lambda — the exact opposite of the pre-registered signature of a
mis-specified higher-order component.** Removing `E_ge2` removes information, not harm.
lambda=0 - lambda=1 = **+0.726 [+0.605, +0.816]**, median +0.663, 34W/92L. **No member of the
family beats the coordinate average.** The effect is uniform rather than concentrated (drop-top-10
sits at the 50-52nd percentile of a uniform-effect null) and holds in every length stratum and every
fold.

**Validity gate**: lambda=1 reproduces `s17/refine_full` on all 126 targets to mean |delta| =
**0.00001 A**; gradients verified against central differences of MATH's own callables (3.3e-10 rel).
EXPERIMENT consumed MATH's object throughout and retired its own provisional implementation.

**THE MECHANISM, from a zero-cost ORACLE diagnostic pre-registered before results.** Fed the
native's own distances, the full objective reaches **1.152 A** and degree-1 only **3.769 A** —
**+0.721 [+0.498, +0.945] worse than Control A.** Perfect distances buy the full objective 2.458 A
and degree-1 only 0.566 A.

> **The truncation destroys the structural information. No distogram improvement could rescue it.**
> This AMENDS the coordinator's "the objective is mis-aimed": the form is sound (L8), and reducing
> its expressive power is not a repair — it is the damage.

**Falsifiers, EXPERIMENT's verdicts**: **F1 FIRED** (4.335, not ~3.6) · **F2** does not fire only
because degree-1 is *decisively worse* than the full objective rather than equal · **F3 FIRED** (no
advantage to die; degree-1 does not beat a matched-magnitude random move at 4.190) · **F4 FIRED**
(sign reverses at n=126) · F5 is MATH's.

**PHASE 9 — the answer is *nowhere*.** Selection by either objective is worse than a random pool
pick (3.600 and 3.531 against 3.551).

**A HAZARD RAISED AND THEN RETIRED BY MEASUREMENT (EXACT), which CORRECTS L7/B5.** The residues
degree-1 cannot see are **exactly** the residues the metric cannot see — identical on **126/126 to
2e-14**. MATH's corollary that the objective is blind to `f_0` and `f_{n-1}` is true, but the
coordinator recorded it as "no arm can place the two termini the frozen metric scores", and **that
over-statement is withdrawn**: those torsion coordinates do not move any Ca, so the metric does not
score them either. Nothing is lost. **This is the second coordinator over-reading of a workstream
result in this sprint.**

**A DISSOCIATION INSIDE DEGREE-1.** MATH's mesh argmin is **not** the objective's optimum — L-BFGS
beats it on 85/126 — and reaching the true continuous optimum lowers the objective on **126/126**
while moving RMSD by **-0.016 [-0.057, +0.032]**. Objective quality and structural quality dissociate
*within* the truncated object too.

**H8 LEVERAGE (n = 74, PARTIAL) — a clean null WITH a mechanism, and the premise was backwards.**
`rho(leverage, |residual|) = -0.43`: the distogram's errors are *smaller* where geometric leverage
is higher. The leverage-weighting lead the coordinator issued (and retracted at L4) was built on an
assumption that measures the wrong sign.

**AND A RESULT THAT INDEPENDENTLY EXPLAINS L10.** `rho(1/sd^2, |residual|) = -0.476` — **the
distogram's confidences are rank-ordered in the correct direction.** This is an independent,
differently-derived explanation of why `wflat` lost by +0.153: the weights are informative because
the confidence ordering tracks the error magnitude correctly.

**Still PARTIAL and labelled as such**: Control B (n=72) and H8 (n=74), both holding sign and
reading from n~45 to n~74. Control B is a **null by construction** — greedy certifies a 1-opt
optimum on 59/72 within budget, 4x budget improves RMSD on neither objective, and degree-1's optimum
is closed-form at 106 calls, so **no optimiser at any budget changes any number here.**

---

## L12 — THE FOURTH ZERO-INFORMATION CONTROL MATCHES, and a methodological error worth more than the arm (2026-09-06, EXPERIMENT)

**1. The reference-measure control answers.** `exp_helixmu` (n = 14, PARTIAL) replaces the pool
torsion marginal with a constant ideal alpha-helix **as the reference measure the degree-1 object is
defined against** — MATH's estimator untouched, one added name to `mu_samples`.

| arm | pool mu | helix mu | difference |
|---|---|---|---|
| refine on `E_le1` | 4.346 | 4.306 | **-0.041 [-0.788, +0.607]** |
| certified argmin | 4.137 | 4.242 | **+0.106 [-0.424, +0.565]** |

**The zero-information measure matches the conditioned one.** Degree-1's information is **generic
backbone plausibility, not the pool's target conditioning.** This is the **fourth** zero-information
control in Sprint 18 to match its informative arm, and it lands inside the objective's own
definition rather than beside it. The falsification is untouched either way; this settles the
interpretation.

*(The coordinator launched the same module independently on a 30-target subset before this landed —
duplicated effort, the coordinator's scheduling error. The n = 30 run supersedes n = 14 when it
completes.)*

**2. A CORRECTION EXPERIMENT MADE AGAINST ITSELF, and it generalises beyond this sprint.**
EXPERIMENT had earlier reported the uniform-mu arm as weak evidence in the same direction. At n = 7
it was; **by n = 17 it reversed** (pool 3.918 vs uniform 4.465). The reading was withdrawn, and the
helix-mu result explains the contradiction:

> **Uniform-on-the-torus is not a valid zero-information measure.** It places mass on *impossible*
> backbone conformations, so it is a **worse** measure rather than an uninformative one.
> **Uniform is not the same as zero-information-but-plausible.**

This is the same class of error as reading "beats random" on an instrument where a constant helix
beats random — the programme's oldest recorded trap, met again in a new place. **Any future
zero-information control in this project must be plausible-but-uninformative, never uniform.**
Logged by EXPERIMENT as E13b.

**3. State frozen and self-consistent.** EXPERIMENT stopped its partial arms because they were
drifting away from the write-up as they accumulated, then regenerated findings, `exp_report.txt` and
the JSON artefacts once so all describe a single state: **complete at n = 126** for the ladder and
its controls; **PARTIAL at Control B n = 85, H8 n = 88, H6 n = 17, helix-mu n = 14**, each labelled
at its own n and each resumable. Control B and H8 held sign and reading unchanged from n ~ 45 through
n ~ 88; **H6 did not, and is flagged rather than smoothed.**

Control B remains a **null by construction**: greedy certifies a 1-opt optimum on 69/85 within
budget, 4x budget improves RMSD on neither objective, and degree-1's optimum is closed-form at 105
calls. **No optimiser at any budget changes a number in this branch.**

**Still open**: H6's magnitude, and helix-mu beyond n = 14.

---

## L13 — CORRECTS L12: the helix-mu control is NOT a match. It is NOT MEASURED (2026-09-06, coordinator's run; EXPERIMENT withdrew it first)

`s18/results/exp_helixmu.json`, **n = 30, `complete: true`**, paired against the same targets' pool-mu
arm in `exp_main_pool_512.json`.

| arm | helix mu | pool mu | helix - pool | W/L |
|---|---|---|---|---|
| refine on `E_le1` | 4.405 | 4.141 | **+0.264 [-0.341, +0.911]** | 13/17 |
| certified argmin | 4.551 | 4.105 | **+0.446 [-0.077, +1.053]** | 14/16 |

**The drift EXPERIMENT flagged is real, and it stabilises the other way from L12:**

    n=14   refine -0.041   argmin +0.106
    n=16   refine +0.154   argmin +0.331
    n=20   refine +0.256   argmin +0.493
    n=25   refine +0.219   argmin +0.334
    n=30   refine +0.264   argmin +0.446      <- COMPLETE

**L12's headline — "the zero-information measure matches the conditioned one", the sprint's fourth
zero-information control — is WITHDRAWN.** Both intervals span zero and both point estimates favour
the *conditioned* measure. The honest label is **NOT MEASURED**, not "matched": a zero-spanning CI
on a point estimate that moved monotonically from -0.041 to +0.264 is absence of evidence.

**EXPERIMENT caught this against itself before the coordinator did, one section after describing the
exact trap for another arm.** The coordinator then propagated the n = 14 reading into `LEDGER.md`
L12, `CLAIMS.md` D12 **and a project memory file** before the correction arrived. **That is the
third coordinator over-reading of a workstream result this sprint, and the first that escaped the
repository into durable memory.** The memory file is corrected.

**WHAT SURVIVES, and it is the properly-powered part.** The zero-information-control pattern rests
on the two arms measured at **n = 126**, not on this one:

  * a constant alpha-helix **START** matches the coordinate average at alpha=1: **-0.124
    [-0.366, +0.109]** (L8, ADVERSARIAL, n = 126)
  * a constant alpha-helix **PRIOR** matches the retrieval prior at every lambda: **-0.024
    [-0.138, +0.088]** at lam=10 (L9, coordinator, n = 126)

Those stand. The claim is therefore **two properly-powered zero-information controls matched their
informative arms**, plus the programme's pre-existing "a constant helix beats random in torsion
space" — **not four.** The reference-measure question inside the objective's own definition is
**OPEN at n = 30** and would need n = 126 to settle.

**THE DURABLE METHODOLOGICAL RULE, which is what this arm actually produced.** On this instrument a
small-n point estimate with a zero-spanning CI is **"not measured"**, never **"matched"**. A match
claim requires either a CI tight enough to exclude the effect size of interest, or n = 126. This is
now the second distinct way the sprint has been bitten by reading an underpowered arm as a null —
the first being uniform-vs-plausible reference measures (L12 item 2).

---

## L14 — TEMPERING THE DISTOGRAM'S WEIGHTS: the deployed exponent is at the optimum. No free win. (2026-09-06, ADVERSARIAL)

`s18/results/q_temper.json`, n = 126 on the main arms. weight = sd^-p, **every arm native-free**
(real dhat, real sd, no native distance read). sd dispersion: CV 0.623, min 0.336, max 4.284.

| arm | weight | RMSD | median | vs deployed p=2 | W/L | folds |
|---|---|---|---|---|---|---|
| p0.0 | sd^-0.0 (uniform) | 3.763 | 3.654 | +0.153 [+0.045, +0.272] | 50/76 | 5/5 |
| p0.5 | sd^-0.5 | 3.732 | 3.567 | +0.122 [+0.038, +0.219] | 47/79 | 5/5 |
| p1.0 | sd^-1.0 | 3.725 | 3.522 | +0.116 [+0.042, +0.202] | 47/79 | 5/5 |
| p1.5 | sd^-1.5 | 3.664 | 3.481 | +0.055 [+0.009, +0.115] | 44/82 | 5/5 |
| **p2.0 (DEPLOYED)** | sd^-2.0 | **3.610** | 3.370 | +0.000 | -- | -- |
| p3.0 | sd^-3.0 | 3.584 | 3.506 | **-0.026 [-0.108, +0.056]** | 66/60 | **3/5** |
| p4.0 | sd^-4.0 | 3.678 | 3.663 | +0.068 [-0.053, +0.209] | 61/65 | 4/5 |
| wperm | 1/sd^2 permuted | 3.752 | 3.606 | +0.142 [-0.005, +0.288] | 49/77 | 4/5 |

**1. THE DEPLOYED EXPONENT IS AT THE OPTIMUM AND THERE IS NO DEPLOYABLE WIN.** The response is
monotone in p up to the shipped value, then flat: p = 3 is **-0.026 [-0.108, +0.056]** — *inside the
instrument's 0.084 A MDE*, CI spanning zero, and **folds 3/5** rather than 5/5. p = 4 is worse.
**Nothing here justifies changing the exponent**, and the argmin-on-the-mean (p = 3) must not be
quoted as an improvement: it is a hyperparameter argmin on the tuning instrument, below the MDE,
without fold consistency.

**2. THIS IS THE THIRD INDEPENDENT CONFIRMATION THAT THE INCUMBENT'S WEIGHTING IS RIGHT**, and the
strongest form of it — a full curve rather than a point comparison. Deleting the weights costs
+0.153 [+0.045, +0.272] (L10's `wflat`, reproduced here at 5/5 folds), and every intermediate
tempering is monotonically worse than the shipped value. Together with EXPERIMENT's
`rho(1/sd^2, |residual|) = -0.476` (L11), the shipped confidence weighting is now the best-supported
native-free design decision in the programme.

**3. A LOGIC SLIP IN THE ARM'S OWN PRE-REGISTERED READING, flagged before it reaches a write-up.**
The pre-registration says: *"If uniform merely MATCHES the permuted weights, sd carries no usable
weighting information either way."* The measurement is uniform - permuted = **+0.011 [-0.120,
+0.134]** — they do match. **But the stated conclusion does not follow.** Uniform (+0.153) and
permuted (+0.142) are *both* worse than the deployed weighting; they agree with each other because
they are two ways of destroying the same information, not because there is no information to
destroy. The correct reading is that **sd carries real usable weighting information (destroying it
costs ~0.15 A) and the two null constructions are equivalent** — which is a clean consistency check
on the null, not evidence of absence. Reading 1 and Reading 2 are both refuted; the answer is
"the deployed form is right".

**STILL RUNNING**: an extension arm (`q_temper_ext`, 50/126 at time of writing).

---

## L15 — PHASE 6: `leg_contact` FAILS AS A TERM IN THE OBJECTIVE, exactly as it failed as a selector (2026-09-06, PHYSICS)

`s18/phys_lambda.py` -> `s18/results/lam.json`, **n = 126, `complete` = true**, config hash
`ae936e1ecaad4e10`, degree-1 source pinned at driver start and re-checked at every target.

**THE PRE-REGISTERED PRIMARY TEST FIRES ITS FALSIFIER.** `lam_c = +1` against `lam_c = 0`:

> **+0.108 A [-0.166, +0.377], median +0.199 (worse), 56W/70L.**

Interval includes zero, median positive, win/loss losing. **H-C is REFUTED.** Sprint 17's
`leg_contact` — the one term carrying **+0.080 [+0.032, +0.128]** of real in-band information the
distance model lacks — does not succeed as a term in the objective any more than it succeeded as a
selector. Branch closed: no ladder extension, no re-normalisation, no sign flip sold as a rescue.

**THE SHARPEST RESULT IS THE ZERO-INFORMATION CONTROL.** Against the identical objective with the MJ
table rebuilt on a **permutation of the residue labels** — functional form, pair set, switch and
magnitude distribution all preserved — the primary arm buys **-0.116 [-0.455, +0.164], 67W/59L**.

> **The MJ table's amino-acid identities are worth nothing measurable. What `leg_contact`
> contributes is the SHAPE of a pair term, not the sequence information in it.**

This is the sprint's **third** zero-information control to match its informative arm at n = 126, and
it lands on the Miyazawa-Jernigan matrix itself.

**THE PRE-REGISTERED SIGN WAS RIGHT AND NO RESCUE WAS EVER AVAILABLE.** `lam_c < 0` — the sign the
global correlation rho = -0.176 would have wanted — is catastrophic: **+1.314 [+0.958, +1.728],
27W/99L**. Declaring the physical sign before the run cost nothing and made a post-hoc flip
indefensible.

**ADDING `leg_contact` DESTROYS THE OBJECTIVE'S GLOBAL ORDERING, MONOTONICALLY IN lam_c**:
rho 0.351 -> 0.282 -> 0.124 -> **-0.060**, while in-band ordering stays flat and insignificant.
**The combined objective inherits `leg_contact`'s anti-ranking, not its information.**

**A FOURTH INDEPENDENT CONFIRMATION OF THE SPRINT-17 INSTRUMENT.** `leg_contact` alone reproduces
`s17/phys_FINDINGS.md` §3 **to three decimals on a rebuilt instrument**: rho global -0.176, rho
in-band +0.071, global argmin 5.634 A.

**INDEPENDENT CORROBORATION OF THE DEGREE-1 CLOSURE, from a fifth direction.** On PHYSICS's own
module with the contact machinery attached, `full` beats `lam_c = 0` by **-0.948 [-1.160, -0.753],
88W/38L**. PHYSICS states plainly that **the base of its own pre-registered primary arm is a closed
object.**

**L30 REPRODUCES.** `full` is +0.394 [+0.255, +0.523] worse than the projected start, and
3.606 - 3.048 = +0.558 against the coordinate average (against L30's +0.561). **Every refined arm in
the ladder is worse than the structure the pipeline already builds.**

**A DEFECT PHYSICS FOUND IN ITS OWN PRE-REGISTRATION, AND REFUSED TO REPAIR POST-HOC.** The ladder
matches the two terms on their **pool sd**, so `lam_c = 1` reads as "one sd against one sd". The
optimiser drives the terms over completely different ranges, so the **realised share** of the
combined objective's drop supplied by the contact term at the primary rung is **0.85-0.94**.

> **The pre-registered ladder never tested a SMALL contact correction.** Its smallest positive rung
> is already contact-dominated; what it measured is *"replace the objective with the contact term"*,
> not *"add a contact correction"*.

**And PHYSICS declined to re-normalise**, on the grounds that a scale chosen after seeing that the
first one gave an unfavourable answer is a tuned parameter wearing a methodological costume. The
realised share is printed beside every arm, and *"a correctly-scaled small contact perturbation is
untested"* is recorded as **OPEN** rather than quietly run. **That is the correct call and it is the
single best piece of methodological discipline in the sprint.**

---

## L16 — QUANTUM/ADVERSARIAL LANE CLOSED; and the `control_D` flag is BOUNDED TIGHTER THAN THE FLAG CLAIMED (2026-09-06)

`s18/quantum_FINDINGS.md` written, `s18/results/quantum_report.txt` regenerates every number, all
artefacts flagged complete. Leakage audit clean: no native quantity in any objective, projection,
labelling, weight, threshold or stopping rule; no bare `hash()` or `np.random` anywhere in
`s18/q_*.py` (the bootstrap was converted to `stable_rng` mid-sprint and every interval regenerated);
the sealed 60-target benchmark was not read, probed or derived.

**THE LANE FIXED MY LOGIC-SLIP CATCH AT THE SOURCE, NOT IN PROSE.** `q_temper.py` was *printing* the
invalid sentence beside its own table, so any future reader regenerating the report would have hit it
again. The module now prints the correction instead. **Fixing the generator rather than the document
is the right response and should be the standing pattern.**

**A PROCESS ERROR THE LANE CAUGHT ITSELF, and it nearly published the wrong conclusion.** Its first
`q_temper` launch ran p in {0..2}; a `pkill` did not take, the old process finished, and the rungs
past the deployed exponent were missing from the first report. It was caught by reading the printed
table against the module's own `PS` constant. **Had it not been, the curve would have been reported
as monotone with the shipped value at its edge** — precisely the wrong conclusion, in the direction
both the lane and the coordinator had flagged as the interesting one. `q_temper_ext` finished; p = 3
and p = 4 are on all 126 rows and nothing is PARTIAL.

### The `control_D` flag, and its real blast radius

The lane flagged that `math_anova`'s `control_D` self-check reports `PASS_machine_precision: false`
with **`E0_rel_err` 0.875-0.905, not converging in S**, and cautioned that *"anything that mixes
`E_le1` with another term at a fixed lambda IS affected"*.

**Measured, not argued** (`s18/math_iface.py`, target 1A13, n = 14, E0 = 74.82):

    E_res      E0 += 1000  ->  value shift +1000.0000   gradient max|delta| 0.000e+00
    E_lambda   E0 += 1000  ->  value shift  +500.0000   gradient max|delta| 0.000e+00

`E0` enters every objective **purely additively** (`float(t.E0 + val.sum())`) and contributes
**exactly zero** to the gradient. At a fixed lambda the mixture shifts by the constant
`(1 - lambda) * dE0`, which is still constant in theta.

> **Therefore NO argmin, NO optimisation trajectory and NO RMSD anywhere in this sprint is affected —
> the lambda ladder included.** The caution is correct in spirit and too broad as written: mixing at
> a fixed lambda is safe, because a lambda-dependent *constant* is still a constant.

**WHAT IS GENUINELY AFFECTED**, and where it must be stated: any **absolute objective value** quoted
from a truncated object (the error is ~0.89 x |E0|, about **67 objective units** on 1A13), and any
quantity that compares objective **magnitudes** across objects rather than **differences**. Every
difference-based quantity is safe by cancellation — which includes PHYSICS's `share` statistic (a
fraction of a *drop*) and every "objective fell X%" computed from two evaluations of the same object.

**Standing instruction for the report**: quote objective values from truncated objects as differences
or not at all, and do not compare an absolute `E_le1` magnitude with an absolute `E_full` magnitude.

---

## L17 — PHYSICS CLOSES: `leg_contact` is HARMFUL on the deployed objective, and every physics filter loses to a random gate (2026-09-06)

Five artefacts, **all `complete` at n = 126**. Instrument checks: G0 contact term = genuine Legacy to
**1.78e-15**; G1 gradient to 3.5e-10; G2b start **bit-identical** to the coordinator's (0.00e+00),
`full` 3.606 vs `objceil` 3.610; **L30 reproduces at +0.558 [+0.400, +0.680], 31W/95L against L30's
+0.561 and an *identical* 31/95**; AMBER k30 **+0.133 [+0.112, +0.165], 24W/99L**, matching s17 with
the **same three gate exclusions (`1D6X 2NB7 7BX2`) for the fifth time.** A mid-run pin
(`_pin_check`) prevented MATH's module landing mid-flight from mixing two degree-1 definitions inside
one artefact.

**1. ON THE DEPLOYED OBJECTIVE, `leg_contact` IS NOT MERELY UNSUPPORTED — IT IS HARMFUL.** The
declared extension gives **+0.722 A [+0.541, +0.893], 38W/88L**, interval excluding zero: **1.28 A
worse than doing nothing, losing on 109/126.** The primary rung's null (L15) becomes a signed
negative once the base is the objective actually deployed.

**2. THE COORDINATOR'S ORACLE GO/NO-GO IS ANSWERED, AND IT IS ZERO.** At the refinement start,
`corr(contact's descent direction, the correction the distogram's error needs)` = **-0.006, median
-0.012**, against its MJ-shuffled null **+0.007 [-0.014, +0.029]**, n = 126.

> **The ~1.0 A available from not trusting the distogram's error direction is UNREACHABLE through
> this term.** The go/no-go returned no-go with the interval tightened to +/-0.03.

PHYSICS also **corrected the statistic the coordinator specified**: the energy contribution and the
*force* have different supports — the force is exactly zero where the switch saturates — and they
disagree in sign at the native (+0.175 vs -0.149). The coordinator's formulation would have measured
the wrong object.

**3. THE HARDEST NEGATIVE IN THE SPRINT, AND IT IS NOT ABOUT ANYONE'S HYPOTHESIS.** Phase 8, on
identical structures: **every physics filter loses to a random gate of the same size.**

| gate | cost vs a RANDOM gate of the same size |
|---|---|
| AMBER single point | **+0.036 [+0.006, +0.061]** |
| `leg_torsion` | **+0.041 [+0.012, +0.068]** |
| `leg_contact` | **+0.053** |
| Legacy (total) | **+0.063 [+0.008, +0.125]** |

**Five scores, every CI excluding zero, every W/L losing.** Halving the candidate set costs
**+0.013 A**; halving it *by any physics score* costs **4-6x that**.

> **The cost is the ORDERING, not the truncation.** This is the most general statement the programme
> has produced about physics-based selection, and it holds for a force field and a knowledge-based
> potential alike.

**4. `leg_torsion` IS CLOSED BY ITS OWN WORKSTREAM.** Its near-native recall gain **reproduces**
(+0.083 [+0.019, +0.135]) and — as `s17` §8 demanded be tested — **it does not convert. It converts
negatively**: +0.055 worse than applying no filter at all. The one gate that raised recall without
damaging the set best does not turn that recall into accuracy.

**5. A RECORDED LAW FAILS UNDER A SELECTIVE INTERVENTION.** The operator law
`d_out = 1.16 * d_set_mean + 0.04 * d_set_best` predicts Legacy's gate should *help*. **The premise
held exactly** (set mean -0.024, set best +0.255) and **the output got +0.076 worse.** Miss
**+0.094 [+0.063, +0.132]** for Legacy and **+0.079 [+0.060, +0.098]** for `leg_torsion`, while the
law **holds for a random gate** (+0.006 [-0.003, +0.014]).

> **The law is descriptive over random subsets and breaks under score-based selection.** It must not
> be used to predict the effect of any gate that orders candidates. Project memory is corrected.

**6. PHASE 10 REFUTES THIS WORKSTREAM'S OWN s17 PREDICTION.** The spacing tax is a **displacement**
tax: a minimum-displacement correction restores 3.80 A spacing at *less* displacement than the
projection (0.801 vs 0.840) and **still costs +0.023 [+0.006, +0.037] more**; a *random* move of the
projection's magnitude costs the same as the projection (+0.009 [-0.014, +0.033]). **The +0.164 A
repair tax is irreducible, and the AMBER pass is never spent.**

**7. A SIXTH INDEPENDENT CONFIRMATION OF THE DEGREE-1 CLOSURE, at full scale.** `E_res` 4.265 vs
`E_full` 3.606 = **+0.658 [+0.483, +0.828], 38W/88L at n = 126**, where the 19-target instrument gave
-0.004. The gauge-invariant object that was null on the small instrument is **decisively worse** on
the real one.

**8. THE PRE-REGISTRATION DEFECT, reported and not repaired** (see L15): pool-sd matching made
`lam_c = 1` a 0.85-0.94 contact-dominated objective, so a *small* correction was never tested. Filed
**OPEN** rather than quietly re-run.

**ONE CHEAP QUESTION THIS LANE OPENED AND DID NOT CLOSE**: *why* score-ordering damages an averaged
set beyond its mean and its best — the survivors' **error covariance** under a score gate. Diversity
is consistent with it (Legacy 2.213 vs random 2.484 at the same m) but was not measured.
