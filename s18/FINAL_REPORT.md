# SPRINT 18 — FINAL REPORT
## The Degree-1 Objective Test

**126 cluster-disjoint targets · full-chain Cα-RMSD, frozen implementation · sealed 60-target
benchmark untouched · four workstreams · n = 126 unless stated**

---

# 0. EXECUTIVE VERDICT

> ## The degree-1 branch is closed. F1, F2, F3, F4 and F5 all fired, from six independent directions, in four independent lanes.
>
> ## The 19-target result that motivated this sprint was an artefact of which 2-bit code names which torsion state.

The brief said: *"If any F1–F5 fires: close the degree-1 branch"* and *"A clean falsification is a
successful sprint."* Every falsifier fired. This is that sprint.

**Nothing was rescued.** No tuned weights, no re-normalisation after an unfavourable answer, no
ladder extended past its pre-registration, no sign flipped post hoc, no easier metric substituted,
no benchmark consulted.

**The headline number.** The incumbent stands at **3.204 Å**; the best structure the pipeline builds
is the coordinate average at **3.048 Å**. The sprint did not move it. It measured, precisely, why
three plausible routes to moving it are closed — and it produced the programme's **first measured
native-free vindication of a deployed design choice.**

## 0.1 The six independent closures

| # | lane | argument | decisive number |
|---|---|---|---|
| 1 | ADVERSARIAL | **gauge** — W1 is not invariant under relabelling the torsion states | shipped labelling at gauge percentile **0.00**, p<0.0001; encoding luck **−0.440 [−0.844, −0.133]**, folds 5/5 |
| 2 | MATH | **wrong object** — the brief's bridge derives the residue-additive object, not W1 | RA vs full **−0.004 [−0.380, +0.307]**, 3W/5L/11 ties |
| 3 | MATH | **mechanism** — the degree-1 object *is* the retrieval torsion prior, and the continuum objective has no prior term | ρ(RA(hamil), torsion prior) = **0.986**, argmins 2.657 vs 2.615 |
| 4 | EXPERIMENT | **expressive power** — truncation destroys structural information | λ ladder **monotone decreasing**; degree-1 **4.335 Å**, +1.287 [+1.124, +1.473] |
| 5 | EXPERIMENT | **ORACLE mechanism** — perfect distances cannot rescue it | full **1.152 Å** vs degree-1 **3.769 Å** |
| 6 | PHYSICS | **full-scale reproduction on a fifth module** | `E_res` 4.265 vs `E_full` 3.606 = **+0.658 [+0.483, +0.828]**, 38W/88L |

The 19-target instrument gave −0.004 for the same object that the 126-target instrument prices at
**+0.658**. That reversal, on the gauge-invariant object, is the whole story of F4.

## 0.2 What actually survived, and it is short

1. **The distogram's errors are worse than random errors of the same magnitude.** Five distinct
   controls plus an independent reimplementation in a second lane.
2. **The shipped 1/sd² confidence weighting is validated, native-free.** Three independent
   confirmations. **The first deployed design choice this programme has ever confirmed rather than
   refuted.**
3. **Where an error improvement lands beats how large it is.** At *identical* residual RMS the
   outcome spans 2.145–2.697 Å.
4. **The cost of a physics filter is its ordering, not its truncation.** Every physics score loses
   to a random gate of the same size.

---

# 1. THE MATHEMATICS

## 1.1 The theorem (EXACT — verified, not asserted)

The residue-additive first-order functional ANOVA object equals the projection onto Walsh
coefficients supported **inside one residue** — weight 0, all weight 1, **and all intra-residue
weight 2**. Verified on 19 targets, raw and rank-uniformised:

| check | residual |
|---|---|
| RA ≡ within-residue Walsh projection | **8.9e-16** |
| qubit-ANOVA order 1 ≡ strict Walsh weight-≤1 | **1.1e-15** |
| orthogonality | 5.5e-17 |
| variance budget `Var(E) = Var(PE) + Var(E−PE)` | 5.0e-16 |

Retained variance: **W1 0.613, RA 0.913.** The strict truncation discards a third of the per-residue
field.

## 1.2 A coordinator error, carried into three workstreams

`s18/BRIEF.md` §4 asserted that under a uniform lattice measure the first-order ANOVA object **is**
the Walsh weight-≤1 projection, and instructed the sprint to verify it numerically. **It is false**;
the two objects differ by 1.4 objective sd. §4's own point 1 already said per-qubit ≠ per-residue and
the next sentence contradicted it. Corrected in place, labelled as the coordinator's error rather
than a workstream's. **The sprint spent effort porting the wrong object because of it.**

## 1.3 The gauge argument — the sprint's best single result

Relabelling which 2-bit code names which torsion state is **pure bookkeeping**: no structure, energy,
ordering or RMSD changes. RA is exactly invariant (sd 0.000 on 19/19). **W1 is not.**

| object | argmin RMSD (n=19) | vs full |
|---|---|---|
| full objective | 2.661 | — |
| **W1, the codebase's labelling** | **2.411** | −0.249 [−0.667, +0.147] |
| **W1, a random relabelling** | **2.852 / 2.862** | **+0.191** |
| RA (gauge-invariant) | 2.657 | −0.004 [−0.393, +0.317] |

Over the exact (S₄)ⁿ orbit the observed labelling sits at **percentile 0.00**, p < 0.0001; 23 of 24
relabellings are *worse* than the full objective. The encoding-luck term is **−0.440 [−0.844,
−0.133], folds 5/5** — **the only interval in the entire degree-1 story that excludes zero, and what
it measures is the arbitrariness of a bit encoding.**

Two lanes reached this independently (2.852 vs 2.862, RA invariant 19/19), by different routes.

The concrete failure: of the three ways to split 4 states into pairs, the encoding privileges two as
weight-1 and discards the third; only 8 of 24 permutations preserve that.

## 1.4 The mechanism — what the degree-1 object actually was

`hamil` is **75% an exactly-additive torsion prior** (prior term 1.0000 residue-additive, distogram
term 0.404, blend 0.9495), and

> **ρ(RA(hamil), the pure retrieval-pool torsion prior) = 0.986**, argmins 2.657 vs 2.615.

**The degree-1 object on the enumerated instrument is, to three nines, the torsion prior.** And
`s17/refine.py`'s continuum objective contains **no prior term at all** — so the truncation had
nothing to reduce to. F4 and F5 are not "it didn't transfer"; **the thing being truncated does not
exist in the continuum objective.**

## 1.5 The control that cuts the other way — quoted beside the closure, not under it

Zero-information 4.003 · ORACLE ceiling 1.030 · **matched-random projection at the same retained
variance 3.552**, which W1 beats by **−1.140 [−1.714, −0.561]**.

**Truncating by degree is genuinely not truncating at random.** The branch closes on *transfer* and
*gauge* — **not** on "degree structure is meaningless."

## 1.6 μ-sensitivity, kept as sensitivity rather than rescue

RA under uniform / pool / rama: −0.004 / −0.066 / −0.021, every interval covering zero, argmin
identical on 14/19. The largest estimate any μ gives is still below the 0.084 Å MDE. Pool and rama
**do not factorise over qubits on any target**, so W1 is not even an ANOVA under any informative μ.

## 1.7 An exact corollary, and the coordinator's over-reading of it

The objective is a function of the **n−2 interior residues**: `f_0 = f_{n−1} ≡ 0`. The coordinator
recorded this as *"no arm can place the two termini the frozen metric scores"* — **withdrawn.**
EXPERIMENT measured it: the residues the objective cannot see are **exactly** the residues the metric
cannot see, identical on **126/126 to 2e-14**. Those torsions move no Cα. **Nothing is lost.**

---

# 2. THE MAIN EXPERIMENT — the λ ladder at n = 126

| λ | 0 (degree-1) | 0.25 | 0.5 | 0.75 | 1 (full) | **Control A** |
|---|---|---|---|---|---|---|
| mean RMSD | **4.335** | 3.742 | 3.683 | 3.619 | 3.610 | **3.048** |
| vs Control A | **+1.287 [+1.124, +1.473]** | +0.694 | +0.635 | +0.570 | +0.561 | — |
| W/L | 18/108 | 23/103 | 27/99 | 32/94 | 31/95 | — |

**Monotone *decreasing* in λ — the exact opposite of the pre-registered signature of a mis-specified
higher-order component.** Removing `E_ge2` removes **information**, not harm. λ=0 − λ=1 = **+0.726
[+0.605, +0.816]**, 34W/92L. **No member of the family beats the coordinate average.** The effect is
uniform, not concentrated (drop-top-10 at the 50–52nd percentile of a uniform-effect null), and holds
in every length stratum and every fold.

**Validity gate**: λ=1 reproduces `s17/refine_full` on all 126 targets to mean |Δ| = **0.00001 Å**;
gradients match central differences of MATH's callables to 3.3e-10.

**The ORACLE mechanism, pre-registered before results.** Fed the native's own distances, the full
objective reaches **1.152 Å** and degree-1 only **3.769 Å** — worse than not refining at all. Perfect
distances buy the full objective 2.458 Å and degree-1 only 0.566 Å.

> **The truncation destroys the structural information. No distogram improvement could rescue it.**

**This amends the coordinator's own framing.** "The objective is mis-aimed" is wrong. The *form is
sound*; **reducing its expressive power is not a repair — it is the damage.**

**Phase 9 — the answer is *nowhere*.** Selection by either objective is worse than a random pool pick
(3.600, 3.531 vs 3.551).

**Control B is a null by construction.** Greedy certifies a 1-opt optimum on 69/85 within budget, 4×
budget improves RMSD on neither objective, degree-1's optimum is closed-form at 105 calls. **No
optimiser at any budget changes a number in this branch.**

---

# 3. THE OBJECTIVE'S ERROR STRUCTURE

## 3.1 The functional form is sound — and the ceiling is 0.083 Å, not 1.152

The parameterisation floor — projecting the **native itself** into the same ideal-geometry torsion
parameterisation — is **0.083 Å** (median 0.060, max 0.435). The parameterisation is not the
constraint anywhere.

| start at α=1 | start RMSD | end RMSD | median | objective at end |
|---|---|---|---|---|
| coordinate average (deployed) | 3.213 | 1.152 | 0.307 | **3.04** |
| member | 3.571 | 1.253 | 0.554 | 3.86 |
| **helix (zero-information)** | 4.070 | **1.028** | 0.412 | 3.68 |
| randtors | 5.258 | 2.080 | 2.278 | 18.38 |
| nativeproj (ORACLE) | 0.083 | **0.082** | 0.033 | **0.30** |

**The objective column decides it**: from the ORACLE start the fit reaches 0.30; from the deployed
start it stops at 3.04 — ten times higher, same objective, same weights. **The optimum is not at
1.152; the optimiser stops there.**

**But the median must be read beside it.** The paired median difference against the ORACLE start is
only **−0.094**. On more than half the targets the deployed start already reaches the native; **the
1.152 mean reports a minority that fall into a distant basin.** Anyone quoting 1.152 without that
sentence is quoting a minority.

**Three coordinator figures withdrawn here**: 1.152 Å as a ceiling; the α ladder as a *requirement
curve* (each rung mixes distance quality with basin reachability and cannot price the former alone);
and "halve the residual → 2.5 Å" as a point figure.

**A control the coordinator's arm did not carry**: a zero-information constant α-helix start reaches
**1.028**, indistinguishable from the coordinate average (−0.124 [−0.366, +0.109]). **The coordinate
average is not a privileged start.** This bears on every arm in the programme that begins from it.

## 3.2 The distogram's errors are worse than random errors of the same magnitude

The sprint's most robust result. All arms carry the **correct** weights; the real distogram is 3.610.

| arm | RMSD | vs the deployed objective |
|---|---|---|
| `shuffled` (residuals permuted) | **2.609** | −1.001 [−1.226, −0.784] |
| `shuf_strat` (permuted *within* separation bins) | **2.560** | −1.050 [−1.257, −0.860] |
| `isotropic` | **2.573** | −1.037 [−1.273, −0.802] |
| `shufr_wkeep` — **independent reimplementation, second lane** | **2.635** | −0.975 [−1.170, −0.783], 102W/24L |

Survived: separation-stratified permutation, isotropic replacement, weight permutation, weight
flattening, and an independent reimplementation on a different start.

**A retraction inside this result.** The coordinator's `shuf_paired` arm (2.072 Å) passed `sd[pi]`
into the fit, permuting the weights as well as the residuals — it was not the deployed functional.
Caught by the adversarial lane **reading the source rather than the table**, and retracted before it
reached the workstream it had been sent to as a lead. The correct reading of the `shuffled` /
`shuf_paired` gap is **residual–weight alignment**: `shufr_wperm − shufr_wkeep = −0.467 [−0.598,
−0.338]`, i.e. **orphaning a residual from its weight is bad.**

## 3.3 Where an improvement lands beats how large it is

A family of error models **all at identical residual RMS** (1.604, the "halve the residual" rung):

| error model | RMSD | vs uniform |
|---|---|---|
| **fix_confident** | **2.145** | **−0.308 [−0.462, −0.150]** |
| **fix_short** | **2.159** | **−0.294 [−0.424, −0.167]** |
| fix_big | 2.324 | −0.129 [−0.267, +0.004] |
| unif (the deployed ladder's rung) | 2.453 | +0.000 |
| **fix_long** | 2.618 | **+0.165 [+0.028, +0.303]** |
| **fix_unconfident** | **2.697** | **+0.244 [+0.097, +0.386]** |

**Spread 0.552 Å at identical residual RMS, straddling 2.5.** Residual RMS does not determine the
outcome. **The direction inverts engineering instinct**: a predictor that reduces RMS by fixing its
worst, longest-range, least-confident pairs lands at 2.62–2.70, not 2.45.

## 3.4 The shipped 1/sd² weighting is validated — three independent ways

**The first measured native-free vindication of a deployed design choice in this programme.**

| p | weight | RMSD | vs deployed p=2 | folds |
|---|---|---|---|---|
| 0.0 | uniform | 3.763 | **+0.153 [+0.045, +0.272]** | 5/5 |
| 0.5 | sd^−0.5 | 3.732 | +0.122 [+0.038, +0.219] | 5/5 |
| 1.0 | sd^−1.0 | 3.725 | +0.116 [+0.042, +0.202] | 5/5 |
| 1.5 | sd^−1.5 | 3.664 | +0.055 [+0.009, +0.115] | 5/5 |
| **2.0 (DEPLOYED)** | sd^−2.0 | **3.610** | — | — |
| 3.0 | sd^−3.0 | 3.584 | −0.026 [−0.108, +0.056] | **3/5** |
| 4.0 | sd^−4.0 | 3.678 | +0.068 [−0.053, +0.209] | 4/5 |

1. **Point comparison**: deleting the confidences costs **+0.153 [+0.045, +0.272]**, 5/5 folds.
2. **The full curve**: monotone below the deployed exponent with 5/5 folds at every such rung, flat
   above it. **p=3's −0.026 is below the MDE, spans zero, folds 3/5 — not quotable as an
   improvement**; it is an argmin-on-the-mean over the tuning instrument.
3. **A different diagnostic entirely**: ρ(1/sd², |residual|) = **−0.476** — the confidences are
   rank-ordered in the correct direction.

**Internal consistency check**: with perfect distances the weighting is irrelevant (`a1_wperm` 1.162
vs `a1_wkeep` 1.152), exactly as it must be.

---

# 4. THE PRIOR ARM — the term the mechanism said was missing

MATH showed the enumerated instrument's degree-1 object *is* the torsion prior and the continuum
objective lacks one. The coordinator added it back: distance + λ × a native-free per-residue von
Mises prior fitted to the retrieved windows. Reproduction gate passes exactly (`prior0` = 3.610).

| arm | RMSD | vs the coordinate average |
|---|---|---|
| coordinate average | **3.048** | — |
| prior λ=0 (= s17 refine_full) | 3.610 | +0.561 [+0.405, +0.715] |
| **prior λ=10 (best)** | **3.387** | **+0.339 [+0.241, +0.449]** |
| prioronly (no distance term) | 3.886 | +0.838 [+0.567, +1.122] |

**The primary hypothesis is refuted**: the best rung recovers 40% of the refinement penalty and still
ends **+0.339 worse than doing nothing.** No λ beats the coordinate average.

**The pre-registered falsifier fired.** At every λ the retrieval-pool torsion prior — **the best
torsion channel measured anywhere in this project** — is statistically indistinguishable from a
constant ideal α-helix (λ=10: −0.024 [−0.138, +0.088]), and at λ=0.03 it is significantly *worse*.
The 0.222 Å recovered is **regularisation**, reported as that smaller claim.

**One precise secondary result**: the `shuf` control (same (μ,κ) tuples permuted across residues) *is*
significantly worse than the true prior at large λ (−0.108 [−0.211, −0.006]).

> **A wrong positional assignment hurts; the right one is worth nothing over making no positional
> claim at all.** The prior's information is detectable and unusable in the same measurement.

**The coordinator's predicted mechanism was refuted by the coordinator's own arm**: displacement moved
only 5.034 → 4.960 rad (1.5%) while 0.222 Å was recovered.

---

# 5. LEGACY

## 5.1 `leg_contact` fails as a term exactly as it failed as a selector

| test | result |
|---|---|
| pre-registered primary (λ_c=+1 vs λ_c=0) | **+0.108 Å [−0.166, +0.377]**, median +0.199, 56W/70L — falsifier fires |
| **on the DEPLOYED objective** | **+0.722 Å [+0.541, +0.893]**, 38W/88L — **1.28 Å worse than doing nothing, losing 109/126** |
| λ_c = −1 (the sign global ρ would want) | **+1.314 [+0.958, +1.728]**, 27W/99L — catastrophic |

Not merely unsupported: **harmful, with a sign.**

## 5.2 The MJ table's amino-acid identities are worth nothing measurable

Against a zero-information control with the MJ table rebuilt on **permuted residue labels** —
functional form, pair set, switch and magnitude distribution all preserved — the primary arm buys
**−0.116 [−0.455, +0.164], 67W/59L.**

> **What `leg_contact` contributes is the *shape* of a pair term, not the sequence information in
> it.**

## 5.3 Adding it destroys the objective's ordering, monotonically

ρ global: 0.351 → 0.282 → 0.124 → **−0.060** as λ_c rises, while in-band ordering stays flat and
insignificant. **The combined objective inherits `leg_contact`'s anti-ranking, not its information.**

`leg_contact` alone reproduces Sprint 17 **to three decimals on a rebuilt instrument** (ρ global
−0.176, ρ in-band +0.071, argmin 5.634 Å) — a fourth independent confirmation of that instrument.

## 5.4 The ORACLE go/no-go returns zero

`corr(contact's descent direction, the correction the distogram's error needs)` = **−0.006**, median
−0.012, against an MJ-shuffled null of +0.007 [−0.014, +0.029], n=126.

> **The ≈1.0 Å available from not trusting the distogram's error direction is unreachable through
> this term.**

PHYSICS also corrected the coordinator's specification of this statistic: the energy contribution and
the **force** have different supports (the force is exactly zero where the switch saturates) and
disagree in sign at the native (+0.175 vs −0.149).

## 5.5 `leg_torsion` is closed by its own workstream

Its near-native recall gain **reproduces** (+0.083 [+0.019, +0.135]) and — as Sprint 17 §8 demanded
be tested — **it does not convert. It converts negatively** (+0.055 worse than no filter).

---

# 6. AMBER

## 6.1 The sprint's hardest negative: every physics filter loses to a random gate

Phase 8, identical structures, n = 126:

| gate | cost vs a **random** gate of the same size |
|---|---|
| AMBER single point | **+0.036 [+0.006, +0.061]** |
| `leg_torsion` | **+0.041 [+0.012, +0.068]** |
| `leg_contact` | **+0.053** |
| Legacy (total) | **+0.063 [+0.008, +0.125]** |

Five scores, **every CI excluding zero, every W/L losing.** Halving the candidate set costs
**+0.013 Å**; halving it *by any physics score* costs **4–6× that**.

> **The cost is the ORDERING, not the truncation** — and it holds for an all-atom force field and a
> knowledge-based potential alike.

## 6.2 The repair tax is irreducible and the AMBER pass is never spent

The spacing tax is a **displacement** tax. A minimum-displacement correction restores 3.80 Å spacing
at *less* displacement than the projection (0.801 vs 0.840) and **still costs +0.023 [+0.006, +0.037]
more**; a *random* move of the projection's magnitude costs the same as the projection (+0.009
[−0.014, +0.033]). **The +0.164 Å tax is irreducible.** This refutes PHYSICS's own Sprint 17
prediction.

## 6.3 AMBER's role is unchanged and reproduces exactly

AMBER k30: **+0.133 [+0.112, +0.165], 24W/99L**, matching Sprint 17 with the **same three gate
exclusions (`1D6X 2NB7 7BX2`) for the fifth time.** AMBER remains **stereochemical repair only**.

---

# 7. CVaR-VQE — the quantum pillar, closed honestly

**The answer was forced before the arms were run, and the pre-registration records that order.**

| object | weight-1 | inside 1 residue | **INTER-RESIDUE** | mean Pauli wt |
|---|---|---|---|---|
| full (deployed) | 0.6133 | 0.9132 | **0.0868** | 1.538 |
| RA | 0.6706 | 1.0000 | **0.0000** | 1.329 |
| W1 | 1.0000 | 1.0000 | **0.0000** | 1.000 |

The truncation **raises** weight-1 above the full objective's 0.613 and drives inter-residue coupling
variance to **exactly zero** — so the sprint's conditional quantum test **cannot fire on this
objective**. This was predicted by the coordinator in advance of the measurement.

Confirmed operationally: greedy 1-opt certifies the truncations' global optimum in **100% of cells at
36 evaluations** (one coordinate pass) against 63.2% on the full objective; both truncations' optima
are closed forms at ~38 table reads.

**The entangled vs CNOT-free arms were run anyway** — null on the consumed readout at every α, on RA,
W1 and full. Running them after the answer was forced, with that order pre-registered, is the
difference between a null and a null that could be accused of having been arranged.

**Q1, Q2 and Q3 all fired. The quantum branch is closed — on a genuine CVaR-VQE that was neither
replaced by a classical approximation nor asked to win.**

---

# 8. ARCHITECTURE DECISION

**Recommendation: keep the incumbent architecture unchanged, and stop spending on objective
engineering, physics ranking, and torsion-channel conditioning.**

| route | status after Sprint 18 |
|---|---|
| Re-engineer the distance objective's **form** | **Closed.** The form is sound — with true distances it reaches the 0.083 Å parameterisation floor. |
| **Truncate** the objective (degree-1, RA, W1) | **Closed**, six ways. |
| Add a **torsion prior** term | **Closed.** +0.339 worse than doing nothing; indistinguishable from a constant helix. |
| Add **`leg_contact`** as a term | **Closed.** Harmful with a sign (+0.722) on the deployed objective. |
| **Physics ranking / gating** (AMBER, Legacy, `leg_torsion`) | **Closed.** Every filter loses to a random gate of the same size. |
| **Search harder** (classical or quantum) | **Closed.** Null by construction; no optimiser at any budget changes a number. |
| **Re-weight** the objective (tempering) | **Closed.** The deployed exponent is at the optimum. |
| **Improve the distogram itself** | **The only route left open** — and §3.3 says *how*. |

**The one open direction, stated precisely.** Every closed route above is a way of *consuming* the
distogram better. The measurements say the consumption is already near-optimal: the weighting is
right, the functional form is right, the readout is right, and the optimiser is not the constraint.
What is wrong is **the distogram's errors themselves** — worse than random errors of the same
magnitude (§3.2), and unrescuable by any downstream operation.

**§3.3 is the actionable instruction and it inverts the natural instinct**: improve **confident,
short-range** pairs. Improving long-range or low-confidence pairs is *significantly worse than a
uniform improvement*. A predictor effort targeted the intuitive way lands at 2.62–2.70 Å; targeted
the measured way, 2.15 Å.

**Sub-2.5 Å assessment, honestly.** No arm in this sprint reached it, and the ORACLE ladder can no
longer be quoted as a requirement curve. What can be said: at identical residual RMS the achievable
band is 2.145–2.697 Å, so **2.5 Å is reachable on distogram improvement alone if and only if the
improvement is concentrated on the right pairs.** That is a claim about a predictor that does not
yet exist, and it is the honest statement of the gap.

---

# 9. BENCHMARK DECISION

> ## The sealed 60-target benchmark remains SEALED. Architecture-freeze conditions are not met.

It was not read, probed, derived from, or tuned against at any point in this sprint; the adversarial
lane's leakage audit confirms it. **No result here justifies unsealing**: the sprint moved no number,
and the one candidate that could have (the sd^−p tempering) came back saying the deployed value is
already optimal. Unsealing to confirm "we changed nothing" would spend a non-renewable resource for
no information — and the record already shows **no fresh benchmark exists**.

---

# 10. ERRORS AND CORRECTIONS — the coordinator's own, first

The brief requires errors be reported, not buried. **Four of the coordinator's calls were wrong and
every one was caught by a workstream rather than by the coordinator.**

| # | error | caught by | disposition |
|---|---|---|---|
| 1 | `BRIEF.md` §4's ANOVA/Walsh equivalence — **false**, and it sent three lanes after the wrong object | ADVERSARIAL + MATH | corrected in place, labelled as coordinator's |
| 2 | `shuf_paired` passed `sd[pi]` — not the deployed functional; the "error assignment" conclusion and the leverage lead built on it | ADVERSARIAL, reading the source | retracted before the receiving lane built on it |
| 3 | over-read MATH's `f_0 ≡ 0` corollary as "no arm can place the termini" | EXPERIMENT, by measurement (126/126 to 2e-14) | withdrawn |
| 4 | published a helix-μ "match" from **n=14** that drifted to the opposite sign by n=30 — **and it reached a durable memory file before correction** | EXPERIMENT, against itself | ledger, claims **and memory** corrected |
| 5 | mis-specified the ORACLE contact statistic (energy contribution vs force — different supports, opposite signs at the native) | PHYSICS | corrected in their implementation |

**Workstream errors, each self-caught and preserved in place**: ADVERSARIAL predicted `wflat` would
win and it lost (recorded as its own); a logic slip in its own pre-registration (*"if uniform matches
permuted, sd carries no information"* — false, both are worse than deployed) which it then **fixed in
the module that printed it, not just in prose**; a failed `pkill` that nearly published a truncated
tempering curve, caught by reading the table against the module's own constant. EXPERIMENT withdrew a
uniform-μ reading that reversed between n=7 and n=17, and killed a sibling's process assuming it was
its own stray. PHYSICS found its own pre-registration made the primary rung 0.85–0.94
contact-dominated so a *small* correction was never tested — **and declined to re-normalise**, on the
grounds that a scale chosen after seeing an unfavourable answer is a tuned parameter wearing a
methodological costume. **That is the single best piece of methodological discipline in the sprint,
and it cost that lane its headline.**

## 10.1 Two methodological rules this sprint produced

1. **A small-n point estimate with a zero-spanning CI is "not measured", never "matched."** A match
   claim needs a CI tight enough to exclude the effect size of interest, or n=126.
2. **A zero-information control must be plausible-but-uninformative, never uniform.** Uniform-on-the-
   torus places mass on impossible backbone conformations, making it a *worse* measure rather than an
   uninformative one.

Both are now in project memory, as is the correction that **the operator law `d_out = 1.16·d_set_mean
+ 0.04·d_set_best` breaks under score-based selection** — it holds for random gates (+0.006) and
points the wrong way for any gate that orders candidates (Legacy miss +0.094 [+0.063, +0.132]).

---

# 11. WHAT IS STILL OPEN

| # | question | why it matters | cost |
|---|---|---|---|
| 1 | **Why score-ordering damages an averaged set beyond its mean and its best** — the survivors' error covariance under a score gate | the mechanism behind §6.1, the sprint's hardest negative; diversity is consistent (Legacy 2.213 vs random 2.484) but unmeasured | cheap |
| 2 | A **correctly-scaled small** `leg_contact` perturbation | PHYSICS's ladder never tested it, and it declined to re-normalise post hoc | cheap |
| 3 | helix-μ at **n=126** | currently NOT MEASURED at n=30 (+0.264 [−0.341, +0.911]) | moderate |
| 4 | H6's magnitude | held at n=17 and did **not** hold its reading | cheap |

**A caveat that must travel with any reported objective value.** `math_anova`'s `control_D` reports
`E0_rel_err` 0.875–0.905, not converging in S. Measured, not argued: `E0` enters every objective
**purely additively** and contributes **exactly zero** to the gradient, so **no argmin, no
optimisation trajectory and no RMSD in this sprint is affected — the λ ladder included** (at fixed λ
the mixture shifts by a constant in θ). What *is* affected is any **absolute objective value** quoted
from a truncated object (~67 units on 1A13) and any comparison of objective **magnitudes** rather
than differences. Every difference-based quantity cancels.

---

# 12. INSTRUMENT INTEGRITY

| check | result |
|---|---|
| sealed benchmark | **not read, probed, derived or tuned against**; leakage audit clean |
| metric | full-chain Cα-RMSD, frozen implementation, unchanged; no alternative reported as primary |
| native information | evaluation and labelled ORACLE diagnostics only; no native quantity in any objective, projection, labelling, weight, threshold or stopping rule |
| seeding | `stable_rng` throughout; no bare `hash()` or `np.random` in any sprint module (one lane converted its bootstrap mid-sprint and **regenerated every interval**) |
| Legacy weights | `DEFAULT_WEIGHTS`, never fitted |
| AMBER discipline | convergence gate declared before use; same three exclusions (`1D6X 2NB7 7BX2`) for the fifth time; gated arms compared to their own gated input |
| cross-lane reproduction | contact term = genuine Legacy to **1.78e-15**; start bit-identical across lanes (**0.00e+00**); `full` 3.606 vs 3.610; L30 reproduces at +0.558 with an **identical 31/95** W/L; λ=1 reproduces s17 to **0.00001 Å** |
| completion flags | every quoted artefact `complete`; PARTIAL arms labelled at their own n and never quoted as results |

---

## Closing

The sprint asked whether a degree-1 truncation solves the objective-alignment problem. **It does not,
and the result that suggested it does was an artefact of an arbitrary bit encoding.** Six independent
closures, four lanes, every falsifier fired, nothing rescued.

What the sprint bought is smaller than a number and more useful than one: **the incumbent's
consumption of the distogram is now measured to be near-optimal in every dimension that was
suspected** — weighting, functional form, readout, optimiser, and prior — **and the remaining error
is in the distogram itself, with a measured instruction for where to attack it.**

*Ledger: `s18/LEDGER.md` (L1–L17). Claims register: `s18/CLAIMS.md` (91 rows, blocks A–H). Workstream
findings: `s18/{math,exp,quantum,phys}_FINDINGS.md`.*
