# SPRINT 16 — FINAL REPORT

**Programme**: quantum-assisted peptide structure prediction, 9–16 residues, full-chain Cα-RMSD
**Date**: 2026-09-06 · **Instrument**: 126 cluster-disjoint targets, five folds
**Sealed benchmark**: 60 targets, **untouched and unspent** — no module in this sprint reads it

> *This report is written to be checkable, not to persuade. Every positive claim carries its
> controls; every failed experiment is here; every retraction is preserved with the evidence that
> forced it. Where a number is a ceiling or an oracle diagnostic it says so in the same sentence.
> All nine workstreams have reported and are integrated.*

---

## 0. What a reader should take away

**Sprint 16 set out to build a steering mechanism, refuted it, refuted its own replacement, and then
refuted five of its own published explanations.** It produced **no method that beats the incumbent**,
its 9-target integration result did not survive the target as the unit, and its last candidate
accuracy effect — AMBER relaxation at −0.0234 Å — was **removed by its own final workstream**, which
showed that the baseline it beats is 0.155 Å worse than doing nothing and that a **matched-magnitude
random displacement is at least as accurate as AMBER's**. **Sprint 16 ends with no surviving
accuracy effect on any instrument.** What it produced instead is a much sharper and much less
flattering picture of where the programme's accuracy actually comes from, four closed directions
that were previously open and plausible, and a methodological result that may outlast the
scientific one.

The headline results, each with the correction that was applied to it:

1. **The flagship is refuted.** Every native-free steering arm chose "do nothing" on every fold
   (+0.000 [+0.000, +0.000], n = 126). Stepping the *exact* native torsion error half-way realises
   7.7% [−0.0, +14.9] of the achievable gain and makes 49 of 126 targets *worse* than not moving.
   **Corrected mechanism:** this is **not** a special property of torsion space — a
   zero-information control interpolating between two arbitrary retrieval windows is *more*
   non-monotonic. It is a **large-displacement regime**, and the displacement is large because the
   torsion-space fit is **87% of a uniformly random torsion assignment** (90° per-angle RMS).

2. **The coordinate-space replacement is refuted, and my first explanation of why was wrong.** The
   original claim — "from the best start every native-free direction has cosine indistinguishable
   from zero" — quoted the arithmetic mean of a *signed* cosine as evidence for a law that is
   **quadratic and two-sided** in it. Read correctly, the direction's **magnitude is 2–3× the
   isotropic null** (|c| 0.305 vs 0.128, W/L 100/26) and what is missing is the **per-target sign** —
   one bit — **and the step magnitude `‖r‖`, which is the RMSD to the native and which nothing in
   this sprint estimates native-free.** **A zero-information constant β-strand carries more of it (|c| 0.398) than any channel
   this sprint built**, with an ORACLE two-sided ceiling of 2.491 Å. That bit is **not obtainable
   native-free**: four references tested at n = 126 in Sprint 15 predict its sign at 0.357–0.413
   against a **constant-guess baseline of 0.516** — every one worse than guessing.

3. **Two of three inherited premises did not survive audit.** The native-free direction's signal is
   reproduced by a zero-information constant α-helix (the pool adds +0.033 [−0.007, +0.075]); the
   quiet-subspace identification is reproduced by a constant α-helix and by a random pool window;
   and the "parameter-free fusion law validated to 0.162 Å" is an **exact algebraic identity**
   (Krogh–Vedelsby, 1995) whose residual was 54% our own arithmetic-mean error and 46% our own
   frame-convention mismatch — **0% the identity failing**.

4. **Legacy contributes nothing to accuracy, and AMBER's accuracy contribution does not exist
   either — the −0.022 Å is the gap between two ways of making the structure worse.** On the shipped
   126-target ensemble Legacy's filtering effect (+0.024) is indistinguishable from dropping the
   same number of windows **at random** (+0.012), because the filter improves the set *mean* while
   destroying the set *best*. AMBER's relaxation effect against the projection is real as
   arithmetic — −0.0234 [−0.0375, −0.0099], n = 126, gated — but the **projection is 0.155 Å worse
   than not repairing at all** (27W/99L), AMBER is **0.133 Å worse** (24W/99L), and a
   **matched-magnitude random displacement is at least as accurate as AMBER's** (+0.0225 [−0.0138,
   +0.0571], 43W/83L). AMBER's displacement even points slightly **away** from the truth: ORACLE
   cos = −0.052 [−0.092, −0.013], and against the random control −0.049 [−0.100, −0.001] — an
   interval excluding zero in the wrong direction. **There is no (k, depth) at which the accuracy
   effect survives on the frame-reproducible subset**, and the frame null fails at every converged
   setting (max 0.137 Å, six times the claimed effect).
   **What is real is the repair.** AMBER reaches Ramachandran 0.874 and zero clashes **while moving
   only 26% of the residual** — but stated alone that is not evidence about a force field, because a
   **constant α-helix scores 1.000 and zero clashes by construction** and beats it 0W/65L. And **k
   is an accuracy knob, not a validity knob**: unrestrained AMBER reaches the same validity and costs
   +0.274 Å. **Legacy's earned role is detection** (AUROC 0.93–0.99 on a non-tautological
   component), and **fitting its weights makes it worse** (+0.78 Å [+0.04, +1.50]), which closes the
   standing escape that it merely lacked training.

5. **The energies' failure is a property of the potential class, not of short peptides.** On 19
   *exhaustively enumerated* spaces (1.28 × 10⁷ configurations, so every argmin is a **certified**
   optimum), a distance-restraint objective's argmin is **−0.926 Å [−1.440, −0.431]** better than a
   random feasible member, on 15 of 19 targets, with the native at the 29.1st percentile and a CI
   excluding chance — while **neither Legacy nor AMBER can be distinguished from random** on the
   same configurations. This is the sprint's only interval excluding zero on an argmin readout, and
   it is the control that makes "neither energy ranks the native" a statement about **contact-style
   potentials** rather than about 9–16-mers.

6. **The quantum component loses to uniform random sampling at matched budget** (+0.209 to +0.301;
   three of four arms with target-level CIs excluding zero, the fourth touching it; 7 of 9 targets
   worse on the largest arm), reconfirming the standing result on a terminal operator it had never
   been tested against. **A sibling workstream measured it against a stronger control and got a
   stronger answer**: the VQE is significantly worse than **greedy 1-opt with restarts at 10 of 10
   objective-quality rungs**, is never significantly better than simulated annealing at matched
   budget, and is indistinguishable from annealing at **a quarter** of its budget at 9 of 10 rungs
   (§5.8). α shows no monotone trend. **Corrected
   mechanism:** the penalty is **roughly half** lost ensemble diversity and half worse conformers
   (42–52% member error depending on the arm) — the original "it is diversity, not conformer
   quality" overstated it about twofold by reading the wrong superposition frame.

7. **The mechanism has no lever.** Diversity-aware selection on the real 126-target instrument is a
   clean null: **+0.015 [−0.019, +0.049], 61W/65L**. Raising the diversity weight lifted diversity
   22% and member error 4.5%, and the readout moved ≤ 0.008 Å — the two terms rise together and
   cancel, exactly as an identity with equal and opposite weights requires.

8. **One architectural claim survived everything**: reading an ensemble out in **coordinate** space
   rather than torsion space is worth **+0.446 [+0.193, +0.706]** at m = 75, on 8 of 9 targets, and
   the advantage grows with ensemble size.

**Against the sprint's RMSD goals**: the incumbent stands unchanged at **3.204 Å**, and **no arm in
Sprint 16 has a surviving accuracy effect on any instrument**. Nothing approaches 2.5 Å, let alone
2.0 Å, and no arm was allowed near the sealed benchmark because none earned it. The one number that
moved in the programme's favour is a **reference point rather than a method**: the raw coordinate
average, before any repair stage, is **3.0498 Å** — 0.155 Å better than what the pipeline ships —
which relocates the question from *how do we improve the prediction* to *why does making it a legal
peptide cost 0.155 Å, and can that be bought more cheaply.*

**A third reporting defect was found in the inherited record, by the workstream it damages.**
`core.project.lam_path` returns unwrapped torsions and the Ramachandran check compares them to fixed
windows; corrected, **AMBER's stereochemical advantage is +0.171, not +0.408 — the record overstates
it 2.39×**. Two previously unrecorded defects in the repair path were also found: the ideal-geometry
projection *degrades* Ramachandran quality (0.836 → 0.703), and AMBER **introduces cis peptide bonds
on 44 of 126 targets**.

**The methodological result may outlast the scientific one.** Across two independent audits, **every
number in both sprints reproduced exactly from its artefact** — and roughly a third of the
conclusions built on those numbers were still wrong, because a quantity was measured correctly and
then read as if it were a different quantity. Sprint 15's own verification pass had *passed* the
fusion-law row. **A transcription audit cannot see a formula error.** The two rules adopted in
consequence are in §7.

---

## 1. The problem, and where this sprint started

The programme predicts backbone structure for 9–16 residue peptides and scores full-chain Cα-RMSD
under a frozen definition (Cα only; every residue including both termini; no trimming; index
correspondence; uniform weights; proper rotations only; model 1 of the native; chain breaks scored,
not rejected). The shipped pipeline is:

    BLOSUM62 retrieval (K = 500) → distogram filter (top-75) → coordinate average
    → ideal-geometry projection (`ramah`, λ = 0.3) → AMBER relaxation      = 3.204 Å

Sprint 15 ran a full generative cascade against it and **failed**: F = 3.321 Å against 3.204,
+0.117 [+0.050, +0.188]. **Zero of thirteen** candidate improvements qualified. The protocol was
frozen with the incumbent unchanged and the 60-target benchmark deliberately left sealed by a
pre-registered rule, because nothing had earned the right to spend it.

Sprint 15 nevertheless left two quantities that had each been measured but never tested against
each other:

- a **native-free direction** `u = θ_pool − θ_fit` with |cos(u, e)| ≈ 0.39 against the true torsion
  error `e`, and
- a **Jacobian spectral split** of torsion space into a small "loud" subspace carrying nearly all
  the RMSD and a large "quiet" one carrying nearly none, with five exactly-null directions.

That pairing — a weak direction plus a way to discard the part of it that cannot matter — was the
largest unexplained ceiling in the programme, and it is what Sprint 16 was built to test.

---

## 2. Provenance: what was inherited, and what survived contact with an audit

Two workstreams re-examined the inheritance **before** the sprint's own results were read. Both
found the record numerically correct and the *interpretation* wrong. This is the single most
important methodological fact of the sprint and it is stated before any result.

### 2.1 Every recorded number reproduced; the comparisons were what failed

AUDIT re-derived Sprint 15's load-bearing claims from source. **0.390, 0.157, +0.499
[+0.326, +0.645], 0.827, 0.735, +0.980, 1.855, −1.822 [−2.037, −1.613] W/L 121/5, 3.677, 1.443 —
all reproduce exactly. There is no transcription error anywhere in the record.** What fails is what
those numbers were compared against:

| claim as recorded | what the audit found |
|---|---|
| the emitted structure's quiet subspace matches the native's at **0.827** | a **zero-information constant α-helix** reaches **0.829**, a random pool window **0.833**; the emitted structure is *significantly worse* at the quarter (−0.016 [−0.029, −0.003]). The 0.499 "null" is exactly k/m = 0.500, i.e. "half the space". Four of the five null directions are fixed coordinate axes. |
| an ORACLE ceiling of **1.855 Å** from rotating error into the quiet subspace | its own zero-information control — the same error norm along a **random** direction in the same quiet half — emits 2.164 Å. **Direction is worth −0.308 [−0.435, −0.172], 16.9% of the recorded effect.** The operation is also not "steer into the quiet subspace": it is a *loud* move of 52.2° RMS. |
| the native-free surrogate has \|cos\| **0.390** against a **0.157** null | a zero-information α-helix reference gives **0.357**; the pool adds **+0.033 [−0.007, +0.075]**, a CI including zero. The correct magnitude-matched null is **0.216, not 0.157**. |
| RMSD is blind to quiet motion, so quiet directions are free | **half-refuted.** RMSD is blind (ratio 0.0018–0.0024) but **the objective is not** (0.018–0.031): a 20° quiet step raises the objective by 53% of its value at the optimum. Ratio of ratios 10–17, > 1 on 99.2% of targets. |

**Consequence adopted for the rest of the sprint, and it should be adopted permanently:** every arm
must carry **both** a zero-information reference control and a matched random-direction control.
Without those two, this programme's evidence read roughly **five times stronger than it was**.

### 2.2 The fusion law is a theorem, and its "validation" was two of our own errors

RETRACT re-derived Sprint 15's "parameter-free fusion law `d_avg ≈ √(r² − (s/2)²)`, validated to
0.162 Å on 126 targets".

- It is the **Krogh–Vedelsby (1995) ambiguity decomposition** for two members — 31 years old,
  standard ensemble-learning material — and is **exact in any inner-product space**, with none of
  the caveats Sprint 15 attached to it.
- Sprint 15 computed `r` as an **arithmetic** mean where the identity requires the **quadratic**
  mean. Correcting that alone moves the residual +0.1619 → **+0.0748**, paired
  −0.0871 [−0.1185, −0.0605] fold-clustered [−0.1285, −0.0556], **W/L 126/0, 5/5 folds same sign**.
  Per target |r₁ − r₂| has mean 0.980 Å and 38% exceed 1 Å, so this is not a rounding matter.
- The remaining +0.0748 decomposes **exactly** into +0.1738 (averaging in the medoid frame rather
  than the target's) − 0.0990 (Kabsch-minimised `s`). In one common frame with the quadratic mean
  the residual is **2.65e−15**. The final superposition contributes exactly zero, provably.
- **So the 0.162 Å was 54% our arithmetic-mean error and 46% our frame-convention mismatch, and 0%
  the identity failing.**
- **Novelty surviving: none.** The proposed survivor "`s` is native-free so it works as an a-priori
  fusion screen" is **refuted**: `s` alone ranks "does fusion win here" at **AUC 0.401 against a
  0.500 null**, and a supervised leave-fold-out threshold degenerates to "always fuse" on 5/5 folds.

**A methodological finding that outlives the number.** Sprint 15's own verification pass **passed**
this row — every figure traced to its artefact exactly. **A transcription audit cannot see a formula
error.** Three hostile reviewer passes in Sprint 15 missed it for the same reason.

### 2.3 What the literature already owned

- **Native-free estimation of structural error direction is prior art, twice**: ATOMRefine (2023)
  predicts coordinate shifts between an initial model and the native and adds them directly;
  DeepAccNet (2021) predicts signed per-pair distance errors. **The number that should have set the
  sprint's prior: ATOMRefine buys GDT-HA 69.84 → 70.04.** Not the ORACLE ceiling.
- **The quiet subspace is prior art in substance**: manipulability ellipsoids, torsion-space normal
  mode analysis (iMod/iMODS), protein inverse-kinematics self-motion manifolds, KGS nullspace
  sampling, and concerted-rotation move sets built to minimise distant-atom displacement. The
  native-free half is the **Gauss–Newton assumption**; the mechanism is **Hansen's discrete Picard
  condition**.
- **Only the COMBINATION was unclaimed** across four literatures — and the combination is what this
  sprint refutes. The contribution is therefore a **clean refutation with a mechanism**, not a
  method.
- A prior-art correction in the programme's favour is also recorded: an AMBER threat LIT reported as
  reinstated was traced by RETRACT to a conflation of two sources, and Sprint 15's risk
  localisation was correct. Two natural defences of our AMBER position nonetheless fail and must
  not be used — their 1E0Q is a 17-residue β-hairpin, one residue above our band; and "free energy
  is a weaker observable than potential-energy ranking" is backwards, since it is Boltzmann-weighted
  and therefore stronger. **"AMBER anti-ranks peptide natives" moves from a finding to a
  confirmation of documented expected behaviour.**

---

## 3. The hypothesis, and what was fixed before any number was seen

**Hypothesis.** Projecting the weak native-free direction `u` onto the loud subspace removes the
part of it that cannot affect RMSD, raising its usable fraction and converting a ≈0.39 cosine into
an RMSD improvement.

**A falsifiable worry was stated in the module's own docstring before it was run** (ledger L2):
quiet directions are *defined* as barely moving the superposed Cα trace, so motion along them
should barely change RMSD by construction — which would make the Sprint 15 observation that
motivated the sprint a counterfactual requiring `θ_native` rather than a reachable move.

**A pre-registration was written while the n = 126 run was in flight** (`s16/PREREG_steer.md`),
fixing four decision rules and the follow-up diagnostics before any full-instrument number
returned. It also recorded, in advance, the three ways the n = 8 smoke was already hostile to the
hypothesis.

**What is native-free and what is not, stated exactly** (ledger L8). Native-free without
qualification: the fit target, the multi-start set, the direction `u`, the consensus direction, the
Jacobian and its spectral split. **Fold-honest supervised, not native-free**: the step size, chosen
by minimising RMSD on four training folds and applied to the held-out fold — the same convention
the programme already uses for its separation debias, but a *trained* quantity. ORACLE and labelled
as such everywhere: the error vector `e`, every `ORACLE_*` arm, the loud/quiet decomposition table,
and the curvature ratio.

---

## 4. The architecture under test

The preferred architecture, as specified:

    native-free error-direction estimation
      + Jacobian quiet-subspace steering
      + probabilistic torsion-space CVaR-VQE
      → ensemble generation
      → structure-level readout
      → Legacy auxiliary / detection branch
      → AMBER all-atom repair / refinement

Sprint 16 tested it as a deep investigation of one architecture rather than a shallow sweep, and
permitted an alternative to displace it only on evidence. One displacement occurred (§5.2) and the
evidence for it is given rather than asserted.

**The four mandated pillars, and where each is genuinely exercised:**

| pillar | where | what makes it genuine |
|---|---|---|
| VQE | `s16/integrate.py`, `α = 1.0` | real MPS ansatz, shot sampling, Adam, real qubit Hamiltonian, hard objective-evaluation budget |
| CVaR-VQE | same module, `α ∈ {0.1, 0.25, 0.5}`; boundary in `s16/qphase.py` | the same circuit optimised against a CVaR tail — arms of one ladder, not separate claims |
| Legacy / Miyazawa–Jernigan | `s16/integrate.py`, `s16/energy_ablate.py` | the genuine eleven-component model at `DEFAULT_WEIGHTS`, **never fitted** |
| AMBER ff14SB/GBn2 | `s16/integrate.py` (single points on demand), `s16/repair_*` (minimisation) | real force field and solvent model; the precomputed table covers 1.1% of the enumerated space and is deliberately not used |

**Legacy and AMBER are made explicitly comparable** by being applied as *the same operation* — drop
the worst quartile of an ensemble by that energy — so that the only difference between them is
which energy does the ranking. Comparing an energy against an *operation* (a filter against a
minimisation) would not compare the energies at all.

---

## 5. The controlled experiments

Full detail lives in `s16/STEER_FINDINGS.md`, `s16/CSTEER_FINDINGS.md`, `s16/INTEGRATE_FINDINGS.md`,
`s16/audit_FINDINGS.md`, `s16/lit_FINDINGS.md`, `s16/retract_FINDINGS.md`, `s16/verify_FINDINGS.md`,
with a one-row-per-claim register in `s16/CLAIMS.md` and a dated decision trail in `s16/LEDGER.md`.

### 5.1 The flagship — Jacobian quiet-subspace steering (`s16/steer.py`, n = 126)

**The mechanism's first step worked.** Restricting the native-free direction to the loud subspace
raised its alignment with the loud part of the true error from **0.318 to 0.536** at rank 4 and
**0.712** at rank 2. The projection did exactly what it was designed to do.

**And it bought nothing.** Every native-free arm — `full`, `loud4`, `loud10`, `loud13`, `consensus`
— chose step 0 on every one of five folds: **+0.000 [+0.000, +0.000]** against the plain fit. Since
`s = 0` is on the ladder, a training fold can always decline to move, so this is the selector
saying, with labels, that the direction is not worth using.

Two pre-registered rules fired:

- **Rule 2** — `loud*` does not beat `full` and both are ≥ 0: *hypothesis refuted at the mechanism
  level.*
- **Rule 4** — `ORACLE_loud4` bought **−0.028 [−0.106, +0.051]** against a rank-4 first-order RMSD
  share of **0.875**, while the full true error bought −3.604. Therefore **no rank-r first-order
  share may be quoted as a finite-step budget** anywhere in this sprint or in the Sprint 15
  provenance sections.

Rule 3 did **not** fire, and that matters: the n = 8 smoke's `quiet_half` anomaly (−0.106 with a CI
**excluding zero**) evaporated at n = 126 (−0.003 [−0.011, +0.007]). This is the sixth small-n
reversal in the programme and the first with a small-n interval that excluded zero. **A confidence
interval at n = 8 is not protection.**

### 5.2 Why — the corrected mechanism

Stepping the **exact** native torsion error, so that `s = 1` lands on the native by construction:

| step s | 0.0 | 0.1 | 0.2 | 0.35 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|---|---|---|
| mean RMSD (Å) | 3.647 | 3.599 | 3.641 | 3.663 | **3.466** | 2.596 | 0.043 |
| per-target fraction of gain | 0 | 1.6% | 0.8% | 1.1% | **7.7% [−0.0, +14.9]** | 32.3% | 100% |
| improved / worsened | — | 82/44 | 75/51 | 73/53 | **77/49** | 104/22 | 126/0 |

Half-way to the correct answer in the correct direction is worth about 8% of the answer, median
16.1%, and **makes 49 of 126 targets worse than not moving at all**. A linear map would give 50% on
every target.

**The first explanation published for this was wrong and is retracted.** The claim that the fit's
torsion errors are *compensating* had no control. VERIFY built one — geodesic interpolation between
two arbitrary retrieval windows, native-free, matched on ‖Δθ‖ — and it bulges **more** than the fit
(−0.024 against +0.014 at s = 0.5). The effect tracks ‖Δθ‖ at Spearman **−0.571**, and **below about
4 rad both are linear**.

**The correct account is worse for the programme.** This is a large-displacement regime that would
afflict any torsion parameterisation, and the displacement is large because **the torsion-space
distance-geometry fit is 87% of a uniformly random torsion assignment** (90° per-angle RMS). The
curvature ratio, ‖J e‖/√n = 9.879 Å against a realised 3.647 Å — **2.709** — is the same fact in the
Jacobian's language.

One objection was tested and **refuted**: `A.wrap` *is* the exact per-angle geodesic on the torus
(agreeing with an independent great-circle construction to 1.8e−15 rad), so the non-monotonicity is
not an artefact of interpolating wrapped angles.

### 5.3 The pivot — coordinate-space steering (`s16/csteer.py`, n = 126)

The flagship's mechanism named its successor. In coordinate space the payoff is closed-form:
`t* = c‖r‖` gives `‖r_new‖ = ‖r‖√(1 − c²)`, which turns the sprint's target into arithmetic — from
3.204 Å, **2.5 Å needs c = 0.625 and 2.0 Å needs c = 0.781**. (`s16/csteer.py` hardcoded 0.615, which
is the value for a 3.170 Å base; it was corrected at source on 2026-09-06.)

**The operational answer is a null**, and it is unambiguous: with a leave-fold-out step and the
channels available, nothing is deployable and 2.5 Å is not approached.

**The published explanation of that null was wrong and is retracted.** The law is quadratic and
two-sided in `c`, so the functional it consumes is **|c|**; the findings quoted the arithmetic mean
of the *signed* `c` and concluded the direction was absent. Corrected, from the best start
(the pool coordinate average, 3.048 Å), against the module's own random control, whose sampled |c|
is 0.136 and whose analytic isotropic null √(2/(π·3n)) at n = 13 is 0.128 — the paired subtractions
below use the 0.128 control mean:

| direction | signed c (as published) | **\|c\|** | **\|c\| − random** | W/L | ORACLE two-sided ceiling |
|---|---|---|---|---|---|
| toward `fit` | +0.012 | **0.305** | **+0.177 [+0.142, +0.212]** | **100/26** | 2.843 Å |
| **constant ideal β-strand (ZERO-INFO)** | −0.012 | **0.398** | **+0.262 [+0.216, +0.310]** | — | **2.491 Å** |
| pure isotropic breathing (ZERO-INFO) | −0.067 | 0.296 | +0.167 [+0.123, +0.214] | — | 2.708 Å |
| random control | −0.016 | 0.136 | +0.000 | — | 3.009 Å |

**A constant ideal β-strand — which knows only the chain length — carries more usable direction than
any channel this sprint built.** The residual `avg → native` is dominated by a **compactness
scalar**: the magnitude is real and the **sign flips per target** (58% need contraction, 42%
expansion).

Two further defects in the same module, both recorded:

- **The ladder had no negative rung.** The fraction of targets whose optimal step is negative is
  0.45–0.67 depending on the channel, so on roughly half the instrument no arm could have helped
  whatever the direction carried. The published null is **partly a property of the ladder**.
- **Two arms that looked like large improvements (−0.599 and −0.164 Å) were endpoint substitutions.**
  Both chose step 1.0 on every fold — they walk the whole way to the pool coordinate average, a
  structure the incumbent already computes. Caught internally, before VERIFY. **Rule adopted: a
  leave-fold-out step selector that lands on the boundary of its ladder is choosing an endpoint, not
  a step.**

**Is the missing bit obtainable native-free? No — and the answer already existed.** Rather than build
a new module, the record was checked: Sprint 15's `s15/expand.py` tested four native-free scale
references at n = 126 with the closed-form optimal isotropic scale. Read at the level of the sign:

| native-free reference | sign accuracy vs the RMSD-optimal scale | correlation |
|---|---|---|
| distogram | **0.357** [0.278, 0.444] | +0.021 |
| separation-debiased | 0.389 [0.302, 0.476] | +0.011 |
| pool-matched | 0.413 [0.325, 0.500] | −0.037 |
| both | 0.389 [0.310, 0.476] | −0.019 |
| **zero-information constant-sign baseline** | **0.516** | — |

**All four are worse than always guessing the majority sign**, and Sprint 15 measured the cost of
using them: +0.080 to +0.123 Å, every interval excluding zero.

A sharper fact falls out of the same table. Against the scale that matches the *predicted distances*
the same references reach 0.54–0.61 — above chance. Against the **RMSD-optimal** scale they are
below a constant guess. **Matching the predicted distance scale is not the same as minimising RMSD,
and the two disagree in sign.**

### 5.4 The integrated architecture, and Legacy against AMBER inside it (`s16/integrate.py`)

The enumerated instrument: 9 exhaustively enumerated targets, 3 seeds, budget 6144 objective
evaluations. Legacy and AMBER are made comparable by being applied as **the same operation** — drop
the worst quartile of an ensemble by that energy — so the only difference is which energy ranks.
AMBER is genuine ff14SB/GBn2 computed **on demand**, because the precomputed table covers 1.1% of
the space and a sampled ensemble essentially never intersects it.

| contrast (m = 75) | mean | target-level CI | median | W/L |
|---|---|---|---|---|
| AMBER − random drop | −0.054 | **[−0.148, +0.040]** | −0.045 | **4 targets worse / 5 better** (101/61 rows) |
| Legacy − random drop | +0.005 | [−0.034, +0.044]\* | −0.025 | 85/77 rows |
| Legacy − random drop, m = 5 | +0.068 | **[−0.039, +0.176]** | +0.036 | 70/90 rows |
| AMBER − Legacy | −0.059 | direction only | −0.043 | 97/65 |
| Legacy→AMBER − random | −0.053 | direction only | −0.051 | 99/63 |

\* row-level intervals; the target-level recomputation widens them.

- **Legacy contributes nothing** as a ranker in this architecture; the point estimate on small
  ensembles is worse than dropping members at random, but its target-level interval
  (+0.068 [−0.039, +0.176]) contains zero, so it is not a significant harm either.
- **AMBER ranks better than Legacy on the filtering readouts** — but not everywhere: on the 19
  exhaustively enumerated spaces the direction reverses and neither energy is distinguishable from
  the other (§5.7, final bullet). And AMBER's own contribution as a ranking filter here
  **does not survive the target as the unit**: the interval widens 2.79×, 4 of 9 targets are worse,
  dropping any of 8 of the 9 destroys it, and seed 2 alone gives −0.0003.
- **The composition adds nothing over AMBER alone** (interaction −0.0044 [−0.1040, +0.0966],
  target level).
- The tie trap that once cost this programme a whole table was censused here: **zero ties** in
  either energy.

**Both energies are dwarfed by the choice of generator**, which moves 0.30 Å.

### 5.5 What the quantum component did, and why (`s16/integrate.py`, `s16/diversity.py`)

At m = 75, against the two mandatory controls at matched budget — positive means the VQE is worse:

| arm | vs untrained-circuit best-of-N (target level) | vs uniform random (target level) |
|---|---|---|
| `cvar0.1` | +0.115 (row-level CI only) | **+0.212 [+0.003, +0.404]** |
| `cvar0.25` | **+0.204 [+0.028, +0.384]** | **+0.301 [+0.115, +0.485]**, 7 targets worse / 2 better |
| `cvar0.5` | +0.137 (row-level CI only) | **+0.234 [+0.011, +0.458]** |
| `cvar1.0` (plain VQE) | +0.112 (row-level CI only) | +0.209 **[−0.002, +0.410]** |

**Three of the four CVaR-VQE arms lose to uniform random sampling at the same budget with a
target-level interval excluding zero; plain VQE's interval touches zero.** These contrasts are
paired within (target, seed), so the row bootstrap was *nearly* right here — the widening that
destroyed the AMBER interval is small for them, but it is not zero, and the target-level intervals
are the ones quoted. This reconfirms the standing project
result — *running the VQE is worse than not running it* — on a **terminal operator it had never been
tested against**: the earlier finding was measured through an argmin, this one through an ensemble →
coordinate-average readout. **α does not rescue it**: no monotone trend, and plain VQE is among the
better arms.

**Two qualifications, both from the sibling QPHASE workstream (§5.8), that this table needs.** First,
**uniform random is a weak control** — QPHASE labels it "proves nothing" and runs the strong ones,
where the result is larger and cleaner. Second, **the sign of the untrained-circuit contrast is
readout-dependent**: QPHASE's argmin readout on the real objective gives the VQE *better* by −0.309
(null), where this 75-member coordinate average gives it worse by +0.112. That is consistent with
the diversity account below — the two readouts consume different quantities — but it means the
untrained-circuit column above must not be read as a general statement. **The general statement is
§5.8's**: against greedy 1-opt and against annealing, at matched budget, the VQE loses on every
readout and at every rung.

**Why.** The terminal operator obeys the Krogh–Vedelsby ambiguity decomposition exactly —
‖avg − native‖² = mean member error² − diversity² — verified through this operator at machine
precision (max relative 5.5e−15; it is the parallel-axis theorem, and is labelled as a theorem
check, not a discovery). Because the decomposition is exact, attributing each generator's readout
penalty between its two terms is not a hypothesis but arithmetic:

| generator | share from member error | share from lost diversity |
|---|---|---|
| `cvar0.25` | 52% | 48% |
| `cvar1.0` | 51% | 49% |
| `cvar0.5` | 42% | 58% |

**The VQE's penalty is roughly half concentration and half worse conformers.** The first published
version of this — "they differ in diversity, not in conformer quality" — read the
free-superposition member error (spread 0.056 Å) where the identity consumes the common-frame one
(spread 0.122 Å), and **overstated the diversity share about twofold**. The recorded project law
that *the terminal operator consumes the set mean* is therefore **incomplete rather than
overturned**: it consumes the set mean *and* the set diversity, at equal and opposite weights.

**The one architectural claim that survived everything.** Reading the same ensemble out in
coordinate space rather than torsion space is worth **+0.214 / +0.297 / +0.446** at m = 5 / 20 / 75,
**8 of 9 targets** at target level (**+0.446 [+0.193, +0.706]**), and **the advantage grows with
ensemble size** — which is what the decomposition predicts, since more members means more diversity
to cancel and torsion averaging cannot cash it.

### 5.6 The one constructive experiment the mechanism suggested (`s16/divselect.py`, n = 126)

If the operator's gain includes a diversity term with equal and opposite weight, the incumbent —
which selects its 75 windows by distogram score alone — optimises one term and ignores the other.
The falsifier was written before the run.

**Instrument check**: the greedy selector at λ = 0 reproduces the shipped top-75 exactly (3.048 Å,
+0.000), so arms differ only in selection.

**The trade is real, measurable, and cancels.** As λ rises at m = 75, diversity rises 1.902 → 2.316
(+22%) and member error rises 3.615 → 3.778 (+4.5%), and **the readout moves by at most 0.008 Å**.

**The deployable arm is a clean null**: **+0.015 [−0.019, +0.049]**, median +0.000, **61W/65L**. The
module's own trap flag fired — three of five folds pinned λ at the ladder end, the second
boundary-pinning instance in the sprint.

**One narrow tie, inseparable from its caveat — and it is a tie, not a positive.** Selecting 75
windows by **pure diversity, ignoring the distogram score entirely**, gives 3.051 Å against the
shipped 3.048 Å — very slightly *worse*, with no interval computed — while matched random gives
3.090 Å. But the candidate set is the top-200 *by score*, so this is **not** a claim
that the score is worthless — it locates the score's selective work in the 500 → 200 cut rather than
the 200 → 75 one, consistent with the recorded result that top-24 filtering caps any reranker.


### 5.7 The causal ablation on the shipped ensemble — Legacy against AMBER at n = 126

`s16/energy_ablate.py`. One **fixed** ensemble (the shipped top-75 windows, a pure cache read with
no RNG), a 2×2 factorial plus a **matched random-filter control**, every AMBER row reported both
gated and ungated with the excluded count printed. This is the full-instrument counterpart of §5.4,
and it is the more trustworthy of the two.

| contrast | gated mean | 95% CI | median | W/L | n |
|---|---|---|---|---|---|
| **AMBER** (`amb − ctrl`) | **−0.0234** | **[−0.0375, −0.0099]** | −0.0074 | 66/57 | 123 |
| Legacy (`leg − ctrl`) | +0.0243 | [−0.0120, +0.0615] | +0.0034 | 59/67 | 126 |
| **Legacy vs its matched random filter** | +0.0118 | [−0.0242, +0.0492] | +0.0024 | 60/66 | 126 |
| random filter alone (`rnd − ctrl`) | +0.0124 | [−0.0022, +0.0291] | — | 57/69 | 126 |
| interaction | +0.0073 | [−0.0073, +0.0218] | +0.0020 | 57/65 | 122 |

**Legacy contributes nothing to accuracy, and the reason is the set-mean trap.** Its entire effect
(+0.024) is indistinguishable from the cost of dropping 19 of 75 windows **at random** (+0.012). The
filter *improves* the ORACLE set **mean** (3.551 → 3.531) while *destroying* the ORACLE set
**best** (2.306 → 2.416) and contracting diversity by 6% — which is the same identity from §5.5 acting
against the pipeline. This is the trap the ledger has warned about, caught in the act.

**AMBER's accuracy effect exists at aggregate and is not resolvable per target.** −0.0234
[−0.0375, −0.0099] with the convergence gate applied (3/126 excluded canonical, 4/126 rotated) —
and the workstream's own **pre-declared frame null FAILS**: +0.0146 Å mean, max 0.140 Å on 11
converged targets. **The band was not loosened.** The consequence is stated exactly: *no per-target
AMBER RMSD claim is resolvable*, while the aggregate clears the 123-target floor by about 10 SE. The
median-vs-mean early warning fired on this row and the null-calibrated concentration check returned
**PASS at the 58.5th percentile**, so the effect is not one or two targets carrying the mean.

**Legacy's real role is detection, and it is good at it.** AUROC **0.93–0.99** on the `torsion`
component — and that one is *not* tautological — and 0.98–0.99 on `steric`, which partly is. Legacy's
only intervals excluding zero anywhere in the sprint are on the validity axis: **+0.032
Ramachandran-favoured and +0.069 Å minimum separation** over its matched random control. That is
precisely the "auxiliary / detection branch" role the architecture assigns it, and it is the role it
earns.

**Three findings that damage the energies further:**

- **The eleven-term Legacy total is worse than its own best component on all ten axes.** MJ contact
  is at chance; `compactness` is an *anti*-detector; two components are identically zero on all 126
  inputs. Removing MJ changes the argmin defect hardly at all (+0.082 → +0.090), so the defect is
  **distributed**, not one bad term — but MJ alone is +0.375 Å worse than random.
- **Fitting Legacy's weights makes it worse**, decisively: leave-fold-out fitting moves the certified
  optimum **+0.7785 Å [+0.0404, +1.4975]** in the wrong direction (8W/11L) on weight vectors that
  gain substantially *in sample* (ρ 0.105 → 0.213). The workstream's own hypothesis — that fitting
  would rescue Legacy — is **refuted**, and with it the standing escape that Legacy underperforms
  only because it was never trained.
- **On 19 exhaustively enumerated spaces, neither energy ranks, and they are indistinguishable from
  each other.** On the identical binding-masked population Legacy's certified argmin is −0.184 Å
  [−0.523, +0.184] against a random feasible member and AMBER's is −0.075 Å [−0.540, +0.358]; head
  to head, Legacy − AMBER = −0.109 Å [−0.666, +0.462]. **Every one of those intervals contains
  zero**, and the ORACLE native percentiles (Legacy 0.531, AMBER 0.407) do not exclude chance
  either. On the *full* enumeration Legacy's argmin is +0.082 Å [−0.259, +0.401] **worse** than
  random, reproducing Roget et al., to whom it must be attributed. **The one arm that breaks the
  pattern is not an energy**: the distance-restraint objective's certified argmin is **−0.926 Å
  [−1.440, −0.431]** better than random (15W/4L, native at the 29.1st percentile, CI excluding
  chance), beating AMBER by −0.815 [−1.357, −0.310] on the identical masked set. That localises the
  pathology to the **potential class**, not to short-chain structure prediction. A zero-information
  constant α-helix reaches an RMSD 0.469 Å below Legacy's certified optimum ([−0.006, +1.041],
  5W/6L with **8 exact ties** — on 8 of 19 targets Legacy's optimum *is* the constant helix), which
  is why "beats uniform random" proves nothing on this instrument.

**Two new defects in the repair path itself**, neither previously recorded:

- **The ideal-geometry projection *degrades* Ramachandran quality** (0.836 → 0.703) and minimum
  separation. AMBER repairs both — but **introduces cis peptide bonds on 44 of 126 targets.**
- Cost at the right number: **13.88 s** for the minimisation against **5.01 s** for the projection,
  about 1,500× an AMBER single point.

**And a correction that runs against this workstream's own case.** A **third** Sprint-14/15 reporting
defect: `core.project.lam_path` returns *unwrapped* torsions (values as extreme as +1022.72°) and
`rama_ok` compares them against fixed windows. Corrected, the baseline arm moves 0.4658 → 0.7029
(+0.2371 [+0.1946, +0.2796], 0W/85L, 85 of 126 targets affected). **AMBER's Ramachandran advantage is
+0.171, not +0.408 — the record overstates it by 2.39×.** Six shared-document lines need the fix.

**One free result worth recording**: all 126 control-arm targets reproduce
`s15/results/phys_repl_canon0.json` at **max |Δ| = 0.000e+00 Å** — across sprints, interpreters, a
code change to `core/amber.py`, and four worker processes.


### 5.8 The CVaR-VQE phase boundary — and the finding that there is no boundary

`s16/qphase.py`, `s16/qphase_FINDINGS.md`. **19 targets, fully enumerated** — nine at 4⁹ and ten at
4¹⁰, **1.28 × 10⁷ exactly labelled structures** — so the certified global optimum of every objective
is known exactly at every rung. The statistical unit is the TARGET; the three seeds are averaged
*within* a target before any interval is formed.

The question was: *at what objective quality does CVaR-VQE begin to beat its classical controls, and
where does the real folding objective sit relative to that boundary?* The answer is that **the
premise is wrong**.

**The control ladder is the part §5.5 should have led with.** §5.5's control was uniform random
sampling, which this workstream labels **"WEAK CONTROL, proves nothing"** — and it is right to. On
the same instrument, against real classical optimisers at matched evaluations:

| control | cost | result |
|---|---|---|
| **greedy 1-opt with restarts** | matched (8,192 evals) | **the VQE is significantly WORSE at 10 of 10 objective-quality rungs** |
| simulated annealing | matched | the VQE is **never** significantly better; significantly worse at 6 of 10 rungs |
| **simulated annealing at ¼ the budget** | 0.25× | **statistically indistinguishable at 9 of 10 rungs**, every effect ≤ the 0.08 Å floor |
| best-of-N from the untrained circuit | matched | the VQE wins only above ρ\* = +0.355 |
| uniform random | matched | *(weak control; the VQE wins above ρ\* = +0.355)* |

**Greedy 1-opt reaches the certified global optimum in 68–100% of cells; the CVaR-VQE reaches it in
0–32%.** On the structured ladder the VQE's returned structure sits 90% of the way from the
untrained circuit's best-of-N to the objective's certified argmin (`diff = +0.051 + 0.904·Δ`,
R² = 0.859). **The VQE is a partial, expensive implementation of "find this objective's optimum",
and the classical arms implement it completely.**

**There is no ρ boundary, and two independent quality knobs prove it.** The `blend` and `noise`
families give *different* ρ\* for the same control (+0.355 against +0.140), a pooled fit of the
paired difference on ρ has **R² = 0.014**, and the implied crossings (−2.74 and −1.10) are outside
the possible range of a correlation. The disagreement between the two parameterisations *is* the
finding, and the brief pre-registered that reading.

**What does govern it, measured.** The one coordinate both families agree on is how much better the
objective's own certified argmin is than what the control already returns:

| candidate axis | Spearman | R² pooled | R² blend |
|---|---|---|---|
| ρ global | −0.059 | 0.014 | 0.007 |
| ρ in-tail (best 1%) | −0.301 | 0.066 | 0.020 |
| argmin percentile | +0.408 | 0.038 | 0.070 |
| **`RMSD(argmin E) − RMSD(control)`** | **+0.636** | **0.427** | **0.859** |

And decisively: **the VQE-versus-classical gap does not depend on objective quality at all.** Across
the whole ten-rung ladder (ρ = +0.264 → +1.000), `vqe − greedy` is **+0.098 ± 0.037**, `vqe − anneal`
**+0.064 ± 0.035**, `vqe − anneal at ¼ budget` **−0.043 ± 0.041**. Flat lines that never cross zero.
**Improving the objective moves both arms together; it does not move the comparison.**

**Where the real objectives sit** (native-free objectives, not oracle diagnostics):

| objective | ρ global | argmin percentile | **RMSD of its own argmin** | vqe − greedy |
|---|---|---|---|---|
| `S14_disto_bayes` | **+0.485** | 0.282 | 3.077 | −0.010 |
| **`hamil`** (deployed) | +0.264 | **0.231** | **2.661** | +0.121 |
| `legacy` | +0.105 | 0.538 | 4.085 | +0.056 |
| `prior` | **+0.005** | 0.389 | 3.463 | +0.131 |
| *space reference* | | | best **1.030**, mean 4.003 | |

**The two axes disagree, and that is the most important row.** `S14_disto_bayes` has nearly twice
`hamil`'s global ρ and a *worse* optimum (3.077 against 2.661 Å). `prior` has ρ = +0.005 — no
ordering skill at all — and yields the *largest* VQE gain over the untrained control. Across all 95
(objective × target) cells the paired difference correlates **−0.038 with ρ**, **+0.519 with the
argmin percentile**, and **+0.867 with `RMSD(argmin) − RMSD(control)`**.

**So the problem does not sit on one side of a phase boundary, because there is no phase boundary.**
It sits at `RMSD(argmin of the best real objective) = 2.661 Å` against a space best of 1.030 Å — and
**a classical 1-opt with restarts reaches that 2.661 Å exactly, at the same budget.**

**A directional tension between workstreams, stated rather than hidden.** On the *real* objective,
QPHASE's argmin readout gives `vqe − untrained = −0.309` (the VQE **better**, CI [−0.649, +0.037],
null) and its ensemble readout −0.241 (better, null), where §5.5's 75-member coordinate-average
readout gives +0.112 to +0.204 (the VQE **worse**). The readouts and controls differ — best-of-N
argmin, a 5-member average, a 75-member average — which is a legitimate explanation and is
consistent with §5.5's own finding that the terminal operator consumes diversity. But **the sign of
the untrained-circuit contrast is readout-dependent**, and the sprint's conclusion should be stated
as: *against the strongest classical controls the VQE loses at every rung and on every readout;
against the weak untrained-circuit control the sign depends on how the ensemble is consumed.*

**And the programme's last standing positive VQE result is retired.** Sprint 15 recorded that *the
one place a VQE wins is an ensemble consumed without a ranker* (+0.36 to 0.57 Å, itself flagged
low-power and awaiting replication). On the full enumerated register, with annealing and greedy in
the table, **both classical searches beat the VQE on the ensemble readout as well — greedy at every
rung, +0.27 to +0.62 Å.** The exception does not survive, and it was the only surviving positive
quantum result in the programme.

**What CVaR's α actually is: a temperature, and a classically reproducible one.** A
**diversity-matched classical thermostat reproduces α's entire ensemble effect** (Pearson +0.93 and
+0.98) and moves the readout **1.6–1.7× further** than α does. The recorded 150× diversity range
attributed to α is a property of the **optimisation budget**, not of α — at this budget it is 1.35×.
This supersedes the programme's "α is a diversity dial" reading: the *effect* is on diversity, but
the *mechanism* is thermal and a classical thermostat does it better.

**And the weak-optimiser escape is closed.** The standing hypothesis was that CVaR's known
gradient-baseline defect made it a weak optimiser and that fixing it would rescue the arm. It is
**refuted three ways**: the defect is null at the target unit and significantly *worse* on the
ensemble readout; its recorded "de-facto step-size reduction" mechanism is **impossible under Adam**
(a 2.6× smaller gradient gives a 7% smaller step); and cutting the learning rate 16× costs +0.156 Å
on the argmin readout and +0.676 Å on the ensemble — **more optimisation is better here, not worse.**

**What the VQE does do, measured on four axes.** It **samples** (the expected energy under its
distribution moves from the 50th percentile to the 15th–25th) and it **represents** (3.6× near-native
mass enrichment at the real objective; the distribution's mode moves 3.792 → 3.014 Å). It does **not
optimise to completion** (0.5% of its mass sits on the certified optimum) and it explicitly does
**not explore** (3,551–4,020 distinct configurations against annealing's 7,265 at the same budget).
That is a coherent and rather specific description of a device, and it is not the description of the
thing this pipeline needs.

**One number here is directly on the "what next" axis.** The deployed objective's certified argmin
is **2.661 Å on the enumerated instrument, better than the incumbent's 3.204 Å**, and a classical
1-opt reaches it exactly in 100% of cells. The gap between 3.204 and 2.661 is not a search problem.


### 5.9 AMBER as a repair *operator* — and the removal of the sprint's last accuracy claim

`s16/repair.py`, `s16/repair_FINDINGS.md`. Every previous section consumes AMBER as a **scalar** —
a number that ranks structures. The mandated architecture ends with "AMBER all-atom
repair/refinement", which is an operator that **moves atoms**, and that operator had never been
swept, never been frame-nulled at more than one setting, and **its displacement had never been
looked at**. n = 126, target as the unit, convergence gate declared before use and reported with
its exclusions every time.

**The pipeline reproduces bit-identically** against `s15/results/phys_repl_canon0.json` in a fresh
interpreter on a different day — max |Δ| = **0.000e+00 Å** on both arms, final energies agreeing to
0.000e+00 kcal/mol — so these numbers sit on the same instrument as the Sprint 15 record and can be
compared to it directly.

#### The reference points, which is where the result is

| arm | RMSD | vs *doing nothing* |
|---|---|---|
| **do nothing** — the raw all-atom coordinate average | **3.0498 Å** | — |
| ideal-geometry projection (the incumbent's repair stage) | 3.2052 Å | **+0.1554 [+0.1217, +0.1911]**, 27W/99L |
| **AMBER k = 30** (gated, n = 123) | 3.1831 Å | **+0.1333 [+0.1031, +0.1645]**, 24W/99L |
| unrestrained AMBER (k = 0) | 3.4775 Å | +0.4309 |
| ZERO-INFORMATION constant ideal α-helix | 4.0696 Å | +1.020 |

> **The programme's −0.022 Å is the gap between two ways of being 0.13–0.16 Å worse than the
> operator that does nothing at all**, and each loses to it on 99 of 123 targets.

#### And the direction is worth nothing over a random one

A matched-magnitude isotropic random displacement of the same input — the control AUDIT made
mandatory — reaches **3.1606 Å against AMBER's 3.1831** (AMBER − random **+0.0225 [−0.0138,
+0.0571]**, 43W/83L) and 3.1826 against the projection's 3.2052 (+0.0226 [−0.0194, +0.0640]).
**Neither operator's direction is worth anything over a random direction of its own size.**

Worse than that, and this is the answer to the question the workstream was set:

> **ORACLE.** cos(v_AMBER, r) = **−0.0521 [−0.0924, −0.0129]**, median −0.073, positive on only
> **37.3%** of targets, against a matched random control at −0.003 / +0.017 / +0.018. Paired,
> **AMBER − random = −0.0491 [−0.0996, −0.0009]** — *a confidence interval excluding zero in the
> wrong direction.* The projection is −0.0371 [−0.0814, +0.0054].

**The physics moves the structure slightly away from the truth**, and significantly so relative to a
random direction.

**And it moves where RMSD cannot see.** The true residual at this input structure is **60.2% in the
loud half** of the Jacobian spectrum, 21.4% quiet, 18.4% non-torsional (isotropic null
0.302/0.333/0.365) — so the error *is* mostly in the directions that move the score. AMBER's
displacement is the opposite: **55.1% non-torsional (1.5× enriched) and only 23.5% loud (depleted)**,
and in torsion space it occupies loud 0.242 against a null of 0.381, quiet 0.488 against 0.421, and
the **exactly-null directions 0.270 against 0.199**. It is adjusting bonds and angles, which is what
a force field is for — and those are precisely the coordinates a Cα-RMSD is blind to. In torsion
space its cosine with the true error is +0.007 [−0.037, +0.051]: nothing.

**The algebra closes exactly**, which is how we know the effect is a magnitude and not a direction:
predicted against realised RMSD agrees to **0.0044 Å**, and the orthogonal-move cost is **+0.1177 Å
for AMBER against +0.1515 Å for the projection**. **AMBER wins by moving 0.09 Å less, not by moving
better.**

This converges with §5.3 from the opposite side. There, no native-free *geometric* channel had a
usable direction from this same starting structure. Here, the *physics* channel does not either —
and it is the only channel in the programme that was never tested for one.

#### No setting rescues it

Five (k, depth) settings, all reported, on the frame-reproducible subset:

| setting | vs projection, frame-reproducible subset | n |
|---|---|---|
| k = 30, full depth | −0.0095 [−0.0242, +0.0056], 39W/48L | 87 |
| k = 30, depth 400 | −0.0061 [−0.0223, +0.0083], 30W/41L | 71 |
| k = 30, depth 100 | +0.0117 [−0.0165, +0.0367], 9W/22L | 31 |
| k = 30, depth 25 | +0.0109 [−0.0136, +0.0353], 19W/28L | 47 |
| k = 0, full depth | **+0.1602 [+0.0502, +0.2724]**, 21W/49L | 70 |

**There is none.** Every setting either contains zero with a losing or near-even record, or is
significantly worse.

**And the frame null fails at every converged setting.** The rotated-lab-frame difference must be
zero by rigid invariance. Gated, its *mean* at k = 30 is small (−0.00095 [−0.0060, +0.0042],
reproducing Sprint 15's ≈0.004 Å floor) — but its **maximum is 0.1373 Å, six times the effect being
claimed**, and 122 of 122 converged targets move by more than 1e−6 Å. Against the standing band
(|mean| ≤ 0.005 **and** max ≤ 0.05) **every full-depth arm fails**. The only arm that passes is
depth 25, and the convergence gate throws away **76 of 126 targets** there, so it is not a usable
operator. **A frame null must be reported with its maximum, not only its mean**: a mean-only reading
makes this pipeline look like it has a 0.004 Å floor when its tail is 0.14 Å at k = 30 and 2.65 Å
unrestrained.

#### The validity claim survives, but only as a conjunction

AMBER k = 30 against the projection: Ramachandran-favoured **0.466 → 0.874**
(−0.4077 [−0.4529, −0.3628], **116W/2L**), clashes **1.397 → 0.000** (**63W/0L**). Sprint 15's
numbers reproduce exactly.

**But a constant ideal α-helix scores rama 1.000 and 0.000 clashes by construction**, and beats
AMBER on Ramachandran **0W/65L**. A validity statistic that a zero-information reference *maximises*
is not evidence about a force field. The defensible claim is the conjunction, and it is a real one:

> restrained relaxation reaches Ramachandran-favoured 0.874 and zero clashes **while moving only
> 0.725 Å — 26% of the residual — from its input**, where the α-helix reaches the same validity by
> moving 0.883 of the residual, i.e. by discarding the structure.

The k-shape was checked on a **pre-registered n = 30 subsample** at k ∈ {5, 15, 60}, fixed before any
RMSD was seen and used to select nothing: it is **monotone toward "do not move"** (RMSD 3.700 → 3.604
as the displacement falls 0.971 → 0.735), and its limit is the do-nothing structure at 3.506.

**And the restraint constant is an accuracy knob, not a validity knob.** Unrestrained AMBER reaches
essentially the same validity as k = 30 (rama 0.868 vs 0.874, clashes 0.000 vs 0.000, geometry
deviation 0.0170 vs 0.0192) while costing **+0.2742 Å [+0.1785, +0.3823]** against the projection.
The force field supplies the stereochemistry; **the restraint supplies nothing but the decision not
to move.** The k-ladder is monotone toward "do not move", and its limit is the do-nothing structure.

#### What this does to the rest of the report

**It removes the sprint's last surviving accuracy claim.** §5.7's AMBER effect,
−0.0234 [−0.0375, −0.0099] at n = 126, is measured against the **projection**, and the projection is
itself 0.155 Å worse than not repairing at all and no better than a random displacement of its own
size. The effect is real as an arithmetic difference and it is **not an accuracy result**: it is the
gap between two ways of degrading the structure. Mechanically, the orthogonal-move cost is +0.1177 Å
for AMBER against +0.1515 Å for the projection, a gap of 0.0338 Å against a realised 0.0221–0.0234 Å
— **the −0.022 Å is a displacement-magnitude difference, not a physics result.**

**What survives, and it is worth keeping.** The mandated pillar is satisfiable and should be kept —
but as a **terminal validity operator with a declared accuracy price**, never as a refinement step.
On this instrument that price is **+0.133 Å against not running it**, bought for **12.14 s per
target** (≈2.7× the projection it replaces, ≈10³× an AMBER single point). If the architecture wants
buildable, clash-free, Ramachandran-favoured output, restrained AMBER at k = 30 is the cheapest
thing measured here that delivers it *while staying near the input*, and the ideal-geometry
projection it would replace is **strictly worse on every axis except wall clock** — rama 0.466,
1.40 clashes, and a larger accuracy price of +0.155 Å.

**AMBER's repair role is real and its accuracy role is not, and the two must never again be quoted
as one number.**


---

## 6. Mechanism: what the sprint now believes, and how much of it is arithmetic

Three results in this sprint are **exact** rather than empirical, and it is worth separating them
from the measurements, because exactness is what makes the measurements interpretable and it is also
what makes them easy to over-read.

**Exact (theorems; measuring them validates an implementation and nothing more):**

- `‖avg − native‖² = mean member error² − diversity²` — the Krogh–Vedelsby ambiguity decomposition,
  verified through the terminal operator at machine precision. This is the parallel-axis theorem.
- `‖r − tv‖` is minimised at `t* = c‖r‖` giving `‖r‖√(1 − c²)` — so coordinate-space RMSD is exactly
  linear in motion toward a superposed target. Reproducible on random point clouds with no protein
  in them (4.6e−07 Å).
- Both survive Kabsch superposition, and provably so: a sum of two symmetric PSD cross-covariances
  is PSD, so the optimal rotation does not change along the segment.

**Measured (and this is where the science is):**

1. **The torsion-space fit is very nearly noise.** 90° per-angle RMS, 87% of a uniformly random
   torsion assignment. Everything else follows from this. Because the displacement to the native is
   that large, the torsion→Cartesian map is far outside its linear regime (curvature ratio 2.709),
   the Jacobian's loud/quiet split describes a neighbourhood the error does not live in, and partial
   motion toward the truth does not pay partially. **No steering scheme in that space can work, and
   this is not a property of the space — a matched zero-information control behaves the same or
   worse.**

2. **The accuracy the pipeline has comes from averaging, and the repair stages spend it.** The raw
   all-atom coordinate average is **3.0498 Å**; the ideal-geometry projection emits 3.2052 Å
   (+0.155, **27W/99L**) and restrained AMBER 3.1831 Å (+0.133, **24W/99L**). Both repair operators
   are worse than doing nothing on four targets in five, **neither beats a matched-magnitude random
   displacement**, and **AMBER's displacement has a negative cosine with the true residual**
   (−0.052 [−0.092, −0.013]; against random, −0.049 [−0.100, −0.001]). The programme's accuracy is
   produced entirely by retrieval and averaging, and 0.13–0.16 Å of it is spent making the result a
   legal peptide.

3. **Averaging works because of the diversity term, and optimisation destroys it.** Generators with
   near-equal conformer quality produce readouts spanning 0.30 Å, and 42–58% of that is the
   diversity term. This is why concentrating the search hurts: **concentration is the destruction of
   the diversity term.** It explains "searching harder on a bad objective loses" without needing the
   objective to be bad. **On α the sprint goes further than the recorded "diversity dial" reading**:
   a diversity-matched **classical thermostat reproduces α's whole ensemble effect** (Pearson +0.93,
   +0.98) and moves the readout 1.6–1.7× further, so α's effect is thermal and is done better
   classically (§5.8).

4. **But the diversity term cannot be bought.** Selecting for it raises member error at almost
   exactly the rate it raises diversity (+22% against +4.5%, net ≤ 0.008 Å), because the two enter
   the identity with equal and opposite weights. The mechanism is exact and **has no lever**.

5. **And the physics moves in the coordinates the score cannot see.** The true residual is **60.2%
   loud**; AMBER's displacement is **55.1% non-torsional and only 23.5% loud**, over-occupying even
   the exactly-null torsion directions (0.270 against a 0.199 null). A force field adjusts bonds and
   angles; a Cα-RMSD is blind to exactly those. The −0.022 Å it appeared to buy is a
   **displacement-magnitude difference** — orthogonal-move cost +0.118 Å against the projection's
   +0.152 Å — and the algebra closes to 0.0044 Å.

6. **What is left of the residual is a scalar, and its sign is not native-free.** From the best
   start, the residual to the native is dominated by a compactness scalar whose magnitude any
   extended reference recovers (a constant β-strand reaches |c| = 0.398 against a 0.128 null) and
   whose **sign flips per target**. Four native-free estimators predict that sign **worse than a
   constant guess**. The ORACLE ceiling behind that one bit is **2.491 Å** — below the sprint's own
   primary goal — which is exactly why it must not be quoted as anything but a ceiling.

7. **The failure of the energies is a property of the potential class, not of short peptides.** On
   19 exhaustively enumerated spaces, where every argmin is certified rather than searched, neither
   Legacy nor AMBER can be distinguished from a random feasible member — but a **distance-restraint
   objective's argmin is −0.926 Å [−1.440, −0.431]** better than random, on 15 of 19 targets, with
   the native at the 29.1st percentile. The problem is not that 9–16-mers are unrankable; it is that
   **contact-style potentials do not rank them.**

8. **And the quantum question has no phase boundary to sit on either side of.** Two independent
   objective-quality knobs give different crossing points for the same control, a fit of the paired
   difference on ρ has R² = 0.014, and the VQE-versus-classical gap is **flat across the whole
   quality ladder** (`vqe − greedy = +0.098 ± 0.037`, never crossing zero). What governs the
   comparison is not objective quality but **how much better the objective's own certified argmin is
   than what the control already returns** (R² 0.427 pooled, 0.859 structured). The VQE lands 90% of
   the way to that argmin; **greedy 1-opt lands on it exactly**, in 68–100% of cells against the
   VQE's 0–32%. The quantum arm is a partial, expensive implementation of something a classical
   1-opt implements completely.

**The one architectural conclusion that survived every control**: ensembles must be read out in
**coordinate** space, not torsion space (+0.446 [+0.193, +0.706] at m = 75, 8/9 targets, growing
with ensemble size). It follows from (3): torsion averaging cannot cash the diversity term.

---

## 7. Methodology: the failure mode this sprint found in itself

This section exists because the sprint's methodological result is more transferable than its
scientific one.

**Two independent audits, on two sprints, reached the same verdict: every number reproduced exactly,
and a large fraction of the conclusions built on those numbers were still wrong.** AUDIT reproduced
eleven load-bearing Sprint 15 figures to the digit. VERIFY reproduced every published Sprint 16
number from its artefact. In both cases the errors were elsewhere.

**The failure mode has a name now:** *a quantity is measured correctly and then read as if it were a
different quantity.* Eight instances, in this sprint's own work and in the record it inherited:

| what was measured | what it was read as | cost |
|---|---|---|
| the arithmetic mean of a **signed** cosine | the input to a law that is **quadratic and two-sided** in it | turned a missing sign bit into a phantom absence of signal; the sprint's most consequential error |
| the **free-superposition** member error | the **common-frame** term the identity consumes | the diversity attribution overstated ≈2× |
| the **arithmetic** mean of two channels' RMSD | the **quadratic** mean the identity requires | 54% of a "prediction error" that was never a prediction |
| a **row** bootstrap over 162 cells | a target-level interval over 9 targets | the sprint's only positive accuracy effect |
| an **algebraic identity** | an empirical validation | a 31-year-old theorem presented as a discovered law, twice |
| a frame null's **mean** | the null's behaviour | a 0.004 Å floor reported where the tail is 0.137 Å — six times the effect being claimed |
| a validity statistic **in isolation** | evidence about a force field | every such statistic here is *maximised* by a constant α-helix that knows only the chain length |
| an operator's **effect against a bad baseline** | the operator's value | −0.022 Å against a projection that is itself 0.155 Å worse than doing nothing |

**Two rules adopted, in order of power.**

> **AUDIT's rule.** Every arm must carry **both** a zero-information reference control (a constant
> ideal α-helix, a random pool window) and a matched random-direction control. Without these two,
> this programme's evidence read about **five times stronger** than it was.

> **VERIFY's rule, which is stronger, because a zero-information control alone would not have caught
> three of the five above.** *Before quoting a statistic as evidence for a law, derive which
> functional of it the law consumes, and quote that functional.* Where a law is quadratic,
> sign-invariant, or frame-dependent, the arithmetic mean of the raw per-target quantity is not it.
> **A mechanism sentence needs its own control. An exhibit derivable on random point clouds is not
> an exhibit.**

**And the rules were broken in this document, after being written into it.** The first draft of §5.5
quoted **row-level** bootstrap intervals and row-level win/loss counts for the quantum table, and
asserted they had survived the target-level recomputation — the fourth of the five failure modes
listed immediately above, committed two pages after it was catalogued, in the section that catalogues
it. One of the four intervals (`cvar1.0`, +0.209) touches zero at the target level and had been
printed as excluding it. A fourth workstream, reviewing the report as a document rather than the
science, caught it. **This is recorded rather than quietly fixed, because a rule that its own author
breaks inside the same document is evidence about how weak the rule is on its own** — the defence
that worked was not the rule but the separate reader.

**A third observation, which is not a rule but should be.** Sprint 15's verification pass **passed**
the fusion-law row, because every figure traced to its artefact exactly. **A transcription audit
cannot see a formula error**, and three hostile reviewer passes missed it for the same reason. An
audit must re-derive, not re-check.

**On small n.** The programme has now had six conclusions reversed by n ≤ 8 reads, and this sprint
supplied the first where the small-n interval **excluded zero** (−0.106 [−0.185, −0.031] at n = 8
against −0.003 [−0.011, +0.007] at n = 126). **A confidence interval at small n is not protection**,
because the failure is in the sample, not the interval.

---

## 8. Limitations

- **Two instruments, different powers.** The steering, coordinate-steering and diversity-selection
  results are at **n = 126**. The integrated architecture, the Legacy/AMBER comparison and the
  diversity decomposition are on the **9 enumerated targets** — the price of an exhaustively
  enumerable space with a certified optimum. Nine targets is inside the range where this programme
  has repeatedly had conclusions reverse, and the target-level recomputation of the AMBER result is
  exactly that risk materialising.
- **The step size is trained, not native-free.** Every steering arm's reported effect is *a
  direction plus a step fitted on labels leave-fold-out*. That is the same convention the programme
  uses for its separation debias and it is fold-honest, but it is not a native-free quantity and is
  never described as one.
- **The ladder shaped one null.** `csteer.STEPS` had no negative rung; on 45–67% of targets no arm
  could have helped regardless of the direction. That null is partly an artefact of the design and
  is reported as such.
- **The compactness bit was closed by re-reading an existing artefact, not by a new experiment.**
  Sprint 15's four native-free scale references are a family — distance-matching estimators. Whether
  some estimator *outside* that family could supply the sign is **untested**, and is stated as an
  open question rather than a closed one.
- **`s16/req.py` was written and not run.** The flagship answered its question by a shorter route.
  The module is kept unrun so the reasoning is inspectable.
- **Everything in §5.9 is measured at ONE input structure** — the top-75 all-atom coordinate
  average. Whether the repair operators behave the same way on a different input (a VQE ensemble
  member, a better-selected average) is untested, and the conclusion "worse than doing nothing"
  is a statement about repairing *this* structure.
- **The frame null fails at every converged setting**, so no per-target AMBER accuracy number
  anywhere in this report is resolvable. Only aggregates are quoted, and they are quoted as
  arithmetic differences rather than as physics.
- **Two seeding-rule violations were found in the Sprint 15 tree, one fixed and one not.**
  `s15/align_jac.py:136` seeded the random-subspace null with `hash()` (salted per process, then
  collapsed onto 7 seeds); REPAIR fixed it to `stable_rng` and verified that the **aggregate does
  not move** (0.4961 → 0.5005 against an analytic 0.5000) while the **per-target value moved by up
  to 0.128** depending on which interpreter drew the null — the aggregate was safe, the per-target
  values were not. `s15/info_regime.py:59` carries the same violation and is **not yet fixed**.
- **Power was not stated anywhere until review, and on the 9-target instrument it is low.** Standard
  error 0.051 Å, minimum detectable effect at 80% power ≈ **0.16 Å**, power to detect 0.05 Å ≈ **9%**.
  Every null quoted on that instrument below about 0.16 Å is uninformative rather than negative, and
  should be read that way — including the Legacy nulls in §5.4.

---

## 9. The next question

**The sprint closes three directions, leaves one open, and relocates the question.**

Closed: torsion-space steering (§5.1–5.2), coordinate-space steering by distance-matching estimators
(§5.3), the search hypothesis for the quantum arm (§5.8), and **AMBER relaxation as an accuracy
operator** (§5.9) — the last of which also closes the only channel in the programme that had never
been tested for a native-free direction: the physics supplies one, and it points the wrong way. Open: one bit and one
magnitude (below). Relocated: the gap that matters is no longer between the pipeline and a better
optimiser, it is between the pipeline and its own objective's certified optimum.

**Closed.** Steering a torsion-space fit by any native-free direction estimate — because the fit is
87% noise and the map is nowhere near linear at that displacement. And steering the coordinate
average by any distance-matching scale estimator — because the residual is a compactness scalar
whose sign those estimators predict worse than a constant guess.

**Open, and precisely stated.** From the pool coordinate average at 3.048 Å, a **zero-information**
constant β-strand direction has an ORACLE two-sided ceiling of **2.491 Å**. That ceiling is computed
with **two** native quantities given: the per-target sign of `c`, and the step magnitude `t* = c‖r‖`,
in which `‖r‖` is the distance to the native. The open question is therefore two questions, and only
the first is a one-bit classification: **(i)** can any native-free signal *outside* the
distance-matching family — sequence composition, predicted secondary structure, the retrieval pool's
own dispersion, a classifier trained fold-honestly on the sign alone — beat 0.516 at predicting the
sign; and **(ii)** what does the 2.491 Å ceiling degrade to when the magnitude is *estimated* rather
than given? **The second is untested and the ceiling is not meaningful without it** — a
leave-fold-out constant step is the cheapest first answer and was not run.

**Why it is worth asking and why it may not be.** It is worth asking because it is the only line in
the programme whose ceiling sits below the primary 2.5 Å goal, and because it is a *classification*
problem on one bit rather than a regression problem on a structure — though the easier-looking
half, the step magnitude, is a regression and is untested. It may not be worth asking
because the recorded native-free compactness proxies reach only 0.24–0.37 against an oracle's 0.909,
and because the same estimators that predict the *distance* scale correctly (0.54–0.61) are
anti-informative about the *RMSD-optimal* scale — which suggests the two quantities are genuinely
different and that distance-based information may be structurally incapable of supplying this bit.

**A second open question, and it is cheaper than the first.** The raw coordinate average is
**3.0498 Å** and everything the pipeline does after it makes it worse — the projection by 0.155 Å,
restrained AMBER by 0.133 Å, both losing on 99 of 123 targets, neither beating a random displacement
of its own size. That is a **0.13–0.16 Å validity tax**, and it is the largest single number in the
programme that is not a prediction error. The question is whether buildable, clash-free geometry can
be bought for less — a smaller displacement, a repair restricted to the non-torsional coordinates it
already occupies, or a different operator entirely. Nothing in this sprint tried; the measurement
that says it is worth trying is new.

**A third thing the sprint learned about where to look, which belongs here rather than in a
footnote.** On the enumerated instrument the **deployed objective's certified argmin is 2.661 Å**,
against the incumbent's 3.204 Å and a space best of 1.030 Å — and a classical 1-opt with restarts
reaches that 2.661 Å **exactly, in 100% of cells, at the same budget.** So the 3.204 → 2.661 gap is
**not a search problem**, and no better optimiser, quantum or classical, will close it. It is a gap
between the objective's optimum and what the pipeline actually returns, which makes it a question
about the terminal operator and the selection stage — the two places this sprint's exact
decomposition applies. And the 2.661 → 1.030 gap is an objective problem, which is where the
distance-restraint result above says to look.

**What must not happen next.** The sealed 60-target benchmark must stay sealed. None of the five
unlock conditions is met, and the first fails outright: **the architecture is not frozen, because
nothing in this sprint earned a place in it.**

**§10 turns the questions above into a staged plan** with measured ceilings, pre-registered
falsifiers, and an explicit statement of which stage the evidence supports and which it does not.


## 10. The path to sub-2.5 Å, and to sub-2.0 Å, keeping all four pillars

This section is a plan, not a result. Every number in it is measured; every inference from those
numbers to a proposed experiment is labelled as an inference. It is written to be attacked.

### 10.1 The gap ladder, measured

All on the 126-target instrument, ORACLE quantities except where marked.

| level | Å | what it is |
|---|---|---|
| shipped argmin readout | 3.454 | what the pipeline returns through its own selector |
| **incumbent (deployed)** | **3.204** | coordinate average → projection → AMBER |
| projection removed | 3.205 | the projection stage alone |
| **raw all-atom coordinate average — "do nothing"** | **3.050** | the best structure the pipeline actually builds |
| deployed objective's **certified** argmin *(19 enumerated targets, coarse k = 4 lattice — does not transfer directly)* | 2.661 | reachable **exactly** by classical 1-opt in 100% of cells |
| **ORACLE best inside the shipped top-75** | **2.306** | the ceiling of any reranker that consumes the current filter |
| **ORACLE best inside the K = 500 retrieval pool** | **1.711** | the ceiling of any reranker that consumes the whole pool |
| perfect distance knowledge (recorded) | ≈1.95–2.0 | what an oracle distogram buys through the current geometry |
| enumerated space best *(19 targets)* | 1.030 | the lattice contains much better than any objective finds |

**Four gaps, and only two of them are search problems.**

| gap | size | character |
|---|---|---|
| 3.204 → 3.050 | **0.155 Å** | **a validity tax.** Both repair operators are worse than doing nothing (§5.9) and neither beats a random displacement of its own size. Recoverable if a cheaper validity operator exists. |
| 3.050 → 2.306 | **0.744 Å** | **a readout/selection gap.** The pipeline averages where the best member is far better. Not a search problem — the objective's own certified argmin is reachable classically (§5.8). |
| **2.306 → 1.711** | **0.595 Å** | **the filter throws away the answer.** The distogram top-75 cut costs 0.595 Å of *ceiling*. On the 18 zero-recall targets it costs 2.393 Å (4.677 against 2.284). |
| 1.711 → 1.030 | 0.681 Å | a retrieval/library gap; out of scope for a selection architecture. |

**The arithmetic that matters: perfect selection inside the existing K = 500 pool gives 1.711 Å.**
Sub-2.0 Å is reachable *in principle without any new structure generation at all* — the answers are
already in the pool on most targets. 53 of 126 targets have a pool best above 2.0 Å, so the mean is
carried by the rest; that is a real limitation and it is why 1.711 is a mean, not a promise.

### 10.2 The one structural change the evidence demands

**The architecture must move from an averaging readout to a selection readout, and the two cannot be
mixed.** This follows from the sprint's own exact identity rather than from preference:

- An **averaging** readout obeys `readout² = mean member error² − diversity²`. Both terms move
  together under any selection pressure — that is precisely what `s16/divselect.py` measured and
  why it was a null (+0.015 [−0.019, +0.049]). **The averaging readout is capped near 3.05 Å and no
  amount of better selection into it will help**, because better members raise the first term as
  fast as they lower the second.
- A **selection** readout has a ceiling of 1.711 Å and is limited by ranking skill alone.

**This is the single largest decision on the table, and the sprint's exact decomposition is what
makes it decidable rather than a matter of taste.**

### 10.3 The blocker, named precisely

The programme's own record says: *nothing ranks within the pool.* In-band ordering is learnable
**within** a target (ρ = 0.986, zero overfitting) but only **0.600 across** targets, and reaching
2.0 Å requires **0.638**. The recorded diagnosis is that capacity is saturated by a *linear* model
and that *the only leverage supplies the per-target sign at inference*.

**Sprint 16 found the same shape in a different place.** The compactness residual has real magnitude
(|c| ≈ 0.3–0.4, well above an isotropic null) and an unpredictable per-target **sign**; four
native-free estimators predict that sign worse than a constant guess. **Two independent lines now
say the missing quantity is per-target calibration, not per-target prediction.** That convergence is
the most actionable thing the sprint produced, and it is what the plan below attacks.

### 10.4 Stage 1 — recover the validity tax (target: 3.204 → ≈3.05 Å)

**The cheapest 0.155 Å in the programme, and it is not a modelling problem.**

`s16/repair_FINDINGS.md` shows the repair stage costs 0.155 Å (projection) or 0.133 Å (AMBER),
loses to doing nothing on 99 of 123 targets, and that AMBER's displacement is **55.1%
non-torsional** — it moves in bond and angle coordinates that Cα-RMSD is blind to, and only 23.5%
in the loud directions where the error actually lives.

**Experiment 1a.** A repair operator restricted to the coordinates AMBER already prefers: minimise
under the force field with the Cα trace **held fixed or tightly restrained**, letting only
non-Cα degrees of freedom relax. The measurement above predicts this keeps essentially all of the
validity (which lives in the non-torsional coordinates) at a fraction of the accuracy price.
**Pre-registered success: validity within 0.02 of AMBER k = 30 on Ramachandran and clashes, with an
accuracy price below 0.05 Å against doing nothing, n = 126, frame null reported with its maximum.**
**Falsifier: if the accuracy price stays above 0.10 Å, the tax is intrinsic and Stage 1 is closed.**

**This is where AMBER stays in the architecture, and it is an honest place for it** — a terminal
validity operator with a declared and now much smaller price, which is exactly what §5.9 concluded
it is good for.

### 10.5 Stage 2 — switch the readout, and make Legacy and AMBER earn their place in it (target: ≈2.5 Å)

**Experiment 2a — the filter is the ceiling; widen it.** The top-75 cut costs 0.595 Å of ceiling and
2.393 Å on the 18 zero-recall targets. Re-run the pipeline with candidate sets of 75, 150, 300 and
the full 500, measuring the **ORACLE ceiling** and the **realised** result separately at each width.
**The expected result is that the ceiling improves and the realised result does not** — because
ranking is the blocker, not recall. **That is the point:** it separates the two failure modes that
the current architecture confounds, and it is cheap.

**Experiment 2b — the three-way readout comparison, on the same candidate set.** This is where the
mandated Legacy/AMBER comparison belongs, and it is the one mode in which they have *never* been
compared. Every arm consumes an identical candidate set; only the readout differs:

| arm | readout |
|---|---|
| `avg` | coordinate average (the incumbent's readout, the control) |
| `sel-dist` | select the argmin of the distance-restraint objective |
| `sel-legacy` | select the argmin of the genuine 11-component Legacy total |
| `sel-amber` | select the argmin of a genuine ff14SB/GBn2 single point |
| `sel-legacy→amber` | Legacy gates, AMBER ranks the survivors |
| `sel-random` | **the control every arm must beat** |
| `sel-ORACLE` | the ceiling, labelled |

**The prior from this sprint is specific and unflattering to the energies**: on 19 exhaustively
enumerated spaces the distance objective's certified argmin is −0.926 Å [−1.440, −0.431] better than
random while **neither Legacy nor AMBER is distinguishable from random** (§5.7). **So the
pre-registered expectation is that `sel-dist` wins and both energy arms fail.** If that holds, the
architecture's ranking stage is structural and the energies are auxiliary — which is what the
evidence already says. If an energy arm *does* win here, that would be the first time in the
programme, and it would be a genuine discovery about the readout rather than about the potential.

**Legacy's earned role is the gate, not the ranker**: AUROC 0.93–0.99 on its non-tautological
`torsion` component (§5.7). **Experiment 2c**: use Legacy purely as a *rejection* gate before
ranking — drop candidates it flags as unbuildable — and price it against dropping the same number at
random. This is the one Legacy use the sprint's evidence actively supports, and it has never been
run in a selection architecture.

### 10.6 Stage 3 — the ranking problem, which is where sub-2.0 Å lives (target: <2.0 Å)

Stage 2 is expected to land near 2.5 Å, not below 2.0. **Sub-2.0 Å requires closing the in-band
ranking gap: cross-target ordering 0.600 → 0.638.** The sprint's contribution is to say what kind of
quantity is missing.

**Experiment 3a — per-target calibration, the one bit.** Two independent lines (the recorded in-band
result and this sprint's compactness result) say the missing quantity is a **per-target scalar
supplied at inference**, not a better global model. Concretely: learn the in-band ranker as usual,
then learn a *separate, fold-honest* per-target calibration from native-free features only —
sequence composition, predicted secondary-structure content, the retrieval pool's own dispersion and
recall proxies, chain length. **Pre-registered falsifier, and it is strict: the calibrated ranker
must beat the uncalibrated one on held-out folds with an interval excluding zero, and must beat a
zero-information constant calibration.** Sprint 16's most-repeated lesson is that a
zero-information reference reproduces more than anyone expects.

**Experiment 3b — where CVaR-VQE genuinely belongs, and it is not the search.** `s16/qphase` settled
what the quantum arm does and does not do, on 1.28 × 10⁷ certified structures:

- it does **not** optimise to completion (0.5% of its mass on the certified optimum) and does **not**
  explore (3,551–4,020 distinct configurations against annealing's 7,265);
- it **does** sample (expected energy percentile 0.50 → 0.15–0.25) and it **does represent** —
  **3.6× near-native mass enrichment**, with the distribution's mode moving 3.792 → 3.014 Å;
- classical 1-opt beats it at finding an argmin at every rung, so **the search framing is closed.**

That leaves one framing in which the pillar is scientifically live, and this sprint's identity makes
it precise. A candidate-set readout is governed by two terms — member quality and diversity. The
VQE **improves the first** (it moves the mode 0.78 Å toward the native) and **damages the second**
(it concentrates). A classical thermostat does the reverse: it holds diversity but does not enrich.

> **Experiment 3b, stated as a falsifiable question: can a CVaR-VQE sampler be operated so that it
> enriches near-native mass *without* paying the diversity it currently pays — and if so, does the
> resulting candidate set beat a classical thermostat matched on both terms?**

The concrete arms: α annealed low→high across the run; a mixture of VQE samples with uniform samples
at matched budget; and the mandatory controls — greedy 1-opt, simulated annealing at matched and at
quarter budget, and a diversity-matched classical thermostat, which `s16/qphase` showed reproduces
α's whole ensemble effect and moves the readout 1.6–1.7× further. **Success is measured on the
(member error, diversity) plane, not on RMSD**, because the identity says RMSD is a function of
those two and reporting only RMSD is what hid the mechanism for two sprints.

**Pre-registered falsifier, and it should be stated before anyone runs it: if no VQE configuration
reaches the Pareto frontier that the classical thermostat plus 1-opt already occupies on that plane,
the quantum pillar has no accuracy role in this architecture and should be retained only as a
studied object, not as a component.** On the evidence in §5.8 that is the more likely outcome, and
saying so now is what makes the experiment worth running.

### 10.7 What the arithmetic supports, and what it does not

| stage | mechanism | expected | confidence |
|---|---|---|---|
| 1 | restrict repair to non-Cα coordinates | 3.204 → ≈3.05 | **high** — the displacement decomposition is measured, and the operator already prefers those coordinates |
| 2 | widen the filter + switch to a selection readout with the distance objective | ≈3.05 → **2.4–2.6** | **moderate** — the ceiling is 1.711 and the distance objective is the only ranker with a CI excluding zero, but no ranker has yet transferred across targets |
| 3 | per-target calibration of the in-band ranker | 2.4–2.6 → **<2.0** | **low** — it requires closing 0.600 → 0.638 on a problem where five previous attempts have failed, and 53 of 126 targets have a pool best above 2.0 Å regardless |

**The honest summary.** Sub-2.5 Å is a credible target with two structural changes, both of which are
justified by measurements in this report rather than by hope. **Sub-2.0 Å requires solving a ranking
problem the programme has failed at five times**, and the only new information Sprint 16 brings to it
is a sharper description of the missing quantity — a per-target scalar, supplied at inference, which
is the same shape as the compactness bit that closed §5.3. That is a real lead and it is not a plan
that should be assumed to work.

**All four mandated pillars survive this plan, with their roles decided by evidence rather than by
mandate:** genuine VQE and genuine CVaR-VQE as a studied sampler on the (error, diversity) plane
against classical controls that currently beat them; genuine Legacy as a rejection gate, which is
the one role where it has an AUROC of 0.93–0.99; and genuine AMBER as a terminal validity operator
with a declared price, plus a head-to-head ranking comparison against Legacy on identical candidate
sets — the one mode in which the two have never been compared.

**And the benchmark stays sealed** until an architecture is frozen, a control is frozen, a success
criterion is pre-registered, the expected result would genuinely change the conclusion, and no
cheaper internal instrument can answer it. Stage 2 is the first point at which those conditions could
plausibly be met.
