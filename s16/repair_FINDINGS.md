# SPRINT 16 — REPAIR workstream findings

**Subject: AMBER as an OPERATOR THAT CHANGES COORDINATES, which nobody had tested.**

The ENERGY workstream tests AMBER as a *ranking filter* over a fixed retrieval ensemble; the
coordinator tests it as a filter over VQE-generated ensembles. Both consume AMBER as a scalar. The
user's mandated architecture ends with "AMBER all-atom repair/refinement" — an operator that moves
atoms — and that operator had never been swept, never been frame-nulled at any setting but one, and
its **displacement had never been looked at**.

Tiering: **DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED.** Every native-derived quantity
is labelled **ORACLE** in the tables and in the prose. The statistical unit is the **TARGET**
(n = 126). Intervals are `I.paired`'s paired bootstrap with a **fold-clustered twin** printed beside
it; with only five clusters the fold-clustered interval is the *weaker* of the two, not the stronger
one, and both are shown (AUDIT's convention, kept).

**Modules** `s16/repair.py` (runs AMBER), `s16/repair_report.py` (all analysis, no AMBER),
`s16/repair_seedfix.py` (the BRIEF §3.7 violation).
**Artefacts** `s16/results/repair_{A,B,C,D}.json`, `repair_report.json`, `repair_report_*.txt`,
`repair_verify.json`, and the logs beside them.

---

## 0. METHOD, and the three things that make the numbers comparable

### 0.1 The pipeline is PHYS's pipeline, and it reproduces BIT-IDENTICALLY

**DEMONSTRATED.** `s16/repair.py` rebuilds the arm from the same inputs PHYS used
(`s14.avgspace.top75_windows` → `s15.phys_repl.averaged_backbone_from` → `core.amber.refine_coords`
at `steps = 0, tolerance = 1.0, threads = 1, memo = False`; arm A via
`s15.phys_repl.project_both`). Run in a **fresh interpreter on a different day**, against
`s15/results/phys_repl_canon0.json`:

| arm | max \|Δ\| over the checked targets |
|---|---|
| the ideal-geometry projection (arm A) | **0.000e+00 Å** |
| AMBER k = 30 (arm B) | **0.000e+00 Å**, and the final energies agree to **0.000e+00 kcal/mol** |

`python -m s16.repair --mode verify` → `s16/results/repair_verify.json`. Every number below is
therefore on the same instrument as `s15/phys_FINDINGS.md` and can be compared to it directly.

**Consequence for the noise floor.** The AMBER path contains **no RNG** and reproduces bit-for-bit
across interpreters, so the **0.08 Å multi-start false-positive floor does not transfer to it** —
its mechanism (a per-process-salted `hash()` seeding a random multi-start) cannot act here. The
floor that does apply is this path's own: **≈0.004 Å with the convergence gate applied and up to
0.038 Å without it** (`s15/phys_FINDINGS.md` §3.1b). Both are used in their right places below.

### 0.2 The convergence gate

Declared in `core.amber.convergence_flags` before it was applied: a restrained minimisation is
**converged iff its final potential energy, restraint switched off, is finite and ≤ 1000 kcal/mol**.
One threshold, no target-specific tuning, no native-derived quantity. It is **reported, never
silently applied**: every table below prints the number of targets it excludes and their names, and
gives the effect both gated and ungated.

### 0.3 The rotated-lab-frame null, at every setting

ff14SB, GBn2, the positional restraint (anchored at the same, rotated, input) and every RMSD in this
project are rigid-invariant, so relaxing the **same structure in a rotated lab frame is zero by
construction**. Proper rotations only (`det = +1`; a reflection is not a symmetry of ff14SB), drawn
under `stable_rng("s16repair", "frame", 1, pdb)`. The pre-declared PASS band is taken unchanged from
`s16/energy_lib.py` so the two workstreams' nulls are judged by one rule: **|mean ΔRMSD| ≤ 0.005 Å
and max |ΔRMSD| ≤ 0.05 Å on the converged subset.** *An arm whose frame null fails is not
reportable.*

The **frame-reproducible ("stable") stratum** is PHYS's definition, reused unchanged so the strata
are comparable across sprints: converged in both frames **and** |E(frame B) − E(frame A)| ≤ 0.5
kcal/mol. It is native-free.

### 0.4 The two controls that AUDIT showed every arm needs

AUDIT established that this programme's evidence read about five times stronger than it was because
arms lacked a zero-information reference and a matched random-direction control. Both are carried
in every claim here.

| control | what it is | which axis it prices |
|---|---|---|
| **ZERO-INFO `none`** | *do nothing* — the raw all-atom coordinate average, the operators' own input | ACCURACY |
| **ZERO-INFO `helix`** | a constant ideal α-helix (φ = −57°, ψ = −47°) of the same length, superposed on the input. It has perfect Ramachandran statistics and zero clashes **by construction** | VALIDITY (and, as a direction, ACCURACY) |
| **RANDOM matched-magnitude** | the input CA displaced by an isotropic random direction (rigid components removed) scaled to the **same RMS displacement the operator applied**, 3 draws under `stable_rng` | ACCURACY / direction quality |

### 0.5 Cost is quoted at the right number

Every AMBER arm here is a **full restrained minimisation**, not a single point. Measured wall clock
is printed per setting in the validity table. A single point is 8–14 ms; these arms are ~10³ times
that. **No single-point cost is quoted for any minimisation arm anywhere in this file.**

### 0.6 Pre-registration of the pass-D shape check

The sweep was re-scoped mid-run by the coordinator (the original 8-setting × 126-target × 2-frame
design was on a ~4.2 h path and was starving sibling agents). The re-scope is recorded here because
it changes what may be claimed:

- **Pass A** — k = 30, both frames, plus the projection and both zero-information references, at
  full n = 126. This alone settles the gate and the frame null.
- **Pass B** — k = 0 (unrestrained), both frames, full n = 126. With A this is the contrast the
  record turns on.
- **Pass C** — minimisation **depth** at k = 30 (25 / 100 / 400 iterations against "until
  converged"), both frames, full n = 126.
- **Pass D** — the **shape** check in k ∈ {5, 15, 60}, on a **pre-registered, native-free,
  fold × length stratified subsample of 30 targets**, fixed by `s16.repair.shape_subsample` and
  printed *before any pass-D RMSD existed*:

```
1D6X 1I6Y 1JBF 1KYC 1KZ2 1MF6 1NIZ 1S9Z 1TOR 2LER 2LNG 2MAI 2MD2 2MK7 2MQ2
2MSA 2NB7 2NBC 2O0S 3BTB 6OQP 6Q08 6QAX 7JS6 7LCW 8FLP 8IS3 8T63 8UN8 9L1M
```

**Nothing is selected on pass D.** It is a shape check; no k it favours is re-quoted at full n as an
out-of-sample choice. (Note in passing: the stratified draw happens to include three of the four
targets PHYS found non-converged — 1D6X, 1MF6, 2NB7 — which is a reason to read its gated column,
not its ungated one.)

A 5-target partial of the original 8-setting design is preserved as
`s16/results/repair_ladder5_superseded.json` and is **quoted nowhere**.

### 0.7 One BRIEF §3.7 violation fixed, and what the fix actually changes

**DEMONSTRATED.** AUDIT found the sprint's one surviving `hash()` seeding at
`s15/align_jac.py:136` — the random-subspace null that produces the 0.499 figure, drawn from
`np.random.default_rng(abs(hash(p)) % 7 + 5)`. `hash()` is salted per process **and** the modulus
collapses 126 targets onto seven seeds. Fixed to
`SD.stable_rng("align_jac", "subspace_null", p)` and verified by `s16/repair_seedfix.py`, which
measures the quantity under every interpreter salt the old code could ever have drawn, under the new
seeding, and against the analytic expectation:

| | value |
|---|---|
| analytic expectation `k/m` | **0.500000** |
| OLD `hash()` seeding, mean over the 7 reachable salts | **0.4961** |
| NEW `stable_rng` seeding | **0.5005** |
| per-target spread **across interpreter salts** under the old seeding | mean **0.0759**, max **0.1280** |

**The reported aggregate does not change** (0.4961 → 0.5005 against an analytic 0.5000), which is
what AUDIT predicted. What the fix removes is a **per-target dependence on which interpreter drew
the null of up to 0.128** — invisible in the mean, and exactly the class of bug BRIEF §3.7 exists to
prevent. Recorded as: *the aggregate was safe; the per-target values were not.*

---

---

## 1. VERDICT — leading with what most damages the programme's existing claims

> **1. Both "repair" operators make the structure WORSE than doing nothing, and the
> zero-information reference wins on accuracy.** The raw all-atom coordinate average is
> **3.0498 Å**. The ideal-geometry projection emits **3.2052 Å** (+0.1554 [+0.1217, +0.1911],
> **27W/99L**). AMBER k = 30 emits **3.1831 Å** (+0.1333 [+0.1031, +0.1645], **24W/99L**, gated
> n = 123). **The programme's −0.022 Å is the gap between two ways of being 0.13–0.16 Å worse
> than the operator that does nothing at all,** and each loses to it on 99 of 123 targets.
>
> **2. A RANDOM displacement of the SAME MAGNITUDE is at least as accurate as AMBER's.**
> Matched-magnitude isotropic control: **3.1606 Å against AMBER's 3.1831** (AMBER − random
> **+0.0225 [−0.0138, +0.0571], 43W/83L**), and **3.1826 against the projection's 3.2052**
> (+0.0226 [−0.0194, +0.0640], 48W/78L). Neither operator's *direction* is worth anything over a
> random direction of its own size.
>
> **3. AMBER's displacement points slightly AWAY from the truth, and significantly so relative
> to random.** ORACLE: cos(v_AMBER, r) = **−0.0521 [−0.0924, −0.0129]**, median −0.073, positive
> on only **37.3%** of targets, against a matched random control at −0.003 / +0.017 / +0.018.
> Paired, **AMBER − random = −0.0491 [−0.0996, −0.0009]** (i.i.d.), fold-clustered
> [−0.0876, −0.0124] — **a confidence interval excluding zero in the wrong direction.** The
> projection is −0.0371 [−0.0814, +0.0054], reproducing the coordinator's independent −0.0677.
>
> **4. There is NO (k, depth) at which the accuracy effect survives on the frame-reproducible
> subset.** Five settings, all reported: k = 30 full **−0.0095 [−0.0242, +0.0056] 39W/48L**
> (n = 87); depth 400 −0.0061; depth 100 **+0.0117**; depth 25 **+0.0109**; k = 0 **+0.1602
> [+0.0502, +0.2724], 21W/49L**. Every one either includes zero with a losing or near-even
> win/loss record, or is significantly worse. **The answer to the brief's question 3 is: there is
> none.**
>
> **5. The frame null FAILS at every converged setting.** The rotated-lab-frame difference is
> zero by rigid invariance. Gated, its *mean* is small at k = 30 (−0.00095 [−0.0060, +0.0042])
> and reproduces PHYS's ≈0.004 Å floor — but its **maximum is 0.1373 Å**, six times the effect
> being claimed, and **122 of 122 converged targets move by more than 1e−6 Å**. Against the
> standing PASS band (|mean| ≤ 0.005 **and** max ≤ 0.05) **every full-depth arm FAILS.** The only
> arm that PASSES is depth 25 (mean +0.00014, max **0.0082**, 20/50 non-zero) — and the
> convergence gate throws away **76 of 126 targets** there, so it is not a usable operator.
>
> **6. The VALIDITY claim is real, large, and entirely reproduced by a ZERO-INFORMATION constant
> α-helix.** AMBER k = 30 versus the projection: Ramachandran-favoured 0.466 → 0.874
> (−0.4077 [−0.4529, −0.3628], **116W/2L**), clashes 1.397 → 0.000 (**63W/0L**) — PHYS's numbers
> reproduce exactly. But a constant ideal α-helix scores **rama 1.000 and 0.000 clashes by
> construction**, and beats AMBER on Ramachandran **0W/65L** (+0.1264 [+0.1008, +0.1536]). A
> validity statistic that a zero-information reference maximises is not evidence about the force
> field; it is only meaningful **jointly** with staying near the input, and that joint statement
> is the one this workstream can defend.
>
> **7. The restraint constant k is an ACCURACY knob, not a VALIDITY knob.** Unrestrained AMBER
> (k = 0) reaches essentially the same validity as k = 30 — rama 0.868 vs 0.874, clashes 0.000 vs
> 0.000, geometry deviation 0.0170 vs 0.0192 — while costing **+0.2742 Å [+0.1785, +0.3823]
> (30W/93L)** against the projection and **+0.4309 Å** against doing nothing. The force field
> supplies the stereochemistry; the restraint supplies nothing but the decision not to move.

**Tiering.** 1, 2, 5, 6, 7 **DEMONSTRATED**. 3 and 4 **DEMONSTRATED** (3's cosines are **ORACLE
DIAGNOSTICS** — they read the native and are labelled so everywhere). Nothing here is a
hypothesis.

---

## 2. Q1 — THE (k, depth) TABLE, with gate exclusions and frame nulls

**n = 126.** ACCURACY is CA-RMSD to native (**ORACLE** evaluation). The gate excludes a target
when its frame-0 minimisation ends above 1000 kcal/mol; the excluded names are printed, not just
counted. `vs PROJECTION` and `vs DO-NOTHING` are gated; the ungated twins are in
`repair_report.json`. The FRAME NULL is measured on the converged-in-both-frames subset.

| setting | k | depth | AMBER (Å) | gate excl | vs PROJECTION | vs DO-NOTHING | wall s | moved (Å) |
|---|---|---|---|---|---|---|---|---|
| `k0_full` | 0 | full | 3.4775 | 3 | **+0.2742 [+0.1785, +0.3823] 30W/93L** | +0.4309 [+0.3328, +0.5374] 24W/99L | 13.51 | 1.526 |
| `k30_d25` | 30 | 25 it | 3.1192 | **76** | +0.0114 [−0.0125, +0.0331] 20W/30L | +0.0707 [+0.0478, +0.0968] 8W/42L | (see note) | 0.465 |
| `k30_d100` | 30 | 100 it | 3.1864 | 3 | −0.0196 [−0.0387, −0.0004] 63W/60L | +0.1371 [+0.1055, +0.1700] 20W/103L | 0.65 | 0.748 |
| `k30_d400` | 30 | 400 it | 3.1866 | 3 | −0.0191 [−0.0329, −0.0058] 62W/61L | +0.1376 [+0.1082, +0.1688] 22W/101L | 2.46 | 0.727 |
| **`k30_full`** | 30 | full | 3.1831 | 3 | **−0.0234 [−0.0375, −0.0099] 66W/57L** | +0.1333 [+0.1031, +0.1645] 24W/99L | 12.14 | 0.725 |
| *reference* | — | 0 it | **3.0498** | — | — | — | 0 | 0 |

The gate excludes the **same three targets — 1D6X, 2NB7, 7BX2 —** at every converged setting.
(PHYS listed four; the fourth, 1MF6, exceeds 1000 kcal/mol only in a rotated frame, and does so in
this run's frame-B arm too.) At depth 25 the gate excludes **76 of 126**, because 25 iterations is
simply not enough to leave the strained input.

**The k = 30 row reproduces PHYS exactly** — PHYS's gated figure was −0.02339 at n = 123; this run
returns **−0.0234 [−0.0375, −0.0099]**, fold-clustered [−0.0325, −0.0164], median **−0.0074**,
66W/57L, concentration **PASS** at mean/sd −0.296. The median-versus-mean early warning fires here
exactly as the ledger says it should: **the mean is 3.2× the median with a near-even win/loss
record**, and the null-calibrated concentration check has almost no power at that mean/sd, so it
settles nothing either way.

### 2.1 The frame null — zero by construction, and it FAILS

| setting | gated n | mean ΔRMSD | 95% CI | sd | **max \|Δ\|** | n with \|Δ\| > 1e−6 | verdict |
|---|---|---|---|---|---|---|---|
| `k0_full` | 122 | −0.00418 | [−0.0815, +0.0724] | 0.432 | **2.6502** | 122 | **FAIL** |
| `k30_d25` | 50 | +0.00014 | [−0.0001, +0.0005] | 0.0012 | **0.0082** | 20 | **PASS** |
| `k30_d100` | 122 | −0.00010 | [−0.0025, +0.0019] | 0.0125 | **0.0886** | 122 | **FAIL** |
| `k30_d400` | 122 | −0.00298 | [−0.0072, +0.0005] | 0.0215 | **0.1639** | 122 | **FAIL** |
| `k30_full` | 122 | −0.00095 | [−0.0060, +0.0042] | 0.0282 | **0.1373** | 122 | **FAIL** |

The band is `|mean| ≤ 0.005 AND max ≤ 0.05` on the converged subset, taken unchanged from
`s16/energy_lib.py` so that ENERGY's and REPAIR's nulls are judged by one rule. **Every arm passes
on the mean and every converged arm fails on the maximum.** This is the honest reading of a
quantity that is mathematically zero: its *central tendency* is a 0.004 Å floor, exactly as PHYS
found, and its *tail* is 0.14 Å at k = 30 and **2.65 Å unrestrained** — because the minimiser
chooses a different basin in a rotated frame on a minority of targets. The effect being claimed is
0.023 Å.

**Consequence, stated plainly.** By the pre-declared rule the brief gave me — *an arm whose frame
null is non-zero is not reportable* — **no full-depth AMBER accuracy number in this file is
reportable**, and the one arm whose null passes is one the convergence gate destroys. I report the
numbers anyway, with the failure attached, because suppressing them would hide the reason.

### 2.2 The frame-reproducible subset — the brief's question 3

Stratum: converged in both frames **and** |E(frame B) − E(frame A)| ≤ 0.5 kcal/mol (PHYS's
definition, native-free, unchanged).

| setting | n stable | effect vs projection | fold-clustered CI | median | W/L |
|---|---|---|---|---|---|
| `k0_full` | 70 | **+0.1602 [+0.0502, +0.2724]** | [+0.1116, +0.2061] | +0.0703 | 21W/49L |
| `k30_d25` | 47 | +0.0109 [−0.0136, +0.0353] | [−0.0049, +0.0243] | +0.0241 | 19W/28L |
| `k30_d100` | 31 | +0.0117 [−0.0165, +0.0367] | [−0.0162, +0.0367] | +0.0211 | 9W/22L |
| `k30_d400` | 71 | −0.0061 [−0.0223, +0.0083] | [−0.0271, +0.0173] | +0.0094 | 30W/41L |
| **`k30_full`** | 87 | **−0.0095 [−0.0242, +0.0056]** | [−0.0266, +0.0101] | **+0.0057** | **39W/48L** |

**There is no (k, depth) at which the accuracy effect survives.** The k = 30 row reproduces PHYS's
−0.009 [−0.024, +0.006] 38W/46L (n = 84) on an **independent frame draw** — mine is
−0.0095 [−0.0242, +0.0056] 39W/48L at n = 87. Every stable-subset median is **positive** (worse)
including k = 30's. Five settings were tested; at α = 0.05 one would expect a ~23% chance of a
spurious hit and none occurred, so no multiplicity correction is needed to reach the conclusion.

### 2.3 Depth is a validity–accuracy trade, and the accuracy optimum is "do not move"

| depth | RMSD | rama | clashes | geom dev | converged |
|---|---|---|---|---|---|
| 0 iterations (= do nothing) | **3.0498** | 0.836 | 5.611 | 0.2480 | — |
| 25 | 3.1192 | 0.639 | 1.706 | 0.1358 | **50/126** |
| 100 | 3.1864 | 0.806 | 0.000 | 0.0236 | 123/126 |
| 400 | 3.1866 | 0.860 | 0.000 | 0.0196 | 123/126 |
| full | 3.1831 | **0.874** | 0.000 | 0.0192 | 123/126 |

Deeper minimisation buys validity monotonically and costs accuracy from the do-nothing limit.
Depth 25 looks best on RMSD purely because it moves least (0.465 Å against 0.725 Å) — and it is
not a valid structure: 1.7 clashes, rama 0.639, and it fails the convergence gate on 60% of the
instrument.

### 2.4 Cost, quoted at the right number

Measured on this box, `threads = 1`, single process, contended (so these are upper bounds):

| operator | mean wall |
|---|---|
| ideal-geometry projection (`lam_path`, nine L-BFGS runs) | **4.44 s** |
| AMBER k = 30, depth 100 | 0.65 s |
| AMBER k = 30, depth 400 | 2.46 s |
| **AMBER k = 30, full** | **12.14 s** |
| AMBER k = 0, full | 13.51 s |

The k = 30 full arm is **≈2.7× the projection** and **≈10³× an AMBER single point (8–14 ms)**.
*Note on the depth-25 timing:* its 2.69 s figure is contaminated — it is the first setting run per
target in pass C and therefore pays the OpenMM builder construction; the ordering-robust depth
costs are the 0.65 s and 2.46 s rows. **No single-point cost is quoted for any minimisation arm in
this file.**

---

## 3. Q2 — WHAT THE OPERATORS ACTUALLY DO TO THE GEOMETRY

**Every quantity in this section reads the native and is an ORACLE DIAGNOSTIC.** None of it enters
a parameter, a threshold, a stopping rule or a hyperparameter. It is measured at the input the
operators are handed — the all-atom coordinate average, `||r||/√n = 3.0498 Å`.

### 3.1 Where the truth lives

| | loud half | quiet half | non-torsional | rigid |
|---|---|---|---|---|
| **the true residual r (ORACLE)** | **0.602** | 0.214 | 0.184 | 0.0000 |
| isotropic-direction null (analytic) | 0.302 | 0.333 | 0.365 | — |

The residual is **twice as concentrated in the loud half of the superposed-CA Jacobian as a random
direction is**, and only 18% of it is outside the reach of torsion motion. **The truth is a loud,
mostly-torsional displacement.**

### 3.2 What the operators do — CA space

Displacement `v` = operator output superposed onto its input, minus the input; the rigid part is
removed by the superposition (measured rigid share 0.0000).

| operator | cos(v, r) **ORACLE** | median | frac > 0 | \|v\|/\|r\| | loud | quiet | non-torsional |
|---|---|---|---|---|---|---|---|
| **AMBER k = 30 full** | **−0.0521 [−0.0924, −0.0129]** | −0.073 | 0.373 | 0.263 | 0.235 | 0.214 | **0.551** |
| AMBER k = 30, depth 100 | −0.0551 [−0.0954, −0.0141] | −0.078 | 0.325 | 0.280 | 0.280 | 0.256 | 0.463 |
| AMBER k = 30, depth 25 | −0.0502 [−0.0898, −0.0091] | −0.057 | 0.405 | 0.216 | 0.319 | 0.316 | 0.365 |
| AMBER k = 0 (unrestrained) | −0.0281 [−0.0880, +0.0294] | −0.043 | 0.476 | 0.602 | **0.567** | 0.226 | 0.207 |
| **ideal-geometry projection** | −0.0371 [−0.0814, +0.0054] | −0.068 | 0.397 | 0.290 | 0.237 | 0.244 | 0.518 |
| ZERO-INFO constant α-helix | −0.0423 [−0.1160, +0.0298] | −0.058 | 0.437 | 0.883 | 0.626 | 0.184 | 0.189 |
| RANDOM matched-magnitude ×3 | −0.003 / +0.017 / +0.018 | ≈0 | ≈0.5 | 0.263 | 0.297 | 0.279 | 0.385 |
| **isotropic-direction null (analytic)** | 0 | 0 | 0.5 | — | **0.302** | **0.333** | **0.365** |

**Three facts.**

1. **Direction.** Both repair operators have a **negative** cosine with the residual — they move
   slightly away from the native. AMBER's excludes zero on its own CI and is **significantly worse
   than the matched random control** (−0.0491 [−0.0996, −0.0009]; fold-clustered
   [−0.0876, −0.0124]). Against the ZERO-INFORMATION α-helix direction it is indistinguishable
   (−0.0098 [−0.0908, +0.0727], 64W/62L).
2. **Magnitude.** AMBER moves **0.725 Å RMS, 26% of the residual**; the projection moves 0.813 Å,
   29%.
3. **Where.** AMBER's motion is **enriched 1.5× in the non-torsional directions** (0.551 against
   an isotropic 0.365) and **depleted in the loud half** (0.235 against 0.302) — it moves where
   the residual is *not*. Unrestrained AMBER is the mirror image: 0.567 loud, and it is precisely
   the arm that loses 0.43 Å.

### 3.3 The same displacement in torsion space (right singular vectors, loud first)

| operator | cos(δθ, e_true) **ORACLE** | \|δθ\|/\|e\| | loud | quiet | exactly-null |
|---|---|---|---|---|---|
| AMBER k = 30 full | +0.0074 [−0.0373, +0.0513] | 0.534 | 0.242 | **0.488** | **0.270** |
| AMBER k = 30, depth 100 | −0.0051 [−0.0431, +0.0331] | 0.576 | 0.237 | 0.500 | 0.263 |
| AMBER k = 0 | +0.0278 [−0.0221, +0.0780] | 0.711 | 0.252 | 0.425 | 0.323 |
| ideal-geometry projection | +0.0301 [−0.0119, +0.0749] | 1.078 | **0.111** | 0.391 | **0.498** |
| ZERO-INFO constant α-helix | **+0.1380 [+0.0855, +0.1893]** | 0.767 | 0.223 | 0.274 | 0.504 |
| **isotropic-direction null (analytic)** | 0 | — | **0.381** | **0.421** | **0.199** |

**AMBER's torsional motion is depleted in the loud half (0.242 vs 0.381) and enriched in the
RMSD-quiet half (0.488 vs 0.421) and in the five exactly-null directions (0.270 vs 0.199).** It
moves, by preference, in directions Cα-RMSD cannot see. That is exactly what a stereochemistry
operator ought to do — it is a positive characterisation of what the force field is for — and it
is also why it cannot fix accuracy.

Note the one row whose cosine excludes zero: the **zero-information constant α-helix**, at
+0.1380. A fixed helix points at the truth in torsion space better than the force field does. This
is the same phenomenon AUDIT recorded (a constant α-helix reproduces 0.357 of the pool surrogate's
0.390) and it is a warning about torsion-space cosines, not a method.

### 3.4 The RMSD change is fully explained by the displacement alone

In coordinate space `n·RMSD_after² = ||r||² − 2 v·r + ||v||²` **exactly**, so a displacement
orthogonal to the residual can only cost, and the cost is fixed by its magnitude.

| operator | predicted | realised | \|pred − real\| | cost of \|v\| alone | corr(\|v\|, ΔRMSD) |
|---|---|---|---|---|---|
| AMBER k = 30 full | 3.1875 | **3.1831** | **0.0044** | +0.1177 | +0.464 |
| AMBER k = 0 | 3.5233 | 3.4775 | 0.0458 | +0.5157 | +0.476 |
| ideal-geometry projection | 3.2101 | **3.2052** | **0.0049** | +0.1515 | +0.507 |
| ZERO-INFO α-helix | 4.2040 | 4.0696 | 0.1343 | +1.7885 | +0.689 |

**The whole of the −0.022 Å is a magnitude difference.** If both operators' displacements were
exactly orthogonal to the residual, the projection would cost +0.1515 Å and AMBER +0.1177 Å — a
gap of **0.0338 Å**, against a realised gap of 0.0221–0.0234 Å. **AMBER "wins" because it moves
0.09 Å less, not because it moves better.** The two small negative cosines account for the rest.

### 3.5 The direction the quoted effect actually lives along

The −0.022 Å figure is a **difference of two displacements**, `v_AMBER − v_projection`. Its cosine
with the true residual is **−0.0258 [−0.0602, +0.0091]**, median −0.0309, positive on 41.3% of
targets. Its magnitude is 0.4471 Å RMS, and the two operators' displacements agree at
cos = +0.7317.

> **So the quantity the programme has been quoting as its only positive physics result is a
> 0.45 Å move in a direction with no measurable relation to the native, whose entire effect on
> RMSD is that it is 0.09 Å shorter than the alternative.**

### 3.6 Cross-check against an independently written module

`s16/csteer.py` measured the same projection direction from the CA-only coordinate average and
reported **cos = −0.0677**. Recomputed here from that same start: **−0.0737 [−0.1132, −0.0331]**,
n = 126. From the *all-atom* average's CA block — the object AMBER is actually handed, which
differs from it by ~0.16 Å — the same direction reads **−0.0371**. Two independently written code
paths agree on sign and order of magnitude; the residual difference is the start object, and it is
recorded here rather than left as a silent discrepancy.

**This also answers the coordinator's question directly.** `s16/CSTEER_FINDINGS.md` found every
native-free direction from the pool coordinate average to have a cosine with the true residual
indistinguishable from zero, and asked whether AMBER's displacement is a genuinely new native-free
direction. **It is not.** Its cosine is −0.052, on the wrong side of zero, worse than a matched
random direction, and indistinguishable from a zero-information α-helix. **That closes the
question cleanly rather than opening one.**

---

## 4. THE SHAPE CHECK IN k — pre-registered subsample, n = 30, NOT a selection

k ∈ {5, 15, 60} was run only on the **pre-registered, native-free, fold × length stratified
30-target subsample** listed in §0.6, which was fixed and printed before any pass-D RMSD existed.
k = 30 is shown on the same 30 targets so the rows are comparable. **Nothing is selected here and
no k below is re-quoted at full n.**

| k | AMBER (Å) | gate excl | moved (Å) | vs PROJECTION | stable subset | rama | clashes | geom dev | wall s |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 3.6998 | 2 | 0.971 | +0.0390 [−0.0000, +0.0799] 10W/18L | −0.0005 [−0.0406, +0.0375] n=15 | 0.794 | 0.000 | 0.0150 | 11.05 |
| 15 | 3.6547 | 2 | 0.880 | −0.0016 [−0.0302, +0.0269] 11W/17L | +0.0059 [−0.0266, +0.0348] n=16 | 0.812 | 0.000 | 0.0171 | 7.56 |
| 30 | 3.6306 | 2 | 0.817 | −0.0292 [−0.0628, +0.0024] 14W/14L | +0.0020 [−0.0406, +0.0377] n=17 | 0.822 | 0.000 | 0.0234 | 13.69 |
| 60 | **3.6037** | 1 | **0.735** | **−0.0533 [−0.0879, −0.0194] 18W/11L** | −0.0270 [−0.0819, +0.0187] n=14 | 0.843 | 0.000 | 0.0353 | 6.59 |
| *do nothing* | **3.5061** | — | 0 | — | — | 0.781 | 7.367 | 0.2885 | 0 |

**The shape is monotone and it is the magnitude law, not a physics optimum.** RMSD falls
monotonically in k (3.700 → 3.655 → 3.631 → 3.604) and the RMS displacement falls monotonically
with it (0.971 → 0.880 → 0.817 → 0.735). Every rung is explained by §3.4: the arm that moves least
loses least. **The limit of the trend is k → ∞, which is exactly the do-nothing structure at
3.5061 Å** — better than every rung. So "raise k" is not a route to accuracy; it is a slow approach
to not running the operator.

On the frame-reproducible subset, even k = 60 returns −0.0270 [−0.0819, +0.0187] at n = 14 — a
confidence interval including zero on a sample small enough that this programme's own record
(n ≤ 8 has reversed a conclusion six times, once with a CI that excluded zero) says it settles
nothing. **The frame null fails at every rung here too** (max |Δ| 0.05–0.14).

Validity rises slightly with k (rama 0.794 → 0.843, clashes 0 throughout) but the **residual
geometric strain rises with it as well** (rms relative bond/angle deviation 0.0150 → 0.0353),
which is the expected trade: a stiffer restraint leaves more of the input's broken geometry in
place.

---

## 5. WHAT THIS CHANGES IN THE PROGRAMME'S RECORD

1. **`s15/LEGACY_VS_AMBER.md` §5 and `s15/ARCHITECTURE.md`** still describe restrained relaxation
   at k = 30 as "worth −0.022 Å at valid geometry". PHYS already flagged this. It now needs three
   further qualifications, all measured here at n = 126: the comparison's own **zero-information
   reference beats both arms by 0.13–0.16 Å on 99 of 123 targets**; a **matched-magnitude random
   displacement is at least as accurate** as AMBER's; and the operator's displacement has a
   **negative** cosine with the true residual that is significantly worse than that random
   control.
2. **The "AMBER repairs stereochemistry" claim survives but must be stated jointly.** In
   isolation, every validity statistic quoted for it is maximised by a **constant α-helix that
   knows nothing about the target**. The defensible claim is the conjunction:
   *restrained relaxation reaches Ramachandran-favoured 0.874 and zero clashes **while moving only
   0.725 Å (26% of the residual) from its input** — the α-helix reaches the same validity at
   0.883 of the residual, i.e. by discarding the structure.* Validity alone is not an axis on
   which a force field can be credited.
3. **The convergence gate must be applied everywhere, and the exclusion is stable and small.**
   Three targets (1D6X, 2NB7, 7BX2) at frame 0, a fourth (1MF6) in a rotated frame. Gating changes
   the k = 30 effect from −0.0221 to −0.0234, i.e. it is free of direction risk, exactly as PHYS
   predicted.
4. **The frame null must be reported with its MAXIMUM, not only its mean.** A mean-only reading
   makes this pipeline look like it has a 0.004 Å floor. The tail is 0.14 Å at k = 30 and 2.65 Å
   unrestrained. The standing band in `s16/energy_lib.py` already encodes both, and **every
   converged arm fails it**.
5. **The restraint constant is not a validity knob.** Sprint 14 selected k = 30 on development-set
   RMSD from {0, 2, 10, 30, 100, 300} and PHYS recorded that the non-circular version of the rule
   picks k = 10 and emits +0.001 Å. The present result is stronger and simpler: **k trades
   displacement magnitude against nothing else that matters**, validity is achieved at every k
   including k = 0, and the accuracy-optimal end of the ladder is not running the operator.
6. **`s15/align_jac.py:136` no longer violates BRIEF §3.7.** The aggregate is unchanged; the
   per-target values moved by up to 0.128 depending on the interpreter (§0.7).

### 5.1 What I would say to the architecture

The mandated pillar is satisfiable and should be kept — but as a **terminal validity operator with
a declared accuracy price**, not as a refinement step. On this instrument that price is
**+0.133 Å against not running it**, bought for **12.14 s per target** (≈2.7× the projection it
replaces, ≈10³× an AMBER single point). If the architecture wants buildable, clash-free,
Ramachandran-favoured output, restrained AMBER at k = 30 is the cheapest thing measured here that
delivers it *while staying near the input*, and the ideal-geometry projection is strictly worse on
every axis except wall clock (rama 0.466, 1.40 clashes, and a **larger** accuracy price of
+0.155 Å). **The honest headline for the final report is that AMBER's repair role is real and its
accuracy role is not — and that the two must never again be quoted as one number.**

---

## 6. WHAT I REFUTED, INCLUDING MY OWN EXPECTATIONS

- **My own prior that a shallower minimisation would be the frame-reproducible sweet spot with
  intact validity.** Depth 25 does pass the frame null — the only arm that does — but it fails the
  convergence gate on **76 of 126 targets** and carries 1.7 clashes and rama 0.639. There is no
  depth at which the operator is both frame-reproducible and valid.
- **My own prior that AMBER's displacement would be uninformative but neutral.** It is
  *significantly worse than random at matched magnitude* on the direction axis
  (−0.0491 [−0.0996, −0.0009]). I expected a null and found a small negative.
- **The hypothesis that raising k would find an accuracy optimum.** The k-shape is monotone toward
  "do not move", and its limit is the do-nothing structure.
- **The framing in which the −0.022 Å is a physics result.** It is a magnitude difference: the
  orthogonal-move cost is +0.1177 Å for AMBER against +0.1515 Å for the projection, a gap of
  0.0338 Å against a realised 0.0221–0.0234 Å (§3.4).

---

## 7. LIMITATIONS AND WHAT REMAINS OPEN

- **Everything here is measured at ONE input structure** — the top-75 all-atom coordinate average
  at 3.0498 Å. Whether AMBER's displacement has a different geometry when applied to a *better* or
  a *differently wrong* start (a VQE-generated conformer, say) is untested here; the coordinator's
  `s16/integrate.py` operates on such ensembles and is the place that question belongs.
- **One frame draw.** PHYS ran three; I ran one (plus PHYS's three at k = 30, which agree). The
  frame null's *tail* is a heavy-tailed quantity and a single draw understates its variability —
  the programme's own record contains a 3-target probe of exactly this quantity that measured the
  wrong thing.
- **Pass D is n = 30.** Its k-shape is monotone and mechanistically explained, but no rung of it
  has been measured at full n, and it must not be used to select k.
- **The frame-instability stratum is still unexplained and unexploited.** PHYS left open whether
  membership is native-free-predictable; it still is. Nothing here fits a selector, deliberately —
  a selector fitted on which targets AMBER happens to help is a native-labelled hyperparameter.
- **The 60-target benchmark was not touched.** No arm, control or diagnostic in this workstream
  reads it.
- **No claim here rests on `I.FAIL18` membership**, and no rank-r first-order share is quoted as a
  finite-step budget.

---

## 8. REPRODUCTION

```
python -m s12.instrument                       # the five pinned constants
python -m s16.repair --mode verify --n 6       # bit-identity against s15/results/phys_repl_canon0.json
python -m s16.repair --mode prereg             # the pass-D subsample, printed before pass D runs
python -m s16.repair --mode pass --pass A      # k = 30, both frames, projection + controls, n = 126  (~60 min)
python -m s16.repair --mode pass --pass B      # k = 0,  both frames, n = 126                          (~40 min)
python -m s16.repair --mode pass --pass C      # depth 25/100/400 at k = 30, both frames, n = 126      (~20 min)
python -m s16.repair --mode pass --pass D      # k = 5/15/60 on the pre-registered n = 30 subsample    (~25 min)
python -m s16.repair_report A,B,C              # every Q1/Q2 table  -> results/repair_report_ABC.txt
python -m s16.repair_report A,D                # the shape check    -> results/repair_report_shape.txt
#   (repair_report.json is written by whichever pass set ran last; the canonical
#    n = 126 copy is repair_report.json and the n = 30 shape copy is
#    repair_report_shape.json)
python -m s16.repair_seedfix                   # the BRIEF S3.7 fix, verified
```

`s16/run_repair.sh` chains B → C → D behind A so that **only one heavy process is ever on the
box**. Threads pinned to 2 (`OMP/MKL/OPENBLAS`), OpenMM at `threads = 1` on the CPU platform with
`DeterministicForces`, `memo = False`, checkpointed every 5 targets with resume. All seeding via
`s15.seed.stable_rng`; **no `hash()` anywhere in this workstream.** Load was sampled before each
launch and the process count held at one; CPU stayed below the 97% ceiling throughout.
