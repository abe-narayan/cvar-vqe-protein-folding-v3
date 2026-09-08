# PHYSICS — STATUS to the coordinator (2026-09-06, in flight)

## COMPLETE at n = 126, and it is the answer to the sprint's physics question

**PHASE 6 — `leg_contact` as a term in the objective. MY OWN HEADLINE HYPOTHESIS IS REFUTED.**
`s18/results/lam.json`, `complete = true`, n = 126. Full tables in `s18/phys_FINDINGS.md` §1.

> **The pre-registered primary test:** `lam_c = +1` vs `lam_c = 0` (degree-1 distance only),
> identical start, fold-clustered — **+0.108 Å [−0.166, +0.377], median +0.199, 56W/70L**.
> Interval spans zero, median worse, losing record. **The pre-registered falsifier FIRES; the
> branch is closed.** No ladder extension, no re-normalisation, no sign flip.

Four things beside it that are worth more than the headline:

1. **The sequence information in MJ is worth nothing measurable here.** Against the
   zero-information control (identical objective, MJ table rebuilt on a permutation of the residue
   labels — same form, same pair set, same magnitudes) the primary arm is
   **−0.116 [−0.455, +0.164], 67W/59L**. What the term contributes is the *shape of a pair term*,
   not the amino acids in it.
2. **The pre-registered sign was right and no negative-λ rescue existed.** `lam_c = −1` — the sign
   the *global* ρ = −0.176 would have wanted — is **+1.314 [+0.958, +1.728], 27W/99L**.
3. **Adding the contact term destroys global ordering monotonically**: ρ(objective, RMSD) over the
   K = 500 pool goes 0.351 → 0.282 → **0.124** → **−0.060** as `lam_c` rises. The combination
   inherits `leg_contact`'s anti-ranking, not its information.
4. **Independent corroboration of your degree-1 closure, from my instrument.** `full` beats the
   degree-1 base by **−0.948 [−1.160, −0.753], 88W/38L**. And L30 reproduces: `full` is
   **+0.394 [+0.255, +0.523]** worse than the projected start (+0.558 against the raw average,
   against L30's +0.561).

**Cross-instrument check (G2b), which you can use.** My `full` arm and your `objceil` α = 0 arm are
the same objective, `dhat`, and start, in different code: **starts bit-identical (max \|Δ\| =
0.00e+00 Å)**, optima **3.606 vs 3.610, mean \|Δ\| 0.0066 Å**. The two instruments are one
instrument.

**A defect in my own pre-registration, reported rather than fixed.** The pool-sd normalisation says
"one sd against one sd", but the optimiser drives the two terms over different ranges: the
**realised share** of the objective's drop supplied by the contact term is **0.66 at `lam_c = +0.5`
and 0.85 at `lam_c = +1`**. The ladder therefore **never tested a small contact correction** — it
tested contact-dominated objectives. I am **not** re-normalising, because a scale chosen after
seeing an unfavourable answer is a tuned parameter in a methodological costume. It is filed OPEN.

## IN FLIGHT (box is at ~30%, so three of mine now run in parallel; BLAS pinned to 1 thread)

| | state | ETA |
|---|---|---|
| **`phys_decorr`** — the ORACLE go/no-go you asked for | **80 / 126**, ~10 s/target | **~8 min** |
| **`phys_down`** — PHASE 8, the causal downstream comparison | resumed at 18 / 126 | **~70–90 min** |
| **`phys_lambda --bases`** — the contact term on the **DEPLOYED** objective `E_full` and on MATH's `E_res` | started, 126 targets, ANOVA cached | **~30 min** |
| **`phys_space`** — PHASE 10 pass 1 (no AMBER) | queued last | ~25 min |

**Preliminary read on `phys_decorr` from its n = 12 smoke — and I am labelling it exactly as you
ask, "not measured", not "matched".** `corr(contact push, correction needed)` at the refinement
start was **mean +0.070, median +0.018**, and the difference from its MJ-shuffled null had a CI
spanning zero at n = 12. If that holds at n = 126 the answer to your question is the **middle**
case: `leg_contact` is **neither corrective nor amplifying — it is decorrelated from the
distogram's error**, which bounds what any λ could have bought and is consistent with the Phase 6
refutation. **I will not state it until n = 126 lands.**

**Preliminary read on `phys_down` from its n = 12 partial, same labelling caveat.** Every physics
filter placed in front of the averaging operator came out *worse* than no filter, and the
matched-random gate of the same count came out slightly *better*. If that survives to n = 126 it
falsifies the prediction I registered before the run — that Legacy's gate should HELP through an
averaging operator, because it improves the set MEAN (weight 1.16) and only damages the set BEST
(weight 0.04). **Not stated until n = 126.**

## Your two reminders, and where they are enforced

* **Own gated input**: `s18/phys_lib.gated` returns the gated comparison, its ungated twin and the
  exclusion count as **one object** — a caller cannot obtain the first without the other two — and
  `s18/phys_down.py` stores each arm's own input structure beside the arm. The trap cannot be
  sprung by inattention here.
* **`DEFAULT_WEIGHTS`, never fitted**: `leg_contact` enters at weight 1.0 exactly as
  `core.energy.DEFAULT_WEIGHTS` defines it, verified bit-for-bit against
  `s16.energy_lib.legacy_components_of_windows` by gate G0 (max \|ΔE\| = **1.78 × 10⁻¹⁵**).
* **Rotated-frame null with its maximum**: the only new AMBER operator in my lane is the incumbent
  k = 30 protocol used *as found*, whose frame null is already on the record; if I report a new
  AMBER operator it will carry its maximum. Phase 10's AMBER pass is gated on Pass 1 leaving an arm
  alive and is the lowest priority in my queue.

## Falsifier verdicts in my own words, so far

| falsifier | verdict |
|---|---|
| **H-C** (`leg_contact` succeeds as a term in the objective where it failed as a selector) | **REFUTED** at n = 126, on the arm and by the sign I pre-registered before running |
| the MJ table's sequence content matters in that role | **REFUTED** — indistinguishable from a residue-label permutation |
| a negative `lam_c` would rescue it | **REFUTED** — catastrophically worse, and it was pre-declared a control, never a rescue |
| my pool-sd normalisation means what its sentence says | **REFUTED by my own diagnostic** — realised share 0.66–0.85; a small contact correction is **untested**, filed OPEN, not quietly run |

---

# UPDATE — two more n = 126 artefacts complete

## A. The go/no-go you asked for: **DECORRELATED**, not amplifying — `s18/results/decorr.json`, n = 126, complete

**ORACLE DIAGNOSTIC.** At the structure the refinement starts from — the only place the two
gradients meet — the correlation between `leg_contact`'s per-pair **descent direction** and the
correction the distogram's error needs is

    corr(push, need)  =  -0.006   median -0.012   vs its MJ-shuffled null  +0.007 [-0.014, +0.029]   62/60
    over 50 pool candidates       -0.001

That is **zero on 126 targets with the interval tightened to +/-0.03** — not a small effect with a
wide interval. Your instrument cross-checks exactly: residual RMS **3.205 A**, MAE **2.347 A**.

**So the answer is your middle case.** `leg_contact` neither corrects nor amplifies the distogram's
error; it is orthogonal to it. **The ~1.0 A that your `shuffled` control shows is available from not
trusting the error's direction is NOT reachable through this term.** Whatever its
+0.080 [+0.032, +0.128] of in-band information is about, it is not about *where the distogram is
wrong*. Measured before the lambda ladder was read, and it predicts the ladder's result.

**One correction to the statistic, which changes the sign of the answer at the native.** The literal
"per-pair contribution" `c_p = MJ_p*switch(d_p)` and the force `push_p = MJ_p*|switch'(d_p)|` have
**different supports**: `c_p` is largest where the cosine switch is saturated and **the force there
is exactly zero**. Correlating the energy contribution with the residual weights pairs the term
cannot move. At the native the two disagree in sign — `corr(c,r) = +0.175` (reads as "amplifying")
against `corr(push,need) = -0.149`, and restricted to switch-active pairs +0.033. Both are reported;
the force statistic at the start is the one that governs an optimum, and it is zero. Also: at the
native the **MJ-shuffled null is MORE amplifying than the real table** (+0.035 [+0.015, +0.059]), so
the amplification there is a property of the pair geometry, not of the amino acids.

## B. The sharper form of my hypothesis, on the DEPLOYED objective — **significantly harmful**

`s18/results/lam_bases.json`, n = 126, complete. Same contact term, same pre-registered sign, same
native-free pool-sd normalisation, ladder cut to {0, +1}. Declared in the module header after your
alpha-ladder and MATH's correction, before it ran.

| base | `lam_c=0` | `lam_c=+1` | vs `lam_c=0` [fold CI] | W/L |
|---|---|---|---|---|
| **`E_full` (DEPLOYED)** | **3.606** | **4.329** | **+0.722 [+0.541, +0.893]** | **38/88** |
| `E_res` (MATH's object) | 4.265 | 4.808 | +0.543 [+0.125, +0.835] | 48/78 |

**On the objective that is actually deployed there is no ambiguity: adding `leg_contact` costs
+0.722 A with an interval excluding zero and a 38/88 record.** The pre-registered arm's +0.108 could
only be called "not supported"; this one is a positive finding of harm. It is not a property of one
base — the same intervention on `E_res` is harmful too.

**And the sequence information is worth nothing, now measured three times at n = 126.** Real MJ
against a residue-label permutation: **-0.033 [-0.151, +0.117]** (66/60) on `E_full`,
**+0.161 [-0.073, +0.428]** (62/64) on `E_res`, **-0.116 [-0.455, +0.164]** (67/59) on the
degree-1 base. All three span zero with near-even records.

**L30 reproduces at the level of the win/loss record.** `E_full` alone vs the coordinate average:
**+0.558 [+0.400, +0.680], 31W/95L**, against L30's +0.561 [+0.407, +0.712], **31W/95L**. Same mean
to 0.003 A and the identical split.

**One number for MATH/EXPERIMENT, not a PHYSICS claim.** Same start, same optimiser, n = 126:
**`E_res` 4.265 vs `E_full` 3.606 — +0.658 [+0.483, +0.828], 38W/88L.** On the enumerated 19-target
instrument the residue-additive object was -0.004 [-0.380, +0.307]; on the real continuous
instrument it is **significantly worse**. A fifth independent direction on the degree-1 closure,
from a module built for something else.

## STILL IN FLIGHT

| | state | ETA |
|---|---|---|
| **`phys_down`** — PHASE 8 | **36 / 126**, ~63 s/target | **~90 min** |
| **`phys_space`** — PHASE 10 pass 1 | just started | ~25 min |

If you need to close the sprint before Phase 8 lands, label it PARTIAL at its own n and I will send
the n = 126 table when it exists. Nothing from it is quoted anywhere yet.

---

# FINAL — ALL FIVE PHYSICS ARTEFACTS COMPLETE AT n = 126

`decorr.json` · `lam.json` · `lam_bases.json` · `down.json` · `space.json` — every one with
`complete = true`, a row count and a config hash. Smoke files renamed `_SMOKE_*` so none can be
read as a result. Full write-up: **`s18/phys_FINDINGS.md`** (verdict first, claim table, then the
five phases).

## PHASE 8 landed, and it is the hardest negative in my lane

**Every physics filter, placed in front of the operator that emits the answer, loses to a random
gate that removes the same number of candidates.** Identical structures, same top-75, f = 0.50.

| filter | vs NO FILTER | vs the MATCHED-RANDOM gate of the same count | W/L |
|---|---|---|---|
| AMBER single point | +0.050 [+0.011, +0.081] | **+0.036 [+0.006, +0.061]** | 52/74 |
| `leg_torsion` | +0.055 [+0.026, +0.081] | **+0.041 [+0.012, +0.068]** | 49/77 |
| `leg_contact` | +0.067 [+0.014, +0.119] | **+0.053 [+0.008, +0.097]** | 56/70 |
| Legacy (11-term) | +0.076 [+0.021, +0.138] | **+0.063 [+0.008, +0.125]** | 53/73 |
| **matched random** [CONTROL] | +0.013 [+0.003, +0.024] | 0 | — |
| constant α-helix [ZERO-INFO] | +0.139 [+0.069, +0.193] | +0.125 [+0.062, +0.181] | 44/82 |

Five scores, every interval excluding zero, every record losing. **Halving the set costs +0.013 Å
by itself; halving it by any physics score costs four to six times that. The cost is the ordering,
not the truncation.** Same sign and same ranking at f = 0.25.

**`leg_torsion` is closed by its own workstream.** Its recall gain reproduces
(**+0.083 [+0.019, +0.135]**, n = 55 targets that have a sub-2 Å member, one denominator) and
`s17/phys_FINDINGS.md` §8 left it open with *"has not been shown to convert into any downstream
accuracy"*. It is now measured: **it converts negatively.**

## A recorded law is wrong under a selective intervention — found by testing my own prediction

I registered before running that Legacy's gate should HELP, because
`d_out = 1.16·d_set_mean + 0.04·d_set_best` says the operator reads what the gate improves. **The
premise held exactly** (set mean 3.551 → 3.527, set best 2.306 → 2.561). **The conclusion was wrong
by the sign.**

| gate | law PREDICTS | OBSERVED | miss [fold CI] |
|---|---|---|---|
| `leg_torsion` | **−0.024** (better) | **+0.055** (worse) | **+0.079 [+0.060, +0.098]** |
| `legacy` | **−0.017** (better) | **+0.076** (worse) | **+0.094 [+0.063, +0.132]** |
| matched random | +0.008 | +0.013 | +0.006 [−0.003, +0.014] |

**The law holds for a random gate and fails for every score gate.** It is a correlational fit
across observed sets and does not survive a *selective* intervention. The natural mechanism is that
score-ordering **correlates the survivors' errors** — and the averaging operator's accuracy is
error cancellation (L30) — but I have not measured that, and the diversity column (`legacy` 2.213
vs `rand` 2.484 at the same m) is only consistent with it. **That is the one genuinely new question
this lane opened, and it is cheap.**

## AMBER: the instrument reproduces to the win/loss record

Against **its own gated input**, unfiltered arm: **+0.133 [+0.112, +0.165], 24W/99L**, against
Sprint 17's +0.1333 [+0.1031, +0.1645], **24W/99L** — same mean to 0.0003 Å and the identical
split. Convergence gate excludes the **same three targets, `1D6X 2NB7 7BX2`**, for the fifth
independent time. Validity on the same structures: clash < 2.0 Å **5.39 → 0.00**, < 2.6 Å
15.28 → 0.03, min heavy 2.116 → **2.809**, Ramachandran-favoured 0.837 → 0.875, cis 0.000 → 0.074.
**AMBER repair does not reorder the filters**: `legacy`→AMBER stays +0.078 [+0.016, +0.145] worse
than `none`→AMBER; against random-gate→AMBER the gaps shrink to null but never reverse.

## PHASE 10: closed, and it refutes this workstream's own Sprint-17 prediction

`s17` §8 item 4 predicted that repairing the coordinate average's 22.4% contraction before AMBER
would be cheap. Four native-free corrections, n = 126. The best — minimum-displacement projection
onto the 3.80 Å constraint set — restores the spacing **exactly** at **less** displacement than the
incumbent (0.801 vs 0.840) and still costs **+0.023 [+0.006, +0.037] more**. **Why**: a *random*
displacement of the projection's own magnitude costs ≈+0.154 Å and the projection costs +0.163 — a
null (+0.009 [−0.014, +0.033]). **The bill is for moving the Cα ~0.8 Å at all, not for choosing a
direction.** The incumbent projection's +0.163 [+0.132, +0.205] tax (L30's +0.164) is irreducible
by these means. Pass 2 (AMBER) was never spent, which is what the staging was for.

*(One arm did beat its matched random — `torsion_rebuild`, −0.154 [−0.209, −0.111] — and still lands
+0.785 Å worse than the average, because its move is 2.17 Å. A purposeful move that is three times
too big.)*

**Nothing further is in flight. My lane is complete.**
