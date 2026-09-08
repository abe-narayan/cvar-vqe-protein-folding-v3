# SPRINT 16 / REVIEWER — the report as a document

**Subject** `s16/FINAL_REPORT.md` as of 2026-09-06 07:16.
**Sources checked** `CLAIMS.md`, `LEDGER.md`, `PREREG_steer.md`, `STEER_FINDINGS.md`,
`CSTEER_FINDINGS.md`, `INTEGRATE_FINDINGS.md`, `audit_FINDINGS.md`, `lit_FINDINGS.md`,
`retract_FINDINGS.md`, `verify_FINDINGS.md`, `energy_FINDINGS.md`, `qphase_FINDINGS.md`,
`repair_FINDINGS.md`, artefacts in `s16/results/`.
**My own compute** `s16/review_checks.py` (three arithmetic/power checks, no heavy compute, no
benchmark read).

> **Verdict up front.** Almost every digit in the report traces to its source. But the report
> commits, in §5.5 and twice in §5.4, **the exact error it catalogues as failure mode #4 in §7** —
> it quotes row-level bootstrap intervals over 162 cells as though they were target-level intervals
> over 9 targets, and asserts in the same sentence that those intervals survived VERIFY. One of the
> four quantum arms does not survive. That single fact is enough for a hostile reader to conclude
> that the report's methodology section is decoration rather than practice, and it will be the first
> thing they find, because the report itself tells them where to look. **Not fit to hand over until
> R1–R4 are fixed.** The rest is repairable wording.

---

## R1 — WORST. §5.5 and §0.5 quote ROW-level intervals and claim they are target-level

**The report, §5.5:**

| arm | vs uniform random (as printed) |
|---|---|
| `cvar0.1` | **+0.212 [+0.059, +0.365]**, 8W/19L |
| `cvar0.25` | **+0.301 [+0.163, +0.437]**, 5W/22L |
| `cvar0.5` | **+0.234 [+0.087, +0.382]**, 8W/19L |
| `cvar1.0` | **+0.209 [+0.046, +0.357]**, 8W/19L |

followed by: *"These contrasts are paired within target, so they survived VERIFY's recomputation."*

**`verify_FINDINGS.md` §6.6, target as the unit, m = 75:**

    cvar0.25  +0.301 [+0.115, +0.485]  (7/2)
    cvar0.5   +0.234 [+0.011, +0.458]
    cvar0.1   +0.212 [+0.003, +0.404]
    cvar1.0   +0.209 [-0.002, +0.410]   <-- NOW TOUCHES ZERO

Three mismatches, in ascending order of seriousness:

1. **Every interval in the table is the row-level interval from `INTEGRATE_FINDINGS` §3**, not the
   target-level one VERIFY computed. `cvar0.25`'s is [+0.163, +0.437] in the report against
   [+0.115, +0.485] in VERIFY. Same for all four rows.
2. **Every W/L is row-level** (27 rows = 9 targets × 3 seeds). VERIFY's target-level W/L for
   `cvar0.25` is 7/2, not 5W/22L. The report's "5–8 wins against 19–22 losses" in §0.5 is a count
   of rows presented as a count of targets, on an instrument whose stated unit is the target.
3. **The claim of survival is false for `cvar1.0`.** §0.5 says "*+0.209 to +0.301, CIs excluding
   zero*"; at the target level `cvar1.0` is +0.209 **[−0.002, +0.410]**. VERIFY's own wording is
   "*the row bootstrap was **nearly** right here*" — the report upgraded "nearly right" to
   "survived".

This matters beyond the digits because §7 lists *"a **row** bootstrap over 162 cells → read as a
target-level interval over 9 targets"* as one of the five named failure modes of the sprint, and
§5.5 sits two pages later.

**Exact replacement, §5.5 table and the sentence under it:**

> | arm | vs untrained-circuit best-of-N (target-level) | vs uniform random (target-level) |
> |---|---|---|
> | `cvar0.1` | +0.115 (row-level CI only) | **+0.212 [+0.003, +0.404]** |
> | `cvar0.25` | **+0.204 [+0.028, +0.384]** | **+0.301 [+0.115, +0.485]**, 7 targets worse / 2 better |
> | `cvar0.5` | +0.137 (row-level CI only) | **+0.234 [+0.011, +0.458]** |
> | `cvar1.0` (plain VQE) | +0.112 (row-level CI only) | +0.209 **[−0.002, +0.410]** |
>
> **Three of the four CVaR-VQE arms lose to uniform random sampling at the same budget with a
> target-level interval excluding zero; plain VQE's interval touches zero.** These contrasts are
> paired within (target, seed), so the row bootstrap was *nearly* right here — the widening that
> destroyed the AMBER interval is small for them, but it is not zero, and the target-level intervals
> are the ones quoted.

**Exact replacement, §0 bullet 5, first sentence:**

> **The quantum component loses to uniform random sampling at matched budget** (+0.209 to +0.301;
> three of four arms with target-level CIs excluding zero, the fourth touching it; 7 of 9 targets
> worse on the largest arm), reconfirming the standing result on a terminal operator it had never
> been tested against.

---

## R2 — §0's "no surviving positive accuracy effect at all" is wrong, and the report says so itself

§0, paragraph 1: *"after its internal verification pass, **no surviving positive accuracy effect at
all**."*

§0, bullet 4, twelve lines later: *"**AMBER does contribute, but only in aggregate**: −0.0234
[−0.0375, −0.0099] at n = 126."*

§0, closing paragraph: *"The single accuracy effect that survives its own statistics at full
instrument size is **AMBER's −0.023 Å**."*

**These are not consistent, and the second and third are right.** The first is inherited verbatim
from `INTEGRATE_FINDINGS`'s correction box — *"the sprint has no surviving positive accuracy
effect"* — which was written by VERIFY about **the 9-target enumerated instrument only**. VERIFY's
own limitation §8 says so: *"The AMBER question is open, not answered negatively — the right
response is to measure it on the 126-target instrument."* ENERGY then did exactly that, at n = 126,
with a matched random-filter control, a pre-declared convergence gate reported both gated and
ungated, and a null-calibrated concentration check that passed at the 58.5th percentile. That is a
surviving positive accuracy effect. Carrying VERIFY's sentence out of its instrument and into the
report's opening line is the same category error the report is about.

Note also that §0's blanket phrasing erases §0 bullet 7 (the coordinate readout, +0.446
[+0.193, +0.706], 8/9 targets, which survived every control VERIFY threw at it).

**Exact replacement, §0 paragraph 1, second sentence:**

> It produced **no method that beats the incumbent**, and its 9-target integration result did not
> survive the target as the unit. **One accuracy effect survives at full instrument size** —
> AMBER's −0.0234 Å [−0.0375, −0.0099] at n = 126 — and it is 0.7% of the gap to 2.5 Å, not
> resolvable per target, and sits inside a rigid-invariance null that fails its own pre-declared
> band. What the sprint produced instead is a much sharper and much less flattering picture of where
> the programme's accuracy actually comes from, and two closed directions that were previously open
> and plausible.

---

## R3 — the correction that created §9 was not carried through: "one bit" is one bit **plus a magnitude**

`verify_FINDINGS.md` §1.5 writes the replacement text for `CSTEER_FINDINGS` §2 and ends it:

> *"The missing ingredient is not the direction. It is the **per-target sign** — one bit, **plus the
> step magnitude**."*

and §1.5's preceding paragraph is explicit: *"the sign of `c` **and the magnitude `||r||`** are both
native-derived, so every 'ORACLE two-sided ceiling' above is a ceiling."* `CSTEER_FINDINGS` §2
repeats it: *"Every 'ORACLE ceiling' above uses the native to place **both the sign and the
magnitude** of the step."*

The report drops the second half in both places it matters:

- §0 bullet 2: *"what is missing is the **per-target sign** — one bit."*
- §9: *"Everything needed to reach it is known except **one bit per target**: whether this peptide's
  averaged structure needs contracting or expanding."*

The optimal step is `t* = c‖r‖`. `‖r‖` is the RMSD to the native. **Nothing in this sprint estimates
it native-free, and the 2.491 Å ceiling is computed with it given.** Sizing the open direction as a
one-bit classification problem — which §9 explicitly leans on (*"it is a **classification** problem
on one bit rather than a regression problem on a structure"*) — is the single most consequential
overclaim left in the report, because it is the sentence that determines what the programme does
next.

**Exact replacement, §0 bullet 2, final clause:** replace *"and what is missing is the **per-target
sign** — one bit."* with:

> and what is missing is the **per-target sign** — one bit — **and the step magnitude `‖r‖`, which
> is the RMSD to the native and which nothing in this sprint estimates native-free.**

**Exact replacement, §9, "Open, and precisely stated" paragraph:**

> **Open, and precisely stated.** From the pool coordinate average at 3.048 Å, a **zero-information**
> constant β-strand direction has an ORACLE two-sided ceiling of **2.491 Å**. That ceiling is
> computed with **two** native quantities given: the per-target sign of `c`, and the step magnitude
> `t* = c‖r‖`, in which `‖r‖` is the distance to the native. The open question is therefore two
> questions, and only the first is a one-bit classification: (i) can any native-free signal *outside*
> the distance-matching family beat 0.516 at predicting the sign, and (ii) what does the 2.491 Å
> ceiling degrade to when the magnitude is estimated rather than given? **The second is untested and
> the ceiling is not meaningful without it** — a leave-fold-out constant step is the cheapest first
> answer and was not run.

**Exact replacement, §9, "Why it is worth asking" first clause:** replace *"and because it is a
*classification* problem on one bit rather than a regression problem on a structure"* with:

> and because its harder half is a *classification* problem on one bit rather than a regression
> problem on a structure — though the easier-looking half, the step magnitude, is a regression and
> is untested

---

## R4 — QPHASE has reported and the report contains none of it, including the control that would strengthen its own quantum claim

`s16/qphase_FINDINGS.md` is on disk (28 KB, mtime 07:03) and the report was written at 07:16. §5.8
is still an italic placeholder, §8 forward-references a §5.9 that does not exist, and the preamble
promises *"Sections marked **PENDING**"* while **no section in the document carries the word
PENDING**. Beyond the housekeeping, three substantive consequences:

1. **The report leads its quantum section with a control its own sibling workstream disqualifies.**
   `qphase_FINDINGS` §0 table: *"uniform random (**WEAK CONTROL, proves nothing**)"*. The report's
   §0.5 and §5.5 headline is "loses to uniform random sampling".
2. **The strong version of the result is sitting unused.** QPHASE §2b: the VQE is significantly
   worse than **greedy 1-opt with restarts at 10 of 10 objective-quality rungs**, is never
   significantly better than **simulated annealing at matched budget**, and is **statistically
   indistinguishable from annealing at a quarter of its objective budget at 9 of 10 rungs**. §2c:
   *"Greedy 1-opt reaches the certified global optimum in 68–100% of cells at 8,192 evaluations; the
   CVaR-VQE reaches it in 0–32%."* §3: greedy also beats the VQE on the **ensemble** readout at
   every rung, which retires Sprint 15's *"the one place a VQE wins is an ensemble consumed without
   a ranker"*. This is a materially stronger and better-controlled refutation than the one the
   report prints, on 19 targets rather than 9. **Underclaimed by omission.**
3. **There is a directional tension the report should not discover in review.** On the *real*
   objective (`signal = 0`), QPHASE's argmin readout gives vqe − untrained = **−0.309** (VQE
   *better*, CI [−0.649, +0.037], null) and its ensemble readout gives **−0.241** (VQE better,
   null), where `INTEGRATE` gives +0.112 to +0.204 (VQE *worse*). The readouts and controls differ —
   best-of-N argmin vs a 5-member coordinate average vs a 75-member one — and that is a legitimate
   explanation. But the report currently states the conclusion universally (*"running the VQE is
   worse than not running it"*) with no acknowledgement that a sibling workstream measured the sign
   the other way on the untrained-circuit control at n = 19. **State the readout dependence
   explicitly or a reader who opens `qphase_FINDINGS.md` will find it for you.**

Also missing and relevant to §9: QPHASE's *"the deployed `hamil` objective's **certified argmin is
2.661 Å** against a space best of 1.030 Å, and a classical 1-opt reaches it exactly in 100% of
cells"* — a ceiling on the enumerated instrument better than the incumbent's 3.204 Å, and directly
on the report's "what next" axis.

**Minimum fix:** either write §5.8 from `qphase_FINDINGS.md` before handing over, or add to the
placeholder line: *"`s16/qphase_FINDINGS.md` has reported and is not yet integrated; its headline —
the VQE loses to greedy 1-opt at 10 of 10 objective-quality rungs and is matched by annealing at a
quarter of its budget — is stronger than §5.5's and supersedes §5.5's choice of uniform random as
the control."* Also delete the "Sections marked PENDING" promise or actually mark them.

---

## R5 — §5.7's enumerated-space bullet compares three numbers across two populations and calls them matched; on the matched set the Legacy/AMBER direction flips

**Report §5.7, final bullet:** *"**On 19 complete enumerations with matched populations**, Legacy's
certified argmin is **+0.082 Å worse than random** (native at the 53.5th percentile); AMBER's is
−0.075 (40.7th); and the **distance objective is −0.926 [−1.440, −0.431]** (15W/4L, native at the
29.1st percentile), beating AMBER by −0.815 [−1.357, −0.310] on the identical masked set."*

`energy_FINDINGS.md` §5, which is careful about this and prints both rows:

| arm | population | n_pop | argmin − mean | CI 95% | W/L |
|---|---|---|---|---|---|
| Legacy | **full enumeration** | 676,056 | **+0.082** | [−0.259, +0.401] | 5/14 |
| Legacy | **binding-masked** | 933 | **−0.184** | [−0.523, +0.184] | 12/7 |
| AMBER | **binding-masked** | 933 | **−0.075** | [−0.540, +0.358] | 7/12 |
| distance objective | binding-masked | 933 | −0.890 | [−1.316, −0.464] | 14/5 |
| distance objective | full enumeration | 676,056 | **−0.926** | [−1.440, −0.431] | 15/4 |

Three failures in one sentence:

1. **"Matched populations" is false for the numbers quoted.** Legacy's +0.082 and the distance
   objective's −0.926 are on the full enumeration; AMBER's −0.075 is on the 933-configuration
   binding-masked stratum.
2. **The direction reverses when you actually match.** On the identical masked set Legacy's
   certified argmin is −0.184 and AMBER's is −0.075: **Legacy's optimum is better**, and
   `energy_FINDINGS` prints the head-to-head as *Legacy − AMBER = −0.109 [−0.666, +0.462]*. The
   report's §0 bullet 4 headline is *"Legacy versus AMBER — **measured twice, on two instruments,
   with the same answer**"* and *"AMBER ranks better than Legacy"*. There is a **third** instrument
   in the report's own §5.7 on which the point estimate goes the other way, and the report cites
   that instrument only for the halves that agree with it.
3. **Selective CI reporting.** Every interval printed in that bullet excludes zero; every interval
   omitted includes it. Legacy [−0.259, +0.401]; AMBER [−0.540, +0.358]; and — for the very next
   sentence — the constant-α-helix advantage over Legacy's optimum, **+0.469 [−0.006, +1.041], W/L
   5/6 with 8 exact ties**, i.e. a coin-flip on the 11 non-tied targets. `energy_FINDINGS` §5.1
   states plainly *"The native sits at chance for both. Legacy 0.535, AMBER 0.407 — **neither CI
   excludes 0.5**."* The report prints both percentiles as if they were results.

**Exact replacement for the whole bullet:**

> - **On 19 exhaustively enumerated spaces, neither energy ranks, and they are indistinguishable
>   from each other.** On the identical binding-masked population Legacy's certified argmin is
>   −0.184 Å [−0.523, +0.184] against a random feasible member and AMBER's is −0.075 Å
>   [−0.540, +0.358]; head to head, Legacy − AMBER = −0.109 Å [−0.666, +0.462]. **Every one of those
>   intervals contains zero**, and the ORACLE native percentiles (Legacy 0.531, AMBER 0.407) do not
>   exclude chance either. On the full enumeration Legacy's argmin is +0.082 Å [−0.259, +0.401]
>   *worse* than random, reproducing Roget et al., and it must be attributed to them. **The one arm
>   that breaks the pattern is not an energy**: the distance-restraint objective's certified argmin
>   is **−0.926 Å [−1.440, −0.431]** better than random (15W/4L, native at the 29.1st percentile,
>   CI excluding chance), beating AMBER by −0.815 [−1.357, −0.310] on the identical masked set. That
>   localises the pathology to the **potential class**, not to short-chain structure prediction.
>   A zero-information constant α-helix reaches an RMSD 0.469 Å below Legacy's certified optimum
>   ([−0.006, +1.041], 5W/6L with **8 exact ties** — i.e. on 8 of 19 targets Legacy's optimum *is*
>   the constant helix), which is why "beats uniform random" proves nothing on this instrument.

---

## R6 — three internal contradictions and one unfootnoted row-level interval in §5.4/§5.7

1. **§5.4** bullet: *"on small ensembles it [Legacy] is **significantly** worse than dropping members
   at random"* (from the bolded **+0.068 [+0.012, +0.125]** row). VERIFY §4.1: at target level this
   becomes **+0.0684 [−0.0389, +0.1757]**. The asterisk footnote (*"row-level intervals; the
   target-level recomputation widens them"*) is present, but the word **significantly** is used
   anyway, and the row is bolded as the only positive Legacy finding. Delete "significantly" and
   print the target-level interval in the table.
2. **§5.4** bullet: *"The composition adds nothing over AMBER alone (interaction −0.004
   [−0.042, +0.032])."* That is the row-level interval and it carries **no asterisk**. Target level
   (VERIFY §4.1): **−0.0044 [−0.1040, +0.0966]**. Replace.
3. **§5.4 table**, AMBER row: the CI is target-level (**[−0.148, +0.040]**) but the **median −0.061
   is row-level**; the target-level median is −0.045. Mixing units inside one row of one table.
   Replace `−0.061` with `−0.045`.
4. **§5.4 vs §5.7 contradiction.** §5.4 bolds a Legacy accuracy CI that excludes zero; §5.7 states
   *"**Legacy's only intervals excluding zero anywhere in the sprint** are on the validity axis"*.
   After fix (1) the contradiction disappears; without it the report contradicts itself across three
   pages.

Also §5.7: *"improves the **ORACLE** set mean (3.551 → 3.531) while destroying the set best
(2.306 → 2.416)"* — the set **best** is equally an ORACLE quantity and is not labelled. Add the
label to both.

---

## R7 — §8's seeding limitation is stale, incomplete, and understates its own consequence

**Report §8:** *"**One rule violation remains in the Sprint 15 tree**: `s15/align_jac.py:136` seeds
with `hash()`. The affected value was verified correct analytically, so **no number changes**."*

Three problems, all in `audit_FINDINGS.md` §5.2 and `repair_FINDINGS.md` §0.7:

- **It does not remain.** REPAIR fixed it to `SD.stable_rng("align_jac", "subspace_null", p)` and
  verified the fix (`s16/repair_seedfix.py`): OLD mean over the 7 reachable salts 0.4961, NEW
  0.5005, analytic 0.5000.
- **There were two, not one.** AUDIT §5.2 also flags `s15/info_regime.py:59`
  (`np.random.default_rng(hash(p) % 2**31)`), and records that it is unfixed.
- **"No number changes" is true of the aggregate and false of the per-target values.** REPAIR
  measured the per-target spread across interpreter salts under the old seeding at **mean 0.0759,
  max 0.1280**, and recorded the finding as *"the aggregate was safe; the per-target values were
  not."* The report keeps only the reassuring half.

**Exact replacement:**

> - **Two seeding-rule violations were found in the Sprint 15 tree, one fixed and one not.**
>   `s15/align_jac.py:136` seeded the random-subspace null with `hash()` (salted per process, then
>   collapsed onto 7 seeds); REPAIR fixed it to `stable_rng` and verified that the **aggregate does
>   not move** (0.4961 → 0.5005 against an analytic 0.5000) while the **per-target value moved by up
>   to 0.128 depending on which interpreter drew the null** — the aggregate was safe, the per-target
>   values were not. `s15/info_regime.py:59` carries the same violation and is **not yet fixed**; it
>   draws only the "one random top-75 window" channel, which `align_jac`'s reported list does not
>   include.

---

## R8 — two small arithmetic defects, one of them in a headline sentence

1. **§5.3: `c = 0.615` for 2.5 Å is wrong; it is 0.625.** `‖r_new‖ = ‖r‖√(1 − c²)` from 3.204 Å
   gives c = **0.6254** for 2.5 Å and 0.7812 for 2.0 Å (`s16/review_checks.py`). 0.615 is the value
   for a **3.170 Å** base. The error originates in `csteer.py`'s docstring and print statement,
   where it is hardcoded rather than computed, and propagated to `CSTEER_FINDINGS` §2 and the
   report. It is small, but it sits in the sentence *"which turns the sprint's target into
   arithmetic"* — an expert who checks one number in this report will check that one.
   **Fix in `s16/csteer.py:22`, `s16/csteer.py:205`, `CSTEER_FINDINGS.md` §2 and `FINAL_REPORT.md`
   §5.3: `2.5 Å needs c = 0.625`.**
2. **§5.3 conflates the analytic null with the sampled control.** The prose says the random control
   is *"at 0.128 — exactly the isotropic null √(2/(π·3n)) at n = 13"*; the table in the same section
   prints the random arm's `|c|` as **0.136**. √(2/(π·3·13)) = 0.1278; 0.136 is the sampled value.
   **Replace the prose with:** *"against the module's own random control, whose sampled |c| is 0.136
   and whose analytic isotropic null √(2/(π·3n)) at n = 13 is 0.128 — the paired subtractions below
   use the 0.128 control mean."*

---

## R9 — UNDERCLAIMED: the strongest positive result in the sprint is a clause inside a bullet

The **distance-restraint objective's certified argmin, −0.926 Å [−1.440, −0.431], 15W/4L, ρ = +0.485
positive on 89% of targets, native at the 29.1st percentile with a CI excluding chance**, measured
on 19 *complete* enumerations (12.8 × 10⁶ configurations, so the argmin is certified, not searched),
is:

- the only result in the sprint with a CI excluding zero on an **argmin** readout;
- the control that turns "both energies fail" from a possible property of 9–16-mers into a
  demonstrated property of **contact-style potentials** — `energy_FINDINGS` §5.1 calls it *"the
  strongest energy result available"*;
- consistent with, and an independent confirmation of, the standing law *a structural objective
  beats both energies*;
- and it appears in the report **only** as a subordinate clause in the third bullet of §5.7, absent
  from §0 entirely, absent from §6 ("what the sprint now believes"), and absent from §9.

This is the report's clearest instance of the opposite failure to overclaiming. It is also the one
result a hostile reader would most want promoted, because it is the only place the sprint
*distinguishes* two hypotheses rather than refuting one.

**Suggested addition, §0, new bullet after bullet 4:**

> 4b. **The energies' failure is a property of the potential class, not of short peptides.** On 19
>    *exhaustively enumerated* spaces (12.8 × 10⁶ configurations, so every argmin is a certified
>    optimum), a distance-restraint objective's argmin is **−0.926 Å [−1.440, −0.431]** better than a
>    random feasible member, 15 of 19 targets, with the native at the 29.1st percentile and a CI
>    excluding chance — while neither Legacy nor AMBER can be distinguished from random on the same
>    configurations. This is the sprint's only interval excluding zero on an argmin readout, and it
>    is the control that makes "neither energy ranks the native" a statement about contact-style
>    potentials rather than about 9–16-mers.

Two smaller under-claims: **Legacy's detection AUROC of 0.93–0.99 on the non-tautological `torsion`
component** is the best-performing measurement in the sprint and gets one sentence; and §5.6's
*"One narrow positive"* labels **3.051 Å against a shipped 3.048 Å** a positive when it is a tie
(and slightly worse), with no interval given — call it a tie, not a positive.

---

## R10 — `CLAIMS.md` carries two stale row-level intervals that contradict the report

`CLAIMS.md` is the register the report points readers to in §5. Two rows were not updated when
VERIFY recomputed:

| row | `CLAIMS.md` says | should say (VERIFY §4.1, §6.5) |
|---|---|---|
| **C7** | +0.446 **[+0.356, +0.537]**, **124W/38L** at m = 75 | +0.446 **[+0.193, +0.706]**, **8/9 targets** |
| **C2** | +0.068 **[+0.012, +0.125]** at m = 5, 70W/90L | +0.068 **[−0.039, +0.176]**, target level; 70W/90L is a row count |

C7 is the sprint's one surviving accuracy-relevant architectural claim, and the register quotes the
interval the sprint has already disowned — while the report quotes the correct one. Fix `CLAIMS.md`.

---

## R11 — no statement of power anywhere, on an instrument the report calls a null five times

The report declares NULL on the 9-target enumerated instrument for: Legacy at m = 75 (§5.4), the
Legacy→AMBER interaction (§5.4), α's monotone trend (§0.5, §5.5), and — in `CLAIMS.md` — C4 and C6.
§8 says only *"Nine targets is inside the range where this programme has repeatedly had conclusions
reverse."* That is a mood, not a number.

From VERIFY's own per-target AMBER effects (`verify_FINDINGS.md` §4.2), recomputed in
`s16/review_checks.py`:

    n = 9, sd = 0.153 A, SE = 0.051 A
    minimum detectable effect at 80% power:  ~0.161 A
    smallest effect whose CI excludes zero:  ~0.118 A
    power to detect a true 0.05 A effect:    ~9%
    power to detect a true 0.10 A effect:    ~37%

**Every "null" on this instrument is uninformative about effects below roughly 0.16 Å**, which is
larger than every effect the instrument was built to measure except the generator contrast (0.30 Å)
and the readout contrast (0.446 Å) — the two that *did* come out significant, which is exactly what
you would expect if power, not physics, were doing the sorting.

**Suggested addition, §8, first limitation, after the existing text:**

> Quantified: on the nine enumerated targets the per-target sd of the AMBER contrast is 0.153 Å, so
> the standard error of a target-level mean is 0.051 Å, the minimum detectable effect at 80% power
> is **≈0.16 Å**, and the power to detect a true 0.05 Å effect is **≈9%**. **Every null reported on
> this instrument — Legacy at m = 75, the Legacy×AMBER interaction, the absence of an α trend — is
> therefore an absence of evidence at effect sizes below ≈0.16 Å, not evidence of absence.** The two
> contrasts that did clear it (0.30 Å generator, 0.446 Å readout) are the only two large enough for
> it to clear.

---

## R12 — smaller labelling and housekeeping items

- **§0 bullet 2** presents `|c| = 0.305 / 0.398` as measurements without noting that `c` is a cosine
  **with the true residual**, i.e. an ORACLE evaluation of a native-free direction. §5.3's table
  labels only the *ceilings* ORACLE. Add "(ORACLE evaluation of a native-free direction)" to the
  §5.3 table header and to §0 bullet 2.
- **§5.1** quotes *"the full true error bought −3.604"* without an ORACLE label in that sentence;
  `STEER_FINDINGS` labels the arm `ORACLE_true`. Add it.
- **Preamble** promises *"Sections marked **PENDING**"*; no section is so marked. Either mark §5.8
  and §5.9 **PENDING** or delete the promise.
- **§8** forward-references *"the REPAIR workstream (§5.9)"*. `repair_FINDINGS.md` currently contains
  method and verification only — *"(Results sections follow as each pass completes.)"* — so §5.9 has
  no results to import yet. Say that rather than referring to a section that does not exist.
- **§2.1** says AUDIT found *"no transcription error anywhere in the record"*, which is right; but
  AUDIT §1.7's most quotable number — *"the pool channel is worth 0.095 Å of a 3.360 Å ORACLE
  ceiling, **2.8%**, with the step and the sign already given away, while an ORACLE line search along
  a **random** direction is worth eight times more"* — is not in the report. It is the sharpest
  single sentence in the provenance section. Consider adding it to the §2.1 table's third row.
- **AUDIT §7 left one control unrun** and the report does not record it: *"The 0.497 loud-half |cos|
  is measured against a magnitude-matched null but **not** against the zero-information α-helix
  reference. That control is the first thing the flagship should run."* The flagship closed the line
  a different way (Rule 4), which is defensible, but the un-run control belongs in §8.

---

# Q3 — Is the narrative honest about its own shape?

**Mostly load-bearing, with one instance of a rule confessed and not applied, and three
under-followed-through retractions.** The distinction that matters to a hostile reader is: *did the
correction change what the programme believes, or does it only change what the programme says about
itself?*

## Load-bearing — these changed beliefs, and would be defensible under attack

| correction | what it changed |
|---|---|
| **AMBER at n = 9 recomputed target-level** (E3 → C1 NOT ESTABLISHED) | withdrew the sprint's headline positive result; changed `CLAIMS.md` C1's status; changed §0. Cost real. |
| **the signed-cosine → \|c\| retraction** (E7) | reversed the sprint's *closing question*. Before: "the direction does not exist, the line is closed." After: "the magnitude is 2–3× the null, one bit is missing, ceiling 2.491 Å." **§9 exists only because of this correction.** This is the most load-bearing correction in the document. |
| **the fusion law → Krogh–Vedelsby** (B6) | destroyed a claimed novelty, produced a citation obligation, and produced the transferable rule in §7 that a transcription audit cannot see a formula error. |
| **the diversity 52/48 re-attribution** (E11, D2 → D4) | **reinstated a project law the sprint had declared overturned**. That is the shape of a correction that costs something: it took back a *win* the sprint had claimed over the existing record. |
| **AUDIT's zero-information controls** (B2, B4, B3) | killed three inherited premises and set the sprint's standing control requirement. Changed what may be built next. |
| **§5.7's `lam_path` correction** | ENERGY found a defect that shrinks **its own** headline 2.39×, and reported it. Straightforwardly costly and correctly foregrounded. |

## Decorative — one, and it is the one that matters

**§7's failure-mode table is a confession the report then fails to act on.** Row 4 of that table is
*"a **row** bootstrap over 162 cells → a target-level interval over 9 targets → the sprint's only
positive accuracy effect."* The report then quotes row-level intervals **four more times** — the
whole of §5.5's table (R1), §5.4's interaction bullet, §5.4's Legacy m = 5 bullet, and §5.4's median
column (R6) — three of them without the asterisk, one of them under an explicit assertion that
VERIFY's recomputation was survived. A hostile reader will read §7, turn back two pages, and find
the rule broken in the section it was derived from. **That single juxtaposition converts the
methodology section from an asset into the report's biggest liability**, and it is the reason for the
overall verdict.

A second, milder instance: **E8 (the coordinate-linearity "cleanest exhibit" demoted to a
tautology)** cost nothing, because the exhibit was never load-bearing for any conclusion — the
operational verdict of §5.3 is unchanged by it. It is a correct and honest demotion, but listing it
as one of five instances in §7 inflates the apparent price of the sprint's self-correction. It should
be kept and *labelled* as costless: "*a free confession — nothing rested on it*".

## Under-followed-through — retractions whose consequences did not reach the conclusions

1. **R3, the largest.** VERIFY's replacement sentence said "one bit **plus the step magnitude**". The
   report kept "one bit" and built §9's forward plan on it. The retraction was accepted; its
   consequence for the *next sprint* was not.
2. **R7.** The seeding retraction was fixed by REPAIR and the report still lists it as outstanding,
   while dropping REPAIR's actual finding (per-target values moved up to 0.128).
3. **R5.** The Roget/enumeration result was accepted, but the report then quotes the pre-matching
   numbers, which flip the Legacy-vs-AMBER direction relative to the matched comparison it also
   cites.

**Net judgement.** The corrections are real and several of them cost the sprint its best claims —
this is not a report performing candour over an empty result. But the display of self-correction is
currently *more* polished than the practice, and R1 is the proof. Fix R1 and the §7 material becomes
the strongest section in the document; leave it and §7 reads as a rhetorical device, which is
precisely the failure the reader is primed for.

---

# Q4 — What the hostile expert asks that the report cannot answer

### 1. "On nine targets with three seeds, what is your minimum detectable effect, and what is your power to detect the 0.05 Å you are calling a null for Legacy?"

**The answer exists and was merely omitted.** It is computable in three lines from VERIFY's own
per-target table: SE = 0.051 Å, MDE at 80% power ≈ **0.16 Å**, power at a true 0.05 Å effect ≈ **9%**
(`s16/review_checks.py`). Every §5.4 null — Legacy, the interaction, α — is an absence of evidence at
effect sizes smaller than roughly three times the effects being discussed, and none of the nulls is
labelled as such. The report says only "nine targets is inside the range where this programme has had
conclusions reverse." That is the weakest possible form of the correct answer. **Fix: R11.**

### 2. "Your own QPHASE workstream calls uniform random 'a WEAK CONTROL, proves nothing' and reports that the VQE loses to greedy 1-opt at 10 of 10 objective-quality rungs and is matched by annealing at a quarter of its budget. Why does your report's quantum section lead with uniform random?"

**The answer exists and was omitted.** `qphase_FINDINGS.md` was on disk thirteen minutes before the
report was written. There is no scientific reason for the omission — it is a sequencing accident —
but the omitted material is *stronger* than what is printed, on twice the targets, against controls
that cannot be dismissed. It also contains the one thing that would blunt the objection: the
directional disagreement at `signal = 0` on the untrained-circuit control (QPHASE −0.309, INTEGRATE
+0.112), which has a readout-dependence explanation the report should give rather than have extracted
from it. **Fix: R4.**

### 3. "Your step is `t* = c‖r‖`. `‖r‖` is the RMSD to the native. What supplies it at inference, and what is the 2.491 Å ceiling when it is estimated rather than given?"

**This answer does not exist at all.** No arm in Sprint 16 estimates the step magnitude native-free
from the coordinate average, no arm re-computes the β-strand ceiling under a leave-fold-out constant
step, and `csteer`'s ladder — the only place a step could have been fitted — has no negative rung and
so could never have expressed the sign in the first place (VERIFY §1.4: 45–67% of targets need one).
The report's §9 states the open direction as one bit *because* the magnitude question was never asked.
This is the only one of the three that requires new work: a single leave-fold-out constant-step arm
from `avg` along the constant β-strand direction, priced against the same random control, at n = 126.
Until it is run, **2.491 Å is a ceiling on a two-parameter oracle, not on a one-bit classifier**, and
§9 should say so. **Fix: R3.**

*(Runner-up, answerable from the artefacts: "On the identical binding-masked population Legacy's
certified optimum is 0.109 Å better than AMBER's. Why does §5.7 quote Legacy's full-enumeration
number against AMBER's masked one, and why is §0's headline 'AMBER ranks better than Legacy'?" —
see R5.)*

---

# VERDICT

**Not yet fit to hand to a hostile expert. Two of the twelve findings are disqualifying as written;
the rest are wording.**

- **R1 is disqualifying** because it is the report's own catalogued failure mode, committed in the
  same document, one section after it is catalogued, with an explicit false claim of survival
  attached. An expert who reads §7 and then re-reads §5.5 will discount every methodological claim
  in the document, and will be right to.
- **R2 is disqualifying** because the executive summary's first paragraph contradicts its own fourth
  bullet and its own closing paragraph, and the contradiction is on the single question a reader most
  wants answered ("did anything work?"). The answer is "one thing, barely, at n = 126" and the report
  should say that in the first paragraph rather than deny it there and concede it twice later.
- **R3 and R4** are not disqualifying for correctness but are the two places the report will be
  caught doing less than it claims: an accepted correction whose consequence was dropped, and a
  sibling workstream's stronger result left on the floor.
- **R5–R8** are real defects of the kind an expert finds by spot-checking, and R8's `c = 0.615` is
  the kind of number a hostile reader checks first *because* the sentence around it claims to turn
  the goal into arithmetic.
- **R9 is the honest counterweight**: the report under-sells the one result that distinguishes
  hypotheses rather than refuting one, and reads as more purely negative than the evidence requires.

**After R1–R8 and the R9 promotion are applied, this is a strong document** — it is checkable, its
controls are present, its ORACLE quantities are labelled almost everywhere, its retractions are
preserved with the evidence that forced them, and its §7 is genuinely transferable. The scientific
content is a clean refutation with a mechanism, and it survives the reading. The problem is not the
science; it is that a report about reading quantities correctly must not misread four of its own.

---

## Modules written by this workstream

`s16/review_checks.py` — the three checks above (closed-form `c` requirement; n = 9 power from
VERIFY's per-target table; the isotropic-null value). No heavy compute; `stable_rng` from
`s15/seed.py`; the sealed 60-target benchmark was not read, listed or referenced.
