# SPRINT 16 — DECISION AND RETRACTION LEDGER

Written as the sprint runs, not reconstructed afterwards. Every entry carries the date, who made
the call, what the evidence was **at the time**, and — where it later changed — what overturned it.
Sprint 15 closed with eleven coordinator retractions, all preserved; the same standard applies here.

The ledger is the place a reader checks whether a number in the final report was **pre-registered
or discovered**. If a claim is not traceable to an entry here, it is a discovered claim and must be
labelled as one.

---

## L1 — The flagship was chosen before any Sprint 16 measurement (2026-09-06, coordinator)

**Decision.** The first major investigation is whether the two separately estimable native-free
quantities from Sprint 15 — the weak direction `u = θ_pool − θ_fit` (|cos(u, e)| ≈ 0.36–0.39) and
the Jacobian's loud/quiet spectral split — combine into a steering mechanism.

**Evidence at the time.** Sprint 15 established both quantities independently and neither had been
tested against the other. It is the largest demonstrated unexplained ceiling in the programme.

**Alternatives not taken.** A shallow sweep of architectures was explicitly rejected in favour of a
deep investigation of one; alternatives are permitted to displace it only on evidence of materially
better expected payoff, not on novelty.

**Status.** Executed as `s16/steer.py`.

---

## L2 — A falsifiable worry was stated in the module BEFORE the experiment (2026-09-06, coordinator)

`s16/steer.py`'s own docstring records, ahead of any run, that quiet directions are *defined* as
barely moving the superposed CA trace, so motion along them should barely change RMSD **by
construction** — which means the ORACLE rotation that motivated the sprint may be a counterfactual
requiring `θ_native` rather than a reachable move.

This matters for the final report: the worry is **pre-registered**, not a post hoc rationalisation
of a null.

---

## L3 — RETRACTED IN ADVANCE: the consensus direction was computed with the wrong mean (2026-09-06, coordinator)

**Defect.** The `consensus` arm averaged the multi-start ensemble's torsions with an arithmetic
mean. Torsions are circular; an arithmetic mean of angles is wrong near the branch cut and places
the consensus at a torsion no member holds.

**Caught** before the n = 126 run, in review of the n = 6 smoke, not by a test.

**Fix.** Circular mean, `atan2(mean sin, mean cos)`.

**Consequence for the record.** The n = 6 smoke's `consensus` and `c_loud*` numbers are void. They
are not quoted anywhere. The n = 8 re-smoke used the corrected mean.

---

## L4 — The n = 8 smoke was hostile to the flagship, and this was written down before n = 126 (2026-09-06, coordinator)

Three results in the n = 8 re-smoke contradict the hypothesis:

- `loud4` (+0.272) did not beat `full` (+0.382); `loud10` (+0.495) was worse than both. Projecting
  the native-free direction onto the loud subspace did **not** raise its usable fraction.
- `ORACLE_loud4` was **+0.320** — stepping the *true* error's loud component made RMSD worse —
  while `ORACLE_true` reached 0.060 Å.
- the inert control moved: `quiet_half` −0.106 [−0.185, −0.031], CI excluding zero, while
  `ORACLE_quiet_half` did nothing (−0.008).

**Action.** `s16/PREREG_steer.md` was written **while the n = 126 run was in flight** and fixes the
decision rule and the follow-up diagnostics before any full-instrument number returned. The
follow-up module `s16/curv.py` was written from that pre-registration and not modified after the
n = 126 numbers were seen.

**The general lesson, stated here because it is likely to recur:** a first-order spectral
decomposition (`rank-4 carries 0.80 of the RMSD`) is **not** a finite-step budget. The measured
curvature ratio at n = 8 — first-order prediction 6.183 Å against an actual 2.725 Å, ratio 2.269 —
says the linearisation is not valid at this error magnitude, and that is the mechanism by which an
ORACLE loud step can lose ground while the decomposition looks favourable.

---

## L5 — Panel B of `s16/curv.py` selects its step in sample, and says so (2026-09-06, coordinator)

The follow-up's step ladder is read at its argmin, which is an in-sample selection and makes every
arm optimistic. It is kept because the diagnostic question is **arm versus arm** (does a random
matched-norm quiet step gain what the native-free quiet step gains?), and the same rule applies to
every arm. The module prints the caveat itself. No panel-B number may be quoted as an
out-of-sample effect; `s16/steer.py`'s leave-fold-out ladder is the only place those live.

---

## L6 — Agents and compute (2026-09-06, coordinator)

Four agents run continuously — AUDIT (`s16/audit_align*.py`), QPHASE (`s16/qphase*.py`), ENERGY
(`s16/energy_*.py`), LIT — with the coordinator holding the flagship. The cap of eight is not
approached. Utilisation was checked at 87% CPU / 71% memory before the flagship was launched
alongside them; the ceiling of 97% is treated as a hard limit, not a target.

---

## L7 — The 60-target benchmark is untouched (2026-09-06, coordinator)

No arm, control, or diagnostic in this sprint reads it. It remains sealed and unspent, as it was
left at the end of Sprint 15. Unlocking requires all five pre-registered conditions, and none has
been met: the architecture is not frozen.

---

## L8 — What "native-free" means in the steering ladder, stated exactly (2026-09-06, coordinator)

A reader will ask this, so it is answered here rather than in a footnote.

**Native-free**, with no qualification: the fit target `dhat` (distogram, leave-fold-out separation
debias), the multi-start set (`s15.distgeo.starts`, seeded by `stable_rng`), the direction
`u = θ_pool − θ_fit` (circular-mean torsions of the top-75 retrieval windows), the consensus
direction (circular mean of the ensemble's own fits), the Jacobian and its loud/quiet split. None
of these reads a native coordinate.

**Fold-honest supervised, not native-free**: the **step size**. It is chosen by minimising mean
RMSD on the four training folds and applied to the held-out fold. Native labels therefore enter a
hyperparameter — legitimately, by the same convention the programme already uses for the
separation debias (`s15.distcal.fit_correction`, itself fitted on native distances leave-fold-out),
but it is a *trained* quantity and must not be described as native-free.

**Consequence.** An arm's reported effect is the effect of *a direction plus a step trained on
labels*. Since `s = 0` is on the ladder, the training folds can always choose "do nothing", so a
positive (worse) test-fold effect — `full` at +0.382 in the n = 8 smoke — is a genuine failure of
generalisation and not an artefact of the selection rule.

**ORACLE, and labelled as such everywhere it appears**: the error vector `e = θ_native − θ_fit`
and every arm whose name begins `ORACLE_`, the loud/quiet decomposition table (panel A), and the
curvature ratio. These are diagnostics and ceilings. None of them may be quoted as a result.

---

## L9 — FLAGSHIP VERDICT: refuted at n = 126, and the pre-registered rules decided it (2026-09-06, coordinator)

`s16/steer.py`, n = 126. Every native-free arm chose "do nothing" on every fold; `loud4`, `loud10`
and `loud13` all returned +0.000 [+0.000, +0.000] against the plain fit. Rule 2 of
`s16/PREREG_steer.md` fires: **hypothesis refuted at the mechanism level.** Rule 4 also fires:
`ORACLE_loud4` bought −0.028 [−0.106, +0.051] against a rank-4 first-order RMSD share of 0.875, so
**no rank-r share may be quoted as a finite-step budget** anywhere in this sprint or in the Sprint 15
provenance sections.

Rule 3 did **not** fire: `quiet_half` came back −0.003 [−0.011, +0.007]. The n = 8 anomaly that drove
the whole follow-up plan was a false positive **whose CI excluded zero at n = 8**. Recorded as the
sixth small-n reversal in the programme and the first with a small-n CI that excluded zero.

**The mechanism, which was not anticipated.** Stepping the *exact* native error in torsion space
realises 5.0% of the achievable gain at half-way and 29.2% at seven-tenths, and is non-monotonic
over the first third. RMSD is strongly convex in torsion displacement toward the truth, so partial
motion toward the truth does not pay partially and no improvement to the direction estimator could
rescue the mechanism. Full account in `s16/STEER_FINDINGS.md`.

`s16/req.py` (the ORACLE requirement curve, written before the verdict) is **not run**: the flagship
answered its question by a shorter route. Sweeping direction quality is pointless once the exact
direction is shown to pay nothing at partial step. The module is kept, unrun, so the reasoning is
inspectable.

---

## L10 — LIT: two of the sprint's ingredients are prior art, and one carries a numerical error (2026-09-06, LIT workstream)

- The **fusion law** `d_avg ≈ √(r² − (s/2)²)` is the **Krogh–Vedelsby (1995) ambiguity
  decomposition** for two members, and is exact in any inner-product space with none of the
  caveats Sprint 15 attached. Sprint 15 also computes `r` as an **arithmetic** mean where the
  identity needs the **quadratic** mean — an error of the same size as the reported prediction
  error. Assigned to the RETRACT workstream.
- The **quiet subspace** is prior art in substance: manipulability ellipsoids, torsion-space NMA
  (iMod/iMODS), protein IK self-motion manifolds, KGS nullspace sampling, and "Proteins Wriggle".
  The native-free half is the **Gauss–Newton assumption**; the mechanism is **Hansen's discrete
  Picard condition**.
- **Native-free estimation of structural error direction is taken twice**: ATOMRefine (2023) and
  DeepAccNet (2021). **The number that should have set the prior: ATOMRefine buys GDT-HA
  69.84 → 70.04.**
- **Only the COMBINATION was unclaimed** — and the combination is what L9 just refuted. The
  sprint's contribution is therefore a **clean refutation with a mechanism**, not a method.
- The Shao & Zhu (2018) AMBER threat is **reinstated**: Sprint 15 downgraded it on a condition that
  is false. Assigned to RETRACT.

---

## L11 — ARCHITECTURAL PIVOT, on evidence rather than novelty (2026-09-06, coordinator)

The flagship's mechanism names its own successor. In **coordinate** space the payoff from a
direction is exact algebra rather than an empirical hope:

    ||r − t v||² = ||r||² − 2t(v·r) + t²   →   t* = c||r||,   ||r_new|| = ||r||√(1 − c²)

so any direction with c > 0 pays, the optimal step is closed-form, and the whole question reduces
to one measurable native-free quantity, **c**. It also converts the sprint's target into
arithmetic: from the incumbent's 3.204 Å, **2.5 Å needs c = 0.625** and **2.0 Å needs c = 0.781** *(0.615 as first written; the
figure was hardcoded in `s16/csteer.py` for a 3.170 Å base and is corrected here)*.
Sprint 15's fusion law is the fixed-half-step special case of the same algebra.

This is a pivot the user's brief permits — an alternative displaces the preferred architecture only
on evidence of materially better expected payoff, and the evidence here is that the preferred
architecture's step does not pay at all in the space it was posed in.

**The falsifier was written before the run**, in the module's docstring: if every native-free
channel's `c` is at or below the torsion-space 0.318, the mechanism buys at most ≈ 0.19 Å from the
fit and cannot reach the incumbent. Executed as `s16/csteer.py`.

**The prior is low and is stated as low.** LIT's ATOMRefine margin (GDT-HA +0.20) is the right
reference class for native-free coordinate-shift estimation, not the ORACLE ceiling.

---

## L12 — Agent roster after LIT completed (2026-09-06, coordinator)

LIT finished; RETRACT was launched into the freed slot with the two correction tasks LIT generated.
Roster remains four — AUDIT, QPHASE, ENERGY, RETRACT — with the coordinator holding `csteer`. Load
was sampled before launching, not assumed.

---

## L13 — PIVOT VERDICT: coordinate-space steering refuted too, by a stronger falsifier (2026-09-06, coordinator)

`s16/csteer.py`, n = 126. The pre-registered falsifier was "every native-free channel's cosine is at
or below the torsion-space 0.318". What fired was stronger: **from the best available start (the
top-75 coordinate average, 3.048 Å) every native-free direction has cosine with the true residual
indistinguishable from zero** — toward the fit 0.012 [−0.052, +0.075], toward the projection −0.068,
toward the medoid −0.010, toward the ensemble mean −0.004, against a random control at 0.005. The
requirement is 0.625 for 2.5 Å.

**A trap was caught inside this module and is recorded rather than buried.** Two arms showed large
effects (−0.599 and −0.164 Å) and both chose **step 1.0 on every fold** — they walk the entire way to
the pool coordinate average and stop. They are endpoint substitutions, not steers, and the structure
they reach is one the incumbent pipeline already computes. **General rule adopted: a leave-fold-out
step selector that lands on the boundary of its ladder is choosing an endpoint, not a step, and the
arm must be read as a substitution until shown otherwise.**

**What the module did establish, and it is the sprint's cleanest exhibit:** coordinate-space RMSD is
**exactly linear** in motion toward the native — 3.647·(1 − s) to three decimals at every rung,
through Kabsch superposition — while torsion space on the same targets and the same fits is
non-monotonic and reaches 5% of the gain at half-way. Full account in `s16/CSTEER_FINDINGS.md`.

---

## L14 — AUDIT: the flagship's premise was ~5x weaker than the record said (2026-09-06, AUDIT workstream)

Every recorded Sprint 15 number **reproduces exactly** — 0.390, 0.157, +0.499 [+0.326, +0.645],
0.827, 0.735, 1.855, −1.822, 3.677, 1.443. There is no transcription error anywhere in the record.
**What fails is what those numbers were compared against.**

- **The quiet-subspace identification is refuted.** A zero-information constant α-helix Jacobian
  overlaps the native's bottom-half subspace at **0.829**, a random pool window at **0.833**, against
  the emitted structure's **0.827** — and the emitted structure is *significantly worse* at the
  quarter. The 0.499 "null" is exactly k/m = 0.500, i.e. "half the space". Four of the five null
  directions are fixed coordinate axes, so ≥ 0.800 of the null-5 overlap is automatic.
- **The 1.855 Å ceiling is 83% direction-free.** Its own zero-information control — the same error
  norm along a *random* direction in the same quiet half — emits 2.164 Å. Direction is worth
  −0.308 [−0.435, −0.172], 16.9% of the recorded 1.822.
- **The surrogate is not a pool signal.** A zero-information α-helix reference gives |cos| 0.357 and
  tracking +0.441; the pool adds **+0.033 [−0.007, +0.075]**, CI including zero. The correct
  magnitude-matched null is **0.216, not 0.157**.
- **The theory check is half-refuted.** RMSD is blind to quiet motion (ratio 0.0018–0.0024), but the
  **objective is not** (0.018–0.031): a 20° quiet step raises the objective by 53% of its value at
  the optimum. Ratio of ratios 10–17, > 1 on 99.2% of targets.
- **One rule violation found**: `s15/align_jac.py:136` seeds with `hash()`. The value was verified
  correct analytically, so no number changes; assigned to REPAIR to fix.
- Leakage and reproducibility clean: benchmark 0/60 overlap and read by no `s15/` or `s16/` module,
  `seed.py` at all five multi-start sites, five pinned constants exact.

**The standing methodological consequence, adopted for the rest of the sprint:** every arm must
carry **both** a zero-information reference control and a matched random-direction control. Without
those two, this programme's evidence read about **five times stronger than it was**.

AUDIT's one surviving item — the surrogate's alignment with the **loud** half, 0.497 against a
matched null of 0.304, W/L 101/25 — is **closed by the flagship**, not left open: stepping the exact
loud component bought −0.028 [−0.106, +0.051].

---

## L15 — Roster (2026-09-06, coordinator)

AUDIT and LIT complete. Running: QPHASE, ENERGY, RETRACT, REPAIR — four, the standing default, with
the cap of eight never approached. The coordinator holds `s16/integrate.py`, which was queued to
start only when `csteer` released the CPU rather than launched alongside it.

---

## L16 — RETRACT: the fusion law is a theorem, and the 0.162 Å was two of our own errors (2026-09-06, RETRACT workstream)

**Sprint 15's "parameter-free fusion law, validated to 0.162 Å on 126 targets" predicts nothing.**
It is an exact algebraic identity. Reproduced bit-identically (max |Δ| = 0.000e+00), then decomposed:
the residual falls +0.1619 → **+0.0748** on switching to the quadratic mean (paired
−0.0871 [−0.1185, −0.0605] i.i.d., [−0.1285, −0.0556] fold-clustered, **W/L 126/0, 5/5 folds**), and
the remaining +0.0748 decomposes **exactly** into +0.1738 (averaging in the medoid frame rather than
the target's) − 0.0990 (Kabsch-minimised `s`). **In one common frame with the quadratic mean the
residual is 2.65e−15.** The final superposition contributes exactly zero, provably (the
cross-covariance is linear in the moving structure and a sum of symmetric PSD matrices is PSD),
verified on 126/126 and on 3000 random trials.

**So the 0.162 Å was 54% our arithmetic-mean error and 46% our frame-convention mismatch, and 0% the
identity failing.**

**Novelty surviving: none.** LIT's proposed survivor (a) — "it holds through Kabsch on 126 real
targets" — is forced, and Sprint 15 never measured it. LIT's survivor (b) — "`s` is native-free so
it works as an a-priori fusion screen" — is **REFUTED**: `s` alone ranks "does fusion win here" at
**AUC 0.401 against a 0.500 null**, and a supervised leave-fold-out threshold degenerates to "always
fuse" on 5/5 folds. What is left is one measurement (+0.602 per-pair error correlation at 3.138 Å
structural disagreement) and one reportable negative — a consensus signal that is *anti*-informative
about fusion.

**A methodological finding that matters more than the number.** `s15/verify_FINDINGS.md` **passed**
this row: every figure traced to the artefact exactly. **A transcription audit cannot see a formula
error.** Sprint 15 ran three hostile reviewer passes and none of them caught this, because they were
all checking that numbers matched their sources.

**RETRACT also corrected LIT, in Sprint 15's favour.** LIT conflated two sources: Sprint 15
downgraded *Maffucci & Contini* (where the ground is true) and separately flagged the Shao & Zhu
PCCP study as UNVERIFIED with the consequence pre-registered. Sprint 15's risk localisation was
correct. Two natural defences of our AMBER position nonetheless fail and must not be used: their
1E0Q is a 17-residue β-hairpin, one residue above our band; and "free energy is a weaker observable"
is backwards — it is Boltzmann-weighted and therefore *stronger*, so their result is upstream of
ours. **"AMBER anti-ranks peptide natives" moves from a finding to a confirmation of documented
expected behaviour.** What survives is the observable, the statistics, the AMOEBA extension, and —
untouched by this literature — **stereochemical repair**.

---

## L17 — INTEGRATION: AMBER contributes, Legacy does not, and the VQE loses to random (2026-09-06, coordinator)

`s16/integrate.py`, 9 enumerated targets × 3 seeds. Full account in `s16/INTEGRATE_FINDINGS.md`.

- **AMBER as a ranker inside the architecture: −0.054 [−0.088, −0.020] against a matched random drop
  at m = 75**, 101W/61L; −0.051 at m = 20. ~~The only positive accuracy effect in Sprint 16.~~
  **[Both the interval and this description were withdrawn — see L21 (VERIFY) and L23 (ENERGY). The
  arm does not survive the target as the unit; the sprint's one surviving accuracy effect is a
  different measurement, AMBER as a relaxation operator on the shipped ensemble.]**
- **Legacy: +0.005 [−0.034, +0.044] at m = 75** — nothing — and **+0.068 [+0.012, +0.125] at
  m = 5**, i.e. significantly *worse* than random. **AMBER − Legacy = −0.059 [−0.089, −0.030].**
- The composition adds nothing over AMBER alone (interaction −0.004 [−0.042, +0.032]).
- **Every CVaR-VQE arm loses to uniform random sampling at matched budget** (+0.209 to +0.301, CIs
  excluding zero, 5–8W/19–22L). α shows no monotone trend and plain VQE is among the better arms.
  This reconfirms the standing result on a **new terminal operator** (ensemble → coordinate average
  rather than argmin).

**SELF-DECLARED DEFECT, recorded before VERIFY reports.** `s16/integrate.py`'s bootstrap resamples
**rows**, and 162 rows sit on only **9 targets**. The programme's rule is that the TARGET is the
unit of analysis, so **every interval in that file is too narrow**. Point estimates and win/loss
counts stand; the intervals are provisional pending VERIFY's recomputation. The caveat is printed at
the head of `s16/INTEGRATE_FINDINGS.md`, not buried.

---

## L18 — The mechanism behind the quantum result, and a correction to a project law (2026-09-06, coordinator)

`s16/diversity.py`. At m = 75 the six generators' **mean member error spans 0.056 Å** while their
**readouts span 0.301 Å**; the difference is **diversity** (2.393 to 2.626), whose ordering
reproduces the readout ordering. The Krogh–Vedelsby decomposition holds through the operator to
machine precision (residual −0.0000, max 0.0000, all three ensemble sizes), so naming which term the
generators differ in is a **complete causal account**: the VQE does not make worse conformers, it
makes **less various** ones, and the terminal operator's entire gain is the diversity term.

**Panel A is a theorem check, not a discovery, and is labelled as such** — RETRACT established the
identity's exactness independently and at machine epsilon on the same day.

**This corrects a recorded project law.** "The terminal operator consumes the set MEAN, not the set
BEST" does not account for this instrument: the set means are equal to 0.056 Å. The missing term is
diversity, and with it the account is exact. It also names the mechanism behind "concentration is
wrong when discrimination binds" — concentration *is* the destruction of the diversity term, which
is why α reads as a diversity dial rather than a tail parameter.

**Native-free diversity is a weak heuristic, not a signal.** Within-target rank correlation with the
readout is only −0.145 to −0.270. It orders the six generators; it does not order individual runs.
RETRACT independently refuted the related native-free screen at AUC 0.401 against a 0.500 null.

---

## L19 — VERIFY launched against Sprint 16's own claims (2026-09-06, coordinator)

Because AUDIT showed Sprint 15's evidence read ~5x stronger than it was, and RETRACT showed a
transcription audit cannot see a formula error, a VERIFY workstream was launched to attack Sprint
16's own claims — the coordinator's first. It is asked specifically whether the flagship's
convexity result is a property of the error or of interpolating **wrapped angles** off the circle;
whether the coordinate-space linearity exhibit is a tautology; whether the csteer cosines survive a
zero-information α-helix reference; whether the integration CIs survive using the target as the
unit; and whether the diversity panel is a theorem check dressed as a discovery. Roster: QPHASE,
ENERGY, REPAIR, VERIFY.

---

## L20 — The mechanism has no lever: diversity-aware selection is a null at n = 126 (2026-09-06, coordinator)

`s16/divselect.py`, the real 126-target instrument. Deployable leave-fold-out arm **+0.015
[−0.019, +0.049], median +0.000, 61W/65L** against the shipped top-75 selection read through the
identical operator. The pre-registered falsifier fires.

The null is informative because both Krogh–Vedelsby terms were printed per arm: raising the mixing
weight lifted diversity from 1.902 to 2.316 (+22%) and member error from 3.615 to 3.778 (+4.5%), and
the readout moved by at most 0.008 Å. **The terms rise together and cancel** — which is what an
identity with equal and opposite weights requires. This is not "the intervention did not fire"; it
is "the intervention fired and was paid for at par".

**The module's own pre-registered trap fired**: three of five folds chose λ = 1.5, the end of the
ladder. Recorded as the second instance in this sprint of a selector pinning to a boundary
(`s16/csteer.py` was the first), and the flag was written into the module before the run.

**Instrument validation**: the greedy selector at λ = 0 reproduces the shipped top-75 exactly
(3.048 Å, +0.000), so every arm differs only in selection.

**One narrow positive, with its caveat inseparable from it**: within the top-200 by score, selecting
by **pure diversity** matches the score-based selection (3.051 vs 3.048) while matched-random gives
3.090. Because the candidate set is itself the top-200 *by score*, this says only that the score's
selective work is done in the 500 → 200 cut, not in the 200 → 75 one. It is **not** a claim that the
score is worthless, and must never be quoted as one.

---

## L21 — VERIFY: THREE OF THE COORDINATOR'S FIVE CLAIMS FAIL, AND FAIL THE SAME WAY (2026-09-06, VERIFY workstream)

**Every published Sprint 16 number reproduces from its artefact exactly.** As in Sprint 15, that was
never the question. Three of five assigned claims fail in the way RETRACT named earlier the same
day: *a quantity is measured correctly and then read as if it were a different quantity*.

**RETRACTED — `s16/CSTEER_FINDINGS.md` §2, "the direction does not exist."** The module derives its
own law, `‖r_new‖ = ‖r‖√(1 − c²)` at `t* = c‖r‖`, and implements it correctly. That law is
**two-sided and quadratic in c**. The findings then quoted the **arithmetic mean of the signed c**
and concluded the direction was absent. Read as the law reads it, from `avg`, against the module's
own random control (0.128, which is exactly the isotropic null √(2/(π·3n)) at n = 13):

| direction | signed c (as published) | \|c\| | \|c\| − random [95% CI] | W/L |
|---|---|---|---|---|
| toward `fit` | +0.012 | **0.305** | **+0.177 [+0.142, +0.212]** | **100/26** |
| toward `ens` | −0.004 | 0.301 | +0.173 [+0.135, +0.210] | 98/28 |
| **constant ideal β-strand (ZERO-INFORMATION)** | −0.012 | **0.398** | **+0.262 [+0.216, +0.310]** | — |
| pure isotropic breathing (ZERO-INFORMATION) | −0.067 | 0.296 | +0.167 [+0.123, +0.214] | — |

**A constant ideal β-strand carries more usable direction from the best start than any channel
Sprint 16 built.** The residual `avg → native` is dominated by a **compactness scalar**: the
magnitude is there and the **sign flips per target** (58% need contraction, 42% expansion). And
`csteer.STEPS` has **no negative rung**, so on 45–67% of targets *no arm on that ladder could have
helped whatever the direction carried* — the observed null is partly a property of the ladder.

**The operational conclusion survives; the stated reason does not.** Nothing is deployable and
2.5 Å is not reached, because both the sign of `c` and `‖r‖` are native-derived. But the correct
statement is *"the magnitude is 2–3× the isotropic null and the per-target sign — one bit — is
missing"*, not *"the direction does not exist"*.

**RETRACTED — `s16/STEER_FINDINGS.md` §3's causal sentence.** The numbers are right (7.66%, median
16.08%, 77/49 all reproduce) and the per-target aggregation was the correct choice. **The geodesic
objection raised in the VERIFY brief is refuted**: `A.wrap` *is* the exact per-angle geodesic
(agrees with an independent great-circle construction to 1.8e−15 rad). But the sentence "the fit's
torsion errors are *compensating*" had **no control**. VERIFY supplied one: geodesic interpolation
between two arbitrary retrieval windows, native-free, matched on ‖Δθ‖, gives **−0.024 at s = 0.5
against the flagship's +0.014** — *the control bulges more*. The whole effect tracks ‖Δθ‖
(Spearman −0.571), and **below ≈4 rad both are linear**. The fit's error is **90° per-angle RMS,
87% of a uniformly random torsion assignment.** So the convexity is **a magnitude regime, not a
property of torsion space**, and the correct reading is that *the torsion-space fit is very nearly
a random torsion assignment*.

**DEMOTED — `s16/CSTEER_FINDINGS.md` §1, the "cleanest exhibit".** It is a tautology: max deviation
from `(1 − s)·rmsd` on **400 random point-cloud pairs with no protein in them is 4.6e−07 Å**. Sum of
two symmetric PSD cross-covariances is PSD, so Kabsch never leaves the identity along the segment —
the same argument RETRACT used to retract the fusion law, and `csteer.py`'s own docstring already
writes the algebra out. **An exhibit derivable on random point clouds is not an exhibit.**

**DOWNGRADED — the AMBER filter.** With the TARGET as the unit, as the programme's rule requires:
**−0.054 [−0.148, +0.040]**, CI **2.79× wider**, 4 targets worse / 5 better, **dropping any of 8 of
9 targets puts the CI over zero**, and **seed 2 alone gives −0.0003**. No per-generator CI excludes
zero. This is the recomputation of the defect the coordinator self-declared in L17, and the point
estimate **does not carry** the description "the only positive accuracy effect in Sprint 16". The
tie trap was censused and **does not apply** — zero ties in Legacy or AMBER.

**CORRECTED — the diversity attribution.** Panel A is the parallel-axis theorem at machine epsilon
(max relative 5.5e−15) and was correctly labelled. But the table printed the **free-superposition**
member error (spread 0.056 Å) while the identity consumes the **common-frame** term (spread
**0.122 Å**). Exact per-target attribution: **cvar0.25 is 52% member error, cvar1.0 51%, cvar0.5
42%.** "They differ in diversity, not conformer quality" — and L18's correction of the set-mean
law — are **wrong by about a factor of two**. It is roughly half and half.

**SURVIVED EVERYTHING**: the coordinate-over-torsion readout (**+0.446 [+0.193, +0.706], 8/9
targets** at target level — the only accuracy claim in `INTEGRATE_FINDINGS` that holds); the VQE
losing to uniform sampling (those contrasts are paired within target, so the row bootstrap was
nearly right); the endpoint-substitution trap, correctly caught and correctly barred; the
projection's 0.165 Å cost; no tie leak.

**STANDING RULE ADOPTED, and it is stronger than the one AUDIT gave us** — a zero-information
control alone would not have caught three of these:

> **Before quoting a statistic as evidence for a law, derive which functional of it the law
> consumes, and quote that functional.** Where a law is quadratic, sign-invariant, or
> frame-dependent, the arithmetic mean of the raw per-target quantity is not it. A mechanism
> sentence needs its own control. An exhibit derivable on random point clouds is not an exhibit.

---

## L22 — VERIFY's open question is CLOSED by an existing Sprint 15 artefact (2026-09-06, coordinator)

VERIFY left one question open: the compactness bit has an ORACLE two-sided ceiling of **2.491 Å**,
below the sprint's primary 2.5 Å goal, and *"whether that bit is obtainable native-free is not
tested here"*. Before building `s16/signbit.py` to test it, the coordinator checked the record —
and **Sprint 15 had already run exactly this experiment** (`s15/expand.py`, `s15/results/expand.json`,
n = 126), testing four native-free scale references with the closed-form optimal isotropic scale
`s* = Σ w d dhat / Σ w d²`.

Re-read at the level of the **sign**, which is what VERIFY's one bit actually asks for:

| native-free reference | sign accuracy vs the RMSD-optimal scale | correlation with it |
|---|---|---|
| distogram `s_pred` | **0.357** [0.278, 0.444] | +0.021 |
| separation-debiased `s_deb` | 0.389 [0.302, 0.476] | +0.011 |
| pool-matched `s_pool` | 0.413 [0.325, 0.500] | −0.037 |
| both `s_both` | 0.389 [0.310, 0.476] | −0.019 |
| **zero-information constant-sign baseline** | **0.516** | — |

**All four native-free references predict the sign WORSE than always guessing the majority sign.**
Sprint 15 also measured the cost of using them: +0.080 to +0.123 Å, every CI excluding zero.

**A sharper finding falls out of the same table.** Against `s_true` — the scale that matches the
*predicted distances* — the references reach 0.54–0.61, i.e. above chance. Against the
**RMSD-optimal** scale they reach 0.36–0.41, below a constant guess. **Matching the predicted
distance scale is not the same as minimising RMSD, and the two disagree in sign.** That is why the
references look informative on their own axis and are anti-informative on the axis that matters.

**Decision: `s16/signbit.py` was NOT built and NOT run.** The experiment exists, at n = 126, and is
negative. Recorded here so that the decision not to spend the compute is visible and reversible
rather than silent. What remains genuinely untested is whether some *other* native-free estimator —
one not of the distance-matching family — could supply the bit; that is new science and is left as
the sprint's stated open question, not claimed as tested.


---

## L23 — ENERGY: the causal ablation, and a third inherited reporting defect (2026-09-06, ENERGY workstream)

Full account in `s16/energy_FINDINGS.md` (886 lines); claims registered as G1-G14 in `s16/CLAIMS.md`.

**The answer the sprint was commissioned to give, at n = 126 on the shipped ensemble.**
**AMBER contributes 0.023 A of accuracy in aggregate and stereochemical repair; Legacy contributes
detection and nothing to accuracy; they do not interact.**

- AMBER `amb - ctrl` **-0.0234 [-0.0375, -0.0099]**, median -0.0074, 66W/57L, n = 123 gated.
- Legacy `leg - ctrl` +0.0243 [-0.0120, +0.0615]; **against its matched random filter +0.0118
  [-0.0242, +0.0492]**, where the random filter alone costs +0.0124. Legacy's whole effect is the
  cost of dropping 19 of 75 windows at random.
- Interaction +0.0073 [-0.0073, +0.0218].

**The mechanism of Legacy's null is the set-mean trap, caught in the act**: the filter improves the
ORACLE set MEAN (3.551 -> 3.531) while destroying the set BEST (2.306 -> 2.416) and contracting
diversity 6%. That is the identity of L18/L20 acting against the pipeline.

**The workstream's own pre-declared frame null FAILS and was not loosened**: +0.0146 A mean, max
0.140 A on 11 converged targets. Consequence stated exactly: **no per-target AMBER RMSD claim is
resolvable**; the aggregate still clears the 123-target floor by ~10 SE, and the null-calibrated
concentration check returns PASS at the 58.5th percentile.

**Legacy's earned role is DETECTION**: AUROC **0.93-0.99** on `torsion` (not tautological), 0.98-0.99
on `steric` (partly tautological). Its only intervals excluding zero anywhere in the sprint are
+0.032 Ramachandran and +0.069 A minimum separation over a matched random control.

**Two escapes closed.**
- **Fitting Legacy's weights makes it worse**: leave-fold-out fitting moves the certified optimum
  **+0.7785 A [+0.0404, +1.4975]** the wrong way (8W/11L) while gaining substantially in sample. The
  workstream's own hypothesis is refuted and the "it was never trained" escape is gone.
- **The eleven-term total is worse than its own best component on all ten axes**; MJ is at chance,
  `compactness` is an anti-detector, two components are identically zero on all 126 inputs, and
  removing MJ barely moves the argmin defect (+0.082 -> +0.090) -- the defect is distributed.

**A zero-information constant alpha-helix beats Legacy's certified optimum by 0.469 A and IS that
optimum on 8 of 19 targets.** The distance objective beats random by -0.926 [-1.440, -0.431] and
beats AMBER by -0.815 [-1.357, -0.310] on the identical masked set.

**THIRD INHERITED REPORTING DEFECT, found by the workstream it damages.** `core.project.lam_path`
returns *unwrapped* torsions (values to +1022.72 deg) and `rama_ok` compares them against fixed
windows. Corrected, the baseline arm moves 0.4658 -> 0.7029 (+0.2371 [+0.1946, +0.2796], 0W/85L, 85
of 126 targets affected). **AMBER's Ramachandran advantage is +0.171, not +0.408 -- the record
overstates it 2.39x.** Six shared-document lines need the fix. This is the third case in two days of
a number that reproduces exactly and means something other than what was written.

**Two new defects in the repair path, neither previously recorded**: the ideal-geometry projection
**degrades** Ramachandran quality (0.836 -> 0.703) and minimum separation; and AMBER **introduces cis
peptide bonds on 44 of 126 targets**. Cost at the right number: 13.88 s minimisation vs 5.01 s
projection, ~1500x an AMBER single point.

**Free result**: all 126 control-arm targets reproduce `s15/results/phys_repl_canon0.json` at
**max |delta| = 0.000e+00 A**, across sprints, interpreters, a change to `core/amber.py`, and four
worker processes.

---

## L24 — REVIEWER launched against the report itself (2026-09-06, coordinator)

Sprint 15 closed with three hostile reviewer reports. Sprint 16 has had AUDIT against its
inheritance and VERIFY against its results; REVIEWER is launched against the **document**, and is
asked specifically whether the report's foregrounding of its own error-correction is load-bearing or
decorative, whether anything is UNDERclaimed as well as overclaimed, and whether the executive
summary's "no surviving positive accuracy effect" is consistent with ENERGY's -0.0234 [-0.0375,
-0.0099] at n = 126. Roster: QPHASE, REPAIR, REVIEWER.

---

## L25 — QPHASE: there is no phase boundary, and the axis is not rho (2026-09-06, QPHASE workstream)

19 fully enumerated targets, 1.28e7 exactly labelled structures, TARGET as the unit. Claims I1-I11
in `s16/CLAIMS.md`; full account in `s16/qphase_FINDINGS.md`.

**The sprint's decisive quantum question was posed as "at what objective quality does CVaR-VQE begin
to beat its classical controls" and the answer is that the premise is wrong.** Two independent
quality knobs give different rho* for the same control (+0.355 blend, +0.140 noise); a cell-level fit
of the paired difference on rho has **R2 = 0.014** with implied crossings at -2.74 and -1.10, outside
the possible range of a correlation. The disagreement between the parameterisations IS the finding,
and the brief pre-registered that reading.

**What governs it instead**, agreed by both families: `RMSD(argmin E) - RMSD(control)`, Spearman
+0.636, R2 0.427 pooled and **0.859** structured, +0.867 across all 95 real (objective x target)
cells against **-0.038** for rho.

**Against real classical optimisers at matched budget the VQE simply loses.** Significantly worse
than **greedy 1-opt at 10 of 10 rungs**; never significantly better than annealing; **indistinguishable
from annealing at a QUARTER of its budget at 9 of 10 rungs**. Greedy reaches the certified global
optimum in **68-100%** of cells against the VQE's **0-32%**, and the VQE lands 90% of the way there
(`diff = +0.051 + 0.904 D`). **The VQE-minus-classical gap is FLAT across the whole quality ladder**
(+0.098 +/- 0.037 vs greedy) and never crosses zero: improving the objective moves both arms together.

**This retires the control the coordinator's own integration section led with.** QPHASE labels
uniform random a "WEAK CONTROL, proves nothing" and ran the strong ones. Report section 5.5 has been
qualified accordingly.

**A directional tension, stated rather than buried**: on the real objective QPHASE's argmin readout
gives vqe - untrained = **-0.309** (VQE better, null) where `s16/integrate.py`'s 75-member
coordinate-average readout gives **+0.112** (VQE worse). The readouts consume different quantities,
which is consistent with L18/L20, but **the sign of the untrained-circuit contrast is
readout-dependent** and the report now says so.

**Two escapes closed.**
- **alpha is a temperature, and a classically reproducible one.** A diversity-matched classical
  thermostat reproduces alpha's whole ensemble effect (Pearson +0.93/+0.98) and moves the readout
  **1.6-1.7x further**. The recorded 150x diversity range is a property of the optimisation BUDGET,
  not of alpha (1.35x here). This **supersedes** the programme's "alpha is a diversity dial" reading:
  the effect is real, the mechanism is thermal, and a classical thermostat does it better.
- **The weak-optimiser hypothesis is REFUTED three ways.** The CVaR gradient-baseline defect is null
  at the target unit and significantly *worse* on the ensemble readout; its recorded "de-facto
  step-size reduction" mechanism is **impossible under Adam** (a 2.6x smaller gradient gives a 7%
  smaller step); and cutting the learning rate 16x costs +0.156 A (argmin) and +0.676 A (ensemble).
  **More optimisation is better here, not worse.**

**What the VQE does, characterised on four axes**: it **samples** (expected energy percentile
0.50 -> 0.15-0.25) and **represents** (3.6x near-native mass enrichment; mode 3.792 -> 3.014 A); it
does **not optimise to completion** (0.5% of mass on the certified optimum) and does **not explore**
(3,551-4,020 distinct configurations against annealing's 7,265 at the same budget).

**And one number for the next sprint**: the deployed objective's **certified argmin is 2.661 A**
against the incumbent's 3.204 A and a space best of 1.030 A, and a classical 1-opt reaches it exactly
in 100% of cells. **The 3.204 -> 2.661 gap is not a search problem.**

---

## L26 — REVIEWER: the report broke its own rule, inside the section that states it (2026-09-06, REVIEWER workstream)

`s16/review_FINDINGS.md`. Twelve findings; the two disqualifying ones are recorded here because they
are the coordinator's.

**R1, and it is the worst thing in the sprint.** The report's section 5.5 quoted **row-level**
bootstrap intervals and row-level win/loss counts for the quantum table and then asserted *"these
contrasts are paired within target, so they survived VERIFY's recomputation"*. They did not, in one
case: `cvar1.0` is +0.209 **[-0.002, +0.410]** at target level and had been printed as
[+0.046, +0.357]. **This is the fourth of the five failure modes the report's own section 7 lists,
committed two pages after it is catalogued, in the document that catalogues it.** Fixed; and section
7 now records that it happened, because a rule its own author breaks inside the same document is
evidence about how weak the rule is without a separate reader.

**R2.** "No surviving positive accuracy effect at all" was carried out of VERIFY's 9-target context
into the report's opening line, contradicting two later statements in the same section. VERIFY's own
limitations say *"the AMBER question is open, not answered negatively"*. Fixed; the later statements
were right.

**Also fixed**: R3 (the "one bit" framing dropped `t* = c||r||`'s magnitude, which nothing in the
sprint estimates native-free -- this changed section 9's forward plan from one question to two);
R4 (QPHASE had reported and was absent; now section 5.8); R5 (the enumerated-space bullet compared
three numbers across two populations and called them matched -- **on the matched set the
Legacy/AMBER direction reverses**, and every printed CI excluded zero while every omitted one
contained it); R6 (four unit-mixing defects in one table); R7 (the seeding limitation was stale,
incomplete and kept only the reassuring half); R8 (**c = 0.615 for 2.5 A is wrong; it is 0.625** --
hardcoded in `s16/csteer.py` for a 3.170 A base, now corrected at source); R10 (two stale row-level
intervals in `CLAIMS.md`); R11 (no power statement anywhere -- on the 9-target instrument SE is
0.051 A, MDE at 80% is **~0.16 A**, and power at 0.05 A is **~9%**, so every null there below 0.16 A
is uninformative rather than negative).

**R9, the opposite failure.** The **distance-restraint objective's certified argmin, -0.926 A
[-1.440, -0.431], 15W/4L, native at the 29.1st percentile on 19 COMPLETE enumerations** -- the only
interval in the sprint excluding zero on an argmin readout, and the control that turns "neither
energy ranks" into a statement about **contact-style potentials** rather than about short peptides --
appeared only as a subordinate clause. Promoted to the executive summary and to section 6.

**REVIEWER's verdict on the corrections themselves**: *mostly load-bearing* -- the |c| retraction
created section 9, the diversity re-attribution took back a win over the existing record, and the
`lam_path` fix shrank ENERGY's own headline 2.39x -- with **one decorative in effect**, section 7's
row-vs-target rule, confessed and then broken four times.

---

## L27 — REPAIR: the last accuracy claim in the programme is removed (2026-09-06, REPAIR workstream)

`s16/repair_FINDINGS.md`; claims J1–J11 in `s16/CLAIMS.md`. The pipeline reproduces PHYS's arm
**bit-identically** in a fresh interpreter on a different day (max |Δ| = 0.000e+00 Å on both arms,
final energies to 0.000e+00 kcal/mol), so every number below sits on the Sprint 15 instrument.

**BOTH REPAIR OPERATORS ARE WORSE THAN DOING NOTHING.** Raw all-atom coordinate average **3.0498 Å**;
ideal-geometry projection **3.2052** (+0.1554 [+0.1217, +0.1911], **27W/99L**); AMBER k = 30
**3.1831** (+0.1333 [+0.1031, +0.1645], **24W/99L**, gated n = 123).

> **The programme's −0.022 Å is the gap between two ways of being 0.13–0.16 Å worse than the
> operator that does nothing at all**, and each loses to it on 99 of 123 targets.

**A matched-magnitude RANDOM displacement is at least as accurate as AMBER's** — 3.1606 against
3.1831 (AMBER − random **+0.0225 [−0.0138, +0.0571]**, 43W/83L) — so neither operator's *direction*
is worth anything over a random direction of its own size.

**AMBER's displacement points slightly AWAY from the truth, significantly so relative to random.**
ORACLE cos(v, r) = **−0.0521 [−0.0924, −0.0129]**, median −0.073, positive on only **37.3%** of
targets; paired against the matched random control **−0.0491 [−0.0996, −0.0009]**, fold-clustered
[−0.0876, −0.0124]. **This closes the coordinator's question directly: the physics does not supply a
new native-free direction** — it is the only channel in the programme that had never been tested for
one, and it is on the wrong side of zero.

**And it moves where RMSD cannot see.** The true residual is **60.2% loud**, 21.4% quiet, 18.4%
non-torsional. AMBER's displacement is **55.1% non-torsional (1.5× enriched) and 23.5% loud
(depleted)**, over-occupying even the **exactly-null** torsion directions (0.270 against a 0.199
null). The mechanism of the −0.022 Å is therefore a **magnitude** difference and not a direction:
orthogonal-move cost **+0.1177 Å (AMBER) against +0.1515 Å (projection)**, and the algebra closes to
**0.0044 Å**. **AMBER wins by moving 0.09 Å less, not by moving better.**

**There is no (k, depth) at which the accuracy effect survives on the frame-reproducible subset**:
k30-full −0.0095 [−0.0242, +0.0056] 39W/48L; d400 −0.0061; d100 +0.0117; d25 +0.0109; k0
**+0.1602 [+0.0502, +0.2724]**. The pre-registered n = 30 shape check at k ∈ {5, 15, 60} is
**monotone toward "do not move"** and was used to select nothing.

**The frame null FAILS at every converged setting, on its maximum.** Gated mean at k = 30 is
−0.00095 — which is where the ≈0.004 Å floor comes from — but the **maximum is 0.1373 Å, six times
the claimed effect**, and 122 of 122 converged targets move by more than 1e−6 Å. **A frame null must
be reported with its maximum, not only its mean.** The only passing arm (depth 25) loses 76 of 126
targets to the convergence gate.

**The VALIDITY claim survives but only as a conjunction.** rama 0.466 → 0.874 (116W/2L) and clashes
1.397 → 0.000 (63W/0L) reproduce Sprint 15 exactly — **but a zero-information constant α-helix scores
rama 1.000 and zero clashes by construction and beats AMBER 0W/65L**. A validity statistic a
zero-information reference maximises is not evidence about a force field. The defensible claim is:
*reaches 0.874 and zero clashes **while moving only 0.725 Å, 26% of the residual**, where the helix
reaches the same validity by moving 0.883 of it — i.e. by discarding the structure.*

**k is an ACCURACY knob, not a VALIDITY knob.** Unrestrained AMBER reaches the same validity
(rama 0.868, 0 clashes, geometry 0.0170) while costing **+0.2742 Å [+0.1785, +0.3823]**. The force
field supplies the stereochemistry; **the restraint supplies nothing but the decision not to move.**

**Two of the workstream's own priors were refuted and are preserved**: that a shallower minimisation
would be the frame-reproducible sweet spot with intact validity (depth 25 passes the null but fails
the gate on 76/126 and carries 1.7 clashes), and that AMBER's displacement would be uninformative but
*neutral* (it is significantly worse than random).

**What this leaves for the architecture.** The mandated pillar is satisfiable and should be kept —
as a **terminal validity operator with a declared accuracy price**, never as a refinement step. The
price is **+0.133 Å against not running it**, for **12.14 s** per target (≈2.7× the projection,
≈10³× a single point), and the projection it would replace is strictly worse on every axis except
wall clock. **AMBER's repair role is real and its accuracy role is not, and the two must never again
be quoted as one number.**

`s15/align_jac.py:136`'s `hash()` violation is **fixed**: aggregate unchanged (0.4961 → 0.5005
against an analytic 0.5000), but per-target values had been swinging by up to **0.128** with the
interpreter salt.

---

## L28 — SPRINT CLOSE (2026-09-06, coordinator)

**Sprint 16 ends with no surviving accuracy effect on any instrument.** The flagship was refuted, its
replacement was refuted, three of the coordinator's published explanations were refuted by an
internal verification pass, and the last candidate effect — AMBER relaxation — was removed by the
final workstream. The incumbent stands unchanged at **3.204 Å**. The **60-target benchmark is sealed
and unspent**; the first unlock condition fails outright, because the architecture is not frozen —
nothing earned a place in it.

**Nine workstreams reported**: AUDIT, LIT, RETRACT, ENERGY, QPHASE, REPAIR, VERIFY, REVIEWER, and the
coordinator's own flagship line (`steer`, `csteer`, `integrate`, `diversity`, `divselect`). Four
agents ran continuously throughout; the cap of eight was never approached; utilisation was sampled
before every launch and the 97% ceiling was treated as a limit rather than a target.

**Deliverables**: `s16/FINAL_REPORT.md` (the report), `s16/CLAIMS.md` (one row per claim, with
status, provenance and whether it is native-free, trained or ORACLE), this ledger, and nine
workstream findings files with their artefacts in `s16/results/`.

**The single most transferable output is not a result.** Across two sprints and four independent
audits, **every published number reproduced exactly from its artefact**, and roughly a third of the
conclusions built on those numbers were still wrong — because a quantity was measured correctly and
then read as if it were a different quantity. The rules adopted in consequence are in the report's
§7, and the strongest of them is VERIFY's: *before quoting a statistic as evidence for a law, derive
which functional of it the law consumes, and quote that functional.*


---

## L29 — The programme's last standing positive quantum result is retired (2026-09-06, coordinator)

Added after a stale Sprint 15 workstream reported in and restated *"the one positive VQE result --
the unranked ensemble, +0.36 to 0.57 A -- is flagged low-power and needs replication"* as still
open. **It is not open; Sprint 16 closed it, and the report had not said so.**

`s16/qphase_FINDINGS.md` section 3, on the full enumerated register: **both classical searches beat
the VQE on the ensemble readout as well, greedy at every rung, +0.27 to +0.62 A.** Sprint 15's
exception -- *the one place a VQE wins is an ensemble consumed without a ranker* -- **does not
survive** once annealing and greedy are in the table rather than only the untrained circuit.

Recorded as claim I12. The report's section 5.8 now states it; it was an omission, caught because a
sibling workstream restated the superseded claim as current.

**The consequence is worth stating plainly**: with I12 retired, **the programme has no surviving
positive quantum result on any readout** -- argmin, ensemble, or coordinate average -- against any
control stronger than an untrained circuit.
