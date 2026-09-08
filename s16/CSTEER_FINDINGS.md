# SPRINT 16 / THE PIVOT — COORDINATE-SPACE STEERING

**Module** `s16/csteer.py` · **artefact** `s16/results/csteer.json` · **log** `s16/results/csteer.log`
**Ledger** L11 (the pivot and its falsifier, both written before the run) · **n = 126**

---

> ## ⚠ CORRECTED 2026-09-06 BY THE VERIFY WORKSTREAM — READ THIS FIRST
>
> **§1 is demoted and §2's stated reason is RETRACTED.** The numbers below all reproduce from the
> artefact exactly; the error is in what they were read as. Both corrections are kept in place with
> the originals struck through, per the programme's standard. See ledger L21 and
> `s16/verify_FINDINGS.md`.
>
> - **§1 is a tautology, not an exhibit.** Max deviation from `(1 − s)·rmsd` on 400 random
>   point-cloud pairs *with no protein in them* is 4.6e−07 Å. Kabsch cannot leave the identity along
>   the segment because a sum of two symmetric PSD cross-covariances is PSD. The module's own
>   docstring already writes the algebra out. **An exhibit derivable on random point clouds is not
>   an exhibit.**
> - **§2 measured the wrong functional.** The module's law is **two-sided and quadratic in c**, so
>   the quantity it consumes is **\|c\|**, not the signed mean. Read correctly, the direction is
>   *present*: \|c\| = 0.305 against a 0.128 isotropic null (+0.177 [+0.142, +0.212], W/L 100/26).
>   What is missing is the **per-target sign** — one bit — **and the step magnitude `‖r‖`, which is
>   the RMSD to the native and which nothing in this sprint estimates native-free.** §2 is
>   rewritten below.
> - The **operational verdict is unchanged**: nothing here is deployable and 2.5 Å is not reached.

## The verdict in one line

**Also refuted — but not for the reason first published.** The pivot's stated falsifier was "every
native-free channel's cosine is at or below the torsion-space 0.318". ~~What actually happened is
worse: from the best available starting structure, every native-free direction has a cosine with the
true residual that is statistically indistinguishable from zero.~~ **[RETRACTED — see the box
above.]** The corrected statement: from the best start the direction's **magnitude** is 2–3× the
isotropic null, its **sign flips per target**, and a **zero-information constant β-strand carries
more of it than any channel this sprint built**. Every arm that *appeared* to win did so by walking
the whole way to a structure the pipeline already had.

---

## 1. ~~The mechanism claim is confirmed exactly, and it is the sprint's cleanest exhibit~~ — DEMOTED TO A THEOREM CHECK

The pivot rested on the claim that coordinate space is a space in which partial motion toward the
truth pays proportionally, where torsion space is not. Stepping the **exact** residual, same 126
targets, same fits, same superposition:

| step s | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.65 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|
| **coordinate space** (from the fit) | 3.647 | 3.282 | 2.918 | 2.553 | 2.188 | 1.824 | 1.277 | 0.729 | 0.000 |
| exact prediction 3.647·(1 − s) | 3.647 | 3.282 | 2.918 | 2.553 | 2.188 | 1.824 | 1.277 | 0.729 | 0.000 |
| **torsion space**, same targets | 3.647 | 3.599 | 3.641 | — | — | 3.466 | — | — | 0.043 |

Coordinate space is **exactly linear**, to three decimals at every rung, and the identity survives
Kabsch superposition. Torsion space is non-monotonic over the first third and reaches 5% of the gain
at half-way. The two rows are the same error, on the same structures, measured the same way; the
only difference is the space the step is taken in.

~~This is a real and, as far as the programme's own record goes, previously unmeasured fact.~~
**[DEMOTED 2026-09-06.]** The coordinate row is a **tautology**: stepping exactly toward a target
and re-superposing must be linear, because a sum of two symmetric PSD cross-covariances is PSD, so
Kabsch never leaves the identity along the segment. VERIFY reproduced it on **400 random point-cloud
pairs with no protein in them**, max deviation **4.6e−07 Å**. The module's own docstring writes the
algebra out, which makes publishing the table as an exhibit the same error RETRACT retracted the
fusion law for, on the same day.

**What content remains in the two rows together?** Only the torsion row, and it is a measurement:
that the *same* error, on the *same* structures, behaves non-monotonically when the step is taken in
torsion space. §3 of `s16/STEER_FINDINGS.md` carries that result and its own — later corrected —
mechanism claim.

## 2. ~~But the direction does not exist~~ — RETRACTED; the magnitude is there, the SIGN is not

Cosine with the true residual, per starting structure, native-free directions only:

| start | mean RMSD | direction | cos [95% CI] | LFO step | effect |
|---|---|---|---|---|---|
| **`avg`** — top-75 coordinate average | **3.048** | toward the fit | **0.012 [−0.052, +0.075]** | 0.0 | +0.006 |
| | | toward the projection | **−0.068 [−0.107, −0.027]** | 0.0 | +0.000 |
| | | toward the medoid | **−0.010 [−0.064, +0.043]** | 0.0 | +0.003 |
| | | toward the ensemble mean | **−0.004 [−0.066, +0.058]** | 0.0 | +0.000 |
| | | *random control* | 0.005 [−0.022, +0.033] | 0.0 | +0.000 |
| `fit` — torsion-space fit | 3.647 | toward `avg` | 0.479 [0.425, 0.531] | **1.0** | −0.599 |
| `proj` — projected average | 3.213 | toward `avg` | 0.316 [0.286, 0.347] | **1.0** | −0.164 |

**The requirement, from the incumbent's 3.204 Å: c = 0.625 for 2.5 Å, c = 0.781 for 2.0 Å.**
*(Corrected 2026-09-06: 0.615 as first written, which is the value for a 3.170 Å base. The figure
was hardcoded in `s16/csteer.py` rather than computed; it has been corrected at source.)*

~~**Available from the best start: c = 0.012, with a confidence interval containing zero and a
random control at 0.005.**~~ **[RETRACTED 2026-09-06.]** The law `‖r_new‖ = ‖r‖√(1 − c²)` is
**quadratic and two-sided in c**, so the functional it consumes is **\|c\|**. The signed mean is not
that functional, and quoting it turned a sign problem into a phantom absence of signal.

### §2, corrected. The magnitude is there; the sign is not.

From `avg` (3.048 Å), n = 126, fold-clustered paired CIs, against the module's own random control —
whose sampled \|c\| is **0.136** and whose analytic isotropic null √(2/(π·3n)) at n = 13 is
**0.128** — the paired subtractions below use the 0.128 control mean:

| direction | signed c (as published) | **\|c\|** | **\|c\| − random [95% CI]** | W/L | ORACLE two-sided ceiling |
|---|---|---|---|---|---|
| toward `fit` | +0.012 | **0.305** | **+0.177 [+0.142, +0.212]** | **100/26** | 2.843 Å |
| toward `ens` | −0.004 | 0.301 | +0.173 [+0.135, +0.210] | 98/28 | 2.836 Å |
| toward `med` | −0.010 | 0.251 | +0.123 [+0.087, +0.160] | 83/43 | 2.894 Å |
| toward `proj` | −0.068 | 0.187 | +0.059 [+0.033, +0.087] | 70/56 | 2.970 Å |

**And the zero-information references AUDIT made mandatory beat every channel we built.** These know
nothing about the target but its chain length:

| direction | kind | **\|c\|** | \|c\| − random | ORACLE ceiling |
|---|---|---|---|---|
| **constant ideal β-strand** | **ZERO-INFORMATION** | **0.398** | **+0.262 [+0.216, +0.310]** | **2.491 Å** |
| constant ideal PPII | ZERO-INFORMATION | 0.382 | +0.246 [+0.200, +0.293] | 2.561 Å |
| constant ideal α-helix | ZERO-INFORMATION | 0.335 | +0.199 [+0.155, +0.243] | 2.733 Å |
| **pure isotropic breathing** (scale the start about its own centroid — no input at all) | **ZERO-INFORMATION** | 0.296 | +0.167 [+0.123, +0.214] | 2.708 Å |
| random | control | 0.136 | +0.000 | 3.009 Å |

**The mechanism.** The residual `avg → native` is dominated by a **compactness scalar**. Every
extended reference aligns with it in magnitude, and **58% of targets need the contracting sign, 42%
the expanding one.** This is the standing project record reappearing exactly — *averaging contracts
the backbone 25.8%*, and *the only leverage supplies the per-target SIGN at inference*.

**A second, separable defect in the same module.** `csteer.STEPS` has **no negative rung**. The
fraction of targets whose optimal step from `avg` is negative is 0.45 (`fit`), 0.67 (`proj`), 0.52
(`med`), 0.48 (`ens`) — so on roughly half the instrument **no arm on that ladder could have helped
whatever the direction carried**. A leave-fold-out selector offered only non-negative steps for a
direction whose sign flips per target must return "do nothing", and it did. **The published null is
partly a property of the ladder, not only of the channels.**

**What is and is not deployable.** Every "ORACLE ceiling" above uses the native to place both the
sign and the magnitude of the step, so none of them is an effect. Coordinate-space steering with a
leave-fold-out step and these channels buys nothing, and 2.5 Å is not reached. That verdict is
unchanged. What changes is the reason, and the reason matters: **the missing ingredient is one bit
per target, not the direction** — and the ceiling behind that bit, 2.491 Å from a zero-information
β-strand, sits below the sprint's own primary goal.

**Whether that bit is obtainable native-free: CLOSED, by an existing artefact.** Sprint 15's
`s15/expand.py` already tested four native-free scale references at n = 126. Read at the level of
the sign, their accuracy against the RMSD-optimal scale is **0.357 / 0.389 / 0.413 / 0.389 against a
zero-information constant-sign baseline of 0.516** — *every one of them worse than a constant
guess*, with correlations of ±0.04. Sprint 15 also measured the cost of using them: +0.080 to
+0.123 Å, all CIs excluding zero. See ledger L22.

## 3. Every apparent win is a re-derivation, not a steer

The two arms with large effects both chose **step 1.0 on all five folds** — they walk the entire way
to the target structure and stop there. `fit → avg` at step 1.0 lands at 3.048 Å, which is exactly
the pool coordinate average; `proj → avg` at step 1.0 lands at 3.048 Å, the same object. Neither is
steering. Both say "discard this structure and use the coordinate average", which is what the
incumbent pipeline already does at stage 3.

Reported as effects they would have looked like a −0.599 Å improvement. They are not. **A
leave-fold-out step selector that lands on the boundary of the ladder is a warning that the arm is
choosing an endpoint, not a step, and every such arm here was.**

## 4. One incidental measurement worth keeping

`proj` (3.213 Å) is the ideal-geometry projection of `avg` (3.048 Å), so on this instrument the
**projection costs 0.165 Å of accuracy**. It is not run for accuracy — it buys buildable geometry,
and the REPAIR workstream is pricing that trade properly. But the number belongs in the record
beside the incumbent's 3.204 Å, because it says the incumbent's projection stage is a validity
purchase, not an accuracy one. This is consistent with Sprint 15's `fig16`, which already showed
selection and projection losing ground while aggregation recovers part of it.

## 5. Convergence with the other workstreams

Three independent lines reached the same place on the same day, and none of them was constructed
knowing the others' answers:

- **Flagship** (`s16/steer.py`): stepping the exact torsion error pays 7.7% at half-way, so no
  direction estimate can be steered with.
- **This module**: in the space where stepping *does* pay proportionally, the direction has
  cosine 0.012 [−0.052, +0.075] from the best start.
- **AUDIT** (`s16/audit_FINDINGS.md`): the native-free surrogate's signal is reproduced by a
  **zero-information constant α-helix** (|cos| 0.357 vs the pool's 0.390; the pool adds
  +0.033 [−0.007, +0.075], CI including zero), the correct matched null is **0.216 rather than
  0.157**, and the quiet-subspace identification is refuted — a constant α-helix Jacobian overlaps
  the native's bottom-half subspace at 0.829 against the emitted structure's 0.827.

AUDIT's one surviving item — the surrogate's alignment with the **loud** half of the error, 0.497
against a matched null of 0.304, W/L 101/25 — is closed by the flagship rather than left open:
stepping the **exact** loud component bought −0.028 [−0.106, +0.051].

## 6. What may be quoted

**May be quoted.** The exact linearity of coordinate-space RMSD under motion toward the native, and
its contrast with torsion space, both at n = 126. The cosines in §2, with their random controls. The
projection's 0.165 Å accuracy cost.

**May not be quoted.** The −0.599 Å and −0.164 Å "effects" as improvements; they are endpoint
substitutions (§3) and are reported here only so that the trap is on the record.
