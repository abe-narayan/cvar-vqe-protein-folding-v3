# PRE-REGISTERED PROTOCOL AND FREEZE

**Written 2026-09-05, before the outcome of the remaining in-flight runs was known.** That timing is
the entire point of this document: a decision rule chosen after seeing which arm won is not a
decision rule. The rule is fixed in §2; §3 records which arms it selected once the runs landed, and
§3 is the only section written afterwards.

---

## 1. WHAT THE FREEZE IS FOR

The **60-target benchmark has not been read in this sprint.** No arm, no agent, and no exploratory
run has touched it. It is a **single-use instrument**: once a number is taken from it, any further
tuning against it destroys its value, so it is spent exactly once, on a protocol fixed in advance.

Everything measured so far is **exploratory**, on the 126-target tuning instrument. The programme has
evaluated a large number of arms without family-wise error control, which is the sharpest criticism
the adversarial reviews produced and the one with no adequate answer. A held-out confirmatory
instrument is the only real defence, and it only works if the protocol is frozen first.

---

## 2. THE DECISION RULE, FIXED IN ADVANCE

### 2.1 What may enter the frozen protocol

An intervention is included **if and only if** all of the following hold on the **126-target tuning
instrument**:

1. it is **native-free at inference** — the predictive path is structurally incapable of reading
   native coordinates, and no hyperparameter of it was chosen using the native except through the
   pinned leave-fold-out folds;
2. it beats its own **matched control** with a **paired 95% confidence interval excluding zero**
   (`I.paired`, fold-aware, 4,000 bootstrap resamples);
3. the comparison is **like-for-like on physical validity** — a projected structure is compared to a
   projected structure, never an unprojected consensus to a projected incumbent;
4. it is measured on the **full 126 targets**, not a subset;
5. any hyperparameter it carries was chosen **leave-fold-out**, and the fold being reported never
   contributed to the choice;
6. where it aggregates a set, the set's **diversity is reported** and is not collapsed — the set-mean
   law is out of domain on a collapsed set.

An intervention that fails any one of these is **excluded**, however promising it looks.

### 2.2 The default, if nothing qualifies

**The frozen protocol is the incumbent pipeline, unchanged**: BLOSUM62 retrieval at K = 500 → learned
leave-fold-out distogram filter to m = 75 → coordinate average → projection onto the ideal-geometry
manifold (`ramah`, λ = 0.3, multi-start, exact gradient) → AMBER ff14SB/GBn2 relaxation.

This is a real possible outcome and it is stated in advance so that it cannot later be presented as a
disappointment that was hedged against. **If nothing qualifies, the benchmark run reports the
incumbent's number and the sprint's contribution is entirely mechanistic.**

### 2.3 What is explicitly NOT eligible, and why

- **Every ORACLE arm.** They are ceilings and diagnostics. `ORACLE_scale_true` is worth
  −0.275 [−0.397, −0.174] on the production pipeline's final stage and is **not eligible**, because
  the scale it uses is derived from native coordinates.
- **Any arm whose only support is a smoke read.** Specifically excluded on these grounds unless their
  full-instrument runs land and qualify: the robust and redescending losses (8 targets), the
  surrogate-destruction arms (6 targets), the per-target scale corrections (8 targets), and the
  fusion law's calibration (24 targets). **⚠ Sprint 16, 2026-09-06 (RETRACT): the fusion law's
  "calibration" is withdrawn entirely, at every n. The expression is the Krogh–Vedelsby (1995)
  ambiguity decomposition and is an exact identity (2.65e−15 in a common frame on 126 targets);
  there is nothing to calibrate. See `s16/retract_FINDINGS.md`.**
- **Anything selected by RMSD among multi-starts.** Selection uses the objective only.
- **Any arm that requires removing a target.** The four self-BLOSUM-hit targets and the 16 verbatim
  own-fold windows are disclosed and **retained**.

### 2.4 What will be reported from the benchmark, also fixed in advance

Fixed now so that no field can be chosen after seeing the result:

- **mean and median full-chain Cα-RMSD**, under the frozen RMSD definition (Cα only, every residue
  including both termini, no trimming, correspondence by index, uniform weights, proper rotations
  only, model 1 of the native);
- the fraction under **2.0 Å** and under **2.5 Å**;
- the **paired difference against the incumbent** with a 95% CI, if any intervention qualified;
- the **per-target table in full**, so nothing can be aggregated selectively;
- the count of any structure excluded by the validity gate, and the reason.

**A single PASS/FAIL verdict** is declared against the pre-registered success criterion below. No
other statistic is promoted to a headline afterwards.

### 2.5 The pre-registered success criterion

| outcome | criterion, on the 60-target benchmark |
|---|---|
| **major success** | mean < 2.5 Å |
| **success** | a paired improvement over the incumbent with a CI excluding zero |
| **null** | no significant difference from the incumbent |
| **failure** | significantly worse than the incumbent |

**A mean in the 3.0–3.2 Å range is not a structural success**, per the programme's stated goals,
regardless of how it compares to the incumbent — unless accompanied by a mechanistic result that
stands on its own.

**And an anticipated outcome, recorded now:** on the tuning instrument the generative cascade emits
**3.321 Å against the incumbent's 3.204 Å**, a significant loss, and the measured distance-accuracy
requirement for 2.5 Å is a **22% reduction in distogram RMSE** that nothing in this sprint achieves.
**The expected verdict is `null`.** Writing that down in advance is what makes the eventual number
worth anything.

### 2.6 WHETHER TO SPEND THE BENCHMARK AT ALL — decided in advance

The benchmark is a **single-use instrument**. Spending it has a cost, and that cost is only worth
paying if the run can change what anyone believes. Fixing the rule now, before the last arms land:

**If at least one intervention qualifies under §2.1**, the benchmark is run once, on the frozen
protocol, and reported under §2.4 and §2.5. That is the case the instrument exists for.

**If nothing qualifies, the benchmark is NOT run**, and that decision is itself reported. The reason
is not caution but arithmetic: if the frozen protocol is the incumbent pipeline **unchanged**, then a
confirmatory run measures a system that has not been modified, against a criterion (§2.5) whose
"success" and "major success" tiers are both unreachable by construction — there is no intervention to
be significantly better than the incumbent, because the protocol *is* the incumbent. The run could
only return `null`, which is already known, and it would consume an instrument that a future sprint
with a real candidate will need.

**Spending a single-use instrument to re-measure an unchanged pipeline is not rigour; it is waste
dressed as rigour.** The honest output in that case is: the protocol is frozen, no intervention
cleared the bar on the tuning instrument, the benchmark remains **unspent and sealed**, and the
sprint's contribution is mechanistic.

**What would make this the wrong call**, stated so a reader can disagree with it on the evidence: if
the mechanistic results implied the incumbent's benchmark number should have *moved* — but they do
not; nothing in this sprint changes the incumbent. Or if the benchmark had never been measured on the
incumbent at all — but the protected set exists precisely because it has been reserved, not because it
is unmeasurable. Or if a reader would update on a `null` — but `null` is the pre-registered
expectation in §2.5, so its information content is close to zero.

---

## 3. WHICH ARMS QUALIFIED

*(this section is completed only after the in-flight runs land, by applying §2.1 mechanically)*

Applying §2.1 mechanically as each result landed:

| candidate | control it must beat | verdict |
|---|---|---|
| the generative cascade as a whole | the incumbent | **FAILED** — +0.117 [+0.050, +0.188] |
| consensus rescaling, distogram-referenced | `project_raw` = the production `fit_ca` | **FAILED** — +0.123 [+0.065, +0.184] |
| consensus rescaling, pool-referenced | same | **FAILED** — +0.095 [+0.033, +0.158]; every native-free scale estimator is uncorrelated with the truth (|r| ≤ 0.106) |
| consensus rescaling, pool + debiased | same | **FAILED** — +0.080 [+0.035, +0.129], and **at the ~0.08 Å false-positive floor** |
| robust / redescending losses, incl. graduated non-convexity | `squared` | **FAILED** — −0.011 [−0.102, +0.077] and −0.025 at n = 126, firing the module's own falsification condition |
| alignment-engineered fitting (22 arms) | the unmodified fit | **FAILED** — best leave-fold-out arm −0.007 [−0.074, +0.055] |
| terminal-relaxed fitting | the unmodified fit | **FAILED** — dropping terminal restraints costs +0.210 [+0.055, +0.369]; even *perfect native* termini are worth −0.000 |
| isotropic Tikhonov regularisation (`tik_LFO`) | the unmodified fit | **clears its own control** at −0.216 [−0.345, −0.094], W/L 80/46, 5/5 folds — **but fails §2.1 clause 2 against the protocol's control**: it improves a *losing* arm from 3.677 Å to ≈3.46 Å, still short of the incumbent's 3.204 Å, and its λ grid is not converged (λ = 30 chosen at the boundary on 5/5 folds) |
| restraint-level channel fusion | the better single channel | **FAILED** — −0.071 [−0.156, +0.014] |
| coordinate-level channel fusion | the better single channel | **clears its own control** at −0.255 [−0.370, −0.145] — **but is a comparison between two generative arms**, both worse than the incumbent (3.367 Å against 3.204 Å), so it improves nothing the pipeline emits |
| pool augmentation of the retrieval set | `incumbent_avg` | **FAILED** — +0.004 [−0.001, +0.011] leave-fold-out, and the procedure **chose w = 0 on 4 of 5 folds**. Its premise was false: the fits are 0.108 Å *worse* than the top-75 set mean they would join (3.659 against 3.551), not better. Even the ORACLE ceiling is −0.024 Å, below the floor |
| the quantum ensemble readout | best-of-N, matched budget **and** matched diversity | **FAILED** — demoted by replication: null at the target level, concentration check failed, control charged 2,048 evaluations against the arm's 819,200 |

Already excluded by §2.1 or §2.3: separation debiasing on the ranking axis, `sd` recalibration
(provably inert), every ORACLE arm, and every arm supported only by a smoke read.

---

## 3.1 THE VERDICT

**No intervention qualified. Not one.** Thirteen candidates were tested against their own controls on
the full 126-target instrument; twelve failed outright, and the two that cleared *their own* control
(`tik_LFO` at −0.216 Å and coordinate-level channel fusion at −0.255 Å) both improve a **losing** arm
without reaching the incumbent, so neither clears §2.1 clause 2 against the protocol's control.

> ### THE FROZEN PROTOCOL IS THE INCUMBENT PIPELINE, UNCHANGED
>
> BLOSUM62 retrieval at **K = 500** → learned leave-fold-out distogram filter to **m = 75** →
> coordinate average → projection onto the ideal-geometry manifold (`ramah`, λ = 0.3, multi-start,
> exact gradient) → AMBER ff14SB/GBn2 relaxation.
>
> Emitting **3.204 Å** on the tuning instrument.

## 3.2 AND THE BENCHMARK IS NOT SPENT

By §2.6, fixed in advance: **the 60-target benchmark remains unread and sealed.**

The frozen protocol is the incumbent **unmodified**, so a confirmatory run would measure a system
this sprint did not change, against a criterion whose "success" and "major success" tiers are
unreachable by construction — there is no intervention to be significantly better than the incumbent,
because the protocol *is* the incumbent. It could only return the **`null`** that §2.5 pre-registered
as the expected verdict. Spending a single-use instrument to re-measure an unchanged pipeline is not
rigour; it is waste dressed as rigour, and a future sprint with a real candidate will need it.

**This decision was written down before the last three candidates returned** (§2.6, timestamped in the
document header), which is the only thing that makes it a decision rather than a rationalisation.

**What would reverse it:** any single intervention clearing §2.1 against the incumbent on the tuning
instrument. The nearest miss is `tik_LFO`, which needs to be re-based on the incumbent rather than on
the generative fit, and needs a converged λ grid (λ = 30 was selected at the boundary on 5/5 folds).
That is a well-defined next experiment, not a hope.

---

## 4. THE ENVIRONMENT, FROZEN

- `PYTHONHASHSEED=0`; `OMP_NUM_THREADS = MKL_NUM_THREADS = OPENBLAS_NUM_THREADS = 2`;
  `torch.set_num_threads(2)`. The thread count is stated because the distogram is **not** bit-stable
  across thread counts (`prob` up to 2.0e-6, `risk` up to 1.3e-4), though no headline constant moves.
- All multi-start seeding via `s15/seed.py`'s `blake2b`-based `stable_seed`, verified identical
  across separate interpreters. **The frozen protocol is re-run under stable seeding before the
  benchmark is touched**; exploratory numbers produced before that fix are correct but not
  bit-reproducible, and are not quoted as reproducible constants.
- The instrument must emit `shipped 3.4540004952559396, pool_best 1.7108244199364904,
  top75_best 2.3061526409453816, synthesis_fit 3.2040761603809194, n_zero_recall 18` immediately
  before the benchmark run. (Noting that this is a cache check, not a reproduction.)
- Binding data rule for any AMBER-labelled tail, in-band, argmin or top-k statistic:
  `amber_kind == 0 AND amber_idx != snap_index`, with per-target n printed.

---

## 5. WHAT WOULD MAKE THIS FREEZE INVALID

Stated so that a reader can check rather than trust:

- selecting an arm after seeing any benchmark number;
- reporting a statistic from §2.4 selectively, or promoting one not listed there;
- re-running the benchmark after a change of any kind;
- excluding a target for any reason;
- quoting an ORACLE arm as a predictive result;
- comparing an unprojected structure to a projected incumbent.

If any of these occurs, the benchmark is spent and the resulting number must be reported as
exploratory, not confirmatory.
