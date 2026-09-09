# WORKSTREAM D — FINDINGS, SPRINT 23

Role: adversary, auditor, statistician. Everything below is independently reproduced from source
(scripts in `s23/`, artefacts in `s23/results/`) or independently derived. Findings that damage a
claim already made (the coordinator's n=40 scoping result, or Workstream A's `gscale.py`) are
marked **DAMAGE**. Findings that damage MY OWN hypothesis are marked **SELF-REFUTED**.

**Context at start of session:** Workstream A had already run `s23/gscale.py` at n=126
(`s23/results/gscale.json`, complete). Its headline: incumbent 3.0483, per-target grid oracle
2.9110 (−0.1374, "CEILING"), global scale nested-CV 3.0486 (**+0.0002, NOT MEASURED — a null**).
This document audits that run and extends it to answer BRIEF's duties 1 and 2 in full.

---

## 1. REPRODUCTION AT n=126 (BRIEF duty 1, first clause)

`s23/d_scale_closedform.py` independently rebuilds the incumbent (score-filter → top-75 coordinate
average, `s12.instrument`) from the raw pool and universe files, not from A's cached curves.

    incumbent, s=1.0 (independent repro)        3.0483   <- matches gscale.json to 4 decimals
    A's oracle, grid [0.90,1.10]                2.9110   <- reproduced EXACTLY

**Reproduction is clean — no arithmetic bug in `gscale.py`.** The defect found below is a *design*
choice (grid width), not a coding error.

---

## 2. THE GRID-TRUNCATION DEFECT — **DAMAGE**, load-bearing

`gscale.py`'s oracle grid is `np.linspace(0.90, 1.10, 81)`, centred on the coordinator's n=40
observation of `s* = 0.975`. Checking the returned `s_star` array directly:

    targets with s_star == 0.90 (low boundary)   42/126
    targets with s_star == 1.10 (high boundary)  21/126
    TOTAL AT THE GRID BOUNDARY                   63/126  (50%)

Half the panel's "oracle" scale is clipped at the edge of the search range, not a true optimum.
Replacing the grid with the EXACT closed-form scalar-Procrustes solution (derived and verified
below) gives:

    oracle, EXACT closed-form (unconstrained s>0)   2.7080
    delta vs incumbent                              -0.3403   <- the TRUE ceiling

**This is 2.5× A's reported −0.1374, and over 2× the coordinator's original n=40 figure of
−0.158.** `gscale.py`'s headline oracle number for this sprint is an undercount and should not be
quoted without this correction. (It does not change §5's verdict on the *global* arm — if anything
it sharpens it, since the true ceiling is bigger and the global capture fraction correspondingly
smaller in relative terms.)

**The derivation, checked not asserted (BRIEF's explicit instruction).** Kabsch superposition finds
rotation `R` minimising `Σ|R(C_i−C̄)−(T_i−T̄)|²` via SVD of `H = Σ(C_i−C̄)(T_i−T̄)ᵀ`. Scaling `C` by a
positive scalar `s` scales `H` by `s` but leaves its singular VECTORS unchanged (only singular
VALUES scale) — **so the Kabsch rotation is scale-invariant.** Fixing `R` once at `s=1`, with
`C'=R(C−C̄)`, `T'=T−T̄`:

    SSE(s) = Cc·s² − 2·Ct·s + Tt     (an exact PARABOLA in s; Cc=Σ|C'|², Ct=ΣC'·T', Tt=Σ|T'|²)
    s*_analytic = Ct / Cc            (closed form, no grid needed)

Verified numerically for all 126 targets: 107/126 matched a brute-force check to `rtol=1e-6`
directly; the other 19 were traced individually and are ALL explained either by the sanity-check
grid itself being too narrow (13 targets whose true `s*` lies outside `[0.5,1.5]`, up to `1.92`) or
by sub-grid-resolution rounding (6 targets, agreeing to `<2e-4` in `s`) — **no bug**.

---

## 3. THE SELECTION-ARTEFACT QUESTION — DERIVED, NOT ASSERTED (BRIEF duty 1, second clause)

**The question posed:** is the continuous per-target oracle's near-universal win a selection
artefact, and — critically — is the right null the SAME as the project's existing min-of-K
framework (S22 M3/M4), or something different?

**Answer: something different, and the difference matters.**

`SSE(1) − SSE(s*) = (Cc−Ct)² / Cc ≥ 0` **always**, with equality iff `s*=1` exactly (a
measure-zero event for continuous data). This is **not** a min-of-K multiple-comparisons problem —
there is no discrete menu of K candidate scales being searched, so the S22 M3/M4 apparatus
("permutation null on a minimum is mis-specified because arms correlate through difficulty," "a
router-construction search carries its own best-of-K exposure") does not directly apply. **The
correct analogy is a single-parameter ordinary-least-squares fit evaluated in-sample**: fitting one
scalar by minimising the very residual that is then reported as the "result" mechanically improves
that residual, exactly as any 1-parameter linear regression's in-sample R² is upward-biased.
`126W/0L` in `gscale.py` is therefore **expected under the null of zero real per-atom
correspondence** and is not, by itself, evidence of anything. **The correct correction is NOT a
permutation null on the minimum (S22 M3 already established that construction is the wrong
instrument for a deterministic per-target minimum) — it is a held-out/transfer test**, exactly what
BRIEF asks for next, and exactly S22 M3's own conclusion ("held-out folds are the correct
construction, not a patch") applied here to a continuous parameter instead of a discrete arm.

---

## 4. THE TRANSFER TEST — THE DECISIVE RESULT (BRIEF duty 1, third clause)

`s23/d_scale_transfer.py`, mirroring `s22/mreal.py` exactly (same `stable_rng` seeding convention,
same half-pool-of-250 split, same 8-repeat/2-direction/fold-clustered-bootstrap design). Each
target's K=500 pool splits into two disjoint halves of 250; within each half, apply the identical
score-filter + top-75 construction; fit the closed-form oracle scale `s*_A` on half A **against
native**, freeze it, apply it **unrefit** to half B's independently-built average, and score against
native. Compared to fixed `s=1.0` on the same half-B data (matched-data null).

**I registered a directional hypothesis that this would NOT fully transfer** (by analogy with the
m-ladder's partial 65% transfer). **It transferred almost completely — the falsifier did not fire:**

| panel | fixed (held out) | transfer s* (held out) | in-sample oracle | PRIMARY (transfer−fixed) | fold CI | W/L | % of ceiling that transfers |
|---|---|---|---|---|---|---|---|
| **full (n=126)** | 3.0712 | 2.7514 | 2.7477 | **−0.3198** SE 0.0485 MDE 0.136 (2.35×) | fold [−0.392,−0.235] | 115W/11L | **99%** |
| FAIL18 (n=18) | 5.7157 | 4.5035 | 4.4981 | **−1.2122** SE 0.1942 MDE 0.544 (2.23×) | fold [−1.643,−0.795] | 17W/1L | **100%** |
| non-FAIL18 (n=108) | 2.6304 | 2.4594 | 2.4560 | **−0.1710** SE 0.0275 MDE 0.077 (2.22×) | fold [−0.214,−0.127] | 98W/10L | **98%** |

Every effect clears its own MDE by more than 2×, on BOTH iid and fold-clustered CIs, on the full
panel and on the non-FAIL18 subset separately. **This is a materially cleaner transfer than S22's
m-ladder (65%) or any router result last sprint (0%, several significantly harmful).**
**CONCLUSION: the per-target propensity for a scale correction is REAL, not an in-sample fitting
artefact — it survives this project's own strongest evidentiary bar.**

**Caveat, stated plainly.** This is still an oracle-information test (both halves use native RMSD
to fit or score), exactly like `mreal.py`. It proves the SIGNAL is real and reusable across
independent draws from a target's own pool; it does **not** by itself supply a native-free way to
find `s*` at inference. See §7.

---

## 5. THE GLOBAL-SCALE CAPTURE FRACTION (BRIEF duty 2) — A'S NULL, INDEPENDENTLY CONFIRMED AND
## SHARPENED

**A's `gscale.py` nested-CV result, independently recomputed from the raw JSON (not the log):**
`+0.00022`, matching the printed `+0.0002` to the last reported digit — **no leakage bug in A's
nested construction.** The nested-CV loop (`s_hat` chosen only from `Cv[tr]`, applied to `Cv[te]`)
is correctly nested; the "in-sample, leaky" arm is correctly and separately labelled.

**Duty 2 asks for the fraction a global constant can be EXPECTED to capture, derived from the
spread of s\*, not just the one achieved number.** `s23/d_scale_globalcapture.py` computes the
population-optimal global scale directly from the exact closed-form quadratics (no grid at all,
search range `[0.3,2.2]` — comfortably covers every observed `s_analytic`, min 0.252, max 1.924),
against the TRUE (non-truncated) ceiling from §2:

| panel | incumbent | best-POSSIBLE global s (in-sample upper bound) | achieved | true oracle ceiling | max capture |
|---|---|---|---|---|---|
| full (n=126) | 3.0483 | s=0.9735 → 3.0424 | −0.0059 | −0.3403 | **1.7%** |
| FAIL18 (n=18) | 5.8319 | s=0.6770 → 5.4817 | −0.3502 | −1.3594 | 25.8% |
| non-FAIL18 (n=108) | 2.5844 | s=0.9905 → 2.5835 | −0.0009 | −0.1705 | **0.6%** |

This is the **best a global constant could EVER do, in-sample, evaluated on all 126 targets at
once** — a strict upper bound on what any nested-CV version (necessarily noisier, trained on 4/5 of
the data) can achieve. **It is 1.7% on the full panel. A's nested-CV +0.0002 (statistically a null)
is therefore not a fitting failure — the ceiling for this arm's own best case is already almost
zero.** Further hyperparameter search, alternate CV schemes, or seed variation on this exact arm
cannot recover more than 1.7% of the ceiling; investment here should stop.

**Mechanism, not just magnitude.** `s*` splits **53/126 wanting expansion (s>1) vs 73/126 wanting
contraction (s<1)** — the corrections partially CANCEL when averaged into one constant. This is the
BRIEF's own predicted failure mode ("if the spread is wide, a global constant captures little") —
confirmed, with the sharper diagnosis that it is not spread alone but **direction heterogeneity**
that kills it.

**A length-conditioned scale (8 buckets, sizes 9–23, in-sample only) does slightly better — 7.6%
capture — but this is explicitly the class the S22 finite-sample bound already forbids at these
bucket sizes** (all buckets far below the `n≈100` threshold for even a single global threshold's own
generalisation gap of 0.39 Å; see BRIEF's standing constraint). **Not pursued further; flagged, not
built.**

---

## 6. THE PLACEBO — A SELF-CAUGHT DEFECT, CORRECTED SAME SESSION

**First attempt (`d_scale_placebo.py`), WITHDRAWN.** Permuting the row order of `nat` within a
target destroys coherent chain geometry and degrades the base RMSD to ~7.5 Å — a different-
difficulty problem from the real 3.05 Å pairing, making any improvement-fraction comparison invalid.
Its headline ("387% mechanical share") is **not used anywhere in this document** and should not be
cited. Recorded here per project convention (S22 M7) rather than silently deleted.

**Corrected construction (`d_scale_placebo2.py`).** Pair each target's own averaged pool structure
with the **native structure of a different, randomly chosen target of the same chain length**
(a real, coherently folded chain — same difficulty class, wrong identity). 20 draws/target,
121/126 targets have a same-length peer.

    REAL pairing     : 3.0483 -> 2.7080   delta -0.3403  (11.2% relative)
    PLACEBO pairing  : 4.1478 -> 3.6146   delta -0.5332  (12.9% relative)

**The true native pairing is NOT distinguishable from a random same-length native by this fit** —
the placebo's absolute AND relative improvement both meet or exceed the real pairing's. **This does
not contradict §4's transfer result** (the signal that transfers between pool-halves of the SAME
target is real and reproducible) — it reframes what that signal IS. §4 proves `s*` is a stable
property that a target's own two independent pool-draws AGREE on. §6 shows that agreement is **not
primarily about matching this target's own specific fold/identity** — a scale fit calibrated
against almost any real protein of the same length recovers as much of the gap. What is stable is
closer to "how badly this target's OWN retrieval/pool construction is size-distorted," which
correlates strongly with length and, even more strongly, with pool quality (§7).

---

## 7. CONCENTRATION — FLAGGED IN THE BRIEF'S OWN WORDS

    FAIL18 (18/126, 14% of panel):     ceiling contribution -1.3594 x 18 = -24.47 (Å.targets)
    non-FAIL18 (108/126, 86%):         ceiling contribution -0.1705 x 108 = -18.41
    FAIL18 SHARE OF THE TOTAL CEILING MAGNITUDE:  57%

**18 already-known zero-recall targets (S10-1's FAIL18, unchanged membership) carry more than half
the true scale-correction ceiling.** This is the same shape the project has hit before (S22 L1's
routing ceiling: 50% in the top 21/126 targets) but sharper: here a SINGLE, PRE-EXISTING,
already-diagnosed failure category (not a post-hoc quantile cut) explains the majority of the
effect. **Practical consequence:** any global or length-keyed scale arm is chasing a signal
overwhelmingly parked in a subset the project already knows is pathological for reasons unrelated
to scale (pool retrieval failure) — which is exactly why no shared constant can serve both FAIL18
and the other 108 without hurting one of them (§5's direction-cancellation mechanism).

---

## 8. LEAKAGE / CORRECTNESS AUDIT OF `gscale.py` (BRIEF duty 3)

- **Nested CV: correctly nested.** `s_hat` is selected only from `Cv[tr].mean(0)` (training folds'
  curves); applied unchanged to `Cv[te]`. Independently recomputed from raw JSON, matches to 5
  decimals (§5). No leakage.
- **In-sample arm correctly labelled** "leaky, for reference" and never used as the achieved number.
- **Native information:** used only to build the training-fold curves (legitimate — training folds
  may use native RMSD) and the CEILING rows (explicitly labelled, never deployed). Not used on the
  held-out fold's own data to pick its own `s`. Clean.
- **Grid-truncation defect: real, load-bearing, documented in §2.** Not a leakage issue, a design
  scope issue — does not invalidate the nested-CV null (§5), but does invalidate the `−0.1374`
  ceiling figure as quoted.
- **Completion flag**, minor process gap: `need = ("curve","rmsd_1","s_star","vbond","rg")` checks
  KEY PRESENCE, not `np.isfinite` (contrast `s22/mreal.py`, which checks both). No bad value was
  actually found in this run — flagged as a process gap for future arms, not a defect in this one.
- **Tie-breaking:** `np.argmin` over an 81-point mean-curve for nested selection and over the
  score-filter's `argsort(..., kind="stable")`. Checked for exact float ties in the training-curve
  argmin across all 5 folds — none found (continuous aggregate over ~100 targets). Low risk, no
  action needed.
- **Basis discipline:** point cloud throughout, on both sides of every comparison in `gscale.py` and
  in this document's extensions. No point-cloud/built-chain mixing found.
- **Frame consistency:** confirmed analytically and numerically that Kabsch rotation is scale-
  invariant (§2), so "scale about the centroid, then re-Kabsch" (`gscale.py`'s method) and "fix the
  rotation once, then scale in that frame" (this document's closed form) are the SAME operation —
  cross-validated by exact numerical agreement (§1).
- **Torsion wrapping / degrees-vs-radians:** not applicable — this arm is pure Cartesian, no
  torsions involved.

**No other lane's code exists yet in `s23/results/` beyond `gscale.py`** (checked at time of
writing); this audit will need to continue as `PREREG_A`'s Experiments 1–3 and any other lane's
arms land.

---

## 9. THE GEOMETRY-COST STANDING CONSTRAINT

The nested-CV global arm's own fitted scale is `s≈0.974` — a further CONTRACTION of an incumbent
already 22% short on virtual Cα–Cα bond length (2.961 Å vs physical 3.812 Å, confirmed independently
in §1's reproduction). **Flagged per BRIEF's explicit instruction, even though the arm is a null and
therefore not being promoted:** had this arm been adopted despite being statistically indistinguishable
from zero, it would have bought nothing in RMSD while making the geometry marginally less physical
still. Since the arm does not ship (§5), no actual trade is being made — recorded so the pattern is
on file before any future arm in this family is evaluated.

---

## 10. STATISTICAL GATEKEEPING SUMMARY (BRIEF duty 4)

Every primary comparison above carries: paired per-target differences (yes, throughout), iid AND
fold-clustered bootstrap CIs (yes, §4; §5's global-CV number reused from A, itself fold-clustered in
`gscale.py`), SE and MDE reported as a multiple (yes: §4's transfers clear 2.2–2.35× their own MDE),
W/L (yes), worst-target degradation (nested global arm: `+0.197 Å` worst case, from A's own log,
independently spot-checked as plausible against the curve data), and an ablation-equivalent (the
FAIL18/non-FAIL18 split throughout, and the grid-vs-closed-form comparison in §2). Seed sensitivity
not separately re-run this pass (both `d_scale_transfer.py` and A's `gscale.py` use
`s15.seed.stable_rng`, this project's standard deterministic-seed convention — not itself a defect,
but a full seed-perturbation sweep was not budgeted this pass and is flagged as not done, not
silently assumed clean).

---

## 11. VERDICT AND RECOMMENDATION TO THE COORDINATOR

1. **Fix `gscale.py`'s grid before the `−0.1374` (or the original n=40 `−0.158`) ceiling number is
   quoted again.** The true unconstrained ceiling is `−0.3403 Å`. This makes the ceiling BIGGER, not
   smaller — it does not rescue the global-scale arm.
2. **The global-constant scale arm (H1 as currently scoped) is REFUTED, not just unmeasured.**
   Independently confirmed at its own theoretical maximum: 1.7% of the true ceiling, in-sample,
   best case. A's nested-CV `+0.0002` was already at that ceiling — there is no room to find more of
   this signal by tuning the same arm. **Recommend closing H1's global-scalar sub-arm entirely.**
3. **Do not pursue a per-length or any richer per-target router for this signal.** 7.6% in-sample
   capture at bucket sizes (9–23) far below the finite-sample bound's own n≈100 threshold — the
   bound BRIEF invokes forbids exactly this construction before it is tried.
4. **The real, transferring signal (§4, 98–100%) is concentrated in FAIL18 (§7, 57% of the ceiling)
   and is NOT primarily an identity-match effect (§6)** — it looks like a symptom of pool-retrieval
   failure for those 18 targets, wearing a "scale" costume. **This points toward diagnosing WHY
   FAIL18's pools are size-distorted (a different investigation, out of scope for this file) rather
   than toward any generic scale-correction architecture.**
5. **If any lane wants to bank the real non-FAIL18 residual (§4: −0.171 Å, transfers at 98%),** that
   requires a native-free ESTIMATOR of `s*` per target — exactly the per-target adaptive-rule class
   the finite-sample bound (S22 C1) forbids at n=126. **This closes the loop the same way S22's
   m-ladder and every router closed: a real, transferable, oracle-confirmed signal that is
   structurally unreachable at this sample size by anything richer than what is already refuted in
   §5.**

**Self-assessment against my own registered hypothesis (PREREG_D):** I expected partial transfer,
by analogy with the m-ladder. **I was wrong — the transfer is nearly total (98–100%), the cleanest
result of this kind in the project's history.** I also built and then withdrew a confounded placebo
mid-session rather than publish it. Both are recorded above rather than smoothed over.
