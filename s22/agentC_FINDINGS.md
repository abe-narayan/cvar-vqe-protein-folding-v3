# SPRINT 22 — WORKSTREAM C FINDINGS
## The routing headroom does not survive honest evaluation; error-placement weighting is null

**Pre-registration**: `s22/PREREG_C.md`, frozen before any number in this file existed; the operator
fork list (Rule 0 clause 2) is in that file's §4, sent before either experiment ran, and its §5
dated addendum records the coordinator's mid-flight correction to `s22/LEDGER.md` L2 and how it
changed the design (blend ladder added; hard difficulty-switch dropped) before data existed.

**Instrument**: 126 targets, pinned 5-fold structure, point-cloud/window basis matching the
incumbent (**3.0483 Å**) except where explicitly marked rebuilt-torsion basis (the four latent
arms). Every arm reused from Sprint 21 was **cross-checked bit-for-bit against its source artefact
before being trusted** — `s22/routerdata.py` raises on any mismatch rather than silently diverging;
all 126 targets passed on every checked arm (`n_checked_vs_s21 = 126/126`).

**Bottom line, stated first because both experiments are negative and the brief asks for that to be
reported as plainly as a positive.** Neither priority experiment moves Cα-RMSD. Experiment 1 (the
routing headroom) fails under held-out-fold evaluation, and a second, independent line of attack
(a nonlinear model, checked in-fold vs. cross-validated) shows *why*: the arm-vs-feature
relationship the router needs is not there for a linear model to miss and not reliably there for a
nonlinear one either — a Random Forest wins hugely in-fold (**-0.276 Å**, i.e. more than half the
ceiling) and **loses** under the same held-out-fold discipline (+0.079 Å). This independently
corroborates the coordinator's concern (`s22/LEDGER.md` L3, and the `mreal.py` split-half test in
flight at hand-off) that a meaningful fraction of the 0.482 Å ORACLE ceiling is not a property of
the targets a deployable model can reach — via a completely different method (cross-target,
feature-based CV rather than within-target pool-split-half). Experiment 2 (error-placement
weighting) is a clean, well-powered null: no structural weight beats the plain-precision control.

---

## EXPERIMENT 1 — CAPTURING THE ROUTING HEADROOM: NEGATIVE, AND WHY

### Arm set and ceiling, cross-checked against the coordinator's independent number

`s22/routerdata.py` reassembles, per target: the m-ladder (`avg_500…avg_1`, matches
`s21/poolgap.py` bit-exactly), the two consensus medoids (same source), the four latent arms
(matches `s21/latentsel.py` bit-exactly on `lat_argmin`/`lat_avg75`), and a **new blend ladder**
(`blend_00…blend_100`, five points superposing the latent's top-75 average onto the pool's top-75
average's frame via Kabsch and mixing per-residue — `blend_100` is checked to equal `avg_75`
exactly, confirming the construction). 16 total candidate arms (the coordinator's independent count
was 13 over a narrower set; my ORACLE ceiling on the superset, **2.5618**, lands within 0.005 Å of
the coordinator's **2.566** — an independent cross-check via a different implementation path, not
an attempt to reproduce their exact number).

| arm set | n arms | ORACLE ceiling | incumbent | headroom |
|---|---|---|---|---|
| m-ladder only | 6 | 2.6985 | 3.0483 | 0.3498 |
| pool_all (World A) | 8 | 2.6628 | 3.0483 | 0.3856 |
| A + latent (World B) | 12 | 2.5664 | 3.0483 | 0.4820 |
| A + blend (new) | 12 | 2.6224 | 3.0483 | 0.4259 |
| **ALL** | **16** | **2.5618** | 3.0483 | 0.4865 |

**Winner-count check reproduces the coordinator's headline finding on the m-ladder exactly**:
`avg_500` (mean 3.396, the *worst* fixed arm) wins the most targets (37/126) — the arm worst on
average wins the most often, which is the signature the coordinator flagged as consistent with
min-of-K noise on a flat, correlated ladder rather than real per-target structure.

### The achieved router, held-out-fold only — NEGATIVE on every arm set tried

Method: one Ridge regression per arm, alpha chosen by generalized CV *within the training folds
only*, fit on four folds, predicting the fifth; routed to `argmin` predicted RMSD; concatenated
out-of-fold. Never scored on a fold that selected it.

| arm set | achieved (CV) | vs incumbent | MDE | capture of headroom |
|---|---|---|---|---|
| m-ladder only | 3.0960 | +0.0477 [-0.0075,+0.1054], 34W/54L | 0.0808 | **-13.6%** |
| pool_all (A) | 3.1078 | +0.0595 [-0.0038,+0.1345], 35W/53L | 0.0982 | **-15.4%** |
| A + latent | 3.1186 | **+0.0703 [+0.0048,+0.1459]**, 34W/55L | 0.0989 | **-14.6%, CI excludes zero — significantly WORSE** |
| A + blend | 3.0806 | +0.0323 [-0.0334,+0.1081], 49W/55L | 0.1040 | **-7.6%** |
| ALL | 3.0917 | +0.0433 [-0.0259,+0.1208], 47W/57L | 0.1040 | **-8.9%** |

**Every arm set is negative** (capture is a negative percentage — the router makes RMSD worse, not
better), and the richest natural set (World A + latent, the set the coordinator's L2 asked me to
treat as first-class) is **significantly** worse than doing nothing. The random-arm control (route
to a matched-count uniformly random arm) is markedly worse still in every set (e.g. ALL: 3.428 vs
incumbent 3.048), which at least confirms the routing PROBLEM is not degenerate — some arms really
are much worse alone — but the learned router does not clear even the much lower bar of "better
than the incumbent," let alone approach the ceiling.

**Ablation**: rg-only (2 features: `rg_z`, `rg_gap`) vs. the full 15-feature set on the ALL arm
set — both negative, both CI-spanning, **not distinguishable from each other**
(rg-only +0.0263 [-0.0455,+0.1104]; full +0.0433 [-0.0239,+0.1213]). The extra engineered features
(score-distribution stats, pool similarity) buy nothing over the two-feature channel that has
cleared every other control in this project, and neither buys anything over the incumbent.

### Attacking the negative before reporting it (`s22/routercv2.py`) — it holds up, and sharpens

**Attack 1 — is this "wins in-fold, loses out-of-fold" (the exact hazard BRIEF names)?** Fit and
predict on the identical 126 targets (unlimited leakage, an upper bound with zero generalisation
gap): **-0.0459 [-0.0943,+0.0002]**, 58W/49L. Even with **no held-out discipline at all**, the
Ridge-per-arm method captures under 10% of the 0.4865 Å ceiling and its CI *still* barely fails to
exclude zero. **The bottleneck is not overfitting alone — the linear signal itself is weak.**

**Attack 2 — does a nonlinear model reveal more signal?** A Random Forest per arm (depth 4, 300
trees), same features, same arm set:

    in-fold (unlimited leakage)     -0.2762 [-0.3681,-0.1943], 80W/29L   -- 57% of the ceiling, "significant"
    held-out fold (proper CV)       +0.0790 [-0.0020,+0.1652], 42W/68L   -- WORSE, direction flips

**This is the textbook in-fold-wins/out-of-fold-loses pattern**, produced deliberately as a
diagnostic and reported in full rather than only the CV number. A tree ensemble with 15 features
over ~100 training targets per fold can memorize enough idiosyncratic structure to look like it
has found a 0.28 Å effect, and that entire effect evaporates — and reverses sign — under the
discipline BRIEF and the coordinator both required.

**Attack 3 — does gating (confidence-thresholded switching) rescue it?** Nested-CV-selected margin
threshold (selected on an inner 70/30 split of the training folds only, never the true held-out
fold), applied to the Ridge router on the ALL set: **+0.0278 [-0.0032,+0.0836]**, and the
nested-selected threshold is **∞ (never switch) on one of five outer folds** and ≥0.3 Å on three of
the other four. The nested procedure itself is telling you the safest policy it can find is close
to "don't route." The gate does not turn the result positive; it mostly disables itself.

### Disposition

**Experiment 1's falsifier — "the achieved held-out-fold router beats the incumbent past its own
MDE with a CI excluding zero" — DOES NOT FIRE, on any of five arm sets, two model classes, two
feature sets, and a gated variant.** One configuration (World A + latent, full features, naive
switching) is **significantly worse** than the incumbent. This is reported as the negative result
BRIEF explicitly asks for, not reframed as a small positive.

**What this adds to the coordinator's own concern.** `s22/LEDGER.md` L3 raised the possibility
(and the in-flight `mreal.py` split-half test is testing directly) that much of the 0.482 Å ORACLE
ceiling is a min-of-K artefact of a flat, correlated arm ladder rather than a real per-target
property. This experiment attacks the same worry from an orthogonal angle — cross-target,
feature-based prediction rather than within-target pool-splitting — and returns the same verdict:
**whatever structure exists in which arm wins per target is not accessible to a model fit on
native-free features of that target, linear or nonlinear, in-fold or out.** Both lines of evidence
should be read together, not as duplicates: a split-half failure would show the *target-level*
signal is unstable; my failure shows that *even a stable target-level signal*, if one exists, is
not recoverable from the features on hand. Either is sufficient on its own to close the router as
currently specified; both failing is a stronger joint statement.

**What is NOT closed.** (i) A held-out-fold discipline with only 5 folds and ~25 targets each is
itself a low-power design for detecting an effect the size of the MDEs above (~0.08–0.12 Å); a
router half this noisy could be hiding under these CIs. (ii) Sequence/ESM features and secondary
structure were declared out of scope for compute/time reasons (PREREG_C.md §1) and are NOT MEASURED,
not tested-and-null. (iii) The blend mechanism itself (§Experiment 1, blend ladder) is not
separately falsified — it was folded into the same router and inherits the same negative, but a
dedicated per-target blend-fraction regression (rather than discrete blend points as arms) was not
attempted and is a smaller, cheaper follow-up than rebuilding the whole router.

### Experiment 3 (readout/window), descriptive only, folded in as pre-registered

Router's chosen arm family by per-target headroom tercile (using the ALL-set CV router's actual
out-of-fold choices — descriptive, since the router itself is a negative result, this shows *what
a losing router does*, not a recommendation):

    low headroom  (n=42): m_ladder 30, blend 7, medoid 4, latent 1
    mid headroom  (n=42): m_ladder 26, blend 7, medoid 7, latent 2
    high headroom (n=42): m_ladder 30, blend 9, latent 2, medoid 1

The router leans on the m-ladder at every headroom level and barely touches the latent arms
regardless of headroom size — consistent with the ablation finding no signal for the router to act
on, not with a real "route to the latent when it's hard" policy it declined to learn.

---

## EXPERIMENT 2 — ERROR-PLACEMENT WEIGHTING `H_D*`: NULL, CLEANLY

`s22/errweight.py`, n=126, four pre-registered native-free weighting schemes multiplying a
Gaussian-approximation precision term, soundness-gated before any structural claim is read.

**Soundness check**: the Gaussian approximation (declared functional fork, PREREG_C.md §2) tracks
the shipped Bayes-risk score reasonably — median Spearman 0.973, argmin agreement 58% — enough to
trust `s_control` as a meaningful null, not a different objective in disguise.

**Certified optimum (pool argmin)**, each structural scheme vs. `s_control`:

    s_local    +0.0301 [-0.0694,+0.1313], MDE 0.1433, 28W/28L
    s_global   +0.0038 [-0.0714,+0.0742], MDE 0.1041, 25W/28L
    s_discrim  +0.0355 [-0.0698,+0.1430], MDE 0.1542, 34W/45L

**Deployed readout (top-75 average)**, the primary per the pre-registration:

    scheme      vs s_control                          vs shipped incumbent (3.048)
    s_local     -0.0108 [-0.0553,+0.0321]              +0.0264 [-0.0327,+0.0857]
    s_global    -0.0165 [-0.0695,+0.0365]              +0.0206 [-0.0350,+0.0787]
    s_discrim   +0.0030 [-0.0600,+0.0672]               +0.0401 [-0.0215,+0.1042]

**None beats `s_control`, and none beats the shipped incumbent — every CI spans zero on both
comparisons.** My own directional prior (`s_discrim`, named in advance per PREREG_C.md §4 clause 5,
precisely so this moment could be checked against it) is numerically the *worst* of the three on
the deployed readout against the incumbent (+0.0401, the largest positive point estimate) — my
stake did not bias the result in its own favour, which is what naming it in advance was for.

### Disposition

**Falsifier — "at least one structural weight beats both `s_control` and the shipped incumbent past
its own MDE on the deployed readout" — DOES NOT FIRE.** Re-weighting the objective's own pairwise
error by chain separation (either direction) or by the pool's own discrimination variance moves
nothing detectably, in either direction, at this sample size. Consistent with the sprint's
standing finding that native-free functionals of the objective's own structure (its score
distribution, and now its pairwise error geometry) do not carry exploitable signal — `the objective
cannot audit itself` (S21 L25/E6) extends, at this power, to a re-weighting of its own error terms.

**What is NOT closed.** The MDEs here (0.06–0.15 Å) are comparable to or larger than the effects
this project's more successful interventions have carried, so a true effect under ~0.1 Å could not
be detected by this design. The four weight forms tested are a small, pre-declared set motivated by
chain separation and pool discrimination; other native-free structural weights (e.g. a genuine
per-pair torsional-Jacobian sensitivity, rather than the chain-separation proxy for it used here)
are NOT MEASURED, not tested-and-null, and were named in the brief as one candidate this file did
not attempt exactly.

---

## ARTEFACTS

All under `s22/results/`, atomic writes (tmp + `os.replace`), config-derived names, completion
flags requiring the full row count and — for `routerdata.json` — an explicit cross-check count
against the two Sprint-21 source artefacts it must reproduce:

    routerdata.json   126/126, complete=true, n_checked_vs_s21=126/126
    routercv.json     complete=true (5 arm-set tables + ablation + descriptive tercile breakdown)
    routercv2.json    complete=true (in-fold and gated-router attacks on the routercv negative)
    errweight.json    126/126, complete=true

Source: `s22/routerdata.py`, `s22/routercv.py`, `s22/routercv2.py`, `s22/errweight.py`.
Pre-registration, including the Rule-0 fork list sent before either experiment ran and the dated
addendum responding to the coordinator's L2 correction: `s22/PREREG_C.md`.

## DATED ADDENDUM — 2026-09-07, AFTER THE ABOVE WAS WRITTEN: L2 IS FULLY WITHDRAWN

The coordinator's `s22/LEDGER.md` L2 (latent-vs-incumbent difficulty interaction) has now been
retracted a second time, algebraically and with a placebo (an artefact-generating shuffle
reproduces an even larger version of the same sign-flip shape from zero real information), and the
"independent stratifier" correction I built §5(b)/(c) of `PREREG_C.md` around is retracted along
with it — the same slope-difference identity (`Cov(stratifier, lat) - Cov(stratifier, inc)`)
applies to *any* difficulty stratifier, `pool_oracle` included. **The one thing that survives is a
single global fact — the latent's RMSD is less difficulty-sensitive than the incumbent's (slope
0.58 vs 1) — which is not a regime, has no crossover, and is not actionable by a native-free router
since exploiting it requires the incumbent's own oracle RMSD to place a target in a regime at all.**

**This does not change anything I built or reported above, and it is worth stating precisely why.**
My router (`s22/routercv.py`) never implemented a difficulty-quartile switch — it fits one Ridge
(and, in the attack, one Random Forest) regression per ARM directly against native-free features
(`rg_z, rg_gap, pool_spread`, score-distribution stats, length), and routes by `argmin` over
predicted RMSD. The four latent arms entered as ordinary candidates in that regression, not through
a hand-built "detect hard targets, switch to latent" rule. So the mechanism now withdrawn was never
load-bearing in my design; PREREG_C.md §5(b)'s blend-ladder addition was motivated by the
now-retracted interaction, but the blend ladder was tested as ordinary router candidates through
the same feature-based regression, not through any difficulty-placement logic — its own negative
result (`A_plus_blend`: +0.032 [-0.033,+0.108], -7.6% of headroom) stands on its own construction
and needs no revision.

**The retraction reinforces rather than undermines my negative.** If anything, a mechanically-
inflated interaction of the kind just retracted is exactly the sort of spurious pattern a flexible
regression (especially the Random Forest attack) could have latched onto and reported as a false
positive; instead both model classes returned a negative under held-out-fold discipline, and the
in-fold/out-of-fold divergence I already reported (RF: -0.276 in-fold vs +0.079 held-out) is now a
second, convergent illustration of the same class of artefact the coordinator's algebra names for a
different statistic.

**Third convergent negative, from the audit lane's own achievable test** (reported in the
coordinator's message, not run by me): a 5-fold held-out router using peptide length to pick the
per-quartile best arm achieved 3.079 against the fixed 3.048 (+0.031 [MDE 0.070], NOT MEASURED,
wrong-signed) — approximately 0% of the ceiling captured, on a *third* independent construction
(length-quartile arm selection) beside my two (Ridge/RF per-arm regression) and the coordinator's
concern about the m-ladder's own min-of-K structure. Three different router constructions, three
negative or null results.

**Action taken in response to the coordinator's explicit instruction**: no new router is being
built, and I am not spending further time hunting for a difficulty regime. Experiment 1's
disposition is unchanged from what is written above — NEGATIVE, held-out-fold, falsifier does not
fire — and is now additionally understood to not depend on any mechanism the coordinator has since
withdrawn. `PREREG_C.md` is left unedited (per its own rule against post-hoc rewriting); this
addendum is the record of the correction's effect on this file's conclusions, and none of Experiment
1's tables or dispositions required updating as a result.

## DATED ADDENDUM 2 — 2026-09-07: `mreal.py` REFUTED THE NOISE HYPOTHESIS; THE COORDINATOR'S
## GEOMETRY-PROFILE FEATURE WAS TRIED AND ALSO FAILS

`s22/mreal.py` (coordinator, split-half, 8 repeats, half-pools of 250, 126/126) refuted the
"m*-headroom is min-of-K noise over a flat correlated ladder" hypothesis: selecting m on one random
half of each target's pool and evaluating it on the disjoint other half transfers **65% of the
apparent headroom** (`-0.239 [-0.309,-0.173]`, 2.4x its own MDE, 79W/47L). **The per-target optimal
m is a real, stable property of the target — not an artefact of finite-pool sampling.** This is an
ORACLE result (the half-A selection uses half-A's own RMSD to native) — it establishes the signal
exists and is stable, not that it is native-free-reachable. Read beside the negatives already in
this file (Ridge on 5 arm sets; RF in-fold −0.276 reversing to +0.079 held-out; the audit lane's
length-quartile router at +0.031), the diagnosis changes from *"the ceiling may be an artefact"* to
*"the ceiling is real and every feature set tried so far is blind to it."*

**The coordinator proposed one further, mechanistically distinct feature family before treating
that blindness as structural**: the **within-window diversity profile** — `spread(m)` = mean
pairwise RMSD among the top-m candidates by score, for m in the ladder, plus its shape (ratios to
the full-pool spread, successive differences, and the score-boundary fraction at each rung). This
is a property of the retained candidate SET's own geometry — mechanistically what `m` actually
controls (how far the averaging benefit persists as the window widens) — and is a different kind of
feature than everything tried before (`rg_*`/`score_*` are target-level distogram-confidence
summaries; this is candidate-set-internal geometry). Built in `s22/mgeom.py` (126/126, complete),
17 features, and tested with the SAME held-out-fold Ridge-per-arm architecture, restricted to the
m-ladder-only arm set to match `mreal.py`'s scope, in `s22/mgeomrouter.py`.

    feature set                          achieved   vs incumbent (full pool, 3.048)          capture
    GEOMETRY ONLY (declared primary)     3.0618     +0.0135 [-0.0512,+0.0790], MDE 0.093       -3.8%
    old features (rg_*/score_*)          3.0960     +0.0477 [-0.0079,+0.1052], MDE 0.081      -13.6%
    combined                             3.0933     +0.0449 [-0.0109,+0.1032], MDE 0.082      -12.8%

**The geometry feature is measurably LESS HARMFUL than the old feature set** (its CI is the
narrowest and most centred on zero of the three, and it is the only one not clearly pointing
positive/worse) **but it still does not capture any of the 65%-transferable signal** — the point
estimate remains slightly negative-of-zero (router very slightly worse, not better) and the
falsifier (beat the incumbent past its own MDE with a CI excluding zero) does not fire.

**A nonlinear model reproduces, a third time, the exact in-fold-wins/held-out-loses pattern already
seen on the old features**, now on the geometry features specifically:

    RF in-fold (unlimited leakage)     -0.1854 [-0.2563,-0.1223], 65W/28L    -- looks like a big win
    RF held-out fold (proper CV)       +0.0599 [+0.0041,+0.1171], 34W/59L   -- significantly WORSE, sign flips

This is the same failure mode documented earlier in this file on the old feature set
(`-0.276` in-fold → `+0.079` held-out), now confirmed on a feature family designed specifically
around the mechanism `mreal.py`'s transfer test pointed to. **Four independent router
constructions now agree**: Ridge/old-features, Ridge/geometry-features, RF/old-features (reversing
sign under CV), RF/geometry-features (reversing sign under CV), plus the audit lane's length
router — **none captures the transferable m* signal `mreal.py` proved is real.**

**Disposition, stated as the coordinator invited**: this is the "one honest attempt" taken rather
than declined, and it closes negatively. The signal `mreal.py` demonstrated is real and stable
across independent halves of the same target's own pool, and is invisible to every native-free
feature family tried by both lanes — target-level distogram confidence, target-level score-
distribution shape, and now candidate-set-internal diversity geometry, under both a linear and a
tree-ensemble model, in-fold and (correctly) held-out. **The honest closing statement is that the
blindness looks structural rather than a feature-engineering shortfall** — consistent with, and now
sharpening, this project's standing finding that native-free functionals of the objective's own
score/candidate structure carry little exploitable rank information (`the objective cannot audit
itself`, S21 L25/E6) — extended here to candidate-set GEOMETRY, which is a genuinely different
channel from anything E6 tested and still returns the same verdict.

**Artefacts added**: `s22/mgeom.py` → `s22/results/mgeom.json` (126/126, complete); `s22/mgeomrouter.py`
→ `s22/results/mgeomrouter.json`. Both atomic writes, config-derived names.

## WHAT THIS LEAVES OPEN (SUPERSEDES THE PRE-`mreal.py` VERSION OF THIS SECTION)

`mreal.py` came back POSITIVE, not negative (Addendum 2): the m* signal is real and 65% transferable
across independent halves of the same pool. So Priority 1 is **not** closed by "the ceiling wasn't
real" — it is closed, more sharply, by **four independent router constructions failing to reach a
signal now proven to exist**. That is a stronger and more specific negative than the one I would
have written if `mreal.py` had also come back null, and it is the version I am reporting.

**What remains genuinely open, in order of how cheap it would be to try next:** (i) the gap between
ORACLE-transferable (65% of 0.350 Å on the m-ladder alone) and native-free-reachable (~0% across
four constructions) is now the sharpest, best-evidenced open question this workstream produced —
worth a dedicated feature search rather than the two families tried here, but I do not have a third
mechanistically-motivated candidate to propose without repeating the same shape of attempt. (ii)
Experiment 2's power limit, independent of all of the above — the MDEs there (0.06–0.15 Å) are wide
enough that a genuinely better structural weight (a real per-pair torsional-sensitivity term, not
the chain-separation proxy used here) is still untried. (iii) Sequence/ESM and secondary-structure
features remain NOT MEASURED for the router, for the same compute/time reasons declared in
`PREREG_C.md` §1, and are now a better-motivated next attempt than they were before `mreal.py`
landed, precisely because two structural/geometric feature families have both failed and a
sequence-level channel would be a genuinely different kind of signal again.
