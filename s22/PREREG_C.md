# SPRINT 22 — PRE-REGISTRATION, WORKSTREAM C
## Capturing the routing headroom, error-placement weighting, and the readout/window question

**Written before any Sprint-22 RMSD number exists on this instrument.** Frozen at first commit;
deviations are recorded as dated addenda in `agentC_FINDINGS.md`, never silently repaired.
Falsifiers registered first, in every block, with an honest prior attached.

I read `s22/BRIEF.md`, `s21/LEDGER.md` (L1–L37, bodies), `s21/CLAIMS.md`, `s12/instrument.py`,
`s21/poolgap.py`, `s21/latentsel.py`, `s21/rgsign.py`, `s21/rgcheck.py`, and — after the
coordinator's message arrived mid-review — `s22/LEDGER.md` L1–L3. This document is written after
that message, so its priors are informed by L1–L3; I say so rather than pretend otherwise, and I
still commit to falsifiers before the router's own numbers exist.

---

## 0. WHAT L1–L3 CHANGE, STATED BEFORE THE DESIGN

1. **L1**: headroom is broadly distributed (50% from top 21 targets, median per-target headroom
   0.317 Å, only 12/126 targets already-optimal at the incumbent). A router does not need to find
   a needle. `corr(headroom, incumbent RMSD) = +0.564`.
2. **L2**: the generative latent's aggregate loss (`lat_avg75 - incumbent = +0.390`, S21 F1) hides
   a monotone, sign-flipping interaction with TRUE difficulty (quartiles +1.469 → +0.448 → +0.152 →
   **−0.455**). Retrieval wins where retrieval works; the latent wins exactly where retrieval fails.
   **Consequence for my arm set: the four latent arms are promoted to first-class routing
   candidates, not a side experiment.**
3. **L3**: two native-free UNIVARIATE difficulty proxies (`score_mean` r=+0.433, `pool_spread`
   r=+0.360) reproduce the WRONG SIGN on the hardest quartile — both say "latent loses" where it
   wins. **A router that estimates difficulty and thresholds on it is closed by L3, in advance.**
   My router must not do that; it must let the model relate *features* to *each arm's own RMSD*,
   because the winning arm depends on more than a difficulty scalar (L2's mechanism is presumably
   about what KIND of information the target's fold requires, not only how hard it is).

**Working hypothesis, stated so it can fail**: a per-arm regression on a small native-free feature
set (rg_z/rg_gap plus pool geometry and score-distribution features), cross-fitted on held-out
folds, captures *some* but not most of the 0.482 Å ceiling — because `rg_z`'s partial ρ with true
difficulty is only +0.126 (per the coordinator) while its ρ with the *bimodal ordering regime*
(E7/L27/L28) was 0.34–0.38 on a DIFFERENT label (`ORACLE_rank_pct` within the latent, n≤13 panel).
Those are different targets of prediction and the transfer is not established. **I expect a
positive but small capture, well under half the ceiling**, consistent with `in-band-ordering-is-
per-target`'s "capacity is saturated by a linear model" and `in-band-signal-limited-not-sample-
limited`'s flat learning curve. I am registering this now so a small positive is not silently
reframed as a win after the fact.

---

## 1. EXPERIMENT 1 — THE ROUTER (PRIORITY 1)

### Question
What fraction of the 0.482 Å ORACLE routing ceiling does a native-free, held-out-fold-scored
router capture, against the 3.048 Å incumbent (pool top-75 coordinate average)?

### Arm set — DECLARED BEFORE FEATURES ARE COMPUTED
**World A** (window/point-cloud basis, `s21/poolgap.py`, bit-identical to the incumbent):
`avg_500, avg_150, avg_75 (INCUMBENT), avg_20, avg_5, avg_1, medoid_all, medoid_75` — 8 arms.
**World B** (rebuilt-torsion basis, `s21/latentsel.py`): `lat_argmin, lat_avg75, lat_medoid,
lat_avg75_rand` — 4 arms, promoted to first-class per L2.
**Total 12 arms.** (The coordinator's L1 states 13; I cannot reconstruct the 13th from the source
files inherited and do not invent one — reported as a discrepancy, not silently matched.)

### Features — DECLARED BEFORE ANY ROUTER IS FIT
`n` (length), `rg_disto, rg_pool_mean, rg_pool_sd, rg_gap, rg_z` (extended to all 126 targets, not
only the n≤13 panel — a scope extension of `rgsign.py`, stated as such), `pool_spread` (native-free
difficulty, `rgcheck.py`'s construction), `score_mean, score_sd, score_margin, score_iqr_over_range`
(functionals of the shipped score's OWN distribution over the K=500 pool — D's exhausted family,
carried anyway because a router combining them with `rg_*` is an untested object even though each
alone is not), `sim_mean_top75` (mean BLOSUM retrieval similarity of the score-selected top-75 —
a cheap sequence-adjacent feature). **NOT computed: ESM embeddings, ideal-geometry secondary
structure, VQE statistics** — declared out of scope for compute/time reasons, not because they are
expected to fail; recorded as NOT MEASURED, not as tested-and-null.

### Method — DECLARED BEFORE RESULTS
Per held-out fold (the pinned 5-fold structure, never touched for model selection): standardize
features on the training folds only; fit one **Ridge regression per arm** (`RMSD_arm ~ features`,
alpha chosen by nested CV within the training folds only) on the training folds; predict every
arm's RMSD on the held-out fold; route each held-out target to `argmin_arm` predicted RMSD.
Concatenate the 5 held-out predictions into one out-of-fold routed series — this is the ACHIEVED
router, and it is the only number allowed to be compared to the incumbent.

Two robustness variants run alongside, not chosen after seeing which wins:
(a) **World-A-only router** (8 arms, single clean basis) — isolates the basis-mixing disclosure.
(b) **rg-only router** (2 features: `rg_z, rg_gap`) vs the **full-feature router** — prices the
    marginal value of the extra engineered features over the one channel that has cleared every
    control.

### Operator forks — RULE 0, all six axes, named before the run
- **functional**: DECLARED the shipped Bayes-risk score wherever an arm depends on one (both
  worlds already use it). NOT TAKEN squared-distance (C2 showed near-interchangeability on
  pool-like structures but D7/C3 showed it does NOT transfer to the latent — so this is not a free
  substitution here and is not attempted).
- **basis**: DECLARED to keep World A and World B in their OWN native bases (window/point-cloud for
  A, rebuilt-torsion for B) rather than force one onto the other. NOT TAKEN a single rebuilt basis
  for everything, which would require rebuilding all pool arms — expensive, and L20/L14 already
  measured the resulting price (+0.016 avg75, +0.076 argmin) as small; disclosed, not corrected
  (correcting an argmin with a mean-shift constant is the M8 error).
- **readout**: THE ROUTER'S OWN OUTPUT — the whole experiment is about selecting a readout/arm.
  NOT TAKEN a single privileged readout; the m-ladder and consensus/latent arms are exactly
  different readouts.
- **normalisation**: DECLARED z-scoring FEATURES only (mean/sd from training folds), never the
  RMSD target — z-scoring RMSD would hide which targets the router fails on (L20's own rule).
- **null**: DECLARED the incumbent itself (`avg_75`, always in the candidate set, so the router can
  choose to do nothing) and the ORACLE ceiling (min over arms, my own 12-arm set, reported beside
  the coordinator's 13-arm 2.566 as a cross-check, not assumed equal).
- **THE LABEL**: no binarised label anywhere in this design — router training regresses continuous
  RMSD per arm, not a thresholded "wins/loses" label, precisely because BRIEF §7's sixth fork
  (L25's `contained512` disaster) was manufactured by a length-dependent threshold. A continuous
  target cannot manufacture a length predictor the same way.

### Falsifier
**Router is a negative result if** the achieved (held-out-fold) router does not beat the incumbent
past the comparison's own MDE (`2.8016×SE`) with a CI excluding zero, **or** if it beats the
incumbent in-fold but the held-out number does not reproduce it (the exact hazard BRIEF names).
Both directions are reported regardless of which one fires; a router that wins in-fold and loses
out-of-fold is written up as the negative result it is.

### Null / matched control
The **oracle ceiling** (min over arms, ORACLE, never achievable) upper-bounds the achievable gain.
The **incumbent** (always-route-to-avg_75) lower-bounds it. A **random arm assignment control**
(route to a uniformly random arm from the 12, matched in count) is also reported — it must be worse
than the incumbent given the sprint's own G2/G3 finding that some arms are individually much worse
alone (medoid_all: +0.658 vs incumbent) — this validates that the routing problem is not vacuous.

### Budget
Full 126-target panel, 5 pinned folds, closed-form Ridge fits (seconds of compute once features
are built). The expensive part is feature extraction (distogram + pool geometry per target),
already largely cached by `s12/instrument.py`.

### Promotion criterion
None from this file alone. A positive result is reported as CAPTURED-FRACTION-OF-CEILING with its
own CI and is a candidate for the sprint's headline only after the coordinator's independent review
per Rule 0 clause 2.

---

## 2. EXPERIMENT 2 — ERROR-PLACEMENT WEIGHTING, `H_D* = Σ s_ij w_ij (Δd_ij)²` (PRIORITY 2)

### Question
Does re-weighting the pairwise squared-distance objective by a NATIVE-FREE structural weight
`s_ij` — attacking discrimination INSIDE the objective, as opposed to post-hoc routing — move the
certified optimum (pool argmin) or the deployed readout (top-75 average) toward the native, on the
full 126-target panel?

### Candidate weights — DECLARED BEFORE ANY RMSD IS COMPUTED, so none is chosen by looking
All four multiply the existing precision weight `1/sd_ij²` (the shipped Bayes-risk score already
carries a distogram-uncertainty term; these are NOT replacements for it but structural ADD-ONS on
top of a Gaussian approximation to it, since the shipped score's own risk-table form is not linearly
decomposable into a re-weightable sum — declared and stated as the functional fork below):

    s_control    1                                   plain precision-weighted squared error (control)
    s_local      1 / (j - i)                          UPWEIGHTS local pairs (torsion-space locality
                                                        theorem: d_ij depends on exactly j-i-1
                                                        torsions -- fewer torsions, more localizable
                                                        per-torsion signal)
    s_global     (j - i)                               UPWEIGHTS long-range pairs (fold-level
                                                        correctness -- RMSD is a GLOBAL statistic and
                                                        local pairs can be satisfied while the global
                                                        fold is wrong)
    s_discrim    Var_pool(d_ij)                         UPWEIGHTS pairs where the retrieval POOL's own
                                                        candidates disagree -- a native-free proxy for
                                                        "this pair carries information that can tell
                                                        candidates apart", as opposed to a pair every
                                                        candidate already agrees on (e.g. i,i+1)

All four are fixed, closed-form, deployment-time-free functions of native-free quantities (chain
separation or the pool's own realised variance) — no fitting against RMSD, so no CV leakage
question arises for this experiment; the four are pre-declared and none is dropped after the run.

### Operator forks — RULE 0
- **functional**: DECLARED a Gaussian approximation to the shipped risk (`Σ (dhat_ij - d_ij)² /
  sd_ij² * s_ij`, sd from the distogram) because the shipped risk table has no closed pairwise-sum
  decomposition to re-weight. NOT TAKEN the shipped Bayes-risk table itself — the fork this
  introduces (Gaussian-approx vs table-lookup) is a NEW functional difference beyond D7/C2's
  squared-vs-Bayes fork, and is why `s_control` (weight ≡ 1) is run first and compared against the
  shipped score's own RMSD as a soundness check before any `s_ij≠1` result is trusted.
- **basis**: DECLARED the pool's own window coordinates (matches the incumbent, World A above).
- **readout**: BOTH reported, never conflated — certified optimum (pool argmin under H_D*) and
  deployed readout (top-75-by-H_D* coordinate average), labelled separately per BRIEF §7.
- **normalisation**: none beyond the existing `1/sd_ij²` precision term, which is already in the
  shipped score and is not this experiment's contribution.
- **null**: `s_control` (weight ≡ 1) IS the null for the three structural weights; if a structural
  weight does not beat `s_control` past its own MDE, it demonstrates nothing about structure.
- **THE LABEL**: none; weights are continuous functions of `(i,j)`.

### Falsifier
**Dead unless at least one of `s_local`, `s_global`, `s_discrim` beats BOTH `s_control` (the
weight-1 null) AND the shipped incumbent past its own MDE, on the DEPLOYED readout** (the certified
optimum is reported for mechanism, not for the primary — Sprint 20/21's repeated lesson that a
certified-optimum win can be a properties-of-the-move artefact, not a discrimination gain).

### Null / matched control
`s_control` as above; and the shipped score's own top-75 average (3.048) as the deployment
anchor, so a structural weight is judged against the thing it would have to replace, not against
a strawman Gaussian.

### Budget
126 targets × 4 weight schemes (+ 1 shipped-score soundness check) × (pool argmin + top-75 average)
— closed-form matrix operations over the cached K=500 pool, no retraining, seconds per target.

### Promotion
None from this file alone; a positive result here is mechanistically important (it would be the
first native-free structural weight to move the objective itself, not merely post-hoc routing) and
must be attacked before promotion, per BRIEF's closing instruction.

---

## 3. EXPERIMENT 3 — THE READOUT AND WINDOW (PRIORITY 3, SCOPED)

**Folded into Experiment 1 rather than run as a separate file**, because the arm set already spans
the readout axis (m-ladder = different window sizes; medoid vs average = different readout
operators; latent arms = different source+readout combination). The router's own choice of arm
per target IS a per-target readout/window selection. I additionally report, as a declared
secondary read of Experiment 1's output (no new fitting): the DISTRIBUTION of which `m` (or
medoid/latent) the router selects, split by the L1 headroom terciles, to see whether the router's
implicit per-target `m` tracks anything sensible. **This is descriptive, not a new falsifiable
claim**, and is labelled as such.

I am NOT attempting "a structural term that makes the selected SET better rather than the selected
member" as its own experiment this sprint — it requires a new differentiable set-level objective,
which is a materially larger build than the remaining budget after Experiments 1–2, and BRIEF's own
scope guidance (search/consensus/source are closed) argues against inventing a tenth consensus-like
mechanism under time pressure. Recorded as **OPEN, not attempted**, rather than silently dropped.

---

## 4. OPERATOR-FORK LIST SENT BEFORE THE RUN (RULE 0, CLAUSE 2)

Sent in this document, dated, before either experiment's numbers exist. Summary for review:

1. Experiment 1 arm-count discrepancy (12 vs the coordinator's 13) — disclosed, not resolved by
   invention.
2. Experiment 1 mixes two RMSD bases (World A window/point-cloud, World B rebuilt-torsion) inside
   one candidate set — disclosed, priced by inheritance from L14/L20, and a single-basis (World-A-
   only) robustness variant is run alongside so a reader can see what the mixing is worth.
3. Experiment 1 uses continuous regression targets, not a binarised label, specifically because of
   BRIEF's sixth fork (the L25 length-dependent-threshold hazard).
4. Experiment 2's functional is a Gaussian approximation to the shipped risk, not the risk table
   itself — a new fork beyond the project's existing squared-vs-Bayes one — and is soundness-gated
   against the shipped score before any weighted variant is trusted.
5. I have a directional hypothesis for Experiment 2 (`s_discrim` is my best guess, on the theory
   that "which pairs can discriminate candidates" is closer to the project's demonstrated
   discrimination bottleneck than either chain-separation weighting) — **naming it now, before the
   run, is the point of clause 2**: if `s_discrim` wins, that alignment with my prior prediction is
   exactly the pattern M2 warns about, and the fork list above exists so a reviewer without that
   stake can check it.

---

## 5. DATED ADDENDUM — 2026-09-07, BEFORE ANY DATA IN THIS FILE, RESPONDING TO THE COORDINATOR'S
## CORRECTION OF `s22/LEDGER.md` L2

The coordinator's L2 (latent wins on the hardest quartile of TRUE difficulty) partially withdrew:
re-stratifying on `pool_oracle` (still ORACLE, but independent of the incumbent's own realisation,
avoiding regression-to-the-mean from binning on the thing being explained) leaves a strong
CONTINUOUS interaction (`b = -0.625`, SE 0.065, 9.6 SE — the latent gains relative to retrieval as
targets get harder) but **the binned claim that the latent significantly wins on the hardest
quartile is WITHDRAWN** (null on the independent stratifier, `-0.155 [-0.388,+0.062]`). At mean
difficulty the latent is still **+0.405 worse**. This landed after §0–§4 above were written but
before any Experiment-1/2 code ran, so it is recorded here as a dated addendum rather than a silent
rewrite, per the file's own rule.

**Three concrete changes to the design above, made before any router number exists:**

**(a) A hard difficulty-threshold switch to the latent is not attempted.** It was never the design
(§1's method is per-arm regression, not a difficulty-then-switch rule), but the coordinator's
correction closes even the possibility that such a rule would have worked — there is no clean
region where the latent wins outright, only a continuous relative gain with no crossing point in
the range tested. Recorded so the router's design is not mistaken for a workaround of a route that
was never open.

**(b) ARM BLENDING is added to the candidate set as a genuinely new mechanism, not a repackaged
switch.** A 9.6-SE monotone interaction with no crossing point is exactly the shape in which a
WEIGHTED COMBINATION of the pool and latent averages can beat both endpoints even though neither
wins alone (the coordinator's point (a)). Implementation: superpose `lat_avg75`'s point cloud onto
`avg_75`'s frame via Kabsch (both are already coordinate-averaged point clouds of the same target,
so this is a superposition, not a basis conversion — no new basis fork), then form
`blend_f = f * avg_75_coords + (1-f) * lat_avg75_superposed`, `f ∈ {0.0, 0.25, 0.5, 0.75, 1.0}`
(`f=1` reproduces the incumbent exactly; `f=0` reproduces `lat_avg75` up to the superposition,
priced separately). These five points are added to the World-A candidate arm set as `blend_00,
blend_25, blend_50, blend_75, blend_100`, subject to the SAME per-arm regression and held-out-fold
routing as every other arm — no new modelling machinery, so no new overfitting axis. **Falsifier
addendum**: if no blend fraction's ORACLE (min over `f` per target) beats the better of the two
pure endpoints (`min(avg_75, lat_avg75)` per target) by more than that comparison's own MDE,
blending is dead and only switching-among-discrete-arms is reported.

**(c) The feature `rg_z` is now read as NOT a difficulty proxy** (its correlation with true
difficulty is +0.126 per the coordinator, far below `score_mean`'s +0.433 and `pool_spread`'s
+0.360, both of which the coordinator has now shown identify the wrong hard targets). It is
retained in the router's feature set exactly BECAUSE it is orthogonal to difficulty — the working
hypothesis is that it may help the regression discriminate WHICH arm wins independent of how hard
the target is, which difficulty alone cannot do (L3). This is not a new falsifier; it sharpens why
`rg_z` is a feature and not a difficulty-based routing rule, and the rg-only-vs-full-feature
ablation in §1 will show whether it is pulling weight in the fitted model or is dead weight.

**No other part of §0–§4 changes.** The arm-count discrepancy, the basis-mixing disclosure, the
continuous (non-binarised) regression target, and Experiment 2 are all unaffected by this
correction — it is entirely about the mechanism by which the latent might help within the router,
not about routing methodology in general.
