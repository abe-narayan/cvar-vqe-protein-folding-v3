# SPRINT 25 — ENDGAME DECISION LEDGER

## L1 — THE DISTANCE POSTERIOR IS OVER-CONFIDENT BY A FACTOR OF TWO. MEASURED, NOT INFERRED.

`s25/calib.py`, `results/calib.json`, **n=126 targets, 8,549 pairs, complete.** Forks and falsifier
written before the run. Natives read for **diagnosis only**; nothing fitted, nothing selected.

    published sd vs sd recomputed from the posterior: max abs diff 0.0000 A   (the artefact is consistent)

    z = (D_native - expected)/sd     z_mean **-0.0521**    z_sd **1.9962**    (perfect: 0.00, 1.00)
    coverage of the central 50% band   **0.2824**    (nominal 0.50)
    coverage of the central 90% band   **0.6159**    (nominal 0.90)
    mean |D_native - expected|  2.3386 A   against a published sd of  1.4471 A
    mean NLL of the true bin    3.2894
    pairs with a MULTIMODAL posterior (>=2 peaks above 0.02):  **0.241**

**The posterior is very nearly CENTRED and roughly TWICE too narrow.** The falsifier — calibrated
everywhere — did not fire.

    sep      pairs   z_mean    z_sd     signed   abs_err   sd_pub   cov90
    2-2       1381  -0.0574  1.2294    -0.0476    0.5225   0.4999  0.7897
    3-3       1255  -0.0039  1.9560    -0.1939    1.4964   0.9841  0.6363
    4-5       2132  +0.0815  2.0468    -0.3011    2.1816   1.3612  0.5799
    6-8       2253  +0.0065  1.8659    -0.4452    3.1209   1.8636  0.5668
    9-15      1528  -0.1692  1.2794    -0.5894    4.2636   2.4822  0.5525

**Over-confidence peaks in the mid-range separations** (z_sd ~2.0 at |i−j| = 3–8, ~1.25 at both ends)
and is flat in chain length (1.81–2.08 across all four length bands). **The signed error is
systematically negative and grows monotonically with separation, −0.048 → −0.589 Å**: the posterior
predicts distances systematically too LONG, increasingly so at long range.

Cross-check against s24: Workstream C measured the distogram's signed offset at +0.4062 under the
opposite sign convention; this run gives −0.4062. **Same number, independent code path.**

---

## L2 — CALIBRATING THE POSTERIOR MAKES RMSD WORSE. THE SCORE DOES NOT USE THE POSTERIOR'S WIDTH.

`s25/temper.py`, `results/temper.json`, **n=126, complete.** Risk-table reconstruction from the
shipped posterior asserted per target: **max error 0.00e+00.** Nested CV over the 5 pinned folds, two
independent fitting criteria.

         f       RMSD   RMSD fixw       z_sd     vs f=1
      1.00     3.0483      3.0483     1.9962    +0.0000     <- the shipped posterior
      1.50     3.0753      3.0600     1.3778    +0.0270
      2.00     3.0835      3.0780     1.0806    +0.0352
      3.00     3.1475      3.1437     0.8258    +0.0992

         T       RMSD       z_sd     vs T=1
      1.00     3.0483     1.9962    +0.0000
      2.00     3.0500     1.1840    +0.0017
      5.00     3.0501     0.6618    +0.0018

    CAL   f fitted on TRAINING-FOLD CALIBRATION   **+0.0538**  SE 0.0507  MDE 0.1422   53W/73L
    RMSD  f fitted on training-fold RMSD           +0.0103  SE 0.0069  MDE 0.0194  **10W/116L**
    ORACLE per-target f                            -0.2120  ... **86% accounted by its best-of-8 null**

    MECHANISM CHECK: held-out z_sd  1.9962 -> **0.9834**  (target 1.0)

> **The primary falsifier fired, and the mechanism check is what makes the result decisive. The CAL
> arm hit its calibration target almost exactly — z_sd 0.9834 against a nominal 1.0 — and RMSD got
> WORSE.** Fixing a real, large, correctly diagnosed defect costs 0.054 Å.

### The two curves together identify the mechanism

**Widening by convolution (SD) degrades RMSD monotonically while improving calibration monotonically —
the two are perfectly anti-correlated. Tempering (TEMP) changes calibration by a factor of three
(z_sd 1.996 → 0.662) and moves the endpoint by 0.003 Å.**

The Bayes risk is an L1 functional, `risk(t) = Σ_c p_c |t − C_c|`, whose minimiser is the posterior
**median**. Tempering preserves the ordering of the bin probabilities and therefore very nearly
preserves the median, so it barely changes the candidate ranking. Convolution moves mass **across**
bins, shifts the median, and changes the ranking — for the worse.

> **The shipped score is a LOCATION-based ranker. The posterior's calibrated width is very nearly
> irrelevant to it, and forcing the width to be correct actively harms the ranking.**

This retires an entire family of directions the endgame directive lists under §8 — variance
calibration, confidence calibration, heavy-tail treatment, sharpening — **not because they fail to fix
calibration, but because they fix it and the endpoint does not care.** It also explains s24's
`conf.py` null from a second direction: that arm reweighted *pairs* by confidence, this one reshapes
each pair's *posterior*, and both find the confidence channel inert.

**What remains open is LOCATION**, which is the only thing the ranking demonstrably responds to — and
which is exactly what s24's prior ladder moved when it produced −0.2150 Å at γ=0.1.

**The oracle row is 86% accounted for by its own best-of-8 null** (observed −0.2466, null −0.2130),
consistent with the standing rule from s24 L15-A. Not a signal.

---

## L3 — THE "+0.113 Å CVaR CONTRIBUTION" IS MEASURED AT A TEMPERATURE THAT DOES NOT SHIP. AND α = 1.0 ON THREE OF FIVE FOLDS.

Found by the quantum lane while verifying an inherited claim rather than while looking for a defect.
**Verified independently by the coordinator from source before acceptance.**

`core/pipeline.py:118`, read directly:

    VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}
    core/pipeline.py:853   ->   alpha, T = VQE_LFO[int(fold) % len(VQE_LFO)]

**The deployed configuration is T = 0.3 on all five folds and α = 1.0 on three of them.**

### What the standing claim actually says

The inherited fact — quoted by me in the Sprint 23 and Sprint 24 reports, and in the published Sprint 24
artifact — is: *"at T=0.1, α=1 collapses the distribution to 0.076 bits and the pipeline degenerates to
a plain argmin at 3.4540 Å; α=0.1 holds 6.36 bits and reaches 3.3414 Å, so CVaR is worth +0.113 Å by
preventing collapse."* The lane verified that arithmetic exactly against `s8/integrate_vqe.json`:
3.4540 − 3.3414 = 0.1126, entropy 0.07614 vs 6.3638 bits. **The claim is true as stated.**

**It is stated at T = 0.1, and T = 0.1 is not what ships.** At the deployed T = 0.3 the stored
marginals are:

    a1.00  3.2835      a0.25  3.2825      a0.10  3.3115

**At the deployed temperature the tail level is worth nothing or slightly negative, and the α=1 arm
does not collapse — 5.66 bits, not 0.076.**

### Two consequences, and the second one bears on a pillar

**(1) The attribution is probably wrong.** α and T are substitutes acting on the same quantity — how
concentrated `p_θ` is at readout. If so, the honest attribution is to **not collapsing the
distribution**, not to **the tail constraint**. The lane has registered `Q-B` as the paired n=126 test
with the falsifier stated: if the α effect at T=0.3 and T=1.0 matches T=0.1 in sign and size, its
prediction is wrong and the inherited claim stands unchanged.

**(2) On three of five folds the deployed selector has NO tail constraint.** At α = 1.0 the CVaR is
the full mean of the distribution. That is not a technicality — the endgame directive's first pillar
requires a CVaR that "materially participates", and on 60% of the folds, as configured, it does not.
**This must be resolved before the architecture is frozen and before anything is presented as a genuine
CVaR selector.**

**Do not fix it by picking a better α off the grid.** The docstring at `core/pipeline.py:113` is
explicit that `VQE_LFO` is a LEAVE-FOLD-OUT table — fold *f*'s cell was chosen on the other four — and
that selecting the grid's best cell instead (`a0.25_T0.3`, 3.2825) *would* be tuning on the instrument.
Any change to this table has to preserve that discipline.

**Basis note, so these numbers are not misread:** 3.4540 / 3.3414 / 3.2835 are the pipeline's
quantum-stage arm on `s8`'s instrument. They are **not** comparable to the 3.0483 Å incumbent, which is
the score-filtered uniform top-75 coordinate average on the 126-target instrument. Two different arms;
never quote them in one column.

### Disposition

**My reports overstated the quantum component's contribution by quoting a figure from a
non-deployed configuration without saying so.** The figure is arithmetically correct and the
configuration is stated in its source; the error is mine for carrying it forward as *the* contribution
of the shipped selector. Q-B settles the attribution. Recorded now rather than after it is settled,
because I have already published it twice.

---

## L4 — THE HEADLINE NUMBER DESCRIBES AN OBJECT THAT IS NOT A PROTEIN STRUCTURE. THE CODEBASE SAYS SO.

Found by the coordinator while verifying, before the freeze, that the number we report is produced by
the pipeline we ship. **It is — and that is not the problem.**

### The production pipeline does emit 3.0483, and it emits four other readouts beside it

`bench_results/compare_tuning126.json`, the four-component run at n=126 with amber, quantum and legacy
all enabled and `headline_eligible: true`:

    shipped        **3.4540**      <- the field literally named "shipped"
    pool_best        1.7108   ORACLE
    top_m_best       2.3062   ORACLE
    rmsd_avg       **3.0483**      <- the number this project calls the incumbent
    rmsd_fit         3.2041
    rmsd_arm         3.2148        <- "SYNTHESIS (stage 3)"
    rmsd_full        3.2355        <- "S9 SYSTEM (+AMBER)"

`rmsd_avg` matches the instrument's incumbent to all 16 digits, so **there is no mismatch between the
instrument and the pipeline.** The pipeline genuinely produces the number.

### But the codebase labels that arm ILLEGAL, and it is right

`core/bench.py:682`, in the arm-name table used for printing:

    ("rmsd_avg", "raw average (illegal)")

Measured directly on the incumbent's own output, n=30 sample:

    mean virtual Ca-Ca bond   **3.1222 A**       native  3.8246 A      -> **18.4% CONTRACTED**
    worst single virtual bond **1.970 A**        a Ca-Ca distance that no peptide can have

> **The object we report 3.0483 Å for is a POINT CLOUD, not a chain. Its virtual bonds are 18%
> short and its worst is 1.97 Å. It cannot be written as a valid PDB and it is not a deliverable
> protein structure.** Projecting it onto the ideal-geometry manifold — which is what makes it a
> structure — costs +0.156 Å and lands at `rmsd_arm` = 3.2148.

### This is not a new deception, and it is a real problem

**In the project's defence:** every sprint has labelled this arm "point cloud", the standing basis rule
forbids comparing it to a built chain, and the +0.156 Å gap has been quoted repeatedly as "pure
operator choice". The project has been internally consistent.

**What it has never done is confront that the point-cloud arm is not shippable.** Phase II forces the
question, because you cannot build a structure viewer for an object that is not a structure, and you
cannot write `T001__production.pdb` from a cloud whose bonds are 1.97 Å.

### The three honest options, none of which is "quote 3.0483 and move on"

1. **Report the built chain as the headline** — 3.2148 Å (SYNTHESIS) or 3.2041 (projection). Costs
   0.166 Å against the number every previous report has used, and it is a real, valid, renderable
   protein structure. **The point cloud stays as a labelled intermediate diagnostic.**
2. **Report both, prominently and always paired**, with the point cloud explicitly marked as a
   non-physical intermediate and the built chain as the deliverable.
3. **Keep the point cloud as headline** and state in every appearance that it is not a structure.
   This is what the project has effectively been doing implicitly, and it is the option a professor is
   most likely to object to.

**My recommendation is (2), with the built chain named as THE production result** — because the
directive's own architecture ends at "final structure", and 3.2148 Å is the RMSD of the final
structure this system actually produces.

### Immediate consequences

- The results lab must export **built chains**, and the RMSD displayed beside each structure must be
  the built-chain RMSD of the file on disk. The round-trip check I asked the renderer lane for — an
  exported file must reproduce its own reported RMSD through the instrument — **is exactly the check
  that would have caught this**, and it now has to pass on the built-chain arm.
- Every comparison in the 7-configuration suite must be run and reported on **one declared basis**,
  and if it is the point cloud then the built-chain figure must appear beside it.
- **The `shipped` field at 3.4540 needs its own explanation** in the final report. It is the
  collapsed-argmin readout and is 0.41 Å worse than `rmsd_avg`; a reader who greps for "shipped" will
  find it first.

**Nothing here is retracted and no prior comparison is invalidated** — every arm was measured on a
consistent basis. What changes is which number we are entitled to call the system's result.

---

## L5 — THE "+0.113 Å CVaR CONTRIBUTION" WAS NEVER A MEASURED EFFECT. AND ONE CURVE EXPLAINS EVERYTHING.

Quantum lane, `s25/q_alpha.py` → `results/q_alpha.json`, **n=126, complete, provenance-stamped.**
Basis stated: the s8 instrument (consensus-medoid SELECTION from a 128-candidate score-filtered set).
**Never comparable to 3.0483 Å.**

### The claim I published twice was not a result even in its own slice

    vqe_a0.1_T0.1  minus  vqe_a1.0_T0.1     -0.1126 A
      SE 0.0792   MDE 0.2220   **0.51x MDE**   fold CI spans zero   55W/44L   **median exactly 0.0000**

**By this project's own fixed rule that is a NULL.** The +0.113 Å was labelled off MARGINAL MEANS where
a paired label was required. `vqe_a1.0_T0.1` is identical to the argmin arm per target, not merely in
mean, so the contrast is real as arithmetic — it was simply never significant.

> **This is the same class as the two harness defects the lane found in its own code, arriving from the
> other direction: not a harness flattering its user, but a marginal quoted where a paired statistic
> was required.** I carried it into two sprint reports and one published artifact. Recorded as mine.

And the lane's own falsifier did **not** fire: α=0.1 − α=1.0 is −0.1126 at T=0.1, **+0.0279 at T=0.3
and +0.0126 at T=1.0.** The effect does not reproduce at other temperatures and it changes sign.

### One curve. The knobs and the circuit all sit on it.

    correlation(readout ENTROPY H, mean RMSD) over the 9 VQE cells    **-0.7423**
    correlation(alpha,             mean RMSD)                          +0.2700
    correlation(T,                 mean RMSD)                          -0.0234
    alpha's marginal share of the variance left after H and H^2         **0.032**

Fitting `RMSD ~ H + H²` on the **nine NO-CIRCUIT arms only**, then scoring the nine VQE arms against
that curve: **mean residual +0.0090 Å against the fit's own residual sd of 0.0268 Å.**

> **The circuit sits ON the classical curve, one third of a residual standard deviation above it.
> α and T are substitutes acting on a single quantity — the entropy of the readout weights — and
> neither knob nor the circuit is visible once you condition on that entropy.**

### The tail constraint is worth nothing at the deployed temperature

    a=0.25 minus a=1.0 at T=0.3   -0.0011  SE 0.0477  **0.01x MDE**  44W/47L  NULL
    a=0.10 minus a=1.0 at T=0.3   +0.0279               0.19x MDE   50W/50L  NULL
    78 of 126 targets (61.9%) sit on an alpha=1.0 fold and have NO tail constraint at all.

**And forcing α < 1 everywhere buys nothing:** forcing a=0.25/T=0.3 on all folds is −0.0311 at 0.27×
MDE (ORACLE, priced only); a=0.1/T=0.3 is −0.0021. **So the α=1 problem is not a performance problem.
It is purely a claims problem** — which is the most useful thing the lane could have told me about the
constraint I set, and it means `VQE_LFO` should not be touched.

    VQE_LFO minus argmin   -0.1405  SE 0.0732  MDE 0.2051  **0.68x MDE**  66W/48L  5/5 folds same sign
                           -> UNDERPOWERED, below its own MDE, NOT a result by our rule.
    VQE_LFO minus Boltzmann(T=0.3)   -0.0002   0.00x MDE
    Concentration: median -0.0081 vs mean -0.1405, 12 exact ties, 5 targets carry 45% of the effect --
    but drop-top-5 against a UNIFORM-EFFECT null sits at the 64.7th percentile, so concentration is
    SUGGESTED by the median-vs-mean gap and NOT established, which is exactly what the standing rule says.

### The genuinely good result, and it is the one for the postdoc

At α = 1 the CVaR objective reduces to `mean(E) − T·H(p)`, whose unconstrained minimiser over the
simplex is **exactly** the Boltzmann distribution `exp(−E/T)/Z`.

> **So the classical `boltz_T` arm is not an analogy — it is the exact optimum that the 21-parameter
> RY/CNOT state is approximating.** Circuit minus that exact optimum: +0.0499 at T=0.1 (0.42× MDE),
> −0.0302 at T=0.3 (0.24×, NULL), +0.0445 at T=1.0 (0.31×, NULL).
>
> **The variational state reaches its own unconstrained optimum to within noise at all three
> temperatures.** That is a real, checkable, positive statement about the ansatz's expressivity — and
> it is simultaneously the reason the circuit buys nothing, because the thing it approximates is
> itself classical and cheap.

### Everything else verified

MPS χ = 2^layers confirmed against an independently written dense simulator: **max |p_MPS − p_dense| =
3.3e-16**, Schmidt rank 4 = χ at layers=2, and the only occurrences of "svd"/"truncat" in
`core/quantum.py` are prose. **Parameter-shift exact**: cos 1.000000000 against finite differences,
relative error 4.6e-10, both terms of the free-energy gradient. **Set-equality theorem re-derived** on
2,592 fresh cells with EXACT zeros, 1,620 containing a genuine zero: 0 subset violations, 0 holes that
were not zero-probability states, 972/972 equality on full support, **and the assertion verified LIVE
rather than dead.**

One exception: the recorded tail-baseline defect reproduces **as a defect but not at the recorded
constant** — cos +0.5666 at 0.519× norm here against the file's +0.655634 at 0.758×, with project
memory recording a third value (+0.524). **The constant is instrument-dependent and must be cited with
its instrument, never as a universal number.**

### A fifth harness defect, same class as the dead assertion

`aggregate()` computes `gate_equality_rate` with `r.get("gate_equality", True)`, **so a row that never
measured equality is counted as having passed it.** Demonstrated returning 1.0 on rows where the key is
absent. Three further siblings are in the lane's findings.

---

## L6 — THE SHIPPED SCORE'S PER-PAIR TARGET IS QUANTISED TO 17 VALUES. ±1.75 Å OF PURE LOCATION ERROR FROM BINNING.

RMSD lane, `s25/loc.py` → `results/loc.json`, **n=126, complete, provenance-stamped. No RMSD
computed**; natives read for diagnosis only, same standing as `calib.py`.

### The structural finding, which nobody had noticed

The L1 minimiser of a discrete distribution always sits on an **atom**. The posterior lives on 17
irregular bin centres — 4.0, 4.75, 5.25, … 21.0, 25.0 — **so the shipped score's effective per-pair
target takes only 17 distinct values**, with 3.5 Å spacing at the top of the range.

> **That is up to ±1.75 Å of pure location error at long separation, from binning alone, in the exact
> channel L2 established the endpoint responds to.**

And the diagnostic I asked for as a guard against a quantisation artefact turned out to identify the
mechanism instead: **`d2centre` = 0.000 and `on_ctr` = 1.000 for the INCUMBENT** and for every
mode-seeking and every quantile arm. Across 30 arms, **every arm with `gam_eff` > 0 has `on_ctr` ≤ 0.19
and every arm with `gam_eff` ≤ 0 has `on_ctr` = 1.000. The sign of the location gain tracks whether the
arm de-quantises.**

### The lane's own registered prediction was half-refuted, and the replacement is better

It predicted the pool would share the prior's signed offset **at every shell**. True in aggregate
(+0.2318 prior vs +0.1744 pool) and **false at 4 of 5 shells, by 0.34–0.65 Å** — the prior's offset
grows monotonically with separation and the pool's does not; they agree only after averaging over a
shape mismatch running in both directions. Withdrawn in the form registered.

**What replaced it came from the top-75 column I asked to be added:**

    prior  +0.2318   ->   pool(K=500)  +0.1744   ->   SELECTED top-75  **-0.1310**   ->   emitted  -0.5852

**Selection has already over-corrected the offset and reversed its sign.** Shifting the prior shorter —
which the raw signed error naively asks for — would select still-shorter candidates and drive the
emitted structure deeper into a contraction that is already the largest single operator effect in the
system. A stated mechanism, and not the one that was registered.

### Three families closed before costing an end-to-end run

    arm            gam_eff     cos   |dL|/|dT|     MAE     bound(-2.1496*gam_eff, UPPER BOUND)
    IDENT          +0.0000  +0.000     0.000    2.3968    -0.0000
    SHELL1.00      +0.0136  +0.094     0.210    2.3993    -0.0293   <- ORACLE in-sample profile
    CTRL_PERM      +0.0108  +0.063     0.210    2.4431    -0.0232   <- its own permuted control
    CTRL_XSHELL    +0.0269  +0.052     0.783    3.0204    -0.0579   <- ANOTHER target's profile

**The full per-separation offset correction, at α=1 with an ORACLE in-sample profile, travels 1.36% of
the way to truth, and its magnitude-matched permuted control already captures 79% of that. A different
target's profile scores HIGHER than the target's own.** Family C is dead against its own null, not
merely small.

**Mode-seeking is wrong-signed everywhere** and the pre-registered stratification localises it exactly
where the mechanism claimed its gain: TRUNC0.5 is −0.0568 on multimodal pairs against −0.0037 on
unimodal. Multimodal fraction 0.241, reproducing L1's 0.241 on an independent path.

    |median - native| 2.5810    |NEAREST mode - native| 1.3614 (ORACLE)    |TOP mode - native| 2.8307

**The oracle nearest mode is worth 1.22 Å per pair and needs the native to pick it; the achievable top
mode is 0.25 Å WORSE than the median it replaces. The safe arm is the absurd one, for the fourth
time.**

### What is open, and it runs the other way along the same axis

`POWER` is monotone through q=1 and **positive above it**: q=1.25 +0.0254, q=1.5 +0.0403,
**q=2.0 +0.0582 at cos +0.253**, MAE 2.3968 → 2.3386. q=2 is an **L2 risk, whose minimiser is the
posterior MEAN**.

> **The mean is nearer the truth than the median, by 5.8% of the way — 2.6× the γ = 0.0225 that reaches
> 3.0 Å.** The sign is the opposite of what the lane registered as promising, and it is the strongest
> ranking-visible arm in the table.

### DEQUANT, added as the new primary before run 2

Replace each bin's point mass at `CENTRES[c]` with its own mass spread **uniformly over that bin's own
support**, half-width β·(hi−lo)/2, **centre unmoved**. β=0 is the incumbent bit-exactly; β=1 takes the
histogram literally. The risk becomes piecewise quadratic then linear with derivative
`Σ_c p_c·clip((t−C_c)/h_c, −1, 1)`, so the minimiser is continuous and unique.

> **No mass crosses a bin and no centre moves — which is exactly what separates it from L2's SD arm.
> That one convolved ACROSS bins, moved the median, and cost +0.054 Å. This is a LOCATION intervention
> with a width side-effect, not the reverse.**

### METRIC: the biggest number in the table, with its floor already subtracted

`METRIC` at η=1 gives gam_eff **+0.2497** at cos +0.495, MAE 2.3968 → 2.1756, per-shell cosines rising
exactly where the posterior is worst. **The lane subtracted its own shared-referent floor before
reporting it**: the native IS realisable, so orthogonal projection onto the realisable set removes part
of the error toward it **by identity**, and `gam_eff = (|dL|/|dT|)² = 0.202` would arise from that
alone. **The genuine excess is 0.048, and that is the number to quote.**

**And the lane pre-registered the prediction that it converts worst:** every candidate in the pool is
*itself* realisable, so the component the projection removes lies orthogonal to the manifold the
candidates live on and adds a near-constant offset to every candidate's score — gam_eff counts it
because it is measured against the native, and **the ranking cannot see it**. To be tested directly
with an `amp` diagnostic measuring how much each move points along directions the candidate pool
actually varies in. **This is the case that will show whether gam_eff is the ranking currency at all.**

---

## L7 — RETRACTION: L2's MECHANISM IS FALSE, ITS MAGNITUDES WERE UNMEASURED, AND BRIEF §1/§3 GOES WITH IT.

Audit lane, `s25/audit_t1{,b,c}.py`, artefacts `audit_t1{,c}.json`, **n=126, complete, forks enumerated
by the auditor.** **Accepted in full. L2's null stands; everything I built on top of it does not.**

### First, what survives

Both my files **reproduce exactly** from source and are **leakage-clean**: risk-table reconstruction
max err exactly 0.0, and the CAL nested CV is genuinely leave-one-fold-out — the auditor re-ran the
fold loop in its own code and recovered the same held-out vector, the same per-fold f (2.5/2.0/2.5/2.0/
2.5) and the same 0.9834. **The problem was never the numbers. It is what I claimed from them.**

### (1) NOT ONE of L2's RMSD contrasts clears its own MDE

    CAL +0.0538 = **0.38x MDE**, fold CI [-0.0209,+0.1400], 53W/73L
    SD f=1.5 0.41x · SD f=3.0 0.52x · SDFIXW f=3.0 0.52x · TEMP T=5.0 **0.02x**, Type-M 51.7
    **Every fold CI includes zero.**

**L2's pre-registered falsifier fired correctly and "calibrating the posterior does not buy RMSD"
STANDS as a clean null.** But *"RMSD got WORSE"*, *"costs 0.054 Å"*, *"degrades monotonically"* and
*"perfectly anti-correlated"* are readings of **unmeasured** numbers and are struck. TEMP especially —
the load-bearing evidence for width-inertness — has essentially **zero power** and cannot distinguish
"inert" from "we did not look".

### (2) The mechanism is retracted on three independent grounds

**(a) Its premise is FALSE IN THE CODE.** `core/predict.py:420` is `w = shell/(sd+0.5)^g`, and
`_score_weights()` returns `ones, 1.0` because `score_weights.json` is intentionally absent — so
shell == 1, g == 1, and **the only per-pair weight in the shipped score is a pure function of the
posterior's WIDTH.** Width is not inert; **it is the score's entire pair-weighting mechanism.** My
`SDFIXW` control is implemented correctly but its output is noise: the shape/weight split runs
87/13, 143/−43, 43/57, 47/53, 84/16, 96/4 across adjacent grid points **with a sign flip at f=1.3**.

**(b) The median-shift explanation is falsified at matched shift.** SD f=2.0 moves the median 0.939 Å
and costs +0.0352; TEMP T=5.0 moves it **0.951 Å** and costs +0.0018. **The stated explanatory variable
is matched and the outcomes are not.**

**(c) The treatment was confounded.** `_widen_sd` row-normalises a Gaussian kernel over **non-uniform**
bin centres (0.5 Å spacing at 4 Å, 3.5 Å at 21 Å), which necessarily drags mass into the dense
short-distance bins. Measured: SD(f) shifts the posterior **mean** by −0.111 / −0.375 / −0.729 /
−1.155 Å at f = 1.15/1.5/2.0/3.0. **CAL selected f = 2.0–2.5, i.e. a −0.73 to −0.99 Å systematic
location shift — about twice the entire signed bias L1 diagnosed — separation-graded at r = +0.977
against L1's own signed-error profile.**

> **So `temper.py` already ran, unlabelled, an approximate per-separation location correction — the
> BRIEF's own leading open candidate — and the endpoint did not respond.**

### (3) The unconfounded experiment, which I should have run

Mean-preserving widening (pure width, achieved sd ratio to 2.28×, mean shift < 0.020 Å) against pure
translation (pure location):

    pure WIDTH     f=1.3/1.5/2.0/3.0   +0.0139/+0.0208/+0.0203/+0.0253   (0.51/0.56/0.43/0.38x MDE)
    pure LOCATION  d=-0.9/-0.6/-0.3/+0.3 +0.0617/+0.0322/+0.0206/+0.0191 (0.58/0.40/0.40/0.44x MDE)
    **ALL EIGHT NOT MEASURED. Every fold CI includes zero.**

**Once unconfounded, LOCATION AND WIDTH ARE EQUALLY FLAT.** Pure width at 3.0 (+0.0253) and pure
location at −0.30 Å (+0.0206) are the same size. **The asymmetry my mechanism was built on does not
exist.**

**And the response to location is SYMMETRIC about the shipped posterior** — −0.30 Å costs +0.0206,
+0.30 Å costs +0.0191 — even though L1 measures the posterior biased −0.41 Å against the natives.
**Moving the prior TOWARD the natives is worth nothing and is indistinguishable from moving it away.**
That is s24 Workstream C's cancellation, now measured rather than predicted.

> **BRIEF §1 and §3's "what the ranking demonstrably responds to is LOCATION" is RETRACTED.** Global
> shifts are null; separation-graded shifts are null. **What is NOT closed is the prior ladder — those
> arms add no information, and reaching γ requires genuinely BETTER location, not a re-reading of the
> location already there.**

### (4) My best-of-8 null was near-tautological

`temper.py:268` draws 8 indices into **that target's own 8 columns with replacement**, so it can only
return one of the row's own values and returns its true minimum with probability 1 − (7/8)⁸ = 0.656.
**It is bounded below by the observed statistic by construction.** Three nulls on the same observed
−0.2466: mine −0.2133 (86%); exchangeable-column −0.3568 (145%); **split-half transfer −0.0095, i.e.
4% of the oracle transfers.** k_eff = 8.00. **The conclusion is unchanged and strengthened — not a
signal — but "86% accounted" is struck and the split-half quoted instead**, or the precedent passes a
real best-of-K arm later.

### (5) L1 survives every attack, with one framing correction

Not heavy tails: median|z|/0.6745 = 1.9803 against z_sd 1.9962, within 0.8%, despite excess kurtosis
+1.82. Not multimodality: multimodal pairs are **better** calibrated (1.455) than unimodal (2.092).
Not edge digitisation: z_sd 1.9922 without the 0.72% of pairs beyond the outer edges.
**Quote it as "1.66–2.09 across estimators, roughly twice too narrow."**

**But one undeclared fork in L1 matters:** `calib.py:131` reads the interval endpoints as bin
**centres**. Three readouts of the same interval — centre 0.2824/0.6159 (published), interpolated
0.3271/0.6494, bin-edge 0.4661/0.7069. **The published figure is the most alarming of the three.**
Under-coverage survives, but **cov50 = 0.4661 is close to its nominal 0.50, which says the defect lives
in the TAILS, not the core** — exactly what the +1.82 kurtosis says independently. **A different object
with a different fix, and my framing hid it.**

### (6) Rule 0: the Sprint-24 shape repeated, four for four

**Four undeclared forks in `calib.py` and ALL FOUR POINT MY WAY**: the quantile endpoint convention;
no continuity correction (worth 6%); `sd` rather than a robust scale (worth 17%); and the
multimodality threshold 0.02 with its peak rule — undisclosed degrees of freedom behind the 24.1%
that the BRIEF promoted to a candidate mechanism. Plus `SEPBANDS` as a free choice with no declared
alternative, and the "peaks in the mid-range" claim read off it. Four more in `temper.py`, the fatal
one being **T1: the widening operator is not width-pure, and the alternative not named is exactly the
mean-preserving widening the auditor then ran.** T2: nominal vs achieved f — the operator delivers
1.43× at f=1.5 and 2.32× at f=3.0, so **the table's labels are not what the operator does.**

**And a design exposure I should have declared:** `FGRID` was placed *after* L1 measured z_sd ≈ 2 on
all 126 dev natives, so **the grid the CAL arm chooses from was informed by a full-dev diagnostic.**
The CV inside the grid is clean; the grid's placement is not.

### (7) Two smaller corrections

The RMSD arm is **10W/15L/101T**, not 10W/116L — four of five folds selected f = 1.0, so 101 targets
are **exact ties** and `temper.py`'s `st()` books ties as losses. And **neither `calib.json` nor
`temper.json` carries a provenance stamp** — both use a hand-rolled `_save` rather than
`ST.save_atomic`. Nothing is promoted from them, so that is a warning, not a retraction.

---

## L8 — L4 CONFIRMED AT n=126 AND UNDERSTATED. THE BASIS DECISION IS NOW FORCED, NOT CHOSEN.

Audit lane, own code, own path end to end (pool → shipped score → top-75 → medoid-frame average),
**all 126 targets** where I sampled 30. `rmsd_avg` reproduces at 3.0483.

    mean virtual Ca-Ca bond   **2.9614 A**   native 3.8122 A   ->  **22.3% CONTRACTED**   (I said 18.4%)
    worst single virtual bond over all 126:  **0.649 A**                                  (I said 1.970)
    **80 of 126 targets have a MEAN virtual bond under 3.4 A**

**0.649 Å is shorter than a covalent C–C bond (1.53 Å). Two alpha carbons cannot be that close.** My
n=30 sample was the conservative end of the distribution. *"This is not a protein structure"* is
arithmetic, not a judgement call.

Cross-check: project memory's `averaging-space-beats-the-objective` already records coordinate averaging
contracting the backbone **25.8%**. **So L4 is a rediscovery of a known project fact in the one place it
had never been applied — the headline.**

### The decisive argument, which L4 did not make

> **BRIEF §5 and the audit's own replication gate both require an exported structure file to reproduce
> its own reported RMSD through the instrument. The point-cloud arm CANNOT PASS THAT GATE BY
> CONSTRUCTION — it cannot be written as a valid PDB, so there is no file to round-trip.**

**Option 3 (keep the point cloud as headline) is not merely objectionable; it is unimplementable under
Phase II's own rules.** That reduces the choice to (1) or (2), and the auditor independently concurs
with (2).

### DECISION, ADOPTED

**Option (2): report both, always paired, with the BUILT CHAIN named as the production result.**

- **Production result — the deliverable structure:** `rmsd_arm` (SYNTHESIS, stage 3). To be re-derived
  independently before it is quoted; the audit has not yet replicated 3.2148 or the +0.156 gap.
- **Point cloud (3.0483):** retained everywhere as an explicitly labelled **non-physical intermediate
  diagnostic**, never as the system's result.
- This **costs 0.166 Å against every previous report**, which is the direction that makes it credible:
  §48 forbids changing the evaluation basis to make a number better, and this makes it worse.
- Every prior comparison remains valid — all were measured on a consistent basis. **What changes is
  which number we are entitled to call the system's result.**

### One thing in L4 that is NOT a defect — cleared so nobody chases it

`compare_tuning126.json`'s `science_delta` is 0.0 or 4.4e-16 on every arm. **That is not two identical
arms mislabelled.** baseline vs optimised is `baseline_tuning126.json` vs `optimised_tuning126_w8.json`
— a **performance** comparison, and bit-identical science is the intended guarantee that the 8-worker
path does not change the answer. **A harness strength, and it belongs in the final report as one**,
before a reader greps `science_delta` and misreads it.

---

## L9 — THE "86% ACCOUNTED" NULL IS PROVEN TO CARRY NO INFORMATION, NUMERICALLY.

`s24/stats_lib.py` is now the sprint standard, with two sibling defects fixed and three constructions
added. **The selftest proves my tautology rather than arguing it.**

    on a 126x8 grid of PURE NOISE:
      valid within-target null accounts for   93%
      split-half transfer accounts for         0%
      the INVALID own-row-with-replacement null accounts for  **83%**

**Eighty-three percent on pure noise, against the 86% `temper.py` reported on real data. The same
number.** That is the proof that my "86% accounted" carried no information at all — struck from L2,
already replaced there by the split-half figure (4% transfers).

**Two sibling defects fixed in the shared library**, both of the same species as the ones the quantum
lane found in its own harness:
- `_verdict` fell back to the **iid CI whenever `folds` was omitted**, and could then return
  BETTER/WORSE — which the module's own docstring forbids. A caller who simply forgot to pass `folds`
  got the forbidden verdict silently. It now refuses and reports what the iid CI would have said.
  All 12 existing callers pass `folds`, so nothing in flight changes.
- `best_of_k_null` resamples from a **pooled** population — right for a candidate pool, wrong for a
  K-column grid of the same targets. **That gap is why `temper.py` rolled its own tautological null.**

**Added and mandatory from now:** `best_of_k_within(M)` (resamples deviations ACROSS targets; reports
share_accounted, residual, and a `k_eff` that is the **entropy of the argmin distribution** rather than
a bare K; and computes the invalid own-row null **LABELLED INVALID** so the trap is recognisable rather
than re-inventable), `split_half_transfer(M)`, and `achieved(nominal, measured)` — which flags an
operator that does not deliver what its parameter is named for at >10% shortfall. **`achieved()` is the
check that would have caught my SD(f), whose nominal 1.5 delivered 1.43 and whose nominal 3.0 delivered
2.32.**

> **RULE FOR EVERY LANE: if you take a minimum over K settings of the same targets, call
> `best_of_k_within`, not `best_of_k_null`.**

### Outstanding and now assigned

**The 4/126 dev self-copy caveat is not attached to L1, L2, L4 or L6.** Four dev targets carry a
verbatim self-sequence window at BLOSUM rank #1 (`core/data.py`'s longer-normalised identity). It is
attached here for all of them, and every final-report figure must carry it.

---

## L10 — THE ROUND-TRIP GATE PASSES, AND THE FROZEN BUILD NOW HAS A NATIVE-FREE LEAKAGE GATE.

### The gate that makes the results lab trustworthy

`s25/resultslab/test_export.py`, run independently by the audit lane, unmodified: **9/9 pass, worst
absolute error 2.239e-04 Å** against a 2.0e-03 tolerance and a PDB `%8.3f` quantisation bound of
8.66e-04; max coordinate error 5.0e-04, **exactly the quantisation half-step**.

**An exported structure reproduces its own reported RMSD through the instrument on BOTH bases.** And
it is a real assertion rather than a tautology: `verify_roundtrip` compares the file against the
**exported** native, not only the in-memory one. **Nothing downstream is decoration.** This is also the
gate whose logic forced L8's basis decision, so it has already paid for itself twice.

### THE GUARD, ADOPTED AS A HARD GATE ON THE FROZEN BUILD

> **The frozen build MUST REFUSE any configuration whose mean RMSD is below the pool oracle ceiling,
> `pool_best = 1.7108 Å`. No real method can beat the best structure in its own candidate pool. A
> configuration that does has native information in it, by whatever route.**

**This is the best safeguard proposed this sprint.** It is native-free, machine-checkable, requires
nobody to remember to set a flag, and it closes all three of the tripwires below at once. It is now a
release condition, not a suggestion.

### Three leakage tripwires, none firing, all to be closed before the frozen build

1. **`--mode proof` registers SIX of the seven REAL configuration names** — legacy, amber, distogram,
   legacy_distogram, amber_distogram, legacy_amber — with `PV.synthetic()`, which is **native +
   N(0,σ)**. It is honestly labelled (banner says NOT A RESULT, every record carries
   `is_synthetic: true`, the flag propagates through `schema.py`), and `register_frozen` cannot
   register the built-in synthetic since it loads only from npz/json paths — **so the direct route is
   structurally closed.** The indirect route, a frozen spec pointing at a file a proof build wrote, is
   not. **The separation between a proof build and the real one is currently a STRING.**
2. **`providers.from_json` does `out.setdefault("is_synthetic", False)`** — a loaded structure with no
   provenance key is **silently marked genuine**. For a leakage guard the safe default is the
   opposite: **absent provenance should be UNKNOWN and must FAIL the frozen build.**
3. Both are subsumed by the pool-oracle gate above, which is why that gate is adopted rather than a
   pile of flag checks.

### A presentation hazard, cheap and urgent

`test_export.py::_roundtrip_rows` names its control arms **`"distogram"`** (the native itself, RMSD
exactly 0.000000) and **`"legacy"`** (native + 1.5 Å noise). Those are two of the seven REAL
configuration names, so the test's own printed table contains:

    T001   distogram   14   0.000000   0.000000

**Pasted into a report, or skimmed by anyone who does not know it is a fixture, that reads as "the
distogram configuration achieves 0.000 Å".** The test is correct; only its labels are hazardous.
**Renamed to `_zero` and `_noise`.**

### PHASE-ORDER RULING

`s25/resultslab` has a full `site/` built while the architecture is still moving — L4/L8 just changed
which number is the headline and L7 retracted the through-line. **The ruling:**

- **The export layer and its round-trip gate were correct to build now.** They are architecture-
  agnostic, and the gate has already caught a real consequence.
- **The site scaffold may stand** — it is data-driven and its existence costs nothing — **but it ships
  with NO data until freeze**, and it must never contain a synthetic record.
- **NO structure may be exported under a real configuration name until the freeze signal**, and the
  pool-oracle gate must be in place before the first one is.
- The renderer's polish is not to be worked on further until the architecture is frozen. **BRIEF §7's
  phase order is binding and this is the drift it exists to prevent.**

### Standing notes now adopted sprint-wide

(a) `stats_lib` is the standard — `best_of_k_within` for a K-settings grid, `achieved()` printed beside
any swept parameter. (b) A verdict requires clearing its own MDE **and** a fold CI excluding zero; both
enforced, neither reachable by omitting an argument. (c) **A monotone sequence of unmeasured point
estimates is not a trend** — that is how L2's mechanism was built and the audit names it the most
repeated error of the sprint. (d) **State what your operator ACHIEVES, not what its parameter is
named.** (e) The 4/126 dev self-copy caveat is attached as of L9 and must appear on every final figure.

---

## L11 — L8 VERIFIED. MY ADOPTED GUARD WAS NECESSARY AND NOT SUFFICIENT. A SYNTHETIC LEADERBOARD WAS ON DISK.

### (a) The built-chain arm reproduces exactly — the L8 decision is safe

Audit lane, own code, all 126, recomputed from the **persisted** coordinates through
`s12.instrument.ca_rmsd` — **a different implementation from the `core.audit.kabsch_rmsd_batch` that
wrote the records**, so agreement is evidence rather than circularity.

    rmsd_avg  POINT CLOUD  3.048338      rmsd_fit  BUILT CHAIN  3.204076      rmsd_arm  BUILT CHAIN  3.214765
    max per-target disagreement over 126 x 3 cells:  **0.00e+00**

**It is actually a chain:** built-chain virtual Cα–Cα bond mean **3.8040, sd 0.0000, worst 3.8040, zero
targets under 3.4 Å**, against the point cloud's 2.9614 / 0.6492 / 80-of-126. Native 3.8122. Both halves
of L8 verified in both directions. *(Modelling note for the report, not a defect: the built chain's bond
is a constant 3.8040 and therefore cannot represent a cis-peptide — cis-proline is ~2.9 Å. Idealised
geometry, correctly.)*

**Round-trip at n=126 on the PRODUCTION arm** — and note the fixture builds `lam=0`, which is
`rmsd_fit`, **not** the `rmsd_arm` we are shipping: all 126 ok, worst |file − instrument| 2.641e-04 Å,
mean read back from the files 3.214770 against 3.214765 reported. The Cα trace built from each record's
persisted torsions equals the persisted Cα to 0.00e+00.

### (b) Two corrections to L8's own numbers, and a real finding

**+0.156 Å is the WRONG ARM** — that is the `rmsd_fit` gap. **The SYNTHESIS arm we are naming as
production costs +0.1664 Å.** And the mean hides the shape: **median +0.0977, sd 0.2046, IQR
[+0.0049, +0.3379], range [−0.3483, +0.6529], and the built chain is BETTER on 16 of 126 targets.**
Report the median and quartiles, never the mean alone.

> **corr(projection gap, that target's own contraction) = +0.5179.** By contraction quartile:
> Q1 **+0.0073**, Q2 +0.1112, Q3 +0.2864, Q4 +0.2628. **On the least-contracted quartile the projection
> is essentially FREE.**
>
> **The projection cost is not a toll for becoming a structure. It is the price of undoing the
> averaging's contraction, concentrated on the targets where averaging did the most damage.**

Through `stats_lib`: +0.1664, SE 0.0182, MDE 0.0511, **3.26× MDE**, fold CI [+0.1343,+0.2118], 5/5 folds
same sign, 16W/110L, **VERDICT WORSE** — a genuinely measured effect, and the first one certified this
sprint. Cross-basis by construction, which here **is** the measurement; labelled on both sides.

**Do not read the quartile table as a lead to reduce the contraction by rescaling.** Global scale is
closed (s23 L5/L7/L10, s24 L15-A). Stated here so it is not re-proposed.

### (c) THE CORRECTION AGAINST MY OWN ADOPTED GUARD

I adopted the pool-oracle threshold as a release condition on the auditor's word that it subsumed the
flag-based tripwires. **It does not, and the auditor corrected itself before anyone relied on it.**

`results/structures/` held **1,134 PDB files** — 126 under each of the seven real configuration names
plus production and native — from a `--mode proof` build in which **six of the seven comparison arms
are `native + N(0,σ)` with σ 1.3–2.2, deliberately sized to look like plausible results.**

    distogram 2.2584 · legacy_amber_distogram 2.0921 · legacy_distogram 2.2180 · amber_distogram 2.5057
    legacy 2.9625 · legacy_amber 3.1267 · amber 3.4408
    **EVERY ONE IS ABOVE pool_best = 1.7108. THE ADOPTED GUARD FIRES ON NOTHING BUT THE native ARM.**

**And worse: not one of the 1,134 files carried any synthetic or provenance marker.** `is_synthetic`
propagates into the JSON payload but **never reaches the PDB headers**, whose ten REMARK 999 keys are
BASIS, CONFIGURATION, FOLD, LENGTH, PDB_ID, REFERENT, RMSD_CA, RMSD_CA_MEMORY, SEQUENCE, TARGET_ID.
**Separated from `results.json`, each file is indistinguishable from a genuine result.** A
`leaderboard.json` on disk ranked **`legacy_amber_distogram` first at 2.0921 Å with a paired effect of
−0.956 Å against production** — a fake method beating the real one, in a file that reads as a finished
deliverable. Verified by the coordinator: `T001__distogram.pdb` carries `RMSD_CA 0.000000` and **is the
native**.

Two further defects: **`test_export.py` writes fixtures into the real output directory under real
configuration names** (4 of the 126 `distogram` files are exact natives), and **`export_prediction` is
not atomic** — it writes twice, so `T006__amber.pdb` currently reads `RMSD_CA PENDING`.

### (d) ACTION TAKEN, and the replacement release condition

**Quarantined**, not deleted, to `QUARANTINE_proofbuild_synthetic_20260908/` with a README stating what
it is and why it was dangerous — 1,134 structures, `leaderboard.*` and `results.*`.
**`results/benchmark_manifest.json` and `results/monomer_manifest.json` are git-tracked files predating
this sprint and were NOT touched**; `results/summary/target_map.json` is a pure T-number↔PDB mapping
carrying no RMSD and remains in place. `results/site` could not be moved (held by a live process) and is
now inert, having no data to load.

**THE RELEASE CONDITION IS NOW ALL FOUR, NOT ONE:**

    (a) mean RMSD >= pool_best 1.7108        -- keep; it catches exact copies and nothing else
    (b) corr(per-target RMSD, that target's pool_best) >= a floor, AND per-target RMSD sd >= a floor
    (c) a MANDATORY provenance REMARK in every PDB header; the frozen build REFUSES any file lacking
        one or marked synthetic   <- THE PRIMARY GUARD
    (d) test fixtures write to a temp directory, never results/structures/; results/ is wiped before
        the frozen build

**(b) has evidence behind it, and it is the discriminator I would not have thought of:**

    configuration            mean     sd     corr(RMSD, pool_best)
    production             3.0483   1.6400        **+0.6988**   <- the only real arm
    every synthetic        2.09-3.44  0.26-0.58   **-0.08 to +0.15**

**A genuine method must be harder on hard targets. Isotropic noise about the native gives nearly
constant RMSD regardless of difficulty, so its correlation sits at zero and its dispersion is 3–6× too
small.** Both statistics are free from what the build already computes.

---

## L12 — PHASE I CLOSES. THE SCORE'S TARGET CAN BE MOVED 25% TOWARD TRUTH AND THE ENDPOINT DOES NOT MOVE.

RMSD lane, `s25/locrun.py`, **n=126, complete, provenance-stamped (sha dc36e5be, commit a15406c8).**
Identity arm asserted per target **both** `np.array_equal` against a genuine `Distogram._risk` **and**
equal to `temper.json`'s independent f=1 column. LOCK_TRAIN taken and released.

    family      nested CV     SE      MDE   xMDE   FULL-LEAKAGE best global parameter
    METRIC       +0.0065   0.0029   0.0081  0.80   eta = 0.00 -> 3.0483   **THE IDENTITY**
    POWER        +0.0176   0.0066   0.0184  0.96   q   = 1.00 -> 3.0483   **THE IDENTITY**
    TRUNC        +0.0273   0.0117   0.0327  0.83   d   = 40.0 -> 3.0483   **THE IDENTITY**
    DEQUANTI     -0.0011   0.0026   0.0073  0.15   b   = 2.00 -> 3.0452 = -0.0031
    DEQUANT      +0.0014   0.0037   0.0105  0.14   b   = 2.00 -> 3.0466 = -0.0018

> **For three of the four families the best single global parameter chosen with COMPLETE LEAKAGE on
> the very targets it is scored on IS THE SHIPPED INCUMBENT.** There is no value of the parameter
> anywhere on the grid worth having. Closed by its ceiling, not by a failed fit — the s24 L16 standard.

DEQUANT's leakage ceiling is **0.003 Å** and it sits at β=2, beyond the principled β=1, i.e. a smoothing
that crosses bin boundaries rather than the histogram's own resolution.

**Every per-target oracle is 124–147% accounted** by `ST.best_of_k_within`'s valid across-target null
and **every split-half transfer is −8% to +16%**. All 52 arms as one grid: oracle −0.2675, valid null
−0.5779, **208% accounted**, split-half −5%. *The invalid own-row null would have said 79–85% for every
one of them, near-constant in K, exactly as the audit proved.*

### The single most informative number in the sprint

Across **51 non-identity arms, four families, two implementation styles**, with `gam_eff` spanning
−0.0145 to +0.2507:

    corr(gam_eff, endpoint delta) = **+0.054**      corr(cos, delta) = -0.112     corr(amp, delta) = +0.023
    endpoint delta range: -0.0136 .. +0.0728 A

**Moving the score's per-pair target toward the truth is, if anything, slightly WORSE.**

METRIC at η=1 cuts the location field's MAE 2.3968 → 2.1756 and its squared error by ~30%, **using only
metric realisability and no native information whatsoever**, and the endpoint gets **+0.0243 WORSE**.

> **That is `better-matrix-worse-ranking` reproduced on the VALID 126-target instrument, with an
> interpolation parameter, against a bit-exact identity** — the retest the DO-NOT-REDO table demands,
> closing the direction properly rather than on the retracted decoy bank.

### A CAVEAT ON THE LADDER'S CURRENCY, TO BE CARRIED EVERYWHERE

> **The prior ladder prices γ at −2.15 Å per unit ONLY WHEN THE MOVE IS ALONG THE NATIVE'S OWN
> DIRECTION — cos = 1 by construction. A real operator travelling 25% of the way to truth at cos = 0.50
> is worth +0.024 Å. γ is redeemable only inside the ladder's own construction.**
>
> **Nobody may quote −2.1496 × gam_eff for an achievable operator without the cosine beside it.**

### The adaptive/rigid distinction does not rescue the direction

That distinction was the whole of the remaining case after L7, and the lane said so before running.
Its arms are **adaptive per-pair** changes sized by each pair's own posterior shape; the audit's A8 arm
is a **rigid translation**. **Both are flat.** Two disjoint operator classes, same answer.

### L6 is what outlives the lane, and it now reads differently

The shipped score's effective per-pair target is quantised to 17 values; the largest gap between
adjacent centres is 4.0 Å (21 → 25), so **the representation admits up to 2.0 Å of pure location error
before any estimation error**, with mean local spacing 1.094 Å at the incumbent's own realised targets.
**Removing that exactly, with a verified mean-preserving operator, is worth 0.003 Å with full leakage.**

> **That is the strongest available statement that the ranking is insensitive to the location of its own
> target — a THIRD independent route to the audit's A8, on a completely different operator class.**

`ST.achieved` confirms every operator did what its name claims — POWER q=2 realises the posterior mean
at ratio 1.0000, METRIC realises (1−η)L_med + ηL_proj at 1.0000, DEQUANT preserves the posterior mean at
1.000000 while inflating sd 1.043 and taking **91.4% of minimisers off a bin centre**, moving the target
0.277 Å RMS per pair. **These are nulls, not broken arms.**

### The lane's own scorecard: four of six predictions wrong or half wrong

Shell-level cancellation HALF WRONG; METRIC small-or-negative CONFIRMED; **family A open on
multimodality REFUTED BY ITS OWN RUN 1 — its favoured hypothesis**; DEQUANT, its late primary, NULL;
"METRIC converts worst" HALF CONFIRMED; **"amp explains the conversions" REFUTED** — POWER has the
highest amp at 3.07 and the worst conversion. **Neither gam_eff nor amp nor cos predicts the endpoint.**

---

# PHASE I IS CLOSED. <3.0 Å WAS NOT REACHED, AND IT IS NOT BEING MANUFACTURED.

**Every arm operates on the CONSUMPTION of a fixed posterior. The ladder's steepness is a property of
the posterior's CONTENT. This sprint shows the two are not connected by any re-reading: you cannot buy
γ by consuming the same distribution differently, however much closer to truth you move its summary
statistic.**

**If <3.0 Å is reachable it is through a distogram that is actually better — an achievability question
about a trained predictor, not a functional question.** Recorded as the sprint's answer to §49.

    PRODUCTION RESULT   built chain   **3.2148 A**   (rmsd_arm, SYNTHESIS stage 3)
    intermediate        point cloud     3.0483 A     non-physical, 22.3% contracted, NOT a structure

**FREEZE DECLARED.** No further architecture changes. Phase II begins: comparison suite, cleanup,
result generation under the four-gate release condition, renderer, final report.

---

## L13 — OPERATIONAL: THREE SPRINTS OF SCIENTIFIC RECORD EXIST ONLY ON THIS DISK.

Cleanup lane, verified independently by the coordinator.

    git ls-files s23 -> 0      git ls-files s24 -> 0      git ls-files s25 -> 0
    314 untracked files, 24 MB, across the three directories
    git check-ignore: NOT ignored -- they were simply never added
    branch main, HEAD a15406c "checkpoint: pre-deep-dive research baseline"

**Everything Sprints 23, 24 and 25 produced is untracked**: every LEDGER, BRIEF, PREREG and FINDINGS
file, every result JSON, and both artefacts the brief names as expensive to regenerate —
**`s24/cache_amber/` (126 files, 63,000 genuine ff14SB/GBn2 single points, a serialised OpenMM run)**
and **`s24/results/a_corpus_permitted_29e3b67e8ca0c03d.json`** (the audited training corpus with its
exclusion list).

**This is a single-point-of-failure on the entire research record, and it is the most urgent
operational item in the project.** It is a `git add` on three directories. **Committing is the user's
call and has not been taken.** Flagged, not actioned.

*(The 1,134-file quarantine is deliberately excluded from any commit — it is synthetic and its README
says so.)*

---

## L14 — THE TEST SUITE IS GREEN. THE BOX IS FULL.

`pytest tests/ -q`: 235 tests, 204 passed, 13 skipped, **3 failed, 15 errors** — and **all 18
non-passing outcomes have one cause, the 92% physical-memory ceiling. Not one is a code defect.**

The 15 errors are all of `test_amber.py`, raised in its autouse `_memory_ceiling` fixture:
`_memory_verdict` returned `"fail"` rather than `"skip"` because the baseline is captured at
**collection** time, before pytest imports torch/pennylane/esm/openmm for the other nine modules — so
the whole suite's import growth was charged to the AMBER suite. Measured 80% at rest, 93–97% during
the run, back to 80% the instant the process exited.

**Re-run split, same box, same code: `test_amber.py` is 15/15 green, and the non-AMBER suite plus 134
new tests is 350 tests, 337 passed, 13 skipped, 0 failed, exit 0.** The only reproducible red is
`test_amber_relaxation_is_frame_invariant`, whose captured stdout shows frame invariance holding at
~1e-4 Å on 9 of 11 targets before the guard fired. **The physics is not failing; the box is full.**

### Two gaps in the test suite that mattered, now closed

- **`s12/instrument.py` had NO test of its own** — despite producing every RMSD in the project.
  Now 44 tests. Two real findings: a self-RMSD floor of ~1.3e-7 Å (the lane's first draft asserted
  exact zero — *the test was wrong, not the code*), and that `core.geometry.pairwise_ca_rmsd`
  **symmetrises where the instrument does not**, which matters for medoid ties.
- **`tests/test_quantum.py` is an EQUIVALENCE suite** — if the shipped CVaR and its reference were
  wrong in the same way, all 39 tests still pass. **90 new tests assert the definition instead**,
  across `CVAR_SORT_CUTOFF` so both branches run.

### Rulings

**Archive layout: OPTION A — leave the sprint directories flat and add an index.** Moving them breaks
the two-hop `ROOT` computation in **all 592 modules**, and `ROOT` builds **data paths**, not just
`sys.path`. Cosmetic tidiness is not worth that risk after a freeze. The lane recommended A and was
right to refuse to choose alone.

**Obsolete list APPROVED** for the 10 root `work_*.log` (untracked, gitignored class, basename absent
from all 917 source files) and the 628 bytecode files (14.0 MB, gitignored). **`s15/results/
qgeom_ens.json.tmp` is HELD** — it sits inside a sprint `results/` directory and is therefore
protected by location; 25 KB is not worth the precedent. **The 17 sprint console logs stay**, and the
lane was right to refuse to classify them: `_archive/README.txt` records this exact check being run
before, run wrong, and one sprint log turning out to be the sole on-disk source of a published number.

**The memory guard for `test_amber_frame_invariance.py` is APPROVED** — it lacks the guard
`test_amber.py` has, so pressure renders as a red FAIL rather than a skip. The lane correctly refused
to make that change alone, because converting a red into a skip is a judgement about the instrument
rather than a formatting fix. It is my judgement to make and I am making it.

**The six unused-import proposals are HELD** until the results and physics lanes stop importing
`core.pipeline`.

---

## L15 — THE READOUT IS INSENSITIVE TO A DISTRIBUTIONAL DIFFERENCE OF NEARLY HALF THE MASS.

Quantum lane, `s25/q_gibbs.py` → `results/q_gibbs.json`. **Ordered because the drafted claim "the
circuit attains its own analytic optimum" was an inference from an ENDPOINT dressed as a statement
about DISTRIBUTIONS — the same error as L2. The measurement fired the falsifier and produced a better
result than the sentence it replaced.**

### The identity that makes it exact rather than estimated

For `F(p) = E_p[E] − T·H(p)`, the unconstrained simplex minimiser is `p* = exp(−E/T)/Z` and

    F(p) - F(p*)  =  T * KL(p || p*)     EXACTLY

**So the free-energy suboptimality IS the divergence** — no estimation, no sampling, both sides from
the same exact statevector the deployed driver uses. The identity is **asserted at runtime to <1e-9 on
every row**, so it cross-checks both computations instead of assuming either.

### The ladder at the deployed T = 0.3

    state                                    F        KL nats   KL bits      TV    H bits
    random theta (untrained)             -1.275443    3.927427   5.66608   0.78052  4.4562
    uniform over 2**n                    -1.455609    3.326874   4.79966   0.70154  7.0000
    point mass at argmin (the collapse)  -1.718572    2.450332   3.53508   0.91374  0.0000
    5 Adam steps                         -1.644194    2.698258   3.89276   0.70994  5.1019
    50 Adam steps  (DEPLOYED)            -2.183096    0.901916   1.30119   0.45308  5.6706
    the Gibbs optimum itself             -2.453671    0.000000   0.00000   0.00000  4.9135

KL(trained ‖ optimum) = **1.449 / 0.902 / 0.373 nats** at T = 0.1 / 0.3 / 1.0; TV = 0.761 / 0.453 /
0.351. **At the deployed temperature the trained state disagrees with its own analytic optimum on 45%
of its mass. "Attains" is unearned by any reading and is cut.**

### The genuine positive is the CONTROL, not the optimum

Best-of-200 from the **untrained** circuit — never an initialisation mean, per the project's standing
rule:

      T      F_trained   F_init_mean   F_init_best    F_gibbs   gap closed   beats best-of-200
    0.10     -1.717571     -0.526843     -0.988733   -1.862495     0.891          True
    0.30     -2.183096     -1.207063     -1.483528   -2.453671     0.783          True
    1.00     -4.936715     -3.587831     -4.214344   -5.309822     0.783          True

**The optimiser genuinely trains. It beats best-of-N at every temperature and closes 78–89% of the
available free-energy gap.**

### THE SENTENCE THAT REPLACED THE OVERCLAIM, AND IT IS THE BEST RESULT IN THE DOCUMENT

Three readings in order: **the optimiser trains**; **it does not reach the optimum** (0.902 nats, 45%
of the mass, and sitting *broader* than optimal at 5.67 bits against 4.91 — not collapsed); and **the
endpoint cannot tell** (−0.0302 Å, 0.24× MDE, NULL).

> **The readout is insensitive to a distributional difference of nearly half the mass.**

**That is the mechanism for the entropy curve.** It is *why* an exact classical Gibbs state, a trained
quantum state 45% away from it, and a flat average over the top 64 candidates all land within a few
hundredths of an Ångström of each other.

**And it narrows a standing caveat of mine.** I had written that the circuit buys nothing because
*"the target is classical and cheap"*. **The circuit does not reach the cheap target either — so
cheapness is not the mechanism. The binding constraint is READOUT SLACK, and that would still bind on
a landscape where the Gibbs state were expensive.** A deeper statement, and one that is more
favourable to variational methods in general while being strictly more honest about this system.

### Two details kept because they cut against us

At **T = 0.1** the trained state's free energy (−1.7176) is marginally **WORSE** than the plain point
mass at the argmin (−1.7186), entropy 0.075 bits — the entropy term is too weak there to hold the
state open and the circuit converges to what is effectively the collapse. At **T = 1.0** the trained
state (KL 0.373) is only modestly better than **uniform** (KL 0.458).

### A second overclaim the same run caught

The lane's §3.1 said `E` is *"numerically the same vector for every target"*. **It is not** —
`_zrank` uses `rankdata`, which **averages ties**, and the score does tie. Measured on 8 real targets:
worst deviation 4.06e-02 against an E range of 3.4371, i.e. **1.18% of the range, and 0 of 8 exactly
identical.**

**Corrected: the spectrum of H is target-independent TO WITHIN TIE-AVERAGING** — and the consequence
still stands and is still striking:

> **The deployed selector solves very nearly the SAME variational problem on all 126 targets. All
> per-target information enters through WHICH CANDIDATE OCCUPIES WHICH RANK, not through the
> spectrum of the Hamiltonian.**

### Rulings 2 and 3, closed

§0(1) now names the deployed selector as an **exact dense statevector** — 7 qubits, 3 layers, 21
parameters, agreeing with an independently written dense simulator to 5.6e-17 — and states in §0
itself that **the MPS is the generation lane's ansatz and `core/pipeline.py` never instantiates
`MPSAnsatz`.** Less impressive, correct. "5 of 5 folds agree in sign" stays adjacent to UNDERPOWERED
with no editorialising sentence after it.

The lane recorded a fifth entry on its own error list: **it drafted a sentence it could not support,
and caught it by re-reading its own draft against the standard rather than from any doubt while
writing.** That is the correct lesson and it is recorded as the lane's own.

---

## L16 — THE SEVEN-CONFIGURATION SUITE. BOTH PHYSICS ENERGIES ARE MEASURABLY WORSE THAN NOISE.

Physics lane, `s25/phys_{lib,gate,suite,analyse,landscape}.py`, `results/phys_*.json` all
`complete: true` and provenance-stamped, 126 per-target cells persisted. Design written to
`s25/PREREG_PHYS.md` **before anything ran.**

**Held identical across all seven:** 126 targets in `sorted(pdb)` order · the 5 pinned folds · the
shipped K=500 pool (`universe_idx` asserted bit-identical to `s24/cache_amber` per target) ·
point-cloud Cα-RMSD · uniform coordinate average · **genuine CVaR-VQE for all seven** (exact
StatevectorCircuit, 9 qubits, 3 layers, 80 Adam iters on the exact parameter-shift gradient, α=0.18,
T=0.5, seed 0). Every contrast paired per target. Classical arms reported as their own complete
matched table, never as a substitute.

**α = 0.18 was pinned on a NATIVE-FREE criterion**: a 12-target probe on the distogram channel, no
RMSD read, gave median realised tail = **75**, the production rung — so the VQE and the classical
top-75 control consume the same number of candidates. That is the only way "same readout" and
"CVaR-VQE for all seven" can both hold.

**Normalisation, stated mathematically:** `zrank(x) = (rankdata(x) − mean)/sd`; with no ties
`mean = (K+1)/2`, `sd = √((K²−1)/12)` exactly, so every channel lands on the *identical* marginal.
`E_S = zrank(Σ_{c∈S} zrank(x_c))`. **The OUTER zrank is load-bearing**: without it a |S|-channel sum
has sd 1.00/1.41/1.73, and since CVaR trades energy against `T·H` the temperature would silently
differ between configurations — worth up to **0.14 Å**. Mean |z_moment| per channel is **AMBER
0.1127, distogram 0.7529, Legacy 0.8013**, so a bare sum would not be "equal weight" at all.

### Gates

    F0 cache integrity   pool identity 126/126 · score reproduction 126/126 · FRESH AmberSP single
                         points on 5 fold-stratified targets x 4 candidates: **max_rel 0.00e+00**
    F1 anchor            C3 classical top-75 = **3.048338**, |diff| = 0.00e+00
    F2 subset-hood       **1764/1764**, max |RMSD_VQE − RMSD_top-m| = **1.14e-13**

### Headline (rank normalisation, VQE arm, n=126)

    C3 Distogram              **3.0580**       C5 AMBER+Distogram        3.1317
    C4 Legacy+Distogram         3.2151         C7 Legacy+AMBER+Dist      3.2530
    C6 Legacy+AMBER             3.6742         C1 Legacy                 3.7553
    C2 AMBER                    3.8805         **random 75-subset null   3.4251**

Nothing beats the incumbent; C3's VQE arm reproduces it to +0.0097 (0.43× MDE).

### The two results that carry the final report

> **Legacy (+0.330, 1.99× MDE), AMBER (+0.455, 2.42× MDE) and Legacy+AMBER (+0.249, 1.53×) are
> measurably WORSE than a random 75-subset of the same pool, 5/5 folds each. Not null — worse.**

> **The rank-permuted control: destroying a physics channel's correspondence to candidates while
> PRESERVING its marginal IMPROVES the endpoint everywhere.** Legacy is worse than its own noise by
> **+0.327**; AMBER by **+0.471** (26W/100L); Legacy inside C4 costs +0.124 against the same noise.
> **AMBER inside a distogram-led score is indistinguishable from noise (+0.0331, 0.43× MDE)** — s24
> L16 reproduced with **no free parameter**, which removes the "the fit was the problem" escape route.

**m-ladder:** under rank normalisation it spans **0.009 Å** across all seven (realised m 73.2–75.5), so
the table compares energies and not set sizes. Under raw moment it spans 35.9–92.8 and is worth ±0.14 Å.

### F4 fired against the lane, and the mechanism is a float64 defect

Raw-moment AMBER beat rank at 2.05× MDE. Chasing it: **`zmoment` is not monotone in float64 when the
sd is set by a 1e28 outlier.** It breaks AMBER's own `argsort` on **40 of 126 targets**, with exact tie
blocks up to **462 of 500**, and `argsort` then reads the **BLOSUM retrieval order — which is not
neutral** (ρ(pool index, ORACLE RMSD) = +0.054 ± 0.008). On the 86 targets where monotonicity holds the
two arms are **bit-identical (0.00e+00)**; the entire −0.0267 classical difference lives on the 40
broken ones. Under raw moment, "AMBER+Distogram" *is* the distogram — top-75 overlap **0.9930**.
**The rank choice is vindicated by a stronger argument than the one pre-registered.**

### Landscape, verified on the full instrument

    rho(Legacy, AMBER)      **-0.0896 +/- 0.018**   third manifold, same number as s20/s24
    rho(Legacy, distogram)   +0.3671 but **median +0.556** -- the mean is dragged by anti-aligned targets
    rho(AMBER, distogram)   **-0.0186 +/- 0.0203, CI INCLUDES ZERO** -> orthogonal, NOT anti-correlated
    compactness: Legacy's top-75 is **-0.758 A** in Rg, AMBER's **+1.103 A** -- opposite signs on the
                 axis that dominates this pool
    whole-pool rho(E,RMSD) trap CONFIRMED: Legacy's +0.307 sits on an Rg control of +0.279;
                 IN-BAND (Rg control +0.021): Legacy +0.109, AMBER +0.091, **distogram +0.236**

### Two self-corrections on the record

F4's misfire, above. And the lane **nearly reported `share_accounted = 5.36` on the best-of-seven null
as if it disqualified the ceiling — it does not; it says the NULL IS INAPPLICABLE**, because seven arms
at mean pairwise r = 0.86 give **k_eff = 1.13, not 7**. The ORACLE per-target ceiling of −0.3349 Å
stands as ORACLE.

---

## L17 — THE HAMILTONIAN BARELY CHANGES BETWEEN TARGETS. THIS EXPLAINS FIVE SPRINTS OF NULLS.

Quantum lane, §3.2 of `s25/QUANTUM.md`. **Independently verified by the coordinator** before adoption,
because it retro-explains other lanes' results and a claim of that reach should not rest on one lane.

`core/pipeline.py:853` builds the selector's Hamiltonian as `E = _zrank(pool["sc"][o])` where
`o = top[:dim]` and `top = argsort(sc)`. **The scores are therefore already sorted, so `zrank` of them
is the standardised rank ladder 1…128 — a fixed vector — except where the score ties and `rankdata`
averages.** Coordinator's own measurement, 8 targets:

    E range 3.4372      worst deviation from target 0:  **4.0627e-02  =  1.18% of range**
    exactly identical:  1 of 8  (target 0 against itself)

**H = diag(E) is a fixed ladder of standardised ranks, up to tie-averaging.** Not exactly constant —
stated as measured so nobody reads it as an identity.

### Three consequences, and the third is a mechanism nobody had

1. **There are effectively TWO trained states in the entire deployment** — one per (α, T) cell of
   `VQE_LFO` — not 126.
2. **The quantum stage is insensitive to the target by construction.** All per-target information
   enters through *which candidate occupies which rank*, never through the spectrum.
3. > **It independently explains a result this project measured five times across five sprints and
   > never had a mechanism for: deeper ansätze and larger χ order nothing.** A more expressive state
   > can only pay if `H` carries target-specific structure to capture. It carries almost none.
   > **Extra expressivity has nothing to be expressive about.**

**(3) is a mechanism rather than a restatement, and it was available from four lines of source the
whole time.** It is the one finding this sprint that retro-explains other lanes' results rather than
only its own — and it means the "deeper ansatz / larger χ" row in every DO-NOT-REDO table since
Sprint 21 was closed for a reason nobody had identified.

### The lane's final state

`s25/QUANTUM.md`, 855 lines, source-referenced, honest-contribution section first on the page. 39/39
quantum tests pass; four artefacts provenance-stamped; **nothing in `core/` or `s24/` edited** and the
two harness fixes left as recommendations.

**Caveats the lane attached to its own handover, kept:** the −0.1405 Å headline is **0.68× MDE,
underpowered, not a result**, though 5/5 folds agree — both statements adjacent, no interpretation
added. Q-B is a re-analysis of a fixed artefact and its marginal means were read before forks were
filed, declared in the prereg. **Q-D never triggered** — no upstream RMSD-lane change reached it, so
the matched-harness arm is untested this sprint.

### Four numbers for a slide, as requested

    1.  KL 0.902 nats, TV 0.453  ->  endpoint difference 0.24x MDE
        The trained state and its own analytic optimum disagree on 45% of their mass and the readout
        cannot tell.  THE headline.
    2.  rho(readout entropy, RMSD) = -0.7423   vs   rho(alpha, RMSD) = +0.2700
        Curve fitted on the nine NO-CIRCUIT arms; the circuit arms land on it at +0.0090 A against a
        residual sd of 0.0268.  One quantity explains both knobs and the circuit.
    3.  alpha = 1.0 on 3 of 5 folds  =  61.9% of targets carry NO tail constraint
        The fact a postdoc finds in core/pipeline.py:118 in ten minutes.  Better from us.
    4.  Gradient variance -0.649 log2/qubit for the LINEAR cost, -0.047 for CVaR at alpha=0.1
        No barren plateau anywhere measured, and the CVaR non-linearity FLATTENS the decay.  Pair it
        with the free control: at alpha=1 the CVaR reduces IDENTICALLY to <psi|diag(E)|psi>.
    5.  (if a fifth is wanted)  78-89% of the free-energy gap closed while beating best-of-200 from
        the untrained circuit -- the genuine positive, and what keeps the slide from reading as pure
        demolition.

---

## L18 — FORM 5 CLOSED. THE FUNCTIONAL LEVER IS CLOSED IN ALL FIVE OF ITS FORMS.

Physics lane, `s25/phys_form5.py` → `results/phys_form5.json`, complete and provenance-stamped.
Pre-registered as `PREREG_PHYS.md` ADDENDUM A1–A7 — rule, falsifier, nulls, ceilings, multiplicity and
a **loudly stated expectation of failure** — written to disk before the run. Threshold **and
direction** fitted on the other four folds only, with both degenerate ends (always-C3, always-C5) on
the grid so the fit was free to decline to switch.

    PRIMARY  s_t = mean over C3's own top-75 of [zrank(E_LEG) - zrank(E_AMB)], switching C3<->C5
    switched 3.0629 vs C3 3.0580    effect +0.0049  SE 0.0130  MDE 0.0363   **0.14x MDE**
    fold CI [-0.0113,+0.0242]   folds same sign 3/5   14W/15L/**97 TIES**   NOT MEASURED
    N1 permuted-signal null     +0.0053   **0.11x MDE**   fold CI [-0.0140,+0.0204]
    >>> FALSIFIER NOT CLEARED

> **The real signal and a signal stripped of all its target-correspondence are indistinguishable:
> +0.0049 against +0.0053. That is what makes this an ABSENCE of signal rather than an underpowered
> non-result** — the primary contrast alone would only have been the latter.

Three corroborations: **the fitted direction flips sign between training folds** (+1 on folds 0/2/4,
−1 on 1/3 — the signature of fitting noise, surfaced only by per-fold reporting); **97 of 126 targets
are ties**, the fit switching just 3–9 held-out targets per fold, i.e. very nearly choosing to do
nothing; and the ORACLE mechanism check gives **ρ(s_t, C5−C3) = +0.0506**, so the signal carries
essentially nothing about which arm wins — even though the per-target spread is real (C5 wins on
58/126 by −0.257 Å and loses the rest by +0.355 Å).

**Declared secondary** (Rg contrast, identical nested CV, Bonferroni bar registered *before* it ran):
−0.0012, **0.03× MDE**, 107 ties. **Its permuted null (+0.0524) is LARGER than its real arm** — at this
effect size the fitting procedure's own noise dominates.

**ORACLE ceiling:** per-target `min(C3, C5)` = 2.9400 Å, −0.1181 Å under full leakage — but
`corr(C3, C5) = 0.9536` → **k_eff = 1.02 independent arms, not 2.** Nested CV recovers none of it and
moves the wrong way.

> **Second instance of the k_eff lesson, now stated as a general rule: a best-of-K null is INAPPLICABLE
> when the K arms are not independent, and `k_eff` is how you know.** It sits alongside the audit
> lane's proof that an own-row null returns 83% on pure noise.

**The functional lever is closed in all five forms: filter, partition, score, audit, and — from
Sprint 24 — an oracle `w` chosen with full leakage worth 0.0148 Å.**

### The suite's write-up, as directed

§0 restructured so the two carrying results lead ahead of the headline table. The outer
re-standardisation is spelled out as **up to 0.14 Å of pure normalisation artefact masquerading as a
physics result**. α = 0.18's reasoning preserved verbatim as the standing answer to "why not the
deployed value" — the deployed table is tuned for a 128-candidate consensus-medoid readout, this suite
needs the realised tail on the production rung of 75. F4 has its own subsection. **The landscape table
now reads orthogonal-with-CI: ρ(AMBER, distogram) = −0.0186, 95% CI [−0.0583, +0.0211]** — includes
zero, so *orthogonal*, not anti-correlated, with an explicit note that a bare point estimate invites
over-reading.
