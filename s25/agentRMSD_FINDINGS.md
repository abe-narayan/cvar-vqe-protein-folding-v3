# SPRINT 25 — RMSD LANE. THE LOCATION QUESTION, ANSWERED NO, WITH THREE MECHANISMS AND ONE STRUCTURAL FINDING.

**Endpoint throughout: mean full-chain Cα-RMSD, 126 cluster-disjoint dev targets, POINT CLOUD basis,
shipped uniform top-75 coordinate average in the medoid frame. Incumbent 3.0483 Å.**
In every arm the pool, the readout and the per-pair weight `w` are IDENTICAL; only the risk table
varies, so the only thing that changes is WHICH 75 candidates are selected. Basis stated on both
sides of every contrast.

Pre-registered in `s25/PREREG_RMSD.md`, fork list sent to the coordinator before each run.
Modules `s25/loc.py` (diagnostic, no RMSD) and `s25/locrun.py` (end-to-end).
Artefacts `s25/results/loc.json`, `s25/results/locrun.json`, both `complete: true`, n=126,
provenance-stamped (`locrun.py`, sha `dc36e5be`, commit `a15406c8`).

**CAVEAT CARRIED ON EVERY FIGURE BELOW:** 4 of the 126 dev targets carry a verbatim self-sequence
window at BLOSUM rank #1, from `core/data.py`'s longer-normalised identity.

---

## THE ANSWER

> **The score's per-pair LOCATION can be moved substantially and demonstrably toward truth from
> information already available — by up to 25% of the way, reducing the location field's MAE by
> 9.2% with no native input at all — and the endpoint does not respond. Across 51 non-identity arms
> spanning `gam_eff` from −0.015 to +0.251, the endpoint stays inside a band of −0.014 to +0.073 Å,
> and `corr(gam_eff, endpoint delta) = +0.054`. Every family is NOT MEASURED against its own
> identity parameter, and for THREE of the four families the best global parameter chosen with
> COMPLETE LEAKAGE on the very targets it is scored on is the incumbent itself.**

The functional side of Phase I's open question is closed, and closed by its ceiling rather than by a
failed fit — the s24 L16 standard.

---

## 1. RUN 1 — THE LOCATION DIAGNOSTIC (`s25/loc.py`, n=126, no RMSD computed)

Natives read for diagnosis only, same standing as `calib.py`. `gam_eff = <ΔL,ΔT>/|ΔT|²` with
`ΔL = L − L_med`, `ΔT = D_native − L_med`: the fraction of the way to truth a location move travels,
in the prior ladder's own currency. **`−2.1496 × gam_eff` is an UPPER BOUND, never a prediction** —
the ladder's γ has `cos = 1` by construction and none of these moves do.

    arm            gam_eff      cos   |dL|/|dT|      MAE   on_bin_centre
    IDENT          +0.0000   +0.000     0.000     2.3968   1.000
    MODE           -0.0189   -0.068     0.258     2.4906   1.000
    POWER q=0.25   -0.0130   -0.060     0.203     2.4554   1.000
    TRUNC d=0.5    -0.0189   -0.068     0.258     2.4906   1.000
    SHELL a=1.0    +0.0136   +0.094     0.210     2.3993   0.189
      CTRL_PERM    +0.0108   +0.063     0.210     2.4431   0.270   <- its own matched null
      CTRL_XSHELL  +0.0269   +0.052     0.783     3.0204   0.079   <- ANOTHER target's profile
    QUANT t=0.60   +0.0613   +0.200     0.268     2.3461   1.000
    MEAN = POWER2  +0.0582   +0.253     0.210     2.3386   0.063
    METRIC e=1.0   +0.2497   +0.495     0.450     2.1756   0.067

### 1a. THE STRUCTURAL FINDING: THE SHIPPED SCORE'S TARGET IS QUANTISED TO 17 VALUES

The coordinator asked for `on_ctr` as a guard against small-q measuring bin quantisation. It caught
something larger. **The L1 minimiser of a discrete distribution sits on an atom**, so the shipped
score's effective per-pair target takes only the 17 values 4.0, 4.75 … 21.0, 25.0. The largest gap
between adjacent centres is **4.0 Å** (21 → 25), so the representation admits **up to 2.0 Å of pure
location error** before any estimation error; the mean local centre spacing at the incumbent's own
realised targets is **1.094 Å**. Every arm with `gam_eff > 0` has `on_ctr ≤ 0.19`; every arm with
`gam_eff ≤ 0` has `on_ctr = 1.000` — across 30 arms the sign of `gam_eff` tracks whether the arm
de-quantises. This is independent of any framing, and it is what motivated family DEQUANT.

### 1b. D1 — THE CANCELLATION TEST. MY PREDICTION WAS HALF WRONG AND I WITHDREW THE HALF THAT FAILED

I pre-registered: `pool_signed[shell] ≈ prior_signed[shell]` to within ~0.1 Å **at every shell**.

    shell  pairs  prior_med  POOL(K=500)   TOP-75   pri-pool  pri-top75
    2-2     1381    +0.0327     -0.0821   -0.1328    +0.1148    +0.1655
    3-3     1255    +0.1827     -0.2555   -0.4706    +0.4382    +0.6533
    4-5     2132    +0.2096     -0.1336   -0.3732    +0.3432    +0.5829
    6-8     2253    +0.2999     +0.3530   -0.0960    -0.0531    +0.3959
    9-15    1528    +0.4692     +1.1158   +0.4517    -0.6467    +0.0175
    ALL             +0.2318     +0.1744   -0.1310    +0.0574    +0.3628

It holds in AGGREGATE (+0.0574, inside the 0.1 Å I named) and **FAILS at 4 of 5 shells, by 0.34 to
0.65 Å, with the sign of the disagreement reversing across the range.** Withdrawn as stated.

**The coordinator's Addition 4 supplied the better mechanism.** The TOP-75 — the candidates that
actually reach the readout — sit at −0.1310, too SHORT, while the prior is too LONG at +0.2318.
With s24 L14's emitted average at −0.5852 the chain is **prior +0.23 → pool +0.17 → SELECTED −0.13 →
emitted −0.59**. Selection has already over-corrected the offset and reversed its sign, so shifting
the prior's target shorter — what the signed error naively asks for — drives the emitted structure
further into a contraction that is already the largest single operator effect in the system.

### 1c. MECHANISM 2 (MULTIMODALITY) PRE-REFUTED, AND IT WAS MY OWN FAVOURED HYPOTHESIS

I registered family A as the open one *because* the 24.1% multimodal pairs live there. Mode-seeking
is wrong-signed at every setting, and the coordinator's pre-registered stratification localises it
exactly where the mechanism claimed its gain: on multimodal pairs `TRUNC(0.5)` is **−0.0568** against
**−0.0037** on unimodal.

    multimodal fraction 0.241   (reproduces s25 L1's 0.241 on an independent code path)
    |median − native|      2.5810
    |NEAREST mode − native| 1.3614   ORACLE choice of mode, needs the native to pick it
    |TOP mode − native|     2.8307   the ACHIEVABLE mode-seeking arm — 0.25 A WORSE than the median

**The safe arm is the absurd one, for the fourth time in this project.** The median-between-modes is
correct; the mode is not reachable without the answer.

### 1d. MECHANISM 1 (PER-SEPARATION OFFSET) DEAD AGAINST ITS OWN NULL

At α=1 with an ORACLE in-sample profile, family SHELL travels `gam_eff` **0.0136**; its
magnitude-matched permuted control gets **0.0108 — 79% of it**; and **another target's profile scores
HIGHER (0.0269) than the target's own.** The shell STRUCTURE, which is the entire content of the
mechanism, is worth 0.0028. Dead before it cost an end-to-end run.

---

## 2. THE FRAME WAS RETRACTED MID-LANE, AND WHAT SURVIVED IT

This lane was commissioned under s25 L2 — "the ranking responds to LOCATION and not to WIDTH".
**s25 L7 retracts that**: the unconfounded experiment finds location and width equally flat, all
eight cells at 0.38–0.58× MDE, and the response to location SYMMETRIC about a posterior biased
−0.41 Å against the natives.

**Nothing measured in run 1 depended on the retracted framing** — it measures where the target is and
how far each functional moves it, which is true either way. What the retraction did change is the
prior on run 2, and **the only thing keeping run 2 alive was one distinction, which I stated before
running it: the audit's location arm is a RIGID TRANSLATION, one constant added to every pair; every
arm here is an ADAPTIVE PER-PAIR change whose size and sign are set by each pair's own posterior
shape.** A rigid shift being flat does not imply an adaptive one is.

**Run 2 settles that: it does. The adaptive/rigid distinction does not rescue the direction.**

On the coordinator's scoping line — genuinely better location (new information) versus a re-reading
of the location already in the posterior — placed before any endpoint was seen:

| arm | side of the line | outcome |
|---|---|---|
| POWER, TRUNC | RE-READING. Median → mean/mode of the same posterior. | null / worse |
| QUANT, SHELL | pure location shifts, DROPPED before running (§4) | closed |
| DEQUANT/DEQUANTI | strictly a RE-READING; what it removes is a REPRESENTATION defect, not an estimation error | null |
| METRIC | genuinely NEW information — the metric constraint is not in the posterior at all | worse; ceiling exactly 0 |

---

## 3. RUN 2 — END TO END (`s25/locrun.py`, n=126, complete)

Identity arm through this module's own path = **3.0483**, asserted per target both `np.array_equal`
against a genuine `core.predict.Distogram._risk` AND equal to `temper.json`'s independently computed
`f=1` column. `w` is bit-identical on both sides of every contrast, so the audit's open question
about the width channel (`w = 1/(sd+0.5)`, `_score_weights()` returning ones) cannot confound
anything here.

    family      param     RMSD   vs ident   gam_eff      cos     amp
    METRIC       0.25   3.0543    +0.0059   +0.0637   +0.497   0.508
    METRIC       1.00   3.0726    +0.0243   +0.2507   +0.496   0.508
    DEQUANTI     1.00   3.0472    -0.0012   +0.0290   +0.159   1.843
    DEQUANTI     2.00   3.0452    -0.0031   +0.0352   +0.180   2.207
    DEQUANT      1.00   3.0474    -0.0010   +0.0333   +0.187   2.151
    POWER        2.00   3.0799    +0.0316   +0.0582   +0.253   3.070
    TRUNC        1.00   3.1211    +0.0728   -0.0145   -0.060   2.000
    MM_TRUNC   mm-only  3.0347    -0.0136   -0.0128   -0.060   1.291

### 3a. NESTED CV, AND THE CEILINGS THAT ACTUALLY CLOSE IT

One global parameter per family, chosen on the training folds, applied held-out, paired against the
family's own identity arm.

    family      nested CV   SE      MDE     x MDE   fold CI            verdict
    METRIC       +0.0065  0.0029  0.0081   0.80   [+0.0000,+0.0198]  NOT MEASURED [TYPE-M]
    DEQUANTI     -0.0011  0.0026  0.0073   0.15   [-0.0064,+0.0024]  NOT MEASURED
    DEQUANT      +0.0014  0.0037  0.0105   0.14   [-0.0003,+0.0043]  NOT MEASURED
    POWER        +0.0176  0.0066  0.0184   0.96   [+0.0000,+0.0384]  NOT MEASURED [TYPE-M]
    TRUNC        +0.0273  0.0117  0.0327   0.83   [+0.0036,+0.0529]  NOT MEASURED [TYPE-M]

**The decisive row is not the fit, it is the ceiling.** The best single global parameter chosen with
FULL LEAKAGE on all 126 targets:

    METRIC     eta = 0.00  ->  3.0483  =  +0.0000    the IDENTITY is the optimum
    POWER      q   = 1.00  ->  3.0483  =  +0.0000    the IDENTITY is the optimum
    TRUNC      d   = 40.0  ->  3.0483  =  +0.0000    the IDENTITY is the optimum
    DEQUANTI   b   = 2.00  ->  3.0452  =  -0.0031    0.4x its own MDE, 0.1% of the incumbent
    DEQUANT    b   = 2.00  ->  3.0466  =  -0.0018

**For three of four families there is no value of the parameter anywhere on the grid worth having,
and the leakage-optimum IS the shipped incumbent.** For DEQUANT the leakage-optimum is 0.003 Å, and
it sits at β=2 — beyond the principled β=1, i.e. a smoothing that crosses bin boundaries rather than
the histogram's own resolution. Every per-target oracle is **124–147% accounted** by
`ST.best_of_k_within`'s valid across-target null, and every split-half transfer is **−8% to +16% of
its own oracle — NOT A SIGNAL**. Taking all 52 arms as one grid: oracle −0.2675, valid null −0.5779,
**208% accounted, split-half transfer −5%.**

*(The invalid own-row null the shared library now flags would have reported 79–85% for every one of
these — a near-constant function of K, exactly as the auditor proved.)*

### 3b. WHAT THE OPERATORS ACTUALLY DELIVERED (`ST.achieved`)

Every arm does what its name claims, so the nulls are nulls and not broken operators:

    POWER q=2 realises the posterior MEAN                    ratio 1.0000
    METRIC e realises (1-e)*L_med + e*L_proj                 ratio 1.0000 at all four e
    DEQUANT/DEQUANTI PRESERVE the posterior mean             ratio 1.000000 at every beta
    DEQUANTI b=1 sd inflation 1.043;  minimiser OFF-centre 0.914  (de-quantisation achieved)
    mean |L_DEQUANT(b=1) - L_median| = 0.2770 A;  DEQUANTI 0.2683 A

**DEQUANT is not a failed operator. It moves the target by 0.28 Å RMS on 91% of pairs, preserves the
posterior mean to 1e-6, and buys 0.003 Å with full leakage.**

### 3c. P1 HALF-CONFIRMED, P2 REFUTED — AND THE REFUTATION IS THE USEFUL PART

Both registered with the coordinator before run 2.

**P1** — METRIC has the largest `gam_eff` and the smallest endpoint movement per unit `gam_eff`.
**First half CONFIRMED**: 0.2507, 4.3× the next arm. **Second half CONFIRMED only against POWER**
(conversion +0.097 vs +0.542) **and FALSE against DEQUANT** (+0.029) — DEQUANT converts less still.
So the mechanism I proposed (the projection removes a component orthogonal to the manifold every
candidate lives on, which is therefore near-constant across candidates and invisible to the ranking)
is consistent with the METRIC/POWER contrast and is not the whole story.

**P2 — REFUTED as stated.** I predicted `amp << 1` for METRIC and `amp ≥ 1` for the others. The
ORDERING is right (METRIC 0.508, lowest by 3.6–6×) but the magnitude claim is wrong, and **`amp` does
not explain the conversions at all**: POWER has the HIGHEST `amp` (3.07) and the WORST conversion,
and `corr(amp, delta) = +0.023` across 51 arms. `amp` is not the missing discount.

### 3d. THE SINGLE MOST INFORMATIVE NUMBER IN THE LANE

Across **51 non-identity arms**, four families, two implementation styles, `gam_eff` spanning
−0.0145 to +0.2507:

    corr(gam_eff, endpoint delta)   = +0.054     positive = toward truth is slightly WORSE
    corr(|gam_eff|, |delta|)        = +0.137
    corr(cos, delta)                = -0.112
    corr(amp, delta)                = +0.023
    endpoint delta range            -0.0136 .. +0.0728 A

**METRIC at η=1 reduces the location field's MAE from 2.3968 to 2.1756 (−9.2%) and its squared error
by ~30%, using ONLY the metric realisability constraint and no native information whatsoever, and
the endpoint gets 0.0243 Å WORSE.** This is `better-matrix-worse-ranking` reproduced on the valid
126-target instrument, with an interpolation parameter, against a bit-exact identity — which is
exactly the retest the DO-NOT-REDO table demands, and it closes the direction properly this time
rather than on the retracted decoy bank.

**Consequence for the prior ladder's currency.** The ladder prices γ at −2.15 Å per unit *when the
move is toward the native along the native's own direction, `cos = 1` by construction*. A move that
travels 25% of the way to truth at `cos = 0.50` is worth **+0.024 Å**. γ is redeemable only in the
ladder's own construction; it is not a general exchange rate for location improvements, and nobody
should quote `−2.15 × gam_eff` for a real operator without the cosine beside it.

### 3e. THE MULTIMODAL STRATIFICATION (pre-registered by the coordinator before run 1)

Target-level split at the median multimodal fraction, plus the pair-level version (the intervention
applied to multimodal pairs ONLY, which is the sharper test):

    arm         param   high-MM(63)  low-MM(63)      MM-only arm, paired vs identity
    METRIC       1.00      +0.0411     +0.0070       —
    DEQUANTI     1.00      -0.0026     +0.0003       MM_DEQUANTI  -0.0004  0.27x MDE  NOT MEASURED
    POWER        2.00      +0.0550     +0.0073       MM_POWER     +0.0663  0.57x MDE  NOT MEASURED
    TRUNC        1.00      +0.0935     +0.0514       MM_TRUNC     -0.0136  0.18x MDE  NOT MEASURED

`MM_TRUNC` is the only arm in the sprint with a negative full-leakage mean (−0.0136) and it is
**0.18× its own MDE with 49W/77L and folds disagreeing 3/5** — noise, and I am not promoting it.
Note it also contradicts §1c's diagnostic, which is another reason to leave it alone.

---

## 4. TWO FAMILIES CLOSED WITHOUT SPENDING A RUN

**QUANT (self-scaled quantile shift) and SHELL (per-separation offset)** are pure location shifts.
Run 1's D1 predicted they would fail end-to-end despite positive `gam_eff`; the audit lane
independently measured exactly that — a global translation is SYMMETRIC and null (−0.30 Å costs
+0.0206, +0.30 Å costs +0.0191) and a separation-graded shift at r=+0.977 to L1's own signed-error
profile is null (+0.0270, 0.41× MDE), and it had already been running unlabelled inside `temper.py`'s
SD arm. Both dropped before running. **SHELL is closed by two independent routes**: dead against its
own permuted null in run 1, and dead on the endpoint in the audit's A3/A8.

---

## 5. SCORECARD ON MY OWN PRE-REGISTERED PREDICTIONS

| # | prediction | outcome |
|---|---|---|
| D1 | pool matches prior at every shell → B/C fail | **HALF WRONG** (aggregate yes, 4/5 shells no); the consequence for B/C was right, via a better mechanism I did not supply |
| 2 | METRIC small or negative | **CONFIRMED** — +0.0243, and its full-leakage ceiling is exactly 0.0000 |
| 3 | family A open because the multimodal pairs live there | **REFUTED by my own run 1** — mode-seeking is wrong-signed, and worst exactly on the multimodal stratum |
| DEQUANT | my late primary, promoted on the quantisation finding | **NULL** — 0.15× MDE nested, −0.0031 Å with full leakage |
| P1 | METRIC largest gam_eff, worst conversion | **HALF CONFIRMED** — largest yes; worst conversion only against POWER |
| P2 | `amp` explains the conversions | **REFUTED** — `corr(amp, delta) = +0.023` |

Four of six wrong or half wrong. The two that landed (METRIC's ceiling, mechanism 2's pre-refutation)
both damaged directions I had argued for.

---

## 6. WHAT THIS LEAVES, STATED WITHOUT PADDING

**<3.0 Å was not reached by this lane and I am not manufacturing it.** The best honest number
remains the incumbent **3.0483 Å**, point cloud, and nothing in `locrun.json` beats it past its own
MDE, in-sample or out.

**What is now closed, with mechanisms:**
1. Re-reading the posterior's location by any risk functional — median → mean, mode, quantile,
   truncated, de-quantised. Three of four families have the incumbent as their leakage-optimum.
2. Metric realisability as a location correction. Retested properly and closed again.
3. Per-separation and global location shifts (two routes each).
4. Mode-seeking / multimodality exploitation — the achievable mode is worse than the median it
   would replace.

**The one thing this lane produced that outlives it** is §1a: the shipped score's effective per-pair
target is quantised to 17 values with a 4.0 Å top gap and a 1.094 Å mean local spacing. Removing
that exactly is worth 0.003 Å, which is itself the strongest available statement that **the ranking
is insensitive to the location of its own target** — a third independent route to the audit's A8,
on a completely different operator class.

**Where I would look next, and it is not in this lane.** Every arm here, and every arm in A8,
operates on the CONSUMPTION of a fixed posterior. The prior ladder's steepness is a property of the
posterior's CONTENT. This lane's result says the two are not connected by any re-reading: you cannot
buy γ by consuming the same distribution differently, no matter how much closer to truth you move
its summary statistic. If <3.0 Å is reachable at all it is through a distogram that is actually
better — which is an achievability question about a trained predictor, is not answered by any ladder,
and is not a functional question.
