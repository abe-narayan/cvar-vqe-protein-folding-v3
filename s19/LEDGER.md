# SPRINT 19 — LEDGER

---

## L1 — THE COORDINATOR'S OPENING HYPOTHESIS IS REFUTED IN ITS MECHANISM (2026-09-06, AGENT D)

**Hypothesis** (`s19/BRIEF.md` §3, coordinator): the objective path collapses a 17-bin predicted
distribution to (mean, sd); on multimodal pairs the mean lands between the modes at an unrealisable
distance, **coherently across correlated pairs** because a region flipping between conformer families
moves many pairs together — a mechanism for Sprint 18's "errors worse than random errors of the same
magnitude."

**Status: REFUTED in mechanism, partially SUPPORTED in premise.** n = 126, artefacts complete
(`s19/results/D_P1P2P3_modality/`, `D_P2P3_hardened/`; code `s19/d_modality.py`, `d_modality2.py`).

| test | result |
|---|---|
| multimodal–multimodal adjacency vs separation-stratified permutation null | 0.072 vs 0.070, **z = +0.27 [+0.08, +0.46]** — scattered, not clustered |
| spatial coherence of the *implied correction* | adjacent 0.540 vs non-adjacent 0.537, **+0.003 [−0.001, +0.008]** — zero |
| spatial coherence of the *real* residual | **+0.0385 [+0.0320, +0.0457]** — REAL |
| …on **unimodal** pairs | **+0.0429 [+0.0336, +0.0537]** |
| …on **multimodal** pairs | **+0.0489 [+0.0268, +0.0731]** |

> **THE DISSOCIATION THAT KILLS IT.** The error's coherence — the thing the Sprint-18 result is
> actually about — is **not carried by multimodality.** It is present in full on pairs with no
> multimodality at all.

**The "aims at unlikely values" harm does not survive its confound.** Unmatched, multimodal pairs
look +0.410 Å worse at the aim point, but carry mean sd 2.16 against 1.25. Matched within
(separation × sd-decile) cells the sign **reverses**: −0.174 [−0.408, +0.072] (**NOT MEASURED**,
W/L 71/49); under regression adjustment on sep × log sd it reverses **significantly**, −0.422
[−0.595, −0.248]. **Conditioned on spread and separation, the mean of a multimodal predictive
distribution is **no different** from the mean of a unimodal one.** *(Amended: AGENT D's
non-parametric nearest-neighbour match gives +0.007 [-0.291, +0.316], n_t = 121, W/L 69/52 -- a
dead null. The significant -0.422 is regression-adjusted ONLY and does not survive the
non-parametric control; it is dropped from the headline at AGENT D's own request. The harm is
REFUTED; the benefit is NOT SUPPORTED.)*

**THE CENSUS REPRODUCES BUT WAS OVERSTATED.** At n = 126: 0.241 multimodal (the n = 12 value 0.217
reproduces), mass-within-1 Å 0.634 (vs 0.641), entropy 1.44 (vs 1.43). D's own bin-artefact attack
**failed** — correcting mass to density for the non-uniform `BIN_EDGES` changes nothing (0.240) —
but **60% of flagged pairs are shoulders, not modes**: 0.097 under a prominence gate, 0.155/0.165
after 2× coarsening. **The honest figure is ~0.10, not 0.22.** The coordinator's brief quoted 0.217
with no prominence criterion.

**WHAT SURVIVES.** ORACLE nearest-mode beats the mean by **−1.068 Å [−1.186, −0.953]** per pair, but
**61% is pure min-of-N** (matched-|offset| random-sign null: −0.654). Mode information proper is
**−0.414 [−0.502, −0.324]** (−0.652 under the prominence definition). **Both native-free mode rules
are WORSE than the deployed mean**: highest-density +0.315, uniform-random +0.529. So the modes carry
real information that no native-free rule yet extracts.

**A COORDINATOR COUNTER-ARGUMENT IS VOID, on an EXACT point.** The brief argued "selection already
consumes the full distribution, so why did it not show there?" `Distogram._risk[·,x] =
Σ_b p_b |x − c_b| w` is a non-negatively-weighted sum of absolute values, hence **convex in x with
minimiser the weighted median of `CENTRES`**. **EXACT — an identity, not a measurement.** The
selection path consumes the distribution in the letter but is **structurally incapable of expressing
multimodality**; its silence is not evidence either way.

**COORDINATOR ERROR, recorded as such.** The mechanism was asserted in a binding shared document with
no measurement behind the coherence claim, and three lanes read it as motivation. Corrected in place
before any lane built on it. **This is the same failure mode as Sprint 18's §4 equivalence error: the
coordinator stating a plausible mechanism in a binding document without measuring it first.** Two
sprints running, the same error class, caught both times by a workstream rather than by the
coordinator.

**THE LEAD IT LEAVES, now the sprint's central target.** The residual's spatial coherence is real and
lives on **unimodal, confident pairs**. What makes it coherent *there* is where the Sprint-18
mechanism actually is. Routed to AGENT A.

**Still running**: D's P5 headroom test on the `s17/refine.py` instrument, `sd` held bit-identical
across arms with only the target vector changed. Validity gate already passes — D's raw arm
reproduces `refine_full` to max |Δ| = 0.0000.

**Coordinator's `s19/distobj.py`** (n = 126, running) is the structural (RMSD) counterpart and is
**pre-committed to close if `mode` does not beat `moment`** — which D's per-pair measurement predicts
it will not.

---

## L2 — THE STRUCTURAL TEST AGREES: the moment-collapse branch is CLOSED (2026-09-06, coordinator)

`s19/distobj.py`, **n = 126, `complete: true`**. Start = the coordinate average (3.048 A). Every arm
native-free; the native is read only to score. `mode` differs from `moment` in the **target value
alone** — same weights, same functional form, same optimiser.

| arm | RMSD | median | vs moment [95% CI] | W/L | vs avg |
|---|---|---|---|---|---|
| coordinate average (start) | **3.048** | 2.837 | — | — | — |
| moment (DEPLOYED) | 3.610 | 3.370 | +0.000 | — | +0.561 |
| **mode** | **3.647** | 3.523 | **+0.037 [-0.009, +0.088]** | **60/66** | +0.598 |
| nll | 3.493 | 3.251 | -0.116 [-0.257, +0.033] | 70/56 | +0.445 |
| risk | 3.600 | 3.372 | -0.010 [-0.078, +0.059] | 67/59 | +0.552 |
| unimodal *(CONTROL)* | 3.627 | 3.454 | +0.018 [-0.063, +0.096] | 63/63 | +0.579 |
| modeshuf *(CONTROL)* | 3.743 | 3.508 | +0.133 [+0.059, +0.208] | 47/79 | +0.695 |
| matched *(CONTROL)* | 4.437 | 4.514 | +0.827 [+0.587, +1.079] | 32/94 | +1.388 |

**THE PRE-REGISTERED PRIMARY FAILS ON BOTH CLAUSES.**

1. `mode - moment = +0.037 [-0.009, +0.088]`, 60W/66L — CI spans zero, magnitude **below the 0.084 A
   MDE**, and the point estimate is in the **wrong direction** (the mode is worse than the mean).
2. The pre-registered secondary — that any gain concentrates on multimodal targets — fails
   **completely**: high-multimodal half **+0.040 [-0.034, +0.130]**, low-multimodal half **+0.034
   [-0.013, +0.090]**. Indistinguishable. **There is no concentration where the mechanism requires
   it.**

**BRANCH CLOSED, as pre-committed before the run.** Two instruments now agree from different
directions: AGENT D's per-pair measurement (native-free mode rules at +0.315 and +0.529, i.e. worse
than the mean — L1) and this structural test. No disagreement to explain.

**THE COMPETING EXPLANATION I PRE-REGISTERED IS THE ONE SUPPORTED.** `unimodal` — the deployed
objective with the multimodal pairs **dropped entirely** — costs +0.018 [-0.063, +0.096], 63W/63L.
**Deleting those pairs costs nothing.** Together with D's finding that, conditioned on separation and
sd, a multimodal mean is at least as close to the truth as a unimodal one, this says the `1/sd^2`
weighting **had already neutralised them**. There was never anything there to recover. Recorded in
the pre-registration as the trap to watch for, and it is what happened.

**The modes are not noise, they are merely not enough.** `modeshuf` (real offsets, permuted across
pairs) is +0.133 [+0.059, +0.208], significantly worse than `moment`, while the true mode assignment
is +0.037. So the mode targets carry **real positional information** — about 0.096 A of it — and it
is **not enough to beat the mean**. This is the third instance in the programme of information that
is present and unusable (cf. `leg_contact`: +0.080 in-band information, selects 0.136 A worse than
random).

**THE ONE FAVOURABLE POINT ESTIMATE, and why it does not reopen anything.** `nll` = 3.493, **-0.116
[-0.257, +0.033]**, 70W/56L — the only arm pointing the right way, and its magnitude exceeds the MDE.
**Label: NOT MEASURED** (§9 — a zero-spanning CI is not a result), and it is confounded relative to
`mode` because it changes the loss *shape* as well as the target. But the decisive point needs no
further compute:

> **Even at its own point estimate, `nll` lands at 3.493 — still +0.445 A WORSE than the 3.048 A
> structure the pipeline already builds.** A powered rerun could only establish a better way of doing
> something harmful. It does not rescue refinement and it does not bear on the mission.

Filed **OPEN but not worth compute**, with that reasoning attached so a future sprint does not spend
on it.

**COORDINATOR NOTE.** The pre-registered expectation was "a real but partial effect, near 3.4-3.5,
concentrated on multimodal targets". The magnitude guess was roughly right for `nll` and the
*direction* and the *concentration* were both wrong for the arm that was designed to be decisive.
The design held: because `mode` changed one thing, the null is interpretable rather than confounded.

---

## L3 — THE MECHANISM IS FOUND: the harm is the residual's CROSS-PAIR SIGN CORRELATION (2026-09-06, AGENT A)

`s19/a_coh.py`, `s19/results/a_coh.json`, **n = 126, COMPLETE**. Instrument is `s18/objceil.py`'s
exact harness; reproduction gate passes exactly (avg 3.048, real 3.610). **ALL ARMS ORACLE** — the
native is read to construct every arm, so nothing here is deployable.

| arm | RMSD | resid RMS | kappa | EDM defect | tri % | vs real |
|---|---|---|---|---|---|---|
| real (deployed distogram) | 3.610 | 3.205 | 0.809 | 0.286 | 4.09 | — |
| signflip *(as first published — magnitude confounded)* | 2.368 | 3.062 | 0.352 | 0.340 | 10.85 | −1.242, 113W/13L |
| **signflip_exact** *(residual RMS matched to machine precision)* | **2.408** | **3.2049** | — | — | — | **−1.202 [−1.408, −1.004]**, 112W/14L, 5/5 |
| shuffled | 2.649 | 3.008 | 0.367 | 0.361 | 11.45 | −0.961 [−1.184, −0.740], 100W/26L |
| iso | 2.655 | 3.006 | 0.367 | 0.368 | 12.02 | −0.955 [−1.190, −0.733], 94W/32L |
| **coherent** | **3.749** | 3.097 | 0.990 | 0.002 | 0.02 | **+0.139 [+0.023, +0.253]**, 56W/70L |
| **coh_scaled** | **3.843** | 3.205 | 0.987 | 0.010 | 0.57 | **+0.233 [+0.122, +0.348]**, 46W/80L |

**`signflip` is the decisive arm.** `dtrue + r * eps` with eps = ±1 i.i.d. per pair: every pair keeps
**its own residual magnitude and its own weight** (no `sd[pi]` orphaning — the Sprint-18 G5 trap),
and the spatial pattern of |r| is untouched. **The only thing destroyed is the cross-pair sign
correlation.**

> ### Sign-flipping alone recovers **128%** of the entire Sprint-18 gap.
>
> real − shuffled = +0.938 · **real − signflip = +1.202** · signflip is significantly **better** than
> shuffled, **−0.264 [−0.448, −0.084]**, 83W/43L, 5/5 folds.

**The effect is entirely the correlation structure.** "Which pair carries which magnitude" — the
thing the coordinator's retracted Sprint-18 L4/G5 lead was about — is worth 0.28 Å in the
**opposite** direction.

**THE POSITIVE CONTROL CLOSES IT.** `coherent` replaces the residual with that of a **real
alternative structure** (the top-75 pool member whose residual RMS is closest to the real one) —
perfectly realisable, matched magnitude. **It is the worst arm on the board**: +0.139 [+0.023,
+0.253] worse than the real distogram, and +0.233 [+0.122, +0.348] rescaled to the exact residual
RMS. A perfectly realisable error is *more* harmful than the distogram's own.

**AT MATCHED MAGNITUDE THE WHOLE FAMILY IS MONOTONE IN COHERENCE.** Every arm sits at residual RMS
3.01–3.21, and RMSD is monotone in kappa (the fraction of the predicted deviation some
ideal-geometry structure can realise): kappa 0.35–0.37 → 2.37–2.65 Å; 0.809 → 3.610; 0.99 →
3.75–3.84. **Span 1.475 Å at constant error size.** The deployed distogram sits at **kappa 0.809** —
four fifths of its error is realisable by a *wrong* structure.

**A CORRECTION TO THE COORDINATOR'S ITEM-3 FRAMING, in the direction he did not expect.** The
predicted field **is** badly non-Euclidean in absolute terms (rank-3 EDM defect 0.286, 4.09% triangle
violations, against the native's 0.002 / 0.02%). **But that is a consequence of error magnitude, not
a mechanism**: a magnitude-matched *incoherent* field is worse on both (0.340 / 10.85%) and lands
1.24 Å **better**.

> **REALISABILITY IS WHAT MAKES AN ERROR HARMFUL** (Pearson +0.984 across the six arm means;
> per-target Spearman of (kappa, RMSD) positive on 116/126). *Phrasing tightened 2026-09-06:
> the original double negative -- "unrealisability is anti-correlated with harm" -- is
> logically identical and correct, but the coordinator inverted it while relaying it to
> another lane. AGENT D caught the relay. Stated positively from here on.* "Independent pairwise marginals are
> not jointly realisable" is measured, real, and **not the mechanism**. Nobody should design against
> it. Kept as a labelled diagnostic only.

**STATUS: this supersedes every prior account of the Sprint-18 finding**, including the coordinator's
refuted moment-collapse mechanism (L1/L2) and Sprint 18's own "which pairs the errors fall on"
(retracted at S18 L4). The answer is the **sign correlation**, and the ceiling it exposes is large:
destroying it alone reaches **2.368 Å**, below the sprint's primary target.

**THE BINDING CAVEAT.** kappa and every arm read `dtrue`. **There is no native-free estimator of
coherence yet.** This is an ORACLE ceiling, not a method.

**NEXT, owned by AGENT A**: (a) which *named* coherent mode owns the gap — per-target offset, global
stretch, separation profile, per-residue additive effect — each surgically removed with a
magnitude-matched twin through the same fit; (b) **SOURCE ATTRIBUTION** by cross-predictor error
alignment: deployed MLP vs a second seed vs PairNet (different architecture, same ESM inputs) vs the
retrieval pool's own mean-distance error. **If the coherent component is shared across independent
families it is an inputs/identifiability limit and no loss or architecture fixes it; if it is
family-specific it is the architecture.** That fork decides whether a predictor intervention is worth
spending on, and it is the single most valuable measurement now open.

---

## L4 — THE DISTRIBUTION-SHAPE DIRECTION IS CLOSED ON ITS MECHANISM, and two numbers are corrected (2026-09-06, AGENT D)

### The extractability bound the coordinator asked for

`s19/d_extract.py`, `s19/results/D_P6_extract/`, **COMPLETE**. Per-pair |aim − true| against the
deployed debiased mean, target-level bootstrap, n = 126, pinned 5 folds. Features: the full 17-bin
vector plus 14 shape summaries.

    ORACLE_bestmode                   -1.095 [-1.201, -0.987]
    min-of-N matched null             -0.611 [-0.710, -0.515]
    mode information proper           -0.484 [-0.570, -0.401]
    ---- what any NATIVE-FREE rule actually gets ----
    leave-fold-out ridge, shape       -0.069 [-0.177, +0.050]    spans zero
    leave-fold-out GBM, shape         -0.045 [-0.158, +0.077]    spans zero
    leave-one-TARGET-out ridge        -0.096 [-0.207, +0.021]    spans zero
    within-fold split-half GBM        +0.057 [-0.075, +0.190]    spans zero, WRONG SIGN
    lfo label-permutation (zero-info) -0.003 [-0.033, +0.027]    control at zero, correctly
    in-sample in-fold GBM             -0.982   ** IN-SAMPLE. NOT a bound. Refused as a ceiling. **

**Zero of the −0.484 Å transfers.** The split-half arm is the one that matters: fitted and evaluated
**inside the same fold**, on the same distribution, out of sample — **still null**. So this is not
cross-fold distribution shift; **the relationship is not there to learn.**

On ALL pairs the only thing that transfers is a better separation/mean/sd map: `lfo_noshape` −0.117
[−0.212, −0.020]; adding the **entire** distribution shape moves it to −0.152 [−0.245, −0.053].
**Shape is worth ~0.035 Å over no shape at all, inside the noise** — and per BRIEF §5 a per-pair MAE
gain is not an RMSD gain.

> **The distribution-shape direction is closed on its MECHANISM, not merely on an outcome.** The
> information exists (−0.484) and is **not extractable** — in-sample-matched, out-of-sample, by a
> linear or a boosted rule, from the full 17-bin vector. That is a stronger closure than an RMSD
> null, because it says why.

### AGENT D corrects its own significant number, and the coordinator's write-up of it

The coordinator wrote *"conditioned on spread and separation, the mean of a multimodal predictive
distribution is **at least as close** to the truth as the mean of a unimodal one"*, leaning on the
regression adjustment — the most model-dependent of the three controls. A non-parametric
nearest-neighbour match (same target, exact same separation, nearest sd, refusing |Δsd| > 0.25):

    matched multi - uni = +0.007 [-0.291, +0.316],  n_t = 121,  W/L 69/52.

Dead null. Four controls now read: unmatched +0.410 (confounded) · cell-matched −0.174 [−0.408,
+0.072] NOT MEASURED · regression-adjusted −0.422 [−0.595, −0.248] · **NN-matched +0.007 [−0.291,
+0.316]**.

> **CORRECTED CLAIM: multimodality has NO measurable effect on the accuracy of the mean once
> separation and spread are controlled. The harm is REFUTED; the benefit is NOT SUPPORTED.**

The significant regression figure is **dropped from the headline** — it does not survive a
non-parametric control. **AGENT D volunteered the retraction of its own significant result rather
than let the brief carry it**, which is the behaviour this programme most needs.

### Three defects AGENT D found in the coordinator's live `s19/distobj.py`

All of the §11 kind: the source contradicted the docstring. **All confirmed by the coordinator
against the source and fixed; the mass-mode artefact voids L2's primary arm and it is being rerun.**

1. **The weighting claim was false.** `_fit`'s docstring said every branch kept `1/sd²`; `nll` and
   `risk` carried **no weight at all**, so they changed the target *and* the weighting at once.
   Sprint 18 priced the weighting at +0.153 [+0.066, +0.247] — **the confound was the same size as
   the effect being looked for.** Matched-weight arms `nllw`/`riskw` added; `nll`/`risk` retained and
   **relabelled as confounded**. `mode` was never affected.
2. **`nll` was a negative log bin MASS, not a density.** On the non-uniform grid `−log w_b` spans
   **2.079 nats** — a spurious reward for wide long-distance bins, pulling outward in exactly the
   direction the shipped separation debias exists to remove. Fixed.
3. **`mode` used the MASS argmax, not the density argmax.** At n = 126 over 8,549 pairs the two
   disagree on **14.0%** of pairs; where they disagree the mass mode sits **+2.72 Å further out**;
   the mean shift is **+0.380 Å outward**. The distogram already over-predicts (+0.509 Å global,
   +1.492 Å at separation 11–15) and the coordinate average contracts 22–25%, so **L2's primary arm
   carried a +0.38 Å expansion confound in the harmful direction.** Fixed to the density mode and
   **the whole n = 126 run is repeating**; the superseded artefact is renamed
   `_SUPERSEDED_distobj_massmode.json`.
4. Minor: peak counting ran over interior bins only, so a mode in the first or last bin was
   invisible. Fixed by padding.

**L2's CONCLUSION IS PROVISIONAL UNTIL THE RERUN LANDS.** Its qualitative reading is still expected
to hold — AGENT D's independent per-pair measurement (native-free mode rules worse than the mean) and
the extractability bound above both point the same way — but the coordinator's primary arm was
measurably confounded and the number must be replaced, not defended.

---

## L5 — L2's "the modes carry real positional information" is DOWNGRADED to INCONCLUSIVE (2026-09-06, AGENT D caught it; coordinator's inference)

L2 read `modeshuf` (+0.133 [+0.059, +0.208]) against `mode` (+0.037) and concluded the mode targets
carry **~0.096 Å of real positional information**. **That inference is not supported by the control
as built.**

`modeshuf` permutes the mode offsets **globally**, which does not only move each offset to a
different pair — it also **orphans each offset from its own pair's weight and separation**. Sprint 18
measured that operation on its own and found it large: `shufr_wperm − shufr_wkeep = −0.467 [−0.598,
−0.338]` (S18 G6a), recorded there as *"orphaning a residual from its weight is bad"*. **The
coordinator built a control containing a confound the programme had already measured and named, one
sprint earlier.**

AGENT D's scaling argument, offered explicitly as an argument and not a measurement: S18 permuted
full residuals of mean magnitude ~2.4 Å at a cost of 0.467 Å; the mode offsets here are ~0.9 Å, so a
roughly linear orphaning tax alone predicts ~0.17 Å — **larger than the +0.133 observed**. The data
are therefore consistent with the mode offsets carrying **zero** positional information.

**STATUS: INCONCLUSIVE.** The claim is withdrawn from L2 pending measurement. It matters beyond L2
because "information present but unusable" is being written into the record as the programme's third
instance of that pattern; the **unusable** half is solid (L4's P6, n = 126, proper min-of-N null,
null even split-half *inside* a fold), and it is only the **present** half that rested on this
control.

**THE FIX IS THE ONE SPRINT 18 ALREADY INVENTED, and it is now running.** A `modeshuf_strat` arm
permutes the offsets **within (separation bin × sd quartile) strata**, so each offset lands on a pair
of comparable weight and separation:

    modeshuf_strat - mode        prices POSITION with alignment held
    modeshuf - modeshuf_strat    prices the ORPHANING

The n = 126 rerun was stopped at 15/126 and relaunched with the arm included, rather than run twice.

**Process note.** The coordinator has now had four separate calls corrected by workstreams in this
sprint alone (the §3 mechanism, the mass-mode confound, the unweighted `nll`/`risk` arms, and this
inference), against zero self-caught. In Sprint 18 the ratio was four to one the same way. **The
adversarial lane is not a luxury in this programme; it is load-bearing.**

---

## L6 — THE IN-MANIFOLD MECHANISM: measured, scaled, and deliberately NOT promoted (2026-09-06, AGENT D)

*Third and final version of this entry. It was built on algebra, contradicted by a first measurement,
and repaired by a corrected one. The sequence is preserved because it is the point: the analytic
reasoning was right, the implementation of the test was wrong, and only running it twice separated
those. Superseded artefact `D_T3_tangent/` (QR defect); canonical `D_T3_tangent2/`, n = 126, complete,
code `s19/d_tangent2.py` (SVD-truncated projection onto the true column space).*

### The geometry — EXACT, and verified numerically rather than assumed

An ideal-geometry Cα trace of n residues has **n−2 virtual angles and n−3 virtual dihedrals**, and
every pair distance is a function of those alone. So the whitened distance Jacobian has rank
**2n − 5**. Measured on all 126 targets:

    mean rank(J_w) = 20.92        mean (2n - 5) = 20.92        mean 2n = 25.92

**Label EXACT.** The generic in-tangent fraction is therefore `rank/npairs`, and the empirical null
lands on it:

    rand_white  (isotropic in the WHITENED space)   0.330 [0.318, 0.341]     <- measurement
    rank(J_w)/npairs                                0.329 [0.318, 0.340]     <- prediction

*The earlier `(2n−2)/P = 0.381` was a genuine over-count and stays struck. The corrected `2n−5` was
right all along; what was wrong was a `np.linalg.qr` projection that returns 2n orthonormal columns
even when `J_w` is rank-deficient, inflating every absolute fraction. The tell was the empirical null
coming out at 2n/npairs instead of rank/npairs.*

### The mechanism, with a scale — ESTABLISHED

    real       0.695 [0.669, 0.722]      coherent    0.846 [0.822, 0.870]
    signflip   0.397 [0.373, 0.420]      rand_white  0.330 [0.318, 0.341]   <- geometric null
    shuffled   0.586 [0.565, 0.608]      rand_raw    0.581 [0.567, 0.595]   <- iso/shuffled's OWN null
    iso        0.577 [0.553, 0.601]

**The clean matched contrast** — `signflip` is the only control preserving each pair's *whitened*
magnitude exactly, so it changes the cross-pair sign pattern and nothing else:

    real - signflip     +0.298 [+0.262, +0.331]   W/L 121/5
    real - rand_white   +0.365 [+0.338, +0.392]   W/L 125/1

> **The distogram's error field puts 0.695 of its whitened energy inside a subspace that a generic
> direction fills to 0.330, and that a sign-randomised version of the SAME field fills to 0.397.**

**Every incoherent control sits at its own correct null.** `iso` 0.577 and `shuffled` 0.586 against
`rand_raw` 0.581 — they are isotropic in **raw** distance space, and after whitening by 1/sd that is
not isotropic, so 0.581 is their reference, not 0.330. The whitening anisotropy is itself measured:
**+0.247 [+0.226, +0.268], W/L 121/5.** *(The withdrawn "shuffled and iso are far above the null" was
a comparison against the wrong null.)*

### THE PRE-COMMITTED DECISION RULE FIRES AGAIN — AGAINST PROMOTION

    across the 5 ARM MEANS, pearson(start tangent, AGENT A's realised kappa)   = +0.889   HIGH
    per-target rho(in-tangent(real) - in-tangent(signflip), gap)
        = +0.134, p = 0.135, 95% CI [-0.057, +0.301]                                     LOW
    per-target rho(in-tangent(real) - rand_white, gap) = +0.160, p = 0.074
    per-target rho(in-tangent(real) raw,          gap) = +0.157, p = 0.080

In AGENT D's own pre-registered words: **the statistic orders arms and not targets, the tangent
projection is a reference and not a model, and the design consequence stays at SUPPORTED.** The
corrected projection moved the per-target correlation from +0.143 to **+0.134** — it did not rescue
it, **and AGENT D states it would not have promoted on +0.16 either.**

### A REAL DISSOCIATION, unexplained — OPEN

AGENT A's *realised* kappa, measured after the fit, **does** predict the per-target gap:

    rho(kappa(real),                   gap) = +0.360   p = 3.4e-05
    rho(kappa(real) - kappa(signflip), gap) = +0.482   p = 1.1e-08

while the start-point linear projection does not (+0.134), **even though the projection tracks kappa
itself well per target (rho ~ +0.66)**.

> **The realised, post-fit quantity carries per-target information about the outcome that the
> start-point geometry does not. Something happens along the ~5-radian trajectory that neither lane
> is measuring.** Recorded **OPEN**; AGENT D states plainly it has no explanation.

**CAVEAT VOLUNTEERED BEFORE ANYONE LEANED ON THE +0.482**: kappa and the gap derive from the same two
fits, so a common cause — targets where the fit converges cleanly having both low objective and low
RMSD — is **not excluded**. Suggestive, not clean. **A clean test needs a native-free predictor of
kappa, which does not exist.**

### THE DESIGN CONSEQUENCE — recorded, and NOT promoted

> If the harm is the in-manifold component, **no change to the loss can fix it**: any objective whose
> argmin lies on the manifold realises whatever part of the error lies in the manifold.

**Held at SUPPORTED.** What makes it credible is the retrodiction — it explains four Sprint-18
results that were never connected (*the functional form is sound*, *the weighting is validated*, *the
optimiser is not the constraint*, *tempering does not help*) as four faces of one fact. What stops it
being promoted is that the linear statistic does not predict per-target outcome, and the per-target
evidence that does exist carries an unexcluded common cause.

**Two live levers, unchanged:** (1) reduce the in-manifold component of the predictor's error — a
training-loss change **in the right basis**, in-manifold vs orthogonal rather than confident vs
unconfident; **gated on AGENT A's cross-family attribution**, because if that component is shared
across independent predictor families it is an identifiability limit and no loss in any basis removes
it. (2) **Do not consume the argmin at all** — nobody is working on this.

### LABELS, as AGENT D asked for them

| | |
|---|---|
| **EXACT** | rank(J_w) = 2n − 5; generic in-tangent fraction = rank/npairs (predicted 0.329, measured 0.330) |
| **ESTABLISHED** | real 0.695 vs matched sign-randomised 0.397, +0.298 [+0.262, +0.331], 121/5; each incoherent control at its own correct null; whitening anisotropy +0.247 [+0.226, +0.268] |
| **SUPPORTED** | the mechanism, and its design consequence |
| **NOT MEASURED** | per-target predictiveness of the linear in-tangent excess, +0.134 [−0.057, +0.301] |
| **OPEN** | why post-fit kappa predicts the per-target gap when start-point geometry does not |
| **STRUCK** | `(2n−2)/P = 0.381`; the claim that AGENT A's kappa values sat on a dimension-ratio null |

### PROCESS

**Five corrections to AGENT D's own output in one day** — the P5 file contention, the "at least as
close" overstatement, the 2n−2 over-count, the QR rank defect, and the claim that rested on it.
**Four of the five were caught by running a measurement rather than by re-reading the argument**, and
the one time algebra was trusted over a control it took two rounds to fix. That is the entire case
for the adversarial lane, and it is also why the design consequence is not being promoted on the
strength of an argument that has now looked airtight three times.

---

---

## L7 — THE P5 HEADROOM TEST: the branch has NO HEADROOM even with an oracle (2026-09-06, AGENT D)

`s19/d_headroom.py`, `s19/results/D_P5_headroom/`, **n = 126, complete**. Instrument is
`s17/refine.py`'s `refine_full` exactly; **only the target vector changes between arms and `sd` is
passed bit-identically to every arm.** Validity gate: `raw` reproduces `s17` refine_full at
max |delta| = **0.0002**, 126/126 within 0.01 A.

| arm | RMSD | median | vs raw [95% CI] | W/L | obj |
|---|---|---|---|---|---|
| avg (start) | **3.048** | 2.837 | | | |
| raw (DEPLOYED) | 3.610 | 3.370 | +0.000 | -- | 56.5 |
| med (Bayes L1) | 3.627 | 3.535 | +0.018 [-0.030, +0.073] | 68/58 | 66.7 |
| mode1 (DENSITY) | 3.644 | 3.431 | +0.035 [-0.017, +0.087] | 54/72 | 72.4 |
| **ORACLE_bestmode** | **3.472** | 3.299 | **-0.138 [-0.213, -0.058]** | 83/43 | 67.5 |
| minofN_null | 3.525 | 3.339 | -0.085 [-0.160, -0.009] | 80/46 | 68.3 |
| rand_flip | 3.629 | 3.503 | +0.019 [-0.052, +0.091] | 58/68 | 78.7 |
| ORACLE_true | **1.148** | 0.307 | -2.462 [-2.770, -2.146] | 117/9 | 3.0 |

> ### THE SENTENCE THAT SETTLES IT, and it needs no control at all.
> **A PERFECT oracle mode-picker lands at 3.472 A. The structure the pipeline already builds is
> 3.048 A. Even with the native in hand to choose every mode, the refinement is +0.424 A WORSE than
> doing nothing.** The branch has no headroom toward the mission whether or not the mode information
> is real.

**THE MODE-DEFINITION DEFECT WAS REAL AND NOT LOAD-BEARING.** D's `mode1` is the **density** mode,
the coordinator's was the **mass** mode. Density: **+0.035 [-0.017, +0.087]**. Mass: **+0.037
[-0.009, +0.088]**. Two definitions differing on 14% of pairs by +2.72 A each, and the outcome moves
by **0.002 A**. The fix was correct and the conclusion never depended on it.

**THE PRE-REGISTERED KILL SHOT, reported exactly as written rather than as it reads best.** F-D5
(*"if `ORACLE_bestmode - raw` does not reach -0.084 with a CI excluding zero, no mode-selection scheme
can move this metric"*) **does NOT fire**: -0.138 [-0.213, -0.058], 83W/43L. There is technically
headroom. **F-D5b, the control clause, splits**: against `rand_flip` as literally pre-registered it
does not fire, but against the tighter `minofN_null` — matched on offset magnitudes **and** the
min-over-k selection — the mode information is **-0.053 [-0.153, +0.044]: spans zero, below MDE, NOT
MEASURED.** AGENT D reported both, named which clause it had written, and argued `minofN_null` is the
correct control because the oracle arm is a minimum over k while `rand_flip` is a single draw.
**Reporting the clause that reads worse, and saying why, is the standard this programme runs on.**

**TWO CROSS-CHECKS WORTH KEEPING.** `ORACLE_true` reaches **1.148 A at objective 3.0**, independently
reproducing Sprint 18's 1.152 basin and G2 (*"perfect distances make this functional point at the
native"*) **on a separately written harness**. And `med` — the Bayes-L1 median, which by AGENT D's
EXACT convexity result is the aim point the **selection** path actually uses — is **+0.018, worth
nothing**. That closes the *"consume the distribution the way selection does"* variant, which neither
lane had tested directly.

**THE BRANCH IS CLOSED ON FOUR INSTRUMENTS**: the coordinator's structural test (L2), D's structural
test (this), D's per-pair aim-point test with a min-of-N null, and D's extractability test (L4).
Final disposition: **the premise stands** (moment collapse is real; 0.641 of mass within 1 A of the
reported mean), **the mechanism is refuted** (coherent conformer flips), and **the information is
real, unextractable out-of-sample, and worth 0.053 +/- 0.10 A structurally even with an oracle.**

---

## L8 — SEEDING DETERMINISM VERIFIED DIRECTLY, and a self-replicate correctly declined (2026-09-06, AGENT D)

AGENT D stopped its own P5 replicate at 70/126 on compute grounds — the machine was at **91% CPU,
7 python processes, 1.88 GB free**, and the replicate was its **second** heavy process against
BRIEF section 15's one-per-workstream rule. It removed `D_P5_headroom_run2/` entirely so no partial
directory could later be mistaken for a result; canonical `D_P5_headroom/` is untouched and complete.

**And it substituted a strictly better check.** A self-replicate of a deterministic pipeline tests
exactly one thing — whether `stable_rng` yields the same stream in a fresh process — because the fits
are deterministic L-BFGS from a fixed start. So that one thing was tested directly, in a fresh
process, at n = 126:

    630 arrays compared (126 targets x 5 seeded/derived arms)   0 mismatches
    checksum over all seeded arms                               176466.6545295715
    distinct streams per target                                 True

**Zero mismatches, and the streams genuinely differ per target** — a constant stream would itself
have been the defect. Four seconds instead of twenty minutes of contended CPU.

**THE STRONGER REPLICATION WAS ALREADY ON THE BOARD, and it is cross-implementation rather than
self**: D's `raw` arm reproduces `s17/refine.py`'s `refine_full` at **max |delta| = 0.0002 A, 126/126
within 0.01** — two independently written harnesses, different modules, different authors, the same
126 numbers. **A self-replicate cannot match that as evidence**, and the coordinator declined the
literal replicate on those grounds rather than spending a core three other lanes were using.

**Coordinator's own compute check at the same moment**: 1.55 GB free of 15.60 GB, **5 python
processes**, only **0.60 GB of it python** — so the memory pressure is not the science jobs, and the
process count is exactly the contract (four lanes plus the coordinator's arm).

---

## L9 — AGENT D's LANE IS COMPLETE, and its most durable output is a methodological rule (2026-09-06)

`s19/agentD_FINDINGS.md`, `s19/PREREG_D.md`, six artefact directories under `s19/results/`, **every
one carrying a COMPLETE flag at n = 126**. `D_T3_tangent/` additionally carries a `SUPERSEDED` file
naming its defect, what in it remains valid (contrasts against `fin_random`, which shared the
inflated subspace) and what does not (any absolute fraction). **Sealed benchmark not read, probed or
derived from; `results/benchmark_manifest.json` never opened.**

### THE RULE THAT SHOULD OUTLIVE THIS SPRINT

> ### A control must be matched in the space the OPERATOR actually works in.

Distinct from "controls are mandatory", which this programme already does well, and distinct from
"a zero-information control must be plausible, not uniform". **Three instances in two sprints, every
one caught by the adversarial lane and none self-caught:**

| # | instance | what the mismatched control actually priced |
|---|---|---|
| 1 | S18 `objceil.py` line 163 — `shuf_paired` passed `sd[pi]` into the fit | the weight permutation, not the error assignment; **0.537 A**, and a lead built on it had already been sent to another lane |
| 2 | S19 `modeshuf` — offsets permuted globally | also **orphaned** each offset from its own pair's weight; S18 had already priced that at -0.467 [-0.598, -0.338] |
| 3 | S19 `iso` / `shuffled` — isotropic in RAW distance space | the fit works in **whitened** space; their correct null is 0.581, not 0.330, and the **0.25** gap looked like a finding until measured |

**The operator transforms the space before acting, and a control built in the untransformed space is
not matched to it.** The mismatch is the same order as the effects being hunted — S18 priced the
weighting at +0.153, exactly the size of the effect it was confounding. **The one correctly matched
control in the sprint, `signflip`, is the one that produced its central result.**

**Companion rule**: *when an analytic null and a measured null disagree, the measurement is the null.*
An algebraic null was offered twice and wrong twice — `2n-2` over-counted the rank, and a
`np.linalg.qr` projection returned 2n columns from a rank-deficient matrix — before the empirical
random-direction null settled it. **A review of the derivation would have passed it both times.**

Both rules are now in project memory, and AGENT D recommends the first gets its own line in the next
brief. **Adopted.**

### LEVER (2) — "do not consume the argmin at all", the one nobody is working on

AGENT D's paragraph, recorded because no other lane touched it. Every stage of this pipeline ends in
an argmin or a top-*m*; T3 says the fit's argmin is **exactly where the in-manifold error gets
realised**; and the recorded operator law says the terminal reads the set **mean**, not the set best.
**A pipeline that optimises a point estimate at every joint is mis-matched at every joint.** Three
forms, increasing in disruption: keep the fit but consume its **posterior** (cheap); keep the
ensemble downstream and report a **weighted set**, which is how NMR peptide structures are deposited
and how these targets are actually defined; or **never form a point estimate anywhere**, so the
objective reweights rather than optimises and the only thing ever chosen is a barycentre's weight
vector.

**The honest obstacle, stated by AGENT D rather than hidden**: the frozen metric is a **point**
metric and section 7 forbids changing it. **But an ensemble's own mean is a point structure and can
be scored**, so measuring what the pipeline produces when no stage ever argmins is a legitimate
diagnostic **under the metric as frozen**. That is the one-experiment version.

### WHAT AGENT D DELIBERATELY DID NOT DO

Block P4 (are the modes an ensembling artefact) is left **OPEN — NOT MEASURED, with its
pre-registration standing rather than deleted**, because P2a and P3 made its outcome unable to change
any conclusion. And it did not touch AGENT A's cross-family attribution, which gates lever (1)
entirely.

### LITERATURE, and a second independent reproduction

`arXiv:2609.02113`'s headline **0.623 A is an ORACLE minimum over 7,332,100 native-scored snapshots
on two targets**; its only native-free numbers are **2.90 A** (chignolin) and **5.39 A** (Trp-cage).
**We are not behind it.** More valuable than the comparison: it **independently reproduces "search
saturates, discrimination binds"** with a 2.18-3.26 A selection gap on a completely different
architecture, and `arXiv:2606.21241` **independently reproduces "the objective does not rank the
native"**. Two external reproductions of this programme's central structural findings, from groups
with no contact with it.

### ONE ENGINEERING ITEM FOR THE NEXT SHARED CONTRACT

`s19/d_tangent.py` recomputed `I.project` per target instead of reading AGENT A's cached
`s19/cache/start_<pdb>.npz`, costing ~8 contended minutes. **That start is deterministic and
identical for every arm and every lane in this sprint.** Any future module needing it should read the
cache; it belongs in the shared contract so the next lane does not rediscover it.

---

## L10 — THE CORRECTED RERUN CONFIRMS L2, resolves L5, and finds ONE arm that is real but useless (2026-09-06, coordinator)

`s19/distobj.py` rerun with the density mode, density-based `nll`, matched-weight arms and AGENT D's
stratified control. **n = 126, `complete: true`.**

| arm | RMSD | median | vs moment [95% CI] | W/L | vs avg |
|---|---|---|---|---|---|
| coordinate average (start) | **3.048** | 2.837 | -- | -- | -- |
| moment (DEPLOYED) | 3.610 | 3.370 | +0.000 | -- | +0.561 |
| **mode** (PRIMARY) | 3.644 | 3.431 | **+0.035 [-0.016, +0.086]** | 54/72 | +0.596 |
| nll *(confounded)* | 3.494 | 3.294 | -0.116 [-0.253, +0.026] | 74/52 | +0.446 |
| nllw | 3.487 | 3.302 | -0.122 [-0.274, +0.027] | 74/52 | +0.439 |
| risk *(confounded)* | 3.600 | 3.372 | -0.010 [-0.078, +0.059] | 67/59 | +0.552 |
| **riskw** | **3.476** | 3.345 | **-0.134 [-0.229, -0.046]** | **77/49** | +0.428 |
| unimodal *(CONTROL, now degraded)* | 3.736 | 3.710 | +0.126 [-0.036, +0.290] | 46/61 | +0.688 |
| modeshuf *(CONTROL)* | 3.689 | 3.455 | +0.079 [+0.012, +0.149] | 47/79 | +0.641 |
| modeshuf_strat *(CONTROL)* | 3.673 | 3.599 | +0.063 [+0.002, +0.126] | 51/75 | +0.625 |
| matched *(CONTROL)* | 4.344 | 4.231 | +0.734 [+0.527, +0.936] | 33/93 | +1.295 |

**1. L2 IS CONFIRMED, AND THE DEFECT WAS NOT LOAD-BEARING.** Density mode **+0.035 [-0.016, +0.086]**
against the superseded mass mode's +0.037, and against AGENT D's independent density-mode
measurement of **+0.035** on a separately written harness. Three numbers within 0.002 A. No
concentration on multimodal targets (high half +0.036, low half +0.033). **The branch stays closed;
nothing reversed.**

**2. L5 IS RESOLVED, AND NEITHER SIDE OF THE ARGUMENT WAS RIGHT.** AGENT D's stratified decomposition:

    position, with alignment held   (strat - mode)     +0.029 [-0.031, +0.088]   W/L 57/69
    the ORPHANING tax               (shuf - strat)     +0.016 [-0.064, +0.096]   W/L 63/63
    both together                   (shuf - mode)      +0.045 [-0.032, +0.125]   W/L 54/72

**Both components are below the 0.084 A MDE and both CIs span zero. NOT MEASURED.** The coordinator's
original claim (+0.096 A of positional information) is **not supported**; and AGENT D's scaling
argument — that orphaning alone could account for the whole penalty, predicting ~0.17 A — is
**also wrong**, by an order of magnitude: the measured orphaning tax is **+0.016**. The honest
statement is that **neither position nor orphaning is measurable here, because both are small.**
*The stratified control was still worth building: it converted one confounded number into two clean
nulls, which is what it was for.*

**3. ONE ARM IS REAL AND STILL USELESS, and both halves must be reported together.** `riskw` — the
Bayes-risk loss carrying the deployed `1/sd^2` weighting — is **-0.134 [-0.229, -0.046] against the
deployed objective, CI excluding zero, above the MDE, 77W/49L.** That is a genuine improvement to the
objective, and it is the **only** arm in either of the coordinator's sprints to beat the deployed
functional on a matched comparison.

> **And it does not matter. `riskw` lands at 3.476 A against the 3.048 A the pipeline already builds
> — still +0.428 A worse than not refining at all.** A better way of doing something harmful.

This is consistent with AGENT D's P5, where the Bayes **median** used as a *target* inside the
squared loss was worth +0.018. The difference is the **loss shape**: `riskw` minimises the risk
itself, an L1-like objective that is robust to the outlying pairs a squared loss chases. Recorded as
a real finding about the objective and an irrelevant one for the mission.

**4. A DETECTOR DEFECT THE COORDINATOR INTRODUCED WHILE FIXING ANOTHER ONE.** Padding the peak
comparison to include the first and last bins — added to fix AGENT D's point (d) — **over-counts
badly**: the multimodal fraction went **0.233 -> 0.857**, against AGENT D's prominence-gated 0.097.
A boundary bin is now scored as a peak whenever it is a local maximum at the edge, including when it
is merely the end of a monotone tail.

**Consequences, stated rather than hidden**: the `unimodal` control now drops 86% of pairs and is no
longer the control it was (its +0.126 is not comparable to the superseded run's +0.018), and the
high/low multimodal split is near-meaningless at this detector setting. **The PRIMARY is unaffected**
— `mode` uses the density argmax and never consults the peak count — which is why the branch's
disposition does not move. **The correct detector is AGENT D's prominence gate at ~0.10, and any
future modality-conditioned arm must use it, not this one.**

**Coordinator's error count for the sprint: five, all caught by workstreams except this one, which
was caught by reading its own output table.**


---

## L11 — THE FORK RESOLVES: the harmful coherent component is SHARED, and two thirds of it is reproduced by ZERO-INFORMATION references (2026-09-06, AGENT A)

`s19/a_source.py`, `s19/results/a_source.json`, **n = 126, COMPLETE**. All arms ORACLE.

**The decomposition is an identity, not a finding**, and is labelled so: each family's error splits
exactly using that family's own terminal fit `X_F` — `r_F = dhat_F - dtrue = r_coh_F + r_inc_F`,
where `r_coh_F = d(X_F) - dtrue` is the part some real structure realises and `r_inc_F` is the part
none can. Cross-family alignment is then read on `r_coh` specifically.

| pair | raw | **COHERENT** | incoherent | RMSD(X_A, X_B) |
|---|---|---|---|---|
| d20 ~ d20_s1 — **same arch, different seed: THE CEILING** | 0.865 | **0.898** | 0.400 | 1.763 |
| deployed ~ d20 — same arch, different regularisation | 0.737 | 0.833 | 0.237 | 2.368 |
| d20 ~ pairnet — **different architecture** | 0.732 | 0.754 | 0.115 | 2.609 |
| deployed ~ pairnet — different architecture, same inputs | 0.630 | 0.731 | 0.097 | 2.690 |
| deployed ~ poolmean — **not a network at all** | 0.709 | **0.784** | 0.090 | 2.522 |
| deployed ~ sepprior — **ZERO-INFO, sequence-blind** | 0.574 | **0.575** | -0.001 | 3.936 |
| deployed ~ helix — **ZERO-INFO, constant alpha-helix** | 0.507 | **0.598** | 0.017 | 3.724 |

As a share of the same-architecture ceiling (0.898): **cross-architecture 0.81, retrieval 0.87,
sequence-blind separation prior 0.64, constant alpha-helix 0.67.**

**And the sharing is specifically in the part that hurts** — exactly the check that was asked for.
The **incoherent** component does **not** share: 0.09-0.19 across families, **0.00** against the
zero-information references.

### Surgery, magnitude-matched (whitening bound: signflip -1.257)

    minus_shared_m  (PairNet-aligned component removed)   2.830   -0.780 [-0.960,-0.612]  102W/24L
    minus_pool_m    (retrieval-aligned removed)           2.623   -0.987 [-1.177,-0.801]  109W/17L
    minus_sep_m     (zero-info-aligned removed)           2.858   -0.752 [-0.934,-0.578]   96W/30L
    only_shared_m   (ONLY the shared component)           3.842   +0.232 [+0.114,+0.360]   55W/71L
    only_pool_m     (ONLY the retrieval-aligned)          4.035   +0.425 [+0.305,+0.550]   39W/87L
    only_sep_m      (ONLY the zero-info-aligned)          3.961   +0.351 [+0.200,+0.502]   46W/80L

Removing the cross-family-shared component recovers **62-79%** of the whitening bound. **Keeping only
the shared component, at matched magnitude, is WORSE than the predictor's actual error** — all three,
every CI excluding zero. **The component the distogram shares with the retrieval pool is the single
most harmful**, and the retrieval pool is both what generates the candidates and where the training
fragments come from.

### THE VERDICT ON THE FORK: **SHARED. Not architectural.**

> **Every predictor, conditioned or not, emits a TYPICAL peptide of that length, and the harmful
> coherent error is the systematic difference between "typical" and this particular native.**
> That is an identifiability / prior-mean statement, not an architecture statement.

**AGENT A refused to overstate it, and the caveat belongs beside the verdict**: sequence conditioning
is **not** worthless — 0.833 (same arch) and 0.784 (retrieval) sit well above 0.598 (helix), so
conditioning does move the error. **But the majority of the harmful coherent component is present in
a reference measure with no sequence information at all, and no loss, architecture or ensemble
reaches that.**

### CORROBORATION FROM THE PREDICTOR ZOO (`s19/a_models.py`, n = 126, COMPLETE)

Eight leave-fold-out predictors, and **not one moves Ca-RMSD past the 0.084 A MDE with a CI excluding
zero**: deployed 3.610, d20 3.574, d20_s1 3.608, sw_none 3.508, sw_lin 3.689, pairnet 3.529,
d20_ens 3.592, combo 3.516. The best, `combo`, is **-0.094 [-0.191, +0.001]** — CI touches zero,
**NOT MEASURED**. Meanwhile **ORACLE MAE falls from 2.347 to 2.073, an 11.7% cut, for nothing.**
**"Better matrix, worse ranking" reproduces on the objective path, eight for eight.**

### TWO CLOSURES

**1. JOINT CONSISTENCY IS CLOSED.** PairNet's triangle update does exactly what it claims
geometrically — triangle violations **4.09% -> 0.68%**, EDM defect **0.286 -> 0.148**, kappa
**0.809 -> 0.914** (the highest of any model) — and buys **-0.080 [-0.242, +0.082]**.
**Making the field more realisable makes it more COHERENTLY WRONG.** This confirms L3's instruction
at the *architecture* level and retires the joint-consistency candidate entirely.

**2. ONE POSITIVE RESULT, AND IT IS A NEGATIVE.** `sw_lin` and `sw_none` share architecture, data and
regularisation and differ **only** in the training separation weight (`min(sep/3, 6)`, up to 6x
toward long-range):

> **sw_lin - sw_none = +0.181 [+0.081, +0.285], 49W/77L, 4/5 folds — while ORACLE MAE moves +0.010
> [-0.039, +0.062].** Weighting the training loss toward long-range pairs costs **0.181 A of RMSD and
> changes MAE by nothing.**

**The Sprint-18 compass reproduced at the PREDICTOR, in the anti-utility direction** — and the
numerically best single MLP on the board, `sw_none` at 3.508, is the one whose loss is **uniform**.

### BLOCKED, NOT DECLINED

The pre-registered positive half of that ladder — retraining with a **short-favouring** loss — is
blocked on memory: **6,429 of 6,788 training sequences are absent from the hot ESM cache**, so it
needs the 1.5 GB bank against 1.90 GB free. Per BRIEF section 15 and the recorded
`machine-fits-two-heavy-jobs` incident, AGENT A declined to launch it. **Recorded OPEN with the
reason, pre-registration standing unedited.**

**Standing caveat, unchanged**: 2.353-2.368 A is an **ORACLE ceiling, not a method**. There is still
**no native-free estimator of coherence**.

---

## L12 — WHY SCORE GATES DAMAGE AN AVERAGED SET: it is ALIGNMENT, not diversity — and three lanes now converge (2026-09-06, AGENT C)

`s19/PREREG_C.md` written before any number. Artefacts in `s19/results/`, all `complete` with row
counts and config hashes.

### The coordinator's briefed hypothesis is REFUTED, by a control AGENT C built to test it

The brief pointed this lane at the **survivors' error covariance** — the Krogh-Vedelsby reading that
a score gate damages an averaged set by destroying diversity. The accounting says exactly that:
Legacy's gate loses `-delta(D^2) = +0.692 A^2 [+0.581, +0.843]`, 14W/112L, **carrying 148% of the
damage.**

**And it is not the cause.** `rand_lowD` — the lowest-`D` of 200 matched-random subsets, containing
**no score, no physics, no sequence** — reaches **D = 1.723, Legacy's D to three decimals**, and
costs **+0.003 A [-0.003, +0.013]** where Legacy costs **+0.063 A [+0.010, +0.124]**. The pure
diversity contrast `rand_lowD - rand_highD` is a null.

> **An accounting identity absorbed the damage; it did not cause it.** Diversity is where the harm
> shows up in the algebra, not where it comes from. **This is the section-10 trap caught by a
> control rather than by review**, and the coordinator briefed the wrong mechanism.

### The mechanism that survives — ALIGNMENT

In a fixed common frame, **exactly**: `readout^2 = ||b_pool||^2 + 2||b_pool||*ALIGN + ||delta||^2`.
Gate damage is a **magnitude** term plus an **alignment** term and nothing else. Pooled, n = 126:

    physics SCORE gates      ALIGN  +0.0209 [+0.0004, +0.0402]
    SCORE-FREE gates         ALIGN  -0.0048 [-0.0124, +0.0020]
    difference                      -0.0256 [-0.0455, -0.0008]

Cross-arm r = +0.838; within-arm per-target r = **+0.951 / +0.967 / +0.955**. Both pooled statistics
pass a drop-top-versus-uniform-effect concentration check.

> ### Score gates push the emitted mean ALONG the direction the pool is already wrong. Score-free gates, at any diversity, do not.

**The zero-information gate is the WORST in the table**: a constant ideal alpha-helix costs
**+0.1253 [+0.0563, +0.1819]**, **twice Legacy's**. **Physics is not being punished for being
physics — any score ordering is.**

### THE CONVERGENCE — three lanes, three instruments, one statement

This is the sprint's synthesis and no single lane could have produced it:

| lane | finding | the shared direction |
|---|---|---|
| **A** (L11) | the harmful coherent error is **shared across all predictor families**, and the **retrieval-pool-aligned** component is the most harmful of all (`only_pool_m` +0.425) | the pool's own systematic error |
| **D** (L6) | the fit **cannot suppress** the in-manifold component; real 0.695 vs a sign-randomised 0.397 | the same direction, in the fit's basis |
| **C** (this) | **score gates move the emitted mean ALONG the pool's error direction**; score-free gates do not | the same direction, at the selection stage |

**The pool has a systematic error direction; the predictor reproduces it; the fit cannot suppress it;
and any score-ordered gate amplifies it.** Every stage of the pipeline is aligned with the same
error, and the retrieval pool is simultaneously the candidate generator, the training-fragment source
and — per AGENT A — the most harmful shared component.

### The design rule works, and is worth zero

`legacy_clust` (diversity-preserving gate design) beats plain Legacy by **-0.0664 [-0.1392, -0.0087],
75W/51L** — recovering the entire Sprint-18 deficit — **and then ties matched-random**
(-0.0033 [-0.0198, +0.0112], **NOT MEASURED**). The native-free half is real: ranking gate *designs*
by `||delta||` predicts damage at **rho = +0.930**.

Sprint 18 reproduces to the third decimal on all five scores **with a 10x tighter null**, as do the
operator-law misses (+0.0937 Legacy, +0.0790 `leg_torsion`, +0.0098 random).

### Q3 — CLOSED, with a sign

Impossible candidates are **real**: mean 2.66 per 75 below a 2.0 A contact, worst pool contact
0.104 A. At the pre-registered r = 10, Legacy `steric` is **worse than matched-random** (+0.0175
[+0.0048, +0.0348]), **worse than a diversity-preserving rejection** (+0.0170 [+0.0091, +0.0286]) and
**worse than not rejecting at all** (+0.0204 [+0.0081, +0.0376]). `minheavy` fails identically.
Dose-response monotone from r = 2. **F-C3 fires.** The last unrefuted Legacy role — hard rejection of
structurally impossible candidates — is closed.

### Q2 — NOT DELIVERED, and reported as a non-delivery

**No Q2 number is quoted at any n**; two partials are `_PARTIAL_*` named and unread. Gates GC2 and
GC3 pass. **F-C2 is UNTESTED and AMBER's standing role is unchanged.** Two blockers:

1. **A DEFECT IN SHARED MACHINERY, found and fixed rather than worked around.** `core.amber`'s
   minimisation is **unbounded (`steps=0`) and does not terminate on some inputs** — over 40 minutes
   at k=3 on 1A1P and over 35 minutes inside the k>=10 ladder on 1D0W, against ~10 s at the incumbent
   k=30. AGENT C added a uniform iteration bound and **certified it bit-inert by gate GC3 at
   0.00e+00 A**. This affects every future sprint that touches `core.amber`.
2. **Memory**: the box sat at 0.6-1.9 GB free against `core.amber.memory_guard`'s 92% ceiling.
   `s19/run_agentC_q2.sh` is queued and resumes per target (~3.5 h on a free box).

### One more self-caught error

Gate GC0's first draft produced a **negative FRAME^2** by mixing free-superposition window RMSDs with
common-frame rebuild quantities — **the exact error class the brief names, and the same one that
broke the Sprint-17 ambiguity identity by 18%.** Caught by the gate, reported rather than silently
repaired.

---

## L13 — A GATE THAT PASSED VACUOUSLY, and the AMBER defect is PATHOLOGICAL rather than merely unbounded (2026-09-06, AGENT C, self-caught)

**AGENT C's first iteration bound was wrong in a way its own gate could not catch.** It set
`steps = 10000`, reasoning from the **708-847 force evaluations** a normal call needs — but
`maxIterations` counts **L-BFGS iterations**, each running a line search over several evaluations.
Measured: a capped k = 10 call on 1CEK still exceeded **seven minutes**
(`_DIAG_1CEK_cap10000_gt7min.log`).

> ### GC3 passed at 10000 because the bound was inert — and it was inert because it was NEVER REACHED.
>
> **A certification that the bound changes nothing is worthless if the bound is never hit.** This is a
> distinct error class from the four in Z1: not a control matched in the wrong space, but a **gate
> that passes vacuously**. Any "this intervention is bit-inert" certificate must also show the
> intervention was **exercised**.

**AND THE DEFECT IS WORSE THAN L12 RECORDED IT.** `core.amber`'s minimisations are not merely
unbounded — they are **pathological**: **more than 10,000 L-BFGS iterations on a thirteen-residue
peptide**, where the incumbent k = 30 on the *same* peptide converges in **4.9 s**. **It is the
SOFTNESS OF THE RESTRAINT, not the target, that triggers it.** That is a property of the operator
this programme deploys, and it has been latent under every AMBER arm in five sprints.

Deployed bound is now **2000**, re-certified by GC3 at **0.00e+00 A / 0.00e+00 kcal/mol** on 5
targets, with `hit_cap` flagged and counted per arm. The cap-10000 partial is preserved as
`_PARTIAL_agentC_pareto_cap10000_n3.json` and **unread**.

### A synthesis contribution from AGENT C, and it unifies its own two results with AGENT A's

C's open item — *"why is a physics score correlated with the pool's error direction at all?"* — has an
answer once AGENT A's L11 is in hand:

> **Any score that prefers compact, well-formed, pool-typical geometry is selecting TOWARD the pool's
> own systematic error.**

That explains **C2** (score gates raise ALIGN, score-free gates do not) and **C4** (the
zero-information constant-alpha-helix gate is the *worst* arm, at twice Legacy's damage) in one
statement — **because a constant helix is the purest possible "prefer pool-typical compact geometry"
rule**, carrying no target information at all and therefore nothing but the bias. Label
**PLAUSIBLE and testable**; it belongs to whichever lane owns the retrieval pool.

**Q2 status**: the n = 126 run is going under the certified bound with a blocking waiter on the
completion flag. **If it lands, both axes with the exclusion count and the rotated-frame maximum; if
not, section 4 stands as written and F-C2 remains UNTESTED.** Q1 and Q3 are final.

**L13 ADDENDUM — AGENT C applied Z6 against its own fix.** GC3 at the deployed `steps = 2000` is
**also never reached** on its five test targets, because the incumbent k = 30 converges there in
3-9 s. That is precisely what GC3 is *for* — the bound must not perturb the incumbent **on inputs
where the incumbent converges** — **but it certifies nothing about an arm that does hit the bound.**
An arm that hits the cap is a **different operator**: *"2000 iterations of restrained relaxation"*,
not *"restrained relaxation"*, and is a legitimate Pareto point only under that label. `hit_cap` is
counted per arm and reported beside both axes.

> **General form, adopted programme-wide: report how many times the bound, guard, threshold or
> fallback actually FIRED. A pass with zero firings is not evidence.**

**AND THE PROVENANCE OF `rand_lowD` MATTERS MORE THAN THE RESULT.** AGENT C states plainly that the
control which refuted the coordinator's briefed diversity mechanism was **not foresight**: falsifier
F-C1 required it to name in advance what would refute its own hypothesis, and *"a score-free gate at
matched D"* is the only honest answer to that question. **C says it would have been very willing to
believe the diversity reading when section 3.4 returned `S_D = 1.48`.**

> **The pre-registration did the work, not the judgement.** This is the strongest single argument in
> the programme's record for pre-registering the FALSIFIER rather than the hypothesis, and it should
> be quoted in the next brief.


---

## L14 — AGENT A CLOSES: the harmful mode is NAMED, and it is NOT ESTIMABLE NATIVE-FREE. This is an identifiability wall. (2026-09-06)

Seven arms, all **n = 126, all `.COMPLETE`**. `s19/agentA_FINDINGS.md`; `s19/PREREG_A.md` unedited.

### 1. Which mode owns the harm — and it is not the one that explains the most variance

`s19/a_struct.py`, ORACLE, **magnitude-matched twins** so an arm cannot earn credit merely for
carrying less error. Reference gap `real - shuffled` = +1.000 A.

| mode removed (magnitude-matched) | RMSD | vs real | share of gap |
|---|---|---|---|
| per-target **OFFSET** | 3.784 | +0.174 [+0.092, +0.265] | **−17.4% — removing it HURTS** |
| per-target **STRETCH** | 3.679 | +0.069 [+0.006, +0.135] | **−6.9% — removing it HURTS** |
| **per-target SEPARATION PROFILE** | **3.085** | **−0.525 [−0.721, −0.343]** | **+52.5% ← the mode** |
| per-residue ADDITIVE | 3.536 | −0.073 [−0.181, +0.036] | +7.3% — NOT MEASURED |

**Variance explained**: offset 0.141, stretch 0.226, separation profile 0.361, per-residue additive
**0.430**.

> **The per-residue mode explains the MOST variance and owns NONE of the harm.**
> **Variance explained does not price damage.** Anyone decomposing a residual and reading the
> largest term is reading the wrong number.

The coherence is a real, **local, share-a-residue** effect: standardised-residual correlation
**+0.087** for pairs sharing a residue, **−0.061** for pairs sharing none, **−0.017** permutation
null. And the offset/stretch signs kill the one component that would have been trivially estimable:
**a global size error is not the problem.**

### 2. Is it estimable native-free? **NO — and the reason is mechanical**

`s19/a_sepfix.py`, n = 126. Five numbers per target is small enough that an estimator might exist;
each arm replaces the ORACLE `mean_b(dtrue)` with a native-free estimate. Only `ORACLE_no_sep` reads
the native.

| arm | RMSD | vs real | W/L | folds |
|---|---|---|---|---|
| *incumbent (coordinate average)* | **3.048** | | | |
| real (deployed) | 3.610 | — | | |
| **ORACLE_no_sep — THE CEILING** | 3.085 | −0.525 [−0.721, −0.343] | 84/42 | 5/5 |
| poolfull (replace `dhat` wholesale) | 3.318 | −0.292 [−0.450, −0.131] | 80/46 | 5/5 |
| poolmed_sep | 3.455 | −0.155 [−0.266, −0.048] | 81/45 | 5/5 |
| pool_sep | 3.512 | −0.098 [−0.195, −0.004] | 71/55 | 4/5 |
| pool_sep_half | 3.554 | −0.056 [−0.120, +0.008] | 77/49 | 4/5 |
| rand_sep *(MATCHED-RANDOM)* | 3.865 | +0.256 [+0.145, +0.372] | 49/77 | 5/5 |
| helix_sep *(ZERO-INFORMATION)* | 3.878 | +0.268 [+0.077, +0.466] | 59/67 | 4/5 |
| global_sep *(ZERO-INFORMATION)* | 4.349 | +0.739 [+0.537, +0.942] | 32/94 | 5/5 |

`pool_sep` **clears every control the brief required**: −0.354 [−0.493, −0.214] against
matched-random and −0.366 [−0.536, −0.199] against the constant-helix zero-information control, both
5/5 folds.

**And it FAILS the trap AGENT A wrote into the module header BEFORE running it**:
`pool_sep − poolfull = +0.193 [+0.070, +0.324]`, 44W/82L. The whole ladder
3.610 → 3.554 → 3.512 → 3.455 → 3.318 is **monotone interpolation toward `poolfull`.** It is *"move
toward the retrieval pool"*, not *"estimate the coherent mode"* — a simpler classical baseline
reproducing the effect, which BRIEF §12 requires be treated as a downgrade. **DOWNGRADED.** Every arm
on the board, `poolfull` included, is **worse than the 3.048 A incumbent**; the native-free route
recovers 19–30% of the 0.525 A ceiling and only by regressing toward what the pipeline already does.

> ### THE WALL, STATED MECHANICALLY
>
> The best native-free estimator of this target's separation profile available anywhere in the
> project **is the retrieval pool** — and L11 measured that the pool **shares 0.87 of the harmful
> component being estimated.**
>
> **The estimator is made of the bias.**

That is not a gap in effort. It is a structural reason why no native-free estimator of coherence has
been found in two sprints of looking.

### 3. A2 REFUTED — the harmful error is not in a detectable subset

`s19/a_subset.py`, n = 126: repair 25% of pairs and ask whether any native-free criterion beats a
size-matched **random** choice, with a **magnitude-matched** column so an arm cannot win by picking
big residuals.

    nf_model (LFO error model)   3.348   +0.232 [+0.103, +0.356]      long_sep    3.210  +0.093 [-0.053,+0.232]
    hi_sd                        3.436   +0.320 [+0.220, +0.423]      short_sep   3.286  +0.170 [+0.053,+0.291]
    lo_sd                        3.277   +0.161 [+0.034, +0.285]      terminal    3.182  +0.066 [-0.042,+0.181]
    multimode                    3.273   +0.157 [+0.054, +0.259]      random      3.116   ---
    ORACLE_absr (perfect)        2.912  -0.205 [-0.364, -0.041]

**Not one native-free criterion beats random**; all seven are worse, five with CIs excluding zero.
And **ORACLE-perfect detection of the largest 25% of errors is worth only −0.205 A** once magnitude
is held fixed. A's own error model ranks the ORACLE error at Spearman **+0.537** and converts to
nothing. **Detecting the error is partly possible; repairing what you detect is worth nothing.**

### Two things AGENT A protected the record against, unprompted

1. **A magnitude confound in its own raw column.** Unmatched, `long_sep` reads −0.264 and `nf_model`
   −0.163, i.e. "beating" random — purely because those subsets carry smaller residuals (resid RMS
   2.14 vs random's 2.79). **Anyone reading that column alone would report a detectable subset.**
   Both columns were run for exactly this reason.
2. **A tension with Sprint 18's G2e that it refused to smooth.** A's matched column has `short_sep`
   +0.170 and `lo_sd` +0.161 (random wins) where S18's G2e has `fix_short` −0.294 and `fix_confident`
   −0.308. **These are different operations, not a contradiction**: G2e *shrinks* residuals within a
   class keeping every pair; A's arm *zeroes* a class and rescales the complement, concentrating the
   surviving residual into a structured subset and thereby **increasing** its coherence. Both agree
   with the sign-coherence mechanism. **They must not be quoted against each other**, and A wrote
   that into the findings rather than leaving the numbers side by side.

### Final disposition of AGENT A's lane

    ESTABLISHED   A3 coherence IS the mechanism (signflip_exact -1.202 = 128% of the gap; coherent +0.139)
    ESTABLISHED   A1 the error is anti-correlated with the utility compass -- and is NOT the mechanism
    ESTABLISHED   A5a the anti-utility training loss costs +0.181 [+0.081,+0.285] at unchanged MAE
    REFUTED       A4 unrealisability as a mechanism (anti-correlated with harm)
    REFUTED       A2 harmful error in a detectable subset
    CLOSED        A6 joint consistency (PairNet); A8 conditional weighting law (premise falsified)
    NOT MEASURED  A7 ensembles (combo -0.094 [-0.191, +0.001])
    DOWNGRADED    A9 native-free estimate of the harmful mode
    OPEN          A5b short-favouring retrain -- NOT RUN, blocked on RAM, pre-registration unedited

**Standing caveat, kept verbatim at AGENT A's request: 2.353–2.368 Å is an ORACLE ceiling, not a
method, and must never be quoted as an achieved number.**

---

## L15 — AGENT A CORRECTS ITS OWN HEADLINE: a clip-induced magnitude confound, and it was pointing AGAINST the claim (2026-09-06)

`s19/a_audit.py`, `s19/results/a_audit.json`, **n = 126, COMPLETE**.

**THE DEFECT.** Every field is clipped at 2.0 A before the fit, exactly as `s18/objceil.py` does.
**Flipping a residual's sign can push a distance below 2.0, so the clip SHRINKS the whitened field.**
`signflip` carried residual RMS **3.062** against `real`'s **3.205** — a 4.5% reduction. On the alpha
ladder a 25% magnitude cut is worth 0.504 A, so 4.5% extrapolates to ~0.091 A, roughly 7% of the
headline. **The arm was not at exactly matched magnitude, which is the one property the entire
argument rests on.**

**THE FIX.** `*_exact` rescales the field **after** the clip to a fixed point where realised residual
RMS equals `real`'s to machine precision — **3.2049 for every arm**, not approximately.

| arm | RMSD | resid RMS | vs real | W/L | folds |
|---|---|---|---|---|---|
| real | 3.610 | 3.2049 | — | | |
| signflip *(as published)* | 2.368 | 3.0621 | −1.242 [−1.458, −1.038] | 113/13 | 5/5 |
| **signflip_exact** | **2.408** | **3.2049** | **−1.202 [−1.408, −1.004]** | 112/14 | 5/5 |
| shuffled *(as published)* | 2.649 | 3.0076 | −0.961 [−1.194, −0.747] | 100/26 | 5/5 |
| **shuffled_exact** | 2.671 | 3.2049 | −0.938 [−1.169, −0.717] | 99/27 | 5/5 |
| **iso_exact** | 2.733 | 3.2049 | −0.877 [−1.106, −0.652] | 93/33 | 5/5 |

    signflip_exact - shuffled_exact = -0.264 [-0.448, -0.084], 83W/43L, 5/5 folds

**THE HEADLINE SURVIVES.** The clip was worth **0.040 A — 3% of 1.242**, less than the 7% the
extrapolation predicted, **because the controls were clipped harder than `signflip` was.** Every
qualitative conclusion is unchanged. The ledger now carries **−1.202 and 128%**, and *"at exactly
matched residual RMS"* is **literally** true rather than approximately true — which for this
particular claim is the difference that matters, since matched magnitude is the whole argument.

**AND THE CONFOUND WAS POINTING AGAINST THE CLAIM, NOT FOR IT.** `shuffled` and `iso` were clipped to
3.008 and 3.006 — a **6.2%** reduction — so the published signflip-minus-shuffled gap of −0.282 was
measured with **the controls carrying LESS error than `signflip` did.** Correcting all three arms to
a common magnitude leaves it at −0.264 [−0.448, −0.084], 5/5 folds. **The claim that beat the
Sprint-18 control was, if anything, understated.**

**THIS IS A FIFTH INSTANCE OF Z1**, and a subtle one: *a control must be matched in the space the
operator actually works in* — and **the clip is part of the operator.** Matching magnitude *before*
the clip is not matching it *after*. Every previous instance was a control built in the wrong space;
this one was built in the right space and then silently moved by a downstream operation.

**Disposition**: the published −1.242 figures are marked **superseded rather than deleted**
throughout `s19/agentA_FINDINGS.md`, so the two can be reconciled by anyone reading the earlier
messages, and item 0 of A's "what I would tell the next sprint" says explicitly to quote −1.202.

**A CORRECTION TO THE COORDINATOR'S CREDIT, from AGENT A.** On "the estimator is made of the bias"
being the sprint's closing sentence, A wrote: *the credit belongs to the measurement, not the
sentence — it is only a result because `pool_sep` cleared matched-random and the zero-information
helix at 5/5 folds first. Had it not cleared those, "we could not find an estimator" would have been
the honest reading.* **Recorded as A stated it.**

Eight arms at n = 126, all with `.COMPLETE` flags: `a_topo`, `a_coh`, `a_struct`, `a_models`,
`a_source`, `a_subset`, `a_sepfix`, `a_audit`. A5b remains OPEN, not run, blocked on RAM,
pre-registration unedited.

---

## L16 — LEVER (2), SOFT CONSUMPTION: falsifier FIRES, and the incumbent's cut is validated on two axes (2026-09-06, coordinator)

`s19/softsel.py`, **n = 126, `complete: true`**. Replaces the shipped hard top-75 cut with a
temperature-weighted barycentre over the **entire** K = 500 pool — nothing selected, nothing
discarded. One knob whose limits contain both endpoints (T -> 0 the argmin, T -> infinity the uniform
pool mean). All arms native-free.

| arm | RMSD | median | vs top75 [95% CI] | W/L | ESS |
|---|---|---|---|---|---|
| **top75 (INCUMBENT)** | **3.048** | 2.837 | -- | -- | 75.0 |
| soft0.002 | 3.420 | 3.441 | +0.372 [+0.250, +0.500] | 33/93 | 1.7 |
| soft0.01 | 3.350 | 3.293 | +0.302 [+0.184, +0.418] | 40/86 | 10.2 |
| soft0.05 | 3.160 | 2.959 | +0.112 [+0.036, +0.187] | 58/68 | 47.2 |
| soft0.1 | 3.061 | 2.849 | +0.012 [-0.034, +0.059] | 64/62 | 80.9 |
| **soft0.25 (best)** | **3.020** | 2.817 | **-0.028 [-0.080, +0.024]** | 66/60 | 178.6 |
| soft0.5 | 3.080 | 2.993 | +0.032 [-0.043, +0.104] | 59/67 | 289.2 |
| soft1 | 3.177 | 3.151 | +0.129 [+0.020, +0.232] | 55/71 | 388.0 |
| poolmean *(ZERO-INFO)* | 3.396 | 3.261 | +0.348 [+0.185, +0.519] | 49/77 | 500.0 |
| rand75 *(MATCHED-RANDOM)* | 3.449 | 3.320 | +0.401 [+0.224, +0.580] | 43/83 | 75.0 |
| hard25 | 3.075 | 2.870 | +0.027 [-0.025, +0.078] | 49/77 | 25.0 |
| hard150 | 3.072 | 2.934 | +0.023 [-0.033, +0.079] | 66/60 | 150.0 |
| hard300 | 3.171 | 3.074 | +0.123 [+0.006, +0.242] | 53/73 | 300.0 |

**1. THE PRE-REGISTERED FALSIFIER FIRES.** The best soft arm is **-0.028 [-0.080, +0.024]** — CI
spans zero, magnitude a third of the 0.084 A MDE. **NOT MEASURED, and it must not be quoted as an
improvement**; it is also an argmin over a temperature grid on the tuning instrument. **Soft
consumption in this form does not beat the hard cut.** The coordinator's pre-registered expectation
("I expect this to fail, and I am running it because it is cheap and unexamined rather than because I
believe it") is what happened, for the reason given: a soft weighting is **ordered by the same
score** — a gentler version of the same operation, not a different one.

**2. THE SECOND FALSIFIER CLAUSE DOES *NOT* FIRE, AND THAT IS THE INFORMATIVE HALF.** The best soft
arm beats the zero-information endpoint by **-0.376 [-0.537, -0.223]**. **The distogram scores
contribute real information through this operator** — the null result is not "the scores are
useless", it is "hard and soft consumption of the same scores are equivalent".

**3. THE INCUMBENT'S CUT IS VALIDATED ON BOTH AXES SIMULTANEOUSLY, which was not the point of the
experiment and is its most useful output.** `hard25` +0.027, `hard150` +0.023, `soft0.1` +0.012,
`soft0.25` -0.028 — **every arm within a third of the MDE of the shipped m = 75.** The response is
flat in cut *size* and flat in cut *hardness* around the deployed setting, and degrades in both
directions away from it (hard300 +0.123, soft1 +0.129, soft0.002 +0.372). **m = 75 sits in the middle
of a broad optimum on a two-dimensional axis nobody had swept.**

**4. AND A NATIVE-FREE POSITIVE FOR THE SHIPPED FILTER: it beats discarding the same NUMBER at random
by +0.401 [+0.224, +0.580], 83/43.**

> **This is an important contrast with L12, and the two are consistent rather than in tension.**
> AGENT C found that **physics** score gates lose to matched-random gates and that a constant-helix
> gate is the worst of all. The **distogram** gate beats matched-random by 0.401 A. Both follow from
> Z7: *any score that prefers compact, well-formed, pool-TYPICAL geometry selects toward the pool's
> own error* — and the distogram score is **target-conditioned**, not a typicality rule. **"Any score
> ordering is harmful" would be the wrong reading of L12; "any POOL-TYPICALITY ordering is harmful"
> is the right one.**

**DISPOSITION.** Lever (2) is **REFUTED in this form** — softening the top-*m* cut. It is **not**
refuted in its strong form, which D stated and nobody has run: *never form a point estimate anywhere*,
so that the objective reweights rather than optimises. This arm softened **one** joint out of several
and left the fit's argmin and the terminal readout untouched. Recorded so the distinction is not lost.

---

## L17 — RETRACTION: the "pathological AMBER minimiser" was CPU STARVATION. Z4 is REFUTED. (2026-09-06, AGENT C self-refuted)

**AGENT C reported this before its repair run completed, because it retracts claims already in this
ledger. The error is C's; the failure to verify it is the coordinator's.**

On a quiet box, uncapped (`steps = 0`, the deployed protocol), the exact calls reported as
non-terminating:

    1A1P k=1    14.7 s      1CEK k=1/3/10   ~6 s each      1D0W k=1/3/10  8.7 / 10.0 / 7.1 s
    1A1P k=3    15.2 s      all converged
    1A1P k=10   12.8 s

**Reported as ">40 minutes of a full core" and ">7 minutes". They take ten to fifteen seconds.**
What was measured was **CPU starvation from six other lanes' processes**, attributed to the
minimiser. **`core.amber` has no pathology and no non-termination defect.**

**THE INDEPENDENT CORROBORATION IS STRONGER THAN THE RETEST**, and it was in C's own completed Q2
table: **max per-call wall 25.7 s, p95 13.6 s, 0 calls over 40 s across 1260 minimisations** —
including the arms called pathological. **C caught it only because that contradicted its own
published claim**, which is the check working.

**WHAT COMES OFF**: Z4 and its "sharpened" form; the *"latent under every AMBER arm for five
sprints"* framing; and `FINAL_REPORT.md` §6's shared-machinery-defect paragraph.

**WHAT SURVIVES, and is now better founded.** **Z6 — a gate can pass vacuously — STANDS**, and is
now demonstrated on **real data rather than on a false premise**: in the completed run the
`steps = 2000` bound fired on **1 of 1260 calls**, so GC3's certificate for that run is vacuous in
exactly the sense the rule defines, and §4.7 reports it as **vacuous rather than as a pass**. *The
rule is right; its origin story was wrong.*

**THE Q2 RESULTS ARE UNAFFECTED.** A call converging before the cap is bit-identical to uncapped, and
essentially none reached it; GC3 showed bit-identity at k30 directly. An exactness check on the
single flagged call (8T61) is running.

**A PRE-REGISTRATION DEVIATION MADE ON A FALSE PREMISE, and being repaired rather than explained.**
C dropped k=1 and k=3 from the pre-registered ladder citing non-termination. **That justification is
void**, so both rungs are being rerun uncapped at n=126 (~45 min) **to restore the ladder rather than
leave a gap C created.** That is the right response and it is more expensive than the alternative.

### COORDINATOR ERROR #6, and it is a different kind from the other five

The other five were wrong hypotheses or defective controls. **This one is a failure of
verification**: the claim was amplified into the ledger, the claims register **and the final report**
as a shared-machinery defect affecting every future sprint, **on one lane's report, with no
independent check** — when a timing test costs seconds and the disconfirming evidence was sitting in
that lane's own completed table. BRIEF §11 says *read the code before believing the claim*.
**It applies to performance claims exactly as it applies to scientific ones**, and I did not apply it.

---

## L18 — Q2 DELIVERS: F-C2 FIRES, and a scalar validity score would have promoted a broken structure (2026-09-06, AGENT C)

> **PROVISIONAL — every number in this entry comes from the run now SUPERSEDED as cap-limited (see
> L19). The science is expected to reproduce and will be re-quoted from the clean run; until then no
> figure here is final.**

**126 rows, complete.** The AMBER Pareto is mapped, and the answer is that the two axes genuinely
dissociate.

**Three arms beat the k=30 incumbent on Cα past the MDE** — `caonly_k300` **−0.110 [−0.136, −0.092]**,
`blend75` **−0.119**, `blend50` **−0.092** — **and all three fail the validity axis.** Every arm whose
validity matches k30 is within **±0.024 Å** of it on accuracy.

> ### The sharpest item in the sprint's physics lane
> **`caonly_k300` is 0.110 Å MORE ACCURATE than the incumbent with a clash count statistically
> indistinguishable from it (−0.008 [−0.050, +0.032]). A scalar clash-based validity score would have
> PROMOTED it.**
>
> **What kills it is `cis_frac +0.4247` — 42% more cis peptide bonds — because pinning only Cα lets
> the peptide plane flip.**

That is a concrete, measured argument for **reporting validity as a vector rather than a scalar**,
and it is exactly the failure mode Sprint 18's cis-peptide finding predicted at the *averaging*
operator now reappearing at the *repair* operator.

**Rotated-frame null max |ΔCα| = 0.014192 Å** (reported as the maximum, per the standing rule).
**Convergence exclusions `1D6X 2NB7 7BX2` for the sixth sprint**, rising to **64/126 at k=1000** —
because a pinned backbone cannot relieve its own strain.

**F-C2 FIRES. AMBER's standing role is unchanged: stereochemical repair only.** There is no
Pareto point that buys Cα accuracy without paying for it in validity.

---

## L19 — THE CAP DID BIND, AND THE RUN IS SUPERSEDED RATHER THAN CAVEATED (2026-09-06, AGENT C)

**The exactness check on the single flagged call FAILED**, and AGENT C reported it immediately
because it changes the disposition of a completed 126-target run:

    8T61  k30  capped vs uncapped :  max dCa 0.1595 A   dE 0.418 kcal/mol

**The `steps = 2000` bound genuinely bound on that call.** The wall-clock proxy was right, C's §4.7
draft claim that *"no arm is cap-limited"* was wrong, and **one row of the completed table is a
different operator from the incumbent on that target.** C caught it **only by running the check
rather than arguing the effect was small.**

### The judgement call, and it is the right one

> *"Arguing that 0.16 A on 1 of 123 targets moves a mean by ~0.001 A would be true and would still be
> the wrong call — the confound is removable, so I am removing it rather than bounding it."*

**Combined with L17: the bound was introduced to fix a defect that does not exist, and it perturbed
at least one measurement by 0.16 A.** The capped run is preserved as
`_SUPERSEDED_agentC_pareto_capped2000_n126.json` — complete and valid data, **superseded and not
quoted.**

**AND ALL THREE COMPUTE DEVIATIONS ARE WITHDRAWN**, because every one rested on the refuted timing
claim. The rerun is **the pre-registration as written**: `STEPS = 0` (the deployed protocol, no bound
anywhere), the **full** ladder `k = 1, 3, 10, 30, 100, 300, 1000`, and `caonly` at the full declared
`k = 30, 100, 300` — 14 arms per target, ~2.6 h on the now-quiet box. **Restoring a pre-registration
is more expensive than explaining a deviation, and C chose the expensive option twice in one night.**

**Q1 and Q3 never touched AMBER and are unaffected.**

### A NEW RULE, and it is the one that caught the error

> ### A bound must be verified on the calls it actually BOUND, not only on the ones it did not.

GC3 certified inertness on **five converging calls** and was **silent about the one call that hit the
bound** — and that call was perturbed by **0.16 A**. This is the operational half of Z6: *a gate can
pass vacuously* says the certificate is empty when the intervention never fires; this says that when
it **does** fire, the firing cases are **the only ones the certificate is about.** C notes it was a
two-minute test it nearly skipped.

---

## L20 — "BEST BUILT 3.048 A" IS A CONTRACTED POINT CLOUD, NOT A STRUCTURE. The bar was wrong. (2026-09-07, WORKSTREAM D audit; coordinator's label)

**Reproduced end to end** from `bench_results/cache/1fc9f2dcf489e2fb`, 126/126 records, no benchmark
manifest opened. Both quoted numbers reproduce exactly:

    rmsd_avg  = 3.0483   avg_ca    coordinate average of the shipped top-75
    rmsd_fit  = 3.2041   fit_ca    lambda=0 projection   <- THE INCUMBENT
    rmsd_arm  = 3.2148   ca        lambda=0.3 ramah arm
    rmsd_full = 3.2355   amber_ca  the pipeline's actual final emission

**INDEPENDENTLY VERIFIED BY THE COORDINATOR** (not taken on the lane's report — the Sprint-19 lesson):

    avg_ca   mean virtual Ca-Ca bond  2.961 A   min over all 126 targets  0.649 A
    fit_ca   mean virtual Ca-Ca bond  3.804 A
    physical ideal 3.804 A  ->  avg_ca is contracted 22.2%

> **`avg_ca` is a shrinkage estimator of Ca positions, not a backbone.** No side chain can be placed
> on it and AMBER cannot score it. **3.048 A is the best POINT CLOUD the pipeline emits; 3.204 A is
> the best STRUCTURE it builds** — projecting the point cloud onto the pipeline's own ideal-geometry
> manifold *gives the incumbent*. The 0.156 A gap is the price of stereochemical realisability
> (1.9x MDE).

**IT IS NOT A METRIC EXPLOIT OTHER ARMS COULD COPY.** The best *constant* global rescale of `fit_ca`
buys only **0.034 A** (3.2041 -> 3.1700 at s = 0.94), below MDE. The averaging operator's contraction
is **per-target adaptive shrinkage**, which is real. **So this is a finding about the LABEL and the
BAR, not about the metric being broken**, and the contraction itself was already recorded (S18: 22.4%).

### What it corrects, all basis mismatches

| where | as published | like-for-like |
|---|---|---|
| `FINAL_REPORT` §4 / `CLAIMS` O6 | ORACLE mode-picker 3.472 vs "3.048 the pipeline already builds" = **+0.424** | 3.472 vs **proj 3.213** = **+0.259**. *Direction unchanged; magnitude overstated by 64%.* |
| `CLAIMS` W4 | "every arm, `poolfull` included, is worse than the 3.048 incumbent" | the arms are `I.build_ca` chains; the comparator must be the built start |
| `s20/BRIEF.md` §1 | "best built 3.048 A" | corrected in place |

**WORKSTREAM D FOUND THIS IN ITS OWN LANE'S PRIOR OUTPUT** — `D_P5_headroom` was its predecessor's
arm — and reported it against itself.

**COORDINATOR ERROR #7, and the second of the "propagated without checking" kind.** The label drifted
between Sprint 12 — which had it right (`s12/dir_FINDINGS.md:241`: *"bayes = 3.048 avg / 3.204
projected"*) — and Sprint 19, where I adopted "best built 3.048" and then used it as the bar in a
report, a claims register and a published artifact. **`rmsd-reporting-basis-mismatch` is a recorded
memory in this project and this is its fourth instance.**

**STANDING RULE, adopted into `s20/BRIEF.md` §1: every arm states its basis** — point cloud / built
chain / repaired emission — and a built chain is never compared against the 3.048 point cloud.

**WHAT IS UNAFFECTED.** Every Sprint-19 coherence result is a comparison *between built chains in one
basis* (`signflip`, `shuffled`, `coherent`, the in-tangent fractions, the source-alignment table, the
ALIGN decomposition). **The mechanism, the source and the wall all stand unchanged.** What changes is
the mission ladder and any comparison that crossed bases.

