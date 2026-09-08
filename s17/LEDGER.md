# SPRINT 17 — DECISION AND RETRACTION LEDGER

Written as the sprint runs. Every entry carries the date, who made the call, the evidence at the
time, and what later overturned it. Sprint 15 closed with eleven coordinator retractions and
Sprint 16 with five; the same standard applies.

---

## L1 — The K = 500 "retrieval ceiling" is a truncation artefact (2026-09-06, coordinator)

The programme has planned against a "K = 500 pool oracle of 1.711 A" as if it were the retrieval
ceiling. It is not. `u["W"]` holds **13,000-27,000 windows per target**; `u["order"]` is the
BLOSUM62 ranking over that universe; `I.pool_idx` truncates at 500. A four-target probe moved 1CB3
from 2.028 A at K = 500 to **1.096 A** over its full universe.

Executed as `s17/oracle_map.py`, the sprint's shared instrument.

> **AMENDED AT L19.** The ceilings stand, but **the obvious reading of them is wrong**. AUDIT's
> matched null -- a random permutation of the identical universe -- gives 0.473 A where the
> observed truncation cost is 0.397 A, so **observed minus null is -0.075 [-0.132, -0.022]** and
> the gain is *below* the order-statistic floor. Widening K reveals **more draws**, not better
> structures the ranking was hiding, and the retrieval apparatus buys **+1.5 targets** of 2.0 A
> recall over a random 500. The caveat below was raised at the time and AUDIT answered it in both
> directions: the deep universe is **not** near-duplicates (d* = 1.355 A, a dense continuum), so
> the ceiling is real -- but the ranking is not extracting it.

**A caveat was raised at the same time it was found, and assigned to AUDIT rather than argued
away**: the minimum of N draws falls with N *even when nothing new is present*. Part of the
apparent ceiling gain from 500 to 13,000 windows may be redundancy and noise-mining rather than
genuinely better structures. AUDIT owns the matched null. **No wide-K ceiling may be leaned on
until it reports.**

---

## L2 — The design equation, and the requirement it produces (2026-09-06, coordinator)

`s17/selection_theory.py`. Given a candidate set whose member-RMSD distribution is measured, and a
selector whose ranking has copula correlation rho with the truth, what does the selector return?
Gaussian copula, one dependence parameter, both limits verified numerically (rho = 1 reproduces the
oracle best exactly, rho = 0 reproduces random selection exactly).

At n = 126 this converts "we need better ranking" into a number, and shows that **widening K does
not merely raise the ceiling -- it lowers the ranking skill required**: the rho needed for 2.0 A
falls from *unreachable* at K = 75 to 0.75 over the full universe, and for 2.5 A from 0.80 to 0.60.

---

## L3 — RETRACTED, and it is the coordinator's own error, in the message warning against it (2026-09-06, coordinator; caught by SELECT)

I wrote to the SELECT workstream that *"the programme's recorded cross-target in-band ordering of
0.600 lands at roughly 2.02 A on the full-universe curve."* **That is wrong.**

The recorded 0.600 is a **pairwise ordering ACCURACY** with a null of 0.500. The copula's axis is a
correlation. They are different quantities:

    tau = 2*acc - 1  and, for a Gaussian copula, tau = (2/pi) arcsin(rho)
    => accuracy 0.600 -> rho 0.309, not rho 0.600
    => the requirement rho 0.65 corresponds to accuracy 0.725, not 0.638

On the corrected axis the programme's measured point sits **far above 2.5 A**, not at 2.02 A. The
strategic conclusion -- widening K lowers the required skill -- is unaffected; the position of our
measured point on that axis moves a long way in the **pessimistic** direction.

**This is BRIEF section 5's failure mode exactly, committed by the coordinator inside the message
instructing another workstream to guard against it.** Fixed at source rather than in prose:
`rho_from_accuracy`, `accuracy_from_rho` and `rho_from_spearman` are now named functions, the first
one's docstring records this error so it cannot recur silently, and every requirement is printed on
both axes. Spearman is a third axis again (rho = 2 sin(pi rho_S / 6), ~3% from rho_S).

**The general lesson, which is worth more than the number:** the failure survives being written
down as a warning. Only an independent reader caught it.

---

## L4 — The global rho is the wrong axis for an argmin (2026-09-06, SELECT workstream)

Measured on identical K = 500 candidate sets, n = 126: the shipped distance objective has **global
Spearman 0.568** but **in-band Spearman 0.131** (best + 1.5 A band). Fed to the copula, the two
bracket the truth: global rho predicts **2.739** against a realized **3.454** (+0.715, 106/126 above
the curve); in-band rho predicts **3.722** (-0.268, 51/126 above).

> **The shipped objective's skill is almost entirely separating garbage from plausible, and nearly
> absent where the argmin is actually decided.**

Consequence adopted: requirements are read off the **band curve**, not the global curve.
`s17/selection_theory.py` now conditions on the objective's own top-B shortlist (B = 5, 10, 25, 50,
100) and reports the required in-band rho on both axes.

---

## L5 — REFUTED: the distance objective cannot be re-engineered into a better selector (2026-09-06, SELECT workstream)

58 functionals of the existing distogram -- pair weighting including sequence separation and
uncertainty; Bayes risk, NLL, absolute, z, log and squared residuals; Huber, Tukey, median, trimmed
and CVaR aggregation; per-pair standardisation and rank-normalisation across the candidate set;
contact-vs-distance mixtures; a distance-to-typicality blend; a spectral/Torgerson form --
leave-fold-out variant selection returns **-0.009 A [-0.128, +0.110]** against the shipped
objective. The pre-registered E2 falsifier fired.

**With L4 attached this reads as a mechanism, not a disappointment: you cannot fix an in-band
problem by re-weighting a global signal.** All 58 are functionals of the same distogram, whose
in-band content is 0.131.

**The sprint's target therefore changes**: it does not need a better global ranker, it needs an
**in-band discriminator**, and those are different objects with different feature requirements.

---

## L6 — The in-band problem is INFORMATION-limited, not calibration-limited — but the information is present and mis-routed (2026-09-06, coordinator)

`s17/inband.py`. Take the top-B band of the shipped objective at K = 500, build native-free
per-candidate features, and fit the same linear scorer three ways so the differences read as
ceilings. **Numbers below are a 12-target smoke and are NOT quotable**; the n = 126 run is queued.

| arm | mean | |
|---|---|---|
| band best | 2.187 | ORACLE |
| **best SINGLE feature per target** | **2.234** | **ORACLE** |
| per-target fit, in-sample | 2.404 | ORACLE |
| per-target fit, leave-one-out | 2.569 | ORACLE |
| distance argmin | 2.661 | |
| GLOBAL fit (deployable) | 2.682 | |
| consensus argmin | 2.758 | |
| random in-band | 2.914 | CONTROL |
| Legacy argmin | 2.948 | **worse than random in-band** |

    CALIBRATION headroom (global - per-target LOO)    +0.114
    INFORMATION deficit  (per-target LOO - band best) +0.381

**The information deficit dominates by about 3x**, so per-target calibration of a fixed feature
combination is worth ~0.11 A and is not the route to sub-2 A. Section 14 is re-scoped, not pursued
as the main line.

**But the third row reframes it.** The ORACLE *best single feature per target* reaches 2.234 against
a band best of 2.187 -- it recovers almost the whole band. **The information is present in the
feature set; what is missing is knowing which feature to trust for this target.** That is neither
calibration (rescaling a fixed combination) nor an information deficit -- it is a
**selector-portfolio routing problem**, exactly sprint sections 35 and 36, now with a measured
ORACLE ceiling instead of a speculation.

The asymmetry that makes it plausible: a *fitted* per-target linear combination (2.569 LOO) is much
worse than *picking one feature* per target (2.234). Fitting p weights on B = 25 candidates
overfits; choosing among ~11 single features is a far lower-capacity decision, and low capacity is
what makes a routing rule learnable from target-level descriptors.

Assigned to SELECT with the falsifier fixed in advance: **a routing policy must beat a constant
policy** (always use the globally-best single selector) on held-out folds with an interval
excluding zero.

---

## L7 — Compute discipline (2026-09-06, coordinator)

Four agents plus three coordinator runs pegged the box at 100% CPU. The in-band run was **queued
behind a bounded CPU hold rather than launched into contention**, and no further coordinator process
will start until the map, the design equation and the readout sweep have drained. Memory held at
72.5%, well inside the limit.

---

## L8 — RETRACTED: the "routing ceiling" was a min-of-N statistic below its own noise floor (2026-09-06, coordinator; caught by SELECT)

L6 read the ORACLE *best single feature per target* (2.234 against a band best of 2.187) as evidence
that **the information is present in the feature set and only the routing is missing**, and
re-scoped section 14 on that basis. **That reading is wrong and is withdrawn.**

`mean_t min_v RMSD(t, v)` is a **min-of-N statistic**: it falls with the number of selectors V
whether or not any per-target structure exists. With V draws from a band of B members the expected
minimum sits near the 1/(V+1) quantile of the band -- for V = 11, B = 25 that is roughly the
second-best of 25, which is approximately where 2.234-against-2.187 sits. **It is where pure noise
lands.**

SELECT measured the matched null on its own instrument, n = 126, K = 500, 25 real native-free
signals:

| | |
|---|---|
| band ORACLE best (shipped top-25) | 2.609 |
| ORACLE min over **25 real signals** | 2.752 |
| **ZERO-INFO min over 25 random picks from the same band** | **2.694** |
| real minus zero-info | **+0.059 [+0.020, +0.100]**, fold [+0.028, +0.102], 63W/63L |

**The real per-target minimum LOSES to drawing the same number of candidates at random from the
same band, with an interval excluding zero.** Real selectors are correlated with one another, so
their per-target minimum does not even reach the independent-draw floor. SELECT hit the identical
trap an hour earlier on its own 58-variant sweep, where an apparent 0.972 A of "calibration
headroom" was beaten by a zero-information min over 58 random selections (2.148 against 2.482).

**Fixed at source**: `s17/inband.py` now computes the zero-information min-of-N null alongside the
routing arm (V matched to the real selector count, 200 repeats, drawn from the same band) and prints
the contrast, so the arm cannot be read without its floor.

**What survives from L6, unaffected**: the calibration-versus-information decomposition, because
`per-target LOO` versus `global` is a held-out comparison and not a min-of-N statistic. The
calibration headroom (+0.114 A on the smoke) stands as measured, as does the observation that
fitting p weights on B = 25 candidates overfits.

**The general lesson, and it is the second time in this sprint that the coordinator has supplied
the failure mode the brief warns about:** an ORACLE arm needs a null of its own construction, not
merely a random-selection control. "Best of V things chosen with labels" is not a ceiling until it
is priced against "best of V things chosen without them".

---

## L9 — Consensus reproduces across a decade of instruments; typicality does not (2026-09-06, SELECT workstream)

n = 126, K = 500, identical candidate sets:

- **consensus medoid of the score's own top-75 vs its MATCHED RANDOM medoid**: 3.282 against 3.726,
  **-0.444 [-0.672, -0.231]**, fold [-0.550, -0.290], 80W/46L.
- **consensus vs distance argmin**: 3.282 against 3.454, **-0.172 [-0.322, -0.031]**, fold
  [-0.282, -0.076], 74W/47L -- reproducing Sprint 8's -0.172 [-0.316, -0.027] **to three decimals on
  a rebuilt instrument**.
- **Legacy argmin vs matched random**: +0.036 [-0.276, +0.353] -- null. **vs distance argmin**:
  +1.036 [+0.694, +1.384], 38W/88L -- decisively worse. Best |in-band Spearman| over its eleven
  components is `leg_torsion` at **+0.044**.

**A conflation in the programme's record, corrected.** "Consensus is the only signal with positive
in-band skill" was carried as one mechanism. It is two, and only one of them holds: the consensus
**medoid** works, while the typicality **score** does not -- in the ORACLE-defined band its in-band
Spearman is **-0.074 [-0.127, -0.021]**, significantly *negative*, and in the top-75 band it is
+0.015 [-0.060, +0.094], null.

---

## L10 — THE IN-BAND RESULT AT n = 126, and it is the sprint's sharpest mechanism (2026-09-06, coordinator)

`s17/inband.py`, n = 126, K = 500, band = the shipped objective's own top-25. All ORACLE arms
labelled. Random-in-band is the control every arm is priced against.

| arm | mean | median | vs random-in-band | W/L | |
|---|---|---|---|---|---|
| band best | 2.609 | 2.515 | -0.859 [-0.978, -0.746] | 126/0 | ORACLE |
| per-target fit, in-sample | 2.756 | 2.629 | -0.711 [-0.841, -0.593] | 121/5 | ORACLE |
| **per-target fit, LOO** | **3.075** | 3.142 | -0.393 [-0.523, -0.268] | 93/33 | ORACLE |
| best single feature | 2.725 | 2.610 | -0.743 [-0.855, -0.632] | 125/1 | ORACLE **min-of-N** |
| *its ZERO-INFO min-of-N null* | *2.757* | *2.705* | *-0.710 [-0.814, -0.609]* | *126/0* | *CONTROL* |
| **GLOBAL fit (deployable)** | **3.369** | 3.265 | **-0.099 [-0.194, -0.005]** | 70/56 | |
| **distance argmin** | **3.454** | 3.478 | **-0.014 [-0.118, +0.094]** | 62/64 | |
| consensus argmin | 3.369 | 3.309 | -0.098 [-0.204, +0.003] | 71/55 | |
| Legacy argmin | 3.415 | 3.235 | -0.053 [-0.160, +0.049] | 60/66 | |
| random in-band | 3.468 | 3.410 | +0.000 | — | CONTROL |

**THE HEADLINE.** *The shipped distance objective has essentially no in-band skill at all.* Its
argmin within **its own top-25** is **-0.014 [-0.118, +0.094], 62W/64L** against picking at random
from that same top-25. It is a garbage filter and nothing more. This is the mechanism behind
SELECT's in-band Spearman of 0.131 and behind their 58-functional null (-0.009 [-0.128, +0.110]):
**you cannot re-weight your way out of a signal that carries no information where the decision is
made.**

**A small in-band signal does exist, but not in the distance objective.** The GLOBAL linear fit over
all features reaches **-0.099 [-0.194, -0.005]**, an interval excluding zero, where the distance
objective alone reaches nothing. The extra content is in the other features -- consensus, radius of
gyration, contact fraction, Legacy terms -- not in the distogram.

**MIN-OF-N CHECK, per L8.** Real best-single-feature minus its zero-information null:
**-0.032 [-0.067, +0.001]** -- the real per-target minimum is marginally on the good side of its own
noise floor but the interval **touches zero**. SELECT's independent instrument with 25 signals gave
**+0.059 [+0.020, +0.100]** (real *loses*). **Both agree there is no established routing headroom**;
mine is null, theirs is significantly negative. L8's retraction stands, refined.

**CALIBRATION vs INFORMATION at n = 126**: calibration headroom (global -> per-target LOO)
**+0.294 A**; information deficit (per-target LOO -> band best) **+0.466 A**. Information binds, but
calibration is not negligible -- both are larger than the 12-target smoke suggested, and the
smoke's ratio (3x) was optimistic; at n = 126 it is 1.6x.

---

## L11 — Readout sweep paused to respect the compute ceiling (2026-09-06, coordinator)

Ten Python processes were running on a box of ~6.43 core-equivalents at a sustained 100% CPU --
past the brief's hard limit of 97%. The coordinator's `s17/readout.py` sweep was **stopped at
15/126** to reduce the coordinator's own footprint to the single shared instrument
(`s17/oracle_map.py`), rather than asking four agents to yield. It checkpoints every 15 targets,
so nothing is lost, and it will be relaunched when the map drains.

Recorded because the decision trades a coordinator experiment for sibling throughput, and the
alternative -- letting the box run oversubscribed at 100% -- would have violated a stated
constraint rather than merely slowed things down.

---

## L12 — WHY WIDENING K HURTS: the shortlist loses the answer (2026-09-06, SELECT workstream)

n = 126, identical candidate sets, target as the unit (`s17/results/sel_bench.json`).

| K | set ORACLE | ORACLE **inside the score's own top-75** | mean of that top-75 | argmin | consensus |
|---|---|---|---|---|---|
| 75 | 2.104 | 2.104 | 4.287 | 3.421 | 3.598 |
| 500 | 1.711 | 2.306 | 3.551 | 3.454 | 3.282 |
| 2000 | 1.504 | 2.406 | 3.526 | 3.520 | 3.344 |
| full | **1.313** | **2.572** | 3.557 | 3.562 | 3.461 |

    top-75 ORACLE, full vs K=75:  +0.468 [+0.304,+0.642]  fold [+0.255,+0.652]   46W/77L
    top-75 MEAN,   full vs K=75:  -0.730 [-0.906,-0.560]  fold [-0.844,-0.623]  104W/22L
    argmin,        full vs K=75:  +0.141 [+0.000,+0.291]  fold [+0.046,+0.253]   51W/68L
    consensus,     full vs K=75:  -0.138 [-0.381,+0.102]  fold [-0.324,+0.142]   70W/56L

**As K grows the shortlist gets BETTER ON AVERAGE and WORSE AT ITS BEST.** The extra windows look
more plausible without being nearer the native, and they **displace genuinely near-native members
out of the top-75**. The 0.468 A of shortlist-ceiling destruction exceeds the 0.141 A of realized
degradation, so it accounts for the whole effect.

**This is Sprint 16's G4 set-mean trap reproduced on the retrieval-width axis** rather than on a
filter's aggressiveness -- same score, same behaviour.

**And it explains the arm split.** Consensus *improves* with K (-0.138) while argmin *degrades*
(+0.141), because **consensus consumes the set MEAN, which is the quantity that improves, and argmin
consumes the set BEST, which is the quantity that is destroyed.** The current architecture pairs a
best-consuming readout with a mean-improving filter. **That mismatch, not ranking skill, is what
makes wide K harmful.**

**THE CONSEQUENCE, and it redirects the sprint.** The wide-K ceiling of 1.313 A is **not reachable
through a score-ranked top-B shortlist at any ranking skill whatsoever**, because the shortlist
ceiling there is 2.572 A. Reaching it requires a **diversity- or coverage-preserving shortlist** --
a **Problem C** question (ensemble construction), not the Problem B ranking question the sprint has
been treating as central.

---

## L13 — SELECT's lane, summarised; three more refutations (2026-09-06, SELECT workstream)

- **E1 identical-candidate instrument** -- built, n = 126 at K = 75/500/2000/full. The distance
  objective beats matched random by **-0.999 [-1.194, -0.799]** at every width; Legacy is null as a
  ranker *and* null as a gate against its own random gate; the **consensus medoid is the only arm
  beating both its matched control and the distance objective**.
- **E2 objective decomposition** -- 58 functionals of the distogram, leave-fold-out
  **-0.009 [-0.128, +0.110]**. **REFUTED.**
- **E3 set-conditional inference** -- answered inside E2; falsifier fired.
- **E4 target-level calibration** -- 27 native-free target features, five arms. **The calibrator
  LOSES to the constant baseline (+0.037 to +0.083) at every shrinkage, and the ladder is
  boundary-pinned at maximum shrinkage.** No feature clears the +-0.175 null band. **REFUTED at
  this feature class** -- which independently confirms L8's retraction from the other direction.
- **E5 hard targets** -- the 53 targets above 2.0 A at K = 500 fall to **24 over the full
  universe** (29 rescued, 55%); above 2.5 A, 27 -> 8; FAIL18 5/18 remain. The residual hard set is
  longer, less helical, less hydrophobic, more polar, with weaker BLOSUM matches and smaller
  libraries. Its **selection gap matches the easy set** (+2.437 vs +2.204), so it is a retrieval
  deficit *on top of* the same selection problem, and it is now worth only **+0.078 A** of the mean
  ceiling -- i.e. the hard-target frontier is much smaller than the sprint assumed.
- **E6 recovered fraction** -- **negative at every width** (-8.4% / -16.4% / -17.9%), which is L12
  measured directly.

---

## L14 — COVERAGE-PRESERVING SHORTLISTS: the ceiling is recoverable, the readout cannot use it (2026-09-06, coordinator)

`s17/coverage.py`, built as the direct consequence of L12. If a score-ranked top-B shortlist
destroys the shortlist ceiling, choose the shortlist for **coverage** of the candidate set's
structural variety instead. Diversity is measured in **pair-distance space**, not by pairwise RMSD,
making it O(K*B*npairs) rather than O(K^2) superpositions -- usable at the full universe and, not
incidentally, deployable.

**This is NOT a repeat of Sprint 16's diversity null, and the distinction is exact.** `divselect`
showed diversity-aware selection is a null for the AVERAGING readout, because that readout obeys
`readout^2 = member error^2 - diversity^2` and the two terms cancel at par under selection pressure.
**That identity does not govern an argmin readout.** Here the quantity of interest is the
shortlist's MINIMUM, not its mean. Same word, different functional -- which is exactly the check
BRIEF section 5 demands before reusing a refuted idea.

Six-target smoke, **NOT quotable**; the n = 126 run is queued behind a CPU hold:

| shortlist (K = 500, B = 75) | CEILING | mean | argmin | medoid |
|---|---|---|---|---|
| `score_topB` (incumbent) | 1.978 | 2.764 | 2.589 | 2.497 |
| `random_B` | 1.728 | 3.757 | 2.830 | 2.615 |
| `maxmin` | 1.949 | 4.476 | 2.589 | 3.897 |
| `strata` | 1.921 | 3.835 | 2.589 | 2.659 |
| **`kmeans_best`** | **1.703** | 3.781 | 2.589 | 2.675 |

    kmeans_best ceiling vs score_topB:  -0.275 [-0.473, -0.091]  5W/1L

**Coverage-preserving shortlists do retain near-native members the score discards** -- the L12
mechanism reversed, and the interval excludes zero even at n = 6.

**And the realized argmin is IDENTICAL across every arm (2.589)**, because the score's own best
candidate is present in all of them and the argmin is score-driven. This is precisely the case the
module's pre-registration named: *a ceiling gain is worthless unless a realized readout can use
it.* With the distance objective having **no in-band skill** (L10: -0.014 [-0.118, +0.094] against
random inside its own top-25), a wider shortlist hands the selector more of exactly what it cannot
discriminate.

> **Problems B and C must be solved together or not at all.** Coverage alone raises a ceiling no
> current selector can reach; ranking alone cannot recover what the shortlist no longer contains.

Note `maxmin` and `maxmin_gated` coincide at K = 500, B = 75 because the gate prefix
(8 x 75 = 600) exceeds K; they separate at larger K. Note also that pure `maxmin` has by far the
worst shortlist mean (4.476) and medoid (3.897) -- unconstrained diversity admits garbage, which is
why `kmeans_best` (coverage **and** score, one representative per cell) is the arm that works.

---

## L15 — THE ORACLE MAP AT n = 126 (2026-09-06, coordinator)

`s17/oracle_map.py`. **Independently cross-checked**: SELECT's separately written `sel_bench.py`
agrees with it to **3.2e-7 A** on 120 targets at all four shared widths.

### A. The ceiling as a function of candidate width [ORACLE]

| K | best | median | mean | coord-avg | diversity | <2.0 A present | <2.5 A present |
|---|---|---|---|---|---|---|---|
| 25 | 2.350 | 4.062 | 4.177 | 3.337 | 2.993 | 42.1% | 54.8% |
| 75 | 2.104 | 4.151 | 4.287 | 3.293 | 3.220 | 46.0% | 62.7% |
| 500 | 1.711 | 4.356 | 4.453 | 3.396 | 3.399 | 57.9% | 78.6% |
| 2000 | 1.504 | 4.564 | 4.596 | 3.519 | 3.487 | 69.0% | 88.1% |
| **full** | **1.313** | 4.880 | 4.816 | -- | -- | **81.0%** | **93.7%** |

### B. Realized vs ceiling, and vs random selection

| K | ORACLE | dist | **gap** | legacy | random | dist - random | W/L |
|---|---|---|---|---|---|---|---|
| 75 | 2.104 | 3.421 | **1.317** | 4.144 | 4.199 | -0.778 [-1.003, -0.555] | 88/38 |
| 500 | 1.711 | 3.454 | **1.743** | 4.490 | 4.381 | -0.927 [-1.181, -0.680] | 92/34 |
| full | 1.313 | 3.562 | **2.249** | -- | 4.755 | -1.193 [-1.475, -0.927] | 102/24 |

**The distance objective is a real filter at every width** (-0.78 to -1.19 against random, all
intervals excluding zero) **and the gap to the ceiling grows monotonically with K**, from 1.317 A at
top-75 to 2.249 A over the full universe.

### C. Where the angstroms are

    incumbent (deployed)                            3.204
    ORACLE best in top-75    [truncation ceiling]   2.104   the top-75 cut costs 0.791 A
    ORACLE best in K = 500   [truncation ceiling]   1.711   the K=500 cut costs 0.397 A
    ORACLE best in FULL pool [retrieval ceiling]    1.313

    targets whose FULL-pool best exceeds 2.0 A:  24/126   (at K = 500: 53)
    targets whose FULL-pool best exceeds 2.5 A:   8/126   (at K = 500: 27)

**A sub-2.0 A candidate is present for 81.0% of targets over the full universe** against 46.0% in
the shipped top-75. The answer is in the pool far more often than the programme believed -- and
**more than half of the "53 hard targets" were a truncation artefact**.

### D. Section 9's central question, answered

| K | extra ceiling | realized change | **recovered** |
|---|---|---|---|
| 500 | -0.393 | +0.033 [-0.086, +0.150] | **-8.4%** |
| 2000 | -0.600 | +0.099 [-0.025, +0.224] | **-16.4%** |
| full | -0.791 | +0.141 [-0.008, +0.293] | **-17.9%** |

> **Widening K buys ceiling and gives back accuracy.** The recovered fraction is *negative at every
> width*. Ranking is the blocker, recall is not, and more recall actively hurts -- for the mechanism
> in L12.

---

## L16 — SELECT closed; FEATURES opened on the last untested class (2026-09-06, coordinator)

SELECT's lane is complete (`s17/select_FINDINGS.md`, 621 lines). Its own summary of the closures:
selection is shut at three independent levels -- 58 distogram functionals, 29 native-free in-band
signals, 27 target-level calibration features -- and the readout triangle contradicts the sprint's
central architectural move, with the coordinate average dominating the medoid in **39/39 (K, m)
cells** and beating the distance argmin by 0.36-0.46 A at every K.

**Every signal tested so far is geometric or energetic** -- computed from a candidate's own
coordinates or from a physical potential. The one class never tested is **learned sequence
representations**, and there is a specific recorded reason to test it: the record carries an
overturned result that **ESM buys 0.288 A over one-hot on SELECTION (p = 0.005, n = 126)**,
measured globally and never in-band.

FEATURES spawned with the falsifier fixed in advance: *if no representation-derived feature beats a
zero-information constant alpha-helix in-band with an interval excluding zero, the last open feature
class is closed and the sprint's selection answer is final.* Roster: QUANTUM, PHYSICS, AUDIT,
FEATURES -- four, the standing default, cap of eight not approached.

**Why this is worth an agent rather than a footnote.** If the 0.288 A is in-band content it is the
most valuable quantity in the programme, because the in-band gap is 2.249 A and nothing else touches
it. If it is more garbage-filtering, the sprint's answer is complete rather than merely current.

---

## L17 — THE AVERAGING READOUT HAS NO HIDDEN HEADROOM (2026-09-06, coordinator)

`s17/avgceil.py`. Every ceiling the programme has quoted is a **best-member** ceiling -- 2.104 A in
the shipped top-75, 1.711 A at K = 500, 1.313 A over the full universe. Those bound a *selection*
readout. But selection is not the readout that wins: the coordinate average dominates the medoid in
39/39 cells and beats the argmin by 0.36-0.46 A. **The programme has been quoting the ceiling of
the wrong operator.**

The question matters because a coordinate average **can be better than any of its members** -- the
Krogh-Vedelsby identity subtracts the diversity term, so with partly cancelling errors the average
lands below the best member. If the ORACLE-best subset average sat far below 1.313 A, the sprint's
target would become *choose a better subset to average*, which is a **weaker requirement than
ranking**: a subset criterion need not order candidates, only partition them.

**It does not.** Five-target smoke, **NOT quotable**, n = 126 queued:

| m | ORACLE subset average | greedy-RANDOM control | score top-m |
|---|---|---|---|
| 2 | 1.735 | 2.824 | |
| **6** | **1.634** | 2.604 | |
| 10 | 1.674 | 2.694 | 2.587 |
| 16 | 1.736 | 2.657 | |

    best ORACLE subset average (m = 6)   1.634
    vs the best-member ceiling (1.615)   +0.019 [-0.532, +0.521]   3W/2L
    vs the greedy-RANDOM control         -0.970 [-1.588, -0.553]

**The pre-registered falsifier fires.** The ORACLE-greedy subset average never beats the best
member -- it starts at the best member and rises monotonically from m = 2. Adding *even the most
helpful available member, chosen with the native*, makes the average worse. The errors do not
cancel usefully, which reproduces the recorded law that *decorrelated errors exist and cannot be
exploited because the gain goes as the square of the weaker channel's skill.*

**Consequence.** The ceiling for **every** readout is the best-member ceiling. Averaging cannot
exceed it and selection cannot reach it, so the programme's ceiling is **1.313 A over the full
universe and the entire 2.249 A gap is ranking.** Problem C's "hidden headroom" question is closed
by measurement.

**Limitation, stated rather than buried.** Greedy forward selection seeded at the best member is a
*lower bound* on the true optimal subset, and it is a biased one -- a subset not containing the best
member could in principle average better. The result is that a strong ORACLE heuristic cannot beat
the best member, not that no subset can. The random-greedy control (-0.970) confirms the search
itself is working, so the null is not a failure of the procedure.

---

## L18 — CORRECTION TO L17: at n = 126 the averaging ceiling DOES beat the best-member ceiling (2026-09-06, coordinator)

`s17/avgceil.py`, n = 126, K = 500, greedy pool = the score's top 200, m up to 24.

| m | ORACLE subset average | greedy-RANDOM control | score top-m |
|---|---|---|---|
| 2 | 1.733 | 3.530 | |
| 5 | 1.600 | 3.235 | 3.231 |
| **6** | **1.598** | 3.194 | |
| 10 | 1.628 | 3.183 | 3.146 |
| 17 | 1.711 | 3.145 | |
| 20 | 1.743 | 3.147 | 3.091 |
| 24 | 1.785 | 3.151 | |

    best ORACLE subset average, m = 6      1.598 A
      vs the best-member ceiling (1.711)   -0.113 [-0.240, +0.023]   90W/36L
      vs the greedy-RANDOM control         -1.597 [-1.759, -1.441]

**L17's conclusion is wrong in SIGN and is corrected.** The 5-target smoke had the ORACLE subset
average rising monotonically from the best member (+0.019, 3W/2L); at n = 126 it **falls to 1.598 A
at m = 6**, beating the best-member ceiling by 0.113 A on **90 of 126 targets**. The random-greedy
control at -1.597 confirms the search is real rather than a min-of-N artefact.

**Labelled SUPPORTED, not ESTABLISHED**: the interval [-0.240, +0.023] **touches zero** despite the
90W/36L win rate. It is a real but modest headroom, not the large one that would redirect the
architecture.

**What it changes.** The true programme ceiling at K = 500 is **1.598 A, not 1.711 A** -- averaging
a well-chosen six-member subset beats the best single member. The optimum sits at **m ~ 6** and the
average is *worse* than the best member beyond m ~ 17, so the useful subset is small. The gap from
the deployable `score top-5` (3.231) to the ORACLE subset average at the same size (1.600) is
**1.63 A**, so subset choice has the same enormous gap as ranking does.

**THE PATTERN, and it is the third instance this sprint.** A small-n read reversed at full n --
after the n = 8 CI-excluding-zero false positive in Sprint 16, the min-of-N routing ceiling (L8),
and the units error (L3). **Every one of this sprint's coordinator errors has been a reading, not a
computation.** The rule that would have caught all three: *do not write a conclusion from a smoke
run at all, even a directional one -- record the smoke as a smoke and wait.*

---

## L19 — AUDIT: the sprint's opening premise is read backwards, and the coordinator's instrument carries three defects (2026-09-06, AUDIT workstream)

Every number re-derived from `s8/generate_univ/*.npz`, never read out of another workstream's
artefact. Two independent RMSD implementations at two dtypes agree to **4.8e-7 A** over ~2.1M
windows, and `oracle_map`'s ceilings all reproduce exactly (best@500 = 1.7108, best@full = 1.3134,
sel@500 = 3.4540, top-75 = 2.3062).

### A1 — the truncation "cost" is BELOW its own order-statistic null

`min` of N draws falls with N even when nothing new is present, so the 0.397 A truncation cost
needs a matched null. Under a **random permutation of the identical universe** the same statistic
is **0.473 A**:

    observed  best(500) - best(full)      0.397  [0.362, 0.429]
    MATCHED NULL, random order            0.473  [0.430, 0.511]
    observed - null                      -0.075  [-0.132, -0.022]   median -0.021, W/L 56/70

**The gain is significantly BELOW the null.** BLOSUM has already front-loaded the quality, so the
marginal return on the next 12,500 windows is *worse* than i.i.d. draws from the same library.
Stratified: on the 53 targets with K=500 best > 2.0 A, obs - null = **+0.046**; on FAIL18,
**-0.013** -- exactly the null and nothing more. **The whole retrieval apparatus buys +1.5 targets
of 2.0 A recall over a random 500 of the same library.** Uniform-effect null passed; fresh
interpreter delta = 0.00.

**Consequence for L1.** The ceiling numbers stand, but "the K = 500 ceiling is a truncation
artefact" must not be read as "the deep universe contains better structures the ranking was
hiding". It contains **more draws**, and the ranking extracts them slightly *less* efficiently than
chance would. L1 is corrected accordingly.

### A2 — the coordinator's budget table inverts the sprint's central comparison

`oracle_map` section C prints `incumbent 3.204` beside `raw coordinate average, K = 500 -> 3.396`.
**The programme's standing 3.050 is the average over the shipped TOP-75**, which AUDIT recomputes
at **3.0483**; the map's row is a different operator on a different set (13.2 of 75 members
shared). So the map makes averaging read **0.19 A worse** than the incumbent when it is
**0.156 A better** -- and **3.048 + the 0.155 A repair tax = 3.203, i.e. the incumbent IS that
average plus repair.** Any argument that averaging should be replaced by selection has to be made
against 3.048, not 3.396.

### A2b — section D is stronger than its own sentence, and is missing a baseline

"Recovered" is not "near zero"; it is **-8% to -18%**, and widening K makes the deployed answer
**worse** (+0.141 A at full vs K=75). But there is no composition-drift baseline: the exact random
pick degrades **+0.529 A** over the same widening, so **net of drift the selector gains -0.388 A
[-0.486, -0.268]**. Both halves must be stated. (Fold-clustered CI excludes zero; the plain target
bootstrap does not -- labelled accordingly.)

### A3 — the Sprint-16 frame trap is live in the shared instrument

Table A prints the **free-superposition** `mean` beside the **common-frame** `diversity`. The
identity `readout^2 = <err^2> - div^2` is **exact** with the common-frame error (residual
7e-14 A^2) and off by **-2.10 A^2** -- about 18% of readout^2 -- with the printed column; using its
arithmetic rather than quadratic mean, up to **8.9 A^2**. `s17/readout.py` already does this
correctly; `oracle_map` does not.

### Further defects

- **A4** the map's random control is 3 draws (per-target SE 0.763 A) when the exact expectation is
  free.
- **A6** `corr(U, chain length) = -0.988` -- "13,000-27,000 windows per target" is **1/n**, and the
  1CB3 exhibit in L1 is simply the target with the most draws.
- **A10** `oracle_map.json` carries **no completion flag** and was read by three workstreams at
  60/126 rows, with fold 3 under-represented 2.7x.
- **A11** the instrument's minimum detectable effect at 80% power is **0.084 A**.
- **A14** two different objects share the label "top-75".

### Clean survivals, stated with equal weight

`u["order"]` is provably native-free. The 1.711 A oracle is robust to the 11.4% of the pool set by
corpus-order tie-break (1.708 A over random tie-breaks, sd 0.018). **No `s17/` module touches the
sealed benchmark; no bare `hash()` anywhere.** `selection_theory.rho_from_accuracy` consumes the
right functional, and s14's 0.600/0.638 really are pairwise accuracies -- with one caveat, that
s14's null is 0.505 rather than 0.500. `readout.py` and `sel_bench.py`'s matched controls hold.

### AUDIT's own falsified hypotheses, preserved

Two of its priors were wrong and are reported as such: **BLOSUM ordering does have skill** (small,
decaying), and **the deep universe is NOT near-duplicates** -- d* = 1.355 A, a dense continuum whose
distinct-conformation count grows 12.2x while K grows 36.9x. That answers the redundancy question
the coordinator raised at L1: the wide-K ceiling is real structure, not duplicates. AUDIT also made
a section-5-shaped error itself (its first recall control compared two things equal by
construction), preserved in place in `audit_strat.py`.

### The highest-value unrun experiment, per AUDIT

Widening K takes the pool from **27.3% to 20.8% peptide**, and the programme's own record says the
**peptide corpus carries ~7x the sequence-structure channel** of the protein fragments. **Every
K-effect measured this sprint is confounded with that composition change and nobody has controlled
it.**

---

## L20 — PHYSICS: the Legacy null reproduced independently, and one result that damages the shortlist rather than the energy (2026-09-06, PHYSICS workstream)

`s17/phys_ident.py`, K = 500 BLOSUM pool, all 126 targets, target as the unit, fold-clustered CIs,
`argmin_tied` throughout, in-band defined as `d <= pool_best + 1.5 A`.

### Independent reproduction of the two nulls

- **Legacy in-band vs matched-random in-band**: **+0.0572 [-0.0348, +0.1547]**, 55W/69L -- null,
  sign of harm.
- **Distance objective in-band vs matched random**: **+0.0587 [-0.0395, +0.1606]**, 57W/62L --
  also null, and slightly *worse* than random. **The coordinator's floor claim (L10) reproduces on
  an independently built instrument.**
- Every Legacy gate is null against its matched random gate on selected RMSD (worst |mean| 0.046).
- `leg_torsion` is the strongest component in-band (rho 0.101 against the distance score's 0.126)
  and the **only** score in a 17-row table with a negative in-band selection delta,
  -0.0446 [-0.1378, +0.0482], 66W/57L -- interval includes zero.

### Legacy DOES preserve coverage, and destroys the best candidate while doing it

On the 73/126 targets that have a sub-2 A candidate, recall and set-best on **one** denominator:

| Legacy gate, f = 0.50 | gate | random | delta |
|---|---|---|---|
| near-native **recall** | 0.613 | 0.501 | **+0.112 [+0.023, +0.198]**, 49W/21L |
| set **BEST** | 1.427 | 1.208 | **+0.219 [+0.031, +0.444]** -- significantly worse |

**Mechanism, measured directly rather than inferred**: in Legacy's own normalised order the pool's
single best candidate sits at **0.450** while the **median** sub-2 A candidate sits at **0.410**
(random = 0.500). **Legacy likes near-natives, and likes the BEST one systematically less than the
typical one.** That is Sprint 16's G4 set-mean trap resolved to the level of an individual
candidate, which is as sharp as that mechanism has ever been stated.

**`leg_torsion` is the exception and is the one physics arm worth keeping**: recall
**+0.159 [+0.073, +0.239]** at f = 0.50 with set-best damage **null** (+0.011 [-0.079, +0.128]).
It is the only gate measured anywhere in this programme that buys coverage without destroying the
best member.

### THE RESULT THAT DAMAGES THE SHORTLIST, NOT THE ENERGY

Shipped distance top-M against a **matched random shortlist of the same size**, n = 126:

| M | P(has sub-2 A), score | random | delta | shortlist ORACLE best, score - random |
|---|---|---|---|---|
| 25 | 0.373 | 0.357 | +0.016 [-0.045, +0.077] null | **+0.195 [+0.036, +0.359]** |
| 75 | 0.437 | 0.452 | -0.016 [-0.066, +0.034] null | **+0.242 [+0.110, +0.382]** |
| 150 | 0.476 | 0.508 | -0.032 [-0.090, +0.024] null | **+0.179 [+0.063, +0.305]** |
| 300 | -- | -- | -- | **+0.097 [+0.016, +0.193]** |

> **The score-ranked shortlist is no better than a random shortlist at KEEPING the answer, and is
> measurably WORSE at keeping the BEST one -- at every width.**

This is L12's K-widening finding reproduced **at fixed K**, and it sharpens what the distance
objective is worth: it raises the shortlist **mean** (SELECT: the average over the score's top-300
is 3.171 against 3.402 over a matched random 300) while actively **degrading the shortlist's
minimum**. Win/loss is near-even (59/67, 63/63, 72/53, 75/26) so the median-vs-mean warning applies
and the effect is mean-carried -- but the direction is identical at all four widths.

### Sprint section 34 -- disagreement as information -- FALSIFIED cleanly

Every distance/Legacy disagreement feature has leave-fold-out skill against the selector's realized
error (rho_dl 0.449, overlap 0.397, rank-of-pick 0.428), and **all of it vanishes once
`dist_score` is partialled out** -- the trivial native-free feature "how good the objective thinks
its own pick is", which alone has LFO 0.598. Partial LFO ranges -0.148 to +0.168 with no consistent
sign. **Disagreement is subsumed by the objective's own confidence.** Route closed.

### Ca-PRESERVING REPAIR -- the pre-registered hypothesis is HALF true (n = 55 of 126, interim)

**Ca-FIXED AMBER** (particle mass 0, epsilon exactly 0) at **identically zero Ca cost**:

- heavy-atom clashes below 2.0 A: **3.40 -> 0.00**; below 2.6 A: 11.02 -> 0.03
- minimum heavy separation: **2.210 -> 2.857 A**
- **but**: fails the convergence gate on **17/55** targets, **degrades** Ramachandran-favoured
  **0.849 -> 0.722** (below the do-nothing input), and introduces **cis peptide bonds at 0.313 mean
  fraction** with 60.6 deg mean omega deviation.

> **Sterics are free. Torsion and peptide-bond plausibility are not.**

### A band-definition discrepancy, recorded so the two are never quoted side by side

SELECT reports a constant ideal alpha-helix at **+0.053** in the *shipped top-25* band,
indistinguishable from the distance objective there. PHYSICS reports the same reference at
**+0.2913 [+0.1966, +0.3899]**, 37W/83L, with in-band Spearman **-0.135**, on the
*pool_best + 1.5 A* band -- i.e. much **worse** than random. **These are two different bands, not a
contradiction**, but they must never be quoted together without saying which.

**And a units bug worth propagating**: `core.geometry.build_backbone` takes **RADIANS**. Passing
degrees silently builds a different constant conformation that scores Ramachandran 0.000. Any
workstream constructing a constant secondary-structure reference must check this.

---

## L21 — Coverage sweep stopped at 20/126; the decision it was built to make is already made (2026-09-06, coordinator)

`s17/coverage.py` was crawling on a saturated box (20/126 after 1809 s, four agents running) and is
stopped. Recorded because it trades a coordinator experiment away, and the reasoning should be
visible rather than silent.

**The strategic question it was built to answer has been answered at n = 126 by PHYSICS on a
cleaner control.** The premise -- that a score-ranked shortlist discards near-native members a
coverage-aware shortlist would keep -- is confirmed from the other direction: a **matched random
shortlist** has a significantly *better* ORACLE best than the score's own shortlist at **every**
width (+0.195 / +0.242 / +0.179 / +0.097 A at M = 25/75/150/300, L20). A random shortlist is the
crudest possible coverage-preserving construction, and it already beats the score.

What the full sweep would have added is **which** coverage construction is best -- `kmeans_best`
versus farthest-point versus strata. That is a refinement of a decision already taken, not the
decision itself.

**What is kept from it**, labelled as what it is: a **6-target smoke, not quotable**, in which
`kmeans_best` raised the shortlist ceiling by -0.275 [-0.473, -0.091] against `score_topB` while
**every arm's realized argmin was identical** -- the case the module's own pre-registration named,
where a ceiling gain exists and no current readout can use it. With the distance objective having
no in-band skill (L10), a wider net hands the selector more of exactly what it cannot discriminate.

> **Problems B and C must be solved together or not at all.** Coverage alone raises a ceiling no
> current selector can reach; ranking alone cannot recover what the shortlist no longer contains.

---

## L22 — QUANTUM: every matched-budget comparison in three sprints was run past classical saturation (2026-09-06, QUANTUM workstream)

### The finding that invalidates the framing, including the coordinator's

**Greedy 1-opt reaches the CERTIFIED GLOBAL OPTIMUM of the deployed objective in 100% of 152
(target x seed) cells at 1,024 objective evaluations, and 62.5% at 36** -- one coordinate-descent
pass. **The CVaR-VQE budget used throughout Sprints 15-17 is 8,192.**

> **Every "matched-budget" quantum-versus-classical comparison this programme has run was conducted
> at least 8x past the point where the classical control saturates. A budget the control does not
> need is not a budget.**

This invalidates the premise of the Pareto framing the coordinator briefed, and **no audit caught it
because every individual number was correct** -- the failure is in what "matched" meant, not in any
computation. It is the sprint's sixth instance of a correctly measured quantity read as a different
quantity, and the largest in scope.

### The mechanism, exact, from the Walsh/Pauli-Z spectrum (19 targets, full registers)

The deployed objective `hamil` puts **61.3% of its variance at Pauli weight 1**, **93.0% at weight
<= 2**, and **95.7% of its weight-2 mass is intra-residue**. It is very nearly a **separable
per-residue field**. Its weight-<=1 truncation -- closed form, **zero search** -- has a certified
argmin of **2.411 A against the full objective's 2.661 A**, i.e. *better than the objective it
truncates*. The ORACLE truth has mean Pauli weight **3.72**, with only 31.4% at weight <= 2.

> **There is no search problem here for a quantum device. The objective is nearly separable and the
> answer is not.** The existing higher-degree content actively points the wrong way.

### Every pre-registered falsifier fired

- **P1, the Pareto test.** All 10 VQE configurations are epsilon-dominated by the 21-point classical
  set on the readout the pipeline consumes, every interval excluding zero (-0.129 to -0.190 A).
  **But the classical leave-one-out null looks identical** (-0.139 to -0.174), so the honest
  statement is *the VQE is one more classical sampler, slightly worse* -- **not uniquely bad**.
- **P3, local VQE.** At 1/4 the exhaustive cost greedy and Metropolis certify the local optimum in
  100% of windows; VQE reaches 68.5-72.2%, and its proposals are worse at set-best than a random
  draw.
- **P4, reweighting.** VQE weights are indistinguishable from **a permutation of their own weights**
  (3.299 vs 3.287 A) and worse than doing nothing (3.135). Entropy-matched Boltzmann reaches 2.865.
- **P5, mixtures.** Best quantum arm contributes 2.07 unique elite structures per target; best
  classical, 5.64.
- **THE ENTANGLEMENT CONTROL, which Sprint 16 never ran.** Deleting the CNOTs gives `mps2fn`, a
  classical product-Bernoulli variational model at identical budget, estimator, optimiser and seed.
  It is **null on member error and on the readout at every alpha**, and the entangled circuit is
  **significantly less diverse**. Adding those three classical points collapses every remaining
  non-dominated cell.

### The one live positive, and it is bounded

`log q_theta` beats both the objective it was trained on (+0.120 in-band) **and** the exact
mean-field classical model (+0.083 [+0.001, +0.179] at alpha = 0.02). Labelled **OPEN, not
SUPPORTED**: the interval clears zero by 0.001, the median (+0.054) sits far below the mean, W/L is
5/3, fold sign 2/4, n = 8 -- the ledger's own concentration warning -- and **it does not convert
through any readout**.

**alpha-as-temperature, refined**: `T_eff` falls 0.668 -> 0.190 monotonically in physical units, but
`KL(q_theta || Boltzmann)` is **4.55-5.99 bits**. **Alpha sets a temperature; q_theta is not that
Boltzmann law.**

### The sprint's one positive lever, and it is not selection

> **Local refinement of a retrieval candidate is worth ~0.7 A** -- 3.636 -> 2.95 A over 6-residue
> windows, **certified by exhaustive enumeration** at 4,096 evaluations per window, with greedy
> reaching the same optimum at 1,024. **Cheap, exact, native-free.**

With selection closed at four levels, this is the only measured mechanism in the sprint that moves
accuracy, and it attacks **Problem A/C** -- improve the candidates -- rather than Problem B.

### Errors QUANTUM caught in its own work, preserved

A **sign error in the primary Pareto statistic** that manufactured large quantum wins (-0.48 to
-0.89 A with intervals excluding zero), and a second defect where the conditional margin was
undefined exactly on the cells the VQE could have won -- replaced by an always-defined
epsilon-dominance indicator. Both recorded in the findings rather than quietly fixed.

### Limitations, stated by the workstream

Pareto n = 10 targets (not 19) and representation n = 8; **both runs were stopped by the workstream
itself when free RAM hit 0.30/0.68 GB**, which the brief forbids exceeding. Checkpoints are on disk.
The box was held at 100% CPU by siblings throughout, so all wall-clock numbers are contended; the
blocking CPU gate never opened and was replaced by a non-blocking contention sampler that records
CPU and RAM into the artefacts.

---

## L23 — FEATURES: the last open feature class is closed, and the binding constraint is the SHORTLIST (2026-09-06, FEATURES workstream)

### The lane-level falsifier fired

**No representation-derived feature has in-band skill.** Nothing in six classes of learned sequence
representation beats a constant ideal alpha-helix **and** random-in-band on selected RMSD with an
interval excluding zero, in either deployable band.

### The recorded 0.288 A is not in-band content

At **matched architecture**, varying only the per-residue block:

- ESM vs one-hot, top-25 band: **+0.036 [-0.063, +0.139]** -- ESM *worse*
- ESM vs one-hot, top-75 band: **+0.116 [-0.004, +0.243]**, fold [+0.019, +0.174] -- ESM *worse*
- target-conditioned pseudo-likelihood (ESM-2-650M, masked, 126 targets): does not separate from a
  **composition-only** control
- residue-aligned window-embedding cosine: counted ordering accuracy **0.503** against a 0.500 null

> **ESM buys a better FILTER, not a DISCRIMINATOR** -- the same object L10 showed the distance
> objective to be. The programme's recorded "ESM buys 0.288 A on selection" is global skill, and
> global skill is not what decides an argmin.

### One real ESM-specific signal, on the correlation axis only

ESM-2's **contact head** reaches in-band Spearman **+0.116 [+0.047, +0.184]** at top-75 (accuracy
0.541), and against its own **ESM-free twin** the increment is **+0.050 [+0.005, +0.097]** -- the
first ESM-specific in-band content measured in this programme. But the twin control shows most of
it is **compactness** (choosing pairs by ESM rather than by the distogram: +0.000 [-0.078, +0.080]),
and **none of it reaches the argmin in any band**.

### THE EXACT RESULT, and it is the sprint's central conclusion arriving from a third direction

Fed through `s17/selection_theory.py`'s band curve:

- rho_S = +0.116 delivers **3.393 A at B = 25** against a band mean of 3.502 -- **0.109 A** against
  a 2.249 A gap.
- **At B <= 25 the 2.5 A target is UNREACHABLE AT rho = 1**: a *perfect* ranker inside the shipped
  top-25 returns **2.609 A**.

> **The shortlist, not the ranker inside it, is the binding constraint.** Labelled **EXACT** -- it
> follows from the measured band contents, not from a model.

This is L12 (widening removes the answer from the shortlist) and L14 (coverage raises a ceiling no
readout can use) arriving independently from the feature side. **Three directions now agree.**

### F7 -- the ESM channel has no shortlist-construction value either

Its shortlist is **significantly worse than matched random** (+0.315 to +0.449). And it
**independently reproduced the PHYSICS result** that the shipped score's top-B shortlist has an
ORACLE best worse than matched random: **+0.209 / +0.226 / +0.170** at B = 25/75/150 against
PHYSICS's +0.195 / +0.242 / +0.179 -- two workstreams, two instruments, same numbers.

A `mix` arm beat the incumbent's ceiling with **both CIs excluding zero** -- and its
**pre-registered matched-random control killed it**: random does as well or better at every size.
Realized argmin identical to three decimals across arms.

### Two things flagged and NOT claimed

- `cf_topd_uniform` -- an **ESM-free** control, the mean realised distance over the *n* pairs the
  distogram predicts closest -- is the battery's best arm (top-75 **-0.135 [-0.252, -0.022]** vs
  random, beating both constants and the shipped score). But it is **1 of ~75 uncorrected tests**
  and **flips sign on the ORACLE band**. **PLAUSIBLE**; belongs to SELECT's E2 lane.
- The ESM contact signal is **length-gated** (rho +0.188 with -0.271 selected at length 14-16;
  rho -0.047 with +0.203 at 9-11). **PLAUSIBLE, post-hoc.**

### A third refinement to the min-of-N rule (L8)

**The min-of-N contrast changes sign with the band** on one instrument (+0.020 / -0.061 / +0.083).
So the rule adopted at L8 needs a further clause: **a min-of-N contrast must name its band.**

### Compliance

The coordinator's radians warning (L20) was **checked empirically before use** -- alpha-helix
d(i, i+3) = 5.20 A, rise 1.58 A/residue, against the degrees bug's 7.49 -- and pinned as an
assertion in `feat_lib.check_constants`. No benchmark reference in any module; all randomness via
`stable_rng`.

---

## L24 — PHYSICS: the validity frontier splits in two, and the incumbent's restraint set wins (2026-09-06, PHYSICS workstream)

Pre-registered, native-free, fold x length stratified 30-target subsample reused **verbatim** from
Sprint 16 so no new selection freedom is created. **SHAPE only** -- no rung is re-quoted at full n
as an out-of-sample choice. Do-nothing on this subsample is 3.5061 A.

| rung | gate excl | Ca disp | RMSD cost vs do-nothing | ramaFav | clash<2.0 | geomDev | **cis** | omega dev |
|---|---|---|---|---|---|---|---|---|
| none | 0 | 0.000 | -- | 0.781 | 7.37 | 0.2885 | **0.000** | 0.0 |
| **`cafix`** (eps = 0) | 15 | **0.000** | **-0.0000** (tautology) | 0.691 | **0.00** | 0.0476 | 0.313 | 59.2 |
| `ca1000` | 6 | 0.185 | **+0.0081 [-0.0033, +0.0187]** | 0.647 | **0.00** | 0.0706 | **0.539** | 89.2 |
| `ca300` | 2 | 0.317 | +0.0305 [+0.0085, +0.0519] | 0.649 | **0.00** | 0.0709 | **0.574** | 96.4 |
| `ca100` | 2 | 0.499 | +0.0662 [+0.0331, **+0.0996**] | 0.717 | **0.00** | 0.0546 | 0.488 | 89.1 |
| `ca30` | 2 | 0.782 | +0.1402 [+0.0843, +0.1968] | 0.814 | **0.00** | 0.0275 | 0.201 | 47.9 |
| `ca3` | 2 | 1.111 | +0.2682 [+0.1588, +0.3833] | **0.847** | **0.00** | 0.0162 | **0.084** | 22.5 |
| `free` | 2 | 1.878 | +0.4586 [+0.2374, +0.6876] | 0.845 | **0.00** | 0.0172 | 0.087 | 21.6 |
| **`k30`** (N/CA/C, incumbent) | 2 | 0.824 | +0.1217 [+0.0688, +0.1769] | 0.828 | **0.00** | 0.0221 | **0.086** | 25.0 |

### The five readings

1. **H1a HOLDS, with no margin.** Every rung whose realised Ca displacement is <= 0.5 A costs less
   than 0.10 A of accuracy -- but `ca100`'s interval upper bound is **+0.0996 against a 0.10
   threshold**. A pass by 0.0004 A.
2. **THE FRONTIER SPLITS IN TWO, and this is the finding.** **Steric validity is FREE**: *every*
   rung, `cafix` included, reaches **0.00** clashes below 2.0 A at a Ca cost anywhere from 0 to
   1.88 A -- the curve is flat. **Torsional and peptide-bond validity is a STRICT TRADE**:
   Ramachandran rises 0.647 -> 0.847 and cis falls 0.574 -> 0.084 **monotonically in the realised
   displacement**, and no rung buys them cheaply.
3. **H1c is REFUTED.** No epsilon reaches the unrestrained arm's validity at < 0.05 A of Ca cost.
   The rungs under 0.05 A return Ramachandran 0.647/0.649 against `free`'s 0.845 and cis
   0.539/0.574 against 0.087. **There is no free lunch on the torsional axis.**
4. **A NULL AGAINST THE LANE'S OWN PREMISE.** The framing was that freeing N and C -- restraining
   *only* Ca -- would let the force field repair in coordinates Ca-RMSD is blind to. At matched
   displacement it does the **opposite**: `k30` (N/CA/C, disp 0.824) beats `ca30` (Ca only, disp
   0.782) on accuracy (+0.122 vs +0.140), Ramachandran (0.828 vs 0.814), geometry (0.0221 vs
   0.0275) and cis (0.086 vs 0.201). **The incumbent's restraint set is better than the Ca-only set
   on essentially every axis.** Sprint 16's 55.1% non-torsional share of AMBER's displacement (J9)
   is **not** an opportunity to be harvested by loosening N and C.
5. **The tight-Ca regime is the WORST place on the ladder for peptide geometry, and it is not
   monotone.** cis peaks at **0.574 at `ca300`** -- worse than `cafix`'s 0.313 and far worse than
   `free`'s 0.087 -- with omega deviation peaking at **96.4 deg**. Mechanism: a Ca-only restraint
   at high k pins the contracted spacing while leaving N and C free to absorb it by flipping omega.
   **A restraint that holds Ca and frees the peptide plane is the one thing not to do.**

### And two pieces of method discipline worth more than the numbers

**The falsifier the coordinator handed this lane could not fire, and the workstream said so before
the run rather than claiming a pass afterwards.** At eps = 0 the Ca coordinates do not move, so
Ca-RMSD equals the do-nothing RMSD **exactly, by construction** -- the brief's falsifier ("accuracy
cost above 0.10 A against doing nothing") is untestable at that rung. It is recorded in
`PREREG_phys.md` in advance, and the eps = 0 rung is judged on validity alone. **That is the
Sprint-16 failure mode caught in the brief itself rather than in a result.**

**And the workstream found the radians units defect in its OWN zero-information control**:
`core.geometry.build_backbone` takes radians, its first constant-alpha-helix reference passed
degrees, silently building a different conformation (recomputed torsions -25.9 / -172.9 deg) scoring
Ramachandran **0.000** instead of 1.000. Every affected number was recomputed and nothing built on
the broken control survives into the document.

### The consequence for the sprint's Level-1 deliverable

**The repair tax cannot be eliminated.** Sterics are free at any restraint strength; torsional and
peptide-bond plausibility are bought strictly with Ca displacement; and the restraint set the
pipeline already ships is the best of the ten tested. **Level 1 -- recover ~3.05 A with valid
all-atom structures -- is not achieved, and the incumbent's repair stage is validated rather than
improved.**

---

## L25 — The one target-level signal that works, and it is the trivial one (2026-09-06, coordinator, from PHYSICS §6)

Buried in the refutation of sprint section 34 is the strongest target-level predictor the sprint
found, and it is not a disagreement feature at all:

> **`dist_score` -- how good the distance objective thinks its own pick is -- predicts the
> selector's realized error at leave-fold-out Spearman 0.598.**

Every disagreement feature (rho_dl 0.449, rank-of-pick 0.428, overlap 0.397) is **entirely
subsumed** by it: partialled, the six features span -0.148 to +0.168 with no consistent sign across
two responses.

**Why this matters despite being trivial.** It does not help *select* -- it is a per-target scalar,
constant across candidates, so it cannot order them. But target-level **difficulty estimation** is
the premise of sprint sections 36 and 39 (adaptive generation: easy targets get selection effort,
hard targets get generation effort), and this is a native-free signal at LFO 0.598 that needs **no
second energy model and no new machinery**. The calibration lane's 27 engineered features could not
beat a constant baseline; the objective's own confidence does this on its own.

**Recorded as the practical residue of a refutation**: sprint section 34's route to calibration is
closed through the physics channel, and what survives is *a selector that scores its own pick badly
is on a hard target.*

---

## L26 — SPRINT CLOSE (2026-09-06, coordinator)

Six workstreams reported: SELECT, QUANTUM, PHYSICS, AUDIT, FEATURES, and the coordinator's own
instruments (`oracle_map`, `selection_theory`, `readout`, `inband`, `avgceil`, `coverage`). Four
agents ran continuously, peaking at five; the cap of eight was never approached.

**Nothing beat the incumbent, and two of its stages were validated rather than left alone.** The
sealed benchmark is unspent and audit-verified untouched. Twenty-six ledger entries, six preserved
coordinator errors, five independent closures of the selection problem, and one exact result --
*at B <= 25 the 2.5 A target is unreachable at rho = 1* -- that redirects the programme from
building a better ranker to building a better shortlist.

**The sprint's own biggest weakness, stated plainly**: every K-effect is confounded with a
composition change nobody controlled (27.3% -> 20.8% peptide as K widens, against a recorded ~7x
sequence-structure channel in the peptide corpus). It is the cheapest outstanding experiment and it
bears on three headline findings.

---

## L27 — CORRECTION: "steric validity is free" was quoted against the WRONG BASELINE, and the workstream caught it in its own draft (2026-09-06, PHYSICS)

**The convergence gate selects a different subset at every rung.** Quoting a gated arm's validity
against the *ungated* instrument-wide input is the exact failure mode BRIEF section 5 warns about --
a quantity measured correctly and compared to a different quantity.

`cafix` converges on the **78 least broken targets**, whose input carries **0.31** sub-2.0 A clashes,
**not** the instrument-wide 5.61. Corrected, against **its own gated input**:

| | uncorrected (what the coordinator published) | **corrected, gated input** |
|---|---|---|
| sub-2.0 A clashes | 7.37 -> 0.00 | **0.31 -> 0.00** (13 better / 0 worse) |
| sub-2.6 A clashes | 11.0 -> 0.03 | **2.03 -> 0.01** (31 better / 0 worse) |

**The effect is real and one-sided, and it is an order of magnitude smaller than the number nearly
published.** Every rung now carries an `in@rung` row (`s17/phys_ca_report.py`).

**The coordinator's report and artifact carried the uncorrected comparison and are amended.** The
qualitative statement -- steric repair is available at zero Ca cost -- survives; the magnitude does
not.

Two further full-n numbers that sharpen the same rung: at eps = 0 the Ca move by **7.1e-15 A**, so
`cafix` is **the first AMBER arm in this programme ever to PASS the rotated-frame null** (max
0.00000) -- *because it buys nothing*. Its Ramachandran cost at full n is **0.906 -> 0.734**
(1 better / 44 worse), cis **0.000 -> 0.313**, omega deviation 60.2 deg, and the convergence gate
fails on **48/126** against the incumbent's 3.

---

## L28 — Sprint 16's cis-peptide defect is EXPLAINED, and it exonerates the force field (2026-09-06, PHYSICS)

Sprint 16 recorded (claim G12) that AMBER "introduces cis peptide bonds on 44 of 126 targets" and
filed it as a new defect of the repair operator. **It is not a defect of AMBER.**

**The coordinate average contracts the backbone by 22.4%** -- Ca-Ca spacing **2.949 A against an
ideal 3.80 A** -- and the minimiser is handed a chain whose peptide planes cannot be satisfied at
that spacing:

    Spearman(input Ca-Ca spacing, cis fraction)                 -0.949
    the same, conditioned on convergence                        +0.807

> **Fix the averaging operator, not AMBER.** The force field is responding correctly to an
> unphysical input.

This also explains L24's non-monotone finding directly: a Ca-only restraint at high k **pins the
contracted spacing** while leaving N and C free to absorb it by flipping omega, which is why cis
peaks in the tight-Ca regime (0.574 at `ca300`) rather than at either end.

---

## L29 — THE DECISIVE AMBER-vs-LEGACY COMPARISON, at full scale (2026-09-06, PHYSICS)

Sprint section 32's mandated experiment, run properly: **63,000 genuine ff14SB/GBn2 single points,
all 126 targets, identical K = 500 candidate sets, no shortlist.**

**AMBER lost every ranking role it was given.**

| | AMBER | Legacy |
|---|---|---|
| global Spearman | **-0.027** | +0.307 |
| in-band AUROC | **0.499** | -- |
| argmin | **5.191 A** (random: 4.427) | 4.490 A |
| incremental value over the distance features | **-0.004 [-0.037, +0.030]** | -- |

`angle` and `torsion` are **significantly anti-informative**. An all-atom force field, scored
exactly, ranks peptide candidates **worse than chance** and worse than random selection.

**What Legacy knows that AMBER does not**: garbage rejection -- and one specific term.
**`leg_contact` (the Miyazawa-Jernigan contact term), which Sprint 16 called useless, reaches
+0.080 [+0.032, +0.128] against a +0.010 null.** It carries real information and **still selects
0.136 A worse than random**, which is the cleanest statement yet of the gap between *having
information* and *being able to use it as a ranker*.

**What AMBER knows that Legacy does not**: where the atoms go. That is a statement about
stereochemistry, not about ranking, and it is the only role either model has earned.

**Instrument reproduction**: an independently built OpenMM System reproduces
`s16/results/repair_A.json` at **max |dRMSD| = 0.000e+00 A** and max |dE| = 2.7e-13 kcal/mol,
n = 126, excluding the same three gate targets.

**One operational hazard recorded**: two `phys_ident` processes racing on a single output file let a
**9-row partial overwrite a complete 126-row artefact**. The `complete` flag adopted for
`oracle_map.json` at L19/A10 should be adopted programme-wide.

---

## L30 — THE OBJECTIVE IS NOT WEAK, IT IS WRONG (2026-09-06, coordinator)

`s17/refine.py`, n = 126. Start from the best structure the pipeline actually builds -- the
coordinate average of the shipped top-75 -- and move toward the distance objective.

| arm | RMSD | median | vs the start | W/L | objective | disp (rad) |
|---|---|---|---|---|---|---|
| **avg (start)** | **3.048** | 2.837 | -- | -- | 192.8 | -- |
| ideal-geometry projection | 3.213 | 2.966 | +0.164 [+0.131, +0.199] | 15/111 | 192.8 | -- |
| **refine_full** | **3.610** | 3.370 | **+0.561 [+0.407, +0.712]** | 31/95 | **56.5** | 5.035 |
| refine_win4 | 3.645 | 3.463 | +0.597 [+0.448, +0.745] | 26/100 | 57.8 | 5.318 |
| refine_win6 | 3.616 | 3.431 | +0.567 [+0.412, +0.730] | 28/98 | 57.2 | 5.365 |
| refine_win8 | 3.611 | 3.438 | +0.563 [+0.415, +0.714] | 29/97 | 56.3 | 5.227 |
| rand_obj (shuffled distogram) | 5.395 | 5.192 | +2.346 [+2.001, +2.688] | 13/113 | -- | CONTROL |
| rand_move (matched magnitude) | 4.547 | 4.490 | +1.499 [+1.217, +1.786] | 18/108 | -- | CONTROL |

> **The objective falls by 71% (192.8 -> 56.5) and the RMSD rises by 0.561 A.** The
> pre-registered falsifier fired on every arm.

**And the move shape does not matter.** Windowed refinement at w = 4, 6 and 8 lands within 0.035 A
of the full fit. The enumerated instrument's 0.7 A local-refinement gain (L22) **does not transfer
to the real instrument**, and the difference is the *objective*, not the move.

**Both controls are far worse** (+2.346 and +1.499), so this is not random damage. The refinement is
purposeful, competent, and aimed at the wrong place.

### The unified diagnosis

This is not "the objective is a weak ranker". It is:

> **The distance objective's optimum is 0.56 A WORSE than the structure the pipeline already
> builds. The pipeline works BECAUSE it does not optimise its own scoring function.**

Everything else in the sprint follows from it and is now one fact rather than five:

- the objective cannot rank in-band (-0.014 [-0.118, +0.094]) -- *its gradient points away from the
  native, so its ordering near the native is meaningless*;
- 58 functionals of it gained nothing (-0.009 [-0.128, +0.110]) -- *re-weighting a signal that
  points the wrong way cannot make it point the right way*;
- its weight-<=1 truncation has a **better** certified argmin than the full objective (2.411 vs
  2.661, L22) -- *the pairwise content is the harmful part*;
- and the recorded law *"MAE is unusable; the prior points the wrong way -- calibration slope
  +0.376, corr +0.283: full amplitude, wrong direction"* is exactly this, measured on the
  predictions rather than on their optimum.

**The coordinate average's accuracy comes from error cancellation across retrieved fragments -- a
statistical mechanism -- and has nothing to do with the objective.** The programme has two signal
sources that do not compose, because one of them is aimed wrongly.

### What this does to the forward plan

**Refinement cannot be the foundation.** Any architecture that moves a structure toward this
objective loses 0.56 A. The plan must **replace or repair the objective first**, and only then is
refinement -- or selection, or a quantum sampler over it -- worth anything.
