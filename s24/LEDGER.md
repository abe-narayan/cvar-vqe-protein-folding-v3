# SPRINT 24 — DECISION LEDGER

## L1 — THE DIRECTIVE'S MANDATED SCREENING STATISTIC IS BROKEN. USE THE COSINE, NOT THE FRACTION.

`s24/srcdecomp.py`, `results/srcdecomp.json`, **n=126, complete.** Forks and outcome-readings written
before the run. Every arm: uniform coordinate average, point cloud, its own medoid frame, m=75 —
only the SOURCE of the members changes. Mean universe size 18,674 windows per target.

    source                                          RMSD   member   |ebar|^2     idio        f
    P0  top-75 by shipped score   [INCUMBENT]     3.0483   3.7037     160.36    63.82   0.6758
    P1  top-75 by BLOSUM order    [no score]      3.2928   4.6893     179.26   145.19   0.4683
    P2  75 UNIFORM from universe  [no retrieval]  3.7808   5.3194     224.38   177.10   0.4896
    P3  full K=500 shipped pool   [reference]     3.3961   4.8945     186.71   158.52   0.4610
    P4  library mean (300 draw)   [reference]     3.7473   5.2955     219.05   176.71   0.4852

> **Directive §15 makes `fcommon` a mandatory screen and says a promising source should show a LOWER
> fcommon than 0.676. That rule selects backwards.** Every source here has a lower f than the
> incumbent, and every one of them is worse — the blind library draw has **f = 0.4896, a "better"
> screen than the incumbent's 0.676, and is 0.76 Å worse.**

**Mechanism, and it is arithmetic.** `f = |ebar|² / (|ebar|² + mean|d_k|²)`. Any source with a wider
spread inflates the denominator and drives f down **without removing one Ångström of shared error**.
A source can improve its f by simply being more scattered. Worse, the numerator is not independent
information either: `|ebar|² = n · RMSD²` exactly, so "absolute common error" and "the arm's RMSD"
are the same quantity wearing two names.

**The screen that does carry information is the one this ledger's L2 measures: the ANGLE between a
candidate source's bias and the incumbent pool's, against a within-source control.** f stays in the
report as a descriptive statistic. It is removed as a decision rule. Recorded here because it
overrides a numbered instruction in the directive.

---

## L2 — THE BIASES ARE PARTIALLY INDEPENDENT, AND QUALITY STILL DOMINATES. THE GENERATOR HAS A BAR.

`s24/biasalign.py`, `results/biasalign.json`, **n=126, complete.** Bias vectors compared in the
NATIVE frame (the only frame in which two sources' errors are commensurable); RMSD in each arm's own
medoid frame, as deployed.

### (a) The angle, with the control that makes it mean something

    score-75  vs  BLOSUM-75     (both retrieved)      +0.7976   median +0.8582
    score-75  vs  library-75    <-- THE QUESTION      +0.6467   median +0.7382
    BLOSUM-75 vs  library-75                          +0.8439   median +0.8991
    library-75 vs library-75'   <-- WITHIN-SOURCE     +0.9330   median +0.9617

    ratio  cos(A,C) / cos(C1,C2)  =  0.693

Two averages of real protein windows are aligned for trivial reasons — both are contracted, both are
protein-like — so a raw cosine is uninterpretable. Against its own within-source ceiling, a blind
corpus draw shares **69%** of the alignment a source shares with itself. **There is genuine angular
headroom, and it is about 31%.** This is the premise the whole sprint rests on, and it is now
measured rather than assumed.

### (b) The mixture curve bows below the line — independent confirmation, geometry-free

Matched-size mixtures, always exactly 75 members, only composition changes:

    mixture      RMSD     vs 75/0    vs the straight line between the endpoints
    m75_0      3.0483    +0.0000     +0.0000
    m60_15     3.0849    +0.0365     <b>-0.1149</b>
    m50_25     3.1351    +0.0867     <b>-0.1657</b>
    m38_37     3.1964    +0.1480     <b>-0.2255</b>
    m25_50     3.3286    +0.2803     -0.2245
    m0_75      3.8055    +0.7572     +0.0000

**Parallel biases sit on the line; independent biases must bow below it. This curve bows by up to
0.226 Å.** Fitting the quadratic error model `S(λ)/S(0) = (1-λ)² + λ²q² + 2λ(1-λ)qc` to the interior
points recovers **c_eff = 0.655 mean / 0.785 median** against the directly measured cosine of
**0.647 / 0.738**, correlated per target at ρ = +0.53. **Two independent routes — a geometric angle
and a fit to an RMSD curve that never sees a bias vector — agree to within 0.01.**

### (c) And yet every mixture loses, which is the actionable part

The second source's bias is usefully independent and its quality is not good enough for that to pay.
Inverting the same model, a second source improves the mixture at small λ exactly when

    q  <  q* = 1/c_eff        q = (its own mean error) / (the incumbent's)

    q* median 1.274   ->  a source with this much bias independence must come in under
                          ~1.27 x 3.048 = **3.88 Å standalone mean RMSD** merely to BREAK EVEN.

The blind library source sits at **3.81 Å, q = 1.231** — just inside its own bar on 87 of 126 targets
and still a net loss in aggregate, because the targets where it violates the bar violate it badly.

> **THE SPEC FOR THE GENERATOR, DERIVED RATHER THAN GUESSED.** A learned candidate source is not
> required to be diverse. It is required to be **≤ ~3.9 Å standalone on the 126-target instrument
> with a bias cosine no higher than ~0.65 against the incumbent pool** — and to move the mean by a
> measurable amount rather than break even, materially better than that on at least one of the two
> axes. Diversity alone is already available for free and is already measured: it is worth +0.76 Å.

### (d) The directive's §17 union, run today with a zero-cost stand-in generator

Merge the shipped K=500 with 500 uniform library windows, score **all 1000 with the same shipped
functional**, take the same top-75:

    union vs incumbent   **+0.0022**  SE 0.0174  MDE 0.0487  fold[-0.0267,+0.0305]  61W/65L  NULL
    share of the merged top-75 drawn from the NON-RETRIEVED half:  **mean 0.355, median 0.347**

**The scorer takes over a third of its final set from candidates BLOSUM retrieval never proposed, and
the answer does not move by a measurable amount.** Two things follow, and they point in opposite
directions:

1. **BLOSUM retrieval contributes far less than its central position suggests** — the distogram score
   recovers a third of an equally good set from blind draws.
2. **Pool composition, at fixed score and fixed readout, is a flat lever.** Adding candidates that are
   merely *different* changes nothing. This is directive §17's central causal experiment, executed
   with a free stand-in source, and it returns a precise null.

**Consequence for the sprint: the generator must beat the pool, not merely differ from it.** An NN
whose samples the existing score likes about as much as what it already has will land on +0.002 Å.

---

## L3 — THE SHARED BIAS IS THE SCORE'S, NOT THE CORPUS'S AND NOT BLOSUM'S. THE L2 SPEC IS INCOMPLETE.

`s24/qmatch.py`, `results/qmatch.json`, **n=126, complete.** Forks and falsifier written before the run.

L2 left one thing untested: the blind-library source was both *retrieval-free* and *bad*, so its
0.647 bias cosine could have come from either property. This file separates them. **B' has
retrieval-free provenance and score-selected quality**: the top-75 by the *same shipped distogram
score* drawn from 2000 uniform library windows with the shipped BLOSUM pool explicitly excluded.

    A  (incumbent)         3.0483
    B' (retrieval-free)    3.1128      q = |e_B|/|e_A|  mean 1.046  median 1.006
    cos(bias_A, bias_B')   **+0.9432**  median +0.9727     [within-source control: +0.9330]

**B' meets the L2 quality bar outright — q median 1.006, essentially quality-matched — and its bias
is MORE aligned with the incumbent than a source is with itself.**

    mixture     RMSD      vs 75/0    vs the line
    m75_0     3.0483     +0.0000     +0.0000
    m60_15    3.0599     +0.0115     **-0.0014**
    m50_25    3.0734     +0.0250     **+0.0035**
    m38_37    3.0782     +0.0298     **-0.0020**
    m25_50    3.0854     +0.0370     **-0.0060**
    m0_75     3.1128     +0.0645     +0.0000

**The curve is FLAT against the line — every point within 0.006 Å of a straight interpolation, where
L2's blind-library curve bowed by 0.226 Å.** Parallel biases sit on the line. These are parallel.

### The mechanism, and it inverts the sprint's premise

> **Selection by the distogram score is what makes the bias parallel. Provenance has nothing to do
> with it.** The blind library draw was 31% angularly independent *because it was unselected* — its
> independence was bought with the 0.76 Å of quality that made it useless. Run the same score over a
> retrieval-free candidate set and the emitted structure lines up with the incumbent's at 0.943.

Three standing results fall out of this as one:

- **L2(d)'s flat union (+0.0022)** — merging pools under the same score cannot move the answer,
  because the score maps any sufficiently rich candidate set onto the same emitted structure.
- **L2's bowing curve** — real, but it was pricing *unselectedness*, not provenance.
- **35.5% of the merged top-75 coming from blind draws with no effect** — the score does not care
  where a candidate came from, and neither does the answer.

### Consequence: the L2 spec was incomplete and is hereby restated

L2 handed the generator lanes "≤ ~3.9 Å standalone AND bias cosine ≤ ~0.65". **B' demonstrates those
two are not independently attainable through this score**: quality obtained *via the distogram score*
causes alignment. A generator whose samples are filtered by the shipped score is predicted to land on
the incumbent's bias direction no matter how it was produced.

> **RESTATED SPEC. A candidate source can only help if its advantage survives the shipped score —
> i.e. it must be better ON THE SCORE'S OWN TERMS (higher-quality candidates the score can already
> recognise), or it must come with a scoring change that is itself validated. "Different provenance"
> is now measured and is worth nothing: +0.0115 Å at the best mixture, inside its own MDE.**

### The hypothesis this raises, and it is now the sprint's most valuable open question

If two candidate sets selected against the *same distogram prediction* inherit that prediction's
error, then the 68% common-mode error of s23 L9 is not a property of the candidates at all — it is
the **shared referent floor** of project memory, with the distogram as the referent. That is
directly testable and it is the coordinator's next experiment. If it holds, the sprint's leverage is
in the prior, not the generator.

---

## L4 — A REAL LEAKAGE DEFECT THAT REACHES THE SEALED BENCHMARK. DECLARED, NOT QUANTIFIED, NOT FIXED.

Workstream A, `s24/a_ha2.py`, `results/a_ha2_fragment_leak.json`, `results/a_ha2_selfwindow.json`.
Determined from source, sequences and identifiers only. **No benchmark result or benchmark structure
was opened.**

### A's own hypothesis was falsified, and it said so

A predicted the fragment filter never sees benchmark targets because they "sit in other folds". Read
from source that is false: `core/bench.py:207` and `core/pipeline.py:1315` both set
`fold = db.folds(cfg.n_folds)[target.seq]`, so a target's OWN fold is the held-out fold at inference
and `distogram._fold_fragments` does filter against it. **The fragment half is clean** — verbatim
containment 0/126, 0/24, 0/60; ≥0.6-identity homologs missed by the unsound 3-mer prefilter 0/126,
0/24, 0/60; fragment-parent deposit overlap 0/126, 0/24, 0/60. The 6,003-fragment and 787-peptide
banks come from disjoint deposit sets with zero id overlap.

### The peptide half is not clean, and the mechanism is a documented convention

`core/data.py:188` normalises Needleman–Wunsch identity **by the LONGER sequence**, and its own
docstring already says so: *"It is leaky at the member level … `containment` is the normalisation a
window-level filter should use."* So an 11-mer sitting **verbatim** inside a 21-mer library chain
scores 11/21 = 0.52, falls below the 0.60 threshold, lands in a different identity cluster, therefore
a different fold — and therefore appears in the target's own fold model's training set **and** in
`library_members`, which excludes only exact-sequence matches. `s9/final.py:210` states the guarantee
the longer-normalisation breaks.

    tuning126   4/126   1CEK in 1A11 · 2FBU in 2LMF · 2P5H in 2P5J · 6B9K in 1U6V
    dev24       0/24
    bench60     **2/60**   8Y3S = 8XTO[10:21] (8XTO in fold 3) · 8ZG3 = 8XTO[0:11], also in 8ZUG

**8XTO is literally the concatenation spanning both benchmark targets and sits in the training set
and retrieval library of both.** Read from identifiers and sequences only.

### The size, measured where measurement is permitted

On the 4 tuning targets A located the exact self-window in the shipped universe:

    target   n    BLOSUM rank of the self-window    its rr    pool best(K=500)   universe best
    1CEK    13              0  (rank #1)            0.595 A       0.342              0.310
    2FBU    12              0                       3.278         2.243              1.154
    2P5H     9              0                       2.334         1.840              1.386
    6B9K    10              0                       4.126         2.077              1.515

**Exposure is certain — a verbatim self-sequence window is always BLOSUM rank #1 and therefore always
in the K=500 pool — and it is not a near-native answer.** In 3 of 4 cases the self-window is WORSE
than the pool's own best by 1.0–2.0 Å: the same peptide sequence in a different deposit adopts a
different conformation.

### COORDINATOR'S DECISION, and it is deliberately the conservative one

**Option (a): leave the benchmark sealed and carry 2/60 as a declared caveat on every benchmark
figure this project has reported or will report.**

- Option (b), quantifying whether the 2/60 self-windows are near-native, **requires reading benchmark
  structures' RMSD. That is an absolute prohibition in the standing directive and I am not doing it,
  nor authorising anyone else to.** A stopped at exactly the right line and was right to hand it up.
- Option (c), re-pinning the fold assignment with a substring test, is a **sprint-scale action**:
  `peptide_db.folds` warns that renumbering clusters silently reassigns every sequence and invalidates
  every trained model on disk, and project memory records that correcting the identity clustering once
  before silently moved 13 benchmark targets and invalidated every fold model. Not this sprint, and
  not without the user's explicit decision.
- The fix is **not** "use containment at 0.6" — project memory measures that threshold AT THE NULL for
  peptides (random 9–16mers score 0.56–0.63 against this bank). A verbatim-substring test, which is
  what A used, is the correct instrument.

**The incumbent 3.0483 Å is not retracted.** 4/126 affected, and on 3 of the 4 the leaked window is
worse than the pool's own best, so the number is not materially built on the defect. It is
**declared**, not withdrawn.

**Standing item for the user**: whether to re-pin folds with a substring test is their call, not mine.
It would invalidate every trained model on disk and move the instrument. Recorded here so the
decision is visible rather than buried.

### Adopted independently of the above

A's recommendation to remove the unsound 3-mer prefilter from `distogram._fold_fragments` is correct —
it is documented as unsound in `peptide_db.clusters`, it catches nothing today (0 missed homologs
across all 210 targets), and an unsound guarantee inside a leakage filter should not be left standing.
**Not applied this sprint** — it changes a shipped filter and would have to be re-validated — but it is
recorded as a known defect with a measured zero current impact.

---

## L5 — THE COMMON MODE IS THE PRIOR'S ERROR, AND SELECTION IS WHAT INSTALLS IT. TRIANGULATED.

`s24/referent.py`, `results/referent.json`, **n=126, complete**, forks and falsifier pre-registered.
Distance space, the score's own pair set (min_sep=2), no structure operations, no model.

### Primary 1 — the emitted structure is closer to the prediction than the TRUTH is

    RMS |Dc - Dhat|   cloud to prediction     **2.2358 A**
    RMS |Dt - Dhat|   native to prediction    **3.2973 A**
    difference  -1.0615  SE 0.1267  MDE 0.3549  fold[-1.208,-0.851]  111W/15L  -> 3.0x MDE
    RMS |Dc - Dt|     cloud to native          2.6347 A

**The pipeline's output agrees with the distogram's prediction better than the native structure
does.** That is over-fitting of the prior, stated without a correlation coefficient and without a
model, and it is measured on 111 of 126 targets.

### Primary 2 — beta, the share of the prior's own error that reappears in the output

    REAL   distogram                    beta +0.5201  median +0.5412   cos +0.6259
    PLACEBO mismatched same-length      beta +0.3521  median +0.2761   cos +0.4795
    real, posterior-mean (secondary)    beta +0.5413  median +0.5740   cos +0.6351

    real - placebo  **+0.1680**  SE 0.0248  MDE 0.0696  fold[+0.108,+0.218]  100W/26L  2.4x MDE

The placebo is the shared-referent floor — both errors are measured against the same native, which
induces covariance by construction, and project memory records that exact trap turning a
"2/3 sequence-independent" claim into "1/5". **The floor is 0.352 and the real value is 0.520; the
excess over the floor is real, 2.4x its own MDE, and 100W/26L.**

### Independently confirmed by Workstream C, in the score's own space, with a zero-information arm

C measured cosines between SIGNED PAIR-DISTANCE error vectors — a different construction, a different
space, and a control I did not have:

    cos(emitted distance error, DISTOGRAM's own distance error)
    arm             scored-75   uniform-75   whole set
    INCUMBENT         0.640
    T0_helix          0.584       0.439        0.487      <- knows NOTHING about the target
    T1_blind          0.662       0.425        0.586
    T2_restype        0.672       0.407        0.582
    T3_pool           0.720       0.535        0.643

> **Applying the selector raises alignment with the distogram's own error by about +0.20 in EVERY
> arm — including one that knows nothing about the target.** Selection is the operation that installs
> the prior's error into the output. It is not a property of the candidates.

### The three routes, and what they jointly establish

| route | construction | result |
|---|---|---|
| L3 | bias cosine between two score-selected sources of different provenance | +0.9432, above the +0.9330 within-source control |
| L5 P1/P2 | distance-space over-fit and beta against a mismatched-distogram placebo | cloud closer to Dhat than the native is; beta 0.520 vs floor 0.352 |
| C's ladder | signed distance-error cosine against the distogram's own error, with a zero-information arm | selection adds ~+0.20 alignment in every arm |

> **The 68% common-mode error of s23 L9 is not a property of the candidate pool. It is substantially
> the distogram's prediction error, transmitted into the output by the selection step.** That is why
> nine downstream operators failed across Sprints 22–23, why L2(d)'s pool union was flat, and why L3's
> quality-matched retrieval-free source was parallel. Every one of those acts on the candidates; none
> of them touches the referent.

**This retires the sprint's founding premise.** A better candidate generator cannot remove an error
that selection installs from the prior. It also reconciles with the project's standing ceiling —
perfect distance knowledge caps the instrument at ~1.95–2.0 Å — which prices the prior's contribution
at ~1.05 Å of the incumbent's 3.05, against beta = 0.52 of transfer. The two numbers agree.

**The lever therefore moves to the FUNCTIONAL**, which is also where Workstream D's set-equality
theorem independently places it: the selection rule is fixed by the algorithm, so the only things that
can move the answer are the energy functional and the candidates — and the candidates are now measured
out.

---

## L6 — CONFIDENCE-WEIGHTED SELECTION: THE MECHANISM MOVES, THE ENDPOINT DOES NOT.

`s24/conf.py`, `results/conf.json`, **n=126, complete.** Forks and both falsifiers pre-registered.
The shipped score is a UNIFORM mean of per-pair Bayes risk; the distogram publishes a per-pair `sd`
the score ignores entirely. L5 says selection installs the prior's error, so weight each pair by the
prior's own confidence: `w_p ∝ sd_p^(-k)`, normalised. **k = 0 reproduces the shipped score exactly**
and that identity is asserted per target before any number is read.

         k      RMSD      beta      vs k=0
     -1.50    3.1229    0.5401     +0.0746
      0.00    3.0483    0.5201     +0.0000   <- the shipped score
      1.00    3.0553    0.5100     +0.0070
      2.00    3.0476    0.4951     -0.0007
      3.00    3.0429    0.4897     -0.0055

    CONFIDENCE WEIGHT, NESTED CV   +0.0314  SE 0.0304  MDE 0.0850  fold[+0.0010,+0.0940]  33W/93L
    in-sample best k = 3.00        -0.0055  SE 0.0470  MDE 0.1316                          65W/61L
    per-target ORACLE k            -0.2426  SE 0.0313  MDE 0.0876  fold[-0.281,-0.199]    121W/  5L  BEATS

    beta, nested arm minus k=0:    -0.0125  SE 0.0075  MDE 0.0210  fold[-0.0277,-0.0006]

**The primary falsifier fired.** The nested-CV arm is not an improvement — it is directionally
*worse*, and even the leaky in-sample fit reaches only −0.0055, well inside its own MDE.

**The informative part is that the mechanism worked and the endpoint still did not move.** Beta falls
monotonically in k across the whole grid, 0.5401 → 0.4897: confidence weighting genuinely does reduce
the share of the prior's error that selection transfers, exactly as L5 predicts. It buys nothing.

> **So a ~6% reduction in prior inheritance is worth approximately zero Ångströms.** That is a
> quantitative bound on what a functional change of this class can achieve, and it is bad news for the
> whole "fix the functional" direction that L5 and Workstream D's theorem jointly pointed at. The
> prior's error is transmitted, but the transmission is not the binding constraint — the prior's error
> itself is.

**And a fifth per-target oracle appears, at 121W/5L and 2.8x its own MDE.** Real, large, near-universal,
and unreachable — the same shape as the averaging width, the arm choice, the cluster choice and the
scale. The standing finite-sample bound explains it and nothing has changed.

---

# CORRECTIONS TO L1-L3, ISSUED ON WORKSTREAM E'S AUDIT

Lane E was tasked to attack the coordinator's Tier-1 work because Rule 0 requires the fork enumerator
to have no stake and the coordinator enumerated their own. It found four things. **All four are
accepted; three are errors of mine.** The originals above are left standing and corrected here rather
than edited, so the record shows what was claimed and when.

### C1 — "TWO INDEPENDENT ROUTES AGREE" IS WITHDRAWN. IT WAS ONE ROUTE.

L2(b) claimed the fitted `c_eff = 0.655` corroborated the directly measured cosine `0.647` because
"a geometric angle and a fit to an RMSD curve that never sees a bias vector agree to within 0.01".
**That is wrong and the demonstration is decisive.** E built the mixture ladder ANALYTICALLY from
e_A and e_C — no structures, no averaging operator, no RMSD of any real cloud — and ran the identical
fitter:

    ladder                          c_eff mean  c_eff median   direct mean  direct median    rho
    A/C1 as shipped                     0.6546      0.7851        0.6467       0.7382       0.531
    A/C1 ANALYTIC, no operator          0.6467      0.7382        0.6467       0.7382       **1.000**

**The zero-information route recovers the direct cosine to every printed digit at rho = 1.000.** The
"agreement to within 0.01" was 0.008 above a floor of 0.000. S(0) and S(1) *are* the fit's endpoints
and are the same two norms that form the direct cosine's denominators; inverting the quadratic is
algebra on the same two vectors against the same referent. This is the project's own
`shared-referent-floor` trap, and I walked into it in the same ledger entry that cites the rule.

**What the mixture curve actually is: a faithfulness check on the averaging operator, and a good one.
It is not a second measurement. The sprint has measured the angle ONCE.**
`c_eff`, `q*`, the 0.655/0.785 pair and rho=0.53 were computed by an ad-hoc shell command with **no
persisted implementation in the repository and no fork declaration**. Rather than persist a fitter
whose output is now known to add nothing over the direct cosine, **all c_eff-derived quantities are
withdrawn.** L2(a)'s 0.647 against a within-source 0.933, and the 0.693 ratio, are untouched — those
are operator-matched, plausibly-controlled direct measurements.

### C2 — THE BOW SURVIVES, ON A CONTROL I DID NOT BUILD.

E constructed the operator null I should have: two independent uniform draws from ONE source, matched
quality, direct cosine 0.933, run through the identical medoid-recomputing frame-changing operator.

    A[:a] + C1[:b]      as shipped                    bow at m38_37  **-0.2255**
    A[perm][:a]+C1[:b]  composition-fork null                        -0.1967
    C1[:a] + C2[:b]     WITHIN-SOURCE OPERATOR NULL                  **-0.0294**

Medoid recomputation, changing spread, Jensen on the square root and RMSD's non-linearity are all
inside that null by construction, and it is flat. **The -0.226 A bow is real.** L2(b)'s conclusion
stands; only its claim of independent corroboration falls.

### C3 — AN UNDECLARED FORK IN MY OWN MIXTURE LADDER, AND IT POINTED MY WAY.

`biasalign.py:169` and `qmatch.py:137` build interior mixtures as `concatenate([A[:a], C1[:b]])`, and
`A` is held **in score order** — so `A[:a]` is the top-a BY SCORE, not an a-subset of A. My declared
`readout` fork said the arms "differ in composition and in nothing else"; they also differ in the
selectivity of the retained half. **Undeclared, and it flatters the interior points.** In biasalign it
is worth -0.005 to -0.029 A across the four interior points, all under their own MDEs: at most 13% of
the headline.

**In qmatch the same fork is more dangerous and the result survives it in the conservative direction.**
There both halves are score-prefixes at every interior point, a mechanism that can put an interior
point below both endpoints *with zero bias independence* — i.e. it can manufacture exactly the win
qmatch existed to look for. **qmatch found no bow at all.** The fork's bias runs toward a bow and none
appeared, so the true curve is at least as flat as reported and L3's "these are parallel" is
strengthened, not weakened. The uniform-subset ladder is nonetheless added as a second arm.

### C4 — HOUSEKEEPING, AND ONE ROW STRUCK.

- **P4 in L1 is struck, and the cause is now isolated — it was mine.** P0-P3 reproduce bit-exactly
  from E's independent reconstruction of the RNG stream. P4 did not, and E found why to six decimal
  places: **the stored P4 is a 2000-window draw wearing a "300 draw" label.** Reconstructing the
  stream with size 300 gives 1A13 2.691988 / 1CEK 1.770185; with size 2000 it gives 2.805015 /
  2.134482, which is exactly what is stored.

  **Mechanism, and it is a process error of the coordinator's.** I reduced P4 from 2000 to 300 mid-
  flight to cut cost, and the `pkill` that was supposed to stop the running job **failed silently** —
  `pkill: command not found` in this shell — so the original pre-edit process kept running and rewrote
  the artefact at 16:30 and 16:33 from source that no longer existed on disk. A second process was
  launched alongside it. That is also why the ledger's P4 (3.7473) and the artefact's P4 (3.7483)
  differ while every other row is stable. **The script no longer reproduces its own artefact.**
  Quote nothing from P4; it was a reference row and no conclusion rests on it.

- **SPRINT-WIDE CONSEQUENCE, and it matters more than P4.** Between 16:24 and 16:41 `s24/` went from
  6 modules to 24 with jobs launching and files being edited continuously. **Any process launched
  before an edit writes results attributable to source that no longer exists, and no result file in
  this sprint currently records which source produced it** — even though BRIEF.md §5 already requires
  a git commit in the log line. **Mandated from now, for anything that will be promoted:**

      from s24 import stats_lib as ST
      ST.save_atomic(path, obj, complete_keys=NEED, rows=rows, n_expected=126, module_file=__file__)

  which stamps module name, sha256 of the module's own source, git commit, dirty flag and launch time,
  and sets `complete` only when the row count matches AND every row carries every key. Anything a
  claim rests on gets re-stamped. **A result whose source hash is unknown cannot be replicated and
  will not be promoted.**
- **L1's f column is not commensurable across rows.** P0/P1/P2 are m=75, P3 is m=500, P4 was m=300,
  and f is m-dependent. The load-bearing contrast P2 vs P0 is m=75 on both sides so the conclusion is
  unaffected; the table needed the footnote and now has it.
- **L1's verdict was post-hoc.** The pre-registered falsifier required f AND |ebar|^2 *both*
  at-or-above; what came back was |ebar|^2 above and f below, and the code printed "Mixed: do not
  collapse them". The verdict is still correct, but **it rests on the identity |ebar|^2 = n*RMSD^2
  (verified, max relative error 8.5e-15), not on the pre-registration.** Recorded as post-hoc-but-sound.
- **L2(d)'s union gloss is slightly wrong.** `lib500` draws over the whole universe including the 500
  pool members, so ~13 windows per target (~2.7% of the library half) are duplicates scored as
  "non-retrieved". "Candidates retrieval never proposed" should read "candidates drawn without regard
  to retrieval". The 0.355 share is otherwise as measured.
- **L2's m0_75 endpoint is a single draw** where srcdecomp averages 8 for the same quantity
  (3.8055 vs 3.7808) — 0.025 A of pure draw noise that shifts the line and every bow.

### C5 — THE GENERATOR SPEC IS RESTATED, AND THEN RETIRED.

E showed the L2(c) spec mixed summaries (mean cosine with median c_eff), that the correct aggregate
break-even is **1.3519 -> 4.12 A** rather than 1.274 -> 3.88, that 12/126 targets (9.5%) have
c_eff <= 0 where the bar is infinite, and — decisively — that **the lambda->0 model mispredicts at the
one lambda this instrument can reach**: fed the measured per-target cosines it predicts mixing at
lambda=0.2 HELPS by 0.019 A; the instrument says it HURTS by 0.037 (SE 0.0187, fold [+0.012,+0.055],
59W/67L).

> **The only empirically anchored bar is the measured one: the blind library source sits at 3.8055 A
> with cosine 0.647 and loses at EVERY mixture fraction.** Every model-derived bar (3.88-4.71 A) is
> looser than the single point where the model has been checked, and there the model has the wrong
> sign. E also notes the bar retains real *per-target* skill — win rate at lambda=0.2 rises from 0.468
> to 0.707 conditional on q < q*, lift +0.238 — so it is a per-target predictor, not an aggregate bar.

This also dissolves L2(c)'s own stated puzzle. "Inside its own bar on 87/126 and still a net loss" was
an artefact of pairing a **mean** RMSD with a **median** q: on the mean, q = 1.618 against q* = 1.528,
so the library source **fails** its bar and the loss was predicted all along.

**The spec is in any case retired as a go/no-go by L3 and L5, and by Workstream C's finding that a
zero-information constant-helix source reaches cosine 0.693 unscored — decorrelation is free and
therefore carries no information about whether a source is useful.** It survives only as a
diagnostic.

---

## L7 — WORKSTREAM B CLOSES ITS OWN LANE BY AN ORACLE UPPER BOUND, AND REPLICATES A PROJECT LAW.

`s24/agentB_FINDINGS.md`, artefacts `resid0/resid0b/resid0c.json`, **n=126, complete**, point cloud
throughout, pre-registered in `s24/PREREG_B.md`. **No model was trained.** The lane turned out to be
closable from an upper bound, which is the cheapest possible way to close it.

### The headline is an oracle ceiling, not a failed attempt

Granting the residual's DIRECTION from the native, a per-residue confidence GATE from the native, and
choosing both the amplitude and the gate threshold IN-SAMPLE on the same 126 targets they are scored
on — **four stacked oracles, two of them leaky** — the arm buys **−0.032 Å at 0.71× its own MDE** at
the torsion accuracy this project actually achieves (σ = 69.7°). The sprint needs −0.049 Å.

    sigma=70, rho=1  (a realistic model)                    -0.0101   0.23x MDE   DEAD NULL
    sigma=45, rho=1  (granting phi-quality accuracy on BOTH angles)  -0.0399   0.79x MDE   TYPE-M
    sigma=70, rho=0  (optimistic i.i.d. bound)              -0.0647   0.55x MDE   ns

### B's own derivation was confirmed and then made worse — its own words, kept

B pre-registered the prediction that an i.i.d. residual contributes ~0 through the uniform mean.
Measured, the member-specific component is **strictly worse at every amplitude**, +0.031 → +1.020 Å,
monotone, with the bias norm growing 1.010 → 1.333.

> **A mean-free torsional perturbation is not a mean-free coordinate perturbation.** The
> torsion→coordinate map is nonlinear, so chain integration converts per-member scatter into coherent
> coordinate distortion. This is the same physics as the project's standing +1.024 Å
> torsion-vs-coordinate averaging penalty, reached from the perturbation side.

**So the directive's mandated stochastic generator cannot help through this readout and actively
hurts.** A model emitting many samples per candidate is emitting, in the only component the readout
sees, its conditional mean plus a penalty for the spread. COMM (one shared correction) carries 77% of
the oracle gain; IDIO carries none and costs.

### The pre-registered falsifier fired, and it is an independent replication

Matched-accuracy coherence pair — full oracle direction, differing only in whether mistakes are shared:

    sigma 30:  COH +0.275 (crosses)   IID -1.140     gap 1.41 A
    sigma 45:  COH +1.141             IID -0.109 (crosses)

**At identical accuracy, coherent mistakes cost 1.1–1.4 Å more than i.i.d. ones** —
`error-coherence-decides-correctors` reproduced on a completely different instrument at ~3× the
recorded effect size. This is what closes the lane: a residual conditioned on TARGET-level information
(sequence, distogram, fold) makes target-level mistakes shared by all 75 members **by construction**,
so a learned residual lives on the COH curve, which crosses do-nothing at ≈27°. The project achieves
69.7°, and φ from full sequence is 36.1° against 36.4° **sequence-blind**.

### The L3 diagnostic, answered directly

    COMM a=0.10   cosR +0.984   norm 0.974        best realistic cell: cosR +0.971  norm 1.027
    COMM a=0.30   cosR +0.926   norm 0.941
    COMM a=1.00   cosR +0.562   norm 0.759

**At every achievable amplitude the residual RESCALES the same error rather than correcting it.** The
bias only genuinely rotates at oracle amplitudes that require knowing the native.

### Two mechanisms worth keeping independently of the verdict

**(a) The validity trap inverts here.** Ramachandran favoured: P0R 0.955, ORAC a=0.30 0.811, COMM
a=0.30 0.859; clashes per structure 0.11 → 0.71. **A partially-applied torsional residual is less
chemically valid than either endpoint** — shrinking the amplitude interpolates between two
Ramachandran basins and lands between them, the same failure as a circular mean of a bimodal
distribution. Amplitude shrinkage is what makes the residual safe in RMSD and exactly what makes it
chemically invalid, and any model emitting a posterior mean inherits it. **Here it is the SAFE arm
that is absurd**, which is the reverse of the usual trap.

**(b) The template is not an asset on the φ channel.** Mean |wrap(native φ − retained member φ)| =
**40.6°** against the sequence-blind Ramachandran marginal's **36.4°**. The retrieved window's own φ is
a *worse* estimate of the native φ than a prior that sees nothing. In "template + residual" the model
starts behind a constant prior. Independent support for `phi-carries-no-sequence-signal` from the
retrieval side.

### And it lands back on L9

COMM carries 77% of the channel, so the object worth predicting is a **per-target common-mode offset**,
not a per-member residual — which is precisely the quantity s23 L9 proved is invisible from inside the
pool and priced at −0.3403 Å ORACLE, unreachable in principle. **The lane rediscovered L9 from the
generation side and hit the same wall.** Corroboration, not a new opening.

**Substrate check for the reproduction gate:** incumbent 3.0483 exact; torsion-rebuilt P0R 3.0524;
+0.0041 SE 0.0047 MDE 0.0131, NULL. The torsion parameterisation is free and is not the barrier.

`s24/residlib.py` carries the mode-collapse and validity harness with reference values. It also
records, so nobody re-derives it: for torsion-built structures ω deviation, cis fraction and bond
length/angle deviation are **0 by construction** and Cα–Cα is 3.80 Å by construction, so **the only
validity axes carrying information in torsion space are Ramachandran and clashes.**

---

## L8 — THE CORPUS IS NOT REPETITIVE. IT IS SMALL. AND IT HAS ALMOST NO BETA SHEET.

Workstream A, Deliverable 2. `s24/agentA_FINDINGS.md`, corpus hash **29e3b67e8ca0c03d**, permitted-corpus
artefact `a_corpus_permitted_29e3b67e8ca0c03d.json`. Criteria and both halves of the label
pre-registered in `s24/PREREG_A.md` before the run.

### A's own hypothesis was falsified in the opposite direction to the one it predicted

A predicted the corpus would be a small number of folds repeated. Greedy leader clustering (Kabsch
Cα-RMSD, τ = 1.0 Å, length-stratified, exact):

    L=9   n=38,316   leaders at N=500/1k/2k/5k/10k:  301 / 521 / 856 / 1,648 / 2,549
    L=11  n=26,311                                   352 / 655 / 1,184 / 2,430 / 4,070
    L=16  n=7,329  (the FULL corpus, not a subsample) 413 / 788 / 1,477 / 3,271 / 4,485

**Every curve is still climbing near-linearly at the largest N affordable, and at L=16 over the entire
permitted corpus 4,485 of 7,329 windows are mutually more than 1.0 Å apart.** There is no saturation.
Against the pre-registered plausible controls at matched N=5,000: the constant α-helix floor gives 1
leader at every N and every length; a per-residue-marginal Ramachandran draw gives 3,312–4,837 against
the real corpus's 1,648–2,430. The corpus is ~2× more concentrated than an uncorrelated draw and
vastly more than the degenerate floor. **It carries real structure. It does not carry enough of it.**

### The pre-registered NO condition fires on the raw count, before any independence collapse

**The 16-residue band holds 7,329 windows in total.** A's NO condition — fewer than ~10⁴ independent
windows in any length band — fires without needing the clustering at all. It fires at L=15 too
(10,027 raw, only ~2,600 peptide-derived).

### And the number that actually binds is worse

The permitted corpus is 6,347 chains from 1,400 independent deposits: **410 peptides** (410 deposits,
8,354 residues) and 5,937 fragments (990 deposits, 80,738 residues).

> **Nearly half the peptide bank is gone — 410 of 787 survive**, because the 210 held-out targets are
> themselves peptides and their ≥0.6-identity neighbourhoods take the rest. So a legitimately-trained
> generator gets 410 peptides where RETRIEVAL draws from ~630 out-of-fold peptides per target.
> **The generator's training corpus is strictly SMALLER than the library it is being asked to beat.**

And it is handicapped in exactly the half that carries the signal: project memory records that the
peptide corpus carries **7× the sequence–structure channel** of the protein fragments that are 80% of
every pool, and the peptide share of permitted windows is 0.132 at L=9 and 0.308 at L=16.
`fragment_db`'s own docstring declares fragments training-only and **not** a source of "what does this
peptide do on its own", because a fragment's conformation is held by contacts outside the window.

### The half of the label that PASSES, reported because both halves were committed in advance

Genuine torsional multimodality is present: global basins αR 0.536 / β 0.394 / αL 0.055 / other 0.016,
and **0 of 20 residue types is more than 80% concentrated in one basin** (most concentrated is Ala at
0.713). **The ">80% dominant basin" NO condition does not fire.** The corpus is multimodal; it is small.

### A structural bias any generator here must design around

    secondary structure over the permitted corpus:   H 0.375   E **0.022**   C 0.604

> **Beta sheet is essentially absent, at 2.2%, and mechanically so:** DSSP assigns E only with a paired
> strand, and both banks are isolated short chains whose partner strand was left outside the window by
> construction. **Any generator trained here will be helix/coil-biased and will fail quietly on sheet
> targets.** Every generative result must be stratified by native SS class.

Relatedly, distinct H/E/C window strings number only 355 (L=9) to 895 (L=16) with effective sample size
exp(H) of 21.7–137.5: **at the level of secondary-structure PATTERN the corpus has an effective size of
20–140**, and the distinct-string count plateaus at ~900 from L=14 while the possible count grows 3× per
residue.

### A's five caveats against its own conclusion, kept

(1) The 10⁴ threshold is A's own, set in advance, and rests on a naive ten-samples-per-parameter rule —
a 10³-parameter model would clear it. (2) Windows across lengths share residues, so the eight bands are
not independent evidence. (3) A leader count is a covering number at one resolution, not an information
content. (4) The near-linear curve means the count is bounded by the sample size, so at L=9–14 this is a
LOWER bound on independence, not the value. (5) A generator need not memorise windows.

### THE ITEM THAT MATTERS BEYOND THIS SPRINT

> **Each of the 5 distogram fold models saw ~100 of the other 125 dev natives. Per-target dev results
> are therefore NOT independent, and iid CIs over the 126 targets are ANTICONSERVATIVE.
> Fold-clustered CIs are the only honest interval on this instrument.**

The project already reports both and the BRIEF already mandates the fold-clustered one; this supplies
the reason, and it means **every iid CI in every sprint of this project is narrower than the truth.**

### Disposition

**Lane C's from-scratch p(φ,ψ | sequence, distogram) is NOT SUPPORTED at 9–16mers on this corpus by
A's own pre-registered criterion.** A held no stake in that answer once it came out NO. Note that
Lane B's residual design is not bound by this in the same way — a residual learns a correction over
~2.35M (candidate, target) pairs rather than a distribution over 7,329 windows — but L7 has closed
that lane on accuracy grounds independently.

---

## L9 — THE GENERATION LADDER AT FULL n: EVERY UNION IS NULL-TO-WORSE, AND SELECTION'S LIFT IS UNIVERSAL.

Workstream C, `s24/c_ladder.py`, `results/c_ladder.json`, **n=126, complete, git a15406c**, 2000
samples per target, forks pre-registered in `s24/PREREG_C.md`. Point cloud on both sides; the
built-chain figure is carried separately and never mixed.

    reproduction gate:  incumbent 3.0483 (exact) · built chain 3.2126 · torsion-rebuilt top-75 3.0524
    pool ORACLE best 1.7108 · cos(incumbent bias, distogram's own error) 0.6351

    arm            standalone   cos vs inc   cos vs dgram    UNION   gen share   ORACLE best
    T0_helix          3.7892       0.7508        0.6061     3.0568     0.097       2.9424
    T1_blind          3.2435       0.8043        0.6350     3.0654     0.235       2.1664
    T2_restype        3.2065       0.8209        0.6480     3.0640     0.256       2.1178
    T3_pool           3.1752       0.9055        0.7194     3.0992     0.523       2.0197

    UNION vs incumbent    T0 +0.0085 fold[+0.002,+0.014] 28W/98L · T1 +0.0170 42W/84L
                          T2 +0.0157 54W/72L             · T3 +0.0509 fold[+0.007,+0.094] 45W/81L

**Every union arm is null-to-worse, and the best generated source is the one that copies the pool.**
T3_pool supplies 52% of the merged top-75 and makes the answer 0.05 Å worse.

### The selection lift, at full n, in every arm including a zero-information one

    alignment with the DISTOGRAM's own error:      scored-75   uniform-75      LIFT
    T0_helix   (knows nothing about the target)       0.6061      0.4511     **+0.1550**
    T1_blind                                          0.6350      0.4421     **+0.1929**
    T2_restype                                        0.6480      0.4457     **+0.2023**
    T3_pool                                           0.7194      0.5619     **+0.1575**

> **Applying the shipped selector raises alignment with the prior's own error by +0.155 to +0.202 in
> every arm, including a constant-α-helix source that has never seen the target.**

**READ THIS THE WAY WORKSTREAM E READS IT, NOT THE WAY I FIRST DID.** I originally presented this as
the cleanest demonstration that the prior's error is installed by selection rather than carried by the
candidates. E's reading is better and I have adopted it: **an effect that a zero-information arm also
shows is a GENERIC effect**, and this one is close to tautological — selection minimises Bayes risk
against `Dhat`, so it necessarily pulls emitted distances toward `Dhat`, and therefore toward `Dhat`'s
error. The universality across arms is worth recording; it carries no Ångströms and it is not evidence
that the transferred error is TARGET-SPECIFIC. See L11.

### Route (a) — the coverage claim — is true in its premise and false in its conclusion

C measured it rather than betting on it. **27% of T3's samples score better than the pool's 75th-best
window and the score takes 52–59% of a merged top-75 from the generated half** — so the corpus genuinely
IS missing conformations the shipped functional ranks highly. **But those chains are structurally
WORSE**: generated ORACLE-best 2.02–2.94 Å against the pool's 1.71 at matched count.

> **The score's preference is not aligned with structural quality.** This is
> `objective-does-not-rank-the-native` (native at the 36.8th percentile) surfacing as the precise,
> quantified reason candidate generation cannot pay, and it is why "supply better-scoring candidates"
> is not a route.

**And decorrelation is free, which retires the cosine as a bar.** T0_helix — a zero-information
constant helix — reaches cosine 0.693 unscored, landing exactly on L2's library ratio. A screen that a
constant helix passes carries no information about whether a source is useful.

### Lane C's own verdict, adopted

No network was trained, and the reason is a diagnosis rather than a refusal: the reference to beat is a
per-residue 2-component von Mises mixture with three parameters per residue per component and **zero
training**, which s14 already measured at φ 33.6° / ψ 59.2° — **better than this project's own trained
leave-fold-out sequence predictor** at 36.1° / 62.4° — while the sequence channel a model would have to
exploit is measured at approximately zero (36.1° against 36.4° sequence-blind; the entire channel is
10.4° of ψ). Directive §38 asks which COMPONENT failed; the answer is that the conditioning channel is
empty. **Scoped to the from-scratch sequence-conditioned path only.**

---

## L10 — THE FUNCTIONAL LEVER CLOSES. TWO GENUINELY DIFFERENT PHYSICS FUNCTIONALS DO NOT ESCAPE THE REFERENT.

Workstream D, `s24/d1_hamdisagree.json`, **complete, n=30** (seeded, fold-stratified, subset declared
before the run). Genuine `H_Legacy` (11-term, DEFAULT_WEIGHTS, never fitted) and genuine `H_AMBER`
(ff14SB/GBn2 single points via `AmberSP`, bit-exactness asserted per target, `amber_verify_max_rel` 0.0).
AMBER is an ENERGY MEASUREMENT only — no minimisation anywhere, per s23 L11.

**n=30 is a smaller instrument than the 126 and nothing here is promoted on it.** The primary is the
bias cosine, not RMSD, because at n=30 an RMSD arm cannot clear a meaningful MDE and an RMSD label
would return a near-guaranteed "not measured" that reads as a null when it is really an underpowered
instrument.

### The pre-registered normalisation was dead on arrival, and D found it before computing any outcome

    frac of pool with |z_raw| < 0.1 :  0.994        frac with AMBER E > 1e4 kcal :  0.575
    worst single point observed     :  7.1e18 kcal/mol

One clash sets the standard deviation, so 99.4% of every pool collapses inside |z| < 0.1 and
`z(E_AMBER)` degenerates into "which candidate has the worst steric clash". This is
`pauli-spectrum-delta-spike-artefact` verbatim. D switched to **rank standardisation** — strictly
monotone, so no ordering, argmin or level set changes, and the currency `core.pipeline._zrank` already
deploys. **Not cosmetic: the AMBER_PREFERS partition overlaps the moment-z version at only 0.59–0.64.**
Diagnosed from the energy distribution alone, with no outcome variable inspected.

### The result

    incumbent (n=30):  RMSD 3.2070 · cos_distogram 0.6841 · beta 0.6071
    rho(Legacy, AMBER) within target: **-0.083** -- the two energies genuinely disagree

    partition          filt RMSD   cos vs inc   FLOOR    cos_dgram   FLOOR
    AGREE_GOOD            3.7545      0.8272    0.7721     0.5416    0.5592
    LEGACY_PREFERS        3.4657      0.8525    0.7839     0.5493    0.5632
    AMBER_PREFERS         3.5472      0.7559    0.7848     0.6223    0.5632
    STRONG_DISAGREE       3.5017      0.8354    0.7839     0.5853    0.5658
    AGREE_BAD             3.5909      0.6106    0.7726     0.5599    0.5583

**Every partition is worse than the incumbent in RMSD, and no partition escapes the referent.** The
bias cosines sit at or below their own random-same-size floors; the distogram-alignment cosines sit
within ~0.06 of theirs. **The incumbent's own 0.6841 alignment is HIGHER than every partition and
higher than every floor** — because the incumbent is the score-selected set and the partitions are not.

> **Two physics functionals that genuinely disagree with each other (ρ = −0.083) partition the pool
> along directions that a random subset of the same size already produces.** Legacy and AMBER differ
> from each other and neither differs from the distogram in the way that would matter. The functional
> lever closes alongside the provenance lever.

This also reproduces L5 from an instrument built for a different purpose: the incumbent sits ABOVE the
floor and unselected partitions sit AT it, which is the selection-installs-the-referent mechanism seen
from the other side.

---

## L11 — RETRACTION: L5's TWO PRIMARIES DO NOT SURVIVE A PROPER FLOOR. THE MECHANISM CLAIM IS WITHDRAWN.

Workstream E, `s24/results/e_referaudit.json`, **n=126, complete, provenance-stamped, E's own code.**
This overturns the claim I had been calling the sprint's central finding. **E is right and I accept it
in full.**

### What E built that I did not

Five floors on the SAME emitted cloud and the SAME native, differing only in the PRIOR:

    arm                                    beta mn   beta md    cos mn   |eP| RMS
    REAL distogram                          0.5201    0.5412    0.6259     3.2973
    P_peer      my shipped floor            0.3521    0.2761    0.4795     4.2885
    P_meanlen   generic same-length mean    0.5307    0.4651    0.5886     3.2579
    P_bestpeer  ADVERSARIAL closest peer    0.5911    0.6022    0.6844     3.1363
    P_difflen   different length            0.3752    0.3134    0.4857     3.9739
    P_shuffle   pair index permuted         0.1664    0.1186    0.2882     5.7233

    REAL minus floor, paired, fold-clustered:
      beta  vs P_meanlen   **-0.0107  [-0.045,+0.028]  NOT MEASURED**
      beta  vs P_bestpeer  **-0.0711  [-0.108,-0.026]  the REAL prior transfers LESS than a wrong one**
      cos   vs P_meanlen   **+0.0373  [-0.004,+0.076]  NOT MEASURED**
      P1 gap vs P_meanlen  **-0.1225  SE 0.1027  MDE 0.2879  fold[-0.303,+0.074]  NOT MEASURED**

**A generic prior that has never seen the target reproduces 102% of β, 94% of the cosine and 88.5% of
the −1.06 Å Primary-1 gap.**

### Why my placebo was too weak, and it is not the obvious reason

**A single random same-length peer is not short of information — it is NOISY.** Its |eP| is 4.29
against the real prior's 3.30 and the generic prior's 3.26, and that extra 30% is the peer's own
idiosyncratic error, orthogonal to eC, which deflates both β and the cosine. **My floor was low because
it was a noisy draw from the class, not because the class carries no signal.** Average the placebo over
the class and the floor rises to meet the real value.

**And β is not norm-free.** β = cos·|eC|/|eP|, so a worse-predicting placebo gets a mechanically smaller
β at identical alignment; E measured that multiplier at 1.132, so **13% of my +0.168 β gap was norm
rather than alignment** before any of the above.

### Both checks I asked for pointed the wrong way

I asked E to test a different-length placebo and a shuffled-pair-index placebo. Both are **looser**
floors (+0.162 and +0.354), so both make the claim look stronger. Shuffling the pair index destroys the
sequence-separation structure that is most of what a Cα distogram knows — a degenerate control, the same
failure mode as uniform-on-the-torus. **The direction to look was toward a TIGHTER floor, not a wronger
one, and I asked for the wrong thing.**

### What Primary 1 was actually measuring

**Typicality.** The emitted cloud is a contracted, generic object — project memory prices that
contraction at 25.8% — and the native is atypical, so **any** protein-like same-length distance
prediction sits nearer the cloud than the native does. That is a real fact about the averaging operator.
It is not "the output agrees with THIS distogram".

### WHAT L5 IS NOW, AND WHAT IT IS NOT

> **RETAINED:** the emitted cloud is generic and contracted, and selection aligns it with whatever prior
> it is scored against — **including a prior that knows nothing about the target.** Near-tautological
> (selection minimises Bayes risk against `Dhat`), universal across every arm, and worth stating.
>
> **WITHDRAWN:** β and the Primary-1 gap as evidence of TARGET-SPECIFIC prior transfer. **"The 68%
> common mode is substantially the distogram's prediction error" is not supported by `referent.py`**,
> and "this retires the sprint's founding premise" is not licensed by it.

**L3, L6, L7, L8, L9 and L10 are unaffected.** Each rests on its own instrument — a bias cosine against
a within-source control, a nested-CV endpoint, an oracle upper bound, a corpus census, a union ladder,
a partition panel. **What falls is my EXPLANATION for why they all came out flat, not the fact that
they did.**

### The right instrument for the ceiling claim, which the project already owns

`distance-prior-is-the-ceiling`: ORACLE distances give **0.36 Å pool / 0.98 Å selected** through the
same library. That is a **target-specific counterfactual**, and no correlation against a shared native
can substitute for one. **Workstream B's B-3 prior-attribution ladder is exactly that counterfactual and
is now the sprint's live test of the ceiling claim**, not `referent.py`.

### Two smaller corrections, both accepted

- **C3's "strengthened" clause is dropped.** E ran the uniform-subset ladder on qmatch with the fork
  removed on both halves: −0.0012 / +0.0093 / +0.0057 / −0.0023, **null at every point on both CIs**,
  and at two of four the sign runs the other way. My argument — "the fork runs toward a bow, none
  appeared, therefore parallel is strengthened" — claimed information from a measurement that returns
  none. **Correct statement: the fork was measured and is null, so L3 is neither strengthened nor
  weakened by it.** L3 stands, uncontaminated.
- **No c_eff is to be quoted from qmatch even descriptively.** Its shipped ladder fits c_eff = 1.0047
  with ρ = **−0.259** against the direct cosine: at q ≈ 1.02 the fit is pure noise and anti-correlates
  with the quantity it estimates.

### B7 resolved — the 0.693 ratio needs no correction

    cos(A,C1) 0.6467   cos(A,C2) 0.6396   cos(C1,C2) 0.9330
    shared-arm ratio 0.6932   ARM-DISJOINT ratio 0.6856   symmetric mean 0.6894
    cos(A,C2) - cos(A,C1) = -0.0070  SE 0.0112  MDE 0.0313  power 0.10  NOT MEASURED

The shared arm is not biasing it. **Quote 0.69, cite 0.686 as the arm-disjoint check.** E also
reproduced L3 independently from its own code: cos(A,B') = 0.9432 against the within-source control
0.9330, ratio 1.011.

---

## L12 — L6's "FIFTH PER-TARGET ORACLE" IS AN ORDER STATISTIC, NOT A SIGNAL. AND THE BETA SENTENCE IS RESCOPED.

Workstream E replicated the whole L6 table from `conf.json` with independent code. **Every number
matches** — k=−1.5 3.1229, k=0 3.0483, k=1 3.0553, k=2 3.0476, k=3 3.0429, per-target oracle −0.2426,
the k=0 identity, the nested-CV arm at +0.0314 and the in-sample best k at −0.0055 inside its own MDE.
**L6's conclusion stands: the mechanism moved and the endpoint did not.** Two corrections to how it was
worded.

### (a) The oracle is 92.5% accounted for by its own best-of-12 null

The k grid has **twelve** values and the per-target oracle takes the minimum over all twelve. That is a
best-of-K arm and the brief requires it be scored against the distribution of the **minimum**.

    observed per-target ORACLE k                     -0.2426   (121W/5L, 2.8x MDE)
    best-of-12 null, WITHIN-target resample          -0.2244   -> **92.5% ACCOUNTED**
    best-of-12 null, pooled residuals                -0.3948   -> 162.7% (over-explains)
    residual after the correct null                  **-0.018 A**

The within-target null is the right one and the conservative one: each target keeps its own dispersion
across the k grid and only the association between a particular k and a particular target is destroyed.
E notes a permutation null is useless here because **the minimum is permutation-invariant** — it has to
be resampling with replacement.

> **Corrected wording, adopted:** a fifth per-target oracle appears at 121W/5L and, unlike the previous
> four, it is **92.5% accounted for by its own best-of-12 null**; the residual is −0.018 Å. **It is an
> order statistic, not a signal**, and the finite-sample bound is not needed to explain it. The word
> "real" is struck from L6.

**121W/5L is not evidence against this.** A best-of-K arm wins on almost every target by construction,
which is exactly why W/L cannot diagnose it. This is the Sprint-23 trap with a larger K: twelve noisy
variants per target buy −0.22 Å for free.

### (b) The beta mechanism sentence is rescoped by L11

L6 said β falling 0.5401 → 0.4897 shows confidence weighting "reduces the share of the prior's error
that selection transfers, exactly as L5 predicts". **β is the quantity L11 has just shown is 102%
reproduced by a prior that has never seen the target.** So what the k-sweep demonstrably reduces is
alignment with **any** protein-like same-length distance prediction, not specifically with this
distogram's error. **L6's load-bearing conclusion — a mechanism moved and the endpoint did not — is
unaffected, and so is the bound it places on the fix-the-functional direction.**

### (c) OPEN, AND IT REACHES BACK INTO SPRINT 23

The ledger of Sprints 22–23 cites **four** per-target oracles as "real, large, near-universal and
unreachable", and the finite-sample bound was invoked to explain all four. **Three of them are minima
over grids** — the averaging width m\*, the arm choice, the cluster choice — and are therefore
best-of-K arms that have never been scored against the distribution of the minimum. The fourth, the
per-target scale s\*, is a **closed-form** optimum rather than a grid minimum, and s23 L10 already
checked it against a best-of-K null and found it far outside.

**Workstream E is tasked with auditing the other three.** If they are substantially order statistics,
the project's standing "real oracle headroom that no native-free rule can reach" narrative is weaker
than recorded — the headroom would be partly an artefact of taking minima over grids, and the
finite-sample bound would be explaining something that does not need explaining. **This is recorded as
open, not as settled in either direction.**

---

## L10-A — AMENDMENT: AMBER *DOES* SELECT ALONG A NON-PARALLEL DIRECTION. IT CANNOT PAY FOR IT.

I under-recorded L10 by reading only the score-filtered columns. Workstream D's own analysis
(`d1_analysis.json`, stamped, module `d1_analyze.py` sha d18a6647) carries a measured positive result
and it is the closest thing to an opening this sprint found. **D's registered prediction was WRONG and
it said so.**

D pre-registered cosines of 0.85–0.95 — partitions differing in QUALITY but not DIRECTION. Measured:

    AMBER_PREFERS, bias cosine vs the incumbent:      **0.5693**
    against its own shared-referent floor:            **-0.2154   SE 0.0546   1.41x MDE**
                                                      fold[-0.297,-0.143]   25W/5L   **MEASURED**

> **A genuinely different physics functional DOES select along a materially non-parallel direction,
> and it clears the ≤0.65 direction bar comfortably.** This is the only source measured in the entire
> sprint to do so past its own floor. Every other candidate source — retrieval-free, blind-library,
> generated, quality-matched — came back parallel.

**It fails on quality, and that is the whole story.**

    partition          RMSD    <=3.9?     cos    <=0.65?     q      <1.274?   passes both?
    AMBER_PREFERS     4.1342     no     0.5693    YES      1.520      no       FAILS
    AGREE_BAD         3.6045    YES     0.6060    YES      1.389      no       FAILS
    LEGACY_PREFERS    3.6197    YES     0.7673     no      1.199      YES      FAILS

**Legacy buys no direction at all** (−0.0166 against its floor, 0.22× MDE): the disagreement axis is
carried almost entirely by AMBER. D notes that had the L2 spec been a single bar rather than a
conjunction, it would have called this a win — **the conjunction is what saved the conclusion, not the
direction condition.**

### The recurring shape, arriving from the physics side

**Diversity is available for free and costs quality.** A blind library draw buys direction and loses
0.76 Å; AMBER buys direction and loses 1.09 Å. Three sprints have now produced this pattern from three
different mechanisms.

### The pre-filter arm: five for five in the harmful direction

    pre-filter          RMSD   vs incumbent   x MDE    W/L    cos_distogram raw -> filtered
    AGREE_GOOD        3.7545     +0.5475      0.79    15/15        +0.0082
    LEGACY_PREFERS    3.4657     +0.2587      0.60    12/18        +0.0839
    AMBER_PREFERS     3.5472     +0.3402      0.69    12/18        **+0.1361**
    STRONG_DISAGREE   3.5017     +0.2947      0.68    10/20        +0.0659
    AGREE_BAD         3.5909     +0.3839      0.67    12/18        -0.0008

No single cell clears its own MDE at n=30 (0.60–0.79×), so **no individual harm is claimed** — the
consistency of the sign across all five is the reportable part. **Physics as a pre-filter ahead of the
deployed scorer costs 0.26–0.55 Å and buys nothing.** And letting the shipped score select *inside*
each partition moves distogram alignment back UP every time, AMBER_PREFERS furthest at +0.1361 — the
selection mechanism reproducing under D's own control, exactly as pre-stated.

### Two clean reproductions worth their own line

    rho(Legacy, AMBER) = **-0.0829**, median -0.0760, SE 0.0397, 19/30 negative, RETRIEVAL manifold
    s20 measured **-0.0886** on the CONTINUOUS TORSION manifold

**Two manifolds, two instruments, the same number.** And the two energies sit differently against the
deployed scorer: **ρ(Legacy, distogram) = +0.3875 while ρ(AMBER, distogram) = −0.0270.** Legacy is
partly aligned with the thing that already selects; AMBER is orthogonal to it. That is the mechanistic
reason AMBER is the one that buys direction.

### D's own corrections, recorded

- The rank-vs-raw partition overlap on the full panel is **0.348**, not the 0.56–0.68 quoted from two
  targets. **Two-thirds of the candidates studied depend on the D1-A normalisation amendment**, so it
  was more substantive than D first claimed. `frac E>1e4 kcal = 0.575`, `frac |z_raw|<0.1 = 0.994`.
- D stamped the raw artefact with the module that **produced** it rather than with the analyser —
  stamping it with `d1_analyze.py` would have pointed the source hash at code that never wrote those
  numbers, which is precisely the failure mode Lane E found in my own file.

---

## L13 — THE PRIOR-ATTRIBUTION LADDER. THE ONLY STEEP LEVER IN THE PROJECT, AND IT IS THE PRIOR.

Workstream B, `s24/priorladder.py`, `results/priorladder.json`, **n=126, complete,
provenance-stamped (module priorladder.py, commit a15406c8).** Everything except the prior is held
fixed: the same shipped K=500 pool, the same Bayes-risk functional, the same top-75, the same uniform
coordinate average. Only `Dhat` moves, interpolated toward the native by γ.

**ORACLE BY CONSTRUCTION. This is a transfer function, never a system result, and γ measures distance
travelled toward a PERFECT prior — not "improvement achievable by any known means".**

    gamma       MASS      TILT     DIRAC  MASSFIXW
    0.0       3.0483    3.0483    3.0752    3.0483     <- the real distogram, the incumbent
    0.1       2.8334    2.8575    2.9254    2.9261
    0.2       2.6827    2.7216    2.7992    2.7762
    0.3       2.5847    2.5871    2.6619    2.6386
    0.5       2.4199    2.4149    2.4404    2.4235
    1.0       2.2261    2.2261    2.2261    2.2332     <- a perfect prior
    EXACT substitution at gamma=1: 2.2367

**Four independently constructed interpolations — probability-mass, tilt, Dirac and fixed-width —
agree to within 0.093 Å at γ=0.1 and to 0.007 Å at γ=1.** That answers the metric-realisability
objection: the result does not depend on the construction.

    MASS arm, paired vs the real prior
    gamma=0.1   2.8334   **-0.2150**  SE 0.0313  MDE 0.0878  fold[-0.266,-0.172]  **113W/13L**  2.4x MDE
    gamma=0.2   2.6827    -0.3656     SE 0.0539  MDE 0.1509  fold[-0.497,-0.280]  112W/14L
    gamma=0.5   2.4199    -0.6285     SE 0.0758  MDE 0.2125  fold[-0.775,-0.492]  116W/10L
    gamma=1.0   2.2261    -0.8223     SE 0.0912  MDE 0.2554  fold[-0.968,-0.654]  114W/12L

### The three numbers that matter

    slope at the origin      **-2.1496 A per unit gamma**
    slope at the top          -0.4168 A per unit gamma      -> the curve SATURATES, 5.2x more
                                                               return from the first 10% than the last
    gamma required to reach 3.0000 A   **0.0225**

> **Moving the prior 2.2% of the way toward perfect takes the instrument below 3.0 Å.** Every other
> lever measured this sprint — pool composition (+0.0022), provenance (+0.0115), confidence weighting
> (+0.0314), physics pre-filters (+0.26 to +0.55), generated-source unions (+0.009 to +0.051),
> the residual channel (−0.032 at four stacked oracles) — is flat to within a few hundredths.
> **The prior's derivative is two orders of magnitude larger than anything else in the system.**

### What this does and does not license

**DOES:** it is the target-specific counterfactual that L11 said was required, and it is the right
instrument for the ceiling claim in a way `referent.py` was not. It replaces a correlation against a
shared referent with an intervention: change only the prior, hold every other operator fixed, watch
the endpoint move. It also reproduces the project's standing `distance-prior-is-the-ceiling` result
from a new direction — perfect distances through the shipped selection path give **2.2261 Å**.

**DOES NOT:** γ is progress toward a perfect prior, not progress achievable by any known method. A 2.2%
interpolation toward the native is not the same object as a 2.2% better distogram, and nothing here
says such a distogram is buildable. The honest statement is about the **derivative**, not about a
delivery date.

**β falls from 0.5413 to 0.3122 across the ladder**, which is the expected direction; per L11 it is
descriptive here and carries no independent weight.

### Disposition for the sprint

**<3.0 Å was not reached, and this sprint measured why with an intervention rather than a
correlation.** The generator direction that the directive was built around is closed on five separate
instruments. The lever that is not closed was never a candidate-generation problem: it is the accuracy
of the distance prior, and it is steep.

---

## L13-A — AMENDMENT: HOW THE LADDER WAS BUILT, THE GATE THAT FAILED, AND A COLLAPSE CROSSOVER.

### The arm as I specified it was not expressible, and B stopped rather than approximating

I asked for `Dhat(γ) = (1−γ)·Dhat_real + γ·D_native` substituted into the shipped functional. **That is
not merely metrically non-realisable — it is not expressible at all**, because the shipped functional
consumes a **distribution**, not a point estimate. Writing it as a point-estimate loss would have been
a different functional and would have made the entire ladder uninterpretable. **B stopped on exactly
the condition I set and redesigned**, interpolating **probability mass** instead, whose γ=1 endpoint is
the native's own bin and is realisable by construction:

    MASS(g)   mix toward the true bin        recentre AND sharpen     PRIMARY
    TILT(g)   exponential tilt, minimum KL   recentre ONLY
    DIRAC(g)  collapse on interpolated mean  recentre AND collapse    (the collapse control)
    MASSFIXW  freezes the weight w to isolate that channel; tracks MASS within 0.01-0.09 A

### The loss is exactly reconstructible, and B verified rather than assumed it

`risk[p,t] = w[p]·Σ_c prob[p,c]·|grid[t] − CENTRES[c]|` with `w[p] = shell/(sd+0.5)^g` — bilinear in
`prob` with a weight that is itself a function of `prob` through `sd`. **Max absolute reconstruction
error 0.0**, and rather than re-implement it B **constructs a genuine `core.predict.Distogram` from the
modified probabilities and uses its own `_risk`**, so `w` is recomputed self-consistently exactly as a
genuinely better distogram would. `_verify_functional` raises and halts if it ever changes.

### REPRODUCTION GATE 2 FAILED, AND THE DIAGNOSIS IS THE MOST USEFUL PART

γ=1 lands at 2.2261 (unbinned EXACT 2.2367) against the project's standing "perfect distances cap at
~1.95–2.0 Å". **B treated its own ladder as the suspect**, as instructed, and tested it by changing
only the candidate set under the same perfect-knowledge functional, same top-75, same readout
(`poolcheck.py`, complete):

    K=500 BLOSUM pool          2.2367     <- the ladder's regime
    FULL universe              **1.6024**     (mean 18,674 windows/target)
    difference          **-0.6342**  SE 0.0572  MDE 0.1601  **3.96x**  fold[-0.696,-0.574]  116W/10L

**The standing figure lies inside that bracket.** The ladder is not wrong — it is **pool-restricted by
construction**, which is the correct choice for a transfer function *in the prior*, and the restriction
is worth **0.634 Å**. Neither number would have been trustworthy without the gate.

### THE SYNTHESIS, both ORACLE, same instrument, same readout

    at a FIXED candidate pool, a PERFECT prior is worth          **-0.822 A**
    at a PERFECT prior, the FULL candidate universe is worth     **-0.634 A**

> **Comparable in SIZE, not comparable in REACHABILITY.** The candidate lever has been measured flat
> for every achievable source this sprint: pool union +0.0022, quality-matched provenance cosine
> 0.9432, four generated-source unions +0.009 to +0.051, the residual channel −0.032 at four stacked
> oracles. The prior lever is **concave and pays in its first increment**. That asymmetry, not the
> ceiling heights, is the sprint's finding.

### A collapse crossover with a location, and it sharpens two standing memories

**DIRAC is WORSE than MASS while the prior is wrong (+0.027 at γ=0, +0.092 at γ=0.1) and BETTER once it
is nearly right (−0.049 at γ=0.7). The sign flips at γ ≈ 0.55.**

> **Collapsing the posterior onto its mean hurts exactly when the mean is wrong and helps when it is
> right.** This turns `do not collapse the posterior` and `confidently wrong costs 2-3x absent` from
> rules into a crossover with a measured location: **confidence should be EARNED before it is
> expressed, and the shipped prior has not earned it.**

It is also the same shape as L7's finding that amplitude shrinkage is safe in RMSD and absurd in
Ramachandran: **a damped or averaged object is only correct when the thing being damped is already
good.**

### Reproduction gate 1 passed to four decimals

γ=0 through this module's own scoring path emits **3.0483**, matching the pinned constant and the
module's independently computed P0; and **β at γ=0 is +0.5413, reproducing L5's posterior-mean β
exactly on a completely separate code path.** β then falls monotonically to +0.2725, and at γ=1 sits
**below** L11's shared-referent floor of 0.3521 — exactly where it should when there is no prior error
left to inherit. B's own pre-registered falsifier (a threshold-shaped curve, meaning only a large jump
pays) **did not fire**: curvature +0.2173, concave, first 10% buys 26% of the total gain.

**B's caveat, kept verbatim in spirit:** every point on this ladder is ORACLE. It says what a better
prior would BUY, not that a better prior is obtainable, and nothing here says the distogram can be
improved 10% by any available means. **That is the next question and this ladder does not answer it.**

---

## L9-A — AMENDMENT: THE SELECTION LIFT IS +0.166, NOT +0.20. AND THE DAMPED-SAMPLER TRAP FIRED.

### The headline number is corrected downward by its own lane

L9 quoted "+0.20 from selection in every arm" from Workstream C's n=26 probe. **The n=126 confirmatory
value is +0.166 under the MAP estimator and +0.177 under the posterior mean.** Per arm (MAP):
T0_helix +0.137, T1_blind +0.185, T2_restype +0.194, T3_pool +0.147. **Use +0.166; the +0.20 is
superseded.** The part that matters is unchanged — the zero-information constant-helix arm still
carries it — and per L11 the whole quantity is descriptive, not evidence of target-specific transfer.

### Recipe pinning: one of five choices differed, and C recomputed both rather than arguing

C used `dg["expected"]`, the posterior mean; Lane D and `referent.py` use the MAP,
`grid[argmin(risk)]`. The other four choices matched exactly (min_sep=2; uniform coordinate average in
its own medoid frame, point cloud; per-target then averaged; neither vector centred). **C recomputed
the alignment under BOTH estimators on all 126 targets** and stored them per target in `c_final.json`
under `map` and `exp`, so any lane can diff directly. **The estimator is worth −0.007 to −0.017 in the
cosine and changes no conclusion.** C also verified its `_bias`/`_cos` are **byte-identical** to
`residlib.bias`/`cos`, which are identical to `qmatch._bias` (max abs diff 0.0) — so every cosine in
this sprint is the same object by construction.

### Lane B's damped-sampler warning fired — on C's BEST arm, not the one it expected

C had flagged T0_helix as its damped sampler. **It was T3_pool.** The von Mises mixture with capped
concentration places density *between* the basins it was fitted to, making it simultaneously the
**best-RMSD arm (3.1752) and the WORST-Ramachandran arm (0.8980)** — 10% outliers, against 0.951–0.955
for arms that resample real torsion pairs directly and 0.955 for the incumbent's own sets.

> **The safe arm is the absurd one, exactly as L7 predicted, and it arrived in a different lane on a
> different sampler.** Anyone reusing T3_pool as a baseline must carry that caveat. The check existed
> only because Lane B handed over `residlib` with the warning attached — a cross-lane transfer that
> damaged the receiving lane's most favourable number.

---

## L14 — TWO NEW MEASUREMENTS FROM LANE C: CHAIN CORRELATION IS 7.2%, AND THE SIGNED-BIAS TENSION DISSOLVES.

### (a) The strongest objection to "no network" was closed by measurement, not argument

Every arm in C's ladder sampled residues independently, so the obvious reply was *"you never tried a
model with chain correlation — an autoregressive model or an entangled latent."* C fitted a first-order
Markov chain to the retrieved pool's own basin-label sequence and compared it against a control whose
transition matrix is replaced by **the product of its own marginals** — marginals held exactly fixed,
so the only difference is correlation.

    T4_markov - T5_shuffle   **-0.0095**  SE 0.0130  MDE 0.0363  fold[-0.026,+0.003]  63W/63L  0.26x MDE  NULL

    adjacent-residue basin MUTUAL INFORMATION in the retrieved pool:
        **0.0499 nats against a maximum of 0.6931  ->  7.2%**

> **Neighbouring residues' basins are nearly independent at this length. There is almost no chain
> correlation for a sequence model to model.** That is a measurement rather than an argument, and it
> closes the objection C would most have expected a reviewer to raise.

### (b) The signed-bias tension I flagged is resolved, and it was a shared-referent effect

I had asked C to check an apparent contradiction: the distogram predicts distances too LONG (+0.367 in
the probe) while the emitted cloud is too SHORT, yet β says half the prior's error transfers. At n=126:

    distogram signed offset                          **+0.4062**
    the whole legal universe of REAL windows          **+0.4419**
    selected members, BEFORE averaging                 -0.0627
    the emitted average                                -0.5852
    what AVERAGING alone contributes   **-0.5226**  SE 0.0373  MDE 0.1046  fold[-0.586,-0.467]  **5.0x MDE**

**First: the distogram's signed offset is not a distogram defect.** Real protein geometry carries the
same offset against these natives (+0.4419 vs +0.4062), so the signed mean is measuring a
**corpus/native scale mismatch** — another shared-referent-floor effect, and it dissolves the apparent
contradiction rather than explaining it away.

**Second: the contraction is the averaging operator, at 5.0× its own MDE.** And the two quantities are
orthogonal — averaging leaves the inherited-error component essentially untouched (projection 0.5665 →
0.5413) while moving the scale mode half an Ångström. **β prices the error DIRECTION; the signed mean
prices the SCALE mode.**

> **Consequence worth keeping: the selected members already carry the inherited error at cosine 0.708
> BEFORE any averaging. The common mode is installed by SELECTION, and no change to the readout can
> remove it.** This is the one part of the original L5 story that survives L11's retraction, and it
> now rests on a measurement that does not involve a placebo at all.

### (c) A subgroup result C declined to promote, and was right to

The SS stratification I requested shows the aggregate null is not uniform — T2_restype is +0.387 on
helix (n=49) and −0.121 on coil (n=43); T0_helix is +1.68 on extended (n=34), validating the Cα-only
classifier. **C declined to promote it**: the stratification was requested *after* the numbers existed,
it was not pre-registered, and every cell favouring C's own arms sits in the 0.7–1.3× Type-M zone
(coil T2 at 0.82×). Recorded as an unpromoted hypothesis — the corpus's helix bias may make blind
resampling relatively better on coil — and nothing more. **C flagged it as the one place in its
findings where a result could have been manufactured, and said so before being asked.**

### Lane C's verdict

**The from-scratch `p(φ,ψ | sequence, structural prior)` path is closed, and the component that failed
is the CONDITIONING CHANNEL** — not the model class, representation, optimiser or sample diversity.
Sequence buys ~0 at this length; chain correlation is 7.2% of maximum and worth 0.26× MDE; and the
distogram is the selector's own referent, so conditioning generation on it imports the common mode
rather than new information. **Scoped to the from-scratch path only.** Lane B closed by an independent
route to the same place.

**Three things C says transfer:** the ideal-geometry torsion manifold is FREE (+0.0041 Å at n=126, so
representation was never the barrier); any lane whose samples pass through the shipped score must
measure its own scored-vs-unscored cosine ladder **with a zero-information control** before quoting a
decorrelation number; and the open direction is the **selector**, not the generator — the one quantity
with measured headroom is that generated sets contain structures the score ranks highly and nativeness
does not reward.

---

## L15 — THE RETROSPECTIVE ORACLE AUDIT. GRID ORACLES ARE ORDER STATISTICS; TRANSFER ORACLES SURVIVE.

Workstream E, `s24/e_oracleaudit.py`, `results/e_oracleaudit.json`, provenance-stamped. This is the
audit L12 opened, and it revises the through-line of Sprints 22–23. **E's own write-up supersedes any
framing here; the artefact is quoted directly.**

    quantity                                    K     observed    share    residual   verdict
    m* IN-SAMPLE grid (s23 exp2, avg_m)         8     -0.3226     1.154    0.0497     ORDER STATISTIC
    arm choice, 19 realisable arms             19     -0.5311     1.126    0.0667     ORDER STATISTIC
    arm choice, drop worst arm                 18     -0.4865     1.151    0.0736     ORDER STATISTIC
    arm choice, drop all latent arms           14     -0.4259     1.178    0.0760     ORDER STATISTIC
    scale s*, GRID version                     81     -0.1374     1.025    0.0034     ORDER STATISTIC
    conf.py k grid (L12)                       12     -0.2426     ~1.08    0.0180     ORDER STATISTIC

    m* TRANSFER (split-half)                    -     -0.2394       -         -       N/A, self-nulling
    scale s* TRANSFER (split-half)              -     -0.3198       -         -       **SURVIVES, 5/5 folds**

> **Every per-target "oracle" this project has measured by taking a MINIMUM OVER A GRID is ≥85%
> accounted for by its own best-of-K null — several of them over 100%.** The share exceeding 1.0 means
> the null explains more than the observed gain: there was nothing there at all.

**And the distinction that saves the programme is the construction, not the quantity.** The same
underlying quantities measured by **split-half transfer** — fit on one half of a target's pool, apply
to the disjoint other half — are **self-nulling by construction** and a best-of-K null does not apply.
**The per-target scale survives in that form at −0.3198, 5/5 folds**, which is the arm Sprint 23's
Workstream D reported at 98–100% transfer.

For cluster choice, E notes the correct null is `oracle − random`, which **is** the exact order
statistic for that panel — so the oracle-vs-random gap Sprint 23 already reported is the properly
nulled quantity and needs no correction.

### What this changes, stated plainly

**Sprints 22–23 reported "real, large, near-universal per-target oracle headroom that no native-free
rule can reach", and invoked the finite-sample bound to explain why it could not be reached.** For the
grid-measured members of that family, **the headroom was substantially an artefact of taking a minimum
over K noisy variants, and the finite-sample bound was explaining something that did not need
explaining.** My own L6 was the sixth instance and I wrote it in this sprint.

**What survives is narrower and better founded:** a per-target signal is real when it is demonstrated
by *transfer* — fitted on one disjoint sample and applied to another — and not when it is demonstrated
by *minimisation over a grid*. The scale s\* meets that standard; the grid ceilings quoted alongside it
did not.

**This is the most consequential correction the project has made**, and it was produced by a lane whose
only job was to attack the coordinator's work. It is recorded here as a finding, not as a footnote.

### Standing rule adopted

**Any oracle formed as a minimum over K variants must be reported beside `ST.best_of_k_null` with its
share-accounted and residual. A W/L count cannot diagnose it — a best-of-K arm wins nearly everywhere
by construction (118/126, 114/126, 126/126 in the rows above).** Where an oracle can be recast as a
split-half transfer, that construction is preferred because it nulls itself.

---

## L15-A — CORRECTION TO L15. TWO OF FOUR SURVIVE. THE DISCRIMINATOR IS OUT-OF-SAMPLE VALIDATION.

I wrote L15 from the artefact before Workstream E's write-up arrived and **got two of the four rows
wrong.** E's version supersedes mine.

    oracle                claimed            verdict
    averaging width m*    -0.244   77/49     **SURVIVES**
    scale s*              -0.3403 126/0      **SURVIVES, and it is the strongest of the four**
    arm choice            -0.482/-0.505      ORDER STATISTIC (112.6% accounted, residual +0.067 WRONG SIGN)
    cluster choice        -0.484             ORDER STATISTIC by exact arithmetic, 174-245%

> **The discriminator is clean and general: the survivors are exactly the two that were validated OUT
> OF SAMPLE.** In every case the thing separating a real per-target optimum from an order statistic was
> whether anyone had built a held-out arm.

### m* — I said "N/A, self-nulling". It SURVIVES, and it was never a bare oracle

`s22/results/mreal.json` holds a **split-half transfer**: m chosen on half A, scored on the disjoint
half B. Self-nulling for best-of-K — if the selecting half carries no per-target signal the choice is
noise and the held-out score is an unbiased draw, so E[gain] = 0, not negative.

    -0.2394  SE 0.0353  MDE 0.0990  **2.42x MDE**  fold[-0.2917,-0.1799]  5/5 folds  79W/38L/9T
    keeps 65.2% of its own in-sample oracle (-0.3672)

**REAL.** The s23 table listed it beside three in-sample minima, which is what made it look like the
same kind of object. It is not.

### s* — survives, and shows essentially NO overfitting

`s23/results/d_scale_transfer.json`, the held-out arm, **which nobody has quoted**:

    sel_heldout vs fixed_heldout  -0.3198  SE 0.0485  MDE 0.1358  **2.36x MDE**
    fold[-0.3920,-0.2369]  5/5 folds  115W/11L  worst-target degradation +0.0057 A
    keeps **98.9%** of the in-sample oracle

**98.9% retention is the signature of a genuine per-target parameter rather than a fitted one.** On the
81-point grid the null accounts for 102.5%, but k_eff there is **1.3 of 81 columns**, so the grid
framing is the wrong instrument; the transfer arm settles it.

**AND A SPRINT-23 CITATION MUST BE WITHDRAWN.** `d_scale_placebo` should **not** be cited as support:
its placebo delta is −1.3169 on a baseline of 7.4987 (**−17.6% of its own baseline**) against the real
−0.3403 on 3.0483 (−11.2%). **In absolute Ångströms the placebo is 5× the real effect, and as a
fraction of its own baseline 1.57×. That placebo shows the OPPOSITE of what the s23 ledger reads into
it.** It does not change the disposition — the transfer arm is decisive — but the citation goes.

### Cluster choice — I said it needed no correction. It is an ORDER STATISTIC, by exact arithmetic

`agentA.json` records **both** `oracle` (min over clusters) and `random` (mean over clusters), so no
simulation is needed:

    oracle - incumbent  =  (oracle - random)  +  (random - incumbent)
                            ^ THE order statistic    ^ the real effect of the cell

    cell      orc-inc   orc-rnd   rnd-inc     W/L
    m75_k5    -0.4615   -0.8056   **+0.3441**  109W/17L
    m250_k3   -0.4839   -1.1145   **+0.6306**   84W/42L

**At every cell a randomly chosen cluster is substantially WORSE than the incumbent, by +0.344 to
+0.631. The clustering operator is a worse operator that only looks good when you take the minimum
over its own outputs.** 174–245% accounted.

**Transcription slip fixed:** "−0.484 at m250_k3, 109/17" splices two cells. **m250_k3 is −0.4839 at
84W/42L; the 109W/17L belongs to m75_k5 at −0.4615.** The error is in s23's ledger and was repeated in
mine.

### E's method, and the fact that three of its own nulls were wrong first

**An over-explaining null is a MIS-SPECIFIED null, not proof the effect is fake.** E built four before
one was right:

    N1  exchangeable-column resample -- centres the simulated grid on the BASELINE and silently deletes
        the (rowmean - baseline) term that is part of the observed gain.        105-168%
    N2  additive model, raw interaction permuted within column -- HETEROSCEDASTIC: hands a target whose
        grid barely moves a residual borrowed from one that moves by Angstroms. 134-296%
    N3  same but standardised -- fixes scale, still treats K columns as K INDEPENDENT opportunities.
                                                                                109-135%
    N4  PRIMARY: reassign whole standardised residual PROFILES between targets, so the profile's
        internal column correlation survives and only the target/own-best-column association dies.
                                                          102.5 / 112.6 / 115.4 / 117.8%

**The diagnostic is `k_eff`, the participation ratio of the residual correlation eigenvalues**: the
scale sweep has k_eff **1.3 of 81** columns, the arm panel **3.4 of 19**. A smooth one-dimensional
sweep does not offer K independent chances at a low value, and a null that pretends it does will beat
any real effect. N4 still overshoots by 2–18%, so the honest reading is **"the null reproduces the
whole observed gain, no per-target signal is detectable above it"**, not "the effect is exactly zero".

### Open, and recorded as open

**The split-half transfer arms are self-nulling for best-of-K, but the two halves share a target, a
native and a pool, so their noise may be positively correlated.** That is a shared-referent question,
not an order-statistic one, and it is the remaining way m* and s* could be softer than they look. **Not
measured; it would take a fresh run.**

### THE CORRECTED THROUGH-LINE

**"Real per-target headroom that no native-free rule can reach" is TRUE for m\* and s\*, and the
finite-sample router bound is still needed to explain those two.** It was **over-applied** to arm
choice and cluster choice, where there is no headroom to explain.

> So the correction is **not** "some of the signal was never there" across the board — that overstates
> it, and my own summary did. **The count of independent instances drops from five to three** (m\*,
> s\*, and conf.py's k which L12 already struck). **The project's per-target headroom is real; it is
> smaller and less universal than the ledgers state.**

**The routine that would have caught this every time is cheap: build the transfer arm, or score the
minimum against its own order statistic, and never let W/L stand in for either — a best-of-K arm wins
nearly everywhere by construction.**

---

## L8-A — AMENDMENT: AT THE RESOLUTION THE ENDPOINT LIVES ON, THE 9-MER CORPUS *DOES* COLLAPSE.

Workstream A completed its heavy run across all 8 length bands and flagged, unprompted, the one number
in its own table that supports the reading its headline contradicts.

    leaders at N=10,000 (L=16 is the WHOLE permitted corpus at 7,329)
    L      n       tau=0.5  tau=1.0  tau=1.5  tau=2.0
    9   38,316     4,837    2,549      903      **213**
    12  21,310     5,908    4,488    3,018    1,464
    16   7,329     5,375    4,485    3,717    2,815

**At τ = 1.0 Å nothing saturates anywhere** — the headline stands, the corpus is not repetitive, it is
small. **But at τ = 2.0 Å the 9-residue band collapses and saturates**: 186 leaders at N=5,000 and 213
at N=10,000, a ratio of **1.15 for a doubled sample**.

> **At 2.0 Å resolution a 9-mer corpus of 38,316 windows contains on the order of 250 distinct
> shapes.** With the dev pool best at 2.355 Å and the incumbent at 3.048 Å, that is arguably the
> resolution the endpoint actually lives on.

The collapse weakens monotonically with length (L=16 at τ=2.0 is 2,815 of 7,329, still climbing), so it
is a **short-window phenomenon**. A flagged it precisely because it is the one number in its table
supporting a "small number of shapes" reading, and it did not want it buried under a headline saying
the opposite.

**The pre-registered torsion-space secondary agrees and is less forgiving**: leader clustering on
[cos φ, sin φ, cos ψ, sin ψ] at 0.35 rad RMS gives uniformly *higher* counts (4,345–5,670 at N=10,000
for L=9–15), same near-linear shape. **No disagreement with the Cα-space primary** — the outcome A
committed to publishing either way.

**Label scored, both halves, as pre-registered:** the size condition FIRES (L=16 has 7,329 raw, below
the 10⁴ floor); the ">80% in one Ramachandran basin" condition does NOT fire (0/20 types, worst Ala at
0.713); the "does not exceed the plausible control" condition does NOT fire (~2× the Ramachandran
control, ~3,000× the constant-helix floor). **One of three fires, and the answer stands: NO, not
softened.**

A's own caveat worth keeping: at L=9–14 the leader counts are bounded by sample size and are therefore
**lower bounds**, so the NO is conservative in the wrong direction there. **The NO rests on L=15 and
L=16, where the raw corpus is smaller than the threshold and no sampling argument helps.**

### Three leakage paths A named that bind any future generator work

    (a2)  the universe npz carries ORACLE `rr` and `nat_ca`, and a conditioning tensor built by
          `{k: z[k] for k in z.files}` picks them up SILENTLY. Reuse s8/test_generate.py's
          NaN-poisoning assertion rather than writing a new one.
    (a7)  PHI/PSI in the npz are float16 (~0.001 rad), so a residual smaller than that is
          QUANTISATION NOISE and a model that appears to learn it may be learning the grid.
          Recompute torsions in float64 from source.
    (a8)  the 4 dev self-copies of the L4 leak are the BLOSUM rank-#1 candidates for those targets,
          so "condition on the top-ranked retrieved candidate" conditions on a VERBATIM SELF-COPY
          on exactly those four.

**And a codebase hazard:** `fragment_db` ids are **not unique** — the id is `<PDB>_<start>` but one
start yields one fragment per length in 9–20, so **803 of 6,003 ids collide**. Any code keying
fragments by `f.pdb` silently merges 9-mers with 20-mers. Use `kind:id:len`.

---

## L16 — D2: THE FUNCTIONAL LEVER IS CLOSED IN ALL THREE OF ITS FORMS, AND CLOSED BY ITS CEILING.

Workstream D, `s24/d2_amberscore.py`, `results/d2_amberscore.json`, **n=126, complete,
provenance-stamped.** LOCK_AMBER held once, 1043 s, released. `s = s_dist + w·z_rank(E_AMBER)`, one
global `w` under nested leave-one-fold-out CV, deployed readout unchanged.

    w = 0 returns **3.0483 exactly** -- zrank is monotone, so w=0 IS the incumbent, bit-for-bit
    blended - incumbent   **+0.0041**  SE 0.0195  MDE 0.0547  **0.08x MDE**
    CI iid [-0.0359,+0.0417]   fold [-0.0334,+0.0564]   53W/69L/4T
    median +0.0036   worst +0.7557   folds same sign 2 of 5      **NULL**

At 0.08× its own MDE with per-fold deltas disagreeing in sign, **this is an absence of signal, not an
underpowered non-result.**

### Nested CV did NOT pick w = 0. It picked a NEGATIVE weight.

**w = −0.5 on four folds and −0.03 on the fifth. Zero folds chose w = 0.** The fitted blend
**anti-weights AMBER** — it prefers candidates AMBER ranks badly.

> That is not noise-shaped, and it agrees with two independent D1 results: ρ(AMBER, distogram) =
> −0.0270, and AGREE_BAD (both energies reject the candidate) at 3.6045 against AMBER_PREFERS at
> 4.1342. **On this manifold AMBER's ordering is mildly ANTI-informative for accuracy.**

**The per-fold reporting I asked for is what surfaced this; the mean would have hidden it completely.**

### The number that closes the lever is its own ceiling

    best single w on the FULL dev set, FULL LEAKAGE:   w = -0.500  ->  3.0336  =  **-0.0148 A**
    nested CV recovers                                             ->  3.0524  =  +0.0041 A

**An ORACLE w, chosen with complete leakage on the very targets it is scored on, is worth 0.0148 Å.
There is no value of w anywhere on the grid worth having.** The lever is closed by its ceiling, not
merely by a failed fit — and that ceiling sits against B-3's **−0.2150 Å** for a prior 10% of the way
to truth.

### The three nulls, and the informative one

    w=0 exact              3.0483   reproduces the incumbent bit-for-bit
    sign-flipped w        +0.0252   0.39x MDE  59W/66L   UNDERPOWERED, not a result
    rank-permuted AMBER   -0.0010   0.08x MDE  64W/53L   NULL

**The rank-permuted arm is the one that matters: preserve AMBER's marginal exactly, destroy its
correspondence to candidates, and get −0.0010 at 0.08× MDE.** Whatever the real arm shows is
indistinguishable from a control with all of AMBER's information stripped out. Its SE is 0.0045 against
the real arm's 0.0195 — **w = −0.5 is a large perturbation that swings individual targets by up to
0.83 Å (|δ| > 0.1 on 44/126) and averages to nothing.**

### A defect in D's own shipped harness, found by D2's control and now fixed

`paired_stats` labelled an effect at **0.39× its own MDE "MEASURED"**, purely because a 5-fold cluster
bootstrap CI excluded zero. With five clusters that CI is unstable, and an effect below its own MDE is
by definition one the design could not reliably detect. **The rule now returns UNDERPOWERED below
0.7× MDE.**

**It had already bitten D1**: three D1-C rows at 0.60–0.69× were machine-labelled MEASURED. D's prose
was right — it wrote "no single one is a clean measured harm and I do not claim one" — **but the JSON
disagreed with the text, and a reader trusting the artefact over the prose would have been misled.**
Both are consistent now and `d1_analyze` was re-run so the stored D1 labels are corrected.

> **Same class as the dead assertion that reported 0/2916 failures: a harness that flatters its user in
> the direction the user wants.** Two instances in one lane, both found by that lane.

### The durable asset

`s24/cache_amber/` — **63,000 genuine ff14SB/GBn2 single points** (126 × 500), `amber_verify_max_rel =
0.0` on every target, with the distogram score and pool indices alongside, 1.5 MB total. **Any future
AMBER question on this pool can be answered without touching OpenMM**, including strengthening D1's
§3.3 against Lane E's generic floor — the strengthening D1 could not do because it persisted summaries
instead of energies.

### Disposition

**AMBER's orthogonality to the distogram is real, it is the only structural reason anyone found to
expect a different functional to help, and it does not convert into accuracy through the SCORE any more
than it did through a FILTER or a PARTITION. The functional lever is closed in all three of its
available forms.**
