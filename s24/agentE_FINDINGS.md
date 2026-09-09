# SPRINT 24 — WORKSTREAM E. STATISTICS, REPLICATION AND ADVERSARIAL AUDIT.

Running log. Newest audit at the top of each section. Every number here was produced by Lane E's
own code (`s24/e_bowaudit.py`, `s24/stats_lib.py`) from the persisted artefacts, never by
re-running the claiming lane's reduction.

---

## AUDIT 1 — THE COORDINATOR'S TWO TIER-1 FILES. 2026-09-08 16:24–16:40.

Targets: `s24/srcdecomp.py` + `results/srcdecomp.json` (L1), `s24/biasalign.py` +
`results/biasalign.json` (L2). Instrument: `s12/instrument.py`, seeds `s15/seed.py`.
Lane E's own artefact: `s24/results/e_bowaudit.json`, `e_bowaudit.log`, n=126 complete.

### VERDICT IN ONE LINE

**L1 stands. L2(a) stands. L2(b) — the bow — SURVIVES a measured operator null and is NOT an
artefact. L2(c) — "two independent routes agree" — is WITHDRAWN: the two routes are the same
route, and the agreement is a shared-referent floor I measured at ρ = 1.000 exactly.
L2(d) — the spec handed to Lanes B and C — is internally inconsistent and needs restating.**

---

### (a) RE-ENUMERATION OF THE SIX FORK AXES, BY A PARTY WITH NO STAKE

I read both files from source. I did not find the Sprint-22 shape (three unstated forks all
pointing the builder's way, 89% of the effect). The undeclared forks I did find are real, are
listed below, and the largest of them is worth **≤13% of the headline bow and does not clear its
own MDE**. That is the honest answer: the coordinator's enumeration was substantially complete.

#### `srcdecomp.py` — forks the docstring does not declare

| # | axis | undeclared fork | direction | priced |
|---|---|---|---|---|
| S1 | **readout / normalisation** | The declared readout is "the uniform mean of exactly m=75 members for P0–P2, so the arms differ in SOURCE and in nothing else." **P3 is m=500 and P4 is m=300.** `f_common` and `member_rmsd` are both functions of m (`S_common` falls with m, `S_idio` does not), so P3's and P4's `f` are **not commensurable with P0–P2's**, yet L1 prints all five in one column. | inflates the apparent generality of "every source has a lower f" | L1's load-bearing contrast is P2 vs P0, **both m=75**, so the conclusion is unaffected. The TABLE is. |
| S2 | **THE LABEL** | Three different statements of what P4 is: docstring line 31 "**all** universe windows, averaged"; inline comment "**2000** is far past convergence"; code `min(N, **300**)`. | **not cosmetic — see S5** | the persisted numbers are the **2000** version; the code and both labels say 300 |
| S3 | **null** | P2 averages the statistic over `NDRAW=8` independent draws; every other arm is a single deterministic set. The alternative not taken (pool 8×75 = 600 windows, or report the single-draw sd as the arm's own noise) is not named. `P2_sd` is recorded but never reported. | reduces P2's estimator noise only, not the arm's — the estimand is right | not load-bearing |
| S4 | **the pre-registered reading** | The docstring's declared falsifier is `|ebar|²(P2) ≥ |ebar|²(P0) **and** f(P2) ≥ f(P0) → H-lib`. Measured: `|ebar|²` **above**, `f` **below** — the code's own verdict branch printed **"Mixed: report both, do not collapse them."** L1 nonetheless reports a definite verdict. | L1's verdict is a **post-hoc** reading of a pre-registered *mixed* outcome | the post-hoc argument is nonetheless **sound and checkable**: `\|ebar\|² = n·RMSD²` is an exact identity (I verified it, max relative error 8.5e-15), so f's numerator carries no information the arm's RMSD does not. L1's conclusion survives on arithmetic, not on the pre-registration. **Say so in the ledger.** |
| **S5** | **reproducibility — ISOLATED, AND IT GENERALISES** | **The persisted artefact was not written by the persisted script.** P0/P1/P2/P3 reproduce bit-exactly from my own reconstruction of the RNG stream. P4 does not — and the cause is now nailed: `srcdecomp.py` on disk draws `min(N, 300)`, but the numbers in `results/srcdecomp.json` are reproduced **to 6 decimal places** by `min(N, 2000)` (1A13 2.805015, 1CEK 2.134482, exact). A long-running process launched **before** the 16:06 edit is still writing results from the pre-edit source. | the artefact, its script, its own `ARMS` label ("300 draw") and the ledger all disagree | **P4 is a 2000-window arm, not 300.** Anyone re-running `srcdecomp.py` today gets a different P4. And S1 is worse than stated: the f column is m = 75, 75, 75, **500**, **2000**. |

`srcdecomp.json` was `complete:false` with 90/126 rows at 16:23 while a re-run was in flight, and
healed to `complete:true` 126/126 by 16:30 (and again at 16:33). **The completion flag worked
exactly as designed** — a reader in that window was told not to trust the file. Recorded as the
discipline functioning, not as a defect.

> ### SPRINT-WIDE HAZARD RAISED BY S5 — EVERY LANE, NOW
>
> A Python process loads its source **once, at launch**. Between 16:24 and 16:41 this repo went
> from 6 to 24 modules in `s24/`, with jobs launched and files edited continuously. **Any process
> launched before an edit keeps writing results attributable to source that no longer exists on
> disk**, and nothing in the result file records which source produced it. That is exactly how
> `srcdecomp.json`'s P4 row came to be a 2000-window arm sitting under a "300 draw" label.
>
> `BRIEF.md` §5 already requires **git commit** in the log line and **no result file in
> `s24/results/` currently carries one.** Minimum fix, three lines, and `stats_lib` will grow a
> helper for it: stamp every result file with the **git commit**, the **sha256 of the module's own
> source** (`hashlib.sha256(open(__file__,'rb').read()).hexdigest()[:16]`), and the **launch
> timestamp**. Until that exists, no artefact in this sprint can be attributed to a script, and
> a re-run is the only way to know what a number means.

**L1 REPLICATED INDEPENDENTLY** (my code, `stats_lib.compare`, point cloud both sides):

    P2 library-75 vs P0 incumbent, RMSD  +0.7325  SE 0.1142  MDE 0.3200  effect/MDE +2.29
                                         fold CI [+0.5643,+0.8683]  5/5 folds same sign  37W/89L
                                         worst degradation +4.5523 (8TXS)   WORSE
    P2 vs P0, f_common                   -0.1862  SE 0.0254  MDE 0.0711  effect/MDE -2.62
                                         fold CI [-0.2323,-0.1314]  90W/36L   "BETTER"

**The screen and the endpoint point in opposite directions on the same pair, past both MDEs, on
5/5 folds. L1's central claim is confirmed.**

#### `biasalign.py` — forks the docstring does not declare

| # | axis | undeclared fork | direction | priced by Lane E |
|---|---|---|---|---|
| **B1** | **readout / composition — THE BIG ONE** | Every interior mixture is `concatenate([A[:a], C1[:b]])`. `A = Wpool[argsort(score)[:75]]` is **held in score order**, so `A[:a]` is the **top-a by score**, not an a-subset of A. The interior arms therefore differ from the endpoints in composition **and** in the selectivity of the retained score half, while the declared `readout` fork asserts size-matching makes them differ "in composition and in nothing else". The alternative not taken — a uniform a-subset of the same 75 — is named nowhere. | **toward the builder**: it adds bow | **measured. −0.005 to −0.029 Å of the bow, i.e. ≤13% of the −0.226 Å headline, and NOT MEASURED at every λ (all four \|effect\| < own MDE).** |
| B2 | null | The `m0_75` endpoint is a **single** 75-window draw, where `srcdecomp` averages 8 for the same quantity (3.8055 here vs 3.7808 there — 0.025 Å of pure draw noise). The line's right endpoint, and therefore every bow, inherits it. | either way | ≤0.025·λ Å |
| B3 | readout | Mixtures take **nested prefixes of the same C1**, so the library half at b<75 is a subset of the endpoint's set and its mean is a noisier estimate of the source mean. | **against** the builder: inflates interior points, **understates** the bow | conservative; the bow is if anything a lower bound on this axis |
| B4 | functional | `lib500 = rng.choice(N, 500, replace=False)` over the **whole** universe, which **includes the 500 shipped pool members**. The §17 union can therefore contain duplicates of retrieved windows, and `union_lib_frac` scores a duplicate as "non-retrieved". Expected overlap ≈ 500·500/18,674 ≈ 13 windows/target. | slightly inflates the 0.355 | ~2.7% of the library half; the BRIEF's gloss "candidates retrieval never proposed" is not exactly what was measured |
| **B5** | **all six** | **`c_eff`, `q*`, the 0.655/0.785 fit and ρ=0.53 are computed by code that is not in the repository.** `biasalign.py`'s `report()` prints no c_eff, no q*, no ρ. The single most-quoted derived number in the sprint has **no persisted implementation and no fork declaration at all** — no statement of the estimator (closed-form LS in squared-ratio space? interior points only? weighted?), no λ-grid fork, no alternative not taken. | unknown until persisted | **I reproduced all four numbers exactly** (mean 0.6546, median 0.7851, ρ 0.531, per-target median q 1.2309) with an unweighted closed-form LS over the interior points, so the arithmetic is right. **The forks were still never declared. Persist the fitter.** |
| B6 | normalisation | The λ grid (0, 0.200, 0.333, 0.493, 0.667, 1) is non-uniform and the fit weights interior points equally in squared-ratio space. Not declared. | small | not isolated |
| B7 | shared referent | `cos(A,C1)` and the control `cos(C1,C2)` **share the arm e_C1** as well as the native referent, so the ratio 0.693 is not a ratio of two independent estimates. The within-source control is still the right idea and still the best available; the caveat belongs in the ledger. | unknown | not isolated (needs cos(A,C2), queued) |

**URGENT CONSEQUENCE FOR `s24/qmatch.py`, WHICH IS RUNNING NOW.** Line 137 repeats B1 —
`concatenate([A[:a], B[:b]])` — and in qmatch **B′ is itself top-75 by the shipped score**, so
**both** halves get the "keep the best" prefix at every interior point while the endpoints keep all
75. In biasalign only one half did. An interior point that beats both endpoints is exactly what
qmatch is looking for, and that fork is a mechanism that can produce one without any bias
independence. **Add the uniform-subset ladder as a second arm before reading qmatch's result.**
It is a two-line change (`perm = rng.permutation(75); A[perm][:a]`, likewise for B′) and it is the
difference between "mixing wins" and "the prefix wins".

---

### (b) ATTACK ON THE BOW. IT SURVIVES.

The ledger reads the −0.226 Å bow as proof of non-parallel bias. Convexity of a norm makes that
true **in the ideal model** — `|(1−λ)u + λv|` is convex and touches the chord iff u ∥ v — but the
deployed operator is not the ideal model: `_avg` recomputes the medoid and the superposition frame
at every mixture, the mixed set's spread changes with composition, and the sub-blocks are not
scaled copies of their parents. I built three nulls and measured each on all 126 targets.

    ladder                                                       bow at m38_37   bow/S(0)
    A[:a] + C1[:b]        AS SHIPPED                                  -0.2255     -0.1511
    A[perm][:a] + C1[:b]  COMPOSITION-FORK NULL (uniform a-subset)    -0.1967     -0.1446
    C1[:a] + C2[:b]       WITHIN-SOURCE OPERATOR NULL                 -0.0294     -0.0063
    analytic |(1-l)eA + l eC1|, NO OPERATOR                           -0.2367     -0.1167

**The operator null is decisive.** `C1[:a] + C2[:b]` mixes two independent uniform draws from one
source: matched quality (q≈1.013), prefixes of a random draw so B1 cannot act, and a directly
measured cosine of 0.933. It runs the **identical** medoid-recomputing, frame-changing, spread-
changing operator. Its bow is **−0.029 to −0.047 Å (bow/S(0) = −0.006 to −0.009)**, and its fitted
`c_eff = 0.979` against its own direct cosine of **0.933** — the operator recovers **slightly LESS**
bow than the geometry implies. **The operator does not manufacture bow; if anything it destroys a
little.** The shipped ladder's bow/S(0) is 10–24× the operator null's.

The composition fork B1 is worth −0.005/−0.010/−0.029/−0.028 Å at the four interior points
(paired, per target): |effect| < its own MDE at every one of them. It is a real undeclared fork
pointing the builder's way, and it is **13% of the headline at most**.

Jensen on the square root, changing spread, and RMSD's non-linearity in the coordinates are all
inside the within-source null by construction, and that null is flat.

> **THE OBSERVED −0.226 Å BOW SURVIVES. Attack (b) fails. The ledger's L2(b) claim is correct as a
> statement about the geometry.** The correction it needs is not to the bow — it is (c).

### (c) THE TWO ROUTES ARE ONE ROUTE. WITHDRAW THE CORROBORATION CLAIM.

L2(b) says: *"Two independent routes — a geometric angle and a fit to an RMSD curve that never
sees a bias vector — agree to within 0.01."* They are not independent, and the fit does see the
bias vectors.

`S(0) = |e_A|/√n` and `S(1) = |e_C|/√n` **are** the fit's endpoints — the same two norms that are
the direct cosine's denominators — and any interior `S(λ)` is, up to the operator, just
`|(1−λ)e_A + λe_C|/√n`. Inverting `S(λ)²/S(0)² = (1−λ)² + λ²q² + 2λ(1−λ)qc` for `c` is therefore
**algebra on the same two vectors**, against the **same native referent**, not a second measurement.

I measured the floor exactly as `shared-referent-floor` prescribes: build the ladder analytically
from e_A and e_C with **no candidate structures, no averaging operator and no RMSD of any real
cloud**, and run the identical fitter on it.

    ladder                              c_eff mean   c_eff median   direct mean   direct median      rho
    A/C1  as shipped                        0.6546         0.7851        0.6467         0.7382    0.531
    A/C1  uniform subset (B1 removed)       0.6514         0.8086        0.6467         0.7382    0.632
    C1/C2 within-source null                0.9787         0.9762        0.9330         0.9617    0.100
    A/C1  ANALYTIC, no operator             0.6467         0.7382        0.6467         0.7382    1.000

**The no-information route recovers the cosine to every printed digit with ρ = 1.000. That is the
floor. The claimed "agreement to within 0.01" is 0.008 ABOVE a floor of 0.000.**

What the shipped fit actually contains beyond the floor: the operator moves the mean c_eff by
**+0.008** and drops the per-target correlation from **1.000 to 0.531**. Read correctly, the
mixture curve is a **faithfulness check on the averaging operator** — it says the operator is
unbiased to +0.008 in the mean while adding a lot of per-target noise. It is a good check and worth
keeping. **It is not a second measurement of bias independence, and the sprint has measured the
angle exactly once.**

L2(a) is untouched by this: the direct cosine 0.647 against the within-source control 0.933 is a
single, real measurement with a plausible, operator-matched control, and the ratio 0.693 stands
(subject to B7).

### (d) THE SPEC. IT NEEDS RESTATING, AND TWO LANES ARE WORKING TO IT.

**Defect 1 — the spec mixes summaries.** The two halves of the sentence handed to Lanes B and C are
different summaries of the same quantity:

    "bias cosine <= ~0.65"   is  MEAN(cos_A_C)     = 0.6467
    "q* = 1.274 -> <=3.9 A"  is  1 / MEDIAN(c_eff) = 1/0.7851 = 1.2737

Used consistently, the λ→0 bar is anywhere in **3.88 – 4.71 Å**:

    1/mean(c_eff)     1.5277  ->  4.657 A        1/median(c_eff)     1.2737  ->  3.884 A
    1/mean(cos_A_C)   1.5464  ->  4.714 A        1/median(cos_A_C)   1.3547  ->  4.129 A

**The spec as issued is the strictest of the four.** It may reject a generator that the model says
would help.

**Defect 2 — `q* = 1.274` is not the median q\*.** It is `1/median(c_eff)`. **12 of 126 targets
(9.5%) have `c_eff ≤ 0`**, where the bar is infinite; the median of the actual per-target q\* over
the 114 targets where it is finite is **1.2115**. The heavy tail (mean 5.7 with a max in the tens of
thousands) is not a nuisance to be medianed away — it is 12 targets on which the model says *any*
second source helps. Report `frac(c_eff ≤ 0) = 0.095` beside the median or the number is misleading.

**Defect 3 — and this is the one that matters — an aggregate bar cannot be built by summarising
each side of a per-target inequality separately.** The sprint's endpoint is `mean_t S_t`. Its
derivative is

    d/dλ mean_t S_t(λ) |_0  =  mean_t [ S_t(0) · (q_t·c_t − 1) ]

so with a common q the aggregate break-even is `q_bar = mean(S0) / mean(S0·c)` = **1.3519 →
4.121 Å**, not 1.274 → 3.884 Å.

**Defect 4 — the model mispredicts the only place it has been tested, and the sign flips.** Fed the
**directly measured** per-target cosines, the model predicts the aggregate at λ=0.2:

    lambda   model prediction   instrument   model - instrument
    0.000            3.0483        3.0483             +0.0000
    0.200            3.0298        3.0849             -0.0551
    0.333            3.0769        3.1351             -0.0582
    0.493            3.1851        3.1964             -0.0112
    0.667            3.3556        3.3286             +0.0270

**The model says adding the library source HELPS the mean by 0.019 Å. The instrument says it HURTS
by 0.037 Å** (SE 0.0187, MDE 0.0524, fold CI [+0.0123,+0.0550], 59W/67L). Per-target sign agreement
is 94/126 (75%). The λ→0 break-even condition that generates q\* is therefore **not validated at any
λ this instrument can actually reach.**

**Defect 5 — what the bar IS good for.** Per target the bar has genuine skill: base win rate at
λ=0.2 is 0.468; conditional on `q < q*` it is 0.707, **lift +0.238**. So `q < q*` is a usable
per-target predictor. It is not an aggregate bar.

> **RESTATED SPEC, LANE E'S RECOMMENDATION.**
>
> 1. **The only empirically anchored bar is `< 3.81 Å standalone mean`.** The blind library source
>    sits at 3.8055 Å with cosine 0.647 and **loses at every mixture fraction measured**. That is a
>    measurement, not a model. Every model-derived bar (3.88–4.71 Å) is *looser* than the one point
>    where the model has been checked, and at that point the model has the wrong sign.
> 2. **Quote the range, not a point.** "The λ→0 model bar is 3.88–4.71 Å depending on the summary;
>    the aggregate-endpoint-correct version is 4.12 Å; the instrument has falsified the model's
>    aggregate sign at the smallest λ it can reach."
> 3. **Keep the two-axis deliverable** the BRIEF already asks Lane C for — standalone mean **and**
>    bias cosine, both on the 126-target instrument, point-cloud basis. That is right and unchanged.
> 4. **Do not let a lane target 3.9 Å as a pass mark.** A generator at 3.9 Å with cosine 0.65 is
>    predicted by a model that mispredicts this exact case. The honest instruction is: *get under
>    3.8 Å standalone and hold the cosine at or below 0.65, and then MEASURE the mixture ladder —
>    with the uniform-subset null (B1) — rather than trusting the bar.*

---

## TASK 2 — `s24/stats_lib.py` IS LIVE. ALL LANES USE IT.

    from s24 import stats_lib as ST
    r = ST.compare(new, old, folds=ST.pinned_folds(pdbs), names=pdbs, label="...")
    print(ST.fmt(r))

One call returns, in one block: paired mean **and median**, SE, **MDE = 2.8016 × SE for that
comparison**, **effect/MDE**, iid CI95 **beside** the 4000-resample fold-clustered CI over the 5
pinned folds, per-fold means and folds-same-sign, W/L/T, **worst-target degradation with the target
named**, p90 degradation, post-hoc power, **Type-M exaggeration factor with a flag at 0.7–1.3×
MDE**, and a concentration block whose drop-top-10 statistic is scored **against a uniform-effect
null** rather than a bare threshold (`median-vs-mean-is-the-free-warning` records that a raw
threshold misfired here once).

`VERDICT` is deliberately strict: an effect must exceed **its own** MDE **and** have the **fold** CI
exclude zero. An iid-only result is not a result.

Also in the module, and mandatory where they apply:

* **`best_of_k_null(values, k, minimise=True)`** — the distribution of the **order statistic**, not
  the mean. Any arm that reports "the best of K" is scored against this.
  **`best_of_k_accounted(arm, baseline, values, k)`** returns `share_accounted`; Sprint 23's
  "−0.077 Å oracle" scored 101% on exactly this.
* **`argmin_tied(score, outcome)`** — mean outcome over the **argmin set**. `np.argmin` on a tied
  signal reads the array's order, which on a sorted pool is the oracle order.
* **`retrodesign(effect, se)`** — power, Type-S, Type-M.
* **`pinned_folds(pdbs)`** — reads the fold assignment from the instrument. Never recompute folds.
* **`provenance(__file__)`** and **`save_atomic(path, obj, complete_keys=..., rows=..., n_expected=...,
  module_file=__file__)`** — added after S5. Atomic tmp + `os.replace`; `complete` set only when the
  row count matches **and** every row carries every key in `complete_keys`; stamps module name,
  **sha256 of the module's own source**, git commit, dirty flag and launch time. Mandatory
  sprint-wide for anything that will be promoted.

`python s24/stats_lib.py` runs a selftest.

---

## AUDIT 2 — L5 (`referent.py`), THE SPRINT'S HEADLINE. 2026-09-08 16:45–16:55.

Lane E artefact: `s24/results/e_referaudit.json` + `.log`, n=126 complete, **provenance-stamped**
(`source_sha256 61aee97c1f830a1d`, commit `a15406c82245`). Script `s24/e_referaudit.py`.

### VERDICT: NEITHER PRIMARY IS SPECIFIC TO THE REAL DISTOGRAM.

Five floors, all evaluated against the **same** emitted cloud `eC` and the **same** native `Dt`,
differing only in the prior:

    arm                                    beta mn   beta md    cos mn   |eP| RMS
    REAL distogram                          0.5201    0.5412    0.6259     3.2973
    P_peer      the SHIPPED floor           0.3521    0.2761    0.4795     4.2885
    P_meanlen   generic same-length mean    0.5307    0.4651    0.5886     3.2579
    P_bestpeer  ADVERSARIAL closest peer    0.5911    0.6022    0.6844     3.1363
    P_difflen   different length            0.3752    0.3134    0.4857     3.9739
    P_shuffle   pair index permuted         0.1664    0.1186    0.2882     5.7233

Real minus floor, paired, fold-clustered (positive = L5 surviving):

    beta   P_peer +0.1680 [+0.108,+0.218] 2.41x MDE  SURVIVES
           P_meanlen -0.0107 [-0.045,+0.028]         NOT MEASURED
           P_bestpeer -0.0711 [-0.108,-0.026]        REAL TRANSFERS LESS THAN A WRONG PRIOR
    cos    P_peer +0.1464 [+0.075,+0.217] 1.89x MDE  SURVIVES
           P_meanlen +0.0373 [-0.004,+0.076]         NOT MEASURED
           P_bestpeer -0.0585 [-0.095,-0.019]

**A prior that has never seen the target reproduces 102.1% of beta and 94.0% of the cosine.**

**Primary 1 fails the same way.** `P1 gap = RMS|Dc−Dhat| − RMS|Dt−Dhat|`:

    REAL -1.0615 | P_meanlen -0.9390 | P_bestpeer -1.2017 | P_peer -0.8019 | P_shuffle -0.3358
    real - P_meanlen  -0.1225  SE 0.1027  MDE 0.2879  fold[-0.3025,+0.0743]  3/5 folds same sign
                              power 0.22  Type-M 2.13   NOT MEASURED

**88.5% of "the emitted structure is closer to the prediction than the truth is" is reproduced by a
prediction carrying no target-specific information.** P1 measures **typicality**: the emitted cloud
is a contracted, generic object (project memory prices that contraction at 25.8%) and the native is
atypical, so any protein-like same-length distance prediction sits nearer the cloud than the native.

### WHY THE SHIPPED PLACEBO WAS TOO WEAK — the transferable part

Not lack of information. **Noise.** One random peer's `|eP|` is **4.29** against the real prior's
3.30 and the generic prior's 3.26; that extra 30% is the peer's own idiosyncratic error, orthogonal
to `eC`, and it deflates both beta and the cosine. The floor was low because it was a *noisy draw
from the class*, not because the class carries no signal. **Average the placebo over the class.**
Separately, beta is not norm-free — `beta = cos·|eC|/|eP|` — so a worse-predicting placebo gets a
mechanically smaller beta at identical alignment; measured norm-only multiplier **1.132** for
`P_peer`, i.e. 13% of the +0.168 beta gap is norm rather than alignment before any of the above.

### BOTH REQUESTED CHECKS POINTED THE WRONG WAY

`P_difflen` (+0.162) and `P_shuffle` (+0.354) are **looser** floors than the shipped one, so they
make the claim look stronger. Shuffling the pair index destroys the sequence-separation structure
that is most of what a CA distogram knows — that is close to a degenerate control, the same failure
mode as uniform-on-the-torus. **The direction to look was a TIGHTER floor, not a wronger one.**

Nuisance checked: 24/126 shipped placebo peers share the target's fold; beta 0.3181 same-fold vs
0.3601 cross-fold. A cross-fold peer's model may have seen this target, which would make the floor
*more* informative and the excess conservative. Not the explanation.

### WHAT THIS DOES AND DOES NOT OVERTURN

It does **not** show the prior is not the ceiling — the project already has far better evidence for
that (`distance-prior-is-the-ceiling`: ORACLE distances give 0.36 Å pool / 0.98 Å selected through
the same library — a **target-specific counterfactual**, which is the right instrument). It shows
that `referent.py`'s two primaries cannot carry "the 68% common mode is substantially the
distogram's prediction error" and cannot license "this retires the sprint's founding premise".
L5's third route — C's ladder showing selection adds ~+0.20 alignment in **every** arm including the
target-blind helix — read straight, supports the generic reading: **an effect a zero-information arm
also shows is a generic effect.** It is also near-tautological (selection minimises Bayes risk
against `Dhat`, so it pulls emitted distances toward `Dhat` and toward `Dhat`'s error) and carries
no Ångströms.

---

## AUDIT 3 — L6 (`conf.py`). 2026-09-08 16:56.

Every number in L6's table replicated from `conf.json` with my own code: k=−1.5 3.1229, k=0 3.0483,
k=1 3.0553, k=2 3.0476, k=3 3.0429, per-target oracle −0.2426, k=0 identity holds. **L6's
conclusion — the falsifier fired, the endpoint does not move — is correct and honestly stated.**

**But the "fifth per-target oracle" is not real.** The k-grid has **twelve** values and the oracle
takes the minimum over all twelve, so it is a best-of-12 and the brief requires the distribution of
the **minimum**:

    observed per-target ORACLE k                          -0.2426   (121W/5L, 2.8x its own MDE)
    best-of-12 null, WITHIN-TARGET resample               -0.2244   -> 92.5% ACCOUNTED
    best-of-12 null, pooled residuals across targets      -0.3948   -> 162.7% (over-explains)
    residual after the conservative null                  -0.018 A

The within-target null keeps each target's own dispersion across the grid and destroys only the
association between a particular k and a particular target. A **permutation** null is useless here
because the minimum is permutation-invariant — it has to be resampling with replacement. The
121W/5L is not evidence against the null: a best-of-K arm wins on nearly every target *by
construction*, which is precisely why W/L cannot diagnose this class.

Interpretive caveat inherited from Audit 2: L6's mechanism sentence ("confidence weighting reduces
the share of the prior's error that selection transfers") rests on **beta**, which Audit 2 shows is
102% reproduced by a target-blind prior. What the k-sweep demonstrably reduces is alignment with
*any* protein-like same-length prediction. L6's endpoint conclusion is unaffected.

---

## AUDIT 4 — C3 (the qmatch fork) AND B7 (the 0.693 ratio). 2026-09-08 16:50.

Lane E artefact: `s24/results/e_audit2.json` + `.log`, n=126 complete, provenance-stamped.

**C3 — the coordinator's L3 reasoning does not hold, but the conclusion survives.** I ran the
uniform-subset ladder on qmatch with the score prefix removed on **both** halves:

    point    shipped MINUS fork-removed
    m60_15   -0.0012  SE 0.0106  MDE 0.0297  fold[-0.0151,+0.0099]  NOT MEASURED
    m50_25   +0.0093  SE 0.0142  MDE 0.0399  fold[-0.0125,+0.0331]  NOT MEASURED
    m38_37   +0.0057  SE 0.0156  MDE 0.0437  fold[-0.0261,+0.0413]  NOT MEASURED
    m25_50   -0.0023  SE 0.0141  MDE 0.0395  fold[-0.0283,+0.0287]  NOT MEASURED

The fork is **null at every point on both CIs**, so L3's "A and B′ are parallel" stands uncontaminated.
But the ledger's argument — *the fork runs toward a bow, none appeared, therefore parallel is
strengthened* — claims information from a measurement that returns none, and at two of four points
the sign is the other way. Correct statement: **the fork was measured and is null; L3 is neither
strengthened nor weakened by it.** Also: never quote a c_eff from qmatch even descriptively — the
shipped ladder's fit gives c_eff 1.0047 with **ρ = −0.259** against the direct cosine (at q ≈ 1.02
the fit is pure noise and *anti*-correlates with what it estimates); fork-removed gives 0.9822,
ρ 0.501; the analytic floor gives 0.9432, ρ 1.000.

**B7 — no correction needed.** `cos(A,C1)/cos(C1,C2) = 0.6932` (quoted), arm-disjoint
`cos(A,C2)/cos(C1,C2) = 0.6856`, symmetric mean 0.6894. The contrast that tests the shared arm,
`cos(A,C2) − cos(A,C1)`, is **−0.0070, SE 0.0112, MDE 0.0313, fold [−0.0259,+0.0140], power 0.10 —
NOT MEASURED**. Quote 0.69, cite 0.686 as the arm-disjoint check. Byproduct confirming L3 from my
own code: `cos(A,B′) = 0.9432` against the within-source control 0.9330, ratio **1.011**.

---

## AUDIT 5 — THE RETROSPECTIVE ORACLE AUDIT. SPRINTS 22–23. 2026-09-08 17:05–17:25.

Lane E artefact: `s24/results/e_oracleaudit.json` + `.log`, provenance-stamped. Script
`s24/e_oracleaudit.py`. Sources: `s22/results/mreal.json`, `s22/results/routerdata.json`,
`s23/results/agentA.json`, `s23/results/gscale.json`, `s23/results/d_scale_closedform.json`,
`s23/results/d_scale_placebo.json`, `s23/results/d_scale_transfer.json`.

### THE VERDICT: TWO OF FOUR SURVIVE, TWO ARE ORDER STATISTICS — AND THE DISCRIMINATOR IS CLEAN

| oracle | claimed | Lane E verdict |
|---|---|---|
| **averaging width m\*** | −0.244 Å, 77/49 | **SURVIVES.** It was never a bare oracle — it is a **split-half transfer**, self-nulling for best-of-K. −0.2394, SE 0.0353, **2.42× MDE**, fold [−0.292,−0.180], **5/5 folds same sign**, 79W/38L/9T. Keeps **65.2%** of its own in-sample oracle. |
| **scale s\*** | −0.3403 Å, 126/0 | **SURVIVES, and it is the strongest of the four.** Held-out transfer −0.3198, SE 0.0485, **2.36× MDE**, fold [−0.392,−0.237], **5/5 folds**, 115W/11L, **worst-target degradation +0.0057 Å**. Keeps **98.9%** of the in-sample oracle — essentially no overfitting, the signature of a genuine per-target parameter. |
| **arm choice** | −0.482 / −0.505 Å | **ORDER STATISTIC.** N4 accounts for **112.6%** (K=19, k_eff 3.4); residual **+0.067 Å, wrong sign**. Robust to the arm-set fork: 112.6% / 115.1% / 117.8% at K=19/18/14. **No transfer arm exists for this claim.** |
| **cluster choice** | −0.484 Å | **ORDER STATISTIC, by exact arithmetic — no simulation needed.** The artefact records both `oracle` (min over clusters) and `random` (mean over clusters), so `oracle − incumbent = (oracle − random) + (random − incumbent)`. The first term IS the order statistic. The second is **+0.344 to +0.631 Å**: a randomly chosen cluster is **substantially WORSE than the incumbent** at every cell. The clustering operator is a worse operator that only looks good when you take the minimum over its own outputs. 174–245% accounted. |

> **The discriminator is not the size of the oracle, and it is not the W/L. It is whether the claim
> was ever validated OUT OF SAMPLE.** The two that survive are exactly the two that were built as
> split-half transfers. The two that fall are exactly the two that were only ever in-sample minima.

### THE NULL, AND WHY THE FIRST THREE VERSIONS OF IT WERE WRONG

I report this because an over-explaining null is a **mis-specified** null, not proof an effect is
fake, and three of mine were mis-specified before the fourth was right.

* **N1** exchangeable-column resample. Mis-specified whenever the baseline column is not the row
  mean — it centres the simulated grid on the *baseline*, deleting the `(rowmean − baseline)` term
  that is part of the observed gain. Safe on `conf.py`'s tight k grid, wrong here. Gave 105–168%.
* **N2** additive model, raw interaction permuted within column. Mis-specified by
  **heteroscedasticity** — hands a target whose grid barely moves a residual borrowed from a target
  whose grid moves by Ångströms, and the minimum is where that lands. Gave 134–296%.
* **N3** same, but standardised and rescaled to the receiving target's own sd. Fixes scale, still
  treats the K columns as **K independent opportunities**. Gave 109–135%.
* **N4, THE PRIMARY.** Reassign whole standardised residual **profiles** between targets. The
  profile's internal column-correlation — and therefore `k_eff` — survives; only the association
  between a target and **its own best column** is destroyed, which is exactly what "a per-target
  optimum exists" asserts. Gave **102.5 / 112.6 / 115.4 / 117.8%** across four independent panels.

`k_eff` (participation ratio of the residual correlation eigenvalues) is reported beside every
result and is the reason N1–N3 fail: **the scale sweep has k_eff 1.3 of 81 columns and the arm
panel k_eff 3.4 of 19.** A smooth 1-D sweep does not offer K independent chances at a low value.

N4 still over-shoots by 2–18%, so the correct reading is **"the null reproduces the whole observed
gain; no per-target signal is detectable above it"**, not "the effect is exactly zero".

### CONSEQUENCE FOR THE PROJECT'S THROUGH-LINE

"Real per-target headroom that no native-free rule can reach" is **true for m\* and s\***, and the
finite-sample router bound is still needed to explain those two. It was **over-applied** to arm
choice and cluster choice, where there is no headroom to explain. The corrected count is **three**
independent instances (m\*, s\*, and whatever else has a transfer arm), not five.

### TWO SMALLER CORRECTIONS FOUND IN PASSING

1. **A W/L transcription slip.** The cluster-choice claim is quoted as "−0.484 Å at m250_k3,
   109/17". Measured: `m250_k3` is −0.4839 at **84W/42L**; the **109W/17L** belongs to `m75_k5`,
   which is −0.4615. The two halves of that citation come from different cells.
2. **`d_scale_placebo` is read the wrong way round.** Its placebo delta is −1.3169 on a baseline of
   7.4987 (**−17.6%** of its own baseline) against the real −0.3403 on 3.0483 (**−11.2%**). In
   absolute Å the placebo is 5× the real effect and as a fraction of its own baseline it is
   **1.57×** — so that placebo does not show the real effect is "far outside", it shows the
   opposite. It does not matter, because `d_scale_transfer` settles s\* independently and
   decisively, but the placebo should not be cited as support.

---

## STANDING WATCHES

### Leakage (Task 4) — swept 2026-09-08 16:58, NOTHING FOUND

* **No trained model exists yet.** No `torch`, no `sklearn`, no `checkpoint`, no `state_dict`
  anywhere in `s24/`. There is as yet nothing for a native to leak into. Re-sweep when Lane B or C
  lands a model.
* **ORACLE labelling is being followed rigorously** across A, B, C and D — every native-consuming
  quantity is suffixed `ORACLE`, and `resid0b.py` even declares "the alpha reported is chosen on
  the same 126 targets it is scored on -> ORACLE, an upper bound". That is the discipline working.
* **Sealed benchmark.** Exactly two touches, both Lane A: `a_corpus.py:123` and `a_ha2.py:84`.
  Both read `p.pdb`, `p.seq`, `p.n` only — no benchmark structure, no `rr`, no result file. Within
  the identifiers-for-exclusion permission; sequences are strictly necessary for a homology audit
  and cannot carry an RMSD. **PASS.** Condition: the audit's per-benchmark-target leak flags are an
  exclusion list and must never select or tune anything.
* **Selector provenance.** `c_ladder.py`'s `gi` is a native-free random permutation. Clean.
* **Basis discipline.** Both Tier-1 files are point cloud on both sides throughout. `c_ladder.py`
  declares its point-cloud/built-chain split explicitly and the argument (an average of 75 built
  chains is not itself an ideal-geometry chain, so it is a point cloud) is sound. Noted but not
  flagged: its mixture ladder averages real windows together with built chains in one set, so the
  two endpoints differ in member provenance as well as in composition.
* **For Lane C.** `c_ladder.py:220` is the **third** file carrying the undeclared score-prefix
  composition fork (`A[:a]` and `B[:b]` both held in score order). My qmatch measurement says it is
  probably null, but it should be declared and ideally nulled. And `c_ladder.json` was
  `complete:false` at 90/126 when I looked — nothing from it should be quoted yet.

### Open

* Re-run the best-of-K null against the **other four** per-target oracles the ledger cites, if any
  of them minimise over a grid. Not yet done.
* Re-sweep for leakage the moment a trained model lands.
