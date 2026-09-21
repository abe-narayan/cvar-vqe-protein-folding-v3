# Sprint 32 — Where the RMSD is lost, and what it would take to get it back

**Status: FINAL — 2026-09-21 11:38.** All six lanes closed; no section is pending. Assembled as
results landed, so nothing here is reconstructed from memory at the end. Every number is
asserted against the artefact that produced it by `s32/s32_verify.py` at **143 recomputed,
0 mismatched, 0 flagged, 18/18 self-tests** across **eleven audits** — each built from the defect
that motivated it, and each carrying a self-test that must trip on that defect.

Branch `s26` · instrument `tuning126`, 126 targets, 9–16 aa · endpoint **mean built-chain Cα RMSD**
Charter `s32/BRIEF.md` (verbatim, 2,087 lines) · Ledger `s32/LEDGER.md` · State `s32/STATE.md`
Contract `s32/S32_CONTRACT.md` · Causal map `s32/CAUSAL_MAP.md` · Verifier `s32/s32_verify.py`
Multiplicity `s32/MULTIPLICITY.md` · Theory `s32/THEORY_Q.md`
Published: <https://claude.ai/artifact/BiAF9fiysaHGTHH1CYNQmL> · rebuild with `python s32/build_report_page.py`
(the page renders this file client-side, so the artifact and the repository cannot drift)

---

## 0. The answer, up front

**The endpoint did not move. It is 3.2105 Å.** Walking every `effect_over_mde` in all 56 result
artefacts:

> ### Improvements clearing 1.0× MDE on the BUILT CHAIN, deployable: **ZERO.**
> Thirty arms clear 1.0× on *some* basis — about twenty read the native, four are filter-vs-random on
> the pool mean, six are in-band diagnostics. **Every arm that clears MDE either reads the native, is
> measured on a basis that is not the endpoint, or points the wrong way.** That is exhaustive over the
> artefacts rather than over memory.

The charter's primary target of < 3.00 Å and its ambition of < 2.50 Å are both unmet.

**What the sprint produced instead is a causal account that makes the failure necessary rather than
unlucky, and it rests on four results that were not available before.**

### 0.1 The readout is a hull projection, and that closes a family by derivation

`‖Σ w_x W_x − t‖² = ⟨w,a⟩ − ½ w'Bw` **is** the squared distance from the native to a point of the
candidate hull. So the readout's optimisation is the **Euclidean projection of the native onto that
hull**, `∂x*/∂t` is the projector onto the active candidates' affine hull — **gain exactly 1 inside,
0 outside** — and the binding constraint is the **hull floor `d = 1.8290`**.

> **A structure estimate good enough to make the readout worth solving is already good enough to
> emit.** S31 measured that solving the CVaR objective exactly *"reshuffles the answer on 66 of 126
> targets and buys nothing"* and recorded it as a surprising null. **It is forced.**

### 0.2 `a` and `μ` are one object, and the requirement is 3n − 6 real numbers

Per-candidate quality `a` and the pool's common mode `μ` are related by a **native-free affine
bijection** — established independently by two lanes working in opposite directions, one deriving
`a` from `μ`, the other recovering `μ` from `{a_k, d_k}` at relative residual **2e-14 on 126/126**.
**S31 called this the sharpest question it produced and left it open. It is closed.**

> **A per-candidate quality estimator with in-band skill IS a structure predictor.** And since
> `rank(d) = 3n − 6` **exactly on 126/126**, the requirement is **exactly ~39 real numbers with no
> compression available.**

### 0.3 The 2.10 Å of "headroom in the pool" is a FIT, not a RETRIEVAL

The control that had never been run: **500 fragments drawn from a *different protein's* universe
reach 1.1626 Å on this target, against the retrieved pool's 1.1167 Å — NOT A RESULT (0.49× MDE).**

> **Do not call a pool "containing the answer" when a pool assembled for a different protein contains
> it equally well.** At 9–16 residues `3n ≈ 39`, and a hull of 500 fragments is close to saturated.
> **The ladder's top rung prices a fit.** Retrieval's share of the 2.10 Å is **2.2%**; what retrieval
> actually buys is the **pool mean**, which is what a uniform average consumes.

### 0.4 "In-band skill is zero" has been read wrong for four sprints — and it still does not help

`Var(ρ) > 0` with `E[ρ] = 0`: the ordering information is **present on every scorer** — AMBER at
exactly +0.0000 mean skill still carries |ρ| at **1.88×** a within-band permutation null — and **the
per-target SIGN is what is missing.** The sign is a property of the **target**, transferring across a
split half at 2.6–4.2× MDE under a *cross-target* null, with a pure-noise arm at 0.15× proving the
test can fail.

**And the ceiling closes it anyway.** 2.0 Å needs in-band `ρ ≈ 0.638`; a **free, perfect** sign gives
the best scorer **0.3228**.

> **Even a free, perfect per-target sign leaves the best in-band scorer 2–3× short. The sign is a real
> finding and a closed route.** Worse for the architecture: **positive in-band skill does exist on the
> deployable band — it is typicality — and the terminal operator is already its argmin.** *The only
> native-free channel with in-band skill is redundant with the operator already deployed.*

### 0.5 The quantum question, answered

**Charter §14: no, on this instrument, and the binding reason is chain length.** `n_res` is 9–16, so
`2^n_res ≤ 65536` on 126/126 and the enumerability condition fails outright. The two properties a
problem would need — **a decision space that grows with the target**, and **a genuinely stochastic
energy** — must hold together, **and this project has never had either.** The last formulation
standing, sparse `s`-of-`K`, escapes all three of S31's obstructions and then falls to monotonicity:
***the hard instances are exactly the ones whose optimum is worse.***

### 0.6 The pattern in this sprint's own errors

**Eleven defects were found in work already written down, and the coordinator wrote most of them.**
Three of mine: a ladder increment quoted on the **member** basis inside a **chain** ladder; a
convexity check that **initialised its accumulator at the pass threshold** and could not fail, written
one hour after I wrote the rule forbidding exactly that; and an inverted gloss on a column named
after the property it is *inversely* related to.

**Every one is the same shape as S31's: a number re-used across a boundary its definition does not
cross.** What is different this sprint is that **most were caught by machine rather than by memory** —
the verifier grew from 58 checks to 134 with 18 self-tests, and each new audit was built from the
defect that motivated it, including one that flags *any aggregate whose FAIL18 and non-FAIL18 strata
have opposite signs*. That generalisation caught its third instance in ten minutes instead of a
retraction.

**And a whole class is now named**: *any comparison of the score's prefix against an alternative
prefix, scored on the BEST MEMBER, is entirely FAIL18* — because FAIL18 is **defined** as the targets
whose prefix excludes the good band. Three instances, each with the aggregate clearing MDE and the
non-circular 108 NOT A RESULT.

---

## 1. What the endpoint is, exactly

The charter's Step 4 asked for the 3.2105 Å endpoint to be verified from artefacts rather than taken
from a report. It was, and the answer is more precise than the number.

**3.2105 Å** is *the λ = 0.3 multi-start projection arm of `s12/instrument.project` applied to
`s29/results/s29_O_structs/<pdb>.npz["prod"]`, CA-RMSD to native, meaned over `tuning126`.*
**It is not a cached scalar.** Five distinct objects live near it and they are **not five estimates
of one**:

| value | what it is |
|---|---|
| 3.048338 | CA point cloud — the top-75 coordinate average, unprojected (`rmsd_avg`) |
| 3.204076 | the **λ = 0** arm — nearest ideal geometry, no Ramachandran penalty (`rmsd_fit`) |
| 3.214765 | the **λ = 0.3** arm — **the chain production itself emits** (`rmsd_arm` / `ca`) |
| 3.235460 | that chain **after AMBER relaxation** (`rmsd_full`) — AMBER costs **+0.0207 Å** |
| **3.210534** | **the canonical endpoint** — a *re-projection* of the stored cloud, λ = 0.3 |

**The canonical endpoint is 0.0043 Å better than the chain production emits.** A reader who assumes
the endpoint is the pipeline's output is wrong by more than several historical claims are large.

### 1.1 Everything upstream of the chain reproduces bit-exactly

Targets 126; folds 25/23/25/23/30; the top-75 **set** reproduced 126/126 from an independent
rescoring; the recomputed coordinate average matches the stored one to **1.42e-14**; cloud
**3.048338**; set mean **3.550683**; pool best K=500 **1.710824**; top-75 best **2.306153**; shipped
argmin **3.454000** to six decimals.

### 1.2 The projection operator is bit-reproducible; its *input* is not uniquely determined

| input cloud | mean chain | mean `d` | max \|d\| | bit-identical |
|---|---|---|---|---|
| `s29_O_structs/<pdb>.npz['prod']` | 3.210533995 | +0.000000 | 0.000000 | **126/126** |
| production cache `avg_ca` | 3.214765154 | +0.004231 | 0.416765 | 0/126 |
| recomputed coordinate average | 3.212625220 | +0.002091 | 0.517410 | 0/126 |

All three are **the same top-75 coordinate average, agreeing to 5.7e-14**, with cloud RMSDs identical
to 9 dp. `s12.instrument.project` reproduces **bit-for-bit across processes, across BLAS thread
counts 1/2/4/8, and across the gap between S29 and today**. Lane R reproduced five separate ladder
rungs bit-for-bit from a different job, script and process.

> **The operator is bit-reproducible. The input is what is fragile.** A **one-ULP (7.1e-15 Å)**
> change in the cloud moves the built chain by **0.10–0.15 Å** — deterministic given identical bits,
> **discontinuous** in them.

**Consequence:** the anchor carries ~**±0.003 Å** of pure arithmetic noise — the draw-to-draw sd of the
endpoint over five 1e-14 Å perturbations (3.207688, 3.207625, 3.206217, 3.203375, 3.211812); the
three independent honest recomputations of the same average span **0.0042**. The charter's targets stay
checkable (< 3.00 is 0.21 Å away, < 2.50 is 0.71 Å). **Unpaired cross-job chain claims below ~0.03 Å
are not resolvable**, and "the same cloud value" does not license a comparison — it must be the same
float64 bits.

---

## 2. Where the RMSD is actually lost

**Built chain, all ORACLE rungs ORACLE / NOT DEPLOYABLE:**

```
best sparse convex combination, K=500, s=10      1.1139        2.10 A of headroom
best single member, K=500                        1.7078        1.50 A of headroom
  + distogram SCORE prefix, 500 -> 128          +0.4357   ->   2.1435
  + prefix 128 -> 75                            +0.1620   ->   2.3055
  + selection / readout (uniform average)       +0.9051   ->   3.2105   PRODUCTION
```

**The pool is not the bottleneck** — but the reason is not the one this ladder first suggested, and
the correction is the deepest result of the sprint.

### 2.0 The top rung prices a FIT, not a RETRIEVAL

The obvious reading — *"the K=500 pool supports 1.1139 Å, so everything downstream destroys 2.10 Å
that retrieval already found"* — attributes the headroom to retrieval. **A size-matched control that
had never been run says retrieval's share of it is 2.2%.** Three arms, same medoid frame, same NNLS
solver, same alternation rounds, same native-pose seeding, n = 126, CA point cloud, **ORACLE / NOT
DEPLOYABLE**:

```
A  BLOSUM500   the shipped pool                        1.1167   (median 1.1002, support 9.9)
B  RAND500     same universe, RETRIEVAL-blind          1.1495   +0.0328   0.38x   63W/63L  NOT A RESULT
C  DONOR500    ANOTHER TARGET'S universe, TARGET-blind 1.1626   +0.0459   0.49x   47W/79L  NOT A RESULT
```

Arm A reproduces the S29 ladder's `hull_pool` = 1.1167 to four decimals, which is what certifies the
implementation is the ladder's operator rather than a lookalike. (Support **9.9** also explains why
the ladder's sparse `s=10` ≈ the full hull: **the hull is already sparse.**)

> ### 500 fragments taken from a DIFFERENT protein reach 1.1626 Å on this target, against the retrieved pool's 1.1167 Å, and the difference is NOT A RESULT.

**So the 2.10 Å is not information retrieval supplied. It is the expressive capacity of 500
fragments in a ~39-dimensional space** — a convex hull of 500 points in `R^{3n}` with `3n ≈ 39` comes
within ~1.16 Å of essentially any target of that length, whether or not the points were chosen for
it. **And realising it requires the ~39 real numbers that ARE the answer.**

Charter §27 warns: *do not call a bad candidate pool a selection problem just because the best
candidate exists somewhere inside it.* **This is that warning one level up: do not call a pool
"containing the answer" when a pool assembled for a different protein contains it equally well.**

**The corrected statement, which is stronger than the one it replaces:**

> **At 9–16 residues fragment space is close to saturated, so the ladder's top rung prices a FIT, not
> a RETRIEVAL.**

**What this does NOT say**, and the distinction is load-bearing: **retrieval is not useless.** It sets
the **pool mean** — 4.4533 against a random universe draw's 4.8155 — and production **averages 75
members rather than fitting a hull**, so mean quality is what production actually consumes.
***The ORACLE top rung and the deployed pipeline consume different properties of the pool, and only
the second is retrieval-sensitive.***

It also makes three other results cohere rather than sitting beside them: lane P's best-member rung
`universe → 500` is **−0.0718 at 0.85×, NOT MEASURED** — retrieval's contribution to the *best member*
is unmeasured and of the same size; the readout's hull floor `d = 1.8290` is for the **top-75** hull,
and the K=500 hull is 1.1167 **because it has 500 points, not because they are the right ones**; and
charter §41's five-bit question is answered by the same fact — **the 2.10 Å "prize" costs ~39 real
numbers of oracle information to collect, which is why every attempt to collect a few bits of it has
returned nothing.**

Charter §56 asks for the earliest irreversible loss. **It is not candidate generation** — and now for
a better reason than "the pool already contains the answer".

### 2.1 Most of the "filter loss" is a bare order statistic, and the rest is 18 circular targets

The 500 → 128 step is **not** retrieval — it is the **distogram Bayes-risk score prefix**, which the
S29 code calls *"the quantum field of view"*: **128 = 2⁷, the VQE register width.** Retrieval is the
earlier arrow. Priced against a **size-matched random subset**, 2000 draws per target (CA cloud,
ORACLE):

```
                              TOTAL      set-size order statistic     ordering effect
K=500 -> 128                 +0.4350            +0.2477                   +0.1872
128  ->  75                  +0.1604            +0.0907                   +0.0696
```

**57% of the 500 → 128 loss is bare set size** — a minimum over 500 is lower than a minimum over 128
for *any* subset. The ordering effect is formally past MDE (1.06×, fold CI excluding zero, 4/5 folds),
and **it is entirely 18 targets**:

```
                 ALL 126              FAIL18 (n=18)        OTHER 108
500 -> 128   +0.1872 (med -0.0380)   +1.4879  0W/18L    -0.0296   0.30x   NOT A RESULT
128 ->  75   +0.0696 (med -0.0045)   +0.3611  3W/15L    +0.0210   0.38x   NOT A RESULT
```

The ten worst targets are **all ten in FAIL18**; the score loses **18 of 18** there. Drop the 10 worst
and the aggregate falls to +0.0294; drop 20 and it goes **negative**. And **FAIL18 is defined in
`s12/instrument.py::selfcheck` as the targets where no pool member within 1.5 A of the pool best
survives into the top-75** — so a contrast asking *"does the score's prefix retain the good
members?"* is **near-circular on precisely those 18**, and *directly* circular for the 128 → 75 arm.
The filter-independent control (split by chain length, median 13) shows **no gradient**: +0.1826 short
against +0.1940 long. **It is not a broad property. It is the 18.**

> ### The score is not a general anti-ordering. On 108 of 126 targets its top-128 retains a marginally *better* best-member than a random 128 — NOT A RESULT. Its failure is catastrophic and total on 14% of targets, and those are exactly the targets already named FAIL18. **The problem is not that it orders badly everywhere; it is that on one target in seven it places its window in the wrong part of the pool entirely.**

This was **independently reproduced by two lanes from different raw artefacts** — the first time in
this project that a retraction has been confirmed that way.

**And the earlier arrow is fine.** `universe → 500` is **−0.0718**, with FAIL18 contributing **0%**
and the other 108 at −0.0833 (0.98×, NOT MEASURED but the right sign). **Retrieval is not the
problem.**

**Also settled and not re-opened:** the prefix length `m` does not transfer — ORACLE global `m*` = 72
is −0.0044 (0.19×) and the leave-fold-out `m` is **+0.0075, 63W/63L, a literal coin flip** on the
endpoint; and a matched **random-subset** family reaches **141%** of the prefix family's gain on the
chain, because prefix variants are nested (lag-1 autocorrelation 0.917 against 0.112) so a
less-correlated family has a larger per-target minimum.

**The corrected increment sentence:** *narrowing 500 → 128 → 75 costs **0.598 A on the BUILT CHAIN** (0.595 on the member/cloud basis —
name the basis) of oracle-best headroom; 0.338 is the set-size order statistic; the remaining 0.257 is 18 targets' worth of the
score placing its window wrongly, and is NOT MEASURED on the other 108.*

**Rank moments, both true and different:** the K=500 best member sits at mean rank **170.3**, median
**134.0** of 500; it survives into the top-128 on **63/126** targets and into the top-75 on
**41/126**. *On half the targets the best available candidate is gone before the readout ever sees
it.*

### 2.3 The ranking cell, measured on the BUILT CHAIN

The S29 ladder's *best single member* changes the ranker **and** the readout at once and reports a
structure production never emits. The cell nobody had isolates the **ranker alone**, leaving the
shipped uniform average in place. Lane P computed it on the cloud and left the chain job queued; it
completed after the lane closed and is **measured, not interpolated** — all three arms projected in
the same process per target, n = 126, **ORACLE / NOT DEPLOYABLE**:

```
                                                   chain    vs PROD    xMDE    fold CI
PRODUCTION (in job)                               3.2126
ORACLE top-75 through the DEPLOYED average        2.1186   -1.0941    4.17x   [-1.195, -0.995]
ORACLE top-5  through the DEPLOYED average        1.5703   -1.6424    5.31x   [-1.778, -1.509]
                                          other 108:       -0.8283
```

> ### **−1.09 Å from ranking alone, with the readout untouched** — and **−1.64 Å** at the ORACLE
> global `m = 5`. *The prefix curve is non-monotone: averaging the five best beats naming the single
> best by 0.23 Å.* **Error cancellation, a second independent instance.**

**This is the largest single ORACLE lever in the sprint, and §3 is why it cannot be collected:**
ranking and readout are **the same requirement**, and knowing which 75 to rank means knowing `a` on
the active subspace, which by the affine bijection means knowing the structure.

### 2.2 The decisive deployable test: replacing the score's window with a random one, at the endpoint

The one intervention in this sprint that needs **no native information at all** — and therefore the
cheapest possible deployable gain — is to throw the score's prefix away. Lane P ran it end to end:
**random 128 of the 500, then the score's top-75 within it, BUILT CHAIN, every arm projected in the
same process as `PROD` for that target, 8 draws, n = 126.**

```
                      effect    xMDE    W/L      draws better    verdict
ALL 126              +0.1648    0.92    52/74      0 of 8        NOT MEASURED, and worse in sign
FAIL18 (circular)    -0.4149    0.70    15/3                     NOT MEASURED -- random HELPS here
OTHER 108            +0.2614    1.52    37/71                    WORSE
```

**The aggregate is two opposite effects cancelling, and the stratification inverts.** Randomising the
window **helps on the 18 targets where the score misplaces it** (15 of 18 win) and **hurts on the
other 108** (1.52x MDE, WORSE). The draw-to-draw sd is **0.0187** against an effect of 0.1648 and
**0 of 8 draws beat production**, so the *direction* is not in doubt even though the paired
per-target MDE says NOT MEASURED.

> ### The score's value is real and it is MEAN CANDIDATE QUALITY FOR AN AVERAGE. Its failure is WINDOW PLACEMENT on one target in seven. Randomising discards the first to fix the second, and the first is worth more.

**This is charter §29's trap, confirmed at the endpoint rather than argued.** A filter can retain a
*worse best member* than random and still produce a *better prediction*, because production averages
75 candidates and never takes a best member. **Every conclusion in §2.1 is about the best member and
none of it transfers to the endpoint** — which is why it was run.

**My registered prediction, committed before the arm reported, was "unchanged or worse".** It holds.
The prediction was made *less* safe mid-flight by lane V's discovery that the score also halves the
pool's spread — a term S32-L2 shows the readout *rewards* — and it survived that anyway: the
mean-quality gain (pool mean 4.4533 → 3.5847) outweighs the spread loss.

**What it leaves open, and it is the only live route in this arrow:** keep the score, and **detect
the 18 targets whose window is misplaced.** That detector must be native-free, and it is the same
missing quantity as §4's per-target sign — not an independent opportunity.

---

## 3. The readout is a hull projection, and that closes the family by derivation

For any weights with `Σw = 1`, `a_x = ‖W_x − t‖²` (**ORACLE**) and `B_xy = ‖W_x − W_y‖²` (native-free):

```
|| sum_x w_x W_x - t ||^2  =  <w, a>  -  0.5 * w' B w
```

`B` is a Euclidean squared-distance matrix, so `w'Bw` is **concave** on the simplex and the program is
**convex and tuning-free**. Reading the terms: minimising wants **low `⟨w,a⟩`** (good candidates) **and
high `w'Bw`** (spread — the variance-cancellation term). ***Spread is rewarded.*** That is why
quality-blind dispersion maximisation — this program with `a` constant — is **+0.1436 Å at 1.22× MDE,
WORSE**: it picks garbage. **`a` is load-bearing and it is the only unknown.**

### 3.1 The left-hand side is a distance to the hull, and that fixes everything

`‖Σw W_x − t‖²` **is** the squared distance from the native to a point of the candidate hull, so the
readout's optimisation is the **Euclidean projection of the native onto the convex hull of the
candidates**. Therefore `∂x*/∂t` is the orthogonal projector onto the **active** candidates' affine
hull: **gain 1 inside it, 0 outside** (finite differences 1.0000000, sd 3.1e-09; orthogonal 4.0e-09;
re-derived independently by the coordinator at 0.9991 / 5.9e-03).

- **`a` is needed only along `|S| − 1` ≈ 5.25 directions**, out of an ambient ~38.9. Everything else in
  `a` is **exactly invisible** to the emitted structure.
- **Gain 0 on ~33 directions is TOTAL suppression, not none** — a generic error is ~87% annihilated,
  and at ε = 4.0 the program turns a 3.65 Å estimate into a 2.23 Å emission.
- **The binding constraint is the LOWER bound**: `d ≤ ‖P_C(t̂) − t‖ ≤ d + ε` with the **hull floor
  `d = 1.8290 ± 0.1178`** (CA cloud, ORACLE). The crossover against direct emission is at **ε ≈ 2.2**.

> ### A structure estimate good enough to make the readout worth solving is already good enough to emit.

**That turns S31's measurement into a theorem.** S31 found that solving the CVaR objective exactly
*"reshuffles the answer on 66 of 126 targets and buys nothing"* (−0.0112 Å, 0.19× MDE) and recorded it
as a surprising null. **It is forced.**

### 3.2 `a` and `μ` are ONE object — S31's sharpest open question, closed twice

`a` is **affine in `t`** and `Σw = 1` makes `‖t‖²` an additive constant, so the readout-relevant part of
`a` is `P_aff{W} t`. Substituting it moves `w` by **7.0e-12** and the structure by **2.9e-12 Å** on
126/126. Since `t = X̄ − μ`, the map **`μ ↔ t ↔ a` is a native-free affine bijection.**

> **A per-candidate quality estimator with in-band skill IS a structure predictor, and a common-mode
> corrector IS a per-candidate quality estimator. One missing channel in two vocabularies.**

**Two lanes established this independently and in opposite directions** — one derived `a` from `μ`,
the other **recovered `μ` from `{a_k, d_k}` by least squares at relative residual 2e-14 on 126/126**
(shuffled-`d` control 0.373, `cos = 1.0000`). *Contract rule 7 satisfied without coordination.*

**And there is no compression hiding in the pool's geometry.** Superposing members *and* native on the
medoid removes 3 translations and 3 rotations, so `rank(d) = 3n − 6` **EXACTLY on 126/126** (median
33) and 128 fragments span it fully. **`μ` lies inside that span necessarily, and the requirement is
exactly `3n − 6` numbers.**

### 3.3 The last quantum formulation falls by monotonicity

Sparse `s`-of-`K` **escapes all three of S31's obstructions** — the subset basis is diagonal, `⟨E⟩` is
linear in `p` because the weights are solved classically inside `E(x)`, and the dimension is `2^K`
not `K`. ***Those obstructions are properties of the candidate-index register, not of CVaR-VQE.***

It falls to a different argument: **`f*(s)` is constant for `s ≥ s*` and strictly worse for
`s < s*`.** So `s ≥ s*` **is** the convex program and `s < s*` is hard **and worse**.

> ### The hard instances are exactly the ones whose optimum is worse.

`s*` at K = 500: mean **10.06**, median 10, p90 13, max 23; **61.1% ≤ 10.** The closure does not depend
on where `s*` falls.

### 3.4 Priced in reals rather than bits

```
r reals     0        1        3        6       10      ~33
CA cloud  3.0532   2.5804   2.2638   2.0312   1.9110  1.8290      ORACLE / NOT DEPLOYABLE
```

**Six reals buy 83% of what thirty-three buy.** *A bit count prices a selection alphabet; this
decision is continuous, so the two are never differenced.* That is charter §41's five-bit question
answered in the right currency.

---

## 4. In-band skill is not zero — the per-target SIGN is missing, and it is still not enough

Four sprints have read *"in-band skill is zero"* as **the information is absent**. It is not.

```
                in band      xMDE   W/L        mean|rho|   perm null   ratio   xMDE   folds
AMBER           +0.0000      0.00   62/64        0.1779      0.0948    1.88    2.36    5/5
DIS             +0.0652      0.83                0.2496      0.0939    2.66    3.05    5/5
LEG_total       +0.0376      0.42                0.2819      0.0908    3.10    3.52    5/5
LEG_torsion     +0.0444      0.65                0.2142      0.0947    2.26    2.77    5/5
```

> ### `Var(ρ) > 0` with `E[ρ] = 0`. The ordering information is present on every scorer — including the one with exactly zero mean skill — and what is missing is the per-target SIGN.

### 4.1 The sign is a property of the target, and the right null proves it

Estimated on half A of a band, applied to **held-out half B**, 16 splits/target, **deduplicated**
(7.7% of band members are exact coordinate duplicates on 122/126 targets), from a committed script
with a pinned seed:

```
scorer        transfer  nullPERM  nullXTGT  globalSGN   xMDE  folds  verdict
AMBER          +0.1163   +0.0058   +0.0004    -0.0725   2.61   5/5   PER-TARGET
DIS            +0.1901   +0.0010   +0.0088    +0.0649   3.19   5/5   PER-TARGET
LEG_total      +0.2180   +0.0042   +0.0082    +0.0426   3.10   5/5   PER-TARGET
LEG_torsion    +0.1492   +0.0046   +0.0017    +0.0460   2.68   5/5   PER-TARGET
RG             +0.3228   +0.0014   +0.0053    +0.0489   4.15   5/5   PER-TARGET
NOISE          +0.0056   -0.0004   +0.0016    +0.0123   0.15   2/5   NOT A RESULT  <- self-test
```

**`nullPERM` — the first null used — could not test the claim.** It permutes `rr` inside half B and
destroys *all* structure, so it is ≈0 for every scorer **including pure noise**. The null the claim
needs is **cross-target**: apply another target's sign. Under it the transfers are genuinely
per-target — and **`NOISE` at 0.15×, 2/5 folds is the self-test proving the audit can fail.**
`DIS` replicates across two independently coded implementations at **+0.1901 / +0.1911**.

**Two framing corrections the lane imposed on its own result.** These are **not** unrelated
Hamiltonians — they all load on compactness, and **plain `Rg`, one line of numpy, beats every
Hamiltonian measured.** *The bit is very likely "is the native more or less compact than its own
band."* And this is **largely a confirmation on a new instrument, not a discovery**: the project
already recorded native-free compactness proxies at 0.24–0.37, and `Rg`'s 0.3228 lands inside that.

### 4.2 The ceiling, which closes the route

2.0 Å requires in-band `ρ ≈ 0.638`. A **perfect, free** per-target sign gives AMBER 0.1163 (18% of
it), `DIS` 0.1901 (30%), `LEG_total` 0.2180 (34%), **`Rg` 0.3228 (51%)**.

> ### Even a free, perfect per-target sign leaves the best in-band scorer 2–3× short of the useful range. **The sign is CLOSED AS A ROUTE and kept as a FINDING.**

**And it is not supplied natively.** Leave-fold-out ridge over 20 native-free per-target features,
scored **against the marginal** rather than 50%: AMBER 0.381 vs marginal 0.508 (**worse than a
constant**); `DIS` 0.548 vs 0.548 (**exactly**); `LEG_total` 0.595 vs 0.587 (**one target in 126**);
per-fold accuracy swings 0.30–0.70. The one mechanism-motivated single feature,
`sign(rg_pred − rg_pool)`, fails with a named mechanism: **it is positive on 81% of targets, so it is
a constant wearing a label. A one-bit feature that is 81/19 cannot carry a 59/41 label.**

**In Ångströms, BUILT CHAIN, n = 126, all arms in one process:** the price of **one ORACLE bit** is
**−0.2101 Å, 1.66× MDE, 5/5 folds — ORACLE / NOT DEPLOYABLE** (`s32_D5_signchain_LEG_total.json`;
*the cloud-basis figure in `s32_D5_signprice.json` is −0.1287 at 1.52× and the two are never
differenced*). Against S31's −0.8102 Å for the
common-mode *direction*, **one bit is 26% of it.** (Median exactly 0.0000 with 74 ties, and that is
**structural**: where the bit is +1, ORACLE ≡ constant. The effect is carried by 52 targets at
≈ −0.51 Å each.)

### 4.3 And the one native-free channel that DOES have in-band skill is redundant

S31's *"in-band skill is zero or the wrong sign"* is a property of the **ORACLE band**, which
conditions on the quantity being predicted. On the band a deployed selector actually ranges over —
the score's own top-24, no native conditioning — consensus has **ρ = +0.2531 ± 0.0434** on the 128
and **+0.2829 ± 0.0322** on the production 75, **positive on ~73% of targets.**

> ### Positive in-band skill EXISTS. It is typicality — small `|d_k|`. And the terminal operator is already the argmin of exactly that: the medoid-superposed average sits at `V = 0`. **The only native-free channel with in-band skill is redundant with the operator already deployed.** That is the wall, and it is mechanical rather than statistical.

---

## 5. The projection is not the earliest irreversible loss — it faithfully transmits an upstream defect

The reconstruction costs **+0.1622 A, 5.05% of the endpoint**, and the charter (§31) asks whether it
is destroying structural improvements. **It is not.**

### 5.1 The cheap projection price was ORACLE-induced

The obvious reading of the S29 ladder — *sparse combinations project for free, dense averages do
not* — is **wrong**, and the control that shows it is the sharpest in the sprint. Same sparsity,
same averaging operator, same projection, same job, three pinned draws (n = 126):

```
arm                          cloud   chain      d  |   price   orthog null |   cos   | vs null
ORACLE sparse s=10          1.1136  1.1139  0.7023 |  +0.0002      +0.2124 |  +0.380 | 4.65x BETTER
RANDSPARSE s=10, 3 draws    3.5541  3.7529  1.1701 |  +0.1988      +0.2088 |  +0.019 | 0.22x NOT MEASURED
SCORESPARSE s=10 (top-10)   3.1455  3.2826  0.6274 |  +0.1371      +0.0905 |  -0.076 | 1.46x WORSE
PROD s=75                   3.0483  3.2105  0.8150 |  +0.1622      +0.1504 |  -0.061 | 0.24x NOT MEASURED
```

> ### The same s = 10 combination pays +0.0002 when the native chose its members and +0.1988 when it did not. **Every object the native did not touch sits at or above its own orthogonal null**, and the deployable score-top-10 is on the *wrong* side of it at 1.46x MDE.

*The `price` and `orthog null` columns are **per-target means, averaged independently**, so they do
not close against each other by arithmetic — for PROD, `sqrt(3.0483² + 0.8150²) − 3.0483 = 0.1071`
against the table's 0.1504. **That gap is Jensen, not an error.**

**So neither sparsity nor alignment is an available intervention.** They are diagnostics of an
object that was already good.

### 5.2 The price is set upstream, by a native-free quantity

Production's off-manifold distance `d` is **rank-determined by the pool's own disagreement**:
`spearman(d, mean pairwise RMSD of the 75 members)` = **+0.9646**, partial on chain length **+0.967**,
partial on the native error **+0.953**, per fold 0.935 / 0.971 / 0.982 / 0.961 / 0.967. **Both sides
are native-free.** The shared-referent floor was measured *first* — permuting the spread within
chain-length strata gives ρ ≈ **+0.10**, max **+0.466** over 4000 draws. Honest limit: the *ratio* has
cv 0.467, so the relation is **monotone, not a proportionality** — quote ρ, never a coefficient.

The whole price therefore decomposes with **only the last link reading the native**:

> **pool disagreement → (ρ +0.965) → `d` → (cos −0.061, orthogonal to slightly adverse) → price
> = +0.1622 = 5.05% of the endpoint.**

**Charter §56, answered for this stage: the projection is bit-reproducible, it reproduces its own
historical numbers exactly, and it adds error in quadrature at a rate fixed upstream. It transmits a
retrieval defect; it does not create one.**

### 5.3 Scalar dilation is closed in both calibrations, and it reconciles two long-quoted numbers

> **✓ PROVENANCE CLEARED.** `s32/s32_R_dilation.py` is committed at `63608f9d` with
> `ST.provenance(__file__)` and the grid pinned in source; **both artefacts were re-emitted
> and every number reproduced exactly.** *It should not have rested on a heredoc — and the
> audit that caught it, asking whether any `.py` actually WRITES the filename, catches a
> failure mode that prose-level checks miss entirely.*

Registered P1.3, **falsified by its own falsifier**: dilating the cloud to ideal virtual-bond length
is **WORSE by +1.0425 at 2.79x MDE, 5/5 folds**; the better-motivated Rg-matched dilation is also
worse, **+0.0622 at 1.47x**. The reason is that the contraction is **separation-dependent**:

> **22.15% at |i−j| = 1, and only 5.40% in the radius of gyration.**

*One scalar matched to one moment is wrong at the others.* **That reconciles two figures this project
has been quoting past each other** — the "3.5% contraction" in project memory and
`core/project.py`'s "2.96 against 3.80" are **both right and measure different separations. Neither
may be substituted for the other.**

### 5.5 The branch-selection hypothesis was mine, and it is falsified on the endpoint

**I opened lane R on a specific mechanism**: `core/project.py` documents two ideal-geometry torsion
solutions *"one Ramachandran-plausible and one not"*, with the reference disagreeing with itself by
up to 1.6 A; every native-free ranker this project has tested is a distance-map function and
therefore **achiral** by theorem G1; **so a chiral criterion should be able to pick the branch where
an achiral one provably cannot.** Lane V tested it as a deployable argmin arm, chain basis, paired
in the same job, tie-averaged, independently of lane R (n = 126; **production in that job is exactly
3.2105**, so these arms sit on the canonical realisation):

**The decision is real.** 126/126 targets have ≥ 2 distinct branches (median 149 of 158, union-find
at 1e-3 Å); the RMSD spread across distinct branches is mean **0.5542**, **> 0.3 Å on 86 targets**;
and the λ = 0.3 objective gap between the best two is **below 1e-3 on 102/126**. *Structurally
distinct branches, chosen by a near-degenerate objective.* **A self-test that could have failed:**
production's own selection rule, re-executed on the recorded branch scalars, reproduces `prod_chain`
on **126/126 with max |error| exactly 0.0** — the branch set provably contains production's answer.

```
                                      chain     vs PROD    xMDE     W/L        ORACLE / NOT DEPLOYABLE
ORACLE best over GEN4      (K=4)     3.1751    -0.0355    1.73x   103W/8L
ORACLE best over GEN4+GEN4D (K=8)    3.1168    -0.0938    2.92x   114W/0L    <- replicates S31-D exactly
ORACLE best over ALL       (K=158)   3.0941    -0.1165    3.43x   124W/0L
  SPLIT-HALF TRANSFER                          -0.0046   =  4% of the oracle  ->  NOT A SIGNAL
  RAND0 zero-signal control (arbitrary column index by construction)    3%
random branch, 5 draws               3.2106    (draw sd 0.0028)
best of 29 deployable arms           3.2006    -0.0100    0.61x              NOT MEASURED
PRODUCTION (lane R, in job)          3.2105    -- exactly the canonical value
```

**The K = 8 row replicates S31-D from a different job and script**: S31 recorded 3.1167631446217534,
−0.09377085032154638, 2.924×, 114W/0L; this run gives 3.1168, −0.0938, 2.92×, 114W/0L. *A cross-sprint
number reproduced rather than quoted.*

> ### Production's multi-start argmin is worth **0.0001 Å** over a coin. Any branch rule competes for 0.117 Å of ORACLE headroom — **96% of it an order statistic** — against an incumbent indistinguishable from random.

> ### Production's multi-start argmin is worth 0.005 A over picking a branch with a coin. Any branch-selection rule is competing for 0.11 A of oracle headroom against an incumbent that is 0.005 A better than random — and **95% of that 0.11 A is an order statistic that does not survive a split half.**

And the sixteen native-free criteria, tested in **BOTH DIRECTIONS** on the full branch set and on
production's own four starts — **64 comparisons, n = 126, paired in-job.** *Testing one direction only
would have been an unregistered choice that halves the apparent multiplicity, and it is what produced
the inverted gloss an earlier draft of this section carried.* **Lane R's in-job production is exactly
3.2105** — it projects the canonical input — so these arms sit on the canonical realisation.

```
                 ARGMIN (best-looking)             ARGMAX (worst-looking)
d_to_C        -0.0091  0.57x  not a result     +0.3324  2.59x  RESULT (worse)  34W/92L
obj1          -0.0068  0.35x  not a result     +0.3408  2.69x  RESULT (worse)  30W/96L
rg            -0.0011  0.04x  not a result     +0.3261  2.76x  RESULT (worse)  30W/95L
typicality    -0.0054  0.29x  not a result     +0.3314  2.58x  RESULT (worse)  35W/91L
rama_nlp      -0.0037  0.20x  not a result     +0.2829  2.28x  RESULT (worse)  38W/88L
legacy        +0.0055  0.21x  not a result     +0.2813  2.29x  RESULT (worse)
ramah         -0.0004  0.02x  not a result     +0.1672  1.88x  RESULT (worse)
disto_risk    +0.0571  1.13x  RESULT (worse)   +0.2923  2.36x  RESULT (worse)
disto_mae     +0.0549  1.03x  RESULT (worse)   +0.2715  2.18x  RESULT (worse)
vbond_mean/sd +0.0041 == max   DEGENERATE: 204 of 204 branches tied (ideal geometry makes these constant)
```

> ### Every criterion identifies the WORST branch at 1.9–2.8× MDE with 5/5 folds. Not one identifies a BETTER branch than production already picks. **This is the third time in the sprint that exact asymmetry appears** — after the score prefix and after consensus. *Native-free observables have real skill at avoiding disasters and none at finding winners.*

**And the compute-matched control removes even the asymmetry.** On **GEN4-only** — production's own
four starts — **nothing is a result in either direction**, largest 0.46×. *Production's four starts
contain no branch bad enough for any criterion to be punished for. The large ARGMAX effects are a
property of having ~200 branches to find a bad one among, not of the criteria.*

**And the chiral hypothesis this lane was opened on is falsified — along with the ARGUMENT behind
it.** In-band ρ across the branch set:

```
d_to_C      +0.1122        obj1  +0.0742        disto_risk  +0.0702        obj0  +0.0603
rama_nlp    +0.0153   fold CI [-0.0461, +0.0608]   median -0.0740   positive on 54/126
```

> ### The ACHIRAL distance-map channels beat the chiral one FOUR TO ONE. Registered prediction P3.3 is falsified, and so is the reasoning that produced it: *"every native-free ranker is a distance-map function and therefore achiral, so a chiral criterion should see what they cannot."* **At this stage they see more, not less.**

**And lane R's first explanation of *why* was refuted by its own data**, which is the better result.
It proposed that the λ = 0.3 branches must all already be Ramachandran-plausible, leaving the
criterion nothing to discriminate. They are not: positive-φ spread **0.632**, **31.2% of branches
above the 17.5% rate of a completely unconstrained fit**, **126/126** targets carrying one, against a
real-library rate of 5.40%. ***The criterion has an enormous amount to discriminate and carries no
information anyway — Ramachandran plausibility and native proximity are orthogonal among the branches
this projection admits.*** That is S9-2's *"validity for free, no accuracy"* re-derived one level
down.

**An independent confirmation arriving from a different stage:** `disto_risk` and `disto_mae` ARGMIN
are **RESULTs at 1.03–1.13× in the WORSE direction** — *picking the branch the distogram likes best is
measurably worse than production's objective argmin.* That is §2's finding about the deployed score,
reproduced at the reconstruction stage by a job not looking for it.

**The search, accounted** (baseline **production**, not the grid mean — Appendix C): choosing the
criterion out of sample is worth **+0.0026 Å, CI [−0.0067, +0.0141]**. ORACLE per-target best criterion
−0.1000; best single arm 0.57× its own nominal MDE. **Honest comparison count ~60, not 64** —
`vbond_mean` and `vbond_sd` are **degenerate by construction**, 204 of 204 branches tied.

**Two checks came out in the lane's favour and are recorded as such.** The branches are **genuinely
distinct structures, not arithmetic noise** — 207.5 per target, 152.7 distinct at 1e-3 A, with
within-cluster `rmsd_nat` spread of median **0.00e+00** and p95 4.1e-04 — so the ULP discontinuity
does *not* dissolve the arm. The adversary went looking for that and did not find it.

> **The hypothesis was specific, the mechanism was named in advance, the falsifier was the right one,
> and it fired. That is the cleanest negative in the sprint and it is the coordinator's.**

### 5.4 Bit-exactness, the strong form

**All 630 per-target chain RMSDs (126 targets × 5 ladder rungs) reproduce S29's recorded values
bit-for-bit — max |Δ| exactly 0.000e+00** — from a different job, script and process. **The mean over
the 126 agrees to 1e-12**, not to the last digit: *an earlier draft of this section said production
"returns 3.210533994943299 to sixteen digits", and lane R's own verifier asserted that and FAILED* —
recomputing the mean in a different summation order gives `...300`. **The per-target identity is
exact and order-independent; a reduction over 126 float64s is not.**

> That sharpens rather than weakens the pairing with §1.2: **the operator is bit-reproducible per
> structure; what is fragile is the input's last bits — and now also any reduction taken over 126 of
> them.**

---

## 6. The quantum question, answered

**Charter §14 asks: can a genuine CVaR-VQE be designed whose quantum state and objective actually
contain information that can improve RMSD? On this instrument, NO** — and the reason is not a
failure of imagination, it is **chain length**.

Five conditions were scored on every discrete decision in the pipeline: **A** a space too big to
enumerate, **B** per-shot, target-dependent and native-free, **C** endpoint-relevant, **D** no cheap
classical algorithm, **E** a genuinely stochastic energy.

> **`n_res` is 9–16, mean 12.96, so `2^n_res ≤ 65536` on 126/126. Condition A fails by chain length
> alone.**

And the deployed stage does not satisfy the others either. **S31's target-invariance result survives
a falsification attempt built to break it** — a self-test that fires on a deliberately tied vector
(0.1221) and a deliberately unsorted one (3.3830) *before* touching data, unlike S31's own, which ran
on random floats that never tie. On the real 126: `sc[o]` non-decreasing 126/126; the **block model
reproduces `E` bit-for-bit, max error 0.0e+00, on all 126**; `max|E − ramp| = 0.0406840`. Priced at
**0.26× MDE**, with the selection readout **identically tied 126/126**, and an operator-matched
control — keep the block *sizes*, relocate the blocks — at **0.38×, the same size.** ***The residual
channel is the existence of duplicates, not which candidates duplicate: one integer per target.***

### The two properties a problem would need

> **P1 — a decision space that GROWS WITH THE TARGET.** Fragment assembly and per-residue branch
> selection are the candidates, and **both exist only where one fragment no longer spans the
> target.**
>
> **P2 — a genuinely STOCHASTIC energy.** A sampled free energy is the only candidate.
>
> **They must hold together, and this project has never had either.**

**That is the charter's "test on longer proteins" arriving as a DERIVED REQUIREMENT rather than a
suggestion** — see §7.

---

## 7. Longer proteins — the requirement is derived, and the instrument now exists

§6 turned the charter's *"maybe test on longer proteins"* into a **derived requirement**: the quantum
question fails condition A because `2^n_res ≤ 65536` at 9–16 residues. §2.0 sharpens it from the
other side: at `3n ≈ 39` a hull of 500 fragments is close to saturated, so the top rung prices a fit.
Both predict that something changes when the chain gets long.

### 7.1 The instrument

**`long40` — 45 targets, 44–60 residues (mean 54.71)**, five frozen folds of nine, assigned by greedy
balanced packing of whole identity-clusters **before the first arm ran**. Funnel: 293 in band → 170
monomers → 90 after leakage → **45** after redundancy clustering at identity 0.40. **Zero overlap with
`tuning126`; `benchmark60` never opened.**

It had to be built from `prots/`, and *why* is itself a finding: **the deployed library serves exactly
ZERO windows at `n ≥ 26`** — `peptide_db` caps at 25 and `fragment_db` at 20 — while `prots/` serves
**1.71 M** at `n = 45`.

**The MDE it can resolve**, built chain, realised: **1.29 Å** pool headroom, **1.45 Å** selection,
**0.48 Å** retrieval, **0.46 Å** sparse gain, **0.82 Å** filter skill — roughly **4× wider** than the
canonical instrument's — **and it varies per rung, so it is quoted per comparison rather than as one
constant.** **Anything under ~0.46 Å here is unresolvable and is reported as such.**

### 7.2 The ladder, BUILT CHAIN, both lengths, complete

One ladder, length-portable at every rung, run through the **same code** at both lengths; basis
**built chain (projector output), pre-AMBER** at both, so every rung *difference* is unaffected by
stage 4's absence.

```
rung                      L ~ 13 (n=126)     L ~ 55 (n=45)
pool_best     ORACLE          1.7078             5.0182
sparse_s10    ORACLE          1.1192             3.6648
top75_best    ORACLE          2.1004             5.8691
avg75                         3.4822             9.7450
avg75_random  control         3.6244            10.0622
```

**Independent reproduction:** `pool_best` **1.7078 matches contract rule 17 exactly** and `sparse_s10`
1.1192 matches its 1.1139 to 0.005 — *from different code, from coordinates*; the cloud arm returns
1.7108, the value pinned in `s12/instrument.py`'s own selfcheck.

```
difference (LOWER IS BETTER)          L ~ 13                    L ~ 55
headroom   avg75 - pool_best      +1.7744  5.50x  2/124     +4.7268  3.67x  1/44
retrieval  top75_best - pool_best  +0.3926  3.25x            +0.8509  1.79x
selection  avg75 - top75_best      +1.3818  4.74x            +3.8760  2.68x
sparse     sparse_s10 - pool_best  -0.5887  6.31x            -1.3534  2.97x
filter     avg75 - avg75_random    -0.1421  1.27x  74/52     -0.3171  0.39x  NOT A RESULT
```

**ALL THREE REGISTERED PREDICTIONS HOLD.** **P1** (pool headroom ≥ 0.75 Å): **+4.73 Å at 3.67× on
44/45 targets** — *generation is **less** the bottleneck at length, not more.* **P2** (selection is the
largest deployable rung): selection 3.876 against retrieval 0.851, a **4.6 : 1** ratio at L ≈ 55
against 3.5 : 1 at L ≈ 13. **P3** (projection cost tracks non-physicality): a real member projects for
**+0.0385** and a dense 75-member average for **+0.3356**.

> ### The readout's share of the recoverable loss is **77.9% at L ≈ 13 and 82.0% at L ≈ 55**. The shape is preserved and **the readout's share GROWS.** This sprint's conclusions are not a peptide-length artefact.

**But say *which* shape is preserved, because one rung does not replicate.** `filter_skill` — the
shipped `avg75` against a matched **random-75** control — is **−0.1421 at 1.27× MDE (a RESULT)** at
L ≈ 13 and **−0.3171 at 0.39× MDE (NOT A RESULT)** at L ≈ 55. ***BLOSUM's endpoint filter value does
not replicate at length.*** The shape is preserved **in the ORACLE rungs and the readout rung**, not
universally.

**And one absolute comparison that must NOT be made.** An earlier draft of this section said
*"averaging buys more in absolute terms at length (−2.58 Å against −1.16 Å)."* **That is
arithmetically right and misleading twice.** Everything at length is ~2.9× larger, so a 2.2× larger
absolute gain is **relatively smaller** (1.35 → 1.27) — *a number crossing a boundary its definition
does not cross.* And more seriously, **at length `m = 75` is past the optimum**: the cloud m-curve
minimises at **m = 3** (8.130) against 9.409 at m = 75, so the comparison is an artefact of the
`pool_mean` reference. **Correctly-sized averaging would buy much more, and m = 75 specifically is
1.28 Å worse than the best prefix.**

**But contract rule 16's CONSTANTS do not survive, only its ORDERING.** *"A real deposited member
projects for free, −0.0007 to −0.0030"* becomes **+0.0385 at L ≈ 55 — 13–55× larger.** The rule is
about the *object*, and that is what transfers; the numbers are peptide-length numbers.

### 7.3 Why the readout rung grows — measured, not asserted

```
                                          L ~ 13      L ~ 55
mean pairwise CA-RMSD among the 75         4.031      10.247
virtual CA-CA bond of the average          2.440       1.802
  ... as a fraction of the native bond     0.640       0.474
Rg of the average / Rg of the native       0.939       0.881
```

> **The object handed to the projector at L ≈ 55 has less than half a real backbone's bond length,
> because the set it averages is spread 2.5× further apart.** That is why the readout rung grows — and
> it is a property of **the operator meeting a wider set**, not of the pool.

### 7.4 Two things that ARE length-scoped, and one that is not

**`m = 75` is a peptide-length constant** (CLOUD basis, **a LEAD, not a result**). The curve is flat
over m = 30/50/75 at L ≈ 13 and has its minimum at **m = 3** at L ≈ 55, where the shipped m = 75 sits
**1.28 Å past it**; every fold's out-of-fold choice at length is m ∈ {3, 5}. Reported as a lead
because it is **0.72× MDE with W/L 21/24 beside a large mean** — the concentration warning — and it is
on the cloud, uncarried to the chain. **On `tuning126` the same arm says m = 75 is already right
(0.47×).** *It is a statement about what would have to change to deploy at length, not a proposal to
change the canonical pipeline.*

**The deployed distance prior is hard-capped at peptide length.** `core/predict.py` sets
`MAXLEN = 26` and `SEP_BINS` tops out at 24, so **every pair with |i−j| ≥ 24 collapses into one
terminal bin** that in training held only |i−j| ∈ {24, 25} — at n = 55 that is **27% of all pairs**.
`sep/26.0` and `n/26.0` reach 1.7–2.3, and the MLP carries raw `n` and raw `j−i`, fitted only on
n ∈ [8, 26]. ***It is evaluated outside its fitted support by construction.*** This is why the
canonical ladder's distogram-defined rungs cannot be evaluated at length at all, and **retraining it
is the single largest named piece of work a long deployment needs.**

**The common mode is NOT length-scoped.** `f` = 0.4683 (SE 0.0163) at L ≈ 13 against **0.5303**
(SE 0.0269) at L ≈ 55 — difference **0.70× MDE, NOT MEASURED**, and if anything it *rises*. The
bias-variance identity verifies to **1.6e-15**. *Definition, because this is exactly the boundary that
gets crossed:* both are the **BLOSUM** top-75; the 0.676 on record is the **DISTOGRAM** top-75, a
different object, never differenced against it.

### 7.5 What the lane falsified, including its own hypothesis

**L-H1 is FALSE: the representation is not the obstacle at 40–60 residues.** The lane pre-registered
that the ideal-geometry representation would be the barrier at length, with the falsifier *"< 1.0 Å at
L = 45"*. The **deployed projector** run on the native itself (ORACLE / NOT DEPLOYABLE):

```
arm                n     L        mean     median    SE       p90     max
PROJ_NAT_SHORT    126   9-16     0.0429    0.0228   0.0053   0.104   0.431
PROJ_NAT_LONG      60   41-60    0.7002    0.7444   0.0295   0.956   1.129
```

**Mean, median *and* p90 all sit below the falsifier.**

> **And the by-product outlives the hypothesis.** The motivating quantity — the **native-torsion
> rebuild** — rises as `0.0299 · L^1.239` (R² 0.9958, paired within 371 molecules) and **is not the
> representability floor.** It overstates the distance to the emittable set by **3.9× at L ≈ 52** and
> **8.1× on the canonical 126** (0.347 against 0.0429). **`fragment_db.REBUILD_TOL = 1.0` and
> `data.REBUILD_TOL = 1.5` gate library admission on exactly this quantity**, on the rationale that
> above it *"the deposited geometry carries something the representation cannot express."* **That does
> not follow — members are being excluded that the projector can express ~8× better than the gate
> assumes.** Not acted on, because changing library admission would change `tuning126`'s pools.

**Two more of the lane's own hypotheses died.** Its exploratory claim that *the sequence channel
strengthens with length* is **NOT SUPPORTED** — in-band Spearman inside the K=500 pool is −0.0464 at
L ≈ 13 and **+0.0261** at L ≈ 55, the *wrong sign* at length, contrast 0.83× NOT MEASURED; so
`structure-and-sequence-are-decoupled` **survives at length**. And its pre-registered admission filter
(`identity ≥ 0.4` against the banks) **rejected 140 of 170 monomers and left ZERO targets** — caught
**by its own null before any RMSD existed**, which showed it rejects **100% of real and 100% of
shuffled** sequences alike.

### 7.6 What was NOT run

**The donor-pool control at long length was not run.** §2.0's finding — that 500 fragments from a
*different protein* reach the target's hull at 9–16 residues — predicts that the same control should
**FAIL** at `3n ≈ 164`. The script exists (`s32/s32_L_hull_capacity.py`) with the prediction registered
in its header. **It ships as a prediction, not a measurement**, and §10 names it as the first thing
S33 should run.

---

## 8. What was falsified, including by its own author

**Registered falsifiers that fired against the lane that wrote them:**

| prediction | registered by | outcome |
|---|---|---|
| a **chiral** criterion can pick the projection branch where achiral ones provably cannot | **the coordinator** | **`rama_nlp` — 0.0031 at 0.16× MDE.** Nothing of 64 comparisons reaches 0.7× except in the wrong direction |
| the random gate improves the endpoint | coordinator (predicted *unchanged or worse*) | **held** — +0.1648 at 0.92×, 0 of 8 draws better |
| dilating the cloud to ideal bond length helps (P1.3) | lane R | **FALSIFIED**, +1.0425 at 2.79× WORSE |
| relaxation displacement points at the target, median `cos ≥ +0.10` | lane D | **FAILED**, 0.09–0.73×; indistinguishable from a matched random direction |
| chiral in-band variance share < 15% | lane D | **FAILED at 23.6%** |
| the 9–16-residue scope excuse (one basin / frozen) | lane D | **FALSIFIED on both horns** — the excuse is *not available* |
| READOUT is the largest decomposition cell | lane P | **partly failed** — readout and ranking are **tied** (−1.05 vs −1.09) |
| the operator law's set-best ratio ≥ 5× | lane P | **FAILED** — measured **2.85×**; the set best *does* reach the output |
| the negative in-band ρ is a collider artefact | lane P | **FAILED** — the collider alone gives the **wrong sign** |
| spread-maximisation above a score floor helps | coordinator | **closed by its own control** — `SPREAD − RANDFLOOR ≈ 0` |
| the chirality/filter synthesis | **coordinator** | **CIRCULAR — controlled and refuted before it was written down.** Lane D's 0.271 and lane V's 0.272 are one number, not two |

**Three theorems closed more than any measurement did.**

- **D-E**: an achiral force field's stochastic propagator is `O(3)`-equivariant, so any achiral
  invariant averaged over it is itself an achiral functional of the seed — a distance-map reading.
  ***Charter §23 entire*** (MD, ensembles, basin populations, metastability, transition rates,
  autocorrelation, dynamic modes) **and §24's ensemble reweighting and temperature response are closed
  as scalar rankers.** *Stochasticity is not an escape: noise is variance, not information.*
- **D-F**: `⟨F(x), G(x)⟩` for two equivariant fields is invariant, hence achiral, hence closed — **every
  physics/prior consistency scalar**, including force-vs-prior-gradient alignment.
- **And AMBER is not a chiral scorer.** ff14SB+GBn2 is **reflection-invariant**. ***This project has
  no chiral scorer at all***, which retires the premise the coordinator opened lane R on.

**Also closed with a named mechanism, not merely "did not work":** AMBER has **no resolution
advantage** — a candidate stores only a Cα trace and φ/ψ, sidechains are modal rotamers, so every atom
AMBER sees is `Ψ(seq, φ, ψ)`; converged relaxed AMBER in band is **−0.0089** against the single
point's +0.0006, so **the zero is not a clash artefact and minimisation does not rescue it**; and
relaxation's 45,867 → −483 kcal/mol for **0.094 Å** of Cα motion is a deterministic rotamer-placement
artefact carrying no candidate information.

> **Every use of a force field on this pipeline is now closed — as a RANKER by theorem and by
> measurement, and as a MOVER by measurement.**

**One free by-product worth more than the arm that produced it:** `cos(e_prod, e_pool75) = +0.9443`.
**The production chain's own error is 94% aligned with the pool's common mode** — S31 §20.1 confirmed
from a different object by a job that was not looking for it.

---

## 9. Is it possible to lower the RMSD?

The charter asks for a direct answer. **Yes — but not by any route this pipeline's architecture makes
available, and the sprint can now say precisely why.**

### 9.1 What would have to be true

Every stage-level intervention reduces to the **same requirement**, and that collapse is this
sprint's central result rather than a coincidence:

- The readout's optimisation is the **projection of the native onto the candidate hull**, with gain
  exactly 1 on the active subspace (§3.1). To use it you must know `a` on ~5 active directions.
- `a` and the common mode `μ` are **one object** under a native-free affine bijection, proved
  independently by two lanes in opposite directions (§3.2). *A per-candidate quality estimator with
  in-band skill **is** a structure predictor.*
- The deviation space has `rank(d) = 3n − 6` **exactly on 126/126**, so the requirement is **exactly
  `3n − 6` ≈ 39 real numbers and there is no compression hiding in the pool's geometry** (§3.2).
- And the hull's 2.10 Å of apparent headroom is **not something retrieval found** — a pool built for
  a *different protein* reaches it equally well (§2.0).

> ### The missing quantity is ~39 real numbers that ARE the answer. Every "channel" this project has looked for is a way of obtaining some of them, and every one measured is worth a few hundredths of an Ångström.

### 9.2 What was actually tried, and what each cost

```
intervention                                  result                         why it failed
random window instead of the score's          +0.1648  0.92x  NOT MEASURED   the score buys the set MEAN
per-target sign, ORACLE, one bit              -0.2101  1.66x  ORACLE only    26% of the direction's value
per-target sign, native-free                  at or below the marginal       a 81/19 feature can't carry 59/41
branch selection, 64 criteria both ways       best 0.51x, search-adj +0.0028 skill at finding the worst only
chiral Ramachandran branch criterion          -0.0031  0.16x                 the achiral twin matches it
scalar dilation (ideal bond)                  +1.0425  2.79x  WORSE          contraction is separation-dependent
scalar dilation (Rg-matched)                  +0.0622  1.47x  WORSE          one scalar, one moment
relaxation displacement, applied              +0.02 to +0.10 at every step   indistinguishable from random
spread-maximisation above a score floor       SPREAD - RANDFLOOR ~ 0         the dispersion term does no work
exact solution of the CVaR objective          -0.0112  0.19x  (S31)          forced by gain 1 (section 3.1)
```

**Not one arm reached its own MDE in the helpful direction. The endpoint is unmoved at 3.2105 Å.**

### 9.3 The three things that would actually move it

**1. A structure predictor** — which by §3.1 makes the pool redundant rather than better. If you can
estimate `t` to better than the hull floor `d = 1.8290`, emit it; the readout adds nothing. *This is
not a criticism of the readout; it is what gain 1 means.*

**2. A longer instrument, where the geometry changes** (§7). At `3n ≈ 164` the hull no longer
saturates, so retrieval starts to determine the *ceiling* and not merely the pool mean; and
`2^n_res` becomes non-enumerable, which is the one condition the quantum question failed on. **Both
of the sprint's two central negatives are length-conditional, and both were derived rather than
guessed.**

**3. Averaging more, or better.** It is the only operator here that demonstrably works — pool mean
4.4533 → emitted 3.2928 on the short instrument, and **−2.58 Å on the long one.** The readout program
says it wants **low `⟨w,a⟩` and high `w'Bw`**; production has no control over the first and destroys
the second at the filter. *That is the one place where the sprint's own theory points at an untested
arm rather than away from one* — and lane P's re-measured operator law gives the coefficient:
`out = −0.9934 + 0.9219·set_mean + 0.3232·set_best`, **ratio 2.85**, R² 0.9162. **The set best is not
unreachable — it reaches the output about a third as hard as the set mean**, which contradicts the
lane's own registered prediction of ≥5× and is the one arm this sprint's data argues *for*.

### 9.4 The honest bottom line

**On `tuning126`, with this architecture, at this chain length: no.** The 2.10 Å is a fit, not a
retrieval; collecting it costs the answer itself; and every native-free channel measured across four
sprints returns a few hundredths.

**The primary target of < 3.00 Å is 0.21 Å away and it is not 0.21 Å of engineering.** It is 0.21 Å
of information that nothing in the pool, the physics, the geometry or the quantum stage has been
shown to contain.

---

## 10. The next bottleneck

**It is the same one, stated more precisely than before, and the sprint's contribution is that it is
now a quantity rather than a direction.**

> ### `3n − 6` ≈ 39 real numbers per target, which are the answer. Not a channel, not a bit, not a score.

Four sprints have looked for *a signal*. This one showed that every candidate signal is a partial
observation of the same ~39-dimensional object, that the readout transcribes it at gain exactly 1
without amplification, and that the pool's apparent 2.10 Å of headroom is the **capacity of fragment
space at this length** rather than anything retrieval discovered.

**What follows for S33, in priority order:**

1. **Go longer.** Both central negatives are **length-conditional and derived, not guessed**: the
   quantum question fails on `2^n_res ≤ 65536`, and the hull saturates because `3n ≈ 39`. At 54.7
   residues `3n ≈ 164`, both change. **The instrument exists** (45 targets) and carries a sharp
   falsifiable prediction: **the donor-pool control of §2.0 should FAIL there.** Run that first — it
   is cheap, and it decides whether this project's architecture has a regime where it is not
   saturated.

   **And the one named engineering prerequisite:** the deployed distance prior is **hard-capped
   at peptide length** — `core/predict.py` sets `MAXLEN = 26`, `SEP_BINS` tops out at 24, so at
   n = 55 **27% of all pairs collapse into a single terminal bin** that in training held only
   separations |i−j| ∈ {24, 25}, and the MLP carries raw `n` and raw `j−i`, fitted only on
   n ∈ [8, 26].
   ***It is evaluated outside its fitted support by construction.*** **Retraining it is the
   single largest named piece of work a long deployment needs**, and until it is done the
   distogram-defined rungs of the ladder cannot be evaluated at length at all.
2. **Weight toward the set best.** The re-measured operator law is
   `out = −0.9934 + 0.9219·set_mean + 0.3232·set_best`, **ratio 2.85, R² 0.9162** — against a
   registered prediction of ≥ 5×. **The set best is not unreachable; it reaches the output about a
   third as hard as the set mean.** *This is the one arm the sprint's own data argues for rather than
   against*, and it is untested.
3. **Stop looking for native-free in-band skill as posed.** It exists, it is typicality, and the
   deployed operator is already its argmin. A new channel must beat the hull floor `d = 1.8290` — at
   which point it makes the pool redundant, which is a different architecture, not an improvement to
   this one.
4. **Do not re-open**: the prefix length `m`; scalar dilation in any calibration; branch selection by
   any of the sixteen criteria; force fields as rankers or movers; the per-target sign as a route;
   ensemble/dynamics scalars (theorem D-E).

---

## Appendix A — every claim withdrawn this sprint

**Sixteen defects were found in work already written down. The coordinator wrote the largest share,
and most were caught by machine rather than by memory.**

| claim | whose | how it died |
|---|---|---|
| *"the retrieval filter costs +0.4357"* | coordinator | **Mislabelled.** 500→128 is the **distogram score prefix** — `s29_O_ladder.py`'s own *"the quantum field of view"*, 128 = 2⁷. Retrieval is the earlier arrow |
| *"the score is worse than random at retaining the best candidate"* | coordinator, from lane V | **Entirely 18 circular targets.** On the other 108 the score is marginally **better**; `n_better = 72` means the score wins on 72 of 126 and I read the W/L backwards; the median **−0.0380** also favours the score and I did not report it |
| the ladder line *"+0.4350 → 2.1435"* | coordinator | **Basis error** — +0.4350 is the **member** basis, the chain increment is **+0.4357**. Caught by AUDIT 9 on its first run: *every ladder block's increments must close on ONE basis* |
| *"gain exactly 1 means no noise suppression"* | coordinator + lane Q | **Backwards.** Gain 0 on ~33 dimensions is **TOTAL** suppression — ~87% of a generic error annihilated; at ε = 4.0 a 3.66 Å estimate emits at 2.23 Å. The headline survives on the **hull floor** `d = 1.8290` |
| the convexity check for S32-L2 | coordinator | **VACUOUS** — accumulator initialised at the pass threshold, max over always-negative quantities, **could not fail on any input.** Written one hour after I wrote the rule forbidding exactly that |
| *"sparse combinations barely leave the manifold"* | coordinator | **False.** Sparse sits **farther** off-manifold (0.8503 vs 0.8150) and pays 800× less. The tax is **direction**, not distance |
| *"the final rung runs at λ = 0, the prior is inactive"* | coordinator | **False.** The ladder is (0.0 → 0.3), the canonical arm is **λ = 0.3**, and the prior is already selecting among branches |
| *"λ = 0.3 costs +0.0107 at the endpoint"* | coordinator | **A cross-job difference, not an effect.** Measured in one job on one cloud: **+0.0055, 0.35× MDE, NOT MEASURED** |
| *"typicality is a RESULT in the wrong direction"* | coordinator, from lane V | **Inverted.** `typicality` is a **distance**, so the arm was the **least** typical branch. Most typical: −0.0054, 0.29×, a null |
| *"a CHIRAL criterion can pick the branch where achiral ones cannot"* | **coordinator — the hypothesis the lane was opened on** | **FALSIFIED, and the argument itself fails.** In-band ρ: `d_to_C` **+0.1122**, `obj1` +0.0742, `disto_risk` +0.0702 — against `rama_nlp` **+0.0153, fold CI [−0.0461, +0.0608], median −0.0740.** ***The achiral distance-map channels beat the chiral one four to one*** |
| *"the chirality/filter synthesis"* | **coordinator** | **CIRCULAR — controlled and refuted before publication.** Lane D's 0.271 and lane V's 0.272 are **one number, not two** |
| *"the λ=0.3 branches must all already be Ramachandran-plausible"* | lane R — its own explanation of the above | **Refuted by its own data**: positive-φ spread **0.632**, **31.2%** of branches above the 17.5% unconstrained rate, **126/126** targets carrying one. *The criterion has enormous variation to work with and carries no information anyway* |
| *"3.210533994943299 to sixteen digits"* | lane R | **Its own verifier asserted it and FAILED.** Per-target is exact (630 RMSDs, max Δ = 0); a **reduction over 126 float64s** depends on summation order |
| `n_distinct` as a structural claim | lane R, caught by lane V | **A greedy-seed count, not connected components.** Fixed by union-find; **39 of 126 targets changed** |
| `split_half_transfer` on a 16-criterion grid | lanes V **and** R, independently | **Baseline artefact** — it centres on the grid's column mean. Read −0.0254 with a CI excluding zero; against production **+0.0026, CI centred on zero** |
| the READOUT cell is the largest; set-best ratio ≥ 5× | lane P | **Both failed.** Readout and ranking are **tied** (−1.05 vs −1.09), and the measured ratio is **2.85×** — *the set best is not unreachable* |

**Plus the reproducibility class, which was invisible to every prose-level check.** AUDIT 11 asks
whether any `.py` actually **writes** a given artefact. It found **seven artefacts with no provenance
block and no producing script**, two of which carried §5.3. Lane R then **ran the same audit over its
own lane and found two more of its own** — four in one lane, all first computed in inline shell
heredocs. Lane D's `s32_D1_signtransfer.json` was the **third** instance of this shape in the
project's history and **the first caught by machine**; `s32_R_gen4d_start.json` had a committed script
that **had never actually run**, dying on a `KeyError` before writing.

> **Every one is the same shape as S31's: a number re-used across a boundary its definition does not
> cross** — a basis, an object, a moment, a stratum, a direction, a baseline, a job.

---


## Appendix B — multiplicity and the search that was run

**Six lanes, 24 ledger entries, ten pre-registrations, 56 result artefacts.** Comparisons emitted:
lane P **38**, lane D **92**, lane V **32 on the chain basis plus 8 in-band**, lane R **64** in the
branch family alone.

**The sprint-level accounting, walked over every `effect_over_mde` in every artefact:**

```
improvements clearing 1.0x MDE on the BUILT CHAIN, deployable:        0
improvements clearing 1.0x MDE on any basis:                         30
   of which ORACLE (read the native to select):                     ~20
   of which filter-vs-random on the POOL MEAN (a diagnostic):         4
   of which in-band Spearman diagnostics:                             6
```

**The three largest ORACLE ceilings**, all **NOT DEPLOYABLE**: ORACLE top-5 through the deployed
operator **−1.567 (5.17×, cloud)**; ORACLE per-target scale **−0.1352 (3.84×, chain)** — but **122%
accounted by the across-target order statistic**, and its leave-fold-out factor is +0.0044 at 0.14×;
ORACLE best branch **−0.1178 (chain)** — **3.3% survives a split half**.

**Searches priced rather than quoted:**

- The branch family: **64 comparisons** (honestly ~60 — four are degenerate with 204 of 204 tied).
  Best single arm **0.57× its nominal MDE**; **out-of-sample search value +0.0026, CI [−0.0067,
  +0.0141]**.
- The prefix family: every best-member contrast against an alternative prefix is **entirely FAIL18**,
  three instances — random-128 (+1.4879 / −0.0296), random-75-of-128 (+0.3611 / +0.0210), BLOSUM-75
  (−1.5871 / +0.0288). **In each the aggregate clears or approaches MDE and the non-circular 108 is
  NOT A RESULT.** This is now a **registered general property of the instrument**, not three
  coincidences, and the verifier flags any line quoting one without its strata.
- Draw distributions, never best draws: random-128 **0 of 8 draws better** (draw sd 0.0187);
  random-75 **0 of 8** (0.0381); random branch **300 draws**, ±0.0076.

---

## Appendix C — methodology this sprint had to learn

Six rules were added to the contract mid-sprint, each paid for by a defect. Three are worth carrying
beyond this project.

**1. `split_half_transfer` nulls the wrong question unless you give it production as the baseline.**
It centres on the grid's **column mean**, not the incumbent. On a 16-criterion grid containing one
implausible column it read **−0.0254 with a CI excluding zero** — an effect nearly 3× the best single
arm, entirely an artefact of the baseline. Done against production: **+0.0004, CI [−0.0091, +0.0116]**,
a CI centred on zero. *Caught by the lane that wrote it, before it reached the report.*

**2. Bit-identity, not value-identity, licenses a cross-job chain comparison.** A **one-ULP
(7.1e-15 Å)** change in the input cloud moves the built chain **0.10–0.15 Å**. The operator is
perfectly deterministic and **discontinuous** in its input. Measured irreducible sd on an unpaired
built-chain mean: **±0.003 Å**. Paired contrasts from bit-identical clouds in one job are unaffected;
**nothing else below ~0.03 Å is resolvable.** And a *reduction* over 126 float64s is not exact either
— the per-target identity is, the mean is not, which a lane found by having its own verifier fail.

**2b. “Does a script write this file?” is not enough — the quoted number must be INSIDE it.** AUDIT 11 asks whether any `.py` writes a given artefact, which caught seven unreproducible files. But **two of lane D's five were worse than unbacked: their headline numbers were not in the file at all.** `s32_D0X_circularity.json` and `s32_D4X_rgsign.json` stored only per-target rows, and the cited variance-retention table and the *“positive on 81% of targets”* figure existed **only in stdout** — so a reader following the citation would have found a 126-row file that did not contain the cited number. **AUDIT 11 as written would have PASSED on both once a script existed.** The standing form of the rule is therefore: ***a `.py` writes this name AND the quoted number is a numeric leaf inside it.*** *The lane caught itself twice, and the second catch is the one the weaker form misses.*

All five were then re-emitted from `s32/s32_D_reemit.py` with seeds and RNG draw **order** pinned, and **3,985 shared numeric leaves reproduce at max |diff| exactly 0.000e+00.** Each file carries a `reemit_check` block recording that diff, so the check is re-runnable rather than a claim. **No number in this report changed.**

**3. An audit that cannot fail is decoration; an audit whose scope shrinks silently is worse.** Both
happened here. A convexity check initialised its accumulator at the pass threshold and reported
`+0.000e+00` in every regime including ones where the true value is −0.26. And a global string
replace narrowed the path audit from 71 paths to 35, **visible only as the matched count falling** —
*an audit whose scope narrows without failing is the same class of defect as one that cannot fail.*
Both now carry positive controls; the path audit asserts its own coverage by file extension.

**4. Retractions in place versus phrase-matching audits.** Contract rule 13 requires struck claims to
stay visible, so a phrase-matching audit fires forever on correctly-handled retractions. The fix is a
**positive obligation, not an exemption**: a struck claim is exempt **only if its replacement is
stated in the same block**, with a self-test proving the same claim asserted *live* below the block is
still caught. That turns the audit from *"stops nagging"* into *"starts catching retractions that
strike a claim without saying what replaces it"*.

**5. A job name is not a lock, and killing a wrapper does not kill its child.** Duplicate launchers
raced on two shards and wrote 24 duplicated rows; a lane found the **orphaned workers** after the
wrappers were killed. **Check the pid, not the name** — and assert **both** `len(rows) == 126` **and**
`len(set(pdbs)) == 126`, because the second alone would have passed.

**6. Staging is not lane-isolated.** Concurrent commits collided on `.git/index.lock`, one lane's
commit swept in another's staged files, and a lane's commit **silently reverted a coordinator edit
that had not yet been committed — twice.** Write and commit in one call, and **verify the edit is in
`HEAD`, not just on disk.**

