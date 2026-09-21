# Sprint 32 — Where the RMSD is lost, and what it would take to get it back

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. Assembled as
results land so nothing is reconstructed from memory at the end.

Branch `s26` · instrument `tuning126`, 126 targets, 9–16 aa · endpoint **mean built-chain Cα RMSD**
Charter `s32/BRIEF.md` (verbatim, 2,087 lines) · Ledger `s32/LEDGER.md` · State `s32/STATE.md`
Contract `s32/S32_CONTRACT.md` · Causal map `s32/CAUSAL_MAP.md` · Verifier `s32/s32_verify.py`
Multiplicity `s32/MULTIPLICITY.md` · Theory `s32/THEORY_Q.md`

---

## 0. The answer, up front

[PENDING — written last.]

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
**−0.2101 Å, 1.66× MDE, 5/5 folds — ORACLE / NOT DEPLOYABLE.** Against S31's −0.8102 Å for the
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

> **⚠ PROVENANCE PENDING — do not finalise this section until it clears.** Lane V's AUDIT 11
> finds `s32_R_dilation_cloud.json` and `s32_R_dilation_rg_cloud.json` have **no provenance block and
> no script in the repository that writes them**. These are two RESULT-grade endpoint numbers and a
> registered falsification that **nobody can currently re-run**. The 22.15% / 5.40% reconciliation
> below is a genuinely valuable output and must not rest on an artefact with no producer. *Lane D hit
> the same defect on D1-T and fixed it by committing the producing script; this needs the same.*

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
in the same job, tie-averaged, independently of lane R (n = 103 at the time of writing; production
in that job is **3.1865** over the subset, and every arm is paired to *it*):

```
ORACLE best branch                    3.0768   -0.1097 vs production   ORACLE / NOT DEPLOYABLE
  SPLIT-HALF TRANSFER                          -0.0051  =  5.0% of the oracle  ->  NOT A SIGNAL
random branch, 300 draws              3.1914   +0.0049  (draw sd 0.0079)
```

> ### Production's multi-start argmin is worth 0.005 A over picking a branch with a coin. Any branch-selection rule is competing for 0.11 A of oracle headroom against an incumbent that is 0.005 A better than random — and **95% of that 0.11 A is an order statistic that does not survive a split half.**

And the sixteen native-free criteria, tested in **BOTH DIRECTIONS** on the full branch set and on
production's own four starts — **64 comparisons**. *Testing one direction only would have been an
unregistered choice that halves the apparent multiplicity, and it is what produced the inverted
gloss an earlier draft of this section carried.* (n = 123; production in-job 3.2198.)

```
                 ARGMIN (best-looking)            ARGMAX (worst-looking)
d_to_C        -0.0082  0.51x  not a result     +0.3417  2.62x  RESULT (worse)
rg            -0.0037  0.12x  not a result     +0.3371  2.82x  RESULT (worse)
typicality    -0.0038  0.20x  not a result     +0.3406  2.61x  RESULT (worse)
obj1          -0.0058  0.30x  not a result     +0.3503  2.72x  RESULT (worse)
rama_nlp      -0.0031  0.16x  not a result     +0.2908  2.31x  RESULT (worse)
legacy        +0.0079  0.30x  not a result     +0.2887  2.30x  RESULT (worse)
disto_risk    +0.0577  1.12x  RESULT (worse)   +0.2986  2.36x  RESULT (worse)
ramah         +0.0009  0.04x  not a result     +0.1731  1.91x  RESULT (worse)
```

> ### Every criterion's ARGMAX is a large, 5/5-fold, 2–3× MDE RESULT in the WORSE direction, and every criterion's ARGMIN is nothing. **The criteria carry real information about branch quality — they reliably identify disasters — and none of it converts into finding a winner.**

**And the compute-matched control kills even that.** On **GEN4-only** — production's own four starts —
**nothing is a RESULT in either direction**; the largest is 0.46×. *So the big ARGMAX effects are a
property of having ~200 branches to find a bad one among, not of the criteria. Production's four
starts contain no branch bad enough for any criterion to be punished by.*

**The chiral Ramachandran criterion the lane was opened for is −0.0031 at 0.16×.**

**The search, accounted** (baseline **production**, not the grid mean — see the note on
`split_half_transfer` in Appendix C): choosing the criterion out of sample is worth **+0.0028 Å, CI
[−0.0054, +0.0128]** — *a CI centred on zero.* The ORACLE per-target best criterion is −0.0983 and
the best single arm is 0.51× its own nominal MDE.

**Honest comparison count: ~60, not 64.** `vbond_mean` and `vbond_sd` are **degenerate by
construction** — every branch has ideal geometry, so 203 of 203 branches tie — and `ramah` min (58
tied) and `posphi_frac` min (92 tied) are heavily degenerate. Tie-averaging handles them correctly.

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

## 7. Longer proteins

[PENDING — lane L, and why the requirement is derived rather than suggested.]

---

## 8. What was falsified, including by its own author

[PENDING — the registered falsifiers that fired, and the corrections.]

---

## 9. Is it possible to lower the RMSD?

[PENDING — the charter's closing question, answered directly.]

---

## 10. The next bottleneck

[PENDING]

---

## Appendix A — every claim withdrawn this sprint

[PENDING]

## Appendix B — multiplicity and the search that was run

[PENDING]
