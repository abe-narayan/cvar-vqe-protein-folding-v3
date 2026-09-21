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

**Consequence:** the anchor carries ~**±0.002 Å** of pure arithmetic noise. The charter's targets stay
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

**The pool is not the bottleneck.** The existing K=500 pool supports **1.1139 A** through a sparse
convex combination of about ten members. Everything downstream of retrieval destroys **2.10 A that is
already present**. Charter §56 asks for the earliest irreversible loss; **it is not candidate
generation.**

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

## 3. The readout is a hull projection, and that closes a family

[PENDING — S32-L5, L(Q1)–L(Q3): gain, the sufficient statistic, the hull floor, monotonicity.]

---

## 4. In-band skill is not zero — the per-target SIGN is missing

[PENDING — S32-L7, L(D2), and whether anything native-free reads the bit.]

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

And the sixteen native-free criteria, x {all branches, GEN4-only} = **32 comparisons**:

```
best of 32:   d_to_C      -0.0095   0.56x   5/5   58W/43L    NOT A RESULT
              obj0        -0.0071   0.42x        rg        +0.0052  0.18x
              rama_nlp    +0.0006   0.03x        ramah     -0.0007  0.04x
              disto_risk  +0.0562   0.98x  NOT MEASURED, wrong direction
              typicality  +0.3434   2.35x  RESULT -- in the WRONG direction, 27W/75L
```

**`rama_nlp` — the chiral Ramachandran criterion the lane was opened for — lands at 0.03x MDE.**
Nothing reaches even the 0.7x band except in the wrong direction, and that is *before* charging the
search: the best of 32 comparisons has an achieved MDE well above its nominal one. On the
**compute-matched** GEN4-only subset — production searches 4 starts, so an arm searching ~200 and
winning is a bigger search rather than a better selector — the best is 0.42x.

***The `typicality` arm is load-bearing as a positive control***: the machinery *can* emit a RESULT,
so the sixteen nulls are a measurement rather than a broken pipeline.

**But read its direction carefully — `typicality` in these rows is a DISTANCE** (mean Cα-RMSD of a
branch to the production top-75), so **lower means more typical**, the sense the consensus medoid
uses. The adversary maximised it. Both directions, from the raw branch rows in one script (n = 111,
tie-averaged argmin, paired to production in-job):

```
ARGMIN typicality = the MOST typical branch    3.1499   -0.0041   0.21x MDE   62W/48L   NOT MEASURED
ARGMAX typicality = the LEAST typical branch   3.4957   +0.3417   2.42x MDE   31W/80L   WORSE
branch-set MEAN (zero-skill reference)         3.1589   +0.0049
PRODUCTION                                     3.1540
```

**So the RESULT is that the most ATYPICAL branch is much worse** — the *expected* direction, and it
**supports** consensus-as-outlier-avoidance rather than contradicting it. *The most typical branch is
a plain null at 0.21x and belongs with the other fifteen.* An earlier draft of this section called it
*"a RESULT in the wrong direction"*; **that was an inverted reading and is corrected here.** The trap
was a column named after the property it is *inversely* related to — lane R has stated the direction
explicitly in every artefact rather than renaming a column other lanes are already reading.

**Two checks came out in the lane's favour and are recorded as such.** The branches are **genuinely
distinct structures, not arithmetic noise** — 207.5 per target, 152.7 distinct at 1e-3 A, with
within-cluster `rmsd_nat` spread of median **0.00e+00** and p95 4.1e-04 — so the ULP discontinuity
does *not* dissolve the arm. The adversary went looking for that and did not find it.

> **The hypothesis was specific, the mechanism was named in advance, the falsifier was the right one,
> and it fired. That is the cleanest negative in the sprint and it is the coordinator's.**

### 5.4 Bit-exactness, the strong form

**All 630 chain RMSDs (126 targets x 5 ladder rungs) reproduce S29's recorded values bit-for-bit**,
mean and max |Δ| exactly **0.000e+00**, from a different job, script and process. Production returns
**3.210533994943299** to sixteen digits; the cloud **3.048338093879531**.

---

## 6. The quantum question, answered

[PENDING — charter §14, the five conditions, P1 and P2, and why chain length is the binding one.]

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
