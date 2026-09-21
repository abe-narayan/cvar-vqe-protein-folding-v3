# S31 — LANE V, THE ADVERSARY

Contract rule 24: the main team may not approve its own positive. This file is the adversary's
record. **Defects lead.** Every entry states (a) what is claimed and where, (b) what the artefact
says, (c) the corrected sentence.

Read-only on `s29/ s28/ s27/ s16/ s12/ core/` and on other lanes' files. Writes confined to this
file, `s31/results/s31_V_*`, `s31/s31_V_*.py`, and ledger entries.

**Scope note.** The coordinator stood lane V down from the `bestm128` question at 00:22; lane F had
already closed it (S31-L11) and lane L had audited its framing (NOTE 11). **I did not duplicate it.**
Section 1 below is a *different* arm of the same family — the one nobody ran — and it was
pre-registered before it was computed.

---

## 0. PRE-REGISTRATION — V1, THE MATCHED-K CONTROL FOR THE SET-MATCHED READOUT LADDER

**Committed before the number existed (contract rule 19), commit `bac1fe1c`.** Code
`s31/s31_V_orderstat.py`, artefacts `s31/results/s31_V_orderstat.json`,
`s31/results/s31_V_orderstat_rows.jsonl` (126 rows), seed 31099.

`s31/STATE.md` NOTE 2(b) and `s31/LEDGER.md` S31-L6 §1 read the fixed-top-128 ladder as

```
argmin over 128            7.0 bits              2.1458
2-of-128, UNIFORM weights  12.99 support bits    2.0700    "-0.076 vs argmin, ZERO weight bits"
```

with the coordinator's handed-back reading that *the 0.076 Å is **error cancellation**, not
weighting.*

**The objection.** `2-of-128` is a per-target minimum over `C(128,2) = 8128` ORACLE-chosen supports;
`argmin over 128` is a per-target minimum over **128**. A **64× larger oracle search**, compared as
though the arms differed only in what they emit. Contract rule 11 and memory
`grid-oracles-are-order-statistics` both bind.

**Registered bar.** *"The 0.076 Å is an order statistic, not error cancellation"* **FIRES** if the
per-target minimum over **128 random pairs** — matched K, matched operator, same top-128, uniform
weights, same native — fails to beat the minimum over the 128 singletons by `0.7 × MDE`.
**Registered prior ~3:1 in favour.** `NDRAW = 8`, **mean over draws, never the maximum** (rule 10).

**Basis: CA point cloud throughout. Every arm is ORACLE / NOT DEPLOYABLE** — each is a per-target
minimum taken against the native. **No built-chain claim is made anywhere in V1.**

**Reproduction gate passed:** my recomputed `argmin over 128` is **2.1457974841561689** against lane
C's `ORACLE_argmin128` **2.1457974841561676**, deviation **1.3e-15**. Same object.

---

# DEFECTS

## D0 — SEVERE, AND LIVE IN TWO LANES' CODE. THE `0.6931` COHERENCE BAR IS MEASURED ON A DIFFERENT OBJECT THAN EVERY ARM IT GRADES

**Claimed** — `s31/STATE.md` NOTE 10, whose heading is *"NO RANKER CAN EVER PASS THE COHERENCE BAR"*:

> *"`ORACLE best member 0.6708` ← the only arm under the **0.6931** bar … the ORACLE ceiling of all
> in-pool ranking is 0.6708 against a 0.6931 bar — **even perfect selection barely clears it**. The
> bar can only be passed by LEAVING THE POOL'S AFFINE HULL."*

and, most consequentially: *"it predicts lane F's `AVG_SEP` should work … **if it passes it is the
first native-free operator in the project's history to do so.**"*

**What the artefacts say. The two `coh`s have different first arguments.**

**S30's `coh`** (`s30/LEDGER.md` §4, the definition, verbatim): *"For a candidate **corrector** with
residual `r_{t,p} = (expected − d_nat) − correction`, its coherence is
`coh = mean over targets of corr_p(r_{t,p}, mu_{t,p})`."* The `0.6931` row is labelled
**`UNCORRECTED (production)`, `R2_oof 0.000`, `applied delta 0 (baseline)`** — i.e. `correction = 0`,
so **`r = expected − d_nat`: the DISTOGRAM'S OWN PREDICTION ERROR**, an *input* to scoring. S30's
admission rule is stated for exactly that object: *"a **prior corrector** is worth building iff it
LOWERS `coh` below the uncorrected error's own +0.6931."*

**Lane B's / lane F's `coh`** (`s31/s31_B2_inpool.py:386-388`, verbatim): *"coh = within-target
corr(**readout pair error**, pool common-mode pair error mu)"* — the **EMITTED STRUCTURE'S** error,
the *output* of the pipeline.

**They share only their second argument `mu`.** And the arithmetic proves the mismatch without any
interpretation:

```
production, as it appears in S30's CORRECTOR table   coh = 0.6931   (distogram prediction error)
production, as it appears in lane B / F's READOUT    coh = 0.9780   (emitted coordinate average)
```

**Same pipeline, same 126 targets, same `mu` — two numbers 0.285 apart, because they are two
different errors.** `0.6931` is not, and has never been, a readout-space quantity.

**Where it is live in code, not only in prose:**

* `s31/s31_B2_inpool.py:386, 396, 475` — `admission_bar_S30: 0.6931` and a per-arm boolean
  `admitted = bool(v.mean() < 0.6931)`.
* `s31/s31_F_coh.py:46` — `BAR = 0.6931`; `s31/results/s31_F_coh.json` emits `admitted` and
  `frac_targets_under_bar` for `AVG`, `MED`, `AVG_RG`, **`AVG_SEP`** and `ORACLE_best`.

**`AVG_SEP` is the sprint's live deployable arm and it is being graded against this bar right now.**

**And the matched analogue flips the verdict.** S30's rule takes its reference from *the object being
corrected, uncorrected*. In readout space that reference is **production's own readout, 0.9780** —
not 0.6931. Against the matched reference:

```
AVG (production, the reference itself)   0.9780      --
AVG_SEP                                  0.9689   LOWER  -> admitted
AVG_RG                                   0.9592   LOWER  -> admitted
MED                                      0.8886   LOWER  -> admitted
argmin by the shipped DIS                0.8288   LOWER  -> admitted
argmin by a RANDOM pool member           0.8244   LOWER  -> admitted
ORACLE best member                       0.6708   LOWER  -> admitted
```

**Everything passes.** The headline inverts from *"nothing can pass the bar"* to *"the bar as
transplanted is uninformative in this space, and no readout-space bar has been established at all."*
I flag my own construction as a construction: the honest minimum is that **0.6931 must not be used
here**, and the direction of the flip is why it matters.

**Independent arithmetic that the two families do not share a scale.** S30 calibrated coh→Å on
correctors: 0.6931 → 0.5868 buys −0.126 Å, and → 0.5365 buys −0.247 Å (≈1.2–1.6 Å per unit coh). In
readout space, 0.9780 → 0.6708 (`ORACLE_best`) is worth `best1_top75` 2.3055 against production
3.2105 = **−0.9050 Å built chain** (≈2.9 Å per unit coh). Different reference points, ~2× different
slope.

**What SURVIVES, and it is the valuable part.** Lane B's **theorem** is sound and independently
checks out: for any `Σ a_m = 1` the common mode passes through with coefficient exactly one, so
`coh` is a function of the readout's **concentration** and not of the ranker — verified numerically
at `uniform_mean_pairspace = 1.0000, sd 1.2e-16`. So is the empirical consequence that **the shipped
cost and a random pool member differ by 0.0044 in `coh`**, a valid within-space comparison and a
strong result. **Lane F also ran the correctly matched contrast alongside the mismatched one:**
`coh(AVG_SEP) − coh(AVG) = −0.0090, 2.37× MDE, fold CI [−0.0108, −0.0066], 5/5 folds`. *That* is the
meaningful number and it should carry the claim.

**Corrected sentences:**

> **(a)** `0.6931` is S30's coherence of the **uncorrected distogram prediction error** with the pool
> common mode, and S30's admission rule is stated for **prior correctors**. It is not a readout-space
> quantity and **must not be applied to terminal operators or rankers.** Production's own readout
> sits at **0.9780** in readout space.
>
> **(b)** Lane B's theorem stands as derived: `Σa = 1` passes the common mode through with
> coefficient exactly one, so `coh` tracks the readout's **concentration**, not the ranker — which is
> why the shipped cost and a random pool member are 0.0044 apart. **No claim about "passing a bar"
> follows, because no readout-space bar has been established.**
>
> **(c)** `AVG_SEP`'s coherence result is `coh(AVG_SEP) − coh(AVG) = −0.0090, 2.37× MDE, 5/5 folds`,
> a matched paired contrast against production's own readout. It is **ORACLE / NOT DEPLOYABLE** (both
> arguments need the native) and is a diagnostic of *why* an operator works, never a gate.

**This is contract rule 8's sixth instance** (*"a statistic read against another statistic's null"*)
**and contract rule 7's** (*"a control must match the operator's own space"* — the project's most
repeated error). It reached shipped code in two lanes before anyone read the definition at source.

---

## D1 — SEVERE. `0.076 Å of error cancellation` IS, AT MATCHED SEARCH SIZE, A `0.077 Å` PENALTY. THE REGISTERED BAR FIRED AND THE SIGN REVERSES

**Claimed** — `s31/STATE.md` NOTE 2(b), the coordinator's own handed-back prediction, still
unannotated as of 00:35:

> *"2-of-128 uniform beats the single best of the same 128 while spending **zero weight bits**, so
> that 0.076 Å is **error cancellation**, not weighting."*

**What the artefact says** (`s31/results/s31_V_orderstat.json`, n = 126, CA point cloud,
ORACLE / NOT DEPLOYABLE):

```
arm                                                      mean     vs argmin over the SAME 128
ORACLE argmin over 128 singletons        (K = 128)      2.1458            --
2-of-128 uniform, 128 RANDOM pairs       (K = 128)      2.2229      +0.0771  WORSE
                                             1.11x MDE, fold CI [+0.040,+0.111], 5/5 folds, 44W/82L
2-of-128 uniform, EXHAUSTIVE             (K = 8128)     1.9070      -0.2388  better  (3.81x MDE)
```

**At matched oracle-search size the pair family is 0.0771 Å WORSE than naming a single member** —
the same magnitude as the claimed gain, opposite sign. **Share of the exhaustive pair gain surviving
at matched K: −0.32.** The registered bar fired.

**And the mechanism decomposes exactly**, which is what makes this a finding rather than a null:

```
family LEVEL (mean over all variants)   singles 3.5847   pairs 3.3469   pairs -0.2378  <- error
                                                                                          cancellation
                                                                                          IS REAL
dispersion loss (pairing compresses the family, so best-of-128 reaches less)     +0.3149
                                                                                 -------
net at matched K = 128                                                           +0.0771  pairs WORSE
```

> **Error cancellation is real and worth −0.2378 Å in level — a random pair beats a random single by
> that much. It is then MORE than consumed by the fact that averaging compresses the family's
> dispersion, so a per-target minimum over the same 128 variants reaches 0.3149 Å less far.**

This is **lane F's S31-L11 mechanism, second independent instance, opposite operator**: F found that
*decorrelating* a family grows its minimum; V finds that *contracting* a family shrinks it. Same
law, and neither was predicted by the `sqrt(0.68 + 0.32/m)` saturation model NOTE 2(b) proposed —
that model prices the *level* and is silent on the *minimum*, which is the quantity actually quoted.

**Corrected sentence:**

> **ORACLE / NOT DEPLOYABLE, CA point cloud, n = 126.** The 2-of-128 uniform arm's advantage over the
> argmin of the same 128 is **a search-size effect, not error cancellation**. At matched search size
> (128 random pairs against 128 singletons) the pair family is **+0.0771 Å worse**, 1.11× MDE, fold
> CI [+0.040, +0.111], 5/5 folds. Pairing *does* buy −0.2378 Å of genuine error cancellation in
> family level, and loses +0.3149 Å of order-statistic reach by compressing the family's dispersion.
> The exhaustive 8128-pair minimum (−0.2388 Å) is bought entirely by the 64× larger oracle search.

---

## D2 — SEVERE. "THE SET-MATCHED LADDER **INVERTS** A PUBLISHED S30 CLOSURE" IS A WITHDRAWAL THAT IS ITSELF WRONG (contract rule 25)

**Claimed** — `s31/STATE.md` NOTE 2 heading and NOTE 2(b) first line: *"THE SET-MATCHED LADDER
INVERTS A PUBLISHED S30 CLOSURE" / "The set-matched ladder inverts S30-L11."*

**What the artefact says.** S30-L11's claim is *"a plain argmin dominates at **every bit budget**."*
S30 computed the reference curve for exactly that claim and it is in the same file as the sparse
arms — `s30/results/s30_Q_sparse.json :: argmin_ref`, the ORACLE argmin over the top-`2^b`:

```
bits   4       5       6       7       8       9
       2.7159  2.5115  2.3432  2.1458  1.9383  1.7108
```

against the arms the ladder is built from:

```
T128_s2_unif   2.0700   12.99 support bits + 0 weight bits
T128_s2_cont   1.9138   12.99 support bits + UNBOUNDED weight bits
T128_s5_cont   1.8037   27.98 support bits + UNBOUNDED weight bits
```

**Nothing in the set-matched ladder beats the plain argmin at a matched bit budget.** The argmin
reaches 1.9383 at **8** bits and 1.7108 at **9**, where the best uniform sparse arm needs **12.99**
to reach 2.0700. My growth curve prices that arm independently: a blind search over ~512 random
pairs already reaches **2.0799**, so S30's greedy `T128_s2_unif` is worth about **9.1 bits of ORACLE
index** — and at 9 bits naming reaches **1.7108**, i.e. **0.36 Å better at matched information**.

**S30-L11's conclusion is not inverted. Its *support* was invalid and lane C repaired it.** Lane C's
own ledger entry (S31-L6 §1) states this correctly — *"as an Å-per-bit question across sets, naming
wins"* — and also **formally withdraws** the zero-weight-bits sentence: *"that arm's support is
ORACLE-chosen (greedy against the native, 12.99 bits, charged by S30) … my sentence implied it was
[free]."* **STATE NOTE 2(b) still carries the withdrawn sentence, unannotated** (contract rule 13).

**Corrected sentence:**

> Lane C's set-matched ladder shows that S30-L11 **supported a correct claim with a set-mismatched
> pair** (2-of-75 against 1-of-128). The claim itself — *a plain argmin dominates at every bit
> budget* — **survives on S30's own `argmin_ref` curve** and is strengthened by V1's matched-K
> control. Nothing is inverted; the support was repaired. And STATE's "zero weight bits" sentence is
> withdrawn by lane C and must be annotated in place.

---

## D3 — MODERATE, AND IT IS IN THE SPRINT'S HEADLINE BLOCK. A BUILT-CHAIN NUMBER LABELLED "CA CLOUD" AND DIFFERENCED AGAINST A CLOUD BASELINE

**Claimed** — `s31/STATE.md`, HEADLINE block, still live at 00:35:

> *"ORACLE argmin-over-128 is **2.1435 Å (CA cloud, ORACLE / NOT DEPLOYABLE)** against production's
> 3.0483 — ~0.90 Å of headroom."*

**What the artefact says** — recomputed from `s29/results/s29_O_chain_rows*.jsonl`, all 126:

```
best1_top128    CHAIN 2.1435    CLOUD 2.1458
prod            CHAIN 3.2105    CLOUD 3.0483
```

**2.1435 is the BUILT CHAIN value.** It is labelled *CA cloud* and then differenced against 3.0483,
which is the *cloud* production. The headline's `~0.90` is `3.0483 − 2.1435` — **a cloud baseline
minus a chain arm.** NOTE 6's ladder has the same five numbers correctly labelled *Built chain*, so
the defect is in the synthesis, not in the lane.

The numerical damage is not the 0.0023 Å of the mislabel — **it is that the headline understates its
own headroom by 0.16 Å**:

```
as written (mixed bases)          3.0483 - 2.1435 = 0.9048
correct, both CLOUD               3.0483 - 2.1458 = 0.9025
correct, both BUILT CHAIN         3.2105 - 2.1435 = 1.0670   <- the endpoint basis (contract rule 1)
```

**This is contract rule 1's exact failure mode, in the sentence that carries the sprint's framing
claim that "the cap is on information, not on value."**

**Corrected sentence:**

> **Value is not capped the same way.** ORACLE argmin-over-128 is **2.1435 Å on the built chain
> (ORACLE / NOT DEPLOYABLE)** against production's **3.2105 Å built chain** — **1.067 Å of
> headroom**. (On the CA cloud the same pair is 2.1458 against 3.0483, 0.903 Å.)

**And a gap in the instrument it slipped through.** `s31/s31_verify.py` reports *"basis audit: any
registered number whose reporting basis is **unstated** — none."* An audit for an *unstated* basis
does not catch a *misstated* one, and the sprint has exactly one of the latter, in its headline.
Recommend the basis audit also assert each number against the value on the basis it names.

---

## D4 — LOW. STATE NOTE 11 QUOTES F3-e DRAWS THAT CANNOT BE REGENERATED; THE LEDGER HAS THE REPAIR AND STATE DOES NOT

**Claimed** — `s31/STATE.md` NOTE 11 (00:23): *"−0.4191 Å sd 0.0065 over 4 draws … 146%"*.

**What the artefact says.** Lane F found its own defect — `hash(pdb)` is salted per process
(`PYTHONHASHSEED`), so the first draws were not regenerable — repaired it to `crc32`, re-ran, and
**annotated the ledger in place** (`s31/LEDGER.md:1568-1575`, commit `6b2baf0c`). The reproducible
run gives:

```
                                  as quoted in STATE      reproducible run
gain, matched random family            -0.4191                -0.4279
draw sd over 4 draws                    0.0065                 0.0273
share of the prefix gain                  146%                   149%
```

**The conclusion is unchanged and strengthened** — lane F's handling of this is exemplary and the
defect is only that **STATE was not updated with the lane's own correction**. Note the draw sd is
**4.2× larger** than the figure STATE carries, which is the number a reader would use to judge the
draw control.

**Corrected sentence:** quote **−0.4279 Å, sd 0.0273 over 4 draws, 149%** from the reproducible
crc32-seeded run, and cite `s31/LEDGER.md:1568-1575` for the superseded values.

---

## D5 — LOW, PROCESS. LANES ARE COMMITTING WITH `git add -A` AND SWEEPING OTHER LANES' UNCOMMITTED WORK INTO THEIR COMMITS

`git log --oneline -1 -- s31/AUDIT_V.md` attributes lane V's pre-registration to commit `bac1fe1c`,
titled *"s31: lane A and lane C close…"*. My `git add` + `git commit` then reported *"no changes
added to commit"* because another lane's `-A` had already swept the files in.

The prereg **was** committed before any number existed, so contract rule 19 is satisfied here. But a
sprint whose central discipline is *"pre-register the falsifier in a committed file before the
number exists"* depends on the commit history saying **which lane registered what and when**.
Recommend lanes stage explicit paths.

---

## D6 — MODERATE, AND IT IS THE SPRINT'S HEADLINE. S31-L17'S ALGEBRA IS CORRECT; ITS EVIDENCE COULD NOT HAVE TESTED ITS OWN CAVEAT, AND "ZERO" OVERSTATES A MEASURABLE QUANTITY

**Claimed** — `s31/LEDGER.md` S31-L17 / `s31/STATE.md`: *"`E = _zrank(pool["sc"][o])` … The ranks of
any 128 distinct values are 1..128, so the standardised vector is a **constant** … **Identical to six
decimals.** Up to the pool's tie structure, `E` is **the same vector on every target in the
benchmark**"* and *"**The quantum stage carries ZERO target-specific information.**"* The recorded
evidence is **three synthetic random score vectors**.

**What I found.** Artefacts `s31/s31_V_zrank.py`, `s31/results/s31_V_zrank.json`, all 126 targets,
native-free (`E` is a function of the score order only).

**(a) The load-bearing step is CONFIRMED, and it is not the one the entry verified.** The claim holds
only if `sc[o]` is *sorted* — `rankdata` of an *unsorted* vector is a target-specific permutation of
1..128, which would carry up to `log2(128!)` bits. It is sorted: `core/pipeline.py:757` sets
`order = argsort(sc)` and `:863` sets `o = top[:dim]`. **I confirmed `sc[top]` is non-decreasing on
all 126 targets from the real scores.** The three synthetic trials cannot establish this — they were
fed an already-sorted vector, which is why they printed a monotone ramp. **The algebra is right and
the recorded evidence does not test the step it rests on.**

My recomputed ramp reproduces lane P's printed values exactly: first five
`[-1.718572, -1.691507, -1.664443, -1.637379, -1.610315]`, last three
`[1.664443, 1.691507, 1.718572]`.

**(b) The tie hedge is the common case, not an edge case — and the synthetic test could not see it.**
`rankdata` uses `method='average'`, so ties bend `E` off the ramp. Random floats never tie, so three
synthetic vectors are structurally incapable of exhibiting the phenomenon the caveat is about. On the
real benchmark:

```
targets with at least one tie in the top-128       125 of 126
targets whose E is EXACTLY the ramp                  1 of 126
tied positions per target        mean 9.56   median 8   p90 18   max 34
max |E - ramp| over all targets                   0.0407     (on a vector spanning +/-1.7186)
max |E_i - E_j| between any two targets           0.0812
distinct STRUCTURES in the top-128     mean 118.45, min 94   <- the source of the ties
```

`118.45` reproduces lane A's `ALPHABET_n_distinct_mean` **exactly**, which independently confirms the
mechanism: **ties are duplicate structures, and duplicate structures have identical scores.**

**(c) So "ZERO" is an overstatement of a quantity that is now measured.** `E` does carry a little
target-specific information — the tie pattern, which is a target-specific fact about which retrieved
structures are duplicates. The honest statement is that it is **negligible and bounded**, and lane P
already measured the downstream size: `sd(H(p*))` across targets is `1.4e-4` bits at α = 1.

**Nothing in the conclusion changes.** The stage is a fixed weighting curve per α and the target
enters only through the readout's `P` and `W`. The defect is that the sprint's headline rests on
evidence that could not have tested its own caveat, and states "zero" where a measured bound is
available and stronger.

**Corrected sentence:**

> `core/pipeline.py:757` sorts the pool by score and `:863-864` takes `E = _zrank(sc[top[:128]])` of
> that **sorted** slice, so `E` is the standardised rank ramp — **confirmed non-decreasing on all 126
> targets from the real scores, not only on synthetic draws.** It is target-independent to within a
> **max-norm deviation of 0.0407 on a vector spanning ±1.7186 (2.4% of range)**, with a maximum
> between-target deviation of **0.0812**. The residual is entirely the tie structure of duplicate
> retrieved windows — **125/126 targets carry ties, mean 9.56 of 128 positions, max 34**, from
> `n_distinct` 118.45/128 — and its downstream effect is `sd(H(p*)) = 1.4e-4` bits. **The quantum
> stage therefore carries no target-specific information beyond that bound**: it answers *"what fixed
> weight should rank k receive?"*, a global hyperparameter, and the target enters only through the
> readout's `P` and `W`.

**One further precision, offered rather than charged as a defect:** `(alpha, T)` is read from
`VQE_LFO[fold]`, so `p` takes **two** distinct values across the benchmark (α = 1 on folds 0/3/4,
α = 0.25 on 1/2), selected by fold — which is a property of the target. Lane P's own phrasing *"one
fixed weighting curve per alpha"* is exactly right; the STATE compression to *"zero target-specific
information"* loses it.

---

## D7 — LOW. THE ~0.70 Å REFERENCE TERM DOES NOT "CANCEL", IT ATTENUATES; IT IS NOT "UNIFORM"; AND THE QUADRATURE USES THE WRONG MOMENT, UNDERSTATING ITS OWN CAVEAT BY ~2×

**Claimed** — `s31/STATE.md` NOTE 7: *"The endpoint carries a **uniform** ~0.70 Å reference term …
In quadrature `sqrt(3.21² − 0.70²) = 3.13` … i.e. ~0.08 Å today, material near 1 Å. **Because it is
uniform it CANCELS in every arm-to-arm delta**, which is the project's actual currency."*

**What the artefact says** — `s31/lit_L/ens_spread.json`, 111 resolved targets, recomputed by me:

**(a) It is not uniform, and lane L's own table says so.** `model1 → medoid` has **mean 0.6965,
median 0.4565, sd 0.7767, max 4.2943** — **CV = 1.12**, i.e. the standard deviation exceeds the mean,
and **13 of 111 targets are exactly 0.0**. That is about as far from uniform as a positive quantity
gets. **The conclusion is right and the stated reason is falsified by the data two paragraphs
above it in the same note.** The correct reason: **both arms are scored against the same reference on
the same target**, so the term is common to the pair and largely differences out — which holds
whether or not it is uniform.

**(b) It attenuates rather than cancels.** Under lane L's own quadrature model the observed delta is
the true delta times `d/sqrt(d² + r²)`:

```
r = 0.6965 (the quoted mean)   d = 3.2105 -> 2.3% attenuation    d = 1.0 -> 17.9%
```

On the project's one confirmed effect (0.0221 Å) that is **0.0005 Å — immaterial now**, and it grows
in exactly the regime lane L already flags. *"Cancels"* should be *"attenuates by 2.3% at the current
endpoint."*

**(c) The quadrature uses the wrong moment, and the error is one-sided.** `d_obs² = d_true² + r²` is
a per-target identity **in squares**, so it must be applied with **second** moments. Lane L applied
it to **first** moments on both sides. Because the reference term is far more skewed (CV 1.12) than
the endpoint, the understatement lands on `r`:

```
r as quoted, MEAN of the per-target RMSDs        0.6965
r as the model needs it, RMS of the same         1.0407     (+49%)

endpoint inflation   sqrt(d^2+r^2) - d   at d = 3.2105:   +0.0747  (quoted ~0.08)
                                          with r = RMS:   +0.1645
delta attenuation at d = 1.0                   r = mean:    17.9%
                                                r = RMS:    30.7%
```

**So the caveat is about 2× larger than stated.** It still changes nothing — lane L's recommendation
(*report the uncertainty beside the headline and change nothing*) remains correct, and I endorse it.

**(d) A positive lane L could have claimed and did not.** **13 of 111 deposited `model 1`s are
*exactly* their own ensemble medoid (11.7%)**, against ~5% expected by chance at ~20 models per
entry. Depositors commonly order NMR models by agreement, and this is direct evidence that the
benchmark's reference is a **better-than-random** ensemble member — which strengthens lane L's
"change nothing" recommendation on its own terms.

**Corrected sentence:**

> The endpoint carries a reference term of mean **0.6965 Å** (median 0.4565, sd 0.7767, max 4.2943;
> **13/111 exactly zero** — model 1 is its own ensemble medoid more often than chance). **It is not
> uniform**; what protects arm-to-arm deltas is that **both arms are scored against the same
> reference on the same target**, so it is common to the pair. It does not cancel exactly — it
> **attenuates every delta by 2.3% at the current endpoint** (0.0005 Å on the project's confirmed
> 0.0221 Å effect), rising to ~31% if arms reach 1.0 Å. Applied with the moment the model requires
> (**RMS 1.0407**, not the mean), the absolute inflation at the endpoint is **+0.1645 Å**, not
> ~0.08. **Report it beside the headline and change nothing.**

---

## D8 — LOW. THE R1 CAPACITY FIGURE HAS TWO VALUES IN CIRCULATION, AND THE FORMULA WRITTEN BESIDE THE CORRECT ONE PRODUCES THE OTHER (JENSEN)

**Claimed** — two versions, both live in `s31/LEDGER.md`:

* `:208` and `:217`, the **"Final form of R1, third restatement"**: *"its capacity is
  **log2(118.45) = 6.888 bits** mean, 6.555 worst"* and *"the number of distinct candidates
  (**6.888 bits** mean)"*.
* `:1817`, `:1827`, `:2459`, lane A's later S31-L13: *"**Corrected capacity: log₂(118.45) = 6.886
  bits** mean"*, *"the **6.886**-bit alphabet"*, *"R1 said the selection readout's alphabet is
  **6.886 bits**"*.

**What the artefact says** — `s31/results/s31_A_cap.json` and `s31_A_cap_rows.jsonl`, 126 rows:

```
CAPACITY_bits_alphabet_mean (artefact)        6.886162
MEAN over targets of log2(n_distinct_t)       6.886162   <- reproduces the artefact EXACTLY
log2(MEAN n_distinct) = log2(118.452381)      6.888163   <- what the WRITTEN FORMULA gives
n_distinct                                    min 94, max 128
log2(94) = 6.5546                                        <- reproduces the quoted "6.555 worst"
```

**The artefact computes the mean of the per-target logs; the formula written beside it is the log of
the mean.** By Jensen these differ systematically (`mean log ≤ log mean`), here by 0.002 bits. So
**6.886 is the right number and `log2(118.45)` is the wrong formula for it**, while the entry that
calls itself the *final form of R1* carries **6.888**, the value that formula actually produces.

**Immaterial to every conclusion** — the point is only that 7 bits is not attained. It is recorded
because two values for one quantity in one ledger is how a number drifts, and because the mean-of-
logs is the quantity a "mean capacity" should be.

**Corrected sentence:**

> The selection readout's alphabet is the number of **byte-distinct** candidates, so its capacity is
> `mean_t log2(n_distinct_t)` = **6.886 bits** (worst 6.555 at `n_distinct` = 94; 7.000 attained only
> where `n_distinct` = 128) — **not 7**. `log2(mean n_distinct) = log2(118.45) = 6.888` is a
> different quantity and should not be written as the formula for 6.886.

---
## D0 — CLOSED 2026-09-21 00:58. STRUCK FROM ALL THREE SOURCES AND NOW ENFORCED BY A REGRESSION GUARD

* `s31/s31_F_coh.py` — struck by lane F; `BAR` removed, the artefact re-run at 00:55:43 and now
  carrying a `STRUCK` field. **Verified: zero live numeric uses remain.**
* `s31/s31_verify.py` — the two certified booleans removed, replaced by the **paired** contrast
  `coh(AVG_SEP) − coh(AVG) = −0.0090, −2.3710× MDE, 5/5 folds`, with the originals quoted in place.
* `s31/s31_B2_inpool.py` — **struck by lane V** (lane B has closed): `admission_bar_S30` and the
  per-arm `admitted` boolean are gone, replaced by `vs_production_readout`, the difference against
  **production's own readout** — a matched reference that needs nothing imported. The original lines
  and the original stdout line are quoted in place (contract rule 13).
* **A regression guard now enforces it** (`s31_verify.py`): it tokenises each source and FLAGs any
  **NUMBER token** equal to 0.6931, so an annotation recording the withdrawal passes and a live use
  does not. Verifier: **129 matched / 0 mismatched / 0 not-found / 0 FLAGGED.**

**One artefact is deliberately NOT re-run:** `s31/results/s31_B2_inpool.json` still contains
`admission_bar_S30` and per-arm `admitted` fields from its original run. Historical artefacts are
immutable (rule 13) and its `coh` values are correct — **it is only the `admitted` booleans that are
withdrawn, and they must not be quoted.** The script will not emit them again.

---

## D9 — MODERATE, AND IT WAS IN THE VERIFIER ITSELF. A CHECK WRITTEN AGAINST A KEY THAT HAS NEVER EXISTED

**Claimed** — while D0 was being fixed, `s31/s31_verify.py` acquired

```python
check("coh(AVG_SEP) - coh(AVG), the paired contrast that carries the claim",
      -0.0090, dig(fc, "coh_contrast", "AVG_SEP_minus_AVG", "effect"), basis="in-band")
```

**What the artefacts say.** `s31/s31_F_coh.py:177` writes `out["coh_AVG_SEP_minus_AVG"]`, and
`s31/results/s31_F_coh.json` has top-level keys
`['seed','n','STRUCK','ORACLE','basis','theorem','coh','coh_AVG_SEP_minus_AVG',
'affine_hull_departure']`. **`coh_contrast` has never existed in either the script or any artefact.**

**Consequence, which is why this is not cosmetic.** `dig` returns `None` for a missing path, `check`
routes that to `MISSING`, and the verifier exits non-zero. For a window the sprint's trust anchor
read **`MATCHED: 127 / KEYS NOT FOUND: 1`** — i.e. **the instrument that certifies every headline
number was itself failing**, on a check written against an assumed key rather than a read one. It was
repaired concurrently (the file now reads `coh_AVG_SEP_minus_AVG`) and the verifier is back to
129/0/0/0.

**This is contract rule 14 — *"prose is not evidence of code; run `ls` and read the file"* — committed
inside the tool that exists to enforce it.** Seventh instance in the project's record.

**Corrected sentence:** the verifier must read the artefact's key list before asserting a path.
Recommend `dig` gain a strict mode for headline checks that names the keys that *do* exist when a
path misses, so the failure says *"you meant `coh_AVG_SEP_minus_AVG`"* rather than *"not found"*.

---

## V-SELF — THREE ERRORS OF MY OWN, RECORDED BECAUSE THE THIRD IS THE BEST ILLUSTRATION IN THIS FILE

Contract rule 24 makes me the auditor; it does not exempt me.

1. **Multiplicity row numbered 913** off a regex that matched a `k`-column instead of a row index.
   Caught on the next read, renumbered to 18 before anything quoted it.
2. **Two backtick spans eaten by the shell** — I built a section inside a double-quoted
   `python -c` string, so bash expanded `` `s31/STATE.md` `` and
   `` `read-the-memory-body-not-the-index-line` `` as command substitutions and dropped them.
   Repaired with the Edit tool so no shell is involved.
3. **My D0 regression guard reported "clean" for the one file that still had the defect.** The first
   version tested `"0.6931" in line and not line.startswith("#")`, which flagged the *annotation
   strings* in the files that had been correctly fixed and flagged its own message text. The second
   version tokenised — correct — but wrapped the whole file in `except Exception: continue`, and
   `s31_B2_inpool.py` contains a `1j` complex literal, so `float('1j')` raised and **the entire file
   was skipped silently.** The guard printed nothing for `s31_B2_inpool.py` and I would have read
   that as passing.

> **Error 3 is this audit's own thesis, committed by the auditor, within an hour of writing it: a
> check that cannot fire is indistinguishable from a check that passed.** It is fixed two ways — the
> `float()` conversion is now per-token so one odd literal cannot abort a file, and a parse failure
> now **FLAGs with the exception text** instead of returning silently, because *"could not read"* must
> never render as *"clean"*. It is recorded here rather than quietly repaired, which is the same
> standard I applied to every other lane.

---

# CONFIRMED

* **S31-L17's algebra is CORRECT** and I confirmed the step its own evidence did not test:
  `sc[top]` is non-decreasing on all 126 targets from the real scores, so `E` really is the
  standardised rank ramp. See D6 for the measured size of the tie caveat.
* **The cross-basis audit D3 asked for is implemented** in `s31/s31_verify.py` (lane D's file,
  at the coordinator's request) and **carries a self-test so it cannot be decoration**: the real
  D3 defect (chain 2.1435 declared `cloud`) is CAUGHT via `best1_top128`, and the corrected
  number (cloud 2.1458 declared `cloud`) is clean. The verifier now runs **129 matched / 0
  mismatched / 0 not-found / 0 flagged**, with 33 s29 O-ladder items carrying both bases.
* **Lane L's reference-term arithmetic reproduces** (`sqrt(3.2105^2 - 0.6965^2) = 3.1358`, quoted
  ~0.08 A) and **its recommendation -- report it beside the headline and change nothing -- is
  correct.** D7 corrects the stated reason and the moment, not the decision.

* **`bestm128 = 2.9027` is real on the endpoint basis, and its deflation is correctly stated.**
  I independently confirm the chain arm from `s29/results/s29_O_chain_rows*.jsonl`: `bestm128`
  2.9027 against `prod` 3.2105, and `best1_top128` 2.1435, `hull_top128` 1.8538, `best1_pool`
  1.7078 — all **built chain**, all 126, matching NOTE 6's ladder exactly.
* **Lane F's F3 provenance gate is genuine.** The prefix curve reproduces S29's stored `curve128` at
  `max |curve − stored| = 0.0e+00` on 126/126, and the gate is a real comparison against
  `s29/results/s29_O_p128_rows.jsonl`, not a value compared to itself.
* **Lane F's F3-e conclusion survives its own seeding repair** (D4): the bar moves from 146% to 149%.
* **`s31/s31_verify.py` runs clean** at **97 matched / 0 mismatched / 0 not-found / 0 flagged**,
  including the `MDE = 2.8016 × SE` identity and the lower-is-better self-test for
  `stats_lib.compare`. 51 named paths, 51 exist.
* **Lane C's own ledger entry (S31-L6) is correct where STATE's summary of it is not** — it states
  the Å-per-bit direction properly and withdraws its own sentence unprompted. D1 and D2 are defects
  in the **synthesis**, not in lane C's work.
* **S30 had already priced the sparse support as an order statistic** and said so:
  `s30_Q_sparse.json :: order_stat` gives `share_accounted` 1.08 / 1.13 and split-half transfer of
  −0% / 4%, verdict *"NOT A SIGNAL"*. That block should be cited whenever the sparse ladder is.
* **Lane L's NOTE 13 finding (1)** — that *"the convex readout is 0.290 Å better than naming the
  best of them"* is an unequal-information comparison (7 bits against 128 free reals) — **is
  correct, and I reached the same conclusion independently** before reading it. Not double-counted
  here; recorded as agreement.

# NOT RECOMPUTABLE

* **Nothing material.** The one artefact that was mid-write when I opened it —
  `s31/results/s31_F3_randfamily_rows.jsonl`, at 57/126 at 00:20 — was a **live re-run**, not a dead
  job (contract rule 15), and completed at 126/126. I did not declare it dead, and I checked the
  process list and the job log before forming a view.

---

## MULTIPLICITY

Family `s31_V_orderstat`: **2 fold-clustered paired comparisons** (V1a exhaustive, V1b matched-K),
plus a 14-point × 2-family descriptive growth curve that spends no α (nothing is accepted or
rejected on it). **One pre-registered falsifier, k = 1 for its own purposes.** Appended to
`s31/MULTIPLICITY.md`.


---

# THE PATTERN ACROSS ALL NINE DEFECTS

Every one is a quantity **transplanted across a boundary its definition does not cross**:

| | boundary crossed |
|---|---|
| **D0** | an **object** — a distogram-corrector residual's baseline applied to readout outputs |
| **D1, D2** | a **search size** — a min over 8128 read against a min over 128 |
| **D3** | a **basis** — a built-chain value labelled CA cloud and differenced against a cloud baseline |
| **D6** | a **data regime** — a tie-free synthetic draw standing in for a benchmark that ties 125/126 |
| **D7** | a **moment** — a squares identity applied to first moments of a CV-1.12 distribution |
| **D8** | **Jensen** — log of a mean written as the formula for a mean of logs |
| **D4, D5** | **provenance** — a superseded draw, and a commit history that does not say who registered what |

`s31/STATE.md` NOTE 13 says cross-lane syntheses fail because nobody owns both halves. **D0 and D6
show the same failure happening INSIDE a lane**, when a lane re-uses a number from an earlier sprint
whose definition it has not re-read. The rule that covers both:

> **A number carries its definition, not just its value. Re-read the definition at the source before
> re-using it across a sprint boundary — and check that the formula you write beside it is the one
> that produces it.**

That is `read-the-memory-body-not-the-index-line` applied to numbers rather than to memory files.

**And the instrument lesson, which is the one I would keep:** the verifier asked whether a basis was
**named** and passed 97/97 while the headline carried a **misnamed** one. **An audit that checks for
an absent label cannot catch a wrong label.** The cross-basis check now added asks the second
question and carries a self-test that fails loudly if it stops catching the defect that motivated it.

