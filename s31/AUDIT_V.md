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

# CONFIRMED

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
