# S32 MULTIPLICITY LEDGER

**Owner: LANE V.** Charter §45, contract rule 9. Every emitted comparison is written here **as
it is emitted**, not at the end. A per-target minimum over K variants is mostly an order
statistic — price best-of-K before calling it a lead, and prefer split-half transfer
(`s24/stats_lib.py::best_of_k_within`, `::split_half_transfer`).

Status legend: **REG** pre-registered · **EXP** exploratory · **CONF** confirmatory
(independent replication of an already-emitted arm) · **AUDIT** verification, no hypothesis.

Rules that bind every row:
- `MDE = 2.8016 x SE`, per comparison. `< 0.7x` NOT A RESULT · `0.7-1.0x` NOT MEASURED ·
  `>= 1.0x` + fold CI excluding zero + >= 4/5 folds agreeing = RESULT.
- Every row names its **basis** (chain / cloud / set-mean / member / selection / in-band).
- A chain contrast must be projected **in the same job from bit-identical clouds** — see D1.

---

## Running totals

### The sprint-level fact, walked over every artefact rather than recalled

`s32/results/*.json`, every `effect_over_mde` in every nested block:

| | count |
|---|---|
| improvements clearing 1.0x MDE on the **BUILT CHAIN, deployable** | **0** |
| improvements clearing 1.0x MDE on **any** basis | 30 |
| … of which ORACLE (read the native to select) | ~20 |
| … of which filter-vs-random on the **pool mean** (a diagnostic) | 4 |
| … of which **in-band Spearman** diagnostics | 6 |

**Not one deployable built-chain improvement anywhere in the sprint.** Every arm that clears MDE
either reads the native, or is measured on a basis that is not the endpoint, or points the wrong
way. The three largest ORACLE ceilings: `ORACLE top-5 through the deployed operator` **−1.567**
(5.17x, cloud), `ORACLE per-target scale` **−0.1352** (3.84x, chain), `ORACLE best branch`
**−0.1072** (2.93x, chain, split-half transfer **3%**).

---

*(Lane-V rows only; the per-lane sections below carry each lane's own count.)*

| | count |
|---|---|
| audit arms (AUDIT, no hypothesis, no multiplicity cost) | 10 (V-A1 … V-A10) |
| confirmatory / independent replication arms (CONF) | 2 (V-A7 `cos`, V-A8 DIS row of D1-T) |
| defects found and reported when found | 7 (D1 … D7) |
| **lane-V comparisons emitted against the endpoint** | **0** |

**Lane V emits no endpoint comparison by design.** Its arms are audits and replications, which
carry no multiplicity cost against the sprint's discovery budget — but the arms they audit do,
and two are logged here explicitly:

- **V-A5** re-gated 2 comparisons already emitted by S32-L3 (it did not add new ones).
- **V-A8** emitted **8** in-band comparisons (4 scorers x {full band, deduplicated}). All are
  DIAGNOSTIC, on the in-band Spearman basis, **never differenced against a chain RMSD**, and one
  of the 8 is a pure-noise falsifier included precisely so the family can fail. Reported in full,
  not best-of.

---

## Emitted comparisons

| # | lane | arm | kind | basis | n | effect | SE | MDE | x MDE | folds | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| — | — | *(none yet)* | — | — | — | — | — | — | — | — | — |

---

## Audit arms (LANE V) — no hypothesis, no multiplicity cost

| id | what | artefact | outcome |
|---|---|---|---|
| V-A1 | Charter Step 4: rebuild the endpoint from artefacts (targets -> pool -> score -> top-75 -> coordinate average -> projection -> chain) | `s32/results/s32_V_step4_endpoint.json` | cloud / set-mean / pool / top-75 confirmed exactly; chain to 0.0021 only — see **D1** |
| V-A2 | Is 3.2105 bit-reproducible from its stated input? | `s32/results/s32_V_chain_bitexact.json` | **YES, 126/126 bit-identical** from `s29_O_structs['prod']`; the two other float64 representations of the same average give 0/126 — see **D1** |
| V-A3 | Size-matched random null for the ladder's narrowing increments (2000 draws/target) | `s32/results/s32_V_ladder_orderstat.json` | see **D2** |
| V-A4 | rr spread of the DIS top-128 vs the pool (mechanism for V-A3) | `s32/results/s32_V_top128_spread.json` | mean 4.4533 -> 3.5847, sd 1.3206 -> 0.6384, **p5 2.6098 -> 2.6184 (unchanged)** |
| V-A5 | The full contract-rule-1 gate + rule-10 draw distribution + rule-12 strata on V-A3 | `s32/results/s32_V_orderstat_gate.json`, `s32_V_orderstat_strata.json` | see **D3** — the aggregate is an outcome-defined stratum |
| V-A6 | The ARITHMETIC-NOISE distribution of the endpoint: eps = 1e-14 A per coordinate, 5 draws x 126 | `s32/results/s32_V_ulp_distribution.json` | **draw-to-draw sd 0.0030 A**; per-target &#124;delta&#124; mean 0.0127 / p90 0.0248 / max 0.5006; cloud unchanged to 1.1e-13. *The irreducible sd on an unpaired built-chain mean.* |
| V-A7 | Is lane R's `cos_align` independent evidence, and is it the cosine it is named for? | `s32/results/s32_V_cos_identity.json` | **NOT independent** — price rebuilt from (e,d,cos) to max error **0.000e+00**, a bijection. **But it IS a real cosine**: direct measurement in a common frame agrees to 0.0038 mean / 0.0642 max. Production mean cos **-0.052**, on the wrong side of the orthogonal null. |
| V-A8 | Adversarial replication of lane D's D1-T sign transfer, with the three controls its artefact lacks | `s32/results/s32_V_D_signadversary.json` | see **D4** |
| V-A9 | **Is the ladder's top rung retrieval, or fragment-space capacity?** Three size-matched ORACLE hull arms | `s32/results/s32_V_hull_capacity.json` | BLOSUM500 **1.1167** (reproduces s29 `hull_pool` to 4 dp); RAND500 **1.1495** (+0.0328, 0.38x, 63W/63L); DONOR500, *another target's fragments*, **1.1626** (+0.0459, 0.49x) — **neither a result.** Retrieval's share of the 2.10 A top rung is **2.2%**. |
| V-A10 | BLOSUM top-75 vs the deployed SCORE top-75 on the best member, stratified | `s32/results/s32_V_prefix_strata.json` | −0.2020 overall (1.14x, RESULT) but **−1.5871 on FAIL18, 18W/0L**, and **+0.0288 on the other 108 (0.26x, NOT A RESULT)**; median over 126 **exactly 0.0000** — the third instance of D3 |

---

## Defects found (reported to the coordinator when found, not at the end)

### D1 — the canonical endpoint is a re-projection, and the chain rung is not continuous in its input

`3.2105` is **not** any cached scalar. The production cache holds four distinct chain-like
means over the same 126 targets:

| value | object |
|---|---|
| 3.048338 | CA point cloud (`avg_ca`), unprojected |
| 3.204076 | lambda=0 projection arm (`fit_ca`) — no Ramachandran penalty |
| 3.214765 | lambda=0.3 arm (`ca`) — **the production pipeline's own emitted chain** |
| 3.235460 | lambda=0.3 arm after AMBER relaxation (`amber_ca`) |

`3.210534` is the mean of the 126 `item="prod"` rows of `s29/results/s29_O_chain_rows*.jsonl`,
written by `s29/s29_O_ladder.py::chain_item`, which **re-projects the stored cloud** through
`s12/instrument.project`. It is 0.0043 A better than the chain production itself emits.

**The mechanism.** The projection is perfectly deterministic given bit-identical input (3
in-process repeats identical to 9 dp; BLAS thread count 1/2/4/8 makes no difference) and is
**discontinuous** in it:

```
2LNG   max|C_recomputed - C_stored| = 7.1e-15   (cloud RMSD identical: 4.249160348 both)
         stored  cloud -> chain 4.778535407
         recomputed   -> chain 4.879180560      0.1006 A apart
6QAX   max|dC| = 3.6e-15    4.226988500 vs 4.078647277      0.1483 A apart
```

Contract rule 3 already prices the amplification at ~1e13; 7e-15 x 1e13 = 0.07-0.15 A, which is
what is measured. **Consequence for every lane: a chain contrast is only meaningful if both
sides are projected in the same job from bit-identical (not merely equal-valued) clouds.**
Unpaired cross-job chain claims below ~0.03 A are not resolvable.

### D2 — the ladder's narrowing increments are unpriced order statistics, and one is mislabelled

The quoted ladder writes `+ retrieval filter K=500 -> 128  +0.4357` and
`+ prefix 128 -> 75  +0.1620`.

(a) **Mislabelled.** Retrieval is the BLOSUM prefix (universe -> K=500). The 128 is the
**distogram Bayes-risk score prefix**, which `s29_O_ladder.py` itself calls *"rung 9: the
top-128 prefix (the quantum field of view)"* — 2^7, the VQE register width. The loss belongs to
the score stage, not the retrieval stage.

(b) **Unpriced.** Both increments are differences of per-target minima over nested sets
(`overlap_top75_in_128 = 1.0` on all 126, so nesting holds and the increments are non-negative
by construction). Size-matched random null, 2000 draws/target, member basis (projection price
-0.0007..-0.0030 A for this object class, contract rule 16):

```
K=500 -> 128 TOTAL                        +0.4350
  size effect (uniform random 128 of 500) +0.2477     57% is bare set size
  ordering effect (DIS 128 - random 128)  +0.1872     [see D3: mean only]
128 -> 75 TOTAL                           +0.1604
  size effect (random 75 of 128)          +0.0907
  ordering effect (DIS 75 - random 75)    +0.0696     [see D3: mean only]
```

Mechanism (V-A4): the score improves the pool's mean rr 4.4533 -> 3.5847 and halves the spread
(sd 1.3206 -> 0.6384) while leaving the 5th percentile **unchanged** (2.6098 -> 2.6184). It
concentrates on the mode and buys nothing in the good tail.

> **RETRACTED IN PART, 2026-09-21, by D3 below, annotated in place per contract rule 13.**
> The original wording of this paragraph read *"The deployed score's ordering is **worse than
> chance** at retaining oracle headroom"* and the closing sentence read *"...and 0.257 is the
> score ordering performing worse than random truncation."* **Both are true of the MEAN only.**
> The median runs the other way (-0.0380), the score wins on 72 of 126 targets, and the entire
> mean is the 18 outcome-defined FAIL18 targets. See D3.

Corrected sentence: *"narrowing 500 -> 128 -> 75 costs 0.595 A of oracle-best headroom, of
which 0.338 is the set-size order statistic; the remaining 0.257 is 18 targets' worth of the
score placing its window in the wrong part of the pool, and is NOT MEASURED on the other 108."*

### D3 — the ordering effect is entirely FAIL18, and FAIL18 is defined by it

`s32/results/s32_V_orderstat_gate.json`, `s32/results/s32_V_orderstat_strata.json`,
`s32/s32_V_orderstat_gate.py` (regenerates the 2000 draws from the pinned seed and asserts it
reproduces `s32_V_ladder_orderstat.json` per-target to 0.00e+00 before computing anything).

The full contract-rule-1 gate, of which S32-L3 quoted one third:

```
SCORE top-128 vs RANDOM 128 of 500
  effect +0.1872  SE 0.0630  MDE 0.1764  1.06x   MEDIAN -0.0380   72W/54L
  folds {0: -0.1087, 1: +0.2621, 2: +0.2938, 3: +0.2566, 4: +0.2344}
  fold CI95 [+0.0383, +0.2743] excludes zero YES; folds same sign 4/5 YES; type_m 1.10 (flag)
  >>> GATE: RESULT
SCORE top-75 vs RANDOM 75 of 128
  effect +0.0696  SE 0.0243  MDE 0.0681  1.02x   MEDIAN -0.0045   65W/61L   >>> GATE: RESULT
```

`stats_lib.compare` is LOWER IS BETTER, so `n_better = 72` means **the score is better on 72 of
126 targets** while the mean runs against it. Mean +0.187 / median -0.038 with a near-even W/L
is the median-vs-mean gap the project calls the free early warning, and it fired.

Contract rule 10, on the aggregate rather than within target: the 126-target draw mean has sd
0.0232 (500->128) and 0.0167 (128->75); **only 69.5% and 53.6% of SINGLE random draws clear
their own comparison's MDE.** The second is a coin flip.

Contract rule 12, which decides it:

```
                 ALL 126            FAIL18 (n=18)      OTHER 108
500 -> 128   +0.1872 (med -0.0380)  +1.4879  0W/18L  -0.0296 (med -0.0743, 72W/36L)  -0.30x NOT A RESULT
128 ->  75   +0.0696 (med -0.0045)  +0.3611  3W/15L  +0.0210 (med -0.0090, 62W/46L)  +0.38x NOT A RESULT
```

The ten worst targets are all in FAIL18. Drop the 10 worst and the aggregate falls +0.1872 ->
+0.0294; drop 20 and it is **negative**, -0.0567. **FAIL18 is defined** (`s12/instrument.py::
selfcheck`) as the targets where no pool member within 1.5 A of the pool best survives into the
production top-75 — so a contrast asking *"does the score's prefix retain the good members?"* is
near-circular there, and directly circular for the 128->75 arm. Filter-independent control
(length <= 13 vs >): +0.1826 vs +0.1940, **no gradient** — it is not a broad property, it is the 18.

**What is true:** *the distogram score is not a general anti-ordering. On 108 of 126 targets its
top-128 retains a marginally better best-member than a random 128 (NOT A RESULT). Its failure is
catastrophic and total on 14% of targets — +1.49 A, 0 wins of 18 — and those are exactly FAIL18.*

Audit 4 of `s32/s32_verify.py` now flags any document line quoting +0.1872 without its strata;
it currently flags 6 live lines across LEDGER, STATE, MULTIPLICITY (this file, above) and
CAUSAL_MAP. Self-tests ST6a/ST6b guard it.

---

### D4 — lane D's D1-T survives the attack, but a one-line geometric scalar beats all four scorers

`s32/results/s32_V_D_signadversary.json`, `s32/s32_V_D_signadversary.py`. Independent
replication from the raw universe, lane D's protocol (16 splits/target, shipped top-75 band),
plus three controls lane D's artefact does not contain.

```
scorer            transfer  nullPERM  nullXTGT  globalSGN  mean rho  mean|rhoB|  folds  per-target?
DIS                +0.1952   -0.0015   +0.0086   +0.0642    +0.0652    0.2673     5/5   PER-TARGET 3.05x
TYPICALITY         +0.3230   -0.0034   +0.1016   +0.2205    +0.2237    0.3731     5/5   PER-TARGET 2.93x
RG                 +0.3323   -0.0031   +0.0016   +0.0487    +0.0501    0.3779     5/5   PER-TARGET 4.44x
NOISE (falsifier)  -0.0012   -0.0017   +0.0006   +0.0073    +0.0100    0.1323     3/5   NOT A RESULT 0.07x
DIS_DEDUP          +0.1856   +0.0052   +0.0095   +0.0708    +0.0677    0.2657     5/5   PER-TARGET 2.85x
RG_DEDUP           +0.3210   +0.0039   +0.0057   +0.0451    +0.0488    0.3787     5/5   PER-TARGET 4.12x
```

- **Replication:** DIS +0.1952 (full band) / +0.1856 (deduplicated) against lane D's **+0.1890**.
- **`nullPERM` cannot test the claim being made.** Permuting `rr` inside B destroys *all*
  structure, so it is ~0 for every scorer **including pure noise**. It asks "is there any
  relation", not "is the relation's SIGN a property of the target".
- **`nullXTGT` is the null the claim needs** — another target's sign applied to this target's
  held-out half, 200 shuffles; its expectation is exactly what a single GLOBAL sign delivers.
  DIS +0.0086 and RG +0.0016, so those transfers really are per-target. **TYPICALITY +0.1016
  against a leave-fold-out global-sign baseline of +0.2205: two thirds of its large transfer
  needs no oracle at all.**
- **Band duplicates:** mean **7.7%** of the 75 members are exact coordinate duplicates, on
  **122/126** targets (mean 69.2 distinct). Deduplicating costs at most **0.012** of transfer —
  **the split-half independence attack fails; the finding survives it.**

**What does not survive is the framing.** *"A one-bit latent that four unrelated Hamiltonians
all respond to"* overstates the independence: **Rg**, one line of numpy, native-free and not a
Hamiltonian, has the largest latent of anything tested. The four all load on compactness and
bare compactness loads harder. And `in-band-ordering-is-per-target` already records native-free
compactness proxies at **0.24–0.37** — 0.3210 and 0.3204 land inside it. D1-T is a confirmation
on a new instrument, not a discovery.

**The ceiling worth quoting instead:** the same memory entry prices 2.0 Å at in-band rho ≈ 0.638.
A *perfect, free* per-target sign takes the best scorer to 0.2263 and Rg to 0.3323. **Even a
perfect sign oracle leaves the best in-band scorer 2–3x short of the useful range.**

### D5 — the sprint's most-quoted result has no producing script and no provenance

`s32/results/s32_D1_signtransfer.json` contains only `{note, n_splits, res}`: **no `provenance`
block** — no module, no `git_commit`, no `source_sha256`, no seed — where every other artefact in
this sprint carries `ST.provenance(__file__)`. And **no script in the repository produces it**:
`grep -rn "signtransfer|D1-T" --include=*.py` returns only `s32/s32_D4_signpred.py` and
`s32/s32_D5_signchain.py`, which *cite* its numbers as a ceiling.

So D1-T is four numbers, quoted as a ceiling by two downstream scripts and as a ledger headline,
that **cannot be re-run by anyone**. Charter §61 requires verifying experiments, artefacts and
seeds; this fails all three. Project memory `findings-prose-is-not-evidence-of-code` records this
exact shape twice already — **this is the third instance.** V-A8 covers the DIS row and all the
controls; AMBER, LEG_total and LEG_torsion must be re-emitted from a committed script with a
provenance block and a pinned seed before the report quotes them.

### D6 — `split_half_transfer` centres on the arm-family mean, not on production (caught in my own code)

`s32/s32_V_R_adversary.py`. My first accounting of lane R's 16-criterion search called
`ST.split_half_transfer(E)` on the per-criterion effects and reported

> *oracle-over-criteria −0.1168, **SPLIT-HALF TRANSFER −0.0254, CI [−0.0356, −0.0153], 22% of
> the oracle*** — a tight CI excluding zero and nearly 3x the best single arm.

**It is an artefact of the baseline.** `split_half_transfer` centres on `M.mean(1)`, the mean
**over criteria**, and that mean includes `typicality` at **+0.3434**. The number measured *"the
chosen criterion beats the average criterion"*, and nobody deploys the average of 16 criteria.
**A control matched to a different arm's magnitude is not a control** — contract rule 6, the
project's most repeated error, and the same shape as rule 4's `CTRL_SHRINK_ORACLE_Y`.

Corrected: choose `argmin` of the per-criterion mean effect **vs production** on one half of the
targets, evaluate that criterion **vs production** on the other half, 400 repeats:
**+0.0004, CI [−0.0091, +0.0116]** — centred on zero. The invalid version is kept in the
artefact under `criterion_search_INVALID_criterion_mean_baseline`, labelled and never quoted,
the same treatment `stats_lib.best_of_k_within` gives the `s25/temper.py` null.

**Generalisable rule this earns:** `split_half_transfer` nulls itself for the question *"does
WHICH setting wins transfer?"* and **not** for *"does the winner beat production?"*. Any lane
applying it to a grid whose columns are not all plausible arms gets an inflated effect with a
tight CI. State the baseline in the same sentence as the transfer.

---

## Lane V — comparisons emitted (charter §45)

**32 comparisons**, all on the **built-chain basis**, all paired to the production chain
**recomputed in lane R's own job** (rule 3), all tie-averaged over the argmin set, all
EXPLORATORY (lane V registered no hypothesis — these are the adversarial re-scoring of lane R's
arms, emitted so the search has an owner):

| family | comparisons | best single arm | its x MDE | the search accounted |
|---|---|---|---|---|
| 16 criteria x **2 directions** x {all branches, GEN4-only} | **64** | `d_to_C` min −0.0082 | 0.51x | **+0.0028, CI [−0.0054, +0.0128]** |

**Why 64 and not 32 — D7, an error of mine that the coordinator caught.** My first version fixed a
direction per criterion (`MAXIMISE = {"typicality"}`), assuming the name meant "how typical". Lane
R's `typicality` is a **distance**, so lower is more typical; I selected the *least* typical branch
and reported "+0.3434, a RESULT in the wrong direction", which the report glossed as *"the most
typical branch is much worse"*. **It is the most ATYPICAL branch that is much worse** — the expected
direction, supporting consensus-as-outlier-avoidance.

**The fix is not a corrected sign, it is no sign at all.** For a criterion with no a priori
direction, testing one direction is an **unregistered choice that halves the apparent
multiplicity**. Every criterion now runs both ways, both are reported, and the comparison count
doubles — which the out-of-sample accounting then charges for. Doubling the search did not move the
verdict: **+0.0004 → +0.0028, both CIs straddling zero.**

**Degeneracy note:** `vbond_mean` and `vbond_sd` return identical min and max with **203 of 203
branches tied** — every branch has ideal geometry, so those statistics are constant. **4 of the 64
comparisons are degenerate by construction**; `ramah` min (58 tied) and `posphi_frac` min (92 tied)
are heavily degenerate. Tie-averaging handles them, but the honest count is ~60.

**And the two-direction table changed the interpretation, which is why it was worth doing.** Every
criterion's ARGMAX is a **2–3x MDE, 5/5-fold RESULT in the WORSE direction** (+0.17 to +0.35) while
every ARGMIN is nothing. *The criteria reliably identify disasters and never identify winners* —
the same mechanism as D2/D3 and as consensus. On the **compute-matched GEN4-only** subset, where
production's own four starts are the whole branch set, **nothing is a result in either direction**
(max 0.46x): production's four starts contain no branch bad enough to be punished for.

Supporting arms, not endpoint comparisons: the ORACLE best branch (**−0.11**, split-half
transfer **5%**, NOT A SIGNAL), and the zero-information control (a uniformly random branch from
the same set, 300 draws, **+0.0049 ± 0.0079** vs production — *production's multi-start argmin is
worth 0.005 Å over a coin*).

**Positive control, which is what makes the sixteen nulls a measurement rather than a broken
pipeline:** `typicality` reports **+0.3434 at 2.35x MDE, 5/5 folds, 27W/75L** — a RESULT, in the
wrong direction. The family can emit one.

Lane V's in-band family (V-A8) emitted **8** further comparisons on the in-band Spearman basis,
DIAGNOSTIC, never differenced against a chain RMSD, one of which is a pure-noise falsifier.

---

## Standing attack surfaces for this sprint (charter §53)

| lane | expected claim | the attack it must survive |
|---|---|---|
| R | ORACLE-best-over-N-restarts on the projection | **order statistic.** Price as best-of-N with `best_of_k_within`, and quote `split_half_transfer`, not the oracle. Given D1, N restarts of a 1e-13-amplifying map is a lottery: the right null is N *perturbed* restarts of the SAME branch structure. |
| Q | a formulation that escapes "classical ordering in quantum notation" | **classical equivalent** (charter §16 item 13). If a cheap classical algorithm solves it exactly, the quantum framing is decoration. |
| P | in-band selection skill | in-band skill has been zero or wrong-signed every time (D2 is the newest instance). Demand the null, the fold structure, and a plausible zero-information control (a constant, not a uniform draw). |
| D | a physical observable escaping theorem G1 | name **which hypothesis of G1** it violates. A deterministic minimiser started at a structure is still a function of that structure. |
| coordinator | any cross-lane synthesis | contract rule 7: in S31 every single-lane result held and **every** cross-lane synthesis failed, 4/4. Re-derive from both lanes' raw artefacts in one script before it is written down. |

---

## Lane D — PHYSICAL RESPONSE AND DYNAMICS (prereg `s32/PREREG_S32_D.md`, commit `34973b1b`)

All lane-D comparisons, logged as emitted. Basis is named on every row. `rho` rows are **in-band
Spearman against true Cα-RMSD inside each target's shipped top-75 band, per target, then aggregated
over n = 126** — *not* a chain RMSD, and never differenced against one. ORACLE `rr` is used as an
evaluation label only; no arm reads it.

Note on `stats_lib.compare` for the `rho` rows: it is LOWER IS BETTER, so its BETTER/WORSE verdict
labels are **inverted** for a correlation and are not quoted. Only effect, se, mde, ×MDE, the
fold-clustered CI and `folds_same_sign` are read.

### Registered — D1-E / D1-L / D0-X (`s32/results/s32_D1_inband.json`, 2026-09-21)

| # | reg? | comparison | basis | effect | se | ×MDE | folds | verdict |
|---|---|---|---|---|---|---|---|---|
| D-1 | R | **AMBER** ff14SB/GBn2 single point, in-band ρ vs 0 | top-75 band, n=126 | **+0.0000** | 0.0202 | **0.00** | 3/5 | **NOT A RESULT — exactly zero skill** |
| D-2 | R | AMBER, in-band ρ, `Rg` partialled | top-75 band | +0.0077 | — | — | — | NOT A RESULT |
| D-3 | R | AMBER, whole-pool ρ (K=500) | 500 pool | −0.0266 | 0.0213 | 0.45 | 4/5 | NOT A RESULT, **wrong sign** |
| D-4 | R | `DIS` (shipped cost), in-band ρ | top-75 band | +0.0652 | 0.0282 | 0.83 | 4/5 | **NOT MEASURED** (reproduces S31 §9's 0.83×) |
| D-5 | R | `LEG_total` (Legacy, 11 terms, default weights), in-band ρ | top-75 band | +0.0376 | 0.0316 | 0.42 | 3/5 | NOT A RESULT |
| D-6 | R | `LEG_torsion`, in-band ρ | top-75 band | +0.0444 | 0.0243 | 0.65 | 5/5 | NOT A RESULT (reproduces S31 §9's 0.65×) |
| D-7 | R | `RG` (compactness control), in-band ρ | top-75 band | +0.0501 | 0.0403 | 0.44 | 3/5 | NOT A RESULT |
| D-8 | E | `DIS`, whole-pool ρ | 500 pool | +0.5678 | 0.0327 | 6.20 | 5/5 | the global axis is outlier detection |
| D-9 | E | `LEG_total`, whole-pool ρ | 500 pool | +0.3073 | 0.0339 | 3.23 | 5/5 | ditto |
| D-10 | R | **D0-X** chiral Cα pseudo-torsion scalar, in-band ρ | top-75 band | +0.0607 | 0.0375 | 0.58 | 4/5 | **NOT MEASURED** |
| D-11 | R | **D0-X** same scalar, whole-pool ρ | 500 pool | **+0.3302** | 0.0413 | **2.85** | 5/5 | real skill globally, lost in band |
| D-12 | R | **D0-X** in-band share of chiral variance | top-75 band | **23.6%** | — | — | — | **registered prediction (<15%) FAILED** |
| D-13 | R | **D0-X** in-band share of chiral variance | 500 pool | 90.9% | — | — | — | — |

Exploratory per-term Legacy in-band ρ (8 of 11 terms are **constant in band** on 18–98 of the 126
targets and are reported as such, not as nulls): `hbond_local` +0.0308 (0.38×), `contact` −0.0279
(0.37×), `solvation` +0.0239 (0.31×). Count these 3 toward multiplicity.

### Registered — D3-M, physics as MOVER (`s32/results/s32_D3_mover_A_0_1.jsonl`, running)

4 restraint rungs × {3 ORACLE cosines + 1 built-chain endpoint + 1 magnitude-matched control}.
**Registered before the numbers: primary = median `cos(Δd, −e_prod)` ≥ +0.10; endpoint = built-chain
CA-RMSD paired against the production chain rebuilt in the same job.** Control draws: 32 per target
per rung, **own distribution reported (mean, sd, max)**, never a best draw (contract rule 10).

## LANE P (pool, filter, selection) — prereg `s32/PREREG_S32_P.md` @ `33dfe0d3`

Basis is named on every row. `OB` = ORACLE band (top-24 by the native label `a`, a collider).
`NFB` = native-free band (the score's own top-24, not a collider).

| # | comparison | basis | registered? | emitted |
|---|---|---|---|---|
| P-01 | R1 universe→500 BLOSUM, BEST vs matched random-500 (200 draws) | cloud label `rr`, ORACLE | REGISTERED (H-P2a) | −0.0735, 0.89× — **NOT MEASURED** |
| P-02 | R1 universe→500 BLOSUM, MEAN vs matched random-500 | cloud label, ORACLE | REGISTERED (H-P2a) | −0.3617, 4.00×, 5/5, 108W/18L — **BETTER** |
| P-03 | R2 500→128 score, BEST vs matched random-128 | cloud label, ORACLE | REGISTERED (H-P2a) | +0.1879, 1.07×, 4/5, 72W/54L — **WORSE** |
| P-04 | R2 500→128 score, MEAN vs matched random-128 | cloud label, ORACLE | REGISTERED (H-P2a) | −0.8687, 4.30×, 5/5 — **BETTER** |
| P-05 | R3 128→75 score, BEST vs matched random-75-of-128 | cloud label, ORACLE | REGISTERED (H-P2a) | +0.0684, 0.98× — **NOT MEASURED** |
| P-06 | R3 128→75 score, MEAN vs matched random-75-of-128 | cloud label, ORACLE | REGISTERED (H-P2a) | −0.0349, 0.71× — **NOT MEASURED** |
| P-07 | R2b 500→75 score (THE DEPLOYABLE FILTER), BEST vs random-75-of-500 | cloud label, ORACLE | REGISTERED (H-P2a) | +0.2286, 1.18×, 4/5, 63W/63L — **WORSE** |
| P-08 | R2b 500→75 score, MEAN vs random-75-of-500 | cloud label, ORACLE | REGISTERED (H-P2a) | −0.9030, 4.01×, 5/5, 110W/16L — **BETTER** |
| P-09 | BLOSUM top-75 vs score top-75, BEST | cloud label, ORACLE | EXPLORATORY | −0.2021, 1.06× TYPE-M, 62W/53L |
| P-10 | BLOSUM top-75 vs score top-75, MEAN | cloud label, ORACLE | EXPLORATORY | +0.7360, 3.55×, 5/5 — **WORSE** |
| P-11 | P3-1 `Var(U)/Var(V)` over the 128 > 1 | cloud, ORACLE | REGISTERED | median 2.62, 78% of targets — **PASS** |
| P-12 | P3-2 in-band ρ(CONS, a) negative on ≥80% | OB, ORACLE | REGISTERED | −0.2745 ± 0.0366, 78.6% — **FAIL on the 80% bar, sign confirmed** |
| P-13 | P3-2 shuffled-U collider control | OB, ORACLE | REGISTERED (control) | +0.1387 ± 0.0214 — collider alone gives the WRONG sign |
| P-14 | P3-3 μ recovered by least squares from {a, d} | cloud, ORACLE | REGISTERED | rel. resid 2e-14, 126/126; shuffled-d control 0.373 — **PASS** |
| P-15 | P3-4 zero-μ control, ρ(CONS, a′) in band | OB, ORACLE | REGISTERED (control) | −0.2745 → **+0.7046 ± 0.0096** — **PASS** |
| P-16 | P3-5 price curve, cos\* for positive in-band ρ | OB, ORACLE | REGISTERED | cos\* = **0.2137** (6 draws/target × 126) |
| P-17 | P3-6 cos(μ̂, μ) for 7 native-free directions | cloud, ORACLE-scored | REGISTERED | max +0.0568 ± 0.0434 (ARGMIN); all **below cos\*** |
| P-18 | in-band ρ of CONS on the NATIVE-FREE band | NFB, ORACLE label only | EXPLORATORY | **+0.2531 ± 0.0434**, 73% positive — positive in-band skill EXISTS |
| P-19 | `Var(U)/Var(V)` on the native-free band | NFB, ORACLE | EXPLORATORY | median 4.89 (mean 19.0) — worse than the ORACLE band |
| P-20 | P5 strata `Var(U)/Var(V)`: tail vs other108 | cloud, ORACLE | REGISTERED (H-P5) | 13.60 / 7.43 (filter-independent) vs 4.71 — **PASS** |
| R-1 | R | PROD(cache cloud) vs PROD(s29 cloud) | R | built chain n=126 | per-target abs mean 0.0133, p90 0.0245, max 0.4168 | n/a (a floor, not an effect) | >1e-3 on 71/126, >0.1 on 3/126; independent replication of lane V's ULP result |
| R-2 | R | lam=0.3 arm vs lam=0 arm (same cloud, same job) | R | built chain n=126 | +0.0055 (SE 0.0057) | 0.35x | NOT MEASURED -- the +0.0107 from the cache is a cross-job difference |
| R-3 | R | PROD chain vs orthogonal null sqrt(e^2+d^2) | R | built chain n=126 | +0.0118 | 0.24x | NOT MEASURED -- P1.2 HOLDS for production |
| R-4 | R | ORACLE best branch (7 family subsets) vs PROD | R | built chain n=126 | GEN4 -0.0355 / PROD8 -0.0938 / MEM75 -0.0720 / RAND0 -0.0783 / ALL -0.1165 | 1.73 / 2.92 / 2.70 / 2.92 / 3.43 | ALL ORACLE. PROD8 replicates S31-D's -0.0938 at 2.92x exactly |
| R-5 | R | order-statistic curve (8 sizes x 24 draws) + best_of_k_within + split-half on 6 families | R | built chain n=126 | ALL: oracle -0.1210, split-half -0.0046 | 4% of oracle | 96% of the ceiling is an order statistic; RAND0 control transfers 3% |
| R-6 | R | in-band rho, 11 criteria | R | within-target rank, n=126 | d_to_C +0.1122 ... rama_nlp +0.0153 | 5 of 11 have fold CI excluding zero | P3.1, P3.2, P3.3 ALL FALSIFIED |
| R-7 | R | SEL_<criterion>_<subset>, 11 criteria x 2 subsets = 22 arms | R | built chain n=126 | best SEL_d_to_C_GEN4+MEM75 -0.0100 | 0.61x | P3.4 FALSIFIED; nothing clears MDE, before charging the search over 22 |
| R-8 | R | BRANCHMEAN, RANDBRANCH (5 draws), MEDOID_ONLY, OBJARGMIN x4 | R | built chain n=126 | RANDBRANCH 3.2106 (draw sd 0.0028) vs PROD 3.2105 | n/a | production's argmin is worth 0.0001 A over a coin |
| R-9 | R | MEDOID_EXTRA, SCALE_NF, SCALE_NF_MED, SCALE_GRID x8 | R | built chain n=126 | MEDOID_EXTRA -0.0064 / SCALE_NF +0.7219 | 0.35x / 2.24x | NOT MEASURED / WORSE 5/5 folds |

### Registered — D1-N, the matched permutation null (`s32/results/s32_D1_signrandom.json`)

The null is **within-band label permutation**: the score vector and the band size are held fixed and
the ORACLE `rr` labels are permuted inside each target's own top-75 band, 24 draws per target, own
distribution reported. This is the control matched to *this* operator's space (contract rule 6) and
it answers a question the mean ρ cannot: **is the per-target |ρ| bigger than chance?**

| # | reg? | scorer | mean\|ρ\| in band | matched null | ratio | ×MDE | folds |
|---|---|---|---|---|---|---|---|
| D-14 | R | **AMBER** | **0.1779** | 0.0948 | **1.88** | **2.36** | **5/5** |
| D-15 | R | `DIS` | 0.2496 | 0.0939 | 2.66 | 3.05 | 5/5 |
| D-16 | R | `LEG_total` | 0.2819 | 0.0908 | 3.10 | 3.52 | 5/5 |
| D-17 | R | `LEG_torsion` | 0.2142 | 0.0947 | 2.26 | 2.77 | 5/5 |

**Reading, and it changes what D-1 means.** AMBER's mean in-band ρ is exactly `+0.0000` (62W/64L)
*while* its per-target |ρ| is **1.88× the permutation null at 2.36× MDE, 5/5 folds**. The two are
not in conflict: `Var(ρ_true) > 0` with `E[ρ_true] = 0`. **In-band ordering content exists on every
scorer measured, including the one with exactly zero mean skill; what is missing is the per-target
SIGN.** This reproduces memory `in-band-ordering-is-per-target` on a channel it had never been
measured on. It is **not** a claim that the sign is recoverable — memory
`in-band-signal-limited-not-sample-limited` records a flat learning curve for exactly that.

### Registered caveat on D-1, measured not asserted

The cached AMBER energies are **single points on the ideal-geometry rebuild and are
clash-dominated**: median energy over the 126 targets' pool medians is **29,668 kcal/mol**, the p90
is **6.83e6**, the median per-target maximum is **1.80e16**, and **58.6% of every 500-pool sits above
1e4 kcal/mol**. D-1 therefore prices **unminimised** AMBER, not AMBER. Rung **D2-R** scores the
relaxed energy on the same band and is the arm that settles it.

### Registered follow-up — D1-S: is the missing per-target in-band SIGN shared across Hamiltonians?

Motivated by D-14…D-17: ordering content exists on every scorer, the sign does not. If two
*independent* Hamiltonians agreed on the sign, a native-free vote would recover it. **They do not,
and where they depart from chance they depart the wrong way.** `s32/results/s32_D1_signshare.json`.
Ten pairwise tests; at 10 comparisons the Bonferroni bar is |z| > 2.81.

| # | reg? | pair | sign agreement (chance 0.50) | z |
|---|---|---|---|---|
| D-18 | E | AMBER vs `LEG_total` | **0.365** | **−3.0** (anti-agreement survives Bonferroni) |
| D-19 | E | AMBER vs `RG` | **0.317** | **−4.1** |
| D-20 | E | `LEG_total` vs `RG` | **0.683** | **+4.1** |
| D-21 | E | AMBER vs `LEG_torsion` | 0.587 | +2.0 (does not survive) |
| D-22 | E | AMBER vs `DIS` | 0.532 | +0.7 |
| D-23..27 | E | the other five pairs | 0.421–0.571 | \|z\| ≤ 1.8 |

**Mechanism, and it is the compactness confound again.** `LEG_total`'s in-band sign agrees with
plain `Rg`'s on **68.3%** of targets while AMBER's agrees on **31.7%** — Legacy tracks compactness in
band and AMBER anti-tracks it, which is *why* the two disagree. A native-free vote therefore fails:

| # | reg? | arm | in-band ρ | se | ×MDE | folds |
|---|---|---|---|---|---|---|
| D-28 | E | AMBER oriented by `LEG_total`'s sign (native-free) | **−0.0508** | 0.0197 | 0.92 | 4/5 |
| D-29 | E | AMBER oriented by `DIS`'s sign (native-free) | −0.0108 | 0.0202 | 0.19 | 4/5 |

Both are NOT MEASURED or NOT A RESULT, and D-28 points the **wrong way** — orienting AMBER by
Legacy is worse than leaving it unoriented.

**The ORACLE price of the sign, labelled ORACLE / NOT DEPLOYABLE** (it is `mean|ρ|`, i.e. what the
scorer would be worth if an oracle supplied the per-target sign — a price, never a method):
`RG` +0.3709, `LEG_total` +0.2819, `DIS` +0.2496, `LEG_torsion` +0.2142, **AMBER +0.1779**, all at
4.9–5.7× MDE. Memory `in-band-ordering-is-per-target` prices 2.0 Å at ρ = 0.638, so **even a perfect
sign oracle on the best of these leaves the in-band channel short of the crossing price by more than
a factor of two.**

**Limitation, stated rather than discovered later:** sign agreement is dominated by the targets where
|ρ| sits in the noise, and the two ρ's being compared are estimated against a **common** referent
(`rr`), so memory `shared-referent-floor` applies to their correlation. The *anti*-agreements are the
robust half of this table; the near-chance rows are uninformative either way.
| P-21 | sign(⟨pc1, μ⟩): 9 native-free sign rules, LFO polarity | cloud, ORACLE bit | EXPLORATORY (9-rule family; read with that multiplicity) | best accuracy 0.5952 ± 0.0226 (1.51×); best payoff **+0.0032 Å** |
| P-22 | the ORACLE bit: sign + ONE global LFO step | cloud, ORACLE | EXPLORATORY | **+0.1835 Å** of the 3.0776 Å common mode; random sign −0.0016 ± 0.0030 (200 draws) |
| P-23 | the ORACLE bit + ORACLE per-target step (1 bit + 1 real) | cloud, ORACLE | EXPLORATORY | +0.4589 Å — **not one bit; quote P-22 for one bit** |
| P-24 | PLACEBO: cos²(pc1, c−t′) for mismatched same-length natives | cloud, ORACLE | EXPLORATORY (control) | 0.2024 vs true 0.2277 vs random 0.0316 — **87.1% of the alignment is placebo** |

### Registered — D0-X-C: the circularity control the coordinator asked for. **THE SYNTHESIS IS CIRCULAR.**

The coordinator proposed that lane D's *"the top-75 filter removes 79% of the pool's chiral
variance"* and lane V's *"the score prefix halves the pool's quality spread"* are **the same event
from two sides**, and that the arm implied is to apply the chiral observable at 500 → 128 where the
chiral variance still exists. He asked to be told if the two halves are one variance measured twice.
**They are.** `s32/results/s32_D0X_circularity.json`.

**Control: the variance each coordinate RETAINS under the shipped top-75 filter, with a size-matched
random-75 control beside it, and with the EXACT achiral twin of the chiral coordinate** — `mean cos τ`
against `mean sin τ` over the same Cα pseudo-torsions, same functional family, same scale, differing
in nothing but parity under reflection (S30's own twin design).

| coordinate | variance retained, score-75 | random-75 control |
|---|---|---|
| `DIS` — *the filter's own criterion* | **0.037** | 0.964 |
| `RG` | 0.111 | 0.979 |
| **`CHI_even_cos` — the ACHIRAL TWIN** | **0.196** | 0.978 |
| `E2E` | 0.244 | 0.957 |
| **`CHI_odd_sin` — the CHIRAL coordinate** | **0.271** | 0.975 |
| **`RR_ORACLE` — lane V's quantity** | **0.272** | 0.965 |
| `LEG_total` | 0.323 | 0.973 |

| # | reg? | comparison | effect | se | ×MDE | folds | W/L |
|---|---|---|---|---|---|---|---|
| D-30 | R | chiral retention **minus its own achiral twin's** | **+0.0755** | 0.0114 | **2.35** | 5/5 | 53/73 |

**Three readings, and all three go against the proposal.**

1. **The filter shrinks everything, not chirality.** The random-75 control retains 0.96–0.98 on every
   coordinate, so the shrinkage is selection and not subsetting — but it is *general* selection. The
   chiral coordinate is the **second-least** shrunk thing measured.
2. **Chirality is shrunk LESS than its own achiral twin**, by +0.0755 at 2.35× MDE, 5/5 folds. *"The
   filter removes the chiral variance"* is exactly backwards as a chirality-specific claim.
3. **Lane D's 0.271 and lane V's 0.272 are the same number to three decimals.** They are not two
   facts. They are one fact — *a score filter narrows every coordinate correlated with the score* —
   and chirality and ORACLE quality happen to correlate with `DIS` at similar strength. **The
   synthesis is circular and must not be written down.**

### D0-X-P: the proposed 500 → 128 arm's precondition, and why it is NOT RUN

| # | reg? | global ρ vs `rr`, K=500 pool | effect | se | ×MDE | folds |
|---|---|---|---|---|---|---|
| D-31 | R | `CHI_odd_sin` (chiral) | +0.3302 | 0.0413 | 2.85 | 5/5 |
| D-32 | R | **`CHI_even_cos` (the ACHIRAL TWIN)** | **−0.3359** | 0.0440 | 2.72 | 5/5 |
| D-33 | R | chiral, `DIS` partialled | +0.1302 | 0.0324 | 1.44 | 5/5 |
| D-34 | R | chiral, `Rg` partialled | +0.2574 | 0.0498 | 1.84 | 5/5 |
| D-35 | E | `DIS`, chiral partialled | +0.4657 | 0.0289 | 5.75 | 5/5 |

rank-corr(chiral, `Rg`) = **+0.464**; rank-corr(chiral, `DIS`) = **+0.374**.

**The arm is not run, and the reason is a derivation, not a budget.** The chiral coordinate's achiral
twin reaches **the same magnitude with the opposite sign** (−0.3359 against +0.3302). By S30's own
twin criterion — *"whatever WRITHE can do, its own reflection-invariant shadow already does"* — **the
G1 escape is not what is doing the work here, for the third time in this project** (WRITHE S30-L26,
XTWIST S31 §9, this). What is left is a weakly-independent global coordinate (+0.1302 past `DIS`,
44% of it shared with `Rg`) proposed as an addition to a global filter — which is **S30-L21's closed
question**: the available fields span ~2 directions, and the best ORACLE global combination of 21 of
them reaches 0.1693, leave-fold-out **0.0948**. Spending AMBER-free CPU on it is allowed; spending
the lane's remaining AMBER window on it is not, and the mechanism that motivated it is refuted above.

### The falsifier that fired, kept in the record rather than overshadowed

D-12: the registered prediction *"in-band share of chiral variance < 15%"* **FAILED at 23.6%.** The
79%/21% split it was meant to explain is, per D-30, **not chirality-specific**. Both statements stand.

### Registered — D1-T: the per-target in-band sign is a real latent and it TRANSFERS

The adversarial check on D-14…D-17. `Var(ρ) > 0` could in principle be dispersion without a stable
per-target *cause*. **Split-half within each target**: the sign is estimated on a random half A of
the top-75 band and applied to the held-out half B, 16 random splits per target, own distribution
averaged (contract rule 10). **Matched null: the same sign, with the labels on B permuted.**
`s32/results/s32_D1_signtransfer.json`.

| # | reg? | scorer | transfer `sign(ρ_A)·ρ_B` | matched null | excess | se | ×MDE | folds |
|---|---|---|---|---|---|---|---|---|
| D-36 | R | **AMBER** | **+0.1125** | −0.0056 | **+0.1181** | 0.0164 | **2.57** | **5/5** |
| D-37 | R | `DIS` | +0.1890 | −0.0033 | +0.1923 | 0.0226 | 3.04 | 5/5 |
| D-38 | R | `LEG_total` | **+0.2263** | −0.0003 | +0.2266 | 0.0237 | **3.41** | **5/5** |
| D-39 | R | `LEG_torsion` | +0.1499 | −0.0008 | +0.1507 | 0.0196 | 2.75 | 5/5 |

**A sign learned on one half of a target's band predicts the ordering on the other half.** The latent
is a property of the **target**, not of the sample — which is what `Var(ρ) > 0` alone could not
establish. AMBER goes from **+0.0000 to +0.1125** and Legacy from **+0.0376 to +0.2263** once the
sign is supplied.

**ORACLE / NOT DEPLOYABLE.** Estimating the sign on half A requires native labels on half A. This
prices *what a partial oracle buys*; it is not a method. And memory
`in-band-signal-limited-not-sample-limited` records that a set-transformer over the full signed
deviation map has a **flat learning curve** when asked to supply this sign natively — so the latent
existing and the latent being *nativly* recoverable are different claims, and only the first is
established here.

### Registered — D2-R: the 9–16-residue reality check. **BOTH HORNS MISSED. My own hypothesis is falsified.**

`s32/results/s32_D2_basin_0_1.jsonl`. 10 targets (**native-free, deterministic**: the first 2 of each
of the 5 frozen folds in pinned pdb order), 24 candidates each from the head of the shipped top-75
band, free AMBER relaxation (`k=0`, `steps=200`, `tol=5.0`), **97.5% converged**, 4.2 s/structure.

| # | reg? | registered horn | threshold | measured | fired? |
|---|---|---|---|---|---|
| D-40 | R | **CONTRACTION** — the pool collapses to a common attractor | median spread ratio ≤ **0.70** | **0.972** (1/10 targets below) | **NO** |
| D-41 | R | **FROZEN** — the trajectory carries no new bit | median per-candidate move ≤ **0.30 Å** | **0.5006 Å** | **NO** |

**H-D2 is FALSIFIED and I registered it as my own expectation.** Relaxation moves each candidate
~0.50 Å (max 1.11 Å) and leaves the band's spread essentially untouched (1.4842 Å → 1.3972 Å). *At
9–16 residues the candidates are not in one basin, and they are not frozen.* **The scope excuse this
lane was entitled to use is not available**, and any future negative on dynamics must name a
different reason.

### Registered — D2-E: minimisation does **not** rescue AMBER's in-band ranking. D-1's caveat is answered.

n = 10 targets, 24 candidates each — **DIAGNOSTIC, no MDE verdict is claimed on this subset.**

| # | reg? | AMBER energy | in-band ρ | se |
|---|---|---|---|---|
| D-42 | R | **single point** (as in D-1) | +0.0006 | 0.1250 |
| D-43 | R | **converged relaxed** | **−0.0089** | 0.1062 |
| D-44 | R | single point, `Rg` partialled | −0.0036 | 0.1179 |
| D-45 | R | relaxed, `Rg` partialled | +0.0254 | 0.1030 |

The relaxed energies are **physically well-behaved** — median −887 to −293 kcal/mol, per-target sd
**1–9 kcal/mol on 8 of 10 targets** — so D-1's zero is **not a clash artefact**. And the per-target
pattern is the headline again, on a completely independent computation: relaxed-energy ρ per target
is {+0.056, +0.092, −0.090, +0.543, −0.686, +0.263, −0.125, +0.117, −0.353, +0.093}, **mean −0.0089,
mean |ρ| 0.242.** *Zero mean, real dispersion, random sign* — reproduced on converged all-atom
physics rather than cached single points.

### The mechanism, and it is the sharpest thing in rung D2

```
energy, median over candidates      3.18e9  ->  -558 kcal/mol
Cα motion required to get there                 0.50 Å
mean CA-RMSD of the band to native  3.1931  ->  3.1849   (-0.0083, nothing)
BEST member of the band             2.5358  ->  2.5741   (+0.0382, WORSE)   ORACLE / NOT DEPLOYABLE
rank preservation rho(rr_in, rr_out)            0.9208
closest pair in the band            0.0643  ->  0.0990 Å
```

**Nine orders of magnitude of energy are relieved by half an Ångström of Cα motion.** The clash lives
in the **sidechains and hydrogens**, which by §0.2c are a *deterministic* function of `(seq, φ, ψ)` —
modal rotamers and frozen local frames, identical construction for every candidate. **AMBER's entire
dynamic range on this pool is spent on a rotamer-placement artefact that carries no candidate
information, and the 0.5 Å of Cα motion it buys leaves the band's quality unchanged (−0.008 Å) and
its best member slightly worse (+0.038 Å).** That is *why* the force field cannot rank here, stated
as a mechanism rather than as a null.
| R-10 | R | ladder rungs x5: observed price vs ORTHOGONAL isotropic null | R | built chain n=126 | prod +0.0118 / bestm +0.0067 / best1 -0.0076 / sp10 -0.2121 / sp20 -0.2186 | 0.24 / 0.13 / 1.16 / 4.65 / 4.54 | P1.2 HOLDS on prod+bestm, FALSIFIED on both sparse rungs (5/5 folds) |
| R-11 | R | ladder rungs x5: cos_align + on-manifold null sqrt(e^2-d^2) | **E** | built chain n=126 | cos: prod -0.061, bestm -0.004, best1 +0.055, sp10 +0.380, sp20 +0.391 | SE 0.016-0.022 | EXPLORATORY, added after 10 of 126 rows |
| P-25 | SPREAD (score floor + max w'Bw over survivors) vs PROD, 4 floors | **CLOUD SCREEN**, gated | EXPLORATORY (coordinator arm) | +0.6447 / +0.1796 / +0.0436 / +0.0076 at floors 25/50/60/75% — 1.96× WORSE to 0.11× NOT MEASURED |
| P-26 | SPREAD vs RANDFLOOR (random 75 from the SAME survivor set) — the control that decides whether spread does work | CLOUD SCREEN | EXPLORATORY (control) | +0.4196 / +0.0864 / −0.0104 / −0.0039 — **spread does no work at any floor** |
| P-27 | SPREAD with the floor chosen leave-fold-out vs PROD | CLOUD SCREEN | EXPLORATORY | +0.0076, 0.11×, 59W/67L — **NOT MEASURED**; LFO picks the tightest floor (0.75) on 5/5 folds, i.e. converges on production |

### LANE Q, continued — Q1/Q2/Q4

| # | lane | family | arms | registered? | basis | emitted | notes |
|---|---|---|---|---|---|---|---|
| Q-4 | Q | Q1-T1 sufficiency | identity / `a`-affine-in-`t` / substitute `P_aff{W}t` for `t` / KKT | REGISTERED | derivation + numerical falsifier, n=126 | 1.3e-11 / 5.8e-13 / Δw 7.0e-12, Δx 2.9e-12 Å / KKT 2.0e-12 | **Not a comparison.** A theorem with a falsifier that can fail; falsifier did not fire. Resolves S31 §20.3. |
| Q-5 | Q | Q1-T2 sensitivity | in-hull gain / orthogonal gain / `dx = P_aff(S)dt` | REGISTERED | derivation + finite differences, n=126×4 | 0.9999999999 (dev 1.2e-08) / 4.0e-09 / rel 3.2e-07 | Step size set by an h-sweep (residual scales as 1/h = roundoff), not by the outcome. Active set moved 0/504. |
| Q-6 | Q | Q1-T3 cardinality | support of the UNCONSTRAINED convex optimum, K=500 | REGISTERED | **CA cloud**, **ORACLE / NOT DEPLOYABLE** | `s*` mean 10.06, median 10, p90 13, max 23; 61.1% ≤ 10 | Closure is by MONOTONICITY and does not depend on where `s*` falls. Value 1.1535 ± 0.0669 is **not** differenced against S32-L1's 1.1139 built chain (different basis and frame). |
| Q-7 | Q | Q4 pricing curve | emitted RMSD vs `r` retained principal components of `P_aff{W}t`, r = 0…33 | REGISTERED | **CA cloud**, **ORACLE / NOT DEPLOYABLE** | 3.0532 (r=0) → 2.0312 (r=6) → 1.8290 (r≈33) | A **ceiling curve in real numbers**, never differenced against a bit count. Each rung is a per-target ORACLE quantity, not a best-of-K. |
| Q-8 | Q | hull projection vs shrinkage | `PROJ − SHRINK` (norm-matched to PROJ's own displacement), 7 noise levels × 8 draws | REGISTERED (control) | **CA cloud**, **ORACLE / NOT DEPLOYABLE** | ε 0.5/1.0/1.5: **+0.370 (4.71×) / +0.265 (3.17×) / +0.170 (1.96×) WORSE**; ε 2.0/2.5/3.0: 0.89×/0.01×/0.77× NOT MEASURED; ε 4.0: −0.137 (1.59×) better, **median −0.074 vs mean −0.137** | **The control fired and killed the arm.** 7 comparisons emitted in one family; the one that clears MDE has a mean/median ratio of 1.85 (concentration warning) and sits under an isotropic error model that is the most favourable geometry, not a neutral one. |

### Registered — D4-S: can the per-target sign be supplied NATIVE-FREE?  **No.** `s32/results/s32_D4_signpred.json`

The deployable form of D1-T. Leave-one-**fold**-out ridge (inner leave-fold-out alpha selection, so
no test fold touches the hyperparameter) over **20 native-free per-target features** — band spread
and its sd, `Rg` mean/sd, the predicted-vs-realised `Rg` disagreement, the chiral summary and its
in-band sd, distogram posterior sd, in-band pair-distance sd, `expected` minus band mean, sequence
composition (hydrophobic/Gly/Pro), peptide-DB fraction, BLOSUM similarity, per-scorer relative score
spread, and `n`. **No feature reads `rr` or `nat_ca`.** Readout is directly comparable to D1-T:
`oriented ρ = mean over targets of sign(ρ̂_lfo) · ρ_true`.

| # | reg? | scorer | oriented (native-free) | **CONSTANT +1** (plausible zero-info control) | label-shuffled (32 draws) | excess over shuffled | ×MDE | folds |
|---|---|---|---|---|---|---|---|---|
| D-46 | R | **AMBER** | **−0.0355** | +0.0000 | −0.0645 | +0.0290 | 0.54 | 3/5 |
| D-47 | R | `DIS` | +0.0674 | +0.0652 | +0.0624 | +0.0050 | 0.10 | 4/5 |
| D-48 | R | `LEG_total` | +0.0373 | +0.0376 | +0.0252 | +0.0121 | 0.22 | 3/5 |
| D-49 | R | `LEG_torsion` | +0.0444 | +0.0444 | +0.0392 | +0.0052 | 0.59 | 4/5 |

**Nothing beats the plausible zero-information control, and AMBER is beaten by it.** Every excess
over the matched label-shuffled control is **0.10–0.59× MDE — NOT A RESULT.** Leave-fold-out R² on
ρ itself: AMBER **−0.0563**, `DIS` +0.0256, `LEG_total` −0.0228, `LEG_torsion` −0.0158. Sign accuracy
0.381 / 0.548 / 0.595 / 0.556 against a chance of 0.50 *and against a constant-`+1` baseline that is
already right on 50–58% of targets.*

**The ORACLE gap this leaves, stated as the lane's closing number:** D1-T's partial-oracle ceiling is
`LEG_total` **+0.2263**, AMBER **+0.1125**; the best native-free arm here is `DIS` **+0.0674**, which
is its own constant control. **The sign is real, it transfers inside a target, and no per-target
native-free summary feature carries it.** This is memory
`in-band-signal-limited-not-sample-limited`'s flat learning curve, reproduced on an explicit
20-feature tabular channel rather than a set transformer.

**A mechanism I guessed and the data did NOT support, recorded because I looked.** I expected the
model to get the sign right only where |ρ| is small. Splitting each scorer at its own median |ρ|:
AMBER 0.365/0.397 (low/high), `DIS` **0.460/0.635**, `LEG_total` **0.619/0.571**, `LEG_torsion`
0.508/0.603. **The split is scorer-dependent and goes the opposite way on `DIS` and `LEG_total`, so
the guess is not supported and is not claimed.** (`DIS`'s 0.635 on its high-|ρ| half is a **post-hoc
subgroup selected by an ORACLE quantity** — |ρ| is not knowable natively — and is recorded as an
observation, not a lead.)
| P-28 | ORACLE top-75 through the **DEPLOYED** uniform average vs production | CLOUD, ORACLE | REGISTERED (H-P1) | **−1.0853, 4.26× MDE, 125W/1L, 5/5** |
| P-29 | ORACLE top-m (m=5, ORACLE-global) through the deployed average vs production | CLOUD, ORACLE | REGISTERED (H-P1) | **−1.5670, 5.17× MDE, 126W/0L, 5/5** |
| P-30 | operator law `out ~ set_mean + set_best`, both prefix families, m=3..75 | CLOUD | REGISTERED (H-P2b-1) | 0.6730 / 0.3256, **ratio 2.07 — the registered ≥5× PREDICTION FAILS** |
| R-12 | R | RANDSPARSE s=10 (3 draws) + SCORESPARSE: price vs own orthogonal null | **E** | built chain n=126 | RANDSPARSE price +0.1988 (draw sd 0.0075), cos +0.019 (draw sd 0.0094); SCORESPARSE price +0.1371, cos -0.076 | 0.22 / 1.46 | THE DECISIVE CONTROL: same s=10 pays +0.0002 (ORACLE members) vs +0.1988 (native-free members). Sparsity explains nothing |
| R-13 | R | spearman(d, top-75 member spread), both sides NATIVE-FREE | **E** | n=126 | rho +0.9646 (partial|n +0.967, partial|e +0.953, 5/5 folds) | null: within-n permutation mean +0.099, p99.9 +0.377, max +0.466 over 4000 draws | ratio cv 0.467 -> MONOTONE, not a proportionality |
| R-14 | R | s_nf dilation (ideal bond) at CLOUD basis | R (P1.3) | CLOUD n=126 | +1.0425 | 2.79x | WORSE, 5/5 folds -- P1.3 FALSIFIED |
| R-15 | R | Rg-matched dilation at CLOUD basis | **E** | CLOUD n=126 | +0.0622 | 1.47x | WORSE, 5/5 folds |
| P-31 | cos(μ̂, μ) for μ̂ = c − t′ (a MISMATCHED same-length deposited native — native-free at inference) | cloud, ORACLE-scored | EXPLORATORY | **−0.0088 ± 0.0235**, 52% positive — zero |
| P-32 | shared-referent FLOOR for P-31: the pure scale direction c−mean(c) | cloud, ORACLE-scored | EXPLORATORY (control, measured first) | +0.0512 ± 0.0342 — the floor is also ~zero, so P-31 is not floor-inflated, it is empty |
| P-33 | cos(μ̂, μ) for μ̂ = c − mean(same-length deposited natives) | cloud, ORACLE-scored | EXPLORATORY | −0.0080 ± 0.0419 |

### Registered — D4-M: the MARGINAL, which is the baseline a one-bit classifier must beat

A constant sign predictor gets `max(p, 1−p)` where `p` is the bit's marginal. **Reported first, as it
must be** — D-46…D-49's accuracies mean nothing without it. `s32/results/s32_D5_signprice.json`.

| # | reg? | scorer | **marginal** (constant predictor) | leave-fold-out accuracy | beats it? | per-fold LFO accuracy |
|---|---|---|---|---|---|---|
| D-50 | R | **AMBER** | 0.508 | **0.381** | **no — worse than constant** | 0.44 0.39 0.36 0.30 0.40 |
| D-51 | R | `DIS` | 0.548 | 0.548 | no — **exactly** the marginal | 0.64 0.39 0.60 0.52 0.57 |
| D-52 | R | `LEG_total` | 0.587 | 0.595 | **+0.008 = one target in 126** | 0.56 0.52 0.64 0.52 0.70 |
| D-53 | R | `LEG_torsion` | 0.556 | 0.556 | no — **exactly** the marginal | 0.60 0.70 0.44 0.43 0.60 |

**No scorer's sign is predictable above its own marginal.** The coordinator's registered prediction
was 55–65% accuracy; the accuracies *are* in that range and it means nothing, because **the marginal
is already there**. Per-fold accuracy swings 0.30–0.44 (AMBER) and 0.52–0.70 (`LEG_total`) — the
`LEG_total` mean is one fold's 0.70 against two folds at 0.52. *This is exactly the case the
coordinator flagged: "if the bit is 70/30, a 70%-accurate classifier is nothing."*

### Registered — D5-C: the Ångström price of one ORACLE bit per target. **CLOUD BASIS — a SCREEN.**

`s32/results/s32_D5_signprice.json`. **CA point cloud, 3.0483 reproduced exactly on PROD** (an
instrument check). This basis is a screen for whether the bit is worth projecting; **it is never
differenced against a built-chain number.** Arm = uniform coordinate average of the best `m` of the
shipped top-75 by the (signed) score. All ORACLE rows **ORACLE / NOT DEPLOYABLE**.

| # | reg? | arm (CLOUD) | mean | Δ vs PROD | se | ×MDE | folds |
|---|---|---|---|---|---|---|---|
| D-54 | R | `LEG_total` m=25 **ORACLE sign** | **2.9196** | **−0.1287** | 0.0303 | 1.52 | 4/5 |
| D-55 | R | `LEG_total` m=25 **const +1** (zero-info control) | 3.1399 | **+0.0916** | 0.0328 | 1.00 | 5/5 |
| D-56 | R | `DIS` m=25 ORACLE | 2.9308 | −0.1176 | 0.0229 | 1.83 | 5/5 |
| D-57 | R | `DIS` m=25 const | 3.0749 | +0.0266 | 0.0269 | 0.35 | 3/5 |
| D-58 | R | **AMBER** m=25 ORACLE | 2.9602 | −0.0881 | 0.0194 | 1.62 | 5/5 |
| D-59 | R | AMBER m=25 const | 3.1377 | +0.0894 | 0.0263 | 1.21 | 5/5 |
| D-60..63 | R | the four m=50 arms | 2.9740–3.0866 | −0.0743…+0.0383 | — | 0.67–1.89 | 4–5/5 |

**Every zero-information (`const +1`) arm is WORSE than the production top-75 average** — selecting
the best 25 by an unsigned scorer costs +0.027 to +0.092 Å. This is lane V's *"the score prefix is
worse than random at retaining the best candidate"* arriving independently through a different
operator. **The ORACLE bit turns a loss into a gain**, and the paired price of the bit itself is
measured on the built chain by `s32_D5_signchain.py` (running; the arms are projected in one process
so the pairing is within-job).

**Order-of-magnitude context, and it is the honest framing:** S31 §20.1 prices the ORACLE common-mode
direction at **−0.8102 Å on the BUILT CHAIN**. One oracle bit per target is worth **~0.13 Å on the
CLOUD**. *The bit is real and it is roughly an order of magnitude smaller than the direction* — which
is charter §41's question answered in the currency it asked for: **the five-bit result is a pricing
clue, not a limit, and per-bit the price is small.**
| P-34 | R2/R2b/R1 BEST-axis deficits, stratified FAIL18 vs other-108 | cloud, ORACLE | REGISTERED (H-P5) | 500→128: mean +0.1853, **median −0.0371**, 72 of 126 favour the score; FAIL18 +1.4773 = **114% of the total**, other-108 −0.0300 at 0.30× — **NOT A RESULT off the circular stratum** |
| P-35 | same for the DEPLOYABLE 500→75 filter | cloud, ORACLE | REGISTERED (H-P5) | mean +0.2260, median +0.0046, 62W/64L; FAIL18 +1.7713 = 112%, other-108 −0.0316 at 0.27× |
| P-36 | same for universe→500 BLOSUM | cloud, ORACLE | REGISTERED (H-P5) | mean −0.0718; FAIL18 contributes **0%**, other-108 −0.0833 at 0.98× — NOT MEASURED, and a different SHAPE from the score's |

### Exploratory but mechanism-motivated — D4-X: `sign(rg_pred − rg_pool)` as the bit. **Also fails, and the reason is that the feature is nearly constant.**

D-20 showed `LEG_total`'s in-band sign agrees with `Rg`'s on 68.3%. That suggests a *physical*
native-free candidate for the bit: **if the predicted structure is more compact than the pool, then
more-compact candidates are better, so a compactness-tracking scorer should be oriented `+`.** The
distogram supplies `rg_pred` without touching a native. One feature, one test, motivated by a
measured mechanism — logged as exploratory. `s32/results/s32_D4X_rgsign.json`.

| # | reg? | scorer | marginal | accuracy of `sign(rg_disagree)` | oriented ρ | vs its own constant control | ×MDE |
|---|---|---|---|---|---|---|---|
| D-64 | E | `RG` | 0.540 | 0.587 | +0.0879 | +0.0378 | 0.33 |
| D-65 | E | `LEG_total` | 0.587 | **0.540 — below its own marginal** | +0.0410 | +0.0034 | 0.04 |
| D-66 | E | **AMBER** | 0.508 | **0.460 — below chance** | −0.0256 | −0.0256 | 0.64 |

**The mechanism for the failure, and it is the feature's own distribution:** `rg_pred − rg_pool_mean`
is **positive on 81.0% of targets**. `sign(rg_disagree)` is therefore *nearly a constant `+1`*, and
its apparent accuracy is its target's marginal wearing a different label. Its correlation with the
ORACLE quantity it is meant to predict is **+0.199**. **A one-bit feature that is 81/19 cannot carry
a 59/41 label.** Every arm is under 0.7× MDE — NOT A RESULT.

## LANE L (length generalisation) — prereg `s32/PREREG_S32_L.md` @ 88f2da39

| # | comparison | type | basis | verdict |
|---|---|---|---|---|
| L-1 | PROJ_NAT_LONG (n=60, L 41–60) vs PROJ_NAT_SHORT (n=126, L 9–16) — deployed projector on the native, λ=0 | **registered** (L-H1) | CA cloud vs native, ORACLE / NOT DEPLOYABLE | 0.700 (SE 0.029) vs 0.0429 (SE 0.0053). **L-H1 FALSIFIED**: the registered falsifier was "< 1.0 Å at L=45"; measured mean 0.700, p90 0.956, max 1.129. The representation is NOT the obstacle at 40–60 residues. |
| L-2 | native-torsion rebuild vs deployed-projector residual, same 60 chains | registered control (L-H1) | CA cloud, ORACLE | 2.742 vs 0.700 — the rebuild bound overstates the manifold distance by **3.9×** at L≈52 and by **8.1×** on the canonical 126 (0.347 vs 0.0429). |
| L-3 | leakage filter `identity(norm="shorter")` vs composition-shuffled null, 60 targets × 3 shuffles | **registered control** (instrument admission) | rejection rate | **At the null**: rejects 100% of real AND 100% of shuffled at every threshold ≤ 0.9. Filter discarded. `verbatim` separates (real 0.330, shuffled 0.000) and is used instead. |
| L-3b | common-mode fraction f, BLOSUM top-75, L 40-60 (n=45) vs L 9-16 (n=126) | **exploratory** (L-H3a) | squared CA error in the medoid frame; bias-variance identity verified to 1.6e-15 | 0.530 (SE 0.027) vs 0.468 (SE 0.016); diff +0.062, SE_diff 0.031, **0.70x MDE -> NOT MEASURED**. The fraction does not measurably move over a 4x length change. **DEFINITION WARNING: these are BLOSUM top-75 subsets. The 0.676 on record is the DISTOGRAM top-75 and is a different object; no claim is made about it.** |
| L-3c | provenance audit of the ORACLE pool-best member, long40 (n=45) | **registered adversarial check**, threshold fixed before reading | sequence identity | Window identity to target: mean 0.231, **median 0.189**, p90 0.400, max 0.519. Longest verbatim common substring: median 3, max 7. Withdrawal threshold was median > 0.40. **HOLDS** -- the long pool's headroom is retrieval reach, not PDB redundancy. ORACLE pool-best cloud RMSD over all 45: mean 4.980, median 5.475. |
| L-3d | in-band Spearman(BLOSUM sim, ORACLE CA-RMSD) inside the K=500 pool, L~55 (n=45) vs L~13 (n=126) | **EXPLORATORY**, found by looking at the ladder's filter_skill column | rank correlation, ORACLE labels | **HYPOTHESIS NOT SUPPORTED.** rho -0.0464 (SE 0.0074) at L~13 vs **+0.0261** (SE 0.0304) at L~55 -- the WRONG SIGN at length, 44.4% of targets negative vs 67.5%. Contrast +0.072, MDE 0.088, **0.83x -> NOT MEASURED**. The global in-band sequence channel does NOT strengthen with length; `structure-and-sequence-are-decoupled` survives at 55 residues as a statement about GLOBAL ordering. What does grow is the HEAD-of-order advantage (top-1 vs pool mean: 8% at L~13, 25% at L~55) -- a different statistic, reported as such. |
| L-3e | readout prefix size: ONE global m fitted OUT-OF-FOLD vs the shipped m=75, at L~55 (n=45) and L~13 (n=126) | **EXPLORATORY**, and the only lane-L arm proposing a deployable change | **CA CLOUD, not the built chain** | **L~13: +0.0201, MDE 0.0424, 0.47x -> NOT A RESULT.** m=75 is right at peptide length (curve flat 3.304/3.293/3.293 over m=30/50/75). **L~55: -1.0127, MDE 1.4156, 0.72x -> NOT MEASURED**, W/L 21/24, per-fold +0.200/-1.111/-0.540/-1.265/-2.347 (4/5 improve). The cloud curve minimises at **m=3** (8.130) and every fold's out-of-fold choice is m in {3,5}; the shipped m=75 sits 1.28 A past the minimum. A LEAD, not a result: near-even W/L with a large mean is the concentration warning, and it is on the CLOUD. |
| L-3f | ORACLE per-target m vs LFO global m, L~55 | ORACLE, priced as an order statistic | CA cloud | 6.557 vs 8.397. The per-target argmin over a 13-point grid is mostly best-of-13 (memory: `grid-oracles-are-order-statistics`) and is quoted only to bound the arm. |
| L-5 | CLOUD ladder, FULL instrument, both lengths (n=126 and n=45) | **registered** (L-H2, P1/P2) | **CA point cloud** -- diagnostic, NOT the endpoint | See the table below. **P1 and P2 both HOLD at full instrument size, every rung a RESULT.** Implementation check: the recomputed canonical `pool_best` is **1.7108**, bit-matching the value pinned in `s12/instrument.py`'s selfcheck. |
| L-4 | `long40` ladder: pool_best / sparse_s10 / top75_best / avg75 / avg75_random, at L 40–60 and L 9–16 | **registered** (L-H2, P1–P3) | built chain, pre-AMBER | in progress |

| P-37 | rank of the pool's deviation matrix vs 3n−6 | cloud, native-free | EXPLORATORY (structural) | **exact on 126/126**; μ lies inside that span on 126/126, which is why μ is exactly recoverable |

### Registered — D3-M: PHYSICS AS MOVER. Interim read at n = 72 (full n = 126 in the ledger entry).

The lane's primary registered arm, and the only escape from G1 that survived D0 (Corollary D-G:
physics is closed as a ranker and open as a mover). AMBER relaxes the production built chain; the
move is converted to pair-distance space, `Δd = d(x_relaxed) − d(x_production)`, and contrasted with
the ORACLE error vectors. **Endpoint arms are BUILT CHAIN, projected in the same process as their
paired production rows.** Registered prediction: **median `cos(Δd, −e_prod)` ≥ +0.10.**

Interim, n = 72 of 126 — final numbers, W/L and fold CIs go in the ledger:

| rung | `cos(Δd, −e_prod)` (ORACLE) | ×MDE | magnitude-matched random control | ENDPOINT Δ (built chain) | ×MDE | folds |
|---|---|---|---|---|---|---|
| k=100 | −0.0018 | 0.02 | −0.0019 | **+0.0119** | 2.46 | 5/5 |
| k=10 (production) | +0.0180 | 0.21 | −0.0052 | **+0.0310** | 1.80 | 5/5 |
| k=1 | +0.0190 | 0.19 | −0.0000 | **+0.0758** | 1.78 | 5/5 |
| k=0 (free) | +0.0247 | 0.23 | −0.0010 | **+0.1045** | 1.89 | 5/5 |

**Registered prediction FAILED, and the falsifier fired on its own terms.** `cos` is 0.02–0.23× MDE
and **indistinguishable from a magnitude-matched random direction in the same space**. The scale
ladder does not rescue it: **every α ∈ {0.25, 0.50, 0.75, 1.0} is positive (worse)** at every rung,
best case α=0.25 at +0.0027. Cα displacement from production: 0.083 Å (k=100) → 0.627 Å (k=0), and
the endpoint damage grows monotonically with it.

**Mechanism: the physical relaxation displacement is ORTHOGONAL to the production structure's error.**
Physics-as-mover cannot reach the prize however it is scaled, because it is not pointing at it —
and applying it costs **+0.012 to +0.105 Å** at the endpoint.

**A cross-check that came free and is worth more than the arm.** `cos(e_prod, e_pool75) = +0.9519`:
the production chain's own pair-distance error and the pool's common-mode error point in **almost
exactly the same direction**. S31 §20.1's common-mode story, confirmed from a different object
(the production chain's error, not the prior's) in a different lane.
| P-38 | the whole H-P3 decomposition repeated on the **production top-75** set (not the quantum 128) | cloud, ORACLE | REGISTERED (H-P3, basis check) | identity 4.6e-13; f_common **0.6780** (matches S23's 0.676); Var(U)/Var(V) median 3.45 / 4.73 in the score's band; ρ_inband(CONS) **−0.2362 ± 0.0145** ORACLE band, **+0.2829 ± 0.0322** score band; zero-μ +0.7828; μ recovery 1e-14 |

### Registered BEFORE the endpoint rows were complete (n=123 of 126 at the time of writing, chain values NOT read)

The operator law fitted in `s32_P_decomp.json` on the SCORE and ORACLE prefix families
(`out = −0.2536 + 0.6730·set_mean + 0.3256·set_best`, cloud, R² 0.9467) is used here to predict a
family it was **not** fitted on. Set statistics of the random-gate arms (ORACLE labels, but the
chain numbers were not looked at):

```
          set_mean   set_best    wBw
PROD        3.5520     2.3056    127.8
R128        3.8576     2.1165    193.4      Δ mean +0.3056  Δ best −0.1891
R75         4.4541     2.0520    316.2      Δ mean +0.9021  Δ best −0.2536
```

| # | PREDICTION (registered) |
|---|---|
| P-39 | `R128 − PROD` cloud effect = **+0.144 Å** (0.673·(+0.3056) + 0.326·(−0.1891)) |
| P-40 | `R75 − PROD` cloud effect = **+0.524 Å** |

Both predict the random gate is **WORSE**, agreeing with the coordinator's own registered
prediction. A chain effect within ±0.10 Å of these confirms the law transfers to a gate family it
was not fitted on; a miss falsifies it for random gates.

### The random-gate endpoint arm — BUILT CHAIN, n = 126, all arms projected in the same process per target

| # | comparison | registered? | emitted |
|---|---|---|---|
| P-41 | `R128` (random 128 of 500, then the score's top-75 within) − PROD | REGISTERED (H-P2b) | **+0.1648, 0.92× MDE — NOT MEASURED**; median +0.0734, 52W/74L, 4/5; **0 of 8 draws better**, draw sd 0.0187 |
| P-42 | `R128` − PROD on the **other 108** (filter-independent) | REGISTERED | **+0.2614, 1.52× — WORSE**, 37W/71L |
| P-43 | `R128` − PROD on **FAIL18 (CIRCULAR)** | REGISTERED, diagnostic only | −0.4149, 0.70× — NOT MEASURED, 15W/3L |
| P-44 | `R75` (random 75 of 500, score never consulted) − PROD | REGISTERED (H-P2b) | **+0.4086, 1.65× MDE, 5/5 folds, 43W/83L — WORSE**; 0 of 8 draws better, draw sd 0.0381 |
| P-45 | `R75` − PROD on the other 108 | REGISTERED | **+0.5627, 2.23× — WORSE** |
| P-46 | `R75` − PROD on FAIL18 (CIRCULAR) | REGISTERED, diagnostic only | −0.5160, 0.94× — NOT MEASURED, 14W/4L |
| P-47 | operator law on the **BUILT CHAIN**, all 17 arms | REGISTERED (H-P2b-1) | `out = −0.9934 + 0.9219·set_mean + 0.3232·set_best`, **ratio 2.85**, R² 0.9162 — the registered ≥5× **FAILS**, and S18's 1.16/0.04 (ratio 29) does not hold under gates that order candidates |
| P-39 | scored: predicted `R128 − PROD` cloud **+0.144**, measured **+0.1369** | REGISTERED | **HIT** (error 0.007) |
| P-40 | scored: predicted `R75 − PROD` cloud **+0.524**, measured **+0.3740** | REGISTERED | **MISS** (error 0.150) — the law over-predicts a large gate change; it is a local linearisation |
| P-48 | PROD reproduction against `s29_O_chain_rows` | REGISTERED (rule 3 check) | mean \|Δ\| 0.0123 (within the 0.0134 floor), p90 0.0260, **max 0.5174 on 2LNG — ABOVE the contract's 0.2285 max floor** |

---

## D1-T RE-EMITTED, and the old artefact SUPERSEDED — `s32/s32_D1_signtransfer.py` → `s32_D1_signtransfer_v2.json`

**A reproducibility defect of mine, recorded in place rather than by deletion (contract rule 13).**
D1-T's first four numbers (D-36…D-39) were produced by an **inline `python -c` one-liner**. The
resulting `s32/results/s32_D1_signtransfer.json` carries **no provenance block, no module, no git
commit, no source hash and no pinned seed, and no script in the repository produced it** — while two
downstream scripts quoted its numbers as a ceiling and the ledger quoted them as a headline. Charter
§61 requires verifiable experiments, artefacts and seeds; it failed all three. **Project memory
records this exact shape twice already (`findings-prose-is-not-evidence-of-code`); this is the
third.** The old artefact is **SUPERSEDED and must not be quoted.** D-36…D-39 stand as originally
written, annotated here, and the table below replaces them.

### Three nulls, not one — and the first one could not test the claim

`nullPERM` (mine) permutes `rr` inside half B, destroying **all** structure, so it sits at ~0 for
every scorer **including pure noise**. It answers *"is there any relation?"*, not *"is the
relation's SIGN a property of the target?"* **`nullXTGT` — another target's sign applied to this
target's held-out half, 200 draws, own distribution — is the null the claim needs**, and
`globalSGN` (one leave-one-fold-out majority sign for all test targets) prices the global
alternative. `NOISE` is the self-test: **an instrument that scored noise as PER-TARGET would be
broken** (contract rule 5).

**DEDUPLICATED numbers are the ones to quote.** 7.7% of band members are exact coordinate
duplicates (**122/126 targets, mean 69.2 distinct of 75**); a duplicate in both halves buys sign
agreement for free.

| # | reg? | scorer | transfer (DEDUP) | nullPERM | **nullXTGT** | globalSGN | ×MDE vs nullXTGT | folds | verdict |
|---|---|---|---|---|---|---|---|---|---|
| D-67 | R | **AMBER** | **+0.1163** | +0.0058 | +0.0004 | −0.0725 | **2.61** | 5/5 | **PER-TARGET** |
| D-68 | R | `DIS` | +0.1901 | +0.0010 | +0.0088 | +0.0649 | 3.19 | 5/5 | PER-TARGET |
| D-69 | R | `LEG_total` | +0.2180 | +0.0042 | +0.0082 | +0.0426 | 3.10 | 5/5 | PER-TARGET |
| D-70 | R | `LEG_torsion` | +0.1492 | +0.0046 | +0.0017 | +0.0460 | 2.68 | 5/5 | PER-TARGET |
| D-71 | R | **`RG`** — one line of numpy, not a Hamiltonian | **+0.3228** | +0.0014 | +0.0053 | +0.0489 | **4.15** | 5/5 | PER-TARGET |
| D-72 | R | **`NOISE`** — the self-test | **+0.0056** | −0.0004 | +0.0016 | +0.0123 | **0.15** | 2/5 | **NOT A RESULT** (correct) |

Raw (non-dedup) numbers move by at most 0.012 and no verdict changes, so **duplicates were not
driving it.** `DIS` at **+0.1901** replicates lane V's independently coded **+0.1911** to 0.001.

### Two corrections to how this must be framed

**(i) They are not unrelated Hamiltonians — they all load on compactness, and plain `Rg` loads
hardest.** `RG` (+0.3228) beats every Hamiltonian measured: `LEG_total` +0.2180, `DIS` +0.1901,
`LEG_torsion` +0.1492, AMBER +0.1163. My own D-20 already said it — `LEG_total`'s in-band sign
agrees with `Rg`'s on **68.3%**. ***The bit is very likely "is the native more or less compact than
its own band", and the Hamiltonians are reading it through their compactness terms.***

**(ii) This is largely a CONFIRMATION on a new instrument, not a discovery.** Memory
`in-band-ordering-is-per-target` already records *"the only leverage supplies the per-target SIGN at
inference — and native-free compactness proxies reach 0.24–0.37 (all CIs exclude zero) against the
oracle's 0.909."* **`RG` at 0.3228 lands inside that interval.** What is new here is the
**cross-target-null verdict on four physical scorers including AMBER**, the **NOISE self-test**, and
the ceiling below.

### THE CEILING, which is the lane's real deliverable

Memory `in-band-ordering-is-per-target` prices **2.0 Å at in-band ρ ≈ 0.638**. A **perfect, free**
per-target sign oracle delivers:

```
AMBER        +0.1163   = 18% of the 0.638 crossing price
LEG_torsion  +0.1492   = 23%
DIS          +0.1901   = 30%
LEG_total    +0.2180   = 34%
RG           +0.3228   = 51%
NOISE        +0.0056   =  1%      (the self-test)
```

> ### Even a free, perfect per-target sign leaves the best in-band scorer 2–3x short of the useful range.
> ### The sign is closed as a ROUTE and kept as a FINDING.

### Registered — D5-B: the Ångström price of the bit on the **BUILT CHAIN**, n = 126

`s32/s32_D5_signchain.py` → `s32_D5_signchain_LEG_total.json`. All three arms projected **in one
process**, so the pairing is within-job (contract rule 3). `PROD` rebuilds to **3.2126** against the
canonical **3.2105** — a **0.0021** gap, inside the 0.0107 mean projection floor. Geometry
secondaries on row one: **virtual-bond mean 3.80395 and sd about 1e-15 on all three arms** — every
arm is a valid ideal-geometry chain.

| # | reg? | arm | mean (chain) | effect | median | se | ×MDE | folds | W/L |
|---|---|---|---|---|---|---|---|---|---|
| D-73 | R | **ORACLE sign**, m=25, vs PROD | 3.0921 | **−0.1206** | −0.0485 | 0.0322 | 1.34 | 4/5 | 86/40 |
| D-74 | R | **const +1** (zero-info control), m=25, vs PROD | 3.3022 | **+0.0896** | +0.0270 | 0.0355 | 0.90 | 4/5 | 49/77 |
| D-75 | R | **PRICE OF ONE ORACLE BIT** = ORACLE vs its own matched const | — | **−0.2101** | **+0.0000** | 0.0453 | **1.66** | **5/5** | **43/9** |

**ORACLE / NOT DEPLOYABLE.** The median is **exactly +0.0000 and W/L is 43/9 with 74 ties — and that
is structural, not concentration**: the bit's marginal is 0.587, so on the 58.7% of targets where the
sign is `+1` the ORACLE arm *is* the const arm and the difference is identically zero. **The whole
effect is carried by the 52 targets where the bit is `−1`, on which it is worth about −0.51 Å.**
(Memory `median-vs-mean-is-the-free-warning` asks for this to be checked; here the explanation is
the marginal, and it is checkable from the artefact's `bit_marginal_frac_positive = 0.587`.)

**Order-of-magnitude context, same basis:** S31 §20.1 prices the ORACLE common-mode **direction** at
**−0.8102 Å on the built chain**. One ORACLE **bit** per target is **−0.2101 Å on the built chain**
— **26% of the direction**, and, per D4-S / D4-M, supplied by nothing native-free.

### Closing notes for LANE P

**The duplicate-row class, resolved.** `s32_P_rand_rows.s*of3.jsonl` reached 150 raw lines over 126
distinct pdbs because two launchers raced on shards 1 and 2 (the coordinator killed the surplus
*wrappers*; lane P then found and killed their orphaned *workers*, pids 25196 and 25960 — killing a
`jobrun` wrapper does **not** kill its child). `_rows()` in `s32_P_rand.py` deduped by pdb, asserted
126 rows **and** 126 distinct pdbs, and checked every duplicate pair for **differing content** (0
conflicting over all 150 lines) before any number was computed, so **no reported figure was ever
computed over 150 rows.** The shard files were then deduped in place — lossless, because the pairs
are byte-identical — and the analysis re-run: `PROD 3.2126`, `R128 +0.1648 (0.92×)`,
`R75 +0.4086 (1.65×)`, **identical to four decimals.** `s32_verify.py` audit 8 no longer flags them.

**Is lane P's `μ`-sign bit the same object as lane D's Rg sign?** Evidence from lane P's own
artefact, offered so two vocabularies are not mistaken for two findings:
- lane P's `RG_probe` sign rule — a two-sided compactness-vs-predicted-Rg probe along `pc1` —
  predicts lane P's `sign(⟨pc1, μ⟩)` bit at **1.32× MDE**, above chance. The compactness axis
  carries information about the `pc1` sign, so they are **not independent**.
- both carry the identical signature: substantial |cos| (`pc1` cos² 0.2277, `RADIAL` cos² 0.1339,
  against a measured same-space random null of 0.0316) with **mean cos indistinguishable from zero**
  — the magnitude is generic, the per-target sign is the missing thing. That is S14's
  `in-band-ordering-is-per-target` in a third vocabulary.
- **Lane P prices the same bit in Ångströms rather than in ρ:** one *perfect* bit per target, with a
  single global leave-fold-out step, is worth **+0.1835 Å** of the 3.0776 Å common mode. That is the
  same conclusion as lane D's "a free perfect sign still leaves you 2–3× short of ρ = 0.638": **a
  real finding and a closed route.**

### LANE L — SELF-DECLARED DEFECT (found by lane L, before any ladder number was quoted)

`s32_L_ladder.run_target` seeds the `avg75_random` zero-information control with
`RNG_SEED + abs(hash(pdb)) % 10000`. **Python's string `hash()` is randomised per process**
(`PYTHONHASHSEED` is unset here; two interpreters return 3567782200507830974 and
4032705310864098348 for the same string). Therefore:

- the control's draws are **not bit-reproducible** across runs, and differ between the two
  worker processes within the same run;
- the drawn subsets are still uniform and the arm is **unbiased**; the reported quantity is
  the **draw mean over N_DRAW = 3** with the **draw-to-draw sd** beside it, as contract
  rule 10 requires, so no conclusion rests on a single lucky draw;
- but a re-run will not reproduce the control column exactly, and the `L2_ladder_*.jsonl`
  artefacts are the only record of the draws actually used (each row carries its three
  per-draw values under `draws`, so the arm can be audited even though it cannot be
  replayed).

Every other lane-L arm is deterministic. The fix (a stable digest of the pdb code in place
of `hash`) is applied to the module **after** this run so that code and artefact are never
silently inconsistent; the artefact predating the fix is flagged here.
| R-16 | R | GEN4D named-start choice: ORACLE global + LEAVE-FOLD-OUT vs PROD | **E** | built chain n=126 | LFO +0.0081 | 0.28x, 2/5 folds | NOT MEASURED. The only family that transfers (60%) and it is redundant: production's argmin already reaches the alpha-helix start |
| R-17 | R | MEDOID_EXTRA / SCALE_NF / SCALE_NF_MED vs PROD (CHAIN basis) | R (R4) | built chain n=126 | -0.0064 / +0.7219 / +0.7228 | 0.35 / 2.24 / 2.24 | MEDOID_EXTRA's fold CI excludes zero and 4/5 agree, but the MDE gate binds |
| R-18 | R | SCALE_GRID per-target ORACLE, priced best-of-K + LEAVE-FOLD-OUT | R (R1) | built chain n=126 | ORACLE -0.1352; LFO +0.0044 | 3.84x ORACLE; 0.14x LFO | 122% accounted by the across-target null. Its split-half (-0.0497) is centred on the GRID MEAN and must NOT be quoted vs production |

### Registered — D3-M FINAL, n = 126. **The registered prediction FAILED; the endpoint is a regression.**

Supersedes the n = 72 interim above. `s32/results/s32_D3_mover_A_0_1.jsonl`, `s32_D_analyse.json`.
**BUILT CHAIN, all arms produced in one process from the production chain they are paired against.**
`PROD` = **3.2126** (canonical 3.2105; gap 0.0021, inside the 0.0107 mean floor) and reproduced to
the same 3.2126 independently in `s32_D5_signchain`.

| # | reg? | rung | `cos(Δd,−e_prod)` ORACLE | ×MDE | matched-random ctrl | **surplus over ctrl** | ×MDE | ENDPOINT Δ (chain) | ×MDE | folds | W/L |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D-76 | R | k=100 | +0.0057 | 0.09 | −0.0019 | +0.0076 | 0.12 | **+0.0204** | 0.84 | 5/5 | 23/103 |
| D-77 | R | k=10 (production) | +0.0412 | 0.63 | −0.0036 | +0.0448 | 0.69 | **+0.0336** | 1.26 | 5/5 | 40/86 |
| D-78 | R | k=1 | +0.0555 | 0.69 | +0.0002 | +0.0553 | 0.69 | **+0.0726** | 1.96 | 5/5 | 33/93 |
| D-79 | R | k=0 (free) | +0.0651 | 0.73 | +0.0014 | +0.0637 | 0.71 | **+0.0985** | 2.21 | 5/5 | 35/91 |
| D-80..91 | R | the 12 scale-ladder arms, α ∈ {0.25, 0.50, 0.75} × 4 rungs | — | — | — | — | — | **+0.0039 … +0.0601, ALL POSITIVE** | — | — | — |

**Registered prediction (median `cos ≥ +0.10`): FAILED.** Against its own magnitude-matched random
direction the surplus is **0.12–0.71× MDE — NOT A RESULT at three rungs, NOT MEASURED at the
fourth.** *The relaxation displacement is indistinguishable from a random direction of the same
length.* The endpoint regression grows monotonically with the freedom given to the physics, 5/5
folds at every rung. **The scale ladder does not rescue it** — every α at every rung is worse, best
case **+0.0039 Å**. *A direction that is not pointing at the target cannot be fixed by shortening
the step.*

**Geometry secondaries, row one (contract rule 15).** Production's virtual bond is **3.80395 Å,
sd 9.2e-16** — an exact ideal-geometry chain. Relaxed: **3.84913 Å sd 2.78e-02** (k=100),
**3.86914 Å sd 4.66e-02** (k=0). **A systematic +1.2% to +1.7% Cα–Cα dilation with real scatter,
growing with the same knob as the damage.** A dilation-only control was **not** run, so this is a
named candidate mechanism, **not a demonstrated one**.

**Free by-product, worth more than the arm:** `cos(e_prod, e_pool75) = +0.9443` and
`cos(e_prod, e_disto) = +0.6621`. **The production chain's own pair-distance error is 94% aligned
with the pool's common-mode error** — S31 §20.1's account confirmed from a different object, in a
job that was not looking for it.

### Lane D — closing tally

**92 comparisons emitted** (D-1…D-91 plus the 3 exploratory Legacy per-term rows), of which
**registered: D-1…D-17, D-30…D-59, D-64 not, D-67…D-91**; exploratory: D-8, D-9, D-18…D-29, D-35,
D-64…D-66. **Registered predictions that FIRED AGAINST THE LANE: three** — the `<15%` chiral in-band
variance share (D-12, measured 23.6%), the H-D2 basin horns (D-40/D-41, neither fired at 0.972 and
0.5006 Å), and the D3-M `cos ≥ +0.10` (D-76…D-79, measured 0.09–0.73× MDE). **One reproducibility
defect of my own, found by lane V and fixed in place** (D1-T's missing script and provenance).
**Zero lane-D arms reached the endpoint as an improvement**, and the one ORACLE arm that did
(D-75, −0.2101 Å) is **ORACLE / NOT DEPLOYABLE** with nothing native-free to supply it.
| R-19 | R | typicality ARGMIN vs ARGMAX (both directions), vs PROD | **E** | built chain n=126 | ARGMIN -0.0054 / ARGMAX +0.3314 | 0.29x / 2.58x | cross-lane sign defect resolved; the MOST typical branch is a null, the LEAST typical is much worse |

### L-5 — the cloud ladder at full instrument size (basis: CA POINT CLOUD, not the endpoint)

| rung | L ~ 13 (n=126) | L ~ 55 (n=45) |
|---|---|---|
| `pool_best` ORACLE / NOT DEPLOYABLE | 1.7108 (SE 0.0791) | 4.9797 (SE 0.4329) |
| `top75_best` ORACLE / NOT DEPLOYABLE | 2.1041 (SE 0.0994) | 5.8330 (SE 0.5055) |
| `avg75` deployable-shaped | 3.2928 (SE 0.1390) | 9.4094 (SE 0.3403) |
| `pool_mean` (a typical member) | 4.4533 | 11.9866 |

| difference (LOWER IS BETTER) | L ~ 13 | L ~ 55 |
|---|---|---|
| headroom `avg75 - pool_best` | **+1.5820**, MDE 0.3045, **5.20x**, W/L 1/125 | **+4.4297**, MDE 1.3164, **3.36x**, W/L 2/43 |
| retrieval `top75_best - pool_best` | +0.3933, MDE 0.1207, 3.26x, W/L 0/102 | +0.8534, MDE 0.4815, 1.77x, W/L 0/25 |
| readout `avg75 - top75_best` | **+1.1887**, MDE 0.2777, **4.28x**, W/L 6/120 | **+3.5764**, MDE 1.4805, **2.42x**, W/L 8/37 |

**Share of the recoverable loss:** retrieval 24.9% / readout 75.1% at L~13; retrieval 19.3%
/ readout **80.7%** at L~55. The ladder's shape is preserved over a 4x length change and the
readout's share GROWS. `top75_best` here is the **BLOSUM** top-75, not the distogram top-75
(2.3062 on record) — a different object, never differenced against it.
| R-20 | R | rama headroom: is there plausibility variation left among lam=0.3 branches? | **E** | diagnostic, n=126 | posphi within-target spread 0.632; 31.2% of branches above the 17.5% unconstrained rate; 126/126 targets carry one | n/a | EXPLORATORY; REFUTED my own first mechanism for R3's negative |
