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

| | count |
|---|---|
| hypotheses opened | 0 |
| pre-registered arms (REG) | 0 |
| exploratory arms (EXP) | 2 |
| confirmatory / replication arms (CONF) | 0 |
| audit arms (AUDIT, no hypothesis, no multiplicity cost) | 4 |
| **comparisons emitted against the endpoint** | **0** |
| headline claims | 0 |

No arm has yet emitted a comparison against the 3.2105 endpoint. The multiplicity budget is
untouched. The moment a lane emits one, it lands here in the same edit as its number.

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
| V-A2 | Is 3.2105 bit-reproducible from its stated input? | `s32/results/s32_V_chain_bitexact.json` | see **D1** |
| V-A3 | Size-matched random null for the ladder's narrowing increments (2000 draws/target) | `s32/results/s32_V_ladder_orderstat.json` | see **D2** |
| V-A4 | rr spread of the DIS top-128 vs the pool (mechanism for V-A3) | `s32/results/s32_V_top128_spread.json` | mean 4.4533 -> 3.5847, sd 1.3206 -> 0.6384, **p5 2.6098 -> 2.6184 (unchanged)** |

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
| R-1 | R | PROD(cache cloud) vs PROD(s29 cloud) | R | built chain n=126 | pending | pending | the 1e-14 cloud floor |
| R-2 | R | lam=0.3 arm vs lam=0 arm (same cloud, same job) | R | built chain n=126 | pending | pending | the ladder's own cost |
| R-3 | R | PROD chain vs isotropic null sqrt(e^2+d^2) | R | built chain n=126 | pending | pending | R1 direction test |
| R-4 | R | ORACLE best branch (9 family subsets) vs PROD | R | built chain n=126 | pending | pending | ORACLE ceiling |
| R-5 | R | order-statistic curve, 9 subset sizes x 24 draws | R | built chain n=126 | pending | pending | best-of-K pricing |
| R-6 | R | in-band rho, 11 criteria | R | within-target rank | pending | pending | R3a |
| R-7 | R | SEL_<criterion>_<subset>, 11 criteria x 2 subsets | R | built chain n=126 | pending | pending | R3b directional |
| R-8 | R | BRANCHMEAN, RANDBRANCH (5 draws), MEDOID_ONLY, OBJARGMIN x5 | R | built chain n=126 | pending | pending | controls + R4 |
| R-9 | R | MEDOID_EXTRA, SCALE_NF, SCALE_NF_MED, SCALE_GRID x8 | R | built chain n=126 | pending | pending | R1/R4 repair job |

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
| R-10 | R | ladder rungs x5: observed price vs ORTHOGONAL isotropic null | R | built chain n=126 | pending | pending | R1 registered test |
| R-11 | R | ladder rungs x5: cos_align + on-manifold null sqrt(e^2-d^2) | **E** | built chain n=126 | pending | pending | EXPLORATORY, added after seeing 10 of 126 ladder rows; the registered test said only "orthogonal or not", the first rows said not, in a direction |
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
| R-12 | R | RANDSPARSE s=10 (3 draws) + SCORESPARSE vs PROD: price, d, cos | **E** | built chain n=126 | pending | pending | EXPLORATORY; the control L9 demands -- separates "sparse is cheap" from "aligned is cheap". Matched to the operator's space: same s, same averaging operator, same projection, same job; only the native's involvement in choosing members differs |
| R-13 | R | d(prod) vs top-75 member spread; d vs cos | **E** | cloud+chain n=126 | pending | pending | EXPLORATORY; is production's off-manifold deviation the AVERAGING artefact (native-free both sides)? |
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
| L-4 | `long40` ladder: pool_best / sparse_s10 / top75_best / avg75 / avg75_random, at L 40–60 and L 9–16 | **registered** (L-H2, P1–P3) | built chain, pre-AMBER | in progress |

