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
