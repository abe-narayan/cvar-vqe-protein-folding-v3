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
  ordering effect (DIS 128 - random 128)  +0.1872     WRONG SIGN
128 -> 75 TOTAL                           +0.1604
  size effect (random 75 of 128)          +0.0907
  ordering effect (DIS 75 - random 75)    +0.0696     WRONG SIGN
```

The deployed score's ordering is **worse than chance** at retaining oracle headroom.
Mechanism (V-A4): the score improves the pool's mean rr 4.4533 -> 3.5847 and halves the spread
(sd 1.3206 -> 0.6384) while leaving the 5th percentile **unchanged** (2.6098 -> 2.6184). It
concentrates on the mode and buys nothing in the good tail, so a random 128 samples the pool's
fat left tail and the score's 128 does not.

Correct sentence: *"narrowing 500 -> 128 -> 75 costs 0.595 A of oracle-best headroom, of which
0.338 is the set-size order statistic and 0.257 is the score ordering performing worse than
random truncation."*

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
