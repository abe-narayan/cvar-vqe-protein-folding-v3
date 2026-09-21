# S32 LEDGER

Newest entries appended at the bottom. Every entry names its artefacts, its basis, and whether it is
ORACLE. Numbers here are asserted by `s32/s32_verify.py` against the artefact that produced them.

**Endpoint: mean built-chain Cα RMSD, `tuning126`, n = 126. Production 3.2105 Å.**
CA point cloud 3.0483 and set mean 3.5507 are *different objects* and are never differenced against
the chain.

---

## S32-L0 -- **THE CHARTER IS STORED VERBATIM AND VERIFIED BY READ-BACK** (2026-09-21 08:00, coordinator)

`s32/BRIEF.md`, 2,087 lines, 51,021 bytes, captured on branch `s26` at commit `373d9176`. Stored
with a delimited provenance header and **verified by reading it back**: all **71** section headers
present with none missing, and the charter's own landmarks confirmed — `1. REPOSITORY`,
`71. BEGIN IMMEDIATELY`, `THE FINAL SCORE IS RMSD`, `LOWER BUILT-CHAIN`, `MDE = 2.8016`,
`Never exceed 8`, the closing `I BELIEVE IN YOU` line, and the trailing instruction about testing on
longer proteins.

**Why this is entry zero.** In S30 the coordinator paraphrased the charter instead of storing it and
the sprint then argued about what had been asked for. S31 fixed it and verified the read-back. This
is the third sprint running the fixed procedure. *A charter that exists only in a summary is not a
specification.*

---

## S32-L1 -- **THE POOL IS NOT THE BOTTLENECK, AND THE PROJECTION COST IS A PROPERTY OF THE OBJECT BEING PROJECTED, NOT A CONSTANT** (2026-09-21 08:00, coordinator)

Extracted from `s29/results/s29_O_ladder_table.json` (S29's O-ladder, 33 rungs carrying both bases).
**This is a cross-sprint quotation and is therefore exactly the class of statement that has been
wrong before** (contract rule 4) — lane V is re-deriving every number from the artefacts that
produced it, and this entry is provisional until it does.

**BUILT CHAIN, all ORACLE rungs ORACLE / NOT DEPLOYABLE:**

```
best sparse convex combination, K=500, s=10   1.1139      <- 2.10 A of headroom
best single member, K=500                     1.7078      <- 1.50 A of headroom
  + distogram SCORE prefix 500 -> 128            +0.4357  ->  2.1435   [see the S32-L3 correction: 0.2477 of this is a pure ORDER STATISTIC and the
                                                        remaining 0.1872 has the SCORE performing WORSE THAN RANDOM] (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)
  + prefix 128 -> 75                         +0.1620  ->  2.3055
  + selection / readout (uniform average)    +0.9051  ->  3.2105   PRODUCTION
```

**Charter §56 asks for the earliest irreversible loss. On this evidence it is not candidate
generation.** The existing K=500 pool supports **1.1139 Å** through a sparse convex combination of
about ten members. Everything downstream of retrieval destroys **2.10 Å that is already present**.
Contract rule 17 is written from this: any claim that generation is the bottleneck must first
explain why 2.10 Å of already-present headroom is not.

**And the second half, which I did not expect and which opens lane R:**

```
PROJECTION COST (built chain minus CA cloud)
  real deposited member (a valid chain)       -0.0007 to -0.0030   ~FREE
  sparse convex combination, K=500, s=10      +0.0002              ~FREE
  sparse convex combination, K=500, s=20      -0.0052              ~FREE
  dense prefix average, K=500                 +0.1701
  PRODUCTION (75-member uniform average)      +0.1622
```

**The projection is free for anything near the manifold of valid chains and costs 0.16 Å for a dense
average, which is not one.** A mean of 75 structures has a contracted backbone (S31: contracts at
*every* separation, least at long range), so the projection must re-expand it, and that repair is
**5% of the endpoint**. Contract rule 16 is written from this: *"the projection costs 0.16 Å" is not
a constant and must not be quoted as one* — every readout proposal states which regime it lands in.

**What this does NOT yet establish.** That the 0.1622 is recoverable. S31 already refuted the
obvious lever: per-target `m` transfers at **+0.0075, 0.22× MDE, 63W/63L — a literal coin flip** on
the endpoint, and the ORACLE global `m` is −0.0044 at 0.19×. So the cost cannot be dodged by
averaging fewer structures. Lane R owns whether it can be dodged at all.

---

## S32-L2 -- **THE READOUT PROGRAM IS CONVEX AND REWARDS SPREAD; AND MY FIRST VERIFICATION OF THAT WAS VACUOUS** (2026-09-21 08:07, coordinator)

Artefacts: `s32/s32_readout_identity.py`, `s32/results/s32_readout_identity.json`. Seed 320001.

**The claim.** For any `Σw = 1`, with `a_x = ‖W_x − t‖²` (**ORACLE**) and `B_xy = ‖W_x − W_y‖²`
(**native-free**):

```
|| sum_x w_x W_x - t ||^2  =  <w, a>  -  0.5 * w' B w
```

`B` is a Euclidean squared-distance matrix, so for `Σv = 0`, `v'Bv = −2‖Σ v_x u_x‖² ≤ 0` with
`u_x = W_x − t`. Hence `w'Bw` is **concave** on the simplex, `−½w'Bw` is **convex** there, and the
program is linear + convex = **convex, with no hyperparameter**.

**Reading the two terms, which is the part worth having in one sentence:** minimising wants **low
`⟨w,a⟩`** (individually good candidates) and **high `w'Bw`** (maximally *spread* candidates — the
variance-cancellation term). ***Spread is rewarded, not penalised.*** That is why quality-blind
dispersion maximisation — this same program with `a` held constant — picks garbage: S31 measured it
at **+0.1436 Å, 1.22× MDE, WORSE**. **`a` is load-bearing and it is the only unknown.**

Verified in four regimes chosen because they are what could break it: a generic pool, a pool **with
duplicate candidates** (real pools have them), a **near-collinear** family, and **every candidate
identical** so that `B = 0`. Identity holds to **6.2e-16** worst-case; conditional negative
semidefiniteness holds in all four.

> **AND THE VERIFICATION I RAN FIRST WAS VACUOUS, one hour after I wrote the rule that forbids it.**
> Contract rule 5 says *a verification must be able to fail*. My first convexity check initialised
> the accumulator at the pass threshold (`worst = 0.0`) and then took a **max over quantities that
> are always negative** — so it reported `+0.000e+00` in every regime, including regimes where the
> true value is `−2.6e-01`, and **it could not have failed on any input.** It is the same shape as
> S31's tie claim "verified" on random floats, which never tie.
>
> **The fix is not a better threshold; it is a positive control.** The shipped version reports
> **both ends** of the range — so a vacuous pass is visible on its face — and runs a symmetric
> non-negative matrix that is **not** a squared-distance matrix, asserting the test **catches** it
> (`+2.685e-01`, caught). *An audit ships with a self-test on the defect that motivated it.*

**Recorded against myself deliberately.** S31's finding was that the coordinator's own claims are
the ones no lane audits. This is the first S32 instance and it was caught by re-reading my own
output, not by a lane.

---

## S32-L3 -- **THE SCORE FILTER IS WORSE THAN RANDOM AT RETAINING THE BEST CANDIDATE, AND MOST OF WHAT I CALLED A "FILTER LOSS" IN S32-L1 WAS A PURE ORDER STATISTIC** (2026-09-21 08:10, coordinator, on lane V's artefact)

> **RETRACTED IN PART at 09:01 by lane V's stratification (S32-L8), and the headline of this entry is wrong as written.** The +0.1872 is **entirely 18 targets**. On the other 108 the score is *marginally better* than random (**−0.0296, 0.30× MDE, NOT A RESULT**, 72W/36L), and it loses on **18 of 18 FAIL18** targets (+1.4879) — a stratum **defined** as targets where no pool member within 1.5 Å of the pool best survives into the top-75, which makes this contrast near-circular on exactly those 18. **Two further errors of mine in this entry:** `n_better = 72` means the score is **BETTER** on 72 of 126 — I wrote “WORSE, 72W/54L”, which reads the W/L backwards — and the **median is −0.0380**, also favouring the score, which I did not report. *Mean +0.187, median −0.038, 57% of targets won by the score: the median-vs-mean gap is the project's free early warning and it fired here.* Original wording left standing per rule 13.*

Lane V's control, `s32/results/s32_V_ladder_orderstat.json`, **2000 draws**, **CA POINT CLOUD**, all
**ORACLE / NOT DEPLOYABLE**. Aggregate re-derived by the coordinator in one script from lane V's raw
rows, per contract rule 7 (a cross-lane claim gets re-derived once before anyone acts on it).

```
best member of K=500                        1.7108
best member of a RANDOM 128 of the 500      1.9586
best member of the SCORE-selected top-128   2.1458
```

**Decomposition of the +0.4350 I called a "filter loss" in S32-L1:**

```
total   best500 -> best128                  +0.4350
  pure ORDER STATISTIC (any 128 of 500)     +0.2477   57% of it -- a random 128 loses this too
  attributable to the SCORE                 +0.1872   the score is WORSE THAN RANDOM (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)
```

`compare(score128, random128)` = **+0.1872, MDE 0.1764, 1.06× — WORSE, 72W/54L**, with `stats_lib` (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)
flagging **`TYPE-M ZONE: magnitude inflated ~1.10×`** — so the honest effect is ~0.17 and it sits
barely past MDE. Top-75 against a random 75 of the 128: **+0.0696 at 1.02×, 65W/61L**, same
direction, barely measured. **The K=500 best member survives into the top-128 on 63/126 targets — a
literal coin flip — at mean rank 170/500.**

> **On half the targets the best available candidate is already gone before the readout ever sees
> it.** That is a far sharper statement of the selection wall than the in-band correlations, and it
> is measured rather than inferred.

**My S32-L1 framing was wrong and is annotated in place.** I called the whole +0.4350 a filter loss.
**57% of it is a pure order statistic** — a minimum over 500 is lower than a minimum over 128 for any
subset — and that part is not attributable to the score at all. *This is contract rule 9 firing
against the coordinator on his own ledger entry, one entry after he wrote the rule.*

**What is NOT established, and it is charter §29's exact trap.** That production would be better with
a random 128. **Production does not take the best member; it averages 75.** By S32-L2 the readout
program rewards low `⟨w,a⟩` **and high spread `w'Bw`**, and a score-selected set is *less diverse than
a random one by construction* — so a score whose *best* member is worse may still select a set that
*averages* better. Consensus has ρ_global **+0.43** and in-band **−0.28**; rejecting outliers for an
average is a different job from finding the best member, and this measurement cannot separate them.

**Registered prediction, before the arm is run** (lane P owns it): replacing the top-128 with a
random 128 and measuring the **built-chain endpoint** leaves it **unchanged or worse**, because the
score's value is outlier rejection for the average rather than best-member selection. **If the
endpoint improves, that is deployable immediately — random selection needs no native information,
which is the wall every other route has hit.** ≥ 8 draws, draw mean and draw-to-draw sd reported,
never the best draw.

---

## S32-L4 -- **THE ENDPOINT'S DEFINITION, RECOVERED: 3.2105 IS NOT A STORED QUANTITY, AND THE CHAIN RUNG IS DISCONTINUOUS IN ITS INPUT AT ONE ULP** (2026-09-21 08:14, lane V)

Lane V's charter Step 4. Artefacts `s32/results/s32_V_step4_endpoint.json`, `s32/s32_step4_rebuild.py`.

**Everything upstream of the chain reproduces bit-exactly.** Targets 126; folds 25/23/25/23/30; the
top-75 **set** reproduced 126/126 from an independent rescoring; the recomputed coordinate average
matches the stored one to **1.42e-14**; CA cloud **3.048338** (want 3.0483); set mean **3.550683**;
pool best K=500 **1.710824**; top-75 best **2.306153**; shipped argmin **3.454000** to six decimals.

**The chain does not, and the reason is structural. Five different objects, not five estimates of
one:**

| value | what it actually is |
|---|---|
| 3.048338 | the CA point cloud — top-75 coordinate average, unprojected (`rmsd_avg`) |
| 3.204076 | the **λ = 0** arm — nearest ideal geometry, no Ramachandran penalty (`rmsd_fit`) |
| 3.214765 | the **λ = 0.3** arm — **the production pipeline's own emitted chain** (`rmsd_arm` / `ca`) |
| 3.235460 | that chain **after AMBER relaxation** (`rmsd_full`) — AMBER costs **+0.0207 Å** |
| **3.210534** | **the canonical 3.2105**: a **re-projection of the stored cloud**, λ = 0.3 arm, the 126 `item="prod"` rows of `s29/results/s29_O_chain_rows*.jsonl` |

> **The canonical endpoint is not a cached scalar.** Its full definition is *the λ=0.3 multi-start
> projection arm of `s12/instrument.project` applied to the production top-75 coordinate average,
> CA-RMSD to native, averaged over tuning126*. On its own terms it reproduces exactly (3.210533995).
> **It is 0.0043 Å better than the chain production itself emits (3.2148).**

**And the mechanism behind the residual, which is the part that binds future measurement:**

```
2LNG   max|C_rebuilt - C_stored| = 7.1e-15   (one float64 ULP; cloud RMSD identical to 9 dp)
         stored cloud  -> chain 4.778535      rebuilt cloud -> chain 4.879181    0.1006 A apart
6QAX   max|dC| = 3.6e-15                      stored 4.226989 / rebuilt 4.078647  0.1483 A apart
```

**The projection is perfectly deterministic given bit-identical input — three in-process repeats
identical to 9 dp, thread count 1/2/4/8 irrelevant — and it is DISCONTINUOUS in that input.**
7e-15 × the ~1e13 amplification = 0.07–0.15 Å, which is exactly what is measured. An independent
re-projection therefore disagrees with the S29 canonical on **126/126** targets: mean **+0.0021**,
p90 **0.026**, max **0.5174** (2LNG — *the same target and magnitude S31's projection pin recorded
as its s27-vs-s29 defect*).

**Consequence, and contract rule 3 is hereby strengthened:** *"same cloud value" is not enough — it
must be the same float64 bits.* Both sides of any chain contrast must be projected in one job from
bit-identical clouds. **Unpaired cross-job chain claims below ~0.03 Å are not resolvable.** Paired
comparisons are untouched: SE of the chain mean is 0.1543, so an unpaired MDE against production is
**0.4324 Å**, and the anchor's 0.002–0.03 wobble cannot threaten a properly paired contrast.

**Ties at the top-75 cut: 22/126 targets (20 with 2, 3 with 3), and they do NOT leak.** Production
uses a seeded random tie key (`rng_for(pdb,"tiekey")`); lane V's independent stable argsort agrees on
the top-75 **set** for 126/126. *But the tie key exists, and an `argmin` on a tied score in new lane
code would leak the pool order* — the defect that once invented a 1.386 Å winner.

**`ok: False` was exactly two assertions, and they are the same one:** built-chain mean 3.212625 vs
3.2105 (+0.0021, tol 5e-4) and the projection price 0.164287 vs 0.1622 (+0.0021, inherited, since
the cloud is exact).

---

## S32-L5 -- **THE READOUT IS THE EUCLIDEAN PROJECTION OF THE NATIVE ONTO THE CANDIDATE HULL, WITH GAIN EXACTLY ONE. SO A QUALITY ESTIMATE GOOD ENOUGH TO MAKE IT WORTH SOLVING IS ALREADY GOOD ENOUGH TO EMIT -- AND S31's SHARPEST OPEN QUESTION IS CLOSED: `a` AND `mu` ARE ONE OBJECT** (2026-09-21 08:22, lane Q; re-derived independently by the coordinator)

Lane Q's derivation. **Re-derived by the coordinator in a separate script before publication**
(contract rule 7): support size 5 of 14, the two objective forms agreeing to **9 decimals**, gain
**0.9991** inside the active affine hull and **5.9e-03** outside, with the check able to distinguish
them. Lane Q's own finite differences: **1.0000000 (sd 3.1e-09)** and **4.3e-08**.

**The identity, read for what it means.** By S32-L2, `‖Σ w_x W_x − t‖² = ⟨w,a⟩ − ½ w'Bw`. The
left-hand side *is* the squared distance from the native to a point of the candidate hull, so

> **the readout's optimisation is the Euclidean projection of the native onto the convex hull of the
> candidates**, `min_{w∈Δ} ‖U'w − t‖²`.

**Consequence 1 — the sufficient statistic.** `∂x*/∂t` is the orthogonal projector onto the
**active** candidates' affine hull: **gain exactly 1 inside it, exactly 0 outside.** So `a` must be
known only along `|S|−1` directions — **about 5 numbers locally, ~33 globally** — and *everything
else in `a` is exactly invisible to the emitted structure*. That is the minimal-information subspace
charter §12 asks for, and it is far smaller than 128.

**Consequence 2 — and this is the closure. [CORRECTED 09:05 by lane V; my original reason was
false.]** I first wrote *"gain exactly 1 means no noise suppression."* **That is backwards.** Gain is
1 on the active affine hull (dimension `|S| − 1` ≈ **5.25**) and **0 on its orthogonal complement —
≈ 33 of ≈ 39 dimensions** (`gain_out_mean` 3.87e-06). **Gain zero is TOTAL suppression, not none**:
the readout annihilates ~87% of a generic error vector by construction, and lane Q's own
`noise_ORACLE` table shows it turning a 3.66 Å estimate into a 2.23 Å emission at ε = 4.0.

**The closure stands, but on the LOWER bound, not on an absence of suppression:**
`d ≤ ‖P_C(t̂) − t‖ ≤ d + ε`, with the **hull floor `d = 1.8290 ± 0.1178`** (CA cloud, ORACLE). The
readout cannot beat the hull floor however good the estimate is, and the measured crossover against
direct emission sits at **ε ≈ 2.2**.

> ### A structure estimate good enough to make the readout worth solving is already good enough to emit.

**(The headline is unaffected: it rests on the measured crossover at ε ≈ 2.2 against the hull floor d = 1.8290, not on the struck claim.) That turns S31's measurement into a theorem.** S31 found that solving the CVaR objective exactly
*"reshuffles the answer on 66 of 126 targets and buys nothing"* (−0.0112 Å at 0.19× MDE) and
recorded it as a surprising null. **It is not surprising and it is not a coincidence — it is forced.**

**Consequence 3 — S31 §20.3 is RESOLVED.** `a` is affine in `t` (3.0e-15 relative), and `Σw = 1`
makes the `‖t‖²` term an additive constant, so the readout-relevant part of `a` is `P_aff{W} t`.
Replacing `t` by `P_aff{W} t` and re-solving moves `w` by **2.5e-14** and the structure by
**6.1e-14 Å**. Since `t = X̄ − μ`, the map **`μ ↔ t ↔ a` is a native-free affine bijection**:

> **A per-candidate quality estimator with in-band skill IS a structure predictor, and a common-mode
> corrector IS a per-candidate quality estimator. One missing channel described in two vocabularies
> — not two requirements.**

S31 left this explicitly open as *"the sharpest question this sprint produces, stated as open rather
than resolved by assertion."* **It is now closed, by derivation.**

**A near-miss lane Q killed with its own control, recorded so nobody re-finds it.** Projecting a
*noisy* structure estimate onto the pool hull looks like a large win — a 3.0 Å estimate emits at
2.15 Å. It is **entirely shrinkage**: against a shrinkage toward the pool mean **norm-matched to the
projection's own displacement**, hull projection is WORSE at every noise level up to 3.0 Å (+0.23 at
ε = 2.0, **6.10×**) and NOT A RESULT at 4.0 (−0.14, 0.41×). *The hull adds nothing over shrinkage.*
(Smoke test n = 3; full-instrument run in flight.)

**What this does NOT close.** Arms that leave the simplex — signed weights, arms that re-embed, or a
cardinality constraint that genuinely binds — are outside the derivation's hypotheses (`Σw = 1`,
unconstrained simplex). **Sparse `s`-of-`K` is the last formulation standing**, and the number that
decides it is whether the cardinality constraint **binds**: if the unconstrained convex optimum's
support is already ≤ s, then the "combinatorial" problem *is* the convex program and the 2.10 Å
sparse headroom is a **convex prize, not a quantum one**. Lane Q is measuring exactly that.

**Q0 separately CONFIRMED S31 rather than falsifying it**, with a self-test that fires on a
deliberately tied vector (0.1221) and a deliberately unsorted one (3.3830) before touching real data
— *unlike S31's, it can fail*. On the real 126: `sc[o]` non-decreasing 126/126; the block model
reproduces `E` **bit-for-bit, max error 0.0e+00, on all 126**; `max|E − ramp| = 0.0406840`. Priced
at **0.26× MDE** (cloud), with the selection readout **identically tied 126/126**, and the
operator-matched control — keep the block *sizes*, relocate the blocks — at **0.38×, the same size**.
**The residual channel is the existence of duplicates, not which candidates duplicate: one integer
per target.** S31 is sharper, not wrong.

---

## S32-L6 -- **THE ENDPOINT IS CLOSED, AND IT IS A POSITIVE RESULT: THE PROJECTION OPERATOR IS BIT-REPRODUCIBLE. WHAT VARIES IS THE INPUT'S LAST BITS, NOT THE STAGE** (2026-09-21 08:40, lane V)

Artefact `s32/results/s32_V_chain_bitexact.json`. Three inputs, **all the same top-75 coordinate
average**, agreeing to 5.7e-14 Å with cloud RMSDs identical to 9 dp:

```
input cloud                            n    mean chain    mean d     max|d|   bit-identical
s29_O_structs/<pdb>.npz['prod']       126   3.210533995  +0.000000  0.000000    126/126
production cache avg_ca               126   3.214765154  +0.004231  0.416765      0/126
recomputed coordinate average         126   3.212625220  +0.002091  0.517410      0/126
```

**This reverses the natural reading of S32-L4.** `s12.instrument.project` is **not flaky**: the
canonical endpoint reproduces from its stated input **bit-for-bit on 126/126**, across processes,
across BLAS thread counts 1/2/4/8, and across the gap between S29 and today, with mean `d` and
max |d| both **exactly 0.0**. That is a stronger reproducibility statement than the project had
before. And re-projecting the production cache's own `avg_ca` reproduces production's stored `ca[]`
**bit-identically on 126/126** — so **production's 3.2148 and S29's 3.2105 are the same operator on
the same mathematical object; S29 recomputed the average and drew different last bits.**

> ### The endpoint, in one sentence, for the report
> **3.2105 Å** is the λ=0.3 multi-start projection arm of `s12/instrument.project` applied to
> `s29/results/s29_O_structs/<pdb>.npz["prod"]`, CA-RMSD to native, meaned over `tuning126`. **The
> operator is bit-reproducible. The *input* is not uniquely determined** — three honest computations
> of the same 75-member average differ in the last bits and give **3.2105, 3.2126 and 3.2148** — so
> the anchor carries roughly **±0.002 Å of pure arithmetic noise**, and the chain production itself
> emits is **3.2148**.

**The charter's targets stay checkable**: < 3.00 is 0.21 Å away and < 2.50 is 0.71 Å away, both two
to three orders of magnitude above the noise. **What is not checkable is any unpaired cross-job chain
claim below ~0.03 Å.** Lane V has queued `s32_V_ulp_distribution.py` (ε = 1e-14 Å per coordinate,
5 draws × 126, asserting the cloud RMSD is unchanged per draw) to convert rule 20 from a caution into
a quotable sd.

### `s32/s32_verify.py` is live: 58 matched, 0 mismatched, 0 flagged, 11/11 self-tests passing

It carries S31's cross-basis audit and adds two that S31 did not have:

- **AUDIT 2, the OBJECT audit.** L4's six-object table is **recomputed live** from
  `bench_results/cache/*` and `s29_O_chain_rows*` at every run, so it cannot rot, and a chain number
  attributed to the wrong object is flagged **even though its basis word ("chain") is correct**.
  ***A basis is not an object*** — that is a genuine advance on S31's audit, which could only catch a
  missing or wrong *basis*.
- **AUDIT 3, cross-job chain deltas.** Any line quoting a chain delta below 0.03 Å that names neither
  same-job pairing nor a bit-identity check is flagged.
- **The path audit now distinguishes MISSING from *pending*** — a path a document itself marks as
  running is not a defect. *Missing, unfinished and crashed look identical to `ls`*, and this is the
  first sprint whose verifier knows the difference.

Self-tests, each fed the real historical defect **and** its corrected form: a chain number declared
cloud (CAUGHT) and the true cloud number (CLEAN); the **emitted** chain 3.2148 called the endpoint
(CAUGHT), the true 3.2105 (CLEAN), the λ=0 arm 3.2041 (CAUGHT), the cloud 3.0483 (CAUGHT); an absent
path (CAUGHT) and a present one (CLEAN); an unsupported cross-job delta (CAUGHT) and the same delta
with a same-job statement (CLEAN); and **ST5, which asserts its own input can exhibit a tie** —
S31's "verified a tie claim on random floats" failure, closed by construction.

---

## S32-L(D1) -- **PHYSICS IS CLOSED AS A RANKER AND OPEN AS A MOVER; AND THE IN-BAND PROBLEM IS THE SIGN, NOT THE SIGNAL** (2026-09-21 09:05, lane D)

Pre-registration `s32/PREREG_S32_D.md`, committed **`34973b1b`** before the first number.
Artefacts: `s32/results/s32_D1_inband.json`, `s32_D1_signrandom.json`, `s32_D1_signshare.json`.
**Basis: in-band Spearman ρ against true Cα-RMSD inside each target's shipped top-75 band, per
target, aggregated over n = 126. This is a diagnostic, not a chain RMSD, and is never differenced
against one.** ORACLE `rr` is an evaluation label; no arm reads it.

### (a) What lane D closed by derivation, before spending anything

Charter §13 permits reopening a closed direction with a named mechanism; **D0 requires the same
burden in reverse — an observable may not be built until it names the hypothesis of G1 it breaks.**
Doing that first closed more than it opened:

- **Theorem D-E.** An achiral potential's stochastic propagator is `O(3)`-equivariant
  (`−∇U` equivariant, thermal noise isotropic), so `⟨A⟩_{x,T,t}` for any achiral invariant `A` is
  **itself an achiral invariant single-structure functional of the seed** — a distance-map reading by
  G1. **Charter §23 (MD, short trajectory ensembles, conformational covariance, basin transitions,
  state populations, metastability, transition rates, autocorrelation, dynamic modes) and §24's
  ensemble reweighting and temperature-dependent response are closed as scalar rankers.** S31 §8
  closed only the single-structure half and recorded that *no thermodynamics was computed*; this
  closes the stochastic half on the same ground.
- **D-E1.** Stochasticity is not an escape: the estimand is deterministic, the thermostat adds
  variance. *Noise is not information.*
- **D-E3.** A deterministic minimiser trajectory is a function of its start. **"Relax it and score
  it" does not escape G1.**
- **Theorem D-F.** `⟨F(x), G(x)⟩` for two `O(3)`-equivariant fields is `O(3)`-**invariant**, hence
  achiral, hence a distance-map reading. **Force-vs-prior-gradient alignment, physics/prior
  consistency and every other contraction of equivariant objects are closed.**
- **Corollary D-G.** The only surviving escape is to **emit a displacement and apply it**.
  *Physics is closed as a ranker and open as a mover* — which is also the shape of the prize, since
  S31 §20.1 prices a **direction**.
- **§0.2c, from the code not from the record.** A pool candidate stores **only** a Cα trace and
  `φ/ψ`; sidechains are **modal rotamers**, hydrogens are **frozen local frames**, and there is no
  PDBFixer in the repository. Every atom AMBER sees is `Ψ(seq, φ, ψ)` for a deterministic `Ψ`.
  **AMBER has no resolution advantage over the torsions, so the "all-atom is a finer map" door is
  closed and no all-atom static score was built.**

**Two claims about chirality corrected against the record before measuring, not after.** The shipped
`ff14SB + gbn2` is **reflection-invariant** (S31 Lemma B1, verified: 1340 torsion phases at distance
0.000e+00 from {0, π}, no CMAP) — **AMBER is not a chiral scorer; this project has none.** And
`core/project.py`'s *"the distance objective is exactly mirror-blind"* is a warning about a
counterfactual objective: its **next sentence** says the shipped one *"is a COORDINATE distance and
is chirality-sensitive"*, with an L-handedness assertion in `test_project.py`.

### (b) AMBER and Legacy measured in band on the shipped instrument — charter §§32, 33

```
                    IN BAND (top-75)          WHOLE POOL (K=500)
scorer            rho      se    xMDE       rho      se    xMDE   folds
AMBER         +0.0000  0.0202    0.00    -0.0266  0.0213    0.45   3/5   62W/64L
DIS           +0.0652  0.0282    0.83    +0.5678  0.0327    6.20   4/5
LEG_total     +0.0376  0.0316    0.42    +0.3073  0.0339    3.23   3/5
LEG_torsion   +0.0444  0.0243    0.65    +0.1676  0.0305    1.96   5/5
RG (control)  +0.0501  0.0403    0.44    +0.2794  0.0382    2.61   3/5
```

**AMBER's mean in-band ranking skill is exactly zero**, and it is the *only* scorer here whose
whole-pool ρ is also wrong-signed. `DIS` at 0.83× and `LEG_torsion` at 0.65× **reproduce S31 §9's
own two figures exactly** — an independent replication of that table on a separately written
instrument.

### (c) The control that changes what the zero means — D1-N

Within-band label permutation (score vector and band size held fixed, 24 draws per target, own
distribution reported):

```
scorer        mean|rho|    null    ratio   xMDE  folds
AMBER            0.1779  0.0948     1.88   2.36   5/5
DIS              0.2496  0.0939     2.66   3.05   5/5
LEG_total        0.2819  0.0908     3.10   3.52   5/5
LEG_torsion      0.2142  0.0947     2.26   2.77   5/5
```

`Var(ρ) > 0` with `E[ρ] = 0`. **In-band ordering content exists on every scorer measured — including
the one with exactly zero mean skill — and what is missing is the per-target SIGN.** This reproduces
memory `in-band-ordering-is-per-target` on a channel it had never been measured on, and it means
*"AMBER cannot rank"* is the wrong sentence: **AMBER ranks, sign-ambiguously.**

### (d) The sign is not shared between independent Hamiltonians, and the reason is compactness

AMBER and `LEG_total` agree on the in-band sign on **36.5%** of targets (z = **−3.0**, *anti*-agreement,
survives Bonferroni at 10 tests). Mechanism: **`LEG_total`'s in-band sign agrees with plain `Rg`'s on
68.3% (z = +4.1) and AMBER's on 31.7% (z = −4.1)** — Legacy tracks compactness in band, AMBER
anti-tracks it, and that is *why* they disagree. Orienting AMBER by Legacy's sign is **−0.0508,
0.92× MDE — NOT MEASURED and pointing the wrong way.**

**ORACLE / NOT DEPLOYABLE price of the sign** (`mean|ρ|`, what the scorer would be worth given a
perfect per-target sign): `RG` +0.3709, `LEG_total` +0.2819, `DIS` +0.2496, **AMBER +0.1779**, all
4.9–5.7× MDE. Memory `in-band-ordering-is-per-target` prices 2.0 Å at ρ = 0.638, so **even a perfect
sign oracle leaves the best of these short of the crossing price by more than 2×.**

### (e) Chirality — D0-X, with a registered prediction that FAILED

The chiral Cα pseudo-torsion scalar has real whole-pool skill (**+0.3302, 2.85× MDE, 5/5 folds**)
that is **lost in band (+0.0607, 0.58× — NOT MEASURED)**. The registered prediction *"in-band share
of chiral variance < 15%"* **failed at 23.6%** and is recorded as failed. What survives is the
gradient: the top-75 filter removes **79% of the pool's chiral variance** (sd 0.3774 → 0.1678),
which is a mechanism for two sprints' *"empty at 9–16 residues"* scope note rather than a restatement
of it.

### (f) The caveat on (b), measured rather than asserted

The cached AMBER energies are **unminimised single points on the ideal-geometry rebuild and are
clash-dominated**: median per-target pool median **29,668 kcal/mol**, p90 **6.83e6**, median
per-target max **1.80e16**, and **58.6% of every 500-pool above 1e4 kcal/mol**. **(b) prices
unminimised AMBER, not AMBER.** Rungs D2-R (relaxed energies in band) and D3-M (physics as mover)
are running and settle it.

## S32-L7 -- **"IN-BAND SKILL IS ZERO" HAS BEEN READ WRONG FOR FOUR SPRINTS. THE ORDERING INFORMATION IS PRESENT ON EVERY SCORER; WHAT IS MISSING IS THE PER-TARGET SIGN** (2026-09-21 08:47, lane D)

Prereg `34973b1b`. Basis: **in-band Spearman ρ against true CA-RMSD inside each target's shipped
top-75 band**, per target, aggregated over n = 126. **Not a chain RMSD** — a diagnostic.

```
                in band       xMDE    W/L
AMBER           +0.0000       0.00x   62/64     <- EXACTLY zero
DIS             +0.0652       0.83x             (reproduces S31 section 9's 0.83x)
LEG_total       +0.0376       0.42x
LEG_torsion     +0.0444       0.65x             (reproduces S31 section 9's 0.65x)
```

**And then the control that changes what the zero means** — a matched **within-band label-permutation
null**, 24 draws per target:

```
scorer        mean|rho|    null     ratio    xMDE    folds
AMBER           0.1779    0.0948    1.88     2.36    5/5
DIS             0.2496    0.0939    2.66     3.05    5/5
LEG_total       0.2819    0.0908    3.10     3.52    5/5
LEG_torsion     0.2142    0.0947    2.26     2.77    5/5
```

> ### `Var(ρ) > 0` with `E[ρ] = 0`. In-band ordering content EXISTS on every scorer — including the one with exactly zero mean skill — and what is missing is the per-target SIGN.

**This is a different problem from the one the project has been trying to solve.** "In-band skill is
zero" has been read as *the information is absent*; it means *the information is present and
unoriented*. It also confirms, by direct measurement on the shipped instrument, an inference
on record from an earlier sprint: *"in-band ordering is learnable but per-target, and the only
leverage supplies the per-target SIGN at inference."*

**D1-S closes the obvious route to the sign.** It is **not shared between independent Hamiltonians**:
AMBER and Legacy agree on **36.5%** of targets — *anti*-agreement at **z = −3.0**. The mechanism is
compactness: Legacy's in-band sign agrees with plain Rg on **68.3%** (z = +4.1) while AMBER's agrees
on **31.7%** (z = −4.1). Orienting AMBER by Legacy's sign is **−0.0508, 0.92× MDE — the wrong way.**

**D0-X, and this is the mechanism that matters for the rest of the sprint.** The chiral CA
pseudo-torsion has **real global skill (+0.3302, 2.85× MDE, 5/5)** that is **lost in band (+0.0607,
0.58×, NOT MEASURED)**. The lane's registered `<15%` in-band-variance-share prediction **FAILED at
23.6%** — recorded as a fired falsifier. But the finding underneath it is:

> **the top-75 filter removes 79% of the pool's chiral variance (sd 0.3774 → 0.1678)** — the
> mechanism behind two sprints' *"chiral channels are empty at 9–16 residues"* scope note. **They are
> not empty; they are filtered out before the readout sees them.**

**Caveat measured rather than asserted:** the cached AMBER energies are **unminimised single points**
on the ideal-geometry rebuild and are clash-dominated — median **29,668 kcal/mol**, p90 **6.83e6**,
**58.6% of every 500-pool above 1e4**. **D-1 prices UNMINIMISED AMBER only.** Rung D2-R settles it.

### The cross-lane proposal this generates — NOT YET ESTABLISHED, and it carries rule 7

Lane D's *"the top-75 filter removes 79% of the chiral variance"* and lane V's *"the score prefix
concentrates on the mode and buys nothing in the good tail — spread halved, 5th percentile
unchanged, and it is worse than random at retaining the best member"* **may be the same event seen
from two sides**: the filter discarding exactly the variance that carries signal. By S32-L2 the
readout **rewards** spread, so the filter would be destroying the chiral signal, the term the readout
wants, and the best member simultaneously.

**The arm:** apply the chiral observable at the **500 → 128** stage, where the chiral variance still
exists, rather than in band where 79% of it is gone. Native-free, deployable, leave-fold-out.

**Neither half has been re-derived by one person in one script, so this is a hypothesis and is
labelled as one.** If lane D's 79% and lane V's spread-halving turn out to be the *same variance
measured twice* rather than two facts, the synthesis is circular. **Every cross-lane synthesis in
S31 failed, four for four**; this one is written down as unaudited on purpose.

---

## S32-L8 -- **THE SCORE IS NOT A GENERAL ANTI-ORDERING. ITS FAILURE IS TOTAL ON ONE TARGET IN SEVEN AND NOT MEASURED ON THE OTHER 108 -- AND THOSE 18 ARE THE STRATUM DEFINED BY THE THING BEING MEASURED** (2026-09-21 09:01, lane V; retracts part of S32-L3)

Artefacts `s32/results/s32_V_orderstat_gate.json`, `s32_V_orderstat_strata.json`, script
`s32/s32_V_orderstat_gate.py` — which regenerates the 2000 draws from the same pinned seed and
**asserts it reproduces `s32_V_ladder_orderstat.json` per-target to 0.00e+00 before computing
anything.** CA POINT CLOUD, ORACLE / NOT DEPLOYABLE.

**Both contrasts formally clear the gate:**

```
SCORE top-128 vs RANDOM 128 of 500
  effect +0.1872  SE 0.0630  MDE 0.1764  1.06x   MEDIAN -0.0380   72W/54L
  folds {-0.1087, +0.2621, +0.2938, +0.2566, +0.2344}   CI95 [+0.0383, +0.2743]   4/5   type_m 1.10
SCORE top-75 vs RANDOM 75 of 128
  effect +0.0696  SE 0.0243  MDE 0.0681  1.02x   MEDIAN -0.0045   65W/61L   4/5
```

**And then the stratification, which is what the entry is about:**

```
                 ALL 126              FAIL18 (n=18)      OTHER 108
500 -> 128   +0.1872 (med -0.0380)   +1.4879  0W/18L   -0.0296 (med -0.0743, 72W/36L)  0.30x  NOT A RESULT (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)
128 ->  75   +0.0696 (med -0.0045)   +0.3611  3W/15L   +0.0210 (med -0.0090, 62W/46L)  0.38x  NOT A RESULT
```

The ten worst targets in the 500→128 contrast — **2BFI, 2NB7, 2JN5, 8T63, 9KAR, 3SGO, 1LB7, 5W52,
3BTB, 2N5C — are all ten in FAIL18**, and the score loses **18 of 18**. Drop the 10 worst and the
aggregate falls from +0.1872 to **+0.0294**; drop 20 and it goes **negative, −0.0567**. (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)

**FAIL18 is defined in `s12/instrument.py::selfcheck` as the targets where no pool member within
1.5 Å of the pool best survives into the production top-75.** A contrast asking *"does the score's
prefix retain the good members?"* is therefore **near-circular on precisely those 18** — and
*directly* circular for the 128→75 arm, which is about the top-75 itself. S31 measured effect size
rising monotonically with a stratum's circularity and named that gradient *the stratum's definition
doing the work*; **this is the same gradient with the stratum left inside the aggregate instead of
named.** Filter-independent control as rule 12 requires: split by chain length (median 13) gives
**+0.1826 short vs +0.1940 long — no gradient.** It is not a length effect. It is the 18.

> ### The score is not a general anti-ordering. On 108 of 126 targets its top-128 retains a marginally *better* best-member than a random 128 — NOT A RESULT. Its failure is catastrophic and total on 14% of targets, and those are exactly the targets already named FAIL18. **The score's problem is not that it orders badly everywhere; it is that on one target in seven it places its window in the wrong part of the pool entirely.**

That localises the loss instead of diffusing it, and it fits lane V's own spread mechanism better
than the sentence it replaces: **the score compresses the same way on every target** (sd 1.32 → 0.64,
5th percentile unchanged 2.6098 → 2.6184) — **but on 18 targets the compressed window is centred in
the wrong place.**

**Contract rule 10, properly applied.** The random arm is a per-target mean over 2000 draws; the
draw-to-draw sd of the **126-target aggregate** is **0.0232** (500→128) and **0.0167** (128→75). A
single random draw — which is what a deployment actually gets — gives +0.1872 ± 0.0232, and **only (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)
69.5% of single draws clear this comparison's own MDE**; for 128→75, **53.6% — a coin flip.**

**The corrected increment sentence, which supersedes S32-L1's:** *narrowing 500 → 128 → 75 costs
0.595 Å of oracle-best headroom; **0.338 is the set-size order statistic**; the remaining **0.257 is
18 targets' worth of the score placing its window wrongly**, and is **NOT MEASURED** on the other
108.*

**Rank moments, both right and different:** mean **170.3**, median **134.0** of 500; the best member
is in the top-128 on **63/126** and in the top-75 on **41/126**. **Name the moment in the same
sentence as the number.**

**Lane V retracted its own earlier sentence to me as part of this** (*"the deployed score's ordering
is worse than chance"*, stated twice without stratification or median) and I propagated it into
S32-L3 with the W/L read backwards on top. **Same failure shape as the whole S31 audit: a quantity
re-used across a boundary its definition does not cross — here a moment and a stratum.**

**Verifier:** now asserts all six stratified numbers, and adds a fourth audit — **any ledger line
quoting an aggregate effect whose FAIL18 and non-FAIL18 strata have opposite signs is flagged
automatically.** That is the generalisable form of this defect. State: 64 matched, 0 mismatched,
0 missing, 0 flagged, 11/11 self-tests.

---

---

## S32-L(D2) -- **THE MISSING PER-TARGET SIGN IS A REAL LATENT AND IT TRANSFERS; AND THE 9–16-RESIDUE SCOPE EXCUSE IS NOT AVAILABLE** (2026-09-21 09:40, lane D)

Pre-registration `s32/PREREG_S32_D.md`, commit **`34973b1b`**. Artefacts
`s32/results/s32_D1_signtransfer.json`, `s32_D2_basin_0_1.jsonl`, `s32_D0X_circularity.json`.
**Basis: in-band Spearman ρ against true Cα-RMSD inside each target's shipped top-75 band. A
diagnostic, never differenced against a chain RMSD.**

### (a) The sign transfers — split-half inside each target, matched null

Sign estimated on a random half A of the band, applied to held-out half B, **16 splits per target**,
own distribution averaged; **matched null = same sign, labels on B permuted**:

```
scorer         transfer   null    excess    se     xMDE  folds
AMBER           +0.1125  -0.0056  +0.1181  0.0164  2.57   5/5
DIS             +0.1890  -0.0033  +0.1923  0.0226  3.04   5/5
LEG_total       +0.2263  -0.0003  +0.2266  0.0237  3.41   5/5
LEG_torsion     +0.1499  -0.0008  +0.1507  0.0196  2.75   5/5
```

**The latent is a property of the TARGET, not of the sample** — which `Var(ρ) > 0` alone could not
establish. AMBER goes **+0.0000 → +0.1125** and Legacy **+0.0376 → +0.2263** once the sign is
supplied. **ORACLE / NOT DEPLOYABLE**: estimating the sign requires native labels on half the band.
Memory `in-band-signal-limited-not-sample-limited` records a **flat learning curve** for supplying
this sign natively, so *the latent existing* and *the latent being natively recoverable* are
different claims and only the first is established here.

### (b) The reality check — both registered horns MISSED, and the hypothesis was mine

10 targets chosen **native-free and deterministically** (first 2 of each frozen fold, pinned order),
24 band members each, free AMBER relaxation, 97.5% converged:

```
CONTRACTION horn   median spread ratio <= 0.70   measured 0.972    DID NOT FIRE
FROZEN horn        median move <= 0.30 A         measured 0.5006   DID NOT FIRE
```

**H-D2 is falsified.** At 9–16 residues the candidates are neither collapsed into one basin nor
frozen. **The scope excuse this lane was entitled to use is not available**, and any future negative
on dynamics must name a different reason.

### (c) Minimisation does not rescue AMBER, and the mechanism is a rotamer artefact

Converged relaxed energy in band **−0.0089** (n = 10, DIAGNOSTIC) against the single point's
**+0.0006**; the relaxed energies are physical (median −887…−293 kcal/mol, per-target sd **1–9
kcal/mol on 8 of 10 targets**), so **D-1's zero is not a clash artefact.** Per-target ρ mean
**−0.0089**, mean |ρ| **0.242** — *zero mean, real dispersion, random sign*, reproduced on converged
all-atom physics.

```
energy, median over candidates      3.18e9  ->  -558 kcal/mol
CA motion required to get there                 0.50 A
band mean CA-RMSD to native         3.1931  ->  3.1849  (-0.0083)   ORACLE / NOT DEPLOYABLE
BEST member of the band             2.5358  ->  2.5741  (+0.0382, WORSE)   ORACLE
rank preservation rho(rr_in,rr_out)             0.9208
```

**Nine orders of magnitude of energy relieved by half an Ångström of Cα motion.** By §0.2c the
sidechains are **modal rotamers** and the hydrogens **frozen local frames** — a deterministic
function of `(seq, φ, ψ)`. **AMBER's dynamic range on this pool is spent on a rotamer-placement
artefact that carries no candidate information**, and the motion it buys leaves the band's quality
unchanged and its best member slightly worse. *That is why the force field cannot rank here.*

### (d) The cross-lane chirality synthesis is CIRCULAR — controlled before it was written down

The coordinator proposed that lane D's *"the top-75 filter removes 79% of the chiral variance"* and
lane V's *"the score prefix halves the pool's quality spread"* are one event seen twice, and asked to
be told if they are one variance measured twice. **They are.** Variance retained under the shipped
filter, with a size-matched **random-75** control and the **exact achiral twin** (`mean cos τ` against
`mean sin τ`, same family, same scale, differing only in parity):

```
                      score-75   random-75
DIS (the criterion)      0.037      0.964
RG                       0.111      0.979
CHI_even_cos  (TWIN)     0.196      0.978
E2E                      0.244      0.957
CHI_odd_sin   (CHIRAL)   0.271      0.975
RR_ORACLE  (lane V's)    0.272      0.965
LEG_total                0.323      0.973
```

The random-75 control retains 0.96–0.98 everywhere, so the shrinkage is **selection**, not
subsetting — but it is **general** selection. Chirality is the **second-least** shrunk coordinate and
is shrunk **less than its own achiral twin** (+0.0755, **2.35× MDE, 5/5 folds**). **Lane D's 0.271
and lane V's 0.272 are the same number to three decimals: one fact, not two.** The proposed
500 → 128 chiral arm is **NOT RUN by derivation** — the achiral twin reaches the same global
magnitude with the opposite sign (**−0.3359** against **+0.3302**), so the G1 escape is not doing the
work, for the **third** time (WRITHE S30-L26, XTWIST S31 §9, this); and what remains is a
weakly-independent global coordinate (+0.1302 past `DIS`, rank-corr **+0.464** with `Rg`) added to a
global filter, which is **S30-L21's closed question** (21 fields span ~2 directions; best
leave-fold-out combination **0.0948**).

**Contract rule 7 in action:** this is a cross-lane claim, it was re-derived from the raw artefacts
in one script by one person before anyone wrote it down, and it did not survive.
