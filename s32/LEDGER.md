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
                                                        remaining 0.1872 has the SCORE performing WORSE THAN RANDOM]
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
  attributable to the SCORE                 +0.1872   the score is WORSE THAN RANDOM
```

`compare(score128, random128)` = **+0.1872, MDE 0.1764, 1.06× — WORSE, 72W/54L**, with `stats_lib`
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
