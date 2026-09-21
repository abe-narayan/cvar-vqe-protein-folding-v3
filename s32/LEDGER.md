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
  + retrieval filter K=500 -> 128            +0.4357  ->  2.1435
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
