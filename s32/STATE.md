# S32 — RUNNING STATE

Maintained per charter §52/§55. Newest notes first, below the headline. **With six lanes an
unwritten research state is how contradictory conclusions coexist unnoticed.**

Endpoint: **mean built-chain Cα RMSD, `tuning126`, n = 126 — production 3.2105 Å.**
CA point cloud 3.0483 and *set mean* 3.5507 are different objects.
Primary target < 3.00 Å. Ambitious < 2.50 Å.

---

## NOTE 2 (2026-09-21 09:37, coordinator): **LANE Q CLOSED. THE QUANTUM QUESTION IS ANSWERED — NO ON THIS INSTRUMENT — AND THE REASON IS CHAIN LENGTH**

`S32-L(Q1)–L(Q3)`. Three things settled:

1. **Sparse `s`-of-`K` falls by monotonicity**, and it is the cleanest closure of the sprint. It
   **escapes all three of S31's obstructions** — *those were properties of the candidate-index
   register, not of CVaR-VQE* — and then dies to a different argument: `f*(s)` is constant for
   `s ≥ s*` and strictly worse for `s < s*`, so ***the hard instances are exactly the ones whose
   optimum is worse.*** `s*` mean 10.06, 61.1% ≤ 10.
2. **Charter §14: NO, on this instrument.** `n_res` is 9–16, so `2^n_res ≤ 65536` on 126/126 and
   the enumerability condition fails **on chain length alone**. The two properties a problem would
   need: **P1 a decision space that grows with the target** (exists only where one fragment no
   longer spans it) and **P2 a genuinely stochastic energy** — *they must hold together and this
   project has never had either.* **That is the user's “test on longer proteins” arriving as a
   DERIVED requirement rather than a suggestion**, and it is what lane L is now pricing.
3. **Six reals buy 83% of what thirty-three buy** (3.0532 → 2.0312 → 1.8290, CA cloud, ORACLE).
   *A bit count prices a selection alphabet; this decision is continuous.*

**The correction that changed the sprint's central theorem, and it was against me.** I wrote that
*“gain exactly 1 means no noise suppression”* and made it the reason S32-L5's headline holds. **Lane
V refuted it from lane Q's own table**: gain is 1 on the active hull (≈ 5.25 dims) and **0 on its
≈ 33-dimensional complement** — *gain zero is TOTAL suppression*, ~87% of a generic error
annihilated, and at ε = 4.0 a 3.65 Å estimate emits at **2.23 Å**. **The headline survives on the
LOWER bound instead: the hull floor `d = 1.8290`, crossover at ε ≈ 2.2.** Charter §4 names *“a
theorem was stronger than the empirical conclusion”* as an S31 failure mode; this was one.

**Live now:** lane R (branch pricing at 32/126, plus the alignment result), lane P (random-128
endpoint arm, 84 distinct of 126), lane D (**is the per-target sign predictable native-free? one
bit, 126 targets** — the sprint's clearest remaining route), lane L (the longer-protein instrument,
now with a derived reason to exist), lane V (ULP distribution, then the R and D adversary passes).

---

## NOTE 1 (2026-09-21 08:25, coordinator): **THE SPRINT'S SHAPE AFTER THREE HOURS — ONE THEOREM CLOSED THE READOUT, AND MY OWN LADDER FRAMING WAS WRONG TWICE**

**What is now established (L2–L5), and it is more than I expected this early:**

1. **The readout is the Euclidean projection of the native onto the candidate hull, with gain exactly
   1** (lane Q, re-derived by me). **`a` is needed only along `|S|−1` ≈ 5.25 directions** — the
   *active* support — out of an ambient `d_mean` ≈ 38.9; everything else is *exactly invisible*.
   Gain is 1 on that support's affine hull and **0 on its ≈ 33-dimensional complement — total
   suppression, not none** [corrected: I first wrote `|S|−1 ≈ 33`, conflating the *support*
   dimension with the *annihilated complement* dimension. They are the two halves of the same
   split, 5.25 + 33 ≈ 38.9, and only the first is what `a` must be known on]. The closure is the **hull floor `d = 1.8290`**, which the
   readout cannot beat: ***a structure estimate good enough to make the readout worth solving is
   already good enough to emit.*** **Arrow 4 is closed by
   derivation**, and S31's "solving it exactly buys nothing" is now forced rather than surprising.
2. **`a` and `μ` are ONE object** — `μ ↔ t ↔ a` is a native-free affine bijection. **S31's sharpest
   open question is resolved.** In-band skill and common-mode correction are one missing channel in
   two vocabularies.
3. **The endpoint is not a cached scalar** (lane V). 3.210534 is a *re-projection* of the stored
   cloud at λ=0.3, and it is **0.0043 Å better than what production itself emits (3.2148)**. Five
   distinct objects live near it.
4. **The chain rung is discontinuous at one ULP**: 7e-15 Å in → **0.10–0.15 Å out**, deterministic
   given identical bits. Contract rules 20–22 written from this.
5. **The score prefix is worse than random** at retaining the best candidate (+0.1872, 1.06× MDE, (median −0.0380; **entirely FAIL18** — −0.0296 at 0.30× on the other 108, NOT A RESULT)
   Type-M ~1.10×), and **57% of what I called a filter loss is a bare order statistic**. Mechanism:
   *the score concentrates on the mode and buys nothing in the good tail* — 5th percentile unchanged
   (2.6098 → 2.6184) while the spread halves.

**Two corrections to me, both annotated in place.** I called the 500→128 step "the retrieval filter"
— it is the **distogram score prefix**, *"the quantum field of view"*, 128 = 2⁷. And I told lane R
the final rung is λ=0 with the Ramachandran prior inactive — **it is λ=0.3 and the prior is already
selecting among branches**, exactly as `lam_path`'s docstring intends. *A penalty is not a selector,
though, and post-hoc branch ranking remains open.*

**And one against myself with no lane involved:** my first verification of L2's convexity was
**vacuous** — accumulator initialised at the pass threshold, max over always-negative quantities,
could not fail on any input. Contract rule 5, violated one hour after I wrote it. The shipped check
carries a positive control.

**Where that leaves the sprint.** Arrows 2–4 are one information problem and it is now bounded by a
theorem. **The live routes are: (i) lane P's random-128 endpoint arm, which needs no `a` at all and
is therefore the cheapest possible deployable gain; (ii) lane R's branch degeneracy, which is the
only loss whose missing quantity may be computable today; (iii) whether the sparse cardinality
constraint binds — if it does not, the 2.10 Å sparse headroom is a convex prize, not a quantum one.**

---

## HEADLINE HYPOTHESIS (provisional, hours in): THE LOSS IS DOWNSTREAM OF RETRIEVAL, AND PART OF IT IS ARITHMETIC

The S29 O-ladder, read on the **built chain**, says the pool already contains the answer:

```
best sparse convex combination, K=500, s=10   1.1139      2.10 A of headroom
best single member, K=500                     1.7078      1.50 A of headroom
PRODUCTION                                    3.2105
```

**So candidate generation is not the earliest irreversible loss.** (Contract rule 17.)

And the projection cost is **a property of the object projected**: ~0 for a real member or a sparse
combination, **+0.1622 for production's 75-member dense average**, which is not a valid chain.

### The lead that opens lane R, from `core/project.py`'s own docstring

> **"A CA trace admits two ideal-geometry torsion solutions at near-equal objective distance, ONE
> RAMACHANDRAN-PLAUSIBLE AND ONE NOT."** … **"THE REFERENCE DISAGREES WITH ITSELF, by up to 1.6 Å**
> … the structures this stage returns on those targets are **not determined by the objective; they
> are determined by the arithmetic."**

~~The final rung runs at **λ = 0**, where the Ramachandran penalty has zero weight, so on degenerate
targets a ~1e-7 coordinate-distance gap — i.e. rounding — picks the branch.~~
**STRUCK — this was wrong.** The λ-ladder is `(0.0 → 0.3)` and **the canonical arm is λ = 0.3**
(`s12/instrument.py:139`), so **the Ramachandran prior is already active and already selecting among
branches** — `lam_path`'s docstring states that as its purpose. See contract rule 22. *A penalty is
not a selector, though: post-hoc ranking of converged branches remains untested.* Original wording
left standing per rule 13; it was sent to lane R before it was caught and is corrected there too.

**Why this might be the first native-free signal with a reason to work.** Every ranking signal this
project has tested is a distance-map function and therefore **achiral** (theorem G1). A torsion
branch choice is exactly what an achiral observable cannot see and **a Ramachandran score can** —
the Ramachandran surface is strongly chiral, native-free, target-specific through `res_classes(seq)`,
and already implemented in that module.

**What is NOT established:** that the ORACLE-best branch is far from production, that any of it is
capturable in band, or that the best-of-N inflation leaves anything after it is priced. Lane R owns
all three. S31's branch-carry at λ=0.3 was **−0.0055 at 0.44× MDE, a NULL** — a different operation,
but it means the obvious version is already tried.

---

## LANES

| lane | question | status |
|---|---|---|
| **R** reconstruction | Is the +0.1622 projection cost recoverable? Is the branch degeneracy an *opportunity*? Does a chiral criterion have in-band skill over branches? | running |
| **Q** CVaR-VQE | Falsify S31's target-invariance first. Then: how precisely must `a` be known for `w*` to be useful? That sensitivity question closes or opens the whole family. | running |
| **P** pool/selection | Decompose pool quality / diversity / common-mode / ranking / readout / reconstruction independently. Is positive in-band skill obtainable at all, or is there a theorem forbidding it? | running |
| **D** physical response | Does a **chiral** native-free observable have in-band skill where achiral ones provably cannot? Cheap first; name the G1 hypothesis violated. | running |
| **V** verification/adversary | Independently rebuild 3.2105 from artefacts. Own the verifier, multiplicity, and the attack on every promising result. | running |

---
