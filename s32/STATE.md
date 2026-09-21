# S32 — RUNNING STATE

Maintained per charter §52/§55. Newest notes first, below the headline. **With six lanes an
unwritten research state is how contradictory conclusions coexist unnoticed.**

Endpoint: **mean built-chain Cα RMSD, `tuning126`, n = 126 — production 3.2105 Å.**
CA point cloud 3.0483 and *set mean* 3.5507 are different objects.
Primary target < 3.00 Å. Ambitious < 2.50 Å.

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

The final rung runs at **λ = 0**, where the Ramachandran penalty has zero weight, so on degenerate
targets a ~1e-7 coordinate-distance gap — i.e. rounding — picks the branch.

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
