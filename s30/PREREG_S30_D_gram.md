# PREREG S30-D-gram — the field class's Gram matrix, effective rank, and combination ceiling

**Written before the measurement completed.** The field vectors were being regenerated when this
was written; the analysis (`s30/s30_D_gram.py analyse`) had not been run on more than 2 targets.
Posted so the predictions below can be wrong in public.

## The question (coordinator, 2026-09-20, following S30-L5)

S30-L5 established that the 21 native-free displacement fields are individually real (11 of 21
with fold CIs excluding zero against the correct signed null; best +7.7 σ) and individually worth
0.0195 Å, because the transfer function `RMSD_prod·√(1−ρ²)` is quadratic near zero. **What ρ does
the best combination of the 21 reach?** If they were orthogonal, √(Σρᵢ²) with eleven fields at
ρ ≈ 0.10–0.11 is ≈ 0.33–0.36, and 3.00 Å needs ρ = 0.358.

## The dictionary, stated so the units are fixed

Lane T's bits↔ρ map is `bits = −((3n−6)/2)·log₂(1−ρ²)`, the Gaussian direction-location capacity
in the rigid-body-removed space. With the dev set's mean ambient dimension **3n−6 = 32.9**:
ρ = 0.14 → 0.470 bits, ρ = 0.358 → 3.25 bits. Both reproduce lane T's quoted values to 3 s.f.,
which is how I know the dictionary I am using is the one the record uses.

## THE PREDICTION, and it is the reason the control is not optional

The fields live in a space of mean dimension 32.9. **Projecting any direction onto a random
21-dimensional subspace of a 32.9-dimensional space captures E[ρ²] = 21/32.9 = 0.638 by dimension
counting alone — ρ ≈ 0.80, which in the dictionary is 24 bits.** That is obviously absurd as an
information statement, and it is what an unconstrained ORACLE combination ceiling will return.

So I predict, before looking:

1. **The unconstrained ORACLE ρ\* will be ≈ 0.78–0.85** and will be **statistically
   indistinguishable from a matched random 21-dimensional subspace.** If it is, the number is an
   artefact of peptide length, not a property of the fields, and must not be quoted as a ceiling.
2. **The Gram will be strongly collinear** — the fields are all differences of averages over
   overlapping subsets of one pool, ranked by scores that all contain DIS. I predict a stable rank
   (trace/λmax) **below 4** and a top-3 eigenvalue share **above 70%**.
3. **The LFO global weighting will land near the best single field** (ρ ≈ 0.11), not near 0.358,
   because `decorrelated-errors-exist-but-are-unusable` prices fusion gain as the square of the
   weaker channel's skill, and because a collinear set has one direction to give.

## The falsifier

The interesting outcome, and the one that would reopen the bound, is:

> **the rank curve's EXCESS over a rank-matched random subspace is materially positive at some
> rank r, AND the leave-fold-out global weighting reaches ρ ≥ 0.20 with a fold CI excluding the
> best single field's ρ.**

ρ = 0.20 implies 2.988 Å through the bound and 0.97 bits — it would not by itself reach the
charter's target, but it would triple the best measured native-free cosine and make the class
worth pursuing. **If the LFO ρ does not clear the best single field by ≥ 1.0× MDE, the combination
question is closed and I will say so.**

## What is ORACLE and what is not

- the Gram, the effective rank, the rank curve, ρ\*, the ORACLE global weighting: **all ORACLE**
  (`u` is the direction to the native). They are ceilings, not methods.
- **only the leave-fold-out global weighting is native-free**, and it is the only number that may
  be read as a capability.
- the ORACLE weighting is fitted with the native and will flatter itself; it is an order statistic
  over weightings and is priced against matched random subspaces for exactly that reason.

## Multiplicity

The rank curve emits 21 comparisons and the rest emit 4; 25 for this file. Under the sprint's
running count a single 1× MDE positive in a 21-point curve is expected and will not be reported
as a result.

Artefacts: `s30/s30_D_gram.py`, `s30/results/s30_D_gram.json`,
`s30/results/s30_D_gram/<pdb>.npz`. Regeneration of the fields is asserted equal to
`s29/results/s29_D_fields_rows.jsonl` cosines to 1e-9, so this measures the same fields S30-L5
priced.
