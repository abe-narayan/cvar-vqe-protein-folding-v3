# PRE-REGISTRATION — the steering experiment (`s16/steer.py`)

Written **while the n = 126 run is in flight** and before any full-instrument number has
returned, so that the reading of the result cannot be chosen after seeing it. The n = 8 smoke
below is the only evidence in hand at the time of writing, and it is explicitly *not* a result:
this project has had a conclusion reversed by a small-n read four times.

## The question

Sprint 15 left two separately estimable **native-free** quantities:

1. a direction `u = θ_pool − θ_fit` that has `|cos(u, e)| ≈ 0.36` with the true torsion error `e`;
2. a Jacobian `J = ∂(superposed CA)/∂θ` whose spectrum splits torsion space into a small **loud**
   subspace that carries almost all of the RMSD and a large **quiet** subspace that carries almost
   none.

The hypothesis is that projecting the weak native-free direction onto the loud subspace removes
the part of `u` that cannot help, raising the usable fraction and turning a 0.36 cosine into an
RMSD improvement.

## What the smoke (n = 8) actually said

```
native-free |cos(u, e)| = 0.318      top-4 singular share of J = 0.883
rank-4 loud carries 0.804 of the RMSD with 0.218 of the error norm   (ORACLE decomposition)
first-order prediction ||J e||/sqrt(n) = 6.183 A against an actual 2.725 A   (ratio 2.269)

full          +0.382 [-0.043,+0.898]
loud4         +0.272 [-0.008,+0.659]
loud10        +0.495 [-0.010,+1.240]
quiet_half    -0.106 [-0.185,-0.031]      <-- the control moved
ORACLE_quiet_half  -0.008 [-0.101,+0.086]
ORACLE_loud4  +0.320 [+0.000,+0.961]
ORACLE_true   -2.665 [-3.093,-2.236]
```

Three things in that block are hostile to the hypothesis and are recorded here **before** the
full run so they are not quietly dropped:

- **the projection made it worse, not better.** `loud4` (+0.272) is no better than `full`
  (+0.382), and `loud10` is worse than both. Filtering the native-free direction through the
  loud subspace did not raise its usable fraction at n = 8.
- **even the ORACLE loud step fails.** `ORACLE_loud4` is **+0.320** — stepping the *true* error's
  loud component makes RMSD worse — while `ORACLE_true` reaches 0.060 Å. A decomposition that says
  loud carries 0.80 of the RMSD, whose loud step then loses ground, is a **first-order statement
  being read as a finite-step promise**. The curvature ratio 2.269 says the linearisation is not
  valid at this error magnitude, which is the mechanism.
- **the control moved.** `quiet_half` is −0.106 with a CI excluding zero. The framing requires the
  quiet arm to do *nothing*. Its ORACLE counterpart does do nothing (−0.008), so this is not the
  quiet subspace carrying hidden RMSD — it is the native-free step doing something that is not
  steering.

## Pre-registered decision rule

Read at **n = 126**, target as the unit, paired bootstrap, step size chosen leave-fold-out.

| outcome | verdict |
|---|---|
| `loud4` or `loud10` beats `full` **and** its CI excludes zero on the improving side | the projection is doing the claimed work — **hypothesis supported**, proceed to integrate steering into the flagship architecture |
| `loud*` does not beat `full`, and both are ≥ 0 | **hypothesis refuted at the mechanism level**: the loud/quiet split is real but does not make the weak native-free direction usable |
| `quiet_half` improves with a CI excluding zero | **the framing is wrong, whatever the loud arms do.** A control that must be inert is not inert; no loud-arm number may be quoted until the quiet arm is explained |
| `ORACLE_loud4` is worse than plain fit while the rank-4 RMSD share is > 0.8 | **the first-order decomposition may not be quoted as a finite-step budget** anywhere in the sprint, including in the Sprint 15 provenance sections |

**No step-size, rank, or arm may be added after the n = 126 numbers are seen** for the purpose of
finding an arm that works. A new arm is permitted only as a *diagnostic of a failure*, and must be
labelled as post hoc in the final report.

## If the quiet control moves at n = 126

The follow-up is fixed in advance (`s16/curv.py`): decide between

- **(a) curvature** — the quiet subspace is defined at `θ_fit`; a finite step rotates `J`, so a
  direction quiet at the start is not quiet along the path. Test: shrink the step until the
  measured RMSD change matches `||J·δ||/√n`; if the gain vanishes as the step shrinks, it is
  curvature and the loud/quiet language survives *locally* only.
- **(b) a non-steering side effect** — the step toward `θ_pool` is a pull toward retrieval torsions
  that improves RMSD for reasons unrelated to the error direction. Test: replace `u` with a random
  direction of the same norm projected into the same quiet subspace. If a random quiet step gains
  the same amount, the arm measures the step, not the direction.

Both tests are native-free except where an ORACLE arm is labelled.
