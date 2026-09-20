# PREREG S30-F2 — filter width vs averaging width, and whether a wider filter is route (b)

Lane F, Sprint 30. Committed before any number for this experiment exists.
Code `s30/s30_F_width.py`; outputs `s30/results/s30_F_width_cloud.json`,
`s30/results/s30_F_width_chain.json`.

## Why this experiment and not filter-on/filter-off

S30-L2 established that the distogram Bayes-risk filter (500 -> 75) is the stage that fails on
hard pools: against a random 75 of the same pool it is +1.7674 Å worse on FAIL18 (0W/18L, random-18
null p = 0) and NOT MEASURED on the other 108, while on the set mean — the quantity the terminal
operator consumes — it buys −1.0855 Å on the 108 (5.11× MDE, 5/5, 100W/8L) and that benefit
vanishes on the tail. Turning the filter off would surrender the 108's benefit to rescue 18
targets and should lose. **The real variable is filter WIDTH**, and a wider filter needs no
detector, which puts it in route (b).

## The confound this experiment exists to remove

In production the top-75 is **both** the filter and the averaging set, so sweeping one number moves
two mechanisms. S29 swept m as the **averaging width** and found it flat (−0.00047 Å per unit m,
ORACLE global value −0.0018, leave-fold-out +0.0079). That flatness is a property of the averaging
width and says nothing about the filter width. Here they are separated:

```
filter width  k : keep the top k of 500 by the shipped Bayes-risk score
averaging width m : uniformly coordinate-average a RANDOM m-subset of those k
```

`k = m = 75` reproduces production. `k = 500, m = 75` is the record's "no filter" arm. Holding
`m = 75` and sweeping `k` moves the filter alone — the measurement nobody here has made.

Grids fixed now: `k in {75, 100, 128, 150, 200, 250, 300, 400, 500}`,
`m in {25, 50, 75, 100, 150}` with `m <= k`, 8 random subsets per cell, seed 30003, common random
numbers across k so the curves are paired.

## Phases

1. **Set-member statistics, ORACLE, point cloud, free.** `mean(rr[top-k])` and `min(rr[top-k])`
   per target per k. Answers the saturation question directly, since `operator-consumes-set-mean`
   gives d_out = 1.16·d_set_mean + 0.04·d_set_best (R² 0.89).
2. **The actual operator, point cloud.** RMSD of the coordinate average of the sampled m-subset.
   Coordinate averaging contracts the backbone ~25.8%, so the set mean is a prediction and this
   is the test of it.
3. **The endpoint, built chain**, on the cells phase 2 selects — the only basis that decides
   anything (contract rule 1).

## Falsifiers, registered before the numbers

**F2a — the coordinator's hypothesis, that width is a free lunch.** Fires only if some `k > 75`
satisfies **both**: (i) on the other 108, the set-mean benefit versus `k = 500` retains **>= 80%**
of what `k = 75` retains; and (ii) on FAIL18 the excess best-member harm versus `k = 500` falls by
**>= 50%** against its value at `k = 75`. Measured at phase 1, confirmed at phase 3.

**F2b — the deployable claim.** `k` is a single global scalar, so it may not be chosen on the 126
natives (contract rule 9). The deployable arm picks `k` **leave-fold-out** — on four folds'
natives, applied to the fifth — and is compared to production on the **built chain**. Fires only at
**>= 0.7× MDE** better with a fold-clustered CI excluding zero. The ORACLE global best `k` is
computed and reported **separately and labelled ORACLE**, never as the deployable number.

**Registered expectation, so that a null cannot be spun afterwards: I expect F2b to FAIL.** S29
priced four single global scalars; two had an ORACLE global optimum of *exactly* zero and all four
were on the wrong side of zero leave-fold-out (`s29/s29_O_FINDINGS.md` U1), and the project's memory
says "tune one number better" is not a strategy this instrument rewards. The reason to run it
anyway is that all four of those scalars were readout or displacement parameters and **none was a
filter width**; F2a can be true while F2b is false, and a true F2a with a false F2b is still a
finding about where the harm lives, reported as such and not as a win.

**F2c — the separation itself.** If the `k` curve at fixed `m = 75` is flat to within 0.7× MDE at
every `k`, then filter width behaves like averaging width, the coordinator's prior is wrong, and
the no-gate branch of route (b) closes. That outcome is registered as a result, not a failure.

## Controls and labelling

- The random m-subset is the control in the operator's own space; its expectation over 8 seeds is
  the priced quantity, and the seed sd is reported beside it.
- A deterministic deployable instantiation (every `k/m`-th member in score rank) is run beside the
  randomized arm so the family has an emittable member.
- Every ORACLE row is labelled ORACLE in every sentence. `rr` and averaged-cloud figures are
  POINT CLOUD and are never mixed with built-chain rows.
- MDE per comparison (2.8016 × SE), fold-clustered CIs, folds-same-sign reported (rules 2–4).
- A FAIL18 claim carries a random-18 null (20,000 draws, seed 30002), and the two
  filter-independent tails from S30-L2 (worst-18 by pool mean; worst-18 by ORACLE best-in-pool)
  are carried alongside, because S30-L2 measured that FAIL18's magnitude is inflated ~2.7× by its
  own definition.
