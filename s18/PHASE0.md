# SPRINT 18 — PHASE 0: THE 19-TARGET RESULT REPRODUCES, AND WAS NEVER SIGNIFICANT

Independently recomputed from `s17/results/quantum_theory.json`, 19 exhaustively enumerated
targets, certified argmins (no search error).

## The numbers reproduce exactly

    full objective certified argmin        mean 2.661 A
    degree-1 (Walsh weight <= 1) argmin    mean 2.411 A
    space best                             mean 1.030 A
    Walsh: weight-1 variance fraction      0.613
    Walsh: cumulative <= weight 2          0.930

Every figure Sprint 17 quoted is confirmed. **Phase 0's gate is passed on reproduction.**

## But the inference drawn from them was never statistically supported

    deg1 - full   mean -0.249   median +0.000   W/L 7/5   95% CI [-0.650, +0.093]

- **The interval includes zero.**
- **The median is exactly 0.000** -- on 7 of 19 targets the two objectives have the *identical*
  argmin, and on those the difference is not small, it is nil.
- The mean is carried by three targets: 7VI4 (3.936 -> 1.396), 1CS9 (5.256 -> 3.655), 7T3H
  (3.156 -> 1.791). Against it: 8HVS (1.824 -> 2.678), 1N9U (2.782 -> 3.495), 5V5B (3.810 -> 4.398).
- **W/L 7/5 with 7 ties** is a coin flip.

## And a tension that points the other way

    rho(objective, RMSD)   full +0.264    degree-1 +0.153

**The degree-1 truncation is a WORSE global correlate of RMSD than the full objective**, even
though its argmin is better on average. Whatever is happening is a property of the *argmin*, not
of the objective's alignment -- which is exactly the distinction Phase 5 of the brief demands be
kept separate.

## Consequence for the sprint

The degree-1 branch is **not** entering Sprint 18 as an established result that must be
translated. It enters as a **point estimate at the noise floor of its own instrument**, with a
plausible mechanism and a median of zero.

**Falsifier F4 ("the effect exists only on the 19-target exhaustive instrument") is therefore
live before the 126-target test begins**, and the honest prior is much weaker than Sprint 17's
report implied. The 126-target experiment is still worth running -- the point estimate is
favourable, the mechanism is derivable, and no cheaper instrument can answer it -- but it is now
a test of a hypothesis, not a confirmation of a finding.

**Sprint 17's report is amended**: it quoted "2.411 vs 2.661" as a lead without its interval.
That is the sprint's own recurring failure mode -- a number quoted without the functional the
claim needs -- and it is recorded here rather than in a footnote.
