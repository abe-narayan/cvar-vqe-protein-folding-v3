# SPRINT 19 — PRE-REGISTRATION, AGENT A (distogram / structure prediction)

Written **before** any arm was run. Not to be edited after seeing results; a mis-specification is
recorded as such and the untested regime is left OPEN.

Instrument: 126 cluster-disjoint tuning targets, `s12/instrument.py`, pinned folds.
Primary metric: **mean full-chain Cα-RMSD** of the structure emitted by the deployed refinement
(`s15/align_lib.fit` on `((d − dhat)/sd)²`, start = ideal-geometry projection of the coordinate
average of the shipped top-75), exactly `s18/objceil.py`'s `a0.0` arm.
Reference constants that must reproduce: `avg` 3.048, `a0.0` 3.610, `shuffled` 2.609,
`shuf_strat` 2.560, `isotropic` 2.573, `a1.0` 1.152.
MDE 0.084 Å. Seeds via `s15/seed.stable_rng`. BLAS threads = 1.

---

## Q1 — PAIRWISE ERROR TOPOLOGY (diagnostic, ORACLE label axis)

**Hypothesis A1.** The predictor's ORACLE error |r| = |dhat − dtrue| is concentrated in the
**low-utility** region identified by the Sprint-18 compass (long-separation, low-confidence pairs),
and is therefore *anti-correlated with utility*.

**Pre-specified predictions.**
- ρ(|r|, sep) > +0.3 and ρ(|r|, 1/sd²) < −0.3 (the latter is a reproduction of D9's −0.476).
- Regressing |r| on native-free features, sep and sd carry most of the explained variance.

**What would refute A1.** ρ(|r|, sep) ≤ 0.1, or the error being uniformly distributed over the
utility axes.

**Note on interpretation, fixed in advance.** A1 being TRUE does **not** by itself explain the
Sprint-18 result (permuted residuals of the same size are better). A1 concerns *where the error is*;
Sprint 18 concerns *what shape it has*. Both are reported, and I will not let a confirmed A1 be
read as a mechanism for S18.

## Q1b — IS THE HARMFUL ERROR IN A NATIVE-FREE-DETECTABLE SUBSET?

**Hypothesis A2.** There exists a native-free pair subset S, |S| ≤ 25% of pairs, such that setting
those pairs to the truth (ORACLE repair, matched pair-count and matched residual-RMS removed)
recovers materially more than repairing a random subset of the same size.

**Primary statistic.** ΔRMSD(repair S) − ΔRMSD(repair random-matched-size), paired over 126 targets,
bootstrap CI.

**Falsifier.** If no native-free criterion beats the size-matched random repair by more than the MDE
with a CI excluding zero, then the harmful error is **not** native-free localisable and every
"fix the bad pairs" intervention is closed.

## Q2 — IS THE ERROR COHERENT?

**Hypothesis A3 (the central claim of my lane).** The residual field r is *structurally coherent*:
`dhat` is much closer to being a realisable Euclidean distance matrix of **some (wrong) structure**
than a magnitude-matched incoherent field is. The fit therefore converges confidently onto a wrong
structure, whereas an incoherent field of the same size has no consistent structure to converge to
and the least-squares compromise stays nearer the truth.

**Primary statistic (ORACLE diagnostic).** The *coherent fraction*
`κ = 1 − ||dhat − d(X*)||_w / ||dhat − dtrue||_w`, where X* is the weighted best-fitting structure.
Equivalently, the terminal objective value f* of the deployed fit, normalised by the residual's
weighted norm.

**Pre-specified prediction.** κ(real) − κ(permuted) > 0.15, and per-target κ correlates positively
with per-target RMSD failure (ρ > +0.25).

**The decisive falsifiable experiment (A3b) — the whitening test.** Construct a residual field that
preserves each pair's residual **magnitude and its pair assignment** (so weights stay aligned,
avoiding the G5/G6a trap) but destroys the *cross-pair correlation structure* by randomising signs
only. Compare with the real field.
- If sign-randomisation alone recovers most of the 1.00 Å `shuffled` gap → the harm is **coherence**,
  and A3 is SUPPORTED.
- If sign-randomisation recovers little → the harm is *which pair carries which magnitude*, not the
  correlation structure, and **A3 is REFUTED**. I will report that as the headline.

I explicitly record in advance that sign-flipping is only a partial whitening (it preserves |r|'s own
spatial pattern); a second arm replaces r by a magnitude-matched draw from an independent field.

## Q3 — METRIC REALISABILITY / MULTIMODALITY

**Hypothesis A4.** The predicted distance field violates Euclidean realisability measurably
(triangle-inequality violations, negative Gram eigenvalues, rank-3 defect), and the *degree* of
violation predicts per-target failure.

**Pre-specified statistic.** Per target: fraction of triples violating the triangle inequality;
`edm_defect` = Σ|λ_k| for k>3 over Σ|λ| of the double-centred −½ J D² J; ρ(edm_defect, a0.0 RMSD).

**Falsifier.** ρ ≤ 0.1 → the realisability defect does not price failure and Q3 is closed as a
diagnostic. NOTE, fixed in advance: A3 and A4 pull in **opposite** directions (A3 says the field is
*too* realisable — coherently wrong; A4 says it is *not* realisable — incoherently wrong). They
cannot both be the mechanism. I pre-commit to reporting whichever is falsified.

## Q4 — THE INTERVENTION: a utility-shaped training loss

**Observation motivating it (checked in source, not in prose).** `core/predict.py::train_fold` calls
`MLP.fit` with `sample_w=None`, so the deployed loss is **uniform over pairs**, and the module's own
comment argues capacity should be moved toward **long-range** pairs ("an unweighted loss spends
capacity where the model is already right"). The Sprint-18 compass says the opposite: accuracy on
long-range/low-confidence pairs converts **worse than uniform** (+0.165, +0.244), while accuracy on
short-range/confident pairs converts best (−0.294, −0.308).

**Hypothesis A5.** Retraining the predictor with a loss weighted **toward short separations**
(the utility direction) lowers Cα-RMSD through the deployed fit, even though it *raises* aggregate
MAE.

**Primary outcome, pre-registered.** Mean Cα-RMSD of `a0.0` (deployed fit, leave-fold-out predictor,
n=126) for `sw=util` versus the deployed `sw=uniform` predictor. Paired bootstrap CI, fold-aware,
median and W/L reported beside the mean.

**Success.** Mean improvement > 0.084 Å with a CI excluding zero.
**Falsifier.** CI spans zero, or the sign is positive → A5 REFUTED and the "reshape the training
loss" branch is closed for this direction.

**Controls, mandatory.**
1. `sw=long` — the *opposite* weighting (`separation_weights('lin')`, the anti-utility direction).
   A5 predicts this is worse than uniform. If util and long are indistinguishable, the training
   weight does nothing and the arm is uninterpretable.
2. Seed control: the uniform arm is **retrained from scratch under the identical script** so the
   comparison is not against a differently-produced checkpoint.
3. MAE is reported for every arm and is **not** an outcome. A lower MAE without lower RMSD is
   recorded as a negative result.

**Explicit anti-goals.** No arm is selected on the benchmark-60. No hyperparameter is chosen on the
outcome; the weight shapes are fixed a priori as `min(sep/3, 6)` (long, the existing implementation)
and its reciprocal-shaped mirror (util). If I sweep shapes, every rung is reported and the argmin is
labelled "argmin-on-the-mean over the tuning instrument", not an improvement.

---

## Reporting rules I bind myself to

- Every number carries n, median, W/L, and a paired CI. No mean alone.
- A CI spanning zero at insufficient power is **NOT MEASURED**, never "matched".
- Every arm that uses `dtrue` is labelled **ORACLE** in the same sentence as its number.
- If Q4 is negative I lead the report with it.
