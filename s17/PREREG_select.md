# PRE-REGISTRATION — SELECT workstream, Sprint 17

Written before any number in `s17/sel_*.py` was produced. Every experiment below states
hypothesis · expected outcome · strongest control · success criterion · **falsifier**.
The findings file records, per experiment, **which rule fired**.

Problem labels are BRIEF §3: **A** retrieval · **B** selection · **C** ensemble · **D** validity.
This workstream is almost entirely **B**.

---

## Standing methodology for every experiment here

1. **TARGET is the unit.** n = 126. Never a row, never a (target, K) cell.
2. **Two intervals, both reported**: a target bootstrap (comparable with the programme's
   history) and a **fold-clustered** bootstrap over the 5 folds (conservative; 5 clusters is
   coarse and that is stated, not hidden).
3. **Both mandatory controls on every arm**: matched random selection of the same count from
   the *same* candidate set, and a zero-information reference.
4. **Tie-safe selection.** `argmin` on a tied score reads the candidate ORDER, which is
   BLOSUM-informative. Every selection statistic averages the outcome over the **whole tied
   argmin set** (S8's `sel_of` trap, which cost a whole table). Enforced in `sel_lib.sel_of`.
5. **Variant selection is leave-fold-out.** A best-of-N over 126 targets is an in-sample
   number and is quoted only as an overfitting upper bound, never as a result. The honest
   number is: choose the variant on 4 folds, evaluate on the 5th, pool the 5 held-out sets.
   This is the exact guard the `score-axis-does-not-transfer` post-mortem says was missing.
6. **ORACLE / REALIZED / GAP are printed separately in every table.**
7. No native quantity enters any score, weight, threshold, or architecture choice.
   `u["rr"]` and `nat_ca` are evaluation only.

---

## E1 — The identical-candidate-set instrument  (Problem B; BRIEF §3B, sprint §10)

**Hypothesis.** When every selector consumes the *same* structures, the distance objective is
the only arm that beats matched random with an interval excluding zero; Legacy is
indistinguishable from random as a ranker; Legacy-as-a-gate is neutral-to-positive.

**Expected outcome.** `sel-dist` < `sel-random` with the CI excluding zero at every K.
`sel-legacy` ≈ `sel-random`. All arms far above the ORACLE ceiling.

**Strongest control.** Matched random of the same count, drawn with `stable_rng`, averaged
over enough draws that its own SE is small relative to the effect being priced.

**Success criterion.** The instrument reproduces the pinned constants (pool best 1.7108,
top-75 best 2.3062) and separates the arms with the two controls attached.

**Falsifier.** If `sel-dist − sel-random` has a CI containing zero at n = 126, *"the distance
objective transfers across targets"* is REFUTED and the sprint's strongest selection signal is
not a selection signal at all.

---

## E2 — Deep decomposition of the distance objective  (Problem B; sprint §12)

**Hypothesis.** The shipped Bayes-risk mean over pairs is not the best functional of the same
predicted distribution. Specifically the shipped score (a) weights all pairs equally although
their predictions differ by an order of magnitude in confidence, (b) averages a per-pair loss
whose scale is pair-dependent, so a few high-dynamic-range pairs dominate, and (c) is a
*point-wise* loss that ignores what the candidate set itself says.

**What is different from `score-axis-does-not-transfer` (11 refuted variants).** Those variants
changed how the prior is consumed *in generation* and were screened best-of-N on 10 targets
that sit inside the benchmark. This experiment (i) holds the candidate set **fixed** so it is a
pure Problem-B measurement, (ii) runs at n = 126, and (iii) selects the variant **leave-fold-out**.
The prior from that memory is therefore *pessimistic*, and it is recorded here as the expected
outcome, not as a reason not to run.

**Variant families.** pair weighting (sequence separation, uncertainty, entropy) · residual
form (Bayes risk, NLL, absolute, z-scored, log) · robust aggregation (trimmed, median, CVaR
over pairs, Huber, Tukey) · per-pair standardisation across the candidate set · contact-vs-
distance mixtures · spectral / low-rank distance-matrix comparison · and the
**centre-blend family** below.

**The centre-blend family, stated separately because it is the one novel idea.**
Typicality — the mean deviation from the candidate set's *own* mean distance profile — is the
only signal in the programme's record with positive in-band skill. It is algebraically the
shipped residual objective with the predicted centre replaced by the set's centre. So define

        score_beta = mean_p  w_p * rho( (D_p - m_p(beta)) / s_p ),
        m_p(beta)  = (1-beta) * expected_p + beta * setmean_p

with `beta = 0` the distance objective and `beta = 1` typicality. **Both endpoints are known
quantities and the interior has never been measured.** This is a candidate-set-conditional
objective, i.e. sprint §13's framing in its cheapest form.

**Expected outcome.** Most variants are null. Uncertainty weighting and per-pair
standardisation are the most likely to move. Prior probability that any variant survives
leave-fold-out selection: honestly low (~30%), because of the standing refutation.

**Strongest control.** (i) the shipped objective on the identical candidate set; (ii) matched
random; (iii) the in-sample best-of-N number printed beside the LFO number so the overfitting
gap is visible.

**Success criterion.** The LFO-selected variant beats the shipped objective by an interval
excluding zero at n = 126, target as the unit, *and* the fold-clustered interval also excludes
zero, *and* the win/loss is not carried by fewer than 10 targets (drop-top-10 keeps the sign).

**Falsifier.** If the LFO-selected variant's advantage over shipped has a CI containing zero,
**the space of functionals of the existing distogram is closed for selection as well as for
generation**, and the `score-axis-does-not-transfer` result is extended from generation to
selection at n = 126. That is a publishable null and it is reported as the headline if it fires.

---

## E3 — Selection as probabilistic inference  (Problem B; sprint §13)

**Hypothesis.** `P(near-native | candidate set)` estimated by a *listwise, set-conditional*
model beats `argmin score`, because the useful information is the candidate's position in the
set's own score/geometry distribution rather than its absolute score.

**Design.** Features per candidate are **set-relative and native-free**: z-score of the score
within the set, rank fraction, distance to the set medoid, typicality, agreement with the set's
top-scoring quantile, per-pair agreement statistics. Labels are native RMSD **inside training
folds only**. Leave-fold-out at the fold level; inference native-free.

**Expected outcome.** Positive but small; the record says a linear model already saturates
within-target capacity and that the loss is entirely in transfer.

**Strongest control.** (i) the same model with **only** the shipped score as a feature
(so any gain is attributable to set-conditioning, not to learning); (ii) matched random;
(iii) a zero-information constant predictor.

**Success criterion.** Beats the shipped objective on held-out folds with a CI excluding zero
AND beats the score-only model.

**Falsifier.** If the set-conditional model does not beat the score-only model, **set
conditioning adds nothing** and sprint §13's framing is closed at this feature class.

---

## E4 — Target-level calibration  (Problem B; sprint §14)

**Hypothesis.** The missing quantity is a **per-target scalar supplied at inference**. Two
independent lines say so (in-band 0.986 within / 0.600 across; the compactness residual's
unpredictable per-target sign).

**Design.** A per-target native-free feature vector — length, sequence composition,
hydrophobicity/charge patterning, predicted secondary-structure content, distogram-predicted
Rg and its uncertainty, retrieval dispersion and clustering, score variance / skew / entropy,
best-vs-median score gap, candidate agreement. Used to select or blend *among* objective
variants per target, or to set a per-target beta in E2's centre-blend family.

**Baselines that the model must beat, all three:**
1. **uncalibrated** — one global variant for every target;
2. **constant** — the single best constant calibration (a zero-information reference);
3. **linear** — a linear model on the same features, so any nonlinear model must earn it.

**Success criterion.** Beats all three on held-out folds with a CI excluding zero.

**Falsifier.** If a per-target calibrator does not beat the *constant* baseline — the exact
failure mode Sprint 15's four estimators hit (0.357–0.413 against a 0.516 constant) — then
per-target calibration is REFUTED at this feature class and sprint §14 is closed for it.

**Boundary trap, declared in advance.** If a leave-fold-out selector lands on the boundary of
its ladder it is choosing an endpoint, not a parameter, and is reported as a substitution.

---

## E5 — The hard targets against the FULL-universe ceiling  (Problems A and B; sprint §38)

**Hypothesis.** A substantial part of the "53 targets whose pool best exceeds 2.0 Å" is a
truncation artefact of `I.pool_idx`'s K = 500 cut, and the residual hard set after the full
universe is consulted is a smaller, structurally distinct population.

**Expected outcome.** The hard set shrinks; what remains is enriched in long chains and in
targets whose predicted secondary structure is neither pure helix nor pure strand.

**Control.** Length-matched and fold-matched comparison against the easy set; every
"structural class" statement carries the class's own base rate.

**Success criterion.** Descriptive — a characterisation with base rates, not a p-value.

**Falsifier.** If the hard set is unchanged by the full universe, **recall is genuinely the
blocker on those targets** and they belong to the retrieval workstream, not to selection.

---

## E6 — Does the extra ceiling reach the selector?  (Problems A vs B; the sprint's pivot)

**Hypothesis.** Widening K buys ceiling and not realized accuracy — ranking is the blocker.

**Success criterion / falsifier are symmetric and both informative.** If the recovered
fraction is near zero at n = 126, **recall is not the blocker and this workstream is the whole
sprint**. If it is large, the pivot is wrong and retrieval width is the cheapest angstrom on
the table. Either answer redirects the sprint, so this is run first among the analyses that
depend on the coordinator's map.
