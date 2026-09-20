# PREREG S30-F3 — what is the distogram confidently wrong about on the pools where its filter inverts?

Lane F, Sprint 30. Committed before any number for this experiment exists.
Code `s30/s30_F_score.py`; output `s30/results/s30_F_score.json`.

## The question, and why it is not a ninth router

S30-L2 located the tail's damage in the distogram Bayes-risk filter; S30-L3 dissolved the stratum
table into one quantity (the filter's ORACLE set-mean benefit, ρ = +0.81 with the effect of
widening) and closed the global-width family by ceiling. What is left is the *cause*: why does a
score that adds 1.09 Å of set-mean quality on 108 targets add nothing on 18?

This experiment measures the score's **in-pool ranking skill per target** and asks what predicts
it. It builds **no router and no detector.** Every predictor it tests is measured against the
native (ORACLE) precisely so that the question is "what is the score wrong about", not "can we
detect it". If a native-free instantiation is ever proposed, its feature correlation with the
closed families below must be reported first.

**The closed router feature families, enumerated from the record so that this lane cannot
accidentally rebuild one** (S22-L7, S23-L7, S26-L110/L115, S28-L6; note the record's own tallies
are inconsistent — S26-L115 calls its own the "sixth through seventeenth" — so "eight routers" is
cited as a tally, not an enumeration):

```
F1  target length (n, length quartile)
F2  the objective's own score distribution / retrieval confidence (score shape, sim entropy, gaps)
F3  posterior entropy (distogram bin entropy, "dg entropy")
F4  candidate-set geometry (spread(m), pax_*, top-75 spread, pool consensus spread)
F5  compactness disagreement (rg_z, rg_gap: predicted Rg vs realised Rg) -- CLOSED at 2% of its MDE
F6  ESM contact-map statistics (con_*)
F7  pool statistical-potential distribution and its rank agreement with the distogram
```

## The measurement

Per target, over that target's **full 500-member BLOSUM pool**:

```
rho_pool   Spearman(shipped Bayes-risk score, ORACLE CA-RMSD to native)   -- positive = skill
```

The score order is recomputed from the posterior with the S30-L3 reproduction gate (recomputed
top-75 must equal the production record's `sub` set-wise on 126/126). The S29-L19 caveat that
deep order below the top-75 cut may differ by ~2-in-126 from a fresh posterior is carried.

**This quantity does not exist in the record.** The pool-ranking measurements (S29-L33, S29-L50,
S14 via S29-L19) are unstratified; the FAIL18-stratified measurements (S28-L23b's gradient cosine,
S29-L2's ladder rho, S29-L20's axis cosine, S29-L35's field cosines) are structure-level cosines
or rank correlations over hand-built ladders, not within-pool ranking skill.

## Falsifiers, registered before the numbers

**F3a — the score's in-pool skill collapses on the tail.** Fires only if `rho_pool` on the other
108 is positive with a fold CI excluding zero, **and** the FAIL18 mean is at most **half** it,
**and** the FAIL18 mean sits outside the central 95% of a 20,000-draw random-18 null (seed 30002).
The two filter-independent tails of S30-L2 (worst-18 by pool mean; worst-18 by ORACLE best-in-pool)
are reported in the same table, because S30-L3 measured that a FAIL18-only effect which does not
replicate on them is a statement about production's definition and not about hardness.

**F3b — the collapse is the distogram's own error.** ORACLE predictor: the distogram's error
against the **native's** true pair distances. Fires only if Spearman(`rho_pool`, distogram error)
over the 126 is **≤ −0.40** with a fold-clustered CI excluding zero.

**F3c — lane L's scale hypothesis, registered as theirs and tested as stated.** Lane L derives that
every fixed-reference pair channel contains a separable term that is a pure function of scale. If
that is what drives the inversion, then the distogram's **scale error** (predicted mean pair
distance minus the native's mean pair distance) should predict `rho_pool` **over and above** its
shape error. Fires only if the partial Spearman of `rho_pool` with scale error, controlling for
shape error, is ≤ −0.30 with a fold CI excluding zero. **Registered prediction: I do not know, and
I am recording that I have no prior here**, because lane L's derivation was relayed to me second
hand and I could not reach lane L directly to get the exact form of the term. If F3c fires on a
statistic I guessed rather than the one lane L would have specified, the result is provisional
until lane L confirms the statistic.

**F3d — confident wrongness, not wrongness.** The project's record says "confidently wrong costs
2-3× absent" (`torsion-restraints-reach-the-target`). Fires only if the interaction
(error × confidence) predicts `rho_pool` better than error alone, by ≥ 0.10 of Spearman, where
confidence is the posterior's own concentration. **Note this uses F3-family (posterior entropy)
material and therefore cannot become a detector**; it is diagnostic only, and will be labelled so.

## Discipline

- Every predictor here is ORACLE (it reads the native's distances). No deployable object is
  produced, and no parameter is tuned on the native (contract rule 9).
- Correlation of every statistic used with the closed families F1–F7 is computed and reported in
  the same table, whether or not it is flattering (contract rule 6's spirit; the charter requires
  it explicitly for anything detector-shaped).
- `n` (length) is residualised out of every headline correlation and both versions are reported,
  because length was the dominant confound in S28-L6 (its two best singles were length proxies at
  corr −0.43 and +0.75).
- MDE per comparison, fold-clustered CIs, folds-same-sign (rules 2–4). Below 0.7× MDE is NOT
  MEASURED. Multiplicity: this experiment runs 4 registered falsifiers plus the F1–F7 correlation
  table, and the count is carried into the sprint-wide tally (rule 26).
