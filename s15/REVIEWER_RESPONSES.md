# HOSTILE REVIEWER REPORTS, AND RESPONSES

Three adversarial reviews, written to be as damaging as the evidence honestly allows, followed by
responses. Where a criticism lands, it is **conceded** and the paper is changed. Concessions are
marked **CONCEDED** and the resulting change is stated; there are eleven of them.

The reviews were written before the responses, and no criticism has been softened after the fact.

---

# REVIEWER 1 — a structural bioinformatician

*"I have read this carefully and I do not think it should be published in its current form. Let me
explain why."*

**R1.1 — The method does not work.** The headline architecture — restraint-driven generation — emits
**3.321 Å** where the authors' own baseline emits **3.204 Å**. That is a loss, with a confidence
interval excluding zero (+0.117 [+0.050, +0.188]). Fifteen sprints of effort have produced a method
that is worse than BLOSUM retrieval plus averaging. Why is this a paper?

**R1.2 — The 0.611 Å number is going to be misread and the authors know it.** An ORACLE fit to the
true distance matrix is not a result about structure prediction. It is a statement that a
20-parameter model can fit 55 numbers. Every abstract that leads with an oracle ceiling is doing
publicity, not science.

**R1.3 — The comparison to prior work is not a comparison.** The only external number quoted is
Zhang et al. at 4.89 Å on 75 fragments. Different targets, different lengths, different RMSD
convention. The authors cannot claim to beat it and should not imply they do.

**R1.4 — 126 targets is small, and they were used for everything.** `K` and `m` were fitted on the
same 126 targets that are reported on. The five folds are reused for the distogram, the bias
correction, every hyperparameter, and the reported evaluation.

**R1.5 — 9–16 residues is a regime nobody cares about**, chosen because it is tractable. Peptides
this short are conformationally heterogeneous; a single native structure is arguably the wrong
target, and RMSD to model 1 of an NMR ensemble is close to meaningless.

**R1.6 — There is no experimental validation of anything.**

## Response to Reviewer 1

**R1.1 — CONCEDED, and it is the paper's central claim rather than a hidden weakness.** The
generative architecture **does not beat the incumbent**, and the abstract says so. The paper's
contribution is not a method that wins; it is a **mechanism for why nothing does**, established by
measurements that a working method would not have produced. Specifically: the ensemble *contains*
structures better than the incumbent emits (G = 2.760, −0.444 [−0.573, −0.325]), and the entire loss
is in the readout — selection costs +0.751, aggregation recovers −0.360, projection costs +0.170. A
paper reporting a 0.1 Å win would have told the field far less.

**R1.2 — Partly conceded, and the framing is now constrained by a rule.** The 0.611 Å is labelled
ORACLE in every table row and in prose, and **it is never stated without the 3.644 Å predicted arm
beside it**. But the reviewer's dismissal is wrong on the substance. The number's function is
**attribution**: because the ORACLE and predicted arms share optimiser, starts, selection rule and
manifold, the 3.03 Å between them is *entirely* restraint error. And it overturns a standing belief —
the field's and this project's — that perfect distance knowledge caps this problem at ~1.95 Å. That
figure was measured through a retrieval library; removing the library moves it to 0.611 Å. Every
prior estimate of the value of better distance prediction was computed against the wrong floor and
understated the floor by 3.19× and the addressable headroom by 2.07×. That is a result about where effort should go, and it is not publicity.

**R1.3 — CONCEDED entirely.** The paper states that Zhang et al.'s 4.89 Å is on a different target
set and that **our 3.204 Å is not a like-for-like win**. We claim no comparison to prior predictive
work, and we say why one is not available.

**R1.4 — CONCEDED, and stated in the methods rather than buried.** `K = 500` and `m = 75` were fitted
on the reported targets and were not cross-validated. The sensitivity is published: `K` moves
`pool_best` 1.970 → 1.504 across 100 → 2000 but moves the headline `shipped` only 3.425 → 3.520 —
**flat, so it is not load-bearing for the reported number**; `m` moves `top75_best` 2.609 → 2.106
across 25 → 150 and **is** load-bearing. A separate 60-target benchmark has been held out and **not
read**, and the protocol is frozen before it is touched.

**R1.5 — Partly conceded.** The regime was chosen for tractability and we say so. But the objection
cuts the other way for the paper's actual claims: the mechanistic findings — that scoring functions
calibrated on large proteins invert at short chain length, that predicted-distance errors are
geometrically realizable, that the variational comparison is not robust to how it is reported — are
*about* the short-chain regime and would be untestable in a regime where a single native is unambiguous. On
NMR model 1: the RMSD convention is frozen and published in full, including that full-chain versus
`[1:-1]` differs by up to 2.4 Å on one structure here.

**R1.6 — CONCEDED without qualification.** No structure in this work has been experimentally
validated. Nothing in the paper claims otherwise.

---

# REVIEWER 2 — a quantum algorithms researcher

*"The word 'quantum' is doing a lot of work in this title."*

**R2.1 — There is no quantum computation here.** Everything is statevector or MPS simulation of 12–18
qubits. That is a classical linear-algebra calculation, and calling it a quantum result is a category
error.

**R2.2 — Your own results say the quantum method is worse than doing nothing.** VQE loses to
best-of-N from its untrained initialisation on every structural readout. What exactly is the quantum
contribution?

**R2.3 — CVaR-VQE for peptide folding is already published.** QuPepFold, PLoS One 2026. This is not
new.

**R2.4 — 12-qubit sub-registers of nine targets at k = 4.** Nine targets. Sub-registers. A
configuration space of 262,144 that you could enumerate exhaustively — and do. Nothing here says
anything about a regime where a quantum computer would be needed.

**R2.5 — The negative result is uninteresting because the ansatz is arbitrary.** Try a better ansatz.

## Response to Reviewer 2

**R2.1 — CONCEDED, and the paper claims no quantum hardware result.** All results are exact
simulation, stated as such throughout. We claim **no quantum advantage** and the word does not appear
in that sense.

**R2.2 — The reviewer is describing a claim we no longer make, and the corrected one is better.**
An extended run (19 targets, a budget ladder, and Sprint 14's own selective objective included as an
arm) **reverses the simple negative**: VQE *beats* the untrained control at a budget of 8,192
(−0.202, 13/19, CI excluding zero). But the win is **objective-independent** — its strongest cell is
the selective objective with the *worst* in-tail ordering measured anywhere, so it is not exploiting
the objective; it is **absent at 32,768** (0 of 6 cells); it goes the *other* way on
diversity-respecting set readouts; and it **never beats a classical greedy search** (0 of 16).

So the contribution is not "VQE loses". It is that **the comparison is not robust to budget, readout,
or n — axes the literature routinely leaves unstated** — and that the quantity dominating every arm,
a selection gap of 1.47–2.14 Å, is **identical in the quantum and control arms**, flat across a 16×
budget range, collapsing to 0.29–0.44 Å under oracle restraints. The sampler is not what is being
measured. The control itself remains, to our knowledge, unreported: the nearest published comparison
samples **uniformly at random**, strictly weaker because it does not share the ansatz's support. We also eliminate the standard explanations rather than invoking them: it is **not** barren
plateaus (measured); **not** ill-conditioning (the gradient's null-space share is 0 to 1e-30 and the
bottom eigenvalue decile carries 0.0016–0.0043 of its squared norm where uniform would be 0.10);
**not** a broken CVaR (verified to 4.2e-14 against Rockafellar–Uryasev); and **not** shot noise (the
headline used exact statevector gradients). What remains is a property of the objective, which we
then test by changing the objective class.

**R2.3 — CONCEDED, and the positioning has been changed because of it.** In-loop CVaR for peptide
folding **is taken**. We do not claim it. What remains ours is the **CVaR defect triple** — findings
about the estimator, not applications of it — of which the most consequential is that the gradient
baseline is centred on the tail only, giving cosine **−0.023** with the true gradient. Plus a
measurement we have not seen elsewhere: **α is a diversity dial**, moving the final distribution from
2.7 distinct configurations out of 4,096 draws at α = 1 to 408 at α = 0.01, with entropy 0.73 → 6.58
bits.

**R2.4 — CONCEDED on scope, and it is stated in the methods.** Conclusions are about a
262,144-configuration discrete space at n = 9, k = 4, and are connected to the continuous problem
only by explicit ceiling measurements. We note the exhaustive enumeration is a *feature* of the
experimental design: because the space is enumerable, the certified optimum is known, so the
optimisation axis and the accuracy axis can be reported separately and never substituted. That is
what let us show they come apart.

**R2.5 — Rejected, with a reason.** "Try a better ansatz" is unfalsifiable unless the control is
fixed, and fixing the control is the point. Our claim is not "this ansatz is bad" but "**training
this circuit is worse than not training it**", which is ansatz-relative by construction: the control
*is* the same ansatz, untrained. A better ansatz raises both arms. The claim would be refuted by an
ansatz whose trained version beats its own untrained best-of-N at matched budget — a specific,
runnable experiment we invite, and one that no published work has yet reported.

**One thing the reviewer did not raise, which we volunteer.** Qubit accounting in this literature is
routinely wrong by a constant: four torsions per chain are **inert** for the Cα trace, so the live
register is `(n − 2)·log₂k`, not `n·log₂k`.

---

# REVIEWER 3 — a statistician

*"I am going to be blunt about the inference."*

**R3.1 — You report a great many comparisons and I see no multiplicity control.** Dozens of arms,
several channels, several architectures, and confidence intervals throughout at 95%.

**R3.2 — Several headline effects rest on n = 6, n = 8 or n = 24.** The realizability result is
n = 24; the sign-destruction effect is n = 6; the scale-correction ceiling is n = 8.

**R3.3 — Your own record shows a high retraction rate.** Three retractions in one sprint, two of them
by the coordinator. Why should I believe the rest?

**R3.4 — The fusion "law" is fitted and validated on the same 24 targets.**

**R3.5 — Paired bootstrap CIs on 126 non-independent targets.** Targets share a library, share folds,
and share a distogram.

**R3.6 — You report FAIL18 in every table and then tell me it is meaningless.**

## Response to Reviewer 3

**R3.1 — CONCEDED, and this is the most serious criticism in all three reviews.** There is no
family-wise error control across the programme's arms. Three mitigations are in place and none is
sufficient on its own: every pre-specified statistic is reported as a **single PASS/FAIL verdict**
rather than as fields a reader may select after seeing them; a near-even win/loss record with a CI
excluding zero triggers a **null-calibrated concentration check** before the result is quoted; and
the **60-target benchmark is held out and unread**, so the confirmatory claims are made on data no
arm has touched. **Change made:** the paper now separates *exploratory* from *confirmatory* claims
explicitly, reports the number of arms evaluated per family, and confines multiplicity-sensitive
language to the confirmatory set.

**R3.2 — CONCEDED, and the small-n reads are labelled as such.** The programme's own history makes
this criticism sharp rather than pedantic: an 8-target read put the predicted arm at 2.792 Å where
the full instrument gives **3.644 Å and reverses the sign of the comparison**, and small samples from
the enumerated nine have reversed a conclusion **four times**. **Change made:** every small-n result
carries its n in the sentence that states it, the claim it supports is stated qualitatively
(direction), and full-instrument runs are reported wherever they exist. For the realizability result
specifically, the qualitative claim is supported by three independent measurements — the surrogate
ladder, the null control (real 0.879 Å versus an i.i.d. null of 2.058 Å, ratio 0.427), and the
channel-disagreement geometry — which is why it is stated as a mechanism rather than a coefficient.

**R3.3 — Rejected as an inference, accepted as a warning.** A visible retraction rate is evidence
about the *reporting process*, not about the surviving results, and the alternative — a clean
narrative — is what should worry a reviewer. Two of the three retractions were caught by checking
machinery against its own stated assumptions (a uniform-grid lookup on a non-uniform grid; a
per-process salted hash), not by scientific intuition, and in both cases the **wrong result was more
plausible than the right one**. The paper reports the retractions in full, including a 25.8% figure
this coordinator reused in four places without checking, which measurement puts at 3.5%.

**R3.4 — CONCEDED.** ~~The `d_avg ≈ √(r² − (s/2)²)` law is validated on all 126 targets, having first been examined
on 24, with one channel pair, using single fits rather than ensembles. **Change made:** it is reported
as a **geometric identity with an empirical calibration check**, not as a fitted model — it has no
free parameters, which is why a 0.162 Å prediction error on unseen quantities, across all 126 targets, is meaningful — and the
paper states explicitly that it must be checked on a second channel pair before being quoted as
general.~~

> **⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT workstream. The response above conceded too little,
> and its remaining defence is now void.** The reviewer's concern was under-stated in three ways.
> (i) The expression is the **two-member Krogh–Vedelsby (1995) ambiguity decomposition** — prior art
> and not ours. (ii) `r` was computed as the **arithmetic** mean where the identity requires the
> **quadratic** mean; corrected, the residual is **+0.075 Å**, not 0.162 (paired −0.087
> [−0.119, −0.061] i.i.d.-target, [−0.129, −0.056] fold-clustered, W/L 126/0, n = 126). (iii) "a
> 0.162 Å prediction error on unseen quantities" is the sentence that must go: **the identity
> predicts nothing unseen.** With both structures in one common frame and the correct mean it is
> exact to **2.65e−15** on all 126 targets, so what was reported as a prediction error was a
> mixture of the wrong mean (54%) and a frame-convention mismatch (46%). It does not need checking
> on a second channel pair; it is an algebraic rearrangement and holds for any pair.
> `s16/retract_law.py`, `s16/retract_exact.py`, `s16/retract_FINDINGS.md`.

**R3.5 — CONCEDED, and partly addressed.** Targets are not independent: they share a fragment
library, five folds, and one distogram. Intervals are computed **fold-aware** where folds apply, and
per-fold effects are reported so a reader can see heterogeneity. This does not fully solve the
dependence and the paper says so. It also discloses the concrete leakage that exists: 16 of 126
universes contain a verbatim own-fold window and 4 targets' top BLOSUM hit is the target's own
sequence — with a **measured impact on `pool_best` of exactly zero** on all four, since they are
different structural determinations 0.59–4.13 Å from native. Those targets are **retained**, because
removing them is the selective removal this programme forbids.

**R3.6 — CONCEDED, and the reviewer has identified a real inconsistency.** `I.FAIL18` is a threshold
artefact: `BAND = 1.5 Å` has no derivation, the zero-recall count runs **45 → 2** across BAND 0.5 →
3.0, and on a 99-combination sweep **exactly one of the 18 targets survives every threshold** (mean
Jaccard 0.498). It is reported only for cross-sprint comparability. **Change made:** the caveat is
attached at first use rather than in a footnote, no argument rests on set membership, and the set is
never described as "the hard targets".

---

# WHAT THE REVIEWS CHANGED

1. Exploratory and confirmatory claims are separated, with arm counts per family (R3.1).
2. Every small-n result carries its n inline and is stated as a direction (R3.2).
3. ~~The fusion law is presented as a parameter-free geometric identity with a calibration check, not
   a fitted model (R3.4).~~ **SUPERSEDED, Sprint 16, 2026-09-06:** it is cited as the Krogh–Vedelsby
   ambiguity decomposition, the mean is corrected, and the "calibration check" is withdrawn — the
   identity is exact and the measured residual was a frame convention (R3.4).
4. The FAIL18 caveat moves to first use (R3.6).
5. No comparison to Zhang et al. is claimed (R1.3).
6. The ORACLE ceiling never appears without the predicted arm beside it (R1.2).
7. CVaR positioning moves off "CVaR for peptide folding" and onto the defect triple and the diversity
   measurement (R2.3).
8. Quantum scope is stated as exact simulation of a 262,144-configuration space (R2.1, R2.4).
9. Constant-fitting on reported targets is stated in the methods with published sensitivities (R1.4).
10. The absence of experimental validation is stated plainly (R1.6).
11. Fold dependence is acknowledged as only partly addressed (R3.5).
12. **The quantum claim is reframed from "VQE loses" to "the comparison is not robust to budget,
    readout or n", after an extended run reversed the simple negative** (R2.2). This concession was
    forced by our own data after the reviews were written, and it strengthens the paper: a
    dependence is a more useful result than a sign.

**A partial answer to R3.1 that arrived after these reviews were written.** An **empirical
false-positive floor** has since been measured directly, by running a comparison that is zero by
construction: it returned **+0.081 [+0.014, +0.169]** on one start draw and −0.003 on another. So a
paired interval can exclude zero on a null effect at ~0.08 Å in this machinery. That does not solve
multiplicity, but it gives a **calibrated magnitude below which no claim is made**, applied
retrospectively to four of this programme's own results including its one positive physics finding.
It is a floor obtained by measurement rather than by correction, and it is reported as such.

**The criticism with no adequate answer** is R3.1. The programme has evaluated a large number of arms
without family-wise error control, and the only real defence — a held-out confirmatory instrument —
has not yet been spent. That is the reason the 60-target benchmark stays sealed until the protocol is
frozen, and it is the honest limit of everything claimed before then.
