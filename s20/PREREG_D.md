# PRE-REGISTRATION — AGENT D (theory / literature / adversarial review), Sprint 20

Written **before** the corresponding result was inspected, in the order the blocks were run.
Per `s20/BRIEF.md` §7 this file is **not edited after results are seen**; where a
pre-registration turns out mis-specified the defect is recorded beneath it, unedited, and the
untested regime is left **OPEN**.

**TARGET is the unit** throughout (§8). Pair-level numbers are labelled diagnostics and never
reported alone. Every RMSD is full-chain Cα against model 1 of the native, frozen implementation.
**Every arm states its basis** — point cloud / built chain / repaired emission — per §1 as
corrected this sprint.

The sealed 60-target benchmark is not read, probed, derived from, or inferred about, and
`results/benchmark_manifest.json` is not opened. Its SHA-256 is recorded in §A2 as an integrity
certificate; a digest of the bytes carries no information about the file's contents.

---

## Block A — THE STARTING AUDIT (directive §4). Not an experiment; no falsifier applies.

Reproduce, from the code, and report every discrepancy against the reports however small:
the incumbent and the best built structure; the continuous-torsion encoding; the VQE
Hamiltonian, the CVaR implementation, the SPSA implementation and its optimizer state; Legacy
and AMBER; the candidate-generation, coordinate-building and final-RMSD paths; seeds and the
target split; and every experimental branch of Sprints 16–19 with a status label.

**A2 — benchmark integrity without opening the file.** Record `sha256`, size and mtime of
`results/benchmark_manifest.json`, plus the same for the pinned `peptide_clusters.json` and
`peptide_folds.json`. This is the certificate future sprints can check byte-identity against
without ever reading a target name.

---

## Block B — CLOSE SPRINT 19 §5 (CVaR-VQE), reported as NOT DELIVERED at freeze

`s19/results/qb_main.json` is complete on disk and was never analysed. **No new computation.**
Execute `s19/PREREG_B.md` §§6–9 exactly as its author wrote them — P1 realised, P2
eps-dominance with the leave-one-out classical control, P3 generation ceiling — and apply every
one of its nine §8 kill rules mechanically.

**Falsifier of my own handling:** if the (M,D) identity `readout² = M² − D²` does not assert to
< 1e-9 Å², the artefact is not sound and nothing is reported from it.

**I do not get to choose the arms after seeing the table.** "Best quantum" and "best classical"
are the arms with the lowest mean `realised_ORACLE` within their pre-declared sets, full stop.

---

## Block C — ATTACK: is the shared-bias alignment a finding, or the shared target?

**Target of attack.** `s19/CLAIMS.md` S1/S2/S4 and `s20/BRIEF.md` §2: the coherent component of
the distogram's error is *shared* — 0.898 same-seed ceiling, 0.754 cross-architecture, 0.784
retrieval pool — and **"two thirds of it is reproduced by ZERO-INFORMATION references"**
(sequence-blind separation prior 0.575, constant α-helix 0.598).

**Hypothesis (mine, adversarial).** The alignment statistic
`corrC(A,B) = corr(d(X_A) − d_nat, d(X_B) − d_nat)` has a large positive null that was never
subtracted, and the null is an algebraic consequence of the two families being compared against
**the same native**, not evidence of shared bias.

**Expected mechanism.** Write `d(X) − d_nat = (d(X) − d_typ) − (d_nat − d_typ)` for any
reference field `d_typ` that both families share (a separation profile, say). The second term is
**identical for every family by construction**. So `corrC` carries a floor set by the ratio of
the native's own deviation from typicality to the families' deviations. **Any two realisable
length-n peptide structures scored against the same native will correlate**, and if that floor
is near 0.6 then the zero-information references are AT their null and S2 is an identity being
read as a discovery (§9).

**Primary endpoint.** `corrC_null` — the same statistic between **two independent structures
that share nothing but the target**, computed three ways of increasing strictness:

| null | construction | what it controls |
|---|---|---|
| `N1 pool×pool` | two windows drawn independently from this target's own K=500 retrieval pool | realisable length-n peptide geometry, no fit, no predictor |
| `N2 fit×fit` | the **same** `s15/align_lib.fit` from the **same** cached start with the **same** `1/sd²` weights, driven to the distance field of a randomly drawn pool window — two independent draws | matched in the space the operator works in (§6.1): same start, same optimiser, same weights, same manifold |
| `N3 cross-target` | a window drawn from a **different** target of the same length, superposed and scored against this native | removes even this target's pool |

**Falsifier (F-D1), fixed now.** If `corrC_null` for **N2** is **below 0.40** at n=126 with a CI
excluding 0.40, my attack fails: the zero-information references at 0.575/0.598 sit well clear
of the null, S2 survives as stated, and I report it as **survived a serious attack** rather than
burying it. If `corrC_null(N2)` is **at or above 0.575** with a CI containing it, S2's "two
thirds of the ceiling from zero information" is **an artefact of the shared target** and must be
restated as an excess over the null.

**Null.** Stated above; the whole block *is* a null measurement.

**Matched control.** N2 is the matched one: same operator, same start, same weights, same
whitening — this is the §6.1 rule applied to the statistic rather than to an arm. N1 and N3 are
reported beside it because a disagreement between them is informative and, per §6.2, the
measurement is the null.

**Budget.** ≤ 300 L-BFGS fits (126 targets × 2, plus smoke). One process, BLAS threads capped
at 1.

**Promotion criterion.** None — this block can only downgrade an existing claim, never promote
one. If F-D1 fires I report S1/S2/S4 as needing an excess-over-null restatement, with the
excesses computed, and I do **not** claim the mechanism is refuted: a reduced excess is still an
excess.

**Recorded in advance, so it cannot be re-narrated:** I expect the null to be **large**
(0.5–0.7) and I expect S1's *ordering* (ceiling > cross-arch > pool > zero-info) to survive it,
because the ordering is what the mechanism actually needs. What I expect to break is the
**"two thirds"** framing and, with it, `s20/BRIEF.md` §2's sentence that the harmful component
"is largely sequence-independent."

---

## Block E — ATTACK: can `H_Legacy ↔ H_AMBER` be swapped with nothing else changing?

**Target of attack.** `s20/BRIEF.md` §4 Q1: *"Hold constant: sequence, torsion representation,
candidate set, ansatz, qubit count, initialisation, optimiser, evaluation budget, CVaR rule,
measurement budget, convergence threshold. **Change only `H_Legacy` ↔ `H_AMBER`.**"*

**Hypothesis (mine).** That instruction is **not executable**, because the two objects are not
two energy functions on one landscape. `core.quantum.FoldingHamiltonian.energy` evaluates a
knowledge-based energy at the built coordinates. `amber_hamiltonian.AmberHamiltonian.energy`
runs `openmm.LocalEnergyMinimizer.minimize(ctx, tol=2.0, maxIterations=50)` under a
`k = 100 kcal/mol/Å²` positional restraint, then sets `k_rest = 0` and reports the
**unrestrained** energy at the **minimised** coordinates. `H_AMBER = E ∘ Relax₅₀`. Swapping the
Hamiltonians also swaps in a relaxation operator, so every landscape metric compared between
them — roughness, Hessian condition number, gradient noise, basin width — is confounded by it.

**Primary endpoint.** On a small pinned set of targets and a pinned set of bitstrings drawn
identically for both arms: (a) the energy shift `E∘Relax₅₀ − E∘Relax₁`; (b) the **Spearman
rank correlation** between the two orderings — the ordering is what CVaR consumes; (c) the
**fraction of calls at which the 50-iteration cap binds** (Z6/Z6b: a bound must be reported on
the calls it actually bound); (d) `n_collapsed`, the `+inf` guard's firing count.

**Falsifier (F-D2), fixed now.** If Spearman(E∘Relax₅₀, E∘Relax₁) ≥ **0.95** and the cap binds
on **< 10%** of calls, the relaxation is a near-monotone reparameterisation, Q1 is executable
roughly as briefed, and I withdraw the objection to the brief's wording. If Spearman < 0.95 or
the cap binds on ≥ 10% of calls, Q1 as briefed compares `E_legacy` against `E_amber ∘ Relax₅₀`
and any lane running it must either (i) add the missing arm `E_amber ∘ Relax₁` — AMBER's energy
with no relaxation, the only true "change only H" arm — or (ii) label every landscape conclusion
**relaxation-confounded**.

**Null / matched control.** The matched control is the *same* Hamiltonian at a different
iteration cap, not a different Hamiltonian: this isolates the operator from the physics, which is
the whole point. `Relax₁` rather than "no relaxation" because OpenMM's `maxIterations = 0` means
*unbounded*, so 1 is the smallest honest setting and the arms then differ only in the cap.

**Budget.** ≤ 3 targets × 64 bitstrings × 2 caps = 384 OpenMM evaluations, one process.

**Promotion criterion.** None; this is a methodological warning to the other lanes, delivered
before they spend a sprint on a confounded comparison.

**Recorded in advance.** I expect the cap to bind on **most** calls and the Spearman to be high
but below 0.95 — that is, I expect F-D2 to fire on the cap clause and possibly not on the
correlation clause. If it fires on the cap alone, the honest statement is that `H_AMBER` is
"50 iterations of restrained relaxation then evaluate", a *different operator* per Z6a, and it
must be named that way rather than "AMBER".

---

## Block F — the in-manifold κ puzzle (Sprint 19 M7, OPEN)

**Not run as an experiment this sprint unless Blocks C and E leave budget.** Recorded so the
scope is on the record: the open question is why **post-fit** κ predicts the per-target
real−signflip gap (ρ +0.360; +0.482 for the difference) while start-point tangent geometry does
not (+0.134). Sprint 19 volunteered the confound itself — κ and the gap are derived from the
same two fits, so a common cause is not excluded — and **a clean test needs a native-free
predictor of κ**, which does not exist. Until one does, M7 stays **OPEN** and no explanation
offered without one should be believed. I pre-commit to *not* offering a story for it.

---

## Statistics, for every block above

Paired target-level bootstrap (4000 resamples) **and** a fold-clustered bootstrap, medians and
W/L beside every mean, the 5-fold sign vector, stratified by length where n allows. **A
zero-spanning CI without the power to exclude the effect of interest is NOT MEASURED, never
"matched."** A result below the 0.084 Å MDE is not a validated improvement. No min-of-N ceiling
without its min-of-N null, naming its band. Concentration reported against a uniform-effect
null, never as a raw drop-top threshold.

**Labels:** EXACT / ORACLE / ESTABLISHED / SUPPORTED / PLAUSIBLE / OPEN / INCONCLUSIVE /
NOT MEASURED / REFUTED / RETRACTED. For every statistic I ask first whether it is an identity,
a theorem, an implementation consequence, an observation, a causal hypothesis, or a learned
relationship — and I derive the operator before interpreting its statistic.
