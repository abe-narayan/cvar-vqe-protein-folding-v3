# SPRINT 17 / QUANTUM — pre-registration

Written **before** any arm was run. Every rule below is reported as FIRED / NOT FIRED in
`s17/quantum_FINDINGS.md`.

The search framing is closed (s16 I1–I12). The only live framing is:

> Can a CVaR-VQE sampler occupy a point on the **(member error, diversity)** Pareto frontier
> that is inaccessible to classical samplers at matched computational budget?

The plane is forced by the exact identity `readout² = mean member error² − diversity²`, so a
candidate-set readout is a function of exactly those two coordinates and nothing else.

**Frame discipline (BRIEF §5, instance 2).** `member error` is the **common-frame** term the
identity consumes, not the free-superposition RMSD, and the mean is **quadratic**, not
arithmetic. Both were read wrongly in Sprint 16 and cost ≈2× and 54% respectively. The
identity is verified to machine precision on every set before any claim (`P0`).

---

## P0 — VERIFICATION (no claim precedes it)

* the identity `readout² = M² − D²` holds to < 1e-9 Å² on every constructed set
* `prior` is exactly additive over residues (so its Boltzmann is an exact product law)
* the fixed-temperature Metropolis thermostat respects its evaluation budget exactly
* `stable_rng` everywhere; `hash()` nowhere
* the s16 instrument constants reproduce

## P1 — THE PARETO TEST (mandatory, sprint §20)

**Hypothesis.** No CVaR-VQE configuration reaches an (M, D) point that the classical
attainable set — a fixed-temperature Metropolis thermostat ladder at **matched objective
budget**, plus greedy 1-opt, plus annealing, plus the retrieval product law, plus uniform —
does not already dominate.

**Expected outcome.** VQE dominated everywhere. s16 §6b already showed a diversity-matched
classical thermostat reproduces α's whole ensemble effect (Pearson +0.93/+0.98) and moves the
readout **1.6–1.7× further**; a thermostat that moves support at matched budget should be
strictly stronger than the `tilt_exact` used there.

**Strongest control.** Fixed-T Metropolis at matched budget (a genuine thermostat that moves
support, which s16 never ran — `tilt_samples` cannot leave its support and `tilt_exact` costs
32–128× the budget), swept over a temperature ladder so it traces a *curve*, not a point.

**Success criterion for the quantum arm.** At least one (α, schedule, init) VQE configuration
whose (M, D) point is **not dominated** by any classical point at matched budget, on ≥ 10 of
19 targets, with the paired per-target dominance margin's CI excluding zero.

**Falsifier (pre-registered, from the brief).** If no VQE configuration reaches a frontier
point the classical thermostat plus 1-opt do not already dominate at matched budget, **the
quantum pillar has no accuracy role in this architecture and is retained as a studied object
only.** On existing evidence this is the likely outcome and will be reported as such.

**Dominance rule, fixed in advance.** Point A dominates B iff `M_A ≤ M_B − τ` and
`D_A ≥ D_B + τ` is false for B over A, i.e. A is at least as good on both coordinates and
strictly better on one by more than the tolerance `τ = 0.02 Å` (a quarter of the programme's
0.08 Å false-positive floor). Lower M is better; higher D is better. Reported both ways so
the reader can apply the opposite diversity preference.

## P2 — MATCHED-BUDGET FAIRNESS (sprint §21)

Every arm reports: objective evaluations · circuit samples drawn · gradient passes ·
distinct configurations touched · wall clock · CPU seconds. No arm is placed at nominal
parity without its true cost printed beside it.

## P3 — LOCAL VQE (sprint §23)

**Hypothesis.** A CVaR-VQE run on a *w*-residue local neighbourhood of a candidate produces a
proposal distribution whose (M, D) point beats what classical local search produces on the
same neighbourhood.

**Expected outcome, and the reason this is a trap.** A neighbourhood small enough for a local
VQE to be tractable (w ≤ 6 residues → 4⁶ = 4,096 configurations) is small enough to
**exhaustively enumerate**, at which point the exact optimum and the exact Boltzmann law of
the subproblem are both available classically for ≤ 4,096 evaluations. The reformulation
therefore *removes its own justification*, and I expect to measure exactly that.

**Control.** Exhaustive enumeration of the same neighbourhood (certified), and the exact
Boltzmann law of the subproblem at matched entropy/diversity.

**Success criterion.** Local VQE reaches an (M, D) point outside the exhaustively-certified
classical frontier of the *same* neighbourhood at a cost below 4^w.

**Falsifier.** If exhaustive enumeration of the neighbourhood costs ≤ the VQE's budget, the
local reformulation is closed by construction and reported as closed.

## P4 — ENSEMBLE REWEIGHTING (sprint §24)

**Hypothesis.** The trained circuit's own probability `q_θ(x)` carries ranking information
about true RMSD **beyond** the objective `E(x)` it was trained on — i.e. the ansatz's limited
correlation structure acts as a useful regulariser of the objective.

**Expected outcome.** No. `q_θ` is a low-bond-dimension approximation to a Boltzmann law of
`E`; anything it adds should be reproduced by a **mean-field (fully factorised) Boltzmann
approximation of `E`**, which is exactly computable here.

**Controls, both mandatory.** (i) the mean-field Boltzmann law of `E` at matched entropy —
the correct classical analogue of a low-correlation variational family; (ii) `E` itself as the
ranker (the zero-information-added reference).

**Success criterion.** Spearman(`log q_θ`, RMSD) > Spearman(`−E`, RMSD) *and* > the mean-field
control, on ≥ 13/19 targets with a CI excluding zero, on the in-band subset that a selector
actually consumes.

**Falsifier.** If the mean-field control matches or beats `q_θ`, the variational family
contributes nothing beyond factorised smoothing and §24 is closed.

## P5 — MIXTURES (sprint §19)

**Hypothesis.** The VQE contributes **unique** high-quality configurations that no classical
component supplies.

**Metric.** Not a mean: the count and near-native mass of configurations in the true top 0.1%
that the VQE's visited set contains and the union of all classical arms' visited sets does
not — and the reverse, which is the control.

**Falsifier.** If the reverse count exceeds the forward count, the VQE is a strict subset
contributor and mixing it in is a dilution, not an addition.

## P6 — THEORY (sprint §22)

**Exact measurement, not a sweep.** The Walsh (Pauli-Z) spectrum of the deployed objective on
the full enumerated register, resolved by Pauli weight, raw and rank-conditioned. The memory
ledger records that a Pauli spectrum of an *unconditioned* energy measures its worst clash, so
the conditioned spectrum is the primary and the raw one is reported beside it.

**Two separate questions, never conflated.** (a) Is there a structural reason this problem
class is hostile to current VQE? (b) Is there a reformulation in which VQE becomes useful?

---

## Standing rules for this workstream

* TARGET is the unit; seeds averaged within target before any interval.
* Paired bootstrap with per-fold means and sign consistency printed; medians beside means;
  W/L beside every mean; the ≈0.08 Å false-positive floor printed on every row.
* Both mandatory controls on every accuracy arm.
* No native information in any predictor. Every ORACLE quantity labelled.
* The sealed 60-target benchmark is not read, probed or derived.
* Claim labels: ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN · REFUTED · EXACT · ORACLE.
