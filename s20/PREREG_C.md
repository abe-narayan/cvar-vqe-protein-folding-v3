# SPRINT 20 — PRE-REGISTRATION, WORKSTREAM C (LEGACY / AMBER / PHYSICS)

Written **before any Sprint-20 physics number existed**. Format per `s20/BRIEF.md` §7. Per §6 rule 5
and s19 claim Z8, **the falsifier is registered first in every block** and the hypothesis second.

Genuine Legacy (`core.energy`, eleven components, `DEFAULT_WEIGHTS`, **never fitted**) and genuine
AMBER (ff14SB/GBn2 via `core.amber`/`s17.phys_lib.ConstrainedBox`, CPU, `Threads=1`). No learned
surrogate stands in for either anywhere in this workstream. Seeding via `s15.seed.stable_rng` only.
The sealed 60-target benchmark is not read, probed, derived from, or tuned against.

Instrument: 126 cluster-disjoint tuning targets, full-chain Cα-RMSD (`s12.instrument.ca_rmsd`),
**MDE at 80% power = 0.084 Å**. TARGET is the unit; every comparison is a paired, fold-aware
bootstrap with median and W/L beside the mean.

**ADDENDUM, added 2026-09-07 after the coordinator's first-hour correction to `BRIEF.md` §1 landed
and BEFORE any Sprint-20 number existed.** The brief now requires every arm to state its **basis**:
point cloud (3.048 Å coordinate average — *not* a buildable structure, mean virtual Cα–Cα 2.961 Å
against the physical 3.804 Å) / built chain (3.204 Å, the incumbent) / repaired emission (3.236 Å).
Every table in this workstream carries its basis in the header. Specifically:

* **Q1's ablation table is on the POINT-CLOUD basis** — it is a re-run of the Sprint-19 gate
  instrument, whose readout is the coordinate average, and it is comparable to s19's 3.050 Å column
  and to nothing else. **Its validity vector is therefore reported on the MEMBERS (built chains),
  not on the contracted average**, because a 22% contracted point cloud has no meaningful bond
  strain and quoting one would be the basis error the brief names.
* **Q3's Pareto is on the repaired-emission basis** throughout: input = the coordinate average,
  output = the AMBER-relaxed all-atom structure, and each arm's Δ is against **its own input**.
* **Q2 is on the built-chain basis**: every θ evaluated is an ideal-geometry rebuild.

---

## 0. WHAT I INHERIT AND WILL NOT RE-DERIVE

Sprint 19 closed Legacy in all four **decision** roles (objective term, ranker, `leg_torsion` gate,
steric rejector) and established that gate damage is exactly
`readout² = ‖b_pool‖² + 2‖b_pool‖·ALIGN + ‖Δ‖²`, with score gates raising ALIGN and score-free gates
not. Sprint 19's Z4 **RETRACTED** the "pathological non-terminating minimiser": it was CPU
starvation. **I do not reintroduce an iteration bound.** All AMBER here runs `STEPS = 0`, the
deployed unbounded protocol.

**My mission is not "which model wins."** Both are closed as selectors. The open question is
**what each one knows**, which does not require either to be a good ranker.

---

## 1. MY HONEST PRIOR, STATED BEFORE THE DATA

I expect `rho(E_Legacy, E_AMBER)` to be **weakly positive and highly variable across targets**
(median in 0.2–0.5), because both contain a steric term and both punish expansion, but Legacy's
eleven terms are dominated by contact/compactness statistics that ff14SB does not contain and
ff14SB's electrostatics and GBn2 solvation have no Legacy analogue.

I expect the **disagreement partitions to separate on validity but NOT on Cα-RMSD**: AMBER-preferred
candidates should be stereochemically cleaner (that is what a force field is), Legacy-preferred
candidates should be more compact/pool-typical, and **neither cell should be closer to the native**,
because Sprint 17–19 measured both as non-rankers and s19 C2 says any pool-typicality ordering
selects *toward* the pool's error. If that is what I find, **the honest headline is "complementary
about geometry, jointly uninformative about accuracy"** — and I will report it that way rather than
dressing a validity separation as a step toward RMSD.

For Q2 I expect AMBER's torsion-space landscape to be **locally far rougher** (larger gradient norms,
a much wider Hessian spectrum, a larger negative-curvature fraction) and I expect **that roughness
to have no predictive relationship to final Cα-RMSD**, because s19/s18 established that the
*objective*, not the optimiser, is what binds. **If roughness predicted RMSD I would be surprised,
and that surprise is the thing worth measuring.**

---

## 2. BLOCK Q1 — DO LEGACY AND AMBER CONTAIN COMPLEMENTARY INFORMATION?

**Falsifier F-C1 (registered first).** The complementarity claim is **REFUTED** if *either*:

* **(a) redundancy** — the per-target Spearman `rho(E_Legacy, E_AMBER)` has a median above **0.80**,
  i.e. the two orderings are substantially the same ordering; **or**
* **(b) structural indistinguishability** — the four partitions (both-prefer / Legacy-only /
  AMBER-only / neither) are **not separable** on any measured structural or validity axis beyond a
  **matched-random partition of identical cell sizes**: formally, no axis on which the
  `Legacy-only − AMBER-only` contrast has a fold-aware CI excluding zero **after** the same
  contrast is computed on 200 random partitions of the same sizes and the observed value must lie
  outside the central 95% of that null.

**Hypothesis H-C1.** Legacy and AMBER order the same 75 candidates differently, and the candidates
they disagree about form **structurally coherent classes**: AMBER-only candidates are
stereochemically valid but geometrically atypical; Legacy-only candidates are compact and
pool-typical but stereochemically strained.

**Expected mechanism.** ff14SB/GBn2 scores *local covalent and non-bonded* physics of an all-atom
structure; Legacy scores *coarse tertiary statistics* (contact, compactness, cooperativity) of a
backbone. Those are different functions of the same coordinates, so the disagreement should be
readable as a validity-versus-typicality axis.

**Primary endpoint.** (i) the distribution over 126 targets of per-target `rho(E_Legacy, E_AMBER)`;
(ii) the `Legacy-only − AMBER-only` contrast on each axis of `s16.energy_lib.panel` **kept as a
vector, never fused into a scalar**, plus radius of gyration, plus ORACLE Cα-RMSD of the rebuild,
plus the ORACLE coherent-error alignment `⟨e_i, b_pool⟩/‖b_pool‖` in the s19 common frame.

**Null / matched control.** 200 random partitions per target at the identical cell sizes
(`s15.seed.stable_rng`). Every partition contrast is reported as an excess over that null.

**Data.** `s18/results/down.json` (`complete`, 126 rows): the eleven genuine Legacy components and
the genuine AMBER single point on the **identical** ideal-geometry rebuilds of the shipped top-75.
Gate **GC20a** must reproduce Legacy totals and the ORACLE labels from that artefact at 0.00e+00
before any Q1 number is quoted.

**Budget.** ≤ 45 min, no new AMBER minimisation.

**Promotion criterion.** A structural class on which the two models disagree is reportable as a
**finding about what each model knows** at CI-excluding-zero *and* outside the matched-random null.
It is **not** promoted to any decision role: per Sprint 19 both are closed as selectors, and I
pre-commit that **no Q1 result will be used to build a gate, a ranker or a rejector.**

### 2.1 The mandated ablation table (directive §23)

On the **identical** K = 75 candidates, keep m = 38, emit the deployed coordinate average:

    Random · Distogram · Legacy · AMBER · Legacy→AMBER · AMBER→Legacy

`Legacy→AMBER` = Legacy keeps the best 56 (75·√(38/75) ≈ two equal-ratio stages), AMBER keeps the
best 38 of those; `AMBER→Legacy` is the same two stages in the opposite order, so both composites
apply **exactly the same total selection pressure**. Reported per arm: emitted Cα-RMSD, best member
RMSD, diversity `D`, and the **full validity vector** (clashes at 2.0 Å and 2.6 Å, min heavy
distance, bond strain, angle strain, Ramachandran three-way, cis fraction, chirality) plus both
energies. **Every arm is compared to a matched-random ordering of the same count and the same
number of stages.**

**Registered expectation, so it cannot be moved later:** I expect **every** ordered arm to be at or
worse than matched-random on Cα-RMSD (Sprint 18/19, reproduced twice), and I expect the composites
to be no better than their first stage. If an ordered arm **beats** matched-random past the MDE that
contradicts two sprints and I will say so.

---

## 3. BLOCK Q2 — WHAT MAKES THE AMBER LANDSCAPE DIFFERENT (the physics half)

Continuous torsion space `θ = (φ, ψ) ∈ R^{2n}`, **no lattice, no binary encoding** (BRIEF §5). Both
potentials are evaluated as `E(θ) = E(rebuild(θ))` through the *same* ideal-geometry rebuild, so the
parameterisation is held identical and only `H_Legacy ↔ H_AMBER` changes.

**Falsifier F-C2 (registered first).** The claim "the AMBER landscape is *harder in a way that
matters*" is **REFUTED** if the landscape metrics that differ most between the two potentials have
**no predictive relationship to final Cα-RMSD**: specifically if, across targets, the per-target
Spearman between each landscape metric and the Cα-RMSD reached by an identical optimiser from an
identical start has a CI spanning zero for **every** metric. In that case the correct report is
**"the landscapes differ and the difference does not price the outcome"**, and every spectral number
is a labelled diagnostic and nothing more.

**Second falsifier F-C2b.** If the two potentials' landscape metrics are computed in different units
(kcal/mol vs Legacy's arbitrary weighted score) any raw comparison of gradient norm or eigenvalue is
**meaningless**. Every cross-potential comparison must therefore be made on **scale-invariant**
quantities only — condition number, negative-curvature fraction, spectral-gap ratio, participation
ratio, anisotropy, and the *rank* correlation of gradients — or it is not reported. I register this
as a falsifier of my own analysis because the naive version of this experiment is a unit error.

**Hypothesis H-C2.** AMBER's torsion landscape has (i) a larger negative-curvature fraction,
(ii) a far larger condition number, and (iii) narrower basins, than Legacy's on the same θ; and
**none of the three predicts final Cα-RMSD.**

**Primary endpoint.** For each of `{Legacy, AMBER}` × each start: `‖∇E‖` (reported only within a
potential), the eigenvalue spectrum of the `2n × 2n` torsion Hessian — density, spectral gap,
negative-curvature fraction, near-zero-mode count, condition number over `|λ|`, participation ratio
— then, from multiple initialisations, the number of distinct minima, basin width (torus distance
from start to the minimum it reaches), and straight-line barrier estimates between minima pairs.
**Then: the per-target Spearman between each metric and the final Cα-RMSD from the same start.**

**Distinguishing local from global (registered method).** A single Hessian is a **local** statement
and will be labelled as one. Global claims come only from the multi-start basin/connectivity
statistics, never from a spectrum.

**Null / matched control.** (i) a **rotated-frame null** — the same θ with the structure rigidly
rotated; every quantity here is rigid-invariant by construction so this is a numerical-floor
measurement and is reported with its **maximum**, never its mean; (ii) a **matched-random torsion
perturbation** of the same magnitude as each optimiser step; (iii) for the correlation-to-RMSD test,
a label-permutation null across targets.

**Gate GC20b (declared before use).** The finite-difference torsion gradient must reproduce a
chain-rule gradient (`Jᵀf`, `f` the analytic OpenMM force) to a stated tolerance on ≥ 3 targets, and
the step size `h` must be shown to sit on the finite-difference plateau (energy-difference metric
flat across a decade of `h`). **If GC20b fails, the Hessian block does not run.** Per s19 Z6 I will
report **how many times any tolerance, guard or fallback actually fired**; a pass with zero firings
is recorded as such.

**Budget.** ≤ 3 h AMBER wall on a quiet box, checkpointed per target, on a declared 30-target
subset drawn deterministically from `I.targets()` (`stable_rng`) — **not** the first 30 in list
order, and the subset list is persisted before the run.

**Promotion criterion.** A landscape metric is reportable as *explanatory* only if it (a) separates
the two potentials on a scale-invariant axis and (b) has a CI-excluding-zero relationship to final
Cα-RMSD. A metric passing only (a) is a **labelled diagnostic** and will be stated as
over-interpretable.

---

## 4. BLOCK Q3 — THE AMBER REPAIR PARETO, FINISHED

Sprint 19 left the uncapped full-ladder run at **28 of 126 rows, `complete: false`**. It is
continued unchanged as `s20.c_pareto_finish` into `s20/results/c_pareto.json`; the s19 rows are
merged **only** after their persisted `cfg_hash` is checked against the current module config, and
s19's own file is never written to.

**Falsifier F-C3 (registered first).** "A differently-constrained AMBER buys the same stereochemical
repair for less Cα displacement" is **REFUTED** if **no** arm is Pareto-non-dominated against the
incumbent `k = 30` on the two axes — i.e. every arm that improves Cα-RMSD against its **own gated
input** is worse on at least one axis of the validity vector, with the comparison made on the
**vector** and never on a fused scalar. Sprint 19 already recorded the shape of this trap:
`caonly_k300` is 0.110 Å more accurate at an indistinguishable clash count and **42% more cis
peptide bonds**, and a scalar validity score would have promoted it. **I pre-commit that no arm is
promoted on an aggregate validity number.**

**Hypothesis H-C3.** No such arm exists; AMBER's +0.164 Å repair tax is a property of the operator,
not of the restraint's parameterisation.

**Primary endpoint.** Two axes, per arm, n = 126: Cα-RMSD (and Δ against the arm's own input) and
the full `EL.panel` vector.

**Null / matched control.** The rotated-frame null, reported with its **maximum**. Every gated arm
against **its own gated input** with the exclusion count stated. The AMBER-free `blend` arms are the
zero-physics control for "step back toward the input."

**Convergence gate (declared before use).** `core.amber.convergence_flags`' `converged` flag, per
call, per arm, with the **exclusion count** printed beside every arm mean; arms are reported both
with and without exclusions. Three targets have failed this gate for six sprints and I expect them
again.

**Budget.** ~3.5 h, resumable per target. **Promotion criterion.** Pareto-non-dominance on the
vector, at n = 126, with a CI excluding zero on the accuracy axis and no axis of the validity vector
degraded past its own matched control. Anything less is reported as "no Pareto point exists."

---

## 5. WHAT WOULD MAKE ME ABANDON EACH BLOCK

| block | abandon when |
|---|---|
| Q1 | F-C1 fires (redundant orderings, or partitions indistinguishable from a matched-random split) |
| Q1 ablation | a composite beats its own first stage only on the tuning subset, or only before the matched-random control |
| Q2 | F-C2 fires (no landscape metric prices RMSD) → every spectral number becomes a labelled diagnostic and the block stops |
| Q2 | F-C2b fires (a cross-potential comparison turns out to be a unit artefact) → that comparison is withdrawn, not rescaled |
| Q3 | F-C3 fires → AMBER's standing role is confirmed unchanged and the ladder closes |

**I will not move a goalpost after seeing the data.** If a pre-registration turns out mis-specified I
will say so, keep it unedited, and record the untested regime as OPEN.
