# SPRINT 18 — QUANTUM/ADVERSARIAL PRE-REGISTRATION

Written before any variational arm was run. Timestamped by its commit; the arm code
(`s18/q_cond.py`) did not exist when this file was written.

---

## 1. What is already known, and is therefore NOT part of the pre-registration

The following were measured before this file and are stated here so that no reader mistakes
them for predictions:

* Sprint 17, `s17/quantum_FINDINGS.md` §0: greedy 1-opt certifies the **full** deployed
  objective's global optimum in 100 % of 152 cells at **1,024** objective evaluations against
  a CVaR-VQE budget of **8,192**; deleting the CNOTs from the ansatz is null at every α.
* This sprint, `s18/q_report anova`: the deployed objective carries **8.7 %** of its Walsh
  variance in inter-residue couplings. The strict Walsh weight-≤1 object (`W1`) and the
  residue-additive ANOVA object (`RA`) each carry **exactly 0.0 %** — they are separable by
  construction, and their global optima are closed-form per-qubit / per-residue argmins at
  ~38 table reads.

**So the conditional's answer is already forced on the objectives this sprint is testing.**
The arms below are run anyway, because a forced answer that is never measured is an assertion,
and because the run is what makes the statement checkable. This is recorded as a prior, not
discovered as a result.

---

## 2. THE PRE-REGISTERED CONDITIONAL

> If the new objective contains useful higher-order correlations, then
> **(a)** greedy 1-opt should **cease** to certify its optimum at tiny budget, and
> **(b)** the CNOT-free product ansatz should become **measurably worse** than the entangled
> VQE at matched budget, with a **target-level interval excluding zero**.

**Decision rules, fixed now:**

| clause | measured as | fires FOR the quantum pillar if | fires AGAINST if |
|---|---|---|---|
| (a) | % of (target × seed) cells where greedy 1-opt reaches the CERTIFIED global optimum, vs budget | certification at 1,024 evals is **< 90 %** on the new objective *and* materially lower than on the full objective | certification at 1,024 is ≥ 90 %, or **higher** than on the full objective |
| (b) | `mps2f − mps2fn` paired over targets, matched budget / estimator / optimiser / seed, on M75, D75 and the readout | mean **negative** (entangled better) with a fold-clustered 95 % CI excluding zero | CI spans zero, or the sign favours the product ansatz |

**The quantum pillar earns a role only if** it reaches a Pareto point in objective quality /
diversity / near-native coverage that classical controls cannot reach at matched budget; **or**
(b) fires in its favour; **or** quantum candidates improve downstream RMSD where matched
classical controls do not. Absent all three, it has no role. **A manufactured advantage is
worse than no result.**

---

## 3. Arms and controls

**Objectives.** `full` (deployed, rank-uniformised), `RA` (residue-additive ANOVA — the object
the brief's §4 defines and the only one that ports to continuous torsions), `W1` (strict Walsh
weight-≤1 — the object that produced 2.411 Å).

**Search arms for (a).** greedy 1-opt · fixed-temperature Metropolis · uniform random, all
budget-capped by the same `Counter`, over budgets {20, 36, 64, 128, 256, 512, 1024, 2048,
8192}, 19 targets × 4 seeds.

**Variational arms for (b).** `mps2f` (entangled, CNOT chain) and `mps2fn` (**identical circuit
with the CNOTs removed** — a product Bernoulli model, i.e. a classical algorithm) at α ∈ {0.05,
0.25, 1.00}, identical budget (8,192), shots (512), learning rate (0.15), `const` baseline,
Adam, and identical seeds. Both consumption modes reported (uniform draw and lowest-objective
selection), never substituted.

**Both mandatory controls on every arm.** Zero-information reference: uniform random over the
register. Matched-random: an equal-count uniform draw from the same register.

---

## 4. Falsifiers for this workstream

**Q1** the new objective's greedy certification is ≥ 90 % at ≤ 1,024 evaluations → there is no
optimisation role. **Q2** `mps2f − mps2fn` spans zero → entanglement has earned nothing.
**Q3** the new objective is *more* separable than the full one → the truncation makes the
quantum case worse and that must be stated plainly, not omitted.

**If Q1 and Q2 both fire, the quantum pillar stays closed and no further arm is run.**

---

## 5. Discipline

TARGET is the unit; paired fold-clustered bootstrap CIs with medians and win/loss beside every
mean. No min-of-N ceiling without its min-of-N null, and the null names its band. Pairwise
ordering accuracy, Spearman, Pearson and copula ρ are distinguished by name and never written
as "ρ = accuracy". No native information in any predictor, parameter, threshold or stopping
rule; every RMSD is an ORACLE post-hoc score. `s15/seed.py::stable_rng`, salt `s18quantum`;
`hash()` nowhere. Unique result directory, completion flags, persisted coefficients and config
hashes. **The 60-target sealed benchmark is not read, not probed and not derived.**
