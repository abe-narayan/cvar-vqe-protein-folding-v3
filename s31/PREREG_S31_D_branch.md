# PREREG — S31 lane D — THE λ=0 BRANCH SELECTION IS MADE ON NOISE. IS FIXING IT WORTH ANYTHING?

**Registered 2026-09-21 00:00 (lane D), before any arm of this experiment was computed.**
The only numbers in this document at the time of writing are the *diagnostic* ones from
`s31/results/s31_D_projection_pin.json` (defect D-B), which measure the operator's
conditioning and say nothing about accuracy.

---

## 1. The finding this rests on (D-B, already established)

There is **no RNG** on the cloud→built-chain projection path. `core.project.fit_multi` loops
over four fixed starts and takes a strict argmin; reprojecting the same cloud twice is
**bit-identical** (max |ΔCA| exactly 0.0). The instrument's irreproducibility is not
stochasticity — it is **conditioning**:

- `core.project.lam_path` solves at **λ=0** from four generic starts and takes the argmin.
- Those four starts converge to objectives that agree to ~**1e-7** — numerical noise — while
  sitting on **different torsion branches**. (1A13: 0.5091924865 / 0.5091925170 /
  0.5091924488 / 0.5091924188.)
- The winner of that noise-level argmin becomes the **warm start** for the λ=0.3 rung, where
  the same four branches are **1e-1 apart** (1A13: 0.947 / 0.519 / 0.562 / 2.396).
- So a decision taken where the objective **cannot discriminate** determines an outcome where
  it **can**. Measured over the targets scored so far: median λ=0 margin ~7e-7, and a 1e-14
  relative cloud perturbation moves the emitted chain by a median of 1.6e-3 Å with a
  0.5 Å tail.

**The λ=0 branch selection is a selection made on noise.** That is a defect in the operator
independent of whether fixing it improves accuracy.

## 2. Hypothesis

**H1 (accuracy).** Carrying all four λ=0 branches forward and taking the argmin at λ=0.3 —
where the objective separates them by ~1e-1 instead of ~1e-7 — lowers the mean built-chain
CA-RMSD over the 126 targets.

**H0.** It does not. A lower value of the λ=0.3 objective is not a better structure, so
removing the noise-selection changes which branch is emitted without improving accuracy.

**H0 is the position I expect to defend.** This project's standing finding is that *the
objective does not rank the native* (native is the distance objective's argmin on 3/126
targets, 36.8th percentile) and that *searching harder on a bad objective hurts*. The B4 arm
is, by construction, **a harder search on the same objective** — so the prior from this
project's own record is that it is neutral-to-harmful for accuracy. The **conditioning**
argument for B4 stands regardless of the accuracy outcome, and the two must not be conflated.

## 3. Arms

All arms are projected **in the same job, from the same stored clouds**
(`s29/results/s29_O_structs/<pdb>.npz`, key `prod`), per the D-B rule. The cloud is the
input and is bit-identical across arms; only the projection differs.

| arm | definition | deployable? |
|---|---|---|
| **PROD** | `s12.instrument.project` exactly as production: `lam_path(pen, (0.0, 0.3), multi=True, grad="exact")`. Candidate set at λ=0.3 = {continuation from the λ=0 argmin branch} ∪ {4 fresh starts at λ=0.3}. | incumbent |
| **B4** | Continue **each of the four λ=0 branch solutions** independently to λ=0.3, union the 4 fresh λ=0.3 starts, take the **argmin of the λ=0.3 objective** over all 8. | **native-free, deterministic, DEPLOYABLE** |
| **ORACLE_B4** | Same 8-candidate set; take the argmin of **CA-RMSD to the native**. | **ORACLE / NOT DEPLOYABLE** — a ceiling, never a method |

### A structural fact registered before the result

**B4's candidate set is a strict superset of PROD's**, and both take the argmin of the *same*
λ=0.3 objective. Therefore **B4's objective is ≤ PROD's objective on every target, by
construction.** B4 cannot lose on the objective. It can lose on RMSD, and whether it does is
exactly the open question — which is what makes this a real test rather than a guaranteed win
dressed as one. Any report of this experiment must state that B4 dominates on the objective by
construction, so that "B4 achieves a lower objective" is never presented as evidence for H1.

## 4. Primary comparison, decided now

**Primary:** `B4 − PROD`, **mean built-chain CA-RMSD, n = 126**, paired, via
`s24.stats_lib.compare` (**lower is better; d = a − b, negative = B4 better**), fold-clustered
CI on the pinned folds, **MDE = 2.8016 × SE**.

**Decision rule, fixed before the number:**

- effect ≤ −1.0 × MDE → **B4 improves accuracy.** Propose adoption; do **not** adopt
  unilaterally, because it changes the endpoint and the canonical 3.2105 Å.
- −1.0 × MDE < effect ≤ −0.7 × MDE → **NOT MEASURED.** Not a result. Report as such.
- |effect| < 0.7 × MDE → **H0 holds on accuracy.** B4 is then a **conditioning** fix only, and
  must be described as one, in the same sentence as its null accuracy result.
- effect ≥ +0.7 × MDE → **B4 is WORSE.** Report it; it would be a clean confirmation that a
  harder search on this objective hurts, which is a result worth having.

**Ceiling first.** `ORACLE_B4 − PROD` is computed and reported **before** B4 is interpreted.
**If the ORACLE ceiling is itself below 0.7 × MDE, the direction is closed** and B4's accuracy
number is not pursued further — there is nothing in the branch set to capture. The ORACLE arm
is labelled **ORACLE / NOT DEPLOYABLE** in the same sentence as its number, every time.

## 5. Secondary, pre-specified

1. **Conditioning, the point of the exercise:** the decision margin of the winning candidate
   over the runner-up, PROD (at λ=0) vs B4 (at λ=0.3). Prediction: median margin rises from
   ~1e-7 to ~1e-2 or better. This is a *measurement of the operator*, not a hypothesis test.
2. **The tail:** per-target |Δ| between B4 and PROD, and how many targets change branch at all.
   Registered expectation: most targets do not change, and the effect — in either direction —
   is carried by the minority that do. Report median vs mean and the drop-top null, per the
   standing median-vs-mean rule.
3. **FAIL18 vs the other 108**, reported whichever way it falls.

## 6. What would make me wrong, and what I will not do

- I will **not** adopt B4 into production in this sprint on my own authority. It changes the
  endpoint; that is the coordinator's call on the evidence.
- I will **not** tune λ, the start set, `maxiter`, or the tie tolerance on native RMSD. The
  arms above are the whole experiment; there is no grid.
- I will **not** report B4's lower objective as support for H1 (see §3).
- If B4 is worse, that is the entry's headline.

## 7. Multiplicity

k = 2 primary comparisons on one basis (B4 vs PROD, ORACLE_B4 vs PROD; the built chain — the
CA cloud is the *input* and is identical across arms, so there is no second basis here).
Registered in `s31/MULTIPLICITY.md`.
