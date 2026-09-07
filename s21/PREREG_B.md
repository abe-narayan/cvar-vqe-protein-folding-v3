# PREREG_B — Sprint 21, WORKSTREAM B (CVaR-VQE MECHANICS)

**Written before any Sprint-21 arm of this lane ran. Not to be edited after seeing results.**
If a specification here is mis-specified it STAYS, is flagged in the findings, and the untested
regime is recorded OPEN (BRIEF §8).

Salt: `s21qb`. Seeding: `s15/seed.py::stable_rng` only; `hash()` appears nowhere.
The sealed 60-target benchmark and `results/benchmark_manifest.json` are not read, probed or
derived from. `dev24` is not run. Target set is `s20.qb2_lib.subset(20)` — the **same 20 tuning
targets** Sprint 20's lane used, so every arm here is paired against a Sprint-20 arm on the same
target.

Representation is **continuous torsions** `z = (phi_1..phi_n, psi_1..psi_n)` in `T^{2n}` throughout.
**No lattice, no binary encoding of an angle.** Hamiltonians are Sprint 20's, imported unchanged
from `s20/qb2_lib.py`: `LEG` (DEFAULT_WEIGHTS, never fitted), `AMB` (ff14SB/GBn2 bare single
point, bit-exactness gate on every target), `AMBc` (`sign(E)log1p|E|` of the same call — a
strictly monotone reparameterisation), `DIST` (the deployed distogram functional).

**Readout, stated because the BRIEF was corrected on this point mid-sprint.** Every arm in this
lane reports the Ca-RMSD of the **argmin over everything evaluated** (`Field.best_z`), the same
`best_seen` readout `core/quantum.py` implements. **No arm here uses an averaged readout.**
Basis: **built chain** from continuous torsions via `core.project.build_ca_exact` — never a point
cloud, never compared to the 3.048 A cloud number.

**4 seeds on every variational arm** (s20 L14: ansatz-seed sensitivity 0.200 A within-target sd =
2.4x MDE). MDE at 80% power = **0.084 A**. Unit of analysis is the **target**; paired target-level
bootstrap CI, median and W/L beside every mean, fold behaviour reported.

---

## 0. WHAT IS SETTLED AND IS NOT RE-ASKED

Priority 2 satisfied (the VQE genuinely trains, -0.210 [-0.339, -0.082], 5/5 folds vs best-of-N
from the *untrained* circuit). CVaR's minimiser is a **face** (EXACT). The deployed ansatz is
**bond dimension 4** (EXACT). The tail parameter alpha is worth nothing and alpha=1 wins. SPSA/central-FD
mismatch is REFUTED. Entanglement is NOT MEASURED with a benefit >0.09 A excluded.
**None of these is re-run.** This lane asks the three open mechanics questions.

---

## 1. BLOCK E — THE ENCODING. FALSIFIER FIRST.

Sprint 20 measured, over 12 (objective, optimiser) cells at n=10, 2 seeds: the embedding
`u = (cos t, sin t)` with an `arctan2` retraction beats plain `t` on final Ca-RMSD in 11/12 cells,
sign-test p = 0.006, up to -0.649 A. Both parameterisations are **physically identical**.

### The hypotheses, named

* **H-E1 (the Sprint-20 lane's reading, and the one I am trying to kill):** the embedding is a
  genuine **optimisation-geometry** effect — smooth periodic chart, no wraparound discontinuity,
  better gradient conditioning. If true it is a real and cheap RMSD lever.
* **H-E0 (mine):** it is a **step-size / probe-scale artefact**. Both `t` and `u` admit a
  one-parameter family of *physically identical* reparameterisations — `t = s*x` for the angle
  chart, `u = r(cos t, sin t)` for the embedding — that change nothing about the problem and
  everything about how far one optimiser step moves in radians. If the effect lives on that
  family, "encoding" is a **gauge**, not physics.

### FALSIFIERS, pre-registered

> **F-E1 (the gauge falsifier — PRIMARY).** The embedding at radius `r` and the angle chart at
> scale `s` are **exact gauge orbits**: `arctan2` is scale-invariant and `t = s*x` is a linear
> chart of the same torus, so the *problem, its solution set and every one of its level sets are
> unchanged* along both. **If final Ca-RMSD depends on `r` within the embedding family, or on `s`
> within the angle family, with a paired CI excluding zero, then outcome is NOT gauge-invariant,
> the Sprint-20 encoding claim RESTS ON AN ARBITRARY ENCODING (BRIEF §6, §12), and H-E1 is
> REFUTED as stated.** The claim is then relabelled **step-size dependence**, and the only
> defensible lever is a step-size rule.

> **F-E2 (the dominance falsifier).** **If the best `t`-family scale `s` in {0.5, 1, 2} matches or
> beats the best embedding radius `r` in {0.5, 1, 2} on RMSD (paired difference CI containing zero
> or favouring `t`), the embedding contributes nothing beyond a scalar already available inside the
> plain angle chart, and H-E1 is REFUTED.** Conversely, **if the embedding at its best radius still
> beats every `t` scale by >= 0.084 A with a CI excluding zero, H-E0 is REFUTED and the embedding is
> a genuine chart effect.** I commit to that reading in advance.

> **F-E3 (the wrap falsifier).** `Field.wrap` re-wraps the iterate to (-pi, pi] after every step for
> `spsa`/`adam_fd` and **never fires** for `nelder`/`lbfgs_fd` (scipy owns the iterate). Wrapping is
> the identity on the objective but **not** on a stateful optimiser's momentum/curvature memory.
> **If `t`-no-wrap recovers the embedding's advantage on the wrapping arms, the mechanism is
> optimiser-state corruption, not periodicity of the encoding.** **If the advantage persists on
> `nelder`, where no wrap ever occurs, the wrap mechanism is REFUTED for that cell** and cannot be
> the general explanation.

> **F-E4 (the retraction falsifier).** Adam/SPSA steps in `u` leave the circle, so the radius
> drifts and the *effective angular step size becomes state-dependent* — an implicit adaptive
> schedule. **If `emb_norm` (radius retracted to `r` after every step, same chart, same
> dimension, same smoothness) loses the advantage, the mechanism is radius drift, i.e. an
> accidental step-size schedule, not the chart.**

> **F-E5 (the harness gate — must fire, not pass vacuously).** `best_of_N` draws from the same
> von Mises basin mixtures in every chart and is mapped back through the same retraction, so its
> emitted torsion set must be **bit-identical** across all 8 parameterisations at the same seed.
> **If it is not bit-identical the harness is not matched and no Block-E number may be read.**
> Reported as a measured max |delta|, not asserted.

> **F-E6 (the RMSD connection).** rho(objective gained, RMSD gained) = +0.003 [-0.068, +0.072] over
> 36 cells. **Any parameterisation that lowers the objective without lowering RMSD is reported as
> a NEGATIVE result.** RMSD is the endpoint; the objective column is secondary and is reported
> with its iteration-count confound attached.

### Design

8 parameterisations, all provably physically identical, all from the identical `z0`:

| id | chart | optimiser coordinate | retraction / post-step | dim |
|---|---|---|---|---|
| `th` | angle | `x = t` | wrap to (-pi, pi] | 2n |
| `th_nowrap` | angle | `x = t` | none | 2n |
| `th_s05` | angle, scale 0.5 | `x = t/0.5` | wrap | 2n |
| `th_s2` | angle, scale 2 | `x = t/2` | wrap | 2n |
| `emb_r1` | circle, radius 1 | `x = 1*(cos t, sin t)` | none (arctan2 at readout) | 4n |
| `emb_r05` | circle, radius 0.5 | `x = 0.5*(cos t, sin t)` | none | 4n |
| `emb_r2` | circle, radius 2 | `x = 2*(cos t, sin t)` | none | 4n |
| `emb_norm` | circle, radius 1 | as `emb_r1` | renormalise each pair to `r` | 4n |

`emb_r1` is Sprint 20's arm exactly. Optimisers `spsa`, `adam_fd` (both wrap, both have explicit
step sizes) on objectives `AMBc` and `DIST`; `nelder` on `AMBc` and `DIST` in a separate block
because Nelder-Mead is **comparison-based and therefore exactly invariant** under both gauge
families and under any strictly monotone transform of the objective (asserted with a measured
gate, not assumed). Budget **B = 512 evaluations, finite-difference probes included**, identical
for every arm. **4 seeds.** n = 20 targets.

**Primary endpoint:** paired target-level mean difference in final Ca-RMSD of the best-evaluated
configuration, seeds averaged within target, negative = the first-named arm better.
**Null:** `best_of_N` at the same budget (mandatory; never an initialisation mean).
**Promotion rule:** a parameterisation is promoted only if it (i) beats `th` by >= 0.084 A with a
CI excluding zero on >= 2 objectives, (ii) is **not** matched by any member of the `t`-scale gauge
family, and (iii) also beats `best_of_N`. **A chart that beats `th` but loses to `best_of_N` is
not a lever; it is a less-bad way of losing to doing nothing.**

---

## 2. BLOCK A — THE ANSATZ AGAINST ITS OWN EXACT LIMIT

`core.quantum.MPSAnsatz(n, layers=L)` sets `chi = 2**L` **by construction in the source** (`D`
doubles once per entangling layer), so the deployed `L=2` is `chi=4`, and depth is the only knob
that moves it. The question is **not** "is deeper better".

### FALSIFIERS

> **F-A1 (the register falsifier — stated first because it damages my own lane's framing most).**
> Every target here has `n <= 16` residues, so the full statevector is <= 65,536 amplitudes and
> **every ansatz at this problem size is classically samplable by brute force, independent of
> bond dimension.** If that holds, "escapes bond dimension 4" is **not sufficient** for a
> quantum-resource claim at this scale and I must say so rather than report chi as if it were.
> Verified by construction and reported as EXACT, not as a discovery.

> **F-A2 (trainability).** Bond dimension is raised by depth. **If the gradient norm of the CVaR
> estimator decays with depth while RMSD does not improve, "expressivity" and "trainability"
> trade and no ansatz in reach is both.** Measured: per-iteration gradient norm, its variance
> across seeds, latent entropy, and RMSD.

> **F-A3 (the endpoint).** **If no ansatz beats `mps_L2` on RMSD by >= 0.084 A with a CI excluding
> zero, ansatz expressivity is closed on this encoding** — a third independent closure after
> Sprint 20's entanglement null and `x_deep` result. This is the outcome I expect.

**Measured, per ansatz, exactly:** bond dimension at every bipartition, by SVD of the exact
amplitude tensor at the trained parameters (n <= 16 makes this exact, not estimated); parameter
count; latent entropy; gradient norm; final Ca-RMSD at **4 seeds**; and the mandatory
`vqe_untrained` and `best_of_N` controls at the same budget.

Ansatz zoo: `prod` (entangler none, chi 1) · `mps_L1` (chi 2) · `mps_L2` (chi 4, **deployed**) ·
`mps_L3` (chi 8) · `mps_L4` (chi 16) · `mps_L2_nofinal` (no trailing RY) · `share_L2` (one RY
angle shared across all residues per layer — `nblocks` parameters instead of `n*nblocks`) ·
`sv_ring_L3` (`StatevectorCircuit`, ring entangler, exact statevector — a topology the open MPS
chain cannot reach).

---

## 3. BLOCK Q — OPTIMIZER, INCLUDING QNG

SPSA-vs-FD mismatch is REFUTED, so the open question is **geometry-aware preconditioning** on the
badly-conditioned AMBER landscape (participation ratio 0.0745, anisotropy 18.8).

The variational parameters are the ansatz angles; the sampled objective is the CVaR of `E(b)`
under `b ~ p_theta`. The right metric for a **sampling** objective optimised by a score-function
estimator is the **classical Fisher information of `p_theta(b)`**, `F = E_p[grad log p grad log p^T]`
— estimated from **the same shots**, costing **zero extra objective evaluations**. This is stated
as a definition. It is **not** claimed to be the Fubini-Study metric of the state; the two coincide
only for the diagonal (computational-basis) part, and that is recorded rather than glossed.

Arms: `adam` (deployed baseline) · `sgd` (plain gradient, matched lr) · `qng` (Fisher-
preconditioned, `(F + lam*I)^{-1} g`) · `qng_diag` (diagonal Fisher) · `spsa_ansatz`.
Objectives `AMBc` and `DIST`. **4 seeds.** Controls: `vqe_untrained`, `best_of_N`.

### FALSIFIERS

> **F-Q1.** **If `qng` lowers the CVaR objective relative to `adam` but does not lower RMSD, it is
> reported as a NEGATIVE result** (BRIEF §3; rho = +0.003). This is the outcome I expect.

> **F-Q2 (landscape vs optimizer — do not conflate).** If **no** optimiser in the battery beats
> `best_of_N` on either Hamiltonian, the limitation is the **landscape/discrimination**, not the
> optimiser, and no optimiser claim may be made. If some optimiser beats `best_of_N` on one
> Hamiltonian and not the other, the limitation is **Hamiltonian-specific** and must be reported
> that way.

> **F-Q3.** If `qng` helps on `DIST` and not on `AMBc`, or vice versa, that is the answer to
> "landscape or optimizer" and must be reported as such rather than averaged away.

---

## 4. WHAT WOULD MAKE ME REPORT A NULL

Any of: an arm loses to `best_of_N`; an effect lives on a gauge orbit; an effect appears only on
the objective and not on RMSD; a CI spans the 0.084 A MDE; an effect dissolves under target-level
analysis; fewer than 4 seeds. **A clean falsification is a successful result.**

## 5. ARTEFACTS

`s21/results/qb3_*.json` with `_COMPLETE` flags that require the **full** configuration
(target list x objectives x arms x seeds), not the subset a call happened to run.
`_PARTIAL_` prefix otherwise. Nothing in `s20/` or earlier is overwritten.

---

## 6. ADDENDUM, 2026-09-07 — BLOCK S: the shots/iterations trade at fixed budget

**Added after Blocks A and Q reported on the OpenMM-free objectives and BEFORE any Block-S arm
ran. Sections 1-5 above are unedited. The falsifier below was fixed before the data existed.**

**Why.** The coordinator's audit lane proposes that the Sprint-20 encoding effect is "take fewer
steps", and notes that this programme's most reproduced law is that optimising the deployed
objectives harder makes the structure worse. That hypothesis has never been tested **directly on
the variational arm**, where the number of gradient steps is a free knob that changes nothing
else: at fixed budget `B = 512`, `shots x iterations = B`, so `shots` in {16, 32, 64, 128} buys
{32, 16, 8, 4} gradient steps of the SAME optimiser on the SAME ansatz from the SAME starts. It is
an 8x sweep of optimisation effort with the encoding, the chart, the dimension, the probe size and
the budget all held exactly fixed. **This is the clean version of the control the encoding
question keeps needing.** It is also the CVaR schedule, which is in this lane's remit.

**Design.** `mps_L2`, `adam`, `alpha = 0.25`, B = 512, **4 seeds**, n = 20 targets, objectives
`DIST` and `LEG` (and `AMBc` if the box permits). Control: `best_of_N` at the same budget.

**Prediction (mine, before the data).** The CVaR training loss falls monotonically with the number
of gradient steps, and Ca-RMSD does not move.

> **F-S1.** **If Ca-RMSD improves monotonically as the step count FALLS (128 shots best, 16 shots
> worst) with a CI excluding zero, the "fewer steps helps" mechanism is SUPPORTED on the
> variational arm**, and the encoding result's step-count explanation gains a direct, independent
> instrument.

> **F-S2.** **If the CVaR loss falls with step count while Ca-RMSD is flat across the whole 8x
> sweep, then "fewer steps" does NOT act through RMSD on this arm**, the step-count explanation
> is not sufficient on its own, and this is one more instrument for
> `search-saturates-discrimination-binds`. I expect F-S2.

> **F-S3.** If the 4-step and 32-step arms differ on RMSD by less than the 0.084 A MDE with a
> zero-spanning CI, the sweep is **NOT MEASURED** and neither F-S1 nor F-S2 may be claimed from
> it; the shot count is then simply not a lever and is reported as such.

**Promotion rule.** None. This block cannot promote anything; it exists to price a mechanism.
