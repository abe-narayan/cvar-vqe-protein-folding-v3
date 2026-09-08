# PREREG_B — Sprint 20, WORKSTREAM B (CVaR-VQE / QUANTUM)

**Written before any Sprint-20 arm ran. Not to be edited after seeing results.**
If a specification here turns out mis-specified it STAYS, is flagged, and the untested regime is
recorded OPEN (BRIEF §7).

Salt: `s20qb`. Seeding: `s15/seed.py::stable_rng` only. `hash()` appears nowhere.
The sealed 60-target benchmark and `results/benchmark_manifest.json` are not read, probed or
derived from. `dev24` is not run.

---

## 0. WHAT SPRINT 19'S LANE ALREADY ANSWERED — recovered, not repeated

`s19/results/qb_main.json` + `qb_report_main.json` hold a **complete n = 126 sweep** whose analysis
was never delivered. It is recovered and credited in `s20/agentB_FINDINGS.md` §1 rather than
re-run. Its verdict (all four pre-registered outcomes of `s19/PREREG_B.md` §6) is **CLOSED against
the sampler hypothesis**, and its kill rules fired. **This lane therefore does not re-ask
"is the VQE a better generator on the distogram objective".** It asks Sprint 20's two questions.

---

## 1. THE OBJECT UNDER TEST

Configuration `z = (phi_1..phi_n, psi_1..psi_n) in T^{2n}` — **continuous torsions**, the same
object `core/project.build_ca_exact` builds from and `s19/qb_lib` used. **No k=4 lattice, no
binary encoding of an angle anywhere in Q1.** (Sprint 18's headline was that a 2-bit labelling
carried its own result; there is no "which bit is which" here.)

Three genuine Hamiltonians, evaluated on the *same* `z`:

| symbol | definition | genuineness |
|---|---|---|
| `LEG` | `core.energy.components_batch` summed with **`DEFAULT_WEIGHTS`**, never fitted | genuine Legacy |
| `AMB` | ff14SB / GBn2 **single point**, no minimisation, on the built all-atom structure | genuine AMBER |
| `AMBc` | `sign(E)*log1p(\|E\|)` of **the same** `AMB` call | a **strictly monotone** reparameterisation |
| `DIST` | the deployed distogram functional (reference axis only, never a Q1 claim) | deployed |

`AMBc` is the load-bearing control of this pre-registration. **A strictly monotone transform
leaves the argmin, the complete ranking of configurations, and every level set unchanged.** It
changes only the *geometry* of the scalar field: gradient magnitudes, Hessian conditioning,
dynamic range. So:

> **any difference between `AMB` and `AMBc` is, by construction, NOT a property of the
> optimisation problem's solution set. It is a property of the objective's parameterisation.**
> This is an identity, not a finding, and is labelled EXACT.

`AMB` bit-exactness against the codebase entry point (`core.amber.refine_coords(k_restraint=0,
steps=-1)`) is asserted to `0.00e+00` relative difference on ≥ 4 structures per target before any
number is read; the lean context call exists only to remove ~9x of Python overhead.

---

## 2. Q1 — THE HYPOTHESIS, AND THE FALSIFIER FIRST

### Falsifier (pre-registered, BRIEF §6 rule 5)

> **F1. If, after monotone conditioning that equalises the two energies' marginal distributions on
> a common reference ensemble, `AMB` STILL shows (a) significantly worse optimizer success rate,
> (b) significantly higher negative-curvature fraction or condition number, (c) significantly
> shorter 1-D autocorrelation length in radians — AND (d) at least one of those landscape metrics
> predicts final Cα-RMSD across targets with a rank correlation whose CI excludes zero — then
> "the AMBER landscape is genuinely harder in a way that matters for structure" STANDS and my
> hypothesis is REFUTED.**

> **F2. If `AMB` and `AMBc` produce statistically indistinguishable optimizer outcomes, the
> conditioning explanation is REFUTED** and roughness/ordering must carry the difference instead.

> **F3. If any landscape metric predicts final Cα-RMSD with |rho| CI excluding zero, I may not
> report "landscape geometry is disconnected from structural outcome."**

### Hypothesis (stated so F1/F2 can kill it)

**H1.** The `AMB`↔`LEG` difference in VQE behaviour is dominated by **objective conditioning**,
not by physics: raw `AMB` on ideal-geometry-built backbones is an r⁻¹² clash spike with a
dynamic range of many orders of magnitude, so gradient estimates, CVaR tails and every quasi-Newton
step are set by the single worst contact. Remedy = condition monotonically, not "use better physics".

**H2.** **No landscape metric of either energy predicts final Cα-RMSD**, because (LEDGER,
`objective-does-not-rank-the-native`, `structural-objective-beats-the-energies`) *both* energies'
optima are in the wrong place. Under H2 the correct headline is *"AMBER's landscape is measurably
harder AND that fact is structurally irrelevant"* — which is a different claim from either
"AMBER is harder" or "AMBER is fine".

### Primary endpoints (Q1)

1. **P1-opt** — per (target, optimizer, seed): best-seen objective reached, and **final Cα-RMSD of
   the best-seen configuration** (ORACLE, post-hoc). Contrast `AMB` vs `AMBc` vs `LEG`.
2. **P1-land** — the scale-free landscape panel (§3), per (target, objective, start).
3. **P1-link** — Spearman(landscape metric, final Cα-RMSD) across targets, **reported for every
   metric including the ones that fail**, with bootstrap CIs. This is the BRIEF §4 requirement.

### Nulls and matched controls (Q1)

* **`best_of_N`** — *N* i.i.d. draws from the same per-residue von Mises basin mixtures used to
  initialise every optimizer, at **exactly the arm's evaluation budget**, scored by the same
  objective. **Mandatory.** No arm in this lane is ever compared to an initialisation mean
  (LEDGER: `concentration-is-wrong-when-discrimination-binds`).
* **`c_metro`, `c_anneal`** — classical thermostat / simulated annealing in continuous torsion
  space at matched budget.
* **`greedy`** — coordinate-descent greedy over single-residue proposals, matched budget.
* Zero-information reference for RMSD: the constant ideal α-helix and the matched empirical
  torsion marginals (BRIEF §6 rule 4 — uniform-on-the-torus is **not** used as a control).

---

## 3. THE LANDSCAPE PANEL — every metric declared before measurement

All computed on the **standardised** field `Ê = (E − mu_ref)/sd_ref`, where `mu_ref, sd_ref` are the
mean and sd of that objective over the target's own K = 500 BLOSUM retrieval pool. **This
standardisation is part of the operator** (BRIEF §6 rule 1): without it a gradient-norm comparison
between an energy in units of 10⁸ and one in units of 10¹ is a units comparison, not a physics one.

Central finite differences, step `h = 0.02 rad`, on `z`:

* `gnorm` ‖∇Ê‖, its distribution over starts, and `gvar` its variance across starts.
* Hessian spectrum: `frac_neg` (λ<0), `cond` = λ_max/λ_min over |λ| ≥ 1e-8 λ_max, `gap`
  (λ₂−λ₁)/|λ₁|, `frac_nearzero` (|λ| < 1e-3 λ_max), `aniso` = λ_max / mean|λ|.
  **A single local Hessian is never used to infer global topology** (BRIEF §4) — every Hessian
  statistic is reported as a distribution over ≥ 2 starts × ≥ 8 targets and labelled *local*.
* 1-D ruggedness on 3 random unit directions, 65 points over ±π: `n_localmin` per 2π,
  `acorr_len` (lag at which the detrended autocorrelation first crosses 1/e, in radians),
  `tv_over_range` = total variation / (max−min) — a scale-free roughness index.
* `barrier` — max of Ê along the straight line between two independent starts, minus the larger
  endpoint, in standardised units.
* Basin width / connectivity: fraction of the 6 starts whose local descent reaches within 0.5 Å
  Cα-RMSD of the same terminus (**diagnostic**, ORACLE-free clustering, RMSD used only to define
  "same").

**Every one of these is reported against P1-link. A metric with no predictive relationship to
Cα-RMSD is reported as such and is not interpreted further.**

---

## 4. THE THREE SEPARATION EXPERIMENTS

**S1 — optimizer.** SPSA · Adam on central-FD gradients · momentum SGD on central-FD gradients ·
L-BFGS-B (FD) · Nelder–Mead (derivative-free) · Powell · plus the mandatory `best_of_N`, `c_metro`,
`c_anneal`, `greedy`. **Matched evaluation budget B, counted identically for every arm: one budget
unit = one evaluation of the objective at one continuous configuration, FD probes included.**
Reported: failure rate (arms that do not beat `best_of_N`), best-seen objective, final Cα-RMSD,
seed sensitivity (sd over seeds), SPSA gradient quality = cosine between the SPSA estimate and the
central-FD gradient at the same point, and update-norm distributions.

**S2 — CVaR α.** Genuine CVaR-VQE (`core.quantum`, `MPSAnsatz`, sampled score-function gradient,
`baseline="const"` — the corrected estimator, NOT the recorded `tail` defect) on the basin latent,
at α ∈ {0.05, 0.25, 1.00}, under `LEG`, `AMB`, `AMBc`. Reported: effective sample size
`ESS = (Σw)²/Σw²` of the CVaR tail, latent entropy / state concentration, gradient variance and
bias against the α = 1 full-mean gradient, basin occupancy, candidate diversity `D`, and Cα-RMSD.
**Recorded prior to respect: concentration is the wrong objective when discrimination binds.**

**S3 — encoding.** Identical optimizers, identical budget, identical starts, run in `theta` versus
the smooth periodic embedding `u = (cos theta, sin theta)` with a unit-norm retraction. Physically
equivalent, geometrically different. **No discretisation is used to make any experiment easier.**
Any Q1 conclusion that differs between the two parameterisations is reported as
**encoding-dependent** and is not claimable.

---

## 5. Q2 — DECORRELATION

For each target, each generator family `F`, form the **realisable (coherent) distance error** of the
structure that family actually emits:

    r_coh_F = d(X_F) − d_true              over CA pairs |i−j| >= 2

(For a *generator* this is the whole error — a structure realises its own distances exactly, so
the `s19/a_source` split `r = r_coh + r_inc` has `r_inc ≡ 0`. **That is an identity, stated so it is
not mistaken for a finding.**) Families: the VQE arms (α ∈ {0.05, 1.0}), `best_of_N`, `c_metro`,
`pool500` (retrieval), the deployed distogram fit, and the zero-information helix.

**Primary Q2 statistic:** `rho(r_coh_VQE, r_coh_pool)` per target, pooled with a target-level
bootstrap CI, **against the same-family-different-seed ceiling** measured here, never against zero.

> **Q2 FALSIFIER (from the brief, adopted verbatim): if the VQE candidates align with the pool at
> the cross-architecture level (≈0.75+), quantum sampling is another view of the same bias and the
> generator role closes too.**

Reported **beside** it and never collapsed with it: geometric diversity (mean pairwise Cα-RMSD)
and member quality. **Error diversity is measured separately from geometric diversity** (BRIEF §4).

---

## 6. INSTRUMENT, BUDGET, STATISTICS

**Targets.** A deterministic fold-balanced, length-spread subset of the 126-target tuning
instrument, drawn once with `stable_rng("s20qb", "subset")` and **persisted in the artefact**.
Declared size: **n = 20** for the landscape and optimizer panels (compute-bound; AMBER single point
is ~15 ms/eval and the panel is ~3·10⁵ AMBER evaluations). **Every claim carries its n.** The
instrument's MDE at 80% power is 0.084 Å **at n = 126**; at n = 20 it is ~2.5x larger and **no
Cα-RMSD claim at n = 20 is presented as a validated pipeline improvement.** Landscape statistics
(gradient norms, spectra, autocorrelation) have far higher per-target precision and are reported
with their own CIs.

**Budget.** Every budgeted arm in S1/S2/S3 gets the identical `B` evaluations, FD probes counted.
`B` is fixed in the artefact config before the run and not changed afterwards.

**Statistics.** Target is the unit. Paired target-level bootstrap CI (4000 resamples), fold-aware,
median and W/L beside every mean, no mean alone. Median-vs-mean printed together with the
uniform-effect null percentile (LEDGER: a raw drop-top threshold is not a valid test).
No directional conclusion is drawn from a smoke.

**Labels.** EXACT / ORACLE / ESTABLISHED / SUPPORTED / PLAUSIBLE / OPEN / INCONCLUSIVE /
NOT MEASURED / REFUTED. Every Cα-RMSD is ORACLE, post-hoc, and never enters an inference decision.

**Artefacts.** `s20/results/qb2_*.json` with `_COMPLETE` flags, persisted config and target list,
deterministic seeds. Partials `_PARTIAL_*`, superseded runs `_SUPERSEDED_*`. No Sprint-19 file is
overwritten.

---

## 7. KILL RULES

Downgrade to CLOSED immediately if an arm: loses to `best_of_N` · loses to `c_metro`/`c_anneal`
(a trivial thermostat reproduces it) · loses to a zero-information null · has a CI spanning the
required effect · dissolves under target-level analysis · differs between `theta` and `(cos,sin)`
(encoding-dependent) · needs native information in an inference decision · improves an objective
while worsening Cα-RMSD.

**A clean falsification is the deliverable if that is what the data says.**

---

## 8. WHAT I EXPECT (recorded so the outcome cannot be re-narrated)

I expect **H1 SUPPORTED** (the AMB↔AMBc gap is large and the AMBc↔LEG gap is small), **H2
SUPPORTED** (P1-link null across the panel), and **Q2 to close** (VQE candidates are built from
the same retrieval pool's basins, so their coherent error should sit at or above the 0.75 line).

If instead `AMBc` does not recover `LEG`-like behaviour, F2 fires, conditioning is refuted, and the
roughness/higher-order-coupling branch is the live one.
