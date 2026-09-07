# SPRINT 21 — PRE-REGISTRATION, WORKSTREAM C
## Hamiltonian design, normalisation, preconditioning and continuation

**Written before any Sprint-21 physics number existed. Falsifier registered first in every block,
with my honest prior attached. Kept unedited; deviations are recorded in `agentC_FINDINGS.md` §D
with their reason and their timing, never repaired silently.**

Timestamp of freeze: see `git log` for this file's first commit; artefacts carry `cfg_hash`.

---

## 0. THE THING THAT MOST DAMAGES MY OWN PRIORITY-ONE EXPERIMENT, STATED FIRST

The directive flags **Legacy → AMBER continuation** as potentially the highest-value idea in the
project. Before I build it I record the three measured results that argue it cannot deliver
Priority 1, so that a null is not a surprise and a positive has to survive them:

1. **`s20` L12 / `c_land_null`.** Minimising either potential in continuous torsion space is a
   *move*, and **72% of the RMSD damage of AMBER's move is its SIZE, not its direction**
   (AMBER +0.6198 raw; matched-magnitude realisable null +0.4438; residual +0.1761
   [−0.1248, +0.5741] **NOT MEASURED**). Continuation changes the *path*. A shorter or
   better-conditioned path is a smaller move, and **a smaller move scores better on this metric
   for reasons that have nothing to do with physics.** Any continuation result is therefore
   worthless without a **move-size-matched null**, and I register that null as the operative
   control, not as a footnote.
2. **`s20` L8 / the programme ledger.** Neither potential ranks candidates by accuracy —
   every `rho(Legacy component, ORACLE d)` is within ±0.07 of zero, and the two ORACLE axes are
   exactly the two that do not separate the potentials. **Search saturates; discrimination
   binds.** Continuation is a *search* intervention. At best it makes a bad objective easier to
   minimise, which the programme has measured as **neutral-to-harmful** on a bad objective.
3. **The λ dial may not be a dial at all.** `E_Legacy ≈ −16` and `E_AMBER ≈ +370` at a good
   start, and `E_AMBER` reaches **+1.9e8** at a real clashing start with a gradient to match.
   In raw units `(1−λ)E_L + λE_A` is AMBER-dominated for every λ above the crossover
   `λ* = ‖∇E_L‖/(‖∇E_L‖+‖∇E_A‖)`, which is **O(10⁻²) at a clean start and can be O(10⁻⁹) at a
   steric one**. A uniform λ grid in raw units is then not an interpolation, it is a step
   function, and **any "abrupt transition" it shows is a units artefact rather than a discovery.**
   This is falsifier F-C2b's class (Sprint 20's own instrument defect) in a new place.

**My honest prior: C1 (the mechanism) partly succeeds, C3 (the RMSD endpoint) fails against its
matched null, and the most durable output of this lane will be the normalisation and the exact
statements about what a transform can and cannot change.** I would rather record that now than
discover it after the fact.

---

## 1. BLOCK N — NORMALISATION. **DECLARED BEFORE ANY RMSD IS READ.**

Nothing may be combined before the components are measured. `s21/results/c_norm.json` records,
per target, on the shipped **top-75 real rebuilds** and on the c_land start set:

    mean, median, sd, MAD, IQR, min, max, p1, p99, skew, kurtosis      (numerical range and tails)
    ||grad|| at the medoid start and its pool spread                    (gradient scale)
    ||H||_2 and the scale-free spectrum vector                          (Hessian scale)
    fraction above 1e4 / 1e6                                            (outlier behaviour)

for **each potential separately** (`BRIEF` §7: separably evaluable at all times).

### The four candidate normalisations, all computed, one declared PRIMARY here

| tag | form | properties |
|---|---|---|
| `raw` | `E` | control; **units mismatch, reported so the artefact is visible** |
| `Nz` | `(E − med_pool)/(1.4826·MAD_pool)` | affine ⇒ gradient exactly rescaled; robust centre/scale; **cannot bound a 1e23 tail** |
| `Ng` | `E / ‖∇E(θ_start)‖` per target | makes λ the true *directional* dial at the start; affine; no tail control |
| **`Nt`** | **`s·asinh((E − med_pool)/s)`, `s = 1.4826·MAD_pool`** | **PRIMARY.** smooth (C^∞), **strictly monotone**, agrees with `Nz` to O(z³) in the bulk, and grows only **logarithmically** in the tail, so the steric singularity is bounded without being deleted |
| `Nr` | within-target rank→normal (`s21/tailprice.py`) | **DECLARED AND REJECTED for continuation, with a reason given before use**: it is a discrete rank on a finite pool, has **zero gradient almost everywhere**, and is **undefined off the pool**. It is a valid *pool-selection* transform and is kept for the selection arms only. |

**`Nt` is declared PRIMARY now, before any RMSD is seen, on these physical grounds and no other:**
AMBER's numerical range on real rebuilds spans **−1170 to +5.5e23** with **53.5% above 1e4**
(`s20` L8); a linear scale therefore reports its worst clash, and this project has already been
burned once by exactly that (`pauli-spectrum-delta-spike-artefact`: an unconditioned energy's
spectrum is its worst clash). `asinh` is the standard bounded-but-signed transform, it needs **no
fitted constant** beyond a robust pool scale, and — unlike a clip — it is **differentiable
everywhere and strictly monotone**, which is what the derivation in §4 requires.

`Nz`, `Ng` and `raw` are computed **alongside** on every arm so the choice is auditable and so a
reader can see whether the conclusion is normalisation-dependent. **If the primary conclusion
flips between `Nt`, `Nz` and `Ng`, the conclusion is reported as NORMALISATION-DEPENDENT and no
arm is promoted.**

**F-N.** *The normalisation is unusable if* `Nt` fails to be monotone on any evaluated pair
(checked exactly), *or if* the pool MAD is zero or non-finite for any target, *or if* the four
normalisations disagree in the sign of the primary endpoint. Firing counts reported.

---

## 2. BLOCK C1 — THE MECHANISM, TESTED BEFORE THE EXPENSIVE ENDPOINT

**Hypothesis H-C1.** Legacy is a compactness model (`s20` L8, 124W/2L). If it carries the state
out of the steric singularity, then at the Legacy-minimised point `θ_L` AMBER should be
**less singular** than at the raw start `θ_0`.

**Prediction.** `E_AMBER(θ_L) ≪ E_AMBER(θ_0)`; `‖∇E_AMBER(θ_L)‖ ≪ ‖∇E_AMBER(θ_0)‖`; and the
**participation ratio of AMBER's Hessian rises** at `θ_L` (curvature spread over more modes).

**Primary endpoint.** `part_ratio(H_AMBER at θ_L) − part_ratio(H_AMBER at θ_0)`, paired over the
30 c_land targets, fold-clustered CI (`s18.phys_lib.paired.ci_fold`).

**FALSIFIER F-C1.** *If AMBER's participation ratio at the Legacy-relaxed point is not higher
than at the raw start with a CI excluding zero, the stated mechanism of the continuation
hypothesis is dead* — Legacy relaxation does not deliver AMBER a better-conditioned region — **and
that verdict stands regardless of what C3's RMSD does.** A C3 gain with F-C1 fired would be an
unexplained empirical effect and must be labelled as such, not as "continuation works because it
de-singularises".

**Matched control.** The move-size-matched `toward_member` geodesic of `s20/c_land_null.py`,
carried to the same torus magnitude as the Legacy minimisation, evaluated on the same three
quantities. **A generic move of the same size may de-singularise AMBER just as well** — that is
the null, and the mechanism claim requires beating it.

**Budget.** 30 targets × (1 Legacy minimisation + 3 AMBER Hessians + nulls). ~1 h.

---

## 3. BLOCK C2 — THE λ SWEEP. **NO ASSUMPTION OF SMOOTHNESS.**

`H(λ) = (1−λ)Ĥ_L + λĤ_A` under each normalisation, `λ ∈ {0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65,
0.8, 0.9, 0.95, 1}` plus, **for `raw` only**, a logarithmic sub-grid `10^{−9..−1}` because §0.3
predicts the whole transition lives there.

Measured at every λ, from the **identical** starts, optimiser and budget as `s20/c_land.py`:

* `‖∇H(λ)‖` and the **gradient variance** across the 5 starts and across the pool;
* the **Hessian spectrum**: `part_ratio`, `neg_frac`, `nearzero_frac`, `cond_med`, `aniso`,
  `spec_skew` (scale-free only — `BRIEF` F-C2b);
* **basin movement**: `θ_moved`, and the torus distance between the λ-minimiser and the
  λ′-minimiser for adjacent λ, which is the **basin-hopping** diagnostic;
* **CVaR concentration**: over the target's top-75 pool, the α = 0.15 tail set of `H(λ)`, its
  Jaccard overlap with the λ = 0 and λ = 1 tail sets, and its mean RMSD;
* **RMSD, first and last.**

**Prediction (registered).** Under `Nt` the sweep is *smooth* in the spectrum metrics and the tail
sets change gradually; under `raw` it is a **step at λ ≈ λ\***, and λ* varies by orders of
magnitude across targets.

**FALSIFIER F-C2.** *If the `raw` and `Nt` sweeps produce the same λ-profile,* my §0.3 claim that
raw λ is a units artefact is **wrong** and I say so. *And if no λ-profile metric moves outside its
own noise band under `Nt`,* the continuation family contains nothing to schedule and C3 is not
worth running at full budget.

**Matched control.** Every λ arm is compared to the **move-size-matched** null at that arm's own
`θ_moved`, because a λ that simply moves less will look better on RMSD for no physical reason.

---

## 4. AN EXACT DERIVATION, STATED BEFORE ITS STATISTICS (`BRIEF` §9)

Let `f` be `C¹` and strictly increasing, `Ĥ = f ∘ H`. Then

1. **Selection is invariant.** `argmin Ĥ = argmin H`, and for every α the CVaR α-tail **SET** of
   `Ĥ` is the α-tail set of `H`. (The CVaR *value* changes; the selected face does not.)
2. **Critical points and their index are invariant.** `∇Ĥ = f′(H)∇H` with `f′ > 0`, so the
   critical set is identical, and at a critical point `∇²Ĥ = f′(H)∇²H`, a positive scalar
   multiple — **so every scale-free spectrum metric, including the participation ratio, is
   EXACTLY unchanged at critical points.**
3. **Away from critical points it is a rank-one update.**
   `∇²Ĥ = f′(H)∇²H + f″(H)∇H∇Hᵀ`. For `asinh` with `z=(H−m)/s`, `f″ = −z/(s(1+z²)^{3/2}) < 0`
   for `z>0`, so at a steric singularity the transform **subtracts** a large rank-one component
   aligned with the clash direction. Whether that spreads or further concentrates the spectrum
   is **an empirical question**, and it is measured, not assumed.

**Consequence, registered now.** *A single-Hamiltonian preconditioner of this class cannot change
what a perfect selector would pick.* Block P is therefore about **optimisation**, and any Block-P
RMSD movement must come from the optimiser reaching a different point — never from "better
ranking". A Block-P arm that claims a ranking gain is mis-specified by this derivation.

---

## 5. BLOCK C3 — THE PRIMARY ENDPOINT: THE STAGED SCHEDULE

**Arms**, identical starts, **identical total function-evaluation budget** (the staged arms split
`MAXITER` across stages; this is checked and the realised counts are reported):

    direct_L        minimise Ĥ_L only                                  (λ = 0)
    direct_A        minimise Ĥ_A only                                  (λ = 1)   THE COMPARATOR
    stage_LA        Ĥ_L then Ĥ_A                                       2 stages
    stage_L5A       Ĥ_L then H(0.5) then Ĥ_A                           3 stages
    ramp            λ = 0 → 1 over 5 equal stages                       5 stages
    NULL_move       toward_member geodesic matched to each arm's own total θ_moved

**Primary endpoint.** `RMSD_end(stage_LA) − RMSD_end(direct_A)`, paired over 30 targets,
fold-clustered CI, median and W/L beside the mean.

**Secondary, and the one that decides promotion.** The same difference **after subtracting each
arm's own move-size-matched null**: `[RMSD_end − NULL_move]_staged − [RMSD_end − NULL_move]_direct`.

**FALSIFIER F-C3.** *Continuation is dead as a Priority-1 intervention if the move-size-corrected
staged-minus-direct difference has a CI spanning zero,* or if it is positive. A raw staged-minus-
direct gain that vanishes under the null correction is reported as **MOVE-SIZE, NOT PHYSICS** —
the exact verdict Sprint 20 issued against its own headline.

**Power, stated in advance.** n = 30 targets, not 126. The programme MDE at 80% power on the full
instrument is 0.084 Å; at n = 30 with the observed between-target sd this block is powered for
roughly **0.3–0.5 Å**, so **a null here is NOT MEASURED, not "matched"**, and I will say so.

**Promotion rule.** No promotion to the pipeline from this block under any circumstances. A
surviving positive is escalated to the coordinator as a *candidate* for an n = 126 run by
Workstream A, with the move-size null attached.

---

## 6. BLOCK P — PRECONDITIONING. `H_AMBER` IS `E ∘ Relax` AND STAYS SEPARABLE.

**The constitutive fact, restated so no arm forgets it.** The deployed `H_AMBER` is
`E ∘ Relax₅₀`; the cap binds **192/192**; unrelaxed AMBER is not finite on **42%** of the k=8
register and is *finite-but-meaningless* on **53.5%** of the real top-75 domain. This lane's
landscape work uses the **bare single point** (`ConstrainedBox.energy_point`, no minimisation,
the `core.amber._run(steps<0)` path verified bit-exact over 150 comparisons in Sprint 20), and
**every table says which object it is measuring.**

**Arms** (all on AMBER alone, so §4 applies — these are optimisation interventions):

    P0  raw AMBER                                       control
    P1  Nt-transformed AMBER                            bounded barrier, monotone, EXACTLY rank-preserving
    P2  trust region: step capped at r rad/coord        directly attacks "72% of the damage is move size"
    P3  term-staged AMBER: bonded groups first,         a genuine SUBSET of the AMBER force groups.
        then all terms                                  **Named `amber_bonded`, NEVER "AMBER".**
    P4  Legacy-preconditioned (= C3's stage_LA)         the continuation, re-read as a preconditioner

**FALSIFIER F-P.** *A preconditioning arm is rejected if it does not beat P0 on its own stated
objective* (final `E_AMBER` at matched evaluations) — a preconditioner that does not precondition
is nothing — **and separately** if it beats P0 on RMSD only through a smaller move (the same null
as C3 applies).

**Surrogate rule, binding.** If any arm in this lane evaluates something that is not genuine
ff14SB/GBn2 on the assembled structure, it is named with a `surrogate_` prefix **in the code, in
the artefact keys and in every table**, and it is validated against true AMBER (rank correlation
and paired energy error on ≥ 500 real rebuilds) before any conclusion uses it. **No arm named
`amber` may be anything but genuine AMBER.** P3 is a *force-group subset*, not a surrogate, and is
named `amber_bonded` for exactly that reason.

---

## 7. GATES — RUN AND REPORTED BEFORE ANY SCIENTIFIC NUMBER

| gate | check | rule |
|---|---|---|
| **GC21a** | `s20.c_land.Pot` reproduces `s20/results/c_land.json`'s `start` block bit-for-bit on 3 targets | if it does not, the instrument has drifted and nothing downstream is comparable to Sprint 20 |
| **GC21b** | the FD steps `H_GRAD=1e-3`, `H_HESS=1e-2` still sit on their measured plateaus **for `H(λ)`**, not only for the endpoints — a mixed Hamiltonian has a different noise floor | plateau table printed; **firings counted** |
| **GC21c** | `Nt` monotonicity on every evaluated pair; MAD finite and positive on every target | exact check, firing count |
| **GC21d** | the benchmark seal: SHA-256 of `results/benchmark_manifest.json` recorded **without reading it** | integrity certificate, never a read |

**Every gate reports how many times its tolerance actually FIRED** (`BRIEF` §7 rule 3). A pass with
zero firings certifies the chosen setting, not the tolerance, and is labelled that way.

---

## 8. WHAT WOULD MAKE ME ABANDON THIS LANE

F-C1 fires **and** F-C2 fires: the mechanism is absent and the λ family is featureless, so the
priority-one experiment has no object. In that case the lane's deliverable is the normalisation
measurement, the §4 derivation and a clean falsification — **which `BRIEF` §12 calls a successful
result, and I will report it as the headline rather than burying it.**
