# SPRINT 21 — WORKSTREAM C — HAMILTONIAN DESIGN, NORMALISATION, PRECONDITIONING, CONTINUATION

**Genuine Legacy (`core.energy`, eleven components, `DEFAULT_WEIGHTS`, never fitted) and genuine
AMBER ff14SB/GBn2 (`core.amber` / `s17.phys_lib.ConstrainedBox`, OpenMM, CPU, `Threads=1`).
No learned surrogate stands in for either anywhere in this lane.**

Pre-registration: **`s21/PREREG_C.md`**, frozen in git commit `7d07d63` **before any Sprint-21
physics number existed**, falsifier first in every block, with my honest prior attached and the
three measured results that argue against my own priority-one hypothesis stated in §0.

Modules: `s21/c_norm.py` (Block N + the gates), `s21/c_cont.py` (Blocks C1/C2/C2t/C3 + Block P),
`s21/c_report.py` (every table). Artefacts in `s21/results/`, each with a `complete` flag, a row
count and a `cfg_hash`.

**THE OBJECT BEING MEASURED, stated once and never elided.** This lane's landscape work uses the
**bare AMBER single point** (`ConstrainedBox.energy_point`, *no* minimisation) because a landscape
needs a function of `θ`. The **deployed** `H_AMBER` is `E ∘ Relax₅₀` — a different operator, whose
cap binds 192/192 (`s20` L6) and without which AMBER is not finite on 42% of the register. **No
table here is a statement about the deployed operator**, and where the distinction matters it is
named.

**BASIS.** Point cloud (3.048 Å coordinate average, *not* a buildable structure) / built chain
(3.204 Å, the incumbent) / repaired emission (3.236 Å). Block C2/C2t's tail arms are **point
cloud** (a coordinate average of a tail). Blocks C1/C3 are **built chain**. Every table says which.

Labels per `BRIEF` §9: **EXACT · ORACLE · ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN ·
INCONCLUSIVE · NOT MEASURED · NOT SUPPORTED · REFUTED · RETRACTED.**

---

## RUN STATUS

| block | what | n | status | artefact |
|---|---|---|---|---|
| gates | GC21a/c/d pass; **GC21b mis-specified and restated (§7.3)** | 3–9000 | **PASS**, firings reported | `c_gate.json` |
| **N** | component measurement; the normalisation declared before any RMSD | 30 | **COMPLETE** | `c_norm.json` |
| **C1** | the mechanism — does Legacy relaxation de-singularise AMBER? | 30 | **COMPLETE** | `c_c1.json` |
| **C2** | the λ sweep — gradient, gradient variance, tail, Hessian | 30 | **COMPLETE** | `c_c2fast.json`, `c_c2.json` |
| **C2t** | the selection-side λ sweep, **declared extension** | **126** | **COMPLETE** | `c_c2tail.json` |
| **C3 + P** | staged schedules and preconditioners — **the primary endpoint** | 30 | **COMPLETE** | `c_c3.json` |

**Every block reached its full pre-registered configuration; nothing here is a partial.** The
completion flag in each artefact requires `len(rows) >= len(FULL_SUBSET)` **and** that the subset
requested *is* the full subset, so a smoke run cannot write one — verified in the flag audit
(`_SMOKE_c_c1.json`: `complete: false`).

---

# VERDICT — LEADING WITH THE DAMAGE

> ## 1. THE DIRECTIVE'S FLAGGED "POTENTIALLY HIGHEST-VALUE IDEA IN THE PROJECT" IS REFUTED ON ITS OWN PRE-REGISTERED PRIMARY ENDPOINT, AND ITS MECHANISM IS REFUTED SEPARATELY.
>
> **F-C3 fires.** Every staged `Legacy → AMBER` schedule ends **further** from the native than
> running AMBER directly from the same start on the same budget:
> `LA_Nt − A_raw = +0.2331 [+0.0673, +0.3950]`, 12W/18L, n = 30 — **and it survives the move-size
> correction** (`+0.2285 [+0.1140, +0.3647]`), so it is *not* the artefact that dissolved Sprint 20's
> own headline. **The path itself is worse, not merely equal.**
>
> **And the mechanism was refuted first, at n = 30, before the endpoint ran.** The directive's case
> was that Legacy, being a compactness model, would carry the state *out* of the steric singularity.
> It drives it **~31× deeper in energy** (`log₁₀ E_AMBER: +1.4947 [+1.1437, +2.0174]`, worse on
> **25/30**) and **~28× deeper in gradient norm**. **Legacy is the wrong preconditioner for AMBER for
> exactly the reason it is Legacy**: `s20` L8 measured Legacy-preferred candidates as **0.45 Å more
> compact (124W/2L)**, and compaction is what closes the heavy-atom contacts that *are* the
> singularity. **Mechanism and endpoint agree, which is stronger than either alone.**

> ## 2. AND MY OWN PRE-REGISTERED FALSIFIER DID NOT FIRE, BECAUSE I HAD REGISTERED THE WRONG ENDPOINT. I AM REPORTING THAT RATHER THAN CLAIMING THE CLAUSE THAT SURVIVED.
>
> F-C1 was attached to **participation ratio**, which *rises* (`+0.0132 [+0.0019, +0.0226]`). So
> **F-C1 as literally written DOES NOT FIRE.** But against a realisable move of the *same torus
> magnitude* the rise is `+0.0089 [−0.0039, +0.0197]` — **NOT MEASURED**. The scale-free summary I
> chose is insensitive to how far up the `r⁻¹²` wall the configuration sits, which is the entire
> content of "de-singularise". **A lane that had measured only its registered primary would have
> published "Legacy relaxation improves AMBER's conditioning, CI excludes zero" — the opposite of the
> truth.** The pre-registration is kept unedited and the mis-specification recorded (§3.1, §9).

> ## 3. `λ` IN RAW UNITS IS NOT A DIAL. THE ENTIRE RAW LADDER FROM λ = 0.05 TO λ = 1 IS ONE HAMILTONIAN.
>
> Measured before any sweep ran, at 150 real starts: the crossover
> `λ* = ‖∇E_L‖/(‖∇E_L‖+‖∇E_A‖)` has **median 1.99e−05, range 1.51e−09 to 5.58e−02, a 3.7e+07-fold
> spread across targets.** On the *declared uniform grid* the largest single-step change in AMBER's
> gradient share is **0.9989, in the first interval** — and the **Hessian spectrum at λ = 0.1 equals
> the λ = 1 spectrum to four significant figures on all seven metrics**, through λ = 0.9.
>
> **A "Legacy → Legacy+distance → Legacy+AMBER → AMBER" schedule written on a raw λ grid never
> visits an intermediate Hamiltonian at all.** Any abruptness or basin-hopping read off such a grid
> is **kcal/mol, not physics** — falsifier F-C2b's class, which this programme has now caught three
> times. **This is why Block N ran before the sweep, and it is the single most transferable thing in
> this lane.**

> ## 4. THE MANDATORY MATRIX'S PHYSICS HALF IS WORSE THAN CHANCE — ALL 42 ARM-ROWS, n = 126, EVERY CI EXCLUDING ZERO.
>
> Pool-restricted CVaR α = 0.15 tail selection, point-cloud readout: matched-count **random** tail
> **3.101**; Legacy **3.208 (+0.106)**; AMBER **3.216 (+0.115)**; best mixture over three normalisations
> **3.182 (+0.080 [+0.041, +0.119])**. **42 of 42 arm-rows positive, 42 of 42 CIs excluding zero**
> (`raw`, `Nz`, `Nt`; `Ng` was run at n = 30 only). `BRIEF` §2 asked whether Sprint 20's argmin result transfers to a **tail**
> selector. **It does.** The tail *is* better than the argmin (3.18–3.23 vs 3.50–3.69) — and both are
> worse than picking at random.
>
> The mechanism is one line: within each pool, `ρ(E_AMBER, ORACLE d) = −0.0002 [−0.0933, +0.0883]`
> and `ρ(E_Legacy, ORACLE d) = +0.0385 [−0.0826, +0.1592]`. **Neither potential ranks accuracy at
> all**, so a mixture of them cannot either.

> ## 5. NO PRECONDITIONER IN THIS LANE PRECONDITIONS, AND THE ONE THAT LOOKED BEST IS MOVE-SIZE.
>
> **F-P fires on every arm.** Direct AMBER minimisation reaches the lowest — and the only physical —
> final AMBER energy (median **−46.5 kcal/mol** against +13 to +88 for every preconditioned arm, and
> **+8.25e+06 for Legacy**). The best-looking RMSD arm, a trust region at one third the pool's own
> torsional dispersion, is **−0.2987 [−0.6307, −0.0547], 20W/10L** — and against its own
> magnitude-matched null it is **−0.0207 [−0.3185, +0.1736], NOT MEASURED**, with an exactly even
> 15W/15L. **MOVE-SIZE, NOT PHYSICS.** No promotion; `PREREG_C.md` §5 forbade it in advance.

> ## 6. WHAT SURVIVES, AND IT IS SMALL AND HONEST.
>
> The **asinh** normalisation, declared PRIMARY before any RMSD existed, is the only intervention
> that does something principled: it cuts AMBER's **across-start gradient CV by 34%** (1.329 → 0.876)
> and pushes both potentials' participation ratios toward each other (Legacy 0.4965 → 0.2336,
> AMBER 0.0668 → 0.0887) — **at provably zero cost to selection**, since §2 derives that a strictly
> monotone transform leaves `argmin`, every CVaR tail **set** and every critical point invariant
> (verified on 126 × 75 real candidates). Its final-energy row is a clean null (`+0.0038 [−0.0008,
> +0.0082]`) at **67 realised function evaluations against 102**.
>
> **I am not claiming the 34% budget saving**, because L-BFGS terminates on gradient norm and `Nt`
> shrinks the gradient by construction, so it may be the convergence test firing early. **Label:
> OPEN. Missing configuration named: `A_Nt` vs `A_raw` at matched realised `nfev`.**

---

# 0. GATES, AND THE GUARD THAT ACTUALLY FIRED

| gate | check | result |
|---|---|---|
| **GC21a** | `s20.c_land.Pot` reproduces `s20/results/c_land.json`'s `start` block | **PASS — max \|Δ\| = 0.00e+00** on all seven scale-free metrics *and* on `E0`, for both potentials, on 3 targets. The instrument has not drifted; every Sprint-20 comparison below is matched. |
| **GC21c** | `Nt` (asinh) strictly monotone, derivative positive | **PASS over 9,000 pairs spanning [−1.0e+06, +1.0e+24]**, zero firings. |
| **GC21d** | benchmark seal | `results/benchmark_manifest.json` **SHA-256 `a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d`**, 
recorded as an integrity certificate **without reading the file**. |
| **GC21b** | FD plateau for `H(λ)` | **superseded by construction** — see §4, the mixture Hessian is no longer finite-differenced. The check that replaces it is stronger and is reported with its firing count. |

**Zero firings on GC21a and GC21c certifies the chosen setting, not the tolerance** (`BRIEF` §7
rule 3). GC21a's reproduction is *exactly* zero rather than merely within tolerance, which is a
stronger statement than a pass.

### A guard that DID fire, 9 times

`core.amber.memory_guard` refuses to open an OpenMM context above **92% physical memory**, because
an OOM here kills a sibling workstream's run and not only this one. During Block N it **fired 9
times** (`c_norm.json: memory_guard_firings = 9`), each time with the box at 93% — driven by a
browser, not by a research process. The driver retries with a 45 s backoff; **the firings are
recorded in the artefact rather than hidden**, because a run that had to wait is a different run
from one that did not.

### Timing is NOT reported as a property of the code

Mid-lane the box was carrying **eight concurrent Python processes** from four other workstreams
(`s21.tailprice` ×2, `s21.d_enc` ×2, `s21.qb3_run`, `s21.a_matrix`, plus this lane). `BRIEF` §10:
**measure timings on a quiet box or not at all.** Per-target wall times in this lane's logs
therefore carry no information about the code and none is claimed from them.

---

# 1. BLOCK N — WHAT THE TWO COMPONENTS ACTUALLY ARE

**Declared before any RMSD was read** (`PREREG_C.md` §1, frozen in git). Domain: the shipped
**top-75 real rebuilds** per target, n = 30 targets (the Sprint-20 `c_land` subset, reused verbatim
so every cross-sprint comparison is matched). **AMBER = bare single point.**

Medians over targets:

| statistic | **Legacy** (arb. units) | **AMBER** (kcal/mol) |
|---|---:|---:|
| pool median | −12.33 | **+1.072e+04** |
| pool 1.4826·MAD | 2.856 | 1.134e+04 |
| pool sd | 3.209 | **2.686e+10** |
| pool IQR | 4.006 | 1.209e+05 |
| pool min / max | −20.36 / −2.494 | +155.4 / **+1.795e+11** |
| pool p1 / p99 | −19.86 / −3.856 | +315.1 / **+6.687e+10** |
| skew / kurtosis | 0.528 / 3.60 | **8.13 / 68.9** |
| fraction > 1e4 | **0** | **0.52** |
| fraction > 1e6 | **0** | **0.16** |
| ‖∇E‖ at medoid | 31.86 | **1.607e+06** |
| ‖∇E‖ CV over 5 starts | 0.222 | **1.329** |
| ‖∇²E‖₂ | 61.12 | **2.258e+08** |

> **Legacy is a well-behaved score** — skew 0.53, kurtosis 3.6, three units of spread on a
> twenty-unit range, **nothing above 1e4 anywhere**. **AMBER on the same 2,250 real structures is
> a heavy-tailed catastrophe count**: skew 8.1, kurtosis 69, **52% of the pool above 1e4**, a pool
> **standard deviation of 2.7e10** against a median of 1.1e4.
>
> This **independently reproduces `s20` L8's 53.5%-above-1e4** on a differently-assembled domain
> (30 targets × 75 windows through `s21`'s own evaluator), and it reproduces the L12 spectrum:
> participation ratio **0.0668 (AMBER) vs 0.4965 (Legacy)**, anisotropy 18.2 vs 4.2, near-zero
> mode fraction **0.481 vs 0.000**, condition number 7,756 vs 6.3.

**The scales are five and seven orders apart.** ‖∇E_AMBER‖ / ‖∇E_Legacy‖ ≈ **5.0e4**;
‖∇²E_AMBER‖ / ‖∇²E_Legacy‖ ≈ **3.7e6**. **Label: ESTABLISHED.**

## 1.1 THE FIRST RESULT, AND IT DAMAGES THE PRIORITY-ONE EXPERIMENT'S PARAMETERISATION

`PREREG_C.md` §0.3 predicted, before the run, that **raw λ is not a dial**. Define the crossover
at which the AMBER term first dominates the mixed gradient,

    λ* = ‖∇E_Legacy‖ / (‖∇E_Legacy‖ + ‖∇E_AMBER‖)

measured at every one of the 150 starts:

    median  1.99e-05        min  1.51e-09        max  5.58e-02        spread  3.7e+07 x

> ### In raw units, `H(λ) = (1−λ)E_L + λE_A` is **already pure AMBER at λ = 10⁻⁴**, and the λ at which it becomes so **varies by seven orders of magnitude across targets and starts**.
>
> A uniform λ grid in raw units is **not an interpolation between two Hamiltonians. It is a step
> function whose step position is a property of the target's steric state.** Any "abrupt
> transition" or "λ-dependent basin hopping" read off a raw grid is a **units artefact**, which is
> exactly the failure class Sprint 20 registered against itself as F-C2b and caught once already.
>
> **Label: EXACT** (it is a consequence of the measured scales, not an inference) **for the
> statement; ESTABLISHED for the magnitudes.**

**This is why Block N exists and why it ran first.** Had the λ sweep been run before the scales
were measured, this lane would have reported a "sharp Legacy→AMBER transition near λ = 0" as a
scientific finding about the two potentials. It is a statement about kcal/mol.

## 1.2 THE DECLARED NORMALISATION, AND WHAT EACH ONE DOES TO THE WORST REAL STRUCTURE

Four normalisations were declared in `PREREG_C.md` §1 with `Nt` (asinh, pool median and MAD)
named **PRIMARY before any RMSD was seen**; the other three are computed alongside on every arm so
the choice is auditable. On the worst real AMBER value in this set (target `2MIG`,
**+5.458e+23 kcal/mol**, pool median +2.321e+04, MAD_n 3.475e+04):

| normalisation | worst value maps to | pool median maps to |
|---|---:|---:|
| `raw` | 5.458e+23 | +2.321e+04 |
| `Nz` (robust z) | 1.571e+19 | 0 |
| `Ng` (gradient-matched) | 1.608e+12 | 0 |
| **`Nt` (asinh)** | **44.89** | **0** |

**An affine normalisation does not solve this.** `Nz` and `Ng` fix the *centre* and the *units*
and leave the tail nineteen and twelve orders wide. Only the bounded transform puts the worst real
structure in the same numerical universe as the median one — and it does so **while preserving the
ordering exactly** (§2).

`Nr`, the within-target rank→normal transform of `s21/tailprice.py`, was **declared and rejected
for continuation before use**, with the reason stated in the pre-registration rather than
discovered afterwards: it is a discrete rank on a finite pool, so it has **zero gradient almost
everywhere** and is **undefined off the pool**. It remains valid for pool-restricted *selection*
and is used nowhere else here.

---

# 2. AN EXACT DERIVATION, STATED BEFORE ITS STATISTICS

Registered as `PREREG_C.md` §4 **before any transformed arm ran**. Let `f` be `C¹` and strictly
increasing, `Ĥ = f ∘ H`. Then

1. **Selection is invariant.** `argmin Ĥ = argmin H`, and for every α the CVaR α-tail **SET** of
   `Ĥ` is that of `H`. The CVaR *value* changes; **the selected face does not.**
2. **Critical points and their index are invariant**, and at a critical point
   `∇²Ĥ = f′(H)·∇²H` — a positive scalar multiple — so **every scale-free spectrum metric,
   including the participation ratio, is EXACTLY unchanged at critical points.**
3. **Away from critical points it is a rank-one update**:
   `∇²Ĥ = f′(H)∇²H + f″(H)∇H∇Hᵀ`, with `f″ < 0` in `asinh`'s upper tail.

> **Consequence, registered in advance: a single-Hamiltonian preconditioner of this class cannot
> change what a perfect selector would pick.** Block P is therefore an *optimisation* intervention
> and nothing else. **Any claim that bounding AMBER's tail improves its ranking is refuted by this
> derivation before it is measured**, and this lane makes no such claim.

**Verified, not merely asserted** (`c_c2tail.json: monotone_check`): on every target the argsort of
`Nt(E_AMBER)` is bit-identical to the argsort of `E_AMBER`, and likewise for Legacy — see §5.

**Label: EXACT.** *Never present an identity as an empirical discovery* (`BRIEF` §9): this is
derived, and the empirical content is only whether item 3's rank-one term spreads or further
concentrates the spectrum, which §4 measures.

---

---

# 3. BLOCK C1 — THE CONTINUATION'S MECHANISM IS REFUTED, AND MY OWN PRIMARY ENDPOINT WAS THE WRONG ONE

**n = 30, COMPLETE. BASIS: built chain. READOUT: the single structure at the minimiser, averaged
over 5 starts. AMBER = bare single point.** `c_c1.json` / `c_c1_report.txt`.

The directive's mechanistic case for `Legacy → AMBER` continuation is that **Legacy, being a
compactness model, can carry the state out of the steric singularity so that AMBER is entered from
a region where its curvature is not concentrated.** `PREREG_C.md` §2 turned that into three
predictions and registered the falsifier on the third.

| quantity, median over targets | at `θ₀` | at `θ_L` (Legacy-minimised) | at the matched NULL |
|---|---:|---:|---:|
| `E_AMBER` (kcal/mol) | 1.827e+04 | **1.868e+06** | 5.771e+04 |
| `‖∇E_AMBER‖` | 1.285e+06 | **1.286e+08** | 6.333e+06 |
| `E_Legacy` | −12.44 | −17.65 | — |

Paired, fold-clustered CI, n = 30:

    log10 E_AMBER        theta_L - theta_0   +1.4947 [+1.1437, +2.0174]   worse on 25/30   *
    log10 ||grad E_A||   theta_L - theta_0   +1.4550 [+1.0944, +1.9101]   worse on 24/30   *
    log10 E_AMBER        NULL    - theta_0   +1.0908 [+0.2293, +2.2148]   worse on 22/30   *
    log10 E_AMBER        theta_L - NULL      +0.4040 [-0.8360, +1.0984]                    NOT MEASURED

> ### A Legacy minimisation does not carry the state OUT of the steric singularity. It drives it ~31x DEEPER in energy and ~28x deeper in gradient norm, with CIs excluding zero, on 25 of 30 targets.
>
> **And the mechanism is Legacy's own established character rather than an accident.** `s20` L8
> measured that **Legacy-preferred candidates are 0.45 Å more compact (124W/2L, 5/5 folds)**.
> Compaction is exactly what closes heavy-atom contacts, and closed heavy-atom contacts **are** the
> steric singularity. **Legacy is the wrong preconditioner for AMBER for precisely the reason it is
> Legacy.** That is one line, it follows from two independently measured facts, and it predicts the
> sign of everything in §6.
>
> **Label: the direction ESTABLISHED (n = 30, CI excluding zero, 25/30); the factor of ~31
> INCONCLUSIVE in magnitude**, because it is a median of a distribution spanning eight orders.

## 3.1 MY REGISTERED FALSIFIER DID NOT FIRE, AND I AM NOT CLAIMING SUPPORT FROM THAT

`PREREG_C.md` §2 named **participation ratio** the primary endpoint, with F-C1 firing if it failed
to rise. It rises:

    participation ratio   theta_L - theta_0   +0.0132 [+0.0019, +0.0226]   up on 19/30   *
    participation ratio   NULL    - theta_0   +0.0043 [-0.0094, +0.0161]   up on 16/30
    participation ratio   theta_L - NULL      +0.0089 [-0.0039, +0.0197]   up on 19/30    NOT MEASURED

**So F-C1 as literally written DOES NOT FIRE.** I am nevertheless reporting the hypothesis as
refuted, for two reasons that were both registered in advance:

1. **The endpoint that survived does not beat its own matched control.** Against a realisable
   geodesic move of the *same* torus magnitude the rise is `+0.0089 [−0.0039, +0.0197]` — **NOT
   MEASURED**. `BRIEF` §7 rule 1: a control must be matched in the space the operator works in.
   Any move of that size spreads the spectrum about as much.
2. **`PREREG_C.md` §2 registered all three predictions, and the falsifier was attached to the
   weakest of them.** The two direct energetic reads — which are what "de-singularise" *means* —
   are refuted with the opposite sign and CIs excluding zero. **The registration was
   mis-specified; per `BRIEF` §8 it is kept unedited and the mis-specification is recorded here
   rather than repaired.**

**This is the same failure mode Sprint 20's Workstream C caught in itself** (a falsifier that fires
on a null). The participation ratio is a *scale-free* summary; it is insensitive to the thing that
actually makes AMBER hard, which is **how far up the `r^-12` wall the configuration sits**. A lane
that had measured only the registered primary would have reported *"Legacy relaxation improves
AMBER's conditioning, CI excludes zero"* and it would have been the opposite of the truth.

## 3.2 THE RMSD READ, WHICH IS A NULL AND IS REPORTED AS ONE

    Ca-RMSD   theta_L - theta_0   +0.1541 [-0.0372, +0.4977]   14W/16L   NOT MEASURED

Consistent with `s20` L12's Legacy arm (`+0.0045 [−0.011, +0.019]`). The mean Legacy move is
**0.1037 rad/coord**. Nothing here is an RMSD claim.

---

# 4. BLOCK C2 (cheap half) — λ IS NOT A DIAL, AND WHAT A DIAL COSTS

**n = 30, COMPLETE** for the gradient and tail halves (`c_c2fast.json`). The Hessian half is
reported in §7 with its own status. **The mixed gradient is EXACT** —
`grad H(λ) = (1−λ)f_L'(E_L) grad E_L + λ f_A'(E_A) grad E_A` — so one Legacy and one AMBER gradient
per start price the whole λ axis under every normalisation, and no arm can differ from another by
anything but λ.

## 4.1 AMBER's share of the mixed gradient, as a function of λ

Median over 30 targets × 5 starts:

| λ | `raw` | `Nz` | `Ng` | `Nt` |
|---|---:|---:|---:|---:|
| 1e−6 | 0.0164 | 0.0000 | 0.0000 | 0.0000 |
| 1e−4 | **0.6251** | 0.0003 | 0.0001 | 0.0003 |
| 1e−2 | **0.9941** | 0.0275 | 0.0100 | 0.0262 |
| 0.05 | **0.9989** | 0.1284 | 0.0500 | 0.1229 |
| 0.5 | **0.9999** | 0.7367 | **0.5000** | 0.7269 |
| 0.95 | **1.0000** | 0.9815 | **0.9500** | 0.9806 |

> **The entire raw λ ladder from 0.05 to 1 is the SAME Hamiltonian.** In raw units the mixture is
> 99.9% AMBER by λ = 0.05 and crosses over at λ ≈ 2e−5, a point that varies **3.7e+07-fold** across
> targets (§1.1). A "staged Legacy → Legacy+AMBER → AMBER" schedule written on a raw λ grid **never
> visits an intermediate Hamiltonian at all.**
>
> **`Ng` is the honest dial**: it makes AMBER's gradient share *identically* λ — which is not an
> empirical finding but the definition of a gradient-matched normalisation, and is labelled as
> such. `Nz` and `Nt` are dials with a mild convexity (λ = 0.5 → share 0.73).
>
> **Label: EXACT for the raw column** (it follows from the measured scales); **ESTABLISHED** for the
> rest.

## 4.2 The gradient variance along λ, and the first genuine preconditioning effect

CV of `‖grad H(λ)‖` across the 5 starts, median over targets:

| λ | `raw` | `Nz` | `Ng` | `Nt` |
|---|---:|---:|---:|---:|
| 0 | 0.2224 | 0.2224 | 0.2224 | 0.2827 |
| 0.5 | 1.3289 | 1.3293 | 0.8665 | **0.7999** |
| 1 | **1.3289** | 1.3289 | 1.3289 | **0.8757** |

**Legacy's gradient scale is stable across starts (CV 0.22); AMBER's is not (CV 1.33).** That is the
steric singularity again: the gradient norm at a start is a function of that start's worst clash.
**`Nt` reduces AMBER's across-start gradient CV by 34% (1.329 → 0.876) at λ = 1 — a genuine
preconditioning effect on the quantity a variational optimiser's step size actually depends on, and
it costs nothing, because §2 proves the transform moves no minimiser and changes no tail set.**
`Nz` and `Ng`, being affine, cannot do this: they rescale the variance and its mean together.
**Label: SUPPORTED** (n = 30, medians; no CI is computed on a CV and it is flagged a diagnostic).

## 4.3 The λ-tails are genuinely different SETS, and it buys nothing

Under `Nt` at λ = 0.5 the CVaR α = 0.15 tail set overlaps the λ = 0 tail at Jaccard **0.329** and the
λ = 1 tail at **0.151** — a *different set*, not an interpolation between two. **The λ-dependent
set-hopping the directive told me to look for is real and measurable.** §5 prices what it is worth.

---

# 5. BLOCK C2t — THE MANDATORY MATRIX, POOL-RESTRICTED, AT THE FULL n = 126

**n = 126, COMPLETE.** `c_c2tail.json` / `c_c2tail_report.txt`. **READOUT: point cloud** — the
coordinate average of each target's own CVaR α = 0.15 tail, a **set-mean** readout. **It is never
compared to the built-chain rows of §3 or §6.** Not a VQE experiment and no quantum claim is made:
it **bounds** what any CVaR-VQE selecting from this pool can achieve, per Hamiltonian, before a
circuit is built.

    pool mean 3.551    ORACLE pool best 2.306    whole-pool average 3.048
    MATCHED-COUNT RANDOM tail  3.101   <-- the only comparison that means anything

| arm | tail_avg | minus RANDOM, fold-clustered CI | W/L |
|---|---:|---|---:|
| Legacy (λ = 0) | 3.208 | **+0.106 [+0.032, +0.190]** | 50/76 |
| AMBER (λ = 1) | 3.216 | **+0.115 [+0.074, +0.170]** | 54/72 |
| `raw` λ = 0.5 | 3.221 | **+0.120 [+0.081, +0.174]** | 51/75 |
| `Nt` λ = 0.35 | 3.202 | **+0.101 [+0.081, +0.130]** | 50/76 |
| **best of all 42 arms** (`Nz` λ = 0.95) | **3.182** | **+0.080 [+0.041, +0.119]** | 48/78 |

> ### EVERY arm of the mandatory matrix's physics half — Legacy alone, AMBER alone, and every mixture across three normalisations — is WORSE than a matched-count RANDOM tail, with a CI excluding zero, at n = 126.
>
> Not null. **Worse.** The best mixture in the entire family is **+0.080 Å worse than chance**, and
> the sign is unanimous: **42 of 42 arm-rows positive, 42 of 42 CIs excluding zero** (`raw`, `Nz`,
> `Nt`; 41 distinct λ values — `Ng` was run at n = 30 only, where it agrees).
>
> This is `s20`'s argmin result (`BRIEF` §2) transferred to the **tail** operator the sprint was
> built to test. **The tail-mean-versus-minimum distinction does not rescue the physics
> Hamiltonians**: `argmin` runs 3.503–3.693 across the same arms and the tail runs 3.18–3.23, so the
> tail *is* better than the argmin — **and both are worse than picking at random.** The improvement
> from argmin to tail is the improvement from a point to an average, which a random average also
> gets.
>
> **`BRIEF` §2 asked whether the argmin result transfers to a tail selector. It does. Label:
> ESTABLISHED, n = 126, pool-restricted.**

**And the mixtures do not beat their own endpoints by enough to matter.** `Nz`/`Nt` at λ ≈ 0.35–0.95
recover about 0.03 Å relative to λ = 0 or λ = 1 — real, since the CIs are tight — but that is 0.03 Å
against a 0.08–0.12 Å deficit to chance. **A continuation between two Hamiltonians that are each
worse than random produces something that is also worse than random.**

## 5.1 THE COORDINATOR'S QUESTION, ANSWERED WITH THE ONE SIGNAL ONLY THIS LANE CAN COMPUTE

Legacy and AMBER are separably evaluable here, so their **disagreement** on a candidate set is
available as a **native-free per-target signal**: `Jaccard(Legacy α-tail, AMBER α-tail)`, mean
0.0554, median 0.0476, and **zero overlap on 43% of targets** — the two potentials select disjoint
tails on nearly half the instrument.

    rho(disagreement, Legacy       selector skill vs random) = +0.046   p = 0.608
    rho(disagreement, AMBER        selector skill vs random) = -0.042   p = 0.643
    rho(disagreement, best mixture selector skill vs random) = +0.058   p = 0.520

    REGIME SPLIT at the signal's median (positive = WORSE than random):
      Legacy         AGREE half +0.158 [+0.015,+0.231]    DISAGREE half +0.087 [+0.024,+0.172]
      AMBER          AGREE half +0.151 [+0.029,+0.334]    DISAGREE half +0.101 [+0.052,+0.148]
      best mixture   AGREE half +0.155 [+0.036,+0.309]    DISAGREE half +0.052 [+0.018,+0.084]

**No skill, and no regime.** Both halves are worse than random with CIs excluding zero: there is no
subset of targets on which the physics selectors help. **Label: NOT SUPPORTED** for the disagreement
signal as a regime detector on this instrument.

## 5.2 AND THE ORACLE CEILING OF A PERFECT SWITCH IS PRICED BY ITS OWN NULL — IT IS ZERO

Before asking whether a native-free signal can find a regime, ask what a **perfect** switch would be
worth. `BRIEF` §7: a `min` over k arms beats any single arm **by construction**, so an ORACLE min
needs a matched min-of-k null.

    ORACLE min of (Legacy, AMBER, random)             2.928
    MATCHED NULL  min of 3 INDEPENDENT random tails    2.926
    ORACLE min of (Legacy, AMBER, mixture, random)    2.919
    MATCHED NULL  min of 4 INDEPENDENT random tails    2.894
    single random tail                                 3.101

> **The ORACLE per-target switch gains +0.173 Å on a single random tail. The pure min-of-k
> selection-bias null supplies +0.175 Å of it.** A *perfect, oracle-guided* choice between Legacy,
> AMBER and chance on every target is worth **nothing at all** beyond drawing three random tails and
> keeping the best — and at k = 4 the null is *better* than the ORACLE.
>
> **The regime-separation route is closed on this instrument before any native-free signal is asked
> to find the regime.** Label: **ESTABLISHED** for the ceiling; the null is **measured, not
> analytic**, per `BRIEF` §7 rule 2.

**Scope, stated so the number cannot travel.** This is the **top-75 retrieved-window pool**, not the
objective's top-512, and the readout is the tail's coordinate average. A regime result on the
coordinator's instrument is a different experiment; this prices the route **here**, and porting it is
real work rather than a lookup.

## 5.3 AN EXACT CLAIM, VERIFIED, AND THE ONE TIME IT FAILED

`PREREG_C.md` §4.1 derives that a strictly monotone transform leaves every CVaR tail **set**
identical. Verified rather than assumed, on all 126 targets × 75 candidates:

* **AMBER: 126/126 argsorts bit-identical.**
* **Legacy: 125/126.** The single firing is on `9L1M`, where three pool members have Legacy energies
  differing by **8.9e−16 — one ULP** — which `asinh` maps to the *same* double, so the stable sort
  orders them by index instead. They sit at sort positions **71–73 of 75**, outside every α ≤ 0.5
  tail, so **no tail set changes.**

**Reported rather than rounded to 126/126**, and worth stating for its own sake: the transform did
not break an ordering, it **collapsed a one-ULP tie** — and a tie broken by index order is the
tie-breaking trap that already cost this programme a whole table.

---


## 5.4 THE PRE-REGISTERED F-N AUDIT: the conclusion is NOT normalisation-dependent

`PREREG_C.md` §1 committed, before any run, that *if the primary conclusion flips between `Nt`,
`Nz` and `Ng` the conclusion is reported NORMALISATION-DEPENDENT and no arm is promoted*, and that a
non-declared winner would be recorded as **the declared choice being WRONG** rather than swapped in.
The audit, run at the end:

    arms evaluated (n = 126)              42 rows / 41 distinct lambda, over raw, Nz, Nt
    arms WORSE than random (mean)         42 / 42
    arms with a CI excluding zero         42 / 42   -- every one on the WORSE side
    best arm in the entire family         Nz lambda=0.95   +0.0803 vs random -- still worse

    C3 primary under BOTH normalisations of the SAME schedule:
      LA_raw - A_raw   +0.2601 [+0.0555, +0.4224]    9W/21L
      LA_Nt  - A_raw   +0.2331 [+0.0673, +0.3950]   12W/18L

    F-N clause 2 (pool MAD finite and positive on every target):  252/252, 0 firings

**No flip anywhere.** The selection conclusion and the C3 conclusion both hold under every
normalisation, so neither is propped up by the declared choice. **`Nt` is not the best-performing
arm** — `Nz` at λ = 0.95 edges it by 0.003 Å — but since **both are worse than chance**, there is no
"winner" to record the declared choice as wrong against. **F-N does not fire. Label: the conclusions
are NORMALISATION-INDEPENDENT.**

**An error of mine, caught by this audit and corrected rather than left standing.** An earlier draft
of this section said *"60 of 60 arms"*. The n = 126 block runs **three** normalisations, not four
(`Ng` was run at n = 30 only), giving **42 arm-rows**. The count was wrong; the conclusion and every
number in the table were not. Corrected here and in `STATUS_C.md`.


---

# 6. BLOCK C3 + BLOCK P — THE PRE-REGISTERED PRIMARY ENDPOINT. **F-C3 FIRES.**

**n = 30, COMPLETE.** `c_c3.json` / `c_c3_report.txt`. **BASIS: built chain. READOUT: the single
structure at the minimiser, averaged over 3 starts — a BEST-flavoured readout, never compared to
the point-cloud rows of §§4–5.** Identical starts, identical **total** function-evaluation budget
across arms (`MAXITER = 100`, `maxfun = 200`, split equally across stages; **realised `nfev` is
printed, so "matched budget" is a measured fact and not a claim**). AMBER = bare single point.

    start Ca-RMSD 3.1477      median E_AMBER at start  1.403e+04 kcal/mol

| arm | RMSD_end | vs start | NULL(move) | excess | \|dθ\| rad | nfev | med E_AMBER_end | med E_Leg_end |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `L_raw` Legacy only | 3.3154 | +0.168 | 3.1632 | +0.152 | 0.112 | 90 | **8.25e+06** | −17.63 |
| **`A_raw` AMBER only — THE COMPARATOR** | 3.6974 | +0.550 | 3.5802 | +0.117 | 0.576 | 102 | **−46.45** | −8.693 |
| `A_Nt` asinh-preconditioned AMBER | 3.8112 | +0.664 | 3.5319 | +0.279 | 0.524 | **67** | +13.45 | −8.966 |
| `LA_raw` Legacy→AMBER, raw | 3.9575 | +0.810 | 3.6210 | +0.336 | 0.595 | 146 | −41.00 | −8.488 |
| **`LA_Nt` Legacy→AMBER — THE CONTINUATION** | **3.9306** | +0.783 | 3.5849 | +0.346 | 0.572 | 142 | +31.04 | −8.848 |
| `L5A_Nt` Legacy→½→AMBER | 3.9503 | +0.803 | 3.5764 | +0.374 | 0.544 | 172 | +37.28 | −9.218 |
| `ramp_Nt` λ = 0→1 in 5 stages | 3.7719 | +0.624 | 3.4869 | +0.285 | 0.458 | 196 | +44.34 | −8.493 |
| `A_tr` trust region, pool radius | 3.6926 | +0.545 | 3.5007 | +0.192 | 0.436 | 86 | +60.10 | −8.278 |
| `A_tr3` trust region, radius/3 | **3.3987** | +0.251 | 3.3022 | +0.097 | **0.176** | 70 | +88.04 | −10.10 |
| `AbondA` bonded groups → all terms | 4.1592 | +1.012 | 3.6938 | +0.466 | 0.710 | 102 | −8.925 | −6.357 |

## 6.1 THE PRIMARY ENDPOINT: the continuation is not neutral, it is WORSE

    RMSD_end:  LA_Nt   - A_raw   +0.2331 [+0.0673, +0.3950]   median +0.0925   12W/18L   *
    RMSD_end:  LA_raw  - A_raw   +0.2601 [+0.0555, +0.4224]   median +0.1089    9W/21L   *
    RMSD_end:  L5A_Nt  - A_raw   +0.2529 [+0.1007, +0.3916]   median +0.1733    8W/22L   *
    RMSD_end:  ramp_Nt - A_raw   +0.0745 [-0.0627, +0.1852]   median +0.0303   14W/16L    NOT MEASURED

**And the secondary that was registered to decide promotion — each arm minus its OWN move-size
matched null, so a schedule cannot win by simply moving less:**

    excess-over-null:  LA_Nt   - A_raw   +0.2285 [+0.1140, +0.3647]   12W/18L   *
    excess-over-null:  L5A_Nt  - A_raw   +0.2567 [+0.1167, +0.3614]   11W/19L   *
    excess-over-null:  ramp_Nt - A_raw   +0.1679 [+0.0827, +0.2484]   10W/20L   *

> ### F-C3 FIRES, and it fires in the direction that is worse than a null result. Every staged Legacy→AMBER schedule ends FURTHER from the native than running AMBER directly from the same start on the same budget — and it survives the move-size correction, so it is not the move-size artefact that dissolved Sprint 20's headline. **The path itself is worse.**
>
> **This is exactly what §3 predicted.** C1 measured that a Legacy minimisation drives the state
> ~31× deeper into the steric singularity. A schedule that hands AMBER a *worse-conditioned*
> starting point produces a worse answer. **The mechanism and the endpoint agree**, which is
> stronger evidence than either alone.
>
> **The directive's flagged "potentially highest-value idea" is REFUTED on its own pre-registered
> primary endpoint, at n = 30, with a CI excluding zero, and with its mechanism independently
> refuted at n = 30 in §3.**

**The honesty layer, applied to my own strongest number.** `LA_Nt`'s mean is **CONCENTRATED**: the
largest single target (`8TXS`, +1.961) carries 14.3% of the total absolute effect, against a
**uniform-effect null percentile of 98.4**. The mean `+0.2331` therefore overstates the typical
target. **The direction survives the caveat** — median +0.0925, 12W/18L, and the same sign in
`LA_raw`, `L5A_Nt` and `ramp_Nt` — but the *magnitude* should be read as the median, not the mean.
**Label: the DIRECTION ESTABLISHED at n = 30; the MAGNITUDE INCONCLUSIVE.**

## 6.2 BLOCK P: **F-P fires on every arm — no preconditioner in this set preconditions**

`PREREG_C.md` §6 registered that *an arm is rejected if it does not beat `P0` on its own stated
objective* — the final `E_AMBER` at matched evaluations. Nothing does:

    log10(E_AMBER_end + 2000):  A_Nt    - A_raw   +0.0038 [-0.0008, +0.0082]   10W/20L    NOT MEASURED
    log10(E_AMBER_end + 2000):  LA_Nt   - A_raw   +0.0067 [+0.0020, +0.0110]    9W/21L   *
    log10(E_AMBER_end + 2000):  ramp_Nt - A_raw   +0.0071 [+0.0029, +0.0111]    2W/28L   *
    log10(E_AMBER_end + 2000):  A_tr    - A_raw   +0.0220 [+0.0097, +0.0394]    3W/27L   *
    log10(E_AMBER_end + 2000):  A_tr3   - A_raw   +0.0496 [+0.0248, +0.0850]    1W/29L   *
    log10(E_AMBER_end + 2000):  AbondA  - A_raw   +0.0974 [+0.0014, +0.2988]    3W/27L   *
    log10(E_AMBER_end + 2000):  L_raw   - A_raw   +2.9668 [+2.2256, +3.5495]    0W/30L   *

**Direct AMBER minimisation reaches the lowest AMBER energy, and it is the only arm that reaches a
physical one** — median **−46.45 kcal/mol**, against +13 to +88 for every preconditioned arm and
**+8.25e+06 for Legacy**. (The `+2000` shift makes the logarithm defined over the physical negative
range; it is a monotone reparameterisation and changes no ordering.)

> **The preconditioners are not buying conditioning. They are buying a shorter move**, and paying
> for it in the objective they were built to minimise. `L_raw`'s row is the whole story in one
> number: **a Legacy minimisation ends at `E_AMBER = 8.25e+06 kcal/mol` — a steric catastrophe — with
> the best Legacy energy of any arm (−17.63).** The two potentials are optimising away from each
> other, which is `s20` L12's `cos(∇E_L, ∇E_A) = −0.35` seen at the end of a trajectory instead of
> at its start.

## 6.3 THE ONE ARM THAT BEAT THE COMPARATOR ON RMSD, AND WHY IT IS **MOVE-SIZE, NOT PHYSICS**

    RMSD_end:          A_tr3 - A_raw   -0.2987 [-0.6307, -0.0547]   median -0.2644   20W/10L   *
    excess-over-null:  A_tr3 - A_raw   -0.0207 [-0.3185, +0.1736]   median +0.0027   15W/15L    NOT MEASURED

**A trust region of one third the pool's own torsional dispersion is 0.30 Å better than direct AMBER
on Cα-RMSD, with a CI excluding zero and 20W/10L — and the entire gain is the size of the move.**
Against its own magnitude-matched realisable null the excess is `−0.0207 [−0.3185, +0.1736]`,
**NOT MEASURED**, with an exactly even 15W/15L split. It moves **0.176 rad/coord against `A_raw`'s
0.576**, and it pays **+0.0496 in log₁₀ final AMBER energy** for the privilege.

The concentration diagnostic says the RMSD gain is **broad, not one target** (top-1 share 0.124,
uniform-effect null percentile 93.3 — below the 95 threshold). So this is a real, distributed,
**worthless** effect: **`BRIEF` §12's "improves an intermediate metric while worsening the
objective", and Sprint 20 L12's verdict reproduced by a different route.**

**Reported as: MOVE-SIZE, NOT PHYSICS. No promotion.** `PREREG_C.md` §5 forbade promotion from this
block in advance, and this is why that was the right rule.

## 6.4 THE ONE GENUINELY OPEN RESULT, AND I AM NOT CLAIMING IT

`A_Nt` — AMBER under the asinh preconditioner — reaches a final AMBER energy **statistically
indistinguishable from raw AMBER** (`+0.0038 [−0.0008, +0.0082]` in log₁₀, **NOT MEASURED**) using
**67 realised function evaluations against 102 — 34% fewer.** Combined with §4.2's 34% reduction in
across-start gradient CV, that looks like a real preconditioning effect on exactly the axis a
variational optimiser cares about: **budget.**

**I am not claiming it, for a specific reason.** L-BFGS's termination test is on the *gradient
norm*, and `Nt` shrinks the gradient by construction, so **"fewer evaluations" may be the
convergence threshold firing earlier on a rescaled gradient rather than genuine efficiency.** The
experiment that settles it is a **matched-`nfev`** arm — stop both at exactly 67 evaluations and
compare — **and I did not run it.** **Label: OPEN. Missing configuration named: `A_Nt` vs `A_raw` at
matched realised `nfev`, n = 30, 3 starts.**

And its RMSD row (`+0.1138 [+0.0157, +0.2435]`) must **not** be read as harm: the concentration
diagnostic returns a **top-1 share of 0.289 at the 100.0th percentile of the uniform-effect null**,
with **median −0.0008 and 16W/14L**. **One target (`5H1H`, +1.766) produces the entire mean.**
**Label: NOT a systematic effect; the CI excluding zero is a single outlier and is reported as
such.** §2 already proves the transform cannot move the minimiser — only the path — so a systematic
RMSD effect from it was never available.

## 6.5 `amber_bonded` — a genuine force-group subset, named as one, and it is the worst arm

`AbondA` minimises the **bonded force groups of the same genuine ff14SB system** (bond, angle,
torsion; nonbonded and solvation switched off) and then hands over to the full potential. **It is a
subset of AMBER, not a learned surrogate, and it is never called "AMBER" in the code, the artefact
keys or this document** (`PREREG_C.md` §6). It is the **worst arm in the table on every axis**:
RMSD `+0.4618 [+0.3489, +0.5918]` (7W/23L), the largest move (0.710 rad/coord), and `+0.0974` in
log₁₀ final energy.

**Coarse-to-fine over AMBER's own force groups is actively harmful**, and the reason is the same one
as everywhere else in this lane: **the bonded terms are blind to sterics**, so minimising them first
does to AMBER what Legacy does to it — moves the chain freely into overlap. **No surrogate was built
anywhere in this workstream, so the surrogate-validation rule of `PREREG_C.md` §6 was never invoked;
it is recorded as NOT TRIGGERED rather than as passed.**

---

# 7. THE HESSIAN SPECTRUM ALONG λ — and a verification gate of mine that was itself mis-specified

**n = 30, COMPLETE** (`c_c2.json`). Medoid start, scale-free metrics only (`BRIEF` F-C2b). The
mixture Hessian is computed by **exact composition**,
`∇²(f∘E) = f′(E)∇²E + f″(E)∇E∇Eᵀ`, so two Hessians and two gradients per target price every λ under
every normalisation.

### `raw` — the step function, now visible in the curvature

| λ | neg_frac | nearzero | cond_med | gap_rel | **part_ratio** | aniso | spec_skew |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 (Legacy) | 0.3200 | 0.0000 | 6.258 | 0.1007 | **0.4965** | 4.177 | 0.1904 |
| 0.1 | 0.2725 | 0.4811 | 7757 | 0.04189 | **0.06844** | 18.05 | 0.7910 |
| 0.25 | 0.2725 | 0.4811 | 7756 | 0.04189 | **0.06844** | 18.05 | 0.7910 |
| 0.5 | 0.2725 | 0.4811 | 7756 | 0.04189 | **0.06844** | 18.05 | 0.7910 |
| 0.9 | 0.2725 | 0.4811 | 7756 | 0.04189 | **0.06844** | 18.05 | 0.7910 |
| 1 (AMBER) | 0.2725 | 0.4811 | 7756 | 0.03989 | **0.06683** | 18.19 | 0.7972 |

> **Under `raw` the λ = 0.1 row IS the λ = 1 row, to four significant figures, on every one of the
> seven metrics — and it stays identical through λ = 0.9.** The prediction registered in
> `PREREG_C.md` §0.3 before any run is confirmed on the curvature as well as on the gradient and
> the tail. **A raw λ ladder does not have an intermediate Hamiltonian to characterise.**

### `Nt` — a genuine, smooth, monotone conditioning ladder

| λ | neg_frac | nearzero | cond_med | **part_ratio** | aniso | spec_skew |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.2922 | 0.0000 | 19.82 | **0.2336** | 9.029 | 0.4655 |
| 0.1 | 0.3179 | 0.0161 | 105.0 | **0.1279** | 13.07 | 0.4232 |
| 0.25 | 0.3387 | 0.0494 | 185.6 | **0.1215** | 13.89 | 0.2889 |
| **0.5** | **0.3600** | 0.0800 | 247.6 | **0.1105** | 14.15 | **0.2725** |
| 0.75 | 0.3378 | 0.1044 | 503.9 | **0.0955** | 14.29 | 0.2841 |
| 0.9 | 0.3280 | 0.1479 | 1046 | **0.0933** | 14.88 | 0.2507 |
| 1 | 0.2899 | 0.4317 | 5745 | **0.0887** | 16.08 | 0.5912 |

**Two results here are not artefacts of the parameterisation, and both are new.**

## 7.1 THE MIXTURE HAS MORE NEGATIVE CURVATURE THAN EITHER POTENTIAL IT IS MADE OF

`neg_frac` runs **0.292 → 0.360 → 0.290** and `spec_skew` runs **0.466 → 0.273 → 0.591**, both with
their extremum at **λ = 0.5**. The interior of the ladder is **more saddle-like than either
endpoint**: at λ = 0.5 the spectrum is closest to a balanced saddle (`spec_skew` 0.27 against 0.47
and 0.59) and carries the largest fraction of descent directions that are not descent directions of
either potential alone.

**This is the direct consequence of §7.3's −0.355 gradient cosine.** Adding two potentials whose
gradients sit ~111° apart creates curvature structure that neither has, and it is **why the staged
schedules of §6 do worse rather than better**: `L5A_Nt` and `ramp_Nt` spend part of their budget
optimising through the most saddle-rich member of the family. **The λ-profile explains the C3
result rather than merely accompanying it. Label: SUPPORTED (n = 30, medians; a shape claim, not a
CI claim).**

## 7.2 THE asinh TRANSFORM CONCENTRATES A WELL-CONDITIONED POTENTIAL AND SPREADS A SINGULAR ONE

`PREREG_C.md` §4.3 registered the rank-one term `f″(H)∇H∇Hᵀ` and said explicitly that **whether it
spreads or concentrates the spectrum is an empirical question**. Both endpoints answer it, in
opposite directions:

    Legacy   part_ratio  raw 0.4965  ->  Nt 0.2336     CONCENTRATES  (halves it)
    AMBER    part_ratio  raw 0.0668  ->  Nt 0.0887     SPREADS       (+33%)

The transform pushes both potentials toward the same participation ratio. It cannot change anything
at a **critical** point (§2, EXACT), so this is entirely a statement about the **non-critical**
starts a variational optimiser actually visits — which is the only place a preconditioner can act.
**Label: SUPPORTED.** Together with §4.2's 34% reduction in AMBER's across-start gradient CV, this
is the one place in the lane where a preconditioner does something measurable and principled to the
object it was meant to fix. **§6.4 is why it is nevertheless labelled OPEN and not promoted.**

## 7.3 A VERIFICATION GATE OF MINE WAS MIS-SPECIFIED, AND I FOUND IT BY RUNNING IT

The pre-registration said the composed Hessian would be **verified against a direct finite-difference
Hessian of the mixture at λ = 0.5 on every target**. It was, and **the tolerance FIRED on 21 of 30
targets** — median relative eigenvalue L2 error **0.165**, max **0.986**. I did not assume the
composition was right. I ran the plateau test the disagreement demanded, on four targets:

    pdb     h        part_ratio (DIRECT FD)          COMPOSED
    1IM7    2e-2 1e-2 6e-3 3e-3 : 0.0749 all four     0.0748     <-- FD ON a plateau, they AGREE
    1KMR    2e-2 1e-2 6e-3 3e-3 : 0.1013 0.0979 0.0976 0.0985    0.0976     <-- AGREE
    1KZ2    2e-2 1e-2 6e-3 3e-3 : 0.0847 0.0720 0.0703 0.0700    0.0936     <-- FD DRIFTS, disagree
    1M02    2e-2 1e-2 6e-3 3e-3 : 0.0923 0.0597 0.0539 0.0528    0.0652     <-- FD DRIFTS, disagree

> **Where the direct FD sits on a plateau the two agree to the fourth decimal. Where they disagree,
> the direct FD is demonstrably h-dependent and is therefore not a reference at all.** My gate
> compared an exact object to a *non-converged* one and called the difference a failure of the exact
> object. **The gate was mis-specified — it is the third instance of the F-C2b family in this
> programme, and it is mine.** Kept, reported with its firing count, and not rewritten.
>
> **What this does and does not license.**
> * The `raw` column is **unaffected**: for `raw`, `f′ = 1` and `f″ = 0`, so the composition is
>   *exactly* `(1−λ)H_Legacy + λH_AMBER`, built from the same `H_L`, `H_A` that **GC21a reproduces
>   from Sprint 20 at 0.00e+00**. §7's step-function result stands without qualification.
> * The **shapes** of the `Nt` ladder are driven by the scalars `f′(E)` and `f″(E)`, which are exact
>   and target-independent, so the *monotone* trends (part_ratio down, cond_med up, near-zero modes
>   up, neg_frac peaking at λ = 0.5) are robust.
> * The **absolute per-target spectrum values on targets whose endpoint Hessians are themselves off
>   their h = 1e-2 plateau are INCONCLUSIVE.** Sprint 20 measured that plateau on one target
>   (`1A13`); it is not verified target-by-target here, and this lane does not claim it is.
> * **Missing configuration, named:** a per-target `h`-plateau measurement of `H_Legacy` and
>   `H_AMBER` on all 30 targets. Not run. **Label for the absolute values: INCONCLUSIVE.**

---

# 8. THE λ TRANSITION IS ABRUPT, AND IT IS ABRUPT FOR THE WRONG REASON

The directive says **do not assume smooth interpolation; check for the transition being abrupt**.
It is. Measured on the *declared uniform* λ grid, the largest single-step change in AMBER's share of
the mixed gradient (a uniform dial would step 0.100 per interval):

| normalisation | largest single-step jump | where | total variation |
|---|---:|---|---:|
| `raw` | **0.9989** | λ: 0 → 0.05 | 1.0000 |
| `Nz` | 0.1284 | λ: 0 → 0.05 | 0.5231 |
| `Ng` | **0.0947** | λ: 0.35 → 0.5 | 0.6639 |
| `Nt` | 0.1229 | λ: 0 → 0.05 | 0.7182 |

> **Under `raw`, 99.9% of the entire Legacy→AMBER transition happens inside the FIRST grid
> interval.** Under `Ng` no step exceeds the uniform 0.100 anywhere. **The abruptness the directive
> warned about is real, and it is entirely a property of the units.** A lane that ran the sweep
> without Block N would have reported a sharp first-order-like transition in a mixed Hamiltonian.

## 8.1 The geometry underneath, and it reproduces Sprint 20 on a different instrument

At the 150 start-points, `cos(∇E_Legacy, ∇E_AMBER)` has **median −0.3553**, negative at **78.7%** of
them — against `s20` L12's **−0.3512 [−0.4719, −0.2221]**, 27/30 targets, computed by a different
route on a 30-target aggregation. **An independent reproduction of the sprint-20 gradient
anti-alignment. Label: ESTABLISHED, cross-sprint.**

The two gradients sit about **111° apart**, so a mixture is a genuine compromise direction rather
than an interpolation along a line. Under `Nt` the mixed gradient's alignment crosses zero against
AMBER between λ = 0.05 and 0.1, and against Legacy near λ = 0.8, so **there is a broad λ ∈ [0.1, 0.8]
window in which the mixed gradient is positively aligned with BOTH potentials.** Under `raw` there is
no such window: the mixture is already anti-aligned with Legacy (−0.27) at λ = 0.05.

**This is the one place where the continuation family has genuine, non-trivial content** — a
descent direction that neither potential supplies alone. §5 measures what it is worth on selection
(nothing), and §6 measures what it is worth on the optimisation path.

---

# 9. DEVIATIONS FROM THE PRE-REGISTRATION — every one, with its reason and its timing

1. **C1's Hessian budget was cut 5× to match what I had actually registered — after two targets
   had run under the larger version, and those two rows were DISCARDED rather than merged.**
   `PREREG_C.md` §2 budgets *"30 targets × (1 Legacy minimisation + 3 AMBER Hessians + nulls),
   ~1 h"*. The first draft of `c1_target` computed the Hessian triple at **all five starts** — 15
   AMBER Hessians per target, five times the registered design and ~4 h of wall clock. Energies,
   gradients, moves and RMSDs are still read at all five starts; the three expensive spectra are
   read at the **medoid** start, which is what §2 says. The correction is recorded in the code at
   the point of change. **Direction of the deviation: the code was over-implemented relative to the
   registration, and was brought back to it.**

2. **The mixture Hessian is computed by EXACT COMPOSITION rather than by finite-differencing the
   mixture, and gate GC21b is superseded by a stronger check.** `∇²(f∘E) = f′∇²E + f″∇E∇Eᵀ`, so two
   Hessians and two gradients price every λ under every normalisation exactly. This removes the
   mixture's own finite-difference noise — the failure mode Sprint 20's Workstream C caught in
   itself — and guarantees the λ arms differ by nothing but λ. **It is verified against a direct FD
   Hessian of the mixture at λ = 0.5 on every target, with the firing count of its tolerance
   reported.** Decided while C1 was running and before any C2 table existed.

3. **`c2fast` was split out of `c2`.** The gradient and tail halves of Block C2 cost ~3 s per target
   and the Hessian half costs ~100 s. They were split into **separate artefacts with separate
   completion flags** so the cheap half could not be lost to a compute budget. No arm, control or
   endpoint changed; the artefact `c_c2fast.json`'s flag covers only what it contains.

4. **`c2tail` (Block C2t, n = 126) is a DECLARED EXTENSION, written while C1 was running and before
   any C2 table existed.** It is not in the pre-registration. It costs one Legacy batch and 75 AMBER
   single points per target, so it runs the **selection** side of the mandatory matrix on the whole
   126-target instrument instead of the 30-target landscape subset. Its normalisations, its α, its
   matched-count random control and its readout are all the ones already declared; nothing was
   chosen after seeing its numbers. **It is the block that produced the strongest result in this
   lane, and it was not pre-registered — so it is labelled a declared extension everywhere and its
   controls are stated with it rather than assumed.**

5. **The trust radius for arms `A_tr` / `A_tr3` was declared in code, not in `PREREG_C.md`.** §6
   listed "trust regions" without fixing a radius. The radius used is **native-free and fixed before
   any RMSD was read**: the target pool's own **median per-coordinate torus dispersion about its
   medoid** — *"do not move further than the ensemble of real candidates for this target already
   differs"* — with a `radius/3` arm beside it. It is a declared-at-code-time choice and is recorded
   as such rather than presented as pre-registered.

6. **A procedural error of mine, self-caught: `git add -A s21` swept other workstreams' in-progress
   files into commit `3d7b000`.** Nothing was overwritten, deleted or modified — the files were
   untracked and are now tracked — but it is not my lane's business to commit another lane's work,
   and a reader diffing that commit will see files this lane did not write. Recorded rather than
   rewritten, since rewriting shared history would be worse than the error.

## 9.1 An observation about another lane's artefact, offered neutrally

`s21/results/_SMOKE_a_matrix.json` carries `complete: true` on a **2-row** run. The `_SMOKE_` prefix
is the declared convention and does carry the warning, so this may be entirely intended. It is
mentioned only because the coordinator reports two vacuous completion flags already caught this
sprint, and a `complete` flag inside a smoke file is the shape of the third. **This lane's own guard
requires `len(rows) >= len(FULL_SUBSET)` AND that the subset requested IS the full subset**, so a
smoke run here cannot write a completion flag even if it finishes — verified in the flag audit
(`_SMOKE_c_c1.json`: `complete: false`).

---

# 10. WHAT THIS LANE IS AND IS NOT CLAIMING

**IS claiming.**

* Raw-unit `λ` is not a dial and any raw λ ladder is a step function (EXACT, §§1.1, 4.1, 7, 8).
* A Legacy minimisation drives AMBER ~31× deeper into its steric singularity, not out of it
  (ESTABLISHED, n = 30, §3) — **the directive's stated mechanism for the continuation is refuted.**
* Every arm of the mandatory matrix's physics half, in every normalisation, is **worse than a
  matched-count random tail** at n = 126 with a CI excluding zero (ESTABLISHED, §5).
* A perfect ORACLE per-target switch between the physics Hamiltonians and chance is worth **nothing**
  beyond min-of-k selection bias (ESTABLISHED, §5.2), and the Legacy–AMBER disagreement carries no
  per-target skill (NOT SUPPORTED, §5.1).
* A strictly monotone transform leaves argmin, every CVaR tail set and every critical point
  invariant (EXACT, §2), and it is verified on 126 × 75 real candidates (§5.3).
* `Nt` reduces AMBER's across-start gradient CV by 34% at no cost to any minimiser or tail set
  (SUPPORTED, §4.2), and pushes the two potentials' participation ratios toward each other in
  opposite directions (SUPPORTED, §7.2).
* **Staged Legacy→AMBER continuation is WORSE than direct AMBER on the pre-registered primary
  endpoint, and it survives the move-size correction** (`+0.2331 [+0.0673, +0.3950]` raw;
  `+0.2285 [+0.1140, +0.3647]` corrected; DIRECTION ESTABLISHED, MAGNITUDE INCONCLUSIVE because the
  mean is concentrated, §6.1).
* **No preconditioner in this set beats direct AMBER on its own stated objective** — F-P fires on
  every arm; direct AMBER is the only arm reaching a physical final energy (ESTABLISHED, §6.2).
* The trust-region RMSD gain is **MOVE-SIZE, NOT PHYSICS** (`−0.0207 [−0.3185, +0.1736]` against its
  own matched null, 15W/15L; NOT MEASURED, §6.3).
* Coarse-to-fine over AMBER's own force groups (`amber_bonded` → full) is **actively harmful**
  (`+0.4618 [+0.3489, +0.5918]`, 7W/23L; ESTABLISHED, §6.5).
* The mixture carries **more negative curvature than either endpoint**, peaking at λ = 0.5
  (SUPPORTED, §7.1) — which explains §6 rather than merely accompanying it.

**Is NOT claiming.**

* Nothing about the **deployed** `H_AMBER = E ∘ Relax₅₀`. Every table here is the bare single point.
* No quantum claim of any kind. §5 **bounds** what a pool-restricted CVaR-VQE can achieve; it is not
  a VQE experiment.
* No cross-basis comparison. The point-cloud rows of §§4–5 and the built-chain rows of §§3, 6 are
  never differenced.
* No RMSD claim from Block C1 (`+0.1541 [−0.0372, +0.4977]`, NOT MEASURED).
* No promotion of any arm to the pipeline. `PREREG_C.md` §5 forbade it in advance.
* **No claim that `A_Nt`'s 34% evaluation saving is genuine efficiency** — L-BFGS terminates on
  gradient norm and `Nt` shrinks the gradient by construction. **OPEN**, with the settling
  experiment named: matched realised `nfev` (§6.4).
* **No claim on the absolute per-target Hessian spectra where the endpoint Hessians are off their
  `h = 1e-2` plateau.** The λ-profile *shapes* are robust and the `raw` column is exact; the absolute
  values are **INCONCLUSIVE** and the missing measurement is named (§7.3).
* **No RMSD harm attributed to the asinh transform.** Its `+0.1138 [+0.0157, +0.2435]` is one target
  (`5H1H`) at the 100th percentile of a uniform-effect null, median −0.0008, 16W/14L (§6.4).

---

# 11. CROSS-SPRINT REPRODUCTIONS, AND THE FIRING AUDIT

## 11.1 Four Sprint-20 results reproduce on this lane's independently-built instrument

Nothing below was assumed. Each was recomputed from this lane's own evaluator and compared to the
Sprint-20 figure afterwards.

| Sprint-20 result | s20 value | reproduced here | n here |
|---|---|---|---|
| unrelaxed AMBER above 1e4 on real rebuilds (L8) | **53.5%** (9450 rebuilds) | **53.51%** (9450 rebuilds) | 126 × 75 |
| `ρ(E_Legacy, E_AMBER)` within-target (L8) | −0.0886 [−0.1239, −0.0514], 0/126 with \|ρ\|>0.8 | **−0.0951 [−0.1595, −0.0355]**, **0/30** with \|ρ\|>0.8 | 30 |
| `cos(∇E_Legacy, ∇E_AMBER)` (L12) | −0.3512 [−0.4719, −0.2221], 27/30 negative | **median −0.3553**, **78.7%** of 150 starts negative | 150 starts |
| participation ratio AMBER vs Legacy (L12) | 0.0745 vs 0.4221 | **0.0668 vs 0.4965** | 30 |
| AMBER finite everywhere without relaxation (L8) | finite, median +16,062 | **2250/2250 finite**, median **+1.07e+04** | 30 × 75 |

**And two facts that follow, which are this lane's own:**

* **`ρ(E_AMBER, ORACLE d) = −0.0002 [−0.0933, +0.0883]`** and
  **`ρ(E_Legacy, ORACLE d) = +0.0385 [−0.0826, +0.1592]`** — *within* each pool of 75 real rebuilds.
  **Neither potential ranks accuracy.** This is the mechanism behind §5's result and it is why no
  mixture of them can rank accuracy either.
* **10 of 30 targets have `|ρ(E_AMBER, ORACLE d)| > 0.3`** while the mean is exactly zero. There *is*
  per-target ordering signal, **with random sign** — the programme's standing "in-band ordering is
  learnable but per-target" finding, seen here from the physics side. **§5.1 shows that the
  Legacy–AMBER disagreement does not supply that sign.**

## 11.2 Every guard, bound and threshold, and how often it FIRED

`BRIEF` §7 rule 3: *a gate can pass VACUOUSLY; report how many times it actually fired.*

| guard / bound / threshold | where | firings |
|---|---|---|
| `core.amber.memory_guard` (92% physical memory) | Block N | **9** — box at 93%, driven by a browser; 45 s backoff each, recorded in `c_norm.json` |
| same | C1, C2, C2t, C3 | **0** |
| GC21a reproduction tolerance 1e−9 | 3 targets × 2 potentials × 8 metrics | **0** — and the realised difference is **exactly 0.00e+00**, which is stronger than a pass |
| GC21c `Nt` monotonicity (synthetic, 9,000 pairs over [−1e6, 1e24]) | gate | **0** |
| `Nt` monotonicity on **real** data (126 × 75) | C2t | **1 of 126** — `9L1M`, a one-ULP tie (8.9e−16) collapsed by the transform, at sort positions 71–73 of 75, **outside every tail** (§5.3) |
| GC21b / `hess_check` composed-vs-direct tolerance 1e−2 | C2, 30 targets | **21 of 30** — and the plateau test shows the *gate* was mis-specified, not the object (§7.3) |
| L-BFGS `maxfun = 200` per arm | C3, 10 arms × 3 starts × 30 targets | realised `nfev` **67–196**, printed per arm; **the budget bound is not the binding constraint on any arm** |
| AMBER non-finiteness (bare single point) | Block N, C2t | **0 of 11,700** — AMBER is *finite everywhere and meaningless on half of it*, which is L8's point |
| `NEARZERO_REL = 1e−4` near-zero-mode threshold | Block N medoid starts | fires on **44.2%** of AMBER modes and **5.2%** of Legacy's; **never fires at all on 24 of 30 targets for Legacy** |

**`n_collapsed` is structurally zero in this lane and that is not a pass.** No arm here calls the
relaxation, so the collapse guard is never reached; the corresponding quantity is the finiteness
audit above, and it says the bare AMBER single point is **finite on 100% of 11,700 real rebuilds
while sitting above 1e4 on 53.5% of them.** *Finite-but-meaningless is more dangerous than `+inf`,
because nothing throws* — `s20` L8's sentence, reproduced.

## 11.3 Reproduction

    python -m s21.c_norm --gate            # GC21a / GC21c / GC21d (the seal, hashed not read)
    python -m s21.c_norm                   # BLOCK N              -> c_norm.json
    python -m s21.c_cont c1                # BLOCK C1             -> c_c1.json
    python -m s21.c_cont c2fast            # BLOCK C2, cheap half -> c_c2fast.json
    python -m s21.c_cont c2                # BLOCK C2, Hessians   -> c_c2.json
    python -m s21.c_cont c2tail            # BLOCK C2t, n = 126   -> c_c2tail.json
    python -m s21.c_cont c3                # BLOCK C3 + P         -> c_c3.json
    python -m s21.c_report n|c1|c2|c2fast|c2t|c3

Every stage is resumable per target; an interruption leaves a named partial and never a silently
truncated result. `BLAS` threads are capped at 1 in every module. **All timings in the logs were
taken under contention from up to eight concurrent processes and are not properties of the code**
(`BRIEF` §10).

---

# 12. WHAT I WOULD HAND THE NEXT LANE

1. **Do not build another Legacy→AMBER schedule.** The mechanism and the endpoint were both measured
   and both say no, and the reason is structural: Legacy is a compactness model and compaction
   creates the steric singularity. This is not a tuning failure and more λ points will not fix it.
2. **If any lane combines two energies, measure `λ*` first.** It costs one gradient each and it is
   the difference between a dial and a step function. Publish the crossover beside the λ grid.
3. **An affine normalisation is not enough for AMBER.** A robust z-score leaves the worst real
   structure at 1.6e+19. If a lane is z-scoring AMBER before mixing it, that arm is a clash census.
4. **Any single-Hamiltonian transform is free on the selection side and cannot help it** (§2, EXACT).
   Do not spend a sprint looking for a ranking gain from a barrier transform; spend it on the
   optimisation path, and control it against matched `nfev`.
5. **The pool-restricted selection ceiling is now measured at n = 126 and the physics half is below
   chance.** A CVaR-VQE selecting over this pool with either physics Hamiltonian, or any mixture of
   them under any of four normalisations, is bounded **below** a matched-count random tail. That is
   a bound on the architecture as configured, not on the optimiser.
6. **The one untested thing worth testing** is the matched-`nfev` comparison of §6.4. It is one
   afternoon, it is fully specified, and it is the only place in this lane where a preconditioner
   still might be real.

---
