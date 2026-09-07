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

| block | what | status | artefact |
|---|---|---|---|
| gates | GC21a/b/c/d | **PASS** | `s21/results/c_gate.json` |
| **N** | component measurement, normalisation declared | **COMPLETE 30/30** | `c_norm.json` |
| **C1** | the mechanism: does Legacy relaxation de-singularise AMBER? | *(see §3)* | `c_c1.json` |
| **C2** | the λ sweep — gradient, tail, Hessian | *(see §4)* | `c_c2.json` |
| **C2t** | the selection-side λ sweep at **n = 126** | *(see §5)* | `c_c2tail.json` |
| **C3+P** | staged schedules and preconditioners | *(see §6)* | `c_c3.json` |

*(This table is updated as each block lands; a partial is quoted at its own n or not at all.)*

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

*(§§3–7 are written as each block lands. The caveats above were written before those tables
existed and cannot be softened after the fact.)*
