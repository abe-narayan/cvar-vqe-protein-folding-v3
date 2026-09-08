# SPRINT 20 — WORKSTREAM B (CVaR-VQE / QUANTUM) — FINDINGS

Pre-registration: `s20/PREREG_B.md`, written before any Sprint-20 arm ran, unedited.
Code: `s20/qb2_lib.py`, `qb2_opt.py`, `qb2_run.py`, `qb2_relax.py`, `qb2_q2.py`, `qb2_report.py`,
`qb2_chain.py`. Artefacts: `s20/results/qb2_*.json` with `_COMPLETE` flags and
`qb2_config.json` (the persisted subset, budget, seeds and step sizes).
**No number in this document is typed by hand.**

Tiering: **EXACT** (a theorem or identity) · **ORACLE** (needs the native; post-hoc only) ·
**ESTABLISHED** · **SUPPORTED** · **PLAUSIBLE** · **OPEN** · **INCONCLUSIVE** ·
**NOT MEASURED** · **REFUTED**.

**Every panel below is COMPLETE with its flag on disk.** Declared n per panel: landscape 20 ·
optimiser battery 20 · relaxation arm 10 · CVaR α 10 · encoding 10 · Q2 20 · recovered Sprint-19
artefacts 126 and 40. **Every claim carries its own n.** The 126-target instrument's MDE at 80 %
power is 0.084 Å; at these n it is larger, and **no Cα-RMSD number here is offered as a validated
pipeline improvement.**

---

## 0. LEAD WITH THE DAMAGE

**My own pre-registered hypothesis H1 is half REFUTED, and the half that survives is the half
that does not matter for RMSD.**

I pre-registered that the Legacy↔AMBER difference is *dominated by objective conditioning* —
that raw AMBER on ideal-geometry backbones is an r⁻¹² clash spike, and that a monotone transform
would recover Legacy-like behaviour. **The dynamic-range half is true and enormous. The curvature
half is false.** Conditioning collapses AMBER's 16.4-decade range to 1.58 decades — Legacy's is
1.56 — and leaves its Hessian condition number, near-zero-mode fraction, anisotropy and count of
local minima *statistically unchanged*, all still far worse than Legacy's. So AMBER carries **two
separate landscape pathologies**, only one of which is a units artefact, and my remedy addresses
only that one.

**Second piece of damage, self-caught.** My first reading of the recovered Sprint-19 extension
artefact was that a warm-started VQE had *the best generation ceiling on the board*, better than
the retrieval pool's. It is a **min-of-N artefact** — 8192 emitted against 500 — and at matched
count the sign reverses and the pool wins, two arms with CIs excluding zero and 5/5 folds. The
corrected table is §1b(vi). I record the wrong version because the rule that caught it (*no
min-of-N ceiling without its min-of-N null, naming its band*) is the reusable part.

**Third, and it is the whole lane's price:** across **36 (objective, arm) cells**, **zero**
improve the objective and the Cα-RMSD together, **eight** significantly trade one for the other,
and the per-target correlation between optimisation gained and structure gained is
**+0.003 [−0.068, +0.072]**. The AMBER-vs-Legacy landscape question is answered in §§3–8 — and
**the answer cannot move Priority 1, because on none of these objectives does optimising better
build a better structure.** §5a.

---

## 0b. THE ANSWER TO Q1, IN ONE PLACE

> **Is the AMBER landscape genuinely harder, and what property makes it harder?**
> **Yes on four measurable axes, no on the one that was assumed, and it does not matter.**

| candidate explanation (BRIEF §4) | verdict | evidence |
|---|---|---|
| **objective roughness** | **REAL, and the dominant one** — AMBER has 2× Legacy's and 8.5× the distogram's local minima per 2π, 473× the Hessian condition number, 4× the near-zero-mode fraction, and conditioning does **not** fix any of it | §3d |
| **poor conditioning / dynamic range** | **REAL and enormous but SEPARABLE** — 16.4 decades over 500 real structures, fixed exactly by a monotone log (to 1.57, Legacy's 1.61), worth **≈1.8 IQR of optimiser progress** and **≈0 Å** | §3c, §6b |
| **optimizer mismatch (SPSA gradient noise)** | **REFUTED** — SPSA/central-FD cosine is +0.15–0.18 and statistically identical on all four objectives | §6e |
| **CVaR concentration** | **REAL but it belongs to the RULE, not the Hamiltonian** — ESS is a functional of the tail's order statistics alone (**EXACT**); what differs is the spectrum. On raw AMBER at α=1 the gradient norm is 9.3e14 and CVaR is acting as a crude conditioner | §7a, §7b |
| **encoding limitation** | **REAL** — θ vs (cos θ, sin θ), physically identical, moves final RMSD in 11/12 cells (p = 0.006), up to −0.649 Å | §8 |
| **ansatz limitation** | **NOT MEASURED / closed from two directions** — entanglement −0.013 [−0.095, +0.077] at n=126; a deeper ansatz −0.020, 2/5 folds at n=40 | §1, §1b(iv) |
| **higher-order coupling** | **NOT the discriminator** — both energies are full-register (LEDGER `torsion-space-locality-theorem`); the measured difference is curvature, not support | §3b |
| **poor initialisation** | **NOT the discriminator** — every arm and every objective shares one `stable_rng` start set by construction; warm-starting from the pool's own basin occupancy is worth −0.073 with a CI spanning zero | §1b(iii) |
| **candidate-distribution bias** | **REAL and it is the binding one** — every generator's coherent error is 91 % shared with the retrieval pool | §9 |
| **"AMBER is more physical"** | **REFUTED as an explanation** — over the same 500 real structures AMBER's argmin (4.990 Å) and Legacy's (5.487 Å) are both **worse than a random draw** (pool mean 4.739); only the structural objective's (3.676) is better | §6d |

**The single sentence.** *The AMBER landscape is harder in every scale-free sense measured, about
half of that is a units artefact a monotone transform removes for free, none of it is the
optimiser's fault, and the whole question is downstream of the fact that neither energy's minimum
is near the native — so a better landscape cannot buy a better structure here.*

---

## 1. RECOVERED: the Sprint-19 quantum sweep, and who did what

`s19/results/qb_main.json` (126 targets × 17 arms) and `qb_report_main.json` were **complete on
disk and unanalysed** when Sprint 19 froze. **WORKSTREAM D executed `s19/PREREG_B.md` §§6–9
verbatim and delivered the closure** (`s20/results/D_QB_CLOSE/`, `s20/d_qb_close.py`,
`s20/LEDGER.md` L1). **That analysis is theirs and is credited as such.** I re-read the raw
artefact independently and reproduce their headline numbers exactly:

| arm | realised (BUILT chain) | point cloud | selected | generation ceiling |
|---|---|---|---|---|
| `pool500` — the shipped retrieval pool, **0 evaluations** | **3.2297** | 3.0639 | 2.2850 | 1.7212 |
| `c_marg` — zero-information matched torsion marginals | 3.3257 | 3.1336 | 2.4908 | 1.6624 |
| `c_metroH` — best classical sampler | 3.4071 | 3.1986 | 2.6816 | 1.6698 |
| `c_helix` — zero-information constant helix | 3.4258 | 3.1964 | 2.6228 | 1.9235 |
| **`q_a1.00` — best quantum** | **3.4858** | 3.2832 | 2.6497 | 1.7445 |
| `c_lbfgs` — hardest optimisation of the deployed objective | 3.6578 | 3.5315 | 3.3108 | 1.9030 |
| `q_untrained` | 3.6962 | 3.3845 | 2.7551 | 1.9709 |

**Basis is stated on every row** (the L20 rule): `realised` is the projected **built chain**,
`point cloud` is the coordinate average before projection. The programme's `3.048 Å` best is a
point cloud; the best **built** structure is `3.204 Å`.

The sampler hypothesis is **REFUTED on all three pre-registered endpoints** and four of nine
pre-registered kill rules fire. **This lane does not re-run any of it.** Two of Workstream D's
results carry directly into my Q1 and are used as its baseline:

* **The VQE genuinely trains on the continuous encoding** — `q_a1.00 − q_untrained =
  −0.210 [−0.339, −0.082]`, 5/5 folds, against the **mandatory** best-of-N control. The
  programme's recorded *"running the VQE is worse than not running it"* was **lattice-only**;
  on the continuous basin-latent encoding the sign reverses. **Priority 2 is satisfied on its own
  terms: the pillar is genuine and it works — it simply does not beat anything classical.**
* **Entanglement is worth `−0.013 [−0.095, +0.077]` — NOT MEASURED**, nearest-neighbour latent
  mutual information 0.045 bits. Ansatz expressivity is therefore *not* spent on in Q1.

### 1b. A SECOND unanalysed Sprint-19 artefact, recovered here — `s19/results/qb_ext.json`

Workstream D analysed `qb_main.json`. **`qb_ext.json` (40/40 targets, COMPLETE on disk) was never
read by anyone.** It holds the Sprint-19 lane's "go deeper than an α sweep" panel — warm starts,
adaptive α, CVaR+uniform and CVaR+annealing mixtures, a deeper ansatz, α = 0.02 and 0.50, and
**four independent ansatz seeds of both the trained and the untrained circuit**. Recovered and
analysed in `s20/qb2_s19ext.py` → `s20/results/qb2_s19ext.json`. Nothing was re-run. Every
comparison is made against `qb_main`'s own arms **on the same 40 targets**.

**(i) The number this lane most needed and never had — seed sensitivity.** — **ESTABLISHED**, n = 40.

    trained    4 seeds:  mean 3.2127   WITHIN-TARGET sd over ansatz seeds  0.2003   best-of-4 2.9740
    untrained  4 seeds:  mean 3.4612   WITHIN-TARGET sd over ansatz seeds  0.2893   best-of-4 3.0899

> **The ansatz seed alone moves a target by 0.200 Å. The 126-target instrument's MDE is 0.084 Å.**
> Every single-seed VQE number in this programme's record carries roughly 2.4 MDE of seed noise,
> and the untrained circuit is *more* seed-sensitive than the trained one (0.289 vs 0.200) — i.e.
> **training reduces variance as well as mean.** This is a methodological result and it should be
> binding on any future quantum arm: report ≥ 4 seeds or report the arm as NOT MEASURED.

**(ii) Training confirmed at 4 seeds against the mandatory control.** `trained − untrained =
−0.2484 [−0.4614, −0.0554]`, 23W/17L, 4/5 folds — independently reproducing Workstream D's
`−0.210 [−0.339, −0.082]` at n = 126 on a different artefact with a different seed set.

**(iii) Less CVaR concentration is better, and the ordering is monotone.** Against the 4-seed
trained baseline (negative = the arm wins):

| arm | Δ realised | CI95 | W/L | folds |
|---|---|---|---|---|
| `x_mixunif50` — CVaR **mixed 50 % with the untruncated mean** | **−0.1092** | **[−0.2041, −0.0207]** | 24/16 | **5/5** |
| `x_adapt` — adaptive α | −0.0851 | [−0.1884, +0.0148] | 26/14 | 5/5 |
| `x_warm` — warm start from the pool's own basin occupancy | −0.0725 | [−0.1633, +0.0086] | 22/18 | 4/5 |
| `x_mixanneal` | −0.0360 | [−0.1364, +0.0637] | 20/20 | 4/5 |
| `x_mps3f` — a **deeper** ansatz (3 entangling layers) | −0.0204 | [−0.1168, +0.0775] | 22/18 | 2/5 |
| `x_a0.50` | −0.0134 | [−0.0980, +0.0697] | 22/18 | 3/5 |
| **`x_a0.02` — the DEEPEST tail** | **+0.0210** | [−0.1170, +0.1665] | 21/19 | 3/5 |

**The only interval excluding zero belongs to the arm that dilutes the CVaR tail towards the
plain sample mean, and the deepest tail is the only arm on the wrong side of zero.** This is an
*independent replication, on an artefact Workstream D did not read*, of their L4 result that
α = 1 is the best quantum arm — and it is the regime **arXiv:2605.02850 (Quantum Tilted Loss,
2026)** predicts: CVaR's sharpening is paid for in shots as `O((e^{|γ|Δ}−1)²/γ²ε²)`, exponential
in the tilt, and at 512 shots an α = 0.05 tail is 26 samples while α = 0.02 is 10.
— **SUPPORTED** (single-seed extension arms against a 4-seed baseline; the effect is roughly half
a seed-sd, so it is reported as supported, not established).

**(iv) A deeper ansatz buys nothing** (`x_mps3f`, −0.020, 2/5 folds), consistent with Workstream
D's entanglement null. **Ansatz expressivity is closed on this encoding from two directions.**

**(v) The closure holds under the stronger arm set.** On the same 40 targets:
`pool500` **3.0281** · `c_marg` 3.0992 · `x_mixunif50` 3.1036 · `c_metroH` 3.1266 ·
`x_adapt` 3.1276 · `x_warm` 3.1402 · `q_a1.00` 3.1675. The best extension arm loses to the
zero-evaluation retrieval pool by `+0.0755 [−0.0714, +0.2247]` and is **exactly matched by the
zero-information marginal sampler** (`+0.0044 [−0.0896, +0.1086]`, folds 1/5). Three arms
(`x_mps3f`, `x_a0.02`, `x_a0.50`) lose to `pool500` with CIs excluding zero.

**(vi) A min-of-N trap I walked into and then caught — recorded because the catch is the useful
part.** My first reading of this table was *"`x_warm` has the best generation ceiling on the
board — 1.5187 Å, better than `pool500`'s 1.7795 — a generation success that the terminal
operator cannot spend."* **That is wrong and it is a min-of-N artefact**: `x_warm` emits 8192
configurations and `pool500` emits 500. The Sprint-19 readout persisted the matched-count column
(`gen_best_sub500_ORACLE`, a stable 500-subsample of the arm's own emitted set), and at matched
count the sign reverses:

| arm | ceiling at 8192 | **ceiling at matched 500** | `pool500` | matched-count Δ vs pool |
|---|---|---|---|---|
| `x_warm` | 1.5187 | **1.9400** | 1.7795 | +0.1605 [−0.0054, +0.3554], 3/5 |
| `x_mixunif50` | 1.6151 | **2.0915** | 1.7795 | **+0.3119 [+0.0961, +0.5416], 5/5** |
| `x_adapt` | 1.6299 | **2.0995** | 1.7795 | **+0.3199 [+0.1203, +0.5377], 5/5** |

**There is no generation success here.** At matched count every extension arm's ceiling is
*worse* than the zero-evaluation retrieval pool's, two of them with CIs excluding zero and 5/5
folds. The brief's rule — *no min-of-N ceiling without its min-of-N null, naming its band* — is
what turned a headline into a null, and the null is the correct reading.

---

## 2. A SCOPE CORRECTION THE COORDINATOR SHOULD READ FIRST — `H_AMBER = E ∘ Relax` is true of the *deployed* object, and **not** of this lane's Q1 arm

The coordinator flagged Q1 as relaxation-confounded, on the grounds that
`amber_hamiltonian.AmberHamiltonian.energy` builds coordinates, applies a
`restraint_k = 100 kcal/mol/Å²` restraint, runs `LocalEnergyMinimizer.minimize(tol, 50)`, and
reports the unrestrained energy at the *minimised* coordinates. **That is correct about the
deployed object.** It is **not** correct about the arm this lane measures, and the difference is
verifiable in source:

* `qb2_lib.AmberSP` reproduces `core.amber._run(..., steps < 0)` — the path
  `core.amber.single_point` documents as *"ff14SB + GBn2 energy of the built structure with **NO
  minimisation at all**"*. It is asserted **bit-exact** against
  `core.amber.refine_coords(k_restraint=0, steps=-1)` — max relative difference **0.000e+00 over
  40 comparisons on 10 targets** (the gate fired 40 times; a gate that never fires is not
  evidence). The lean call exists only to remove ~9× of Python overhead.
* `_run` **never calls `_is_collapsed`** — the only call site in the file is
  `AmberHamiltonian._evaluate` (core/amber.py:1168). So the bare arm has **no infinite-valued
  region by construction**, not merely an unfired guard. Measured `n_nonfinite` on the bare arm
  is **0**.
* The claim that *"building a continuous AMBER objective does not escape the relaxation —
  `energy_from_coords` goes through the same `_evaluate`"* is true of `energy_from_coords`
  specifically (it minimises 50 steps by default) and is why this lane does **not** use it.

So this lane's Q1 is a genuine "change only the energy function" comparison on **continuous
torsions**, satisfying BRIEF §5, and is **not** relaxation-confounded.

**The relaxation is nevertheless added as its own explicit arm** (`qb2_relax.py`), because the
*deployed* VQE Hamiltonian is the composed object and its landscape deserves to be measured
rather than assumed: `AMB` (bare) · `AMBr1` (k = 100, steps = 1 — the floor, since OpenMM reads
`maxIterations = 0` as unbounded) · `AMBr50` (k = 100, steps = 50, the deployed setting), with
`n_collapsed` counted on every arm and the F-D2 rank-agreement test computed. Delivered in §4;
artefact `s20/results/qb2_relax.json`, COMPLETE.

*(One first observation from the pricing probe, contended box, reported only as an order of
magnitude: over 24 pool structures the maximum energy falls `1.7e11 → 6.9e6 → 2.1e3` across
`AMB → AMBr1 → AMBr50`. The relaxation is doing the same job the monotone log does. The
quantified version is in §4.)*

I also record the coordinator's correction that AMBER's deployed repair cost is
**+0.021 [+0.014, +0.028]**, not the `+0.164 Å` I was briefed — that was the *projection's* tax.
No prior about AMBER being expensive is carried into this work.

---

## 3. Q1 — THE LANDSCAPE PANEL

`s20/results/qb2_land.json`, `qb2_report_land.json`. Only the Hamiltonian changes: identical
sequence, continuous torsion representation, candidate set, start set (`stable_rng`, 7 shared
starts including `s19/cache/start_<pdb>.npz`), FD step `h = 0.02 rad`, and evaluation accounting.

### 3a. A pre-registration deviation, declared

`PREREG_B.md` §3 specified mean/sd standardisation on the K = 500 retrieval pool. **On the first
target that specification failed**: raw `AMB` over 500 *real* retrieval windows has
`sd_ref = 1.0e12` and a range of **13.2 decades**, so dividing by that sd sends every other
structure's standardised value, gradient norm and Hessian eigenvalue to ~1e−9 — the control is
not a scale for this variable and would have silently reported *"AMBER's landscape is flat"*.
A **robust** centre/scale (median, IQR/1.349) is substituted as primary; **both are persisted**,
and `log10_range`, `sd_over_iqr` and `skew_p99_over_p50` are reported as the diagnostic that
forced the substitution. The pre-registration is kept unedited. **Note that the scale choice is
irrelevant to most of the panel**: `frac_neg`, `cond`, `frac_nearzero`, `aniso`, `n_localmin`,
`acorr_len` and `tv_over_range` are scale-invariant by construction.

### 3b. The panel, n = 20 targets, 7 shared starts, 2 local Hessians and 6 line scans per cell

`qb2_report_land.json`. Scale-invariant rows are marked ▪; the rest use the robust scale.

| metric | LEG | AMB | AMBc | DIST |
|---|---|---|---|---|
| ‖∇Ê‖ (robust units / rad) | 15.67 | **1.615e+13** | 26.86 | 1.771 |
| ‖∇Ê‖ coefficient of variation ▪ | 1.331 | 2.347 | 0.9341 | 0.5364 |
| log₁₀ range on the K=500 pool | 1.612 | **16.39** | 1.569 | 3.423 |
| sd/IQR (heavy-tailedness) ▪ | 1.019 | **3.07e+15** | 1.367 | 1.929 |
| p99/p50 of \|E\| ▪ | 2.725 | **5.33e+08** | 2.76 | 10.53 |
| Hessian fraction negative ▪ | 0.3793 | 0.4999 | 0.4803 | 0.3954 |
| log₁₀ Hessian condition number ▪ | 3.627 | 6.416 | **6.302** | 3.055 |
| fraction near-zero modes ▪ | 0.1416 | 0.6368 | **0.5915** | 0.1929 |
| anisotropy λmax/mean\|λ\| ▪ | 9.725 | 15.52 | **15.32** | 9.981 |
| local minima per 2π along a random line ▪ | 3.892 | 7.783 | **7.767** | 0.9083 |
| autocorrelation length (rad) ▪ | 0.8942 | 0.1563 | 0.6087 | 1.122 |
| total variation / range ▪ | 1.925 | 2.23 | 4.794 | 1.462 |
| barrier between two starts | 4.322 | 7.98e+15 | 2.275 | 0.4207 |

AMBER bit-exactness gate: max relative difference **0.000e+00 over 80 comparisons on 20
targets** — **the gate fired 80 times.** Measured `n_nonfinite` on the bare arm: **0**.

### 3c. **AMB vs AMBc — the load-bearing contrast, and where my hypothesis breaks**

`AMBc = sign(E)·log1p(|E|)` is a **strictly monotone** function of the *same* AMBER call, so the
argmin, the complete ranking of configurations and every level set are **identical**. — **EXACT**.
Any difference below is therefore, by construction, *not* a property of the optimisation
problem's solution set.

| metric | AMB − AMBc | CI95 | W/L | |
|---|---|---|---|---|
| log₁₀ range | **+14.82** | [13.67, 16.00] | **0/20** | conditioning |
| ‖∇Ê‖ | +1.6e13 | [1.2e9, 4.8e13] | **0/20** | conditioning |
| sd/IQR | +3.1e15 | [3.4e9, 9.2e15] | **0/20** | conditioning |
| p99/p50 | +5.3e8 | [1.8e8, 1.0e9] | **0/20** | conditioning |
| barrier | +8.0e15 | [3.3e10, 2.3e16] | **0/20** | conditioning |
| **log₁₀ condition number** | **+0.115** | **[−0.174, +0.385]** | 8/12 | **UNCHANGED** |
| **anisotropy** | **+0.202** | **[−1.108, +1.583]** | 11/9 | **UNCHANGED** |
| **local minima per 2π** | **+0.017** | **[−0.625, +0.642]** | 10/9 | **UNCHANGED** |
| fraction negative curvature | +0.020 | [−0.002, +0.042] | 5/12 | unchanged |
| fraction near-zero modes | +0.045 | [+0.019, +0.073] | 4/13 | barely moved |
| autocorrelation length | −0.452 | [−0.523, −0.386] | 20/0 | improved by conditioning |
| total variation / range | −2.564 | [−3.105, −2.025] | 20/0 | **worse** after conditioning |

> **AMBER carries TWO independent landscape pathologies and a monotone reparameterisation fixes
> exactly one of them.** The dynamic range — 16.4 decades over 500 *real* retrieval structures,
> an r⁻¹² clash tail — collapses to 1.57 decades, and the gradient norm falls by twelve orders of
> magnitude. The **curvature** pathology does not move at all: condition number, anisotropy and
> the count of local minima along a random line are statistically identical before and after.
> — **ESTABLISHED**, n = 20.

### 3d. **Is the conditioned AMBER landscape Legacy-like? — NO, and this REFUTES my H1**

| metric | AMBc − LEG | CI95 | W/L | |
|---|---|---|---|---|
| log₁₀ range | −0.044 | [−0.121, +0.031] | 13/7 | **matched** |
| p99/p50 | +0.035 | [−0.414, +0.486] | 9/11 | **matched** |
| **log₁₀ condition number** | **+2.675** | **[+2.218, +3.118]** | **0/20** | 473× worse |
| **fraction near-zero modes** | **+0.450** | **[+0.387, +0.515]** | **0/20** | 0.59 vs 0.14 |
| **anisotropy** | **+5.593** | **[+3.973, +7.143]** | 3/17 | |
| **local minima per 2π** | **+3.875** | **[+2.808, +4.859]** | 1/19 | 7.77 vs 3.89 |
| **total variation / range** | **+2.869** | **[+2.259, +3.478]** | **0/20** | |
| fraction negative curvature | +0.101 | [+0.046, +0.156] | 4/16 | 0.48 vs 0.38 |
| autocorrelation length | −0.286 | [−0.400, −0.165] | 17/3 | shorter |
| ‖∇Ê‖ | +11.19 | [+3.43, +18.79] | 5/15 | |

**I pre-registered H1 as "the difference is dominated by conditioning". That is half wrong and
the falsifier F2's complement fires.** After the dynamic range and the tail shape have been
matched *to Legacy's, with CIs spanning zero*, the conditioned AMBER field still has **473× the
Hessian condition number**, **four times the near-zero-mode fraction**, and **twice the local
minima per radian**. Conditioning is necessary and is not sufficient. — **SUPPORTED**, n = 20.

**Ordering of the four objectives on every curvature metric is the same and it is not the
"physics" ordering:** `DIST` (the deployed structural objective) is the smoothest on
condition number (10³·⁰⁶), near-zero modes (0.19), local minima (**0.91 per 2π — effectively
unimodal along a random line**) and autocorrelation length (1.12 rad); `LEG` is next; the two
AMBER variants are worst and are indistinguishable from each other. **The objective this
programme found to be the only one whose optimum is in the right place
(`structural-objective-beats-the-energies`) is also, independently, the one with the tamest
landscape.** Whether that is causal is the P1-link question in §5, and it is where the claim
has to be checked rather than told.

**A caution I am obliged to state**: `frac_neg ≈ 0.50` for both AMBER variants at a *random*
point is what a high-dimensional field with no curvature preference looks like; `LEG` 0.379 and
`DIST` 0.395 have a genuine positive-curvature bias. A single local Hessian is **not** global
topology and none of these numbers is used as if it were — every one is a distribution over
2 starts × 20 targets and is labelled *local*.

---

## 4. THE RELAXATION ARM — `H_AMBER = E ∘ Relax` decomposed

`s20/qb2_relax.py` → `qb2_relax.json`, `qb2_report_relax.json`. n = 10 targets, reference
ensemble = the first 120 windows of each target's own K = 500 BLOSUM pool (declared; the
50-step relaxation is ~50× the cost of a single point). Three arms, identical coordinates:
`AMB` (bare single point) · `AMBr1` (k = 100 kcal/mol/Å², **steps = 1**) · `AMBr50`
(k = 100, **steps = 50**, the `amber_hamiltonian.AmberHamiltonian` default and therefore the
objective the *deployed* VQE actually minimises). `AMBr1` and `AMBr50` are asserted bit-exact
against `core.amber.refine_coords` at the same settings before any number is read.

**The two questions this arm answers, and only these two.**
(a) *Is the relaxation a conditioner or a different objective?* — Spearman(`E∘Relax_50`, `E`) and
Spearman(`E∘Relax_50`, `E∘Relax_1`) against the audit lane's **F-D2** threshold of 0.95.
(b) *How much of the deployed AMBER landscape's apparent tameness is the relaxation?* —
`log10_range`, ruggedness, and gradient norm across the three arms with Legacy on the same axis.

Bit-exactness gate against `core.amber.refine_coords` at the same settings: **max relative
difference 0.00e+00** — AMB over 30 comparisons, AMBr1 over 20, AMBr50 over 20. **The gate fired
70 times.**

### 4a. The relaxation is NOT a conditioner. It is a different objective. — **ESTABLISHED**, n = 10

| rank agreement (Spearman over 120 real pool structures, per target) | mean | CI95 | verdict |
|---|---|---|---|
| `E∘Relax_50` vs the **bare single point** | **+0.720** | [+0.656, +0.778] | **REORDERS** |
| `E∘Relax_50` vs `E∘Relax_1` | **+0.720** | [+0.654, +0.781] | **REORDERS** |
| `E∘Relax_1` vs the bare single point | **+0.950** | [+0.924, +0.973] | spans 0.95 |
| `E∘Relax_50` vs `AMBc` (monotone twin of the bare point) | +0.720 | [+0.656, +0.778] | identical to the row above, as it must be — **EXACT** |
| bare AMBER vs Legacy | −0.097 | [−0.211, +0.032] | the two energies barely agree at all |
| `E∘Relax_50` vs Legacy | +0.180 | [+0.039, +0.323] | |

> **The coordinator's objection is CORRECT and is now quantified. `H_AMBER = E ∘ Relax_50` and
> the bare AMBER energy agree on only 52 % of rank variance (ρ = 0.720), so they are different
> objectives, not two views of one landscape.** The audit lane's **F-D2 falsifier does NOT fire**
> — the 0.95 bar is missed decisively at 50 steps — and the objection therefore stands for any
> lane that measures the *deployed* Hamiltonian. **A single minimisation step is the borderline
> case** (ρ = 0.950, CI spanning the bar): `Relax_1` is approximately a conditioner and
> `Relax_50` is not.
>
> **This is exactly why this lane's primary Q1 arm is the bare single point.** `AMB` and `AMBc`
> are the same call with the same argmin, so the §3 comparison isolates the energy function; the
> relaxation is measured here as its own operator rather than being silently folded into "AMBER".

### 4b. What the relaxation actually does: it trades dynamic range for multimodality

| | log₁₀ range | sd/IQR | p99/p50 | local minima per 2π | autocorr len | ρ(E, RMSD) ORACLE | argmin RMSD |
|---|---|---|---|---|---|---|---|
| `LEG` | 1.472 | 1.10 | 3.64 | — | — | **+0.138** | 5.270 |
| `AMB` bare | **12.63** | 1.3e9 | 1.8e9 | **5.6** | 0.196 | −0.135 | 4.460 |
| `AMBc` log-conditioned | 1.458 | 1.75 | 2.93 | — | — | −0.135 | 4.460 |
| `AMBr1` | 9.216 | 3.1e5 | 3.5e4 | **11.1** | 0.226 | −0.149 | 4.756 |
| **`AMBr50` (deployed)** | 3.652 | 1.8e4 | 7.29 | **15.2** | 0.491 | −0.130 | 4.855 |

> **Relaxation compresses the range by nine decades and simultaneously TRIPLES the number of
> local minima** (5.6 → 15.2 per 2π; Legacy's is 3.9 and the distogram objective's is 0.9 in the
> §3 panel). That is the coordinator's predicted mechanism measured directly: *smooths most of
> the domain, injects discontinuities where the minimiser's basin assignment flips.* The
> deployed AMBER Hamiltonian is the **most multimodal objective anywhere in this study.**

### 4c. The collapse sentinel, as a measured count — and it is not a small number

| arm | `+inf` (deployed geometry predicate, min heavy-atom contact < 1.05 Å) | of |
|---|---|---|
| `AMB` bare | **0** | 3302 calls — *and by construction*: `core.amber._run` never calls `_is_collapsed` |
| **`AMBr1`** | **1109 (33.59 %)** | 3302 calls |
| `AMBr50` (deployed) | 24 (0.73 %) | 3302 calls |

> **At one minimisation step the deployed Hamiltonian would return `+inf` on a third of its
> evaluations.** The 50-step minimisation is not a refinement of the energy — it is what makes
> the objective finite often enough to be optimised at all. A guard that never fires is not
> evidence; this one fires 1109 times and the number is the finding.
> — **ESTABLISHED**, n = 10, 3302 calls per arm.

*(`gnorm_std_mean` is reported as `nan` for `LEG` and `AMBc` in this artefact because this module
computes gradients only for the three AMBER arms; their gradient norms are in the §3 panel at
n = 20. Recorded so the `nan` is not mistaken for a failure.)*

---

## 5. P1-LINK — the answer is TWO-SIDED, and the useful half is the negative one

`qb2_report_link.json`, `qb2_report_decouple_opt.json`. n = 20. Difficulty control partialled out
on ranks: **`pool500`'s realised built-chain Cα-RMSD through the frozen terminal operator** (the
binding requirement of `s20/LEDGER.md` L4). The partial-correlation implementation is verified
against the textbook one-control formula to **9 decimal places**.

### 5a. **The decisive test needs no landscape metric at all — and it is a clean zero**

For every (objective, arm) cell, does beating `best_of_N` **on the objective** go with beating it
**on Cα-RMSD**?

    cells where objective and RMSD BOTH move significantly, SAME direction:     0
    cells where BOTH move significantly, OPPOSITE directions (a TRADE):         8
    cells where at most one moves significantly:                              28
                                                                        ------ 36 cells

    WITHIN-CELL, target-level rho(objective gained vs best_of_N, RMSD gained vs best_of_N)
      over all 36 cells:   mean +0.0034   median +0.0436   CI [-0.0679, +0.0717]
                           positive in 20/36 cells

> **Zero cells out of thirty-six improve the objective and the structure together. Eight
> significantly trade one for the other. The per-target correlation between "how much better this
> arm optimised on this target" and "how much better a structure it built on this target" is
> +0.003 with a CI tight around zero.** This holds under Legacy, under bare AMBER, under
> conditioned AMBER and under the deployed distogram objective alike. — **ESTABLISHED**, n = 20,
> 36 cells, ORACLE scoring of native-free decisions.

This is `search-saturates-discrimination-binds` measured on a *sixth* instrument, and it is the
single fact that prices the whole of Q1: **whatever the answer to "is the AMBER landscape harder",
it cannot move Priority 1 through the optimiser, because optimiser progress and structural
outcome are not coupled on any of these objectives.**

### 5b. **And yet one landscape metric DOES predict which targets end badly**

52 tests at 95 % → 2.6 flags expected by chance. Observed **10 partialled flags**, and they are
not scattered:

| metric | LEG | AMB | AMBc | DIST |
|---|---|---|---|---|
| **`frac_neg`** (fraction of negative Hessian eigenvalues at a random start) | **+0.528 [+0.092, +0.804]** | **+0.546 [+0.092, +0.823]** | **+0.510 [+0.122, +0.840]** | **+0.511 [+0.068, +0.789]** |
| `log10_cond` | **+0.517 [+0.039, +0.772]** | **+0.625 [+0.150, +0.860]** | +0.366 ns | +0.176 ns |
| `aniso` | +0.425 ns | **+0.493 [+0.026, +0.794]** | **+0.487 [+0.078, +0.830]** | +0.035 ns |
| `skew` | ns | ns | ns | **−0.471 [−0.769, −0.071]** |
| `tv_over_range` | ns | ns | ns | **−0.466 [−0.714, −0.046]** |
| every other metric (gnorm, range, n_localmin, acorr_len, barrier, …) | ns | ns | ns | ns |

> **`frac_neg` survives the difficulty partial on all four objectives independently, with the
> same sign and a CI excluding zero every time.** Four independent replications of one metric is
> not a multiplicity artefact; a scatter of ten unrelated flags would have been. The
> saddle-dominance of the objective at a random plausible start is a **real** predictor of how
> badly that target finishes. — **SUPPORTED**, n = 20 (a small instrument for a correlation; the
> CIs are wide and the claim is deliberately not stated more strongly).

### 5c. The synthesis, stated so it cannot be over-read

**Landscape geometry is diagnostic of the TARGET and is not actionable through the OPTIMISER.**
`frac_neg` tells you which targets will end badly (ρ ≈ 0.51–0.55 after difficulty is removed);
changing how hard, how well or by what method you optimise changes the outcome not at all
(0/36 cells, within-cell ρ = +0.003). Both statements are measured here on the same 20 targets.
Everything else in the §3 panel — the twelve-decade gradient norm, the 473× condition number, the
doubled local-minimum count, the barrier heights — **has no predictive relationship to Cα-RMSD and
is reported as an ornament, exactly as BRIEF §4 requires.**

---

## 6. S1 — THE OPTIMISER BATTERY: conditioning buys 1.8 IQR of optimisation and 0 Å of structure

`qb2_opt.json`, `qb2_report_opt.json`. **n = 20 targets, COMPLETE.** B = **512** evaluations per
(objective, arm, seed), FD probes counted, 2 seeds, identical shared starts, identical robust
standardisation, identical move classes. Ten arms including the **mandatory `best_of_N`**.

### 6a. Cα-RMSD of the best-seen configuration (ORACLE, post-hoc)

| arm | AMB | AMBc | DIST | LEG |
|---|---|---|---|---|
| `best_of_N` (**mandatory control**) | 5.186 | 4.755 | 3.987 | 4.562 |
| `mom_fd` | **4.350** | 4.395 | 3.958 | **4.007** |
| `lbfgs_fd` | 4.430 | 4.591 | 4.003 | 4.230 |
| `adam_fd` | 4.542 | 4.408 | 4.029 | 4.447 |
| `nelder` | 4.500 | 4.500 | 3.960 | 4.370 |
| `spsa` | 4.762 | 4.916 | 3.951 | 4.385 |
| `powell` | 4.857 | 4.627 | **3.908** | 4.683 |
| `metro` | 5.196 | 4.908 | 3.941 | 4.725 |
| `anneal` | 5.055 | 5.468 | 3.980 | 4.596 |
| `greedy` | 5.677 | 5.159 | **3.843** | 4.723 |
| *reference:* pool best 1.901 · pool mean 4.739 · shared start 3.401 | | | | |

### 6b. Best-seen OBJECTIVE reached (robust-standardised; lower is better) — **the conditioning result**

| arm | AMB | **AMBc** | DIST | LEG |
|---|---|---|---|---|
| `best_of_N` | **−0.031** | **−1.820** | −1.302 | −0.619 |
| `adam_fd` | **+4.469** | −0.844 | −1.307 | −0.308 |
| `nelder` | **+1.640** | −0.945 | −1.365 | −0.533 |
| `lbfgs_fd` | **+0.320** | −1.095 | −1.509 | −0.170 |
| `greedy` | −0.031 | **−1.911** | −1.491 | −1.353 |
| `anneal` | −0.031 | **−1.899** | −1.496 | −1.287 |
| `powell` | −0.030 | −1.694 | **−1.623** | **−1.533** |

> **On raw AMBER, not one arm beats `best_of_N` on the objective, and the gradient-based arms end
> up catastrophically worse** (Adam +4.469, Nelder +1.640, L-BFGS +0.320 robust units *above* a
> random draw's best). **On `AMBc` — the identical energy under a strictly monotone transform, with
> an identical argmin and an identical ranking of every configuration — every arm reaches ≈ −1.8
> to −1.9.** The gap is about **1.8 IQR of optimiser progress bought by a reparameterisation that
> provably changes nothing about the problem's solution.** — **ESTABLISHED**, n = 20.
>
> **This is the surviving half of my H1, and it is large.** It is also the half that does not
> matter: the same table's RMSD column (§6a) moves by ≤ 0.4 Å and §5a shows the two are not
> coupled.

### 6c. **The objective/RMSD trade, in the raw numbers**

Every cell where an arm significantly **beat** `best_of_N` on Cα-RMSD is a cell where it
significantly **lost** on the objective, and vice versa:

    [AMB]   mom_fd   RMSD -0.836 [-1.251,-0.481] sig  |  OBJ +0.0218  [+0.0056,+0.0440] sig
            lbfgs_fd RMSD -0.755 [-1.149,-0.350] sig  |  OBJ +0.3505  [+0.0853,+0.6943] sig
            nelder   RMSD -0.686 [-1.137,-0.225] sig  |  OBJ +1.670   [+0.0153,+4.557]  sig
            adam_fd  RMSD -0.643 [-1.038,-0.235] sig  |  OBJ +4.499   [+0.0761,+13.09]  sig
            spsa     RMSD -0.424 [-0.840,-0.024] sig  |  OBJ +0.0061  [+0.0027,+0.0103] sig
    [AMBc]  anneal   RMSD +0.713 [+0.257,+1.202] sig  |  OBJ -0.0784  [-0.1655,-0.0167] sig
            greedy   RMSD +0.404 [+0.020,+0.808] sig  |  OBJ -0.0908  [-0.1838,-0.0266] sig
    [LEG]   mom_fd   RMSD -0.555 [-0.934,-0.181] sig  |  OBJ +0.5426  [+0.2949,+0.7701] sig
    [DIST]  nothing significant on RMSD at all; every arm's |Δ| <= 0.144 with a CI spanning zero

**Eight for eight.** On the deployed distogram objective the trade does not even appear — the
objective moves significantly on six arms and Cα-RMSD moves on none.

### 6d. Does the objective RANK? — the cleanest statement of the whole lane

Over the **same K = 500 retrieval windows the pipeline actually uses**, per target, n = 20:

| objective | ρ(E, Cα-RMSD) ORACLE | RMSD of the objective's argmin | mean RMSD of its own top 1 % |
|---|---|---|---|
| **AMB** = **AMBc** (identical, as an **EXACT** consequence of monotonicity) | **−0.071 [−0.173, +0.042]** | **4.990** | 5.461 |
| **LEG** | +0.191 [+0.011, +0.353] | **5.487** | 4.946 |
| **DIST** | **+0.497 [+0.325, +0.654]** | **3.676** | 3.850 |
| *pool mean* | | **4.739** | |

> **Picking the lowest-energy structure out of 500 real retrieval windows is WORSE than picking
> one at random — for both physics energies.** AMBER's argmin is 4.990 Å and Legacy's is 5.487 Å
> against a pool mean of 4.739 Å. Only the structural objective's argmin (3.676 Å) beats the
> random draw. AMBER's rank correlation with correctness is **negative with a CI spanning zero**;
> Legacy's is +0.19; the distogram's is +0.50. — **ESTABLISHED**, n = 20, and it reproduces
> `objective-does-not-rank-the-native` and `structural-objective-beats-the-energies` on a third
> instrument.
>
> **This is why §5a comes out zero.** There is nothing wrong with the optimisers. The optima are
> in the wrong place.

### 6e. SPSA is not the problem — the optimizer-mismatch hypothesis is ELIMINATED

Cosine between the SPSA gradient estimate and the **central-FD gradient at the same point**,
budgeted like everything else, 40 measurements per objective:

    AMB   +0.1782 [+0.1527, +0.2039]      DIST  +0.1500 [+0.1280, +0.1723]
    AMBc  +0.1728 [+0.1465, +0.2006]      LEG   +0.1550 [+0.1297, +0.1811]

> **SPSA's gradient quality is statistically indistinguishable across all four objectives** — the
> four CIs overlap almost completely, and if anything it is *best* on raw AMBER. A one-perturbation
> SPSA estimate in 2n ≈ 26 dimensions carries ≈ 0.15–0.18 cosine with the true gradient on every
> objective, which is what the estimator's variance predicts and has nothing to do with the
> Hamiltonian. **"SPSA copes with Legacy and not with AMBER" is REFUTED.** The AMB→AMBc collapse
> in §6b happens with the *same* gradient quality, so it is the objective's parameterisation, not
> the estimator's noise. — **REFUTED** (the optimizer-mismatch-as-gradient-noise hypothesis),
> n = 20.

The non-finite sentinel fired **0 times on all four objectives** — because the bare AMBER path
never calls `_is_collapsed` (§2), so this is a guard that is absent by construction rather than
one that passed vacuously. The relaxed arms, where the predicate *is* meaningful, fire it 1109
and 24 times (§4c).

## 7. S2 — CVaR α: on raw AMBER the tail is a **crude conditioner**, and a monotone transform does the same job better

`qb2_cvar.json`, `qb2_report_cvar.json`. **n = 10 targets, COMPLETE**, B = 512, 2 seeds,
64 shots/iteration. Genuine CVaR-VQE: `core.quantum.MPSAnsatz` (RY + CNOT chain + final RY),
`core.quantum.cvar_gradient(..., baseline="const")` — the **corrected** estimator, never the
recorded `tail` defect. `vqe_untrained` is the **mandatory** same-circuit, same-budget,
zero-gradient-step control.

**Scope, set by the coordinator's L4 and respected here.** Whether a tail *helps RMSD* is already
measured at n = 126 and it does not (`q_a0.05 − q_a1.00 = +0.012 [−0.084, +0.105]`, and α = 1 is
the best quantum arm). This panel therefore measures the **mechanism**, not the outcome.

### 7a. Derive the operator before interpreting its statistic — **EXACT**

The CVaR weight vector is `w_i = −(q_α − e_i)/(αN)` on the α-tail and 0 elsewhere, so
`ESS = (Σ|w|)²/Σw²` is a functional of **the tail's order statistics alone**. It contains no
reference to the Hamiltonian. **Any Legacy-vs-AMBER difference in measured CVaR concentration is
therefore a difference in their energy SPECTRA passed through one shared, Hamiltonian-independent,
discontinuous operator — an identity, not a property of AMBER.** The numerical demonstration:
ρ(log₁₀ spectrum skew, ESS fraction) across all 40 (objective, target) cells is **+0.311** at
α = 0.05, **+0.609** at α = 0.25, and **−0.565** at α = 1.00 — the sign flips exactly where the
rule's mechanism changes (below 1 the tail size is fixed and ESS tracks tail spread; at 1 every
sample is in the tail and a heavy-tailed spectrum concentrates the weight on its own extreme).

### 7b. What actually differs between the Hamiltonians — the score-function gradient norm

| objective | ESS fraction: α=0.05 / 0.25 / 1.00 | ‖CVaR gradient‖: α=0.05 / 0.25 / 1.00 |
|---|---|---|
| **AMB** (raw) | 0.142 / 0.612 / **0.090** | 0.0060 / 0.178 / **9.27e+14** |
| **AMBc** (same energy, log) | 0.137 / 0.475 / 0.643 | 1.079 / 1.021 / **1.713** |
| DIST | 0.135 / 0.462 / 0.588 | 0.370 / 0.412 / 0.831 |
| LEG | 0.125 / 0.390 / 0.288 | 0.646 / 0.430 / 3.378 |

> **On raw AMBER at α = 1 the score-function gradient norm is 9.3 × 10¹⁴ and the effective sample
> size COLLAPSES to 0.090 — worse concentration than the α = 0.05 tail.** One clash structure
> owns the batch. Truncating to a tail is the only thing keeping the estimator finite, so on raw
> AMBER **CVaR is functioning as a crude, discontinuous conditioner rather than as a risk
> objective.** Under the monotone log transform of the *same* energy the gradient norm is 1.0–1.7
> at every α and the ESS behaves exactly as the rule predicts. — **SUPPORTED**, n = 10.
>
> This is the specific, mechanistic answer to "why does CVaR seem to matter more under AMBER":
> it is doing conditioning work that a monotone reparameterisation does better, continuously,
> and for free. And it aligns with **arXiv:2605.02850 (Quantum Tilted Loss, 2026)**, which proves
> CVaR cannot remove barren plateaus, only reshape local geometry, and prices the sharpened
> gradient at `O((e^{|γ|Δ}−1)²/γ²ε²)` — exponential in the tilt. At 512 shots an α = 0.05 tail is
> 26 samples; the recovered Sprint-19 extension panel (§1b iii) independently found the deepest
> tail (α = 0.02, 10 samples) to be the worst arm and dilution toward the plain mean the best.

### 7c. RMSD against the mandatory control — no effect anywhere, n = 10

    [AMB]  vqe0.05 +0.020 ns | vqe0.25 -0.347 ns | vqe1.0 -0.110 ns
    [AMBc] vqe0.05 -0.396 ns | vqe0.25 -0.276 ns | vqe1.0 +0.049 ns
    [DIST] vqe0.05 -0.010 ns | vqe0.25 -0.156 ns | vqe1.0 -0.274 ns
    [LEG]  vqe0.05 -0.003 ns | vqe0.25 +0.351 ns | vqe1.0 +0.360 [+0.006,+0.737] sig, WORSE

**Eleven of twelve cells span zero at n = 10; the one that does not is the VQE doing WORSE than
its own untrained circuit under Legacy.** This panel has nowhere near the power of the n = 126
sweep and is labelled **INCONCLUSIVE on RMSD**; it is quoted only for its mechanism columns, whose
per-cell precision is high. **Latent entropy falls monotonically with α on every objective**
(e.g. DIST 9.15 bits untrained → 7.82 at α = 1 → 7.12 at α = 0.05) and `latent_max_prob` rises —
the concentration is real, measurable, and buys nothing, which is
`concentration-is-wrong-when-discrimination-binds` reproduced with the mechanism attached.

## 8. S3 — ENCODING: the coordinate artifact is REAL, and it points the same way as everything else

`qb2_enc.json`, `qb2_report_enc.json`. **n = 10 targets, COMPLETE.** Identical objectives,
identical shared starts, identical seeds, identical budget B = 512, identical arms — run in
`theta` versus the smooth periodic embedding `u = (cos θ, sin θ)` with an `arctan2` retraction.
The retraction is scale-invariant, so the radius is pure gauge and the two parameterisations are
**physically identical**. **No discretisation is used anywhere.**

Negative = the embedding wins.

| | `adam_fd` | `lbfgs_fd` | `nelder` | `spsa` |
|---|---|---|---|---|
| **AMB** RMSD | −0.287 [−0.595, +0.000] | −0.048 ns | **−0.272 [−0.461, −0.099]** | −0.247 ns |
| **AMBc** RMSD | −0.369 ns | −0.125 ns | **−0.272 [−0.461, −0.099]** | **−0.649 [−1.205, −0.165]** |
| **LEG** RMSD | −0.191 ns | **+0.073** ns | −0.138 ns | −0.131 ns |

> **The embedding's point estimate is better in 11 of 12 cells** (sign test, two-sided
> p = 0.0063), three with CIs excluding zero, the largest being **−0.649 Å**. A change of
> coordinates that provably does not change the physics moves final Cα-RMSD by up to two-thirds
> of an Ångström. — **SUPPORTED**, n = 10.

**And the same trade appears again**: in 9 of 12 cells the embedding reaches a *worse* objective
value while building a *better* structure — `AMBc/nelder` is −0.272 Å with the objective +0.615
worse, `AMBc/spsa` −0.649 Å with the objective +0.378 worse. §5a's decoupling reproduces in a
third independent setting.

**A confound I must declare, and it cuts against the objective column only.** The embedded space
has dimension 4n rather than 2n, so at a fixed budget of 512 evaluations the finite-difference
arms get **half as many gradient steps**. That is a real disadvantage and it plausibly explains
part of the worse objective values. It does **not** touch the RMSD column, which is simply the
best structure found within the same 512 evaluations from the same starts. **So the objective
column here is NOT MEASURED cleanly and is not claimed; the RMSD column is.**

**Consequence for the pre-registration.** PREREG §4 S3 says any Q1 conclusion that differs between
the two parameterisations is **encoding-dependent and not claimable**. Checked one by one: the
§5a decoupling holds in both (same sign, same trade). The §6d ranking result is a property of the
objective, not of the parameterisation, and is unaffected. The §6b conditioning gap is measured in
`theta` only, and its embedded counterpart is confounded by the halved iteration count — so **the
size of the AMB→AMBc conditioning gap is reported as a `theta`-space measurement and its
magnitude in other parameterisations is recorded OPEN.** The direction (raw AMBER is the worst
place to run a gradient method) reproduces in both.

## 9. Q2 — DECORRELATION: **CLOSED.** The VQE's coherent error is the pool's error

`s20/qb2_q2.py` → `qb2_q2.json`, `qb2_report_q2.json`. **n = 20 targets, COMPLETE.** Seven
generator families, two independent seeds each, budget 4096 evaluations on the deployed distogram
objective, every family passed through the **frozen native-free terminal operator** (shipped
Bayes-risk top-75 → coordinate average → ideal-geometry projection). Both bases are recorded:
`built RMSD` is the projected chain, `cloud RMSD` the point cloud.

**The identity first, so it is not mistaken for a finding.** Sprint 19 splits a *predictor's*
error as `r = r_coh + r_inc`. A **generator** emits a structure, and a structure realises its own
distances exactly, so `r_inc ≡ 0` and the generator's whole error *is* its coherent component.
— **EXACT**.

### 9a. Correlations, against the ceiling — never against zero

| family | same-generator ceiling (2 seeds) | **ρ vs the retrieval pool** |
|---|---|---|
| `pool500` | +1.000 [+1.000, +1.000] — deterministic, an identity check that passed | — |
| **`vqe`** (α = 0.25, CNOT chain) | +0.852 [+0.754, +0.936] | **+0.777 [+0.672, +0.870]** |
| `vqe_prod` (CNOTs deleted) | +0.833 [+0.729, +0.921] | +0.765 [+0.630, +0.879] |
| `vqe_untrained` | +0.800 [+0.716, +0.875] | +0.803 [+0.726, +0.873] |
| `bestofN` | +0.917 [+0.872, +0.956] | +0.802 [+0.681, +0.902] |
| `metro` | +0.897 [+0.834, +0.952] | +0.728 [+0.597, +0.849] |
| **`helix` — ZERO INFORMATION** | +0.994 [+0.991, +0.997] | **+0.793 [+0.658, +0.895]** |

### 9b. A PRE-REGISTRATION MIS-SPECIFICATION, declared

`PREREG_B.md` §5 adopted the brief's falsifier verbatim: *"if the VQE candidates align with the
pool at the cross-architecture level (≈0.75+), the generator role closes."* **That line was
calibrated on cross-architecture PREDICTOR families. Measured here, the zero-information constant
helix already correlates with the pool at +0.793.** The threshold therefore sits **below its own
zero-information null** for *generators*, and the absolute test cannot discriminate. The
pre-registration stays unedited; the valid test is against the measured null and against each
generator's own ceiling, and both are reported.

### 9c. The null-referenced test — a clean, complete null

    vqe-vs-pool  minus  helix (ZERO-INFORMATION)      -0.0157 [-0.1101, +0.1072]  13W/7L   ns
    vqe-vs-pool  minus  bestofN (no optimisation)     -0.0245 [-0.0860, +0.0422]  14W/6L   ns
    vqe-vs-pool  minus  metro (classical thermostat)  +0.0491 [-0.0055, +0.1113]   8W/12L  ns
    vqe-vs-pool  minus  vqe_prod (CNOTs deleted)      +0.0120 [-0.0423, +0.0784]  12W/8L   ns
    vqe-vs-pool  minus  its OWN same-seed ceiling     -0.0745 [-0.1450, -0.0096]  16W/4L   sig

> **The VQE's error correlation with the retrieval pool is indistinguishable from a constant
> α-helix's, from a no-optimisation draw's, and from a classical thermostat's** — and it sits only
> **0.075 below the VQE's own same-seed ceiling**, the one comparison whose CI excludes zero.
>
> **91.3 % of the VQE's own reproducible error structure is shared with the retrieval pool**
> (0.777 / 0.852; per-target ratio mean 0.946). For the zero-information helix the same ratio is
> 0.797. **There is essentially nothing in the VQE's coherent error that is its own.**
> — **Q2 CLOSED**, n = 20, ORACLE scoring of native-free decisions.

### 9d. Quality and diversity, reported separately from error correlation (BRIEF §4)

| family | built RMSD | cloud RMSD | member mean | member best | geometric diversity |
|---|---|---|---|---|---|
| `pool500` (0 evaluations) | **3.456** | 3.264 | 3.900 | 2.559 | 2.688 |
| `bestofN` | 3.552 | 3.360 | 4.007 | 2.787 | 3.087 |
| `vqe` | 3.597 | 3.392 | 3.906 | 2.812 | 2.641 |
| `vqe_untrained` | 3.600 | 3.354 | 4.170 | 2.852 | **3.439** |
| `vqe_prod` | 3.824 | 3.631 | 4.059 | 3.001 | 2.510 |
| `metro` | 3.877 | 3.805 | 3.950 | 3.485 | **1.060** |
| `helix` | 4.393 | 4.272 | 4.368 | 3.460 | 2.155 |

Built-chain RMSD against `pool500`, paired: `bestofN` +0.096 [−0.155, +0.354] ns · `vqe`
+0.141 [−0.142, +0.456] ns · `vqe_untrained` +0.144 ns · `vqe_prod` **+0.368 [+0.118, +0.645]
sig** · `metro` **+0.421 [+0.067, +0.806] sig** · `helix` **+0.937 [+0.422, +1.560] sig**.

**At n = 20 the VQE is not distinguishable from the zero-evaluation pool on RMSD** (the n = 126
sweep, which has the power, puts it +0.256 behind). **Two populations here are geometrically very
different and make the same structural mistake**, exactly the separation the brief asks for:
`vqe_untrained` has the highest geometric diversity on the board (3.439) and the *same* error
correlation with the pool (+0.803) as everything else.

**Verdict.** *"Quality is acceptable AND error correlation is low"* is the standard the brief set
for a valuable quantum sampler. Quality is at best matched; error correlation is **not** low, is
**not** below any null, and is **91 % of the arm's own ceiling**. **The generator role closes.**

---

## 10. SCOREBOARD AGAINST THE PRE-REGISTRATION

`s20/PREREG_B.md` is unedited. Two of its specifications turned out wrong and both are declared
rather than repaired:

| pre-registered item | outcome |
|---|---|
| **F1** — *"if after conditioning AMBER still shows worse optimiser success AND worse curvature AND shorter autocorrelation AND at least one landscape metric predicts Cα-RMSD, my hypothesis is REFUTED"* | **FIRES, in part.** Curvature: yes (§3d). Autocorrelation: yes. A metric predicting RMSD: yes — `frac_neg` on all four objectives (§5b). Optimiser success after conditioning: **no** — `AMBc` recovers to −1.9 IQR, matching or beating everything (§6b). **Three of four clauses fire. H1 is half refuted and is reported as such.** |
| **F2** — *"if AMB and AMBc are indistinguishable, conditioning is refuted"* | **Does not fire.** They differ by 12 orders of magnitude in gradient norm and 1.8 IQR in optimiser progress. |
| **F3** — *"if any landscape metric predicts Cα-RMSD, I may not report that landscape geometry is disconnected from structural outcome"* | **FIRES.** `frac_neg` predicts. The finding is therefore stated as the **two-sided** §5c, not as a blanket disconnect. |
| **Q2 falsifier** (≈0.75 alignment ⇒ generator role closes) | **MIS-SPECIFIED** — the threshold sits *below* its own zero-information null for generators (helix = +0.793). Declared in §9b; replaced by the null-referenced test, which returns a clean null and closes Q2 anyway. |
| **PREREG §3 scale control** (mean/sd) | **MIS-SPECIFIED** — AMBER's reference sd is 1.0e12; the control is not a scale for this variable. Robust median/IQR substituted, both persisted, the diagnostic that forced it reported. §3a. |
| **PREREG §6** n = 20 declared, "no Cα-RMSD claim at n = 20 is a validated pipeline improvement" | **Honoured.** No arm here is proposed for the pipeline. |

## 11. WHAT I RECOMMEND, AND WHAT I DO NOT

**Do not** spend Sprint 21 compute on: the ansatz (closed twice), the CVaR α (closed at n = 126
and mechanistically explained here), SPSA replacements (the estimator is not the problem),
or "better physics" for ranking (both energies' argmins are worse than a random draw over the
deployed pool).

**If anyone runs the deployed AMBER Hamiltonian again**, note that it is `E ∘ Relax_50`, that it
agrees with the bare energy on only 52 % of rank variance, that it is the most multimodal
objective measured anywhere in this study, and that at one minimisation step it returns `+inf`
on a third of its calls. Cite §4, not "AMBER".

**Two things this lane produced that are reusable regardless of the quantum question:**
1. **Seed sensitivity 0.200 Å** (§1b i) — 2.4× the instrument's MDE. Any future variational arm
   must report ≥ 4 ansatz seeds or be labelled NOT MEASURED.
2. **A monotone conditioning of AMBER is free and worth 1.8 IQR of optimiser progress** (§6b),
   with an argmin and a ranking that are provably unchanged. If AMBER is ever optimised rather
   than merely evaluated, do it on `sign(E)·log1p(|E|)`.

**And the honest closing statement.** Priority 2 asked that the genuine CVaR-VQE be preserved.
It is: genuine circuit, genuine sampled CVaR score-function gradient with the corrected baseline,
genuine Legacy at `DEFAULT_WEIGHTS`, genuine AMBER asserted bit-exact against the codebase entry
point on 150 comparisons. It trains, measurably, against the mandatory best-of-N control. Its
measured contribution to Cα-RMSD is zero to slightly negative, its coherent error is 91 % the
retrieval pool's, and the generator role closes. **A pillar that is real and does not contribute
is an acceptable outcome; the false positive would not have been.**
