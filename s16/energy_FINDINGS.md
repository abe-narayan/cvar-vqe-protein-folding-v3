# SPRINT 16 — ENERGY findings

**The two mandated energy pillars, and the causal question of what each actually contributes.**
Genuine Legacy / Miyazawa–Jernigan (`core.energy`, eleven separately-measured components,
reconstructed from those components to **1.47 × 10⁻⁵ over 12.8 × 10⁶ structures**) and genuine
all-atom AMBER ff14SB / GBn2 (`core.amber`, OpenMM). **No learned surrogate stands in for either
anywhere in this workstream.**

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Modules: `s16/energy_lib.py` (gate, standing frame null, monotone conditioning, binding rule,
Legacy detectors, stereochemistry panel), `s16/energy_gate.py`, `s16/energy_ablate.py`,
`s16/energy_report.py`, `s16/energy_repair.py`, `s16/energy_weights.py`,
`s16/energy_contrast.py`, `tests/test_amber_frame_invariance.py`.
Artefacts in `s16/results/`.

---

## 0. VERDICT, AND THE CLAIM LEDGER

> **AMBER contributes stereochemical repair and 0.023 Å of accuracy. Legacy contributes detection
> and nothing to accuracy. They do not interact. Neither ranks a native, and fitting Legacy's
> weights leave-fold-out makes it worse — so the ranking failure is not a training-status artefact.
> A distance-restraint objective on the same configurations at the same chain lengths is −0.926 Å
> better than random, which locates the pathology in the potential class rather than in the
> problem.**
>
> Two corrections run **against** this workstream and are reported in full: the AMBER minimiser's
> rigid-invariance null is **+0.0146 Å live on 11 converged targets**, larger than the accuracy
> effect it certifies; and the published Ramachandran gain is **2.39× overstated** because the
> incumbent scores unwrapped torsions.

| claim | tier |
|---|---|
| the four fixes are in place: gate, standing frame null, monotone conditioning, binding rule | **DEMONSTRATED** (§2) |
| AMBER improves stereochemical validity, hugely, on every steric and Ramachandran axis | **DEMONSTRATED** (§1.2, §3.1) |
| AMBER improves accuracy by −0.0234 Å [−0.0375, −0.0099] as an aggregate, gated, n = 123 | **DEMONSTRATED**, with §2.2 attached |
| no per-target AMBER RMSD claim is resolvable — the exact frame null's residual is 5–7× the effect | **DEMONSTRATED** (§2.2) |
| **AMBER introduces cis peptide bonds on 44 / 126 targets** | **DEMONSTRATED** (new defect, §3.1) |
| Legacy contributes nothing to accuracy; its effect equals the cost of a matched random filter | **DEMONSTRATED** (§1.1) |
| Legacy and AMBER do not interact (+0.007 [−0.007, +0.022]) | **DEMONSTRATED** (§1.1) |
| Legacy detects stereochemical defects at AUROC 0.93–0.99 through 2 of its 11 components | **DEMONSTRATED** (§3.2) |
| the eleven-term Legacy total is worse than its own best component on all 10 detection axes | **DEMONSTRATED** (§3.2) |
| leave-fold-out fitting of Legacy's weights makes its certified optimum +0.7785 Å WORSE | **DEMONSTRATED / REFUTED** (my own hypothesis, §4) |
| neither energy ranks the native; they are indistinguishable from each other | **DEMONSTRATED** (§5) |
| the distance objective's certified argmin is −0.926 Å better than random on the same spaces | **DEMONSTRATED** (§5) |
| the published arm-A Ramachandran fraction is an unwrapped-torsion artefact | **DEMONSTRATED** (new defect, §3.4) |
| my own pre-declared PASS band for the standing frame null | **FAILED, reported as failed** (§2.2) |
| Legacy's ranking failure is a training-status artefact | **REFUTED** (§4) |
| the MJ contact term is the removable seat of the anti-ranking | **REFUTED** (§5.4) |
| a Legacy of `steric` + `hbond_longrange` + `coop_helix` minus `compactness` would rank better | **HYPOTHESIS**, ORACLE-informed, untested (§5.4) |

---

## 1. THE CAUSAL ABLATION — one fixed ensemble, four arms, two axes

**DEMONSTRATED**, n = 126 targets, `s16/energy_ablate.py` → `s16/results/energy_report.json`.

Every arm consumes the **same 75 windows per target** — the shipped distogram filter's top-75 set,
a pure cache read with no RNG. Nothing about the ensemble changes between arms; only the readout
does. **LEGACY-ON** and the **matched random control** were both pre-registered in the module
docstring before any filtered arm's RMSD existed (drop the worst quartile, 19 of 75, by the genuine
eleven-component Legacy total at the unfitted `DEFAULT_WEIGHTS`; the control drops a uniformly
random 19 of 75, two independent `stable_rng` draws). Every AMBER row is gated **and** ungated with
the excluded count printed.

| arm | construction | mean CA-RMSD *(ORACLE)* | median |
|---|---|---|---|
| **ctrl** — neither | average → ideal-geometry projection | 3.2052 | 2.9661 |
| **amb** — AMBER only | average → ff14SB/GBn2 k = 30 | **3.1831** | 2.9588 |
| **leg** — Legacy only | Legacy filter → average → projection | 3.2294 | 3.0518 |
| **legamb** — Legacy → AMBER | Legacy filter → average → AMBER | 3.2126 | 3.0118 |
| **rnd** — matched random filter (2 draws) | random filter → average → projection | 3.2176 | 3.0029 |
| **rndamb** — random filter → AMBER | random filter → average → AMBER | 3.1888 | 2.9629 |

### 1.1 ACCURACY — CA-RMSD, ORACLE label. Negative is better.

| contrast | **gated** mean | CI 95% | median | W/L | n | ungated (excl) |
|---|---|---|---|---|---|---|
| **AMBER main effect** `amb − ctrl` | **−0.0234** | **[−0.0375, −0.0099]** | −0.0074 | 66/57 | 123 | −0.0221 (3) |
| **LEGACY main effect** `leg − ctrl` | **+0.0243** | [−0.0120, +0.0615] | +0.0034 | 59/67 | 126 | — (0) |
| **LEGACY vs its RANDOM control** `leg − rnd` | **+0.0118** | [−0.0242, +0.0492] | +0.0024 | 60/66 | 126 | — (0) |
| random filter alone `rnd − ctrl` | +0.0124 | [−0.0022, +0.0291] | +0.0015 | 57/69 | 126 | — (0) |
| LEGACY given AMBER on `legamb − amb` | +0.0326 | [−0.0039, +0.0718] | +0.0012 | 59/63 | 122 | +0.0295 (4) |
| AMBER given LEGACY on `legamb − leg` | **−0.0181** | **[−0.0343, −0.0026]** | −0.0015 | 63/62 | 125 | −0.0168 (1) |
| both vs neither `legamb − ctrl` | +0.0054 | [−0.0359, +0.0490] | −0.0038 | 66/59 | 125 | +0.0074 (1) |
| both vs random + AMBER `legamb − rndamb` | +0.0238 | [−0.0122, +0.0626] | +0.0040 | 57/64 | 121 | +0.0239 (5) |
| **INTERACTION** `(legamb−amb) − (leg−ctrl)` | **+0.0073** | [−0.0073, +0.0218] | +0.0020 | 57/65 | 122 | +0.0053 (4) |

**The median-vs-mean early warning fires on the AMBER row** (mean 3.2× the median, near-even 66/57,
CI excluding zero), so the ledger's mandated null-calibrated concentration check was run — a raw
drop-top threshold is not a test. **One verdict:** observed drop-top-10 −0.00795 against a
uniform-effect null of −0.00943 [−0.02344, +0.00477]; observed at the **58.5th percentile** of that
null ⇒ **PASS, not concentrated**. Same verdict for `legamb − leg` (60.4th percentile).

**Reading the accuracy axis:**

1. **AMBER contributes −0.023 Å, and only that.** It is 0.7% of a 3.2 Å baseline, near-even on
   wins, and — per §2.2 — larger than the aggregate but **smaller than the per-target floor of a
   null that must be zero**.
2. **LEGACY contributes NOTHING to accuracy, and if anything hurts.** `leg − ctrl` is **+0.024 Å**
   (worse), and against its own matched control it is **+0.012 [−0.024, +0.049]**, W/L 60/66. Since
   dropping 19 of 75 windows *at random* already costs +0.012 Å, **Legacy's entire measured
   accuracy effect is the cost of throwing windows away, not the identity of the windows it
   throws.** This replicates the standing ledger entry *"Legacy in GENERATION hurts the answer"*
   under a matched control it had never been given.
3. **They do not interact.** The interaction term is **+0.007 [−0.007, +0.022]** with a losing
   57/65. AMBER's effect is the same whether Legacy filtered the ensemble or not (−0.023 vs
   −0.018), and Legacy's is the same either way. **The two energies are additive and each is
   separately measurable — which is exactly what the sprint's pillar requirement asks for, and it
   is a null result, not a synergy.**
4. **Nothing survives the classical control.** `legamb − rndamb` = +0.024 [−0.012, +0.063]: the
   full Legacy → AMBER pipeline is *not* better than random-filter → AMBER.

### 1.2 VALIDITY — a separate axis, scored on the SAME structure whose RMSD is quoted

| axis | ctrl | **amb** | leg | legamb | rnd | **AMBER main effect** `amb − ctrl`, gated | improved/degraded |
|---|---|---|---|---|---|---|---|
| Ramachandran favoured ↑ | 0.7029 | **0.8736** | 0.7234 | 0.8846 | 0.6918 | **+0.1698 [+0.1398, +0.2021]** | **103 / 5** |
| Ramachandran outlier ↓ | 0.1204 | **0.0410** | 0.1102 | 0.0448 | 0.1274 | **−0.0806 [−0.1027, −0.0600]** | **68 / 9** |
| heavy clashes < 2.0 Å ↓ | 1.3968 | **0.0000** | 1.3651 | 0.0397 | 1.5000 | **−1.3740 [−1.7561, −1.0325]** | **60 / 0** |
| heavy clashes < 2.6 Å ↓ | 6.3333 | **0.0317** | 6.1825 | 0.2381 | 6.4921 | **−6.1545 [−7.3740, −5.0323]** | **112 / 0** |
| min heavy separation (Å) ↑ | 1.9629 | **2.8061** | 2.0225 | 2.7828 | 1.9536 | **+0.8320 [+0.7482, +0.9124]** | **120 / 3** |
| bond strain ↓ | **0.0000** | 0.0272 | 0.0000 | 0.0332 | 0.0000 | +0.0262 [+0.0229, +0.0298] | 0 / 123 |
| angle strain ↓ | **0.0000** | 0.0257 | 0.0000 | 0.0261 | 0.0000 | +0.0243 [+0.0229, +0.0258] | 0 / 123 |
| **cis-peptide fraction ↓** | **0.0000** | **0.0760** | 0.0000 | 0.0640 | 0.0000 | **+0.0737 [+0.0525, +0.0979]** | **0 / 44** |
| ω deviation from 180° ↓ | **0.0** | 21.88 | 0.0 | 19.73 | 0.0 | +21.45 [+17.77, +25.45] | 0 / 123 |
| rms rel. bond/angle deviation ↓ | **0.0000** | 0.0192 | 0.0000 | 0.0194 | 0.0000 | +0.0188 [+0.0167, +0.0210] | 0 / 123 |

*(`rama_allowed` is a residual of the other two and is not independently interpretable; it is in the
JSON. The projection's zeros on the bottom five rows are **exact by construction** — it builds from
ideal internal coordinates — so those rows price what AMBER PAYS, not what it fails at.)*

**Legacy's validity contribution, against its matched random control** — the only place in the
entire ablation where a Legacy contrast has a CI excluding zero:

| axis | `leg − ctrl` | `leg − rnd` | improved/degraded |
|---|---|---|---|
| Ramachandran favoured ↑ | +0.0205 [−0.0047, +0.0463] | **+0.0316 [+0.0064, +0.0580]** | 48 / 37 |
| min heavy separation ↑ | +0.0596 [−0.0024, +0.1216] | **+0.0689 [+0.0155, +0.1248]** | 62 / 64 |
| clashes < 2.0 Å ↓ | −0.0317 [−0.1984, +0.1429] | −0.1349 [−0.3095, +0.0516] | 38 / 21 |
| clashes < 2.6 Å ↓ | −0.1508 [−0.5952, +0.3016] | −0.3095 [−0.7143, +0.1111] | 63 / 45 |

> ### The causal answer, in one paragraph
>
> **AMBER contributes stereochemical repair and 0.023 Å of accuracy; Legacy contributes a small
> amount of stereochemical pre-filtering and nothing to accuracy; they do not interact; and the
> accuracy half of AMBER's contribution is inside a null that must be zero.** On validity AMBER is
> overwhelming — 103W/5L on Ramachandran, 112W/0L on 2.6 Å clashes, 120W/3L on minimum heavy
> separation — and it **pays for it** in internal-coordinate strain and, newly measured here, in
> **cis peptide bonds on 44 of 126 targets where neither the input nor the projection has any**.
> On accuracy AMBER's −0.023 Å [−0.0375, −0.0099] is real as an aggregate and unresolvable
> per-target; Legacy's +0.024 Å is indistinguishable from the cost of discarding 19 windows at
> random. **Legacy's only CI-excluding-zero effect anywhere in this ablation is +0.032 Ramachandran
> and +0.069 Å minimum separation over its matched random control — a detection effect, on the
> validity axis, exactly where the record says its defensible role is.**

### 1.3 What the Legacy filter did to the ensemble — the set-mean law, and its trap

The ledger requires set **diversity** beside any aggregate, because a method that concentrates the
set improves the proxy while degrading the output. Measured (ORACLE labels):

| filter | set diversity (mean pairwise RMSD) | ORACLE set **mean** | ORACLE set **best** | raw average RMSD |
|---|---|---|---|---|
| none (75 windows) | **2.453** | 3.551 | **2.306** | 3.050 |
| **Legacy, worst quartile dropped** | **2.307** (−6.0%) | **3.531** (better) | **2.416** (+0.110, worse) | 3.078 |
| random 1 (19 dropped) | 2.444 | 3.546 | 2.334 | 3.056 |
| random 2 | 2.436 | 3.548 | 2.348 | 3.053 |

**The Legacy filter does exactly what the set-mean law predicts and warns about**: it improves the
set MEAN (3.551 → 3.531, better than either random draw) while **destroying the set BEST**
(2.306 → 2.416, worse than either random draw) and contracting diversity 6%. The proxy it optimises
moves the right way; the thing that would actually matter moves the wrong way. This is the
mechanism behind the +0.024 Å in §1.1, measured rather than inferred.

**Why the pre-registered LEGACY-ON arm is a rank filter and not a pure clash gate**, recorded before
any RMSD was consulted: Legacy's steric term is non-zero on a **mean 10.2%** of the 75 windows
(max 44%, and **exactly zero on 20 of 126 targets**), so a pure clash gate cannot move a 75-window
coordinate average on a sixth of the instrument. The clash gate is measured instead as a
**detector**, in §3, which is the axis where it has a defensible role — and where it turns out to
be near-perfect.

### 1.4 Convergence and cost

Non-converged minimisations excluded by the gate, per arm: `ctrlamber` 3, `legamber` 1,
`rnd1amber` 3, `rnd2amber` 1. Cost measured over 126 targets under four-shard contention:
projection **5.01 s**, AMBER restrained minimisation **13.88 s mean / 13.48 s median** — consistent
with the uncontended 4.46 s / 12.57 s. **Not the 8–14 ms single point.**

---

## 2. THE FOUR FIXES

Nothing in §1, §3, §4 or §5 is quotable without these, and all four are in place.

### 2.1 The convergence gate — DEMONSTRATED, and it is now in the force field's own module

`core.amber.refine_coords` asked `LocalEnergyMinimizer` to run to convergence and never checked
that it had. The rule is **declared in `core/amber.py` itself, in the code, before it was applied
anywhere**:

> a restrained minimisation is **CONVERGED** iff its final potential energy, with the restraint
> switched off, is finite and **≤ 1000.0 kcal/mol** (`CONVERGE_MAX_KCAL`).

One threshold. No target-specific tuning. No native-derived quantity. The same 1000 kcal/mol the
codebase already used for its bond+angle strain gate. It is **reported, never silently applied**:
`refine_coords` still returns the structure and now also returns `converged` / `converge_reason`,
so the *consumer* decides. `energy_lib.gated_paired` makes it impossible to report a gated number
without its ungated twin and the excluded count.

The gate adds **no arithmetic to the minimisation**, which was checked rather than assumed:
`tests/test_amber.py` passes unchanged, and all **126** Sprint-16 control-arm targets reproduce the
Sprint-15 artefact at **max |Δ| = 0.000e+00 Å on both arms** (§2.5).

**On the pinned 126-target instrument the gate excludes 3 targets in the canonical frame** —
1D6X (1785 kcal/mol), 2NB7 (3685), 7BX2 (1493) — and **4 in the rotated frame**, where 1MF6
reaches **8.9 × 10⁸ kcal/mol**. Both figures are reported because the gate's bite is
frame-dependent, which is itself the defect.

| the k = 30 accuracy effect (AMBER − projection) | n | mean | CI 95% | W/L |
|---|---|---|---|---|
| ungated, as published | 126 | −0.02207 | [−0.03591, −0.00857] | 67/59 |
| **GATED** | **123** | **−0.02339** | [−0.03745, −0.00995] | 66/57 |

Gating moves the effect *away* from zero, so the fix carries no direction risk — and it does not
rescue the claim, for the reason in §2.2.

### 2.2 Frame invariance as a standing null — DEMONSTRATED, and it FAILS its own pre-declared band

ff14SB, GBn2, the positional restraint (`½k|x − x₀|²` with `x₀` from the same input) and every
RMSD here are rigid-invariant, so **relaxing the same structure in a rotated lab frame is zero by
construction**. Reflections are excluded (`det = +1` enforced) because a reflection is *not* a
symmetry of ff14SB.

The PASS band was declared in `energy_lib` before any Sprint-16 draw: **|mean| ≤ 0.005 Å and
max |Δ| ≤ 0.05 Å on the converged subset.**

| frame draw | n | ungated mean | sd | max\|Δ\| | **gated** mean | sd | **max\|Δ\|** | PASS |
|---|---|---|---|---|---|---|---|---|
| 1 | 126 → 122 | +0.01173 [−0.0034,+0.0383] | 0.137 | **1.5100** (1MF6) | **−0.00049** [−0.0046,+0.0035] | 0.024 | **0.1242** | **FAIL** |
| 2 | 126 → 123 | −0.00432 [−0.0089,+0.0002] | 0.026 | 0.1381 | −0.00370 [−0.0083,+0.0007] | 0.025 | **0.1381** | **FAIL** |
| 3 | 126 → 123 | +0.00021 [−0.0049,+0.0052] | 0.029 | 0.1723 | +0.00017 [−0.0052,+0.0053] | 0.029 | **0.1723** | **FAIL** |

**The gate fixes the mean of the null and does not fix its tail.** With the non-converged
minimisations removed the mean collapses from +0.0117 Å to −0.0005 Å — a genuine floor of about
0.004 Å, five times below the quoted effect — but on **converged** targets the AMBER minimiser
still moves individual structures by **0.124–0.172 Å** under a transformation that is
mathematically no transformation at all, on three independent frames.

> **The consequence, stated plainly. The MEAN effect survives its own null (the null's per-target
> sd is 0.024–0.029 Å, so the mean of 123 targets has a null standard error of ≈0.0023 Å and the
> −0.0234 Å effect is ≈10 of them). NO PER-TARGET AMBER RMSD CLAIM SURVIVES IT: the residual is
> 5–7× the entire aggregate effect on individual targets.** Any future arm that selects, weights
> or conditions on a per-target AMBER relaxation outcome is operating inside this floor.

**The check is permanent, it was run live, and it FAILED.** `s16/energy_gate.py --frame 12` and
`tests/test_amber_frame_invariance.py` run the null on a fixed 12-target slice with an
*independent* rotation draw (`stable_rng("s16energy", "frame", 1, pdb)`, 24 fresh restrained
minimisations):

| live standing null, 12 targets | n | mean | sd | max \|Δ\| | non-zero | verdict |
|---|---|---|---|---|---|---|
| ungated | 12 | +0.01314 | 0.0403 | 0.14008 | **12 / 12** | — |
| **gated** (1D6X excluded, E = 1785 kcal/mol) | 11 | **+0.01464** | 0.0419 | **0.14008** (1DEP) | **11 / 11** | **FAIL** on both criteria |

> **A comparison that is mathematically zero returns +0.0146 Å on eleven converged targets —
> LARGER IN MAGNITUDE than the −0.0234 Å accuracy effect the pipeline exists to measure.** On 123
> targets the same null averages to −0.0005 Å, so the *aggregate* effect still clears it by ≈10
> null standard errors. At n ≈ 12 it does not clear it at all.

The band has **not been loosened to make the test pass**. The test now prints the pre-declared
verdict (`PASS = False`) and asserts only a *regression guard* set from the measured floor
(|mean| ≤ 0.05, max ≤ 0.25) — it catches the minimiser getting worse, which is what a regression
check is for, and it leaves the failure visible in its own output. `1DEP` is frame-sensitive here
and again on cached frame 3 (0.135 Å), so the sensitivity is a reproducible target property.

### 2.3 Monotone conditioning — DEMONSTRATED, and winsorising at the 99th percentile is not enough

Pooled over the **17,735 binding-masked** AMBER single points in the 19 enumerated caches:

| | value |
|---|---|
| fraction > 1e6 kcal/mol | **35.8%** |
| fraction > 1e12 | **3.14%** |
| maximum | **2.65 × 10²⁰** |
| 99th percentile | **9.4 × 10¹³** |

| conditioning map | top-10 share of variance | kurtosis | order-preserving? |
|---|---|---|---|
| raw | **0.99900** | 1.71 × 10⁴ | — |
| winsorise at p99 | **0.04873** | **79.4** | **NO** (constant above the clip) |
| signed log, `sign(x−m)·log1p(|x−m|)` | 0.00450 | 1.49 | **yes**, ρ = 1.000000 |
| rank | 0.00169 | 1.80 | **yes**, ρ = 1.000000 |

Ten points out of 17,735 should carry 0.056% of the variance under any well-behaved distribution.
Raw AMBER gives them **99.9%**; winsorising at the 99th percentile still leaves them **4.9%,
87× the uniform expectation**, with a kurtosis of 79. **Per target** the winsorised figure is worse
still (mean top-10 share **0.751**), because within one target the ten worst configurations are a
larger fraction of the sample. Only the two monotone maps condition the distribution, and they
condition it *without* changing a single ordering statement — verified, not assumed
(`spearman_with_raw` = 1.000000 exactly).

**Every rank-based statistic in this document (Spearman, argmin, percentile, in-band mean) is
therefore invariant to the choice of map and immune to the tail.** The maps matter only for
moment- or spectrum-based claims, and this workstream makes none without them.

### 2.4 The binding data rule — DEMONSTRATED, enforced by construction

`energy_lib.binding_mask(z)` returns `amber_kind == 0 AND amber_idx != snap_index` and there is no
other route to an AMBER array in `s16/`. Across all 19 enumerated caches:

| | |
|---|---|
| AMBER single points cached | **43,871** |
| `kind == 0` ("unbiased") | 17,751 |
| **binding (`kind == 0` minus the forced snap index)** | **17,735** |
| files where the ORACLE snap index sits inside `kind == 0` | **16 / 19** |

Per-target n is printed beside every statistic in `s16/results/energy_contrast.json` (933 mean;
1190–1195 at n = 9, 698–700 at n = 10). **`s14/results/obj_floor.json` was not opened by any
module in `s16/`**, per the retraction.

### 2.5 A determinism result that came free — DEMONSTRATED

The ablation's control arms were recomputed from scratch rather than read from Sprint 15's cache.
**On all 126 targets** they reproduce `s15/results/phys_repl_canon0.json` at
**max |Δ| = 0.000 × 10⁰ Å on both arms** — a different sprint, a different interpreter, a different
day, `PYTHONHASHSEED` unset, and across four independent worker processes. This extends Sprint 15's
cross-process determinism check across a *code change to `core/amber.py`* and confirms the
sprint's rule that **the ≈0.08 Å multi-start false-positive floor cannot reach this pipeline** —
it has no RNG at all. It also means the ≈0.08 Å floor **does** reach the ablation's random-filter
control arms, which do have one, and those are reported over two independent draws for that reason.

---

## 3. THE REPAIR BENCHMARK — the comparison worth keeping

**DEMONSTRATED**, n = 126, `s16/energy_repair.py` → `s16/results/energy_repair.json`.

**One input structure per target, identical for every arm:** the raw all-atom coordinate average of
the shipped top-75 window set. It is genuinely broken — **5.6 heavy-atom clashes below 2.0 Å,
15.9 below 2.6 Å, 0.315 rms relative bond strain, 0.185 angle strain** — which is what makes this a
repair benchmark rather than a polishing benchmark. Pure cache read, no RNG.

**Legacy enters as the DETECTOR** (its eleven components from `core.energy.energy_components` on the
input's *real* coordinates, not on an ideal-geometry rebuild), **AMBER as the REPAIRER**, and **the
ideal-geometry projection as the classical control repairer**. Both repairers start from the same
input, and every statistic is scored on the same structure whose RMSD is quoted.

### 3.1 The repair panel

| axis | input | projection | **AMBER** | PROJECTION − input | **AMBER − input** |
|---|---|---|---|---|---|
| clashes < 2.0 Å ↓ | 5.611 | 1.397 | **0.000** | −4.214 [−6.048, −2.619] 44/16 | **−5.390 [−7.520, −3.553] 54/0** |
| clashes < 2.6 Å ↓ | 15.881 | 6.333 | **0.032** | −9.548 [−13.318, −6.246] 62/53 | **−15.252 [−19.951, −10.992] 76/0** |
| min heavy sep. (Å) ↑ | 2.096 | 1.963 | **2.806** | **−0.133 [−0.229, −0.034] (WORSE)** | **+0.692 [+0.567, +0.821] 91/32** |
| Rama favoured ↑ | 0.836 | 0.703 | **0.874** | **−0.134 [−0.163, −0.106] (WORSE), 99 targets** | **+0.038 [+0.011, +0.066]** |
| Rama outlier ↓ | 0.081 | 0.120 | **0.041** | **+0.040 [+0.015, +0.065] (WORSE)** | **−0.039 [−0.058, −0.021]** |
| bond strain ↓ | 0.315 | **0.000** | 0.027 | **−0.315 [−0.351, −0.281] 126/0** | −0.284 [−0.316, −0.251] 122/1 |
| angle strain ↓ | 0.185 | **0.000** | 0.026 | **−0.185 [−0.210, −0.160] 126/0** | −0.158 [−0.183, −0.133] 91/32 |
| rms rel. geom. dev. ↓ | 0.248 | **0.000** | 0.019 | **−0.248 [−0.281, −0.217] 126/0** | −0.225 [−0.256, −0.194] 114/9 |
| **cis-peptide fraction ↓** | **0.000** | **0.000** | **0.076** | 0.000 (no change) | **+0.074 [+0.053, +0.098] on 44 targets** |
| ω deviation from 180° ↓ | **0.0** | **0.0** | 21.88 | 0.0 | **+21.45 [+17.77, +25.45] 0/123** |

**Three findings, one of them new and inconvenient.**

1. **AMBER wins every steric and Ramachandran axis, decisively.** It removes *all* sub-2.0 Å
   clashes on 54 targets to 0, cuts sub-2.6 Å clashes 15.881 → 0.032, and gains +0.692 Å of minimum
   heavy separation.
2. **The classical control repairer makes the structure WORSE on the axes AMBER is best at.** The
   ideal-geometry projection *degrades* Ramachandran-favoured from **0.836 to 0.703** (99 targets
   worse to 9 better) and *reduces* the minimum heavy separation from 2.096 Å to 1.963 Å, while
   perfecting internal coordinates by construction. **The projection buys ideal bonds and angles by
   spending clashes and torsion plausibility.** So the repair comparison is not "AMBER vs nothing";
   it is a genuine trade between two repairers, and they trade in opposite directions.
3. **NEW DEFECT, and it counts against AMBER: restrained ff14SB/GBn2 relaxation introduces
   cis peptide bonds on 44 of 126 targets** (mean cis fraction 0.000 → 0.076), where neither the
   input nor the projection has a single one, together with a 21.9° mean ω deviation. This has not
   been reported anywhere in this project. It does not change AMBER's ranking on the clash and
   Ramachandran axes, and it means the phrase *"AMBER produces a valid structure"* is too strong:
   **AMBER produces a sterically valid structure with a peptide-bond geometry defect the incumbent
   does not have.** `k_restraint = 30` pins N/CA/C but not ω directly, which is the plausible
   mechanism; it was not tested.

### 3.2 Legacy as the DETECTOR — AUROC of each component on the INPUT structure

Label = "this target is worse than the instrument median on this axis". `hbond_longrange` and
`coop_sheet` are **identically zero on all 126 inputs** and have no detector at all — the two
largest empirical weights in the potential (3.0 and 2.0) contribute nothing to this class of
structure.

| defect | **steric** | contact (MJ) | hbond_local | coop_helix | solvation | electrostatic | aromatic | **torsion** | compactness | **total (UNFITTED)** |
|---|---|---|---|---|---|---|---|---|---|---|
| clashes < 2.0 Å | **0.989** | 0.600 | 0.785 | 0.779 | 0.492 | 0.496 | 0.467 | 0.901 | 0.362 | 0.938 |
| clashes < 2.6 Å | **0.991** | 0.565 | 0.785 | 0.754 | 0.503 | 0.457 | 0.442 | 0.902 | 0.371 | 0.916 |
| min heavy separation | **0.984** | 0.593 | 0.833 | 0.796 | 0.485 | 0.464 | 0.477 | 0.920 | 0.357 | 0.925 |
| Rama outlier | 0.797 | 0.485 | 0.883 | 0.858 | 0.556 | 0.417 | 0.466 | **0.925** | 0.521 | 0.855 |
| Rama favoured | 0.793 | 0.444 | 0.929 | 0.901 | 0.580 | 0.438 | 0.446 | **0.976** | 0.488 | 0.894 |
| bond strain | 0.905 | 0.470 | 0.915 | 0.893 | 0.552 | 0.421 | 0.454 | **0.991** | 0.446 | 0.949 |
| angle strain | 0.838 | 0.441 | 0.958 | 0.942 | 0.561 | 0.419 | 0.438 | **0.993** | 0.481 | 0.926 |
| rms rel. geom. dev. | 0.875 | 0.457 | 0.933 | 0.917 | 0.558 | 0.418 | 0.449 | **0.995** | 0.453 | 0.936 |
| D-chirality centres | 0.799 | 0.508 | 0.768 | 0.720 | 0.509 | 0.467 | 0.356 | **0.887** | 0.516 | 0.801 |
| ω deviation | 0.628 | 0.352 | 0.766 | 0.704 | 0.614 | 0.394 | 0.435 | 0.761 | 0.613 | 0.682 |

**Five readings, and the first is a caveat against my own headline.**

1. **`steric`'s ≈0.99 on the clash axes is largely TAUTOLOGICAL** — the Legacy steric term is
   computed *from* heavy-atom distances, so "the clash term detects clashes" is a consistency
   check on the implementation, not a discovery. It is reported as a confirmation that the
   detector works, not as a capability.
2. **`torsion` is the genuinely informative detector, and it is not tautological.** It is a
   Ramachandran-basin score with no knowledge of bond lengths, yet it separates bond strain at
   **AUROC 0.991**, angle strain at **0.993** and total geometric deviation at **0.995** — better
   than the steric term and better than the eleven-term total. A coordinate average whose bonds
   have collapsed lands in implausible torsion space, and Legacy sees that.
3. **`hbond_local` (0.77–0.96) and `coop_helix` (0.70–0.94) are strong secondary detectors** of
   general deformation.
4. **The MJ contact term is useless as a detector** (0.35–0.60, at or below chance on eight of ten
   axes) and **`compactness` is an ANTI-detector** on every steric axis (0.357–0.371, i.e. it
   points systematically the wrong way — the same term §5.4 shows is anti-correlated with
   correctness and whose removal moves the native from the 53.5th to the 33.6th percentile).
   `solvation`, `electrostatic` and `aromatic` sit at chance (0.42–0.61).
5. **The eleven-term total is worse than its own best single component on every one of the ten
   axes.** Combining the components at the unfitted `DEFAULT_WEIGHTS` *destroys* detection skill —
   the same shape as §4 (fitting them destroys ranking skill) and §5.4 (`only_steric` ranks better
   than the total). **Legacy is a good instrument in pieces and a poor one assembled.**

> **The comparison that survives.** On an axis where both models have a defensible role, they are
> genuinely complementary and separately measurable: **Legacy detects the defects (AUROC 0.98–0.99
> for steric, 0.93–0.99 for torsion) and AMBER removes them (all clashes, +0.69 Å separation,
> +0.038 Ramachandran-favoured against the input and +0.170 against the projection).** Neither can
> do the other's job, and neither of these claims requires either model to rank a native.

### 3.3 Cost, at the right number

`refine_coords(k = 30, steps = 0)` is a full restrained **minimisation**: **13.88 s mean / 13.48 s
median** measured here over 126 targets under four-shard contention, consistent with the
uncontended **12.57 s**. The projection is **5.01 s** measured (4.46 s uncontended). An AMBER
**single point** is 8–14 ms — **≈1,500× cheaper, and the wrong number for this arm.** Every earlier
budget-matched claim in this project used the single-point figure and was therefore wrong by three
orders of magnitude. **AMBER is the expensive repairer, by ≈2.8×**, and Legacy's detection pass is
essentially free (the eleven components on 126 structures cost seconds in total).

### 3.4 A THIRD Sprint-14/15 reporting defect, found here, and it runs AGAINST this workstream

**DEMONSTRATED, n = 126.** `core.project.lam_path` optimises torsions with an unconstrained L-BFGS
and returns them **unwrapped** — this instrument contains φ values of **+1022.72°** and ψ values of
**−1115.93°**, meaning −57.28° and −35.93°. `s14.ener_avgrefine.rama_ok` (and every caller,
including `s15.phys_repl`) compares those raw degrees against fixed Ramachandran windows, so an
unwrapped *favoured* residue is scored as an **outlier**.

| arm A (the ideal-geometry projection), Ramachandran-favoured | value |
|---|---|
| **as published** (unwrapped torsions) | **0.4658** |
| **corrected** (torsions wrapped to (−180°, 180]) | **0.7029** |
| paired difference | **+0.2371 [+0.1946, +0.2796]**, **0 W / 85 L** — never in the other direction |
| targets affected | **85 / 126** |
| AMBER arm's own figure | **unchanged, max \|Δ\| = 0.000e+00** (its torsions come from coordinates) |

> **Consequence: AMBER's Ramachandran advantage over the projection is +0.171, not +0.408 — the
> published figure overstates it by 2.39×.** Sprint 15 corrected this same statistic upward
> (0.734 → 0.466 for the baseline) after finding it was scored on the wrong structure; that
> correction was right in *kind* and its magnitude was inflated by this second, independent bug in
> the same line of code. The Sprint 15 headline **"0.466 → 0.874 Ramachandran-favoured, 116W/2L"**
> should read **"0.703 → 0.874, 103W/5L"**.
>
> **The validity claim survives this correction and is still large and one-sided** — 103 targets to
> 5, with clashes and minimum separation untouched by the bug — but it is **not as large as the
> record says**, and the record should be fixed. `s16.energy_lib.wrap_deg` fixes it for this
> workstream; `s14/ener_avgrefine.py` and `s15/phys_repl.py` are flagged, not edited.

---

## 4. LEGACY'S ELEVEN WEIGHTS, RESOLVED

**The path taken: fit them leave-fold-out, measure the ablation, and let the measurement decide.**
It decided against the fit.

`DEFAULT_WEIGHTS` is documented in the repository as *"a variance-balanced starting point, not a
fit"*; `WEIGHT_ORIGIN` marks seven of the eleven `empirical … MUST be calibrated`, and
`calibrate_weights` was never run. Comparing a never-fitted statistical potential against a force
field fitted over decades compares **training status**, not physics — so the comparison had to be
made fair before it could be made at all.

**Instrument.** The 19 exhaustively enumerated targets — 9 at n = 9 (262,144 configurations) and
10 at n = 10 (1,048,576) — **12,845,056 configurations**, each carrying all eleven Legacy
components and an exact ORACLE CA-RMSD. The spaces are *complete*, so an argmin is a **certified
optimum**, not a search result.

**Leakage.** The fit consumes ORACLE CA-RMSD as a **training label**, which is permitted only
leave-fold-out and only declared. `fit_fold(f)` trains on targets whose pinned fold ≠ f and scores
only targets whose fold = f. No held-out target's native touches its own weight vector. Every
fitted row is tagged **FITTED_LOFO** and every unfitted row **UNFITTED**; the two are never
presented as equally trained.

**Pre-registered before the first fold returned:** OLS/ridge of the eleven *per-target
standardised* components onto the *per-target standardised* ORACLE RMSD, pooled over training
targets; ridge fixed at `1e-3 · tr(XᵀX)/p`, one value, not tuned; 40,000 configurations sampled
per training target. Success = fitted beats unfitted on **both** mean ρ(E, RMSD) **and** the
certified-argmin RMSD, with a paired CI excluding zero on the argmin.

### 4.1 The result — REFUTED (my own hypothesis that fitting would rescue Legacy)

| statistic (19 targets, 5 pinned folds) | **UNFITTED** | **FITTED_LOFO** | fitted − unfitted | CI 95% | W/L |
|---|---|---|---|---|---|
| Spearman ρ(E, ORACLE RMSD) | +0.1054 | +0.0370 | **−0.0684** | [−0.2935, +0.1428] | 9/10 |
| **certified-argmin RMSD (Å)** | 4.0851 | 4.8636 | **+0.7785** | **[+0.0404, +1.4975]** | 8/11 |
| argmin − space mean (Å) | +0.0820 | +0.8605 | +0.7785 | [+0.0404, +1.4975] | 8/11 |
| ORACLE native percentile | 0.5350 | 0.5035 | −0.0315 | [−0.1949, +0.1304] | 8/11 |

Per fold, the argmin penalty is +1.373 / +0.050 / +2.660 / +0.733 / −0.140; median +0.647.
**Leave-fold-out fitting makes Legacy's certified optimum significantly WORSE**, and does not
improve its rank correlation.

### 4.2 The fitter is not broken — the generalisation is. ORACLE DIAGNOSTIC.

Fitted on *all nineteen* targets and scored on the same nineteen (an in-sample ceiling, not a
claim):

| statistic | UNFITTED | FITTED (in-sample) | Δ |
|---|---|---|---|
| ρ | +0.1054 | **+0.2130** | +0.108 |
| certified-argmin RMSD | 4.0851 | **3.6720** | −0.413 |
| argmin − space mean | +0.0820 | **−0.3311** | −0.413 |
| ORACLE native percentile | 0.5350 | **0.2848** | −0.250 |

In-sample the fit **doubles** the rank correlation, moves the certified optimum from 0.08 Å *worse*
than a random member to 0.33 Å *better*, and moves the native from the 53.5th to the 28.5th
percentile. So the eleven-component basis **does** contain a weighting that ranks — and it does not
transfer between folds.

**And the weight vectors themselves are stable**, which rules out the obvious explanation. Across
the five leave-fold-out fits the standardised coefficients move very little:

| term | fold 0 | 1 | 2 | 3 | 4 | in-sample |
|---|---|---|---|---|---|---|
| **steric** | +0.232 | +0.237 | +0.242 | +0.235 | +0.172 | **+0.227** |
| **contact (MJ)** | +0.116 | +0.097 | +0.119 | +0.051 | +0.154 | +0.104 |
| solvation | −0.080 | −0.068 | −0.066 | −0.091 | −0.031 | −0.070 |
| coop_helix | +0.057 | +0.047 | +0.050 | +0.034 | +0.037 | +0.046 |
| aromatic | −0.026 | −0.072 | −0.045 | −0.038 | −0.041 | −0.045 |
| **compactness** | **−0.069** | **+0.078** | **−0.087** | **−0.083** | **+0.122** | −0.013 |
| electrostatic | +0.001 | +0.038 | +0.007 | +0.033 | +0.071 | +0.028 |
| hbond_longrange | −0.021 | −0.019 | −0.021 | −0.019 | −0.016 | −0.019 |
| torsion, coop_sheet, hbond_local | ≈0 | ≈0 | ≈0 | ≈0 | ≈0 | ≈0 |

Nine of the eleven coefficients are essentially identical across folds. The fit is *reproducible*
and still does not transfer, so the failure is not weight-estimation noise: **the direction the
eleven components jointly point in is target-specific.** That is the same shape as the standing
ledger entry *in-band ordering is learnable but per-target* (0.986 within a target, 0.600 across),
reaching Legacy's own basis. Only `compactness` and `electrostatic` change sign between folds.

### 4.3 The decision, and the limit it puts on every Legacy claim

> **Legacy's eleven weights are CHARACTERISED AS FIXED.** The alternative — fitting them — was
> attempted under a pre-registered protocol, is reported here in full, and **loses**. Every Legacy
> claim in this document is therefore a claim about the *unfitted, variance-balanced*
> `DEFAULT_WEIGHTS` potential, and is limited accordingly:
>
> * Legacy may not be described as a *trained* ranker, and its ranking failures may not be excused
>   as "it was never fitted" — because fitting it, honestly, makes it worse.
> * Legacy versus AMBER on ranking is a comparison between an **unfitted** statistical potential
>   and a **fitted** force field, and every such row in §5 says so.
> * A future Legacy weight fit must beat +0.7785 Å on this instrument to be reportable at all.

**What this does NOT show.** It does not show that no fit of any Legacy-like potential can
generalise — only that a leave-fold-out linear fit on the eleven existing components does not, on
19 targets at n = 9–10, against RMSD as the label. n = 19 is small; the argmin CI's near bound is
+0.04 Å. A pairwise-ranking loss, a nonlinear combiner, or a different label were not tried.

---

## 5. THE THREE-WAY ENERGY CONTRAST — one instrument, identical chain lengths

**DEMONSTRATED.** All three objectives scored on the *same* 19 exhaustively enumerated spaces
(n = 9 and n = 10, k = 4), against the same ORACLE RMSD labels, with the same statistics. Because
the spaces are complete, the Legacy and distogram argmins are **certified optima**: no search, no
budget, no optimiser to blame. AMBER exists only on the cached subsample, so it is scored on the
**binding-masked** population, and Legacy and the distogram are *also* reported restricted to
exactly those positions — three arms on one identical configuration set, not three different ones.

The control is `population mean RMSD` = the expected RMSD of **a randomly chosen feasible member
of the enumerated space**, which is precisely the control Roget et al. use. A zero-information
constant α-helix is reported beside it, because on this instrument uniform random is not a
sufficient control.

| arm | population | n_pop | argmin (Å) | pop. mean | **argmin − mean** | CI 95% | median | W/L | ρ | ρ > 0 | in-band-100 | ORACLE native %ile |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Legacy/MJ (UNFITTED)** | full enumeration | 676,056 | 4.085 | 4.003 | **+0.082** | [−0.259, +0.401] | +0.168 | **5/14** | +0.105 | 74% | 3.983 | **0.535** |
| Legacy/MJ (UNFITTED) | binding-masked | 933 | 3.822 | 4.006 | −0.184 | [−0.523, +0.184] | −0.481 | 12/7 | +0.110 | 79% | 4.024 | 0.531 |
| **AMBER ff14SB/GBn2 (FITTED force field)** | binding-masked | 933 | 3.931 | 4.006 | **−0.075** | [−0.540, +0.358] | +0.079 | **7/12** | +0.148 | 68% | 3.781 | **0.407** |
| **distance-restraint objective** | binding-masked | 933 | 3.116 | 4.006 | **−0.890** | **[−1.316, −0.464]** | −0.765 | 14/5 | +0.483 | 89% | 3.355 | — |
| **distance-restraint objective** | full enumeration | 676,056 | 3.077 | 4.003 | **−0.926** | **[−1.440, −0.431]** | −0.457 | **15/4** | +0.485 | 89% | 3.141 | **0.291** |

Head to head, on the identical binding-masked set, paired and fold-aware:

| contrast | argmin (Å) | CI 95% | W/L | Δρ | CI 95% | W/L |
|---|---|---|---|---|---|---|
| distogram − Legacy | **−0.706** | [−1.274, −0.148] | 13/5 | **+0.373** | [+0.186, +0.564] | 16/3 |
| distogram − AMBER | **−0.815** | [−1.357, −0.310] | 14/5 | **+0.335** | [+0.186, +0.491] | 14/5 |
| **Legacy − AMBER** | −0.109 | [−0.666, +0.462] | 11/8 | −0.038 | [−0.263, +0.181] | 9/10 |

ORACLE native percentile against chance (0.5), paired over 19 targets:

| | mean | vs 0.5 | CI 95% | W/L |
|---|---|---|---|---|
| AMBER (binding-masked) | 0.407 | −0.093 | [−0.226, +0.039] | 10/9 |
| Legacy (binding-masked) | 0.531 | +0.031 | [−0.146, +0.198] | 8/11 |
| Legacy (full enumeration) | 0.535 | +0.035 | [−0.141, +0.201] | 8/11 |
| **distance-restraint (full enumeration)** | **0.291** | **−0.209** | **[−0.337, −0.068]** | 13/6 |

### 5.1 What the table says

1. **Neither energy ranks the native, and they are indistinguishable from each other.** Legacy's
   certified optimum is 0.082 Å *worse* than a random feasible member and loses on 14 of 19
   targets; AMBER's argmin is 0.075 Å better with a median of +0.079 and a *losing* 7/12 record;
   the head-to-head between them is +0.109 Å with a CI four times its width and ρ difference
   −0.038. **A never-fitted statistical potential and a decades-fitted all-atom force field are
   at the same place**, which is the cleanest available evidence that the failure is not a
   training-status artefact.
2. **The native sits at chance for both.** Legacy 0.535, AMBER 0.407 — neither CI excludes 0.5.
3. **The distance-restraint objective breaks the pattern at the same chain lengths, on the same
   configurations, through the same statistics.** Certified argmin **−0.926 Å [−1.440, −0.431]**
   better than random, 15/4, ρ +0.485 positive on 89% of targets, native at the **29th percentile**
   with a CI excluding chance. It beats *both* energies on the identical masked set with CIs
   excluding zero on both the argmin and ρ.

> **The pathology is a property of CONTACT-STYLE POTENTIALS, not of short-chain structure
> prediction.** At n = 9–10 an objective built from predicted pairwise distances finds structures
> that are almost a full ångström better than random, while both physics energies do not. This is
> the strongest energy result available and it needs both energies to exist as controls.

### 5.2 The zero-information control, which makes the Legacy row worse

| | mean RMSD (Å) | vs the enumerated space mean (4.003) |
|---|---|---|
| constant α-helix (nearest library state to −57°/−47° at every residue) | **3.616** | −0.387 [−0.995, +0.143], 9W/10L |
| constant extended / PPII (−120°/130°) | 4.092 | +0.089 [−1.012, +1.141], 8W/11L |
| enumerated space best (ORACLE ceiling) | 1.030 | — |

| certified optimum vs the constant α-helix | Δ (Å) | CI 95% | W/L |
|---|---|---|---|
| **Legacy (full enumeration)** | **+0.469** | [−0.006, +1.041] | 5/6 (**8 exact ties**) |
| Legacy (binding-masked) | +0.206 | [−0.307, +0.755] | 12/7 |
| AMBER (binding-masked) | +0.315 | [−0.454, +1.156] | 10/9 |
| distance-restraint (full enumeration) | **−0.539** | [−1.152, +0.017] | 10/7 |

**On 8 of the 19 targets Legacy's certified optimum *is* the constant α-helix**, to within float
error. A zero-information sequence-blind constant beats the certified optimum of an eleven-term
statistical potential by 0.469 Å. This is why "beats uniform random" proves nothing on this
instrument, and it is reported here because it makes the Legacy row *worse*, not better.

### 5.3 Prior art — this must be cited, not claimed

Roget et al., **arXiv:2606.21241**, enumerate every peptide up to length 15 across 12,446 PDB
structures and report that the minimum-cost conformation has on average a **larger** RMSD than a
randomly chosen feasible one, with the mechanism that the Miyazawa–Jernigan contact potential was
parameterised on proteins above 50 residues, so ρ(cost, RMSD) is negative at these lengths.

**The Legacy row above is their result.** It reproduces at +0.082 Å on 19 complete spaces and must
be attributed to them.

**What is ours** is (a) the same measurement for an **all-atom force field**, whose parameters have
no length provenance and which their mechanism cannot cover, arriving at the same place; (b) the
demonstration that **fitting the statistical potential's weights leave-fold-out makes it worse**
(§4), which removes the training-status escape from their result; and (c) the **counterpart arm**
— a distance-restraint objective at the *same* chain lengths on the *same* configurations that is
−0.926 Å better than random — which localises the pathology to the potential class rather than to
the problem.

### 5.4 Where inside Legacy the pathology lives — a full eleven-term ablation

**DEMONSTRATED**, on all 19 complete enumerations, 12.8 × 10⁶ configurations. Each term is removed
(`no_X`, weight → 0) and each is run alone (`only_X`), with the argmin averaged over its **full
tied set** — the ledger's recorded trap, and it bites hard here: `only_steric` has a *median of
620,109 tied minima*, so `np.argmin` on it would have read the cache's ORACLE sort order and
invented an optimum. The `only_X` rows with large tie counts are **filters**, not rankers, and are
labelled as such.

| arm | certified argmin (Å) | argmin − space mean | CI 95% | W/L | ρ | ORACLE native %ile | median ties |
|---|---|---|---|---|---|---|---|
| **full Legacy (UNFITTED)** | 4.085 | **+0.082** | [−0.26, +0.40] | 5/14 | +0.105 | 0.535 | 1 |
| no_steric | 5.748 | +1.745 | [+1.13, +2.38] | 2/17 | −0.083 | 0.537 | 1 |
| **only_steric** *(a FILTER)* | 3.870 | **−0.133** | **[−0.21, −0.06]** | **14/5** | **+0.215** | **0.176** | 620,109 |
| no_contact | 4.093 | +0.090 | [−0.25, +0.41] | 5/14 | +0.074 | 0.553 | 1 |
| **only_contact (MJ alone)** | 4.378 | **+0.375** | [−0.22, +0.96] | 9/10 | +0.161 | 0.465 | 4 |
| no_hbond_local | 4.621 | +0.618 | [+0.08, +1.20] | 4/15 | +0.101 | 0.529 | 1 |
| only_hbond_local *(filter)* | 3.846 | −0.157 | [−0.50, +0.18] | 11/8 | −0.011 | 0.237 | — |
| no_hbond_longrange | 3.968 | **−0.035** | [−0.37, +0.28] | 8/11 | +0.108 | 0.534 | 1 |
| only_hbond_longrange | 6.026 | +2.023 | [+1.40, +2.60] | 2/17 | −0.119 | **0.016** | 4 |
| no_coop_helix | 4.943 | +0.940 | [+0.42, +1.44] | 3/16 | +0.104 | 0.539 | 1 |
| only_coop_helix *(filter)* | 3.839 | −0.164 | [−0.38, +0.05] | 13/6 | +0.030 | **0.038** | 3,916 |
| no_coop_sheet | 4.085 | +0.082 | — | 5/14 | +0.105 | 0.535 | 1 |
| no_solvation | 4.251 | +0.248 | [−0.11, +0.59] | 5/14 | +0.159 | 0.525 | 1 |
| only_solvation | 5.124 | +1.121 | [+0.66, +1.60] | 2/17 | −0.143 | 0.538 | — |
| no_electrostatic | 4.093 | +0.090 | [−0.25, +0.41] | 5/14 | +0.099 | 0.528 | 1 |
| no_aromatic | 4.085 | +0.082 | [−0.26, +0.40] | 5/14 | +0.109 | 0.532 | 1 |
| no_torsion | 4.131 | +0.128 | [−0.22, +0.44] | 5/14 | +0.106 | 0.545 | 1 |
| only_torsion | 3.557 | −0.446 | [−0.94, +0.07] | 13/6 | +0.003 | 0.194 | 1 |
| **no_compactness** | 4.085 | +0.082 | [−0.26, +0.40] | 5/14 | **+0.140** | **0.336** | 1 |
| only_compactness *(filter)* | 4.651 | +0.647 | [+0.29, +1.02] | 5/14 | **−0.279** | 0.554 | 99,216 |

*(`only_coop_sheet`, `only_electrostatic` and `only_aromatic` are constant on at least one target
and are omitted rather than reported at a meaningless argmin.)*

**Four readings.**

1. **The argmin defect is DISTRIBUTED — the standing ledger entry survives at k = 4 on complete
   spaces.** No single removal repairs it: the best is `no_hbond_longrange` at −0.035 Å, whose CI
   spans zero. **Removing the MJ contact term changes essentially nothing** (+0.082 → +0.090).
   There is no term to delete.
2. **And yet MJ *alone* carries the phenomenon**, exactly as Roget et al. predict: `only_contact`
   is **+0.375 Å worse than a random feasible member**, loses 9/10, and puts the native at the
   46.5th percentile. It is the mechanism, and it is buried under ten other terms in the total.
3. **Legacy's ONE statistic with a confidence interval excluding zero is its steric term used as a
   FILTER.** Restricting to the zero-steric set is worth **−0.133 Å [−0.21, −0.06], 14W/5L**, with
   ρ = +0.215 (double the full potential's) and the native at the **17.6th percentile**. This is
   precisely and only the "clash gate" role the record already assigns it — now measured, on a
   complete space, with the tie set averaged. Removing steric costs **+1.745 Å**: it is the one
   load-bearing component.
4. **`compactness` is actively harmful and is a single-term defect for the NATIVE, if not for the
   argmin.** Alone it is *anti*-correlated with correctness (ρ = **−0.279**); deleting it raises
   the full potential's ρ from +0.105 to +0.140 and moves the native from the **53.5th to the
   33.6th** percentile — the largest single-term movement in the table. `only_hbond_longrange`
   (native %ile **0.016**) and `only_coop_helix` (**0.038**) score the native beautifully and have
   terrible argmins, which is the textbook shape of a term that favours the right answer *and* a
   great many wrong ones.

> **HYPOTHESIS (not tested here):** a Legacy restricted to `steric` as a hard filter plus
> `hbond_longrange` and `coop_helix` as native-favouring terms, with `compactness` deleted, would
> rank better than the eleven-term total. Every number needed to construct it is above, and every
> one of them is an **ORACLE DIAGNOSTIC** — the native percentile is a native-derived quantity, so
> building a potential from this table is a leak unless the construction is redone leave-fold-out.
> §4 is the warning: a leave-fold-out fit on these same components *lost*.

### 5.5 The same decomposition for AMBER — and the ledger's "interaction-only" recipe, retested

**DEMONSTRATED.** The five ff14SB/GBn2 energy groups, on the binding-masked population (mean
n = 933 per target, 19 targets), argmin averaged over its tied set.

| AMBER arm | argmin (Å) | argmin − pop. mean | CI 95% | ρ | W/L |
|---|---|---|---|---|---|
| **total** | 3.931 | −0.075 | [−0.540, +0.358] | +0.148 | 7/12 |
| **nonbonded + solvation** *(the ledger's "interaction-only")* | 3.843 | **−0.163** | [−0.607, +0.250] | +0.148 | 8/11 |
| **solvation only (GBn2)** | 3.878 | −0.128 | [−0.399, +0.182] | **+0.244** | **13/6** |
| nonbonded only | 4.150 | **+0.144** | [−0.268, +0.517] | +0.147 | 7/12 |
| bonded only (bond + angle + torsion) | 4.095 | +0.089 | [−0.217, +0.397] | **−0.012** | 11/8 |
| no_solvation | 4.194 | +0.188 | [−0.214, +0.550] | +0.147 | 7/12 |
| no_torsion | 3.788 | −0.217 | [−0.660, +0.201] | +0.148 | 9/10 |
| only_bond | 4.006 | **+0.000** | [+0.000, +0.000] | — | 0/0 |
| only_angle | 3.993 | −0.013 | [−0.054, +0.025] | — | 4/4 |

**Three readings.**

1. **`bond` is exactly constant across the whole configuration space and `angle` almost so** —
   these are single points of *ideal-geometry* builds, so the bonded terms carry no information
   by construction. AMBER's entire discriminating content on this instrument is
   `torsion + nonbonded + solvation`, and the torsion group's rank correlation is **−0.012**.
2. **The ledger's "converged interaction-only AMBER (nonbonded + solvation)" recipe does improve
   on the total** (−0.163 vs −0.075) — but the improvement is carried by **solvation**, not by the
   interaction part: `nonbonded` alone is **+0.144 Å worse than random**, and removing solvation
   costs +0.188 Å. The right statement is **"GBn2 implicit solvation is the only AMBER group with
   above-total rank skill"** (ρ +0.244, 13W/6L), not "the interaction terms are what rank".
3. **No AMBER arm's argmin CI excludes zero.** The single-term picture is the same as Legacy's:
   there is a best component, it is not the total, and it still does not rank.

---

## 6. WHAT I REFUTED, INCLUDING MY OWN

- **My own hypothesis that Legacy's ranking failure is a training-status artefact that fitting
  would repair.** A pre-registered leave-fold-out fit on the eleven components makes the certified
  optimum **+0.7785 Å [+0.0404, +1.4975] worse**, on stable weight vectors, with the same fit
  gaining substantially in-sample. **REFUTED** (§4).
- **My own pre-declared PASS band for the standing frame-invariance null.** I declared
  |mean| ≤ 0.005 Å *and* max |Δ| ≤ 0.05 Å before taking a draw. The mean passes on all three
  frames once gated; **the max fails on all three (0.124 / 0.138 / 0.172 Å)**. I have not
  loosened the band. **FAILED, reported as failed** (§2.2).
- **The expectation that the MJ contact term is the removable seat of the anti-ranking.** Deleting
  it moves the certified argmin from +0.082 to +0.090 Å — nothing. The defect is distributed, as
  the standing ledger says. MJ *alone*, however, does carry it (+0.375 Å worse than random),
  which is the Roget mechanism visible in isolation. **Partly refuted, partly confirmed** (§5.4).
- **The published arm-A Ramachandran-favoured fraction, and with it part of AMBER's own validity
  claim.** `core.project.lam_path` returns unwrapped torsions and the incumbent `rama_ok` compares
  them to fixed windows, so the incumbent scores favoured residues as outliers on **85 of 126
  targets**. Corrected, arm A is 0.4658 → **0.7029** and AMBER's Ramachandran advantage shrinks
  from +0.408 to **+0.171, a factor of 2.39**. **This runs AGAINST the workstream's own headline
  and is reported anyway** (§3.4).
- **My own assumption that Legacy's clash gate would be its usable generative role.** Its steric
  term fires on a mean 10.2% of pool windows and on *zero* windows for 20 of 126 targets, so it
  cannot move a coordinate average; and the rank filter that replaced it is
  **indistinguishable from a matched random filter** on accuracy (§1.1, §1.3).
- **The expectation that Legacy's components would combine.** The eleven-term total at
  `DEFAULT_WEIGHTS` is worse than its own best single component on **all ten** detection axes
  (§3.2) and on the certified argmin (§5.4), and a leave-fold-out refit is worse still (§4).
  **Assembling this potential destroys the skill its pieces have.**

## 7. WHAT REMAINS OPEN

- **The frame-null tail.** Gating fixes the null's mean and not its 0.12–0.17 Å per-target
  residual. Where that residual comes from inside `LocalEnergyMinimizer` is not diagnosed here.
  Until it is, no per-target AMBER RMSD claim is resolvable.
- **A non-linear or ranking-loss Legacy fit.** §4 rules out a leave-fold-out *linear* fit against
  an RMSD label on 19 targets. It does not rule out a pairwise-ranking objective, a per-target
  adaptive weighting (which §4.2's stability result points at), or a larger fitting set.
- **A Legacy rebuilt from §5.4** — steric as a hard filter, `hbond_longrange` and `coop_helix`
  retained, `compactness` deleted — is an ORACLE-informed construction and would have to be
  redone leave-fold-out before it could be quoted.
- **Two of Legacy's eleven components are identically zero on every one of the 126 raw coordinate
  averages** (`hbond_longrange`, `coop_sheet` — the two largest empirical weights, 3.0 and 2.0).
  Whether that is a property of coordinate averaging or of 9–16-residue peptides generally is not
  established here.
- **Why restrained relaxation introduces cis peptide bonds on 44 of 126 targets.** `k = 30` pins
  N/CA/C positionally and does not restrain ω directly; that is a plausible mechanism and it was
  not tested. A cis-ω penalty, or restraining ω, is the obvious next experiment and is cheap.
- **Whether AMBER's repair advantage over the projection survives a repairer that does both** —
  ideal internal coordinates *and* clash removal. Neither existing operator achieves both, and the
  panel in §3.1 shows exactly which axes each one spends.
- **Nothing in this workstream has been near the 60-target benchmark**, which stays sealed.

---

## 8. WHAT MUST CHANGE IN THE SHARED DOCUMENTS

Flagged, not edited — these are the coordinator's files.

1. **`s15/LEGACY_VS_AMBER.md` §5.2** and **`s15/phys_FINDINGS.md` §5.3** quote arm A's
   Ramachandran-favoured fraction as **0.466** and the gain as **0.466 → 0.874, 116W/2L**. Both are
   **unwrapped-torsion artefacts**. The corrected figures are **0.703 → 0.874, 103W/5L**, a gain of
   **+0.171 rather than +0.408** (§3.4). The clash and minimum-separation rows are unaffected.
2. **`s15/LEGACY_VS_AMBER.md` §5** describes AMBER's product as *"clash-free"* geometry. It should
   add that the relaxed structure carries **cis peptide bonds on 44/126 targets and a 21.9° mean ω
   deviation**, neither of which the projection has (§3.1).
3. **`s15/LEGACY_VS_AMBER.md` §5 table** lists Legacy's surviving role as *"clash gate, and nothing
   else"*. §1.3 shows the clash gate cannot act on a 75-window average (10.2% of windows flagged,
   zero on 20 targets); §3.2 shows the defensible role is **detection**, at AUROC 0.98–0.99
   (`steric`) and 0.93–0.99 (`torsion`), on structures rather than on pool members.
4. **`s15/LEGACY_VS_AMBER.md` §1** says Legacy's weights *"were never fitted"* and leaves the
   implication that fitting them would help. §4 fits them leave-fold-out and they get **+0.7785 Å
   worse**. The document should record that the escape has been closed by measurement.
5. **The ledger entry "Converged interaction-only Amber is the right way to score Amber"** should
   be narrowed: on 19 complete enumerations the gain over the total is carried by **GBn2 solvation**
   (ρ +0.244, 13W/6L), while `nonbonded` alone is **+0.144 Å worse than random** (§5.5).
6. **`s14/ener_avgrefine.rama_ok` and `s15/phys_repl`** need `wrap_deg` (or an equivalent) before
   the Ramachandran window test. `s16.energy_lib.wrap_deg` is the fix.

---

## REPRODUCTION

```
python -m s12.instrument                      # pinned constants, 126 targets
python -m s16.energy_gate                     # the four fixes, from cached artefacts
python -m s16.energy_gate --frame 12          # the STANDING null, live (~5 min of AMBER)
pytest tests/test_amber_frame_invariance.py   # the permanent regression check
python -m s16.energy_weights                  # Legacy's weights, leave-fold-out (~10 min)
python -m s16.energy_contrast                 # the three-way contrast (~5 min)
python -m s16.energy_ablate --shard 0/4       # the causal ablation, 4 cooperating workers
python -m s16.energy_ablate --shard 1/4       # ~1 h wall clock in parallel on this box
python -m s16.energy_ablate --shard 2/4
python -m s16.energy_ablate --shard 3/4
python -m s16.energy_report                   # the ablation's read-out
python -m s16.energy_repair                   # the repair benchmark
```

Auxiliary analyses reported in §5.4, §5.2 and §5.5 are inline scripts over the cached enumerations;
their outputs are `s16/results/energy_termablate.json`, `energy_constant_controls.json`,
`energy_vs_helix.json` and `energy_amberterms.json`.

Environment: `OMP_NUM_THREADS = MKL_NUM_THREADS = OPENBLAS_NUM_THREADS = 2`, OpenMM CPU platform
at `Threads = 1` with `DeterministicForces`, **checkpointed every target**. Workers cooperate
through `_done_everywhere()` (each re-reads every shard file before each target) so extra workers
can be added mid-run; `load_all()` de-duplicates by pdb.

**Two operational failures worth recording**, because both cost real work: (1) a mid-run
`MemoryError` from `core.amber.memory_guard` at a 92% physical-memory ceiling killed a shard and
cost eight targets under the original every-tenth-target checkpoint — hence per-target
checkpointing; (2) an *unbounded* CPU hold at a 96% threshold deadlocked both workers indefinitely
while sibling agents sat at 100%, producing zero progress for ~25 minutes. The hold is now
**bounded** (yield up to 30 s per target, then proceed) with a separate memory hold, which is the
gate that actually binds. All seeding via `s15.seed.stable_rng`; **no `hash()` anywhere**.

**Determinism.** The control arms reproduce `s15/results/phys_repl_canon0.json` at
**max |Δ| = 0.000e+00 Å**, across sprints, interpreters and a code change to `core/amber.py`.
The only RNG in this workstream is the random-filter control (`stable_rng("s16energy",
"rndfilter", draw, pdb)`), reported over two independent draws for exactly that reason.
