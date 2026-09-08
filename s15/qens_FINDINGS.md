# Sprint 15 — QENS workstream findings

Agent: QENS. One question: **is the project's only positive quantum result real, or is it noise?**

Tiering: **DEMONSTRATED** / **ORACLE DIAGNOSTIC** / **HYPOTHESIS** / **REFUTED**.

Instrument reproduced **at start and at end** (`python -m s12.instrument`, a cache read per the
Phase 0 audit, not a pipeline reproduction): `shipped 3.4540004952559396, pool_best
1.7108244199364904, top75_best 2.3061526409453816, synthesis_fit 3.2040761603809194,
n_zero_recall 18`.

**Status: COMPLETE.** 1,672 new paired cells (R1 456, R3 760, R4 456), plus the 108-cell
bit-level reproduction of the original and 135 enrichment-control cells. Every module
checkpoints after every cell and resumes by re-running the same command.

---

## 0. VERDICT

> **CONFIRMED as an effect, REFUTED as a quantum result.** The ensemble win replicates at
> power — **−0.385 Å [−0.748, −0.035] at alpha = 1 and −0.295 Å [−0.495, −0.069] at
> alpha = 0.25**, 19 targets × 8 seeds, target-level unit, concentration DIFFUSE, and it is
> *stronger* on the ten targets that were never used to find it — but it is a pure **location**
> effect belonging to the **objective**, not to the circuit: a one-line classical Boltzmann
> reweighting of the same untrained circuit by the same objective at matched entropy, using
> **1/200th** of the objective budget, **beats the CVaR-VQE by +0.25 to +0.48 Å on every
> readout at every alpha**, classical simulated annealing at 1/400th of the budget matches or
> beats it, and the per-target sign of the "win" is set by how good the objective's own
> certified optimum is on that target (**rho = +0.52 to +0.69**).

Three subsidiary verdicts, all DEMONSTRATED:

* **The alpha non-monotonicity is REFUTED.** Mapped at ten alphas × 19 targets × 4 seeds the
  response is smooth and monotone decreasing in diversity, significant at nine of ten alphas.
  The reported "win at 1 and 0.05 but not at 0.25" was two noisy three-seed cells.
* **The original claim did not survive its own data at the correct unit of analysis** — 27
  (target, seed) cells were treated as 27 independent units; at the target unit (n = 9) its
  alpha = 1 headline is NULL and its alpha = 0.05 cell fails the concentration check. The
  replication rescues the effect that the unit correction had removed, and then the classical
  controls take it away again for a different reason.
* **The QRESTRAINT sub-2 Å enrichment cell replicates exactly (5.20/6.16/0.98/6.70/18.64×) and
  collapses to the same explanation**: classical simulated annealing at the identical
  8,192-evaluation budget is statistically indistinguishable on all four native-free
  objectives and significantly better on the ORACLE one; an entropy-matched classical use of
  the same objective doubles it.

**Nothing in this workstream clears the programme's declared bar.** No arm beats an untrained
circuit at matched budget in a way a cheaper classical method does not beat more.

---

## 1. THE CLAIM UNDER TEST, STATED EXACTLY

From `s15/qgeom_FINDINGS.md` §E2, relayed as `coord_FINDINGS.md` Q3.1:

> Consumed as an **ensemble with no ranker anywhere**, a CVaR-VQE beats best-of-N from its own
> untrained start by **0.36–0.57 Å** on a coordinate-average readout, CIs excluding zero, at
> alpha = 1 and alpha = 0.05. n = 27 (9 targets × 3 seeds). Flagged **LOW POWER**
> (mean/sd = −0.40); replicate before building on it.

Its recorded numbers, for reference:

| readout | alpha | VQE | control | diff | CI95 | W/L | verdict |
|---|---|---|---|---|---|---|---|
| rand-5 coordinate average | 1.00 | 2.722 | 3.289 | **−0.567** | [−1.080, −0.068] | 18/9 | SIG |
| rand-5 coordinate average | 0.05 | 2.930 | 3.289 | **−0.359** | [−0.708, −0.033] | 16/11 | SIG |
| rand-75 coordinate average | 1.00 | 2.673 | 3.110 | **−0.437** | [−0.836, −0.032] | 15/12 | SIG |

---

## 2. VERIFICATION BEFORE ANY CLAIM — `s15/qens_verify.py`, `s15/qens_repl.py R0`

| check | result |
|---|---|
| instrument's four pinned constants | all four reproduce to 1e-9 |
| `qens_lib.Struct` vs `qgeom_cvar.Struct` on the nine n=9 targets | `full` index vector **identical**; ORACLE rmsd column max diff **0.0** |
| precomputed CA table vs the stored ORACLE rmsd column, all 19 targets | max err **2.8e-07 – 6.7e-07 Å** (float32 storage of the cache) |
| `qens_lib.hamil_sub` vs `uniformise(combine(*tabulate(pdb), 0.25)[full])`, nine targets | **rank-identical, max\|dE\| = 0.0, rho = 1.000000000** on all nine |
| CVaR-VQE repeated from the same seed | max \|Δp_final\| **0.0** |
| **R0: bit-level reproduction of all 108 original E2 cells** | **max absolute difference 0.000e+00 over every readout and every cell** |

**DEMONSTRATED.** The original result is exactly reproducible and this workstream's machinery
is the original machinery. Nothing that follows can be a re-implementation artefact.

### 2b. A DEFECT IN THE ORIGINAL READOUT, found here — DEMONSTRATED

`qgeom_ens.unranked_readouts` computes

    out["set_coordavg_rmsd"] = st.coord_avg_rmsd(u[: min(len(u), 400)])

where `u = np.unique(idx)`. `np.unique` returns configuration indices **sorted**, so `u[:400]`
is the 400 **lowest-indexed** configurations — a systematic slice of the big-endian register
(the most significant sub-register residue held in its low states), not a random subsample.
It truncates arms with more than 400 distinct draws and leaves narrower arms untouched, so it
**bites unequally across arms** and is not a valid cross-arm comparison. The affected number in
the original table is `whole-set coordinate average` (recorded −0.414, NULL). The `rand-m`
readouts, including the headline `rand-5`, use a proper `rng.choice` and are unaffected.

This workstream keeps `set_coordavg_rmsd` under the original definition (so R0 reproduces
exactly) and adds `setrand_coordavg_rmsd`, the unbiased version, which is the one it reads.

---

## 3. THE UNIT OF ANALYSIS — the original claim does not survive its own data

**DEMONSTRATED**, from `qgeom_ens.json` alone, before a single new cell was run.

The original treated **27 (target, seed) cells as 27 independent units**. Three seeds on the
same target are repeated measurements of the same instance, not independent replicates of a
target effect. Recomputed with seeds averaged first (n = 9 targets), on exactly the same data:

| readout | alpha | **cell-level (the original convention, n=27)** | **target-level (n=9)** |
|---|---|---|---|
| rand-5 coordavg | 1.00 | −0.567 [−1.080, −0.068] 18/9 **SIG** | −0.567 [−1.234, **+0.085**] 8/1 **NULL** |
| rand-5 coordavg | 0.25 | −0.138 [−0.543, +0.238] NULL | −0.138 [−0.625, +0.317] NULL |
| rand-5 coordavg | 0.05 | −0.359 [−0.708, −0.033] 16/11 **SIG** | −0.359 [−0.797, −0.008] 7/2 SIG, **CONCENTRATED (FAIL)** |
| rand-75 coordavg | 1.00 | −0.437 [−0.836, −0.032] 15/12 **SIG** | −0.437 [−1.056, **+0.160**] 6/3 **NULL** |
| rand-75 coordavg | 0.05 | −0.126 [−0.310, +0.060] NULL | −0.126 [−0.305, +0.071] NULL |
| set mean | 1.00 | −0.847 [−1.240, −0.450] 22/5 SIG | −0.847 [−1.430, −0.275] 8/1 **SIG** |
| set best | 1.00 | +1.046 [+0.725, +1.395] 1/24 SIG | +1.046 [+0.605, +1.563] 0/9 **SIG** |

**Three readings.**

1. **At the honest unit of analysis the alpha = 1 headline is NULL**, on both coordinate-average
   readouts, with intervals crossing zero. The point estimate is unchanged (it must be — the
   design is balanced), but the interval widens by roughly √3 and swallows the effect.
2. **The alpha = 0.05 cell survives the interval and fails the null-calibrated concentration
   check at the target level** (p_share < 0.05): it is carried by two or three of the nine
   targets. The brief's rule says such a result must not be quoted as a general effect.
3. **The two rows that survive at both units are the ones that were never the positive claim** —
   the set MEAN (VQE better by 0.847) and the set BEST (VQE **worse** by 1.046). Those are the
   concentration trade, not an ensemble result.

*This is not yet a refutation.* n = 9 targets has very little power in either direction, which
is exactly why the replication below runs 19 targets and 12 seeds. It does establish that the
original claim was significant only under a convention that overstates n by a factor of three.

---

## 4. THE REPLICATION AT POWER — design

`s15/qens_repl.py`. **Original: 9 targets × 3 seeds = 27 cells per alpha. Here: 19 targets ×
12 seeds = 228 cells per alpha, and 19 independent targets against 9.**

* **19 enumerated targets.** The nine n=9 (`s13/results/qarch_enum_*.npz`) plus the ten n=10
  (`s14/cache/obj_enum_*.npz`). Every target uses the same 12-qubit sub-register (residues
  1–6), so the instance size is identical across all 19 and only the target changes.
* **Ansatz, optimiser, budget, objective: unchanged** — `ring` L=2, 200 Adam iterations at
  lr 0.10, exact statevector gradient, 2,048 readout draws, native-free `hamil` at w=0.25.
* **Six arms**, all consumed with no ranker anywhere:

| arm | what it is | objective evaluations | budget class |
|---|---|---|---|
| `vqe` | the CVaR-VQE's final distribution | **819,200** | — |
| `untrained` | the same circuit at `random_theta(seed)` — the incumbent control | 0 | A |
| `uniform` | uniform random over the sub-register — **WEAK**, see below | 0 | A |
| `anneal` | classical SA on the same objective | 2,048 | A (equal draws) |
| `anneal_cost` | classical SA, converged tail of the trajectory | **819,200** | B (equal evaluations) |
| `tilt` | classical Boltzmann tilt of `p_init`, temperature solved so its entropy equals the VQE's | 4,096 | B, diagnostic |

* **BUDGET ACCOUNTING, stated because the original convention flatters the VQE.** An exact
  statevector gradient reads **all 4,096 objective values on every one of the 200 iterations**,
  i.e. **819,200 objective evaluations**, while the original control was charged **2,048
  draws**. That is a **400× advantage** and it is the convention the positive result was
  measured under. Class A reproduces it; class B charges the classical arms the VQE's real
  count. Both are reported.
* **`uniform` proves nothing.** The brief records that a zero-information constant α-helix
  beats uniform random by 0.457 Å. It is reported for calibration only.
* **Matched-diversity readout `mdN_coordavg`**: every arm's DISTINCT configuration set is
  subsampled to a common size `D = min(distinct(vqe), distinct(untrained))` and `N` members are
  drawn uniformly from that. Equal diversity, equal budget, so only *which* configurations each
  distribution favours can move it.
* **Diversity beside every RMSD**: distinct-configuration count, draw entropy, and mean pairwise
  CA-RMSD, with a `COLLAPSED` flag. The set-mean law is out of domain on a collapsed set.

ORACLE labelling: every RMSD, every `set_*` / `rand*` / `md*` readout and every `Ep_rmsd` is
computed from native coordinates **post hoc** and is **ORACLE**. No native quantity enters any
objective, circuit, temperature solve or sampling step.

### 4a. THE FALSE-POSITIVE FLOOR, applied to every row

The RESTRAINT workstream measured an empirical false-positive floor for single-draw paired
comparisons in this machinery: a comparison that is **zero by construction** (maximum
likelihood against a constant-width Gaussian *is* least squares) returned
**+0.081 [+0.014, +0.169]** on one start draw and −0.003 on another, with absolute arm means
varying at **sd 0.132 Å across four draws**. A 95% interval excluded zero on a null effect,
from the start draw alone, and the null-calibrated concentration test flagged that row — and
only that row — as FAIL.

**Every paired mean in this workstream is therefore read against a floor of ≈ 0.08 Å.** Rows
at or below it are marked `<=FLOOR` and are reported as **UNRESOLVED**, whatever their
interval says. The cell under replication is 0.36–0.57 Å, comfortably above the floor, so the
floor does not threaten the original claim; it does bind several rows of the decomposition.
The structural defence is the one this design already uses: **eight independent start draws
per target**, averaged before the target-level test, rather than one.

### 4c. RESOURCE DISCIPLINE, recorded because it changed the design

The box ran at 91–100% CPU with five to eight other agent processes throughout. Two design
parameters were cut for resource reasons, both **before** any result from the affected cells
was read, and both applied uniformly:

* the equal-evaluation SA control (`anneal_cost`, 819,200 single-move evaluations, ~60% of
  the module's runtime) runs on **seeds 0–3 only** — 19 targets × 4 seeds = 76 paired cells;
* R1's seed count was reduced from 12 to **8**, and R3's from 8 to **4**, after five targets.
  The discarded cells were seeds 8–11, chosen by index and not by outcome; the R1 checkpoint
  was truncated so the panel is balanced.

At the coordinator's instruction the workstream was reduced to **one concurrent job** partway
through; `qens_enrich` was paused (it checkpoints after every cell) and resumed afterwards.

---

## 4b. THE SECOND POSITIVE CELL — the classical control for the QRESTRAINT enrichment

`s15/qens_enrich.py`. Added at the coordinator's request mid-workstream.

**The claim under test** (`s15/qrestraint_FINDINGS.md` §7): a trained CVaR-VQE state puts
**5.2× / 6.2× / 6.7×** more probability mass below 2.0 Å than its own untrained
initialisation on three native-free objectives (**18.6×** on an ORACLE objective), with a
clean internal control — the weakest objective (`E_ls_pool`) gives **0.98×**, no enrichment
at all. Nine n=9 targets, 3 seeds, alpha = 0.25, a hard budget of **8,192** objective
evaluations enforced by `s14.vqe_lib.Counter`, and the whole 2^18 register enumerated so the
probabilities are computed rather than sampled.

**The gap, in the finding agent's own words:** *"whether a classical sampler at matched budget
reaches the same enrichment is untested here. Without it, this is a measurement of what the
VQE does, not a claim that only a VQE does it."*

Arms, every one scored with the identical statistic — the probability mass its own
distribution places below 2.0 / 2.5 / 3.0 Å — with entropy and effective support printed
beside it, because an enrichment bought by collapsing the distribution is not an enrichment:

| arm | what it is | objective evaluations |
|---|---|---|
| `vqe` | the trained state, **re-run here per seed** rather than read as a 3-seed mean | **8,192** |
| `uniform` | the space's own base rate | 0 |
| `untrained` | the ansatz at `theta0` — the VQE's own control | 0 |
| `untr_matchdiv` | draw from `theta0` until the trained state's **effective support** distinct configurations have been seen, weight them uniformly — the coordinator's literal request | 0 |
| `anneal` | classical SA at the VQE's **exact** 8,192-evaluation budget; distribution = its empirical visit frequency | **8,192 — MATCHED** |
| `anneal_tail` | the same annealer, converged tail only — the analogue of reading the VQE's *final* distribution rather than its whole history | **8,192 — MATCHED** |
| `topM` | uniform over the M lowest-objective configurations, M = the trained state's effective support | 262,144 |
| `boltz` | `p ∝ exp(−E/T)`, T solved so the entropy equals the trained state's | 262,144 |
| `boltz_p0` | `p ∝ p0·exp(−E/T)` at the same matched entropy — the untrained circuit reweighted by the same objective | 262,144 |

`topM`, `boltz` and `boltz_p0` read the whole tabulated objective, which only an enumerated
instrument allows; they are **diagnostics of what the objective makes available**, not
budget-matched competitors. `anneal` and `anneal_tail` **are** budget-matched and they are the
arms that decide the question.

**A reporting correction applied here.** The recorded 5.2–6.7× is a **ratio of means over
targets**, which one target with large absolute near-native mass can dominate. This workstream
prints the **per-target ratio distribution** (median, geometric mean, and the fraction of
targets above 1×) beside it.

---

## 5. R1 — THE REPLICATION AT POWER. The effect is REAL and it is NOT the VQE's.

**DEMONSTRATED.** 456 cells (19 targets × 8 seeds × 3 alphas). Negative = the VQE ensemble is
better than the untrained circuit. `rand-5 coordinate average` is the original headline
readout; every number is ORACLE post-hoc.

| readout | alpha | **cell (n=152)** | **target (n=19, primary)** | W/L | concentration |
|---|---|---|---|---|---|
| **rand-5 coordavg** | 1.00 | −0.385 [−0.557, −0.218] SIG | **−0.385 [−0.748, −0.035] SIG** | 14/5 | DIFFUSE (PASS) |
| **rand-5 coordavg** | 0.25 | −0.295 [−0.441, −0.154] SIG | **−0.295 [−0.495, −0.069] SIG** | 15/4 | DIFFUSE (PASS) |
| **rand-5 coordavg** | 0.05 | −0.161 [−0.291, −0.031] SIG *(conc. FAIL)* | −0.161 [−0.328, **+0.008**] NULL | 13/6 | DIFFUSE (PASS) |
| rand-20 coordavg | 1.00 | −0.326 SIG | −0.326 [−0.691, +0.014] NULL | 13/6 | PASS |
| rand-20 coordavg | 0.25 | −0.216 SIG | **−0.216 [−0.385, −0.043] SIG** | 15/4 | PASS |
| rand-20 coordavg | 0.05 | −0.169 SIG | **−0.169 [−0.317, −0.017] SIG** | 13/6 | PASS |
| rand-75 coordavg | 1.00 | −0.320 SIG | −0.320 [−0.672, +0.015] NULL | 13/6 | PASS |
| rand-75 coordavg | 0.25 | −0.146 SIG | **−0.146 [−0.274, −0.032] SIG** | 12/7 | PASS |
| rand-75 coordavg | 0.05 | −0.122 SIG | **−0.122 [−0.202, −0.041] SIG** | 16/3 | PASS |
| whole-set coordavg (unbiased) | 1.00 | −0.326 SIG | −0.326 [−0.712, +0.032] NULL | 13/6 | PASS |
| **drawn-set MEAN** | 1.00 | −0.836 SIG | **−0.836 [−1.146, −0.531] SIG** | 17/2 | PASS |
| **drawn-set BEST** | 1.00 | +1.201 SIG | **+1.201 [+0.922, +1.513] SIG (WORSE)** | 0/19 | PASS |

**5a. THE HEADLINE REPLICATES, AND IT SURVIVES THE UNIT-OF-ANALYSIS CORRECTION.** At the
target level with 19 targets and 8 seeds the alpha = 1 win is **−0.385 [−0.748, −0.035]**,
14/5, concentration DIFFUSE — where the original's own data gave a NULL at that unit. The
alpha = 0.25 cell, which the original recorded as nothing, is **−0.295 [−0.495, −0.069]**,
15/4. Both are far above the 0.08 Å false-positive floor. **The original point estimates were
right and its interval was over-narrow at the cell level and over-wide at the target level;
at power the effect is a little smaller (−0.385 against −0.567) and much better resolved.**

**5b. IT IS NOT A PROPERTY OF THE ORIGINAL NINE TARGETS.** Split by target set, at the target
unit on rand-5:

| subset | alpha=1.0 | alpha=0.25 | alpha=0.05 |
|---|---|---|---|
| the ORIGINAL nine | −0.354 [−1.044, +0.329] NULL | −0.226 [−0.583, +0.199] NULL | −0.110 NULL |
| **the TEN NEW targets** | **−0.413 [−0.682, −0.193] SIG** 9/1 | **−0.357 [−0.533, −0.189] SIG** 8/2 | **−0.206 [−0.382, −0.046] SIG** 7/3 |

The effect is **larger and cleaner on the ten targets that were never used to find it**. That
is the strongest single argument that it is not a selection artefact.

**5c. SEED VARIANCE IS LARGE AND THE ORIGINAL'S THREE SEEDS WERE LUCKY-ISH.** Per-seed mean
difference at alpha = 1 (n = 19 targets each): **−0.123, −0.464, −0.431, −0.460, −0.442,
−0.426, −0.650, −0.088** — sd ≈ 0.19 Å across start draws, entirely consistent with the
RESTRAINT workstream's sd 0.132 Å figure. Two of eight seeds land at or near the false-positive
floor. **A three-seed study could have reported anything from −0.09 to −0.65.**

**5d. THE DIVERSITY COLUMN, mandatory beside every RMSD.**

| alpha | arm | distinct / 2048 | median | entropy (bits) | mean pairwise RMSD | cells collapsed (≤5 distinct) |
|---|---|---|---|---|---|---|
| 1.00 | **vqe** | **2.3** | **2** | 0.640 | 1.733 | **0.97 — COLLAPSED** |
| 1.00 | untrained | 471.5 | 451 | 7.543 | 3.471 | 0.00 |
| 1.00 | tilt (classical) | **5.8** | 4 | 0.636 | 1.592 | **0.62 — COLLAPSED** |
| 1.00 | anneal | 1025.6 | 1026 | 9.629 | 3.395 | 0.00 |
| 1.00 | anneal_cost | 439.6 | 442 | 8.163 | 3.292 | 0.00 |
| 0.25 | **vqe** | 128.5 | 131 | 4.088 | 3.407 | 0.00 |
| 0.25 | tilt | 102.0 | 98 | 4.120 | 2.997 | 0.00 |
| 0.05 | **vqe** | 250.3 | 250 | 5.774 | 3.459 | 0.00 |
| 0.05 | tilt | 240.3 | 239 | 5.771 | 3.260 | 0.00 |

> **At alpha = 1 the VQE arm is COLLAPSED in 97% of cells** (2.3 distinct structures in 2,048
> draws) and its "5-member coordinate average" is one structure. **The alpha = 1 row is not an
> ensemble result and must never be quoted as one.** At alpha = 0.25 and 0.05 the ensemble is
> genuinely diverse (mean pairwise RMSD 3.41 and 3.46 against the control's 3.47) and those
> rows *are* ensemble results.

---

## 6. THE DECOMPOSITION — what carries it: LOCATION, and nothing else

**DEMONSTRATED.**

**6a. It is a pure LOCATION effect.** The expected RMSD under each distribution, computed
exactly with no sampling noise (ORACLE post-hoc):

| alpha | E_p[RMSD] VQE | E_p[RMSD] untrained | diff | W/L |
|---|---|---|---|---|
| 1.00 | 3.083 | 3.918 | **−0.834 [−0.979, −0.690]** | 125/27 |
| 0.25 | 3.444 | 3.918 | **−0.474 [−0.558, −0.396]** | 127/25 |
| 0.05 | 3.627 | 3.918 | **−0.291 [−0.352, −0.232]** | 117/32 |

**6b. The coordinate-average effect is entirely the set-mean effect, at the coefficient the
recorded operator law predicts.** OLS of the paired coordinate-average difference on the
paired set-mean and set-best differences, cell level:

| alpha | fitted | R² |
|---|---|---|
| 1.00 | `d_coordavg = +0.285 + 0.943·d_set_mean + 0.099·d_set_best` | **0.747** |
| 0.25 | `d_coordavg = +0.185 + 1.008·d_set_mean − 0.035·d_set_best` | 0.325 |
| 0.05 | `d_coordavg = +0.201 + 1.204·d_set_mean − 0.157·d_set_best` | 0.281 |

The recorded law is `d_out = 1.16·d_set_mean + 0.04·d_set_best`. Measured here the set-mean
coefficient is **0.94 to 1.20**. **There is no residual "shape" term: everything the ensemble
readout sees is where the distribution sits, not how it is shaped.**

**6c. AT MATCHED DIVERSITY THE EFFECT DISAPPEARS wherever the arms are genuinely diverse.**
`mdN` subsamples every arm's distinct set to the common size `D = min(distinct(vqe),
distinct(untrained))` and averages `N` members drawn **uniformly** from it — equal diversity,
equal budget, so only *which* configurations each distribution favours can move it.

| readout | alpha | target (n=19) | verdict |
|---|---|---|---|
| md-5 coordavg | 1.00 | **−0.760 [−1.117, −0.406]** 17/2 SIG | but D ≈ 2, see below |
| md-5 coordavg | 0.25 | −0.069 [−0.161, +0.032] NULL | **≤ FLOOR — UNRESOLVED** |
| md-5 coordavg | 0.05 | −0.025 [−0.123, +0.078] NULL | **≤ FLOOR — UNRESOLVED** |
| md-20 coordavg | 0.25 | −0.051 NULL | **≤ FLOOR — UNRESOLVED** |
| md-20 coordavg | 0.05 | −0.016 NULL | **≤ FLOOR — UNRESOLVED** |
| md-75 coordavg | 0.25 | −0.126 [−0.193, −0.061] SIG | above floor |
| md-75 coordavg | 0.05 | −0.053 [−0.115, +0.005] NULL | **≤ FLOOR — UNRESOLVED** |

**Read this correctly.** At alpha = 1 the match forces `D ≈ 2`, so the "matched" comparison is
*the VQE's single converged structure against the average of two random draws from the
untrained circuit* — the match makes the control worse, and the +0.76 there is a consequence of
the matching, not evidence of shape. **At alpha = 0.25 and 0.05, where both arms are genuinely
diverse and the match is meaningful, re-weighting each arm's own distinct configurations
uniformly removes the effect on four of five rows and pushes it below the false-positive
floor.** The VQE's advantage lives in **how much probability mass it puts on its good
configurations**, not in which configurations it can reach.

---

## 7. THE CLASSICAL CONTROLS — and this is where the claim dies

**DEMONSTRATED.** Target unit, n = 19. **POSITIVE = the VQE is WORSE.**

| comparison | readout | alpha=1.0 | alpha=0.25 | alpha=0.05 |
|---|---|---|---|---|
| **vqe − tilt** (classical Boltzmann reweight of the untrained circuit at matched entropy) | rand-5 | **+0.398 [+0.061,+0.735] SIG** 6/13 | **+0.304 [+0.128,+0.487] SIG** 5/14 | **+0.275 [+0.139,+0.410] SIG** 2/17 |
| **vqe − tilt** | rand-75 | **+0.392 [+0.097,+0.683] SIG** | **+0.429 [+0.215,+0.637] SIG** | **+0.291 [+0.129,+0.435] SIG** |
| **vqe − tilt** | set mean | **+0.470 [+0.128,+0.811] SIG** | **+0.350 [+0.214,+0.487] SIG** 1/18 | **+0.250 [+0.137,+0.359] SIG** 2/17 |
| vqe − anneal (equal draws, 2,048 evals) | rand-5 | −0.239 [−0.518,+0.019] NULL | −0.149 [−0.310,+0.024] NULL | −0.015 NULL ≤FLOOR |
| vqe − anneal | rand-75 | −0.034 NULL ≤FLOOR | **+0.141 [+0.058,+0.228] SIG (worse)** | **+0.165 [+0.049,+0.284] SIG (worse)** |
| vqe − anneal_cost (equal evaluations, 819,200) | rand-5 | +0.068 NULL ≤FLOOR | +0.141 [−0.027,+0.319] NULL | **+0.367 [+0.107,+0.633] SIG (worse)** |
| vqe − anneal_cost | rand-75 | +0.097 [−0.078,+0.254] NULL | **+0.325 [+0.184,+0.459] SIG (worse)** | **+0.302 [+0.120,+0.477] SIG (worse)** |
| vqe − uniform (**WEAK CONTROL, proves nothing**) | rand-5 | −0.541 SIG | −0.450 SIG | −0.316 SIG |

**7a. THE DECISIVE ROW.** `tilt` is `p ∝ p_untrained · exp(−E/T)` with `T` solved so its
entropy equals the VQE's — a **one-line classical reweighting of the same untrained circuit's
own distribution by the same objective**, matched on entropy and therefore on diversity, using
**4,096 objective evaluations against the VQE's 819,200 (a 200× cheaper arm)**. It beats the
CVaR-VQE **on every readout at every alpha**, by **+0.25 to +0.48 Å**, with intervals excluding
zero, W/L up to 1/18, and DIFFUSE concentration passes at the target unit. Every one of these
margins is 3–6× the false-positive floor.

**7b. AT MATCHED BUDGET, CLASSICAL SIMULATED ANNEALING IS AT LEAST AS GOOD.** At equal readout
draws (`anneal`, 2,048 evaluations — a **400× smaller** objective budget than the VQE's) the
comparison is NULL on rand-5 at every alpha and **significantly against the VQE on rand-75 at
alpha = 0.25 and 0.05**. At equal objective evaluations (`anneal_cost`, 819,200) the VQE is
significantly worse on three of six cells and never significantly better.

**7c. THE ONE ARM THE VQE BEATS IS THE ONE THE BRIEF SAYS PROVES NOTHING.** It beats uniform
random by 0.32–0.54 Å. The brief records that a zero-information constant alpha-helix beats
uniform random by 0.457 Å.

> **The ensemble effect is real, it is a location effect, and it belongs to the OBJECTIVE.
> Any method that concentrates probability mass on this objective's low-energy region gets it,
> and two classical methods get more of it for 1/200th and 1/400th of the objective budget.**

---

## 8. THE MECHANISM, measured — it is the objective's quality, target by target

**ORACLE DIAGNOSTIC.** On the same 12-qubit sub-registers, the objective's **certified global
minimiser** (known exactly because the space is enumerated) sits at ORACLE RMSD **2.537 Å**
against a space mean of **3.992 Å** and a space best of **1.384 Å** — mean percentile
**0.208**. On this sub-register the native-free structural objective is a *good* objective:
its own optimum is 1.455 Å better than a random draw.

Per target, the paired VQE-minus-untrained difference tracks that quality:

| alpha | rho(per-target diff, objective's argmin **percentile**) | rho(per-target diff, objective's argmin **RMSD**) |
|---|---|---|
| 1.00 | **+0.524** | +0.396 |
| 0.25 | **+0.635** | +0.600 |
| 0.05 | **+0.690** | +0.667 |

The extremes make it concrete. On **2MK7** and **7N2I**, where the objective's certified argmin
is at percentile 0.000 and 0.001, the VQE "wins" by **−2.016** and **−1.825 Å**. On **6F3V**
and **1CS9**, where the objective's argmin is at percentile 0.939 and 0.938, it **loses** by
**+1.250** and **+0.742 Å**. **The sign of the project's only positive quantum result is set by
whether the objective is any good on that target** — a property of the objective, observable
without any circuit.

This also explains why the result appears *here* and not in the rest of the programme: the
12-qubit sub-register freezes seven of nine or ten residues at a fixed configuration, and on
that reduced space this objective's optimum is genuinely good, where on the full register and
on real pools it is not.

---

## 9. R3 — THE ALPHA MAP. The non-monotonicity is REFUTED.

**DEMONSTRATED.** 760 cells (10 alphas × 19 targets × 4 seeds). Target unit, rand-5
coordinate average, ORACLE post-hoc.

| alpha | 1.00 | 0.75 | 0.50 | 0.35 | 0.25 | 0.15 | 0.10 | 0.05 | 0.025 | 0.01 |
|---|---|---|---|---|---|---|---|---|---|---|
| **diff** | −0.442 | **−0.565** | −0.464 | −0.310 | −0.334 | −0.302 | −0.273 | −0.295 | −0.188 | −0.158 |
| CI excludes 0 | yes | yes | yes | yes | yes | yes | yes | yes | yes | **no** |
| distinct configs | 2.5 | 24 | 57 | 93 | 117 | 166 | 202 | 242 | 269 | 324 |
| entropy (bits) | 0.70 | 1.55 | 2.50 | 3.52 | 4.02 | 4.90 | 5.40 | 5.91 | 6.16 | 6.79 |

**The response is smooth and monotone decreasing in diversity, significant at nine of ten
alphas.** The reported "win at alpha = 1 and alpha = 0.05 but not in between" was **two noisy
cells either side of a monotone decline** — at three seeds the 0.25 cell happened to land at
−0.138 and at eight seeds it is −0.295 with an interval excluding zero. **REFUTED: there is no
non-monotonicity to explain.**

The set-mean row is monotone over the whole range too (−0.807 → −0.129, significant at all ten
alphas), which is the same location effect seen without the coordinate-average operator.

---

## 10. THE COORDINATOR'S PREDICTION — the win does NOT track distinct-configuration count

**REFUTED.** The proposal was that alpha's effect is really an *effective-budget* effect —
alpha controls concentration, concentration controls how many distinct configurations a fixed
budget yields, so the win should track diversity rather than alpha, and "the arms that win may
simply be the ones that concentrate least".

Measured on 1,216 pooled cells (R1 + R3):

| statistic | value |
|---|---|
| rho(paired diff, log2 distinct) | +0.096 |
| rho(paired diff, log alpha) | −0.103 |
| rho(log2 distinct, log alpha) | **−0.904** — the two predictors are 0.90 collinear |
| OLS `diff = −0.411 − 0.0071·log2(distinct) − 0.081·log(alpha)` | R² = **0.012** |
| **mean WITHIN-alpha rho(diff, distinct)** | **+0.004** |

**Within a fixed alpha — the only place the two can be separated — the distinct-configuration
count carries nothing** (rho +0.004, running −0.15 to +0.25 across ten alphas with no trend).
And the direction of the aggregate is the opposite of the prediction: **the arms that win are
the ones that concentrate MOST** (−0.44 at 2.5 distinct configurations, −0.16 at 324).

That is consistent with §6 and §8 rather than in tension with them: concentrating on a good
objective moves the distribution's *location*, and location is the whole effect. It does not
contradict the QRESTRAINT budget observation, which is about an **argmin** readout, where a
better search genuinely can find a worse structure; it says the two are different phenomena on
different readouts.

---

## 11. R4 — THE OBJECTIVE AXIS, reported separately as the brief requires

**DEMONSTRATED.** 456 cells. The objective is rank-uniformised over the 4,096-configuration
sub-register, so a value is the mean rank-percentile of what a distribution samples: **0 is the
certified global optimum, 0.5 is a uniform random draw.**

| alpha | E_p[E] VQE | E_p[E] untrained | E_p[E] annealer's tail | annealer's **best visited E** | p(exact argmin) under the VQE | VQE mode E | VQE mode RMSD (ORACLE) |
|---|---|---|---|---|---|---|---|
| 1.00 | **0.0266** | 0.5196 | 0.1737 | **0.000000** | 0.0587 | 0.0292 | 3.101 |
| 0.25 | 0.3219 | 0.5196 | 0.1737 | **0.000000** | 0.1764 | 0.0365 | 2.754 |
| 0.05 | 0.4220 | 0.5196 | 0.1737 | **0.000000** | 0.1425 | 0.1112 | 2.822 |

**Two readings, and they must not be substituted for each other.**

1. **On the optimisation axis the VQE genuinely concentrates** — at alpha = 1 its distribution
   sits at the 2.7th percentile of the objective against the untrained circuit's 52nd, which is
   what "the optimiser works" means. But **classical simulated annealing at 2,048 objective
   evaluations reaches the objective's certified global optimum in 100% of cells** (best
   visited E = 0.000000 exactly), where the VQE with **819,200** exact-gradient evaluations
   places only **5.9%** of its mass on it. On the objective axis the classical arm wins by a
   distance at 1/400th the budget.
2. **On the structural axis that concentration is worth +0.385 Å against the untrained circuit
   and −0.398 Å against a classical reweighting.** The two axes agree in direction here (unlike
   the rest of the programme) precisely because on this sub-register the objective is good
   (§8) — but the classical arms convert the same objective into the same or better structures
   far more cheaply.

---

## 12. THE SECOND POSITIVE CELL — the QRESTRAINT enrichment, replicated and controlled

**DEMONSTRATED.** `s15/qens_enrich.py`, 135 cells (9 n=9 targets × 5 objectives × 3 seeds),
exact enumeration of the 2^18 register, alpha = 0.25, hard budget **8,192** objective
evaluations. Every number is `P(ORACLE CA-RMSD < 2.0 Å)` under the arm's own distribution.

**12a. THE ENRICHMENT REPLICATES EXACTLY.** Re-running the VQE independently per seed
reproduces the recorded 3-seed means to the printed digits: **5.20× / 6.16× / 0.98× / 6.70× /
18.64×** on `E_ls_pred` / `E_ml_pred` / `E_ls_pool` / `E_combined` / `E_ORACLE_true`. The
internal control holds — the weakest objective gives no enrichment.

**12b. THE CLASSICAL CONTROL, at the VQE's exact budget.** `× untrained`, mean over 27 cells;
entropy and effective support printed because an enrichment bought by collapse is not one.

| arm | evals | `E_ls_pred` | `E_ml_pred` | `E_ls_pool` | `E_combined` | `E_ORACLE_true` (ORACLE) | entropy (bits) |
|---|---|---|---|---|---|---|---|
| **vqe** | 8,192 | **5.20** | **6.16** | 0.98 | **6.70** | **18.64** | 9.4–9.9 |
| untrained | 0 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 12.74 |
| **untr_matchdiv** (matched distinct count) | 0 | **1.09** | **1.08** | 1.04 | **1.05** | **1.11** | 9.6–10.1 |
| uniform | 0 | 1.18 | 1.18 | 1.18 | 1.18 | 1.18 | 18.00 |
| **anneal** (whole trajectory) | **8,192** | 3.53 | 5.50 | 1.25 | 5.19 | 10.99 | 12.0–12.5 |
| **anneal_tail** (converged tail) | **8,192** | **5.34** | **8.07** | 0.68 | **8.00** | **24.31** | 8.5–9.0 |
| topM (objective's best M, M = the VQE's eff. support) | 262,144 | **10.73** | **14.16** | 1.14 | **13.08** | **51.69** | matched |
| boltz (Gibbs at matched entropy) | 262,144 | **10.47** | **13.61** | 1.03 | **15.25** | **55.65** | matched |
| boltz_p0 (untrained reweighted, matched entropy) | 262,144 | 5.57 | 9.08 | 1.24 | 7.20 | 24.27 | matched |

Paired at the target unit (n = 9), **POSITIVE = the VQE puts more mass below 2.0 Å**:

| comparison | `E_ls_pred` | `E_ml_pred` | `E_ls_pool` | `E_combined` | `E_ORACLE_true` |
|---|---|---|---|---|---|
| vqe − untrained | +0.053 SIG | +0.065 SIG | −0.000 NULL | +0.072 SIG | +0.223 SIG 9/0 |
| **vqe − untr_matchdiv** | **+0.052 SIG** | +0.064 NULL | −0.001 NULL | **+0.072 SIG** | **+0.222 SIG 9/0** |
| **vqe − anneal** | +0.021 NULL | +0.008 NULL | −0.003 NULL | +0.019 NULL | +0.097 SIG |
| **vqe − anneal_tail** | −0.002 NULL | −0.024 NULL | +0.004 NULL | −0.016 NULL | **−0.072 SIG (worse) 1/8** |
| vqe − topM | −0.070 NULL | **−0.101 SIG (worse)** | −0.002 NULL | −0.081 NULL | **−0.418 SIG (worse) 0/9** |
| vqe − boltz | −0.067 NULL | −0.094 NULL | −0.001 NULL | **−0.108 SIG (worse)** | **−0.468 SIG (worse) 0/9** |

**Three findings.**

1. **The enrichment is NOT free diversity.** The untrained circuit sampled to the trained
   state's own distinct-configuration count gives **1.04–1.11×** — nothing. The coordinator's
   literal control confirms the enrichment requires the objective, and the VQE beats it
   significantly.
2. **A classical sampler at the identical 8,192-evaluation budget reaches the same
   enrichment.** `anneal_tail` gives **5.34 / 8.07 / 0.68 / 8.00 / 24.31×** against the VQE's
   **5.20 / 6.16 / 0.98 / 6.70 / 18.64×**, and the paired difference is **NULL on all four
   native-free objectives** and **significantly against the VQE on the ORACLE objective**
   (−0.072, 1/8). It reproduces the internal control too (0.68× on `E_ls_pool`), so it is not
   an artefact of the readout either.
3. **At exactly matched diversity, a one-line classical use of the same objective roughly
   doubles the enrichment.** `topM` — uniform over the M lowest-objective configurations with
   M set to the trained state's own effective support — gives **10.7 / 14.2 / 1.1 / 13.1 /
   51.7×**, and the entropy-matched Gibbs distribution gives **10.5 / 13.6 / 1.0 / 15.3 /
   55.7×**. Both read the whole tabulated objective (262,144 evaluations, only possible on an
   enumerated instrument), so they are diagnostics of what the objective makes available rather
   than budget-matched competitors — but they bound the VQE from above by a factor of two.

**12c. A REPORTING CORRECTION: the headline ratio is a ratio of means and the typical target
sees far less.** Per-target ratio distribution on the native-free objectives:

| objective | ratio of means (headline) | **median per-target ratio** | targets with any enrichment |
|---|---|---|---|
| `E_ls_pred` | 5.20 | **1.13** | **5 of 9** |
| `E_ml_pred` | 6.16 | **2.07** | **5 of 9** |
| `E_combined` | 6.70 | **2.66** | **5 of 9** |
| `E_ORACLE_true` (ORACLE) | 18.64 | 14.84 | **9 of 9** |

**On the native-free objectives the enrichment is carried by five of nine targets and the
median target sees 1.1–2.7×, not 5–7×.** The null-calibrated concentration check agrees: the
`vqe − untrained` rows are **CONCENTRATED (FAIL)** on three of five objectives. On the ORACLE
objective the effect is uniform across all nine targets and DIFFUSE — which is exactly the
signature of an effect that is really about how good the objective is.

> **Both positive cells collapse to the same explanation.** The unranked-ensemble win (§5–§8)
> and the sub-2 Å enrichment are two readouts of one thing: **a distribution concentrated on a
> good objective's low-energy region**. Classical simulated annealing at the same budget gets
> the enrichment; a classical entropy-matched reweighting gets twice it, and beats the VQE on
> the ensemble readout as well. Neither cell survives as evidence that *only* a VQE does it.

---

## 13. WHAT I REFUTED, INCLUDING MY OWN EXPECTATIONS

**Of the record I inherited:**

1. **"The win is non-monotone in alpha."** REFUTED — smooth and monotone over ten alphas
   (§9). The original's non-monotonicity was seed noise at n = 3.
2. **"The original result is significant."** It is, at the cell unit; at the target unit on
   its own data the alpha = 1 headline is NULL and the alpha = 0.05 row fails the
   null-calibrated concentration check (§3). The effect itself survives replication at power,
   so this is a correction to the *inference*, not to the phenomenon.
3. **`qgeom_ens.unranked_readouts`'s `set_coordavg_rmsd` is computed on a SORTED, not random,
   subsample** of the distinct configurations, so it is biased differently for arms with
   different distinct counts (§2b). The unbiased version is added and used here.
4. **The coordinator's relayed prediction that the win tracks distinct-configuration count
   rather than alpha** — REFUTED. Within a fixed alpha the correlation is **+0.004** on 1,216
   pooled cells, and the aggregate direction is the opposite of the prediction (§10).
5. **The QRESTRAINT enrichment as evidence of what a VQE uniquely does** — REFUTED by the
   matched-budget classical control that workstream itself named as missing (§12).

**My own, kept in place with the evidence:**

6. **I expected the replication to WEAKEN or fail.** It did not: at power the effect is
   cleanly significant at two of three alphas at the target unit and at nine of ten in the
   alpha map, and it is *larger on the held-out ten targets*. I was wrong about the
   phenomenon and right only about its interpretation.
7. **I expected the matched-diversity control to be the decisive test.** It is informative
   (§6c) but not decisive: at alpha = 1 the match degenerates to D approximately 2 and
   *manufactures* a larger apparent win, which is an artefact of the matching and not
   evidence. **The decisive control turned out to be the classical entropy-matched tilt**,
   which I added as a secondary arm. Recorded so the next agent does not repeat the design
   error.
8. **I expected `anneal_cost` (equal objective evaluations) to be the binding budget
   convention.** It is not the tightest arm — `anneal` at 1/400th of the budget is already
   enough, and `tilt` at 1/200th is stronger still. The equal-cost convention was the less
   informative of the two here.

---

## 14. HONEST LIMITATIONS

* **Everything is on 12-qubit sub-registers of 19 enumerated targets at k = 4**, with seven of
  nine or ten residues frozen. §8 shows this matters: on the reduced space the native-free
  structural objective's certified optimum is 2.537 Å against a space mean of 3.992 Å, i.e.
  genuinely good, which is *not* the regime the full instrument is in. The result is a
  statement about a 4,096-configuration space, not about the folding problem.
* **The `tilt` arm is given the VQE's final entropy** as its temperature target, so it is a
  diagnostic and not a standalone competitor. A native-free rule for choosing that entropy was
  not measured. `anneal` and `anneal_tail` are standalone and budget-matched, and they carry
  the same conclusion more weakly.
* **`topM` / `boltz` / `boltz_p0` read the whole tabulated objective** (4,096 or 262,144
  values), which only an enumerated instrument allows. They bound what the objective makes
  available; they are not deployable methods.
* **`anneal` and `anneal_tail` are not exactly entropy-matched** (12.0–12.5 and 8.5–9.0 bits
  against the VQE's 9.4–9.9). Their entropies bracket the VQE's, which is why both are
  reported.
* **The enrichment cell is nine targets, three seeds, one alpha** (0.25), inherited from the
  design being controlled. Its per-target ratio distribution is bimodal and its
  `vqe − untrained` rows fail the concentration check on three of five objectives.
* **The `geo-mean ratio` column in the enrichment ratio table is degenerate** wherever an arm
  has exactly zero mass below 2.0 Å on some target (it prints 0.00). Read the median ratio and
  the fraction of targets above 1× instead; both are reported.
* **Two design parameters were cut for resource reasons mid-run** (§4c), by index and not by
  outcome.

---

## 15. WHAT REMAINS OPEN

1. **Does the location effect survive on the full register, where the objective's optimum is
   NOT good?** §8 predicts it reverses. That is a falsifiable prediction and it was not run.
2. **A native-free rule for choosing the tilt temperature.** If one exists, `tilt` becomes a
   deployable method rather than a diagnostic — and on this instrument it is worth 0.25–0.48 Å
   over the VQE at 1/200th the cost.
3. **Whether an entropy-matched classical tilt of a retrieval pool** (rather than of a circuit
   distribution) beats the production consensus operator. That is the same mechanism aimed at
   the pipeline the project actually ships, and it is cheap.
4. **`anneal_tail` beat the VQE significantly on the ORACLE objective and only there.** Whether
   that is objective quality (§8's law extended to the enrichment readout) or a property of
   annealing was not separated.

---

## 16. REPRODUCTION

    python -m s12.instrument            # pinned constants, run at start and at end

    python -m s15.qens_verify           # V1-V5: instrument, Struct, objective, determinism, timing
    python -m s15.qens_repl R0          # bit-level reproduction of the original 108 E2 cells
    python -m s15.qens_repl R1          # the replication at power + matched diversity (456 cells)
    python -m s15.qens_repl R3          # the alpha map, 10 alphas (760 cells)
    python -m s15.qens_repl R4          # the objective axis (456 cells)
    python -m s15.qens_report           # every table in sections 3, 5-11
    python -m s15.qens_enrich           # the classical control for the QRESTRAINT enrichment
    python -m s15.qens_enrich report    # its tables only, from the checkpoint

Every module checkpoints after every cell via `qens_lib.ck` (merge-on-write, retries the
Windows `os.replace` lock) and skips cells already present, so an interrupted run resumes by
re-running the same command. Seeds are explicit and process-stable: `s15.seed.stable_rng`
(blake2b) everywhere except the R0 reproduction, which deliberately uses the original's
`np.random.default_rng(1000*s+7)` / `default_rng(s)` convention so that it reproduces
bit-for-bit. `OMP_NUM_THREADS = MKL_NUM_THREADS = OPENBLAS_NUM_THREADS = 2` are set at module
import and on every launch.

Artefacts: `s15/results/qens_verify.json`, `qens_repl.json` (R0/R1/R3/R4), `qens_report.json`,
`qens_enrich.json`, and the console logs `qens_report.log` and `qens_enrich_report.log`.

### Leakage audit

No benchmark file was read or written. Grepping the four `qens_*` modules for
`benchmark_manifest`, `peptide_folds`, `peptide_clusters`, `dev24` and `dev_set` returns
**nothing**. The 60-target protected benchmark was never touched, looked for, or referenced.

Native quantities enter only as: `Enum.rmsd` and `s12.instrument.load_univ(...)["nat_ca"]` for
**post-hoc scoring**, and `E_ORACLE_true` in §12, which is an inherited **ORACLE DIAGNOSTIC**
objective labelled as such in every row it appears in. No objective, circuit, temperature
solve, annealing schedule, subsample or seed used a native quantity. No `qens_*` module
imports `LogPTable` / `ShiftedLogP` / `SumLogP`; the only route to them is through
`s15.qrestraint`'s own objective construction in §12, which uses the **corrected** class (its
own §0.2 checks this). No AMBER quantity enters any measurement here, so the
`amber_kind == 0 AND amber_idx != snap_index` rule is not engaged.

Nothing under `core/`, `s5/ s7/ s8/ s9/ s12/ s13/ s14/` or `tests/` was modified.
`s15/qgeom_*.py` and `s15/qrestraint.py` were imported and not changed.
