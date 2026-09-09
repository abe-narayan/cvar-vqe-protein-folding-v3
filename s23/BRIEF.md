# SPRINT 23 — FINAL. DRIVE THE DEV MEAN BELOW 3.0 Å. PRESERVE THE ARCHITECTURE.

## 0. THE THREE PILLARS — NON-NEGOTIABLE

1. **CVaR-VQE remains the central selector.** Classical methods may be diagnostics, proposals,
   refinement, diversity, generation, post-processing or controls — **not the selector**.
2. **`H_Legacy` genuine and independently evaluable** — the 11-term potential at `DEFAULT_WEIGHTS`,
   never fitted, never replaced by a proxy.
3. **`H_AMBER` genuine and independently evaluable** — ff14SB/GBn2. Distinguish its ROLES:
   Hamiltonian / training objective / readout objective / refinement / landscape probe.

**PRIMARY ENDPOINT: mean full-chain Cα-RMSD, 126 cluster-disjoint dev targets. Incumbent 3.048 Å
(point cloud, score-filter → top-75 coordinate average, K=500 pool). TARGET < 3.0 Å.**
**The 60-target benchmark is SEALED.** No inspection, no tuning, no threshold selection against it.

---

## 1. DO NOT REDO — closed by prior sprints, with evidence

| closed | evidence |
|---|---|
| More VQE search | whole `2**n` latent enumerated; exact argmin **ties** a zero-eval pool; 0/122 reversals |
| Deeper ansatz / larger χ | layers 1→20 swept; on AMBER **best-of-N wins the ladder**; χ orders nothing |
| Different optimiser | none beats best-of-N on RMSD; training the circuit at all = +0.003 Å |
| Encoding θ vs (sin,cos) | 5,760 cells, 5 legs; **radius ≡ scale, the same knob** |
| Changing the TRAINING Hamiltonian | +0.056 [−0.020,+0.129], null |
| Replacing the candidate SOURCE | generative latent adds **−0.103 Å at the ORACLE of the union** |
| The consensus family | 9 native-free arms, **0%** of the in-pool gap |
| Legacy→AMBER continuation | staged **worse** than direct; raw λ never visits an intermediate H |
| Naive torsion averaging | coordinate averaging beat it by **1.024 Å** — *but see §3, the circular-mean question* |
| Hedging across m | `hedge_all` −0.0002, `hedge_core` −0.0117 vs MDE 0.027 — falsifier fired |
| Native-free routers over arms | 4 constructions, ~0% capture, one **significantly harmful** |

**A finite-sample bound from this project's own σ ≈ 0.41 Å:** at n≈100/fold even a *one-global-threshold*
router has a generalisation gap of 0.39 Å; anything richer needs 280–1045 targets. **Do not build a
multi-feature per-target router. A single GLOBAL scalar fitted on folds is the only router class the
sample size supports.**

---

## 2. THE COORDINATOR'S SCOPING RESULT — start here

Coordinate averaging **contracts** the structure: mean virtual Cα–Cα bond **3.12 Å** against a native
**3.82 Å**. Kabsch fits rotation and translation only — **it does not fit scale** — so a systematic
scale error propagates straight into RMSD. Measured, n=40:

    incumbent (scale 1.0)                 2.8547
    ORACLE best per-target scale          2.6964   mean s* = 0.975   <- CEILING, -0.158 [-0.255,-0.085], 38W/2L
    native-free: match predicted Rg       2.9697   mean s  = 1.084   +0.115  WORSE
    native-free: least-squares fit to dhat 2.9556  mean s  = 1.079   +0.101  WORSE

**Read this carefully, it is counter-intuitive.** The optimal scale is **0.975 — slightly MORE
contraction**, not less. Both distogram-derived estimators want to **expand** (s ≈ 1.08) and both
make RMSD worse. Restoring the physical bond length (s = 1.277) is catastrophic: **+0.82 Å**.

**So: a scale correction is worth up to 0.158 Å, it is nearly universal (38W/2L), and every
physically-motivated estimator of it points the WRONG WAY.** A single global constant fitted on
training folds is the obvious first arm and is within the only router class the bound allows.

---

## 3. PRIORITY HYPOTHESES

**H1 — GLOBAL SCALE CORRECTION.** One scalar, fitted per training fold, applied at inference.
Nested CV. Also test: scale about the centroid vs about the medoid frame; scale applied before vs
after any repair.

**H2 — CLUSTER-RESTRICTED AVERAGING.** If the top-m spans two structural basins, their average lies
between both and resembles neither. `medoid_75` (pick one member, 3.282) and `avg_75` (average all,
3.048) are both tested; **"cluster the top-m, average within the dominant cluster" is NOT.** This is
the clearest untested gap in the aggregation family.

**H3 — CIRCULAR TORSION AVERAGING, DONE PROPERLY.** Prior sprints found torsion averaging loses by
1.024 Å. **Verify from source whether that test used a CIRCULAR mean** (`atan2(Σsin, Σcos)`) or a
naive arithmetic mean of angles. A naive mean across the ±π branch cut is a *bug*, not a result, and
would explain a loss of exactly that size. If it was naive, redo it properly with von Mises weighting
and basin conditioning. If it was circular, the direction stays closed.

**H4 — rg_z-GATED AGGREGATION.** The compactness-disagreement signal (partial ρ 0.34–0.38, survives
length, difficulty and realised-extension controls) is the only native-free signal that has cleared
every control. It is **not** a difficulty proxy (r = +0.126 with true difficulty). Test it as a
**gate on aggregation width or weighting**, not as a router over many arms.

**H5 — Cα-PRESERVING REPAIR.** Ideal-geometry projection and AMBER relaxation both improve validity
and can worsen RMSD. Test a repair that **fixes bonds/angles/clashes while minimising Cα
displacement** — a restrained minimisation with a strong positional restraint on Cα only.

**H6 — ADAPTIVE α.** α was swept {1.0, 0.5, 0.30, 0.15, 0.05, annealed, adaptive} crossed with three
Hamiltonians. **A fixed α is not the explanation.** Only pursue conditioning α on a native-free
quantity (rg_z, ensemble entropy) if H1–H5 leave budget.

---

## 4. HARD RULES

- **MDE = 2.8016 × SE, per comparison.** Report SE and the effect as a multiple of its own MDE.
- **Paired, fold-clustered CIs beside iid.** W/L/tie. Report worst-target degradation explicitly.
- **Basis discipline:** point-cloud and built-chain RMSD are never compared. State the basis on both sides.
- **ORACLE vs ACHIEVED vs FLOOR** labelled at every appearance.
- **Nested CV for anything fitted.** A hyperparameter fitted on the folds it is scored on is leakage.
- **Native information: evaluation only.** Never at inference.
- **Rule 0 — operator forks.** For any directional hypothesis, enumerate SIX axes in the module
  docstring and NAME THE ALTERNATIVE NOT TAKEN: functional, basis, readout, normalisation, null,
  **and THE LABEL**. Send the fork list to the coordinator BEFORE the run.
- **Best-of-K nulls are the distribution of the MAXIMUM.**
- **Completion flags demand the full key set**, not a row count. Atomic writes, config-derived filenames.
- **Ablate anything that wins.** A win without an ablation is not promoted.

## 5. LOGGING
experiment id · git commit · seed · target ids · config · Hamiltonian · ansatz · qubits · optimiser ·
α · shots/exact · iterations · ensemble size · post-processing · AMBER settings · repair settings ·
mean/median/worst RMSD · CI · per-target delta vs incumbent.
