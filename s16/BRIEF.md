# SPRINT 16 — SHARED AGENT BRIEF

Binding on every agent. Read it before touching anything. Sprint 15's brief
(`s15/BRIEF.md`) remains in force except where this supersedes it.

---

## 0. THE SPRINT'S QUESTION

> **Can native-free information about the DIRECTION of structural error be combined with
> native-free information about structurally QUIET torsional directions to steer genuine
> CVaR-VQE toward conformations with substantially lower Cα-RMSD?**

**Targets:** primary **< 2.5 Å** mean held-out Cα-RMSD; preferred **< 2.0 Å**. These are goals,
not outcomes. **A 2.55 Å result that survives every control is scientifically preferable to an
artificial 1.99 Å.**

**Four mandatory pillars, all genuine and all separately measurable:** VQE · CVaR · Legacy/MJ ·
all-atom AMBER (ff14SB/GBn2). Everything else is negotiable.

---

## 1. WHY THIS HYPOTHESIS — the Sprint 15 evidence it rests on

| finding | value |
|---|---|
| ORACLE rotation of the fit's own error out of the loud half, **at fixed magnitude** | **1.855 Å** — −1.822 [−2.037, −1.613], W/L 121/5 |
| best native-free arm, of 22 tried | −0.007 [−0.074, +0.055] — **≈4% of the oracle gain** |
| native-free error-direction surrogate `θ_fit − θ_pool` | \|cos\| **0.390** with the true error vs a **0.157** null; alignment corr **+0.499 [+0.33, +0.65]** |
| native-free quiet-subspace estimate | bottom-half subspace overlap cos² **0.827** vs a **0.499** null; spectra correlated **+0.980** |
| **the two have never been combined** | this sprint's flagship |

**The prior is not high, and say so.** \|cos\| = 0.390 explains ≈15% of the error direction's
variance, and the standing ledger records that gain from a noisy channel goes as the **square** of
its skill. The pooled partial correlation of ΔRMSD with Δalignment given Δraw error is **+0.087**,
against β = +0.586 for plain distance MAE. **An honest negative here is a real result** — it would
close the largest unexplained ceiling in the programme with a mechanism rather than with 22 failed
arms.

---

## 2. ESTABLISHED SPRINT 15 RESULTS — treat as given unless your audit disproves them

- incumbent **3.204 Å**; Sprint 15 generative cascade **3.321 Å** (+0.117 [+0.050, +0.188])
- ensemble ORACLE best **2.760**; selection **+0.751**, aggregation **−0.360**, projection **+0.170**
- true-distance continuous-torsion ORACLE **≈0.6 Å** (sd 0.132 across start draws); predicted **3.644 Å**
- **error shape beats magnitude**: outlier-shaped error at 4.12 Å RMS emits 1.993 Å; i.i.d. at
  3.00 Å emits 2.561 Å; the real distogram at 3.70 Å emits 3.644 Å
- **the errors are geometrically realizable** — 2.4× closer to realizable than matched noise, ratio
  0.413, paired −1.391 [−1.623, −1.184]
- fusion law `d_avg ≈ √(r² − (s/2)²)`, **parameter-free**, predicts to 0.162 Å on 126 targets, and
  **s is native-free**
- **perfect decorrelation bottoms out at 2.320 Å** — 2.5 Å needs 86% of the coherence destroyed and
  2.0 Å needs *both* a smaller error and an incoherent one
- **five** exactly-null torsion directions, not four (the fifth a distributed whole-chain crank);
  participation ratio **3.47** effective directions out of ~25.9
- CVaR-VQE: every effect **dissolves** under a classical Boltzmann reweighting of the same untrained
  circuit at **1/200th** the budget; annealing finds the certified optimum in 100% of cells at 2,048
  evaluations where the VQE puts 5.9% of its mass using 819,200
- the per-target VQE win correlates **+0.52 to +0.69** with where the objective's certified argmin
  sits in the true RMSD distribution — **concentration pays when the objective's optimum is the
  answer**
- AMBER's demonstrated role is **stereochemical repair** (0.466 → 0.874 Ramachandran-favoured,
  116W/2L; 1.397 → 0.000 clashes, 63W/0L), **not** native ranking; its accuracy claim is WEAKENED
- Legacy's defensible role is **detection**, not ranking; its **eleven weights were never fitted**

---

## 3. HARD RULES

### 3.1 Leakage — absolute

Native structures may be used **only** for evaluation, clearly-labelled ORACLE diagnostics,
auditing, and post-hoc interpretation. **No** native-derived quantity may enter predictive model
parameters, architecture decisions after benchmark exposure, candidate selection, target-specific
tuning, thresholds, learned weights, hyperparameters, or stopping rules.

Maintain a visible separation between **predictive / tuning / oracle / diagnostic / benchmark**.
Label every ORACLE arm in **every table row and in prose**.

### 3.2 The 60-target benchmark stays SEALED

It was not read in Sprint 15 and is not to be read now. Unlock requires **all** of: the architecture
frozen, the control frozen, the success criterion pre-registered, the expected result would
genuinely change the scientific conclusion, and no cheaper internal instrument can answer it.
Architecture discovery happens on the 126-target instrument.

### 3.3 RMSD, frozen and unchanged

Cα only · every residue including both termini · no trimming · index correspondence · uniform
weights · **proper rotations only** · model 1 of the native · chain breaks scored, not rejected.
Full-chain versus `[1:-1]` differs by up to 2.4 Å on one structure here.

### 3.4 Statistics

Paired, fold-aware bootstrap 95% intervals via `I.paired`. Report **n**, **median beside mean**,
per-fold detail, and the W/L. **The statistical unit is the TARGET** — treating (target × seed)
cells as independent overstated n threefold in Sprint 15 and turned a null into a "significant"
result. Run the null-calibrated concentration check on any near-even W/L with an interval excluding
zero. Report pre-specified statistics as **one PASS/FAIL verdict**, never as fields to select from.

### 3.5 Controls — no single budget convention

Every VQE claim reports against, as applicable: **best-of-N from the untrained circuit**;
**matched diversity** (equal distinct-configuration count); **matched cost** in objective
evaluations *and* wall-clock; **classical Boltzmann reweighting of the same circuit's samples**;
simulated annealing; greedy; exact enumeration where the space allows.

**Uniform random is not a control** — a zero-information constant α-helix beats it on this
instrument. **Never report one budget convention alone when the other reverses the conclusion** —
equal-iterations gave QNG 5/5 wins where equal-cost gave 0/7.

### 3.6 The empirical false-positive floor, and its limit

A comparison that is **zero by construction** returned **+0.081 [+0.014, +0.169]** on one start
draw and −0.003 on another. **Effects at or below ≈0.08 Å are unresolvable without replication
across draws** — in pipelines that *have* a random multi-start.

**The floor transfers only if its MECHANISM does.** It is a per-process-salted random multi-start.
The AMBER k = 30 pipeline has **no RNG** and reproduces bit-identically; applying the floor there
was wrong. **Check the mechanism before importing the number.**

### 3.7 Seeding and reproducibility

Use `s15/seed.py`'s `stable_rng` / `stable_seed` (blake2b). **Never `hash()`** — Python salts string
hashes per process, which made no multi-start constant in Sprint 15 reproducible. Bit-for-bit
reproducibility is preferred wherever the computation is deterministic; document any trade.
**A determinism check that never starts a second interpreter cannot see this class of bug.**

### 3.8 AMBER data rules

- **Add and use a convergence gate.** `refine_coords` has none: 1D6X, 1MF6, 2NB7, 7BX2 end
  minimisation above 1000 kcal/mol (1MF6 at **8.9 × 10⁸**) and are silently scored. Gated, the
  exact frame null collapses from +0.0117 to −0.0005 [−0.0046, +0.0035].
- **Binding rule** for every tail, in-band, argmin or top-k statistic:
  `amber_kind == 0 AND amber_idx != snap_index`, with per-target n printed.
- **`s14/results/obj_floor.json` is RETRACTED.** Do not read it.
- 37.6% of cached single points exceed 1e6 kcal/mol, 3.2% exceed 1e12, max 2.6e20; the 99th
  percentile is itself 1e12–1e15, so **winsorising there is not sufficient conditioning**.
- Cost: a single point is 8–14 ms, but the k = 30 arm is a **full minimisation at 12.57 s mean**,
  ≈1,500× a single point and ≈2.8× the projection. Budget-match at the *right* number.
- **Preserve frame invariance** and test it: the physics is rigid-invariant, so a rotated-frame null
  must be zero.

### 3.9 Legacy

Genuine Legacy/MJ, identifiable and separately measurable. **Do not silently substitute a learned
surrogate.** Its eleven weights were **never fitted** — resolve this honestly by either fitting them
leave-fold-out (and reporting the ablation) or characterising them as fixed and limiting the claim.
**Never compare a fitted and an unfitted system in a way that implies equal training status.**

### 3.10 Known traps

- **`I.FAIL18` is a threshold artefact** — `BAND = 1.5 Å` has no derivation, the zero-recall count
  runs 45 → 2 across BAND 0.5 → 3.0, and only 1 of the 18 targets survives every threshold. Report
  the column for continuity; **no argument may rest on set membership**.
- **The set-mean law is out of domain on a collapsed set.** Report the set's **diversity** beside
  any aggregate; a method that concentrates improves the proxy while degrading the output.
- **A raw drop-top threshold is not a test of concentration** — it needs a uniform-effect null.
- **`np.argmin` on a tied signal reads the ORACLE sort order.** Average over the tied argmin set.
- **Band/component decompositions must print the JOINT value** — single-band effects summed to
  +0.785 against a joint +0.381, a factor of two.
- **`I.write`'s completeness flag has a 12/12 false-alarm rate.** Check row counts, not the flag.
- **`n ≤ 8` reads have reversed conclusions five times.** A 3-target probe of a heavy-tailed
  quantity measures the wrong thing.

---

## 4. MACHINERY THAT EXISTS — use it, do not rebuild it

**`s12/instrument.py`** — `I.targets()` (126, pinned), `I.load_univ`, `I.build_ca`, `I.ca_rmsd`,
`I.kabsch_rmsd_batch`, `I.project`, `I.coordinate_average`, `I.pairwise_rmsd`, `I.medoid`,
`I.superpose_batch`, `I.distogram`, `I.pair_index`, `I.paired`, `I.summary`, `I.write`, `I.FAIL18`.
Must reproduce `shipped 3.4540004952559396 · pool_best 1.7108244199364904 · top75_best
2.3061526409453816 · synthesis_fit 3.2040761603809194 · n_zero_recall 18`. **This is a cache read,
not a pipeline reproduction.**

**`s15/align_lib.py`** — the flagship's machinery. `sup_jacobian` (exact analytic
∂(superposed CA)/∂(φ,ψ), verified cosine 1.000000000000000), `spectrum`, `alignment(J, e)`,
`energy_split`, `rigid_basis`, and **`fit(..., basis=...)`** which optimises over
`θ = θ₀ + basis @ c` with an exact chain rule — the subspace-constrained arm is already built.

**`s15/`** — `distgeo` (restraint fits), `distml` (`LogPTable`, corrected non-uniform grid),
`pooldist` (pool distance channel), `robust` (`fit_robust`, `fit_gnc`, the `rho` family),
`cascade` (G→S→A→F), `coherence` (realizability + fusion law), `decorr` (the coherence ladder),
`feasible`, `expand`, `seed`, `figures`.

**`s14/vqe_lib.py`** — `Enum` (fully enumerated targets), `Counter` (hard budget),
`search_random/greedy/anneal/ga`, `summarise`, `spearman`, `wait_for_memory`, and crucially
**`blend_objective(base, truth, signal)`** and **`noisy_truth`** — rank-space objective-quality
knobs written for the phase-boundary question, with no scale or tail-shape confound.
**`s14/vqe_run.py`** — `make_ansatz`, `init_theta`, `run` (CVaR-VQE with hard budget).

**`core/quantum.py`** — `StatevectorCircuit`, `MPSAnsatz`, `cvar_exact`,
`grad_cvar_paramshift/fd/score`, `Adam`, `mem_gate`.

**Enumerated spaces** — `s13/results/qarch_enum_*.npz` (9 targets, n = 9, k = 4, 262,144 configs)
plus 10 more n = 10 targets in `s14/cache/` — **19 targets, 1.28e7 exactly labelled structures**.

---

## 5. THE MACHINE

15.6 GB RAM; 8 cores that are 4 fast + 4 slow ≈ **6.43 core-equivalents**. Target ~95% CPU,
**never exceed 97%**. Two heavy jobs fit comfortably; more contend. Set
`OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2` at the top of every module. Checkpoint
long runs to JSON every 10 targets. Before launching anything heavy:

```
powershell -NoProfile -Command "$c=Get-CimInstance Win32_Processor|Measure-Object -Property LoadPercentage -Average; $o=Get-CimInstance Win32_OperatingSystem; '{0}% CPU, {1:N1} GB free' -f [int]$c.Average,($o.FreePhysicalMemory/1MB)"
```

Hold if CPU > 90%. `s14.vqe_lib.wait_for_memory(min_gb, tag)` blocks until RAM frees.

---

## 6. HOW TO WORK

For every candidate: **state the hypothesis → identify the causal mechanism → define the control →
define the evaluation unit → define the success criterion → run the ablations → inspect the failures
→ check it replicates → audit it adversarially → promote, kill, or quarantine.**

Before implementing anything, check the Sprint 5–15 record and state whether it is new, a variation
of a killed idea, already tested, or useful only as a diagnostic. **Do not rename a dead idea.**

**A refutation is as valuable as a confirmation.** Preserve your own errors and retractions in
place — never edit them away. Do not soften an inconvenient result. Report negatives plainly.

Write findings to `s16/<name>_FINDINGS.md`, tiering every claim **DEMONSTRATED / ORACLE DIAGNOSTIC /
HYPOTHESIS / REFUTED**. Do not send interim reports; work until done, then write the file.

**The final report must let an expert determine:** the hypothesis and why it follows from Sprint 15,
what changed and what was held fixed, the controls, what the quantum component actually did, what
CVaR actually changed, what Legacy and AMBER each contributed, which effects were causal, which were
oracle-only, which generalised, which failed, how uncertainty was handled, what was pre-registered
versus discovered, what the literature already established, what is genuinely new, and what remains
unresolved. **Not marketing copy. No inflated novelty. No hidden failures.**
