# REPRODUCIBILITY

What reproduces, what does not, and exactly how far each claim of reproduction extends. Written to
be checked rather than believed. Where a check was run and passed, the tolerance is stated; where a
check was never run, that is stated too.

Status: **living**, completed before the protocol freeze.

---

## 1. WHAT REPRODUCES BIT-IDENTICALLY FROM PRIMARY INPUTS

Verified in the Phase 0 audit by rebuilding from source rather than reading cache.

| layer | result | scope |
|---|---|---|
| `peptide_db.npz` from 1,524 PDB files | **787/787 records bit-identical**; `ca`, `phi`, `psi` max diff exactly 0.0 | 45 s |
| all 126 window universes | **bit-identical**, 9 arrays × 126 targets | 2.35 M windows, 31 s |
| `pool_best = 1.7108244199364904` | **reproduced from the PDB files**, not from cache | 126 targets |
| `top75_best`, `n_zero_recall` | reproduced **without the pipeline cache**; the production `sub` is the score's top-75 on 126/126 | 126 targets |
| `shipped = 3.4540004952559396` | reproduced at torch **1, 2 and 8 threads**, zero argmin flips | 126/126 |
| `synthesis_fit = 3.2040761603809194` | bit-identical end-to-end — `sub`, `ca`, `fit_ca` and `amber_ca` all max-diff 0.0 | full run |
| AMBER golden −489.9138948277905 | **bit-exact**; ff14SB/GBn2 parameters verified against an independently constructed ForceField at max diff 0.0 | — |
| Legacy total | equals the `DEFAULT_WEIGHTS`-weighted sum of its 11 components, worst 2.7e-5 | **1.28e7 structures** |
| CVaR | verified against Rockafellar–Uryasev, max error **4.2e-14**; `dCVaR/dp` median FD error 0.0 | 400 cases |
| RMSD | Horn quaternion reimplementation agrees with `kabsch_rmsd_batch` to **2.04e-13 Å** | 63,000 real structures |
| analytic torsion gradient | **cosine 1.000000000000000** against central differences | real targets |
| `LogPTable` derivative (after correction) | **cosine 1.000000000000000**, max abs error 5e-09 | — |

**Checked again at the end of the sprint**, after all module patches, the seeding fix, the
`LogPTable` correction and every arm above. The instrument emits
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`.

---

## 2. WHAT DOES NOT REPRODUCE

### 2.1 The multi-start draw was not reproducible across processes — FIXED

Every multi-start module seeded its RNG with `hash(pdb)`. **Python's string hash is salted per
process** unless `PYTHONHASHSEED` is set, and the environment lock records it as `null`. Two
consecutive interpreters return `hash('1A13') = 217586290314588545` and `2408026022170661001`.

**Scope of the damage, stated precisely.**

- **Within one process: harmless.** All arms in a run share one start set, so every paired
  comparison inside a single result file is start-matched and valid. Every within-file conclusion of
  this sprint stands.
- **Across processes: not harmless.** Cross-module comparisons used different start draws, and
  multi-start constants — including the ORACLE 0.611 Å and predicted 3.644 Å — were **not
  bit-reproducible on re-running**. They are correct measurements carrying an undeclared
  run-to-run variance.

**Fixed** in `s15/seed.py` with a `blake2b`-based `stable_seed`, verified identical across two
separate interpreters (`stable_seed('1A13') = 3025285571` both times). Seven modules use it.

**Policy adopted, and recorded so it is never presented as a convenience:** exploratory runs need
not be bit-reproducible; the **frozen protocol is re-run under stable seeding before the benchmark
is touched**, and only those numbers are quoted as reproducible constants.

**Why every existing check missed it.** Determinism had been verified repeatedly, including at three
thread counts — but always *within* a process or against a *cached* artefact. **A reproducibility
audit that never starts a second interpreter cannot detect the most common source of
irreproducibility in Python.**

### 2.2 The distogram is not bit-stable across BLAS thread counts

`torch.set_num_threads` at 1, 2 and 8 gives three distinct byte patterns: `prob` differs by up to
**2.0e-6** and `risk` by up to **1.3e-4**. The cache was written at torch's default of 8; the sprint
brief mandates 2, so **following the brief changes the bits**.

No headline constant moves, and no argmin flips. But any claim of bit-identical reproduction must
**state the thread count**. All Sprint 15 results are produced at 2 threads.

### 2.3 `python -m s12.instrument` is a cache read, not a reproduction

In its 2.2 s it recomputes only the Kabsch step: `pool_best`, `top75_best` and `n_zero_recall`
recompute nothing. It also **asserts only three of its five constants** — `shipped` and
`synthesis_fit` are printed and unchecked. Running it is a consistency check on the cache, and
describing it as verifying the pipeline would be wrong.

---

## 3. ARTEFACT HYGIENE

Problems a reviewer would find, listed so they do not have to.

- **325 of 363 result JSONs carry no completeness flag.** Of the 38 that do, **all 12 "incomplete"
  flags are false alarms** — the writer takes `len()` of a dictionary of *arms*, not of targets. A
  flag with a 12/12 false-alarm rate is worse than none. **Do not cite it as evidence a run
  finished**; check the row count.
- **Five Sprint 14 results exist in two places** (`s14/results/X.json` and `s12/results/s14_X.json`)
  with no canonical copy designated.
- **`verify/amber_audit.json`'s `5_translation_invariant: 0.5446` measures the minimiser, not the
  physics.** Separated: single-point 1.3e-6, minimised 0.148 kcal/mol. The field name asserts one
  quantity and reports another, with no PASS/FAIL.
- **Sprint 14's flagship ablation ladder (4.072 → 3.485 → 3.314 → 3.204) traces cleanly but to four
  different files**, one of them misleadingly named `vqe_classical_limit.json`, with no artefact
  holding the ladder itself.
- **One result artefact is retracted outright**: every AMBER row in `s14/results/obj_floor.json`,
  computed with no `amber_kind` mask yet flagged `"complete": true`.

---

## 4. NUMERICAL CAUTIONS THAT AFFECT REPRODUCTION

- **37.6% of cached AMBER single points exceed 1e6 kcal/mol**, 3.2% exceed 1e12, maximum 2.6e20. The
  99th percentile is itself 1e12–1e15, so **winsorising at the 99th percentile is not sufficient
  conditioning**. Any spectral analysis of the unconditioned energy measures its worst steric clash.
- **The `amber_kind == 0` stratum is unbiased on the mean and not in the tails**, because the index
  array force-includes the ORACLE snap index. Binding rule: `amber_kind == 0 AND
  amber_idx != snap_index` for every tail, in-band, argmin or top-k statistic, with per-target n
  printed.
- **Plug-in mutual information is unusable on this data.** Measured directly: the plug-in estimator
  reports +0.40 to +0.88 bits where held-out predictive cross-entropy reports −0.11 to −0.23. Use
  held-out cross-entropy.
- **AMBER cost is 8.3 ms at 127 atoms rising to 23.3 ms at 237 atoms**, i.e. 8–14 ms on the
  enumerated targets — not the 28 ms previously assumed. Any budget-matched comparison computed with
  28 ms over-charged AMBER by 2–3×.

---

## 5. LOAD-BEARING CONSTANTS AND THEIR SENSITIVITY

| constant | status | sensitivity |
|---|---|---|
| `K = 500` | fitted on the reported targets | `pool_best` 1.970 → 1.504 across 100 → 2000; `shipped` 3.425 → 3.520 (**flat**) |
| `m = 75` | fitted on the reported targets | `top75_best` 2.609 → 2.106 across 25 → 150 |
| `BAND = 1.5 Å` | **no derivation exists** | `n_zero_recall` 45 → 2 across 0.5 → 3.0; only 1 of 18 FAIL18 targets survives every threshold |
| `IDENTITY_THRESHOLD = 0.6` | reused unjustified | at the null for peptides — random sequences score 0.56–0.63 |
| `min_sep = 2` | undocumented | **no sweep exists**; all pair statistics inherit it |
| `torsion_window = 8` | undocumented | not swept |
| AMBER tolerance 1.0 | undocumented | not swept |
| `bond + angle > 1000` strain gate | undocumented | not swept |
| Legacy's 11 weights | **never fitted** | documented as "a variance-balanced starting point, not a fit" |
| projection penalty | **pre-registration mismatch** | `s8/project_devarm.json` pins `phip@0.03`; production ships `ramah@0.3`. Worth 0.004 Å |

---

## 6. HOW TO REPRODUCE THIS SPRINT

1. Set `PYTHONHASHSEED=0` and `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2`, and pin
   `torch.set_num_threads(2)`. State the thread count in anything you report.
2. Run `python -m s12.instrument` and confirm the five constants. Note this checks the cache, not
   the pipeline.
3. Rebuild from primary inputs with the `s15/audit_*.py` modules if a genuine reproduction is
   wanted; expect ~45 s for the database and ~31 s for the universes, and ~90 minutes for the full
   end-to-end pipeline check.
4. Run the sprint modules in dependency order: `distacc` → `distcal` → `distgeo` → `distml` →
   `pooldist` → `feasible` → `robust` → `cascade` → `expand` → `augment`, then `quant`,
   `qrestraint`, `qgeom_*`, `info_*`, `align_*`.
5. Every module checkpoints to `s15/results/<name>.json` every 10 targets and prints its own arm
   table. Compare arm **orderings** first and absolute values second: orderings are start-matched
   within a run, absolute values carry the start-draw variance for anything produced before the
   seeding fix.

---

## 7. WHAT IS NOT REPRODUCIBLE BY DESIGN

The **60-target protected benchmark** has not been read in this sprint. Reproducing the tuning
results requires none of it. The confirmatory benchmark run happens once, after the protocol freeze,
and it is a **single-use instrument** — a reproduction that re-tunes against it is not a
reproduction of this work.
