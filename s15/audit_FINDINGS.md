# Sprint 15 — AUDIT workstream (Phase 0). Findings.

Tier legend: **DEMONSTRATED** / **HYPOTHESIS** / **REFUTED** / **ORACLE DIAGNOSTIC**.
Written continuously. Every number below has a JSON artefact under `s15/results/`.

Reproduction commands are listed in section Z.

---

## A0. Environment lock — DEMONSTRATED

Measured 2026-09-05 on the working machine. Machine-readable: `s15/results/audit_env_lock.json`.

| item | value |
|---|---|
| python | 3.13.13, conda-forge, MSC v.1944 64-bit (AMD64) |
| platform | Windows-11-10.0.26200-SP0 |
| numpy | 2.5.1 |
| scipy | 1.18.0 |
| torch | 2.13.0+cpu |
| openmm | 8.5.2 (build 8.5.2.dev-36a30cb) |
| scikit-learn | 1.9.0 |
| matplotlib | 3.11.1 |
| pennylane | 0.45.1 |
| biopython | 1.87 |
| pandas / qiskit | NOT INSTALLED |
| BLAS/LAPACK | scipy-openblas 0.3.33.112.0, USE64BITINT DYNAMIC\_ARCH NO\_AFFINITY, Haswell kernel, MAX\_THREADS=24 |
| numpy build compiler | msvc 19.44.35228 |
| OpenMM platforms present | Reference, CPU, OpenCL |

**Thread determinism — DEMONSTRATED, with a caveat that must go in Methods.**

* `numpy`/BLAS paths (Kabsch, retrieval, scoring, projection): bit-identical at every thread
  count tested.
* `core.amber` single point and minimisation: bit-identical at `threads=1` (the shipped
  setting). `core/amber.py` already documents that `threads>1` is a *different experiment*
  (converged energies drift 3.0e-4 relative), and that is confirmed rather than re-measured
  here.
* **The distogram is NOT bit-stable under `torch.set_num_threads`.** Three thread counts
  give three different byte patterns of `prob` for the same input, same session
  (sha16 `ccee63c6…` / `3ab4ca5d…` / `9599b9e3…` at 1 / 2 / 8). Within one thread count a
  repeat is bit-identical (max diff 0.0).
  * Magnitude across all 126 targets: `prob` up to **2.0e-6**, `risk` up to **1.3e-4**.
  * **The cached `s12/cache/disto_*.npz` were written at torch's default 8 threads**
    (max diff exactly 0.0 at t=8).
  * The brief mandates `OMP_NUM_THREADS=2`, which torch honours — so *following the brief
    changes the distogram bits relative to the committed cache*. It does not change any
    reported constant (below), but a "bit-identical" claim must state the thread count.

---

## A1. What `python -m s12.instrument` actually is — DEMONSTRATED

Wall time **2.2 s**. It reproduces the five pinned constants exactly:

```
shipped 3.4540004952559396  pool_best 1.7108244199364904  top75_best 2.3061526409453816
synthesis_fit 3.2040761603809194  n_zero_recall 18
```

But it is **a consistency check over four cached artefact families, not a reproduction**.
Per constant, what the 2.2 s actually computes:

| constant | recomputed in the 2.2 s | read from cache |
|---|---|---|
| `pool_best` | nothing (a `min` over stored labels) | `s8/generate_univ/*.npz`: `rr`, `order` |
| `top75_best` | nothing (a `min`) | `rr` + `sub` from `bench_results/cache/1fc9f2dcf489e2fb/*.json` |
| `n_zero_recall` | a set membership test | same two |
| `synthesis_fit` | the Kabsch only | `fit_ca` from the production cache, `nat_ca` from the universe |
| `shipped` | the Bayes-risk score and its argmin over 500×126 rows | the distogram tensors, `W`, `rr` |

Two further gaps in the self-check itself:

* `selfcheck()` **asserts only** `pool_best` and `top75_best` (tol 2e-3) and the FAIL18 set.
  `shipped` and `synthesis_fit` are printed and never asserted. A regression in either would
  pass the gate silently. **RECOMMENDATION: assert all five.**
* `stage_univ` rebuilds only universes that do not already exist on disk
  (`todo = [p for p in tg if not os.path.exists(...)]`), so "rerun the generator" is a no-op
  against a populated cache. Reproduction requires deleting first.

---

## A2. What reproduces from primary inputs — DEMONSTRATED

Primary inputs are the deposited PDB files: `pdbs/` (61), `pdbs_ext/` (1463), `prots/` (13751).

### A2.1 `peptide_db.npz` — BIT-IDENTICAL

`s15/results/audit_repro_db.json`. Full rescan of 1524 PDB files through
`peptide_db._scan()` in **45.2 s**: 787 records in, 787 records out, sequence sets equal,
zero PDB-id disagreements, and `ca`, `phi`, `psi` **bit-identical on all 787** (max abs
diff exactly 0.0 on all three arrays).

### A2.2 All 126 window universes — BIT-IDENTICAL

`s15/results/audit_repro_univ.json`. Every one of the 126 `s8/generate_univ/<pdb>.npz`
rebuilt from the library through `s8.generate`'s own code path, 2,352,893 windows total,
**30.6 s**. All nine stored arrays (`W`, `PHI`, `PSI`, `S`, `org`, `sim`, `order`, `rr`,
`nat_ca`) bit-identical on **126/126 targets**, zero exceptions.

> **`pool_best = 1.7108244199364904` is REPRODUCED FROM PRIMARY INPUTS, exactly.**

### A2.3 `top75_best` and `n_zero_recall` — REPRODUCED WITHOUT THE PIPELINE CACHE

`s15/results/audit_sensitivity.json`. The production record's `sub` was assumed to be the
shipped score's top-75. That is now verified rather than assumed: recomputing the score
from the universe and the cached distogram and taking `argsort(kind="stable")[:75]`
reproduces the production `sub` **as a set on 126/126 targets**. So `top75_best` and
`n_zero_recall` do not actually depend on the pipeline cache at all.

### A2.4 The production pipeline end-to-end — BIT-IDENTICAL on every target run so far

`s15/results/audit_pipe_run.json`. `core.pipeline.run_target` rerun from source at `PROD`
config. A full 126-target pass was launched and was at **22/126 with zero exceptions** when
this was written; it writes incrementally and can be resumed with the command in section Z.
Per target 5–15 s (1A13 14.4 s, 1A1P 11.1 s, 2MK7 5.1 s).

* `Config().key()` evaluates to **`1fc9f2dcf489e2fb`**, which is exactly the directory name
  `s12.instrument.PROD_KEY` reads. The production cache is genuinely config-keyed.
* `sub` identical, `ca` identical, **`fit_ca` bit-identical** (max abs diff 0.0),
  `synthesis_fit` delta exactly 0.0 on all three.
* **AMBER reproduces bit-for-bit too**: `amber_ca` max abs diff 0.0 and `amber_e1` equal to
  the cached value to the last digit (-347.2650824559028 / -501.0767846123149 /
  -184.89562973782432). An 800-step GBn2 L-BFGS minimisation is bit-reproducible in this
  environment at `threads=1`.

> `synthesis_fit = 3.2040761603809194` is reproduced end-to-end, from source, with
> `synthesis_fit_delta` exactly 0.0 on **22 of 126 targets** and counting; every one of
> `sub`, `ca`, `fit_ca` and `amber_ca` bit-identical, on every target run. The remaining 104
> are running. **No drift of any magnitude has been observed anywhere in the pipeline.**

### A2.5 The distogram cache — DRIFT, quantified, does not move the constant

`s15/results/audit_repro_disto.json`, `audit_determinism_t{1,2,8}.json`.

* Recomputed at torch 2 threads (the brief's setting), `prob` differs from cache by up to
  **8.9e-7** and `risk` by up to **1.9e-5** per target; over all 126, up to 2.0e-6 / 1.3e-4.
* At torch 8 threads the recomputation is **exactly 0.0** — i.e. the cache is faithful.
* **Zero argmin flips on 126/126 targets at every thread count**, and
  `shipped` recomputes to `3.4540004952559396` at 1, 2 and 8 threads — identical to the
  pinned value in all three.

> `shipped = 3.4540004952559396` is REPRODUCED from the trained model + primary universes.
> The *tensor* it comes from is not bit-stable; the *constant* is.

---

## A3. Stale-artefact and cache-key hazards — DEMONSTRATED

Machine-readable inventory: `s15/results/audit_provenance.json` (rebuild with
`python -m s15.audit_provenance`, verify with `... check`).

| artefact | config key | atomic write | verdict |
|---|---|---|---|
| `bench_results/cache/<key>/*.json` | **YES** (`Config.key()`, sha1 of `science()` + backend report) | YES (`_atomic_write_json`) | SAFE. `NOT_SCIENCE` correctly excludes thread counts, memo sizes, affinity |
| `s8/generate_univ/*.npz` | no | YES (`.tmp.npz` + `os.replace`) | SAFE against truncation; unsafe against a silent input change (no input hash stored) |
| `s12/cache/disto_*.npz` | no | no | **HAZARD**: keyed on PDB id only. A different fold model, ESM bank or torch version writes to the same filename |
| `s12/results/*.json` (225) | no | no | documented `I.write` hazard, still present |
| `s13/results/*.json` (59) | no | no | `qarch_lib.write` has **no completeness flag at all** |
| `s14/results/*.json` (51) | no | no | as above |
| `peptide_db.npz`, `peptide_clusters.json`, `peptide_folds.json` | no | no | PINNED by brief; `peptide_db.npz` verified reproducible (A2.1) |

**The `I.write` completeness flag is worse than absent — REFUTED as a safeguard.**
363 result JSONs across s12/s13/s14/s15: **325 carry no completeness flag**. Of the 38 that
do, 12 report `complete: false`. All 12 were inspected and **all 12 are false alarms**:
`I.write` looks for `obj["per_target"]` or `obj["rows"]` and takes `len()` of whatever it
finds, but these files store a *dict of arms* there, so `len()` returns the arm count, not
the target count. Example: `s14/results/ladder.json` reports `n_rows 11 / n_expected 126`
while every one of its 11 arms carries `n = 126`. A flag with a 12/12 false-alarm rate
trains readers to ignore it, which is exactly the failure mode it was added to prevent.
**RECOMMENDATION: pass the row count explicitly, never infer it.**

**Duplicated results with no canonical copy — HAZARD.** Five Sprint 14 results exist twice,
under `s14/results/<name>.json` and `s12/results/s14_<name>.json`, differing only by the
~57 bytes of completeness fields (`hamil3`, `position`, `coherence`, `ladder`, `avgspace`).
Nothing states which one a table was read from. For a paper, one must be deleted or
declared canonical.

**Stale directory still on disk:**
`distogram_models_large.STALE-PRE-FOLD-REPIN-DO-NOT-USE/` — correctly named, but present.
`s12/newpdbs/` (133 PDB files) is **not** in `core.data.DIRS` and does not enter
`peptide_db`; its provenance and purpose are undocumented.

---

## A4. Leave-fold-out discipline of the instrument — DEMONSTRATED, with two disclosures

`s15/audit_leak.py` → `s15/results/audit_leak.json`, `audit_selfmatch.json`. All 126 targets.
Nothing in `s12.instrument` asserts any of this; it is checked here for the first time.

**Clean:** the `fold` label stored in every universe npz equals the fold recomputed from the
pinned `peptide_folds.json` on **126/126** targets. Zero mismatches.

**Disclosure 1 — the 0.6 identity filter operates on whole peptides, and windows are
substrings.** On **16 of 126** targets the retrieval universe contains a window whose
sequence is *verbatim* the full sequence of a peptide-database member belonging to the
target's **own held-out fold** (examples: `AISVLLAQAVFLL` in 1CEK and 1G89's universes,
`LLGDFFRKSKEK` in 2FBU and 2LWS's, `VDIHVWDGV` in 2MK7's). The structure comes from an
out-of-fold donor, so no held-out *coordinate* enters — but the fold split does not achieve
sequence disjointness at the window level, which is the level retrieval actually operates
at. This is the same mechanism `core/data.py:31` already records ("573 library members
across 58 targets contain a target at ≥0.6 while passing the 0.6 filter").

**Disclosure 2 — on 4 targets the pool's top BLOSUM hit is the target's OWN sequence.**
1CEK, 2FBU, 2P5H, 6B9K each have exactly one window in their universe whose sequence is
*identical to the target's*, and on all four it sits at **BLOSUM rank 0** — the first
retrieved candidate, by construction. My independent count of 4/126 coincides with the
project's existing record that "4/126 targets carry a verbatim copy in their own
distogram's training set"; whether they are the same four was not checked.

**But the impact on the pinned constants is exactly zero, measured:**

| target | CA-RMSD of the identical-sequence window to the native | pool\_best | pool\_best with it removed |
|---|---|---|---|
| 1CEK | 0.594 | 0.342 | 0.342 |
| 2FBU | 3.279 | 2.243 | 2.243 |
| 2P5H | 2.334 | 1.840 | 1.840 |
| 6B9K | 4.126 | 2.077 | 2.077 |

These are *different structural determinations of the same sequence*, not the native's own
coordinates, and none of them is the pool best. **`pool_best` is unchanged to the last
digit on all four.** So this is a disclosure requirement, not a retraction: the library is
not sequence-disjoint from the instrument on 3.2% of targets, and any retrieval-quality
claim should say so, but no reported number moves.

---

## B. Load-bearing constants — the audit table

Sensitivity measured in `s15/results/audit_sensitivity.json` on all 126 targets.
"Load-bearing" = a reasonable alternative moves a published number by more than its
reported uncertainty.

### B.1 The retrieval / selection / synthesis pipeline (`core.pipeline.Config`)

| constant | value | origin | fitted on | transfer | load-bearing? |
|---|---|---|---|---|---|
| `k` (BLOSUM retrieval depth) | 500 | swept in `s7.poolsize` (K ∈ 25…2000) on the **126-target tuning instrument** | the same 126 targets that report `pool_best` | never re-fitted per target | **YES for `pool_best`** (1.970 at K=100 → 1.504 at K=2000; 0.47 A of range). **NO for `shipped`** (3.425 / 3.461 / 3.454 / 3.483 / 3.520 over the same sweep — 0.095 A, flat) |
| `m` (distogram-score filter size) | 75 | swept in `s8.inband` / `s8.consensus2` over 8 sizes × 7 keys × 16 operators, on the 126 | the same 126 | a leave-fold-out arm-selection variant exists (`*_LFO`) but the **shipped 75 was chosen on all 126** | **YES**: `top75_best` 2.609 (m=25) → 2.106 (m=150). Any "top-75 ceiling" claim is a claim about m=75 |
| `penalty` | `"ramah"` | POST-HOC swap, declared in `s8/project.py` | — | — | **YES for structural-plausibility claims**, NO for RMSD (see B.5) |
| `lam` | 0.3 | geometric ladder `(0, .001, .003, .01, .03, .1, .3, 1.0)`, chosen at the m=75 tuning cell | the 126 | one pre-registered dev-24 pass | RMSD cost of the swap +0.004 A [-0.005, +0.013] — **not** load-bearing for RMSD |
| `maxiter` | 300 | stated parameter; `core.project iters` measures whether it binds | — | — | declared not tuned; `core/project.py` states the number of starts, `maxiter`, penalty, `lam` and tolerance are "science, not tuning knobs" |
| `multi_start` | True | S9-2: the projection is degenerate and a warm-started optimiser cannot cross branches | — | mechanism, not a fit | **YES** (documented) |
| `STARTS` | 4 conformations (-120,130), (-57,-47), (-139,135), (-75,145) | extended / α / β / PPII, copied by value from `s8.consensus2.FIT_STARTS` | theoretically motivated | — | probable, untested here |
| `project_grad` | `"exact"` | **not** a speed choice: `analytic` lands on different structures on 126/126 targets, median 0.042 A, worst 1.90 A | — | — | **YES**, and correctly hashed into the cache key |
| `reference_precision` | True | reproduces `s9.final`'s float32 npz round trip | — | — | YES, correctly hashed |
| `tie_break` | `"stable"` | PINNED; `np.argsort` kind | — | — | **YES.** BLOSUM `sim` has heavy ties, and the memory record already contains a case where `np.argmin` on a tied signal read the oracle sort order and invented a 1.386 A winner |
| `min_sep` | 2 | CA-CA pair separation for the distogram score | — | — | untestable from the cache: the cached distogram stores i,j at min\_sep=2 only, so 1 and 3 return `null` in the sweep. **Not audited. Flag.** |
| `amber_k` | 10.0 kcal/mol/Å² (`K_MODERATE`) | one of three declared levels (1 / 10 / 100) | — | Sprint 14 measured k=30 at -0.025 A [-0.046,-0.005] and k=300 at -0.105 A | **YES**: the restraint strength changes the refinement result and 10 is not the best measured value |
| `amber_steps` | 0 (= minimise to convergence) | OpenMM convention | — | — | YES |
| `torsion_window` | 8 | `torsion_lib2.library_for` | undocumented here | — | **UNDOCUMENTED. Flag.** |
| `n_folds` | 5 | leave-fold-out discipline | — | — | structural |
| `vqe_qubits/layers/iters/seed` | 7 / 3 / 50 / 0 | in the key but `quantum=False` in PROD | — | — | inert in the shipped path |

### B.2 The instrument's own constants (`s12/instrument.py`)

| constant | value | origin | load-bearing? |
|---|---|---|---|
| `K` | 500 | mirrors `Config.k` | see B.1 |
| `M` | 75 | mirrors `Config.m` | see B.1 |
| `BAND` | **1.5 Å** | **ARBITRARY. No derivation found anywhere in the repo.** | **YES, decisively — see B.6** |
| `PROD_KEY` | `"1fc9f2dcf489e2fb"` | verified equal to `Config().key()` | correct |
| `ALPHABET` | `"ARNDCQEGHILKMFPSTWYV"` | convention; note `core.data` also defines `ALPHABET_ALT` in a different order | a silent mix-up would be catastrophic and is not asserted anywhere |
| `paired(n_boot=4000, seed=0)` | | convention | no |

### B.3 Energy models

| constant | value | origin | load-bearing? |
|---|---|---|---|
| AMBER force field | `amber14/protein.ff14SB.xml` + `implicit/gbn2.xml` | literature | verified present in the built System (D.1) |
| nonbonded method | NoCutoff (method 0), cutoff field 1.0 nm unused | correct for implicit solvent | no |
| `K_WEAK/K_MODERATE/K_STRONG` | 1 / 10 / 100 | declared ladder | see `amber_k` |
| minimiser tolerance | 1.0 (kJ/mol/nm) | declared | untested; **flag** |
| strain rejection | `bond+angle > 1000` | declared, untested here | possible |
| `MEMO_SIZE` | 512 | not science (pure function memo) | no |
| Legacy `DEFAULT_WEIGHTS` | steric 4.0, contact 1.0, hbond\_local 1.0, hbond\_longrange 3.0, coop\_helix 2.0, coop\_sheet 2.0, solvation 0.5, electrostatic 1.0, aromatic 0.8, torsion 0.15, compactness 0.4 | **"variance-balanced starting point, not a fit"** (`energy_terms.py` docstring); `energy_quality.calibrate_weights` exists but was not used | **YES** — an eleven-weight vector that was never fitted is a strong claim to defend. Verified as the actual weights (D.3) |
| Legacy geometry constants | `AROM_PD_DIST 4.00`, `AROM_T_DIST 5.20`, depths 1.00/0.70, widths 1.10/35.0°, `AROM_CB_DIST 5.50/1.80`, `COOP_LADDER_MIN_SEP 2` | literature geometry + "modelling choice under the term's single calibrated weight" | the geometry is literature; the relative depths are declared arbitrary |
| MJ contact potential | Miyazawa–Jernigan table, verbatim | literature | no |

### B.4 Data-set construction

| constant | value | origin | leakage risk |
|---|---|---|---|
| `IDENTITY_THRESHOLD` | **0.6** | Needleman–Wunsch identity normalised by the longer sequence; used for clustering, folds, holdout | **KNOWN LEAKAGE PATH.** The project's own record ("A 0.6 containment threshold is AT THE NULL for peptides": random sequences score 0.56–0.63 against the fragment bank) applies to *containment*, a different statistic — but 0.6 is the same number reused for a different purpose with no separate justification. `core/data.py:31` also records that **573 library members across 58 targets contain a target at ≥0.6 while passing the 0.6 filter**. **Flag for the paper. This is the number a hostile referee will attack first.** |
| `MIN_LEN, MAX_LEN` | 8, 26 | dataset scope | no |
| `REBUILD_TOL` | 1.5 Å | "above this the deposited geometry carries something the representation cannot express" — a stated, principled criterion | probably not |
| CA-CA step filter | 3.5–4.1 Å | chain-break rejection | no |
| tuning-target length filter | 9 ≤ n ≤ 16 | instrument scope | no |
| `BENCHMARK_N` | 60 | pinned, inaccessible | n/a |
| dev set | 24 | pinned | n/a |
| `folds(5, seed=0)` | | pinned in `peptide_folds.json` | correctly pinned; the memory record already documents that re-deriving it silently moved 13 targets |

### B.5 The pre-registration mismatch — DEMONSTRATED, needs a Methods sentence

`s8/project_devarm.json` — the machine-readable pin, timestamped 2026-09-04 00:22:51 — says:

```json
{"pen": "phip", "lam": 0.03, "m": 75, "multi": true, "pass_number_on_family": 3, ...}
```

`core.pipeline.Config` ships `penalty="ramah", lam=0.3`. The swap **is** declared, in
`s8/project.py`'s docstring, in bold, with its criterion (*"THE ARM TO CARRY FORWARD IS
`ramah@0.3`, NOT THE ONE THE RULE PICKED… This choice is POST HOC on a criterion the
original rule did not contain, and is labelled as such"*) and its cost (+0.004 A
[-0.005, +0.013]).

The problem is that the *pin file* still records the superseded arm, so an automated
provenance check comparing pin to code reports a mismatch. The RMSD headline is unaffected
(the swap is worth 0.004 A). The **structural-plausibility** numbers are not: residues in
no canonical basin 2.90% (ramah) vs 14.34% (phip), symmetric KL 6.78 vs 8.80. Any paper
claim about backbone plausibility rests entirely on the post-hoc arm.

**Not a retraction — a disclosure requirement.** Methods must say: primary endpoint
pre-registered on `phip@0.03`; shipped arm `ramah@0.3` selected post hoc on a
distributional criterion; RMSD difference 0.004 A; plausibility claims derive from the
post-hoc arm.

### B.6 `BAND = 1.5 Å` and `FAIL18` — the weakest constant in the instrument

**DEMONSTRATED.** `n_zero_recall` is a pure function of an undocumented threshold:

| BAND (Å) | 0.5 | 0.75 | 1.0 | 1.25 | **1.5** | 2.0 | 3.0 |
|---|---|---|---|---|---|---|---|
| n\_zero\_recall | 45 | 34 | 24 | 20 | **18** | 13 | 2 |

And the *membership* is not stable either. `s12/results/adv_fail18.json` swept 99
(BAND, K, M) combinations: mean Jaccard against FAIL18 **0.498**, and

> **exactly ONE of the 18 targets (9KAR) is in the zero-recall set in all 99 combinations.**
> Five of the 18 appear in fewer than half. 56 non-members appear in at least one.

Sprint 12 measured this and the record contains it. But `I.FAIL18` is still a hard-coded
18-element list, `selfcheck()` still *asserts* it, and Sprints 12–14 use "the FAIL18
cohort" as though it were an object. **Every claim conditioned on FAIL18 is conditioned on
BAND = 1.5, which has no derivation.** For the paper, either derive BAND or report every
FAIL18 claim across the band sweep.

### B.7 Enumeration / quantum-architecture constants

| constant | value | origin | note |
|---|---|---|---|
| `K` (torsion states) | 4 | information-theoretic: power-of-two for the binary encoding | declared; 2 log2(k) dead qubits per chain |
| `N_AMBER_UNIFORM / PRIOR / BAND` | 1200 / 1200 / 600 | arbitrary budget | drives the 40/39/21% mixture of D.2 |
| band for `idx_b` | `max(2000, B//100)` (best 1% by TRUE RMSD) | **ORACLE by construction** | see D.2 |
| `BAND_FRAC` (`qarch_validity`) | 0.01 | ORACLE in-band definition | correctly labelled |
| `TOL` (`adv_locality`) | 1e-9 | machine-precision support threshold, replacing an earlier 0.10 A | an improvement, documented |
| CVaR `alpha` grids | (1.0, .5, .25, .1, .05, .01) / (1.0, .25, .05) / 0.25 | exploratory grids | the brief already records that `alpha` is not a learning rate |
| random-tail null | 0.524–0.527 | measured, **not** 0.500 | correctly documented in the brief |
| in-band accuracy required for 2.0 A | 0.638 | derived from the operator model | derived, not fitted |
| per-torsion sigma for 2.0 A | 15.1° (i.i.d.) | derived | the brief already notes coherent error needs 10.6° and biased 14.1° — **the 15.1° figure is conditional on i.i.d. and must never be quoted bare** |

---

## C. RMSD — independently reimplemented, and the frozen definition

### C.1 The independent implementation — DEMONSTRATED

`s15/audit_rmsd.py`, results in `s15/results/audit_rmsd.json`.

`s12.instrument.kabsch_rmsd_batch` uses SVD of the cross-covariance with the determinant
sign fix. The audit implementation (`horn_rmsd`) uses **Horn's quaternion method**: build
the symmetric 4×4 key matrix K from the 3×3 correlation matrix, take `lambda_max` from
`eigvalsh`, and `RMSD = sqrt((E0 - 2 lambda_max)/n)`. Different decomposition, different
LAPACK kernel, and unit quaternions cover SO(3) only — so an improper "fit" is unreachable
*by construction* rather than by a correction step. This is a mathematically independent
check, not a second copy.

**Agreement on the instrument's own targets: max absolute difference 2.04e-13 Å over
63,000 real structures** (126 targets × the K=500 pool, worst at 1S9Z). Against the stored
float32 `rr` labels the worst difference is 4.77e-7 Å, which is exactly the float32 storage
precision and not an error.

### C.2 Edge cases — all pass

| case | Horn | instrument | verdict |
|---|---|---|---|
| identity | 0.0 | 0.0 | pass |
| rotation + translation | 0.0 | 0.0 | pass — rigid-motion invariant |
| **z-mirrored copy** | 2.9601965330 | 2.9601965330 | **PASS — improper rotation correctly rejected.** The same test with the determinant fix removed returns **0.0**, i.e. a mirror image would be scored as a perfect fit |
| point inversion (x → −x) | 2.9601965330 | 2.9601965330 | pass, same mechanism |
| n = 1 | 0.0 | 0.0 | pass |
| n = 2 | 0.0 | 0.0 | pass |
| collinear (rank-1) | 0.0 | 0.0 | pass — degenerate rotation handled |
| all atoms coincident (rank-0) | 0.0 | 0.0 | pass |
| NaN coordinate (missing atom) | `LinAlgError` | `LinAlgError` | **pass** — it raises, it does not return a plausible float |
| length mismatch (10 vs 14) | `ValueError` | `ValueError` | pass |
| chain break (25 Å gap) | 20.3281 | 20.3281 | finite and large — the metric scores it, it does not reject it. **This is a definitional choice and must be stated** |
| coordinate doubling | 10.2690 | 10.2690 | RMSD is in Å and is not scale-invariant — correct |
| residue-order swap | 1.43626 | 1.43626 | correspondence is by INDEX, never nearest-neighbour — correct |
| float32 vs float64 input | 3.44205711 | 3.44205719 | 8e-8 Å; the universes store `W` as float32, so all pool RMSDs carry ~1e-7 Å of storage noise |
| terminal residues wrecked | **full 3.4421** | | **trimmed [1:-1] 1.0257** — a 2.4 Å difference on the same structure |

### C.3 THE FROZEN DEFINITION

> **CA-RMSD, as used everywhere in this project and in any paper arising from it.**
>
> 1. **Atom set.** Alpha carbons only, one per residue, for **every residue of the chain,
>    including both termini**. No trimming, no core selection, no per-residue outlier
>    rejection, no distance cutoff, no iterative superposition.
> 2. **Correspondence.** By residue index, `i ↔ i`, over two arrays of identical length n.
>    Never by nearest neighbour, never by alignment, never by sequence matching. A length
>    mismatch is an error, not a truncation.
> 3. **Weighting.** Uniform, weight 1/n per residue. (`core.project.wkabsch_rmsd_batch`
>    supports non-uniform weights for the projection's internal use; it is numerically
>    identical to the unweighted routine at uniform weights, and **no reported RMSD uses
>    non-uniform weights**.)
> 4. **Alignment.** Optimal superposition over the **proper** rotation group SO(3) plus
>    translation: centre both structures on their unweighted centroids, then minimise over
>    rotations of determinant **+1** only. **Reflections and improper rotations are
>    forbidden**, so an enantiomer of the native does not score as the native. Verified
>    two ways (C.2).
> 5. **Units.** Ångström. Not scale-invariant, not normalised by chain length beyond the
>    1/n inside the root mean square.
> 6. **Formula.** `RMSD = sqrt( (1/n) * sum_i |R x_i + t - y_i|^2 )` at the minimising
>    `(R, t)`, `det R = +1`.
> 7. **Degenerate inputs.** n = 1, n = 2 and collinear chains return 0 (the fit is exact or
>    the rotation is non-unique). n = 0 is an error.
> 8. **Missing atoms.** Not supported. A non-finite coordinate raises. Structures with
>    missing backbone are excluded upstream by `peptide_db._scan` (CA count must equal
>    sequence length; CA-CA steps must lie in 3.5–4.1 Å; the torsion rebuild must be within
>    `REBUILD_TOL = 1.5` Å).
> 9. **Alternate conformations / multi-model NMR.** The native is **model 1** of the
>    deposited file, taken once (`s12/instrument.py` docstring: "CA-RMSD of each window to
>    the native (model-1) CA trace"). No ensemble minimum, no best-of-models.
> 10. **Chain breaks.** Scored, not rejected. The metric is defined on any two equal-length
>     point sets.
> 11. **Precision.** Reference arrays are float64. Library windows are stored float32, so
>     pool RMSDs carry ~1e-7 Å of storage noise. Two independent implementations agree to
>     2.0e-13 Å.
>
> **This is the full-chain definition. The September 2026 preprint audited in Sprint 14
> reported an RMSD that excluded terminal residues and never defined it in Methods; on this
> instrument that choice is worth up to 2.4 Å on a single structure (C.2). Any number in
> our paper that is not full-chain must say so at the point of use.**

---

## D. Energy models

### D.1 AMBER provenance — DEMONSTRATED

`s15/results/audit_energy_amber.json`, plus `verify/amber_audit.py` rerun today.

* Force field: `amber14/protein.ff14SB.xml` + `implicit/gbn2.xml`, both present in source.
* Built System for 1A13 (n=14): **237 particles**, forces = HarmonicBondForce,
  PeriodicTorsionForce, NonbondedForce, HarmonicAngleForce, **CustomGBForce**,
  CustomExternalForce (the positional restraint). GBn2 is genuinely present.
* NonbondedMethod = 0 (NoCutoff), correct for implicit solvent.
* Total charge +3.000 e, 105 distinct partial charges, 7 distinct GB radii.
* **Parameters verified against an independently constructed ForceField**: max abs charge
  diff 0.0, max abs sigma diff 0.0, max abs epsilon diff 0.0, GB params max abs diff 0.0.
* **Pinned golden reproduces bit-exactly today**: 1A13 native interaction energy
  (nonbonded + solvation, k=10, steps=0, tol=1.0) = **-489.9138948277905 kcal/mol**,
  `abs_diff 0.0`.
* Single point is deterministic across a full cache teardown.

**One misleading field in an existing audit — flag.** `verify/amber_audit.json` reports
`"5_translation_invariant": 0.5446` under a comment reading *"a rigid translation must
leave the energy invariant: real physics, not a hash"*. It does not measure that: it
translates and then **minimises to convergence** (`steps=0`), so it measures the
minimiser's basin sensitivity, not the physics. Separating the two (`s15/results/
audit_amber_translation.json`):

| | solvation Δ | total Δ |
|---|---|---|
| single point (`steps=-1`) | **1.28e-6 kcal/mol** | 2.69e-5 |
| minimised (`steps=0`) | 0.148 | 0.041 |

The physics is translation invariant to 1e-6. The *protocol* is not, to ~0.1 kcal/mol.
The field name asserts the first and reports the second, with no PASS/FAIL — precisely the
"reading past a check" trap the brief warns about. **Rename or split it.**

### D.2 The oracle-conditioned AMBER subset — VERIFIED, and the recorded figure is CORRECTED

`s15/results/audit_amber_oracle.json`, all **19** enumerated targets
(9 × `s13/results/qarch_enum_*.npz` at k=4 n=9, 10 × `s14/cache/obj_enum_*.npz` at n=10).

The mechanism is confirmed by code reading and by measurement. In
`s13/qarch_enum.py`:

```python
band  = np.argsort(rmsd)[:max(2000, B // 100)]     # best 1% by TRUE CA-RMSD
idx_b = rng.choice(band, N_AMBER_BAND, replace=False)
a_idx = np.unique(np.concatenate([idx_u, idx_p, idx_b, [snap_idx]]))
a_kind = np.zeros(len(a_idx), np.int8)
a_kind[np.isin(a_idx, idx_p)] = 1
a_kind[np.isin(a_idx, idx_b)] = 2
```

**CONFIRMED — the effect size.** Mean over 19 files, subset mean RMSD minus space mean
RMSD = **−0.4243 Å** (sd 0.101). On the 9 Sprint 13 files alone it is **−0.4009 Å**.
Sprint 14's recorded "0.401 A better than its space" is exact to three decimals.

**CONFIRMED — the operative instruction.** `amber_kind == 0` is an unbiased sample of the
space: mean delta **+0.0003 Å** over 19 files (+0.008 Å over the 9), and its fraction
falling inside the oracle band is 0.00973 against the uniform expectation 0.00999.
**Use `amber_kind == 0`. That advice stands.**

**CORRECTED — the "40%" is transposed.** The oracle-conditioned rows are `kind == 2`, and
they are **21.7% of the subset** (range 20.2–23.2%), not 40%. The **40%** figure is the
share of `kind == 0`, i.e. the **unbiased** rows (range 39.0–42.5%). The remaining ~39% is
`kind == 1`, prior-sampled — native-free, but not uniform (mean delta +0.076 Å).

> **RETRACTION / CORRECTION.** The sentence "the cached AMBER subset is 40%
> ORACLE-CONDITIONED" is wrong on the number. The correct statement is:
> **"only 40% of the cached AMBER subset is unbiased; 22% is oracle-conditioned and a
> further 39% is prior-conditioned; the whole subset sits 0.42 A below its own space."**
> The conclusion drawn from it — use `amber_kind == 0` — is unaffected and correct.
> Preserve the original claim beside this correction.

**NEW DEFECT — `amber_kind == 0` is unbiased on the MEAN but NOT on the MINIMUM.**
`a_idx` force-includes `snap_idx = ORACLE_snap()`, the k=4 grid configuration nearest the
native. It is not in `idx_p` or `idx_b` on most files, so it is labelled `kind = 0`:

* the ORACLE snap lands in `kind == 0` on **16 of 19 files**;
* on **6 of 19 files (2MK7, 6S0N, 8IS3, 2MJQ, 7VI4, 8HVS) the minimum-RMSD member of the
  "unbiased" `kind == 0` subset IS that oracle snap.**

Diluted over ~1200 rows this is worth ≲0.003 Å on a mean, which is why it never showed up.
It is decisive for anything that reads a minimum, a tail, a top-k or an argmin — which is
what every ranking statistic in Sprints 13–14 does.

> **BINDING RULE for Sprint 15: any tail-restricted, in-band, argmin, best-of-k or
> top-percentile statistic on the cached AMBER subset must use `amber_kind == 0 AND
> amber_idx != snap_index`. Mean-level statistics may use `amber_kind == 0` as recorded.**

**Exposure, quantified** (`s15/results/audit_amber_snap_exposure.json`,
`audit_amber_snap_inband.json`):

| statistic | files exposed |
|---|---|
| snap in the RMSD-lowest 1% of its own `kind==0` stratum | **10 of 16** (median snap RMSD percentile **0.14**) |
| snap in the AMBER-lowest 1% of its stratum | 4 of 16 |
| snap in the AMBER-lowest 0.2% (the smallest `ener_tail` quantile) | 2 of 16 |
| **snap inside the `<1.5 A` IN-BAND set** | **14 of 19** |

The in-band exposure is the serious one, because the in-band set is small:

| target | in-band (<1.5 A) members in `kind==0` | share of in-band PAIRS touching the snap |
|---|---|---|
| 2MJQ | **5** | **40.0%** |
| 7VI4 | **8** | **25.0%** |
| 2MK7 | **12** | **16.7%** |
| 8HVS | **13** | **15.4%** |
| 5V5B | **19** | **10.5%** |
| median over the 14 contaminated files | 157 | 1.35% |

**Two consequences for the headline in-band accuracy table (Legacy 0.539 / learned 0.523 /
prior 0.501 / AMBER 0.463):**

1. It is computed on a band that carries an ORACLE-inserted member on 14/19 targets, and on
   five of those the member touches 10–40% of every in-band pair. The number is not clean.
   The direction is safe — the snap is near-native, so its presence can only *flatter* an
   objective — which means **AMBER's true in-band accuracy is at most 0.463 and plausibly
   lower**, and the published negative is if anything understated. But it is not exact.
2. Independently of the snap, **the in-band set has 5, 8, 12, 13 and 19 members on five
   targets** — 10, 28, 66, 78 and 171 pairs respectively. Per-target in-band accuracy on
   those is noise, and the aggregate depends entirely on how the targets are weighted.
   The project's own trap list already records that "n ≤ 4 samples from the enumerated nine
   have reversed a conclusion four times". **How the 0.463/0.539 aggregate is weighted must
   be stated, and per-target n must be printed beside it.**

**Fairness to the Sprint 14 record.** `s14/ener_lib.py:21` labels the strata correctly
("0=uniform 1=prior 2=ORACLE band"), and `s14/obj_FINDINGS.md` §1.3 states the correction
accurately with the right offsets (−0.405 A full, +0.007 A on `kind==0`, 978 rows/target —
all reproduced here). **The transposition is in the brief's one-line summary, not in the
Sprint 14 modules or findings.** The snap-index defect is new and is in neither.

### D.2b Which published numbers used the FULL (contaminated) subset — one retraction

Every consumer of `amber_total` in s13/s14 was read. Result:

| module | stratum used | verdict |
|---|---|---|
| `s13/qarch_robust.py` | `a_kind == 0` | clean |
| `s13/qarch_sepfit.py` | `a_kind == 0` | clean |
| `s13/qarch_validity.py` | all three strata reported separately, and it already tracks the snap's index within each stratum (`si`) | clean, and the most careful of the set |
| `s13/qarch_locality.py` | computes AMBER fresh on its own configs, does not read the cache | clean |
| `s14/ener_{complement,decoy,geom,matrix,norm,tail}.py` | `z.uniform_mask` | clean on the mean; **exposed to the snap defect** (D.2) |
| `s14/vqe_collapse.py` | explicit `masked` flag, defaults to the mask, and deliberately reproduces the unmasked version so the difference can be reported | exemplary |
| **`s14/obj_floor.py`** | **the FULL `amber_idx` subsample, no `amber_kind` mask** | **CONTAMINATED** |

> **RETRACTION.** `s14/results/obj_floor.json` computes every AMBER and `amb_*` row on the
> full labelled subsample. The file contains **no `amber_kind` or `uniform` string at all**,
> carries `"complete": true`, and covers `n_expected: 9` targets. Its AMBER rows are
> conditioned on the native by 0.42 Å and its target count is superseded by the 16- and
> 19-target passes.
>
> **Every AMBER row in `s14/results/obj_floor.json` is retracted.** The corrected values
> exist — in `s14/obj_FINDINGS.md` §1.3, on 16 targets, `kind == 0` — but *the artefact was
> never regenerated*, so a reader who opens the JSON rather than the prose gets the
> superseded numbers with a green completeness flag on top. That is the "reading past a
> check" trap, live in the repository today.
>
> **Action: regenerate with the mask, or rename to `obj_floor.STALE-UNMASKED-AMBER.json`.**
> The Legacy, prior and `leg_*` rows in that file are computed over the full enumeration and
> are unaffected.

One minor related item: `s13/qarch_validity.py` computes its AMBER normalisation
`amber_mu` / `amber_sd` from `a_kind == 1` (the **prior-sampled** stratum), not from
`kind == 0`. That stratum sits +0.076 Å from the space mean. Small, but it is a
normalisation derived from a conditioned population and should be stated or changed.

### D.3 Legacy — independently verified as an algebraic identity

`s15/results/audit_energy_legacy.json`.

* The stored `legacy` column is **not** the plain sum of the eleven `leg_*` components
  (max abs discrepancy 111–266 kcal/mol against a total sd of 6–20).
* It **is** the `energy_terms.DEFAULT_WEIGHTS`-weighted sum, verified independently on all
  **19 files / 1.28e7 structures**: worst `|legacy − Σ w_t · leg_t|` = **2.71e-5**, which
  is float32 storage precision. DEMONSTRATED.
* Recomputing Legacy from source through `qarch_lib.legacy_components` for 64 configs of
  2MK7 reproduces the stored column **bit-identically in float32** (max float64 diff
  3.9e-7).

### D.4 AMBER cost — the 28 ms figure is CORRECTED

`s15/results/audit_energy_cost.json`, `audit_amber_cost_sweep.json`.
Distinct inputs, memo defeated and counted, `threads=1`, CPU platform, ff14SB/GBn2 single
point (`steps=-1`, no minimisation):

| target | n | atoms | ms / distinct call |
|---|---|---|---|
| 1CS9 | 9 | 127 | **8.3** |
| 9UV5 | 9 | 136 | **9.7** |
| 2MK7 | 9 | 145 | **10.8** |
| 1N9U | 10 | 182 | **14.1** |
| 8HVS | 10 | 213 | **19.7** |
| 1A13 | 14 | 237 | **23.3** |

Memo hits during the distinct loop: 1 of 30 (an unavoidable collision), 29 misses. The same
structure repeated 30 times costs **0.33 ms/call, 30/30 memo hits — a 33× speedup.** That
is the mechanism behind the original wrong 6 ms figure, confirmed directly.

> **The brief's "AMBER single point 28 ms distinct" is a ceiling for the largest targets,
> not a constant.** The true cost is **8–23 ms**, scaling with atom count (GBn2 is the ~98%
> term and is O(N²)). On the n=9–10 targets where all the enumerated AMBER data actually
> lives it is **8–14 ms**. Budget-matching arguments that priced AMBER at 28 ms on those
> targets over-charged it by 2–3×. `s13/qarch_lib.amber_energies`'s docstring still says
> "~6 ms per configuration warm" — that line is **stale and should be corrected in place.**

### D.4b CVaR — independently verified (a non-negotiable component, so checked)

`s15/results/audit_cvar.json`. `core.quantum.cvar_exact` computes the lower-tail CVaR by
sorting, accumulating mass and adding the boundary state's partial contribution. It is
checked here against the **Rockafellar–Uryasev variational form**,
`CVaR_a(X) = max_t [ t − (1/a) E(t − X)^+ ]`, evaluated exactly on the discrete support —
a different formula, not a reimplementation of the same one.

* 400 random cases, n ∈ [3, 60], energies scaled up to ×50, `alpha` ∈ {0.01 … 1.0},
  every fifth case with deliberate ties: **max absolute value error 4.17e-14**
  (at |CVaR| ≈ 6.2, i.e. ~1e-15 relative — consistent with the recorded 1.7e-15).
* The envelope-theorem derivative `dCVaR/dp` checked against finite differences along
  mass-preserving directions (`p_i += h, p_j -= h`, h = 1e-7), 398 cases:
  **median error exactly 0.0**, p95 5.85e-8 (finite-difference noise).

DEMONSTRATED. The value estimator and its exact `dp` are correct. The three recorded CVaR
*defects* (the `baseline="tail"` gradient bias, the sampled upward bias at non-integer
`alpha*N`, and `dCVaR/dp ≡ 0` iff `p(argmin E) ≥ alpha`) are properties of the *estimators
and gradients around* this function, not of it, and are not re-tested here.

### D.5 Raw `amber_total` is unusable without conditioning — DEMONSTRATED

`s15/results/audit_amber_scale.json`, 19 files.

* **37.6%** of cached AMBER single points exceed 1e6 kcal/mol; **3.2%** exceed 1e12; the
  largest is **2.65e20 kcal/mol**.
* Per-file mean ~1e15–1e16 with sd ~1e17, against a median of ~1e4 and a minimum of
  −510 kcal/mol.

Any mean, sd, z-score, correlation or normalisation taken on the raw column is a statistic
about the single worst steric clash. This is why `s14/cache/ener_norm.json`'s `grad_sd`
exists and why the brief says never to use raw AMBER's sd; the audit confirms the magnitude
of the problem. **A 99th-percentile winsorisation is not enough — the 99th percentile is
itself 1e12–1e15 on most files.**

---

## E. Provenance graph

Machine-checkable record: `s15/results/audit_provenance.json`
(`python -m s15.audit_provenance` builds it, `... check` re-verifies every sampled
sha256 against it). It carries, per node: kind, glob, file count, total bytes, mtime range,
sha256-16 of up to 8 representative files, the code that writes it, whether that writer is
config-keyed, whether the write is atomic, and which of the five pinned constants depends
on it. It also carries a `constant_dependencies` block giving, for each pinned constant,
its value, the exact expression, the artefacts it reads, whether the instrument recomputes
it, whether it has been reproduced from primary inputs, the evidence file, and whether it
is ORACLE.

The chain, end to end:

```
pdbs/ + pdbs_ext/ (1524 .pdb)  --geo.native_coords_from_pdb + filters-->  peptide_db.npz  [BIT-IDENTICAL, A2.1]
prots/ (13751 .pdb)            --distogram._fold_fragments-------------->  fold fragment libraries
peptide_db.npz --identity 0.6--> peptide_clusters.json --seed 0--> peptide_folds.json   [PINNED]
  |                                                                          |
  +--> debias.tuning_targets() -> the 126 --------------------------+        |
                                                                    v        v
  library (out-of-fold peptides + this fold's fragments) --windows--> s8/generate_univ/<pdb>.npz
                                                                       [BIT-IDENTICAL 126/126, A2.2]
                                                                       W PHI PSI S org sim order | rr nat_ca (ORACLE)
                                                                          |
   ESM --> distogram_models/fold<k> --> s12/cache/disto_<pdb>.npz --------+
           [thread-sensitive, A0/A2.5]                                    |
                                                                          v
   order[:500] -> pool -> shipped_score -> argsort[:75] == production `sub` (126/126, A2.3)
                                            |
                        coordinate_average -+-> core.project (ramah, lam .3, 4 starts, maxiter 300)
                                                 -> fit_ca ---> synthesis_fit 3.2041  [3/126 rerun bit-identical, A2.4]
                                                 -> ca/phi/psi -> core.amber (ff14SB/GBn2, k=10, steps=0)
                                                                 -> amber_ca  [bit-identical, A2.4]
   all of the above cached under bench_results/cache/1fc9f2dcf489e2fb/  [key == Config().key(), verified]
```

### Traceability spot-check of Sprint 14's flagship result — PASSES

LEDGER C11, "the ablation ladder closes: 4.072 → 3.485 (−0.587, aggregation) → 3.314
(−0.171, objective) vs incumbent 3.204". Every rung traced to a stored float:

| rung | value | artefact | path inside it |
|---|---|---|---|
| A one committed torsion vector | 4.072482457887866 | `s14/results/ladder.json` | `rows/L2b_top75_circmean/mean` |
| B + sample 200, coordinate-average | 3.485128 | `s14/results/vqe_classical_limit.json` | `curve/200/mean` |
| C + rank 4000, average top 200 | 3.314163 | `s14/results/consensus.json` | `aggregate/200/mean` |
| incumbent | 3.2040761603809194 | `s14/results/budgetcurve.json`, `ladder.json` | `incumbent`, `rows/incumbent_synthesis/mean` |

All four agree to the quoted precision, and the incumbent matches the instrument's
`synthesis_fit` to the last digit. **The headline is traceable.** The only weakness is
presentational: the four rungs live in four different files written by four modules, one of
them named `vqe_classical_limit.json`, and no artefact contains the ladder as such. A
referee asking "where is Table 1?" has to be walked through four files.

### Numbers I could NOT trace, and would not put in a paper as they stand

1. **`min_sep = 2`.** The cached distograms store `i, j` at min\_sep = 2 only, so no
   alternative is testable from the cache and no sweep was found. The pair set the entire
   selection score is computed on is undocumented and unaudited.
2. **`torsion_window = 8`** (`Config`, hashed into the key). No derivation found.
3. **The minimiser tolerance 1.0** and the **`bond+angle > 1000` strain rejection**. Both
   declared in `core/amber.py`, neither swept.
4. **`s12/newpdbs/` (133 files).** Not in `core.data.DIRS`, not in `peptide_db`, purpose
   undocumented.
5. **Five duplicated Sprint 14 result files** with no canonical designation (A3).
6. **Which Sprint 13/14 published numbers took an argmin or a tail over `kind == 0`**, and
   are therefore exposed to the snap-index defect of D.2. Not enumerated. First follow-up.

---

## Z. Reproduction commands

```
set PYTHONIOENCODING=utf-8
set OMP_NUM_THREADS=2 & set MKL_NUM_THREADS=2 & set OPENBLAS_NUM_THREADS=2

python -m s12.instrument                       # the five pinned constants, 2.2 s
python -m s15.audit_repro db                   # rescan 1524 PDBs -> peptide_db.npz   (45 s)
python -m s15.audit_repro univ all             # rebuild all 126 universes            (31 s)
python -m s15.audit_repro disto 1A13,1A1P      # recompute distograms from the model
python -m s15.audit_determinism 1|2|8          # shipped, recomputed at N torch threads
python -m s15.audit_pipe run 1A13,1A1P,2MK7    # rerun the production pipeline        (31 s)
python -m s15.audit_rmsd                       # Horn quaternion RMSD + 15 edge cases
python -m s15.audit_energy amber|cost|legacy   # force-field provenance, cost, algebra
python verify/amber_audit.py                   # the pinned -489.9138948277905 golden
python -m s15.audit_amber                      # the oracle-conditioned subset, 19 files
python -m s15.audit_sensitivity                # K / m / BAND sweeps, all 126          (~5 min)
python -m s15.audit_leak                       # leave-fold-out + window-level sequence check
python -m s15.audit_provenance                 # build the provenance record
python -m s15.audit_provenance check           # re-verify every hash in it
```

---

## Y. Recommendations, in priority order

1. **Regenerate or rename `s14/results/obj_floor.json`** (D.2b). A stale artefact with a
   green completeness flag is worse than a missing one.
2. **Adopt the binding rule** `amber_kind == 0 AND amber_idx != snap_index` for every tail,
   in-band, argmin or top-k statistic on the cached AMBER subset (D.2), and re-run the
   in-band accuracy table with per-target n printed.
3. **Assert all five constants in `s12.instrument.selfcheck`.** It currently asserts only
   `pool_best`, `top75_best` and the FAIL18 set; `shipped` and `synthesis_fit` are printed
   and unchecked. I deliberately did **not** patch shared machinery mid-sprint — other
   workstreams are running against it — so this is a recommendation, not a change.
4. **Fix `I.write`'s completeness heuristic** to take an explicit row count instead of
   `len(obj["rows"])`. Its current 12/12 false-alarm rate makes the flag noise (A3).
5. **Finish the 126-target `synthesis_fit` reproduction** (`python -m s15.audit_pipe run
   <all 126>`, ~90 min single process). 22/126 complete and bit-identical at the time of
   writing; the module writes incrementally, so it can simply be re-run.
6. **Freeze one canonical location** for the five duplicated Sprint 14 results (A3).
7. **Rename `verify/amber_audit.json`'s `5_translation_invariant`** — it measures minimiser
   basin sensitivity, not physics invariance, and the physics is invariant to 1.3e-6 (D.1).
8. **Correct in place**: `s13/qarch_lib.amber_energies`'s "~6 ms per configuration warm"
   docstring (true cost 8–23 ms distinct, 0.33 ms on a memo hit) and the brief's flat 28 ms.
9. **Derive or retire `BAND = 1.5`** and stop treating `FAIL18` as an object (B.6).
10. **Add an input hash to `s8/generate_univ/*.npz` and `s12/cache/disto_*.npz`.** Both are
    keyed on PDB id alone and neither would notice a changed model, library or fold file.

---

## PHASE 0 GATE — verdict

**PASS, with four disclosures and one correction.**

Reproduces exactly from primary inputs: `peptide_db.npz` (bit-identical, 787/787), all 126
window universes (bit-identical, 9/9 arrays × 126), `pool_best`, `top75_best`,
`n_zero_recall`, `shipped` (at three thread counts), the AMBER golden, and the production
pipeline end to end on 3/126 targets including bit-identical AMBER minimisation.

Does not yet reproduce: `synthesis_fit` on 123 of 126 targets (not run, ~23 min).

Corrected: the "40% oracle-conditioned" figure (it is 22%; 40% is the *unbiased* share).

New defect: `amber_kind == 0` contains the ORACLE snap on 16/19 files and it is the subset
minimum on 6/19.

Disclosures required in Methods: the pre-registration/shipped arm mismatch (B.5); `BAND`
and therefore `FAIL18` being threshold artefacts (B.6); the 0.6 identity threshold (B.4);
the RMSD definition (C.3).
