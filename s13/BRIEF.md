# Sprint 13 — shared agent brief: TORSION-CONSTRAINED VQE/CVaR

READ FULLY. You are one of several agents on a sprint whose goal is a genuinely new
architecture **and** a rigorous variational-trainability study. Obey every rule here.

## 0. What this sprint is

Two deliverables, equally weighted:

1. **A torsion-constrained VQE/CVaR peptide-folding architecture** with a credible path to
   < 2.0 Å mean CA-RMSD, keeping genuine Legacy and genuine AMBER ff14SB/GBn2.
2. **A controlled study of how molecular energy-model complexity reshapes the geometry and
   trainability of the variational quantum problem.** Legacy vs AMBER, everything else
   matched: same representation, ansatz, initialisation, CVaR α, optimiser, seed, and
   **objective-evaluation budget** (not wall time — AMBER is ~20× costlier per call).

The audience is a professor working on variational quantum algorithms, barren plateaus,
quantum information geometry, and trainability. The protein problem is the laboratory; the
VQE is the object of study.

## 1. Hard rules (non-negotiable)

1. **benchmark60 is PROTECTED.** Never read `results/benchmark_manifest.json` targets, never
   open `s9/final_cache/*`, never tune on it. One pre-registered pass at the very end, only
   if the coordinator authorises it.
2. **dev24 is reserved.** Do not run on it without the coordinator's say-so.
3. **NO NATIVE INFORMATION AT INFERENCE.** Native coordinates, torsions, distances, RMSD,
   energies, contact maps and secondary structure may be used ONLY for (a) training a
   supervised predictor under leave-fold-out discipline, (b) explicitly labelled ORACLE
   diagnostics, (c) post-hoc evaluation. Name every oracle arm `o_*` or `ORACLE_*`.
4. **Tier every claim** in your findings: DEMONSTRATED (measured on held-out targets) /
   ORACLE DIAGNOSTIC (needs native info, not a predictive result) / HYPOTHESIS /
   LITERATURE-SUPPORTED / PROPOSED. The most dangerous failure mode in this repository is an
   oracle diagnostic quietly becoming a headline.
5. **Do not modify** `core/`, root modules, `s5/ s7/ s8/ s9/ s12/`, `tests/`, `verify/`,
   `pdbs/`, `prots/`, `pdbs_ext/`, `peptide_*.json`, root `*.npz`, `bench_results/`. No git
   writes. Your code goes in `s13/<yourname>_*.py`, results in `s13/results/<yourname>_*.json`,
   findings in `s13/<yourname>_FINDINGS.md`.
6. **No fake components.** Do not call a surrogate "AMBER". Do not call enumeration "VQE".
   Do not call a classical optimiser a quantum result. If you use an AMBER surrogate, name it
   a surrogate everywhere.
7. **COMPUTE, and this is live:** the user is in a meeting for roughly the first two hours, so
   until the coordinator says otherwise keep to **≤ 2 concurrent heavy processes and ≤ 1.2 GB
   per process**. Call `s12.instrument.free_gb()` before every heavy step and wait if under
   1.5 GB. `OMP_NUM_THREADS=2`. After the coordinator lifts the cap, target ~95% CPU/RAM but
   never exceed 97%.
8. **Statistics on every headline**: n, mean, median, sd, paired comparison, bootstrap CI,
   W/L, per-fold, drop-top-10/20, fraction under 2.0 Å and under 1.5 Å, FAIL18 vs other-108.
   A mean gain carried by two targets is not a result.
9. **Budget parity is a correctness condition.** Track exact objective-evaluation counts.
   `budget.BudgetedEnergyModel` does this; use it. A VQE that quietly got fewer evaluations
   than its classical control is a broken experiment, not a negative result.

## 2. What already exists — do not rebuild it

The repository **already contains** the torsion-space VQE machinery. Sprint 12 concluded the
retrieval architecture is information-limited; the torsion oracle says the representation the
project abandoned was the right one. Read these before writing anything:

- `torsion_lib2.library_for(seq, k, exclude_seq=seq)` → a **sequence-conditioned, leakage-safe**
  torsion library: k states per residue. `PerResidueTorsion(seq, tab, chi_bits=False)` wraps it
  and reports `n_bits = n_res * log2(k)`.
- `core.quantum` (= `foldvqe` + `qansatz` + `vqe` + `objective` + `hamiltonian` consolidated):
  `StatevectorCircuit` (exact, n ≲ 20), `MPSAnsatz` (exact, no qubit ceiling, bond dim capped
  at 2^layers by the chain topology), PennyLane `lightning.qubit` (n ≤ 30),
  `run_global_cvar_vqe`, `run_cvar_vqe`, `cvar`, `cvar_gradient`, `FoldObjective`,
  `FoldingHamiltonian`, `AmberHamiltonian`, `_spsa`, `Reservoir`, `fold`.
- `core.energy` — the genuine Legacy 11-term knowledge-based model. `energy_components` /
  `components_batch` (0.3 ms per structure).
- `core.amber` — genuine ff14SB/GBn2. **CORRECTED TIMING (an earlier figure of 6 ms in this
  brief was wrong — it hit the result memo by re-evaluating the SAME state).** Measured on
  DISTINCT configurations, one warm builder, threads=1:

  | call | cost | note |
  |---|---|---|
  | `single_point`, distinct states | **28 ms** | the real per-evaluation cost |
  | `single_point`, repeated state | 0.8 ms | memo hit — never time a loop this way |
  | `refine_coords(k=10, steps=0)` | **9.1 s** | full restrained minimisation |
  | `core.energy.components_batch` (Legacy) | **0.3 ms** | |

  So **AMBER is ~90× Legacy per objective evaluation.** A 19,200-evaluation VQE run costs
  ~9 minutes per target under AMBER and ~6 seconds under Legacy. That is affordable, and it
  is exactly why budget parity must be counted in EVALUATIONS, never in wall time.
- `s12/instrument.py` — the shared measurement instrument: `targets()` (126 tuning targets),
  `load_univ`, `ca_rmsd`, `kabsch_rmsd_batch`, `build_ca(phi, psi)` (accepts `(n,)` or `(B,n)`),
  `project`, `paired`, `summary`, `write`, `FAIL18`, `free_gb`.

**THE CVaR GRADIENT DEFECT MUST NOT COME BACK.** `core.quantum.cvar_gradient` defaults to
`baseline="const"`, which is correct. `baseline="tail"` reproduces a measured bias (cosine
+0.58 with the exact gradient at 0.56× norm). Never use `"tail"` for a result.

**AMBER's 92% MEMORY GUARD.** `core.amber.builder_for` raises `MemoryError` above 92%
physical memory and the agent fleet trips this routinely. Poll `core.amber.memory_percent()`
and WAIT for headroom. Do not catch the exception and continue — a NaN column silently
becomes "AMBER has no signal", which would be a fabricated result.

## 3. The empirical record you are building on

Sprint 12 (`s12/SPRINT12_DOSSIER.md`, 898 lines) established, on the 126-target tuning set:

| quantity | value |
|---|---|
| shipped retrieval pipeline | 3.213 Å (FAIL18 6.019) |
| learned set decoder over the full signed deviation map | flat learning curve — in-band selection is **signal-limited** |
| perfect distance objective, m=75 / co-optimised m=5 | 2.395 / ~1.99 Å |
| whole architecture's sequence information | **0.776 Å** total; **negative on FAIL18** (blind 5.425 beats shipped 6.019) |
| **ORACLE continuous torsions, σ=12°, full coverage** | **1.486 Å, 86% under 2 Å, FAIL18 1.80** |
| same torsions used as a pool *filter* | 2.972 Å — **the value is in BUILDING, not filtering** |
| coverage sensitivity at σ=12° | 90% → 2.14, 75% → 2.79, 50% → 3.29 (worse than incumbent) |
| terminal operator law | output ≈ 1.16·(set mean) + 0.04·(set best); a perfect rank-1 is worth 1.74 Å through argmin and 0.03 Å through the shipped average |
| a leave-fold-out **4-state torsion-bin predictor** | 0.690 accuracy overall, **0.517 on FAIL18** (majority baseline 0.562) — it failed as a retrieval key |
| 10 of the 18 failures are fibril or lasso peptides | Fisher p = 1.2e-06; lasso is unreachable by linear-window retrieval |
| fresh benchmark supply | 16 targets, 10 of them amyloid fibrils — **no adequately-powered fresh instrument exists** |

**Representation ceiling measured at the start of this sprint** (`s13/ceiling.py`, ORACLE):
on 1A13, ideal-geometry rebuild of the native torsions is 0.284 Å; the discrete library
reaches 1.89 Å at k=4 (28 qubits), 1.63 at k=8 (42), 0.89 at k=16 (56). The full 126-target
sweep is running; read `s13/results/ceiling_report.json` before assuming any ceiling.

## 4. The central design tension you must respect

In torsion space the chain builds sequentially, so **a distance between residues i and j
depends on every torsion between them**. A pairwise distance term is therefore NOT 2-local in
the torsion variables — it is |i−j|-local. This is the mechanism behind the whole trainability
question and it is the sprint's most promising theoretical thread:

- the torsion-prior term is exactly **1-local**;
- Legacy's terms are pairwise in CA/CB distance → high-weight in torsion space;
- AMBER's are pairwise over ~230 atoms → higher-weight still, and denser.

So do not assume the Hamiltonian decomposes. Two honest routes: (a) evaluate the energy on
**sampled bitstrings** and optimise CVaR of the sampled distribution — this is what
`run_global_cvar_vqe` already does and it needs no decomposition; (b) build a genuinely
few-body Hamiltonian and **verify** it reproduces the true objective to numerical precision.
If you take route (b), the verification is not optional.

## 5. Objective validity comes before optimisation

For **every** Hamiltonian you propose, measure before optimising it:

- the native's percentile under it,
- Spearman ρ(energy, CA-RMSD) over a real candidate population,
- whether the argmin is nearer the native than a random sample,
- the energy-vs-RMSD scatter.

Sprint 12's hardest-won lesson: **if the objective does not rank the native, optimising it
harder makes the structure worse.** The certified exact optimum of the previous sprint's
assembly Hamiltonian emitted 0.169 Å *worse* than its own anchor. Do not repeat that.

## 6. Reporting

Write `s13/<yourname>_FINDINGS.md` incrementally. Per experiment: hypothesis, motivation,
exact config (seed, split, energy model, representation, ansatz, optimiser, α, shots, budget),
result table, CI, W/L, concentration, leakage audit, interpretation, limitations. Classify
every finding PROVEN / STRONGLY SUPPORTED / SUPPORTED / OPEN / WEAK / REFUTED / INVALIDATED /
UNKNOWN. Never silently replace an earlier result — record the correction.

Your final message to the coordinator: ≤ 60 lines, the numbers that matter, the file paths,
and an explicit statement of what you did NOT establish.
