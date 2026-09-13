# State of the repository — brief, 2026-09-12

A one-sitting summary of what this repository is, what it has measured, what it has concluded,
and where it stands operationally. Written from the tree at commit `ae86a124` (2026-09-09,
"refactor: deep clean and reorganize codebase"), working tree clean. Every number below is
quoted from a persisted artefact or a sprint ledger; where a number has been retracted the
retraction is stated, not the original.

Read in this order if you have time for more: `README.md` → `ARCHITECTURE.md` →
`results/summary/professor_brief.md` → `s25/LEDGER.md` → `docs/FINDINGS.md`.

---

## 1. What the project is

A structure predictor for short peptides (9–16 residues), and — more importantly — the measuring
instrument built around it. The predictor retrieves candidate backbones from a library, filters
them with a learned distogram, has a genuine CVaR-VQE select a subset, coordinate-averages the
survivors, projects the average onto an ideal-geometry backbone manifold, and optionally relaxes
under AMBER ff14SB/GBn2. Two genuine physics Hamiltonians (an 11-term "Legacy" potential and
all-atom AMBER via OpenMM) are independently evaluable at the scoring stage.

**Endpoint:** mean full-chain Cα-RMSD to deposited model 1 after optimal superposition.
**Instrument:** 126 cluster-disjoint development targets (`tuning126`), 5 pinned folds, target as
the unit of analysis, paired comparisons with bootstrap and fold-clustered CIs.
**Sealed benchmark:** 60 targets, spent exactly once (S9-10) on pre-registered constants. The CLI
refuses to run it without `--i-am-spending-the-benchmark`.

The programme ran 25 sprints. The last ~14 (S12–S25) were multi-agent research campaigns with
pre-registered falsifiers, per-sprint ledgers, and a standing rule that no superseded claim is
ever deleted. **The findings, not the folding, are the project's actual output.**

---

## 2. The headline result

| | value | basis |
|---|---|---|
| **Production result** | **3.2148 Å** (leaderboard rebuild: 3.2126) | ideal-geometry built chain — a real, valid, renderable structure |
| intermediate | 3.0483 Å | raw coordinate average — 22.3% contracted (mean virtual Cα–Cα 2.961 Å vs 3.812 native, worst bond 0.649 Å), **not a protein structure**, labelled `"raw average (illegal)"` in `core/bench.py` |
| deployed emission | 3.236 Å | after restrained AMBER relaxation |

The basis decision (S25 L4/L8/L11) is that **the built chain is the result** and the point cloud
is a labelled intermediate. It costs 0.166 Å against every earlier report; it was adopted because
it is the honest direction.

**On the sealed benchmark (S9-10, the only read):** synthesis architecture 2.9610 vs the shipped
distogram-argmin baseline 2.9507 — paired +0.0103, 95% CI [−0.160, +0.180], 31W/29L. A tuning
gain of −0.250 did not transfer. **There is no validated accuracy improvement over the baseline.**
Neither the <3.0 Å (dev) nor the <2.5 Å / <2.0 Å targets set in S15–S25 were reached.

For scale, on the same 126 targets (ORACLE arms read the native and are diagnostics only):

| arm | mean Cα-RMSD |
|---|---|
| annealing on Legacy energy | 4.624 |
| random 75-subset of the K=500 pool | 3.425 |
| constant α-helix, zero information | 4.065 |
| sequence-only torsion predictor | 3.770 |
| **shipped pipeline (built chain)** | **3.215** |
| ORACLE perfect distance prior, same pipeline (S24 γ=1) | 2.226 |
| ORACLE best member of K=500 pool | 1.711 |
| ORACLE best member of the full retrieval universe | 1.313 |
| ORACLE torsion-space ceiling, ~24 qubits | 1.594 |
| ORACLE torsion-space distance geometry from true distances | 0.611 |

---

## 3. Repository state

### Layout (flat by design — see README "Layout, and why it is flat")

| what | where |
|---|---|
| production path | `core/` — 11 modules (`data, geometry, predict, energy, amber, quantum, cache, pipeline, bench, project, __init__`) |
| reference arm | 24 root modules (`protein_geometry.py`, `peptide_db.py`, `distogram.py`, `foldvqe.py`, `energy_terms.py`, `amber_refine.py`, …) plus `s5/ s7/ s8/ s9/` — **the baseline arm every equivalence claim is measured against; not leftovers** |
| research log | `s12/` … `s25/` — one directory per sprint: BRIEF, PREREG_*, LEDGER, per-agent FINDINGS, scripts, results JSON |
| research record | `docs/FINDINGS.md` (9,081 lines, S5–S11, corrections ledger first), `docs/CONDENSED_REPORT.md` (S5–S13), `ARCHITECTURE.md` (frozen spec, S25) |
| results lab | `s25/resultslab/` builds `results/structures/` (1,135 PDBs: 126 targets × 9 configurations incl. native, + 1,008 chain variants), `results/summary/` (leaderboard, per-target results, `professor_brief.md`), `results/site/` (static 3D viewer, regenerated) |
| tests | `tests/` — 12 files, 368 tests |
| verification | `verify/` — standalone audits and their JSON; also an output directory of `core/project.py` |
| data | `pdbs/` (61 structures, tracked, **do not add/remove files — sort order feeds BLOSUM tie-breaking**), pinned caches at root |

Git: 4,560 tracked files, 10 commits on `main` (the pre-consolidation tree with all 372 Python
files is at `5fa05cd`). The last commit removed 5,395 tracked non-source files (per-target caches,
a 1,139-file synthetic quarantine, RCSB downloads) and verified every `.py` AST-identical to its
parent once docstrings are stripped.

### Pinned things you must not casually touch

- `peptide_folds.json`, `peptide_clusters.json`, `catrace_prior.npz` are **write-on-first-use
  caches**. Deleting or moving them silently re-derives folds; doing that once moved 13 benchmark
  targets and invalidated every trained model.
- `results/benchmark_manifest.json` — the sealed 60.
- `distogram_models/` (5 leave-fold-out checkpoints). `distogram_models_large.STALE-…-DO-NOT-USE/`
  is quarantined by rename so `FRAG_LARGE=1` fails loudly.
- `pdbs/` file set (see above).

### Not in git (rebuilt on demand, or at risk)

- `prots/` 5.5 GB, `pdbs_ext/` 555 MB, `esm_cache.npz` 1.5 GB, `peptide_db.npz`, `fragment_db*.npz`,
  `*_models/` — owners named in `.gitignore`.
- `_archive/logs/` — 99 console logs from S5–S9; one (`s8/predictor_report.log`) is the sole
  source of a published number.
- **`s24/cache_amber/` (126 `.npz`, 1.5 MB, 63,000 genuine ff14SB/GBn2 single points) is caught by
  the blanket `*.npz` ignore rule and is NOT tracked.** S25 L13 flagged this as the most urgent
  operational item; the S23–S25 text and JSON were committed since (324 files tracked) but this
  cache was not. Cheap to whitelist.

### Declared defects that stand

- **Identity leak** (`core/data.py:188`): Needleman–Wunsch identity normalised by the *longer*
  sequence, so a verbatim k-mer inside a longer fragment can score <0.6 and land in another fold.
  4/126 dev targets and **2/60 sealed-benchmark targets** carry a verbatim self-copy in their own
  fold model's training set. Priced at +0.0004 Å (dev) / +0.0030 Å (benchmark, S10-4); the 2/60
  benchmark effect was deliberately never quantified because that would open the benchmark.
- Ideal-geometry projection uses a constant 3.804 Å virtual bond → cannot represent cis-peptides.
- The production leaderboard row shows `pool_gate = WARN` with 2 per-target pool violations
  (targets where the emitted structure beats the pool's own ORACLE best member — possible because
  averaging and projection can land off the pool). Declared, not a defect.
- Fold models each saw ~100 of the other 125 dev natives → per-target dev results are not
  independent; **iid CIs are anticonservative; fold-clustered CIs are the honest ones.**

---

## 4. The production pipeline (frozen, `ARCHITECTURE.md`, commit `a15406c`)

```
sequence → ESM-2 embedding → distogram (17-bin distance posterior, leave-fold-out)
        → BLOSUM62 retrieval, K=500 real protein windows (787 peptides + 6,003 fragments)
        → L1 Bayes-risk score s(w) = Σ_p w_p · E_post|t − d_p(w)|,  w_p = 1/(sd_p+0.5)
        → top-128 → E = zrank(scores) → H = diag(E)
        → CVaR-VQE: 7 qubits, 3 layers, 21 params, exact statevector, Adam on exact
          parameter-shift gradient, F(p)=E_p[E] − T·H(p), VQE_LFO table (T=0.3; α=1.0 on 3/5 folds)
        → CVaR tail → uniform coordinate average (top-75 rung)   [3.048 Å, point cloud]
        → multi-start ideal-geometry projection, ramah@0.3        [3.215 Å, built chain]
        → optional restrained AMBER ff14SB/GBn2 relaxation, k=10  [3.236 Å]
```

`core.backend()` resolves each operation to `core.*` or to the root reference module;
`CORE_BACKENDS=legacy` forces the reference arm. Results are checkpointed per target under a
SHA-1 config key that includes the live backend set, so an optimised result can never be served
from a baseline cache. `bench_results/compare_tuning126.json` shows `science_delta` = 0 between
the 1-worker baseline and 8-worker optimised runs — bit-identical science across the two arms is
the intended guarantee.

---

## 5. Results so far — what has been established

### 5.1 The binding constraint is the distance prior, and it is not reachable by re-reading it

- **The prior's derivative is the only steep lever** (S24 L13): interpolating the posterior toward
  truth with everything else fixed moves the endpoint −2.15 Å per unit γ; γ=0.0225 would reach
  3.0 Å. Four independent interpolation constructions agree. But that slope holds only along the
  native's own direction; a real operator travelling 25% at cos 0.5 is worth +0.024 Å.
- **Every re-reading of the existing posterior is closed** (S25 L2/L6/L7/L12): width calibration
  (the posterior is ~2× over-confident, and calibrating it makes RMSD *worse*), location shifts,
  mode-seeking, alternative risk functionals, metric projection, de-quantisation (the score's
  per-pair target takes only 17 values, gaps up to 4 Å). For 3 of 4 families the best parameter
  chosen *with full leakage* is the shipped default. Across 51 arms,
  corr(progress-toward-truth, endpoint) = +0.054.
- **The pool's error is 68% common-mode** (S23 L9, exact identity): a bias shared by every pool
  member is invisible to any within-pool statistic, so every operator downstream of retrieval is
  capped.
- **The distogram's errors are worse than random errors of the same magnitude** (S18, five
  controls + independent reimplementation; S19 mechanism): every predictor emits a "typical
  peptide of that length"; the coherent error is the gap between typical and this native, and
  score-ordered gates *amplify* it.

### 5.2 Selection / ranking inside the pool is closed

- Pool contains a sub-2 Å candidate for 81% of targets (full universe) / 46% (top-75), yet
  widening K makes the realised answer *worse* — extra windows displace near-native members from
  the shortlist (S17). A perfect ranker inside the shipped top-25 returns 2.609 Å.
- The terminal operator consumes the set **mean**, not the set best
  (d_out ≈ 1.16·mean + 0.04·best, R² 0.89). The last untested ranker class — a set-transformer
  over the full signed deviation map — has a flat learning curve (S13).
- Exhaustive enumeration of the whole 2ⁿ latent (S21): exact argmin ties a zero-evaluation pool;
  ORACLE over the same set is 1.72 Å better with 0/122 reversals. Search is saturated;
  discrimination binds.
- Consensus/medoid family (9 arms): 0% of the in-pool gap. Diversity-maximising selection: dead.
- Per-target optimal set size m is real and transfers across pool halves (−0.239 Å, S22 L4), but
  **five independent router constructions fail to predict it native-free**, two significantly
  harmful held-out (S22 L7, S23 L7).

### 5.3 The physics energies do not rank nativeness on this pool

- **Seven-configuration suite, genuine CVaR-VQE for all seven, matched readout** (S25 L16):
  Distogram 3.058 · AMBER+Dist 3.132 · Legacy+Dist 3.215 · L+A+Dist 3.253 · random-75 3.425 ·
  Legacy+AMBER 3.674 · Legacy 3.755 · AMBER 3.881. **Legacy (+0.330) and AMBER (+0.455) are
  measurably worse than a random subset, 5/5 folds.** Rank-permuting a physics channel while
  preserving its marginal *improves* the endpoint.
- Legacy is a compactness model (top-75 Rg −0.758 Å); AMBER has the opposite sign (+1.103 Å);
  ρ(Legacy, AMBER) = −0.09; ρ(AMBER, distogram) = −0.019, CI includes zero (orthogonal).
- The functional lever is closed in all five forms: filter, partition, additive score, per-target
  audit sign, and an oracle weight chosen with full leakage worth 0.015 Å (S24 L16, S25 L18).
- AMBER relaxation's apparent −0.023 Å gain was against a projection that is itself 0.155 Å worse
  than doing nothing; a matched-magnitude random displacement is as accurate (S16). In torsion
  space, Legacy's certified global optimum is +0.139 Å *worse* than random sampling (S13).
- Raw moment standardisation of AMBER is not monotone in float64 (1e28 outliers), breaking
  argsort on 40/126 targets — the rank normalisation is the correct one.

### 5.4 What the quantum component is and does

- Verified genuine (S25): exact statevector to 5.6e-17 against an independent simulator; exact
  parameter-shift gradient, cos 1.000000000 vs finite differences; real CVaR. 39 equivalence
  tests + 90 definition tests pass.
- **Its selection is provably classical:** the CVaR tail's support is always a subset of an
  initial prefix of the energy order (2,592 adversarial cells, 0 violations). It can delete a
  top-m member, never add one.
- **The readout is insensitive to it:** the trained state sits 0.902 nats / 45% of its mass from
  its own analytic Gibbs optimum and the endpoint difference is 0.24× MDE. RMSD tracks readout
  entropy (ρ −0.74), not α or T; the circuit arms land on the curve fitted to the nine
  no-circuit arms at +0.009 Å.
- **The Hamiltonian barely changes between targets** (S25 L17): `E = zrank(sorted scores)` is the
  standardised rank ladder up to tie-averaging (worst deviation 1.18% of range) → effectively two
  trained states in the whole deployment. This is the mechanism for why deeper ansätze and larger
  χ ordered nothing, five times over five sprints.
- The one genuine positive: the optimiser trains — beats best-of-200 from the untrained circuit
  at every temperature, closing 78–89% of the free-energy gap. No barren plateau at any width
  measured; the CVaR non-linearity flattens the gradient-variance decay.
- **Withdrawn:** "+0.113 Å CVaR contribution" (measured at a non-deployed temperature; 0.51× MDE
  as a paired contrast).

### 5.5 Trainability result (S13, the publishable half)

- Exact locality theorem in torsion space: d_ij depends on exactly the j−i−1 residues between i
  and j (agreement 1.0000); no 2-local Ising form of a distance objective exists in this encoding.
- The Pauli spectrum of a force field predicts its gradient variance with no free parameter
  (Legacy mean weight 2.236, measured/predicted 1.006; AMBER 3.015, 1.001) — after monotone
  rank-preserving conditioning; raw AMBER's spectrum is a δ-spike artefact of its worst clash.
- QNG has nothing to fix (metric full rank, g_ii = 0.25 exactly); the metric contains no
  Hamiltonian (bit-identical across energy models).

### 5.6 Other closed routes

Chemical-shift torsion restraints (54/126 coverage; ORACLE-perfect torsions still 2.021 Å);
sequence-only torsion prediction (φ carries no sequence signal: 36.1° vs 36.4° blind); degree-1
truncation of the distance objective (S18, five falsifiers fired); native-free error-direction
steering (S16, every arm chose "do nothing"); learned residual / from-scratch candidate
generation (S24, closed by oracle upper bound and corpus census: the corpus is small and has
almost no β-sheet); Legacy→AMBER continuation and preconditioning (S21); probability-weighted
readout (S23); per-target scale correction (S23: real, universal, unreachable native-free);
growing the fragment library (~0.43 Å per decade); a fresh benchmark (all 204 identity clusters
of 9–16-mers are spent; the world supply of containment-fresh targets is 16, 10 of them amyloid).

### 5.7 What is open

1. Whether a genuinely better distance predictor is *obtainable* — the project measured what one
   would buy, not whether one exists. ESM-pca32 is worth −0.288 Å over one-hot on selection.
2. The 2/60 benchmark self-copy leak (declared, unquantified by design).
3. Publishing the trainability half (S13 + S25 QUANTUM.md).
4. The ansatz's dynamical Lie algebra (not measured; scoped as "nowhere near a 2-design").

---

## 6. Tests

### 6.1 Unit / invariant suite (`pytest tests/`)

368 tests across 12 files: `test_amber` 15 · `test_amber_frame_invariance` 3 · `test_cvar` 90
(asserts the CVaR *definition*, both `CVAR_SORT_CUTOFF` branches) · `test_data` 41 ·
`test_energy` 9 · `test_equivalence` 14 (core vs reference arm) · `test_geometry` 24 ·
`test_instrument` 43 (added S25 — `s12/instrument.py` had no test despite producing every RMSD) ·
`test_integration` 25 · `test_pipeline` 37 (NaN-poisons every native quantity and asserts emitted
coordinates are bit-identical) · `test_project` 28 · `test_quantum` 39 (equivalence against
`lightning.qubit`, `==` on probabilities).

Last recorded status: 368 pass, 0 warnings (commit `ae86a124`). S25 L14 notes the AMBER tests
carry a memory-ceiling guard: on this 15.6 GB box a full run can push past 92% RAM, and the guard
turns those into skips/errors that are not code defects.

**S26 governed run (2026-09-13, branch `s26`, commit `601a39c7`, under `s26/governor.py`; per-file counts, every skip reason and peak RSS in `s26/TEST_RUN.md` / `s26/results/test_run.json`): 370 tests, 357 passed, 13 skipped, 0 failed, 0 errors, 0 memory-guard skips**, split into three jobs so no more than one AMBER job from the lane was live: the 10 non-AMBER files as job `pytest_core_post` (351 tests, 338 passed, 13 skipped; 240 s, peak RSS 1.692 GB), `tests/test_amber.py` as AMBER job `pytest_amber` (16 passed; 261 s, 0.872 GB) and `tests/test_amber_frame_invariance.py` as AMBER job `pytest_amber_frame` (3 passed; 256 s, 0.324 GB). Passed per file: `cvar` 90 · `data` 42 · `energy` 9 · `equivalence` 11 (+3 skipped) · `geometry` 24 · `instrument` 43 · `integration` 17 (+8 skipped) · `pipeline` 35 (+2 skipped) · `project` 28 · `quantum` 39 · `amber` 16 · `amber_frame_invariance` 3. The 13 skips are real, none from the memory guard: 3 in `test_equivalence.py` and 8 in `test_integration.py` are `set VERIFY_SLOW=1` opt-ins, 2 in `test_pipeline.py` are absent artefacts (`bench_results/optimised_tuning126_w6.json`; no smoke result on disk). The suite is 370 rather than 368 because S26 added one test to `test_amber.py` (the memory-guard message) and one to `test_data.py` (the identity flag ships dark). The same files ran before any S26 edit (job `pytest_core`, tree `a4db170c`): 337 passed, 13 skipped, 0 failed, identical skip list.

### 6.2 Standalone audits (`verify/`)

Ansatz, CVaR, AMBER platform/threads/affinity, determinism, gradient key collision, headline,
leak, legacy, projection equivalence/exactness/stability, VQE leave-fold-out, containment — each
with its JSON. `verify/run_equiv2.sh` is flagged stale.

### 6.3 Results-lab gates (`s25/resultslab/`)

Four release gates on every exported structure set: pool-oracle gate, difficulty-correlation
gate, per-target sd floor, mandatory provenance REMARK in every PDB header (build refuses files
without one). Every exported structure reproduces its own RMSD through the instrument to within
PDB quantisation (worst 2.6e-04 Å). A synthetic native-plus-noise leaderboard found on disk in
S25 was quarantined (its README survives as `docs/quarantine-synthetic-proofbuild.md`).

### 6.4 Experimental discipline — the standing rules

Every experiment since S15 is pre-registered (`PREREG_*.md`) with falsifiers before it runs;
artefacts are saved atomically with module sha256, git commit and dirty flag; `complete: true`
is gated on a full key set. Comparison statistics: paired per target, bootstrap CI, fold-clustered
CI beside iid, W/L, drop-top-10 concentration check, and **MDE = 2.8016 × SE per comparison**
(the "0.084 Å constant" was wrong by up to 84× in both directions). Controls must match the
operator's space; a zero-information control must be plausible (constant α-helix), not uniform;
a best-of-K null is inapplicable when the K arms are correlated (report k_eff); tied argmins are
averaged, never broken by array order (that trap once invented a 1.386 Å winner); grid oracles
are order statistics — only split-half transfer arms count.

**Retraction count:** 23 corrections across S5–S13; five S25 claims retracted by internal audit
(L5, L7, L9, L11, L15). Nothing superseded is deleted; ledgers carry both.

---

## 7. Sprint index (one line each)

| sprint | question | outcome |
|---|---|---|
| S5–S9 | build the retrieval + distogram + synthesis architecture | tuning gain −0.250; benchmark +0.010 (null); geometry better, accuracy not |
| S10 | leakage adjudication, opportunity ladders | non-replication is concentration, not contamination; recognition is the barrier |
| S11 | consolidation into `core/` | two-arm equivalence, bit-identical science |
| S12 | adversarial/aggregation/assembly | terminal operator consumes the set mean; multi-piece assembly closed |
| S13 | torsion space, quantum architecture | locality theorem; Pauli-spectrum → gradient variance; neither energy ranks |
| S14 | chemical-shift restraints, VQE redesign | shift route closed by arithmetic; structural objective's optimum is fine, energies' is not |
| S15 | paper-driven programme | 3.321 Å (worse); error *shape* beats magnitude; 0.611 Å true distance-geometry floor |
| S16 | native-free steering | every arm chose "do nothing"; AMBER relax gain dissolves under matched-random control |
| S17 | averaging → selection readout | shortlist, not ranker, binds; widening K hurts; selection closed at four levels |
| S18 | degree-1 objective | all five falsifiers fired; 1/sd² weighting validated native-free |
| S19 | why errors are maladaptive | pool has a systematic error direction the predictor reproduces; shared-referent floor |
| S20 | landscape geometry, Legacy vs AMBER | Legacy is compactness; AMBER difficulty is a steric singularity; sampler hypothesis refuted |
| S21 | CVaR-VQE as selector | exhaustive latent enumerated; argmin ties the pool; MDE is per-comparison; χ orders nothing |
| S22 | final campaign, routing | per-target m transfers but no router finds it; readout-H lever = A2b null |
| S23 | drive dev < 3.0 | per-target scale real but unreachable; error 68% common-mode |
| S24 | learned candidate generation | closed by oracle bound and corpus census; prior-attribution ladder is the only steep lever |
| S25 | endgame: RMSD, freeze, clean, results lab | basis decision (built chain); physics worse than noise; two trained states; `ARCHITECTURE.md` frozen |

---

## 8. Operational items outstanding

1. **Whitelist `s24/cache_amber/*.npz`** in `.gitignore` and commit it (1.5 MB, 63k AMBER single
   points, a serialised OpenMM run; only copy is this disk).
2. `verify/run_equiv2.sh` is stale (`PROJECT_GRAD` no longer read from env).
3. Six unused-import proposals in `core/` held until the results/physics lanes stop importing
   `core.pipeline`.
4. `_archive/logs/` is untracked by decision; one log is a cited source.
5. Machine budget: 15.6 GB RAM, ~11 GB usable, two heavy jobs max; 8 cores = 4 fast + 4 slow
   ≈ 6.43 core-equivalents.
