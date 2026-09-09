<!--
DRAFT — S25 CLEANUP LANE. NOT INSTALLED. The live README.md at the repository root is
untouched. This becomes README.md after the architecture freeze, once the placeholders
below are filled from the frozen pipeline and the Phase-II results.

Every section marked <!-- PLACEHOLDER: FREEZE --> is waiting on a specific fact that does
not exist yet. Each one names the fact it needs.
-->

# Peptide structure prediction, and the instrument that measures it

A structure predictor for short peptides (8–20 residues), and — more to the point — the
measuring instrument built around it. The predictor retrieves candidate backbones from a
library, filters them with a learned distance posterior, selects among them with a
CVaR-VQE, synthesises one structure from the survivors, and relaxes it under AMBER
ff14SB/GBn2. The instrument prices every one of those stages against held-out targets, and
**its findings are the actual output of the project.**

> <!-- PLACEHOLDER: FREEZE -->
> **The headline.** One paragraph: the final mean full-chain Cα-RMSD on the 126
> cluster-disjoint dev targets against the incumbent 3.0483 Å, the basis (point-cloud vs
> built-chain — they differ by 0.156 Å of pure operator choice and are never compared),
> the paired CI both iid and fold-clustered, and the benchmark-60 number if and only if
> the benchmark is spent. If < 3.0 Å was not reached, this paragraph says so plainly and
> states what remains limiting. *Needs: the frozen architecture's Phase-II run.*

---

## Quick start

```bash
git clone <repo> && cd Protein-Folding-Algorithm
pip install -e ".[test]"
pytest                                              # ~5 min on a quiet box

python -m core.pipeline run   --manifest smoke8     # 8 targets, a few minutes
python -m core.pipeline stats --manifest smoke8
```

Manifests, smallest first: `smoke8`, `smoke24`, `dev24`, `tuning126`, `benchmark60`.

**`benchmark60` is sealed.** It was spent once, on pre-registered constants. The CLI
refuses to run it without `--i-am-spending-the-benchmark`. Everything else is developed on
`tuning126`; `dev24` is the cluster-disjoint dev split.

### Running the suite on a busy machine

This box is 15.6 GB and `core.amber` refuses to open OpenMM contexts above 92% physical
memory. Running the whole suite in one process imports torch, pennylane, esm and openmm
together and can cross that line on its own, which surfaces as errors in the AMBER tests
rather than as an out-of-memory. If you see that, run the AMBER files separately:

```bash
pytest tests/ --ignore=tests/test_amber.py --ignore=tests/test_amber_frame_invariance.py
pytest tests/test_amber.py tests/test_amber_frame_invariance.py
```

---

## Where everything is

| if you are looking for | it is in |
|---|---|
| **production code** | `core/` — 11 modules, one implementation per operation |
| **the pipeline** | `core/pipeline.py` |
| **the quantum component** | `core/quantum.py` (reference arm: `foldvqe.py`, `qansatz.py`, `vqe.py`, `objective.py`, `hamiltonian.py`) |
| **physics** | `core/energy.py` (H_Legacy), `core/amber.py` (H_AMBER) |
| **the distance posterior** | `core/predict.py` |
| **the evaluation instrument** | `s12/instrument.py` — every RMSD in the project comes out of here |
| **tests** | `tests/` — 12 files, equivalence + scientific invariants |
| **reproducibility evidence** | `verify/*.json`, `bench_results/*.json` |
| **the research record** | `FINDINGS.md`, and each sprint's `LEDGER.md` |
| **results, structures and renders** | `results/` |
| **sprint experiments** | `s12/` … `s25/` — the record of how every conclusion was reached |

**Read [`FINDINGS.md`](FINDINGS.md) if you want to know what this project learned.** It is
the consolidated research record with a corrections ledger at the top, because several
findings withdraw or revise earlier ones. Most of the substantive results in this project
are negative results, and they are the part worth reading.

---

## The pipeline

> <!-- PLACEHOLDER: FREEZE -->
> The stage table below is the pre-freeze architecture. Component ORDER may change if an
> experiment justifies it, and the four pillars fix what must be present, not the order.
> *Needs: the frozen component order and the final K/M/λ/α constants.*

| stage | what it does | where |
|---|---|---|
| **retrieve** | K=500 BLOSUM62-scored windows over fold-disjoint peptides plus fold fragments | `core.data`, `s5.lib` |
| **filter** | a learned distance posterior scores the pool by Bayes risk; top-75 survive | `core.predict` |
| **evaluate** | H_Legacy (11 terms at `DEFAULT_WEIGHTS`) and H_AMBER (ff14SB/GBn2), rank-standardised | `core.energy`, `core.amber` |
| **select** | CVaR-VQE over the Hamiltonian-derived objective; the CVaR tail is the ensemble | `core.quantum` |
| **synthesise** | coordinate-average the survivors, then multi-start project onto the ideal-geometry manifold | `core.geometry`, `core.project` |
| **relax** | AMBER ff14SB/GBn2 restrained relaxation, converged | `core.amber` |

`ARCHITECTURE.md` has the equations, the Hamiltonians and the data flow.

### Resuming

Every target writes one atomic JSON checkpoint under
`bench_results/cache/<config-key>/<pdb>.json`. The config key is a SHA-1 over every
parameter that can change a number — K, M, penalty, λ, the AMBER schedule, the fold count,
the pair separation, the optimiser iteration cap, the tie-break rule, **and the live
backend set**. Resuming is just re-running the same command; a run with a changed
parameter gets a different key and cannot silently inherit the old answer.

---

## How to reproduce a result

1. **A single number from `FINDINGS.md` or a sprint `LEDGER.md`.** Each entry names the
   module that produced it and the result JSON it wrote. Re-run that module.
2. **The whole benchmark.**
   ```bash
   python -m core.bench --arm optimised --manifest tuning126 --workers 6 --fresh
   ```
   `--fresh` deletes that arm's checkpoints first, which is what makes it a cold
   measurement.
3. **The reference arm.** `core/__init__.py` does not import root modules directly; it
   asks `core.backend("geometry")`, which returns the consolidated module if it exposes
   the whole contract and otherwise the root module it replaces.
   ```bash
   CORE_BACKENDS=legacy python -m core.bench --arm baseline --manifest tuning126 --fresh
   ```
   That is the pre-consolidation implementation, and it is the arm every equivalence claim
   is proven against. `core.backend_report()` records which implementation each backend
   resolved to; that set goes into the results file **and into the cache key**, so an
   optimised result can never be served out of a baseline cache.
4. **The equivalence audits.** `verify/run_equiv.sh`, `verify/run_equiv2.sh`.

---

## Where structures and renders live, and how to view results

> <!-- PLACEHOLDER: FREEZE -->
> *Needs: the results-lab lane's final layout.* Expected shape:
>
> ```
> results/
> ├── summary/     results.json, results.csv, leaderboard.csv,
> │                final_report.md, professor_brief.md
> ├── structures/  T001__method.pdb — deterministic naming, one per target per method
> └── site/        the 3D renderer: native overlay, RMSD shown prominently
> ```
>
> This section gets: how to open the renderer, what the overlay shows, which
> configuration each leaderboard row is, and the exact command that regenerates all of it.

The seven-configuration comparison — Legacy, AMBER, Distogram and the four combinations —
is run on the same targets, the same folds, the same metric, the same basis, with the
CVaR-VQE as the matched selector throughout.

---

## The instrument discipline

These rules are why the numbers here mean anything. Breaking one does not produce a wrong
number; it produces a number that cannot be interpreted.

- **`benchmark60` is held out.** All development is on `tuning126`; `dev24` is the
  cluster-disjoint dev split.
- **Basis discipline.** Point-cloud RMSD and built-chain RMSD are never compared. The gap
  is 0.156 Å of pure operator choice. State the basis on both sides of every comparison.
- **ORACLE / ACHIEVABLE / PRODUCTION** are labelled at every appearance. An oracle is
  never promoted to a result.
- **Nested CV for anything fitted.** An in-sample fit is reported as the optimism, never
  as the result.
- **MDE = 2.8016 × SE, per comparison**, and SE is reported beside every mean. An effect
  at 0.7–1.3× MDE is a Type-M zone and is not a result.
- **Fold-clustered CIs beside iid ones.** iid CIs on this instrument are anticonservative:
  each of the 5 distogram fold models saw ~100 of the other 125 dev natives.
- **Folds and clusters are pinned** in `peptide_folds.json` and `peptide_clusters.json`,
  and the benchmark in `results/benchmark_manifest.json`. Correcting the identity
  clustering once silently moved 13 benchmark targets and invalidated every fold model.
- **Natives are reporting labels only.** The deployable path sees a `deployable_view`;
  `label()` is the only function that opens a native and it runs after the structure is
  final. `tests/test_pipeline.py` NaN-poisons every native quantity and asserts every
  emitted coordinate is bit-identical.
- **A known filter defect stands** (S9-10, S10-4): `identity` normalises by the *longer*
  sequence, so a long fragment can carry a target's exact k-mer and still align below 0.6.
  Priced at +0.0030 Å on the benchmark, +0.0004 Å on tuning. It is documented rather than
  fixed because fixing it would re-pin the folds.

---

## Setup

```bash
pip install -e .            # runtime deps
pip install -e ".[test]"    # plus pytest
```

Nine runtime dependencies, each annotated in `pyproject.toml` with the module that needs
it: numpy, scipy, torch, openmm, pennylane, pennylane-lightning, fair-esm, scikit-learn,
biopython. `lightning.qubit` is the device `core.quantum` builds circuits for and the one
`tests/test_quantum.py` asserts probabilities against with `==`, so it is not an
interchangeable backend. Likewise there is no fallback for OpenMM: AMBER is the validity
stage, not an option.

### Data this repository does not carry

`pdbs/` (61 structures) is tracked — it is the reporting set. The rest is large, ignored,
and rebuilt on demand by the module that owns it: `prots/` (5.5 GB) and `pdbs_ext/`
(555 MB) from the RCSB queries, `esm_cache.npz` (1.5 GB) from `esm_features.compute`,
`peptide_db.npz` and `fragment_db*.npz` from their `build` entry points, the `*_models/`
checkpoint directories from `train_fold`, and the per-sprint `cache/` directories from the
sprint module beside them. `.gitignore` names the owner of each and carries the exact
`git show` commands for the four RCSB fetch scripts, which live in history rather than in
the tree.

> **Two traps worth knowing.** `peptide_folds.json`, `peptide_clusters.json` and
> `catrace_prior.npz` are **write-on-first-use caches**: if the file is missing it is not
> an error, it is silently re-derived. Moving them therefore does not crash, it re-pins
> the folds. And `pdbs/` is globbed and sorted by three modules — the sort is stable given
> the same set of files, but BLOSUM similarity sums are small integers with large tie sets
> resolved by enumeration order, so adding or removing one `.pdb` moves up to 47 of the
> 500 pool members. **Do not add or remove files there.**

---

## Layout

```
README.md            this file
ARCHITECTURE.md      data flow, equations, Hamiltonians, VQE/CVaR, known limitations
FINDINGS.md          the research record — corrections ledger first
pyproject.toml       the only config; every dependency annotated with its consumer
core/                the production path (11 modules)
tests/               the equivalence and scientific-invariant suite (12 files)
verify/              standalone audits — and a runtime input/output directory of core/project.py
s12/instrument.py    the evaluation instrument
<24 root modules>    the CORE_BACKENDS=legacy reference arm
s5/ s7/ s8/ s9/      the reference arm's sprint packages
s12/ … s25/          the research record: experiments, LEDGERs, PREREGs, FINDINGS, results
pdbs/                the 61-structure reporting set — order-sensitive, do not modify
results/             manifests, structures, renders, the summary tables
bench_results/       golden results, whitelisted past a blanket ignore
```

### Why the tree is flat

Four things pin it, and each was verified rather than assumed.

1. **The root modules resolve data relative to themselves.** `peptide_db.py` sets
   `BASE = dirname(abspath(__file__))` and assumes it *is* the repository root, then
   resolves `pdbs/`, `results/benchmark_manifest.json`, `peptide_folds.json` and
   `peptide_clusters.json` from it.
2. **The sprint packages hard-code their hop count.** Every module under `s5/ s7/ s8/ s9/`
   and `s12/` does `sys.path.insert(0, dirname(dirname(abspath(__file__))))`. One level
   deeper and that points at the wrong directory.
3. **`core/bench.py` builds a subprocess preamble** that inserts the repository root and
   imports `s9.final`. A shim in a parent directory would not reach the child process.
4. **`verify/` is not just scripts.** `core/project.py` writes its evidence there and
   reads `project_inputs.json` and `project_equiv.json` back at runtime;
   `tests/test_project.py` globs `verify/project_stability*.json`;
   `tests/test_integration.py` imports three `verify/` scripts as modules.

---

## History

The pre-consolidation tree — 372 Python files, ten sprint packages, every one-off
experiment — is at commit `5fa05cd`. `git checkout 5fa05cd` recovers it whole, and
`git log --diff-filter=D --name-only` finds any single deleted file by name. The Sprint 11
consolidation reduced the production import closure to the 11 modules in `core/` **without
changing a single scientific decision**; the root modules survived because they are the
arm that claim is proven against, and `tests/test_equivalence.py` is where it is proven.

`_archive/` holds console logs moved out of the Sprint 5–9 packages in the 2026-09-04
consolidation. Nothing was deleted — see `_archive/README.txt`, which also records the one
log that is cited by path in `FINDINGS.md`.

`distogram_models_large.STALE-PRE-FOLD-REPIN-DO-NOT-USE/` is quarantined by rename, on
purpose, so that `FRAG_LARGE=1` fails loudly instead of silently loading a model that
predates the fold repin.
