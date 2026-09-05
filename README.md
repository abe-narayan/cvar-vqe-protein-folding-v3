# Peptide structure prediction, and the instrument that measures it

A structure predictor for short peptides (roughly 8-20 residues), and — more to the point —
the measuring instrument built around it. The predictor retrieves candidate backbones from a
library, filters them with a learned distogram, synthesises one structure from the survivors,
and relaxes it under AMBER. The instrument prices every one of those stages against a
held-out benchmark, and its findings are the actual output of the project.

**Read [`FINDINGS.md`](FINDINGS.md) first.** It is the consolidated research record: 72
findings across seven sprints, with a corrections ledger at the top because several of them
withdraw or revise earlier ones. The headline is there and it is a negative result:

> On the 60-target held-out benchmark, the synthesis architecture scores **2.9610** against
> the shipped distogram baseline's **2.9507** — a paired gain of **+0.0103, 95% CI
> [-0.1596, +0.1803]** — against a tuning-instrument gain of -0.250. The mean transferred;
> the effect did not. **There is no validated accuracy improvement over the baseline.**

The architecture is physically better-behaved than the baseline (exact 3.804 Å virtual
bonds, native-like torsion statistics, strain removed) and no more accurate. The binding
constraint is *recognition* — ranking within the near-native band — and S9-8 measured that
closure as informational rather than architectural.


## Where everything is

| if you are looking for | it is in |
|---|---|
| **production code** | `core/` — 11 modules, one implementation per operation |
| **the four pipeline stages** | `core/pipeline.py` |
| **quantum code** | `core/quantum.py`; reference arm `foldvqe.py` `qansatz.py` `vqe.py` `objective.py` `hamiltonian.py` `budget.py` `refine2.py` |
| **retrieval** | `core/data.py`, `s5/lib.py`, `s7/debias.py`, `s7/poolsize.py` |
| **physics** | `core/energy.py`, `core/amber.py`; reference arm `energy_terms.py` `legacy_field.py` `legacy_refine.py` `amber_refine.py` `amber_hamiltonian.py` `protein_geometry.py` |
| **tests** | `tests/` — 9 files, equivalence and scientific invariants |
| **experiments** | `s5/ s7/ s8/ s9/` — sprint packages, and still load-bearing (see below) |
| **reproducibility evidence** | `verify/*.json`, `bench_results/*.json`, and the sprint packages' own `*.json` |
| **verification scripts** | `verify/*.py`, `verify/run_equiv*.sh` |
| **data** | `pdbs/` (61 structures), `results/` (two frozen manifests), and the pinned caches at root |
| **the research record** | `FINDINGS.md` |

The 24 modules at the repository root are the **reference arm** — the implementations the
consolidated `core/` path is measured against. Grouped by role:

| role | modules |
|---|---|
| geometry, and the root of the import graph | `protein_geometry.py` (41 importers) |
| data, folds, clusters, manifests | `peptide_db.py`, `fragment_db.py` |
| prediction | `distogram.py`, `pairnet.py`, `priors.py`, `esm_features.py` |
| energy and physics | `energy_terms.py`, `legacy_field.py`, `legacy_refine.py`, `amber_refine.py`, `amber_hamiltonian.py` |
| quantum | `foldvqe.py`, `qansatz.py`, `vqe.py`, `objective.py`, `hamiltonian.py`, `budget.py`, `refine2.py` |
| representation and geometry helpers | `representations.py`, `sidechains.py`, `torsion_lib2.py`, `catrace.py`, `floor.py` |

`esm_features.py`, `torsion_lib2.py`, `legacy_refine.py`, `representations.py`,
`sidechains.py`, `floor.py`, `budget.py`, `catrace.py`, `priors.py`, `pairnet.py` and
`refine2.py` are reached by `core/` on **both** arms — they are production, not archive.


## The pipeline

Four stages, in `core.pipeline`:

| stage | what it does | where |
|---|---|---|
| **retrieve** | K=500 BLOSUM62-scored windows over fold-disjoint peptides plus fold fragments | `core.data`, `s5.lib` |
| **filter** | a learned distogram scores the pool; top-75 survive | `core.predict` |
| **synthesise** | coordinate-average the survivors, then multi-start project onto the ideal-geometry manifold with `ramah@0.3` | `core.geometry`, `core.project` |
| **relax** | AMBER ff14SB/GBn2 restrained relaxation, k=10, converged | `core.amber` |

Two optional components sit outside that path and are measured separately, because they add
work the reference path does not do: a **CVaR-VQE hypothesis selector** (`core.quantum`) and
a **Legacy late refiner** (`core.energy`). Enable both with `--components`.

`core/` is the consolidated production closure — one implementation per operation. Each
module is cheap to import: nothing loads a database, a model or the ESM bank at import time.

    core.data       the peptide database, folds, identity clusters, manifests
    core.geometry   backbone build / Kabsch / the manifold projection
    core.predict    the distogram predictor
    core.energy     the classical energy terms
    core.amber      the ff14SB/GBn2 restrained relaxation
    core.quantum    the VQE / CVaR optimiser
    core.cache      the on-disk caches, and subset extraction out of the big banks
    core.pipeline   the four stages, parallel and resumable
    core.bench      the harness that measures them, per stage

### The legacy backends, and why the root modules are still here

`core/__init__.py` does not import root modules directly. It asks `core.backend("geometry")`,
which returns `core.geometry` if that module exists *and* exposes the whole contract, and
otherwise the root module it replaces (`protein_geometry`). That fallback is not vestigial:
`CORE_BACKENDS=legacy` is how `core.bench --arm baseline` guarantees it is measuring the
reference path rather than the consolidated one, so the root modules and the surviving
`s5/`, `s7/`, `s8/`, `s9/` packages **are** the baseline arm. `core.backend_report()` records
which implementation every backend resolved to, and that set goes into the results file and
into the cache key — so an optimised result can never be served out of a baseline cache.

    python -c "import core; print(core.backend_report())"
    CORE_BACKENDS=legacy python -c "import core; print(core.backend_report())"

The second prints `amber_refine, peptide_db, energy_terms, protein_geometry, legacy_field,
s7.audit, distogram, foldvqe, s8.project`. **Those nine modules and their transitive closure
are why the root modules and the sprint packages survived consolidation.** They are not
leftovers; they are the arm every equivalence claim in this project is proven against, and
deleting them would not tidy the repository, it would make those claims unreproducible.

The full pre-consolidation tree — all 372 Python files, ten sprint packages, every one-off
experiment — is at commit **`5fa05cd`**. `git checkout 5fa05cd` recovers it whole;
`git log --diff-filter=D --name-only` finds any single file by name.


## Running a full-quality experiment

    python -m core.pipeline run --manifest tuning126 --workers 6
    python -m core.pipeline stats --manifest tuning126

Manifests: `smoke8`, `smoke24`, `tuning126`, `dev24`, `benchmark60`.

To measure rather than just run — per-stage clocks, both arms, a real speedup:

    python -m core.bench --arm baseline  --manifest tuning126 --fresh
    python -m core.bench --arm optimised --manifest tuning126 --workers 6 --fresh
    python -m core.bench --compare bench_results/baseline_tuning126.json \
                                   bench_results/optimised_tuning126.json

`--fresh` deletes that arm's checkpoints first, which is what makes it a cold measurement;
a run that inherited work from disk is not one, and the results file says so. `--components`
runs the full four-component system, which is *more* work than the reference path, so its
wall clock answers a different question and is labelled as such.

### Resuming

Every target writes one atomic JSON checkpoint under
`bench_results/cache/<config-key>/<pdb>.json`. The config key is a SHA-1 over every parameter
that can change a number — K, M, penalty, λ, the AMBER schedule, the fold count, the pair
separation, the optimiser iteration cap, the tie-break rule, and the live backend set.

**So resuming is just re-running the same command.** A run that dies at 99% resumes at 99%.
A run with a changed parameter gets a different key and cannot silently inherit the old
answer. `stats` reports how many targets were resumed, and the results file records it,
because a resumed run is not a cold timing measurement.


## The instrument discipline

These rules are why the numbers in `FINDINGS.md` mean anything. Breaking one does not produce
a wrong number; it produces a number that cannot be interpreted.

- **`benchmark60` is held out.** It was spent exactly once, by S9-10, on pre-registered
  constants. The CLI will not run it without `--i-am-spending-the-benchmark`. All
  optimisation runs on `tuning126`; `dev24` is the cluster-disjoint dev split.
- **Folds and clusters are pinned** in `peptide_folds.json` and `peptide_clusters.json`, and
  the benchmark in `results/benchmark_manifest.json`. Correcting the identity clustering once
  silently moved 13 benchmark targets and invalidated every trained fold model.
- **Natives are reporting labels only.** The deployable path sees a `deployable_view`;
  `label()` is the only function that opens a native and it runs after the structure is
  final. `tests/test_pipeline.py` NaN-poisons every native quantity and asserts every emitted
  coordinate is bit-identical.
- **A known filter defect stands** (S9-10, S10-4): `identity` normalises by the *longer*
  sequence, so a long fragment can carry a target's exact k-mer and still align below 0.6.
  Priced at +0.0030 Å on the benchmark and +0.0004 Å on tuning. It applies to the
  member-level filter everywhere in this repository.


## Setup

    pip install -e .          # runtime deps
    pip install -e ".[test]"  # plus pytest
    pytest                    # tests/

Nine runtime dependencies, each annotated in `pyproject.toml` with the module that needs it:
numpy, scipy, torch, openmm, pennylane, pennylane-lightning, fair-esm, scikit-learn,
biopython. `lightning.qubit` is the device `core.quantum` builds circuits for and the one
`tests/test_quantum.py` asserts probabilities against with `==`, so it is not an
interchangeable backend.

### Data this repository does not carry

`pdbs/` (61 structures) is tracked — it is the reporting set. The rest is large, ignored, and
rebuilt on demand by the module that owns it: `prots/` (5.5 GB) and `pdbs_ext/` (555 MB) from
the RCSB queries, `esm_cache.npz` (1.5 GB) from `esm_features.compute`, `peptide_db.npz` and
`fragment_db*.npz` from their `build` entry points, and the `*_models/` checkpoint
directories from `train_fold`. `.gitignore` names the owner of each, and carries the exact
`git show` commands for the four RCSB fetch scripts, which live in history rather than in the
tree. A checkpoint is valid only for the fold definition committed alongside it.


## Layout, and why it is flat

    FINDINGS.md          the research record — 72 findings, corrections ledger first
    README.md            this file
    pyproject.toml       the only config; every dependency annotated with its consumer
    core/                the consolidated production path (11 modules)
    <24 root modules>    the reference arm — see "Where everything is" above
    s5/ s7/ s8/ s9/      the reference arm's sprint packages, plus their evidence JSON
    tests/               the equivalence and scientific-invariant suite (9 files)
    verify/              standalone audits — AND an output directory of core/project.py
    pdbs/                the 61-structure reporting set
    results/             the two frozen manifests, and nothing else
    bench_results/       golden results, whitelisted past a blanket ignore
    core_cache/          a .gitignore sentinel; the cache itself is not tracked

The flatness is deliberate. Four things pin it, and each was verified rather than assumed:

**The root modules resolve data relative to themselves.** `peptide_db.py` sets
`BASE = dirname(abspath(__file__))` and assumes it *is* the repository root, then resolves
`pdbs/`, `results/benchmark_manifest.json`, `peptide_folds.json` and `peptide_clusters.json`
from it. Move the module and those resolve somewhere else.

**The sprint packages hard-code their hop count.** All eleven modules under `s5/ s7/ s8/ s9/`
do `sys.path.insert(0, dirname(dirname(abspath(__file__))))`. One level deeper and that points
at the wrong directory. `core/bench.py` also builds a subprocess preamble that inserts the
repository root and then imports `s9.final` — a shim in the parent would not reach the child.
And `FINDINGS.md` cites these paths by name several hundred times.

**`verify/` is not just scripts.** `core/project.py` *writes* its evidence there (`_outpath`)
and *reads* `project_inputs.json` and `project_equiv.json` back at runtime;
`tests/test_project.py` globs `verify/project_stability*.json`; `tests/test_integration.py`
imports three `verify/` scripts as modules. Renaming the directory is a code change.

> **The trap worth knowing about.** `peptide_folds.json`, `peptide_clusters.json` and
> `catrace_prior.npz` are **write-on-first-use caches**. If the file is missing they are not
> an error — they are silently re-derived. Moving them therefore does not crash, it re-pins
> the folds, and correcting the identity clustering once silently moved 13 benchmark targets
> and invalidated every model on disk. This is the most dangerous operation in the repository
> and it looks like a trivial tidy-up.

`pdbs/` is globbed and sorted by three modules. The sort makes the order stable *given the
same set of files*, but adding, removing or renaming one `.pdb` changes the set — and BLOSUM
similarity sums are small integers with large tie sets resolved by library enumeration order,
so an unstable order moves up to 47 of the 500 pool members. **Do not add or remove files
there.**

Everything else — 372 Python files across ten sprints, ~1 GB of caches, four sprints of
one-off experiment scripts — was deleted during the Sprint 11 consolidation and is in git
history. The full pre-consolidation tree is at commit `5fa05cd`;
`git log --diff-filter=D --name-only` finds any single file by name.

Untracked and deliberately so: `_archive/logs/` holds 99 console logs moved out of the sprint
packages during the 2026-09-04 consolidation (nothing deleted — see `_archive/README.txt`),
and `distogram_models_large.STALE-PRE-FOLD-REPIN-DO-NOT-USE/` is quarantined by rename so
that `FRAG_LARGE=1` fails loudly rather than loading a model that predates the fold repin.
