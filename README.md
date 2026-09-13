# Peptide structure prediction, and the instrument that measures it

A structure predictor for short peptides (roughly 8-20 residues), and - more to the point -
the measuring instrument built around it. The predictor retrieves candidate backbones from a
library, filters them with a learned distogram, synthesises one structure from the survivors,
and relaxes it under AMBER. The instrument prices every one of those stages against a
held-out benchmark, and its findings are the actual output of the project.

**Read [`docs/FINDINGS.md`](docs/FINDINGS.md) first.** It is the consolidated research record: 72
findings across seven sprints, with a corrections ledger at the top because several of them
withdraw or revise earlier ones. The headline is there and it is a negative result:

> On the 60-target held-out benchmark, the synthesis architecture scores **2.9610** against
> the shipped distogram baseline's **2.9507** - a paired gain of **+0.0103, 95% CI
> [-0.1596, +0.1803]** - against a tuning-instrument gain of -0.250. The mean transferred;
> the effect did not. **There is no validated accuracy improvement over the baseline.**

The architecture is physically better-behaved than the baseline (exact 3.804 Å virtual
bonds, native-like torsion statistics, strain removed) and no more accurate. The binding
constraint is *recognition* - ranking within the near-native band - and S9-8 measured that
closure as informational rather than architectural.


## Where everything is

| if you are looking for | it is in |
|---|---|
| **production code** | `core/` - 11 modules, one implementation per operation |
| **the four pipeline stages** | `core/pipeline.py` |
| **quantum code** | `core/quantum.py`; reference arm `foldvqe.py` `qansatz.py` `vqe.py` `objective.py` `hamiltonian.py` `budget.py` `refine2.py` |
| **retrieval** | `core/data.py`, `s5/lib.py`, `s7/debias.py`, `s7/poolsize.py` |
| **physics** | `core/energy.py`, `core/amber.py`; reference arm `energy_terms.py` `legacy_field.py` `legacy_refine.py` `amber_refine.py` `amber_hamiltonian.py` `protein_geometry.py` |
| **tests** | `tests/` - 12 files, equivalence and scientific invariants |
| **experiments** | `s5/ s7/ s8/ s9/` - sprint packages, still load-bearing (see below) |
| **the research log** | `s12/` … `s25/` - one directory per sprint: scripts, their results JSON, and a findings write-up. Not imported by `core/`. |
| **reproducibility evidence** | `verify/*.json`, `bench_results/*.json`, and the sprint packages' own `*.json` |
| **verification scripts** | `verify/*.py`, `verify/run_equiv*.sh` |
| **data** | `pdbs/` (61 structures), `results/` (frozen manifests, exported structures, the results-lab summary), and the pinned caches at root |
| **the research record** | `docs/FINDINGS.md` |

The 24 modules at the repository root are the **reference arm** - the implementations the
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
`refine2.py` are reached by `core/` on **both** arms - they are production, not archive.


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

`core/` is the consolidated production closure - one implementation per operation. Each
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
into the cache key - so an optimised result can never be served out of a baseline cache.

    python -c "import core; print(core.backend_report())"
    CORE_BACKENDS=legacy python -c "import core; print(core.backend_report())"

The second prints `amber_refine, peptide_db, energy_terms, protein_geometry, legacy_field,
s7.audit, distogram, foldvqe, s8.project`. **Those nine modules and their transitive closure
are why the root modules and the sprint packages survived consolidation.** They are not
leftovers; they are the arm every equivalence claim in this project is proven against, and
deleting them would not tidy the repository, it would make those claims unreproducible.

The full pre-consolidation tree - all 372 Python files, ten sprint packages, every one-off
experiment - is at commit **`5fa05cd`**. `git checkout 5fa05cd` recovers it whole;
`git log --diff-filter=D --name-only` finds any single file by name.


## Running a full-quality experiment

    python -m core.pipeline run --manifest tuning126 --workers 6
    python -m core.pipeline stats --manifest tuning126

Manifests: `smoke8`, `smoke24`, `tuning126`, `dev24`, `benchmark60`.

To measure rather than just run - per-stage clocks, both arms, a real speedup:

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
that can change a number - K, M, penalty, λ, the AMBER schedule, the fold count, the pair
separation, the optimiser iteration cap, the tie-break rule, and the live backend set.

**So resuming is just re-running the same command.** A run that dies at 99% resumes at 99%.
A run with a changed parameter gets a different key and cannot silently inherit the old
answer. `stats` reports how many targets were resumed, and the results file records it,
because a resumed run is not a cold timing measurement.


## The instrument discipline

These rules are why the numbers in `docs/FINDINGS.md` mean anything. Breaking one does not produce
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

`pdbs/` (61 structures) is tracked - it is the reporting set. The rest is large, ignored, and
rebuilt on demand by the module that owns it: `prots/` (5.5 GB) and `pdbs_ext/` (555 MB) from
the RCSB queries, `esm_cache.npz` (1.5 GB) from `esm_features.compute`, `peptide_db.npz` and
`fragment_db*.npz` from their `build` entry points, and the `*_models/` checkpoint
directories from `train_fold`. `.gitignore` names the owner of each, and carries the exact
`git show` commands for the four RCSB fetch scripts, which live in history rather than in the
tree. A checkpoint is valid only for the fold definition committed alongside it.


## The S26 resource governor

This is a 16.75 GB, 8-core box shared with the user's own sessions; the baseline load with no
campaign work is about 67% RAM. Sprint 26 added a small governor so that several agents can run
heavy jobs without pushing the box past its ceiling or running more than two OpenMM jobs at once.
Three scripts, all under `s26/`, pure `psutil`:

| script | what it does |
|---|---|
| `s26/governor.py` | samples every 5 s; **above 93% RAM or CPU suspends the newest registered job**; above 95% for 15 s kills the newest (CTRL_BREAK first so it can checkpoint, 25 s grace, then terminate) and puts its spec back on the queue; resumes suspended jobs once the box is back under 90%; **launches the next queued job when the box has sat under 88% for 60 s**; caps AMBER-tagged jobs at two (a third is suspended on sight until a slot frees). `--once` prints one sample, `--status` the last snapshot. |
| `s26/jobrun.py` | runs one command as a registered job: waits for headroom (and an AMBER slot if tagged AMBER), writes `s26/jobs/<name>.json`, streams stdout+stderr to `s26/logs/<name>.log`, samples the process tree's RSS, and on exit writes `s26/jobs_done/<name>.json` with the exit code, wall time and **peak RSS** (the honest memory figure to quote). |
| `s26/enqueue.py` | puts a job spec on `s26/queue/` with a priority; the governor launches it through `jobrun.py` when the band allows. |

The band: **88% low water** (launch queued work), **90%** (resume suspended work below it),
**93% ceiling** (suspend), **95% for 15 s** (kill and requeue). Tags: `CPU` (default), `AMBER`
(anything that imports OpenMM or calls `core.amber`; at most two at once), `ESM` (loads the ESM
bank), `TEST` (pytest). Live state is in `s26/governor_state.json` (rewritten every sample; read
it before launching anything heavy); every action is appended to `s26/governor.log`, which is
whitelisted past the blanket `*.log` ignore because it is a deliverable.

    python s26/governor.py                                   # foreground; Ctrl-C stops it
    python s26/jobrun.py --agent I --tag CPU --name probe --est-ram 0.5 -- python some_script.py
    python s26/enqueue.py --agent P --tag AMBER --name relax --priority 20 --est-ram 2.5 -- python s26/p_relax.py

The test suite runs under it split three ways so at most one AMBER job is live per lane: the
non-AMBER files as one `TEST` job, then `tests/test_amber.py` and
`tests/test_amber_frame_invariance.py` as two sequential `AMBER` jobs. The record of each run
(per-file counts, every skip reason, peak RSS) is `s26/TEST_RUN.md` / `s26/results/test_run.json`.
Both AMBER files carry an autouse memory guard: above `core.amber`'s 92% ceiling they skip (or
fail, if this suite's own working set is what filled the box) with a message that names the
ceiling and quotes the governor's last reading when one is running.

`python s26/examine.py` (or `examine.sh` / `examine.bat` at the root; there is no Makefile)
regenerates the module map `s26/results/module_map.json` and re-reads every number in
`s26/results/claims.json` from the artefact it is claimed from, reporting `OK`, `MISMATCH` or
`ABSENT`.


## Layout, and why it is flat

    README.md            this file
    ARCHITECTURE.md      the frozen production spec
    docs/                the research record, four files:
                           FINDINGS.md          72 findings, corrections ledger first
                           CONDENSED_REPORT.md  the sprint 5-13 summary
                           consolidation-2026-09-04.json   what the last structural pass did
                           quarantine-synthetic-proofbuild.md  a leakage post-mortem
    pyproject.toml       the only config; every dependency annotated with its consumer
    core/                the consolidated production path (11 modules)
    <24 root modules>    the reference arm - see "Where everything is" above
    s5/ s7/ s8/ s9/      the reference arm's sprint packages, plus their evidence JSON
    tests/               the equivalence and scientific-invariant suite (12 files)
    verify/              standalone audits - AND an output directory of core/project.py
    pdbs/                the 61-structure reporting set
    results/             frozen manifests, exported PDBs, and the results-lab summary
    s12/ … s25/          the research log - one directory per sprint, not imported by core/
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
repository root and then imports `s9.final` - a shim in the parent would not reach the child.
And `docs/FINDINGS.md` cites these paths by name several hundred times.

**`verify/` is not just scripts.** `core/project.py` *writes* its evidence there (`_outpath`)
and *reads* `project_inputs.json` and `project_equiv.json` back at runtime;
`tests/test_project.py` globs `verify/project_stability*.json`; `tests/test_integration.py`
imports three `verify/` scripts as modules. Renaming the directory is a code change.

> **The trap worth knowing about.** `peptide_folds.json`, `peptide_clusters.json` and
> `catrace_prior.npz` are **write-on-first-use caches**. If the file is missing they are not
> an error - they are silently re-derived. Moving them therefore does not crash, it re-pins
> the folds, and correcting the identity clustering once silently moved 13 benchmark targets
> and invalidated every model on disk. This is the most dangerous operation in the repository
> and it looks like a trivial tidy-up.

`pdbs/` is globbed and sorted by three modules. The sort makes the order stable *given the
same set of files*, but adding, removing or renaming one `.pdb` changes the set - and BLOSUM
similarity sums are small integers with large tie sets resolved by library enumeration order,
so an unstable order moves up to 47 of the 500 pool members. **Do not add or remove files
there.**

Everything else - 372 Python files across ten sprints, ~1 GB of caches, four sprints of
one-off experiment scripts - was deleted during the Sprint 11 consolidation and is in git
history. The full pre-consolidation tree is at commit `5fa05cd`;
`git log --diff-filter=D --name-only` finds any single file by name.

Untracked and deliberately so: `_archive/logs/` holds 99 console logs moved out of the sprint
packages during the 2026-09-04 consolidation (nothing deleted - see `_archive/README.txt`),
and `distogram_models_large.STALE-PRE-FOLD-REPIN-DO-NOT-USE/` is quarantined by rename so
that `FRAG_LARGE=1` fails loudly rather than loading a model that predates the fold repin.
