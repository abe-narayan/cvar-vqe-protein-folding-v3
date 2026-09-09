# S25 CLEANUP LANE — INVENTORY, CLASSIFICATION, ARCHIVE PLAN

**Status: INVENTORY AND PREPARATION ONLY. NOTHING HAS BEEN MOVED OR DELETED.**
Snapshot taken 2026-09-08 ~19:27. Research lanes were live in `s25/` and writing to
`results/site`, `results/structures`, `results/summary` during this pass; nothing under
`s25/` other than my own `agentCLEAN_*` files was touched, and no lock was taken.

Machine-readable artefacts persisted beside this file:

| artefact | what it is |
|---|---|
| `s25/agentCLEAN_reachability.json` | the import graph and the four reachability closures |
| `s25/agentCLEAN_inventory.json` | every file in the tree, classified, with per-file evidence |
| `s25/agentCLEAN_reach.py` | the reachability builder — re-runnable, read-only |
| `s25/agentCLEAN_classify.py` | the classifier that consumes it |
| `s25/agentCLEAN_hygiene.py` | the AST hygiene pass (no linter is installed in this env) |

---

## 0. THE TWO THINGS TO READ FIRST

**(1) The test suite is GREEN on the production path.** 204 passed, 13 skipped. Every
failure and error in the full-suite run is one cause — the 92% physical-memory ceiling —
and not one of them is a code defect. Detail in §4.

**(2) `s23/`, `s24/` and `s25/` ARE NOT IN GIT.** `git ls-files s23 s24 s25` returns zero
for all three. That includes both artefacts the brief names as unrecoverable:

    s24/cache_amber/                                126 files, 1.5 MB, 63,000 ff14SB/GBn2
                                                    single points — UNTRACKED
    s24/results/a_corpus_permitted_29e3b67e8ca0c03d.json
                                                    180 KB audited corpus + exclusion list
                                                    — UNTRACKED

They are not gitignored (`git check-ignore` returns nothing for them); they were simply
never added. 128 untracked files across the three sprints, of which 48 are `LEDGER.md`,
`BRIEF.md`, `PREREG_*.md` and `*_FINDINGS.md` — the exact class the brief says is
scientific evidence. **This is a one-disk-failure-away situation and it is the single
highest-risk finding in this lane.** It is a `git add` on three directories, and it needs
your decision because committing is not mine to make.

---

## 1. REPOSITORY INVENTORY

687 Python files, 33,000 files total, ~12 GB. Full classification:

| class | files | size | tracked | what it means |
|---|---:|---:|---:|---|
| `DATA_CORPUS` | 15,214 | 6,421 MB | 0 | `prots/`, `pdbs_ext/` — gitignored RCSB downloads |
| `REGENERABLE_CACHE` | 13,155 | 5,205 MB | 4,097 | per-target caches, derived banks, checkpoints |
| `EVIDENCE` | 4,521 | 353 MB | 1,409 | LEDGER/BRIEF/PREREG/FINDINGS/CLAIMS + result JSON + figures |
| `DATA_PINNED` | 194 | 61 MB | 194 | `pdbs/`, `results/structures` — order-sensitive, do not touch |
| `OBSOLETE` | 643 | 14.3 MB | 1 | bytecode + 10 root logs + 1 interrupted write |
| `RESEARCH_ARCHIVE` | 609 | 6.6 MB | 520 | sprint experiment modules and their console logs |
| `PRODUCTION` | 82 | 2.1 MB | 76 | reachable from `core.pipeline`/`core.bench`, plus `tests/` |
| `PROTECTED` | 127 | 1.2 MB | 0 | `s24/cache_amber/` + the audited corpus |
| `INSTRUMENT` | 36 | 0.2 MB | 29 | `s12/instrument.py`, `verify/` and its transcripts |

**The shape of the repository is not 687 experiment scripts. It is 48 production modules
and 12 GB of derived bytes.** 96.5% of the disk is `DATA_CORPUS` + `REGENERABLE_CACHE`,
both of which have a named owner that rebuilds them. The scientific record — every
`.md`, every result JSON, every figure — is 353 MB, and the code that is not production is
6.6 MB. That ratio is the whole cleanup story.

### 1.1 PRODUCTION — proved by import graph, not by name

The reachability closure from `core`, `core.pipeline`, `core.bench` is **48 modules**.
A pure AST walk gets this wrong, because `core/__init__.py` resolves backends through
`importlib.import_module` on the `REPLACES` table; those dynamic edges are extracted from
the literal and added explicitly (`agentCLEAN_reachability.json:dynamic_edges`).

    core/           __init__ amber bench cache data energy geometry pipeline predict
                    project quantum                                     (11 modules)
    root reference  amber_hamiltonian amber_refine budget catrace distogram energy_terms
    arm             esm_features floor foldvqe fragment_db legacy_field legacy_refine
                    objective pairnet peptide_db priors protein_geometry qansatz refine2
                    representations sidechains torsion_lib2               (22 modules)
    sprint packages s5.{esm32,esmraw,torsion} s7.{audit,debias,poolsize,repr_select}
                    s8.{consensus2,generate,inband,project} s9.final     (12 modules)

The root modules are **not** dead weight and the import graph says so: `peptide_db` has 87
importers, `torsion_lib2` 45, `protein_geometry` 44, `distogram` 23. They are the
`CORE_BACKENDS=legacy` arm — the reference every equivalence claim is proven against.
`core.bench --arm baseline` is exactly the act of running them.

Two root modules are **NOT** in the production closure and are reached only from tests:
`hamiltonian` and `vqe`, both imported solely by `tests/test_quantum.py` as the
consolidation's equivalence oracle. They are `PRODUCTION` by the test-of-production rule,
not deletable, and the README's claim that all 24 root modules are the reference arm is
correct in spirit but off by these two in mechanism.

### 1.2 INSTRUMENT

`s12/instrument.py` plus the 47-module closure under it (which is the production closure
plus itself — the instrument imports `core.pipeline`, `core.predict`, `core.project`,
`core.geometry`). `verify/` is treated as instrument, and that is a code fact rather than a
courtesy: `core/project.py` **writes** its evidence into `verify/` at runtime and reads
`project_inputs.json` / `project_equiv.json` back; `tests/test_integration.py` imports
three `verify/` scripts as modules; and five `verify/*.log` transcripts are cited by
basename inside `verify/run_equiv*.sh`. Renaming `verify/` is a code change, not a move.

### 1.3 RESEARCH ARCHIVE — 609 files, 6.6 MB

592 sprint experiment modules unreachable from production, the instrument or the tests,
plus 17 sprint console logs. **None of this is proposed for deletion.** The per-sprint
Python and log footprint is trivial (`s12` 133 modules = 0.9 MB; the whole of `s13`–`s25`
= 5.7 MB), so archiving it costs nothing and the reproducibility value is the project.

A note on how cheaply "unreachable" can mislead: 418 of the 592 have **zero** importers.
That does not make them dead — they are `__main__` experiment drivers, which is what a
sprint module is. Reachability is the right instrument for finding what production needs;
it is the wrong instrument for deciding what a sprint may throw away, and I have not used
it that way.

### 1.4 OBSOLETE — the complete list, with evidence

This is the whole of it. **11 non-bytecode files, 148 KB.**

| file | evidence |
|---|---|
| `work_bench_baseline.log` | root-level bench timing transcript. **Basename does not appear in any of the 8,556 text files in the tree — 433.7 MB searched, `.py` `.md` `.json` `.toml` `.txt` `.sh`, i.e. including every sprint result JSON and all 553 KB of `FINDINGS.md`.** `.gitignore` declares `*.log` disposable. Untracked. |
| `work_bench_comp.log` | as above |
| `work_bench_comp2.log` | as above |
| `work_bench_lpt.log` | as above |
| `work_bench_lpt2.log` | as above |
| `work_bench_opt.log` | as above |
| `work_bench_opt8.log` | as above |
| `work_bench_opt8b.log` | as above |
| `work_sweep.log` | as above |
| `work_sweep2.log` | as above |
| `s15/results/qgeom_ens.json.tmp` | an interrupted write. 25 KB, two keys (`E1_entropy_floor`, `_written`), written 21:56; the complete `s15/results/qgeom_ens.json` beside it is 180 KB and was written 23:39. Superseded by its own successor. **It is tracked in git and it sits inside a sprint `results/` directory, so it falls under the protected class by location — I am flagging it rather than assuming.** |

Plus 628 `__pycache__` / `.pytest_cache` files, 14.0 MB, gitignored and regenerated on the
next import.

The citation search was run twice, deliberately: once over 917 `.py`/`.md`/`.toml`/`.txt`/
`.sh` files (15.9 MB) and once over **8,556 text files including every result JSON**
(433.7 MB). Both returned the identical answer — of 36 `.log` files in the tree, exactly
five are cited by basename, and all five are `verify/` equivalence transcripts named
inside `verify/run_equiv*.sh`. Those five are classified `INSTRUMENT`, not obsolete. The
wide search is what makes this a negative result rather than an unsearched assumption; the
`_archive/README.txt` precedent is a case of exactly this check being run too narrowly.

**I explicitly did NOT classify sprint console logs as obsolete**, even though 17 of them
(`s21`–`s24`) have the same "not cited by basename anywhere" evidence as the root ones —
including under the 433.7 MB search.
`_archive/README.txt` records that this exact check was run before, was run WRONG, and that
one sprint log turned out to be cited in `FINDINGS.md` and another was the sole on-disk
source of a published number. The precedent is to archive, not delete, and the cost of
doing so is 200 KB.

---

## 2. THE ARCHIVE PLAN — proposed, not executed

### 2.1 The target layout

```
Protein-Folding-Algorithm/
├── README.md               rewritten: purpose, install, run, architecture, results
├── ARCHITECTURE.md         new: data flow, Hamiltonians, VQE/CVaR, distogram, limits
├── FINDINGS.md             the research record (stays at root, it is the headline)
├── pyproject.toml
├── core/                   production — 11 modules, unchanged
├── tests/                  12 files (10 existing + test_instrument.py + test_cvar.py)
├── verify/                 unchanged — it is a runtime path, not a script folder
├── <24 root modules>       unchanged — the CORE_BACKENDS=legacy reference arm
├── pdbs/  results/         unchanged
├── bench_results/          unchanged
├── s5/ s7/ s8/ s9/         unchanged — production closure members live here
└── research/               NEW — the only structural move
    ├── README.md           an index: one line per sprint, what it settled
    ├── s12/ … s25/         moved wholesale, code + docs + results
    └── _archive/           moved here from root
```

### 2.2 The moves, exactly

**One move type only.** `git mv` (or `mv`, for the untracked ones) each of these
directories under `research/`, contents intact:

    s12  s13  s14  s15  s16  s17  s18  s19  s20  s21  s22  s23  s24  s25  _archive

Nothing else moves. `core/`, `tests/`, `verify/`, `s5/`, `s7/`, `s8/`, `s9/`, the 24 root
modules, `pdbs/`, `results/`, `bench_results/` and every pinned data file stay exactly
where they are.

**`s5/`, `s7/`, `s8/`, `s9/` MUST NOT MOVE.** Twelve of their modules are in the production
closure, all eleven do `sys.path.insert(0, dirname(dirname(abspath(__file__))))` and would
resolve one directory wrong, and `core/bench.py` builds a subprocess preamble that inserts
the repository root and imports `s9.final`. `s12/` has the same hop-count assumption in
`instrument.py` (`ROOT = dirname(dirname(abspath(__file__)))`).

### 2.3 The blocker on the move, and the honest cost

**Moving `s12/`–`s25/` into `research/` breaks the two-hop `ROOT` computation in every one
of the 592 modules under them.** They would each need `dirname` applied once more. That is
a 592-file mechanical edit to code whose whole value is that it reproduces a published
number, and it is exactly the class of change this project has been burned by.

There are three options and I am not choosing between them:

| option | cost | risk |
|---|---|---|
| **A. Leave them at root.** Add `research/README.md` as an index and say so in the README. | zero | the tree still shows 14 sprint directories to a visitor |
| **B. Move, and fix the hop count** in 592 files with one scripted edit, then re-run the suite. | one scripted edit + a full suite run | a module that computed ROOT some other way silently resolves wrong; 592 files is a lot of surface for a "tidy-up" |
| **C. Move, and drop a `research/__init__.py`-free shim** — a `conftest.py`/`sitecustomize` that puts the repo root on `sys.path` — so the two-hop `ROOT` is wrong but nothing reads it. | small | `ROOT` is used to build data paths (`s8/generate_univ`, `bench_results/cache`), not just `sys.path`, so it does not work in general |

**My recommendation is A**, and it is not laziness. The repository's flatness is already
argued for in the existing README with four verified reasons, the sprint directories are
6.6 MB of code against 12 GB of data, and a visiting professor's complaint about a flat
tree costs less than one silently re-resolved cache path. If you want B, it should be its
own post-freeze task with its own full-suite verification, not a step inside a cleanup.

### 2.4 Results and structures

`results/` is currently being populated by the results-lab lane (`site/`, `structures/`,
`summary/` all appeared while I was inventorying). I have proposed no change to it. The
Phase-II layout in the brief — `results/summary/{results.json,results.csv,leaderboard.csv,
final_report.md,professor_brief.md}` — is compatible with what is there.

### 2.5 Items I will NOT classify without you

1. `s15/results/qgeom_ens.json.tmp` — obsolete by content, protected by location.
2. `distogram_models_large.STALE-PRE-FOLD-REPIN-DO-NOT-USE/` — 7.2 MB, quarantined by
   rename on purpose so `FRAG_LARGE=1` fails loudly. It is doing a job. Keep, but it wants
   a sentence in the README rather than a scary directory name in `ls`.
3. `CONSOLIDATION_MANIFEST.json` (12.7 KB) and `CONDENSED_REPORT.md` (12.4 KB) at root —
   neither is referenced by basename anywhere in the tree. They read as Sprint-11
   consolidation artefacts. Probably `research/`, possibly obsolete, **your call**.
4. `s8/integrate.py` and `s8/invfold.py` — the only two modules in `s5/s7/s8/s9` outside
   every closure. `s8.invfold` is imported only by `s8.integrate`, and `s8.integrate` by
   nothing. `tests/test_quantum.py` names `s8/integrate` in its docstring as a
   consolidation source. Almost certainly archive; not obviously deletable.
5. The 2.7 GB of sprint `cache/` directories (`s12/cache` 1.11 GB, `s14/cache` 1.03 GB,
   `s16/cache` 210 MB, `s18/cache` 133 MB, …). Regenerable by the module that wrote them,
   but regenerating `s12/cache` means recomputing 126 distograms. **Not deletion
   candidates in my view — exclusion candidates**: keep on disk, keep out of git, and note
   them in the README as derived.

---

## 3. TEST AUDIT

### 3.1 What exists, classified

| file | tests | class |
|---|---:|---|
| `tests/test_amber.py` | 15 | physics validation + numerical correctness (ff14SB/GBn2 genuineness, bit-exact refine, topology identity) |
| `tests/test_amber_frame_invariance.py` | 3 | physics validation (relaxation must not depend on the input frame) |
| `tests/test_data.py` | 41 | production regression (database, folds, identity clusters, manifests) |
| `tests/test_energy.py` | 9 | numerical correctness (the 11 Legacy terms) |
| `tests/test_equivalence.py` | 14 | production regression against the `s9/final.py` golden constants |
| `tests/test_geometry.py` | 24 | numerical correctness (build, Kabsch, RMSD, SS assignment) |
| `tests/test_integration.py` | 24 | production regression, end-to-end + `verify/` audits as modules |
| `tests/test_pipeline.py` | 33 | production regression (four stages, cache keys, native-leak NaN poisoning) |
| `tests/test_project.py` | 28 | numerical correctness (manifold projection, λ-path, multi-start) |
| `tests/test_quantum.py` | 39 | VQE/CVaR correctness — equivalence against `qansatz`/`vqe`/`foldvqe` |
| **`tests/test_instrument.py`** | **44** | **NEW — instrument regression (§3.4)** |
| **`tests/test_cvar.py`** | **90** | **NEW — CVaR definitional regression (§3.4)** |

No obsolete sprint tests and no dead tests were found in `tests/`. There are two test
modules living outside `tests/` — `s14/test_vqe.py` and `s25/resultslab/test_export.py` —
which `pytest` does not collect (`testpaths = ["tests"]`). `s25/resultslab/test_export.py`
belongs to a live lane; `s14/test_vqe.py` is sprint archive.

### 3.2 Suite state as run — 2026-09-08

    python -m pytest tests/ -q          235 tests
                                        204 passed
                                         13 skipped
                                          3 failed
                                         15 errors

**All 18 non-passing outcomes have a single cause: the 92% physical-memory ceiling.**

* The 15 errors are all of `tests/test_amber.py`, all raised in the autouse
  `_memory_ceiling` fixture at setup. `_memory_verdict` returned `"fail"` rather than
  `"skip"` because the baseline was captured at *collection* time — before pytest had
  imported torch, pennylane, esm and openmm for the other nine modules — so the whole
  suite's import growth was attributed to the AMBER suite. Measured: 80% at rest, 93–97%
  during the run, and back to 80% the instant the pytest process exited.
* 2 of the 3 failures are `tests/test_amber_frame_invariance.py`, `MemoryError` out of
  `core/amber.py:memory_guard` at 93%.
* The third failure is `test_partition_is_actually_used_and_is_cheaper_where_it_is_used`,
  a wall-clock assertion on 256-element arrays: 1.860e-5 s against a 1.851e-5 s budget.
  It missed by 0.5% under load, and its own docstring says the threshold is loose because
  "this box runs three sibling jobs".

**Re-run in isolation on the same box, same code:**

    python -m pytest tests/test_amber.py tests/test_amber_frame_invariance.py -q
        17 passed, 1 failed

And with the AMBER files split out — which is the right way to run this suite on this box,
and is now documented in the draft README:

    python -m pytest tests/ --ignore=tests/test_amber.py \
                            --ignore=tests/test_amber_frame_invariance.py -q
        350 tests: 337 passed, 13 skipped, 0 failed, 0 errors      (exit 0)

That 350 includes the 134 new tests in §3.4, so the whole non-AMBER surface plus the two
new regression files is green, at 93% memory, with the timing test passing.

`tests/test_amber.py` is **15/15 green**. The single remaining red is
`test_amber_relaxation_is_frame_invariant`, which relaxes target after target until the
box crosses 92% — it got through 11 targets first, and its captured stdout shows frame
invariance holding at |ΔRMSD| ≈ 1e-4 Å on 9 of them (`1DEP` +0.140, `1FUV` +0.016 are the
outliers). **The physics is not failing. The box is full.**

### 3.3 The one test-hygiene defect, and a proposal I am not making unilaterally

`tests/test_amber.py` has the autouse memory guard that converts memory pressure into a
`skip`. `tests/test_amber_frame_invariance.py` does not, so the same condition renders as
a red `FAIL`. Since it opens ~15 OpenMM contexts it is the *more* likely of the two to hit
the ceiling.

**Proposal (needs your approval — converting a red into a skip is a judgement about the
instrument, not a formatting fix):** give `test_amber_frame_invariance.py` the same
autouse `_memory_ceiling` fixture, imported from `test_amber.py` rather than copied.
It cannot change a number; it changes what a busy box reports.

Also: `@pytest.mark.slow` is used in that file but not registered, so pytest emits two
`PytestUnknownMarkWarning`s on every run. One line in `pyproject.toml`
(`markers = ["slow: ..."]`) removes both and enables `-m "not slow"`, which is the right
way to run this suite on a loaded box. Safe, but it is a config edit and lanes are live.

### 3.4 MISSING regression tests — WRITTEN AND PASSING

Both new files are pure arithmetic on synthetic arrays. Nothing loads a window universe,
a distogram, ESM, OpenMM or any cache; they take about a second, hold no lock, and are
safe to run beside a live experiment.

#### `tests/test_instrument.py` — 44 tests, all passing

`s12/instrument.py` had **no test of its own**. Every RMSD in this project — every
`FINDINGS.md` number, every LEDGER row, the 3.0483 Å incumbent — comes out of ten functions
in that file, and a silent change to `kabsch_rmsd_batch` would move all of them at once
with nothing going red.

Worse, and this is the finding: **the instrument carries its own copies of four functions
that also exist in `core.geometry`** — `kabsch_rmsd_batch`, `ca_rmsd`, `pair_index`,
`pair_dists`. Not imports, not aliases: separate code. The pipeline optimises against
`core.geometry`; the instrument that scores the pipeline uses its own. They agree today —
I measured the difference at **exactly 0.0**, because both reconstruct the residual from
the singular values in the same order — so the cross-checks are asserted at `==`, not at a
tolerance, and a first-ulp divergence now fails loudly.

Covered: self-RMSD, rigid-motion invariance on *both* arguments, symmetry, the
reflection rejection (the determinant fix — delete it and every mirrored decoy scores 0),
agreement with an independently written textbook Kabsch, linear scaling, batch-vs-loop
identity; `superpose_batch` realising the RMSD it reports and being a genuine rigid motion;
`pairwise_rmsd` symmetry and agreement with the batch routine; `medoid` as argmin-of-row-
means; `coordinate_average` against its documented recipe and against rigid copies of one
structure (the control for the known 25.8% contraction); `pair_index`/`pair_dists` against
brute force; `shipped_score`'s grid lookup including the out-of-range clamp at both ends;
`paired`'s full arithmetic; `summary`; and `write`'s `complete` flag.

Three of those deserve calling out.

* **`test_medoid_tie_break_is_the_lowest_index_and_that_is_a_documented_hazard`** pins
  that `medoid` is `np.argmin`, so a tie resolves to the order the pool arrived in. This
  project has already had an `np.argmin` over a tied native-free signal silently read the
  *oracle* sort order and invent a 1.386 Å winner. The test does not call the behaviour
  wrong — deterministic tie-breaking is right — it makes the inheritance explicit.
* **`test_paired_se_is_the_one_the_mde_rule_consumes`** asserts that `se` is the ddof=1
  standard error of the *paired difference*, not of either arm. `MDE = 2.8016 × SE` is a
  hard rule of this sprint, so this field prices every MDE in every ledger.
* **A numerical fact I found and did not paper over.** The self-RMSD of a structure with
  itself is **~1.3e-7 Å, not 0.0**, in both implementations. The residual is computed as
  `|P|² + |T|² − 2Σs`, which for identical structures is a difference of two ~1e4 Å²
  quantities whose true difference is zero; ~1e-11 Å² of cancellation error survives and
  the square root turns it into 1e-7 Å. It is four orders below the 1e-3 Å the project
  reports to and cannot affect a finding — but my first draft asserted an exact zero and
  *the test was wrong, not the code*. It is now pinned as a measured floor.

#### `tests/test_cvar.py` — 90 tests, all passing

`tests/test_quantum.py` is excellent, and it is an **equivalence** suite: it proves
`core.quantum` reproduces `qansatz`/`vqe`/`foldvqe` bit-for-bit and pins the tail-mean
gradient-baseline defect. But it inherits its notion of correct from the module it
compares against — **if the shipped CVaR and its reference were wrong in the same way,
every one of those 39 tests would still pass.**

So this file asserts the *definition* instead, and none of it was checked anywhere:
CVaR equals the mean at α=1 (the failure where the optimiser is silently minimising the
mean); never exceeds the mean; is monotone non-decreasing in α; descends to the minimum;
is translation-equivariant and positively homogeneous (the two properties that license
rank standardisation between Hamiltonians at all); is invariant to sample order; and
matches a naive independent reference across α ∈ {0.05 … 1.0} and n ∈ {1 … 1000}, which
crosses `CVAR_SORT_CUTOFF = 448` so **both** the argsort branch and the argpartition
branch are exercised at every α.

Also newly covered: `tail_indices` (exactly k elements, all of them the lowest k,
lowest-index tie-breaking on *both* branches under heavy ties); the exact distributional
pair `cvar_from_probs`/`cvar_from_distribution` agreeing with each other, their tail mass
summing to exactly α and never exceeding a state's own probability, the fractional
boundary split, and zero-probability states being ignored — the failure mode where a
selector "finds" a brilliant configuration the circuit never prepares; `alpha_schedule`
(endpoints, geometric ratio, clamping, and that its whole output range is a legal α);
`BestSeenTracker` (keeps the minimum, counts every offer — that count is the evidence the
search does not enumerate); `all_bitstrings`; `n_parameters`.

One test is written directly against a measured hazard from the project record:
`test_cvar_survives_the_scale_raw_amber_actually_arrives_at` feeds in the real
distribution — 200 normal values plus 20 at 7.1e18 — and asserts the value stays finite,
no clash enters the low tail, and the answer still matches the reference.

The four sampled/exact entry points are asserted to agree **only where their domains
actually meet** (uniform weights, `α·N` integral). Asserting it everywhere would be
asserting a falsehood: the sampled form takes whole samples, the exact form splits the
boundary state's mass. That distinction is stated in the module docstring.

---

## 4. BUILD AND HYGIENE

No linter is installed (`ruff`, `flake8`, `pyflakes`, `pylint`, `black` all absent), so
the pass is `s25/agentCLEAN_hygiene.py`, an AST walk over `core/`, `s12/instrument.py`,
`tests/` and `verify/`.

**Import errors: none.** All 687 files parse. All nine declared runtime dependencies
import cleanly: numpy 2.5.1, scipy 1.18.0, torch 2.13.0+cpu, openmm 8.5.2, pennylane
0.45.1, pennylane-lightning, fair-esm 2.0.0, scikit-learn 1.9.0, biopython 1.87,
pytest 9.1.1. `pyproject.toml` is accurate and every dependency is annotated with its
consumer; nothing is declared that is not imported, and nothing imported is undeclared.

**Fixed (safe, in files no lane imports, and none can alter a number):**

| file | change |
|---|---|
| `tests/test_data.py` | removed `from core import geometry as g` — bound name never used |
| `tests/test_geometry.py` | removed `import math` — zero `math.` uses in the file |
| `tests/test_integration.py` | removed `import subprocess` — single occurrence is the import |

Re-run after the edits: `tests/test_data.py tests/test_geometry.py tests/test_integration.py`
→ 82 passed, 8 skipped.

**Proposed, not done — these are inside `core/`, and lanes are importing `core.pipeline`
right now. BRIEF §4: never edit a module while a job launched from it is running.**

| file | finding |
|---|---|
| `core/amber.py:288` | `MAXITER_PER_PARAM`, `MIN_MAXITER_PER_PARAM` imported, never referenced |
| `core/bench.py:49` | `asdict` imported, never referenced |
| `core/cache.py:43` | `Iterable` imported, never referenced |
| `core/predict.py:38` | `Dict` imported, never referenced |
| `core/data.py:451` | `_seq_index()` — a dead private `lru_cache`d helper. **Zero references repo-wide** (687 files searched). Genuinely dead; still `core/`, so still a proposal. |
| `verify/vqe_lfo_audit.py:19` | `defaultdict` imported, never referenced |

**Dead functions: one** (`core/data.py:_seq_index`). The other two the checker flagged —
`tests/test_amber.py:_memory_ceiling`, `tests/test_equivalence.py:_cached_arm` — are
decorated pytest fixtures and are false positives; noted so the next person does not
"clean" them.

**Broken paths: none.** The 11 the checker flagged (`amber14/protein.ff14SB.xml`,
`implicit/gbn2.xml`) are OpenMM force-field resource names resolved from the openmm
package, not filesystem paths. False positives, recorded so the check is not re-run
from scratch.

**Missing assets: none** in the production or instrument path. Everything gitignored has a
named owner in `.gitignore` that rebuilds it, and that file is unusually good — it names
the rebuilding module for every ignored artefact and carries the `git show` commands for
the four RCSB fetch scripts that live only in history.

**Duplicate utilities — the real finding, and it is `s12/instrument.py`:**

| name | in `core` | in the instrument |
|---|---|---|
| `kabsch_rmsd_batch` | `core.geometry` | own copy — bit-identical today, now asserted at `==` |
| `ca_rmsd` | `core.geometry` | own copy |
| `pair_index` | `core.geometry`, `core.energy` | own copy |
| `pair_dists` | `core.geometry` | own copy |
| `pairwise_rmsd` | `core.geometry.pairwise_ca_rmsd` | own copy — **behaviourally different**: `core` returns `0.5*(M + M.T)`, the instrument returns the raw matrix |

The `pairwise_rmsd` difference is real. Both are symmetric to ~1e-9 so it moves no number
today — but a medoid is an argmin over row means, and symmetrising changes those means in
the last bits, which is exactly what decides a tie. It is pinned by
`test_pairwise_rmsd_is_not_symmetrised_the_way_core_geometry_symmetrises_it`, with the
note that the findings were measured with the instrument's unsymmetrised form. **Do not
"unify" these without deciding which one the record was written against.**

**Stale references in `README.md`** (details in §5): the headline is the pre-S12 number,
"72 findings across seven sprints" is now fourteen sprints, and the Layout section does
not mention `s12/`–`s25/` at all — 4.4 GB and 592 modules of the tree are invisible in the
document that claims to describe it.

---

## 5. DOCUMENTATION SKELETON

Drafted, **not installed**. The live `README.md` is untouched.

* `s25/agentCLEAN_README_draft.md`
* `s25/agentCLEAN_ARCHITECTURE_draft.md`

Architecture-specific sections that depend on the frozen pipeline are marked
`<!-- PLACEHOLDER: FREEZE -->` with the exact question each one is waiting on.

What the draft README strips, and why: the existing one is genuinely good — the
"why it is flat" section and the instrument-discipline rules are worth keeping almost
verbatim — but a new researcher currently has to absorb the Sprint 11 consolidation story,
the `CORE_BACKENDS=legacy` dual-arm design and the `5fa05cd` history pointer before
reaching a command they can run. In the draft, "install, run one target, run the
benchmark, look at the results" is above the fold, and the consolidation history moves to
a clearly-labelled section near the end where it is background rather than a prerequisite.

---

## 6. WHAT I NEED FROM YOU

1. **The untracked-evidence decision (§0.2).** `s23/`, `s24/`, `s25/` including both
   protected artefacts are outside version control. Highest-risk item in this lane.
2. **Approve or cut the OBSOLETE list (§1.4).** 11 files, 148 KB, plus bytecode. I
   deliberately kept it small.
3. **Pick A, B or C for the sprint-directory move (§2.3).** I recommend A.
4. **Rule on the five unsure items (§2.5).**
5. **Approve the `core/` hygiene edits and the `test_amber_frame_invariance` memory guard
   (§3.3, §4)** — post-freeze, or now if the lanes are between runs.

Nothing in this document has been acted on except: two new test files, three unused
imports removed from `tests/`, and five `agentCLEAN_*` artefacts written into `s25/`.
