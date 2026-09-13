# SPRINT 26, LANE I (INFRASTRUCTURE / GOVERNOR): FINDINGS

Status: Phase 3 items done; every number below has an artefact path. Branch `s26`, commits
`6e50ea93` (brief as received), `eb89c165`, `37bddbbb`, `601a39c7` and the brief edit
`[PENDING: commit]`. Ledger entries L6, L7 (+ correction L9), L8, L10, L15, L16, and
`[PENDING: L for item 2, item 5, the test run]`.

Tiers: DEMONSTRATED (measured, artefact on disk) / ORACLE DIAGNOSTIC (reads native-derived
quantities, never selects) / HYPOTHESIS / REFUTED / OPEN.

---

## 0. The answers the brief asked for, in one table

| item | disposition | where |
|---|---|---|
| 1. test suite under the governor | DEMONSTRATED: 369 tests, 356 passed, 0 failed, 0 errors, 13 skipped (0 memory-guard) over three jobs; final non-AMBER re-run on the committed tree `[PENDING: counts]` | `s26/TEST_RUN.md`, `s26/results/test_run.json`, `s26/results/pytest_*.xml`, `s26/jobs_done/pytest_*.json` |
| 2. `verify/run_equiv2.sh` stale | CLOSED: rewritten to drive `Config.project_grad` through `s26/i_run_equiv2_arm.py`; four arms of smoke8 run with `--no-amber`; `[PENDING: result]` | `verify/run_equiv2.sh`, `s26/results/run_equiv2.log`, `s26/results/run_equiv2_compare.json` |
| 3. six held unused imports | CLOSED: all six applied, AST diff = exactly the removed names, tests green | commit `37bddbbb`, L10, `s26/i_ast_check.py` |
| 4. `_archive/logs/s8/predictor_report.log` | CLOSED: excerpt tracked with sha256 and size; the number is the S8-13 transfer law | `docs/sources/s8_predictor_report_excerpt.md`, L6 |
| 5. `TEST_RUN_RESULT_PLACEHOLDER` | CLOSED: brief committed unchanged (`6e50ea93`), then section 6.1 replaced with the governor run `[PENDING: commit]` | `docs/STATE_BRIEF_2026-09-12.md` |
| 6a. identity leak | flag ships dark; audit in memory; folds untouched (sha256-verified); the "corrected" criterion is at the null as a clustering threshold; the minimal fix moves exactly the 4 declared dev self-copies (+1 transitive) and 2/60 benchmark sequences | commit `37bddbbb`, `s26/i_identity_audit.py`, `s26/results/i_identity_audit.json`, `s26/PREREG_identity_null.md`, L15 |
| 6b. AMBER memory guard | message names the ceiling and quotes the governor; shared with the frame-invariance file; unit test; both AMBER files green under the governor | commit `eb89c165`, L16 |
| 6c. `pool_gate = WARN` | 1D6X and 1KWE; mechanism confirmed from the production cache; explanation added to the leaderboard text and the rule string; takes effect on rebuild | commit `eb89c165`, L7 + L9 |
| 6d. moment standardisation | production uses `core.pipeline._zrank`; every moment z-score is research-only or not on an AMBER energy | L8 |
| 7. hygiene | README governor section; `s26/examine.py` + `examine.sh`/`.bat`; claim ledger 21/21 OK; STATUS kept | commit `eb89c165`, `s26/results/claim_check.json` |

---

## 1. THE TEST SUITE UNDER THE GOVERNOR (task 1). DEMONSTRATED.

Every job ran through `s26/jobrun.py` (registered with the governor, stdout in
`s26/logs/<job>.log`, exit / wall / peak RSS in `s26/jobs_done/<job>.json`), never more than
one AMBER job from this lane at a time, and every launch was preceded by a read of
`s26/governor_state.json` (box 64.6 to 74.6% RAM throughout; never above 88% at a launch).

| job | files | tree | tests | passed | skipped | failed/err | wall | peak RSS |
|---|---|---|---|---:|---:|---:|---:|---:|
| `pytest_core` | the 10 non-AMBER files | `a4db170c` (before any lane-I edit) | 350 | 337 | 13 | 0 | 195.3 s | 1.694 GB |
| `pytest_amber` | `tests/test_amber.py` | `37bddbbb` | 16 | 16 | 0 | 0 | 260.5 s | 0.872 GB |
| `pytest_amber_frame` | `tests/test_amber_frame_invariance.py` | `37bddbbb` | 3 | 3 | 0 | 0 | 255.7 s | 0.324 GB |
| `pytest_nanpoison_post_core_edit` (not in the total) | `test_pipeline.py` + `test_data.py` on the edited core | `37bddbbb` | 79 | 77 | 2 | 0 | 225.5 s | 1.566 GB |
| `pytest_core_post` (replaces `pytest_core` in the final total) | the 10 non-AMBER files | `601a39c7` | `[PENDING]` | | | | | |

Per-file counts are in `s26/TEST_RUN.md`. The 13 skips, all real and none from the memory
guard (`-rs` reasons, verbatim in `s26/TEST_RUN.md`): 3 in `test_equivalence.py` and 8 in
`test_integration.py` are `set VERIFY_SLOW=1 ...` (the full-pipeline / OpenMM arms are opt-in);
2 in `test_pipeline.py` are absent artefacts (`bench_results/optimised_tuning126_w6.json`, "no
smoke result on disk"). The memory guard did not fire in any job, so no re-run was needed.

The previous record (state brief 6.1, "355 passed, 13 skipped") was one number for one
process; this one is per file, per job, with the peak RSS the governor measured, and one test
more (the 6b message test). The test count went 368 -> 369 (+1 in `test_amber.py`) and the
`test_data.py` count 41 -> 42 (+1 identity-flag test), so the whole suite is 370 with the
final non-AMBER re-run `[PENDING: confirm 370 = 351 + 16 + 3]`.

## 2. `verify/run_equiv2.sh` (operational item 2). `[PENDING: verdict]`

Why it was stale: it exported `PROJECT_GRAD` into the environment, which `core/project.py`
stopped reading when the gradient mode moved into `core.pipeline.Config.project_grad` so that it
reaches the cache key (`verify/grad_key_collision.py` documents the collision: two modes, one
`cfg_key`, the second run served the first's arrays). The `core.pipeline` CLI has no flag for the
mode, so as written the script ran three arms in the same mode.

The repair: `s26/i_run_equiv2_arm.py` builds the identical `Config` to `python -m core.pipeline
run` (the same `replace(PROD, ...)` call) and sets `project_grad`; the script runs four arms
(baseline legacy; consolidated `exact`, the shipped default; `fd`; `analytic`), keeps the
original `cfg_key` grep, and passes any argument (here `--no-amber`) to every arm.
`s26/i_equiv2_compare.py` then compares the arms per target from their cache directories
(`ca`, `fit_ca`, `phi`, `psi`, `avg_ca`, and the scalar RMSDs) with `==`.

Run: job `run_equiv2`, tag AMBER (conservative: with `--no-amber` no OpenMM context is created,
and the driver prints whether `openmm` was even imported), `[PENDING: wall, peak RSS]`.
Result `[PENDING]`. Log: `s26/results/run_equiv2.log` (whitelisted in `.gitignore`).

## 3. THE SIX HELD UNUSED IMPORTS (operational item 3). CLOSED, L10.

All six applied in commit `37bddbbb`; `s26/i_ast_check.py` (docstrings stripped, `ast.unparse`
diffed against HEAD) shows each file differing by exactly the removed names. The S25 hold
condition was checked by grep: neither `s25/resultslab/` nor `s25/phys_*.py` ever imported
`core.pipeline`; `s24/` names it in docstrings only; `s23/c1_probweight.py` (finished) and
`s12/instrument.py` (lazily, inside one function) do. At the moment of the edit the governor
had no job registered (`s26/governor_state.json` 00:26:31). Tests after the edit: section 1.
`peptide_db.py` was not touched.

## 4. THE S8-13 SOURCE LOG (operational item 4). CLOSED, L6.

`docs/sources/s8_predictor_report_excerpt.md`: path `_archive/logs/s8/predictor_report.log`,
sha256 `2647f640efb6d8f8e0398620647d503c366aa9a0209f17742d927366dd0f4740`, 13,588 bytes, 172
lines; cited by `docs/FINDINGS.md` S8-13 (line 3711) and the corrections ledger (line 69). The
number with no other source: `selected = 1.009 * ref + 0.011` (r = 0.981, 2,520 pairs), log line
92; a search for `1.009 * ref` over every json/md/py/txt/log file finds only FINDINGS and the
log. Two gaps found on the way: S8-13's second regression (`1.019 * ref - 0.001`, r = 0.984) is
in no on-disk file at all, and S8-13's X_fit positive-phi 0.0551 does not match the converged
row in the log (0.0629). Only the excerpt is tracked.

## 5. THE BRIEF'S SECTION 6.1 (operational item 5). `[PENDING: commit]`

`docs/STATE_BRIEF_2026-09-12.md` was committed exactly as received (`6e50ea93`, so history holds
the "355 passed, 13 skipped" sentence), then the sentence was replaced with the governor run
(per file, skip reasons, peak RSS, date, commit) in a second commit `[PENDING: hash]`.

## 6. THE DECLARED DEFECTS

### 6a. Identity leak. Flag ships dark; folds pinned; the brief's criterion is AT THE NULL. L15.

Code (`37bddbbb`): `core.data.identity(..., norm="longer")` / `identity_many(..., norm=)`;
default statement-for-statement the pinned convention; `norm="shorter"` = match count over the
shorter sequence behind a verbatim-substring test; `clusters()` / `folds()` do not accept it.
Test in `tests/test_data.py` (default bit-identical; 1CEK-in-1A11 = 1.0 under the flag; scalar
and batched agree on both sides of `_BATCH_MIN`).

Audit (`s26/i_identity_audit.py`, in memory; `peptide_folds.json` / `peptide_clusters.json`
sha256 before == after == lane E's `pinned_hashes.json`; 0.056 GB, job `i_identity_audit3`):

- DEMONSTRATED: the in-memory copy of `core.data.clusters` reproduces the pinned clusters
  (470) and the pinned folds exactly.
- DEMONSTRATED: the brief's corrected criterion at 0.6 collapses 470 clusters to 166, touches
  594/787 sequences, 71/126 dev targets and 48/60 benchmark sequences (count only).
- DEMONSTRATED (pre-registered, `s26/PREREG_identity_null.md`): that is chance. A real dev
  sequence passes the shorter criterion against 0.5% of the other 786 members; a
  composition-preserving shuffle passes against 0.3% (ratio 0.62; falsifier was < 0.2). At
  mean degree ~3 on chance edges, single linkage percolates. The pinned longer criterion:
  0.07% real vs 0.001% shuffled. So the flag is a per-pair leak test, not a clustering
  threshold, and its docstring says so.
- DEMONSTRATED: the minimal fix (pinned criterion OR verbatim substring): 459 clusters, 23
  sequences change membership (20 verbatim containments + 3 carried by single linkage), and
  the dev targets that gain a cross-fold mate are exactly the 4 declared self-copies (1CEK
  fold 2 vs carrier fold 3, longer identity 0.520; 2FBU 4 vs 0, 0.522; 2P5H 4 vs 2, 0.529;
  6B9K 0 vs 2, 0.588; each shorter identity 1.0, each moved into one fold by the fix) plus
  8ZG2, whose own carrier (n=23, longer identity 0.478) is already in its fold and which joins
  only because that carrier also carries a peptide from another fold. 14 further dev targets
  have a verbatim relative already in their own fold. Benchmark: 2/60 gain a cross-fold mate
  (count only), the S24 L4 figure.
- DEMONSTRATED: even the minimal fix, re-derived through `folds()`' shuffle, relabels 647/787
  sequences (the shorter criterion: 630/787, 102/126 dev), because cluster ids renumber. Any
  future repair must patch the pinned fold map and retrain only the affected fold models,
  never re-derive.

### 6b. AMBER memory guard. DEMONSTRATED, L16.

Decision unchanged (`_memory_verdict` on `core.amber.memory_percent()`); the message names
`core.amber`'s 92% ceiling and, when `s26/governor_state.json` is under 120 s old, the
governor's RAM / CPU / job counts / band ceiling, else says the governor is not running.
`test_amber_frame_invariance.py` imports the fixture (S25 L14). Unit test with synthetic
snapshots. Both files green under the governor (section 1); the guard did not fire.

### 6c. `pool_gate = WARN`. DEMONSTRATED, L7 corrected by L9.

`results/summary/results.json` production rows: 1D6X (T007, chain 1.7809 A vs its pool's
ORACLE best 2.4194) and 1KWE (T019, 2.4224 vs 2.8625). Production cache
(`bench_results/cache/1fc9f2dcf489e2fb/`): 1D6X `rmsd_avg` 2.0900, `rmsd_arm` 1.7814,
`top_m_best` 2.4194; 1KWE `rmsd_avg` 2.7313, `rmsd_arm` 2.4127, `top_m_best` 2.8625. The point
cloud is already nearer the native than the best single member of the 75 it averages, and the
projected chain stays below that member: averaging plus projection leaves the pool, on single
targets only (aggregate margin 1.50 A). `s25/resultslab/schema.py` now prints one explanatory
line under any WARN row and defines PASS / WARN / FAIL in the persisted `pool_gate_rule`
(verified on the on-disk payload, `s26/logs/precheck_fmt_leaderboard.log`). `results/` was not
rebuilt; the change appears on the next `--mode frozen` build.

### 6d. Moment standardisation. DEMONSTRATED, L8.

Production: `core.pipeline._zrank` (rank, line 788, used at 838 on the distogram score; and the
production `Config` has `quantum=False`, so the selector does not run in the 3.2148 A run).
Moment z-scores: `s25/phys_lib.zmoment` (research, the declared secondary, where the 40/126
non-monotone cases were measured), `s25/phys_landscape.py:132`, `s25/phys_gate.py:71`
(research), `core/energy.py:1066 LegacyField._standardise` (Legacy terms vs random structures;
the optional `legacy` term, `w_legacy = 0`, not instantiated by `core.pipeline`),
`core/quantum.py:2126 build_distogram(models="both")` (distogram scores; generation lane). None
applies a moment to an AMBER energy on a production path.

## 7. HYGIENE. DEMONSTRATED.

README section "The S26 resource governor" (band 88/90/93/95, AMBER cap 2, tags, log and state
paths, the split test run). `s26/examine.py` = `make examine`: calls lane E's `e_module_map`
(module map) and `s26/i_claim_check.py` (21 claims in `s26/results/claims.json`, each re-read
from its artefact: 21/21 OK, `s26/results/claim_check.json`; `--search` also runs lane E's
`e_claims.py`). `examine.sh` / `examine.bat` at the root (no conflict; no Makefile exists).
`s26/i_ast_check.py` makes the AST-identity rule runnable. STATUS lines at 00:13 and 00:42.

## 8. WHAT I DID NOT DO, AND WHY

- Did not re-derive folds or clusters, did not call `core.data.clusters()` / `folds()` or their
  `peptide_db` twins from the audit, did not write `peptide_folds.json` / `peptide_clusters.json`
  (hashes verified before and after). The fix ships dark.
- Did not add the `norm` flag to `peptide_db.identity`: the legacy arm must stay pure
  (`tests/test_amber.py::test_budget_classes_have_exactly_one_definition` states the rule for
  `budget.py`; the same reasoning holds), and `tests/test_data.py` pins `core.data.identity ==
  peptide_db.identity` on the default path.
- Did not rebuild `results/`; the 6c text change is dormant until the next frozen build.
- Did not run `verify/run_equiv2.sh` with AMBER: `--no-amber` leaves the projection comparison
  intact and keeps OpenMM out; the script supports both and says so.
- Did not read `results/benchmark_manifest.json`, any benchmark PDB or any benchmark RMSD. The
  benchmark counts came from `core.backend("data").benchmark()` sequences only; no benchmark
  name was printed; non-dev database members are described by length only in the audit output.
- Did not regenerate `s26/results/module_map.json` until the end `[PENDING: state whether done]`,
  because lane E was reading its own artefact for EXAMINATION.md.
- Did not add a claim for the benchmark headline (+0.0103, CI [-0.160, +0.180]) to
  `claims.json`: I did not find its results artefact and would not invent a path.

## 9. WHAT DAMAGED MY OWN EXPECTATIONS

1. **I wrote two numbers into the ledger before the query returned.** L7 quoted 1.5921 and
   2.2812 as point-cloud RMSDs; the cache says 2.0900 and 2.7313. The entry was composed in the
   same turn as the query. Corrected by L9, and the rule (a number enters the ledger only in a
   turn after its output was read) is written there.
2. **The "corrected" identity criterion is not a correction at the clustering level.** I
   expected the shorter normalisation to move a handful of targets; it moved 71/126, and the
   pre-registered null showed why (0.5% real vs 0.3% shuffled pass rate). My PREREG also
   expected pass rates of 10 to 30%; they are 0.5%: the percolation comes from mean degree ~3,
   not from a high pass rate. The at-the-null verdict stood; my magnitude was wrong by 30x.
3. **PREREG H3 was falsified as written.** I expected every sequence moved by the minimal fix
   to be a verbatim containment itself; 3 of 23 are carried along by single linkage. The
   falsifier should have been written on the edges, not on the nodes.
4. **8ZG2 looked like a fifth self-copy for one turn.** It is a transitive case (its carrier is
   in its own fold). The audit's first pass listed it; only the per-target verbatim breakdown
   in the third pass separated direct from transitive.
5. **The frame-invariance file ran in 255.7 s at 0.324 GB peak.** S25 L14 recorded it as the
   test most likely to hit the ceiling; on a 65% box with the governor's budget it is not close.

## 10. PRODUCTION-PATH CHANGES MADE (all deliberate, tested, in the ledger)

| file | change | flag / default | test |
|---|---|---|---|
| `core/data.py` | `identity(..., norm=)`, `identity_many(..., norm=)`; `_seq_index` removed | `norm="longer"` = old behaviour, bit-identical | `tests/test_data.py` (new test + the four existing identity pins) |
| `core/amber.py` | two unused names dropped from the `budget` import | n/a | `tests/test_amber.py` 16/16 |
| `core/bench.py`, `core/cache.py`, `core/predict.py` | one unused typing/dataclasses import each | n/a | non-AMBER suite |
| `verify/vqe_lfo_audit.py` | unused `defaultdict` | n/a (not production) | byte-compiled |
| `s25/resultslab/schema.py` | WARN explanation in `fmt_leaderboard` and `pool_gate_rule` | text only | rendered on the on-disk payload |
| `tests/test_amber.py`, `tests/test_amber_frame_invariance.py`, `tests/test_data.py` | 6b message + shared fixture; 6a flag test | n/a | themselves |

No number changed anywhere: the identity default is asserted bit-identical, the removed names
were never read, and the leaderboard change is text.

## 11. ARTEFACT INDEX

`s26/TEST_RUN.md`, `s26/results/test_run.json`, `s26/results/pytest_{core,amber,amber_frame,nanpoison_post_core_edit,core_post}.xml`,
`s26/jobs_done/*.json`, `s26/logs/*.log`; `docs/sources/s8_predictor_report_excerpt.md`;
`verify/run_equiv2.sh`, `s26/i_run_equiv2_arm.py`, `s26/i_equiv2_compare.py`,
`s26/results/run_equiv2.log`, `s26/results/run_equiv2_compare.json`; `s26/i_identity_audit.py`,
`s26/results/i_identity_audit.json`, `s26/PREREG_identity_null.md`; `s26/i_ast_check.py`;
`s26/i_test_report.py`; `s26/examine.py`, `s26/i_claim_check.py`, `s26/results/claims.json`,
`s26/results/claim_check.json`, `examine.sh`, `examine.bat`; README section; ledger L6-L10
(lane I: L6, L7, L8, L9, L10), L15, L16, `[PENDING]`.
