# SPRINT 26, LANE I (INFRASTRUCTURE / GOVERNOR): FINDINGS

Status: Phase 3 complete. Every number below has an artefact path. Branch `s26`. Lane-I commits:
`6e50ea93` (the brief as received), `eb89c165`, `37bddbbb`, `601a39c7`, `83305549` (the brief's
6.1), `52d15a2e` (the test record), and the closing commit that carries this file. Ledger entries
(lane I): L6, L7 with its correction L9, L8, L10, L15, L16, L18, L19, L20, L21.

Tiers: DEMONSTRATED (measured, artefact on disk) / ORACLE DIAGNOSTIC (reads native-derived
quantities, never selects) / HYPOTHESIS / REFUTED / OPEN.

---

## 0. CLOSURE TABLE

| item | disposition | commit | artefacts | ledger |
|---|---|---|---|---|
| 1. test suite under the governor | DONE. 370 tests, 357 passed, 13 skipped (0 memory-guard), 0 failed, 0 errors, three jobs, peak RSS 1.692 / 0.872 / 0.324 GB | `52d15a2e` (record); tree tested `601a39c7` | `s26/TEST_RUN.md`, `s26/results/test_run.json`, `s26/results/pytest_{core,core_post,amber,amber_frame,nanpoison_post_core_edit}.xml`, `s26/jobs_done/pytest_*.json`, `s26/logs/pytest_*.log` | L16, L20 |
| 2. `verify/run_equiv2.sh` stale | CLOSED. Repaired to set `Config.project_grad`; run with `--no-amber`; four distinct cfg_keys; baseline and the shipped `exact` mode bit-identical on 8/8 targets | `eb89c165` (script + driver); `52d15a2e` (transcript) | `verify/run_equiv2.sh`, `s26/i_run_equiv2_arm.py`, `s26/i_equiv2_compare.py`, `s26/results/run_equiv2.log`, `s26/results/run_equiv2_compare.json`, `s26/jobs_done/run_equiv2.json` | L19 |
| 3. six held unused imports | CLOSED. All six applied; AST diff = exactly the removed names; tests green after | `37bddbbb` | `s26/i_ast_check.py`, `s26/jobs_done/pytest_nanpoison_post_core_edit.json`, `pytest_amber.json`, `pytest_core_post.json` | L10 |
| 4. `_archive/logs/s8/predictor_report.log` | CLOSED. The S8-13 transfer law is its only-source number; excerpt tracked with sha256 and size | `eb89c165` | `docs/sources/s8_predictor_report_excerpt.md` | L6 |
| 5. `TEST_RUN_RESULT_PLACEHOLDER` | CLOSED. Brief committed unchanged, then 6.1 replaced with the governed run | `6e50ea93`, then `83305549` | `docs/STATE_BRIEF_2026-09-12.md` section 6.1 | L20 |
| 6a. identity leak | DONE, ships dark. `identity(..., norm=)` flag, default bit-identical; in-memory audit; pinned files sha256-unchanged; the brief's criterion is at the null as a clustering threshold; the minimal fix moves exactly the 4 declared dev self-copies (+1 transitive) and 2/60 benchmark sequences | `37bddbbb` (flag + test), `601a39c7` (audit) | `core/data.py`, `tests/test_data.py`, `s26/i_identity_audit.py`, `s26/results/i_identity_audit.json`, `s26/PREREG_identity_null.md`, `s26/jobs_done/i_identity_audit3.json` | L15, L18 |
| 6b. AMBER memory guard | DONE. Skip/fail text names the 92% ceiling and quotes the governor; fixture shared with the frame-invariance file; unit test; both files green | `eb89c165` | `tests/test_amber.py`, `tests/test_amber_frame_invariance.py`, `s26/results/pytest_amber.xml`, `pytest_amber_frame.xml` | L16 |
| 6c. `pool_gate = WARN` | DONE. 1D6X and 1KWE; mechanism confirmed from the production cache; one explanatory line in the leaderboard text and the persisted rule string; takes effect on the next build, `results/` not rebuilt | `eb89c165` | `s25/resultslab/schema.py`, `s26/logs/precheck_fmt_leaderboard.log` | L7, L9 |
| 6d. moment standardisation | DONE. Production uses `core.pipeline._zrank`; every moment z-score is research-only or not on an AMBER energy; table per path | none (analysis) | L8 (the table) | L8 |
| 7. hygiene | DONE. README governor section; `python s26/examine.py` (+ `examine.sh` / `examine.bat`) = module map + pinned-hash drift + claim ledger; STATUS lines at 00:13, 00:42, 08:40, 08:4x | `eb89c165`, closing commit | `README.md`, `s26/examine.py`, `s26/i_claim_check.py`, `s26/results/claims.json`, `s26/results/claim_check.json`, `s26/logs/i_examine_full.log` | L21 |

---

## 1. THE TEST SUITE UNDER THE GOVERNOR (task 1). DEMONSTRATED.

Every job ran through `s26/jobrun.py` (registered with the governor, stdout in
`s26/logs/<job>.log`, exit / wall / peak RSS in `s26/jobs_done/<job>.json`), never more than
one AMBER job from this lane at a time, and every launch was preceded by a read of
`s26/governor_state.json` (box 64.6 to 74.6% RAM during the lane's jobs; never above 88% at a
launch; the governor never suspended or killed anything of mine, `s26/governor.log`).

| job | files | tree | tests | passed | skipped | failed/err | wall | peak RSS | in total |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| `pytest_core` | the 10 non-AMBER files | `a4db170c`, before any lane-I edit | 350 | 337 | 13 | 0 | 195.3 s | 1.694 GB | superseded |
| `pytest_amber` | `tests/test_amber.py` | `37bddbbb` | 16 | 16 | 0 | 0 | 260.5 s | 0.872 GB | yes |
| `pytest_amber_frame` | `tests/test_amber_frame_invariance.py` | `37bddbbb` | 3 | 3 | 0 | 0 | 255.7 s | 0.324 GB | yes |
| `pytest_nanpoison_post_core_edit` | `test_pipeline.py` + `test_data.py` on the edited core | `37bddbbb` | 79 | 77 | 2 | 0 | 225.5 s | 1.566 GB | no (re-run of counted files) |
| `pytest_core_post` | the 10 non-AMBER files | `601a39c7` | 351 | 338 | 13 | 0 | 240.4 s | 1.692 GB | yes |

**Suite total on the committed tree: 370 tests, 357 passed, 13 skipped, 0 failed, 0 errors,
0 memory-guard skips.** Per-file counts are in `s26/TEST_RUN.md`. The 13 skips are real and
none is from the memory guard (`-rs` reasons verbatim in `s26/TEST_RUN.md`): 3 in
`test_equivalence.py` and 8 in `test_integration.py` are `set VERIFY_SLOW=1 ...` opt-ins (the
full-pipeline / OpenMM arms); 2 in `test_pipeline.py` are absent artefacts
(`bench_results/optimised_tuning126_w6.json`; "no smoke result on disk"). The memory guard did
not fire in any job, so no re-run was needed. The suite is 370 rather than 368 because S26
added one test to `test_amber.py` (the 6b message) and one to `test_data.py` (the 6a flag).

The previous record (state brief 6.1, "355 passed, 13 skipped") was one number for one
process; the S26 record is per file, per job, with the peak RSS the governor measured, and the
same files were run before and after the production-path edits with an identical skip list.

## 2. `verify/run_equiv2.sh` (operational item 2). DEMONSTRATED, L19.

Why it was stale: it exported `PROJECT_GRAD` into the environment, which `core/project.py`
stopped reading when the gradient mode moved into `core.pipeline.Config.project_grad` so that it
reaches the cache key (`verify/grad_key_collision.py`: two modes, one `cfg_key`, the second run
served the first's arrays). The `core.pipeline` CLI has no flag for the mode, so as written the
script ran three arms in the same mode.

Repair: `s26/i_run_equiv2_arm.py` builds the identical `Config` to `python -m core.pipeline
run` (the same `replace(PROD, ...)` call) and sets `project_grad`; the script runs four arms
(baseline legacy; consolidated `exact`, the shipped default; `fd`; `analytic`), keeps the
original `cfg_key` grep, and passes any argument (here `--no-amber`) to every arm.
`s26/i_equiv2_compare.py` compares the arms per target from their cache directories with `==`.

Run: job `run_equiv2`, tag AMBER (conservative: with `--no-amber` no OpenMM context is created,
but the driver's diagnostic shows `openmm` is still imported because `core.pipeline` resolves
the amber backend at start), exit 0, 145.3 s, peak RSS 0.596 GB. Per arm on smoke8: baseline
67.4 s, exact 33.0 s, fd 20.2 s, analytic 8.8 s.

| arm | cfg_key | against `opt_exact`, 8 targets |
|---|---|---|
| `baseline_legacy` | `65ec272db3d31f05` | `ca`, `fit_ca`, `phi`, `psi`, `avg_ca`, `rmsd_avg`, `rmsd_fit`, `rmsd_arm`, `n_windows`, `n_top`: bit-identical 8/8 |
| `opt_exact` | `66050f6daae4ca07` | the shipped default |
| `opt_fd` | `85faafd84d76a827` | `avg_ca`, `rmsd_avg` identical; `rmsd_fit` differs 8/8, max 6.9e-4 A; `rmsd_arm` max 1.2e-2 A |
| `opt_analytic` | `f4e137a48586bd4b` | `avg_ca`, `rmsd_avg` identical; `rmsd_fit` max 1.0e-1 A; `rmsd_arm` max 4.8e-2 A, 8/8 |

Four distinct keys, so the collision is closed; baseline against the shipped mode bit-identical
on every array and scalar, so the consolidation is faithful; `fd` and `analytic` differ by the
builder and the gradient as `core/project.py` states. Raw coordinate max|d| between modes (22 to
32 A) and phi/psi differences near 2 pi are lab-frame and wrap differences of un-superposed
chains, not the science. Transcript: `s26/results/run_equiv2.log` (whitelisted in `.gitignore`);
comparison: `s26/results/run_equiv2_compare.json`. Not run: the AMBER-inclusive form; stage 4
is downstream of the projection the script compares.

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
row in the log (0.0629). Only the excerpt is tracked; `python s26/examine.py` checks the
transfer law as claim `s8_13_transfer_law` with the excerpt as fallback.

## 5. THE BRIEF'S SECTION 6.1 (operational item 5). CLOSED, L20.

`docs/STATE_BRIEF_2026-09-12.md` was committed exactly as received (`6e50ea93`, so history holds
the "355 passed, 13 skipped" sentence), then the live-run paragraph was replaced with the
governed run (totals, three jobs with peak RSS, passed count per file, the skip reasons, the
pointer to `s26/TEST_RUN.md` / `s26/results/test_run.json`, and the pre-edit baseline beside it)
in commit `83305549`.

## 6. THE DECLARED DEFECTS

### 6a. Identity leak. Flag ships dark; folds pinned; the brief's criterion is AT THE NULL. L15, L18.

Code (`37bddbbb`): `core.data.identity(..., norm="longer")` / `identity_many(..., norm=)`;
default statement-for-statement the pinned convention; `norm="shorter"` = match count over the
shorter sequence behind a verbatim-substring test; `clusters()` / `folds()` do not accept it;
nothing on the production path passes it. Test in `tests/test_data.py` (default bit-identical;
1CEK-in-1A11 = 1.0 under the flag; scalar and batched agree on both sides of `_BATCH_MIN`).

Audit (`s26/i_identity_audit.py`, in memory; `peptide_folds.json` / `peptide_clusters.json`
sha256 before == after == lane E's `pinned_hashes.json`; job `i_identity_audit3`, 180.6 s,
0.056 GB; `s26/results/i_identity_audit.json`):

- DEMONSTRATED: the in-memory copy of `core.data.clusters` reproduces the pinned clusters
  (470) and the pinned folds exactly.
- DEMONSTRATED: the brief's corrected criterion at 0.6 collapses 470 clusters to 166, touches
  594/787 sequences, 71/126 dev targets and 48/60 benchmark sequences (count only; how the
  benchmark sequences were obtained is L18).
- DEMONSTRATED (pre-registered, `s26/PREREG_identity_null.md`): that is chance. A real dev
  sequence passes the shorter criterion against 0.5% of the other 786 members; a
  composition-preserving shuffle passes against 0.3% (ratio 0.62; the falsifier was < 0.2). At
  mean degree ~3 on chance edges, single linkage percolates. The pinned longer criterion: 0.07%
  real vs 0.001% shuffled. So the flag is a per-pair leak test, not a clustering threshold, and
  its docstring says so.
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

## 7. HYGIENE. DEMONSTRATED, L21.

README section "The S26 resource governor" (band 88/90/93/95, AMBER cap 2, tags, log and state
paths, how to start it, the split test run). `python s26/examine.py` = `make examine` (no
Makefile, no `make` on this box; `examine.sh` / `examine.bat` at the root, no conflict): lane
E's `e_module_map` (module map), lane E's `e_hashes --check` (pinned-hash drift), lane I's
`i_claim_check` (claim ledger), `--search` for lane E's `e_claims`. Full run at 08:41 (job
`i_examine_full`, 15 s, 0.299 GB, `s26/logs/i_examine_full.log`): 725 modules mapped, backend
reports default and legacy as in the README, **no drift across 101 pinned entries** (the
benchmark manifest matches the S20 record), **21/21 claims OK** (`s26/results/claim_check.json`).
`s26/i_ast_check.py` makes the contract's AST-identity rule runnable; `s26/i_test_report.py`
turns junit XML plus `jobs_done` into the test record. STATUS lines at 00:13, 00:42, 08:40, 08:4x.

## 8. WHAT I DID NOT DO, AND WHY

- Did not re-derive folds or clusters, did not call `core.data.clusters()` / `folds()` or their
  `peptide_db` twins from the audit, did not write `peptide_folds.json` / `peptide_clusters.json`
  (hashes verified before and after, and again by `examine.py` at 08:41). The fix ships dark.
- Did not add the `norm` flag to `peptide_db.identity`: the legacy arm must stay pure, and
  `tests/test_data.py` pins `core.data.identity == peptide_db.identity` on the default path.
- Did not rebuild `results/`; the 6c text change is dormant until the next frozen build.
- Did not run `verify/run_equiv2.sh` with AMBER: `--no-amber` leaves the projection comparison
  intact; the script supports both and says so.
- Did not open `results/benchmark_manifest.json` directly, any benchmark PDB, native or RMSD.
  The benchmark counts came from `core.backend("data").benchmark()` sequences only, which is the
  route the brief named and which itself parses the manifest for its ids (stated plainly in
  L18); no benchmark name was printed; non-dev database members appear by length only. Nothing
  further will be computed on benchmark targets by this lane.
- Did not add a claim for the benchmark headline (+0.0103, CI [-0.160, +0.180]) to
  `claims.json`: I did not find its results artefact and would not invent a path.
- Did not add a Makefile: no `make` on this box; the two launchers and `s26/examine.py` are
  the target.

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
   in its own fold). Only the per-target verbatim breakdown in the third audit pass separated
   direct from transitive.
5. **The frame-invariance file ran in 255.7 s at 0.324 GB peak.** S25 L14 recorded it as the
   test most likely to hit the ceiling; on a 65% box under the governor's budget it is not close.
6. **`--no-amber` still imports openmm.** I expected the flag to keep OpenMM out of the process;
   `core.pipeline` resolves the amber backend at start, so only the contexts are avoided. The
   AMBER tag on that job was conservative and correct.

## 10. PRODUCTION-PATH CHANGES MADE (all deliberate, tested, in the ledger)

| file | change | flag / default | test |
|---|---|---|---|
| `core/data.py` | `identity(..., norm=)`, `identity_many(..., norm=)`; `_seq_index` removed | `norm="longer"` = old behaviour, bit-identical | `tests/test_data.py` (new test + the four existing identity pins) |
| `core/amber.py` | two unused names dropped from the `budget` import | n/a | `tests/test_amber.py` 16/16 |
| `core/bench.py`, `core/cache.py`, `core/predict.py` | one unused typing/dataclasses import each | n/a | non-AMBER suite, `pytest_core_post` |
| `verify/vqe_lfo_audit.py` | unused `defaultdict` | n/a (not production) | byte-compiled |
| `s25/resultslab/schema.py` | WARN explanation in `fmt_leaderboard` and `pool_gate_rule` | text only | rendered on the on-disk payload |
| `tests/test_amber.py`, `tests/test_amber_frame_invariance.py`, `tests/test_data.py` | 6b message + shared fixture; 6a flag test | n/a | themselves |

No number changed anywhere: the identity default is asserted bit-identical, the removed names
were never read, the leaderboard change is text, and the four-arm equivalence run shows the
shipped mode bit-identical to the legacy baseline after the edits.

## 11. ARTEFACT INDEX

`s26/TEST_RUN.md`, `s26/results/test_run.json`, `s26/results/pytest_{core,core_post,amber,amber_frame,nanpoison_post_core_edit}.xml`,
`s26/jobs_done/*.json`, `s26/logs/*.log`; `docs/sources/s8_predictor_report_excerpt.md`;
`verify/run_equiv2.sh`, `s26/i_run_equiv2_arm.py`, `s26/i_equiv2_compare.py`,
`s26/results/run_equiv2.log`, `s26/results/run_equiv2_compare.json`; `s26/i_identity_audit.py`,
`s26/results/i_identity_audit.json`, `s26/PREREG_identity_null.md`; `s26/i_ast_check.py`;
`s26/i_test_report.py`; `s26/examine.py`, `s26/i_claim_check.py`, `s26/results/claims.json`,
`s26/results/claim_check.json`, `examine.sh`, `examine.bat`; README section; ledger L6, L7, L8,
L9, L10, L15, L16, L18, L19, L20, L21.

---

# EXTENDED WINDOW (L77, 22:05 to 04:30): EVERY TEST AND DIAGNOSTIC THAT COULD BE RUN

Closure table for the six extended-window items (commit hashes are lane-I commits unless noted).

| item | disposition | commit | artefacts | ledger |
|---|---|---|---|---|
| (1) opt-in test tier, `VERIFY_SLOW=1` on the 11 skips | DONE: 11 run, 11 passed; equivalence 3/3 (40.2 s, 0.313 GB), integration 8/8 (165.6 s, 1.139 GB). Two defects found on the way: the L92 relaunch had dropped the flag from my shell (one null run set aside, L98), and `tests/test_equivalence.py`'s summary parser could never have run its three arms (fixed, L102). Suite with the tier folded in: **370 tests, 368 passed, 0 failed, 2 skipped** (absent artefacts) | `7be8e4b0`, `7e08b968`, `7ad4ef68` (+ `a6ce3ab6`, swept in by PH) | `s26/TEST_RUN.md`, `s26/results/test_run.json`, `s26/results/pytest_slow_{equivalence,integration}.xml`, `s26/logs/pytest_slow_*`, the `.NULLRUN_*` / `.PARSER_SKIP_*` files set aside | L98, L102, L104, L113 |
| (2) every `verify/` audit re-run, JSON beside the tracked one, diffed | VERIFY_DISPOSITION. Nothing tracked was overwritten (sha256 asserted per audit; deletion guards probed) | `b0912487`, `2505e575`, closing commit | `s26/i_verify_rerun.py`, `s26/results/verify/REPORT.md`, `REPORT.json`, `<audit>.rerun.json`, `verify__<name>.json`, `s26/logs/verify_*.log`, `s26/jobs_done/verify_*.json` | VERIFY_TABLE_LEDGER |
| (3) frozen results-lab rebuild | DONE, **REPRODUCED**: 2016/2016 per-target RMSDs exact on both bases, every leaderboard number / gate / verdict identical (3.2126 / 3.0483, WARN 2, PASS), 2142/2142 PDB ATOM records identical, L7 WARN text present; `results/summary` committed, structures restored; the tracked leaderboard's three descriptive columns were written by uncommitted code and do not reproduce | `083c9b95` | `s26/i_resultslab_rebuild.py`, `s26/results/resultslab_rebuild/{snapshot_before.json,post_verdict.json,post_stdout.txt,summary_before/}`, `s26/logs/resultslab_rebuild.log` | L127 |
| (4) final AST gate against `ae86a124` | DONE: 55 production modules, 50 identical, the 5 differing files exactly the L10 / L15 edits | `858a1be2` | `s26/results/ast_gate_ae86a124.txt` | L98 |
| (5) `python s26/examine.py` on the final tree | DONE (job `i_examine_final`, 15 s, peak RSS 0.452 GB): 787 modules mapped; no drift across the 101 pinned entries, the sealed benchmark manifest still hashes to the S20 record; 21/21 claims re-read OK from their artefacts (production 3.2148 / 3.0483 / 3.2355 / 3.2041 means of the 126 cache records, the compare deltas, the leaderboard row and gate, the pinned hashes, the S8-13 transfer law, the S26 test totals) | closing commit | `s26/logs/i_examine_final.log`, `s26/results/claim_check.json`, `module_map.json` | closing entry |
| (6) ledger entries and this table | DONE | closing commit | this section | L98, L102, L104, L113, L127, L128, VERIFY_TABLE_LEDGER, closing entry |

## E1. The opt-in tier (item 1)

Tiers: DEMONSTRATED. Job `pytest_slow_equivalence3` (AMBER, flag inside the command, tree `7e08b968`): 3 passed in 40.2 s, peak RSS 0.313 GB; both arms served from the on-disk caches (production `1fc9f2dcf489e2fb`, baseline `464a0ddb5f283e04`, mtimes 2026-09-04), which the test's docstring allows, so the tier certifies the equivalence of the two recorded arms up to the projection (`avg_ca` and seven scalars `==` on 8/8) and that the projection's divergence stays inside its pinned band. Job `pytest_slow_integration2` (AMBER, registered 01:33 after 4,620 s in the queue): 8 passed in 165.6 s, peak RSS 1.139 GB: Legacy terms column-for-column on a real pool; ff14SB/GBn2 parameters against a System built from the same XML; the pinned 1A13 interaction energy bit-exact; rigid-translation invariance; NaN-poisoning through `run_target` with stage 4 on 1CS9 and 1CB3; cross-process and 4-thread bit-identity through child processes.

What damaged my expectations here: (a) the relaunch by another lane dropped `VERIFY_SLOW=1` because it lived in my shell, not in the command (L98); an environment a job needs belongs in its command. (b) The three equivalence arms had never run, in any recorded suite: `_run_arm` located the summary with `rfind("{")`, which on an `indent=2` dump is a nested brace, and the `except` turned the parse failure into a skip (L102). Every "13 skipped" on the record included three tests that could not run.

## E2. The `verify/` audits (item 2)

Runner: `s26/i_verify_rerun.py` executes each audit in-process with `builtins.open`, `os.replace`, `os.rename`, `os.remove`, `os.unlink`, `shutil.move`, `shutil.rmtree` wrapped: writes under `verify/` or `bench_results/` are redirected to `s26/results/verify/`, deletions and move-outs there are refused (six operations probed), the tracked file's sha256 is asserted unchanged after every run, and the fresh JSON is compared with the tracked one leaf by leaf. `determinism_audit` runs D1/D3/D4/D5 by function and not D2, which as written moves the last two PRODUCTION cache records to a temp dir, re-runs smoke8 (which does not contain them) and rmtree's the backup.

VERIFY_TABLE_FINDINGS

## E3. The rebuild (item 3)

DEMONSTRATED, L127. Governed job `resultslab_rebuild` (the documented command unchanged; 4,472.9 s under seven concurrent jobs, peak RSS 0.11 GB, after 4,765 s in the queue). Bracketed by `s26/i_resultslab_rebuild.py`: `pre` snapshotted the tracked `results/summary` and the sha256 + ATOM-record sha256 of all 2,142 tracked PDBs; `post` compared. Verdict REPRODUCED on every number (`s26/results/resultslab_rebuild/post_verdict.json`). Byte-wise, the committed `results/summary` changed in: provenance blocks; per-record timestamp / git_commit / module_hash; `results.csv`'s `hamiltonians` (code -> name) and `distogram_used` ('' -> true/false), now consistent with `results.json`; `pool_gate_rule` now carries the L7 text. `leaderboard.json` lost three descriptive per-row keys (`selector`, `hamiltonians`, `distogram_used`) that no committed `schema.py` produces: the tracked build ran from a dirty working copy. Not hand-patched. `results/structures` restored to the tracked bytes (only `REMARK GIT_COMMIT` / `MODULE_SHA` differed).

## E4. The AST gate (item 4)

DEMONSTRATED, L98. `s26/i_ast_check.py --ref ae86a124` over `core/` (11), the 24 root modules and `s5/ s7/ s8/ s9/` (20): 50 AST-identical modulo docstrings; `core/amber.py`, `core/bench.py`, `core/cache.py`, `core/predict.py` differ by one removed import each (L10) and `core/data.py` by the `norm` keyword on `identity` / `identity_many` and the removed dead `_seq_index` (L15). `s26/results/ast_gate_ae86a124.txt`.

## E5. Governor defect found on the way (L128)

A job whose root process blocks in `subprocess.run` reads "running" to psutil after `suspend()`, so the governor re-suspended my `verify_grad_key_collision` tree on every hot sample (18 SUSPEND lines, 0 RESUME) and, Windows counting suspends per thread, the child needed 16 `resume()` calls; 28 minutes lost. `s26/i_tree_watchdog.py` guards lane I's own jobs only (resumes a stopped descendant when the root is running; never touches a governor-suspended root or another lane's job). Two-line fix proposed in L128 for the governor itself.

## E6. Hygiene list for the next sprint (from what broke in this one)

1. **A launcher's cap is read at every tick from the day it is written.** Seven waiters started before jobrun v2.3 read the cap of four at import and starved for two hours under a file cap of six (L92); the fix was one line, the cost was a lane's evening.
2. **A file-based stop request instead of CTRL_BREAK.** The governor's kill path sends `CTRL_BREAK_EVENT` to a process group and waits 25 s for a checkpoint; on Windows the signal does not reach a grandchild reliably and a suspended tree cannot handle it at all. A job that polls `s26/stop/<name>` every few seconds and checkpoints on its own is deterministic and testable.
3. **The supervisor's own state is the truth for suspensions.** The governor decided "suspended" from psutil's status of the ROOT process; a root blocked in `subprocess.run` reads running after `suspend()`, so it was re-suspended 18 times, never resumed, and the child accumulated 16 stacked suspends (L128). Keep the suspended set in the governor's state, act on that, and on resume loop until the tree runs.
4. **A test's summary parser must fail, not skip.** `tests/test_equivalence.py` turned a parse failure into a skip, and the three opt-in arms had never run in any recorded suite (L102). A skip is a statement that the test does not apply; a parse failure is a defect.
5. **A job's environment belongs in its command.** `VERIFY_SLOW=1` lived in my shell; the relaunch by another lane re-ran the command faithfully and the tier ran as eight skips in five seconds (L98).
6. **Stage and commit in one call.** Eight lanes share one git index; four staged lane-I files were swept into another lane's commit (L104). Or give each lane a worktree.
7. **A number enters the ledger only in a turn after its output was read.** L7's two wrong point-cloud numbers, corrected in L9.
8. **A generated artefact is not hand-patched; a column no committed code produces is a finding.** The tracked leaderboard's three descriptive columns came from a dirty working copy (L127); the repair is the one-function schema fix and a rebuild, not an edit.
9. **No audit under `verify/` may mutate a production cache.** `determinism_audit` D2 moves and then deletes two records of `bench_results/cache/1fc9f2dcf489e2fb` as written (E2); re-runs must redirect writes and refuse deletions under evidence directories (`s26/i_verify_rerun.py` does both).
10. **An audit whose arms are cache keys must record the keys it needs and refuse silently vacuous comparisons.** `project_arms` with a repeated key reports identity by construction; `equiv_compare` recomputes keys that move whenever `Config` gains a field.
11. **Estimate memory from a measured peak, never from a guess, and record it beside the result.** Every job in this lane carries its `peak_rss_gb` from `s26/jobs_done/`; the largest lane-I job of the window was the integration tier at 1.139 GB against an estimate of 1.8.

## E7. What I did not do in the window, and why

- Did not run the audits past the 03:30 deadline; the heavy 126-target projection audits that did not complete are recorded as NOT RUN with their cost class.
- Did not hand-patch the rebuilt leaderboard's three missing descriptive columns: a generated artefact is not edited by hand; the one-function schema fix and a 75-minute rebuild are recommended instead (L127).
- Did not run `determinism_audit` D2 (it would delete two production-cache records) nor `projection_divergence` on invented arms (its two smoke8 caches are gone; on baseline vs production it trips on 7 S11 baseline records whose AMBER declined at the 92% ceiling; `project_arms` fresh answers its question: 126/126 bit-identical).
- Did not run the AMBER-inclusive `verify/run_equiv2.sh` (L19 stands).
