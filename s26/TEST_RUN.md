# S26 TEST RUN (under the governor)

Every job below ran through `s26/jobrun.py` (registered with the governor, stdout in `s26/logs/<job>.log`, exit code / wall / PEAK RSS in `s26/jobs_done/<job>.json`). Counts come from the junit XML each job wrote to `s26/results/<job>.xml`; the machine-readable form is `s26/results/test_run.json`. A skip whose message names the memory ceiling is a memory-guard skip and is counted apart from real skips.

## job `pytest_core`

- command: `C:/Users/abena/miniforge_3/python.exe -m pytest tests/ -q -rs -p no:cacheprovider --deselect tests/test_amber.py --deselect tests/test_amber_frame_invariance.py --ignore=tests/test_amber.py --ignore=tests/test_amber_frame_invariance.py --junitxml=s26/results/pytest_core.xml`
- commit: `a4db170c`  note: the non-AMBER files BEFORE any lane-I edit (tree a4db170c); superseded in the suite total by pytest_core_post, which re-ran the same files on the committed tree; kept as the pre-edit baseline (337 passed / 13 skipped / 0 failed, identical skip list)
- start 2026-09-13T00:11:10  end 2026-09-13T00:14:25  wall 195.3 s  exit 0  peak RSS 1.694 GB  tag TEST  est 2.5 GB
- totals: 350 tests, 337 passed, 0 failed, 0 errors, 13 skipped (0 memory-guard)

| file | tests | passed | failed | errors | skipped | memory-guard skips |
|---|---:|---:|---:|---:|---:|---:|
| `test_cvar.py` | 90 | 90 | 0 | 0 | 0 | 0 |
| `test_data.py` | 41 | 41 | 0 | 0 | 0 | 0 |
| `test_energy.py` | 9 | 9 | 0 | 0 | 0 | 0 |
| `test_equivalence.py` | 14 | 11 | 0 | 0 | 3 | 0 |
| `test_geometry.py` | 24 | 24 | 0 | 0 | 0 | 0 |
| `test_instrument.py` | 43 | 43 | 0 | 0 | 0 | 0 |
| `test_integration.py` | 25 | 17 | 0 | 0 | 8 | 0 |
| `test_pipeline.py` | 37 | 35 | 0 | 0 | 2 | 0 |
| `test_project.py` | 28 | 28 | 0 | 0 | 0 | 0 |
| `test_quantum.py` | 39 | 39 | 0 | 0 | 0 | 0 |

skip reasons (`-rs`):

- `test_equivalence.py::test_everything_up_to_the_projection_is_bit_identical` [real]: set VERIFY_SLOW=1 to run the full pipeline arms
- `test_equivalence.py::test_the_projection_selects_a_different_degenerate_branch` [real]: set VERIFY_SLOW=1 to run the full pipeline arms
- `test_equivalence.py::test_no_stage_is_skipped_in_the_consolidated_arm` [real]: set VERIFY_SLOW=1 to run the full pipeline arms
- `test_integration.py::test_legacy_terms_match_the_shipped_module_on_a_real_pool` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_amber_is_real_ff14sb_gbn2_with_the_shipped_parameters` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_the_1a13_native_interaction_energy_is_bit_exact` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_single_point_energy_is_invariant_under_rigid_translation` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_nan_poisoning_the_native_changes_no_deployable_quantity[1CS9]` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_nan_poisoning_the_native_changes_no_deployable_quantity[1CB3]` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_the_same_target_is_bit_identical_across_processes` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_the_result_does_not_depend_on_the_ambient_thread_count` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_pipeline.py::test_tuning_instrument_reference_numbers[optimised_tuning126_w6]` [real]: bench_results/optimised_tuning126_w6.json not present -- run the harness first
- `test_pipeline.py::test_the_harness_refuses_to_call_a_development_run_a_headline` [real]: no smoke result on disk

## job `pytest_amber`

- command: `C:/Users/abena/miniforge_3/python.exe -m pytest tests/test_amber.py -q -rs -p no:cacheprovider --junitxml=s26/results/pytest_amber.xml`
- commit: `37bddbbb`  note: tests/test_amber.py with the S26 memory-guard message (defect 6b) and its new unit test (16 = 15 + 1); core/amber.py already carried the import removal committed as 37bddbbb; the guard did not fire
- start 2026-09-13T00:29:14  end 2026-09-13T00:33:34  wall 260.5 s  exit 0  peak RSS 0.872 GB  tag AMBER  est 2.0 GB
- totals: 16 tests, 16 passed, 0 failed, 0 errors, 0 skipped (0 memory-guard)

| file | tests | passed | failed | errors | skipped | memory-guard skips |
|---|---:|---:|---:|---:|---:|---:|
| `test_amber.py` | 16 | 16 | 0 | 0 | 0 | 0 |

## job `pytest_nanpoison_post_core_edit`

- command: `C:/Users/abena/miniforge_3/python.exe -m pytest tests/test_pipeline.py tests/test_data.py -q -rs -p no:cacheprovider --junitxml=s26/results/pytest_nanpoison_post_core_edit.xml`
- commit: `37bddbbb`  note: tests/test_pipeline.py (the NaN-poison test) + tests/test_data.py re-run on the edited core (identity flag, import removals, dead helper); these files are already counted through pytest_core, so this job is recorded but not added to the suite total
- start 2026-09-13T00:29:58  end 2026-09-13T00:33:44  wall 225.5 s  exit 0  peak RSS 1.566 GB  tag TEST  est 1.7 GB
- totals: 79 tests, 77 passed, 0 failed, 0 errors, 2 skipped (0 memory-guard)

| file | tests | passed | failed | errors | skipped | memory-guard skips |
|---|---:|---:|---:|---:|---:|---:|
| `test_data.py` | 42 | 42 | 0 | 0 | 0 | 0 |
| `test_pipeline.py` | 37 | 35 | 0 | 0 | 2 | 0 |

skip reasons (`-rs`):

- `test_pipeline.py::test_tuning_instrument_reference_numbers[optimised_tuning126_w6]` [real]: bench_results/optimised_tuning126_w6.json not present -- run the harness first
- `test_pipeline.py::test_the_harness_refuses_to_call_a_development_run_a_headline` [real]: no smoke result on disk

## job `pytest_amber_frame`

- command: `C:/Users/abena/miniforge_3/python.exe -m pytest tests/test_amber_frame_invariance.py -q -rs -p no:cacheprovider --junitxml=s26/results/pytest_amber_frame.xml`
- commit: `37bddbbb`  note: tests/test_amber_frame_invariance.py, now carrying tests/test_amber.py's autouse memory guard by import (S25 L14 approval, S26 defect 6b); the guard did not fire; 3 = the file's own count
- start 2026-09-13T00:35:34  end 2026-09-13T00:39:50  wall 255.7 s  exit 0  peak RSS 0.324 GB  tag AMBER  est 1.5 GB
- totals: 3 tests, 3 passed, 0 failed, 0 errors, 0 skipped (0 memory-guard)

| file | tests | passed | failed | errors | skipped | memory-guard skips |
|---|---:|---:|---:|---:|---:|---:|
| `test_amber_frame_invariance.py` | 3 | 3 | 0 | 0 | 0 | 0 |

## job `pytest_core_post`

- command: `C:/Users/abena/miniforge_3/python.exe -m pytest tests/ -q -rs -p no:cacheprovider --deselect tests/test_amber.py --deselect tests/test_amber_frame_invariance.py --ignore=tests/test_amber.py --ignore=tests/test_amber_frame_invariance.py --junitxml=s26/results/pytest_core_post.xml`
- commit: `601a39c7`  note: the 10 non-AMBER files on the committed tree 601a39c7 (identity flag + its test, import removals, dead helper removed); this is the suite-total run; the 13 skips are the same VERIFY_SLOW / absent-artefact skips as before the edits, none from the memory guard
- start 2026-09-13T00:41:56  end 2026-09-13T00:45:56  wall 240.4 s  exit 0  peak RSS 1.692 GB  tag TEST  est 2.0 GB
- totals: 351 tests, 338 passed, 0 failed, 0 errors, 13 skipped (0 memory-guard)

| file | tests | passed | failed | errors | skipped | memory-guard skips |
|---|---:|---:|---:|---:|---:|---:|
| `test_cvar.py` | 90 | 90 | 0 | 0 | 0 | 0 |
| `test_data.py` | 42 | 42 | 0 | 0 | 0 | 0 |
| `test_energy.py` | 9 | 9 | 0 | 0 | 0 | 0 |
| `test_equivalence.py` | 14 | 11 | 0 | 0 | 3 | 0 |
| `test_geometry.py` | 24 | 24 | 0 | 0 | 0 | 0 |
| `test_instrument.py` | 43 | 43 | 0 | 0 | 0 | 0 |
| `test_integration.py` | 25 | 17 | 0 | 0 | 8 | 0 |
| `test_pipeline.py` | 37 | 35 | 0 | 0 | 2 | 0 |
| `test_project.py` | 28 | 28 | 0 | 0 | 0 | 0 |
| `test_quantum.py` | 39 | 39 | 0 | 0 | 0 | 0 |

skip reasons (`-rs`):

- `test_equivalence.py::test_everything_up_to_the_projection_is_bit_identical` [real]: set VERIFY_SLOW=1 to run the full pipeline arms
- `test_equivalence.py::test_the_projection_selects_a_different_degenerate_branch` [real]: set VERIFY_SLOW=1 to run the full pipeline arms
- `test_equivalence.py::test_no_stage_is_skipped_in_the_consolidated_arm` [real]: set VERIFY_SLOW=1 to run the full pipeline arms
- `test_integration.py::test_legacy_terms_match_the_shipped_module_on_a_real_pool` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_amber_is_real_ff14sb_gbn2_with_the_shipped_parameters` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_the_1a13_native_interaction_energy_is_bit_exact` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_single_point_energy_is_invariant_under_rigid_translation` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_nan_poisoning_the_native_changes_no_deployable_quantity[1CS9]` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_nan_poisoning_the_native_changes_no_deployable_quantity[1CB3]` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_the_same_target_is_bit_identical_across_processes` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_integration.py::test_the_result_does_not_depend_on_the_ambient_thread_count` [real]: set VERIFY_SLOW=1 for the OpenMM / pipeline checks
- `test_pipeline.py::test_tuning_instrument_reference_numbers[optimised_tuning126_w6]` [real]: bench_results/optimised_tuning126_w6.json not present -- run the harness first
- `test_pipeline.py::test_the_harness_refuses_to_call_a_development_run_a_headline` [real]: no smoke result on disk

## Combined (jobs marked as counting toward the suite total)

**370 tests: 357 passed, 0 failed, 0 errors, 13 skipped (0 memory-guard skips).**

Rendered 2026-09-13 08:38 by `s26/i_test_report.py`.
