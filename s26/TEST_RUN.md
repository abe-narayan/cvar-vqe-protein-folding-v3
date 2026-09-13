# S26 TEST RUN (under the governor)

Every job below ran through `s26/jobrun.py` (registered with the governor, stdout in `s26/logs/<job>.log`, exit code / wall / PEAK RSS in `s26/jobs_done/<job>.json`). Counts come from the junit XML each job wrote to `s26/results/<job>.xml`; the machine-readable form is `s26/results/test_run.json`. A skip whose message names the memory ceiling is a memory-guard skip and is counted apart from real skips.

## job `pytest_core`

- command: `C:/Users/abena/miniforge_3/python.exe -m pytest tests/ -q -rs -p no:cacheprovider --deselect tests/test_amber.py --deselect tests/test_amber_frame_invariance.py --ignore=tests/test_amber.py --ignore=tests/test_amber_frame_invariance.py --junitxml=s26/results/pytest_core.xml`
- commit: `a4db170c`  note: the non-AMBER files, run BEFORE any S26 lane-I edit (tree = a4db170c plus the coordinator's untracked s26/ files); 13 skips are all VERIFY_SLOW or absent-artefact skips, none from the memory guard
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

## Combined (jobs marked as counting toward the suite total)

**350 tests: 337 passed, 0 failed, 0 errors, 13 skipped (0 memory-guard skips).**

Rendered 2026-09-13 00:29 by `s26/i_test_report.py`.
