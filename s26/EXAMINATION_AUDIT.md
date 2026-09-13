# S26 EXAMINATION AUDIT (lane A, Adversary)

Written 2026-09-13, 08:55 to 09:15, on branch `s26`, HEAD `76287a33` at audit time (the working
tree differs from HEAD only in shared ledger/status files and other lanes' untracked files; no
production module differs). Object of the audit: `s26/EXAMINATION.md`, `s26/BRIEF.md`,
`s26/agentE_FINDINGS.md`, all at commit `80159ae9` (2026-09-13 08:52:46), file mtimes 08:52:09,
08:51:50 and 08:52:10, unchanged from the start of the audit to the writing of this file
(re-checked with `git log -1 -- s26/EXAMINATION.md` and `stat` before writing). Ledger read L0
to L29. Every number below carries the artefact it was read from. Verdict at the end.

Labels: MATERIAL blocks Phase 1 until fixed and re-checked; MINOR is recorded and does not block.

## 1. The three reproductions (Part 2.3)

What I ran. `s26/e_reproduce.py` hard-codes its output path to lane E's committed artefact
`s26/results/e_reproduce.json`, so re-running it as written would overwrite another lane's file
(contract section 0). `s26/a_reproduce_head.py` executes `e_reproduce.main()` unchanged and
redirects its single `ST.save_atomic` call to `s26/results/a_reproduce_head.json`; the
computation is the Examiner's code byte for byte (the output's provenance `source_sha256`
4bd3111ea1007fb7 is the sha of `s26/e_reproduce.py`, equal at HEAD and at `e3570aed`). Job
`a_reproduce_head` through `s26/jobrun.py` (`s26/jobs_done/a_reproduce_head.json`: exit 0, wall
5.0 s, sampler peak RSS 0.005 GB, an under-read of a 5 s process, the same caveat E records; log
`s26/logs/a_reproduce_head.log`). Lane E's file was not touched: sha256 8b8eb38aceb6fc25...
before and after.

Provenance: stored file `git_commit` 74e073e26d5c dirty (E's run, 07:12 UTC); fresh file
`git_commit` 76287a339c6b dirty (HEAD, 15:56 UTC). HEAD carries lane I's `37bddbbb`
production-path edits, which E's section I.1 says its numbers predate.

| quantity | E quotes (EXAMINATION H) | `e_reproduce.json :: summary/*/mean_fresh` | `a_reproduce_head.json :: summary/*/mean_fresh` | max per-target abs disagreement fresh vs stored, both files |
|---|---|---|---|---|
| rmsd_avg | 3.048338093879532 | same | same | 0.0 |
| rmsd_fit | 3.2040761603809194 | same | same | 0.0 |
| rmsd_arm | 3.2147651542109985 | same | same | 0.0 |
| rmsd_full | 3.2354598538973844 | same | same | 0.0 |

T030 is 1S9Z (`results/summary/target_map.json :: t_of_pdb/1S9Z` = "T030").
`summary/T030_1S9Z/rmsd_arm` 0.18198112330908295 in both files, equal to `stored_rmsd_arm`; avg
0.19588790009124518, fit 0.1819434462106372, full 0.3173791640496416, as quoted. Across the 126
rows and seven columns (`rmsd_avg`, `rmsd_fit`, `rmsd_arm`, `rmsd_full`, `bond_mean_avg`,
`bond_mean_chain`, `bond_mean_native`) the max abs difference between the stored and the HEAD
file is 0.0. `fold_mismatch_vs_pinned_folds` is [] in both; `cfg_keys` ['1fc9f2dcf489e2fb'];
`amber_errors` []. The virtual-bond block is identical (avg 2.9613938229671057, native
3.812176007269521, contraction 22.317%, chain 3.803954938363982, chain sd max 1.8e-15). The
stored aggregates E quotes (shipped 3.4540004952559396, pool_best 1.7108244199364904, top_m_best
2.3061526409453816) are `summary/stored_shipped/mean`, `stored_pool_best/mean`,
`stored_top_m_best/mean`.

The 2.6e-04 A tolerance is met with zero disagreement on HEAD. E's caveat I.1 is closed by this
re-run. No finding.

## 2. The claim ledger

Method: a script of mine read each cited key path from each cited JSON file (nothing under
`s9/final_*` and not `results/benchmark_manifest.json`). Sample: `random.Random(26).sample(...)`
over the sourced ids excluding C06, C24, C26, C27, C34, C35 gives C02, C05, C08, C15, C17, C19,
C21, C23, C28, C32. I opened those (C05 and C28 through the test status in the junit files,
section 4), every UNSOURCED entry, every markdown-sourced entry, and 17 others: 29 artefact
leaves in total.

Sourced entries, value found at the cited key, stored precision; every one equals E's quote:

| id | artefact :: key | value found |
|---|---|---|
| C01 | `bench_results/baseline_tuning126.json :: science/rmsd_arm/mean`; `optimised_tuning126_w8.json :: science/rmsd_arm/mean` | 3.214765154210998; 3.2147651542109985 |
| C02 | `baseline_tuning126.json :: science/rmsd_avg/mean` | 3.048338093879532 |
| C03 | `baseline_tuning126.json :: science/rmsd_full/mean`; `bench_results/cache/1fc9f2dcf489e2fb/9KAR.json :: amber_e1`; `2BP4.json :: amber_e1`, `amber_strain_after` | 3.2354598538973844; 1262.4103583126937; 845.3467373422843, 1172.6967039637698 |
| C07 | `baseline_tuning126.json :: science/pool_best/mean`; `results/summary/pool_best.json :: mean` | 1.7108244199364904 both |
| C08 | `s25/results/phys_suite.json :: random_null_rank/mean`, `random_null_moment/mean` | 3.4251256499338507 both |
| C09 | `s12/results/s14_ladder.json :: rows/L0_constant_helix/mean`; `s14_coherence.json :: real_emitters/L0_constant_helix/rmsd` | 4.064753929494389 both |
| C10 | `s13/results/tors_eval.json :: reports/a_pepPos/build/mean` | 3.7704579282796167 |
| C12 | `bench_results/recon_library_saturation.json :: universe_best/small` | 1.3134468768690186 |
| C13 | `s13/results/ceiling_report.json :: k/4/descent/mean`; `adv_ceiling_k4_report.json :: arms/real/summary/mean` | 1.5940389994797128 both |
| C14 | `s15/results/distgeo.json :: arms/ORACLE_true_distances/mean` | 0.6108221854859448 |
| C18 | `s25/results/q_gibbs.json :: results/divergence_ladder/1/ladder/"50 Adam steps (DEPLOYED)"/kl_nats` (T 0.3, `deployed_at_this_T` true) | 0.9019158325027764 |
| C19 | `s25/results/q_alpha.json :: results/expressivity_vs_boltzmann/"circuit - exact Boltzmann @T=0.3"` mean, mde, eff_over_mde | -0.03024176839325163, 0.12605272002164683, 0.2399136519073786 (W 43, L 44, ties 39, verdict NULL) |
| C20 | `q_gibbs.json :: results/spectrum_target_independence` worst, E_range, worst_frac_of_range | 0.04063087303875024, 3.4371432161567106, 0.011821117271971639 (n_checked 8) |
| C21 | `s25/results/q_verify.json :: results/"max \|p_statevector - p_dense\|"` | "5.551e-17" |
| C22 | `q_verify.json :: results/"cos(param-shift, exact FD), min over 12"`, `"cos(free-energy grad, FD), min over 4 (a,T)"` | "1.000000000" both |
| C23 | `q_verify.json :: results/"rel \|g_ps - g_fd\| / \|g_fd\|, max over 12"`, `"rel err of free-energy grad, max over 4"` | "4.597e-10", "4.663e-10" |
| C24 | `q_verify.json :: results/"cos(exact-expectation, TAIL baseline)"`, `"\|g_tail\| / \|g_exact\|"`; `core/quantum.py:42`; `tests/test_quantum.py:59` | "0.566586", "0.519"; +0.655634 at 0.758; `REC_TAIL_COS = 0.655634` (asserted at `:327`, passed) |
| C25 | `s25/results/q_plateau.json :: results/{linear_alpha1_T0, cvar_alpha025_T0, cvar_alpha01_T0, deployed_a1_T03}/log2_slope_per_qubit` | -0.6491874127317976, -0.2522012761595383, -0.047158929969141276, -0.31053500868274697 |
| C29 | `phys_suite.json :: configs_rank/{3,5,4,7,6,1,2}/mean` | 3.058016219425174, 3.131686249878013, 3.215116303141773, 3.2530157270854128, 3.6742023446060776, 3.7553415439547573, 3.880496706249537 |
| C30 | `phys_suite.json :: controls_rank/{1,2}/vs_random/{effect, mde, effect_over_mde, ci95_fold, folds_same_sign}` | 0.3302158940209069, 0.1662403218173309, 1.9863766528541527, [0.23753999909000342, 0.4158172931555582], 5; 0.45537105631568636, 0.18834102878507575, 2.4178006207841753, [0.328398404491368, 0.5355771582778813], 5 |
| C31 | `q_alpha.json :: results/cell_correlations/H` | -0.7422836835050723 |
| C32 | `q_verify.json :: results/n_cells`, `"SUBSET-hood violations  <- THE THEOREM"`, `"full-support cells with EXACT value equality"` | 2592, 0, "972 / 972" |
| C33 | `q_gibbs.json :: results/training_control/{0,1,2}/frac_gap_closed` | 0.8914955163069174, 0.782951322968052, 0.7833283350457241 |
| C04, C05, C28 | `s26/results/pytest_core.xml`, `pytest_core_post.xml` | `test_the_committed_benchmark_report_still_holds_its_reference_numbers` and `test_benchmark60_cached_record_reproduces_and_is_not_rerun` passed in both; `tests/test_data.py:74` `ok_canon >= 4` inside a passed test |

Derived entries reproduced: C11 mean of `s24/results/priorladder.json :: rows[*]/MASS1.0` =
2.226080063107832 (TILT1.0 and DIRAC1.0 identical, MASSFIXW1.0 2.2331859074632265); C15 slope
(MASS0.1 2.8333818616053317 minus MASS0.0 3.0483380938795324) / 0.1 = -2.149562; C16 gamma
0.022487; C17 `s23/results/errdecomp.json :: rows[*]/f_common` mean 0.675770, median 0.673693,
n 126. All as E states.

Precision: E quotes stored values at stored precision; no number is quoted to more digits than
its artefact holds. The two last-digit pairs E lists separately are real (C01 3.214765154210998
against 3.2147651542109985; C03 `rmsd_fit` 3.20407616038092 in `baseline_tuning126.json`
against 3.2040761603809194 in `e_reproduce.json`): float summation order, not a discrepancy.

Basis: `s25/results/phys_suite.json :: basis` = "point_cloud". The artefact declares its own
basis, which is stronger than the inference in L28 item 3 (section 9).

The five UNSOURCED rulings:

- **C26 (36.1 / 36.4 deg): now DERIVED.** `s13/cache/tors_rows.npz` holds one JSON string per
  arm; arms `n_marg` (sequence-blind corpus marginal) and `p_grid` (full sequence context plus
  properties) each carry 126 rows with `err_phi` and `err_psi`, per-residue absolute errors
  (1,507 residues in total). Pooled over residues: phi MAE `n_marg` 36.416 deg, `p_grid` 36.133
  deg; psi MAE 72.772 and 62.355 deg. The dossier's 36.4 / 36.1 and 72.8 / 62.4
  (`s13/SPRINT13_DOSSIER.md:261-263`) are these numbers at one decimal. The mean of per-target
  means is 37.246 / 36.559; the record quotes the pooled statistic. Other arms reproduce the same
  way (`a_pepPos` 34.526 / 61.597, `p_point` 35.176 / 62.670). MINOR; the number is sourced now
  and may be quoted with this path.
- **C27 (+0.0004 / +0.0030): not absent, in git history.** `git ls-tree -r --name-only 5fa05cd`
  (the pre-consolidation tree named in the state brief) lists `s10/idaudit_price.json`,
  `s10/idaudit_audit.json`, `s10/idaudit_cause.json`, `s10/idaudit_overlap.json`,
  `s10/idaudit_pools.json`, `s10/idaudit.py`, `s10/test_idaudit.py`; the cited commit `49fc708`
  (`docs/FINDINGS.md:4948`) exists. I opened none of them: the price stage covers the sealed
  benchmark (the +0.0030 half), so the file may hold benchmark per-target values (Rule 1, rule
  14). Ruling: the +0.0030 half can never be re-derived in S26 and stays document-only, quoted
  only as "S10-4, historical"; the +0.0004 dev half is re-derivable on the dev set and sits inside
  lane W's mandatory direction (L25). MINOR.
- **C34 (|z_moment| 0.7529 / 0.8013 / 0.1127): confirmed document-only.** A grep of the three
  literals over `s25/` (py, log, txt, json, md) hits only `s25/LEDGER.md:1097` and
  `s25/agentPHYS_FINDINGS.md:115`; the `z_*` fields of `s25/results/calib.json` and
  `audit_t1.json` are posterior-calibration z-scores, a different quantity. Not quoted in the
  state brief. Re-derivable from the S25 channel energies by lane PH if a deliverable needs it.
  MINOR: not on a slide otherwise.
- **C35 (355 / 13): document-only, superseded.** The replacement is 370 / 357 / 13
  (`s26/results/test_run.json :: combined`), not the 369 / 356 that E writes; section 4. MINOR
  here; the examination's own restatement is the finding.
- **C24's 0.524: confirmed document-only.** `s8/integrate.py:1162` `stage_cvarcheck` writes
  `s8/integrate_cvarcheck.json` (absent). The three `_archive/logs/s8/*.log` hits for "0.524"
  are unrelated numbers (`chirality.log:16` `iv_geom_mirror_worse_frac` 0.524375,
  `diffuse_report.log:109` `gen_rama` +0.524, `predictor_report.log:73` `G_llfit` +0.524); no S8
  cvarcheck console log survives. The two sourced values for the same defect stand: 0.655634
  (`core/quantum.py:42`, `tests/test_quantum.py:59`) and 0.566586 (`q_verify.json`). MINOR: quote
  those, never 0.524; a re-run of `stage_cvarcheck` would give a seed-dependent sampled value.

C06 (+0.0103): not verified by me either (Rule 1). The two means beside it are pinned by
`tests/test_pipeline.py:338-347` (`BENCH60` constants, tolerance `RMSD_TOL`), passed in every
governed junit file (section 4).

## 3. The dataflow trace against the code

Three stages checked at HEAD:

1. **Score.** `core/predict.py:416` `self.w = shell / np.maximum((self.sd + 0.5) ** g, 1e-6)`
   with `(ws, g)` from `_score_weights()` (`:389-395`), which returns `np.ones(5), 1.0` because
   `score_weights.json` is intentionally absent (`:378-385`). The production weight is therefore
   `1/(sd_p + 0.5)`, as the state brief writes; E's module-map form `shell/(sd+0.5)^g` is the
   general code at its production values. `_risk` is an attribute, `(npairs, 760)` float32 at
   `:422`, on `np.arange(2.0, 40.0, 0.05)` (`:419`); the trace's `(105, 760)` and `(91, 760)`
   follow.
2. **Selector.** `core/pipeline.py:836-848`: `o = top[:128]`, `E = _zrank(sc[o])` (`_zrank` at
   `:788`, rankdata then standardise), `(alpha, T) = VQE_LFO[fold % 5]` (`:113`: alpha 1.0 on
   folds 0, 3, 4), readout `sel = o[consensus_medoid(block, p)]`, `sel_uniform =
   o[consensus_medoid(block)]`. The trace's `E_first5`, `E_last5`, the per-fold (alpha, T)
   (0.25 on fold 1, 1.0 on fold 3) and `sel` follow this code. The count behind "78 of 126" is
   the pinned fold sizes 25 + 23 + 30 (`s26/results/e_reproduce.json :: rows[*]/fold`: folds 0
   to 4 hold 25, 23, 25, 23, 30 targets).
3. **Average.** `core/pipeline.py:934` `b = int(np.argmin(cc._medoid(Ps)))`, `_medoid` = mean
   pairwise CA-RMSD to the other members (`s8/consensus2.py:293-295`), then superpose on member
   b and take the mean. The trace reports `medoid_local` and `medoid_pool_index` = `sub[b]`.

The identical pool index 449 on both targets (average stage and readout) looked like a trace
bug and is not: `bench_results/cache/1fc9f2dcf489e2fb/1S9Z.json :: sub[44]` = 449 and
`9KAR.json :: sub[34]` = 449; 449 appears in each `sub` once, at the position the trace names.
A coincidence of retrieval indices, verified against the records. The p-weighted medoid over the
128-block lands on the same member as the uniform medoid over the 75 on both targets, which is
what the identical `sel`, `sel_uniform` and `medoid_pool_index` say.

The AMBER row is read from the record, as E states, and matches the cache: 1S9Z `amber_e0`
5370.044288286983, `amber_e1` -1290.59114035452, `amber_moved` 0.19126831218030266,
`amber_strain_after` 51.98570245549598, `amber_ca` bond mean 3.8618725; 9KAR
6020473359349.543, 1262.4103583126937, 0.6096586566244366, 1166.1306874049387, 3.9338333. 9KAR
is the largest `rmsd_arm` minus `pool_best` gap in the cache (6.038; next 2NDM 5.236, 2MQ2
4.837; `e_reproduce.json` rows), as E says. Note, not a finding: `average()` breaks a tie in the
medoid criterion by array order (`np.argmin`); the criterion is a float mean of pairwise RMSDs,
a tie is measure-zero, and the code is pinned production from before S26.

The trace describes what the code does. No finding.

## 4. The test run

Facts from `s26/jobs_done/pytest_*.json` and the junit files (parsed with `xml.etree`):

| job | tree | junit file: tests / passed / failed / errors / skipped | wall s | peak RSS GB | exit |
|---|---|---|---|---|---|
| `pytest_core` | a4db170c | `s26/results/pytest_core.xml` 350 / 337 / 0 / 0 / 13 | 195.3 | 1.694 | 0 |
| `pytest_amber` | 37bddbbb | `pytest_amber.xml` 16 / 16 / 0 / 0 / 0 | 260.5 | 0.872 | 0 |
| `pytest_nanpoison_post_core_edit` | 37bddbbb | `pytest_nanpoison_post_core_edit.xml` 79 / 77 / 0 / 0 / 2 (not in the total) | 225.5 | 1.566 | 0 |
| `pytest_amber_frame` | 37bddbbb | `pytest_amber_frame.xml` 3 / 3 / 0 / 0 / 0 | 255.7 | 0.324 | 0 |
| `pytest_core_post` | 601a39c7 | `pytest_core_post.xml` 351 / 338 / 0 / 0 / 13 | 240.4 | 1.692 | 0 |

`s26/TEST_RUN.md` (rendered 08:38) and `s26/results/test_run.json :: combined` say 370 tests,
357 passed, 0 failed, 0 errors, 13 skipped, 0 memory-guard skips = `pytest_core_post` +
`pytest_amber` + `pytest_amber_frame`, with `pytest_core` marked superseded. L20 and the state
brief section 6.1 say the same. The 13 skips in `pytest_core_post.xml`: 3 in
`test_equivalence.py` and 8 in `test_integration.py` with message "set VERIFY_SLOW=1 ...", and 2
absent-artefact skips in `test_pipeline.py` (`bench_results/optimised_tuning126_w6.json`; "no
smoke result on disk"); no skip message names the memory ceiling. E's account of the skips is
correct. The tests E relies on in section C (C04, C05, the NaN-poison guard, the
pre-registration pin, the four-component reproduction, the tail-baseline regression) are
`passed` in `pytest_core.xml` and in `pytest_core_post.xml`.

**FINDING A1, MATERIAL.** `s26/EXAMINATION.md` section D reads "Combined: 369 tests, 356
passed, 0 failed, 0 errors, 13 skipped", built from `pytest_core` (350 / 337) plus the two AMBER
jobs; it lists four junit files and omits `pytest_core_post.xml`; it names `s26/TEST_RUN.md`
"rendered 2026-09-13 00:41", the render that predates `pytest_core_post` (00:41:56 to 00:45:56).
Row C35 then states "the governed S26 run (`s26/TEST_RUN.md`) is 369 tests, 356 passed, 13
skipped"; `s26/agentE_FINDINGS.md` X1 and item 4 of "Claims I could not source" repeat it;
`s26/LEDGER.md` L28 item 2 repeats it ("369 / 356 / 13"). The artefact named holds 370 / 357 /
13. By the audit rule (a path that exists but does not contain the number) this is material. It
is a stale read of a superseded render, not a wrong measurement; the correct count is already in
L20, `s26/TEST_RUN.md` and the state brief; nothing in Phase 1 depends on it. Fix: lane E (or the
coordinator, if E has finished its turn) appends a correction to `s26/EXAMINATION.md` (section D
and row C35) and to `s26/agentE_FINDINGS.md`, and posts a ledger entry naming L28 item 2; I
re-check on the append. `s26/BRIEF.md` does not carry the count and needs no change.

## 5. The declared defects

All five are restated with file and line, and every citation resolves at HEAD (`grep -n`):

1. `core/data.py:182` `def identity(..., norm="longer")`; `:216` `containment`; `:321`
   `identity_many`; `peptide_db.py:119` `def identity`, `:170` `shared / max(len(a), len(b))`;
   `tests/test_data.py:192-211` (1CEK in 1A11, 13/25 = 0.520) and `:74` `ok_canon >= 4`. The
   state brief's `core/data.py:188` points inside the same function (at `ae86a124` line 188 is
   `n, m = len(a), len(b)`).
2. `core/project.py:149-151, 155` and `core/geometry.py:57-65` constants; `OMEGA_TRANS =
   math.pi` at `project.py:155` and `geometry.py:65`; the step gate is `core/data.py:439`
   (`step.min() < 3.5 or step.max() > 4.1`) and `:725` for windows. Lane PH's L22 cites
   `core/data.py:406-407` and `697-698`, the line numbers before `37bddbbb`; same content, noted
   for the L22 check, not a finding here.
3. `s25/resultslab/exportlib.py:598` `pool_oracle_gate` (`difficulty_gate` at `:653`);
   `schema.py:309-310`; `results/summary/leaderboard.json` production row `pool_gate` WARN,
   `n_pool_violations` 2, `pool_margin` 1.501796830456039.
4. `core/predict.py:224` `_fold_fragments`, `:345` `train_fold`; the five model hashes are in
   `s26/results/pinned_hashes.json :: files/distogram_models/*`.
5. `s25/phys_lib.py:98` `zmoment`; `core/pipeline.py:788` `_zrank`, applied at `:838`;
   `core/amber.py:1282` `CONVERGE_MAX_KCAL = 1000.0`; `s8/integrate.py:292` `if cm["bond"] +
   cm["angle"] > 1000.0`.

No finding.

## 6. The pinned-file hash table

Recomputed with `sha256sum` and compared to `s26/results/pinned_hashes.json ::
files/<path>/{sha256, bytes}`:

| file | recomputed sha256 | bytes | table |
|---|---|---|---|
| `peptide_folds.json` | 89bc58edc59c6d794cd83f093cdafe74f05075aa43251f9ed3adfbe8115e4683 | 18956 | equal |
| `peptide_clusters.json` | 9ae3d8f3afc281df7d37319fb04420cd22976032696f7c27f517de79549c75b9 | 20217 | equal |
| `catrace_prior.npz` | d93a5c627d5ca4dfb8b6c358249415266c9f320b8cbb35caa7faf6a7449a1bf8 | 7172 | equal |
| `distogram_models/fold1_esm_frag.pt` | 23691cc531435766c88c394bee1650515067249a804fa4ce04b5e6947f9d29fd | 1497519 | equal |

Four of four. I did not re-hash `results/benchmark_manifest.json` (rule 14); E's byte hash
equals the S20 record per `benchmark_manifest_matches_S20_record` = true, and
`s26/e_hashes.py:55-60` reads bytes in chunks and never parses. No finding.

## 7. What is not in git

`du -sb`: `bench_results/cache/1fc9f2dcf489e2fb` 715,001 bytes, 126 files (E: 0.7 MB / 126);
`s8/generate_univ` 117,402,802 bytes, 126 files (E: 117.4 MB / 126); `esm_small.npz` 22,062,826
bytes (E: 22.1 MB; the hash table's 22062826). `git check-ignore -v`: `.gitignore:19`
(`bench_results/*`), `:90` (`s8/generate_univ/`), `:113` (`*.npz`); `git ls-files` returns 0
files for the first two. `s26/models/` is no longer empty (lane P's checkpoints); trivial. No
finding.

## 8. Rule 1 closure

`grep -rn "final_report\|benchmark_manifest\|benchmark(" s26/*.py` hits three scripts:
`s26/e_claims.py` (docstring `:11`, comment `:28-30`, `SKIP_FILES` `:32-33` = the manifest,
`results/monomer_manifest.json`, `s9/final_report.json`, `s9/final_synth.json`; `SKIP_SUBSTR`
`:34` = "benchmark", "s9/final_cache", "bench60", "benchmark60"; applied at `:139`; committed in
`e3570aed`, unchanged at HEAD), `s26/e_hashes.py` (`:31` the manifest in the pinned list, `:92`
its entry; bytes only, read at `:55-60`), and `s26/i_identity_audit.py` (`:262`, `:375`
`core.backend("data").benchmark()`, the sequence-only count read that L18 declares and L18b
rules conservative). No other S26 script reads a benchmark file. The incident is closed in the
committed script.

Two observations for the coordinator, not findings of this audit: (a)
`tests/test_pipeline.py:338-347` parses `s9/final_report.json` (three aggregate means and `n`),
and `test_reproduces_the_committed_benchmark_synthesis_bit_for_bit` reads `s9/final_synth.json`
and `s9/final_cache` (coordinates, no RMSD); these pinned regression tests ran three times under
lane I before L29's rule, and any further full `pytest tests/` in S26 parses them again. Say
whether that is accepted or the test is `--deselect`ed for the rest of the sprint. (b)
`core.data.benchmark()` is reached indirectly by `dev_set()` (`core/data.py:633`) and
`s7/debias.py:110`: every construction of the 126-target instrument reads the manifest's `pdb`
list to exclude it. That is the instrument's pre-existing definition, not an S26 read.

## 9. Other checks (no findings)

- L28 item 3 against the artefacts: `s25/results/phys_suite.json :: configs_rank/{3,5,4,7,6,1,2}/mean`
  as in the C29 row, and `:: basis` = "point_cloud"; `results/summary/leaderboard.json :: rows[*]`
  for distogram, amber_distogram, legacy_distogram, legacy_amber_distogram, legacy_amber, legacy,
  amber: `mean` (basis `built_chain_bb`) 3.218708055115966, 3.309961759812365,
  3.3731520305329474, 3.4220825100678876, 3.8248241126995124, 3.8843706034949403,
  4.101517931701575; `mean_secondary` (basis `point_cloud_ca`) 3.0580168439307194,
  3.1316921198294936, 3.2150989596380146, 3.253030982555741, 3.6742099693114274,
  3.7553479973755235, 3.8805028264654857 (the suite's means after the PDB round trip). The
  random null `random_null_rank/mean` 3.4251256499338507 has no built-chain counterpart on disk.
  L28 item 3 stands and is understated: the artefact declares its basis. The state brief's
  section 2 table and section 5.3 (`docs/STATE_BRIEF_2026-09-12.md:57, 190-191`) put these
  beside the built-chain 3.215 without a basis label; L29 item 2's rule covers every S26
  deliverable.
- `s26/BRIEF.md` numbers: A1's reference `s26/results/q_mde_reference.json ::
  contrasts/rmsd_q_synth-rmsd_u_synth` effect -0.01346251788977831, se 0.03420495699920276, mde
  0.09582860752896645, 0.14x (its `source_cache` `bench_results/cache/464a0ddb5f283e04` is an
  S25 cache dated 2026-09-04, so no S26 endpoint ran for it); B3 from L14; C1 from L26; C4's
  m-router MDE 0.055 = 2.8016 x 0.0196; C5's f_common 0.6758; C3's +0.0207 = 3.2355 - 3.2148.
  All consistent with their sources.
- H2 / "78 of 126": `s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint`
  0.6190476190476191 x 126 = 78.0, equal to the pinned count of folds 0, 3 and 4.

## Verdict

One material finding (A1, section 4: the examination's test-suite total 369 / 356 is not the
370 / 357 its named artefact holds; documentary fix, then re-check). Minor: C27's artefacts are
in git history at `5fa05cd`, not absent; C26 is now derived from `s13/cache/tors_rows.npz` and
matches; C34, C35 and the 0.524 stay document-only with the rulings above.

AUDIT: MATERIAL FINDINGS 1, LISTED ABOVE
