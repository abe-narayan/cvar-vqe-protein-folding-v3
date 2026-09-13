# agentE_FINDINGS (Sprint 26, lane E: Examiner / Librarian)

Format as S12 to S25: tiers DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN,
every number with its artefact path, then "what damaged my own expectations" and "what I did
not do and why". Companion documents: `s26/EXAMINATION.md` (sections A to I) and `s26/BRIEF.md`.
Scripts and artefacts committed in `e3570aed`.

## DEMONSTRATED

D1. **The production result reproduces exactly from the persisted records through an
independent RMSD implementation.** `s26/results/e_reproduce.json` (`s26/e_reproduce.py`, job
`e_reproduce`, exit 0): fresh vs stored means 3.048338093879532 / 3.2040761603809194 /
3.2147651542109985 / 3.2354598538973844 on `rmsd_avg` / `rmsd_fit` / `rmsd_arm` / `rmsd_full`,
max per-target disagreement 0.0 on all four, n = 126, natives from `core.backend("data").load()`
through `core.pipeline._q`, RMSD through `s12.instrument.ca_rmsd`. T030 = 1S9Z `rmsd_arm`
0.18198112330908295 both ways. Fold labels match the pinned folds on all 126 records.

D2. **A fresh `run_target` + `label` reproduces the record on both traced targets to 0.0.**
`s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json`: shipped, pool_best, top_m_best,
rmsd_avg, rmsd_fit, rmsd_arm all diff 0.0; `sub` identical; emitted `ca` max abs diff 0.0;
every intermediate (avg_ca, fit_ca, ca, phi, psi) matches the record at 0.0.

D3. **The selector's Hamiltonian is the standardised rank ladder to 0.394% of range on 1S9Z
(2 ties) and 9KAR (4 ties)** (`stages/hamiltonian/dev_frac_of_range` 0.003938280770824405 and
0.0039392382627650765), consistent with S25's worst case of 1.18% over 8 targets
(`s25/results/q_gibbs.json :: results/spectrum_target_independence/worst_frac_of_range`
0.011821117271971639).

D4. **The sealed benchmark manifest is unchanged.** `s26/results/pinned_hashes.json`:
`results/benchmark_manifest.json` sha256 a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d,
8002 bytes, equal to the S20 record; hashed as bytes, never parsed. All other pinned files hashed
(table in EXAMINATION F); lane I's audit found the folds and clusters equal before and after its
in-memory re-clustering (L15).

D5. **The claim search sourced 30 of the 35 listed numbers to an artefact leaf, a passing test, or
a stated derivation from stored rows**, with file and key (EXAMINATION C;
`s26/results/claim_search.json`, `.txt`). Examples at stored precision: 1.7108244199364904
(`bench_results/baseline_tuning126.json :: science/pool_best/mean`), 3.4251256499338507
(`s25/results/phys_suite.json :: random_null_rank/mean`), -0.7422836835050723
(`s25/results/q_alpha.json :: results/cell_correlations/H`), 0.9019158325027764
(`s25/results/q_gibbs.json :: results/divergence_ladder/1/ladder/"50 Adam steps (DEPLOYED)"/kl_nats`),
2592 and 0 (`s25/results/q_verify.json :: results/n_cells`, `"SUBSET-hood violations  <- THE THEOREM"`).

D6. **The prior-ladder slope and the gamma for 3.0 A are arithmetic on stored rows, not stored
leaves.** `s24/results/priorladder.json` rows: mean MASS0.0 3.048338, MASS0.1 2.833382, MASS1.0
2.226080; slope on the first rung -2.149562 A per unit gamma; linear interpolation gives gamma
0.022487 for 3.0 A. Both match the ledger's -2.1496 and 0.0225.

D7. **The 68% common-mode figure is the mean of a per-row leaf.** `s23/results/errdecomp.json ::
rows[*]/f_common`: mean 0.675770, median 0.673693 (n = 126).

## ORACLE DIAGNOSTIC

O1. **9KAR shows the selection failure mode in numbers**: pool best 1.3997329473495483, top-75
best 5.343628883361816 (the score filter keeps nothing near the native), coordinate average with
virtual bonds of mean 1.9696 A and minimum 1.1556 A, projection to 3.804 A, AMBER from
6.020473359349543e12 to 1262.4103583126937 kcal/mol (above the 1000 kcal/mol convergence gate),
`rmsd_full` 7.543710007423382 (`s26/results/e_trace_9KAR.json`).

O2. **One of the 126 production relaxations ends above `core.amber.CONVERGE_MAX_KCAL` and one
more fails the S8 strain-rejection rule; both are inside the 3.2355 mean**: 9KAR (e1
1262.4103583126937 > 1000, bond+angle strain 1166.1306874049387, moved 0.6096586566244366) and
2BP4 (e1 845.3467373422843, converged by the gate; strain 1172.6967039637698 > 1000, moved
0.503434471788075); `bench_results/cache/1fc9f2dcf489e2fb/{2BP4,9KAR}.json`; lane PH's L24
counts 125 of 126 converged. Mean `amber_moved` over 126 is 0.2332 A, max 0.6097 (9KAR).

O3. **Point cloud contraction re-measured**: mean virtual bond 2.9613938229671057 vs native
3.812176007269521 (22.317%), worst bond 0.6491799200244281, 80/126 under 3.4 A; built chain
3.803954938363982 with sd at most 1.802075856660821e-15 (`e_reproduce.json :: summary/virtual_bond`).

## HYPOTHESIS

H1. The point-cloud basis of the seven-configuration suite, of the random-75 null and of the
Legacy/AMBER "worse than noise" effects (+0.330 / +0.455; `s25/results/phys_suite.json ::
controls_rank/{1,2}/vs_random/effect` 0.3302158940209069 / 0.45537105631568636) is not carried
into the presenter documents; on the built chain the same leaderboard rows read 3.2187 to 4.1015
(`results/summary/leaderboard.json :: rows[*]/mean`). The ranking would survive a basis change;
the magnitudes would not. To be checked by the Adversary before any slide quotes them.

H2. Because alpha = 1.0 on folds 0, 3 and 4 (`core/pipeline.py:113`; 78 of 126 targets,
`s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint` 0.6190476190476191),
any Proposal A statement about "the CVaR tail" describes 48 targets unless the arm re-derives the
(alpha, T) table.

## REFUTED

None. No number I re-measured disagreed with its record.

## OPEN

X1. **Four listed numbers have no artefact on disk** (EXAMINATION C): 36.1 / 36.4 deg (C26;
cited `s13/SPRINT13_DOSSIER.md:261-263`, `s13/tors_FINDINGS.md:373, 385`); +0.0004 / +0.0030
(C27; `docs/FINDINGS.md:4573, 5026`; the cited `s10/idaudit_*.json` do not exist, there is no
`s10/`); the |z_moment| triple 0.7529 / 0.8013 / 0.1127 (C34; `s25/LEDGER.md:1097`,
`s25/agentPHYS_FINDINGS.md:115`); 355 passed / 13 skipped (C35; `docs/STATE_BRIEF_2026-09-12.md:274`,
superseded by `s26/TEST_RUN.md`: 369 tests, 356 passed, 13 skipped). One component of C24, the
0.524 cosine of the sampled tail-only estimator (`docs/FINDINGS.md:3161`), has no artefact
(`s8/integrate_cvarcheck.json` is absent). `s13/cache/tors_rows.npz` may hold the phi MAE values
of C26; I did not open npz caches in the claim search.

X2. **The production cache `bench_results/cache/1fc9f2dcf489e2fb/` is gitignored** (EXAMINATION
G). Every S12 to S26 reproduction, the results lab's incumbent and the lane pre-registrations
read it. Regenerating it needs the untracked inputs (`pdbs_ext/`, `prots/`, the ESM caches, the
pinned models) and one AMBER pass.

X3. **The benchmark paired delta +0.0103 is named to `s9/final_report.json` but not verified by
this lane** (Rule 1). The two means it sits between are asserted to 5e-4 by two passing tests.

## What damaged my own expectations

- I expected the claim search to be a grep. Numbers at the quoted precision are cited in
  documents far more often than they are stored, and stored leaves round to the same value by
  coincidence in the per-target arrays of every sprint (C09 4.065 has 497 leaves that round to
  it, one of which is the constant helix). Only the aggregate-leaf filter and the targeted
  per-sprint lookups separated citation from artefact; the first pass would have marked several
  sourced claims as unsourced and one unsourced claim (C34) as sourced.
- I expected the physics-suite and the state brief to share a basis. They do not (H1).
- I expected the `jobrun` peak RSS to be usable for every job. For a 5-second job the 5-second
  sampler read 0.004 GB; the in-process `pl.proc_rss()` peak (46.09 MB) is the honest figure and
  jobs shorter than the sampling interval must report it themselves.
- I expected the working tree to be frozen during Phase 0. Lane I's ledgered edits to `core/`
  (L10, L15, L16; commit `37bddbbb`) landed between my traces and this write-up. The identity
  path they touch is not on `run_target` and the tests passed on the edited tree, but my
  reproduction stands on the earlier tree and should be re-run on HEAD (5 s).
- My first `e_claims.py` run parsed `s9/final_report.json` and printed its aggregate
  `dist/shipped/mean` (2.9507235775391263), a number already in `README.md:14`. I had excluded
  the manifest and the per-target cache but not the report. The exclusion is now in the script
  (`SKIP_FILES`, `SKIP_SUBSTR`), the artefacts were regenerated, and no per-target benchmark
  value was printed or read. Recorded in ledger L28.

## What I did not do and why

- Did not run `pytest` (contract section 3; lane I owns the run, `s26/TEST_RUN.md`).
- Did not re-run AMBER on either traced target (RAM was under the box's budget but the brief
  allowed reading `amber_ca` from the record; the trace scores the stored relaxation).
- Did not open `s9/final_report.json`, `s9/final_synth.json`, `s9/final_cache/` or the manifest
  (Rule 1); C04 to C06 are sourced by the tests that assert them.
- Did not open `.npz` caches in the claim search (C26 may live in `s13/cache/tors_rows.npz`).
- Did not compute built-chain means for the random-75 null (it would need 126 projections; an
  endpoint computation, gated).
- Did not commit `s26/STATUS.md` (shared, coordinator-owned); my hourly lines are appended in
  the working tree.
- Did not write the lane Q pre-registration rows from their files (not on disk at assembly);
  BRIEF section 1 is written from the `s26/q_*.py` docstrings and `s26/results/q_mde_reference.json`.

## Claims I could not source (the list the coordinator asked for)

1. 36.1 deg vs 36.4 deg phi MAE (C26): document-only.
2. +0.0004 / +0.0030 identity-leak prices (C27): document-only; cited artefacts absent.
3. mean |z_moment| 0.7529 / 0.8013 / 0.1127 (C34): document-only.
4. 355 passed / 13 skipped (C35): document-only; superseded by `s26/results/test_run.json`.
5. 0.524 sampled tail-only cosine (part of C24): document-only.


---

## ADDENDUM (appended 2026-09-13 09:05, coordinator, after L31 item A1)

X1 above quotes 369 / 356 / 13 for the governed test run; the committed artefact holds
370 / 357 / 13 (`s26/results/test_run.json :: combined`). See the addendum in `s26/EXAMINATION.md`.
