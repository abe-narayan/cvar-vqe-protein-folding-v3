# S26 STATUS

Each lane appends two lines under its heading every 60 minutes: what it is running, what it
will do next. Timestamp every entry (local time). Never edit another lane's lines.

## Coordinator
- 2026-09-12 23:50 running: governor started (pid 36196); cache_amber committed; reading the record.
- 2026-09-13 01:30 next: spawn lanes E, I, Q, P; continue reading core modules; adversary after EXAMINATION.md lands.
  (clock correction: the line above was written at about 00:03, not 01:30; lanes E, I, Q, P were spawned 00:03 to 00:06.)
- 2026-09-13 00:20 running: lanes E, I, Q, P active; coordinator finished reading predict/data/walsh/geo/S25 physics; committed contract+ledger (74e073e2); spawning PH.
- 2026-09-13 00:20 next: Adversary when s26/EXAMINATION.md lands; Presentation lane when the first proposal verdict is in; Wildcard after that.

- 2026-09-13 08:40 running: all five lanes resumed after the 00:46 session-limit cut (L17); governor continuous; no job in flight at the cut.
- 2026-09-13 08:40 next: Adversary when EXAMINATION.md lands; PR lane on the first proposal verdict; commit shared ledger/status hourly.

## E (Examiner / Librarian)
- 2026-09-13 00:12 running: e_reproduce (126 records through s12.instrument), e_trace_1S9Z then e_trace_9KAR under jobrun; module_map.json (699 modules) and pinned_hashes.json (benchmark manifest matches S20) written.
- 2026-09-13 00:12 next: finish reading core/energy, core/amber, the 24 root modules, s5/s7/s8/s9, tests, verify, resultslab; then write EXAMINATION.md sections A-H.
- 2026-09-13 08:42 running: resumed after the 00:46 session cut; scripts and artefacts committed; reading ledger L6-L16b and the PREREG files; assembling EXAMINATION.md (A-H), then BRIEF.md and agentE_FINDINGS.md.
- 2026-09-13 08:42 next: ledger entry for the unsourced claims and the incident notes; commit the three deliverables; two-line summary for the coordinator.

## I (Infrastructure / Governor)

- 2026-09-13 00:13 running: read contract/brief/code; job pytest_core (non-AMBER suite, tag TEST, est 2.5 GB) launched under the governor at ~00:11; brief committed unchanged.
- 2026-09-13 00:13 next: AMBER test files as two sequential AMBER jobs; S8-13 source excerpt; identity audit script (in memory only); run_equiv2.sh fix.
- 2026-09-13 00:42 running: run_equiv2 (AMBER slot, --no-amber, 4 arms of smoke8) and pytest_core_post (final non-AMBER suite on 601a39c7). Done: pytest_core 337/13/0, pytest_amber 16/0, frame 3/0; identity audit (in memory) L15; items 3, 4 closed (L6, L10); 6b/6c/6d entered; 3 commits (eb89c165, 37bddbbb, 601a39c7).
- 2026-09-13 00:42 next: compare the equiv2 arms, write s26/results/run_equiv2.log, brief 6.1 edit + commit, examine.py full run (regenerates s26/results/module_map.json, lane E: it will be rewritten once), agentI_FINDINGS.md.
- 2026-09-13 08:40 running: resumed after the 00:46 session cut (nothing restarted). Recorded pytest_core_post (final suite: 370 tests, 357 passed, 13 skipped, 0 failed) and run_equiv2 (4 arms, baseline = exact bit-identical 8/8); brief 6.1 replaced and committed; ledger L17 (benchmark provenance), L18 (item 2).
- 2026-09-13 08:40 next: examine.py now also runs e_hashes --check; full examine run (rewrites s26/results/module_map.json once); closure table in agentI_FINDINGS.md; ledger notes for item 5 and hygiene; final commit; then finish.
- 2026-09-13 08:43 running: nothing; examine.py full run passed (725 modules, no drift on 101 pinned entries, 21/21 claims); ledger L20 (item 5), L21 (hygiene); findings closure table written.
- 2026-09-13 08:43 next: closing commit (findings, examine.py, claim_check.json), then the lane finishes its turn. Items 1-5 and defects 6a-6d closed; no production number changed.

## Q (Quantum lane)
- 2026-09-13 08:50 running: resumed after the session cut at 00:46; q_*.py, probes and reference MDEs committed (2fb2998b); PREREG_A2/A4 on disk; launching A2 (DLA, s26/q_dla.py) and A4 (variance, s26/q_var.py) under the governor per the coordinator ruling.
- 2026-09-13 08:50 next: PREREG_A1/A3, IDEA files, TRAINABILITY_PAPER_OUTLINE.md, agentQ_FINDINGS.md; A1/A3 endpoints wait for PHASE 0 SIGNED OFF.

- 2026-09-13 09:25 running: A4 (s26/q_var.py, job a4_var) grown-circuit sweep in progress; A2 complete (s26/results/q_dla.json, ledger L27: full so(2^n) from depth 2 at n=7, my H2b falsified). PREREG_A1-A4, three IDEA files, TRAINABILITY_PAPER_OUTLINE.md, agentQ_FINDINGS.md on disk and committed (05f39d2f).
- 2026-09-13 09:25 next: A4 results into the findings and a ledger entry; A1/A3 builds and labels wait for PHASE 0 SIGNED OFF (harness ready, --label self-gated).
## P (Prior / Learning lane)
- 2026-09-13 00:20 running: read LANE_CONTRACT, state brief, S24 ledger (all), S25/S23/S22/S17/S19 sections, S7/S8/S9/S10 findings, S12 dossier+agg, core/predict, data, distogram, esm_features, pairnet, s7/repr_select, priorladder, instrument, stats_lib, errdecomp; probes p_probe_esm (CPU, 1.0 GB) then p_probe_esmcache (ESM, 2.5 GB) launched under jobrun.
- 2026-09-13 00:20 next: write PREREG_C1..C5, PREREG_B1..B3, three IDEA files, s26/p_ladder.py + synthetic unit test, C2 one-fold training-time probe, B1 feasibility json, agentP_FINDINGS.md; no endpoint runs before sign-off.
- 2026-09-13 01:00 running: nothing (phase gate). Written: PREREG_C1..C5, PREREG_B1..B3, IDEA_better_prior_inputs / window_ensembling / coherence_penalised_training, p_ladder.py (11 rungs, tests pass), b1_feasibility.json (B1 infeasible: RAM 8.76 GB vs 4.4 GB headroom, openfold/omegaconf absent), probes done (esm_cache 1.79 GB; pca32 fold-0 retrain 453 s / 1.25 GB reproduces the pinned model, argmin 25/25); ledger L11-L14.
- 2026-09-13 01:00 next: awaiting sign-off (question in L12: may the 5-fold ladder TRAINING run before sign-off?); then eval shipped -> noesm -> conly -> pca32 -> mix -> pairnet -> pca32f -> pca128 -> raw -> wide -> esm8m, C1 write-up, B3, C4, C5, PROPOSAL_C/B/B_REPLACEMENT.
- 2026-09-13 08:45 running: resumed after the 00:46 session cut (no training job had launched; only the probe model exists). Training order fixed in PREREG_C2 addendum 2; launching the 8-rung training chain one governor job at a time (est-ram 1.5 GB); PREREG_C3 addendum 2 and IDEA_amber_prior_partner.md written.
- 2026-09-13 08:45 next: C1 closure write-up; B3/C4/C5 code with synthetic tests; agentP_FINDINGS sections; ledger entry for the launch (re-reading the tail first); commit; hourly STATUS.
- 2026-09-13 09:25 running: C2 training chain (p_train_noesm in progress, then conly, pca32, wide, pca32f, pca128, esm8m, raw; one governor job at a time; ~6-10 h). Committed 19ca2ea3: C1 write-up, p_stats/p_b3/p_c4/p_deliver, features cached, IDEA_amber_prior_partner, PREREG_C2/C3 addenda; ledger L26.
- 2026-09-13 09:25 next: turn ends once p_c5 selftest is green and committed; after "PHASE 0 SIGNED OFF": eval shipped then each rung as its checkpoints complete, B3/C4/C5 runs, PROPOSAL_C/B/B_REPLACEMENT, p_deliver for the best rung (ledger + finish turn for PH).

## A (Adversary) -- not yet spawned

## PH (Physics lane)
- 2026-09-13 08:41 running: nothing (gate closed). Done: reading list; PREREG_amber_reject/cis/c3_control; five IDEA files; ph_lib/ph_cis/ph_reject/ph_c3 + synthetic tests; census jobs ph_cis_census (0.038 GB), ph_reject_census (0.117 GB), ph_c3_nativefree; ledger L22-L24; findings and Part IV notes; commit 5dc7a3a6.
- 2026-09-13 08:41 next: READY TO LAUNCH the minute PHASE 0 SIGNED OFF lands: `ph_cis.py floor` (1 min), `ph_c3.py stage1` (2 min), `ph_reject.py cloud` (5 min) then `ph_reject.py chain` (CPU, ~3 h, per-target cells); C3 stage 2 when P delivers (AMBER, probe first). Awaiting tournament ranking for IDEA_branch_select / rotamer_relief part B.
- 2026-09-13 09:35 running: nothing (gate closed). Done since 08:41: PREREG_branch_select.md + ph_branch.py and PREREG_rotamer_relief.md + ph_relief.py (probe modes, per-target cells, 7/7 synthetic tests pass); reject report gains the moved-subset secondary and the empty-draw fallback; findings section 2.3 states the R-first reading; Part IV notes checked against the six prompt items, Walsh numbers and L23 added, speaker-notes table (section 8).
- 2026-09-13 09:35 next: on sign-off: floor, C3 stage 1, reject cloud then chain (as above). If ranked: `jobrun --tag AMBER --name ph_branch_probe --est-ram 1.0 -- python s26/ph_branch.py probe --pdb 1A13` after `ph_branch.py solutions` (CPU, 10 min); `jobrun --tag AMBER --name ph_relief_probe --est-ram 0.8 -- python s26/ph_relief.py probe --pdb 1A13 --m 5`. C3 stage 2 on s26/results/p_best_rung_chains.json when P delivers.

## PR (Presentation lane) -- not yet spawned

## W (Wildcard)
