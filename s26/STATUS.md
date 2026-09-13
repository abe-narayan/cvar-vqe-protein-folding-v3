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

- 2026-09-13 08:54 running: Adversary auditing EXAMINATION.md (spawned 08:55); lane E moved to report writing; P training rung noesm; Q A4 running; PH and W pre-gate work done or in progress; lane I finished.
- 2026-09-13 08:54 next: PHASE 0 SIGNED OFF on a clean audit, then release Q (A1/A3), P (ladder eval), PH (cis floor, C3 stage 1, reject) and W; PR lane on the first verdict.

- 2026-09-13 09:19 running: PHASE 0 SIGNED OFF (L33) after the audit's one material item was fixed (L32); Q, P, PH, W released to endpoints; A re-checks then attacks L27/L22-L24/L30 and ranks the tournament; E drafting the report.
- 2026-09-13 09:19 next: spawn PR on the first proposal verdict; relay P's best-rung file to PH; commit shared files hourly; watch the four-job launch cap.

## E (Examiner / Librarian)
- 2026-09-13 00:12 running: e_reproduce (126 records through s12.instrument), e_trace_1S9Z then e_trace_9KAR under jobrun; module_map.json (699 modules) and pinned_hashes.json (benchmark manifest matches S20) written.
- 2026-09-13 00:12 next: finish reading core/energy, core/amber, the 24 root modules, s5/s7/s8/s9, tests, verify, resultslab; then write EXAMINATION.md sections A-H.
- 2026-09-13 08:42 running: resumed after the 00:46 session cut; scripts and artefacts committed; reading ledger L6-L16b and the PREREG files; assembling EXAMINATION.md (A-H), then BRIEF.md and agentE_FINDINGS.md.
- 2026-09-13 08:42 next: ledger entry for the unsourced claims and the incident notes; commit the three deliverables; two-line summary for the coordinator.
- 2026-09-13 09:26 running: report writer; s26/REPORT.md Parts I to VI, IX, Appendix A and the B rows committed part by part (stub commit 5aaea13c was a path error, content followed in 3312996a); addenda for L31 to L34 applied (C26 artefact, Phase 0 signed off).
- 2026-09-13 09:26 next: fill Part VII / VIII / Appendix C slots and the Part V S26 additions as verdicts land from L33 on; keep the hourly STATUS; the 09:40 stamp on the previous pair was a clock error, the lines were written about 09:15.
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
- 2026-09-13 09:35 running: gate open (L33); A1 launched as two governed shards (a1_build_s0/s1, est-ram 0.4, threads 1, per-target checkpoints s26/results/a1/); A3 chained to start when both A1 shards finish (s26/q_chain_a3.py); A4 (a4_var) on its last cell.
- 2026-09-13 09:35 next: A2 figure (s26/figures/), PROPOSAL_B_REPLACEMENT.md drafted; A4 findings + figure when it lands; A1 label/stats, ledger entries with ST.fmt, replications, PROPOSAL_A.md.
- 2026-09-13 09:45 running: A1 shards a1_build_s0/s1 (63 targets each) building under the governor; A3 chained behind them. A4 complete (s26/results/q_var.json, ledger L35, figure s26/figures/a4_variance_slopes.png).
- 2026-09-13 09:45 next: findings section 4 and outline C9/F9 from A4; then A1 label + stats + ledger (ST.fmt verbatim) when the shards finish; PROPOSAL_A.md.
## P (Prior / Learning lane)
- 2026-09-13 00:20 running: read LANE_CONTRACT, state brief, S24 ledger (all), S25/S23/S22/S17/S19 sections, S7/S8/S9/S10 findings, S12 dossier+agg, core/predict, data, distogram, esm_features, pairnet, s7/repr_select, priorladder, instrument, stats_lib, errdecomp; probes p_probe_esm (CPU, 1.0 GB) then p_probe_esmcache (ESM, 2.5 GB) launched under jobrun.
- 2026-09-13 00:20 next: write PREREG_C1..C5, PREREG_B1..B3, three IDEA files, s26/p_ladder.py + synthetic unit test, C2 one-fold training-time probe, B1 feasibility json, agentP_FINDINGS.md; no endpoint runs before sign-off.
- 2026-09-13 01:00 running: nothing (phase gate). Written: PREREG_C1..C5, PREREG_B1..B3, IDEA_better_prior_inputs / window_ensembling / coherence_penalised_training, p_ladder.py (11 rungs, tests pass), b1_feasibility.json (B1 infeasible: RAM 8.76 GB vs 4.4 GB headroom, openfold/omegaconf absent), probes done (esm_cache 1.79 GB; pca32 fold-0 retrain 453 s / 1.25 GB reproduces the pinned model, argmin 25/25); ledger L11-L14.
- 2026-09-13 01:00 next: awaiting sign-off (question in L12: may the 5-fold ladder TRAINING run before sign-off?); then eval shipped -> noesm -> conly -> pca32 -> mix -> pairnet -> pca32f -> pca128 -> raw -> wide -> esm8m, C1 write-up, B3, C4, C5, PROPOSAL_C/B/B_REPLACEMENT.
- 2026-09-13 08:45 running: resumed after the 00:46 session cut (no training job had launched; only the probe model exists). Training order fixed in PREREG_C2 addendum 2; launching the 8-rung training chain one governor job at a time (est-ram 1.5 GB); PREREG_C3 addendum 2 and IDEA_amber_prior_partner.md written.
- 2026-09-13 08:45 next: C1 closure write-up; B3/C4/C5 code with synthetic tests; agentP_FINDINGS sections; ledger entry for the launch (re-reading the tail first); commit; hourly STATUS.
- 2026-09-13 09:25 running: C2 training chain (p_train_noesm in progress, then conly, pca32, wide, pca32f, pca128, esm8m, raw; one governor job at a time; ~6-10 h). Committed 19ca2ea3: C1 write-up, p_stats/p_b3/p_c4/p_deliver, features cached, IDEA_amber_prior_partner, PREREG_C2/C3 addenda; ledger L26.
- 2026-09-13 09:25 next: turn ends once p_c5 selftest is green and committed; after "PHASE 0 SIGNED OFF": eval shipped then each rung as its checkpoints complete, B3/C4/C5 runs, PROPOSAL_C/B/B_REPLACEMENT, p_deliver for the best rung (ledger + finish turn for PH).
- 2026-09-13 09:26 running: PHASE 0 SIGNED OFF (L33). Training chain: noesm done (5 folds, 29 min, 0.63 GB), conly training. Eval chain launched (s26/p_eval_chain.sh): shipped anchor -> noesm -> B3 run -> conly -> pca32 -> mix -> pairnet -> C4 -> C5 -> wide/pca32f/pca128/esm8m/raw as checkpoints appear; one lane-P job beside the training job.
- 2026-09-13 09:26 next: ledger entry per rung with ST.fmt verbatim (s26/p_rung_report.py); PROPOSAL_B.md after B3; p_deliver for the best rung; PROPOSAL_C.md when PH's C3_RESULT.md exists.

## A (Adversary)
- 2026-09-13 09:15 running: nothing. Done: `s26/EXAMINATION_AUDIT.md` and ledger L31 (reproduction exact on HEAD; one MATERIAL: the 369/356 test count in EXAMINATION D / C35 / L28 against the artefact's 370/357; C26 re-derived from `s13/cache/tors_rows.npz`; C27 in git history at 5fa05cd).
- 2026-09-13 09:15 next: re-check lane E's correction when it lands; then section 2 checks of L27 (DLA, `s26/q_dla.py`) and L22 to L24 (PH censuses); discipline reading (stats_lib docstrings, S25/S16/S21 ledgers, FINDINGS corrections).

## PH (Physics lane)
- 2026-09-13 08:41 running: nothing (gate closed). Done: reading list; PREREG_amber_reject/cis/c3_control; five IDEA files; ph_lib/ph_cis/ph_reject/ph_c3 + synthetic tests; census jobs ph_cis_census (0.038 GB), ph_reject_census (0.117 GB), ph_c3_nativefree; ledger L22-L24; findings and Part IV notes; commit 5dc7a3a6.
- 2026-09-13 08:41 next: READY TO LAUNCH the minute PHASE 0 SIGNED OFF lands: `ph_cis.py floor` (1 min), `ph_c3.py stage1` (2 min), `ph_reject.py cloud` (5 min) then `ph_reject.py chain` (CPU, ~3 h, per-target cells); C3 stage 2 when P delivers (AMBER, probe first). Awaiting tournament ranking for IDEA_branch_select / rotamer_relief part B.
- 2026-09-13 09:35 running: nothing (gate closed). Done since 08:41: PREREG_branch_select.md + ph_branch.py and PREREG_rotamer_relief.md + ph_relief.py (probe modes, per-target cells, 7/7 synthetic tests pass); reject report gains the moved-subset secondary and the empty-draw fallback; findings section 2.3 states the R-first reading; Part IV notes checked against the six prompt items, Walsh numbers and L23 added, speaker-notes table (section 8).
- 2026-09-13 09:35 next: on sign-off: floor, C3 stage 1, reject cloud then chain (as above). If ranked: `jobrun --tag AMBER --name ph_branch_probe --est-ram 1.0 -- python s26/ph_branch.py probe --pdb 1A13` after `ph_branch.py solutions` (CPU, 10 min); `jobrun --tag AMBER --name ph_relief_probe --est-ram 0.8 -- python s26/ph_relief.py probe --pdb 1A13 --m 5`. C3 stage 2 on s26/results/p_best_rung_chains.json when P delivers.
- 2026-09-13 09:21 running: gate open (L33, gate_open() True). ph_cis_floor and ph_c3_stage1 launched 09:20 under jobrun, HELD at the 4/4 job cap (box 73.9% RAM, cpu 75-93%); they start when a slot frees. Then ph_reject cloud, then chain (~3 h, cells).
- 2026-09-13 09:21 next: ledger entries with ST.fmt verbatim per result; C3_RESULT.md after stage 1; answer Adversary checks on L22-L24 in the ledger.

## PR (Presentation lane) -- not yet spawned

## W (Wildcard)
- 2026-09-13 09:15 running: w_selfcopy_retrieval (126 targets, channel A native-free, est 0.5 GB) and w_selfcopy_envelope (126 x 4 leaked fold models, est 0.6 GB) under jobrun. Done: reading list; PREREG_selfcopy_bound; IDEA x4 (selfcopy_proxy_bound, tiebreak_noise_floor, conformational_identity_floor, window_provenance); w_selfcopy.py + synthetic tests (ALL OK); census (0.064 GB) and two 1-target probes (0.104 / 0.298 GB); ledger entry.
- 2026-09-13 09:15 next: commit; agentW_FINDINGS.md skeleton; await the coordinator's ruling on Part B's retrains (PREREG section 8); on sign-off run `w_selfcopy.py endpoint`, `floor`, `report` (seconds each); write the tie-break floor script if the Adversary ranks it.
