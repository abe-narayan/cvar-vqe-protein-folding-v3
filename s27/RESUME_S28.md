# RESUME_S28 (pause notes per lane; append-only)

## D (paused 2026-09-14 23:28)

Lane D, the Adversary and the test suite. Everything of mine is committed on `s26` (last commit
`7008a145` plus the S28-L36 commit `fafbc5bf`; `git status` shows only `s26/governor.log`, not
mine). Nothing of mine is running; the persistent file/ledger watch in my session ends with it.

### 1. Every S28 entry I checked, with its verdict and the open caveats

| entry checked | my entry | verdict | open caveats (still to be answered by the lane) |
|---|---|---|---|
| PREREG_S28_A (amplitude readout) | S28-L1 | STANDS WITH CAVEAT | all five answered by A in S28-L12 (addendum 3) |
| PREREG_S28_B (hopping) | S28-L2 | STANDS WITH CAVEAT | (a) production as comparator: adopted; (b) middle leg: representability fit of the J = 3 ground state still not run (asked again in S28-L26) |
| PREREG_S28_C (FAIL18 detector; readouts) | S28-L3 | STANDS WITH CAVEAT | (a) inter-quantile clause fixed by C's addendum 3; rest answered |
| S28-L6 (C: detector null) | S28-L8 | STANDS | `s28_C_fail18.json` carries no `complete` flag |
| lane B J = 0 anchor | S28-L9 | STANDS (bit-exact, 2 targets x 2 seeds) | none |
| suite + 3 test findings | S28-L10 | recorded | lane A's own poison test is a determinism test (covered by `tests/test_s28_D.py`); `run_chain_target` raises on a NaN native (scorer, not a leak) |
| S28-L8b (B: F5 trainability) | S28-L11 | STANDS WITH CAVEAT | (a)-(e) all accepted and reworded by B in S28-L21 |
| S28-L1b (A: ORACLE ceiling) | S28-L13 | STANDS WITH CAVEAT | quote the like-for-like median contrast (-0.156), not -0.320 |
| S28-L15 (C: readouts, point cloud) | S28-L16 | STANDS (intermediate) | harm magnitudes are Type-M zone: say the sign |
| instrument floor of the built chain | S28-L18, scoped in S28-L27b | recorded; R1, R2 in RETRACTIONS_S28 | tail 0.5 A on 2LNG at 126 targets; two sides of a chain contrast must share a code path |
| S28-L18b (A: recognition, point cloud) | S28-L20 | STANDS (intermediate negative) | quote the matched-budget a500 beside the converged one; "least bad vs random subspace" is Type-M |
| S28-L21 (B: hopping, point cloud) | S28-L22 | STANDS (intermediate negative) | the leading cell (R3 s1 J0.1) is dead: seed-0 twin, PERM control, +0.25 vs production; nine chain cells are a best-of-nine; never quote GS R2's W/L |
| PREREG addenda A2 / B2 | S28-L23 | STANDS WITH CAVEAT | A2: cosine null, e-grid pricing, like-for-like chain control, circuit baseline residual (all answered in S28-L23b); B2: Perron overlap, disconnected graphs (answered in S28-L25/L29) |
| S28-L23b (A2.1: ORACLE cosine) | S28-L24 | STANDS | none |
| S28-L25 (B2: F5-B2 first clause) | S28-L26 | STANDS WITH CAVEAT | representability fit of the J = 3 GS (expressivity vs optimisation) still owed before "optimisation quality" is quoted |
| S28-L26b (A: built-chain verdict, REFUTED) | S28-L27b (+ S28-L28 numbering) | STANDS | lam 1 is Type-M zone with the same sign; no seed-1 owed |
| S28-L29 (B2: second clause fails) and S28-L30 (A2.2 ladder, point cloud) | S28-L31 | STANDS WITH CAVEAT / STANDS (intermediate) | B2's "within 30x" rule is a hypothesis fitted to 2 graphs x 3 J; the A2 random best-of-8 is NOT A SIGNAL (priced) |
| S28-L32 (C: Part 2 chain verdict, closed) | S28-L33 | STANDS | none |
| PREREG_S28_C2 (recognition audit) | S28-L34 | STANDS WITH CAVEAT | ties at 0.5 and geometry beside preferences (adopted in S28-L35); chain multiplicity over all 31 scorers |
| S28-L35 (C2: CAGEO candidate) | S28-L36 | **VETOED** (CAGEO); closure claim STANDS; CONTACT marginal (Type-M) | C2 accepted the veto in S28-L37 and reproduced the pool-member control; it is carried into the chain entry for all 31 scorers |

Retractions index: `s27/RETRACTIONS_S28.md` R1 (scope of C's 1e-5 floor), R2 (scope of my own
S28-L18, 22 -> 126 targets), R3 (C2's "falsified at the registered bar by CAGEO", vetoed).

### 2. Entries NOT yet checked or pending

- Lane B's built-chain verdict (job `s28B_chain`, `s27/results/s28_B_chain_rows.jsonl` at
  1,764 of 2,268 rows at 23:27; 18 arms). To check: production as the comparator (S28-L2(a)),
  the nine 0.7x cells priced as a best-of-nine (`ST.best_of_k_within` over the nine chain
  columns), FAIL18 / 108, the shared-code-path floor (S28-L18/L27b), the three-way split with
  the representability fit still owed (S28-L26), never quote GS R2's 79W/47L as evidence.
- Lane C2's built-chain entry (job `s28C2_chain2`, `s28_C2_chain_rows.jsonl` at 70 of 126):
  expect CAGEO's preference to collapse on ideal-geometry chains (this control's prediction,
  S28-L36 note a); the pool-member control for all 31 scorers (C2 promised it in S28-L37);
  the multiplicity null over all 31.
- Lane A2's chain entry (job `s28A2_chain`, `s28_A2_chain_rows.jsonl` not yet started): the
  step ladder on the built chain, production re-projected in the same job, the random control
  as the mean of the SAME two projected draws (S28-L23(c)), e-grid priced.
- Lane B2's endpoint (k = 10, J = 3 only; gated on the B chain verdict and my check of it).
- S28-L37 (C2's acceptance of the veto): read, consistent with S28-L36; no check entry needed
  unless the chain entry re-opens CAGEO.

### 3. Suite state

- Light non-AMBER files (`tests/test_cvar.py test_data test_energy test_equivalence
  test_geometry test_instrument test_project test_quantum`): 286 pass / 3 skip / 0 fail
  (`s28D_pytest_light`, peak RSS 0.909 GB; S28-L4).
- S28 lane files (`tests/test_s28_A.py A2 B B2 C C2 D`): 70 pass / 0 fail
  (`s28D_pytest_lanes_v3`, 15 s, peak 0.329 GB). Total green gate 356 pass / 3 skip / 0 fail.
- DEFERRED by the coordinator (S28-L5): `tests/test_pipeline.py` (forks 2 workers, tree 1.46
  GB; killed twice by the governor at 95.7% / 98.2% RAM; its kill took lane B's train job as
  collateral), `tests/test_integration.py`, and the two AMBER files (`tests/test_amber.py`,
  `test_amber_frame_invariance.py`, one TEST job each, one at a time). Quiet window: "after
  the B chain lands" (coordinator status 22:47). Not re-queued; nothing launched.
- `python s26/examine.py` not yet run (planned at the close).

### 4. Reproductions of S27 from its artefacts (one per hour, seed stated)

Seeds 101 (pool row 2L7T / LEG_steric), 102 (chain row 2MP9 / DIS_MEAN), 103 (vqe row 3BTB /
DIS+LEG_steric, m 79), 104 (pool row 6MBM / DIS+0.25*LEG_steric), 105 (chain row 2EFZ / DIS):
all exact to the last digit (`s27/results/s28_D_reproduce_{pool,chain,vqe}_seed10{1..5}.json`).
Next: seed 106, kind `vqe`.

### 5. The next three steps, as commands (run from the repository root, after the coordinator lifts the pause)

1. Check lane B's chain verdict when it posts (replace `<arm>` names by the entry's; production
   comparator; then the best-of-nine pricing over the nine 0.7x chain cells named in the entry):
   `C:/Users/abena/miniforge_3/python.exe s27/s28_D_attack.py --rows s27/results/s28_B_chain_rows.jsonl --arm-expr "source|seed|graph|J|readout" --value-col rmsd_chain --basis chain --family "vqe\|" --out s27/results/s28_D_attack_B_chain.json`
   (if the rows use a different arm column, `--arm-col arm`; the kit prints `ST.fmt`, the
   FAIL18 / 108 strata and the family's `best_of_k_within`); then append
   `## S28-L<next> -- ADVERSARY CHECK OF <B's chain entry>` to `s27/LEDGER.md` after re-reading
   the tail, and commit with the trailer lines.
2. Hour-6 reproduction and the C2 / A2 chain checks as they post:
   `C:/Users/abena/miniforge_3/python.exe s26/jobrun.py --agent S28D --tag CPU --name s28D_reproduce_seed106 --est-ram 0.35 -- C:/Users/abena/miniforge_3/python.exe s27/s28_D_reproduce.py --seed 106 --kind vqe`
   then, for C2's chain entry, re-run the pool-member control on the chain rows once C2 names
   its file (`s27/s28_D_c2_poolmember.py` reads `s28_C2_ca_rows.jsonl`; point `ROWS` at the
   chain rows file and the chain channels in `s27/cache`), and for A2's chain entry
   `s27/s28_D_attack.py --rows s27/results/s28_A2_chain_rows.jsonl ...` with production from
   the same job.
3. In the coordinator's quiet window, the deferred files one at a time as TEST jobs, then
   examine:
   `C:/Users/abena/miniforge_3/python.exe s26/jobrun.py --agent S28D --tag TEST --name s28D_pytest_pipeline_q --est-ram 1.6 -- C:/Users/abena/miniforge_3/python.exe -m pytest tests/test_pipeline.py -q -p no:cacheprovider -ra`
   then the same for `tests/test_integration.py` (`--est-ram 1.2`), then `tests/test_amber.py`
   and `tests/test_amber_frame_invariance.py` with `--tag TEST --est-ram 1.0` (one at a time,
   no other AMBER job live), then `C:/Users/abena/miniforge_3/python.exe s26/examine.py`;
   record the counts in a `## S28-L<next> -- SUITE STATUS (close)` entry and finish
   `s27/s28_D_FINDINGS.md` (qualifier table rows for the B chain, C2 chain, A2 chain entries).
## B (paused 2026-09-14 23:30)

Lane B (Gaussian hopping Hamiltonian H = diag(E) - J A) and B2 (kNN spread-spectrum graph).
Branch `s26`; last commit at the pause `afe1f808`. Prereg `s27/PREREG_S28_B.md` (base +
addendum 0; addendum 1 = B2; addendum 2 = B2 endpoint scope). Findings draft
`s27/s28_B_FINDINGS.md` (section 2.2, the built-chain verdict, is NOT written).

### My ledger entries
S28-L8b (F5 trainability vs J, rank-one mechanism, departure diagnostics partial), S28-L21
(point cloud 126/126, intermediate, departure table at n = 126, S28-L11 rewording accepted),
S28-L25 (B2 kNN trainability; three-way split middle leg; m-ladder identity), S28-L29 (B2
share table, second clause failed, kNN rank-one decomposition, connectivity). Lane D's checks
of them: S28-L2 (prereg), S28-L9 (J = 0 anchor bit-for-bit), S28-L11 (F5), S28-L22 (point
cloud), S28-L23 (B2 prereg). No lane D caveat is unanswered: S28-L11 (a) to (e) reworded in
S28-L21; S28-L22's three (production comparator; best-of-nine pricing; no W/L as evidence)
are carried into the chain entry; S28-L23 (a) and (b) answered in S28-L29; (c), (d) are the
gate below.

### Artefacts and completeness (all committed at afe1f808)
- `s27/results/s28_B_rows.jsonl`: 4,914 / 4,914 rows (126 targets x 39 arms), COMPLETE.
- `s27/results/s28_B_summary.json`: point-cloud statistics, COMPLETE for the point cloud; its
  `chain` section is absent until re-run after the chain lands.
- `s27/results/s28_B_chain_rows.jsonl`: PARTIAL, 1,782 rows = 99 of 126 targets x 18 arms
  (`s28_B_chain_arms.txt` lists the 18) at 23:28; job `s28B_chain` (pid in
  `s26/jobs/s28B_chain.json`) is RUNNING and checkpoints per (arm, target); the governor
  re-queues it on a kill (`s26/governor.log`); if it is ever absent from `s26/jobs/` and
  `s26/jobs_done/s28B_chain.json` has exit != 0, relaunch with the same command and it
  resumes: `python s26/jobrun.py --agent S28B --tag CPU --name s28B_chain_r2 --est-ram 0.4 --
  python s27/s28_B_hop.py --chain --arms "$(cat s27/results/s28_B_chain_arms.txt)"`.
- `s27/results/s28_B_split.json` + `s28_B_split_rows.jsonl`: 1,638 rows = 126 x 13, COMPLETE.
- `s27/results/s28_B_mladder.json`: COMPLETE (eps <= 1.1e-13 on every arm).
- `s27/results/s28_B_train.json` + `s28_B_train_rows.jsonl` (504 rows, 12 targets): COMPLETE.
- `s27/results/s28_B_rank1.json`: COMPLETE.
- `s27/results/s28_B2_train.json` + `s28_B2_train_rows.jsonl` (720 rows): COMPLETE.
- `s27/results/s28_B2_share.json` + `s28_B2_share_rows.jsonl` (144 rows): COMPLETE.
- `s27/results/s28_B2_rank1.json` + `s28_B2_rank1_rows.jsonl` (144 rows): COMPLETE.
- Job records: `s26/jobs_done/s28B_{probe_graph,probe_1A13,profile,run,run2,run3,train,
  train2,rank1,split,split2,mladder,chain}.json`, `s28B2_{train,share,rank1}.json`.

### When s28B_chain lands (next three steps, as commands)
1. `python s27/s28_B_analyse.py` (no flags; it reads `s28_B_rows.jsonl`, `s28_B_chain_rows.jsonl`
   and, if present, `s28_B_chain_rows_p1.jsonl`, plus S27's `chain_rows.jsonl :: DIS` for the
   production chain; writes `s28_B_summary.json` with a `chain` section: per arm `mean_chain`,
   `fail18_chain`, `other_chain`, `vs_J0_R1` (paired vs `vqe|s0|NONE|J0|R1` or the seed-1 twin),
   `vs_J0_R1_nonfail18`, `vs_J0_R1_fail18`, `vs_production`; every `fmt` block verbatim).
   Then price the nine 0.7x cells as a best-of-nine: `ST.best_of_k_within` on the (126 x 9)
   matrix of chain differences vs J = 0 R1 (S28-L22).
2. Write ledger entry `## S28-L<n> -- VERDICT ON THE BUILT CHAIN ... (date, B)`: per arm
   `ST.fmt` verbatim from `summary.json :: chain[arm].vs_J0_R1.fmt` and `.vs_production.fmt`;
   answer S28-L2 (a) production as the comparator, (b) the three-way split (Hamiltonian
   quality = the GS rows; optimisation quality = `s28_B_split.json`, S28-L25; emitted
   structure = the chain), (c) sign coherence beside every hopping value (S28-L21 table),
   (d) GS-R3 at J = 0 degenerate wherever it appears, (e) `best_of_k_within` on the J grid and
   on the nine cells; contract addendum 1 (built chain only; FAIL18 vs 108; what is new vs
   S26/S27: cite S27 L6 (VQE tracks the classical prefix on 80 diagonal energies) and S25
   L17); the registered prior (WORSE or null). Then complete `s27/s28_B_FINDINGS.md` section
   2.2 and remove its DRAFT header. Commit both.
3. ONLY after lane D posts its check of that entry: the scoped B2 endpoint (prereg addendum
   2: k = 10, J = 3, both seeds, VQE R1/R3, GS R3, PERM, point cloud; chain only at 0.7x MDE):
   `python s26/jobrun.py --agent S28B --tag CPU --name s28B2_run --est-ram 0.4 -- python
   s27/s28_B2_knn.py --run --k 10` (NOTE: `run_main` currently runs the full J_GRID_B2
   {0, 0.3, 1, 3} and all three graphs; before launching, restrict it to J in {0, 3} and
   graphs (REAL, PERM) per addendum 2 by setting `B.J_GRID = (0.0, 3.0)` and `B.GRAPHS =
   ("REAL", "PERM")` inside `run_main`, add a test, commit), then a B2 endpoint entry `(date,
   B2)` with production as the comparator.

### Gate and scope
The B2 endpoint runs only after the S28B built-chain verdict entry is posted AND lane D has
checked it (S28-L23(d), contract addendum 1 item 15). Scope is addendum 2 exactly; no k or J
beyond k = 10, J = 3. Nothing else is queued or planned for this lane.

## C (paused 2026-09-14 23:30; appended after a write race with lane D's section, whose content stands)

### Ledger entries that are mine
S28-L6 (Part 1, FAIL18 detector: refuted at F1; ORACLE switch ceiling), S28-L7 and S28-L7b
(reproduction check, correction of its wall time), S28-L15 (Part 2, point cloud, intermediate),
S28-L32 (Part 2 verdict on the built chain; D: STANDS at S28-L33), S28-L35 (C2 CA-level audit;
its "closure falsified by CAGEO" is withdrawn), S28-L37 (accepts D's S28-L36 veto, pool-member
control reproduced). Pre-registrations: `s27/PREREG_S28_C.md` (base + addenda 1 to 3),
`s27/PREREG_S28_C2.md` (base + addendum 1). Findings: `s27/s28_C_FINDINGS.md` (Parts 1 and 2
complete; C2 sections C2.0 to C2.2b written; C2.3 is a draft header, no chain number read).

### Artefacts and completeness
| path | state |
|---|---|
| `s27/results/s28_C_features.json` | 126/126, complete |
| `s27/results/s28_C_fail18.json` | complete (blocks, singles, singles_max_null, oracle_switch) |
| `s27/results/s28_C_null_pred_*.npy` (5) | per-permutation predictions; only needed for a switched arm, which was not run |
| `s27/results/s28_C_reproduce.json` | 12 values, max deviation 0.0 |
| `s27/results/s28_C_per_target.md` | the per-target table (FAIL18 first); `python s27/s28_C_fail18.py table` regenerates it |
| `s27/results/s28_C_readout_cloud_rows.jsonl` / `_summary.json` | 29 arms x 126, complete |
| `s27/results/s28_C_readout_chain_rows.jsonl` / `_summary.json` | 12 arms x 126 (+1 probe row), complete, no duplicates |
| `s27/results/s28_C2_ca_rows.jsonl` / `s28_C2_ca_summary.json` | 126/126 (control draws 0-3, with geometry), complete; the summary has ties-at-0.5, geometry, circ_s0 contrasts, head-to-head, the linear combination, the pool-member control |
| `s27/results/s28_C2_ca_rows_seed2.jsonl` / `s28_C2_ca_seed2_summary.json` | 126/126 (control draws 4-7), complete (its summary predates the pool-member-control code; `python s27/s28_C2_recog_audit.py analyse_ca --tag seed2` regenerates it with that control) |
| `s27/results/s28_C2_chain_rows.jsonl` | LIVE, written by job `s28C2_chain2` (70/126 at 23:27); not committed while the job runs |
| `s27/results/s28_C2_chain_rows_checkpoint_70.jsonl` | committed copy of the first 70 rows |
| `s27/results/s28_C2_chain_rows_prepatch_partial11.jsonl` | 11 rows from the first chain launch (module without the geometry fields), superseded, kept for the record |
| `tests/test_s28_C.py` (11), `tests/test_s28_C2.py` (6) | all pass |

### C2 CA-audit result state
POSTED: S28-L35 (seed 0 = draws 0-3; the table is from `s28_C2_ca_summary.json` after the
ties / geometry / circ_s0 patch; seed-2 contrasts beside it from `s28_C2_ca_seed2_summary.json`).
ANSWERED: S28-L36 (D's pool-member veto) in S28-L37, the control reproduced to the third decimal
from `s28_C2_ca_summary.json :: pool_member_control`. Standing reading: the closure claim ("no
scorer in the S27 library recognises the ORACLE structures") stands; CAGEO vetoed; CONTACT /
CONTACT_LL marginal (Type-M on the pool-member control, which was not registered); the learned
linear combination is anti-production (prefers random signed combinations on 0.89; head-to-head
ORACLE vs random 0.53).

### s28C2_chain2
Job `s28C2_chain2` (jobrun, agent S28C, tag CPU, est 0.6 GB; the CA jobs peaked at 0.35 GB),
8 projections per target, per-target checkpoint in `s27/results/s28_C2_chain_rows.jsonl`,
resumable. If it is killed, relaunch with
`python s26/jobrun.py --agent S28C --tag CPU --name s28C2_chain3 --est-ram 0.6 -- python s27/s28_C2_recog_audit.py chain`
(it skips the pdbs already in the rows file). When 126 rows exist:
`python s27/s28_C2_recog_audit.py analyse_chain`
writes `s27/results/s28_C2_chain_summary.json` and prints the CA-level tables for 31 scorers
(the 16 backbone scorers on the projected chains + the 15 CA scorers re-evaluated on them,
suffix `@chain`; one max-over-31 sign-flip null; the pool-member control for every scorer with
a pool channel, the chain ones against the pool members' real-torsion values in `s27/cache`;
geometry of every projected structure; ties at 0.5; circ_s0 contrasts). It feeds the C2
built-chain ledger entry (next number after the ledger tail; "C2 RECOGNITION AUDIT, BUILT
CHAIN ...") and section C2.3 of `s27/s28_C_FINDINGS.md`. Registered expectation (S28-L36(a)):
CAGEO's preference collapses on the projected chains.

### Lane D's S28-L34 caveats and their state
(a) ties at 0.5 with tie counts: DONE (S28-L35 table carries tie counts).
(b) geometry beside every preference: DONE at the CA level (bond / Rg of every ladder point and
    control in S28-L35); the chain rows carry `geom_chain` / `geom_cloud` and `analyse_chain`
    prints them.
(c) closure decided on circ_s0 as well: DONE (circ_s0 contrasts and clause check in the summary
    and the entry).
(d) FAIL18 stratum as k of 18: DONE (`fail18_k`; the entry's table).
(e) chain multiplicity over all 31 scorers: IMPLEMENTED in `analyse_chain` (the `@chain`
    merge), not yet run (waits for the job).
S28-L36(b) (pool-member control on the chain for all 31): IMPLEMENTED in
`pool_member_control`, runs inside `analyse_chain`, not yet run.

### Next three steps, as commands
1. `python s27/s28_C2_recog_audit.py analyse_chain` (after `wc -l s27/results/s28_C2_chain_rows.jsonl` reads 126; if the job died, relaunch with the jobrun line above first).
2. Post the built-chain C2 entry: `grep -n '^## S28-L' s27/LEDGER.md | tail -1` for the number; compose from `s27/results/s28_C2_chain_summary.json` (pref table for 31 scorers, pool-member control, max-over-31 null, the linear combination and its controls, geometry, FAIL18 k/18; `ST.fmt` verbatim for any contrast that clears a clause); append; `git add s27/LEDGER.md s27/results/s28_C2_chain_rows.jsonl s27/results/s28_C2_chain_summary.json && git commit` with the two trailer lines.
3. Fill `s27/s28_C_FINDINGS.md` section C2.3 from the same summary (the closure verdict on the built chain; what damaged expectations; what was not done: AMB deferred, the seed-2 pool-member control), update `s27/STATUS.md` under `## C`, and commit.

## note (lane D, 2026-09-14 23:36): a write race on this file
Lane C's `## C` section (commit 139198a3's message names it) is NOT in the file at that commit
(that commit carries only my `## D` section: my Write landed between C's write and C's `git add`),
and my subsequent restore duplicated `## D` (removed here). Lane C's section must be
re-appended by lane C from its own notes when the pause lifts. Nothing else was changed.
Correction (lane D, 23:38): lane C's `## C` section IS present above (lane C re-appended it at
23:30 after the race); the note's "must be re-appended" is void. The file now carries D, B, C.

## A (paused 2026-09-14 23:35)

Lane A (the amplitude readout) and A2 (the objective's local behaviour at the production point).
Everything of mine is committed on `s26` (last commit "s28 A: pause save"); the only file that
keeps changing is `s27/results/s28_A_chain_rows.jsonl`, appended per target by the running job.

### Ledger entries (mine) and their state
- S28-L1b (headed "S28-L1 -- ORACLE EXPRESSIVITY ...", line ~613): the ORACLE point-cloud ceiling, DONE; lane D S28-L13 STANDS WITH CAVEAT (median contrast), accepted in L26b and the findings.
- S28-L12: numbering correction + answers to D's S28-L1 caveats (a) to (e), DONE.
- S28-L17: PARTIAL notes (superseded by L26b and L18b), DONE.
- S28-L18b (headed "S28-L18 -- RECOGNITION ON THE POINT CLOUD ..."): point cloud 126/126, intermediate, DONE; D S28-L20 STANDS.
- S28-L19: numbering correction, DONE.
- S28-L23b: A2.1 ORACLE cosine diagnostic 126/126 + answers to D's S28-L23 caveats (a) to (e), DONE.
- S28-L26b (headed "S28-L26 -- VERDICT ON THE BUILT CHAIN ..."): THE VERDICT, F1 silent, every lam arm WORSE, DONE; D S28-L27b STANDS.
- S28-L27: numbering correction, DONE.
- S28-L30: A2.2 step ladder on the point cloud (intermediate) + the ORACLE objective diagnostic, DONE.
- NOT YET WRITTEN: the A2.2 built-chain verdict entry (feeds from the running job, below).

### Artefacts and completeness
- `s27/results/s28_A_oracle_rows.jsonl` 126/126 (ORACLE ceilings, production frame check).
- `s27/results/s28_A_recog_rows.jsonl` 126/126 (every recognition arm, point cloud, F parts).
- `s27/results/s28_A_chain_rows.jsonl`: primary arms (prod, oracle_circ, circ_l0/0.3/1/3_i80) 126/126; oracle_aff500 22/126 (dropped by decision, S28-L17); A2 arms (step_e*, rand0/1_e*, circP, circ_e*) 90/126 at the pause and growing (the file holds one merged row per target; the LAST row per pdb is the complete one).
- `s27/results/s28_A2_cosine_rows.jsonl` 126/126; `s27/results/s28_A2_ladder_rows.jsonl` 126/126.
- `s27/results/s28_A_summary.json`, `s28_A2_summary.json`, `s28_A_objdiag.json` (complete; the A2 summary's chain block is partial until the job lands).
- `s27/results/s28_A_structs/<pdb>{,_recog,_a2}.npz` 379 files, committed (whitelisted in `.gitignore`, 12 MB): the emitted clouds the chain phase projects.
- `s27/s28_A_FINDINGS.md`: DRAFT, complete except section 4.3.
- Code: `s27/s28_A_amp.py`, `s27/s28_A_analyse.py`, `s27/s28_A_objdiag.py`, `s27/s28_A2_local.py`, `s27/s28_A2_analyse.py`; tests `tests/test_s28_A.py` (23), `tests/test_s28_A2.py` (5), all green at the pause.

### Running / queued
- `s28A2_chain` (jobrun, agent S28A, tag CPU, est 0.5 GB, peak so far 0.31 GB): projects the 13 A2 arms per target; 90/126 at the pause; per-target checkpoint into `s27/results/s28_A_chain_rows.jsonl`; leave it running. If it is killed, relaunch from the checkpoint with EXACTLY:
  `python s26/jobrun.py --agent S28A --tag CPU --name s28A2_chain_r2 --est-ram 0.5 -- python s27/s28_A_amp.py chain --arms step_e0.1,rand0_e0.1,rand1_e0.1,step_e0.3,rand0_e0.3,rand1_e0.3,step_e1,rand0_e1,rand1_e1,circP,circ_e0.1,circ_e0.3,circ_e1`
- When it lands (`s26/jobs_done/s28A2_chain.json`, exit 0, 126 targets carrying `step_e1`): run
  `python s27/s28_A2_analyse.py` (prints the built-chain block "A2.2 BUILT CHAIN" with `ST.fmt` for every A2 arm vs production, the FAIL18/108 split, and the step-vs-random contrast paired on the SAME two projected draws; writes `s27/results/s28_A2_summary.json`), then append the ledger entry `## S28-L<next> -- A2.2 THE STEP LADDER ON THE BUILT CHAIN (verdict) (date, A2)` with those blocks verbatim, the falsifier decision (some e beats production beyond MDE, fold CI, 5/5 AND beats the random mean beyond MDE; prior: does not fire), `ST.best_of_k_within` over the three e for the step arm and for the circuit arm (D's S28-L23 caveat (b)), and the circP residual beside its own RMSD (caveat (e)); read the ledger tail and append in ONE python process (number collisions happened three times).

### Next three steps to resume
1. `python s27/s28_A2_analyse.py` (after the job lands; if not landed, `wc -l s27/results/s28_A_chain_rows.jsonl` and `tail -1 s26/logs/s28A2_chain.log` first).
2. Append the A2.2 chain verdict entry (numbered from the tail in the same process), then `git add s27/LEDGER.md s27/results/s28_A2_summary.json s27/results/s28_A_chain_rows.jsonl && git commit` with the trailer lines.
3. Fill section 4.3 of `s27/s28_A_FINDINGS.md` from that entry, remove the DRAFT header, update `s27/STATUS.md` under `## A`, and commit. Then `python -m pytest tests/test_s28_A.py tests/test_s28_A2.py -q -p no:cacheprovider` to confirm 28 pass.

### Open lane D caveats
None unanswered: S28-L1 (a) to (e) answered in S28-L12 and addendum 3; S28-L13 accepted in L26b and the findings; S28-L20 (a), (b) answered in L26b; S28-L23 (a) to (e) answered in L23b, with (b), (c), (e) to be honoured in the A2.2 chain entry when written; S28-L27b STANDS with one wording request ("worse at every lam; the size at lam 1 is in the Type-M zone"), adopted in the findings.
