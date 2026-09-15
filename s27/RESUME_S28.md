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
