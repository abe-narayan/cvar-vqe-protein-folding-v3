# s28_D_FINDINGS (Sprint 28, lane D: the Adversary and the test suite)

S12 to S25 format: tiers DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN, every
number with its artefact path, then "what damaged my own expectations" and "what I did not do
and why". Companion: `s27/RETRACTIONS_S28.md` (append-only). Scripts `s27/s28_D_*.py`, results
`s27/results/s28_D_*.json`. Ledger entries `## S28-L<n> -- ADVERSARY CHECK OF S28-L<m>` in
`s27/LEDGER.md`; hourly lines under `## D` in `s27/STATUS.md`.

Method (inherited from `s26/agentA_FINDINGS.md`): every positive is attacked within the hour of
its ledger entry, every finalist whether or not positive; the checklist is leakage (grep for
`nat_ca`, `oracle_rr`, `rr`, `native`, `rmsd` in the lane's script, every hit traced; the
NaN-poison test run by me), tie-breaking by array order, iid vs fold-clustered CI disagreement,
concentration in the top-10 targets, k_eff for any minimum over K variants or a grid, regression
to the mean on any tuned parameter, the control's match to the operator's space (set size,
magnitude, basis), a second seed and reversed fold order for anything positive, power on every
null. Verdicts: STANDS / STANDS WITH CAVEAT / VETOED.

## DEMONSTRATED

D1. **Lane B's J = 0 loop is S27's VQE arm bit-for-bit.** `s27/results/s28_D_anchor_B.json`
(job `s28D_anchor_B`, 10 s, peak RSS 0.317 GB): on 1A13 and 2N9M at seeds 0 and 1 the
512-vector p is identical (`==`, max |diff| 0.0), the tail identical (m 72/67/77/71), the
point-cloud RMSD equal to `s27/results/vqe_rows.jsonl :: DIS` to the last digit (abs diff 0.0
on all four). Ledger S28-L9.

D2. **S27 reproduces from its artefact at random.** Hour 1, seed 101 (`s27/s28_D_reproduce.py
--seed 101 --kind pool`): row 7935 of 25578, 2L7T / LEG_steric, stored 2.9099985141915305,
recomputed 2.9099985141915305, abs diff 0.0; the DIS top-75 anchor of the same target 3.3032989
exact (`s27/results/s28_D_reproduce_pool_seed101.json`).

D3. **Lane A's deployable path is native-free, end to end.** `tests/test_s28_D.py ::
test_lane_A_recognition_phase_is_bit_identical_under_nan_poison`: the whole recognition phase
on a poisoned pool (nat_ca / oracle_rr NaN through a monkeypatched `load_pool`) emits every
structure, objective value and weight diagnostic bit-identical while every ORACLE `rmsd_cloud`
is NaN (`s26/logs/s28D_pytest_D_v3.log`, 5 pass). Lane A's own poison test is a determinism
test (S28-L10 item 2). Lane C's and B's poison tests are genuine (they poison `cand` /
monkeypatch `channels_for`).

D4. **Lane B's F5 numbers recomputed (S28-L11).** From `s27/results/s28_B_train.json`,
`s28_B_rank1.json`: J* = 85.7 at n = 9 (the coupling at which the hopping gradient variance
equals the CVaR term's); hop-only slope -1.844 over n = 4..9 and -1.700 over the padding-free
4..8 (per-qubit steps -2.48, -1.37, -2.26, -0.58, -3.25); the diagonal same-spectrum control
decays at -1.856 over 4..8 (the entry's -1.27 is the n = 9 outlier 7.18e-5); the VQE's hop
share at J = 3 equals its sign coherence (0.331 vs 0.321 on 30 targets), the rank-one
mechanism read back.

## ORACLE DIAGNOSTIC

(none of my own)

## HYPOTHESIS

H1. **Lane A's ORACLE expressivity ceilings are dimension counting.** From
`s27/results/s28_A_oracle_rows.jsonl` (126 rows, before lane A's entry): the ORACLE affine-75
least squares is 0.0000 on every target because 75 weights exceed 3n <= 48 (a rank statement);
the ORACLE random 27-dim affine subspace is 0.000 at n = 9, 10 (3n - 6 <= 24 < 26 dof), 0.19 at
n = 11 and rises to 0.97 to 1.09 at n = 15, 16 (Spearman with n 0.80); the ORACLE circuit
family rises 0.12 -> 0.46 over the same n (Spearman 0.62). What survives as a non-trivial
ORACLE reading is the circuit family beating a random linear family at the same parameter
count by 0.32 A; the rest is the count of free weights against 3n - 6. To be put to lane A when
its entry posts (not pre-empted here).

## REFUTED

(none yet)

## OPEN (my queue)

- Lane A's first result entry (ORACLE ceiling, then the recognition arms): H1 above, the
  lam*S~ vs CVaR scale (S28-L1 caveat a), the denominator failures, the nested lam.
- Lane C's Part 2 built chain when it posts (my point-cloud pass with `s27/s28_D_attack.py`
  found every ranker-informed cell null-to-WORSE; the permuted controls are inside their MDE).
- Lane B's endpoint arms F1 to F4 (production comparator, S28-L2 caveat a; the three-way split).
- Hour-2 reproduction (seed 102, a chain row); hour-3 (seed 103, a vqe row).
- The deferred suite files (pipeline, integration, the two AMBER files) in the coordinator's
  quiet window; `python s26/examine.py` at the end.

## What damaged my own expectations

- I expected the light suite to run as one job as it did in S26 (peak 1.69 GB); on this box at
  84 to 87% baseline the full 10-file job peaked at 2.03 GB and was killed, then
  `test_pipeline.py` alone forked two workers (tree 1.46 GB) and took lane B's train job down as
  collateral (S28-L5). The test suite is not free on a loaded box; it is a job like any other.
- I expected lane B's "-1.27 for the diagonal control" to support its own reading weakly; it
  supports it strongly once the n = 9 outlier is dropped (the two controls decay at the same
  rate on 4..8).

## What I did not do and why

- Did not re-queue `tests/test_pipeline.py`, `test_integration.py` or the AMBER test files
  after the coordinator's 19:33 decision (S28-L5); they wait for the announced quiet window.
- Did not run any endpoint experiment of my own; every number here is a recomputation from a
  lane's artefact, a reproduction of S27, or a synthetic-data test.
- Did not post on lane A's ORACLE ceilings before lane A's own entry (H1 is held here).

## Qualifier table (every number that would reach a report, and what is said beside it)

| claim | say with it | source |
|---|---|---|
| (none yet) | | |
