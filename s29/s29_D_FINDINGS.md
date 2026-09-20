# s29_D_FINDINGS (Sprint 29, lane D: the Adversary, the cost-RMSD meter, the suite)

S12 to S25 format: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN, every number
with its artefact path, then "what damaged my own expectations" and "what I did not do and why",
then the qualifier table for every number that would reach the report. Companion:
`s29/RETRACTIONS_S29.md` (append-only). Scripts `s29/s29_D_*.py`, results
`s29/results/s29_D_*.json`. Ledger entries `## S29-L<n> -- ADVERSARY CHECK OF S29-L<m>` in
`s29/LEDGER.md`; hourly lines under `## D` in `s29/STATUS.md`.

Method (inherited from `s27/s28_D_FINDINGS.md`, itself from `s26/agentA_FINDINGS.md`): every
positive attacked within the hour of its ledger entry, every finalist whether or not positive;
the checklist is leakage (grep for `nat_ca`, `oracle_rr`, `rr`, `native`, `rmsd` in the lane's
script, every hit traced; the NaN-poison test run by me), tie-breaking by array order, iid vs
fold-clustered CI disagreement, concentration, k_eff and split-half transfer for any grid or
minimum, regression to the mean, the control's match to the operator's space, a second seed and
reversed fold order for anything positive, power on every null, and for quantum claims the ten
controls and nine questions of charter section 11. Every check states the sprint's comparison
count (`s29/STATE.md`) and prices a 1x-MDE positive against the max-over-K null. Verdicts:
STANDS / STANDS WITH CAVEAT / VETOED.

## DEMONSTRATED

D1. **The cost-RMSD meter reproduces every S28 anchor through its own code path (S29-L2).**
`s29/s29_D_cost_audit.py`, `tests/test_s29_D.py` (7 pass), results
`s29/results/s29_D_cost_audit_{DIS_ca,DIS_chain-s28rows,DIS_SURR_ca}.json`. Ladder Spearman
-0.4023 (built chain, S28's -0.402) and -0.1818 (CA, S28's -0.182); gradient cosine -0.0339
(S28-L23b's -0.034) with the random-direction reference at 0.140 and Spearman(cos, production
RMSD) -0.372; the native's pool percentile 0.3688 on the surrogate (S28-L30's 0.369, with
production below the native on 99/126 and below the pool's best on 6/126) and 0.3676 under the
shipped lookup; pref(ORACLE circ_best vs PROD) 0.0714 on the chain and 0.2063 at CA level with
the pool-member control at 0.0201 / 0.1265 and contrasts +0.0513 (0.80x MDE) / +0.0799 (0.83x).
ORACLE throughout. The ladder cache regenerates C2's structures and asserts them against lane
A's `per_start[0]` / `per_sub[0]` (< 1e-6) and C2's stored DIS scores (< 1e-5) on 126/126.

D2. **The sign of the ladder rho depends on which rungs are in the ladder, and it flips
(S29-L2; new, not in S28).** For the shipped cost: S28's five rungs (production plus four
near-native structures) give -0.182 CA / -0.402 chain; the charter's six rungs (which include a
random signed combination at 3.88 A and a matched Gaussian at 4.12 A) give +0.260 CA / -0.092
chain; all nine give +0.118 / -0.236. The shipped cost orders the BULK correctly and
anti-orders the near-native half. Consequence: "rho = -0.40" must always name its ladder, a
positive charter rho with a negative S28 rho is the expected failure mode of a new objective
rather than a pass, and a cost proposed to this sprint must beat the shipped cost on the S28
ladder. This is `decoy-bank-not-a-pool-proxy` ("garbage rejection reads as skill") measured on
the objective itself.

D3. **Suite green at the first gate (S29-L5).** 17 light files, 381 tests: 378 passed, 3
skipped (the `VERIFY_SLOW` opt-ins), 0 failed, exit 0, 115.2 s, peak RSS 0.93 GB
(`s26/jobs_done/s29D_pytest_light.json`). `test_pipeline.py`, `test_integration.py` and the two
AMBER files are deferred to the coordinator's quiet window, as in S28-L5.

## ORACLE DIAGNOSTIC

OD1. The four meter numbers are ORACLE by construction (the ladder is scored against the
native, the gradient cosine points at the native, the percentile places the native in the pool,
the preference is for a structure chosen against the native). The meter adds no deployable arm
and tunes nothing; passing it is necessary and never sufficient.

OD2. The shipped cost prefers lane A's own native-free optimum (`circ_l1_i80`, mean 3.385 A,
0.34 A worse than production) to production on 0.786 of targets, and scores it below 0.969 of
real pool traces: S28-L18b's "the objective's minimiser is away from the native" read as a
preference, with the mechanism (the minimiser is a structure no real trace resembles) attached.

## HYPOTHESIS

H1. **Rung 6 of lane O and P3 of lane X are both confounded with CONTRACTION until controlled
(S29-L3(a), S29-L4(c)).** Production is a 22%-contracted trace; the shipped cost prefers it to
87% of real pool members on that alone (S28-L36). Lane O's step displaces production along a
difference of two contracted averages; lane X's R2 - R1 differs in effective set size and hence
in how contracted the emitted structure is. Predicted: both effects, if any, will track a
scale-only control. Registered as a prediction here so it can be wrong.

## REFUTED

(nothing yet)

## OPEN (my queue)

- Lane O's rung-6 numbers, when posted: the three required controls of S29-L3 (random direction
  of matched displacement, scale-only, and the blind-difference cosine null) and F6a's measured
  signed-mean reference.
- Lane X's D1 and P1 to P5, when posted: the scrambled-chimera control, the participation-ratio
  control on P3, power and Type-M beside any GO.
- The heavy test files in the quiet window.
- A reproduction of one S28 number from its artefact every two hours (seeds stated).
- Lane M's harness audit, lane T's derivations and lane L's imports as they land.

## What damaged my own expectations

- I expected the meter's ladder rho to be a property of the cost. It is a property of the cost
  AND the ladder: the same cost is +0.26 on the charter's rungs and -0.18 on S28's at the same
  CA level. I had read S28's "-0.40" as a statement about the objective; it is a statement
  about the objective on the near-native half of a ladder, and the brief's six rungs are a
  weaker test than the five they replace, not a stronger one.
- I expected the shipped table and its linear surrogate to give the same four numbers. They
  agree to 1e-3 on the cosine and the preferences but differ on the percentile in the third
  decimal (0.3676 vs 0.3688) and on the CA ladder rho in the third (−0.1818 vs −0.1857),
  because the table is piecewise constant and creates ties the surrogate breaks. Any anchor
  quoted to three decimals must say which of the two it came from.
- I expected finite differences on a piecewise-linear cost to be safe at h = 0.01 A. They are
  not: the distogram's risk grid has 0.05 A spacing, so h = 0.01 straddles knots and the cosine
  drifts by up to 1e-2. h = 1e-3 holds to 3e-4.

## What I did not do and why

- Did not run the heavy test files (pipeline, integration with VERIFY_SLOW=1, the two AMBER
  files): they wait for the coordinator's quiet window (S28-L5's pattern), and no file outside
  `s29/` and `tests/` has changed on this branch since S28-L45's full-tree run.
- Did not meter the other 30 S27 scorers. The meter takes them (`--f CAGEO`, `--f LEG`, ...)
  but a table of 31 costs nobody proposes is a table, not a check; they are metered on request.
- Did not build the `--basis chain` cache (28 projections per target would be an hour-long job).
  The shipped cost's chain numbers come from C2's stored rows (`--basis chain-s28rows`), which
  is the artefact S28 verdicts were written on; a new cost that needs the real chain basis gets
  the cache built for it in one job, on request.
- Did not post on lane O's or lane X's expected results, only on their preregs.

## Qualifier table (every number that would reach a report, and what is said beside it)

| claim | say with it | source |
|---|---|---|
| the cost-RMSD meter reproduces S28's four anchors (rho -0.402 / -0.182, cos -0.034, pctile 0.369, pref 0.071 / 0.206) | every number ORACLE; the rho must name its ladder AND its basis (chain vs CA), the percentile must name the surrogate S~ (0.3688) or the shipped table (0.3676), and the cosine is the SURROGATE's gradient because the shipped table is piecewise constant (as in S28-L23b) | S29-L2, `s29/results/s29_D_cost_audit_*.json` |
| the shipped cost's ladder rho is +0.26 on the charter's rungs and -0.18 on S28's at the same CA level | the charter's ladder contains two garbage rungs (RAND_SIGNED 3.88 A, GAUSS_MATCHED 4.12 A) that any ordering cost rejects; garbage rejection reads as skill (`decoy-bank-not-a-pool-proxy`); the S28 ladder is the binding one and a new cost must beat the shipped cost there | S29-L2 |
| the meter refuses a cost that reads the native | the poison is the CONTEXT (`ctx.nat_ca`, `ctx.oracle_rr` NaN), so a cost that reads the native through another path (its own `load_univ`, a cached `oracle_rr`) is NOT caught: the leakage grep of the lane's own source stays part of the check | `tests/test_s29_D.py`, S29-L2 |
| lane O's prereg: accepted for rungs 1 to 5 and 7, held for rung 6 | the three controls of S29-L3 must be in the same entry as any rung-6 positive; F6a's second clause compares a mean over 126 against a single-draw magnitude and cannot fire as written | S29-L3 |
| lane X's prereg: accepted, probe may run | D1's ORACLE best is a minimum over up to 262,144 structures against a minimum over 8 and needs the scrambled-chimera control; P3 needs a participation-ratio-matched random-weight arm; any GO on 12 targets carries its power and Type-M and the words "GO, not a result" | S29-L4 |
| the suite is green (378 / 3 / 0 on 17 light files) | the heavy files and the two AMBER files are NOT in this count and wait for the quiet window; the carry-over from S28-L45 covers only files neither lane has touched | S29-L5 |
