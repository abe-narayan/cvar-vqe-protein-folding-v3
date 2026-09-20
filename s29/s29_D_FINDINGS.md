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

D4. **The meter's first customer, lane X's pair log-score (S29-L6).** On the binding S28 ladder
it is +0.200 rho above the shipped cost (fold CI [+0.100, +0.343], 5/5 folds, 1.45x MDE, power
0.98), which makes it the first cost in the record that is UNINFORMATIVE (+0.018, CI straddling
zero) rather than ANTI-INFORMATIVE (-0.182, CI below zero) on the ladder from production to the
native. On the charter's rungs the difference is 0.43x MDE and on all nine 0.93x: the gain is
specifically on the near-native half. Recognition is unchanged: the native sits at the 37.8th
percentile of its own pool (shipped 36.8th, difference 0.25x MDE), production scores below the
native on 78/126, and the pool-member contrast (+0.092, 0.88x) is in the same Type-M zone as
the shipped cost's. Its gradient is undefined (a bin lookup; 99.2% zero components, cosine NaN
on 114/126), so the +0.110 printed on the surviving 12 is not a measurement and is not quoted.
`s29/results/s29_D_cost_audit_X_cost_nll_ca.json`.

D5. **Addendum 20 is mechanical, and the shipped cost's descent EXPANDS (S29-L10).** Every
cosine the meter prints now carries its shrink signature: descend 0.3 A along -grad f from
production and report the bond and Rg ratios with the count of contracting targets. Shipped
cost: bond x1.0438, Rg x1.0248, 18/126 contract. Lane X's pair log-score: x1.0232 / x1.0065,
3/12. Neither buys its cosine by shrinking, so "the cost just prefers contracted things" is
measured FALSE at the gradient level for both -- while remaining true of the shipped cost's
RANKING (production beats 87% of real traces, S28-L36). Two defects in my own meter were found
and fixed in the same patch (an all-NaN cosine axis crashed the renderer after the run had
completed; a partially defined cosine was printed as a measurement), with regression tests.

D6. **The within-realism-band measurement: F1 fires, F2 fails (S29-L33).** 126 targets, 500 pool
members, 12 CA scorers, 3 pre-registered realism statistics, 2 arms. Positive within-band
ordering exists on 50 of 70 cells (DIS +0.5075 [+0.462, +0.547] under R1; +0.5397 under R2;
max-over-scorers p_max 0.000), so information orthogonal to realism is present and substantial.
But rho_in > rho_across fails on 58 of 70: DIS is +0.568 globally and +0.508 / +0.540 / +0.361
in-band under R1 / R2 / R3, with the paired fold CI below zero. **Conditioning on realism removes
ordering skill rather than revealing any.** The one sign change is SS_MATCH (-0.009 -> +0.060
[+0.045, +0.077], replicated on two realism definitions and both arms), worth ~1e-3 A by the
square law. Compactness loading, measured: rho(R, Rg) = +0.460 (R1), -0.015 (R2, but a FOLDED
function of Rg), +0.404 (R3). The registered Gaussian partial correlation OVERSTATES what
survives narrow banding (DIS/R1 predicted +0.596, measured +0.508, width curve to +0.418).
`s29/results/s29_D_band_ca.json`.

D7. **Assumption B2 of the achievable bound survives 21 new displacement fields (S29-L35), and
the limiting quantity is the SIGN, not the alignment.** No field's signed mean cosine clears the
0.140 random reference with a 2-SE margin (best CHAN_DISTPOT +0.1128 [+0.088, +0.137]); the
largest ORACLE gain through a global step is 0.019 A. But every field's PER-TARGET |cos| is
0.251 to 0.325, about twice a random draw's, so with a perfect per-target sign the same fields
reach 2.708 A (EXPAND) to 2.894 A (MEDOID) on the point cloud -- below lane T's 2.98 A margin,
and reproducing the bound's structure across an operator class rather than one direction. EXPAND
(pure de-contraction) has the largest unsigned alignment and a sign at 0.468: production is 22%
contracted and expanding it points away from the native as often as toward it.
`s29/results/s29_D_fields.json`.

D8. **The shrink grid: the cosine half of the rule-20 mechanism is refuted, the percentile half
confirmed (S29-L37).** Over s in 1.0 down to 0.3 the cosine goes -0.034 -> -0.056 (s = 0.6) ->
-0.033: it never rises and never crosses zero. The native percentile worsens monotonically,
0.369 -> 0.491, 9 of 9 steps. Over the same grid the S28 ladder rho IMPROVES (-0.186 -> -0.054)
and pref(circ_best) nearly doubles -- a trade-off between meter numbers that a single-number gate
would have been gamed by. Contract addendum 20 stands as a rule; its stated mechanism does not.

D9. **Lane B's non-prefix subset optimum reproduces independently (S29-L38).** All 12 best pairs
and values reproduce from my own single-pass code; the search is exhaustive (124,750 = C(500,2)
on every target); the frame is the deployed readout's own at max |dev| 0.00e+00; and the verdict
survives the EXACT shipped lookup as well as the piecewise-linear surrogate, with 0 sign flips
(two small-gap targets change WHICH pair wins and neither becomes a prefix). 7JGX's +0.0000 gap
is a genuine argmin at the prefix, not a code path. One defect recorded: lane B's tie rule takes
the first tied index in ARRAY order despite its comment, which biases toward the prefix and is
therefore conservative for its own claim; tie sets were size 1 or 2 and it changed nothing.
`s29/results/s29_D_bcheck_pairs.json`.

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

H2. **A proper scoring rule of the same posterior is not adversarial where its Bayes risk is
(S29-L6).** The shipped cost is an expected-L1 Bayes risk under a ~2x over-confident posterior
and its minimiser contracts; the pair log score is a proper scoring rule of the same posterior
and its minimiser does not. Same information, different functional, +0.200 of ladder rho.
Prediction, so it can be wrong: the difference will NOT survive to an endpoint RMSD through a
fixed readout, because ordering the ladder is not recognising the native and the native's
percentile did not move (`better-matrix-worse-ranking`; `averaging-space-beats-the-objective`
prices the whole objective channel at 0.171 A).

H1. **Rung 6 of lane O and P3 of lane X are both confounded with CONTRACTION until controlled
(S29-L3(a), S29-L4(c)).** Production is a 22%-contracted trace; the shipped cost prefers it to
87% of real pool members on that alone (S28-L36). Lane O's step displaces production along a
difference of two contracted averages; lane X's R2 - R1 differs in effective set size and hence
in how contracted the emitted structure is. Predicted: both effects, if any, will track a
scale-only control. Registered as a prediction here so it can be wrong.

## REFUTED

(nothing yet)

## OPEN (my queue)

- Lane O's rung-6 numbers, when posted (its artefact `s29/results/s29_O_lfo.json` already shows
  the leave-fold-out step selecting t = 0 on all five folds, so the arm IS production and the
  three controls of S29-L3 are moot for a null arm -- as S29-L3 said they would be; what still
  needs saying in that entry is that the per-target ORACLE argmin sits at the grid's LEFT edge
  (t = -1, the blind average) on 29 of 126 targets, i.e. the grid truncates the operator in the
  direction OPPOSITE to H1): the three required controls of S29-L3 (random direction
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
- I expected my own meter to be the one thing in the sprint that could not lose a result. It
  lost a 126-target CAGEO run to a print statement (an all-NaN cosine axis, `ci95_fold = None`,
  TypeError in `render` AFTER every target had been computed). The lesson is the S24 stats_lib
  one restated: the axis that is undefined is the axis nobody wrote a test for.
- I expected the shipped cost's descent direction to contract the structure, because its
  ranking prefers contraction (S28-L36) and its minimiser is a contracted average (S25 L2). It
  expands: bond x1.044, Rg x1.025, 18/126 contracting. Ranking preference and gradient direction
  are different statements about the same cost and I had them fused.
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
| the band experiment: F1 fires on 50/70 cells, F2 fails on 58/70 | say BOTH halves: positive within-band ordering exists AND conditioning removes skill rather than revealing it; name the realism definition in the same sentence (a result on one does not transfer, and R1/R3 are compactness-loaded at rho(R,Rg) +0.46/+0.40 while R2 is linearly clean but a FOLDED function of Rg); the design is the incidental-parameter elimination remedy and is SILENT BY CONSTRUCTION about the per-target sign | S29-L33 |
| SS_MATCH goes from -0.009 global to +0.060 in-band | replicated on two realism definitions and both arms, so it is not a single cell; and at that size the square law prices it at ~1e-3 A -- a mechanism result, never an operator | S29-L33 |
| B2 survives: no displacement field's signed mean cosine reaches 0.140 | **state the FAIL18 half in the same breath**: 10 of 21 fields clear a 20,000-draw random-18 null on the FAIL18 stratum where 1 is expected (CONS_TRIM +0.348 vs +0.046, p 0.0001), AND it does not move the bound, because at n = 18 the intervals are wide and even CONS_TRIM's lower bound of 0.128 sits under the 0.140 reference. A critic who reads only the first half has a case; a reader who gets both has the truth. Also say that eight of the 21 have signed means whose fold CI excludes ZERO -- they are small real alignments, not noise | S29-L35 |
| the sign-oracle ceiling is 2.708 A | ORACLE in the sign AND the step; it is the bound's own formula fed the per-target \|cos\|, not a reachable arm; it sits BELOW lane T's 2.98 A margin and is reported as a refinement of the bound's description, not a challenge to it | S29-L35 |
| the shrink grid refutes the cosine half of rule 20's mechanism | the percentile half is confirmed (0.369 -> 0.491, 9/9 steps), so the RULE stands and only its justification is withdrawn; and the same grid shows the ladder rho improving while the percentile degrades, which is why the meter has four numbers | S29-L37 |
| lane B's claim 1 stands | reproduced under a SECOND functional (the exact shipped lookup) as well as the surrogate, because three of the twelve gaps are smaller than the two readings' typical difference; and the tie-rule defect is recorded as conservative-for-B, not as a flaw in the result | S29-L38 |
| lane B's claim 2 (+0.4278 A at 0.59x MDE, n = 12) | read it in BOTH directions: not evidence the lift helps, and at power 0.38 not strong evidence it hurts either (Type-M 1.61 says the magnitude is inflated ~60% if it is noise-significant); n = 35 clears its own MDE and n = 126 costs 16 minutes | S29-L38 |
| the deployed VQE arm vs a fixed target-independent prefix | -0.0097 A at 0.43x MDE on the point cloud (NOT MEASURED); the built chain is the registered primary and decides; T's "equals on 120/126" is false at the structure level (43/126 within the chain floor) and true only as an endpoint statement | S29-L26 |
