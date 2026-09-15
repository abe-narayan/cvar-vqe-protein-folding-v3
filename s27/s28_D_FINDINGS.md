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

(none yet)

## ORACLE DIAGNOSTIC

(none of my own)

## HYPOTHESIS

(none yet)

## REFUTED

(none yet)

## OPEN (my queue)

- The three lanes' preregs as they land: check each falsifier is falsifiable and each control
  matches the operator's space; one ledger entry per prereg.
- Every positive within the hour; every finalist.
- `tests/test_s28_<lane>.py` for any operator a lane leaves untested.
- One S27 reproduction per hour (seed stated).
- `python s26/examine.py` once at the end.

## What damaged my own expectations

(none yet)

## What I did not do and why

(none yet)

## Qualifier table (every number that would reach a report, and what is said beside it)

| claim | say with it | source |
|---|---|---|
| (none yet) | | |
