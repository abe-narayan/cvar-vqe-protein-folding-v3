# PREREG_D — WORKSTREAM D, AUDITOR/ADVERSARY/STATISTICIAN (Sprint 23)

Never edited after data; addenda are dated and appended. Filed 2026-09-08, after reading
`s23/BRIEF.md`, `s22/LEDGER.md`, `s22/CLAIMS.md`, `s21/LEDGER.md` in full, and after discovering
that Workstream A had already landed `s23/gscale.py` / `s23/results/gscale.json` (n=126, complete)
before this file was written — that run is the object of Duty 1/2 below, not a hypothetical.

**Instrument.** `s12/instrument.py`, 126 cluster-disjoint dev targets, 5 pinned folds, K=500 pool,
shipped Bayes-risk distogram score, top-75 coordinate average (incumbent, point cloud, 3.0483 Å —
matches BRIEF's 3.048). Point-cloud basis throughout, matching Workstream A.

## PLAN

**D1 (audits BRIEF duty 1).** Independently reproduce the n=126 scoping result (already run by A
as `gscale.py`); audit whether the per-target continuous oracle's −0.158/−0.137 Å is a selection
artefact, and derive (not assert) the correct null for a continuous one-parameter oracle, contrasted
with the project's existing min-of-K framework (S22 M3/M4). Run the split-half transfer test BRIEF
names explicitly, mirroring `s22/mreal.py`'s construction exactly (same seeding convention, same
half-pool-of-250 design, same fold-clustered bootstrap).

**D2 (audits BRIEF duty 2).** Independently verify A's nested-CV global-scale number from the raw
JSON array (not the printed log). Derive the population-optimal global scale directly from
closed-form per-target quadratics (no grid dependency) and report the captured fraction against the
TRUE (not grid-truncated) ceiling.

**D3 (audits BRIEF duty 3).** Read `gscale.py` end to end for leakage, tie-breaking, completion-flag
and basis-discipline defects, per BRIEF's named checklist.

**D4 (BRIEF duty 4, applied throughout).** Every number below carries paired per-target differences,
iid AND fold-clustered bootstrap CIs, SE, MDE = 2.8016×SE reported as a multiple, W/L, and worst-
target degradation where relevant.

**Falsifiers, declared before running anything beyond what A had already produced:**
- D1: if the split-half transfer fails to beat fixed s=1.0 past its own MDE with a fold-clustered CI
  excluding zero, the oracle ceiling is confirmed to be a selection artefact and A's H1 should be
  closed entirely, not just the global-constant sub-arm.
- D2: if the exact (non-grid) population-optimal global scale captures materially more than A's
  measured ~0%, A's own verdict was itself an artefact of grid resolution and should be revisited.

Six-axis operator forks are stated per-script (`d_scale_closedform.py`, `d_scale_transfer.py`,
`d_scale_globalcapture.py`, `d_scale_placebo.py`/`d_scale_placebo2.py`), each before its own RMSD
was read, in the same convention as A's `gscale.py`.

---

## ADDENDUM 2026-09-08 (same session) — ALL FOUR DUTIES COMPLETE; ONE SELF-CAUGHT DEFECT IN MY OWN
## FIRST PLACEBO CONSTRUCTION, WITHDRAWN AND REPLACED

Full results and derivations in `s23/agentD_FINDINGS.md`. Headline dispositions, dated at time of
writing, not retroactively edited:

1. **A defect found in `gscale.py`, load-bearing.** The oracle grid `[0.90,1.10]` (81 points)
   truncates the per-target optimum at the boundary for **63/126 targets (50%)**. The TRUE
   unconstrained closed-form ceiling is **−0.3403 Å**, not the reported **−0.1374 Å** — 2.5× larger.
   **A's own headline oracle number for this sprint is an undercount and should be corrected before
   further use.**
2. **D1 falsifier did NOT fire — the opposite of my registered directional expectation.** The
   split-half transfer test recovers **98–100%** of the in-sample ceiling out of sample (full panel
   −0.320 held-out vs −0.323 in-sample; non-FAIL18 −0.171 vs −0.174), far cleaner than S22's m-ladder
   (65%). **The per-target propensity for a scale correction is real, not a selection artefact, on
   this project's own strongest evidentiary standard (the transfer test that saved a claim last
   sprint).** I hold my own directional hypothesis wrong and say so.
3. **But a self-caught defect in my own first placebo, corrected same-session.**
   `d_scale_placebo.py`'s within-target atom-permutation control was confounded (it degrades the
   base RMSD to ~7.5 Å, a different-difficulty problem) and its "387%" headline is WITHDRAWN, not
   used anywhere below. Replaced by `d_scale_placebo2.py` (cross-target, same-length real natives),
   which is fair-difficulty and shows the true native pairing is **not distinguishable from a random
   same-length native** by this fit (placebo delta −0.53 Å ≥ real delta −0.34 Å) — reframing, not
   reversing, finding 2: the signal that transfers is real but is NOT primarily an identity-match
   effect.
4. **D2 falsifier did NOT fire — A's near-zero global capture is CONFIRMED, not weakened, by using
   the exact (non-grid) closed form.** Best-possible single global scalar captures **1.7%** of the
   true ceiling (0.6% on non-FAIL18); direction of s* splits 53/126 wants-expansion vs 73/126 wants-
   contraction — the corrections partially cancel, which is the mechanism, not just "spread is wide."
   **Do not invest further in a global scalar; the theoretical maximum is already priced.**
5. A per-length-bucket constant (8 buckets, in-sample) captures 7.6% — better than global, still
   small, and the standing finite-sample bound (S22 C1) applies directly at these bucket sizes
   (9–23 per bucket, far below n≈100) — **not pursued further, per the BRIEF's own standing
   constraint.**
6. **Concentration, flagged in the BRIEF's own words.** 18 FAIL18 targets (14% of the panel) carry
   **57%** of the true ceiling's total magnitude. Any global or length-keyed arm is chasing a signal
   overwhelmingly parked in a small, already-known-pathological subset.

No arm is promoted from this file. This prices a ceiling and audits an arm; it does not ship one.
