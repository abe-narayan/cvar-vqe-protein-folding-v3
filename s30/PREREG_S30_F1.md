# PREREG S30-F1 — is the tail pool-limited or selection-limited?

Lane F, Sprint 30. Committed before any aggregate number for this experiment exists.
Code: `s30/s30_F_stagegap.py`. Output: `s30/results/s30_F_stagegap.json`.

## The question

Production emits mean 3.2105 Å built chain (126 targets); the worst 18 (`s12.instrument.FAIL18`,
pinned, not re-derived here) average 6.0200 and the other 108 average 2.7423. Capping the worst
targets is worth more to the mean than improving every target. Before asking *how* to fix the
tail, decide *what is broken* there:

- **pool-limited** — the candidate pool on those targets contains nothing good, so no selection
  or readout change can help, and the fix must be in generation/retrieval; or
- **selection-limited** — the pool contains good members and the pipeline discards them, so the
  fix is in the filter/readout.

These point at completely different work, so this is measured first and nothing else in the lane
starts until it reads out.

## The decomposition (ORACLE; every row below is ORACLE and is labelled so wherever quoted)

The shipped path has three narrowing stages per target. Write `rr` for the ORACLE point-cloud
CA-RMSD of a window to the native (`s8/generate_univ/<pdb>.npz`, field `rr`).

```
U   universe of windows (9.8k-21.5k per target)          ORACLE best = min rr over all
P   BLOSUM pool, order[:500]                             ORACLE best = min rr over pool
T   score-selected top-75, pool[sub] from the            ORACLE best = min rr over top-75
    production record bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json
E   the emitted structure (the m=75 average, chained)    production RMSD
```

Gaps, per target:

```
G_retr   = best(P) - best(U)     what BLOSUM retrieval gives up
G_filt   = best(T) - best(P)     what the score's 500 -> 75 filter gives up
G_read   = prod    - best(T)     what the readout/averaging gives up
```

## Controls, each in the operator's own space (contract rule 6)

- `G_filt` is compared against **best-of-random-75 drawn from the same 500 pool**, 2,000 draws
  per target, common random numbers across targets, seed 30001. Best-of-k is an order statistic
  (`grid-oracles-are-order-statistics`), so the raw drop from 500 to 75 is *expected* and is not
  evidence of anything; only the excess over the random-75 control is.
- `G_retr` is compared against **best-of-random-500 drawn from the same universe**, 2,000 draws
  per target, same seed stream. Same reasoning.
- Skill is reported both as an Ångström excess and as the **per-target percentile** of the
  realised best within its own random distribution (0.5 = no skill, > 0.5 = the stage is worse
  than random, < 0.5 = better than random).

## Falsifiers, registered before the numbers

**F1a — the "pool is not the limit on the tail" claim.**
Fires (claim SURVIVES) only if ORACLE best-of-pool on FAIL18 is **below 3.00 Å**, the cap used in
the sprint-open counterfactual. If ORACLE best-of-pool on FAIL18 is at or above 3.00 Å the tail is
pool-limited at the cap and this lane pivots to retrieval/generation. No null needed: this is a
one-sided factual threshold on an ORACLE quantity.

**F1b — the "the tail is filter-limited, distinctively" claim.**
Registered statistic: `share_filt = G_filt / G_total`, `G_total = prod - best(P)`, computed as a
ratio of means within each stratum (not a mean of ratios). Claim fires only if **all three** hold:

1. `share_filt(FAIL18) >= 2 x share_filt(other 108)`;
2. the FAIL18 `G_filt` mean lies **outside** the central 95% of a **random-18 null** (20,000 draws
   of 18 targets from the 126, seed 30002) — the null S29 ran for exactly this purpose;
3. the effect survives a **difficulty control**: FAIL18 is selected on production RMSD, so
   `G_filt` is regressed on production RMSD over the 108 and the FAIL18 `G_filt` is compared to
   that regression's extrapolated prediction. The claim fires only if the FAIL18 residual is
   positive and outside the 95% prediction band. If `G_filt` is simply a monotone function of
   difficulty, this is a difficulty statement, not a regime statement, and will be reported as one.

**F1c — the "the score's filter has negative skill on the tail" claim.**
Fires only if, on FAIL18, ORACLE best-of-top-75 is **worse** than the mean ORACLE best-of-random-75,
with median per-target percentile > 0.5, and the fold-clustered CI of the paired difference
excluding zero. Prediction registered now, from `sequence-conditioning-hurts-the-failures` (BLOSUM's
deployable channel is rho = +0.001 on FAIL18 and the blind pipeline beats the shipped one there):
**I expect the filter to be better than random on the 108 and no better than random, possibly
worse, on FAIL18.** If instead the filter beats random on both strata, the mechanism in that memory
does not reach the filter stage and I will say so.

## Basis and labelling

`rr` is a **point cloud** quantity. The built-chain figures for `prod`, `best1_top75`,
`best1_top128` and `best1_pool` already exist from S29 lane O
(`s29/results/s29_O_ladder_table.json`) and are quoted from there; every chain/cloud basis is
named in the same sentence (contract rule 1). Nothing here tunes a deployable parameter on the
native (contract rule 9): this measurement produces no deployable object at all, only a diagnosis.

MDE is computed per comparison as 2.8016 x SE from that comparison's own paired differences
(contract rule 3); fold-clustered CIs govern (rule 4). Below 0.7x MDE is NOT MEASURED (rule 2).

## What I will NOT conclude

That a fix exists. Every number here is ORACLE and bounds a prize rather than delivering one. A
selection-limited verdict says the information is present in the pool, not that any native-free
rule can find it — and S29 closed recognition three independent ways. The output of this
experiment is a direction for the rest of the lane, not a result.
