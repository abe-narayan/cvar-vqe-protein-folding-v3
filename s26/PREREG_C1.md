# PREREG_C1 -- CLOSURE REPRODUCTION FROM PERSISTED ARTEFACTS (lane P, Sprint 26)

Written 2026-09-13. This is a reproduction protocol over persisted artefacts, not an endpoint
experiment: no structure is built, no model is trained, nothing is selected. It is allowed before
sign-off (lane brief, deliverable 5) and its numbers were read while locating the artefacts.

## Hypothesis

**H_C1.** The two closures that Proposal C rests on are reproducible from the artefacts on disk to
three decimals: (a) the set-transformer over the full signed deviation map has a FLAT learning
curve (S12 agg_FINDINGS section 6b: 3.043 -> 3.026 from n = 8 to n = 75 training targets, against
a leaked label at 2.534 from n = 8), and (b) a perfect ranker inside the shipped top-25 returns
2.609 A (S17 L10/L23, `s17/inband.py`).

## Exact falsifier

If any quoted number is not reproduced from the named artefact to within 0.001 A, the closure
is DOWNGRADED to "cited, not reproduced" in the Proposal C write-up, and the ledger records the
discrepancy. If the artefact is missing, the closure is marked "artefact lost" (as `s7/repr_tune.json`
already is for S7-11).

## Comparison arm, basis

No comparison arm: a reproduction. Basis, as in the originals: S12 decoder numbers are RAW point
clouds (mean weighted coordinate average, `mean_raw`); the S17 in-band numbers are selected
single-candidate CA-RMSDs inside the shipped objective's own top-25 (selection basis). Neither is
the built chain, and neither will be compared to it.

## What was read (2026-09-13), artefact paths, and the result

`s12/results/agg_dec_v2_n8.json` mean_raw 3.0433 · `agg_dec_v2_n16.json` 3.0491 ·
`agg_dec_v2_n32.json` 3.0456 · `agg_dec_v2_n64.json` 3.0358 · `agg_dec_v2.json` 3.0258 ·
`agg_dec_v2_101.json` 3.1471 (fixed 40 epochs, no early stop; confounded, not read as a decline) ·
leaked label: `agg_dec_v2_oracle_n8.json` 2.5342 · `agg_dec_v2_oracle_n32.json` 2.2004 ·
`agg_dec_v2_oracle.json` 2.1460. v1 harness: `agg_dec_weight_n8.json` 3.0561 · `agg_dec_weight.json`
3.0640 · `agg_dec_oracle_n8.json` 2.5335 · `agg_dec_oracle.json` 2.1466. Each file carries
`raw` as a 126-entry per-target dict. All reproduce the S12 table (3.043 -> 3.026; 2.534) exactly.

`s17/results/inband.json`, 126 rows, cells['25']: band_best mean 2.6087 (median 2.5155),
random-in-band 3.4676, distance argmin 3.4540, consensus 3.3692, Legacy 3.4149, band mean 3.5016;
band_best - random: -0.8589, SE 0.0600, MDE 0.1680, 126W/0L. Reproduces S17 L10's table
(2.609 / 3.468 / 3.454 / 3.369 / 3.415) exactly.

## Expected effect vs MDE

Not applicable (no effect is estimated). Recorded for the presenter: the in-band contrast's own
MDE is 0.168 A, and the learning-curve range (3.026 to 3.049) is 0.023 A wide, inside the S12
seed-to-seed spread (seed 0 3.0640 vs seed 1 3.0363 on v1; 3.0258 vs 3.0337 on v2).

## Memory, agent-hours

Reads JSON only; < 100 MB; 0.5 agent-hour including the write-up paragraph and table in
`s26/PROPOSAL_C.md`.

## Rule 0 forks (a reproduction still has them)

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | the artefacts' own stored means (`mean_raw`, `cells['25']`) | recomputing from per-target rows with a different aggregation |
| basis | as in the originals (raw cloud; in-band selection) | translating either to the built chain |
| readout | the stored learning-curve points at n = 8/16/32/64/75 | interpolating a slope |
| normalisation | none | none |
| null | the leaked-label control, as stored | a fresh permutation null |
| label | mean CA-RMSD as stored | median, W/L |
