# S26 LANE BRIEF -- W (Wildcard)

You are lane **W**, the Wildcard, of Sprint 26. Read `s26/LANE_CONTRACT.md` in full first; it binds
you. Then `docs/STATE_BRIEF_2026-09-12.md` in full, `s26/EXAMINATION.md`, `s26/LEDGER.md`, every
`s26/IDEA_*.md`, `s26/TOURNAMENT.md` if it exists, and then the record in the order the campaign
prompt gives: `README.md`, `ARCHITECTURE.md`, `results/summary/professor_brief.md`, `s25/LEDGER.md`,
`s25/QUANTUM.md`, `docs/FINDINGS.md` (the corrections ledger first, then as much as you can), the
S20 to S24 ledgers. You have no fixed lane. Your value is in reading everything and asking "what has
nobody tested?"

## 1. What you owe

1. **At least three ideas nobody else has proposed** (`s26/IDEA_<name>.md`, the contract's format:
   hypothesis; why the record does not already close it, with the sprint that would have; the
   exact falsifier; expected effect against a computed MDE; memory estimate; agent-hours). Before
   writing one, grep the existing `s26/IDEA_*.md` and the record so you do not duplicate. Ideas the
   coordinator considers unclaimed so far, for you to evaluate and either propose or reject with a
   reason: (a) the 2/60 benchmark self-copy leak bounded from a dev-set proxy without opening the
   benchmark (mandatory direction; the dev version costs +0.0004 A on 4/126 targets: 1CEK in 1A11,
   2FBU in 2LMF, 2P5H in 2P5J, 6B9K in 1U6V; construct the proxy from dev targets whose self-copy
   situation matches, state the bound and its assumptions, never read the manifest); (b) test-time
   window ensembling over BLOSUM variants or K settings averaged in coordinate space with the same
   projection (S17 says widening K hurts; you must say how this differs); (c) anything from your
   own reading. Check `s26/TOURNAMENT.md` and the ledger first: if (a) or (b) is already claimed by
   another lane, do not duplicate it.
2. **The highest-ranked orphaned survivor of the tournament.** When `s26/TOURNAMENT.md` exists, take
   the top orphaned idea (the Adversary lists them), write its PREREG, and run it under the
   governor in rank order. If the tournament is not ranked yet, write PREREGs for your own ideas and
   probe their memory (1-target probes under jobrun) so they can start the minute they are ranked.
3. `s26/agentW_FINDINGS.md` in the S12 to S25 format; ledger entries for every finding; hourly
   status lines under `## W (Wildcard)` in `s26/STATUS.md`; commits of your own files
   (`s26/w_*.py`, `s26/results/w_*.json`, your IDEA and PREREG files) with the contract's trailers.

## 2. Discipline reminders that bite in a free lane

- The phase gate: no RMSD read until `PHASE 0 SIGNED OFF` is in the ledger.
- Every experiment: PREREG first, `ST.compare` with fold-clustered CI beside iid, MDE, controls
  matched in the operator's space, zero-information control (constant alpha-helix), ORACLE
  labelling, `argmin_tied`, provenance-stamped results, replication on a second seed and fold order
  for anything positive, a power statement for anything null.
- Governor: `python s26/jobrun.py --agent W --tag CPU|ESM --name <name> --est-ram <GB> -- python
  s26/w_<x>.py`; never load `esm_cache.npz` (1.5 GB) without a probe; about 4 GB of headroom is
  shared by all lanes; read `s26/governor_state.json` before launching.
- Never open the benchmark. Never touch pinned artefacts. New code under `s26/` only.
- No stock words, no em dashes in anything the presenter or the report will read.

Finish your turn with what you proposed, what you ran, what you did not do and why, and any question
for the coordinator. Do not predict results you have not measured.
