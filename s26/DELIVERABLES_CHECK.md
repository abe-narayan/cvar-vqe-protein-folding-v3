# S26 DELIVERABLES CHECK (lane A, the Adversary) -- FINAL PASS

Walks the Part 10 table of the campaign prompt (14 rows). Final pass started 2026-09-14 04:00
on the committed tree (ledger at L133); the report row is re-checked at about 04:40 on lane E's
committed final text and its result is appended below. Each row: exists / meets the stated
standard / what is missing. Rule 13 greps (banned words; U+2014; U+2013) were run over every
presenter-facing file and the deck dump at 04:02; the results are in the second table.

| # | deliverable | exists | standard met | missing / note |
|---|---|---|---|---|
| 1 | governor + log (`s26/governor.py` 471 lines, `governor.log` 1,946 lines, `jobrun.py`, `enqueue.py`) | OK | v2 band; v2.1-v2.4 fixes each ledgered (L18b, L36, L59, L74, L83, L92, L111, L128/L129); host kill at 10:10 and the stall breaker recorded; log whitelisted and continuous | -- |
| 2 | EXAMINATION + AUDIT | OK | `EXAMINATION.md` (A-I + addendum), `EXAMINATION_AUDIT.md` (L31); one MATERIAL fixed (L32), re-checked (L34); Phase 0 signed off (L33) | -- |
| 3 | BRIEF | OK | `BRIEF.md`: one-sentence falsifiers per item, basis named | -- |
| 4 | PREREG_* (29 files) | OK | every arm that ran has a prereg written before compute (checked on every entry the Adversary examined, L45-L97, L120-L121); addenda appended, never edited | -- |
| 5 | PROPOSAL_A / B / C + B_REPLACEMENT | OK | all four FINAL; A REPLACE (L69/L70/L75), B REPLACE (L120: notes sourced, isolations persisted `s26/results/a_ladder_isolations.json`), C KEEP WITH EDITS (L120; the MATERIAL document defect fixed by L124's addendum 1; C5 final in addendum 2 per L132); scripts 192 / 214 words | B's note [5] rho -0.638 cites `p_b3.json` but lives in `agentP_FINDINGS.md` section 10 (L124 adopted the note); MINOR |
| 6 | IDEA_* (18) + TOURNAMENT | OK | ranked with kill attempts, run order, memory-deferred and orphans (L51 accepted); every ranked survivor ran or was deferred with its reason (L133 for raw) | -- |
| 7 | agent*_FINDINGS (E, I, Q, P, PH, W, A, PR) | OK | all eight on disk in the S12-S25 format; A's carries the per-slide qualifier table | -- |
| 8 | LEDGER + RETRACTIONS | OK | 137 entries, append-only, collisions suffixed; `RETRACTIONS.md` R1-R10 plus the prior-sprint disposition table (R2 the depth-3 scope; the S13 Pauli weights sourcing gap; S7-11 -0.288 replaced by L62; S10-4 re-derived; the suite's basis; S16 sharpened; routers extended; L38 scoped by L89) | -- |
| 9 | repository fixes on `s26` | OK | items 1-5 and defects 6a-6d closed (L1, L6, L7/L9, L8, L10, L15, L16, L19, L20, L21); `s24/cache_amber` tracked (126 files); AST identity checked; the opt-in tier 11/11 (L113, after L102/L104); the frozen lab rebuild reproduces 2016/2016 RMSDs (L127) | L127's three descriptive leaderboard columns written by uncommitted code is lane I's item, recorded there |
| 10 | the deck + PRESENTATION_CHANGES | OK | `vqe_research_overview.pptx` (900 KB, built 03:55 from artefacts, L61/L73/L76/L126); `PRESENTATION_CHANGES.md` (584 lines) and `pr_notes.md` map every number to its path; the dump `pr_verify_dump.txt` is rule-13 clean; the Adversary's qualifiers are on the slides (L46, L54, L55, L70, L75/L82, L122, L123) | the presenter says each qualifier with its number (table in `agentA_FINDINGS.md`) |
| 11 | TRAINABILITY_PAPER_OUTLINE | OK | on disk (184 lines), rule-13 clean; the plateau claim is placed against the depth-3 regime (line 127) and "a barren plateau from a small gradient" is in the claims-not-to-make list; "positive" reads "exact" (L122) | MINOR: put the literal "at depth 3" in the F5 caption |
| 12 | REPORT (`s26/REPORT.md`) | PARTIAL at 04:00 | 3,423 lines at commit `439ce5f8` (02:21); rule-13 clean at 04:02; lane E's final pass (Parts VII, VIII, Appendix B, Appendix C) runs 04:15 to 04:45 | re-checked at ~04:40 on the committed text; result below |
| 13 | STATUS | OK | hourly lines per lane, timestamped (230 lines) | -- |
| 14 | docs/FINDINGS.md corrections ledger | OK | "Sprint 26 additions" block at line 117 (L130, commit `5f841706`); rule-13 clean; R6-R10 added to RETRACTIONS at 04:05, so the coordinator appends matching rows (L130's own clause) | PENDING: the R6-R10 rows |

## Rule 13 pass (banned words; U+2014; U+2013), 2026-09-14 04:02

| file | banned | U+2014 | U+2013 |
|---|---|---|---|
| `s26/REPORT.md` (at `439ce5f8`) | 0 | 0 | 0 |
| `s26/PRESENTATION_CHANGES.md`, `s26/pr_notes.md`, `s26/pr_verify_dump.txt` (the deck dump), `s26/pr_verify.txt` | 0 | 0 | 0 |
| `s26/PROPOSAL_A.md`, `PROPOSAL_B.md`, `PROPOSAL_C.md`, `PROPOSAL_B_REPLACEMENT.md` | 0 | 0 | 0 |
| `s26/TRAINABILITY_PAPER_OUTLINE.md`, `TOURNAMENT.md`, `RETRACTIONS.md`, `EXAMINATION.md`, `EXAMINATION_AUDIT.md`, `BRIEF.md`, `agentA_FINDINGS.md` | 0 | 0 | 0 |
| `docs/FINDINGS.md`, the S26 block | 0 | 0 | 0 |

## Numbers that must not reach a slide or the report without their qualifier (final list)

- C34 (|z_moment| triple), C24's 0.524 and the S13 Pauli mean weights 2.236 / 3.015: off every slide (L32, L122).
- The seven-configuration suite, the random-75 null 3.425, Legacy +0.330 / AMBER +0.455: point-cloud basis on both sides (L28/L29).
- C3: "+0.0111 worse than random" is Type-M; lean on +0.0207 and +0.0385; "validity step on 124 of 126" (L46); replicated (L87).
- Steric reject: +0.228 / +0.248 and +0.167 / +0.157 are Type-M, tail-carried; direction measured at 1e3 and by S@1e4, on both bases (L54, L95).
- The 2/60 bound: 0.028 (mean-CI envelope) beside 0.151 (worst target), MINOR under every reading, A2 named (L55, L58).
- A1: "twelve of twelve" is one observation (corr 0.955); 0.23x / 0.36x MDE, 0.06 A resolution; growth at alpha = 1 is inert, not absent (L70, L75, L82).
- The tie-break floor 0.024 A applies across runs that do not share the tie-break, not to paired contrasts (L71).
- The omega floor: 0.35 A is the own-torsion upper bound, 0.083 A the tight floor (L89).
- The strain flag: the pool's own spread predicts the error (rho +0.45); `moved` is its proxy and adds nothing given it; never "the first native-free quantity above 0.4" (L121, L123).
- The toward-member 0.02 A is a partial reversal of the projection's displacement, not a lever (L95); the branch degeneracy 0.08 A is an ORACLE order statistic (L96).
- "No exponential plateau": always "at depth 3" (R2, L119).
- 3.2148 (cache) vs 3.2126 (rebuild): quote one and name it; the C2 anchor is 3.2126 (L57).
- C5: null (GLOBAL) to harmful (RIDGE +0.164, 1.79x, 5/5); the -1.87 / -3.13 A ceilings are ORACLE (L132).
- The raw rung was not run (4 of 5 folds trained, L133): the ladder is nine of ten rungs; say so.

## Report re-check (appended at about 04:40)

(pending)
