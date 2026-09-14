# S26 DELIVERABLES CHECK (lane A, the Adversary) -- LIVING FILE

Walks the Part 10 table of the campaign prompt (14 rows). Updated as artefacts land so the
final pass before close is short. Each row: exists / meets the stated standard / what is
missing. Banned-word and em-dash greps (rule 13) are run over every presenter-facing file each
pass: `grep -nE "genuinely|honestly|leverage|robust|delve|underscore"` and a count of U+2014.

Last pass: 2026-09-14 00:15 (after L97).

| # | deliverable | exists | standard met | missing / note |
|---|---|---|---|---|
| 1 | governor + log (`s26/governor.py`, `s26/governor.log`, `jobrun.py`, `enqueue.py`) | OK | OK: v2 band, whitelisted log, jobrun v2.2 (L18b, L36, L41); host-kill at 10:10 recorded (L40) | -- |
| 2 | EXAMINATION + AUDIT (`s26/EXAMINATION.md`, `s26/EXAMINATION_AUDIT.md`) | OK | OK: audit L31, one MATERIAL fixed by addenda (L32), re-checked L34; Phase 0 signed off L33 | -- |
| 3 | BRIEF (`s26/BRIEF.md`) | OK | OK: one-sentence falsifiers per item, basis named | -- |
| 4 | PREREG_* (22 files) | OK | OK for the arms that ran; survivors still without a PREREG (L51): strain_difficulty, window_ensembling, window_provenance, amber_prior_partner, product_state_optimum (unless inside PREREG_A1) | PENDING: those five before compute |
| 5 | PROPOSAL_A / B / C (+ B_REPLACEMENT) | PARTIAL | all four on disk; A's verdict REPLACE accepted (L69) and supported by the Adversary (L70); B_REPLACEMENT per L13 | C's verdict waits on the ladder eval (L62-L67 so far: noesm +0.208 Type-M, conly / wide / pca32f not measured) |
| 6 | IDEA_* + TOURNAMENT | OK | 16 IDEA files; `s26/TOURNAMENT.md` ranked (L51 accepted); banned-word clean | -- |
| 7 | agent*_FINDINGS | PARTIAL | A, E, I, P, PH, Q, W on disk in the S12-S25 format; PR lane (spawned L42) has none yet | PENDING: agentPR_FINDINGS; each lane's "what I did not do" section present (checked A, E, W) |
| 8 | LEDGER + RETRACTIONS | OK | append-only, numbered, collisions suffixed (L16b, L18b); `s26/RETRACTIONS.md` R1-R5 (R5 = L75's inert-growth correction), R2 marked FOR docs/FINDINGS.md | keep current: no S26 endpoint has contradicted a prior sprint yet (L39, L43, L44 confirm the record) |
| 9 | repository fixes on `s26` | OK | lane I's L6, L7/L9, L8, L10, L15, L16, L19, L20, L21 with AST-identity checks and governed tests (370/357/13); no production number changed | -- |
| 10 | presentation + PRESENTATION_CHANGES | MISSING | the pptx was never on this machine (L2); the PR lane builds the deck from the prompt's structure with python-pptx and writes `s26/PRESENTATION_CHANGES.md`; neither on disk at this pass | PENDING (PR lane). Every slide number must carry an artefact path; C34 and the 0.524 stay off every slide (L32); C26 now has `s26/results/a_c26_phi_mae.json`; the seven-configuration suite and +0.330/+0.455 carry the point-cloud basis (L29); the C3 "+0.0111 vs random" carries its Type-M flag (L46); the 2/60 caveat reads per L55 |
| 11 | TRAINABILITY_PAPER_OUTLINE | OK | on disk (lane Q); the "no plateau" claim must carry "at depth 3" (RETRACTIONS R2) | check the outline text for "at depth 3" at the final pass |
| 12 | REPORT (`s26/REPORT.md`) | PARTIAL | 2,239 lines, Parts I-VI, IX, Appendix A drafted (L37); banned words 0, em dashes 0 at this pass. A second report `docs/REPORT_S26.md` + `docs/REPORT_S26_SUMMARY.md` exists that no S26 lane wrote (L42) | PENDING: Parts VII, VIII, Appendix B; the docs/ pair must be reconciled or quarantined by the coordinator before close |
| 13 | STATUS (`s26/STATUS.md`) | OK | hourly lines per lane, timestamped | -- |
| 14 | docs/FINDINGS.md corrections ledger | PENDING | untouched since `ae86a124`; the coordinator makes the edit at close from RETRACTIONS.md (R2 the depth-3 scope; the C27 sourcing correction is a note, not a retraction) | PENDING at close |

## Banned-word / em-dash pass (rule 13), presenter-facing files

| file | banned words | em dashes (U+2014) | pass |
|---|---|---|---|
| `s26/REPORT.md` | 0 | 0 | 2026-09-13 19:45 |
| `s26/TOURNAMENT.md` | 0 (one rule-13 word removed, `e3ea9f33`) | 0 | 19:45 |
| `s26/RETRACTIONS.md` | 0 | 0 | 19:45 |
| `s26/EXAMINATION_AUDIT.md` | 0 | 0 | 19:45 |
| `s26/agentA_FINDINGS.md` | 0 | 0 | 19:45 |
| `s26/LEDGER.md` (A's entries) | 0 (one rule-13 word removed) | 0 | 19:45 |
| deck / PRESENTATION_CHANGES | built (L61), slide 8 rebuilt (L73), L75 wording applied (L76); grep at the final pass | -- | -- |

## Numbers that must not reach a slide or the report without their qualifier

- C34 (|z_moment| triple) and C24's 0.524: document-only, off every slide (L32).
- The seven-configuration suite 3.058..3.881, the random-75 null 3.425, Legacy +0.330 / AMBER
  +0.455: point-cloud basis, name it on both sides (L28/L29); built-chain twins in
  `results/summary/leaderboard.json :: rows[*]/mean`.
- C3: "+0.0111 worse than random" is Type-M (1.12x); lean on +0.0207 vs do-nothing and +0.0385
  vs member; "validity step on 124 of 126" (L46).
- Steric reject: "+0.228" and "+0.167" are Type-M; the harm is tail-carried (median +0.003); the
  direction is measured at 1e3 and by S@1e4 (L54).
- The 2/60 bound: "0.028" is the mean-CI envelope; the worst single target gives 0.151; MINOR
  under every reading; quote both with A2 named (L55).
- A1: "twelve of twelve arms negative" is one correlated observation (mean pairwise corr 0.955); quote it with the 0.23x / 0.36x MDE and the 0.06 A resolution (L70).
- The tie-break floor 0.024 A applies across runs that do not share the tie-break; paired within-run contrasts (C3 etc.) are not inside it (L71).
- The cis/omega floor: "0.35 A" is the own-torsion upper bound (L38); the tight floor is 0.083 A (L89); the projection's 0.166 A cost is displacement, not representation.
- The toward-member 0.02 A improvement (L87) is a measured reversal of part of the projection's displacement, not a physics result and not a proposal (L95); the branch-degeneracy 0.08 A is an ORACLE order statistic (L96).
- "No exponential plateau": always "at depth 3" (R2).
- 3.2148 (cache) vs 3.2126 (leaderboard rebuild): quote one and name it (EXAMINATION H).
