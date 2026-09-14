# agentA_FINDINGS (Sprint 26, lane A: the Adversary)

S12 to S25 format: tiers DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN, every
number with its artefact path, then "what damaged my own expectations" and "what I did not do
and why". Companion documents: `s26/EXAMINATION_AUDIT.md`, `s26/TOURNAMENT.md`,
`s26/RETRACTIONS.md`. My scripts: `s26/a_reproduce_head.py`, `s26/a_c26_phi_mae.py`,
`s26/a_dla_check.py`, `s26/a_append_checks.py`.

## DEMONSTRATED

D1. **Phase 0 audit: reproduction exact on HEAD; one MATERIAL documentary finding, fixed.**
`s26/EXAMINATION_AUDIT.md` (ledger L31). `s26/e_reproduce.py` re-run on HEAD through
`s26/a_reproduce_head.py` (job `a_reproduce_head`, exit 0; `s26/results/a_reproduce_head.json`)
equals the stored `s26/results/e_reproduce.json` at 0.0 on all four bases (3.048338 / 3.204076 /
3.214765 / 3.235460) and all 126 rows; T030 = 1S9Z rmsd_arm 0.18198112330908295. 29 claim-ledger
leaves opened; every sourced value sits at its cited key at stored precision. The one MATERIAL
finding (the test total 369/356 against the artefact's 370/357) was fixed by L32 and re-checked
in L34 (STANDS, no veto).

D2. **C26 (phi 36.1 vs 36.4 deg) is now sourced.** `s26/results/a_c26_phi_mae.json` (provenance
stamped): pooled over 1,507 residues of `s13/cache/tors_rows.npz`, phi MAE p_grid 36.133 /
n_marg 36.416 deg, psi 62.355 / 72.772; matches `s13/SPRINT13_DOSSIER.md:261-263` at one decimal.
The document-only number now has an artefact path for the deck.

D3. **L27 (the DLA) reproduces independently.** `s26/a_dla_check.py` ->
`s26/results/a_dla_check.json` (88 s), sharing no closure code with `s26/q_dla.py` and using a
different numeric rank rule, agrees with `s26/results/q_dla.json` at every fixed-ansatz cell
(n = 4..7, L = 1..4), across all four conjugation/gate conventions (8128 at n = 7 from L = 2),
and on the numeric route at n = 4, 5 (6/6). Ledger L45: STANDS.

D4. **L22, L23, L24, L38, L39, L35 recomputed from their artefacts.** Ledger L46-L48. Every
number in the PH censuses, the C3 stage-1 contrasts and the cis floor reproduces; the AMBER-cache
pool-identity assertion holds on three spot-checked targets; the C3 `rand` control matches
S16's `s16/repair.py` construction exactly.

## ORACLE DIAGNOSTIC

None of my own; my checks re-scored other lanes' ORACLE artefacts (L38, L39) without adding one.

## HYPOTHESIS

H1. The `rand` control in L39 beats AMBER, and a move TOWARD a random pool member beats
do-nothing (-0.0178, fold CI [-0.0255, -0.0124], 5/5, `s26/results/ph_c3_stage1.json`). A
zero-information move toward a member improving RMSD is a striking native-free signal; it is a
control here, not a proposal, and it is small. Noted, not pursued (it would be an endpoint arm).

## REFUTED

None of my own claims. I recorded four other lanes' or sprints' claims as retracted or scoped in
`s26/RETRACTIONS.md` (R1 the S26 test count; R2 the S25 "no plateau" scoped to depth 3; R3 lane
Q's H2b; R4 lane PH's cis-ensemble prediction).

## OPEN (my queue after this turn)

- Adversary checks of the endpoint results that landed while I ranked: L43 (the steric-reject
  point-cloud arm, +0.228 worse -- a confirmed-harmful, so leakage/control discipline to verify,
  not a gain to attack) and L44 (the 2/60 bound, MINOR; verify the triangle bounds and the
  envelope fold CI against the L30 prereg).
- The remaining endpoint results as they post (A1, A3, the C2 ladder rungs, C4, C5): each gets a
  `## L<n> -- ADVERSARY CHECK` within the hour, the Part-4.3 checklist answered item by item.
- The deliverables check (brief section 5) before the sprint closes, including the banned-word /
  em-dash grep over the report and the deck.

## What damaged my own expectations

- The identical pool index 449 on both traced targets read as a trace bug; it is a coincidence
  verified in both cache records (`sub[44]` = 449 on 1S9Z, `sub[34]` = 449 on 9KAR).
- I expected C27's S10 artefacts to be lost (the examination read them absent); they are in git
  history at `5fa05cd`, and lane W re-derived the +0.0004 dev half exactly (L44).
- I expected an independent DLA closure to disagree at n = 6 or 9 (the "divisible by 3" pattern);
  it agrees at every cell, including the intermediate 510 / 1023 at n = 6.
- I expected no clean kill in the tournament and found none: every closing experiment the record
  offers measured a different operator or space, so all 16 ideas survive as scored hypotheses.

## What I did not do and why

- Did not open `s10/idaudit_*.json` from git history: its price stage covers the sealed benchmark
  (Rule 1, rule 14).
- Did not re-hash `results/benchmark_manifest.json` (rule 14); the Examiner's byte hash equals
  the S20 record and `s26/e_hashes.py` never parses it.
- Did not run any endpoint experiment of my own; every check re-read an artefact or re-derived a
  property. The only governed job I launched is `a_reproduce_head` (5 s); `a_dla_check.py` and
  `a_c26_phi_mae.py` are tiny and ran directly / through jobrun.
- Did not compute a slope CI for L35's grown-vs-fixed comparison myself; I flagged that the
  artefact carries none and that the "within error" claim is therefore asserted, not computed.


## Claims that reached a slide, and the qualifier the presenter says with each (as of 2026-09-13 22:20; slides per `s26/PRESENTATION_CHANGES.md`, L61 / L73 / L76)

| slide | claim on the slide | say with it | source of the qualifier |
|---|---|---|---|
| 2 | phi MAE 36.1 (full sequence context) vs 36.4 deg (sequence-blind) | pooled over 1,507 residues of 126 targets; the 0.3 deg gap is the whole measurable sequence-to-phi channel at this length; psi 62.4 vs 72.8 | `s26/results/a_c26_phi_mae.json`, L31 |
| 2 (notes) | the dev identity-leak price +0.0004 A | on the lam = 0 chain (`fit`), 0.33x its MDE; +0.0018 on the built chain with a CI spanning zero; the +0.0030 benchmark half is S10-4's historical figure and is not on the slide | L44, L50 |
| 3, 4 | 3.2148 A built chain (production) | the cache's mean; the results-lab rebuild reads 3.2126 (PDB round trip); quote one and name it | EXAMINATION H, L57 |
| 3, 4 | 3.0483 A raw average | a point cloud, 22.3% contracted, not a structure; never compared with a built chain | E section H, L28 |
| 4 | random-75 3.425, constant helix 4.065 | point-cloud and ladder bases; the built-chain twins of the physics rows are 3.2187 to 4.1015 | L28 item 3, L29 |
| 4 | benchmark +0.0103, CI [-0.160, +0.180] | the artefact is `s9/final_report.json`, never opened in S26; the two means beside it are pinned by passing tests; no validated improvement | C04-C06, L32 item 2 |
| 4 | the 2/60 self-copy bound 0.028 A | the mean-CI envelope under assumption A2; the same envelope's worst single target gives 0.151 A; MINOR under every reading; dev-proxy price 0.002 A; cannot move the benchmark verdict either way | L55, L58 |
| 5 | 0.902 nats from the Gibbs optimum; 5.6e-17; cos 1.000000000 | on the 78 alpha = 1 targets the optimum is a product state (KL to the product of marginals at most 7.9e-4, per-target rows); the readout is insensitive at 0.24x MDE | L68, L70 |
| 6, 9 | the S25 slopes (-0.649 to -0.047 log2 Var per qubit) and "no exponential plateau" | always "at depth 3"; the DLA is the full so(128) from depth 2, so the slopes are a shallowness statement, not a favourable-algebra statement; never "barren plateau" either way | R2, L45 |
| 6 | grown circuits give no width-scaling argument (A4) | at alpha = 1 the grown circuits are product circuits; the one non-trivial family (alpha = 0.25, L2) decays at -0.302 against the fixed -0.243 with no persisted slope CI, so "within error" is asserted | L47 |
| 7 | Legacy +0.330 / AMBER +0.455 worse than a random subset | point-cloud basis, both sides; 5/5 folds | L28, L29 |
| 7, 10 | "refine with physics is a validity step, not an accuracy step"; +0.0207 worse than doing nothing | +0.0207 at 2.16x MDE is clean; "+0.0111 worse than a random move" is Type-M (1.12x): say the sign, not the size; "validity step on 124 of 126" (2BP4 and 9KAR break a virtual bond; 9KAR does not converge) | L46 |
| 7 | AMBER as a steric reject filter is harmful (+0.228 at 1e4) | Type-M magnitude (1.17x); the direction is measured at 1e3 and by the shrink arm; the harm is tail-carried (median +0.003, worst +4.15); point cloud so far, the built chain pending | L54 |
| 7 (notes) | strain predicts error, Spearman +0.433 | a calibration flag, never a gain; the signal is the pool's own disagreement (the top-75 spread, native-free, rho +0.45 before AMBER runs); `moved` tracks it at rho 0.76 and adds nothing given it (partial +0.08, CI spans zero); "CIs replicated", not "result replicated"; never "the first native-free quantity above 0.4" (vetoed as worded) | L81, `s26/results/a_strain_vs_spread.json`, follow-up entry |
| 10 (notes) | C5 "not run to completion" | C5 was RUNNING at 02:06 with PROPOSAL_C stamped FINAL 03:30; the outcome is whatever the close finds, reported in the present tense; not a result either way | L120 |
| 9 | the S13 Pauli mean weights 2.236 / 3.015 | typed from the claim ledger by PR (L61), not read from an artefact; re-read from `s13/results/geo_pauli.json` before the slide says them, or say "every S26 figure" | L120 |
| 8 | A1: ADAPT is null, -0.014 / -0.022 A at 0.23x / 0.36x MDE | the comparison resolves 0.06 A; "twelve of twelve arms negative" is one observation (pairwise correlation 0.955); the 21-parameter match is a budget with realised counts 7 to 21; growth at alpha = 1 is inert, not absent (L75) | L70, L75, L82 |
| 8 | the fixed circuit stops 0.90 nats short; ADAPT reaches the Gibbs state | on the 78 alpha = 1 targets; on the 48 alpha = 0.25 targets the Gibbs state is not the CVaR optimum (gibbs_T is +0.067 there) | L70 caveat 3 |
| 10 | the C2 anchor 3.2126 and the rung contrasts | rebuild basis (L57); noesm +0.208 is Type-M with the sign at 5/5; every other rung so far is under its MDE ("not measured", MDE beside it); any best rung at the close is an order statistic (`best_of_k_within`) | L79 |
| 7 (notes) | the constant-omega representation floor 0.35 A (L38) | the own-torsion UPPER bound; the tight floor (native projected through the production projection) is 0.083 A, 0.043 with the prior off; the projection's 0.166 A cost is its displacement, not representation | L89, L97 |
| 7, 10 (notes) | a zero-information move toward a pool member improves the built chain by 0.02 A (L87) | measured and replicated (1.77x MDE, 5/5), but a partial reversal of the projection's own 0.166 A displacement toward the cloud the members surround; not physics, not a proposal; the same step toward the cloud itself was not run | L95 |
| 7 (notes) | the projection-branch degeneracy is worth 0.08 A to a perfect chooser (L88) | an ORACLE order statistic; the 56% split-half transfer is against the random pick and lands on the production choice, so no fixed-start rule beats the objective | L96 |
| 11 | -- | no number of mine; the tie-break floor 0.024 A, if quoted, applies across runs that do not share the tie-break, not to paired within-run contrasts | L71 |

## OPEN (queue as of 2026-09-13 22:20)

- The pool-spread control for L53 (`s26/a_strain_vs_spread.py`, job queued behind the four-job cap); L81 is provisional until it lands.
- Every new entry within the hour: the reject chain, A3, the remaining rungs (esm8m after L78's kill, mix, pairnet, raw), branch_select, rotamer_relief B, the ensembling, W's survivors, lane I's slow test tier and the results-lab rebuild.
- RETRACTIONS.md (R1-R5) and DELIVERABLES_CHECK.md kept current; the final Part 10 pass at about 04:00 with the rule-13 grep over the report, the deck dump and every presenter-facing file.

## Final pass (2026-09-14 04:00 to 04:12; ledger entry "DELIVERABLES CHECK (FINAL PASS)")

`s26/DELIVERABLES_CHECK.md` final: 13 of 14 rows OK; the report row re-checked at about 04:40
on lane E's committed text (result in the ledger). Rule 13 clean on every presenter-facing file
and the deck dump. `s26/RETRACTIONS.md` R1-R10 plus the prior-sprint disposition table for the
`docs/FINDINGS.md` corrections block (L130). L131 and L132 checked: both STAND.

## What I did not do (final)

- L115 (C4, eighteen routers) and L112 / L116 (the delivery file) were not given their own
  check entries; both are nulls or reproductions and nothing positive stands unchecked.
- The S13 locality theorem and the S13 Pauli-spectrum prediction were cited, not re-derived,
  this sprint (L120, L122); the S13 Pauli mean weights are off the slides for that reason.
- The `attn` rung (3 to 3.5 GB) and the raw rung's evaluation (L133) did not run on this box;
  they are the two prior inputs the ladder could not price.
