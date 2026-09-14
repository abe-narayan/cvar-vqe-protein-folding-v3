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
