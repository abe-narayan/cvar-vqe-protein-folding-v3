# LANE PR (PRESENTATION), SPRINT 26: FINDINGS

Deliverable: `vqe_research_overview.pptx` at the repository root, built by `s26/pr_build_deck.py`
from artefacts; the record `s26/PRESENTATION_CHANGES.md` (written by `s26/pr_changes.py`); the
spoken text `s26/pr_notes.md`; the number registry `s26/pr_values.py` -> `s26/pr_values.json`;
the figures `s26/pr_figures.py` -> `s26/figures/pr_*.png`; the verification records
`s26/pr_verify.txt` and `s26/pr_verify_dump.txt`. Tiers as in S12 to S25. Every number carries
its artefact path; the full table is `s26/PRESENTATION_CHANGES.md`.

## 0. The case I was in

The deck named in the campaign prompt does not exist on this machine or in git history (ledger
L2), so I BUILT it from the eleven-slide structure the prompt describes rather than editing a
file. `s26/PRESENTATION_CHANGES.md` maps every slide's content and every number's artefact by
title, so the edits can be carried onto the real file if the presenter supplies one.

## 1. What is built. DEMONSTRATED (a file, reopened and checked).

- Eleven slides, dark theme (RGB 18, 18, 24; off-white text; one amber accent; figures on white
  plates at 190 dpi; Calibri; title + body + notes on every slide). Reopened with python-pptx:
  11 slides, 0 U+2014, 0 U+2013, 0 occurrences of any of the six banned words, the only
  non-ASCII character U+2022 (`s26/pr_verify.txt`).
- Slides 1 to 7 and 11 (step a): every number read from an artefact at build time by
  `s26/pr_values.py` (237 registered tokens, `s26/pr_values.json`), 30 of the 35 claims of
  `s26/EXAMINATION.md` section C re-read at their leaves. Slide 4 carries the two
  ORACLE-superposed CA overlays: 1S9Z (T030, n = 16, built chain 0.181981 A) and 9KAR (n = 15,
  7.437696 A), the built chain `ca` from `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`
  against `nat_ca` from `s8/generate_univ/<pdb>.npz`, superposed with
  `s12.instrument.superpose_batch`; the RMSD recomputed by `s26/pr_figures.py` equals the
  record to 1e-9 on both (`s26/pr_figures.py`, assertion). Slide 6 carries my rebuild of the S25
  width sweep from `s25/results/q_plateau.json` (five cells, fitted slopes, the S13 depth-8
  reference base 0.5039 from `s13/results/geo_kernel.json`) beside lane Q's A2 figure; slide 7
  the accuracy ladder with the basis on every bar.
- Slide 10 (step b) from `s26/PROPOSAL_C.md` and `s26/C3_RESULT.md` with the Adversary's L46
  wording; slide 9 (step c) from `s26/PROPOSAL_B.md` and `s26/PROPOSAL_B_REPLACEMENT.md`;
  slide 8 (step d, the "until then" state) a PENDING placeholder with the measured A2 / A4
  facts and lane Q's A4 figure; the builder swaps it for the proposal slide when
  `s26/PROPOSAL_A.md` exists. Spoken word counts: 234 / 248 / 244 on slides 8 / 9 / 10, all
  under 250; the ALSO and SOURCES blocks in the notes are not spoken.
- Slide 11's notes carry the direction line, marked DRAFT (step e): the evidence favours
  direction C in its learn-the-prior form as the accuracy study (the prior is the only lever
  measured as steep, `s24/results/priorladder.json`, -2.15 A per unit), paired with direction A
  as a scientific study rather than an accuracy study (the product-circuit result,
  `s26/results/q_dla.json` adapt sets and `s26/results/q_var.json`); direction B as stated
  cannot be measured on this machine (`s26/results/b1_feasibility.json`).
- Basis rule (L28 item 3, L29 item 2) applied: the seven-configuration suite, the random-75
  null 3.4251 and the Legacy +0.330 / AMBER +0.455 verdicts are named point-cloud everywhere
  they appear; the production 3.2148 is named built chain; the C2 anchor is named the
  leaderboard-rebuild basis 3.2126 (L57); the S25 quantum contrasts are named the selection
  basis; every RMSD in the notes' SOURCES list carries a `[basis]` tag.

## 2. Numbers I could not source to a results artefact, and what I did with them

| number | status | what the deck does |
|---|---|---|
| +0.0103 [-0.1596, +0.1803], 31W/29L (benchmark, claim C06) | the artefact `s9/final_report.json` is not opened by any S26 lane (Rule 1) | on slide 4, typed from the claim ledger, labelled NAMED_NOT_OPENED in the registry; the two means 2.9507 / 2.9610 are SOURCED_BY_TEST (`tests/test_pipeline.py:338-347`) |
| 4.4 GB of campaign headroom at the B1 measurement | a governor reading recorded in L13, no JSON leaf | on slide 9 with the ledger entry as its source |
| the +0.0030 A benchmark leak price (C27) | artefact absent from the tree (in git history at `5fa05cd`, not opened) | NOT on any slide; the benchmark caveat uses lane W's bound (`s26/results/w_selfcopy_bound.json`, L58) instead |
| the abs z_moment triple 0.7529 / 0.8013 / 0.1127 (C34) and the 0.524 cosine (C24) | document-only (L28, L31, L32) | NOT on any slide or note |
| ESM -0.288 A [-0.484, -0.092] on selection (S7-11) | per-target artefact lost (L11) | NOT on slide 9; named as document-only in its notes |
| 355 / 13 test count (C35) | superseded | replaced by 370 / 357 / 13 (`s26/results/test_run.json`), registered and not on a slide |
| the S13 "2.236 / 3.015" Pauli mean weights and the "78 cells" of the outline | I could not reproduce them from `s13/results/geo_pauli.json` in the time I had (the cell means there are 2.93 / 5.14 under a different conditioning) | NOT on slide 9; the ratio I did reproduce (median measured/predicted 0.9969 over 104 cells, `cells[*]/legacy/ratio_meas_over_pred_exact`) is what the slide quotes |
| the operator law d_out = 1.16 mean + 0.04 best (S12) | a memory / dossier number, no leaf named in the claim ledger | slide 3 states the mechanism ("the output tracks the set mean") without the coefficients |

## 3. What I changed against the record's older wording, and why

Listed in `s26/PRESENTATION_CHANGES.md` ("Stale or unqualified numbers"): the point-cloud basis
named on the physics suite; "no barren plateau" replaced by "gradient variance at depth 3" with
the five slopes, the DLA and the scope; "AMBER refines" replaced by the L39 / L46 sentence with
the Type-M flag on +0.0111 and "a validity step on 124 of 126"; the C2 basis per L57; the
product-circuit reading grounded in `q_dla.json` per L47 and no CI claimed on the -0.302 vs
-0.243 difference; the steric reject quoted with the Type-M flag and "tail-carried" per L54; the
strain signal quoted as a calibration flag, never a gain, per L53; the identity-leak price now
carrying lane W's artefact per L44 and L50. Nothing softened: every proposal slide carries the
verdict form and its PENDING state.

## 4. Where the deck's numbers agree with the ledger by recomputation (a check on both)

- L14 (tors minus arm +0.5557, median +0.2314, SE 0.1318, MDE 0.3694, 37W/89L) recomputed from
  `s13/cache/tors_rows.npz['a_pepPos']` against the production cache: identical to four decimals.
- L35's matched-P slopes (-0.079 / +0.006 / +0.035 / -0.008 / -0.246 / -0.302) recomputed from
  the `P_matched` rows of `s26/results/q_var.json`: identical to three decimals; the grown/fixed
  ratio range 1.63 to 370 reproduced over 39 matched rows.
- L26's C1 closures (3.0433 -> 3.0258, 2.5342 -> 2.1460, 2.6087 / 3.4676 / 3.4540) read from
  `s12/results/agg_dec_v2*.json` and `s17/results/inband.json`: identical.
- L53's quartile means 2.286 / 2.936 / 3.758 / 3.923 recomputed from `s26/results/ph_strain.json`.
- L46's three C3 contrasts and the cosine -0.0491 read from `s26/results/ph_c3_stage1.json`.
- The overlays' RMSDs (0.181981, 7.437696) recomputed through `s12.instrument.ca_rmsd`.

## 5. What damaged my own expectations

1. I expected the ladder's C2 built-chain anchor to be the production 3.2148. It is the rebuild
   3.2126 (L57): the same cloud, fed to the projection at 1e-14, lands in a different local
   optimum on 120 of 126 targets. The deck names the basis; it would have been wrong by 0.002 A
   and, worse, by a basis, without L57.
2. I expected the +0.0111 (AMBER vs random) to be the headline of the C3 sentence. The Adversary
   (L46) is right that it is a Type-M number; the two clean contrasts (+0.0207, +0.0385) carry the
   slide and the +0.0111 is quoted with its flag.
3. I expected "product circuit" to be a q_var.json fact. It is a q_dla.json fact (the adapt
   sets' abelian closure and `n_distinct_ops`); L47 caught the attribution before I wrote it.
4. A verification pass that prints the banned words by name puts them into every file that
   pastes it. The checks now print them masked.

## 6. What I did not do and why

- I did not open `s9/final_report.json`, `results/benchmark_manifest.json` or any benchmark
  record (Rule 1); the four benchmark numbers on slide 4 are typed from claim C06 and labelled.
- I did not put the four document-only numbers (C24's 0.524, C26 before L32, C27's +0.0030,
  C34) or the lost ESM -0.288 on a slide.
- I did not regenerate lane Q's A2 and A4 figures; they are already at 190 dpi on white from the
  results JSON (`s26/figures/a2_dla_dimension.png`, `a4_variance_slopes.png`, 1368 x 874 and
  2040 x 843 at 189.99 dpi) and the task names them as such; my own figures (overlays, width
  sweep, ladder) are regenerated on every build.
- I did not build slide 8's proposal form: `s26/PROPOSAL_A.md` does not exist (checked 19:22 and
  at every build). The placeholder carries only measured facts and is marked PENDING.
- I did not read `docs/REPORT_S26.md` beyond its Part X presentation guide and Appendix D, and
  copied nothing from it: no external citation, no personal detail; every number on a slide is
  traced by me to an artefact or a ledger entry.
- I did not run anything through jobrun: the build peaks at 0.12 GB (`psutil` sampling of the
  child process at 50 ms; figure step alone 0.10 GB), under the 200 MB line.

## 7. Artefacts and memory

    vqe_research_overview.pptx            the deck (11 slides)
    s26/pr_values.py -> s26/pr_values.json the registry (237 tokens: value, path, basis, status, note)
    s26/pr_figures.py -> s26/figures/pr_overlay_1S9Z.png, pr_overlay_9KAR.png, pr_width_sweep.png, pr_ladder.png
    s26/pr_notes.md                        spoken text with {TOKEN} placeholders
    s26/pr_build_deck.py                   the build and the verification (pr_verify.txt, pr_verify_dump.txt, pr_manifest.json)
    s26/pr_changes.py -> s26/PRESENTATION_CHANGES.md
    peak RSS: build 0.12 GB, figures 0.10 GB (measured; both under the jobrun line)
