"""PR lane, Sprint 26: write `s26/PRESENTATION_CHANGES.md` from the build's own records.

Reads `s26/pr_manifest.json` (slide -> title, tokens, spoken words), `s26/pr_values.json`
(token -> value, path, basis, status) and `s26/pr_verify.txt` (the verification checks), and
combines them with the slide descriptions below into the record the brief asks for: every
slide's content, the artefact behind every number, the stale numbers that were fixed, and the
verification dump. Run after `s26/pr_build_deck.py`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s26", "PRESENTATION_CHANGES.md")

SLIDES = {
    1: ("Title", "Title, subtitle, four bullets (targets and folds; the built-chain result with the point cloud named as an "
        "intermediate; the quantum component's size and the null contribution; what survives), audience line.", []),
    2: ("Why 9 to 16 residue peptides are hard", "Five labelled points: the phi sequence signal; the common-mode error; the prior "
        "as the only steep lever (ORACLE, point-cloud basis); scale against the constant helix; the instrument.", []),
    3: ("The pipeline, stage by stage", "Eight boxes in two rows (sequence/ESM-2; distogram; retrieval; score and cut; CVaR-VQE "
        "selector marked OFF in production; coordinate average marked not a structure; ideal-geometry projection marked THE RESULT; "
        "optional AMBER relaxation marked a validity step) plus three notes (set mean; the rank-ladder Hamiltonian; reproduction).", []),
    4: ("Results on 126 held-out targets, with two overlays", "Two ORACLE-superposed CA overlays (1S9Z = T030, 9KAR) and the summary "
        "block: mean / median / best / worst / fractions under 2 and 3 A on the built chain; point cloud and relaxed means; the random-75 "
        "and helix controls with their bases; ORACLE pool best and top-75 best; the 9KAR mechanism; the sealed-benchmark delta.",
        ["s26/figures/pr_overlay_1S9Z.png", "s26/figures/pr_overlay_9KAR.png"]),
    5: ("The quantum component: verified, then measured", "Two columns. Verified: register, exactness vs an independent simulator, "
        "parameter-shift vs finite differences, no sampling, the rank-ladder Hamiltonian, the set-equality theorem. Measured: the "
        "optimiser trains; distance from the Gibbs optimum and the endpoint's insensitivity; selector vs argmin (UNDERPOWERED); vs "
        "Boltzmann; alpha = 1 folds; off in production; no advantage claimed.", []),
    6: ("Gradient variance at depth 3, and the algebra (not a barren-plateau claim)", "Left: the S25 width sweep (five cells, fitted "
        "slopes, the depth-8 reference) and the depth sweep at n = 7, rebuilt from q_plateau.json. Right: lane Q's A2 figure (dim(DLA) "
        "vs n). Text: the five slopes, what they say, what they do not say, the draw counts.",
        ["s26/figures/pr_width_sweep.png", "s26/figures/a2_dla_dimension.png"]),
    7: ("Where the remaining accuracy lives", "Left: the accuracy ladder with the basis on every bar (controls, production, ORACLE). "
        "Right: closed by measurement (selection, physics as selector, physics as relaxer with the Type-M flag, physics as steric "
        "filter, search) and open (the prior; the common mode). A basis line under the figure.", ["s26/figures/pr_ladder.png"]),
    8: ("Direction A: ADAPT-VQE on this Hamiltonian (PENDING)", "PLACEHOLDER until s26/PROPOSAL_A.md exists. Left: lane Q's A4 figure "
        "with the matched-P slopes and the grown/fixed ratio. Right: the A2 algebra, the product-state target, A1 status, verdict PENDING.",
        ["s26/figures/a4_variance_slopes.png"]),
    9: ("Direction B: a learned folding model as the prior, and its replacement", "Left: what the proposal says; B1 final (three "
        "independent blockers); B2 pending; B3 pending with the persisted arms; verdict PENDING B2/B3, REPLACE if null. Right: the "
        "replacement paper's claim chain, the two claims it will not make, the venue statement.", []),
    10: ("Direction C: learn a better distance prior; physics for validity, not accuracy", "Left: what the proposal says; C1 final "
         "(closures reproduced); C2 running (rungs, falsifier, anchor and its basis). Right: C3 final (the seven-contrast result with "
         "the Type-M flag and the 124-of-126 validity wording); C4 and C5 pending; verdict PENDING.", []),
    11: ("Goal and ask", "Goal (two threads), ask (the study and one question), what will not be claimed. The direction line is in "
         "the notes, marked DRAFT until the coordinator's verdict entry.", []),
}

FIXES = """
## Stale or unqualified numbers, and what the deck does instead

| item | where it stood | what the deck does | artefact |
|---|---|---|---|
| the seven-configuration suite, the random-75 null 3.4251 and the Legacy +0.330 / AMBER +0.455 verdicts quoted beside the built-chain 3.2148 without a basis | `docs/STATE_BRIEF_2026-09-12.md:57, 190-191`, `results/summary/professor_brief.md:115-125` (L28 item 3, L29 item 2) | every appearance names the point-cloud basis; the built-chain means of the same rows are in the slide 7 notes | `s25/results/phys_suite.json :: basis` = point_cloud; `results/summary/leaderboard.json :: rows[*]/mean` |
| phi MAE 36.1 / 36.4 deg (document-only, C26) | `s13/SPRINT13_DOSSIER.md`, the state brief | on slide 2 with the Adversary's derivation as its artefact (L32) | `s26/results/a_c26_phi_mae.json :: summary/arms/{p_grid,n_marg}/pooled_mae_phi_deg` |
| the identity-leak prices +0.0004 / +0.0030 (C27, artefact absent) | `docs/FINDINGS.md` | the +0.0004 dev half appears in the slide 2 notes with lane W's re-derivation as its artefact; the +0.0030 benchmark half is NOT quoted; the benchmark caveat uses lane W's bound instead (L55, L58) | `s26/results/w_selfcopy_endpoint.json :: A/ge06_mean_delta_all126`; `s26/results/w_selfcopy_bound.json :: verdict/arm` |
| the mean abs z_moment triple 0.7529 / 0.8013 / 0.1127 (C34, document-only) | `s25/LEDGER.md` L16 | not on any slide or note | none |
| the 0.524 sampled tail-only cosine (inside C24, document-only) | `docs/FINDINGS.md:3161` | not on any slide; the sourced 0.566586 is registered (Q_TAIL_COS) and not spoken | `s25/results/q_verify.json :: results/'cos(exact-expectation, TAIL baseline)'` |
| 355 passed / 13 skipped (C35, superseded) | the state brief's earlier text | replaced by the governed run 370 / 357 / 13 (registered, not on a slide) | `s26/results/test_run.json :: combined` |
| "+0.113 A CVaR contribution" | `core/pipeline.py` docstring; older documents | withdrawn (S25 L5); named as withdrawn in the slide 5 notes and never quoted as a gain | `s25/results/q_alpha.json` (the paired contrast is 0.51x MDE) |
| "no barren plateau" | earlier wording | slide 6 is titled "Gradient variance at depth 3, and the algebra (not a barren-plateau claim)"; the five slopes and the DLA are given with scope | `s25/results/q_plateau.json`, `s26/results/q_dla.json` |
| "AMBER refines the structure" | earlier wording | "refine with physics is a validity step, not an accuracy step" with the L46 qualifiers (Type-M flag on +0.0111; 124 of 126) | `s26/results/ph_c3_stage1.json`, `ph_c3_nativefree.json` |
| ESM -0.288 A on selection (S7-11, per-target artefact lost, L11) | `s26/PROPOSAL_B.md` | not on slide 9; named as document-only in its notes | none on disk (`s7/repr_tune.json` lost) |
| the C2 ladder's built-chain basis | L56's title | slide 10 names the rebuild basis 3.2126 for the anchor and every C2 contrast (L57) | `s26/results/p_ladder_report_shipped_s0.json :: arm/mean` |
| "the same within error" for the grown -0.302 vs fixed -0.243 slopes | L35 | slide 8 says both are far from a 2-design's -1 and that no CI on the difference is on disk (L47) | `s26/results/q_var.json` |
| the product-circuit reading | L35 attributed to q_var.json | slide 8 grounds it in the adapt sets of `s26/results/q_dla.json`, as L47 requires | `s26/results/q_dla.json :: results/adapt_sets` |
| the steric reject +0.228 / +0.167 | L43 | quoted with the Type-M flag and "tail-carried, median +0.003" (L54) | `s26/results/ph_reject_report.json :: report/point_cloud/1e4` |
| the strain signal (Spearman +0.433) | L53 | quoted in the notes as a calibration flag, never as a gain, as L53 asks | `s26/results/ph_strain.json` |

## Numbers on a slide or in a note that are not read from a results artefact

| token | value | why | where it is traced |
|---|---|---|---|
| BENCH_DELTA, BENCH_CI_LO, BENCH_CI_HI, BENCH_WL | +0.0103, -0.1596, +0.1803, 31W/29L | the artefact is `s9/final_report.json`, which no S26 lane opens (Rule 1) | claim C06 of `s26/EXAMINATION.md`; `docs/FINDINGS.md:4505`; `README.md:14` |
| BENCH_FULL, BENCH_SHIPPED | 2.9610, 2.9507 | the same file; the two means are asserted within 5e-4 by a passing test | `tests/test_pipeline.py:338-347` (passed, `s26/TEST_RUN.md`); claims C04, C05 |
| B1_HEADROOM_GB | 4.4 | a governor reading at the time of the B1 measurement, not a JSON leaf | `s26/LEDGER.md` L13 (64.8% of 16.75 GB used, ceiling 93%) |
| C2_N_RUNGS | 11 | the ladder's design, not a measurement | `s26/PREREG_C2.md`, `s26/PROPOSAL_C.md` |
| MDE_FACTOR | 2.8016 | a rule of the instrument | `s26/LANE_CONTRACT.md` section 2; `s25/results/q_verify.json :: 'MDE == 2.8016*SE exactly'` |
| SELFCOPY_DEV, SELFCOPY_BENCH | 4/126, 2/60 | counts | `s26/results/i_identity_audit.json` (L15, L18); `tests/test_data.py:72-74` |
| the "124 of 126" validity wording, "seven rungs trained", "33 of 126 A1 targets" | text | ledger and status facts | L46; L41 and `s26/models/p_ladder/`; `s26/STATUS.md` (Q, 19:20) |

Everything else in `s26/pr_values.json` is read from the artefact named in its `path` at build
time (status SOURCED, SOURCED_BY_TEST or DERIVED with the arithmetic stated).

## DRAFT and PENDING items

- Slide 8 is a PLACEHOLDER marked PENDING: `s26/PROPOSAL_A.md` did not exist at build time; the builder replaces the
  placeholder with the proposal slide when the file appears (`s26/pr_build_deck.py`, `builders[8]`).
- Slide 9's verdict is PENDING B2/B3 (lane P); the replacement form is shown as lane Q wrote it.
- Slide 10's verdict is PENDING C2 to C5 (lane P); C1 and C3 are final.
- Slide 11's note naming the direction the evidence favours is DRAFT until the coordinator's verdict entry in
  `s26/LEDGER.md`; it was drafted from `s26/PROPOSAL_B.md`, `PROPOSAL_B_REPLACEMENT.md`, `PROPOSAL_C.md` and the ledger.
"""


def main():
    man = json.load(open(os.path.join(ROOT, "s26", "pr_manifest.json"), encoding="utf-8"))
    vals = json.load(open(os.path.join(ROOT, "s26", "pr_values.json"), encoding="utf-8"))
    verify = open(os.path.join(ROOT, "s26", "pr_verify.txt"), encoding="utf-8").read()
    L = []
    L.append("# PRESENTATION CHANGES (PR lane, Sprint 26)\n")
    L.append(f"Generated by `s26/pr_changes.py` at {datetime.now().strftime('%Y-%m-%d %H:%M')} from `s26/pr_manifest.json`, "
             "`s26/pr_values.json` and `s26/pr_verify.txt`, which `s26/pr_build_deck.py` writes.\n")
    L.append("## The case\n")
    L.append("`vqe_research_overview.pptx` did not exist on this machine or in git history (ledger L2), so the deck was BUILT from "
             "the eleven-slide structure the campaign prompt describes (title; six 'what I built' slides; three 'where I want to go' "
             "slides; goal and ask) with python-pptx 1.0.2 by `s26/pr_build_deck.py`. If the presenter supplies the real file later, "
             "this record maps every slide's content and every number's artefact onto it by title. Dark theme: background RGB(18, 18, 24), "
             "off-white text, one amber accent, figures on white plates at 190 dpi, Calibri throughout, title + body + notes on every slide.\n")
    L.append("Build chain: `s26/pr_values.py` reads every number from its artefact (the registry, written to `s26/pr_values.json`); "
             "`s26/pr_figures.py` regenerates the overlays, the width sweep and the ladder from artefacts; `s26/pr_notes.md` holds the "
             "spoken text with `{TOKEN}` placeholders; the build fills them and appends a SOURCES list to every notes frame; then it "
             "reopens the deck and writes `s26/pr_verify_dump.txt` (every text frame and notes frame) and `s26/pr_verify.txt`.\n")
    L.append("## Verification (pasted from `s26/pr_verify.txt`)\n")
    L.append("```\n" + verify.strip() + "\n```\n")
    L.append("The only non-ASCII character in the deck is U+2022 (the bullet). Word counts are of the SPOKEN block only; the "
             "ALSO and SOURCES blocks in the notes are not spoken.\n")
    L.append("## Slide by slide\n")
    for no in range(1, 12):
        title, desc, figs = SLIDES[no]
        m = man[str(no)]
        L.append(f"### Slide {no}: {title}\n")
        L.append(desc + "\n")
        if figs:
            L.append("Figures: " + ", ".join(f"`{f}`" for f in figs) + "\n")
        L.append(f"Spoken words in the notes: {m['spoken_words']}. Tokens used on the slide's notes ({len(m['tokens'])}):\n")
        L.append("| token | value | artefact (file :: key) | basis | status |")
        L.append("|---|---|---|---|---|")
        for tok in m["tokens"]:
            d = vals[tok]
            v = d["value"]
            if isinstance(v, float):
                sv = format(v, ".6g")
            elif isinstance(v, list):
                sv = ", ".join(format(x, ".4g") if isinstance(x, float) else str(x) for x in v)
            else:
                sv = str(v)
            L.append(f"| {tok} | {sv} | `{d['path']}` | {d['basis'] or ''} | {d['status']} |")
        L.append("")
    L.append(FIXES.strip() + "\n")
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    text = open(OUT, encoding="utf-8").read()
    from s26.pr_build_deck import BANNED
    bad = [w[0] + "*" + w[2:] for w in BANNED if w in text.lower()]
    print(f"wrote {os.path.relpath(OUT, ROOT)}: {len(text.splitlines())} lines; em dashes {text.count(chr(0x2014))}, "
          f"en dashes {text.count(chr(0x2013))}, banned words present: {bad}")


if __name__ == "__main__":
    main()
