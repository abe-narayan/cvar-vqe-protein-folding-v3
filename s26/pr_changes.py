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
    8: ("Direction A: let the circuit grow (qubit-ADAPT-VQE). Verdict: REPLACE", "What we learned by letting the circuit grow (the "
        "diagnosis), from s26/PROPOSAL_A.md (L68, accepted L69). Left: lane Q's A4 figure and the A1 endpoint block (the two primaries "
        "with SE, MDE, x MDE, fold CI, W/L; NOT MEASURED at the registered threshold; the 0.06 A resolution; all 12 arms; the exact "
        "Gibbs state as a control). Right: what the proposal says; the product-state diagnosis (KL to the product of marginals, the "
        "7-rotation state, the fixed circuit's 0.90 nats, the inert L-BFGS growth on the 78 alpha = 1 targets, the full DLA, the "
        "product circuits of A4, the per-step algebra of L138); the A1 block quotes both seeds (L139) as NOT MEASURED on either, with the "
        "L140 qualifier (the deployed circuit moves 0.033 A between seeds, the grown circuits are seed-stable); verdict REPLACE and "
        "what replaces it (slide 9).", ["s26/figures/a4_variance_slopes.png"]),
    9: ("Direction B (verdict REPLACE), and what we publish: the trainability paper", "Final, from s26/PROPOSAL_B.md (a88ea259) and "
        "PROPOSAL_B_REPLACEMENT.md; verdict REPLACE (L117, Adversary L120). Left: what the proposal says; B1 final (three blockers); B2 final "
        "(noesm, conly, esm8m against the shipped prior with the isolations); B3 final (the sign classifier at its permutation null, the size "
        "of the gain carried by the pool's strand content). Right: the paper's claim chain with its Sprint 26 additions, the two claims it "
        "will not make, the venue statement; the notes scope the Adversary's checking to the S26 additions (L122).", []),
    10: ("Direction C (verdict KEEP WITH EDITS): learn a better prior; physics for validity only", "Final, from s26/PROPOSAL_C.md "
         "(0a85323f with addendum 1) and C3_RESULT.md addenda 1 to 4; verdict KEEP WITH EDITS (L117, L120). Left: the four edits; C1; C2 "
         "final on nine rungs (none beats the shipped prior; the ESM-2 650M channel; the five null rungs against their MDEs; raw 3 of 5 "
         "folds, not evaluated); C4 final (12 m* and 6 s* routers). Right: C3 final with the L46 flags, the replication (L87), the validity "
         "axis (L100: 34 to 1 clash targets, 1.3% bond and 2.5% angle strain, a validity step on 124 of 126), stage 2 = stage 1 (L112, "
         "L118); C5 closed (addendum 2, L132: GLOBAL null, RIDGE harmful, the ORACLE ceilings); raw not run (L133).", []),
    11: ("Goal and ask", "Goal (two threads), ask (the study and one question), what will not be claimed. The direction line in the "
         "notes is the coordinator's L117 ruling verbatim, amended by L122 ('exact and complete'; the S26 additions checked, the S13 inputs "
         "cited without an S26 re-check; the Pauli mean weights off the slides).", []),
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
| "ADAPT ... selects no entangling operator on any of the 78 alpha = 1 targets under L-BFGS" / "declines them on all 78 targets" / L68's "stops with no operator selected on 78 of 78" | `s26/PROPOSAL_A.md` sections 3 and 5 (as first written); L68 | RETRACTED by lane Q in L75 on this lane's L73 flag and corrected by the PROPOSAL_A.md addendum: the growth is inert, not absent. The deck carries L75's final wording on slide 8 and in its notes: operators are appended on 60 (pool V) and 68 (pool L2) of the 78 alpha = 1 targets under L-BFGS, every one multi-qubit, and on 78 of 78 under Adam; they lower the free energy by at most 1.2e-4 nats under L-BFGS and 8.6e-4 under Adam, their angles stay at or below 0.018 rad under L-BFGS, and the state remains a product state to KL 4.1e-4; the 7-rotation state is already within 7.9e-4 nats of the Gibbs optimum. The endpoint numbers and the verdict REPLACE are unchanged. The ideal-ladder statements of L27 and L35 stand. | `s26/results/a1/*.json :: adapt/{V,L2}_{lbfgs,adam_best}_zrank/{stopped, sequence, trace, theta}; arms/adapt*_P21/kl_to_product` |
| the A1 basis | "built chain" in L68 | slide 8 names it `rmsd_q_synth` (the production projection of the WEIGHTED average over the 128 candidates, 3.2280 A for the deployed selector) and distinguishes it from the production top-75 arm 3.2148 | `s26/results/a1_stats.json :: means/fixed_zrank_it50/rmsd_q_synth`; `s26/results/q_mde_reference.json` |
| the steric reject +0.228 / +0.167 | L43 | quoted with the Type-M flag and "tail-carried, median +0.003" (L54) | `s26/results/ph_reject_report.json :: report/point_cloud/1e4` |
| the strain signal (Spearman +0.433, "how far the relaxation moves the chain predicts its error") | L53; the deck's 21:55 build (slide 7 and slide 10 notes) | RETRACTED as a framing by lane PH in L123 on the Adversary's L121: the pool's own disagreement (the top-75's pairwise CA-RMSD spread, native-free) predicts the error at Spearman +0.452 partial on n and Rg, fold CI [+0.280, +0.609], 5/5 folds; the relaxation's displacement tracks that spread at rho 0.756 and adds +0.082 given it (iid CI [-0.103, +0.259], permutation p 0.39). The deck's slide 7 and slide 10 notes carry that wording; the quartile means stand as a presentable form of the phenomenon; "the first native-free quantity above 0.4" was never on a slide | `s26/results/a_strain_vs_spread.json :: summary/{partial_n_rg/spread_mean, rho_moved_spread, partial_n_rg_spread/moved}` |
| L68's "null at the registered threshold" for A1 (the deck's 21:56 to 03:55 builds said "no change" / "null") | L68; `s26/PROPOSAL_A.md` section 3 | RETRACTED in scope by lane Q's L139 (seed 1: -0.045 / -0.052 A at 0.71x / 0.79x MDE, fold CIs excluding zero) and the Adversary's L140 (R11): slide 8 and its notes now say NOT MEASURED on either seed (underpowered at seed 0, Type-M at seed 1), quote both seeds with their MDE multiples, and carry the qualifier that the deployed fixed circuit moves 0.033 A between seeds (3.228 / 3.261) while the grown circuits are seed-stable (3.214 / 3.216); slide 5 carries the same qualifier; slide 4 states 3.2148 as the seed-free production number (Config quantum = False); the verdict REPLACE stands on seed-independent facts | `s26/results/a1s1_stats.json`, `a1_stats.json`; L139, L140 |
| the slide 11 line's "positive and complete" and "the Adversary has checked each" | L117 (the DRAFT line was replaced by L117 verbatim at the 02:1x build) | per the Adversary's L120 and the coordinator's L122: "exact and complete"; the checking is scoped to the S26 additions (A2 L45, A4 L47 with the L119 intervals, the product-state fact L70) and the S13 inputs are cited from their artefacts without an S26 re-check; the S13 Pauli mean weights 2.236 / 3.015 stay off every slide | `s26/results/q_var_boot.json :: results/ci` (L119); L120, L122 |
| PROPOSAL_C.md's pre-written C5 outcome under a future FINAL stamp (af05d987) | the 02:01 file | the Adversary's L120 MATERIAL item; lane P's 0a85323f rewrote it (live status, raw fold count, addendum 1); slide 10 and its notes are built from that file and read the live C5 checkpoint at build time (rows in `s26/results/p_c5.json`, `complete` not set: not a result) | `s26/results/p_c5.json`; `s26/models/p_ladder/raw_fold*_s0.pt` |

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

- Slide 8 is in its proposal form (verdict REPLACE, L68, accepted L69 subject to the Adversary's check of L68); the
  builder falls back to the PENDING placeholder only if `s26/PROPOSAL_A.md` is absent (`s26/pr_build_deck.py`, `builders[8]`).
- Slide 9's verdict is REPLACE, final (`s26/PROPOSAL_B.md` a88ea259; L117; Adversary L120 STANDS).
- Slide 10's verdict is KEEP WITH EDITS, final (`s26/PROPOSAL_C.md` 0a85323f; L117; L120 STANDS with the document fix applied);
  C5 completed at 03:12 and lane P's addendum 2 (03:54, L132) is applied: slide 10 and its notes state C5 as CLOSED (null to
  harmful; the ORACLE ceilings) with every number read from `s26/results/p_c5.json :: summary`; the raw rung is NOT RUN
  (4 of 5 folds trained, L133). Nothing on the deck is live or pending.
- Slide 11's direction line is the coordinator's ruling (L117) verbatim, with L122's amendment; nothing on the deck is DRAFT.
- Slide 8 is unchanged from the 21:56 build (L75 wording) apart from two additions to its notes: the L119 intervals and one
  sentence on A3 (L125: 124 of 126 distinct trained states under the entropy-matched raw score against 40 under the rank ladder;
  +0.0034 A, 0.04x MDE; `s26/results/a3_stats.json`, `a3_property.json`).
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
