# S26 LANE BRIEF -- PR (Presentation lane)

You are lane **PR** of Sprint 26. Read `s26/LANE_CONTRACT.md` in full first; it binds you. Then
`docs/STATE_BRIEF_2026-09-12.md` in full, `s26/EXAMINATION.md` (the claim ledger in it is the
source of every number on slides 2 to 7), `s26/LEDGER.md`, every `s26/PROPOSAL_*.md` that exists,
`s26/TRAINABILITY_PAPER_OUTLINE.md` if it exists, `results/summary/professor_brief.md`, and
`docs/CONDENSED_REPORT.md`. Then the S25 quantum findings (`s25/QUANTUM.md`, `s25/agentQ_FINDINGS.md`)
because slides 5 and 6 rest on them.

## 1. The situation you inherit (ledger L2)

The campaign prompt names `vqe_research_overview.pptx`, an 11-slide dark-theme deck: title; six
"what I built" slides (why peptides are hard; the pipeline; results with structure overlays; the
verified quantum component; the barren-plateau result; where the remaining accuracy lives); three
"where I want to go" slides (ADAPT-VQE; AlphaFold gaps; learned model plus physics); goal and ask.
**That file does not exist on this machine or in git history** (ledger L2). You build it from that
structure with python-pptx (`import pptx`, version 1.0.2 is installed; matplotlib 3.11 for figures),
at the repository root as `vqe_research_overview.pptx`, and `s26/PRESENTATION_CHANGES.md` records
every slide's content and the artefact behind every number, so the edits map onto the real file if
the presenter supplies it later. If the file appears in the repository root before you start (the
user may have added it), read it with python-pptx, keep its slides 1 to 7 and 11 as the base, and
edit rather than rebuild; say which case you were in.

## 2. What the deck must satisfy (campaign prompt Part 6)

- Slides 8, 9, 10 carry the edited or replaced proposals exactly as `s26/PROPOSAL_A.md`,
  `s26/PROPOSAL_B.md` (or `s26/PROPOSAL_B_REPLACEMENT.md`) and `s26/PROPOSAL_C.md` state them. Each
  must be sayable by a sixteen-year-old in under two minutes (about 250 spoken words; put the
  spoken text in the speaker notes and count the words). Every number in the notes carries an
  artefact path.
- Slides 2 to 7: every number verified against the Phase 0 claim ledger in `s26/EXAMINATION.md`.
  A number with no artefact does not go on a slide. Stale numbers are fixed and the fix is listed
  in `PRESENTATION_CHANGES.md` with the artefact.
- Figures: if lane Q produced a clean A2 figure (DLA dimension versus n) or A4 figure (gradient
  variance for grown versus fixed ansatz), it goes on slide 6 or 8; if lane P produced a C2 ladder
  figure, it goes on slide 10. Regenerate every figure yourself from the results JSON at 190 dpi on
  a white background (`fig.savefig(path, dpi=190, facecolor="white")`), so it sits on a white
  plate inside the dark theme. Figures under `s26/figures/`. Never copy a number by hand into a
  figure; read it from the artefact.
- Slide 11 speaker notes get one plain line saying which of the three directions the evidence
  favours and why, so the presenter can answer "which would you do first?" without hesitating.
  Take that line from the coordinator's ledger entry on the proposal verdicts; if it is not posted
  yet, draft it from the PROPOSAL files and mark it DRAFT in `PRESENTATION_CHANGES.md`.
- Style: no em dashes anywhere (U+2014; also avoid U+2013 as a substitute). No "genuinely",
  "honestly", "leverage", "robust", "delve", "underscore". Short sentences. Formal but human. The
  presenter should sound like a sixteen-year-old who did the work. Never claim a quantum advantage.
  Never call a small gradient a barren plateau; slide 6 says what the width sweep measured
  (log2 Var per qubit -0.649 / -0.252 / -0.047 / -0.311 / -0.243 across the cells it was run on)
  and what it does not say.
- Dark theme: background near-black (for example RGB 18,18,24), text off-white, one accent colour,
  figures on a white plate. Title, body, notes on every slide. Consistent fonts.

## 3. Structure overlays on slide 4

Slide 4 shows results with structure overlays. Draw CA traces (predicted versus native, after
Kabsch superposition) for T030 (the 0.182 A target) and one hard target above 5 A, from the
production cache `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` (`ca` is the built chain) and
the native from `s8/generate_univ/<pdb>.npz` `nat_ca`. Reading a dev-target native for a FIGURE is
an ORACLE use and is allowed; label the figure ORACLE-superposed in the notes. Never touch a
benchmark target. Use `core.numerics` or `s12/instrument.py` for the superposition (do not write
your own Kabsch).

## 4. Process

1. Wait for nothing you do not need: start with slides 1 to 7 and 11 as soon as `s26/EXAMINATION.md`
   exists, then do slides 8 to 10 as each PROPOSAL file lands (the coordinator will message you
   when one does; check `ls s26/PROPOSAL_*.md` yourself each hour as well).
2. Build with a script `s26/pr_build_deck.py` that reads the artefacts and writes the deck, so the
   deck is regenerable; keep the spoken text in a companion `s26/pr_notes.md` that the script
   reads. Run the build script directly (it is light) but if it loads any results JSON larger than
   50 MB, go through jobrun.
3. Verify the output: reopen the .pptx with python-pptx, count slides (11), dump every text frame
   and every notes frame, grep the dump for the banned words and for U+2014/U+2013, count notes
   words per proposal slide, and paste the checks into `s26/PRESENTATION_CHANGES.md`.
4. Write `s26/agentPR_FINDINGS.md` (what you built, what you could not source, what you changed
   and why) and hourly status lines under `## PR (Presentation lane)` in `s26/STATUS.md`.
5. Commit your own files (`vqe_research_overview.pptx`, `s26/pr_*.py`, `s26/pr_notes.md`,
   `s26/figures/*.png`, `s26/PRESENTATION_CHANGES.md`, `s26/agentPR_FINDINGS.md`) with the trailer
   lines from the contract. `*.png` under `s26/figures/` may need a `.gitignore` whitelist; check
   `git check-ignore -v` and add `!s26/figures/*.png` with a ledger note if needed.

Finish your turn with: what is built, which slides still wait on which PROPOSAL file, and any
number you could not source. Do not invent a number for a slide that has no artefact; leave the
placeholder text "NUMBER PENDING: <artefact needed>" and list it.
