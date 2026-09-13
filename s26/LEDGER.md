# SPRINT 26 -- DECISION LEDGER

Append-only. One entry per finding or decision. Nobody edits another lane's entry; a
retraction is a new entry naming the entry it retracts. Every number carries an artefact path.

---

## L0 -- THE GOVERNOR IS RUNNING; THE BASELINE LOAD LEAVES ABOUT 4 GB (2026-09-12, coordinator)

`s26/governor.py` started at 23:50:41 (pid 36196), commit `a4db170c`. Self-test: a 0.3 GB job
registered through `s26/jobrun.py`, ran, and reported peak RSS 0.328 GB
(`s26/jobs_done/selftest_alloc.json`). Every action is logged to `s26/governor.log`.

Measured before any campaign process existed (`python s26/governor.py --once`):

    RAM total 16.75 GB (15.6 GiB)   used 11.16 GB = 66.6%   available 5.59 GB   CPU 10.5%

The 66.6% is the user's own environment (several VS Code and claude sessions, Defender;
`psutil` process table, 2026-09-12 23:47). The campaign's 93% ceiling therefore leaves about
4.4 GB for all lanes together. The "90 to 93% band" instruction is read as: fill up to, never
above, 93%; the governor launches queued work only when the box has sat under 88% for 60 s.

---

## L1 -- OPERATIONAL ITEM 1 CLOSED: `s24/cache_amber/*.npz` IS TRACKED (2026-09-12, coordinator)

Commit `6b494780`, the first on branch `s26`: `.gitignore` whitelist `!s24/cache_amber/*.npz`
placed after the blanket `*.npz` rule; 126 files, 1.5 MB, `git ls-files s24/cache_amber`
returns 126. S25 L13's single-point-of-failure on the 63,000 ff14SB/GBn2 single points is
closed.

---

## L2 -- THE PRESENTATION FILE DOES NOT EXIST ON THIS MACHINE (2026-09-12, coordinator)

`vqe_research_overview.pptx` is named in the campaign prompt as the 11-slide deck to tune. It
is not in the repository, not in git history (`git log --all --diff-filter=A -- '*.pptx'`
returns nothing), and not anywhere under `C:\Users\abena` to depth 7 excluding `AppData`
(`find` over Desktop, Documents, Downloads, OneDrive and the profile root; the only `.pptx`
files are six unrelated lecture decks in Downloads and python-pptx's template).

Decision: Phase 4 will build the deck from the structure the prompt describes (title; six
"what I built" slides; three proposal slides; goal and ask; dark theme) using python-pptx,
and `s26/PRESENTATION_CHANGES.md` will record every slide's content and its artefact so the
edits map onto the real file if the user supplies it. No number will be invented for a slide
that has no artefact.

---

## L3 -- `TEST_RUN_RESULT_PLACEHOLDER` IS NOT IN THE TREE (2026-09-12, coordinator)

The prompt says the state brief contains the literal string `TEST_RUN_RESULT_PLACEHOLDER`.
`grep -rn` over every `.md`, `.py` and `.txt` finds no occurrence. The brief on disk
(`docs/STATE_BRIEF_2026-09-12.md` section 6.1) already carries a live run: "355 passed, 13
skipped, 0 failed, 0 errors, exit 0" dated 2026-09-12. Operational item 5 is therefore "replace
that line with the S26 run under the governor, with per-file counts and the skip reasons", and
the Infrastructure lane owns it.

---

## L4 -- BRANCH AND READING ORDER (2026-09-12, coordinator)

Branch `s26` created off `ae86a124`. Coordinator read, in this order and in full: the state
brief, README, ARCHITECTURE, professor brief, S25 LEDGER (1,270 lines), S25 QUANTUM.md (857),
S25 agentQ_FINDINGS, S20 agentD_FINDINGS (1,006), docs/FINDINGS.md (all 9,081 lines, corrections
ledger first), docs/CONDENSED_REPORT, the S24, S23, S22, S21, S20, S19, S18, S17, S16, S15 and
S14 ledgers, the S13 and S12 dossiers, S13 qarch_FINDINGS part 1, the S25 BRIEF and PREREG_Q,
s12/instrument.py, s24/stats_lib.py, core/pipeline.py, and the deployed-selector region of
core/quantum.py. The remaining code read (core modules, root modules, tests, verify,
resultslab) is assigned to the Examiner and Infrastructure lanes, whose module map is the
Phase 0 deliverable.

---
## L5 -- PHASE-GATE READING FOR THE CIS-PEPTIDE CENSUS, AND LANE PH SPAWNED (2026-09-13, coordinator)

The phase gate forbids any endpoint experiment (anything reading an RMSD to a native) before
"PHASE 0 SIGNED OFF". The cis-peptide census reads omega angles of the 126 dev natives and CA-CA
distances of pool windows; it computes no RMSD, selects nothing and changes nothing. It is an ORACLE
DIAGNOSTIC and is allowed before the gate. The projection-floor cost on cis targets (an RMSD) waits
for the gate. Lane PH (Physics) spawned at 00:20 with brief `s26/briefs/PH.md`: AMBER as a steric
reject filter, the cis-peptide gap, the steric singularity for Part IV, and the C3 matched-random
control (stage 1 on the production cache, stage 2 on the best C2 rung when lane P delivers it).
Lanes active: E, I, Q, P, PH (`s26/lanes.json`).

---

