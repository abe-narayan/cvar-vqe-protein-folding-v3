# SPRINT 30 — STATUS

## Setup notes (before any science)

**Two corrections to the charter's reading list, made by checking rather than assuming**
(the habit S29 paid for twice — see project memory `findings-prose-is-not-evidence-of-code`):

- charter says `s27/REPORT_S29.md` -> the S29 report is at **`s29/REPORT_S29.md`** (1,619 lines)
- charter says `s27/LEDGER.md through the S29 close` -> `s27/LEDGER.md` is **S28's** ledger
  (4,215 lines, ends at S28-L50). S29's ledger is **`s29/LEDGER.md`** (5,823 lines, S29-L0..L57).
  Both are read; they are different documents.

`s27/RETRACTIONS_S28.md` exists and is 11 lines. There is no separate S29 retractions file --
S29's withdrawals live in `s29/REPORT_S29.md` section 14(e) (18 entries, 11 of them the
coordinator's) and are annotated in place in the ledger.

**Section 7's meter is already built.** `s29/s29_D_cost_audit.py` (767 lines) implements the
cost-RMSD meter with the six diagnostics S29 specified. S30 VERIFIES and EXTENDS it rather than
rebuilding it -- the charter's "build this first" is satisfied by making it run and adding the
diagnostics section 7 lists that S29's version lacks (per-target distributions, FAIL18/108 split
as standard output).

**One inherited defect is still open and will bite this sprint immediately:**
`s26/jobrun.py` `CPU_START = 85.0` refuses to launch while smoothed CPU exceeds 85%, while
`s26/governor.py` `CPU_CEILING = 101.0` deliberately holds the box near 100%. Under the charter's
own 94-95% utilisation target the launch gate is closed most of the time. See project memory
`paired-thresholds-move-together`.
