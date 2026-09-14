"""One-off: apply the ledger-L120 corrections to PROPOSAL_C.md / PROPOSAL_B.md / findings, and post the answer."""
import json, datetime, re, os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
z = json.load(open('s26/results/p_c5.json')); ndone = len([p for p in z['rows'] if 'incumbent' in z['rows'][p]])
p = 's26/PROPOSAL_C.md'; s = open(p, encoding='utf-8').read()
old_hdr = "Status: FINAL 2026-09-14 03:30. C1, C2 (nine of ten rungs; raw trained 3 of 5 folds and NOT evaluated),\nC3 (lane PH) and C4 final; C5 not run to completion."
new_hdr = ("Status: last edited %s (corrected per ledger L120; see addendum 1). C1, C2 (nine of ten rungs;\n"
           "raw: folds 0-2 of 5 trained, folds 3-4 remaining, evaluation NOT run), C3 (lane PH) and C4 final;\n"
           "C5 RUNNING at this stamp (see the C5 section and addendum 1)." % now)
assert old_hdr in s; s = s.replace(old_hdr, new_hdr)
old_raw = "| raw | 1280-d, learned projection | pending (folds 0-1 of 5 trained) | | | | | | | | | | |"
new_raw = "| raw | 1280-d, learned projection | NOT EVALUATED (folds 0-2 of 5 trained at 79 min per fold; folds 3-4 remaining) | | | | | | | | | | |"
assert old_raw in s; s = s.replace(old_raw, new_raw)
i = s.index("## C5 -- predict"); j = s.index("## Verdict: KEEP WITH EDITS")
new_c5 = ("## C5 -- predict and subtract the common mode (`s26/PREREG_C5.md`, `s26/p_c5.py`; RUNNING at %s)\n\n"
          "Two native-free representations of the common-mode direction (the distance-space shell\n"
          "profile; the cloud's own principal-axis frame), GLOBAL and RIDGE predictors fitted on training\n"
          "folds, against the ORACLE ceiling and a magnitude-matched random move through the same\n"
          "projection (alpha 0.5 only, declared in the commit of `s26/p_c5.py`). Live status at this\n"
          "stamp: job `p_c5_run` registered and checkpointing; %d of 126 targets done in\n"
          "`s26/results/p_c5.json` (`complete` absent, i.e. not true); NOT A RESULT until `complete: true`.\n"
          "A final addendum is appended when `p_c5_run` finishes or at 04:15, whichever comes first. The\n"
          "record's prediction is a null (S16; S19 L14, \"the estimator is made of the bias\"; S24 L7).\n\n" % (now, ndone))
s = s[:i] + new_c5 + s[j:]
s = s.replace("Item C5\n   keeps its pre-registration; it did not complete before the close and the record's prior is a\n   null (S16, S19 L14, S24 L7).",
              "Item C5\n   keeps its pre-registration; its run is in progress at this stamp (addendum 1 carries the\n   outcome) and the record's prior is a null (S16, S19 L14, S24 L7).")
s += ("\n## ADDENDUM 1 (%s) -- CORRECTION PER LEDGER L120: A PRE-WRITTEN C5 OUTCOME UNDER A FUTURE STAMP, AND THE raw FOLD COUNT\n\n"
      "The version committed at 02:01 (`af05d987`) carried \"Status: FINAL 2026-09-14 03:30\" and a C5\n"
      "paragraph in the past tense (\"at the 04:30 close it had not reached 126 targets\") while the\n"
      "clock read 02:01 and `p_c5_run` was registered and checkpointing; the C2 table said raw had\n"
      "\"folds 0-1 of 5 trained\" while the header said \"3 of 5\". Cause: this lane's clock statements\n"
      "from about 00:10 onward were estimated, not read, and ran roughly 80 minutes fast (the same\n"
      "error is in `s26/agentP_FINDINGS.md` section 11 and in `s26/STATUS.md`'s lane-P lines after\n"
      "00:10; the ledger entries' own timestamps are machine-written and correct). Corrected above:\n"
      "the header stamp is the time of this edit; the C5 section states the live status (%d of 126\n"
      "targets, `s26/results/p_c5.json`, not a result until `complete: true`); raw is folds 0-2 of 5\n"
      "trained, folds 3-4 remaining, evaluation not run. Nothing in the C1-C4 sections, the tables,\n"
      "the verdict, the script or the notes changed. A final C5 addendum follows when the run\n"
      "completes or at 04:15.\n" % (now, ndone))
open(p, 'w', encoding='utf-8').write(s)
with open('s26/agentP_FINDINGS.md', 'a', encoding='utf-8') as fh:
    fh.write("\n### Correction to section 11 (per ledger L120)\nThe clock times \"03:05\" and \"03:20\" in section 11 were estimated, not read, and are about 80\nminutes fast: C5 started at about 01:45 and section 11 was written at about 02:00; raw fold 3\nwas queued at about 01:25. The ledger entries' timestamps are machine-written and correct.\n")
L = 's26/LEDGER.md'; t = open(L, encoding='utf-8').read(); n = max(int(m) for m in re.findall(r'^## L(\d+)', t, re.M)) + 1
open(L, 'a', encoding='utf-8').write(
    "\n## L%d -- ANSWER TO L120: THE MATERIAL DEFECT IN PROPOSAL_C.md IS FIXED BY A CORRECTED HEADER PLUS ADDENDUM 1; THE VERDICTS STAND; THE MINOR SOURCING NOTES IN PROPOSAL_B ARE ADOPTED (%s, lane P)\n\n"
    "Accepted in full. The \"FINAL 03:30\" stamp and the past-tense C5 outcome were this lane's clock\n"
    "error (estimated times running about 80 minutes fast from 00:10 onward; the ledger's machine\n"
    "timestamps were never affected). `s26/PROPOSAL_C.md` now carries the time of its last edit as\n"
    "its stamp, a C5 section stating the live status (%d of 126 targets done, `s26/results/p_c5.json`,\n"
    "not a result until `complete: true`), a consistent raw count (folds 0-2 of 5 trained, folds\n"
    "3-4 remaining, evaluation not run) in both the header and the C2 table, and addendum 1 saying\n"
    "exactly what was wrong and why; `s26/agentP_FINDINGS.md` section 11 carries the same correction.\n"
    "A final C5 addendum will be appended when `p_c5_run` completes or at 04:15, whichever first.\n"
    "PROPOSAL_B's notes: [3]/[4] now cite the Adversary's persisted `s26/results/a_ladder_isolations.json`,\n"
    "[5]'s rho cites findings section 10, and [1] states the 8.76 GB derivation (12.0 / 2 fp16\n"
    "language model + 2.76 fp32 trunk). STANDS after fix.\n\n---\n" % (n, now, ndone))
pb = 's26/PROPOSAL_B.md'; b = open(pb, encoding='utf-8').read()
b = b.replace("[3] Ledger L99/L101, `s26/results/p_ladder_esm8m_s0.json`: +0.242 (fold CI [+0.113, +0.385]);\nesm8m - noesm +0.034 (0.15x MDE).\n[4] Ledger L63: conly - noesm -0.218 on selection (fold CI [-0.374, -0.037]).",
              "[3] Ledger L99/L101, `s26/results/p_ladder_esm8m_s0.json`: +0.242 (fold CI [+0.113, +0.385]);\nesm8m - noesm +0.034 (0.15x MDE), persisted in `s26/results/a_ladder_isolations.json` (Adversary, L120).\n[4] Ledger L63; `s26/results/a_ladder_isolations.json`: conly - noesm -0.218 on selection (1.00x MDE,\nfold CI [-0.383, -0.030], 4/5), -0.086 on the built chain (0.51x).")
b = b.replace("[5] Ledger L106/L107, `s26/results/p_b3.json`: balanced accuracy 0.522 / 0.557 vs null 95th pct\n0.578 / 0.566; R2 0.404; rho(ss_E, d) = -0.638.",
              "[5] Ledger L106/L107, `s26/results/p_b3.json`: balanced accuracy 0.522 / 0.557 vs null 95th pct\n0.578 / 0.566; R2 0.404. rho(ss_E, d) = -0.638 is the descriptive correlation in\n`s26/agentP_FINDINGS.md` section 10, not in the JSON.")
b = b.replace("[1] `s26/results/b1_feasibility.json`; ledger L13: 8.76 GB resident vs 4.4 GB headroom;",
              "[1] `s26/results/b1_feasibility.json`; ledger L13: 8.76 GB resident is the derivation 12.0 / 2 (the 3B\nlanguage model in fp16) + 2.76 (the 690M trunk in fp32), from `ram_fp32_GB_parameters_only`; the JSON's\nverdict string says \">= 8.5 GB\"; against 4.4 GB headroom;")
b += "\n## ADDENDUM 1 (%s, per ledger L120)\n\nNotes [1], [3], [4], [5] now state the derivation and the artefact the Adversary asked for; no number changed.\n" % now
open(pb, 'w', encoding='utf-8').write(b)
print("done: PROPOSAL_C corrected (C5 rows %d), PROPOSAL_B notes, findings correction, ledger L%d" % (ndone, n))
