"""s26/p_ledger_rung.py -- post one rung's report to s26/LEDGER.md as the next free L<n>.

Re-reads the ledger tail immediately before appending (five lanes append concurrently; LEDGER
L16b rule 4).  The block is `s26/p_rung_report.py`'s text verbatim (ST.fmt for arm / cloud /
sel, the gamma-equivalent with the cos caveat), preceded by the one-line reading.

    python s26/p_ledger_rung.py --rung noesm --title "..." [--seed 0] [--note "..."]
"""
from __future__ import annotations

import argparse
import datetime
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

from s26 import p_rung_report as R                   # noqa: E402

LEDGER = os.path.join(HERE, "LEDGER.md")


def next_number(text):
    nums = [int(m) for m in re.findall(r"^## L(\d+)", text, re.M)]
    return max(nums) + 1


def post(rung, title, seed=0, note=""):
    rep = R.report(rung, seed)
    text = open(LEDGER, encoding="utf-8").read()
    n = next_number(text)
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    body = "\n## L%d -- LANE P, C2 RUNG %s: %s (%s, lane P)\n\n" % (n, rung.upper(), title, stamp)
    body += ("Artefacts: `s26/results/p_ladder_%s_s%d.json` (126 rows, complete), "
             "`s26/results/p_ladder_report_%s_s%d.json`; anchor `s26/results/p_ladder_shipped_s%d.json`. "
             "Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior "
             "in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame "
             "average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior "
             "through the same path. Negative = the rung is better.\n\n" % (rung, seed, rung, seed, seed))
    if note:
        body += note.rstrip() + "\n\n"
    body += "```\n" + rep["text"] + "\n```\n\n---\n"
    # re-read right before the append so a concurrent entry cannot be overwritten or double-numbered
    text2 = open(LEDGER, encoding="utf-8").read()
    n2 = next_number(text2)
    if n2 != n:
        body = body.replace("## L%d --" % n, "## L%d --" % n2, 1); n = n2
    with open(LEDGER, "a", encoding="utf-8") as fh:
        fh.write(body)
    print("posted L%d for rung %s" % (n, rung))
    return n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rung", required=True); ap.add_argument("--title", required=True)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--note", default="")
    a = ap.parse_args(); post(a.rung, a.title, a.seed, a.note)
