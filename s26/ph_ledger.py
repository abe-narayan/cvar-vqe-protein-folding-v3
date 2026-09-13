"""s26/ph_ledger.py -- collision-safe ledger append for lane PH.

    python s26/ph_ledger.py --title "TITLE (2026-09-13, PH)" --body s26/results/_entry.md [--status "..."]

Re-reads `s26/LEDGER.md` immediately before appending (L16b rule 4), takes max(L<n>) + 1, writes
`## L<n> -- TITLE` followed by the body file verbatim and a `---` line.  Prints the number taken.
Never edits an existing entry.  Optionally appends one STATUS line under `## PH`.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LED = os.path.join(HERE, "LEDGER.md")
STATUS = os.path.join(HERE, "STATUS.md")


def next_number() -> int:
    with open(LED, encoding="utf-8") as fh:
        txt = fh.read()
    nums = [int(m) for m in re.findall(r"^## L(\d+)", txt, flags=re.M)]
    return max(nums) + 1


def append_entry(title: str, body: str) -> int:
    n = next_number()
    with open(LED, encoding="utf-8") as fh:
        txt = fh.read()
    entry = f"## L{n} -- {title}\n\n{body.rstrip()}\n\n---\n"
    with open(LED, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(txt.rstrip("\n") + "\n\n" + entry)
    return n


def append_status(lines) -> None:
    with open(STATUS, encoding="utf-8") as fh:
        rows = fh.read().split("\n")
    # insert after the last line of the ## PH block
    start = next(i for i, l in enumerate(rows) if l.startswith("## PH"))
    end = start + 1
    while end < len(rows) and not rows[end].startswith("## "):
        end += 1
    while end > start + 1 and rows[end - 1].strip() == "":
        end -= 1
    rows[end:end] = list(lines)
    with open(STATUS, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--title", required=True)
    ap.add_argument("--body", required=True)
    ap.add_argument("--status", action="append", default=[])
    a = ap.parse_args()
    with open(a.body, encoding="utf-8") as fh:
        body = fh.read()
    n = append_entry(a.title, body)
    print(f"L{n}")
    if a.status:
        now = time.strftime("%Y-%m-%d %H:%M")
        append_status([f"- {now} {s}" for s in a.status])
