#!/usr/bin/env python
"""s27/s28_D_leakgrep.py -- lane D: the leakage grep over a lane's scripts.

Greps every hit of the native-bearing tokens (`nat_ca`, `oracle_rr`, `\\brr\\b`, `native`,
`rmsd`, `Dnat`, `oracle`) in the named files, prints them with line numbers so each can be
traced by hand, and writes the hit list to `s27/results/s28_D_leakgrep_<tag>.json`.  A hit is
not a leak; an untraced hit is.  The trace (which function the hit lives in, and whether that
function's output reaches a deployable arm) is written into the ledger entry, not here.

Usage:  python s27/s28_D_leakgrep.py --tag A s27/s28_A_amp.py [more files]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

TOKENS = {
    "nat_ca": r"nat_ca",
    "oracle_rr": r"oracle_rr",
    "rr": r"\brr\b",
    "native": r"native",
    "rmsd": r"rmsd",
    "Dnat": r"Dnat",
    "oracle": r"oracle",
}


def _enclosing_def(lines, i):
    """Name of the nearest `def`/`class` above line i (0-based), or '<module>'."""
    for j in range(i, -1, -1):
        m = re.match(r"^\s*(def|class)\s+([A-Za-z_][A-Za-z0-9_]*)", lines[j])
        if m:
            return m.group(2)
    return "<module>"


def grep(path):
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    hits = []
    for i, line in enumerate(lines):
        for name, pat in TOKENS.items():
            if re.search(pat, line):
                hits.append(dict(file=os.path.relpath(path, ROOT), line=i + 1, token=name,
                                 fn=_enclosing_def(lines, i), text=line.strip()[:160]))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    out = {"tag": a.tag, "files": [], "hits": []}
    for f in a.files:
        p = f if os.path.isabs(f) else os.path.join(ROOT, f)
        if not os.path.exists(p):
            print("MISSING", f)
            out["files"].append(dict(file=f, exists=False))
            continue
        h = grep(p)
        out["files"].append(dict(file=f, exists=True, n_hits=len(h)))
        out["hits"].extend(h)
        print(f"== {f}: {len(h)} hits")
        for x in h:
            print(f"  {x['line']:5d}  [{x['token']:9s}]  in {x['fn']:30s}  {x['text']}")
    by_fn = {}
    for x in out["hits"]:
        by_fn.setdefault((x["file"], x["fn"]), []).append(x["token"])
    out["by_function"] = [dict(file=k[0], fn=k[1], tokens=sorted(set(v)), n=len(v))
                          for k, v in sorted(by_fn.items())]
    print("\nfunctions carrying hits (trace each):")
    for r in out["by_function"]:
        print(f"  {r['file']} :: {r['fn']}  {r['tokens']}  ({r['n']})")
    op = os.path.join(HERE, "results", f"s28_D_leakgrep_{a.tag}.json")
    os.makedirs(os.path.dirname(op), exist_ok=True)
    with open(op, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print("wrote", op)


if __name__ == "__main__":
    main()
