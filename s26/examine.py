#!/usr/bin/env python
"""s26/examine.py -- the repository's `make examine` (there is no Makefile; `examine.sh` and
`examine.bat` at the root call this).

Two steps, both cheap (no database, model or ESM bank is loaded):

  1. regenerate the module map, `s26/results/module_map.json`, by calling lane E's
     `s26/e_module_map.py` (imported, not duplicated);
  2. re-read every claimed number in `s26/results/claims.json` from the artefact it is
     claimed from (`s26/i_claim_check.py`) and report OK / MISMATCH / ABSENT.

    python s26/examine.py               # both
    python s26/examine.py --no-map      # claims only
    python s26/examine.py --no-claims   # map only
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-map", action="store_true")
    ap.add_argument("--no-claims", action="store_true")
    ap.add_argument("--claims", default=os.path.join(HERE, "results", "claims.json"))
    ap.add_argument("--search", action="store_true",
                    help="also run lane E's s26/e_claims.py, which finds WHERE each number "
                         "of EXAMINATION.md is cited and stored (slower: a tree-wide grep)")
    a = ap.parse_args(argv)
    rc = 0
    if not a.no_map:
        import e_module_map                                  # lane E's script, called as-is
        print("== module map (s26/e_module_map.py) ==")
        rc |= int(e_module_map.main() or 0)
    if not a.no_claims:
        import i_claim_check
        print("\n== claim ledger (s26/i_claim_check.py) ==")
        rc |= int(i_claim_check.main(["--claims", a.claims]) or 0)
    if a.search:
        import e_claims                                      # lane E's search, called as-is
        print("\n== claim search (s26/e_claims.py) ==")
        rc |= int(e_claims.main() or 0)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
