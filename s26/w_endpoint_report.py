#!/usr/bin/env python
"""s26/w_endpoint_report.py -- the gated half of s26/PREREG_selfcopy_bound.md in one governed
job: `posterior` (Part B, on whatever channel-B models exist on disk; no training), then
`endpoint` (Parts A, B, C signed against the natives), then `report` (Part D, the bound).  A
driver so that `s26/w_selfcopy.py` is not edited.  Refuses to run before "PHASE 0 SIGNED OFF".

    python s26/jobrun.py --agent W --tag CPU --name w_endpoint_report --est-ram 0.5 -- python s26/w_endpoint_report.py
"""
from __future__ import annotations

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import w_selfcopy as W                               # noqa: E402


def main():
    W.require_signoff("w_endpoint_report")
    t0 = time.time()
    print("== posterior (Part B, existing models only; nothing is trained here)", flush=True)
    W.posterior()
    print("== endpoint (%.0f s)" % (time.time() - t0), flush=True)
    W.endpoint()
    print("== report (%.0f s)" % (time.time() - t0), flush=True)
    W.report()
    print("done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
