"""s23/c1_explore_t0_analyze.py -- analysis for the DATED, NOT-PRE-REGISTERED T=0 addendum.

Same aggregation as `c1_analyze.py`, applied to `c1_explore_t0_raw.json`. Reported separately
and labelled EXPLORATORY throughout -- never merged into the pre-registered primary/secondary
verdicts.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s23 import qc_lib as QL          # noqa: E402
from s23 import c1_analyze as CA      # noqa: E402


def run():
    raw = json.load(open(os.path.join(QL.RESULTS, "c1_explore_t0_raw.json")))
    assert raw["complete"] and raw["n_expected"] == 126
    d = CA.analyze_T(raw, "0.0")
    out = {"config": raw["config"], "by_T": {"0.0": d}, "complete": True,
          "EXPLORATORY": True, "dated": "2026-09-08", "pre_registered": False}
    QL.save_json("c1_explore_t0_analysis.json", out)
    CA.report(out)
    return out


if __name__ == "__main__":
    run()
