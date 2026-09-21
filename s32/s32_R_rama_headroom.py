#!/usr/bin/env python
"""S32 LANE R -- WHY the chiral Ramachandran criterion has no in-band skill.

EXPLORATORY, logged as R-20 in `s32/MULTIPLICITY.md`.

I WROTE THE WRONG MECHANISM FIRST AND MY OWN DATA REFUTED IT.  The first version of this
check asserted: *"every branch that survives the lam = 0.3 objective is already
Ramachandran-plausible, so the dense log-density has almost nothing left to separate."*
**That is false.**  The measurement below shows the branches differ ENORMOUSLY in
Ramachandran plausibility -- the within-target spread of the positive-phi rate is 0.63, and
124 of 124 targets carry at least one branch above the 17.5% rate of a completely
unconstrained fit.  The criterion has plenty to discriminate.  It simply carries no
information about which branch is nearer the native.

SO THE ACTUAL MECHANISM IS THE STRONGER ONE:

    Ramachandran plausibility and native proximity are ORTHOGONAL among the branches this
    projection admits.  `core/project.py`'s docstring is RIGHT that a CA trace admits
    solutions "one Ramachandran-plausible and one not" -- the branches really do differ that
    way -- but **the plausible one is not the one closer to the native.**

That is lane R re-deriving S9-2's conclusion ("a torsion prior buys physical validity for
free, and no accuracy") one level down, at the level of individual branches, from a
different measurement.

REFERENCE RATES (S7 / S9-2, measured, not assumed): the real positive-phi rate over library
windows is 5.40% (band 2.90-5.66%); an unconstrained lam = 0 fit sits at 17.5%; the shipped
`ramah` prior pulls the emitted structure to 5.5%.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
UNCONSTRAINED = 0.175          # S9-2's measured lam=0 positive-phi rate
REAL = 0.0540                  # S9-2's measured real-library rate


def main():
    rows = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_R_branches_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    rows[r["pdb"]] = r
    P = list(rows)
    pf = [np.asarray(rows[p]["posphi_frac"], float) for p in P]
    rn = [np.asarray(rows[p]["rama_nlp"], float) for p in P]
    spread = np.array([v.max() - v.min() for v in pf])
    worst = np.array([v.max() for v in pf])
    bad = np.array([float((v > UNCONSTRAINED).mean()) for v in pf])
    rspread = np.array([v.max() - v.min() for v in rn])

    out = {
        "n": len(P), "INCOMPLETE": len(P) < 126,
        "reference_rates": {"real_library": REAL, "unconstrained_lam0": UNCONSTRAINED,
                            "shipped_ramah_emitted": 0.055,
                            "source": "S9-2, measured"},
        "posphi_mean_over_branches": float(np.mean([v.mean() for v in pf])),
        "posphi_within_target_spread_mean": float(spread.mean()),
        "posphi_within_target_spread_median": float(np.median(spread)),
        "posphi_worst_branch_mean": float(worst.mean()),
        "frac_branches_above_unconstrained_rate": float(bad.mean()),
        "n_targets_with_any_branch_above_unconstrained": int((bad > 0).sum()),
        "rama_nlp_within_target_spread_mean": float(rspread.mean()),
        "rama_nlp_within_target_spread_median": float(np.median(rspread)),
        "inband_rho_rama_nlp_see": "s32_R_analysis.json :: R3a_inband_rho.rama_nlp",
        "MECHANISM": ("The branch set is NOT uniformly plausible -- the within-target "
                      "positive-phi spread is %.3f and %d of %d targets carry a branch "
                      "above the 17.5%% unconstrained rate -- so the criterion has plenty "
                      "to discriminate. It has no in-band skill anyway. Ramachandran "
                      "plausibility and native proximity are ORTHOGONAL among the branches "
                      "this projection admits."
                      % (spread.mean(), int((bad > 0).sum()), len(P))),
        "RETRACTED_FIRST_MECHANISM": ("An earlier draft of this check claimed the branches "
                                      "were all already plausible so there was nothing to "
                                      "separate. The data above refutes that; it is "
                                      "recorded rather than deleted."),
    }
    ST.save_atomic(os.path.join(RESULTS, "s32_R_rama_headroom.json"), out,
                   module_file=__file__)

    print("n = %d%s\n" % (out["n"], "  (INCOMPLETE)" if out["INCOMPLETE"] else ""))
    print("IS THERE RAMACHANDRAN VARIATION LEFT AMONG THE lam=0.3 BRANCHES?  YES, A LOT.")
    print("  reference: real library %.1f%%, unconstrained lam=0 fit %.1f%%, shipped %.1f%%"
          % (100 * REAL, 100 * UNCONSTRAINED, 5.5))
    print("  mean positive-phi over all branches           %.4f  (%.1f%%)"
          % (out["posphi_mean_over_branches"], 100 * out["posphi_mean_over_branches"]))
    print("  within-target spread (max-min)  mean %.4f  median %.4f"
          % (out["posphi_within_target_spread_mean"],
             out["posphi_within_target_spread_median"]))
    print("  worst branch per target, mean                 %.4f"
          % out["posphi_worst_branch_mean"])
    print("  branches above the 17.5%% unconstrained rate   %.1f%%"
          % (100 * out["frac_branches_above_unconstrained_rate"]))
    print("  targets with ANY branch above 17.5%%           %d of %d"
          % (out["n_targets_with_any_branch_above_unconstrained"], out["n"]))
    print("  rama_nlp within-target spread  mean %.4f nats  median %.4f"
          % (out["rama_nlp_within_target_spread_mean"],
             out["rama_nlp_within_target_spread_median"]))
    print("\n" + out["MECHANISM"])
    print("\nwrote", os.path.join(RESULTS, "s32_R_rama_headroom.json"))


if __name__ == "__main__":
    main()
