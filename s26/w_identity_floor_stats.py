#!/usr/bin/env python
"""s26/w_identity_floor_stats.py -- the two paired contrasts registered in
s26/PREREG_identity_floor.md, from s26/results/w_selfcopy_floor.json (Part E of the self-copy
bound, already measured).  ORACLE DIAGNOSTIC on both sides.  Reads no native: the per-row RMSDs
are in the artefact.

    python s26/jobrun.py --agent W --tag CPU --name w_identity_floor_stats --est-ram 0.2 -- python s26/w_identity_floor_stats.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s24 import stats_lib as ST                      # noqa: E402

SRC = os.path.join(HERE, "results", "w_selfcopy_floor.json")
OUT = os.path.join(HERE, "results", "w_identity_floor.json")


def main():
    z = json.load(open(SRC))
    assert z.get("complete"), "the Part E artefact is not complete"
    rows = z["rows"]
    per = {}
    for r in rows:
        d = per.setdefault(r["pdb"], {"fold": r["fold"], "n": r["n"], "copy": [], "pool_best": r["pool_best_rr"],
                                      "pool_mean": r["pool_mean_rr"], "partners": 0})
        d["copy"].append(r["cross_deposit_rmsd"]); d["partners"] += 1
    pdbs = sorted(per)
    copy = np.array([np.mean(per[p]["copy"]) for p in pdbs])
    best = np.array([per[p]["pool_best"] for p in pdbs])
    mean = np.array([per[p]["pool_mean"] for p in pdbs])
    folds = ST.pinned_folds(pdbs)
    out = {"label": "identity floor: paired contrasts (ORACLE on both sides)", "n_targets": len(pdbs),
           "n_pairs": len(rows), "per_target": {p: {**per[p], "copy_mean": float(np.mean(per[p]["copy"]))} for p in pdbs},
           "median_over_targets": float(np.median(copy)), "median_over_pairs": float(z["summary"]["median"]),
           "frac_targets_below_1.0": float((copy < 1.0).mean()), "frac_targets_below_1.5": float((copy < 1.5).mean()),
           "source": os.path.relpath(SRC, ROOT), "prereg": "s26/PREREG_identity_floor.md"}
    for name, other in (("copy_minus_pool_best", best), ("copy_minus_pool_mean", mean)):
        r = ST.compare(copy, other, folds, names=pdbs, label="%s (ORACLE both sides; single window vs native), n=%d targets" % (name, len(pdbs)))
        out[name] = r
        print(ST.fmt(r))
    print("median over targets %.3f (over the 22 pairs %.3f); below 1.0 A: %d/%d targets; below 1.5 A: %d/%d"
          % (out["median_over_targets"], out["median_over_pairs"], int((copy < 1.0).sum()), len(copy), int((copy < 1.5).sum()), len(copy)))
    ST.save_atomic(OUT, out, module_file=__file__)
    print("  ->", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
