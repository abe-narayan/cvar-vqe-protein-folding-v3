#!/usr/bin/env python
"""S32 LANE R -- both directions of `typicality`, persisted, because the name inverted a
cross-lane claim.

`typicality` in `s32_R_branches*.jsonl` is a DISTANCE: the mean Ca-RMSD of a branch to the
production top-75 members.  **Lower means MORE typical** -- the sense `s8/inband.py` and the
consensus medoid use.  `s32/s32_V_R_adversary.py:56` sets `MAXIMISE = {"typicality"}` and so
selected the LEAST typical branch, and the draft report glossed that arm as "the most typical
branch is much worse".  It is the most ATYPICAL branch that is much worse, which is the
EXPECTED direction and supports consensus-as-outlier-avoidance rather than contradicting it.

This is a contract rule 7 defect -- a claim spanning two lanes, audited by neither -- and it
is re-derived here in ONE script from ONE raw artefact so the two directions can never again
be quoted as one number.  The naming is lane R's fault; the column is not renamed because
rows are already written and other lanes read them.

CONTRACT RULE 5: this check can fail.  If the two directions agreed, the sign would not
matter and the script would say so.  They differ by 0.35 A.

BASIS: built-chain Ca-RMSD, tuning126.  Every arm is a SELECTION over branches built in one
job from one cloud, paired to that job's own PRODUCTION row.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST                                          # noqa: E402
from s12 import instrument as I                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")


def tied_mean(score, outcome, sign):
    """Mean outcome over the TIED argmin set of `sign * score` (S12's standing lesson:
    np.argmin on a tied criterion reads the array order, which was the ORACLE sort)."""
    s = np.asarray(score, float) * sign
    o = np.asarray(outcome, float)
    m = s.min()
    t = np.isclose(s, m, atol=1e-12, rtol=1e-9)
    return float(o[t].mean()), int(t.sum())


def main():
    rows = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_R_branches_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    rows[r["pdb"]] = r
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in rows]
    n = len(pdbs)
    F = np.array([rows[p]["fold"] for p in pdbs], int)
    prod = np.array([rows[p]["prod_chain"] for p in pdbs])
    bmean = np.array([float(np.mean(rows[p]["rmsd_nat"])) for p in pdbs])

    out = {"n": n, "INCOMPLETE": n < 126,
           "definition": ("typicality = mean CA-RMSD of the branch to the production "
                          "top-75 members; a DISTANCE, so LOWER = MORE TYPICAL"),
           "why": ("s32_V_R_adversary.py:56 MAXIMISE={'typicality'} selects the LEAST "
                   "typical branch; the draft report glossed it as the most typical"),
           "PROD_chain_mean": float(prod.mean()),
           "branch_set_mean_zero_skill_reference": float(bmean.mean()),
           "arms": {}}

    for sign, name, gloss in ((+1, "ARGMIN_typicality", "the MOST typical branch"),
                              (-1, "ARGMAX_typicality", "the LEAST typical branch")):
        vals, ties = zip(*[tied_mean(rows[p]["typicality"], rows[p]["rmsd_nat"], sign)
                           for p in pdbs])
        v = np.array(vals, float)
        out["arms"][name] = {
            "gloss": gloss, "sign": sign,
            "chain_mean": float(v.mean()),
            "chain_se": float(v.std(ddof=1) / math.sqrt(n)),
            "mean_tie_size": float(np.mean(ties)),
            "vs_prod": ST.compare(v, prod, F, names=pdbs, label="%s - PROD" % name),
        }
    a, b = out["arms"]["ARGMIN_typicality"], out["arms"]["ARGMAX_typicality"]
    out["directions_differ_by"] = abs(a["chain_mean"] - b["chain_mean"])
    out["SELFTEST_can_fail"] = ("if the two directions agreed the sign would be immaterial; "
                                "they differ by %.4f A" % out["directions_differ_by"])

    ST.save_atomic(os.path.join(RESULTS, "s32_R_typicality_sign.json"), out,
                   module_file=__file__)
    print("n = %d%s\n%s\n" % (n, "  (INCOMPLETE)" if n < 126 else "", out["definition"]))
    for k, r in out["arms"].items():
        c = r["vs_prod"]
        print("%-18s %-28s %.4f  eff %+.4f  %.2fxMDE  %dW/%dL  %d/5 folds  %s"
              % (k, r["gloss"], r["chain_mean"], c["effect"], abs(c["effect_over_mde"]),
                 c["n_better"], c["n_worse"], c["folds_same_sign"], c["verdict"][:34]))
    print("%-18s %-28s %.4f  eff %+.4f"
          % ("BRANCH_MEAN", "zero-skill reference", bmean.mean(), (bmean - prod).mean()))
    print("%-18s %-28s %.4f" % ("PRODUCTION", "", prod.mean()))
    print("\nthe two directions differ by %.4f A" % out["directions_differ_by"])
    print("wrote", os.path.join(RESULTS, "s32_R_typicality_sign.json"))


if __name__ == "__main__":
    main()
