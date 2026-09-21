#!/usr/bin/env python
"""s31/s31_E7_contrasts.py -- S31 lane E: persist every paired contrast S31-L22 quotes.

Contract rule 17: every meaningful number carries an artefact path.  The cross-file contrasts
(E2 arms against E5 control arms) live in neither aggregate, so they are computed and stored here.
Endpoint basis: BUILT CHAIN, paired against PROD in the same rows; CA cloud reported beside it.
Every arm involved is ORACLE / NOT DEPLOYABLE except the four marked native-free.

    python s31/s31_E7_contrasts.py
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

OUT = os.path.join(HERE, "results", "s31_E7_contrasts.json")

PAIRS = [
    ("ORACLE_Y_ALONG_MU", "ORACLE_FULL", "the along half vs the PERFECT prior -- expect a TIE"),
    ("ORACLE_Y_ALONG_MU", "ORACLE_Y_PERP_MU", "the decomposition asymmetry"),
    ("ORACLE_Y_PERP_MU", "CTRL_SHRINK_ORACLE_Y", "perp vs its NORM-matched control"),
    ("ORACLE_Y_PERP_MU", "ENERGYMATCHED_PERP", "perp vs its ENERGY-MATCHED-DIRECTION control"),
    ("ORACLE_Y_ALONG_MU", "ENERGYMATCHED_ALONG", "along vs its ENERGY-MATCHED-DIRECTION control"),
    ("ORACLE_Y_ALONG_MU", "SHRINK_Y_090", "along vs norm-matched shrinkage -- the honest caveat"),
    ("ENERGYMATCHED_ALONG", "SHRINK_Y_090", "energy-matched along vs the same shrinkage"),
    ("ENERGYMATCHED_PERP", "CTRL_SHRINK_ORACLE_Y", "energy-matched perp vs norm-matched shrinkage"),
    ("ORACLE_Y_PERP_MU", "PROD", "is the orthogonal correction harmful on its own?"),
    ("RANDDIR_PERP", "ORACLE_FULL", "the degenerate isotropic split, for scale"),
    ("E2o_PROJ_MU", "FIT_N3", "ORACLE projection vs not projecting"),
    ("E3o_CONSTR_MU", "FIT_N3", "ORACLE constrained refit vs not constraining"),
    ("E2_PROJ_MUHAT", "FIT_N3", "NATIVE-FREE projection vs not projecting"),
    ("E2_PROJ_MUHAT", "CTRL_SHRINK_E2", "NATIVE-FREE projection vs plain shrinkage"),
    ("E3_CONSTR_MUHAT", "CTRL_SHRINK_E2", "NATIVE-FREE constrained refit vs plain shrinkage"),
]
NATIVE_FREE = {"FIT_N3", "E2_PROJ_MUHAT", "E3_CONSTR_MUHAT", "CTRL_SHRINK_E2", "CTRL_PROJ_RAND",
               "PROD"}


def main():
    R = {}
    for p in glob.glob(os.path.join(HERE, "results", "s31_E2_rows.s*.jsonl")) + \
             glob.glob(os.path.join(HERE, "results", "s31_E5_rows.s*.jsonl")):
        for ln in open(p):
            r = json.loads(ln)
            R.setdefault(r["pdb"], {}).update(r)
    tg = I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.array([int(t["fold"]) for t in tg])
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    assert all(p in R for p in pdbs), "missing targets"

    def arm(nm, basis):
        return np.array([R[p][nm + "_" + basis] for p in pdbs])

    out = {"n": len(pdbs), "basis_note": "primary = built chain; CA cloud beside it",
           "prod_chain_this_run": float(arm("PROD", "chain").mean()),
           "prod_chain_canonical_S29": 3.2105,
           "prod_cloud_this_run": float(arm("PROD", "cloud").mean()),
           "contrasts": {}}
    for a_, b_, why in PAIRS:
        rec = {"why": why,
               "status": ("both arms native-free at inference"
                          if a_ in NATIVE_FREE and b_ in NATIVE_FREE
                          else "ORACLE / NOT DEPLOYABLE (at least one arm reads the native)")}
        for basis in ("chain", "cloud"):
            c = ST.compare(arm(a_, basis), arm(b_, basis), folds=folds, names=pdbs,
                           label="%s vs %s|%s" % (a_, b_, basis))
            d = arm(a_, basis) - arm(b_, basis)
            rec[basis] = {"effect": c["effect"], "median_effect": c["median_effect"],
                          "se": c["se"], "mde": c["mde"], "x_mde": c["effect_over_mde"],
                          "ci95_fold": c["ci95_fold"], "per_fold": c["per_fold"],
                          "folds_same_sign": c["folds_same_sign"],
                          "W": c["n_better"], "L": c["n_worse"],
                          "type_m_flag": c["type_m_flag"],
                          "FAIL18_delta": float(d[f18].mean()),
                          "other108_delta": float(d[~f18].mean())}
        out["contrasts"]["%s vs %s" % (a_, b_)] = rec
    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, OUT)
    print("PROD chain this run %.4f (canonical 3.2105)   n = %d"
          % (out["prod_chain_this_run"], out["n"]))
    print("%-48s %9s %7s %4s %8s  %s" % ("contrast (chain, negative = first better)", "effect",
                                         "xMDE", "ss", "W/L", "status"))
    for k, v in out["contrasts"].items():
        c = v["chain"]
        print("%-48s %+9.4f %7.2f %4d %4d/%-4d %s"
              % (k, c["effect"], abs(c["x_mde"]), c["folds_same_sign"], c["W"], c["L"],
                 "ORACLE" if v["status"].startswith("ORACLE") else "native-free"))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
