#!/usr/bin/env python
"""S32 LANE R -- R4 analysis: the repair arms, and the ORACLE scale grid priced as the
order statistic it is.

Reads `s32/results/s32_R_repair_shard*.jsonl` -> `s32/results/s32_R_repair.json`.

BASIS: built-chain CA-RMSD, tuning126, n = 126.  PROD and every arm are produced by
`core.project.lam_path` in the SAME process per target, so every contrast is paired in-job.

CONTRACT RULE 9.  `SCALE_GRID` takes a per-target minimum over 8 dilation factors.  That is
a best-of-K arm and it is priced as one: `stats_lib.best_of_k_within` for the across-target
null and, the number to quote, `split_half_transfer` -- choose the factor on a random half
of the targets, score it on the other half, both directions, repeated.  If the winning
factor is target-specific noise the transfer is ~0 and no separate null is needed.

CONTRACT RULE 15.  Every arm carries its emitted chain's virtual-bond mean and sd and its
displacement from the production chain, from row one.
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
ARMS = ["MEDOID_EXTRA", "SCALE_NF", "SCALE_NF_MED"]


def main():
    R = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "s32_R_repair_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    R[r["pdb"]] = r
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in R]
    n = len(pdbs)
    folds = np.array([R[p]["fold"] for p in pdbs], int)
    prod = np.array([R[p]["PROD"]["rmsd"] for p in pdbs])
    out = {"n": n, "INCOMPLETE": n < 126, "prereg_commit": "02754f5a",
           "basis": "built chain CA-RMSD, tuning126; every arm via lam_path in the same process",
           "PROD_chain_mean": float(prod.mean()),
           "PROD_chain_se": float(prod.std(ddof=1) / math.sqrt(n)),
           "s_nf_mean": float(np.mean([R[p]["s_nf"] for p in pdbs])),
           "arms": {}}

    print("n = %d%s\nPROD chain mean %.6f\n" % (n, "  (INCOMPLETE)" if n < 126 else "",
                                                prod.mean()))
    for k in ARMS:
        v = np.array([R[p][k]["rmsd"] for p in pdbs])
        rec = {
            "chain_mean": float(v.mean()),
            "chain_se": float(v.std(ddof=1) / math.sqrt(n)),
            "chain_median": float(np.median(v)),
            "d_to_cloud_mean": float(np.mean([R[p][k]["d_to_C"] for p in pdbs])),
            "vbond_mean": float(np.mean([R[p][k]["vbond_mean"] for p in pdbs])),
            "vbond_sd": float(np.mean([R[p][k]["vbond_sd"] for p in pdbs])),
            "displacement_from_PROD_mean": float(np.mean([R[p][k]["d_to_prod"] for p in pdbs])),
            "displacement_from_PROD_max": float(np.max([R[p][k]["d_to_prod"] for p in pdbs])),
            "n_identical_to_PROD": int(sum(R[p][k]["d_to_prod"] < 1e-9 for p in pdbs)),
            "ORACLE_NOT_DEPLOYABLE": False,
            "vs_prod": ST.compare(v, prod, folds, names=pdbs, label="%s - PROD" % k),
        }
        out["arms"][k] = rec
        c = rec["vs_prod"]
        print("%-14s %.4f  eff %+.4f  %.2fxMDE  %dW/%dL/%dT  foldCI[%+.4f,%+.4f] %d/5"
              % (k, v.mean(), c["effect"], abs(c["effect_over_mde"]), c["n_better"],
                 c["n_worse"], c["n_tied"], c["ci95_fold"][0], c["ci95_fold"][1],
                 c["folds_same_sign"]))
        print("               %s" % c["verdict"][:88])
        print("               vbond %.4f +- %.4f | moved from PROD: mean %.4f max %.4f | "
              "identical to PROD on %d/%d"
              % (rec["vbond_mean"], rec["vbond_sd"], rec["displacement_from_PROD_mean"],
                 rec["displacement_from_PROD_max"], rec["n_identical_to_PROD"], n))

    # ------------------------------------------------- the ORACLE grid, priced as best-of-K
    keys = sorted(R[pdbs[0]]["SCALE_GRID"].keys(), key=float)
    M = np.array([[R[p]["SCALE_GRID"][s]["rmsd"] for s in keys] for p in pdbs])
    orc = M.min(1)
    bok = ST.best_of_k_within(M, n_boot=400, seed_parts=("s32R", "scalegrid"))
    glob_s = keys[int(np.argmin(M.mean(0)))]
    gl = M[:, keys.index(glob_s)]
    out["SCALE_GRID"] = {
        "ORACLE_NOT_DEPLOYABLE": True, "grid": keys,
        "oracle_per_target_mean": float(orc.mean()),
        "vs_prod_ORACLE": ST.compare(orc, prod, folds, names=pdbs,
                                     label="ORACLE per-target scale - PROD"),
        "best_global_factor": glob_s, "global_mean": float(gl.mean()),
        "vs_prod_global_ORACLE_CHOICE": ST.compare(gl, prod, folds, names=pdbs,
                                                   label="global best scale - PROD"),
        "best_of_k": bok,
    }
    c = out["SCALE_GRID"]["vs_prod_ORACLE"]
    print("\nSCALE_GRID, ORACLE / NOT DEPLOYABLE, %d factors %s" % (len(keys), keys))
    print("  per-target ORACLE  %.4f  eff %+.4f  %.2fxMDE"
          % (orc.mean(), c["effect"], abs(c["effect_over_mde"])))
    cg = out["SCALE_GRID"]["vs_prod_global_ORACLE_CHOICE"]
    print("  best GLOBAL factor %s -> %.4f  eff %+.4f  %.2fxMDE  (the factor itself chosen "
          "on native RMSD: ORACLE)" % (glob_s, gl.mean(), cg["effect"],
                                       abs(cg["effect_over_mde"])))
    print("  BEST-OF-K PRICING (contract rule 9): observed per-target gain %+.4f"
          % bok["observed_gain"])
    print("    across-target null      %+.4f  (%.0f%% of the observed gain accounted)"
          % (bok["null_across_targets"], 100 * bok["share_accounted"]))
    print("    SPLIT-HALF TRANSFER     %+.4f  (%.0f%% of oracle)  <- the number to quote"
          % (bok["split_half"], 100 * bok["split_half_frac"]))
    print("    k_eff %.2f of %d   |   %s" % (bok["k_eff"], len(keys), bok["verdict"]))

    ST.save_atomic(os.path.join(RESULTS, "s32_R_repair.json"), out, module_file=__file__)
    print("\nwrote", os.path.join(RESULTS, "s32_R_repair.json"))


if __name__ == "__main__":
    main()
