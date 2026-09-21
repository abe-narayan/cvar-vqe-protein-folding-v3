#!/usr/bin/env python
"""S32 LANE R -- R1: does `price = sqrt(e^2 + d^2) - e` predict the S29 O-ladder?

Reads `s32/results/s32_R_laddernull_shard*.jsonl` and emits
`s32/results/s32_R_ladder_null.json`.

WHAT IS BEING TESTED, and what would falsify it
For each rung: `e` = RMSD(cloud, native), `d` = RMSD(projected chain, cloud), and the
observed price is RMSD(chain, native) - e.  The ISOTROPIC NULL is the price you get if the
projection displacement is statistically orthogonal to the cloud's native error.  The claim
is that the object-dependence contract rule 16 records -- prices from -0.0052 to +0.1701 --
is ENTIRELY through `d`.

Falsified per rung if the paired difference (observed price - null price) clears its own MDE:
NEGATIVE means the projection moves TOWARD the native there and some of the price is
recoverable; POSITIVE means it moves away and the price is worse than orthogonal.

REPRODUCTION CHECK: each rung's chain RMSD is compared to the value S29 recorded for the
same (pdb, item).  `prod` must reproduce the canonical 3.210533994943299.
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
CONTRACT_PRICE = {"prod": +0.1622, "bestm": +0.1701, "best1_pool": -0.0030,
                  "sparse_pool_s10": +0.0002, "sparse_pool_s20": -0.0052}


def main():
    rows = {}
    for p in sorted(glob.glob(os.path.join(RESULTS, "s32_R_laddernull_shard*.jsonl"))):
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    rows[r["pdb"]] = r
    pdbs = [t["pdb"] for t in I.targets() if t["pdb"] in rows]
    n = len(pdbs)
    folds = np.array([rows[p]["fold"] for p in pdbs], int)
    out = {"n": n, "prereg_commit": "02754f5a",
           "basis": "built chain CA-RMSD, tuning126; clouds from s29_O_structs",
           "INCOMPLETE": n < 126, "rungs": {}}

    for key, cprice in CONTRACT_PRICE.items():
        have = [p for p in pdbs if key in rows[p]["rungs"]]
        if len(have) < 2:
            continue
        f = np.array([rows[p]["fold"] for p in have], int)
        g = lambda k: np.array([rows[p]["rungs"][key][k] for p in have], float)  # noqa: E731
        e, d = g("cloud_rmsd"), g("d_to_cloud")
        chain = g("chain_rmsd")
        obs = chain - e
        nullp = np.hypot(e, d) - e
        ref = np.array([rows[p]["rungs"][key].get("s29_chain_rmsd") or np.nan for p in have])
        ok = np.isfinite(ref)
        rec = {
            "name": rows[have[0]]["rungs"][key]["name"],
            "ORACLE_NOT_DEPLOYABLE": bool(rows[have[0]]["rungs"][key]["ORACLE_NOT_DEPLOYABLE"]),
            "n": len(have),
            "cloud_mean": float(e.mean()), "chain_mean": float(chain.mean()),
            "chain_se": float(chain.std(ddof=1) / math.sqrt(len(have))),
            "d_mean": float(d.mean()), "d_median": float(np.median(d)),
            "d_p90": float(np.percentile(d, 90)), "d_max": float(d.max()),
            "cloud_vbond_mean": float(g("cloud_vbond_mean").mean()),
            "chain_vbond_mean": float(g("chain_vbond_mean").mean()),
            "price_observed_mean": float(obs.mean()),
            "price_observed_se": float(obs.std(ddof=1) / math.sqrt(len(have))),
            "price_isotropic_null_mean": float(nullp.mean()),
            "price_contract_rule16": cprice,
            "abs_err_vs_contract": abs(float(obs.mean()) - cprice),
            "s29_reproduction": {
                "n_compared": int(ok.sum()),
                "mean_abs_diff": float(np.abs(chain[ok] - ref[ok]).mean()) if ok.any() else None,
                "max_abs_diff": float(np.abs(chain[ok] - ref[ok]).max()) if ok.any() else None,
                "n_bit_identical": int((chain[ok] == ref[ok]).sum()) if ok.any() else 0,
                "s29_chain_mean": float(ref[ok].mean()) if ok.any() else None},
            "observed_minus_null": ST.compare(obs, nullp, f,
                                              names=have,
                                              label="%s: observed price - isotropic null" % key),
        }
        out["rungs"][key] = rec

    ST.save_atomic(os.path.join(RESULTS, "s32_R_ladder_null.json"), out, module_file=__file__)

    print("n = %d targets%s\n" % (n, "  (INCOMPLETE)" if n < 126 else ""))
    hdr = ("%-18s %8s %8s %8s | %9s %9s %9s | %8s %6s"
           % ("rung", "cloud", "chain", "d", "obs price", "iso null", "rule16", "xMDE", "verdict"))
    print(hdr); print("-" * len(hdr))
    for k, v in out["rungs"].items():
        c = v["observed_minus_null"]
        print("%-18s %8.4f %8.4f %8.4f | %+9.4f %+9.4f %+9.4f | %8.2f %s"
              % (k, v["cloud_mean"], v["chain_mean"], v["d_mean"],
                 v["price_observed_mean"], v["price_isotropic_null_mean"],
                 v["price_contract_rule16"], abs(c["effect_over_mde"]),
                 "BELOW null" if c["effect"] < 0 else "above null"))
    print("\nreproduction of S29's own chain RMSD per rung:")
    for k, v in out["rungs"].items():
        s = v["s29_reproduction"]
        print("  %-18s n=%3d  mean|d| %s  max|d| %s  bit-identical %d"
              % (k, s["n_compared"],
                 ("%.3e" % s["mean_abs_diff"]) if s["mean_abs_diff"] is not None else "-",
                 ("%.3e" % s["max_abs_diff"]) if s["max_abs_diff"] is not None else "-",
                 s["n_bit_identical"]))
    print("\nwrote", os.path.join(RESULTS, "s32_R_ladder_null.json"))


if __name__ == "__main__":
    main()
