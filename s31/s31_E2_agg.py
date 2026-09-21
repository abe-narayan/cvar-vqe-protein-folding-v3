#!/usr/bin/env python
"""s31/s31_E2_agg.py -- S31 lane E: aggregate the applied arms on BOTH bases.

Endpoint basis, stated on every row: BUILT CHAIN (production 3.2105 A, the S29 `prod` row; this
lane's own re-projection is printed beside it and every delta is PAIRED against that same
re-projection, so the unpinned-projection spread cancels).  CA POINT CLOUD (production 3.0483 A)
is an intermediate and is reported beside it.

Registered verdict (`s31/PREREG_S31_E.md` §3, unchanged by AMENDMENT 1):
    RESULT      delta < 0 and |delta| >= 1.0x MDE, fold CI excluding zero, >= 4/5 folds same sign
    NOT MEASURED    0.7x <= |delta| < 1.0x MDE
    FALSIFIED   |delta| < 0.7x MDE, or delta > 0

    python s31/s31_E2_agg.py
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s31.s31_E2_deltas import ORACLE_ARMS  # noqa: E402

OUT = os.path.join(HERE, "results", "s31_E2_applied.json")
PROD_CHAIN_CANON, PROD_CLOUD_CANON = 3.2105, 3.0483


ALL_ORACLE = False


def status(nm):
    return ("baseline" if nm == "PROD" else
            "ORACLE / NOT DEPLOYABLE" if (ALL_ORACLE or nm in ORACLE_ARMS) else
            "LFO-supervised, native-free at inference")


def verdict(c):
    x = abs(c["effect_over_mde"])
    if c["effect"] > 0:
        return "HARMFUL/FALSIFIED (positive delta)" if x >= 0.7 else "FALSIFIED (|d| < 0.7x MDE, positive)"
    if x < 0.7:
        return "FALSIFIED (|d| %.4f < 0.7x MDE %.4f)" % (abs(c["effect"]), 0.7 * c["mde"])
    if x < 1.0:
        return "NOT MEASURED (0.7-1.0x MDE)"
    if c["ci95_fold"] is not None and c["ci95_fold"][1] >= 0:
        return "NOT MEASURED (fold CI includes zero)"
    if c["folds_same_sign"] < 4:
        return "NOT MEASURED (%d/5 folds same sign)" % c["folds_same_sign"]
    return "RESULT"


def main():
    global ALL_ORACLE
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows-glob", default="s31_E2_rows*.jsonl")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--all-oracle", action="store_true")
    a = ap.parse_args()
    ALL_ORACLE = a.all_oracle
    outpath = a.out
    R = {}
    for p in sorted(glob.glob(os.path.join(HERE, "results", a.rows_glob))):
        for ln in open(p):
            try:
                r = json.loads(ln)
            except Exception:
                continue
            R[r["pdb"]] = r
    tg = I.targets()
    rows = [R[t["pdb"]] for t in tg if t["pdb"] in R]
    missing = [t["pdb"] for t in tg if t["pdb"] not in R]
    folds = np.array([r["fold"] for r in rows]); pdbs = [r["pdb"] for r in rows]
    f18 = np.array([bool(r["fail18"]) for r in rows])
    names = sorted(k[:-6] for k in rows[0] if k.endswith("_cloud"))

    res = {"n": len(rows), "n_missing": len(missing), "missing": missing,
           "prod_chain_canonical_S29": PROD_CHAIN_CANON,
           "prod_cloud_canonical_S29": PROD_CLOUD_CANON,
           "prod_chain_this_run": float(np.mean([r["PROD_chain"] for r in rows])),
           "prod_cloud_this_run": float(np.mean([r["PROD_cloud"] for r in rows])),
           "note": ("every delta is PAIRED against PROD computed in the same process; the "
                    "unpinned multi-start projection seed therefore cancels in the delta"),
           "arms": {}}
    for basis in ("chain", "cloud"):
        base = np.array([r["PROD_" + basis] for r in rows])
        for nm in names:
            if nm == "PROD":
                continue
            v = np.array([r[nm + "_" + basis] for r in rows])
            c = ST.compare(v, base, folds=folds, names=pdbs, label="%s|%s" % (nm, basis))
            d = {"status": status(nm), "basis": basis, "mean": float(v.mean()),
                 "median": float(np.median(v)), "prod": float(base.mean()),
                 "delta": c["effect"], "median_delta": c["median_effect"], "se": c["se"],
                 "mde": c["mde"], "x_mde": c["effect_over_mde"], "ci95_fold": c["ci95_fold"],
                 "ci95_iid": c["ci95_iid"], "per_fold": c["per_fold"],
                 "folds_same_sign": c["folds_same_sign"], "W": c["n_better"], "L": c["n_worse"],
                 "type_m": c["type_m"], "type_m_flag": c["type_m_flag"],
                 "stats_lib_verdict": c["verdict"], "PREREG_VERDICT": verdict(c),
                 "FAIL18_delta": float((v - base)[f18].mean()),
                 "other108_delta": float((v - base)[~f18].mean())}
            res["arms"].setdefault(nm, {})[basis] = d
    for nm in res["arms"]:
        ch, cl = res["arms"][nm]["chain"]["delta"], res["arms"][nm]["cloud"]["delta"]
        res["arms"][nm]["cloud_to_chain_transfer"] = float(ch / cl) if abs(cl) > 1e-9 else None

    tmp = outpath + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, outpath)

    print("n = %d  (missing %d)   PROD chain %.4f (canonical 3.2105)  cloud %.4f (canonical 3.0483)"
          % (res["n"], res["n_missing"], res["prod_chain_this_run"], res["prod_cloud_this_run"]))
    for basis in ("chain", "cloud"):
        print("\n=== BUILT CHAIN ===" if basis == "chain" else "\n=== CA POINT CLOUD ===")
        print("%-22s %-40s %7s %8s %6s %7s %-38s"
              % ("arm", "status", "mean", "delta", "xMDE", "W/L", "PREREG verdict"))
        order = sorted(res["arms"], key=lambda k: res["arms"][k][basis]["delta"])
        for nm in order:
            a = res["arms"][nm][basis]
            print("%-22s %-40s %7.4f %+8.4f %6.2f %3d/%-3d %-38s"
                  % (nm, a["status"], a["mean"], a["delta"], abs(a["x_mde"]), a["W"], a["L"],
                     a["PREREG_VERDICT"]))
    print("\nwrote", outpath)


if __name__ == "__main__":
    main()
