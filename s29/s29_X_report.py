#!/usr/bin/env python
"""s29/s29_X_report.py -- lane X: print the probe in LANE T's LADDER ORDER (S29-L15).

The ladder is read in this order and the first failure stops the reading:

    GATE 1  TV(p_Gamma, p_0) on the SAMPLED distribution must exceed 0.45, because below
            that the readout provably cannot resolve the difference (S25 L15).
    GATE 2  the correctly-named classical counterpart -- for a LOCAL mixer that is a
            classical thermal sampler / simulated annealing over the SAME configuration
            space at matched evaluations, NOT an eigensolver.
    GATE 3  the endpoint contrast on the BUILT CHAIN.

Beside them, the lane's own pre-registered falsifiers: D1 (does recombination contain
anything, against the SCRAMBLED matched null lane D required) and P1 to P5.

    python s29/s29_X_report.py [--json s29/results/s29_X_probe.json]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

DEFAULT = os.path.join(HERE, "results", "s29_X_probe.json")

LADDER = [
    ("GATE 2 (T's classical counterpart: SA over the same space at matched evaluations)",
     ["P2 VQE R2 - SA R2 (chain)", "P2 VQE R2 - SA R2 (cloud)"]),
    ("GATE 3 (the endpoint, built chain)",
     ["P1 VQE_g1_s0|R2 - PRODUCTION (chain)", "P1 VQE_g1_s0|R1 - PRODUCTION (chain)",
      "P1 VQE_g0_s0|R2 - PRODUCTION (chain)", "P1 EXACT_argmin - PRODUCTION (chain)",
      "P1 EXACT_top8 - PRODUCTION (chain)", "D2 EXACT top-75 - PRODUCTION (chain)",
      "P1 CVAROPT|R2 - PRODUCTION (chain)", "P1 SA_s0|R2 - PRODUCTION (chain)"]),
    ("D1 (does recombination contain anything: the SCRAMBLED matched null first)",
     ["D1 chimera ORACLE best - SCRAMBLED ORACLE best (cloud)",
      "D1u chimera ORACLE best - best PARENT (cloud, UNMATCHED order statistic)",
      "D1p chimera ORACLE best - POOL ORACLE best K=500 (cloud)"]),
    ("P3 / P3c (the quantum stage beyond m, with lane D's PR-matched control)",
     ["P3 VQE R2 - VQE R1 (chain)", "P3c VQE R2 - PR-matched random weights (chain)"]),
    ("P4 / P5 and the charter's ten controls (built chain)",
     ["P4 VQE R2 - GIBBS(matched) R2 (chain)", "P5 VQE gamma - VQE gamma=0, R2 (chain)",
      "M6 VQE(gamma=0) R2 - CVaR-OPTIMAL LAW R2 (chain)",
      "M6g VQE(gamma) R2 - CVaR-OPTIMAL LAW R2 (chain)",
      "C-untrained VQE R2 - UNTRAINED R2 (chain)",
      "C-product VQE R2 - PRODUCT R2 (chain)",
      "C-seed VQE R2 s0 - s1 (chain)",
      "C-diagonalised VQE R2 - GS R2 (chain)",
      "C-perm VQE R2 - PERM VQE R2 (chain)",
      "C-ordstat VQE R2 - EXACT top-m (chain)"]),
]

ARMS = ["ORACLE_best_chimera", "EXACT_argmin", "EXACT_top8", "EXACT_top75",
        "EXACT_topm_matched", "CVAROPT|R1", "CVAROPT|R2", "GIBBS_T1|R2", "GIBBS_match|R2",
        "SA_s0|R2", "VQE_g0_s0|R1", "VQE_g0_s0|R2", "VQE_g1_s0|R1", "VQE_g1_s0|R2",
        "VQE_g1_s0|R3", "VQE_g1_s1|R2", "VQE_g05_s0|R2", "VQE_g2_s0|R2",
        "VQE_prod_s0|R2", "UNTRAINED_s0|R2", "GS_g1|R2", "GS_g0|R2", "PERM_VQE_g1_s0|R2",
        "PERM_EXACT_top75"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=DEFAULT)
    a = ap.parse_args()
    o = json.load(open(a.json))
    n = o["n"]
    print(f"S29 LANE X PROBE -- n = {n} targets: {', '.join(o['pdbs'])}")
    print(f"artefact: {a.json}\n")

    print("=" * 78)
    print("GATE 1 (lane T, S29-L15): TV(p_Gamma, p_0) on the SAMPLED distribution > 0.45")
    print("=" * 78)
    tg = o.get("tv_gate", {})
    if not tg:
        print("  (absent)")
    for k in ("VQE_g05_s0", "VQE_g1_s0", "VQE_g2_s0"):
        if k in tg:
            v = tg[k]
            print("  %-12s mean TV %.4f  median %.4f  [%.4f, %.4f]  above 0.45 on %d/%d"
                  % (k, v["mean"], v["median"], v["min"], v["max"],
                     v["n_above_gate"], v["n"]))
    passed = any(v["mean"] > 0.45 for v in tg.values()) if tg else False
    print("  GATE 1 VERDICT: %s" % ("PASSES (a grid point's mean TV clears 0.45)" if passed
                                    else "FAILS -- the cell is empty for this instrument"))

    print("\n" + "=" * 78)
    print("MEANS BY ARM (A; cloud = point cloud, chain = the reporting basis)")
    print("=" * 78)
    print("  %-24s %8s %8s %8s %8s" % ("arm", "cloud", "med", "chain", "med"))
    for arm in ARMS:
        c = o["arms"].get(f"{arm}|rmsd_cloud")
        h = o["arms"].get(f"{arm}|rmsd_chain")
        if c is None and h is None:
            continue
        print("  %-24s %8.4f %8.4f %8s %8s"
              % (arm, c["mean"] if c else float("nan"), c["median"] if c else float("nan"),
                 ("%8.4f" % h["mean"]) if h else "--",
                 ("%8.4f" % h["median"]) if h else "--"))
    prod = o.get("production")
    if prod:
        print("  %-24s %8.4f %8s %8.4f" % ("PRODUCTION (DIS s0)", prod["cloud"], "", prod["chain"]))

    print("\n" + "=" * 78)
    print("ORACLE DIAGNOSTICS (labelled ORACLE; they choose nothing)")
    print("=" * 78)
    for k, v in o.get("oracle", {}).items():
        print("  %-36s mean %8.4f  median %8.4f" % (k, v["mean"], v["median"]))

    print("\n" + "=" * 78)
    print("MECHANISM (contract rule 18; addendum 20(c); lane T's M5)")
    print("=" * 78)
    for k, v in o.get("meta", {}).items():
        if isinstance(v, dict):
            print("  %-28s mean %10.4f  [%.4f, %.4f]" % (k, v["mean"], v["min"], v["max"]))
        else:
            print("  %-28s %10.4f" % (k, v))

    for title, labels in LADDER:
        print("\n" + "=" * 78)
        print(title)
        print("=" * 78)
        for lab in labels:
            if lab in o.get("fmt", {}):
                print(o["fmt"][lab])
                print()
    missing = [l for _, ls in LADDER for l in ls if l not in o.get("fmt", {})]
    if missing:
        print("  (absent from this artefact: %s)" % ", ".join(missing))


if __name__ == "__main__":
    main()
