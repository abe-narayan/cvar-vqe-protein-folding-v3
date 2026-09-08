"""s17/sel_hard.py -- THE HARD TARGETS, AGAINST THE FULL-UNIVERSE CEILING.  Sprint 17 section 38, PREREG E5.

The programme's "53 of 126 targets have a pool best above 2.0 A" is a statement about the
FIRST 500 entries of a BLOSUM ranking, not about the library.  This module recomputes the hard
set against the full window universe (13,000-27,000 windows per target) and then asks what,
if anything, the residual hard targets have in common.

Reads `s17/results/sel_bench.json` (ORACLE ceilings at K = 75 / 500 / 2000 / full, all on the
identical candidate sets) and `s17/results/sel_cal.json` (native-free per-target features).
Nothing here is a predictor; it is a characterisation, and every class statement carries the
class's own base rate so a "3 of 4 hard targets are X" cannot be read as evidence when X is
also 70% of the easy set.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I         # noqa: E402
from s17 import sel_lib as L            # noqa: E402
from s17 import sel_obj as O            # noqa: E402


def report():
    B = {r["pdb"]: r for r in json.load(open(os.path.join(L.RESULTS, "sel_bench.json")))["rows"]}
    C = {r["pdb"]: r for r in json.load(open(os.path.join(L.RESULTS, "sel_cal.json")))["rows"]}
    pdbs = [p for p in B if p in C]
    n = len(pdbs)
    o75 = np.array([B[p]["cells"]["75"]["oracle"] for p in pdbs])
    o500 = np.array([B[p]["cells"]["500"]["oracle"] for p in pdbs])
    ofull = np.array([B[p]["cells"]["full"]["oracle"] for p in pdbs])
    d500 = np.array([B[p]["cells"]["500"]["dist"] for p in pdbs])
    dfull = np.array([B[p]["cells"]["full"]["dist"] for p in pdbs])
    nn = np.array([B[p]["n"] for p in pdbs])
    fold = np.array([B[p]["fold"] for p in pdbs])
    U = np.array([B[p]["U"] for p in pdbs])

    print(f"\n{'='*96}\nTHE HARD TARGETS AGAINST THE FULL-UNIVERSE CEILING   n = {n}\n{'='*96}")
    print("All ORACLE quantities.  'hard' = the ORACLE best in the candidate set exceeds 2.0 A,")
    print("i.e. no window in the set is within 2.0 A of the native.\n")
    for th in (2.0, 2.5):
        print(f"  above {th} A:  top-75 {int((o75>th).sum()):3d}/{n}   "
              f"K=500 {int((o500>th).sum()):3d}/{n}   FULL {int((ofull>th).sum()):3d}/{n}")
    hard500 = set(np.array(pdbs)[o500 > 2.0])
    hardfull = set(np.array(pdbs)[ofull > 2.0])
    print(f"\n  the K=500 hard set has {len(hard500)} targets; {len(hardfull)} survive the full universe")
    print(f"  RESCUED by widening (hard at 500, easy at full): {len(hard500 - hardfull)}")
    print(f"  hard at full but not at 500: {len(hardfull - hard500)}  (must be 0 by nesting)")
    print(f"  FAIL18 (zero-recall) still above 2.0 A over the full universe: "
          f"{len([p for p in I.FAIL18 if p in hardfull])}/18")

    h = np.array([p in hardfull for p in pdbs])
    print(f"\nA. WHAT THE RESIDUAL HARD SET LOOKS LIKE   ({h.sum()} hard / {(~h).sum()} easy)")
    print(f"  {'quantity':<20}{'hard':>10}{'easy':>10}{'CI on the difference':>28}")
    fk = sorted(C[pdbs[0]]["feat"].keys())
    X = np.array([[C[p]["feat"][k] for k in fk] for p in pdbs], float)
    rowsx = [("chain length", nn.astype(float)), ("universe size", U.astype(float)),
             ("ORACLE best, full", ofull), ("ORACLE best, K=500", o500),
             ("realized dist, full", dfull)]
    rowsx += [(k, X[:, q]) for q, k in enumerate(fk)]
    for name, v in rowsx:
        a, b = v[h], v[~h]
        rng = np.random.default_rng(0)
        bs = np.array([a[rng.integers(0, len(a), len(a))].mean()
                       - b[rng.integers(0, len(b), len(b))].mean() for _ in range(4000)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        star = "*" if lo > 0 or hi < 0 else " "
        print(f" {star}{name:<20}{a.mean():>10.3f}{b.mean():>10.3f}"
              f"   {a.mean()-b.mean():+.3f} [{lo:+.3f},{hi:+.3f}]")
    print("  * = the bootstrap CI on the hard-easy difference excludes zero.")

    print(f"\nB. BASE RATES BY FOLD AND LENGTH  (a class statement needs its own base rate)")
    for f in np.unique(fold):
        m = fold == f
        print(f"    fold {f}:  hard {int(h[m].sum())}/{int(m.sum())} = {h[m].mean():.1%}")
    for lo, hi in ((9, 11), (12, 13), (14, 16)):
        m = (nn >= lo) & (nn <= hi)
        print(f"    length {lo}-{hi}:  hard {int(h[m].sum())}/{int(m.sum())} = {h[m].mean():.1%}")

    print(f"\nC. IS THE HARD SET A SELECTION PROBLEM OR A RETRIEVAL PROBLEM?")
    print(f"  On the HARD set   ORACLE full {ofull[h].mean():.3f}   realized {dfull[h].mean():.3f}"
          f"   gap {dfull[h].mean()-ofull[h].mean():+.3f}")
    print(f"  On the EASY set   ORACLE full {ofull[~h].mean():.3f}   realized {dfull[~h].mean():.3f}"
          f"   gap {dfull[~h].mean()-ofull[~h].mean():+.3f}")
    print("  If the GAP is similar on both, the hard set is a retrieval problem on top of the")
    print("  same selection problem, not a different selection problem.")

    print(f"\nD. WHAT THE WHOLE-INSTRUMENT MEAN WOULD BE IF THE HARD SET WERE SOLVED")
    print(f"  current ORACLE full-universe ceiling         {ofull.mean():.3f}")
    print(f"  same, with every hard target capped at 2.0   "
          f"{np.minimum(ofull, np.where(h, 2.0, ofull)).mean():.3f}")
    print(f"  contribution of the {int(h.sum())} hard targets to the ceiling: "
          f"{(ofull[h].sum() - 2.0*h.sum())/n:+.3f} A of the mean")


if __name__ == "__main__":
    report()
