"""s21/d_null.py -- D3: WHICH MATCHED-DISPLACEMENT NULL?  The "72%" depends entirely on it.

    python -m s21.d_null

WHAT IS UNDER ATTACK.  `s21/BRIEF.md` and my own commission both quote Sprint 20 L12:
*"72% of AMBER's damage is move size, not direction"*, and instruct me to hold the lambda
continuation to "the matched-displacement null, the one that matters".  There are TWO
matched-displacement nulls in `s20/results/c_land_null.json`, they are matched to the same
per-coordinate RMS torus magnitude from the same starts, and **they disagree about the sign of
AMBER's direction.**

    d_real(amber)           +0.6198   the minimiser's move
    d_rand_iso(amber)       +0.9977   an ISOTROPIC move of the same torus magnitude
    d_toward_member(amber)  +0.4438   a move of the same magnitude along the geodesic toward
                                      another randomly chosen pool member

So against the isotropic null AMBER's direction is BETTER than noise; against the
toward-member null it is worse.  The "72%" is 0.4438 / 0.6198 -- a ratio against the second.

THE LABEL IS THE DEFECT, NOT THE NUMBER.  `s20/c_land_null.py:24`, `s20/agentC_FINDINGS.md` §4
and `s20/LEDGER.md` L12 each describe `toward_member` as *"realisable, native-free, **zero
information**"*.  It is realisable and native-free.  It is not zero-information: a move of fixed
size toward a uniformly random member of the pool is, in expectation, a move toward the POOL'S
CENTROID, and pool consensus is the one native-free signal this programme has established has
positive in-band skill (S12's score-filter + consensus medoid, -0.172 [-0.316, -0.027]; project
memory `consensus-is-the-only-in-band-discriminator` and `consensus-is-outlier-avoidance`).

BRIEF section 7 rule 4 asks for a **plausible-but-uninformative** control rather than a uniform
one, and `toward_member` is exactly that -- the right control, wearing the wrong label.  But
"uninformative" and "zero information" are different bars, and the difference is what "72%"
means.  This module measures the difference instead of asserting it.

    H0:  toward_member and rand_iso are the same null.
    If they are not, then "72% of the damage is move size" is a statement about a control that
    itself carries skill, and the honest form of the sentence names its null.

No new structure is built and no energy is evaluated: every number is recomputed from Workstream
C's own Sprint-20 artefact, target-paired, with a FOLD-CLUSTERED interval quoted beside the
i.i.d. one (`s12.instrument.paired` is i.i.d. over targets across 293 call sites; `PL.paired`
gives both).

CONSEQUENCE, registered before the numbers below are read: if the two nulls differ, then any
Sprint-21 lambda-continuation arm must be reported against **both**, and "beats toward_member"
is the bar while "beats rand_iso" is not.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s18 import phys_lib as PL      # noqa: E402

SRC = os.path.join(ROOT, "s20", "results", "c_land_null.json")


def run():
    o = json.load(open(SRC))
    rows = o["rows"]
    folds = np.array([int(r["fold"]) for r in rows])
    n = len(rows)
    print(f"\n=== D3  WHICH MATCHED-DISPLACEMENT NULL?  n = {n} targets "
          f"(Workstream C's Sprint-20 subset), 30 draws each ===")
    print("Every arm is matched to the SAME per-coordinate RMS torus magnitude from the SAME")
    print("starts.  Basis: BUILT CHAIN (ideal-geometry CA trace).  CI = fold-clustered, with")
    print("the i.i.d.-over-targets interval beside it.\n")
    out = {"n": n, "src": "s20/results/c_land_null.json", "cells": {}}
    for pot in ("legacy", "amber"):
        mag = np.array([r["null"][pot]["theta_moved"] for r in rows])
        start = np.array([r["null"][pot]["rmsd_start"] for r in rows])
        real = np.array([r["null"][pot]["rmsd_end"] for r in rows])
        iso = np.array([r["null"][pot]["rand_iso"] for r in rows])
        twd = np.array([r["null"][pot]["toward_member"] for r in rows])
        print(f"  [{pot}]  matched magnitude {mag.mean():.3f} rad/coordinate, "
              f"start {start.mean():.3f} A")
        for tag, a, b in (("real   - start", real, start),
                          ("iso    - start", iso, start),
                          ("toward - start", twd, start),
                          ("real   - iso", real, iso),
                          ("real   - toward", real, twd),
                          ("THE NULLS DIFFER: toward - iso", twd, iso)):
            st = PL.paired(a, b, folds=folds)
            ci = st.get("ci_fold", st["ci"])
            out["cells"][f"{pot}|{tag.strip()}"] = {
                "mean": st["mean"], "ci95_fold": list(ci), "ci95_iid": list(st["ci"]),
                "median": st["median"], "W": st["W"], "L": st["L"]}
            print(f"    {tag:<32}{st['mean']:+7.4f} "
                  f"[{ci[0]:+.4f},{ci[1]:+.4f}] fold  "
                  f"[{st['ci'][0]:+.4f},{st['ci'][1]:+.4f}] iid   "
                  f"W/L {st['W']}/{st['L']}  med {st['median']:+.4f}")
        frac = (twd - start).mean() / max(1e-12, (real - start).mean())
        frac_iso = (iso - start).mean() / max(1e-12, (real - start).mean())
        out["cells"][f"{pot}|share_toward"] = float(frac)
        out["cells"][f"{pot}|share_iso"] = float(frac_iso)
        print(f"    share of the move's cost reproduced by the null: "
              f"toward_member {100*frac:.0f}%   isotropic {100*frac_iso:.0f}%")
        print()
    print("READ.  'X% of the damage is move size' is not a property of the potential; it is a")
    print("ratio against a chosen null, and the two matched-magnitude nulls here differ.  The")
    print("honest form of the sentence names its null.  Any lambda-continuation arm this sprint")
    print("must be reported against BOTH, and beating `toward_member` is the bar.")
    json.dump(out, open(os.path.join(RESULTS, "d_null.json"), "w"), indent=1)
    ok = (n >= 30 and all(f"{p}|{t}" in out["cells"] for p in ("legacy", "amber")
                          for t in ("real   - start", "iso    - start", "toward - start",
                                    "real   - iso", "real   - toward",
                                    "THE NULLS DIFFER: toward - iso")))
    p = os.path.join(RESULTS, "d_null.COMPLETE")
    if ok:
        with open(p, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"n={n} potentials=legacy,amber nulls=rand_iso,toward_member "
                     f"comparisons=6 per potential, fold-clustered and iid CIs both recorded, "
                     f"source={out['src']}\n")
        print("\nCOMPLETE.")
    elif os.path.exists(p):
        os.remove(p)
    return out


if __name__ == "__main__":
    run()
