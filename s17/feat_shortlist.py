"""s17/feat_shortlist.py -- F7: is the ESM contact channel a SHORTLIST CONSTRUCTOR?

WHY THIS EXISTS.  F1 found the largest in-band RANK CORRELATION any native-free signal has shown
in this sprint (ESM contact head, rho_S +0.116 [+0.047, +0.184] in the shipped top-75) while its
ARGMIN stayed null.  Those two facts together say the signal orders the band slightly but cannot
pick its best.  L12 and L14 identify an unoccupied role for exactly such a signal: **choosing the
shortlist** rather than choosing within it.  L12 measured that the score-ranked top-B DESTROYS the
shortlist ceiling as K grows (top-75 ORACLE degrades 0.468 A from K = 75 to the full universe),
and PHYSICS independently measured the shipped shortlist's ORACLE best as significantly WORSE than
a matched random shortlist of the same size.

WHAT IS MEASURED, per target at K = 500, for shortlist sizes B = 25, 75, 150:

    ceiling   ORACLE best inside the shortlist          <- the quantity L12 says is destroyed
    mean      shortlist mean                            <- the quantity that improves instead
    argmin    realized pick by the shipped score inside the shortlist
    medoid    realized consensus medoid inside the shortlist
    p_sub2    fraction of targets whose shortlist contains a sub-2.0 A candidate

ARMS.  `dist` (the incumbent top-B), `esmcon` (top-B by the best F1 score), `rand` (MATCHED
RANDOM of the same size, the control without which none of this is readable), and `mix` (B/2 from
each score, deduplicated and topped up by the distance score to EXACTLY B, so the union arm is
never priced against a larger shortlist -- that free win is the trap this arm exists to avoid).

A CEILING GAIN IS NOT A RESULT.  L14's pre-registration made the point and its own smoke test
then demonstrated it: `kmeans_best` improved the shortlist ceiling by 0.275 A and the realized
argmin was IDENTICAL across every arm, because the argmin is score-driven and the score's own
favourite is in every shortlist.  So the realized readouts are reported beside every ceiling here.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s17 import sel_lib as L               # noqa: E402
from s17 import feat_lib as F              # noqa: E402

K = 500
SIZES = (25, 75, 150)
OUTJSON = os.path.join(F.RESULTS, "feat_shortlist.json")


def esmcon_score(t, p):
    """The F1 arm with the largest in-band rho: contact-probability-weighted realised contact."""
    con = F.esm_contacts()[t["seq"]]
    i, j = np.asarray(p["i"], int), np.asarray(p["j"], int)
    pc = np.clip(np.asarray(con, float)[i, j], 1e-6, 1 - 1e-6)
    s = F.soft_contact(np.asarray(p["D"], float))
    return -(s * pc[None, :]).sum(1) / max(pc.sum(), 1e-9)


def one_target(t):
    p = L.pack(t["pdb"], K, want=("D", "W"))
    rr = p["rr"]
    W = p["W"]
    sc = L.shipped(p)
    ec = esmcon_score(t, p)
    rng = SD.stable_rng(t["pdb"], "s17featshort")
    rk_d = np.argsort(sc, kind="stable")
    rk_e = np.argsort(ec, kind="stable")
    row = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
           "pool_best": float(rr.min()), "arms": {}}
    for B in SIZES:
        half = B // 2
        mix = list(dict.fromkeys(list(rk_d[:half]) + list(rk_e[:half])))
        q = 0
        while len(mix) < B and q < len(rk_d):          # top up to EXACTLY B, never more
            if rk_d[q] not in mix:
                mix.append(int(rk_d[q]))
            q += 1
        arms = {"dist": rk_d[:B], "esmcon": rk_e[:B],
                "rand": rng.choice(len(rr), B, replace=False), "mix": np.array(mix[:B], int)}
        cell = {}
        for nm, idx in arms.items():
            idx = np.asarray(idx, int)
            sub = rr[idx]
            #: realized readouts, so a ceiling is never quoted alone (L14)
            P = I.pairwise_rmsd(W[idx])
            cell[nm] = {"ceiling": float(sub.min()), "mean": float(sub.mean()),
                        "argmin": L.sel_of(sc[idx], sub),
                        "medoid": float(sub[int(np.argmin(P.mean(1)))]),
                        "sub2": int(sub.min() < 2.0)}
        row["arms"][str(B)] = cell
    return row


def run(targets=None, verbose=True):
    tg = targets if targets is not None else I.targets()
    rows, t0 = [], time.time()
    for q, t in enumerate(tg):
        rows.append(one_target(t))
        if verbose and (q + 1) % 20 == 0:
            print(f"  {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(OUTJSON, "w"))
    json.dump({"rows": rows}, open(OUTJSON, "w"))
    print(f"wrote {OUTJSON}", flush=True)
    return rows


def report(rows=None):
    rows = rows if rows is not None else json.load(open(OUTJSON))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    print(f"\nF7  SHORTLIST CONSTRUCTION   n = {len(rows)}, K = 500")
    print("   `ceiling` is an ORACLE quantity (the shortlist's best member).  A ceiling gain is")
    print("   worthless unless a realized readout can consume it, so argmin and medoid sit beside it.")
    for B in SIZES:
        b = str(B)
        print(f"\n  --- shortlist size B = {B}")
        print(f"     {'arm':<10}{'CEILING':>9}{'mean':>9}{'argmin':>9}{'medoid':>9}{'P(sub-2A)':>11}")
        for nm in ("dist", "esmcon", "rand", "mix"):
            g = lambda k: np.array([r["arms"][b][nm][k] for r in rows], float)   # noqa: E731
            print(f"     {nm:<10}{g('ceiling').mean():>9.3f}{g('mean').mean():>9.3f}"
                  f"{g('argmin').mean():>9.3f}{g('medoid').mean():>9.3f}{g('sub2').mean():>11.3f}")
        print("     paired contrasts on the CEILING (negative = better):")
        for a, c in (("esmcon", "dist"), ("esmcon", "rand"), ("dist", "rand"), ("mix", "dist")):
            x = np.array([r["arms"][b][a]["ceiling"] for r in rows], float)
            y = np.array([r["arms"][b][c]["ceiling"] for r in rows], float)
            print("       " + L.fmt_pair(L.report_pair(f"{a} vs {c}", x, y, fold)))
        print("     paired contrasts on the REALIZED medoid:")
        for a, c in (("esmcon", "dist"), ("mix", "dist")):
            x = np.array([r["arms"][b][a]["medoid"] for r in rows], float)
            y = np.array([r["arms"][b][c]["medoid"] for r in rows], float)
            print("       " + L.fmt_pair(L.report_pair(f"{a} vs {c}", x, y, fold)))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
        report()
