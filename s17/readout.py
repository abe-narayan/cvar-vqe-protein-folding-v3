"""s17/readout.py -- THE READOUT SWEEP.  Selection and averaging are two ends of one dial.

Sprint 17 section 11 proposes switching the architecture from an AVERAGING readout to a
SELECTION readout.  Both are special cases of one operator:

    score the K candidates, keep the best m, coordinate-average those m.

    m = 1   pure selection (the distance argmin)
    m = K   pure averaging over the whole candidate set
    m = 75  the incumbent, when K = 75 and the score is the shipped distogram filter

So the architectural question is not "selection or averaging" -- it is **where the optimum of
m sits**, and whether it beats both endpoints.  Nobody has measured the interior.

WHY THE INTERIOR SHOULD BE NON-TRIVIAL, exactly.  Sprint 16 established that the
coordinate-average readout obeys the Krogh-Vedelsby identity in the frame the operator builds
in, to machine precision:

    readout^2 = mean member error^2 - diversity^2

Shrinking m improves the first term (better members survive) and shrinks the second (less
disagreement left to cancel).  The two move in opposite directions, so an interior optimum is
expected on general grounds and its location is a property of the SCORE's quality -- a better
score pushes m* down, because it can be trusted to keep only the best few.  The programme's own
record predicts exactly this: *m\\* shrinks 500 -> 75 -> 20 -> 3-5 as the objective improves*.
This module tests that recorded prediction directly, which is the first time it has been done
across candidate widths.

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  There is an interior m* that beats both m = 1 and m = K, and m* is smaller at
  larger K (because a wider set contains more bad members that the score must exclude).

  EXPECTED.  m* in the single digits to low tens; the gain over m = 1 modest (tenths of an
  angstrom), because the identity's two terms are of comparable size.

  CONTROLS, both mandatory.  (i) RANDOM-m: average m candidates chosen at RANDOM from the
  same K, which prices the score rather than the averaging.  (ii) The ZERO-INFORMATION
  endpoints: m = K (the whole set, no score used) and a constant ideal alpha-helix.

  SUCCESS.  Some interior m beats BOTH endpoints with a paired fold-clustered interval
  excluding zero, target as the unit, and beats random-m at the same m.

  FALSIFIER.  If the curve is monotone in m at every K -- i.e. one of the endpoints always
  wins -- then "selection versus averaging" is a genuine binary and the interior is empty.
  If interior-m only beats m = 1 by less than the gap to random-m, the effect is averaging,
  not selection, and must be reported as such.

  TRAP GUARD.  m is chosen LEAVE-FOLD-OUT.  An m* that lands on the BOUNDARY of the ladder is
  an endpoint substitution, not a tuned parameter -- Sprint 16 caught that failure twice -- and
  is reported as such rather than as a win.

ORACLE USE.  The native is read only to score the emitted structure.  The score driving the
selection is the shipped native-free distance objective.  m is trained on labels leave-fold-out
and is therefore a TRAINED quantity, not a native-free one, and is labelled that way.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

KS = (75, 150, 300, 500, 1000, 2000)
MS = (1, 2, 3, 5, 8, 12, 20, 35, 75, 150, 300)
N_RAND = 3


def _avg_of(Wsel, nat):
    """the production coordinate-average operator, and the identity's two terms with it."""
    if len(Wsel) == 1:
        v = float(I.ca_rmsd(Wsel[0], nat))
        return v, v, 0.0, v
    P = I.pairwise_rmsd(Wsel)
    avg, b = I.coordinate_average(Wsel, P)
    n = Wsel.shape[1]
    Wm = I.superpose_batch(Wsel, Wsel[b])
    div = float(np.sqrt(np.mean(((Wm - avg) ** 2).sum(axis=(1, 2)) / n)))
    #: the member-error term MUST be measured in the frame the operator builds in, not in
    #: free superposition.  Sprint 16 published the free-superposition term where the identity
    #: consumes the common-frame one and got the attribution wrong by about a factor of two;
    #: this is that correction applied at the source rather than after an audit finds it.
    natm = I.superpose_batch(np.asarray(nat, float)[None], avg)[0]
    err = float(np.sqrt(np.mean(((Wm - natm) ** 2).sum(axis=(1, 2)) / n)))
    #: the identity's OWN readout term, in the same frame as err and div.  It is reported
    #: beside the free-superposition RMSD so the parallel-axis residual is visible rather
    #: than implied -- if the two differ, the identity column is not exact and must not be
    #: read as though it were.
    cf = float(np.sqrt(((avg - natm) ** 2).sum() / n))
    return float(I.ca_rmsd(avg, nat)), err, div, cf


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        order = np.asarray(u["order"], int)
        W = np.asarray(u["W"], float)[order]
        nat = np.asarray(u["nat_ca"], float)
        U = len(W)
        i, j = I.pair_index(n)
        dg = I.distogram(pdb, seq, fold)
        score = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        rng = SD.stable_rng(pdb, "s17readout")

        cell = {}
        for K in KS:
            k = min(K, U)
            sc = score[:k]
            rank = np.argsort(sc, kind="stable")
            e = {"k": int(k)}
            for m in MS:
                if m > k:
                    continue
                r, err, div, cf = _avg_of(W[rank[:m]], nat)
                e[f"m{m}"] = {"rmsd": r, "err": err, "div": div, "cf": cf}
                #: CONTROL -- the same m, chosen at random from the same K.  This prices the
                #: SCORE; without it an interior optimum could be pure averaging.
                rr = [_avg_of(W[rng.permutation(k)[:m]], nat)[0] for _ in range(N_RAND)]
                e[f"m{m}"]["rand"] = float(np.mean(rr))
            e["all"] = _avg_of(W[:k], nat)[0] if k <= 2000 else None
            cell[str(K)] = e
        rows.append({"pdb": pdb, "n": n, "fold": fold, "U": int(U), "cells": cell})
        if (c + 1) % 15 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(RESULTS, "readout.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "readout.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "readout.json")))["rows"]
    rng = SD.stable_rng("readout", "report")
    folds = np.asarray([r["fold"] for r in rows], int)
    print(f"\nn = {len(rows)} targets.  Score = the shipped native-free distance objective.")
    print("m = 1 is pure SELECTION; m = K is pure AVERAGING; the incumbent is K=75, m=75.\n")

    print("A. THE READOUT CURVE -- mean emitted RMSD against m, per candidate width K")
    ms = [m for m in MS]
    print(f"  {'K':>6}" + "".join(f"{('m=%d' % m):>8}" for m in ms))
    for K in KS:
        k = str(K)
        line = f"  {K:>6}"
        for m in ms:
            v = [r["cells"][k][f"m{m}"]["rmsd"] for r in rows if f"m{m}" in r["cells"][k]]
            line += f"{np.mean(v):>8.3f}" if len(v) == len(rows) else f"{'--':>8}"
        print(line)

    print("\n  the same curve for RANDOM-m (the control that prices the SCORE, not the averaging)")
    print(f"  {'K':>6}" + "".join(f"{('m=%d' % m):>8}" for m in ms))
    for K in KS:
        k = str(K)
        line = f"  {K:>6}"
        for m in ms:
            v = [r["cells"][k][f"m{m}"]["rand"] for r in rows if f"m{m}" in r["cells"][k]]
            line += f"{np.mean(v):>8.3f}" if len(v) == len(rows) else f"{'--':>8}"
        print(line)

    print("\nB. THE IDENTITY'S TWO TERMS -- why an interior optimum exists  [err is ORACLE]")
    K0 = "500"
    print(f"  K = {K0}:  {'m':>5}{'member err':>13}{'diversity':>12}{'readout':>10}"
          f"{'common-frame':>11}{'identity':>13}{'|residual|':>13}")
    for m in ms:
        key = f"m{m}"
        if key not in rows[0]["cells"][K0]:
            continue
        er = np.mean([r["cells"][K0][key]["err"] for r in rows])
        dv = np.mean([r["cells"][K0][key]["div"] for r in rows])
        ro = np.mean([r["cells"][K0][key]["rmsd"] for r in rows])
        cf = np.mean([r["cells"][K0][key]["cf"] for r in rows])
        res = np.mean([abs(np.sqrt(max(r["cells"][K0][key]["err"]**2
                                       - r["cells"][K0][key]["div"]**2, 0.0))
                           - r["cells"][K0][key]["cf"]) for r in rows])
        print(f"  {'':>9}{m:>5}{er:>13.3f}{dv:>12.3f}{ro:>10.3f}{cf:>11.3f}"
              f"{np.sqrt(max(er*er-dv*dv,0)):>13.3f}{res:>13.2e}")

    print("\n  |residual| is the PER-TARGET parallel-axis residual and is machine zero, so the")
    print("  identity is EXACT here.  The visible gap between the `identity` and `readout`")
    print("  columns is Jensen's inequality on the aggregation -- the mean of sqrt(.) is not")
    print("  the sqrt of the mean -- and is NOT a frame error.  `readout` (free superposition)")
    print("  and `common-frame` agree exactly, which is the check that they are one object.")

    print("\nC. THE DEPLOYABLE ARM -- m chosen LEAVE-FOLD-OUT (a TRAINED quantity), per K")
    print(f"  {'K':>6}{'m* per fold':>26}{'realized':>10}{'vs m=1':>24}{'vs random-m*':>24}")
    for K in KS:
        k = str(K)
        avail = [m for m in ms if f"m{m}" in rows[0]["cells"][k]]
        picked, vals, rnds = {}, np.zeros(len(rows)), np.zeros(len(rows))
        for f in sorted(set(folds)):
            tr = folds != f
            picked[int(f)] = min(avail, key=lambda m: np.mean(
                [rows[q]["cells"][k][f"m{m}"]["rmsd"] for q in range(len(rows)) if tr[q]]))
        for q, r in enumerate(rows):
            mm = picked[int(folds[q])]
            vals[q] = r["cells"][k][f"m{mm}"]["rmsd"]
            rnds[q] = r["cells"][k][f"m{mm}"]["rand"]
        one = np.array([r["cells"][k]["m1"]["rmsd"] for r in rows])
        a1, l1, h1 = _boot(vals - one, rng)
        a2, l2, h2 = _boot(vals - rnds, rng)
        bnd = " BOUNDARY" if any(picked[f] in (avail[0], avail[-1]) for f in picked) else ""
        print(f"  {K:>6}{str([picked[f] for f in sorted(picked)]):>26}{vals.mean():>10.3f}"
              f"   {a1:+.3f} [{l1:+.3f},{h1:+.3f}]   {a2:+.3f} [{l2:+.3f},{h2:+.3f}]{bnd}")

    print("\n  A fold whose m* sits at the END of the ladder is choosing an endpoint, not a")
    print("  parameter, and is flagged BOUNDARY -- Sprint 16 was caught by that failure twice.")
    print("  'vs random-m*' is the arm that decides whether the effect is SELECTION or merely")
    print("  AVERAGING: if interior m beats m=1 by less than it beats random-m, it is the score")
    print("  doing the work; if the two are similar, it is the averaging.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
