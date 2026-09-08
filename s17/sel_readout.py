"""s17/sel_readout.py -- THE READOUT TRIANGLE, on one instrument and identical candidate sets.

THE QUESTION.  Three readouts consume the same shortlist and the programme has never priced
them against each other on one instrument:

    argmin      return the top-1 by the distance objective          (a SELECTION readout)
    medoid      return the member minimising mean RMSD to the rest  (a SELECTION readout)
    average     superpose on the medoid and take the coordinate mean (an AVERAGING readout)

The numbers currently in circulation come from three different places -- argmin 3.454 from the
shipped pipeline, medoid 3.282 from Sprint 8, average ~3.05 from Sprint 16's all-atom stage --
and they are not comparable.  This module puts all three on the same candidate sets, the same
n = 126, the same CA-only RMSD definition, at every (K, m).

WHY IT MATTERS ARCHITECTURALLY.  Sprint 17's central hypothesis is that the readout should move
from AVERAGING to SELECTION, because averaging obeys `readout^2 = mean member error^2 -
diversity^2` and is capped near 3.05 A while selection has a 1.711 A ceiling.  If averaging in
fact BEATS both selection readouts at every m on identical candidates, the hypothesis is not
refuted -- the ceilings are still what they are -- but the claim that switching readouts is an
immediate improvement is, and that distinction has to be measured rather than argued.

CONTROLS.  Every readout is computed on the score's top-m AND on a matched random m drawn from
the same K-set, so the score's contribution to the shortlist is separated from the readout's
contribution to the answer.  Ties are handled by `sel_lib.sel_of` wherever an argmin is taken.

COST NOTE.  The top-m sets are NESTED in m, so the O(m^2) pairwise matrix is computed once per
(target, K) at the largest m and every smaller m is a sub-block.  Recomputing it per m would be
13x the work for identical numbers.
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

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from s17 import sel_lib as L             # noqa: E402

KS = (75, 500, 2000)
MS = (1, 2, 3, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300)
MMAX = max(MS)
NRAND = 4          # matched-random shortlists per (target, K)


def _readouts(W, rr, P, order):
    """argmin-free readouts on the nested prefix `order`; P is the pairwise matrix over it."""
    out = {}
    for m in MS:
        if m > len(order):
            continue
        sub = order[:m]
        Pm = P[:m, :m]
        b = int(np.argmin(Pm.mean(1)))
        out[m] = {
            "medoid": float(rr[sub[b]]),
            "oracle": float(rr[sub].min()),
            "mean": float(rr[sub].mean()),
        }
        if m == 1:
            out[m]["avg"] = float(rr[sub[0]])
            out[m]["div"] = 0.0
        else:
            A = I.superpose_batch(W[sub], W[sub[b]])
            C = A.mean(0)
            out[m]["avg"] = float(I.ca_rmsd(C, _readouts.nat))
            n = W.shape[1]
            out[m]["div"] = float(np.sqrt(np.mean(((A - C) ** 2).sum(axis=(1, 2)) / n)))
    return out


def one_target(t):
    pdb = t["pdb"]
    u = I.load_univ(pdb)
    o = np.asarray(u["order"], int)
    Wall = np.asarray(u["W"], float)[o]
    rrall = np.asarray(u["rr"], float)[o]
    nat = np.asarray(u["nat_ca"], float)
    _readouts.nat = nat
    U = len(o)

    p = L.pack(pdb, min(max(KS), U))
    dg_i, dg_j = p["i"], p["j"]
    rows_p = np.arange(len(dg_i))[None, :]

    rng = SD.stable_rng(pdb, "s17readout")
    cells = {}
    for K in KS:
        k = min(K, U)
        W = Wall[:k]
        rr = rrall[:k]
        D = np.linalg.norm(W[:, dg_i, :] - W[:, dg_j, :], axis=-1)
        g = L.gbin(D, p["grid"])
        sc = (p["risk_raw"][rows_p, g] * p["w_ship"][None, :]).mean(1)

        top = np.argsort(sc, kind="stable")[:min(MMAX, k)]
        P = I.pairwise_rmsd(W[top])
        cell = {"k": int(k), "argmin": L.sel_of(sc, rr), "oracle_K": float(rr.min()),
                "random_K": float(rr.mean()), "score": _readouts(W, rr, P, top)}

        acc = []
        for _ in range(NRAND):
            pick = rng.choice(k, min(MMAX, k), replace=False)
            acc.append(_readouts(W, rr, I.pairwise_rmsd(W[pick]), pick))
        cell["rand"] = {m: {q: float(np.mean([a[m][q] for a in acc if m in a]))
                            for q in acc[0][min(acc[0].keys())]}
                        for m in acc[0]}
        cells[str(K)] = cell
    return {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]), "U": int(U), "cells": cells}


def run(targets=None, verbose=True):
    tg = targets if targets is not None else L.targets()
    rows, t0 = [], time.time()
    for q, t in enumerate(tg):
        rows.append(one_target(t))
        if verbose and (q + 1) % 10 == 0:
            print(f"  {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(L.RESULTS, "sel_readout.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(L.RESULTS, "sel_readout.json"), "w"))
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(L.RESULTS, "sel_readout.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    print(f"\n{'='*104}\nTHE READOUT TRIANGLE   n = {len(rows)} targets, identical candidate sets, "
          f"CA-only RMSD\n{'='*104}")
    print("Every row: the SAME shortlist consumed three ways.  ORACLE = best member of the")
    print("shortlist, a ceiling.  'rand' columns use a matched random shortlist of the same size.\n")

    for K in KS:
        k = str(K)
        avail = [r for r in rows if k in r["cells"]]
        arg = np.array([r["cells"][k]["argmin"] for r in avail])
        okk = np.array([r["cells"][k]["oracle_K"] for r in avail])
        rkk = np.array([r["cells"][k]["random_K"] for r in avail])
        print(f"--- K = {K}   ORACLE over the whole set {okk.mean():.3f}   "
              f"random-of-one {rkk.mean():.3f}   distance argmin {arg.mean():.3f}")
        print(f"  {'m':>5}{'ORACLE@m':>10}{'medoid':>9}{'average':>9}{'div':>7}"
              f"{'| rand medoid':>14}{'rand avg':>10}{'| medoid-avg [95% CI]':>26}{'W/L':>8}")
        for m in MS:
            sub = [r for r in avail if str(m) in {str(x) for x in r["cells"][k]["score"]}]
            sub = [r for r in avail if m in r["cells"][k]["score"] or str(m) in r["cells"][k]["score"]]
            if not sub:
                continue
            gg = lambda r, d, q: (d[m] if m in d else d[str(m)])[q]      # noqa: E731
            om = np.array([gg(r, r["cells"][k]["score"], "oracle") for r in sub])
            md = np.array([gg(r, r["cells"][k]["score"], "medoid") for r in sub])
            av = np.array([gg(r, r["cells"][k]["score"], "avg") for r in sub])
            dv = np.array([gg(r, r["cells"][k]["score"], "div") for r in sub])
            rmd = np.array([gg(r, r["cells"][k]["rand"], "medoid") for r in sub])
            rav = np.array([gg(r, r["cells"][k]["rand"], "avg") for r in sub])
            f = np.array([r["fold"] for r in sub], int)
            pr = L.report_pair("", md, av, f)
            print(f"  {m:>5}{om.mean():>10.3f}{md.mean():>9.3f}{av.mean():>9.3f}{dv.mean():>7.3f}"
                  f"{rmd.mean():>14.3f}{rav.mean():>10.3f}"
                  f"   {pr['diff']:+.3f} [{pr['lo']:+.3f},{pr['hi']:+.3f}]{pr['W']:>5}/{pr['L']}")
        print()

    print("THE TRIANGLE, at the best m for each readout (chosen LEAVE-FOLD-OUT, not in-sample)")
    for K in KS:
        k = str(K)
        avail = [r for r in rows if k in r["cells"]]
        f = np.array([r["fold"] for r in avail], int)
        gg = lambda r, d, m, q: (d[m] if m in d else d[str(m)])[q]      # noqa: E731
        ms = [m for m in MS if all((m in r["cells"][k]["score"]) or (str(m) in r["cells"][k]["score"])
                                   for r in avail)]
        for tag in ("medoid", "avg"):
            M = np.array([[gg(r, r["cells"][k]["score"], m, tag) for m in ms] for r in avail])
            lfo, ch, ins = L.lfo_variant(M, f)
            print(f"  K={K:>5} {tag:<7} LFO m {dict((a, ms[b]) for a, b in ch.items())}  "
                  f"-> {lfo.mean():.3f}   (in-sample best m = {ms[ins]} at {M[:, ins].mean():.3f})")
        arg = np.array([r["cells"][k]["argmin"] for r in avail])
        print(f"  K={K:>5} argmin                                             -> {arg.mean():.3f}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
        report()
