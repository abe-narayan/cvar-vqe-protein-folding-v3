"""s17/avgceil.py -- THE CEILING OF THE READOUT THAT ACTUALLY WINS.

Every ceiling the programme has ever quoted is a **best-member** ceiling: 2.104 A in the shipped
top-75, 1.711 A at K = 500, 1.313 A over the full universe.  Those bound a SELECTION readout.

But selection is not the readout that wins.  The SELECT workstream resolved the readout triangle
at n = 126 on identical shortlists: the coordinate average (3.056) beats the medoid (3.282) beats
the distance argmin (3.454), the average dominating in **39 of 39 (K, m) cells** with every
interval excluding zero.  And the sprint has closed selection at three independent levels.

**So the programme has been quoting the ceiling of the wrong operator.**

THE POINT, and it is not a technicality.  A coordinate average can be BETTER THAN ANY OF ITS
MEMBERS -- that is the entire reason averaging works, and the Krogh-Vedelsby identity says so
exactly:

    readout^2 = mean member error^2 - diversity^2

The diversity term is subtracted.  With members whose errors partly cancel, the average lands
below the best member.  So the ORACLE-best SUBSET AVERAGE is not bounded by the best-member
ceiling and could be far below 1.313 A.  **Nobody has measured it**, and it is the ceiling that
governs the architecture the evidence actually supports.

WHAT THIS DECIDES.  If the best achievable subset average is far below the best-member ceiling,
then the sprint's target should be *choose a better subset to average*, which is a Problem C
question with a different and possibly easier signal requirement than ranking -- a subset
criterion need not order candidates, only partition them.  If it is no better than the
best-member ceiling, averaging has no hidden headroom and the programme's ceiling is what it
has always been.

METHOD.  Greedy forward selection under the ORACLE: start from the single best member, then
repeatedly add whichever remaining candidate most reduces the coordinate average's RMSD, up to
m_max.  Greedy is a LOWER BOUND on the true optimum -- the real ceiling can only be better --
and that direction is the safe one for a ceiling claim.  The running average is maintained
incrementally so each candidate trial costs one superposition rather than a re-average.

CONTROLS, both mandatory.

  greedy_random   the identical greedy procedure driven by a RANDOM criterion instead of the
                  native, which prices the SEARCH rather than the information.  A greedy
                  forward search over K candidates will improve a set even with no signal, and
                  without this control the whole result is a min-of-N artefact -- exactly the
                  trap that already cost this sprint one retracted ceiling (ledger L8).
  score_topm      the incumbent behaviour at the same m, so the gap is legible.

PRE-REGISTRATION.

  HYPOTHESIS.  The ORACLE subset average lands materially below the best-member ceiling,
  because the diversity term is subtracted and errors partly cancel.

  SUCCESS.  ORACLE subset average beats the best-member ceiling at the same K, target as the
  unit, interval excluding zero -- AND beats `greedy_random` by more than `greedy_random`
  beats `score_topm`, so the gain is information rather than search.

  FALSIFIER.  If the ORACLE subset average is no better than the best-member ceiling, averaging
  has no hidden headroom and this direction is closed.  If it beats the best-member ceiling but
  not `greedy_random`, the apparent headroom is a search artefact and must be reported as one.

ORACLE USE.  This module reads the native to drive the greedy criterion.  Everything it
produces is a CEILING and none of it is deployable or may be quoted as a result.  It exists to
tell the architecture where to aim.
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

KS = (500, 2000)
M_MAX = 40
POOL_CAP = 500          # candidates considered at each greedy step, by score, for tractability


def _greedy(W, d, nat, order_pool, m_max, criterion, rng=None):
    """forward greedy on the coordinate average.  `criterion` is 'oracle' or 'random'.

    Members are superposed onto the current medoid-free reference (the first pick) so the
    running mean is well defined; this is the same frame convention the production operator
    uses when it superposes on a fixed member.
    """
    ref = int(order_pool[int(np.argmin(d[order_pool]))]) if criterion == "oracle" \
        else int(order_pool[int(rng.integers(len(order_pool)))])
    sel = [ref]
    Wr = I.superpose_batch(W[order_pool], W[ref])          # all candidates in the ref frame
    pos = {int(c): q for q, c in enumerate(order_pool)}
    acc = Wr[pos[ref]].copy()
    traj = [float(I.ca_rmsd(acc, nat))]
    used = {ref}
    for _ in range(1, m_max):
        best, bestv = None, None
        for c in order_pool:
            c = int(c)
            if c in used:
                continue
            cand = (acc * len(sel) + Wr[pos[c]]) / (len(sel) + 1)
            v = rng.random() if criterion == "random" else float(I.ca_rmsd(cand, nat))
            if bestv is None or v < bestv:
                bestv, best = v, c
        if best is None:
            break
        acc = (acc * len(sel) + Wr[pos[best]]) / (len(sel) + 1)
        sel.append(best); used.add(best)
        traj.append(float(I.ca_rmsd(acc, nat)))
    return sel, traj


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
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        d = I.kabsch_rmsd_batch(W, nat)                    # ORACLE
        rng = SD.stable_rng(pdb, "s17avgceil")

        cell = {}
        for K in KS:
            k = min(K, U)
            #: the greedy pool is the score's top POOL_CAP, so the ceiling is one a deployable
            #: shortlist could in principle contain rather than an unreachable global optimum.
            pool = np.argsort(sc[:k], kind="stable")[:min(POOL_CAP, k)]
            _s, tro = _greedy(W, d, nat, pool, M_MAX, "oracle", rng)
            _s2, trr = _greedy(W, d, nat, pool, M_MAX, "random", rng)
            rkm = np.argsort(sc[:k], kind="stable")
            topm = {}
            for m in (5, 10, 20, 40):
                sel = rkm[:m]
                P = I.pairwise_rmsd(W[sel])
                avg, _b = I.coordinate_average(W[sel], P)
                topm[str(m)] = float(I.ca_rmsd(avg, nat))
            cell[str(K)] = {"best_member": float(d[:k].min()),
                            "oracle_traj": tro, "random_traj": trr, "score_topm": topm}
        rows.append({"pdb": pdb, "n": n, "fold": fold, "cells": cell})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(RESULTS, "avgceil.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "avgceil.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "avgceil.json")))["rows"]
    rng = SD.stable_rng("avgceil", "report")
    print(f"\nn = {len(rows)} targets.  EVERY arm here is an ORACLE CEILING except `score_topm`.")
    print("Greedy is a LOWER bound on the true optimum, so the real ceiling can only be better.\n")
    for K in KS:
        k = str(K)
        if k not in rows[0]["cells"]:
            continue
        bm = np.array([r["cells"][k]["best_member"] for r in rows])
        L = min(len(r["cells"][k]["oracle_traj"]) for r in rows)
        tro = np.array([r["cells"][k]["oracle_traj"][:L] for r in rows])
        trr = np.array([r["cells"][k]["random_traj"][:L] for r in rows])
        print(f"=== K = {K}   best-member ceiling {bm.mean():.3f} A")
        print(f"  {'m':>4}{'ORACLE subset avg':>20}{'greedy-RANDOM ctrl':>21}"
              f"{'score top-m':>14}")
        tm = rows[0]["cells"][k]["score_topm"]
        for m in range(1, L):
            sm = ""
            if str(m + 1) in tm:
                sm = f"{np.mean([r['cells'][k]['score_topm'][str(m+1)] for r in rows]):>14.3f}"
            print(f"  {m+1:>4}{tro[:, m].mean():>20.3f}{trr[:, m].mean():>21.3f}{sm}")
        j = int(np.argmin(tro.mean(0)))
        mu, lo, hi = _boot(tro[:, j] - bm, rng)
        mu2, lo2, hi2 = _boot(tro[:, j] - trr[:, j], rng)
        w = int((tro[:, j] < bm).sum()); l = int((tro[:, j] > bm).sum())
        print(f"\n  best ORACLE subset average at m = {j+1}: {tro[:, j].mean():.3f} A")
        print(f"    vs the best-member ceiling   {mu:+.3f} [{lo:+.3f},{hi:+.3f}]  {w}W/{l}L")
        print(f"    vs the greedy-RANDOM control {mu2:+.3f} [{lo2:+.3f},{hi2:+.3f}]")
        print(f"    (the second line is the one that decides whether the headroom is")
        print(f"     INFORMATION or merely the min-of-N of a forward search)\n")
    print("READ.  If the ORACLE subset average sits well below the best-member ceiling AND well")
    print("below its own random-greedy control, then the architecture should aim at CHOOSING A")
    print("SUBSET TO AVERAGE rather than at ranking candidates -- a subset criterion need not")
    print("order candidates, only partition them, which is a weaker requirement than the in-band")
    print("ranking the sprint has closed at three levels.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
