"""s17/coverage.py -- COVERAGE-PRESERVING SHORTLISTS.  Problem C, and the lever L12 exposed.

WHY THIS EXISTS.  The SELECT workstream isolated, at n = 126, why widening the candidate set
hurts: as K grows the score's own top-75 shortlist gets **better on average** (-0.730
[-0.906, -0.560]) and **worse at its best** (+0.468 [+0.304, +0.642]).  The extra windows look
more plausible without being nearer the native, and they DISPLACE the genuinely near-native
members out of the shortlist.  The shortlist-ceiling destruction (0.468 A) exceeds the realized
degradation (0.141 A), so it accounts for the whole effect.

The consequence is sharp: **the full-universe ceiling of 1.313 A is unreachable through a
score-ranked top-B shortlist at ANY ranking skill**, because the shortlist ceiling there is
2.572 A.  No amount of Problem-B ranking work can recover something the shortlist no longer
contains.  That makes this a Problem C question -- ensemble construction -- and it is the one
the sprint has not attacked.

THE HYPOTHESIS.  A shortlist chosen for COVERAGE of the candidate set's structural variety,
rather than purely for score, retains near-native members that a score ranking discards, and so
raises the shortlist ceiling at the same B.

THIS IS NOT A REPEAT OF SPRINT 16's DIVERSITY NULL, and the distinction is exact.  Sprint 16's
`divselect` showed diversity-aware selection is a null for the AVERAGING readout, because that
readout obeys `readout^2 = member error^2 - diversity^2` and the two terms move together under
selection pressure, cancelling at par.  **That identity does not govern an argmin readout.**
Here the quantity of interest is the shortlist's MINIMUM, not its mean, and preserving coverage
protects the minimum precisely where averaging gains nothing.  Same word, different functional.

METHOD.  Diversity is measured in PAIR-DISTANCE space rather than by pairwise RMSD, which makes
it O(K*B*npairs) instead of O(K^2) superpositions and therefore usable at the full universe --
and, not incidentally, makes any resulting shortlist deployable.

ARMS at matched B, all consuming the same candidate set:

    score_topB     top-B by the shipped distance objective         the incumbent behaviour
    random_B       B drawn uniformly                               CONTROL (coverage without score)
    maxmin         farthest-point sampling seeded by the best      coverage, score only as seed
    maxmin_gated   farthest-point sampling restricted to the       coverage INSIDE a score prefix
                   score's top-P (P >> B)
    strata         equal counts from each of B score bins          the cheapest coverage proxy
    kmeans_best    k-means into B cells, best-scoring per cell     coverage + score, per cell

MEASURED for each arm: the shortlist ORACLE best (the ceiling -- the quantity L12 says is being
destroyed), the shortlist mean, and two realized readouts (score argmin, consensus medoid).

PRE-REGISTRATION.

  SUCCESS.  Some coverage arm raises the shortlist ORACLE ceiling above `score_topB` at matched
  B, with a paired fold-clustered interval excluding zero, target as the unit.

  FALSIFIER.  If no coverage arm beats `score_topB` on the ceiling, the near-native members are
  not merely mis-ranked but genuinely unreachable by any B-sized shortlist of this candidate
  set, and Problem C is closed by this route.

  THE TRAP TO AVOID.  A ceiling gain is worthless if the realized readout cannot use it.  Both
  are reported side by side, and a coverage arm that raises the ceiling while degrading the
  realized readout is reported as what it is -- a wider net that the selector cannot exploit,
  which would say the two problems must be solved together or not at all.

ORACLE USE.  Shortlist ceilings read the native and are labelled ORACLE.  Every shortlist
construction is native-free.
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

KS = (500, 2000, 0)          # 0 = the full universe
BS = (25, 75, 200)
GATE_MULT = 8                # maxmin_gated searches the score's top GATE_MULT*B


def _farthest_point(D, B, seed_idx, cand=None):
    """farthest-point sampling in pair-distance space; returns B indices including the seed."""
    n = len(D)
    pool = np.arange(n) if cand is None else np.asarray(cand, int)
    sel = [int(seed_idx)]
    d2 = ((D[pool] - D[sel[0]]) ** 2).sum(1)
    while len(sel) < min(B, len(pool)):
        a = int(np.argmax(d2))
        nxt = int(pool[a])
        if nxt in sel:
            d2[a] = -1.0
            continue
        sel.append(nxt)
        d2 = np.minimum(d2, ((D[pool] - D[nxt]) ** 2).sum(1))
    return np.asarray(sel, int)


def _kmeans_best(D, score, B, rng, iters=12):
    """k-means in pair-distance space; keep the best-scoring member of each cell."""
    n = len(D)
    c = D[rng.choice(n, size=min(B, n), replace=False)].copy()
    lab = np.zeros(n, int)
    for _ in range(iters):
        d2 = ((D[:, None, :] - c[None, :, :]) ** 2).sum(2) if n * len(c) < 4_000_000 else None
        if d2 is None:
            lab = np.array([int(np.argmin(((c - D[q]) ** 2).sum(1))) for q in range(n)])
        else:
            lab = np.argmin(d2, axis=1)
        for k in range(len(c)):
            m = lab == k
            if m.any():
                c[k] = D[m].mean(0)
    out = []
    for k in range(len(c)):
        m = np.where(lab == k)[0]
        if len(m):
            out.append(int(m[int(np.argmin(score[m]))]))
    return np.asarray(out, int)


def _strata(score, B):
    """equal counts from each of B equal-width score bins -- the cheapest coverage proxy."""
    r = np.argsort(score, kind="stable")
    chunks = np.array_split(r, B)
    return np.asarray([int(c[0]) for c in chunks if len(c)], int)


def _readouts(Wsel, dsel, score_sel, nat_ref):
    """shortlist ceiling (ORACLE) and two native-free realized readouts."""
    out = {"ceiling": float(dsel.min()), "mean": float(dsel.mean()),
           "argmin": float(dsel[int(np.argmin(score_sel))])}
    if len(Wsel) > 1:
        P = I.pairwise_rmsd(Wsel)
        out["medoid"] = float(dsel[int(I.medoid(P))])
        #: THE COORDINATE AVERAGE.  Added after the SELECT workstream resolved the readout
        #: triangle at n = 126 on identical shortlists: average 3.056 < medoid 3.282 < argmin
        #: 3.454, with the average winning all 39 (K, m) cells and every interval excluding
        #: zero.  A coverage module that omits the readout that actually wins would be
        #: measuring the wrong operator, and averaging consumes the set MEAN -- exactly the
        #: quantity these coverage arms move most.
        avg, _b = I.coordinate_average(Wsel, P)
        out["average"] = float(I.ca_rmsd(avg, nat_ref[0]))
    else:
        out["medoid"] = out["argmin"]
        out["average"] = out["argmin"]
    return out


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
        Dall = I.pair_dists(W, i, j)
        score_all = np.asarray(I.shipped_score(dg, Dall), float)
        d_all = I.kabsch_rmsd_batch(W, nat)          # ORACLE
        rng = SD.stable_rng(pdb, "s17cover")

        cell = {}
        for K in KS:
            k = U if K == 0 else min(K, U)
            D, sc, d = Dall[:k], score_all[:k], d_all[:k]
            #: standardise the distance space so no single pair dominates the geometry
            Z = (D - D.mean(0)) / np.maximum(D.std(0), 1e-9)
            rk = np.argsort(sc, kind="stable")
            ck = {}
            for B in BS:
                if B > k:
                    continue
                arms = {}
                arms["score_topB"] = rk[:B]
                arms["random_B"] = rng.permutation(k)[:B]
                arms["maxmin"] = _farthest_point(Z, B, int(rk[0]))
                P = min(GATE_MULT * B, k)
                arms["maxmin_gated"] = _farthest_point(Z, B, int(rk[0]), cand=rk[:P])
                arms["strata"] = _strata(sc, B)
                arms["kmeans_best"] = _kmeans_best(Z, sc, B, rng)
                ck[str(B)] = {a: _readouts(W[s], d[s], sc[s], (nat,))
                              for a, s in arms.items()}
                ck[str(B)]["_set_ceiling"] = float(d.min())
            cell["full" if K == 0 else str(K)] = ck
        rows.append({"pdb": pdb, "n": n, "fold": fold, "U": int(U), "cells": cell})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(RESULTS, "coverage.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "coverage.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


ARMS = ("score_topB", "random_B", "maxmin", "maxmin_gated", "strata", "kmeans_best")


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "coverage.json")))["rows"]
    rng = SD.stable_rng("coverage", "report")
    print(f"\nn = {len(rows)} targets.  Shortlist CEILING is ORACLE.  `argmin` and `medoid` are")
    print("native-free realized readouts on the same shortlist.  Every shortlist is native-free.\n")
    ks = [k for k in rows[0]["cells"]]
    for K in ks:
        for B in BS:
            b = str(B)
            if b not in rows[0]["cells"][K]:
                continue
            setc = np.mean([r["cells"][K][b]["_set_ceiling"] for r in rows])
            print(f"=== K = {K}, B = {B}   full-set ORACLE ceiling {setc:.3f}")
            print(f"  {'shortlist':<16}{'CEILING':>9}{'mean':>8}{'argmin':>9}"
                  f"{'medoid':>9}{'AVERAGE':>10}"
                  f"   {'ceiling vs score_topB [95% CI]':>32}")
            base = np.array([r["cells"][K][b]["score_topB"]["ceiling"] for r in rows])
            for a in ARMS:
                cl = np.array([r["cells"][K][b][a]["ceiling"] for r in rows])
                mn = np.mean([r["cells"][K][b][a]["mean"] for r in rows])
                am = np.mean([r["cells"][K][b][a]["argmin"] for r in rows])
                md = np.mean([r["cells"][K][b][a]["medoid"] for r in rows])
                av = np.mean([r["cells"][K][b][a].get("average", np.nan) for r in rows])
                mu, lo, hi = _boot(cl - base, rng)
                w = int((cl < base).sum()); l = int((cl > base).sum())
                tag = "" if a != "score_topB" else "  <- incumbent"
                print(f"  {a:<16}{cl.mean():>9.3f}{mn:>8.3f}{am:>9.3f}{md:>9.3f}"
                      f"{av:>10.3f}   {mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}{tag}")
            print()
    print("READ.  A coverage arm that lowers the CEILING column has retained near-native members")
    print("the score ranking discards -- that is the L12 mechanism reversed.  But a ceiling gain")
    print("is worthless unless a realized readout can use it, so the argmin and medoid columns")
    print("decide whether Problem C alone is enough or whether B and C must be solved together.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
