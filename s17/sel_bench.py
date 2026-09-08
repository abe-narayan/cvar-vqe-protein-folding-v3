"""s17/sel_bench.py -- THE IDENTICAL-CANDIDATE-SET INSTRUMENT.  BRIEF section 3B, sprint section 10, PREREG E1.

Every selector in this programme has been priced on a different candidate set at some point,
which makes the comparisons void.  This module fixes one candidate set per (target, K) and
runs every arm through it.  Nothing here may consume anything the other arms cannot see.

THE ARMS, all on the identical set:

    ORACLE            best member                                   [CEILING, never a result]
    dist              argmin of the shipped Bayes-risk objective
    legacy            argmin of the genuine 11-component Legacy total
    legacy_gate_dist  Legacy drops the worst q of the set, dist ranks the survivors
    rand_gate_dist    THE CONTROL FOR THE GATE: drop the same COUNT at random, dist ranks
    consensus         S8's deployable arm: medoid of dist's own top-m
    rand_consensus    THE CONTROL FOR CONSENSUS: medoid of a random m of the same size
    random            matched random selection of one -- EXACT expectation, noise-free
    bestof_m          random best-of-m: the control for any arm allowed to look at m members

WHY EVERY GATE ARM CARRIES ITS OWN RANDOM GATE.  Removing candidates changes the argmin
distribution even when the removal is uninformative -- a smaller set has a worse expected
minimum under a noisy score and a better one under a good score.  Sprint 16's G3 measured
exactly this: Legacy's filter looked worth +0.024 A and its matched random filter was worth
+0.012 A of the same.  A gate priced against "no gate" is priced against the wrong thing.

RHO IS REPORTED FOR EVERY ARM AT EVERY K.  The sprint's design equation
(`s17/selection_theory.py`) consumes the selector's rank correlation with the true ordering,
not its RMSD; RMSD alone cannot distinguish a selector that is on a path to the target from
one that is not.  Three rho are reported and they are different quantities:

    rho_within   mean per-target Spearman(score, true RMSD).  THIS is the copula's parameter.
    rho_band     the same restricted to the near-native band (best + 1.5 A).  This is the
                 quantity comparable with the programme's recorded 0.600 / 0.638.
    rho_pool     Spearman over all (target, candidate) pairs pooled.  Measures whether the
                 score is comparable BETWEEN targets -- irrelevant to argmin-per-target, and
                 decisive for anything that thresholds a score globally.

THE COPULA GAP, the diagnostic the coordinator asked for.  Given a selector's measured
rho_within and the target's own member distribution, the design equation predicts what that
rho should deliver under exchangeable errors.  `gap = realized - predicted`.  A positive gap
means the selector is worse than its own rho implies, i.e. its errors are STRUCTURED rather
than exchangeable -- a different and more fixable failure than low rho.
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

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from s17 import sel_lib as L              # noqa: E402
from s17 import sel_obj as O              # noqa: E402
from s17 import selection_theory as TH    # noqa: E402

KS = (75, 500, 2000, 0)          # 0 = the whole universe
K_LEGACY = 500                   # Legacy's component builder is the expensive part
GATES = (0.25, 0.50, 0.75)
M_CONS = 75                      # the S8 deployable arm's filter width
NREP = 64                        # repeats for every matched-random control


def _kname(K):
    return "full" if K == 0 else str(K)


def _consensus(W):
    """Medoid of a candidate set: S8's consensus readout, tie-safe by construction."""
    P = I.pairwise_rmsd(W)
    return int(np.argmin(P.mean(1)))


def one_target(t, verbose=False):
    pdb, n, fold, seq = t["pdb"], int(t["n"]), int(t["fold"]), t["seq"]
    u = I.load_univ(pdb)
    order = np.asarray(u["order"], int)
    U = len(order)
    Wall = np.asarray(u["W"], float)[order]
    rrall = np.asarray(u["rr"], float)[order]
    PHI = np.asarray(u["PHI"], float)[order]
    PSI = np.asarray(u["PSI"], float)[order]

    dg = I.distogram(pdb, seq, fold)
    i, j = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
    prob = np.asarray(dg["prob"], float)
    centres = np.asarray(dg["centres"], float)
    grid = np.asarray(dg["grid"], float)
    risk_raw = (prob[:, None, :] * np.abs(grid[None, :, None] - centres[None, None, :])).sum(2)
    wship = L._wship(np.asarray(dg["sd"], float))
    rows_p = np.arange(len(i))[None, :]

    Dall = np.linalg.norm(Wall[:, i, :] - Wall[:, j, :], axis=-1)
    g = L.gbin(Dall, grid)
    score_all = (risk_raw[rows_p, g] * wship[None, :]).mean(1)

    #: Legacy on the first K_LEGACY candidates only
    from s16 import energy_lib as EL
    kL = min(K_LEGACY, U)
    comp = EL.legacy_components_of_windows(seq, PHI[:kL], PSI[:kL])
    leg_all = np.asarray(EL.legacy_total_from(comp), float)

    rng = SD.stable_rng(pdb, "s17bench")
    out = {"pdb": pdb, "n": n, "fold": fold, "U": int(U), "cells": {}}
    for K in KS:
        k = U if K == 0 else min(K, U)
        rr = rrall[:k]
        sc = score_all[:k]
        e = {"k": int(k), "oracle": float(rr.min()), "mean": float(rr.mean()),
             "median": float(np.median(rr)), "random": float(rr.mean())}
        e["dist"] = L.sel_of(sc, rr)
        e["rho_dist"] = O.spearman(sc, rr)
        m = rr <= rr.min() + 1.5
        e["rho_band_dist"] = O.spearman(sc[m], rr[m]) if m.sum() >= 8 else None
        e["frac_band"] = float(m.mean())

        # --- best-of-m random: the control for any arm allowed to look at m members
        for mm in (5, 25, 75):
            e[f"bestof{mm}"] = L.topm_random(rr, mm, pdb, f"bo{mm}k{K}", n=NREP)

        # --- consensus: medoid of dist's own top-m, and its matched random control
        if k >= M_CONS:
            top = np.argsort(sc, kind="stable")[:M_CONS]
            e["consensus"] = float(rr[top[_consensus(Wall[:k][top])]])
            acc = []
            for _ in range(8):
                pick = rng.choice(k, M_CONS, replace=False)
                acc.append(rr[pick[_consensus(Wall[:k][pick])]])
            e["rand_consensus"] = float(np.mean(acc))
            e["top75_oracle"] = float(rr[top].min())
            e["top75_mean"] = float(rr[top].mean())

        # --- Legacy arms, only where Legacy has been computed
        if k <= kL:
            leg = leg_all[:k]
            e["legacy"] = L.sel_of(leg, rr)
            e["rho_legacy"] = O.spearman(leg, rr)
            e["rho_band_legacy"] = O.spearman(leg[m], rr[m]) if m.sum() >= 8 else None
            for q in GATES:
                keep = max(2, int(round((1.0 - q) * k)))
                srv = np.argsort(leg, kind="stable")[:keep]           # Legacy keeps the best
                e[f"legacy_gate{int(q*100)}"] = L.sel_of(sc[srv], rr[srv])
                acc = []
                for _ in range(NREP):
                    p2 = rng.choice(k, keep, replace=False)
                    acc.append(L.sel_of(sc[p2], rr[p2]))
                e[f"rand_gate{int(q*100)}"] = float(np.mean(acc))
        out["cells"][_kname(K)] = e
    return out


def run(targets=None, verbose=True):
    tg = targets if targets is not None else L.targets()
    rows, t0 = [], time.time()
    for q, t in enumerate(tg):
        rows.append(one_target(t))
        if verbose and (q + 1) % 10 == 0:
            print(f"  {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(L.RESULTS, "sel_bench.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(L.RESULTS, "sel_bench.json"), "w"))
    return rows


# --------------------------------------------------------------------------- reporting
def _col(rows, k, f):
    return np.array([r["cells"][k][f] for r in rows if r["cells"][k].get(f) is not None], float)


def _both(rows, k, f, ref):
    """paired vectors, restricted to targets where both exist"""
    a, b = [], []
    for r in rows:
        c = r["cells"][k]
        if c.get(f) is not None and c.get(ref) is not None:
            a.append(c[f]); b.append(c[ref])
    return np.array(a, float), np.array(b, float)


def _fold(rows, k, f, ref):
    return np.array([r["fold"] for r in rows
                     if r["cells"][k].get(f) is not None and r["cells"][k].get(ref) is not None], int)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(L.RESULTS, "sel_bench.json")))["rows"]
    n = len(rows)
    print(f"\n{'='*104}\nTHE IDENTICAL-CANDIDATE-SET INSTRUMENT   n = {n} targets\n{'='*104}")
    print("Every arm in a row consumes THE SAME structures.  ORACLE reads the native and is a")
    print("ceiling, never a result.  'random' is the EXACT expectation of selecting one member.\n")

    print("A. ORACLE CEILING / REALIZED / GAP, by candidate width")
    print(f"  {'K':>6}{'ORACLE':>9}{'dist':>9}{'gap':>8}{'legacy':>9}{'consensus':>11}"
          f"{'random':>9}{'bestof75':>10}{'rho_w':>8}{'rho_b':>8}")
    for K in KS:
        k = _kname(K)
        o = _col(rows, k, "oracle"); d = _col(rows, k, "dist"); rd = _col(rows, k, "random")
        lg = _col(rows, k, "legacy"); cs = _col(rows, k, "consensus")
        b75 = _col(rows, k, "bestof75")
        rw = _col(rows, k, "rho_dist"); rb = _col(rows, k, "rho_band_dist")
        print(f"  {k:>6}{o.mean():>9.3f}{d.mean():>9.3f}{d.mean()-o.mean():>8.3f}"
              + (f"{lg.mean():>9.3f}" if len(lg) else f"{'--':>9}")
              + (f"{cs.mean():>11.3f}" if len(cs) else f"{'--':>11}")
              + f"{rd.mean():>9.3f}{b75.mean():>10.3f}{rw.mean():>8.3f}{rb.mean():>8.3f}")

    print("\nB. EVERY ARM AGAINST ITS OWN MATCHED CONTROL   (negative = arm is better)")
    for K in KS:
        k = _kname(K)
        print(f"  --- K = {k} ---")
        pairs = [("dist", "random"), ("dist", "bestof5"), ("legacy", "random"),
                 ("consensus", "rand_consensus"), ("consensus", "dist"),
                 ("legacy", "dist")]
        for q in GATES:
            pairs.append((f"legacy_gate{int(q*100)}", f"rand_gate{int(q*100)}"))
            pairs.append((f"legacy_gate{int(q*100)}", "dist"))
        for f, ref in pairs:
            a, b = _both(rows, k, f, ref)
            if len(a) == 0:
                continue
            print("    " + L.fmt_pair(L.report_pair(f"{f} vs {ref}", a, b, _fold(rows, k, f, ref))))

    print("\nC. DOES WIDENING K REACH THE SELECTOR?   ceiling gained vs realized gained")
    k0 = "75"
    o0, d0 = _col(rows, k0, "oracle"), _col(rows, k0, "dist")
    for K in KS[1:]:
        k = _kname(K)
        o, d = _col(rows, k, "oracle"), _col(rows, k, "dist")
        dc = o.mean() - o0.mean()
        r = L.report_pair("", d, d0, np.array([x["fold"] for x in rows], int))
        rec = 0.0 if abs(dc) < 1e-9 else r["diff"] / dc
        print(f"  K={k:>5}  ceiling {dc:+.3f}   realized {r['diff']:+.3f} "
              f"[{r['lo']:+.3f},{r['hi']:+.3f}]  {r['W']}W/{r['L']}L   recovered {rec:>6.1%}")

    print("\nD. RHO AS A FUNCTION OF K  -- is the ranker's skill stable as the set widens?")
    print(f"  {'K':>6}{'rho_within':>12}{'rho_band':>10}{'frac in band':>14}{'rho_pool':>10}")
    for K in KS:
        k = _kname(K)
        rw = _col(rows, k, "rho_dist"); rb = _col(rows, k, "rho_band_dist")
        fb = _col(rows, k, "frac_band")
        print(f"  {k:>6}{rw.mean():>12.3f}{rb.mean():>10.3f}{fb.mean():>14.3f}{'n/a':>10}")

    print("\nE. THE COPULA GAP  -- realized minus what the measured rho should deliver")
    print("   positive = the selector is worse than its own rho implies, i.e. its errors are")
    print("   STRUCTURED rather than exchangeable.")
    rng = SD.stable_rng("s17bench", "copula")
    for K in KS:
        k = _kname(K)
        gaps, preds, reals = [], [], []
        for r in rows:
            c = r["cells"][k]
            if c.get("rho_dist") is None:
                continue
            u = I.load_univ(r["pdb"])
            o = np.asarray(u["order"], int)
            d = np.asarray(u["rr"], float)[o][:c["k"]]
            z = TH._norm_scores(d)
            pr = TH._expected(d, z, max(0.0, float(c["rho_dist"])), rng, ndraw=200)[0]
            preds.append(pr); reals.append(c["dist"]); gaps.append(c["dist"] - pr)
        g = np.array(gaps)
        print(f"  K={k:>5}  realized {np.mean(reals):>7.3f}   copula-predicted "
              f"{np.mean(preds):>7.3f}   gap {g.mean():+7.3f}  "
              f"median {np.median(g):+.3f}  {(g>0).sum()}/{len(g)} above curve")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
        report()
