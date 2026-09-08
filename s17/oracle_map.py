"""s17/oracle_map.py -- THE ORACLE MAP.  Where every angstrom actually is, before any modelling.

Sprint 17 section 7 requires a complete internal map of attainable accuracy BEFORE any
sophisticated selection model is built.  This module is that map, and it is the shared
instrument every other workstream in the sprint reads.

THE FACT THAT MOTIVATES IT.  The programme has been quoting a "K = 500 pool oracle" of
1.711 A as the retrieval ceiling.  That is not the retrieval ceiling -- it is the ceiling of
the first 500 entries of a BLOSUM62 ranking over a universe of 13,000-27,000 windows per
target.  Over the full universe the ceiling is 1.313 A.

**AND THE OBVIOUS READING OF THAT IS WRONG.  A1 (AUDIT).**  The minimum of N draws falls with
N even when nothing new is present, so the 0.397 A "truncation cost" needs a matched null.
Under a RANDOM PERMUTATION of the identical universe the same statistic is 0.473 A, so
observed minus null is **-0.075 [-0.132, -0.022]** -- the gain is significantly BELOW the
order-statistic null.  BLOSUM has already front-loaded the quality and the marginal return on
the next 12,500 windows is *worse* than i.i.d. draws.  The whole retrieval apparatus buys
**+1.5 targets** of 2.0 A recall over a random 500 of the same library.

So the ceilings below are correct and the deep universe is real structure rather than
near-duplicates (AUDIT measured d* = 1.355 A, a dense continuum) -- but "widening K reveals
better structures the ranking was hiding" is NOT what the numbers say.  It reveals **more
draws**, extracted slightly less efficiently than chance.  Note also corr(U, chain length)
= -0.988: the universe size is essentially 1/n, so any per-target exhibit is confounded with
chain length.

WHAT IS MEASURED, per target and per candidate width K:

    best / median / mean member RMSD          ORACLE -- ceilings and diagnosis only
    coordinate-average RMSD                   the incumbent's readout operator
    diversity (the Krogh-Vedelsby term)       NATIVE-FREE
    distance-objective-selected RMSD          NATIVE-FREE selector (the shipped objective)
    Legacy-selected RMSD                      NATIVE-FREE selector
    random-selected RMSD                      the control every selector must beat
    recall of < 1.0 / 1.5 / 2.0 / 2.5 A       ORACLE -- how often the answer is present at all

AMBER selection is deliberately NOT here: a single point is ~6 ms, so scoring 13k windows x
126 targets is ~2.7 CPU-hours and belongs to the physics workstream on a shortlist, not to
the map.

CONVENTIONS THAT ARE NOW STANDING LAW IN THIS PROGRAMME.

* Every table separates **ORACLE CEILING** from **REALIZED** and prints the **GAP**.  A 1.7 A
  oracle with a 3 A realized prediction is not a 1.7 A method.
* The unit of analysis is the TARGET.  Means are printed with medians, win/loss and
  fold-clustered bootstrap intervals wherever a comparison is made.
* Every selector is priced against **random selection of the same count from the same
  candidate set** -- not against a weaker selector.
* Diversity is measured as the Krogh-Vedelsby ambiguity term in the frame the averaging
  operator itself builds in, because that is the functional the identity consumes.  Sprint 16
  established that reading the free-superposition term instead changes the answer by ~2x.

COST NOTE.  The coordinate-average operator superposes on the set medoid, which is O(K^2).
It is therefore computed up to K = MAX_AVG and skipped above, with the omission printed rather
than hidden.  Everything else is O(K) or O(K log K) and runs over the full universe.
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
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

#: candidate widths.  0 is the sentinel for "the whole universe".
KS = (25, 50, 75, 100, 150, 300, 500, 1000, 2000, 5000, 0)
MAX_AVG = 2000          # coordinate averaging is O(K^2); above this it is skipped, not faked
THRESH = (1.0, 1.5, 2.0, 2.5)
N_RAND = 3              # random-selection draws, seeded through s15.seed


def _kname(K):
    return "full" if K == 0 else str(K)


def _diversity(Wsel, P=None):
    """The Krogh-Vedelsby ambiguity term, in the frame the averaging operator builds in.

    Returns (coordinate-average structure, diversity in A, medoid index).
    """
    if P is None:
        P = I.pairwise_rmsd(Wsel)
    avg, b = I.coordinate_average(Wsel, P)
    Wm = I.superpose_batch(Wsel, Wsel[b])
    n = Wsel.shape[1]
    div = float(np.sqrt(np.mean(((Wm - avg) ** 2).sum(axis=(1, 2)) / n)))
    return np.asarray(avg, float), div, int(b)


def run(targets=None, verbose=True):
    from s16 import energy_lib as EL
    tg = targets if targets is not None else I.targets()
    rows = []
    t0 = time.time()

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        order = np.asarray(u["order"], int)
        W = np.asarray(u["W"], float)[order]          # universe in BLOSUM-ranked order
        PHI = np.asarray(u["PHI"], float)[order]
        PSI = np.asarray(u["PSI"], float)[order]
        nat = np.asarray(u["nat_ca"], float)
        U = len(W)

        #: ORACLE -- per-window RMSD to the native.  Diagnosis and ceilings only.
        d = I.kabsch_rmsd_batch(W, nat)

        #: NATIVE-FREE -- the shipped distance objective on every window.
        i, j = I.pair_index(n)
        dg = I.distogram(pdb, seq, fold)
        score = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)

        rng = SD.stable_rng(pdb, "s17map")
        cell = {}
        for K in KS:
            k = U if K == 0 else min(K, U)
            dk = d[:k]
            e = {
                "k_actual": int(k),
                "best": float(dk.min()), "median": float(np.median(dk)),
                "mean": float(dk.mean()),
                "sel_dist": float(dk[int(np.argmin(score[:k]))]),
                "sel_rand": float(np.mean([dk[int(rng.integers(k))] for _ in range(N_RAND)])),
            }
            for th in THRESH:
                e[f"recall_{th}"] = float((dk < th).mean())
                e[f"has_{th}"] = bool((dk < th).any())
            if k <= MAX_AVG:
                avg, div, b_ = _diversity(W[:k])
                e["coord_avg"] = float(I.ca_rmsd(avg, nat))
                e["diversity"] = div
                #: A3 (AUDIT).  The identity readout^2 = <err^2> - div^2 consumes the
                #: COMMON-FRAME quadratic member error.  Table A previously printed the
                #: free-superposition arithmetic `mean` beside the common-frame `diversity`,
                #: which breaks the identity by ~18% of readout^2.  Both are now stored and
                #: the identity residual is printed so the two can never be read as one.
                Wm = I.superpose_batch(W[:k], W[:k][b_])
                natm = I.superpose_batch(np.asarray(nat, float)[None], avg)[0]
                e["err_cf"] = float(np.sqrt(np.mean(((Wm - natm) ** 2).sum(axis=(1, 2)) / n)))
                e["identity_resid"] = float(abs(np.sqrt(max(e["err_cf"] ** 2 - div * div, 0.0))
                                                - e["coord_avg"]))
            cell[_kname(K)] = e

        #: Legacy on the top-500 only.  The eleven-component model is ~O(K) but the
        #: component builder is the expensive part, and Legacy has never been shown to rank;
        #: 500 is enough to price it as a selector against its own random control.
        kL = min(500, U)
        comp = EL.legacy_components_of_windows(seq, PHI[:kL], PSI[:kL])
        leg = np.asarray(EL.legacy_total_from(comp), float)
        for K in KS:
            k = U if K == 0 else min(K, U)
            if k <= kL:
                cell[_kname(K)]["sel_legacy"] = float(d[:k][int(np.argmin(leg[:k]))])

        rows.append({"pdb": pdb, "n": n, "fold": fold, "U": int(U), "cells": cell})
        if verbose and (c + 1) % 20 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False, "n_expected": len(tg)},
                      open(os.path.join(RESULTS, "oracle_map.json"), "w"))

    #: A10 (AUDIT).  A partial artefact was read by three workstreams at 60/126 rows with fold 3
    #: under-represented 2.7x.  The completion flag makes an incomplete read detectable.
    json.dump({"rows": rows, "complete": len(rows) == len(tg), "n_expected": len(tg)},
              open(os.path.join(RESULTS, "oracle_map.json"), "w"))
    report(rows)
    return rows


# ------------------------------------------------------------------ reporting
def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "oracle_map.json")))["rows"]
    rng = SD.stable_rng("oraclemap", "report")
    ks = [_kname(K) for K in KS]
    g = lambda k, f: np.array([r["cells"][k][f] for r in rows if f in r["cells"][k]])  # noqa: E731
    U = np.array([r["U"] for r in rows])

    print(f"\nn = {len(rows)} targets.  Universe size: min {U.min()}  median "
          f"{int(np.median(U))}  max {U.max()}.")
    print("Every RMSD in angstroms.  ORACLE columns read the native and are ceilings, never results.\n")

    print("A. THE CEILING AS A FUNCTION OF CANDIDATE WIDTH   [ORACLE]")
    print(f"  {'K':>6}{'best':>9}{'median':>9}{'mean':>9}{'coordavg':>10}{'diversity':>11}"
          f"{'<2.0 A present':>16}{'<2.5 A present':>16}")
    for k in ks:
        b = g(k, "best")
        ca = g(k, "coord_avg"); dv = g(k, "diversity")
        h2 = np.array([r["cells"][k]["has_2.0"] for r in rows]).mean()
        h25 = np.array([r["cells"][k]["has_2.5"] for r in rows]).mean()
        print(f"  {k:>6}{b.mean():>9.3f}{g(k,'median').mean():>9.3f}{g(k,'mean').mean():>9.3f}"
              + (f"{ca.mean():>10.3f}{dv.mean():>11.3f}" if len(ca) else f"{'--':>10}{'--':>11}")
              + f"{h2:>15.1%}{h25:>16.1%}")

    print("\nB. REALIZED SELECTORS vs THE ORACLE CEILING, and vs RANDOM SELECTION")
    print(f"  {'K':>6}{'ORACLE':>9}{'dist':>9}{'gap':>8}{'legacy':>9}{'random':>9}"
          f"{'   dist - random [95% CI]':>28}{'  W/L':>8}")
    for k in ks:
        b = g(k, "best"); sd_ = g(k, "sel_dist"); sr = g(k, "sel_rand")
        sl = g(k, "sel_legacy")
        mu, lo, hi = _boot(sd_ - sr, rng)
        w = int((sd_ < sr).sum()); l = int((sd_ > sr).sum())
        print(f"  {k:>6}{b.mean():>9.3f}{sd_.mean():>9.3f}{sd_.mean()-b.mean():>8.3f}"
              + (f"{sl.mean():>9.3f}" if len(sl) == len(b) else f"{'--':>9}")
              + f"{sr.mean():>9.3f}   {mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}")

    print("\nC. WHERE THE ANGSTROMS ARE  (means over targets)")
    fullb = g("full", "best"); k500 = g("500", "best"); k75 = g("75", "best")
    ca500 = g("500", "coord_avg")
    print(f"  incumbent (deployed)                                3.204")
    print(f"  raw coordinate average over the SHIPPED TOP-75      3.048   <- the right comparison")
    print(f"    (+ the 0.155 A repair tax = 3.203 ~ the incumbent: the incumbent IS that")
    print(f"     average plus repair.  A2, AUDIT.)")
    print(f"  coordinate average over ALL of K = 500              {ca500.mean():>6.3f}   "
          f"<- a DIFFERENT operator on a DIFFERENT set; not comparable to 3.204")
    print(f"  ORACLE best in top-75      [truncation ceiling]    {k75.mean():>6.3f}")
    print(f"  ORACLE best in K = 500     [truncation ceiling]    {k500.mean():>6.3f}")
    print(f"  ORACLE best in FULL pool   [retrieval ceiling]     {fullb.mean():>6.3f}")
    print(f"\n  the K=500 truncation costs   {k500.mean()-fullb.mean():>6.3f} A of ceiling")
    print(f"  the top-75 truncation costs  {k75.mean()-fullb.mean():>6.3f} A of ceiling")
    print(f"  targets whose FULL-pool best is still above 2.0 A: "
          f"{int((fullb>2.0).sum())}/{len(fullb)}  (K=500: {int((k500>2.0).sum())})")
    print(f"  targets whose FULL-pool best is still above 2.5 A: "
          f"{int((fullb>2.5).sum())}/{len(fullb)}  (K=500: {int((k500>2.5).sum())})")

    print("\nD. THE CENTRAL QUESTION OF SECTION 9 -- does widening K help the SELECTOR,")
    print("   or only the ceiling?")
    b75, s75 = g("75", "best"), g("75", "sel_dist")
    for k in ("150", "300", "500", "1000", "2000", "5000", "full"):
        b = g(k, "best"); s = g(k, "sel_dist")
        mu, lo, hi = _boot(s - s75, rng)
        print(f"  K={k:>5}  ceiling {b.mean()-b75.mean():+.3f}   realized {mu:+.3f} "
              f"[{lo:+.3f},{hi:+.3f}]   recovered "
              f"{(0.0 if abs(b.mean()-b75.mean())<1e-9 else mu/(b.mean()-b75.mean())):>6.1%}")
    print("\n  'recovered' is the fraction of the extra CEILING that the distance selector")
    print("  actually converts into realized accuracy.  Near zero means widening K buys a")
    print("  bigger ceiling and nothing else -- i.e. ranking is the blocker, not recall.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
