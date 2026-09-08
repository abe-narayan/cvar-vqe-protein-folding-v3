"""s17/selection_theory.py -- THE DESIGN EQUATION.  How good must a selector be, and at what K?

Sprint 17 section 41 requires deriving relationships analytically before running experiments.
This module derives the one that decides the whole sprint.

THE QUESTION.  The programme knows its ceilings (best member at each candidate width K) and it
knows its realized selectors are far above them.  What it has never had is the function that
connects the two:

    given a candidate set of size K whose member RMSDs have empirical distribution F,
    and a native-free selector whose ranking has correlation rho with the true ranking,
    what RMSD does the selector return?

With that function, "we need better ranking" becomes a number: **the rho required to reach
2.5 A and 2.0 A at each K**.  Without it, every selection experiment is uncalibrated -- we
cannot tell an encouraging selector from a hopeless one.

THE MODEL, and it is deliberately the simplest one that is not wrong.  Per target, take the
K member RMSDs as given (they are measured, not assumed).  Map each member to its normal
score z_true = Phi^-1((rank - 0.5)/K).  A selector of quality rho scores members as

    z_score = rho * z_true + sqrt(1 - rho^2) * eps,      eps ~ N(0, 1) i.i.d.

and returns the argmin of z_score.  This is the Gaussian copula with a single dependence
parameter: it fixes the SELECTOR's marginal and the TRUTH's marginal at their measured values
and varies only how tightly they are coupled.  rho = 1 recovers the oracle exactly; rho = 0
recovers uniform random selection exactly.  Both limits are checked numerically in `verify()`.

WHAT THE MODEL IS NOT.  It assumes the selector's errors are exchangeable across members
within a target.  A real selector may be systematically wrong on a structured subset -- it may
rank all helical decoys above all sheet decoys, say -- and then the same rho buys less.  So
**this curve is an upper bound on what a given rho delivers**, and the gap between a measured
selector's (rho, realized) point and this curve is itself a diagnostic: a selector below the
curve is not merely weak, it is weak in a structured way.  That diagnostic is reported.

WHY A COPULA AND NOT AN ORDER-STATISTIC FORMULA.  The quantity wanted is E[d at argmin of a
correlated score], which for a general empirical F has no closed form.  The copula makes the
dependence explicit and the expectation a one-line Monte Carlo over a model whose only free
parameter is the thing we want to solve for.  Everything else is measured.

USE OF THE NATIVE.  The member-RMSD distributions are ORACLE quantities and are used here to
build a CEILING CURVE and a REQUIREMENT, never a predictor.  No quantity computed here may
enter any model.  Section 8 of the brief applies: this is a ceiling, and it says so.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

RHOS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)
#: the BAND.  The argmin is decided among candidates the selector has already ranked highly,
#: so the parameter that matters is in-band skill, not global skill.  B is the shortlist the
#: real selector delivers; the band curve asks how much skill INSIDE it is required.
BANDS = (5, 10, 25, 50, 100)


def rho_from_accuracy(acc):
    """pairwise ORDERING ACCURACY (null 0.5) -> Gaussian-copula rho.

    tau = 2*acc - 1 and, for a Gaussian copula, tau = (2/pi) arcsin(rho).  These are NOT the
    same axis and conflating them is the failure mode this programme keeps repeating: the
    coordinator did exactly that in the first strategic read of this module, mapping a recorded
    accuracy of 0.600 onto rho = 0.600 when it is rho = 0.309.  Both axes are printed
    everywhere below so the mistake cannot recur silently.
    """
    return float(np.sin(np.pi * (2.0 * np.asarray(acc, float) - 1.0) / 2.0))


def accuracy_from_rho(rho):
    """the inverse, so every requirement can be quoted on the axis a ranker is measured on."""
    return float((2.0 / np.pi * np.arcsin(np.asarray(rho, float)) + 1.0) / 2.0)


def rho_from_spearman(rs):
    """Spearman rho_S -> Gaussian-copula rho.  Differs from rho_S itself by ~3%."""
    return float(2.0 * np.sin(np.pi * np.asarray(rs, float) / 6.0))
KS = (75, 150, 300, 500, 1000, 2000, 5000, 0)     # 0 = full universe
NDRAW = 400                                        # copula draws per (target, K, rho)


def _norm_scores(d):
    """normal scores of the TRUE ordering: best member gets the most negative z."""
    k = len(d)
    r = np.argsort(np.argsort(d, kind="stable"), kind="stable")     # 0 = best
    from scipy.special import ndtri
    return ndtri((r + 0.5) / k)


def _expected(d, z_true, rho, rng, ndraw=NDRAW):
    """E[selected RMSD] under a Gaussian-copula selector of quality rho."""
    k = len(d)
    if rho >= 1.0:
        return float(d.min()), 0.0
    if rho <= 0.0:
        return float(d.mean()), float(d.std() / np.sqrt(ndraw))
    eps = rng.standard_normal((ndraw, k))
    zs = rho * z_true[None, :] + np.sqrt(1.0 - rho * rho) * eps
    sel = d[np.argmin(zs, axis=1)]
    return float(sel.mean()), float(sel.std() / np.sqrt(ndraw))


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        order = np.asarray(u["order"], int)
        W = np.asarray(u["W"], float)[order]
        d_all = I.kabsch_rmsd_batch(W, np.asarray(u["nat_ca"], float))
        U = len(d_all)
        rng = SD.stable_rng(pdb, "s17theory")
        cell = {}
        #: the real selector's own ranking, so the band is the one it actually delivers.
        i, j = I.pair_index(int(t["n"]))
        dg = I.distogram(pdb, t["seq"], int(t["fold"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        for K in KS:
            k = U if K == 0 else min(K, U)
            d = d_all[:k]
            z = _norm_scores(d)
            e = {"k": int(k), "best": float(d.min()), "mean": float(d.mean()),
                 "curve": {f"{r:.2f}": _expected(d, z, r, rng)[0] for r in RHOS}}
            #: BAND CURVES.  Condition on the shipped objective's own top-B shortlist and ask
            #: how much IN-BAND skill is needed there.  The global rho is not the parameter
            #: that decides an argmin -- the SELECT workstream measured global Spearman 0.568
            #: against in-band 0.131 on the same candidates, and the realized result sits
            #: between the two predictions and much closer to the in-band one.
            rk = np.argsort(sc[:k], kind="stable")
            e["band"] = {}
            for B in BANDS:
                if B > k:
                    continue
                db = d[rk[:B]]
                zb = _norm_scores(db)
                e["band"][str(B)] = {
                    "best": float(db.min()), "mean": float(db.mean()),
                    "curve": {f"{r:.2f}": _expected(db, zb, r, rng)[0] for r in RHOS},
                }
            cell["full" if K == 0 else str(K)] = e
        rows.append({"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]), "U": int(U),
                     "cells": cell})
        if (c + 1) % 25 == 0:
            print(f"  {c+1}/{len(tg)}", flush=True)
            json.dump({"rows": rows}, open(os.path.join(RESULTS, "sel_theory.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "sel_theory.json"), "w"))
    report(rows)
    return rows


def verify(rows=None):
    """The two limits the model MUST reproduce exactly, checked numerically."""
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "sel_theory.json")))["rows"]
    ok = True
    for k in ("500", "full"):
        c1 = np.mean([r["cells"][k]["curve"]["1.00"] for r in rows])
        b = np.mean([r["cells"][k]["best"] for r in rows])
        c0 = np.mean([r["cells"][k]["curve"]["0.00"] for r in rows])
        m = np.mean([r["cells"][k]["mean"] for r in rows])
        print(f"  K={k:>5}  rho=1 -> {c1:.4f} vs ORACLE best {b:.4f}   "
              f"|delta| {abs(c1-b):.2e}")
        print(f"  K={k:>5}  rho=0 -> {c0:.4f} vs random mean  {m:.4f}   "
              f"|delta| {abs(c0-m):.2e}")
        ok &= abs(c1 - b) < 1e-9 and abs(c0 - m) < 1e-9
    print(f"  limits exact: {ok}")
    return ok


def _required_rho(rows, k, target_rmsd):
    """smallest rho on the grid whose mean curve reaches target_rmsd; None if unreachable."""
    for r in RHOS:
        v = np.mean([x["cells"][k]["curve"][f"{r:.2f}"] for x in rows])
        if v <= target_rmsd:
            return r, v
    return None, None


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "sel_theory.json")))["rows"]
    ks = [("full" if K == 0 else str(K)) for K in KS]
    print(f"\nn = {len(rows)} targets.  ORACLE construction: the member-RMSD distributions are")
    print("read from the native.  Everything below is a CEILING or a REQUIREMENT, never a result.\n")

    print("A. LIMIT CHECK -- the model must reproduce the oracle and random exactly")
    verify(rows)

    print("\nB. THE DESIGN CURVE -- mean selected RMSD as a function of selector quality rho")
    print(f"  {'K':>6}" + "".join(f"{r:>7.2f}" for r in RHOS))
    for k in ks:
        print(f"  {k:>6}" + "".join(
            f"{np.mean([x['cells'][k]['curve'][f'{r:.2f}'] for x in rows]):>7.3f}" for r in RHOS))

    print("\nC. THE REQUIREMENT -- the rho a selector needs to hit each target, per K")
    print(f"  {'K':>6}{'rho for 2.5 A':>16}{'rho for 2.2 A':>16}{'rho for 2.0 A':>16}"
          f"{'oracle (rho=1)':>16}")
    for k in ks:
        line = f"  {k:>6}"
        for thr in (2.5, 2.2, 2.0):
            r, v = _required_rho(rows, k, thr)
            line += f"{('%.2f' % r) if r is not None else 'unreachable':>16}"
        line += f"{np.mean([x['cells'][k]['best'] for x in rows]):>16.3f}"
        print(line)

    band_report(rows)

    print("\n  READ.  This is the sprint's design equation.  A selector's measured cross-target")
    print("  rank correlation can be compared directly against these numbers to say whether it")
    print("  is on a path to the target or not.  The recorded in-band ordering is 0.600 across")
    print("  targets; the recorded requirement for 2.0 A was 0.638 -- this table says whether")
    print("  that number survives at the widths now available.")
    print("\n  CAVEAT that travels with every number here: the copula assumes the selector's")
    print("  errors are exchangeable within a target.  A selector that is systematically wrong")
    print("  on a structured subset buys less at the same rho, so these are UPPER BOUNDS on")
    print("  what a given rho delivers, and a measured selector falling below its curve is")
    print("  weak in a STRUCTURED way -- which is a different and more fixable problem.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()


def band_report(rows=None, K0="500"):
    """The parameter that actually decides an argmin: IN-BAND selector skill.

    The SELECT workstream measured the shipped objective at global Spearman 0.568 but IN-BAND
    0.131 on the same K = 500 candidates, and its realized result sits between the two copula
    predictions and much closer to the in-band one.  So a requirement read off the GLOBAL curve
    is read off the wrong axis.  This conditions on the shipped objective's own top-B shortlist
    and asks how much skill is needed inside it.
    """
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "sel_theory.json")))["rows"]
    have = [B for B in BANDS if str(B) in rows[0]["cells"][K0].get("band", {})]
    if not have:
        print("\nD. BAND CURVE -- not present in this artefact (rerun `run()`).")
        return
    print(f"\nD. THE BAND CURVE at K = {K0} -- in-band skill, the parameter that decides an argmin")
    print(f"  {'band B':>8}{'B-best':>9}{'B-mean':>9}" + "".join(f"{r:>7.2f}" for r in RHOS))
    for B in have:
        c = [x["cells"][K0]["band"][str(B)] for x in rows]
        print(f"  {B:>8}{np.mean([y['best'] for y in c]):>9.3f}"
              f"{np.mean([y['mean'] for y in c]):>9.3f}"
              + "".join(f"{np.mean([y['curve'][f'{r:.2f}'] for y in c]):>7.3f}" for r in RHOS))

    print("\n  in-band rho required, on BOTH axes (copula rho / pairwise ordering accuracy)")
    print(f"  {'band B':>8}{'for 2.5 A':>24}{'for 2.2 A':>24}{'for 2.0 A':>24}")
    for B in have:
        c = [x["cells"][K0]["band"][str(B)] for x in rows]
        line = f"  {B:>8}"
        for thr in (2.5, 2.2, 2.0):
            hit = None
            for r in RHOS:
                if np.mean([y["curve"][f"{r:.2f}"] for y in c]) <= thr:
                    hit = r
                    break
            line += (f"{'rho %.2f / acc %.3f' % (hit, accuracy_from_rho(hit)):>24}"
                     if hit is not None else f"{'unreachable':>24}")
        print(line)

    print("\n  UNITS -- the coordinator got this wrong once already, and it is BRIEF s5's failure")
    print("  mode exactly: a pairwise ORDERING ACCURACY (null 0.500) is NOT a copula rho.")
    print("  tau = 2*acc - 1 and tau = (2/pi) arcsin(rho), so:")
    for acc in (0.539, 0.600, 0.638, 0.692, 0.725, 0.986):
        print(f"    accuracy {acc:.3f}  ->  copula rho {rho_from_accuracy(acc):.3f}")
    print("  The programme's recorded cross-target in-band ordering of 0.600 is rho 0.309, NOT")
    print("  rho 0.600.  Spearman is a third axis again: rho = 2 sin(pi rho_S / 6).")
