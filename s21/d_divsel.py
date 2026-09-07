"""s21/d_divsel.py -- R1 PILOT: SELECT THE SET AS A SET, NOT AS A TOP-K OF INDEPENDENT SCORES.

    python -m s21.d_divsel run | report

WHY (the mechanism comes out of D2, not out of the air).  `s21/d_cvarop.py` measured, at n = 126
and against matched-count random controls, that a Hamiltonian's sign FLIPS with the readout:
Legacy BEATS its control on a single-member readout (-0.407 [-0.557, -0.259]) and LOSES to it on
the medoid and the coordinate average (+0.302, +0.329, CIs excluding zero).  The mechanism is
s20 L2c -- Legacy prefers compact, pool-typical geometry, so its tail is a tightly clustered,
LOW-DIVERSITY set, and that destroys the error diversity the averaging operator lives on.

The deployed terminal operator is a coordinate average over a selected set, and project memory
`averaging-space-beats-the-objective` prices WHERE you average at 1.024 A against WHAT you rank
with at 0.171 A.  Yet **every selector in this project ranks candidates independently and takes a
top-k.**  Nothing has ever selected the set AS A SET.

    THE QUESTION: within a native-free score band, does choosing the 75 members to maximise
    geometric DIVERSITY beat choosing the 75 best-scoring ones, through the same averaging
    readout?

PRE-REGISTRATION (`s21/PREREG_D.md` R1; this docstring fixes the arms before the run).

  Hypothesis:  the averaging operator's value is set-diversity-driven, so a diverse in-band set
               averages better than the top-k of the same band.
  Primary:     divmax(band) - top75, and divmax(band) - rand_in_band(band), paired over 126
               targets, bootstrap CI over targets, fold-clustered interval beside it.
  Falsifier:   R1 is DEAD if divmax does not beat BOTH the incumbent top75 AND the
               rand_in_band control by more than the MDE (0.084 A) with a CI excluding zero.
               It is ALSO dead if `divmin`, the sign control, does not lose -- a criterion that
               does not hurt when reversed is not measuring what it says.
  Nulls:       rand_in_band is the OPERATIVE control (a diverse set and a random set are easy to
               confuse); divmin is the SIGN control; top75 is the incumbent.
  Budget:      n = 126, pool-restricted, no optimisation, no energy evaluation beyond the
               shipped distogram score.  Reuses the one 500x500 pairwise-RMSD matrix per target.
  Promotion:   none from a pilot.  A surviving effect earns a dev-set replication, not a claim.

  PRIOR AGAINST, STATED FIRST.  `decorrelated-errors-exist-but-are-unusable` shows fusion gains
  going as the SQUARE of the weaker channel's skill.  That is about fusing two PREDICTORS and
  this is about choosing an averaging SET -- a different operator -- but it is the reason to
  expect a small effect and to kill this quickly if the first CI spans zero.

BASIS.  Every RMSD here is a **POINT CLOUD** (a coordinate average is not a buildable backbone,
s20 X1).  The incumbent on this basis is 3.0483 A, which `d_cvarop.py` reproduced BIT-EXACTLY as
`tail_avg(0.15) | disto`, so the comparison is anchored to the production number and not to a
re-implementation of it.  Nothing here is comparable to the 3.204 built-chain incumbent.

NO NATIVE INFORMATION enters any selection.  The native is read only to score, and every ORACLE
row is labelled.
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

from s12 import instrument as I      # noqa: E402
from s15 import seed as SD           # noqa: E402
from s18 import phys_lib as PL       # noqa: E402

M = I.M                 # 75, the shipped set size
BANDS = (0.30, 0.50, 1.00)          # fraction of the K=500 pool the set may be drawn from
R_DRAW = 8                          # repeats for the random-in-band control
CFG = {"M": M, "BANDS": list(BANDS), "R_DRAW": R_DRAW, "K": I.K,
       "score": "shipped Bayes-risk distogram", "basis": "POINT CLOUD"}


def _greedy(P, cand, k, seed_idx, sign=+1.0):
    """Greedy max-min (sign=+1) or min-max (sign=-1) selection of k members from `cand`.

    max-min: repeatedly add the candidate whose MINIMUM pairwise RMSD to the already-chosen set
    is largest -- the standard farthest-point / max-min-diversity construction.  sign=-1 gives
    the mirror-image criterion and is the SIGN CONTROL: if diversity is what matters, reversing
    it must hurt.
    """
    cand = np.asarray(cand, int)
    chosen = [int(seed_idx)]
    rest = [int(c) for c in cand if int(c) != int(seed_idx)]
    d = P[np.ix_(rest, chosen)].min(1)
    while len(chosen) < k and rest:
        j = int(np.argmax(sign * d))
        chosen.append(rest[j])
        rest.pop(j)
        d = np.delete(d, j)
        if rest:
            d = np.minimum(d, P[np.ix_(rest, [chosen[-1]])].ravel())
    return np.asarray(chosen, int)


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows, t0 = [], time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        pool = I.pool_idx(u)
        W = np.asarray(u["W"], float)[pool]
        m = len(W)
        rng = SD.stable_rng(pdb, "s21_d_divsel")
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        s = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        order = np.argsort(s, kind="stable")
        P = I.pairwise_rmsd(W)
        true = I.kabsch_rmsd_batch(W, nat)

        def avg_rmsd(idx):
            idx = np.asarray(idx, int)
            a, _b = I.coordinate_average(W[idx], P[np.ix_(idx, idx)])
            return float(I.ca_rmsd(np.asarray(a, float), nat))

        e = {"pdb": pdb, "n": n, "fold": fold, "m": m,
             "incumbent_top75": avg_rmsd(order[:M]),
             "pool_mean_ORACLE": float(true.mean()),
             "ORACLE_pool_best": float(true.min())}
        for B in BANDS:
            nb = max(M, int(round(B * m)))
            band = order[:nb]
            e[f"divmax{B:g}"] = avg_rmsd(_greedy(P, band, M, band[0], +1.0))
            e[f"divmin{B:g}"] = avg_rmsd(_greedy(P, band, M, band[0], -1.0))
            e[f"randband{B:g}"] = float(np.mean(
                [avg_rmsd(rng.choice(band, M, replace=False)) for _ in range(R_DRAW)]))
            #: diagnostics, native-free: how diverse each set actually is
            dm = _greedy(P, band, M, band[0], +1.0)
            e[f"spread_divmax{B:g}"] = float(P[np.ix_(dm, dm)].mean())
            e[f"spread_top75"] = float(P[np.ix_(order[:M], order[:M])].mean())
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            _save(rows, len(tg))
    _save(rows, len(tg))
    return rows


def _save(rows, n_expected):
    json.dump({"rows": rows, "cfg": CFG, "n_expected": int(n_expected)},
              open(os.path.join(RESULTS, "d_divsel.json"), "w"))
    need = (["incumbent_top75", "spread_top75"]
            + [f"{a}{B:g}" for B in BANDS
               for a in ("divmax", "divmin", "randband", "spread_divmax")])
    ok = (len(rows) == n_expected == 126
          and all(np.isfinite(r.get(k, np.nan)) for r in rows for k in need))
    p = os.path.join(RESULTS, "d_divsel.COMPLETE")
    if ok:
        with open(p, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"n=126 M={M} bands={list(BANDS)} R_DRAW={R_DRAW} "
                     f"arms=divmax,divmin,randband,incumbent_top75 keys={len(need)}\n")
    elif os.path.exists(p):
        os.remove(p)


def report():
    d = json.load(open(os.path.join(RESULTS, "d_divsel.json")))
    rows = d["rows"]
    g = lambda k: np.array([r[k] for r in rows])       # noqa: E731
    folds = g("fold")
    print(f"\n=== R1 PILOT  set-diversity selection, n = {len(rows)}, "
          f"basis POINT CLOUD ===")
    print(f"  incumbent top-75 by score : {g('incumbent_top75').mean():.4f} A "
          f"(production `rmsd_avg` = 3.0483, reproduced bit-exactly by d_cvarop)")
    print(f"  pool mean (ORACLE)        : {g('pool_mean_ORACLE').mean():.4f}")
    print(f"  mean pairwise spread of the top-75 set: {g('spread_top75').mean():.3f} A\n")

    def PP(a, b):
        st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
        return st, st.get("ci_fold", st["ci"])

    out = {"n": len(rows), "cfg": d["cfg"], "cells": {}}
    for B in BANDS:
        nb = int(round(B * I.K))
        print(f"  band = top {B:g} of the pool ({max(M, nb)} candidates, choose {M})   "
              f"divmax spread {g(f'spread_divmax{B:g}').mean():.3f} A")
        for arm in ("divmax", "divmin", "randband"):
            print(f"    {arm:<10}{g(f'{arm}{B:g}').mean():>8.4f}", end="")
            for ref, refname in ((g("incumbent_top75"), "vs top75"),
                                 (g(f"randband{B:g}"), "vs rand_in_band")):
                st, ci = PP(g(f"{arm}{B:g}"), ref)
                out["cells"][f"{arm}|{B:g}|{refname}"] = {
                    "mean": st["mean"], "ci95_fold": list(ci), "ci95_iid": list(st["ci"]),
                    "W": st["W"], "L": st["L"], "median": st["median"]}
                print(f"   {refname} {st['mean']:+7.4f} [{ci[0]:+.4f},{ci[1]:+.4f}] "
                      f"{st['W']}/{st['L']}", end="")
            print()
        print()
    print("PRE-REGISTERED FALSIFIER: R1 is DEAD unless divmax beats BOTH top75 AND")
    print("rand_in_band by more than the MDE (0.084 A) with a CI excluding zero, AND divmin")
    print("(the sign control) loses.  A residual inside the MDE is NOT MEASURED, not 'matched'.")
    best = min(((abs(out['cells'][f'divmax|{B:g}|vs top75']['mean']), B) for B in BANDS))
    for B in BANDS:
        a = out["cells"][f"divmax|{B:g}|vs top75"]
        b = out["cells"][f"divmax|{B:g}|vs rand_in_band"]
        s = out["cells"][f"divmin|{B:g}|vs top75"]
        alive = (a["mean"] < -0.084 and a["ci95_fold"][1] < 0
                 and b["mean"] < -0.084 and b["ci95_fold"][1] < 0
                 and s["mean"] > 0)
        print(f"  band {B:g}: {'ALIVE' if alive else 'DEAD'}")
    json.dump(out, open(os.path.join(RESULTS, "d_divsel_report.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
