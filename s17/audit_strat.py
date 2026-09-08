"""s17/audit_strat.py -- STRATIFICATION, CONCENTRATION AND FRESH-INTERPRETER REPRODUCTION
of this audit's own headline (the order-statistic null on the ceiling map).

BRIEF section 4 requires stratification and medians beside means; the programme's standing
memory requires that a near-even W/L with a CI excluding zero be checked for concentration
against a UNIFORM-EFFECT NULL, not against a raw drop-top threshold.  This module applies
both to the audit's own claim, so the audit is held to the standard it is enforcing.

It also recomputes the headline in a fresh interpreter from the artefact and reports
max |delta| against the first run, and re-derives it from the raw npz files without
going through `audit_ordernull.py` at all.
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

from s12 import instrument as I     # noqa: E402
from s15 import seed as SD          # noqa: E402

UNIV = os.path.join(ROOT, "s8", "generate_univ")


def _boot(d, rng, B=8000, folds=None):
    d = np.asarray(d, float)
    if folds is None:
        m = d[rng.integers(0, len(d), size=(B, len(d)))].mean(1)
    else:
        folds = np.asarray(folds); uf = np.unique(folds)
        idx = [np.where(folds == f)[0] for f in uf]
        m = np.array([d[np.concatenate([idx[q] for q in rng.integers(0, len(uf), len(uf))])].mean()
                      for _ in range(B)])
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def independent_headline(R=200):
    """Re-derive observed-vs-null from the npz files, not from audit_ordernull.json."""
    tg = I.targets()
    obs, nul, pdbs, folds, ns = [], [], [], [], []
    for t in tg:
        z = np.load(os.path.join(UNIV, f"{t['pdb']}.npz"), allow_pickle=True)
        rr = np.asarray(z["rr"], float); order = np.asarray(z["order"], int)
        U = len(rr); k = min(500, U)
        obs.append(float(rr[order[:k]].min() - rr.min()))
        rng = SD.stable_rng(t["pdb"], "s17auditstrat")
        v = np.array([rr[rng.permutation(U)[:k]].min() for _ in range(R)])
        nul.append(float(v.mean() - rr.min()))
        pdbs.append(t["pdb"]); folds.append(int(t["fold"])); ns.append(int(t["n"]))
    return (np.array(obs), np.array(nul), pdbs, np.array(folds), np.array(ns))


def recall_ctrl(rows, rng, R=200):
    """The matched-K recall control: BLOSUM top-500 vs a RANDOM 500 of the same universe."""
    print("  MATCHED-K RECALL CONTROL (BLOSUM 500 vs random 500, R=200 draws/target):")
    tg = {t["pdb"]: t for t in I.targets()}
    hits_b = {th: 0.0 for th in (1.5, 2.0, 2.5)}
    hits_r = {th: 0.0 for th in (1.5, 2.0, 2.5)}
    for r in rows:
        z = np.load(os.path.join(UNIV, f"{r['pdb']}.npz"), allow_pickle=True)
        rr = np.asarray(z["rr"], float); order = np.asarray(z["order"], int)
        U = len(rr); k = min(500, U)
        g = SD.stable_rng(r["pdb"], "s17auditrecall")
        mins = np.array([rr[g.permutation(U)[:k]].min() for _ in range(R)])
        bmin = rr[order[:k]].min()
        for th in hits_b:
            hits_b[th] += float(bmin < th)
            hits_r[th] += float((mins < th).mean())
    for th in (1.5, 2.0, 2.5):
        print(f"    t={th} A:  BLOSUM-500 {hits_b[th]:.0f}/{len(rows)}"
              f"   random-500 {hits_r[th]:.1f}/{len(rows)}"
              f"   retrieval buys {hits_b[th]-hits_r[th]:+.1f} targets")


def main():
    rng = SD.stable_rng("s17audit", "strat")
    rows = json.load(open(os.path.join(RESULTS, "audit_ordernull.json")))["rows"]
    pdbs = [r["pdb"] for r in rows]
    fold = np.array([r["fold"] for r in rows])
    n = np.array([r["n"] for r in rows])
    b500 = np.array([r["blosum"]["500"] for r in rows])
    bfull = np.array([r["blosum"]["full"] for r in rows])
    r500 = np.array([r["rand_mean"]["500"] for r in rows])
    rfull = np.array([r["rand_mean"]["full"] for r in rows])
    obs = b500 - bfull
    nul = r500 - rfull
    d = obs - nul

    print("A. FRESH-INTERPRETER, INDEPENDENT RE-DERIVATION OF THE HEADLINE")
    o2, n2, p2, f2, _ = independent_headline()
    assert p2 == pdbs
    print(f"  max |delta| observed gain : {np.abs(o2 - obs).max():.2e} A")
    print(f"  max |delta| null gain     : {np.abs(n2 - nul).max():.3f} A"
          f"   (different permutation seed by design; mean delta "
          f"{np.abs(n2-nul).mean():.3f})")
    m1, l1, h1 = _boot(o2 - n2, rng, folds=f2)
    print(f"  observed - null, independent route: {m1:+.3f} [{l1:+.3f},{h1:+.3f}]")

    print("\nB. THE AUDIT'S OWN CLAIM, HELD TO THE SPRINT'S REPORTING STANDARD")
    m, lo, hi = _boot(d, rng, folds=fold)
    print(f"  mean {m:+.3f} [{lo:+.3f},{hi:+.3f}]   median {np.median(d):+.3f}"
          f"   W/L {int((d<0).sum())}/{int((d>0).sum())}   sd {d.std(ddof=1):.3f}")
    order = np.argsort(d)
    print(f"  drop the 10 most negative targets: {d[order[10:]].mean():+.3f}"
          f"   drop 20: {d[order[20:]].mean():+.3f}")
    # uniform-effect null for the drop-top statistic
    B = 4000
    sim = rng.normal(d.mean(), d.std(ddof=1), size=(B, len(d)))
    sim.sort(1)
    q10 = np.percentile(sim[:, 10:].mean(1), [2.5, 97.5])
    print(f"  UNIFORM-EFFECT NULL for drop-top-10: [{q10[0]:+.3f},{q10[1]:+.3f}]"
          f"  -> observed {d[order[10:]].mean():+.3f} is "
          f"{'INSIDE (no concentration)' if q10[0] <= d[order[10:]].mean() <= q10[1] else 'OUTSIDE'}")
    print("  per fold: " + "  ".join(f"f{f}:{d[fold==f].mean():+.3f}(n={int((fold==f).sum())})"
                                     for f in np.unique(fold)))

    print("\nC. STRATIFICATION")
    for lab, mask in (("n 9-11", n <= 11), ("n 12-13", (n >= 12) & (n <= 13)),
                      ("n 14-16", n >= 14)):
        print(f"  {lab:9} n={int(mask.sum()):3d}  observed gain {obs[mask].mean():.3f}"
              f"   null gain {nul[mask].mean():.3f}   obs-null {d[mask].mean():+.3f}")
    hard = b500 > 2.0
    print(f"  K=500 best > 2.0 A  n={int(hard.sum()):3d}  observed gain {obs[hard].mean():.3f}"
          f"   null {nul[hard].mean():.3f}   obs-null {d[hard].mean():+.3f}")
    print(f"  K=500 best <= 2.0 A n={int((~hard).sum()):3d}  observed gain {obs[~hard].mean():.3f}"
          f"   null {nul[~hard].mean():.3f}   obs-null {d[~hard].mean():+.3f}")
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    print(f"  FAIL18              n={int(f18.sum()):3d}  observed gain {obs[f18].mean():.3f}"
          f"   null {nul[f18].mean():.3f}   obs-null {d[f18].mean():+.3f}")

    print("\nD. WHAT THE WIDER CEILING ACTUALLY BUYS IN RECALL TERMS")
    # ERROR PRESERVED IN PLACE.  The first version of this block printed the recall of the
    # RANDOM-ORDER full universe as a "reference" for the BLOSUM full universe.  That is
    # vacuous: at K = full the set is the whole universe under EVERY permutation, so the
    # two are equal by construction (they printed 77/77, 102/102, 118/118).  It is a
    # tautology, not a control, and it is exactly the shape of mistake this audit exists to
    # catch.  The control that means something is at MATCHED K: random-500 vs BLOSUM-500.
    for th in (1.5, 2.0, 2.5):
        a = int((b500 < th).sum()); bq = int((bfull < th).sum())
        print(f"  targets with a window better than {th} A:  K=500 {a}/126"
              f"  ->  full universe {bq}/126   (+{bq-a})")
    print("  MATCHED-K control for the retrieval ordering (mean over the random draws of")
    print("  the per-target probability that a random 500 contains a window under t):")
    for th in (1.5, 2.0, 2.5):
        rec = np.array([np.mean(np.asarray(r["rand_all_500"], float) < th) for r in rows]) \
            if "rand_all_500" in rows[0] else None
        if rec is None:
            print(f"    t={th}: not stored by audit_ordernull.py; see audit_strat.recall_ctrl()")
            break
    print(f"  targets with the FULL-universe best still > 2.0 A: {int((bfull>2.0).sum())}/126")
    recall_ctrl(rows, rng)

    json.dump({"obs": obs.tolist(), "null": nul.tolist(), "pdb": pdbs,
               "fold": fold.tolist(), "n": n.tolist()},
              open(os.path.join(RESULTS, "audit_strat.json"), "w"))


if __name__ == "__main__":
    main()
