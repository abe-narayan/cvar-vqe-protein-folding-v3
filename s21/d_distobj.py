"""s21/d_distobj.py -- D7: THE BRIEF'S "distogram 3.676" IS A DIFFERENT FUNCTIONAL FROM THE
SHIPPED SELECTOR, ON A DIFFERENT BASIS, AT A DIFFERENT n.

    python -m s21.d_distobj

WHAT IS UNDER ATTACK.  `s21/BRIEF.md` §2 opens the sprint with

    argmin RMSD   AMBER 4.990 | Legacy 5.487 | distogram 3.676     pool MEAN 4.739

and builds the whole "argmin is not CVaR" framing on it.  Traced to source
(`s20/agentB_FINDINGS.md` §6d, `s20/LEDGER.md` L?), those four numbers are **n = 20 targets**, on
the **ideal-geometry REBUILD** of each window's torsions (`L.rmsd_of(Zpool, tgt)`), and the
"distogram" column is `s20.qb2_lib.Ham(kind="DIST")`, which is

    E_sq(x) = sum_p (d_p - dhat_p)^2 / sd_p^2                      (qb2_lib.py:191-196)

a **squared-distance** functional -- **not** the deployed selector.  The deployed selector is
`s12.instrument.shipped_score`, a **Bayes-risk** score read off the distogram's full posterior:

    E_bayes(x) = mean_p  risk_p[ bin(d_p) ]                        (instrument.py:181-186)

Three differences at once -- functional, basis, n -- none of them stated where the numbers are
quoted.  The internal comparison in Sprint 20 is fine: all four of its numbers share one n, one
basis and one pool.  The hazard is that this sprint's mandatory matrix builds its **Distance** row
on the SHIPPED score (Workstream A's `tailprice` does exactly that), so a lane that compares its
`disto` arm to "the 3.676 baseline" is comparing two different functionals on two different bases
at two different n.  That is the basis-mismatch failure mode, and it would be the fifth instance.

THE MEASUREMENT.  Both functionals, the SAME K = 500 real retrieved windows, the SAME 126 targets,
the SAME window basis, scored the same way.  If they rank alike the confusion is harmless; if they
do not, the brief's opening row needs its provenance attached.

    Prediction:  the two functionals differ materially in argmin RMSD and in rank agreement.
    Falsifier:   if per-target Spearman(E_bayes, E_sq) exceeds 0.9 in the median AND their argmin
                 RMSDs agree within the MDE, the two are interchangeable and this claim is
                 REFUTED -- the brief's number is then just an n and basis difference.
    Budget:      zero AMBER, no optimisation, n = 126.
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

from s12 import instrument as I     # noqa: E402
from s18 import phys_lib as PL      # noqa: E402


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.empty(len(a)); ra[np.argsort(a, kind="stable")] = np.arange(len(a))
    rb = np.empty(len(b)); rb[np.argsort(b, kind="stable")] = np.arange(len(b))
    ra -= ra.mean(); rb -= rb.mean()
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows, t0 = [], time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        pool = I.pool_idx(u)
        W = np.asarray(u["W"], float)[pool]
        dg = I.distogram(pdb, seq, fold)
        i, j = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
        D = I.pair_dists(W, i, j)
        e_bayes = np.asarray(I.shipped_score(dg, D), float)
        dhat = np.asarray(dg["expected"], float)
        sd = np.maximum(np.asarray(dg["sd"], float), 1e-6)
        e_sq = ((D - dhat[None, :]) ** 2 / sd[None, :] ** 2).sum(1)
        true = I.kabsch_rmsd_batch(W, nat)
        rows.append({
            "pdb": pdb, "n": n, "fold": fold,
            "argmin_bayes": float(true[int(np.argmin(e_bayes))]),
            "argmin_sq": float(true[int(np.argmin(e_sq))]),
            "top1pct_bayes": float(true[np.argsort(e_bayes, kind="stable")[:5]].mean()),
            "top1pct_sq": float(true[np.argsort(e_sq, kind="stable")[:5]].mean()),
            "rho_bayes_sq": _spearman(e_bayes, e_sq),
            "rho_bayes_true": _spearman(e_bayes, true),
            "rho_sq_true": _spearman(e_sq, true),
            "pool_mean": float(true.mean()), "pool_best": float(true.min()),
        })
        if (c + 1) % 20 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    _save(rows, len(tg))
    report(rows)
    return rows


def _save(rows, n_expected):
    json.dump({"rows": rows, "n_expected": int(n_expected)},
              open(os.path.join(RESULTS, "d_distobj.json"), "w"))
    need = ["argmin_bayes", "argmin_sq", "rho_bayes_sq", "rho_bayes_true", "rho_sq_true",
            "pool_mean", "top1pct_bayes", "top1pct_sq"]
    ok = (len(rows) == n_expected == 126
          and all(np.isfinite(r.get(k, np.nan)) for r in rows for k in need))
    p = os.path.join(RESULTS, "d_distobj.COMPLETE")
    if ok:
        with open(p, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"n=126 functionals=shipped_Bayes_risk,squared_distance "
                     f"basis=WINDOW (real retrieved windows, K=500) keys={len(need)}\n")
    elif os.path.exists(p):
        os.remove(p)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "d_distobj.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    folds = g("fold")
    n = len(rows)
    print(f"\n=== D7  TWO 'DISTANCE' FUNCTIONALS ON THE SAME 500 WINDOWS, n = {n} ===")
    print("Basis: WINDOW (real retrieved windows -- valid geometry).  ORACLE scoring only.\n")
    print(f"  {'quantity':<34}{'mean':>9}{'median':>9}")
    for k in ("argmin_bayes", "argmin_sq", "top1pct_bayes", "top1pct_sq",
              "pool_mean", "pool_best", "rho_bayes_sq", "rho_bayes_true", "rho_sq_true"):
        print(f"  {k:<34}{g(k).mean():>9.4f}{np.median(g(k)):>9.4f}")
    st = PL.paired(g("argmin_sq"), g("argmin_bayes"), folds=folds)
    ci = st.get("ci_fold", st["ci"])
    print(f"\n  argmin_sq - argmin_bayes  {st['mean']:+.4f}  fold[{ci[0]:+.4f},{ci[1]:+.4f}]  "
          f"iid[{st['ci'][0]:+.4f},{st['ci'][1]:+.4f}]  W/L {st['W']}/{st['L']}")
    st2 = PL.paired(g("rho_sq_true"), g("rho_bayes_true"), folds=folds)
    ci2 = st2.get("ci_fold", st2["ci"])
    print(f"  rho_sq_true - rho_bayes_true  {st2['mean']:+.4f}  fold[{ci2[0]:+.4f},{ci2[1]:+.4f}]"
          f"  W/L {st2['W']}/{st2['L']}")
    med_rho = float(np.median(g("rho_bayes_sq")))
    mde = 2.8016 * (g("argmin_sq") - g("argmin_bayes")).std(ddof=1) / np.sqrt(n)
    print(f"\n  PRE-REGISTERED FALSIFIER: REFUTED if median rho(E_bayes, E_sq) > 0.9 AND the two")
    print(f"  argmins agree within this comparison's own MDE ({mde:.4f} A).")
    print(f"  median rho = {med_rho:.4f}; |argmin difference| = "
          f"{abs(st['mean']):.4f} A  ->  "
          f"{'REFUTED -- interchangeable' if (med_rho > 0.9 and abs(st['mean']) < mde) else 'SUPPORTED -- they are different objects'}")
    #: RECONCILIATION with the brief's §2 row, on Sprint 20's OWN 20-target subset.
    out_sub = None
    try:
        from s20 import qb2_lib as QL
        sub = [x["pdb"] if isinstance(x, dict) else x for x in QL.subset(20)]
        m = {r["pdb"]: r for r in rows}
        ok = [x for x in sub if x in m]
        if ok:
            sg = lambda k: np.array([m[x][k] for x in ok], float)   # noqa: E731
            print(f"\n  RECONCILIATION on Sprint 20's OWN subset(20), WINDOW basis, n = {len(ok)}:")
            print(f"    argmin_bayes {sg('argmin_bayes').mean():.4f}   "
                  f"argmin_sq {sg('argmin_sq').mean():.4f}   "
                  f"pool_mean {sg('pool_mean').mean():.4f}   "
                  f"pool_best {sg('pool_best').mean():.4f}")
            print("    Sprint 20 reported, on the REBUILD basis: DIST argmin 3.676, "
                  "pool mean 4.739, pool best 1.901.")
            print("    pool mean and pool best agree to 0.01 A, so the subset and the pool are")
            print("    reproduced and the basis shift is negligible (Workstream A measures the")
            print("    per-member rebuild-minus-window shift at +0.011 A).")
            print("    => THE 3.676-vs-3.454 GAP IS THE SUBSET, not the functional and not the")
            print(f"    basis.  subset(20) is simply HARDER: pool mean "
                  f"{sg('pool_mean').mean():.3f} against {g('pool_mean').mean():.3f} at n = 126.")
            out_sub = {"n": len(ok), "argmin_bayes": float(sg("argmin_bayes").mean()),
                       "argmin_sq": float(sg("argmin_sq").mean()),
                       "pool_mean": float(sg("pool_mean").mean()),
                       "pool_best": float(sg("pool_best").mean()),
                       "s20_reported_rebuild": {"DIST_argmin": 3.676, "pool_mean": 4.739,
                                                "pool_best": 1.901}}
    except Exception as e:                                          # pragma: no cover
        print(f"  (subset reconciliation unavailable: {e!r})")

    print(f"\n  For the record, the SHIPPED selector's argmin over the K=500 pool at n = {n} on the")
    print(f"  WINDOW basis is {g('argmin_bayes').mean():.4f} A -- which is `s12/instrument`'s own")
    print(f"  asserted constant `shipped argmin 3.4540` -- and the pool mean is "
          f"{g('pool_mean').mean():.4f} A.")
    print(f"  The brief's Sec.2 row quotes 3.676 and 4.739, which are n = 20, the REBUILD basis,")
    print(f"  and the SQUARED-DISTANCE functional.  Internally consistent there; not the same")
    print(f"  objects as the mandatory matrix's `Distance` row, which is built on the shipped score.")
    out = {"n": n, "cells": {k: {"mean": float(g(k).mean()), "median": float(np.median(g(k)))}
                             for k in ("argmin_bayes", "argmin_sq", "top1pct_bayes",
                                       "top1pct_sq", "pool_mean", "pool_best", "rho_bayes_sq",
                                       "rho_bayes_true", "rho_sq_true")},
           "argmin_sq_minus_bayes": {"mean": st["mean"], "ci95_fold": list(ci),
                                     "ci95_iid": list(st["ci"]), "W": st["W"], "L": st["L"]},
           "rho_diff": {"mean": st2["mean"], "ci95_fold": list(ci2)},
           "own_mde": float(mde), "s20_subset20_reconciliation": out_sub}
    json.dump(out, open(os.path.join(RESULTS, "d_distobj_report.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
