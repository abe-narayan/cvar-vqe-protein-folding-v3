"""FAIL18 forensics, step 15: the DEPLOYABLE analogue of the 10-pair oracle.

E14 (Q1): the distogram's own per-pair sd correlates r = 0.49 with |expected - true| on
FAIL18 and r = 0.48 on the controls, and the ten worst pairs sit at the 71st percentile of
its sd.  So the model has a partial, honest handle on which pairs it got wrong.  Does
using it help?

Arms (deployable; the shipped Bayes-risk score with a per-pair weight):
  shipped      uniform weight                                       (control)
  invsd^p      weight proportional to sd^-p, p in {0.5, 1, 2}
  droptop-k    drop the k pairs with the largest sd, k in {3, 10, 20}
  droplong     drop every pair with |i-j| >= 9 (the shell E2 shows is worst)
Metric: top-75 best and recall (no projection), all 126, paired vs shipped.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I


def main():
    rows = []
    for k, t in enumerate(I.targets()):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]; band = rr <= rr.min() + I.BAND
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n); sep = j - i
        grid = np.asarray(dg["grid"], float); risk = np.asarray(dg["risk"], float)
        sd = np.asarray(dg["sd"], float)
        D = I.pair_dists(Wc, i, j).astype(np.float32).astype(float)
        g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
        R = risk[np.arange(risk.shape[0])[None, :], g]           # (b, npairs)

        def ev(w):
            w = np.asarray(w, float)
            if w.sum() <= 0:
                w = np.ones_like(w)
            sc = (R * w[None, :]).sum(1) / w.sum()
            sub = np.argsort(sc, kind="stable")[:I.M]
            return dict(top75_best=float(rr[sub].min()), top75_mean=float(rr[sub].mean()),
                        recall=int(band[sub].any()))

        arms = {"shipped": ev(np.ones(len(sd)))}
        for p_ in (0.5, 1.0, 2.0):
            arms[f"invsd{p_}"] = ev(1.0 / np.maximum(sd, 1e-3) ** p_)
        for kk in (3, 10, 20):
            w = np.ones(len(sd))
            if kk < len(sd):
                w[np.argsort(-sd)[:kk]] = 0.0
            arms[f"dropsd{kk}"] = ev(w)
        arms["droplong"] = ev((sep < 9).astype(float))
        rows.append(dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18, arms=arms))
    I.write("fail_sdweight", rows)

    f = [r for r in rows if r["fail18"]]; o = [r for r in rows if not r["fail18"]]
    folds = np.array([r["fold"] for r in rows])
    base = np.array([r["arms"]["shipped"]["top75_best"] for r in rows])
    print(f"{'arm':12s}{'all126':>9s}{'d':>8s}{'ci95':>19s}{'FAIL18':>9s}{'other108':>10s}"
          f"{'recF':>6s}{'drop10':>8s}")
    summ = {}
    for a in rows[0]["arms"]:
        A = np.array([r["arms"][a]["top75_best"] for r in rows])
        st = I.paired(A, base, folds=folds)
        summ[a] = dict(all126=float(A.mean()),
                       fail18=float(np.mean([r["arms"][a]["top75_best"] for r in f])),
                       other108=float(np.mean([r["arms"][a]["top75_best"] for r in o])),
                       recall_fail=float(np.mean([r["arms"][a]["recall"] for r in f])),
                       paired=st)
        s = summ[a]
        print(f"{a:12s}{s['all126']:9.3f}{st['mean_diff']:8.3f} "
              f"[{st['ci95'][0]:6.3f},{st['ci95'][1]:6.3f}]{s['fail18']:9.3f}"
              f"{s['other108']:10.3f}{s['recall_fail']:6.2f}"
              f"{(st['drop_top10_mean_diff'] or 0):8.3f}")
    I.write("fail_sdweight_summary", summ)


if __name__ == "__main__":
    main()
