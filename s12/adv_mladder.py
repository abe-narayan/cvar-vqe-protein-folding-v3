"""COORDINATOR RE-PRIORITISATION, T3 -- the DEPLOYABLE m ladder through the REAL projection.

Is a FIXED small m better than the shipped m = 75, under the SHIPPED objective, with no
router and no oracle?  This is the only change in the sprint that would be deployable
today if it holds, so it gets the full production terminal (lambda = 0 `fit` arm, which
A3 showed reproduces production to 7e-4 A -- the lambda arm's ~0.02 A noise floor would
swamp the effect sizes at stake here).

Per target: shipped Bayes-risk score -> top-m -> superpose on the subset medoid ->
coordinate average -> `core.project` -> CA-RMSD to native, for every m in MS.
Then a LEAVE-FOLD-OUT choice of m: pick the m that minimises the mean on the OTHER four
folds and apply it to this fold, so the reported number contains no in-sample tuning.
FAIL18 / other-108 split on every row, per the coordinator's constraint (2).
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I

MS = [1, 3, 5, 10, 25, 50, 75, 150]


def main(lo=0, hi=None, shard=None):
    tg = I.targets()[lo:hi]
    rows = []
    t0 = time.time()
    for t in tg:
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        W = np.asarray(u["W"][p], float); rr = np.asarray(u["rr"][p], float); nat = u["nat_ca"]
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(n)
        sc = I.shipped_score(dg, I.pair_dists(W, i, j))
        order = np.argsort(sc, kind="stable")
        P = I.pairwise_rmsd(W[order[:max(MS)]])
        row = dict(pdb=pdb, n=n, fold=t["fold"], fail=int(pdb in I.FAIL18))
        for m in MS:
            idx = np.arange(m)
            if m == 1:
                C = W[order[0]]
            else:
                C, _ = I.coordinate_average(W[order[:m]], P[np.ix_(idx, idx)])
            row[f"avg_m{m}"] = float(I.ca_rmsd(C, nat))
            row[f"fit_m{m}"] = float(I.ca_rmsd(I.project(C, t["seq"], t["fold"])["fit_ca"], nat))
            row[f"setmean_m{m}"] = float(rr[order[:m]].mean())
            row[f"setbest_m{m}"] = float(rr[order[:m]].min())
        rows.append(row)
        print(f"{len(rows):3d}/{len(tg)} {pdb} " +
              " ".join(f"m{m}={row[f'fit_m{m}']:.2f}" for m in MS) +
              f"  {time.time()-t0:.0f}s", flush=True)

    fail = np.array([r["fail"] for r in rows], bool)
    folds = np.array([r["fold"] for r in rows])
    names = [r["pdb"] for r in rows]
    agg = {"n": len(rows), "MS": MS}
    for tag in ("fit", "avg"):
        agg[tag] = {m: dict(mean=float(np.mean([r[f"{tag}_m{m}"] for r in rows])),
                            fail=float(np.mean([r[f"{tag}_m{m}"] for r in rows if r["fail"]])),
                            other=float(np.mean([r[f"{tag}_m{m}"] for r in rows if not r["fail"]])))
                    for m in MS}
    for m in MS:
        if m == 75:
            continue
        agg[f"paired_fit_m{m}_vs_m75"] = I.paired(
            np.array([r[f"fit_m{m}"] for r in rows]),
            np.array([r["fit_m75"] for r in rows]), folds=folds, names=names)
        for sub, msk in (("FAIL18", fail), ("other108", ~fail)):
            agg[f"paired_fit_m{m}_vs_m75_{sub}"] = I.paired(
                np.array([r[f"fit_m{m}"] for r in rows])[msk],
                np.array([r["fit_m75"] for r in rows])[msk])
    # ---- LEAVE-FOLD-OUT choice of m (no in-sample tuning)
    lfo = np.empty(len(rows)); chosen = {}
    for f in np.unique(folds):
        tr = folds != f
        means = {m: float(np.mean([r[f"fit_m{m}"] for r, k in zip(rows, tr) if k])) for m in MS}
        mstar = min(MS, key=lambda m: means[m])
        chosen[int(f)] = mstar
        lfo[folds == f] = [r[f"fit_m{mstar}"] for r, k in zip(rows, folds == f) if k]
    agg["lfo_chosen_m"] = chosen
    agg["lfo_mean"] = float(lfo.mean())
    agg["lfo_fail"] = float(lfo[fail].mean())
    agg["lfo_other"] = float(lfo[~fail].mean())
    agg["paired_lfo_vs_m75"] = I.paired(lfo, np.array([r["fit_m75"] for r in rows]),
                                        folds=folds, names=names)
    print("\n m   fit_mean  FAIL18   other108")
    for m in MS:
        e = agg["fit"][m]
        print(f"{m:4d}  {e['mean']:7.3f}  {e['fail']:7.3f}  {e['other']:7.3f}")
    print(f"LFO   {agg['lfo_mean']:7.3f}  {agg['lfo_fail']:7.3f}  {agg['lfo_other']:7.3f}  chosen {chosen}")
    for m in MS:
        if m == 75:
            continue
        v = agg[f"paired_fit_m{m}_vs_m75"]
        print(f"m{m} vs m75: d={v['mean_diff']:+.4f} CI{[round(x,4) for x in v['ci95']]} "
              f"W/L {v['n_better']}/{v['n_worse']} drop10 {v['drop_top10_mean_diff']:+.4f}")
    v = agg["paired_lfo_vs_m75"]
    print(f"LFO vs m75: d={v['mean_diff']:+.4f} CI{[round(x,4) for x in v['ci95']]} "
          f"W/L {v['n_better']}/{v['n_worse']}")
    I.write(f"adv_mladder_{shard}" if shard is not None else "adv_mladder", dict(agg=agg, rows=rows))


if __name__ == "__main__":
    a = sys.argv
    main(int(a[1]), int(a[2]), a[3]) if len(a) > 3 else main()
