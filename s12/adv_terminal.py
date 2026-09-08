"""COORDINATOR RE-PRIORITISATION -- the TERMINAL x OBJECTIVE-QUALITY surface.

Section D established that the shipped terminal (uniform coordinate average of the top-75)
is a near-linear map of the SET MEAN and almost blind to the SET BEST.  If that is right,
every in-band ranking experiment this project has run was scored through an operator
structurally unable to pay for what it measured.  Three deliverables:

  T1  RE-PRICING.  Score the SAME ranking signals through many terminals (argmin = m1, and
      uniform averages at m = 1,2,3,5,10,25,50,75,150,300,500) on the SAME pools.  If a
      ranking contrast that is ~0 at m=75 is large at m=5, S9-8's "in-band discrimination is
      informationally impossible" was measured through the wrong instrument.
      Signals: shipped Bayes risk / consensus centrality / score+consensus fusion /
      random (3 seeds, the null) / ORACLE rr (the ceiling).

  T2  THE SURFACE.  Objective quality x m.  Quality is swept two ways, because they answer
      different questions:
        (a) MAE ladder -- score = L1 against (native matrix + i.i.d. noise) bisected to a
            requested MAE.  This is the obj agent's `iid` family, the well-behaved one.
        (b) RANK-SKILL ladder -- score = a*z(rr) + (1-a)*z(noise), which fixes Spearman
            skill directly.  This is the axis the m question actually lives on.
      Deliverable: argmin(m) as a function of quality, and where the shipped distogram sits.

  T3  THE DEPLOYABLE FALLBACK.  Is a FIXED small m better than m=75 under the SHIPPED
      objective?  Reported here on the coordinate average; `s12/adv_mladder.py` runs the
      same ladder through the real projection for the headline.

  T4  A NULL FOR MY OWN LAW.  The set-mean law was fitted on oracle-insertion perturbations.
      Here it is re-fitted, out of design, on the m-ladder sets (which vary cardinality
      rather than composition) and on random subsets.  If the coefficients move, the law is
      an artefact of how I perturbed the set.

ORACLE/DIAGNOSTIC: the rr-based scorers and the corrupted-native objectives are oracle arms
and are labelled as such.  The shipped/consensus/random arms are deployable.
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I

MS = [1, 2, 3, 5, 10, 25, 50, 75, 150, 300, 500]
MAES = [0.0, 0.25, 0.5, 1.0, 1.5, 2.339, 3.5, 6.0]
ALPHAS = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]


def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / max(x.std(), 1e-9)


def avg_of(W, idx, Pfull=None):
    """Uniform coordinate average of W[idx] after superposing on the subset's medoid."""
    if len(idx) == 1:
        return np.asarray(W[idx[0]], float)
    sub = np.asarray(W[idx], float)
    P = Pfull[np.ix_(idx, idx)] if Pfull is not None else None
    C, _ = I.coordinate_average(sub, P)
    return C


def main(limit=None, seed=0):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    rng = np.random.default_rng(seed)
    rows = []
    t0 = time.time()
    for t in tg:
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        W = np.asarray(u["W"][p], float); rr = np.asarray(u["rr"][p], float); nat = u["nat_ca"]
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(n)
        D = I.pair_dists(W, i, j)
        dtrue = I.pair_dists(nat[None], i, j)[0]
        sc = I.shipped_score(dg, D)
        # consensus centrality over the K=500 pool (deployable, native-free)
        P500 = I.pairwise_rmsd(W)
        cons = P500.mean(1)

        scorers = {
            "shipped": sc,
            "consensus": cons,
            "fuse_sc_cons": z(sc) + z(cons),
            "ORACLE_rr": rr,
        }
        for s in range(3):
            scorers[f"random{s}"] = rng.permutation(len(W)).astype(float)
        # (a) MAE ladder: L1 against a corrupted native matrix, i.i.d. noise at target MAE
        for mae in MAES:
            e = rng.normal(size=len(dtrue))
            e = e * (mae / max(np.abs(e).mean(), 1e-9))
            tgt = dtrue + e
            scorers[f"ORACLE_mae{mae}"] = np.abs(D - tgt[None, :]).mean(1)
        # (b) rank-skill ladder
        zn = rng.normal(size=len(W))
        for a in ALPHAS:
            scorers[f"ORACLE_a{a}"] = a * z(rr) + (1 - a) * z(zn)

        row = dict(pdb=pdb, n=n, fold=t["fold"], fail=int(pdb in I.FAIL18),
                   pool_best=float(rr.min()), pool_mean=float(rr.mean()))
        for nm, s in scorers.items():
            order = np.argsort(np.asarray(s, float), kind="stable")
            row[f"rho_{nm}"] = float(spearmanr(s, rr).statistic)
            for m in MS:
                idx = order[:m]
                row[f"{nm}_m{m}"] = float(I.ca_rmsd(avg_of(W, idx, P500), nat))
                row[f"{nm}_m{m}_setmean"] = float(rr[idx].mean())
                row[f"{nm}_m{m}_setbest"] = float(rr[idx].min())
        # T4: random subsets of size 75 (out-of-design perturbation for the law)
        for s in range(6):
            idx = rng.choice(len(W), 75, replace=False)
            row[f"rand75_{s}"] = float(I.ca_rmsd(avg_of(W, idx, P500), nat))
            row[f"rand75_{s}_setmean"] = float(rr[idx].mean())
            row[f"rand75_{s}_setbest"] = float(rr[idx].min())
        rows.append(row)
        print(f"{len(rows):3d}/{len(tg)} {pdb} shipped m1 {row['shipped_m1']:.2f} "
              f"m10 {row['shipped_m10']:.2f} m25 {row['shipped_m25']:.2f} "
              f"m75 {row['shipped_m75']:.2f} | orc m75 {row['ORACLE_rr_m75']:.2f} "
              f"{time.time()-t0:.0f}s", flush=True)

    names = [r["pdb"] for r in rows]
    fail = np.array([r["fail"] for r in rows], bool)
    folds = np.array([r["fold"] for r in rows])
    keys = [k[:-3] for k in rows[0] if k.endswith("_m75") and not k.endswith("setmean_m75")]
    keys = sorted({k.rsplit("_m", 1)[0] for k in rows[0]
                   if "_m" in k and k.rsplit("_m", 1)[-1].isdigit()})
    agg = {"n": len(rows), "MS": MS}
    table = {}
    for nm in keys:
        e = {}
        for m in MS:
            v = np.array([r[f"{nm}_m{m}"] for r in rows])
            e[m] = dict(mean=float(v.mean()), fail=float(v[fail].mean()),
                        other=float(v[~fail].mean()))
        best_m = min(MS, key=lambda m: e[m]["mean"])
        e["best_m"] = best_m
        e["best_m_fail"] = min(MS, key=lambda m: e[m]["fail"])
        e["best_m_other"] = min(MS, key=lambda m: e[m]["other"])
        e["rho"] = float(np.mean([r[f"rho_{nm}"] for r in rows])) if f"rho_{nm}" in rows[0] else None
        e["gain_best_vs_m75"] = e[best_m]["mean"] - e[75]["mean"]
        table[nm] = e
    agg["table"] = table
    # deployable m ladder paired against m=75
    for m in MS:
        if m == 75:
            continue
        agg[f"paired_shipped_m{m}_vs_m75"] = I.paired(
            np.array([r[f"shipped_m{m}"] for r in rows]),
            np.array([r["shipped_m75"] for r in rows]), folds=folds, names=names)
    # T1 re-pricing: contrast between two deployable rankers at each m
    for a, b in (("fuse_sc_cons", "shipped"), ("consensus", "shipped"),
                 ("shipped", "random0"), ("ORACLE_rr", "shipped")):
        agg[f"contrast_{a}_vs_{b}"] = {
            m: I.paired(np.array([r[f"{a}_m{m}"] for r in rows]),
                        np.array([r[f"{b}_m{m}"] for r in rows]), names=names) for m in MS}
    I.write("adv_terminal", dict(agg=agg, rows=rows))
    # console summary
    print("\n=== emitted (coordinate average) vs m ===")
    hdr = "scorer".ljust(20) + "rho".rjust(7) + "".join(f"m{m}".rjust(8) for m in MS) + "  best_m"
    print(hdr)
    for nm in keys:
        e = table[nm]
        print(nm.ljust(20) + (f"{e['rho']:+.3f}" if e["rho"] is not None else "   -  ").rjust(7)
              + "".join(f"{e[m]['mean']:8.3f}" for m in MS) + f"   {e['best_m']}")
    return agg


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
