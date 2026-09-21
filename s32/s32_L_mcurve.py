"""LANE L / L4 -- is m = 75 simply the wrong prefix size at 40-60 residues?

NOT PRE-REGISTERED.  Declared EXPLORATORY, and it is the one lane-L arm that proposes a
DEPLOYABLE change, so it is held to the stricter standard of the two.

THE MECHANISM IT TESTS.  `avg75` uniformly averages the 75 retained members.  At L ~ 13
those 75 are spread 4.03 A apart (mean pairwise CA-RMSD) and the average is a usable
object.  At L ~ 55 they are spread **10.25 A** apart and the average has a mean virtual
CA-CA bond of 1.80 A against a native 3.81 -- it has collapsed to less than half a
backbone.  Averaging is only a sane operator on a concentrated set, and concentration
degrades with length.  Project memory already records that the optimal prefix shrinks as
the objective improves (`operator-consumes-set-mean`: m* 500 -> 75 -> 20 -> 3-5).  m = 75
was fixed at peptide length and has never been revisited at any other.

TWO ARMS, AND THE DIFFERENCE BETWEEN THEM IS THE POINT.

    m_oracle   per-target argmin over the m grid.  **ORACLE / NOT DEPLOYABLE**, and it is
               mostly an ORDER STATISTIC over |grid| draws (memory: `grid-oracles-are-
               order-statistics`).  Reported only to bound the arm, never as the result.
    m_lfo      ONE global m chosen on four folds and applied to the held-out fifth, for
               all five folds.  **DEPLOYABLE.**  This is the number that counts.

BASIS: **CA point cloud**, not the built chain.  Stated on every number here.  The cloud is
a diagnostic (contract rule 0) and a cloud gain is NOT an endpoint gain: it must be carried
to the chain before it means anything.  The chain carry is named as the next step rather
than assumed, because at L ~ 55 the projection cost is NOT the ~0 it is at peptide length
(measured: +0.29 A for a dense average) and it could eat a cloud gain whole.

    python -m s32.s32_L_mcurve run
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
K = 500
GRID = (1, 2, 3, 5, 8, 12, 20, 30, 50, 75, 128, 250, 500)
SHIPPED_M = 75


def _curve(W, nat):
    """Cloud CA-RMSD of the uniform medoid-frame average of the top-m, for m in GRID."""
    from s12 import instrument as I
    P = I.pairwise_rmsd(W)
    out = []
    for m in GRID:
        m = min(m, len(W))
        C, _ = I.coordinate_average(W[:m], P[:m, :m])
        out.append(I.ca_rmsd(C, nat))
    return np.array(out)


def run(verbose=True):
    from core import geometry as geo
    from s12 import instrument as I
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    out = {"basis": "CA point cloud (NOT the built chain)", "grid": list(GRID),
           "shipped_m": SHIPPED_M}
    for kind in ("short", "long"):
        rows = []
        if kind == "short":
            tg = I.targets()
        else:
            tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
            cp.load()
        for k, t in enumerate(tg):
            if kind == "short":
                u = I.load_univ(t["pdb"])
                W = np.asarray(u["W"], float)[I.pool_idx(u, K)]
                nat = np.asarray(u["nat_ca"], float)
                fold = int(t["fold"])
            else:
                W, _, _ = LD._bank_long(t)
                p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
                if not os.path.exists(p):
                    p = glob.glob(os.path.join(BASE, "prots",
                                               t["pdb"].lower() + ".pdb"))[0]
                _, coords, _, _ = geo.native_coords_from_pdb(p)
                nat = np.asarray(coords["CA"], float)
                fold = int(t["fold"])
            c = _curve(W, nat)
            rows.append(dict(pdb=t["pdb"], n=int(len(nat)), fold=fold,
                             curve=[float(x) for x in c]))
            if verbose and (k + 1) % 20 == 0:
                print(f"  {kind} {k+1}/{len(tg)}", flush=True)
        C = np.array([r["curve"] for r in rows])
        folds = np.array([r["fold"] for r in rows])
        gi = list(GRID).index(SHIPPED_M)
        shipped = C[:, gi]
        # ORACLE per-target argmin -- an order statistic over len(GRID) draws.
        m_oracle = C.min(1)
        # DEPLOYABLE: global m fitted out-of-fold.
        lfo = np.empty(len(rows)); chosen = {}
        for f in sorted(set(folds.tolist())):
            tr = folds != f
            best = int(np.argmin(C[tr].mean(0)))
            chosen[int(f)] = GRID[best]
            lfo[~tr] = C[~tr, best]
        # A global m fitted on ALL targets, reported separately as the in-sample version.
        gbest = int(np.argmin(C.mean(0)))
        out[kind] = dict(
            n=len(rows), mean_len=float(np.mean([r["n"] for r in rows])),
            curve_mean=[float(x) for x in C.mean(0)],
            shipped_m75=float(shipped.mean()),
            m_oracle=float(m_oracle.mean()),
            m_lfo=float(lfo.mean()),
            m_chosen_per_fold=chosen,
            m_global_insample=GRID[gbest],
            m_global_insample_value=float(C[:, gbest].mean()),
            folds=[int(x) for x in folds],
            shipped_rows=[float(x) for x in shipped],
            lfo_rows=[float(x) for x in lfo],
            rows=rows)
    # the deployable contrast, with its MDE
    from s24 import stats_lib as ST
    out["contrast"] = {}
    for kind in ("short", "long"):
        d = out[kind]
        c = ST.compare(np.array(d["lfo_rows"]), np.array(d["shipped_rows"]),
                       folds=np.array(d["folds"]), label=f"{kind}: m_lfo vs m=75")
        out["contrast"][kind] = {k: c[k] for k in
                                 ("n", "mean_a", "mean_b", "effect", "median_effect",
                                  "se", "mde", "effect_over_mde", "ci95_iid",
                                  "n_better", "n_worse")}
        if "per_fold" in c:
            out["contrast"][kind]["per_fold"] = c["per_fold"]
    path = os.path.join(RESULTS, "L4_mcurve.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        for kind in ("short", "long"):
            d = out[kind]
            print(f"\n== {kind}  n={d['n']}  L~{d['mean_len']:.1f}   "
                  f"BASIS: CA CLOUD, not the built chain")
            print("   m:      " + "".join(f"{m:>8}" for m in GRID))
            print("   cloud:  " + "".join(f"{x:8.3f}" for x in d["curve_mean"]))
            print(f"   shipped m=75 {d['shipped_m75']:.4f} | LFO global m "
                  f"{d['m_lfo']:.4f} (per fold {d['m_chosen_per_fold']}) | "
                  f"ORACLE per-target {d['m_oracle']:.4f} (ORDER STATISTIC, not a result)")
            c = out["contrast"][kind]
            band = ("RESULT" if abs(c["effect_over_mde"]) >= 1 else
                    "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
            print(f"   LFO m vs m=75: {c['effect']:+.4f}  MDE {c['mde']:.4f}  "
                  f"{abs(c['effect_over_mde']):.2f}x  W/L {c['n_better']}/{c['n_worse']}"
                  f"  {band}   [CLOUD BASIS -- must be carried to the chain]")
        print("\nwrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
