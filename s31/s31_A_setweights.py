"""S31 lane A -- is the derived readout's loss the SET or the WEIGHTS?

Production is: choose the SET by the DIS score (top-75), then weight UNIFORMLY.
The derived QP is: take the top-128 SET, then choose the WEIGHTS by <w,a_hat> - 1/2 w'Bw.

Those differ in two places at once, so the +0.13 A deficit is unattributed.  This decomposes it
on a 2x2: {top-75, top-128} x {uniform, derived QP}, with a_hat = const (the free half alone,
zero parameters, native-free).  a_hat = const is used so the comparison isolates the WEIGHTING
RULE and cannot be confounded by the quality model.

Basis: CA POINT CLOUD.  All four cells are native-free and deployable.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s8 import consensus2 as cc                                        # noqa: E402
from s31.s31_A_r1 import kabsch_rmsd_batch, medoid_uniform             # noqa: E402
from s31.s31_A_readout import min_quad_simplex                         # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s31", "results", "s31_A_setweights.json")


def cell(W, nat, n):
    """Return (uniform RMSD, MEB-weighted RMSD) for a candidate set, in its own medoid frame."""
    m = len(W)
    P = np.empty((m, m), float)
    for a in range(m):
        P[a] = kabsch_rmsd_batch(W, W[a])
    P = P.astype(np.float32).astype(float)
    b = medoid_uniform(P)
    X = cc.superpose_batch(W, W[b]).reshape(m, -1)
    sq = (X ** 2).sum(1)
    B = sq[:, None] + sq[None, :] - 2.0 * (X @ X.T)
    np.fill_diagonal(B, 0.0)
    B = np.maximum(B, 0.0)

    def rms(w):
        return float(kabsch_rmsd_batch(
            ((np.asarray(w, float) @ X).reshape(n, 3))[None], nat)[0])

    w_meb = min_quad_simplex(-B, np.zeros(m))
    return rms(np.full(m, 1.0 / m)), rms(w_meb), int((w_meb > 1e-6).sum())


def main(limit=None):
    files = UNIV if limit is None else UNIV[:limit]
    rows = []
    for i, f in enumerate(files):
        z = np.load(f)
        pdb, n, fold = str(z["pdb"]), int(z["n"]), int(z["fold"])
        pool = z["order"].astype(np.int64)[:500]
        nat = z["nat_ca"].astype(np.float64)
        W500 = z["W"][pool].astype(np.float64)
        DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
        order = np.argsort(DIS, kind="stable")
        u75, m75, s75 = cell(W500[order[:75]], nat, n)
        u128, m128, s128 = cell(W500[order[:128]], nat, n)
        rows.append({"pdb": pdb, "fold": fold, "u75": u75, "meb75": m75,
                     "u128": u128, "meb128": m128, "sup75": s75, "sup128": s128})
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(files)}", flush=True)

    folds = np.array([r["fold"] for r in rows], int)
    gl = sorted(set(folds.tolist()))
    g = lambda k: np.array([r[k] for r in rows], float)                 # noqa: E731

    def cmp(x, tag):
        mu = float(x.mean())
        se = float(np.array([x[folds == f].mean() for f in gl]).std(ddof=1) / np.sqrt(len(gl)))
        mde = 2.8016 * se
        return {"tag": tag, "mean": mu, "SE": se, "MDE": mde,
                "ratio_to_MDE": float(abs(mu) / mde) if mde > 0 else 0.0,
                "W": int((x < 0).sum()), "L": int((x > 0).sum()),
                "verdict": ("NOT A RESULT" if abs(mu) < 0.7 * mde else
                            "NOT MEASURED" if abs(mu) < mde else
                            ("BETTER" if mu < 0 else "WORSE"))}

    summ = {"n": len(rows),
            "cells": {k: float(g(k).mean()) for k in ("u75", "meb75", "u128", "meb128")},
            "meb_support": {"top75": float(g("sup75").mean()),
                            "top128": float(g("sup128").mean())},
            "comparisons": [
                cmp(g("u128") - g("u75"), "SET effect at uniform:  top128 - top75"),
                cmp(g("meb75") - g("u75"), "WEIGHT effect on top75: MEB - uniform"),
                cmp(g("meb128") - g("u128"), "WEIGHT effect on top128: MEB - uniform"),
                cmp(g("meb128") - g("u75"), "both:  MEB@128 - PROD75"),
                cmp(g("meb75") - g("meb128"), "SET effect at MEB: top75 - top128"),
            ]}
    with open(OUT, "w") as fh:
        json.dump(summ, fh, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
