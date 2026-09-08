"""SPRINT 13 -- ADVERSARIAL AUDIT 6b: is the gradient-variance ratio reproducible?

`s13/results/geo_pauli.json` reports, for 1A13 L=4 k=4 layers=3, raw AMBER
`ratio_meas_over_pred_exact = 0.9804`.  Recomputing the same quantity here with a different
theta draw gave 0.6605.  The ratio is an estimate from `n_theta` random parameter vectors,
and AMBER's table is spike-dominated, so a heavy-tailed variance estimate is exactly what
should be unstable.  If the ratio moves by 30 % between theta draws, "median ratio 0.993 /
0.915" needs a confidence interval before it is a headline.

This runs the identical computation at S independent theta seeds per cell and reports the
spread, for all four variants, and at two values of n_theta so the convergence rate is
visible.

    python -m s13.adv_pauli2
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s13 import geo_common as G                     # noqa: E402
from s13 import geo_pauli as GP                     # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")


def ratios(pdb, L, k, layers, n_theta, seed):
    n = int(round(np.log2(k))) * L
    N = 1 << n
    circ = G.circuit(n, layers=layers)
    P = circ.n_params()
    w = GP.popcount(n)
    rng = np.random.default_rng(seed)
    TH = rng.uniform(-np.pi, np.pi, size=(n_theta, P))
    s1 = np.zeros((P, N)); s2 = np.zeros((P, N))
    Js = np.zeros((n_theta, P, N))
    for t in range(n_theta):
        PR = circ.probs_batch(circ._shift_grid(TH[t], np.pi / 2))
        J = (PR[0::2] - PR[1::2]) / 2.0
        Js[t] = J
        Jh = GP.fwht(J)
        s1 += Jh; s2 += Jh ** 2
    VarS = s2 / n_theta - (s1 / n_theta) ** 2
    vw = np.zeros(n + 1)
    m = VarS.mean(0)
    for ww in range(n + 1):
        sel = (w == ww)
        vw[ww] = float(m[sel].mean()) if sel.any() else 0.0
    out = {}
    for v in G.VARIANTS:
        E = G.variant_table(pdb, L, k, v)
        c2 = (GP.fwht(E) / N) ** 2
        pred_x = float((c2[None, :] * VarS).sum(1).mean())
        pred_w = float((c2 * vw[w]).sum())
        g = np.einsum("tpx,x->tp", Js, E)
        meas = float(g.var(0).mean())
        out[v] = {"ratio_exact": meas / max(pred_x, 1e-300),
                  "ratio_weight_kernel": meas / max(pred_w, 1e-300)}
    return out


def main():
    cells = [("1A13", 4, 4, 3, 128), ("2BFI", 5, 4, 3, 96), ("1A1P", 4, 8, 3, 64)]
    seeds = [0, 1, 2, 3, 4, 5, 6, 7]
    res = []
    for pdb, L, k, layers, nt in cells:
        if not os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")):
            continue
        for nth in (nt, 2 * nt):
            got = {v: {"exact": [], "wk": []} for v in G.VARIANTS}
            for s in seeds:
                r = ratios(pdb, L, k, layers, nth, 7717 * (L * int(np.log2(k))) + 31 * layers + s)
                for v in G.VARIANTS:
                    got[v]["exact"].append(r[v]["ratio_exact"])
                    got[v]["wk"].append(r[v]["ratio_weight_kernel"])
            row = {"pdb": pdb, "L": L, "k": k, "layers": layers, "n_theta": nth,
                   "n_seeds": len(seeds)}
            for v in G.VARIANTS:
                a = np.array(got[v]["exact"]); b = np.array(got[v]["wk"])
                row[v] = {"ratio_exact_mean": float(a.mean()), "ratio_exact_sd": float(a.std(ddof=1)),
                          "ratio_exact_min": float(a.min()), "ratio_exact_max": float(a.max()),
                          "ratio_wk_mean": float(b.mean()), "ratio_wk_sd": float(b.std(ddof=1)),
                          "ratio_wk_min": float(b.min()), "ratio_wk_max": float(b.max()),
                          "ratio_exact_all": [float(x) for x in a]}
                print(f"  {pdb} L={L} k={k} nth={nth} {v:12s} exact {a.mean():.3f} +- {a.std(ddof=1):.3f} "
                      f"[{a.min():.3f},{a.max():.3f}]   wkernel {b.mean():.3f} +- {b.std(ddof=1):.3f}",
                      flush=True)
            res.append(row)
            json.dump({"what": "theta-seed stability of the gradient-variance ratio",
                       "cells": res}, open(os.path.join(RESULTS, "adv_pauli2.json"), "w"), indent=1)
    print("wrote adv_pauli2.json")


if __name__ == "__main__":
    main()
