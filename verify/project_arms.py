"""Three arms, to separate a faithful consolidation from a changed optimisation.

With `core.project` live the consolidated pipeline ships an ANALYTIC gradient where the
reference used one-sided finite differences at eps=1e-5. Those are not the same
optimisation -- L-BFGS-B takes different steps and stops at a different point of the same
basin -- so a plain "baseline vs consolidated" diff conflates two very different claims:

    baseline  vs  opt_fd    IS the consolidation faithful?          (expect bit-identical)
    opt_fd    vs  opt_an    what does the analytic gradient cost?   (expect nonzero)
    baseline  vs  opt_an    what does the SHIPPED default emit?     (the headline)

Usage:  python verify/project_arms.py BASE_KEY FD_KEY AN_KEY
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(_ROOT, "bench_results", "cache")

ARRAYS = ("ca", "phi", "psi", "avg_ca", "fit_ca", "amber_ca")
SCALARS = ("shipped", "pool_best", "top_m_best", "top_m_mean",
           "rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full")


def load(key):
    d = os.path.join(CACHE, key)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "*.json"))):
        r = json.load(open(f))
        out[r["pdb"]] = r
    return out


def compare(A, B, label):
    common = sorted(set(A) & set(B))
    rep = {"label": label, "n": len(common)}
    if not common:
        rep["VACUOUS"] = True
        return rep
    sub_same = sum(list(A[p].get("sub", [])) == list(B[p].get("sub", []))
                   for p in common)
    rep["sub_order_identical"] = f"{sub_same}/{len(common)}"
    for k in ARRAYS:
        d, ident = [], 0
        for p in common:
            if k not in A[p] or k not in B[p]:
                continue
            x, y = np.asarray(A[p][k], float), np.asarray(B[p][k], float)
            if x.shape != y.shape:
                continue
            d.append(float(np.max(np.abs(x - y))))
            ident += int(np.array_equal(x, y))
        if d:
            rep[k] = {"worst_max_abs": max(d), "bit_identical": f"{ident}/{len(d)}"}
    for k in SCALARS:
        d, ident, va, vb = [], 0, [], []
        for p in common:
            x, y = A[p].get(k), B[p].get(k)
            if x is None or y is None:
                continue
            d.append(abs(float(x) - float(y)))
            ident += int(x == y)
            va.append(float(x))
            vb.append(float(y))
        if d:
            rep[k] = {"worst_abs": max(d), "bit_identical": f"{ident}/{len(d)}",
                      "mean_a": float(np.mean(va)), "mean_b": float(np.mean(vb)),
                      "mean_diff": float(np.mean(vb) - np.mean(va))}
    return rep


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    base, fd, an = sys.argv[1], sys.argv[2], sys.argv[3]
    A, F, N = load(base), load(fd), load(an)
    out = {"keys": {"baseline": base, "opt_fd": fd, "opt_analytic": an},
           "n": {"baseline": len(A), "opt_fd": len(F), "opt_analytic": len(N)},
           "consolidation_faithful__baseline_vs_opt_fd": compare(A, F, "baseline vs opt_fd"),
           "gradient_cost__opt_fd_vs_opt_analytic": compare(F, N, "opt_fd vs opt_analytic"),
           "shipped_default__baseline_vs_opt_analytic": compare(A, N, "baseline vs opt_an")}
    print(json.dumps(out, indent=2, sort_keys=True, default=float))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "project_arms.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=float)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
