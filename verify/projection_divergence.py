"""The consolidated projection selects a different degenerate branch on some targets.

`core/project.py` states that it "changes only how fast one L-BFGS-B iteration is
computed", and every science knob it names -- starts, maxiter, penalty, lam, tolerance --
is indeed untouched. But the ARITHMETIC of the objective did change: the inner loop uses a
scan builder that reproduces the reference to 5.3e-14 A rather than bit-exactly.

That is normally harmless. Here it is not, because the module's own docstring says the
problem is DEGENERATE: "a CA trace admits two ideal-geometry torsion solutions at
near-equal objective distance". When two minima are near-equal, a 1e-14 difference in the
objective is enough to change which one L-BFGS-B reports, and the emitted structure then
moves by angstroms rather than by epsilon.

This script measures that, per target, on the two arms actually run:

    baseline    every backend legacy    (s8.project)
    optimised   core.project, PROJECT_GRAD=fd  -- the REFERENCE's own gradient, so the
                gradient formula is not the variable here

Reported per target: how far apart the two arms' own structures are (Kabsch RMSD of one
against the other), and what that costs against the native.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
CACHE = os.path.join(_ROOT, "bench_results", "cache")


def load(key):
    out = {}
    for f in glob.glob(os.path.join(CACHE, key, "*.json")):
        r = json.load(open(f))
        out[r["pdb"]] = r
    return out


def main(base_key=None, opt_key=None):
    import protein_geometry as geo
    base_key = base_key or sys.argv[1]
    opt_key = opt_key or sys.argv[2]
    A, B = load(base_key), load(opt_key)
    common = sorted(set(A) & set(B))

    rows = []
    for p in common:
        a, b = A[p], B[p]
        row = {"pdb": p, "n": a.get("n")}
        for k in ("avg_ca", "fit_ca", "ca", "amber_ca"):
            if k in a and k in b:
                x, y = np.asarray(a[k], float), np.asarray(b[k], float)
                row[f"{k}_arms_apart_rmsd"] = float(geo.rmsd(x, y))
                row[f"{k}_raw_max_abs"] = float(np.max(np.abs(x - y)))
        for k in ("rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full"):
            if k in a and k in b:
                row[f"{k}_base"] = float(a[k])
                row[f"{k}_opt"] = float(b[k])
                row[f"{k}_delta"] = float(b[k] - a[k])
        rows.append(row)

    def col(name):
        return np.array([r[name] for r in rows if name in r], float)

    summary = {
        "baseline_key": base_key, "optimised_key": opt_key, "n_targets": len(rows),
        "upstream_of_projection_bit_identical": bool(
            np.all(col("avg_ca_raw_max_abs") == 0.0)),
        "targets_where_the_arms_disagree_by_more_than_1A": [
            r["pdb"] for r in rows if r.get("ca_arms_apart_rmsd", 0) > 1.0],
        "worst_arms_apart_rmsd_ca": float(col("ca_arms_apart_rmsd").max()),
        "rmsd_full_mean_base": float(col("rmsd_full_base").mean()),
        "rmsd_full_mean_opt": float(col("rmsd_full_opt").mean()),
        "rmsd_full_mean_delta": float(col("rmsd_full_delta").mean()),
        "rmsd_full_worst_abs_delta": float(np.abs(col("rmsd_full_delta")).max()),
        "rmsd_fit_mean_delta": float(col("rmsd_fit_delta").mean()),
        "n_rmsd_full_bit_identical": int((col("rmsd_full_delta") == 0.0).sum()),
    }
    out = {"summary": summary, "per_target": rows}

    print(f"{'pdb':7} {'arms apart (A)':>15} {'rmsd_full base':>15} "
          f"{'rmsd_full opt':>14} {'delta':>10}")
    for r in rows:
        print(f"{r['pdb']:7} {r.get('ca_arms_apart_rmsd', float('nan')):>15.4f} "
              f"{r['rmsd_full_base']:>15.6f} {r['rmsd_full_opt']:>14.6f} "
              f"{r['rmsd_full_delta']:>+10.6f}")
    print()
    print(json.dumps(summary, indent=2, sort_keys=True))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "projection_divergence.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=float)
    return out


if __name__ == "__main__":
    main()
