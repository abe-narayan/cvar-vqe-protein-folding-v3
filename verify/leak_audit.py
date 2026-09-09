"""Leakage guard: NaN-poison the native and assert every deployable quantity survives.

The project's standard guard. If any native coordinate reaches a deployable stage, the
NaN propagates and the output changes (or becomes NaN); if the separation is real, every
emitted array is bit-identical to the unpoisoned run.

Poisoned: `target.ca`, `target.phi`, `target.psi` -- every native quantity the target
record carries. Compared: the FULL deployable record out of `run_target`, which is
retrieval, filtering, synthesis, projection and AMBER.
"""
from __future__ import annotations

import copy
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

NATIVE_FIELDS = ("ca", "phi", "psi")


def poison(target):
    """A copy of `target` with every native array replaced by NaN."""
    t = copy.deepcopy(target)
    n_poisoned = []
    for f in NATIVE_FIELDS:
        if not hasattr(t, f):
            continue
        v = np.asarray(getattr(t, f), float)
        bad = np.full(v.shape, np.nan)
        try:
            object.__setattr__(t, f, bad)
            n_poisoned.append(f)
        except Exception:
            try:
                setattr(t, f, bad)
                n_poisoned.append(f)
            except Exception:
                pass
    return t, n_poisoned


def flatten(rec):
    """Every emitted numeric quantity, as a flat name -> ndarray map."""
    out = {}
    for k, v in rec.items():
        if k == "amber" and isinstance(v, dict):
            for k2, v2 in v.items():
                if isinstance(v2, (int, float, np.floating)):
                    out[f"amber.{k2}"] = np.asarray([v2], float)
                elif isinstance(v2, np.ndarray):
                    out[f"amber.{k2}"] = np.asarray(v2, float)
            continue
        if isinstance(v, np.ndarray):
            out[k] = np.asarray(v, float)
        elif isinstance(v, (int, float, np.floating)) and not isinstance(v, bool):
            out[k] = np.asarray([v], float)
    return out


def main(pdbs=("1CS9", "1CB3", "1CEK")):
    import core
    from core import pipeline as P

    db = core.backend("data")
    folds = db.folds(P.PROD.n_folds)
    report = {"backends": core.backend_report(), "targets": {}}

    for pdb in pdbs:
        tg = [t for t in P.manifest("smoke8") if t.pdb == pdb]
        if not tg:
            tg = [db.by_pdb(pdb)] if hasattr(db, "by_pdb") else []
        if not tg:
            report["targets"][pdb] = {"error": "target not found"}
            continue
        t = tg[0]
        fold = int(folds[t.seq])

        clean, _, _ = P.run_target(t, fold, P.PROD)
        pois_t, poisoned = poison(t)
        dirty, _, _ = P.run_target(pois_t, fold, P.PROD)

        a, b = flatten(clean), flatten(dirty)
        keys = sorted(set(a) & set(b))
        diffs, nan_leaks = {}, []
        for k in keys:
            x, y = a[k], b[k]
            if x.shape != y.shape:
                diffs[k] = "shape mismatch"
                continue
            if np.isnan(y).any() and not np.isnan(x).any():
                nan_leaks.append(k)
            if not np.array_equal(x, y):
                d = np.abs(x - y)
                diffs[k] = {"max_abs": float(np.nanmax(d)),
                            "n_differing": int((x != y).sum())}
        report["targets"][pdb] = {
            "poisoned_fields": poisoned,
            "n_quantities_compared": len(keys),
            "quantities": keys,
            "n_bit_identical": len(keys) - len(diffs),
            "differences": diffs,
            "nan_leaked_into": nan_leaks,
            "CLEAN": not diffs and not nan_leaks,
            "only_baseline_keys": sorted(set(a) - set(b)),
            "only_poisoned_keys": sorted(set(b) - set(a)),
        }
        print(f"{pdb}: {len(keys)} quantities, "
              f"{'CLEAN' if report['targets'][pdb]['CLEAN'] else 'LEAK'}", flush=True)

    report["ALL_CLEAN"] = all(v.get("CLEAN") for v in report["targets"].values())
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "leak_audit.json"), "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True, default=str)
    print(json.dumps({"ALL_CLEAN": report["ALL_CLEAN"],
                      "per_target": {k: {kk: vv for kk, vv in v.items()
                                         if kk in ("poisoned_fields", "CLEAN",
                                                   "n_quantities_compared",
                                                   "n_bit_identical", "differences",
                                                   "nan_leaked_into")}
                                     for k, v in report["targets"].items()}},
                     indent=2, sort_keys=True, default=str))
    return report

if __name__ == "__main__":
    main()
