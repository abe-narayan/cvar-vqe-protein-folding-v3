"""Per-target, per-stage equivalence: the legacy root modules vs the consolidated `core`.

Compares the two arms' persisted per-target records field by field. Deterministic stages
are expected to be BIT-IDENTICAL; anything else is localised to the stage that introduced
it and quantified in the units the stage reports in.

Stages compared, in pipeline order:
    sub        the filtered top-m membership (set and ORDER)
    ca/phi/psi the retrieved pool's best member
    avg_ca     the coordinate average          (synthesis)
    fit_ca     the projected coordinates       (manifold projection)
    amber_ca   the post-AMBER coordinates
    rmsd_*     the final per-target numbers
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

ARRAYS = ("ca", "phi", "psi", "avg_ca", "fit_ca", "amber_ca", "amber_phi", "amber_psi")
SCALARS = ("shipped", "pool_best", "pool_mean", "top_m_best", "top_m_mean",
           "rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full", "d_rmsd",
           "amber_e0", "amber_e1", "amber_moved", "amber_strain_after", "n_windows")


def keys_for_arms():
    """The two arms' cache keys.

    Recomputing them from the CURRENT `Config` is unreliable while siblings are still
    adding fields -- a field added after a run makes the key move and the comparison
    silently finds nothing. So the keys may be passed on the command line, and when they
    are not, they are recomputed and the caller is warned if either directory is empty.
    """
    if len(sys.argv) >= 3:
        return sys.argv[1], sys.argv[2]
    import core
    from core.pipeline import PROD
    base = PROD.key({"forced": "legacy"})
    opt = PROD.key(core.backend_report())
    return base, opt


def load_arm(key):
    d = os.path.join(_ROOT, "bench_results", "cache", key)
    recs = {}
    if not os.path.isdir(d):
        return recs, d
    for fn in sorted(os.listdir(d)):
        if fn.endswith(".json"):
            with open(os.path.join(d, fn)) as fh:
                r = json.load(fh)
            recs[r["pdb"]] = r
    return recs, d


def main():
    base_key, opt_key = keys_for_arms()
    B, bdir = load_arm(base_key)
    O, odir = load_arm(opt_key)

    report = {"baseline_key": base_key, "baseline_dir": bdir, "n_baseline": len(B),
              "optimised_key": opt_key, "optimised_dir": odir, "n_optimised": len(O)}
    common = sorted(set(B) & set(O))
    report["n_common"] = len(common)
    report["only_baseline"] = sorted(set(B) - set(O))
    report["only_optimised"] = sorted(set(O) - set(B))

    stage = {}
    per_target = {}
    for pdb in common:
        b, o = B[pdb], O[pdb]
        row = {}

        # filtered membership: order AND set
        sb, so = list(b.get("sub", [])), list(o.get("sub", []))
        row["sub_order_identical"] = sb == so
        row["sub_set_identical"] = set(sb) == set(so)
        row["sub_symdiff"] = len(set(sb) ^ set(so)) // 2
        row["sub_len"] = (len(sb), len(so))

        for k in ARRAYS:
            if k not in b or k not in o:
                continue
            x, y = np.asarray(b[k], float), np.asarray(o[k], float)
            if x.shape != y.shape:
                row[k] = {"shape_mismatch": [list(x.shape), list(y.shape)]}
                continue
            d = np.abs(x - y)
            den = np.maximum(np.abs(x), 1e-300)
            row[k] = {"max_abs": float(d.max()),
                      "max_rel": float((d / den).max()),
                      "bit_identical": bool(np.array_equal(x, y))}

        for k in SCALARS:
            if k not in b or k not in o:
                continue
            x, y = b[k], o[k]
            if x is None or y is None:
                continue
            row[k] = {"baseline": x, "optimised": y,
                      "abs": abs(float(x) - float(y)),
                      "rel": abs(float(x) - float(y)) / max(abs(float(x)), 1e-300),
                      "bit_identical": x == y}
        per_target[pdb] = row

    # roll up worst-case per stage across targets
    for k in ARRAYS + SCALARS:
        vals = [per_target[p][k] for p in common
                if k in per_target[p] and isinstance(per_target[p][k], dict)
                and "max_abs" in per_target[p][k] or
                (k in per_target[p] and "abs" in per_target[p].get(k, {}))]
        if not vals:
            continue
        if "max_abs" in vals[0]:
            stage[k] = {"worst_max_abs": max(v["max_abs"] for v in vals),
                        "worst_max_rel": max(v["max_rel"] for v in vals),
                        "n_bit_identical": sum(v["bit_identical"] for v in vals),
                        "n": len(vals)}
        else:
            stage[k] = {"worst_abs": max(v["abs"] for v in vals),
                        "worst_rel": max(v["rel"] for v in vals),
                        "n_bit_identical": sum(v["bit_identical"] for v in vals),
                        "n": len(vals)}
    stage["sub"] = {
        "n_order_identical": sum(per_target[p]["sub_order_identical"] for p in common),
        "n_set_identical": sum(per_target[p]["sub_set_identical"] for p in common),
        "worst_symdiff": max([per_target[p]["sub_symdiff"] for p in common] or [0]),
        "n": len(common)}

    report["stage_rollup"] = stage
    report["per_target"] = per_target
    # A verdict over an empty intersection is VACUOUS, not a pass. Say so loudly.
    if not common:
        report["ALL_BIT_IDENTICAL"] = None
        report["VERDICT"] = ("VACUOUS: the two arms share no targets. Pass the two cache "
                             "keys explicitly: python verify/equiv_compare.py BASE OPT")
    else:
        ok = all((v.get("n_bit_identical") == v.get("n")) for k, v in stage.items()
                 if k != "sub" and "n_bit_identical" in v) \
            and stage["sub"]["n_order_identical"] == stage["sub"]["n"]
        report["ALL_BIT_IDENTICAL"] = bool(ok)
        report["VERDICT"] = ("BIT-IDENTICAL on all %d targets" % len(common) if ok
                             else "DIFFERENCES FOUND on %d targets" % len(common))

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "equivalence.json"), "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True, default=float)

    print(json.dumps({k: v for k, v in report.items() if k != "per_target"},
                     indent=2, sort_keys=True, default=float))
    return report

if __name__ == "__main__":
    main()
