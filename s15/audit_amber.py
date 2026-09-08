"""S15 AUDIT / Part D -- independent verification of the ORACLE-CONDITIONED AMBER subset.

Sprint 14 recorded: "the cached AMBER subset in `s13/results/qarch_enum_*.npz` is 40%
ORACLE-CONDITIONED (0.401 A better than its space); only `amber_kind == 0` rows are an
unbiased sample."  This module verifies that from the npz files themselves and from the
sampling code in `s13/qarch_enum.py`, and then asks the follow-up question the record does
not: is `amber_kind == 0` REALLY unbiased?

    python -m s15.audit_amber
"""
from __future__ import annotations
import os, sys, json, glob
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)

# Declared in s13/qarch_enum.py
N_AMBER_UNIFORM, N_AMBER_PRIOR, N_AMBER_BAND = 1200, 1200, 600


def files():
    out = sorted(glob.glob(os.path.join(ROOT, "s13", "results", "qarch_enum_*.npz")))
    out += sorted(glob.glob(os.path.join(ROOT, "s14", "cache", "qarch_enum_*.npz")))
    out += sorted(glob.glob(os.path.join(ROOT, "s14", "cache", "obj_enum_*.npz")))
    return out


def main():
    rows = []
    for f in files():
        z = np.load(f, allow_pickle=True)
        if "amber_kind" not in z.files:
            rows.append({"file": os.path.basename(f), "skipped": "no amber_kind",
                         "keys": list(z.files)[:12]})
            continue
        pdb = str(z["pdb"]); rmsd = np.asarray(z["rmsd"], np.float64)
        ai = np.asarray(z["amber_idx"], int); ak = np.asarray(z["amber_kind"], int)
        snap = int(z["snap_index"])
        r_sub = rmsd[ai]
        n = len(ai)
        frac = {int(k): float((ak == k).mean()) for k in np.unique(ak)}
        # 'oracle-conditioned' = drawn from the true-RMSD-sorted band (kind 2)
        f_oracle = float((ak == 2).mean())
        d_all = float(r_sub.mean() - rmsd.mean())
        d_k0 = float(r_sub[ak == 0].mean() - rmsd.mean())
        d_k1 = float(r_sub[ak == 1].mean() - rmsd.mean())
        d_k2 = float(r_sub[ak == 2].mean() - rmsd.mean())
        # is snap (the ORACLE nearest-native config) inside the 'unbiased' kind-0 rows?
        snap_pos = np.where(ai == snap)[0]
        snap_kind = int(ak[snap_pos[0]]) if len(snap_pos) else None
        # band size the code used
        band_size = max(2000, len(rmsd) // 100)
        band_thresh = float(np.sort(rmsd)[band_size - 1])
        in_band_k0 = float((r_sub[ak == 0] <= band_thresh).mean())
        rows.append({
            "file": os.path.basename(f), "pdb": pdb, "n_configs": int(len(rmsd)),
            "n_amber": n, "kind_fractions": frac,
            "frac_oracle_conditioned_kind2": f_oracle,
            "space_mean_rmsd": float(rmsd.mean()),
            "subset_mean_rmsd": float(r_sub.mean()),
            "delta_full_subset_minus_space": d_all,
            "delta_kind0_minus_space": d_k0,
            "delta_kind1_minus_space": d_k1,
            "delta_kind2_minus_space": d_k2,
            "subset_min_rmsd": float(r_sub.min()),
            "kind0_min_rmsd": float(r_sub[ak == 0].min()),
            "space_min_rmsd": float(rmsd.min()),
            "snap_index_kind": snap_kind,
            "snap_rmsd": float(rmsd[snap]),
            "snap_is_space_argmin": bool(snap == int(np.argmin(rmsd))),
            "band_size_used": int(band_size), "band_rmsd_threshold": band_thresh,
            "kind0_frac_inside_band": in_band_k0,
            "kind0_expected_frac_inside_band": float(band_size / len(rmsd)),
            "n_amber_finite": int(np.isfinite(np.asarray(z["amber_total"], float)).sum()),
        })
        print(json.dumps(rows[-1]), flush=True)
        with open(os.path.join(OUT, "audit_amber_oracle.json"), "w") as fh:
            json.dump({"rows": rows, "complete": False}, fh, indent=1)

    good = [r for r in rows if "pdb" in r]
    agg = {
        "n_files": len(good),
        "mean_frac_oracle_conditioned": float(np.mean([r["frac_oracle_conditioned_kind2"] for r in good])),
        "mean_delta_full_subset_minus_space": float(np.mean([r["delta_full_subset_minus_space"] for r in good])),
        "mean_delta_kind0_minus_space": float(np.mean([r["delta_kind0_minus_space"] for r in good])),
        "mean_delta_kind1_minus_space": float(np.mean([r["delta_kind1_minus_space"] for r in good])),
        "mean_delta_kind2_minus_space": float(np.mean([r["delta_kind2_minus_space"] for r in good])),
        "snap_kinds": sorted({r["snap_index_kind"] for r in good}),
        "n_files_with_snap_in_kind0": int(sum(r["snap_index_kind"] == 0 for r in good)),
        "n_files_where_snap_is_space_argmin": int(sum(r["snap_is_space_argmin"] for r in good)),
        "kind0_frac_inside_band_mean": float(np.mean([r["kind0_frac_inside_band"] for r in good])),
        "kind0_expected_frac_inside_band_mean": float(np.mean([r["kind0_expected_frac_inside_band"] for r in good])),
    }
    print(json.dumps(agg, indent=1))
    with open(os.path.join(OUT, "audit_amber_oracle.json"), "w") as fh:
        json.dump({"rows": rows, "aggregate": agg, "complete": True,
                   "declared_constants": {"N_AMBER_UNIFORM": N_AMBER_UNIFORM,
                                          "N_AMBER_PRIOR": N_AMBER_PRIOR,
                                          "N_AMBER_BAND": N_AMBER_BAND}}, fh, indent=1)
    print("wrote", os.path.join(OUT, "audit_amber_oracle.json"))


if __name__ == "__main__":
    main()
