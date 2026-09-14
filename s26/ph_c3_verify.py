"""s26/ph_c3_verify.py -- C3 stage 2 reduction check (lane PH, coordinator's ruling of 2026-09-14).

Lane P's best C2 rung is the shipped prior itself (L112), so `s26/results/p_best_rung_chains.json`
should carry the production emission. Verify per target, against `bench_results/cache/
1fc9f2dcf489e2fb/<pdb>.json`: `ca` equal to the cache's `ca`; `phi`/`psi`, wrapped to (-pi, pi]
(L114: delivered unwrapped), equal to the cache's `phi`/`psi` wrapped the same way; `rmsd_arm`
equal to the cache's. Report the max deviation and every target beyond 1e-6. Native-free (the
cache's `rmsd_arm` is compared as a stored number, not recomputed). CPU, seconds.

    python s26/ph_c3_verify.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s26 import ph_lib as L                                   # noqa: E402

DELIVERY = os.path.join(L.RESULTS, "p_best_rung_chains.json")
OUT = os.path.join(L.RESULTS, "ph_c3_stage2_verify.json")
TOL = 1e-6


def wrap(x):
    x = np.asarray(x, float)
    return ((x + math.pi) % (2 * math.pi)) - math.pi


def main() -> dict:
    with open(DELIVERY, encoding="utf-8") as fh:
        d = json.load(fh)
    rows = {r["pdb"]: r for r in d["rows"]}
    out_rows, bad = [], []
    for t in L.targets():
        p = t["pdb"]
        rec = L.prod_record_nativefree(p)
        r = rows[p]
        dca = float(np.abs(np.asarray(r["ca"], float) - np.asarray(rec["ca"], float)).max())
        dphi = float(np.abs(wrap(r["phi"]) - wrap(rec["phi"])).max())
        dpsi = float(np.abs(wrap(r["psi"]) - wrap(rec["psi"])).max())
        # the cached record's rmsd_arm is a stored number; compared as such, not recomputed
        with open(os.path.join(L.PROD_DIR, f"{p}.json"), encoding="utf-8") as fh:
            stored = json.load(fh)["rmsd_arm"]
        drm = abs(float(r["rmsd_arm"]) - float(stored))
        row = {"pdb": p, "rung": r.get("rung"), "d_ca_max": dca, "d_phi_max": dphi, "d_psi_max": dpsi,
               "d_rmsd_arm": drm, "max_abs_torsion_delivered": float(r.get("max_abs_torsion", np.nan)),
               "identical": bool(max(dca, dphi, dpsi, drm) < TOL)}
        out_rows.append(row)
        if not row["identical"]:
            bad.append(p)
    summ = {"n": len(out_rows), "n_identical": sum(r["identical"] for r in out_rows), "differing": bad,
            "max_d_ca": max(r["d_ca_max"] for r in out_rows), "max_d_phi": max(r["d_phi_max"] for r in out_rows),
            "max_d_psi": max(r["d_psi_max"] for r in out_rows), "max_d_rmsd_arm": max(r["d_rmsd_arm"] for r in out_rows),
            "rungs": sorted(set(str(r["rung"]) for r in out_rows)),
            "max_abs_torsion_delivered": max(r["max_abs_torsion_delivered"] for r in out_rows),
            "delivery_provenance": d.get("provenance"), "delivery_mean_rmsd_arm": d.get("mean_rmsd_arm"),
            "tol": TOL}
    print(json.dumps(summ, indent=1, default=str))
    L.save(OUT, {"what": "C3 stage 2 reduction check: lane P's best-rung delivery vs the production cache",
                 "rows": out_rows, "summary": summ}, rows=out_rows,
           complete_keys=("pdb", "d_ca_max", "d_phi_max", "d_psi_max", "d_rmsd_arm", "identical"),
           n_expected=126, module_file=__file__)
    return summ


if __name__ == "__main__":
    main()
