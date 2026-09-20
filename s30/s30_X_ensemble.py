"""S30 lane X, H-X2: the deposited "native" is model 1 of an NMR ensemble. How much of the
endpoint -- and how much of the FAIL18 tail -- is the ensemble's own spread?

Pre-registration: `s30/PREREG_S30_X.md` section 2 H-X2, committed at d4305d17 BEFORE this
file existed. My OWN registered prediction there is that this kills my own hypothesis:
mean floor ~0.6 A, explaining < 5% of a 6.28 A tail.

`core/geometry.py:808` `parse_pdb(..., model_index=0)` -- the benchmark scores against
model 1 only.  `ca_rmsd_to_ensemble` exists in two modules and is called by nothing.

The FLOOR, defined so it is not over-read: a method that predicted the ensemble's MEDOID
perfectly would still score `RMSD(medoid, model 1)` against the benchmark's native. That is
the irreducible cost of the model-1 convention, and it is ORACLE by construction (it uses
the deposited native), so every row is labelled ORACLE.

Reads PDB files and one stored bench artefact. No pipeline compute.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE)
OUT = os.path.join(BASE, "s30", "results", "s30_X_ensemble.json")
N_TAIL = 18


def ca_rmsd(a, b):
    a = a - a.mean(0)
    b = b - b.mean(0)
    u, s, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    s = s.copy()
    s[-1] *= d
    e = (a * a).sum() + (b * b).sum() - 2 * s.sum()
    return float(np.sqrt(max(e, 0.0) / len(a)))


def main():
    from core import geometry as geo
    from s7 import debias

    tg = list(debias.tuning_targets())
    bench = json.load(open(os.path.join(BASE, "bench_results",
                                        "baseline_tuning126.json")))["per_target"]
    end = {r["pdb"]: r["rmsd_arm"] for r in bench}      # BUILT CHAIN (contract rule 1)
    pool_best = {r["pdb"]: r["pool_best"] for r in bench}

    import glob
    rows, missing = [], []
    for p in tg:
        hits = (glob.glob(os.path.join(BASE, "pdbs", p.pdb + ".pdb"))
                + glob.glob(os.path.join(BASE, "pdbs_ext", p.pdb + ".pdb")))
        if not hits:
            missing.append(p.pdb)
            continue
        try:
            models = geo.parse_pdb_ensemble(hits[0])
        except Exception as e:                                    # pragma: no cover
            missing.append(p.pdb + ":" + type(e).__name__)
            continue
        cas = [np.asarray(m[2], float) for m in models
               if len(m[2]) == len(p.ca)]
        if len(cas) < 1:
            missing.append(p.pdb + ":len")
            continue
        X = np.stack(cas)
        nm = len(X)
        # pairwise spread
        pw = [ca_rmsd(X[i], X[j]) for i in range(nm) for j in range(i + 1, nm)]
        # medoid: the model minimising mean RMSD to the rest
        D = np.zeros((nm, nm))
        for i in range(nm):
            for j in range(i + 1, nm):
                D[i, j] = D[j, i] = ca_rmsd(X[i], X[j])
        med = int(np.argmin(D.mean(1))) if nm > 1 else 0
        rows.append({
            "pdb": p.pdb, "n": int(p.n), "n_models": nm,
            "spread_pairwise_ORACLE": float(np.mean(pw)) if pw else 0.0,
            "spread_to_model1_ORACLE": float(D[0].sum() / (nm - 1)) if nm > 1 else 0.0,
            "floor_medoid_to_model1_ORACLE": float(D[0, med]),
            "endpoint_builtchain": end.get(p.pdb),
            "pool_best_ORACLE": pool_best.get(p.pdb),
        })

    rows = [r for r in rows if r["endpoint_builtchain"] is not None]
    e = np.array([r["endpoint_builtchain"] for r in rows])
    fl = np.array([r["floor_medoid_to_model1_ORACLE"] for r in rows])
    sp = np.array([r["spread_pairwise_ORACLE"] for r in rows])
    nm = np.array([r["n_models"] for r in rows])

    order = np.argsort(-e)
    tail, rest = order[:N_TAIL], order[N_TAIL:]

    d_end = float(e[tail].mean() - e[rest].mean())
    d_floor = float(fl[tail].mean() - fl[rest].mean())
    d_spread = float(sp[tail].mean() - sp[rest].mean())

    out = {
        "prereg": "s30/PREREG_S30_X.md @ d4305d17 (H-X2)",
        "n_targets_resolved": len(rows), "missing": missing,
        "n_multi_model": int((nm > 1).sum()), "n_single_model": int((nm == 1).sum()),
        "n_models_mean": float(nm.mean()), "n_models_median": float(np.median(nm)),
        "floor_medoid_to_model1_ORACLE": {
            "mean": float(fl.mean()), "median": float(np.median(fl)),
            "sd": float(fl.std(ddof=1)), "max": float(fl.max()), "min": float(fl.min()),
            "frac_above_1A": float((fl > 1.0).mean()),
        },
        "spread_pairwise_ORACLE": {
            "mean": float(sp.mean()), "median": float(np.median(sp)),
            "max": float(sp.max()),
        },
        "endpoint_builtchain": {"mean": float(e.mean()), "median": float(np.median(e)),
                                "tail%d" % N_TAIL: float(e[tail].mean()),
                                "rest": float(e[rest].mean())},
        "corr_floor_endpoint": float(np.corrcoef(fl, e)[0, 1]),
        "corr_spread_endpoint": float(np.corrcoef(sp, e)[0, 1]),
        "corr_floor_poolbest_ORACLE": float(np.corrcoef(
            fl, [r["pool_best_ORACLE"] for r in rows])[0, 1]),
        "tail_enrichment": {
            "d_endpoint_tail_minus_rest": d_end,
            "d_floor_tail_minus_rest_ORACLE": d_floor,
            "d_spread_tail_minus_rest_ORACLE": d_spread,
            "share_of_tail_excess_explained_by_floor": float(d_floor / d_end),
            "share_of_tail_excess_explained_by_spread": float(d_spread / d_end),
        },
        "DECISION_RULE": ("PRICED AND DEAD at this endpoint if mean floor < 1.0 A AND the "
                          "floor explains < 25% of the FAIL18 excess; LIVE otherwise "
                          "(registered before the numbers)"),
        "rows": rows,
    }
    out["VERDICT"] = ("PRICED AND DEAD at this endpoint"
                      if (out["floor_medoid_to_model1_ORACLE"]["mean"] < 1.0
                          and abs(d_floor / d_end) < 0.25) else "LIVE")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1)
    os.replace(tmp, OUT)
    small = {k: v for k, v in out.items() if k != "rows"}
    json.dump(small, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
