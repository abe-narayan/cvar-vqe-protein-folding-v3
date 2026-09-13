"""s26/e_reproduce.py -- reproduce the production means from the persisted records.

Lane E (Examiner), Sprint 26, EXAMINATION.md section H (campaign Part 2.3).  Reads the 126
production records in `bench_results/cache/1fc9f2dcf489e2fb/` and re-scores every stored
structure against its native through `s12.instrument.ca_rmsd` -- a different RMSD
implementation from `core.geometry.kabsch_rmsd_batch`, which wrote the records.  Natives
come from `core.backend("data").load()` and are put through the same float32 round trip as
`core.pipeline._q` (Config.reference_precision=True), because that is what `label()` scored
against.

Reports: the built-chain mean (`ca`, rmsd_arm), the point-cloud mean (`avg_ca`, rmsd_avg),
the lam=0 projection (`fit_ca`, rmsd_fit), the AMBER emission (`amber_ca`, rmsd_full), the
max per-target disagreement against the stored values, T030's (1S9Z) per-target rmsd_arm,
and the virtual Ca-Ca bond statistics behind the "22.3% contracted" claim.  Also reads the
results-lab leaderboard for the production row so the 3.2126 vs 3.2148 pair is stated.

Reads natives: this is a REPORTING script, not an operator; nothing here selects anything.

    python s26/jobrun.py --agent E --tag CPU --name e_reproduce --est-ram 0.8 -- \
        python s26/e_reproduce.py

Writes `s26/results/e_reproduce.json` (provenance-stamped, complete gated on 126 rows).
"""
from __future__ import annotations

import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np                                                          # noqa: E402

PROD_KEY = "1fc9f2dcf489e2fb"
CACHE = os.path.join(ROOT, "bench_results", "cache", PROD_KEY)
NEED = ("pdb", "n", "fold", "rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full",
        "stored_rmsd_avg", "stored_rmsd_fit", "stored_rmsd_arm", "stored_rmsd_full",
        "bond_mean_avg", "bond_min_avg", "bond_mean_chain", "bond_mean_native")


def bonds(ca):
    return np.linalg.norm(np.diff(np.asarray(ca, float), axis=0), axis=1)


def main():
    import core
    from core import pipeline as pl
    from s12 import instrument as I
    from s24 import stats_lib as ST

    db = core.backend("data")
    natives = {p.pdb: p for p in db.load()}
    folds = db.folds(5)
    cfg = pl.PROD
    files = sorted(glob.glob(os.path.join(CACHE, "*.json")))
    print(f"records: {len(files)} under {os.path.relpath(CACHE, ROOT)}")
    rows = []
    fold_mismatch = []
    for f in files:
        with open(f) as fh:
            r = json.load(fh)
        p = natives[r["pdb"]]
        nat = pl._q(np.asarray(p.ca, float), cfg)          # the reference's float32 round trip
        if int(folds[p.seq]) != int(r["fold"]):
            fold_mismatch.append(r["pdb"])
        row = {"pdb": r["pdb"], "n": int(r["n"]), "fold": int(r["fold"]), "seq": r["seq"],
               "rmsd_avg": I.ca_rmsd(r["avg_ca"], nat),
               "rmsd_fit": I.ca_rmsd(r["fit_ca"], nat),
               "rmsd_arm": I.ca_rmsd(r["ca"], nat),
               "rmsd_full": (I.ca_rmsd(r["amber_ca"], nat) if r.get("amber_ca") is not None
                             else None),
               "stored_rmsd_avg": r["rmsd_avg"], "stored_rmsd_fit": r["rmsd_fit"],
               "stored_rmsd_arm": r["rmsd_arm"], "stored_rmsd_full": r["rmsd_full"],
               "stored_shipped": r["shipped"], "stored_pool_best": r["pool_best"],
               "stored_top_m_best": r["top_m_best"],
               "bond_mean_avg": float(bonds(r["avg_ca"]).mean()),
               "bond_min_avg": float(bonds(r["avg_ca"]).min()),
               "bond_mean_chain": float(bonds(r["ca"]).mean()),
               "bond_sd_chain": float(bonds(r["ca"]).std()),
               "bond_mean_native": float(bonds(p.ca).mean()),
               "bond_mean_amber": (float(bonds(r["amber_ca"]).mean())
                                   if r.get("amber_ca") is not None else None),
               "cfg_key": r["cfg_key"], "dev_mode": r["dev_mode"], "amber_err": r.get("amber_err")}
        rows.append(row)
    rows.sort(key=lambda x: x["pdb"])
    n = len(rows)

    def col(k):
        return np.array([x[k] for x in rows if x[k] is not None], float)

    summary = {}
    for k in ("rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full"):
        fresh, stored = col(k), col("stored_" + k)
        summary[k] = {"n": int(len(fresh)), "mean_fresh": float(fresh.mean()),
                      "mean_stored": float(stored.mean()),
                      "max_abs_disagreement": float(np.max(np.abs(fresh - stored))),
                      "median": float(np.median(fresh)), "sd": float(fresh.std(ddof=1)),
                      "frac_under_2": float((fresh < 2.0).mean())}
    for k in ("stored_shipped", "stored_pool_best", "stored_top_m_best"):
        summary[k] = {"mean": float(col(k).mean())}
    gap = col("rmsd_arm") - col("rmsd_avg")
    summary["projection_gap_arm_minus_avg"] = {
        "mean": float(gap.mean()), "median": float(np.median(gap)), "sd": float(gap.std(ddof=1)),
        "q1": float(np.percentile(gap, 25)), "q3": float(np.percentile(gap, 75)),
        "min": float(gap.min()), "max": float(gap.max()),
        "n_chain_better": int((gap < 0).sum())}
    bavg = col("bond_mean_avg")
    summary["virtual_bond"] = {
        "avg_ca_mean_bond": float(bavg.mean()), "avg_ca_worst_bond": float(col("bond_min_avg").min()),
        "avg_ca_targets_mean_under_3.4": int((bavg < 3.4).sum()),
        "chain_mean_bond": float(col("bond_mean_chain").mean()),
        "chain_bond_sd_max": float(col("bond_sd_chain").max()),
        "native_mean_bond": float(col("bond_mean_native").mean()),
        "contraction_pct": float(100 * (1 - bavg.mean() / col("bond_mean_native").mean())),
        "amber_mean_bond": float(col("bond_mean_amber").mean())}
    t030 = [x for x in rows if x["pdb"] == "1S9Z"][0]
    summary["T030_1S9Z"] = {k: t030[k] for k in ("rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full",
                                                 "stored_rmsd_arm", "stored_pool_best")}
    summary["fold_mismatch_vs_pinned_folds"] = fold_mismatch
    summary["cfg_keys"] = sorted(set(x["cfg_key"] for x in rows))
    summary["dev_modes"] = sorted(set(x["dev_mode"] for x in rows))
    summary["amber_errors"] = [x["pdb"] for x in rows if x["amber_err"]]

    # the results-lab leaderboard, for the 3.2126 rebuild figure
    lb_path = os.path.join(ROOT, "results", "summary", "leaderboard.json")
    lb_prod = None
    try:
        with open(lb_path) as fh:
            lb = json.load(fh)
        cand = lb.get("rows") if isinstance(lb, dict) else lb
        if isinstance(cand, list):
            for r in cand:
                if isinstance(r, dict) and r.get("configuration") == "production":
                    lb_prod = {k: v for k, v in r.items() if not isinstance(v, (list, dict))}
        if lb_prod is None and isinstance(lb, dict):
            lb_prod = {"_keys": list(lb.keys())}
    except Exception as exc:                                                 # noqa: BLE001
        lb_prod = {"error": str(exc)}
    summary["leaderboard_production_row"] = lb_prod

    print(json.dumps(summary, indent=1, default=float))
    print()
    for k in ("rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full"):
        s = summary[k]
        print(f"{k:10} fresh mean {s['mean_fresh']:.6f}  stored mean {s['mean_stored']:.6f}  "
              f"max per-target |diff| {s['max_abs_disagreement']:.3e}  (n={s['n']})")
    print(f"T030 1S9Z rmsd_arm fresh {t030['rmsd_arm']:.6f} stored {t030['stored_rmsd_arm']:.6f}")
    out = {"summary": summary, "rows": rows, "n_records": n, "prod_key": PROD_KEY,
           "rmsd_implementation": "s12.instrument.ca_rmsd",
           "native_source": "core.backend('data').load() through core.pipeline._q (float32 round trip)"}
    path = os.path.join(ROOT, "s26", "results", "e_reproduce.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ST.save_atomic(path, out, complete_keys=NEED, rows=rows, n_expected=126, module_file=__file__)
    print(f"wrote {os.path.relpath(path, ROOT)}  rss {pl.proc_rss()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
