#!/usr/bin/env python
"""s27/analyse_pool.py -- paired statistics of every configuration in `pool_rows.jsonl`
against the shipped score (DIS top-75), on the point-cloud endpoint, with the S24 statistics
library; the ranked table; the diagnostics; the nested leave-fold-out choice for the weighted
family; `best_of_k_within` for every grid.  Writes `s27/results/pool_summary.json`.
"""
from __future__ import annotations

import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "pool_rows.jsonl")
OUT = os.path.join(RESULTS, "pool_summary.json")


def load():
    by = defaultdict(dict)
    with open(ROWS, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            by[r["config"]][r["pdb"]] = r
    return by


def compact(o):
    return dict(effect=o["effect"], se=o["se"], mde=o["mde"], x_mde=o.get("effect_over_mde"),
                ci_iid=o["ci95_iid"], ci_fold=o.get("ci95_fold"), folds_same_sign=o.get("folds_same_sign"),
                W=o["n_better"], L=o["n_worse"], T=o["n_tied"], median=o["median_effect"], verdict=o["verdict"],
                conc_pctile=(o.get("concentration") or {}).get("pctile_in_null"), power=o.get("power"),
                type_m=o.get("type_m"))


def main():
    by = load()
    pdbs = sorted(by["DIS"])
    folds = ST.pinned_folds(pdbs)
    dis = np.array([by["DIS"][p]["rmsd"] for p in pdbs])
    rnd = np.array([by["DIS"][p]["rand_mean"] for p in pdbs])
    out = dict(n=len(pdbs), dis_mean=float(dis.mean()), rand_mean=float(rnd.mean()), configs={})
    for cfg, rows in by.items():
        if len(rows) != len(pdbs):
            continue
        a = np.array([rows[p]["rmsd"] for p in pdbs])
        c = ST.compare(a, dis, folds=folds, names=pdbs, label=f"{cfg} - DIS")
        cr = ST.compare(a, rnd, folds=folds, names=pdbs, label=f"{cfg} - random75")
        kind = rows[pdbs[0]]["kind"]
        d = dict(kind=kind, mean=float(a.mean()), median=float(np.median(a)), vs_dis=compact(c), vs_random=compact(cr),
                 overlap_dis=float(np.mean([rows[p]["overlap_dis"] for p in pdbs])),
                 rho_pool=float(np.nanmean([rows[p]["rho_pool"] for p in pdbs])),
                 rho_inband=float(np.nanmean([rows[p]["rho_inband"] for p in pdbs])),
                 tie_frac_top=float(np.mean([rows[p]["tie_frac_top"] for p in pdbs])),
                 n_distinct=float(np.mean([rows[p]["n_distinct"] for p in pdbs])),
                 per_fold={str(f): float(np.mean(a[folds == f] - dis[folds == f])) for f in sorted(set(folds))})
        if kind == "single":
            d["rho_dis"] = float(np.mean([rows[p]["rho_dis"] for p in pdbs]))
            d["cost_ms"] = float(np.nanmean([rows[p].get("cost_ms", float("nan")) for p in pdbs]))
        out["configs"][cfg] = d
    # permuted-control contrast for the pairs: real minus permuted
    for cfg in list(out["configs"]):
        if out["configs"][cfg]["kind"] == "pair" and (cfg + "~perm") in by:
            a = np.array([by[cfg][p]["rmsd"] for p in pdbs]); b = np.array([by[cfg + "~perm"][p]["rmsd"] for p in pdbs])
            out["configs"][cfg]["vs_perm"] = compact(ST.compare(a, b, folds=folds, names=pdbs, label=f"{cfg} - perm"))
    # nested leave-fold-out choice of the weight per complement (H3): weights {0.25, 0.5, 1.0}
    nested = {}
    comps = sorted({c[len("DIS+"):] for c in by if c.startswith("DIS+") and "*" not in c and "~" not in c and "adapt" not in c})
    for c in comps:
        arms = {w: f"DIS+{w}*{c}" if w != 1.0 else f"DIS+{c}" for w in (0.25, 0.5, 1.0)}
        if not all(k in by and len(by[k]) == len(pdbs) for k in arms.values()):
            continue
        Mx = np.column_stack([[by[arms[w]][p]["rmsd"] for p in pdbs] for w in (0.25, 0.5, 1.0)])
        held = np.empty(len(pdbs))
        chosen = {}
        for f in sorted(set(folds)):
            tr, te = folds != f, folds == f
            wbest = int(np.argmin(Mx[tr].mean(0) - dis[tr].mean()))
            held[te] = Mx[te, wbest]
            chosen[str(f)] = [0.25, 0.5, 1.0][wbest]
        cn = ST.compare(held, dis, folds=folds, names=pdbs, label=f"nested(w) {c} - DIS")
        bk = ST.best_of_k_within(Mx - dis[:, None])
        nested[c] = dict(held_out=compact(cn), chosen_w_per_fold=chosen,
                         full_leak_best=float((Mx - dis[:, None]).mean(0).min()),
                         best_of_k_within={k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                                           for k, v in bk.items() if isinstance(v, (int, float, np.floating, str))})
    out["nested_weight"] = nested
    # ranked table by mean RMSD
    ranked = sorted(out["configs"].items(), key=lambda kv: kv[1]["mean"])
    out["ranked"] = [dict(config=k, kind=v["kind"], mean=v["mean"], effect=v["vs_dis"]["effect"],
                          x_mde=v["vs_dis"]["x_mde"], ci_fold=v["vs_dis"]["ci_fold"], verdict=v["vs_dis"]["verdict"])
                     for k, v in ranked]
    ST.save_atomic(OUT, out, module_file=__file__)
    print(f"n={len(pdbs)}  DIS top-75 mean {dis.mean():.4f}  random-75 {rnd.mean():.4f}")
    print(f"{'rank':>4} {'config':30s} {'kind':11s} {'mean':>7} {'effect':>8} {'xMDE':>6} {'fold CI':>20} {'W/L':>8} verdict")
    for i, r in enumerate(out["ranked"][:60]):
        v = out["configs"][r["config"]]["vs_dis"]
        ci = v["ci_fold"] or [float('nan')] * 2
        print(f"{i+1:4d} {r['config']:30s} {r['kind']:11s} {r['mean']:7.4f} {r['effect']:+8.4f} {r['x_mde'] if r['x_mde'] is not None else 0:6.2f} "
              f"[{ci[0]:+.3f},{ci[1]:+.3f}] {v['W']:3d}/{v['L']:<3d} {r['verdict']}")
    print("...")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
