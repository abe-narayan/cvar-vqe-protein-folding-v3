#!/usr/bin/env python
"""s27/analyse_arms.py -- statistics of the CVaR-VQE arm (`vqe_rows.jsonl`) and the built-chain
arm (`chain_rows.jsonl`) against DIS, plus the seed replication.  Writes `arms_summary.json`.
"""
from __future__ import annotations

import json
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


def compact(o):
    return dict(effect=o["effect"], se=o["se"], mde=o["mde"], x_mde=o.get("effect_over_mde"), ci_iid=o["ci95_iid"],
                ci_fold=o.get("ci95_fold"), folds_same_sign=o.get("folds_same_sign"), W=o["n_better"], L=o["n_worse"],
                T=o["n_tied"], verdict=o["verdict"])


def main():
    out = {}
    # ---- VQE arm
    by = defaultdict(dict)
    for line in open(os.path.join(RESULTS, "vqe_rows.jsonl"), encoding="utf-8"):
        r = json.loads(line); by[(r["config"], r["seed"])][r["pdb"]] = r
    pdbs = sorted(by[("DIS", 0)])
    folds = ST.pinned_folds(pdbs)
    dis0 = np.array([by[("DIS", 0)][p]["rmsd_vqe"] for p in pdbs])
    v = {}
    for (cfg, seed), rows in by.items():
        if len(rows) != len(pdbs):
            continue
        a = np.array([rows[p]["rmsd_vqe"] for p in pdbs])
        d = dict(seed=seed, mean_vqe=float(a.mean()),
                 mean_top75=float(np.mean([rows[p]["rmsd_top75"] for p in pdbs])),
                 mean_m=float(np.mean([rows[p]["m"] for p in pdbs])),
                 eps_max=float(np.max(np.abs([rows[p]["eps"] for p in pdbs]))),
                 mladder_mean=float(np.mean([rows[p]["mladder"] for p in pdbs])),
                 entropy_bits=float(np.mean([rows[p]["entropy_bits"] for p in pdbs])),
                 ess=float(np.mean([rows[p]["ess"] for p in pdbs])),
                 gate_pass_all=bool(all(rows[p]["gate_pass"] for p in pdbs)),
                 gate_equality_frac=float(np.mean([rows[p]["gate_equality"] for p in pdbs])),
                 secs_mean=float(np.mean([rows[p]["secs"] for p in pdbs])))
        d["vs_dis_vqe"] = compact(ST.compare(a, dis0, folds=folds, names=pdbs, label=f"vqe {cfg} s{seed} - DIS vqe s0"))
        v.setdefault(cfg, {})[str(seed)] = d
    # seed replication: seed-1 effect inside seed-0's fold CI?
    for cfg, ds in v.items():
        if "0" in ds and "1" in ds:
            e1 = ds["1"]["vs_dis_vqe"]["effect"]; ci0 = ds["0"]["vs_dis_vqe"]["ci_fold"]
            ds["replicates"] = bool(ci0[0] <= e1 <= ci0[1])
    out["vqe"] = v
    # ---- chain arm
    byc = defaultdict(dict)
    if os.path.exists(os.path.join(RESULTS, "chain_rows.jsonl")):
        for line in open(os.path.join(RESULTS, "chain_rows.jsonl"), encoding="utf-8"):
            r = json.loads(line); byc[r["config"]][r["pdb"]] = r
        pc = sorted(byc["DIS"]); fc = ST.pinned_folds(pc)
        disc = np.array([byc["DIS"][p]["rmsd_chain"] for p in pc])
        disk = np.array([byc["DIS"][p]["rmsd_cloud"] for p in pc])
        c = {}
        for cfg, rows in byc.items():
            if len(rows) != len(pc):
                continue
            a = np.array([rows[p]["rmsd_chain"] for p in pc]); k = np.array([rows[p]["rmsd_cloud"] for p in pc])
            c[cfg] = dict(n=len(pc), mean_chain=float(a.mean()), mean_cloud=float(k.mean()),
                          vs_dis_chain=compact(ST.compare(a, disc, folds=fc, names=pc, label=f"chain {cfg} - DIS")),
                          vs_dis_cloud=compact(ST.compare(k, disk, folds=fc, names=pc, label=f"cloud {cfg} - DIS")))
        out["chain"] = c
    ST.save_atomic(os.path.join(RESULTS, "arms_summary.json"), out, module_file=__file__)
    print("VQE ARM (genuine CVaR-VQE, 9 qubits, alpha 0.18, T 0.5), point cloud, vs DIS's own VQE arm (seed 0)")
    print(f"{'config':28s} {'vqe s0':>7} {'eff':>8} {'xMDE':>6} {'fold CI':>18} {'vqe s1':>7} {'eff s1':>8} {'rep':>4} {'m':>5} {'eps':>8} {'H bits':>6} {'gate':>5}")
    for cfg, ds in sorted(v.items(), key=lambda kv: kv[1]["0"]["mean_vqe"]):
        d0 = ds["0"]; a = d0["vs_dis_vqe"]; ci = a["ci_fold"]; d1 = ds.get("1")
        print(f"{cfg:28s} {d0['mean_vqe']:7.4f} {a['effect']:+8.4f} {a['x_mde'] or 0:6.2f} [{ci[0]:+.3f},{ci[1]:+.3f}] "
              f"{d1['mean_vqe'] if d1 else float('nan'):7.4f} {d1['vs_dis_vqe']['effect'] if d1 else float('nan'):+8.4f} "
              f"{'yes' if ds.get('replicates') else 'no':>4} {d0['mean_m']:5.1f} {d0['eps_max']:8.1e} {d0['entropy_bits']:6.2f} {'ok' if d0['gate_pass_all'] else 'FAIL':>5}")
    if "chain" in out:
        print("\nBUILT CHAIN (production projection) vs DIS, and the same sets on the point cloud")
        print(f"{'config':28s} {'n':>4} {'chain':>7} {'eff':>8} {'xMDE':>6} {'fold CI':>18} {'cloud':>7} {'eff':>8} {'xMDE':>6}")
        for cfg, d in sorted(out["chain"].items(), key=lambda kv: kv[1]["mean_chain"]):
            a = d["vs_dis_chain"]; b = d["vs_dis_cloud"]; ci = a["ci_fold"]
            print(f"{cfg:28s} {d['n']:4d} {d['mean_chain']:7.4f} {a['effect']:+8.4f} {a['x_mde'] or 0:6.2f} [{ci[0]:+.3f},{ci[1]:+.3f}] {d['mean_cloud']:7.4f} {b['effect']:+8.4f} {b['x_mde'] or 0:6.2f}")


if __name__ == "__main__":
    main()
