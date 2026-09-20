#!/usr/bin/env python
"""s29/s29_M_F1_analyse.py -- the statistics for F1, exactly as `s29/PREREG_S29_M_F1.md` registers.

Primary: LOG - PROD on the BUILT CHAIN, paired per target, ST.compare with the pinned folds,
decided on the fold-clustered CI. Everything else is a control or a decomposition and says so.

    python s29/s29_M_F1_analyse.py [--probe]
Artefacts: s29/results/s29_M_F1_summary.json, s29/results/s29_M_F1_fmt.txt
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s29_M_F1_rows.jsonl")
ARMS = ("PROD", "LOG", "LOGW", "L2RISK", "LOGPERM")


def load():
    rows = [json.loads(l) for l in open(ROWS, encoding="utf-8")]
    by = {}
    for r in rows:
        by.setdefault(r["arm"], {})[r["pdb"]] = r
    common = set.intersection(*[set(v) for v in by.values()]) if by else set()
    pdbs = sorted(common)
    return by, pdbs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    a = ap.parse_args()
    by, pdbs = load()
    n = len(pdbs)
    folds = ST.pinned_folds(pdbs)
    out = {"n": n, "probe": bool(a.probe), "pdbs": pdbs}
    lines = []

    def col(arm, key):
        return np.array([by[arm][p][key] for p in pdbs], float)

    out["means"] = {arm: dict(chain=float(col(arm, "rmsd_chain").mean()),
                              cloud=float(col(arm, "rmsd_cloud").mean()),
                              bond_cloud=float(col(arm, "bond_cloud").mean()),
                              rg_cloud=float(col(arm, "rg_cloud").mean()),
                              overlap_prod=float(col(arm, "overlap_prod").mean()),
                              tie_frac=float(col(arm, "tie_frac_at_cut").mean()),
                              cos_err_vs_prod=float(np.nanmean(col(arm, "cos_err_vs_prod"))))
                    for arm in ARMS if arm in by}
    out["rg_native"] = float(col("PROD", "rg_native").mean())
    out["floor_rate"] = float(col("LOG", "floor_rate").mean())

    # ---- the endpoint contrasts. PRIMARY FIRST, and labelled.
    out["contrasts"] = {}
    for arm, label, tag in (("LOG", "PRIMARY  LOG - PROD (built chain)", "primary"),
                            ("LOGW", "control  LOGW - PROD (built chain)", "control"),
                            ("L2RISK", "control  L2RISK - PROD (built chain)", "control"),
                            ("LOGPERM", "control  LOGPERM - PROD (built chain)", "control")):
        if arm not in by:
            continue
        c = ST.compare(col(arm, "rmsd_chain"), col("PROD", "rmsd_chain"), folds=folds,
                       names=pdbs, label=label, seed_parts=("s29M", "F1"))
        out["contrasts"][arm + "_chain"] = c
        lines.append(ST.fmt(c))
    for arm in ("LOG",):
        c = ST.compare(col(arm, "rmsd_cloud"), col("PROD", "rmsd_cloud"), folds=folds, names=pdbs,
                       label="intermediate  %s - PROD (POINT CLOUD, not the result)" % arm,
                       seed_parts=("s29M", "F1"))
        out["contrasts"][arm + "_cloud"] = c
        lines.append(ST.fmt(c))
    if "LOGW" in by:
        c = ST.compare(col("LOG", "rmsd_chain"), col("LOGW", "rmsd_chain"), folds=folds, names=pdbs,
                       label="decomposition  LOG - LOGW (built chain): what the shipped weight costs",
                       seed_parts=("s29M", "F1"))
        out["contrasts"]["LOG_minus_LOGW"] = c
        lines.append(ST.fmt(c))

    # ---- MECHANISM (rule 18): contraction, measured beside the outcome
    c = ST.compare(col("LOG", "bond_cloud"), col("PROD", "bond_cloud"), folds=folds, names=pdbs,
                   label="MECHANISM  LOG - PROD mean virtual CA-CA bond of the emitted cloud (A; "
                         "native 3.8122; POSITIVE = LOG contracts LESS)", seed_parts=("s29M", "F1"))
    out["contrasts"]["bond_cloud"] = c
    lines.append(ST.fmt(c))

    # ---- the ORACLE parallel-bias check (S24 L2/L3)
    cosv = col("LOG", "cos_err_vs_prod")
    nn = np.array([by["PROD"][p]["n"] for p in pdbs], int)
    rand_ref = float(np.mean(1.0 / np.sqrt(3 * nn - 6)))
    out["parallel_bias"] = dict(mean_cos=float(np.nanmean(cosv)), median_cos=float(np.nanmedian(cosv)),
                                min_cos=float(np.nanmin(cosv)), max_cos=float(np.nanmax(cosv)),
                                random_reference=rand_ref,
                                note="ORACLE DIAGNOSTIC: cos between LOG's and PROD's error vectors "
                                     "against the native, rigid body removed. Near 1 = the same answer.")

    # ---- FAIL18 / 108 with a random-18 null (ORACLE partition, labelled)
    d = col("LOG", "rmsd_chain") - col("PROD", "rmsd_chain")
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    if f18.sum() >= 2:
        rng = np.random.default_rng(29001)
        null = np.array([d[rng.choice(len(d), size=int(f18.sum()), replace=False)].mean()
                         for _ in range(2000)])
        obs = float(d[f18].mean())
        out["fail18"] = dict(n_fail18=int(f18.sum()), effect_fail18=obs,
                             effect_rest=float(d[~f18].mean()),
                             null_p10=float(np.percentile(null, 10)),
                             null_p50=float(np.percentile(null, 50)),
                             null_p90=float(np.percentile(null, 90)),
                             percentile_of_obs=float((null < obs).mean()),
                             note="FAIL18 is an ORACLE partition (s29/CONVENIENCE_CHOICES C8); "
                                  "it selects nothing and is reported as a split only.")

    # ---- the registered withdrawal diagnostics
    out["diagnostics"] = dict(
        mean_top75_overlap_LOG_vs_PROD=float(col("LOG", "overlap_prod").mean()),
        posterior_floor_rate=float(col("LOG", "floor_rate").mean()),
        mean_tie_frac_at_cut_LOG=float(col("LOG", "tie_frac_at_cut").mean()))

    txt = "\n\n".join(lines)
    open(os.path.join(RESULTS, "s29_M_F1_fmt.txt"), "w", encoding="utf-8").write(txt)
    ST.save_atomic(os.path.join(RESULTS, "s29_M_F1_summary.json"), out, module_file=__file__)
    print("n =", n, " (probe)" if a.probe else " (full)")
    print("means (built chain):", {k: round(v["chain"], 4) for k, v in out["means"].items()})
    print("means (point cloud):", {k: round(v["cloud"], 4) for k, v in out["means"].items()})
    print("mean cloud bond    :", {k: round(v["bond_cloud"], 4) for k, v in out["means"].items()},
          " native 3.8122")
    print("top-75 overlap with PROD:", {k: round(v["overlap_prod"], 3) for k, v in out["means"].items()})
    print("error-direction cos vs PROD (ORACLE):",
          {k: round(v["cos_err_vs_prod"], 3) for k, v in out["means"].items()},
          " random ref %.3f" % rand_ref)
    print()
    print(txt)


if __name__ == "__main__":
    sys.exit(main())
