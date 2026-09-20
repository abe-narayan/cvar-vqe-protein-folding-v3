#!/usr/bin/env python
"""s27/s28_D_attack_B_chain.py -- lane D: the standing checks on lane B's built-chain rows.

Reads `s27/results/s28_B_chain_rows.jsonl` (18 arms x targets, `arm|pdb|rmsd_chain|rmsd_cloud`)
and prints, with `ST.fmt` verbatim:

  1. every arm vs PRODUCTION on the built chain (S27 `chain_rows.jsonl :: DIS`, the reporting
     basis; S28-L2 caveat (a)), whole set and FAIL18 / other-108 (contract addendum 1, 14a).
     These contrasts CROSS a code path (no production arm is re-projected in B's chain job), so
     they carry the S28-L18/L27b floor (mean 0.006 A, tail 0.5 A on one target); said beside them.
  2. F1 on the chain: each of the nine 0.7x point-cloud cells vs the SAME readout at J = 0, same
     seed (R1 and R3 cells; the two R2 cells have no J = 0 R2 chain arm and are paired with
     J = 0 R1, as lane B's own analysis does), shared code path, whole set and strata.
  3. the nine cells priced as a best-of-nine: `ST.best_of_k_within` on the (targets x 9) matrix
     of chain differences vs the J = 0 R1 arm of the same seed (S28-L22), and on the same-readout
     matrix; the per-target minimum is an order statistic, never a result.
  4. the anchor: `rmsd_cloud` of `vqe|s{0,1}|NONE|J0|R1` vs `vqe_rows.jsonl :: DIS rmsd_vqe`.
  5. the largest point-cloud cell (s1 REAL J0.1 R3, S28-L21/L22) with its seed-0 twin and its
     PERM controls on the chain, same readout at J = 0 (a cell its PERM twin matches is the
     presence of a hopping term, not the pool graph's structure).
  GS R2 / R3 vs production are in 1 (W/L printed by `ST.fmt`, never quoted as evidence).

ORACLE: every value is an RMSD to a native (the endpoint); nothing is chosen by it.
Writes `s27/results/s28_D_attack_B_chain.json`.
"""
from __future__ import annotations

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
ROWS = os.path.join(RESULTS, "s28_B_chain_rows.jsonl")
OUT = os.path.join(RESULTS, "s28_D_attack_B_chain.json")
FAIL18 = set(I.FAIL18)

NINE = ["vqe|s0|REAL|J0.1|R1", "vqe|s0|REAL|J3|R2", "vqe|s0|PERM|J0.1|R1", "vqe|s0|PERM|J0.3|R3",
        "vqe|s0|PERM|J1|R3", "vqe|s0|RAND|J3|R2", "vqe|s1|REAL|J0.1|R3", "vqe|s1|PERM|J0.1|R3",
        "vqe|s1|PERM|J0.3|R3"]


def load_jsonl(path):
    with open(path, encoding="utf-8") as f:
        return [json.loads(l) for l in f if l.strip()]


def compact(o):
    d = {k: v for k, v in o.items() if k != "concentration"}
    d["fmt"] = ST.fmt(o)
    return d


def contrast(name, a_map, b_map, label):
    pdbs = sorted(p for p in a_map if p in b_map and np.isfinite(a_map[p]) and np.isfinite(b_map[p]))
    a = np.array([a_map[p] for p in pdbs]); b = np.array([b_map[p] for p in pdbs])
    folds = ST.pinned_folds(pdbs)
    out = {"n": len(pdbs)}
    r = ST.compare(a, b, folds, names=pdbs, label=f"{label} (n={len(pdbs)})")
    print(ST.fmt(r)); out["all"] = compact(r)
    for tag, m in (("FAIL18", np.array([p in FAIL18 for p in pdbs])),
                   ("other108", np.array([p not in FAIL18 for p in pdbs]))):
        if m.sum() < 5:
            continue
        rs = ST.compare(a[m], b[m], folds[m], names=[p for p, s in zip(pdbs, m) if s],
                        label=f"    {label}, {tag} (n={int(m.sum())})")
        print(ST.fmt(rs)); out[tag] = compact(rs)
    return out


def main(rows_path=ROWS, out_path=OUT):
    rows = load_jsonl(rows_path)
    T = {}
    for r in rows:
        T.setdefault(r["arm"], {})[r["pdb"]] = float(r["rmsd_chain"])
    arms = sorted(T)
    complete = sorted(set.intersection(*(set(T[a]) for a in arms)))
    print(f"{len(rows)} rows, {len(arms)} arms, {len(complete)} targets complete on every arm")
    prod = {r["pdb"]: float(r["rmsd_chain"]) for r in load_jsonl(os.path.join(RESULTS, "chain_rows.jsonl"))
            if r["config"] == "DIS"}
    res = {"rows": os.path.relpath(rows_path, ROOT).replace(os.sep, "/"), "n_targets": len(complete),
           "arms": arms, "vs_production": {}, "F1_chain": {}, "best_of_nine": {}, "anchor": {}}

    # 4. anchor
    vq = {}
    for r in load_jsonl(os.path.join(RESULTS, "vqe_rows.jsonl")):
        if r["config"] == "DIS":
            vq[(r["pdb"], int(r["seed"]))] = float(r["rmsd_vqe"])
    for s in (0, 1):
        cl = {r["pdb"]: float(r["rmsd_cloud"]) for r in rows if r["arm"] == f"vqe|s{s}|NONE|J0|R1"}
        d = [abs(cl[p] - vq[(p, s)]) for p in cl if (p, s) in vq]
        res["anchor"][f"s{s}"] = {"n": len(d), "max_abs_diff": float(max(d))}
        print(f"anchor seed {s}: J0 R1 rmsd_cloud vs S27 vqe_rows DIS: n {len(d)} max |diff| {max(d):.3g}")

    # 1. vs production
    print("\n=== 1. every arm vs PRODUCTION (S27 chain_rows :: DIS; crosses a code path, S28-L18/L27b floor) ===")
    for a in arms:
        res["vs_production"][a] = contrast(a, T[a], prod, f"{a} vs production (chain)")

    # 2. F1 on the chain, same readout at J = 0, same seed
    print("\n=== 2. F1 on the chain: the nine 0.7x cells vs the same readout at J = 0, same seed (shared code path) ===")
    diffs_same, diffs_r1 = {}, {}
    for a in NINE:
        src, seed, graph, J, R = a.split("|")
        comp_same = f"vqe|{seed}|NONE|J0|{R}"
        comp_r1 = f"vqe|{seed}|NONE|J0|R1"
        if comp_same not in T:
            print(f"  {a}: no {comp_same} chain arm; paired with {comp_r1} (as lane B does)")
            comp_same = comp_r1
        res["F1_chain"][a] = {"comparator": comp_same,
                              **contrast(a, T[a], T[comp_same], f"F1 {a} - {comp_same} (chain)")}
        diffs_same[a] = {p: T[a][p] - T[comp_same][p] for p in complete}
        diffs_r1[a] = {p: T[a][p] - T[comp_r1][p] for p in complete}

    # 3. best-of-nine pricing
    print("\n=== 3. the nine cells priced as a best-of-nine (ST.best_of_k_within) ===")
    for tag, D in (("vs_same_readout_J0", diffs_same), ("vs_J0_R1", diffs_r1)):
        M = np.array([[D[a][p] for a in NINE] for p in complete])
        bk = ST.best_of_k_within(M)
        res["best_of_nine"][tag] = {"arms": NINE, **bk}
        print(f"  [{tag}] per-target min - mean {bk['observed_gain']:+.4f}   across-target null {bk['null_across_targets']:+.4f}"
              f"   share accounted {bk['share_accounted']:.2f}   k_eff {bk['k_eff']:.2f}")
        print(f"    split-half transfer {bk['split_half']:+.4f} ({100*bk['split_half_frac']:.0f}% of oracle)   {bk['verdict']}")
        print("    argmin counts:", ", ".join(f"{a} ({c:.0f})" for a, c in zip(NINE, bk["argmin_counts"])))
        means = M.mean(0)
        print("    column means:", ", ".join(f"{a} {m:+.4f}" for a, m in zip(NINE, means)))
        res["best_of_nine"][tag]["column_means"] = {a: float(m) for a, m in zip(NINE, means)}

    # 5. the largest point-cloud cell (S28-L21/L22: s1 REAL J0.1 R3): its seed-0 twin and its PERM
    #    control on the chain, same readout at J = 0, shared code path; a cell that its PERM twin
    #    matches is the hopping term's presence, not the pool graph's structure.
    print("\n=== 5. the largest point-cloud cell's seed-0 twin and PERM control on the chain ===")
    res["largest_cell"] = {}
    for name, a, b in (("cell s1 REAL J0.1 R3", "vqe|s1|REAL|J0.1|R3", "vqe|s1|NONE|J0|R3"),
                       ("seed-0 twin", "vqe|s0|REAL|J0.1|R3", "vqe|s0|NONE|J0|R3"),
                       ("PERM control s1", "vqe|s1|PERM|J0.1|R3", "vqe|s1|NONE|J0|R3"),
                       ("PERM control s0 (J0.3)", "vqe|s0|PERM|J0.3|R3", "vqe|s0|NONE|J0|R3")):
        if a in T and b in T:
            res["largest_cell"][name] = {"a": a, "b": b, **contrast(name, T[a], T[b], f"{name}: {a} - {b} (chain)")}
    ST.save_atomic(out_path, res, complete_keys=["vs_production", "F1_chain", "best_of_nine", "anchor", "largest_cell"],
                   module_file=__file__)
    print("wrote", out_path)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default=ROWS)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    main(a.rows, a.out)
