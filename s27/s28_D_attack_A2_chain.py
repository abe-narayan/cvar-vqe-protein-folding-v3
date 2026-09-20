#!/usr/bin/env python
"""s27/s28_D_attack_A2_chain.py -- lane D: the standing checks on lane A2's built-chain arms.

Reads `s27/results/s28_A_chain_rows.jsonl` (one merged record per pdb, `arms: {arm: {rmsd_chain,
rmsd_cloud, rg, bond}}`; the LAST record per pdb is the complete one) and prints, `ST.fmt`
verbatim, on the built chain:

  1. step_e vs prod (production re-projected by the same script; S28-L23(c)), e in 0.1, 0.3, 1;
     whole set and FAIL18 / other-108;
  2. step_e vs the MEAN of the SAME two projected random draws (rand0_e, rand1_e), the registered
     control (PREREG_S28_A addendum 4; S28-L23(c)); the random draws' best-of-2 is printed as an
     order statistic and never as a control;
  3. the circuit family: circP vs prod (the nearest point), circ_e vs circP (one step in theta);
  4. e-grid pricing: `ST.best_of_k_within` over the three e for the step arm and for the circuit
     arm (S28-L23(b)), on the (targets x 3) matrix of chain differences vs prod / vs circP;
  5. the point-cloud twin of every contrast (rmsd_cloud) beside the chain one, so the projection's
     shrinkage of a harm or a gain is visible;
  6. geometry: mean virtual bond and Rg of every arm (the chain has ideal bonds by construction).

The falsifier (addendum 4): some e beats production beyond MDE with the fold CI, 5/5, AND beats
the random mean beyond MDE. ORACLE: every value is an RMSD to a native; nothing is chosen by it.
Writes `s27/results/s28_D_attack_A2_chain.json`.
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
ROWS = os.path.join(RESULTS, "s28_A_chain_rows.jsonl")
OUT = os.path.join(RESULTS, "s28_D_attack_A2_chain.json")
FAIL18 = set(I.FAIL18)
ES = ("0.1", "0.3", "1")


def compact(o):
    d = {k: v for k, v in o.items() if k != "concentration"}
    d["fmt"] = ST.fmt(o)
    return d


def contrast(a_map, b_map, label, strata=True):
    pdbs = sorted(p for p in a_map if p in b_map and np.isfinite(a_map[p]) and np.isfinite(b_map[p]))
    a = np.array([a_map[p] for p in pdbs]); b = np.array([b_map[p] for p in pdbs])
    folds = ST.pinned_folds(pdbs)
    out = {"n": len(pdbs)}
    r = ST.compare(a, b, folds, names=pdbs, label=f"{label} (n={len(pdbs)})")
    print(ST.fmt(r)); out["all"] = compact(r)
    if strata:
        for tag, m in (("FAIL18", np.array([p in FAIL18 for p in pdbs])),
                       ("other108", np.array([p not in FAIL18 for p in pdbs]))):
            if m.sum() < 5:
                continue
            rs = ST.compare(a[m], b[m], folds[m], names=[p for p, s in zip(pdbs, m) if s],
                            label=f"    {label}, {tag} (n={int(m.sum())})")
            print(ST.fmt(rs)); out[tag] = compact(rs)
    return out


def main(rows_path=ROWS, out_path=OUT):
    last = {}
    with open(rows_path, encoding="utf-8") as f:
        for l in f:
            if l.strip():
                r = json.loads(l); last[r["pdb"]] = r
    need = ["prod", "circP"] + [f"{k}_e{e}" for e in ES for k in ("step", "rand0", "rand1", "circ")]
    pdbs = sorted(p for p, r in last.items() if all(a in r["arms"] for a in need))
    print(f"{len(last)} targets in the file, {len(pdbs)} complete on all {len(need)} A2 arms")

    def col(arm, key="rmsd_chain"):
        return {p: float(last[p]["arms"][arm][key]) for p in pdbs}

    res = {"rows": os.path.relpath(rows_path, ROOT).replace(os.sep, "/"), "n_targets": len(pdbs),
           "chain": {}, "cloud": {}, "grid": {}, "geometry": {}}
    for basis, key in (("chain", "rmsd_chain"), ("cloud", "rmsd_cloud")):
        print(f"\n########## basis: {basis} ({key}) ##########")
        B = res[basis]
        prod = col("prod", key)
        for e in ES:
            step = col(f"step_e{e}", key)
            r0, r1 = col(f"rand0_e{e}", key), col(f"rand1_e{e}", key)
            rmean = {p: 0.5 * (r0[p] + r1[p]) for p in pdbs}
            rbest = {p: min(r0[p], r1[p]) for p in pdbs}
            print(f"\n=== 1. step e={e} vs prod ({basis}) ===")
            B[f"step_e{e}_vs_prod"] = contrast(step, prod, f"step e={e} vs prod ({basis})")
            print(f"\n=== 2. step e={e} vs the mean of the same two projected random draws ({basis}) ===")
            B[f"step_e{e}_vs_randmean"] = contrast(step, rmean, f"step e={e} vs rand mean-of-2 ({basis})")
            B[f"randmean_e{e}_vs_prod"] = contrast(rmean, prod, f"rand mean-of-2 e={e} vs prod ({basis})", strata=False)
            B[f"randbest_e{e}_mean"] = float(np.mean([rbest[p] for p in pdbs]))
            print(f"    [order statistic, not a control] rand best-of-2 e={e} mean {B[f'randbest_e{e}_mean']:.4f}"
                  f" vs rand mean-of-2 {np.mean([rmean[p] for p in pdbs]):.4f} vs prod {np.mean([prod[p] for p in pdbs]):.4f}")
        print(f"\n=== 3. the circuit family ({basis}) ===")
        circP = col("circP", key)
        B["circP_vs_prod"] = contrast(circP, prod, f"circP (nearest point) vs prod ({basis})", strata=False)
        for e in ES:
            B[f"circ_e{e}_vs_circP"] = contrast(col(f"circ_e{e}", key), circP, f"circuit one-step e={e} vs circP ({basis})", strata=False)
            B[f"circ_e{e}_vs_prod"] = contrast(col(f"circ_e{e}", key), prod, f"circuit one-step e={e} vs prod ({basis})", strata=False)
        print(f"\n=== 4. e-grid pricing ({basis}) ===")
        G = res["grid"].setdefault(basis, {})
        for fam, base in (("step", prod), ("circ", circP)):
            M = np.array([[col(f"{fam}_e{e}", key)[p] - base[p] for e in ES] for p in pdbs])
            bk = ST.best_of_k_within(M)
            G[fam] = {"es": list(ES), **bk, "column_means": {e: float(m) for e, m in zip(ES, M.mean(0))}}
            print(f"  {fam} over e {ES}: per-target min - mean {bk['observed_gain']:+.4f}  null {bk['null_across_targets']:+.4f}"
                  f"  share {bk['share_accounted']:.2f}  k_eff {bk['k_eff']:.2f}  split-half {bk['split_half']:+.4f}"
                  f" ({100*bk['split_half_frac']:.0f}%)  {bk['verdict']}")
            print("    column means:", ", ".join(f"e={e} {m:+.4f}" for e, m in zip(ES, M.mean(0))))
    print("\n=== 6. geometry of every arm (mean over targets; chain bonds are ideal by construction) ===")
    for arm in need:
        g = {k: float(np.mean([last[p]["arms"][arm][k] for p in pdbs])) for k in ("rg", "bond") if k in last[pdbs[0]]["arms"][arm]}
        g["rmsd_chain_mean"] = float(np.mean([last[p]["arms"][arm]["rmsd_chain"] for p in pdbs]))
        g["rmsd_cloud_mean"] = float(np.mean([last[p]["arms"][arm]["rmsd_cloud"] for p in pdbs]))
        res["geometry"][arm] = g
        print(f"  {arm:12s} " + "  ".join(f"{k} {v:.4f}" for k, v in g.items()))
    ST.save_atomic(out_path, res, complete_keys=["chain", "cloud", "grid", "geometry"], module_file=__file__)
    print("wrote", out_path)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default=ROWS)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    main(a.rows, a.out)
