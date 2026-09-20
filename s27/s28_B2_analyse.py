#!/usr/bin/env python
"""s27/s28_B2_analyse.py -- statistics for the B2 ENDPOINT at the scope of prereg addendum 2.

Reads `s27/results/s28_B2_rows_k10.jsonl` (k = 10, J in {0, 3}, graphs REAL and PERM, VQE seeds
0 and 1 plus the exact ground state; point cloud). Registered cells (addendum 2): VQE R1 and
R3 at J = 3 REAL (both seeds), GS R3 at J = 3 REAL (DEGENERATE where the k = 10 graph has more
than one component; the count is on every row), PERM at the same cell, and J = 0.
Comparators: the SAME readout at J = 0 (F1), PRODUCTION (DIS top-75 uniform, `rmsd_dis75`;
S28-L2(a)), PERM (F3); both seeds; FAIL18 / 108. Every contrast is `s24.stats_lib.compare`
with its `ST.fmt` block verbatim. The built chain is run only if a registered cell reaches
0.7x MDE on the point cloud against production or J = 0.

Anchors: the J = 0 rows must equal `s28_B_rows.jsonl :: vqe|s*|NONE|J0` bit-for-bit (the kNN
graph is multiplied by zero at J = 0, so the run is `run_cvar_vqe` line for line).
Writes `s27/results/s28_B2_summary.json`.
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

from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402
from s27.s28_B_analyse import compact      # noqa: E402

ROWS = os.path.join(B.RESULTS, "s28_B2_rows_k10.jsonl")
OUT = os.path.join(B.RESULTS, "s28_B2_summary.json")
SCREEN_X = 0.7


def main():
    rows = B.load_rows(ROWS)
    by = {}
    for r in rows:
        by.setdefault(B.arm_of(r), {})[r["pdb"]] = r
    pdbs = sorted(by["vqe|s0|NONE|J0"])
    n = len(pdbs)
    folds = ST.pinned_folds(pdbs)
    fail = np.array([p in B.FAIL18 for p in pdbs])
    arms = sorted(a for a, d in by.items() if len(d) == n)
    prod = np.array([by["vqe|s0|NONE|J0"][p]["rmsd_dis75"] for p in pdbs])

    def vec(arm, R):
        return np.array([by[arm][p][R]["rmsd"] for p in pdbs])

    out = dict(n=n, pdbs=pdbs, arms={}, contrasts={}, anchors={}, screen=[], scope="prereg addendum 2: k 10, J {0, 3}, REAL/PERM")
    # ---- anchors: J = 0 bit-identical to the Gaussian run's J = 0 (same code path, A x 0)
    gauss = {}
    for r in B.load_rows(B.ROWS):
        if r["graph"] == "NONE" and r["J"] == 0.0 and r["source"] == "vqe":
            gauss[(r["seed"], r["pdb"])] = r
    for s in B.SEEDS:
        a = f"vqe|s{s}|NONE|J0"
        same = [by[a][p]["R1"]["rmsd"] == gauss[(s, p)]["R1"]["rmsd"] and by[a][p]["m"] == gauss[(s, p)]["m"]
                and by[a][p]["F"] == gauss[(s, p)]["F"] for p in pdbs if (s, p) in gauss]
        out["anchors"][f"J0_bit_identical_s{s}"] = dict(n_checked=len(same), n_identical=int(sum(same)))
    out["anchors"]["dis75_mean"] = float(prod.mean())
    # ---- per-arm means and native-free diagnostics (incl. the component count)
    for arm in arms:
        rs = [by[arm][p] for p in pdbs]
        d = dict(n=n, m_mean=float(np.mean([r["m"] for r in rs])), pr_median=float(np.median([r["pr"] for r in rs])),
                 jac75_mean=float(np.mean([r["jac75"] for r in rs])),
                 n_components_max=int(max(r.get("n_components", 1) for r in rs)),
                 n_targets_disconnected=int(sum(r.get("n_components", 1) > 1 for r in rs)),
                 degenerate_any=bool(any(r.get("degenerate", False) for r in rs)),
                 lam2_over_lam1_median=float(np.median([r.get("lam2_over_lam1", np.nan) for r in rs])),
                 perron_uniform_overlap_median=float(np.median([r.get("perron_uniform_overlap", np.nan) for r in rs])))
        for k in ("tv_vs_J0", "hop", "hop_abs", "sign_coh", "F"):
            vals = [r[k] for r in rs if r.get(k) is not None]
            d[f"{k}_mean"] = float(np.mean(vals)) if vals else None
        if arm.startswith("vqe") and "J3" in arm:
            coh = np.array([r["sign_coh"] for r in rs])
            d["n_coherent_over_0.5"] = int((coh > 0.5).sum()); d["n_coh_between_0.1_0.5"] = int(((coh > 0.1) & (coh <= 0.5)).sum())
        for R in ("R1", "R2", "R3"):
            v = vec(arm, R)
            d[f"mean_{R}"] = float(v.mean()); d[f"fail18_{R}"] = float(v[fail].mean()); d[f"other_{R}"] = float(v[~fail].mean())
        out["arms"][arm] = d
    C = out["contrasts"]

    def cmp(key, a, b, label):
        C[key] = compact(ST.compare(a, b, folds=folds, names=pdbs, label=label))
        C[key + "|108"] = compact(ST.compare(a[~fail], b[~fail], folds=folds[~fail], names=[p for p, f in zip(pdbs, fail) if not f], label=label + " (108 non-FAIL18)"))
        C[key + "|FAIL18"] = compact(ST.compare(a[fail], b[fail], folds=folds[fail], names=[p for p, f in zip(pdbs, fail) if f], label=label + " (FAIL18)"))
        return C[key]

    # ---- the registered cells
    for s in B.SEEDS:
        base = f"vqe|s{s}|NONE|J0"
        for g in ("REAL", "PERM"):
            arm = f"vqe|s{s}|{g}|J3"
            if arm not in out["arms"]:
                continue
            for R in ("R1", "R3"):
                c = cmp(f"F1|{arm}|{R}", vec(arm, R), vec(base, R), f"B2 F1 {arm} {R} - {base} {R} (point cloud, kNN k10)")
                p_ = cmp(f"PROD|{arm}|{R}", vec(arm, R), prod, f"B2 PROD {arm} {R} - DIS top-75 uniform (point cloud)")
                if g == "REAL" and (abs(c["x_mde"]) >= SCREEN_X or abs(p_["x_mde"]) >= SCREEN_X):
                    out["screen"].append(dict(arm=f"{arm}|{R}", x_mde_F1=c["x_mde"], x_mde_PROD=p_["x_mde"]))
            if g == "REAL":
                perm = f"vqe|s{s}|PERM|J3"
                if perm in out["arms"]:
                    for R in ("R1", "R3"):
                        cmp(f"CTRL|{arm}|PERM|{R}", vec(arm, R), vec(perm, R), f"B2 CTRL {arm} {R} - {perm} {R} (point cloud)")
        for R in ("R1", "R3"):
            cmp(f"PROD|{base}|{R}", vec(base, R), prod, f"B2 PROD {base} {R} - DIS top-75 uniform (point cloud)")
    garm = "gs|s-1|REAL|J3"
    if garm in out["arms"]:
        cmp(f"GSPROD|{garm}|R3", vec(garm, "R3"), prod, f"B2 GSPROD {garm} R3 - DIS top-75 uniform (point cloud)")
        cmp(f"GS|{garm}|R3", vec(garm, "R3"), vec("vqe|s0|NONE|J0", "R1"), f"B2 GS {garm} R3 - vqe|s0|NONE|J0 R1 (point cloud)")
        for s in B.SEEDS:
            cmp(f"F2|vqe|s{s}|REAL|J3|R3", vec(f"vqe|s{s}|REAL|J3", "R3"), vec(garm, "R3"), f"B2 F2 vqe|s{s}|REAL|J3 R3 - {garm} R3 (point cloud)")
    for key in [k for k in C if k.startswith("F1|vqe|s0|") and not k.endswith(("|108", "|FAIL18"))]:
        k1 = key.replace("|s0|", "|s1|")
        if k1 in C:
            C[key]["effect_seed1"] = C[k1]["effect"]
            C[key]["replicates_seed1"] = bool(C[key]["ci_fold"][0] <= C[k1]["effect"] <= C[key]["ci_fold"][1])
    ST.save_atomic(OUT, out, module_file=__file__)
    A = out["anchors"]
    print(f"n={n}  DIS top-75 {A['dis75_mean']:.6f}  J0 bit-identical to the Gaussian run: {A['J0_bit_identical_s0']} {A['J0_bit_identical_s1']}")
    print(f"\n{'arm':20s} {'R1':>7} {'R3':>7} {'m':>5} {'PR':>6} {'jac75':>6} {'hop':>6} {'coh':>5} {'ncoh':>4} {'comp':>4} {'F':>8}")
    for arm in arms:
        d = out["arms"][arm]
        f = lambda x, w=6, p=3: (f"{x:{w}.{p}f}" if x is not None else " " * w)
        print(f"{arm:20s} {d['mean_R1']:7.4f} {d['mean_R3']:7.4f} {d['m_mean']:5.1f} {d['pr_median']:6.0f} {d['jac75_mean']:6.3f} {f(d['hop_mean'])} {f(d['sign_coh_mean'],5,2)} {d.get('n_coherent_over_0.5', ''):>4} {d['n_targets_disconnected']:4d} {f(d['F_mean'],8,4)}")
    print("\ncontrasts (point cloud):")
    for key, c in C.items():
        if key.endswith(("|108", "|FAIL18")):
            continue
        k108 = C.get(key + "|108"); kf = C.get(key + "|FAIL18")
        print(f"  {key:30s} eff {c['effect']:+.4f} x{c['x_mde']:+.2f} fold [{c['ci_fold'][0]:+.4f},{c['ci_fold'][1]:+.4f}] {c['folds_same_sign']}/5 {c['W']}W/{c['L']}L/{c['T']}T "
              f"| 108 {k108['effect']:+.4f} x{k108['x_mde']:+.2f} | FAIL18 {kf['effect']:+.4f} x{kf['x_mde']:+.2f}  {c['verdict'].split(' [')[0]}"
              + (f"  rep={c.get('replicates_seed1')}" if "replicates_seed1" in c else ""))
    print("\nSCREEN (REAL cells at |effect| >= 0.7x MDE vs J0 or vs production; the chain runs only for these):")
    for s_ in out["screen"]:
        print("  ", s_)
    if not out["screen"]:
        print("   none: no built chain is run (prereg addendum 2)")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
