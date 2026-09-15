#!/usr/bin/env python
"""s27/s28_B_analyse.py -- statistics for the S28-B run. Supersedes `s28_B_hop.py --analyse`.

Reads `s28_B_rows.jsonl` (point cloud) and, when present, `s28_B_chain_rows.jsonl` (built
chain) and S27's `chain_rows.jsonl` (the production built chain, config DIS). Every contrast is
`s24.stats_lib.compare` (paired, MDE, iid CI beside the fold CI, W/L, concentration null) and
its `ST.fmt` block is stored verbatim under `fmt`. Writes `s28_B_summary.json`.

Contrast families (`s27/PREREG_S28_B.md` section 6; S28-L2 caveats (a), (d), (e)):
  F1    VQE hopping arm (readout R, seed s, graph g, J) vs the SAME readout at J = 0, same seed
  F1b   the same arm vs the DEPLOYED readout at J = 0 (R1)
  PROD  the same arm vs PRODUCTION (DIS top-75 uniform; point cloud `rmsd_dis75`)
  CTRL  REAL vs PERM and REAL vs RAND at the same J, readout, seed
  F2    VQE vs the exact ground state at the same J and graph, R2 and R3 (R3-GS at J = 0 is
        degenerate and is never compared)
  GS    the eigensolver's readouts vs the deployed J = 0 R1
  GRID  the minimum over the four non-zero J rungs priced as an order statistic
  CHAIN the built chain: each projected arm vs J = 0 R1 (chain), vs production (S27 DIS chain),
        and on the 108 non-FAIL18 targets
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

S27_VQE = os.path.join(B.RESULTS, "vqe_rows.jsonl")
S27_CHAIN = os.path.join(B.RESULTS, "chain_rows.jsonl")
SCREEN_X = 0.7


def compact(o):
    return dict(effect=o["effect"], se=o["se"], mde=o["mde"], x_mde=o["effect_over_mde"],
                ci_iid=o["ci95_iid"], ci_fold=o["ci95_fold"], folds_same_sign=o.get("folds_same_sign"),
                W=o["n_better"], L=o["n_worse"], T=o["n_tied"], verdict=o["verdict"],
                mean_a=o["mean_a"], mean_b=o["mean_b"], median_effect=o["median_effect"],
                fmt=ST.fmt(o))


def main():
    rows = B.load_rows()
    by = {}
    for r in rows:
        by.setdefault(B.arm_of(r), {})[r["pdb"]] = r
    pdbs = sorted(by["vqe|s0|NONE|J0"])
    n = len(pdbs)
    folds = ST.pinned_folds(pdbs)
    fail = np.array([p in B.FAIL18 for p in pdbs])
    arms_full = [a for a, d in by.items() if len(d) == n]

    def vec(arm, R):
        return np.array([by[arm][p][R]["rmsd"] for p in pdbs])

    prod = np.array([by["vqe|s0|NONE|J0"][p]["rmsd_dis75"] for p in pdbs])
    out = dict(n=n, pdbs=pdbs, arms={}, contrasts={}, anchors={}, screen=[])

    # ---- anchors, bit-for-bit against S27
    out["anchors"]["dis75_mean"] = float(prod.mean())
    s27 = {}
    if os.path.exists(S27_VQE):
        for r in B.load_rows(S27_VQE):
            if r["config"] == "DIS":
                s27[(r["seed"], r["pdb"])] = r
    for s in B.SEEDS:
        a = f"vqe|s{s}|NONE|J0"
        out["anchors"][f"vqe_J0_R1_s{s}"] = float(vec(a, "R1").mean())
        if s27:
            same = [by[a][p]["R1"]["rmsd"] == s27[(s, p)]["rmsd_vqe"] and by[a][p]["m"] == s27[(s, p)]["m"]
                    for p in pdbs if (s, p) in s27]
            out["anchors"][f"s27_bit_identical_s{s}"] = dict(n_checked=len(same), n_identical=int(sum(same)),
                                                            s27_mean=float(np.mean([s27[(s, p)]["rmsd_vqe"] for p in pdbs if (s, p) in s27])))

    # ---- per-arm means and native-free diagnostics
    for arm in arms_full:
        rs = [by[arm][p] for p in pdbs]
        d = dict(n=n, m_mean=float(np.mean([r["m"] for r in rs])),
                 pr_median=float(np.median([r["pr"] for r in rs])),
                 jac75_mean=float(np.mean([r["jac75"] for r in rs])),
                 mass_top75_mean=float(np.mean([r["mass_top75"] for r in rs])),
                 entropy_bits_mean=float(np.mean([r["entropy_bits"] for r in rs])),
                 gate_pass_all=bool(all(r["gate_pass"] for r in rs)),
                 gate_equality_frac=float(np.mean([r["gate_equality"] for r in rs])),
                 holes_max=int(max(r["gate_n_holes"] for r in rs)),
                 secs_mean=float(np.mean([r["secs"] for r in rs])),
                 pad_mass_max=float(max(r["pad_mass"] for r in rs)),
                 degenerate=bool(rs[0].get("degenerate", False)))
        for R in ("R1", "R2", "R3"):
            v = vec(arm, R)
            d[f"mean_{R}"] = float(v.mean()); d[f"fail18_{R}"] = float(v[fail].mean()); d[f"other_{R}"] = float(v[~fail].mean())
        for k in ("tv_vs_J0", "hop", "hop_abs", "sign_coh", "F", "cvar", "e0", "gap"):
            vals = [r[k] for r in rs if r.get(k) is not None]
            d[f"{k}_mean"] = float(np.mean(vals)) if vals else None
        for R in ("R1", "R3"):
            d[f"{R}_mean_pair_rmsd"] = float(np.mean([r[R]["mean_pair_rmsd"] for r in rs]))
            d[f"{R}_overlap_dis75"] = float(np.mean([r[R]["overlap_dis75"] for r in rs]))
        out["arms"][arm] = d

    C = out["contrasts"]

    def cmp(key, a, b, label, mask=None):
        if mask is None:
            o = ST.compare(a, b, folds=folds, names=pdbs, label=label)
        else:
            o = ST.compare(a[mask], b[mask], folds=folds[mask], names=[p for p, m in zip(pdbs, mask) if m], label=label)
        C[key] = compact(o)
        return C[key]

    # ---- F1 / F1b / PROD / CTRL / F2
    for s in B.SEEDS:
        base = f"vqe|s{s}|NONE|J0"
        for g in B.GRAPHS:
            for J in B.J_GRID[1:]:
                arm = f"vqe|s{s}|{g}|J{J:g}"
                if arm not in out["arms"]:
                    continue
                for R in ("R1", "R2", "R3"):
                    c = cmp(f"F1|{arm}|{R}", vec(arm, R), vec(base, R), f"F1 {arm} {R} - {base} {R} (point cloud)")
                    if R != "R1":
                        cmp(f"F1b|{arm}|{R}", vec(arm, R), vec(base, "R1"), f"F1b {arm} {R} - {base} R1 (point cloud)")
                    cmp(f"PROD|{arm}|{R}", vec(arm, R), prod, f"PROD {arm} {R} - DIS top-75 uniform (point cloud)")
                    if abs(c["x_mde"]) >= SCREEN_X:
                        out["screen"].append(dict(arm=f"{arm}|{R}", x_mde=c["x_mde"], effect=c["effect"]))
                if g == "REAL":
                    for ctrl in ("PERM", "RAND"):
                        carm = f"vqe|s{s}|{ctrl}|J{J:g}"
                        if carm in out["arms"]:
                            for R in ("R1", "R2", "R3"):
                                cmp(f"CTRL|{arm}|{ctrl}|{R}", vec(arm, R), vec(carm, R), f"CTRL {arm} {R} - {carm} {R} (point cloud)")
                garm = f"gs|s-1|{g}|J{J:g}"
                if garm in out["arms"]:
                    for R in ("R2", "R3"):
                        cmp(f"F2|{arm}|{R}", vec(arm, R), vec(garm, R), f"F2 {arm} {R} - {garm} {R} (point cloud)")
    # ---- the eigensolver's readouts vs the deployed J = 0 R1 and vs production
    for g in B.GRAPHS:
        for J in B.J_GRID[1:]:
            garm = f"gs|s-1|{g}|J{J:g}"
            if garm in out["arms"]:
                for R in ("R1", "R2", "R3"):
                    cmp(f"GS|{garm}|{R}", vec(garm, R), vec("vqe|s0|NONE|J0", "R1"), f"GS {garm} {R} - vqe|s0|NONE|J0 R1 (point cloud)")
                    cmp(f"GSPROD|{garm}|{R}", vec(garm, R), prod, f"GSPROD {garm} {R} - DIS top-75 uniform (point cloud)")
    # ---- production anchors: J = 0 R1/R2/R3 vs DIS top-75
    for s in B.SEEDS:
        base = f"vqe|s{s}|NONE|J0"
        for R in ("R1", "R2", "R3"):
            cmp(f"PROD|{base}|{R}", vec(base, R), prod, f"PROD {base} {R} - DIS top-75 uniform (point cloud)")
    # ---- seed replication
    for key in [k for k in C if k.startswith("F1|vqe|s0|")]:
        k1 = key.replace("|s0|", "|s1|")
        if k1 in C:
            ci0 = C[key]["ci_fold"]; e1 = C[k1]["effect"]
            C[key]["effect_seed1"] = e1
            C[key]["replicates_seed1"] = bool(ci0[0] <= e1 <= ci0[1])
    # ---- the J grid as an order statistic (min over 4 rungs), per readout, seed, graph
    for s in B.SEEDS:
        base = f"vqe|s{s}|NONE|J0"
        for g in B.GRAPHS:
            for R in ("R1", "R2", "R3"):
                arms = [f"vqe|s{s}|{g}|J{J:g}" for J in B.J_GRID[1:]]
                if all(a in out["arms"] for a in arms):
                    M = np.stack([vec(a, R) - vec(base, R) for a in arms], 1)
                    C[f"GRID|{g}|s{s}|{R}"] = ST.best_of_k_within(M)

    # ---- built chain
    if os.path.exists(B.CHAIN_ROWS):
        cr = B.load_rows(B.CHAIN_ROWS)
        byc = {}
        for r in cr:
            byc.setdefault(r["arm"], {})[r["pdb"]] = r
        prod_chain = {}
        if os.path.exists(S27_CHAIN):
            for r in B.load_rows(S27_CHAIN):
                if r["config"] == "DIS":
                    prod_chain[r["pdb"]] = r["rmsd_chain"]
        out["chain"] = {}
        full = [a for a, d in byc.items() if all(p in d for p in pdbs)]
        ref = "vqe|s0|NONE|J0|R1"
        pc = np.array([prod_chain.get(p, np.nan) for p in pdbs])
        out["anchors"]["production_chain_mean"] = float(np.nanmean(pc)) if prod_chain else None
        for arm in full:
            a = np.array([byc[arm][p]["rmsd_chain"] for p in pdbs])
            k = np.array([byc[arm][p]["rmsd_cloud"] for p in pdbs])
            d = dict(n=n, mean_chain=float(a.mean()), mean_cloud=float(k.mean()),
                     fail18_chain=float(a[fail].mean()), other_chain=float(a[~fail].mean()))
            comp = ref.replace("|s0|", "|s1|") if "|s1|" in arm else ref
            if comp in byc and arm != comp:
                b = np.array([byc[comp][p]["rmsd_chain"] for p in pdbs])
                d["vs_J0_R1"] = compact(ST.compare(a, b, folds=folds, names=pdbs, label=f"CHAIN {arm} - {comp} (built chain)"))
                d["vs_J0_R1_nonfail18"] = compact(ST.compare(a[~fail], b[~fail], folds=folds[~fail],
                                                             names=[p for p, f in zip(pdbs, fail) if not f],
                                                             label=f"CHAIN {arm} - {comp} (built chain, 108 non-FAIL18)"))
                d["vs_J0_R1_fail18"] = compact(ST.compare(a[fail], b[fail], folds=folds[fail],
                                                          names=[p for p, f in zip(pdbs, fail) if f],
                                                          label=f"CHAIN {arm} - {comp} (built chain, FAIL18)"))
            if prod_chain and np.isfinite(pc).all():
                d["vs_production"] = compact(ST.compare(a, pc, folds=folds, names=pdbs, label=f"CHAIN {arm} - production DIS top-75 (built chain, S27 chain_rows DIS)"))
            out["chain"][arm] = d
    ST.save_atomic(B.SUMMARY, out, module_file=__file__)

    # ---- print
    A = out["anchors"]
    print(f"n={n}  DIS top-75 {A['dis75_mean']:.6f}  VQE J0 R1 s0 {A['vqe_J0_R1_s0']:.6f} s1 {A['vqe_J0_R1_s1']:.6f}  S27 bit-identical: {A.get('s27_bit_identical_s0')} {A.get('s27_bit_identical_s1')}")
    print(f"\n{'arm':22s} {'R1':>7} {'R2':>7} {'R3':>7} {'m':>5} {'PR':>7} {'jac75':>6} {'mass75':>6} {'Hbits':>6} {'TV':>6} {'hop':>6} {'coh':>5} {'eq':>4} {'F':>8}")
    for arm in sorted(out["arms"], key=lambda a: (a.split("|")[0], a.split("|")[2], float(a.split("|")[3][1:]), a.split("|")[1])):
        d = out["arms"][arm]
        f = lambda x, w=6, p=3: (f"{x:{w}.{p}f}" if x is not None else " " * w)
        print(f"{arm:22s} {d['mean_R1']:7.4f} {d['mean_R2']:7.4f} {d['mean_R3']:7.4f} {d['m_mean']:5.1f} {d['pr_median']:7.1f} "
              f"{d['jac75_mean']:6.3f} {d['mass_top75_mean']:6.3f} {d['entropy_bits_mean']:6.2f} {f(d['tv_vs_J0_mean'])} "
              f"{f(d['hop_mean'])} {f(d['sign_coh_mean'],5,2)} {d['gate_equality_frac']:4.2f} {f(d['F_mean'],8,4)}")
    print("\nF1 (VQE hopping arm vs its own J=0, same readout and seed), point cloud:")
    for key, c in C.items():
        if key.startswith("F1|"):
            print(f"  {key:30s} eff {c['effect']:+.4f} x{c['x_mde']:+.2f} fold [{c['ci_fold'][0]:+.4f},{c['ci_fold'][1]:+.4f}] "
                  f"{c['W']}W/{c['L']}L/{c['T']}T rep={c.get('replicates_seed1')}  {c['verdict']}")
    print("\nPROD (arm vs DIS top-75 uniform), point cloud:")
    for key, c in C.items():
        if key.startswith("PROD|"):
            print(f"  {key:30s} eff {c['effect']:+.4f} x{c['x_mde']:+.2f} fold [{c['ci_fold'][0]:+.4f},{c['ci_fold'][1]:+.4f}] {c['W']}W/{c['L']}L  {c['verdict']}")
    print("\nCTRL (REAL vs PERM / RAND), F2 (VQE vs GS), GS (eigensolver vs J0 R1):")
    for key, c in C.items():
        if key.startswith(("CTRL|", "F2|", "GS|", "GSPROD|")):
            print(f"  {key:36s} eff {c['effect']:+.4f} x{c['x_mde']:+.2f} fold [{c['ci_fold'][0]:+.4f},{c['ci_fold'][1]:+.4f}] {c['W']}W/{c['L']}L  {c['verdict']}")
    print("\nGRID (min over J priced as an order statistic):")
    for key, c in C.items():
        if key.startswith("GRID|"):
            print(f"  {key:22s} {json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in c.items() if k != 'per_target'})}")
    print("\nSCREEN (|effect| >= 0.7x MDE on F1, to the chain):")
    for s_ in out["screen"]:
        print(f"  {s_['arm']:34s} eff {s_['effect']:+.4f} x{s_['x_mde']:+.2f}")
    if "chain" in out:
        print(f"\nCHAIN (production chain mean {A.get('production_chain_mean')}):")
        for arm, d in out["chain"].items():
            v = d.get("vs_J0_R1"); pr_ = d.get("vs_production"); nf = d.get("vs_J0_R1_nonfail18")
            print(f"  {arm:30s} chain {d['mean_chain']:.4f} cloud {d['mean_cloud']:.4f} fail18 {d['fail18_chain']:.3f} other {d['other_chain']:.3f}"
                  + (f" | vsJ0 {v['effect']:+.4f} x{v['x_mde']:+.2f} [{v['ci_fold'][0]:+.4f},{v['ci_fold'][1]:+.4f}] {v['verdict']}" if v else "")
                  + (f" | non-FAIL18 {nf['effect']:+.4f} x{nf['x_mde']:+.2f}" if nf else "")
                  + (f" | vsPROD {pr_['effect']:+.4f} x{pr_['x_mde']:+.2f} {pr_['verdict']}" if pr_ else ""))
    print("wrote", B.SUMMARY)


if __name__ == "__main__":
    main()
