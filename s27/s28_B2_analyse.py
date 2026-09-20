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

    # ---- PREREG addendum 3: the coherence-class split, REGISTERED before the job. Class from
    # the trained state's own sign coherence: COHERENT > 0.5, INCOHERENT <= 0.1, MIXED between.
    # Class counts per seed and graph; the registered cells contrasted on each class separately
    # beside the all-126 contrast. Changes no verdict.
    def cls(arm):
        coh = np.array([by[arm][p]["sign_coh"] for p in pdbs])
        return np.where(coh > 0.5, "coherent", np.where(coh <= 0.1, "incoherent", "mixed"))

    out["classes"] = {}
    for s in B.SEEDS:
        for g in ("REAL", "PERM"):
            arm = f"vqe|s{s}|{g}|J3"
            if arm not in out["arms"]:
                continue
            c_ = cls(arm)
            rs = [by[arm][p] for p in pdbs]
            d = dict(n_coherent=int((c_ == "coherent").sum()), n_mixed=int((c_ == "mixed").sum()),
                     n_incoherent=int((c_ == "incoherent").sum()))
            for name in ("coherent", "incoherent", "mixed"):
                m_ = c_ == name
                if m_.sum():
                    d[f"{name}_hop_mean"] = float(np.mean([r["hop"] for r, k in zip(rs, m_) if k]))
                    d[f"{name}_bound_mean"] = float(np.mean([r["hop_abs"] for r, k in zip(rs, m_) if k]))
                    d[f"{name}_coh_mean"] = float(np.mean([r["sign_coh"] for r, k in zip(rs, m_) if k]))
                    d[f"{name}_F_mean"] = float(np.mean([r["F"] for r, k in zip(rs, m_) if k]))
                    d[f"{name}_m_mean"] = float(np.mean([r["m"] for r, k in zip(rs, m_) if k]))
            out["classes"][arm] = d
            base = f"vqe|s{s}|NONE|J0"
            for name in ("coherent", "incoherent"):
                m_ = c_ == name
                if m_.sum() < 8:
                    continue
                sub = [p for p, k in zip(pdbs, m_) if k]
                for R in ("R1", "R3"):
                    a_ = vec(arm, R)[m_]
                    C[f"CLASS|{arm}|{R}|{name}|vs_J0"] = compact(ST.compare(a_, vec(base, R)[m_], folds=folds[m_], names=sub,
                                                                              label=f"B2 CLASS {name} (n {m_.sum()}) {arm} {R} - {base} {R} (point cloud; registered addendum 3)"))
                    C[f"CLASS|{arm}|{R}|{name}|vs_PROD"] = compact(ST.compare(a_, prod[m_], folds=folds[m_], names=sub,
                                                                                label=f"B2 CLASS {name} (n {m_.sum()}) {arm} {R} - DIS top-75 uniform (point cloud; registered addendum 3)"))
    if all(f"vqe|s{s}|REAL|J3" in out["arms"] for s in B.SEEDS):
        c0, c1 = cls("vqe|s0|REAL|J3"), cls("vqe|s1|REAL|J3")
        out["classes"]["REAL_targets"] = dict(coherent_both=int(((c0 == "coherent") & (c1 == "coherent")).sum()),
                                              coherent_either=int(((c0 == "coherent") | (c1 == "coherent")).sum()),
                                              coherent_neither=int(((c0 != "coherent") & (c1 != "coherent")).sum()))

    # ---- the built chain (addendum 2: only for a cell at 0.7x MDE on the point cloud; paired
    # against production on the chain). Production = S27 `chain_rows.jsonl :: DIS`, re-projected
    # in this lane's process and identical on 126/126 (`s28_B_prodcheck.json`). The J = 0 twin
    # on the chain is the Gaussian run's `vqe|s*|NONE|J0|R3` row (`s28_B_chain_rows.jsonl`):
    # the J = 0 state is bit-identical between the two runs (anchor above), so it is reused.
    chain_path = os.path.join(B.RESULTS, "s28_B2_chain_rows_k10.jsonl")
    if os.path.exists(chain_path):
        byc = {}
        for r in B.load_rows(chain_path):
            byc.setdefault(r["arm"], {})[r["pdb"]] = r
        gch = {}
        for r in B.load_rows(B.CHAIN_ROWS):
            gch.setdefault(r["arm"], {})[r["pdb"]] = r
        prod_chain = {r["pdb"]: r["rmsd_chain"] for r in B.load_rows(os.path.join(B.RESULTS, "chain_rows.jsonl")) if r["config"] == "DIS"}
        pc = np.array([prod_chain[p] for p in pdbs])
        out["chain"] = dict(production_chain_mean=float(pc.mean()))
        for arm, d_ in byc.items():
            if not all(p in d_ for p in pdbs):
                out["chain"][arm] = dict(n=len(d_), partial=True)
                continue
            a_ = np.array([d_[p]["rmsd_chain"] for p in pdbs])
            k_ = np.array([d_[p]["rmsd_cloud"] for p in pdbs])
            src, seed, graph, J, R = arm.split("|")
            e = dict(n=n, mean_chain=float(a_.mean()), mean_cloud=float(k_.mean()),
                     fail18_chain=float(a_[fail].mean()), other_chain=float(a_[~fail].mean()))
            e["vs_production"] = compact(ST.compare(a_, pc, folds=folds, names=pdbs, label=f"B2 CHAIN {arm} - production DIS top-75 (built chain, S27 chain_rows DIS)"))
            e["vs_production_nonfail18"] = compact(ST.compare(a_[~fail], pc[~fail], folds=folds[~fail], names=[p for p, f in zip(pdbs, fail) if not f], label=f"B2 CHAIN {arm} - production (built chain, 108 non-FAIL18)"))
            e["vs_production_fail18"] = compact(ST.compare(a_[fail], pc[fail], folds=folds[fail], names=[p for p, f in zip(pdbs, fail) if f], label=f"B2 CHAIN {arm} - production (built chain, FAIL18)"))
            twin = f"{src}|{seed}|NONE|J0|{R}"
            if twin in gch and all(p in gch[twin] for p in pdbs):
                b_ = np.array([gch[twin][p]["rmsd_chain"] for p in pdbs])
                e["vs_same_readout_J0"] = compact(ST.compare(a_, b_, folds=folds, names=pdbs, label=f"B2 CHAIN F1 {arm} - {twin} (built chain; the J = 0 twin from s28_B_chain_rows.jsonl, bit-identical state)"))
            ref = f"vqe|{seed}|NONE|J0|R1"
            if ref in gch and all(p in gch[ref] for p in pdbs):
                b_ = np.array([gch[ref][p]["rmsd_chain"] for p in pdbs])
                e["vs_J0_R1"] = compact(ST.compare(a_, b_, folds=folds, names=pdbs, label=f"B2 CHAIN {arm} - {ref} (built chain)"))
            # the registered class split on the chain (addendum 3)
            c_ = cls(f"{src}|{seed}|{graph}|{J}")
            for name in ("coherent", "incoherent"):
                m_ = c_ == name
                if m_.sum() >= 8:
                    sub = [p for p, k in zip(pdbs, m_) if k]
                    e[f"class_{name}_vs_production"] = compact(ST.compare(a_[m_], pc[m_], folds=folds[m_], names=sub, label=f"B2 CHAIN CLASS {name} (n {m_.sum()}) {arm} - production (built chain; registered addendum 3)"))
            out["chain"][arm] = e
    ST.save_atomic(OUT, out, module_file=__file__)
    A = out["anchors"]
    print(f"n={n}  DIS top-75 {A['dis75_mean']:.6f}  J0 bit-identical to the Gaussian run: {A['J0_bit_identical_s0']} {A['J0_bit_identical_s1']}")
    print(f"\n{'arm':20s} {'R1':>7} {'R3':>7} {'m':>5} {'PR':>6} {'jac75':>6} {'hop':>6} {'coh':>5} {'ncoh':>4} {'comp':>4} {'F':>8}")
    for arm in arms:
        d = out["arms"][arm]
        f = lambda x, w=6, p=3: (f"{x:{w}.{p}f}" if x is not None else " " * w)
        print(f"{arm:20s} {d['mean_R1']:7.4f} {d['mean_R3']:7.4f} {d['m_mean']:5.1f} {d['pr_median']:6.0f} {d['jac75_mean']:6.3f} {f(d['hop_mean'])} {f(d['sign_coh_mean'],5,2)} {d.get('n_coherent_over_0.5', ''):>4} {d['n_targets_disconnected']:4d} {f(d['F_mean'],8,4)}")
    print("\ncoherence classes (addendum 3; class rule coherent > 0.5, incoherent <= 0.1, mixed between):")
    for arm, d in out["classes"].items():
        if arm == "REAL_targets":
            print(f"  REAL targets: coherent on both seeds {d['coherent_both']}, either {d['coherent_either']}, neither {d['coherent_neither']}")
            continue
        print(f"  {arm:18s} coherent {d['n_coherent']:3d}  mixed {d['n_mixed']:3d}  incoherent {d['n_incoherent']:3d}"
              + "".join(f" | {nm}: hop {d.get(nm + '_hop_mean', float('nan')):.3f} bound {d.get(nm + '_bound_mean', float('nan')):.3f} coh {d.get(nm + '_coh_mean', float('nan')):.3f} F {d.get(nm + '_F_mean', float('nan')):.3f} m {d.get(nm + '_m_mean', float('nan')):.1f}"
                        for nm in ("coherent", "incoherent", "mixed") if d.get(nm + "_hop_mean") is not None))
    print("\ncontrasts (point cloud):")
    for key, c in C.items():
        if key.endswith(("|108", "|FAIL18")):
            continue
        if key.startswith("CLASS|"):
            n_sub = c["fmt"].split("n=")[1].split()[0]
            print(f"  {key:52s} n {n_sub:>3} eff {c['effect']:+.4f} x{c['x_mde']:+.2f} MDE {c['mde']:.4f} fold [{c['ci_fold'][0]:+.4f},{c['ci_fold'][1]:+.4f}] {c['folds_same_sign']}/5  {c['verdict'].split(' [')[0]}")
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
    if "chain" in out:
        print(f"\nBUILT CHAIN (production chain mean {out['chain']['production_chain_mean']:.4f}):")
        for arm, e in out["chain"].items():
            if arm == "production_chain_mean":
                continue
            if e.get("partial"):
                print(f"  {arm:24s} PARTIAL {e['n']}/126")
                continue
            pr_ = e["vs_production"]; pn = e["vs_production_nonfail18"]; pf = e["vs_production_fail18"]
            sr = e.get("vs_same_readout_J0"); v = e.get("vs_J0_R1")
            print(f"  {arm:24s} chain {e['mean_chain']:.4f} cloud {e['mean_cloud']:.4f} fail18 {e['fail18_chain']:.3f} other {e['other_chain']:.3f}"
                  f" | vsPROD {pr_['effect']:+.4f} x{pr_['x_mde']:+.2f} [{pr_['ci_fold'][0]:+.4f},{pr_['ci_fold'][1]:+.4f}] {pr_['folds_same_sign']}/5 {pr_['verdict'].split(' [')[0]}"
                  f" (108 {pn['effect']:+.4f} x{pn['x_mde']:+.2f}; FAIL18 {pf['effect']:+.4f} x{pf['x_mde']:+.2f})"
                  + (f" | F1 same-readout J0 {sr['effect']:+.4f} x{sr['x_mde']:+.2f} [{sr['ci_fold'][0]:+.4f},{sr['ci_fold'][1]:+.4f}] {sr['folds_same_sign']}/5" if sr else "")
                  + (f" | vsJ0R1 {v['effect']:+.4f} x{v['x_mde']:+.2f}" if v else "")
                  + "".join(f" | class {nm} vsPROD {e['class_' + nm + '_vs_production']['effect']:+.4f} x{e['class_' + nm + '_vs_production']['x_mde']:+.2f}" for nm in ("coherent", "incoherent") if f"class_{nm}_vs_production" in e))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
