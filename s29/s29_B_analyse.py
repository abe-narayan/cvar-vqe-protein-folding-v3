#!/usr/bin/env python
"""s29/s29_B_analyse.py -- statistics for S29 lane B (`s29/PREREG_S29_B.md`).

analyse1  measurement 1: the hop-only gradient variance of A, A_c and G at n = 4..9, this lane's
          INDEPENDENT rows beside lane T's (S29-L11) and S28's (`s28_B_train.json`), the fitted
          log2 slopes on both register ranges, the law check Var = r_stable / D^2, J*, and the
          verdicts on the pre-registered bright line B1 and on lane T's predictions T1 and T2.

analyse2  measurement 2: what the exact ground state of diag(E) - J M selects, ORACLE-scored,
          paired against each target's own DIS top-75 uniform average, with the grid priced as an
          order statistic, the PERM and SPEC controls, the pole-cancellation check (S29-L11
          prediction 3), the one-parameter collapse fit and the SIGN statistics (prereg A3).

Every RMSD here is ORACLE and is labelled so in the entry; none of it tunes anything.
"""
from __future__ import annotations

import json
import os
import sys
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s29 import s29_B_compat as C                                    # noqa: E402

T_GRAD = os.path.join(HERE, "results", "s29_T_grad_rows.jsonl")
S28_TRAIN = os.path.join(ROOT, "s27", "results", "s28_B_train.json")


def _load(path: str) -> List[Dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(l) for l in fh if l.strip()]


def _slope(ns, vals):
    ns, vals = np.asarray(ns, float), np.asarray(vals, float)
    ok = vals > 0
    if ok.sum() < 3:
        return float("nan")
    return float(np.polyfit(ns[ok], np.log2(vals[ok]), 1)[0])


# ------------------------------------------------------------------------------- measurement 1
def analyse1(out_path: str = None) -> Dict:
    rows = C.load_all(C.GRAD_ROWS)
    out = dict(kind="s29_B_grad", lane="B", sprint=29, ns=C.NS,
               note="property measurement; no native, no RMSD; contract rule 9 applies",
               n_rows=len(rows), n_targets=len(set(r["pdb"] for r in rows)),
               var_diag_n9=C.VAR_DIAG_N9, settings=dict(n_theta=C.N_THETA, seed=C.SEED,
                                                        init_sd=C.INIT_SD, layers=C.LAYERS))
    per = {}
    for mat in ("A", "A_c", "G"):
        cells, vs = {}, []
        for n in C.NS:
            rs = [r for r in rows if r["matrix"] == mat and r["n"] == n]
            if not rs:
                continue
            v = float(np.median([r["var_g0"] for r in rs]))
            cells[str(n)] = dict(
                var_g0_median=v,
                var_g0_min=float(min(r["var_g0"] for r in rs)),
                var_g0_max=float(max(r["var_g0"] for r in rs)),
                r_stable_median=float(np.median([r["r_stable"] for r in rs])),
                pred_var_median=float(np.median([r["pred_var"] for r in rs])),
                meas_over_pred=float(np.median([r["var_g0"] / r["pred_var"] for r in rs])),
                lam2_over_lam1_median=float(np.median([r["spec_lam2_over_lam1"] for r in rs])),
                perron_uniform_overlap2_median=float(
                    np.median([r["spec_perron_uniform_overlap2"] for r in rs])),
                rank_median=float(np.median([r["spec_rank"] for r in rs])),
                ratio_to_diag=float(C.VAR_DIAG_N9 / v),
                n_targets=len(rs))
            vs.append((n, v))
        ns_ = [x[0] for x in vs]
        vv = [x[1] for x in vs]
        rs9 = [r for r in rows if r["matrix"] == mat and r["n"] == 9]
        r_st9 = float(np.median([r["r_stable"] for r in rs9])) if rs9 else float("nan")
        per[mat] = dict(per_n=cells,
                        slope_4_9=_slope(ns_, vv),
                        slope_4_8=_slope([n for n in ns_ if n <= 8],
                                         [v for n, v in vs if n <= 8]),
                        r_stable_n9=r_st9,
                        J_star_n9=float(512.0 * np.sqrt(C.VAR_DIAG_N9 / r_st9)) if r_st9 > 0 else None,
                        fd_rel_max=float(max((r.get("fd_rel_max", 0.0) for r in rows
                                              if r["matrix"] == mat), default=0.0)))
    out["matrices"] = per

    # --- cross-check against lane T's independent rows, paired by (pdb, matrix, n)
    if os.path.exists(T_GRAD):
        t = {(r["pdb"], r["matrix"], r["n"]): r["var_g0"] for r in _load(T_GRAD)}
        diffs, pairs = [], 0
        for r in rows:
            k = (r["pdb"], r["matrix"], r["n"])
            if k in t and t[k] != 0:
                diffs.append(abs(r["var_g0"] - t[k]) / abs(t[k]))
                pairs += 1
        out["cross_check_lane_T"] = dict(
            n_pairs=pairs,
            max_rel_diff=float(max(diffs)) if diffs else None,
            median_rel_diff=float(np.median(diffs)) if diffs else None,
            note="same 120 draws at seed 1009; an exact match is expected and a difference is a bug")

    # --- cross-check A against S28's own hop-only row (the historical anchor)
    if os.path.exists(S28_TRAIN):
        s28 = json.load(open(S28_TRAIN, encoding="utf-8"))["summary"].get("hop_only|J1", {})
        pn = s28.get("per_n", {})
        out["cross_check_S28_hop_only_A"] = {
            str(n): dict(s28=pn.get(str(n), {}).get("var_g0_median"),
                         mine=per["A"]["per_n"].get(str(n), {}).get("var_g0_median"))
            for n in C.NS if str(n) in pn}

    # --- the pre-registered verdicts
    ac, g = per.get("A_c", {}), per.get("G", {})
    v9 = lambda m: m.get("per_n", {}).get("9", {}).get("var_g0_median", float("nan"))
    b1_clause_i = bool(ac.get("slope_4_9", -9) > -1.0 and g.get("slope_4_9", -9) > -1.0)
    b1_clause_ii = bool(v9(ac) >= C.VAR_DIAG_N9 / 30.0 and v9(g) >= C.VAR_DIAG_N9 / 30.0)
    out["verdicts"] = dict(
        B1_threshold_var_n9=float(C.VAR_DIAG_N9 / 30.0),
        B1_clause_i_slopes_shallower_than_minus1=b1_clause_i,
        B1_clause_ii_within_30x_of_diagonal=b1_clause_ii,
        B1="SUPPORTED" if (b1_clause_i and b1_clause_ii) else
           ("REFUTED" if not (b1_clause_i or b1_clause_ii) else "PARTIAL"),
        T1_falsified=bool(ac.get("slope_4_9", -9) > -1.3 or g.get("slope_4_9", -9) > -1.3
                          or v9(ac) > 3e-5 or v9(g) > 3e-5),
        T2_max_meas_over_pred_n_ge_7=float(max(
            [per[m]["per_n"][str(n)]["meas_over_pred"] for m in per for n in (7, 8, 9)
             if str(n) in per[m]["per_n"]], default=float("nan"))),
    )
    out["verdicts"]["T2_falsified"] = bool(out["verdicts"]["T2_max_meas_over_pred_n_ge_7"] > 3.0
                                           or out["verdicts"]["T2_max_meas_over_pred_n_ge_7"] < 1 / 3.0)
    from s24 import stats_lib as ST
    ST.save_atomic(out_path or C.GRAD_OUT, out)
    print(json.dumps({k: out[k] for k in ("n_rows", "n_targets", "verdicts",
                                          "cross_check_lane_T")}, indent=1))
    for m in ("A", "A_c", "G"):
        p = per[m]
        print(f"{m:5s} slope 4..9 {p['slope_4_9']:+.3f}  4..8 {p['slope_4_8']:+.3f}  "
              f"var(n=9) {p['per_n']['9']['var_g0_median']:.4e}  "
              f"r_stable(n=9) {p['r_stable_n9']:.3f}  ratio/diag {p['per_n']['9']['ratio_to_diag']:.0f}x  "
              f"J* {p['J_star_n9']:.1f}")
    return out


# ------------------------------------------------------------------------------- measurement 2
def _cells_of(rows: List[Dict]):
    return sorted({(r["matrix"], r["J"]) for r in rows}, key=lambda x: (x[0], x[1]))


def analyse2(targets: str = "12", out_path: str = None) -> Dict:
    from s24 import stats_lib as ST
    rows = C.load_all(C.GS_ROWS)
    pdbs = sorted({r["pdb"] for r in rows})
    if targets == "12":
        keep = set(C.picks_12())
        rows = [r for r in rows if r["pdb"] in keep]
        pdbs = sorted({r["pdb"] for r in rows})
    by = {}
    for r in rows:
        by.setdefault((r["matrix"], r["J"]), {})[r["pdb"]] = r
    cells = _cells_of(rows)
    folds = ST.pinned_folds(pdbs)
    prod = np.array([by[cells[0]][p]["rmsd_prod"] for p in pdbs])
    out = dict(kind="s29_B_gs", lane="B", sprint=29, n=len(pdbs), pdbs=pdbs,
               label="ORACLE DIAGNOSTIC -- every RMSD reads the native; nothing here is deployable",
               j_grid=list(C.J_GRID), matrices=list(C.MAT_ARMS),
               anchors=dict(prod_mean=float(prod.mean()), prod_se=float(prod.std(ddof=1) / np.sqrt(len(prod))),
                            s10_5_ladder_constant=2.954,
                            s10_5_note="the brief's 2.954 is the S10-5 ladder figure, a different "
                                       "pool era; the paired comparator here is each target's own "
                                       "DIS top-75 uniform average"),
               arms={}, contrasts={})

    def vec(cell, field):
        return np.array([by[cell][p][field[0]][field[1]] if isinstance(field, tuple)
                         else by[cell][p][field] for p in pdbs], float)

    for cell in cells:
        d = by[cell]
        if len(d) != len(pdbs):
            continue
        name = f"{cell[0]}|J{cell[1]:g}"
        rs = [d[p] for p in pdbs]
        arm = dict(n=len(pdbs))
        for R in ("R2", "R3", "R4"):
            v = np.array([r[R]["rmsd"] for r in rs], float)
            fin = np.isfinite(v)
            arm[f"mean_{R}"] = float(v[fin].mean()) if fin.any() else float("nan")
            arm[f"se_{R}"] = float(v[fin].std(ddof=1) / np.sqrt(fin.sum())) if fin.sum() > 1 else float("nan")
            arm[f"n_defined_{R}"] = int(fin.sum())
            arm[f"rg_{R}"] = float(np.nanmean([r[R]["rg"] for r in rs]))
            arm[f"bond_{R}"] = float(np.nanmean([r[R]["bond"] for r in rs]))
        for f in ("pr", "sign_coh", "gap", "jac75", "mass_top75", "mean_pair_rmsd_R3",
                  "pole_imbalance", "hop", "d_R3_prod", "frac_neg_w", "ess_w", "eta", "pc1_r2"):
            vals = np.array([r.get(f, np.nan) for r in rs], float)
            arm[f + "_median"] = float(np.nanmedian(vals)) if np.isfinite(vals).any() else float("nan")
        # S29-L11 prediction 3: does the p-readout emit the pool mean?
        arm["frac_R3_within_floor"] = float(np.mean([r["d_R3_prod"] <= C.FLOOR_A for r in rs]))
        # prereg A3: the sign
        sc = [r.get("sign_correct") for r in rs if r.get("sign_correct") is not None]
        arm["n_sign_defined"] = len(sc)
        arm["sign_correct_frac"] = float(np.mean(sc)) if sc else float("nan")
        bs = np.array([r.get("rmsd_pc1_bestsign", np.nan) for r in rs], float)
        arm["mean_pc1_bestsign"] = float(np.nanmean(bs))
        arm["n_pc1_defined"] = int(np.isfinite(bs).sum())
        out["arms"][name] = arm
        # the paired contrast that the gate reads
        for R in ("R3", "R4", "R2"):
            v = np.array([d[p][R]["rmsd"] for p in pdbs], float)
            if not np.isfinite(v).all():
                out["contrasts"][f"{name}|{R}"] = dict(skipped="undefined on some targets",
                                                       n_defined=int(np.isfinite(v).sum()))
                continue
            o = ST.compare(v, prod, folds=folds, names=pdbs,
                           label=f"gs {name} {R} - PROD (point cloud, ORACLE)")
            out["contrasts"][f"{name}|{R}"] = dict(
                effect=o["effect"], se=o["se"], mde=o["mde"], x_mde=o["effect_over_mde"],
                ci_fold=o["ci95_fold"], ci_iid=o["ci95_iid"], folds_same_sign=o["folds_same_sign"],
                W=o["n_better"], L=o["n_worse"], T=o["n_tied"], verdict=o["verdict"],
                mean_a=o["mean_a"], mean_b=o["mean_b"], fmt=ST.fmt(o))

    # --- the grid as an order statistic (prereg gate condition 2), REAL arms only, R3 and R4
    for R in ("R3", "R4"):
        real = [c for c in cells if "|" not in c[0] and c[0] in ("A_c", "G")]
        M = []
        ok = True
        for p in pdbs:
            row = []
            for c in real:
                v = by[c][p][R]["rmsd"]
                row.append(v if np.isfinite(v) else np.nan)
            M.append(row)
        M = np.array(M, float)
        if np.isfinite(M).all() and M.shape[1] > 1:
            d = M - prod[:, None]
            bk = ST.best_of_k_within(d)
            out[f"best_of_k_{R}"] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                                     for k, v in bk.items() if k != "per_draw"}
            out[f"best_of_k_{R}"]["cells"] = [f"{a}|J{b:g}" for a, b in real]
        else:
            out[f"best_of_k_{R}"] = dict(skipped="undefined cells", n_defined=int(np.isfinite(M).sum()),
                                         n_total=int(M.size))

    # --- the gate
    gate = dict(condition_1=[], condition_4=[], note="prereg section 5.5 + addendum A3")
    for cell in cells:
        if "|" in cell[0] or cell[0] == "A":
            continue
        name = f"{cell[0]}|J{cell[1]:g}"
        for R in ("R3", "R4"):
            c = out["contrasts"].get(f"{name}|{R}")
            if not c or "effect" not in c:
                continue
            if c["effect"] <= -0.7 * c["mde"]:
                gate["condition_1"].append(dict(cell=name, readout=R, effect=c["effect"],
                                                x_mde=c["x_mde"], pr=out["arms"][name]["pr_median"],
                                                sign_correct_frac=out["arms"][name]["sign_correct_frac"]))
    gate["opened"] = bool(gate["condition_1"])
    out["gate"] = gate
    ST.save_atomic(out_path or C.GS_OUT, out)
    print(json.dumps(dict(n=out["n"], anchors=out["anchors"], gate=gate,
                          best_of_k_R3=out.get("best_of_k_R3"),
                          best_of_k_R4=out.get("best_of_k_R4")), indent=1)[:4000])
    print("\narm                    R3      R4      R2    |  pr   coh  polelim  d(R3,prod)  R4def  sign%  bestsign")
    for name, a in out["arms"].items():
        print(f"{name:22s} {a['mean_R3']:.4f} {a['mean_R4'] if np.isfinite(a['mean_R4']) else float('nan'):7.4f} "
              f"{a['mean_R2']:.4f} | {a['pr_median']:6.1f} {a['sign_coh_median']:.3f} "
              f"{a['pole_imbalance_median']:+.3f} {a['d_R3_prod_median']:9.4f} "
              f"{a['n_defined_R4']:3d}/{out['n']:d} "
              f"{a['sign_correct_frac'] if np.isfinite(a['sign_correct_frac']) else float('nan'):.3f} "
              f"{a['mean_pc1_bestsign'] if np.isfinite(a['mean_pc1_bestsign']) else float('nan'):.4f}")
    print(f"\nPROD (DIS top-75 uniform, ORACLE) {out['anchors']['prod_mean']:.4f} "
          f"+- {out['anchors']['prod_se']:.4f} (SE), n = {out['n']}")
    return out


# ================================================= measurement 5b: the endpoint (prereg addendum 3)
def analyse_end(basis: str = "cloud", out_path: str = None) -> Dict:
    """Every pre-registered contrast of measurement 5b, on the point cloud or the built chain.

    basis = "cloud"  reads `rmsd_Ralpha` / `rmsd_tailset` from `s29_B_tta_end_rows*.jsonl`
    basis = "chain"  reads `rmsd_chain` from `s29_B_tta_chain_rows*.jsonl`, keyed by readout
    """
    from s24 import stats_lib as ST
    from s29 import s29_B_tta as TT
    if basis == "cloud":
        rows = C.load_all(TT.END_ROWS)
        fields = {"R_alpha": "rmsd_Ralpha", "tailset": "rmsd_tailset"}
    else:
        rows = C.load_all(TT.CHAIN_ROWS)
        fields = {"R_alpha": "rmsd_chain", "tailset": "rmsd_chain"}
    by = {}
    for r in rows:
        if basis == "chain" and r.get("readout") not in (None, "C_Ralpha", "C_tailset"):
            continue
        k = r["arm"] if basis == "cloud" else (r["arm"], r.get("readout"))
        by.setdefault(k, {})[r["pdb"]] = r
    pdbs = sorted(set(r["pdb"] for r in rows))
    folds = ST.pinned_folds(pdbs)
    out = dict(kind="s29_B_tta_end", lane="B", sprint=29, basis=basis, n=len(pdbs), pdbs=pdbs,
               label="ACHIEVABLE SELECTIONS, ORACLE-SCORED; nothing here tunes anything",
               lam_grid=list(TT.LAM_GRID), arms={}, contrasts={})

    def vec(arm, field):
        d = by[arm]
        return np.array([d[p][field] for p in pdbs], float)

    full = [a for a in by if len(by[a]) == len(pdbs)]
    for a in sorted(full, key=str):
        d = by[a]
        row = dict(n=len(pdbs))
        for nm, f in fields.items():
            if f in next(iter(d.values())):
                v = vec(a, f)
                row["mean_" + nm] = float(v.mean())
                row["se_" + nm] = float(v.std(ddof=1) / np.sqrt(len(v)))
        for f in ("m_tail", "pr", "jac75", "lam_w_ratio", "rg_Ralpha", "bond_Ralpha",
                  "tail_is_prefix", "pad_mass", "F", "cvar", "H_nats", "s_val", "secs"):
            vals = [d[p].get(f) for p in pdbs if d[p].get(f) is not None]
            if vals:
                row[f + "_mean"] = float(np.mean([float(x) for x in vals]))
        out["arms"][str(a)] = row

    def add_contrast(label, a, b, field):
        if a not in by or b not in by:
            return
        if len(by[a]) != len(pdbs) or len(by[b]) != len(pdbs):
            return
        va, vb = vec(a, field), vec(b, field)
        if not (np.isfinite(va).all() and np.isfinite(vb).all()):
            return
        o = ST.compare(va, vb, folds=folds, names=pdbs, label=label)
        out["contrasts"][label] = dict(
            effect=o["effect"], se=o["se"], mde=o["mde"], x_mde=o["effect_over_mde"],
            ci_fold=o["ci95_fold"], ci_iid=o["ci95_iid"], folds_same_sign=o["folds_same_sign"],
            W=o["n_better"], L=o["n_worse"], T=o["n_tied"], verdict=o["verdict"],
            mean_a=o["mean_a"], mean_b=o["mean_b"], fmt=ST.fmt(o))

    if basis == "cloud":
        for nm, f in fields.items():
            for lam in TT.LAM_GRID:
                for sd in (0, 1):
                    arm = "vqe|lam%g|s%d" % (lam, sd)
                    add_contrast("%s %s - production" % (arm, nm), arm, "production", f)
                    if lam != 0.0:
                        add_contrast("%s %s - its own lam=0 (same seed)" % (arm, nm),
                                     arm, "vqe|lam0|s%d" % sd, f)
                    add_contrast("%s %s - untrained best-of-16 (same lam)" % (arm, nm),
                                 arm, "untr16|lam%g" % lam, f)
                    add_contrast("%s %s - fixed profile M6" % (arm, nm), arm,
                                 "fixed_profile|M6", f)
            add_contrast("fixed profile M6 %s - production" % nm, "fixed_profile|M6",
                         "production", f)
            add_contrast("untrained best-of-16 (lam=0) %s - production" % nm, "untr16|lam0",
                         "production", f)
        # the lam grid as an order statistic, per seed and readout (prereg B2.6)
        for nm, f in fields.items():
            for sd in (0, 1):
                arms = ["vqe|lam%g|s%d" % (lam, sd) for lam in TT.LAM_GRID]
                if all(a in by and len(by[a]) == len(pdbs) for a in arms):
                    M = np.column_stack([vec(a, f) for a in arms])
                    base = vec("production", f)
                    bk = ST.best_of_k_within(M - base[:, None])
                    out["best_of_k_%s_s%d" % (nm, sd)] = {
                        kk: (float(vv) if isinstance(vv, (int, float, np.floating)) else vv)
                        for kk, vv in bk.items() if kk != "per_draw"}
                    out["best_of_k_%s_s%d" % (nm, sd)]["arms"] = arms
    else:
        for lam in TT.LAM_GRID:
            for sd in (0, 1):
                arm = "vqe|lam%g|s%d" % (lam, sd)
                for ro in ("C_Ralpha", "C_tailset"):
                    add_contrast("CHAIN %s [%s] - production" % (arm, ro), (arm, ro),
                                 ("production", ro), "rmsd_chain")
                    if lam != 0.0:
                        add_contrast("CHAIN %s [%s] - its own lam=0" % (arm, ro), (arm, ro),
                                     ("vqe|lam0|s%d" % sd, ro), "rmsd_chain")
    ST.save_atomic(out_path or os.path.join(HERE, "results", "s29_B_tta_end_%s.json" % basis), out)
    print("n =", out["n"], "basis", basis)
    print("%-24s %9s %9s %7s %7s %7s %7s" % ("arm", "R_alpha", "tailset", "m", "PR", "jac75",
                                             "prefix"))
    for a in sorted(out["arms"], key=str):
        r = out["arms"][a]
        print("%-24s %9.4f %9.4f %7.1f %7.1f %7.3f %7.2f" % (
            a, r.get("mean_R_alpha", float("nan")), r.get("mean_tailset", float("nan")),
            r.get("m_tail_mean", float("nan")), r.get("pr_mean", float("nan")),
            r.get("jac75_mean", float("nan")), r.get("tail_is_prefix_mean", float("nan"))))
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyse1", action="store_true")
    ap.add_argument("--analyse2", action="store_true")
    ap.add_argument("--analyse-end", action="store_true")
    ap.add_argument("--basis", default="cloud")
    ap.add_argument("--targets", default="12")
    a = ap.parse_args()
    if a.analyse1:
        analyse1()
    if a.analyse2:
        analyse2(a.targets)
    if a.analyse_end:
        analyse_end(a.basis)
