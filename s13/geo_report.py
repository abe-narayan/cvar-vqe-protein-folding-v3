"""SPRINT 13 GEO -- pooled tables and paired statistics over every geo_* result file.

Reads `s13/results/geo_{audit,anova,pauli,metric,metric_arm1,cvar,land,opt}.json` and emits
`s13/results/geo_report.json` plus a console digest.  Every headline gets n, mean, median,
paired difference, bootstrap CI and W/L where a pairing exists.

    python -m s13.geo_report
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402


def load(name):
    p = os.path.join(G.RESULTS, f"geo_{name}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def paired(a, b, n_boot=4000, seed=0):
    """b - a, paired.  Positive = b larger."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    d = b - a
    if d.size < 2:
        return {"n": int(d.size)}
    lo, hi = G.boot_ci(d, n_boot, seed)
    return {"n": int(d.size), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "d": float(d.mean()), "median_d": float(np.median(d)),
            "ci95": [lo, hi], "W": int((d < 0).sum()), "L": int((d > 0).sum()),
            "tie": int((d == 0).sum())}


def sec_pauli(out):
    d = load("pauli")
    if not d:
        return
    cells = d["cells"]
    key = lambda r: (r["pdb"], r["L"], r["k"], r["layers"])
    mw = {m: [r[m]["spectrum"]["mean_weight"] for r in cells] for m in G.VARIANTS}
    s2 = {m: [r[m]["spectrum"]["share_w_le2"] for r in cells] for m in G.VARIANTS}
    out["pauli"] = {
        "n_cells": len(cells),
        "mean_weight": {m: {"mean": float(np.mean(v)), "median": float(np.median(v)),
                            "min": float(np.min(v)), "max": float(np.max(v))}
                        for m, v in mw.items()},
        "share_w_le2": {m: {"mean": float(np.mean(v)), "median": float(np.median(v))}
                        for m, v in s2.items()},
        "paired_mean_weight_amber_minus_legacy": paired(mw["legacy"], mw["amber"]),
        "paired_share_le2_amber_minus_legacy": paired(s2["legacy"], s2["amber"]),
        "pred_ratio_median": {m: float(np.median(
            [r[m]["ratio_meas_over_pred_exact"] for r in cells])) for m in G.VARIANTS},
    }
    # gradient-variance scaling in n, at fixed depth, per model
    sc = {}
    for lay in sorted({r["layers"] for r in cells}):
        for m in G.VARIANTS:
            for a in ("a1.0", "a0.25", "a0.05"):
                pts = defaultdict(list)
                for r in cells:
                    if r["layers"] != lay:
                        continue
                    s = r[m]["scales"]["iqr"] or 1.0
                    pts[r["n_qubits"]].append(r[m]["var_grad_measured"][a] / s ** 2)
                if len(pts) >= 3:
                    ns = np.array(sorted(pts))
                    vs = np.array([np.mean(pts[n]) for n in ns])
                    sc[f"lay{lay}_{m}_{a}"] = {
                        "n": [int(x) for x in ns],
                        "var_iqr_normalised": [float(x) for x in vs],
                        **G.loglog_fit(ns, vs)}
    out["pauli"]["grad_var_scaling"] = sc


def sec_cvar(out):
    d = load("cvar")
    if not d:
        return
    rows = d["rows"]
    res = {}
    for m in sorted({r["model"] for r in rows}):
        for a in sorted({r["alpha"] for r in rows}):
            sel = [r for r in rows if r["model"] == m and r["alpha"] == a]
            if not sel:
                continue
            res[f"{m}|alpha={a}"] = {
                "n": len(sel),
                "grad_var_random_theta": float(np.mean(
                    [r["grad_var_random_theta"]["var_mean"] for r in sel])),
                "entropy_bits": float(np.mean([r["entropy_bits"] for r in sel])),
                "entropy_max_bits": float(np.mean([r["entropy_max_bits"] for r in sel])),
                "p_top": float(np.mean([r["p_top"] for r in sel])),
                "low_energy_tail_mass": float(np.mean(
                    [r["low_energy_tail_mass"] for r in sel])),
                "cvar_drop_frac_span": float(np.mean(
                    [r["cvar_drop"] / max(abs(r["cvar_start"]), 1e-12) for r in sel])),
                "effective_sample_size": float(np.mean(
                    [r["estimability"]["effective_sample_size"] for r in sel])),
                "cvar_shot_sd_relative": float(np.mean(
                    [r["estimability"]["relative_sd"] for r in sel])),
                "ORACLE_rmsd_argmax_p": float(np.mean(
                    [r["ORACLE_vqe"]["rmsd_argmax_p"] for r in sel])),
                "ORACLE_rmsd_best_in_top16": float(np.mean(
                    [r["ORACLE_vqe"]["rmsd_best_in_top16"] for r in sel])),
                "ORACLE_rmsd_p_weighted": float(np.mean(
                    [r["ORACLE_vqe"]["rmsd_p_weighted_top256"] for r in sel])),
                "ORACLE_gibbs_rmsd_p_weighted": float(np.mean(
                    [r["ORACLE_gibbs_matched_entropy"]["rmsd_p_weighted_top256"]
                     for r in sel])),
                "paired_vqe_minus_gibbs_pweighted": paired(
                    [r["ORACLE_gibbs_matched_entropy"]["rmsd_p_weighted_top256"]
                     for r in sel],
                    [r["ORACLE_vqe"]["rmsd_p_weighted_top256"] for r in sel]),
                "paired_vqe_minus_gibbs_top16": paired(
                    [r["ORACLE_gibbs_matched_entropy"]["rmsd_best_in_top16"] for r in sel],
                    [r["ORACLE_vqe"]["rmsd_best_in_top16"] for r in sel]),
            }
    out["cvar"] = res


def sec_opt(out):
    d = load("opt")
    if not d:
        return
    rows = d["rows"]
    res = {}
    for m in sorted({r["model"] for r in rows}):
        for b in sorted({r["budget"] for r in rows}):
            base = [r for r in rows if r["model"] == m and r["budget"] == b
                    and r["arm"] == "random"]
            bk = {(r["pdb"], r["L"], r["seed"]): r for r in base}
            for arm in sorted({r["arm"] for r in rows}):
                sel = [r for r in rows if r["model"] == m and r["budget"] == b
                       and r["arm"] == arm]
                if not sel:
                    continue
                pr = [bk[(r["pdb"], r["L"], r["seed"])] for r in sel
                      if (r["pdb"], r["L"], r["seed"]) in bk]
                res[f"{m}|budget={b}|{arm}"] = {
                    "n": len(sel),
                    "evals_spent_mean": float(np.mean(
                        [r["unique_evals_spent"] for r in sel])),
                    "best_energy_percentile": float(np.mean(
                        [r["best_energy_percentile"] for r in sel])),
                    "ORACLE_rmsd_of_best_seen": float(np.mean(
                        [r["ORACLE_rmsd_of_best_seen"] for r in sel])),
                    "metric_cond_mean": (float(np.mean(
                        [r["metric_cond_mean"] for r in sel]))
                        if sel[0].get("metric_cond_mean") else None),
                    "paired_vs_random_energy_pct": paired(
                        [r["best_energy_percentile"] for r in pr],
                        [r["best_energy_percentile"] for r in sel]) if pr else {},
                    "paired_vs_random_rmsd": paired(
                        [r["ORACLE_rmsd_of_best_seen"] for r in pr],
                        [r["ORACLE_rmsd_of_best_seen"] for r in sel]) if pr else {},
                }
    out["opt"] = res


def sec_land(out):
    d = load("land")
    if not d:
        return
    res = {}
    for m in G.VARIANTS:
        rows = [r for r in d["rows"] if m in r]
        if not rows:
            continue
        res[m] = {
            "n_cells": len(rows),
            "frac_negative_eigs_random_theta": float(np.mean(
                [np.mean([h["frac_negative"] for h in r[m]["hessian_random_theta"]])
                 for r in rows])),
            "frac_negative_eigs_at_endpoint": float(np.mean(
                [np.mean([h["frac_negative"] for h in r[m]["hessian_at_endpoint"]])
                 for r in rows])),
            "hessian_cond_at_endpoint": float(np.mean(
                [np.mean([h["cond_abs"] for h in r[m]["hessian_at_endpoint"]])
                 for r in rows])),
            "hessian_spectral_norm_at_endpoint": float(np.mean(
                [np.mean([h["spectral_norm"] for h in r[m]["hessian_at_endpoint"]])
                 for r in rows])),
            "n_distinct_minima_tol0.01": float(np.mean(
                [r[m]["multimodality"]["n_distinct_tol0.01"] for r in rows])),
            "frac_starts_within_1pct_of_best": float(np.mean(
                [r[m]["multimodality"]["frac_within_1pct_of_best"] for r in rows])),
            "final_value_sd_over_starts": float(np.mean(
                [r[m]["multimodality"]["val_sd"] for r in rows])),
        }
    res["geometry_vs_structure"] = {
        "corr_l2_vs_structural": float(np.mean(
            [r["geometry_vs_structure"]["corr_l2_vs_structural"] for r in d["rows"]])),
        "corr_fs_vs_structural": float(np.mean(
            [r["geometry_vs_structure"]["corr_fs_vs_structural"] for r in d["rows"]])),
        "corr_l2_vs_fs": float(np.mean(
            [r["geometry_vs_structure"]["corr_l2_vs_fs"] for r in d["rows"]])),
    }
    out["land"] = res


def main():
    out = {"what": "pooled s13 trainability-and-geometry report"}
    sec_pauli(out)
    sec_cvar(out)
    sec_opt(out)
    sec_land(out)
    G.write("geo_report", out)
    print(json.dumps(out, indent=1, default=float)[:12000])


if __name__ == "__main__":
    main()
