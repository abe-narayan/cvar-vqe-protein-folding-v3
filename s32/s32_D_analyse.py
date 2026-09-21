#!/usr/bin/env python
"""s32/s32_D_analyse.py -- S32 lane D aggregation for D3-M (mover) and D2-R (basin).

Registered in s32/PREREG_S32_D.md, commit 34973b1b.

BASIS.  Every RMSD is a BUILT CHAIN.  The production chain is rebuilt inside the D3-M job from
the shipped top-75 indices and the relaxed chain is produced from it in the SAME job, so the
pairing is within-job and contract rule 3 is satisfied without a bit-identity check.

SIGNS.  `stats_lib.compare` is LOWER IS BETTER.  For the ENDPOINT rows a negative effect is an
improvement and the verdict labels are read as printed.  For the COSINE rows the quantity is a
correlation, so the verdict labels would be inverted and are NOT quoted -- only effect, se, mde,
x_mde, the fold-clustered CI and folds_same_sign are read.

    python s32/s32_D_analyse.py
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s24 import stats_lib as ST          # noqa: E402

RES = os.path.join(HERE, "results")


def load(pat):
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, pat))):
        with open(f, encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if ln:
                    rows.append(json.loads(ln))
    rows.sort(key=lambda r: r["pdb"])
    return rows


def brief(c):
    return {"mean": c["effect"], "median": c["median_effect"], "se": c["se"], "mde": c["mde"],
            "x_mde": abs(c["effect"]) / c["mde"] if c["mde"] else float("nan"),
            "ci95_fold": c["ci95_fold"], "folds_same_sign": c["folds_same_sign"],
            "n_better": c["n_better"], "n_worse": c["n_worse"], "verdict": c["verdict"]}


def analyse_mover(out):
    rows = [r for r in load("s32_D3_mover_A_*.jsonl") if "k0_rmsd" in r]
    if not rows:
        return
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    n = len(rows)
    zero = np.zeros(n)
    prod = np.array([r["rmsd_prod"] for r in rows])
    out["D3M"] = {"n": n, "mean_rmsd_prod_built_chain": float(prod.mean()),
                  "median_rmsd_prod_built_chain": float(np.median(prod)),
                  "se_rmsd_prod": float(prod.std(ddof=1) / np.sqrt(n)),
                  "rungs": {}, "geometry": {}}
    for rung in ("k100", "k10", "k1", "k0"):
        if "%s_rmsd" % rung not in rows[0]:
            continue
        g = {}
        # --- (a) ORACLE / NOT DEPLOYABLE directional diagnostic
        for ref in ("eprod", "epool75", "edisto"):
            v = np.array([r["%s_cos_%s" % (rung, ref)] for r in rows], float)
            g["cos_%s" % ref] = brief(ST.compare(v, zero, folds=folds, names=pdbs,
                                                 label="%s cos vs -e_%s" % (rung, ref)))
            g["cos_%s" % ref]["ORACLE"] = "ORACLE / NOT DEPLOYABLE"
        cm = np.array([r["%s_ctrl_cos_mean" % rung] for r in rows], float)
        csd = np.array([r["%s_ctrl_cos_sd" % rung] for r in rows], float)
        cmx = np.array([r["%s_ctrl_cos_max" % rung] for r in rows], float)
        g["ctrl_magnitude_matched_random"] = {
            "draw_mean_over_targets": float(cm.mean()), "draw_sd_mean": float(csd.mean()),
            "draw_max_mean": float(cmx.mean()), "n_draws_per_target": 32,
            "note": "own distribution, not a best draw (contract rule 10)"}
        v = np.array([r["%s_cos_eprod" % rung] for r in rows], float)
        g["cos_eprod_minus_control"] = brief(ST.compare(v, cm, folds=folds, names=pdbs,
                                                        label="%s cos - matched random" % rung))
        # --- (b) NATIVE-FREE ENDPOINT, built chain, paired, same job
        rm = np.array([r["%s_rmsd" % rung] for r in rows], float)
        g["ENDPOINT_built_chain"] = brief(ST.compare(rm, prod, folds=folds, names=pdbs,
                                                     label="%s relaxed vs production chain" % rung))
        g["ENDPOINT_built_chain"]["mean_arm"] = float(rm.mean())
        g["ENDPOINT_built_chain"]["mean_prod"] = float(prod.mean())
        # --- the SCALE ladder.  Each alpha is its own registered arm; no per-target maximum.
        # BASIS: an interpolate of two valid chains, NOT reprojected; the virtual-bond mean and
        # sd travel with it so the deviation from a valid chain can be priced.
        g["SCALE_ladder"] = {}
        for al in (25, 50, 75):
            k = "%s_a%02d_rmsd" % (rung, al)
            if k not in rows[0]:
                continue
            v = np.array([r[k] for r in rows], float)
            e = brief(ST.compare(v, prod, folds=folds, names=pdbs,
                                 label="%s alpha=%.2f vs production" % (rung, al / 100)))
            e["mean_arm"] = float(v.mean())
            e["vb_mean"] = float(np.mean([r["%s_a%02d_vb_mean" % (rung, al)] for r in rows]))
            e["vb_sd"] = float(np.mean([r["%s_a%02d_vb_sd" % (rung, al)] for r in rows]))
            e["basis"] = "interpolated chain, NOT reprojected"
            g["SCALE_ladder"]["alpha_%.2f" % (al / 100)] = e
        # --- (c) the G1-CLOSED twin and the geometry secondaries (contract rule 15, row one)
        nd = np.array([r["%s_norm_dd" % rung] for r in rows], float)
        mv = np.array([r["%s_move_ca" % rung] for r in rows], float)
        g["G1_closed_twin_norm_dd"] = {"mean": float(nd.mean()), "median": float(np.median(nd))}
        g["move_ca_from_production"] = {"mean": float(mv.mean()), "median": float(np.median(mv)),
                                        "max": float(mv.max())}
        g["virtual_bond"] = {
            "mean": float(np.mean([r["%s_vb_mean" % rung] for r in rows])),
            "sd": float(np.mean([r["%s_vb_sd" % rung] for r in rows])),
            "prod_mean": float(np.mean([r["vb_mean_prod"] for r in rows])),
            "prod_sd": float(np.mean([r["vb_sd_prod"] for r in rows]))}
        g["energy"] = {"e0_median": float(np.median([r["%s_e0" % rung] for r in rows])),
                       "e1_median": float(np.median([r["%s_e1" % rung] for r in rows])),
                       "frac_converged": float(np.mean([r["%s_converged" % rung] for r in rows]))}
        out["D3M"]["rungs"][rung] = g
    out["D3M"]["reference_geometry"] = {
        "cos_eprod_epool75": float(np.mean([r["cos_eprod_epool75"] for r in rows])),
        "cos_eprod_edisto": float(np.mean([r["cos_eprod_edisto"] for r in rows])),
        "norm_e_prod": float(np.mean([r["norm_e_prod"] for r in rows])),
        "norm_e_pool75": float(np.mean([r["norm_e_pool75"] for r in rows])),
        "note": "ORACLE / NOT DEPLOYABLE -- reference error vectors, recomputed in this job"}


def analyse_basin(out):
    rows = [r for r in load("s32_D2_basin_*.jsonl") if "spread_in" in r]
    if not rows:
        return
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    ratio = g("spread_out") / g("spread_in")
    out["D2R"] = {
        "n_targets": len(rows), "k_per_target": int(rows[0]["k"]),
        "REGISTERED_horn_contraction_ratio_le_0.70": {
            "median_ratio": float(np.median(ratio)), "mean_ratio": float(ratio.mean()),
            "n_below_0.70": int((ratio <= 0.70).sum()), "fired": bool(np.median(ratio) <= 0.70)},
        "REGISTERED_horn_frozen_move_le_0.30A": {
            "median_move": float(np.median(g("move_median"))),
            "mean_move": float(g("move_mean").mean()),
            "fired": bool(np.median(g("move_median")) <= 0.30)},
        "spread_in_median": float(np.median(g("spread_in"))),
        "spread_out_median": float(np.median(g("spread_out"))),
        "min_pair_in_median": float(np.median(g("min_pair_in"))),
        "min_pair_out_median": float(np.median(g("min_pair_out"))),
        "rho_e_in_mean": float(g("rho_e_in").mean()),
        "rho_e_out_mean": float(g("rho_e_out").mean()),
        "rho_e_in_rgpart_mean": float(g("rho_e_in_rgpart").mean()),
        "rho_e_out_rgpart_mean": float(g("rho_e_out_rgpart").mean()),
        "rho_e_in_se": float(g("rho_e_in").std(ddof=1) / np.sqrt(len(rows))),
        "rho_e_out_se": float(g("rho_e_out").std(ddof=1) / np.sqrt(len(rows))),
        "rr_in_mean": float(g("rr_in_mean").mean()),
        "rr_out_mean": float(g("rr_out_mean").mean()),
        "rr_in_best_mean": float(g("rr_in_best").mean()),
        "rr_out_best_mean": float(g("rr_out_best").mean()),
        "rho_rrin_rrout_mean": float(g("rho_rrin_rrout").mean()),
        "wall_mean_s": float(g("wall_mean").mean()),
        "frac_converged": float(g("frac_converged").mean()),
        "note": ("ORACLE rr used to SCORE only.  rho columns are in-band Spearman over the first "
                 "40 members of the shipped top-75 band, per target, averaged over targets.  "
                 "n_targets is small: these are DIAGNOSTIC, reported with se, and no MDE verdict "
                 "is claimed on a 15-target subset.")}


def main():
    out = {"prereg_commit": "34973b1b",
           "note": ("Lane D.  RMSD rows are BUILT CHAIN, paired within the same job.  Cosine rows "
                    "are ORACLE / NOT DEPLOYABLE diagnostics.  stats_lib.compare is LOWER IS "
                    "BETTER: verdicts are quoted for the ENDPOINT rows and NOT quoted for the "
                    "cosine rows, where they would be inverted.")}
    analyse_mover(out)
    analyse_basin(out)
    p = os.path.join(RES, "s32_D_analyse.json")
    ST.save_atomic(p, out, module_file=__file__)
    print(json.dumps(out, indent=2, default=float))
    print("\nwrote", p)


if __name__ == "__main__":
    main()
