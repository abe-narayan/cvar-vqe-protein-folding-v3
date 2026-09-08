"""SPRINT 15, INFO, PART B -- resolve the per-target-skill / compactness discrepancy.

Sprint 14 left one number in two incompatible states:

    OBJ           per-target skill correlates +0.909 with the native Rg Z-SCORED
                  AGAINST THE ENUMERATED SPACE, over 12 targets, leave-one-target-out,
                  skill = `rho_global` on the FULL space.
    coordinator   the same chain's ORACLE ceiling measures +0.416 (perm p 0.177), over
                  the same 12 targets, skill = `cross_rho` from the CEILING experiment
                  (fit on the other targets' RMSD<=2.5 BANDS, scored on this target's
                  band), compactness LENGTH-RESIDUALISED rather than space-z-scored.

Five things differ between the two.  This module varies them ONE AT A TIME on a common
target set, so the discrepancy is decomposed rather than adjudicated:

    axis 1  training/eval region      full space   vs  RMSD<=2.5 band
    axis 2  skill statistic           rho_global   vs  rho(band)  vs  in-band pair acc
                                                   vs  -d_top100 (emitted RMSD)
    axis 3  compactness normalisation space z-score vs length-residual vs raw
    axis 4  correlation               Spearman     vs  Pearson
    axis 5  target set                12           vs  19

ORACLE DIAGNOSTIC throughout: every compactness variable reads the native trace.  This
measures a CEILING; nothing here is an inference-time procedure.

    python -m s15.info_sign
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s15 import info_lib as L            # noqa: E402

HEADLINE = os.path.join(ROOT, "s12", "results", "s14_obj_headline.json")
CEILING = os.path.join(ROOT, "s12", "results", "s14_obj_ceiling.json")


def load_skills():
    with open(HEADLINE) as fh:
        hl = json.load(fh)
    with open(CEILING) as fh:
        ce = json.load(fh)
    ceil = {r["pdb"]: r for r in ce["rows"]}
    return hl, ceil


def per_target():
    hl, ceil = load_skills()
    rows = []
    for p in L.enum_targets():
        g = L.geom_cache(p)
        n = int(g["n"])
        rg_nat = float(g["rg_native_ORACLE"])
        u_rg, u_r = np.asarray(g["u_rg"], float), np.asarray(g["u_rmsd"], float)
        b_rg, b_r = np.asarray(g["b_rg"], float), np.asarray(g["b_rmsd"], float)
        r = {
            "pdb": p, "n": n, "fold": int(g["fold"]),
            # ---- ORACLE compactness, three normalisations
            "rg_nat_ORACLE": rg_nat,
            "rg_space_mean": float(u_rg.mean()), "rg_space_sd": float(u_rg.std()),
            "z_space_ORACLE": float((rg_nat - u_rg.mean()) / u_rg.std()),
            # ---- model-free per-target ordering axes (the flip diagnostic, structural)
            "rho_space_rg_rmsd": L.spearman(u_rg, u_r),
            "rho_band_rg_rmsd": L.spearman(b_rg, b_r),
            "rho_space_c8_rmsd": L.spearman(np.asarray(g["u_c8"], float), u_r),
            # ---- reference objectives on the same uniform sample
            "rho_space_legacy": L.spearman(np.asarray(g["u_legacy"], float), u_r),
            "rho_space_prior": L.spearman(np.asarray(g["u_prior"], float), u_r),
            "band_n": int(len(b_r)), "space_mean_rmsd": float(u_r.mean()),
            "space_min_rmsd": float(g["rmsd_min_ORACLE"]),
        }
        a = hl["arms"]["learned"][p]
        r["hl_rho_global"] = a["rho_global"]
        r["hl_inband_lt15"] = a["inband_lt15"]
        r["hl_negd_top100"] = -a["d_top100"]
        r["hl_negd_decile"] = -a["d_decile"]
        for b in ("legacy", "prior"):
            r["base_%s_rho_global" % b] = hl["baselines"][b][p]["rho_global"]
            r["base_%s_inband" % b] = hl["baselines"][b][p]["inband_lt15"]
        if p in ceil:
            r["ceil_cross_rho"] = ceil[p]["cross_rho"]
            r["ceil_cross_pairacc"] = ceil[p]["cross_pairacc"]
            r["ceil_cross_negdt100"] = -ceil[p]["cross_dt100"]
            r["ceil_heldout_pairacc"] = ceil[p]["heldout_pairacc"]
        rows.append(r)
    return rows


SKILLS = [
    ("hl_rho_global", "LFO learned, rho_global on full space"),
    ("hl_inband_lt15", "LFO learned, in-band pair accuracy (<1.5 A)"),
    ("hl_negd_top100", "LFO learned, -d_top100 (emitted RMSD gain)"),
    ("rho_space_rg_rmsd", "model-free: rho(Rg, RMSD) over the space"),
    ("rho_band_rg_rmsd", "model-free: rho(Rg, RMSD) inside the band"),
    ("base_legacy_rho_global", "Legacy, rho_global"),
    ("ceil_cross_rho", "CEILING cross_rho (band-trained, band-scored)"),
    ("ceil_cross_pairacc", "CEILING cross-target in-band pair accuracy"),
    ("ceil_cross_negdt100", "CEILING -d_top100 inside the band"),
]

COMPACT = [
    ("z_space_ORACLE", "native Rg z-scored vs the enumerated space"),
    ("resid_n", "native Rg residualised on chain length"),
    ("rg_nat_ORACLE", "native Rg, raw"),
]


def grid(rows, subset=None, tag=""):
    rs = [r for r in rows if subset is None or r["pdb"] in subset]
    nvec = np.asarray([r["n"] for r in rs], float)
    rgn = np.asarray([r["rg_nat_ORACLE"] for r in rs], float)
    comp = {"z_space_ORACLE": np.asarray([r["z_space_ORACLE"] for r in rs], float),
            "resid_n": L.residualise(rgn, nvec),
            "rg_nat_ORACLE": rgn}
    out = {}
    for sk, sdesc in SKILLS:
        if not all(sk in r for r in rs):
            continue
        y = np.asarray([r[sk] for r in rs], float)
        for ck, cdesc in COMPACT:
            x = comp[ck]
            out[f"{sk}|{ck}"] = {
                "skill": sdesc, "compact": cdesc, "n": len(rs),
                "spearman": L.spearman(x, y), "pearson": L.pearson(x, y),
                "perm_p_spearman": L.perm_p(x, y, fn=L.spearman, n_perm=10000),
                "ci95_pearson": L.boot_ci(x, y)}
    return {"tag": tag, "targets": [r["pdb"] for r in rs], "cells": out}


def main():
    rows = per_target()
    all19 = [r["pdb"] for r in rows]
    twelve = sorted([r["pdb"] for r in rows if "ceil_cross_rho" in r])
    g19 = grid(rows, None, "all 19 enumerated targets")
    g12 = grid(rows, set(twelve), "the 12 targets of the ceiling experiment")
    out = {"per_target": rows, "grid_19": g19, "grid_12": g12,
           "targets_19": all19, "targets_12": twelve}
    L.jwrite("info_sign_grid", out)

    print("PART B.1 -- DECOMPOSING THE +0.909 vs +0.416 DISCREPANCY")
    print(f"19 enumerated targets: {' '.join(all19)}")
    print(f"12 ceiling targets:    {' '.join(twelve)}\n")
    for g in (g19, g12):
        print(f"--- {g['tag']} (n={len(g['targets'])}) ---")
        print(f"{'skill statistic':<46}{'compactness':<40}"
              f"{'spear':>8}{'pears':>8}{'permp':>8}")
        for kk, v in g["cells"].items():
            print(f"{v['skill']:<46}{v['compact']:<40}"
                  f"{v['spearman']:>8.3f}{v['pearson']:>8.3f}"
                  f"{v['perm_p_spearman']:>8.3f}")
        print()

    print("--- per-target table (ORACLE compactness vs the skill axes) ---")
    hdr = ("pdb", "n", "natRg", "spRg", "z", "rhoRg_sp", "rhoRg_bd", "hl_rhoG",
           "hl_inb", "ceil_rho")
    print("".join(f"{h:>10}" for h in hdr))
    for r in sorted(rows, key=lambda r: -r["z_space_ORACLE"]):
        print(f"{r['pdb']:>10}{r['n']:>10}{r['rg_nat_ORACLE']:>10.2f}"
              f"{r['rg_space_mean']:>10.2f}{r['z_space_ORACLE']:>10.2f}"
              f"{r['rho_space_rg_rmsd']:>10.3f}{r['rho_band_rg_rmsd']:>10.3f}"
              f"{r['hl_rho_global']:>10.3f}{r['hl_inband_lt15']:>10.3f}"
              + (f"{r['ceil_cross_rho']:>10.3f}" if "ceil_cross_rho" in r else f"{'-':>10}"))
    return out


if __name__ == "__main__":
    main()
