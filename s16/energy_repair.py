"""SPRINT 16, ENERGY -- THE REPAIR BENCHMARK.  The comparison worth keeping.

Both energies have a defensible role on ONE axis, and it is not ranking.  Legacy detects
stereochemical defects; AMBER removes them.  This module makes that comparison explicit
on IDENTICAL input structures.

THE INPUT, identical for every arm: the raw all-atom coordinate average of the shipped
top-75 window set.  It is genuinely broken -- ~26% bond contraction, ~1.4 heavy-atom
clashes below 2.0 A per target -- which is what makes it a repair benchmark rather than a
polishing benchmark.  It is a pure cache read with no RNG.

THE PANEL, scored on the SAME structure whose RMSD is quoted (`energy_lib.panel`):
Ramachandran favoured / allowed / outlier, heavy-atom clashes below 2.0 and 2.6 A,
minimum heavy separation, bond strain, angle strain, cis-peptide fraction, omega
deviation, chirality (CB improper, D-centre count).

THE TWO ROLES
  DETECTORS  Legacy's eleven components, evaluated by `core.energy.energy_components`
             on the INPUT's real coordinates (not on an ideal-geometry rebuild), scored
             against every defect by Spearman rho and by AUROC for "this target has the
             defect".  This is the genuine eleven-term potential; no surrogate.
  REPAIRERS  AMBER ff14SB/GBn2 restrained relaxation at k = 30, against the incumbent
             ideal-geometry projection as the CLASSICAL CONTROL repairer.  Both start
             from the same input.

COST, at the RIGHT number.  Not the 8-14 ms AMBER single point: the k = 30 arm is a full
restrained MINIMISATION, measured at 12.57 s mean on this box -- ~1500x a single point
and ~2.8x the 4.46 s projection.  Every earlier budget-matched claim in this project used
the single-point number and was therefore wrong by three orders of magnitude.  The wall
times are re-measured here from the ablation's own timings and reported beside the panel.

    python -m s16.energy_repair
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                      # noqa: E402
from s14.avgspace import top75_windows               # noqa: E402
from s15.phys_repl import averaged_backbone_from      # noqa: E402
from s16 import energy_lib as L                      # noqa: E402
from s16 import energy_ablate as A                   # noqa: E402

#: The defect axes.  `sign` says which direction is WORSE, so a detector that correlates
#: POSITIVELY with `sign * defect` is flagging the defect.
DEFECTS = {"n_clash_2A": +1, "n_clash_2p6A": +1, "min_heavy": -1,
           "rama_outlier": +1, "rama_favoured": -1, "bond_strain": +1,
           "angle_strain": +1, "cis_frac": +1, "omega_dev": +1,
           "n_D_centres": +1, "geom_rms_rel_dev": +1}


def detectors(cache=True):
    """Legacy's eleven components on the raw all-atom coordinate average, per target.

    Cheap: no AMBER, no projection.  ~0.5 s per target.  Cached, because the input is a
    deterministic cache read and recomputing it under box contention is pure waste.
    """
    import json
    from core import energy as et
    p = os.path.join(L.CACHE, "legacy_detectors.json")
    if cache and os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    rows = []
    for t in I.targets():
        pdb, seq = t["pdb"], t["seq"]
        W, PHI, PSI, u = top75_windows(pdb)
        avg, C_ca, _ = averaged_backbone_from(W, PHI, PSI)
        c = {a: np.asarray(avg[a], float) for a in L.ATOMS if a in avg}
        ph, ps = L.torsions_of(c)
        comp = et.energy_components(seq, c, ph, ps)
        tot = float(sum(float(et.DEFAULT_WEIGHTS[k]) * float(comp[k])
                        for k in L.LEG_TERMS))
        rows.append(dict(pdb=pdb, fold=int(t["fold"]), n=len(seq),
                         legacy_total_UNFITTED=tot,
                         **{"leg_" + k: float(comp[k]) for k in L.LEG_TERMS}))
        print(f"  detect {pdb} legacy {tot:+.2f} steric {comp['steric']:+.3f}", flush=True)
    if cache:
        with open(p, "w") as fh:
            json.dump(rows, fh, indent=1)
    return rows


def run():
    det = {r["pdb"]: r for r in detectors()}
    abl = {r["pdb"]: r for r in A.load_all()}
    pdbs = sorted(set(det) & set(abl))
    print(f"joined {len(pdbs)} targets", flush=True)

    per = []
    for p in pdbs:
        a = abl[p]
        row = dict(pdb=p, fold=a["fold"], n=a["n"], **{k: v for k, v in det[p].items()
                                                       if k.startswith("leg")})
        row["input"] = a["ctrl_rawavg_panel"]
        row["projected"] = a["ctrl"]["panel"]
        row["amber"] = a["ctrlamber"].get("panel")
        row["input_rmsd_ORACLE"] = a["ctrl_rawavg_rmsd"]
        row["projected_rmsd_ORACLE"] = a["ctrl"]["rmsd"]
        row["amber_rmsd_ORACLE"] = a["ctrlamber"].get("rmsd")
        row["amber_converged"] = a["ctrlamber"].get("converged")
        row["amber_energy"] = a["ctrlamber"].get("energy")
        row["wall_projection"] = a["ctrl"]["wall"]
        row["wall_amber"] = a["ctrlamber"].get("wall")
        per.append(row)

    fold = np.array([r["fold"] for r in per])
    names = [r["pdb"] for r in per]
    E = np.array([r["amber_energy"] if r["amber_energy"] is not None else np.inf
                  for r in per])

    # ---------------- REPAIR: the panel, on identical inputs
    repair = {}
    for key in DEFECTS:
        inp = np.array([r["input"].get(key, np.nan) for r in per], float)
        pro = np.array([r["projected"].get(key, np.nan) for r in per], float)
        amb = np.array([(r["amber"] or {}).get(key, np.nan) for r in per], float)
        ok = np.isfinite(inp) & np.isfinite(pro) & np.isfinite(amb)
        if ok.sum() < 5:
            continue
        repair[key] = {
            "n": int(ok.sum()),
            "input_mean": float(inp[ok].mean()), "input_median": float(np.median(inp[ok])),
            "projected_mean": float(pro[ok].mean()),
            "amber_mean": float(amb[ok].mean()),
            "AMBER_vs_input": L.gated_paired(amb[ok], inp[ok], E[ok],
                                             folds=fold[ok],
                                             names=[names[i] for i in np.where(ok)[0]]),
            "AMBER_vs_projection": L.gated_paired(amb[ok], pro[ok], E[ok],
                                                  folds=fold[ok],
                                                  names=[names[i] for i in np.where(ok)[0]]),
            "PROJECTION_vs_input": I.paired(pro[ok], inp[ok], folds=fold[ok],
                                            names=[names[i] for i in np.where(ok)[0]])}

    # ---------------- DETECT: how well does each Legacy component flag each defect?
    detect = {}
    for key, sgn in DEFECTS.items():
        y = np.array([r["input"].get(key, np.nan) for r in per], float)
        rem = np.array([(r["input"].get(key, np.nan)
                         - (r["amber"] or {}).get(key, np.nan)) for r in per], float)
        ok = np.isfinite(y)
        if ok.sum() < 5 or np.nanstd(y[ok]) == 0:
            continue
        # "has the defect" = worse than this instrument's own median on that axis
        thr = np.median(y[ok])
        lab = (sgn * y[ok]) > (sgn * thr)
        d = {}
        for term in L.LEG_TERMS + ("total_UNFITTED",):
            s = np.array([r["leg_" + term] if term != "total_UNFITTED"
                          else r["legacy_total_UNFITTED"] for r in per], float)[ok]
            if np.std(s) == 0:
                d[term] = {"rho": None, "auroc": None, "rho_removed": None}
                continue
            d[term] = {"rho": L.spearman(sgn * s, sgn * y[ok]) if False else
                              L.spearman(s, sgn * y[ok]),
                       "auroc": L.auroc(s, lab),
                       "rho_removed": (L.spearman(s[np.isfinite(rem[ok])],
                                                  (sgn * rem[ok])[np.isfinite(rem[ok])])
                                       if np.isfinite(rem[ok]).sum() > 4 else None)}
        detect[key] = {"n": int(ok.sum()), "sign": sgn, "median": float(thr),
                       "frac_labelled_defective": float(lab.mean()), "by_term": d}

    cost = {"projection_s": float(np.mean([r["wall_projection"] for r in per])),
            "amber_minimisation_s": float(np.mean([r["wall_amber"] for r in per
                                                   if r["wall_amber"]])),
            "amber_minimisation_median_s": float(np.median([r["wall_amber"] for r in per
                                                            if r["wall_amber"]])),
            "note": "measured under 2-shard contention; the s15 uncontended figures are "
                    "12.57 s (AMBER) and 4.46 s (projection); an AMBER SINGLE POINT is "
                    "8-14 ms, ~1500x cheaper, and is the WRONG number for this arm"}
    return {"per_target": per, "repair": repair, "detect": detect, "cost": cost,
            "n_targets": len(per)}


if __name__ == "__main__":
    t0 = time.time()
    out = run()
    print("\nwrote", L.write("energy_repair", out), f"({time.time()-t0:.0f}s)")
    print(f"\n{'defect':18s} {'input':>9s} {'projected':>10s} {'AMBER':>9s} "
          f"{'AMBER-proj':>11s} {'CI(gated)':>22s} {'W/L':>9s}")
    for k, v in out["repair"].items():
        g = v["AMBER_vs_projection"].get("gated") or v["AMBER_vs_projection"]["ungated"]
        print(f"{k:18s} {v['input_mean']:9.3f} {v['projected_mean']:10.3f} "
              f"{v['amber_mean']:9.3f} {g['mean_diff']:+11.4f} "
              f"[{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}] "
              f"{g['n_better']:4d}/{g['n_worse']:<4d}")
    print("\nBEST LEGACY DETECTOR per defect (AUROC on the INPUT structure):")
    for k, v in out["detect"].items():
        best = sorted(((d["auroc"], t) for t, d in v["by_term"].items()
                       if d["auroc"] is not None and np.isfinite(d["auroc"])),
                      key=lambda x: -abs(x[0] - 0.5))[:3]
        print(f"  {k:18s} " + "  ".join(f"{t}={a:.3f}" for a, t in best))
    print("\ncost:", {k: (round(v, 3) if isinstance(v, float) else v)
                      for k, v in out["cost"].items()})
