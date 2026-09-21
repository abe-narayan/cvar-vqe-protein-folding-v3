#!/usr/bin/env python
"""s32/s32_D2_basin.py -- S32 lane D rung D2-R: is basin physics OPERATIVE at 9-16 residues,
and does MINIMISATION make AMBER a ranker in band?

Registered in s32/PREREG_S32_D.md, commit 34973b1b.  Two questions, one set of rows.

(1) THE REALITY CHECK (prereg H-D2).  These are 9-16 residue peptides.  Conformational basins,
    metastability and folding cooperativity may not be well defined at that length, and either
    horn kills the dynamics hypothesis WITH A NAMED MECHANISM:
      CONTRACTION  median pairwise CA-RMSD among relaxed candidates <= 0.7x its pre-relaxation
                   value  ->  relaxation destroys the discrimination it was invoked to supply;
      FROZEN       median per-candidate CA displacement <= 0.3 A
                   ->  the trajectory carries no bit the starting structure did not carry.
    Both thresholds were fixed in the prereg before this file produced a number, so neither
    horn can be chosen after seeing the outcome.  The lane survives only if NEITHER fires.

(2) THE COMPLETION OF D1-E.  `s32_D1_inband.py` measured AMBER's in-band rho at +0.0000 using
    the cached SINGLE POINTS on the ideal-geometry rebuild, which are clash-dominated (e0 is
    routinely 1e3-1e8 kcal/mol).  That is not AMBER's best shot and the prereg does not let it
    stand as one.  The project's one positive physics-ranking result
    (memory `physics-ranks-real-geometry-not-lattice`, later overturned on other grounds) used
    CONVERGED `refine_coords`, not single points.  This file scores the RELAXED energy in band
    on the same band, with the Rg control, and reports both.

TARGET SELECTION IS NATIVE-FREE AND DETERMINISTIC: the first `PER_FOLD` targets of each of the
5 frozen folds, in pinned pdb order.  No native quantity enters the choice.
CANDIDATE SELECTION: the first `NCAND` members of the shipped top-75 band, in the shipped
order.  ORACLE `rr` is read only to SCORE, never to choose.

    python s32/s32_D2_basin.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I          # noqa: E402
from s32.s32_D1_inband import spear, spear_partial, rg   # noqa: E402

RES = os.path.join(HERE, "results")
PER_FOLD = 1
NCAND = 24
K_REST = 0.0          # free relaxation: the largest move physics can make
STEPS = 200
TOL = 5.0


def chosen_targets():
    tg = I.targets()
    out = []
    for f in sorted({int(t["fold"]) for t in tg}):
        out += [t for t in tg if int(t["fold"]) == f][:PER_FOLD]
    return out


def pair_rmsd_median(W):
    return float(np.median(I.pairwise_rmsd(W)[np.triu_indices(len(W), 1)]))


def one_target(t):
    from core import geometry as geo
    from core import amber as ar
    import torsion_lib2 as tl2

    u = I.load_univ(t["pdb"])
    p = I.pool_idx(u, k=500)
    sub = np.asarray(I.shipped_record(t["pdb"])["sub"], int)[:NCAND]
    W = np.asarray(u["W"][p], float)[sub]
    PHI = np.asarray(u["PHI"][p], float)[sub]
    PSI = np.asarray(u["PSI"][p], float)[sub]
    rr = np.asarray(u["rr"][p], float)[sub]                 # ORACLE label, scoring only
    nat = np.asarray(u["nat_ca"], float)

    tab = tl2.library_for(t["seq"], 8, t["seq"])
    rep = tl2.PerResidueTorsion(t["seq"], tab, chi_bits=False)
    BB = geo.build_backbone_batch(PHI, PSI)

    # The pre-relaxation object is the IDEAL-GEOMETRY REBUILD, not the deposited W, because that
    # is what AMBER is handed.  Basis stated: every RMSD below is on rebuilt chains.
    ca_in = np.asarray(BB["CA"], float)
    e_in, e_out, ca_out, moved, conv, wall = [], [], [], [], [], []
    for b in range(len(sub)):
        cd = {k: v[b] for k, v in BB.items()}
        t0 = time.time()
        try:
            r = ar.refine_coords(t["seq"], rep, cd, k_restraint=K_REST, steps=STEPS,
                                 tolerance=TOL, components=False, memo=False)
        except Exception as exc:                                  # noqa: BLE001
            return {"pdb": t["pdb"], "error": repr(exc)[:200]}
        wall.append(time.time() - t0)
        e_in.append(float(r["energy_initial"]))
        e_out.append(float(r["energy"]))
        ca_out.append(np.asarray(r["ca"], float))
        conv.append(bool(r["converged"]))
    ar.clear_cache()
    ca_out = np.stack(ca_out)
    e_in = np.array(e_in); e_out = np.array(e_out)
    moved = np.array([I.ca_rmsd(ca_out[b], ca_in[b]) for b in range(len(sub))])
    rr_in = np.array([I.ca_rmsd(ca_in[b], nat) for b in range(len(sub))])    # ORACLE
    rr_out = np.array([I.ca_rmsd(ca_out[b], nat) for b in range(len(sub))])  # ORACLE

    row = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "k": len(sub),
           "wall_mean": float(np.mean(wall)), "frac_converged": float(np.mean(conv)),
           # --- (1) the reality check
           "spread_in": pair_rmsd_median(ca_in),
           "spread_out": pair_rmsd_median(ca_out),
           "move_median": float(np.median(moved)),
           "move_mean": float(np.mean(moved)),
           "move_max": float(np.max(moved)),
           "min_pair_out": float(np.min(I.pairwise_rmsd(ca_out)[np.triu_indices(len(sub), 1)])),
           "min_pair_in": float(np.min(I.pairwise_rmsd(ca_in)[np.triu_indices(len(sub), 1)])),
           # --- energetics
           "e_in_median": float(np.median(e_in)), "e_out_median": float(np.median(e_out)),
           "e_out_sd": float(np.std(e_out)), "e_out_range": float(np.ptp(e_out)),
           # --- (2) does relaxation make AMBER a ranker IN BAND?
           "rho_e_in": spear(e_in, rr_in),
           "rho_e_out": spear(e_out, rr_out),
           "rho_e_out_vs_rrin": spear(e_out, rr_in),
           "rho_e_in_rgpart": spear_partial(e_in, rr_in, rg(ca_in)),
           "rho_e_out_rgpart": spear_partial(e_out, rr_out, rg(ca_out)),
           # --- does relaxation itself improve the candidates?  (ORACLE)
           "rr_in_mean": float(rr_in.mean()), "rr_out_mean": float(rr_out.mean()),
           "rr_in_best": float(rr_in.min()), "rr_out_best": float(rr_out.min()),
           "rho_rrin_rrout": spear(rr_in, rr_out),
           "rr_deposited_mean": float(rr.mean())}
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--per-fold", type=int, default=PER_FOLD)
    ap.add_argument("--ncand", type=int, default=NCAND)
    a = ap.parse_args()
    globals()["PER_FOLD"] = a.per_fold
    globals()["NCAND"] = a.ncand
    tg = chosen_targets()
    os.makedirs(RES, exist_ok=True)
    jl = os.path.join(RES, "s32_D2_basin_%d_%d.jsonl" % (a.shard, a.of))
    done = set()
    if os.path.exists(jl):
        with open(jl, encoding="utf-8") as fh:
            for ln in fh:
                try:
                    done.add(json.loads(ln)["pdb"])
                except Exception:                                 # noqa: BLE001
                    pass
    t0 = time.time()
    for i, t in enumerate(tg):
        if i % a.of != a.shard or t["pdb"] in done:
            continue
        row = one_target(t)
        with open(jl, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(row) + "\n")
        print("[D2 %d/%d] %s  elapsed %.0fs  %s" % (
            i + 1, len(tg), t["pdb"], time.time() - t0,
            json.dumps({k: round(v, 4) for k, v in row.items()
                        if k in ("spread_in", "spread_out", "move_median", "rho_e_in",
                                 "rho_e_out")})), flush=True)
    print("done", jl)


if __name__ == "__main__":
    main()
