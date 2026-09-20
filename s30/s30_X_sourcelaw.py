"""S30 lane X, H-X1: does the SOURCE LAW predict the endpoint of every enlargement /
generation family in the record?

Pre-registration: `s30/PREREG_S30_X.md` section 2, committed at d4305d17 BEFORE this file
existed and before any number below was computed (contract rule 12, 14).

The law, coefficients FROZEN at the Sprint 20 wide-set values (2,142 (arm,target) cells,
17 generator arms at 8,192 draws each; project memory `operator-consumes-set-mean`):

    d_out  =  0.803 * d_set_mean  +  0.298 * d_set_best

Nothing here is fitted for the PRIMARY arm. The SECONDARY arm fits one slope and says so.

Reads only stored artefacts. No new pipeline compute, no box pressure.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A_MEAN, B_BEST = 0.803, 0.298          # FROZEN, Sprint 20 wide-set fit. Never refitted here.
OUT = os.path.join(BASE, "s30", "results", "s30_X_sourcelaw.json")

# ---------------------------------------------------------------- the K-ladder, published
# `docs/FINDINGS.md` section S7-6, `s7/poolsize.py`, `s7/poolsize_kcurve.json`, n=126.
# Transcribed from the FINDINGS table and verified in place before transcription.
# columns: K, selected endpoint, pool best (ORACLE), pool mean
KLADDER = [
    (25,   3.439, 2.350, 4.177),
    (50,   3.399, 2.177, 4.232),
    (100,  3.425, 1.970, 4.316),
    (250,  3.461, 1.808, 4.393),
    (500,  3.454, 1.711, 4.453),   # the shipped reference rung
    (1000, 3.483, 1.568, 4.522),
    (2000, 3.520, 1.504, 4.596),
]
REF_K = 500
DISCLOSED_CELL = (25, 2000)   # computed by hand before the prereg; excluded from the bar


def kladder():
    """PRIMARY. Frozen-coefficient prediction of every rung against the K=500 reference."""
    ref = [r for r in KLADDER if r[0] == REF_K][0]
    _, d_ref, best_ref, mean_ref = ref
    cells = []
    for K, d, best, mean in KLADDER:
        if K == REF_K:
            continue
        dmean, dbest = mean - mean_ref, best - best_ref
        pred = A_MEAN * dmean + B_BEST * dbest
        obs = d - d_ref
        cells.append({"K": K, "d_mean": dmean, "d_best": dbest,
                      "pred": pred, "obs": obs, "resid": pred - obs,
                      "sign_ok": bool(np.sign(pred) == np.sign(obs)) if obs != 0 else None})

    # the disclosed, NOT-out-of-sample cell, reported separately and excluded from the bar
    a = [r for r in KLADDER if r[0] == DISCLOSED_CELL[0]][0]
    b = [r for r in KLADDER if r[0] == DISCLOSED_CELL[1]][0]
    disc = {"cell": "K%d_to_K%d" % DISCLOSED_CELL,
            "pred": A_MEAN * (b[3] - a[3]) + B_BEST * (b[2] - a[2]),
            "obs": b[1] - a[1]}
    disc["resid"] = disc["pred"] - disc["obs"]

    res = np.array([c["resid"] for c in cells])
    signs = [c["sign_ok"] for c in cells if c["sign_ok"] is not None]

    # identifiability: across the ladder, mean and best are near-perfectly anti-correlated,
    # so this family CANNOT identify the two coefficients separately.  Registered in advance.
    means = np.array([r[3] for r in KLADDER])
    bests = np.array([r[2] for r in KLADDER])
    collin = float(np.corrcoef(means, bests)[0, 1])

    # the zero-information control in the operator's own space (contract rule 6, 7): the
    # BEST-ONLY law, i.e. the gate lane X was handed -- "the ceiling predicts the endpoint".
    best_only = [{"K": c["K"], "pred": c["d_best"], "obs": c["obs"],
                  "resid": c["d_best"] - c["obs"]} for c in cells]
    mean_only = [{"K": c["K"], "pred": c["d_mean"], "obs": c["obs"],
                  "resid": c["d_mean"] - c["obs"]} for c in cells]

    return {
        "cells": cells,
        "disclosed_cell_NOT_out_of_sample": disc,
        "n_cells": len(cells),
        "mean_abs_resid": float(np.abs(res).mean()),
        "max_abs_resid": float(np.abs(res).max()),
        "sign_correct": int(sum(signs)), "sign_n": len(signs),
        "FALSIFIER_mean_abs_resid_bar": 0.10,
        "FALSIFIER_sign_wrong_frac_bar": 1.0 / 3.0,
        "verdict": ("PASS" if (np.abs(res).mean() <= 0.10
                               and (len(signs) - sum(signs)) / max(len(signs), 1) < 1 / 3)
                    else "REFUTED"),
        "collinearity_mean_vs_best_r": collin,
        "control_best_only_mean_abs_resid": float(
            np.abs([c["resid"] for c in best_only]).mean()),
        "control_mean_only_mean_abs_resid": float(
            np.abs([c["resid"] for c in mean_only]).mean()),
        "control_best_only_sign_correct": int(sum(
            1 for c in best_only if np.sign(c["pred"]) == np.sign(c["obs"]))),
    }


# ------------------------------------------------- S24 lane C: 504 (arm, target) cells
ARMS = ("T0_helix", "T1_blind", "T2_restype", "T3_pool")


def s24_cells():
    """SECONDARY. The law's dominant claim -- the endpoint tracks the SELECTED SET's mean --
    on 126x4 cells from a different sprint, a different corpus and different samplers.

    One slope IS fitted here (the S20 coefficients were fitted on other data); the arm is
    labelled as a fit, not a prediction."""
    p = os.path.join(BASE, "s24", "results", "c_ladder.json")
    rows = json.load(open(p))["rows"]
    x, y, fold, tag = [], [], [], []
    union = []
    for r in rows:
        for a in ARMS:
            d = r[a]
            x.append(d["sel_member_mean_ORACLE"])
            y.append(d["rmsd"])
            fold.append(int(r["fold"]))
            tag.append(a)
            # the UNION ceiling S24 never formed: min(pool best, generated best)
            union.append({
                "pdb": r["pdb"], "arm": a, "fold": int(r["fold"]),
                "pool_best": r["pool_best_ORACLE"],
                "union_best_matched": min(r["pool_best_ORACLE"],
                                          d["gen_best_ORACLE_matched"]),
                "union_best_full": min(r["pool_best_ORACLE"], d["gen_best_ORACLE_full"]),
                "union_endpoint": d["union"], "pool_endpoint": r["incumbent"],
                "beats": int(d["gen_best_ORACLE_matched"] < r["pool_best_ORACLE"]),
            })
    x, y, fold = np.array(x), np.array(y), np.array(fold)
    # WITHIN-TARGET identification (contract rule 6: the control lives in the operator's own
    # space).  A pooled cross-target slope is contaminated by "hard targets are hard"; the
    # operator's response is the slope AFTER removing each target's mean.  4 arms/target.
    pdbs = np.array([r["pdb"] for r in rows for _ in ARMS])
    xw, yw = x.astype(float).copy(), y.astype(float).copy()
    for t in set(pdbs.tolist()):
        m = pdbs == t
        xw[m] -= xw[m].mean()
        yw[m] -= yw[m].mean()
    wslope = float(np.polyfit(xw, yw, 1)[0])
    wpred = wslope * xw
    wr2 = float(1 - ((yw - wpred) ** 2).sum() / (yw ** 2).sum())
    rngw = np.random.default_rng(1)
    tl = sorted(set(pdbs.tolist()))
    wb = []
    for _ in range(2000):
        pk = rngw.choice(tl, len(tl), replace=True)
        ii = np.concatenate([np.where(pdbs == t)[0] for t in pk])
        wb.append(np.polyfit(xw[ii], yw[ii], 1)[0])
    wlo, whi = np.percentile(wb, [2.5, 97.5])

    slope, icpt = np.polyfit(x, y, 1)
    pred = slope * x + icpt
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()

    # fold-clustered bootstrap on the slope (contract rule 4)
    rng = np.random.default_rng(0)
    folds = sorted(set(fold.tolist()))
    boots = []
    for _ in range(2000):
        pick = rng.choice(folds, len(folds), replace=True)
        idx = np.concatenate([np.where(fold == f)[0] for f in pick])
        if len(set(x[idx])) < 2:
            continue
        boots.append(np.polyfit(x[idx], y[idx], 1)[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])

    per_arm = {}
    for a in ARMS:
        m = np.array([t == a for t in tag])
        s, ic = np.polyfit(x[m], y[m], 1)
        pr = s * x[m] + ic
        per_arm[a] = {"slope": float(s), "intercept": float(ic),
                      "r2": float(1 - ((y[m] - pr) ** 2).sum()
                                  / ((y[m] - y[m].mean()) ** 2).sum()),
                      "set_mean": float(x[m].mean()), "endpoint": float(y[m].mean())}

    # the union ceiling, aggregated
    ub = {}
    for a in ARMS:
        rs = [u for u in union if u["arm"] == a]
        ub[a] = {
            "n_targets_with_a_better_structure_than_the_whole_pool_ORACLE":
                int(sum(u["beats"] for u in rs)),
            "pool_best_ORACLE": float(np.mean([u["pool_best"] for u in rs])),
            "union_best_matched_ORACLE": float(np.mean([u["union_best_matched"] for u in rs])),
            "union_best_full_ORACLE": float(np.mean([u["union_best_full"] for u in rs])),
            "d_ceiling_matched": float(np.mean([u["union_best_matched"] - u["pool_best"]
                                                for u in rs])),
            "d_endpoint": float(np.mean([u["union_endpoint"] - u["pool_endpoint"] for u in rs])),
        }
        # what Delta(set mean) the law REQUIRES to explain this endpoint, given the ceiling
        ub[a]["implied_d_set_mean"] = float(
            (ub[a]["d_endpoint"] - B_BEST * ub[a]["d_ceiling_matched"]) / A_MEAN)
    return {"n_cells": int(len(x)), "slope_pooled": float(slope), "intercept": float(icpt),
            "r2_pooled": float(r2), "fold_CI95_slope_pooled": [float(lo), float(hi)],
            "n_boot": len(boots), "per_arm": per_arm, "union_ceiling": ub,
            "WITHIN_TARGET_slope": wslope, "WITHIN_TARGET_r2": wr2,
            "WITHIN_TARGET_CI95_target_clustered": [float(wlo), float(whi)],
            "WITHIN_TARGET_note": ("target-demeaned, 4 arms per target: this is the "
                                   "operator's own response, free of cross-target difficulty")}


def gate_scorecard(kl, s24):
    """The head-to-head the lane exists for: the CEILING gate lane X was handed against the
    SOURCE LAW, scored on sign over every cell in the record where both are computable."""
    cells = []
    for c in kl["cells"]:
        cells.append({"family": "S7-6 K-ladder", "cell": "K%d" % c["K"],
                      "obs": c["obs"], "law_pred": c["pred"],
                      "ceiling_gate_pred": c["d_best"]})
    for a, u in s24["union_ceiling"].items():
        cells.append({"family": "S24 generated-source union", "cell": a,
                      "obs": u["d_endpoint"],
                      "law_pred": None,          # set mean not stored; sign not predictable
                      "ceiling_gate_pred": u["d_ceiling_matched"]})

    def sgn_ok(p, o):
        return None if (p is None or o == 0) else bool(np.sign(p) == np.sign(o))

    law = [sgn_ok(c["law_pred"], c["obs"]) for c in cells]
    gate = [sgn_ok(c["ceiling_gate_pred"], c["obs"]) for c in cells]
    return {
        "cells": cells,
        "ceiling_gate_sign_correct": int(sum(1 for g in gate if g)),
        "ceiling_gate_sign_n": int(sum(1 for g in gate if g is not None)),
        "source_law_sign_correct": int(sum(1 for g in law if g)),
        "source_law_sign_n": int(sum(1 for g in law if g is not None)),
        "READING": ("the ceiling gate predicts the DIRECTION of the endpoint move; a gate "
                    "that is wrong on the sign is worse than no gate, because it rejects "
                    "exactly the spaces it should admit and admits the ones it should reject"),
    }


# ---------------------------------------------------- S29-L56 chimera, approximate cell
def chimera():
    import glob
    fs = [f for f in sorted(glob.glob(os.path.join(BASE, "s29", "results",
                                                   "s29_X_probe_*.json")))
          if "fmt" not in f]
    sm, cb, pb = [], [], []
    for f in fs:
        o = json.load(open(f))["oracle"]
        sm.append(o["space_mean_rmsd"])
        cb.append(o["best_chimera_rmsd_cloud"])
        pb.append(o["pool_best_rmsd"])
    return {
        "n": len(fs),
        "space_mean_rmsd": float(np.mean(sm)),
        "chimera_best_ORACLE": float(np.mean(cb)),
        "pool_best_ORACLE": float(np.mean(pb)),
        "NOTE": ("the pool's MEAN on these 12 targets is not stored in any artefact, so the "
                 "chimera cell is an INVERSION, not a prediction: given the observed endpoint "
                 "deficit and the measured ceiling deficit, what Delta(set mean) does the law "
                 "require? Labelled approximate wherever it is quoted."),
        "d_ceiling": float(np.mean(cb) - np.mean(pb)),
    }


def main():
    kl, s24 = kladder(), s24_cells()
    out = {"frozen_coefficients": {"a_set_mean": A_MEAN, "b_set_best": B_BEST,
                                   "source": "Sprint 20 wide-set fit, 2142 cells"},
           "prereg": "s30/PREREG_S30_X.md @ d4305d17",
           "H_X1_primary_kladder": kl,
           "H_X1_secondary_s24": s24,
           "GATE_SCORECARD": gate_scorecard(kl, s24),
           "chimera_inversion_APPROXIMATE": chimera()}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tmp = OUT + ".%d.tmp" % os.getpid()          # contract rule 19: pid in the temp name
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1)
    os.replace(tmp, OUT)
    json.dump(out, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
