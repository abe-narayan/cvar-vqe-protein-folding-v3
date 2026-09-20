"""S30 lane X, H-X3: is a "better typical member" a BIAS improvement or a SPREAD reduction?

Pre-registration: `s30/PREREG_S30_X.md` ADDENDUM 1, committed at 539de5a8 BEFORE this file
existed and before any per-target number below was computed.

The coordinator asked: can a generator be built whose TYPICAL member beats the pool's, rather
than one whose BEST member does?  Measure the set mean first.

Measuring it shows it is two quantities.  For a coordinate-average terminal over m members,
the common-mode identity (exact to 2.7e-14, `pool-error-is-68-percent-common-mode`) gives

    set_mean^2  ~=  B^2 + S^2          B = RMSD(set average, native) = THE ENDPOINT
    endpoint     =  B                  S = the set's spread about its own centroid

so a set mean improved purely by CONCENTRATION (S down, B fixed) is worth exactly zero.
Only B converts, and B is the endpoint.

`set_mean` here is the arithmetic mean of member RMSDs, not the root-mean-square, so S is a
Jensen-biased estimate of the spread.  The direction of that bias is stated with the result and
it cannot manufacture the effect: it makes S too SMALL, i.e. it understates concentration.

504 (arm, target) cells from S24 lane C plus the shipped pool.  Stored artefacts only.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "s30", "results", "s30_X_typicalgood.json")
ARMS = ("T0_helix", "T1_blind", "T2_restype", "T3_pool")


def spread(set_mean, b):
    return math.sqrt(max(set_mean * set_mean - b * b, 0.0))


def fold_ci(vals, folds, n_boot=4000, seed=0):
    vals, folds = np.asarray(vals, float), np.asarray(folds)
    rng = np.random.default_rng(seed)
    fs = sorted(set(folds.tolist()))
    bs = []
    for _ in range(n_boot):
        pick = rng.choice(fs, len(fs), replace=True)
        idx = np.concatenate([np.where(folds == f)[0] for f in pick])
        bs.append(vals[idx].mean())
    lo, hi = np.percentile(bs, [2.5, 97.5])
    per = {int(f): float(vals[folds == f].mean()) for f in fs}
    same = sum(1 for f in fs if np.sign(per[f]) == np.sign(vals.mean()))
    return {"mean": float(vals.mean()), "median": float(np.median(vals)),
            "se": float(vals.std(ddof=1) / math.sqrt(len(vals))),
            "mde": float(2.8016 * vals.std(ddof=1) / math.sqrt(len(vals))),
            "fold_CI95": [float(lo), float(hi)],
            "folds_same_sign": "%d/%d" % (same, len(fs)), "per_fold": per,
            "n": int(len(vals))}


def main():
    rows = json.load(open(os.path.join(BASE, "s24", "results", "c_ladder.json")))["rows"]
    bench = {r["pdb"]: r for r in json.load(open(
        os.path.join(BASE, "bench_results", "baseline_tuning126.json")))["per_target"]}

    cells, pool = [], []
    for r in rows:
        b = bench.get(r["pdb"])
        if b is None:
            continue
        p_mean, p_B = b["top_m_mean"], r["incumbent"]     # SAME set: the shipped top-75
        p_S = spread(p_mean, p_B)
        pool.append({"pdb": r["pdb"], "fold": int(r["fold"]),
                     "set_mean": p_mean, "B_endpoint": p_B, "S_spread": p_S,
                     "avg_gain": p_mean - p_B})
        for a in ARMS:
            d = r[a]
            g_mean, g_B = d["sel_member_mean_ORACLE"], d["rmsd"]
            g_S = spread(g_mean, g_B)
            cells.append({
                "pdb": r["pdb"], "arm": a, "fold": int(r["fold"]),
                "set_mean": g_mean, "B_endpoint": g_B, "S_spread": g_S,
                "avg_gain": g_mean - g_B,
                "d_set_mean": g_mean - p_mean,
                "d_B": g_B - p_B,
                "d_S": g_S - p_S,
                "wins_set_mean": int(g_mean < p_mean),
                "wins_endpoint": int(g_B < p_B),
                "wins_BOTH": int(g_mean < p_mean and g_B < p_B),
            })

    fold = np.array([c["fold"] for c in cells])
    out = {"prereg": "s30/PREREG_S30_X.md ADDENDUM 1 @ 539de5a8",
           "n_cells": len(cells), "n_pool": len(pool)}

    out["pool"] = {k: float(np.mean([p[k] for p in pool]))
                   for k in ("set_mean", "B_endpoint", "S_spread", "avg_gain")}
    out["per_arm"] = {}
    for a in ARMS:
        cs = [c for c in cells if c["arm"] == a]
        f = np.array([c["fold"] for c in cs])
        out["per_arm"][a] = {
            "set_mean": float(np.mean([c["set_mean"] for c in cs])),
            "B_endpoint": float(np.mean([c["B_endpoint"] for c in cs])),
            "S_spread": float(np.mean([c["S_spread"] for c in cs])),
            "avg_gain": float(np.mean([c["avg_gain"] for c in cs])),
            "d_set_mean": fold_ci([c["d_set_mean"] for c in cs], f),
            "d_B": fold_ci([c["d_B"] for c in cs], f),
            "d_S": fold_ci([c["d_S"] for c in cs], f),
            "n_wins_set_mean": int(sum(c["wins_set_mean"] for c in cs)),
            "n_wins_endpoint": int(sum(c["wins_endpoint"] for c in cs)),
            "n_wins_BOTH": int(sum(c["wins_BOTH"] for c in cs)),
        }

    # ---- PRIMARY falsifier: wins on BOTH, and whether it transfers split-half
    both = np.array([c["wins_BOTH"] for c in cells])
    frac_both = float(both.mean())
    rng = np.random.default_rng(7)
    pdbs = sorted({c["pdb"] for c in cells})
    tr = []
    for _ in range(400):                       # split-half transfer on the BOTH-winning cells
        perm = rng.permutation(pdbs)
        h1, h2 = set(perm[:len(perm) // 2]), set(perm[len(perm) // 2:])
        # pick the arm that wins on both most often in half 1, score it in half 2
        best, bv = None, -1
        for a in ARMS:
            v = np.mean([c["wins_BOTH"] for c in cells
                         if c["arm"] == a and c["pdb"] in h1])
            if v > bv:
                best, bv = a, v
        tr.append(np.mean([c["d_B"] for c in cells
                           if c["arm"] == best and c["pdb"] in h2]))
    out["PRIMARY_falsifier"] = {
        "frac_cells_winning_BOTH_set_mean_and_endpoint": frac_both,
        "n_cells_winning_BOTH": int(both.sum()),
        "REFUTED_if_frac_ge": 0.25,
        "split_half_transfer_d_B_mean": float(np.mean(tr)),
        "split_half_transfer_d_B_CI95": [float(np.percentile(tr, 2.5)),
                                         float(np.percentile(tr, 97.5))],
        "verdict": ("REFUTED" if frac_both >= 0.25 and np.mean(tr) < 0 else "H-X3 STANDS"),
    }

    # ---- SECONDARY: within-target corr(S, B).  Does concentration buy bias?
    xs, ys = [], []
    for t in pdbs:
        cs = [c for c in cells if c["pdb"] == t]
        p = [q for q in pool if q["pdb"] == t][0]
        S = np.array([c["S_spread"] for c in cs] + [p["S_spread"]])
        B = np.array([c["B_endpoint"] for c in cs] + [p["B_endpoint"]])
        xs.append(S - S.mean())
        ys.append(B - B.mean())
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    out["SECONDARY_within_target_corr_S_B"] = {
        "r": float(np.corrcoef(xs, ys)[0, 1]), "n": int(len(xs)),
        "registered_prediction": "near zero or negative",
        "reading": ("positive r would mean a more concentrated set also has a smaller bias, "
                    "i.e. concentration buys the endpoint; negative means the opposite"),
    }

    # ---- TERTIARY: does the averaging gain scale with S?  (the terminal extracts spread)
    allS = np.array([c["S_spread"] for c in cells] + [p["S_spread"] for p in pool])
    allG = np.array([c["avg_gain"] for c in cells] + [p["avg_gain"] for p in pool])
    sl, ic = np.polyfit(allS, allG, 1)
    out["TERTIARY_avg_gain_vs_S"] = {
        "slope": float(sl), "intercept": float(ic),
        "r": float(np.corrcoef(allS, allG)[0, 1]), "n": int(len(allS)),
        "reading": ("if the averaging gain is a function of S alone, the terminal's whole value "
                    "is spread extraction and a CONCENTRATED source hands it nothing to extract"),
    }

    # ---- the decisive contrast: the most CONCENTRATED source in the record
    conc = min(out["per_arm"], key=lambda a: out["per_arm"][a]["S_spread"])
    out["MOST_CONCENTRATED_SOURCE"] = {
        "arm": conc, "S_spread": out["per_arm"][conc]["S_spread"],
        "pool_S_spread": out["pool"]["S_spread"],
        "its_endpoint": out["per_arm"][conc]["B_endpoint"],
        "pool_endpoint": out["pool"]["B_endpoint"],
        "is_it_the_worst_endpoint": bool(
            out["per_arm"][conc]["B_endpoint"]
            == max(v["B_endpoint"] for v in out["per_arm"].values())),
    }

    out["JENSEN_NOTE"] = ("set_mean is the arithmetic mean of member RMSDs, not the RMS, so S is "
                          "biased SMALL -- it understates concentration and cannot manufacture "
                          "the effect reported here.")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump({**out, "cells": cells, "pool": pool}, fh, indent=1)
    os.replace(tmp, OUT)
    json.dump(out, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
