#!/usr/bin/env python
"""s26/w_tiebreak_report.py -- the ledger-facing numbers of the tie-break floor from
s26/results/w_selfcopy_tiebreak_endpoint.json (written by `w_tiebreak.py endpoint`): the
draw-0-vs-draw-1 paired contrast through ST.fmt (verbatim in the ledger), the recorded effects
placed against the floor, and the native-free membership statistics from the draws artefact.
Reads no native (every RMSD is already in the endpoint artefact).

    python s26/w_tiebreak_report.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s24 import stats_lib as ST                      # noqa: E402

RES = os.path.join(HERE, "results")
#: recorded hundredths-level effects (A, built chain unless stated), with their source
RECORDED = [
    ("production AMBER relaxation minus do-nothing (L39, C3 stage 1)", 0.0207),
    ("AMBER minus matched random move (L39)", 0.0111),
    ("restrained AMBER relaxation of the average, S24 (memory)", 0.022),
    ("ORACLE functional weight, S24 L16", 0.015),
    ("all-atom AMBER reranking, S17 (memory)", 0.004),
    ("C27 dev leak price, lam = 0 chain (S10-4, L44)", 0.0004),
    ("C27 dev leak price, built chain (L44)", 0.0018),
    ("projection mode fd vs exact, rmsd_arm max (L19)", 0.012),
    ("ph_c3 predicted AMBER-vs-random (L24)", 0.013),
]


def main():
    E = json.load(open(os.path.join(RES, "w_selfcopy_tiebreak_endpoint.json")))
    Z = json.load(open(os.path.join(RES, "w_selfcopy_tiebreak_draws.json")))
    rows = E["rows"]; pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    out = {"label": "tie-break floor, ledger numbers", "n": len(rows), "source": ["w_selfcopy_tiebreak_endpoint.json", "w_selfcopy_tiebreak_draws.json"]}
    # ---- membership statistics (native-free)
    rep = np.array([[d["top75_replaced"] for d in r["draws"]] for r in Z["rows"]], float)
    same = np.array([[d["argmin_same"] for d in r["draws"]] for r in Z["rows"]], float)
    n_in = np.array([r["n_in"] for r in Z["rows"]], float); n_tied = np.array([r["n_tied"] for r in Z["rows"]], float)
    out["membership"] = {"top75_replaced_mean": float(rep.mean()), "top75_replaced_median": float(np.median(rep)),
                         "top75_replaced_max": float(rep.max()), "argmin_same_frac": float(same.mean()),
                         "n_in_pool_median": float(np.median(n_in)), "n_tied_median": float(np.median(n_tied)),
                         "targets_with_zero_replacement_on_every_draw": int((rep.max(1) == 0).sum())}
    print("membership:", json.dumps(out["membership"]))
    # ---- per basis: the two headline contrasts through ST.fmt
    for b in ("arm", "cloud", "sel", "fit"):
        M = np.array([[e[b] for e in r["draws"]] for r in rows])
        prod = np.array([r["production"][b] for r in rows])
        r01 = ST.compare(M[:, 0], M[:, 1], folds, names=pdbs, label="tie-break draw 0 minus draw 1 [%s], n=126" % b)
        rpm = ST.compare(prod, M.mean(1), folds, names=pdbs, label="production convention minus mean over 8 draws [%s], n=126" % b)
        print(ST.fmt(r01)); print(ST.fmt(rpm))
        s = E["stats"][b]
        out[b] = {"draw0_vs_draw1": r01, "production_vs_draws": rpm, "summary": {k: v for k, v in s.items() if k != "production_vs_draws"}}
        print("  [%s] m_tie (sd of the 126-mean over draws) %.4f | mean over draws %.4f production %.4f | s_tie median %.4f p90 %.4f | paired MDE between draws median %.4f max %.4f | max |paired effect| %.4f"
              % (b, s["m_tie_sd_of_mean_over_draws"], s["mean_over_draws"], s["production_mean"], s["s_tie_median"], s["s_tie_p90"],
                 s["paired_mde_between_draws_median"], s["paired_mde_between_draws_max"], s["paired_effect_between_draws_maxabs"]))
        # per-target sd on the built chain: how many targets move more than 0.1 / 0.5 A
        if b == "arm":
            st = np.array([r["s_tie"][b] for r in rows])
            rng_t = M.max(1) - M.min(1)
            out["arm_per_target"] = {"s_tie_gt_0.1": int((st > 0.1).sum()), "s_tie_gt_0.5": int((st > 0.5).sum()),
                                     "range_over_draws_median": float(np.median(rng_t)), "range_over_draws_p90": float(np.percentile(rng_t, 90)),
                                     "range_over_draws_max": float(rng_t.max()), "worst_target": pdbs[int(np.argmax(rng_t))]}
            print("  arm per target: s_tie > 0.1 A on %d/126, > 0.5 A on %d/126; range over 8 draws median %.3f p90 %.3f max %.3f (%s)"
                  % (out["arm_per_target"]["s_tie_gt_0.1"], out["arm_per_target"]["s_tie_gt_0.5"], out["arm_per_target"]["range_over_draws_median"],
                     out["arm_per_target"]["range_over_draws_p90"], out["arm_per_target"]["range_over_draws_max"], out["arm_per_target"]["worst_target"]))
    # ---- the recorded effects against the floor (built chain)
    s = E["stats"]["arm"]
    floor_mde = s["paired_mde_between_draws_median"]; floor_m = s["m_tie_sd_of_mean_over_draws"]
    table = []
    for name, val in RECORDED:
        table.append({"effect": name, "value_A": val, "below_paired_mde_between_draws": bool(val < floor_mde),
                      "below_2_m_tie": bool(val < 2 * floor_m)})
    out["recorded_effects_vs_floor"] = {"paired_mde_between_draws_median_arm": floor_mde, "m_tie_arm": floor_m, "table": table}
    print("recorded effects against the floor (built chain): paired MDE between draws %.4f, 2 x m_tie %.4f" % (floor_mde, 2 * floor_m))
    for t in table:
        print("  %-62s %.4f   below draw MDE: %s   below 2 m_tie: %s" % (t["effect"], t["value_A"], t["below_paired_mde_between_draws"], t["below_2_m_tie"]))
    ST.save_atomic(os.path.join(RES, "w_tiebreak_report.json"), out, module_file=__file__)
    print("  ->", os.path.join(RES, "w_tiebreak_report.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
