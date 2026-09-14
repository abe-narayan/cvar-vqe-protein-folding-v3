#!/usr/bin/env python
"""s26/w_bound_addendum.py -- the two L55 caveats applied to s26/results/w_selfcopy_bound.json
without editing s26/w_selfcopy.py (a job that imports it is running).

1. Adds the envelope's paired-gain row (leaked-model mean of arm - sel, minus the clean one,
   `ST.compare`, fold CI) to `signed_bounds_gated/C_envelope_fold_ci/gain`, and recomputes the
   `verdict/gain` class from it, so the artefact agrees with L44's pre-registered class (MINOR).
2. Adds the per-target readings of the same envelope: (2/60) x |leaked-mean minus clean| at the
   worst single target, at the p95 target, and over (target, model) pairs, on every basis, so
   that "0.028" is never quoted without "0.151 at the worst target (2BP4), under assumption A2".

The original artefact's provenance is kept under `provenance_original`; the file is re-saved with
this script's provenance.  Reads no native: every number comes from
s26/results/w_selfcopy_endpoint.json (gated, already on disk).

    python s26/jobrun.py --agent W --tag CPU --name w_bound_addendum --est-ram 0.2 -- python s26/w_bound_addendum.py
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
BOUND = os.path.join(RES, "w_selfcopy_bound.json")
ENDPOINT = os.path.join(RES, "w_selfcopy_endpoint.json")
N_BENCH, N_LEAKED = 60, 2
IMMATERIAL, CI_HALF = 0.017, 0.170


def materiality(x):
    if not np.isfinite(x):
        return "not measured"
    if x < IMMATERIAL:
        return "IMMATERIAL (< %.3f A, one tenth of the benchmark CI half-width %.3f)" % (IMMATERIAL, CI_HALF)
    if x < CI_HALF:
        return "MINOR (%.3f to %.3f A)" % (IMMATERIAL, CI_HALF)
    return "MATERIAL (>= %.3f A)" % CI_HALF


def main():
    B = json.load(open(BOUND)); E = json.load(open(ENDPOINT))
    rows = E["C"]["rows"]; pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    f = N_LEAKED / N_BENCH
    # ---- 1. the envelope's paired-gain row
    g_leak = np.array([r["leaked_mean"]["arm"] - r["leaked_mean"]["sel"] for r in rows])
    g_clean = np.array([r["clean"]["arm"] - r["clean"]["sel"] for r in rows])
    res = ST.compare(g_leak, g_clean, folds, names=pdbs,
                     label="Part C envelope: paired gain (arm - sel), leaked-model mean minus clean (ORACLE)")
    print(ST.fmt(res))
    env = B["signed_bounds_gated"]["C_envelope_fold_ci"]
    env["gain"] = {"effect": res["effect"], "mde": res["mde"], "ci95_fold": res["ci95_fold"], "verdict": res["verdict"],
                   "bound": float(f * max(abs(x) for x in res["ci95_fold"]))}
    # ---- 2. per-target readings of the envelope (L55 caveat 2), every basis
    per = {}
    for b in ("arm", "cloud", "sel", "fit", "gain"):
        if b == "gain":
            d = g_leak - g_clean
            dm = np.array([[(r["leaked"][j]["arm"] - r["leaked"][j]["sel"]) - (r["clean"]["arm"] - r["clean"]["sel"]) for j in r["leaked"]] for r in rows], dtype=object)
        else:
            d = np.array([r["leaked_mean"][b] - r["clean"][b] for r in rows])
            dm = np.array([[r["leaked"][j][b] - r["clean"][b] for j in r["leaked"]] for r in rows], dtype=object)
        flat = np.array([x for row in dm for x in row], float)
        k = int(np.argmax(np.abs(d)))
        per[b] = {"median_signed": float(np.median(d)), "p05_signed": float(np.percentile(d, 5)),
                  "worst_target": pdbs[k], "worst_abs": float(np.abs(d).max()), "p95_abs": float(np.percentile(np.abs(d), 95)),
                  "worst_abs_over_target_model_pairs": float(np.abs(flat).max()),
                  "bound_worst_target": float(f * np.abs(d).max()), "bound_p95_target": float(f * np.percentile(np.abs(d), 95)),
                  "bound_worst_pair": float(f * np.abs(flat).max()), "bound_mean_ci": env[b]["bound"] if b in env else None}
    B["signed_bounds_gated"]["C_envelope_per_target"] = per
    # ---- the verdict, recomputed with the gain row present; the per-target readings beside it
    for b in ("arm", "gain", "sel", "cloud"):
        cands = []
        sb = B["signed_bounds_gated"]
        if "both_removed4" in sb and b in sb["both_removed4"]:
            cands.append(sb["both_removed4"][b])
        elif "A_real4" in sb and b in sb["A_real4"]:
            cands.append(sb["A_real4"][b])
        if b in env:
            cands.append(env[b]["bound"])
        x = max(cands)
        B["verdict"][b] = {"bound_A": x, "class": materiality(x), "source": "signed",
                           "envelope_readings_A": {"mean_ci_limit": env[b]["bound"], "p95_target": per[b]["bound_p95_target"],
                                                   "worst_target": per[b]["bound_worst_target"], "worst_target_pdb": per[b]["worst_target"],
                                                   "worst_target_model_pair": per[b]["bound_worst_pair"]},
                           "class_under_every_reading": materiality(max(env[b]["bound"], per[b]["bound_worst_pair"])),
                           "assumption": "A2 of s26/PREREG_selfcopy_bound.md: the benchmark 2 are no worse than the dev envelope"}
    B["addendum_L55"] = {"applied": "2026-09-13", "what": "gain envelope row added; per-target envelope readings added; verdict recomputed",
                         "ledger": "L55 (Adversary), addendum entry by lane W"}
    B["provenance_original"] = B.get("provenance")
    ST.save_atomic(BOUND, B, module_file=__file__)
    print(json.dumps({"verdict": B["verdict"], "per_target_arm": per["arm"], "per_target_gain": per["gain"]}, indent=1))
    print("  ->", BOUND)
    return 0


if __name__ == "__main__":
    sys.exit(main())
