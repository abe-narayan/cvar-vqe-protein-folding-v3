"""Do the RECORDED scientific numbers reproduce on the integrated path?

Two instruments:
  * the 126-target tuning instrument, against `s9/synth_cache/full_sc_75_*.json`, which is
    the s9 reference the record was written from
  * the 60-target benchmark, against `s9/final_report.json` -- READ ONLY. The single
    pre-registered pass is spent; this file never runs it.

Per-target, per-arm, so a mean that matches for two offsetting reasons cannot pass.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

#: recorded -> (harness field, s9 reference field)
ARMS = {"shipped": ("shipped", "score"),
        "pool_best": ("pool_best", "pool_best"),
        "top_m_best": ("top_m_best", "sub_best"),
        "top_m_mean": ("top_m_mean", "sub_mean"),
        "rmsd_avg": ("rmsd_avg", "avg"),
        "rmsd_fit": ("rmsd_fit", "fit")}

RECORDED_126 = {"shipped": 3.4540, "pool_best": 1.7108, "top_m_best": 2.306,
                "rmsd_fit": 3.2005, "rmsd_avg": 3.0483}
RECORDED_60 = {"full": 2.9610, "shipped": 2.9507, "pool_best": 1.4666}


def load_s9_126():
    out = {}
    for f in glob.glob(os.path.join(_ROOT, "s9", "synth_cache", "full_sc_75_*.json")):
        d = json.load(open(f))
        out[d["pdb"]] = d
    return out


def load_harness_126():
    p = os.path.join(_ROOT, "bench_results", "baseline_tuning126.json")
    if not os.path.exists(p):
        return {}, None
    d = json.load(open(p))
    return {r["pdb"]: r for r in d["per_target"]}, d


def main():
    out = {}
    s9 = load_s9_126()
    hr, meta = load_harness_126()
    out["n_s9_reference"] = len(s9)
    out["n_harness"] = len(hr)
    out["harness_backends"] = (meta or {}).get("backends")
    out["harness_cfg"] = (meta or {}).get("config")

    common = sorted(set(s9) & set(hr))
    out["n_common"] = len(common)

    per_arm = {}
    for name, (hk, sk) in ARMS.items():
        a, b, pdbs = [], [], []
        for p in common:
            hv = hr[p].get(hk)
            sv = s9[p].get(sk, s9[p].get("rmsd", {}).get(sk)
                           if isinstance(s9[p].get("rmsd"), dict) else None)
            if isinstance(s9[p].get("rmsd"), dict) and sk in s9[p]["rmsd"]:
                sv = s9[p]["rmsd"][sk]
            if hv is None or sv is None:
                continue
            a.append(float(hv))
            b.append(float(sv))
            pdbs.append(p)
        if not a:
            continue
        a, b = np.array(a), np.array(b)
        d = np.abs(a - b)
        worst = int(np.argmax(d))
        per_arm[name] = {
            "n": len(a),
            "harness_mean": float(a.mean()),
            "s9_reference_mean": float(b.mean()),
            "mean_diff": float(a.mean() - b.mean()),
            "max_abs_per_target": float(d.max()),
            "worst_target": pdbs[worst],
            "worst_harness": float(a[worst]),
            "worst_s9": float(b[worst]),
            "n_bit_identical": int((a == b).sum()),
            "n_within_1e-9": int((d < 1e-9).sum()),
        }
    out["arms_126"] = per_arm

    # the recorded constants
    rec = {}
    for k, want in RECORDED_126.items():
        got = per_arm.get(k, {}).get("harness_mean")
        rec[k] = {"recorded": want, "harness": got,
                  "abs_diff": (None if got is None else abs(got - want)),
                  "matches_to_4dp": (None if got is None
                                     else abs(got - want) < 5e-4)}
    out["recorded_126"] = rec

    # -- the 60-target benchmark, READ ONLY -------------------------------
    fr = os.path.join(_ROOT, "s9", "final_report.json")
    if os.path.exists(fr):
        d = json.load(open(fr))
        pt = d["per_target"]
        got = {"full": float(np.mean([r["full"] for r in pt])),
               "shipped": float(np.mean([r["shipped"] for r in pt])),
               "pool_best": float(np.mean([r["pool_best"] for r in pt]))}
        out["benchmark60"] = {
            "n": d["n"], "source": "s9/final_report.json (cached, NOT re-run)",
            "preregistration": d.get("preregistration"),
            "arms": {k: {"recorded": RECORDED_60[k], "cached": v,
                         "abs_diff": abs(v - RECORDED_60[k]),
                         "matches_to_4dp": abs(v - RECORDED_60[k]) < 5e-4}
                     for k, v in got.items()}}

    print(json.dumps(out, indent=2, sort_keys=True, default=float))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "headline_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=float)
    return out


if __name__ == "__main__":
    main()
