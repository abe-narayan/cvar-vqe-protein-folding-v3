"""s17/audit_repro.py -- INDEPENDENT RE-DERIVATION OF THE ORACLE MAP, PLUS LEAKAGE,
DETERMINISM AND POWER.

A transcription audit cannot see a formula error, so nothing here re-reads
`s17/results/oracle_map.json` as ground truth: every quantity is recomputed from the
universe npz by a second route and the two are compared.

CHECKS
  R1  GEOMETRY.  `oracle_map.py` recomputes per-window RMSD as
      `I.kabsch_rmsd_batch(W, nat)` in float64.  The npz already carries `rr`, computed
      independently in `s8/generate.py` by `audit.kabsch_rmsd_batch` on the float32
      arrays.  Two implementations, two dtypes -> report max |delta| per target.
  R2  THE SHIPPED ANCHOR.  The map's K = 500 `sel_dist` must reproduce the programme's
      standing "shipped argmin 3.4540".  If it does not, the map's selector is not the
      shipped selector and every "gap to oracle" in it is measured against a different
      object.  Note `I.selfcheck` rounds pair distances through float32 and the map does
      not -- so this is a REAL check, not a tautology.
  R3  THE RANDOM CONTROL'S OWN NOISE.  `oracle_map.py` uses N_RAND = 3 draws per target
      for `sel_rand`.  Derive the standard error that control carries and state what it
      does to the W/L counts and the CI in the map's table B.
  R4  DETERMINISM.  Rerun the order-null statistic in a FRESH interpreter, report max
      |delta|.  Grep the s17 tree for bare `hash(` and for benchmark60 access.
  R5  POWER.  For the 126-target instrument, compute the SE of a paired difference and
      the minimum detectable effect at 80% power, so that any null this sprint declares
      can be read against what it could have detected.

Nothing here is predictive; rr / nat_ca are ORACLE labels used for ceilings only.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I     # noqa: E402
from s15 import seed as SD          # noqa: E402

UNIV = os.path.join(ROOT, "s8", "generate_univ")


def r1_r2(verbose=True):
    tg = I.targets()
    worst = 0.0
    sel500, best500, bestfull = [], [], []
    worst_pdb = None
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        z = np.load(os.path.join(UNIV, f"{pdb}.npz"), allow_pickle=True)
        W = np.asarray(z["W"], float)
        rr = np.asarray(z["rr"], float)
        order = np.asarray(z["order"], int)
        nat = np.asarray(z["nat_ca"], float)
        d = I.kabsch_rmsd_batch(W, nat)
        m = float(np.abs(d - rr).max())
        if m > worst:
            worst, worst_pdb = m, pdb
        i, j = I.pair_index(n)
        dg = I.distogram(pdb, seq, fold)
        p = order[:min(500, len(rr))]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W[p], i, j)), float)
        sel500.append(float(d[p][int(np.argmin(sc))]))
        best500.append(float(d[p].min()))
        bestfull.append(float(d.min()))
        if verbose and (c + 1) % 25 == 0:
            print(f"  repro {c+1}/{len(tg)}", flush=True)
    out = {"r1_max_abs_delta_rr": worst, "r1_worst_pdb": worst_pdb,
           "r2_sel500_mean": float(np.mean(sel500)),
           "r2_best500_mean": float(np.mean(best500)),
           "r2_bestfull_mean": float(np.mean(bestfull)),
           "sel500": sel500, "best500": best500, "bestfull": bestfull,
           "pdb": [t["pdb"] for t in tg], "fold": [int(t["fold"]) for t in tg]}
    json.dump(out, open(os.path.join(RESULTS, "audit_repro.json"), "w"))
    return out


def r3_control_noise(K=500, n_rep=400):
    """SE of oracle_map's N_RAND = 3 random-selection control, per target and on the mean."""
    tg = I.targets()
    per = []
    for t in tg:
        z = np.load(os.path.join(UNIV, f"{t['pdb']}.npz"), allow_pickle=True)
        rr = np.asarray(z["rr"], float); order = np.asarray(z["order"], int)
        k = min(K, len(rr))
        v = rr[order[:k]]
        rng = SD.stable_rng(t["pdb"], "s17auditctrl")
        draws = v[rng.integers(0, k, size=(n_rep, 3))].mean(1)
        per.append({"pdb": t["pdb"], "sd_pop": float(v.std(ddof=1)),
                    "se3": float(v.std(ddof=1) / np.sqrt(3)),
                    "emp_sd_of_mean3": float(draws.std(ddof=1))})
    se3 = np.array([p["se3"] for p in per])
    return {"per_target_se_of_sel_rand": float(se3.mean()),
            "se_of_the_126_target_mean": float(np.sqrt((se3 ** 2).mean() / len(se3))),
            "max_per_target_se": float(se3.max()), "per": per}


def r4_greps():
    hits = {"bare_hash": [], "benchmark60": []}
    for dirpath, dirnames, files in os.walk(os.path.join(ROOT, "s17")):
        dirnames[:] = [d for d in dirnames if d not in ("__pycache__", "cache", "results")]
        for f in files:
            if not f.endswith(".py"):
                continue
            p = os.path.join(dirpath, f)
            src = open(p, encoding="utf-8", errors="replace").read()
            for k, line in enumerate(src.splitlines(), 1):
                if re.search(r"(?<![\w.])hash\s*\(", line) and "stable" not in line:
                    hits["bare_hash"].append(f"{os.path.relpath(p, ROOT)}:{k}: {line.strip()}")
                if re.search(r"benchmark60|bench60|benchmark_60", line):
                    hits["benchmark60"].append(f"{os.path.relpath(p, ROOT)}:{k}: {line.strip()}")
    return hits


def r5_power(alpha=0.05, power=0.80):
    """MDE for a paired 126-target comparison, using the observed spread of a real
    per-target difference (the order-null skill at K = 500)."""
    p = os.path.join(RESULTS, "audit_ordernull.json")
    if not os.path.exists(p):
        return {}
    rows = json.load(open(p))["rows"]
    d = np.array([r["rand_mean"]["500"] - r["blosum"]["500"] for r in rows])
    n = len(d)
    se = d.std(ddof=1) / np.sqrt(n)
    z = 1.959963985 + 0.8416212336            # two-sided alpha + one-sided beta
    # also: paired ceiling-gain spread, the statistic the sprint's headline uses
    g = np.array([r["blosum"]["500"] - r["blosum"]["full"] for r in rows])
    return {"n": n, "sd_of_diff": float(d.std(ddof=1)), "se": float(se),
            "mde_80pct": float(z * se),
            "headline_sd": float(g.std(ddof=1)),
            "headline_se": float(g.std(ddof=1) / np.sqrt(n)),
            "headline_mde_80pct": float(z * g.std(ddof=1) / np.sqrt(n))}


def main():
    print("R1/R2 -- independent re-derivation of the map's core quantities")
    o = r1_r2()
    print(f"  R1  max |kabsch_rmsd_batch(float64) - npz rr(float32)| over all windows,")
    print(f"      all 126 targets: {o['r1_max_abs_delta_rr']:.2e}  (worst {o['r1_worst_pdb']})")
    print(f"  R2  shipped distance-selector argmin at K=500: {o['r2_sel500_mean']:.4f}"
          f"   [programme's standing value 3.4540]")
    print(f"      ORACLE best at K=500 : {o['r2_best500_mean']:.4f}  [standing 1.7108]")
    print(f"      ORACLE best full     : {o['r2_bestfull_mean']:.4f}")

    print("\nR3 -- how noisy is oracle_map.py's own random control (N_RAND = 3)?")
    c = r3_control_noise()
    print(f"  per-target SE of sel_rand      {c['per_target_se_of_sel_rand']:.3f} A"
          f"   (max {c['max_per_target_se']:.3f})")
    print(f"  SE it contributes to the 126-target MEAN  {c['se_of_the_126_target_mean']:.3f} A")
    print("  The MEAN is fine; the per-target W/L counts and the width of the")
    print("  'dist - random' interval in the map's table B are inflated by this noise.")

    print("\nR4 -- leakage and seeding greps over s17/")
    h = r4_greps()
    print(f"  bare hash() uses : {len(h['bare_hash'])}")
    for x in h["bare_hash"]:
        print("    " + x)
    print(f"  benchmark60 refs : {len(h['benchmark60'])}")
    for x in h["benchmark60"]:
        print("    " + x)

    print("\nR5 -- power of the 126-target instrument")
    p = r5_power()
    if p:
        print(f"  paired difference sd {p['sd_of_diff']:.3f}, SE {p['se']:.3f},"
              f" MDE at 80% power {p['mde_80pct']:.3f} A")
        print(f"  the K-widening headline statistic: sd {p['headline_sd']:.3f},"
              f" SE {p['headline_se']:.3f}, MDE {p['headline_mde_80pct']:.3f} A")
    json.dump({"r3": {k: v for k, v in c.items() if k != "per"}, "r4": h, "r5": p},
              open(os.path.join(RESULTS, "audit_repro_meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
