"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 8 -- "optimise harder, get worse", measured.

Two of my own measurements disagreed and the disagreement is the result.

    section 3  (9 targets, n=9, k=4, COMPLETE 262,144-configuration enumeration):
               Legacy's CERTIFIED GLOBAL OPTIMUM is +0.139 A WORSE than a random draw.
    section 8  (126 targets, k=4, 12,001 RANDOM configurations each):
               Legacy's best-of-12,001 is -0.617 A BETTER than a random draw.

Both are correct. They differ in how hard the objective was optimised. This file measures the
whole curve on the enumerated spaces, where the sample-budget argmin and the certified global
optimum are both available exactly: for m = 10, 100, ..., 262144, draw m configurations
uniformly, take the best by each objective, and record the CA-RMSD reached, averaged over
`reps` independent draws. If the curve falls and then RISES, "optimising harder makes the
structure worse" is not a slogan, it is a measured non-monotonicity, and the optimisation
budget is a hyperparameter of the ANSWER rather than of the search.

Output: `s13/results/qarch_budget.json`.

    python -m s13.qarch_budget
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402

BUDGETS = (10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 262144)
REPS = 40


def curve(F, rmsd, budgets=BUDGETS, reps=REPS, seed=0):
    rng = np.random.default_rng(seed)
    B = len(F)
    out = {}
    for m in budgets:
        if m > B:
            continue
        vals = []
        for _ in range(reps if m < B else 1):
            idx = rng.choice(B, m, replace=False) if m < B else np.arange(B)
            E = F[idx]
            tie = np.flatnonzero(E == E.min())
            vals.append(float(rmsd[idx][tie].mean()))
        out[m] = {"mean_rmsd": float(np.mean(vals)),
                  "sd": float(np.std(vals)), "reps": len(vals)}
    return out


def main():
    pdbs = sorted(f[len("qarch_enum_"):-4] for f in os.listdir(Q.RESULTS)
                  if f.startswith("qarch_enum_") and f.endswith(".npz"))
    rows = []
    for p in pdbs:
        z = np.load(os.path.join(Q.RESULTS, f"qarch_enum_{p}.npz"))
        rmsd = np.asarray(z["rmsd"], float)
        objs = {"legacy_total": np.asarray(z["legacy"], float),
                "prior_empirical": np.asarray(z["prior"], float),
                "leg_steric": np.asarray(z["leg_steric"], float),
                "ORACLE_rmsd": rmsd}
        rec = {"pdb": p, "n": int(z["n"]), "k": int(z["k"]),
               "pool_mean": float(rmsd.mean()), "best": float(rmsd.min()),
               "curves": {nm: curve(F, rmsd) for nm, F in objs.items()}}
        rows.append(rec)
        c = rec["curves"]["legacy_total"]
        print(f"  {p}: " + " ".join(f"{m}:{c[m]['mean_rmsd']:.2f}" for m in sorted(c)),
              flush=True)
    Q.write("qarch_budget", {"what": "CA-RMSD reached vs number of objective evaluations, "
                                    "on the complete enumerated spaces",
                             "budgets": list(BUDGETS), "reps": REPS, "rows": rows})
    print(f"\n{'budget':>8s} " + " ".join(f"{nm:>16s}" for nm in rows[0]["curves"]))
    for m in BUDGETS:
        line = f"{m:8d} "
        for nm in rows[0]["curves"]:
            v = [r["curves"][nm][m]["mean_rmsd"] for r in rows if m in r["curves"][nm]]
            line += f"{np.mean(v):16.3f} " if v else " " * 17
        print(line)
    print(f"{'pool mean':>8s} " + f"{np.mean([r['pool_mean'] for r in rows]):16.3f}")
    print(f"{'space best':>8s} " + f"{np.mean([r['best'] for r in rows]):16.3f}")
    return rows


if __name__ == "__main__":
    main()
