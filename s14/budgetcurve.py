"""SPRINT 14, coordinator -- does searching a GOOD objective harder help or hurt?

THE DECISIVE QUESTION FOR THE WHOLE SPRINT.  Sprint 13 established that optimising the
physical energies harder makes structures worse: on nine fully enumerated targets, RMSD
against evaluation budget runs 3.764 at 10 evaluations -> 3.667 at 300 -> 3.920 at the
CERTIFIED GLOBAL OPTIMUM, and the turn is governed by the fraction of space explored rather
than the evaluation count.  The stated consequence was that a VQE at n=12-16 with 1e4
evaluations sits on the improving limb of a curve whose limit is known to be bad, and will
LOOK like it is working.

That result was obtained on objectives with essentially no structural signal -- Legacy's
rank correlation inside its own low-energy decile is +0.043.  So it was confounded: it may
show that search hurts, or merely that searching a bad objective hurts.

`s14/hamil.py` removed the confound by building a native-free STRUCTURAL objective whose
in-decile rank correlation is **+0.370**, roughly eight times Legacy's and the best measured
anywhere in this project.  This file asks the question again on that objective:

    if the objective is good, does more search produce better structures?

If the curve descends monotonically, search is worth doing and a variational optimiser has a
real job.  If it turns over even here, then no optimiser -- classical or quantum -- can help
in this representation, and that is the sprint's answer independent of any ansatz.

METHOD.  Sample a large native-free configuration pool ONCE per target, then read the
argmin over increasing prefixes of it.  The budget curve therefore comes from a single pass
and every point on it is a strict subset of the next, which removes sampling noise between
budgets -- a paired design rather than independent runs.

Reported on BOTH axes always, because the project has established that within-target
correlation between energy reached and RMSD is -0.006 across six optimisers:
  * objective axis: the best objective value found
  * structural axis: the CA-RMSD of the configuration that value belongs to

Run:
    python -m s14.budgetcurve
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s14 import hamil as H                 # noqa: E402
from s14 import ladder as L                # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

BUDGETS = [10, 30, 100, 300, 1000, 3000, 10000, 20000]
W = 0.25              # the weight with the best in-decile rank correlation
N_MAX = 20000


def curve_for_target(t, w=W, n_max=N_MAX, seed=0):
    T = H.Terms(t["pdb"], t["seq"], t["n"], t["fold"])
    rng = np.random.default_rng(hash((t["pdb"], seed)) % (2 ** 32))

    # native-free proposal: half prior, half uniform, exactly as in hamil.py
    S1 = H._sample(T.P, n_max // 2, rng)
    S2 = rng.integers(0, T.k, size=(n_max - n_max // 2, T.n))
    S = np.vstack([S1, S2])
    rng.shuffle(S)                                   # so a prefix is an unbiased subsample

    ep, ed = H._z(T.e_prior(S)), H._z(T.e_disto(S))
    E = (1.0 - w) * ep + w * ed
    rm = T.ORACLE_rmsd(S)                            # ORACLE, post hoc only

    row = {"pdb": t["pdb"], "n": T.n, "fold": int(t["fold"]),
           "state_space": float(T.k ** T.n),
           "best_rmsd_in_pool": float(rm.min()), "points": {}}
    for B in BUDGETS:
        if B > len(S):
            continue
        k = int(np.argmin(E[:B]))
        row["points"][str(B)] = {
            "objective": float(E[k]),                        # objective axis
            "rmsd": float(rm[k]),                            # structural axis
            "rmsd_of_best_available": float(rm[:B].min()),   # ORACLE, the selection gap
            "frac_space_seen": float(B / (T.k ** T.n)),
        }
    return row


def run(targets=None, w=W):
    tg = targets if targets is not None else I.targets()
    rows = [curve_for_target(t, w) for t in tg]

    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    fail = np.isin(pdbs, I.FAIL18)
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)

    agg = {}
    for B in BUDGETS:
        k = str(B)
        got = [r for r in rows if k in r["points"]]
        if not got:
            continue
        rm = np.asarray([r["points"][k]["rmsd"] for r in got], float)
        ob = np.asarray([r["points"][k]["objective"] for r in got], float)
        av = np.asarray([r["points"][k]["rmsd_of_best_available"] for r in got], float)
        fl = np.asarray([r["points"][k]["frac_space_seen"] for r in got], float)
        f2 = np.isin([r["pdb"] for r in got], I.FAIL18)
        agg[k] = {"n": len(got),
                  "rmsd_mean": float(rm.mean()), "rmsd_median": float(np.median(rm)),
                  "frac_under_2": float((rm < 2.0).mean()),
                  "FAIL18": float(rm[f2].mean()) if f2.any() else None,
                  "objective_mean": float(ob.mean()),
                  "oracle_best_available": float(av.mean()),
                  "selection_gap": float(rm.mean() - av.mean()),
                  "mean_frac_space_seen": float(fl.mean())}

    # is the curve monotone in the structural axis?
    ks = [k for k in map(str, BUDGETS) if k in agg]
    ys = [agg[k]["rmsd_mean"] for k in ks]
    best_i = int(np.argmin(ys))
    turned = best_i < len(ys) - 1
    agg_meta = {"monotone_improving": not turned,
                "best_budget": int(ks[best_i]),
                "rmsd_at_best_budget": float(ys[best_i]),
                "rmsd_at_max_budget": float(ys[-1]),
                "degradation_past_minimum": float(ys[-1] - ys[best_i])}

    # split by fraction of space explored, the variable Sprint 13 identified as governing
    small = [r for r in rows if r["points"][ks[-1]]["frac_space_seen"] >= 0.01]
    large = [r for r in rows if r["points"][ks[-1]]["frac_space_seen"] < 0.001]

    def sub(rs):
        if not rs:
            return None
        y = [float(np.mean([r["points"][k]["rmsd"] for r in rs])) for k in ks]
        b = int(np.argmin(y))
        return {"n": len(rs), "curve": dict(zip(ks, y)), "best_budget": int(ks[b]),
                "degradation_past_minimum": float(y[-1] - y[b])}

    out = {"n": len(rows), "w_distogram": w, "budgets": BUDGETS,
           "aggregate": agg, "meta": agg_meta,
           "by_space_explored": {"seen_ge_1pct": sub(small), "seen_lt_0.1pct": sub(large)},
           "incumbent": float(ref.mean()),
           "per_target": rows}

    rm_max = np.asarray([r["points"][ks[-1]]["rmsd"] for r in rows], float)
    out["max_budget_vs_incumbent"] = I.paired(rm_max, ref, folds=folds, names=pdbs)

    with open(os.path.join(RESULTS, "budgetcurve.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_budgetcurve", out, n_expected=len(rows))

    print(f"native-free structural Hamiltonian, w(distogram) = {w}")
    print(f"incumbent {ref.mean():.3f}\n")
    print(f"{'budget':>8}{'objective':>12}{'RMSD':>9}{'<2A':>7}{'FAIL18':>9}"
          f"{'ORACLE best avail':>20}{'selection gap':>15}{'space seen':>12}")
    for k in ks:
        a = agg[k]
        f18 = f"{a['FAIL18']:.3f}" if a["FAIL18"] is not None else "  -  "
        print(f"{k:>8}{a['objective_mean']:>12.3f}{a['rmsd_mean']:>9.3f}"
              f"{a['frac_under_2']:>7.2f}{f18:>9}{a['oracle_best_available']:>20.3f}"
              f"{a['selection_gap']:>15.3f}{a['mean_frac_space_seen']:>12.2e}")
    m = agg_meta
    print(f"\nstructural axis {'DESCENDS MONOTONICALLY' if m['monotone_improving'] else 'TURNS OVER'}"
          f" -- minimum at budget {m['best_budget']} ({m['rmsd_at_best_budget']:.3f} A),"
          f" max budget {m['rmsd_at_max_budget']:.3f} A,"
          f" degradation {m['degradation_past_minimum']:+.3f} A")
    for tag, s in out["by_space_explored"].items():
        if s:
            print(f"  {tag:<18} n={s['n']:<4} min at {s['best_budget']:<6}"
                  f" degradation {s['degradation_past_minimum']:+.3f} A")
    v = out["max_budget_vs_incumbent"]
    print(f"\nargmin at max budget vs incumbent: {v['mean_diff']:+.3f} "
          f"[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]  {v['n_better']}W/{v['n_worse']}L")
    return out


if __name__ == "__main__":
    run()
