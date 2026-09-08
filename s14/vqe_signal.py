"""SPRINT 14 / VQE -- THE DECISIVE EXPERIMENT.

    At what OBJECTIVE QUALITY does VQE/CVaR start to beat classical search, and does it
    ever beat it on STRUCTURE rather than on ENERGY?

Sprint 13 showed VQE would lose, but only ever tested objectives with no structural signal,
and on such an objective better search makes structures WORSE (the certified Legacy optimum
is +0.139 A worse than a random draw).  "VQE lost" was therefore confounded with "the
objective was bad".  This deconfounds it by making objective quality a controlled variable.

THE DESIGN
----------
Nine fully enumerated targets (n=9, k=4, 262,144 configurations, true CA-RMSD on every
one).  A signal-tunable objective family blends a realistic objective toward the true RMSD
in RANK space, so that objective quality moves and nothing else does; the realised rank
correlation is MEASURED at every setting, never assumed.  At each signal level, on the same
space, at MATCHED budgets counted in objective evaluations:

    random | greedy 1-opt | simulated annealing | GA | EXACT ENUMERATION
    | expectation-value VQE (alpha=1) | CVaR-VQE at alpha = 0.25, 0.1, 0.05

Two axes, always separate and never substituted for one another:
    objective gap    E_found - E_exact          (did the optimiser optimise?)
    structural gap   RMSD_returned - RMSD_best  (did the structure get better?)

THE TRAPS THIS IS BUILT TO AVOID
--------------------------------
* "Beats random" proves nothing: a zero-information constant alpha-helix beats uniform
  random by 0.457 A.  Every arm is therefore reported against the CERTIFIED OPTIMUM of its
  own objective and against the space best, not against random.
* RMSD vs budget is non-monotone: 3.764 at 10 evaluations, 3.667 at 300, 3.920 at the
  certified optimum.  A VQE on the improving limb looks like it works.  Every arm here is
  run at several budgets AND the exact optimum is computed, so the limb and the limit are
  both visible in the same table.
* An initialisation that hands over the answer is not a result: `init_rmsd_mean` is
  recorded on every VQE run and printed next to the outcome.

    python -m s14.vqe_signal            # full grid, ~2-3 h, one heavy process
    python -m s14.vqe_signal --quick    # smoke test
"""
from __future__ import annotations

import sys
import time

import numpy as np

from s14 import vqe_lib as V
from s14 import vqe_run as R

SIGNALS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)
BUDGETS = (3000, 30000)
ALPHAS = (1.0, 0.25, 0.1, 0.05)
SEEDS = (0, 1, 2)
ANSATZ = "mps2f"
SHOTS = 512


def cell(e: V.Enum, E: np.ndarray, budget: int, seeds=SEEDS, alphas=ALPHAS,
         ansatz=ANSATZ, shots=SHOTS) -> dict:
    """Every arm on one (target, objective, budget) cell."""
    exact_i = int(np.argmin(E))
    out = {"exact": {"evals": e.N, "objective_gap": 0.0,
                     "rmsd_returned": float(e.rmsd[exact_i]),
                     "structural_gap": float(e.rmsd[exact_i] - e.rmsd.min()),
                     "best_e": float(E[exact_i])},
           "rho_obj_rmsd": V.spearman(E, e.rmsd),
           "space_best_rmsd": float(e.rmsd.min()),
           "space_mean_rmsd": float(e.rmsd.mean())}
    for name, fn in V.CLASSICAL.items():
        rs = []
        for s in seeds:
            c = fn(E, e.n, e.k, budget, np.random.default_rng(1000 * s + 7))
            rs.append(V.summarise(c, e.rmsd, E, exact_i))
        out[name] = _agg(rs)
    for a in alphas:
        rs = []
        for s in seeds:
            # The exact final distribution (entropy, near-native probability mass) costs an
            # MPS contraction over all 4^n basis states -- 3 s uncontended, far more on a
            # loaded box, and it dominates the run.  It is a property of the converged state,
            # so one seed prices it; the other seeds still contribute every sampled quantity.
            r = R.run(E, e.n_qubits, a, budget=budget, shots=min(shots, budget),
                      ansatz=ansatz, seed=s, rmsd=e.rmsd,
                      exact_dist=(s == seeds[0]))
            r["objective_gap"] = float(r["best_e"] - E[exact_i])
            r["structural_gap"] = float(r["rmsd_returned"] - e.rmsd.min())
            if "mode_rmsd" in r:            # only the seed that priced the exact state
                r["mode_structural_gap"] = float(r["mode_rmsd"] - e.rmsd.min())
            rs.append(r)
        out[f"vqe_a{a}"] = _agg(rs)
    return out


_AGG_KEYS = ("objective_gap", "structural_gap", "rmsd_returned", "rmsd_mean_seen",
             "rmsd_best_seen", "n_distinct", "diversity", "entropy_bits",
             "mode_rmsd", "mode_structural_gap", "mean_rmsd_under_p",
             "pmass_below_2.0", "pmass_below_2.5", "pmass_below_3.0",
             "frac_below_2.0", "frac_below_2.5", "frac_below_3.0",
             "init_rmsd_mean", "evals", "grad_var", "max_prob", "eff_support")


def _agg(rs):
    o = {"n_seeds": len(rs)}
    for k in _AGG_KEYS:
        v = [r[k] for r in rs if k in r and np.isfinite(r[k])]
        if v:
            o[k] = float(np.mean(v))
            o[k + "_sd"] = float(np.std(v))
    o["per_seed_structural_gap"] = [r.get("structural_gap") for r in rs]
    o["per_seed_objective_gap"] = [r.get("objective_gap") for r in rs]
    return o


def main(argv):
    quick = "--quick" in argv
    targets = V.ENUM_TARGETS[:2] if quick else V.ENUM_TARGETS
    signals = (0.0, 0.3, 1.0) if quick else SIGNALS
    budgets = (3000,) if quick else BUDGETS
    seeds = (0,) if quick else SEEDS
    # Load-shedding switches.  The box is shared and was carrying six to seven Python
    # processes from five agents while this ran; these trade grid density for completion
    # on the axis that matters (budget 30,000 is the one past the s13 1%-of-space turn).
    if "--b30" in argv:
        budgets = (30000,)
    if "--seeds2" in argv:
        seeds = (0, 1)
    base = "prior" if "--prior" in argv else "legacy"

    V.wait_for_memory(1.2, "vqe_signal")
    res = {"config": {"targets": list(targets), "signals": list(signals),
                      "budgets": list(budgets), "alphas": list(ALPHAS),
                      "seeds": list(seeds), "ansatz": ANSATZ, "shots": SHOTS,
                      "base": base, "space": "s13 qarch_enum n=9 k=4 262144"},
           "cells": {}}
    if "--resume" in argv:
        # Cells already computed are kept verbatim.  They were produced by identical code
        # on identical seeds; only the exact-distribution book-keeping changed, and that is
        # recorded per cell, so resuming does not mix two configurations of the SEARCH.
        import json
        import os
        p0 = os.path.join(V.RESULTS, f"vqe_signal_{base}.json")
        if os.path.exists(p0):
            res["cells"] = json.load(open(p0))["cells"]
            print(f"resuming: {len(res['cells'])} cells already present", flush=True)
    t0 = time.time()
    for p in targets:
        e = V.Enum(p)
        b = getattr(e, base)
        for sig in signals:
            E = V.blend_objective(b, e.rmsd, sig)
            for bud in budgets:
                key = f"{p}|{base}|{sig}|{bud}"
                if key in res["cells"]:
                    continue
                res["cells"][key] = cell(e, E, bud, seeds=seeds)
                r = res["cells"][key]
                print(f"[{time.time()-t0:7.0f}s] {key:28s} rho={r['rho_obj_rmsd']:+.3f} "
                      f"exact={r['exact']['rmsd_returned']:.3f} "
                      f"rnd={r['random']['rmsd_returned']:.3f} "
                      f"sa={r['anneal']['rmsd_returned']:.3f} "
                      f"vqe1={r['vqe_a1.0']['rmsd_returned']:.3f} "
                      f"vqe.05={r['vqe_a0.05']['rmsd_returned']:.3f}", flush=True)
                V.write(f"vqe_signal_{base}", res)
    V.write(f"vqe_signal_{base}", res)
    print(f"done in {time.time()-t0:.0f}s -> s14/results/vqe_signal_{base}.json")


if __name__ == "__main__":
    main(sys.argv[1:])
