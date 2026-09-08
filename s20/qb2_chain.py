"""SPRINT 20 / WORKSTREAM B -- ONE heavy process, the whole Q1/Q2 programme in order.

Declared scope (compute-bound; every claim carries its own n, and the prefixes are persisted in
`s20/results/qb2_config.json`):

    relax   10 targets   the RELAXATION arm -- H_AMBER = E o Relax, decomposed
    opt     12 targets   S1, the optimiser battery, B = 512 evaluations per arm and seed
    cvar     8 targets   S2, the CVaR alpha panel
    enc      8 targets   S3, theta vs (cos theta, sin theta)
    q2      20 targets   the decorrelation question

`land` (20 targets) is run separately and is expected to be on disk already.
The budget B = 512 was fixed before any optimiser arm ran and is not changed afterwards.
"""
from __future__ import annotations

import time

from s20 import qb2_lib as L
from s20 import qb2_run as R
from s20 import qb2_relax as RX
from s20 import qb2_q2 as Q2

PLAN = [("relax", 10), ("opt", 12), ("cvar", 8), ("enc", 8), ("q2", 20)]


def main():
    ss = [t["pdb"] for t in L.subset(L.N_SUBSET)]
    L.write("qb2_config", {"salt": L.SALT, "subset": ss, "B_OPT": R.B_OPT,
                           "SEEDS": list(R.SEEDS), "R_START": R.R_START, "FD_H": L.FD_H,
                           "N_HESS": R.N_HESS, "ALPHAS": list(R.ALPHAS),
                           "kinds": list(R.KINDS), "plan": PLAN,
                           "relax_n_ref": RX.N_REF, "q2_budget": Q2.BUDGET,
                           "q2_alpha": Q2.ALPHA, "q2_m_sel": Q2.M_SEL}, complete=True)
    for mode, n in PLAN:
        t0 = time.time()
        pdbs = ss[:n]
        print(f"\n########## {mode} on {n} targets ##########", flush=True)
        if mode == "relax":
            RX.run(pdbs)
        elif mode == "opt":
            R.run_opt(pdbs)
        elif mode == "cvar":
            R.run_cvar(pdbs)
        elif mode == "enc":
            R.run_opt(pdbs, arms=["spsa", "adam_fd", "lbfgs_fd", "nelder"], tag="enc",
                      kinds=("LEG", "AMB", "AMBc"), embed=True)
        elif mode == "q2":
            Q2.run(pdbs)
        print(f"########## {mode} done in {time.time()-t0:.0f}s ##########", flush=True)
    print("CHAIN COMPLETE")


if __name__ == "__main__":
    main()
