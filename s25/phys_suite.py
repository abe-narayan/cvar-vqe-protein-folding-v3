"""S25 / PHYSICS LANE -- THE 7-CONFIGURATION COMPARISON SUITE. Pre-registered: `s25/PREREG_PHYS.md`.

    python s25/phys_suite.py --phase run       # 126 targets x 7 configs x 2 normalisations
    python s25/phys_suite.py --phase analyse   # statistics; reads the per-target cache

THE SEVEN CONFIGURATIONS
========================
    1 Legacy   2 AMBER   3 Distogram   4 Legacy+Distogram
    5 AMBER+Distogram     6 Legacy+AMBER      7 Legacy+AMBER+Distogram

WHAT IS HELD IDENTICAL, AND WHY THE CODE RATHER THAN THE AUTHOR HOLDS IT
=======================================================================
Same 126 targets in the same order, same 5 pinned folds, same K=500 pool object, same
point-cloud metric, same uniform coordinate-average readout, same genuine CVaR-VQE selector at
the same (alpha=0.18, T=0.5, 3 layers, 80 iters, seed 0), same normalisation. The
configuration is a set of channel NAMES handed to `phys_lib.combine`, which is the only place
an energy is built -- so no configuration can receive a different pool, readout or
normalisation from another. Every contrast is PAIRED per target.

CVaR-VQE IS THE MATCHED SELECTOR FOR ALL SEVEN
==============================================
The classical arms below are CONTROLS reported as their own complete matched table for every
configuration, never as a substitute for a VQE row. Three of them, each answering a different
question:

    top-75        the PRODUCTION rung -- what the pipeline ships, and the anchor (F1)
    top-m_c       the SIZE-MATCHED bar at the VQE's own realised m; `rmsd_vqe - rmsd_topm` is
                  asserted ~0 by the s24 set-equality theorem, and a non-zero value is a BUG
    permuted      the energy's marginal preserved exactly, its correspondence to candidates
                  destroyed -- the control that caught s24 L16's real arm as indistinguishable
                  from noise. Run on the classical top-75 arm, where it is exact and free.
    random        16 random 75-subsets through the IDENTICAL coordinate average -- the
                  zero-information null matched in the operator's space, with max/min kept so
                  any best-of-K claim can be priced against the distribution of the extremum.

THE m-LADDER IS REPORTED, NOT HIDDEN
====================================
`RMSD_VQE(c) = RMSD_top75(c) + [RMSD_top-m_c(c) - RMSD_top75(c)] + eps_c`. The middle bracket
is the only channel the circuit has (s24 D1 SS1), and it is the one thing that could confound
"which energy" with "how many candidates were averaged" (`operator-consumes-set-mean`). Both
terms are persisted per cell so the decomposition can be read off directly.

NO MINIMISATION ANYWHERE. AMBER is an energy measurement and a feasibility diagnostic; s23 L11
measured 17 of 17 repair settings at or worse than the no-op.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402

CELLS = os.path.join(P.RESULTS, "suite_cells")
os.makedirs(CELLS, exist_ok=True)
NORMS = ("rank", "moment")


# =================================================================== one target, all cells
def run_target(pdb: str) -> Dict:
    cand, ch = P.channels(pdb)
    k = cand.k
    rng = SD.stable_rng(pdb, "randnull75", 0, P.M_PROD, salt=P.SALT)
    rnd = np.array([P.rmsd_of_set(cand, rng.choice(k, size=P.M_PROD, replace=False))
                    for _ in range(16)], float)
    # the permuted marginals, drawn ONCE per target so every configuration that permutes a
    # given channel permutes it the SAME way -- otherwise the null differs between arms.
    prng = SD.stable_rng(pdb, "permnull", 0, salt=P.SALT)
    perm = {c: prng.permutation(k) for c in P.CHANNELS}

    rows = []
    for norm in NORMS:
        for cid, name, subset in P.CONFIGS:
            E = P.combine(ch, subset, norm=norm)
            q = H.arm_vqe(cand, E, alpha=P.ALPHA, T=P.TEMP, layers=P.LAYERS,
                          iters=P.ITERS, seed=P.SEED)
            m = int(q["m"])
            gate = H.gate_set_equality(E, q["cands"], m)
            order = np.argsort(E, kind="stable")
            r_vqe = P.rmsd_of_set(cand, q["cands"])
            r_topm = P.rmsd_of_set(cand, order[:m])
            r_top75 = P.rmsd_of_set(cand, order[:P.M_PROD])
            # NULL (d): permute the PHYSICS channels only. A configuration with no physics
            # channel has no permuted arm, and the field is left as NaN rather than faked.
            phys = [c for c in subset if c != "DIS"]
            if phys:
                chp = dict(ch)
                for c in phys:
                    chp[c] = ch[c][perm[c]]
                Ep = P.combine(chp, subset, norm=norm)
                r_perm75 = P.rmsd_of_set(cand, np.argsort(Ep, kind="stable")[:P.M_PROD])
            else:
                r_perm75 = float("nan")
            rows.append(dict(
                pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(k),
                config=int(cid), config_name=name, subset="+".join(subset), norm=norm,
                basis="point_cloud", label="ACHIEVABLE",
                alpha=P.ALPHA, T=P.TEMP, seed=P.SEED,
                q_m=m, q_tail_size=int(q["tail_size"]), q_pad=int(q["n_pad_in_tail"]),
                q_cvar=float(q["cvar"]), q_entropy_bits=float(q["entropy_bits"]),
                q_ess=float(q["ess"]), q_secs=float(q["secs"]),
                gate_pass=bool(gate["pass_"]), gate_equality=bool(gate["equality"]),
                gate_n_holes=int(gate["n_holes"]),
                rmsd_vqe=r_vqe, rmsd_topm=r_topm, rmsd_top75=r_top75, rmsd_perm75=r_perm75,
                eps_vqe_minus_topm=float(r_vqe - r_topm),
                mladder_topm_minus_top75=float(r_topm - r_top75),
                rand_mean=float(rnd.mean()), rand_sd=float(rnd.std(ddof=1)),
                rand_best=float(rnd.min()), rand_worst=float(rnd.max()),
                E_mean=float(E.mean()), E_sd=float(E.std()),
                oracle_pool_best=float(np.min(cand.oracle_rr)),
                oracle_pool_mean=float(np.mean(cand.oracle_rr)),
            ))
    return dict(pdb=pdb, rows=rows)


def phase_run(pdbs: List[str]) -> int:
    t0 = time.time()
    todo = [p for p in pdbs if not os.path.exists(os.path.join(CELLS, f"{p}.json"))]
    print(f"phase=run  {len(todo)} of {len(pdbs)} targets "
          f"({len(pdbs)-len(todo)} cached)  "
          f"{len(P.CONFIGS)} configs x {len(NORMS)} normalisations", flush=True)
    import json
    for ii, pdb in enumerate(todo):
        r = run_target(pdb)
        tmp = os.path.join(CELLS, f"{pdb}.tmp{os.getpid()}")
        with open(tmp, "w") as fh:
            json.dump(r, fh)
        os.replace(tmp, os.path.join(CELLS, f"{pdb}.json"))
        el = time.time() - t0
        print(f"  [{ii+1}/{len(todo)}] {pdb}  {el:.0f}s  ({el/(ii+1):.1f}s/target, "
              f"eta {(len(todo)-ii-1)*el/(ii+1)/60:.0f}min)", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=("run",), default="run")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    pdbs = P.targets()
    if a.limit:
        pdbs = pdbs[:a.limit]
    return phase_run(pdbs)


if __name__ == "__main__":
    sys.exit(main())
