"""SPRINT 22 / WORKSTREAM A -- A2: CVaR degeneracy and the tail as a distribution.

Pre-registered in `s22/PREREG_A.md` A2, as amended by ADDENDUM 1 (dated 2026-09-07): the
beta-weighted structural term is replaced by the entropy-regularised free energy
`F_alpha,T = CVaR_alpha - T*H(p)`, because pure CVaR-over-p collapses to a single-candidate
delta (ESS=1) regardless of alpha for continuous, generically-distinct energies -- a pilot
finding recorded before this experiment ran, not discovered here.

TWO SUB-EXPERIMENTS.

  2a COLLAPSE DYNAMICS (T=0 throughout).  alpha in {1.0, 0.5, 0.15, 0.05, annealed, adaptive} x
     4 seeds x 16 targets, entangler=cnot only (the deployed default; `entangler=none` is left to
     2b, where T>0 actually exercises non-trivial support -- both entanglers collapse a T=0 run
     to the same single point, so crossing them here would not answer a new question, and that is
     stated rather than silently narrowing scope). Checkpoints at iters
     {5,20,50,100,150,200} expose HOW FAST the tail collapses, not only the converged endpoint.

  2b ENTROPY-REGULARISED DEGENERACY BREAKING.  alpha=0.15 (matches the incumbent's own m=75/500)
     x T in {0, 0.05, 0.2, 0.5} (pool-sd-relative is not needed here since H_DIST is already in
     fixed native units and T multiplies an entropy term in NATS/BITS, not H's own units -- see
     PREREG_A A2 fork 4) x entangler in {cnot, none} x 4 seeds x 16 targets, at final convergence
     (iters=200). Reports the matched-count random control AT THE SAME REALISED TAIL SIZE for
     every cell, so "does the trained tail carry information" is never confounded with "a bigger
     tail averages better regardless of content" (s21 `operator-consumes-set-mean`).

No AMBER. No OpenMM. No shared-budget bookkeeping (E is a lookup table).
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import s22.qcand_lib as QC                      # noqa: E402
from s20 import qb2_lib as L20                  # noqa: E402

ALPHAS_2A = [1.0, 0.5, 0.15, 0.05, "annealed", "adaptive"]
SEEDS = [0, 1, 2, 3]
CHECKPOINTS = [5, 20, 50, 100, 150, 200]
ITERS = 200
TS_2B = [0.0, 0.05, 0.2, 0.5]
ALPHA_2B = 0.15


def snapshot_readouts(enc, pool, snap, alpha):
    p = snap["p"]
    face = QC.exact_face(enc.E, p, alpha)
    _b, cset = QC.argmin_readout(enc, p)
    argmin_rr, _ = QC.tail_average_rmsd(pool, cset)
    cands = QC.tail_candidates(enc, face["tail"])
    r_tailavg, _ = QC.tail_average_rmsd(pool, cands)
    tail_rr = pool["rr"][cands] if len(cands) else np.array([])
    return {
        "tail_size": face["tail_size"], "ess": face["ess"],
        "entropy_bits": snap["entropy_bits"], "argmin_rr": float(argmin_rr),
        "r_tailavg": r_tailavg,
        "tail_rr_sd": float(tail_rr.std()) if len(tail_rr) > 1 else 0.0,
        "tail_rr_range": float(tail_rr.max() - tail_rr.min()) if len(tail_rr) > 1 else 0.0,
        "tail_energy_range": float(enc.scores[cands].max() - enc.scores[cands].min())
        if len(cands) > 1 else 0.0,
        "n_tail_cands": int(len(cands)),
    }


def run_2a(tg):
    rows = []
    t_start = time.time()
    for ti, t in enumerate(tg):
        pool = QC.build_pool(t["pdb"])
        enc = QC.Encoding(pool["score_dist"])
        exact = QC.classical_exact_sort(pool, pool["score_dist"], 0.15)
        cell_rows = []
        for alpha in ALPHAS_2A:
            for seed in SEEDS:
                out = QC.train(enc, alpha_schedule=alpha, seed=seed, iters=ITERS,
                               entangler="cnot", checkpoints=CHECKPOINTS)
                ck = {it: snapshot_readouts(enc, pool, snap, snap["alpha"])
                     for it, snap in out["snapshots"].items()}
                cell_rows.append({"alpha": str(alpha), "seed": seed, "checkpoints": ck,
                                  "final_alpha": out["hist_alpha"][-1],
                                  "grad_norm_mean": out["grad_norm_mean"]})
        rows.append({"pdb": t["pdb"], "n": pool["n"], "n_qubits": enc.n_qubits,
                    "R_score_argmin": exact["r_argmin"], "R_score_tailavg": exact["r_tailavg"],
                    "R_pool": float(pool["rr"].min()), "cells": cell_rows})
        print(f"[2a {ti+1}/{len(tg)}] {t['pdb']} elapsed={time.time()-t_start:.0f}s", flush=True)
    return {"experiment": "A2a_collapse_dynamics", "alphas": [str(a) for a in ALPHAS_2A],
           "seeds": SEEDS, "checkpoints": CHECKPOINTS, "iters": ITERS,
           "n_targets": len(rows), "rows": rows, "wall_time_s": time.time() - t_start}


def run_2b(tg):
    rows = []
    t_start = time.time()
    for ti, t in enumerate(tg):
        pool = QC.build_pool(t["pdb"])
        enc = QC.Encoding(pool["score_dist"])
        exact = QC.classical_exact_sort(pool, pool["score_dist"], ALPHA_2B)
        cell_rows = []
        for entangler in ("cnot", "none"):
            for T in TS_2B:
                for seed in SEEDS:
                    out = QC.train(enc, alpha_schedule=ALPHA_2B, seed=seed, iters=ITERS,
                                   entangler=entangler, T_schedule=T)
                    face = QC.exact_face(enc.E, out["p_final"], ALPHA_2B)
                    _b, cset = QC.argmin_readout(enc, out["p_final"])
                    argmin_rr, _ = QC.tail_average_rmsd(pool, cset)
                    cands = QC.tail_candidates(enc, face["tail"])
                    r_tailavg, _ = QC.tail_average_rmsd(pool, cands)
                    tail_rr = pool["rr"][cands] if len(cands) else np.array([])
                    m_match = max(1, len(cands))
                    rnd = QC.random_subset_control(pool, m_match, n_draws=24, seed=seed)
                    cell_rows.append({
                        "entangler": entangler, "T": T, "seed": seed,
                        "tail_size": face["tail_size"], "ess": face["ess"],
                        "entropy_bits": out["hist_entropy"][-1],
                        "argmin_rr": float(argmin_rr), "r_tailavg": r_tailavg,
                        "tail_rr_sd": float(tail_rr.std()) if len(tail_rr) > 1 else 0.0,
                        "n_tail_cands": int(len(cands)),
                        "random_control_matched_m": rnd,
                        "final_cvar": float(out["hist_val"][-1]),
                    })
        rows.append({"pdb": t["pdb"], "n": pool["n"], "n_qubits": enc.n_qubits,
                    "R_score_argmin": exact["r_argmin"], "R_score_tailavg": exact["r_tailavg"],
                    "R_pool": float(pool["rr"].min()), "cells": cell_rows})
        print(f"[2b {ti+1}/{len(tg)}] {t['pdb']} elapsed={time.time()-t_start:.0f}s", flush=True)
    return {"experiment": "A2b_entropy_regularised_degeneracy", "Ts": TS_2B,
           "alpha": ALPHA_2B, "seeds": SEEDS, "iters": ITERS,
           "n_targets": len(rows), "rows": rows, "wall_time_s": time.time() - t_start}


def main():
    tg = L20.subset(n=16)
    out2a = run_2a(tg)
    req_a = ["pdb", "n_qubits", "R_score_argmin", "R_score_tailavg", "cells"]
    path = QC.write("a2a_collapse", out2a, required_keys=req_a)
    print("WROTE", path)

    out2b = run_2b(tg)
    req_b = ["pdb", "n_qubits", "R_score_argmin", "R_score_tailavg", "cells"]
    path = QC.write("a2b_entropy", out2b, required_keys=req_b)
    print("WROTE", path)


if __name__ == "__main__":
    main()
