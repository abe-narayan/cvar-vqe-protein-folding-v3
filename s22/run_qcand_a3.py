"""SPRINT 22 / WORKSTREAM A -- A3: multi-stage CVaR-VQE on the candidate-basis encoding.

Pre-registered in `s22/PREREG_A.md` A3. One circuit, warm-started parameters, TWO stages with
stated purposes:

  Stage 1 (purpose: enter a useful structural region).  H_Legacy, 60 iterations. Motivated by
  s21 L3's finding that Legacy's own tail is compact/pool-typical -- a plausible "narrow the
  field to plausible structures first" role. This is a DIFFERENT substrate from s21 C7'/C7''
  (continuous torsion, where staging was harmful): tested here, not assumed to inherit that
  result.

  Stage 2 (purpose: discriminate candidates). H_DIST, 140 iterations, continuing the SAME
  circuit parameters and Adam moments from stage 1.

  FINAL selection: genuine CVaR-VQE readout on H_DIST (argmin and tail-average), matching the
  single-stage control's own Hamiltonian so Pillar 1 (CVaR-VQE is the selector, on the
  Hamiltonian actually being scored) is not silently switched to Legacy's answer.

CONTROL: single-stage H_DIST-only training for the SAME total 200 iterations, from the SAME
initial theta per seed (paired design) -- never a free compute allowance for staging.

PROMOTION CRITERION (pre-registered): staged beats single-stage by more than the comparison's own
MDE on a majority of targets, AND the argmin readout still ties the exact H_DIST global optimum
(checked explicitly per cell -- a staged run whose argmin drifts to Legacy's optimum would be
reporting Legacy's answer under a DIST label).

REGISTERED EXPECTATION (before running): staging does NOT help, by analogy with s21 C7'/C7'' on
the continuous encoding -- stated so a negative result is not read as a surprise.
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
from core import quantum as Q                   # noqa: E402

ALPHA = 0.15
SEEDS = [0, 1, 2, 3]
STAGE1_ITERS = 60
STAGE2_ITERS = 140
TOTAL_ITERS = STAGE1_ITERS + STAGE2_ITERS


def readouts(enc, pool, p):
    face = QC.exact_face(enc.E, p, ALPHA)
    _b, cset = QC.argmin_readout(enc, p)
    argmin_rr, _ = QC.tail_average_rmsd(pool, cset)
    cands = QC.tail_candidates(enc, face["tail"])
    r_tailavg, _ = QC.tail_average_rmsd(pool, cands)
    return {"argmin_set": sorted(int(x) for x in cset), "argmin_rr": float(argmin_rr),
           "r_tailavg": r_tailavg, "tail_size": face["tail_size"], "ess": face["ess"]}


def run_seed(pool, enc_dist, enc_leg, seed):
    n = enc_dist.n_qubits
    import numpy as _np
    rng = _np.random.default_rng(seed)
    theta_init = (_np.pi / 2.0) + rng.normal(0.0, 0.8, Q.MPSAnsatz(n, layers=2, final_ry=True,
                                                                   entangler="cnot").n_params())

    # ---- single-stage control: H_DIST only, TOTAL_ITERS, from theta_init
    single = QC.train(enc_dist, alpha_schedule=ALPHA, seed=seed, iters=TOTAL_ITERS,
                      entangler="cnot", theta_init=theta_init.copy())
    r_single = readouts(enc_dist, pool, single["p_final"])

    # ---- staged: Legacy (STAGE1_ITERS) -> Dist (STAGE2_ITERS), same params, same Adam moments
    stage1 = QC.train(enc_leg, alpha_schedule=ALPHA, seed=seed, iters=STAGE1_ITERS,
                      entangler="cnot", theta_init=theta_init.copy())
    stage2 = QC.train(enc_dist, alpha_schedule=ALPHA, seed=seed, iters=STAGE2_ITERS,
                      entangler="cnot", theta_init=stage1["theta"].copy(),
                      opt_state=stage1["opt"])
    r_staged = readouts(enc_dist, pool, stage2["p_final"])

    exact = QC.classical_exact_sort(pool, pool["score_dist"], ALPHA)
    exact_set = set(int(x) for x in exact["argmin_set"])
    argmin_tied_single = (set(r_single["argmin_set"]) == exact_set)
    argmin_tied_staged = (set(r_staged["argmin_set"]) == exact_set)

    return {
        "seed": seed,
        "single": r_single, "staged": r_staged,
        "argmin_tied_single": bool(argmin_tied_single),
        "argmin_tied_staged": bool(argmin_tied_staged),
        "stage1_final_legacy_cvar": float(stage1["hist_val"][-1]),
        "single_final_dist_cvar": float(single["hist_val"][-1]),
        "staged_final_dist_cvar": float(stage2["hist_val"][-1]),
    }


def main():
    tg = L20.subset(n=16)
    rows = []
    t_start = time.time()
    for ti, t in enumerate(tg):
        pool = QC.build_pool(t["pdb"])
        enc_dist = QC.Encoding(pool["score_dist"])
        enc_leg = QC.Encoding(pool["score_leg"])
        exact = QC.classical_exact_sort(pool, pool["score_dist"], ALPHA)
        seed_rows = [run_seed(pool, enc_dist, enc_leg, s) for s in SEEDS]
        rows.append({"pdb": t["pdb"], "n": pool["n"], "n_qubits": enc_dist.n_qubits,
                    "R_score_argmin": exact["r_argmin"], "R_score_tailavg": exact["r_tailavg"],
                    "R_pool": float(pool["rr"].min()), "seeds": seed_rows})
        print(f"[A3 {ti+1}/{len(tg)}] {t['pdb']} elapsed={time.time()-t_start:.0f}s", flush=True)

    # ---- paired analysis: staged tail-average RMSD vs single-stage, target-level mean over seeds
    staged_tv, single_tv = [], []
    staged_arg, single_arg = [], []
    n_tied_single = n_tied_staged = n_cells = 0
    for r in rows:
        staged_tv.append(np.mean([s["staged"]["r_tailavg"] for s in r["seeds"]]))
        single_tv.append(np.mean([s["single"]["r_tailavg"] for s in r["seeds"]]))
        staged_arg.append(np.mean([s["staged"]["argmin_rr"] for s in r["seeds"]]))
        single_arg.append(np.mean([s["single"]["argmin_rr"] for s in r["seeds"]]))
        for s in r["seeds"]:
            n_cells += 1
            n_tied_single += int(s["argmin_tied_single"])
            n_tied_staged += int(s["argmin_tied_staged"])

    ci_tailavg = QC.paired_ci(staged_tv, single_tv)
    ci_argmin = QC.paired_ci(staged_arg, single_arg)

    out = {
        "experiment": "A3_multistage",
        "alpha": ALPHA, "stage1_iters": STAGE1_ITERS, "stage2_iters": STAGE2_ITERS,
        "total_iters": TOTAL_ITERS, "seeds": SEEDS, "n_targets": len(rows),
        "rows": rows,
        "primary_staged_minus_single_tailavg": ci_tailavg,
        "secondary_staged_minus_single_argmin": ci_argmin,
        "argmin_tied_exact_optimum": {
            "single": f"{n_tied_single}/{n_cells}", "staged": f"{n_tied_staged}/{n_cells}"},
        "wall_time_s": time.time() - t_start,
    }
    req = ["pdb", "n_qubits", "R_score_argmin", "R_score_tailavg", "seeds"]
    path = QC.write("a3_multistage", out, required_keys=req)
    print("WROTE", path)
    print("PRIMARY staged - single (tail-average RMSD):", ci_tailavg)
    print("argmin tied to exact optimum:", out["argmin_tied_exact_optimum"])


if __name__ == "__main__":
    main()
