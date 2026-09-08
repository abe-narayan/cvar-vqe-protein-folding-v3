"""SPRINT 17 / QUANTUM -- T6, the two knobs Sprint 16 did not vary: SHOTS-vs-ITERATIONS at
fixed budget, and RESTARTS.  Derivation first, then measurement (sprint sections 17-18).

WHY THESE TWO AND NOT A SWEEP.  Sprint 16 already varied alpha (6 values), the learning rate
(5 values spanning 16x) and the gradient baseline, and its verdict was that the run is
UNDER-stepping: cutting the learning rate 16x costs +0.156 A on the argmin readout and
+0.676 A on the ensemble readout.  Two knobs follow from that verdict and were never turned:

  SHOTS.  At a fixed objective budget `B = shots x iters`, halving `shots` doubles the number
  of Adam steps and doubles the per-step gradient variance.  Adam normalises by a running RMS
  of the gradient, so a uniform rescaling of the gradient does not change the step length
  (this is exactly why Sprint 16's "the biased estimator is a de-facto step-size reduction"
  was mechanically impossible).  What changes is the number of steps and the direction noise.
  Given "more optimisation is better here", the DERIVED PREDICTION is that FEWER SHOTS AND
  MORE ITERATIONS is better, until the CVaR tail `alpha x shots` becomes too small to
  estimate.  At alpha = 0.25 the tail has 32 samples at shots = 128 and 512 at shots = 2048,
  so the prediction is testable across the whole ladder without entering Sprint 15's
  non-integer-`alpha N` bias regime.

  RESTARTS.  A CVaR-VQE has no restart mechanism.  Splitting the budget into `R` independent
  runs of `B/R` and pooling their draws is a pure DIVERSITY purchase: the pooled support is a
  union of `R` partially-converged supports.  The DERIVED PREDICTION is that restarts move
  the arm UP AND RIGHT on the (member error, diversity) plane -- exactly like raising a
  temperature -- and therefore cannot produce a point a temperature ladder does not already
  supply.  That is a falsifiable prediction, not a hope: if restarts produced a point outside
  the classical frontier, alpha-as-temperature would be incomplete.

SUCCESS CRITERION.  Same as P1: a point outside the classical frontier measured in
`q_pareto`.  This module reports the (M, D) points; `q_report` scores them against the same
classical set.

RUN:  python -m s17.q_sched run
"""
from __future__ import annotations

import sys
import time

import numpy as np

from s14 import vqe_run as R
from s15 import seed as SD
from s17 import q_lib as L
from s17 import q_pareto as P

SEEDS = (0,)
ALPHA = 0.25
SHOTS_LADDER = (128, 512, 2048)
RESTARTS = (4,)

#: THE ENTANGLEMENT CONTROL, and the most important arm in this module.
#: `mps2fn` is the identical circuit with the CNOTs REMOVED -- a fully factorised Bernoulli
#: model over the qubits, trained by the IDENTICAL CVaR score-function estimator, the
#: IDENTICAL Adam optimiser and the IDENTICAL 8,192-evaluation budget.  Sampling and
#: gradients of a product Bernoulli model are trivially classical, so this arm IS a classical
#: variational sampler; it is counted in the CLASSICAL attainable set in `q_report`.  If it
#: matches `mps2f`, the entanglement in the ansatz is contributing nothing and the "quantum"
#: arm is a classical variational sampler wearing a circuit.  Sprint 16 never ran it.
PRODUCT_ALPHAS = (0.05, 0.25, 1.00)


def run_cell(pdb, seed):
    ins = L.inst(pdb)
    E = ins.hamil()
    rng = SD.stable_rng(pdb, seed, "readout", salt=L.SALT)
    res = {"_meta": {"pdb": pdb, "seed": seed, "fold": ins.fold, "n": ins.n}}
    for al in PRODUCT_ALPHAS:
        c = L.Cost()
        v = R.run(E, ins.n_qubits, al, L.BUDGET, shots=L.SHOTS, ansatz="mps2fn", seed=seed,
                  rmsd=ins.rmsd, bits_per_res=ins.bits_per_res, exact_dist=False,
                  keep_seen=True, lr=0.15, baseline="const")
        seen = np.asarray(v["seen"], np.int64)
        r, _ = P.arm_readout(ins, E, seen, rng,
                             c.stop(v["evals"], circuit_samples=v["evals"],
                                    grad_passes=v["iters"],
                                    distinct=int(np.unique(seen).size)))
        r["iters"] = int(v["iters"])
        res[f"product_a{al}"] = r
    for sh in SHOTS_LADDER:
        c = L.Cost()
        v = R.run(E, ins.n_qubits, ALPHA, L.BUDGET, shots=sh, ansatz=L.ANSATZ, seed=seed,
                  rmsd=ins.rmsd, bits_per_res=ins.bits_per_res, exact_dist=False,
                  keep_seen=True, lr=0.15, baseline="const")
        seen = np.asarray(v["seen"], np.int64)
        r, _ = P.arm_readout(ins, E, seen, rng,
                             c.stop(v["evals"], circuit_samples=v["evals"],
                                    grad_passes=v["iters"],
                                    distinct=int(np.unique(seen).size)))
        r["iters"] = int(v["iters"])
        r["tail_samples"] = float(ALPHA * sh)
        res[f"shots{sh}"] = r
    for rr in RESTARTS:
        c = L.Cost()
        got, it = [], 0
        for j in range(rr):
            v = R.run(E, ins.n_qubits, ALPHA, L.BUDGET // rr, shots=min(512, L.BUDGET // rr),
                      ansatz=L.ANSATZ, seed=seed * 100 + j, rmsd=ins.rmsd,
                      bits_per_res=ins.bits_per_res, exact_dist=False, keep_seen=True,
                      lr=0.15, baseline="const")
            got.append(np.asarray(v["seen"], np.int64))
            it += int(v["iters"])
        seen = np.concatenate(got)
        r, _ = P.arm_readout(ins, E, seen, rng,
                             c.stop(seen.size, circuit_samples=seen.size, grad_passes=it,
                                    distinct=int(np.unique(seen).size)))
        r["iters"] = it
        r["restarts"] = rr
        res[f"restart{rr}"] = r
    return res


def run(targets=L.TARGETS19, seeds=SEEDS):
    done = L.ck_load("sched")
    for sd in seeds:
        for pdb in targets:
            key = f"{pdb}_{sd}"
            if key in done:
                print(f"  skip {key}", flush=True)
                continue
            L.gate(key)
            t0 = time.time()
            r = run_cell(pdb, sd)
            L.ck("sched", key, r)
            done[key] = r
            print(f"  {key} {time.time()-t0:.0f}s  shots128 M/D "
                  f"{r['shots128']['M75']:.3f}/{r['shots128']['D75']:.3f}  shots2048 "
                  f"{r['shots2048']['M75']:.3f}/{r['shots2048']['D75']:.3f}  restart4 "
                  f"{r['restart4']['M75']:.3f}/{r['restart4']['D75']:.3f}  prod0.25 "
                  f"{r['product_a0.25']['M75']:.3f}/{r['product_a0.25']['D75']:.3f}",
                  flush=True)


if __name__ == "__main__":
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        run(sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS19)
