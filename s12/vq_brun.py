"""s12 QUANTUM-ROLE -- PROBLEM B: subset selection for the terminal operator.

The aggregation agent dumped, for all 126 targets, a 16-binary-variable instance
(`s12/results/agg_subset_instances.json`) with a certified exhaustive optimum:

    H(x) = x'Mx/k^2 - 2 b'x/k + const ,   k = sum(x) >= 2 ,  x in {0,1}^16

That is NOT a QUBO -- the 1/k and 1/k^2 make it a ratio of quadratics -- so the honest
statement of the encoding is: 16 qubits, one per candidate, whole register feasible except
the k < 2 states, which are given a penalty equal to the finite maximum plus the span.  The
`k`-normalisation is what makes the objective non-separable in the way that matters here
(the risk of an average is not the average of risks), and it is kept exactly.

Every arm is run against the SAME `budget.BudgetedEnergyModel`, so quantum and classical
compete on unique-bitstring evaluations.  Emission (uniform coordinate average of the
selected medoid-superposed candidates) is reproduced from the aggregation agent's own cache
and validated against its dumped RMSDs before any number here is used.

Usage: python -m s12.vq_brun [n_targets] [budgets]
"""
from __future__ import annotations
import os, sys, json, time, warnings
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_common as AC
from s12 import vq_classical as C
from s12 import vq_quantum as Qm

INST = os.path.join(ROOT, "s12", "results", "agg_subset_instances.json")


def instances():
    with open(INST) as fh:
        return json.load(fh)["instances"]


def emit(pdb, hypo_idx, sel):
    """Uniform coordinate average of the selected hypotheses, medoid-superposed.

    `hypo_idx` indexes the shipped top-75; the top-75's own medoid is the superposition
    reference, exactly as the production operator uses it.
    """
    d = AC.load(pdb)
    W75 = d["Wp"][d["sub"].astype(int)]
    ref = W75[int(np.argmin(I.pairwise_rmsd(W75).mean(1)))]
    Wh = W75[np.asarray(hypo_idx, int)]
    sel = np.asarray(sel, int)
    A = I.superpose_batch(Wh[sel], ref)
    return A.mean(0), d["nat"]


def validate_emission(inst_list, n=8):
    """Reproduce the aggregation agent's dumped RMSDs before trusting the operator."""
    out = []
    for q in inst_list[:n]:
        hy = np.asarray(q["hypo_idx_into_top75"], int)
        for field in ("uniform_all_hypotheses", "exhaustive_surrogate", "ORACLE_best_subset"):
            e = q[field]
            S = list(range(len(hy))) if field == "uniform_all_hypotheses" else e["S"]
            X, nat = emit(q["pdb"], hy, S)
            out.append(dict(pdb=q["pdb"], field=field, dumped=float(e["rmsd"]),
                            mine=float(I.ca_rmsd(X, nat))))
    d = np.array([abs(r["dumped"] - r["mine"]) for r in out])
    return dict(n=len(out), max_abs_diff=float(d.max()), mean_abs_diff=float(d.mean()),
                rows=out[:6])


def energy_vector(M, b, const):
    """Full 2^16 energy vector with the k<2 states penalised (not left at +inf)."""
    N = len(b)
    X = ((np.arange(1 << N)[:, None] >> np.arange(N - 1, -1, -1)[None, :]) & 1).astype(float)
    E = C.b_energy(np.asarray(M, float), np.asarray(b, float), float(const), X)
    fin = np.isfinite(E)
    pen = float(E[fin].max() + (E[fin].max() - E[fin].min()))
    E = np.where(fin, E, pen)
    return E, X, N


def run(ntarg=None, budgets=(256, 1024, 4096, 16384), layers=3, alpha=0.25,
        shots=16, restarts=1, seed=0, out="vq_stage3_B"):
    from core import quantum as Q
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ins = instances()
    val = validate_emission(ins, 8)
    print("emission validation:", json.dumps({k: v for k, v in val.items() if k != "rows"}), flush=True)
    if ntarg:
        ins = ins[:: max(1, len(ins) // ntarg)][:ntarg]
    rows = []
    t00 = time.time()
    for q, inst in enumerate(ins):
        while I.free_gb() < 1.2:
            print('  waiting on memory', flush=True); time.sleep(20)
        M = np.asarray(inst["M"], float); bv = np.asarray(inst["b"], float)
        const = float(inst["const"])
        E, X, N = energy_vector(M, bv, const)
        hy = np.asarray(inst["hypo_idx_into_top75"], int)
        opt = float(E.min())
        opt_state = int(np.argmin(E))
        # emitted RMSD for every one of the 65536 subsets (evaluation only)
        d = AC.load(inst["pdb"])
        W75 = d["Wp"][d["sub"].astype(int)]
        ref = W75[int(np.argmin(I.pairwise_rmsd(W75).mean(1)))]
        Ah = I.superpose_batch(W75[hy], ref)                    # (16, n, 3)
        cnt = X.sum(1)
        SUM = X @ Ah.reshape(len(hy), -1)
        AVG = (SUM / np.maximum(cnt, 1)[:, None]).reshape(len(X), *Ah.shape[1:])
        rr = I.kabsch_rmsd_batch(AVG, d["nat"])
        rr[cnt < 2] = np.nan
        row = dict(pdb=inst["pdb"], n=int(inst["n"]), fold=int(inst["fold"]),
                   n_qubits=N, dim=int(len(E)), opt=opt,
                   opt_rmsd=float(rr[opt_state]),
                   oracle_best_rmsd=float(np.nanmin(rr)),
                   avg75_rmsd=float(inst["avg75_rmsd"]),
                   uniform16_rmsd=float(inst["uniform_all_hypotheses"]["rmsd"]),
                   dumped_exhaustive_H=float(inst["exhaustive_surrogate"]["H"]),
                   rmsd_percentile_of_optimum=float(np.nanmean(rr < rr[opt_state])),
                   corr_E_rmsd=float(np.corrcoef(E[np.isfinite(rr)], rr[np.isfinite(rr)])[0, 1]),
                   budgets={})
        for B in budgets:
            entry = {}
            ham = Qm.qubo_hamiltonian(E, N, eval_budget=B)
            t0 = time.perf_counter()
            try:
                res = Q.run_global_cvar_vqe(ham, layers=layers, alpha=alpha, shots=shots,
                                            restarts=restarts, seed=seed, optimizer="SPSA",
                                            device="lightning.qubit",
                                            final_shots=max(256, shots * 4), verbose=False)
                bs = res["best_seen_bitstring"]
                entry["vqe"] = dict(best_energy=float(res["best_seen_energy"]),
                                    entropy_bits=float(res["distribution_entropy_bits"]),
                                    spsa_iters=int(res["n_spsa_iterations_total"]),
                                    n_energy_evaluations=int(res["n_energy_evaluations"]),
                                    rmsd=float(rr[int(bs, 2)]) if bs else None,
                                    wall=float(time.perf_counter() - t0))
            except Exception as e:                              # noqa: BLE001
                entry["vqe"] = dict(error=str(e))
            for nm, fn in (("random", Qm.budgeted_random),
                           ("anneal", Qm.budgeted_anneal),
                           ("greedy_ls", Qm.budgeted_greedy_ls)):
                hm = Qm.qubo_hamiltonian(E, N, eval_budget=B)
                r = fn(hm, seed=seed)
                r["rmsd"] = float(rr[int(r["best_bitstring"], 2)]) if r["best_bitstring"] else None
                entry[nm] = r
            row["budgets"][str(B)] = entry
        rows.append(row)
        if q % 5 == 0:
            print(f"  [B] {q+1}/{len(ins)} {inst['pdb']} {time.time()-t00:.0f}s "
                  f"free={I.free_gb():.2f}", flush=True)
            I.write(out, dict(problem="B_subset", validation=val, alpha=alpha, shots=shots,
                              layers=layers, restarts=restarts, rows=rows))
    I.write(out, dict(problem="B_subset", validation=val, alpha=alpha, shots=shots,
                      layers=layers, restarts=restarts, rows=rows))
    return rows


if __name__ == "__main__":
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else None
    bg = tuple(int(x) for x in sys.argv[2].split(",")) if len(sys.argv) > 2 else (256, 1024, 4096, 16384)
    run(nt, bg)
