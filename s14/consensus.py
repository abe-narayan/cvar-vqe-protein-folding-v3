"""SPRINT 14, coordinator -- select a SET with the structural Hamiltonian, then average it.

THIS IS THE ARCHITECTURE THE SPRINT'S OWN MEASUREMENTS POINT AT.  Three results combine.

  C2  Averaging in COORDINATE space beats committing to one torsion vector.  The retrieval
      pool's position-specific torsion distribution is the best torsion channel in this
      project (phi 33.6 deg, psi 59.2 deg) and still emits 4.072 A, while the pipeline that
      coordinate-averages the identical windows emits 3.204 A.  Torsion errors compound
      through chain integration; coordinate errors cancel in place.

  C7  A native-free STRUCTURAL Hamiltonian -- retrieval torsion prior plus the shipped
      distogram's Bayes risk -- orders the space eight times better than any physical
      energy (in-decile rank correlation +0.370 against Legacy's +0.043).  But its ARGMIN
      emits 3.575 A while the same 4,000 sampled configurations contain a 1.683 A structure.
      A 1.9 A selection gap survives an eightfold gain in ranking.

  S12 The terminal operator consumes the set MEAN, not the set BEST:
      `d_out = 1.16*d_set_mean + 0.04*d_set_best`, R^2 0.893.  A perfect rank-1 decision is
      worth -1.743 A through argmin and only -0.029 A through an m=75 average.

Put together they say the argmin is the wrong readout.  A well-ordered objective should be
used to SELECT A SET whose members are individually mediocre but whose errors are
independent, and that set should then be aggregated in coordinate space where the errors
cancel.  That is exactly what the production pipeline does to retrieved windows, applied
instead to configurations the Hamiltonian chose.

It is also the classical limit of the quantum architecture the brief asks for: a variational
state IS a distribution over configurations, CVaR shapes its low-energy tail, and measuring
it B times gives precisely a low-energy SET.  So this file measures what a VQE would have to
beat, with no VQE in it -- which is the honest control the brief's section 23 demands, and
the number any later quantum claim must be differenced against.

The set size m is swept because Sprint 12 established that the optimal set size shrinks
(500 -> 75 -> 20 -> 3-5) as the objective improves, so m* is itself a readout of objective
quality.

NATIVE-FREE throughout; the native is read only to score.

Run:
    python -m s14.consensus
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

SET_SIZES = [1, 3, 5, 10, 20, 50, 100, 200, 500]
N_SAMPLE = 4000
W = 0.25                       # best in-decile rank correlation from s14/hamil.py


def run_target(t, w=W, n_sample=N_SAMPLE, sizes=SET_SIZES, seed=0, project=True):
    T = H.Terms(t["pdb"], t["seq"], t["n"], t["fold"])
    rng = np.random.default_rng(hash((t["pdb"], seed)) % (2 ** 32))

    S1 = H._sample(T.P, n_sample // 2, rng)
    S2 = rng.integers(0, T.k, size=(n_sample - n_sample // 2, T.n))
    S = np.vstack([S1, S2])

    ep, ed = H._z(T.e_prior(S)), H._z(T.e_disto(S))
    E = (1.0 - w) * ep + w * ed
    order = np.argsort(E)
    rm = T.ORACLE_rmsd(S)                                    # ORACLE, post hoc only

    u = I.load_univ(t["pdb"])
    nat = u["nat_ca"]
    rows = np.arange(T.n)
    out = {"pdb": t["pdb"], "fold": int(t["fold"]),
           "oracle_best_sampled": float(rm.min()),
           "argmin_rmsd": float(rm[int(order[0])]), "sets": {}}

    for m in sizes:
        if m > len(S):
            continue
        idx = order[:m]
        CA = I.build_ca(T.sp.PHI[rows[None, :], S[idx]],
                        T.sp.PSI[rows[None, :], S[idx]])
        CA = np.asarray(CA, float).reshape(m, T.n, 3)
        C, _ = I.coordinate_average(CA) if m > 1 else (CA[0], 0)
        rec = {"set_mean_rmsd": float(rm[idx].mean()),
               "set_best_rmsd": float(rm[idx].min()),
               "consensus_raw": float(I.ca_rmsd(C, nat))}
        if project:
            rec["consensus_projected"] = float(
                I.ca_rmsd(I.project(C, t["seq"], int(t["fold"]))["fit_ca"], nat))
        out["sets"][str(m)] = rec
    return out


def run(targets=None, w=W, sizes=SET_SIZES, project=True):
    tg = targets if targets is not None else I.targets()
    rows = [run_target(t, w, N_SAMPLE, sizes, project=project) for t in tg]

    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    fail = np.isin(pdbs, I.FAIL18)
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)

    key = "consensus_projected" if project else "consensus_raw"
    agg = {}
    for m in sizes:
        k = str(m)
        got = [r for r in rows if k in r["sets"]]
        if not got:
            continue
        v = np.asarray([r["sets"][k][key] for r in got], float)
        vr = np.asarray([r["sets"][k]["consensus_raw"] for r in got], float)
        sm = np.asarray([r["sets"][k]["set_mean_rmsd"] for r in got], float)
        sb = np.asarray([r["sets"][k]["set_best_rmsd"] for r in got], float)
        f2 = np.isin([r["pdb"] for r in got], I.FAIL18)
        agg[k] = {"n": len(got), **I.summary(v),
                  "FAIL18": float(v[f2].mean()) if f2.any() else None,
                  "consensus_raw_mean": float(vr.mean()),
                  "set_mean_rmsd": float(sm.mean()),
                  "set_best_rmsd": float(sb.mean()),
                  "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}

    ks = [k for k in map(str, sizes) if k in agg]
    ys = [agg[k]["mean"] for k in ks]
    best = int(np.argmin(ys))

    out = {"n": len(rows), "w_distogram": w, "n_sample": N_SAMPLE,
           "projected": project, "sizes": sizes, "aggregate": agg,
           "best_set_size": int(ks[best]), "best_mean": float(ys[best]),
           "incumbent": float(ref.mean()),
           "oracle_best_sampled": float(np.mean([r["oracle_best_sampled"] for r in rows])),
           "argmin_only": float(np.mean([r["argmin_rmsd"] for r in rows])),
           "per_target": rows}
    with open(os.path.join(RESULTS, "consensus.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_consensus", out, n_expected=len(rows))

    print(f"incumbent {ref.mean():.3f}    argmin only {out['argmin_only']:.3f}"
          f"    ORACLE best sampled {out['oracle_best_sampled']:.3f}\n")
    print(f"{'m':>5}{'set mean':>10}{'set best':>10}{'consensus':>11}{'projected':>11}"
          f"{'<2A':>7}{'FAIL18':>9}{'vs incumbent':>26}")
    for k in ks:
        a = agg[k]
        f18 = f"{a['FAIL18']:.3f}" if a["FAIL18"] is not None else "  -  "
        v = a["vs_incumbent"]
        print(f"{k:>5}{a['set_mean_rmsd']:>10.3f}{a['set_best_rmsd']:>10.3f}"
              f"{a['consensus_raw_mean']:>11.3f}{a['mean']:>11.3f}"
              f"{a['frac_under_2.0']:>7.2f}{f18:>9}"
              f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print(f"\nbest set size m* = {out['best_set_size']} at {out['best_mean']:.3f} A")
    return out


if __name__ == "__main__":
    run()
