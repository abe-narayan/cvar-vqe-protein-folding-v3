"""SPRINT 14 / VQE -- the coordinate-averaged CVaR tail.

The coordinator's finding: on `s14/hamil`, ranking 4,000 configurations and
COORDINATE-AVERAGING the top 200 recovers 0.29 A over the argmin (3.608 -> 3.314) even
though no individual set member improves.  Coordinate-space aggregation is the one operator
that recovers anything at this objective quality.

That is exactly the operator CVaR is built to shape.  Every other VQE readout tested in this
sprint decodes to a SINGLE bitstring, which throws away the distribution -- and a
distribution is the only thing a variational state has that a classical search does not.  So
this is the honest quantum analogue of the coordinator's classical control, and the question
is sharp:

    does a CVaR-shaped tail coordinate-average BETTER than a classical top-m
    coordinate-average at MATCHED objective-evaluation budget?

Arms, all at the same budget, all reported with the same aggregation:
    classical_topm   rank the sampled configurations by the objective, average the top m.
    vqe_bestm        the VQE's m lowest-energy SAMPLED configurations, averaged.
    vqe_tail_p       the alpha-tail of the FINAL VARIATIONAL DISTRIBUTION, averaged with
                     weights p(x).  This one uses the learned state, not the sample, and is
                     the only arm in the sprint for which "the VQE did something" is even a
                     well-posed claim.
    vqe_full_p       the whole final distribution, averaged with weights p(x).

The reference to beat is the argmin of the same arm, and the control is the classical one.

    python -m s14.vqe_tailavg
"""
from __future__ import annotations

import numpy as np

from core import geometry as geo
from core import quantum as Q
from s12 import instrument as I
from s14 import vqe_lib as V
from s14 import vqe_run as R

ALPHAS = (1.0, 0.25, 0.1, 0.05)
MS = (1, 5, 20, 75, 200, 500)


def ca_of(sp, idx, e: V.Enum):
    """CA traces for a list of configuration indices."""
    return sp.ca(e.states(np.asarray(idx, np.int64)))


def coord_average(cas, weights, ref_i=0):
    """Kabsch-superpose every member onto `ref_i` and take the weighted mean."""
    ref = np.asarray(cas[ref_i], float)
    w = np.asarray(weights, float)
    w = w / max(w.sum(), 1e-300)
    stack = np.stack([geo.kabsch_superpose(np.asarray(c, float), ref) for c in cas])
    return (stack * w[:, None, None]).sum(0)


def _rmsd(sp, ca):
    return float(I.ca_rmsd(np.asarray(ca, float), sp.nat_ca))       # post-hoc only


def one_target(pdb, signal, budget=20480, shots=512, seeds=(0, 1, 2)):
    from s13.qarch_lib import Space
    e = V.Enum(pdb)
    sp = Space(pdb, k=4)
    E = V.blend_objective(e.legacy, e.rmsd, signal)
    out = {"pdb": pdb, "signal": signal, "rho": V.spearman(E, e.rmsd),
           "space_best": float(e.rmsd.min()),
           "exact_rmsd": float(e.rmsd[int(np.argmin(E))])}

    # ---- classical control: uniform sampling at the same budget, rank, average top m
    cl = {m: [] for m in MS}
    for s in seeds:
        rng = np.random.default_rng(900 + s)
        idx = rng.integers(0, e.N, budget)
        order = idx[np.argsort(E[idx])]
        for m in MS:
            sel = order[:m]
            cas = ca_of(sp, sel, e)
            cl[m].append(_rmsd(sp, coord_average(cas, np.ones(len(sel)))))
    out["classical_topm"] = {str(m): float(np.mean(v)) for m, v in cl.items()}

    # ---- VQE arms
    for a in ALPHAS:
        bm = {m: [] for m in MS}
        tp, fp, arg = [], [], []
        for s in seeds:
            r = R.run(E, e.n_qubits, a, budget=budget, shots=shots, ansatz="mps2f",
                      seed=s, rmsd=e.rmsd, keep_seen=True)
            # (i) the m lowest-energy configurations the VQE ACTUALLY SAMPLED.  Uses only
            #     energies already paid for, so it is budget-matched to the classical arm.
            seen = np.unique(np.asarray(r["seen"], np.int64))
            order = seen[np.argsort(E[seen])]
            for m in MS:
                sel_m = order[:m]
                bm[m].append(_rmsd(sp, coord_average(ca_of(sp, sel_m, e),
                                                     np.ones(len(sel_m)))))
            # (ii) the alpha-tail of the FINAL VARIATIONAL DISTRIBUTION, weighted by p(x).
            #      This is the only readout that uses the learned state rather than the
            #      sample.  Computing the exact tail over all 4^n states reads the whole E
            #      array, so it is an ORACLE-BUDGET diagnostic, labelled as such: it prices
            #      the ceiling of the operator, not a deployable arm.
            an = R.make_ansatz("mps2f", e.n_qubits)
            p = np.maximum(an.probs(np.array(r["theta"])), 0.0)
            p = p / p.sum()
            _, _, mass = Q.cvar_from_probs(E, p, a)
            sel = np.flatnonzero(mass > 0)
            if sel.size == 0:
                sel = np.array([int(np.argmax(p))])
            if sel.size > 4000:                       # cap the superposition cost
                sel = sel[np.argsort(mass[sel])[::-1][:4000]]
            tp.append(_rmsd(sp, coord_average(ca_of(sp, sel, e), mass[sel])))
            # (iii) the whole final distribution, its top 2,000 states by mass.
            keep = np.argsort(p)[::-1][:2000]
            fp.append(_rmsd(sp, coord_average(ca_of(sp, keep, e), p[keep])))
            arg.append(float(e.rmsd[r["best_i"]]))
        out[f"vqe_a{a}"] = {
            "argmin_rmsd": float(np.mean(arg)),
            "tail_p_avg": float(np.mean(tp)),
            "full_p_avg": float(np.mean(fp)),
            "bestm": {str(m): float(np.mean(v)) for m, v in bm.items()},
            "tail_size": int(sel.size)}
    return out


def main():
    V.wait_for_memory(1.2, "vqe_tailavg")
    rows = []
    for sig in (0.1, 0.3):
        for pdb in V.ENUM_TARGETS:
            r = one_target(pdb, sig)
            rows.append(r)
            print(f"{pdb} sig={sig} rho={r['rho']:+.3f} exact={r['exact_rmsd']:.3f} "
                  f"| classical top-200 {r['classical_topm']['200']:.3f} "
                  f"| a=0.1 argmin {r['vqe_a0.1']['argmin_rmsd']:.3f} "
                  f"tail-p {r['vqe_a0.1']['tail_p_avg']:.3f} "
                  f"full-p {r['vqe_a0.1']['full_p_avg']:.3f}", flush=True)
            V.write("vqe_tailavg", {"rows": rows})
    report(rows)
    V.write("vqe_tailavg", {"rows": rows})


def report(rows):
    print()
    print("=" * 100)
    print("COORDINATE-AVERAGED READOUTS -- mean CA-RMSD over the nine enumerated targets")
    print("=" * 100)
    for sig in sorted({r["signal"] for r in rows}):
        rs = [r for r in rows if r["signal"] == sig]
        print(f"\n--- signal {sig}  (mean rho {np.mean([r['rho'] for r in rs]):+.3f}, "
              f"certified optimum {np.mean([r['exact_rmsd'] for r in rs]):.3f} A, "
              f"space best {np.mean([r['space_best'] for r in rs]):.3f}) ---")
        print(f"  {'arm':22s} " + "".join(f"{f'm={m}':>9s}" for m in MS)
              + f"{'tail-p':>9s}{'full-p':>9s}{'argmin':>9s}")
        v = [np.mean([r["classical_topm"][str(m)] for r in rs]) for m in MS]
        print(f"  {'classical top-m':22s} " + "".join(f"{x:9.3f}" for x in v)
              + f"{'-':>9s}{'-':>9s}{v[0]:9.3f}")
        for a in ALPHAS:
            k = f"vqe_a{a}"
            v = [np.mean([r[k]["bestm"][str(m)] for r in rs]) for m in MS]
            print(f"  {'CVaR-VQE a=' + str(a):22s} " + "".join(f"{x:9.3f}" for x in v)
                  + f"{np.mean([r[k]['tail_p_avg'] for r in rs]):9.3f}"
                  + f"{np.mean([r[k]['full_p_avg'] for r in rs]):9.3f}"
                  + f"{np.mean([r[k]['argmin_rmsd'] for r in rs]):9.3f}")
        # the paired test that matters
        best_cl = min(MS, key=lambda m: np.mean([r["classical_topm"][str(m)] for r in rs]))
        print(f"\n  best classical aggregation: m={best_cl} at "
              f"{np.mean([r['classical_topm'][str(best_cl)] for r in rs]):.3f} A")
        for a in ALPHAS:
            k = f"vqe_a{a}"
            for name, vals in (("tail-p", [r[k]["tail_p_avg"] for r in rs]),
                               ("full-p", [r[k]["full_p_avg"] for r in rs])):
                p = I.paired(np.array(vals),
                             np.array([r["classical_topm"][str(best_cl)] for r in rs]),
                             names=[r["pdb"] for r in rs])
                print(f"    a={a:<5g} {name:7s} vs classical m={best_cl}: "
                      f"{p['mean_diff']:+.3f} "
                      f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
                      f"W/L {p['n_better']}/{p['n_worse']}")


if __name__ == "__main__":
    main()
