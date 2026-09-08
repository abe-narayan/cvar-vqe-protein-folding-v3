"""SPRINT 15 / QGEOM -- PART D continued: CVaR AS AN ENSEMBLE GENERATOR, AND THE ENTROPY FLOOR.

Two things fell out of the D1/D2 trajectory that deserve their own experiment.

E1. THE ENTROPY FLOOR LAW.  Sprint 14 measured that the final distribution BROADENS as alpha
    falls (4.07 -> 10.26 bits) and attributed it to the gradient being noisier at small alpha
    (`alpha * shots` samples carry it).  With an EXACT gradient there is no shot noise at all,
    and the broadening is still there -- so the recorded mechanism cannot be the whole story.
    The mechanism proposed and tested here is Sprint 14's own iff condition seen as a STOPPING
    RULE:

        dCVaR_alpha/dp is identically zero  <=>  p(argmin E) >= alpha

    so an exact-gradient CVaR optimiser **stops concentrating the moment the argmin carries
    mass alpha**.  `alpha` is then not a learning rate and not a tail width but a
    **floor on the entropy of the answer**, with an exactly known stopping condition.
    PREDICTION, pre-registered: at convergence `p(argmin E) ~ alpha` from above, over orders
    of magnitude in alpha.

E2. THE UNRANKED SET READOUT.  D2 finds that the VQE's drawn sample has a significantly LOWER
    MEAN RMSD than the control's (-0.85 A at alpha=1, W/L 22/5) while its objective-SELECTED
    top-m is no better and often worse.  Those two facts are only compatible if concentration
    helps the set and selection then throws the advantage away -- which is the project's
    recorded situation, since nothing ranks within a pool.  So the honest question is:

        as an ENSEMBLE GENERATOR consumed WITHOUT a ranker, does a CVaR-VQE beat
        best-of-N from its own initial distribution?

    Readouts that use no objective at all: whole-set coordinate average, random-m coordinate
    average, set mean.  Diversity and contraction are reported alongside, because a narrow
    set's "coordinate average" is just its mode and must not be sold as an ensemble.

    python -m s15.qgeom_ens
"""
from __future__ import annotations

import time

import numpy as np

from core.quantum import cvar_exact
from s12 import instrument as I
from s15 import qgeom_lib as G
from s15.qgeom_cvar import Struct, make_circ, run_vqe

TAG = "ens"
TARGETS = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")
SUB = (1, 2, 3, 4, 5, 6)


# ==================================================================== E1
def entropy_floor(pdbs=TARGETS, alphas=(1.0, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005),
                  iters=400, seeds=(0, 1), residues=SUB, w=0.25):
    print("=" * 108)
    print("E1. THE ENTROPY FLOOR LAW.  Prediction: p(argmin E) -> alpha from above,")
    print("    because dCVaR/dp == 0 iff p(argmin) >= alpha.  EXACT gradient, no shots.")
    print("=" * 108)
    from s14.vqe_hamil import tabulate, combine
    prev = G.ck_load(TAG)
    out = prev.get("E1_entropy_floor", {})
    print(f"{'alpha':>7s} {'n':>4s} {'p(argmin)':>10s} {'p/alpha':>8s} {'max p':>8s} "
          f"{'H bits':>8s} {'ESSfrac':>9s} {'set-mean RMSD':>14s}")
    for pdb in pdbs:
        st = Struct(pdb, residues)
        pr, ds = tabulate(pdb)
        E = G.V.uniformise(combine(pr, ds, w)[st.full])
        xstar = int(np.argmin(E))
        circ = make_circ(st.n_qubits)
        for a in alphas:
            for s in seeds:
                key = f"{pdb}|a{a}|s{s}"
                if key in out:
                    continue
                r = run_vqe(circ, E, a, iters, seed=s, lowmem=True)
                p = r["p_final"]
                out[key] = {"p_argmin": float(p[xstar]), "alpha": a,
                            "max_p": float(p.max()),
                            "entropy_bits": r["final"]["entropy_bits"],
                            "ess_frac": r["final"]["ess_frac"],
                            "set_mean_rmsd": float(p @ st.rmsd)}
                G.ck(TAG, "E1_entropy_floor", out)
    for a in alphas:
        rows = [v for k, v in out.items() if f"|a{a}|" in k]
        if not rows:
            continue
        pa = float(np.mean([r["p_argmin"] for r in rows]))
        print(f"{a:7.3f} {len(rows):4d} {pa:10.5f} {pa / a:8.2f} "
              f"{np.mean([r['max_p'] for r in rows]):8.4f} "
              f"{np.mean([r['entropy_bits'] for r in rows]):8.3f} "
              f"{np.mean([r['ess_frac'] for r in rows]):9.5f} "
              f"{np.mean([r['set_mean_rmsd'] for r in rows]):14.3f}")
    print()
    print("  `p/alpha` near 1 confirms the stopping rule; the entropy column is then a")
    print("  DIAL and not an accident.  alpha=1 has no such condition and runs to a point.")
    return out


# ==================================================================== E2
def unranked_readouts(st, idx, ms=(5, 20, 75), rng=None):
    """Readouts that use NO objective-based selection at all."""
    idx = np.asarray(idx, np.int64)
    rng = rng or np.random.default_rng(0)
    R = st.rmsd                                          # ORACLE, post-hoc scoring only
    out = {"set_mean_rmsd": float(R[idx].mean()),
           "set_best_rmsd": float(R[idx].min()),
           "n_distinct": int(np.unique(idx).size)}
    u = np.unique(idx)
    out["set_coordavg_rmsd"] = st.coord_avg_rmsd(u[: min(len(u), 400)])
    for m in ms:
        pick = rng.choice(idx, size=m, replace=len(idx) < m)
        out[f"rand{m}_mean_rmsd"] = float(R[pick].mean())
        out[f"rand{m}_coordavg_rmsd"] = st.coord_avg_rmsd(np.unique(pick))
    # diversity: mean pairwise CA-RMSD inside a 40-member subsample
    sub = rng.choice(idx, size=min(40, len(idx)), replace=False)
    W = st.build(np.unique(sub))
    if len(W) > 1:
        d = [float(I.kabsch_rmsd_batch(W[i + 1:], W[i])[0]) for i in range(len(W) - 1)]
        out["mean_pairwise_rmsd"] = float(np.mean(d))
    else:
        out["mean_pairwise_rmsd"] = 0.0
    # contraction of the coordinate average against the native radius of gyration
    A = I.coordinate_average(W)
    A = A[0] if isinstance(A, tuple) else A
    rg = lambda X: float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
    out["rg_ratio_avg_over_native"] = rg(np.asarray(A, float)) / rg(st.nat)
    return out


def generator_test(pdbs=TARGETS, alphas=(1.0, 0.25, 0.05, 0.01), iters=200,
                   seeds=(0, 1, 2), budget=2048, residues=SUB, w=0.25):
    print()
    print("=" * 108)
    print(f"E2. CVaR-VQE AS AN ENSEMBLE GENERATOR, CONSUMED WITHOUT A RANKER.")
    print(f"    Matched budget ({budget} draws).  Control = best-of-N from the arm's own")
    print("    untrained initial distribution, the same draws, the same readouts.")
    print("=" * 108)
    from s14.vqe_hamil import tabulate, combine
    prev = G.ck_load(TAG)
    out = prev.get("E2_generator", {})
    for pdb in pdbs:
        st = Struct(pdb, residues)
        pr, ds = tabulate(pdb)
        E = G.V.uniformise(combine(pr, ds, w)[st.full])
        circ = make_circ(st.n_qubits)
        for a in alphas:
            for s in seeds:
                key = f"{pdb}|a{a}|s{s}"
                if key in out:
                    continue
                t0 = time.time()
                rng = np.random.default_rng(1000 * s + 7)
                r = run_vqe(circ, E, a, iters, seed=s, lowmem=True)
                ci = rng.choice(len(E), size=budget, p=r["p_init"])
                vi = rng.choice(len(E), size=budget, p=r["p_final"])
                out[key] = {
                    "control": unranked_readouts(st, ci, rng=np.random.default_rng(s)),
                    "vqe": unranked_readouts(st, vi, rng=np.random.default_rng(s)),
                    "entropy_bits": r["final"]["entropy_bits"],
                    "secs": time.time() - t0}
                G.ck(TAG, "E2_generator", out)
        print(f"  {pdb} done", flush=True)
    summarise_generator(out, alphas)
    return out


def summarise_generator(out, alphas=(1.0, 0.25, 0.05, 0.01)):
    print()
    print("=" * 108)
    print("E2b. PAIRED.  NEGATIVE = the VQE ensemble is better than best-of-N from its own")
    print("     untrained start, at the same budget, with NO ranker anywhere.")
    print("=" * 108)
    keys = ("set_mean_rmsd", "set_coordavg_rmsd", "rand5_coordavg_rmsd",
            "rand20_coordavg_rmsd", "rand75_coordavg_rmsd", "rand20_mean_rmsd",
            "set_best_rmsd", "mean_pairwise_rmsd", "n_distinct",
            "rg_ratio_avg_over_native")
    res = {}
    for ro in keys:
        print(f"\n  --- {ro} ---")
        print(f"  {'alpha':>6s} {'n':>4s} {'VQE':>9s} {'control':>9s} "
              f"{'diff':>10s} {'CI95':>24s} {'W/L':>7s} {'verdict':>12s}")
        for a in alphas:
            v = [x["vqe"][ro] for k, x in out.items() if f"|a{a}|" in k]
            c = [x["control"][ro] for k, x in out.items() if f"|a{a}|" in k]
            if not v:
                continue
            pr = G.paired(v, c)
            res[f"{ro}|a{a}"] = {"vqe_mean": float(np.mean(v)),
                                 "ctrl_mean": float(np.mean(c)), **pr}
            print(f"  {a:6.2f} {pr['n']:4d} {np.mean(v):9.3f} {np.mean(c):9.3f} "
                  f"{pr['mean']:+10.4f} "
                  f"{f'[{pr['ci_lo']:+.4f}, {pr['ci_hi']:+.4f}]':>24s} "
                  f"{f'{pr['win']}/{pr['loss']}':>7s} {pr['verdict']:>12s}")
    print()
    print("  READ `n_distinct` AND `mean_pairwise_rmsd` BEFORE READING ANY COORDINATE")
    print("  AVERAGE.  A set with one distinct member has a coordinate average equal to")
    print("  that member; calling it an ensemble would be a category error.")
    G.ck(TAG, "E2b_paired", res)
    return res


def main():
    G.wait_mem(0.8, "qgeom_ens")
    G.ck_load(TAG)
    entropy_floor()
    generator_test()
    print("\nwritten -> s15/results/qgeom_ens.json")


if __name__ == "__main__":
    main()
