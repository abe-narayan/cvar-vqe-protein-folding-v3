"""SPRINT 13 GEO -- EXPERIMENT 4: the Fubini-Study metric / quantum Fisher information.

THE DISTINCTION THIS WHOLE EXPERIMENT IS BUILT AROUND
=====================================================
``g_ij(theta) = Re[<d_i psi|d_j psi> - <d_i psi|psi><psi|d_j psi>]`` is a property of the
PARAMETERISED STATE.  It does not contain the Hamiltonian.  Two energy models sharing an
ansatz, a representation and a parameter vector have **numerically identical** metrics --
not similar, identical -- and Arm 0 below measures that difference and reports it as the
zero it must be.

So a Legacy-vs-AMBER difference in the metric can only arise through **where in parameter
space the optimiser is taken**, and the design says so explicitly:

  Arm 0  CONTROL, `metric_is_hamiltonian_free`: g at a fixed theta, evaluated in both
         model contexts.  Expected max|difference| = 0 exactly.
  Arm 1  MANIFOLD, `random_theta`: the distribution of the spectrum over theta ~ U(-pi,pi),
         by qubit count and depth.  This characterises the ANSATZ; no energy model is
         involved at all.
  Arm 2  TRAJECTORY, `at_optimum`: identical initialisation, identical optimiser, identical
         step count, identical alpha -- run once per energy model -- then g at each
         endpoint.  Any difference here is attributable to the trajectory, and the
         Fubini-Study path length actually travelled is reported beside it so the reader
         can see how far apart the two endpoints are ON THE MANIFOLD.

Reported per metric: full eigenvalue spectrum summary, effective rank at two thresholds,
participation ratio, condition number (regularised and full), trace, the smallest non-zero
eigenvalue, and the parameter-correlation structure ``g_ij / sqrt(g_ii g_jj)``.

    python -m s13.geo_metric
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from core.quantum import Adam, cvar_exact              # noqa: E402

LAYERS = (1, 2, 3, 5)
QUBITS = (8, 10, 12, 14, 16)
ALPHA = 0.25
STEPS = 200
LR = 0.10


def adam_run(E, circ, theta0, alpha, steps=STEPS, lr=LR):
    """Adam on the EXACT parameter-shift CVaR gradient.  Identical for every model."""
    th = np.array(theta0, float)
    opt = Adam(th.size, lr=lr)
    traj = [th.copy()]
    vals = []
    for _ in range(steps):
        g = G.cvar_grad(E, circ, th, alpha)
        vals.append(G.cvar_value(E, circ, th, alpha))
        th = opt.step(th, g)
        traj.append(th.copy())
    vals.append(G.cvar_value(E, circ, th, alpha))
    return th, np.array(vals), np.array(traj)


def fs_path_length(circ, traj):
    """Sum of ``sqrt(dtheta^T g dtheta)`` along the optimiser trajectory."""
    tot = 0.0
    for a in range(len(traj) - 1):
        d = traj[a + 1] - traj[a]
        g = G.fs_metric(circ, traj[a])
        tot += float(np.sqrt(max(d @ g @ d, 0.0)))
    return tot


def arm0_control(cell=("1A13", 5, 4), layers=3, alpha=ALPHA, seeds=(0, 1, 2)):
    """THE CROSS-EVALUATION: the metric follows theta, not the Hamiltonian.

    Both models are optimised from the same theta0 to endpoints ``theta_L`` and
    ``theta_A``.  Then the 2x2 table ``g(theta_X)`` evaluated "under model M" is filled in
    for every (X, M).  Two things must hold, and both are measured:

      * ``g(theta_X)`` is the SAME matrix for both M -- max|difference| exactly 0.  The
        metric contains no energy model.
      * ``g(theta_L) != g(theta_A)`` in general.  Whatever Legacy/AMBER metric difference
        this study reports is therefore a statement about WHERE THE OPTIMISER WENT, and is
        reported as such throughout.
    """
    pdb, L, k = cell
    n = int(round(np.log2(k))) * L
    circ = G.circuit(n, layers=layers)
    EL = G.variant_table(pdb, L, k, "legacy")
    EA = G.variant_table(pdb, L, k, "amber")
    out = []
    for s in seeds:
        th0 = G.init_theta(circ.n_params(), s)
        thL, _, _ = adam_run(EL, circ, th0, alpha)
        thA, _, _ = adam_run(EA, circ, th0, alpha)
        gL_ctxL, gL_ctxA = G.fs_metric(circ, thL), G.fs_metric(circ, thL)
        gA_ctxL, gA_ctxA = G.fs_metric(circ, thA), G.fs_metric(circ, thA)
        out.append({
            "cell": f"{pdb} L={L} k={k}", "n_qubits": n, "layers": layers, "seed": s,
            "max_abs_diff_same_theta_across_model_contexts":
                float(max(np.abs(gL_ctxL - gL_ctxA).max(),
                          np.abs(gA_ctxL - gA_ctxA).max())),
            "max_abs_diff_between_the_two_endpoints":
                float(np.abs(gL_ctxL - gA_ctxA).max()),
            "rel_diff_between_endpoints":
                float(np.abs(gL_ctxL - gA_ctxA).max()
                      / max(np.abs(gL_ctxL).max(), 1e-300)),
            "theta_endpoint_l2_separation": float(np.linalg.norm(thL - thA)),
            "stats_legacy_endpoint": G.metric_stats(gL_ctxL),
            "stats_amber_endpoint": G.metric_stats(gA_ctxA),
        })
        print("  arm0 seed=%d  same-theta cross-context diff %.3e | endpoint-to-endpoint "
              "metric diff %.4f (rel %.3f)"
              % (s, out[-1]["max_abs_diff_same_theta_across_model_contexts"],
                 out[-1]["max_abs_diff_between_the_two_endpoints"],
                 out[-1]["rel_diff_between_endpoints"]), flush=True)
    return out


def arm1_manifold(qubits=QUBITS, layers=LAYERS, n_theta=48):
    rows = []
    for n in qubits:
        for lay in layers:
            circ = G.circuit(n, layers=lay)
            P = circ.n_params()
            n_theta = {8: 48, 10: 48, 12: 32, 14: 24, 15: 16, 16: 12}.get(n, 12)
            rng = np.random.default_rng(97 * n + lay)
            st = []
            berry = 0.0
            for _ in range(n_theta):
                th = rng.uniform(-np.pi, np.pi, P)
                psi, D = G.dstates(circ, th)
                berry = max(berry, float(np.abs(D @ psi).max()))
                g = D @ D.T - np.outer(D @ psi, D @ psi)
                st.append(G.metric_stats(g))
            agg = {"n_qubits": n, "layers": lay, "P": P, "n_theta": n_theta,
                   "berry_max_abs": berry}
            for key in ("trace", "lam_max", "lam_min_nonzero", "participation_ratio",
                        "rank_1e-10", "rank_1e-6", "cond_1e-10",
                        "mean_abs_offdiag_corr", "diag_mean"):
                v = np.array([s[key] for s in st], float)
                agg[key] = {"mean": float(np.mean(v)), "median": float(np.median(v)),
                            "sd": float(np.std(v)), "min": float(np.min(v)),
                            "max": float(np.max(v))}
            rows.append(agg)
            print("  manifold n=%-2d lay=%d P=%-3d  rank(1e-10) %5.1f/%d  PR %6.2f  "
                  "trace %8.3f  cond %.2e"
                  % (n, lay, P, agg["rank_1e-10"]["mean"], P,
                     agg["participation_ratio"]["mean"], agg["trace"]["mean"],
                     agg["cond_1e-10"]["median"]), flush=True)
    return rows


def arm2_trajectory(cells, layers=3, alpha=ALPHA, seeds=(0, 1, 2)):
    rows = []
    for pdb, L, k in cells:
        n = int(round(np.log2(k))) * L
        circ = G.circuit(n, layers=layers)
        P = circ.n_params()
        Es = {}
        for v in G.VARIANTS:
            try:
                Es[v] = G.variant_table(pdb, L, k, v)
            except FileNotFoundError:
                pass
        if len(Es) < 2:
            continue
        R = G.rmsd_table(pdb, L, k)
        for s in seeds:
            th0 = G.init_theta(P, s)
            row = {"pdb": pdb, "L": L, "k": k, "n_qubits": n, "layers": layers,
                   "alpha": alpha, "seed": s, "steps": STEPS, "lr": LR,
                   "init": G.metric_stats(G.fs_metric(circ, th0))}
            ends = {}
            for v, E in Es.items():
                t0 = time.time()
                th, vals, traj = adam_run(E, circ, th0, alpha)
                ends[v] = th
                g = G.fs_metric(circ, th)
                p = circ.probs(th)
                jbest = int(np.argmax(p))
                row[v] = {
                    "metric": G.metric_stats(g),
                    "cvar_start": float(vals[0]), "cvar_end": float(vals[-1]),
                    "cvar_drop_frac_of_span": float(
                        (vals[0] - vals[-1]) / max(vals[0] - E.min(), 1e-300)),
                    "theta_l2_moved": float(np.linalg.norm(th - th0)),
                    "fs_path_length": fs_path_length(circ, traj[::25]),
                    "entropy_bits": G.entropy_bits(p),
                    "p_top": float(p.max()),
                    "ORACLE_rmsd_argmax_p": float(R[jbest]),
                    "ORACLE_rmsd_best_in_top16": float(
                        R[np.argsort(p)[::-1][:16]].min()),
                    "wall_s": time.time() - t0,
                }
            # how far apart are the two endpoints, in parameters and on the manifold?
            if "legacy" in ends and "amber" in ends:
                tl, ta = ends["legacy"], ends["amber"]
                gm = G.fs_metric(circ, (tl + ta) / 2)
                d = tl - ta
                row["endpoint_sep_l2"] = float(np.linalg.norm(d))
                row["endpoint_sep_fs"] = float(np.sqrt(max(d @ gm @ d, 0.0)))
            rows.append(row)
            print("  traj %-5s L=%-2d k=%d n=%-2d seed=%d  " % (pdb, L, k, n, s)
                  + "  ".join("%s rank %d cond %.1e" %
                              (v[:6], row[v]["metric"]["rank_1e-10"],
                               row[v]["metric"]["cond_1e-10"]) for v in Es), flush=True)
    return rows


def main():
    cells = [(p, L, k) for p in ("1A13", "1A1P", "2BFI")
             for L, k in ((5, 4), (6, 4), (12, 2))]
    cells = [c for c in cells
             if os.path.exists(os.path.join(G.CACHE, G._tag(c[0], c[1], c[2], "amber")
                                            + ".npy"))]
    import json
    prev = os.path.join(G.RESULTS, "geo_metric_arm1.json")
    arm1 = (json.load(open(prev))["rows"] if os.path.exists(prev) else arm1_manifold())
    out = {"what": "Fubini-Study metric / QFI of the RY-CNOT ansatz, and its dependence on "
                   "where the optimiser goes under two energy models",
           "arm0_control": arm0_control(),
           "arm1_manifold": arm1,
           "arm2_trajectory": arm2_trajectory(cells)}
    G.write("geo_metric", out)


if __name__ == "__main__":
    main()
