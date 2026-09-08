"""SPRINT 14 / VQE -- PART A measurements: what alpha actually does.

The correctness of the estimator is settled in `s14/test_vqe.py`.  This module measures the
four things the brief asks for that a correctness test cannot answer:

  A. THE GRADIENT AUDIT, tabulated.  paramshift vs finite differences vs the sampled score
     estimator vs exact differentiation, across alpha x n x depth x seed, with the
     zero-gradient cells identified rather than averaged into a NaN.

  B. DOES ALPHA CHANGE THE TRAJECTORY OR MERELY RESCALE IT?  The cheap wrong way to ask is
     to compare final energies.  The right way is to evaluate the gradient of DIFFERENT
     alphas AT THE SAME POINT: if alpha only rescaled, every such pair would be parallel
     (cosine 1) and alpha would be a learning rate.  Then the same question dynamically, as
     parameter-space divergence between trajectories that share an initialisation and a
     seed and differ only in alpha.

  C. THE ALPHA SWEEP on a real enumerated target: stability, gradient variance, entropy,
     diversity, near-native probability mass, and both structural readouts.

  D. THE COLLAPSE TRAP.  The brief records that at alpha <= 0.25 two AMBER variants became
     literally the same objective.  Section A of `test_vqe.py` proves the mechanism as an
     identity; here it is MEASURED on the real pairs -- how much of the alpha-tail two
     genuinely different objectives share, and how parallel their CVaR gradients become, as
     a function of alpha.

    python -m s14.vqe_cvar
"""
from __future__ import annotations

import numpy as np

from core import quantum as Q
from s14 import vqe_lib as V
from s14 import vqe_run as R

ALPHAS = (1.0, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01)


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else float("nan")


# ------------------------------------------------------------------ A. the audit
def gradient_audit():
    print("=" * 78)
    print("A. GRADIENT AUDIT -- paramshift vs finite differences, exact statevector")
    print("=" * 78)
    print(f"{'n':>3s} {'L':>2s} {'alpha':>6s} {'cos':>14s} {'rel err':>10s} "
          f"{'zero-grad cells':>16s}")
    rows = []
    for n in (5, 7, 9, 11):
        for layers in (1, 2, 3, 5):
            for a in ALPHAS:
                cs, re, zero = [], [], 0
                for seed in range(4):
                    rng = np.random.default_rng(1000 * seed + 13 * n + layers)
                    circ = Q.StatevectorCircuit(n, layers)
                    th = rng.normal(0, 0.8, circ.n_params())
                    E = rng.normal(0, 1, 1 << n)
                    gp = Q.grad_cvar_paramshift(circ, th, E, a)
                    if np.linalg.norm(gp) == 0.0:
                        zero += 1
                        continue
                    gf = Q.grad_cvar_fd(circ, th, E, a, h=1e-5)
                    cs.append(_cos(gp, gf))
                    re.append(np.linalg.norm(gp - gf) / np.linalg.norm(gf))
                c = np.mean(cs) if cs else float("nan")
                r = np.mean(re) if re else float("nan")
                rows.append({"n": n, "layers": layers, "alpha": a, "cos": c,
                             "relerr": r, "zero_cells": zero, "of": 4})
                print(f"{n:3d} {layers:2d} {a:6.3g} {c:14.12f} {r:10.2e} {zero:10d}/4")
    ok = [x for x in rows if np.isfinite(x["cos"])]
    print(f"\nWORST cosine over {len(ok)} non-degenerate cells: "
          f"{min(x['cos'] for x in ok):.12f}")
    print(f"WORST relative error: {max(x['relerr'] for x in ok):.3e}")
    print(f"cells with an EXACTLY zero gradient: "
          f"{sum(x['zero_cells'] for x in rows)} of {4*len(rows)}")
    return rows


def sampled_estimator_audit():
    print()
    print("=" * 78)
    print("A2. SAMPLED SCORE ESTIMATOR vs the exact gradient, by baseline")
    print("    (4 x 200,000 shots, averaged: what is left is BIAS, not variance)")
    print("=" * 78)
    print(f"{'alpha':>7s} | {'const cos':>10s} {'|g|ratio':>9s} | {'tail cos':>9s} "
          f"{'|g|ratio':>9s} | {'none cos':>9s} {'|g|ratio':>9s}")
    out = []
    n, layers = 7, 3
    for a in ALPHAS:
        rng = np.random.default_rng(4)
        circ = Q.StatevectorCircuit(n, layers)
        th = rng.normal(0, 0.8, circ.n_params())
        E = rng.normal(0, 1, 1 << n)
        gp = Q.grad_cvar_paramshift(circ, th, E, a)
        row = {"alpha": a}
        for bl in ("const", "tail", "none"):
            g = np.mean([Q.grad_cvar_score(circ, th, E, a, shots=200000,
                                           rng=np.random.default_rng(s), baseline=bl)[0]
                         for s in range(4)], axis=0)
            row[bl + "_cos"] = _cos(g, gp)
            row[bl + "_ratio"] = float(np.linalg.norm(g) / np.linalg.norm(gp))
        out.append(row)
        print(f"{a:7.3g} | {row['const_cos']:10.6f} {row['const_ratio']:9.4f} | "
              f"{row['tail_cos']:9.6f} {row['tail_ratio']:9.4f} | "
              f"{row['none_cos']:9.6f} {row['none_ratio']:9.4f}")
    print("\nThe defect's severity is alpha- and instance-dependent: tail-baseline cosine")
    print(f"ranges {min(r['tail_cos'] for r in out):+.3f} to "
          f"{max(r['tail_cos'] for r in out):+.3f} across alpha at a FIXED point.")
    print("The recorded -0.023 is one draw from that spread, not a constant.")
    return out


# -------------------------------------------------- B. rescale or new direction?
def alpha_is_not_a_learning_rate():
    print()
    print("=" * 78)
    print("B1. IS ALPHA A RESCALING?  Gradients of different alpha AT THE SAME POINT.")
    print("    If alpha only rescaled the step, every cosine would be exactly 1.")
    print("=" * 78)
    n, layers = 9, 3
    rng = np.random.default_rng(21)
    M = np.zeros((len(ALPHAS), len(ALPHAS)))
    NR = np.zeros(len(ALPHAS))
    G = []
    for a in ALPHAS:
        gs = []
        for seed in range(6):
            r = np.random.default_rng(500 + seed)
            circ = Q.StatevectorCircuit(n, layers)
            th = r.normal(0, 0.8, circ.n_params())
            E = r.normal(0, 1, 1 << n)
            gs.append(Q.grad_cvar_paramshift(circ, th, E, a))
        G.append(gs)
        NR[ALPHAS.index(a)] = np.mean([np.linalg.norm(g) for g in gs])
    for i in range(len(ALPHAS)):
        for j in range(len(ALPHAS)):
            M[i, j] = np.mean([_cos(G[i][s], G[j][s]) for s in range(6)])
    print("cosine between grad(alpha_i) and grad(alpha_j) at identical theta, 6 seeds")
    print("        " + "".join(f"{a:>8.3g}" for a in ALPHAS))
    for i, a in enumerate(ALPHAS):
        print(f"{a:7.3g} " + "".join(f"{M[i,j]:8.4f}" for j in range(len(ALPHAS))))
    print("\nmean |grad| by alpha: " + "  ".join(f"{a:g}:{v:.4f}"
                                                 for a, v in zip(ALPHAS, NR)))
    print(f"off-diagonal cosine range: {M[~np.eye(len(ALPHAS), dtype=bool)].min():.4f}"
          f" .. {M[~np.eye(len(ALPHAS), dtype=bool)].max():.4f}")
    return {"cos_matrix": M.tolist(), "alphas": list(ALPHAS), "grad_norms": NR.tolist()}


def trajectory_divergence(pdb="1CS9", signal=0.3, budget=20480, shots=512):
    print()
    print("=" * 78)
    print(f"B2. TRAJECTORY DIVERGENCE in parameter space ({pdb}, signal={signal})")
    print("    Same seed, same initialisation, same shot stream: ONLY alpha differs.")
    print("=" * 78)
    e = V.Enum(pdb)
    E = V.blend_objective(e.legacy, e.rmsd, signal)
    ref = None
    print(f"{'alpha':>7s} {'|theta-theta0|':>15s} {'path len':>10s} "
          f"{'|th - th(a=1)|':>15s} {'cos(dtheta, dtheta_a1)':>23s}")
    out = []
    for a in ALPHAS:
        r = R.run(E, e.n_qubits, a, budget=budget, shots=shots, ansatz="mps2f",
                  seed=0, rmsd=e.rmsd)
        th = np.array(r["theta"]); th0 = np.array(r["theta0"])
        if ref is None:
            ref = th
        d = float(np.linalg.norm(th - ref))
        c = _cos(th - th0, ref - th0)
        out.append({"alpha": a, "disp": r["param_displacement"],
                    "path": r["param_path_length"], "dist_to_a1": d, "cos_to_a1": c,
                    "rmsd_returned": r["rmsd_returned"], "mode_rmsd": r["mode_rmsd"],
                    "entropy_bits": r["entropy_bits"]})
        print(f"{a:7.3g} {r['param_displacement']:15.4f} "
              f"{r['param_path_length']:10.3f} {d:15.4f} {c:23.4f}")
    print("\nIf alpha were a learning rate the trajectories would be COLLINEAR")
    print("(cos = 1) and differ only in how far along they got.")
    return out


# ----------------------------------------------------------------- C. the sweep
def alpha_sweep(pdb="1CS9", signal=0.3, budget=20480, shots=512, seeds=(0, 1, 2, 3, 4)):
    print()
    print("=" * 78)
    print(f"C. ALPHA SWEEP on {pdb}, signal={signal} (rho measured below), "
          f"budget={budget} evals")
    print("=" * 78)
    e = V.Enum(pdb)
    E = V.blend_objective(e.legacy, e.rmsd, signal)
    exact_i = int(np.argmin(E))
    print(f"rho(objective, true RMSD) = {V.spearman(E, e.rmsd):+.4f};  "
          f"exact optimum RMSD = {e.rmsd[exact_i]:.3f};  space best = {e.rmsd.min():.3f};"
          f"  random draw = {e.rmsd.mean():.3f}")
    print()
    hdr = (f"{'alpha':>8s} {'objgap':>8s} {'sd':>7s} {'RMSD ret':>9s} {'sd':>6s} "
           f"{'RMSD mode':>10s} {'H bits':>7s} {'div':>6s} {'gradvar':>9s} "
           f"{'P(<2.5A)':>9s} {'maxp':>8s}")
    print(hdr)
    rows = []
    settings = [(a, False) for a in ALPHAS] + [((0.5, 0.05), True), ((1.0, 0.01), True)]
    for a, anneal in settings:
        rs = [R.run(E, e.n_qubits, a, budget=budget, shots=shots, ansatz="mps2f",
                    seed=s, rmsd=e.rmsd, alpha_anneal=anneal) for s in seeds]
        og = [r["best_e"] - E[exact_i] for r in rs]
        rr = [r["rmsd_returned"] for r in rs]
        row = {"alpha": a if not anneal else list(a), "anneal": anneal,
               "objective_gap": float(np.mean(og)), "objective_gap_sd": float(np.std(og)),
               "rmsd_returned": float(np.mean(rr)), "rmsd_returned_sd": float(np.std(rr)),
               "mode_rmsd": float(np.mean([r["mode_rmsd"] for r in rs])),
               "entropy_bits": float(np.mean([r["entropy_bits"] for r in rs])),
               "diversity": float(np.mean([r["diversity"] for r in rs])),
               "grad_var": float(np.mean([r["grad_var"] for r in rs])),
               "pmass_below_2.5": float(np.mean([r["pmass_below_2.5"] for r in rs])),
               "max_prob": float(np.mean([r["max_prob"] for r in rs])),
               "init_rmsd_mean": float(np.mean([r["init_rmsd_mean"] for r in rs]))}
        rows.append(row)
        lab = f"{a}" if not anneal else f"{a[0]}->{a[1]}"
        print(f"{lab:>8s} {row['objective_gap']:8.4f} {row['objective_gap_sd']:7.4f} "
              f"{row['rmsd_returned']:9.3f} {row['rmsd_returned_sd']:6.3f} "
              f"{row['mode_rmsd']:10.3f} {row['entropy_bits']:7.2f} "
              f"{row['diversity']:6.3f} {row['grad_var']:9.2e} "
              f"{row['pmass_below_2.5']:9.4f} {row['max_prob']:8.4f}")
    print(f"\ninitialisation RMSD (mean over sampled starts): "
          f"{rows[0]['init_rmsd_mean']:.3f} -- the VQE starts at a random point, not at "
          f"the answer.")
    return rows


# -------------------------------------------------------------- D. the collapse
def collapse_trap(pdb="1CS9"):
    print()
    print("=" * 78)
    print("D. THE COLLAPSE TRAP, measured on real objective pairs")
    print("   How much of the alpha-tail do two genuinely different objectives share,")
    print("   and how parallel do their CVaR gradients become?")
    print("=" * 78)
    e = V.Enum(pdb)
    # AMBER variants live on the 2,955-configuration labelled subset
    d = np.load(f"{V.ENUM_DIR}/qarch_enum_{pdb}.npz")
    amb = {c[4:]: np.asarray(d[c], float) for c in d.files if c.startswith("amb_")}
    tot = np.asarray(d["amber_total"], float)
    pairs = {
        "AMBER total vs AMBER minus solvation":
            (V.uniformise(tot), V.uniformise(tot - amb["solvation"])),
        "AMBER total vs AMBER nonbonded+solvation":
            (V.uniformise(tot), V.uniformise(amb["nonbonded"] + amb["solvation"])),
        "Legacy total vs Legacy minus steric (subset)":
            (V.uniformise(e.legacy[np.asarray(d["amber_idx"])]),
             V.uniformise((e.legacy - e.leg["steric"])[np.asarray(d["amber_idx"])])),
    }
    M = len(tot)
    rng = np.random.default_rng(0)
    # a realistic, mildly concentrated distribution over the subset
    p = rng.dirichlet(np.full(M, 0.6))
    print(f"\nsubset size {M};  distribution: Dirichlet(0.6), max p = {p.max():.4f}")
    print(f"\n{'pair':46s}" + "".join(f"{a:>8.3g}" for a in ALPHAS))
    out = {}
    for name, (E1, E2) in pairs.items():
        rho = V.spearman(E1, E2)
        jac, gcos = [], []
        for a in ALPHAS:
            _, _, m1 = Q.cvar_from_probs(E1, p, a)
            _, _, m2 = Q.cvar_from_probs(E2, p, a)
            s1, s2 = m1 > 0, m2 > 0
            jac.append(float((s1 & s2).sum() / max(1, (s1 | s2).sum())))
            d1 = Q.cvar_exact(E1, p, a)[2]
            d2 = Q.cvar_exact(E2, p, a)[2]
            gcos.append(_cos(d1, d2))
        out[name] = {"rho": rho, "tail_jaccard": jac, "dcvar_cosine": gcos}
        print(f"{name[:44]:46s}" + "".join(f"{v:8.3f}" for v in jac)
              + f"   (rho between the two objectives = {rho:+.3f})")
        print(f"{'  cosine of dCVaR/dp':46s}" + "".join(f"{v:8.3f}" for v in gcos))
    print("\nA Jaccard of 1.000 and a dCVaR cosine of 1.000 mean the two objectives are")
    print("LITERALLY the same CVaR function at that alpha, however different they are")
    print("elsewhere.  That is the recorded trap, and it is a property of CVaR.")
    return out


def main():
    V.wait_for_memory(1.0, "vqe_cvar")
    out = {}
    out["gradient_audit"] = gradient_audit()
    out["sampled_audit"] = sampled_estimator_audit()
    out["alpha_rescale"] = alpha_is_not_a_learning_rate()
    out["trajectory"] = trajectory_divergence()
    out["sweep_sig0.3"] = alpha_sweep(signal=0.3)
    out["sweep_sig0.0"] = alpha_sweep(signal=0.0)
    out["sweep_sig1.0"] = alpha_sweep(signal=1.0)
    out["collapse"] = collapse_trap()
    V.write("vqe_cvar", out)
    print("\nwritten -> s14/results/vqe_cvar.json")


if __name__ == "__main__":
    main()
