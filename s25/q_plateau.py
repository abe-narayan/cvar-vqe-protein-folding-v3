"""S25 / LANE Q -- GRADIENT-VARIANCE SCALING OF THE CVaR OBJECTIVE ON THE DEPLOYED ANSATZ.

PROPERTY MEASUREMENT. Reads no RMSD, touches no native, has no outcome that can favour a
hypothesis about the endpoint, so Rule 0's fork enumeration does not apply (s24 Lane-D
precedent). It exists because `s25/QUANTUM.md` has a barren-plateau section and the honest
alternative to measuring this was writing NOT MEASURED in it.

THE QUESTION, AND WHY THIS DESIGN ANSWERS IT
============================================
A barren plateau is the statement that Var_theta[d C / d theta_k] decays exponentially in the
qubit count. The standard results are proved for a cost that is LINEAR in the state,
C = <psi|H|psi>. **CVaR_alpha is not linear in p** -- it is a concave, piecewise-linear
functional of the probability vector -- so the standard picture does not transfer by
inspection, and this project has never measured whether it transfers by fact.

The design gives that its own control, at zero extra cost:

    alpha = 1.0   CVaR reduces EXACTLY to the mean of E under p, i.e. <psi|H|psi> with
                  H = diag(E). This is the LINEAR cost -- the textbook barren-plateau setting.
    alpha < 1.0   the same circuit, the same H, the same theta samples, the same estimator,
                  differing ONLY in the non-linearity.

So the alpha=1 row is not an analogy; it is the linear special case of the same expression,
and the alpha<1 rows are the identical experiment with the non-linearity switched on.

T = 0 isolates the CVaR term; T = 0.3 is the DEPLOYED free energy and adds -T H(p), which is
also non-linear in p. Both are reported so the two non-linearities are not conflated.

Gradients are the EXACT parameter-shift gradient (`grad_cvar_paramshift` / `free_energy`), so
every number here is free of shot noise by construction and the variance measured is variance
over THETA alone -- which is what a barren-plateau statement is about.

SCOPE, STATED SO IT IS NOT OVER-READ
====================================
* The deployed register is n = 7. Everything above n = 7 here is an EXTRAPOLATION INSTRUMENT,
  not a deployed configuration, and is labelled so.
* P = layers * n parameters. At layers = 3, P = 3n, which is exponentially smaller than
  dim so(2**n) = 2**(n-1)(2**n - 1). This circuit is nowhere near a 2-design at any n
  measured here, so an observed decay rate is a property of THIS shallow ansatz and must not
  be quoted as a 2-design result.
* Real-amplitude ansatz: RY and CNOT are real, so the reachable manifold lies in SO(2**n),
  not SU(2**n). Stated because it changes which Lie algebra any theory statement would be
  about.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q            # noqa: E402
from s24 import stats_lib as ST          # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
R = {}


def deployed_E(dim: int) -> np.ndarray:
    """The DEPLOYED energy vector's shape: zrank of the score over the filtered set.

    `core/pipeline.quantum_stage` sets `E = _zrank(score[top[:dim]])`, and `top` is already
    the score ORDER, so E is the standardised ranks 1..dim -- the same vector for every
    target. Using it here (rather than iid noise) keeps the spectrum of H exactly the one the
    production objective sees, at every register size.
    """
    r = np.arange(1, dim + 1, dtype=float)
    return (r - r.mean()) / r.std()


def measure(n: int, layers: int, alpha: float, T: float, n_theta: int, seed: int,
            init_sd: float) -> dict:
    circ = Q.StatevectorCircuit(n, layers)
    E = deployed_E(circ.dim)
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    g0, gn = [], []
    for _ in range(n_theta):
        th = rng.normal(0.0, init_sd, P)
        if T == 0.0:
            g = Q.grad_cvar_paramshift(circ, th, E, alpha)
        else:
            _, g, _, _, _ = Q.free_energy(circ, th, E, alpha, T)
        g0.append(float(g[0]))
        gn.append(float(np.dot(g, g)))
    g0 = np.asarray(g0)
    gn = np.asarray(gn)
    return dict(n=n, layers=layers, P=P, dim=circ.dim, alpha=alpha, T=T,
                var_g0=float(g0.var(ddof=1)), mean_g0=float(g0.mean()),
                mean_sq_norm=float(gn.mean()),
                mean_sq_per_param=float(gn.mean() / P), n_theta=n_theta)


#: draws per register size. Variance of a variance goes as 2/(m-1), so m=250 gives the
#: estimate a ~9% relative SE; the big registers are cut to 80 (~16%) because a single
#: parameter-shift gradient there is 2P exact simulations of 2**n amplitudes. The decay
#: SLOPE is read off five decades, so a 16% error bar on the last point cannot move it.
DRAWS = {4: 250, 6: 250, 7: 250, 8: 250, 10: 200, 12: 120, 13: 80}


def sweep(tag, ns, layers, alpha, T, n_theta=None, init_sd=0.6):
    print(f"\n  {tag}   layers={layers} alpha={alpha} T={T} "
          f"theta ~ N(0, {init_sd}^2)")
    print(f"    {'n':>3} {'dim':>7} {'P':>4} {'draws':>6} {'Var[dF/dth_0]':>16} "
          f"{'E|g|^2/P':>14}")
    rows = []
    for n in ns:
        nt = n_theta or DRAWS.get(n, 100)
        r = measure(n, layers, alpha, T, nt, seed=1000 + n, init_sd=init_sd)
        rows.append(r)
        print(f"    {n:3d} {r['dim']:7d} {r['P']:4d} {nt:6d} {r['var_g0']:16.6e} "
              f"{r['mean_sq_per_param']:14.6e}", flush=True)
    ns_ = np.array([r["n"] for r in rows], float)
    v = np.array([r["var_g0"] for r in rows], float)
    ok = v > 0
    slope = float(np.polyfit(ns_[ok], np.log2(v[ok]), 1)[0]) if ok.sum() > 2 else float("nan")
    print(f"    fitted decay: log2 Var ~ {slope:+.4f} * n      "
          f"(a barren plateau from a 2-design would be about -1 per qubit; "
          f"{'consistent' if slope < -0.8 else 'NOT consistent'})")
    return dict(rows=rows, log2_slope_per_qubit=slope)


if __name__ == "__main__":
    NS = [4, 6, 7, 8, 10, 12, 13]
    print("GRADIENT VARIANCE vs REGISTER SIZE -- exact parameter shift, no shot noise")
    print("The deployed register is n=7, layers=3. Larger n is an EXTRAPOLATION INSTRUMENT.")

    # THE CONTROL AND THE TREATMENT, same circuit, same H, same theta law.
    R["linear_alpha1_T0"] = sweep("LINEAR COST (alpha=1, T=0): CVaR == <psi|H|psi>",
                                  NS, 3, 1.0, 0.0)
    R["cvar_alpha025_T0"] = sweep("NON-LINEAR CVaR (alpha=0.25, T=0)", NS, 3, 0.25, 0.0)
    R["cvar_alpha01_T0"] = sweep("NON-LINEAR CVaR (alpha=0.10, T=0)", NS, 3, 0.10, 0.0)
    R["deployed_a1_T03"] = sweep("DEPLOYED free energy (alpha=1.0, T=0.3)", NS, 3, 1.0, 0.3)
    R["deployed_a025_T03"] = sweep("DEPLOYED free energy (alpha=0.25, T=0.3)",
                                   NS, 3, 0.25, 0.3)

    # depth, at the deployed width
    print("\n  DEPTH SWEEP AT THE DEPLOYED WIDTH n=7 (alpha=0.25, T=0.3)")
    print(f"    {'L':>3} {'P':>4} {'Var[dF/dth_0]':>16} {'E|g|^2/P':>14}")
    dep = []
    for L in (1, 2, 3, 4, 6, 8, 12):
        r = measure(7, L, 0.25, 0.3, 250, seed=77 + L, init_sd=0.6)
        dep.append(r)
        print(f"    {L:3d} {r['P']:4d} {r['var_g0']:16.6e} {r['mean_sq_per_param']:14.6e}")
    R["depth_sweep_n7"] = dep

    # the ratio that answers the coordinator's question directly, at every n
    print("\n  DOES THE NON-LINEARITY CHANGE THE VARIANCE PICTURE?")
    print("    ratio Var[grad CVaR_alpha] / Var[grad MEAN], same circuit, same theta law")
    print(f"    {'n':>3} {'alpha=0.25':>14} {'alpha=0.10':>14}")
    lin = {r["n"]: r["var_g0"] for r in R["linear_alpha1_T0"]["rows"]}
    a25 = {r["n"]: r["var_g0"] for r in R["cvar_alpha025_T0"]["rows"]}
    a10 = {r["n"]: r["var_g0"] for r in R["cvar_alpha01_T0"]["rows"]}
    ratios = {}
    for n in NS:
        ratios[n] = (a25[n] / lin[n], a10[n] / lin[n])
        print(f"    {n:3d} {ratios[n][0]:14.4f} {ratios[n][1]:14.4f}")
    R["nonlinearity_variance_ratio"] = {str(k): list(v) for k, v in ratios.items()}

    # the entanglement control the project already has a standing result about
    print("\n  ENTANGLEMENT CONTROL at n=7: chi=2**layers vs a product state (MPS path)")
    for ent in ("cnot", "none"):
        a = Q.MPSAnsatz(7, layers=2, final_ry=True, entangler=ent)
        print(f"    entangler={ent:<5} chi={a.chi}  n_params={a.n_params()}")
    print("    (the DEPLOYED selector is StatevectorCircuit, not MPSAnsatz; the standing")
    print("     project result on chi is s21/s20 -- entanglement worth -0.013 [-0.095,+0.077],")
    print("     NOT MEASURED, not refuted. Not re-derived here.)")

    ST.save_atomic(os.path.join(OUT, "q_plateau.json"),
                   dict(kind="property measurement, gradient variance, no RMSD", lane="Q",
                        sprint=25, results=R), module_file=__file__)
    print("\nwrote", os.path.join(OUT, "q_plateau.json"))
