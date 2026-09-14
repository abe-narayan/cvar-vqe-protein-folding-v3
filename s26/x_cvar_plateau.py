#!/usr/bin/env python
"""Graph 2: gradient variance of the project's CVaR objective versus qubit count, one line per
alpha, on the deployed circuit and objective.

Everything is the project's own code, unchanged:
  circuit    `core.quantum.StatevectorCircuit(n, layers=3)` (RY layers + CNOT chain/ring, real
             amplitudes, exact dense statevector; the deployed selector is n = 7, layers = 3)
  objective  CVaR_alpha(E; p) with p = |psi(theta)|^2 and E the DEPLOYED energy shape
             (standardised ranks 1..2^n, exactly what `core.pipeline.quantum_stage` hands the
             selector on every target, `s25.q_plateau.deployed_E`)
  gradient   the EXACT parameter-shift gradient `core.quantum.grad_cvar_paramshift`
             (no shot noise: the variance is over theta alone)
  theta law  theta ~ N(0, 0.6^2) per parameter, as `core.quantum.run_cvar_vqe` initialises
  sampler    `s25.q_plateau.measure` (the S25 width-sweep instrument, reproduced bit-exactly
             by lane Q in S26 L35), seed 1000 + n per register so the SAME theta draws are used
             for every alpha at a given n

Quantity plotted: Var_theta[dCVaR/dtheta_0] (the first parameter, the S25 convention) with a
bootstrap 95% CI over the theta draws; the mean squared gradient per parameter is stored beside
it in the data file. T = 0 isolates the CVaR term (the deployed free energy adds -T*H(p) at
T = 0.3; those rows are also computed and stored, not plotted). No fit is drawn; the least-
squares slope of log2 Var against n is reported in the data file with its own CI so the reader
can judge whether an exponential law is supported, per the project's rule that a small
gradient is not a barren plateau without the variance-versus-width analysis.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q            # noqa: E402
from s25 import q_plateau as QP          # noqa: E402

ALPHAS = [1.0, 0.5, 0.25, 0.10]
NS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
LAYERS = 3
DRAWS = {4: 250, 5: 250, 6: 250, 7: 250, 8: 250, 9: 200, 10: 200, 11: 150, 12: 120, 13: 80}
INIT_SD = 0.6
OUT_JSON = os.path.join(HERE, "results", "x_cvar_plateau.json")
OUT_PNG = os.path.join(HERE, "figures", "graph2_gradient_variance_vs_qubits.png")


def samples(n, alpha, T, n_theta, seed):
    """The per-draw gradient samples behind `s25.q_plateau.measure`, kept so a bootstrap CI
    can be formed; same circuit, energy, theta law, seed and gradient call as `measure`."""
    circ = Q.StatevectorCircuit(n, LAYERS)
    E = QP.deployed_E(circ.dim)
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    g0, gn = np.empty(n_theta), np.empty(n_theta)
    for k in range(n_theta):
        th = rng.normal(0.0, INIT_SD, P)
        if T == 0.0:
            g = Q.grad_cvar_paramshift(circ, th, E, alpha)
        else:
            _, g, _, _, _ = Q.free_energy(circ, th, E, alpha, T)
        g0[k] = float(g[0])
        gn[k] = float(np.dot(g, g))
    return g0, gn, P, circ.dim


def boot_var(x, rng, nb=2000):
    v = np.array([x[rng.integers(0, len(x), len(x))].var(ddof=1) for _ in range(nb)])
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def main():
    rows = []
    rng_b = np.random.default_rng(26)
    ck = OUT_JSON + ".rows.tmp"
    if os.path.exists(ck):
        rows = json.load(open(ck, encoding="utf-8"))["rows"]
        print(f"  resuming with {len(rows)} cells")
    done = {(r["n"], r["alpha"], r["T"]) for r in rows}
    for T in (0.0, 0.3):
        for alpha in ALPHAS:
            for n in NS:
                if (n, alpha, T) in done:
                    continue
                t0 = time.time()
                g0, gn, P, dim = samples(n, alpha, T, DRAWS[n], seed=1000 + n)
                var = float(g0.var(ddof=1))
                ci = boot_var(g0, rng_b)
                # cross-check against the S25 instrument itself on the same draws
                ref = QP.measure(n, LAYERS, alpha, T, DRAWS[n], seed=1000 + n, init_sd=INIT_SD)
                rows.append(dict(n=n, alpha=alpha, T=T, P=P, dim=dim, draws=DRAWS[n],
                                 var_g0=var, var_g0_ci95=ci,
                                 mean_sq_per_param=float(gn.mean() / P),
                                 mean_sq_per_param_ci95=[float(np.percentile(
                                     [gn[rng_b.integers(0, len(gn), len(gn))].mean() / P
                                      for _ in range(2000)], q)) for q in (2.5, 97.5)],
                                 s25_measure_var_g0=ref["var_g0"],
                                 rel_dev_vs_s25_measure=abs(ref["var_g0"] - var) / max(var, 1e-300),
                                 seconds=time.time() - t0))
                print(f"  T={T} alpha={alpha:<5} n={n:2d} dim={dim:5d} P={P:2d} draws={DRAWS[n]:3d} "
                      f"Var[dF/dth0]={var:.4e} ci[{ci[0]:.2e},{ci[1]:.2e}] "
                      f"(S25 instrument {ref['var_g0']:.4e})  {time.time()-t0:5.1f}s", flush=True)
                with open(ck, "w", encoding="utf-8") as fh:
                    json.dump(dict(rows=rows), fh)

    # slopes of log2 Var against n, per (alpha, T), with a bootstrap CI over the draws' variance
    slopes = {}
    for T in (0.0, 0.3):
        for alpha in ALPHAS:
            rr = sorted([r for r in rows if r["alpha"] == alpha and r["T"] == T], key=lambda r: r["n"])
            ns = np.array([r["n"] for r in rr], float)
            v = np.array([r["var_g0"] for r in rr])
            lo = np.array([r["var_g0_ci95"][0] for r in rr]); hi = np.array([r["var_g0_ci95"][1] for r in rr])
            ok = v > 0
            slope = float(np.polyfit(ns[ok], np.log2(v[ok]), 1)[0]) if ok.sum() > 2 else None
            # slope CI: refit on log2 of values resampled uniformly inside each point's CI
            sl = []
            for _ in range(2000):
                vv = np.exp(np.log(np.maximum(lo, 1e-300)) + rng_b.random(len(v)) * (np.log(np.maximum(hi, 1e-300)) - np.log(np.maximum(lo, 1e-300))))
                if (vv > 0).sum() > 2:
                    sl.append(float(np.polyfit(ns[vv > 0], np.log2(vv[vv > 0]), 1)[0]))
            ratio = float(v[ok][-1] / v[ok][0]) if ok.sum() > 1 else None
            slopes[f"alpha={alpha}_T={T}"] = dict(
                log2_var_slope_per_qubit=slope,
                slope_ci95=[float(np.percentile(sl, 2.5)), float(np.percentile(sl, 97.5))] if sl else None,
                var_ratio_last_over_first=ratio, n_first=int(ns[ok][0]), n_last=int(ns[ok][-1]),
                reading=("decays (CI excludes 0)" if slope is not None and sl and np.percentile(sl, 97.5) < 0 else
                         "no decay resolved (CI includes 0)" if slope is not None else "undefined"))
    out = dict(protocol=__doc__, alphas=ALPHAS, ns=NS, layers=LAYERS, draws=DRAWS, init_sd=INIT_SD,
               rows=rows, slopes=slopes,
               environment=dict(python=sys.version.split()[0], numpy=np.__version__))
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    # ---- Graph 2 (T = 0, the CVaR objective; one line per alpha)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 6.2), dpi=190)
    fig.patch.set_facecolor("white")
    colors = {1.0: "#4c72b0", 0.5: "#55a868", 0.25: "#dd8452", 0.10: "#c44e52"}
    for alpha in ALPHAS:
        rr = sorted([r for r in rows if r["alpha"] == alpha and r["T"] == 0.0], key=lambda r: r["n"])
        ns = np.array([r["n"] for r in rr]); v = np.array([r["var_g0"] for r in rr])
        lo = v - np.array([r["var_g0_ci95"][0] for r in rr]); hi = np.array([r["var_g0_ci95"][1] for r in rr]) - v
        s = slopes[f"alpha={alpha}_T=0.0"]
        lab = f"alpha = {alpha}" + (" (linear cost, <psi|H|psi>)" if alpha == 1.0 else "")
        lab += f"; slope {s['log2_var_slope_per_qubit']:+.2f} log2/qubit [{s['slope_ci95'][0]:+.2f}, {s['slope_ci95'][1]:+.2f}]"
        ax.errorbar(ns, v, yerr=[np.maximum(lo, 0), np.maximum(hi, 0)], marker="o", capsize=3,
                    color=colors[alpha], label=lab)
    ax.set_yscale("log")
    ax.axvline(7, color="grey", ls=":", lw=1)
    ax.annotate("deployed n = 7", xy=(7, 1), xycoords=("data", "axes fraction"), xytext=(3, -12),
                textcoords="offset points", fontsize=8, color="grey", va="top")
    ax.set_xlabel("number of qubits n (register of 2^n candidates; P = 3n parameters, depth 3)")
    ax.set_ylabel("Var over theta of dCVaR / dtheta_0 (log scale)")
    ax.set_title("Gradient variance of the deployed CVaR objective versus qubit count\n"
                 "exact parameter-shift gradient; theta ~ N(0, 0.6^2); deployed rank-ladder H; T = 0\n"
                 "error bars = bootstrap 95% CI over the theta draws; slopes = least squares on log2 Var, with CI",
                 fontsize=9.5)
    ax.set_xticks(NS)
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1, frameon=False)
    ax.grid(which="both", alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, facecolor="white")
    print("\nSLOPES (log2 Var per qubit, T = 0):")
    for alpha in ALPHAS:
        s = slopes[f"alpha={alpha}_T=0.0"]
        print(f"  alpha={alpha:<5} slope {s['log2_var_slope_per_qubit']:+.3f} ci {s['slope_ci95']}  "
              f"Var(n={s['n_last']})/Var(n={s['n_first']}) = {s['var_ratio_last_over_first']:.3g}  -> {s['reading']}")
    print("data:", OUT_JSON)
    print("graph:", OUT_PNG)


if __name__ == "__main__":
    main()
