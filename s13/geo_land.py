"""SPRINT 13 GEO -- EXPERIMENT 8: landscape geometry around matched states.

THREE QUESTIONS, ONE PER SECTION
    1. **Curvature.**  Hessian of ``CVaR_alpha(theta)`` by central differences on the EXACT
       parameter-shift gradient, at (a) random theta and (b) each model's own optimiser
       endpoint reached from the SAME theta0.  Reported: eigenvalue spectrum, fraction of
       negative eigenvalues, spectral norm, trace, and the condition number of |H|.
       Curvature carries the energy's units, so every comparison is on IQR-standardised
       energies and the standardiser is recorded.
    2. **Multimodality.**  M independent Adam descents from random starts; endpoints
       clustered by final objective value.  Reported: number of distinct minima at two
       tolerances, the spread of final objective values, and the fraction of starts
       reaching within 1 % of the best endpoint found.
    3. **Does parameter distance mean structural distance?**  For random theta pairs,
       ``||dtheta||_2``, the Fubini-Study distance ``sqrt(dtheta' g dtheta)`` at the
       midpoint, and the CA-RMSD between the two states' argmax-probability structures.
       If the parameterisation is badly matched to the manifold these decouple, and that
       is the concrete content of "QNG should help".

    python -m s13.geo_land
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                        # noqa: E402
from s12 import instrument as I                        # noqa: E402
from core.quantum import Adam                          # noqa: E402

ALPHA = 0.25
LAYERS = 3
CELLS = [(p, 5, 4) for p in ("1A13", "1A1P", "2BFI")] + [(p, 6, 4) for p in ("1A13", "1A1P", "2BFI")]
N_START = 12
STEPS = 150
LR = 0.10


def std_energy(E):
    """IQR-standardised energy (affine: rankings and argmin unchanged, units removed)."""
    E = G.clean(E)
    s = G.scales(E)["iqr"]
    return (E - np.median(E)) / (s if s > 0 else 1.0), s


def hessian(E, circ, theta, alpha, h=1e-3):
    P = theta.size
    H = np.empty((P, P))
    for i in range(P):
        tp = theta.copy(); tp[i] += h
        tm = theta.copy(); tm[i] -= h
        H[:, i] = (G.cvar_grad(E, circ, tp, alpha)
                   - G.cvar_grad(E, circ, tm, alpha)) / (2 * h)
    H = (H + H.T) / 2.0
    w = np.linalg.eigvalsh(H)
    aw = np.abs(w)
    nz = aw[aw > 1e-12 * max(aw.max(), 1e-300)]
    return {"eigs_top6": [float(x) for x in w[::-1][:6]],
            "eigs_bot6": [float(x) for x in w[:6]],
            "lam_max": float(w.max()), "lam_min": float(w.min()),
            "trace": float(w.sum()), "spectral_norm": float(aw.max()),
            "frac_negative": float((w < 0).mean()),
            "frac_near_zero_1e-3_rel": float((aw < 1e-3 * aw.max()).mean()),
            # A completely flat Hessian (AMBER at random theta on some cells) has every
            # eigenvalue at fp noise, so the "smallest non-negligible" set can be empty.
            "cond_abs": float(aw.max() / max(float(nz.min()) if nz.size else 0.0, 1e-300)),
            "neg_mass_frac": float(aw[w < 0].sum() / max(aw.sum(), 1e-300))}


def descend(E, circ, theta0, alpha, steps=STEPS, lr=LR):
    th = np.array(theta0, float)
    opt = Adam(th.size, lr=lr)
    for _ in range(steps):
        th = opt.step(th, G.cvar_grad(E, circ, th, alpha))
    return th, G.cvar_value(E, circ, th, alpha)


def multimodality(E, circ, alpha, n_start=N_START, seed=0):
    rng = np.random.default_rng(seed)
    ends, vals = [], []
    for _ in range(n_start):
        th0 = rng.uniform(-np.pi, np.pi, circ.n_params())
        th, v = descend(E, circ, th0, alpha)
        ends.append(th); vals.append(v)
    vals = np.array(vals)
    best = vals.min()
    span = max(vals.max() - best, 1e-12)

    def n_distinct(tol):
        s = np.sort(vals)
        c = 1
        for a in range(1, len(s)):
            if s[a] - s[a - 1] > tol * span:
                c += 1
        return c

    return {"n_starts": n_start, "val_best": float(best), "val_worst": float(vals.max()),
            "val_mean": float(vals.mean()), "val_sd": float(vals.std()),
            "n_distinct_tol0.01": n_distinct(0.01),
            "n_distinct_tol0.05": n_distinct(0.05),
            "frac_within_1pct_of_best": float((vals <= best + 0.01 * span).mean()),
            "frac_within_5pct_of_best": float((vals <= best + 0.05 * span).mean())}, ends


def geometry_vs_structure(circ, pdb, L, k, n_pairs=64, seed=0):
    from s13.geo_cvar import ca_coords
    W, nat = ca_coords(pdb, L, k)
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    d2, dfs, drm = [], [], []
    for _ in range(n_pairs):
        a = rng.uniform(-np.pi, np.pi, P)
        b = rng.uniform(-np.pi, np.pi, P)
        g = G.fs_metric(circ, (a + b) / 2)
        d = b - a
        d2.append(float(np.linalg.norm(d)))
        dfs.append(float(np.sqrt(max(d @ g @ d, 0.0))))
        ja = int(np.argmax(circ.probs(a))); jb = int(np.argmax(circ.probs(b)))
        drm.append(float(I.ca_rmsd(W[ja], W[jb])))
    d2, dfs, drm = np.array(d2), np.array(dfs), np.array(drm)
    return {"n_pairs": n_pairs,
            "corr_l2_vs_structural": float(np.corrcoef(d2, drm)[0, 1]),
            "corr_fs_vs_structural": float(np.corrcoef(dfs, drm)[0, 1]),
            "spearman_l2_vs_structural": G.spearman(d2, drm),
            "spearman_fs_vs_structural": G.spearman(dfs, drm),
            "corr_l2_vs_fs": float(np.corrcoef(d2, dfs)[0, 1]),
            "mean_l2": float(d2.mean()), "mean_fs": float(dfs.mean()),
            "mean_structural_rmsd": float(drm.mean())}


def main():
    rows = []
    for pdb, L, k in CELLS:
        if not os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")):
            continue
        n = int(round(np.log2(k))) * L
        circ = G.circuit(n, layers=LAYERS)
        G.gate(tries=1)
        t0 = time.time()
        row = {"pdb": pdb, "L": L, "k": k, "n_qubits": n, "layers": LAYERS,
               "alpha": ALPHA, "P": circ.n_params(),
               "geometry_vs_structure": geometry_vs_structure(circ, pdb, L, k)}
        for v in G.VARIANTS:
            Ez, s = std_energy(G.variant_table(pdb, L, k, v))
            mm, _ends = multimodality(Ez, circ, ALPHA)
            hs_rand, hs_end = [], []
            for seed in (0, 1):
                th0 = G.init_theta(circ.n_params(), seed)
                hs_rand.append(hessian(Ez, circ, th0, ALPHA))
                thE, _ = descend(Ez, circ, th0, ALPHA)
                hs_end.append(hessian(Ez, circ, thE, ALPHA))
            row[v] = {"iqr_standardiser": s, "multimodality": mm,
                      "hessian_random_theta": hs_rand,
                      "hessian_at_endpoint": hs_end}
        row["wall_s"] = time.time() - t0
        rows.append(row)
        print("  %-5s L=%-2d n=%-2d | neg-eig frac (endpoint) leg %.3f amb %.3f soft %.3f "
              "| distinct minima leg %d amb %d | %.0fs"
              % (pdb, L, n,
                 np.mean([h["frac_negative"] for h in row["legacy"]["hessian_at_endpoint"]]),
                 np.mean([h["frac_negative"] for h in row["amber"]["hessian_at_endpoint"]]),
                 np.mean([h["frac_negative"]
                          for h in row["amber_soft"]["hessian_at_endpoint"]]),
                 row["legacy"]["multimodality"]["n_distinct_tol0.01"],
                 row["amber"]["multimodality"]["n_distinct_tol0.01"],
                 row["wall_s"]), flush=True)
        G.write("geo_land", {"what": "landscape curvature, multimodality and the "
                                     "parameter-vs-structure distance relation",
                             "rows": rows})
    G.write("geo_land", {"what": "landscape curvature, multimodality and the "
                                 "parameter-vs-structure distance relation", "rows": rows})


if __name__ == "__main__":
    main()
