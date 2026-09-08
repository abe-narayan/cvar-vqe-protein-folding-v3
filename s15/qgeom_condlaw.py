"""SPRINT 15 / QGEOM -- PART B(iii): the POSITIVE form of the conditioning result.

B2 shows the conditioned point is geometrically indistinguishable from an entropy-matched
scramble of its own prior on every statistic measured.  That is the negative.  The positive
statement -- and the one worth a paper sentence -- is that the geometry is a FUNCTION OF THE
ENTROPY of whatever distribution the state was prepared for.  This module tests that
directly, pooled over every arm and every target:

    (i)  correlate each metric statistic with H(p) across all 7 arms x 9 targets x 3 seeds
    (ii) ask how much variance is left for ARM IDENTITY once H(p) is regressed out.
         If target-specificity mattered, the arm label would carry residual variance.
    (iii) the same for the ORIGINAL prior's entropy, which is what an information-ladder
         story would want the geometry to depend on.

    python -m s15.qgeom_condlaw
"""
from __future__ import annotations

import json
import os

import numpy as np

from s15 import qgeom_lib as G

TAG = "condlaw"
STATS = ("cond", "eff_rank_frac", "log_pseudo_det", "max_dev_from_I4", "grad_var",
         "grad_share_bottom_decile", "max_p")
ARMS = ("uniform", "random", "COND", "scram_state", "scram_res", "cross", "dirichlet")


def collect(key="B1_blockL2"):
    p = os.path.join(G.RESULTS, "qgeom_cond.json")
    d = json.load(open(p))[key]
    rows = []
    for cell, v in d.items():
        if "|s" not in cell:
            continue
        for a in ARMS:
            g = v["arms"][a]["geom"]
            r = {"cell": cell, "arm": a, "H": g["entropy_bits"],
                 "prior_H": v["prior_entropy_bits"]}
            for s in STATS:
                r[s] = np.log10(g[s]) if s == "cond" else g.get(s, np.nan)
            rows.append(r)
    return rows


def r2_of_labels(y, labels):
    """Fraction of variance in y explained by a categorical label (one-way ANOVA R^2)."""
    y = np.asarray(y, float)
    tot = float(((y - y.mean()) ** 2).sum())
    if tot <= 0:
        return 0.0
    ss = 0.0
    for lb in set(labels):
        m = np.array([l == lb for l in labels])
        ss += float(((y[m] - y[m].mean()) ** 2).sum())
    return float(1.0 - ss / tot)


def main(key="B1_blockL2"):
    rows = collect(key)
    print("=" * 104)
    print(f"B3. IS THE GEOMETRY A FUNCTION OF THE ENTROPY OF THE PREPARED DISTRIBUTION?")
    print(f"    {len(rows)} points = {len(set(r['cell'] for r in rows))} cells x "
          f"{len(ARMS)} arms, ansatz key {key}")
    print("=" * 104)
    H = np.array([r["H"] for r in rows])
    labels = [r["arm"] for r in rows]
    lab_noU = [r["arm"] for r in rows if r["arm"] != "uniform"]
    print(f"{'statistic':26s} {'rho(stat,H)':>12s} {'R2 lin(H)':>10s} "
          f"{'R2 arm label':>13s} {'R2 arm | H resid':>17s} {'R2 arm|H, no uniform':>21s} "
          f"{'R2 arm|H, FITTED':>16s}")
    out = {}
    for s in STATS:
        y = np.array([r[s] for r in rows], float)
        ok = np.isfinite(y) & np.isfinite(H)
        yy, hh = y[ok], H[ok]
        lb = [labels[i] for i in range(len(rows)) if ok[i]]
        rho = G.spearman(yy, hh)
        b = np.polyfit(hh, yy, 1)
        res = yy - np.polyval(b, hh)
        r2lin = 1.0 - float((res ** 2).sum() / max(((yy - yy.mean()) ** 2).sum(), 1e-30))
        r2arm = r2_of_labels(yy, lb)
        r2res = r2_of_labels(res, lb)
        m2 = np.array([l != "uniform" for l in lb])
        r2res_noU = r2_of_labels(res[m2], [l for l in lb if l != "uniform"])
        # the cleanest test: only the FITTED arms, all of which were produced by the
        # identical optimiser on a product prior and differ ONLY in whose prior it was
        fit = ("COND", "scram_state", "scram_res", "cross", "dirichlet")
        m3 = np.array([l in fit for l in lb])
        r2res_fit = r2_of_labels(res[m3], [l for l in lb if l in fit])
        out[s] = {"spearman_with_H": rho, "r2_linear_in_H": r2lin,
                  "r2_arm_label": r2arm, "r2_arm_on_H_residual": r2res,
                  "r2_arm_on_H_residual_excl_uniform": r2res_noU,
                  "r2_arm_on_H_residual_fitted_arms_only": r2res_fit,
                  "n": int(ok.sum()), "n_fitted": int(m3.sum())}
        print(f"{s:26s} {rho:+12.3f} {r2lin:10.3f} {r2arm:13.3f} {r2res:17.3f} "
              f"{r2res_noU:21.3f} {r2res_fit:16.3f}")
    print()
    print("  `R2 arm | H resid` is the variance the ARM LABEL still explains AFTER the")
    print("  entropy of the prepared distribution has been regressed out.  A small value")
    print("  means the arm -- i.e. whether the conditioning was the RIGHT target's prior,")
    print("  a scramble of it, another target's, or synthetic -- carries no geometry once")
    print("  entropy is known.  The last column drops `uniform`, whose entropy is an")
    print("  extreme leverage point, so the test cannot be carried by one arm.")
    G.ck(TAG, key, out)

    # the same against the PRIOR's entropy rather than the realised state entropy
    print()
    print("  Against the ORIGINAL prior's entropy (COND / scram arms only):")
    sel = [r for r in rows if r["arm"] in ("COND", "scram_state", "scram_res")]
    ph = np.array([r["prior_H"] for r in sel])
    for s in ("cond", "log_pseudo_det"):
        y = np.array([r[s] for r in sel], float)
        print(f"    rho({s}, H(prior)) = {G.spearman(y, ph):+.3f}   "
              f"rho({s}, H(state)) = "
              f"{G.spearman(y, np.array([r['H'] for r in sel])):+.3f}")
    return out


if __name__ == "__main__":
    main()
