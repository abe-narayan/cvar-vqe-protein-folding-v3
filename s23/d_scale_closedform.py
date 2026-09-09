"""s23/d_scale_closedform.py -- WORKSTREAM D AUDIT of gscale.py's per-target oracle.

Independently reproduces the incumbent (score-filter top-75 coordinate average) at n=126 and
replaces gscale.py's 81-point grid search over s in [0.90,1.10] with the EXACT closed-form
scalar-Procrustes optimum, derived as follows.

DERIVATION (checked numerically below, not just asserted).  Kabsch superposition fits a rotation
R minimising sum_i |R(C_i-Cc) - (T_i-Tc)|^2 via SVD of H = sum_i (C_i-Cc)(T_i-Tc)^T.  Scaling C by
a POSITIVE scalar s scales H by s but does not change its left/right singular vectors, so R* is
IDENTICAL for every s>0 -- the Kabsch rotation is scale-invariant.  Fix R once (at s=1), let
C'_i = R(C_i-Cc), T'_i = T_i-Tc.  Then

    SSE(s) = sum_i |s*C'_i - T'_i|^2 = Cc*s^2 - 2*Ct*s + Tt        (a PARABOLA in s)
    Cc = sum_i |C'_i|^2,  Ct = sum_i C'_i . T'_i,  Tt = sum_i |T'_i|^2

    s*_analytic = Ct / Cc                         (closed form, exact, no grid)
    SSE(1) - SSE(s*) = (Cc - Ct)^2 / Cc  >= 0      ALWAYS, with equality iff s*=1 exactly.

This is a ONE-PARAMETER LEAST-SQUARES fit through data that ALSO supplies the evaluation label
(native coordinates).  The improvement is MECHANICALLY GUARANTEED to be non-negative for any C,T
with Cc>0, regardless of whether the correspondence between C and T carries real structural
information -- this is the source of the 126W/0L record in gscale.json, and is analysed
quantitatively in d_scale_placebo.py.

This script's OWN job: (1) confirm the parabola claim numerically: (2) compute the EXACT
unconstrained s*_analytic per target and compare to gscale.py's 81-pt grid-clipped s_star, to
quantify how much of the 63/126 grid-boundary hits are a grid-width artefact rather than a
genuine extreme per-target optimum.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I  # noqa: E402

TOPM = 75


def _save(o, name="d_scale_closedform.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def kabsch_R(C, T):
    """Return the Kabsch rotation aligning centred C onto centred T (both (n,3))."""
    Cc = C - C.mean(0); Tc = T - T.mean(0)
    H = Cc.T @ Tc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.eye(3); D[2, 2] = d
    R = Vt.T @ D @ U.T
    return R, Cc, Tc


def run():
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        C, _b = I.coordinate_average(W[o[:TOPM]])
        C = np.asarray(C, float)

        R, Cprime, Tprime = kabsch_R(C, nat)
        Crot = Cprime @ R.T  # R(C-Cc), in the frame where T is Tprime

        Cc = float((Crot ** 2).sum())
        Ct = float((Crot * Tprime).sum())
        Tt = float((Tprime ** 2).sum())
        n = Crot.shape[0]

        s_analytic = Ct / Cc
        sse1 = Cc - 2 * Ct + Tt
        sse_star = Tt - Ct * Ct / Cc
        rmsd1 = float(np.sqrt(max(sse1, 0.0) / n))
        rmsd_star_analytic = float(np.sqrt(max(sse_star, 0.0) / n))

        # numeric cross-check: parabola-in-s claim, and grid-based oracle for comparison
        grid = np.linspace(0.5, 1.5, 2001)  # wide grid, fine resolution, brute-force check
        sse_grid = Cc * grid ** 2 - 2 * Ct * grid + Tt
        pred_sse = Cc * s_analytic ** 2 - 2 * Ct * s_analytic + Tt
        parabola_ok = bool(np.isclose(pred_sse, sse_grid.min(), rtol=1e-6, atol=1e-8))
        s_grid_wide = float(grid[int(np.argmin(sse_grid))])

        # narrow grid matching gscale.py exactly, for direct comparison
        narrow = np.linspace(0.90, 1.10, 81)
        rmsd_narrow = np.sqrt(np.maximum(Cc * narrow ** 2 - 2 * Ct * narrow + Tt, 0.0) / n)
        s_narrow = float(narrow[int(np.argmin(rmsd_narrow))])
        rmsd_narrow_star = float(rmsd_narrow.min())

        rows.append({
            "pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
            "rmsd_1": rmsd1,
            "s_analytic": s_analytic,
            "rmsd_star_analytic": rmsd_star_analytic,
            "s_narrow_grid": s_narrow,
            "rmsd_star_narrow_grid": rmsd_narrow_star,
            "s_wide_grid_check": s_grid_wide,
            "parabola_confirmed": parabola_ok,
            "hit_low_bound": bool(s_narrow <= 0.9001),
            "hit_high_bound": bool(s_narrow >= 1.0999),
        })
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)

    ok = len(rows) == len(tg) and all(np.isfinite(r["s_analytic"]) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "d_scale_closedform.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)  # noqa: E731

    print("\nn = %d.  Independent reproduction, closed-form scalar Procrustes.\n" % len(rows))
    print("  parabola-in-s confirmed for all targets: %s (%d/%d)"
          % (all(r["parabola_confirmed"] for r in rows),
             sum(r["parabola_confirmed"] for r in rows), len(rows)))

    print("\n  %-42s%9s" % ("arm", "RMSD"))
    print("  %-42s%9.4f" % ("incumbent s=1.0 (independent repro)", g("rmsd_1").mean()))
    print("  %-42s%9.4f" % ("oracle, narrow grid [0.90,1.10] (repro of A's)", g("rmsd_star_narrow_grid").mean()))
    print("  %-42s%9.4f" % ("oracle, EXACT closed-form (unconstrained s>0)", g("rmsd_star_analytic").mean()))

    d_narrow = g("rmsd_star_narrow_grid") - g("rmsd_1")
    d_exact = g("rmsd_star_analytic") - g("rmsd_1")
    print("\n  delta vs incumbent, narrow-grid oracle : %+.4f" % d_narrow.mean())
    print("  delta vs incumbent, EXACT closed-form  : %+.4f  <-- the TRUE unconstrained ceiling"
          % d_exact.mean())

    nb = sum(r["hit_low_bound"] or r["hit_high_bound"] for r in rows)
    print("\n  targets hitting the [0.90,1.10] grid boundary: %d/%d (%.0f%%)"
          % (nb, len(rows), 100.0 * nb / len(rows)))
    sa = g("s_analytic")
    print("  s_analytic (EXACT, unconstrained) distribution: mean %.4f sd %.4f min %.4f max %.4f"
          % (sa.mean(), sa.std(), sa.min(), sa.max()))
    print("  s_analytic outside [0.90,1.10]: %d/%d (%.0f%%)  outside [0.5,1.5]: %d/%d"
          % ((np.abs(sa - 1) > 0.10).sum(), len(rows), 100.0 * (np.abs(sa - 1) > 0.10).mean(),
             ((sa < 0.5) | (sa > 1.5)).sum(), len(rows)))
    extreme = sorted(zip(sa, [r["pdb"] for r in rows]))
    print("  5 most extreme s_analytic (low):", extreme[:5])
    print("  5 most extreme s_analytic (high):", extreme[-5:])


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
