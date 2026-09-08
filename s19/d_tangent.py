"""AGENT D, Sprint 19 -- T3: WHY realisable error is harmful, and what the analytic null is.

AGENT A measured `kappa` -- the fraction of a distance-error field that a real ideal-geometry
structure can REALISE -- and found real 0.809, signflip 0.352, shuffled/iso 0.367, a real
alternative structure 0.990, with RMSD ordering the arms in the same order (pearson +0.984 across
the six arm means, per-target Spearman > 0 on 116/126).

THIS MODULE SUPPLIES THE MISSING PIECE: **kappa has a free analytic null.**

The weighted fit is, to first order, an orthogonal projection of the whitened residual onto the
TANGENT SPACE of the ideal-geometry manifold at the start point.  That tangent space has dimension
D = 2n - 2 (S18 B5: the objective is a function of the n-2 interior residues' torsions, and the two
terminal torsion coordinates move no CA), embedded in a pair space of dimension P = (n-1)(n-2)/2.
A GENERIC error vector therefore retains a fraction

        kappa_null  ~  D / P

of its energy under the projection, with NO information content whatsoever.  Mean over the 126
targets: (2n-2)/((n-1)(n-2)/2) = 0.381.  AGENT A's incoherent arms measure 0.352-0.367.

    **The incoherent arms are AT the analytic null.  The real field's 0.809 is 2.1x it.**

That converts "the errors are worse than random errors of the same magnitude" into a statement with
a mechanism and a scale: the fit is a DENOISER against the component of the error orthogonal to the
manifold and has NO power at all against the component inside it, and the distogram puts twice as
much of its error inside the manifold as chance would.

What this module measures directly, at the deployed start point and with the deployed weights:

  * the exact in-tangent energy fraction of the real residual, by explicit projection
  * the same for signflip / shuffled / isotropic / a real-alternative-structure residual
  * the analytic null D/P, per target
  * whether the per-target in-tangent fraction predicts the per-target Sprint-18 gap

NO OPTIMISATION IS RUN.  This is linear algebra on a 78 x 26 Jacobian per target.

Run:  python -m s19.d_tangent
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I            # noqa: E402
from s14 import avgspace as AV             # noqa: E402
from s15 import distcal as C               # noqa: E402
from s15 import seed as SD                 # noqa: E402

OUT = os.path.join(HERE, "results", "D_T3_tangent")
os.makedirs(OUT, exist_ok=True)
EPS = 1e-5


def jacobian(phi, psi, i, j):
    """d(pair distance) / d(torsion), by central differences on the EXACT builder."""
    n = len(phi)
    th = np.concatenate([phi, psi])

    def dist(t):
        ca = I.build_ca(t[:n], t[n:])
        return np.sqrt(((ca[i] - ca[j]) ** 2).sum(1))

    J = np.zeros((len(i), 2 * n))
    for k in range(2 * n):
        tp, tm = th.copy(), th.copy()
        tp[k] += EPS
        tm[k] -= EPS
        J[:, k] = (dist(tp) - dist(tm)) / (2 * EPS)
    return J


def frac_in(J_w, r_w):
    """Fraction of ||r_w||^2 lying in the column space of the whitened Jacobian."""
    q, _ = np.linalg.qr(J_w)
    p = q.T @ r_w
    d = float((r_w * r_w).sum())
    return float((p * p).sum() / max(d, 1e-30))


def main():
    tg = I.targets()
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([data[p]["fold"] for p in pdbs], int)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    rows = []
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd = d["i"], d["j"], d["sd"]
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        dtrue = d["dtrue"]
        rng = SD.stable_rng(pdb, "s19D_tangent")

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

        J = jacobian(phi0, psi0, i, j)
        inv = 1.0 / sd
        J_w = J * inv[:, None]

        r = dhat - dtrue
        rms = float(np.sqrt((r ** 2).mean()))
        DW = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1)
        rw_all = DW - dtrue[None]
        k = int(np.argmin(np.abs(np.sqrt((rw_all ** 2).mean(1)) - rms)))
        fields = {
            "real": r,
            "signflip": r * rng.choice([-1.0, 1.0], size=len(r)),
            "shuffled": rng.permutation(r) * rng.choice([-1.0, 1.0], size=len(r)),
            "iso": rng.standard_normal(len(r)) * rms,
            "coherent": rw_all[k],
        }
        npair, dof = len(i), 2 * n - 2
        e = {"pdb": pdb, "n": n, "fold": fold, "npairs": npair, "dof": dof,
             "rank_Jw": int(np.linalg.matrix_rank(J_w, tol=1e-8)),
             "null_DP": dof / npair, "resid_rms": rms}
        for a, v in fields.items():
            e["fin_" + a] = frac_in(J_w, v * inv)
        # a genuinely random direction in pair space -- the empirical version of D/P
        e["fin_random"] = float(np.mean([frac_in(J_w, rng.standard_normal(npair))
                                         for _ in range(20)]))
        rows.append(e)
        if (c + 1) % 25 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}", flush=True)

    rng = SD.stable_rng("D", "T3")

    def boot(v, B=4000):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        b = v[rng.integers(0, len(v), size=(B, len(v)))].mean(1)
        return float(v.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

    g = lambda k: np.asarray([r[k] for r in rows], float)      # noqa: E731
    print(f"\nn = {len(rows)} targets.  In-tangent energy fraction of the WHITENED residual at the "
          f"deployed start point.\n")
    print(f"  {'field':<26}{'in-tangent':>12}{'[95% CI]':>20}")
    for a in ("real", "signflip", "shuffled", "iso", "coherent", "random"):
        m, lo, hi = boot(g("fin_" + a))
        print(f"  {a:<26}{m:>12.3f}   [{lo:.3f},{hi:.3f}]")
    m, lo, hi = boot(g("null_DP"))
    print(f"  {'ANALYTIC NULL  D/P':<26}{m:>12.3f}   [{lo:.3f},{hi:.3f}]")
    print(f"  rank(J_w) mean {g('rank_Jw').mean():.2f} vs D = 2n-2 mean {g('dof').mean():.2f}"
          f"   (rank/npairs = {(g('rank_Jw')/g('npairs')).mean():.3f})")

    m, lo, hi = boot(g("fin_real") - g("fin_random"))
    print(f"\n  real MINUS a random direction in the same space: {m:+.3f} [{lo:+.3f},{hi:+.3f}]"
          f"   W/L {int((g('fin_real') > g('fin_random')).sum())}/"
          f"{int((g('fin_real') < g('fin_random')).sum())}")

    # does the per-target in-tangent excess predict the per-target Sprint-18 gap?
    try:
        A = {r["pdb"]: r for r in json.load(open(os.path.join(ROOT, "s19", "results",
                                                              "a_coh.json")))["rows"]}
        gap = np.asarray([A[r["pdb"]]["real"] - A[r["pdb"]]["signflip"] for r in rows], float)
        exc = g("fin_real") - g("fin_random")
        kap = np.asarray([A[r["pdb"]]["real_kappa"] for r in rows], float)
        from scipy.stats import spearmanr
        print(f"\n  per-target Spearman(in-tangent excess, real - signflip gap) = "
              f"{spearmanr(exc, gap).statistic:+.3f}")
        print(f"  per-target Spearman(in-tangent excess, A's kappa)            = "
              f"{spearmanr(exc, kap).statistic:+.3f}")
        print(f"  per-target Spearman(A's kappa,          gap)                 = "
              f"{spearmanr(kap, gap).statistic:+.3f}")
    except Exception as ex:                                    # noqa: BLE001
        print("  (a_coh cross-check unavailable:", ex, ")")

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "n_expected": len(tg)},
              open(os.path.join(OUT, "tangent.json"), "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={len(rows)}\n")


if __name__ == "__main__":
    main()
