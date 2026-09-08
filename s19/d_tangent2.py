"""AGENT D, Sprint 19 -- T3 CORRECTED.  `s19/d_tangent.py` has a defect and this fixes it.

THE DEFECT.  `frac_in` projected with `np.linalg.qr(J_w)`.  Reduced QR of an (npairs x 2n) matrix
returns **2n** orthonormal columns even when J_w is rank-deficient -- and it IS rank-deficient,
because the map (phi, psi) -> CA trace is not injective (an ideal-geometry CA trace of n residues
has at most 2n-5 internal coordinates).  So every in-tangent fraction was computed against a
2n-dimensional subspace instead of the true column space.  The tell is in the output: the EMPIRICAL
random-direction null came out at **0.417**, which is 2n/npairs, not rank(J_w)/npairs = 0.329.

Contrasts against `fin_random` were still valid, because the random control shared the inflated
subspace -- but no absolute number from that run is quotable, and this module re-measures with an
SVD-truncated projection onto the genuine column space.

It also reuses AGENT A's cached starts (`s19/cache/start_<pdb>.npz`), which is what `d_tangent.py`
should have done: the projection of the coordinate average is deterministic, identical for every
arm and every lane, and costs 2-4 s per target to rebuild.

Run:  python -m s19.d_tangent2
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
from s19 import a_fit as F                 # noqa: E402

OUT = os.path.join(HERE, "results", "D_T3_tangent2")
os.makedirs(OUT, exist_ok=True)
EPS = 1e-5
RTOL = 1e-6                                #: singular-value cutoff, relative to the largest


def jac_basis(phi, psi, i, j, inv):
    """Orthonormal basis of the TRUE column space of the whitened Jacobian, by SVD."""
    n = len(phi)
    th = np.concatenate([phi, psi])

    def dist(t):
        ca = I.build_ca(t[:n], t[n:])
        return np.sqrt(((ca[i] - ca[j]) ** 2).sum(1))

    J = np.empty((len(i), 2 * n))
    for k in range(2 * n):
        tp, tm = th.copy(), th.copy()
        tp[k] += EPS
        tm[k] -= EPS
        J[:, k] = (dist(tp) - dist(tm)) / (2 * EPS)
    J_w = J * inv[:, None]
    U, s, _ = np.linalg.svd(J_w, full_matrices=False)
    r = int((s > RTOL * s[0]).sum())
    return U[:, :r], r


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
        rng = SD.stable_rng(pdb, "s19D_tangent2")
        inv = 1.0 / sd
        phi0, psi0, _avg = F.start(pdb, seq, fold)            # AGENT A's cached start
        U, r = jac_basis(np.asarray(phi0, float), np.asarray(psi0, float), i, j, inv)

        def frac(v):
            vw = v * inv
            p = U.T @ vw
            return float((p * p).sum() / max(float((vw * vw).sum()), 1e-30))

        res = dhat - dtrue
        rms = float(np.sqrt((res ** 2).mean()))
        W = np.asarray(AV.top75_windows(pdb)[0], float)
        DW = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1)
        rw_all = DW - dtrue[None]
        k = int(np.argmin(np.abs(np.sqrt((rw_all ** 2).mean(1)) - rms)))
        fields = {
            "real": res,
            "signflip": res * rng.choice([-1.0, 1.0], size=len(res)),
            "shuffled": rng.permutation(res) * rng.choice([-1.0, 1.0], size=len(res)),
            "iso": rng.standard_normal(len(res)) * rms,
            "coherent": rw_all[k],
        }
        e = {"pdb": pdb, "n": n, "fold": fold, "npairs": len(i), "rank": r,
             "null_rank": r / len(i), "resid_rms": rms}
        for a, v in fields.items():
            e["fin_" + a] = frac(v)
        #: TWO nulls, and the difference between them is itself informative.
        #: `rand_white`  isotropic in the WHITENED space -- the pure geometric null, E = rank/npairs
        #: `rand_raw`    isotropic in RAW distance space then whitened -- what `iso`/`shuffled`
        #:               actually are, and NOT isotropic after whitening
        e["fin_rand_white"] = float(np.mean([
            float(((U.T @ (z := rng.standard_normal(len(i)))) ** 2).sum() / (z * z).sum())
            for _ in range(30)]))
        e["fin_rand_raw"] = float(np.mean([frac(rng.standard_normal(len(i)) * rms)
                                           for _ in range(30)]))
        rows.append(e)
        if (c + 1) % 25 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}", flush=True)

    rng2 = SD.stable_rng("D", "T3b")

    def boot(v, B=4000):
        v = np.asarray(v, float)
        b = v[rng2.integers(0, len(v), size=(B, len(v)))].mean(1)
        return float(v.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

    g = lambda k: np.asarray([r[k] for r in rows], float)      # noqa: E731
    print(f"\nn = {len(rows)}.  SVD-truncated projection onto the TRUE column space of J_w.\n")
    print(f"  {'field':<28}{'in-tangent':>12}{'[95% CI]':>20}")
    for a in ("real", "signflip", "shuffled", "iso", "coherent",
              "rand_white", "rand_raw"):
        m, lo, hi = boot(g("fin_" + a))
        tag = "   <- the geometric null" if a == "rand_white" else (
            "   <- the null iso/shuffled actually live at" if a == "rand_raw" else "")
        print(f"  {a:<28}{m:>12.3f}   [{lo:.3f},{hi:.3f}]{tag}")
    m, lo, hi = boot(g("null_rank"))
    print(f"  {'rank(J_w)/npairs (PREDICTED)':<28}{m:>12.3f}   [{lo:.3f},{hi:.3f}]")
    print(f"  mean rank {g('rank').mean():.2f},  mean 2n-5 = "
          f"{(2*g('n')-5).mean():.2f},  mean 2n = {(2*g('n')).mean():.2f}")
    for lab, a, b in [("real - signflip (MATCHED whitened magnitudes)", "real", "signflip"),
                      ("real - rand_white", "real", "rand_white"),
                      ("iso - rand_white  (the whitening anisotropy)", "iso", "rand_white")]:
        m, lo, hi = boot(g("fin_" + a) - g("fin_" + b))
        print(f"  {lab:<46}{m:+.3f} [{lo:+.3f},{hi:+.3f}]  W/L "
              f"{int((g('fin_'+a) > g('fin_'+b)).sum())}/"
              f"{int((g('fin_'+a) < g('fin_'+b)).sum())}")

    try:
        A = {r["pdb"]: r for r in json.load(open(os.path.join(HERE, "results",
                                                             "a_coh.json")))["rows"]}
        from scipy.stats import spearmanr
        gap = np.asarray([A[r["pdb"]]["real"] - A[r["pdb"]]["signflip"] for r in rows])
        exc = g("fin_real") - g("fin_signflip")
        for lab, v in [("in-tangent(real) - in-tangent(signflip)", exc),
                       ("in-tangent(real) - rand_white", g("fin_real") - g("fin_rand_white")),
                       ("in-tangent(real), raw", g("fin_real"))]:
            r_ = spearmanr(v, gap)
            print(f"  rho({lab:<42}, per-target gap) = {r_.statistic:+.3f}  p={r_.pvalue:.3g}")
        bs = np.array([spearmanr(exc[k], gap[k]).statistic
                       for k in rng2.integers(0, len(gap), size=(2000, len(gap)))])
        print(f"  95% CI = [{np.percentile(bs,2.5):+.3f}, {np.percentile(bs,97.5):+.3f}]")
        kp = np.asarray([A[r["pdb"]]["real_kappa"] for r in rows])
        print(f"  pearson across the 5 arm means (start tangent vs A's realised kappa) = "
              f"{np.corrcoef([g('fin_'+a).mean() for a in ('real','signflip','shuffled','iso','coherent')], [np.mean([A[r['pdb']][a+'_kappa'] for r in rows]) for a in ('real','signflip','shuffled','iso','coherent')])[0,1]:+.3f}")
        print(f"  per-target rho(in-tangent real, A's kappa real) = "
              f"{spearmanr(g('fin_real'), kp).statistic:+.3f}")
    except Exception as ex:                                    # noqa: BLE001
        print("  (a_coh cross-check unavailable:", ex, ")")

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "n_expected": len(tg)},
              open(os.path.join(OUT, "tangent2.json"), "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={len(rows)}\n")


if __name__ == "__main__":
    main()
