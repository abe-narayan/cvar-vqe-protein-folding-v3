#!/usr/bin/env python
"""s31/s31_E2_deltas.py -- S31 lane E: build every correction vector, and the class ceiling.

E1 killed the native-free direction estimate (`s31/results/s31_E1_direction.json`).  What survives
is the coordinator's redirect, registered in `s31/PREREG_S31_E.md` AMENDMENT 1:

    Over ALL corrections orthogonal to the TRUE common mode `mu`, what is the best achievable
    endpoint?  Project the ORACLE error `y` itself onto mu's orthogonal complement and apply it.

That is the CLASS CEILING.  **ORACLE / NOT DEPLOYABLE** -- it reads `d_nat` twice over (`y` and
`mu`).  It is a price and a bound, never an achievement.

Every arm is a per-pair correction `delta` applied as production's posterior translation
`score_shift(dg, D, -delta)` -- identical to `s30/s30_P_lr.py` and `s30/s30_P_prior.py`.

    ARM                     delta                                   status
    PROD                    0                                       baseline
    FIT_N3                  yhat (S30 N3_plus_pool, R2 0.2355)      LFO-supervised, native-free at inference
    E2_PROJ_MUHAT           P_perp(yhat; mu_hat)                    ditto -- the REGISTERED arm whose gate failed
    E3_CONSTR_MUHAT         refit on mu_hat-residualised features    ditto
    CTRL_SHRINK_E2          yhat * ||P_perp(yhat;mu_hat)||/||yhat||  ditto -- matched shrinkage control
    CTRL_PROJ_RAND          P_perp(yhat; perm(mu_hat))              ditto -- matched random-direction control
    E2o_PROJ_MU             P_perp(yhat; mu)                        ORACLE / NOT DEPLOYABLE
    E3o_CONSTR_MU           refit on mu-residualised features        ORACLE / NOT DEPLOYABLE
    CTRL_SHRINK_E2o         yhat * ||P_perp(yhat;mu)||/||yhat||      ORACLE / NOT DEPLOYABLE
    ORACLE_FULL             y                                        ORACLE -- the unconstrained ceiling
    ORACLE_Y_PERP_MU        P_perp(y; mu)                            ORACLE -- THE CLASS CEILING
    ORACLE_Y_ALONG_MU       P_along(y; mu)                           ORACLE -- the discarded half
    ORACLE_Y_PERP_1MU       y residualised on [1, mu]                ORACLE -- ceiling incl. the offset
    ORACLE_Y_PERP_MUHAT     P_perp(y; mu_hat)                        ORACLE magnitude, native-free direction
    CTRL_SHRINK_ORACLE_Y    y * ||P_perp(y;mu)||/||y||               ORACLE -- matched shrinkage for the ceiling

SCOPE OF THE WORD "CEILING": `P_perp(y; mu)` is the L2-optimal approximation to the ideal
correction `y` inside the constrained class.  The endpoint is not a quadratic in delta, so this is
the natural ceiling of the class, NOT a certified endpoint-optimal member of it.  Stated here so the
report cannot overclaim it.

    python s31/s31_E2_deltas.py
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
    os.environ.setdefault(_v, "1")

from s31 import s31_E_lib as L             # noqa: E402

OUT_NPZ = os.path.join(HERE, "results", "s31_E2_deltas.npz")
OUT_JSON = os.path.join(HERE, "results", "s31_E2_deltas.json")
SEED = 310502

ORACLE_ARMS = {"E2o_PROJ_MU", "E3o_CONSTR_MU", "CTRL_SHRINK_E2o", "ORACLE_FULL",
               "ORACLE_Y_PERP_MU", "ORACLE_Y_ALONG_MU", "ORACLE_Y_PERP_1MU",
               "ORACLE_Y_PERP_MUHAT", "CTRL_SHRINK_ORACLE_Y"}


def p_along(v, u):
    uu = float(u @ u)
    return (float(v @ u) / uu) * u if uu > 1e-12 else np.zeros_like(v)


def p_perp(v, u):
    return v - p_along(v, u)


def resid_on(v, B):
    """Residual of v after least squares on the columns of B (per target)."""
    w, *_ = np.linalg.lstsq(B, v, rcond=None)
    return v - B @ w


def per_target_op(vec, direc, tgt, op):
    out = np.empty_like(vec)
    for _, m in L.per_target(vec, tgt):
        out[m] = op(vec[m], direc[m])
    return out


def main():
    d = L.load()
    X, y, mu, mh = d["X"], d["y"], d["mu"], d["mu_hat"]
    tgt, fold, sep = d["tgt"], d["fold"], d["sep"]
    pdbs = [str(p) for p in d["pdbs"]]
    feats = [str(f) for f in d["feats"]]
    rng = np.random.default_rng(SEED)

    # ---- S30's fitted corrector, reproduced bit-for-bit (lfo_ridge has no RNG)
    cols_n3 = [feats.index(f) for f in L.NEST["N3_plus_pool"]]
    yhat = L.lfo_ridge(X[:, cols_n3], y, fold)
    cols_n1 = [feats.index(f) for f in L.NEST["N1_plus_calibration"]]
    yhat1 = L.lfo_ridge(X[:, cols_n1], y, fold)
    ss = lambda a: float((a ** 2).sum())
    r2_n3 = 1 - ss(y - yhat) / ss(y - y.mean())
    r2_n1 = 1 - ss(y - yhat1) / ss(y - y.mean())

    # ---- per-target random direction with mu_hat's own marginal (matched control)
    rand_dir = np.empty_like(mh)
    for _, m in L.per_target(mh, tgt):
        rand_dir[m] = mh[m][rng.permutation(int(m.sum()))]

    # ---- constrained refits: residualise every feature within target against [1, u]
    def constrained_fit(u):
        Xr = np.empty_like(X)
        for _, m in L.per_target(y, tgt):
            B = np.column_stack([np.ones(int(m.sum())), u[m]])
            for k in range(X.shape[1]):
                Xr[m, k] = resid_on(X[m, k], B)
        yh = L.lfo_ridge(Xr[:, cols_n3], y, fold)
        return per_target_op(yh, u, tgt, p_perp)      # enforce the E2 constraint exactly

    def shrink_to(v, target_vec):
        out = np.empty_like(v)
        for _, m in L.per_target(v, tgt):
            nv = np.linalg.norm(v[m])
            out[m] = v[m] * (np.linalg.norm(target_vec[m]) / nv if nv > 1e-12 else 0.0)
        return out

    A = {}
    A["PROD"] = np.zeros_like(y)
    A["FIT_N3"] = yhat.copy()
    A["E2_PROJ_MUHAT"] = per_target_op(yhat, mh, tgt, p_perp)
    A["E2o_PROJ_MU"] = per_target_op(yhat, mu, tgt, p_perp)
    A["CTRL_PROJ_RAND"] = per_target_op(yhat, rand_dir, tgt, p_perp)
    A["CTRL_SHRINK_E2"] = shrink_to(yhat, A["E2_PROJ_MUHAT"])
    A["CTRL_SHRINK_E2o"] = shrink_to(yhat, A["E2o_PROJ_MU"])
    A["E3_CONSTR_MUHAT"] = constrained_fit(mh)
    A["E3o_CONSTR_MU"] = constrained_fit(mu)
    A["ORACLE_FULL"] = y.copy()
    A["ORACLE_Y_PERP_MU"] = per_target_op(y, mu, tgt, p_perp)
    A["ORACLE_Y_ALONG_MU"] = per_target_op(y, mu, tgt, p_along)
    A["ORACLE_Y_PERP_MUHAT"] = per_target_op(y, mh, tgt, p_perp)
    A["CTRL_SHRINK_ORACLE_Y"] = shrink_to(y, A["ORACLE_Y_PERP_MU"])
    yp1 = np.empty_like(y)
    for _, m in L.per_target(y, tgt):
        yp1[m] = resid_on(y[m], np.column_stack([np.ones(int(m.sum())), mu[m]]))
    A["ORACLE_Y_PERP_1MU"] = yp1

    # ---- diagnostics: R^2 as an approximation to y, coherence of the residual, norm fractions
    def coh(r):
        v = []
        for _, m in L.per_target(r, tgt):
            if r[m].std() > 1e-12 and mu[m].std() > 1e-12:
                v.append(float(np.corrcoef(r[m], mu[m])[0, 1]))
        return float(np.mean(v))

    def normfrac(v, u):
        f = []
        for _, m in L.per_target(v, tgt):
            nv = np.linalg.norm(v[m])
            f.append(float(np.linalg.norm(p_along(v[m], u[m])) / nv) if nv > 1e-12 else np.nan)
        return float(np.nanmean(f)), float(np.nanmedian(f))

    diag = {"R2_N3_reproduced": r2_n3, "R2_N1_reproduced": r2_n1,
            "S30_published_R2_N3": 0.23545312895964132,
            "S30_published_coh_N3": 0.9171667941971987,
            "S30_published_coh_uncorrected": 0.6930959406368322,
            "seed": SEED, "arms": {}}
    for k, v in A.items():
        r = y - v
        am, amed = normfrac(v, mu)
        hm, hmed = normfrac(v, mh)
        diag["arms"][k] = {
            "status": "ORACLE / NOT DEPLOYABLE" if k in ORACLE_ARMS else
                      ("baseline" if k == "PROD" else "LFO-supervised, native-free at inference"),
            "R2_vs_y": float(1 - ss(y - v) / ss(y - y.mean())),
            "coh_residual_with_mu_ORACLE": coh(r),
            "norm_frac_along_mu_mean": am, "norm_frac_along_mu_median": amed,
            "norm_frac_along_muhat_mean": hm, "norm_frac_along_muhat_median": hmed,
            "rms_delta": float(np.sqrt((v ** 2).mean()))}
    diag["coh_UNCORRECTED_reproduced"] = coh(y)

    tmp = OUT_NPZ[:-4] + ".%d.tmp.npz" % os.getpid()
    np.savez_compressed(tmp, off=d["off"], pdbs=d["pdbs"], names=np.array(sorted(A)),
                        **{"delta_" + k: v for k, v in A.items()})
    os.replace(tmp, OUT_NPZ)
    tj = OUT_JSON + ".%d.tmp" % os.getpid()
    with open(tj, "w") as fh:
        json.dump(diag, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tj, OUT_JSON)

    print("R2(N3) reproduced %.6f  (S30 published %.6f)  delta %.2e"
          % (r2_n3, diag["S30_published_R2_N3"], abs(r2_n3 - diag["S30_published_R2_N3"])))
    print("coh uncorrected reproduced %.6f  (S30 published %.6f)"
          % (diag["coh_UNCORRECTED_reproduced"], diag["S30_published_coh_uncorrected"]))
    print("%-22s %-40s %8s %8s %8s %8s" % ("arm", "status", "R2", "coh", "|along mu|", "rms"))
    for k in sorted(A):
        a = diag["arms"][k]
        print("%-22s %-40s %+8.4f %+8.4f %8.4f %8.4f"
              % (k, a["status"], a["R2_vs_y"], a["coh_residual_with_mu_ORACLE"],
                 a["norm_frac_along_mu_mean"], a["rms_delta"]))
    print("wrote", OUT_NPZ, OUT_JSON)


if __name__ == "__main__":
    main()
