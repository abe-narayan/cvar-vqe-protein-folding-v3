#!/usr/bin/env python
"""s32/s32_Q1_sufficiency.py -- LANE Q, questions Q1/Q2/Q4.

Registered in `s32/PREREG_S32_Q.md` @ a8f9d6a7.

WHAT THIS MEASURES, and why it is a derivation with a numerical falsifier rather than an arm.

With `Sum_x w_x = 1`, `a_x = ||W_x - t||^2`, `B_xy = ||W_x - W_y||^2`:

    || Sum_x w_x W_x - t ||^2  =  <w,a> - (1/2) w'Bw                                   (I)

`B` is native-free; `a` is the only unknown.  But `a` is NOT a free vector:

    a_x = ||W_x||^2 - 2 <W_x, t> + ||t||^2                                             (II)

so `a` is AFFINE in `t`, and since `Sum w = 1` makes the `||t||^2` term an additive constant that
cannot move an argmin over the simplex, the readout-relevant part of `a` is the linear functional
`g_x = <W_x, t>` -- i.e. `P_aff{W} t`, the projection of the native onto the candidates' affine
hull, of dimension `rank(aff{W_x}) << 128`.

    THEOREM Q1-T1.  The per-candidate quality vector `a` and the pool common mode `mu = Xbar - t`
    are the SAME OBJECT up to a known, native-free affine bijection, and the sufficient statistic
    for the entire sum-to-one readout problem is `P_aff{W} t`.

That resolves S31 section 20.3, which states as OPEN whether `ahat` and `mu` are one requirement
or two.

    THEOREM Q1-T2 (sensitivity; the coordinator's question).  Write `x(t) = argmin_{w in simplex}
    ||U'w - t||^2` evaluated, i.e. the Euclidean projection of `t` onto conv{W_x}.  Then for a
    perturbation `dt` that does not change the active set `S`,

        dx  =  P_aff{W_x : x in S} dt            -- gain EXACTLY 1 inside, EXACTLY 0 outside

    so the readout reads `t` through an `(|S|-1)`-dimensional window and is blind to the rest, and
    globally `x` is 1-Lipschitz in `t` (projection onto a convex set).  Consequences:
      * an estimate error `eps` produces an emitted error in `[d, d+eps]` with `d = dist(t, conv)`;
      * the projection readout BEATS emitting the estimate directly iff the estimate's error
        exceeds `d` AND lies outside the hull's affine span.

Every arm here reads the native and is **ORACLE / NOT DEPLOYABLE**.  Basis is the **CA point
cloud**, a diagnostic; nothing here is differenced against the 3.2105 built-chain endpoint.

    python s32/s32_Q1_sufficiency.py [--limit N] [--draws 4]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import core                                                            # noqa: E402
from core.pipeline import consensus_medoid                             # noqa: E402
from s8 import consensus2 as cc                                        # noqa: E402

CACHE = os.path.join(ROOT, "s31", "results", "s31_C_cache")
OUT = os.path.join(HERE, "results")
NQ = 128


# ------------------------------------------------------------------ simplex QP
def proj_simplex(v):
    """Euclidean projection onto {w >= 0, sum w = 1} (Duchi et al.)."""
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    rho = np.nonzero(u - css / np.arange(1, len(v) + 1) > 0)[0][-1]
    return np.maximum(v - css[rho] / (rho + 1.0), 0.0)


def simplex_qp(U, t, rho=None, tol=1e-9, _G=None):
    """min_w ||U'w - t||^2 over the simplex, by NNLS on the sum-to-one-augmented system.

    WHY NOT FISTA, which is what this function did first.  A 4000-iteration FISTA plus a
    drop-only polish passed a 3-target smoke and then FAILED ITS OWN KKT CERTIFICATE on the
    full instrument -- max residual 1.26 at K=128 and 10.85 at K=500 -- because a drop-only
    polish never ADDS a violated index back, so a support that FISTA got wrong stays wrong.
    The certificate caught it; the smoke did not.  Lawson-Hanson NNLS is a finite exact
    active-set method, so the support is right by construction.

    The sum-to-one constraint is imposed by appending the row `sqrt(rho) * 1'` with target
    `sqrt(rho)`; `|sum w - 1|` is returned so the penalty is CHECKED, not assumed.

    Returns (w, support, kkt_residual, sum_err).  KKT for `min w'Gw - 2w'g` on the simplex is
    `(Gw - g)_i = nu` on the support and `>= nu` off it.
    """
    from scipy.optimize import nnls
    U = np.asarray(U, float)
    t = np.asarray(t, float)
    K = U.shape[0]
    G = U @ U.T if _G is None else _G
    g = U @ t
    if rho is None:
        rho = 1e6 * max(float(np.abs(G).max()), 1.0)
    A = np.vstack([U.T, np.sqrt(rho) * np.ones((1, K))])
    b = np.concatenate([t, [np.sqrt(rho)]])
    w, _res = nnls(A, b, maxiter=50 * K)
    S = np.flatnonzero(w > 1e-12)
    # exact equality-constrained polish on the NNLS support (removes the penalty's bias)
    for _ in range(60):
        s = len(S)
        M = np.zeros((s + 1, s + 1))
        M[:s, :s] = G[np.ix_(S, S)]
        M[:s, s] = -1.0
        M[s, :s] = 1.0
        sol, *_ = np.linalg.lstsq(M, np.concatenate([g[S], [1.0]]), rcond=None)
        wS = sol[:s]
        if (wS >= -1e-13).all():
            w = np.zeros(K)
            w[S] = np.maximum(wS, 0.0)
            break
        S = S[wS > -1e-13]
        if len(S) == 0:
            break
    r = G @ w - g
    S = np.flatnonzero(w > 1e-12)
    nu = float(r[S].mean()) if len(S) else 0.0
    kkt = max(float(np.abs(r[S] - nu).max()) if len(S) else 0.0,
              float(max(0.0, (nu - r).max())))
    kkt /= max(1.0, float(np.abs(r).max()))                 # RELATIVE, so the bar is scale-free
    return w, S, kkt, float(abs(w.sum() - 1.0))


# ------------------------------------------------------------------ per target
def do_target(pdb, draws, rng):
    geo = core.backend("geometry")
    z = np.load(os.path.join(CACHE, "%s.npz" % pdb))
    order = np.asarray(z["order"], int)
    W = np.asarray(z["W"], float)
    nat = np.asarray(z["nat"], float)                                  # ORACLE
    fold, n = int(z["fold"]), int(z["n"])
    o = order[:NQ]
    Wo = W[o]

    # ---- the deployed common frame: superpose everything onto the uniform consensus medoid
    Pt = geo.pairwise_ca_rmsd(Wo)
    b = consensus_medoid(Pt)
    ref = Wo[b]
    Sup = cc.superpose_batch(Wo, ref)
    tt = cc.superpose_batch(nat[None], ref)[0]                         # ORACLE
    U = Sup.reshape(NQ, -1)                                            # (K, d)
    t = tt.reshape(-1)
    d = U.shape[1]
    scale = float(np.sqrt(n))                                          # ||.|| -> RMSD units

    # ---- check 1: the readout identity, simplex AND affine draws
    a = ((U - t) ** 2).sum(1)
    sq = (U * U).sum(1)
    B = sq[:, None] + sq[None, :] - 2.0 * (U @ U.T)
    np.fill_diagonal(B, 0.0)
    err1 = 0.0
    for _ in range(40):
        w = rng.random(NQ)
        w = w / w.sum()                                                # simplex
        lhs = float(((U.T @ w - t) ** 2).sum())
        rhs = float(w @ a - 0.5 * w @ B @ w)
        err1 = max(err1, abs(lhs - rhs) / max(abs(lhs), 1e-12))
    for _ in range(40):
        w = rng.normal(size=NQ)
        w = w / w.sum()                                                # AFFINE, signed
        lhs = float(((U.T @ w - t) ** 2).sum())
        rhs = float(w @ a - 0.5 * w @ B @ w)
        err1 = max(err1, abs(lhs - rhs) / max(abs(lhs), 1e-12))

    # ---- check 2: a is affine in t
    g = U @ t
    a_from_t = sq - 2.0 * g + float(t @ t)
    err2 = float(np.abs(a_from_t - a).max() / max(np.abs(a).max(), 1e-12))

    # ---- check 3: the sufficient statistic.  Replace t by its affine-hull projection.
    Ub = U.mean(0)
    D = U - Ub
    _u, sv, vt = np.linalg.svd(D, full_matrices=False)
    rk = int((sv > sv[0] * 1e-10).sum())
    V = vt[:rk]                                                        # aff direction basis
    coef = V @ (t - Ub)
    t_par = Ub + V.T @ coef                                            # P_aff{W} t
    hull_gap = float(np.linalg.norm(t - t_par)) / scale                # RMSD units, ORACLE

    G = U @ U.T
    w_full, S_full, kkt_full, sum_full = simplex_qp(U, t, _G=G)
    w_par, S_par, kkt_par, _sp = simplex_qp(U, t_par, _G=G)
    err3 = float(np.abs(w_full - w_par).max())
    x_full = U.T @ w_full
    x_par = U.T @ w_par
    err3x = float(np.linalg.norm(x_full - x_par)) / scale

    # emitted values (ORACLE)
    hull_fixed = float(np.linalg.norm(x_full - t)) / scale             # fixed-frame
    hull_kabsch = float(geo.ca_rmsd(x_full.reshape(n, 3), nat))        # re-superposed

    # ---- check 4 (Q1-T2): dx = P_aff(S) dt, gain 1 inside / 0 outside
    A = U[S_full].T                                                    # (d, s)
    if len(S_full) > 1:
        Dc = (U[S_full] - U[S_full].mean(0)).T                         # (d, s) aff directions
        Qa, _r = np.linalg.qr(Dc)
        rs = int(np.linalg.matrix_rank(Dc, tol=1e-9))
        Qa = Qa[:, :rs]
        Pa = Qa @ Qa.T
    else:
        Pa = np.zeros((d, d))
        rs = 0
    h = 1e-7 * scale
    sens_rel, gain_in, gain_out = [], [], []
    n_moved = 0
    set0 = set(S_full.tolist())
    for _ in range(draws):
        dt = rng.normal(size=d)
        dt = h * dt / np.linalg.norm(dt)
        w2, S2, _k2, _s2 = simplex_qp(U, t + dt, _G=G)
        if set(S2.tolist()) != set0:      # Q1-T2 is stated FOR A FIXED ACTIVE SET
            n_moved += 1
            continue
        dx = U.T @ (w2 - w_full)
        pred = Pa @ dt
        sens_rel.append(float(np.linalg.norm(dx - pred) / max(np.linalg.norm(pred), 1e-30)))
        # pure in-hull and pure orthogonal directions
        din = Pa @ rng.normal(size=d)
        nin = float(np.linalg.norm(din))
        if nin > 1e-30:                       # |S| == 1 => the in-hull subspace is EMPTY
            din = h * din / nin               # and the gain is undefined, not zero
            w3, S3, _k3, _s3 = simplex_qp(U, t + din, _G=G)
            if set(S3.tolist()) == set0:
                gain_in.append(float(np.linalg.norm(U.T @ (w3 - w_full)) / np.linalg.norm(din)))
        dout = rng.normal(size=d)
        dout = dout - Pa @ dout
        dout = h * dout / max(np.linalg.norm(dout), 1e-30)
        w4, S4, _k4, _s4 = simplex_qp(U, t + dout, _G=G)
        if set(S4.tolist()) == set0:
            gain_out.append(float(np.linalg.norm(U.T @ (w4 - w_full)) / np.linalg.norm(dout)))

    # ---- Q4: the pricing curve in REAL NUMBERS.  ORACLE / NOT DEPLOYABLE.
    #      truncate t's representation to the top-r candidate-spread directions (the ORDERING
    #      is native-free -- it is the pool's own SVD -- only the coefficients are ORACLE).
    curve = []
    for r in sorted({0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 28, 32, rk}):
        if r > rk:
            continue
        tr = Ub + (V[:r].T @ coef[:r] if r else 0.0)
        wr, _Sr, _kr, _sr = simplex_qp(U, tr, _G=G)
        xr = U.T @ wr
        curve.append({"r": int(r),
                      "fixed": float(np.linalg.norm(xr - t)) / scale,
                      "kabsch": float(geo.ca_rmsd(xr.reshape(n, 3), nat))})

    # ---- the "how accurate must the estimate be" curve.  ORACLE / NOT DEPLOYABLE.
    #      isotropic error of RMSD magnitude eps, 4 draws; direct emission vs projected emission.
    noise = []
    for eps in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0):
        dir_v, prj_v = [], []
        for _ in range(draws):
            e = rng.normal(size=d)
            e = eps * scale * e / np.linalg.norm(e)
            th = t + e
            wn, _Sn, _kn, _sn = simplex_qp(U, th, _G=G)
            xn = U.T @ wn
            dir_v.append(float(geo.ca_rmsd(th.reshape(n, 3), nat)))
            prj_v.append(float(geo.ca_rmsd(xn.reshape(n, 3), nat)))
        noise.append({"eps": float(eps),
                      "direct_mean": float(np.mean(dir_v)), "direct_sd": float(np.std(dir_v, ddof=1)),
                      "proj_mean": float(np.mean(prj_v)), "proj_sd": float(np.std(prj_v, ddof=1))})

    # ---- Q1-C: does the s-of-K cardinality constraint BIND?  ORACLE / NOT DEPLOYABLE.
    #      If the UNCONSTRAINED convex optimum over the full K=500 pool is already sparser than
    #      s, then "best sparse convex combination" IS the convex program and is classically
    #      solved in milliseconds -- charter section 16 item 13 closes it without a circuit.
    W5 = cc.superpose_batch(W, ref)
    U5 = W5.reshape(len(W5), -1)
    G5 = U5 @ U5.T
    w5, S5, kkt5, sum5 = simplex_qp(U5, t, _G=G5)
    x5 = U5.T @ w5

    return {
        "pdb": pdb, "fold": fold, "n": n, "d": int(d),
        "k500_support": int(len(S5)), "k500_kkt": float(kkt5),
        "k500_fixed": float(np.linalg.norm(x5 - t)) / scale,
        "k500_kabsch": float(geo.ca_rmsd(x5.reshape(n, 3), nat)),
        "err_identity_rel": err1, "err_a_affine_in_t_rel": err2,
        "err_sufficient_w": err3, "err_sufficient_x_rmsd": err3x,
        "rank_aff": rk, "hull_gap_rmsd": hull_gap,
        "support": int(len(S_full)), "kkt": kkt_full, "kkt_par": kkt_par,
        "support_aff_dim": rs,
        "hull_fixed": hull_fixed, "hull_kabsch": hull_kabsch,
        "sens_rel_mean": float(np.mean(sens_rel)) if sens_rel else float("nan"),
        "sens_rel_max": float(np.max(sens_rel)) if sens_rel else float("nan"),
        "n_sens": int(len(sens_rel)), "n_activeset_moved": int(n_moved),
        "sum_err": float(sum_full), "k500_sum_err": float(sum5),
        "gain_in_mean": float(np.mean(gain_in)) if gain_in else float("nan"),
        "n_gain_in": int(len(gain_in)),
        "gain_out_mean": float(np.mean(gain_out)) if gain_out else float("nan"),
        "n_gain_out": int(len(gain_out)),
        "curve": curve, "noise": noise,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--draws", type=int, default=4)
    args = ap.parse_args()

    names = sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))
    if args.limit:
        names = names[:args.limit]

    rows, t0 = [], time.time()
    for i, pdb in enumerate(names):
        rng = np.random.default_rng(abs(hash(("s32Q1", pdb))) % (2 ** 32))
        rows.append(do_target(pdb, args.draws, rng))
        if i % 10 == 0:
            print("  %3d/%d %s  %.1fs" % (i + 1, len(names), pdb, time.time() - t0), flush=True)

    os.makedirs(OUT, exist_ok=True)
    tmp = os.path.join(OUT, "s32_Q1_rows.jsonl.%d.tmp" % os.getpid())
    with open(tmp, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    os.replace(tmp, os.path.join(OUT, "s32_Q1_rows.jsonl"))

    def col(k):
        return np.array([r[k] for r in rows], float)

    rs = sorted({c["r"] for r in rows for c in r["curve"]})
    curve = []
    for r in rs:
        v = [c for row in rows for c in row["curve"] if c["r"] == r]
        if len(v) == len(rows):
            curve.append({"r": r, "n": len(v),
                          "fixed_mean": float(np.mean([x["fixed"] for x in v])),
                          "kabsch_mean": float(np.mean([x["kabsch"] for x in v])),
                          "kabsch_se": float(np.std([x["kabsch"] for x in v], ddof=1) / np.sqrt(len(v)))})
    eps_list = sorted({c["eps"] for r in rows for c in r["noise"]})
    noise = []
    for e in eps_list:
        v = [c for row in rows for c in row["noise"] if c["eps"] == e]
        noise.append({"eps": e,
                      "direct_mean": float(np.mean([x["direct_mean"] for x in v])),
                      "proj_mean": float(np.mean([x["proj_mean"] for x in v])),
                      "proj_draw_sd_mean": float(np.mean([x["proj_sd"] for x in v]))})

    summary = {
        "prereg_commit": "a8f9d6a7",
        "label": "ORACLE / NOT DEPLOYABLE -- every arm reads the native",
        "basis": "CA point cloud (diagnostic).  NOT the 3.2105 built-chain endpoint.",
        "n": len(rows),
        "T1": {
            "err_identity_rel_max": float(col("err_identity_rel").max()),
            "err_a_affine_in_t_rel_max": float(col("err_a_affine_in_t_rel").max()),
            "err_sufficient_w_max": float(col("err_sufficient_w").max()),
            "err_sufficient_x_rmsd_max": float(col("err_sufficient_x_rmsd").max()),
            "rank_aff_mean": float(col("rank_aff").mean()),
            "rank_aff_min": int(col("rank_aff").min()),
            "rank_aff_max": int(col("rank_aff").max()),
            "d_mean": float(col("d").mean()),
            "hull_gap_rmsd_mean": float(col("hull_gap_rmsd").mean()),
            "kkt_max": float(col("kkt").max()),
        },
        "T2": {
            "support_mean": float(col("support").mean()),
            "support_median": float(np.median(col("support"))),
            "support_max": int(col("support").max()),
            "support_aff_dim_mean": float(col("support_aff_dim").mean()),
            "sens_rel_mean": float(col("sens_rel_mean").mean()),
            "sens_rel_max": float(col("sens_rel_max").max()),
            "gain_in_mean": float(col("gain_in_mean").mean()),
            "gain_in_sd": float(col("gain_in_mean").std(ddof=1)),
            "gain_out_mean": float(col("gain_out_mean").mean()),
            "gain_out_max": float(col("gain_out_mean").max()),
        },
        "hull": {
            "fixed_mean": float(col("hull_fixed").mean()),
            "kabsch_mean": float(col("hull_kabsch").mean()),
            "kabsch_se": float(col("hull_kabsch").std(ddof=1) / np.sqrt(len(rows))),
        },
        "Q1C_cardinality_ORACLE": {
            "k500_support_mean": float(col("k500_support").mean()),
            "k500_support_median": float(np.median(col("k500_support"))),
            "k500_support_max": int(col("k500_support").max()),
            "k500_support_p90": float(np.percentile(col("k500_support"), 90)),
            "frac_support_le_10": float((col("k500_support") <= 10).mean()),
            "k500_kabsch_mean": float(col("k500_kabsch").mean()),
            "k500_kabsch_se": float(col("k500_kabsch").std(ddof=1) / np.sqrt(len(rows))),
            "k500_kkt_max": float(col("k500_kkt").max()),
        },
        "Q4_curve_ORACLE": curve,
        "noise_ORACLE": noise,
        "elapsed_s": time.time() - t0,
    }
    tmp = os.path.join(OUT, "s32_Q1_sufficiency.json.%d.tmp" % os.getpid())
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    os.replace(tmp, os.path.join(OUT, "s32_Q1_sufficiency.json"))
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
