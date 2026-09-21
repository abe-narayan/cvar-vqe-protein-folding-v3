#!/usr/bin/env python
"""s32/s32_P_inband.py -- S32 LANE P, H-P3: what in-band skill IS, and what it costs.

PREREG: s32/PREREG_S32_P.md @ 33dfe0d3.  Basis of every number here: the CA POINT CLOUD in
S31 lane A's own frame (the 128 score-top candidates superposed on the uniform medoid, the
native superposed likewise).  Nothing here is an endpoint claim; H-P3's only deployable
outcome (P3-6) is carried to the built chain by a separate job if it fires.

THE DECOMPOSITION, exact and asserted per target:

    a_k = |x_k - t|^2 = |mu|^2 + 2<mu, d_k> + |d_k|^2      mu = c - t,  d_k = x_k - c
          ORACLE label   constant   U_k                      V_k
                         in k       UNOBSERVABLE             OBSERVABLE

ORACLE / NOT DEPLOYABLE for every arm except the cos(mu_hat_DIS, mu) diagnostics, which read
the native only to SCORE a native-free direction.

    python s32/s32_P_inband.py run
    python s32/s32_P_inband.py analyse
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

from s12 import instrument as I                                        # noqa: E402
from scipy.stats import spearmanr, rankdata                            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s32_P_inband_rows.jsonl")
OUT = os.path.join(RESULTS, "s32_P_inband.json")

DIM = 128                 # the score-top window the readout acts on (S31 lane A's set)
NBAND = 24                # S31's in-band definition: the 24 truly best by the ORACLE label
SEED = 32_0_777
NPERM = 8                 # shuffle-U control draws
NDRAW = 6                 # price-curve draws per cos
COSGRID = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00]


def zs(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / max(v.std(), 1e-12)


def rho(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.std() < 1e-15 or y.std() < 1e-15:
        return float("nan")
    return float(spearmanr(x, y).statistic)


def dis_gradient(dg, C):
    """d(mean pairwise Bayes risk)/d(coords) at C (n,3), by finite difference on the risk table.

    The shipped score is `risk[p, bin(D_p)].mean(p)`; the table is a lookup on a 0.05 A grid, so
    the derivative wrt D_p is a central difference of the SAME table the score reads -- not a
    different function.  Chain rule to coordinates through D_p = |x_i - x_j|.
    """
    grid = np.asarray(dg["grid"], float); risk = np.asarray(dg["risk"], float)
    ii = np.asarray(dg["i"], int); jj = np.asarray(dg["j"], int)
    dv = C[ii] - C[jj]
    D = np.linalg.norm(dv, axis=1)
    g = np.clip(((D - grid[0]) / 0.05).astype(np.int64), 1, len(grid) - 2)
    ar = np.arange(risk.shape[0])
    drdD = (risk[ar, g + 1] - risk[ar, g - 1]) / (2.0 * 0.05)          # (npairs,)
    unit = dv / np.maximum(D, 1e-9)[:, None]
    G = np.zeros_like(C)
    np.add.at(G, ii, (drdD[:, None] * unit) / len(D))
    np.add.at(G, jj, -(drdD[:, None] * unit) / len(D))
    return G


def one(t, rng):
    pdb = t["pdb"]; n = int(t["n"])
    u = I.load_univ(pdb)
    dg = I.distogram(pdb, t["seq"], t["fold"])
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(float)
    pool = np.asarray(u["order"], int)[:500]
    W500 = np.asarray(u["W"], float)[pool]
    rr500 = np.asarray(u["rr"], float)[pool]
    nat = np.asarray(u["nat_ca"], float)

    order = np.argsort(DIS, kind="stable")
    top = order[:DIM]
    Wt = W500[top]
    Pt = I.pairwise_rmsd(Wt).astype(np.float32).astype(float)
    b = I.medoid(Pt)

    Sup = I.superpose_batch(Wt, Wt[b])
    tt = I.superpose_batch(nat[None], Wt[b])[0]
    Xf = Sup.reshape(DIM, -1); tf = tt.reshape(-1)
    a = ((Xf - tf[None, :]) ** 2).sum(1)                       # ORACLE label, A^2
    c = Xf.mean(0)
    mu = c - tf
    d = Xf - c[None, :]
    U = 2.0 * (d @ mu)
    V = (d ** 2).sum(1)
    mu2 = float(mu @ mu)
    ident = float(np.abs(a - (mu2 + U + V)).max())

    band = np.argsort(a, kind="stable")[:NBAND]
    #: the ORACLE band conditions on `a` and is therefore a collider; the band a DEPLOYED
    #: selector actually ranges over is the score's own top-24, which conditions on nothing
    #: native.  Both are carried so the collider claim can be separated from the deployability
    #: claim.  (`top` is already DIS-sorted, so the score's top-24 is 0..23.)
    nfband = np.arange(NBAND)

    # ------------------------------------------------------------------ P3-1 variance shares
    varU, varV = float(U.var()), float(V.var())
    covUV = float(np.cov(U, V)[0, 1])
    varU_b, varV_b = float(U[band].var()), float(V[band].var())
    covUV_b = float(np.cov(U[band], V[band])[0, 1])

    # ------------------------------------------------------------------ signals
    CONS = Pt.mean(1)                       # lane B's medoid criterion (native-free)
    DISz = zs(rankdata(DIS[top]))           # the deployed score's z-rank (native-free)
    sig = {"CONS": CONS, "V": V, "DIS": DISz, "ORACLE_UV": U + V}
    rhoG = {k: rho(v, a) for k, v in sig.items()}
    rhoIB = {k: rho(v[band], a[band]) for k, v in sig.items()}
    rhoNF = {k: rho(v[nfband], a[nfband]) for k, v in sig.items()}
    varU_n, varV_n = float(U[nfband].var()), float(V[nfband].var())

    # ------------------------------------------------------------------ P3-2 collider control
    # Shuffle U across candidates, breaking any U-V dependence but keeping both marginals.
    # If the negative in-band rho of a V-signal survives, it is the band's CONDITIONING and
    # carries no information about typicality being anti-informative.
    sh_V, sh_C = [], []
    for _ in range(NPERM):
        ash = mu2 + U[rng.permutation(DIM)] + V
        bs = np.argsort(ash, kind="stable")[:NBAND]
        sh_V.append(rho(V[bs], ash[bs])); sh_C.append(rho(CONS[bs], ash[bs]))

    # ------------------------------------------------------------------ P3-4 zero-mu control
    a0 = V.copy()                                    # the same pool with the common mode removed
    b0 = np.argsort(a0, kind="stable")[:NBAND]
    rho0 = {"CONS": rho(CONS[b0], a0[b0]), "DIS": rho(DISz[b0], a0[b0])}
    rho0_full = {"CONS": rho(CONS, a0), "DIS": rho(DISz, a0)}

    # ------------------------------------------------------------------ P3-3 mu recovery
    A = np.column_stack([np.ones(DIM), 2.0 * d])
    y = a - V
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = float(np.linalg.norm(A @ sol - y) / max(np.linalg.norm(y), 1e-12))
    mu_rec = sol[1:]
    # the recoverable part of mu is its projection on span(d); report that honestly
    Q, _ = np.linalg.qr(d.T)                                   # (3n, r) orthonormal basis of span(d)
    mu_par = Q @ (Q.T @ mu)
    rec_cos = float(mu_rec @ mu_par / max(np.linalg.norm(mu_rec) * np.linalg.norm(mu_par), 1e-12))
    rec_err = float(np.linalg.norm(mu_rec - mu_par) / max(np.linalg.norm(mu_par), 1e-12))
    ds = d[rng.permutation(DIM)]
    As = np.column_stack([np.ones(DIM), 2.0 * ds])
    ss, *_ = np.linalg.lstsq(As, y, rcond=None)
    res_shuf = float(np.linalg.norm(As @ ss - y) / max(np.linalg.norm(y), 1e-12))

    # ------------------------------------------------------------------ P3-5 price curve
    uhat = mu_par / max(np.linalg.norm(mu_par), 1e-12)
    nrm = float(np.linalg.norm(mu_par))
    price, price_nf = {}, {}
    for cs in COSGRID:
        vals, vnf = [], []
        for _ in range(NDRAW):
            g = Q @ rng.standard_normal(Q.shape[1])
            g = g - (g @ uhat) * uhat
            g /= max(np.linalg.norm(g), 1e-12)
            mh = nrm * (cs * uhat + np.sqrt(max(1.0 - cs * cs, 0.0)) * g)
            ah = 2.0 * (d @ mh) + V
            vals.append(rho(ah[band], a[band])); vnf.append(rho(ah[nfband], a[nfband]))
        price["%.2f" % cs] = [float(np.mean(vals)), float(np.std(vals, ddof=1)) if NDRAW > 1 else 0.0]
        price_nf["%.2f" % cs] = [float(np.mean(vnf)), float(np.std(vnf, ddof=1)) if NDRAW > 1 else 0.0]

    # ------------------------------------------------------------------ P3-6 native-free mu_hat
    Cn = c.reshape(n, 3)
    cand = {}
    G = dis_gradient(dg, Cn).reshape(-1)                       # descending risk moves toward t
    cand["GRAD"] = G
    cand["ARGMIN"] = c - Xf[int(np.argmin(DIS[top]))]          # best-scoring member -> centroid
    for T in (0.5, 2.0):
        w = np.exp(-zs(DIS[top]) / T); w /= w.sum()
        cand["SOFT%.1f" % T] = c - (w[:, None] * Xf).sum(0)    # softmax mean -> mu estimate
    q = DIM // 4
    cand["CONTRAST"] = Xf[3 * q:].mean(0) - Xf[:q].mean(0)     # worst-scoring - best-scoring
    _, _, Vt = np.linalg.svd(d, full_matrices=False)
    pc1 = Vt[0]
    cand["PC1"] = pc1 * np.sign(rho(d @ pc1, DISz) or 1.0)     # sign fixed by the score, not by mu
    # radial: the pool-vs-prediction Rg disagreement (memory: a native-free contrast)
    ii, jj = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
    Dc = np.linalg.norm(Cn[ii] - Cn[jj], axis=1)
    Dp = np.asarray(dg["expected"], float)
    s_rg = float(np.sum(Dc * Dp) / max(np.sum(Dc * Dc), 1e-12))
    cand["RADIAL"] = (Cn - Cn.mean(0)).reshape(-1) * (1.0 - s_rg)

    nf = {}
    for k, v in cand.items():
        vv = Q @ (Q.T @ v)
        nv = max(np.linalg.norm(vv), 1e-12)
        nf["cos_" + k] = float(vv @ mu_par / (nv * max(nrm, 1e-12)))
        nf["cosraw_" + k] = float(v @ mu / max(np.linalg.norm(v) * max(np.linalg.norm(mu), 1e-12), 1e-12))
        nf["rhoIB_dir_" + k] = rho((d @ vv)[band], a[band])
        nf["rhoIB_scaled_" + k] = rho((2.0 * (d @ (vv / nv * nrm)) + V)[band], a[band])
        nf["rhoNFB_dir_" + k] = rho((d @ vv)[nfband], a[nfband])
        nf["rhoNFB_scaled_" + k] = rho((2.0 * (d @ (vv / nv * nrm)) + V)[nfband], a[nfband])

    return dict(
        pdb=pdb, n=n, fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
        pool_best=float(rr500.min()), pool_mean=float(rr500.mean()),
        rr_top128_best=float(rr500[top].min()), rr_top128_mean=float(rr500[top].mean()),
        identity_max_abs=ident, mu2=mu2, meanV=float(V.mean()), rank_d=int(Q.shape[1]),
        f_common=float(mu2 / max((a).mean(), 1e-12)),
        varU=varU, varV=varV, covUV=covUV, varU_band=varU_b, varV_band=varV_b, covUV_band=covUV_b,
        ratio_all=float(varU / max(varV, 1e-12)), ratio_band=float(varU_b / max(varV_b, 1e-12)),
        varU_nfband=varU_n, varV_nfband=varV_n,
        ratio_nfband=float(varU_n / max(varV_n, 1e-12)),
        rho_global={k: v for k, v in rhoG.items()}, rho_inband={k: v for k, v in rhoIB.items()},
        rho_nfband={k: v for k, v in rhoNF.items()}, price_nfband=price_nf,
        shufU_rhoIB_V=[float(x) for x in sh_V], shufU_rhoIB_CONS=[float(x) for x in sh_C],
        zero_mu_rhoIB=rho0, zero_mu_rho_global=rho0_full,
        mu_recovery=dict(rel_resid=res, rel_resid_shuffled_d=res_shuf, cos_to_mu_par=rec_cos,
                         rel_err_to_mu_par=rec_err, mu_par_frac=float(nrm ** 2 / max(mu2, 1e-12))),
        price=price, native_free=nf,
    )


def run():
    done = set()
    if os.path.exists(ROWS):
        for ln in open(ROWS):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    ts = I.targets(); t0 = time.time()
    for k, t in enumerate(ts):
        if t["pdb"] in done:
            continue
        rng = np.random.default_rng(SEED + abs(hash(t["pdb"])) % 100003)
        r = one(t, rng)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        print("[%3d/%3d %5.0fs] %s ident=%.2e U/V=%5.2f band=%5.2f  rhoIB CONS %+0.3f V %+0.3f "
              "shufU_V %+0.3f  zeroMU CONS %+0.3f  murec %.1e  cosGRAD %+0.3f"
              % (k + 1, len(ts), time.time() - t0, r["pdb"], r["identity_max_abs"],
                 r["ratio_all"], r["ratio_band"], r["rho_inband"]["CONS"], r["rho_inband"]["V"],
                 float(np.mean(r["shufU_rhoIB_V"])), r["zero_mu_rhoIB"]["CONS"],
                 r["mu_recovery"]["rel_resid"], r["native_free"]["cos_GRAD"]), flush=True)
    print("run complete")


def _fold_stat(x, folds):
    x = np.asarray(x, float); m = np.array([np.nanmean(x[folds == f]) for f in sorted(set(folds))])
    se = float(m.std(ddof=1) / np.sqrt(len(m)))
    return float(np.nanmean(x)), se, [float(v) for v in m]


def analyse():
    R = [json.loads(ln) for ln in open(ROWS) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    if len(R) != 126:
        raise RuntimeError("have %d rows, not 126" % len(R))
    folds = np.array([r["fold"] for r in R])
    names = [r["pdb"] for r in R]
    g = lambda k: np.array([r[k] for r in R], float)                      # noqa: E731

    pm, pb = g("pool_mean"), g("pool_best")
    strata = {
        "all126": np.ones(126, bool),
        "worst18_poolmean": np.isin(np.arange(126), np.argsort(-pm)[:18]),
        "worst18_bestpool": np.isin(np.arange(126), np.argsort(-pb)[:18]),
        "FAIL18_CIRCULAR": np.array([r["fail18"] for r in R], bool),
    }
    other108 = ~strata["worst18_poolmean"]

    def agg(vals, mask=None):
        v = np.asarray(vals, float)
        m = np.ones(len(v), bool) if mask is None else mask
        mu_, se, per = _fold_stat(v[m], folds[m])
        return dict(mean=mu_, se=se, median=float(np.nanmedian(v[m])), n=int(m.sum()),
                    per_fold=per, frac_neg=float(np.nanmean(v[m] < 0)))

    ratio_all = g("ratio_all"); ratio_band = g("ratio_band")
    ibC = np.array([r["rho_inband"]["CONS"] for r in R])
    ibV = np.array([r["rho_inband"]["V"] for r in R])
    ibD = np.array([r["rho_inband"]["DIS"] for r in R])
    ibO = np.array([r["rho_inband"]["ORACLE_UV"] for r in R])
    gC = np.array([r["rho_global"]["CONS"] for r in R])
    gD = np.array([r["rho_global"]["DIS"] for r in R])
    shV = np.array([np.mean(r["shufU_rhoIB_V"]) for r in R])
    shV_sd = np.array([np.std(r["shufU_rhoIB_V"], ddof=1) for r in R])
    shC = np.array([np.mean(r["shufU_rhoIB_CONS"]) for r in R])
    z0C = np.array([r["zero_mu_rhoIB"]["CONS"] for r in R])
    z0D = np.array([r["zero_mu_rhoIB"]["DIS"] for r in R])

    def curve(key):
        pr = {}
        for cs in COSGRID:
            k = "%.2f" % cs
            m = np.array([r[key][k][0] for r in R]); s = np.array([r[key][k][1] for r in R])
            pr[k] = dict(mean=float(np.nanmean(m)), se=float(_fold_stat(m, folds)[1]),
                         median=float(np.nanmedian(m)), draw_sd_within_target=float(np.nanmean(s)),
                         frac_pos=float(np.nanmean(m > 0)))
        ys = np.array([pr["%.2f" % c]["mean"] for c in COSGRID]); xs = np.array(COSGRID)
        cs_ = float("nan")
        for i in range(len(xs) - 1):
            if ys[i] <= 0 <= ys[i + 1] and ys[i + 1] != ys[i]:
                cs_ = float(xs[i] + (xs[i + 1] - xs[i]) * (0 - ys[i]) / (ys[i + 1] - ys[i])); break
        return pr, cs_

    price, cstar = curve("price")
    price_nf, cstar_nf = curve("price_nfband")

    nfkeys = [k[4:] for k in R[0]["native_free"] if k.startswith("cos_")]
    nf = {}
    for k in nfkeys:
        cos = np.array([r["native_free"]["cos_" + k] for r in R])
        cosr = np.array([r["native_free"]["cosraw_" + k] for r in R])
        rdir = np.array([r["native_free"]["rhoIB_dir_" + k] for r in R])
        rsc = np.array([r["native_free"]["rhoIB_scaled_" + k] for r in R])
        rdn = np.array([r["native_free"]["rhoNFB_dir_" + k] for r in R])
        rsn = np.array([r["native_free"]["rhoNFB_scaled_" + k] for r in R])
        nf[k] = dict(cos_in_span=agg(cos), cos_raw=agg(cosr),
                     rhoIB_direction_only=agg(rdir), rhoIB_scaled_to_ORACLE_mu_norm=agg(rsc),
                     rhoNFB_direction_only=agg(rdn), rhoNFB_scaled_to_ORACLE_mu_norm=agg(rsn),
                     beats_cstar=bool(np.nanmean(cos) > cstar) if cstar == cstar else None)

    out = dict(
        prereg="s32/PREREG_S32_P.md @ 33dfe0d3", hypothesis="H-P3",
        basis="CA POINT CLOUD, S31 lane A frame (128 score-top, uniform medoid); NOT an endpoint",
        ORACLE="ORACLE / NOT DEPLOYABLE -- a_k, mu, U_k and every rho against a read the native",
        n=126, n_band=NBAND, dim=DIM, seed=SEED,
        identity_max_abs_over_targets=float(g("identity_max_abs").max()),
        P3_1_variance_share=dict(
            ratio_varU_over_varV_all128=agg(ratio_all), ratio_in_band=agg(ratio_band),
            median_ratio_all=float(np.median(ratio_all)), median_ratio_band=float(np.median(ratio_band)),
            frac_targets_ratio_gt_1_all=float((ratio_all > 1).mean()),
            frac_targets_ratio_gt_1_band=float((ratio_band > 1).mean()),
            f_common_top128=agg(g("f_common")),
            VERDICT="PASS" if np.median(ratio_all) > 1 else "FAIL"),
        P3_2_collider=dict(
            rhoIB_V=agg(ibV), rhoIB_CONS=agg(ibC), rhoIB_DIS=agg(ibD),
            rho_global_CONS=agg(gC), rho_global_DIS=agg(gD),
            shuffledU_rhoIB_V=agg(shV), shuffledU_rhoIB_CONS=agg(shC),
            shuffledU_draw_sd_within_target=float(shV_sd.mean()),
            real_minus_shuffled_V=agg(ibV - shV), real_minus_shuffled_CONS=agg(ibC - shC),
            CONS_minus_V_inband=agg(ibC - ibV),
            frac_targets_rhoIB_CONS_neg=float((ibC < 0).mean()),
            VERDICT="PASS" if (np.nanmean(ibC) < 0 and (ibC < 0).mean() >= 0.80) else "FAIL"),
        P3_3_mu_recovery=dict(
            rel_resid=agg(np.array([r["mu_recovery"]["rel_resid"] for r in R])),
            rel_resid_shuffled_d=agg(np.array([r["mu_recovery"]["rel_resid_shuffled_d"] for r in R])),
            frac_under_1e6=float(np.mean([r["mu_recovery"]["rel_resid"] < 1e-6 for r in R])),
            cos_recovered_to_mu_par=agg(np.array([r["mu_recovery"]["cos_to_mu_par"] for r in R])),
            mu_par_frac_of_mu2=agg(np.array([r["mu_recovery"]["mu_par_frac"] for r in R])),
            VERDICT="PASS" if np.mean([r["mu_recovery"]["rel_resid"] < 1e-6 for r in R]) >= 120 / 126 else "FAIL"),
        P3_4_zero_mu_control=dict(
            rhoIB_CONS_with_mu=agg(ibC), rhoIB_CONS_mu_removed=agg(z0C),
            rhoIB_DIS_mu_removed=agg(z0D),
            rho_global_CONS_mu_removed=agg(np.array([r["zero_mu_rho_global"]["CONS"] for r in R])),
            flip=agg(z0C - ibC),
            VERDICT="PASS" if np.nanmean(z0C) > 0.5 else "FAIL"),
        P3_5_price_curve=dict(grid=COSGRID, curve_ORACLE_band=price, cos_star_ORACLE_band=cstar,
                              curve_native_free_band=price_nf, cos_star_native_free_band=cstar_nf,
                              sanity_rhoIB_ORACLE_UV=agg(ibO),
                              VERDICT="PASS" if cstar == cstar and cstar > 0 else "FAIL"),
        NATIVE_FREE_BAND=dict(
            note="band = the SCORE's own top-24, which conditions on nothing native; this is the "
                 "set a deployed selector ranges over and it is NOT a collider",
            ratio_varU_over_varV=agg(g("ratio_nfband")),
            rho=dict(CONS=agg(np.array([r["rho_nfband"]["CONS"] for r in R])),
                     V=agg(np.array([r["rho_nfband"]["V"] for r in R])),
                     DIS=agg(np.array([r["rho_nfband"]["DIS"] for r in R])),
                     ORACLE_UV=agg(np.array([r["rho_nfband"]["ORACLE_UV"] for r in R])))),
        P3_6_native_free=nf,
        P5_strata={s: dict(ratio_all=agg(ratio_all, m), ratio_band=agg(ratio_band, m),
                           f_common=agg(g("f_common"), m), rhoIB_CONS=agg(ibC, m),
                           rhoIB_DIS=agg(ibD, m), mu2=agg(g("mu2"), m), meanV=agg(g("meanV"), m))
                   for s, m in strata.items()},
        P5_other108=dict(ratio_all=agg(ratio_all, other108), ratio_band=agg(ratio_band, other108),
                         f_common=agg(g("f_common"), other108), rhoIB_CONS=agg(ibC, other108)),
        strata_members={s: [names[i] for i in np.nonzero(m)[0]] for s, m in strata.items()
                        if s != "all126"},
    )
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "strata_members"}, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["run", "analyse"])
    a = ap.parse_args()
    (run if a.phase == "run" else analyse)()
