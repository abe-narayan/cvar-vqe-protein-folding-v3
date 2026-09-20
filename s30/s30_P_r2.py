#!/usr/bin/env python
"""s30/s30_P_r2.py -- S30 lane P.

THE ONE MEASUREMENT: can ANY native-free feature set explain the ORACLE error, out of fold?

Pre-registration `s30/PREREG_S30_P.md` (committed 097352c5, before this file existed).

Lane T's restatement of assumption B2:  rho_max = sqrt(R^2(e ~ S)).
    rho 0.14  (B2's ceiling)     <-> R^2  1.96 %
    rho 0.358 (3.00 A primary)   <-> R^2 12.82 %
    rho 0.628 (2.50 A ambitious) <-> R^2 39.44 %

e is ORACLE (the displacement from production's point cloud to the native).  EVERY feature is
native-free.  EVERY fit is leave-fold-out on the pinned folds.  A matched-dimension random-frame
control is refit alongside every arm and reported beside every number.

The basis carries lane L's S30-L18 result: the harmful mode is a per-target SEPARATION PROFILE
(five numbers, +52.5% of the gap, -0.525 A ORACLE), so the five per-separation-bin distogram
descent directions are IN the basis rather than left to the pool's PCs to span.

    python s30/s30_P_r2.py rows        # stage 1: per-target basis, response, features (cached)
    python s30/s30_P_r2.py fit         # stage 2: leave-fold-out ridge, controls, strata

Results: s30/results/s30_P_rows.npz, s30/results/s30_P.json
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A2_local as A2         # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s30_P_rows.npz")
OUT = os.path.join(RESULTS, "s30_P.json")

POOL_K, M, MIN_SEP = 500, 75, 2
NSEP = 5                 # lane L's five-number separation profile
KDIR = 14                # scale + consensus + 5 sep-gradients + 7 pool PCs
N_CTRL = 4               # matched-dimension random frames per target
SEED = 300201

REFS = ["consensus", "bestscore", "medoid", "distograd", "scale", "expand", "provenance",
        "rgmatch", "filterdisp", "bestresid", "scoredir", "rgdir", "typicdir"] + \
       ["sepgrad%d" % b for b in range(NSEP)]
SIGNED = REFS + ["skew", "corr_score", "corr_rg", "corr_sim", "corr_resid"]
DIST_DERIVED = {"bestscore", "distograd", "rgmatch", "filterdisp", "bestresid", "scoredir",
                "corr_score", "corr_resid"} | {"sepgrad%d" % b for b in range(NSEP)}
UNSDIR = ["eigshare", "projsd", "idx", "is_scale", "is_consensus", "is_sepgrad"]
TARGF = ["n", "score_mean", "score_sd", "score_gap", "top75_prmsd", "rg0", "rg_pool", "rg_pred",
         "rg_disagree", "dres_C0", "post_sd", "post_ent", "n_distinct", "sim_top75",
         "stable_rank", "pc1_share", "frac_pep"]


# --------------------------------------------------------------------------- helpers
def _unit(v):
    nv = np.linalg.norm(v)
    return v / nv if nv > 1e-12 else np.zeros_like(v)


def _rg(X):
    c = X - X.mean(0)
    return float(np.sqrt((c ** 2).sum(1).mean()))


def _skew(x):
    x = np.asarray(x, float); s = x.std()
    return float((((x - x.mean()) / s) ** 3).mean()) if s > 1e-12 else 0.0


def _corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.std() < 1e-12 or b.std() < 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _grad_ls(C0, ii, jj, exp, mask=None):
    """Descent direction of  sum_p (d_p - exp_p)^2  at C0, optionally restricted to a pair set."""
    if mask is not None:
        ii, jj, exp = ii[mask], jj[mask], exp[mask]
    v = C0[ii] - C0[jj]
    d = np.linalg.norm(v, axis=1)
    coef = (2.0 * (d - exp) / np.maximum(d, 1e-12))[:, None] * v
    g = np.zeros_like(C0)
    np.add.at(g, ii, coef)
    np.add.at(g, jj, -coef)
    return -g.ravel()


def frame_features(U, ctx):
    """Signed and unsigned per-direction features for an orthonormal frame U (K, 3n).

    Every entry of S flips sign with its direction; every entry of Zd is invariant.
    """
    V500, refs = ctx["V500"], ctx["refs"]
    score, rgm, sim, dres, sigma = ctx["score"], ctx["rgm"], ctx["sim"], ctx["dres"], ctx["sigma"]
    K = len(U)
    P = V500 @ U.T
    S = np.zeros((K, len(SIGNED))); Zd = np.zeros((K, len(UNSDIR)))
    for a, nm in enumerate(REFS):
        S[:, a] = U @ refs[nm]
    var = (P ** 2).mean(0); o = len(REFS)
    for k in range(K):
        p = P[:, k]
        S[k, o] = _skew(p)
        S[k, o + 1] = -_corr(p, score)
        S[k, o + 2] = _corr(p, rgm)
        S[k, o + 3] = _corr(p, sim)
        S[k, o + 4] = -_corr(p, dres)
        Zd[k] = [var[k] / max(var.sum(), 1e-12), np.sqrt(var[k]) / sigma, k / (K - 1.0),
                 float(ctx["tag"][k] == 1), float(ctx["tag"][k] == 2), float(ctx["tag"][k] == 3)]
    return S, Zd


# --------------------------------------------------------------------------- stage 1
def target_row(t, rng):
    pdb = t["pdb"]
    u = I.load_univ(pdb); rec = I.shipped_record(pdb); dg = I.distogram(pdb)
    n = int(u["n"])
    C0 = np.asarray(rec["avg_ca"], float)
    nat = np.asarray(u["nat_ca"], float)
    pool = np.asarray(u["order"], int)[:POOL_K]
    sub = np.asarray(rec["sub"], int)
    W = I.superpose_batch(np.asarray(u["W"], float)[pool], C0)
    B = A2.rigid_basis(C0)
    prj = lambda X: X - (X @ B) @ B.T
    V500 = prj((W - C0).reshape(POOL_K, -1))
    e = prj(A2.oracle_direction(C0, nat).reshape(1, -1))[0]            # ORACLE
    assert abs(np.sqrt((e ** 2).sum() / n) - rec["rmsd_avg"]) < 5e-3, pdb
    Vsub = V500[sub]
    sigma = float(np.sqrt((Vsub ** 2).sum(1).mean()))

    ii, jj = I.pair_index(n, MIN_SEP)
    exp = np.asarray(dg["expected"], float)
    sep = jj - ii
    D500 = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
    score = I.shipped_score(dg, D500)
    dres = np.abs(D500 - exp[None]).mean(1)
    rgm = np.array([_rg(w) for w in W])
    sim = np.asarray(u["sim"], float)[pool]
    org = np.asarray(u["org"], bool)[pool]
    dsq = np.zeros((n, n)); dsq[ii, jj] = exp ** 2; dsq[jj, ii] = exp ** 2
    for a in range(n - 1):
        dsq[a, a + 1] = dsq[a + 1, a] = 3.80 ** 2
    rg_pred = float(np.sqrt(dsq.sum() / (2.0 * n * n)))
    cen = C0.mean(0); rg0 = _rg(C0); rg_pool = float(rgm.mean())
    Pw = I.pairwise_rmsd(W[sub])
    consm = I.pairwise_rmsd(W).mean(1)                       # typicality of each pool member

    refs = {}
    refs["consensus"] = _unit(prj(Vsub.mean(0)[None])[0])
    refs["bestscore"] = _unit(V500[int(np.argmin(score))])
    refs["medoid"] = _unit(Vsub[int(I.medoid(Pw))])
    refs["distograd"] = _unit(prj(_grad_ls(C0, ii, jj, exp)[None])[0])
    Zs = _unit(prj((C0 - cen).ravel()[None])[0]); refs["scale"] = Zs
    refs["expand"] = _unit(prj(((C0 - cen) * (rg_pool / max(rg0, 1e-12) - 1.0)).ravel()[None])[0])
    pep = V500[org].mean(0) if org.any() else np.zeros_like(e)
    frg = V500[~org].mean(0) if (~org).any() else np.zeros_like(e)
    refs["provenance"] = _unit(prj((pep - frg)[None])[0])
    refs["rgmatch"] = _unit(V500[int(np.argmin(np.abs(rgm - rg_pred)))])
    refs["filterdisp"] = _unit(prj((V500.mean(0) - Vsub.mean(0))[None])[0])
    refs["bestresid"] = _unit(V500[int(np.argmin(dres))])
    zs = lambda v: (v - v.mean()) / max(v.std(), 1e-12)
    refs["scoredir"] = _unit(prj((-zs(score) @ V500)[None])[0])
    refs["rgdir"] = _unit(prj((zs(rgm) @ V500)[None])[0])
    refs["typicdir"] = _unit(prj((-zs(consm) @ V500)[None])[0])
    # lane L's five-number separation profile, transported to coordinate space
    edges = np.quantile(sep, np.linspace(0, 1, NSEP + 1))
    edges[0] -= 0.5; edges[-1] += 0.5
    sg = []
    for b in range(NSEP):
        mk = (sep > edges[b]) & (sep <= edges[b + 1])
        g = _unit(prj(_grad_ls(C0, ii, jj, exp, mk)[None])[0]) if mk.any() else np.zeros_like(e)
        refs["sepgrad%d" % b] = g; sg.append(g)

    # ---- the native-free basis
    Vc = Vsub - Vsub.mean(0)
    _, sv, Vt = np.linalg.svd(Vc, full_matrices=False)
    cand = [Zs, refs["consensus"]] + sg + [Vt[k] for k in range(min(7, Vt.shape[0]))]
    tag = [1, 2] + [3] * NSEP + [0] * 7
    Q, tg = [], []
    for v, tv in zip(cand, tag):
        w = v.copy()
        for q in Q:
            w = w - (w @ q) * q
        if np.linalg.norm(w) > 1e-8:
            Q.append(_unit(w)); tg.append(tv)
        if len(Q) == KDIR:
            break
    while len(Q) < KDIR:
        w = prj(rng.normal(size=e.shape)[None])[0]
        for q in Q:
            w = w - (w @ q) * q
        if np.linalg.norm(w) > 1e-8:
            Q.append(_unit(w)); tg.append(0)
    U = np.array(Q[:KDIR]); tg = tg[:KDIR]

    ctx = dict(V500=V500, refs=refs, score=score, rgm=rgm, sim=sim, dres=dres, sigma=sigma, tag=tg)
    S, Zd = frame_features(U, ctx)
    c = U @ e                                                            # ORACLE response

    cS, cZd, cc = [], [], []
    for _ in range(N_CTRL):
        R = prj(rng.normal(size=(KDIR, e.size)))
        Qr = []
        for v in R:
            w = v.copy()
            for q in Qr:
                w = w - (w @ q) * q
            Qr.append(_unit(w))
        Ur = np.array(Qr)
        ctxr = dict(ctx); ctxr["tag"] = [0] * KDIR
        s2, z2 = frame_features(Ur, ctxr)
        cS.append(s2); cZd.append(z2); cc.append(Ur @ e)

    x = [n, float(score.mean()), float(score.std()),
         float(np.sort(score)[M] - score.min()), float(Pw[np.triu_indices(M, 1)].mean()),
         rg0, rg_pool, rg_pred, rg_pred - rg0,
         float(np.abs(np.linalg.norm(C0[ii] - C0[jj], axis=1) - exp).mean()),
         float(np.asarray(dg["sd"], float).mean()),
         float((-(np.asarray(dg["prob"], float)
                  * np.log(np.clip(np.asarray(dg["prob"], float), 1e-12, None))).sum(1)).mean()),
         float(rec["n_distinct"]), float(sim[sub].mean()),
         float((sv ** 2).sum() ** 2 / max((sv ** 4).sum(), 1e-12)),
         float(sv[0] ** 2 / max((sv ** 2).sum(), 1e-12)), float(org[sub].mean())]

    rr = np.asarray(u["rr"], float)[pool]
    return dict(pdb=pdb, n=n, fold=int(t["fold"]), fail18=pdb in I.FAIL18,
                pool_mean=float(rr.mean()), pool_best=float(rr.min()),
                enorm2=float((e ** 2).sum()), sigma=sigma, c=c, S=S, Zd=Zd, x=np.array(x, float),
                cS=np.array(cS), cZd=np.array(cZd), cc=np.array(cc), tag=np.array(tg),
                prod=float(rec["rmsd_avg"]), prod_fit=float(rec["rmsd_fit"]))


def stage_rows(limit=None):
    rng = np.random.default_rng(SEED)
    tg = I.targets()[: limit or None]
    R, t0 = [], time.time()
    for k, t in enumerate(tg):
        R.append(target_row(t, rng))
        if (k + 1) % 20 == 0:
            print("  %d/%d  %.1fs" % (k + 1, len(tg), time.time() - t0), flush=True)
    out = {k: np.array([r[k] for r in R]) for k in
           ["c", "S", "Zd", "x", "cS", "cZd", "cc", "tag"]}
    for k in ["n", "fold", "fail18", "enorm2", "sigma", "prod", "prod_fit",
              "pool_mean", "pool_best"]:
        out[k] = np.array([r[k] for r in R])
    out["pdb"] = np.array([r["pdb"] for r in R])
    np.savez_compressed(ROWS, **out)
    print("wrote %s  (%d targets, %.1fs)" % (ROWS, len(R), time.time() - t0))


# --------------------------------------------------------------------------- stage 2
def design(S, Zd, x, use, inter):
    """Sign-equivariant design for one target.  NO intercept: an intercept is not equivariant."""
    A = S[:, use]
    if not inter:
        return A
    Z = np.concatenate([Zd, np.repeat(x[None], len(S), 0)], 1)
    return np.concatenate([A] + [A * Z[:, [b]] for b in range(Z.shape[1])], 1)


def _ridge(X, y, w, alpha):
    Xw = X * np.sqrt(w)[:, None]
    G = Xw.T @ Xw + alpha * np.eye(X.shape[1])
    return np.linalg.solve(G, Xw.T @ (y * np.sqrt(w)))


def _std(Xtr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd = np.where(sd > 1e-12, sd, 1.0)
    return (Xtr - mu) / sd, (Xte - mu) / sd


def lfo_predict(X, y, w, folds, alphas):
    """Leave-fold-out ridge, penalty chosen by a NESTED leave-fold-out inside the training folds.
    Standardisation is fit on training rows only."""
    F = np.array(sorted(set(folds.tolist())))
    yh = np.zeros_like(y)
    for f in F:
        tr, te = folds != f, folds == f
        best, ba = np.inf, alphas[0]
        for a in alphas:
            err = 0.0
            for g in F[F != f]:
                itr, ite = tr & (folds != g), folds == g
                A, Bx = _std(X[itr], X[ite])
                err += float((w[ite] * (y[ite] - Bx @ _ridge(A, y[itr], w[itr], a)) ** 2).sum())
            if err < best:
                best, ba = err, a
        A, Bx = _std(X[tr], X[te])
        yh[te] = Bx @ _ridge(A, y[tr], w[tr], ba)
    return yh


def run_arm(C, Sv, Zv, Xv, EN, SIG, folds_t, use, inter, weighted):
    nt, K = C.shape
    X = np.concatenate([design(Sv[t], Zv[t], Xv[t], use, inter) for t in range(nt)], 0)
    y = np.concatenate([C[t] / SIG[t] for t in range(nt)])
    w = np.concatenate([np.full(K, SIG[t] ** 2 if weighted else 1.0) for t in range(nt)])
    fd = np.concatenate([np.full(K, folds_t[t]) for t in range(nt)])
    yh = lfo_predict(X, y, w, fd, [1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5])
    ch = (yh.reshape(nt, K).T * SIG).T
    res = EN - (C ** 2).sum(1) + ((C - ch) ** 2).sum(1)
    return (EN - res) / EN, res, ch


def summarise(expl, res, EN, N, folds, pdbs, label, base_expl=None, base_res=None):
    r2 = float(1.0 - res.sum() / EN.sum())
    o = dict(label=label, R2_pooled=r2, R2_mean=float(expl.mean()),
             rho_pooled=float(np.sqrt(max(r2, 0.0))),
             rmsd_after_cloud=float(np.sqrt(res / N).mean()),
             rmsd_before_cloud=float(np.sqrt(EN / N).mean()))
    if base_expl is not None:
        # ST.compare's convention is "negative = a better", i.e. LOWER is better.  Explained
        # variance is HIGHER-is-better, so both arms are negated: `effect` then reads as a
        # deficit and a negative effect means the model beats the control.  Quoting the raw
        # verdict string on the un-negated arms would invert it -- lane D's S30-L5 trap.
        cmp = ST.compare(-np.asarray(expl), -np.asarray(base_expl), folds=folds, names=pdbs,
                         label=label + "|excess(neg: lower=better)")
        o["excess"] = {k: cmp[k] for k in ("effect", "se", "mde", "effect_over_mde", "ci95_fold",
                                           "ci95_iid", "folds_same_sign", "n_better", "n_worse",
                                           "verdict", "per_fold")}
        o["excess_R2_mean"] = float(-cmp["effect"])      # positive = model beats control
        o["R2_pooled_control"] = float(1.0 - base_res.sum() / EN.sum())
        o["excess_R2_pooled"] = r2 - o["R2_pooled_control"]
    return o


def stage_fit():
    z = np.load(ROWS, allow_pickle=True)
    C, Sv, Zv, Xv = z["c"], z["S"], z["Zd"], z["x"]
    cS, cZ, cc = z["cS"], z["cZd"], z["cc"]
    EN, SIG, folds = z["enorm2"], z["sigma"], z["fold"].astype(int)
    N = z["n"].astype(float); pdbs = list(z["pdb"]); nt = len(pdbs)
    f18 = z["fail18"].astype(bool); tag = z["tag"]
    tailA = np.zeros(nt, bool); tailA[np.argsort(-z["pool_mean"])[:18]] = True
    tailB = np.zeros(nt, bool); tailB[np.argsort(-z["pool_best"])[:18]] = True
    strata = {"FAIL18_CIRCULAR": f18, "other108": ~f18,
              "tailA_poolmean": tailA, "tailB_oraclebest": tailB}
    iall = list(range(len(SIGNED)))
    ind = [i for i, s in enumerate(SIGNED) if s not in DIST_DERIVED]
    capt = (C ** 2).sum(1) / EN
    capc = np.array([[(cc[t, r] ** 2).sum() / EN[t] for r in range(cc.shape[1])]
                     for t in range(nt)]).mean(1)
    res = {"meta": {"n_targets": nt, "K": int(C.shape[1]), "signed": SIGNED,
                    "no_dist_features": [SIGNED[i] for i in ind],
                    "bars": {"B2": 0.0196, "primary_3.00A": 0.1282, "ambitious_2.50A": 0.3944},
                    "prod_cloud": float(np.sqrt(EN / N).mean()),
                    "prod_builtchain": float(z["prod_fit"].mean()),
                    "ORACLE_capture_basis_mean": float(capt.mean()),
                    "ORACLE_capture_ctrl_mean": float(capc.mean()),
                    "ORACLE_capture_excess": float(capt.mean() - capc.mean()),
                    "ORACLE_capture_cmp": ST.compare(capt, capc, folds=folds, names=pdbs,
                                                     label="ORACLE capture basis vs ctrl")},
           "arms": {}, "strata": {}}
    arms = [("main", iall, False), ("full", iall, True),
            ("main_nodist", ind, False), ("full_nodist", ind, True)]
    for wtd in (True, False):
        for nm, use, inter in arms:
            key = "%s|%s" % (nm, "wt" if wtd else "unwt")
            ex, rs, ch = run_arm(C, Sv, Zv, Xv, EN, SIG, folds, use, inter, wtd)
            cex, crs = [], []
            for r in range(cc.shape[1]):
                e2, r2, _ = run_arm(cc[:, r], cS[:, r], cZ[:, r], Xv, EN, SIG, folds, use,
                                    inter, wtd)
                cex.append(e2); crs.append(r2)
            cex = np.mean(cex, 0); crs = np.mean(crs, 0)
            a = summarise(ex, rs, EN, N, folds, pdbs, key, cex, crs)
            sc = C[:, 0]; scc = ch[:, 0]
            sg = np.array([np.isin(np.where(tag[t] == 3)[0], np.arange(C.shape[1]))
                           for t in range(nt)])
            a["scale_mode"] = {
                "share_of_error_ORACLE": float((sc ** 2).sum() / EN.sum()),
                "R2_of_scale_coeff": float(1 - ((sc - scc) ** 2).sum() / (sc ** 2).sum()),
                "R2_shape_pooled": float(1 - (rs.sum() - ((sc - scc) ** 2).sum())
                                         / (EN.sum() - (sc ** 2).sum())),
                "mean_scale_coeff": float(sc.mean()), "sd_scale_coeff": float(sc.std())}
            sep_idx = [k for k in range(C.shape[1]) if (tag[:, k] == 3).mean() > 0.5]
            if sep_idx:
                cs, cshat = C[:, sep_idx], ch[:, sep_idx]
                a["sepprofile_mode"] = {
                    "dirs": sep_idx,
                    "share_of_error_ORACLE": float((cs ** 2).sum() / EN.sum()),
                    "R2_of_sep_coeffs": float(1 - ((cs - cshat) ** 2).sum() / (cs ** 2).sum())}
            res["arms"][key] = a
            for sn, mask in strata.items():
                res["strata"].setdefault(key, {})[sn] = {
                    "n": int(mask.sum()),
                    "R2_pooled": float(1 - rs[mask].sum() / EN[mask].sum()),
                    "R2_pooled_control": float(1 - crs[mask].sum() / EN[mask].sum()),
                    "excess_R2_mean": float(ex[mask].mean() - cex[mask].mean()),
                    "ORACLE_capture": float(capt[mask].mean())}
            print("%-18s R2 %+.4f  ctrl %+.4f  excess %+.4f  %s"
                  % (key, a["R2_pooled"], a["R2_pooled_control"], a["excess_R2_pooled"],
                     a["excess"]["verdict"]), flush=True)

    rng = np.random.default_rng(SEED + 7)
    perm = np.arange(nt)
    for nv in np.unique(N):
        idx = np.where(N == nv)[0]
        perm[idx] = rng.permutation(idx)
    exP, rsP, _ = run_arm(C, Sv[perm], Zv[perm], Xv[perm], EN, SIG, folds, iall, True, True)
    res["arms"]["CPERM|full|wt"] = summarise(exP, rsP, EN, N, folds, pdbs, "CPERM")
    print("CPERM             R2 %+.4f" % res["arms"]["CPERM|full|wt"]["R2_pooled"])

    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("wrote", OUT)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["rows", "fit"])
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    stage_rows(a.limit) if a.stage == "rows" else stage_fit()


if __name__ == "__main__":
    main()
