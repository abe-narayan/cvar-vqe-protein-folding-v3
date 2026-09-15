#!/usr/bin/env python
"""s27/s28_C_fail18.py -- S28 lane C, Part 1: a FAIL18 detector on a feature class no router has
seen (`s27/PREREG_S28_C.md` section 1).

NATIVE-FREE features per target (`features_one`): the distribution of the three universe-fitted
statistical potentials (DISTPOT, ENV, CONTACT of `s27/ham_lib.py`) over the shipped K = 500 pool
and their agreement with the shipped distogram score (block SP), plus the control classes that
previous routers used (block CTRL).  Nothing in `features_one` reads `nat_ca`, `oracle_rr` or an
RMSD; the test suite NaN-poisons both and asserts bit-identical features.

ORACLE label (read only inside `run`): `s12.instrument.FAIL18` membership, the 18 zero-recall
targets.  It is the training label of a nested leave-fold-out ridge-logistic classifier and is
never read at inference.

    python s27/s28_C_fail18.py features            # native-free, cached to s27/results
    python s27/s28_C_fail18.py selftest            # synthetic labels, no RMSD read
    python s27/s28_C_fail18.py run [--n-perm 500]  # AUROC vs permutation null; the ORACLE switch
    python s27/s28_C_fail18.py switch              # the switched arm (only if F1 fired)
"""
from __future__ import annotations

import argparse
import json
import math
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
from s27 import ham_lib as HL              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

RESULTS = os.path.join(HERE, "results")
FEAT = os.path.join(RESULTS, "s28_C_features.json")
OUT = os.path.join(RESULTS, "s28_C_fail18.json")
SWITCH_OUT = os.path.join(RESULTS, "s28_C_switch.json")
S26_FEAT = os.path.join(ROOT, "s26", "results", "p_c4_features.json")
CHAIN_ROWS = os.path.join(RESULTS, "chain_rows.jsonl")
M = 75
SP_CHANNELS = ("DISTPOT", "ENV", "CONTACT")
ALPHAS = np.logspace(-2, 4, 13)
N_PERM = 500
N_RAND = 200


# ============================================================================ features
def _skew(x):
    x = np.asarray(x, float)
    c = x - x.mean()
    s = c.std()
    return float((c ** 3).mean() / s ** 3) if s > 0 else 0.0


def _bin_entropy(p):
    p = np.clip(np.asarray(p, float), 1e-12, 1 - 1e-12)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))


def features_one(cand, ch, sim, dg):
    """One target's NATIVE-FREE feature dict.  `cand` supplies W (the pool) and n; `ch` the
    S27 channels; `sim` the pool's BLOSUM sims; `dg` the distogram.  Never `nat_ca` / `rr`."""
    from scipy.stats import spearmanr
    n = int(cand.n)
    k = int(cand.k)
    key = RP.rng_for(cand.pdb, "tiekey").random(k)
    dis = np.asarray(ch["DIS"], float)
    top = RP.topm(dis, M, key)
    n_pairs3 = (n - 3) * (n - 2) // 2
    f = {}
    for X in SP_CHANNELS:
        x = np.asarray(ch[X], float)
        q25, q75 = np.percentile(x, [25, 75])
        size = float(n) if X == "ENV" else float(n_pairs3)
        z = HL.zrank(x)
        tx = RP.topm(x, M, key)
        f[f"{X}_spread"] = float((q75 - q25) / size)
        f[f"{X}_skew"] = _skew(x)
        f[f"{X}_rho_dis"] = float(spearmanr(x, dis).correlation) if np.std(x) > 0 else 0.0
        f[f"{X}_overlap"] = float(len(set(tx.tolist()) & set(top.tolist())) / M)
        f[f"{X}_top75_z"] = float(z[top].mean())
        f[f"{X}_top75_sd"] = float(z[top].std())
    cons = np.asarray(ch["CONS"], float)
    f["cons_mean"] = float(cons.mean())
    f["cons_top75_gap"] = float((cons[top].mean() - cons.mean()) / max(cons.std(), 1e-12))
    W = np.asarray(cand.W, float)
    i, j = I.pair_index(n, 3)
    D = I.pair_dists(W, i, j)
    fr = (D < 8.0).mean(0)
    f["pool_contact_ent"] = float(_bin_entropy(fr).mean())
    # CTRL block
    f["n"] = float(n)
    pr = np.clip(np.asarray(dg["prob"], float), 1e-12, 1)
    f["dg_ent_mean"] = float(-(pr * np.log(pr)).sum(1).mean())
    sim = np.asarray(sim, float)
    zs = sim - sim.max()
    w = np.exp(zs)
    w /= w.sum()
    f["sim_ent"] = float(-(w * np.log(np.clip(w, 1e-300, 1))).sum() / math.log(len(w)))
    f["sim_gap75"] = float((sim[0] - sim[M - 1]) / (sim.std() + 1e-9))
    P = I.pairwise_rmsd(W[top])
    f["top75_spread"] = float(P[np.triu_indices(M, 1)].mean())
    return f


SP_NAMES = [f"{X}_{s}" for X in SP_CHANNELS for s in ("spread", "skew", "rho_dis", "overlap", "top75_z", "top75_sd")] \
    + ["cons_mean", "cons_top75_gap", "pool_contact_ent"]
CTRL_NAMES = ["n", "dg_ent_mean", "sim_ent", "sim_gap75", "top75_spread"]


def build_features(verbose=True, limit=0):
    if os.path.exists(FEAT) and not limit:
        return json.load(open(FEAT))
    from s25 import phys_lib as P
    pdbs = P.targets()
    if limit:
        pdbs = pdbs[:limit]
    rows = []
    t0 = time.time()
    for c, pdb in enumerate(pdbs):
        cand, ch, _ = RP.channels_for(pdb)
        u = I.load_univ(pdb)
        sim = np.asarray(u["sim"][I.pool_idx(u)], float)
        del u
        dg = I.distogram(pdb, cand.seq, cand.fold)
        f = features_one(cand, ch, sim, dg)
        rows.append(dict(pdb=pdb, fold=int(cand.fold), **f))
        if verbose:
            print(f"  features [{c+1}/{len(pdbs)}] {pdb} n={cand.n} ({time.time()-t0:.0f}s)", flush=True)
    names = SP_NAMES + CTRL_NAMES
    out = {"rows": rows, "names": names, "sp": SP_NAMES, "ctrl": CTRL_NAMES}
    if not limit:
        ST.save_atomic(FEAT, out, complete_keys=names, rows=rows, n_expected=126, module_file=__file__)
    return out


def blocks():
    """Feature blocks: the new class, the controls, both, and the S26 comparison blocks."""
    feat = build_features(verbose=False)
    pdbs = [r["pdb"] for r in feat["rows"]]
    folds = np.array([int(r["fold"]) for r in feat["rows"]])
    X = {nm: np.array([float(r[nm]) for r in feat["rows"]]) for nm in feat["names"]}
    s26 = json.load(open(S26_FEAT))
    s26rows = {r["pdb"]: r for r in s26["rows"]}
    assert all(p in s26rows for p in pdbs)
    for nm in s26["names"]:
        X["s26:" + nm] = np.array([float(s26rows[p][nm]) for p in pdbs])
    B = {
        "SP": list(SP_NAMES),
        "CTRL": list(CTRL_NAMES),
        "SP+CTRL": list(SP_NAMES) + list(CTRL_NAMES),
        "old_S22": ["s26:" + nm for nm in s26["names"] if nm.startswith("old_")],
        "s26_new_all": ["s26:" + nm for nm in s26["names"] if not nm.startswith("old_")],
    }
    return pdbs, folds, X, B


# ============================================================================ the classifier
def auroc(y, s):
    """Rank-based AUROC (ties averaged)."""
    from scipy.stats import rankdata
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    pos = y > 0.5
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def logistic_fit(X, y, alpha, iters=50, tol=1e-9):
    """L2-penalised logistic regression, balanced class weights, standardised columns
    (training statistics), unpenalised intercept.  Newton iterations."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    mu = X.mean(0)
    sd = X.std(0) + 1e-9
    Z = (X - mu) / sd
    n, p = Z.shape
    n1 = max(int((y > 0.5).sum()), 1)
    n0 = max(n - n1, 1)
    w = np.where(y > 0.5, n / (2.0 * n1), n / (2.0 * n0))
    A = np.hstack([np.ones((n, 1)), Z])
    beta = np.zeros(p + 1)
    pen = np.full(p + 1, float(alpha))
    pen[0] = 0.0

    def loss(b):
        eta = A @ b
        # weighted log-loss, numerically stable: log(1+exp(eta)) - y*eta
        ll = np.logaddexp(0.0, eta) - y * eta
        return float((w * ll).sum() + 0.5 * (pen * b * b).sum())

    cur = loss(beta)
    for _ in range(iters):
        eta = A @ beta
        pr = 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
        g = A.T @ (w * (pr - y)) + pen * beta
        Hm = (A * (w * pr * (1 - pr))[:, None]).T @ A + np.diag(pen)
        step = np.linalg.solve(Hm + 1e-10 * np.eye(p + 1), g)
        t = 1.0
        nxt = beta - t * step
        nl = loss(nxt)
        while nl > cur + 1e-12 and t > 1e-4:          # damped Newton: halve until it descends
            t *= 0.5
            nxt = beta - t * step
            nl = loss(nxt)
        beta, prev, cur = nxt, cur, nl
        if float(np.abs(t * step).max()) < tol or abs(prev - cur) < 1e-10 * max(1.0, abs(cur)):
            break
    return {"beta": beta, "mu": mu, "sd": sd}


def logistic_decision(m, X):
    Z = (np.asarray(X, float) - m["mu"]) / m["sd"]
    return m["beta"][0] + Z @ m["beta"][1:]


def nested_logistic(X, y, folds, alphas=ALPHAS, fold_order=None):
    """Held-out decision values for every target.  For each outer fold f: alpha chosen by inner
    leave-one-fold-out over the other folds (pooled inner held-out AUROC; ties -> larger alpha),
    refit on all other folds, fold f scored.  Also returns the prevalence-matched threshold of
    each fold (the (1 - training prevalence) quantile of the training decision values)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    folds = np.asarray(folds)
    F = sorted(set(folds.tolist()))
    if fold_order is not None:
        F = list(fold_order)
    dec = np.zeros(len(y))
    thr = np.zeros(len(y))
    chosen = {}
    for f in F:
        tr = folds != f
        inner = [g for g in F if g != f]
        best, best_a = -np.inf, alphas[-1]
        for a in alphas:
            s_in = np.zeros(len(y))
            for g in inner:
                itr = tr & (folds != g)
                ite = folds == g
                m = logistic_fit(X[itr], y[itr], a)
                s_in[ite] = logistic_decision(m, X[ite])
            sc = auroc(y[tr], s_in[tr])
            if sc >= best - 1e-12:
                best, best_a = sc, a
        m = logistic_fit(X[tr], y[tr], best_a)
        dec[folds == f] = logistic_decision(m, X[folds == f])
        d_tr = logistic_decision(m, X[tr])
        prev = float(y[tr].mean())
        thr[folds == f] = float(np.quantile(d_tr, 1.0 - prev))
        chosen[int(f)] = float(best_a)
    return dec, thr, chosen


def nested_single(x, y, folds):
    """The one-feature rule: sign chosen on the training folds by AUROC; held-out decision =
    signed feature.  Threshold: the prevalence-matched quantile of the training rows."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    folds = np.asarray(folds)
    dec = np.zeros(len(y))
    thr = np.zeros(len(y))
    for f in sorted(set(folds.tolist())):
        tr = folds != f
        s = 1.0 if auroc(y[tr], x[tr]) >= 0.5 else -1.0
        dec[folds == f] = s * x[folds == f]
        prev = float(y[tr].mean())
        thr[folds == f] = float(np.quantile(s * x[tr], 1.0 - prev))
    return dec, thr


def per_fold_auroc(y, dec, folds):
    vals = []
    for f in sorted(set(np.asarray(folds).tolist())):
        m = np.asarray(folds) == f
        v = auroc(y[m], dec[m])
        if np.isfinite(v):
            vals.append(v)
    return float(np.mean(vals)) if vals else float("nan")


def balanced_acc(y, pred):
    y = np.asarray(y) > 0.5
    pred = np.asarray(pred) > 0.5
    tpr = float(pred[y].mean()) if y.any() else float("nan")
    tnr = float((~pred[~y]).mean()) if (~y).any() else float("nan")
    return 0.5 * (tpr + tnr)


# ============================================================================ the switch endpoint
def chain_table():
    """Per-target built-chain and point-cloud RMSD of the named configurations, from S27's
    `chain_rows.jsonl` (ORACLE numbers; read only by `run` / `switch`)."""
    by = {}
    with open(CHAIN_ROWS, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            by.setdefault(r["config"], {})[r["pdb"]] = r
    return by


def switched_endpoint(pred_fail, pdbs, by, alt, basis="rmsd_chain"):
    base = np.array([by["DIS"][p][basis] for p in pdbs])
    altv = np.array([by[alt][p][basis] for p in pdbs])
    return np.where(np.asarray(pred_fail, bool), altv, base), base


def reproduce_chain(n_check=6):
    """Re-run the production projection of the DIS and DIS+DISTPOT top-75 averages on the first
    `n_check` targets and compare with `chain_rows.jsonl` (the reproduction check of 1.6)."""
    from s24 import d_harness as H
    from s27 import run_vqe_chain as RV
    from s25 import phys_lib as P
    by = chain_table()
    out = []
    for pdb in P.targets()[:n_check]:
        cand, ch, _ = RP.channels_for(pdb)
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        for cfg in ("DIS", "DIS+DISTPOT"):
            E = RV.energy_for(cfg, ch, pdb, key)
            top = np.lexsort((key, E))[:M]
            C, _ = H.readout_uniform(cand, top)
            ca = H.readout_projected(cand, C)
            rc = float(I.ca_rmsd(C, cand.nat_ca))
            rk = float(I.ca_rmsd(ca, cand.nat_ca))
            out.append(dict(pdb=pdb, config=cfg, rmsd_cloud=rc, rmsd_chain=rk,
                            ref_cloud=by[cfg][pdb]["rmsd_cloud"], ref_chain=by[cfg][pdb]["rmsd_chain"],
                            d_cloud=abs(rc - by[cfg][pdb]["rmsd_cloud"]),
                            d_chain=abs(rk - by[cfg][pdb]["rmsd_chain"])))
    return out


# ============================================================================ run
def run(n_perm=N_PERM, seed=0, fold_order=None, tag=""):
    pdbs, folds, X, B = blocks()
    fail = set(I.FAIL18)                                   # ORACLE label
    y = np.array([1.0 if p in fail else 0.0 for p in pdbs])
    assert int(y.sum()) == 18
    rng = np.random.default_rng(seed)
    res = {"n": len(pdbs), "n_pos": int(y.sum()), "seed": seed, "n_perm": n_perm,
           "fold_order": list(fold_order) if fold_order else None, "blocks": {}, "singles": {}}
    t0 = time.time()
    for bname, names in B.items():
        Xb = np.column_stack([X[nm] for nm in names])
        dec, thr, chosen = nested_logistic(Xb, y, folds, fold_order=fold_order)
        obs = auroc(y, dec)
        pf = per_fold_auroc(y, dec, folds)
        pred05 = dec >= 0.0
        predpm = dec >= thr
        null = np.empty(n_perm)
        null_pred05 = np.zeros((n_perm, len(y)), bool)
        null_predpm = np.zeros((n_perm, len(y)), bool)
        for t in range(n_perm):
            yp = y[rng.permutation(len(y))]
            d_p, t_p, _ = nested_logistic(Xb, yp, folds, fold_order=fold_order)
            null[t] = auroc(yp, d_p)
            null_pred05[t] = d_p >= 0.0
            null_predpm[t] = d_p >= t_p
        p_perm = float((null >= obs).mean())
        res["blocks"][bname] = dict(
            names=names, n_features=len(names), auroc=obs, auroc_per_fold_mean=pf,
            bal_acc_05=balanced_acc(y, pred05), bal_acc_pm=balanced_acc(y, predpm),
            n_pred_pos_05=int(pred05.sum()), n_pred_pos_pm=int(predpm.sum()),
            tp_05=int((pred05 & (y > 0.5)).sum()), tp_pm=int((predpm & (y > 0.5)).sum()),
            alphas=chosen, null_mean=float(null.mean()), null_p50=float(np.percentile(null, 50)),
            null_p95=float(np.percentile(null, 95)), null_p99=float(np.percentile(null, 99)),
            p_perm=p_perm, F1_fires=bool(obs > np.percentile(null, 95)),
            dec=dec.tolist(), thr=thr.tolist(), pred05=pred05.tolist(), predpm=predpm.tolist(),
            null=null.tolist())
        np.save(os.path.join(RESULTS, f"s28_C_null_pred_{bname.replace('+', '_')}{tag}.npy"),
                np.stack([null_pred05, null_predpm]))
        print(f"  block {bname:12s} p={len(names):2d}  AUROC {obs:.3f} (per-fold mean {pf:.3f})  "
              f"null mean {null.mean():.3f} p95 {np.percentile(null, 95):.3f}  p_perm {p_perm:.3f}  "
              f"balacc@0.5 {balanced_acc(y, pred05):.3f} (pred+ {int(pred05.sum())}, TP {int((pred05 & (y > 0.5)).sum())})  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        ST.save_atomic(OUT, res, module_file=__file__)
    # singles: every feature, sign chosen on the training folds; a shared null of the signed rule
    for nm in list(SP_NAMES) + list(CTRL_NAMES):
        dec, thr = nested_single(X[nm], y, folds)
        obs = auroc(y, dec)
        null = np.empty(n_perm)
        for t in range(n_perm):
            yp = y[rng.permutation(len(y))]
            d_p, _ = nested_single(X[nm], yp, folds)
            null[t] = auroc(yp, d_p)
        res["singles"][nm] = dict(auroc=obs, auroc_raw=auroc(y, X[nm]), null_p95=float(np.percentile(null, 95)),
                                  p_perm=float((null >= obs).mean()), bal_acc_pm=balanced_acc(y, dec >= thr))
    ST.save_atomic(OUT, res, module_file=__file__)
    print("  singles (held-out AUROC, sign nested; raw AUROC of the feature; perm p):")
    for nm, v in sorted(res["singles"].items(), key=lambda kv: -kv[1]["auroc"]):
        print(f"    {nm:20s} {v['auroc']:.3f}  raw {v['auroc_raw']:.3f}  p95 {v['null_p95']:.3f}  p {v['p_perm']:.3f}")
    # the ORACLE switch (the ceiling), independent of the detector
    by = chain_table()
    res["oracle_switch"] = {}
    for alt in ("DIS+DISTPOT", "DIS+ENV"):
        for basis in ("rmsd_chain", "rmsd_cloud"):
            a, b = switched_endpoint(y > 0.5, pdbs, by, alt, basis)
            r = ST.compare(a, b, folds, names=pdbs, label=f"ORACLE switch FAIL18->{alt} vs production ({basis})")
            print(ST.fmt(r))
            res["oracle_switch"][f"{alt}|{basis}"] = {k: v for k, v in r.items() if k != "concentration"} | {"concentration": r["concentration"]}
            # strata on this basis
            d = a - b
            res["oracle_switch"][f"{alt}|{basis}"]["fail18_effect"] = float(d[y > 0.5].mean())
            res["oracle_switch"][f"{alt}|{basis}"]["fail18_se"] = float(d[y > 0.5].std(ddof=1) / math.sqrt(18))
    ST.save_atomic(OUT, res, module_file=__file__)
    return res


def switch(block="SP+CTRL", alt="DIS+DISTPOT", rule="05", n_rand=N_RAND, seed=0, res=None):
    """The switched arm (run only if F1 fired for `block`).  Predicted-FAIL -> alt's top-75,
    else DIS.  Controls: random subsets of the same size; the label-permutation endpoints."""
    res = res or json.load(open(OUT))
    pdbs, folds, X, B = blocks()
    by = chain_table()
    fail = set(I.FAIL18)
    y = np.array([1.0 if p in fail else 0.0 for p in pdbs])
    pred = np.asarray(res["blocks"][block]["pred" + rule], bool)
    out = {"block": block, "alt": alt, "rule": rule, "n_pred_pos": int(pred.sum()),
           "tp": int((pred & (y > 0.5)).sum()), "fp": int((pred & (y < 0.5)).sum())}
    rng = np.random.default_rng(seed + 1)
    nullp = np.load(os.path.join(RESULTS, f"s28_C_null_pred_{block.replace('+', '_')}.npy"))[0 if rule == "05" else 1]
    for basis in ("rmsd_chain", "rmsd_cloud"):
        a, b = switched_endpoint(pred, pdbs, by, alt, basis)
        r = ST.compare(a, b, folds, names=pdbs, label=f"switched[{block},{rule}] -> {alt} vs production ({basis})")
        print(ST.fmt(r))
        d = a - b
        rand = []
        for _ in range(n_rand):
            pk = np.zeros(len(y), bool)
            pk[rng.choice(len(y), size=int(pred.sum()), replace=False)] = True
            ar, br = switched_endpoint(pk, pdbs, by, alt, basis)
            rand.append(float((ar - br).mean()))
        perm = []
        for t in range(min(n_rand, len(nullp))):
            ap, bp = switched_endpoint(nullp[t], pdbs, by, alt, basis)
            perm.append(float((ap - bp).mean()))
        rand = np.array(rand)
        perm = np.array(perm)
        out[basis] = {k: v for k, v in r.items()}
        out[basis].update(dict(
            fail18_effect=float(d[y > 0.5].mean()), non_fail18_effect=float(d[y < 0.5].mean()),
            rand_mean=float(rand.mean()), rand_p05=float(np.percentile(rand, 5)), rand_p95=float(np.percentile(rand, 95)),
            rand_share_below=float((rand <= d.mean()).mean()),
            perm_mean=float(perm.mean()), perm_p05=float(np.percentile(perm, 5)), perm_share_below=float((perm <= d.mean()).mean())))
        print(f"    FAIL18 {d[y > 0.5].mean():+.4f}  non-FAIL18 {d[y < 0.5].mean():+.4f}  "
              f"random-subset control mean {rand.mean():+.4f} [p05 {np.percentile(rand, 5):+.4f}]  "
              f"label-perm endpoint mean {perm.mean():+.4f} [p05 {np.percentile(perm, 5):+.4f}]")
    ST.save_atomic(SWITCH_OUT, out, module_file=__file__)
    return out


def singles_max_null(n_perm=N_PERM, seed=0, res=None):
    """The best single feature is an order statistic over the 26 singles (contract rule 6,
    `grid-oracles-are-order-statistics`): price it against the distribution of the MAXIMUM
    held-out AUROC over the same 26 signed rules under label permutation.  Also a labelled
    diagnostic: the top singles residualised on length n (native-free), to see whether they
    are length proxies."""
    res = res or json.load(open(OUT))
    pdbs, folds, X, B = blocks()
    fail = set(I.FAIL18)                                   # ORACLE label
    y = np.array([1.0 if p in fail else 0.0 for p in pdbs])
    names = list(SP_NAMES) + list(CTRL_NAMES)
    rng = np.random.default_rng(seed + 7)
    obs = {nm: res["singles"][nm]["auroc"] for nm in names}
    best = max(obs, key=obs.get)
    mx = np.empty(n_perm)
    for t in range(n_perm):
        yp = y[rng.permutation(len(y))]
        mx[t] = max(auroc(yp, nested_single(X[nm], yp, folds)[0]) for nm in names)
    out = {"best_single": best, "best_auroc": obs[best], "max_null_mean": float(mx.mean()),
           "max_null_p50": float(np.percentile(mx, 50)), "max_null_p95": float(np.percentile(mx, 95)),
           "p_max": float((mx >= obs[best]).mean()), "n_singles": len(names), "n_perm": n_perm}
    print(f"  best single {best} AUROC {obs[best]:.3f}; max-over-{len(names)} null mean {mx.mean():.3f} "
          f"p50 {np.percentile(mx, 50):.3f} p95 {np.percentile(mx, 95):.3f}; p_max {out['p_max']:.3f}")
    # length-residualised diagnostic for the top five
    n = X["n"]
    diag = {}
    for nm in sorted(obs, key=obs.get, reverse=True)[:5]:
        if nm == "n":
            continue
        x = X[nm]
        A = np.column_stack([np.ones(len(n)), n])
        resid = x - A @ np.linalg.lstsq(A, x, rcond=None)[0]
        dec, _ = nested_single(resid, y, folds)
        rho = float(np.corrcoef(x, n)[0, 1])
        diag[nm] = {"auroc_residual_on_n": auroc(y, dec), "corr_with_n": rho}
        print(f"    {nm:20s} corr(x, n) {rho:+.3f}   held-out AUROC after residualising on n {auroc(y, dec):.3f}")
    out["length_residual_diag"] = diag
    res["singles_max_null"] = out
    ST.save_atomic(OUT, res, module_file=__file__)
    return out


def selftest():
    """Synthetic: a planted signal must clear its permutation null; random labels must not.
    No RMSD, no native read."""
    rng = np.random.default_rng(0)
    n, p = 126, 8
    folds = np.repeat(np.arange(5), 26)[:n]
    X = rng.normal(size=(n, p))
    w = rng.normal(size=p)
    lin = X @ w + 0.7 * rng.normal(size=n)
    y_sig = (lin > np.quantile(lin, 1 - 18 / 126)).astype(float)
    y_rnd = np.zeros(n)
    y_rnd[rng.choice(n, 18, replace=False)] = 1.0
    ds, _, _ = nested_logistic(X, y_sig, folds)
    dr, _, _ = nested_logistic(X, y_rnd, folds)
    null = np.array([auroc(yp, nested_logistic(X, yp, folds)[0])
                     for yp in (y_rnd[rng.permutation(n)] for _ in range(30))])
    print(f"  planted: AUROC {auroc(y_sig, ds):.3f}; random labels {auroc(y_rnd, dr):.3f}; null p95 {np.percentile(null, 95):.3f}")
    assert auroc(y_sig, ds) > 0.85 and auroc(y_sig, ds) > np.percentile(null, 95)
    assert auroc(y_rnd, dr) < max(np.percentile(null, 99), 0.7)
    print("  s28_C_fail18 selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["features", "selftest", "run", "switch", "reproduce", "singles_null"])
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--reverse-folds", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--block", default="SP+CTRL")
    ap.add_argument("--alt", default="DIS+DISTPOT")
    ap.add_argument("--rule", default="05")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    if a.mode == "features":
        build_features(limit=a.limit)
    elif a.mode == "selftest":
        selftest()
    elif a.mode == "reproduce":
        rows = reproduce_chain()
        for r in rows:
            print(f"  {r['pdb']} {r['config']:12s} cloud {r['rmsd_cloud']:.6f} (ref {r['ref_cloud']:.6f}, d {r['d_cloud']:.1e})  "
                  f"chain {r['rmsd_chain']:.6f} (ref {r['ref_chain']:.6f}, d {r['d_chain']:.1e})")
        ST.save_atomic(os.path.join(RESULTS, "s28_C_reproduce.json"), {"rows": rows}, module_file=__file__)
    elif a.mode == "run":
        global OUT
        if a.reverse_folds or a.seed:
            OUT = OUT.replace(".json", f"_seed{a.seed}{'_rev' if a.reverse_folds else ''}.json")
        run(n_perm=a.n_perm, seed=a.seed, fold_order=[4, 3, 2, 1, 0] if a.reverse_folds else None,
            tag=f"_seed{a.seed}{'_rev' if a.reverse_folds else ''}" if (a.reverse_folds or a.seed) else "")
    elif a.mode == "switch":
        switch(block=a.block, alt=a.alt, rule=a.rule, seed=a.seed)
    elif a.mode == "singles_null":
        singles_max_null(n_perm=a.n_perm, seed=a.seed)


if __name__ == "__main__":
    main()
